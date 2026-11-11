#pragma once

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <optional>
#include <span>
#include <utility>
#include <vector>

#include <geom/brep/face_query.hpp>
#include <geom/brep/face_uv.hpp>
#include <geom/cgeom/arrangement2.hpp>
#include <geom/cgeom/polygon2.hpp>
#include <geom/cgeom/triangulate2.hpp>
#include <geom/core/error.hpp>
#include <geom/core/result.hpp>
#include <geom/primitives/point.hpp>
#include <geom/primitives/vec.hpp>

namespace geom {
inline namespace v1 {

// One curve along which a curved face is split. The polyline lives in the
// surface parameter domain reported by face_loops_uv() for the same face. A
// closed cut repeats its first point as its last point.
struct FaceCut {
  std::vector<Point2d> uv;
};

// One piece of a split curved face. `loops` holds the 3D rings and `loops_uv`
// the matching parameter rings; index 0 is the outer ring and the rest are
// holes. The rings are implicitly closed, outer is counter-clockwise and holes
// are clockwise. `reversed` is always false here: the boolean driver flips the
// B-side fragments when it needs the opposite material side.
struct CurvedFragment {
  std::vector<std::vector<Point3d>> loops;
  std::vector<std::vector<Point2d>> loops_uv;
  bool reversed = false;
};

namespace detail {

inline bool trim_face_finite2(const Point2d& p) noexcept {
  return std::isfinite(p[0]) && std::isfinite(p[1]);
}

inline bool trim_face_finite3(const Point3d& p) noexcept {
  return std::isfinite(p[0]) && std::isfinite(p[1]) && std::isfinite(p[2]);
}

inline bool trim_face_lex_less(const Point2d& a, const Point2d& b) noexcept {
  if (a[0] != b[0]) return a[0] < b[0];
  return a[1] < b[1];
}

inline bool trim_face_same_xy(const Point2d& a, const Point2d& b) noexcept {
  return a[0] == b[0] && a[1] == b[1];
}

inline Point2d trim_face_project_on_edge(const Point2d& p, const Point2d& a,
                                         const Point2d& b) noexcept {
  const Vec2d edge = b - a;
  const double length2 = dot(edge, edge);
  if (!(length2 > 0.0)) return a;
  double t = dot(p - a, edge) / length2;
  if (!(t > 0.0)) {
    t = 0.0;
  } else if (t >= 1.0) {
    t = 1.0;
  }
  return Point2d{a[0] + edge[0] * t, a[1] + edge[1] * t};
}

// Where a cut endpoint landed relative to the projected face boundary.
struct TrimFaceSnap {
  Point2d point{};
  bool boundary = false;
  bool vertex = false;
  std::size_t ring = 0;
  std::size_t edge = 0;
};

// Snaps a point to the nearest boundary vertex or boundary segment within
// `radius`. A vertex wins a tie with a segment. Points already on a segment are
// projected onto it, so the caller can insert them as exact ring vertices.
inline TrimFaceSnap trim_face_find_snap(
    const std::vector<std::vector<Point2d>>& loops, const Point2d& p,
    double radius) {
  TrimFaceSnap snap;
  snap.point = p;
  bool found = false;
  double best = 0.0;

  for (std::size_t r = 0; r < loops.size(); ++r) {
    for (const Point2d& vertex : loops[r]) {
      const double distance = std::hypot(vertex[0] - p[0], vertex[1] - p[1]);
      if (!(distance <= radius)) continue;
      if (found && !(distance < best)) continue;
      found = true;
      best = distance;
      snap.point = vertex;
      snap.boundary = true;
      snap.vertex = true;
      snap.ring = r;
      snap.edge = 0;
    }
  }

  for (std::size_t r = 0; r < loops.size(); ++r) {
    const std::vector<Point2d>& ring = loops[r];
    for (std::size_t i = 0; i < ring.size(); ++i) {
      const Point2d& a = ring[i];
      const Point2d& b = ring[(i + 1) % ring.size()];
      const double distance = face_query_point_segment_distance(p, a, b);
      if (!(distance <= radius)) continue;
      // A vertex already at the same distance keeps priority.
      if (found && !(distance < best)) continue;
      found = true;
      best = distance;
      snap.point = trim_face_project_on_edge(p, a, b);
      snap.boundary = true;
      snap.vertex = false;
      snap.ring = r;
      snap.edge = i;
      // A projection that lands on an end vertex is that vertex.
      if (trim_face_same_xy(snap.point, a)) {
        snap.point = a;
        snap.vertex = true;
      } else if (trim_face_same_xy(snap.point, b)) {
        snap.point = b;
        snap.vertex = true;
      }
    }
  }
  return snap;
}

inline bool trim_face_snap_better(const TrimFaceSnap& candidate,
                                  const TrimFaceSnap& current) noexcept {
  if (candidate.boundary != current.boundary) return candidate.boundary;
  if (candidate.vertex != current.vertex) return candidate.vertex;
  return trim_face_lex_less(candidate.point, current.point);
}

inline std::size_t trim_face_find_root(std::vector<std::size_t>& parent,
                                       std::size_t index) noexcept {
  while (parent[index] != index) {
    parent[index] = parent[parent[index]];
    index = parent[index];
  }
  return index;
}

struct TrimFaceInsert {
  std::size_t ring = 0;
  std::size_t edge = 0;
  Point2d point{};
};

// Inserts snapped points into the ring edges they landed on. All inserts refer
// to the original ring indices, so each ring is rebuilt once.
inline void trim_face_apply_inserts(std::vector<std::vector<Point2d>>& rings,
                                    std::vector<TrimFaceInsert> inserts) {
  if (inserts.empty()) return;
  std::sort(inserts.begin(), inserts.end(),
            [](const TrimFaceInsert& lhs, const TrimFaceInsert& rhs) noexcept {
              if (lhs.ring != rhs.ring) return lhs.ring < rhs.ring;
              if (lhs.edge != rhs.edge) return lhs.edge < rhs.edge;
              return trim_face_lex_less(lhs.point, rhs.point);
            });
  inserts.erase(std::unique(inserts.begin(), inserts.end(),
                            [](const TrimFaceInsert& lhs,
                               const TrimFaceInsert& rhs) noexcept {
                              return lhs.ring == rhs.ring &&
                                     lhs.edge == rhs.edge &&
                                     trim_face_same_xy(lhs.point, rhs.point);
                            }),
                inserts.end());

  std::size_t index = 0;
  while (index < inserts.size()) {
    const std::size_t ring_index = inserts[index].ring;
    std::size_t group_end = index;
    while (group_end < inserts.size() && inserts[group_end].ring == ring_index) {
      ++group_end;
    }
    std::vector<Point2d>& ring = rings[ring_index];
    const std::size_t vertex_count = ring.size();
    std::vector<std::vector<std::pair<double, Point2d>>> extra(vertex_count);
    for (std::size_t i = index; i < group_end; ++i) {
      const std::size_t edge = inserts[i].edge;
      const Point2d& a = ring[edge];
      const Point2d& b = ring[(edge + 1) % vertex_count];
      const Vec2d edge_vector = b - a;
      const double length2 = dot(edge_vector, edge_vector);
      const double t =
          length2 > 0.0 ? dot(inserts[i].point - a, edge_vector) / length2 : 0.0;
      extra[edge].emplace_back(t, inserts[i].point);
    }
    std::size_t extra_total = 0;
    for (std::vector<std::pair<double, Point2d>>& edge_extra : extra) {
      std::sort(edge_extra.begin(), edge_extra.end(),
                [](const std::pair<double, Point2d>& lhs,
                   const std::pair<double, Point2d>& rhs) noexcept {
                  if (lhs.first != rhs.first) return lhs.first < rhs.first;
                  return trim_face_lex_less(lhs.second, rhs.second);
                });
      extra_total += edge_extra.size();
    }
    std::vector<Point2d> rebuilt;
    rebuilt.reserve(vertex_count + extra_total);
    for (std::size_t i = 0; i < vertex_count; ++i) {
      rebuilt.push_back(ring[i]);
      for (const std::pair<double, Point2d>& entry : extra[i]) {
        rebuilt.push_back(entry.second);
      }
    }
    ring.swap(rebuilt);
    index = group_end;
  }
}

// Walks the single boundary cycle owned by an arrangement face.
inline bool trim_face_walk_cycle(const Arrangement2& arrangement,
                                 std::uint32_t face,
                                 std::vector<Point2d>& ring) {
  const std::span<const Arrangement2::HalfEdge> half_edges =
      arrangement.half_edges();
  const std::span<const Point2d> vertices = arrangement.vertices();
  const std::uint32_t start = arrangement.faces()[face].edge;
  std::uint32_t half_edge = start;
  std::size_t guard = 0;
  do {
    ring.push_back(vertices[half_edges[half_edge].origin]);
    half_edge = half_edges[half_edge].next;
    if (++guard > half_edges.size()) return false;
  } while (half_edge != start);
  return true;
}

// Removes consecutive duplicate vertices, which a bridge cycle produced by a
// keyhole cut can contain.
inline void trim_face_drop_repeats(std::vector<Point2d>& ring) {
  if (ring.size() < 2) return;
  std::vector<Point2d> cleaned;
  cleaned.reserve(ring.size());
  for (const Point2d& point : ring) {
    if (!cleaned.empty() && trim_face_same_xy(cleaned.back(), point)) continue;
    cleaned.push_back(point);
  }
  while (cleaned.size() >= 2 && trim_face_same_xy(cleaned.front(), cleaned.back())) {
    cleaned.pop_back();
  }
  ring.swap(cleaned);
}

// Reverses the ring when needed and rotates it so the lexicographically
// smallest vertex comes first. This makes the ring start independent of the
// arrangement's internal vertex ids.
inline void trim_face_normalize_ring(std::vector<Point2d>& ring, bool want_ccw) {
  trim_face_drop_repeats(ring);
  if (ring.size() < 3) return;
  const int orientation = ring_orientation2(ring);
  if (orientation != 0 && ((orientation > 0) != want_ccw)) {
    std::reverse(ring.begin(), ring.end());
  }
  std::size_t first = 0;
  for (std::size_t i = 1; i < ring.size(); ++i) {
    if (trim_face_lex_less(ring[i], ring[first])) first = i;
  }
  if (first != 0) {
    std::rotate(ring.begin(), ring.begin() + static_cast<std::ptrdiff_t>(first),
                ring.end());
  }
}

// Area scale |dS/du x dS/dv| of the surface at `uv`, computed with a centered
// difference. It converts a parameter-domain area into a surface area. Falls
// back to one when the derivative is not usable.
inline double trim_face_metric_scale(const FaceUvProjection& projection,
                                     const Point2d& uv) noexcept {
  const double h =
      1e-6 * std::max(1.0, std::max(std::fabs(uv[0]), std::fabs(uv[1])));
  const Point3d u_minus = face_uv_eval(projection, Point2d{uv[0] - h, uv[1]});
  const Point3d u_plus = face_uv_eval(projection, Point2d{uv[0] + h, uv[1]});
  const Point3d v_minus = face_uv_eval(projection, Point2d{uv[0], uv[1] - h});
  const Point3d v_plus = face_uv_eval(projection, Point2d{uv[0], uv[1] + h});
  if (!trim_face_finite3(u_minus) || !trim_face_finite3(u_plus) ||
      !trim_face_finite3(v_minus) || !trim_face_finite3(v_plus)) {
    return 1.0;
  }
  const Vec3d du{(u_plus[0] - u_minus[0]) / (2.0 * h),
                 (u_plus[1] - u_minus[1]) / (2.0 * h),
                 (u_plus[2] - u_minus[2]) / (2.0 * h)};
  const Vec3d dv{(v_plus[0] - v_minus[0]) / (2.0 * h),
                 (v_plus[1] - v_minus[1]) / (2.0 * h),
                 (v_plus[2] - v_minus[2]) / (2.0 * h)};
  const Vec3d normal = cross(du, dv);
  const double scale = std::sqrt(dot(normal, normal));
  if (!std::isfinite(scale) || !(scale > 0.0)) return 1.0;
  return scale;
}

struct TrimFaceCycle {
  std::vector<Point2d> ring;
  std::uint32_t face = 0;
  double area = 0.0;
  Point2d rep{};
  bool has_rep = false;
  bool material = false;
  std::size_t parent = static_cast<std::size_t>(-1);
  std::vector<std::size_t> children;
};

inline bool trim_face_same_fragment(const CurvedFragment& lhs,
                                    const CurvedFragment& rhs) noexcept {
  if (lhs.loops_uv.size() != rhs.loops_uv.size()) return false;
  for (std::size_t r = 0; r < lhs.loops_uv.size(); ++r) {
    if (lhs.loops_uv[r].size() != rhs.loops_uv[r].size()) return false;
    for (std::size_t i = 0; i < lhs.loops_uv[r].size(); ++i) {
      if (!trim_face_same_xy(lhs.loops_uv[r][i], rhs.loops_uv[r][i])) {
        return false;
      }
    }
  }
  return true;
}

}  // namespace detail

// Splits `face` along `cuts` and returns the resulting fragments in a
// deterministic order (lexicographic minimum of the outer parameter ring, then
// hole count). Every cut polyline is given in the parameter domain of the
// face's surface, on the same unwrapped u branch as the rings returned by
// face_loops_uv(); a cut that crosses the periodic seam without being unwrapped
// is an open cut endpoint off the boundary and is rejected. An open cut endpoint
// must lie on the face boundary within `tolerance`, weld to another open cut
// endpoint within that distance, or T-junction onto another cut (its projection
// is inserted into that cut). A closed cut (first point repeated as last) may
// lie anywhere. Cut segments shorter than `tolerance` and fragments whose
// surface area is at most tolerance^2 are dropped. Errors returned by
// face_loops_uv() are propagated unchanged.
inline Result<std::vector<CurvedFragment>> split_face_curved(
    const Brep& brep, FaceId face, std::span<const FaceCut> cuts,
    double tolerance) {
  if (face.is_null() || !brep.face_alive(face)) {
    return Error{Errc::invalid_topology, "split_face_curved: face is not alive",
                 face.index(), 0u};
  }
  if (!std::isfinite(tolerance)) {
    return Error{Errc::non_finite, "split_face_curved: tolerance is not finite",
                 face.index(), 0u};
  }
  if (tolerance < 0.0) {
    return Error{Errc::tolerance_violation,
                 "split_face_curved: tolerance is negative", face.index(), 0u};
  }

  const Result<std::vector<FaceUvLoop>> loops_result =
      face_loops_uv(brep, face, tolerance);
  if (!loops_result) return loops_result.error();
  const std::vector<FaceUvLoop>& uv_loops = loops_result.value();
  if (uv_loops.empty()) {
    return Error{Errc::invalid_topology,
                 "split_face_curved: face has no boundary loops", face.index(),
                 0u};
  }
  if (!uv_loops.front().outer) {
    return Error{Errc::invalid_topology,
                 "split_face_curved: face does not start with an outer loop",
                 face.index(), 0u};
  }
  // Defensive: SolidEditor::commit() already rejects a face with more than one
  // outer loop, so this branch is only reachable for a corrupt Brep that was
  // assembled some other way.
  std::size_t outer_loops = 0;
  for (const FaceUvLoop& loop : uv_loops) {
    if (loop.outer) ++outer_loops;
  }
  if (outer_loops != 1) {
    return Error{Errc::unsupported,
                 "split_face_curved: face must have exactly one outer loop",
                 face.index(), static_cast<std::uint32_t>(outer_loops)};
  }

  const Result<detail::FaceUvProjection> projection_result =
      detail::face_uv_projection(brep, brep.face(face), tolerance);
  if (!projection_result) return projection_result.error();
  const detail::FaceUvProjection& projection = projection_result.value();

  std::vector<std::vector<Point2d>> loops;
  loops.reserve(uv_loops.size());
  double min_u = 0.0;
  double max_u = 0.0;
  double min_v = 0.0;
  double max_v = 0.0;
  bool have_extent = false;
  for (const FaceUvLoop& loop : uv_loops) {
    if (loop.points.empty()) {
      return Error{Errc::invalid_topology,
                   "split_face_curved: face has an empty boundary loop",
                   face.index(), 0u};
    }
    loops.push_back(loop.points);
    for (const Point2d& point : loop.points) {
      if (!detail::trim_face_finite2(point)) {
        return Error{Errc::non_finite,
                     "split_face_curved: boundary point is not finite",
                     face.index(), 0u};
      }
      if (!have_extent) {
        min_u = max_u = point[0];
        min_v = max_v = point[1];
        have_extent = true;
      } else {
        min_u = std::min(min_u, point[0]);
        max_u = std::max(max_u, point[0]);
        min_v = std::min(min_v, point[1]);
        max_v = std::max(max_v, point[1]);
      }
    }
  }
  const double extent = std::hypot(max_u - min_u, max_v - min_v);
  if (!std::isfinite(extent) || !(extent > 0.0)) {
    return Error{Errc::degenerate,
                 "split_face_curved: face boundary has no extent",
                 face.index(), 0u};
  }
  const double snap_radius = std::max(tolerance, 1e-12 * extent);

  // Validate cuts and classify them as open or closed.
  std::vector<std::vector<Point2d>> cut_points;
  std::vector<char> cut_closed;
  cut_points.reserve(cuts.size());
  cut_closed.reserve(cuts.size());
  for (std::size_t c = 0; c < cuts.size(); ++c) {
    const std::vector<Point2d>& uv = cuts[c].uv;
    if (uv.size() < 2) {
      return Error{Errc::invalid_topology,
                   "split_face_curved: a cut has fewer than two points",
                   face.index(), static_cast<std::uint32_t>(c)};
    }
    for (const Point2d& point : uv) {
      if (!detail::trim_face_finite2(point)) {
        return Error{Errc::non_finite,
                     "split_face_curved: cut point is not finite", face.index(),
                     static_cast<std::uint32_t>(c)};
      }
    }
    const double closing = std::hypot(uv.front()[0] - uv.back()[0],
                                      uv.front()[1] - uv.back()[1]);
    const bool closed = closing <= snap_radius;
    std::vector<Point2d> points = uv;
    if (closed) points.back() = points.front();
    cut_points.push_back(std::move(points));
    cut_closed.push_back(closed ? 1 : 0);
  }

  // Snap every open cut endpoint. A welded group of endpoints shares one exact
  // coordinate so touching cuts and the boundary meet exactly. A free endpoint
  // (not on the boundary and not welded to another cut) must T-junction onto
  // another cut; the projection is inserted into that cut so the arrangement
  // sees an exact split.
  struct Endpoint {
    Point2d raw{};
    detail::TrimFaceSnap snap{};
    Point2d final{};
    std::size_t cut = 0;
  };
  std::vector<Endpoint> endpoints;
  for (std::size_t c = 0; c < cuts.size(); ++c) {
    if (cut_closed[c] != 0) continue;
    endpoints.push_back(Endpoint{cut_points[c].front(),
                                 detail::trim_face_find_snap(
                                     loops, cut_points[c].front(), snap_radius),
                                 Point2d{}, c});
    endpoints.push_back(Endpoint{cut_points[c].back(),
                                 detail::trim_face_find_snap(
                                     loops, cut_points[c].back(), snap_radius),
                                 Point2d{}, c});
  }
  const std::size_t endpoint_count = endpoints.size();
  std::vector<std::size_t> parent(endpoint_count);
  for (std::size_t i = 0; i < endpoint_count; ++i) parent[i] = i;
  for (std::size_t i = 0; i < endpoint_count; ++i) {
    for (std::size_t j = i + 1; j < endpoint_count; ++j) {
      const double distance = std::hypot(
          endpoints[i].raw[0] - endpoints[j].raw[0],
          endpoints[i].raw[1] - endpoints[j].raw[1]);
      if (!(distance <= snap_radius)) continue;
      const std::size_t root_i = detail::trim_face_find_root(parent, i);
      const std::size_t root_j = detail::trim_face_find_root(parent, j);
      if (root_i != root_j) parent[root_j] = root_i;
    }
  }
  std::vector<detail::TrimFaceSnap> canonical(endpoint_count);
  std::vector<char> has_canonical(endpoint_count, 0);
  for (std::size_t i = 0; i < endpoint_count; ++i) {
    const std::size_t root = detail::trim_face_find_root(parent, i);
    if (has_canonical[root] == 0 ||
        detail::trim_face_snap_better(endpoints[i].snap, canonical[root])) {
      canonical[root] = endpoints[i].snap;
      has_canonical[root] = 1;
    }
  }

  struct CutInsert {
    std::size_t cut = 0;
    std::size_t edge = 0;
    double t = 0.0;
    Point2d point{};
  };
  std::vector<CutInsert> cut_inserts;
  std::vector<char> group_done(endpoint_count, 0);
  for (std::size_t i = 0; i < endpoint_count; ++i) {
    const std::size_t root = detail::trim_face_find_root(parent, i);
    if (group_done[root] != 0) continue;
    group_done[root] = 1;
    const detail::TrimFaceSnap& snap = canonical[root];
    Point2d final = snap.point;
    if (!snap.boundary) {
      // How many different cut polylines meet at this welded point?
      std::size_t distinct_cuts = 0;
      std::size_t sole_cut = 0;
      for (std::size_t j = 0; j < endpoint_count; ++j) {
        if (detail::trim_face_find_root(parent, j) != root) continue;
        if (distinct_cuts == 0) {
          distinct_cuts = 1;
          sole_cut = endpoints[j].cut;
        } else if (endpoints[j].cut != sole_cut) {
          distinct_cuts = 2;
        }
      }
      if (distinct_cuts < 2) {
        // A single free endpoint is a T-junction candidate: it must land in the
        // interior of a segment of another cut. Two welded endpoints of the SAME
        // cut cannot reach this branch: a cut whose first and last point are
        // within snap_radius was classified as a closed loop before endpoint
        // grouping, so this group holds exactly one endpoint.
        bool have_hit = false;
        double best_distance = 0.0;
        CutInsert best{};
        for (std::size_t host = 0; host < cuts.size(); ++host) {
          if (host == sole_cut) continue;
          const std::vector<Point2d>& chain = cut_points[host];
          for (std::size_t e = 0; e + 1 < chain.size(); ++e) {
            const Point2d& a = chain[e];
            const Point2d& b = chain[e + 1];
            const Vec2d edge = b - a;
            const double length2 = dot(edge, edge);
            if (!(length2 > 0.0)) continue;
            const double t = dot(snap.point - a, edge) / length2;
            if (!(t > 0.0) || !(t < 1.0)) continue;
            const Point2d hit_point{a[0] + edge[0] * t,
                                    a[1] + edge[1] * t};
            const double distance =
                std::hypot(snap.point[0] - hit_point[0],
                           snap.point[1] - hit_point[1]);
            if (!(distance <= snap_radius)) continue;
            if (have_hit && !(distance < best_distance)) continue;
            have_hit = true;
            best_distance = distance;
            best = CutInsert{host, e, t, hit_point};
          }
        }
        if (!have_hit) {
          return Error{Errc::numeric_failure,
                       "split_face_curved: an open cut endpoint is not on the "
                       "face boundary",
                       face.index(), static_cast<std::uint32_t>(endpoint_count)};
        }
        final = best.point;
        cut_inserts.push_back(best);
      }
    }
    for (std::size_t j = 0; j < endpoint_count; ++j) {
      if (detail::trim_face_find_root(parent, j) == root) {
        endpoints[j].final = final;
      }
    }
  }

  // Insert T-junction projections into the host cut polylines at their original
  // segment positions. The endpoint and the host then share one exact vertex.
  if (!cut_inserts.empty()) {
    std::vector<std::vector<CutInsert>> per_cut(cuts.size());
    for (const CutInsert& insert : cut_inserts) {
      per_cut[insert.cut].push_back(insert);
    }
    for (std::size_t c = 0; c < cuts.size(); ++c) {
      if (per_cut[c].empty()) continue;
      std::stable_sort(per_cut[c].begin(), per_cut[c].end(),
                       [](const CutInsert& lhs, const CutInsert& rhs) noexcept {
                         if (lhs.edge != rhs.edge) return lhs.edge < rhs.edge;
                         return lhs.t < rhs.t;
                       });
      const std::vector<Point2d>& chain = cut_points[c];
      std::vector<Point2d> rebuilt;
      rebuilt.reserve(chain.size() + per_cut[c].size());
      std::size_t next = 0;
      for (std::size_t e = 0; e + 1 < chain.size(); ++e) {
        rebuilt.push_back(chain[e]);
        while (next < per_cut[c].size() && per_cut[c][next].edge == e) {
          const Point2d& point = per_cut[c][next].point;
          bool present = false;
          for (const Point2d& vertex : chain) {
            if (detail::trim_face_same_xy(vertex, point)) {
              present = true;
              break;
            }
          }
          if (!present && !detail::trim_face_same_xy(rebuilt.back(), point)) {
            rebuilt.push_back(point);
          }
          ++next;
        }
      }
      rebuilt.push_back(chain.back());
      cut_points[c] = std::move(rebuilt);
    }
  }

  // Insert segment snaps into the boundary rings so the cut endpoint and the
  // ring share one exact vertex.
  std::vector<std::vector<Point2d>> rings = loops;
  std::vector<detail::TrimFaceInsert> inserts;
  std::vector<char> inserted(endpoint_count, 0);
  for (std::size_t i = 0; i < endpoint_count; ++i) {
    const std::size_t root = detail::trim_face_find_root(parent, i);
    if (inserted[root] != 0) continue;
    inserted[root] = 1;
    const detail::TrimFaceSnap& snap = canonical[root];
    if (!snap.boundary || snap.vertex) continue;
    bool already_present = false;
    for (const Point2d& vertex : rings[snap.ring]) {
      if (detail::trim_face_same_xy(vertex, snap.point)) {
        already_present = true;
        break;
      }
    }
    if (!already_present) {
      inserts.push_back({snap.ring, snap.edge, snap.point});
    }
  }
  detail::trim_face_apply_inserts(rings, std::move(inserts));

  // Rewrite the open cut chains with the welded endpoint coordinates.
  std::size_t endpoint_index = 0;
  for (std::size_t c = 0; c < cuts.size(); ++c) {
    if (cut_closed[c] != 0) continue;
    cut_points[c].front() = endpoints[endpoint_index].final;
    cut_points[c].back() = endpoints[endpoint_index + 1].final;
    endpoint_index += 2;
  }

  // A sub-tolerance FIRST or LAST segment of an open chain is a snapping
  // artifact: the rewritten endpoint sits within the tolerance of its
  // neighbour. Collapse that neighbour onto the endpoint. Interior segments are
  // never dropped, however short: dropping one would silently disconnect the
  // chain into two open pieces, and the arrangement would see two free interior
  // endpoints instead of a cut crossing the face.
  for (std::size_t c = 0; c < cuts.size(); ++c) {
    if (cut_closed[c] != 0) continue;
    std::vector<Point2d>& chain = cut_points[c];
    while (chain.size() > 2 &&
           !(std::hypot(chain[1][0] - chain[0][0], chain[1][1] - chain[0][1]) >
             tolerance)) {
      chain.erase(chain.begin() + 1);
    }
    while (chain.size() > 2 &&
           !(std::hypot(chain[chain.size() - 2][0] - chain.back()[0],
                        chain[chain.size() - 2][1] - chain.back()[1]) >
             tolerance)) {
      chain.erase(chain.end() - 2);
    }
  }

  // All boundary loops belong to one source, so their contributions cancel
  // inside a hole: the outer ring winds once and the hole ring winds once the
  // other way. Each cut is its own source.
  std::vector<ArrangementSegment> segments;
  std::size_t segment_total = 0;
  for (const std::vector<Point2d>& ring : rings) segment_total += ring.size();
  segments.reserve(segment_total + cuts.size() * 8u);
  for (std::size_t r = 0; r < rings.size(); ++r) {
    const std::vector<Point2d>& ring = rings[r];
    for (std::size_t i = 0; i < ring.size(); ++i) {
      const Point2d& a = ring[i];
      const Point2d& b = ring[(i + 1) % ring.size()];
      if (detail::trim_face_same_xy(a, b)) continue;
      segments.push_back(ArrangementSegment{a, b, 0u});
    }
  }
  for (std::size_t c = 0; c < cuts.size(); ++c) {
    const std::vector<Point2d>& chain = cut_points[c];
    std::size_t used = 0;
    for (std::size_t i = 0; i + 1 < chain.size(); ++i) {
      const Point2d& a = chain[i];
      const Point2d& b = chain[i + 1];
      // Sub-tolerance first/last segments were collapsed above; interior short
      // segments are kept so the chain stays connected. Exact duplicate points
      // still carry no direction and are skipped.
      if (detail::trim_face_same_xy(a, b)) {
        continue;
      }
      ++used;
      segments.push_back(ArrangementSegment{a, b,
                                            static_cast<std::uint32_t>(1u + c)});
    }
    if (used == 0) {
      return Error{Errc::degenerate,
                   "split_face_curved: cut collapses to a point", face.index(),
                   static_cast<std::uint32_t>(c)};
    }
  }

  const detail::SplitResult split = detail::split_segments(segments);
  if (!std::isfinite(split.max_residual)) {
    return Error{Errc::numeric_failure,
                 "split_face_curved: segment intersection is not representable",
                 face.index(), 0u};
  }
  const Arrangement2 arrangement = Arrangement2::build(segments);

  // Every bounded cycle is a candidate fragment boundary. A cycle is material
  // when the boundary source winds around it; hole rings contribute the opposite
  // winding, so a point inside a hole has winding zero and an island inside a
  // hole has winding one again.
  std::vector<detail::TrimFaceCycle> cycles;
  for (std::uint32_t index = 0; index < arrangement.face_count(); ++index) {
    if (!arrangement.faces()[index].bounded) continue;
    detail::TrimFaceCycle record;
    record.face = index;
    if (!detail::trim_face_walk_cycle(arrangement, index, record.ring)) {
      return Error{Errc::invalid_topology,
                   "split_face_curved: arrangement face cycle is not closed",
                   face.index(), index};
    }
    detail::trim_face_drop_repeats(record.ring);
    if (record.ring.size() < 3) continue;
    record.area = std::fabs(signed_area2(record.ring));
    if (!(record.area > 0.0)) continue;

    record.material = arrangement.winding(index, 0u) != 0;

    Polygon2d polygon;
    polygon.vertices = record.ring;
    polygon.ring_starts.push_back(0u);
    polygon.ring_starts.push_back(
        static_cast<std::uint32_t>(record.ring.size()));
    const Triangulation2 triangulation = triangulate_polygon2(polygon);
    if (triangulation.ok && !triangulation.triangles.empty()) {
      const std::array<std::uint32_t, 3>& triangle =
          triangulation.triangles.front();
      record.rep = Point2d{
          (polygon.vertices[triangle[0]][0] + polygon.vertices[triangle[1]][0] +
           polygon.vertices[triangle[2]][0]) /
              3.0,
          (polygon.vertices[triangle[0]][1] + polygon.vertices[triangle[1]][1] +
           polygon.vertices[triangle[2]][1]) /
              3.0};
      record.has_rep = true;
    }
    cycles.push_back(std::move(record));
  }

  // Parent of a cycle is the smallest-area larger cycle that contains its
  // representative point. Immediate children are the fragment's holes.
  for (std::size_t index = 0; index < cycles.size(); ++index) {
    if (!cycles[index].has_rep) continue;
    std::size_t best = static_cast<std::size_t>(-1);
    for (std::size_t candidate = 0; candidate < cycles.size(); ++candidate) {
      if (candidate == index) continue;
      if (!(cycles[candidate].area > cycles[index].area)) continue;
      if (!detail::face_query_point_in_ring(cycles[candidate].ring,
                                            cycles[index].rep)) {
        continue;
      }
      if (best == static_cast<std::size_t>(-1) ||
          cycles[candidate].area < cycles[best].area) {
        best = candidate;
      }
    }
    cycles[index].parent = best;
  }
  for (std::size_t index = 0; index < cycles.size(); ++index) {
    if (cycles[index].parent != static_cast<std::size_t>(-1)) {
      cycles[cycles[index].parent].children.push_back(index);
    }
  }

  std::vector<CurvedFragment> fragments;
  const double area_floor = tolerance * tolerance;
  for (const detail::TrimFaceCycle& record : cycles) {
    if (!record.material || !record.has_rep) continue;

    CurvedFragment fragment;
    fragment.reversed = false;
    std::vector<Point2d> outer = record.ring;
    detail::trim_face_normalize_ring(outer, true);
    if (outer.size() < 3) continue;
    const double outer_area = std::fabs(signed_area2(outer));
    if (!(outer_area > 0.0)) continue;
    double area_uv = outer_area;
    std::vector<std::vector<Point2d>> holes_uv;
    holes_uv.reserve(record.children.size());
    for (const std::size_t child : record.children) {
      std::vector<Point2d> hole = cycles[child].ring;
      detail::trim_face_normalize_ring(hole, false);
      if (hole.size() < 3) continue;
      const double hole_area = std::fabs(signed_area2(hole));
      if (!(hole_area > 0.0)) continue;
      area_uv -= hole_area;
      holes_uv.push_back(std::move(hole));
    }
    // Holes are ordered by area (descending) and then lexicographically.
    std::sort(holes_uv.begin(), holes_uv.end(),
              [](const std::vector<Point2d>& lhs,
                 const std::vector<Point2d>& rhs) noexcept {
                const double lhs_area = std::fabs(signed_area2(lhs));
                const double rhs_area = std::fabs(signed_area2(rhs));
                if (lhs_area != rhs_area) return lhs_area > rhs_area;
                return detail::trim_face_lex_less(lhs.front(), rhs.front());
              });
    if (!(area_uv > 0.0)) continue;

    // A fragment suppressed by the area floor (for example a tiny island) does
    // not remove its cycle from the arrangement: the enclosing fragment keeps
    // the island outline as a hole, so the emitted parent ring has a hole where
    // no fragment was emitted. The hole is left in place rather than patched,
    // because its outline is still an exact piece of the input subdivision.
    const double scale = detail::trim_face_metric_scale(projection, record.rep);
    if (!(area_uv * scale > area_floor)) continue;

    fragment.loops_uv.push_back(std::move(outer));
    for (std::vector<Point2d>& hole : holes_uv) {
      fragment.loops_uv.push_back(std::move(hole));
    }

    bool duplicate = false;
    for (const CurvedFragment& existing : fragments) {
      if (detail::trim_face_same_fragment(existing, fragment)) {
        duplicate = true;
        break;
      }
    }
    if (duplicate) continue;
    fragments.push_back(std::move(fragment));
  }

  std::stable_sort(fragments.begin(), fragments.end(),
                   [](const CurvedFragment& lhs,
                      const CurvedFragment& rhs) noexcept {
                     const Point2d& la = lhs.loops_uv.front().front();
                     const Point2d& ra = rhs.loops_uv.front().front();
                     if (la[0] != ra[0]) return la[0] < ra[0];
                     if (la[1] != ra[1]) return la[1] < ra[1];
                     return lhs.loops_uv.size() < rhs.loops_uv.size();
                   });

  // Lift every ring into 3D only after the fragment set is final.
  for (CurvedFragment& fragment : fragments) {
    fragment.loops.reserve(fragment.loops_uv.size());
    for (const std::vector<Point2d>& ring : fragment.loops_uv) {
      std::vector<Point3d> ring3;
      ring3.reserve(ring.size());
      for (const Point2d& uv : ring) {
        const Point3d point = detail::face_uv_eval(projection, uv);
        if (!detail::trim_face_finite3(point)) {
          return Error{Errc::non_finite,
                       "split_face_curved: surface evaluation is not finite",
                       face.index(), 0u};
        }
        ring3.push_back(point);
      }
      fragment.loops.push_back(std::move(ring3));
    }
  }
  return fragments;
}

}  // namespace v1
}  // namespace geom
