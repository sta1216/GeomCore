#pragma once

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <optional>
#include <span>
#include <utility>
#include <vector>

#include <geom/brep/face_query.hpp>
#include <geom/brep/face_uv.hpp>
#include <geom/brep/trim_domain.hpp>
#include <geom/cgeom/arrangement2.hpp>
#include <geom/cgeom/polygon2.hpp>
#include <geom/cgeom/triangulate2.hpp>
#include <geom/core/error.hpp>
#include <geom/core/result.hpp>
#include <geom/primitives/point.hpp>
#include <geom/primitives/vec.hpp>

namespace geom {
inline namespace v1 {

// One piece of a boundary arc carried by a fragment ring: the circle and the
// parameter interval the ring follows. t0 < t1 and t1 - t0 <= 2*pi; a segment
// that closes a full circle has t1 - t0 == 2*pi. The interval is the geometric
// extent of the fragment's signed arc reference, which records the traversal
// direction, so the ring may follow the interval in either direction.
struct CurvedArcSpan {
  CircleCurve3 circle;
  double t0 = 0.0;
  double t1 = 0.0;
};

// One piece of a split curved face. `loops` holds the 3D rings and `loops_uv`
// the matching parameter rings; index 0 is the outer ring and the rest are
// holes. The rings are implicitly closed, outer is counter-clockwise and holes
// are clockwise. `reversed` is always false here: the boolean driver flips the
// B-side fragments when it needs the opposite material side.
// On a seam-only face (a sphere or torus face whose loop is the seam slit) the
// loops contain only cut geometry, and a seam-closed ring is a closed chain
// whose first point is repeated as last: the repeated first point and the point
// before it are the two parameter copies of one physical seam point, identified
// with the same 1e-12 relative bound the parameter projection uses. Such a ring
// is ordered first as the outer boundary (before any interior closed cut, which
// is a hole) because its uv area is a chord artefact, not a region area.
// `arc_spans[l][i]` describes the segment from point i to point (i+1)%size of
// loop l. It is empty by default and, after split_face_curved, is sized exactly
// like `loops`, with a value only on segments that follow an arc.
struct CurvedFragment {
  std::vector<std::vector<Point3d>> loops;
  std::vector<std::vector<Point2d>> loops_uv;
  bool reversed = false;
  std::vector<std::vector<std::optional<CurvedArcSpan>>> arc_spans;
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

// Inserts snapped boundary points into a ring and splits the arc reference of
// every segment a point lands on. An endpoint snapped onto an arc must become a
// loop point of its circle: the fragment mapping matches ring vertices against
// the loop points exactly, and a ring vertex left inside an arc segment is a
// hard error, while the arc spans are read per segment, so the segment a point
// splits has to carry two exact references.
inline void trim_face_apply_inserts_with_arcs(
    const FaceUvProjection& projection, std::vector<Point2d>& ring,
    std::vector<std::optional<FaceUvArcRef>>& arcs,
    std::vector<TrimFaceInsert> inserts) {
  if (inserts.empty()) return;
  std::sort(inserts.begin(), inserts.end(),
            [](const TrimFaceInsert& lhs, const TrimFaceInsert& rhs) noexcept {
              if (lhs.edge != rhs.edge) return lhs.edge < rhs.edge;
              return trim_face_lex_less(lhs.point, rhs.point);
            });
  inserts.erase(std::unique(inserts.begin(), inserts.end(),
                            [](const TrimFaceInsert& lhs,
                               const TrimFaceInsert& rhs) noexcept {
                              return lhs.edge == rhs.edge &&
                                     trim_face_same_xy(lhs.point, rhs.point);
                            }),
                inserts.end());

  const std::size_t count = ring.size();
  std::vector<std::vector<std::pair<double, Point2d>>> extra(count);
  for (const TrimFaceInsert& insert : inserts) {
    if (insert.edge >= count) continue;
    const Point2d& a = ring[insert.edge];
    const Point2d& b = ring[(insert.edge + 1) % count];
    const Vec2d edge_vector = b - a;
    const double length2 = dot(edge_vector, edge_vector);
    const double t =
        length2 > 0.0 ? dot(insert.point - a, edge_vector) / length2 : 0.0;
    extra[insert.edge].emplace_back(t, insert.point);
  }
  for (std::vector<std::pair<double, Point2d>>& edge_extra : extra) {
    std::sort(edge_extra.begin(), edge_extra.end(),
              [](const std::pair<double, Point2d>& lhs,
                 const std::pair<double, Point2d>& rhs) noexcept {
                if (lhs.first != rhs.first) return lhs.first < rhs.first;
                return trim_face_lex_less(lhs.second, rhs.second);
              });
  }

  // Splits the reference of the segment starting at ring[i] at the arc angle of
  // an inserted point. The point was projected onto this circle, so its angle
  // lies within the segment's sweep; when it does not, the segment keeps its
  // reference and the inserted point stays an un-referenced vertex.
  const auto split_angle = [&projection](const FaceUvArcRef& ref,
                                         const Point2d& point,
                                         double& angle) {
    const Result<Point3d> lifted = face_uv_eval(projection, point);
    if (!lifted.has_value()) return false;
    const double sweep = ref.t1 - ref.t0;
    if (!(std::fabs(sweep) > 0.0)) return false;
    double offset = face_uv_circle_angle(ref.circle, lifted.value()) - ref.t0;
    if (offset > kFaceUvPi) {
      offset -= 2.0 * kFaceUvPi;
    } else if (offset < -kFaceUvPi) {
      offset += 2.0 * kFaceUvPi;
    }
    angle = ref.t0 + offset;
    if (sweep > 0.0) return angle > ref.t0 && angle < ref.t1;
    return angle < ref.t0 && angle > ref.t1;
  };

  std::vector<Point2d> rebuilt;
  std::vector<std::optional<FaceUvArcRef>> rebuilt_arcs;
  rebuilt.reserve(count + inserts.size());
  rebuilt_arcs.reserve(count + inserts.size());
  for (std::size_t i = 0; i < count; ++i) {
    const FaceUvArcRef* source =
        i < arcs.size() && arcs[i].has_value() ? &*arcs[i] : nullptr;
    double start_t = source != nullptr ? source->t0 : 0.0;
    rebuilt.push_back(ring[i]);
    for (const std::pair<double, Point2d>& entry : extra[i]) {
      double angle = start_t;
      if (source != nullptr && split_angle(*source, entry.second, angle)) {
        rebuilt_arcs.push_back(FaceUvArcRef{source->circle, start_t, angle});
        start_t = angle;
      } else {
        rebuilt_arcs.push_back(std::nullopt);
      }
      rebuilt.push_back(entry.second);
    }
    if (source) {
      rebuilt_arcs.push_back(FaceUvArcRef{source->circle, start_t, source->t1});
    } else {
      rebuilt_arcs.push_back(std::nullopt);
    }
  }
  ring = std::move(rebuilt);
  arcs = std::move(rebuilt_arcs);
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

inline constexpr double kTrimFaceTwoPi = 2.0 * kFaceUvPi;

// A point where a cut segment meets a boundary circle. `s` is the fraction along
// the cut segment (0 or 1 when the hit is a cut endpoint).
struct TrimFaceCircleHit {
  double s = 0.0;
  Point3d point{};
  double t = 0.0;
};

// Intersects the segment p0->p1 with `circle`: the unique crossing of the
// circle's plane, or the two intersections inside that plane. A segment that
// only grazes the circle (no penetration deeper than tolerance) reports no
// interior hit, so a tangent cut never splits the boundary.
inline std::vector<TrimFaceCircleHit> trim_face_circle_hits(
    const CircleCurve3& circle, const Point3d& p0, const Point3d& p1,
    double tolerance) {
  std::vector<TrimFaceCircleHit> hits;
  const double radius = circle.radius;
  const double tol = std::max(tolerance, 1e-12 * std::max(1.0, radius));
  const Vec3d normal = cross(circle.x_axis, circle.y_axis);
  const Vec3d d = p1 - p0;
  const Vec3d f = p0 - circle.center;
  const auto on_circle = [&](const Point3d& point) {
    const Vec3d offset = point - circle.center;
    return std::fabs(norm2(offset) - radius) <= tol &&
           std::fabs(dot(offset, normal)) <= tol;
  };
  if (on_circle(p0)) {
    hits.push_back(
        TrimFaceCircleHit{0.0, p0, face_uv_circle_angle(circle, p0)});
  }
  if (on_circle(p1)) {
    hits.push_back(
        TrimFaceCircleHit{1.0, p1, face_uv_circle_angle(circle, p1)});
  }
  const double along_normal = dot(d, normal);
  if (std::fabs(along_normal) <= tol) {
    // The segment lies in (or parallel to) the circle's plane.
    if (std::fabs(dot(f, normal)) > tol) return hits;
    const double start_u = dot(f, circle.x_axis);
    const double start_v = dot(f, circle.y_axis);
    const double delta_u = dot(d, circle.x_axis);
    const double delta_v = dot(d, circle.y_axis);
    const double a = delta_u * delta_u + delta_v * delta_v;
    if (!(a > 0.0)) return hits;
    const double b = 2.0 * (start_u * delta_u + start_v * delta_v);
    const double c = start_u * start_u + start_v * start_v - radius * radius;
    const double discriminant = b * b - 4.0 * a * c;
    if (!(discriminant > 0.0)) return hits;
    // Reject a tangency: the chord must sink below the radius by more than the
    // tolerance at its closest point.
    const double middle_s = -b / (2.0 * a);
    const double middle_u = start_u + delta_u * middle_s;
    const double middle_v = start_v + delta_v * middle_s;
    const double middle_distance =
        std::hypot(middle_u, middle_v);
    if (!(radius - middle_distance > tol)) return hits;
    const double root = std::sqrt(discriminant);
    const double params[2] = {(-b - root) / (2.0 * a),
                              (-b + root) / (2.0 * a)};
    for (const double s : params) {
      if (!(s > 0.0) || !(s < 1.0)) continue;
      const Point3d point{p0[0] + d[0] * s, p0[1] + d[1] * s,
                          p0[2] + d[2] * s};
      hits.push_back(
          TrimFaceCircleHit{s, point, face_uv_circle_angle(circle, point)});
    }
    return hits;
  }
  const double s = -dot(f, normal) / along_normal;
  if (!(s > 0.0) || !(s < 1.0)) return hits;
  const Point3d point{p0[0] + d[0] * s, p0[1] + d[1] * s, p0[2] + d[2] * s};
  if (std::fabs(norm2(point - circle.center) - radius) > tol) return hits;
  hits.push_back(
      TrimFaceCircleHit{s, point, face_uv_circle_angle(circle, point)});
  return hits;
}

// The loop point a fragment ring vertex came from.
using TrimFaceRingRef = std::pair<std::size_t, std::size_t>;

inline void trim_face_drop_repeats_refs(
    std::vector<Point2d>& ring,
    std::vector<std::optional<TrimFaceRingRef>>& refs) {
  if (refs.size() != ring.size()) refs.resize(ring.size());
  std::vector<Point2d> cleaned;
  std::vector<std::optional<TrimFaceRingRef>> cleaned_refs;
  cleaned.reserve(ring.size());
  cleaned_refs.reserve(ring.size());
  for (std::size_t i = 0; i < ring.size(); ++i) {
    if (!cleaned.empty() && trim_face_same_xy(cleaned.back(), ring[i])) continue;
    cleaned.push_back(ring[i]);
    cleaned_refs.push_back(refs[i]);
  }
  while (cleaned.size() >= 2 &&
         trim_face_same_xy(cleaned.front(), cleaned.back())) {
    cleaned.pop_back();
    cleaned_refs.pop_back();
  }
  ring.swap(cleaned);
  refs.swap(cleaned_refs);
}

// trim_face_normalize_ring with the per-vertex loop references kept in lockstep.
inline void trim_face_normalize_ring_refs(
    std::vector<Point2d>& ring,
    std::vector<std::optional<TrimFaceRingRef>>& refs, bool want_ccw) {
  trim_face_drop_repeats_refs(ring, refs);
  if (ring.size() < 3) return;
  const int orientation = ring_orientation2(ring);
  if (orientation != 0 && ((orientation > 0) != want_ccw)) {
    std::reverse(ring.begin(), ring.end());
    std::reverse(refs.begin(), refs.end());
  }
  std::size_t first = 0;
  for (std::size_t i = 1; i < ring.size(); ++i) {
    if (trim_face_lex_less(ring[i], ring[first])) first = i;
  }
  if (first != 0) {
    const auto offset = static_cast<std::ptrdiff_t>(first);
    std::rotate(ring.begin(), ring.begin() + offset, ring.end());
    std::rotate(refs.begin(), refs.begin() + offset, refs.end());
  }
}

// The arc span of ring segment i, when that segment is a boundary arc. The
// fragment's signed arc reference names the segment (refs[i] is the loop index
// the segment starts at) and the reference's own direction gives the interval.
inline std::optional<CurvedArcSpan> trim_face_ring_span(
    const std::vector<std::optional<TrimFaceRingRef>>& refs,
    const std::vector<std::vector<Point2d>>& loops,
    const std::vector<std::vector<std::optional<FaceUvArcRef>>>& loop_arcs,
    std::size_t i) {
  if (refs.empty()) return std::nullopt;
  const std::size_t j = (i + 1) % refs.size();
  if (!refs[i] || !refs[j]) return std::nullopt;
  const TrimFaceRingRef a = *refs[i];
  const TrimFaceRingRef b = *refs[j];
  if (a.first != b.first) return std::nullopt;
  if (a.first >= loops.size()) return std::nullopt;
  const std::size_t count = loops[a.first].size();
  if (count < 2) return std::nullopt;
  if (!((a.second + 1) % count == b.second ||
        (b.second + 1) % count == a.second)) {
    return std::nullopt;
  }
  // The ring may run either way around the loop, so the arc reference that
  // describes this edge starts at whichever endpoint comes first in loop order.
  const std::size_t start_index =
      (b.second + 1) % count == a.second ? b.second : a.second;
  if (start_index >= loop_arcs[a.first].size()) return std::nullopt;
  const std::optional<FaceUvArcRef>& arc = loop_arcs[a.first][start_index];
  if (!arc) return std::nullopt;
  const double sweep = arc->t1 - arc->t0;
  if (!(std::fabs(sweep) > 0.0)) return std::nullopt;
  const double low = std::min(arc->t0, arc->t1);
  const double high = std::max(arc->t0, arc->t1);
  const double shifted =
      low - kTrimFaceTwoPi * std::floor(low / kTrimFaceTwoPi);
  return CurvedArcSpan{arc->circle, shifted, shifted + (high - low)};
}

// True when the uv point lies on the face boundary within `tolerance`.
inline bool trim_face_uv_on_ring(const std::vector<std::vector<Point2d>>& loops,
                                 const Point2d& point, double tolerance) noexcept {
  for (const std::vector<Point2d>& ring : loops) {
    const std::size_t count = ring.size();
    if (count < 2) continue;
    for (std::size_t i = 0; i < count; ++i) {
      if (face_query_point_segment_distance(point, ring[i], ring[(i + 1) % count]) <= tolerance) {
        return true;
      }
    }
  }
  return false;
}

// True when the uv point is on the face: on a ring counts as inside, otherwise
// the ring parity decides (a hole ring flips it back to outside).
inline bool trim_face_uv_inside(const std::vector<std::vector<Point2d>>& loops,
                                const Point2d& point, double tolerance) noexcept {
  if (trim_face_uv_on_ring(loops, point, tolerance)) return true;
  bool inside = false;
  for (const std::vector<Point2d>& ring : loops) {
    if (face_query_point_in_ring(ring, point)) inside = !inside;
  }
  return inside;
}

// Closest point of the boundary to `point`.
inline Point2d trim_face_nearest_ring_point(const std::vector<std::vector<Point2d>>& loops,
                                            const Point2d& point) noexcept {
  Point2d best = point;
  double best_distance = std::numeric_limits<double>::infinity();
  for (const std::vector<Point2d>& ring : loops) {
    const std::size_t count = ring.size();
    if (count < 2) continue;
    for (std::size_t i = 0; i < count; ++i) {
      const Point2d& a = ring[i];
      const Point2d& b = ring[(i + 1) % count];
      const double distance = face_query_point_segment_distance(point, a, b);
      if (!(distance < best_distance)) continue;
      best_distance = distance;
      best = trim_face_project_on_edge(point, a, b);
    }
  }
  return best;
}

// The point where the segment from an inside point to an outside point crosses
// the boundary. Bisection converges on the ring, then the crossing is projected
// onto the crossed ring segment so the caller gets an exact boundary point.
inline Point2d trim_face_boundary_crossing(const std::vector<std::vector<Point2d>>& loops,
                                           const Point2d& inside, const Point2d& outside,
                                           double tolerance) noexcept {
  double low = 0.0;
  double high = 1.0;
  for (int step = 0; step < 60; ++step) {
    const double mid = 0.5 * (low + high);
    const Point2d point{inside[0] + mid * (outside[0] - inside[0]),
                        inside[1] + mid * (outside[1] - inside[1])};
    if (trim_face_uv_inside(loops, point, tolerance)) {
      low = mid;
    } else {
      high = mid;
    }
  }
  const Point2d point{inside[0] + low * (outside[0] - inside[0]),
                      inside[1] + low * (outside[1] - inside[1])};
  return trim_face_nearest_ring_point(loops, point);
}

// Clips one cut chain to the face. The intersector traces a curve over the whole
// parameter patch, which for a trimmed face (a disk cap, say) reaches past the
// real boundary. The arc junctions inserted before this point are exact boundary
// points, so a chain that leaves the face and comes back keeps those exact
// crossings and only the outside pieces are dropped here.
inline std::vector<std::vector<Point2d>> trim_face_clip_cut(
    const std::vector<std::vector<Point2d>>& loops, const std::vector<Point2d>& chain, bool closed,
    double tolerance) {
  const std::size_t count = chain.size();
  if (count < 2) return {};
  std::vector<char> flags(count, 0);
  std::size_t inside_count = 0;
  for (std::size_t i = 0; i < count; ++i) {
    flags[i] = trim_face_uv_inside(loops, chain[i], tolerance) ? 1 : 0;
    inside_count += flags[i] != 0 ? 1u : 0u;
  }
  if (inside_count == count) return {chain};
  if (inside_count == 0) return {};

  std::vector<Point2d> points = chain;
  std::vector<char> state = flags;
  if (closed) {
    // Start at an outside point so an inside run never wraps the array end.
    std::size_t start = 0;
    while (start < points.size() && state[start] != 0) ++start;
    if (start == points.size()) return {};
    const auto offset = static_cast<std::ptrdiff_t>(start);
    std::rotate(points.begin(), points.begin() + offset, points.end());
    std::rotate(state.begin(), state.begin() + offset, state.end());
  }

  std::vector<std::vector<Point2d>> runs;
  std::vector<Point2d> current;
  for (std::size_t i = 0; i < points.size(); ++i) {
    if (state[i] != 0) {
      if (current.empty() && i > 0) {
        current.push_back(trim_face_boundary_crossing(loops, points[i], points[i - 1], tolerance));
      }
      current.push_back(points[i]);
      continue;
    }
    if (current.empty()) continue;
    current.push_back(trim_face_boundary_crossing(loops, points[i - 1], points[i], tolerance));
    runs.push_back(std::move(current));
    current.clear();
  }
  if (!current.empty()) runs.push_back(std::move(current));

  std::vector<std::vector<Point2d>> clipped;
  clipped.reserve(runs.size());
  for (std::vector<Point2d>& run : runs) {
    std::vector<Point2d> piece;
    piece.reserve(run.size());
    for (const Point2d& point : run) {
      if (!piece.empty()) {
        const double du = point[0] - piece.back()[0];
        const double dv = point[1] - piece.back()[1];
        if (du * du + dv * dv <= tolerance * tolerance) continue;
      }
      piece.push_back(point);
    }
    if (piece.size() >= 2) clipped.push_back(std::move(piece));
  }
  return clipped;
}

}  // namespace detail

namespace detail {

// The full parameter rectangle of a seam-only face. `boolean_solid_curved.hpp`
// declares the same domains in curved_boolean_surface_uv_domain(), but that
// header includes this one, so the splitter cannot call it without an include
// cycle. Only Sphere3 and Torus3 are seam-only in this task; Cone is Task 7.
inline std::optional<std::array<double, 4>> trim_face_seam_only_domain(
    const SurfaceGeometry& surface) noexcept {
  if (std::holds_alternative<Sphere3>(surface)) {
    return std::array<double, 4>{0.0, 2.0 * kFaceUvPi, -0.5 * kFaceUvPi,
                                 0.5 * kFaceUvPi};
  }
  if (std::holds_alternative<Torus3>(surface)) {
    return std::array<double, 4>{0.0, 2.0 * kFaceUvPi, 0.0, 2.0 * kFaceUvPi};
  }
  return std::nullopt;
}

// True when the value lies within the standard parameter tolerance of one
// domain side value. The tolerance is the same relative 1e-12 the projection
// uses for parameter equality; a single point may touch a side (an open cut
// endpoint on the boundary), but a whole segment may not.
inline bool trim_face_seam_only_near_side(double value, double side) noexcept {
  return std::fabs(value - side) <= 1e-12 * std::max(1.0, std::fabs(value));
}

// True when the segment lies along one domain side. Such a cut cannot be told
// apart from the synthetic chain on that side after the arrangement, and the
// strip step would delete it together with the synthetic geometry, so the
// caller rejects the cut up front.
inline bool trim_face_seam_only_on_same_side(
    const Point2d& a, const Point2d& b,
    const std::array<double, 4>& domain) noexcept {
  return (trim_face_seam_only_near_side(a[0], domain[0]) &&
          trim_face_seam_only_near_side(b[0], domain[0])) ||
         (trim_face_seam_only_near_side(a[0], domain[1]) &&
          trim_face_seam_only_near_side(b[0], domain[1])) ||
         (trim_face_seam_only_near_side(a[1], domain[2]) &&
          trim_face_seam_only_near_side(b[1], domain[2])) ||
         (trim_face_seam_only_near_side(a[1], domain[3]) &&
          trim_face_seam_only_near_side(b[1], domain[3]));
}

// True when the two endpoints are the two parameter copies of one physical
// seam point: one near each u side with the same v (a periodic u seam), or one
// near each v side with the same u (a periodic v seam). The other coordinate
// and the side distance use the same relative 1e-12 bound the projection uses
// for parameters; comparing them exactly would reject a chain whose seam copies
// drifted by rounding as open. Two different poles never form a pair.
inline bool trim_face_seam_only_seam_pair(
    const Point2d& a, const Point2d& b, const std::array<double, 4>& domain,
    bool periodic_u, bool periodic_v) noexcept {
  const auto same = [](double lhs, double rhs) noexcept {
    return std::fabs(lhs - rhs) <=
           1e-12 * std::max({1.0, std::fabs(lhs), std::fabs(rhs)});
  };
  if (periodic_u && same(a[1], b[1]) &&
      ((trim_face_seam_only_near_side(a[0], domain[0]) &&
        trim_face_seam_only_near_side(b[0], domain[1])) ||
       (trim_face_seam_only_near_side(a[0], domain[1]) &&
        trim_face_seam_only_near_side(b[0], domain[0])))) {
    return true;
  }
  if (periodic_v && same(a[0], b[0]) &&
      ((trim_face_seam_only_near_side(a[1], domain[2]) &&
        trim_face_seam_only_near_side(b[1], domain[3])) ||
       (trim_face_seam_only_near_side(a[1], domain[3]) &&
        trim_face_seam_only_near_side(b[1], domain[2])))) {
    return true;
  }
  return false;
}

// The ordering key of a ring: its lexicographic minimum. Unlike the ring front
// it does not depend on where the arrangement's walk started or which way the
// ring runs, so ring and fragment order are stable under cycle rotation.
inline Point2d trim_face_seam_only_ring_key(
    const std::vector<Point2d>& ring) noexcept {
  Point2d key = ring.front();
  for (const Point2d& point : ring) {
    if (trim_face_lex_less(point, key)) key = point;
  }
  return key;
}

// Splits the synthetic flags of one arrangement cycle into the maximal runs of
// cut edges, as indices local to the cycle. The walk starts at a run head so a
// run that straddles the end of the cycle array is emitted once and whole
// (`{6, 7, 0}` for the pattern `{0, 1, 1, 0, 0, 1, 0, 0}`); a no-wrapping walk
// would lose index 0 there. An all-synthetic cycle has no runs; an all-cut cycle
// is one closed run of its own.
inline std::vector<std::vector<std::size_t>> trim_face_seam_runs(
    std::span<const char> synthetic) {
  std::vector<std::vector<std::size_t>> runs;
  const std::size_t n = synthetic.size();
  if (n == 0) return runs;
  std::size_t cut_count = 0;
  std::size_t head = n;
  for (std::size_t k = 0; k < n; ++k) {
    if (synthetic[k] != 0) continue;
    ++cut_count;
    if (head == n && synthetic[(k + n - 1) % n] != 0) head = k;
  }
  if (cut_count == 0) return runs;
  if (head == n) {
    std::vector<std::size_t> run(n);
    for (std::size_t k = 0; k < n; ++k) run[k] = k;
    runs.push_back(std::move(run));
    return runs;
  }
  std::vector<std::size_t> current;
  for (std::size_t step = 0; step < n; ++step) {
    const std::size_t k = (head + step) % n;
    if (synthetic[k] != 0) {
      if (!current.empty()) {
        runs.push_back(std::move(current));
        current.clear();
      }
    } else {
      current.push_back(k);
    }
  }
  if (!current.empty()) runs.push_back(std::move(current));
  return runs;
}

// True when the edge lies along one domain side. The arrangement has no source
// tag on its half-edges, so this geometric test is how the synthetic chains are
// recognized after the build. Requiring the same side (not merely two endpoints
// that each touch some side) keeps a cut that joins two different sides from
// being misclassified as synthetic.
inline bool trim_face_seam_only_synthetic_edge(
    const Point2d& a, const Point2d& b,
    const std::array<double, 4>& domain) noexcept {
  return (a[0] == domain[0] && b[0] == domain[0]) ||
         (a[0] == domain[1] && b[0] == domain[1]) ||
         (a[1] == domain[2] && b[1] == domain[2]) ||
         (a[1] == domain[3] && b[1] == domain[3]);
}

// The sphere chart is the only seam-only chart whose v range is exactly
// [-pi/2, pi/2]; the cone's is [0, r/tan(theta)] and the torus's [0, 2pi]. On
// that chart a synthetic edge along the top or bottom side is a pole row, which
// collapses to one 3D point even though it spans the whole u range.
inline bool trim_face_seam_only_pole_v(const std::array<double, 4>& domain,
                                       double& pole_v) noexcept {
  constexpr double kHalfPi = 1.57079632679489661923;
  if (domain[2] == -kHalfPi) {
    pole_v = domain[2];
    return true;
  }
  if (domain[3] == kHalfPi) {
    pole_v = domain[3];
    return true;
  }
  return false;
}

// Appends one side of the synthetic rectangle. Interior samples interpolate the
// two corners, so every sample of a horizontal side keeps v exactly and every
// sample of a vertical side keeps u exactly; the side therefore stays on its
// domain line, which is what the synthetic-edge test relies on.
inline void trim_face_seam_only_append_side(
    std::vector<ArrangementSegment>& segments, const Point2d& from,
    const Point2d& to, double step) {
  const double du = to[0] - from[0];
  const double dv = to[1] - from[1];
  const double pieces_d = std::ceil(std::hypot(du, dv) / step);
  const std::size_t pieces =
      static_cast<std::size_t>(pieces_d < 1.0 ? 1.0 : pieces_d);
  Point2d previous = from;
  for (std::size_t i = 1; i <= pieces; ++i) {
    const double t = static_cast<double>(i) / static_cast<double>(pieces);
    const Point2d next =
        i == pieces ? to : Point2d{from[0] + t * du, from[1] + t * dv};
    segments.push_back(ArrangementSegment{previous, next, 0u});
    previous = next;
  }
}

// The public seam-only splitter still routes sphere faces here, so it keeps the
// historical input screening: an interior island, an open chain, a cut that
// leaves the rectangle or one that crosses the seam twice used to be refused
// before any arrangement work. Task 3 moves the dispatch, and only then does the
// public path accept those inputs; this screen keeps the old errors stable until
// it does. The geometric checks do not modify the cuts: the splitter below
// windows them again, on the same descriptors the caller passes.
inline Result<void> trim_face_seam_only_screen(
    FaceId face, const FaceUvProjection& projection,
    std::span<const FaceCut> cuts, const std::array<double, 4>& domain) {
  const double u_lo = domain[0];
  const double u_hi = domain[1];
  const double v_lo = domain[2];
  const double v_hi = domain[3];
  const double turn = 2.0 * kFaceUvPi;
  const bool periodic_v = projection.periodic_v();

  // Normalize every cut into the domain. The u unwrap is the same whole-turn
  // slide the projection uses; the chain is then moved as a whole into
  // [u_lo, u_hi], so a chain that crosses the seam is rejected instead of
  // silently wrapping.
  std::vector<std::vector<Point2d>> cut_points;
  cut_points.reserve(cuts.size());
  for (std::size_t c = 0; c < cuts.size(); ++c) {
    const std::vector<Point2d>& uv = cuts[c].uv;
    const std::uint32_t cut_index = static_cast<std::uint32_t>(c);
    if (uv.size() < 2) {
      return Error{Errc::invalid_topology,
                   "split_face_curved: a cut has fewer than two points",
                   face.index(), cut_index};
    }
    for (const Point2d& point : uv) {
      if (!trim_face_finite2(point)) {
        return Error{Errc::non_finite,
                     "split_face_curved: cut point is not finite", face.index(),
                     cut_index};
      }
    }
    std::vector<Point2d> points = uv;
    for (std::size_t i = 1; i < points.size(); ++i) {
      points[i][0] = face_uv_unwrap(points[i][0], points[i - 1][0]);
      if (periodic_v) {
        points[i][1] = face_uv_unwrap(points[i][1], points[i - 1][1]);
      }
    }
    double min_u = points.front()[0];
    double max_u = min_u;
    for (const Point2d& point : points) {
      min_u = std::min(min_u, point[0]);
      max_u = std::max(max_u, point[0]);
    }
    const double shift = -std::floor((min_u - u_lo) / turn) * turn;
    for (Point2d& point : points) point[0] += shift;
    min_u += shift;
    max_u += shift;
    // An SSI cut endpoint can project a few ulp inside the domain. The seam-pair
    // closure tolerates that drift, but the arrangement does not: the cut would
    // leave a dangling chain and build a zero-width face. Snap a coordinate that
    // close to a side onto that side before the arrangement.
    const double seam_eps =
        1e-12 * std::max({1.0, std::fabs(u_lo), std::fabs(u_hi)});
    for (Point2d& point : points) {
      if (std::fabs(point[0] - u_lo) <= seam_eps) {
        point[0] = u_lo;
      } else if (std::fabs(point[0] - u_hi) <= seam_eps) {
        point[0] = u_hi;
      }
      if (periodic_v) {
        if (std::fabs(point[1] - v_lo) <= seam_eps) {
          point[1] = v_lo;
        } else if (std::fabs(point[1] - v_hi) <= seam_eps) {
          point[1] = v_hi;
        }
      }
    }
    // The bounds check below uses these; the values above predate the snap.
    min_u = points.front()[0];
    max_u = min_u;
    for (const Point2d& point : points) {
      min_u = std::min(min_u, point[0]);
      max_u = std::max(max_u, point[0]);
    }
    if (!(min_u >= u_lo) || !(max_u <= u_hi)) {
      return Error{Errc::unsupported,
                   "split_face_curved: cut leaves the parameter domain",
                   face.index(), cut_index};
    }
    if (periodic_v) {
      double min_v = points.front()[1];
      for (const Point2d& point : points) min_v = std::min(min_v, point[1]);
      const double v_shift = -std::floor((min_v - v_lo) / turn) * turn;
      for (Point2d& point : points) point[1] += v_shift;
      for (const Point2d& point : points) {
        if (!(point[1] >= v_lo) || !(point[1] <= v_hi)) {
          return Error{Errc::unsupported,
                       "split_face_curved: cut leaves the parameter domain",
                       face.index(), cut_index};
        }
      }
    } else {
      for (Point2d& point : points) {
        const double slack = 1e-9 * std::max(1.0, std::fabs(point[1]));
        if (point[1] < v_lo - slack || point[1] > v_hi + slack) {
          return Error{Errc::unsupported,
                       "split_face_curved: cut leaves the parameter domain",
                       face.index(), cut_index};
        }
        // A pole touched within the allowed slack is the pole exactly.
        if (point[1] < v_lo) {
          point[1] = v_lo;
        } else if (point[1] > v_hi) {
          point[1] = v_hi;
        }
      }
    }
    // A cut that runs along a domain side would merge with the synthetic chain
    // on that side; after the strip step its geometry would look synthetic and
    // disappear. Reject it before the arrangement so no partial work is done.
    for (std::size_t i = 1; i < points.size(); ++i) {
      if (trim_face_seam_only_on_same_side(points[i - 1], points[i], domain)) {
        return Error{Errc::unsupported,
                     "split_face_curved: cut lies on the parameter domain "
                     "boundary",
                     face.index(), cut_index};
      }
    }
    // A cut chain is closed when it ends where it began, or when its two ends
    // are the two copies of one seam point (both tests use the same 1e-12
    // relative bound; ring closure below uses it too). Anything else cannot
    // form a loop of cut geometry on this face.
    const bool closed =
        trim_face_same_xy(points.front(), points.back()) ||
        trim_face_seam_only_seam_pair(points.front(), points.back(), domain,
                                      projection.periodic_u(), periodic_v);
    if (!closed) {
      return Error{Errc::unsupported,
                   "split_face_curved: open cuts on a seam-only face",
                   face.index(), cut_index};
    }
    // This splitter only joins cut runs that end on the seam. A closed cut
    // wholly inside the domain would need its two arcs joined into one ring
    // (one side would otherwise be lost), and a closed cut that crosses the
    // seam at two places would need several runs joined; both are rejected
    // here, before any arrangement work, so the caller never gets a partial
    // split. A latitude-like cut has exactly one seam point, realized by its
    // two ends.
    std::vector<double> seam_contacts;
    for (const Point2d& point : points) {
      if (!trim_face_seam_only_near_side(point[0], u_lo) &&
          !trim_face_seam_only_near_side(point[0], u_hi)) {
        continue;
      }
      bool seen = false;
      for (const double v : seam_contacts) {
        if (std::fabs(v - point[1]) <=
            1e-12 * std::max({1.0, std::fabs(v), std::fabs(point[1])})) {
          seen = true;
          break;
        }
      }
      if (!seen) seam_contacts.push_back(point[1]);
    }
    if (seam_contacts.empty()) {
      return Error{Errc::unsupported,
                   "split_face_curved: closed cut does not reach the seam",
                   face.index(), cut_index};
    }
    if (seam_contacts.size() > 1 ||
        trim_face_same_xy(points.front(), points.back()) ||
        !trim_face_seam_only_seam_pair(points.front(), points.back(), domain,
                                       projection.periodic_u(), periodic_v)) {
      return Error{Errc::unsupported,
                   "split_face_curved: closed cut crosses the seam twice",
                   face.index(), cut_index};
    }
    cut_points.push_back(std::move(points));
  }
  for (std::size_t c = 0; c < cut_points.size(); ++c) {
    bool degenerate = true;
    for (const Point2d& point : cut_points[c]) {
      if (!trim_face_same_xy(point, cut_points[c].front())) {
        degenerate = false;
        break;
      }
    }
    if (degenerate) {
      return Error{Errc::degenerate,
                   "split_face_curved: cut collapses to a point", face.index(),
                   static_cast<std::uint32_t>(c)};
    }
  }
  return Result<void>{};
}

// The side of the domain an arrangement edge lies on, in the TrimFaceDomain
// order (0 = u_lo, 1 = u_hi, 2 = v_lo, 3 = v_hi), or -1 for a real cut edge. The
// arrangement carries no source tag on its half-edges, so this geometric test
// is how the synthetic chains are recognized after the build. Requiring the
// same side (not merely two endpoints that each touch some side) keeps a cut
// that joins two different sides from being misclassified as synthetic.
inline int trim_face_domain_edge_side(const Point2d& a, const Point2d& b,
                                      const TrimFaceDomain& domain) noexcept {
  if (a[0] == domain.bounds[0] && b[0] == domain.bounds[0]) return 0;
  if (a[0] == domain.bounds[1] && b[0] == domain.bounds[1]) return 1;
  if (a[1] == domain.bounds[2] && b[1] == domain.bounds[2]) return 2;
  if (a[1] == domain.bounds[3] && b[1] == domain.bounds[3]) return 3;
  return -1;
}

// True when the point lies strictly inside the closed planar ring, using the
// even-odd rule on the exact orientation predicate. Boundary points are not
// inside: a ring does not contain the vertices it shares with its complement
// fragment, and the island keyhole reuses the island vertices exactly.
inline bool trim_face_point_strictly_inside_ring(
    const Point2d& point, const std::vector<Point2d>& ring) noexcept {
  const std::size_t n = ring.size();
  if (n < 3) return false;
  bool inside = false;
  for (std::size_t i = 0; i < n; ++i) {
    const Point2d& a = ring[(i + n - 1) % n];
    const Point2d& b = ring[i];
    if (a[0] == b[0] && a[1] == b[1]) continue;
    if ((a[1] > point[1]) == (b[1] > point[1])) continue;
    const int side = sign_orient2d(a[0], a[1], b[0], b[1], point[0], point[1]);
    if (side == 0) return false;
    if ((side > 0) == (b[1] > a[1])) inside = !inside;
  }
  return inside;
}

// True when the two closed rings cross each other properly. Two rings that do
// not cross are either disjoint or nested; crossing rings are the overlapping
// islands case, which the area partition guard owns.
inline bool trim_face_rings_cross(const std::vector<Point2d>& first,
                                  const std::vector<Point2d>& second) noexcept {
  const std::size_t n = first.size();
  const std::size_t m = second.size();
  for (std::size_t i = 1; i < n; ++i) {
    for (std::size_t j = 1; j < m; ++j) {
      const Point2d& a = first[i - 1];
      const Point2d& b = first[i];
      const Point2d& c = second[j - 1];
      const Point2d& d = second[j];
      if ((a[0] == b[0] && a[1] == b[1]) ||
          (c[0] == d[0] && c[1] == d[1])) {
        continue;
      }
      const int o1 = sign_orient2d(a[0], a[1], b[0], b[1], c[0], c[1]);
      const int o2 = sign_orient2d(a[0], a[1], b[0], b[1], d[0], d[1]);
      const int o3 = sign_orient2d(c[0], c[1], d[0], d[1], a[0], a[1]);
      const int o4 = sign_orient2d(c[0], c[1], d[0], d[1], b[0], b[1]);
      if (o1 * o2 < 0 && o3 * o4 < 0) return true;
    }
  }
  return false;
}

// The transverse coordinates where a closed cut passes through one periodic
// seam: the interpolated value at each u/v wrap, plus the value at any interior
// vertex that lies exactly on the seam. A chain whose two ends are the seam
// copies is closed by that seam itself, so its end vertices are not crossings.
inline std::vector<double> trim_face_seam_crossings(
    const std::vector<Point2d>& points, const TrimFaceDomain& domain,
    int axis) {
  std::vector<double> crossings;
  if (points.size() < 3) return crossings;
  const int transverse = axis == 0 ? 1 : 0;
  const double lo = domain.bounds[2 * axis];
  const double hi = domain.bounds[2 * axis + 1];
  const double period = hi - lo;
  if (!(period > 0.0)) return crossings;
  std::vector<double> unwrapped(points.size(), points.front()[axis]);
  for (std::size_t i = 1; i < points.size(); ++i) {
    double value = points[i][axis];
    while (value - unwrapped[i - 1] > 0.5 * period) value -= period;
    while (unwrapped[i - 1] - value > 0.5 * period) value += period;
    unwrapped[i] = value;
  }
  for (std::size_t i = 0; i + 1 < points.size(); ++i) {
    const double from = unwrapped[i];
    const double to = unwrapped[i + 1];
    if (from == to) continue;
    const double band_from = std::floor((from - lo) / period);
    const double band_to = std::floor((to - lo) / period);
    if (band_from == band_to) continue;
    const double boundary =
        band_to > band_from ? lo + (band_from + 1.0) * period
                            : lo + band_from * period;
    const double t = (boundary - from) / (to - from);
    // A crossing at a chain vertex is that vertex (handled below), not an
    // additional pass through the seam.
    if (!(t > 1e-12) || !(t < 1.0 - 1e-12)) continue;
    crossings.push_back(points[i][transverse] +
                        t * (points[i + 1][transverse] -
                             points[i][transverse]));
  }
  for (std::size_t i = 1; i + 1 < points.size(); ++i) {
    if (trim_face_seam_only_near_side(points[i][axis], lo) ||
        trim_face_seam_only_near_side(points[i][axis], hi)) {
      crossings.push_back(points[i][transverse]);
    }
  }
  return crossings;
}

// True when a closed chain crosses a periodic seam twice with different
// transverse coordinates. Equal coordinates pair up (the tail can be translated
// by one period without changing the shoelace); a twisted pair cannot, because
// the closing chord would cut across the chart instead of following the seam.
inline bool trim_face_closed_cut_twisted(const std::vector<Point2d>& points,
                                         const TrimFaceDomain& domain,
                                         double tolerance) {
  for (int axis = 0; axis < 2; ++axis) {
    const bool periodic = axis == 0 ? domain.periodic_u : domain.periodic_v;
    if (!periodic) continue;
    const int transverse = axis == 0 ? 1 : 0;
    const double lo = domain.bounds[2 * axis];
    const double hi = domain.bounds[2 * axis + 1];
    std::vector<double> crossings =
        trim_face_seam_crossings(points, domain, axis);
    // A chain whose two ends are the two copies of this seam is closed by the
    // seam itself: the closing segment runs along the seam, so it crosses at
    // the tail's transverse coordinate even though no edge does. Without this
    // crossing a ring that enters and leaves through the seam at different
    // transverse coordinates would see fewer than two crossings and slip
    // through to the arrangement, which then fails with an open-chain error.
    if (trim_face_seam_only_seam_pair(points.front(), points.back(),
                                      domain.bounds, domain.periodic_u,
                                      domain.periodic_v) &&
        (trim_face_seam_only_near_side(points.back()[axis], lo) ||
         trim_face_seam_only_near_side(points.back()[axis], hi))) {
      crossings.push_back(points.back()[transverse]);
    }
    // A ring closed by repeating its first point (rather than by the two copies
    // of a seam point) that starts exactly on this seam has its closing pass
    // exactly at that vertex: the edge scan above skips t == 1, and the interior
    // vertex scan skips the duplicated tail, so the pass is counted here. Without
    // it a ring that straddles the seam through its own start point sees a single
    // crossing and is misread as an open chain instead of a twisted one.
    if (trim_face_same_xy_rel(points.front(), points.back(), tolerance) &&
        (trim_face_seam_only_near_side(points.front()[axis], lo) ||
         trim_face_seam_only_near_side(points.front()[axis], hi))) {
      crossings.push_back(points.front()[transverse]);
    }
    if (crossings.size() < 2) continue;
    if (crossings.size() % 2 != 0) return true;
    std::sort(crossings.begin(), crossings.end());
    for (std::size_t i = 0; i < crossings.size(); i += 2) {
      if (std::fabs(crossings[i + 1] - crossings[i]) > tolerance) return true;
    }
  }
  return false;
}

// Chooses the midpoint of the widest free band on a periodic axis. `covered`
// lists the intervals already occupied by cuts as {lo, hi}; an interval with
// hi < lo wraps through the period. Overlaps are merged, the midpoint of the
// widest remaining band is normalized into [0, period), and an exact tie goes to
// the lower midpoint. std::nullopt means the axis is fully covered, the widest
// free band is narrower than 2*`margin` (so no point keeps the requested
// clearance from the cuts), the arguments are not finite, or more than
// `max_intervals` intervals were supplied.
inline std::optional<double> trim_face_periodic_free_midpoint(
    std::span<const std::array<double, 2>> covered, double period,
    double margin, std::size_t max_intervals = 4) {
  if (!std::isfinite(period) || !(period > 0.0) || !std::isfinite(margin) ||
      !(margin >= 0.0)) {
    return std::nullopt;
  }
  if (covered.size() > max_intervals) return std::nullopt;
  std::vector<std::array<double, 2>> bands;
  bands.reserve(2 * covered.size());
  for (const std::array<double, 2>& interval : covered) {
    if (!std::isfinite(interval[0]) || !std::isfinite(interval[1])) {
      return std::nullopt;
    }
    double width = interval[1] - interval[0];
    if (width >= period) return std::nullopt;  // nothing is left free
    // A wrapped interval of a full turn or more covers the circle as well.
    // Without this check width += period turns exactly -period into a
    // zero-width point, which would incorrectly leave the whole circle free.
    if (width <= -period) return std::nullopt;
    if (width == 0.0) continue;               // a point blocks no band
    if (width < 0.0) width += period;         // wrapped interval
    double lo = interval[0] - period * std::floor(interval[0] / period);
    if (lo < 0.0 || lo >= period) lo = 0.0;  // guard against rounding
    const double hi = lo + width;
    if (hi <= period) {
      bands.push_back(std::array<double, 2>{lo, hi});
    } else {
      bands.push_back(std::array<double, 2>{lo, period});
      bands.push_back(std::array<double, 2>{0.0, hi - period});
    }
  }
  std::sort(bands.begin(), bands.end());
  std::vector<std::array<double, 2>> merged;
  merged.reserve(bands.size());
  for (const std::array<double, 2>& band : bands) {
    if (!merged.empty() && band[0] <= merged.back()[1]) {
      merged.back()[1] = std::max(merged.back()[1], band[1]);
    } else {
      merged.push_back(band);
    }
  }
  double best_length = -1.0;
  double best_midpoint = 0.0;
  const auto consider = [&](double gap_start, double gap_end) {
    const double length = gap_end - gap_start;
    if (!(length > 0.0)) return;
    if (length > best_length) {
      best_length = length;
      best_midpoint = 0.5 * (gap_start + gap_end);
    }
  };
  if (merged.empty()) {
    consider(0.0, period);  // the whole circle is free
  } else {
    for (std::size_t i = 0; i < merged.size(); ++i) {
      const double gap_start = merged[i][1];
      const double gap_end =
          i + 1 < merged.size() ? merged[i + 1][0] : merged[0][0] + period;
      consider(gap_start, gap_end);
    }
  }
  if (best_length < 2.0 * margin) return std::nullopt;
  double midpoint = best_midpoint - period * std::floor(best_midpoint / period);
  if (midpoint < 0.0 || midpoint >= period) midpoint = 0.0;
  return midpoint;
}

// Splits a synthetic parameter rectangle along cut chains. The rectangle sides
// are synthetic: a seam copy is dropped (the ring closes across the seam), a
// collapsed side becomes the single pole point, and a real boundary side stays
// in the ring as if it were a cut. Closed chains that touch no side bound
// material on both sides, so they are emitted as an island plus the complement
// of the rectangle around them.
//
// `holes` are the stored inner uv rings of the face (material outside them).
// They are handled together with the other cut geometry in Task 3; this entry
// point only splits on `cuts`, so an empty span is the expected input here.
inline Result<std::vector<CurvedFragment>> trim_face_split_domain(
    FaceId face, const FaceUvProjection& projection,
    std::span<const FaceCut> cuts, const TrimFaceDomain& domain,
    std::span<const std::vector<Point2d>> holes, double tolerance) {
  if (!holes.empty()) {
    return Error{Errc::unsupported,
                 "split_face_curved: stored holes are not supported yet",
                 face.index(), 0u};
  }
  const double u_lo = domain.bounds[0];
  const double u_hi = domain.bounds[1];
  const double v_lo = domain.bounds[2];
  const double v_hi = domain.bounds[3];

  // A closed chain that crosses a periodic seam is windowed into open pieces
  // and rejoined through the seam. The join is only area-preserving when the two
  // crossings keep the same transverse coordinate: translating the tail by one
  // period must not move it. A twisted crossing (different transverse
  // coordinates) would need the chart re-origined, so refuse it before any
  // arrangement work instead of building a chord across the chart.
  for (std::size_t c = 0; c < cuts.size(); ++c) {
    const std::vector<Point2d>& uv = cuts[c].uv;
    if (uv.size() < 2) continue;
    const bool closed =
        trim_face_same_xy_rel(uv.front(), uv.back(), tolerance) ||
        trim_face_seam_only_seam_pair(uv.front(), uv.back(), domain.bounds,
                                      domain.periodic_u, domain.periodic_v);
    if (!closed) continue;
    if (trim_face_closed_cut_twisted(uv, domain, tolerance)) {
      return Error{Errc::unsupported,
                   "split_face_curved: closed cut crosses the seam at "
                   "different transverse coordinates; twisted seam crossings "
                   "are deferred to M3c-3c-4b3c",
                   face.index(), static_cast<std::uint32_t>(c)};
    }
  }

  // The cut polylines are windowed first, so every piece lies inside the
  // rectangle and its open ends sit exactly on a side. A closed island comes
  // back as one closed piece; a latitude-like chain is one piece whose two ends
  // are the two copies of one seam point.
  Result<std::vector<std::vector<Point2d>>> windowed =
      trim_face_window_chains(cuts, domain, tolerance);
  if (!windowed) return windowed.error();
  std::vector<std::vector<Point2d>> pieces = std::move(windowed.value());

  const double slack_u =
      1e-9 * std::max({1.0, std::fabs(u_lo), std::fabs(u_hi)});
  const double slack_v =
      1e-9 * std::max({1.0, std::fabs(v_lo), std::fabs(v_hi)});
  for (std::size_t i = 0; i < pieces.size(); ++i) {
    std::vector<Point2d>& piece = pieces[i];
    for (const Point2d& point : piece) {
      if (!trim_face_finite2(point)) {
        return Error{Errc::non_finite,
                     "split_face_curved: cut point is not finite", face.index(),
                     static_cast<std::uint32_t>(i)};
      }
    }
    // SSI seam endpoints can land a few ulp inside the chart. The seam pairing
    // and the synthetic-side recognition are exact, so snap a coordinate that
    // close to a side onto that side; otherwise the cut leaves a
    // near-degenerate sliver along the seam and its ring no longer closes.
    const double seam_eps_u =
        1e-12 * std::max({1.0, std::fabs(u_lo), std::fabs(u_hi)});
    const double seam_eps_v =
        1e-12 * std::max({1.0, std::fabs(v_lo), std::fabs(v_hi)});
    for (Point2d& point : piece) {
      if (std::fabs(point[0] - u_lo) <= seam_eps_u) {
        point[0] = u_lo;
      } else if (std::fabs(point[0] - u_hi) <= seam_eps_u) {
        point[0] = u_hi;
      }
      if (std::fabs(point[1] - v_lo) <= seam_eps_v) {
        point[1] = v_lo;
      } else if (std::fabs(point[1] - v_hi) <= seam_eps_v) {
        point[1] = v_hi;
      }
    }
    for (Point2d& point : piece) {
      if (!domain.periodic_u) {
        if (point[0] < u_lo - slack_u || point[0] > u_hi + slack_u) {
          return Error{Errc::unsupported,
                       "split_face_curved: cut leaves the parameter domain",
                       face.index(), static_cast<std::uint32_t>(i)};
        }
        if (point[0] < u_lo) {
          point[0] = u_lo;
        } else if (point[0] > u_hi) {
          point[0] = u_hi;
        }
      }
      if (!domain.periodic_v) {
        if (point[1] < v_lo - slack_v || point[1] > v_hi + slack_v) {
          return Error{Errc::unsupported,
                       "split_face_curved: cut leaves the parameter domain",
                       face.index(), static_cast<std::uint32_t>(i)};
        }
        // A pole touched within the allowed slack is the pole exactly.
        if (point[1] < v_lo) {
          point[1] = v_lo;
        } else if (point[1] > v_hi) {
          point[1] = v_hi;
        }
      }
    }
    for (std::size_t k = 1; k < piece.size(); ++k) {
      if (trim_face_seam_only_on_same_side(piece[k - 1], piece[k],
                                           domain.bounds)) {
        return Error{Errc::unsupported,
                     "split_face_curved: cut lies on the parameter domain "
                     "boundary",
                     face.index(), static_cast<std::uint32_t>(i)};
      }
    }
  }

  // One entry per glued chain. A closed chain that touches no side is an
  // interior island: it bounds material on both sides inside the rectangle.
  struct DomainChain {
    std::vector<Point2d> points;
    bool island = false;
  };
  std::vector<DomainChain> chains;
  std::vector<char> paired(pieces.size(), 0);
  for (std::size_t i = 0; i < pieces.size(); ++i) {
    if (paired[i]) continue;
    std::vector<Point2d>& piece = pieces[i];
    if (piece.size() < 2) continue;
    paired[i] = 1;
    const bool closed =
        trim_face_same_xy_rel(piece.front(), piece.back(), tolerance);
    const bool seam_closed =
        trim_face_seam_only_seam_pair(piece.front(), piece.back(),
                                      domain.bounds, domain.periodic_u,
                                      domain.periodic_v);
    if (closed || seam_closed) {
      DomainChain chain;
      chain.points = std::move(piece);
      if (closed) {
        chain.island = true;
        for (const Point2d& point : chain.points) {
          if (trim_face_point_on_domain_side(point, domain, tolerance)) {
            chain.island = false;
            break;
          }
        }
      }
      chains.push_back(std::move(chain));
      continue;
    }
    // Two open pieces that meet at the seam are the two halves of one chain
    // that crosses the seam: join them in order and let the arrangement close
    // the physical loop across the seam.
    std::size_t match = pieces.size();
    for (std::size_t j = i + 1; j < pieces.size(); ++j) {
      if (paired[j] || pieces[j].size() < 2) continue;
      if (trim_face_same_xy_rel(pieces[j].front(), pieces[j].back(),
                                tolerance)) {
        continue;
      }
      if (trim_face_seam_only_seam_pair(piece.back(), pieces[j].front(),
                                        domain.bounds, domain.periodic_u,
                                        domain.periodic_v) &&
          trim_face_seam_only_seam_pair(pieces[j].back(), piece.front(),
                                        domain.bounds, domain.periodic_u,
                                        domain.periodic_v)) {
        match = j;
        break;
      }
    }
    if (match == pieces.size()) {
      // An open piece whose ends lie on two non-seam domain sides cuts the face
      // from boundary to boundary: the rest of the fragment's boundary is the
      // synthetic sides themselves, so the piece needs no partner across the
      // seam. A piece that leaves through a periodic seam and returns elsewhere
      // would need the chart re-origined at the cut, so refuse the mix instead
      // of closing it with a chord across the chart.
      const bool front_on_side =
          trim_face_point_on_domain_side(piece.front(), domain, tolerance);
      const bool back_on_side =
          trim_face_point_on_domain_side(piece.back(), domain, tolerance);
      const auto on_seam = [&domain, tolerance](const Point2d& point) noexcept {
        if (domain.periodic_u &&
            (trim_face_near_side(point[0], domain.bounds[0], tolerance) ||
             trim_face_near_side(point[0], domain.bounds[1], tolerance))) {
          return true;
        }
        return domain.periodic_v &&
               (trim_face_near_side(point[1], domain.bounds[2], tolerance) ||
                trim_face_near_side(point[1], domain.bounds[3], tolerance));
      };
      const bool front_on_seam = on_seam(piece.front());
      const bool back_on_seam = on_seam(piece.back());
      if (front_on_seam != back_on_seam) {
        return Error{Errc::unsupported,
                     "split_face_curved: open cut crosses the seam",
                     face.index(), static_cast<std::uint32_t>(i)};
      }
      if (!front_on_side || !back_on_side || front_on_seam || back_on_seam) {
        return Error{Errc::unsupported,
                     "split_face_curved: open cuts on a seam-only face",
                     face.index(), static_cast<std::uint32_t>(i)};
      }
      DomainChain chain;
      chain.points = std::move(piece);
      // The synthetic sides are sampled at fixed steps, so an open end that
      // stops between two samples must pin the exact side value here; otherwise
      // the side chord misses the cut endpoint and the fragment does not close.
      const auto snap = [u_lo, u_hi, v_lo, v_hi, tolerance](Point2d& point) {
        if (trim_face_near_side(point[0], u_lo, tolerance)) {
          point[0] = u_lo;
        } else if (trim_face_near_side(point[0], u_hi, tolerance)) {
          point[0] = u_hi;
        }
        if (trim_face_near_side(point[1], v_lo, tolerance)) {
          point[1] = v_lo;
        } else if (trim_face_near_side(point[1], v_hi, tolerance)) {
          point[1] = v_hi;
        }
      };
      snap(chain.points.front());
      snap(chain.points.back());
      chains.push_back(std::move(chain));
      continue;
    }
    paired[match] = 1;
    DomainChain chain;
    chain.points = std::move(piece);
    chain.points.insert(chain.points.end(), pieces[match].begin(),
                        pieces[match].end());
    chains.push_back(std::move(chain));
  }
  for (std::size_t c = 0; c < chains.size(); ++c) {
    bool degenerate = true;
    for (const Point2d& point : chains[c].points) {
      if (!trim_face_same_xy(point, chains[c].points.front())) {
        degenerate = false;
        break;
      }
    }
    if (degenerate) {
      return Error{Errc::degenerate,
                   "split_face_curved: cut collapses to a point", face.index(),
                   static_cast<std::uint32_t>(c)};
    }
  }

  // Nested closed islands are deferred: the arrangement keeps one cycle per
  // face, so a fragment cannot carry the inner island as a hole and would
  // instead claim its area. Nested rings do not cross, while overlapping
  // islands (the case the partition guard below covers) do, so containment
  // without a crossing is unambiguous nesting.
  for (std::size_t a = 0; a < chains.size(); ++a) {
    if (!chains[a].island) continue;
    for (std::size_t b = 0; b < chains.size(); ++b) {
      if (a == b || !chains[b].island) continue;
      if (trim_face_rings_cross(chains[a].points, chains[b].points)) continue;
      if (trim_face_point_strictly_inside_ring(chains[a].points.front(),
                                               chains[b].points)) {
        return Error{Errc::numeric_failure,
                     "split_face_curved: nested closed cuts are not supported",
                     face.index(), 0u};
      }
    }
  }

  // Four synthetic sides, counter-clockwise: bottom, right, top, left. The
  // classification below recovers them by their geometry, not by this order.
  const double step = kFaceUvPi / 32.0;
  std::vector<ArrangementSegment> segments;
  segments.reserve(4u * 65u + chains.size() * 64u);
  const Point2d bottom_left{u_lo, v_lo};
  const Point2d bottom_right{u_hi, v_lo};
  const Point2d top_right{u_hi, v_hi};
  const Point2d top_left{u_lo, v_hi};
  trim_face_seam_only_append_side(segments, bottom_left, bottom_right, step);
  trim_face_seam_only_append_side(segments, bottom_right, top_right, step);
  trim_face_seam_only_append_side(segments, top_right, top_left, step);
  trim_face_seam_only_append_side(segments, top_left, bottom_left, step);
  for (std::size_t c = 0; c < chains.size(); ++c) {
    const std::uint32_t source = static_cast<std::uint32_t>(1u + c);
    for (std::size_t i = 0; i + 1 < chains[c].points.size(); ++i) {
      if (trim_face_same_xy(chains[c].points[i], chains[c].points[i + 1])) {
        continue;
      }
      segments.push_back(ArrangementSegment{chains[c].points[i],
                                            chains[c].points[i + 1], source});
    }
  }

  const SplitResult split = split_segments(segments);
  if (!std::isfinite(split.max_residual)) {
    return Error{Errc::numeric_failure,
                 "split_face_curved: segment intersection is not representable",
                 face.index(), 0u};
  }
  const Arrangement2 arrangement = Arrangement2::build(segments);

  const std::span<const Point2d> vertices = arrangement.vertices();
  const std::span<const Arrangement2::HalfEdge> half_edges =
      arrangement.half_edges();
  const std::span<const Arrangement2::Face> faces = arrangement.faces();

  std::vector<CurvedFragment> fragments;
  for (std::uint32_t f = 0; f < arrangement.face_count(); ++f) {
    if (!faces[f].bounded) continue;
    // The synthetic rectangle is source 0, so a nonzero winding of source 0 is
    // exactly the material side of the synthetic domain.
    if (arrangement.winding(f, 0u) == 0) continue;

    const std::uint32_t start = faces[f].edge;
    std::vector<std::uint32_t> cycle;
    std::uint32_t h = start;
    do {
      cycle.push_back(h);
      h = half_edges[h].next;
    } while (h != start && cycle.size() <= half_edges.size());
    if (h != start || cycle.size() < 3) continue;

    // A material cycle that owns part of a collapsed side reaches that pole:
    // the side is one 3D point, so the ring is the rim plus exactly one pole
    // point instead of a seam-pair closure.
    double pole_v = 0.0;
    bool touches_pole = false;
    std::vector<char> synthetic(cycle.size(), 0);
    std::size_t kept_edges = 0;
    for (std::size_t i = 0; i < cycle.size(); ++i) {
      const Arrangement2::HalfEdge& edge = half_edges[cycle[i]];
      const Point2d& from = vertices[edge.origin];
      const Point2d& to = vertices[half_edges[edge.twin].origin];
      const int side = trim_face_domain_edge_side(from, to, domain);
      if (side < 0) {
        ++kept_edges;
        continue;
      }
      const TrimSideKind kind = domain.sides[static_cast<std::size_t>(side)];
      if (kind == TrimSideKind::kRealBoundary) {
        // A real boundary is material: it stays in the ring like a cut edge.
        ++kept_edges;
        continue;
      }
      synthetic[i] = 1;
      if (kind == TrimSideKind::kCollapsed) {
        pole_v =
            domain.bounds[side == 2 ? 2u : 3u];
        touches_pole = true;
      }
    }
    if (kept_edges == 0) continue;

    struct DomainRing {
      std::vector<Point2d> points;
      bool seam = false;
    };
    std::vector<DomainRing> rings;

    // A domain that is periodic on both axes and has no pole is a torus chart.
    // A band across it is already one closed walk: the cycle carries both rims
    // and both seam columns, so walking it once closes the band. Emitting the
    // rims as separate seam rings would put two zero-area boundaries on one
    // region, so the whole cycle becomes a single ring. A chart that is
    // collapsed in v keeps the rim runs: its band is joined below instead.
    bool built_whole_cycle = false;
    if (domain.periodic_u && domain.periodic_v && !domain.collapsed_v()) {
      // Start at a real (cut) edge where one exists, never at the seam column.
      std::size_t first_kept = cycle.size();
      for (std::size_t i = 0; i < cycle.size(); ++i) {
        if (!synthetic[i]) {
          first_kept = i;
          break;
        }
      }
      if (first_kept == cycle.size()) first_kept = 0;
      std::vector<Point2d> whole;
      whole.reserve(cycle.size() + 1);
      bool seam = false;
      for (std::size_t k = 0; k < cycle.size(); ++k) {
        const std::size_t index = (first_kept + k) % cycle.size();
        const Arrangement2::HalfEdge& edge = half_edges[cycle[index]];
        whole.push_back(vertices[edge.origin]);
        seam = seam || (synthetic[index] != 0);
      }
      whole.push_back(vertices[half_edges[cycle[first_kept]].origin]);
      trim_face_normalize_ring(whole, true);
      if (whole.size() >= 3) {
        rings.push_back(DomainRing{std::move(whole), seam});
        built_whole_cycle = true;
      }
    }

    // Split the cycle into maximal runs of kept edges. trim_face_seam_runs
    // starts its walk at a run head, so a run that straddles the end of the
    // cycle array is emitted once and whole (a cycle with no kept edge was
    // skipped above). A whole-cycle ring above already covers the cycle.
    std::vector<std::vector<std::size_t>> runs;
    if (!built_whole_cycle) {
      runs = trim_face_seam_runs(std::span<const char>(synthetic));
    }
    for (const std::vector<std::size_t>& run : runs) {
      const Arrangement2::HalfEdge& first_edge =
          half_edges[cycle[run.front()]];
      std::vector<Point2d> points;
      points.reserve(run.size() + 1);
      points.push_back(vertices[first_edge.origin]);
      for (const std::size_t index : run) {
        const Arrangement2::HalfEdge& edge = half_edges[cycle[index]];
        points.push_back(vertices[half_edges[edge.twin].origin]);
      }
      if (touches_pole) {
        // The rim ends on a pole row, and that row is one 3D point: close the
        // ring through the exact pole between the last rim sample and the
        // repeated first point, i.e. along the meridian chords rim -> pole ->
        // rim. The pole comes from the domain, so v is exactly the panel edge.
        points.push_back(Point2d{points.back()[0], pole_v});
        if (domain.collapsed_row_is_linear_apex &&
            !trim_face_seam_only_near_side(points.front()[1], pole_v)) {
          // The run started on the opposite side row, so the ring must still
          // travel the whole collapsed row to the front's meridian: the row is
          // one 3D point, and its copy at the front's u carries that travel.
          // Without it the ring closes with a uv diagonal from the pole to the
          // front, which under-covers the chart by the front's u extent and
          // leaves the pcurve off the ring.
          points.push_back(Point2d{points.front()[0], pole_v});
        }
        points.push_back(points.front());
        if (points.size() < 3) continue;
        trim_face_normalize_ring(points, true);
        rings.push_back(DomainRing{std::move(points), false});
        continue;
      }
      const Point2d front = points.front();
      const Point2d back = points.back();
      const Point3d front3 = face_uv_eval(projection, front);
      const Point3d back3 = face_uv_eval(projection, back);
      const bool identical = trim_face_finite3(front3) &&
                             trim_face_finite3(back3) &&
                             front3[0] == back3[0] && front3[1] == back3[1] &&
                             front3[2] == back3[2];
      bool seam = false;
      if (!identical) {
        // The two u sides are two parameter copies of one physical seam, so a
        // run that leaves at (u_lo, v) and returns at (u_hi, v) is one
        // physical loop; the v sides of a torus are one seam in the same way.
        // Two different poles never close: a chain from one pole to the other
        // is a meridian, not a loop. Closure uses the same 1e-12 relative bound
        // as the input classification, so a seam copy that drifted by rounding
        // still closes; the ring is then closed by repeating its first point.
        seam = trim_face_seam_only_seam_pair(front, back, domain.bounds,
                                             domain.periodic_u,
                                             domain.periodic_v);
        if (!seam) {
          // A run can also enter and leave through the same parameter copy of a
          // seam, as a bound-from-boundary cut does on a chart whose seam is
          // not touched by the cut. The seam generatrix itself is then the
          // fragment's boundary edge, so the ring closes along that copy
          // instead of across the chart. The chord is a real edge, not a chart
          // artifact, so the ring stays a normal (non-seam) ring.
          const bool same_copy =
              (domain.periodic_u &&
               ((trim_face_seam_only_near_side(front[0], domain.bounds[0]) &&
                 trim_face_seam_only_near_side(back[0], domain.bounds[0])) ||
                (trim_face_seam_only_near_side(front[0], domain.bounds[1]) &&
                 trim_face_seam_only_near_side(back[0], domain.bounds[1])))) ||
              (domain.periodic_v &&
               ((trim_face_seam_only_near_side(front[1], domain.bounds[2]) &&
                 trim_face_seam_only_near_side(back[1], domain.bounds[2])) ||
                (trim_face_seam_only_near_side(front[1], domain.bounds[3]) &&
                 trim_face_seam_only_near_side(back[1], domain.bounds[3]))));
          if (!same_copy) {
            return Error{Errc::numeric_failure,
                         "split_face_curved: synthetic boundary chain is open",
                         face.index(), f};
          }
        }
        points.push_back(front);
      } else {
        points.push_back(front);
        if (points.size() < 3) continue;
        trim_face_normalize_ring(points, true);
      }
      if (points.size() < 3) continue;
      // A seam ring is explicitly closed by the repeated first point, whose
      // other parameter copy is the point before the last. It keeps the
      // arrangement's direction and is not rotated.
      rings.push_back(DomainRing{std::move(points), seam});
    }
    if (rings.empty()) continue;

    // A pole-free fragment of a chart with a pole, bounded by exactly two
    // rings, is the band between two latitude runs: the two rims, joined
    // across the seam column at u = 0. Joining turns them into one closed
    // ring, which is what a band is; leaving them as two rings would put two
    // zero-area boundaries on one region. A cap has one ring (or touches a
    // pole) and is never joined. Both runs are seam rings here: their ends
    // are the two parameter copies of the seam, one whole turn apart, so the
    // seam flag itself cannot decide.
    if (domain.collapsed_v() && !touches_pole && rings.size() == 2) {
      std::vector<Point2d> merged;
      merged.reserve(rings[0].points.size() + rings[1].points.size());
      for (const DomainRing& ring : rings) {
        // Each seam ring repeats its first point at the end; the join closes
        // the merged ring instead, so that repeat is dropped.
        for (std::size_t i = 0; i + 1 < ring.points.size(); ++i) {
          if (!merged.empty() && merged.back()[0] == ring.points[i][0] &&
              merged.back()[1] == ring.points[i][1]) {
            continue;
          }
          merged.push_back(ring.points[i]);
        }
      }
      merged.push_back(merged.front());
      if (merged.size() >= 4) {
        trim_face_normalize_ring(merged, true);
        rings.clear();
        rings.push_back(DomainRing{std::move(merged), false});
      }
    }

    // A seam ring closes across the periodic seam, so its uv shoelace area is a
    // chord artefact (the closing edge spans the whole u period at nearly
    // constant v), not a region area: it must not be used to choose the outer
    // ring. The seam rings are the boundary of the full uv domain, so they are
    // the outer boundary of a fragment and come first; interior closed rings
    // are the holes, ordered by descending |uv area|. Ties among seam rings (a
    // band has two of them) are broken by the lexicographic minimum point,
    // which does not depend on where the arrangement's walk started.
    std::stable_sort(rings.begin(), rings.end(),
                     [](const DomainRing& lhs, const DomainRing& rhs) {
                       if (lhs.seam != rhs.seam) return lhs.seam;
                       if (!lhs.seam) {
                         const double lhs_area =
                             std::fabs(signed_area2(lhs.points));
                         const double rhs_area =
                             std::fabs(signed_area2(rhs.points));
                         if (lhs_area != rhs_area) return lhs_area > rhs_area;
                         if (lhs.points.size() != rhs.points.size()) {
                           return lhs.points.size() > rhs.points.size();
                         }
                       }
                       return trim_face_lex_less(
                           trim_face_seam_only_ring_key(lhs.points),
                           trim_face_seam_only_ring_key(rhs.points));
                     });
    CurvedFragment fragment;
    fragment.reversed = false;
    for (std::size_t l = 0; l < rings.size(); ++l) {
      // Interior rings get the outer/hole orientation. A seam ring keeps its
      // walk direction so its two seam copies stay first and last; its uv area
      // is zero anyway, so the orientation carries no information.
      if (!rings[l].seam) {
        trim_face_normalize_ring(rings[l].points, l == 0);
      }
      fragment.loops_uv.push_back(std::move(rings[l].points));
    }
    if (fragment.loops_uv.empty()) continue;
    fragments.push_back(std::move(fragment));
  }

  // An interior island owns material on both sides, but the arrangement stores
  // its outer side as the unbounded record of the island ring, so only the
  // island itself becomes a fragment above. When every cut is such an island
  // the material is exactly the islands plus one complement region, and that
  // complement is a single keyhole ring: the rectangle walked clockwise, with
  // each island spliced in at one point through a doubled bridge (the two
  // bridge edges cancel in the shoelace, so the ring area is -(area of the
  // rectangle) plus the island areas). Marking it a seam ring keeps the
  // arrangement's clockwise direction: the island and its complement must
  // traverse their shared edge in the same direction to partition the chart,
  // which is the convention the island fragment test asserts.
  bool all_islands = !chains.empty();
  for (const DomainChain& chain : chains) {
    if (!chain.island) all_islands = false;
  }
  if (all_islands && !domain.periodic_v) {
    const Point2d corner{u_lo, v_lo};
    std::vector<Point2d> ring;
    ring.reserve(5u + chains.size() * 5u);
    // Start at a corner, splice every island in through a bridge that leaves and
    // returns to that same corner, then walk the rest of the rectangle. Each
    // bridge is traversed twice in opposite directions, so the two edges cancel
    // in the shoelace and the ring area is -(rectangle) plus the island areas.
    // The return to the corner before the rectangle continues is what keeps the
    // ring free of a zero-width spike: a ring that left and came back through
    // the bridge at its very end would walk the closing edge backwards.
    ring.push_back(corner);
    for (const DomainChain& chain : chains) {
      const Point2d hub = chain.points.front();
      ring.push_back(hub);
      for (std::size_t i = 1; i + 1 < chain.points.size(); ++i) {
        ring.push_back(chain.points[i]);
      }
      ring.push_back(hub);
      ring.push_back(corner);
    }
    ring.push_back(Point2d{u_lo, v_hi});
    ring.push_back(Point2d{u_hi, v_hi});
    ring.push_back(Point2d{u_hi, v_lo});
    ring.push_back(corner);
    CurvedFragment fragment;
    fragment.reversed = false;
    fragment.loops_uv.push_back(std::move(ring));
    fragments.push_back(std::move(fragment));
  }

  // Defensive backstop: the checks above refuse every input that could strip
  // all the cut geometry, so a handled synthetic domain always yields at least
  // one fragment. This guard keeps a future regression from silently returning
  // an empty success.
  if (fragments.empty()) {
    return Error{Errc::numeric_failure, "split_face_curved: no fragments",
                 face.index(), 0u};
  }

  // The fragments must tile the chart exactly: outer loops cover, holes remove.
  // The budget is well defined only for the keyhole decomposition. A seam ring
  // closes across a periodic seam, so its planar shoelace is a chord artefact
  // (see the ring sort) rather than a region area; where every cut is a closed
  // island the rings are true planar polygons. There overlapping islands come
  // apart as several material faces while the island complement still counts
  // each full island, so the covered area can exceed the chart: reject that
  // instead of returning a tiling that double-covers material.
  if (all_islands) {
    const double chart_area = (u_hi - u_lo) * (v_hi - v_lo);
    double covered = 0.0;
    for (const CurvedFragment& fragment : fragments) {
      for (std::size_t l = 0; l < fragment.loops_uv.size(); ++l) {
        const double area = std::fabs(signed_area2(fragment.loops_uv[l]));
        covered += l == 0 ? area : -area;
      }
    }
    if (std::fabs(covered - chart_area) > 1e-9 * chart_area) {
      return Error{Errc::numeric_failure,
                   "split_face_curved: fragments do not partition the chart",
                   face.index(), 0u};
    }
  }

  // An island next to another closed cut is a region with a hole: the fragment
  // that owns the island would need the island ring as a second cycle, which the
  // one-cycle-per-face extraction cannot produce. A planar containment test on
  // the emitted rings cannot see this: a pole chart closes its ring with a chord
  // from the last pole point back to the first rim point, so the ring polygon is
  // not the chart region and the island's representative is not inside it.
  // Refusing every mixed configuration is deterministic and covers exactly the
  // cases that need the second cycle, nested or not.
  if (!all_islands) {
    for (const DomainChain& chain : chains) {
      if (!chain.island) continue;
      return Error{Errc::numeric_failure,
                   "split_face_curved: nested closed cuts are not supported",
                   face.index(), 0u};
    }
  }

  std::stable_sort(fragments.begin(), fragments.end(),
                   [](const CurvedFragment& lhs,
                      const CurvedFragment& rhs) noexcept {
                     const Point2d la =
                         trim_face_seam_only_ring_key(lhs.loops_uv.front());
                     const Point2d ra =
                         trim_face_seam_only_ring_key(rhs.loops_uv.front());
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
        const Point3d point = face_uv_eval(projection, uv);
        if (!trim_face_finite3(point)) {
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

// Splits a seam-only face: a face whose uv loop is the seam slit with zero
// area, so the ordinary pipeline clips every cut away. The surface is opened
// onto the full parameter rectangle and handed to trim_face_split_domain with
// the four sides classified as seam copies (and poles collapsed). The screen
// keeps the historical rejections of the public path; the general splitter does
// not need it, because it can carry an island with its complement.
inline Result<std::vector<CurvedFragment>> trim_face_split_seam_only(
    FaceId face, const FaceUvProjection& projection,
    std::span<const FaceCut> cuts, const std::array<double, 4>& domain) {
  const Result<void> screen =
      trim_face_seam_only_screen(face, projection, cuts, domain);
  if (!screen) return screen.error();
  TrimFaceDomain descriptor;
  descriptor.bounds = domain;
  descriptor.periodic_u = projection.periodic_u();
  descriptor.periodic_v = projection.periodic_v();
  descriptor.sides = {TrimSideKind::kSeamCopy, TrimSideKind::kSeamCopy,
                      TrimSideKind::kSeamCopy, TrimSideKind::kSeamCopy};
  double pole_v = 0.0;
  if (trim_face_seam_only_pole_v(domain, pole_v)) {
    // Both v sides of the sphere chart are pole rows: each collapses to its own
    // 3D point, which is which is decided per cycle by the row it touches.
    descriptor.sides[2] = TrimSideKind::kCollapsed;
    descriptor.sides[3] = TrimSideKind::kCollapsed;
  }
  return trim_face_split_domain(face, projection, cuts, descriptor, {}, 1e-12);
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
// lie anywhere on an ordinary face, where it bounds an interior island. A
// sphere seam-only face (its boundary loop is the zero-area seam slit) has its
// parameter rectangle synthesized and its closed cuts must be latitude-like:
// the two ends must be the two copies of one physical seam point, matched with
// the same 1e-12 relative bound the parameter projection uses (so a seam point
// that drifted by rounding still counts as one point), and no other vertex may
// touch the seam. On the sphere chart a latitude cut therefore meets the seam
// at one physical point; a cap fragment that reaches a pole is emitted with
// exactly one pole point (v = -pi/2 or +pi/2, taken from the domain, not a
// projected approximation), and a pole-free band between two latitudes is
// emitted as one ring joined through the seam columns. That path inspects the u
// sides only. A closed cut that never touches the seam (an island on a
// seam-only face), a closed cut tangent to the seam, a closed cut crossing the
// seam more than once ("crosses the seam twice"), and any cut running along a
// domain side are rejected as unsupported before any arrangement work. A torus
// face whose stored outer loop covers the full parameter rectangle (area
// 4*pi^2) is split on that synthetic domain instead, where both seams are
// periodic: a latitude cut becomes whole-cycle bands and a meridian cut closes
// through the v seam. Closed cuts wholly inside that chart and stored inner
// rings are refused as unsupported, and so is a torus face whose loop covers
// less than the full chart, so a genuinely trimmed patch is never reopened and
// its trim discarded. Cut segments shorter than `tolerance` and fragments whose
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
  // Junction radius for corners where two cuts meet. The intersector traces
  // every cut over the face's parameter patch and stops where the patch ends,
  // so a trace can stop short of (or run past) the crossing with another cut by
  // a fraction of a sample step: measured 2.5e-5 rad on a unit cylinder, about
  // 4e-6 of the lateral face's uv diagonal. Those two endpoints describe one
  // geometric corner, so they are reconciled at this radius instead of the exact
  // weld radius. The constant is a parametrization error, not a model tolerance;
  // it stays far below a real gap, because a cut that stops 1e-2 short of the
  // boundary must still be reported as an open endpoint.
  const double junction_tol = std::max(64.0 * tolerance, 1e-5 * extent);

  const SurfaceGeometry& surface =
      brep.geometry().surfaces[brep.face(face).surface];

  // A canonical cone face stores its base rim plus the single apex point as the
  // outer loop. The ring closes with the implicit chart diagonal from the apex
  // back to the base seam, so its shoelace area is exactly half the chart
  // (pi * v_apex). The area alone is not enough: a band that keeps the base and
  // stops at half the apex height has the same half-chart area, because the
  // shoelace sum of a band is 2*pi*z and the chart is 2*pi*v_apex. The canonical
  // loop is therefore also required to sit entirely on the two v sides, with
  // exactly one point on the collapsed apex side. Any other loop is a genuinely
  // trimmed cone; dispatching it here would reopen the full chart and split away
  // the trim. A synthetic face with no cuts falls through to the ordinary path,
  // which returns it unchanged instead of building an empty synthetic
  // arrangement.
  if (std::holds_alternative<Cone3>(surface)) {
    const Result<detail::TrimFaceDomain> domain =
        detail::trim_face_domain(brep, face, projection, tolerance);
    if (!domain) return domain.error();
    bool synthetic = loops.size() == 1u;
    if (synthetic) {
      const double chart =
          (domain.value().bounds[1] - domain.value().bounds[0]) *
          (domain.value().bounds[3] - domain.value().bounds[2]);
      const double loop_area = std::fabs(signed_area2(loops.front()));
      synthetic = std::fabs(loop_area - 0.5 * chart) <= 1e-9 * chart;
      if (synthetic) {
        const double v_lo = domain.value().bounds[2];
        const double v_hi = domain.value().bounds[3];
        std::size_t apex_points = 0;
        bool on_rows = true;
        for (const Point2d& point : loops.front()) {
          if (detail::trim_face_near_side(point[1], v_hi, tolerance)) {
            ++apex_points;
          } else if (!detail::trim_face_near_side(point[1], v_lo, tolerance)) {
            on_rows = false;
            break;
          }
        }
        synthetic = on_rows && apex_points == 1u;
      }
    }
    if (!synthetic) {
      return Error{Errc::unsupported,
                   "split_face_curved: trimmed cone faces are not supported "
                   "yet",
                   face.index(), 0u};
    }
    if (!cuts.empty()) {
      return detail::trim_face_split_domain(
          face, projection, cuts, domain.value(),
          std::span<const std::vector<Point2d>>{}, tolerance);
    }
  }

  // A full torus face stores its whole parameter rectangle as the outer loop, so
  // its uv area is nonzero and the clipping below would treat the rectangle's
  // four sides as real face boundary: a latitude cut in the chart interior then
  // lands as an island and the arrangement splits it the wrong way. The synthetic
  // domain splitter knows the seam copies, so route a periodic-v torus through it
  // with the same domain its fragments are built from. A sphere keeps its
  // historical seam-only path below (its slit loop has zero area and its v sides
  // collapse to poles); only the torus is dispatched here.
  if (std::holds_alternative<Torus3>(surface)) {
    const Result<detail::TrimFaceDomain> domain =
        detail::trim_face_domain(brep, face, projection, tolerance);
    if (!domain) return domain.error();
    // A stored inner ring is refused by the splitter with the historical message;
    // a face with no cuts and no holes has nothing to split and falls through to
    // the ordinary path, which returns it unchanged instead of building an empty
    // synthetic arrangement.
    const std::span<const std::vector<Point2d>> holes(loops.data() + 1,
                                                      loops.size() - 1);
    if (!cuts.empty() || !holes.empty()) {
      // A closed cut that touches no synthetic side is an island: it bounds
      // material on both sides of the chart, but a fragment cannot carry the
      // complement as a hole. Refuse this scoped input explicitly before the
      // arrangement instead of leaving the generic partition guard to report it.
      bool islands_only = holes.empty() && !cuts.empty();
      for (const FaceCut& cut : cuts) {
        const std::vector<Point2d>& uv = cut.uv;
        if (!islands_only) break;
        if (uv.size() < 4 ||
            !detail::trim_face_same_xy_rel(uv.front(), uv.back(), tolerance)) {
          islands_only = false;
          break;
        }
        for (const Point2d& point : uv) {
          if (detail::trim_face_point_on_domain_side(point, domain.value(),
                                                     tolerance)) {
            islands_only = false;
            break;
          }
        }
      }
      if (islands_only) {
        return Error{Errc::unsupported,
                     "split_face_curved: closed cuts on a periodic face are "
                     "not supported yet",
                     face.index(), 0u};
      }
      // Only the full parameter rectangle is a synthetic chart. A stored loop
      // that covers less of the rectangle is a trimmed torus patch whose
      // removed material is real: dispatching it here would reopen the full
      // chart and split away the trim, so refuse it instead. The bound is
      // relative to the chart area so a canonical torus, whose loop area is
      // exactly 4*pi^2, keeps passing.
      if (holes.empty()) {
        const double chart_area =
            (domain.value().bounds[1] - domain.value().bounds[0]) *
            (domain.value().bounds[3] - domain.value().bounds[2]);
        const double loop_area = std::fabs(signed_area2(loops.front()));
        if (loops.size() != 1u ||
            std::fabs(loop_area - chart_area) > 1e-9 * chart_area) {
          return Error{Errc::unsupported,
                       "split_face_curved: trimmed torus faces are not "
                       "supported yet",
                       face.index(), 0u};
        }
      }
      return detail::trim_face_split_domain(face, projection, cuts,
                                            domain.value(), holes, tolerance);
    }
  }

  // A seam-only face has no loop-derived uv area: its loop is the seam slit
  // itself, so the clipping below would discard every cut and return no
  // fragments. On a sphere the slit's uv area is a rounding artefact (a chart
  // recovered off axis has area around 1e-33, not exactly zero); the relative
  // bound below recognizes it without touching the torus path, whose loop area
  // stays exactly zero and whose near-tangent cuts are not covered by that
  // reasoning. Open the surface onto its full parameter rectangle instead.
  const double loop_area = std::fabs(signed_area2(loops.front()));
  const bool sphere_slit =
      std::holds_alternative<Sphere3>(surface) &&
      loop_area <= 1e-12 * extent * extent;
  if (loop_area == 0.0 || sphere_slit) {
    // The seam-only splitter splits the outer ring only. A synthetic-domain face
    // that already carries a stored inner ring (a band with a hole) would lose
    // that hole silently, so refuse it explicitly until the splitter carries
    // multiple cycles.
    if (sphere_slit && loops.size() > 1) {
      return Error{Errc::unsupported,
                   "split_face_curved: stored holes are not supported yet",
                   face.index(), static_cast<std::uint32_t>(loops.size() - 1)};
    }
    const std::optional<std::array<double, 4>> seam_domain =
        detail::trim_face_seam_only_domain(surface);
    if (seam_domain.has_value()) {
      return detail::trim_face_split_seam_only(face, projection, cuts,
                                               seam_domain.value());
    }
  }

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
  // A cut whose points all coincide carries no direction. Reject it before the
  // clipping step below: clipping would drop it entirely and the face would be
  // returned unsplit instead of reporting the degenerate input. The test is
  // exact, matching the arrangement's duplicate-segment test, so a tiny but
  // genuine closed island still reaches the small-fragment filter.
  for (std::size_t c = 0; c < cut_points.size(); ++c) {
    bool degenerate = true;
    for (const Point2d& point : cut_points[c]) {
      if (!detail::trim_face_same_xy(point, cut_points[c].front())) {
        degenerate = false;
        break;
      }
    }
    if (degenerate) {
      return Error{Errc::degenerate,
                   "split_face_curved: cut collapses to a point", face.index(),
                   static_cast<std::uint32_t>(c)};
    }
  }

  // Exact arc junctions. A cut segment that meets a boundary arc is split at
  // the analytic point on the circle: the loop gains a vertex whose arc
  // reference carries the exact parameter, and the cut gains the same projected
  // vertex, so the arrangement never approximates the crossing with the
  // inscribed chord.
  std::vector<std::vector<std::optional<FaceUvArcRef>>> loop_arcs(loops.size());
  for (std::size_t l = 0; l < loops.size(); ++l) {
    loop_arcs[l] = uv_loops[l].arcs;
    loop_arcs[l].resize(loops[l].size());
  }
  struct ArcJunction {
    double along = 0.0;
    double t = 0.0;
    Point2d uv{};
  };
  struct CutJunction {
    std::size_t cut = 0;
    std::size_t edge = 0;
    double s = 0.0;
    Point2d uv{};
  };
  struct CutRewrite {
    std::size_t cut = 0;
    std::size_t vertex = 0;
    Point2d uv{};
  };
  std::vector<std::vector<std::vector<ArcJunction>>> loop_junctions(loops.size());
  for (std::size_t l = 0; l < loops.size(); ++l) {
    loop_junctions[l].resize(loops[l].size());
  }
  std::vector<CutJunction> cut_junctions;
  std::vector<CutRewrite> cut_rewrites;
  std::vector<std::vector<Point3d>> cut_points3(cut_points.size());
  for (std::size_t c = 0; c < cut_points.size(); ++c) {
    cut_points3[c].reserve(cut_points[c].size());
    for (const Point2d& point : cut_points[c]) {
      cut_points3[c].push_back(detail::face_uv_eval(projection, point));
    }
  }
  for (std::size_t l = 0; l < loops.size(); ++l) {
    const std::size_t count = loops[l].size();
    for (std::size_t i = 0; i < count; ++i) {
      if (!loop_arcs[l][i]) continue;
      const CircleCurve3& circle = loop_arcs[l][i]->circle;
      const double sweep = loop_arcs[l][i]->t1 - loop_arcs[l][i]->t0;
      if (!(std::fabs(sweep) > 0.0)) continue;
      const double direction = sweep > 0.0 ? 1.0 : -1.0;
      const double length = std::fabs(sweep);
      const double slack =
          std::max(1e-12, tolerance / std::max(circle.radius, 1e-12));
      for (std::size_t c = 0; c < cut_points.size(); ++c) {
        for (std::size_t e = 0; e + 1 < cut_points[c].size(); ++e) {
          const std::vector<detail::TrimFaceCircleHit> hits =
              detail::trim_face_circle_hits(circle, cut_points3[c][e],
                                            cut_points3[c][e + 1], tolerance);
          for (const detail::TrimFaceCircleHit& hit : hits) {
            double offset = hit.t - loop_arcs[l][i]->t0;
            if (offset > detail::kFaceUvPi) {
              offset -= detail::kTrimFaceTwoPi;
            } else if (offset < -detail::kFaceUvPi) {
              offset += detail::kTrimFaceTwoPi;
            }
            const double along = offset * direction;
            if (!(along >= 0.0) || !(along <= length)) continue;
            if (along < slack || along > length - slack) {
              // The contact lands on an existing boundary vertex; the ordinary
              // endpoint snap already puts the cut exactly there.
              continue;
            }
            const Result<Point2d> projected = detail::face_uv_project(
                projection, hit.point, loops[l][i][0], true);
            if (!projected) return projected.error();
            Point2d point = projected.value();
            if (projection.periodic_u()) {
              point[0] += detail::kTrimFaceTwoPi *
                          std::round((loops[l][i][0] - point[0]) /
                                     detail::kTrimFaceTwoPi);
            }
            if (projection.periodic_v()) {
              point[1] += detail::kTrimFaceTwoPi *
                          std::round((loops[l][i][1] - point[1]) /
                                     detail::kTrimFaceTwoPi);
            }
            const double t = loop_arcs[l][i]->t0 + direction * along;
            loop_junctions[l][i].push_back(ArcJunction{along, t, point});
            if (hit.s > 0.0 && hit.s < 1.0) {
              cut_junctions.push_back(CutJunction{c, e, hit.s, point});
            } else if (hit.s <= 0.0) {
              cut_rewrites.push_back(CutRewrite{c, e, point});
            } else {
              cut_rewrites.push_back(CutRewrite{c, e + 1, point});
            }
          }
        }
      }
    }
  }
  {
    std::vector<std::vector<Point2d>> rebuilt_loops(loops.size());
    std::vector<std::vector<std::optional<FaceUvArcRef>>> rebuilt_arcs(
        loops.size());
    for (std::size_t l = 0; l < loops.size(); ++l) {
      for (std::size_t i = 0; i < loops[l].size(); ++i) {
        rebuilt_loops[l].push_back(loops[l][i]);
        std::vector<ArcJunction>& junctions = loop_junctions[l][i];
        std::sort(junctions.begin(), junctions.end(),
                  [](const ArcJunction& lhs, const ArcJunction& rhs) {
                    return lhs.along < rhs.along;
                  });
        // Each junction splits one segment reference into two; the traversal
        // direction of every piece is carried by its own signed interval.
        double start_t = loop_arcs[l][i] ? loop_arcs[l][i]->t0 : 0.0;
        for (const ArcJunction& junction : junctions) {
          if (detail::trim_face_same_xy(rebuilt_loops[l].back(), junction.uv)) {
            continue;
          }
          if (loop_arcs[l][i]) {
            rebuilt_arcs[l].push_back(FaceUvArcRef{loop_arcs[l][i]->circle,
                                                   start_t, junction.t});
          } else {
            rebuilt_arcs[l].push_back(std::nullopt);
          }
          rebuilt_loops[l].push_back(junction.uv);
          start_t = junction.t;
        }
        if (loop_arcs[l][i]) {
          rebuilt_arcs[l].push_back(FaceUvArcRef{loop_arcs[l][i]->circle,
                                                 start_t,
                                                 loop_arcs[l][i]->t1});
        } else {
          rebuilt_arcs[l].push_back(std::nullopt);
        }
      }
    }
    loops.swap(rebuilt_loops);
    loop_arcs.swap(rebuilt_arcs);
  }
  for (const CutRewrite& rewrite : cut_rewrites) {
    if (rewrite.vertex < cut_points[rewrite.cut].size()) {
      cut_points[rewrite.cut][rewrite.vertex] = rewrite.uv;
    }
  }
  if (!cut_junctions.empty()) {
    std::sort(cut_junctions.begin(), cut_junctions.end(),
              [](const CutJunction& lhs, const CutJunction& rhs) {
                if (lhs.cut != rhs.cut) return lhs.cut < rhs.cut;
                if (lhs.edge != rhs.edge) return lhs.edge < rhs.edge;
                return lhs.s < rhs.s;
              });
    for (std::size_t c = 0; c < cut_points.size(); ++c) {
      const std::size_t begin = static_cast<std::size_t>(
          std::lower_bound(cut_junctions.begin(), cut_junctions.end(), c,
                           [](const CutJunction& lhs, std::size_t value) {
                             return lhs.cut < value;
                           }) -
          cut_junctions.begin());
      std::size_t next = begin;
      if (next == cut_junctions.size() || cut_junctions[next].cut != c) {
        continue;
      }
      std::vector<Point2d> rebuilt;
      rebuilt.reserve(cut_points[c].size() + 8);
      rebuilt.push_back(cut_points[c].front());
      for (std::size_t e = 0; e + 1 < cut_points[c].size(); ++e) {
        while (next < cut_junctions.size() && cut_junctions[next].cut == c &&
               cut_junctions[next].edge == e) {
          if (!detail::trim_face_same_xy(rebuilt.back(),
                                         cut_junctions[next].uv)) {
            rebuilt.push_back(cut_junctions[next].uv);
          }
          ++next;
        }
        if (!detail::trim_face_same_xy(rebuilt.back(), cut_points[c][e + 1])) {
          rebuilt.push_back(cut_points[c][e + 1]);
        }
      }
      if (cut_closed[c] != 0) rebuilt.back() = rebuilt.front();
      cut_points[c].swap(rebuilt);
    }
  }

  // The intersector traces a cut over the whole parameter patch, which for a
  // trimmed face reaches past the real boundary (a plane through a disk cap is
  // the usual case). Keep only the part of each chain that lies on the face; the
  // arc junctions inserted above are exact boundary points, so a chain that
  // leaves through an arc keeps an exact crossing. Nothing changes for callers
  // whose cuts already end on the boundary.
  {
    std::vector<std::vector<Point2d>> clipped_points;
    std::vector<char> clipped_closed;
    clipped_points.reserve(cut_points.size());
    clipped_closed.reserve(cut_points.size());
    for (std::size_t c = 0; c < cut_points.size(); ++c) {
      std::vector<std::vector<Point2d>> pieces = detail::trim_face_clip_cut(
          loops, cut_points[c], cut_closed[c] != 0, snap_radius);
      for (std::vector<Point2d>& piece : pieces) {
        const bool piece_closed =
            piece.size() >= 3 && detail::trim_face_same_xy(piece.front(), piece.back());
        clipped_points.push_back(std::move(piece));
        clipped_closed.push_back(piece_closed ? 1 : 0);
      }
    }
    cut_points.swap(clipped_points);
    cut_closed.swap(clipped_closed);
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
    bool front = false;
  };
  std::vector<Endpoint> endpoints;
  for (std::size_t c = 0; c < cut_points.size(); ++c) {
    if (cut_closed[c] != 0) continue;
    for (const bool front : {true, false}) {
      const Point2d& raw = front ? cut_points[c].front() : cut_points[c].back();
      detail::TrimFaceSnap snap = detail::trim_face_find_snap(loops, raw, junction_tol);
      // A chord projection onto an arc boundary is not exact: the chord lies
      // inside the arc by its sagitta, and the arrangement only accepts boundary
      // vertices that are loop points or exact circle points. Project the raw
      // endpoint onto the arc circle itself instead, keeping the angle of the
      // trace (radial projection), and fall back to reconciling the endpoint
      // against the other cuts when it does not land on the arc's span.
      if (snap.boundary && !snap.vertex && snap.ring < loop_arcs.size() &&
          snap.edge < loop_arcs[snap.ring].size() && loop_arcs[snap.ring][snap.edge]) {
        bool projected = false;
        const FaceUvArcRef& arc = *loop_arcs[snap.ring][snap.edge];
        const CircleCurve3& circle = arc.circle;
        const double sweep = arc.t1 - arc.t0;
        const Result<Point3d> lifted = detail::face_uv_eval(projection, raw);
        if (lifted.has_value() && std::fabs(sweep) > 0.0) {
          const Vec3d normal = cross(circle.x_axis, circle.y_axis);
          const Vec3d radial = lifted.value() - circle.center;
          const double height = dot(radial, normal);
          const Vec3d in_plane{radial[0] - normal[0] * height,
                               radial[1] - normal[1] * height,
                               radial[2] - normal[2] * height};
          const double length = std::sqrt(norm2(in_plane));
          if (length > 0.0) {
            const Point3d on_arc{
                circle.center[0] + in_plane[0] * (circle.radius / length),
                circle.center[1] + in_plane[1] * (circle.radius / length),
                circle.center[2] + in_plane[2] * (circle.radius / length)};
            double offset =
                detail::face_uv_circle_angle(circle, on_arc) - arc.t0;
            if (offset > detail::kFaceUvPi) {
              offset -= detail::kTrimFaceTwoPi;
            } else if (offset < -detail::kFaceUvPi) {
              offset += detail::kTrimFaceTwoPi;
            }
            const double along = sweep > 0.0 ? offset : -offset;
            if (along >= 0.0 && along <= std::fabs(sweep)) {
              const Result<Point2d> uv = detail::face_uv_project(
                  projection, on_arc, raw[0], true);
              if (uv.has_value()) {
                snap.point = uv.value();
                projected = true;
              }
            }
          }
        }
        if (!projected) {
          snap.boundary = false;
          snap.point = raw;
        }
      }
      endpoints.push_back(Endpoint{raw, snap, Point2d{}, c, front});
    }
  }
  const std::size_t endpoint_count = endpoints.size();
  std::vector<std::size_t> parent(endpoint_count);
  for (std::size_t i = 0; i < endpoint_count; ++i) parent[i] = i;
  for (std::size_t i = 0; i < endpoint_count; ++i) {
    for (std::size_t j = i + 1; j < endpoint_count; ++j) {
      const double distance = std::hypot(
          endpoints[i].raw[0] - endpoints[j].raw[0],
          endpoints[i].raw[1] - endpoints[j].raw[1]);
      if (!(distance <= junction_tol)) continue;
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
        for (std::size_t host = 0; host < cut_points.size(); ++host) {
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
            if (!(distance <= junction_tol)) continue;
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
      } else {
        // Two open ends meet at a corner that neither trace samples exactly.
        // Weld them onto the intersection of their end segments: welding onto
        // either raw point would drag the other chain off its own exact plane
        // (the vertical cut carries the exact coordinate, the patch-limited arc
        // trace is the one that stops short), leaving the two faces a trace step
        // apart when the result is assembled.
        std::size_t first = endpoint_count;
        std::size_t second = endpoint_count;
        double closest = 0.0;
        for (std::size_t j = 0; j < endpoint_count; ++j) {
          if (detail::trim_face_find_root(parent, j) != root) continue;
          for (std::size_t k = j + 1; k < endpoint_count; ++k) {
            if (detail::trim_face_find_root(parent, k) != root) continue;
            if (endpoints[j].cut == endpoints[k].cut) continue;
            const double distance =
                std::hypot(endpoints[j].raw[0] - endpoints[k].raw[0],
                           endpoints[j].raw[1] - endpoints[k].raw[1]);
            if (first == endpoint_count || distance < closest) {
              closest = distance;
              first = j;
              second = k;
            }
          }
        }
        if (first != endpoint_count) {
          const Endpoint& lhs = endpoints[first];
          const Endpoint& rhs = endpoints[second];
          const auto end_direction = [&cut_points](const Endpoint& endpoint) -> Vec2d {
            const std::vector<Point2d>& chain = cut_points[endpoint.cut];
            if (endpoint.front) {
              return Vec2d{chain[1][0] - chain[0][0], chain[1][1] - chain[0][1]};
            }
            const std::size_t last = chain.size() - 1;
            return Vec2d{chain[last - 1][0] - chain[last][0],
                         chain[last - 1][1] - chain[last][1]};
          };
          const Vec2d lhs_edge = end_direction(lhs);
          const Vec2d rhs_edge = end_direction(rhs);
          const double lhs_length = std::hypot(lhs_edge[0], lhs_edge[1]);
          const double rhs_length = std::hypot(rhs_edge[0], rhs_edge[1]);
          // The acceptance bound below is a uv displacement, so the end
          // directions are normalised first: a bound on the raw parameter would
          // scale with the (arbitrary) length of the end segment and let a long
          // final segment move the corner further than junction_tol.
          if (lhs_length > 0.0 && rhs_length > 0.0) {
            const Vec2d lhs_unit{lhs_edge[0] / lhs_length,
                                 lhs_edge[1] / lhs_length};
            const Vec2d rhs_unit{rhs_edge[0] / rhs_length,
                                 rhs_edge[1] / rhs_length};
            const double determinant =
                lhs_unit[0] * rhs_unit[1] - lhs_unit[1] * rhs_unit[0];
            if (std::fabs(determinant) > 1e-12) {
              const Vec2d delta{rhs.raw[0] - lhs.raw[0],
                                rhs.raw[1] - lhs.raw[1]};
              const double along_lhs =
                  (delta[0] * rhs_unit[1] - delta[1] * rhs_unit[0]) /
                  determinant;
              const double along_rhs =
                  (delta[0] * lhs_unit[1] - delta[1] * lhs_unit[0]) /
                  determinant;
              if (std::fabs(along_lhs) <= 2.0 * junction_tol &&
                  std::fabs(along_rhs) <= 2.0 * junction_tol) {
                final = Point2d{lhs.raw[0] + along_lhs * lhs_unit[0],
                                lhs.raw[1] + along_lhs * lhs_unit[1]};
              }
            }
          }
        }
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
    std::vector<std::vector<CutInsert>> per_cut(cut_points.size());
    for (const CutInsert& insert : cut_inserts) {
      per_cut[insert.cut].push_back(insert);
    }
    for (std::size_t c = 0; c < cut_points.size(); ++c) {
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
  // ring share one exact vertex. The inserts are applied to the loop rings
  // themselves (not a copy): a snapped point that lands on an arc segment has
  // to become a loop point of its circle, or the fragment mapping, which
  // matches ring vertices against the loop points exactly, rejects the vertex.
  std::vector<std::vector<detail::TrimFaceInsert>> per_ring(loops.size());
  std::vector<char> inserted(endpoint_count, 0);
  for (std::size_t i = 0; i < endpoint_count; ++i) {
    const std::size_t root = detail::trim_face_find_root(parent, i);
    if (inserted[root] != 0) continue;
    inserted[root] = 1;
    const detail::TrimFaceSnap& snap = canonical[root];
    if (!snap.boundary || snap.vertex) continue;
    bool already_present = false;
    for (const Point2d& vertex : loops[snap.ring]) {
      if (detail::trim_face_same_xy(vertex, snap.point)) {
        already_present = true;
        break;
      }
    }
    if (!already_present) {
      per_ring[snap.ring].push_back({snap.ring, snap.edge, snap.point});
    }
  }
  for (std::size_t l = 0; l < loops.size(); ++l) {
    if (per_ring[l].empty()) continue;
    detail::trim_face_apply_inserts_with_arcs(projection, loops[l], loop_arcs[l],
                                              std::move(per_ring[l]));
  }
  std::vector<std::vector<Point2d>> rings = loops;

  // Rewrite the open cut chains with the welded endpoint coordinates.
  std::size_t endpoint_index = 0;
  for (std::size_t c = 0; c < cut_points.size(); ++c) {
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
  for (std::size_t c = 0; c < cut_points.size(); ++c) {
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
  segments.reserve(segment_total + cut_points.size() * 8u);
  for (std::size_t r = 0; r < rings.size(); ++r) {
    const std::vector<Point2d>& ring = rings[r];
    for (std::size_t i = 0; i < ring.size(); ++i) {
      const Point2d& a = ring[i];
      const Point2d& b = ring[(i + 1) % ring.size()];
      if (detail::trim_face_same_xy(a, b)) continue;
      segments.push_back(ArrangementSegment{a, b, 0u});
    }
  }
  for (std::size_t c = 0; c < cut_points.size(); ++c) {
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

  // Map the arrangement rings back to the boundary loops. The arrangement
  // copies input vertices bit for bit, so a vertex that came from a loop point
  // compares exactly equal to it; a vertex that is still on an arc but is not
  // one of its points means the junction insertion was missed.
  std::vector<std::vector<Point3d>> loop_points3(loops.size());
  for (std::size_t l = 0; l < loops.size(); ++l) {
    loop_points3[l].reserve(loops[l].size());
    for (const Point2d& uv : loops[l]) {
      loop_points3[l].push_back(detail::face_uv_eval(projection, uv));
    }
  }
  auto map_ring = [&](const std::vector<Point2d>& ring)
      -> Result<std::vector<std::optional<detail::TrimFaceRingRef>>> {
    std::vector<std::optional<detail::TrimFaceRingRef>> refs(ring.size());
    for (std::size_t i = 0; i < ring.size(); ++i) {
      const Point3d point = detail::face_uv_eval(projection, ring[i]);
      bool found = false;
      for (std::size_t l = 0; l < loop_points3.size() && !found; ++l) {
        for (std::size_t j = 0; j < loop_points3[l].size(); ++j) {
          const Point3d& candidate = loop_points3[l][j];
          if (candidate[0] == point[0] && candidate[1] == point[1] &&
              candidate[2] == point[2]) {
            refs[i] = detail::TrimFaceRingRef{l, j};
            found = true;
            break;
          }
        }
      }
      if (found) continue;
      for (std::size_t l = 0; l < loops.size(); ++l) {
        const std::size_t count = loops[l].size();
        for (std::size_t j = 0; j < count; ++j) {
          const auto& arc = loop_arcs[l][j];
          if (!arc) continue;
          const CircleCurve3& circle = arc->circle;
          const double sweep = arc->t1 - arc->t0;
          if (!(std::fabs(sweep) > 0.0)) continue;
          const double direction = sweep > 0.0 ? 1.0 : -1.0;
          const double length = std::fabs(sweep);
          const double eps =
              std::max(tolerance, 1e-12 * std::max(1.0, circle.radius));
          const Vec3d radial = point - circle.center;
          const Vec3d normal = cross(circle.x_axis, circle.y_axis);
          if (!(std::fabs(dot(radial, normal)) <= eps)) continue;
          if (!(std::fabs(norm2(radial) - circle.radius) <= eps)) continue;
          double offset =
              detail::face_uv_circle_angle(circle, point) - arc->t0;
          if (offset > detail::kFaceUvPi) {
            offset -= detail::kTrimFaceTwoPi;
          } else if (offset < -detail::kFaceUvPi) {
            offset += detail::kTrimFaceTwoPi;
          }
          const double along = offset * direction;
          const double slack =
              std::max(1e-12, eps / std::max(circle.radius, 1e-12));
          if (along > slack && along < length - slack) {
            // Not covered by a regression test: firing this needs a ring vertex
            // that sits strictly inside an arc segment yet is not bitwise equal
            // to a loop point, which the current arrangement never produces
            // (computed crossings that end up on the arc are inserted into the
            // loop beforehand, so they map exactly).
            return Error{Errc::numeric_failure,
                         "split_face_curved: fragment vertex left its loop",
                         face.index(), 0u};
          }
        }
      }
    }
    return refs;
  };
  auto span_ring = [&](const std::vector<std::optional<detail::TrimFaceRingRef>>&
                           refs,
                       std::size_t count) {
    std::vector<std::optional<CurvedArcSpan>> spans(count);
    for (std::size_t i = 0; i < count; ++i) {
      spans[i] = detail::trim_face_ring_span(refs, loops, loop_arcs, i);
    }
    return spans;
  };

  std::vector<CurvedFragment> fragments;
  const double area_floor = tolerance * tolerance;
  for (const detail::TrimFaceCycle& record : cycles) {
    if (!record.material || !record.has_rep) continue;

    CurvedFragment fragment;
    fragment.reversed = false;
    std::vector<Point2d> outer = record.ring;
    Result<std::vector<std::optional<detail::TrimFaceRingRef>>> outer_refs =
        map_ring(outer);
    if (!outer_refs) return outer_refs.error();
    detail::trim_face_normalize_ring_refs(outer, outer_refs.value(), true);
    if (outer.size() < 3) continue;
    const double outer_area = std::fabs(signed_area2(outer));
    if (!(outer_area > 0.0)) continue;
    double area_uv = outer_area;
    struct HolePiece {
      std::vector<Point2d> ring;
      std::vector<std::optional<detail::TrimFaceRingRef>> refs;
    };
    std::vector<HolePiece> holes;
    holes.reserve(record.children.size());
    for (const std::size_t child : record.children) {
      HolePiece hole;
      hole.ring = cycles[child].ring;
      Result<std::vector<std::optional<detail::TrimFaceRingRef>>> hole_refs =
          map_ring(hole.ring);
      if (!hole_refs) return hole_refs.error();
      hole.refs = std::move(hole_refs.value());
      detail::trim_face_normalize_ring_refs(hole.ring, hole.refs, false);
      if (hole.ring.size() < 3) continue;
      const double hole_area = std::fabs(signed_area2(hole.ring));
      if (!(hole_area > 0.0)) continue;
      area_uv -= hole_area;
      holes.push_back(std::move(hole));
    }
    // Holes are ordered by area (descending) and then lexicographically.
    std::sort(holes.begin(), holes.end(),
              [](const HolePiece& lhs, const HolePiece& rhs) noexcept {
                const double lhs_area = std::fabs(signed_area2(lhs.ring));
                const double rhs_area = std::fabs(signed_area2(rhs.ring));
                if (lhs_area != rhs_area) return lhs_area > rhs_area;
                return detail::trim_face_lex_less(lhs.ring.front(),
                                                  rhs.ring.front());
              });
    if (!(area_uv > 0.0)) continue;

    // A fragment suppressed by the area floor (for example a tiny island) does
    // not remove its cycle from the arrangement: the enclosing fragment keeps
    // the island outline as a hole, so the emitted parent ring has a hole where
    // no fragment was emitted. The hole is left in place rather than patched,
    // because its outline is still an exact piece of the input subdivision.
    const double scale = detail::trim_face_metric_scale(projection, record.rep);
    if (!(area_uv * scale > area_floor)) continue;

    // Hole rings get exactly the same span treatment as the outer ring. The
    // current tests only use straight holes (no engaged arc), so an arc running
    // along a hole boundary is supported by construction but not covered by a
    // regression test.
    const std::size_t outer_count = outer.size();
    fragment.arc_spans.push_back(span_ring(outer_refs.value(), outer_count));
    fragment.loops_uv.push_back(std::move(outer));
    for (HolePiece& hole : holes) {
      const std::size_t hole_count = hole.ring.size();
      fragment.arc_spans.push_back(span_ring(hole.refs, hole_count));
      fragment.loops_uv.push_back(std::move(hole.ring));
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
