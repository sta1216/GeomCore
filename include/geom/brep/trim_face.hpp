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
    const std::optional<FaceUvArcRef> source =
        i < arcs.size() ? arcs[i] : std::nullopt;
    double start_t = source ? source->t0 : 0.0;
    rebuilt.push_back(ring[i]);
    for (const std::pair<double, Point2d>& entry : extra[i]) {
      double angle = start_t;
      if (source && split_angle(*source, entry.second, angle)) {
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
