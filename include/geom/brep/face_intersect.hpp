#pragma once

// Cut segments where two planar B-rep faces meet.
//
// Two planes either cross along one line or are parallel. The crossing line is
// clipped to each face by projecting it into the face frame and keeping the
// intervals whose midpoint lies strictly inside the face. Parallel faces that
// lie on one plane are intersected as 2D polygons and every result ring edge is
// returned. Curved faces are unsupported.
//
// The tolerance decides plane membership and coplanarity; no other epsilon is
// introduced beyond the parallel guard below.

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <span>
#include <vector>

#include <geom/brep/face_query.hpp>
#include <geom/brep/plane_frame.hpp>
#include <geom/brep/topology.hpp>
#include <geom/cgeom/boolean2.hpp>
#include <geom/cgeom/polygon2.hpp>
#include <geom/core/error.hpp>
#include <geom/core/result.hpp>
#include <geom/primitives/point.hpp>
#include <geom/surfaces/plane.hpp>

namespace geom {
inline namespace v1 {

// One cut segment, with endpoints in 3D world coordinates.
struct Segment3 {
  Point3d a;
  Point3d b;
};

// The cut of one face pair. `coplanar` distinguishes a true crossing with no
// visible cut (empty segments) from two faces that lie on the same plane.
struct FacePairCut {
  std::vector<Segment3> segments;
  bool coplanar = false;
};

namespace detail {

// Plane normals whose cross product is this short are treated as parallel.
// This is an angle guard, not a topology decision.
constexpr double kFaceIntersectParallelEpsilon = 1e-12;

inline double face_intersect_cross2(const Vec2d& a, const Vec2d& b) noexcept {
  return a[0] * b[1] - a[1] * b[0];
}

inline bool face_intersect_point_less(const Point3d& a, const Point3d& b) noexcept {
  if (a[0] != b[0]) return a[0] < b[0];
  if (a[1] != b[1]) return a[1] < b[1];
  return a[2] < b[2];
}

// Endpoints are ordered lexicographically so the result is deterministic.
inline Segment3 face_intersect_ordered(const Point3d& a, const Point3d& b) noexcept {
  if (face_intersect_point_less(b, a)) return Segment3{b, a};
  return Segment3{a, b};
}

// One parameter interval on the infinite crossing line x = p0 + t * u, where u
// is the unit line direction.
struct FaceIntersectInterval {
  double lo = 0.0;
  double hi = 0.0;
};

// Builds a Polygon2d from rings that are already projected into one frame.
// Ring orientation is left as stored: boolean2 normalizes outer and hole rings
// itself, so a pre-pass here would duplicate that rule.
inline Polygon2d face_intersect_polygon(
    const std::vector<std::vector<Point2d>>& loops) {
  Polygon2d polygon;
  polygon.ring_starts.push_back(0u);
  for (const std::vector<Point2d>& ring : loops) {
    polygon.vertices.insert(polygon.vertices.end(), ring.begin(), ring.end());
    polygon.ring_starts.push_back(static_cast<std::uint32_t>(polygon.vertices.size()));
  }
  return polygon;
}

// Clips the crossing line to one face. The line is projected into the face
// frame; every boundary edge contributes the parameters where it crosses the
// line, including both endpoints of an edge that lies along the line. Sorted
// crossing parameters bound the candidate intervals, and the midpoint of each
// one is tested against the face interior.
inline Result<std::vector<FaceIntersectInterval>> face_intersect_clip_line(
    const Brep& brep, FaceId id, const PlaneFrame& frame, const Point3d& p0,
    const Vec3d& direction, double tol) {
  const Result<std::vector<std::vector<Point2d>>> loops =
      face_loops_2d(brep, id, frame);
  if (!loops) return loops.error();

  const Point2d q0 = plane_frame_to_2d(frame, p0);
  const Vec2d line_direction{dot(direction, frame.u), dot(direction, frame.v)};
  const double line_length2 = dot(line_direction, line_direction);

  std::vector<double> crossings;
  for (const std::vector<Point2d>& ring : loops.value()) {
    if (ring.size() < 2) continue;
    for (std::size_t i = 0; i < ring.size(); ++i) {
      const Point2d& p1 = ring[i];
      const Point2d& p2 = ring[(i + 1) % ring.size()];
      const Vec2d edge = p2 - p1;
      const Vec2d offset = p1 - q0;
      const double denom = face_intersect_cross2(line_direction, edge);
      if (denom != 0.0) {
        const double s = face_intersect_cross2(offset, edge) / denom;
        const double k = face_intersect_cross2(offset, line_direction) / denom;
        if (k >= 0.0 && k <= 1.0) crossings.push_back(s);
      } else if (face_intersect_cross2(offset, line_direction) == 0.0) {
        // The edge is collinear with the line. Its endpoints bound every part
        // of the edge that can belong to a clipped interval.
        if (line_length2 > 0.0) {
          crossings.push_back(dot(p1 - q0, line_direction) / line_length2);
          crossings.push_back(dot(p2 - q0, line_direction) / line_length2);
        }
      }
    }
  }

  std::sort(crossings.begin(), crossings.end());
  std::vector<FaceIntersectInterval> intervals;
  for (std::size_t i = 1; i < crossings.size(); ++i) {
    const double s0 = crossings[i - 1];
    const double s1 = crossings[i];
    if (!(s1 > s0)) continue;
    const double mid = 0.5 * (s0 + s1);
    const Point3d midpoint = p0 + direction * mid;
    const Result<bool> inside = point_in_face(brep, id, midpoint, tol);
    if (!inside) return inside.error();
    if (!inside.value()) continue;
    const Point3d x0 = p0 + direction * s0;
    const Point3d x1 = p0 + direction * s1;
    const double t0 = dot(x0 - p0, direction);
    const double t1 = dot(x1 - p0, direction);
    FaceIntersectInterval interval;
    interval.lo = t0 < t1 ? t0 : t1;
    interval.hi = t0 < t1 ? t1 : t0;
    intervals.push_back(interval);
  }
  return intervals;
}

// Sorts and joins touching or overlapping intervals. A crossing that passes
// exactly through a boundary vertex produces two intervals that share an
// endpoint, and the result must still be one segment.
inline std::vector<FaceIntersectInterval> face_intersect_merge_intervals(
    std::vector<FaceIntersectInterval> intervals) {
  std::sort(intervals.begin(), intervals.end(),
            [](const FaceIntersectInterval& a, const FaceIntersectInterval& b) {
              return a.lo < b.lo;
            });
  std::vector<FaceIntersectInterval> merged;
  for (const FaceIntersectInterval& interval : intervals) {
    if (merged.empty() || interval.lo > merged.back().hi) {
      merged.push_back(interval);
    } else if (interval.hi > merged.back().hi) {
      merged.back().hi = interval.hi;
    }
  }
  return merged;
}

// Every ring edge of every result polygon, wrapped around and ordered.
inline void face_intersect_emit_polygons(const PolygonSet& set,
                                         const PlaneFrame& frame,
                                         std::vector<Segment3>& segments) {
  for (const Polygon2d& polygon : set.polygons) {
    for (std::uint32_t r = 0; r < polygon.ring_count(); ++r) {
      const std::span<const Point2d> ring = polygon.ring(r);
      if (ring.size() < 2) continue;
      for (std::size_t i = 0; i < ring.size(); ++i) {
        const Point3d a = plane_frame_to_3d(frame, ring[i]);
        const Point3d b = plane_frame_to_3d(frame, ring[(i + 1) % ring.size()]);
        segments.push_back(face_intersect_ordered(a, b));
      }
    }
  }
}

}  // namespace detail

// Returns the cut segments of two planar faces. Both faces must be alive and
// backed by planes; a curved face is unsupported. A pair of coplanar faces is
// handled as a 2D polygon intersection and reports coplanar = true even when
// the intersection is empty.
inline Result<FacePairCut> intersect_planar_faces(const Brep& a, FaceId fa,
                                                  const Brep& b, FaceId fb,
                                                  double tol) {
  const Result<bool> face_a_check = detail::face_query_check_face(a, fa);
  if (!face_a_check) return face_a_check.error();
  const Result<bool> face_b_check = detail::face_query_check_face(b, fb);
  if (!face_b_check) return face_b_check.error();
  const Result<bool> tolerance_check = detail::face_query_check_tolerance(tol);
  if (!tolerance_check) return tolerance_check.error();

  const Result<Plane3> plane_a = face_plane(a, fa);
  if (!plane_a) return plane_a.error();
  const Result<Plane3> plane_b = face_plane(b, fb);
  if (!plane_b) return plane_b.error();
  const Result<PlaneFrame> frame_a_result = make_plane_frame(plane_a.value());
  if (!frame_a_result) return frame_a_result.error();
  const Result<PlaneFrame> frame_b_result = make_plane_frame(plane_b.value());
  if (!frame_b_result) return frame_b_result.error();
  const PlaneFrame& frame_a = frame_a_result.value();
  const PlaneFrame& frame_b = frame_b_result.value();

  const Vec3d& normal_a = frame_a.normal;
  const Vec3d& normal_b = frame_b.normal;
  const Vec3d line_vector = cross(normal_a, normal_b);

  FacePairCut cut;
  if (dot(line_vector, line_vector) <=
      detail::kFaceIntersectParallelEpsilon * detail::kFaceIntersectParallelEpsilon) {
    // Parallel planes: coplanar when their distances from the origin agree
    // after orienting B's normal like A's. Antiparallel normals describe the
    // same plane, so B's distance is negated before the comparison.
    const Plane3 unit_a = plane_frame_plane(frame_a);
    const Plane3 unit_b = plane_frame_plane(frame_b);
    const double oriented_b =
        dot(normal_a, normal_b) < 0.0 ? -unit_b.offset : unit_b.offset;
    if (std::fabs(unit_a.offset - oriented_b) > tol) return cut;
    cut.coplanar = true;

    const Result<std::vector<std::vector<Point2d>>> loops_a =
        face_loops_2d(a, fa, frame_a);
    if (!loops_a) return loops_a.error();
    const Result<std::vector<std::vector<Point2d>>> loops_b =
        face_loops_2d(b, fb, frame_a);
    if (!loops_b) return loops_b.error();
    const std::array<Polygon2d, 2> operands = {
        detail::face_intersect_polygon(loops_a.value()),
        detail::face_intersect_polygon(loops_b.value())};
    const Result<BooleanResult> intersection =
        boolean2(operands, BooleanOp::intersection);
    if (!intersection) return intersection.error();
    detail::face_intersect_emit_polygons(intersection.value().polygons, frame_a,
                                         cut.segments);
    return cut;
  }

  // The crossing line solves nA . x = bA, nB . x = bB and d . x = 0, where d is
  // the cross product of the normals. Cramer's rule on the 3x3 matrix whose
  // columns are the normal components and d.
  const Plane3 unit_a = plane_frame_plane(frame_a);
  const Plane3 unit_b = plane_frame_plane(frame_b);
  const Vec3d column0{normal_a[0], normal_b[0], line_vector[0]};
  const Vec3d column1{normal_a[1], normal_b[1], line_vector[1]};
  const Vec3d column2{normal_a[2], normal_b[2], line_vector[2]};
  const Vec3d rhs{unit_a.offset, unit_b.offset, 0.0};
  const double det = dot(column0, cross(column1, column2));
  if (!(det != 0.0)) {
    return Error{Errc::numeric_failure, "face_intersect: crossing line is degenerate"};
  }
  const Point3d p0{dot(rhs, cross(column1, column2)) / det,
                  dot(column0, cross(rhs, column2)) / det,
                  dot(column0, cross(column1, rhs)) / det};

  // Clipping parameters and the emitted points both use the unit direction, so
  // the interval parameter is a distance along the line and cannot be scaled by
  // the cross product magnitude.
  const double line_length = std::sqrt(dot(line_vector, line_vector));
  const Vec3d line_direction{line_vector[0] / line_length,
                             line_vector[1] / line_length,
                             line_vector[2] / line_length};

  const Result<std::vector<detail::FaceIntersectInterval>> intervals_a =
      detail::face_intersect_clip_line(a, fa, frame_a, p0, line_direction, tol);
  if (!intervals_a) return intervals_a.error();
  const Result<std::vector<detail::FaceIntersectInterval>> intervals_b =
      detail::face_intersect_clip_line(b, fb, frame_b, p0, line_direction, tol);
  if (!intervals_b) return intervals_b.error();

  const std::vector<detail::FaceIntersectInterval> merged_a =
      detail::face_intersect_merge_intervals(intervals_a.value());
  const std::vector<detail::FaceIntersectInterval> merged_b =
      detail::face_intersect_merge_intervals(intervals_b.value());

  std::size_t i = 0;
  std::size_t j = 0;
  while (i < merged_a.size() && j < merged_b.size()) {
    const double lo = std::max(merged_a[i].lo, merged_b[j].lo);
    const double hi = std::min(merged_a[i].hi, merged_b[j].hi);
    if (hi > lo) {
      cut.segments.push_back(detail::face_intersect_ordered(
          p0 + line_direction * lo, p0 + line_direction * hi));
    }
    if (merged_a[i].hi < merged_b[j].hi) {
      ++i;
    } else {
      ++j;
    }
  }
  return cut;
}

}  // namespace v1
}  // namespace geom
