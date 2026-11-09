#pragma once

// Exact point-in-face and point-in-shell queries on a B-rep.
//
// Containment is decided by exact 2D orientation predicates. The tolerance is
// used only to accept a point as lying on a face plane and to measure the
// distance to a boundary segment; it never turns a near miss into a hit.
// Surfaces that are not planes are reported as unsupported.

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <span>
#include <utility>
#include <variant>
#include <vector>

#include <geom/brep/plane_frame.hpp>
#include <geom/brep/topology.hpp>
#include <geom/core/error.hpp>
#include <geom/core/result.hpp>
#include <geom/numeric/predicates2.hpp>
#include <geom/primitives/point.hpp>
#include <geom/surfaces/plane.hpp>

namespace geom {
inline namespace v1 {

namespace detail {

// Exact point-in-ring test. A half-open ray crossing is counted from each
// directed edge; points exactly on an edge may be classified either way,
// which is why boundary membership is a separate query. A ring with fewer
// than three points has no interior.
inline bool face_query_point_in_ring(const std::vector<Point2d>& ring,
                                     const Point2d& p) noexcept {
  const std::size_t count = ring.size();
  if (count < 3) return false;
  int winding = 0;
  for (std::size_t i = 0; i < count; ++i) {
    const Point2d& a = ring[i];
    const Point2d& b = ring[(i + 1) % count];
    const int side = sign_orient2d(a[0], a[1], b[0], b[1], p[0], p[1]);
    if (a[1] <= p[1]) {
      if (b[1] > p[1] && side > 0) ++winding;
    } else if (b[1] <= p[1] && side < 0) {
      --winding;
    }
  }
  return winding != 0;
}

// Distance from p to the segment a-b. A zero-length segment degenerates to
// the distance to its single point.
inline double face_query_point_segment_distance(const Point2d& p, const Point2d& a,
                                                const Point2d& b) noexcept {
  const Vec2d ab = b - a;
  const Vec2d ap = p - a;
  const double length2 = dot(ab, ab);
  if (!(length2 > 0.0)) return norm2(ap);
  const double t = std::clamp(dot(ap, ab) / length2, 0.0, 1.0);
  return norm2(ap - ab * t);
}

inline Result<bool> face_query_check_face(const Brep& brep, FaceId id) {
  if (id.is_null() || !brep.face_alive(id)) {
    return Error{Errc::index_out_of_range, "face_query: face is not alive", id.index(), 0u};
  }
  return true;
}

inline Result<bool> face_query_check_tolerance(double tol) {
  if (!std::isfinite(tol)) {
    return Error{Errc::tolerance_violation, "face_query: tolerance is not finite"};
  }
  if (tol < 0.0) {
    return Error{Errc::tolerance_violation, "face_query: tolerance is negative"};
  }
  return true;
}

inline Result<bool> face_query_check_point(const Point3d& p) {
  if (!std::isfinite(p[0]) || !std::isfinite(p[1]) || !std::isfinite(p[2])) {
    return Error{Errc::non_finite, "face_query: point has a non-finite component"};
  }
  return true;
}

}  // namespace detail

// The face's surface as a plane. A face without a surface, a surface index
// outside the geometry store and any curved surface are unsupported.
inline Result<Plane3> face_plane(const Brep& brep, FaceId id) {
  const Result<bool> face_check = detail::face_query_check_face(brep, id);
  if (!face_check) return face_check.error();
  const Face& face = brep.face(id);
  if (face.surface == kNoGeometry) {
    return Error{Errc::unsupported, "face_query: face has no surface", id.index(), 0u};
  }
  const GeometryStore& store = brep.geometry();
  if (static_cast<std::size_t>(face.surface) >= store.surfaces.size()) {
    return Error{Errc::unsupported, "face_query: face surface is missing", id.index(), 0u};
  }
  const auto* plane =
      std::get_if<Plane3>(&store.surfaces[static_cast<std::size_t>(face.surface)]);
  if (plane == nullptr) {
    return Error{Errc::unsupported, "face_query: face surface is not planar", id.index(), 0u};
  }
  return *plane;
}

// The face boundary projected onto `frame` as 2D rings: the outer loop first,
// then the holes in chain order. Adjacent coedges share their endpoints, so
// only each coedge's start vertex is emitted; a repeated closing point is
// dropped. The frame is usually make_plane_frame(face_plane(brep, id)).
inline Result<std::vector<std::vector<Point2d>>> face_loops_2d(const Brep& brep, FaceId id,
                                                              const PlaneFrame& frame) {
  const Result<bool> face_check = detail::face_query_check_face(brep, id);
  if (!face_check) return face_check.error();
  // A curved face has no single projection plane, so reject it before
  // walking any loop.
  const Result<Plane3> plane = face_plane(brep, id);
  if (!plane) return plane.error();

  struct LoopRing {
    bool outer = false;
    std::vector<Point2d> points;
  };
  std::vector<LoopRing> rings;

  const std::size_t loop_limit = brep.loops().size();
  const std::size_t coedge_limit = brep.coedges().size();
  std::size_t loop_count = 0;
  LoopId loop_id = brep.face(id).first_loop;
  while (!loop_id.is_null()) {
    if (loop_count++ >= loop_limit) {
      return Error{Errc::invalid_topology, "face_query: loop chain is cyclic", id.index(), 0u};
    }
    if (!brep.loop_alive(loop_id)) {
      return Error{Errc::invalid_topology, "face_query: loop is not alive", loop_id.index(), 0u};
    }
    const Loop& loop = brep.loop(loop_id);
    LoopRing ring;
    ring.outer = loop.outer;
    if (loop.start.is_null()) {
      return Error{Errc::invalid_topology, "face_query: loop has no coedges", loop_id.index(), 0u};
    }
    std::size_t coedge_count = 0;
    CoedgeId coedge_id = loop.start;
    do {
      if (coedge_count++ >= coedge_limit) {
        return Error{Errc::invalid_topology, "face_query: coedge cycle is cyclic", loop_id.index(),
                     0u};
      }
      if (!brep.coedge_alive(coedge_id)) {
        return Error{Errc::invalid_topology, "face_query: coedge is not alive", coedge_id.index(),
                     0u};
      }
      const Coedge& coedge = brep.coedge(coedge_id);
      if (!brep.edge_alive(coedge.edge)) {
        return Error{Errc::invalid_topology, "face_query: edge is not alive", coedge.edge.index(),
                     0u};
      }
      const Edge& edge = brep.edge(coedge.edge);
      const VertexId vertex_id = coedge.reversed ? edge.v1 : edge.v0;
      if (!brep.vertex_alive(vertex_id)) {
        return Error{Errc::invalid_topology, "face_query: vertex is not alive", vertex_id.index(),
                     0u};
      }
      ring.points.push_back(plane_frame_to_2d(frame, brep.vertex(vertex_id).point));
      coedge_id = coedge.next;
    } while (!coedge_id.is_null() && coedge_id != loop.start);
    if (coedge_id.is_null()) {
      return Error{Errc::invalid_topology, "face_query: coedge cycle is not closed",
                   loop_id.index(), 0u};
    }
    if (ring.points.size() > 1 && ring.points.front()[0] == ring.points.back()[0] &&
        ring.points.front()[1] == ring.points.back()[1]) {
      ring.points.pop_back();
    }
    rings.push_back(std::move(ring));
    loop_id = loop.next;
  }

  std::vector<std::vector<Point2d>> loops;
  loops.reserve(rings.size());
  for (LoopRing& ring : rings) {
    if (ring.outer) loops.push_back(std::move(ring.points));
  }
  for (LoopRing& ring : rings) {
    if (!ring.outer) loops.push_back(std::move(ring.points));
  }
  return loops;
}

// True when p lies on the face: on its plane within `tol`, inside the outer
// loop and outside every hole.
inline Result<bool> point_in_face(const Brep& brep, FaceId id, const Point3d& p, double tol) {
  const Result<bool> face_check = detail::face_query_check_face(brep, id);
  if (!face_check) return face_check.error();
  const Result<bool> tol_check = detail::face_query_check_tolerance(tol);
  if (!tol_check) return tol_check.error();
  const Result<bool> point_check = detail::face_query_check_point(p);
  if (!point_check) return point_check.error();

  const Result<Plane3> plane = face_plane(brep, id);
  if (!plane) return plane.error();
  const Result<PlaneFrame> frame_result = make_plane_frame(plane.value());
  if (!frame_result) return frame_result.error();
  const PlaneFrame& frame = frame_result.value();
  if (std::fabs(dot(p - frame.origin, frame.normal)) > tol) return false;

  const Result<std::vector<std::vector<Point2d>>> loops = face_loops_2d(brep, id, frame);
  if (!loops) return loops.error();
  const std::vector<std::vector<Point2d>>& rings = loops.value();
  if (rings.empty()) return false;
  const Point2d q = plane_frame_to_2d(frame, p);
  if (!detail::face_query_point_in_ring(rings.front(), q)) return false;
  for (std::size_t i = 1; i < rings.size(); ++i) {
    if (detail::face_query_point_in_ring(rings[i], q)) return false;
  }
  return true;
}

// True when p is on the face plane within `tol` and within `tol` of any
// boundary segment. This is the only query where a tolerance decides a hit.
inline Result<bool> point_on_face_boundary(const Brep& brep, FaceId id, const Point3d& p,
                                           double tol) {
  const Result<bool> face_check = detail::face_query_check_face(brep, id);
  if (!face_check) return face_check.error();
  const Result<bool> tol_check = detail::face_query_check_tolerance(tol);
  if (!tol_check) return tol_check.error();
  const Result<bool> point_check = detail::face_query_check_point(p);
  if (!point_check) return point_check.error();

  const Result<Plane3> plane = face_plane(brep, id);
  if (!plane) return plane.error();
  const Result<PlaneFrame> frame_result = make_plane_frame(plane.value());
  if (!frame_result) return frame_result.error();
  const PlaneFrame& frame = frame_result.value();
  if (std::fabs(dot(p - frame.origin, frame.normal)) > tol) return false;

  const Result<std::vector<std::vector<Point2d>>> loops = face_loops_2d(brep, id, frame);
  if (!loops) return loops.error();
  const Point2d q = plane_frame_to_2d(frame, p);
  for (const std::vector<Point2d>& ring : loops.value()) {
    if (ring.size() < 2) continue;
    for (std::size_t i = 0; i < ring.size(); ++i) {
      const Point2d& a = ring[i];
      const Point2d& b = ring[(i + 1) % ring.size()];
      if (detail::face_query_point_segment_distance(q, a, b) <= tol) return true;
    }
  }
  return false;
}

// Ray-parity containment for a closed shell. Each direction is tried in turn;
// a hit that lands on a face boundary is ambiguous, so that direction is
// discarded. When every direction is ambiguous the query fails instead of
// guessing.
inline Result<bool> point_in_shell(const Brep& brep, ShellId id, const Point3d& p, double tol) {
  if (id.is_null() || !brep.shell_alive(id)) {
    return Error{Errc::index_out_of_range, "face_query: shell is not alive", id.index(), 0u};
  }
  const Result<bool> tol_check = detail::face_query_check_tolerance(tol);
  if (!tol_check) return tol_check.error();
  const Result<bool> point_check = detail::face_query_check_point(p);
  if (!point_check) return point_check.error();

  const std::span<const FaceId> faces = brep.shell_faces(id);
  if (faces.empty()) {
    return Error{Errc::index_out_of_range, "face_query: shell has no faces", id.index(), 0u};
  }

  const Vec3d point{p[0], p[1], p[2]};

  // Four regular-tetrahedron directions. None is axis aligned and any three
  // are independent, so several directions stay usable on a box.
  constexpr std::array<Vec3d, 4> kRawDirections = {{
      Vec3d{0.5773, 0.5774, 0.5775},
      Vec3d{0.5774, -0.5775, -0.5773},
      Vec3d{-0.5775, 0.5773, -0.5774},
      Vec3d{-0.5773, -0.5774, 0.5775},
  }};

  for (const Vec3d& raw : kRawDirections) {
    const double length = std::sqrt(dot(raw, raw));
    const Vec3d direction{raw[0] / length, raw[1] / length, raw[2] / length};
    int parity = 0;
    bool ambiguous = false;
    for (const FaceId face_id : faces) {
      const Result<Plane3> plane = face_plane(brep, face_id);
      if (!plane) return plane.error();
      const Result<PlaneFrame> frame_result = make_plane_frame(plane.value());
      if (!frame_result) return frame_result.error();
      const Plane3 unit_plane = plane_frame_plane(frame_result.value());
      const double den = dot(unit_plane.normal, direction);
      if (!(std::fabs(den) > 1e-12)) continue;
      const double t = (unit_plane.offset - dot(unit_plane.normal, point)) / den;
      if (!(t > tol)) continue;
      const Point3d hit = p + direction * t;
      const Result<bool> boundary = point_on_face_boundary(brep, face_id, hit, tol);
      if (!boundary) return boundary.error();
      if (boundary.value()) {
        ambiguous = true;
        break;
      }
      const Result<bool> inside = point_in_face(brep, face_id, hit, tol);
      if (!inside) return inside.error();
      if (inside.value()) ++parity;
    }
    if (ambiguous) continue;
    return parity % 2 == 1;
  }
  return Error{Errc::numeric_failure, "face_query: no reliable ray direction", id.index(), 0u};
}

}  // namespace v1
}  // namespace geom
