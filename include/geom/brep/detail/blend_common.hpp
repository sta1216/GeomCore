#pragma once

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <span>
#include <string_view>
#include <vector>

#include <geom/brep/topology.hpp>
#include <geom/core/error.hpp>
#include <geom/core/result.hpp>
#include <geom/curves/analytic_curve.hpp>
#include <geom/primitives/point.hpp>
#include <geom/primitives/vec.hpp>
#include <geom/surfaces/plane.hpp>

namespace geom {
inline namespace v1 {
namespace detail {

// Shared planar-convex view of a solid, used by every blend operation that
// replaces one edge or face by a planar cut: a chamfer (M4c-1) and the shell
// wall offset both reduce to "cut this convex planar solid with a plane".
//
// The view is a snapshot of geometry only: it borrows nothing from the Brep,
// so callers may keep cutting copies of the input while this view stays valid.
struct BlendPlanarSolid {
  SolidId solid;                    // the single solid of the input
  std::vector<FaceId> faces;        // shell faces, in shell order
  std::vector<Vec3d> outward;       // unit outward normal per face
  std::vector<double> offset;       // plane offset per face: dot(outward, x) = offset
  double scale = 1.0;               // max(1, max |vertex coordinate|)
};

// Error::message is a static `const char*`, not a std::string: a caller tag
// cannot be interpolated at run time. Each diagnostic therefore keeps one
// complete literal per public entry point, and the tag selects which one.
// This is the reason the tag is threaded through every helper below.
inline Error blend_error(std::string_view caller, Errc code, const char* chamfer_message,
                         const char* shell_message) noexcept {
  return Error{code, caller == "shell_solid" ? shell_message : chamfer_message};
}

// Point3d and Vec3d are distinct types with no implicit conversion; geometry
// code needs the vector form often enough to warrant one named helper.
inline Vec3d blend_as_vec(Point3d p) noexcept {
  return Vec3d{p[0], p[1], p[2]};
}

// max(1, max |coordinate|) over the vertices: the reference length that turns
// the absolute tolerance into a relative one and keeps 1e-12-style floors
// meaningful for models far from the unit cube.
inline double blend_scale_of(const Brep& brep) noexcept {
  double scale = 1.0;
  for (const VertexId id : brep.vertex_ids()) {
    const Point3d p = brep.vertex(id).point;
    scale = std::max({scale, std::fabs(p[0]), std::fabs(p[1]), std::fabs(p[2])});
  }
  return scale;
}

// Unit outward normal and normalized offset of one planar face, with
// Face::reversed applied (the flag flips both the normal and the offset, since
// the stored plane equation is normal . x = offset).
struct BlendFacePlane {
  Vec3d outward;
  double offset = 0.0;
};

inline Result<BlendFacePlane> blend_face_plane(const Brep& brep, FaceId face_id,
                                               std::string_view caller) {
  const Face& face = brep.face(face_id);
  const GeometryStore& geometry = brep.geometry();
  if (face.surface == kNoGeometry ||
      static_cast<std::size_t>(face.surface) >= geometry.surfaces.size()) {
    return blend_error(caller, Errc::unsupported,
                       "chamfer_edges: all faces must be planar",
                       "shell_solid: all faces must be planar");
  }
  const Plane3* plane = std::get_if<Plane3>(&geometry.surfaces[face.surface]);
  if (plane == nullptr) {
    return blend_error(caller, Errc::unsupported,
                       "chamfer_edges: all faces must be planar",
                       "shell_solid: all faces must be planar");
  }
  const Vec3d n = plane->normal;
  const double len = std::hypot(n[0], n[1], n[2]);
  if (!(len > 0.0) || !std::isfinite(len) || !std::isfinite(plane->offset)) {
    return blend_error(caller, Errc::degenerate,
                       "chamfer_edges: face plane is degenerate",
                       "shell_solid: face plane is degenerate");
  }
  Vec3d outward{n[0] / len, n[1] / len, n[2] / len};
  double offset = plane->offset / len;
  if (face.reversed) {
    outward = outward * -1.0;
    offset = -offset;
  }
  return BlendFacePlane{outward, offset};
}

// Boundary vertices of a face in loop order: for each loop the traversal starts
// at Loop::start and follows Coedge::next, taking the edge endpoint the coedge
// starts from. Consecutive points (last wraps to first) are real boundary
// segments, which blend_face_width relies on.
inline Result<std::vector<Point3d>> blend_face_loop_points(const Brep& brep, FaceId face_id,
                                                           std::string_view caller) {
  const Face& face = brep.face(face_id);
  if (face.first_loop.is_null()) {
    return blend_error(caller, Errc::invalid_topology,
                       "chamfer_edges: face has no boundary loop",
                       "shell_solid: face has no boundary loop");
  }
  const std::size_t coedge_limit = brep.coedge_ids().size();
  const std::size_t loop_limit = brep.loop_ids().size();
  std::vector<Point3d> points;
  LoopId loop_id = face.first_loop;
  std::size_t loop_guard = 0;
  while (!loop_id.is_null()) {
    const Loop& loop = brep.loop(loop_id);
    if (loop.start.is_null()) {
      return blend_error(caller, Errc::invalid_topology,
                         "chamfer_edges: face loop has no start coedge",
                         "shell_solid: face loop has no start coedge");
    }
    CoedgeId coedge_id = loop.start;
    std::size_t coedge_guard = 0;
    do {
      const Coedge& coedge = brep.coedge(coedge_id);
      const Edge& edge = brep.edge(coedge.edge);
      points.push_back(brep.vertex(coedge.reversed ? edge.v1 : edge.v0).point);
      coedge_id = coedge.next;
      if (coedge_id.is_null() || ++coedge_guard > coedge_limit) {
        return blend_error(caller, Errc::invalid_topology,
                           "chamfer_edges: face loop is not closed",
                           "shell_solid: face loop is not closed");
      }
    } while (coedge_id != loop.start);
    loop_id = loop.next;
    if (++loop_guard > loop_limit) {
      return blend_error(caller, Errc::invalid_topology,
                         "chamfer_edges: face loop chain is cyclic",
                         "shell_solid: face loop chain is cyclic");
    }
  }
  return points;
}

// Average of the boundary vertices. Every corner of a simple loop appears
// twice, so the average is the vertex centroid of the face; the solid has
// already been proven convex, so that centroid is strictly inside the face.
inline Result<Point3d> blend_face_centroid(const Brep& brep, FaceId face_id,
                                           std::string_view caller) {
  const Result<std::vector<Point3d>> loop =
      blend_face_loop_points(brep, face_id, caller);
  if (!loop.has_value()) return loop.error();
  if (loop.value().empty()) {
    return blend_error(caller, Errc::invalid_topology,
                       "chamfer_edges: face has no boundary",
                       "shell_solid: face has no boundary");
  }
  double sum[3] = {0.0, 0.0, 0.0};
  for (const Point3d p : loop.value()) {
    sum[0] += p[0];
    sum[1] += p[1];
    sum[2] += p[2];
  }
  const double count = static_cast<double>(loop.value().size());
  return Point3d{sum[0] / count, sum[1] / count, sum[2] / count};
}

// Validates and describes the one solid/one closed planar shell this milestone
// supports, and proves convexity: every vertex must lie on the inner side of
// every face plane within max(64 * tolerance, 1e-12 * scale). The 64x factor
// matches split_solid's weld slack so a model that survived assembly is not
// rejected by a tighter convexity probe.
inline Result<BlendPlanarSolid> blend_planar_convex_solid(const Brep& brep, double tolerance,
                                                          std::string_view caller) {
  if (!std::isfinite(tolerance) || !(tolerance > 0.0)) {
    return blend_error(caller, Errc::tolerance_violation,
                       "chamfer_edges: tolerance must be finite and positive",
                       "shell_solid: tolerance must be finite and positive");
  }
  if (brep.solid_count() != 1u || brep.shell_count() != 1u) {
    return blend_error(caller, Errc::unsupported,
                       "chamfer_edges: exactly one solid with one shell is required",
                       "shell_solid: exactly one solid with one shell is required");
  }
  const SolidId solid = brep.solid_ids().front();
  const std::span<const ShellId> solid_shells = brep.solid_shells(solid);
  if (solid_shells.size() != 1u) {
    return blend_error(caller, Errc::unsupported,
                       "chamfer_edges: exactly one solid with one shell is required",
                       "shell_solid: exactly one solid with one shell is required");
  }
  const ShellId shell = solid_shells.front();
  if (!brep.shell(shell).closed) {
    return blend_error(caller, Errc::invalid_topology,
                       "chamfer_edges: solid shell must be closed",
                       "shell_solid: solid shell must be closed");
  }

  BlendPlanarSolid view;
  view.solid = solid;
  view.scale = blend_scale_of(brep);
  const std::span<const FaceId> shell_faces = brep.shell_faces(shell);
  view.faces.assign(shell_faces.begin(), shell_faces.end());
  view.outward.reserve(view.faces.size());
  view.offset.reserve(view.faces.size());
  for (const FaceId face_id : view.faces) {
    const Result<BlendFacePlane> plane = blend_face_plane(brep, face_id, caller);
    if (!plane.has_value()) return plane.error();
    view.outward.push_back(plane.value().outward);
    view.offset.push_back(plane.value().offset);
  }

  const double convex_tolerance = std::max(64.0 * tolerance, 1e-12 * view.scale);
  for (const VertexId vertex_id : brep.vertex_ids()) {
    const Vec3d point = blend_as_vec(brep.vertex(vertex_id).point);
    for (std::size_t i = 0; i < view.outward.size(); ++i) {
      if (dot(view.outward[i], point) - view.offset[i] > convex_tolerance) {
        return blend_error(caller, Errc::unsupported,
                           "chamfer_edges: solid must be convex",
                           "shell_solid: solid must be convex");
      }
    }
  }
  return view;
}

// Unit direction of a straight edge. Only LineCurve3 edges can be chamfered:
// a curved edge has no single in-face direction, and silently flattening it
// would produce a cut that misses the edge it claims to replace.
inline Result<Vec3d> blend_edge_direction(const Brep& brep, EdgeId edge_id,
                                          std::string_view caller) {
  if (!brep.edge_alive(edge_id)) {
    return blend_error(caller, Errc::index_out_of_range,
                       "chamfer_edges: edge does not belong to the solid",
                       "shell_solid: edge does not belong to the solid");
  }
  const Edge& edge = brep.edge(edge_id);
  const GeometryStore& geometry = brep.geometry();
  if (edge.curve == kNoGeometry ||
      static_cast<std::size_t>(edge.curve) >= geometry.curves.size()) {
    return blend_error(caller, Errc::unsupported,
                       "chamfer_edges: only straight edges can be chamfered",
                       "shell_solid: only straight edges can be chamfered");
  }
  const LineCurve3* line = std::get_if<LineCurve3>(&geometry.curves[edge.curve]);
  if (line == nullptr) {
    return blend_error(caller, Errc::unsupported,
                       "chamfer_edges: only straight edges can be chamfered",
                       "shell_solid: only straight edges can be chamfered");
  }
  const Vec3d direction = line->direction;
  const double len = std::hypot(direction[0], direction[1], direction[2]);
  if (!(len > 0.0) || !std::isfinite(len)) {
    return blend_error(caller, Errc::degenerate,
                       "chamfer_edges: edge direction is degenerate",
                       "shell_solid: edge direction is degenerate");
  }
  return Vec3d{direction[0] / len, direction[1] / len, direction[2] / len};
}

// In-face direction perpendicular to the edge and pointing into the face:
// cross(outward, edge_dir) gives one of the two perpendiculars, and the face
// centroid picks the one on the interior side. The edge must actually border
// the face, otherwise the cross product has no geometric meaning.
inline Result<Vec3d> blend_inward_direction(const Brep& brep, FaceId face_id, EdgeId edge_id,
                                            const Vec3d& edge_dir, const Vec3d& edge_mid,
                                            std::string_view caller) {
  bool adjacent = false;
  for (const CoedgeId id : brep.coedge_ids()) {
    const Coedge& coedge = brep.coedge(id);
    if (coedge.edge == edge_id && coedge.face == face_id) {
      adjacent = true;
      break;
    }
  }
  if (!adjacent) {
    return blend_error(caller, Errc::invalid_topology,
                       "chamfer_edges: edge is not adjacent to the face",
                       "shell_solid: edge is not adjacent to the face");
  }
  const Result<BlendFacePlane> plane = blend_face_plane(brep, face_id, caller);
  if (!plane.has_value()) return plane.error();
  const Result<Point3d> centroid = blend_face_centroid(brep, face_id, caller);
  if (!centroid.has_value()) return centroid.error();

  Vec3d inward = cross(plane.value().outward, edge_dir);
  const double len = std::hypot(inward[0], inward[1], inward[2]);
  if (!(len > 0.0) || !std::isfinite(len)) {
    return blend_error(caller, Errc::degenerate,
                       "chamfer_edges: edge is parallel to the face normal",
                       "shell_solid: edge is parallel to the face normal");
  }
  inward = inward * (1.0 / len);
  const Vec3d to_centroid = blend_as_vec(centroid.value()) - edge_mid;
  if (dot(inward, to_centroid) < 0.0) inward = inward * -1.0;
  return inward;
}

// Distance from the edge (given by its midpoint and the in-face direction)
// across the face to the opposite boundary. The boundary is projected into the
// frame (direction, cross(outward, direction)); the +direction ray through the
// origin is intersected with every boundary segment, and the smallest positive
// hit wins. The edge's own endpoints project onto x = 0 and are skipped by the
// positive-distance test, so the opposite side is what gets measured.
inline Result<double> blend_face_width(const Brep& brep, FaceId face_id, const Vec3d& origin,
                                       const Vec3d& direction, std::string_view caller) {
  const Result<BlendFacePlane> plane = blend_face_plane(brep, face_id, caller);
  if (!plane.has_value()) return plane.error();
  const double direction_len = std::hypot(direction[0], direction[1], direction[2]);
  if (!(direction_len > 0.0) || !std::isfinite(direction_len)) {
    return blend_error(caller, Errc::degenerate,
                       "chamfer_edges: width direction is degenerate",
                       "shell_solid: width direction is degenerate");
  }
  const Vec3d u = direction * (1.0 / direction_len);
  Vec3d v = cross(plane.value().outward, u);
  const double v_len = std::hypot(v[0], v[1], v[2]);
  if (!(v_len > 0.0) || !std::isfinite(v_len)) {
    return blend_error(caller, Errc::degenerate,
                       "chamfer_edges: width direction lies in the face plane",
                       "shell_solid: width direction lies in the face plane");
  }
  v = v * (1.0 / v_len);

  const Result<std::vector<Point3d>> loop = blend_face_loop_points(brep, face_id, caller);
  if (!loop.has_value()) return loop.error();
  if (loop.value().size() < 3u) {
    return blend_error(caller, Errc::invalid_topology,
                       "chamfer_edges: face boundary is degenerate",
                       "shell_solid: face boundary is degenerate");
  }
  std::vector<Vec2d> projected;
  projected.reserve(loop.value().size());
  for (const Point3d p : loop.value()) {
    const Vec3d d = blend_as_vec(p) - origin;
    projected.push_back(Vec2d{dot(d, u), dot(d, v)});
  }

  // Anything closer than this to the origin is the selected edge itself, not
  // the opposite side of the face.
  const double eps = 1e-12 * blend_scale_of(brep);
  double best = std::numeric_limits<double>::infinity();
  for (std::size_t i = 0; i < projected.size(); ++i) {
    const Vec2d a = projected[i];
    const Vec2d b = projected[(i + 1) % projected.size()];
    if (a[1] == b[1]) continue;  // segment parallel to the ray
    const double t = a[1] / (a[1] - b[1]);
    if (t < -1e-12 || t > 1.0 + 1e-12) continue;
    const double x = a[0] + t * (b[0] - a[0]);
    if (x > eps && x < best) best = x;
  }
  if (!std::isfinite(best)) {
    return blend_error(caller, Errc::degenerate,
                       "chamfer_edges: face has no width in the given direction",
                       "shell_solid: face has no width in the given direction");
  }
  return best;
}

}  // namespace detail
}  // namespace v1
}  // namespace geom
