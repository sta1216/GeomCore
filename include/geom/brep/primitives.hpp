#pragma once

#include <array>
#include <cmath>
#include <cstdint>
#include <numbers>

#include <geom/brep/solid_editor.hpp>
#include <geom/brep/topology.hpp>
#include <geom/core/error.hpp>
#include <geom/core/result.hpp>
#include <geom/curves/analytic_curve.hpp>
#include <geom/primitives/point.hpp>
#include <geom/primitives/vec.hpp>
#include <geom/surfaces/analytic.hpp>
#include <geom/surfaces/plane.hpp>

namespace geom {
inline namespace v1 {

namespace detail {

inline bool primitive_point_finite(Point3d p) noexcept {
  return std::isfinite(p[0]) && std::isfinite(p[1]) && std::isfinite(p[2]);
}

inline bool primitive_vec_finite(Vec3d v) noexcept {
  return std::isfinite(v[0]) && std::isfinite(v[1]) && std::isfinite(v[2]);
}

inline bool primitive_point_equal(Point3d a, Point3d b) noexcept {
  return a[0] == b[0] && a[1] == b[1] && a[2] == b[2];
}

}  // namespace detail

// Axis-aligned box [lo, hi] as a closed B-rep solid.
//
// The eight corners are placed in the order
//
//   0 (lo.x, lo.y, lo.z)   1 (hi.x, lo.y, lo.z)   2 (hi.x, hi.y, lo.z)   3 (lo.x, hi.y, lo.z)
//   4 (lo.x, lo.y, hi.z)   5 (hi.x, lo.y, hi.z)   6 (hi.x, hi.y, hi.z)   7 (lo.x, hi.y, hi.z)
//
// and the twelve edges form the bottom ring (0-1, 1-2, 2-3, 3-0), the top ring
// (4-5, 5-6, 6-7, 7-4) and the four verticals (1-5, 2-6, 3-7, 0-4). Each edge is
// a unit-speed LineCurve3 from v0 to v1 over [0, |v1 - v0|].
//
// The six planar faces are ordered bottom (-z), top (+z), front (-y), back (+y),
// left (-x), right (+x). Each face's single outer loop runs counter-clockwise
// when viewed from outside along the face normal, so the right-hand rule of the
// loop reproduces the stored outward normal; the plane is
// Plane3{normal, dot(normal, corner)} with a unit axis normal. The shell is
// closed and the solid owns it.
//
// Validation:
//   - non-finite lo or hi                -> Errc::non_finite
//   - non-finite or negative tolerance   -> Errc::tolerance_violation
//   - lo[axis] >= hi[axis] on any axis   -> Errc::degenerate
// The three checks above run before any editor state is allocated. Edge
// construction then checks every length; a length that overflows to a
// non-finite value -> Errc::degenerate, and that edge (with its curve) is not
// staged. Any error from a SolidEditor add_* call is propagated unchanged;
// commit() performs the final validation.
inline Result<Brep> make_box(Point3d lo, Point3d hi, double tol) {
  if (!detail::primitive_point_finite(lo) || !detail::primitive_point_finite(hi)) {
    return Error{Errc::non_finite, "make_box: box corner is not finite", 0u, 0u};
  }
  if (!std::isfinite(tol) || tol < 0.0) {
    return Error{Errc::tolerance_violation,
                 "make_box: tolerance must be finite and non-negative", 0u, 0u};
  }
  for (int axis = 0; axis < 3; ++axis) {
    if (!(lo[axis] < hi[axis])) {
      return Error{Errc::degenerate, "make_box: box is degenerate or inverted", 0u, 0u};
    }
  }

  const double x0 = lo[0];
  const double y0 = lo[1];
  const double z0 = lo[2];
  const double x1 = hi[0];
  const double y1 = hi[1];
  const double z1 = hi[2];
  const Point3d corners[8] = {
      {x0, y0, z0}, {x1, y0, z0}, {x1, y1, z0}, {x0, y1, z0},
      {x0, y0, z1}, {x1, y0, z1}, {x1, y1, z1}, {x0, y1, z1},
  };

  SolidEditor editor;

  VertexId vertices[8];
  for (int i = 0; i < 8; ++i) {
    const Result<VertexId> vertex = editor.add_vertex(corners[i], tol);
    if (!vertex.has_value()) return vertex.error();
    vertices[i] = vertex.value();
  }

  struct EdgeSpec {
    int v0;
    int v1;
  };
  const EdgeSpec edge_specs[12] = {
      {0, 1}, {1, 2}, {2, 3}, {3, 0},  // bottom ring
      {4, 5}, {5, 6}, {6, 7}, {7, 4},  // top ring
      {1, 5}, {2, 6}, {3, 7}, {0, 4},  // verticals
  };

  EdgeId edges[12];
  for (int i = 0; i < 12; ++i) {
    const Point3d from = corners[edge_specs[i].v0];
    const Point3d to = corners[edge_specs[i].v1];
    const Vec3d delta = to - from;
    const double length = std::hypot(delta[0], delta[1], delta[2]);
    if (!std::isfinite(length) || !(length > 0.0)) {
      return Error{Errc::degenerate, "make_box: box edge length is not representable", 0u,
                   0u};
    }
    const Vec3d direction{delta[0] / length, delta[1] / length, delta[2] / length};
    const std::uint32_t curve = editor.add_curve(LineCurve3{from, direction});
    const Result<EdgeId> edge = editor.add_edge(
        vertices[edge_specs[i].v0], vertices[edge_specs[i].v1], curve, 0.0, length, tol);
    if (!edge.has_value()) return edge.error();
    edges[i] = edge.value();
  }

  struct CoedgeUse {
    int edge;
    bool reversed;
  };
  struct FaceSpec {
    Vec3d normal;
    int corner;
    CoedgeUse loop[4];
  };
  // Loop order is counter-clockwise seen from outside along `normal`.
  const FaceSpec face_specs[6] = {
      {{0.0, 0.0, -1.0}, 0, {{3, true}, {2, true}, {1, true}, {0, true}}},    // bottom
      {{0.0, 0.0, 1.0}, 4, {{4, false}, {5, false}, {6, false}, {7, false}}},  // top
      {{0.0, -1.0, 0.0}, 0, {{0, false}, {8, false}, {4, true}, {11, true}}},  // front
      {{0.0, 1.0, 0.0}, 2, {{2, false}, {10, false}, {6, true}, {9, true}}},   // back
      {{-1.0, 0.0, 0.0}, 0, {{11, false}, {7, true}, {10, true}, {3, false}}}, // left
      {{1.0, 0.0, 0.0}, 1, {{1, false}, {9, false}, {5, true}, {8, true}}},    // right
  };

  std::array<FaceId, 6> faces{};
  for (int f = 0; f < 6; ++f) {
    const FaceSpec& spec = face_specs[f];
    const Point3d corner = corners[spec.corner];
    const std::uint32_t surface = editor.add_surface(
        Plane3{spec.normal, dot(spec.normal, Vec3d{corner[0], corner[1], corner[2]})});
    const Result<FaceId> face = editor.add_face(surface, false, tol);
    if (!face.has_value()) return face.error();
    faces[static_cast<std::size_t>(f)] = face.value();

    std::array<CoedgeSpec, 4> coedges{};
    for (int i = 0; i < 4; ++i) {
      coedges[static_cast<std::size_t>(i)] =
          CoedgeSpec{edges[spec.loop[i].edge], spec.loop[i].reversed};
    }
    const Result<LoopId> loop = editor.add_loop(face.value(), coedges, true);
    if (!loop.has_value()) return loop.error();
  }

  const Result<ShellId> shell = editor.add_shell(faces, true);
  if (!shell.has_value()) return shell.error();
  const std::array<ShellId, 1> shells{shell.value()};
  const Result<SolidId> solid = editor.add_solid(shells, tol);
  if (!solid.has_value()) return solid.error();
  return editor.commit();
}

// Right circular cylinder: the circle of `radius` around `base`, extruded along
// the unit `axis` by `height`; `tol` is stored on every entity.
//
// The lateral surface is cut along u = 0, so the bottom and top circles are
// each a single closed edge whose only vertex is the seam point, and the seam
// is a line edge used twice by the lateral loop (up and down). The three faces
// are the bottom cap (plane normal -axis), the top cap (plane normal +axis) and
// the lateral face (Cylinder3), and the shell is closed with V = 2, E = 3,
// F = 3 (chi = 2).
//
// Frame: x is the unit perpendicular to axis chosen by the brief's candidate
// rule, y = cross(axis, x), so (x, y, axis) is right handed.
//
// Validation:
//   - non-finite base, axis, radius or height  -> Errc::non_finite
//   - non-finite or negative tolerance         -> Errc::tolerance_violation
//   - axis not unit length, radius <= 0 or
//     height <= 0                              -> Errc::degenerate
//   - top/seam points that round to base or to
//     each other (zero-length seam)             -> Errc::degenerate
// Any error from a SolidEditor add_* call is propagated unchanged; commit()
// performs the final validation.
inline Result<Brep> make_cylinder(Point3d base, Vec3d axis, double radius, double height,
                                  double tol) {
  if (!detail::primitive_point_finite(base) || !detail::primitive_vec_finite(axis) ||
      !std::isfinite(radius) || !std::isfinite(height)) {
    return Error{Errc::non_finite, "make_cylinder: input is not finite", 0u, 0u};
  }
  if (!std::isfinite(tol) || tol < 0.0) {
    return Error{Errc::tolerance_violation,
                 "make_cylinder: tolerance must be finite and non-negative", 0u, 0u};
  }
  const double axis_length = std::hypot(axis[0], axis[1], axis[2]);
  if (!(std::fabs(axis_length - 1.0) <= detail::kAnalyticSurfaceUnitTolerance)) {
    return Error{Errc::degenerate, "make_cylinder: axis must be unit length", 0u, 0u};
  }
  if (!(radius > 0.0)) {
    return Error{Errc::degenerate, "make_cylinder: radius must be positive", 0u, 0u};
  }
  if (!(height > 0.0)) {
    return Error{Errc::degenerate, "make_cylinder: height must be positive", 0u, 0u};
  }

  // Any candidate works as long as it is not parallel to the axis.
  const Vec3d candidate =
      std::fabs(axis[0]) < 0.9 ? Vec3d{1.0, 0.0, 0.0} : Vec3d{0.0, 1.0, 0.0};
  const Vec3d perpendicular = cross(cross(axis, candidate), axis);
  const double perpendicular_length =
      std::hypot(perpendicular[0], perpendicular[1], perpendicular[2]);
  if (!std::isfinite(perpendicular_length) || !(perpendicular_length > 0.0)) {
    return Error{Errc::degenerate, "make_cylinder: axis frame is degenerate", 0u, 0u};
  }
  const Vec3d x_axis{perpendicular[0] / perpendicular_length,
                     perpendicular[1] / perpendicular_length,
                     perpendicular[2] / perpendicular_length};
  const Vec3d y_axis = cross(axis, x_axis);

  const Point3d top{base[0] + height * axis[0], base[1] + height * axis[1],
                    base[2] + height * axis[2]};
  if (!detail::primitive_point_finite(top)) {
    return Error{Errc::degenerate, "make_cylinder: top center is not representable", 0u, 0u};
  }
  if (detail::primitive_point_equal(top, base)) {
    return Error{Errc::degenerate, "make_cylinder: seam has zero length", 0u, 0u};
  }
  const Point3d p0{base[0] + radius * x_axis[0], base[1] + radius * x_axis[1],
                   base[2] + radius * x_axis[2]};
  const Point3d p1{top[0] + radius * x_axis[0], top[1] + radius * x_axis[1],
                   top[2] + radius * x_axis[2]};
  if (!detail::primitive_point_finite(p0) || !detail::primitive_point_finite(p1)) {
    return Error{Errc::degenerate, "make_cylinder: seam point is not representable", 0u, 0u};
  }
  if (detail::primitive_point_equal(p0, p1)) {
    return Error{Errc::degenerate, "make_cylinder: seam has zero length", 0u, 0u};
  }

  constexpr double kTwoPi = 2.0 * std::numbers::pi_v<double>;

  SolidEditor editor;
  const Result<VertexId> bottom_vertex = editor.add_vertex(p0, tol);
  if (!bottom_vertex.has_value()) return bottom_vertex.error();
  const Result<VertexId> top_vertex = editor.add_vertex(p1, tol);
  if (!top_vertex.has_value()) return top_vertex.error();

  const std::uint32_t bottom_circle =
      editor.add_curve(CircleCurve3{base, x_axis, y_axis, radius});
  const Result<EdgeId> bottom_edge = editor.add_edge(
      bottom_vertex.value(), bottom_vertex.value(), bottom_circle, 0.0, kTwoPi, tol);
  if (!bottom_edge.has_value()) return bottom_edge.error();

  const std::uint32_t top_circle = editor.add_curve(CircleCurve3{top, x_axis, y_axis, radius});
  const Result<EdgeId> top_edge = editor.add_edge(
      top_vertex.value(), top_vertex.value(), top_circle, 0.0, kTwoPi, tol);
  if (!top_edge.has_value()) return top_edge.error();

  const std::uint32_t seam_curve = editor.add_curve(LineCurve3{p0, axis});
  const Result<EdgeId> seam_edge = editor.add_edge(bottom_vertex.value(), top_vertex.value(),
                                                   seam_curve, 0.0, height, tol);
  if (!seam_edge.has_value()) return seam_edge.error();

  const Vec3d down{-axis[0], -axis[1], -axis[2]};
  const Vec3d base_vector{base[0], base[1], base[2]};
  const Vec3d top_vector{top[0], top[1], top[2]};
  const std::uint32_t bottom_surface = editor.add_surface(Plane3{down, dot(down, base_vector)});
  const Result<FaceId> bottom_face = editor.add_face(bottom_surface, false, tol);
  if (!bottom_face.has_value()) return bottom_face.error();
  const std::array<CoedgeSpec, 1> bottom_loop{CoedgeSpec{bottom_edge.value(), true}};
  const Result<LoopId> bottom_loop_id = editor.add_loop(bottom_face.value(), bottom_loop, true);
  if (!bottom_loop_id.has_value()) return bottom_loop_id.error();

  const std::uint32_t top_surface = editor.add_surface(Plane3{axis, dot(axis, top_vector)});
  const Result<FaceId> top_face = editor.add_face(top_surface, false, tol);
  if (!top_face.has_value()) return top_face.error();
  const std::array<CoedgeSpec, 1> top_loop{CoedgeSpec{top_edge.value(), false}};
  const Result<LoopId> top_loop_id = editor.add_loop(top_face.value(), top_loop, true);
  if (!top_loop_id.has_value()) return top_loop_id.error();

  const std::uint32_t lateral_surface =
      editor.add_surface(Cylinder3{base, axis, x_axis, radius});
  const Result<FaceId> lateral_face = editor.add_face(lateral_surface, false, tol);
  if (!lateral_face.has_value()) return lateral_face.error();
  const std::array<CoedgeSpec, 4> lateral_loop{
      CoedgeSpec{bottom_edge.value(), false}, CoedgeSpec{seam_edge.value(), false},
      CoedgeSpec{top_edge.value(), true}, CoedgeSpec{seam_edge.value(), true}};
  const Result<LoopId> lateral_loop_id =
      editor.add_loop(lateral_face.value(), lateral_loop, true);
  if (!lateral_loop_id.has_value()) return lateral_loop_id.error();

  const std::array<FaceId, 3> faces{bottom_face.value(), top_face.value(), lateral_face.value()};
  const Result<ShellId> shell = editor.add_shell(faces, true);
  if (!shell.has_value()) return shell.error();
  const std::array<ShellId, 1> shells{shell.value()};
  const Result<SolidId> solid = editor.add_solid(shells, tol);
  if (!solid.has_value()) return solid.error();
  return editor.commit();
}

// Right circular cone: the circle of `base_radius` around `base`, tapering to
// the apex at base + height * axis. The half angle is derived as
// atan2(base_radius, height) rather than taken as an independent parameter, so
// the apex is exactly where the Cone3 lateral surface closes and the base
// circle is exactly its v = 0 cross-section; no second parameter can contradict
// the first.
//
// The lateral surface is cut along u = 0. The base circle is a single closed
// edge whose only vertex is the seam point base + base_radius * x, the seam is
// the slant line from that point to the apex (used twice, up and down), and the
// two faces are the base cap (plane normal -axis, circle reversed) and the
// lateral face (Cone3). V = 2, E = 2, F = 2 (chi = 2).
//
// Frame: the same candidate rule as make_cylinder, x is the unit perpendicular
// to axis and y = cross(axis, x), so (x, y, axis) is right handed.
//
// Validation:
//   - non-finite base, axis, base_radius or height -> Errc::non_finite
//   - non-finite or negative tolerance            -> Errc::tolerance_violation
//   - axis not unit length, base_radius <= 0 or
//     height <= 0                                 -> Errc::degenerate
//   - a derived half angle that rounds to zero    -> Errc::degenerate
//   - an apex or seam point that is not
//     representable, or a seam whose endpoints
//     round to one point                          -> Errc::degenerate
// Any error from a SolidEditor add_* call is propagated unchanged; commit()
// performs the final validation.
inline Result<Brep> make_cone(Point3d base, Vec3d axis, double base_radius, double height,
                              double tol) {
  if (!detail::primitive_point_finite(base) || !detail::primitive_vec_finite(axis) ||
      !std::isfinite(base_radius) || !std::isfinite(height)) {
    return Error{Errc::non_finite, "make_cone: input is not finite", 0u, 0u};
  }
  if (!std::isfinite(tol) || tol < 0.0) {
    return Error{Errc::tolerance_violation,
                 "make_cone: tolerance must be finite and non-negative", 0u, 0u};
  }
  const double axis_length = std::hypot(axis[0], axis[1], axis[2]);
  if (!(std::fabs(axis_length - 1.0) <= detail::kAnalyticSurfaceUnitTolerance)) {
    return Error{Errc::degenerate, "make_cone: axis must be unit length", 0u, 0u};
  }
  if (!(base_radius > 0.0)) {
    return Error{Errc::degenerate, "make_cone: base radius must be positive", 0u, 0u};
  }
  if (!(height > 0.0)) {
    return Error{Errc::degenerate, "make_cone: height must be positive", 0u, 0u};
  }
  const double half_angle = std::atan2(base_radius, height);
  if (!(half_angle > 0.0)) {
    return Error{Errc::degenerate, "make_cone: half angle is not representable", 0u, 0u};
  }

  // Any candidate works as long as it is not parallel to the axis.
  const Vec3d candidate =
      std::fabs(axis[0]) < 0.9 ? Vec3d{1.0, 0.0, 0.0} : Vec3d{0.0, 1.0, 0.0};
  const Vec3d perpendicular = cross(cross(axis, candidate), axis);
  const double perpendicular_length =
      std::hypot(perpendicular[0], perpendicular[1], perpendicular[2]);
  if (!std::isfinite(perpendicular_length) || !(perpendicular_length > 0.0)) {
    return Error{Errc::degenerate, "make_cone: axis frame is degenerate", 0u, 0u};
  }
  const Vec3d x_axis{perpendicular[0] / perpendicular_length,
                     perpendicular[1] / perpendicular_length,
                     perpendicular[2] / perpendicular_length};

  const Point3d apex{base[0] + height * axis[0], base[1] + height * axis[1],
                     base[2] + height * axis[2]};
  if (!detail::primitive_point_finite(apex)) {
    return Error{Errc::degenerate, "make_cone: apex is not representable", 0u, 0u};
  }
  const Point3d p0{base[0] + base_radius * x_axis[0], base[1] + base_radius * x_axis[1],
                   base[2] + base_radius * x_axis[2]};
  if (!detail::primitive_point_finite(p0)) {
    return Error{Errc::degenerate, "make_cone: seam point is not representable", 0u, 0u};
  }
  if (detail::primitive_point_equal(p0, apex)) {
    return Error{Errc::degenerate, "make_cone: seam has zero length", 0u, 0u};
  }
  // The slant direction is formed from the parameters, not from the rounded
  // endpoints, so the LineCurve3 is exactly the u = 0 generatrix of Cone3.
  const Vec3d delta{height * axis[0] - base_radius * x_axis[0],
                    height * axis[1] - base_radius * x_axis[1],
                    height * axis[2] - base_radius * x_axis[2]};
  if (!detail::primitive_vec_finite(delta)) {
    return Error{Errc::degenerate, "make_cone: seam is not representable", 0u, 0u};
  }
  const double slant = std::hypot(delta[0], delta[1], delta[2]);
  if (!std::isfinite(slant) || !(slant > 0.0)) {
    return Error{Errc::degenerate, "make_cone: seam has zero length", 0u, 0u};
  }
  const Vec3d seam_direction{delta[0] / slant, delta[1] / slant, delta[2] / slant};

  constexpr double kTwoPi = 2.0 * std::numbers::pi_v<double>;

  SolidEditor editor;
  const Result<VertexId> base_vertex = editor.add_vertex(p0, tol);
  if (!base_vertex.has_value()) return base_vertex.error();
  const Result<VertexId> apex_vertex = editor.add_vertex(apex, tol);
  if (!apex_vertex.has_value()) return apex_vertex.error();

  const Vec3d y_axis = cross(axis, x_axis);
  const std::uint32_t base_circle =
      editor.add_curve(CircleCurve3{base, x_axis, y_axis, base_radius});
  const Result<EdgeId> base_edge = editor.add_edge(
      base_vertex.value(), base_vertex.value(), base_circle, 0.0, kTwoPi, tol);
  if (!base_edge.has_value()) return base_edge.error();

  const std::uint32_t seam_curve = editor.add_curve(LineCurve3{p0, seam_direction});
  const Result<EdgeId> seam_edge = editor.add_edge(base_vertex.value(), apex_vertex.value(),
                                                   seam_curve, 0.0, slant, tol);
  if (!seam_edge.has_value()) return seam_edge.error();

  const Vec3d down{-axis[0], -axis[1], -axis[2]};
  const Vec3d base_vector{base[0], base[1], base[2]};
  const std::uint32_t cap_surface = editor.add_surface(Plane3{down, dot(down, base_vector)});
  const Result<FaceId> cap_face = editor.add_face(cap_surface, false, tol);
  if (!cap_face.has_value()) return cap_face.error();
  const std::array<CoedgeSpec, 1> cap_loop{CoedgeSpec{base_edge.value(), true}};
  const Result<LoopId> cap_loop_id = editor.add_loop(cap_face.value(), cap_loop, true);
  if (!cap_loop_id.has_value()) return cap_loop_id.error();

  const std::uint32_t lateral_surface =
      editor.add_surface(Cone3{base, axis, x_axis, base_radius, half_angle});
  const Result<FaceId> lateral_face = editor.add_face(lateral_surface, false, tol);
  if (!lateral_face.has_value()) return lateral_face.error();
  const std::array<CoedgeSpec, 3> lateral_loop{
      CoedgeSpec{base_edge.value(), false}, CoedgeSpec{seam_edge.value(), false},
      CoedgeSpec{seam_edge.value(), true}};
  const Result<LoopId> lateral_loop_id =
      editor.add_loop(lateral_face.value(), lateral_loop, true);
  if (!lateral_loop_id.has_value()) return lateral_loop_id.error();

  const std::array<FaceId, 2> faces{cap_face.value(), lateral_face.value()};
  const Result<ShellId> shell = editor.add_shell(faces, true);
  if (!shell.has_value()) return shell.error();
  const std::array<ShellId, 1> shells{shell.value()};
  const Result<SolidId> solid = editor.add_solid(shells, tol);
  if (!solid.has_value()) return solid.error();
  return editor.commit();
}

// Sphere of `radius` around `center`, cut along the meridian half great circle
// in the plane that contains the canonical x axis (1, 0, 0) and axis
// (0, 0, 1). The seam is CircleCurve3{center, x_axis, axis, radius} over
// [-pi/2, pi/2], running from the south pole to the north pole, and the single
// face's loop uses it once in each direction. V = 2, E = 1, F = 1 (chi = 2).
//
// The frame is canonical because make_sphere has no frame parameter: axis is
// (0, 0, 1) and x_axis is (1, 0, 0), so Sphere3::eval(0, v) is exactly the
// seam curve and the poles are the two vertices.
//
// Validation:
//   - non-finite center or radius      -> Errc::non_finite
//   - non-finite or negative tolerance -> Errc::tolerance_violation
//   - radius <= 0                      -> Errc::degenerate
//   - a pole that is not representable,
//     or two poles rounding to one point -> Errc::degenerate
// Any error from a SolidEditor add_* call is propagated unchanged; commit()
// performs the final validation.
inline Result<Brep> make_sphere(Point3d center, double radius, double tol) {
  if (!detail::primitive_point_finite(center) || !std::isfinite(radius)) {
    return Error{Errc::non_finite, "make_sphere: input is not finite", 0u, 0u};
  }
  if (!std::isfinite(tol) || tol < 0.0) {
    return Error{Errc::tolerance_violation,
                 "make_sphere: tolerance must be finite and non-negative", 0u, 0u};
  }
  if (!(radius > 0.0)) {
    return Error{Errc::degenerate, "make_sphere: radius must be positive", 0u, 0u};
  }

  const Vec3d axis{0.0, 0.0, 1.0};
  const Vec3d x_axis{1.0, 0.0, 0.0};
  const Point3d south{center[0] - radius * axis[0], center[1] - radius * axis[1],
                      center[2] - radius * axis[2]};
  const Point3d north{center[0] + radius * axis[0], center[1] + radius * axis[1],
                      center[2] + radius * axis[2]};
  if (!detail::primitive_point_finite(south) || !detail::primitive_point_finite(north)) {
    return Error{Errc::degenerate, "make_sphere: pole is not representable", 0u, 0u};
  }
  if (detail::primitive_point_equal(south, north)) {
    return Error{Errc::degenerate, "make_sphere: seam has zero length", 0u, 0u};
  }

  constexpr double kHalfPi = std::numbers::pi_v<double> / 2.0;

  SolidEditor editor;
  const Result<VertexId> south_vertex = editor.add_vertex(south, tol);
  if (!south_vertex.has_value()) return south_vertex.error();
  const Result<VertexId> north_vertex = editor.add_vertex(north, tol);
  if (!north_vertex.has_value()) return north_vertex.error();

  const std::uint32_t seam_curve =
      editor.add_curve(CircleCurve3{center, x_axis, axis, radius});
  const Result<EdgeId> seam_edge = editor.add_edge(
      south_vertex.value(), north_vertex.value(), seam_curve, -kHalfPi, kHalfPi, tol);
  if (!seam_edge.has_value()) return seam_edge.error();

  const std::uint32_t sphere_surface =
      editor.add_surface(Sphere3{center, axis, x_axis, radius});
  const Result<FaceId> face = editor.add_face(sphere_surface, false, tol);
  if (!face.has_value()) return face.error();
  const std::array<CoedgeSpec, 2> loop{CoedgeSpec{seam_edge.value(), false},
                                        CoedgeSpec{seam_edge.value(), true}};
  const Result<LoopId> loop_id = editor.add_loop(face.value(), loop, true);
  if (!loop_id.has_value()) return loop_id.error();

  const std::array<FaceId, 1> faces{face.value()};
  const Result<ShellId> shell = editor.add_shell(faces, true);
  if (!shell.has_value()) return shell.error();
  const std::array<ShellId, 1> shells{shell.value()};
  const Result<SolidId> solid = editor.add_solid(shells, tol);
  if (!solid.has_value()) return solid.error();
  return editor.commit();
}

}  // namespace v1
}  // namespace geom
