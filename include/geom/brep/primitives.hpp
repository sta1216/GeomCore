#pragma once

#include <array>
#include <cmath>
#include <cstdint>

#include <geom/brep/solid_editor.hpp>
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

inline bool primitive_point_finite(Point3d p) noexcept {
  return std::isfinite(p[0]) && std::isfinite(p[1]) && std::isfinite(p[2]);
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

}  // namespace v1
}  // namespace geom
