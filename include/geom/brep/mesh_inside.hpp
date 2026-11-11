#pragma once

// Point classification against a closed triangulation.
//
// The boolean driver classifies fragment representative points by ray parity.
// A ray that grazes an edge or a vertex cannot prove inside/outside, so four
// independent directions are tried; if all four are degenerate the query fails
// instead of guessing. A point within `tolerance` of the mesh surface also
// fails: classifying a point on the boundary is caller-ambiguous.

#include <array>
#include <cmath>
#include <cstddef>
#include <limits>
#include <vector>

#include <geom/brep/tessellate.hpp>
#include <geom/core/error.hpp>
#include <geom/primitives/point.hpp>
#include <geom/primitives/vec.hpp>

namespace geom {
namespace v1 {
namespace detail {
namespace mesh_inside_detail {

inline Vec3d mesh_inside_cross(const Vec3d& a, const Vec3d& b) noexcept {
  return Vec3d{a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2],
               a[0] * b[1] - a[1] * b[0]};
}

inline double mesh_inside_dot(const Vec3d& a, const Vec3d& b) noexcept {
  return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
}

// Squared distance from p to the triangle (a, b, c). Ericson, Real-Time
// Collision Detection, 5.1.5.
inline double mesh_inside_triangle_distance2(const Point3d& p, const Point3d& a, const Point3d& b,
                                             const Point3d& c) noexcept {
  const Vec3d ab = b - a;
  const Vec3d ac = c - a;
  const Vec3d ap = p - a;
  const double d1 = mesh_inside_dot(ab, ap);
  const double d2 = mesh_inside_dot(ac, ap);
  if (d1 <= 0.0 && d2 <= 0.0) return mesh_inside_dot(ap, ap);

  const Vec3d bp = p - b;
  const double d3 = mesh_inside_dot(ab, bp);
  const double d4 = mesh_inside_dot(ac, bp);
  if (d3 >= 0.0 && d4 <= d3) return mesh_inside_dot(bp, bp);

  const double vc = d1 * d4 - d3 * d2;
  if (vc <= 0.0 && d1 >= 0.0 && d3 <= 0.0) {
    const double v = d1 / (d1 - d3);
    const Vec3d q = ab * v;
    const Vec3d r = ap - q;
    return mesh_inside_dot(r, r);
  }

  const Vec3d cp = p - c;
  const double d5 = mesh_inside_dot(ab, cp);
  const double d6 = mesh_inside_dot(ac, cp);
  if (d6 >= 0.0 && d5 <= d6) return mesh_inside_dot(cp, cp);

  const double vb = d5 * d2 - d1 * d6;
  if (vb <= 0.0 && d2 >= 0.0 && d6 <= 0.0) {
    const double w = d2 / (d2 - d6);
    const Vec3d q = ac * w;
    const Vec3d r = ap - q;
    return mesh_inside_dot(r, r);
  }

  const double va = d3 * d6 - d5 * d4;
  if (va <= 0.0 && (d4 - d3) >= 0.0 && (d5 - d6) >= 0.0) {
    const double w = (d4 - d3) / ((d4 - d3) + (d5 - d6));
    const Vec3d edge = c - b;
    const Vec3d r = p - (b + edge * w);
    return mesh_inside_dot(r, r);
  }

  const double denominator = 1.0 / (va + vb + vc);
  const double v = vb * denominator;
  const double w = vc * denominator;
  const Vec3d q = ab * v + ac * w;
  const Vec3d r = ap - q;
  return mesh_inside_dot(r, r);
}

}  // namespace mesh_inside_detail
}  // namespace detail

// Returns true when `p` is inside the closed mesh. Fails when `p` lies on the
// surface within `tolerance`, when a probe ray grazes an edge/vertex in all four
// directions, or when the input is invalid.
inline Result<bool> point_in_tessellation(const TriMesh& mesh, const Point3d& p,
                                          double tolerance) {
  if (!std::isfinite(p[0]) || !std::isfinite(p[1]) || !std::isfinite(p[2])) {
    return Error{Errc::non_finite, "point_in_tessellation: point is not finite", 0u, 0u};
  }
  if (!(tolerance > 0.0) || !std::isfinite(tolerance)) {
    return Error{Errc::tolerance_violation,
                 "point_in_tessellation: tolerance must be positive and finite", 0u, 0u};
  }
  if (mesh.positions.empty() || mesh.triangles.empty()) {
    return Error{Errc::degenerate, "point_in_tessellation: mesh is empty", 0u, 0u};
  }

  const double tolerance2 = tolerance * tolerance;
  for (const std::array<std::uint32_t, 3>& triangle : mesh.triangles) {
    if (triangle[0] >= mesh.positions.size() || triangle[1] >= mesh.positions.size() ||
        triangle[2] >= mesh.positions.size()) {
      return Error{Errc::invalid_topology, "point_in_tessellation: triangle index is out of range",
                   0u, 0u};
    }
    const Point3d& a = mesh.positions[triangle[0]];
    const Point3d& b = mesh.positions[triangle[1]];
    const Point3d& c = mesh.positions[triangle[2]];
    if (detail::mesh_inside_detail::mesh_inside_triangle_distance2(p, a, b, c) <= tolerance2) {
      return Error{Errc::numeric_failure, "point_in_tessellation: point is on the mesh surface",
                   0u, 0u};
    }
  }

  const double inv14 = 1.0 / std::sqrt(14.0);
  const std::array<Vec3d, 6> directions = {
      Vec3d{1.0, 0.0, 0.0},
      Vec3d{0.0, 1.0, 0.0},
      Vec3d{0.0, 0.0, 1.0},
      Vec3d{1.0 * inv14, 2.0 * inv14, 3.0 * inv14},
      Vec3d{3.0 * inv14, 1.0 * inv14, 2.0 * inv14},
      Vec3d{2.0 * inv14, 3.0 * inv14, 1.0 * inv14},
  };

  for (const Vec3d& direction : directions) {
    std::size_t crossings = 0;
    bool ambiguous = false;
    for (const std::array<std::uint32_t, 3>& triangle : mesh.triangles) {
      const Point3d& a = mesh.positions[triangle[0]];
      const Point3d& b = mesh.positions[triangle[1]];
      const Point3d& c = mesh.positions[triangle[2]];
      const Vec3d e1 = b - a;
      const Vec3d e2 = c - a;
      const Vec3d pvec = detail::mesh_inside_detail::mesh_inside_cross(direction, e2);
      const double det = detail::mesh_inside_detail::mesh_inside_dot(e1, pvec);
      if (!(std::abs(det) > 0.0)) continue;
      const double inv = 1.0 / det;
      const Vec3d tvec = p - a;
      const double u = detail::mesh_inside_detail::mesh_inside_dot(tvec, pvec) * inv;
      const Vec3d qvec = detail::mesh_inside_detail::mesh_inside_cross(tvec, e1);
      const double v = detail::mesh_inside_detail::mesh_inside_dot(direction, qvec) * inv;
      const double t = detail::mesh_inside_detail::mesh_inside_dot(e2, qvec) * inv;
      if (!(t > 0.0)) continue;

      const double scale = 1.0 / (std::sqrt(detail::mesh_inside_detail::mesh_inside_dot(e1, e1)) +
                                  std::sqrt(detail::mesh_inside_detail::mesh_inside_dot(e2, e2)) +
                                  std::sqrt(detail::mesh_inside_detail::mesh_inside_dot(c - b, c - b)));
      const double eps = tolerance * scale;
      if (u < -eps || v < -eps || u + v > 1.0 + eps) continue;
      if (u <= eps || v <= eps || u + v >= 1.0 - eps) {
        ambiguous = true;
        break;
      }
      ++crossings;
    }
    if (ambiguous) continue;
    return (crossings % 2u) == 1u;
  }
  return Error{Errc::numeric_failure, "point_in_tessellation: every probe ray is degenerate", 0u,
               0u};
}

}  // namespace v1
}  // namespace geom
