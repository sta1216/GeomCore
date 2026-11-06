#pragma once

#include <cmath>

#include <geom/primitives/point.hpp>
#include <geom/primitives/vec.hpp>

namespace geom {
inline namespace v1 {

// Plane n . x = offset. The normal is stored as supplied and is never
// normalized: all queries divide by |n| (or |n|^2 for projection), so they are
// correct for any nonzero normal, and scaling the normal together with the
// offset leaves the plane unchanged.
//
// Precondition: |normal| != 0. Plane3 does not validate the normal, and a zero
// normal makes the queries divide by zero (sign and projected point are then
// unspecified); there is no silent zero-distance fallback.
struct Plane3 {
  Vec3d normal{0.0, 0.0, 1.0};
  double offset = 0.0;  // n . x = offset

  double signed_distance(Point3d p) const noexcept {
    const double n2 = dot(normal, normal);
    return (dot(normal, Vec3d{p[0], p[1], p[2]}) - offset) / std::sqrt(n2);
  }

  int side(Point3d p) const noexcept {
    const double d = signed_distance(p);
    return d > 0.0 ? 1 : (d < 0.0 ? -1 : 0);
  }

  // p - ((n . p - offset) / |n|^2) n, i.e. the closest point on the plane.
  Point3d project(Point3d p) const noexcept {
    const double n2 = dot(normal, normal);
    const double t = (dot(normal, Vec3d{p[0], p[1], p[2]}) - offset) / n2;
    return Point3d{p[0] - t * normal[0], p[1] - t * normal[1], p[2] - t * normal[2]};
  }
};

}  // namespace v1
}  // namespace geom
