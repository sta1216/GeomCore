#pragma once

#include <cmath>

#include <geom/primitives/point.hpp>
#include <geom/primitives/vec.hpp>

namespace geom {
inline namespace v1 {

// Plane n . x = offset. The normal is stored as supplied and is never
// normalized: every query reduces the plane to its unit form (unit = n / |n|,
// c = offset / |n|) so results are independent of the normal's scale, and
// scaling the normal and the offset together leaves the plane unchanged.
//
// |n| is computed with the three-argument std::hypot, which is scale-safe and
// therefore avoids intermediate underflow (dot(n, n) can flush subnormal
// normals to zero) and overflow (dot(n, n) can exceed the double range) for
// extreme but representable normals. offset / |n| can overflow only when the
// true distance magnitude exceeds the double range, in which case an infinite
// distance is the honest result.
//
// Precondition: |normal| != 0. Plane3 does not validate the normal, and a zero
// normal makes the queries divide by zero (sign and projected point are then
// unspecified); there is no silent zero-distance fallback.
struct Plane3 {
  Vec3d normal{0.0, 0.0, 1.0};
  double offset = 0.0;  // n . x = offset

  double signed_distance(Point3d p) const noexcept {
    const double h = length();
    const Vec3d unit{normal[0] / h, normal[1] / h, normal[2] / h};
    return dot(unit, Vec3d{p[0], p[1], p[2]}) - offset / h;
  }

  int side(Point3d p) const noexcept {
    const double d = signed_distance(p);
    return d > 0.0 ? 1 : (d < 0.0 ? -1 : 0);
  }

  // p - (dot(unit, p) - offset / |n|) * unit, i.e. the closest point on the
  // plane for any nonzero normal.
  Point3d project(Point3d p) const noexcept {
    const double h = length();
    const Vec3d unit{normal[0] / h, normal[1] / h, normal[2] / h};
    const double t = dot(unit, Vec3d{p[0], p[1], p[2]}) - offset / h;
    return Point3d{p[0] - t * unit[0], p[1] - t * unit[1], p[2] - t * unit[2]};
  }

 private:
  double length() const noexcept { return std::hypot(normal[0], normal[1], normal[2]); }
};

}  // namespace v1
}  // namespace geom
