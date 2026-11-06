#pragma once

#include <geom/primitives/point.hpp>
#include <geom/primitives/vec.hpp>

namespace geom {
inline namespace v1 {

// Plane n . x = offset with a user-supplied normal. The struct does not
// normalize the normal; callers who want distances in world units must pass a
// unit normal.
struct Plane3 {
  Vec3d normal{0.0, 0.0, 1.0};  // unit normal
  double offset = 0.0;          // n . x = offset

  double signed_distance(Point3d p) const noexcept {
    return dot(normal, Vec3d{p[0], p[1], p[2]}) - offset;
  }

  int side(Point3d p) const noexcept {
    const double d = signed_distance(p);
    return d > 0.0 ? 1 : (d < 0.0 ? -1 : 0);
  }

  Point3d project(Point3d p) const noexcept {
    const double d = signed_distance(p);
    return Point3d{p[0] - d * normal[0], p[1] - d * normal[1], p[2] - d * normal[2]};
  }
};

}  // namespace v1
}  // namespace geom
