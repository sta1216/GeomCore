#pragma once

#include <algorithm>
#include <cmath>

#include <geom/core/result.hpp>
#include <geom/surfaces/plane.hpp>

namespace geom {
inline namespace v1 {

// Orthonormal right-handed frame anchored on an oriented plane: u x v = normal.
// The basis depends only on the normal direction, so positively scaling a
// plane's normal and offset together leaves the frame unchanged; a negative
// scale flips both normal and u while the origin stays put. Every query against
// the same oriented plane shares one frame.
struct PlaneFrame {
  Point3d origin;
  Vec3d u;
  Vec3d v;
  Vec3d normal;
};

namespace detail {

// Normalizes a nonzero finite vector. Scales by the largest magnitude
// component first, so the intermediate length stays in [1, sqrt(3)] and
// cannot overflow or underflow even for denormal or near-maximum inputs.
inline Vec3d plane_frame_normalize(const Vec3d& v) noexcept {
  const double ax = std::fabs(v[0]);
  const double ay = std::fabs(v[1]);
  const double az = std::fabs(v[2]);
  const double m = std::max({ax, ay, az});
  const double x = v[0] / m;
  const double y = v[1] / m;
  const double z = v[2] / m;
  const double h = std::sqrt(x * x + y * y + z * z);
  return Vec3d{x / h, y / h, z / h};
}

// Builds an orthonormal right-handed basis (u, v) for the UNIT normal.
// The reference axis is the coordinate axis LEAST aligned with `normal`
// (smallest |normal[i]|). Using the largest component instead would make
// cross(ref, normal) vanish for axis-aligned normals such as (0, 1, 0),
// whose frame the planar work must support.
inline void plane_frame_basis(const Vec3d& normal, Vec3d& u, Vec3d& v) {
  int axis = 0;
  double smallest = std::fabs(normal[0]);
  for (int i = 1; i < 3; ++i) {
    const double component = std::fabs(normal[i]);
    if (component < smallest) {
      smallest = component;
      axis = i;
    }
  }
  Vec3d ref{};
  ref[axis] = 1.0;
  u = plane_frame_normalize(cross(ref, normal));
  v = cross(normal, u);
}

}  // namespace detail

// Validates the plane and builds its frame. A non-finite normal component or
// a non-finite offset is rejected as non_finite; a zero normal as degenerate.
// The distance offset / |normal| is computed scale-safely and must be finite,
// otherwise the plane is rejected as non_finite. If that distance rounds to
// zero because it lies below the smallest subnormal, the frame is accepted and
// origin is the nearest representable point (the world origin).
inline Result<PlaneFrame> make_plane_frame(const Plane3& plane) {
  const Vec3d& raw = plane.normal;
  if (!std::isfinite(raw[0]) || !std::isfinite(raw[1]) || !std::isfinite(raw[2])) {
    return Error{Errc::non_finite, "plane normal has a non-finite component"};
  }
  if (!std::isfinite(plane.offset)) {
    return Error{Errc::non_finite, "plane offset is not finite"};
  }
  const double ax = std::fabs(raw[0]);
  const double ay = std::fabs(raw[1]);
  const double az = std::fabs(raw[2]);
  const double m = std::max({ax, ay, az});
  if (m == 0.0) {
    return Error{Errc::degenerate, "plane normal is zero"};
  }
  const double x = raw[0] / m;
  const double y = raw[1] / m;
  const double z = raw[2] / m;
  const double h = std::sqrt(x * x + y * y + z * z);
  const Vec3d normal{x / h, y / h, z / h};
  // Divide by h in [1, sqrt(3)] first: offset / m can overflow even when the
  // true distance offset / (m * h) is finite.
  const double distance = (plane.offset / h) / m;
  if (!std::isfinite(distance)) {
    return Error{Errc::non_finite, "plane offset is too large for a finite frame"};
  }

  PlaneFrame frame;
  frame.normal = normal;
  detail::plane_frame_basis(normal, frame.u, frame.v);
  frame.origin = Point3d{normal[0] * distance, normal[1] * distance, normal[2] * distance};
  return frame;
}

inline Point2d plane_frame_to_2d(const PlaneFrame& frame, const Point3d& p) noexcept {
  const Vec3d d = p - frame.origin;
  return Point2d{dot(d, frame.u), dot(d, frame.v)};
}

inline Point3d plane_frame_to_3d(const PlaneFrame& frame, const Point2d& q) noexcept {
  return frame.origin + frame.u * q[0] + frame.v * q[1];
}

// {unit normal, frame distance}, i.e. the same oriented plane rescaled to unit
// normal. dot(normal, origin) equals offset / |normal| rounded to the nearest
// representable double (it may be zero when that distance underflows).
inline Plane3 plane_frame_plane(const PlaneFrame& frame) noexcept {
  const Vec3d& n = frame.normal;
  const Point3d& o = frame.origin;
  return Plane3{n, n[0] * o[0] + n[1] * o[1] + n[2] * o[2]};
}

}  // namespace v1
}  // namespace geom
