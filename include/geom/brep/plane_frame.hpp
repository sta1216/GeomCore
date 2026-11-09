#pragma once

#include <cmath>

#include <geom/core/result.hpp>
#include <geom/surfaces/plane.hpp>

namespace geom {
inline namespace v1 {

// Orthonormal right-handed frame anchored on a plane: u x v = normal.
// The basis depends only on the normal direction, so scaling a plane's normal
// and offset together leaves the frame unchanged and every query against the
// same geometric plane shares one frame.
struct PlaneFrame {
  Point3d origin;
  Vec3d u;
  Vec3d v;
  Vec3d normal;
};

namespace detail {

// |v| must be nonzero; returns v / |v| computed with the scale-safe hypot.
inline Vec3d plane_frame_normalize(const Vec3d& v) noexcept {
  const double h = std::hypot(v[0], v[1], v[2]);
  return Vec3d{v[0] / h, v[1] / h, v[2] / h};
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

// Validates the plane normal and builds its frame. A non-finite normal
// component is rejected as non_finite; a zero normal as degenerate. The plane
// offset never participates in validation: the resulting frame is valid for
// any offset.
inline Result<PlaneFrame> make_plane_frame(const Plane3& plane) {
  const Vec3d& raw = plane.normal;
  if (!std::isfinite(raw[0]) || !std::isfinite(raw[1]) || !std::isfinite(raw[2])) {
    return Error{Errc::non_finite, "plane normal has a non-finite component"};
  }
  const double h = std::hypot(raw[0], raw[1], raw[2]);
  if (h == 0.0) {
    return Error{Errc::degenerate, "plane normal is zero"};
  }
  const Vec3d normal{raw[0] / h, raw[1] / h, raw[2] / h};
  const double distance = plane.offset / h;

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

// {unit normal, offset / |normal|}, i.e. the same geometric plane rescaled to
// unit normal. dot(normal, origin) equals offset / |normal| because origin is
// the normal-foot point nearest the world origin.
inline Plane3 plane_frame_plane(const PlaneFrame& frame) noexcept {
  const Vec3d& n = frame.normal;
  const Point3d& o = frame.origin;
  return Plane3{n, n[0] * o[0] + n[1] * o[1] + n[2] * o[2]};
}

}  // namespace v1
}  // namespace geom
