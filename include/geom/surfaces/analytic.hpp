#pragma once

#include <cmath>
#include <numbers>

#include <geom/primitives/point.hpp>
#include <geom/primitives/vec.hpp>

namespace geom {
inline namespace v1 {

namespace detail {

// Unit-length and orthogonality slack shared by the analytic surface validity
// checks. Like every geometric tolerance in this library it is a local
// constant, not a global epsilon, and it only bounds representation error in
// inputs that are meant to be exactly unit / exactly orthogonal. The analytic
// curves use the same 1e-12 bound (kAnalyticUnitTolerance); the constant is
// repeated here so that this header stays independent of the curve header.
inline constexpr double kAnalyticSurfaceUnitTolerance = 1e-12;

inline bool analytic_surface_point_finite(Point3d p) noexcept {
  return std::isfinite(p[0]) && std::isfinite(p[1]) && std::isfinite(p[2]);
}

inline bool analytic_surface_vec_finite(Vec3d v) noexcept {
  return std::isfinite(v[0]) && std::isfinite(v[1]) && std::isfinite(v[2]);
}

// axis and x_axis must be finite, unit length and orthogonal. The in-plane
// basis is then x_axis and y = cross(axis, x_axis), which is unit and
// orthogonal as well, so (x_axis, y, axis) is a right-handed orthonormal frame
// and every normal written below has unit length without a division.
inline bool analytic_surface_frame_valid(Vec3d axis, Vec3d x_axis) noexcept {
  if (!analytic_surface_vec_finite(axis) || !analytic_surface_vec_finite(x_axis)) {
    return false;
  }
  // std::hypot is scale-safe and returns the norm exactly for axes with zero
  // components; dot(x, x) can underflow or overflow for extreme axes. The
  // tolerance is measured on the norm, not on the squared norm, so an axis with
  // | |x| - 1 | <= 1e-12 is accepted even though |dot(x, x) - 1| can be 2e-12.
  const double n_axis = std::hypot(axis[0], axis[1], axis[2]);
  const double n_x = std::hypot(x_axis[0], x_axis[1], x_axis[2]);
  const double axis_dot_x = dot(axis, x_axis);
  return std::fabs(n_axis - 1.0) <= kAnalyticSurfaceUnitTolerance &&
         std::fabs(n_x - 1.0) <= kAnalyticSurfaceUnitTolerance &&
         std::fabs(axis_dot_x) <= kAnalyticSurfaceUnitTolerance;
}

// radial(u) = cos(u) * x_axis + sin(u) * cross(axis, x_axis), the unit outward
// direction of the frame at angle u. Precondition: the frame is valid.
inline Vec3d analytic_surface_radial(double u, Vec3d axis, Vec3d x_axis) noexcept {
  const Vec3d y_axis = cross(axis, x_axis);
  const double c = std::cos(u);
  const double s = std::sin(u);
  return Vec3d{c * x_axis[0] + s * y_axis[0], c * x_axis[1] + s * y_axis[1],
               c * x_axis[2] + s * y_axis[2]};
}

// d(radial)/du = -sin(u) * x_axis + cos(u) * cross(axis, x_axis), the unit
// angular tangent. Precondition: the frame is valid.
inline Vec3d analytic_surface_tangent(double u, Vec3d axis, Vec3d x_axis) noexcept {
  const Vec3d y_axis = cross(axis, x_axis);
  const double c = std::cos(u);
  const double s = std::sin(u);
  return Vec3d{-s * x_axis[0] + c * y_axis[0], -s * x_axis[1] + c * y_axis[1],
               -s * x_axis[2] + c * y_axis[2]};
}

}  // namespace detail

// Infinite right circular cylinder:
//   p(u, v) = origin + v * axis + radius * radial(u)
// u is the angle around the axis (radians, periodic with 2*pi) and v is the
// signed distance along the axis. axis must be unit length, x_axis must be
// unit length and orthogonal to axis, and radius must be positive; is_valid
// reports violations and never throws. There is no guard on u or v: a
// non-finite parameter simply produces a non-finite point or derivative.
struct Cylinder3 {
  Point3d origin{};
  Vec3d axis{0.0, 0.0, 1.0};
  Vec3d x_axis{1.0, 0.0, 0.0};
  double radius = 1.0;

  bool is_valid() const noexcept {
    return detail::analytic_surface_point_finite(origin) &&
           detail::analytic_surface_frame_valid(axis, x_axis) &&
           std::isfinite(radius) && radius > 0.0;
  }

  Point3d eval(double u, double v) const noexcept {
    return origin + axis * v + detail::analytic_surface_radial(u, axis, x_axis) * radius;
  }

  Vec3d du(double u, double) const noexcept {
    return detail::analytic_surface_tangent(u, axis, x_axis) * radius;
  }

  Vec3d dv(double, double) const noexcept { return axis; }

  // radial(u) is already the unit outward normal, so the normal is independent
  // of v and cannot degenerate for a valid entity.
  Vec3d normal(double u, double) const noexcept {
    return detail::analytic_surface_radial(u, axis, x_axis);
  }
};

// Cone whose v = 0 cross-section is a circle of the given radius:
//   p(u, v) = origin + v * axis + (radius - v * tan(half_angle)) * radial(u)
// The radius shrinks as v grows and reaches zero at v = radius / tan(half_angle)
// (the apex). For v beyond the apex the parameterization continues on the
// opposite nappe of the same double cone: eval() stays on the surface, but the
// (u, v) frame turns over, so normal() flips sign to stay consistent with
// cross(du, dv). Callers that trim the cone at the apex never reach that branch.
// half_angle is the angle between the surface and the axis and must lie in
// (0, pi/2). The frame and radius requirements match Cylinder3.
struct Cone3 {
  Point3d origin{};
  Vec3d axis{0.0, 0.0, 1.0};
  Vec3d x_axis{1.0, 0.0, 0.0};
  double radius = 1.0;
  double half_angle = 0.5;

  bool is_valid() const noexcept {
    return detail::analytic_surface_point_finite(origin) &&
           detail::analytic_surface_frame_valid(axis, x_axis) &&
           std::isfinite(radius) && radius > 0.0 && std::isfinite(half_angle) &&
           half_angle > 0.0 && half_angle < std::numbers::pi_v<double> / 2.0;
  }

  Point3d eval(double u, double v) const noexcept {
    const double rho = radius - v * std::tan(half_angle);
    return origin + axis * v + detail::analytic_surface_radial(u, axis, x_axis) * rho;
  }

  Vec3d du(double u, double v) const noexcept {
    return detail::analytic_surface_tangent(u, axis, x_axis) * (radius - v * std::tan(half_angle));
  }

  Vec3d dv(double u, double) const noexcept {
    return axis - detail::analytic_surface_radial(u, axis, x_axis) * std::tan(half_angle);
  }

  // cross(du, dv) = (radius - v*tan(half_angle)) * (radial + tan(half_angle)*axis),
  // so the unit parameterization normal is
  //   sign(radius - v*tan(half_angle)) *
  //       (cos(half_angle)*radial + sin(half_angle)*axis).
  // The closed form is unit everywhere, including the apex, where du vanishes
  // and cross(du, dv) stops determining a direction; there the limit normal of
  // the v > 0 side is returned.
  Vec3d normal(double u, double v) const noexcept {
    const Vec3d n = detail::analytic_surface_radial(u, axis, x_axis) * std::cos(half_angle) +
                    axis * std::sin(half_angle);
    return radius - v * std::tan(half_angle) < 0.0 ? Vec3d{-n[0], -n[1], -n[2]} : n;
  }
};

// Sphere of the given radius:
//   p(u, v) = center + radius * (cos(v) * radial(u) + sin(v) * axis)
// u is the longitude (periodic with 2*pi), v is the latitude in
// [-pi/2, pi/2]. radius must be positive and the frame must be orthonormal.
struct Sphere3 {
  Point3d center{};
  Vec3d axis{0.0, 0.0, 1.0};
  Vec3d x_axis{1.0, 0.0, 0.0};
  double radius = 1.0;

  bool is_valid() const noexcept {
    return detail::analytic_surface_point_finite(center) &&
           detail::analytic_surface_frame_valid(axis, x_axis) &&
           std::isfinite(radius) && radius > 0.0;
  }

  Point3d eval(double u, double v) const noexcept {
    const double c = std::cos(v);
    const double s = std::sin(v);
    return center + detail::analytic_surface_radial(u, axis, x_axis) * (radius * c) +
           axis * (radius * s);
  }

  Vec3d du(double u, double v) const noexcept {
    return detail::analytic_surface_tangent(u, axis, x_axis) * (radius * std::cos(v));
  }

  Vec3d dv(double u, double v) const noexcept {
    const double c = std::cos(v);
    const double s = std::sin(v);
    return (axis * c - detail::analytic_surface_radial(u, axis, x_axis) * s) * radius;
  }

  // p - center divided by radius, written directly as the unit vector
  // cos(v)*radial + sin(v)*axis. cross(du, dv) is parallel to this vector but
  // vanishes at the poles (as cos(v)); the closed form is the directional
  // limit there, returning +axis at v = +pi/2 and -axis at v = -pi/2. At the
  // poles the radial residual cos(v)*radial(u) is one ulp of 1 (about 6e-17),
  // so equal-to-±axis comparisons need one ulp of slack, not exact equality.
  Vec3d normal(double u, double v) const noexcept {
    const double c = std::cos(v);
    const double s = std::sin(v);
    return detail::analytic_surface_radial(u, axis, x_axis) * c + axis * s;
  }
};

// Torus with major radius R and minor radius r:
//   p(u, v) = center + (R + r*cos(v)) * radial(u) + r*sin(v) * axis
// u goes around the axis, v goes around the tube. Both parameters are angles
// (radians); the surface is smooth and periodic in each. 0 < r < R is required
// so the tube never self-intersects and R + r*cos(v) stays positive.
struct Torus3 {
  Point3d center{};
  Vec3d axis{0.0, 0.0, 1.0};
  Vec3d x_axis{1.0, 0.0, 0.0};
  double major_radius = 2.0;
  double minor_radius = 0.5;

  bool is_valid() const noexcept {
    return detail::analytic_surface_point_finite(center) &&
           detail::analytic_surface_frame_valid(axis, x_axis) &&
           std::isfinite(major_radius) && std::isfinite(minor_radius) && minor_radius > 0.0 &&
           minor_radius < major_radius;
  }

  Point3d eval(double u, double v) const noexcept {
    const double c = std::cos(v);
    const double s = std::sin(v);
    const Vec3d radial = detail::analytic_surface_radial(u, axis, x_axis);
    return center + radial * (major_radius + minor_radius * c) + axis * (minor_radius * s);
  }

  Vec3d du(double u, double v) const noexcept {
    return detail::analytic_surface_tangent(u, axis, x_axis) *
           (major_radius + minor_radius * std::cos(v));
  }

  Vec3d dv(double u, double v) const noexcept {
    const double c = std::cos(v);
    const double s = std::sin(v);
    return axis * (minor_radius * c) -
           detail::analytic_surface_radial(u, axis, x_axis) * (minor_radius * s);
  }

  // cos(v)*radial + sin(v)*axis; unit for a valid entity and never degenerate
  // because R + r*cos(v) > 0 when 0 < r < R.
  Vec3d normal(double u, double v) const noexcept {
    const double c = std::cos(v);
    const double s = std::sin(v);
    return detail::analytic_surface_radial(u, axis, x_axis) * c + axis * s;
  }
};

}  // namespace v1
}  // namespace geom
