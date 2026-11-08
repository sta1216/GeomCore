#pragma once

#include <cassert>
#include <cmath>

#include <geom/primitives/box.hpp>
#include <geom/primitives/point.hpp>
#include <geom/primitives/vec.hpp>

namespace geom {
inline namespace v1 {

namespace detail {

// Unit-length and orthogonality slack shared by the analytic curve validity
// checks. Like every geometric tolerance in this library it is a local
// constant, not a global epsilon, and it only bounds representation error in
// inputs that are meant to be exactly unit / exactly orthogonal.
inline constexpr double kAnalyticUnitTolerance = 1e-12;

inline bool analytic_point_finite(Point3d p) noexcept {
  return std::isfinite(p[0]) && std::isfinite(p[1]) && std::isfinite(p[2]);
}

inline bool analytic_vec_finite(Vec3d v) noexcept {
  return std::isfinite(v[0]) && std::isfinite(v[1]) && std::isfinite(v[2]);
}

}  // namespace detail

// Infinite line p(t) = origin + t * direction, parameterized by arc length:
// direction must be finite and unit length, so eval(t) advances t units along
// the line. There is no is_finite guard on t; an infinite or NaN t simply
// produces an infinite or NaN point.
struct LineCurve3 {
  Point3d origin{};
  Vec3d direction{1.0, 0.0, 0.0};

  bool is_valid() const noexcept {
    if (!detail::analytic_point_finite(origin) || !detail::analytic_vec_finite(direction)) {
      return false;
    }
    const double d2 = dot(direction, direction);
    return std::fabs(1.0 - d2) <= detail::kAnalyticUnitTolerance;
  }

  Point3d eval(double t) const noexcept { return origin + direction * t; }

  // The unit direction is constant, so the derivative does not depend on t.
  Vec3d derivative(double) const noexcept { return direction; }

  // Precondition: t0 <= t1. The arc of a straight line is bounded by its two
  // endpoints, which is exact.
  Box3d bounds(double t0, double t1) const noexcept {
    assert(t0 <= t1);
    const Point3d p0 = eval(t0);
    const Point3d p1 = eval(t1);
    Box3d box;
    for (int i = 0; i < 3; ++i) {
      box.lo[i] = p0[i] < p1[i] ? p0[i] : p1[i];
      box.hi[i] = p0[i] < p1[i] ? p1[i] : p0[i];
    }
    return box;
  }
};

// Circle p(angle) = center + radius * (cos(angle) * x_axis + sin(angle) * y_axis).
// x_axis and y_axis span the circle plane and must be finite, unit length and
// orthogonal; radius must be finite and positive. The parameter is the angle in
// radians, so eval is periodic with 2*pi.
struct CircleCurve3 {
  Point3d center{};
  Vec3d x_axis{1.0, 0.0, 0.0};
  Vec3d y_axis{0.0, 1.0, 0.0};
  double radius = 1.0;

  bool is_valid() const noexcept {
    if (!detail::analytic_point_finite(center) || !detail::analytic_vec_finite(x_axis) ||
        !detail::analytic_vec_finite(y_axis) || !std::isfinite(radius)) {
      return false;
    }
    if (!(radius > 0.0)) return false;
    const double xx = dot(x_axis, x_axis);
    const double yy = dot(y_axis, y_axis);
    const double xy = dot(x_axis, y_axis);
    return std::fabs(1.0 - xx) <= detail::kAnalyticUnitTolerance &&
           std::fabs(1.0 - yy) <= detail::kAnalyticUnitTolerance &&
           std::fabs(xy) <= detail::kAnalyticUnitTolerance;
  }

  Point3d eval(double angle) const noexcept {
    const double c = std::cos(angle);
    const double s = std::sin(angle);
    return Point3d{center[0] + radius * (c * x_axis[0] + s * y_axis[0]),
                   center[1] + radius * (c * x_axis[1] + s * y_axis[1]),
                   center[2] + radius * (c * x_axis[2] + s * y_axis[2])};
  }

  Vec3d derivative(double angle) const noexcept {
    const double c = std::cos(angle);
    const double s = std::sin(angle);
    return Vec3d{radius * (-s * x_axis[0] + c * y_axis[0]),
                 radius * (-s * x_axis[1] + c * y_axis[1]),
                 radius * (-s * x_axis[2] + c * y_axis[2])};
  }

  // Exact bounds of the arc [a0, a1], not a sampled approximation.
  //
  // Along axis i the circle is center[i] + radius * A * cos(t - phi) with
  // A = hypot(x_axis[i], y_axis[i]) and phi = atan2(y_axis[i], x_axis[i]), so
  // its extrema on the arc are attained at a0, at a1, or at a stationary angle
  // phi + k*pi inside the arc. Candidates are enumerated around
  // k0 = floor((a0 - phi)/pi): four consecutive stationary angles span 3*pi, so
  // k0-1 .. k0+2 contain every stationary angle in [a0, a1] for an arc shorter
  // than 2*pi, and for a longer arc they still contain both a maximum and a
  // minimum. A candidate outside the arc is clamped to the nearest endpoint,
  // where it re-evaluates a value that is already a bound of the arc.
  //
  // An axis with x_axis[i] == y_axis[i] == 0 is constant, which the endpoint
  // evaluations already cover. For very large angles the candidate angles are
  // reconstructed at the double spacing of a0, but each candidate is a
  // stationary point, where the circle value varies only quadratically with
  // the angle, so the bounds stay close to the true extrema.
  //
  // Precondition: a0, a1 finite and a0 <= a1.
  Box3d bounds(double a0, double a1) const noexcept {
    assert(std::isfinite(a0) && std::isfinite(a1) && a0 <= a1);

    constexpr double kPi = 3.14159265358979323846;
    Box3d box = Box3d::empty();
    auto expand_at = [&](double t) noexcept {
      double tc = t;
      if (tc < a0) tc = a0;
      if (tc > a1) tc = a1;
      box.expand(eval(tc));
    };

    expand_at(a0);
    expand_at(a1);
    for (int i = 0; i < 3; ++i) {
      const double x = x_axis[i];
      const double y = y_axis[i];
      if (x == 0.0 && y == 0.0) continue;
      const double base = std::atan2(y, x);
      const double k0 = std::floor((a0 - base) / kPi);
      for (int off = -1; off <= 2; ++off) {
        expand_at(base + (k0 + static_cast<double>(off)) * kPi);
      }
    }
    return box;
  }
};

}  // namespace v1
}  // namespace geom
