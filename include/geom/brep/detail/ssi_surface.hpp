#pragma once

#include <algorithm>
#include <cmath>
#include <numbers>
#include <type_traits>
#include <variant>

#include <geom/brep/plane_frame.hpp>
#include <geom/brep/topology.hpp>
#include <geom/core/error.hpp>
#include <geom/core/result.hpp>
#include <geom/surfaces/analytic.hpp>
#include <geom/surfaces/nurbs_surface.hpp>
#include <geom/surfaces/plane.hpp>

namespace geom {
inline namespace v1 {

// One sample of a surface parameterization: the point together with its two
// first partial derivatives. SSI consumes the triple everywhere it needs a
// local frame, so evaluation returns all three at once instead of asking the
// surface to repeat the work.
struct SurfacePoint {
  Point3d p{};
  Vec3d du{};
  Vec3d dv{};
};

// True when every component of the point and both partials is finite.
inline bool surface_point_finite(const SurfacePoint& point) noexcept {
  for (int i = 0; i < 3; ++i) {
    if (!std::isfinite(point.p[i]) || !std::isfinite(point.du[i]) ||
        !std::isfinite(point.dv[i])) {
      return false;
    }
  }
  return true;
}

namespace detail {

// Full turn, the period of every periodic surface parameter in this library.
inline constexpr double kSsiTwoPi = 2.0 * std::numbers::pi_v<double>;

// Folds an angle into [-pi, pi). std::remainder already lands in [-pi, pi];
// the extra branch moves the positive boundary onto the negative one so the
// seam angle has a single representative in the half-open interval.
inline double ssi_wrap_angle(double value) noexcept {
  const double folded = std::remainder(value, kSsiTwoPi);
  return folded == std::numbers::pi_v<double> ? -std::numbers::pi_v<double> : folded;
}

// Validates a surface without evaluating it, using the same rejection rules as
// surface_eval: a plane is judged by make_plane_frame so its Errc passes
// through unchanged, every other surface by its own is_valid().
inline Result<void> ssi_validate_surface(const SurfaceGeometry& surface) {
  return std::visit(
      [](const auto& s) -> Result<void> {
        using Surface = std::decay_t<decltype(s)>;
        if constexpr (std::is_same_v<Surface, Plane3>) {
          const Result<PlaneFrame> frame = make_plane_frame(s);
          if (!frame.has_value()) return frame.error();
          return Result<void>{};
        } else {
          if (!s.is_valid()) {
            return Error{Errc::degenerate, "ssi: surface is not valid"};
          }
          return Result<void>{};
        }
      },
      surface);
}

// A plane has no native (u, v); its patch coordinates are the plane frame's
// orthonormal basis, with u along frame.u and v along frame.v.
inline Result<SurfacePoint> ssi_eval_plane(const Plane3& plane, double u, double v) {
  const Result<PlaneFrame> frame = make_plane_frame(plane);
  if (!frame.has_value()) return frame.error();
  const PlaneFrame& f = frame.value();
  SurfacePoint point;
  point.p = f.origin + f.u * u + f.v * v;
  point.du = f.u;
  point.dv = f.v;
  return point;
}

// Analytic and NURBS surfaces expose eval/du/dv directly. Those members are
// only defined for a valid surface, so the validity check must come first.
template <class Surface>
inline Result<SurfacePoint> ssi_eval_parametric(const Surface& surface, double u, double v) {
  if (!surface.is_valid()) {
    return Error{Errc::degenerate, "ssi: surface is not valid"};
  }
  SurfacePoint point;
  point.p = surface.eval(u, v);
  point.du = surface.du(u, v);
  point.dv = surface.dv(u, v);
  return point;
}

// Unit normal from the partials, for surfaces whose cross product is not
// available as a closed form. Scaling by the largest component first keeps the
// intermediate in [1, sqrt(3)]: cross products of near-degenerate partials can
// underflow, and the scaled hypot cannot overflow either.
inline Result<Vec3d> ssi_normal_from_derivatives(const Vec3d& du, const Vec3d& dv) {
  const Vec3d n = cross(du, dv);
  if (!std::isfinite(n[0]) || !std::isfinite(n[1]) || !std::isfinite(n[2])) {
    return Error{Errc::non_finite, "ssi: surface normal is not finite"};
  }
  const double scale = std::max({std::fabs(n[0]), std::fabs(n[1]), std::fabs(n[2])});
  if (scale == 0.0) {
    return Error{Errc::numeric_failure, "ssi: surface normal is degenerate"};
  }
  const Vec3d scaled{n[0] / scale, n[1] / scale, n[2] / scale};
  const double length = std::hypot(scaled[0], scaled[1], scaled[2]);
  return Vec3d{scaled[0] / length, scaled[1] / length, scaled[2] / length};
}

}  // namespace detail

// Evaluates position and first partials. A plane is evaluated through
// make_plane_frame, so its (u, v) parameterization is the orthonormal frame and
// a rejected plane propagates make_plane_frame's Errc; every other surface uses
// its own parameterization and is rejected as degenerate when invalid.
inline Result<SurfacePoint> surface_eval(const SurfaceGeometry& surface, double u, double v) {
  return std::visit(
      [u, v](const auto& s) -> Result<SurfacePoint> {
        using Surface = std::decay_t<decltype(s)>;
        if constexpr (std::is_same_v<Surface, Plane3>) {
          return detail::ssi_eval_plane(s, u, v);
        } else {
          return detail::ssi_eval_parametric(s, u, v);
        }
      },
      surface);
}

// Unit surface normal, oriented like cross(du, dv) of the parameterization.
// Planes and NURBS surfaces build it from the partials; the analytic surfaces
// use their closed forms, which stay unit at poles and apexes where the cross
// product of the partials vanishes.
inline Result<Vec3d> surface_normal(const SurfaceGeometry& surface, double u, double v) {
  return std::visit(
      [u, v](const auto& s) -> Result<Vec3d> {
        using Surface = std::decay_t<decltype(s)>;
        if constexpr (std::is_same_v<Surface, Plane3>) {
          const Result<PlaneFrame> frame = make_plane_frame(s);
          if (!frame.has_value()) return frame.error();
          return cross(frame.value().u, frame.value().v);
        } else if constexpr (std::is_same_v<Surface, NurbsSurface>) {
          if (!s.is_valid()) {
            return Error{Errc::degenerate, "ssi: surface is not valid"};
          }
          return detail::ssi_normal_from_derivatives(s.du(u, v), s.dv(u, v));
        } else {
          if (!s.is_valid()) {
            return Error{Errc::degenerate, "ssi: surface is not valid"};
          }
          return s.normal(u, v);
        }
      },
      surface);
}

// True when `axis` (0 = u, 1 = v) is periodic with 2*pi: the angular parameter
// of a cylinder, cone, sphere or torus. Planes and NURBS surfaces are bounded
// parameterizations, so neither axis is periodic for them.
inline bool surface_periodic(const SurfaceGeometry& surface, int axis) noexcept {
  if (axis != 0 && axis != 1) return false;
  return std::visit(
      [axis](const auto& s) noexcept -> bool {
        using Surface = std::decay_t<decltype(s)>;
        (void)s;
        if constexpr (std::is_same_v<Surface, Cylinder3> || std::is_same_v<Surface, Cone3> ||
                      std::is_same_v<Surface, Sphere3>) {
          return axis == 0;
        } else if constexpr (std::is_same_v<Surface, Torus3>) {
          return true;
        } else {
          return false;
        }
      },
      surface);
}

// Folds a periodic parameter into [-pi, pi) and returns a non-periodic one
// unchanged. A non-finite parameter is rejected before folding, because
// std::remainder would only turn it into NaN.
inline Result<double> wrap_surface_param(const SurfaceGeometry& surface, int axis, double value) {
  if (!std::isfinite(value)) {
    return Error{Errc::non_finite, "ssi: surface parameter is not finite"};
  }
  if (!surface_periodic(surface, axis)) return value;
  return detail::ssi_wrap_angle(value);
}

}  // namespace v1
}  // namespace geom
