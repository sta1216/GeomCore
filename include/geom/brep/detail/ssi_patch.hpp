#pragma once

#include <array>
#include <cmath>
#include <utility>
#include <variant>

#include <geom/brep/detail/ssi_surface.hpp>
#include <geom/brep/topology.hpp>
#include <geom/core/error.hpp>
#include <geom/core/result.hpp>

namespace geom {
inline namespace v1 {

// Axis-aligned rectangular domain on one surface, in that surface's own (u, v)
// parameterization. Bounds are stored exactly as supplied; periodic axes are
// folded at query time, so a patch may sit anywhere on the parameter circle.
struct SsiPatch {
  SurfaceGeometry surface;
  std::array<double, 2> u{};
  std::array<double, 2> v{};
};

namespace detail {

// Closed-interval test for one axis. A non-periodic axis is a plain comparison.
// A periodic axis is compared on the parameter circle: bounds and query are all
// folded into [-pi, pi), and if the folded bounds invert, the interval is the
// one crossing the seam. A span of one full period covers the whole axis.
inline bool ssi_patch_axis_contains(const SsiPatch& patch, int axis, double value) noexcept {
  if (!std::isfinite(value)) return false;
  const std::array<double, 2>& range = axis == 0 ? patch.u : patch.v;
  const double lo = range[0];
  const double hi = range[1];
  if (!surface_periodic(patch.surface, axis)) {
    return lo <= value && value <= hi;
  }
  if (hi - lo >= kSsiTwoPi) return true;
  const double folded_lo = ssi_wrap_angle(lo);
  const double folded_hi = ssi_wrap_angle(hi);
  const double folded = ssi_wrap_angle(value);
  if (folded_lo <= folded_hi) {
    return folded_lo <= folded && folded <= folded_hi;
  }
  return folded >= folded_lo || folded <= folded_hi;
}

// A NURBS surface has a bounded parameter domain and its eval() clamps
// out-of-domain parameters. A patch wider than that domain would silently
// denote clamped geometry instead of the parameter rectangle its bounds
// describe, so it is rejected as unsupported. Analytic surfaces have unbounded
// parameters (a cylinder's v has no domain) and planes have no native bounds,
// so only NURBS surfaces are checked.
inline Result<void> ssi_patch_domain_check(const SurfaceGeometry& surface,
                                           const std::array<double, 2>& u,
                                           const std::array<double, 2>& v) {
  const NurbsSurface* nurbs = std::get_if<NurbsSurface>(&surface);
  if (nurbs == nullptr) return Result<void>{};
  const std::pair<double, double> domain_u = nurbs->domain_u();
  const std::pair<double, double> domain_v = nurbs->domain_v();
  if (u[0] < domain_u.first || u[1] > domain_u.second || v[0] < domain_v.first ||
      v[1] > domain_v.second) {
    return Error{Errc::unsupported, "ssi: patch domain lies outside the NURBS surface domain"};
  }
  return Result<void>{};
}

}  // namespace detail

// Builds a patch after validating the surface and the domain. Non-finite
// bounds are rejected before the ordering test, because a NaN bound compares
// false against everything and would otherwise look merely inverted. An
// inverted or empty interval is degenerate. A periodic axis wider than 2*pi is
// unsupported: it would cover the same parameterization more than once.
// Non-periodic axes have no width limit. A NURBS patch must also lie inside the
// surface's own domain, because eval() would silently clamp it.
inline Result<SsiPatch> make_ssi_patch(const SurfaceGeometry& surface, std::array<double, 2> u,
                                       std::array<double, 2> v) {
  const Result<void> valid = detail::ssi_validate_surface(surface);
  if (!valid.has_value()) return valid.error();
  if (!std::isfinite(u[0]) || !std::isfinite(u[1]) || !std::isfinite(v[0]) ||
      !std::isfinite(v[1])) {
    return Error{Errc::non_finite, "ssi: patch domain is not finite"};
  }
  if (!(u[0] < u[1]) || !(v[0] < v[1])) {
    return Error{Errc::degenerate, "ssi: patch domain is inverted or empty"};
  }
  const Result<void> in_domain = detail::ssi_patch_domain_check(surface, u, v);
  if (!in_domain.has_value()) return in_domain.error();
  if (surface_periodic(surface, 0) && u[1] - u[0] > detail::kSsiTwoPi) {
    return Error{Errc::unsupported, "ssi: patch spans more than one period in u"};
  }
  if (surface_periodic(surface, 1) && v[1] - v[0] > detail::kSsiTwoPi) {
    return Error{Errc::unsupported, "ssi: patch spans more than one period in v"};
  }
  SsiPatch patch;
  patch.surface = surface;
  patch.u = u;
  patch.v = v;
  return patch;
}

// True when (u, v) lies in the closed patch domain, folding periodic axes.
inline bool patch_contains(const SsiPatch& patch, double u, double v) noexcept {
  return detail::ssi_patch_axis_contains(patch, 0, u) &&
         detail::ssi_patch_axis_contains(patch, 1, v);
}

}  // namespace v1
}  // namespace geom
