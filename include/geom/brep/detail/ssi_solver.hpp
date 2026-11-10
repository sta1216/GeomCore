#pragma once

#include <algorithm>
#include <array>
#include <cmath>
#include <geom/brep/detail/ssi_patch.hpp>
#include <geom/brep/detail/ssi_surface.hpp>
#include <geom/core/error.hpp>
#include <geom/core/result.hpp>
#include <limits>

namespace geom {
inline namespace v1 {

// One verified surface-surface intersection point. The two parameter pairs are
// the refined seeds, p is the point on surface A at (ua[0], ua[1]), and residual
// is ssi_residual for those parameters. Periodic parameters are stored folded
// into [-pi, pi), so ua/va always name the same parameter value the patch
// domain uses.
struct SsiPoint {
  std::array<double, 2> ua{};
  std::array<double, 2> va{};
  Point3d p{};
  double residual = 0.0;
};

namespace detail {

// A 3x3 system whose pivot falls below this fraction of its column's original
// magnitude is treated as singular: at that point the Newton direction is
// dominated by rounding noise and the step would be meaningless.
inline constexpr double kSsiSolverPivotThreshold = 1e-13;

// The step is halved at most this many times while looking for a decrease of
// |F|; 2^-40 shrinks a unit step to below 1e-12, far past any useful motion.
inline constexpr int kSsiSolverMaxHalvings = 40;

// Folds the periodic axes of one parameter pair in place. Returns false only
// for a non-finite value, which callers treat as a contract violation.
inline bool ssi_solver_wrap(const SurfaceGeometry& surface,
                            std::array<double, 2>& params) noexcept {
  for (int axis = 0; axis < 2; ++axis) {
    if (!surface_periodic(surface, axis)) continue;
    if (!std::isfinite(params[static_cast<std::size_t>(axis)])) return false;
    params[static_cast<std::size_t>(axis)] = ssi_wrap_angle(params[static_cast<std::size_t>(axis)]);
  }
  return true;
}

inline double ssi_solver_det3(const std::array<std::array<double, 3>, 3>& m) noexcept {
  return m[0][0] * (m[1][1] * m[2][2] - m[1][2] * m[2][1]) -
         m[0][1] * (m[1][0] * m[2][2] - m[1][2] * m[2][0]) +
         m[0][2] * (m[1][0] * m[2][1] - m[1][1] * m[2][0]);
}

// Gaussian elimination with partial pivoting. Each pivot is compared against
// the largest entry the column had before elimination, so the test is
// scale-invariant. False means the system is numerically singular (or the
// solution overflowed); the caller reports numeric_failure and never uses a
// half-computed direction.
inline bool ssi_solver_solve3(std::array<std::array<double, 3>, 3> matrix, const Vec3d& rhs,
                              std::array<double, 3>& solution) noexcept {
  double column_max[3] = {0.0, 0.0, 0.0};
  for (int column = 0; column < 3; ++column) {
    for (int row = 0; row < 3; ++row) {
      column_max[column] = std::max(
          column_max[column],
          std::fabs(matrix[static_cast<std::size_t>(row)][static_cast<std::size_t>(column)]));
    }
    if (!(column_max[column] > 0.0)) return false;
  }
  double scratch[3] = {rhs[0], rhs[1], rhs[2]};
  for (int column = 0; column < 3; ++column) {
    int pivot_row = column;
    double pivot =
        std::fabs(matrix[static_cast<std::size_t>(column)][static_cast<std::size_t>(column)]);
    for (int row = column + 1; row < 3; ++row) {
      const double candidate =
          std::fabs(matrix[static_cast<std::size_t>(row)][static_cast<std::size_t>(column)]);
      if (candidate > pivot) {
        pivot = candidate;
        pivot_row = row;
      }
    }
    if (!std::isfinite(pivot) || pivot <= kSsiSolverPivotThreshold * column_max[column]) {
      return false;
    }
    if (pivot_row != column) {
      std::swap(matrix[static_cast<std::size_t>(column)],
                matrix[static_cast<std::size_t>(pivot_row)]);
      std::swap(scratch[column], scratch[pivot_row]);
    }
    const double diagonal =
        matrix[static_cast<std::size_t>(column)][static_cast<std::size_t>(column)];
    for (int row = column + 1; row < 3; ++row) {
      const double factor =
          matrix[static_cast<std::size_t>(row)][static_cast<std::size_t>(column)] / diagonal;
      if (!std::isfinite(factor)) return false;
      for (int k = column; k < 3; ++k) {
        matrix[static_cast<std::size_t>(row)][static_cast<std::size_t>(k)] -=
            factor * matrix[static_cast<std::size_t>(column)][static_cast<std::size_t>(k)];
      }
      scratch[row] -= factor * scratch[column];
    }
  }
  for (int row = 2; row >= 0; --row) {
    double accumulated = scratch[row];
    for (int k = row + 1; k < 3; ++k) {
      accumulated -= matrix[static_cast<std::size_t>(row)][static_cast<std::size_t>(k)] *
                     solution[static_cast<std::size_t>(k)];
    }
    const double value =
        accumulated / matrix[static_cast<std::size_t>(row)][static_cast<std::size_t>(row)];
    if (!std::isfinite(value)) return false;
    solution[static_cast<std::size_t>(row)] = value;
  }
  return true;
}

}  // namespace detail

// Residual of a candidate intersection point: the distance between the two
// surface points at the stored parameter pairs, recomputed here so a stale or
// forged p cannot mask a parameter mismatch. This is only the parameter gap.
// It does not certify that the two parameter pairs correspond to the same
// geometric point on a self-intersecting patch, and it does not check patch
// membership: seams and poles can map distinct pairs onto one point, and a
// candidate outside a patch can still show a small gap. Callers that need
// those guarantees must check them separately (patch_contains), exactly as
// refine_ssi_point does. Non-finite input or a failed evaluation yields +inf.
inline double ssi_residual(const SsiPatch& a, const SsiPatch& b, const SsiPoint& point) noexcept {
  std::array<double, 2> ua{point.ua[0], point.ua[1]};
  std::array<double, 2> va{point.va[0], point.va[1]};
  if (!detail::ssi_solver_wrap(a.surface, ua) || !detail::ssi_solver_wrap(b.surface, va)) {
    return std::numeric_limits<double>::infinity();
  }
  const Result<SurfacePoint> on_a = surface_eval(a.surface, ua[0], ua[1]);
  const Result<SurfacePoint> on_b = surface_eval(b.surface, va[0], va[1]);
  if (!on_a.has_value() || !on_b.has_value()) return std::numeric_limits<double>::infinity();
  if (!surface_point_finite(on_a.value()) || !surface_point_finite(on_b.value())) {
    return std::numeric_limits<double>::infinity();
  }
  const double gap = norm2(on_a.value().p - on_b.value().p);
  if (!std::isfinite(gap)) return std::numeric_limits<double>::infinity();
  return gap;
}

// Refines a four-parameter seed on the intersection of two patches with the
// literal 3x4 Newton scheme of the plan: the step lies in the null space of
// F = SA(u, v) - SB(s, t), obtained from the 3x3 minors of the Jacobian, with
// the largest component of the null vector held fixed. The remaining 3x3
// system is solved with partial pivoting; the full step is then halved until
// |F| decreases. Periodic axes are folded every step, a seed outside either
// patch is rejected instead of clamped, and an iterate that cannot be kept
// inside both patches is a numeric_failure.
inline Result<SsiPoint> refine_ssi_point(const SsiPatch& a, const SsiPatch& b,
                                         std::array<double, 4> seed, double tolerance,
                                         int max_iterations = 24) {
  if (!std::isfinite(tolerance) || !(tolerance > 0.0)) {
    return Error{Errc::tolerance_violation,
                 "ssi: Newton tolerance must be finite and greater than zero"};
  }
  if (max_iterations <= 0) {
    return Error{Errc::numeric_failure, "ssi: Newton iteration budget must be positive"};
  }
  for (int index = 0; index < 4; ++index) {
    if (!std::isfinite(seed[static_cast<std::size_t>(index)])) {
      return Error{Errc::non_finite, "ssi: Newton seed is not finite"};
    }
  }
  const Result<SsiPatch> valid_a = make_ssi_patch(a.surface, a.u, a.v);
  if (!valid_a.has_value()) return valid_a.error();
  const Result<SsiPatch> valid_b = make_ssi_patch(b.surface, b.u, b.v);
  if (!valid_b.has_value()) return valid_b.error();
  const SsiPatch& patch_a = valid_a.value();
  const SsiPatch& patch_b = valid_b.value();
  if (!patch_contains(patch_a, seed[0], seed[1]) || !patch_contains(patch_b, seed[2], seed[3])) {
    return Error{Errc::numeric_failure, "ssi: Newton seed lies outside a patch domain"};
  }
  std::array<double, 2> ua{seed[0], seed[1]};
  std::array<double, 2> va{seed[2], seed[3]};
  if (!detail::ssi_solver_wrap(patch_a.surface, ua) ||
      !detail::ssi_solver_wrap(patch_b.surface, va)) {
    return Error{Errc::numeric_failure, "ssi: Newton seed cannot be wrapped"};
  }
  const Result<SurfacePoint> seed_a = surface_eval(patch_a.surface, ua[0], ua[1]);
  if (!seed_a.has_value()) return seed_a.error();
  const Result<SurfacePoint> seed_b = surface_eval(patch_b.surface, va[0], va[1]);
  if (!seed_b.has_value()) return seed_b.error();
  if (!surface_point_finite(seed_a.value()) || !surface_point_finite(seed_b.value())) {
    return Error{Errc::non_finite, "ssi: surface evaluation is not finite"};
  }
  SurfacePoint current_a = seed_a.value();
  SurfacePoint current_b = seed_b.value();
  double gap = norm2(current_a.p - current_b.p);
  if (!std::isfinite(gap)) {
    return Error{Errc::numeric_failure, "ssi: initial separation is not finite"};
  }

  for (int iteration = 0;; ++iteration) {
    if (gap <= tolerance) break;
    if (iteration >= max_iterations) {
      return Error{Errc::numeric_failure, "ssi: Newton iteration did not converge"};
    }
    const Vec3d mismatch = current_a.p - current_b.p;
    std::array<Vec3d, 4> columns{current_a.du, current_a.dv,
                                 Vec3d{-current_b.du[0], -current_b.du[1], -current_b.du[2]},
                                 Vec3d{-current_b.dv[0], -current_b.dv[1], -current_b.dv[2]}};
    double scale = 0.0;
    for (const Vec3d& column : columns) {
      for (int component = 0; component < 3; ++component) {
        scale = std::max(scale, std::fabs(column[component]));
      }
    }
    if (!(scale > 0.0) || !std::isfinite(scale)) {
      return Error{Errc::numeric_failure, "ssi: Jacobian has no usable column"};
    }
    const double inverse_scale = 1.0 / scale;
    for (Vec3d& column : columns) {
      for (int component = 0; component < 3; ++component) column[component] *= inverse_scale;
    }

    // Null direction from the 3x3 minors, scaled by the same factor so the
    // determinants stay O(1) and cannot overflow.
    double null[4] = {0.0, 0.0, 0.0, 0.0};
    for (int drop = 0; drop < 4; ++drop) {
      std::array<std::array<double, 3>, 3> minor{};
      int column_index = 0;
      for (int column = 0; column < 4; ++column) {
        if (column == drop) continue;
        for (int row = 0; row < 3; ++row) {
          minor[static_cast<std::size_t>(row)][static_cast<std::size_t>(column_index)] =
              columns[static_cast<std::size_t>(column)][row];
        }
        ++column_index;
      }
      null[drop] = (drop % 2 == 0 ? 1.0 : -1.0) * detail::ssi_solver_det3(minor);
    }
    int fixed = 0;
    double largest = 0.0;
    for (int index = 0; index < 4; ++index) {
      const double magnitude = std::fabs(null[index]);
      if (!std::isfinite(magnitude)) {
        return Error{Errc::numeric_failure, "ssi: null direction is not finite"};
      }
      if (magnitude > largest) {
        largest = magnitude;
        fixed = index;
      }
    }
    if (!(largest > 0.0)) {
      return Error{Errc::numeric_failure, "ssi: intersection is numerically rank deficient"};
    }

    std::array<std::array<double, 3>, 3> system{};
    std::array<int, 3> variable{};
    for (int slot = 0; slot < 3; ++slot) {
      variable[static_cast<std::size_t>(slot)] = slot < fixed ? slot : slot + 1;
      for (int row = 0; row < 3; ++row) {
        system[static_cast<std::size_t>(row)][static_cast<std::size_t>(slot)] =
            columns[static_cast<std::size_t>(variable[static_cast<std::size_t>(slot)])][row];
      }
    }
    const Vec3d rhs{-mismatch[0], -mismatch[1], -mismatch[2]};
    std::array<double, 3> solution{};
    if (!detail::ssi_solver_solve3(system, rhs, solution)) {
      return Error{Errc::numeric_failure, "ssi: Newton system is numerically singular"};
    }
    std::array<double, 4> delta{};
    for (int slot = 0; slot < 3; ++slot) {
      // The columns were divided by `scale` above, so the solution of the
      // scaled system describes the same direction at 1/scale of its true
      // length: undo that division before applying the step.
      delta[static_cast<std::size_t>(variable[static_cast<std::size_t>(slot)])] =
          solution[static_cast<std::size_t>(slot)] * inverse_scale;
    }

    bool accepted = false;
    std::array<double, 2> accepted_ua{};
    std::array<double, 2> accepted_va{};
    SurfacePoint accepted_a{};
    SurfacePoint accepted_b{};
    double accepted_gap = 0.0;
    double alpha = 1.0;
    for (int halving = 0; halving <= detail::kSsiSolverMaxHalvings; ++halving, alpha *= 0.5) {
      std::array<double, 2> trial_ua{ua[0] + alpha * delta[0], ua[1] + alpha * delta[1]};
      std::array<double, 2> trial_va{va[0] + alpha * delta[2], va[1] + alpha * delta[3]};
      bool usable = std::isfinite(trial_ua[0]) && std::isfinite(trial_ua[1]) &&
                    std::isfinite(trial_va[0]) && std::isfinite(trial_va[1]);
      if (usable) usable = detail::ssi_solver_wrap(patch_a.surface, trial_ua);
      if (usable) usable = detail::ssi_solver_wrap(patch_b.surface, trial_va);
      if (usable) usable = patch_contains(patch_a, trial_ua[0], trial_ua[1]);
      if (usable) usable = patch_contains(patch_b, trial_va[0], trial_va[1]);
      if (!usable) continue;
      const Result<SurfacePoint> trial_a = surface_eval(patch_a.surface, trial_ua[0], trial_ua[1]);
      if (!trial_a.has_value()) return trial_a.error();
      const Result<SurfacePoint> trial_b = surface_eval(patch_b.surface, trial_va[0], trial_va[1]);
      if (!trial_b.has_value()) return trial_b.error();
      if (!surface_point_finite(trial_a.value()) || !surface_point_finite(trial_b.value())) {
        return Error{Errc::non_finite, "ssi: surface evaluation is not finite"};
      }
      const double trial_gap = norm2(trial_a.value().p - trial_b.value().p);
      if (!std::isfinite(trial_gap) || !(trial_gap < gap)) continue;
      accepted = true;
      accepted_ua = trial_ua;
      accepted_va = trial_va;
      accepted_a = trial_a.value();
      accepted_b = trial_b.value();
      accepted_gap = trial_gap;
      break;
    }
    if (!accepted) {
      return Error{Errc::numeric_failure,
                   "ssi: Newton step cannot be kept inside both patch domains"};
    }
    ua = accepted_ua;
    va = accepted_va;
    current_a = accepted_a;
    current_b = accepted_b;
    gap = accepted_gap;
  }

  SsiPoint result;
  result.ua = ua;
  result.va = va;
  result.p = current_a.p;
  result.residual = ssi_residual(patch_a, patch_b, result);
  return result;
}

}  // namespace v1
}  // namespace geom
