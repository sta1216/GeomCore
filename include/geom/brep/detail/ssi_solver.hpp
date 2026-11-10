#pragma once

#include <algorithm>
#include <array>
#include <cmath>
#include <geom/brep/detail/ssi_patch.hpp>
#include <geom/brep/detail/ssi_surface.hpp>
#include <geom/core/error.hpp>
#include <geom/core/result.hpp>
#include <limits>
#include <vector>

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

// One traced polyline segment of an intersection curve. a and b are the two
// refined endpoints, ua/ub are the A-side parameters, va/vb the B-side
// parameters, and residual is the gap-only endpoint residual
// max(ssi_residual(a), ssi_residual(b)).
struct SsiSegment {
  Point3d a{};
  Point3d b{};
  std::array<double, 2> ua{};
  std::array<double, 2> ub{};
  std::array<double, 2> va{};
  std::array<double, 2> vb{};
  double residual = 0.0;
};

// One traced intersection curve. closed is true only when the traced polyline
// could be welded back to its first vertex; the closing vertex is not stored as
// a duplicate and the final segment ends at the first vertex.
struct SsiCurve {
  std::vector<SsiSegment> segments;
  bool closed = false;
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

namespace detail {

// A closing candidate is accepted only when its tangent still points the same
// way as the start tangent. Without this a long step can pass near the start on
// the far side of a loop and close the curve early, silently dropping the arc
// that lies between the two visits.
inline constexpr double kSsiTraceClosureTangent = 0.9;

// Parameter difference on one axis, measured across the periodic seam when the
// surface wraps on that axis.
inline double ssi_trace_param_delta(const SurfaceGeometry& surface, int axis, double first,
                                    double second) noexcept {
  double delta = second - first;
  if (surface_periodic(surface, axis)) delta = ssi_wrap_angle(delta);
  return delta;
}

// Distance between two parameter pairs on one surface, folding periodic axes.
inline double ssi_trace_param_distance(const SurfaceGeometry& surface,
                                       const std::array<double, 2>& first,
                                       const std::array<double, 2>& second) noexcept {
  const double du = ssi_trace_param_delta(surface, 0, first[0], second[0]);
  const double dv = ssi_trace_param_delta(surface, 1, first[1], second[1]);
  return std::sqrt(du * du + dv * dv);
}

// Least-squares parameter velocity of unit surface speed along `tangent`:
// solves [du dv] * velocity = tangent for one surface. The intersection tangent
// lies in both tangent planes, so this is the parameter velocity that moves the
// surface point along the intersection.
inline bool ssi_trace_parameter_velocity(const SurfacePoint& point, const Vec3d& tangent,
                                         std::array<double, 2>& velocity) noexcept {
  const double uu = dot(point.du, point.du);
  const double uv = dot(point.du, point.dv);
  const double vv = dot(point.dv, point.dv);
  const double ut = dot(point.du, tangent);
  const double vt = dot(point.dv, tangent);
  const double determinant = uu * vv - uv * uv;
  if (!std::isfinite(determinant) || !(determinant > 0.0)) return false;
  const double inverse = 1.0 / determinant;
  velocity[0] = (vv * ut - uv * vt) * inverse;
  velocity[1] = (uu * vt - uv * ut) * inverse;
  return std::isfinite(velocity[0]) && std::isfinite(velocity[1]);
}

// Unit intersection tangent plus the parameter velocity on each surface at one
// refined point. False when the normals are parallel (tangent surfaces) or a
// surface derivative is degenerate, in which case the trace cannot continue.
struct SsiTraceFrame {
  Vec3d tangent{};
  std::array<double, 2> velocity_a{};
  std::array<double, 2> velocity_b{};
};

inline bool ssi_trace_frame(const SsiPatch& a, const SsiPatch& b, const SsiPoint& point,
                            SsiTraceFrame& frame) noexcept {
  const Result<SurfacePoint> on_a = surface_eval(a.surface, point.ua[0], point.ua[1]);
  const Result<SurfacePoint> on_b = surface_eval(b.surface, point.va[0], point.va[1]);
  if (!on_a.has_value() || !on_b.has_value()) return false;
  if (!surface_point_finite(on_a.value()) || !surface_point_finite(on_b.value())) return false;
  const Vec3d normal_a = cross(on_a.value().du, on_a.value().dv);
  const Vec3d normal_b = cross(on_b.value().du, on_b.value().dv);
  const Vec3d candidate = cross(normal_a, normal_b);
  const double length = norm2(candidate);
  if (!std::isfinite(length) || !(length > 0.0)) return false;
  frame.tangent = candidate * (1.0 / length);
  return ssi_trace_parameter_velocity(on_a.value(), frame.tangent, frame.velocity_a) &&
         ssi_trace_parameter_velocity(on_b.value(), frame.tangent, frame.velocity_b);
}

// Why one directional walk stopped. left_patch and refinement_failed are normal
// open-curve endings; step_budget means the walk neither closed nor left the
// patch before the shared budget ran out.
enum class SsiTraceEnd { closed, left_patch, refinement_failed, tangent_failed, step_budget };

struct SsiTraceHalf {
  std::vector<SsiPoint> points;
  SsiTraceEnd end = SsiTraceEnd::left_patch;
};

// Walks one direction from `start`. Every accepted vertex is the result of
// refine_ssi_point, and every candidate is checked against both patch domains
// and against the tolerance: the gap residual alone does not certify that the
// parameters stayed inside a patch. A candidate that leaves either patch ends
// the walk at the last accepted vertex. The step budget is shared by both
// directions of one curve.
inline SsiTraceHalf ssi_trace_half(const SsiPatch& a, const SsiPatch& b, const SsiPoint& start,
                                   const SsiTraceFrame& start_frame, double sign, double step_start,
                                   double step_max, double tolerance, int max_steps,
                                   int& steps_used) {
  SsiTraceHalf half;
  half.points.push_back(start);
  SsiPoint current = start;
  double step = step_start;
  for (;;) {
    if (steps_used >= max_steps) {
      half.end = SsiTraceEnd::step_budget;
      return half;
    }
    SsiTraceFrame frame;
    if (!ssi_trace_frame(a, b, current, frame)) {
      half.end = SsiTraceEnd::tangent_failed;
      return half;
    }
    std::array<double, 2> predicted_ua{current.ua[0] + sign * step * frame.velocity_a[0],
                                       current.ua[1] + sign * step * frame.velocity_a[1]};
    std::array<double, 2> predicted_va{current.va[0] + sign * step * frame.velocity_b[0],
                                       current.va[1] + sign * step * frame.velocity_b[1]};
    // refine_ssi_point checks patch membership before folding periodic
    // parameters, so predictions must be folded into the patch interval first.
    if (!ssi_solver_wrap(a.surface, predicted_ua) || !ssi_solver_wrap(b.surface, predicted_va)) {
      half.end = SsiTraceEnd::refinement_failed;
      return half;
    }
    if (!patch_contains(a, predicted_ua[0], predicted_ua[1]) ||
        !patch_contains(b, predicted_va[0], predicted_va[1])) {
      half.end = SsiTraceEnd::left_patch;
      return half;
    }
    const Result<SsiPoint> refined = refine_ssi_point(
        a, b, {predicted_ua[0], predicted_ua[1], predicted_va[0], predicted_va[1]}, tolerance);
    if (!refined.has_value()) {
      half.end = SsiTraceEnd::refinement_failed;
      return half;
    }
    const SsiPoint next = refined.value();
    if (!patch_contains(a, next.ua[0], next.ua[1]) || !patch_contains(b, next.va[0], next.va[1]) ||
        !(next.residual <= tolerance)) {
      half.end = SsiTraceEnd::refinement_failed;
      return half;
    }
    // A return to the start is recognised from the parameters themselves, not
    // from a raw 3D distance: two nearby parameters can sit on opposite sides of
    // a loop. Three accepted steps must precede the candidate so a one- or
    // two-step walk cannot weld itself shut.
    const double closure_limit = std::max(2.0 * step, 8.0 * tolerance);
    const double closure_a = ssi_trace_param_distance(a.surface, next.ua, start.ua);
    const double closure_b = ssi_trace_param_distance(b.surface, next.va, start.va);
    const double closure_3d = norm2(next.p - start.p);
    if (half.points.size() >= 4 && closure_a <= closure_limit && closure_b <= closure_limit &&
        closure_3d <= closure_limit) {
      SsiTraceFrame candidate;
      if (ssi_trace_frame(a, b, next, candidate) &&
          dot(candidate.tangent, start_frame.tangent) >= kSsiTraceClosureTangent) {
        // The candidate is a re-visit of the start; stop before adding it so the
        // caller can close the curve with the start vertex itself.
        half.end = SsiTraceEnd::closed;
        return half;
      }
    }
    half.points.push_back(next);
    current = next;
    step = std::min(step * 1.5, step_max);
    ++steps_used;
  }
}

}  // namespace detail

// Traces the intersection curves of two validated SSI patches from refined
// seeds. Each curve is walked in both null-space directions; a walk that
// returns to its start closes the curve, and the closing vertex is the start
// itself rather than a duplicate. The result is ordered by each curve's first
// vertex so it does not depend on seed order, and duplicated seeds collapse to
// one curve. A walk that crosses a patch boundary ends the curve there; a walk
// that exhausts the shared step budget without closing or leaving a patch is a
// numeric_failure rather than a silently truncated curve.
inline Result<std::vector<SsiCurve>> trace_ssi_curves(const SsiPatch& a, const SsiPatch& b,
                                                      const std::vector<SsiPoint>& seeds,
                                                      double tolerance, int max_steps = 4096) {
  if (!std::isfinite(tolerance) || !(tolerance > 0.0)) {
    return Error{Errc::tolerance_violation,
                 "ssi: trace tolerance must be finite and greater than zero"};
  }
  if (max_steps <= 0) {
    return Error{Errc::numeric_failure, "ssi: trace step budget must be positive"};
  }
  const Result<SsiPatch> validated_a = make_ssi_patch(a.surface, a.u, a.v);
  if (!validated_a.has_value()) return validated_a.error();
  const Result<SsiPatch> validated_b = make_ssi_patch(b.surface, b.u, b.v);
  if (!validated_b.has_value()) return validated_b.error();
  const SsiPatch& patch_a = validated_a.value();
  const SsiPatch& patch_b = validated_b.value();

  Box3d bounds = patch_bounds(patch_a, tolerance);
  if (bounds.is_empty()) {
    return Error{Errc::numeric_failure, "ssi: trace cannot bound patch A"};
  }
  const Box3d other = patch_bounds(patch_b, tolerance);
  if (other.is_empty()) {
    return Error{Errc::numeric_failure, "ssi: trace cannot bound patch B"};
  }
  bounds.expand(other);
  const double diagonal = norm2(bounds.hi - bounds.lo);
  if (!std::isfinite(diagonal) || !(diagonal > 0.0)) {
    return Error{Errc::numeric_failure, "ssi: trace patch bounds are degenerate"};
  }
  const double step_start = std::min(diagonal / 64.0, 0.25);
  const double step_max = diagonal / 8.0;
  if (!std::isfinite(step_start) || !(step_start > 0.0) || !std::isfinite(step_max)) {
    return Error{Errc::numeric_failure, "ssi: trace step size is not finite"};
  }

  // Refine every seed first, then order and deduplicate the refined starts.
  // Ordering by the refined parameters, not by input index, makes the surviving
  // representative of a duplicate cluster independent of seed order. A seed
  // that cannot be refined is skipped: the curve it would have started is a
  // question for the caller, not a reason to drop the other curves.
  std::vector<SsiPoint> origins;
  origins.reserve(seeds.size());
  for (const SsiPoint& seed : seeds) {
    std::array<double, 2> ua = seed.ua;
    std::array<double, 2> va = seed.va;
    if (!detail::ssi_solver_wrap(patch_a.surface, ua) ||
        !detail::ssi_solver_wrap(patch_b.surface, va)) {
      continue;
    }
    const Result<SsiPoint> refined =
        refine_ssi_point(patch_a, patch_b, {ua[0], ua[1], va[0], va[1]}, tolerance);
    if (!refined.has_value()) continue;
    if (!(refined.value().residual <= tolerance)) continue;
    origins.push_back(refined.value());
  }
  std::stable_sort(origins.begin(), origins.end(), [](const SsiPoint& left, const SsiPoint& right) {
    if (left.ua != right.ua) return left.ua < right.ua;
    return left.va < right.va;
  });
  std::vector<SsiPoint> unique_origins;
  unique_origins.reserve(origins.size());
  for (const SsiPoint& origin : origins) {
    bool duplicate = false;
    for (const SsiPoint& kept : unique_origins) {
      const double distance_a =
          detail::ssi_trace_param_distance(patch_a.surface, origin.ua, kept.ua);
      const double distance_b =
          detail::ssi_trace_param_distance(patch_b.surface, origin.va, kept.va);
      if (distance_a <= 8.0 * tolerance && distance_b <= 8.0 * tolerance) {
        duplicate = true;
        break;
      }
    }
    if (!duplicate) unique_origins.push_back(origin);
  }

  std::vector<SsiCurve> curves;
  curves.reserve(unique_origins.size());
  for (const SsiPoint& origin : unique_origins) {
    detail::SsiTraceFrame start_frame;
    if (!detail::ssi_trace_frame(patch_a, patch_b, origin, start_frame)) continue;
    int steps_used = 0;
    detail::SsiTraceHalf forward =
        detail::ssi_trace_half(patch_a, patch_b, origin, start_frame, 1.0, step_start, step_max,
                               tolerance, max_steps, steps_used);
    if (forward.end == detail::SsiTraceEnd::step_budget) {
      return Error{Errc::numeric_failure, "trace_ssi_curves: step budget exhausted"};
    }
    std::vector<SsiPoint> points;
    bool closed = false;
    if (forward.end == detail::SsiTraceEnd::closed) {
      points = forward.points;
      closed = true;
    } else {
      detail::SsiTraceHalf backward =
          detail::ssi_trace_half(patch_a, patch_b, origin, start_frame, -1.0, step_start, step_max,
                                 tolerance, max_steps, steps_used);
      if (backward.end == detail::SsiTraceEnd::step_budget) {
        return Error{Errc::numeric_failure, "trace_ssi_curves: step budget exhausted"};
      }
      points.reserve(backward.points.size() + forward.points.size() - 1);
      for (auto it = backward.points.rbegin(); it != backward.points.rend(); ++it) {
        points.push_back(*it);
      }
      for (std::size_t index = 1; index < forward.points.size(); ++index) {
        points.push_back(forward.points[index]);
      }
    }
    if (points.size() < 2) continue;
    SsiCurve curve;
    curve.closed = closed;
    const std::size_t segment_count = closed ? points.size() : points.size() - 1;
    curve.segments.reserve(segment_count);
    for (std::size_t index = 0; index < segment_count; ++index) {
      const SsiPoint& first = points[index];
      const SsiPoint& second = points[(index + 1) % points.size()];
      SsiSegment segment;
      segment.a = first.p;
      segment.b = second.p;
      segment.ua = first.ua;
      segment.ub = second.ua;
      segment.va = first.va;
      segment.vb = second.va;
      segment.residual =
          std::max(ssi_residual(patch_a, patch_b, first), ssi_residual(patch_a, patch_b, second));
      curve.segments.push_back(segment);
    }
    curves.push_back(curve);
  }

  std::stable_sort(curves.begin(), curves.end(), [](const SsiCurve& left, const SsiCurve& right) {
    const SsiSegment& one = left.segments.front();
    const SsiSegment& two = right.segments.front();
    if (one.ua != two.ua) return one.ua < two.ua;
    if (one.va != two.va) return one.va < two.va;
    for (int axis = 0; axis < 3; ++axis) {
      if (one.a[axis] != two.a[axis]) return one.a[axis] < two.a[axis];
    }
    return left.segments.size() < right.segments.size();
  });
  return curves;
}

}  // namespace v1
}  // namespace geom
