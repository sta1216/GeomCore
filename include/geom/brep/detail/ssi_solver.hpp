#pragma once

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <geom/brep/detail/ssi_patch.hpp>
#include <geom/brep/detail/ssi_surface.hpp>
#include <geom/core/error.hpp>
#include <geom/core/result.hpp>
#include <limits>
#include <type_traits>
#include <variant>
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

// A step whose tangent turns by more than this angle is too coarse for the
// local curvature: the step is halved and retried. This keeps small or strongly
// curved loops from winding past themselves.
inline constexpr double kSsiTraceCurvatureTurn = 0.5;

// No parameter may move by more than this fraction of its patch span in one
// step. The step is a 3D length, so the cap converts the patch geometry into a
// per-step surface length and keeps the trace independent of the model scale.
inline constexpr double kSsiTraceStepFraction = 0.05;

// Upper bound on halvings/retries in the curvature and closing loops.
inline constexpr int kSsiTraceMaxRetries = 60;

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

// Midpoint of two parameters on one axis, measured across the periodic seam
// when the surface wraps on that axis.
inline double ssi_trace_param_midpoint(const SurfaceGeometry& surface, int axis, double first,
                                       double second) noexcept {
  return first + 0.5 * ssi_trace_param_delta(surface, axis, first, second);
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

// True when both parameter pairs are inside their patch domains.
inline bool ssi_trace_inside(const SsiPatch& a, const SsiPatch& b, const std::array<double, 2>& ua,
                             const std::array<double, 2>& va) noexcept {
  return patch_contains(a, ua[0], ua[1]) && patch_contains(b, va[0], va[1]);
}

// Largest 3D step whose parameter motion on either surface stays under
// kSsiTraceStepFraction of every patch span it actually moves along. Only axes
// with non-zero velocity constrain the step, so a line along the long axis of an
// elongated patch is not limited by the patch's narrow direction.
inline double ssi_trace_step_cap(const SsiPatch& a, const SsiPatch& b,
                                 const SsiTraceFrame& frame) noexcept {
  const SsiPatch* patches[2] = {&a, &b};
  const std::array<double, 2>* velocities[2] = {&frame.velocity_a, &frame.velocity_b};
  double cap = std::numeric_limits<double>::infinity();
  for (int surface = 0; surface < 2; ++surface) {
    const std::array<double, 2> spans{patches[surface]->u[1] - patches[surface]->u[0],
                                      patches[surface]->v[1] - patches[surface]->v[0]};
    for (int axis = 0; axis < 2; ++axis) {
      const double speed = std::fabs((*velocities[surface])[axis]);
      if (speed > 0.0) cap = std::min(cap, kSsiTraceStepFraction * spans[axis] / speed);
    }
  }
  return cap;
}

// Refines the patch-boundary endpoint between an accepted inside vertex and a
// rejected prediction that left one of the patches. The parameter interval is
// bisected (circularly on periodic axes) and every midpoint inside both patches
// is refined; the last inside point becomes the endpoint. The loop stops when
// the interval is 1e-12 of a patch span or the endpoint stops moving, so the
// endpoint sits on the boundary instead of up to a full step short of it.
inline SsiPoint ssi_trace_boundary_point(const SsiPatch& a, const SsiPatch& b,
                                         const SsiPoint& inside, std::array<double, 2> outside_ua,
                                         std::array<double, 2> outside_va, double tolerance) {
  std::array<double, 2> low_ua = inside.ua;
  std::array<double, 2> low_va = inside.va;
  SsiPoint endpoint = inside;
  const std::array<double, 2> spans_a{a.u[1] - a.u[0], a.v[1] - a.v[0]};
  const std::array<double, 2> spans_b{b.u[1] - b.u[0], b.v[1] - b.v[0]};
  for (int iteration = 0; iteration < 128; ++iteration) {
    std::array<double, 2> mid_ua{ssi_trace_param_midpoint(a.surface, 0, low_ua[0], outside_ua[0]),
                                 ssi_trace_param_midpoint(a.surface, 1, low_ua[1], outside_ua[1])};
    std::array<double, 2> mid_va{ssi_trace_param_midpoint(b.surface, 0, low_va[0], outside_va[0]),
                                 ssi_trace_param_midpoint(b.surface, 1, low_va[1], outside_va[1])};
    if (!ssi_solver_wrap(a.surface, mid_ua) || !ssi_solver_wrap(b.surface, mid_va)) break;
    bool accepted = false;
    bool converged = false;
    if (ssi_trace_inside(a, b, mid_ua, mid_va)) {
      const Result<SsiPoint> refined =
          refine_ssi_point(a, b, {mid_ua[0], mid_ua[1], mid_va[0], mid_va[1]}, tolerance);
      if (refined.has_value() && ssi_trace_inside(a, b, refined.value().ua, refined.value().va)) {
        // The second stopping rule: an accepted midpoint that moves the endpoint
        // by less than the tolerance is already on the boundary within tolerance.
        converged = norm2(refined.value().p - endpoint.p) <= tolerance;
        endpoint = refined.value();
        low_ua = endpoint.ua;
        low_va = endpoint.va;
        accepted = true;
      }
    }
    if (converged) break;
    if (!accepted) {
      outside_ua = mid_ua;
      outside_va = mid_va;
    }
    double width = 0.0;
    for (int axis = 0; axis < 2; ++axis) {
      width = std::max(
          width, std::fabs(ssi_trace_param_delta(a.surface, axis, low_ua[axis], outside_ua[axis])) /
                     spans_a[axis]);
      width = std::max(
          width, std::fabs(ssi_trace_param_delta(b.surface, axis, low_va[axis], outside_va[axis])) /
                     spans_b[axis]);
    }
    if (width <= 1e-12) break;
  }
  return endpoint;
}

// Why one directional walk stopped. left_patch, refinement_failed and
// tangent_failed are normal open-curve endings; step_budget means the walk
// could not emit its next segment inside the shared budget; closure_failed
// means the near-start refinement could not place a vertex within 8*tolerance
// of the start.
enum class SsiTraceEnd {
  closed,
  left_patch,
  refinement_failed,
  tangent_failed,
  step_budget,
  closure_failed
};

struct SsiTraceHalf {
  std::vector<SsiPoint> points;
  SsiTraceEnd end = SsiTraceEnd::left_patch;
};

// Walks the last part of a closing trace into the start vertex. The curve is
// known to pass through the start, so each trial is a refined step of the
// current 3D gap toward it; trials that do not make clear progress are halved
// and retried from the same vertex. The walk stops as soon as the last accepted
// vertex is within 8*tolerance of the start (the closing weld is then short
// enough), or after pushing one more vertex that close. Every vertex it emits
// is refined, so the closing chord is never O(step).
inline SsiTraceEnd ssi_trace_close(const SsiPatch& a, const SsiPatch& b, SsiPoint base,
                                   SsiTraceFrame base_frame, const SsiPoint& start, double sign,
                                   double step_limit, double tolerance, int max_steps,
                                   int& segments_used, std::vector<SsiPoint>& points) {
  const double target = 8.0 * tolerance;
  double gap = norm2(base.p - start.p);
  if (!(gap > target)) return SsiTraceEnd::closed;
  // Aim for a vertex `target` short of the start rather than at the start
  // itself: a candidate that lands exactly on the start would turn the final
  // step into the closing chord, which is O(step). Stopping one target short
  // keeps the closing weld (the last emitted segment) within 8*tolerance.
  double trial = std::min(step_limit, gap - 0.5 * target);
  for (int attempt = 0; attempt < kSsiTraceMaxRetries; ++attempt) {
    if (segments_used >= max_steps) return SsiTraceEnd::step_budget;
    if (!(trial > 0.0)) return SsiTraceEnd::closure_failed;
    std::array<double, 2> predicted_ua{base.ua[0] + sign * trial * base_frame.velocity_a[0],
                                       base.ua[1] + sign * trial * base_frame.velocity_a[1]};
    std::array<double, 2> predicted_va{base.va[0] + sign * trial * base_frame.velocity_b[0],
                                       base.va[1] + sign * trial * base_frame.velocity_b[1]};
    if (!ssi_solver_wrap(a.surface, predicted_ua) || !ssi_solver_wrap(b.surface, predicted_va)) {
      trial *= 0.5;
      continue;
    }
    if (!ssi_trace_inside(a, b, predicted_ua, predicted_va)) {
      trial *= 0.5;
      continue;
    }
    const Result<SsiPoint> refined = refine_ssi_point(
        a, b, {predicted_ua[0], predicted_ua[1], predicted_va[0], predicted_va[1]}, tolerance);
    if (!refined.has_value()) {
      trial *= 0.5;
      continue;
    }
    const SsiPoint next = refined.value();
    if (!ssi_trace_inside(a, b, next.ua, next.va) || !(next.residual <= tolerance)) {
      trial *= 0.5;
      continue;
    }
    const double next_gap = norm2(next.p - start.p);
    if (next_gap <= target) {
      points.push_back(next);
      ++segments_used;
      return SsiTraceEnd::closed;
    }
    if (next_gap < 0.9 * gap) {
      SsiTraceFrame next_frame;
      if (!ssi_trace_frame(a, b, next, next_frame)) return SsiTraceEnd::tangent_failed;
      points.push_back(next);
      ++segments_used;
      base = next;
      base_frame = next_frame;
      gap = next_gap;
      trial = std::min(step_limit, gap - 0.5 * target);
      continue;
    }
    trial *= 0.5;
  }
  return SsiTraceEnd::closure_failed;
}

// Walks one direction from `start`. Every accepted vertex is the result of
// refine_ssi_point, and every candidate is checked against both patch domains
// and against the tolerance: the gap residual alone does not certify that the
// parameters stayed inside a patch. A step whose prediction leaves either patch
// is refined onto the boundary and ends the walk there. A tangent that turns by
// more than kSsiTraceCurvatureTurn halving-retries the step. The step budget is
// shared by both directions of one curve and counts every emitted segment.
inline SsiTraceHalf ssi_trace_half(const SsiPatch& a, const SsiPatch& b, const SsiPoint& start,
                                   const SsiTraceFrame& start_frame, double sign, double step_start,
                                   double step_max, double tolerance, int max_steps,
                                   int& segments_used, bool allow_closure) {
  SsiTraceHalf half;
  half.points.push_back(start);
  SsiPoint current = start;
  SsiTraceFrame frame = start_frame;
  double step = step_start;
  int retries = 0;
  for (;;) {
    const double cap = ssi_trace_step_cap(a, b, frame);
    if (!std::isfinite(cap) || !(cap > 0.0)) {
      half.end = SsiTraceEnd::tangent_failed;
      return half;
    }
    step = std::min(step, cap);
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
    if (!ssi_trace_inside(a, b, predicted_ua, predicted_va)) {
      // The step leaves a patch: refine the exit point itself so the curve ends
      // on the boundary rather than a whole step short of it.
      const SsiPoint endpoint =
          ssi_trace_boundary_point(a, b, current, predicted_ua, predicted_va, tolerance);
      if (norm2(endpoint.p - current.p) > tolerance) {
        if (segments_used >= max_steps) {
          half.end = SsiTraceEnd::step_budget;
          return half;
        }
        half.points.push_back(endpoint);
        ++segments_used;
      }
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
    // Defensive re-check, not the load-bearing gate: refine_ssi_point already
    // guarantees both patches and a residual within tolerance. It is kept here
    // because this walk depends on those invariants and the cost is one fold.
    if (!ssi_trace_inside(a, b, next.ua, next.va) || !(next.residual <= tolerance)) {
      half.end = SsiTraceEnd::refinement_failed;
      return half;
    }
    SsiTraceFrame next_frame;
    if (!ssi_trace_frame(a, b, next, next_frame)) {
      half.end = SsiTraceEnd::tangent_failed;
      return half;
    }
    if (dot(next_frame.tangent, frame.tangent) < std::cos(kSsiTraceCurvatureTurn)) {
      step *= 0.5;
      if (++retries > kSsiTraceMaxRetries) {
        half.end = SsiTraceEnd::tangent_failed;
        return half;
      }
      continue;
    }
    // Closure is decided in 3D, because a parameter window mixes radians with
    // lengths and breaks as soon as the model scale changes. The tangent gate
    // rejects the far side of a loop, where a long step can pass near the start
    // while still travelling the other way.
    if (allow_closure && half.points.size() >= 4 &&
        norm2(next.p - start.p) <= std::max(2.0 * step, 8.0 * tolerance) &&
        dot(next_frame.tangent, start_frame.tangent) >= kSsiTraceClosureTangent) {
      half.end = ssi_trace_close(a, b, current, frame, start, sign, step, tolerance, max_steps,
                                 segments_used, half.points);
      return half;
    }
    if (segments_used >= max_steps) {
      half.end = SsiTraceEnd::step_budget;
      return half;
    }
    half.points.push_back(next);
    ++segments_used;
    current = next;
    frame = next_frame;
    retries = 0;
    step = std::min(step * 1.5, step_max);
  }
}

}  // namespace detail

// Traces the intersection curves of two validated SSI patches from refined
// seeds. Each curve is walked in both null-space directions; a walk that
// returns to its start closes the curve through a refined vertex within
// 8*tolerance of the start, so the closing chord is never O(step). The result is
// ordered by each curve's first vertex so it does not depend on seed order, and
// duplicated seeds collapse to one curve. A walk that crosses a patch boundary
// ends on the boundary itself; a walk that cannot emit its next segment inside
// the shared budget is a numeric_failure rather than a silently truncated
// curve, and a seed that traces to nothing is reported instead of returning an
// empty success.
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
  const double step_start = diagonal / 64.0;
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
    // The first step is also bounded by the local patch metric (see
    // ssi_trace_step_cap), so a small feature inside a large patch is not
    // stepped over on the first move.
    const double step =
        std::min(step_start, detail::ssi_trace_step_cap(patch_a, patch_b, start_frame));
    if (!std::isfinite(step) || !(step > 0.0)) continue;
    int segments_used = 0;
    detail::SsiTraceHalf forward =
        detail::ssi_trace_half(patch_a, patch_b, origin, start_frame, 1.0, step, step_max,
                               tolerance, max_steps, segments_used, true);
    if (forward.end == detail::SsiTraceEnd::step_budget) {
      return Error{Errc::numeric_failure, "trace_ssi_curves: step budget exhausted"};
    }
    if (forward.end == detail::SsiTraceEnd::closure_failed) {
      return Error{Errc::numeric_failure, "trace_ssi_curves: closure refinement did not converge"};
    }
    std::vector<SsiPoint> points;
    bool closed = false;
    if (forward.end == detail::SsiTraceEnd::closed) {
      points = forward.points;
      closed = true;
    } else {
      detail::SsiTraceHalf backward =
          detail::ssi_trace_half(patch_a, patch_b, origin, start_frame, -1.0, step, step_max,
                                 tolerance, max_steps, segments_used, false);
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
    if (points.size() < 2) {
      // This seed has no extent, but other seeds may still describe real
      // curves; only an all-empty result is reported below.
      continue;
    }
    const bool weld = closed && norm2(points.back().p - points.front().p) > 0.0;
    if (weld) {
      if (segments_used >= max_steps) {
        return Error{Errc::numeric_failure, "trace_ssi_curves: step budget exhausted"};
      }
      ++segments_used;
    }
    SsiCurve curve;
    curve.closed = closed;
    const std::size_t segment_count =
        closed ? (weld ? points.size() : points.size() - 1) : points.size() - 1;
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

  // A refined seed whose tangent frame is degenerate has no traceable extent
  // (for example a sphere and a plane that touch at a pole). Returning success
  // with no curves would hide that from the caller, so it is reported instead.
  if (curves.empty() && !seeds.empty()) {
    return Error{Errc::numeric_failure, "trace_ssi_curves: seed has no traceable extent"};
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

// ---------------------------------------------------------------------------
// Degenerate-input diagnostics
// ---------------------------------------------------------------------------
//
// surfaces_coincident and patches_tangent classify the two cases SSI cannot
// trace: the patches lie on the same infinite carrier (infinitely many
// intersection curves), or the carriers touch with parallel normals (an
// isolated contact or a tangency curve, not a transverse intersection).
//
// Both predicates are noexcept and allocation-free. They answer false, rather
// than reporting an error, for any input the sampling below cannot certify: a
// non-finite coordinate, a non-finite or non-positive tolerance, an invalid or
// inverted patch, or a NURBS patch. NURBS patches are not diagnosed (false for
// both predicates): their exact implicit distance is not available without a
// Newton projection, so neither coincidence nor tangency can be certified from
// a sampled point set alone. A mixed surface-type pair is never coincident (two
// different carrier families are never the same surface) but can be tangent.

namespace detail {

// Sampling resolution shared by both diagnostics. Five samples per axis cover
// the four patch corners and the parameter midpoint.
inline constexpr int kSsiDiagnosticGrid = 5;
inline constexpr int kSsiDiagnosticSampleCount = kSsiDiagnosticGrid * kSsiDiagnosticGrid;

// Matching slack for carrier frames, measured on 1 - |unit axis dot unit axis|.
inline constexpr double kSsiDiagnosticFrameSlack = 1e-12;

// Angular slack for the tangency contact, measured on |n_a x n_b| of unit normals.
inline constexpr double kSsiDiagnosticTangentSlack = 1e-6;

inline constexpr double kSsiDiagnosticInfinity = std::numeric_limits<double>::infinity();

struct SsiDiagnosticSample {
  Point3d point{};
  std::array<double, 2> uv{};
};

inline bool ssi_diagnostic_vec_finite(const Vec3d& value) noexcept {
  return std::isfinite(value[0]) && std::isfinite(value[1]) && std::isfinite(value[2]);
}

inline Vec3d ssi_diagnostic_vec(const Point3d& point) noexcept {
  return Vec3d{point[0], point[1], point[2]};
}

inline bool ssi_diagnostic_unit_parallel(const Vec3d& a, const Vec3d& b) noexcept {
  const double alignment = std::fabs(dot(a, b));
  return std::isfinite(alignment) && 1.0 - alignment <= kSsiDiagnosticFrameSlack;
}

// Unit carrier data for a Plane3: the unit normal and the offset along it. The
// plane frame normalizes the caller's normal exactly once, so planes whose
// normals differ only in length compare equal here.
struct SsiDiagnosticPlane {
  Vec3d normal{};
  double offset = 0.0;
};

inline bool ssi_diagnostic_plane(const Plane3& plane, SsiDiagnosticPlane& out) noexcept {
  const Result<PlaneFrame> frame = make_plane_frame(plane);
  if (!frame.has_value()) return false;
  out.normal = frame.value().normal;
  out.offset = dot(out.normal, ssi_diagnostic_vec(frame.value().origin));
  return std::isfinite(out.offset);
}

// A patch is usable for diagnostics when its box is finite and ordered, does
// not exceed one periodic turn on a periodic axis, and its surface validates.
inline bool ssi_diagnostic_patch_valid(const SsiPatch& patch) noexcept {
  if (!std::isfinite(patch.u[0]) || !std::isfinite(patch.u[1]) ||
      !std::isfinite(patch.v[0]) || !std::isfinite(patch.v[1])) {
    return false;
  }
  if (!(patch.u[0] < patch.u[1]) || !(patch.v[0] < patch.v[1])) return false;
  if (surface_periodic(patch.surface, 0) && patch.u[1] - patch.u[0] > kSsiTwoPi) {
    return false;
  }
  if (surface_periodic(patch.surface, 1) && patch.v[1] - patch.v[0] > kSsiTwoPi) {
    return false;
  }
  return ssi_validate_surface(patch.surface).has_value();
}

inline bool ssi_diagnostic_frames_match(const Vec3d& axis_a, const Vec3d& x_a,
                                        const Vec3d& axis_b, const Vec3d& x_b) noexcept {
  return ssi_diagnostic_unit_parallel(axis_a, axis_b) &&
         ssi_diagnostic_unit_parallel(x_a, x_b);
}

// Same infinite carrier for two patches of one surface type. The caller has
// already checked that the variant indices are equal, so the second get_if
// always succeeds here.
inline bool ssi_diagnostic_same_carrier(const SurfaceGeometry& a, const SurfaceGeometry& b,
                                        double tolerance) noexcept {
  if (const Plane3* plane_a = std::get_if<Plane3>(&a)) {
    const Plane3* plane_b = std::get_if<Plane3>(&b);
    if (plane_b == nullptr) return false;
    SsiDiagnosticPlane unit_a{};
    SsiDiagnosticPlane unit_b{};
    if (!ssi_diagnostic_plane(*plane_a, unit_a) || !ssi_diagnostic_plane(*plane_b, unit_b)) {
      return false;
    }
    if (!ssi_diagnostic_unit_parallel(unit_a.normal, unit_b.normal)) return false;
    // One carrier can be written with opposite normals; fold the sign into B
    // before comparing the offsets.
    const double sign = dot(unit_a.normal, unit_b.normal) < 0.0 ? -1.0 : 1.0;
    return std::fabs(unit_a.offset - sign * unit_b.offset) <= tolerance;
  }
  if (const Sphere3* sphere_a = std::get_if<Sphere3>(&a)) {
    const Sphere3* sphere_b = std::get_if<Sphere3>(&b);
    if (sphere_b == nullptr) return false;
    return norm2(sphere_a->center - sphere_b->center) <= tolerance &&
           std::fabs(sphere_a->radius - sphere_b->radius) <= tolerance &&
           ssi_diagnostic_frames_match(sphere_a->axis, sphere_a->x_axis, sphere_b->axis,
                                       sphere_b->x_axis);
  }
  if (const Cylinder3* cylinder_a = std::get_if<Cylinder3>(&a)) {
    const Cylinder3* cylinder_b = std::get_if<Cylinder3>(&b);
    if (cylinder_b == nullptr) return false;
    if (std::fabs(cylinder_a->radius - cylinder_b->radius) > tolerance) return false;
    if (!ssi_diagnostic_frames_match(cylinder_a->axis, cylinder_a->x_axis, cylinder_b->axis,
                                     cylinder_b->x_axis)) {
      return false;
    }
    // Sliding the origin along the axis does not move a cylinder; only the
    // radial offset between the two axis lines matters.
    const Vec3d delta = cylinder_b->origin - cylinder_a->origin;
    const Vec3d radial = delta - cylinder_a->axis * dot(delta, cylinder_a->axis);
    return norm2(radial) <= tolerance;
  }
  if (const Cone3* cone_a = std::get_if<Cone3>(&a)) {
    const Cone3* cone_b = std::get_if<Cone3>(&b);
    if (cone_b == nullptr) return false;
    if (std::fabs(cone_a->half_angle - cone_b->half_angle) > tolerance) return false;
    if (!ssi_diagnostic_frames_match(cone_a->axis, cone_a->x_axis, cone_b->axis,
                                     cone_b->x_axis)) {
      return false;
    }
    // (origin, radius) is a redundant parametrization of one nappe; the apex is
    // the invariant. Folding it this way also handles a flipped axis, which
    // moves the apex in the opposite direction along the axis.
    const double apex_a = cone_a->radius / std::tan(cone_a->half_angle);
    const double apex_b = cone_b->radius / std::tan(cone_b->half_angle);
    const Point3d point_a = cone_a->origin + cone_a->axis * apex_a;
    const Point3d point_b = cone_b->origin + cone_b->axis * apex_b;
    return norm2(point_a - point_b) <= tolerance;
  }
  if (const Torus3* torus_a = std::get_if<Torus3>(&a)) {
    const Torus3* torus_b = std::get_if<Torus3>(&b);
    if (torus_b == nullptr) return false;
    return norm2(torus_a->center - torus_b->center) <= tolerance &&
           std::fabs(torus_a->major_radius - torus_b->major_radius) <= tolerance &&
           std::fabs(torus_a->minor_radius - torus_b->minor_radius) <= tolerance &&
           ssi_diagnostic_frames_match(torus_a->axis, torus_a->x_axis, torus_b->axis,
                                       torus_b->x_axis);
  }
  return false;  // NurbsSurface: coincidence is not diagnosed.
}

// Exact implicit distance from point to the infinite carrier (zero on the
// surface), used to sample the patch domain. Plane/sphere/cylinder/torus
// distances are exact; the cone expression is exact for the double cone. A
// NURBS carrier has no exact implicit form here, so it reports infinity and
// every caller rejects it.
inline double ssi_diagnostic_carrier_distance(const SurfaceGeometry& surface,
                                              const Point3d& point) noexcept {
  return std::visit(
      [&point](const auto& carrier) noexcept -> double {
        using Surface = std::decay_t<decltype(carrier)>;
        if constexpr (std::is_same_v<Surface, Plane3>) {
          SsiDiagnosticPlane plane{};
          if (!ssi_diagnostic_plane(carrier, plane)) return kSsiDiagnosticInfinity;
          return std::fabs(dot(ssi_diagnostic_vec(point), plane.normal) - plane.offset);
        } else if constexpr (std::is_same_v<Surface, Sphere3>) {
          return std::fabs(norm2(point - carrier.center) - carrier.radius);
        } else if constexpr (std::is_same_v<Surface, Cylinder3>) {
          const Vec3d local = point - carrier.origin;
          const double axial = dot(local, carrier.axis);
          return std::fabs(norm2(local - carrier.axis * axial) - carrier.radius);
        } else if constexpr (std::is_same_v<Surface, Cone3>) {
          const Vec3d local = point - carrier.origin;
          const double axial = dot(local, carrier.axis);
          const double radial = norm2(local - carrier.axis * axial);
          const double from_apex = axial - carrier.radius / std::tan(carrier.half_angle);
          return std::fabs(radial * std::cos(carrier.half_angle) -
                           std::fabs(from_apex) * std::sin(carrier.half_angle));
        } else if constexpr (std::is_same_v<Surface, Torus3>) {
          const Vec3d local = point - carrier.center;
          const double axial = dot(local, carrier.axis);
          const double radial = norm2(local - carrier.axis * axial);
          return std::fabs(
              std::hypot(radial - carrier.major_radius, axial) - carrier.minor_radius);
        } else {
          return kSsiDiagnosticInfinity;
        }
      },
      surface);
}

// Sample the patch on a full 5x5 grid, including both ends of each parameter
// interval. Non-finite evaluation or non-finite resulting points fail the
// sample; callers treat that as "not diagnosed" rather than as an error.
inline bool ssi_diagnostic_sample_grid(
    const SsiPatch& patch,
    std::array<SsiDiagnosticSample, kSsiDiagnosticSampleCount>& samples) noexcept {
  for (int iu = 0; iu < kSsiDiagnosticGrid; ++iu) {
    const double u =
        iu + 1 == kSsiDiagnosticGrid
            ? patch.u[1]
            : patch.u[0] + (patch.u[1] - patch.u[0]) *
                               (static_cast<double>(iu) /
                                static_cast<double>(kSsiDiagnosticGrid - 1));
    for (int iv = 0; iv < kSsiDiagnosticGrid; ++iv) {
      const double v =
          iv + 1 == kSsiDiagnosticGrid
              ? patch.v[1]
              : patch.v[0] + (patch.v[1] - patch.v[0]) *
                                 (static_cast<double>(iv) /
                                  static_cast<double>(kSsiDiagnosticGrid - 1));
      const Result<SurfacePoint> evaluated = surface_eval(patch.surface, u, v);
      if (!evaluated.has_value() || !surface_point_finite(evaluated.value())) return false;
      SsiDiagnosticSample& sample =
          samples[static_cast<std::size_t>(iu) * kSsiDiagnosticGrid + iv];
      sample.point = evaluated.value().p;
      sample.uv = {u, v};
    }
  }
  return true;
}

}  // namespace detail

// True when both patches lie on the same infinite carrier within `tolerance`.
//
// The answer is certified over the sampled patch domain rather than at one
// point: the defining carrier parameters must match (unit frames within
// 1 - |dot| <= 1e-12, and offsets/centres/radii/apex within `tolerance`), and
// all 25 samples of patch a must lie within `tolerance` of carrier b's exact
// implicit distance. NURBS patches answer false because their implicit
// distance is unavailable; a mixed plane/analytic pair never shares a carrier.
// Non-finite coordinates and a non-positive or non-finite tolerance answer
// false as well. The function is noexcept, allocation-free, and deterministic.
inline bool surfaces_coincident(const SsiPatch& a, const SsiPatch& b,
                                double tolerance) noexcept {
  if (!std::isfinite(tolerance) || !(tolerance > 0.0)) return false;
  if (a.surface.index() != b.surface.index()) return false;
  if (std::holds_alternative<NurbsSurface>(a.surface)) return false;
  if (!detail::ssi_diagnostic_patch_valid(a) || !detail::ssi_diagnostic_patch_valid(b)) {
    return false;
  }
  if (!detail::ssi_diagnostic_same_carrier(a.surface, b.surface, tolerance)) return false;
  std::array<detail::SsiDiagnosticSample, detail::kSsiDiagnosticSampleCount> samples{};
  if (!detail::ssi_diagnostic_sample_grid(a, samples)) return false;
  for (const detail::SsiDiagnosticSample& sample : samples) {
    const double distance = detail::ssi_diagnostic_carrier_distance(b.surface, sample.point);
    if (!std::isfinite(distance) || distance > tolerance) return false;
  }
  return true;
}

// True when the two carriers touch with parallel normals, without coinciding.
//
// The sampled contact gap is the minimum over all 25x25 grid pairs of
// max(distance(p_a, carrier_b), distance(p_b, carrier_a)); it must be within
// `tolerance`. At the winning pair the unit normals must satisfy
// |n_a x n_b| <= 1e-6. Coincident carriers, crossing planes, and a sphere cut
// by a plane in a circle therefore answer false: coincidence is excluded
// explicitly, and a transverse intersection has non-parallel normals. NURBS
// patches answer false because their implicit distance is unavailable, so
// tangency cannot be certified (and coincidence cannot be excluded). The same
// non-finite/rejected-input rules as surfaces_coincident apply. The function is
// noexcept, allocation-free, and deterministic.
inline bool patches_tangent(const SsiPatch& a, const SsiPatch& b, double tolerance) noexcept {
  if (!std::isfinite(tolerance) || !(tolerance > 0.0)) return false;
  if (std::holds_alternative<NurbsSurface>(a.surface) ||
      std::holds_alternative<NurbsSurface>(b.surface)) {
    return false;
  }
  if (!detail::ssi_diagnostic_patch_valid(a) || !detail::ssi_diagnostic_patch_valid(b)) {
    return false;
  }
  if (surfaces_coincident(a, b, tolerance)) return false;

  std::array<detail::SsiDiagnosticSample, detail::kSsiDiagnosticSampleCount> samples_a{};
  std::array<detail::SsiDiagnosticSample, detail::kSsiDiagnosticSampleCount> samples_b{};
  if (!detail::ssi_diagnostic_sample_grid(a, samples_a) ||
      !detail::ssi_diagnostic_sample_grid(b, samples_b)) {
    return false;
  }

  double contact = detail::kSsiDiagnosticInfinity;
  std::size_t contact_a = 0;
  std::size_t contact_b = 0;
  for (std::size_t ia = 0; ia < samples_a.size(); ++ia) {
    for (std::size_t ib = 0; ib < samples_b.size(); ++ib) {
      const double to_b =
          detail::ssi_diagnostic_carrier_distance(b.surface, samples_a[ia].point);
      const double to_a =
          detail::ssi_diagnostic_carrier_distance(a.surface, samples_b[ib].point);
      const double gap = std::max(to_a, to_b);
      if (!std::isfinite(gap)) return false;
      if (gap < contact) {
        contact = gap;
        contact_a = ia;
        contact_b = ib;
      }
    }
  }
  if (!(contact <= tolerance)) return false;

  const Result<Vec3d> normal_a =
      surface_normal(a.surface, samples_a[contact_a].uv[0], samples_a[contact_a].uv[1]);
  const Result<Vec3d> normal_b =
      surface_normal(b.surface, samples_b[contact_b].uv[0], samples_b[contact_b].uv[1]);
  if (!normal_a.has_value() || !normal_b.has_value()) return false;
  const Vec3d normal_cross = cross(normal_a.value(), normal_b.value());
  if (!detail::ssi_diagnostic_vec_finite(normal_cross)) return false;
  return norm2(normal_cross) <= detail::kSsiDiagnosticTangentSlack;
}

}  // namespace v1
}  // namespace geom
