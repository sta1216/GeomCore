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
#include <map>
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

// Period-aware distance between two parameter pairs on one patch: a periodic
// axis is compared modulo 2*pi, so a seam pair such as -pi+eps and pi-eps is
// 2*eps apart rather than almost a full period. Non-periodic axes are compared
// as plain differences. The result is in parameter units, not radians alone:
// it is the Euclidean length of the folded (du, dv) vector.
inline double ssi_param_distance(const SsiPatch& patch, std::array<double, 2> p,
                                 std::array<double, 2> q) noexcept {
  std::array<double, 2> delta{p[0] - q[0], p[1] - q[1]};
  for (std::size_t axis = 0; axis < 2; ++axis) {
    if (surface_periodic(patch.surface, static_cast<int>(axis))) {
      delta[axis] = std::remainder(delta[axis], detail::kSsiTwoPi);
    }
  }
  return std::hypot(delta[0], delta[1]);
}

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

// Refines an intersection point pinned to one parameter of patch a. A walk that
// leaves the patch has no unconstrained intersection point past the limit, so
// this solve keeps ua[axis] at the boundary value and finds the remaining
// (u, s, t) that satisfy SA = SB. It is the boundary analogue of
// refine_ssi_point; when the pinned system is singular the caller falls back to
// the parameter bisection.
inline Result<SsiPoint> ssi_trace_boundary_refine(const SsiPatch& a, const SsiPatch& b,
                                                  std::array<double, 2> ua,
                                                  std::array<double, 2> va, int axis,
                                                  double boundary, double tolerance) {
  const std::size_t pinned = static_cast<std::size_t>(axis);
  ua[pinned] = boundary;
  if (!detail::ssi_solver_wrap(a.surface, ua) ||
      !detail::ssi_solver_wrap(b.surface, va) || ua[pinned] != boundary) {
    return Error{Errc::numeric_failure, "ssi: boundary seed cannot be pinned"};
  }
  if (!patch_contains(a, ua[0], ua[1]) || !patch_contains(b, va[0], va[1])) {
    return Error{Errc::numeric_failure, "ssi: boundary seed lies outside a patch domain"};
  }
  const std::size_t other = pinned == 0 ? 1 : 0;
  Result<SurfacePoint> on_a = surface_eval(a.surface, ua[0], ua[1]);
  if (!on_a.has_value()) return on_a.error();
  Result<SurfacePoint> on_b = surface_eval(b.surface, va[0], va[1]);
  if (!on_b.has_value()) return on_b.error();
  double gap = norm2(on_a.value().p - on_b.value().p);
  for (int iteration = 0; iteration < 24 && gap > tolerance; ++iteration) {
    const Vec3d mismatch = on_a.value().p - on_b.value().p;
    const Vec3d column_a = pinned == 0 ? on_a.value().dv : on_a.value().du;
    const Vec3d column_bu{-on_b.value().du[0], -on_b.value().du[1], -on_b.value().du[2]};
    const Vec3d column_bv{-on_b.value().dv[0], -on_b.value().dv[1], -on_b.value().dv[2]};
    std::array<std::array<double, 3>, 3> system{};
    for (int row = 0; row < 3; ++row) {
      const std::size_t at = static_cast<std::size_t>(row);
      system[at][0] = column_a[row];
      system[at][1] = column_bu[row];
      system[at][2] = column_bv[row];
    }
    std::array<double, 3> solution{};
    if (!detail::ssi_solver_solve3(system, Vec3d{-mismatch[0], -mismatch[1], -mismatch[2]},
                                   solution)) {
      return Error{Errc::numeric_failure, "ssi: boundary system is numerically singular"};
    }
    bool accepted = false;
    double alpha = 1.0;
    for (int halving = 0; halving <= detail::kSsiSolverMaxHalvings; ++halving, alpha *= 0.5) {
      std::array<double, 2> trial_ua = ua;
      trial_ua[other] += alpha * solution[0];
      std::array<double, 2> trial_va{va[0] + alpha * solution[1], va[1] + alpha * solution[2]};
      bool usable = std::isfinite(trial_ua[0]) && std::isfinite(trial_ua[1]) &&
                    std::isfinite(trial_va[0]) && std::isfinite(trial_va[1]);
      if (usable) usable = detail::ssi_solver_wrap(a.surface, trial_ua);
      if (usable) usable = detail::ssi_solver_wrap(b.surface, trial_va);
      if (usable) usable = trial_ua[pinned] == boundary;
      if (usable) usable = patch_contains(a, trial_ua[0], trial_ua[1]);
      if (usable) usable = patch_contains(b, trial_va[0], trial_va[1]);
      if (!usable) continue;
      const Result<SurfacePoint> trial_a = surface_eval(a.surface, trial_ua[0], trial_ua[1]);
      if (!trial_a.has_value()) return trial_a.error();
      const Result<SurfacePoint> trial_b = surface_eval(b.surface, trial_va[0], trial_va[1]);
      if (!trial_b.has_value()) return trial_b.error();
      const double trial_gap = norm2(trial_a.value().p - trial_b.value().p);
      if (!std::isfinite(trial_gap) || !(trial_gap < gap)) continue;
      accepted = true;
      ua = trial_ua;
      va = trial_va;
      on_a = trial_a;
      on_b = trial_b;
      gap = trial_gap;
      break;
    }
    if (!accepted) {
      return Error{Errc::numeric_failure, "ssi: boundary Newton step made no progress"};
    }
  }
  if (!(gap <= tolerance)) {
    return Error{Errc::numeric_failure, "ssi: boundary Newton iteration did not converge"};
  }
  const SsiPoint point{ua, va, on_a.value().p, gap};
  return SsiPoint{point.ua, point.va, point.p, ssi_residual(a, b, point)};
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
  bool converged_flag = false;
  // The curve can reach the limit while the parameter-space midpoint refine
  // slides along the curve, so pin the violated coordinate and solve for the
  // boundary intersection directly before falling back to bisection.
  for (int axis = 0; axis < 2; ++axis) {
    const std::size_t at = static_cast<std::size_t>(axis);
    if (ssi_patch_axis_contains(a, axis, outside_ua[at])) continue;
    const double low = axis == 0 ? a.u[0] : a.v[0];
    const double high = axis == 0 ? a.u[1] : a.v[1];
    // The coordinate left the interval, but on a periodic axis "below low" can
    // also mean "past high, a full period over": a folded value just below zero
    // may sit on the far side of a wedge cut from a nearly full-turn patch.
    // Pin the boundary that is nearest to the inside vertex in parameter space.
    double boundary = low;
    double best = std::fabs(ssi_trace_param_delta(a.surface, axis, inside.ua[at], low));
    const double high_distance =
        std::fabs(ssi_trace_param_delta(a.surface, axis, inside.ua[at], high));
    if (high_distance < best) {
      best = high_distance;
      boundary = high;
    }
    const Result<SsiPoint> pinned =
        ssi_trace_boundary_refine(a, b, inside.ua, inside.va, axis, boundary, tolerance);
    if (pinned.has_value()) return pinned.value();
  }
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
    if (converged) { converged_flag = true; break; }
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
// could not emit its next segment inside the shared budget; overrun means it
// travelled more than 1.25 patch parameter diagonals without closing or
// reaching a boundary; closure_failed means the near-start refinement could not
// place a vertex within 8*tolerance of the start.
enum class SsiTraceEnd {
  closed,
  left_patch,
  refinement_failed,
  tangent_failed,
  step_budget,
  overrun,
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
  // A closing vertex must agree with the start in parameters as well as in 3D:
  // on a periodic patch a point can be within `target` of the start while its
  // stored parameter is a full period away. The window is the smaller of a
  // multiple of the 3D tolerance and a relative parameter span, so it scales
  // with both the model and the patch domain.
  const double span_a = std::hypot(a.u[1] - a.u[0], a.v[1] - a.v[0]);
  const double span_b = std::hypot(b.u[1] - b.u[0], b.v[1] - b.v[0]);
  const double param_window =
      std::max(32.0 * tolerance, 1e-6 * std::max(span_a, span_b));
  double gap = norm2(base.p - start.p);
  if (!(gap > target)) {
    // The base already sits inside the weld band. That is only a closure when
    // the parameters agree too: on a periodic patch the same 3D point can be
    // stored a full period away from the start's parameter branch, and a 3D-only
    // test would accept a closing chord the caller cannot weld unambiguously.
    if (ssi_param_distance(a, base.ua, start.ua) <= param_window &&
        ssi_param_distance(b, base.va, start.va) <= param_window) {
      return SsiTraceEnd::closed;
    }
    return SsiTraceEnd::closure_failed;
  }
  // Aim for a vertex `target` short of the start rather than at the start
  // itself: a candidate that lands exactly on the start would turn the final
  // step into the closing chord, which is O(step). Stopping one target short
  // keeps the closing weld (the last emitted segment) within 8*tolerance.
  double trial = std::min(step_limit, gap - 0.5 * target);
  double direction = sign;
  for (int attempt = 0; attempt < kSsiTraceMaxRetries; ++attempt) {
    if (segments_used >= max_steps) return SsiTraceEnd::step_budget;
    if (!(trial > 0.0)) return SsiTraceEnd::closure_failed;
    std::array<double, 2> predicted_ua{base.ua[0] + direction * trial * base_frame.velocity_a[0],
                                       base.ua[1] + direction * trial * base_frame.velocity_a[1]};
    std::array<double, 2> predicted_va{base.va[0] + direction * trial * base_frame.velocity_b[0],
                                       base.va[1] + direction * trial * base_frame.velocity_b[1]};
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
      if (ssi_param_distance(a, next.ua, start.ua) <= param_window &&
          ssi_param_distance(b, next.va, start.va) <= param_window) {
        points.push_back(next);
        ++segments_used;
        return SsiTraceEnd::closed;
      }
      // The 3D gap closed but the parameters did not: the candidate names a
      // different parameter branch of the same point. Back off instead of
      // accepting a closure the caller could not weld unambiguously.
      trial *= 0.5;
      continue;
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
    // The trial overshot the start or stalled short of it: the tangent chord
    // does not line up with the remaining gap any more. Walk the other way for
    // the next trial, halving the step so a stalled approach still converges.
    direction = -direction;
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
// `param_diagonal` bounds the accumulated parameter travel: a walk that winds
// more than 1.25 diagonals without closing or reaching a boundary is not
// progressing around anything, and is stopped instead of returning a many-lap
// arc as if it were a real open curve.
inline SsiTraceHalf ssi_trace_half(const SsiPatch& a, const SsiPatch& b, const SsiPoint& start,
                                   const SsiTraceFrame& start_frame, double sign, double step_start,
                                   double step_max, double tolerance, double param_diagonal,
                                   int max_steps, int& segments_used, bool allow_closure) {
  SsiTraceHalf half;
  half.points.push_back(start);
  SsiPoint current = start;
  SsiTraceFrame frame = start_frame;
  double step = step_start;
  double travel = 0.0;
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
    // A step can start and end inside a patch while crossing a trimmed part of
    // the parameter domain, e.g. the wedge cut out of a nearly full-turn patch:
    // the end points are legal but the segment between them is not. Sample the
    // midpoint as well, so a coarse step cannot jump across the cut.
    const std::array<double, 2> mid_ua{
        detail::ssi_trace_param_midpoint(a.surface, 0, current.ua[0], predicted_ua[0]),
        detail::ssi_trace_param_midpoint(a.surface, 1, current.ua[1], predicted_ua[1])};
    const std::array<double, 2> mid_va{
        detail::ssi_trace_param_midpoint(b.surface, 0, current.va[0], predicted_va[0]),
        detail::ssi_trace_param_midpoint(b.surface, 1, current.va[1], predicted_va[1])};
    if (!ssi_trace_inside(a, b, mid_ua, mid_va)) {
      // Halfway is already outside: the boundary lies between the current
      // vertex and the midpoint, so clip there instead of on the far side.
      const SsiPoint endpoint = ssi_trace_boundary_point(a, b, current, mid_ua, mid_va, tolerance);
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
    // Parameter travel of this accepted step, measured on the better-conditioned
    // of the two surfaces: a walk that keeps going around a loop it cannot close
    // must be reported, not returned as a very long open arc.
    double step_travel = 0.0;
    for (int axis = 0; axis < 2; ++axis) {
      if (surface_periodic(a.surface, axis)) {
        step_travel = std::max(step_travel, std::fabs(detail::ssi_trace_param_delta(
                                                   a.surface, axis, current.ua[static_cast<std::size_t>(axis)],
                                                   next.ua[static_cast<std::size_t>(axis)])));
      }
      if (surface_periodic(b.surface, axis)) {
        step_travel = std::max(step_travel, std::fabs(detail::ssi_trace_param_delta(
                                                   b.surface, axis, current.va[static_cast<std::size_t>(axis)],
                                                   next.va[static_cast<std::size_t>(axis)])));
      }
    }
    travel += step_travel;
    if (travel > 1.25 * param_diagonal) {
      half.end = SsiTraceEnd::overrun;
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
    // The 3D gate below bounds the remaining gap to the start; the closing walk
    // in ssi_trace_close then enforces the parameter window as well, so a point
    // one period away in parameter space is not accepted as a closure. The
    // tangent gate rejects the far side of a loop, where a long step can pass
    // near the start while still travelling the other way.
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
  const double param_span_a =
      std::hypot(patch_a.u[1] - patch_a.u[0], patch_a.v[1] - patch_a.v[0]);
  const double param_span_b =
      std::hypot(patch_b.u[1] - patch_b.u[0], patch_b.v[1] - patch_b.v[0]);
  const double param_diagonal = std::max(param_span_a, param_span_b);
  if (!std::isfinite(param_diagonal) || !(param_diagonal > 0.0)) {
    return Error{Errc::numeric_failure, "ssi: trace parameter domain is degenerate"};
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
                               tolerance, param_diagonal, max_steps, segments_used, true);
    if (forward.end == detail::SsiTraceEnd::step_budget) {
      return Error{Errc::numeric_failure, "trace_ssi_curves: step budget exhausted"};
    }
    if (forward.end == detail::SsiTraceEnd::overrun) {
      return Error{Errc::numeric_failure,
                   "trace_ssi_curves: curve does not close or reach a boundary"};
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
                                 tolerance, param_diagonal, max_steps, segments_used, false);
      if (backward.end == detail::SsiTraceEnd::step_budget) {
        return Error{Errc::numeric_failure, "trace_ssi_curves: step budget exhausted"};
      }
      if (backward.end == detail::SsiTraceEnd::overrun) {
        return Error{Errc::numeric_failure,
                     "trace_ssi_curves: curve does not close or reach a boundary"};
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
// intersection curves), or the carriers meet with parallel normals and do not
// cross (an isolated contact or a tangency curve rather than a transverse
// intersection).
//
// surfaces_coincident compares the exact defining parameters of the two
// carriers (unit frame directions within 1 - |dot| <= 1e-12, positions and
// radii within the caller's tolerance, with the seam parameter x_axis ignored
// for rotationally symmetric carriers) and then certifies the patch domain on
// a 9x9 sample grid of patch a: the 5x5 node grid plus the 4x4 cell centres,
// i.e. the fractions k/8, k = 0..8, of each parameter interval. Every sample
// must lie within `tolerance` of carrier b's exact implicit distance. Because
// the carriers are infinite, this is an exact parameter comparison plus a
// domain check, not a surface-to-surface fitting heuristic.
//
// patches_tangent refines the minimum point-to-point distance between the two
// patches: a 9x9 x 9x9 coarse grid finds a starting pair, then Gauss-Newton on
// the four patch parameters (3x4 Jacobian, damped 4x4 normal equations, pivot
// guard, step halving, periodic wrap, parameters kept in-patch, at most 50
// iterations) converges to the candidate contact. A second alternating pass then
// moves the pair onto the tangency set, because a line or circle contact leaves
// the gap degenerate along the contact: a damped Gauss-Newton step on the normal
// mismatch (finite-difference normal Jacobian) is followed by ordinary contact
// Gauss-Newton that re-closes the gap, and the mismatch at least halves per
// round. Tangency then requires:
//  - the refined gap is within `tolerance`;
//  - the unit normals at the refined contact are parallel: |cross| < 1e-10;
//  - the patches are not coincident (surfaces_coincident is false);
//  - both refined contact parameters are inside their patches (patch_contains);
//  - a global one-sided test over each whole patch domain. A fixed 21x21
//    parameter grid (endpoints included; periodic axes wrap by 2*pi) is
//    evaluated on the patch, and every sample's exact signed distance to the
//    other carrier must stay on one side. Two bands are checked: the tolerance
//    band `tolerance`, and a geometric band max(1e-12 * (1 + max|p|),
//    0.05 * tolerance). A pair is rejected when either band has samples beyond
//    both +band and -band.
//
//    The geometric band is also the detection floor of this predicate, and the
//    floor is deliberate: a cut shallower than
//    max(1e-12 * (1 + max|p|), 0.05 * tolerance) leaves every sample of the cut
//    carrier inside the band, so such a cut is reported as tangent_contact even
//    though the carriers technically cross. The floor is a tolerance-scale
//    cutoff (5% of the caller's tolerance, never below double round-off at the
//    sample scale), not an accident: cuts above the tested floor are rejected,
//    exact contact is accepted, and the trade-off is the same one the caller
//    makes when choosing `tolerance`. In particular the predicate does NOT
//    claim to catch arbitrarily shallow cuts.
//
// NURBS patches are unsupported by both predicates: their exact signed
// implicit distance is not available without a Newton projection, so neither
// coincidence nor tangency can be certified (both answer false). Mixed
// analytic types are never coincident; they can still be tangent (a plane
// touching a sphere), which only the tangency path allows. Both predicates are
// noexcept and allocation-free, ignore non-finite coordinates, reject
// non-finite or non-positive tolerances, and are deterministic (no global
// state).
// ---------------------------------------------------------------------------

namespace detail {

// Coincidence domain sampling: the 5x5 node grid plus the 4x4 cell centres of
// its cells; together these are the fractions k/8 (k = 0..8) of each parameter
// interval, i.e. a 9x9 grid. The nodes are the contract minimum; the centres
// keep a carrier from passing only because the nodes happen to land in band.
inline constexpr int kSsiDiagnosticCoincidenceGrid = 9;

// Contact search: a 9x9 coarse grid per patch, then Gauss-Newton refinement.
inline constexpr int kSsiDiagnosticContactGrid = 9;
inline constexpr int kSsiDiagnosticRefineIterations = 50;
inline constexpr double kSsiDiagnosticPivotSlack = 1e-14;
inline constexpr double kSsiDiagnosticDamping = 1e-10;

// Tangency alignment: the gap refinement stalls on a line or circle contact at
// a parameter drift whose size is the normal tilt, so a second pass alternates
// a damped Gauss-Newton step on the normal mismatch (finite-difference normal
// Jacobian) with contact Gauss-Newton that re-closes the gap. Each round halves
// the mismatch, so 32 rounds take a 1e-6 starting mismatch down to 1e-16.
inline constexpr int kSsiDiagnosticAlignRounds = 32;
inline constexpr int kSsiDiagnosticAlignGaussNewtonIterations = 4;
inline constexpr double kSsiDiagnosticAlignStep = 1e-6;

// Global one-sided test: a fixed 21x21 parameter grid over the whole patch
// domain, and two neutral bands: the caller's tolerance and a geometric band
// that is 0.05 * tolerance above the round-off floor. The first catches a
// crossing that separates the carriers by a full tolerance or more; the second
// catches a cut whose cap is only a fraction of the tolerance above the other
// carrier, which the tolerance band alone would treat as neutral.
inline constexpr int kSsiDiagnosticGlobalTestGrid = 21;
inline constexpr double kSsiDiagnosticGlobalTestFraction = 0.05;

// Carrier frame slack (1 - |unit dot unit|) and refined-contact normal slack.
// The normal check is strict: a crossing at 1e-9 rad must be rejected.
inline constexpr double kSsiDiagnosticFrameSlack = 1e-12;
inline constexpr double kSsiDiagnosticNormalSlack = 1e-10;

// Relative round-off floor of the geometric band: samples whose signed distance
// is within 1e-12 * (1 + max|p|) of the other carrier are treated as neutral,
// so double-precision noise at an exact contact is not read as a crossing side.
inline constexpr double kSsiDiagnosticSignSlack = 1e-12;

inline constexpr int kSsiDiagnosticGridCapacity =
    kSsiDiagnosticCoincidenceGrid * kSsiDiagnosticCoincidenceGrid;

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

// True when two unit directions are parallel or antiparallel within the frame
// slack. Sign folding is intentional: a carrier does not care whether its axis
// or normal points the other way.
inline bool ssi_diagnostic_unit_parallel(const Vec3d& a, const Vec3d& b) noexcept {
  const double alignment = std::fabs(dot(a, b));
  return std::isfinite(alignment) && 1.0 - alignment <= kSsiDiagnosticFrameSlack;
}

// Unit normal and signed offset of a plane, in the canonical form n . x = d.
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

inline bool ssi_diagnostic_patch_valid(const SsiPatch& patch) noexcept {
  if (!(std::isfinite(patch.u[0]) && std::isfinite(patch.u[1]) &&
        std::isfinite(patch.v[0]) && std::isfinite(patch.v[1]))) {
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

// Same infinite carrier within `tolerance`. The variants have already been
// checked to hold the same alternative, so get_if keeps the walk noexcept.
inline bool ssi_diagnostic_same_carrier(const SurfaceGeometry& a, const SurfaceGeometry& b,
                                        double tolerance) noexcept {
  if (const Plane3* plane_a = std::get_if<Plane3>(&a)) {
    const Plane3* plane_b = std::get_if<Plane3>(&b);
    if (plane_b == nullptr) return false;
    SsiDiagnosticPlane carrier_a{};
    SsiDiagnosticPlane carrier_b{};
    if (!ssi_diagnostic_plane(*plane_a, carrier_a) ||
        !ssi_diagnostic_plane(*plane_b, carrier_b)) {
      return false;
    }
    if (!ssi_diagnostic_unit_parallel(carrier_a.normal, carrier_b.normal)) return false;
    const double sign = dot(carrier_a.normal, carrier_b.normal) < 0.0 ? -1.0 : 1.0;
    return std::fabs(carrier_a.offset - sign * carrier_b.offset) <= tolerance;
  }
  if (const Sphere3* sphere_a = std::get_if<Sphere3>(&a)) {
    const Sphere3* sphere_b = std::get_if<Sphere3>(&b);
    if (sphere_b == nullptr) return false;
    // A sphere is fixed by its centre and radius; axis and x_axis are only a
    // pole/seam parameterization and must not participate.
    return norm2(sphere_a->center - sphere_b->center) <= tolerance &&
           std::fabs(sphere_a->radius - sphere_b->radius) <= tolerance;
  }
  if (const Cylinder3* cylinder_a = std::get_if<Cylinder3>(&a)) {
    const Cylinder3* cylinder_b = std::get_if<Cylinder3>(&b);
    if (cylinder_b == nullptr) return false;
    // Parallel (possibly opposite) axis directions, equal radius, and the same
    // axis line: the origin offset may slide along the axis but not across it.
    if (!ssi_diagnostic_unit_parallel(cylinder_a->axis, cylinder_b->axis)) return false;
    const Vec3d delta = cylinder_b->origin - cylinder_a->origin;
    const Vec3d radial = delta - cylinder_a->axis * dot(delta, cylinder_a->axis);
    return std::fabs(cylinder_a->radius - cylinder_b->radius) <= tolerance &&
           norm2(radial) <= tolerance;
  }
  if (const Cone3* cone_a = std::get_if<Cone3>(&a)) {
    const Cone3* cone_b = std::get_if<Cone3>(&b);
    if (cone_b == nullptr) return false;
    // A double cone is fixed by its apex, axis line, and half angle. The apex
    // is invariant under an axis flip to the other nappe. The half angle is an
    // angle, so its comparison is scaled by the carrier's local characteristic
    // radius (the base radius), not by its absolute position: a 1e-10 rad
    // difference is irrelevant on a millimetre-scale cone and significant on a
    // kilometre-scale one, while translating the same cone must not change the
    // verdict.
    if (!ssi_diagnostic_unit_parallel(cone_a->axis, cone_b->axis)) return false;
    const Point3d apex_a =
        cone_a->origin + cone_a->axis * (cone_a->radius / std::tan(cone_a->half_angle));
    const Point3d apex_b =
        cone_b->origin + cone_b->axis * (cone_b->radius / std::tan(cone_b->half_angle));
    if (norm2(apex_a - apex_b) > tolerance) return false;
    const double angular_scale = std::max(std::max(cone_a->radius, cone_b->radius), 1e-300);
    return std::fabs(cone_a->half_angle - cone_b->half_angle) * angular_scale <= tolerance;
  }
  if (const Torus3* torus_a = std::get_if<Torus3>(&a)) {
    const Torus3* torus_b = std::get_if<Torus3>(&b);
    if (torus_b == nullptr) return false;
    return norm2(torus_a->center - torus_b->center) <= tolerance &&
           ssi_diagnostic_unit_parallel(torus_a->axis, torus_b->axis) &&
           std::fabs(torus_a->major_radius - torus_b->major_radius) <= tolerance &&
           std::fabs(torus_a->minor_radius - torus_b->minor_radius) <= tolerance;
  }
  return false;  // NURBS: unsupported.
}

// Signed implicit distance to an analytic carrier. The sign is negative on the
// side that contains the carrier's interior (plane: opposite the normal;
// sphere/cylinder/torus: inside; cone: inside the double cone), and its
// magnitude is the exact distance to the carrier surface. NURBS returns
// +infinity so the callers reject it before sampling.
inline double ssi_diagnostic_signed_distance(const SurfaceGeometry& surface,
                                             const Point3d& point) noexcept {
  return std::visit(
      [&point](const auto& carrier) noexcept -> double {
        using Carrier = std::decay_t<decltype(carrier)>;
        if constexpr (std::is_same_v<Carrier, Plane3>) {
          SsiDiagnosticPlane plane{};
          if (!ssi_diagnostic_plane(carrier, plane)) return kSsiDiagnosticInfinity;
          return dot(ssi_diagnostic_vec(point), plane.normal) - plane.offset;
        } else if constexpr (std::is_same_v<Carrier, Sphere3>) {
          return norm2(point - carrier.center) - carrier.radius;
        } else if constexpr (std::is_same_v<Carrier, Cylinder3>) {
          const Vec3d delta = point - carrier.origin;
          const double axial = dot(delta, carrier.axis);
          return norm2(delta - carrier.axis * axial) - carrier.radius;
        } else if constexpr (std::is_same_v<Carrier, Cone3>) {
          const Vec3d delta = point - carrier.origin;
          const double axial = dot(delta, carrier.axis);
          const double radial = norm2(delta - carrier.axis * axial);
          const double from_apex = axial - carrier.radius / std::tan(carrier.half_angle);
          return radial * std::cos(carrier.half_angle) -
                 std::fabs(from_apex) * std::sin(carrier.half_angle);
        } else if constexpr (std::is_same_v<Carrier, Torus3>) {
          const Vec3d delta = point - carrier.center;
          const double axial = dot(delta, carrier.axis);
          const double radial = norm2(delta - carrier.axis * axial);
          return std::hypot(radial - carrier.major_radius, axial) - carrier.minor_radius;
        } else {
          return kSsiDiagnosticInfinity;
        }
      },
      surface);
}

// Samples a patch on a `grid` x `grid` uniform grid including both parameter
// boundaries. grid must not exceed kSsiDiagnosticCoincidenceGrid.
inline bool ssi_diagnostic_sample_patch(
    const SsiPatch& patch, int grid,
    std::array<SsiDiagnosticSample, kSsiDiagnosticGridCapacity>& samples) noexcept {
  const double denominator = static_cast<double>(grid - 1);
  for (int iu = 0; iu < grid; ++iu) {
    const double fraction_u = static_cast<double>(iu) / denominator;
    const double u = iu == 0 ? patch.u[0]
                             : (iu + 1 == grid ? patch.u[1]
                                               : patch.u[0] + (patch.u[1] - patch.u[0]) * fraction_u);
    for (int iv = 0; iv < grid; ++iv) {
      const double fraction_v = static_cast<double>(iv) / denominator;
      const double v = iv == 0 ? patch.v[0]
                               : (iv + 1 == grid ? patch.v[1]
                                                 : patch.v[0] + (patch.v[1] - patch.v[0]) * fraction_v);
      const Result<SurfacePoint> evaluated = surface_eval(patch.surface, u, v);
      if (!evaluated.has_value() || !surface_point_finite(evaluated.value())) return false;
      SsiDiagnosticSample& sample = samples[static_cast<std::size_t>(iu) * grid + iv];
      sample.point = evaluated.value().p;
      sample.uv = {u, v};
    }
  }
  return true;
}

// Wraps a periodic parameter back into the patch interval and clamps a
// non-periodic one, so the contact search never leaves the patch domain.
inline double ssi_diagnostic_place(double value, double lo, double hi,
                                   bool periodic) noexcept {
  if (!std::isfinite(value)) return value;
  if (periodic) {
    double folded = value - lo;
    folded -= kSsiTwoPi * std::floor(folded / kSsiTwoPi);
    value = lo + folded;
  }
  return std::clamp(value, lo, hi);
}

// Dense 4x4 solve with partial pivoting and a relative pivot guard. The
// normal equations are damped by the caller, so a singular pivot means the
// contact is not resolvable at this precision; false lets tangency answer
// false instead of stepping to garbage.
inline bool ssi_diagnostic_solve4(std::array<std::array<double, 4>, 4> matrix,
                                  std::array<double, 4> rhs,
                                  std::array<double, 4>& solution) noexcept {
  double scale = 0.0;
  for (int index = 0; index < 4; ++index) {
    scale = std::max(scale, std::fabs(matrix[index][index]));
  }
  if (!std::isfinite(scale) || !(scale > 0.0)) return false;
  const double threshold = kSsiDiagnosticPivotSlack * scale;
  for (int column = 0; column < 4; ++column) {
    int pivot = column;
    double best = std::fabs(matrix[column][column]);
    for (int row = column + 1; row < 4; ++row) {
      const double candidate = std::fabs(matrix[row][column]);
      if (candidate > best) {
        best = candidate;
        pivot = row;
      }
    }
    if (!(best > threshold)) return false;
    if (pivot != column) {
      std::swap(matrix[column], matrix[pivot]);
      std::swap(rhs[column], rhs[pivot]);
    }
    const double diagonal = matrix[column][column];
    for (int row = column + 1; row < 4; ++row) {
      const double factor = matrix[row][column] / diagonal;
      for (int index = column; index < 4; ++index) {
        matrix[row][index] -= factor * matrix[column][index];
      }
      rhs[row] -= factor * rhs[column];
    }
  }
  for (int row = 3; row >= 0; --row) {
    double value = rhs[row];
    for (int index = row + 1; index < 4; ++index) {
      value -= matrix[row][index] * solution[index];
    }
    const double diagonal = matrix[row][row];
    if (!(std::fabs(diagonal) > threshold)) return false;
    solution[row] = value / diagonal;
  }
  for (int index = 0; index < 4; ++index) {
    if (!std::isfinite(solution[index])) return false;
  }
  return true;
}

struct SsiDiagnosticContact {
  std::array<double, 2> uv_a{};
  std::array<double, 2> uv_b{};
  Point3d point_a{};
  Point3d point_b{};
  double gap = kSsiDiagnosticInfinity;
};

// Evaluates both patches at `parameters` and stores the contact. Returns false
// only for an unusable (missing or non-finite) evaluation.
inline bool ssi_diagnostic_store_contact(const SsiPatch& a, const SsiPatch& b,
                                         const std::array<double, 4>& parameters,
                                         SsiDiagnosticContact& contact) noexcept {
  const Result<SurfacePoint> eval_a = surface_eval(a.surface, parameters[0], parameters[1]);
  const Result<SurfacePoint> eval_b = surface_eval(b.surface, parameters[2], parameters[3]);
  if (!eval_a.has_value() || !eval_b.has_value()) return false;
  if (!surface_point_finite(eval_a.value()) || !surface_point_finite(eval_b.value())) {
    return false;
  }
  contact.uv_a = {parameters[0], parameters[1]};
  contact.uv_b = {parameters[2], parameters[3]};
  contact.point_a = eval_a.value().p;
  contact.point_b = eval_b.value().p;
  contact.gap = norm2(eval_a.value().p - eval_b.value().p);
  return std::isfinite(contact.gap);
}

// Damped Gauss-Newton on the four patch parameters that minimizes |A - B|^2.
// Line contacts and sphere poles leave a null direction in JtJ, so the step is
// damped; a step is accepted only when it strictly shrinks the gap, and the
// iteration stops as soon as no halved step does.
inline bool ssi_diagnostic_gauss_newton(const SsiPatch& a, const SsiPatch& b,
                                        std::array<double, 4>& parameters,
                                        int iterations) noexcept {
  const std::array<std::array<double, 2>, 4> ranges{a.u, a.v, b.u, b.v};
  const std::array<bool, 4> periodic{surface_periodic(a.surface, 0),
                                     surface_periodic(a.surface, 1),
                                     surface_periodic(b.surface, 0),
                                     surface_periodic(b.surface, 1)};
  for (int iteration = 0; iteration < iterations; ++iteration) {
    const Result<SurfacePoint> eval_a = surface_eval(a.surface, parameters[0], parameters[1]);
    const Result<SurfacePoint> eval_b = surface_eval(b.surface, parameters[2], parameters[3]);
    if (!eval_a.has_value() || !eval_b.has_value()) return false;
    if (!surface_point_finite(eval_a.value()) || !surface_point_finite(eval_b.value())) {
      return false;
    }
    const Vec3d residual = eval_a.value().p - eval_b.value().p;
    const double gap = norm2(residual);
    if (!std::isfinite(gap)) return false;
    const std::array<Vec3d, 4> columns{eval_a.value().du, eval_a.value().dv,
                                       eval_b.value().du * -1.0,
                                       eval_b.value().dv * -1.0};
    std::array<std::array<double, 4>, 4> normal{};
    std::array<double, 4> gradient{};
    for (int row = 0; row < 4; ++row) {
      for (int column = 0; column < 4; ++column) {
        normal[row][column] = dot(columns[row], columns[column]);
      }
      gradient[row] = dot(columns[row], residual);
    }
    const double scale = (normal[0][0] + normal[1][1] + normal[2][2] + normal[3][3]) / 4.0;
    if (!std::isfinite(scale) || !(scale > 0.0)) return false;
    const double damping = kSsiDiagnosticDamping * scale;
    std::array<double, 4> rhs{};
    for (int index = 0; index < 4; ++index) {
      normal[index][index] += damping;
      rhs[index] = -gradient[index];
    }
    std::array<double, 4> step{};
    if (!ssi_diagnostic_solve4(normal, rhs, step)) return false;

    bool accepted = false;
    double factor = 1.0;
    for (int halving = 0; halving < 24; ++halving) {
      std::array<double, 4> candidate{};
      for (int index = 0; index < 4; ++index) {
        candidate[index] = ssi_diagnostic_place(parameters[index] + factor * step[index],
                                                ranges[index][0], ranges[index][1],
                                                periodic[index]);
      }
      const Result<SurfacePoint> candidate_a =
          surface_eval(a.surface, candidate[0], candidate[1]);
      const Result<SurfacePoint> candidate_b =
          surface_eval(b.surface, candidate[2], candidate[3]);
      if (candidate_a.has_value() && candidate_b.has_value() &&
          surface_point_finite(candidate_a.value()) &&
          surface_point_finite(candidate_b.value())) {
        const double candidate_gap =
            norm2(candidate_a.value().p - candidate_b.value().p);
        if (std::isfinite(candidate_gap) && candidate_gap < gap) {
          parameters = candidate;
          accepted = true;
          break;
        }
      }
      factor *= 0.5;
    }
    if (!accepted) break;  // The contact is stationary at this precision.
  }
  return true;
}

// Unit normal at a patch parameter pair; false when the evaluation is unusable.
inline bool ssi_diagnostic_normal_at(const SurfaceGeometry& surface,
                                     const std::array<double, 2>& uv,
                                     Vec3d& normal) noexcept {
  const Result<Vec3d> result = surface_normal(surface, uv[0], uv[1]);
  if (!result.has_value()) return false;
  if (!ssi_diagnostic_vec_finite(result.value())) return false;
  normal = result.value();
  return true;
}

// Moves a contact that already touches within tolerance onto the tangency set.
// The gap refinement of a line or circle contact stalls at a tangential drift
// whose size is the normal tilt, so this pass alternates a damped Gauss-Newton
// step on the normal mismatch m = cross(nA, nB) (finite-difference normal
// Jacobian) with ordinary contact Gauss-Newton that re-closes the gap. The gap
// step splits the correction between the patches, so the mismatch halves per
// round instead of disappearing; 24 rounds take a 1e-6 mismatch below 1e-9.
// Returns true only when the final unit normals are parallel; a pair with
// constant normals (two crossing planes) has an all-zero mismatch Jacobian and
// returns false immediately.
inline bool ssi_diagnostic_align_contact(const SsiPatch& a, const SsiPatch& b,
                                         SsiDiagnosticContact& contact) noexcept {
  std::array<double, 4> parameters{contact.uv_a[0], contact.uv_a[1], contact.uv_b[0],
                                   contact.uv_b[1]};
  const std::array<std::array<double, 2>, 4> ranges{a.u, a.v, b.u, b.v};
  const std::array<bool, 4> periodic{surface_periodic(a.surface, 0),
                                     surface_periodic(a.surface, 1),
                                     surface_periodic(b.surface, 0),
                                     surface_periodic(b.surface, 1)};
  std::array<double, 2> uv_a{parameters[0], parameters[1]};
  std::array<double, 2> uv_b{parameters[2], parameters[3]};
  Vec3d normal_a{};
  Vec3d normal_b{};
  if (!ssi_diagnostic_normal_at(a.surface, uv_a, normal_a) ||
      !ssi_diagnostic_normal_at(b.surface, uv_b, normal_b)) {
    return false;
  }
  for (int round = 0; round < kSsiDiagnosticAlignRounds; ++round) {
    const Vec3d mismatch = cross(normal_a, normal_b);
    if (!ssi_diagnostic_vec_finite(mismatch)) return false;
    const double mismatch_norm = norm2(mismatch);
    if (mismatch_norm < kSsiDiagnosticNormalSlack) {
      return ssi_diagnostic_store_contact(a, b, parameters, contact);
    }
    std::array<Vec3d, 4> columns{};
    for (int index = 0; index < 4; ++index) {
      std::array<double, 2> plus = index < 2 ? uv_a : uv_b;
      std::array<double, 2> minus = plus;
      const int axis = index % 2;
      plus[axis] += kSsiDiagnosticAlignStep;
      minus[axis] -= kSsiDiagnosticAlignStep;
      const SurfaceGeometry& surface = index < 2 ? a.surface : b.surface;
      Vec3d normal_plus{};
      Vec3d normal_minus{};
      if (!ssi_diagnostic_normal_at(surface, plus, normal_plus) ||
          !ssi_diagnostic_normal_at(surface, minus, normal_minus)) {
        return false;
      }
      const Vec3d derivative =
          (normal_plus - normal_minus) * (0.5 / kSsiDiagnosticAlignStep);
      columns[index] = index < 2 ? cross(derivative, normal_b) : cross(normal_a, derivative);
    }
    std::array<std::array<double, 4>, 4> normal{};
    std::array<double, 4> rhs{};
    for (int row = 0; row < 4; ++row) {
      for (int column = 0; column < 4; ++column) {
        normal[row][column] = dot(columns[row], columns[column]);
      }
      rhs[row] = -dot(columns[row], mismatch);
    }
    const double scale = (normal[0][0] + normal[1][1] + normal[2][2] + normal[3][3]) / 4.0;
    if (!std::isfinite(scale) || !(scale > 0.0)) return false;
    const double damping = kSsiDiagnosticDamping * scale;
    for (int index = 0; index < 4; ++index) normal[index][index] += damping;
    std::array<double, 4> step{};
    if (!ssi_diagnostic_solve4(normal, rhs, step)) return false;

    bool accepted = false;
    double factor = 1.0;
    for (int halving = 0; halving < 24; ++halving) {
      std::array<double, 4> candidate{};
      for (int index = 0; index < 4; ++index) {
        candidate[index] = ssi_diagnostic_place(parameters[index] + factor * step[index],
                                                ranges[index][0], ranges[index][1],
                                                periodic[index]);
      }
      const std::array<double, 2> candidate_uv_a{candidate[0], candidate[1]};
      const std::array<double, 2> candidate_uv_b{candidate[2], candidate[3]};
      Vec3d candidate_normal_a{};
      Vec3d candidate_normal_b{};
      if (ssi_diagnostic_normal_at(a.surface, candidate_uv_a, candidate_normal_a) &&
          ssi_diagnostic_normal_at(b.surface, candidate_uv_b, candidate_normal_b)) {
        const Vec3d candidate_mismatch = cross(candidate_normal_a, candidate_normal_b);
        if (ssi_diagnostic_vec_finite(candidate_mismatch) &&
            norm2(candidate_mismatch) < mismatch_norm) {
          parameters = candidate;
          uv_a = candidate_uv_a;
          uv_b = candidate_uv_b;
          normal_a = candidate_normal_a;
          normal_b = candidate_normal_b;
          accepted = true;
          break;
        }
      }
      factor *= 0.5;
    }
    if (!accepted) return false;
    if (!ssi_diagnostic_gauss_newton(a, b, parameters,
                                     kSsiDiagnosticAlignGaussNewtonIterations)) {
      return false;
    }
    uv_a = {parameters[0], parameters[1]};
    uv_b = {parameters[2], parameters[3]};
    if (!ssi_diagnostic_normal_at(a.surface, uv_a, normal_a) ||
        !ssi_diagnostic_normal_at(b.surface, uv_b, normal_b)) {
      return false;
    }
  }
  const Vec3d remainder = cross(normal_a, normal_b);
  if (!ssi_diagnostic_vec_finite(remainder)) return false;
  if (!(norm2(remainder) < kSsiDiagnosticNormalSlack)) return false;
  return ssi_diagnostic_store_contact(a, b, parameters, contact);
}

// Coarse 9x9 x 9x9 search followed by damped Gauss-Newton on the four patch
// parameters. Returns false only when an evaluation, the linear solve, or the
// final gap is unusable (non-finite); a contact that simply does not exist is
// reported through contact.gap.
inline bool ssi_diagnostic_refine_contact(const SsiPatch& a, const SsiPatch& b,
                                          SsiDiagnosticContact& contact) noexcept {
  std::array<SsiDiagnosticSample, kSsiDiagnosticGridCapacity> samples_a{};
  std::array<SsiDiagnosticSample, kSsiDiagnosticGridCapacity> samples_b{};
  if (!ssi_diagnostic_sample_patch(a, kSsiDiagnosticContactGrid, samples_a) ||
      !ssi_diagnostic_sample_patch(b, kSsiDiagnosticContactGrid, samples_b)) {
    return false;
  }
  double best_gap = kSsiDiagnosticInfinity;
  std::array<double, 2> best_uv_a{};
  std::array<double, 2> best_uv_b{};
  const std::size_t sample_count = static_cast<std::size_t>(kSsiDiagnosticContactGrid) *
                                   static_cast<std::size_t>(kSsiDiagnosticContactGrid);
  for (std::size_t ia = 0; ia < sample_count; ++ia) {
    for (std::size_t ib = 0; ib < sample_count; ++ib) {
      const double gap = norm2(samples_a[ia].point - samples_b[ib].point);
      if (!std::isfinite(gap)) return false;
      if (gap < best_gap) {
        best_gap = gap;
        best_uv_a = samples_a[ia].uv;
        best_uv_b = samples_b[ib].uv;
      }
    }
  }

  std::array<double, 4> parameters{best_uv_a[0], best_uv_a[1], best_uv_b[0], best_uv_b[1]};
  if (!ssi_diagnostic_gauss_newton(a, b, parameters, kSsiDiagnosticRefineIterations)) {
    return false;
  }
  return ssi_diagnostic_store_contact(a, b, parameters, contact);
}

// Global one-sided test: the whole domain of `patch` is sampled on a fixed
// 21x21 parameter grid (endpoints included; periodic axes wrap by 2*pi) and
// every sample's signed distance to `other`'s carrier must stay on one side.
// Two bands are checked: the tolerance band and the geometric band
// max(1e-12 * (1 + max|p|), 0.05 * tolerance). Either band showing samples
// beyond both +band and -band rejects the pair. A crossing or cutting carrier
// produces both signs; a touching carrier produces one. A cut shallower than
// the geometric band stays inside the detection floor and is not visible here
// (see the banner comment). The grid is fixed and stack-allocated, so the test
// stays noexcept and allocation-free.
inline bool ssi_diagnostic_sign_test(const SsiPatch& patch, const SurfaceGeometry& other,
                                     double tolerance) noexcept {
  const double denominator = static_cast<double>(kSsiDiagnosticGlobalTestGrid - 1);
  const double geometric_fraction = kSsiDiagnosticGlobalTestFraction * tolerance;
  bool tolerance_positive = false;
  bool tolerance_negative = false;
  bool geometric_positive = false;
  bool geometric_negative = false;
  for (int iu = 0; iu < kSsiDiagnosticGlobalTestGrid; ++iu) {
    const double u = ssi_diagnostic_place(
        patch.u[0] + (patch.u[1] - patch.u[0]) * (static_cast<double>(iu) / denominator),
        patch.u[0], patch.u[1], surface_periodic(patch.surface, 0));
    for (int iv = 0; iv < kSsiDiagnosticGlobalTestGrid; ++iv) {
      const double v = ssi_diagnostic_place(
          patch.v[0] + (patch.v[1] - patch.v[0]) * (static_cast<double>(iv) / denominator),
          patch.v[0], patch.v[1], surface_periodic(patch.surface, 1));
      const Result<SurfacePoint> evaluated = surface_eval(patch.surface, u, v);
      if (!evaluated.has_value() || !surface_point_finite(evaluated.value())) return false;
      const Point3d& point = evaluated.value().p;
      const double signed_distance = ssi_diagnostic_signed_distance(other, point);
      if (!std::isfinite(signed_distance)) return false;
      const double magnitude = std::max(std::max(std::fabs(point[0]), std::fabs(point[1])),
                                        std::fabs(point[2]));
      const double geometric_band =
          std::max(kSsiDiagnosticSignSlack * (1.0 + magnitude), geometric_fraction);
      if (signed_distance > tolerance) tolerance_positive = true;
      if (signed_distance < -tolerance) tolerance_negative = true;
      if (signed_distance > geometric_band) geometric_positive = true;
      if (signed_distance < -geometric_band) geometric_negative = true;
      if ((tolerance_positive && tolerance_negative) ||
          (geometric_positive && geometric_negative)) {
        return false;
      }
    }
  }
  return true;
}

}  // namespace detail

// True when both patches lie on the same infinite carrier within `tolerance`.
// The carrier parameters are compared exactly (see the banner comment) and the
// domain of patch a is certified on the 9x9 sample grid described above;
// patch b's domain is covered by the same carrier, so it needs no separate
// sampling. Non-finite coordinates, non-positive or non-finite tolerances,
// mixed analytic types, and NURBS patches answer false.
inline bool surfaces_coincident(const SsiPatch& a, const SsiPatch& b,
                                double tolerance) noexcept {
  if (!std::isfinite(tolerance) || !(tolerance > 0.0)) return false;
  if (a.surface.index() != b.surface.index()) return false;
  if (std::holds_alternative<NurbsSurface>(a.surface)) return false;
  if (!detail::ssi_diagnostic_patch_valid(a) || !detail::ssi_diagnostic_patch_valid(b)) {
    return false;
  }
  if (!detail::ssi_diagnostic_same_carrier(a.surface, b.surface, tolerance)) return false;

  std::array<detail::SsiDiagnosticSample, detail::kSsiDiagnosticGridCapacity> samples{};
  if (!detail::ssi_diagnostic_sample_patch(a, detail::kSsiDiagnosticCoincidenceGrid, samples)) {
    return false;
  }
  const std::size_t sample_count =
      static_cast<std::size_t>(detail::kSsiDiagnosticCoincidenceGrid) *
      static_cast<std::size_t>(detail::kSsiDiagnosticCoincidenceGrid);
  for (std::size_t index = 0; index < sample_count; ++index) {
    const double signed_distance =
        detail::ssi_diagnostic_signed_distance(b.surface, samples[index].point);
    if (!std::isfinite(signed_distance) || std::fabs(signed_distance) > tolerance) {
      return false;
    }
  }
  return true;
}

// True when the two carriers touch within `tolerance` at a refined contact with
// parallel normals and do not coincide. The refined contact must be inside both
// patches, and the global one-sided test described in the banner comment
// rejects carriers that cross, cut, or interpenetrate. Cuts shallower than the
// documented detection floor max(1e-12 * (1 + max|p|), 0.05 * tolerance) are
// reported as tangent_contact by contract. NURBS patches answer false because
// their tangency cannot be certified without an exact projection.
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

  detail::SsiDiagnosticContact contact;
  if (!detail::ssi_diagnostic_refine_contact(a, b, contact)) return false;
  if (!(contact.gap <= tolerance)) return false;
  if (!detail::ssi_diagnostic_align_contact(a, b, contact)) return false;
  if (!(contact.gap <= tolerance)) return false;
  if (!patch_contains(a, contact.uv_a[0], contact.uv_a[1]) ||
      !patch_contains(b, contact.uv_b[0], contact.uv_b[1])) {
    return false;
  }
  const Result<Vec3d> normal_a = surface_normal(a.surface, contact.uv_a[0], contact.uv_a[1]);
  const Result<Vec3d> normal_b = surface_normal(b.surface, contact.uv_b[0], contact.uv_b[1]);
  if (!normal_a.has_value() || !normal_b.has_value()) return false;
  const Vec3d normal_cross = cross(normal_a.value(), normal_b.value());
  if (!detail::ssi_diagnostic_vec_finite(normal_cross)) return false;
  if (!(norm2(normal_cross) < detail::kSsiDiagnosticNormalSlack)) return false;
  if (!detail::ssi_diagnostic_sign_test(a, b.surface, tolerance)) return false;
  if (!detail::ssi_diagnostic_sign_test(b, a.surface, tolerance)) return false;
  return true;
}

// ---------------------------------------------------------------------------
// Driver helpers
//
// These support the ssi() orchestration in brep/ssi.hpp: they turn
// triangle-triangle endpoints into patch parameters, refine them into
// seeds, and group the seeds into one traceable seed per intersection
// curve. They live here so the public driver header stays a thin
// vocabulary plus orchestration layer.
// ---------------------------------------------------------------------------

namespace detail {

// Barycentric coordinates of p projected onto triangle (q0, q1, q2), clamped to
// the closed triangle. A triangle-triangle endpoint lies on the triangle
// boundary, so the clamp is what keeps the interpolated parameters inside the
// patch domain instead of a few ulps outside it. False means the triangle is too
// thin to interpolate safely or a weight is not finite.
inline bool ssi_driver_barycentric(const Point3d& p, const Point3d& q0, const Point3d& q1,
                                   const Point3d& q2, std::array<double, 3>& weights) noexcept {
  const Vec3d e0 = q1 - q0;
  const Vec3d e1 = q2 - q0;
  const Vec3d offset = p - q0;
  const double d00 = dot(e0, e0);
  const double d01 = dot(e0, e1);
  const double d11 = dot(e1, e1);
  const double det = d00 * d11 - d01 * d01;
  if (!std::isfinite(det) || !(det > 1e-20 * d00 * d11)) return false;
  const double d20 = dot(offset, e0);
  const double d21 = dot(offset, e1);
  const double raw1 = (d11 * d20 - d01 * d21) / det;
  const double raw2 = (d00 * d21 - d01 * d20) / det;
  const double raw0 = 1.0 - raw1 - raw2;
  if (!std::isfinite(raw0) || !std::isfinite(raw1) || !std::isfinite(raw2)) return false;
  const double clamped0 = std::min(1.0, std::max(0.0, raw0));
  const double clamped1 = std::min(1.0, std::max(0.0, raw1));
  const double clamped2 = std::min(1.0, std::max(0.0, raw2));
  const double sum = clamped0 + clamped1 + clamped2;
  if (!(sum > 0.0)) return false;
  weights = {clamped0 / sum, clamped1 / sum, clamped2 / sum};
  return true;
}

// Interpolates the (u, v) parameters stored for one triangle at a 3D point.
inline bool ssi_driver_triangle_uv(const SsiPatchMesh& mesh,
                                   const std::array<std::uint32_t, 3>& triangle, const Point3d& p,
                                   std::array<double, 2>& uv) noexcept {
  std::array<double, 3> weights{};
  if (!ssi_driver_barycentric(p, mesh.points[triangle[0]], mesh.points[triangle[1]],
                              mesh.points[triangle[2]], weights)) {
    return false;
  }
  const std::array<double, 2>& uv0 = mesh.uv[triangle[0]];
  const std::array<double, 2>& uv1 = mesh.uv[triangle[1]];
  const std::array<double, 2>& uv2 = mesh.uv[triangle[2]];
  uv[0] = weights[0] * uv0[0] + weights[1] * uv1[0] + weights[2] * uv2[0];
  uv[1] = weights[0] * uv0[1] + weights[1] * uv1[1] + weights[2] * uv2[1];
  return std::isfinite(uv[0]) && std::isfinite(uv[1]);
}

inline Box3d ssi_driver_triangle_bounds(const SsiPatchMesh& mesh,
                                        const std::array<std::uint32_t, 3>& triangle) noexcept {
  Box3d bounds = Box3d::empty();
  bounds.expand(mesh.points[triangle[0]]);
  bounds.expand(mesh.points[triangle[1]]);
  bounds.expand(mesh.points[triangle[2]]);
  return bounds;
}

inline bool ssi_driver_bounds_overlap(const Box3d& first, const Box3d& second) noexcept {
  for (int axis = 0; axis < 3; ++axis) {
    if (first.lo[axis] > second.hi[axis] || second.lo[axis] > first.hi[axis]) return false;
  }
  return true;
}

// Total order on refined seeds. The parameter pairs dominate so clustering and
// output never depend on mesh traversal order; the remaining keys only break
// exact ties deterministically.
inline bool ssi_driver_seed_less(const SsiPoint& left, const SsiPoint& right) noexcept {
  if (left.ua < right.ua) return true;
  if (right.ua < left.ua) return false;
  if (left.va < right.va) return true;
  if (right.va < left.va) return false;
  if (left.residual < right.residual) return true;
  if (right.residual < left.residual) return false;
  if (left.p[0] != right.p[0]) return left.p[0] < right.p[0];
  if (left.p[1] != right.p[1]) return left.p[1] < right.p[1];
  return left.p[2] < right.p[2];
}

inline std::size_t ssi_driver_seed_root(std::vector<std::size_t>& parent,
                                        std::size_t index) noexcept {
  std::size_t root = index;
  while (parent[root] != root) root = parent[root];
  while (parent[index] != root) {
    const std::size_t next = parent[index];
    parent[index] = root;
    index = next;
  }
  return root;
}

// Stands for a triangle-triangle endpoint that did not survive Newton
// refinement, so a half-refined segment contributes only the endpoint it has.
inline constexpr std::size_t kSsiDriverNoSeed = static_cast<std::size_t>(-1);

// Collapses the refined triangle-triangle endpoints into one seed per traced
// curve. Two seeds name the same intersection point when their 3D positions
// agree within seed_tolerance, and the two endpoints of one segment name two
// points of the same curve. Both relations are unioned, so a connected
// component is one intersection curve: adjacent mesh triangles share endpoints
// and the chain of segments is exactly the curve. One deterministic
// representative (the smallest parameter pair) is returned per component;
// without this, tracing every distinct point of a curve would retrace the same
// curve once per point.
//
// The comparison is in 3D on purpose. Parameter units differ between surfaces
// (a plane's parameters are lengths, a sphere's or torus's are angles), and
// Newton paths from neighbouring triangle pairs leave the shared hit point a
// small but dimension-dependent amount apart; a parameter-space radius would
// grow with the plane extent until it merged two genuinely distinct torus
// branches (~0.6 apart at z=0.4) into one curve.
inline std::vector<SsiPoint> ssi_driver_curve_seeds(
    const SsiPatch& a, const SsiPatch& b, const std::vector<SsiPoint>& seeds,
    const std::vector<std::array<std::size_t, 2>>& segments, double seed_tolerance,
    double max_chord_error_a = 0.0, double max_chord_error_b = 0.0) {
  // The patch pair is kept in the interface for symmetry with the tracing
  // helpers; the union itself only compares the refined 3D endpoints.
  (void)a;
  (void)b;
  const std::size_t count = seeds.size();
  std::vector<std::size_t> parent(count);
  for (std::size_t index = 0; index < count; ++index) parent[index] = index;

  // The tessellation, not the solver tolerance, bounds how far apart two
  // endpoints that name the same intersection point can be: adaptive leaves
  // create hanging nodes, so segments from neighbouring triangle pairs meet
  // about one chord error apart. The cluster radius has to cover that error or
  // one loop is reported as several arcs. The cost is that two genuinely
  // distinct branches closer than 2*chord_error merge; that near-tangency case
  // is diagnosed by the Task 6 predicates and recorded as an M3b-1 limitation.
  seed_tolerance =
      std::max(seed_tolerance, 1.5 * std::max(max_chord_error_a, max_chord_error_b));

  // Near a pole the parameter slide can spread duplicate seeds a few times
  // farther than one chord error even though their 3D points still coincide,
  // so the radius also gets a small fraction of the smaller patch diagonal.
  // It stays well below the spacing of distinct branches (about 0.6 on the
  // test torus).
  const Box3d bounds_a = patch_bounds(a, seed_tolerance);
  const Box3d bounds_b = patch_bounds(b, seed_tolerance);
  if (!bounds_a.is_empty() && !bounds_b.is_empty()) {
    const double diagonal_a = norm2(bounds_a.hi - bounds_a.lo);
    const double diagonal_b = norm2(bounds_b.hi - bounds_b.lo);
    seed_tolerance = std::max(seed_tolerance, 0.01 * std::min(diagonal_a, diagonal_b));
  }

  // A widened radius can also reach across to a different branch that passes
  // close by, so a candidate merge must agree in direction as well: adjacent
  // mesh segments continue the same curve and are nearly collinear, while two
  // crossing branches leave their common neighbourhood at an angle.
  std::vector<std::array<double, 3>> seed_tangent(count);
  std::vector<bool> seed_tangent_valid(count, false);
  for (const std::array<std::size_t, 2>& segment : segments) {
    if (segment[0] >= count || segment[1] >= count || segment[0] == segment[1]) continue;
    const Vec3d direction = seeds[segment[1]].p - seeds[segment[0]].p;
    const double length = norm2(direction);
    if (!(length > 0.0)) continue;
    const std::array<double, 3> unit{direction[0] / length, direction[1] / length,
                                     direction[2] / length};
    seed_tangent[segment[0]] = unit;
    seed_tangent_valid[segment[0]] = true;
    seed_tangent[segment[1]] = unit;
    seed_tangent_valid[segment[1]] = true;
  }
  const auto same_direction = [&](std::size_t first, std::size_t second) {
    if (!seed_tangent_valid[first] || !seed_tangent_valid[second]) return true;
    double dot = 0.0;
    for (int axis = 0; axis < 3; ++axis) {
      dot += seed_tangent[first][static_cast<std::size_t>(axis)] *
             seed_tangent[second][static_cast<std::size_t>(axis)];
    }
    return std::fabs(dot) >= 0.9;
  };

  // Duplicates of one intersection point sit at the same 3D location, so a hash
  // of the refined point finds them without a quadratic scan; the neighbourhood
  // probe keeps a pair from being missed when it straddles a cell wall. The 3D
  // distance below is the actual clustering rule.
  std::map<std::array<long long, 3>, std::vector<std::size_t>> buckets;
  const auto cell_of = [seed_tolerance](const Point3d& point) {
    std::array<long long, 3> cell{};
    for (int axis = 0; axis < 3; ++axis) {
      const double scaled = point[axis] / seed_tolerance;
      const double finite = std::isfinite(scaled) ? scaled : 0.0;
      cell[static_cast<std::size_t>(axis)] = static_cast<long long>(std::floor(finite));
    }
    return cell;
  };
  for (std::size_t index = 0; index < count; ++index) {
    buckets[cell_of(seeds[index].p)].push_back(index);
  }
  for (std::size_t index = 0; index < count; ++index) {
    const std::array<long long, 3> cell = cell_of(seeds[index].p);
    for (int dx = -2; dx <= 2; ++dx) {
      for (int dy = -2; dy <= 2; ++dy) {
        for (int dz = -2; dz <= 2; ++dz) {
          const auto found = buckets.find(
              std::array<long long, 3>{cell[0] + dx, cell[1] + dy, cell[2] + dz});
          if (found == buckets.end()) continue;
          for (const std::size_t other : found->second) {
            if (other <= index) continue;
            if (!(norm2(seeds[index].p - seeds[other].p) <= seed_tolerance)) continue;
            if (!same_direction(index, other)) continue;
            const std::size_t root_first = ssi_driver_seed_root(parent, index);
            const std::size_t root_second = ssi_driver_seed_root(parent, other);
            if (root_first != root_second) parent[root_second] = root_first;
          }
        }
      }
    }
  }
  for (const std::array<std::size_t, 2>& segment : segments) {
    if (segment[0] >= count || segment[1] >= count) continue;
    const std::size_t root_first = ssi_driver_seed_root(parent, segment[0]);
    const std::size_t root_second = ssi_driver_seed_root(parent, segment[1]);
    if (root_first != root_second) parent[root_second] = root_first;
  }
  std::vector<SsiPoint> representatives;
  std::vector<std::size_t> slot(count, kSsiDriverNoSeed);
  for (std::size_t index = 0; index < count; ++index) {
    const std::size_t root = ssi_driver_seed_root(parent, index);
    if (slot[root] == kSsiDriverNoSeed) {
      slot[root] = representatives.size();
      representatives.push_back(seeds[index]);
    } else if (ssi_driver_seed_less(seeds[index], representatives[slot[root]])) {
      representatives[slot[root]] = seeds[index];
    }
  }
  std::stable_sort(representatives.begin(), representatives.end(), ssi_driver_seed_less);
  return representatives;
}

}  // namespace detail

}  // namespace v1
}  // namespace geom
