#pragma once

// Public surface-surface intersection diagnostics and driver.
//
// ssi() is the header-only intersection driver: it classifies degenerate,
// coincident, and tangent input, tessellates both patches, turns the endpoints
// of the triangle-triangle intersection segments into Newton seeds, and traces
// the curves. The predicates and the point/curve tracing helpers it builds on
// are inline functions in detail/ssi_solver.hpp; this header pins the shared
// vocabulary and the orchestration so every translation unit sees the same
// definition.

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstring>
#include <limits>
#include <utility>
#include <vector>

#include <geom/brep/detail/ssi_patch.hpp>
#include <geom/brep/detail/ssi_solver.hpp>
#include <geom/brep/detail/ssi_tri_tri.hpp>
#include <geom/core/result.hpp>

namespace geom {
inline namespace v1 {

// Why an SSI request returned no usable curve, or why its input was rejected
// before tracing. none is the zero value, so an uninitialized issue can never
// masquerade as a classified diagnostic.
enum class SsiDiagnostic : std::uint8_t {
  none = 0,
  tangent_contact,
  coincident,
  degenerate_patch,
  no_convergence,
};

// One diagnostic attached to an SsiResult. params carries the four patch
// parameters most relevant to the issue, in (u_a, v_a, u_b, v_b) order, or
// zeros when the issue is purely structural. message is a static, non-owning
// explanation owned by the implementation.
struct SsiIssue {
  SsiDiagnostic code = SsiDiagnostic::none;
  std::array<double, 4> params{};
  const char* message = "";
};

// Result of one SSI request. curves are the traced intersection curves (empty
// when the request was classified as degenerate before tracing). The counters
// describe the work done and are filled by the driver; the defaults are the
// zero-work state.
struct SsiResult {
  std::vector<SsiCurve> curves;
  std::vector<SsiIssue> issues;
  // Triangle pairs whose axis-aligned boxes overlapped in the broad phase.
  std::uint64_t candidate_pairs = 0;
  // Segments in the returned curves.
  std::uint64_t accepted_segments = 0;
  // Extracted triangle-pair segments that the returned curves do not account
  // for. When the trace resolves a curve more finely than the tessellated seed
  // set, several trace segments can share one extracted segment; the difference
  // then saturates at zero rather than going negative.
  std::uint64_t rejected_segments = 0;
  // Largest residual over the returned segments; zero when nothing was traced.
  double max_residual = 0.0;
};

// Budget for one SSI request. tolerance is the geometric acceptance band,
// max_depth bounds the adaptive patch subdivision, and max_segments bounds the
// traced output. An invalid tolerance or a non-positive budget is rejected by
// ssi() rather than silently clamped. max_depth must be at least 6: the seed
// clustering assumes at least 64 cells per patch side, and coarser meshes
// produce duplicate curves or order-dependent closure failures.
struct SsiOptions {
  double tolerance = 1e-9;
  int max_depth = 12;
  std::uint64_t max_segments = 100000;
};

namespace detail {

// Builds one structural diagnostic. params stays zero for issues that do not
// name a parameter location.
inline SsiIssue ssi_driver_issue(SsiDiagnostic code, const char* message) noexcept {
  SsiIssue issue;
  issue.code = code;
  issue.message = message;
  return issue;
}

}  // namespace detail

// Traces the intersection of two patches. Degenerate, coincident, and tangent
// input is classified before the tessellator runs; otherwise both patches are
// tessellated adaptively, every AABB-overlapping triangle pair is intersected,
// and the segment endpoints become Newton seeds. Seeds that share a refined
// endpoint belong to the same intersection curve, so one seed per connected
// segment component is traced: one curve yields one traced curve instead of one
// per mesh segment. The counters record the broad phase (candidate pairs), the
// traced output (accepted segments and the worst residual), and the extracted
// triangle-pair segments that the traced curves did not consume (rejected).
// Curves are ordered by their first point and diagnostics by code and
// parameters, so two calls on identical input produce bit-identical output.
inline Result<SsiResult> ssi(const SsiPatch& a, const SsiPatch& b, const SsiOptions& options) {
  if (!std::isfinite(options.tolerance) || !(options.tolerance > 0.0)) {
    return Error{Errc::numeric_failure, "ssi: tolerance must be positive and finite"};
  }
  // The seed clustering assumes at least 64 cells per patch side: coarser
  // meshes make neighbouring triangle pair seeds land far enough apart that one
  // loop splits into duplicate, argument-order-dependent curves. The floor is
  // part of this task's contract.
  if (options.max_depth < 6) {
    return Error{Errc::unsupported, "ssi: max_depth must be at least 6"};
  }
  if (options.max_segments < 1) {
    return Error{Errc::numeric_failure, "ssi: max_segments must be positive"};
  }
  SsiResult result;
  const double tolerance = options.tolerance;

  // A patch whose 3D extent is inside the acceptance band has no well-defined
  // parameterization to trace on, so it is diagnosed before any predicate or
  // tessellation work.
  const Box3d bounds_a = patch_bounds(a, tolerance);
  const Box3d bounds_b = patch_bounds(b, tolerance);
  if (bounds_a.is_empty() || bounds_b.is_empty()) {
    result.issues.push_back(detail::ssi_driver_issue(SsiDiagnostic::degenerate_patch,
                                                     "ssi: patch bounds are not finite"));
    return result;
  }
  const double diagonal_a = norm2(bounds_a.diagonal());
  const double diagonal_b = norm2(bounds_b.diagonal());
  if (!std::isfinite(diagonal_a) || !std::isfinite(diagonal_b) || diagonal_a <= tolerance ||
      diagonal_b <= tolerance) {
    result.issues.push_back(detail::ssi_driver_issue(
        SsiDiagnostic::degenerate_patch, "ssi: patch is degenerate within tolerance"));
    return result;
  }

  if (surfaces_coincident(a, b, tolerance)) {
    result.issues.push_back(detail::ssi_driver_issue(
        SsiDiagnostic::coincident, "ssi: surfaces are coincident over the shared domain"));
    return result;
  }
  if (patches_tangent(a, b, tolerance)) {
    result.issues.push_back(detail::ssi_driver_issue(
        SsiDiagnostic::tangent_contact, "ssi: patches touch tangentially"));
    return result;
  }

  // A band below the arithmetic resolution of the coordinates cannot be
  // certified: the evaluator resolves points to a few ulp of the model scale, so
  // a tighter request accepts only bit-identical evaluations and reports a
  // residual that is zero by construction. The floor is two ulp of the largest
  // coordinate magnitude; requests at or above it trace, requests below it are
  // diagnosed instead of being returned as an empty "converged" result.
  double max_abs = 0.0;
  for (const Box3d* bounds : {&bounds_a, &bounds_b}) {
    for (int axis = 0; axis < 3; ++axis) {
      max_abs = std::max(max_abs, std::fabs(bounds->lo[axis]));
      max_abs = std::max(max_abs, std::fabs(bounds->hi[axis]));
    }
  }
  const double resolution =
      2.0 * std::numeric_limits<double>::epsilon() * std::max(1.0, max_abs);
  if (tolerance < resolution) {
    result.issues.push_back(detail::ssi_driver_issue(
        SsiDiagnostic::no_convergence,
        "ssi: tolerance is below the coordinate resolution of the patches"));
    return result;
  }

  Box3d bounds = bounds_a;
  bounds.expand(bounds_b);
  const double diagonal = norm2(bounds.diagonal());
  if (!std::isfinite(diagonal) || !(diagonal > 0.0)) {
    result.issues.push_back(detail::ssi_driver_issue(SsiDiagnostic::degenerate_patch,
                                                     "ssi: patch bounds are degenerate"));
    return result;
  }
  const double chord_tolerance = std::max(1e-3 * diagonal, 64.0 * tolerance);
  const Result<SsiPatchMesh> mesh_a = tessellate_ssi_patch(a, chord_tolerance, options.max_depth);
  if (!mesh_a.has_value()) {
    result.issues.push_back(detail::ssi_driver_issue(SsiDiagnostic::degenerate_patch,
                                                     mesh_a.error().message));
    return result;
  }
  const Result<SsiPatchMesh> mesh_b = tessellate_ssi_patch(b, chord_tolerance, options.max_depth);
  if (!mesh_b.has_value()) {
    result.issues.push_back(detail::ssi_driver_issue(SsiDiagnostic::degenerate_patch,
                                                     mesh_b.error().message));
    return result;
  }
  const SsiPatchMesh& patch_mesh_a = mesh_a.value();
  const SsiPatchMesh& patch_mesh_b = mesh_b.value();
  if (patch_mesh_a.triangles.empty() || patch_mesh_b.triangles.empty()) {
    result.issues.push_back(detail::ssi_driver_issue(
        SsiDiagnostic::degenerate_patch, "ssi: patch tessellation produced no triangles"));
    return result;
  }

  // Broad phase: one AABB per triangle, built in mesh order so the seed order is
  // a pure function of the tessellation.
  std::vector<Box3d> boxes_a;
  boxes_a.reserve(patch_mesh_a.triangles.size());
  for (const std::array<std::uint32_t, 3>& triangle : patch_mesh_a.triangles) {
    boxes_a.push_back(detail::ssi_driver_triangle_bounds(patch_mesh_a, triangle));
  }
  std::vector<Box3d> boxes_b;
  boxes_b.reserve(patch_mesh_b.triangles.size());
  for (const std::array<std::uint32_t, 3>& triangle : patch_mesh_b.triangles) {
    boxes_b.push_back(detail::ssi_driver_triangle_bounds(patch_mesh_b, triangle));
  }

  std::vector<SsiPoint> seeds;
  std::vector<std::array<std::size_t, 2>> segments;
  for (std::size_t index_a = 0; index_a < patch_mesh_a.triangles.size(); ++index_a) {
    const std::array<std::uint32_t, 3>& triangle_a = patch_mesh_a.triangles[index_a];
    for (std::size_t index_b = 0; index_b < patch_mesh_b.triangles.size(); ++index_b) {
      if (!detail::ssi_driver_bounds_overlap(boxes_a[index_a], boxes_b[index_b])) continue;
      ++result.candidate_pairs;
      const std::array<std::uint32_t, 3>& triangle_b = patch_mesh_b.triangles[index_b];
      const TriTriResult hit = tri_tri_intersect(
          patch_mesh_a.points[triangle_a[0]], patch_mesh_a.points[triangle_a[1]],
          patch_mesh_a.points[triangle_a[2]], patch_mesh_b.points[triangle_b[0]],
          patch_mesh_b.points[triangle_b[1]], patch_mesh_b.points[triangle_b[2]], tolerance);
      if (hit.kind != TriTriKind::segment) continue;
      std::array<std::size_t, 2> endpoints{detail::kSsiDriverNoSeed, detail::kSsiDriverNoSeed};
      const Point3d positions[2] = {hit.a, hit.b};
      for (int endpoint = 0; endpoint < 2; ++endpoint) {
        std::array<double, 2> ua{};
        std::array<double, 2> va{};
        if (!detail::ssi_driver_triangle_uv(patch_mesh_a, triangle_a, positions[endpoint], ua)) {
          continue;
        }
        if (!detail::ssi_driver_triangle_uv(patch_mesh_b, triangle_b, positions[endpoint], va)) {
          continue;
        }
        const Result<SsiPoint> refined =
            refine_ssi_point(a, b, {ua[0], ua[1], va[0], va[1]}, tolerance);
        if (!refined.has_value()) continue;
        if (!(refined.value().residual <= tolerance)) continue;
        endpoints[static_cast<std::size_t>(endpoint)] = seeds.size();
        seeds.push_back(refined.value());
      }
      segments.push_back(endpoints);
    }
  }

  const double seed_tolerance =
      std::max(64.0 * tolerance, 1e-9 * std::max(diagonal_a, diagonal_b));
  const std::vector<SsiPoint> representatives =
      detail::ssi_driver_curve_seeds(a, b, seeds, segments, seed_tolerance);
  Result<std::vector<SsiCurve>> traced = trace_ssi_curves(a, b, representatives, tolerance);
  if (!traced.has_value()) {
    // An all-degenerate seed set is a diagnostic at the driver layer, not a
    // failed request: the caller still receives a result that explains why no
    // curve was produced. Every other trace error is a contract failure and is
    // propagated unchanged.
    if (std::strcmp(traced.error().message,
                    "trace_ssi_curves: seed has no traceable extent") != 0) {
      return traced.error();
    }
  } else {
    result.curves = std::move(traced.value());
  }
  if (result.curves.empty() && !segments.empty()) {
    result.issues.push_back(detail::ssi_driver_issue(
        SsiDiagnostic::no_convergence, "ssi: no seed converged to an intersection curve"));
  }

  std::uint64_t accepted_segments = 0;
  double max_residual = 0.0;
  for (const SsiCurve& curve : result.curves) {
    accepted_segments += static_cast<std::uint64_t>(curve.segments.size());
    for (const SsiSegment& segment : curve.segments) {
      max_residual = std::max(max_residual, segment.residual);
    }
  }
  if (accepted_segments > options.max_segments) {
    return Error{Errc::numeric_failure, "ssi: segment budget exhausted"};
  }
  result.accepted_segments = accepted_segments;
  result.max_residual = max_residual;
  const std::uint64_t extracted_segments = static_cast<std::uint64_t>(segments.size());
  result.rejected_segments =
      extracted_segments > accepted_segments ? extracted_segments - accepted_segments : 0;

  // Curve order is a pure function of geometry: first point first, then the
  // parameterization and the segment count to break ties. Diagnostics are
  // ordered by code and parameters so equal inputs report issues identically;
  // ssi() emits at most one diagnostic today, so this sort is defensive.
  std::stable_sort(result.curves.begin(), result.curves.end(),
                   [](const SsiCurve& left, const SsiCurve& right) {
                     const SsiSegment& l = left.segments.front();
                     const SsiSegment& r = right.segments.front();
                     for (int axis = 0; axis < 3; ++axis) {
                       if (l.a[axis] != r.a[axis]) return l.a[axis] < r.a[axis];
                     }
                     if (l.ua != r.ua) return l.ua < r.ua;
                     if (l.va != r.va) return l.va < r.va;
                     return left.segments.size() < right.segments.size();
                   });
  std::stable_sort(result.issues.begin(), result.issues.end(),
                   [](const SsiIssue& left, const SsiIssue& right) {
                     if (left.code != right.code) {
                       return static_cast<int>(left.code) < static_cast<int>(right.code);
                     }
                     return left.params < right.params;
                   });
  return result;
}

}  // namespace v1
}  // namespace geom
