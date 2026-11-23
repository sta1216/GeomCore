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
  // Saturated accounting difference max(0, extracted - accepted). Trace
  // refinement can emit many output segments per extracted seed segment, so
  // accepted may exceed the extracted count; this is therefore not a count of
  // untouched geometry, only of the unconsumed extraction surplus.
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
  // Exact analytic sphere-plane chart. Callers that need the traced
  // parameterization (the boolean driver, whose seam handling is built around
  // it) can opt out; public callers get the exact path by default.
  bool analytic_sphere_plane = true;
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

// Exact sphere/plane fast path. A plane perpendicular to the sphere's axis
// cuts the sphere in one latitude circle, so the intersection can be emitted
// without tessellation or Newton refinement. The path is taken only when the
// whole circle lies strictly inside both patch domains and the result is
// certifiable at the requested tolerance; every other case returns false so
// the generic driver produces the same classification and clipping it always
// did.
inline bool ssi_try_analytic_sphere_plane(const SsiPatch& a, const SsiPatch& b,
                                          const SsiOptions& options,
                                          double chord_tolerance,
                                          SsiResult& result) {
  const Sphere3* sphere = std::get_if<Sphere3>(&a.surface);
  const Plane3* plane = std::get_if<Plane3>(&b.surface);
  bool swapped = false;
  if (sphere == nullptr || plane == nullptr) {
    sphere = std::get_if<Sphere3>(&b.surface);
    plane = std::get_if<Plane3>(&a.surface);
    swapped = true;
  }
  if (sphere == nullptr || plane == nullptr) return false;
  if (!sphere->is_valid()) return false;
  const SsiPatch& sphere_patch = swapped ? b : a;
  const SsiPatch& plane_patch = swapped ? a : b;
  if (sphere_patch.u[1] - sphere_patch.u[0] < kSsiTwoPi - 1e-12) return false;

  const Result<PlaneFrame> frame_result = make_plane_frame(*plane);
  if (!frame_result.has_value()) return false;
  const PlaneFrame& frame = frame_result.value();
  const double axial_tilt = norm2(cross(frame.normal, sphere->axis));
  if (!(axial_tilt <= 1e-12)) return false;

  const double d = dot(frame.normal, sphere->center - frame.origin);
  const double scale = std::max({1.0, std::fabs(sphere->center[0]),
                                 std::fabs(sphere->center[1]),
                                 std::fabs(sphere->center[2])}) +
                       sphere->radius + std::fabs(plane->offset);
  if (options.tolerance <
      2.0 * std::numeric_limits<double>::epsilon() * scale) {
    return false;
  }
  const double radius = sphere->radius;
  if (!(std::fabs(d) < radius)) return false;
  const double gap = radius - std::fabs(d);
  const double tangent_band =
      std::max(1e-12 * (1.0 + radius + std::fabs(d)), 0.05 * options.tolerance);
  if (gap <= tangent_band) {
    result.issues.push_back(ssi_driver_issue(SsiDiagnostic::tangent_contact,
                                             "ssi: patches touch tangentially"));
    return true;
  }

  const double rho = std::sqrt(std::max(0.0, radius * radius - d * d));
  if (!(rho > chord_tolerance)) return false;
  const double ratio = chord_tolerance / rho;
  if (!(ratio < 1.0)) return false;
  const double step = std::acos(1.0 - ratio);
  if (!(step > 0.0)) return false;
  const double count = std::ceil(0.5 * kSsiTwoPi / step);
  if (!std::isfinite(count) || count < 8.0) return false;
  const auto segments = static_cast<std::uint64_t>(count);
  if (segments > 4096u || segments > options.max_segments) return false;

  const double v = std::asin(std::clamp(-d / radius, -1.0, 1.0));
  const double v_margin = 1e-12 * (1.0 + std::fabs(sphere_patch.v[0]) + std::fabs(sphere_patch.v[1]));
  if (!(v > sphere_patch.v[0] + v_margin) || !(v < sphere_patch.v[1] - v_margin)) {
    return false;
  }

  const Vec3d radial_x = frame.u;
  const Vec3d radial_y = frame.v;
  const Point3d center = sphere->center - frame.normal * d;
  const double u_margin = 1e-12 * (1.0 + std::fabs(plane_patch.u[0]) + std::fabs(plane_patch.u[1]));
  const double w_margin = 1e-12 * (1.0 + std::fabs(plane_patch.v[0]) + std::fabs(plane_patch.v[1]));
  const std::uint64_t sample_count = std::max<std::uint64_t>(segments, 64u);
  for (std::uint64_t i = 0; i < sample_count; ++i) {
    const double angle = kSsiTwoPi * static_cast<double>(i) / static_cast<double>(sample_count);
    const Point3d p = center + radial_x * (rho * std::cos(angle)) + radial_y * (rho * std::sin(angle));
    const Point2d uv = plane_frame_to_2d(frame, p);
    if (!(uv[0] > plane_patch.u[0] + u_margin) || !(uv[0] < plane_patch.u[1] - u_margin) ||
        !(uv[1] > plane_patch.v[0] + w_margin) || !(uv[1] < plane_patch.v[1] - w_margin)) {
      return false;
    }
  }

  const Vec3d y_axis = cross(sphere->axis, sphere->x_axis);
  // Start the loop on the u = 0 seam so the parameter sequence is monotone and
  // the closing segment is the only one that crosses the seam, exactly like a
  // traced closed curve.
  const double phi = std::atan2(dot(radial_x, y_axis), dot(radial_x, sphere->x_axis));
  const double angle_start = -phi;
  std::vector<Point3d> points;
  std::vector<Point2d> plane_uv;
  std::vector<std::array<double, 2>> sphere_uv;
  points.reserve(static_cast<std::size_t>(segments));
  plane_uv.reserve(static_cast<std::size_t>(segments));
  sphere_uv.reserve(static_cast<std::size_t>(segments));
  for (std::uint64_t i = 0; i < segments; ++i) {
    const double angle =
        angle_start + kSsiTwoPi * static_cast<double>(i) / static_cast<double>(segments);
    const Point3d p =
        center + radial_x * (rho * std::cos(angle)) + radial_y * (rho * std::sin(angle));
    const Point2d uv = plane_frame_to_2d(frame, p);
    points.push_back(p);
    plane_uv.push_back(uv);
    sphere_uv.push_back({kSsiTwoPi * static_cast<double>(i) / static_cast<double>(segments), v});
  }

  SsiCurve curve;
  curve.closed = true;
  curve.segments.reserve(static_cast<std::size_t>(segments));
  double worst_residual = 0.0;
  for (std::uint64_t i = 0; i < segments; ++i) {
    const std::uint64_t next = (i + 1) % segments;
    const std::array<double, 2> sphere_a = sphere_uv[static_cast<std::size_t>(i)];
    const std::array<double, 2> plane_a{plane_uv[static_cast<std::size_t>(i)][0],
                                        plane_uv[static_cast<std::size_t>(i)][1]};
    std::array<double, 2> sphere_b = sphere_uv[static_cast<std::size_t>(next)];
    if (next == 0) sphere_b[0] = kSsiTwoPi;  // closing endpoint keeps the winding
    const std::array<double, 2> plane_b{plane_uv[static_cast<std::size_t>(next)][0],
                                        plane_uv[static_cast<std::size_t>(next)][1]};
    const SsiPoint point_a = swapped ? SsiPoint{plane_a, sphere_a, points[static_cast<std::size_t>(i)]}
                                     : SsiPoint{sphere_a, plane_a, points[static_cast<std::size_t>(i)]};
    const SsiPoint point_b = swapped ? SsiPoint{plane_b, sphere_b, points[static_cast<std::size_t>(next)]}
                                     : SsiPoint{sphere_b, plane_b, points[static_cast<std::size_t>(next)]};
    SsiSegment segment;
    segment.a = point_a.p;
    segment.b = point_b.p;
    segment.ua = point_a.ua;
    segment.ub = point_b.ua;
    segment.va = point_a.va;
    segment.vb = point_b.va;
    segment.residual = std::max(ssi_residual(swapped ? plane_patch : sphere_patch,
                                            swapped ? sphere_patch : plane_patch, point_a),
                                ssi_residual(swapped ? plane_patch : sphere_patch,
                                            swapped ? sphere_patch : plane_patch, point_b));
    worst_residual = std::max(worst_residual, segment.residual);
    curve.segments.push_back(segment);
  }
  result.curves.push_back(std::move(curve));
  result.candidate_pairs = 0;
  result.accepted_segments = segments;
  result.max_residual = worst_residual;
  return true;
}

}  // namespace detail

// Traces the intersection of two patches. Degenerate, coincident, and tangent
// input is classified before the tessellator runs; otherwise both patches are
// tessellated adaptively, every AABB-overlapping triangle pair is intersected,
// and the segment endpoints become Newton seeds. Seeds that share a refined
// endpoint belong to the same intersection curve, so one seed per connected
// segment component is traced: one curve yields one traced curve instead of one
// per mesh segment. The counters record the broad phase (candidate pairs), the
// traced output (accepted segments and the worst residual), and the saturated
// difference between extracted triangle-pair segments and accepted output.
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

  // The analytic sphere/plane chart is exact, so it runs before patch_bounds
  // and the tangency predicates: those sample grids that cost milliseconds and
  // are only needed by the generic tracer.
  if (options.analytic_sphere_plane) {
  if (const Sphere3* sphere = std::get_if<Sphere3>(&a.surface);
      sphere != nullptr && std::get_if<Plane3>(&b.surface) != nullptr) {
    SsiResult analytic;
    const double chord = std::max(1e-3 * 2.0 * sphere->radius, 64.0 * tolerance);
    if (detail::ssi_try_analytic_sphere_plane(a, b, options, chord, analytic)) {
      return analytic;
    }
  } else if (const Sphere3* sphere_b = std::get_if<Sphere3>(&b.surface);
             sphere_b != nullptr && std::get_if<Plane3>(&a.surface) != nullptr) {
    SsiResult analytic;
    const double chord = std::max(1e-3 * 2.0 * sphere_b->radius, 64.0 * tolerance);
    if (detail::ssi_try_analytic_sphere_plane(a, b, options, chord, analytic)) {
      return analytic;
    }
  }
  }

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

  // A band below the arithmetic resolution of the coordinates cannot be
  // certified: the evaluator resolves points to a few ulp of the model scale, so
  // a tighter request cannot be certified whether or not a particular evaluation
  // happens to be bit-identical. The floor is two ulp of max(1, largest
  // coordinate magnitude); requests at or above it trace, requests below it are
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
  bool saw_contact = false;
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
      if (hit.kind == TriTriKind::touching || hit.kind == TriTriKind::coplanar) {
        saw_contact = true;
      }
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
  if (result.curves.empty() && (saw_contact || !segments.empty())) {
    result.issues.push_back(detail::ssi_driver_issue(
        SsiDiagnostic::no_convergence,
        "ssi: no intersection curve could be certified at this tolerance"));
  }
  // Empty curves cannot contribute segments; drop them before accounting. This
  // is defensive: tracing never emits a curve without segments today.
  result.curves.erase(
      std::remove_if(result.curves.begin(), result.curves.end(),
                     [](const SsiCurve& curve) { return curve.segments.empty(); }),
      result.curves.end());
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
