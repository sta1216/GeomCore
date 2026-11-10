#pragma once

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <limits>
#include <map>
#include <utility>
#include <variant>
#include <vector>

#include <geom/brep/detail/ssi_surface.hpp>
#include <geom/brep/topology.hpp>
#include <geom/core/error.hpp>
#include <geom/core/result.hpp>
#include <geom/primitives/box.hpp>

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

// Piecewise-linear approximation of a patch. Every vertex carries its patch
// parameters and a unit normal oriented like the surface's cross(du, dv), so
// callers can recover where a point came from without another evaluation.
// Triangle indices are CCW in the patch's (u, v) plane.
struct SsiPatchMesh {
  std::vector<Point3d> points;
  std::vector<Vec3d> normals;
  std::vector<std::array<double, 2>> uv;
  std::vector<std::array<std::uint32_t, 3>> triangles;
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

namespace detail {

// patch_bounds samples the same adaptive grid as tessellation but has no depth
// argument; this cap bounds its work for callers who only need a box.
inline constexpr int kSsiPatchBoundsMaxDepth = 12;

inline double ssi_patch_max_abs_component(const Vec3d& vector) noexcept {
  return std::max({std::fabs(vector[0]), std::fabs(vector[1]), std::fabs(vector[2])});
}

// Distance from the cell center to the plane spanned by its four corners.
// Scaling each diagonal by its largest component before the cross product
// keeps the intermediate finite for cells whose edges are tiny or huge; a
// degenerate (zero-area) cell has no plane to deviate from and reports zero.
inline double ssi_patch_plane_deviation(const SurfacePoint& p00, const SurfacePoint& p10,
                                        const SurfacePoint& p11, const SurfacePoint& p01,
                                        const SurfacePoint& center) noexcept {
  Vec3d first = p11.p - p00.p;
  Vec3d second = p10.p - p01.p;
  const double first_scale = ssi_patch_max_abs_component(first);
  const double second_scale = ssi_patch_max_abs_component(second);
  if (!(first_scale > 0.0) || !(second_scale > 0.0)) return 0.0;
  for (int axis = 0; axis < 3; ++axis) {
    first[axis] /= first_scale;
    second[axis] /= second_scale;
  }
  Vec3d normal = cross(first, second);
  const double normal_scale = ssi_patch_max_abs_component(normal);
  if (!(normal_scale > 0.0)) return 0.0;
  const Vec3d scaled{normal[0] / normal_scale, normal[1] / normal_scale, normal[2] / normal_scale};
  const double length = std::hypot(scaled[0], scaled[1], scaled[2]);
  normal = Vec3d{scaled[0] / length, scaled[1] / length, scaled[2] / length};
  return std::fabs(dot(center.p - p00.p, normal));
}

// Distance from a sampled edge midpoint to the straight chord between that
// edge's corners. This is the sagitta the chord approximation ignores.
inline double ssi_patch_chord_deviation(const Point3d& sample, const Point3d& corner_a,
                                        const Point3d& corner_b) noexcept {
  const double x = sample[0] - 0.5 * corner_a[0] - 0.5 * corner_b[0];
  const double y = sample[1] - 0.5 * corner_a[1] - 0.5 * corner_b[1];
  const double z = sample[2] - 0.5 * corner_a[2] - 0.5 * corner_b[2];
  return std::hypot(x, y, z);
}

// A cell is flat enough when its center stays close to the corner plane and
// each edge midpoint stays close to its corner chord.
inline double ssi_patch_cell_flatness(const SurfacePoint& p00, const SurfacePoint& p10,
                                      const SurfacePoint& p11, const SurfacePoint& p01,
                                      const SurfacePoint& bottom, const SurfacePoint& right,
                                      const SurfacePoint& top, const SurfacePoint& left,
                                      const SurfacePoint& center) noexcept {
  double flatness = ssi_patch_plane_deviation(p00, p10, p11, p01, center);
  flatness = std::max(flatness, ssi_patch_chord_deviation(bottom.p, p00.p, p10.p));
  flatness = std::max(flatness, ssi_patch_chord_deviation(right.p, p10.p, p11.p));
  flatness = std::max(flatness, ssi_patch_chord_deviation(top.p, p01.p, p11.p));
  flatness = std::max(flatness, ssi_patch_chord_deviation(left.p, p00.p, p01.p));
  return flatness;
}

struct SsiPatchMeshBuilder {
  const SsiPatch* patch = nullptr;
  SsiPatchMesh mesh;
  // Dyadic cell corners are recomputed by every incident cell from the same
  // (u, v) values, so exact key equality merges shared vertices deterministically.
  std::map<std::pair<double, double>, std::uint32_t> vertices;
  double chord_tolerance = 0.0;
  int max_depth = 0;
};

inline Result<std::uint32_t> ssi_patch_emit_vertex(SsiPatchMeshBuilder& builder,
                                                   const SurfacePoint& sample, double u,
                                                   double v) {
  const std::pair<double, double> key{u, v};
  const auto existing = builder.vertices.find(key);
  if (existing != builder.vertices.end()) return existing->second;
  if (builder.mesh.points.size() >=
      static_cast<std::size_t>(std::numeric_limits<std::uint32_t>::max())) {
    return Error{Errc::numeric_failure, "ssi: patch mesh exceeds the 32-bit vertex index range"};
  }
  const Result<Vec3d> normal = surface_normal(builder.patch->surface, u, v);
  if (!normal.has_value()) return normal.error();
  const auto index = static_cast<std::uint32_t>(builder.mesh.points.size());
  builder.mesh.points.push_back(sample.p);
  builder.mesh.normals.push_back(normal.value());
  builder.mesh.uv.push_back({u, v});
  builder.vertices.emplace(key, index);
  return index;
}

// Recursive quadtree. The four corner samples and their parameters are passed
// in, so a cell only evaluates its center and four edge midpoints. Leaves emit
// four triangles around the center, wound CCW in (u, v).
inline Result<void> ssi_patch_tessellate_cell(SsiPatchMeshBuilder& builder,
                                              const SurfacePoint& p00, const SurfacePoint& p10,
                                              const SurfacePoint& p11, const SurfacePoint& p01,
                                              double u0, double u1, double v0, double v1,
                                              int depth) {
  const double um = u0 + 0.5 * (u1 - u0);
  const double vm = v0 + 0.5 * (v1 - v0);
  if (!(um > u0 && um < u1 && vm > v0 && vm < v1)) {
    return Error{Errc::degenerate, "ssi: patch cell cannot be subdivided further"};
  }
  const Result<SurfacePoint> center = surface_eval(builder.patch->surface, um, vm);
  if (!center.has_value()) return center.error();
  const Result<SurfacePoint> bottom = surface_eval(builder.patch->surface, um, v0);
  if (!bottom.has_value()) return bottom.error();
  const Result<SurfacePoint> right = surface_eval(builder.patch->surface, u1, vm);
  if (!right.has_value()) return right.error();
  const Result<SurfacePoint> top = surface_eval(builder.patch->surface, um, v1);
  if (!top.has_value()) return top.error();
  const Result<SurfacePoint> left = surface_eval(builder.patch->surface, u0, vm);
  if (!left.has_value()) return left.error();
  if (!surface_point_finite(center.value()) || !surface_point_finite(bottom.value()) ||
      !surface_point_finite(right.value()) || !surface_point_finite(top.value()) ||
      !surface_point_finite(left.value())) {
    return Error{Errc::non_finite, "ssi: patch sample is not finite"};
  }

  const double flatness =
      ssi_patch_cell_flatness(p00, p10, p11, p01, bottom.value(), right.value(), top.value(),
                              left.value(), center.value());
  if (flatness > builder.chord_tolerance && depth < builder.max_depth) {
    const Result<void> lower_left =
        ssi_patch_tessellate_cell(builder, p00, bottom.value(), center.value(), left.value(), u0,
                                  um, v0, vm, depth + 1);
    if (!lower_left.has_value()) return lower_left.error();
    const Result<void> lower_right =
        ssi_patch_tessellate_cell(builder, bottom.value(), p10, right.value(), center.value(), um,
                                  u1, v0, vm, depth + 1);
    if (!lower_right.has_value()) return lower_right.error();
    const Result<void> upper_left =
        ssi_patch_tessellate_cell(builder, left.value(), center.value(), top.value(), p01, u0, um,
                                  vm, v1, depth + 1);
    if (!upper_left.has_value()) return upper_left.error();
    const Result<void> upper_right =
        ssi_patch_tessellate_cell(builder, center.value(), right.value(), p11, top.value(), um, u1,
                                  vm, v1, depth + 1);
    if (!upper_right.has_value()) return upper_right.error();
    return Result<void>{};
  }

  const Result<std::uint32_t> i00 = ssi_patch_emit_vertex(builder, p00, u0, v0);
  if (!i00.has_value()) return i00.error();
  const Result<std::uint32_t> i10 = ssi_patch_emit_vertex(builder, p10, u1, v0);
  if (!i10.has_value()) return i10.error();
  const Result<std::uint32_t> i11 = ssi_patch_emit_vertex(builder, p11, u1, v1);
  if (!i11.has_value()) return i11.error();
  const Result<std::uint32_t> i01 = ssi_patch_emit_vertex(builder, p01, u0, v1);
  if (!i01.has_value()) return i01.error();
  const Result<std::uint32_t> ic = ssi_patch_emit_vertex(builder, center.value(), um, vm);
  if (!ic.has_value()) return ic.error();

  builder.mesh.triangles.push_back({i00.value(), i10.value(), ic.value()});
  builder.mesh.triangles.push_back({i10.value(), i11.value(), ic.value()});
  builder.mesh.triangles.push_back({i11.value(), i01.value(), ic.value()});
  builder.mesh.triangles.push_back({i01.value(), i00.value(), ic.value()});
  return Result<void>{};
}

}  // namespace detail

// Adaptively tessellates a patch until every leaf cell is within
// chord_tolerance of its corner bilinear, or until max_depth subdivisions have
// been made. The tolerance is the only flatness authority: an invalid value is
// rejected instead of defaulting. The patch is revalidated with
// make_ssi_patch, so a hand-built patch cannot bypass the surface and domain
// checks. Output depends only on the patch and tolerance: child order, vertex
// emission order, and triangle order are fixed.
inline Result<SsiPatchMesh> tessellate_ssi_patch(const SsiPatch& patch, double chord_tolerance,
                                                 int max_depth) {
  if (!std::isfinite(chord_tolerance) || !(chord_tolerance > 0.0)) {
    return Error{Errc::tolerance_violation,
                 "ssi: chord tolerance must be finite and positive"};
  }
  if (max_depth < 1) {
    return Error{Errc::index_out_of_range, "ssi: patch max depth must be positive"};
  }
  const Result<SsiPatch> valid = make_ssi_patch(patch.surface, patch.u, patch.v);
  if (!valid.has_value()) return valid.error();

  detail::SsiPatchMeshBuilder builder;
  builder.patch = &patch;
  builder.chord_tolerance = chord_tolerance;
  builder.max_depth = max_depth;

  const Result<SurfacePoint> p00 = surface_eval(patch.surface, patch.u[0], patch.v[0]);
  if (!p00.has_value()) return p00.error();
  const Result<SurfacePoint> p10 = surface_eval(patch.surface, patch.u[1], patch.v[0]);
  if (!p10.has_value()) return p10.error();
  const Result<SurfacePoint> p11 = surface_eval(patch.surface, patch.u[1], patch.v[1]);
  if (!p11.has_value()) return p11.error();
  const Result<SurfacePoint> p01 = surface_eval(patch.surface, patch.u[0], patch.v[1]);
  if (!p01.has_value()) return p01.error();
  if (!surface_point_finite(p00.value()) || !surface_point_finite(p10.value()) ||
      !surface_point_finite(p11.value()) || !surface_point_finite(p01.value())) {
    return Error{Errc::non_finite, "ssi: patch corner is not finite"};
  }

  const Result<void> status = detail::ssi_patch_tessellate_cell(
      builder, p00.value(), p10.value(), p11.value(), p01.value(), patch.u[0], patch.u[1],
      patch.v[0], patch.v[1], 0);
  if (!status.has_value()) return status.error();
  return std::move(builder.mesh);
}

// Axis-aligned bounds of the same adaptive grid used by tessellate_ssi_patch,
// grown outward by chord_tolerance to cover the chord error. Invalid input
// produces an empty box: there is no enclosing box to report.
inline Box3d patch_bounds(const SsiPatch& patch, double chord_tolerance) {
  const Result<SsiPatchMesh> mesh =
      tessellate_ssi_patch(patch, chord_tolerance, detail::kSsiPatchBoundsMaxDepth);
  if (!mesh.has_value()) return Box3d::empty();
  Box3d bounds = Box3d::empty();
  for (const Point3d& point : mesh.value().points) bounds.expand(point);
  if (bounds.is_empty()) return bounds;
  for (int axis = 0; axis < 3; ++axis) {
    bounds.lo[axis] -= chord_tolerance;
    bounds.hi[axis] += chord_tolerance;
  }
  return bounds;
}

}  // namespace v1
}  // namespace geom
