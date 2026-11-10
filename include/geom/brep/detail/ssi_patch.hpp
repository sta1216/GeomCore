#pragma once

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <geom/brep/detail/ssi_surface.hpp>
#include <geom/brep/topology.hpp>
#include <geom/core/error.hpp>
#include <geom/core/result.hpp>
#include <geom/primitives/box.hpp>
#include <limits>
#include <map>
#include <utility>
#include <variant>
#include <vector>

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
//
// max_chord_error is the largest sampled leaf flatness that was accepted; it is
// an estimate, not a bound. tolerance_met is false when at least one leaf still
// exceeded the requested chord tolerance at max_depth.
struct SsiPatchMesh {
  std::vector<Point3d> points;
  std::vector<Vec3d> normals;
  std::vector<std::array<double, 2>> uv;
  std::vector<std::array<std::uint32_t, 3>> triangles;
  double max_chord_error = 0.0;
  bool tolerance_met = true;
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

// patch_bounds uses a fixed sampling grid rather than adaptive subdivision, so
// its cost does not depend on chord_tolerance. This many cells are sampled per
// parameter axis, giving a (N + 1) x (N + 1) grid of corner samples.
inline constexpr int kSsiPatchBoundsCells = 16;

// Bilinear interpolation of the cell's four corner positions at (s, t), with
// s and t in [0, 1]. This is the surface two corner-only triangles would
// interpolate, so sample deviations from it measure the error of a cell.
inline Point3d ssi_patch_bilinear(const Point3d& p00, const Point3d& p10, const Point3d& p11,
                                  const Point3d& p01, double s, double t) noexcept {
  const double w00 = (1.0 - s) * (1.0 - t);
  const double w10 = s * (1.0 - t);
  const double w11 = s * t;
  const double w01 = (1.0 - s) * t;
  Point3d point;
  for (int axis = 0; axis < 3; ++axis) {
    point[axis] = w00 * p00[axis] + w10 * p10[axis] + w11 * p11[axis] + w01 * p01[axis];
  }
  return point;
}

// Distance from one parameter sample to the cell's bilinear surface. Measuring
// against the bilinear, rather than the corner plane, makes this sensitive to
// saddles whose error is zero on the edges and at the center.
inline double ssi_patch_sample_deviation(const Point3d& p00, const Point3d& p10, const Point3d& p11,
                                         const Point3d& p01, const Point3d& sample, double s,
                                         double t) noexcept {
  return norm2(sample - ssi_patch_bilinear(p00, p10, p11, p01, s, t));
}

// The nine samples a cell evaluates: center, four edge midpoints, and four
// quarter points at (0.25, 0.25), (0.75, 0.25), (0.75, 0.75), (0.25, 0.75).
struct SsiPatchCellSamples {
  SurfacePoint center;
  SurfacePoint edge_bottom;
  SurfacePoint edge_right;
  SurfacePoint edge_top;
  SurfacePoint edge_left;
  SurfacePoint quarter_00;
  SurfacePoint quarter_10;
  SurfacePoint quarter_11;
  SurfacePoint quarter_01;
};

// Sampled flatness of a cell: the largest distance from any of its nine samples
// to the bilinear surface through its four corners. The quarter points are what
// catch a bicubic saddle, whose error vanishes at all corners, edge midpoints,
// and the center. This is an estimate, not a mathematical bound.
inline double ssi_patch_cell_flatness(const SurfacePoint& p00, const SurfacePoint& p10,
                                      const SurfacePoint& p11, const SurfacePoint& p01,
                                      const SsiPatchCellSamples& samples) noexcept {
  double flatness =
      ssi_patch_sample_deviation(p00.p, p10.p, p11.p, p01.p, samples.center.p, 0.5, 0.5);
  flatness = std::max(flatness, ssi_patch_sample_deviation(p00.p, p10.p, p11.p, p01.p,
                                                           samples.edge_bottom.p, 0.5, 0.0));
  flatness = std::max(flatness, ssi_patch_sample_deviation(p00.p, p10.p, p11.p, p01.p,
                                                           samples.edge_right.p, 1.0, 0.5));
  flatness = std::max(flatness, ssi_patch_sample_deviation(p00.p, p10.p, p11.p, p01.p,
                                                           samples.edge_top.p, 0.5, 1.0));
  flatness = std::max(flatness, ssi_patch_sample_deviation(p00.p, p10.p, p11.p, p01.p,
                                                           samples.edge_left.p, 0.0, 0.5));
  flatness = std::max(flatness, ssi_patch_sample_deviation(p00.p, p10.p, p11.p, p01.p,
                                                           samples.quarter_00.p, 0.25, 0.25));
  flatness = std::max(flatness, ssi_patch_sample_deviation(p00.p, p10.p, p11.p, p01.p,
                                                           samples.quarter_10.p, 0.75, 0.25));
  flatness = std::max(flatness, ssi_patch_sample_deviation(p00.p, p10.p, p11.p, p01.p,
                                                           samples.quarter_11.p, 0.75, 0.75));
  flatness = std::max(flatness, ssi_patch_sample_deviation(p00.p, p10.p, p11.p, p01.p,
                                                           samples.quarter_01.p, 0.25, 0.75));
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
  double max_chord_error = 0.0;
  bool tolerance_met = true;
};

inline Result<std::uint32_t> ssi_patch_emit_vertex(SsiPatchMeshBuilder& builder,
                                                   const SurfacePoint& sample, double u, double v) {
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
// in, so a cell evaluates only its nine interior samples. Leaves emit four
// triangles around the center, wound CCW in (u, v). The root cell is level 1.
//
// Adaptive leaves at different depths are not edge-conforming: a coarse leaf's
// edge is not split where a finer neighbour introduces a midpoint, so hanging
// nodes, and cracks of the order of the coarse leaf's sampled error, are
// possible. At a pole several parameter samples may evaluate to the same
// position; those vertices are kept, with distinct parameters.
inline Result<void> ssi_patch_tessellate_cell(SsiPatchMeshBuilder& builder, const SurfacePoint& p00,
                                              const SurfacePoint& p10, const SurfacePoint& p11,
                                              const SurfacePoint& p01, double u0, double u1,
                                              double v0, double v1, int depth) {
  const double um = u0 + 0.5 * (u1 - u0);
  const double vm = v0 + 0.5 * (v1 - v0);
  if (!(um > u0 && um < u1 && vm > v0 && vm < v1)) {
    return Error{Errc::degenerate, "ssi: patch cell cannot be subdivided further"};
  }
  const double uq0 = u0 + 0.25 * (u1 - u0);
  const double uq1 = u0 + 0.75 * (u1 - u0);
  const double vq0 = v0 + 0.25 * (v1 - v0);
  const double vq1 = v0 + 0.75 * (v1 - v0);
  const std::array<std::pair<double, double>, 9> parameters{{{um, vm},
                                                             {um, v0},
                                                             {u1, vm},
                                                             {um, v1},
                                                             {u0, vm},
                                                             {uq0, vq0},
                                                             {uq1, vq0},
                                                             {uq1, vq1},
                                                             {uq0, vq1}}};
  std::array<SurfacePoint, 9> samples{};
  for (std::size_t i = 0; i < samples.size(); ++i) {
    const Result<SurfacePoint> sample =
        surface_eval(builder.patch->surface, parameters[i].first, parameters[i].second);
    if (!sample.has_value()) return sample.error();
    if (!surface_point_finite(sample.value())) {
      return Error{Errc::non_finite, "ssi: patch sample is not finite"};
    }
    samples[i] = sample.value();
  }
  SsiPatchCellSamples cell;
  cell.center = samples[0];
  cell.edge_bottom = samples[1];
  cell.edge_right = samples[2];
  cell.edge_top = samples[3];
  cell.edge_left = samples[4];
  cell.quarter_00 = samples[5];
  cell.quarter_10 = samples[6];
  cell.quarter_11 = samples[7];
  cell.quarter_01 = samples[8];

  const double flatness = ssi_patch_cell_flatness(p00, p10, p11, p01, cell);
  if (flatness > builder.chord_tolerance && depth < builder.max_depth) {
    const Result<void> lower_left = ssi_patch_tessellate_cell(
        builder, p00, cell.edge_bottom, cell.center, cell.edge_left, u0, um, v0, vm, depth + 1);
    if (!lower_left.has_value()) return lower_left.error();
    const Result<void> lower_right = ssi_patch_tessellate_cell(
        builder, cell.edge_bottom, p10, cell.edge_right, cell.center, um, u1, v0, vm, depth + 1);
    if (!lower_right.has_value()) return lower_right.error();
    const Result<void> upper_left = ssi_patch_tessellate_cell(
        builder, cell.edge_left, cell.center, cell.edge_top, p01, u0, um, vm, v1, depth + 1);
    if (!upper_left.has_value()) return upper_left.error();
    const Result<void> upper_right = ssi_patch_tessellate_cell(
        builder, cell.center, cell.edge_right, p11, cell.edge_top, um, u1, vm, v1, depth + 1);
    if (!upper_right.has_value()) return upper_right.error();
    return Result<void>{};
  }

  builder.max_chord_error = std::max(builder.max_chord_error, flatness);
  if (flatness > builder.chord_tolerance) builder.tolerance_met = false;

  const Result<std::uint32_t> i00 = ssi_patch_emit_vertex(builder, p00, u0, v0);
  if (!i00.has_value()) return i00.error();
  const Result<std::uint32_t> i10 = ssi_patch_emit_vertex(builder, p10, u1, v0);
  if (!i10.has_value()) return i10.error();
  const Result<std::uint32_t> i11 = ssi_patch_emit_vertex(builder, p11, u1, v1);
  if (!i11.has_value()) return i11.error();
  const Result<std::uint32_t> i01 = ssi_patch_emit_vertex(builder, p01, u0, v1);
  if (!i01.has_value()) return i01.error();
  const Result<std::uint32_t> ic = ssi_patch_emit_vertex(builder, cell.center, um, vm);
  if (!ic.has_value()) return ic.error();

  builder.mesh.triangles.push_back({i00.value(), i10.value(), ic.value()});
  builder.mesh.triangles.push_back({i10.value(), i11.value(), ic.value()});
  builder.mesh.triangles.push_back({i11.value(), i01.value(), ic.value()});
  builder.mesh.triangles.push_back({i01.value(), i00.value(), ic.value()});
  return Result<void>{};
}

}  // namespace detail

// Adaptively tessellates a patch until every leaf cell's nine-sample flatness is
// at most chord_tolerance, or until max_depth cell levels have been reached.
// The root cell is level 1, so max_depth == 1 emits a single cell and must be
// at least 1. Flatness is sampled at nine points per cell and measured against
// the cell's bilinear surface; it is a sampled estimate, not a mathematical
// certificate. When the depth cap is reached with a leaf still above tolerance,
// that leaf is emitted, tolerance_met is set false, and max_chord_error records
// the worst accepted sampled flatness. Leaves of different depths are not
// edge-conforming, so cracks of the order of the coarse leaf error are possible.
// The tolerance is the only flatness authority: an invalid value is rejected
// instead of defaulting. The patch is revalidated with make_ssi_patch, so a
// hand-built patch cannot bypass the surface and domain checks. Output depends
// only on the patch and tolerance: child order, vertex emission order, and
// triangle order are fixed.
inline Result<SsiPatchMesh> tessellate_ssi_patch(const SsiPatch& patch, double chord_tolerance,
                                                 int max_depth) {
  if (!std::isfinite(chord_tolerance) || !(chord_tolerance > 0.0)) {
    return Error{Errc::tolerance_violation, "ssi: chord tolerance must be finite and positive"};
  }
  if (max_depth < 1) {
    return Error{Errc::index_out_of_range, "ssi: patch max depth must be positive"};
  }
  const Result<SsiPatch> valid = make_ssi_patch(patch.surface, patch.u, patch.v);
  if (!valid.has_value()) return valid.error();

  detail::SsiPatchMeshBuilder builder;
  builder.patch = &valid.value();
  builder.chord_tolerance = chord_tolerance;
  builder.max_depth = max_depth;

  const Result<SurfacePoint> p00 =
      surface_eval(valid.value().surface, valid.value().u[0], valid.value().v[0]);
  if (!p00.has_value()) return p00.error();
  const Result<SurfacePoint> p10 =
      surface_eval(valid.value().surface, valid.value().u[1], valid.value().v[0]);
  if (!p10.has_value()) return p10.error();
  const Result<SurfacePoint> p11 =
      surface_eval(valid.value().surface, valid.value().u[1], valid.value().v[1]);
  if (!p11.has_value()) return p11.error();
  const Result<SurfacePoint> p01 =
      surface_eval(valid.value().surface, valid.value().u[0], valid.value().v[1]);
  if (!p01.has_value()) return p01.error();
  if (!surface_point_finite(p00.value()) || !surface_point_finite(p10.value()) ||
      !surface_point_finite(p11.value()) || !surface_point_finite(p01.value())) {
    return Error{Errc::non_finite, "ssi: patch corner is not finite"};
  }

  const Result<void> status = detail::ssi_patch_tessellate_cell(
      builder, p00.value(), p10.value(), p11.value(), p01.value(), valid.value().u[0],
      valid.value().u[1], valid.value().v[0], valid.value().v[1], 1);
  if (!status.has_value()) return status.error();
  builder.mesh.max_chord_error = builder.max_chord_error;
  builder.mesh.tolerance_met = builder.tolerance_met;
  return std::move(builder.mesh);
}

// Conservative axis-aligned bounds for a patch, built from a fixed 17 x 17 grid
// of corner samples plus a 3 x 3 subgrid inside every cell. Each subgrid sample
// is compared with its cell's bilinear surface; the largest such sag defines a
// padding of twice that sag plus 1e-12 times the box diagonal, which is added on
// every axis. This is a sampled conservative box, not an exact bound, and its
// size does not depend on chord_tolerance: the tolerance is validated for
// consistency with tessellate_ssi_patch but is otherwise unused. Invalid input,
// including an invalid tolerance, produces an empty box.
inline Box3d patch_bounds(const SsiPatch& patch, double chord_tolerance) {
  if (!std::isfinite(chord_tolerance) || !(chord_tolerance > 0.0)) return Box3d::empty();
  const Result<SsiPatch> valid = make_ssi_patch(patch.surface, patch.u, patch.v);
  if (!valid.has_value()) return Box3d::empty();
  const SsiPatch& target = valid.value();
  const int cells = detail::kSsiPatchBoundsCells;
  const int grid_width = cells + 1;
  std::vector<Point3d> grid(static_cast<std::size_t>(grid_width) * grid_width);
  Box3d bounds = Box3d::empty();
  for (int iu = 0; iu < grid_width; ++iu) {
    for (int iv = 0; iv < grid_width; ++iv) {
      const double u = target.u[0] + (target.u[1] - target.u[0]) *
                                         (static_cast<double>(iu) / static_cast<double>(cells));
      const double v = target.v[0] + (target.v[1] - target.v[0]) *
                                         (static_cast<double>(iv) / static_cast<double>(cells));
      const Result<SurfacePoint> sample = surface_eval(target.surface, u, v);
      if (!sample.has_value()) return Box3d::empty();
      if (!surface_point_finite(sample.value())) return Box3d::empty();
      grid[static_cast<std::size_t>(iu) * grid_width + iv] = sample.value().p;
      bounds.expand(sample.value().p);
    }
  }
  if (bounds.is_empty()) return bounds;

  double max_sag = 0.0;
  for (int iu = 0; iu < cells; ++iu) {
    for (int iv = 0; iv < cells; ++iv) {
      const double u0 = target.u[0] + (target.u[1] - target.u[0]) *
                                          (static_cast<double>(iu) / static_cast<double>(cells));
      const double u1 = target.u[0] + (target.u[1] - target.u[0]) * (static_cast<double>(iu + 1) /
                                                                     static_cast<double>(cells));
      const double v0 = target.v[0] + (target.v[1] - target.v[0]) *
                                          (static_cast<double>(iv) / static_cast<double>(cells));
      const double v1 = target.v[0] + (target.v[1] - target.v[0]) * (static_cast<double>(iv + 1) /
                                                                     static_cast<double>(cells));
      const Point3d& p00 = grid[static_cast<std::size_t>(iu) * grid_width + iv];
      const Point3d& p10 = grid[static_cast<std::size_t>(iu + 1) * grid_width + iv];
      const Point3d& p11 = grid[static_cast<std::size_t>(iu + 1) * grid_width + iv + 1];
      const Point3d& p01 = grid[static_cast<std::size_t>(iu) * grid_width + iv + 1];
      for (int su = 0; su < 3; ++su) {
        for (int sv = 0; sv < 3; ++sv) {
          const double s = 0.5 * static_cast<double>(su);
          const double t = 0.5 * static_cast<double>(sv);
          const double u = u0 + s * (u1 - u0);
          const double v = v0 + t * (v1 - v0);
          const Result<SurfacePoint> sample = surface_eval(target.surface, u, v);
          if (!sample.has_value()) return Box3d::empty();
          if (!surface_point_finite(sample.value())) return Box3d::empty();
          max_sag = std::max(max_sag, detail::ssi_patch_sample_deviation(p00, p10, p11, p01,
                                                                         sample.value().p, s, t));
          bounds.expand(sample.value().p);
        }
      }
    }
  }

  if (!std::isfinite(max_sag)) return Box3d::empty();
  const double diagonal = norm2(bounds.hi - bounds.lo);
  if (!std::isfinite(diagonal)) return Box3d::empty();
  const double padding = 2.0 * max_sag + 1e-12 * diagonal;
  if (!std::isfinite(padding)) return Box3d::empty();
  for (int axis = 0; axis < 3; ++axis) {
    bounds.lo[axis] -= padding;
    bounds.hi[axis] += padding;
  }
  return bounds;
}

}  // namespace v1
}  // namespace geom
