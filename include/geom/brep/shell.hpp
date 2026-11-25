#pragma once

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <span>
#include <utility>
#include <vector>

#include <geom/brep/boolean_solid.hpp>
#include <geom/brep/detail/blend_common.hpp>
#include <geom/brep/primitives.hpp>
#include <geom/brep/validate.hpp>
#include <geom/core/error.hpp>
#include <geom/core/result.hpp>
#include <geom/surfaces/plane.hpp>

namespace geom {
inline namespace v1 {

struct ShellOptions {
  double thickness = 0.0;  // uniform wall thickness, > 0
  double tolerance = 1e-9;
};

struct ShellResult {
  Brep brep;
  std::uint32_t cavity_shells = 0;  // enclosed pockets; always >= 1 on success
};


// Hollows a convex planar solid to a uniform wall of `thickness`, leaving an
// inner cavity. M4c-1 accepts closed, axis-aligned solids only: an open face
// list, an oblique face, a non-convex body, or a curved face is refused with a
// stable `Errc::unsupported` diagnostic rather than returning an order-
// dependent solid. M4c-2 lifts the open-face and oblique restrictions together
// with the splitter/boolean face-order fix.
//
// The cavity is the axis-aligned box inset by `thickness` on every side, then
// the final regularized difference A - cavity leaves the uniform wall. The box
// is built directly instead of being cut out of a seed: every accepted face is
// axis-aligned, so a convex solid with only axis-aligned faces is a box, and
// constructing it directly avoids the order-sensitive splitter path.
//
// A thickness that consumes the body in any axis leaves no interior and is
// reported as degenerate.
//
// Diagnostics are static texts because Error::message is a const char*.
inline Result<ShellResult> shell_solid(const Brep& brep, std::span<const FaceId> open_faces = {},
                                       const ShellOptions& options = {}) {
  if (!std::isfinite(options.thickness)) {
    return Error{Errc::non_finite, "shell_solid: thickness must be finite"};
  }
  if (!(options.thickness > 0.0)) {
    return Error{Errc::degenerate, "shell_solid: thickness must be positive"};
  }
  if (!std::isfinite(options.tolerance)) {
    return Error{Errc::non_finite, "shell_solid: tolerance must be finite"};
  }
  if (!(options.tolerance > 0.0)) {
    return Error{Errc::tolerance_violation, "shell_solid: tolerance must be positive"};
  }

  // Shape, planarity, closedness and convexity are properties of the whole
  // solid and are checked before any individual face is looked at.
  const Result<detail::BlendPlanarSolid> view =
      detail::blend_planar_convex_solid(brep, options.tolerance, "shell_solid");
  if (!view.has_value()) return view.error();

  // Requested open faces are validated against the single shell of the input
  // and deduplicated by marking their shell position; a repeat is a no-op.
  std::vector<char> is_open(view.value().faces.size(), 0);
  for (const FaceId face : open_faces) {
    if (!brep.face_alive(face)) {
      return Error{Errc::index_out_of_range,
                   "shell_solid: open face does not belong to the solid"};
    }
    const auto found = std::find(view.value().faces.begin(), view.value().faces.end(), face);
    if (found == view.value().faces.end()) {
      return Error{Errc::index_out_of_range,
                   "shell_solid: open face does not belong to the solid"};
    }
    is_open[static_cast<std::size_t>(found - view.value().faces.begin())] = 1;
  }
  const bool any_open =
      std::any_of(is_open.begin(), is_open.end(), [](char flag) { return flag != 0; });
  if (any_open) {
    return Error{Errc::unsupported, "shell_solid: open faces are not supported in M4c-1"};
  }

  const double tolerance = options.tolerance;

  // This slice of the heal/blend milestone supports closed, axis-aligned
  // convex solids (boxes and prisms built from box faces). Sequential
  // half-space clipping and the boolean stitch are order sensitive for open
  // pockets and oblique faces, so refuse those inputs explicitly rather than
  // returning an order-dependent result. M4c-2 lifts the restriction together
  // with the splitter/boolean ordering fix.
  std::size_t oblique_faces = 0;
  for (std::size_t i = 0; i < view.value().faces.size(); ++i) {
    int non_axis = 0;
    for (int axis = 0; axis < 3; ++axis) {
      if (std::abs(view.value().outward[i][axis]) > 1e-9) ++non_axis;
    }
    if (non_axis > 1) ++oblique_faces;
  }
  if (oblique_faces > 0u) {
    return Error{Errc::unsupported, "shell_solid: only axis-aligned faces are supported"};
  }

  // Every accepted face is axis-aligned. A closed convex solid whose faces are
  // all axis-aligned is a box: intersecting the half-spaces leaves one interval
  // per axis, independent of the other axes. The cavity is therefore the input
  // bounding box inset by the wall thickness on every side.
  Point3d lo = brep.vertex(brep.vertex_ids().front()).point;
  Point3d hi = lo;
  for (const VertexId id : brep.vertex_ids()) {
    const Point3d& point = brep.vertex(id).point;
    for (int axis = 0; axis < 3; ++axis) {
      lo[axis] = std::min(lo[axis], point[axis]);
      hi[axis] = std::max(hi[axis], point[axis]);
    }
  }
  const Point3d inner_lo{lo[0] + options.thickness, lo[1] + options.thickness,
                         lo[2] + options.thickness};
  const Point3d inner_hi{hi[0] - options.thickness, hi[1] - options.thickness,
                         hi[2] - options.thickness};
  for (int axis = 0; axis < 3; ++axis) {
    if (!(inner_hi[axis] > inner_lo[axis])) {
      return Error{Errc::degenerate, "shell_solid: wall thickness leaves no interior"};
    }
  }
  const Result<Brep> inner_box = make_box(inner_lo, inner_hi, tolerance);
  if (!inner_box.has_value()) {
    return Error{Errc::degenerate, "shell_solid: wall thickness leaves no measurable interior"};
  }

  Brep inner = inner_box.value();

  SolidBooleanOptions boolean_options;
  boolean_options.tolerance = tolerance;
  const Result<SolidBooleanResult> difference = boolean_solid(
      brep, inner, SolidBooleanOp::difference, boolean_options);
  if (!difference.has_value()) {
    // A wall that is only a few tolerances wide cannot be stitched back onto
    // the outer solid; report the geometric limit rather than a stitcher error.
    if (difference.error().code == Errc::non_finite) {
      return Error{Errc::non_finite, "shell_solid: cavity difference failed a numeric check"};
    }
    return Error{Errc::degenerate,
                 "shell_solid: wall thickness leaves no measurable interior at this tolerance"};
  }

  Brep shelled = std::move(difference.value().brep);
  if (shelled.solid_count() != 1u) {
    return Error{Errc::invalid_topology, "shell_solid: difference did not leave one solid"};
  }
  if (!validate_topology(shelled).ok()) {
    return Error{Errc::invalid_topology, "shell_solid: shell produced an invalid solid"};
  }

  ShellResult result;
  if (any_open) {
    result.cavity_shells = 0u;
  } else {
    const std::size_t shells = shelled.solid_shells(shelled.solid_ids().front()).size();
    if (shells < 2u) {
      return Error{Errc::invalid_topology, "shell_solid: closed shell produced no cavity"};
    }
    result.cavity_shells = static_cast<std::uint32_t>(shells - 1u);
  }
  result.brep = std::move(shelled);
  return result;
}

}  // namespace v1
}  // namespace geom
