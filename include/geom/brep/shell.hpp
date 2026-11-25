#pragma once

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <span>
#include <utility>
#include <vector>

#include <geom/brep/boolean_solid.hpp>
#include <geom/brep/detail/blend_common.hpp>
#include <geom/brep/primitives.hpp>
#include <geom/brep/split_solid.hpp>
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
  std::uint32_t cavity_shells = 0;  // enclosed pockets; 0 when any face is open
};

// Hollows a convex planar solid to a uniform wall of `thickness`, leaving an
// inner cavity (or, for requested open faces, a pocket that opens through
// them).
//
// The cavity is built as a half-space intersection, not by cutting the input:
// a bounding box around the input, expanded by twice the wall thickness, is
// cut once per face. A closed face is moved inward by `thickness`; an open
// face is moved outward, so its half-space keeps the cavity open through the
// face and lets the inner tool protrude. The final regularized difference
// A - cavity then keeps every closed wall and removes the material in front of
// an open one. Building the cavity from a seed box (rather than from the input)
// is what makes an open pocket possible: an inward cut can only remove
// material, never extend the cavity past the face.
//
// A half-space that leaves nothing below it means the requested wall consumes
// the body in that direction and is reported as degenerate.
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
  bool any_open = false;
  for (const char flag : is_open) {
    if (flag != 0) any_open = true;
  }

  if (any_open) {
    std::size_t closed_faces = 0;
    for (const char flag : is_open) {
      if (flag == 0) ++closed_faces;
    }
    if (closed_faces == 0u) {
      return Error{Errc::unsupported, "shell_solid: at least one face must stay closed"};
    }
  }

  const double tolerance = options.tolerance;

  // The cavity is the intersection of one half-space per face, seeded by a box
  // that contains the solid with room to spare on every side. A closed face
  // cuts `thickness` inside the input; an open face cuts `thickness` outside it,
  // so the cavity protrudes through the opening and the difference below opens
  // the pocket instead of leaving a coincident cap.
  Point3d lo = brep.vertex(brep.vertex_ids().front()).point;
  Point3d hi = lo;
  for (const VertexId id : brep.vertex_ids()) {
    const Point3d& point = brep.vertex(id).point;
    for (int axis = 0; axis < 3; ++axis) {
      lo[axis] = std::min(lo[axis], point[axis]);
      hi[axis] = std::max(hi[axis], point[axis]);
    }
  }
  const double margin = 2.0 * options.thickness;
  Point3d seed_lo{lo[0] - margin, lo[1] - margin, lo[2] - margin};
  Point3d seed_hi{hi[0] + margin, hi[1] + margin, hi[2] + margin};
  if (!std::isfinite(seed_lo[0]) || !std::isfinite(seed_lo[1]) || !std::isfinite(seed_lo[2]) ||
      !std::isfinite(seed_hi[0]) || !std::isfinite(seed_hi[1]) || !std::isfinite(seed_hi[2])) {
    return Error{Errc::non_finite, "shell_solid: wall thickness moves the seed box out of range"};
  }
  const Result<Brep> seed = make_box(seed_lo, seed_hi, tolerance);
  if (!seed.has_value()) {
    return Error{Errc::degenerate, "shell_solid: could not build the inner seed box"};
  }

  Brep inner = seed.value();
  for (std::size_t i = 0; i < view.value().faces.size(); ++i) {
    const bool open = is_open[i] != 0;
    const Plane3 plane{view.value().outward[i],
                       view.value().offset[i] + (open ? options.thickness : -options.thickness)};

    // A plane at or below every vertex leaves no interior on its below side;
    // cutting there would either be an identity or a tangency error.
    double distance_min = std::numeric_limits<double>::infinity();
    for (const VertexId id : inner.vertex_ids()) {
      distance_min = std::min(distance_min, plane.signed_distance(inner.vertex(id).point));
    }
    if (!(distance_min < -tolerance)) {
      return Error{Errc::degenerate, "shell_solid: wall thickness leaves no interior"};
    }

    const Result<SplitSolidResult> split =
        split_solid_by_plane(inner, inner.solid_ids().front(), plane, tolerance);
    if (!split.has_value()) {
      return Error{Errc::degenerate,
                   "shell_solid: wall thickness leaves no measurable interior at this tolerance"};
    }
    if (split.value().below.solid_count() == 0u) {
      return Error{Errc::degenerate, "shell_solid: wall thickness leaves no interior"};
    }
    if (split.value().below.solid_count() != 1u) {
      return Error{Errc::invalid_topology, "shell_solid: cut split the inner solid"};
    }
    inner = std::move(split.value().below);
  }

  SolidBooleanOptions boolean_options;
  boolean_options.tolerance = tolerance;
  const Result<SolidBooleanResult> difference = boolean_solid(
      brep, inner, SolidBooleanOp::difference, boolean_options);
  if (!difference.has_value()) {
    // A wall that is only a few tolerances wide cannot be stitched back onto
    // the outer solid; report the geometric limit rather than a stitcher error.
    if (difference.error().code == Errc::non_finite) return difference.error();
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
