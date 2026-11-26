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
#include <geom/brep/detail/half_space.hpp>
#include <geom/brep/mass_properties.hpp>
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
  std::uint32_t cavity_shells = 0;  // enclosed pockets; 0 when a face is open
};

// Hollows a convex planar solid to a uniform wall of `thickness` at every face
// that is not listed in `open_faces`. A listed face is left unwalled: its
// cavity half-space is relaxed outward past the input surface, so the pocket
// breaks through and the result is an open shell with no enclosed pocket.
// M4c-2 supports any closed convex planar solid, including oblique faces, and
// any non-empty proper subset of open faces. Duplicate open-face ids are
// accepted once; opening every face is refused because no wall would remain.
//
// A closed solid whose faces are all axis-aligned is a box, and intersecting
// its half-spaces leaves one interval per axis. That case keeps the direct
// inset-box path: the cavity is the input AABB inset by `thickness` on every
// side and no boolean splitter is involved. Every other input, including a
// closed oblique solid and any open shell, uses the general path: the
// inward-offset face half-spaces are intersected with a bounding seed box (the
// input AABB expanded by at least the wall thickness), which keeps the cavity
// bounded when an open face is relaxed outward. The final regularized
// difference A - cavity leaves the uniform wall.
//
// `cavity_shells` is the number of enclosed pockets after the cut and is 0
// when any face is open.
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
  const detail::BlendPlanarSolid& solid = view.value();
  const double tolerance = options.tolerance;

  // Requested open faces are validated against the single shell of the input,
  // then deduplicated: a repeat is accepted once. Opening every face leaves
  // nothing to wall and is refused before any geometry is built.
  std::vector<std::size_t> open_indices;
  open_indices.reserve(open_faces.size());
  for (const FaceId face : open_faces) {
    if (!brep.face_alive(face)) {
      return Error{Errc::index_out_of_range,
                   "shell_solid: open face does not belong to the solid"};
    }
    const auto found = std::find(solid.faces.begin(), solid.faces.end(), face);
    if (found == solid.faces.end()) {
      return Error{Errc::index_out_of_range,
                   "shell_solid: open face does not belong to the solid"};
    }
    open_indices.push_back(static_cast<std::size_t>(found - solid.faces.begin()));
  }
  std::sort(open_indices.begin(), open_indices.end());
  open_indices.erase(std::unique(open_indices.begin(), open_indices.end()),
                     open_indices.end());

  const bool any_open = !open_indices.empty();
  if (any_open && open_indices.size() == solid.faces.size()) {
    return Error{Errc::unsupported,
                 "shell_solid: at least one face must stay closed"};
  }

  // The direct inset-box path is only valid for a closed solid whose unit face
  // normals are coordinate axes. The slack is angular (a dimensionless
  // deviation between unit normals), not a distance: a user tolerance must not
  // widen it, or a rotated solid at a coarse tolerance would be misclassified
  // as axis-aligned and hollowed with the wrong AABB-inset cavity.
  const double axis_slack = 1e-12;
  const auto axis_aligned = [axis_slack](const Vec3d& normal) {
    int dominant = -1;
    for (int axis = 0; axis < 3; ++axis) {
      if (std::fabs(normal[axis]) > axis_slack) {
        if (dominant >= 0) return false;
        dominant = axis;
      }
    }
    return dominant >= 0 &&
           std::fabs(std::fabs(normal[dominant]) - 1.0) <= axis_slack;
  };
  bool all_axis_aligned = true;
  for (const Vec3d& outward : solid.outward) {
    if (!axis_aligned(outward)) {
      all_axis_aligned = false;
      break;
    }
  }

  // A wall that is only a few tolerances wide cannot be measured or stitched
  // back onto the outer solid; both paths report the same geometric limit.
  const char* const no_interior =
      "shell_solid: wall thickness leaves no measurable interior at this tolerance";

  const auto build_cavity = [&]() -> Result<Brep> {
    if (!any_open && all_axis_aligned) {
      // Every accepted face is axis-aligned, so the cavity is the input
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
          return Error{Errc::degenerate,
                       "shell_solid: wall thickness leaves no interior"};
        }
      }
      const Result<Brep> inner_box = make_box(inner_lo, inner_hi, tolerance);
      if (!inner_box.has_value()) {
        return Error{Errc::degenerate, no_interior};
      }
      return inner_box.value();
    }

    // General path. The seed box is the input AABB expanded on every axis by
    // a margin of at least the wall thickness: a relaxed open face must push
    // the pocket past the input surface, and the seed box keeps it bounded.
    Point3d lo = brep.vertex(brep.vertex_ids().front()).point;
    Point3d hi = lo;
    for (const VertexId id : brep.vertex_ids()) {
      const Point3d& point = brep.vertex(id).point;
      for (int axis = 0; axis < 3; ++axis) {
        lo[axis] = std::min(lo[axis], point[axis]);
        hi[axis] = std::max(hi[axis], point[axis]);
      }
    }
    double extent = 0.0;
    for (int axis = 0; axis < 3; ++axis) {
      extent = std::max(extent, hi[axis] - lo[axis]);
    }
    const double margin = std::max({options.thickness, 1e-6 * std::max(1.0, extent),
                                    1000.0 * tolerance});

    std::vector<detail::HalfSpacePlane> planes;
    planes.reserve(6u + solid.faces.size());
    for (int axis = 0; axis < 3; ++axis) {
      Vec3d positive{};
      positive[axis] = 1.0;
      planes.push_back(detail::HalfSpacePlane{positive, hi[axis] + margin});
      Vec3d negative{};
      negative[axis] = -1.0;
      planes.push_back(detail::HalfSpacePlane{negative, -(lo[axis] - margin)});
    }

    // A relaxed open face at an input extreme can land exactly on the matching
    // seed-box plane when the margin equals the thickness. Coincident
    // half-spaces would build the same face twice and fail assembly, so the
    // seed box keeps the plane in that case.
    const double plane_slack = std::max(1000.0 * tolerance, 1e-12 * solid.scale);
    for (std::size_t i = 0; i < solid.faces.size(); ++i) {
      const bool is_open =
          std::binary_search(open_indices.begin(), open_indices.end(), i);
      const Vec3d& normal = solid.outward[i];
      const double offset =
          solid.offset[i] + (is_open ? options.thickness : -options.thickness);
      const bool duplicate = std::any_of(
          planes.begin(), planes.end(),
          [&](const detail::HalfSpacePlane& existing) {
            return dot(existing.normal, normal) >= 1.0 - 1e-12 &&
                   std::fabs(existing.offset - offset) <= plane_slack;
          });
      if (!duplicate) {
        planes.push_back(detail::HalfSpacePlane{normal, offset});
      }
    }

    const Result<Brep> cavity = detail::half_space_solid(planes, tolerance);
    if (!cavity.has_value()) {
      if (cavity.error().code == Errc::non_finite) return cavity.error();
      return Error{Errc::degenerate, no_interior};
    }
    return cavity.value();
  };

  const Result<Brep> cavity = build_cavity();
  if (!cavity.has_value()) return cavity.error();

  SolidBooleanOptions boolean_options;
  boolean_options.tolerance = tolerance;
  const Result<SolidBooleanResult> difference = boolean_solid(
      brep, cavity.value(), SolidBooleanOp::difference, boolean_options);
  if (!difference.has_value()) {
    // A wall that is only a few tolerances wide cannot be stitched back onto
    // the outer solid; report the geometric limit rather than a stitcher error.
    if (difference.error().code == Errc::non_finite) return difference.error();
    return Error{Errc::degenerate, no_interior};
  }

  Brep shelled = std::move(difference.value().brep);
  if (shelled.solid_count() != 1u) {
    return Error{Errc::invalid_topology, "shell_solid: difference did not leave one solid"};
  }
  if (!validate_topology(shelled).ok()) {
    return Error{Errc::invalid_topology, "shell_solid: shell produced an invalid solid"};
  }

  if (any_open) {
    // An open pocket that never enters the solid leaves the input unchanged,
    // which a successful boolean can still do when the relaxed cavity lies
    // wholly outside the input. Require that material was actually removed
    // before reporting an open shell; this is also the geometric limit of an
    // over-thick wall. The floor is a fraction of the input volume, not an
    // absolute volume: a 1e-5-scale model must keep the same fate as its unit
    // counterpart.
    const Result<MassProperties> before = mass_properties(brep);
    if (!before.has_value()) return before.error();
    const Result<MassProperties> after = mass_properties(shelled);
    if (!after.has_value()) return after.error();
    const double removed = before.value().volume - after.value().volume;
    if (removed <= 1e-12 * std::fabs(before.value().volume)) {
      return Error{Errc::degenerate, no_interior};
    }
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
