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
// The inner solid is the input cut by one plane per shell face: a closed face
// is moved inward by `thickness` and the below side is kept. A face listed in
// `open_faces` is instead moved outward: the plane is then strictly outside the
// convex input, so keeping the below side is the identity and the cut is
// skipped; the final difference removes that wall, which opens the pocket.
// Either way the inner solid is never conflated with the outer one, and the
// regularized difference A - B only depends on B inside A.
//
// The cuts are derived from the input, then applied to a working copy one at a
// time. A plane that would leave no interior below it means the requested wall
// consumes the body in that direction and is reported as degenerate.
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

  const double tolerance = options.tolerance;
  Brep inner = brep;
  for (std::size_t i = 0; i < view.value().faces.size(); ++i) {
    const bool open = is_open[i] != 0;
    const Plane3 plane{view.value().outward[i],
                       view.value().offset[i] + (open ? options.thickness : -options.thickness)};

    // The below side has positive volume only when some vertex is strictly
    // below the plane. A plane that is at or above every vertex consumes the
    // whole body, i.e. the requested wall leaves no interior.
    double distance_min = std::numeric_limits<double>::infinity();
    double distance_max = -std::numeric_limits<double>::infinity();
    for (const VertexId id : inner.vertex_ids()) {
      const double distance = plane.signed_distance(inner.vertex(id).point);
      distance_min = std::min(distance_min, distance);
      distance_max = std::max(distance_max, distance);
    }
    if (!(distance_min < -tolerance)) {
      return Error{Errc::degenerate, "shell_solid: wall thickness leaves no interior"};
    }
    if (!(distance_max > tolerance)) {
      // The plane does not split the solid: keeping the below side leaves the
      // working solid unchanged. This is the open-face case, whose outward
      // plane is strictly outside a convex input.
      continue;
    }

    const Result<SplitSolidResult> split = split_solid_by_plane(
        inner, inner.solid_ids().front(), plane, tolerance);
    if (!split.has_value()) return split.error();
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
  if (!difference.has_value()) return difference.error();

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
