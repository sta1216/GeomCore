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

namespace detail {

// True when the half-spaces have a point strictly more than `tolerance` inside
// every plane. Used only after a split failure: a genuine splitter bug must not
// be reported as a wall that is too thin, while a cavity narrower than the
// numeric resolution has no such point and is honestly degenerate.
//
// A triple intersection lies exactly on its three planes, so it can never be
// the witness itself. Instead collect the vertices of the intersection that
// satisfy every other half-space, and test their centroid: a convex
// combination of feasible vertices is strictly inside the intersection, and its
// clearance is the smallest wall half-width for a slab-like cavity.
inline bool shell_half_spaces_have_interior(const std::vector<Plane3>& planes, double tolerance) {
  const std::size_t count = planes.size();
  if (count > 64u) return true;  // too many triples to certify a failure; surface it
  std::vector<Point3d> vertices;
  vertices.reserve(count);
  for (std::size_t i = 0; i < count; ++i) {
    for (std::size_t j = i + 1; j < count; ++j) {
      for (std::size_t k = j + 1; k < count; ++k) {
        const Vec3d& ni = planes[i].normal;
        const Vec3d& nj = planes[j].normal;
        const Vec3d& nk = planes[k].normal;
        const Vec3d cjk = cross(nj, nk);
        const double det = dot(ni, cjk);
        if (!std::isfinite(det) || std::abs(det) < 1e-15) continue;
        const Vec3d numerator = cjk * planes[i].offset + cross(nk, ni) * planes[j].offset +
                                cross(ni, nj) * planes[k].offset;
        const Vec3d point = numerator * (1.0 / det);
        if (!std::isfinite(point[0]) || !std::isfinite(point[1]) || !std::isfinite(point[2])) {
          continue;
        }
        const Point3d candidate{point[0], point[1], point[2]};
        bool feasible = true;
        for (const Plane3& plane : planes) {
          if (plane.signed_distance(candidate) > tolerance) {
            feasible = false;
            break;
          }
        }
        if (feasible) vertices.push_back(candidate);
      }
    }
  }
  if (vertices.empty()) return false;
  Point3d centroid = vertices.front();
  for (std::size_t i = 1; i < vertices.size(); ++i) {
    centroid = centroid + (vertices[i] - centroid) * (1.0 / static_cast<double>(i + 1));
  }
  // Equality means the cavity is exactly two tolerances wide: the assembler
  // cannot stitch that wall, so it counts as unmeasurable rather than interior.
  const double required = tolerance * (1.0 + 4.0 * std::numeric_limits<double>::epsilon());
  for (const Plane3& plane : planes) {
    if (!(plane.signed_distance(centroid) < -required)) return false;
  }
  return true;
}

}  // namespace detail

// Hollows a convex planar solid to a uniform wall of `thickness`, leaving an
// inner cavity. M4c-1 accepts closed, axis-aligned solids only: an open face
// list, an oblique face, a non-convex body, or a curved face is refused with a
// stable `Errc::unsupported` diagnostic rather than returning an order-
// dependent solid. M4c-2 lifts the open-face and oblique restrictions together
// with the splitter/boolean face-order fix.
//
// The cavity is built as a half-space intersection, not by cutting the input:
// a bounding box around the input, expanded by twice the wall thickness, is
// cut once per face, with each face plane moved inward by `thickness`. The
// final regularized difference A - cavity then leaves the uniform wall.
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
  // The margin must also stay well above the numeric resolution: a seed face
  // only `3 * thickness` away from a cut would leave a sliver that the stitcher
  // merges away for very thin walls. A relative floor keeps the seed a real
  // solid at any thickness.
  double extent = 0.0;
  for (int axis = 0; axis < 3; ++axis) extent = std::max(extent, hi[axis] - lo[axis]);
  const double margin =
      std::max({2.0 * options.thickness, 1e-6 * std::max(1.0, extent), 1000.0 * tolerance});
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
  // Sequential half-space clipping is order sensitive when an oblique plane
  // meets a face that an earlier oblique cut already trimmed. This slice of
  // M4c supports at most one oblique face; box-like parts and extrusions with a
  // single slanted face have a deterministic, well-conditioned cut order.
  // Cutting axis-aligned planes first also keeps those cuts on box faces.
  std::vector<std::size_t> cut_order(view.value().faces.size());
  for (std::size_t i = 0; i < cut_order.size(); ++i) cut_order[i] = i;
  const auto oblique = [&view](std::size_t i) {
    int non_axis = 0;
    for (int axis = 0; axis < 3; ++axis) {
      if (std::abs(view.value().outward[i][axis]) > 1e-9) ++non_axis;
    }
    return non_axis > 1;
  };
  std::stable_partition(cut_order.begin(), cut_order.end(),
                        [&oblique](std::size_t i) { return !oblique(i); });
  std::vector<Plane3> planes;
  planes.reserve(cut_order.size());
  for (const std::size_t i : cut_order) {
    const bool open = is_open[i] != 0;
    planes.push_back(Plane3{view.value().outward[i],
                            view.value().offset[i] +
                                (open ? options.thickness : -options.thickness)});
  }
  for (const Plane3& plane : planes) {
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
      if (!detail::shell_half_spaces_have_interior(planes, tolerance)) {
        return Error{Errc::degenerate,
                     "shell_solid: wall thickness leaves no measurable interior at this tolerance"};
      }
      // A real splitter failure: keep its code, but keep the public prefix.
      return Error{Errc::invalid_topology,
                   "shell_solid: cavity intersection could not be assembled"};
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
