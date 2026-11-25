#pragma once

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <span>
#include <utility>
#include <vector>

#include <geom/brep/detail/blend_common.hpp>
#include <geom/brep/split_solid.hpp>
#include <geom/brep/validate.hpp>
#include <geom/core/error.hpp>
#include <geom/core/result.hpp>
#include <geom/surfaces/plane.hpp>

namespace geom {
inline namespace v1 {

struct ChamferOptions {
  double distance = 0.0;  // setback along each adjacent face, > 0
  double tolerance = 1e-9;
};

struct ChamferResult {
  Brep brep;
  std::uint32_t chamfered_edges = 0;  // distinct edges actually cut
};

// Replaces each selected straight edge of a convex planar solid with a planar
// facet: the edge and the strips of width `distance` on its two adjacent faces
// are removed, leaving a symmetric chamfer whose cross-section is a right
// triangle (leg `distance`) swept along the edge.
//
// The cuts are derived from the input solid (one plane per distinct edge, in
// first-occurrence order) and then applied to a working copy one at a time, so
// edges that share a vertex each subtract only their own share. chamfered_edges
// counts distinct edges; a repeated EdgeId is cut once.
//
// Diagnostics are static texts because Error::message is a const char*.
inline Result<ChamferResult> chamfer_edges(const Brep& brep, std::span<const EdgeId> edges,
                                           const ChamferOptions& options = {}) {
  if (!std::isfinite(options.distance)) {
    return Error{Errc::non_finite, "chamfer_edges: distance must be finite"};
  }
  if (!(options.distance > 0.0)) {
    return Error{Errc::degenerate, "chamfer_edges: distance must be positive"};
  }
  if (!std::isfinite(options.tolerance)) {
    return Error{Errc::non_finite, "chamfer_edges: tolerance must be finite"};
  }
  if (!(options.tolerance > 0.0)) {
    return Error{Errc::tolerance_violation, "chamfer_edges: tolerance must be positive"};
  }

  // Shape, planarity, closedness and convexity are properties of the whole
  // solid and are checked before any individual edge is looked at.
  const Result<detail::BlendPlanarSolid> view =
      detail::blend_planar_convex_solid(brep, options.tolerance, "chamfer_edges");
  if (!view.has_value()) return view.error();

  if (edges.empty()) {
    return Error{Errc::unsupported, "chamfer_edges: at least one edge is required"};
  }
  std::vector<EdgeId> distinct;
  distinct.reserve(edges.size());
  for (const EdgeId edge : edges) {
    if (std::find(distinct.begin(), distinct.end(), edge) == distinct.end()) {
      distinct.push_back(edge);
    }
  }

  // Derive every cut plane from the original solid before cutting anything, so
  // a later edge's width or topology is judged against the input, not against
  // an already chamfered neighbour.
  std::vector<Plane3> planes;
  planes.reserve(distinct.size());
  for (const EdgeId edge : distinct) {
    if (!brep.edge_alive(edge)) {
      return Error{Errc::index_out_of_range, "chamfer_edges: edge does not belong to the solid"};
    }
    const Result<Vec3d> direction =
        detail::blend_edge_direction(brep, edge, "chamfer_edges");
    if (!direction.has_value()) return direction.error();

    std::size_t coedge_count = 0;
    std::vector<FaceId> adjacent;
    for (const CoedgeId coedge_id : brep.coedge_ids()) {
      const Coedge& coedge = brep.coedge(coedge_id);
      if (coedge.edge != edge) continue;
      ++coedge_count;
      if (std::find(adjacent.begin(), adjacent.end(), coedge.face) == adjacent.end()) {
        adjacent.push_back(coedge.face);
      }
    }
    if (coedge_count != 2u) {
      return Error{Errc::invalid_topology,
                   "chamfer_edges: edge must have exactly two coedges"};
    }
    if (adjacent.size() != 2u) {
      return Error{Errc::invalid_topology,
                   "chamfer_edges: edge must border two distinct faces"};
    }

    const Edge& edge_data = brep.edge(edge);
    const Vec3d start = detail::blend_as_vec(brep.vertex(edge_data.v0).point);
    const Vec3d end = detail::blend_as_vec(brep.vertex(edge_data.v1).point);
    const Vec3d midpoint = (start + end) * 0.5;

    const Result<Vec3d> inward0 = detail::blend_inward_direction(
        brep, adjacent[0], edge, direction.value(), midpoint, "chamfer_edges");
    if (!inward0.has_value()) return inward0.error();
    const Result<Vec3d> inward1 = detail::blend_inward_direction(
        brep, adjacent[1], edge, direction.value(), midpoint, "chamfer_edges");
    if (!inward1.has_value()) return inward1.error();

    for (const FaceId face : adjacent) {
      const Result<Vec3d> inward =
          face == adjacent[0] ? inward0 : inward1;
      const Result<double> width =
          detail::blend_face_width(brep, face, midpoint, inward.value(), "chamfer_edges");
      if (!width.has_value()) return width.error();
      if (!(options.distance < width.value())) {
        return Error{Errc::unsupported,
                     "chamfer_edges: distance exceeds the adjacent face width"};
      }
    }

    const Vec3d sum = inward0.value() + inward1.value();
    const double sum_len = std::hypot(sum[0], sum[1], sum[2]);
    if (!(sum_len > 1e-12) || !std::isfinite(sum_len)) {
      return Error{Errc::unsupported,
                   "chamfer_edges: adjacent faces are coplanar"};
    }
    const Vec3d normal{sum[0] / sum_len, sum[1] / sum_len, sum[2] / sum_len};
    const Vec3d setback = midpoint + inward0.value() * options.distance;
    planes.push_back(Plane3{normal, dot(normal, setback)});
  }

  Brep current = brep;
  for (const Plane3& plane : planes) {
    const Result<SplitSolidResult> split = split_solid_by_plane(
        current, current.solid_ids().front(), plane, options.tolerance);
    if (!split.has_value()) return split.error();
    if (split.value().above.solid_count() != 1u) {
      return Error{Errc::invalid_topology, "chamfer_edges: cut split the solid"};
    }
    // The edge lies on the negative side of its cut plane; the material to
    // keep is the positive side.
    current = std::move(split.value().above);
  }

  if (!validate_topology(current).ok()) {
    return Error{Errc::invalid_topology, "chamfer_edges: chamfer produced an invalid solid"};
  }

  ChamferResult result;
  result.brep = std::move(current);
  result.chamfered_edges = static_cast<std::uint32_t>(planes.size());
  return result;
}

}  // namespace v1
}  // namespace geom
