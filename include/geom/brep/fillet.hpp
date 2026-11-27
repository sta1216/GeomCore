#pragma once

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <span>
#include <utility>
#include <vector>

#include <geom/brep/detail/fillet_box.hpp>
#include <geom/brep/detail/fillet_prism.hpp>
#include <geom/brep/topology.hpp>
#include <geom/core/error.hpp>
#include <geom/core/result.hpp>

namespace geom {
inline namespace v1 {

struct FilletOptions {
  double radius = 0.0;
  double tolerance = 1e-9;
};

struct FilletResult {
  Brep brep;
  std::uint32_t filleted_edges = 0;
};

namespace detail {

// Box body of fillet_edges(), kept in one place so the box behaviour and its
// messages stay exactly as they were before prism support was added.
inline Result<FilletResult> fillet_box_selection(const Brep& brep,
                                                 std::span<const EdgeId> edges,
                                                 const FilletBox& box, double radius,
                                                 double tolerance) {
  double scale = 1.0;
  for (int axis = 0; axis < 3; ++axis) {
    scale = std::max(scale, std::abs(box.lo[axis]));
    scale = std::max(scale, std::abs(box.hi[axis]));
  }
  const double slack = std::max(64.0 * tolerance, 1e-12 * scale);

  std::vector<EdgeId> unique;
  unique.reserve(edges.size());
  for (const EdgeId id : edges) {
    if (std::find(unique.begin(), unique.end(), id) == unique.end()) unique.push_back(id);
  }

  int family_axis = -1;
  std::vector<int> corners;
  corners.reserve(unique.size());
  for (const EdgeId id : unique) {
    if (!brep.edge_alive(id)) {
      return Error{Errc::index_out_of_range, "fillet_edges: edge does not belong to the solid",
                   id.index()};
    }
    const Edge& edge = brep.edge(id);
    const Point3d& start = brep.vertex(edge.v0).point;
    const Point3d& end = brep.vertex(edge.v1).point;

    int direction_axis = -1;
    int differing = 0;
    for (int axis = 0; axis < 3; ++axis) {
      if (std::abs(end[axis] - start[axis]) > slack) {
        direction_axis = axis;
        ++differing;
      }
    }
    if (differing != 1) {
      return Error{Errc::unsupported, "fillet_edges: only axis-aligned box edges are supported",
                   id.index()};
    }
    if (family_axis < 0) {
      family_axis = direction_axis;
    } else if (family_axis != direction_axis) {
      return Error{Errc::unsupported, "fillet_edges: selected edges must be parallel",
                   id.index()};
    }

    const double start_family = start[family_axis];
    const double end_family = end[family_axis];
    const bool start_low = std::abs(start_family - box.lo[family_axis]) <= slack;
    const bool start_high = std::abs(start_family - box.hi[family_axis]) <= slack;
    const bool end_low = std::abs(end_family - box.lo[family_axis]) <= slack;
    const bool end_high = std::abs(end_family - box.hi[family_axis]) <= slack;
    if (!((start_low && end_high) || (start_high && end_low))) {
      return Error{Errc::unsupported, "fillet_edges: edge is not a box corner edge",
                   id.index()};
    }

    const Point3d& base =
        std::abs(start_family - box.lo[family_axis]) <=
                std::abs(end_family - box.lo[family_axis])
            ? start
            : end;
    const int axis_a = family_axis == 0 ? 1 : 0;
    const int axis_b = family_axis == 2 ? 1 : 2;
    const double coordinate_a = base[axis_a];
    const double coordinate_b = base[axis_b];
    if (std::abs(coordinate_a - box.lo[axis_a]) > slack &&
        std::abs(coordinate_a - box.hi[axis_a]) > slack) {
      return Error{Errc::unsupported, "fillet_edges: edge is not a box corner edge",
                   id.index()};
    }
    if (std::abs(coordinate_b - box.lo[axis_b]) > slack &&
        std::abs(coordinate_b - box.hi[axis_b]) > slack) {
      return Error{Errc::unsupported, "fillet_edges: edge is not a box corner edge",
                   id.index()};
    }
    const bool bit_a =
        std::abs(coordinate_a - box.hi[axis_a]) < std::abs(coordinate_a - box.lo[axis_a]);
    const bool bit_b =
        std::abs(coordinate_b - box.hi[axis_b]) < std::abs(coordinate_b - box.lo[axis_b]);
    corners.push_back((bit_a ? 1 : 0) | (bit_b ? 2 : 0));
  }

  for (std::size_t first = 0; first < unique.size(); ++first) {
    const Edge& left = brep.edge(unique[first]);
    for (std::size_t second = first + 1; second < unique.size(); ++second) {
      const Edge& right = brep.edge(unique[second]);
      if (left.v0 == right.v0 || left.v0 == right.v1 || left.v1 == right.v0 ||
          left.v1 == right.v1) {
        return Error{Errc::unsupported,
                     "fillet_edges: filleted edges must not share a vertex"};
      }
    }
  }

  Result<Brep> filleted = fillet_box_edges(box, corners, family_axis, radius, tolerance);
  if (!filleted.has_value()) return filleted.error();

  FilletResult result;
  result.brep = std::move(filleted.value());
  result.filleted_edges = static_cast<std::uint32_t>(corners.size());
  return result;
}

}  // namespace detail

// Rounds the selected edges of an axis-aligned box or of a convex prism with
// exact cylinder patches. Box edges must belong to one parallel family; prism
// edges must be vertical.
inline Result<FilletResult> fillet_edges(const Brep& brep, std::span<const EdgeId> edges,
                                         const FilletOptions& options = {}) {
  if (edges.empty()) {
    return Error{Errc::unsupported, "fillet_edges: at least one edge is required"};
  }
  if (!std::isfinite(options.radius)) {
    return Error{Errc::non_finite, "fillet_edges: radius must be finite"};
  }
  if (!(options.radius > 0.0)) {
    return Error{Errc::degenerate, "fillet_edges: radius must be positive"};
  }
  if (!std::isfinite(options.tolerance)) {
    return Error{Errc::non_finite, "fillet_edges: tolerance must be finite"};
  }
  if (!(options.tolerance > 0.0)) {
    return Error{Errc::tolerance_violation, "fillet_edges: tolerance must be positive"};
  }

  const Result<detail::FilletBox> box_bounds = detail::fillet_box_bounds(brep, options.tolerance);
  if (box_bounds.has_value()) {
    return detail::fillet_box_selection(brep, edges, box_bounds.value(), options.radius,
                                        options.tolerance);
  }

  const Result<detail::FilletPrism> prism = detail::fillet_prism_bounds(brep, options.tolerance);
  if (!prism.has_value()) {
    return Error{Errc::unsupported, "fillet_edges: only a box or a convex prism is supported"};
  }
  const detail::FilletPrism& shape = prism.value();

  std::vector<std::uint32_t> selected;
  selected.reserve(edges.size());
  for (const EdgeId id : edges) {
    if (!brep.edge_alive(id)) {
      return Error{Errc::index_out_of_range, "fillet_edges: edge does not belong to the solid",
                   id.index()};
    }
    const auto match = std::find(shape.vertical.begin(), shape.vertical.end(), id);
    if (match == shape.vertical.end()) {
      return Error{Errc::unsupported, "fillet_edges: only vertical prism edges are supported",
                   id.index()};
    }
    const std::uint32_t index = static_cast<std::uint32_t>(match - shape.vertical.begin());
    if (std::find(selected.begin(), selected.end(), index) == selected.end()) {
      selected.push_back(index);
    }
  }

  Result<Brep> filleted =
      detail::fillet_prism_edges(shape, selected, options.radius, options.tolerance);
  if (!filleted.has_value()) return filleted.error();

  FilletResult result;
  result.brep = std::move(filleted.value());
  result.filleted_edges = static_cast<std::uint32_t>(selected.size());
  return result;
}

}  // namespace v1
}  // namespace geom
