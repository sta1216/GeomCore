#pragma once

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <numbers>
#include <optional>
#include <span>
#include <utility>
#include <vector>

#include <geom/brep/assemble_curved.hpp>
#include <geom/brep/detail/fillet_box.hpp>
#include <geom/brep/face_uv.hpp>
#include <geom/brep/topology.hpp>
#include <geom/core/error.hpp>
#include <geom/core/result.hpp>
#include <geom/curves/analytic_curve.hpp>
#include <geom/surfaces/analytic.hpp>
#include <geom/surfaces/plane.hpp>

namespace geom {
inline namespace v1 {
namespace detail {

// Bounds of the box whose twelve edges are all rounded. `lo` and `hi` are the
// per-axis minima and maxima of the eight corner vertices.
struct FilletBoxRounded {
  double lo[3];
  double hi[3];
};

namespace fillet_box_rounded_detail {

using fillet_box_detail::axis_direction;
using fillet_box_detail::cross_axis_a;
using fillet_box_detail::cross_axis_b;
using fillet_box_detail::LoopPoints;
using fillet_box_detail::LoopSpans;
using fillet_box_detail::make_face;
using fillet_box_detail::point_at;

// Inset planar face of the rounded box: the slab at `plane_axis = plane_value`
// whose other two coordinates run inside [lo + radius, hi - radius].
inline CurvedAssembledFace plane_face(const FilletBoxRounded& box, int plane_axis, int sign,
                                      double radius) {
  const int p_axis = cross_axis_a(plane_axis);
  const int q_axis = cross_axis_b(plane_axis);
  const double plane_value = sign > 0 ? box.hi[plane_axis] : box.lo[plane_axis];
  const Vec3d normal = axis_direction(plane_axis) * static_cast<double>(sign);
  const double offset = sign > 0 ? box.hi[plane_axis] : -box.lo[plane_axis];
  const double p_lo = box.lo[p_axis] + radius;
  const double p_hi = box.hi[p_axis] - radius;
  const double q_lo = box.lo[q_axis] + radius;
  const double q_hi = box.hi[q_axis] - radius;

  LoopPoints points{point_at(plane_axis, plane_value, p_axis, p_lo, q_axis, q_lo),
                    point_at(plane_axis, plane_value, p_axis, p_hi, q_axis, q_lo),
                    point_at(plane_axis, plane_value, p_axis, p_hi, q_axis, q_hi),
                    point_at(plane_axis, plane_value, p_axis, p_lo, q_axis, q_hi)};
  if (dot(cross(axis_direction(p_axis), axis_direction(q_axis)), normal) < 0.0) {
    std::reverse(points.begin() + 1, points.end());
  }

  std::vector<Point2d> uv;
  uv.reserve(points.size());
  for (const Point3d& point : points) {
    uv.push_back(Point2d{point[p_axis], point[q_axis]});
  }
  return make_face(SurfaceGeometry{Plane3{normal, offset}}, std::move(points), LoopSpans{},
                   std::move(uv));
}

// Quarter-cylinder patch at one rounded cross-section corner. The cylinder
// spans [lo + radius, hi - radius] along the family axis, so both end arcs meet
// the corner octants exactly.
inline CurvedAssembledFace cylinder_face(const FilletBoxRounded& box, int family_axis, int bits_a,
                                         int bits_b, double radius) {
  const int axis_a = cross_axis_a(family_axis);
  const int axis_b = cross_axis_b(family_axis);
  const Vec3d direction_a = axis_direction(axis_a);
  const Vec3d direction_f = axis_direction(family_axis);
  const Vec3d frame_b = cross(direction_f, direction_a);
  const double parity = dot(frame_b, axis_direction(axis_b)) > 0.0 ? 1.0 : -1.0;
  const double sign_a = bits_a != 0 ? 1.0 : -1.0;
  const double sign_b = bits_b != 0 ? 1.0 : -1.0;
  const double middle = std::atan2(sign_b * parity, sign_a);
  const double t0 = middle - std::numbers::pi / 4.0;
  const double t1 = middle + std::numbers::pi / 4.0;
  const double height = box.hi[family_axis] - box.lo[family_axis] - 2.0 * radius;

  Point3d origin{};
  origin[axis_a] = bits_a != 0 ? box.hi[axis_a] - radius : box.lo[axis_a] + radius;
  origin[axis_b] = bits_b != 0 ? box.hi[axis_b] - radius : box.lo[axis_b] + radius;
  origin[family_axis] = box.lo[family_axis] + radius;

  const auto evaluate = [&](double angle, double height_value) {
    const Vec3d radial = direction_a * std::cos(angle) + frame_b * std::sin(angle);
    return origin + direction_f * height_value + radial * radius;
  };

  LoopPoints points{evaluate(t0, 0.0), evaluate(t1, 0.0), evaluate(t1, height),
                    evaluate(t0, height)};
  const CircleCurve3 bottom{origin, direction_a, frame_b, radius};
  const CircleCurve3 top{origin + direction_f * height, direction_a, frame_b, radius};
  LoopSpans spans{std::optional<CurvedArcSpan>{CurvedArcSpan{bottom, t0, t1}}, std::nullopt,
                  std::optional<CurvedArcSpan>{CurvedArcSpan{top, t0, t1}}, std::nullopt};
  std::vector<Point2d> uv{Point2d{t0, 0.0}, Point2d{t1, 0.0}, Point2d{t1, height},
                          Point2d{t0, height}};

  return make_face(SurfaceGeometry{Cylinder3{origin, direction_f, direction_a, radius}},
                   std::move(points), std::move(spans), std::move(uv));
}

// Spherical octant at one box corner. The chart axis points from the sphere
// centre to the corner tangency point, which makes that point the chart pole;
// the three great-circle arcs then project to an exact quarter lune, and the
// loop is ordered so the pole closes the u band without a seam jump.
inline std::optional<CurvedAssembledFace> sphere_face(const FilletBoxRounded& box, int bits_x,
                                                      int bits_y, int bits_z, double radius) {
  const Point3d center{bits_x != 0 ? box.hi[0] - radius : box.lo[0] + radius,
                       bits_y != 0 ? box.hi[1] - radius : box.lo[1] + radius,
                       bits_z != 0 ? box.hi[2] - radius : box.lo[2] + radius};
  const double sign[3] = {bits_x != 0 ? 1.0 : -1.0, bits_y != 0 ? 1.0 : -1.0,
                          bits_z != 0 ? 1.0 : -1.0};
  const Vec3d axis = axis_direction(0) * sign[0];
  const Vec3d x_axis = axis_direction(1) * sign[1];
  const Sphere3 surface{center, axis, x_axis, radius};

  const Vec3d directions[3] = {axis, x_axis, axis_direction(2) * sign[2]};
  const Point3d corners[3] = {center + directions[0] * radius, center + directions[1] * radius,
                              center + directions[2] * radius};
  const bool forward = sign[0] * sign[1] * sign[2] > 0.0;
  const int order[3] = {forward ? 0 : 2, 1, forward ? 2 : 0};

  LoopPoints points;
  LoopSpans spans;
  points.reserve(3);
  spans.reserve(3);
  for (int index = 0; index < 3; ++index) {
    points.push_back(corners[order[index]]);
    const CircleCurve3 arc{center, directions[order[index]], directions[order[(index + 1) % 3]],
                           radius};
    spans.push_back(CurvedArcSpan{arc, 0.0, std::numbers::pi / 2.0});
  }

  std::vector<Point2d> uv;
  uv.reserve(points.size());
  double previous_u = 0.0;
  bool has_previous_u = false;
  for (const Point3d& point : points) {
    const Result<Point2d> projected =
        face_uv_project_sphere(surface, point, previous_u, has_previous_u);
    if (!projected.has_value()) return std::nullopt;
    uv.push_back(projected.value());
    previous_u = projected.value()[0];
    has_previous_u = true;
  }

  return make_face(SurfaceGeometry{surface}, std::move(points), std::move(spans), std::move(uv));
}

}  // namespace fillet_box_rounded_detail

// True only when `edges`, after deduplication, are exactly the twelve live
// edges of `box`, in any order.
inline bool fillet_box_all_edges(const Brep& brep, std::span<const EdgeId> edges,
                                 const FilletBox& box, double tolerance) {
  double scale = 1.0;
  for (int axis = 0; axis < 3; ++axis) {
    scale = std::max(scale, std::abs(box.lo[axis]));
    scale = std::max(scale, std::abs(box.hi[axis]));
  }
  const double slack = std::max(64.0 * tolerance, 1e-12 * scale);

  std::vector<EdgeId> unique;
  unique.reserve(edges.size());
  for (const EdgeId id : edges) {
    if (!brep.edge_alive(id)) return false;
    if (std::find(unique.begin(), unique.end(), id) == unique.end()) unique.push_back(id);
  }
  if (unique.size() != 12) return false;

  const auto is_corner = [&](const Point3d& point) {
    for (int axis = 0; axis < 3; ++axis) {
      if (std::abs(point[axis] - box.lo[axis]) > slack &&
          std::abs(point[axis] - box.hi[axis]) > slack) {
        return false;
      }
    }
    return true;
  };

  std::vector<EdgeId> box_edges;
  box_edges.reserve(12);
  const std::vector<EdgeId> all = brep.edge_ids();
  for (const EdgeId id : all) {
    const Edge& edge = brep.edge(id);
    if (is_corner(brep.vertex(edge.v0).point) && is_corner(brep.vertex(edge.v1).point)) {
      box_edges.push_back(id);
    }
  }
  if (box_edges.size() != 12) return false;

  for (const EdgeId id : unique) {
    if (std::find(box_edges.begin(), box_edges.end(), id) == box_edges.end()) return false;
  }
  return true;
}

// Rounds all twelve edges of a box with six inset planar faces, twelve quarter
// cylinders and eight spherical corner octants.
inline Result<Brep> fillet_box_rounded(const FilletBoxRounded& box, double radius,
                                       double tolerance) {
  // At half the smallest dimension the inset planes and the cylinder heights
  // collapse to zero, so the twenty-six patches can no longer be assembled.
  const double smallest = std::min({box.hi[0] - box.lo[0], box.hi[1] - box.lo[1],
                                    box.hi[2] - box.lo[2]});
  if (!(radius < 0.5 * smallest)) {
    return Error{Errc::unsupported,
                 "fillet_edges: radius must be less than half the smallest box dimension", 0u, 0u};
  }

  std::vector<CurvedAssembledFace> faces;
  faces.reserve(26);

  for (int plane_axis = 0; plane_axis < 3; ++plane_axis) {
    for (int sign = -1; sign <= 1; sign += 2) {
      faces.push_back(fillet_box_rounded_detail::plane_face(box, plane_axis, sign, radius));
    }
  }
  for (int family_axis = 0; family_axis < 3; ++family_axis) {
    for (int bits_a = 0; bits_a < 2; ++bits_a) {
      for (int bits_b = 0; bits_b < 2; ++bits_b) {
        faces.push_back(
            fillet_box_rounded_detail::cylinder_face(box, family_axis, bits_a, bits_b, radius));
      }
    }
  }
  for (int bits_x = 0; bits_x < 2; ++bits_x) {
    for (int bits_y = 0; bits_y < 2; ++bits_y) {
      for (int bits_z = 0; bits_z < 2; ++bits_z) {
        std::optional<CurvedAssembledFace> corner =
            fillet_box_rounded_detail::sphere_face(box, bits_x, bits_y, bits_z, radius);
        if (!corner.has_value()) {
          return Error{Errc::non_finite, "fillet_edges: rounded box corner patch is not finite"};
        }
        faces.push_back(std::move(corner.value()));
      }
    }
  }

  for (std::size_t index = 0; index < faces.size(); ++index) {
    for (const std::vector<Point3d>& loop : faces[index].loops) {
      if (loop.size() < 3) {
        return Error{Errc::invalid_topology, "fillet_edges: rounded box face has too few points",
                     static_cast<std::uint32_t>(index)};
      }
      for (const Point3d& point : loop) {
        if (!std::isfinite(point[0]) || !std::isfinite(point[1]) || !std::isfinite(point[2])) {
          return Error{Errc::non_finite, "fillet_edges: rounded box face is not finite",
                       static_cast<std::uint32_t>(index)};
        }
      }
    }
  }

  Result<Brep> assembled = assemble_curved_solid(faces, tolerance);
  if (!assembled.has_value()) {
    // The public API must not leak the assembler's internal diagnostic.
    return Error{Errc::invalid_topology, "fillet_edges: fillet could not be assembled"};
  }
  return assembled;
}

}  // namespace detail
}  // namespace v1
}  // namespace geom
