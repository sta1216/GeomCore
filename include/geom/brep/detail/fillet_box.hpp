#pragma once

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <numbers>
#include <optional>
#include <span>
#include <utility>
#include <vector>

#include <geom/brep/assemble_curved.hpp>
#include <geom/brep/detail/blend_common.hpp>
#include <geom/brep/topology.hpp>
#include <geom/brep/trim_face.hpp>
#include <geom/brep/validate.hpp>
#include <geom/core/error.hpp>
#include <geom/core/result.hpp>
#include <geom/curves/analytic_curve.hpp>
#include <geom/surfaces/analytic.hpp>
#include <geom/surfaces/plane.hpp>

namespace geom {
inline namespace v1 {
namespace detail {

// Axis-aligned bounds of the box a fillet selection belongs to. `lo` and `hi`
// are the per-axis minima and maxima of the eight corner vertices.
struct FilletBox {
  double lo[3];
  double hi[3];
};

namespace fillet_box_detail {

using LoopPoints = std::vector<Point3d>;
using LoopSpans = std::vector<std::optional<CurvedArcSpan>>;

inline Vec3d axis_direction(int axis) {
  Vec3d direction{0.0, 0.0, 0.0};
  direction[axis] = 1.0;
  return direction;
}

inline Vec3d unit_vector(Vec3d value) {
  const double length = std::hypot(value[0], value[1], value[2]);
  if (!(length > 0.0)) return Vec3d{};
  return value * (1.0 / length);
}

inline Point3d point_at(int axis_a, double value_a, int axis_b, double value_b, int axis_c,
                        double value_c) {
  Point3d point{};
  point[axis_a] = value_a;
  point[axis_b] = value_b;
  point[axis_c] = value_c;
  return point;
}

inline CurvedAssembledFace make_face(SurfaceGeometry surface, LoopPoints points, LoopSpans spans,
                                     std::vector<Point2d> uv) {
  CurvedAssembledFace face;
  face.surface = std::move(surface);
  face.loops.push_back(std::move(points));
  face.loops_uv.push_back(std::move(uv));
  if (!spans.empty()) face.arc_spans.push_back(std::move(spans));
  return face;
}

// The two cross-section axes of a family, in ascending axis order.
inline int cross_axis_a(int family_axis) { return family_axis == 0 ? 1 : 0; }
inline int cross_axis_b(int family_axis) { return family_axis == 2 ? 1 : 2; }

// Cap of the box. The ring is first built counter-clockwise about +e_f; the low
// cap is then reversed so every cap loop runs counter-clockwise about its own
// outward normal.
inline CurvedAssembledFace cap_face(const FilletBox& box, int family_axis,
                                    const std::array<std::array<bool, 2>, 2>& selected,
                                    bool at_high, double radius) {
  const int axis_a = cross_axis_a(family_axis);
  const int axis_b = cross_axis_b(family_axis);
  const Vec3d direction_a = axis_direction(axis_a);
  const Vec3d direction_b = axis_direction(axis_b);
  const Vec3d direction_f = axis_direction(family_axis);
  const double cap = at_high ? box.hi[family_axis] : box.lo[family_axis];

  std::array<std::array<int, 2>, 4> order{};
  if (dot(cross(direction_a, direction_b), direction_f) > 0.0) {
    order[0] = {0, 0};
    order[1] = {1, 0};
    order[2] = {1, 1};
    order[3] = {0, 1};
  } else {
    order[0] = {0, 0};
    order[1] = {0, 1};
    order[2] = {1, 1};
    order[3] = {1, 0};
  }

  const auto corner_point = [&](const std::array<int, 2>& bits) {
    return point_at(axis_a, bits[0] != 0 ? box.hi[axis_a] : box.lo[axis_a], axis_b,
                    bits[1] != 0 ? box.hi[axis_b] : box.lo[axis_b], family_axis, cap);
  };

  LoopPoints points;
  LoopSpans spans;
  points.reserve(8);
  spans.reserve(8);
  for (std::size_t index = 0; index < order.size(); ++index) {
    const std::array<int, 2> bits = order[index];
    const Point3d corner = corner_point(bits);
    if (!selected[static_cast<std::size_t>(bits[0])][static_cast<std::size_t>(bits[1])]) {
      points.push_back(corner);
      spans.push_back(std::nullopt);
      continue;
    }
    const Point3d previous = corner_point(order[(index + 3) % order.size()]);
    const Point3d next = corner_point(order[(index + 1) % order.size()]);
    const Vec3d incoming = unit_vector(corner - previous);
    const Vec3d outgoing = unit_vector(next - corner);
    const Point3d center = corner + (outgoing - incoming) * radius;
    const CircleCurve3 arc{center, outgoing * -1.0, incoming, radius};
    points.push_back(corner - incoming * radius);
    spans.push_back(CurvedArcSpan{arc, 0.0, std::numbers::pi / 2.0});
    points.push_back(corner + outgoing * radius);
    spans.push_back(std::nullopt);
  }

  if (!at_high) {
    LoopPoints reversed_points;
    LoopSpans reversed_spans(points.size());
    reversed_points.reserve(points.size());
    reversed_points.push_back(points.front());
    for (std::size_t index = points.size() - 1; index > 0; --index) {
      reversed_points.push_back(points[index]);
    }
    reversed_spans[0] = spans.back();
    for (std::size_t index = 1; index < spans.size(); ++index) {
      reversed_spans[index] = spans[spans.size() - 1 - index];
    }
    points = std::move(reversed_points);
    spans = std::move(reversed_spans);
  }

  std::vector<Point2d> uv;
  uv.reserve(points.size());
  for (const Point3d& point : points) {
    uv.push_back(Point2d{point[axis_a], point[axis_b]});
  }

  Vec3d normal = direction_f;
  double offset = box.hi[family_axis];
  if (!at_high) {
    normal = normal * -1.0;
    offset = -box.lo[family_axis];
  }
  return make_face(SurfaceGeometry{Plane3{normal, offset}}, std::move(points), std::move(spans),
                   std::move(uv));
}

// Side face parallel to the family axis. `p_axis` runs across the face and the
// family axis spans its full height; `p_lo`/`p_hi` are already inset by the
// selected corner radius.
inline CurvedAssembledFace side_face(const FilletBox& box, int plane_axis, int sign, int p_axis,
                                     double p_lo, double p_hi) {
  const int q_axis = 3 - plane_axis - p_axis;
  const double plane_value = sign > 0 ? box.hi[plane_axis] : box.lo[plane_axis];
  const Vec3d normal = axis_direction(plane_axis) * static_cast<double>(sign);
  const double offset = sign > 0 ? box.hi[plane_axis] : -box.lo[plane_axis];
  const double q_lo = box.lo[q_axis];
  const double q_hi = box.hi[q_axis];

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

// Quarter-cylinder patch at one selected cross-section corner. The angular
// interval is derived from the actual cylinder frame so the same code serves
// every family and corner.
inline CurvedAssembledFace cylinder_face(const FilletBox& box, int family_axis, int bits_a,
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
  const double height = box.hi[family_axis] - box.lo[family_axis];

  Point3d origin{};
  origin[axis_a] = bits_a != 0 ? box.hi[axis_a] - radius : box.lo[axis_a] + radius;
  origin[axis_b] = bits_b != 0 ? box.hi[axis_b] - radius : box.lo[axis_b] + radius;
  origin[family_axis] = box.lo[family_axis];

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

}  // namespace fillet_box_detail

// Recognises one axis-aligned, convex, closed box and returns its bounds. Every
// rejection is the stable box-only message so callers need no extra mapping.
inline Result<FilletBox> fillet_box_bounds(const Brep& brep, double tolerance) {
  const auto reject = [] {
    return Error{Errc::unsupported, "fillet_edges: only an axis-aligned box is supported"};
  };
  if (!std::isfinite(tolerance) || !(tolerance > 0.0)) return reject();
  if (!validate_topology(brep).ok()) return reject();
  const Result<BlendPlanarSolid> shape =
      blend_planar_convex_solid(brep, tolerance, "fillet_edges");
  if (!shape.has_value()) return reject();
  if (brep.vertex_count() != 8u || brep.edge_count() != 12u || brep.face_count() != 6u) {
    return reject();
  }
  if (shape.value().faces.size() != 6u) return reject();

  FilletBox box;
  for (int axis = 0; axis < 3; ++axis) {
    box.lo[axis] = std::numeric_limits<double>::infinity();
    box.hi[axis] = -std::numeric_limits<double>::infinity();
  }
  double scale = 1.0;
  for (const VertexId id : brep.vertex_ids()) {
    const Point3d& point = brep.vertex(id).point;
    for (int axis = 0; axis < 3; ++axis) {
      box.lo[axis] = std::min(box.lo[axis], point[axis]);
      box.hi[axis] = std::max(box.hi[axis], point[axis]);
      scale = std::max(scale, std::abs(point[axis]));
    }
  }
  for (int axis = 0; axis < 3; ++axis) {
    if (!(box.hi[axis] - box.lo[axis] > 0.0)) return reject();
  }
  const double slack = std::max(64.0 * tolerance, 1e-12 * scale);

  std::array<int, 8> corner_uses{};
  for (const VertexId id : brep.vertex_ids()) {
    const Point3d& point = brep.vertex(id).point;
    int pattern = 0;
    for (int axis = 0; axis < 3; ++axis) {
      const double distance_lo = std::abs(point[axis] - box.lo[axis]);
      const double distance_hi = std::abs(point[axis] - box.hi[axis]);
      if (distance_lo > slack && distance_hi > slack) return reject();
      if (distance_hi < distance_lo) pattern |= 1 << axis;
    }
    if (corner_uses[static_cast<std::size_t>(pattern)] != 0) return reject();
    ++corner_uses[static_cast<std::size_t>(pattern)];
  }
  for (const int uses : corner_uses) {
    if (uses != 1) return reject();
  }

  for (const EdgeId id : brep.edge_ids()) {
    const Edge& edge = brep.edge(id);
    const Point3d& start = brep.vertex(edge.v0).point;
    const Point3d& end = brep.vertex(edge.v1).point;
    int differing = 0;
    for (int axis = 0; axis < 3; ++axis) {
      if (std::abs(end[axis] - start[axis]) > slack) ++differing;
    }
    if (differing != 1) return reject();
  }

  for (const FaceId id : shape.value().faces) {
    const Result<std::vector<Point3d>> loop = blend_face_loop_points(brep, id, "fillet_edges");
    if (!loop.has_value() || loop.value().size() != 4u) return reject();
    const std::vector<Point3d>& points = loop.value();
    for (std::size_t index = 0; index < points.size(); ++index) {
      const Point3d& start = points[index];
      const Point3d& end = points[(index + 1) % points.size()];
      int differing = 0;
      for (int axis = 0; axis < 3; ++axis) {
        if (std::abs(end[axis] - start[axis]) > slack) ++differing;
      }
      if (differing != 1) return reject();
    }
  }
  return box;
}

// Builds the filleted solid for the box. Every selected corner of the one
// parallel family gets its own quarter-cylinder patch; the two caps and the
// four side faces are built once and walk all four cross-section corners.
inline Result<Brep> fillet_box_edges(const FilletBox& box, std::span<const int> corners,
                                     int family_axis, double radius, double tolerance) {
  if (!std::isfinite(radius)) {
    return Error{Errc::non_finite, "fillet_edges: radius must be finite"};
  }
  if (!(radius > 0.0)) {
    return Error{Errc::degenerate, "fillet_edges: radius must be positive"};
  }
  if (!std::isfinite(tolerance)) {
    return Error{Errc::non_finite, "fillet_edges: tolerance must be finite"};
  }
  if (!(tolerance > 0.0)) {
    return Error{Errc::tolerance_violation, "fillet_edges: tolerance must be positive"};
  }
  if (family_axis < 0 || family_axis > 2 || corners.empty()) {
    return Error{Errc::unsupported, "fillet_edges: only an axis-aligned box is supported"};
  }

  const int axis_a = fillet_box_detail::cross_axis_a(family_axis);
  const int axis_b = fillet_box_detail::cross_axis_b(family_axis);
  const double extent_a = box.hi[axis_a] - box.lo[axis_a];
  const double extent_b = box.hi[axis_b] - box.lo[axis_b];
  if (!(extent_a > 0.0) || !(extent_b > 0.0)) {
    return Error{Errc::unsupported, "fillet_edges: only an axis-aligned box is supported"};
  }
  if (!(radius < 0.5 * std::min(extent_a, extent_b))) {
    return Error{Errc::unsupported, "fillet_edges: radius exceeds the box cross-section"};
  }

  std::array<std::array<bool, 2>, 2> selected{};
  std::size_t selected_count = 0;
  for (const int corner : corners) {
    if (corner < 0 || corner > 3) {
      return Error{Errc::unsupported, "fillet_edges: edge is not a box corner edge"};
    }
    bool& slot = selected[static_cast<std::size_t>(corner & 1)]
                         [static_cast<std::size_t>((corner >> 1) & 1)];
    if (!slot) {
      slot = true;
      ++selected_count;
    }
  }

  const double a_side_lo = box.lo[axis_b] +
                           (selected[1][0] ? radius : 0.0);
  const double a_side_hi = box.hi[axis_b] -
                           (selected[1][1] ? radius : 0.0);
  const double a_side_lo_rev = box.lo[axis_b] +
                               (selected[0][0] ? radius : 0.0);
  const double a_side_hi_rev = box.hi[axis_b] -
                               (selected[0][1] ? radius : 0.0);
  const double b_side_lo = box.lo[axis_a] +
                           (selected[0][1] ? radius : 0.0);
  const double b_side_hi = box.hi[axis_a] -
                           (selected[1][1] ? radius : 0.0);
  const double b_side_lo_rev = box.lo[axis_a] +
                               (selected[0][0] ? radius : 0.0);
  const double b_side_hi_rev = box.hi[axis_a] -
                               (selected[1][0] ? radius : 0.0);

  std::vector<CurvedAssembledFace> faces;
  faces.reserve(6u + selected_count);
  faces.push_back(fillet_box_detail::cap_face(box, family_axis, selected, true, radius));
  faces.push_back(fillet_box_detail::cap_face(box, family_axis, selected, false, radius));
  faces.push_back(
      fillet_box_detail::side_face(box, axis_a, 1, axis_b, a_side_lo, a_side_hi));
  faces.push_back(
      fillet_box_detail::side_face(box, axis_a, -1, axis_b, a_side_lo_rev, a_side_hi_rev));
  faces.push_back(
      fillet_box_detail::side_face(box, axis_b, 1, axis_a, b_side_lo, b_side_hi));
  faces.push_back(
      fillet_box_detail::side_face(box, axis_b, -1, axis_a, b_side_lo_rev, b_side_hi_rev));
  for (int bits_a = 0; bits_a < 2; ++bits_a) {
    for (int bits_b = 0; bits_b < 2; ++bits_b) {
      if (!selected[static_cast<std::size_t>(bits_a)][static_cast<std::size_t>(bits_b)]) {
        continue;
      }
      faces.push_back(
          fillet_box_detail::cylinder_face(box, family_axis, bits_a, bits_b, radius));
    }
  }

  Result<Brep> assembled = assemble_curved_solid(faces, tolerance);
  if (!assembled.has_value()) return assembled.error();
  if (!validate_topology(assembled.value()).ok()) {
    return Error{Errc::invalid_topology, "fillet_edges: fillet produced an invalid solid"};
  }
  return assembled;
}

}  // namespace detail
}  // namespace v1
}  // namespace geom
