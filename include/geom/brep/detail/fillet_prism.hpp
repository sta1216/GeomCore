#pragma once

#include <algorithm>
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
#include <geom/brep/face_uv.hpp>
#include <geom/brep/plane_frame.hpp>
#include <geom/brep/topology.hpp>
#include <geom/brep/trim_face.hpp>
#include <geom/brep/validate.hpp>
#include <geom/cgeom/polygon2.hpp>
#include <geom/core/error.hpp>
#include <geom/core/result.hpp>
#include <geom/curves/analytic_curve.hpp>
#include <geom/numeric/predicates2.hpp>
#include <geom/primitives/point.hpp>
#include <geom/primitives/vec.hpp>
#include <geom/surfaces/analytic.hpp>
#include <geom/surfaces/plane.hpp>

namespace geom {
inline namespace v1 {
namespace detail {

// A convex prism recognised in a Brep, in the frame of its extrusion axis.
// `ring` is the bottom cap CCW about `axis`; `vertical[i]` is the edge joining
// ring[i] to its top counterpart. Its stored direction is not guaranteed.
struct FilletPrism {
  Point3d bottom{};
  Vec3d axis{0.0, 0.0, 1.0};
  double height = 0.0;
  std::vector<Point3d> ring;
  std::vector<EdgeId> vertical;
};

namespace fillet_prism_detail {

using LoopPoints = std::vector<Point3d>;
using LoopSpans = std::vector<std::optional<CurvedArcSpan>>;

inline double point_distance(const Point3d& a, const Point3d& b) {
  return std::hypot(b[0] - a[0], b[1] - a[1], b[2] - a[2]);
}

inline Point3d point_offset(const Point3d& base, const Vec3d& direction, double amount) {
  return Point3d{base[0] + direction[0] * amount, base[1] + direction[1] * amount,
                 base[2] + direction[2] * amount};
}

inline Vec3d unit_vector(const Vec3d& value) {
  const double length = norm2(value);
  if (!(length > 0.0) || !std::isfinite(length)) return Vec3d{};
  return value * (1.0 / length);
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

// The cap ring as 3D points in face-loop order. Only the point positions are
// used; the cap is known to be planar, so the uv chart introduces no bias.
inline bool read_cap_ring(const Brep& brep, FaceId face_id, double tolerance,
                          std::vector<Point3d>& points) {
  const Result<std::vector<FaceUvLoop>> loops = face_loops_uv(brep, face_id, tolerance);
  if (!loops.has_value() || loops.value().size() != 1u) return false;
  const Result<FaceUvProjection> projection =
      face_uv_projection(brep, brep.face(face_id), tolerance);
  if (!projection.has_value()) return false;
  const std::vector<Point2d>& uv = loops.value().front().points;
  if (uv.size() < 3u) return false;
  points.clear();
  points.reserve(uv.size());
  for (const Point2d& value : uv) points.push_back(face_uv_eval(projection.value(), value));
  return true;
}

}  // namespace fillet_prism_detail

// Recognises a convex prism with planar caps perpendicular to a common axis and
// strictly convex polygonal cross-sections. Returns the bottom ring, the axis
// and the vertical edges of every ring vertex; anything else is unsupported.
inline Result<FilletPrism> fillet_prism_bounds(const Brep& brep, double tolerance) {
  const auto reject = []() -> Error {
    return Error{Errc::unsupported, "fillet_edges: only a box or a convex prism is supported"};
  };

  if (!std::isfinite(tolerance) || !(tolerance > 0.0)) return reject();
  if (!validate_topology(brep).ok()) return reject();
  const Result<BlendPlanarSolid> blended =
      blend_planar_convex_solid(brep, tolerance, "fillet_edges");
  if (!blended.has_value()) return reject();
  const BlendPlanarSolid& shape = blended.value();
  if (shape.faces.size() < 5u) return reject();

  // One boundary loop per face; the cap rings are read through that loop.
  for (const FaceId face_id : shape.faces) {
    const Face& face = brep.face(face_id);
    if (face.first_loop.is_null()) return reject();
    std::size_t loop_count = 0;
    for (LoopId loop_id = face.first_loop; !loop_id.is_null();
         loop_id = brep.loop(loop_id).next) {
      if (++loop_count > 1u) return reject();
    }
    if (loop_count != 1u) return reject();
  }

  // Every face normal is a candidate axis: exactly two faces must be parallel to
  // it (the caps) and every other face perpendicular (the sides).
  for (std::size_t first = 0; first < shape.faces.size(); ++first) {
    const Vec3d axis = fillet_prism_detail::unit_vector(shape.outward[first]);
    if (!(norm2(axis) > 0.0)) continue;

    std::size_t caps[2] = {shape.faces.size(), shape.faces.size()};
    std::size_t cap_count = 0;
    bool perpendicular = true;
    for (std::size_t other = 0; other < shape.faces.size(); ++other) {
      const double alignment = std::abs(dot(shape.outward[other], axis));
      if (alignment > 1.0 - 1e-9) {
        if (cap_count < 2u) caps[cap_count] = other;
        ++cap_count;
      } else if (alignment > 1e-9) {
        perpendicular = false;
        break;
      }
    }
    if (!perpendicular || cap_count != 2u) continue;
    if (!(dot(shape.outward[caps[0]], shape.outward[caps[1]]) < -1.0 + 1e-9)) continue;

    const std::size_t top = dot(shape.outward[caps[0]], axis) > 0.0 ? caps[0] : caps[1];
    const std::size_t bottom = top == caps[0] ? caps[1] : caps[0];
    if (!(dot(shape.outward[top], axis) > 1.0 - 1e-9)) continue;
    if (!(dot(shape.outward[bottom], axis) < -1.0 + 1e-9)) continue;
    const double height = shape.offset[top] + shape.offset[bottom];
    if (!std::isfinite(height) || !(height > 0.0)) continue;

    std::vector<Point3d> bottom_points;
    std::vector<Point3d> top_points;
    if (!fillet_prism_detail::read_cap_ring(brep, shape.faces[bottom], tolerance, bottom_points)) {
      continue;
    }
    if (!fillet_prism_detail::read_cap_ring(brep, shape.faces[top], tolerance, top_points)) {
      continue;
    }
    const std::size_t n = bottom_points.size();
    if (n < 3u || n > 64u || top_points.size() != n) continue;

    // Both caps project onto the same planar chart perpendicular to the axis.
    const Result<PlaneFrame> frame_result = make_plane_frame(Plane3{axis, 0.0});
    if (!frame_result.has_value()) continue;
    const PlaneFrame& frame = frame_result.value();
    std::vector<Point2d> bottom_planar;
    std::vector<Point2d> top_planar;
    bottom_planar.reserve(n);
    top_planar.reserve(n);
    for (const Point3d& point : bottom_points) {
      bottom_planar.push_back(plane_frame_to_2d(frame, point));
    }
    for (const Point3d& point : top_points) top_planar.push_back(plane_frame_to_2d(frame, point));

    const int bottom_orientation = ring_orientation2(bottom_planar);
    const int top_orientation = ring_orientation2(top_planar);
    if (bottom_orientation == 0 || top_orientation == 0) continue;
    if (bottom_orientation < 0) {
      std::reverse(bottom_points.begin(), bottom_points.end());
      std::reverse(bottom_planar.begin(), bottom_planar.end());
    }
    if (top_orientation < 0) {
      std::reverse(top_points.begin(), top_points.end());
      std::reverse(top_planar.begin(), top_planar.end());
    }

    // Strictly convex: every triple turns left.
    bool convex = true;
    for (std::size_t index = 0; index < n && convex; ++index) {
      const Point2d& previous = bottom_planar[(index + n - 1u) % n];
      const Point2d& current = bottom_planar[index];
      const Point2d& next = bottom_planar[(index + 1u) % n];
      convex = sign_orient2d(previous[0], previous[1], current[0], current[1], next[0],
                             next[1]) > 0;
    }
    for (std::size_t index = 0; index < n && convex; ++index) {
      const Point2d& previous = top_planar[(index + n - 1u) % n];
      const Point2d& current = top_planar[index];
      const Point2d& next = top_planar[(index + 1u) % n];
      convex = sign_orient2d(previous[0], previous[1], current[0], current[1], next[0],
                             next[1]) > 0;
    }
    if (!convex) continue;

    // Cyclic alignment by minimising the worst planar distance.
    std::size_t shift = 0;
    double worst = std::numeric_limits<double>::infinity();
    for (std::size_t candidate = 0; candidate < n; ++candidate) {
      double candidate_worst = 0.0;
      for (std::size_t index = 0; index < n; ++index) {
        const Point2d& a = top_planar[(index + candidate) % n];
        const Point2d& b = bottom_planar[index];
        candidate_worst = std::max(candidate_worst, std::hypot(a[0] - b[0], a[1] - b[1]));
      }
      if (candidate_worst < worst) {
        worst = candidate_worst;
        shift = candidate;
      }
    }
    const double alignment_tolerance = 64.0 * tolerance;
    if (!(worst <= alignment_tolerance)) continue;

    std::vector<Point3d> top_aligned(n);
    bool aligned = true;
    for (std::size_t index = 0; index < n; ++index) {
      top_aligned[index] = top_points[(index + shift) % n];
      const Point3d predicted =
          fillet_prism_detail::point_offset(bottom_points[index], axis, height);
      if (fillet_prism_detail::point_distance(top_aligned[index], predicted) >
          alignment_tolerance) {
        aligned = false;
        break;
      }
    }
    if (!aligned) continue;

    // Exactly one edge must join ring[index] to the aligned top vertex.
    std::vector<EdgeId> vertical(n);
    std::vector<int> matches(n, 0);
    for (const EdgeId edge_id : brep.edge_ids()) {
      const Edge& edge = brep.edge(edge_id);
      const Point3d& a = brep.vertex(edge.v0).point;
      const Point3d& b = brep.vertex(edge.v1).point;
      for (std::size_t index = 0; index < n; ++index) {
        const bool forward = fillet_prism_detail::point_distance(a, bottom_points[index]) <=
                                 tolerance &&
                             fillet_prism_detail::point_distance(b, top_aligned[index]) <=
                                 tolerance;
        const bool backward = fillet_prism_detail::point_distance(b, bottom_points[index]) <=
                                  tolerance &&
                              fillet_prism_detail::point_distance(a, top_aligned[index]) <=
                                  tolerance;
        if (forward || backward) {
          vertical[index] = edge_id;
          ++matches[index];
        }
      }
    }
    bool edges_ok = true;
    for (std::size_t index = 0; index < n; ++index) {
      if (matches[index] != 1) {
        edges_ok = false;
        break;
      }
    }
    for (std::size_t left = 0; left < n && edges_ok; ++left) {
      for (std::size_t right = left + 1u; right < n; ++right) {
        if (vertical[left] == vertical[right]) {
          edges_ok = false;
          break;
        }
      }
    }
    if (!edges_ok) continue;

    return FilletPrism{bottom_points.front(), axis, height, std::move(bottom_points),
                       std::move(vertical)};
  }
  return reject();
}

// Rounds the selected vertical edges of a recognised convex prism. Each selected
// corner is replaced by a tangent arc on both caps and a matching cylinder
// patch; the caps and side faces are trimmed to the tangent points.
inline Result<Brep> fillet_prism_edges(const FilletPrism& prism,
                                       std::span<const std::uint32_t> vertices, double radius,
                                       double tolerance) {
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

  const std::size_t n = prism.ring.size();
  if (n < 3u || prism.vertical.size() != n || vertices.empty()) {
    return Error{Errc::unsupported, "fillet_edges: only vertical prism edges are supported"};
  }
  const Vec3d axis = fillet_prism_detail::unit_vector(prism.axis);
  if (!(norm2(axis) > 0.0) || !std::isfinite(prism.height) || !(prism.height > 0.0)) {
    return Error{Errc::unsupported, "fillet_edges: only a box or a convex prism is supported"};
  }

  std::vector<bool> selected(n, false);
  for (const std::uint32_t index : vertices) {
    if (static_cast<std::size_t>(index) >= n) {
      return Error{Errc::unsupported, "fillet_edges: only vertical prism edges are supported"};
    }
    selected[index] = true;
  }

  const Result<PlaneFrame> frame_result = make_plane_frame(Plane3{axis, 0.0});
  if (!frame_result.has_value()) {
    return Error{Errc::unsupported, "fillet_edges: only a box or a convex prism is supported"};
  }
  const PlaneFrame& frame = frame_result.value();

  struct Corner {
    bool selected = false;
    double trim = 0.0;
    double sweep = 0.0;
    Point3d tangent_in{};
    Point3d tangent_out{};
    CircleCurve3 circle{};
  };
  std::vector<Corner> corners(n);

  double centroid_x = 0.0;
  double centroid_y = 0.0;
  double centroid_z = 0.0;
  for (const Point3d& point : prism.ring) {
    centroid_x += point[0];
    centroid_y += point[1];
    centroid_z += point[2];
  }
  const double inverse_count = 1.0 / static_cast<double>(n);
  const Point3d centroid{centroid_x * inverse_count, centroid_y * inverse_count,
                         centroid_z * inverse_count};
  const double pi = std::numbers::pi_v<double>;

  for (std::size_t index = 0; index < n; ++index) {
    if (!selected[index]) continue;
    const Point3d& vertex = prism.ring[index];
    const Point3d& previous = prism.ring[(index + n - 1u) % n];
    const Point3d& next = prism.ring[(index + 1u) % n];

    const Vec3d incoming = fillet_prism_detail::unit_vector(previous - vertex);
    const Vec3d outgoing = fillet_prism_detail::unit_vector(next - vertex);
    if (!(norm2(incoming) > 0.0) || !(norm2(outgoing) > 0.0)) {
      return Error{Errc::unsupported, "fillet_edges: only a box or a convex prism is supported"};
    }
    const Vec3d sine = cross(incoming, outgoing);
    const double theta = std::atan2(std::abs(sine[2]), dot(incoming, outgoing));
    if (!(theta > 0.0) || !(theta < pi)) {
      return Error{Errc::unsupported, "fillet_edges: only a box or a convex prism is supported"};
    }
    const double tangent_length = radius / std::tan(theta / 2.0);
    if (!(tangent_length > 0.0) ||
        !(tangent_length < fillet_prism_detail::point_distance(vertex, previous)) ||
        !(tangent_length < fillet_prism_detail::point_distance(vertex, next))) {
      return Error{Errc::unsupported, "fillet_edges: radius exceeds the prism cross-section"};
    }

    const Point3d tangent_in = fillet_prism_detail::point_offset(vertex, incoming, tangent_length);
    const Point3d tangent_out =
        fillet_prism_detail::point_offset(vertex, outgoing, tangent_length);

    // Inward normals of the two edges meeting at the vertex.
    Vec3d inward_in = cross(axis, incoming * -1.0);
    Vec3d inward_out = cross(axis, outgoing);
    const Vec3d toward = centroid - vertex;
    if (dot(inward_in, toward) < 0.0) inward_in = inward_in * -1.0;
    if (dot(inward_out, toward) < 0.0) inward_out = inward_out * -1.0;
    inward_in = fillet_prism_detail::unit_vector(inward_in);
    inward_out = fillet_prism_detail::unit_vector(inward_out);

    const double a11 = dot(inward_in, frame.u);
    const double a12 = dot(inward_in, frame.v);
    const double a21 = dot(inward_out, frame.u);
    const double a22 = dot(inward_out, frame.v);
    const double determinant = a11 * a22 - a12 * a21;
    if (!(determinant > 0.0)) {
      return Error{Errc::unsupported, "fillet_edges: radius exceeds the prism cross-section"};
    }
    const double offset_u = radius * (a22 - a12) / determinant;
    const double offset_v = radius * (a11 - a21) / determinant;
    const Point3d center = vertex + frame.u * offset_u + frame.v * offset_v;

    const Vec3d radial = fillet_prism_detail::unit_vector(tangent_in - center);
    if (!(norm2(radial) > 0.0)) {
      return Error{Errc::unsupported, "fillet_edges: radius exceeds the prism cross-section"};
    }
    const Vec3d tangent_axis = cross(axis, radial);

    Corner& corner = corners[index];
    corner.selected = true;
    corner.trim = tangent_length;
    corner.sweep = pi - theta;
    corner.tangent_in = tangent_in;
    corner.tangent_out = tangent_out;
    corner.circle = CircleCurve3{center, radial, tangent_axis, radius};
  }

  // Every side face keeps a positive strip between its trimmed ends.
  for (std::size_t index = 0; index < n; ++index) {
    const std::size_t next = (index + 1u) % n;
    const double edge_length = fillet_prism_detail::point_distance(prism.ring[index],
                                                                   prism.ring[next]);
    const double left_trim = selected[index] ? corners[index].trim : 0.0;
    const double right_trim = selected[next] ? corners[next].trim : 0.0;
    if (!(edge_length - left_trim - right_trim > tolerance)) {
      return Error{Errc::unsupported,
                   "fillet_edges: radius leaves a prism face narrower than the tolerance"};
    }
  }

  const double height = prism.height;
  // Plane offsets along `axis`: the recognised ring may sit at any axis
  // coordinate, not necessarily at zero.
  const double base = axis[0] * prism.ring.front()[0] + axis[1] * prism.ring.front()[1] +
                      axis[2] * prism.ring.front()[2];
  const double top_offset = base + height;
  const auto build_cap = [&](bool at_top) {
    fillet_prism_detail::LoopPoints points;
    fillet_prism_detail::LoopSpans spans;
    points.reserve(n + 1u);
    spans.reserve(n + 1u);
    for (std::size_t index = 0; index < n; ++index) {
      const Corner& corner = corners[index];
      if (corner.selected) {
        const Point3d tangent_in =
            at_top ? fillet_prism_detail::point_offset(corner.tangent_in, axis, height)
                   : corner.tangent_in;
        const Point3d tangent_out =
            at_top ? fillet_prism_detail::point_offset(corner.tangent_out, axis, height)
                   : corner.tangent_out;
        CircleCurve3 arc_circle = corner.circle;
        if (at_top) {
          arc_circle.center =
              fillet_prism_detail::point_offset(corner.circle.center, axis, height);
        }
        points.push_back(tangent_in);
        spans.push_back(CurvedArcSpan{arc_circle, 0.0, corner.sweep});
        points.push_back(tangent_out);
        spans.push_back(std::nullopt);
      } else {
        const Point3d base = prism.ring[index];
        points.push_back(at_top ? fillet_prism_detail::point_offset(base, axis, height) : base);
        spans.push_back(std::nullopt);
      }
    }
    return std::make_pair(std::move(points), std::move(spans));
  };

  std::pair<fillet_prism_detail::LoopPoints, fillet_prism_detail::LoopSpans> top_cap =
      build_cap(true);
  const std::pair<fillet_prism_detail::LoopPoints, fillet_prism_detail::LoopSpans> bottom_cap =
      build_cap(false);

  // The bottom cap winds CCW about -axis, which is the reverse of the ring order.
  fillet_prism_detail::LoopPoints bottom_points;
  bottom_points.reserve(bottom_cap.first.size());
  bottom_points.push_back(bottom_cap.first.front());
  for (std::size_t index = bottom_cap.first.size(); index-- > 1u;) {
    bottom_points.push_back(bottom_cap.first[index]);
  }
  fillet_prism_detail::LoopSpans bottom_spans(bottom_cap.second.size());
  bottom_spans[0] = bottom_cap.second.back();
  for (std::size_t index = 1; index < bottom_cap.second.size(); ++index) {
    bottom_spans[index] = bottom_cap.second[bottom_cap.second.size() - 1u - index];
  }

  const auto planar = [&frame](const fillet_prism_detail::LoopPoints& points) {
    std::vector<Point2d> uv;
    uv.reserve(points.size());
    for (const Point3d& point : points) uv.push_back(plane_frame_to_2d(frame, point));
    return uv;
  };

  std::vector<CurvedAssembledFace> faces;
  faces.reserve(2u + n + vertices.size());
  std::vector<Point2d> top_uv = planar(top_cap.first);
  faces.push_back(fillet_prism_detail::make_face(SurfaceGeometry{Plane3{axis, top_offset}},
                                                 std::move(top_cap.first),
                                                 std::move(top_cap.second), std::move(top_uv)));
  std::vector<Point2d> bottom_uv = planar(bottom_points);
  faces.push_back(fillet_prism_detail::make_face(SurfaceGeometry{Plane3{axis * -1.0, -base}},
                                                 std::move(bottom_points),
                                                 std::move(bottom_spans), std::move(bottom_uv)));

  for (std::size_t index = 0; index < n; ++index) {
    const std::size_t next = (index + 1u) % n;
    const Point3d left = corners[index].selected ? corners[index].tangent_out
                                                  : prism.ring[index];
    const Point3d right = corners[next].selected ? corners[next].tangent_in
                                                  : prism.ring[next];
    const Point3d left_top = fillet_prism_detail::point_offset(left, axis, height);
    const Point3d right_top = fillet_prism_detail::point_offset(right, axis, height);
    const Vec3d edge = right - left;
    const Vec3d normal = fillet_prism_detail::unit_vector(cross(edge, axis));
    if (!(norm2(normal) > 0.0)) {
      return Error{Errc::invalid_topology, "fillet_edges: fillet could not be assembled"};
    }
    const double offset = normal[0] * left[0] + normal[1] * left[1] + normal[2] * left[2];
    const Plane3 plane{normal, offset};
    const Result<PlaneFrame> side_frame = make_plane_frame(plane);
    if (!side_frame.has_value()) {
      return Error{Errc::invalid_topology, "fillet_edges: fillet could not be assembled"};
    }
    fillet_prism_detail::LoopPoints points{left, right, right_top, left_top};
    std::vector<Point2d> uv;
    uv.reserve(points.size());
    for (const Point3d& point : points) {
      uv.push_back(plane_frame_to_2d(side_frame.value(), point));
    }
    faces.push_back(fillet_prism_detail::make_face(SurfaceGeometry{plane}, std::move(points),
                                                   fillet_prism_detail::LoopSpans{},
                                                   std::move(uv)));
  }

  for (std::size_t index = 0; index < n; ++index) {
    const Corner& corner = corners[index];
    if (!corner.selected) continue;
    const Point3d in_bottom = corner.tangent_in;
    const Point3d out_bottom = corner.tangent_out;
    const Point3d in_top = fillet_prism_detail::point_offset(in_bottom, axis, height);
    const Point3d out_top = fillet_prism_detail::point_offset(out_bottom, axis, height);
    fillet_prism_detail::LoopPoints points{in_bottom, out_bottom, out_top, in_top};
    fillet_prism_detail::LoopSpans spans;
    spans.reserve(points.size());
    spans.push_back(CurvedArcSpan{corner.circle, 0.0, corner.sweep});
    spans.push_back(std::nullopt);
    CircleCurve3 top_circle = corner.circle;
    top_circle.center = fillet_prism_detail::point_offset(corner.circle.center, axis, height);
    spans.push_back(CurvedArcSpan{top_circle, 0.0, corner.sweep});
    spans.push_back(std::nullopt);
    std::vector<Point2d> uv{Point2d{0.0, 0.0}, Point2d{corner.sweep, 0.0},
                            Point2d{corner.sweep, height}, Point2d{0.0, height}};
    faces.push_back(fillet_prism_detail::make_face(
        SurfaceGeometry{Cylinder3{corner.circle.center, axis, corner.circle.x_axis, radius}},
        std::move(points), std::move(spans), std::move(uv)));
  }

  Result<Brep> assembled = assemble_curved_solid(faces, tolerance);
  if (!assembled.has_value()) {
    return Error{Errc::invalid_topology, "fillet_edges: fillet could not be assembled"};
  }
  if (!validate_topology(assembled.value()).ok()) {
    return Error{Errc::invalid_topology, "fillet_edges: fillet produced an invalid solid"};
  }
  return assembled;
}

}  // namespace detail
}  // namespace v1
}  // namespace geom

