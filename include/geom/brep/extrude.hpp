#pragma once

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <span>
#include <vector>

#include <geom/brep/assemble.hpp>
#include <geom/brep/plane_frame.hpp>
#include <geom/brep/topology.hpp>
#include <geom/cgeom/polygon2.hpp>
#include <geom/core/error.hpp>
#include <geom/core/result.hpp>
#include <geom/surfaces/plane.hpp>

namespace geom {
inline namespace v1 {

namespace detail {

// Outward unit normal of one extrusion side wall. `a` and `b` are consecutive
// points of the (counter-clockwise) bottom ring and `axis` is the extrusion
// direction; cross(b - a, axis) points away from the enclosed material.
// The largest component is scaled out first so the length stays in
// [1, sqrt(3)] for extreme coordinates. A zero or non-finite cross product
// returns a zero vector, which assemble_solid rejects as a degenerate face.
inline Vec3d extrude_side_normal(const Point3d& a, const Point3d& b,
                                 const Vec3d& axis) noexcept {
  const Vec3d edge = b - a;
  const Vec3d raw = cross(edge, axis);
  const double scale =
      std::max({std::fabs(raw[0]), std::fabs(raw[1]), std::fabs(raw[2])});
  if (!(scale > 0.0) || !std::isfinite(scale)) return Vec3d{};
  const double x = raw[0] / scale;
  const double y = raw[1] / scale;
  const double z = raw[2] / scale;
  const double length = std::sqrt(x * x + y * y + z * z);
  return Vec3d{x / length, y / length, z / length};
}

}  // namespace detail

// Extrudes the outer ring of a planar polygon into a closed B-rep solid.
//
// The ring is mapped through the plane's frame, normalized to counter-clockwise
// (a clockwise input is accepted and reversed), and swept along the plane normal
// by `height`. The result has one bottom cap, one top cap and one planar quad
// per ring edge; every face normal points out of the solid. Task 1 accepts only
// a profile without holes.
//
// Errors, in check order:
//   - non-finite or non-positive tolerance -> Errc::tolerance_violation;
//   - make_plane_frame failure -> that error unchanged;
//   - non-finite height -> Errc::non_finite; height <= 0 -> Errc::degenerate;
//   - validate_polygon2 failure -> its Errc::non_finite or Errc::degenerate;
//   - more than one ring -> Errc::unsupported (holes are a later task);
//   - undefined ring orientation -> Errc::degenerate;
//   - consecutive equal ring points, including last == first -> degenerate;
//   - assemble_solid failure -> Errc::invalid_topology (the assembler's
//     internal message is not exposed).
inline Result<Brep> extrude_profile(const Plane3& plane,
                                    const Polygon2d& profile, double height,
                                    double tolerance) {
  if (!std::isfinite(tolerance) || !(tolerance > 0.0)) {
    return Error{Errc::tolerance_violation,
                 "extrude_profile: tolerance must be finite and positive"};
  }

  const Result<PlaneFrame> frame_result = make_plane_frame(plane);
  if (!frame_result) return frame_result.error();
  const PlaneFrame& frame = frame_result.value();

  if (!std::isfinite(height)) {
    return Error{Errc::non_finite, "extrude_profile: height is not finite"};
  }
  if (!(height > 0.0)) {
    return Error{Errc::degenerate, "extrude_profile: height must be positive"};
  }

  const Errc profile_status = validate_polygon2(profile);
  if (profile_status == Errc::non_finite) {
    return Error{Errc::non_finite,
                 "extrude_profile: profile coordinates are not finite"};
  }
  if (profile_status != Errc::none) {
    return Error{Errc::degenerate, "extrude_profile: profile is degenerate"};
  }

  if (profile.ring_count() != 1u) {
    return Error{Errc::unsupported,
                 "extrude_profile: hole rings are not supported yet"};
  }

  const std::span<const Point2d> source = profile.ring(0u);
  std::vector<Point2d> ring;
  ring.reserve(source.size());
  const int orientation = ring_orientation2(source);
  if (orientation == 0) {
    return Error{Errc::degenerate,
                 "extrude_profile: profile ring is degenerate"};
  }
  if (orientation > 0) {
    ring.assign(source.begin(), source.end());
  } else {
    ring.assign(source.rbegin(), source.rend());
  }

  // Repeated points (including a duplicated closing vertex) are not rejected
  // by ring_orientation2, which skips them when it picks the extreme vertex.
  const auto same_point = [](const Point2d& first,
                             const Point2d& second) noexcept {
    return first[0] == second[0] && first[1] == second[1];
  };
  for (std::size_t i = 0; i < ring.size(); ++i) {
    if (same_point(ring[i], ring[(i + 1u) % ring.size()])) {
      return Error{Errc::degenerate,
                   "extrude_profile: profile ring has a repeated vertex"};
    }
  }

  std::vector<Point3d> bottom;
  std::vector<Point3d> top;
  bottom.reserve(ring.size());
  top.reserve(ring.size());
  for (const Point2d& point : ring) {
    const Point3d base = plane_frame_to_3d(frame, point);
    bottom.push_back(base);
    top.push_back(base + frame.normal * height);
  }

  const Vec3d& normal = frame.normal;
  std::vector<AssembledFace> faces;
  faces.reserve(2u + ring.size());

  // The top cap keeps the normalized counter-clockwise ring; the bottom cap is
  // the same ring reversed so it winds counter-clockwise about -normal.
  AssembledFace bottom_cap;
  bottom_cap.outer.resize(bottom.size());
  std::reverse_copy(bottom.begin(), bottom.end(), bottom_cap.outer.begin());
  bottom_cap.normal = Vec3d{-normal[0], -normal[1], -normal[2]};
  faces.push_back(std::move(bottom_cap));

  AssembledFace top_cap;
  top_cap.outer = top;
  top_cap.normal = normal;
  faces.push_back(std::move(top_cap));

  for (std::size_t i = 0; i < ring.size(); ++i) {
    const std::size_t next = (i + 1u) % ring.size();
    AssembledFace side;
    side.outer = {bottom[i], bottom[next], top[next], top[i]};
    side.normal = detail::extrude_side_normal(bottom[i], bottom[next], normal);
    faces.push_back(std::move(side));
  }

  const Result<Brep> assembled = assemble_solid(faces, tolerance);
  if (!assembled) {
    return Error{Errc::invalid_topology,
                 "extrude_profile: solid could not be assembled"};
  }
  return assembled;
}

}  // namespace v1
}  // namespace geom
