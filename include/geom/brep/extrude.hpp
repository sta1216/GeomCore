#pragma once

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <span>
#include <vector>

#include <geom/brep/assemble.hpp>
#include <geom/brep/plane_frame.hpp>
#include <geom/brep/topology.hpp>
#include <geom/cgeom/polygon2.hpp>
#include <geom/cgeom/triangulate2.hpp>
#include <geom/core/error.hpp>
#include <geom/core/result.hpp>
#include <geom/surfaces/plane.hpp>

namespace geom {
inline namespace v1 {

namespace detail {

// Outward unit normal of one extrusion side wall. `a` and `b` are consecutive
// points of a normalized bottom ring and `axis` is the extrusion direction; the
// outer ring is counter-clockwise and holes are clockwise, so cross(b - a, axis)
// points away from the material in both cases.
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

// Triangulation of the bottom cap from the top cap's triangulation. The bottom
// cap is built from every ring reversed, so the vertex at position i of a ring
// moves to position n - 1 - i. Reading the top triangulation through that
// mapping keeps every triangle's cyclic order unchanged, i.e. still
// counter-clockwise about the top normal and therefore clockwise about the
// bottom normal, so the last two indices of each triangle are swapped.
inline Triangulation2 extrude_reverse_triangulation(
    const Triangulation2& source,
    std::span<const std::uint32_t> ring_sizes) noexcept {
  Triangulation2 reversed;
  reversed.ok = source.ok;
  std::vector<std::uint32_t> starts(ring_sizes.size() + 1u, 0u);
  for (std::size_t r = 0; r < ring_sizes.size(); ++r) {
    starts[r + 1u] = starts[r] + ring_sizes[r];
  }
  const auto map_vertex = [&](std::uint32_t index) noexcept -> std::uint32_t {
    for (std::size_t r = 0; r < ring_sizes.size(); ++r) {
      if (index < starts[r + 1u]) {
        return starts[r] + (starts[r + 1u] - 1u - index);
      }
    }
    return index;
  };
  reversed.triangles.reserve(source.triangles.size());
  for (const std::array<std::uint32_t, 3>& triangle : source.triangles) {
    reversed.triangles.push_back({map_vertex(triangle[0]),
                                  map_vertex(triangle[2]),
                                  map_vertex(triangle[1])});
  }
  return reversed;
}

}  // namespace detail

// Extrudes a planar polygon, with optional holes, into a closed B-rep solid.
//
// Every ring is mapped through the plane's frame and normalized: the outer ring
// (ring 0) ends counter-clockwise and every hole ends clockwise, so the material
// lies to the left of each ring; an input ring wound the other way is reversed.
// The normalized rings are swept along the plane normal by `height`. The result
// has one bottom cap, one top cap and one planar quad per edge of every ring;
// every face normal points out of the solid.
//
// Errors, in check order:
//   - non-finite or non-positive tolerance -> Errc::tolerance_violation;
//   - make_plane_frame failure -> that error unchanged;
//   - non-finite height -> Errc::non_finite; height <= 0 -> Errc::degenerate;
//   - validate_polygon2 failure -> its Errc::non_finite or Errc::degenerate;
//   - undefined orientation of any ring -> Errc::degenerate;
//   - consecutive equal points in any normalized ring, including last == first
//     -> Errc::degenerate;
//   - a profile that triangulate_polygon2 cannot triangulate, i.e. is not a
//     simple polygon with holes, -> Errc::degenerate;
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

  std::vector<std::vector<Point2d>> rings;
  rings.reserve(profile.ring_count());
  for (std::uint32_t r = 0; r < profile.ring_count(); ++r) {
    const std::span<const Point2d> source = profile.ring(r);
    const int orientation = ring_orientation2(source);
    if (orientation == 0) {
      return Error{Errc::degenerate,
                   "extrude_profile: profile ring is degenerate"};
    }
    // Outer rings must end counter-clockwise and holes clockwise, so the
    // material is on the left of each ring; reverse an input wound the other
    // way. A clockwise outer ring is therefore accepted.
    const bool wants_ccw = r == 0u;
    std::vector<Point2d> ring;
    ring.reserve(source.size());
    if ((orientation > 0) == wants_ccw) {
      ring.assign(source.begin(), source.end());
    } else {
      ring.assign(source.rbegin(), source.rend());
    }
    rings.push_back(std::move(ring));
  }

  // Repeated points (including a duplicated closing vertex) are not rejected
  // by ring_orientation2, which skips them when it picks the extreme vertex.
  // Check them after normalization because reversing a ring preserves them.
  const auto same_point = [](const Point2d& first,
                             const Point2d& second) noexcept {
    return first[0] == second[0] && first[1] == second[1];
  };
  for (const std::vector<Point2d>& ring : rings) {
    for (std::size_t i = 0; i < ring.size(); ++i) {
      if (same_point(ring[i], ring[(i + 1u) % ring.size()])) {
        return Error{Errc::degenerate,
                     "extrude_profile: profile ring has a repeated vertex"};
      }
    }
  }

  // Certify that the profile is a simple polygon with holes before any face is
  // built. Triangulation is the existing check that sees self-intersections and
  // holes that lie outside the outer ring. It runs on the normalized rings
  // because their order is also the cap order, so the caps reuse this
  // triangulation instead of triangulating the same polygon twice more; the
  // caps still keep the polygon-with-holes rings.
  Polygon2d certificate_profile;
  certificate_profile.ring_starts.reserve(rings.size() + 1u);
  certificate_profile.ring_starts.push_back(0u);
  for (const std::vector<Point2d>& ring : rings) {
    certificate_profile.vertices.insert(certificate_profile.vertices.end(),
                                        ring.begin(), ring.end());
    certificate_profile.ring_starts.push_back(
        static_cast<std::uint32_t>(certificate_profile.vertices.size()));
  }
  const Triangulation2 certificate = triangulate_polygon2(certificate_profile);
  if (!certificate.ok) {
    return Error{Errc::degenerate,
                 "extrude_profile: profile is not a simple polygon with holes"};
  }

  std::size_t point_count = 0;
  for (const std::vector<Point2d>& ring : rings) point_count += ring.size();

  std::vector<std::vector<Point3d>> bottom_rings(rings.size());
  std::vector<std::vector<Point3d>> top_rings(rings.size());
  for (std::size_t r = 0; r < rings.size(); ++r) {
    bottom_rings[r].reserve(rings[r].size());
    top_rings[r].reserve(rings[r].size());
    for (const Point2d& point : rings[r]) {
      const Point3d base = plane_frame_to_3d(frame, point);
      bottom_rings[r].push_back(base);
      top_rings[r].push_back(base + frame.normal * height);
    }
  }

  const Vec3d& normal = frame.normal;
  std::vector<AssembledFace> faces;
  faces.reserve(2u + point_count);

  // The top cap keeps the normalized order: outer counter-clockwise, holes
  // clockwise about +normal. The bottom cap reverses every ring, so it winds
  // the same way about -normal.
  std::vector<std::uint32_t> ring_sizes(rings.size(), 0u);
  for (std::size_t r = 0; r < rings.size(); ++r) {
    ring_sizes[r] = static_cast<std::uint32_t>(rings[r].size());
  }
  AssembledFace bottom_cap;
  AssembledFace top_cap;
  bottom_cap.normal = Vec3d{-normal[0], -normal[1], -normal[2]};
  top_cap.normal = normal;
  top_cap.triangulation = certificate;
  bottom_cap.triangulation =
      detail::extrude_reverse_triangulation(certificate, ring_sizes);
  for (std::size_t r = 0; r < rings.size(); ++r) {
    std::vector<Point3d> reversed(bottom_rings[r].rbegin(),
                                  bottom_rings[r].rend());
    if (r == 0u) {
      bottom_cap.outer = std::move(reversed);
      top_cap.outer = top_rings[r];
    } else {
      bottom_cap.holes.push_back(std::move(reversed));
      top_cap.holes.push_back(top_rings[r]);
    }
  }
  faces.push_back(std::move(bottom_cap));
  faces.push_back(std::move(top_cap));

  // One side quad per edge of every ring. For a hole the ring is clockwise and
  // cross(b - a, normal) points into the hole, i.e. out of the material.
  for (std::size_t r = 0; r < rings.size(); ++r) {
    const std::vector<Point3d>& bottom = bottom_rings[r];
    const std::vector<Point3d>& top = top_rings[r];
    for (std::size_t i = 0; i < rings[r].size(); ++i) {
      const std::size_t next = (i + 1u) % rings[r].size();
      AssembledFace side;
      side.outer = {bottom[i], bottom[next], top[next], top[i]};
      side.normal =
          detail::extrude_side_normal(bottom[i], bottom[next], normal);
      faces.push_back(std::move(side));
    }
  }

  Result<Brep> assembled = assemble_solid(faces, tolerance);
  if (!assembled) {
    return Error{Errc::invalid_topology,
                 "extrude_profile: solid could not be assembled"};
  }
  return assembled;
}

}  // namespace v1
}  // namespace geom
