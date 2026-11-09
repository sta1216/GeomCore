#pragma once

// Point classification against a closed planar B-rep solid.
//
// The solid boundary is treated as closed: a point within `tol` of any face,
// including a face interior, an edge, a vertex or a hole rim, is inside. A
// point that is not on the boundary is classified by ray parity. Each probe
// direction is tried in turn; a ray that hits a face boundary is ambiguous
// because the crossing may be counted twice or not at all, so that direction is
// discarded. When every direction is ambiguous the query fails instead of
// guessing.

#include <array>
#include <cmath>
#include <cstddef>
#include <span>
#include <vector>

#include <geom/brep/face_query.hpp>
#include <geom/brep/plane_frame.hpp>
#include <geom/brep/topology.hpp>
#include <geom/core/error.hpp>
#include <geom/core/result.hpp>
#include <geom/primitives/point.hpp>
#include <geom/surfaces/plane.hpp>

namespace geom {
inline namespace v1 {

namespace detail {

inline Result<bool> inside_check_point(const Point3d& p) {
  if (!std::isfinite(p[0]) || !std::isfinite(p[1]) || !std::isfinite(p[2])) {
    return Error{Errc::non_finite, "point_in_solid: point has a non-finite component"};
  }
  return true;
}

// A probe ray advances by `tol` per unit of `t`, so the tolerance must be
// strictly positive. Errc::degenerate covers zero, negative and non-finite
// tolerances; there is no meaningful distance threshold otherwise.
inline Result<bool> inside_check_tolerance(double tol) {
  if (!std::isfinite(tol) || !(tol > 0.0)) {
    return Error{Errc::degenerate, "point_in_solid: tolerance must be finite and positive"};
  }
  return true;
}

// Every face of every shell of the solid. An unknown solid, a dead shell, a
// solid without shells or a solid whose shells carry no faces is invalid
// topology rather than an outside classification.
inline Result<std::vector<FaceId>> inside_solid_faces(const Brep& brep, SolidId solid) {
  if (solid.is_null() || !brep.solid_alive(solid)) {
    return Error{Errc::invalid_topology, "point_in_solid: solid is not alive", solid.index(),
                 0u};
  }
  const std::span<const ShellId> shells = brep.solid_shells(solid);
  if (shells.empty()) {
    return Error{Errc::invalid_topology, "point_in_solid: solid has no shells", solid.index(),
                 0u};
  }
  std::vector<FaceId> faces;
  for (const ShellId shell : shells) {
    if (shell.is_null() || !brep.shell_alive(shell)) {
      return Error{Errc::invalid_topology, "point_in_solid: shell is not alive", shell.index(),
                   0u};
    }
    const std::span<const FaceId> shell_faces = brep.shell_faces(shell);
    faces.insert(faces.end(), shell_faces.begin(), shell_faces.end());
  }
  if (faces.empty()) {
    return Error{Errc::invalid_topology, "point_in_solid: solid has no faces", solid.index(),
                 0u};
  }
  return faces;
}

}  // namespace detail

// Classifies `p` against `solid` with distance tolerance `tol`.
//
// The solid is closed: a point on any face, including an edge, a vertex or a
// hole rim (within `tol`), is inside. Other points are classified by the parity
// of ray crossings. A probe direction whose face hit is within `tol` of that
// face's boundary is ambiguous, so the next fixed direction is tried.
//
// Errors, in check order:
//   - non-finite point                         -> Errc::non_finite
//   - tol <= 0, NaN or infinite                -> Errc::degenerate
//   - unknown/dead solid, no shells, no faces  -> Errc::invalid_topology
//   - every probe direction ambiguous          -> Errc::numeric_failure
inline Result<bool> point_in_solid(const Brep& brep, SolidId solid, const Point3d& p,
                                   double tol) {
  const Result<bool> point_check = detail::inside_check_point(p);
  if (!point_check) return point_check.error();
  const Result<bool> tol_check = detail::inside_check_tolerance(tol);
  if (!tol_check) return tol_check.error();
  const Result<std::vector<FaceId>> faces_result = detail::inside_solid_faces(brep, solid);
  if (!faces_result) return faces_result.error();
  const std::vector<FaceId>& faces = faces_result.value();

  // Closed-solid convention: every point on a face, including its boundary,
  // belongs to the solid. Checking face interiors as well as rims keeps the
  // answer independent of the probe direction.
  for (const FaceId face : faces) {
    const Result<bool> interior = point_in_face(brep, face, p, tol);
    if (!interior) return interior.error();
    if (interior.value()) return true;
    const Result<bool> boundary = point_on_face_boundary(brep, face, p, tol);
    if (!boundary) return boundary.error();
    if (boundary.value()) return true;
  }

  const Vec3d point{p[0], p[1], p[2]};

  // The first direction is the box diagonal. On a box it can leave through a
  // corner or an edge, which exercises the retry path. The remaining three
  // tetrahedral directions are not axis aligned.
  constexpr std::array<Vec3d, 4> kRawDirections = {{
      Vec3d{1.0, 1.0, 1.0},
      Vec3d{0.5774, -0.5775, -0.5773},
      Vec3d{-0.5775, 0.5773, -0.5774},
      Vec3d{-0.5773, -0.5774, 0.5775},
  }};

  for (const Vec3d& raw : kRawDirections) {
    const double length = std::sqrt(dot(raw, raw));
    if (!(std::isfinite(length) && length > 0.0)) continue;
    const Vec3d direction{raw[0] / length, raw[1] / length, raw[2] / length};
    int parity = 0;
    bool ambiguous = false;
    for (const FaceId face : faces) {
      const Result<Plane3> plane_result = face_plane(brep, face);
      if (!plane_result) return plane_result.error();
      const Result<PlaneFrame> frame_result = make_plane_frame(plane_result.value());
      if (!frame_result) return frame_result.error();
      const Plane3 plane = plane_frame_plane(frame_result.value());
      const double den = dot(plane.normal, direction);
      if (!(std::fabs(den) > 1e-12)) continue;
      const double t = (plane.offset - dot(plane.normal, point)) / den;
      if (!(t > tol)) continue;
      const Point3d hit = p + direction * t;
      if (!std::isfinite(hit[0]) || !std::isfinite(hit[1]) || !std::isfinite(hit[2])) {
        continue;
      }
      const Result<bool> boundary = point_on_face_boundary(brep, face, hit, tol);
      if (!boundary) return boundary.error();
      if (boundary.value()) {
        ambiguous = true;
        break;
      }
      const Result<bool> inside = point_in_face(brep, face, hit, tol);
      if (!inside) return inside.error();
      if (inside.value()) ++parity;
    }
    if (ambiguous) continue;
    return parity % 2 == 1;
  }

  return Error{Errc::numeric_failure, "point_in_solid: every probe ray is degenerate",
               solid.index(), 0u};
}

}  // namespace v1
}  // namespace geom
