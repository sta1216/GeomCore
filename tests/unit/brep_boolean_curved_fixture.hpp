#pragma once

// Shared measurement helpers for the curved boolean tests. `sphere_face_volume`
// sums the signed chart volumes of a solid's spherical faces; each face must be
// a polar cap (a whole-sphere boolean result is exactly two such caps sharing
// one rim). The QA-16 corpus and the QA-17 arbitrary-placement lens both use it,
// so it lives here rather than in either translation unit.

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <numbers>
#include <optional>
#include <variant>
#include <vector>

#include <geom/brep/face_uv.hpp>
#include <geom/brep/topology.hpp>
#include <geom/primitives/point.hpp>
#include <geom/surfaces/analytic.hpp>

namespace geom::test {

// Signed chart volume of one spherical face, or nullopt when the face is not
// spherical or is not a polar cap in its own chart. The face must be a polar cap:
// its single uv loop has to reach a pole (v = +-pi/2), and the latitude range
// gives the cap height. The loop's signed shoelace area carries the orientation
// of the ring itself (+ outward, - reversed). The result does not consult the
// face's `reversed` flag; callers that need the flag's sign apply it themselves.
inline std::optional<double> sphere_face_cap_volume(const Brep& brep, FaceId face_id) {
  constexpr double kTolerance = 1e-9;
  const Face& face = brep.face(face_id);
  const SurfaceGeometry& surface = brep.geometry().surfaces[face.surface];
  const Sphere3* sphere = std::get_if<Sphere3>(&surface);
  if (sphere == nullptr) return std::nullopt;
  const Result<std::vector<FaceUvLoop>> loops =
      face_loops_uv(brep, face_id, kTolerance);
  if (!loops.has_value() || loops.value().size() != 1) return std::nullopt;
  const std::vector<Point2d>& ring = loops.value().front().points;
  if (ring.size() < 3) return std::nullopt;

  double twice_area = 0.0;
  double v_lo = 1e30;
  double v_hi = -1e30;
  for (std::size_t i = 0; i < ring.size(); ++i) {
    const Point2d& point = ring[i];
    const Point2d& next = ring[(i + 1) % ring.size()];
    twice_area += point[0] * next[1] - next[0] * point[1];
    v_lo = std::min(v_lo, point[1]);
    v_hi = std::max(v_hi, point[1]);
  }
  if (twice_area == 0.0) return std::nullopt;

  constexpr double kHalfPi = 0.5 * std::numbers::pi_v<double>;
  // An off-axis chart recovers the pole through asin(dot(offset, axis) / r),
  // which amplifies the |axis|^2 = 1 +/- 1e-16 rounding to about 2e-8 in v.
  // Accept poles within 1e-6; the induced cap-height error is O(1e-12), far
  // below the 1e-6 volume assertions.
  constexpr double kPoleTolerance = 1e-6;
  double height = 0.0;
  if (std::abs(v_lo + kHalfPi) <= kPoleTolerance) {
    // Southern cap: the pole is at v = -pi/2, so its axial extent is
    // r * (sin(v_hi) + 1).
    height = sphere->radius * (1.0 + std::sin(v_hi));
  } else if (std::abs(v_hi - kHalfPi) <= kPoleTolerance) {
    // Northern cap: r * (1 - sin(v_lo)).
    height = sphere->radius * (1.0 - std::sin(v_lo));
  } else {
    return std::nullopt;
  }
  const double sign = twice_area > 0.0 ? 1.0 : -1.0;
  const double radius = sphere->radius;
  return sign * std::numbers::pi_v<double> * height * height *
         (3.0 * radius - height) / 3.0;
}

// Signed chart volume of a solid built only from spherical faces. Each spherical
// face must be a polar cap; `sphere_face_cap_volume` owns the per-face sign and
// the cap geometry. A face without a pole at either v end is not a polar cap and
// is reported as nullopt rather than as zero.
inline std::optional<double> sphere_face_volume(const Brep& brep) {
  double volume = 0.0;
  bool any_sphere = false;
  for (const FaceId face_id : brep.face_ids()) {
    const Face& face = brep.face(face_id);
    const SurfaceGeometry& surface = brep.geometry().surfaces[face.surface];
    if (!std::holds_alternative<Sphere3>(surface)) continue;
    const std::optional<double> cap = sphere_face_cap_volume(brep, face_id);
    if (!cap.has_value()) return std::nullopt;
    volume += *cap;
    any_sphere = true;
  }
  if (!any_sphere) return std::nullopt;
  return volume;
}

}  // namespace geom::test
