#pragma once

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <utility>
#include <vector>

#include <geom/brep/topology.hpp>
#include <geom/core/error.hpp>
#include <geom/core/result.hpp>

namespace geom {
inline namespace v1 {

// Polyline pcurves in a face's parameter domain.
//
// A pcurve runs along its edge from t0 to t1 (not along the coedge's traversal
// direction, which may be reversed), so pcurve_split() cuts it at the fraction of
// its chord length and both halves carry the same computed junction point.
// Degenerate input is rejected instead of silently producing a wrong split:
// fewer than two points and zero total length are Errc::degenerate, non-finite
// coordinates (or a non-finite fraction) are Errc::non_finite, and so is a chord
// sum that overflows to infinity.

// True when the polyline has at least two points and every coordinate is finite.
// A pcurve with coincident consecutive points is still well formed; it merely
// contributes zero length there.
inline bool pcurve_well_formed(const Pcurve2d& pcurve) noexcept {
  if (pcurve.points.size() < 2u) return false;
  for (const Point2d& point : pcurve.points) {
    if (!std::isfinite(point[0]) || !std::isfinite(point[1])) return false;
  }
  return true;
}

namespace detail {

inline bool point2d_is_finite(const Point2d& point) noexcept {
  return std::isfinite(point[0]) && std::isfinite(point[1]);
}

inline bool points_are_equal(const Point2d& a, const Point2d& b) noexcept {
  return a[0] == b[0] && a[1] == b[1];
}

inline double chord_length(const Pcurve2d& pcurve, std::size_t segment) noexcept {
  return std::hypot(pcurve.points[segment + 1u][0] - pcurve.points[segment][0],
                    pcurve.points[segment + 1u][1] - pcurve.points[segment][1]);
}

// Splits the polyline at the normalized chord-length fraction `t`, clamped to
// [0, 1], and returns both halves. The junction point is computed once and stored
// in both halves, so pcurve_merge() can drop the duplicate by exact comparison.
// At a split that lands on an existing vertex the vertex is not duplicated, and a
// clamped endpoint yields a degenerate (zero-length, two-point) other half.
inline Result<std::pair<Pcurve2d, Pcurve2d>> pcurve_split_halves(const Pcurve2d& pcurve,
                                                                 double t) {
  const std::size_t count = pcurve.points.size();
  if (count < 2u) {
    return Error{Errc::degenerate, "pcurve2d: pcurve has fewer than two points", 0u, 0u};
  }
  for (std::size_t i = 0; i < count; ++i) {
    if (!point2d_is_finite(pcurve.points[i])) {
      return Error{Errc::non_finite, "pcurve2d: pcurve has a non-finite point",
                   static_cast<std::uint32_t>(i), 0u};
    }
  }
  if (!std::isfinite(t)) {
    return Error{Errc::non_finite, "pcurve2d: split fraction is not finite", 0u, 0u};
  }

  double total = 0.0;
  for (std::size_t i = 0; i + 1u < count; ++i) total += chord_length(pcurve, i);
  // Every point is finite, but the chord sum can still overflow to +inf; the
  // fraction would then be NaN and the halves silently wrong.
  if (!std::isfinite(total)) {
    return Error{Errc::non_finite, "pcurve_split: chord length is not finite", 0u, 0u};
  }
  if (!(total > 0.0)) {
    return Error{Errc::degenerate, "pcurve2d: pcurve has zero length", 0u, 0u};
  }
  const double target = std::clamp(t, 0.0, 1.0) * total;

  // Walk the chords, skipping zero-length ones (they cannot host the junction).
  // The default is the end of the polyline: it is the answer when rounding leaves
  // `target` equal to the running total at the last segment.
  std::size_t segment = count - 2u;
  double fraction = 1.0;
  double accumulated = 0.0;
  for (std::size_t i = 0; i + 1u < count; ++i) {
    const double length = chord_length(pcurve, i);
    if (!(length > 0.0)) continue;
    if (target <= accumulated + length) {
      segment = i;
      fraction = (target - accumulated) / length;
      break;
    }
    accumulated += length;
  }

  // Snap to an existing vertex when the fraction lands on an endpoint, so the
  // junction is bit-identical to the vertex and a merge can dedupe it exactly.
  Point2d junction = pcurve.points[segment];
  if (fraction >= 1.0) {
    junction = pcurve.points[segment + 1u];
  } else if (fraction > 0.0) {
    junction[0] = pcurve.points[segment][0] +
                  (pcurve.points[segment + 1u][0] - pcurve.points[segment][0]) * fraction;
    junction[1] = pcurve.points[segment][1] +
                  (pcurve.points[segment + 1u][1] - pcurve.points[segment][1]) * fraction;
  }

  Pcurve2d lower;
  lower.points.reserve(segment + 2u);
  lower.points.insert(lower.points.end(), pcurve.points.begin(),
                      pcurve.points.begin() + static_cast<std::ptrdiff_t>(segment) + 1);
  if (!points_are_equal(junction, lower.points.back())) lower.points.push_back(junction);
  // A split at the start clamps to a zero-length lower half, which still needs
  // two points to stay well formed.
  if (lower.points.size() < 2u) lower.points.push_back(junction);

  Pcurve2d upper;
  upper.points.reserve(count - segment);
  if (!points_are_equal(junction, pcurve.points[segment + 1u])) {
    upper.points.push_back(junction);
  }
  upper.points.insert(upper.points.end(),
                      pcurve.points.begin() + static_cast<std::ptrdiff_t>(segment) + 1,
                      pcurve.points.end());
  if (upper.points.size() < 2u) upper.points.push_back(junction);

  return std::make_pair(std::move(lower), std::move(upper));
}

}  // namespace detail

// Total chord length of the polyline; 0 for a malformed pcurve and +inf when the
// chord sum overflows. Allocation-free and noexcept by contract.
inline double pcurve_length(const Pcurve2d& pcurve) noexcept {
  if (!pcurve_well_formed(pcurve)) return 0.0;
  double length = 0.0;
  for (std::size_t i = 0; i + 1u < pcurve.points.size(); ++i) {
    length += detail::chord_length(pcurve, i);
  }
  return length;
}

// Lower half of the polyline up to the normalized chord-length fraction `t`,
// clamped to [0, 1]. The returned polyline ends at the junction point; call
// detail::pcurve_split_halves() when both halves are needed.
inline Result<Pcurve2d> pcurve_split(const Pcurve2d& pcurve, double t) {
  Result<std::pair<Pcurve2d, Pcurve2d>> halves = detail::pcurve_split_halves(pcurve, t);
  if (!halves.has_value()) return halves.error();
  return std::move(halves.value().first);
}

// Concatenates two halves that meet at a junction. When the lower half's last
// point and the upper half's first point are exactly equal, the duplicate
// junction is dropped once; a non-equal junction is kept on both sides.
inline Pcurve2d pcurve_merge(const Pcurve2d& lower, const Pcurve2d& upper) {
  Pcurve2d merged;
  merged.points.reserve(lower.points.size() + upper.points.size());
  merged.points.insert(merged.points.end(), lower.points.begin(), lower.points.end());
  std::size_t first = 0;
  if (!merged.points.empty() && !upper.points.empty() &&
      detail::points_are_equal(merged.points.back(), upper.points.front())) {
    first = 1;
  }
  merged.points.insert(merged.points.end(),
                       upper.points.begin() + static_cast<std::ptrdiff_t>(first),
                       upper.points.end());
  return merged;
}

}  // namespace v1
}  // namespace geom
