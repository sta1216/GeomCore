#pragma once

// Collapsing a known circular cut chain into exact arc spans.
//
// SSI returns cut curves as polylines whose vertices are Newton-refined onto
// both surfaces. When the cut is analytically known to be a circle (for example
// a cylinder cut perpendicular to its axis, or a plane cut of a sphere), the
// polyline can be replaced by one or two exact `CurvedArcSpan` pieces.
// `assemble_curved_solid` then builds exact `CircleCurve3` edges instead of
// chords, which keeps the volume of a split or boolean result exact.
//
// Only segments that carry no arc provenance yet are touched; original boundary
// arcs are left alone. A full-turn run is split into two spans because the
// assembler cannot recover the traversal direction of a single 2*pi span (its
// endpoints coincide).

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <optional>
#include <utility>
#include <vector>

#include <geom/brep/trim_face.hpp>
#include <geom/core/error.hpp>
#include <geom/curves/analytic_curve.hpp>
#include <geom/primitives/point.hpp>
#include <geom/primitives/vec.hpp>

namespace geom {
inline namespace v1 {
namespace detail {
namespace arc_heal_detail {

inline constexpr double kArcHealPi = 3.14159265358979323846;

// Largest angular step accepted between two consecutive vertices of a healed
// chain (28.6 degrees). Cut circles are sampled at 1.4 to 3.2 degrees, so a
// segment spanning more than this between two on-circle vertices is a straight
// chord of the face boundary (for example where a clipped planar fragment meets
// a side face), not a sample of the circle.
inline constexpr double kArcHealMaxStep = 0.5;

inline double arc_heal_wrap(double delta) noexcept {
  const double turn = 2.0 * kArcHealPi;
  delta = std::fmod(delta, turn);
  if (delta <= -kArcHealPi) delta += turn;
  if (delta > kArcHealPi) delta -= turn;
  return delta;
}

inline Vec3d arc_heal_normal(const CircleCurve3& circle) noexcept {
  return cross(circle.x_axis, circle.y_axis);
}

inline double arc_heal_angle(const CircleCurve3& circle, const Point3d& point) noexcept {
  const Vec3d offset = point - circle.center;
  return std::atan2(dot(offset, circle.y_axis), dot(offset, circle.x_axis));
}

struct ArcHealLoop {
  std::vector<Point3d> points;
  std::vector<Point2d> uv;
  std::vector<std::optional<CurvedArcSpan>> spans;
};

// Builds the span for the vertex run [first, last] (inclusive). The sweep
// follows the loop order and is the sum of the wrapped deltas.
inline Result<CurvedArcSpan> arc_heal_make_span(const CircleCurve3& circle,
                                                const std::vector<Point3d>& points,
                                                std::size_t first, std::size_t last,
                                                std::size_t size, double angle_tolerance) {
  double sweep = 0.0;
  int sign = 0;
  std::size_t index = first;
  while (index != last) {
    const std::size_t next = (index + 1u) % size;
    const double delta = arc_heal_wrap(arc_heal_angle(circle, points[next]) -
                                       arc_heal_angle(circle, points[index]));
    int step_sign = 0;
    if (delta > angle_tolerance) {
      step_sign = 1;
    } else if (delta < -angle_tolerance) {
      step_sign = -1;
    } else {
      return Error{Errc::invalid_topology, "collapse_circular_chains: arc run backtracks", 0u, 0u};
    }
    if (sign != 0 && sign != step_sign) {
      return Error{Errc::invalid_topology, "collapse_circular_chains: arc run backtracks", 0u, 0u};
    }
    sign = step_sign;
    sweep += delta;
    index = next;
  }
  if (!(std::abs(sweep) > 0.0) || std::abs(sweep) >= 2.0 * kArcHealPi - angle_tolerance) {
    return Error{Errc::invalid_topology, "collapse_circular_chains: arc run is a full turn", 0u,
                 0u};
  }
  const double t0 = arc_heal_angle(circle, points[first]);
  return CurvedArcSpan{circle, t0, t0 + sweep};
}

inline Result<ArcHealLoop> arc_heal_loop(const std::vector<Point3d>& points,
                                         const std::vector<Point2d>& uv,
                                         const std::vector<std::optional<CurvedArcSpan>>& spans,
                                         const CircleCurve3& circle, double tolerance) {
  const std::size_t size = points.size();
  if (size < 3) {
    return ArcHealLoop{points, uv, {}};
  }
  const double slack = std::max(8.0 * tolerance, 1e-9 * circle.radius);
  const double angle_tolerance =
      std::min(kArcHealPi, std::max(1e-12, 2.0 * tolerance / circle.radius));
  std::vector<bool> on(size, false);
  for (std::size_t i = 0; i < size; ++i) {
    const Vec3d offset = points[i] - circle.center;
    const Vec3d normal = arc_heal_normal(circle);
    const double axial = dot(offset, normal);
    const Vec3d in_plane = offset - normal * axial;
    on[i] = std::abs(axial) <= slack &&
            std::abs(std::sqrt(dot(in_plane, in_plane)) - circle.radius) <= slack;
  }

  // A run continues along a segment only when that segment samples the circle.
  // Two vertices can both lie on the circle while the segment joining them is a
  // straight chord of the face boundary (where a clipped planar fragment meets
  // a side face, for example); healing such a chord as an arc would bend the
  // boundary. Sampling is far finer than kArcHealMaxStep, so a longer step
  // between two on-circle vertices marks a chord.
  std::vector<bool> arc_edge(size, false);
  for (std::size_t i = 0; i < size; ++i) {
    const std::size_t next = (i + 1u) % size;
    if (!on[i] || !on[next]) continue;
    const double delta = arc_heal_wrap(arc_heal_angle(circle, points[next]) -
                                       arc_heal_angle(circle, points[i]));
    arc_edge[i] = std::abs(delta) <= kArcHealMaxStep;
  }

  std::size_t anchor = size;
  for (std::size_t i = 0; i < size; ++i) {
    if (!on[i]) {
      anchor = i;
      break;
    }
  }
  if (anchor == size) {
    for (std::size_t i = 0; i < size; ++i) {
      if (!arc_edge[i]) {
        anchor = i;
        break;
      }
    }
  }

  if (anchor == size) {
    // The whole loop lies on the circle: split at the vertex nearest a half turn.
    double half = 0.0;
    std::size_t split = size;
    double best = 1e300;
    for (std::size_t i = 1; i < size; ++i) {
      half += arc_heal_wrap(arc_heal_angle(circle, points[i]) -
                            arc_heal_angle(circle, points[i - 1]));
      if (std::abs(std::abs(half) - kArcHealPi) < best) {
        best = std::abs(std::abs(half) - kArcHealPi);
        split = i;
      }
    }
    if (split == size || size < 4) {
      return Error{Errc::invalid_topology, "collapse_circular_chains: no half-turn vertex", 0u, 0u};
    }
    Result<CurvedArcSpan> first =
        arc_heal_make_span(circle, points, 0u, split, size, angle_tolerance);
    if (!first) return first.error();
    Result<CurvedArcSpan> second =
        arc_heal_make_span(circle, points, split, 0u, size, angle_tolerance);
    if (!second) return second.error();
    ArcHealLoop out;
    out.points = {points[0], points[split]};
    out.uv = {uv[0], uv[split]};
    out.spans = {first.value(), second.value()};
    return out;
  }

  ArcHealLoop out;
  out.points.reserve(size);
  out.uv.reserve(size);
  out.spans.reserve(size);

  std::size_t index = 0;
  std::size_t current = anchor;
  while (index < size) {
    const std::size_t begin = current;
    const bool engaged = begin < spans.size() && spans[begin].has_value();
    if (on[begin] && arc_edge[begin] && !engaged) {
      // A run may wrap past the end of the ring; the ring has no natural start,
      // so the walk begins at a chord and must be allowed to close the loop.
      std::size_t length = 1;
      while (length < size) {
        const std::size_t next = (begin + length) % size;
        if (!on[next]) break;
        if (!arc_edge[(begin + length - 1u) % size]) break;
        if (next < spans.size() && spans[next].has_value()) break;
        ++length;
      }
      if (length >= 3) {
        const std::size_t end = (begin + length - 1u) % size;
        Result<CurvedArcSpan> span =
            arc_heal_make_span(circle, points, begin, end, size, angle_tolerance);
        if (!span) return span.error();
        out.points.push_back(points[begin]);
        out.uv.push_back(uv[begin]);
        out.spans.push_back(span.value());
        index += length - 1u;
        current = end;
        if (index >= size) break;
        continue;
      }
    }
    out.points.push_back(points[begin]);
    out.uv.push_back(uv[begin]);
    out.spans.push_back(engaged ? spans[begin] : std::nullopt);
    ++index;
    current = (begin + 1u) % size;
  }
  return out;
}

}  // namespace arc_heal_detail

// Replaces unengaged co-circular chains of `fragment` with exact arc spans.
// Returns true when at least one chain was collapsed.
inline Result<bool> collapse_circular_chains(CurvedFragment& fragment, const CircleCurve3& circle,
                                             double tolerance) {
  if (!circle.is_valid() || !(tolerance > 0.0) || !std::isfinite(tolerance)) {
    return Error{Errc::degenerate, "collapse_circular_chains: invalid input", 0u, 0u};
  }
  bool changed = false;
  for (std::size_t loop = 0; loop < fragment.loops.size(); ++loop) {
    const std::size_t size = fragment.loops[loop].size();
    if (size < 3) continue;
    if (loop >= fragment.loops_uv.size() || fragment.loops_uv[loop].size() != size) {
      return Error{Errc::invalid_topology, "collapse_circular_chains: uv ring mismatch", 0u, 0u};
    }
    std::vector<std::optional<CurvedArcSpan>> spans;
    if (loop < fragment.arc_spans.size()) spans = fragment.arc_spans[loop];
    Result<arc_heal_detail::ArcHealLoop> healed = arc_heal_detail::arc_heal_loop(
        fragment.loops[loop], fragment.loops_uv[loop], spans, circle, tolerance);
    if (!healed) return healed.error();
    if (healed.value().points.size() == size) continue;
    changed = true;
    fragment.loops[loop] = std::move(healed.value().points);
    fragment.loops_uv[loop] = std::move(healed.value().uv);
    if (fragment.arc_spans.size() < fragment.loops.size()) {
      fragment.arc_spans.resize(fragment.loops.size());
    }
    fragment.arc_spans[loop] = std::move(healed.value().spans);
  }
  return changed;
}

}  // namespace detail
}  // namespace v1
}  // namespace geom
