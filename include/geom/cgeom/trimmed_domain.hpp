#pragma once

#include <geom/brep/face_uv.hpp>
#include <geom/cgeom/polygon2.hpp>
#include <geom/cgeom/triangulate2.hpp>
#include <geom/core/result.hpp>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <optional>
#include <span>
#include <utility>
#include <vector>

namespace geom {
inline namespace v1 {

// A trimmed face boundary ring in the face's uv parameterization.
//
// The ring is implicitly closed: `points` does not repeat the first point.
// `arcs` is parallel to `points`; arcs[i] describes the segment from points[i]
// to points[(i + 1) % points.size()] and is disengaged when the segment is a
// straight edge (or a numerical correction of a resampled arc endpoint).
struct TrimmedBoundary {
  std::vector<Point2d> points;
  std::vector<std::optional<FaceUvArcRef>> arcs;
  bool outer = false;
};

// A triangulation of the trimmed uv domain. `points` holds every ring vertex,
// outer ring first, and `triangles` indexes into it.
struct TrimmedUvTriangulation {
  std::vector<Point2d> points;
  std::vector<std::array<std::uint32_t, 3>> triangles;
};

namespace detail {

inline constexpr double kTrimmedPi = 3.141592653589793238462643383279502884;
inline constexpr double kTrimmedTwoPi = 2.0 * kTrimmedPi;
inline constexpr std::size_t kTrimmedMaxSamples = 1000000u;
inline constexpr double kTrimmedTolerance = 1e-12;

inline bool trimmed_point_finite(const Point2d& p) noexcept {
  return std::isfinite(p[0]) && std::isfinite(p[1]);
}

inline bool trimmed_points_close(const Point2d& a, const Point2d& b) noexcept {
  return std::abs(a[0] - b[0]) <= kTrimmedTolerance &&
         std::abs(a[1] - b[1]) <= kTrimmedTolerance;
}

// Slid by a multiple of 2*pi so the result is within pi of `reference`. The
// sampled arc points come back in (-pi, pi] from the cylinder/cone/torus
// projections, so each sample is unwrapped against the previous ring point.
inline double trimmed_unwrap_u(double value, double reference) noexcept {
  return value + kTrimmedTwoPi * std::round((reference - value) / kTrimmedTwoPi);
}

inline bool trimmed_same_circle(const CircleCurve3& a, const CircleCurve3& b) noexcept {
  const double eps = 1e-12;
  return std::abs(a.radius - b.radius) <= eps &&
         std::abs(a.center[0] - b.center[0]) <= eps &&
         std::abs(a.center[1] - b.center[1]) <= eps &&
         std::abs(a.center[2] - b.center[2]) <= eps &&
         std::abs(a.x_axis[0] - b.x_axis[0]) <= eps &&
         std::abs(a.x_axis[1] - b.x_axis[1]) <= eps &&
         std::abs(a.x_axis[2] - b.x_axis[2]) <= eps &&
         std::abs(a.y_axis[0] - b.y_axis[0]) <= eps &&
         std::abs(a.y_axis[1] - b.y_axis[1]) <= eps &&
         std::abs(a.y_axis[2] - b.y_axis[2]) <= eps;
}

// A maximal run of segments that share one analytic circle and are contiguous
// in its angle parameter. Merging the upstream (usually very fine) chain keeps
// the resampled count tied to the requested deflection instead of the
// tessellation that produced the face.
struct TrimmedRingRun {
  std::size_t first_segment = 0;
  std::size_t last_segment = 0;
  bool arc = false;
  CircleCurve3 circle{};
  double t0 = 0.0;
  double t1 = 0.0;
};

inline std::vector<TrimmedRingRun> trimmed_ring_runs(
    const std::vector<std::optional<FaceUvArcRef>>& arcs) {
  const std::size_t n = arcs.size();
  std::vector<TrimmedRingRun> runs;
  runs.reserve(n);
  std::size_t i = 0;
  while (i < n) {
    TrimmedRingRun run;
    run.first_segment = i;
    if (!arcs[i]) {
      run.last_segment = i;
      runs.push_back(run);
      ++i;
      continue;
    }
    run.arc = true;
    run.circle = arcs[i]->circle;
    run.t0 = arcs[i]->t0;
    double t_end = arcs[i]->t1;
    std::size_t last = i;
    while (last + 1 < n && arcs[last + 1] &&
           trimmed_same_circle(run.circle, arcs[last + 1]->circle)) {
      const double next_t0 = trimmed_unwrap_u(arcs[last + 1]->t0, t_end);
      if (std::abs(next_t0 - t_end) > 1e-9) break;
      t_end = trimmed_unwrap_u(arcs[last + 1]->t1, next_t0);
      ++last;
    }
    run.last_segment = last;
    run.t1 = t_end;
    runs.push_back(run);
    i = last + 1;
  }
  return runs;
}

// Half of the angular step that keeps a chord within `deflection` of a circle
// of the given radius: 2 * acos(1 - deflection / radius), clamped.
inline double trimmed_arc_step(double radius, double deflection) noexcept {
  if (!(radius > 0.0) || !std::isfinite(radius)) return 0.0;
  double ratio = 1.0 - deflection / radius;
  ratio = std::clamp(ratio, -1.0, 1.0);
  return 2.0 * std::acos(ratio);
}

inline Error trimmed_cap_error(FaceId face) {
  return Error{Errc::numeric_failure,
               "trimmed_face_boundaries: sample cap exceeded on face", face.index(), 0u};
}

// Planned number of output points for one ring, computed before any ring
// storage is allocated. Each arc run contributes its segment samples plus one
// exact-endpoint correction at most; the closing duplicate is dropped later.
inline Result<std::size_t> trimmed_plan_ring(const std::vector<TrimmedRingRun>& runs,
                                             double deflection, FaceId face) {
  std::size_t total = 1u;  // the ring start point
  for (const TrimmedRingRun& run : runs) {
    if (!run.arc) {
      ++total;
      if (total > kTrimmedMaxSamples) return trimmed_cap_error(face);
      continue;
    }
    if (!std::isfinite(run.t0) || !std::isfinite(run.t1)) {
      return Error{Errc::numeric_failure, "trimmed_face_boundaries: arc parameter is not finite",
                   face.index(), 0u};
    }
    const double step = trimmed_arc_step(run.circle.radius, deflection);
    const double span = std::abs(run.t1 - run.t0);
    const double segments = std::ceil(span / step);
    if (!(segments >= 0.0) || !std::isfinite(segments) ||
        segments > static_cast<double>(kTrimmedMaxSamples)) {
      return trimmed_cap_error(face);
    }
    total += static_cast<std::size_t>(segments) + 1u;
    if (total > kTrimmedMaxSamples) return trimmed_cap_error(face);
  }
  return total;
}

inline void trimmed_reverse_ring(TrimmedBoundary& ring) {
  const std::size_t n = ring.points.size();
  std::vector<std::optional<FaceUvArcRef>> reversed(n);
  for (std::size_t j = 0; j < n; ++j) {
    const std::size_t source = (n - 2u - j + n) % n;
    const std::optional<FaceUvArcRef>& arc = ring.arcs[source];
    if (arc) {
      reversed[j] = FaceUvArcRef{arc->circle, arc->t1, arc->t0};
    }
  }
  std::reverse(ring.points.begin(), ring.points.end());
  ring.arcs = std::move(reversed);
}

inline Result<TrimmedBoundary> trimmed_build_ring(const FaceUvLoop& loop,
                                                  const FaceUvProjection& projection,
                                                  double deflection,
                                                  const std::vector<TrimmedRingRun>& runs,
                                                  std::size_t planned, bool outer,
                                                  FaceId face) {
  const std::size_t n = loop.points.size();
  if (n < 3 || loop.arcs.size() != n) {
    return Error{Errc::invalid_topology,
                 "trimmed_face_boundaries: loop buffers are inconsistent", face.index(), 0u};
  }

  TrimmedBoundary ring;
  ring.outer = outer;
  ring.points.reserve(planned + 1u);
  ring.arcs.reserve(planned);
  ring.points.push_back(loop.points[0]);
  double previous_u = loop.points[0][0];

  for (const TrimmedRingRun& run : runs) {
    if (!run.arc) {
      const Point2d end = loop.points[(run.last_segment + 1u) % n];
      ring.points.push_back(end);
      ring.arcs.push_back(std::nullopt);
      previous_u = end[0];
      continue;
    }

    if (!(run.circle.radius > 0.0) || !std::isfinite(run.circle.radius)) {
      return Error{Errc::numeric_failure, "trimmed_face_boundaries: arc radius is not positive",
                   face.index(), 0u};
    }
    const double step = trimmed_arc_step(run.circle.radius, deflection);
    const double span = std::abs(run.t1 - run.t0);
    const double segments = std::ceil(span / step);
    if (!(segments >= 0.0) || !std::isfinite(segments) ||
        segments > static_cast<double>(kTrimmedMaxSamples)) {
      return trimmed_cap_error(face);
    }
    const std::size_t count = static_cast<std::size_t>(segments);
    double t_previous = run.t0;
    for (std::size_t k = 1; k <= count; ++k) {
      const double fraction = static_cast<double>(k) / static_cast<double>(count);
      const double t = run.t0 + (run.t1 - run.t0) * fraction;
      const Result<Point2d> projected = face_uv_project(projection, run.circle.eval(t),
                                                        previous_u, true);
      if (!projected) return projected.error();
      Point2d sample = projected.value();
      sample[0] = trimmed_unwrap_u(sample[0], previous_u);
      ring.points.push_back(sample);
      ring.arcs.push_back(FaceUvArcRef{run.circle, t_previous, t});
      previous_u = sample[0];
      t_previous = t;
    }

    const Point2d end = loop.points[(run.last_segment + 1u) % n];
    if (!trimmed_points_close(ring.points.back(), end)) {
      ring.points.push_back(end);
      ring.arcs.push_back(std::nullopt);
      previous_u = end[0];
    }
  }

  if (ring.points.size() < 2u || !trimmed_points_close(ring.points.back(), ring.points.front())) {
    return Error{Errc::numeric_failure, "trimmed_face_boundaries: ring did not close",
                 face.index(), 0u};
  }
  ring.points.pop_back();  // implicitly closed ring: drop the repeated first point

  // Drop consecutive duplicates; the arc that encoded the zero-length segment
  // is dropped with the duplicate point so arcs stays parallel to points.
  std::size_t i = 1u;
  while (i < ring.points.size()) {
    if (trimmed_points_close(ring.points[i], ring.points[i - 1u])) {
      ring.points.erase(ring.points.begin() + static_cast<std::ptrdiff_t>(i));
      ring.arcs.erase(ring.arcs.begin() + static_cast<std::ptrdiff_t>(i - 1u));
    } else {
      ++i;
    }
  }
  if (ring.points.size() < 3u) {
    return Error{Errc::invalid_topology,
                 "trimmed_face_boundaries: ring has fewer than three points", face.index(), 0u};
  }
  if (geom::signed_area2(ring.points) == 0.0) {
    return Error{Errc::degenerate, "trimmed_face_boundaries: ring has zero area", face.index(),
                 0u};
  }

  // Re-unwrap the ring so no consecutive u step jumps across the periodic seam.
  for (std::size_t k = 1u; k < ring.points.size(); ++k) {
    ring.points[k][0] = trimmed_unwrap_u(ring.points[k][0], ring.points[k - 1u][0]);
  }

  const double area = geom::signed_area2(ring.points);
  if ((outer && area < 0.0) || (!outer && area > 0.0)) trimmed_reverse_ring(ring);

  return ring;
}

}  // namespace detail

// Extracts the trimmed uv boundary rings of `face`. The face must have exactly
// one outer loop. Circular segments are resampled at `deflection` and the
// sample count is capped at 1,000,000 per ring before any allocation.
inline Result<std::vector<TrimmedBoundary>> trimmed_face_boundaries(
    const Brep& brep, FaceId face, const detail::FaceUvProjection& projection,
    double deflection) {
  if (!std::isfinite(deflection) || !(deflection > 0.0)) {
    return Error{Errc::tolerance_violation,
                 "trimmed_face_boundaries: deflection must be finite and positive", face.index(),
                 0u};
  }

  Result<std::vector<FaceUvLoop>> loops =
      face_loops_uv(brep, face, projection.tolerance);
  if (!loops) return loops.error();

  std::size_t outer_count = 0u;
  for (const FaceUvLoop& loop : loops.value()) {
    if (loop.outer) ++outer_count;
  }
  if (outer_count != 1u) {
    return Error{Errc::invalid_topology,
                 "trimmed_face_boundaries: face must have exactly one outer loop", face.index(),
                 static_cast<std::uint32_t>(outer_count)};
  }

  std::vector<TrimmedBoundary> result;
  result.reserve(loops.value().size());

  const auto append_ring = [&](const FaceUvLoop& loop) -> Result<void> {
    const std::vector<detail::TrimmedRingRun> runs = detail::trimmed_ring_runs(loop.arcs);
    const Result<std::size_t> planned =
        detail::trimmed_plan_ring(runs, deflection, face);
    if (!planned) return planned.error();
    Result<TrimmedBoundary> ring = detail::trimmed_build_ring(loop, projection, deflection, runs,
                                                              planned.value(), loop.outer, face);
    if (!ring) return ring.error();
    result.push_back(std::move(ring.value()));
    return {};
  };

  // The outer ring is reported first, matching face_loops_uv.
  for (const FaceUvLoop& loop : loops.value()) {
    if (!loop.outer) continue;
    const Result<void> appended = append_ring(loop);
    if (!appended) return appended.error();
  }
  for (const FaceUvLoop& loop : loops.value()) {
    if (loop.outer) continue;
    const Result<void> appended = append_ring(loop);
    if (!appended) return appended.error();
  }
  return result;
}

// Triangulates the trimmed uv domain described by `boundaries`. The outer ring
// comes first in the output point pool. Duplicate, non-finite, or structurally
// invalid rings are rejected; a failed ear clipping is reported as
// Errc::degenerate.
inline Result<TrimmedUvTriangulation> triangulate_trimmed_domain(
    std::span<const TrimmedBoundary> boundaries) {
  std::size_t outer_index = boundaries.size();
  std::size_t outer_count = 0u;
  for (std::size_t i = 0; i < boundaries.size(); ++i) {
    if (boundaries[i].outer) {
      ++outer_count;
      outer_index = i;
    }
  }
  if (outer_count != 1u) {
    return Error{Errc::invalid_topology,
                 "triangulate_trimmed_domain: exactly one outer ring is required", 0u,
                 static_cast<std::uint32_t>(outer_count)};
  }

  std::size_t total = 0u;
  for (const TrimmedBoundary& ring : boundaries) {
    if (ring.points.size() < 3u) {
      return Error{Errc::degenerate,
                   "triangulate_trimmed_domain: ring has fewer than three points", 0u, 0u};
    }
    if (ring.arcs.size() != ring.points.size()) {
      return Error{Errc::invalid_topology,
                   "triangulate_trimmed_domain: ring buffers are inconsistent", 0u, 0u};
    }
    for (const Point2d& point : ring.points) {
      if (!detail::trimmed_point_finite(point)) {
        return Error{Errc::non_finite,
                     "triangulate_trimmed_domain: ring point is not finite", 0u, 0u};
      }
    }
    total += ring.points.size();
  }
  if (total > static_cast<std::size_t>(std::numeric_limits<std::uint32_t>::max())) {
    return Error{Errc::numeric_failure,
                 "triangulate_trimmed_domain: too many ring points", 0u, 0u};
  }
  for (std::size_t i = 0; i < boundaries.size(); ++i) {
    for (std::size_t j = i + 1u; j < boundaries.size(); ++j) {
      const TrimmedBoundary& a = boundaries[i];
      const TrimmedBoundary& b = boundaries[j];
      if (a.points.size() == b.points.size() && a.points.front()[0] == b.points.front()[0] &&
          a.points.front()[1] == b.points.front()[1]) {
        return Error{Errc::degenerate, "triangulate_trimmed_domain: duplicate ring", 0u, 0u};
      }
    }
  }

  std::vector<Point2d> ordered;
  ordered.reserve(total);
  Polygon2d polygon;
  polygon.ring_starts.push_back(0u);
  polygon.vertices.reserve(total);

  const auto append_ring = [&](const TrimmedBoundary& ring) {
    ordered.insert(ordered.end(), ring.points.begin(), ring.points.end());
    polygon.vertices.insert(polygon.vertices.end(), ring.points.begin(), ring.points.end());
    polygon.ring_starts.push_back(static_cast<std::uint32_t>(polygon.vertices.size()));
  };

  append_ring(boundaries[outer_index]);
  for (std::size_t i = 0; i < boundaries.size(); ++i) {
    if (i == outer_index) continue;
    append_ring(boundaries[i]);
  }

  // Translate by the outer ring's first vertex: keeps the exact predicates in
  // their well-conditioned range and is undone by returning `ordered` verbatim.
  const Point2d origin = polygon.vertices.front();
  for (Point2d& vertex : polygon.vertices) {
    vertex[0] -= origin[0];
    vertex[1] -= origin[1];
  }

  if (const Errc status = validate_polygon2(polygon); status != Errc::none) {
    return Error{status, "triangulate_trimmed_domain: polygon is invalid", 0u, 0u};
  }

  const Triangulation2 triangulation = triangulate_polygon2(polygon);
  if (!triangulation.ok) {
    return Error{Errc::degenerate, "triangulate_trimmed_domain: triangulation failed", 0u, 0u};
  }
  for (std::size_t i = 0; i < triangulation.triangles.size(); ++i) {
    const std::array<std::uint32_t, 3>& triangle = triangulation.triangles[i];
    if (triangle[0] >= polygon.vertices.size() || triangle[1] >= polygon.vertices.size() ||
        triangle[2] >= polygon.vertices.size()) {
      return Error{Errc::degenerate, "triangulate_trimmed_domain: triangle index is out of range",
                   static_cast<std::uint32_t>(i), 0u};
    }
    const Point2d& a = polygon.vertices[triangle[0]];
    const Point2d& b = polygon.vertices[triangle[1]];
    const Point2d& c = polygon.vertices[triangle[2]];
    const double twice = (b[0] - a[0]) * (c[1] - a[1]) - (b[1] - a[1]) * (c[0] - a[0]);
    if (!(twice > 0.0) || !std::isfinite(twice)) {
      return Error{Errc::degenerate,
                   "triangulate_trimmed_domain: triangle is not positively oriented",
                   static_cast<std::uint32_t>(i), 0u};
    }
  }

  TrimmedUvTriangulation out;
  out.points = std::move(ordered);
  out.triangles = triangulation.triangles;
  return out;
}

}  // namespace v1
}  // namespace geom
