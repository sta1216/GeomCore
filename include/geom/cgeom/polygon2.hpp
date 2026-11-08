#pragma once

#include <cassert>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <span>
#include <vector>

#include <geom/core/error.hpp>
#include <geom/numeric/predicates2.hpp>
#include <geom/primitives/point.hpp>

namespace geom {
inline namespace v1 {

// A 2D polygon stored as one shared vertex pool plus a ring index array: ring r
// owns the vertices [ring_starts[r], ring_starts[r + 1]). The first ring is the
// outer boundary, the remaining rings are holes. Rings are not closed: the last
// vertex is implicitly connected back to the first one.
//
// `ring_starts` is always expected to be {0, c0, c0 + c1, ..., vertices.size()};
// validate_polygon2 checks that structure and the finiteness of every coordinate.
struct Polygon2d {
  std::vector<Point2d> vertices;
  std::vector<std::uint32_t> ring_starts;

  std::uint32_t ring_count() const noexcept {
    return ring_starts.empty() ? 0u : static_cast<std::uint32_t>(ring_starts.size() - 1);
  }

  // Precondition: i < ring_count(). The returned span aliases `vertices`, which
  // must outlive it.
  std::span<const Point2d> ring(std::uint32_t i) const noexcept {
    assert(i < ring_count());
    const std::uint32_t first = ring_starts[i];
    const std::uint32_t count = ring_starts[i + 1] - first;
    return std::span<const Point2d>(vertices.data() + first, static_cast<std::size_t>(count));
  }
};

struct PolygonSet {
  std::vector<Polygon2d> polygons;
};

// Structural validation of a Polygon2d. Coordinates are checked first, so a
// non-finite coordinate is reported as Errc::non_finite even when the ring table
// is also malformed. The structure is valid iff it describes at least one ring:
//   - ring_starts.size() >= 2,
//   - ring_starts.front() == 0 and ring_starts.back() == vertices.size(),
//   - ring_starts is strictly increasing,
//   - every ring holds at least 3 vertices.
// Any structural violation is reported as Errc::degenerate. This is a structural
// check only: duplicate vertices, collinear runs and self-intersections are not
// rejected.
inline Errc validate_polygon2(const Polygon2d& poly) noexcept {
  for (const Point2d& p : poly.vertices) {
    if (!(std::isfinite(p[0]) && std::isfinite(p[1]))) return Errc::non_finite;
  }

  const std::size_t n = poly.ring_starts.size();
  if (n < 2) return Errc::degenerate;
  if (poly.ring_starts.front() != 0u) return Errc::degenerate;
  if (static_cast<std::size_t>(poly.ring_starts.back()) != poly.vertices.size()) {
    return Errc::degenerate;
  }
  for (std::size_t i = 1; i < n; ++i) {
    if (poly.ring_starts[i] <= poly.ring_starts[i - 1]) return Errc::degenerate;
    if (poly.ring_starts[i] - poly.ring_starts[i - 1] < 3u) return Errc::degenerate;
  }
  return Errc::none;
}

// Exact orientation of a closed ring: 1 for counter-clockwise, -1 for clockwise,
// 0 when the orientation is undefined (fewer than 3 vertices, every vertex equal,
// or a degenerate turn at the extreme vertex).
//
// The decision is made at the lexicographically smallest vertex (x, then y),
// where the interior angle is convex for a simple ring, so the local turn through
// that vertex has the orientation of the whole ring. Vertices equal to it are
// skipped on both sides, so repeated points (including a duplicated closing
// vertex) do not change the answer. The turn itself is classified by
// sign_orient2d, never by a floating-point cross product; a simple ring is
// assumed, self-intersecting rings are undefined.
inline int ring_orientation2(std::span<const Point2d> ring) noexcept {
  const std::size_t n = ring.size();
  if (n < 3) return 0;

  std::size_t k = 0;
  for (std::size_t i = 1; i < n; ++i) {
    const Point2d& p = ring[i];
    const Point2d& q = ring[k];
    if (p[0] < q[0] || (p[0] == q[0] && p[1] < q[1])) k = i;
  }

  const double kx = ring[k][0];
  const double ky = ring[k][1];
  const auto equals_k = [kx, ky](const Point2d& p) noexcept {
    return p[0] == kx && p[1] == ky;
  };

  std::size_t prev = k;
  for (std::size_t step = 1; step < n; ++step) {
    const std::size_t j = (k + n - step) % n;
    if (!equals_k(ring[j])) {
      prev = j;
      break;
    }
  }
  std::size_t next = k;
  for (std::size_t step = 1; step < n; ++step) {
    const std::size_t j = (k + step) % n;
    if (!equals_k(ring[j])) {
      next = j;
      break;
    }
  }

  // The same neighbour on both sides means no two distinct directions exist.
  if (prev == next) return 0;
  return sign_orient2d(ring[prev][0], ring[prev][1], kx, ky, ring[next][0], ring[next][1]);
}

// Shoelace area of a closed ring, signed by its winding. `long double` keeps the
// accumulation wider than the inputs, but this is still a floating-point
// approximation for reporting and tests; it is not a topological predicate.
inline double signed_area2(std::span<const Point2d> ring) noexcept {
  const std::size_t n = ring.size();
  if (n < 3) return 0.0;

  long double acc = 0.0L;
  for (std::size_t i = 0; i < n; ++i) {
    const Point2d& a = ring[i];
    const Point2d& b = ring[(i + 1) % n];
    acc += static_cast<long double>(a[0]) * static_cast<long double>(b[1]) -
           static_cast<long double>(b[0]) * static_cast<long double>(a[1]);
  }
  return static_cast<double>(acc * 0.5L);
}

// Area of the polygon: the absolute outer area minus the absolute area of every
// hole. Winding of each ring is irrelevant because absolute values are used.
// Precondition: validate_polygon2(poly) == Errc::none.
inline double polygon_area2(const Polygon2d& poly) noexcept {
  if (poly.ring_count() == 0u) return 0.0;

  double area = std::abs(signed_area2(poly.ring(0)));
  for (std::uint32_t r = 1; r < poly.ring_count(); ++r) {
    area -= std::abs(signed_area2(poly.ring(r)));
  }
  return area;
}

// Total area of every polygon in the set.
inline double polygon_set_area2(const PolygonSet& set) noexcept {
  long double total = 0.0L;
  for (const Polygon2d& poly : set.polygons) {
    total += static_cast<long double>(polygon_area2(poly));
  }
  return static_cast<double>(total);
}

}  // namespace v1
}  // namespace geom
