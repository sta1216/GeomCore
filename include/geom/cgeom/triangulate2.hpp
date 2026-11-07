#pragma once

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <span>
#include <utility>
#include <vector>

#include <geom/cgeom/polygon2.hpp>
#include <geom/core/error.hpp>
#include <geom/numeric/predicates2.hpp>
#include <geom/primitives/point.hpp>

namespace geom {
inline namespace v1 {

// Triangulation of a polygon with holes. Every triangle is counter-clockwise and
// stores indices into the original `Polygon2d::vertices` pool (not into an
// internal deduplicated pool). Triangles are emitted in the ear-clipping order;
// they are not sorted, so callers that need a canonical order must sort them.
struct Triangulation2 {
  std::vector<std::array<std::uint32_t, 3>> triangles;
  bool ok = false;
};

namespace detail {

// A ring vertex plus the index of the same vertex in the input Polygon2d. The
// original index survives deduplication and collinear-vertex removal.
struct TriVertex {
  Point2d p;
  std::uint32_t original;
};

using TriRing = std::vector<TriVertex>;

// Internal self-check tolerance for the covered-area invariant. The public
// contract checks 1e-12 relative in tests; this wider bound absorbs the
// floating-point summation error of large triangulations.
constexpr double kTriangulateAreaRelTol = 1e-9;

inline bool same_xy(const Point2d& a, const Point2d& b) noexcept {
  return a[0] == b[0] && a[1] == b[1];
}

inline bool lex_less_xy(const Point2d& a, const Point2d& b) noexcept {
  return a[0] < b[0] || (a[0] == b[0] && a[1] < b[1]);
}

// Drops consecutive duplicate coordinates, including a duplicate across the
// closing edge (e.g. an explicitly repeated first/last vertex). The surviving
// vertex keeps its original pool index.
inline TriRing dedup_ring(std::span<const Point2d> ring, std::uint32_t first) {
  TriRing out;
  out.reserve(ring.size());
  for (std::size_t i = 0; i < ring.size(); ++i) {
    TriVertex v{ring[i], first + static_cast<std::uint32_t>(i)};
    if (out.empty() || !same_xy(out.back().p, v.p)) out.push_back(v);
  }
  if (out.size() > 1 && same_xy(out.front().p, out.back().p)) out.pop_back();
  return out;
}

// Removes every exactly-collinear vertex, preserving the relative order of the
// remaining vertices (and therefore their original indices). A ring that
// collapses below three vertices is reported as degenerate.
inline bool remove_collinear_vertices(TriRing& ring) {
  bool changed = true;
  while (changed && ring.size() > 3) {
    changed = false;
    const std::size_t n = ring.size();
    for (std::size_t i = 0; i < n; ++i) {
      const TriVertex& prev = ring[(i + n - 1) % n];
      const TriVertex& cur = ring[i];
      const TriVertex& next = ring[(i + 1) % n];
      const bool duplicate = same_xy(prev.p, cur.p) || same_xy(cur.p, next.p);
      const bool collinear =
          sign_orient2d(prev.p[0], prev.p[1], cur.p[0], cur.p[1], next.p[0], next.p[1]) == 0;
      if (!duplicate && !collinear) continue;
      ring.erase(ring.begin() + static_cast<std::ptrdiff_t>(i));
      changed = true;
      break;
    }
  }
  return ring.size() >= 3;
}

// Orientation of a work ring, decided by the shared exact predicate through
// ring_orientation2. Returns 0 for a degenerate (zero-area) ring.
inline int work_orientation(const TriRing& ring) {
  std::vector<Point2d> pts;
  pts.reserve(ring.size());
  for (const TriVertex& v : ring) pts.push_back(v.p);
  return ring_orientation2(pts);
}

inline void normalize_orientation(TriRing& ring, int wanted) {
  if (work_orientation(ring) != wanted) std::reverse(ring.begin(), ring.end());
}

// True iff the closed triangle abc (assumed counter-clockwise) contains p; the
// boundary counts as contained. Exact, because all three half-plane tests use
// sign_orient2d.
inline bool point_in_closed_triangle(const Point2d& a, const Point2d& b, const Point2d& c,
                                     const Point2d& p) noexcept {
  return sign_orient2d(a[0], a[1], b[0], b[1], p[0], p[1]) >= 0 &&
         sign_orient2d(b[0], b[1], c[0], c[1], p[0], p[1]) >= 0 &&
         sign_orient2d(c[0], c[1], a[0], a[1], p[0], p[1]) >= 0;
}

// Proper crossing: both endpoints of each segment lie strictly on opposite
// sides of the other segment. Endpoint touches and collinear overlaps are not
// "proper" and are handled separately by the bridge visibility test.
inline bool properly_cross(const Point2d& a, const Point2d& b, const Point2d& c,
                           const Point2d& d) noexcept {
  const int o1 = sign_orient2d(a[0], a[1], b[0], b[1], c[0], c[1]);
  const int o2 = sign_orient2d(a[0], a[1], b[0], b[1], d[0], d[1]);
  if (o1 == 0 || o2 == 0 || o1 == o2) return false;
  const int o3 = sign_orient2d(c[0], c[1], d[0], d[1], a[0], a[1]);
  const int o4 = sign_orient2d(c[0], c[1], d[0], d[1], b[0], b[1]);
  return o3 != 0 && o4 != 0 && o3 != o4;
}

// True iff the open diagonal (a, c) properly crosses a ring edge that does not
// share a vertex with the diagonal. After hole bridges are spliced in, a convex
// corner with an empty ear triangle can still have a diagonal that pierces a
// non-incident edge; clipping such an ear would tear the ring. The check is
// exact (properly_cross uses sign_orient2d only). Edges incident to a or c are
// exempt by vertex index, exactly like their shared endpoints in the diagonal.
inline bool diagonal_crosses_ring(const TriRing& ring, std::size_t ia,
                                  std::size_t ic) noexcept {
  const std::size_t n = ring.size();
  const Point2d& a = ring[ia].p;
  const Point2d& c = ring[ic].p;
  for (std::size_t j = 0; j < n; ++j) {
    const std::size_t jn = (j + 1) % n;
    if (j == ia || j == ic || jn == ia || jn == ic) continue;
    if (properly_cross(a, c, ring[j].p, ring[jn].p)) return true;
  }
  return false;
}

// p lies on the open segment ab (strictly between a and b). The collinearity
// test is exact; the bounding-box test only accepts an already-collinear point.
inline bool on_open_segment(const Point2d& a, const Point2d& b, const Point2d& p) noexcept {
  if (same_xy(p, a) || same_xy(p, b)) return false;
  if (sign_orient2d(a[0], a[1], b[0], b[1], p[0], p[1]) != 0) return false;
  return p[0] >= std::min(a[0], b[0]) && p[0] <= std::max(a[0], b[0]) &&
         p[1] >= std::min(a[1], b[1]) && p[1] <= std::max(a[1], b[1]);
}

// A bridge (b, m) is visible when no ring edge properly crosses it and no ring
// vertex lies in its open interior. Edges incident to b or m are exempt, as in
// the classic bridge construction; duplicate copies of a bridge endpoint are
// recognised by coordinate equality.
inline bool bridge_is_visible(const Point2d& b, const Point2d& m, const TriRing& outer,
                              const std::vector<TriRing>& rings) {
  const auto edge_blocks = [&b, &m](const Point2d& u, const Point2d& v) {
    if (same_xy(u, b) || same_xy(v, b) || same_xy(u, m) || same_xy(v, m)) return false;
    if (properly_cross(b, m, u, v)) return true;
    return on_open_segment(b, m, u) || on_open_segment(b, m, v);
  };

  const std::size_t n = outer.size();
  for (std::size_t i = 0; i < n; ++i) {
    if (edge_blocks(outer[i].p, outer[(i + 1) % n].p)) return false;
  }
  for (const TriRing& ring : rings) {
    const std::size_t rn = ring.size();
    for (std::size_t i = 0; i < rn; ++i) {
      if (edge_blocks(ring[i].p, ring[(i + 1) % rn].p)) return false;
    }
  }
  return true;
}

struct RayHit {
  bool hit = false;
  std::size_t edge = 0;
  double x = 0.0;
};

// Nearest intersection of the +x ray starting at m with the current outer ring.
// The half-open y rule ((a.y > m.y) != (b.y > m.y)) counts a vertex exactly once,
// so a ray through a vertex does not double-count the two incident edges. On an
// exact x tie the first edge in ring order wins, which keeps the choice
// deterministic.
inline RayHit nearest_ray_hit(const TriRing& outer, const Point2d& m) noexcept {
  RayHit best;
  const std::size_t n = outer.size();
  for (std::size_t i = 0; i < n; ++i) {
    const Point2d& a = outer[i].p;
    const Point2d& b = outer[(i + 1) % n].p;
    if ((a[1] > m[1]) == (b[1] > m[1])) continue;
    const double t = (m[1] - a[1]) / (b[1] - a[1]);
    const double x = a[0] + t * (b[0] - a[0]);
    if (!(x > m[0])) continue;
    if (!best.hit || x < best.x) {
      best.hit = true;
      best.edge = i;
      best.x = x;
    }
  }
  return best;
}

// Splices the hole into the outer ring through the bridge (outer[b_index], hole[m_index]).
// Both endpoints appear twice in the merged ring, once on each side of the
// bridge. Hole rings are clockwise, so walking the hole from m keeps the polygon
// interior on the left, matching the outer ring's counter-clockwise winding.
inline void splice_hole(TriRing& outer, std::size_t b_index, const TriRing& hole,
                        std::size_t m_index) {
  TriRing merged;
  merged.reserve(outer.size() + hole.size() + 2);
  merged.insert(merged.end(), outer.begin(),
                outer.begin() + static_cast<std::ptrdiff_t>(b_index) + 1);
  const std::size_t h = hole.size();
  for (std::size_t k = 0; k < h; ++k) merged.push_back(hole[(m_index + k) % h]);
  merged.push_back(hole[m_index]);
  merged.insert(merged.end(), outer.begin() + static_cast<std::ptrdiff_t>(b_index), outer.end());
  outer.swap(merged);
}

// True when the direction from ring[i] to target lies inside the polygon's
// interior wedge at that vertex. After one or more holes are bridged the same
// coordinate can occur several times; every copy has its own pair of incident
// edges and therefore its own wedge.
inline bool interior_cone_contains(const TriRing& ring, std::size_t i,
                                   const Point2d& target) noexcept {
  const std::size_t n = ring.size();
  const std::size_t ip = (i + n - 1) % n;
  const std::size_t in = (i + 1) % n;
  const Point2d& p = ring[ip].p;
  const Point2d& c = ring[i].p;
  const Point2d& q = ring[in].p;
  const int turn = sign_orient2d(p[0], p[1], c[0], c[1], q[0], q[1]);
  const int s1 = sign_orient2d(p[0], p[1], c[0], c[1], target[0], target[1]);
  const int s2 = sign_orient2d(c[0], c[1], q[0], q[1], target[0], target[1]);
  if (turn > 0) return s1 >= 0 && s2 >= 0;
  if (turn < 0) return s1 >= 0 || s2 >= 0;
  return s1 >= 0;
}

// Bridges one clockwise hole into the counter-clockwise outer ring.
//
// Candidate order (deterministic):
//   1. the hole's rightmost vertex (max x, then max y);
//   2. if neither endpoint of the hit outer edge is visible, the same search is
//      repeated for the remaining hole vertices in ring order.
// For a candidate hole vertex, the two endpoints of the nearest +x ray hit edge
// are tried: first the endpoint with the larger x (the classic bridge choice),
// then the other one. An exact x tie between the endpoints is broken by
// lexicographic (x, then y) order, i.e. the lexicographically smaller visible
// partner wins that tie.
//
// An endpoint coordinate can occur more than once after earlier bridges. The
// copies are scanned in ring order and the first copy whose interior wedge
// contains the bridge direction is used; splicing into the wrong copy (for
// example the far side of an existing bridge) produces a locally inverted ring
// that ear clipping cannot resolve. If no copy passes the wedge test the first
// copy is used, preserving the pre-existing fallback.
inline bool bridge_hole(TriRing& outer, const TriRing& hole,
                        const std::vector<TriRing>& rings) {
  const std::size_t h = hole.size();
  if (h < 3 || outer.size() < 3) return false;

  std::size_t rightmost = 0;
  for (std::size_t i = 1; i < h; ++i) {
    const Point2d& p = hole[i].p;
    const Point2d& q = hole[rightmost].p;
    if (p[0] > q[0] || (p[0] == q[0] && p[1] > q[1])) rightmost = i;
  }

  for (std::size_t step = 0; step < h; ++step) {
    const std::size_t m_index = (rightmost + step) % h;
    const Point2d& m = hole[m_index].p;
    const RayHit hit = nearest_ray_hit(outer, m);
    if (!hit.hit) continue;

    const std::size_t n = outer.size();
    const std::size_t e0 = hit.edge;
    const std::size_t e1 = (hit.edge + 1) % n;
    std::size_t first = e0;
    std::size_t second = e1;
    if (outer[e0].p[0] < outer[e1].p[0] ||
        (outer[e0].p[0] == outer[e1].p[0] && lex_less_xy(outer[e1].p, outer[e0].p))) {
      first = e1;
      second = e0;
    }

    const auto splice_at = [&](const Point2d& coord) -> bool {
      std::size_t chosen = n;
      std::size_t fallback = n;
      for (std::size_t i = 0; i < n; ++i) {
        if (!same_xy(outer[i].p, coord)) continue;
        if (fallback == n) fallback = i;
        if (interior_cone_contains(outer, i, m)) {
          chosen = i;
          break;
        }
      }
      if (chosen == n) chosen = fallback;
      if (chosen == n) return false;
      splice_hole(outer, chosen, hole, m_index);
      return true;
    };

    const Point2d endpoints[2] = {outer[first].p, outer[second].p};
    const std::size_t endpoint_count = same_xy(endpoints[0], endpoints[1]) ? 1 : 2;
    for (std::size_t e = 0; e < endpoint_count; ++e) {
      const Point2d& coord = endpoints[e];
      if (!bridge_is_visible(coord, m, outer, rings)) continue;
      if (splice_at(coord)) return true;
    }

    // Fallback: when several holes are present, both endpoints of the ray-hit
    // edge can be blocked because the bridge would graze another hole. Any
    // visible outer vertex is then a valid partner; take the lexicographically
    // smallest visible coordinate so the choice stays deterministic.
    std::size_t best_coord = n;
    for (std::size_t i = 0; i < n; ++i) {
      const Point2d& coord = outer[i].p;
      bool seen = false;
      for (std::size_t j = 0; j < i && !seen; ++j) seen = same_xy(outer[j].p, coord);
      if (seen) continue;
      if (!bridge_is_visible(coord, m, outer, rings)) continue;
      if (best_coord == n || lex_less_xy(coord, outer[best_coord].p)) best_coord = i;
    }
    if (best_coord != n && splice_at(outer[best_coord].p)) return true;
  }
  return false;
}

// Ear clipping on the merged counter-clockwise ring. The lowest valid ear index
// is always removed first, so the output is deterministic. A convex corner is a
// valid ear when no other ring vertex lies in or on its closed triangle and its
// diagonal (prev, next) does not properly cross a non-incident ring edge;
// duplicate bridge vertices are skipped by coordinate equality. The diagonal
// test makes ear validation O(n) per candidate, so pathological rings are
// O(n^3) overall; M1b has no triangulation performance gate.
inline bool ear_clip(TriRing ring, std::vector<std::array<std::uint32_t, 3>>& out) {
  while (ring.size() > 3) {
    const std::size_t n = ring.size();
    bool clipped = false;
    for (std::size_t i = 0; i < n; ++i) {
      const std::size_t ip = (i + n - 1) % n;
      const std::size_t in = (i + 1) % n;
      const Point2d& a = ring[ip].p;
      const Point2d& b = ring[i].p;
      const Point2d& c = ring[in].p;
      if (sign_orient2d(a[0], a[1], b[0], b[1], c[0], c[1]) <= 0) continue;

      bool blocked = false;
      for (std::size_t j = 0; j < n; ++j) {
        if (j == ip || j == i || j == in) continue;
        const Point2d& q = ring[j].p;
        if (same_xy(q, a) || same_xy(q, b) || same_xy(q, c)) continue;
        if (point_in_closed_triangle(a, b, c, q)) {
          blocked = true;
          break;
        }
      }
      if (blocked) continue;

      // The ear triangle can be free of vertices yet still be invalid: after
      // hole bridging, its diagonal may properly cross another edge of the
      // ring. Reject the candidate; if no ear survives the ring is not simple.
      if (diagonal_crosses_ring(ring, ip, in)) continue;

      out.push_back({ring[ip].original, ring[i].original, ring[in].original});
      ring.erase(ring.begin() + static_cast<std::ptrdiff_t>(i));
      clipped = true;
      break;
    }
    if (!clipped) return false;
  }

  if (ring.size() != 3) return false;
  if (sign_orient2d(ring[0].p[0], ring[0].p[1], ring[1].p[0], ring[1].p[1], ring[2].p[0],
                    ring[2].p[1]) <= 0) {
    return false;
  }
  out.push_back({ring[0].original, ring[1].original, ring[2].original});
  return true;
}

inline long double triangle_area2(const Point2d& a, const Point2d& b,
                                  const Point2d& c) noexcept {
  const long double ab = static_cast<long double>(a[0]) * static_cast<long double>(b[1]) -
                         static_cast<long double>(b[0]) * static_cast<long double>(a[1]);
  const long double bc = static_cast<long double>(b[0]) * static_cast<long double>(c[1]) -
                         static_cast<long double>(c[0]) * static_cast<long double>(b[1]);
  const long double ca = static_cast<long double>(c[0]) * static_cast<long double>(a[1]) -
                         static_cast<long double>(a[0]) * static_cast<long double>(c[1]);
  return std::fabs(ab + bc + ca) * 0.5L;
}

}  // namespace detail

// Triangulates a simple polygon with optional holes.
//
// Method: validate -> normalize ring orientation (outer CCW, holes CW) -> drop
// consecutive duplicate points and exactly-collinear points -> bridge every hole
// to the outer ring through a deterministic visible vertex -> ear clipping.
// Convexity and point-in-triangle decisions use the exact sign_orient2d, so no
// global epsilon is involved. Complexity is O(n^2) in the number of vertices
// (each ear scan is O(n) and each hole bridge is O(n^2) worst case).
//
// Preconditions (documented, not repaired):
//   - validate_polygon2(poly) == Errc::none;
//   - every ring is simple and rings do not intersect each other;
//   - every hole lies strictly inside the outer ring and holes are not nested;
//   - ring windings may be arbitrary; adjacent duplicate points are tolerated.
// If validation or any precondition-dependent step fails, the result is
// {{}, false}; the function never indexes outside `poly.vertices` and throws no
// exceptions of its own (standard allocation failures may still propagate).
//
// On success every emitted triangle is counter-clockwise, non-degenerate, and
// indexes `poly.vertices` (indices refer to the original pool even though
// duplicate/collinear vertices were dropped internally). The covered area equals
// polygon_area2(poly) within kTriangulateAreaRelTol relative error.
inline Triangulation2 triangulate_polygon2(const Polygon2d& poly) {
  if (validate_polygon2(poly) != Errc::none) return {{}, false};
  if (poly.ring_count() == 0u) return {{}, false};
  if (poly.vertices.size() > static_cast<std::size_t>(std::numeric_limits<std::uint32_t>::max())) {
    return {{}, false};
  }

  std::vector<detail::TriRing> rings;
  rings.reserve(poly.ring_count());
  for (std::uint32_t r = 0; r < poly.ring_count(); ++r) {
    detail::TriRing ring = detail::dedup_ring(poly.ring(r), poly.ring_starts[r]);
    if (ring.size() < 3) return {{}, false};
    if (!detail::remove_collinear_vertices(ring)) return {{}, false};
    if (ring.size() < 3) return {{}, false};
    if (detail::work_orientation(ring) == 0) return {{}, false};
    detail::normalize_orientation(ring, r == 0u ? 1 : -1);
    rings.push_back(std::move(ring));
  }

  detail::TriRing merged = rings.front();
  for (std::size_t h = 1; h < rings.size(); ++h) {
    if (!detail::bridge_hole(merged, rings[h], rings)) return {{}, false};
  }

  const double expected_area = polygon_area2(poly);
  if (!(expected_area > 0.0)) return {{}, false};

  std::vector<std::array<std::uint32_t, 3>> triangles;
  if (!detail::ear_clip(merged, triangles)) return {{}, false};

  long double covered_area = 0.0L;
  for (const std::array<std::uint32_t, 3>& t : triangles) {
    const Point2d& a = poly.vertices[t[0]];
    const Point2d& b = poly.vertices[t[1]];
    const Point2d& c = poly.vertices[t[2]];
    if (sign_orient2d(a[0], a[1], b[0], b[1], c[0], c[1]) <= 0) return {{}, false};
    covered_area += detail::triangle_area2(a, b, c);
  }

  const long double target = static_cast<long double>(expected_area);
  const long double diff = covered_area > target ? covered_area - target : target - covered_area;
  const long double tolerance =
      static_cast<long double>(detail::kTriangulateAreaRelTol) * target;
  if (diff > tolerance) return {{}, false};

  return {std::move(triangles), true};
}

}  // namespace v1
}  // namespace geom
