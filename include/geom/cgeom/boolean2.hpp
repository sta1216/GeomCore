#pragma once

#include <geom/cgeom/arrangement2.hpp>
#include <geom/cgeom/polygon2.hpp>
#include <geom/core/error.hpp>
#include <geom/core/result.hpp>
#include <geom/numeric/predicates2.hpp>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <span>
#include <utility>
#include <vector>

namespace geom {
inline namespace v1 {

// Planar boolean operations on polygons, computed on an arrangement.
//
// Every ring of every operand is inserted into one Arrangement2 with the operand
// index as the segment source, so the arrangement resolves each crossing,
// overlap and shared edge once, exactly. Each cell then gets one winding per
// operand and is accepted by the requested fill rule. The aggregate
// `Face::winding` is never used for that decision: it cannot distinguish "two
// operands overlap here" from "one ring wraps twice".
//
// The accepted cells are the result region. Its boundary is traced on the DCEL
// by walking `next()` and crossing to `twin()` wherever an edge has accepted
// cells on both sides; loops with positive signed area become outer rings and
// loops with negative area become holes of the smallest outer ring containing
// them. The emitted set is therefore a regularized region even when the input
// rings self-intersect.
enum class BooleanOp {
  union_,               // covered by at least one operand
  intersection,         // covered by every operand
  difference,           // covered by operands[0] and by no other operand
  symmetric_difference  // covered by an odd number of operands
};

// Per-operand fill rule applied to that operand's own winding number.
enum class FillRule {
  non_zero,  // winding != 0
  even_odd,  // winding is odd
  positive   // winding > 0
};

struct BooleanOptions {
  FillRule fill_rule = FillRule::non_zero;
  // When false, two non-adjacent edges of the same input ring that meet make the
  // call fail with Errc::self_intersecting. Adjacent edges always meet at their
  // shared vertex, so they are exempt.
  bool allow_self_intersections = true;
};

struct BooleanReport {
  std::uint32_t input_edges = 0;       // ring segments handed to the arrangement
  std::uint32_t output_polygons = 0;   // accepted connected components
  std::uint32_t merged_vertices = 0;   // exactly collinear vertices removed
  double max_residual = 0.0;           // arrangement intersection residual
};

struct BooleanResult {
  PolygonSet polygons;
  BooleanReport report;
};

namespace detail {

// Lexicographic order of (x, y); Point2d has no ordering of its own because
// geometry code must never compare points by accident.
inline bool boolean_point_less(const Point2d& a, const Point2d& b) noexcept {
  if (a[0] != b[0]) return a[0] < b[0];
  return a[1] < b[1];
}

inline bool boolean_accept(std::int32_t winding, FillRule rule) noexcept {
  if (rule == FillRule::even_odd) return (winding % 2) != 0;
  if (rule == FillRule::positive) return winding > 0;
  return winding != 0;  // non_zero
}

// Signed shoelace area of one traced loop, computed after the exact power-of-two
// scaling 2^-exponent. Each loop is scaled by its own exponent (see the caller),
// so a tiny loop keeps a representable area even beside a huge one; the scaling
// preserves both sign and |area| order, which is all the classification needs.
inline double boolean_ring_area(std::span<const Point2d> ring, int exponent) noexcept {
  const std::size_t n = ring.size();
  if (n < 3) return 0.0;
  long double accumulated = 0.0L;
  for (std::size_t i = 0; i < n; ++i) {
    const Point2d& a = ring[i];
    const Point2d& b = ring[(i + 1) % n];
    const double ax = std::ldexp(a[0], -exponent);
    const double ay = std::ldexp(a[1], -exponent);
    const double bx = std::ldexp(b[0], -exponent);
    const double by = std::ldexp(b[1], -exponent);
    accumulated += static_cast<long double>(ax) * static_cast<long double>(by) -
                   static_cast<long double>(bx) * static_cast<long double>(ay);
  }
  return static_cast<double>(accumulated * 0.5L);
}

// Winding-number point test with exact orientation predicates. `ring` is a
// simple ring (an outer boundary), so a nonzero winding means "inside".
inline bool boolean_point_in_ring(const Point2d& p, std::span<const Point2d> ring) noexcept {
  const std::size_t n = ring.size();
  if (n < 3) return false;
  int winding = 0;
  for (std::size_t i = 0; i < n; ++i) {
    const Point2d& a = ring[i];
    const Point2d& b = ring[(i + 1) % n];
    const int side = sign_orient2d(a[0], a[1], b[0], b[1], p[0], p[1]);
    if (a[1] <= p[1]) {
      if (b[1] > p[1] && side > 0) ++winding;
    } else if (b[1] <= p[1] && side < 0) {
      --winding;
    }
  }
  return winding != 0;
}

// Removes vertices whose neighbours are exactly collinear with them (including
// coincident points) and returns how many went away. Deletion can expose a new
// collinear turn, so the ring is rescanned until it is stable.
inline std::uint32_t boolean_merge_collinear(std::vector<Point2d>& ring) noexcept {
  std::uint32_t removed = 0;
  bool changed = true;
  while (changed && ring.size() > 3) {
    changed = false;
    const std::size_t n = ring.size();
    for (std::size_t i = 0; i < n; ++i) {
      const Point2d& prev = ring[(i + n - 1) % n];
      const Point2d& cur = ring[i];
      const Point2d& next = ring[(i + 1) % n];
      if (sign_orient2d(prev[0], prev[1], cur[0], cur[1], next[0], next[1]) != 0) continue;
      ring.erase(ring.begin() + static_cast<std::ptrdiff_t>(i));
      ++removed;
      changed = true;
      break;
    }
  }
  return removed;
}

// Rotates a ring so it starts at its lexicographically smallest vertex; the
// arrangement offers no canonical start, and output must be reproducible.
inline void boolean_rotate_to_min(std::vector<Point2d>& ring) noexcept {
  if (ring.empty()) return;
  std::size_t best = 0;
  for (std::size_t i = 1; i < ring.size(); ++i) {
    if (boolean_point_less(ring[i], ring[best])) best = i;
  }
  std::rotate(ring.begin(), ring.begin() + static_cast<std::ptrdiff_t>(best), ring.end());
}

struct BooleanLoop {
  std::vector<Point2d> points;
  // |signed area| = area_mantissa * 2^area_exponent, with
  // area_mantissa in [0.5, 1). Keeping the magnitude as a mantissa and an
  // exponent lets hole attachment compare areas across the whole accepted span.
  double area_mantissa = 0.0;
  int area_exponent = 0;
};

// True when |area_a| < |area_b|, both given as mantissa * 2^exponent.
inline bool boolean_area_less(double mantissa_a, int exponent_a, double mantissa_b,
                              int exponent_b) noexcept {
  if (exponent_a != exponent_b) return exponent_a < exponent_b;
  return std::fabs(mantissa_a) < std::fabs(mantissa_b);
}

struct BooleanComponent {
  BooleanLoop outer;
  std::vector<BooleanLoop> holes;
};

}  // namespace detail

// Boolean combination of polygon sets. See BooleanOp for the per-operation
// semantics; `difference` subtracts the union of operands[1..] from operands[0].
//
// Errors:
//   - Errc::degenerate for an empty operand list or a malformed Polygon2d
//     (which includes a ring with fewer than three vertices);
//   - Errc::non_finite for a non-finite coordinate;
//   - Errc::unsupported when the coordinate exponent span exceeds 1314 bits,
//     the same guard delaunay2 uses (the arrangement face areas are scaled by
//     one global power of two and cannot separate such magnitudes);
//   - Errc::self_intersecting when allow_self_intersections is false and two
//     non-adjacent edges of the same input ring meet. Overlaps between distinct
//     rings or distinct operands are not self-intersections.
//
// Scale limit inherited from the arrangement: face winding is seeded from face
// areas scaled by one global power of two, so a feature more than about 1e154
// times smaller than the largest coordinate can still be misclassified even
// though the exponent span is accepted; keep every feature within ~1e154 of the
// largest coordinate.
//
// The result may be empty; an empty PolygonSet is a legal answer.
inline Result<BooleanResult> boolean2(std::span<const Polygon2d> operands, BooleanOp op,
                                     const BooleanOptions& options = {}) {
  const std::size_t operand_count = operands.size();
  if (operand_count == 0) {
    return Error{Errc::degenerate, "polygon boolean needs at least one operand", 0u, 0u};
  }

  // Reject unsupported input in every build, so the guard stays testable.
  constexpr int kMaxSpan = 1314;
  int max_exp = 0;
  int min_exp = 0;
  bool have_exp = false;
  for (std::size_t i = 0; i < operand_count; ++i) {
    const Errc problem = validate_polygon2(operands[i]);
    if (problem != Errc::none) {
      return Error{problem, "invalid polygon boolean operand",
                   static_cast<std::uint32_t>(i), 0u};
    }
    for (const Point2d& q : operands[i].vertices) {
      for (const double c : {q[0], q[1]}) {
        if (c == 0.0) continue;  // non-finite was rejected by validation
        int e = 0;
        std::frexp(std::fabs(c), &e);
        if (!have_exp) {
          max_exp = e;
          min_exp = e;
          have_exp = true;
        } else {
          if (e > max_exp) max_exp = e;
          if (e < min_exp) min_exp = e;
        }
      }
    }
  }
  if (have_exp && max_exp - min_exp > kMaxSpan) {
    return Error{Errc::unsupported, "polygon boolean exponent span exceeds 1314", 0u, 0u};
  }

  // 1. Normalize every ring (outer counter-clockwise, holes clockwise) and emit
  //    its segments with the operand index as the source. The flat segment index
  //    is remembered so input_pairs() can be mapped back to a ring position.
  std::vector<ArrangementSegment> segments;
  std::vector<std::uint32_t> segment_ring;
  std::vector<std::uint32_t> segment_position;
  std::vector<std::uint32_t> ring_sizes;
  segments.reserve(operands[0].vertices.size());
  for (std::size_t i = 0; i < operand_count; ++i) {
    const Polygon2d& poly = operands[i];
    for (std::uint32_t r = 0; r < poly.ring_count(); ++r) {
      const std::span<const Point2d> ring = poly.ring(r);
      std::vector<Point2d> oriented(ring.begin(), ring.end());
      const int orientation = ring_orientation2(ring);
      const bool want_ccw = (r == 0u);
      if ((want_ccw && orientation < 0) || (!want_ccw && orientation > 0)) {
        std::reverse(oriented.begin(), oriented.end());
      }
      const std::uint32_t ring_id = static_cast<std::uint32_t>(ring_sizes.size());
      const std::uint32_t ring_size = static_cast<std::uint32_t>(oriented.size());
      ring_sizes.push_back(ring_size);
      for (std::uint32_t k = 0; k < ring_size; ++k) {
        ArrangementSegment segment;
        segment.a = oriented[k];
        segment.b = oriented[(k + 1u) % ring_size];
        segment.source = static_cast<std::uint32_t>(i);
        segments.push_back(segment);
        segment_ring.push_back(ring_id);
        segment_position.push_back(k);
      }
    }
  }

  // 2. One arrangement over every operand resolves all crossings exactly.
  const Arrangement2 arrangement = Arrangement2::build(segments);
  const std::span<const Arrangement2::HalfEdge> half_edges = arrangement.half_edges();
  const std::span<const Arrangement2::Face> faces = arrangement.faces();
  const std::span<const Point2d> vertices = arrangement.vertices();

  if (!options.allow_self_intersections) {
    for (const std::pair<std::uint32_t, std::uint32_t>& pair : arrangement.input_pairs()) {
      const std::size_t first = pair.first;
      const std::size_t second = pair.second;
      if (segment_ring[first] != segment_ring[second]) continue;
      const std::uint32_t low = std::min(segment_position[first], segment_position[second]);
      const std::uint32_t high = std::max(segment_position[first], segment_position[second]);
      const std::uint32_t ring_size = ring_sizes[segment_ring[first]];
      const bool adjacent = (high - low == 1u) || (low == 0u && high + 1u == ring_size);
      if (!adjacent) {
        return Error{Errc::self_intersecting, "polygon ring self-intersects",
                     static_cast<std::uint32_t>(first), static_cast<std::uint32_t>(second)};
      }
    }
  }

  // 3. Classify each cell with the per-source winding of every operand. Every
  //    face record is classified, including the `!bounded` record that a hole
  //    side of a multi-connected region gets: those records carry the region's
  //    winding, while the true outer face keeps winding 0, which no fill rule
  //    accepts. Filtering on `Face::bounded` here would drop the holes.
  const std::size_t face_total = faces.size();
  std::vector<char> kept(face_total, 0);
  for (std::size_t f = 0; f < face_total; ++f) {
    std::size_t accepted = 0;
    for (std::size_t i = 0; i < operand_count; ++i) {
      if (detail::boolean_accept(
              arrangement.winding(static_cast<std::uint32_t>(f), static_cast<std::uint32_t>(i)),
              options.fill_rule)) {
        ++accepted;
      }
    }
    bool keep = false;
    switch (op) {
      case BooleanOp::union_: keep = accepted != 0; break;
      case BooleanOp::intersection: keep = accepted == operand_count; break;
      case BooleanOp::difference:
        keep = accepted == 1 &&
               detail::boolean_accept(
                   arrangement.winding(static_cast<std::uint32_t>(f), 0u), options.fill_rule);
        break;
      case BooleanOp::symmetric_difference: keep = (accepted % 2u) != 0; break;
    }
    kept[f] = keep ? 1 : 0;
  }

  // 4. Trace the boundary of the kept region. An edge is on the boundary when
  //    exactly one of its two cells is kept, and it is walked in the direction
  //    that has the kept side on the left. Stepping over an edge with kept cells
  //    on both sides moves the walk into the neighbouring cell, so a component
  //    spanning several cells yields one loop.
  const std::size_t half_edge_total = half_edges.size();
  std::vector<char> visited(half_edge_total, 0);
  std::vector<std::vector<Point2d>> loops;
  for (std::size_t start = 0; start < half_edge_total; ++start) {
    if (visited[start] != 0) continue;
    if (kept[half_edges[start].face] == 0) continue;
    if (kept[half_edges[half_edges[start].twin].face] != 0) continue;
    std::vector<Point2d> loop;
    std::uint32_t h = static_cast<std::uint32_t>(start);
    for (std::size_t step = 0; step <= half_edge_total; ++step) {
      visited[h] = 1;
      loop.push_back(vertices[half_edges[h].origin]);
      std::uint32_t edge = half_edges[h].next;
      for (std::size_t hop = 0; hop <= half_edge_total; ++hop) {
        if (kept[half_edges[half_edges[edge].twin].face] == 0) break;
        edge = half_edges[half_edges[edge].twin].next;
      }
      h = edge;
      if (h == start) break;
    }
    if (loop.size() >= 3) loops.push_back(std::move(loop));
  }

  // 5. Classify and canonicalize the loops. Every loop is normalized by its own
  //    power of two: a global exponent would underflow the area of a small loop
  //    next to a huge one and silently drop it.
  std::uint32_t merged_vertices = 0;
  std::vector<detail::BooleanComponent> components;
  std::vector<detail::BooleanLoop> holes;
  for (const std::vector<Point2d>& loop : loops) {
    double max_magnitude = 0.0;
    for (const Point2d& vertex : loop) {
      const double x = std::fabs(vertex[0]);
      const double y = std::fabs(vertex[1]);
      if (x > max_magnitude) max_magnitude = x;
      if (y > max_magnitude) max_magnitude = y;
    }
    const int loop_exp =
        (max_magnitude > 0.0 && std::isfinite(max_magnitude)) ? std::ilogb(max_magnitude) : 0;
    const double signed_area = detail::boolean_ring_area(loop, loop_exp);
    if (signed_area == 0.0) continue;  // zero-area sliver, not part of any region
    detail::BooleanLoop ring;
    ring.points = loop;
    const std::uint32_t removed = detail::boolean_merge_collinear(ring.points);
    if (ring.points.size() < 3) continue;
    merged_vertices += removed;
    detail::boolean_rotate_to_min(ring.points);
    int mantissa_exp = 0;
    const double mantissa = std::frexp(signed_area, &mantissa_exp);
    ring.area_mantissa = mantissa;
    ring.area_exponent = mantissa_exp + 2 * loop_exp;
    if (mantissa > 0.0) {
      detail::BooleanComponent component;
      component.outer = std::move(ring);
      components.push_back(std::move(component));
    } else {
      holes.push_back(std::move(ring));
    }
  }

  const auto by_first_vertex = [](const detail::BooleanLoop& a, const detail::BooleanLoop& b) {
    return detail::boolean_point_less(a.points.front(), b.points.front());
  };
  std::sort(components.begin(), components.end(),
            [&by_first_vertex](const detail::BooleanComponent& a,
                               const detail::BooleanComponent& b) {
              return by_first_vertex(a.outer, b.outer);
            });

  // 6. Attach every hole to the smallest outer ring containing its smallest
  //    vertex; nested outers then keep their own holes.
  for (detail::BooleanLoop& hole : holes) {
    if (components.empty()) break;
    std::size_t fallback = 0;
    for (std::size_t c = 1; c < components.size(); ++c) {
      const detail::BooleanLoop& candidate = components[c].outer;
      const detail::BooleanLoop& current = components[fallback].outer;
      if (detail::boolean_area_less(candidate.area_mantissa, candidate.area_exponent,
                                    current.area_mantissa, current.area_exponent)) {
        fallback = c;
      }
    }
    std::size_t best = fallback;
    bool found = false;
    for (std::size_t c = 0; c < components.size(); ++c) {
      if (!detail::boolean_point_in_ring(hole.points.front(), components[c].outer.points)) continue;
      if (!found) {
        best = c;
        found = true;
        continue;
      }
      const detail::BooleanLoop& candidate = components[c].outer;
      const detail::BooleanLoop& current = components[best].outer;
      if (detail::boolean_area_less(candidate.area_mantissa, candidate.area_exponent,
                                    current.area_mantissa, current.area_exponent)) {
        best = c;
      }
    }
    components[best].holes.push_back(std::move(hole));
  }

  BooleanResult result;
  for (detail::BooleanComponent& component : components) {
    std::sort(component.holes.begin(), component.holes.end(), by_first_vertex);
    Polygon2d poly;
    poly.ring_starts.push_back(0u);
    poly.vertices.insert(poly.vertices.end(), component.outer.points.begin(),
                         component.outer.points.end());
    poly.ring_starts.push_back(static_cast<std::uint32_t>(poly.vertices.size()));
    for (const detail::BooleanLoop& hole : component.holes) {
      poly.vertices.insert(poly.vertices.end(), hole.points.begin(), hole.points.end());
      poly.ring_starts.push_back(static_cast<std::uint32_t>(poly.vertices.size()));
    }
    result.polygons.polygons.push_back(std::move(poly));
  }

  result.report.input_edges = static_cast<std::uint32_t>(segments.size());
  result.report.output_polygons = static_cast<std::uint32_t>(result.polygons.polygons.size());
  result.report.merged_vertices = merged_vertices;
  result.report.max_residual = arrangement.max_residual();
  return result;
}

}  // namespace v1
}  // namespace geom
