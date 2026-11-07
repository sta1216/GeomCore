#pragma once

#include <algorithm>
#include <array>
#include <cassert>
#include <cmath>
#include <cstdint>
#include <limits>
#include <span>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>

#include <geom/cgeom/delaunay2.hpp>
#include <geom/cgeom/triangulate2.hpp>
#include <geom/core/error.hpp>
#include <geom/core/result.hpp>
#include <geom/numeric/predicates2.hpp>
#include <geom/primitives/point.hpp>

namespace geom {
inline namespace v1 {

// A single constraint edge, given as two indices into the input point array.
struct CdtConstraint {
  std::uint32_t a = 0;
  std::uint32_t b = 0;
};

// Constrained Delaunay triangulation of a planar point set.
//
// The output is a triangulation whose vertex set is the input points, whose
// triangles are counter-clockwise and index into the input, sorted by
// (a, b, c), and which contains every constraint edge. Away from the
// constraints the triangulation is Delaunay: for every unconstrained interior
// edge, the opposite vertex of the neighbouring triangle is not strictly
// inside the circumcircle.
struct CdtResult {
  std::vector<std::array<std::uint32_t, 3>> triangles;
};

namespace detail {

inline bool cdt_same_xy(const Point2d& a, const Point2d& b) noexcept {
  return a[0] == b[0] && a[1] == b[1];
}

inline std::uint64_t cdt_directed(std::uint32_t u, std::uint32_t v) noexcept {
  return (static_cast<std::uint64_t>(u) << 32) | v;
}

// Canonical key for the undirected edge {u, v}.
inline std::uint64_t cdt_undirected(std::uint32_t u, std::uint32_t v) noexcept {
  return u < v ? cdt_directed(u, v) : cdt_directed(v, u);
}

inline bool cdt_on_segment(const Point2d& a, const Point2d& b, const Point2d& p) noexcept {
  if (sign_orient2d(a[0], a[1], b[0], b[1], p[0], p[1]) != 0) return false;
  return std::min(a[0], b[0]) <= p[0] && p[0] <= std::max(a[0], b[0]) &&
         std::min(a[1], b[1]) <= p[1] && p[1] <= std::max(a[1], b[1]);
}

// True when the two (undirected, distinct) constraint segments are not allowed
// to coexist: they cross properly, touch without sharing an endpoint, overlap
// collinearly, or are collinear and leave a shared endpoint in the same
// direction (one contains the other's start). Two segments that meet only at a
// shared endpoint are allowed, which is exactly the stated contract.
inline bool cdt_constraints_conflict(const Point2d& a, const Point2d& b, const Point2d& c,
                                     const Point2d& d) noexcept {
  const bool share_ac = cdt_same_xy(a, c);
  const bool share_ad = cdt_same_xy(a, d);
  const bool share_bc = cdt_same_xy(b, c);
  const bool share_bd = cdt_same_xy(b, d);

  if (share_ac || share_ad || share_bc || share_bd) {
    const Point2d* p = nullptr;
    const Point2d* q = nullptr;
    const Point2d* r = nullptr;
    if (share_ac) {
      p = &a;
      q = &b;
      r = &d;
    } else if (share_ad) {
      p = &a;
      q = &b;
      r = &c;
    } else if (share_bc) {
      p = &b;
      q = &a;
      r = &d;
    } else {
      p = &b;
      q = &a;
      r = &c;
    }
    // Non-collinear segments meet only at the shared endpoint.
    if (sign_orient2d((*p)[0], (*p)[1], (*q)[0], (*q)[1], (*r)[0], (*r)[1]) != 0) {
      return false;
    }
    // Collinear: overlap beyond the shared endpoint iff q and r leave p in the
    // same direction. Compare on an axis where q actually moves away from p.
    if ((*q)[0] != (*p)[0]) return ((*q)[0] > (*p)[0]) == ((*r)[0] > (*p)[0]);
    return ((*q)[1] > (*p)[1]) == ((*r)[1] > (*p)[1]);
  }

  const int o1 = sign_orient2d(a[0], a[1], b[0], b[1], c[0], c[1]);
  const int o2 = sign_orient2d(a[0], a[1], b[0], b[1], d[0], d[1]);
  const int o3 = sign_orient2d(c[0], c[1], d[0], d[1], a[0], a[1]);
  const int o4 = sign_orient2d(c[0], c[1], d[0], d[1], b[0], b[1]);

  if (o1 == 0 && o2 == 0) {
    // All four points are collinear; reject any common point.
    double lo1 = 0.0, hi1 = 0.0, lo2 = 0.0, hi2 = 0.0;
    if (a[0] != b[0]) {
      lo1 = std::min(a[0], b[0]);
      hi1 = std::max(a[0], b[0]);
      lo2 = std::min(c[0], d[0]);
      hi2 = std::max(c[0], d[0]);
    } else {
      lo1 = std::min(a[1], b[1]);
      hi1 = std::max(a[1], b[1]);
      lo2 = std::min(c[1], d[1]);
      hi2 = std::max(c[1], d[1]);
    }
    return !(hi1 < lo2 || hi2 < lo1);
  }

  // Proper crossing: each segment strictly separates the other's endpoints.
  if (o1 != 0 && o2 != 0 && o3 != 0 && o4 != 0 && ((o1 > 0) != (o2 > 0)) &&
      ((o3 > 0) != (o4 > 0))) {
    return true;
  }

  // Endpoint-on-segment touches (a T-junction), without a shared index.
  if (o1 == 0 && cdt_on_segment(a, b, c)) return true;
  if (o2 == 0 && cdt_on_segment(a, b, d)) return true;
  if (o3 == 0 && cdt_on_segment(c, d, a)) return true;
  if (o4 == 0 && cdt_on_segment(c, d, b)) return true;
  return false;
}

struct CdtTri {
  std::array<std::uint32_t, 3> v;
  bool alive = true;
};

inline void cdt_make_ccw(std::array<std::uint32_t, 3>& v, std::span<const Point2d> p) {
  const int s = sign_orient2d(p[v[0]][0], p[v[0]][1], p[v[1]][0], p[v[1]][1], p[v[2]][0],
                              p[v[2]][1]);
  if (s < 0) std::swap(v[1], v[2]);
  assert(s != 0 && "cdt2: degenerate triangle");
}

// Directed edge -> owning live triangle.
inline std::unordered_map<std::uint64_t, std::uint32_t> cdt_owners(
    const std::vector<CdtTri>& tris) {
  std::unordered_map<std::uint64_t, std::uint32_t> owner;
  owner.reserve(3 * tris.size());
  for (std::uint32_t t = 0; t < tris.size(); ++t) {
    if (!tris[t].alive) continue;
    for (int e = 0; e < 3; ++e) {
      owner[cdt_directed(tris[t].v[e], tris[t].v[(e + 1) % 3])] = t;
    }
  }
  return owner;
}

// Recovers one constraint edge with a Sloan-style cavity walk: collect every
// triangle crossed by the segment, delete them, split the cavity boundary by
// the segment into two simple polygons, and ear-clip both. The call assumes the
// constraint does not already exist.
inline bool cdt_recover(std::vector<CdtTri>& tris, std::span<const Point2d> p,
                        std::uint32_t a, std::uint32_t b) {
  auto owner = cdt_owners(tris);
  if (owner.find(cdt_directed(a, b)) != owner.end() ||
      owner.find(cdt_directed(b, a)) != owner.end()) {
    return true;
  }

  // Any triangle incident to `a` whose interior wedge at `a` contains the
  // direction towards `b` starts the walk. The wedge is measured with the exact
  // orientation predicate.
  std::uint32_t start = 0xFFFFFFFFu;
  for (std::uint32_t t = 0; t < tris.size() && start == 0xFFFFFFFFu; ++t) {
    if (!tris[t].alive) continue;
    const auto& v = tris[t].v;
    for (int k = 0; k < 3; ++k) {
      if (v[k] != a) continue;
      const std::uint32_t x = v[(k + 1) % 3];
      const std::uint32_t y = v[(k + 2) % 3];
      const int s1 = sign_orient2d(p[a][0], p[a][1], p[x][0], p[x][1], p[b][0], p[b][1]);
      const int s2 = sign_orient2d(p[a][0], p[a][1], p[y][0], p[y][1], p[b][0], p[b][1]);
      if (s1 > 0 && s2 < 0) {
        start = t;
        break;
      }
    }
  }
  if (start == 0xFFFFFFFFu) return false;

  std::vector<std::uint32_t> cavity;
  std::uint32_t cur = start;
  bool reached = false;
  for (std::size_t guard = 0; guard <= tris.size(); ++guard) {
    cavity.push_back(cur);
    bool has_b = false;
    for (int k = 0; k < 3; ++k) {
      if (tris[cur].v[k] == b) has_b = true;
    }
    if (has_b) {
      reached = true;
      break;
    }
    // The segment leaves the current triangle through the unique edge that it
    // properly crosses. Edges incident to `a` are never crossed. Requiring
    // o3 > 0 and o4 < 0 selects the exit edge: `a` lies strictly on the
    // triangle's interior side of it and `b` strictly beyond it. Without the
    // side signs the just-entered edge (a outside, b inside) also passes the
    // opposition test, so the walk could step back and oscillate.
    std::uint32_t nxt = 0xFFFFFFFFu;
    for (int e = 0; e < 3; ++e) {
      const std::uint32_t u = tris[cur].v[e];
      const std::uint32_t w = tris[cur].v[(e + 1) % 3];
      if (u == a || w == a) continue;
      const int o1 = sign_orient2d(p[a][0], p[a][1], p[b][0], p[b][1], p[u][0], p[u][1]);
      const int o2 = sign_orient2d(p[a][0], p[a][1], p[b][0], p[b][1], p[w][0], p[w][1]);
      const int o3 = sign_orient2d(p[u][0], p[u][1], p[w][0], p[w][1], p[a][0], p[a][1]);
      const int o4 = sign_orient2d(p[u][0], p[u][1], p[w][0], p[w][1], p[b][0], p[b][1]);
      if (o1 == 0 || o2 == 0 || o3 == 0 || o4 == 0) continue;
      if (!((o1 > 0) != (o2 > 0)) || o3 <= 0 || o4 >= 0) continue;
      const auto it = owner.find(cdt_directed(w, u));
      if (it == owner.end()) return false;
      nxt = it->second;
      break;
    }
    if (nxt == 0xFFFFFFFFu) return false;
    cur = nxt;
  }
  if (!reached) return false;

  // Cavity boundary: directed edges of the cavity whose reverse is not itself
  // a cavity edge. Their orientation follows the cavity's CCW boundary.
  std::unordered_set<std::uint64_t> cavity_edges;
  cavity_edges.reserve(3 * cavity.size());
  for (const std::uint32_t t : cavity) {
    for (int e = 0; e < 3; ++e) {
      cavity_edges.insert(cdt_directed(tris[t].v[e], tris[t].v[(e + 1) % 3]));
    }
  }
  std::vector<std::pair<std::uint32_t, std::uint32_t>> boundary;
  for (const std::uint32_t t : cavity) {
    for (int e = 0; e < 3; ++e) {
      const std::uint32_t u = tris[t].v[e];
      const std::uint32_t w = tris[t].v[(e + 1) % 3];
      if (cavity_edges.find(cdt_directed(w, u)) == cavity_edges.end()) {
        boundary.emplace_back(u, w);
      }
    }
  }

  std::unordered_map<std::uint32_t, std::uint32_t> next;
  next.reserve(boundary.size());
  for (const auto& [u, w] : boundary) {
    if (!next.emplace(u, w).second) return false;  // boundary is not simple
  }

  // Walk the cycle once, then split it at b. The two halves plus the constraint
  // diagonal are exactly the two replacement polygons.
  std::vector<std::uint32_t> cycle;
  std::uint32_t v = a;
  do {
    cycle.push_back(v);
    const auto it = next.find(v);
    if (it == next.end()) return false;
    v = it->second;
    if (cycle.size() > boundary.size() + 1) return false;
  } while (v != a);

  std::size_t pos_b = cycle.size();
  for (std::size_t i = 0; i < cycle.size(); ++i) {
    if (cycle[i] == b) pos_b = i;
  }
  if (pos_b == cycle.size()) return false;

  std::vector<std::uint32_t> poly1(cycle.begin(), cycle.begin() + static_cast<std::ptrdiff_t>(pos_b) + 1);
  std::vector<std::uint32_t> poly2(cycle.begin() + static_cast<std::ptrdiff_t>(pos_b), cycle.end());
  poly2.push_back(a);

  auto make_ring = [&](std::vector<std::uint32_t>& poly) -> bool {
    std::vector<Point2d> pts;
    pts.reserve(poly.size());
    for (const std::uint32_t i : poly) pts.push_back(p[i]);
    const int orientation = ring_orientation2(pts);
    if (orientation == 0) return false;
    if (orientation < 0) std::reverse(poly.begin(), poly.end());
    return poly.size() >= 3;
  };
  if (!make_ring(poly1) || !make_ring(poly2)) return false;

  std::vector<std::array<std::uint32_t, 3>> fresh;
  for (std::vector<std::uint32_t>* poly : {&poly1, &poly2}) {
    TriRing ring;
    ring.reserve(poly->size());
    for (const std::uint32_t i : *poly) ring.push_back(TriVertex{p[i], i});
    if (!ear_clip(std::move(ring), fresh)) return false;
  }

  for (const std::uint32_t t : cavity) tris[t].alive = false;
  for (auto& t : fresh) {
    CdtTri added{t, true};
    cdt_make_ccw(added.v, p);
    tris.push_back(added);
  }
  return true;
}

// Restores the Delaunay property on every unconstrained interior edge by legal
// flips. Constraint edges are never flipped; a flip is only attempted across a
// strictly convex quadrilateral. The flip count is capped by four times the
// triangle count to turn a numerical pathology into an error rather than a
// livelock.
inline bool cdt_restore_delaunay(std::vector<CdtTri>& tris, std::span<const Point2d> p,
                                 const std::unordered_set<std::uint64_t>& constrained) {
  const std::size_t cap = 4 * tris.size();
  std::size_t flips = 0;
  bool flipped = true;
  while (flipped) {
    flipped = false;
    const auto owner = cdt_owners(tris);
    for (std::uint32_t t = 0; t < tris.size() && !flipped; ++t) {
      if (!tris[t].alive) continue;
      const std::array<std::uint32_t, 3> T = tris[t].v;
      for (int e = 0; e < 3 && !flipped; ++e) {
        const std::uint32_t a = T[e];
        const std::uint32_t b = T[(e + 1) % 3];
        const std::uint32_t c = T[(e + 2) % 3];
        if (a > b) continue;  // each interior edge once, from its CCW side
        if (constrained.find(cdt_undirected(a, b)) != constrained.end()) continue;
        const auto it = owner.find(cdt_directed(b, a));
        if (it == owner.end()) continue;  // hull edge
        const std::array<std::uint32_t, 3> U = tris[it->second].v;
        std::uint32_t d = 0xFFFFFFFFu;
        for (int k = 0; k < 3; ++k) {
          if (U[k] != a && U[k] != b) d = U[k];
        }
        if (d == 0xFFFFFFFFu) continue;
        const int ua = sign_orient2d(p[c][0], p[c][1], p[d][0], p[d][1], p[a][0], p[a][1]);
        const int ub = sign_orient2d(p[c][0], p[c][1], p[d][0], p[d][1], p[b][0], p[b][1]);
        if (ua == 0 || ub == 0 || ((ua > 0) == (ub > 0))) continue;  // not convex
        if (sign_incircle(p[a][0], p[a][1], p[b][0], p[b][1], p[c][0], p[c][1], p[d][0],
                          p[d][1]) <= 0) {
          continue;  // legal edge
        }
        std::array<std::uint32_t, 3> first{c, d, b};
        std::array<std::uint32_t, 3> second{c, a, d};
        cdt_make_ccw(first, p);
        cdt_make_ccw(second, p);
        tris[t].v = first;
        tris[it->second].v = second;
        ++flips;
        if (flips > cap) return false;
        flipped = true;
      }
    }
  }
  return true;
}

}  // namespace detail

// Constrained Delaunay triangulation.
//
// Contract:
//   - finite coordinates; input.size() <= 2^32 - 1 so indices fit uint32_t;
//   - the coordinate exponent span is at most 1314 bits, the largest span the
//     exact predicates can resolve. Wider spans return Errc::unsupported;
//   - duplicate points are skipped, keeping the first occurrence; a constraint
//     that names duplicates is mapped to the surviving representative;
//   - non-constraint points are in general position (no three collinear, no
//     cocircular quadruple), the same precondition as delaunay2. Constraint
//     endpoints may be collinear with one another, but constraint segments may
//     only meet at a shared endpoint. Proper crossings, T-junction touches and
//     collinear overlaps between distinct constraints return
//     Errc::self_intersecting;
//   - duplicate constraints, including a reversed duplicate, are accepted and
//     deduplicated.
//
// Errors: Errc::non_finite, Errc::index_out_of_range, Errc::degenerate (fewer
// than three distinct points, all points collinear, or a constraint with equal
// endpoints), Errc::self_intersecting, Errc::unsupported (exponent span), and
// Errc::numeric_failure if the recovery/flip machinery cannot finish.
//
// Method: build the unconstrained delaunay2 mesh, then recover each missing
// constraint edge by deleting the triangles its segment crosses and
// retriangulating the resulting cavity (split along the constraint by ear
// clipping). Finally, restore the Delaunay property with Lawson flips that skip
// constraint edges. Every geometric decision uses the exact sign_orient2d /
// sign_incircle predicates.
inline Result<CdtResult> cdt2(std::span<const Point2d> points,
                              std::span<const CdtConstraint> constraints) {
  assert(points.size() <= 0xFFFFFFFFull);

  constexpr int kMaxSpan = 1314;
  int max_exp = 0;
  int min_exp = 0;
  bool have_exp = false;
  for (const Point2d& q : points) {
    for (const double c : {q[0], q[1]}) {
      if (!std::isfinite(c)) {
        return Error{Errc::non_finite, "cdt2: non-finite input coordinate", 0u, 0u};
      }
      if (c == 0.0) continue;
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
  if (points.size() < 3) {
    return Error{Errc::degenerate, "cdt2: needs at least three points", 0u, 0u};
  }
  if (have_exp && max_exp - min_exp > kMaxSpan) {
    return Error{Errc::unsupported, "cdt2: coordinate exponent span exceeds 1314", 0u, 0u};
  }

  // Map every input index to the first occurrence of its coordinates, matching
  // delaunay2's duplicate handling.
  std::vector<std::uint32_t> representative(points.size(), 0u);
  for (std::uint32_t i = 0; i < points.size(); ++i) {
    representative[i] = i;
    for (std::uint32_t j = 0; j < i; ++j) {
      if (detail::cdt_same_xy(points[i], points[j])) {
        representative[i] = representative[j];
        break;
      }
    }
  }

  std::vector<std::pair<std::uint32_t, std::uint32_t>> edges;
  std::unordered_set<std::uint64_t> seen;
  edges.reserve(constraints.size());
  seen.reserve(constraints.size());
  for (const CdtConstraint& constraint : constraints) {
    if (constraint.a >= points.size() || constraint.b >= points.size()) {
      return Error{Errc::index_out_of_range, "cdt2: constraint index out of range",
                   constraint.a, constraint.b};
    }
    const std::uint32_t a = representative[constraint.a];
    const std::uint32_t b = representative[constraint.b];
    if (a == b) {
      return Error{Errc::degenerate, "cdt2: constraint endpoints coincide", constraint.a,
                   constraint.b};
    }
    const std::uint64_t key = detail::cdt_undirected(a, b);
    if (seen.find(key) != seen.end()) continue;  // duplicate, including reversed
    for (const auto& [c, d] : edges) {
      if (detail::cdt_constraints_conflict(points[a], points[b], points[c], points[d])) {
        return Error{Errc::self_intersecting, "cdt2: constraints intersect", a, b};
      }
    }
    seen.insert(key);
    edges.emplace_back(a, b);
  }

  const std::vector<std::array<std::uint32_t, 3>> base = delaunay2(points);
  if (base.empty()) {
    return Error{Errc::degenerate, "cdt2: no triangulation (degenerate point set)", 0u, 0u};
  }

  std::vector<detail::CdtTri> tris;
  tris.reserve(base.size() + edges.size() * 2);
  for (const auto& t : base) tris.push_back(detail::CdtTri{t, true});

  std::unordered_set<std::uint64_t> constrained;
  constrained.reserve(edges.size());
  for (const auto& [a, b] : edges) {
    bool present = false;
    for (const detail::CdtTri& t : tris) {
      for (int e = 0; e < 3; ++e) {
        if ((t.v[e] == a && t.v[(e + 1) % 3] == b) ||
            (t.v[e] == b && t.v[(e + 1) % 3] == a)) {
          present = true;
        }
      }
    }
    if (!present) {
      if (!detail::cdt_recover(tris, points, a, b)) {
        return Error{Errc::numeric_failure, "cdt2: constraint recovery failed", a, b};
      }
    }
    constrained.insert(detail::cdt_undirected(a, b));
  }

  if (!detail::cdt_restore_delaunay(tris, points, constrained)) {
    return Error{Errc::numeric_failure, "cdt2: Delaunay restoration did not converge", 0u, 0u};
  }

  CdtResult result;
  result.triangles.reserve(tris.size());
  for (const detail::CdtTri& t : tris) {
    if (!t.alive) continue;
    std::array<std::uint32_t, 3> v = t.v;
    detail::cdt_make_ccw(v, points);
    result.triangles.push_back(v);
  }
  std::sort(result.triangles.begin(), result.triangles.end());
  return result;
}

}  // namespace v1
}  // namespace geom
