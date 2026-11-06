#pragma once

#include <algorithm>
#include <array>
#include <cassert>
#include <cmath>
#include <cstdint>
#include <span>
#include <unordered_map>
#include <utility>
#include <vector>

#include <geom/cgeom/hull2d.hpp>
#include <geom/numeric/predicates2.hpp>
#include <geom/primitives/point.hpp>

namespace geom {
inline namespace v1 {

// 2D Delaunay triangulation: Bowyer-Watson insertion with a linear triangle
// search, O(n^2) worst case; deliberately simple for M1.
//
// Seeding. The implementation does NOT use the classical "super-triangle".
// A super-triangle is only correct while every super vertex stays outside the
// circumcircle of every real triangle; for nearly collinear point sets the
// circumcircles become arbitrarily large compared to the bounding box, so even
// a twenty-times-extent margin leaves super vertices inside them. The result is
// that triangles touching super vertices invade the interior, and dropping them
// at the end punches holes in the output (measured: a perturbed 10x10 grid gave
// T = 183 instead of 184 and lost the hull edge 1-9). Instead the algorithm is
// seeded with the exact triangulation of the convex hull (a fan from one hull
// vertex, see convex_hull2d) which is then made Delaunay by Lawson edge flips.
// The seed has no extra vertices, so nothing has to be dropped and no hole can
// appear; every decision uses the exact sign_incircle / sign_orient2d
// predicates.
//
// Contract:
//   - distinct points in general position: no three collinear and no point
//     exactly on a circumcircle of a triangle of the current triangulation.
//     Degenerate input is outside the M1 contract (no symbolic perturbation
//     yet) and trips a debug assert where it is cheap to detect. Duplicate
//     points are skipped, keeping the first occurrence.
//   - finite coordinates; input.size() <= 2^32 - 1 so indices fit uint32_t.
//   - supported numeric span: exp2(max|coord|) - exp2(min nonzero |coord|) <=
//     1314, the same exponent span sign_incircle can resolve (see
//     predicates2.hpp). Non-finite coordinates and wider spans are outside the
//     supported range and return an empty vector instead of a partial or
//     garbage triangulation; the guard runs in every build.
//
// Output triangles are CCW (verified with sign_orient2d, never a floating cross
// product), index into the input, are strictly non-degenerate, and the list is
// sorted by (a, b, c). The sort makes the geometric result canonical for a
// given point set; the index lists themselves follow the input point order, so
// the same points supplied in a different order yield permuted indices rather
// than identical vectors. Fewer than three distinct points, or all points
// collinear, yield an empty result.
inline std::vector<std::array<std::uint32_t, 3>> delaunay2(std::span<const Point2d> input) {
  assert(input.size() <= 0xFFFFFFFFull);

  // Reject unsupported input in every build, so the guard is testable in Debug.
  constexpr int kMaxSpan = 1314;
  int max_exp = 0;
  int min_exp = 0;
  bool have_exp = false;
  for (const Point2d& q : input) {
    for (const double c : {q[0], q[1]}) {
      if (!std::isfinite(c)) return {};
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
  if (have_exp && max_exp - min_exp > kMaxSpan) return {};

  // Drop exact duplicates, keeping the first occurrence of each point.
  std::vector<std::uint32_t> unique;
  unique.reserve(input.size());
  for (std::uint32_t i = 0; i < input.size(); ++i) {
    bool duplicate = false;
    for (std::uint32_t j : unique) {
      if (input[i][0] == input[j][0] && input[i][1] == input[j][1]) {
        duplicate = true;
        break;
      }
    }
    if (!duplicate) unique.push_back(i);
  }
  if (unique.size() < 3) return {};

  std::vector<Point2d> p;
  p.reserve(unique.size());
  for (std::uint32_t i : unique) p.push_back(input[i]);

  const std::vector<std::uint32_t> hull = convex_hull2d(p);
  if (hull.size() < 3) return {};  // all distinct points are collinear

  struct Tri {
    std::uint32_t v[3];
    bool alive;
  };
  std::vector<Tri> tris;
  tris.reserve(2 * p.size());

  // Make a triangle CCW in place using the exact orientation predicate.
  auto make_ccw = [&](Tri& t) {
    const int s = sign_orient2d(p[t.v[0]][0], p[t.v[0]][1], p[t.v[1]][0], p[t.v[1]][1],
                                p[t.v[2]][0], p[t.v[2]][1]);
    if (s < 0) std::swap(t.v[1], t.v[2]);
    assert(s != 0 && "delaunay2: degenerate triangle");
  };

  // Seed: fan triangulation of the (strictly convex, CCW) hull.
  for (std::size_t k = 1; k + 1 < hull.size(); ++k) {
    Tri t{{hull[0], hull[k], hull[k + 1]}, true};
    make_ccw(t);
    tris.push_back(t);
  }

  // Make the seed Delaunay with Lawson flips. Incremental insertion below is
  // only correct when it starts from a Delaunay triangulation.
  bool flipped = true;
  while (flipped) {
    flipped = false;
    std::unordered_map<std::uint64_t, std::uint32_t> owner;
    owner.reserve(3 * tris.size());
    for (std::uint32_t t = 0; t < tris.size(); ++t) {
      if (!tris[t].alive) continue;
      for (int e = 0; e < 3; ++e) {
        owner[(static_cast<std::uint64_t>(tris[t].v[e]) << 32) | tris[t].v[(e + 1) % 3]] = t;
      }
    }
    for (std::uint32_t t = 0; t < tris.size() && !flipped; ++t) {
      if (!tris[t].alive) continue;
      const Tri T = tris[t];
      for (int e = 0; e < 3 && !flipped; ++e) {
        const std::uint32_t a = T.v[e];
        const std::uint32_t b = T.v[(e + 1) % 3];
        const std::uint32_t c = T.v[(e + 2) % 3];
        const auto it = owner.find((static_cast<std::uint64_t>(b) << 32) | a);
        if (it == owner.end()) continue;  // hull edge, no neighbour
        const std::uint32_t u = it->second;
        if (u <= t) continue;  // each interior edge is considered once
        const Tri U = tris[u];
        std::uint32_t d = U.v[0];
        for (int k = 0; k < 3; ++k) {
          if (U.v[k] != a && U.v[k] != b) d = U.v[k];
        }
        if (sign_incircle(p[a][0], p[a][1], p[b][0], p[b][1], p[c][0], p[c][1], p[d][0],
                          p[d][1]) <= 0) {
          continue;  // legal edge
        }
        Tri t1{{a, c, d}, true};
        Tri t2{{c, b, d}, true};
        make_ccw(t1);
        make_ccw(t2);
        tris[t] = t1;
        tris[u] = t2;
        flipped = true;
      }
    }
  }

  std::vector<bool> inserted(p.size(), false);
  for (std::uint32_t h : hull) inserted[h] = true;

  std::vector<std::uint32_t> bad;
  std::unordered_map<std::uint64_t, std::uint32_t> edge;
  std::vector<std::pair<std::uint32_t, std::uint32_t>> boundary;

  for (std::uint32_t i = 0; i < p.size(); ++i) {
    if (inserted[i]) continue;
    bad.clear();
    for (std::uint32_t t = 0; t < tris.size(); ++t) {
      if (!tris[t].alive) continue;
      const Tri& T = tris[t];
      const int s = sign_incircle(p[T.v[0]][0], p[T.v[0]][1], p[T.v[1]][0], p[T.v[1]][1],
                                  p[T.v[2]][0], p[T.v[2]][1], p[i][0], p[i][1]);
      if (s > 0) bad.push_back(t);
#if !defined(NDEBUG)
      if (s == 0) {
        // General-position precondition violated: the point lies exactly on the
        // triangle's circumcircle, so the insertion order would matter.
        assert(false && "delaunay2: cocircular input is outside the M1 contract");
      }
#endif
    }
    assert(!bad.empty() && "delaunay2: interior point not covered by the triangulation");

    // Directed boundary edges: an edge of a bad triangle without its reverse.
    edge.clear();
    boundary.clear();
    for (std::uint32_t t : bad) {
      const Tri& T = tris[t];
      for (int e = 0; e < 3; ++e) {
        const std::uint64_t key =
            (static_cast<std::uint64_t>(T.v[e]) << 32) | T.v[(e + 1) % 3];
        ++edge[key];
      }
    }
    for (const auto& entry : edge) {
      const std::uint32_t a = static_cast<std::uint32_t>(entry.first >> 32);
      const std::uint32_t b = static_cast<std::uint32_t>(entry.first & 0xFFFFFFFFull);
      const std::uint64_t reverse = (static_cast<std::uint64_t>(b) << 32) | a;
      if (edge.find(reverse) == edge.end()) boundary.emplace_back(a, b);
    }
    // Deterministic triangle creation order (the output is sorted anyway).
    std::sort(boundary.begin(), boundary.end());

    for (std::uint32_t t : bad) tris[t].alive = false;
    for (const auto& [a, b] : boundary) {
      Tri t{{a, b, i}, true};
      make_ccw(t);
      tris.push_back(t);
    }
  }

  std::vector<std::array<std::uint32_t, 3>> out;
  out.reserve(tris.size());
  for (const Tri& T : tris) {
    if (!T.alive) continue;
    const int s = sign_orient2d(p[T.v[0]][0], p[T.v[0]][1], p[T.v[1]][0], p[T.v[1]][1],
                               p[T.v[2]][0], p[T.v[2]][1]);
    if (s == 0) {
      assert(false && "delaunay2: degenerate triangle in output");
      continue;
    }
    if (s > 0) {
      out.push_back({unique[T.v[0]], unique[T.v[1]], unique[T.v[2]]});
    } else {
      out.push_back({unique[T.v[1]], unique[T.v[0]], unique[T.v[2]]});
    }
  }
  std::sort(out.begin(), out.end());
  return out;
}

}  // namespace v1
}  // namespace geom
