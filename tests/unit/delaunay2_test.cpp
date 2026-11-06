#include <doctest/doctest.h>

#include <algorithm>
#include <cmath>
#include <vector>

#include <geom/cgeom/delaunay2.hpp>
#include <geom/cgeom/hull2d.hpp>
#include <geom/numeric/predicates2.hpp>

#include <cstdint>
#include <map>
#include <utility>

using geom::v1::convex_hull2d;
using geom::v1::delaunay2;
using geom::v1::Point2d;

namespace {

double incircle_perm(const std::vector<Point2d>& p, std::uint32_t a, std::uint32_t b,
                     std::uint32_t c, std::uint32_t d) {
  return static_cast<double>(geom::v1::sign_incircle(p[a][0], p[a][1], p[b][0], p[b][1],
                                                     p[c][0], p[c][1], p[d][0], p[d][1]));
}

}  // namespace

TEST_CASE("Delaunay satisfies the empty circumcircle property") {
  std::vector<Point2d> pts;
  for (int y = 0; y < 6; ++y) {
    for (int x = 0; x < 6; ++x) {
      // Deterministic perturbation keeps the input in general position
      // (no exact cocircular quadruples), which is the M1 contract.
      const double jx = 0.001 * static_cast<double>((x * 7 + y * 3) % 5);
      const double jy = 0.001 * static_cast<double>((x * 5 + y * 7) % 3);
      pts.push_back(Point2d{static_cast<double>(x) + jx, static_cast<double>(y) + jy});
    }
  }
  pts.push_back(Point2d{2.3, 3.1});
  pts.push_back(Point2d{4.7, 0.4});

  const auto tris = delaunay2(pts);
  REQUIRE_FALSE(tris.empty());
  for (const auto& t : tris) {
    CHECK(geom::v1::sign_orient2d(pts[t[0]][0], pts[t[0]][1], pts[t[1]][0], pts[t[1]][1],
                                  pts[t[2]][0], pts[t[2]][1]) > 0);
    for (std::uint32_t i = 0; i < pts.size(); ++i) {
      if (i == t[0] || i == t[1] || i == t[2]) continue;
      CHECK(incircle_perm(pts, t[0], t[1], t[2], i) <= 0);
    }
  }
}

TEST_CASE("Delaunay triangle count matches Euler's formula") {
  std::vector<Point2d> pts;
  for (int i = 0; i < 20; ++i) {
    pts.push_back(Point2d{std::cos(0.37 * i) * 10.0 + (i % 3) * 0.01,
                          std::sin(0.53 * i) * 10.0 + (i % 5) * 0.02});
  }
  const auto tris = delaunay2(pts);
  const auto hull = convex_hull2d(pts);
  REQUIRE_FALSE(tris.empty());
  // T = 2n - 2 - h for n points in general position, h hull vertices.
  CHECK(tris.size() == 2 * pts.size() - 2 - hull.size());
}

TEST_CASE("Delaunay skips duplicate points and handles tiny input") {
  std::vector<Point2d> dup{{0, 0}, {0, 0}, {1, 0}, {0, 1}, {1, 2}};
  CHECK_FALSE(delaunay2(dup).empty());
  std::vector<Point2d> two{{0, 0}, {1, 0}};
  CHECK(delaunay2(two).empty());
}

namespace {

// Perturbed square grid: spacing 0.1 with deterministic, incommensurate
// offsets small enough to preserve the grid combinatorics but large enough to
// destroy every exact cocircular quadruple.
std::vector<Point2d> perturbed_grid(int n) {
  std::vector<Point2d> pts;
  pts.reserve(static_cast<std::size_t>(n) * static_cast<std::size_t>(n));
  for (int y = 0; y < n; ++y) {
    for (int x = 0; x < n; ++x) {
      const double ox = 0.037 * std::sin(1.7 * x + 0.9 * y + 0.3);
      const double oy = 0.041 * std::cos(2.3 * x - 1.1 * y + 0.7);
      pts.push_back(Point2d{0.1 * x + ox, 0.1 * y + oy});
    }
  }
  return pts;
}

}  // namespace

TEST_CASE("Delaunay perturbed grids satisfy all structural invariants") {
  for (int n : {5, 10, 15}) {
    CAPTURE(n);
    const auto pts = perturbed_grid(n);
    const std::size_t v = pts.size();
    const auto tris = delaunay2(pts);
    const auto hull = convex_hull2d(pts);
    REQUIRE_FALSE(tris.empty());

    // Euler: T = 2V - 2 - h with h hull vertices (convex_hull2d is the oracle).
    CHECK(tris.size() == 2 * v - 2 - hull.size());

    std::map<std::pair<std::uint32_t, std::uint32_t>, int> edge_uses;
    for (const auto& t : tris) {
      // CCW and strictly non-degenerate.
      CHECK(geom::v1::sign_orient2d(pts[t[0]][0], pts[t[0]][1], pts[t[1]][0], pts[t[1]][1],
                                    pts[t[2]][0], pts[t[2]][1]) > 0);
      // Empty circumcircle: no other input point lies strictly inside.
      for (std::uint32_t i = 0; i < pts.size(); ++i) {
        if (i == t[0] || i == t[1] || i == t[2]) continue;
        CHECK(incircle_perm(pts, t[0], t[1], t[2], i) <= 0);
      }
      for (int e = 0; e < 3; ++e) {
        std::uint32_t a = t[e];
        std::uint32_t b = t[(e + 1) % 3];
        if (a > b) std::swap(a, b);
        ++edge_uses[{a, b}];
      }
    }
    for (const auto& [edge, uses] : edge_uses) {
      CAPTURE(edge.first);
      CAPTURE(edge.second);
      CHECK(uses <= 2);
    }
  }
}

TEST_CASE("Delaunay output is deterministic and canonically ordered") {
  const auto pts = perturbed_grid(7);
  const auto first = delaunay2(pts);
  const auto second = delaunay2(pts);
  REQUIRE_FALSE(first.empty());
  CHECK(first == second);
  CHECK(std::is_sorted(first.begin(), first.end()));
  for (std::size_t i = 1; i < first.size(); ++i) {
    CHECK(first[i - 1] != first[i]);
  }
}

TEST_CASE("Delaunay returns empty for degenerate input") {
  CHECK(delaunay2(std::vector<Point2d>{}).empty());
  CHECK(delaunay2(std::vector<Point2d>{{1.0, 2.0}}).empty());
  // Collinear points have no triangulation with positive area.
  CHECK(delaunay2(std::vector<Point2d>{{0, 0}, {1, 1}, {2, 2}, {3, 3}}).empty());
  // Duplicates only reduce to a single point.
  CHECK(delaunay2(std::vector<Point2d>{{5, 5}, {5, 5}, {5, 5}}).empty());
}
