#include <doctest/doctest.h>

#include <geom/cgeom/cdt2.hpp>
#include <geom/numeric/predicates2.hpp>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <limits>
#include <unordered_map>
#include <vector>

using geom::v1::CdtConstraint;
using geom::v1::cdt2;
using geom::v1::Errc;
using geom::v1::Point2d;

namespace {

bool has_edge(const std::vector<std::array<std::uint32_t, 3>>& tris,
              std::uint32_t a, std::uint32_t b) {
  for (const auto& t : tris) {
    for (int e = 0; e < 3; ++e) {
      const std::uint32_t u = t[e], v = t[(e + 1) % 3];
      if ((u == a && v == b) || (u == b && v == a)) return true;
    }
  }
  return false;
}

int edge_use(const std::vector<std::array<std::uint32_t, 3>>& tris,
             std::uint32_t a, std::uint32_t b) {
  int n = 0;
  for (const auto& t : tris) {
    for (int e = 0; e < 3; ++e) {
      const std::uint32_t u = t[e], v = t[(e + 1) % 3];
      if ((u == a && v == b) || (u == b && v == a)) ++n;
    }
  }
  return n;
}

// Deterministic general-position 6x6 grid. The rows are only slightly
// perturbed, so a long chord such as 0-4 crosses several interior edges.
std::vector<Point2d> perturbed_grid(int n) {
  std::vector<Point2d> pts;
  pts.reserve(static_cast<std::size_t>(n) * static_cast<std::size_t>(n));
  for (int y = 0; y < n; ++y) {
    for (int x = 0; x < n; ++x) {
      const int i = y * n + x;
      pts.push_back(Point2d{static_cast<double>(x) + 0.031 * std::sin(1.0 + 0.37 * i),
                            static_cast<double>(y) + 0.029 * std::cos(0.5 + 0.53 * i)});
    }
  }
  return pts;
}

}  // namespace

TEST_CASE("cdt2: recovers either diagonal of a convex quadrilateral") {
  const std::vector<Point2d> pts{{0, 0}, {4, 0}, {3, 2}, {0, 3}};
  const CdtConstraint c02[1] = {{0, 2}};
  auto r02 = cdt2(pts, c02);
  REQUIRE(r02.has_value());
  CHECK(has_edge(r02.value().triangles, 0, 2));
  CHECK(edge_use(r02.value().triangles, 0, 2) == 2);

  const CdtConstraint c13[1] = {{1, 3}};
  auto r13 = cdt2(pts, c13);
  REQUIRE(r13.has_value());
  CHECK(has_edge(r13.value().triangles, 1, 3));
  CHECK(edge_use(r13.value().triangles, 1, 3) == 2);
}

TEST_CASE("cdt2: constrained edge on the hull and duplicate constraints") {
  const std::vector<Point2d> pts{{0, 0}, {4, 0}, {3, 2}, {0, 3}};
  const CdtConstraint cs[3] = {{0, 1}, {1, 0}, {0, 1}};
  auto r = cdt2(pts, cs);
  REQUIRE(r.has_value());
  CHECK(edge_use(r.value().triangles, 0, 1) == 1);  // hull edge
}

TEST_CASE("cdt2: rejects crossing constraints and bad indices") {
  const std::vector<Point2d> pts{{0, 0}, {4, 0}, {3, 2}, {0, 3}};
  const CdtConstraint cross[2] = {{0, 2}, {1, 3}};
  auto r = cdt2(pts, cross);
  REQUIRE(!r.has_value());
  CHECK(r.error().code == Errc::self_intersecting);

  const CdtConstraint bad[1] = {{0, 9}};
  auto b = cdt2(pts, bad);
  REQUIRE(!b.has_value());
  CHECK(b.error().code == Errc::index_out_of_range);

  const CdtConstraint degenerate[1] = {{2, 2}};
  auto d = cdt2(pts, degenerate);
  REQUIRE(!d.has_value());
  CHECK(d.error().code == Errc::degenerate);
}

TEST_CASE("cdt2: rejects collinear overlap and T-junction constraints") {
  {
    // {0,1} lies inside {0,2}: collinear overlap through the shared endpoint 0.
    const std::vector<Point2d> pts{{0, 0}, {1, 0}, {2, 0}, {0, 2}};
    const CdtConstraint cs[2] = {{0, 2}, {0, 1}};
    auto r = cdt2(pts, cs);
    REQUIRE(!r.has_value());
    CHECK(r.error().code == Errc::self_intersecting);
  }
  {
    // The endpoint 2 of {2,3} touches the interior of {0,1} without sharing an
    // index, so it is not a legal shared-endpoint meeting.
    const std::vector<Point2d> pts{{0, 0}, {4, 0}, {2, 0}, {2, 2}};
    const CdtConstraint cs[2] = {{0, 1}, {2, 3}};
    auto r = cdt2(pts, cs);
    REQUIRE(!r.has_value());
    CHECK(r.error().code == Errc::self_intersecting);
  }
}

TEST_CASE("cdt2: non-constrained edges stay Delaunay") {
  std::vector<Point2d> pts;
  for (int y = 0; y < 6; ++y) {
    for (int x = 0; x < 6; ++x) {
      // NOTE: the brief's x + 0.013*y / y + 0.017*x is an affine shear of the
      // integer grid, so it preserves collinearity and violates delaunay2's
      // general-position precondition (delaunay2 asserts on this very grid in
      // Debug). Use a genuine, deterministic perturbation instead.
      const int i = y * 6 + x;
      const double jx = 0.031 * std::sin(1.0 + 0.37 * i);
      const double jy = 0.029 * std::cos(0.5 + 0.53 * i);
      pts.push_back(Point2d{static_cast<double>(x) + jx, static_cast<double>(y) + jy});
    }
  }
  const CdtConstraint cs[2] = {{0, 7}, {2, 9}};
  auto r = cdt2(pts, cs);
  REQUIRE(r.has_value());
  std::vector<bool> constrained(pts.size() * pts.size(), false);
  constrained[0 * pts.size() + 7] = constrained[7 * pts.size() + 0] = true;
  constrained[2 * pts.size() + 9] = constrained[9 * pts.size() + 2] = true;
  for (const auto& t : r.value().triangles) {
    for (int e = 0; e < 3; ++e) {
      const std::uint32_t a = t[e], b = t[(e + 1) % 3], c = t[(e + 2) % 3];
      if (constrained[a * pts.size() + b]) continue;
      // Opposite vertex across edge (a,b): find the other triangle.
      std::uint32_t d = 0xFFFFFFFFu;
      for (const auto& u : r.value().triangles) {
        for (int k = 0; k < 3; ++k) {
          const std::uint32_t p = u[k], q = u[(k + 1) % 3];
          if ((p == b && q == a) || (p == a && q == b)) d = u[(k + 2) % 3];
        }
      }
      if (d == 0xFFFFFFFFu) continue;  // hull edge
      CHECK(geom::v1::sign_incircle(pts[a][0], pts[a][1], pts[b][0], pts[b][1],
                                    pts[c][0], pts[c][1], pts[d][0], pts[d][1]) <= 0);
    }
  }
}

TEST_CASE("cdt2: output is a manifold triangulation (Euler and edge use)") {
  std::vector<Point2d> pts;
  for (int y = 0; y < 4; ++y) {
    for (int x = 0; x < 4; ++x) {
      const int i = y * 4 + x;
      pts.push_back(Point2d{static_cast<double>(x) + 0.031 * std::sin(1.0 + 0.37 * i),
                            static_cast<double>(y) + 0.029 * std::cos(0.5 + 0.53 * i)});
    }
  }
  const CdtConstraint cs[5] = {{0, 5}, {2, 7}, {5, 9}, {8, 13}, {10, 15}};
  auto r = cdt2(pts, cs);
  REQUIRE(r.has_value());
  const auto& tris = r.value().triangles;

  std::unordered_map<std::uint64_t, int> directed_use;
  std::unordered_map<std::uint64_t, int> undirected_use;
  std::vector<bool> used(pts.size(), false);
  for (const auto& t : tris) {
    CHECK(geom::v1::sign_orient2d(pts[t[0]][0], pts[t[0]][1], pts[t[1]][0], pts[t[1]][1],
                                  pts[t[2]][0], pts[t[2]][1]) > 0);
    for (int e = 0; e < 3; ++e) {
      const std::uint32_t a = t[e], b = t[(e + 1) % 3];
      used[a] = true;
      ++directed_use[(static_cast<std::uint64_t>(a) << 32) | b];
      ++undirected_use[(static_cast<std::uint64_t>(std::min(a, b)) << 32) | std::max(a, b)];
    }
  }
  for (const auto& [key, count] : directed_use) {
    (void)key;
    CHECK(count == 1);  // a directed edge belongs to at most one triangle
  }
  std::size_t hull_edges = 0;
  for (const auto& [key, count] : undirected_use) {
    (void)key;
    CHECK(count >= 1);
    CHECK(count <= 2);
    if (count == 1) ++hull_edges;
  }
  for (const bool u : used) CHECK(u);
  CHECK(has_edge(tris, 0, 5));
  CHECK(has_edge(tris, 2, 7));
  CHECK(has_edge(tris, 5, 9));
  CHECK(has_edge(tris, 8, 13));
  CHECK(has_edge(tris, 10, 15));
  // V - E + F = 2 with F = T + 1 for a triangulation of a disk.
  CHECK(undirected_use.size() == pts.size() + tris.size() - 1);
  CHECK(3 * tris.size() == 2 * undirected_use.size() - hull_edges);
}

// Constraint {0,4} is not a base Delaunay edge on this 6x6 grid, so recovery
// must delete the crossed cavity and retriangulate it. This failed with
// numeric_failure while the walk could step back across its entry edge.
TEST_CASE("cdt2: recovers a multi-triangle constraint chord") {
  const std::vector<Point2d> pts = perturbed_grid(6);
  const CdtConstraint cs[1] = {{0, 4}};
  auto r = cdt2(pts, cs);
  REQUIRE(r.has_value());
  const auto& tris = r.value().triangles;
  CHECK(edge_use(tris, 0, 4) == 2);

  std::unordered_map<std::uint64_t, int> directed;
  for (const auto& t : tris) {
    for (int e = 0; e < 3; ++e) {
      const std::uint32_t a = t[e], b = t[(e + 1) % 3];
      ++directed[(static_cast<std::uint64_t>(a) << 32) | b];
    }
  }
  CHECK(directed[(static_cast<std::uint64_t>(0) << 32) | 4] == 1);
  CHECK(directed[(static_cast<std::uint64_t>(4) << 32) | 0] == 1);
  for (const auto& [key, count] : directed) {
    (void)key;
    CHECK(count == 1);
  }
}

// Non-vacuous Delaunay check: the constraint above really did cross edges, so
// every unconstrained interior edge must still be locally Delaunay after the
// restore pass. A skipped restore leaves illegal edges here (this chord's
// ear-clipped cavity is not Delaunay, unlike the short 0-4 chord).
TEST_CASE("cdt2: constrained triangulation stays Delaunay away from constraints") {
  const std::vector<Point2d> pts = perturbed_grid(6);
  const CdtConstraint cs[1] = {{31, 35}};
  auto r = cdt2(pts, cs);
  REQUIRE(r.has_value());
  const auto& tris = r.value().triangles;
  CHECK(has_edge(tris, 31, 35));

  std::size_t checked = 0;
  for (const auto& t : tris) {
    for (int e = 0; e < 3; ++e) {
      const std::uint32_t a = t[e], b = t[(e + 1) % 3], c = t[(e + 2) % 3];
      if (a > b) continue;  // each undirected edge once
      if ((a == 31 && b == 35) || (a == 35 && b == 31)) continue;
      std::uint32_t d = 0xFFFFFFFFu;
      for (const auto& u : tris) {
        for (int k = 0; k < 3; ++k) {
          if (u[k] == b && u[(k + 1) % 3] == a) d = u[(k + 2) % 3];
        }
      }
      if (d == 0xFFFFFFFFu) continue;  // hull edge
      ++checked;
      CHECK(geom::v1::sign_incircle(pts[a][0], pts[a][1], pts[b][0], pts[b][1],
                                    pts[c][0], pts[c][1], pts[d][0], pts[d][1]) <= 0);
    }
  }
  CHECK(checked > 0);
}

TEST_CASE("cdt2: duplicate coordinates map to the surviving representative") {
  const std::vector<Point2d> pts{{0, 0}, {2, 0}, {0, 2}, {0, 0}};  // 3 duplicates 0
  const CdtConstraint same[1] = {{0, 3}};
  auto bad = cdt2(pts, same);
  REQUIRE(!bad.has_value());
  CHECK(bad.error().code == Errc::degenerate);

  const CdtConstraint mapped[1] = {{3, 1}};
  auto ok = cdt2(pts, mapped);
  REQUIRE(ok.has_value());
  CHECK(has_edge(ok.value().triangles, 0, 1));
}

TEST_CASE("cdt2: rejects exponent spans beyond the accepted range") {
  const std::vector<Point2d> pts{{0.0, 0.0}, {1e200, 0.0}, {0.0, 1e-200}};
  const std::vector<CdtConstraint> none;
  auto r = cdt2(pts, none);
  REQUIRE(!r.has_value());
  CHECK(r.error().code == Errc::unsupported);
}
