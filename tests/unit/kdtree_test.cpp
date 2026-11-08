#include <doctest/doctest.h>

#include <geom/index/kdtree.hpp>

#include <algorithm>
#include <cmath>
#include <limits>
#include <random>
#include <span>
#include <vector>

using geom::v1::Box3d;
using geom::v1::KdTree;
using geom::v1::Point3d;

namespace {

double dist2(const Point3d& a, const Point3d& b) {
  const double dx = a[0] - b[0];
  const double dy = a[1] - b[1];
  const double dz = a[2] - b[2];
  return dx * dx + dy * dy + dz * dz;
}

}  // namespace

TEST_CASE("kdtree: empty tree") {
  KdTree tree;
  double d2 = 0.0;
  CHECK(tree.size() == 0);
  CHECK(tree.bounds().is_empty());
  CHECK(tree.nearest(Point3d{0, 0, 0}, d2) == 0xFFFFFFFFu);
  CHECK(d2 == std::numeric_limits<double>::infinity());
  std::uint32_t out[3]{};
  CHECK(tree.nearest_k(Point3d{0, 0, 0}, 3, out) == 0);
  CHECK(tree.range(Box3d::empty()).empty());
}

TEST_CASE("kdtree: single point and duplicate coordinates") {
  std::vector<Point3d> pts{{1, 1, 1}, {1, 1, 1}};
  KdTree tree;
  tree.build(pts);
  double d2 = 0.0;
  CHECK(tree.nearest(Point3d{0, 0, 0}, d2) == 0);
  CHECK(d2 == doctest::Approx(3.0));
  CHECK(tree.size() == 2);

  std::vector<std::uint32_t> out(2);
  CHECK(tree.nearest_k(Point3d{0, 0, 0}, 2, out) == 2);
  CHECK(out[0] == 0);
  CHECK(out[1] == 1);
}

TEST_CASE("kdtree: nearest, k-nearest and range match brute force") {
  std::mt19937 rng(12345);
  std::uniform_real_distribution<double> uni(-100.0, 100.0);
  std::vector<Point3d> pts(2000);
  for (Point3d& p : pts) p = Point3d{uni(rng), uni(rng), uni(rng)};

  KdTree tree;
  tree.build(pts);
  CHECK(tree.size() == pts.size());

  for (int q = 0; q < 200; ++q) {
    const Point3d query{uni(rng), uni(rng), uni(rng)};
    double best = std::numeric_limits<double>::infinity();
    std::uint32_t bi = 0xFFFFFFFFu;
    for (std::uint32_t i = 0; i < pts.size(); ++i) {
      const double d = dist2(query, pts[i]);
      if (d < best || (d == best && i < bi)) {
        best = d;
        bi = i;
      }
    }
    double got = 0.0;
    CHECK(tree.nearest(query, got) == bi);
    CHECK(got == best);

    std::vector<std::uint32_t> all(pts.size());
    for (std::uint32_t i = 0; i < pts.size(); ++i) all[i] = i;
    std::sort(all.begin(), all.end(), [&](std::uint32_t a, std::uint32_t b) {
      const double da = dist2(query, pts[a]);
      const double db = dist2(query, pts[b]);
      return da != db ? da < db : a < b;
    });
    std::vector<std::uint32_t> out(9);
    const std::uint32_t n = tree.nearest_k(query, 9, out);
    REQUIRE(n == 9);
    for (std::uint32_t i = 0; i < n; ++i) CHECK(out[i] == all[i]);
  }

  const Box3d box{Point3d{-50, -50, -50}, Point3d{50, 50, 50}};
  std::vector<std::uint32_t> want;
  for (std::uint32_t i = 0; i < pts.size(); ++i) {
    const bool in = pts[i][0] >= box.lo[0] && pts[i][0] <= box.hi[0] &&
                    pts[i][1] >= box.lo[1] && pts[i][1] <= box.hi[1] &&
                    pts[i][2] >= box.lo[2] && pts[i][2] <= box.hi[2];
    if (in) want.push_back(i);
  }
  CHECK(tree.range(box) == want);
}

TEST_CASE("kdtree: tie discrimination picks the smallest index") {
  // 16 points: indices 0..7 at (0,0,0), indices 8..15 at (10,0,0). Leaf size 8
  // puts each group in its own leaf, and the smaller-index leaf is visited
  // second. Every point is at distance 25 from the query, so both the pruning
  // equality rule and the index tie-breaks decide the answer.
  std::vector<Point3d> pts;
  for (int i = 0; i < 8; ++i) pts.push_back(Point3d{0, 0, 0});
  for (int i = 0; i < 8; ++i) pts.push_back(Point3d{10, 0, 0});
  KdTree tree;
  tree.build(pts);

  double d2 = 0.0;
  CHECK(tree.nearest(Point3d{5, 0, 0}, d2) == 0);
  CHECK(d2 == 25.0);

  std::vector<std::uint32_t> out(2);
  CHECK(tree.nearest_k(Point3d{5, 0, 0}, 2, out) == 2);
  CHECK(out[0] == 0);
  CHECK(out[1] == 1);
}

TEST_CASE("kdtree: all duplicate coordinates keep index order") {
  std::vector<Point3d> pts(100, Point3d{3, 3, 3});
  KdTree tree;
  tree.build(pts);

  double d2 = 0.0;
  CHECK(tree.nearest(Point3d{3, 3, 3}, d2) == 0);
  CHECK(d2 == 0.0);

  std::vector<std::uint32_t> out(5);
  CHECK(tree.nearest_k(Point3d{3, 3, 3}, 5, out) == 5);
  CHECK(out[0] == 0);
  CHECK(out[1] == 1);
  CHECK(out[2] == 2);
  CHECK(out[3] == 3);
  CHECK(out[4] == 4);
}

TEST_CASE("kdtree: nearest_k k and output-size boundaries") {
  std::vector<Point3d> pts{{0, 0, 0}, {1, 0, 0}, {2, 0, 0}, {3, 0, 0}};
  KdTree tree;
  tree.build(pts);

  std::vector<std::uint32_t> out(4, 99u);
  CHECK(tree.nearest_k(Point3d{0, 0, 0}, 0, out) == 0);
  CHECK(out[0] == 99u);  // k == 0 writes nothing

  CHECK(tree.nearest_k(Point3d{0, 0, 0}, 10, out) == 4);  // k > n clamps to n
  CHECK(out[0] == 0);
  CHECK(out[1] == 1);
  CHECK(out[2] == 2);
  CHECK(out[3] == 3);

  std::vector<std::uint32_t> small(2);
  CHECK(tree.nearest_k(Point3d{0, 0, 0}, 4, small) == 2);  // out.size() < k
  CHECK(small[0] == 0);
  CHECK(small[1] == 1);

  CHECK(tree.nearest_k(Point3d{0, 0, 0}, 4, std::span<std::uint32_t>{}) == 0);
}

TEST_CASE("kdtree: empty build, rebuild and clearing") {
  KdTree tree;
  tree.build(std::span<const Point3d>{});
  CHECK(tree.size() == 0);
  CHECK(tree.bounds().is_empty());
  double d2 = 0.0;
  CHECK(tree.nearest(Point3d{0, 0, 0}, d2) == 0xFFFFFFFFu);
  CHECK(d2 == std::numeric_limits<double>::infinity());
  CHECK(tree.range(Box3d::empty()).empty());

  const std::vector<Point3d> first{{0, 0, 0}, {100, 0, 0}};
  tree.build(first);
  CHECK(tree.size() == 2);
  CHECK(tree.nearest(Point3d{0, 0, 0}, d2) == 0);

  tree.build(std::span<const Point3d>{});
  CHECK(tree.size() == 0);
  CHECK(tree.bounds().is_empty());
  CHECK(tree.nearest(Point3d{0, 0, 0}, d2) == 0xFFFFFFFFu);
  CHECK(d2 == std::numeric_limits<double>::infinity());
  std::uint32_t out[2]{};
  CHECK(tree.nearest_k(Point3d{0, 0, 0}, 2, out) == 0);
  CHECK(tree.range(Box3d::empty()).empty());

  const std::vector<Point3d> second{{1, 0, 0}, {2, 0, 0}};
  tree.build(second);  // rebuild with a different point set
  CHECK(tree.size() == 2);
  CHECK(tree.nearest(Point3d{0, 0, 0}, d2) == 0);
  CHECK(d2 == 1.0);
  std::vector<std::uint32_t> out2(2);
  CHECK(tree.nearest_k(Point3d{0, 0, 0}, 2, out2) == 2);
  CHECK(out2[0] == 0);
  CHECK(out2[1] == 1);
}

TEST_CASE("kdtree: range boundaries") {
  std::vector<Point3d> pts{{-2, -2, -2}, {0, 0, 0}, {2, 2, 2}};
  KdTree tree;
  tree.build(pts);

  CHECK(tree.range(Box3d::empty()).empty());
  CHECK(tree.range(Box3d{Point3d{10, 10, 10}, Point3d{20, 20, 20}}).empty());

  const std::vector<std::uint32_t> all = tree.range(tree.bounds());
  REQUIRE(all.size() == 3);
  CHECK(all[0] == 0);
  CHECK(all[1] == 1);
  CHECK(all[2] == 2);

  const std::vector<std::uint32_t> inside =
      tree.range(Box3d{Point3d{0, 0, 0}, Point3d{2, 2, 2}});
  REQUIRE(inside.size() == 2);
  CHECK(inside[0] == 1);
  CHECK(inside[1] == 2);
}

TEST_CASE("kdtree: bounds and queries outside the root box") {
  std::vector<Point3d> pts{{1, 2, 3}, {-4, 5, -6}, {0, 0, 0}};
  KdTree tree;
  tree.build(pts);
  REQUIRE(tree.size() == 3);

  const Box3d& b = tree.bounds();
  CHECK(b.lo[0] == -4.0);
  CHECK(b.lo[1] == 0.0);
  CHECK(b.lo[2] == -6.0);
  CHECK(b.hi[0] == 1.0);
  CHECK(b.hi[1] == 5.0);
  CHECK(b.hi[2] == 3.0);

  double d2 = 0.0;
  CHECK(tree.nearest(Point3d{-1000, 0, 0}, d2) == 1);
  CHECK(d2 == 992077.0);

  std::vector<std::uint32_t> out(2);
  CHECK(tree.nearest_k(Point3d{1000, 0, 0}, 2, out) == 2);
  CHECK(out[0] == 0);
  CHECK(out[1] == 2);
}
