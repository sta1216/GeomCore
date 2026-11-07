#include <doctest/doctest.h>

#include <geom/index/kdtree.hpp>

#include <algorithm>
#include <cmath>
#include <limits>
#include <random>
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
  CHECK(tree.bounds().lo[0] > tree.bounds().hi[0]);  // empty box
  CHECK(tree.nearest(Point3d{0, 0, 0}, d2) == 0xFFFFFFFFu);
  CHECK(std::isinf(d2));
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
