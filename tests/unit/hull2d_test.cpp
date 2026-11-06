#include <doctest/doctest.h>

#include <cmath>
#include <vector>

#include <geom/cgeom/hull2d.hpp>
#include <geom/numeric/predicates2.hpp>

using geom::v1::convex_hull2d;
using geom::v1::Point2d;

TEST_CASE("hull of a square with interior points is the square") {
  std::vector<Point2d> pts{{0, 0}, {1, 0}, {1, 1}, {0, 1},
                           {0.5, 0.5}, {0.25, 0.9}, {0.1, 0.1}};
  const auto hull = convex_hull2d(pts);
  REQUIRE(hull.size() == 4);
  for (std::uint32_t i = 0; i < hull.size(); ++i) {
    const auto& a = pts[hull[i]];
    const auto& b = pts[hull[(i + 1) % hull.size()]];
    const auto& c = pts[hull[(i + 2) % hull.size()]];
    CHECK(geom::v1::sign_orient2d(a[0], a[1], b[0], b[1], c[0], c[1]) > 0);
  }
}

TEST_CASE("hull handles degenerate inputs") {
  CHECK(convex_hull2d(std::vector<Point2d>{}).empty());
  CHECK(convex_hull2d(std::vector<Point2d>{{1, 2}}).size() == 1);

  std::vector<Point2d> line{{0, 0}, {1, 1}, {2, 2}, {3, 3}};
  const auto h = convex_hull2d(line);
  REQUIRE(h.size() == 2);
  CHECK((h[0] == 0 && h[1] == 3));

  std::vector<Point2d> dup{{0, 0}, {0, 0}, {1, 0}, {0, 1}};
  CHECK(convex_hull2d(dup).size() == 3);
}

TEST_CASE("hull of an all-collinear input returns the endpoints in sorted order") {
  // The sort is by (x, then y), so the returned pair is the leftmost endpoint
  // first and the rightmost endpoint second, independent of input order.
  std::vector<Point2d> shuffled{{3, 3}, {0, 0}, {2, 2}, {1, 1}};
  const auto h = convex_hull2d(shuffled);
  REQUIRE(h.size() == 2);
  CHECK(h[0] == 1);  // (0, 0)
  CHECK(h[1] == 0);  // (3, 3)
}

TEST_CASE("hull ordering ties break on x then y") {
  // Equal x values force the secondary y comparison; the first hull index is the
  // lexicographically smallest input point (0, 0).
  std::vector<Point2d> same_x{{0, 0}, {0, 1}, {1, 0}, {1, 1}};
  const auto hx = convex_hull2d(same_x);
  REQUIRE(hx.size() == 4);
  CHECK(hx[0] == 0);
  CHECK((hx[1] == 2 && hx[2] == 3 && hx[3] == 1));

  // Equal y values with a collinear interior point: only the endpoints survive.
  std::vector<Point2d> same_y{{0, 0}, {1, 0}, {2, 0}, {1, 1}};
  const auto hy = convex_hull2d(same_y);
  REQUIRE(hy.size() == 3);
  CHECK((hy[0] == 0 && hy[1] == 2 && hy[2] == 3));
}

TEST_CASE("hull of identical points is a single index") {
  std::vector<Point2d> same{{2, 2}, {2, 2}, {2, 2}};
  const auto h = convex_hull2d(same);
  REQUIRE(h.size() == 1);
  CHECK(h[0] == 0);
}

TEST_CASE("hull keeps the smallest input index of duplicated coordinates") {
  std::vector<Point2d> dup{{1, 0}, {0, 0}, {1, 0}};
  const auto h = convex_hull2d(dup);
  REQUIRE(h.size() == 2);
  CHECK((h[0] == 1 && h[1] == 0));
}

TEST_CASE("hull classifies near-collinear points exactly") {
  // The third point is one ulp above the line y = x: a float cross product
  // computed in doubles may not see it, sign_orient2d does.
  std::vector<Point2d> bent{{0, 0}, {1, 1}, {2, std::nextafter(2.0, 3.0)}};
  CHECK(convex_hull2d(bent).size() == 3);

  std::vector<Point2d> flat{{0, 0}, {1, 1}, {2, 2}};
  CHECK(convex_hull2d(flat).size() == 2);
}
