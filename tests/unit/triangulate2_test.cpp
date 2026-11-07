#include <doctest/doctest.h>

#include <geom/cgeom/polygon2.hpp>
#include <geom/cgeom/triangulate2.hpp>
#include <geom/numeric/predicates2.hpp>

#include <array>
#include <cmath>
#include <cstdint>
#include <limits>
#include <utility>
#include <vector>

using geom::v1::Point2d;
using geom::v1::Polygon2d;
using geom::v1::Triangulation2;
using geom::v1::polygon_area2;
using geom::v1::triangulate_polygon2;

namespace {

double triangle_area2(const Polygon2d& p, const std::array<std::uint32_t, 3>& t) {
  return std::fabs(geom::v1::signed_area2(
      std::vector<Point2d>{p.vertices[t[0]], p.vertices[t[1]], p.vertices[t[2]]}));
}

double covered_area(const Polygon2d& p, const Triangulation2& tri) {
  double sum = 0.0;
  for (const auto& t : tri.triangles) sum += triangle_area2(p, t);
  return sum;
}

// The area contract is an explicit 1e-12 relative comparison, not an implicit
// default tolerance.
void check_area(const Polygon2d& p, const Triangulation2& tri) {
  const double expected = polygon_area2(p);
  const double actual = covered_area(p, tri);
  CHECK(std::fabs(actual - expected) <= 1e-12 * std::fabs(expected));
}

Polygon2d single(std::vector<Point2d> ring) {
  Polygon2d p;
  p.vertices = std::move(ring);
  p.ring_starts = {0, static_cast<std::uint32_t>(p.vertices.size())};
  return p;
}

}  // namespace

TEST_CASE("triangulate2: square") {
  const Polygon2d p = single({{0, 0}, {1, 0}, {1, 1}, {0, 1}});
  const Triangulation2 tri = triangulate_polygon2(p);
  REQUIRE(tri.ok);
  CHECK(tri.triangles.size() == 2);
  for (const auto& t : tri.triangles) {
    CHECK(geom::v1::sign_orient2d(p.vertices[t[0]][0], p.vertices[t[0]][1],
                                  p.vertices[t[1]][0], p.vertices[t[1]][1],
                                  p.vertices[t[2]][0], p.vertices[t[2]][1]) > 0);
  }
  check_area(p, tri);
}

TEST_CASE("triangulate2: concave L-shape") {
  const Polygon2d p = single({{0, 0}, {2, 0}, {2, 1}, {1, 1}, {1, 2}, {0, 2}});
  const Triangulation2 tri = triangulate_polygon2(p);
  REQUIRE(tri.ok);
  CHECK(tri.triangles.size() == 4);
  check_area(p, tri);
}

TEST_CASE("triangulate2: polygon with a hole") {
  Polygon2d p;
  p.vertices = {{0, 0}, {4, 0}, {4, 4}, {0, 4},
                {1, 1}, {1, 3}, {3, 3}, {3, 1}};
  p.ring_starts = {0, 4, 8};
  const Triangulation2 tri = triangulate_polygon2(p);
  REQUIRE(tri.ok);
  CHECK(tri.triangles.size() == 8);  // n + 2h - 2
  check_area(p, tri);
}

TEST_CASE("triangulate2: reversed rings and collinear vertices") {
  Polygon2d reversed;
  reversed.vertices = {{0, 0}, {0, 1}, {1, 1}, {1, 0}};  // CW outer
  reversed.ring_starts = {0, 4};
  const Triangulation2 tri_cw = triangulate_polygon2(reversed);
  REQUIRE(tri_cw.ok);
  check_area(reversed, tri_cw);

  const Polygon2d collinear = single({{0, 0}, {1, 0}, {2, 0}, {2, 2}, {0, 2}});
  const Triangulation2 tri_col = triangulate_polygon2(collinear);
  REQUIRE(tri_col.ok);
  CHECK(tri_col.triangles.size() == 2);  // collinear vertex (1,0) is dropped
  check_area(collinear, tri_col);
}

TEST_CASE("triangulate2: rejects malformed input") {
  Polygon2d p;
  p.ring_starts = {0, 2};
  p.vertices = {{0, 0}, {1, 0}};
  CHECK(triangulate_polygon2(p).ok == false);

  p.vertices[0][0] = std::numeric_limits<double>::quiet_NaN();
  p.ring_starts = {0, 3};
  CHECK(triangulate_polygon2(p).ok == false);
}
