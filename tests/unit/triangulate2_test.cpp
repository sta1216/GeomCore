#include <doctest/doctest.h>

#include <geom/cgeom/polygon2.hpp>
#include <geom/cgeom/triangulate2.hpp>
#include <geom/numeric/predicates2.hpp>

#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <random>
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

// Every emitted triangle must be counter-clockwise and have strictly positive
// area; a count-only check would accept degenerate or inverted output.
void check_ccw(const Polygon2d& p, const Triangulation2& tri) {
  for (const auto& t : tri.triangles) {
    const Point2d& a = p.vertices[t[0]];
    const Point2d& b = p.vertices[t[1]];
    const Point2d& c = p.vertices[t[2]];
    CHECK(geom::v1::sign_orient2d(a[0], a[1], b[0], b[1], c[0], c[1]) > 0);
    CHECK(triangle_area2(p, t) > 0.0);
  }
}

bool point_in_triangle(const Point2d& a, const Point2d& b, const Point2d& c, double x,
                       double y) {
  return geom::v1::sign_orient2d(a[0], a[1], b[0], b[1], x, y) >= 0 &&
         geom::v1::sign_orient2d(b[0], b[1], c[0], c[1], x, y) >= 0 &&
         geom::v1::sign_orient2d(c[0], c[1], a[0], a[1], x, y) >= 0;
}

// Even-odd crossing test across all rings: inside the outer ring and outside
// every hole. Sample points avoid integer coordinates, so no boundary case is
// expected.
bool point_in_polygon(const Polygon2d& p, double x, double y) {
  bool inside = false;
  for (std::uint32_t r = 0; r < p.ring_count(); ++r) {
    const auto ring = p.ring(r);
    const std::size_t n = ring.size();
    for (std::size_t i = 0, j = n - 1; i < n; j = i++) {
      const Point2d& a = ring[i];
      const Point2d& b = ring[j];
      if ((a[1] > y) != (b[1] > y)) {
        const double x_cross = a[0] + (y - a[1]) * (b[0] - a[0]) / (b[1] - a[1]);
        if (x < x_cross) inside = !inside;
      }
    }
  }
  return inside;
}

// Deterministic low-discrepancy sample: no point may be covered by more than one
// triangle, and a point is covered exactly when it lies inside the polygon.
void check_coverage(const Polygon2d& p, const Triangulation2& tri) {
  constexpr int kSamples = 200;
  for (int k = 1; k <= kSamples; ++k) {
    const double fx = std::fmod(0.5 + static_cast<double>(k) * 0.7548776662466927, 1.0);
    const double fy = std::fmod(0.5 + static_cast<double>(k) * 0.5698402909980532, 1.0);
    const double x = 10.0 * fx;
    const double y = 10.0 * fy;
    int hits = 0;
    for (const auto& t : tri.triangles) {
      if (point_in_triangle(p.vertices[t[0]], p.vertices[t[1]], p.vertices[t[2]], x, y)) {
        ++hits;
      }
    }
    INFO("sample " << k << " at (" << x << ", " << y << ") hit " << hits << " triangle(s)");
    CHECK(hits <= 1);
    CHECK((hits == 1) == point_in_polygon(p, x, y));
  }
}

Polygon2d single(std::vector<Point2d> ring) {
  Polygon2d p;
  p.vertices = std::move(ring);
  p.ring_starts = {0, static_cast<std::uint32_t>(p.vertices.size())};
  return p;
}

// 10x10 square with two disjoint axis-aligned unit-square holes whose
// lower-left corners are (ax, ay) and (bx, by).
Polygon2d two_unit_holes(double ax, double ay, double bx, double by) {
  Polygon2d p;
  p.vertices = {{0, 0},  {10, 0}, {10, 10}, {0, 10},
                {ax, ay}, {ax + 1, ay}, {ax + 1, ay + 1}, {ax, ay + 1},
                {bx, by}, {bx + 1, by}, {bx + 1, by + 1}, {bx, by + 1}};
  p.ring_starts = {0, 4, 8, 12};
  return p;
}

}  // namespace

TEST_CASE("triangulate2: square") {
  const Polygon2d p = single({{0, 0}, {1, 0}, {1, 1}, {0, 1}});
  const Triangulation2 tri = triangulate_polygon2(p);
  REQUIRE(tri.ok);
  CHECK(tri.triangles.size() == 2);
  check_ccw(p, tri);
  check_area(p, tri);
}

TEST_CASE("triangulate2: concave L-shape") {
  const Polygon2d p = single({{0, 0}, {2, 0}, {2, 1}, {1, 1}, {1, 2}, {0, 2}});
  const Triangulation2 tri = triangulate_polygon2(p);
  REQUIRE(tri.ok);
  CHECK(tri.triangles.size() == 4);
  check_ccw(p, tri);
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
  check_ccw(p, tri);
  check_area(p, tri);
}

TEST_CASE("triangulate2: reversed rings and collinear vertices") {
  Polygon2d reversed;
  reversed.vertices = {{0, 0}, {0, 1}, {1, 1}, {1, 0}};  // CW outer
  reversed.ring_starts = {0, 4};
  const Triangulation2 tri_cw = triangulate_polygon2(reversed);
  REQUIRE(tri_cw.ok);
  check_ccw(reversed, tri_cw);
  check_area(reversed, tri_cw);

  const Polygon2d collinear = single({{0, 0}, {1, 0}, {2, 0}, {2, 2}, {0, 2}});
  const Triangulation2 tri_col = triangulate_polygon2(collinear);
  REQUIRE(tri_col.ok);
  CHECK(tri_col.triangles.size() == 2);  // collinear vertex (1,0) is dropped
  check_ccw(collinear, tri_col);
  check_area(collinear, tri_col);
}

TEST_CASE("triangulate2: two holes regression (diagonal must not cross a bridge)") {
  // Reviewer repro: hole A at (1,2), hole B at (3,1). With vertex-only ear
  // validation an ear after the second bridge crossed the first bridge edge and
  // the triangulation failed with ok == false.
  const Polygon2d p = two_unit_holes(1, 2, 3, 1);
  const Triangulation2 tri = triangulate_polygon2(p);
  REQUIRE(tri.ok);
  CHECK(tri.triangles.size() == 14);  // n + 2h - 2 = 12 + 4 - 2
  check_ccw(p, tri);
  check_area(p, tri);  // 100 - 1 - 1 = 98
  check_coverage(p, tri);
}

TEST_CASE("triangulate2: randomized disjoint unit-square hole pairs") {
  std::mt19937 rng(20261107u);
  int generated = 0;
  while (generated < 300) {
    const int ax = 1 + static_cast<int>(rng() % 8u);
    const int ay = 1 + static_cast<int>(rng() % 8u);
    const int bx = 1 + static_cast<int>(rng() % 8u);
    const int by = 1 + static_cast<int>(rng() % 8u);
    const int dx = ax - bx;
    const int dy = ay - by;
    // Unit squares overlap or touch exactly when both coordinate gaps are < 2.
    if (dx > -2 && dx < 2 && dy > -2 && dy < 2) continue;
    ++generated;

    const Polygon2d p = two_unit_holes(static_cast<double>(ax), static_cast<double>(ay),
                                       static_cast<double>(bx), static_cast<double>(by));
    const Triangulation2 tri = triangulate_polygon2(p);
    INFO("holes at (" << ax << "," << ay << ") and (" << bx << "," << by << ")");
    REQUIRE(tri.ok);
    CHECK(tri.triangles.size() == 14);
    check_ccw(p, tri);
    check_area(p, tri);
  }
  CHECK(generated == 300);
}

TEST_CASE("triangulate2: rejects malformed input") {
  Polygon2d p;
  p.ring_starts = {0, 2};
  p.vertices = {{0, 0}, {1, 0}};
  const Triangulation2 short_ring = triangulate_polygon2(p);
  CHECK(short_ring.ok == false);
  CHECK(short_ring.triangles.empty());

  p.vertices[0][0] = std::numeric_limits<double>::quiet_NaN();
  p.ring_starts = {0, 3};
  const Triangulation2 non_finite = triangulate_polygon2(p);
  CHECK(non_finite.ok == false);
  CHECK(non_finite.triangles.empty());
}
