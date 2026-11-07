#include <doctest/doctest.h>

#include <geom/cgeom/polygon2.hpp>

#include <limits>

using geom::v1::Errc;
using geom::v1::Polygon2d;
using geom::v1::Point2d;

namespace {

Polygon2d square_with_hole() {
  Polygon2d p;
  // 4x4 outer CCW, 2x2 hole CW.
  p.vertices = {{0, 0}, {4, 0}, {4, 4}, {0, 4},
                {1, 1}, {1, 3}, {3, 3}, {3, 1}};
  p.ring_starts = {0, 4, 8};
  return p;
}

}  // namespace

TEST_CASE("polygon2: ring iteration, orientation and area") {
  const Polygon2d poly = square_with_hole();
  CHECK(poly.ring_count() == 2);
  CHECK(poly.ring(0).size() == 4);
  CHECK(poly.ring(1).size() == 4);
  CHECK(geom::v1::ring_orientation2(poly.ring(0)) == 1);
  CHECK(geom::v1::ring_orientation2(poly.ring(1)) == -1);
  CHECK(geom::v1::signed_area2(poly.ring(0)) == doctest::Approx(16.0));
  CHECK(geom::v1::signed_area2(poly.ring(1)) == doctest::Approx(-4.0));
  CHECK(geom::v1::polygon_area2(poly) == doctest::Approx(12.0));
}

TEST_CASE("polygon2: area is independent of winding") {
  Polygon2d p;
  p.vertices = {{0, 0}, {0, 1}, {1, 1}, {1, 0}};  // CW outer
  p.ring_starts = {0, 4};
  CHECK(geom::v1::ring_orientation2(p.ring(0)) == -1);
  CHECK(geom::v1::polygon_area2(p) == doctest::Approx(1.0));
  CHECK(geom::v1::ring_orientation2(std::vector<Point2d>{{0, 0}, {1, 1}}) == 0);
}

TEST_CASE("polygon2: duplicate consecutive points do not change orientation") {
  Polygon2d p;
  p.vertices = {{0, 0}, {0, 0}, {1, 0}, {1, 1}, {0, 1}};
  p.ring_starts = {0, 5};
  CHECK(geom::v1::ring_orientation2(p.ring(0)) == 1);
}

TEST_CASE("polygon2: structural validation") {
  Polygon2d p;
  p.vertices = {{0, 0}, {1, 0}, {0, 1}};
  p.ring_starts = {0, 3};
  CHECK(geom::v1::validate_polygon2(p) == Errc::none);

  p.ring_starts = {0, 2, 3};
  CHECK(geom::v1::validate_polygon2(p) == Errc::degenerate);

  p.ring_starts = {0, 4};
  CHECK(geom::v1::validate_polygon2(p) == Errc::degenerate);

  p.ring_starts = {1, 3};
  CHECK(geom::v1::validate_polygon2(p) == Errc::degenerate);

  p.ring_starts = {0, 3};
  p.vertices[0][0] = std::numeric_limits<double>::quiet_NaN();
  CHECK(geom::v1::validate_polygon2(p) == Errc::non_finite);

  p.vertices[0][0] = std::numeric_limits<double>::infinity();
  CHECK(geom::v1::validate_polygon2(p) == Errc::non_finite);
}

TEST_CASE("polygon2: set area sums polygons") {
  geom::v1::PolygonSet set;
  set.polygons.push_back(square_with_hole());
  set.polygons.push_back(square_with_hole());
  CHECK(geom::v1::polygon_set_area2(set) == doctest::Approx(24.0));
}
