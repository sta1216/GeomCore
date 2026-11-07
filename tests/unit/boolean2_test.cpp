#include <doctest/doctest.h>

#include <geom/cgeom/boolean2.hpp>
#include <geom/cgeom/polygon2.hpp>

#include <cmath>
#include <limits>
#include <vector>

using geom::v1::BooleanOp;
using geom::v1::BooleanOptions;
using geom::v1::boolean2;
using geom::v1::Errc;
using geom::v1::FillRule;
using geom::v1::Polygon2d;
using geom::v1::PolygonSet;
using geom::v1::Point2d;
using geom::v1::polygon_area2;
using geom::v1::polygon_set_area2;
using geom::v1::ring_orientation2;

namespace {

Polygon2d rect(double x0, double y0, double x1, double y1) {
  Polygon2d p;
  p.vertices = {{x0, y0}, {x1, y0}, {x1, y1}, {x0, y1}};
  p.ring_starts = {0, 4};
  return p;
}

double area_of(const PolygonSet& set) { return polygon_set_area2(set); }

// One ring appended to a shared vertex pool.
void append_ring(Polygon2d& p, std::vector<Point2d> ring) {
  p.vertices.insert(p.vertices.end(), ring.begin(), ring.end());
  if (p.ring_starts.empty()) p.ring_starts.push_back(0);
  p.ring_starts.push_back(static_cast<std::uint32_t>(p.vertices.size()));
}

}  // namespace

TEST_CASE("boolean2: overlapping squares, four operations") {
  const Polygon2d a = rect(0, 0, 1, 1);
  const Polygon2d b = rect(0.5, 0.5, 1.5, 1.5);
  const Polygon2d ops[2] = {a, b};

  auto uni = boolean2(ops, BooleanOp::union_);
  REQUIRE(uni.has_value());
  CHECK(area_of(uni.value().polygons) == doctest::Approx(1.75));
  CHECK(uni.value().report.output_polygons == 1);

  auto inter = boolean2(ops, BooleanOp::intersection);
  REQUIRE(inter.has_value());
  CHECK(area_of(inter.value().polygons) == doctest::Approx(0.25));

  auto diff = boolean2(ops, BooleanOp::difference);
  REQUIRE(diff.has_value());
  CHECK(area_of(diff.value().polygons) == doctest::Approx(0.75));

  auto xr = boolean2(ops, BooleanOp::symmetric_difference);
  REQUIRE(xr.has_value());
  CHECK(area_of(xr.value().polygons) == doctest::Approx(1.5));
}

TEST_CASE("boolean2: disjoint and nested") {
  const Polygon2d a = rect(0, 0, 1, 1);
  const Polygon2d far = rect(3, 3, 4, 4);
  const Polygon2d ops[2] = {a, far};
  auto uni = boolean2(ops, BooleanOp::union_);
  REQUIRE(uni.has_value());
  CHECK(uni.value().polygons.polygons.size() == 2);
  CHECK(area_of(uni.value().polygons) == doctest::Approx(2.0));
  auto inter = boolean2(ops, BooleanOp::intersection);
  REQUIRE(inter.has_value());
  CHECK(inter.value().polygons.polygons.empty());

  const Polygon2d outer = rect(0, 0, 4, 4);
  const Polygon2d inner = rect(1, 1, 3, 3);
  const Polygon2d nested[2] = {outer, inner};
  auto diff = boolean2(nested, BooleanOp::difference);
  REQUIRE(diff.has_value());
  CHECK(area_of(diff.value().polygons) == doctest::Approx(12.0));
  REQUIRE(diff.value().polygons.polygons.size() == 1);
  CHECK(diff.value().polygons.polygons[0].ring_count() == 2);
}

TEST_CASE("boolean2: shared edge produces one clean rectangle") {
  const Polygon2d a = rect(0, 0, 1, 1);
  const Polygon2d b = rect(1, 0, 2, 1);
  const Polygon2d ops[2] = {a, b};
  auto uni = boolean2(ops, BooleanOp::union_);
  REQUIRE(uni.has_value());
  REQUIRE(uni.value().polygons.polygons.size() == 1);
  const Polygon2d& p = uni.value().polygons.polygons[0];
  CHECK(p.ring_count() == 1);
  CHECK(p.ring(0).size() == 4);  // collinear shared-edge endpoints merged away
  CHECK(area_of(uni.value().polygons) == doctest::Approx(2.0));
}

TEST_CASE("boolean2: fill rules and self-intersection") {
  Polygon2d bowtie;
  bowtie.vertices = {{0, 0}, {2, 2}, {2, 0}, {0, 2}};
  bowtie.ring_starts = {0, 4};
  const Polygon2d ops[1] = {bowtie};
  auto ok = boolean2(ops, BooleanOp::union_);
  REQUIRE(ok.has_value());
  CHECK(std::isfinite(area_of(ok.value().polygons)));

  BooleanOptions strict;
  strict.allow_self_intersections = false;
  auto bad = boolean2(ops, BooleanOp::union_, strict);
  REQUIRE(!bad.has_value());
  CHECK(bad.error().code == Errc::self_intersecting);

  // Fill rules apply per operand, so each coincident unit square has winding
  // +1 inside (odd) and the union of the two is the unit square again.
  const Polygon2d a = rect(0, 0, 1, 1);
  const Polygon2d b = rect(0, 0, 1, 1);
  const Polygon2d same[2] = {a, b};
  BooleanOptions even;
  even.fill_rule = FillRule::even_odd;
  auto x = boolean2(same, BooleanOp::union_, even);
  REQUIRE(x.has_value());
  CHECK(x.value().polygons.polygons.size() == 1);
  CHECK(area_of(x.value().polygons) == doctest::Approx(1.0));
}

TEST_CASE("boolean2: rejects unsupported input") {
  Polygon2d nan = rect(0, 0, 1, 1);
  nan.vertices[0][0] = std::numeric_limits<double>::quiet_NaN();
  const Polygon2d bad[1] = {nan};
  auto r = boolean2(bad, BooleanOp::union_);
  REQUIRE(!r.has_value());
  CHECK(r.error().code == Errc::non_finite);

  Polygon2d huge = rect(0, 0, 1e200, 1);
  huge.vertices[1][1] = 1e-200;
  const Polygon2d wide[1] = {huge};
  auto w = boolean2(wide, BooleanOp::union_);
  REQUIRE(!w.has_value());
  CHECK(w.error().code == Errc::unsupported);

  std::vector<Polygon2d> none;
  auto e = boolean2(none, BooleanOp::union_);
  REQUIRE(!e.has_value());
  CHECK(e.error().code == Errc::degenerate);
}

// Distinct case: every operation is variadic, not just binary.
TEST_CASE("boolean2: three operands and report fields") {
  const Polygon2d a = rect(0, 0, 2, 2);
  const Polygon2d b = rect(1, 0, 3, 2);
  const Polygon2d c = rect(2, 0, 4, 2);
  const Polygon2d ops[3] = {a, b, c};

  auto uni = boolean2(ops, BooleanOp::union_);
  REQUIRE(uni.has_value());
  CHECK(uni.value().polygons.polygons.size() == 1);
  CHECK(uni.value().report.output_polygons == 1);
  CHECK(area_of(uni.value().polygons) == doctest::Approx(8.0));
  CHECK(uni.value().report.input_edges == 12);
  CHECK(uni.value().report.max_residual == 0.0);

  auto inter = boolean2(ops, BooleanOp::intersection);
  REQUIRE(inter.has_value());
  CHECK(inter.value().polygons.polygons.empty());

  auto diff = boolean2(ops, BooleanOp::difference);
  REQUIRE(diff.has_value());
  CHECK(area_of(diff.value().polygons) == doctest::Approx(2.0));

  auto xr = boolean2(ops, BooleanOp::symmetric_difference);
  REQUIRE(xr.has_value());
  CHECK(area_of(xr.value().polygons) == doctest::Approx(4.0));
}

// Distinct case: output rings are oriented outer CCW / holes CW and the report
// counts collinear merges; a lone operand is returned unchanged by difference.
TEST_CASE("boolean2: output orientation, merge count and single operand") {
  const Polygon2d outer = rect(0, 0, 4, 4);
  const Polygon2d inner = rect(1, 1, 3, 3);
  const Polygon2d nested[2] = {outer, inner};
  auto diff = boolean2(nested, BooleanOp::difference);
  REQUIRE(diff.has_value());
  REQUIRE(diff.value().polygons.polygons.size() == 1);
  const Polygon2d& p = diff.value().polygons.polygons[0];
  REQUIRE(p.ring_count() == 2);
  CHECK(ring_orientation2(p.ring(0)) == 1);
  CHECK(ring_orientation2(p.ring(1)) == -1);
  CHECK(polygon_area2(p) == doctest::Approx(12.0));
  CHECK(diff.value().report.output_polygons == 1);
  CHECK(diff.value().report.input_edges == 8);
  CHECK(diff.value().report.merged_vertices == 0);

  const Polygon2d one[1] = {outer};
  auto self = boolean2(one, BooleanOp::difference);
  REQUIRE(self.has_value());
  CHECK(area_of(self.value().polygons) == doctest::Approx(16.0));
  CHECK(self.value().report.output_polygons == 1);
}

// Distinct case: two holes attach to one outer ring, and `difference` removes
// the union of every operand after the first.
TEST_CASE("boolean2: two holes stay with their outer ring") {
  const Polygon2d outer = rect(0, 0, 4, 4);
  const Polygon2d h1 = rect(0.5, 0.5, 1.5, 1.5);
  const Polygon2d h2 = rect(2.5, 2.5, 3.5, 3.5);
  const Polygon2d ops[3] = {outer, h1, h2};
  auto diff = boolean2(ops, BooleanOp::difference);
  REQUIRE(diff.has_value());
  REQUIRE(diff.value().polygons.polygons.size() == 1);
  const Polygon2d& p = diff.value().polygons.polygons[0];
  REQUIRE(p.ring_count() == 3);
  CHECK(ring_orientation2(p.ring(0)) == 1);
  CHECK(ring_orientation2(p.ring(1)) == -1);
  CHECK(ring_orientation2(p.ring(2)) == -1);
  CHECK(polygon_area2(p) == doctest::Approx(14.0));
  CHECK(diff.value().report.input_edges == 12);
}

// Distinct case: `positive` ignores negative winding. Three nested rings of one
// operand give windings +1 (annulus), 0 (between the two holes) and -1 (island);
// `non_zero` and `even_odd` keep the island, `positive` does not.
TEST_CASE("boolean2: positive fill rule drops negative winding") {
  Polygon2d poly;
  append_ring(poly, {{0, 0}, {4, 0}, {4, 4}, {0, 4}});          // CCW outer
  append_ring(poly, {{1, 1}, {3, 1}, {3, 3}, {1, 3}});          // hole, CCW input
  append_ring(poly, {{1.5, 1.5}, {2.5, 1.5}, {2.5, 2.5}, {1.5, 2.5}});  // island
  REQUIRE(poly.ring_count() == 3);
  const Polygon2d ops[1] = {poly};

  auto nz = boolean2(ops, BooleanOp::union_);
  REQUIRE(nz.has_value());
  CHECK(nz.value().polygons.polygons.size() == 2);
  CHECK(area_of(nz.value().polygons) == doctest::Approx(13.0));

  BooleanOptions pos;
  pos.fill_rule = FillRule::positive;
  auto pv = boolean2(ops, BooleanOp::union_, pos);
  REQUIRE(pv.has_value());
  REQUIRE(pv.value().polygons.polygons.size() == 1);
  CHECK(pv.value().polygons.polygons[0].ring_count() == 2);
  CHECK(area_of(pv.value().polygons) == doctest::Approx(12.0));

  BooleanOptions eo;
  eo.fill_rule = FillRule::even_odd;
  auto ev = boolean2(ops, BooleanOp::union_, eo);
  REQUIRE(ev.has_value());
  CHECK(area_of(ev.value().polygons) == doctest::Approx(13.0));
}
