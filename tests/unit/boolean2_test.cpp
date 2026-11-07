#include <doctest/doctest.h>

#include <geom/cgeom/boolean2.hpp>
#include <geom/cgeom/polygon2.hpp>

#include <cmath>
#include <limits>
#include <utility>
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
  CHECK(uni.value().report.merged_vertices == 2);
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
  // Both lobes are accepted with opposite winding and become two components.
  REQUIRE(ok.value().polygons.polygons.size() == 2);
  CHECK(area_of(ok.value().polygons) == doctest::Approx(2.0));

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

// Regression C1: one global area scale loses every loop much smaller than the
// largest operand. Each loop is normalized by its own exponent, so a tiny
// disjoint operand or hole survives next to a 1e150 operand. Areas are not
// asserted here: polygon_area2 overflows or underflows at this scale.
TEST_CASE("boolean2: mixed-scale union keeps the tiny disjoint operand") {
  const Polygon2d big = rect(2e150, 0.0, 3e150, 1e150);
  const Polygon2d tiny = rect(0.0, 0.0, 1e-150, 1e-150);
  const Polygon2d ops[2] = {big, tiny};
  auto uni = boolean2(ops, BooleanOp::union_);
  REQUIRE(uni.has_value());
  REQUIRE(uni.value().polygons.polygons.size() == 2);
  CHECK(uni.value().report.output_polygons == 2);
  for (const Polygon2d& p : uni.value().polygons.polygons) {
    CHECK(p.ring_count() == 1);
    CHECK(p.ring(0).size() == 4);
  }
}

TEST_CASE("boolean2: mixed-scale difference keeps the tiny hole") {
  const Polygon2d big = rect(0.0, 0.0, 1e150, 1e150);
  const Polygon2d tiny = rect(1e-149, 1e-149, 2e-149, 2e-149);
  const Polygon2d ops[2] = {big, tiny};
  auto diff = boolean2(ops, BooleanOp::difference);
  REQUIRE(diff.has_value());
  REQUIRE(diff.value().polygons.polygons.size() == 1);
  const Polygon2d& p = diff.value().polygons.polygons[0];
  REQUIRE(p.ring_count() == 2);
  CHECK(p.ring(0).size() == 4);
  CHECK(p.ring(1).size() == 4);
}

TEST_CASE("boolean2: mixed-scale intersection keeps the tiny operand") {
  const Polygon2d big = rect(0.0, 0.0, 1e150, 1e150);
  const Polygon2d tiny = rect(1e-149, 1e-149, 2e-149, 2e-149);
  const Polygon2d ops[2] = {big, tiny};
  auto inter = boolean2(ops, BooleanOp::intersection);
  REQUIRE(inter.has_value());
  REQUIRE(inter.value().polygons.polygons.size() == 1);
  const Polygon2d& p = inter.value().polygons.polygons[0];
  CHECK(p.ring_count() == 1);
  CHECK(p.ring(0).size() == 4);
}

// A ring's first and last segments meet at the closing vertex; strict mode must
// not report that shared endpoint as a self-intersection.
TEST_CASE("boolean2: strict mode accepts wrap-around adjacency") {
  Polygon2d triangle;
  triangle.vertices = {{0.0, 0.0}, {1.0, 0.0}, {0.0, 1.0}};
  triangle.ring_starts = {0, 3};
  const Polygon2d tri[1] = {triangle};
  BooleanOptions strict;
  strict.allow_self_intersections = false;
  auto t = boolean2(tri, BooleanOp::union_, strict);
  REQUIRE(t.has_value());
  CHECK(t.value().polygons.polygons.size() == 1);
  CHECK(area_of(t.value().polygons) == doctest::Approx(0.5));
}

// Crossings between distinct operands are not self-intersections.
TEST_CASE("boolean2: strict mode accepts cross-operand overlap") {
  const Polygon2d a = rect(0.0, 0.0, 1.0, 1.0);
  const Polygon2d b = rect(0.5, 0.5, 1.5, 1.5);
  const Polygon2d ops[2] = {a, b};
  BooleanOptions strict;
  strict.allow_self_intersections = false;
  auto uni = boolean2(ops, BooleanOp::union_, strict);
  REQUIRE(uni.has_value());
  CHECK(area_of(uni.value().polygons) == doctest::Approx(1.75));
}

// An outer ring and one of its hole rings are different rings; their contact at
// a corner is not a self-intersection even in strict mode.
TEST_CASE("boolean2: strict mode accepts a hole touching its outer at a corner") {
  Polygon2d poly = rect(0.0, 0.0, 4.0, 4.0);
  append_ring(poly, {{0.0, 0.0}, {1.0, 0.5}, {0.5, 1.0}});
  REQUIRE(poly.ring_count() == 2);
  const Polygon2d ops[1] = {poly};
  BooleanOptions strict;
  strict.allow_self_intersections = false;
  auto r = boolean2(ops, BooleanOp::union_, strict);
  REQUIRE(r.has_value());
  CHECK(!r.value().polygons.polygons.empty());
}

// Four nested rings with even_odd: the result has two outers (0..10 and 4..6)
// and both contain the 5..5.5 hole, which must attach to the smaller one.
TEST_CASE("boolean2: hole attaches to the smallest containing outer") {
  Polygon2d poly;
  append_ring(poly, {{0.0, 0.0}, {10.0, 0.0}, {10.0, 10.0}, {0.0, 10.0}});
  append_ring(poly, {{2.0, 2.0}, {8.0, 2.0}, {8.0, 8.0}, {2.0, 8.0}});
  append_ring(poly, {{4.0, 4.0}, {6.0, 4.0}, {6.0, 6.0}, {4.0, 6.0}});
  append_ring(poly, {{5.0, 5.0}, {5.5, 5.0}, {5.5, 5.5}, {5.0, 5.5}});
  REQUIRE(poly.ring_count() == 4);
  const Polygon2d ops[1] = {poly};
  BooleanOptions even;
  even.fill_rule = FillRule::even_odd;
  auto r = boolean2(ops, BooleanOp::union_, even);
  REQUIRE(r.has_value());
  REQUIRE(r.value().polygons.polygons.size() == 2);

  const Polygon2d* big = nullptr;
  const Polygon2d* island = nullptr;
  for (const Polygon2d& p : r.value().polygons.polygons) {
    if (polygon_area2(p) > 10.0) {
      big = &p;
    } else {
      island = &p;
    }
  }
  REQUIRE(big != nullptr);
  REQUIRE(island != nullptr);
  REQUIRE(big->ring_count() == 2);
  REQUIRE(island->ring_count() == 2);
  CHECK(polygon_area2(*big) == doctest::Approx(64.0));
  CHECK(polygon_area2(*island) == doctest::Approx(3.75));
  CHECK(big->ring(1).size() == 4);
  CHECK(island->ring(1).size() == 4);
}

// Regression: the per-loop area exponent must include the 2*loop_exp term.
// All coordinates are <= 10, but the loops sit at different magnitudes. If the
// term is dropped, the small hole's area exponent is understated and it looks
// larger than the small island, so it attaches to the big outer ring instead.
TEST_CASE("boolean2: per-loop area exponent preserves nested hole ownership") {
  Polygon2d poly;
  append_ring(poly, {{0.0, 0.0}, {10.0, 0.0}, {10.0, 10.0}, {0.0, 10.0}});  // outer, CCW
  append_ring(poly, {{0.5, 0.5}, {0.5, 8.0}, {8.0, 8.0}, {8.0, 0.5}});      // hole, CW
  append_ring(poly, {{1.0, 1.0}, {7.0, 1.0}, {7.0, 7.0}, {1.0, 7.0}});      // island, CCW
  append_ring(poly, {{2.0, 2.0}, {2.0, 3.0}, {3.0, 3.0}, {3.0, 2.0}});      // hole, CW
  REQUIRE(poly.ring_count() == 4);
  const Polygon2d ops[1] = {poly};
  BooleanOptions even;
  even.fill_rule = FillRule::even_odd;
  auto r = boolean2(ops, BooleanOp::union_, even);
  REQUIRE(r.has_value());
  REQUIRE(r.value().polygons.polygons.size() == 2);

  const Polygon2d* outer = nullptr;
  const Polygon2d* island = nullptr;
  for (const Polygon2d& p : r.value().polygons.polygons) {
    if (polygon_area2(p) > 40.0) {
      outer = &p;
    } else {
      island = &p;
    }
  }
  REQUIRE(outer != nullptr);
  REQUIRE(island != nullptr);
  REQUIRE(outer->ring_count() == 2);
  REQUIRE(island->ring_count() == 2);
  CHECK(polygon_area2(*outer) == doctest::Approx(43.75).epsilon(1e-12));
  CHECK(polygon_area2(*island) == doctest::Approx(35.0).epsilon(1e-12));
}

// A plain square has no self-intersection; strict mode must accept it.
TEST_CASE("boolean2: strict mode accepts a simple square") {
  const Polygon2d a = rect(0.0, 0.0, 1.0, 1.0);
  const Polygon2d ops[1] = {a};
  BooleanOptions strict;
  strict.allow_self_intersections = false;
  auto r = boolean2(ops, BooleanOp::union_, strict);
  REQUIRE(r.has_value());
  CHECK(r.value().polygons.polygons.size() == 1);
  CHECK(area_of(r.value().polygons) == doctest::Approx(1.0));
}

namespace {

Polygon2d polygon_of(std::vector<Point2d> vertices) {
  Polygon2d p;
  p.vertices = std::move(vertices);
  p.ring_starts = {0, static_cast<std::uint32_t>(p.vertices.size())};
  return p;
}

}  // namespace

// The kite's leftmost vertex owns a lower-half-plane edge. Arrangement winding
// seeding used to invert the interior/exterior, so this simple polygon unioned
// to nothing and its intersection with itself was empty as well.
TEST_CASE("boolean2: kite keeps its topology for single-operand ops") {
  const Polygon2d kite = polygon_of({{0, 0}, {2, 1}, {0, 4}, {-2, 1}});
  const Polygon2d ops[1] = {kite};

  const auto uni = boolean2(ops, BooleanOp::union_);
  REQUIRE(uni.has_value());
  REQUIRE(uni.value().polygons.polygons.size() == 1);
  CHECK(uni.value().polygons.polygons[0].ring_count() == 1);
  CHECK(polygon_area2(uni.value().polygons.polygons[0]) == doctest::Approx(8.0));

  const auto self = boolean2(ops, BooleanOp::intersection);
  REQUIRE(self.has_value());
  REQUIRE(self.value().polygons.polygons.size() == 1);
  CHECK(polygon_area2(self.value().polygons.polygons[0]) == doctest::Approx(8.0));

  const Polygon2d pair[2] = {kite, rect(100.0, 100.0, 101.0, 101.0)};
  const auto far = boolean2(pair, BooleanOp::union_);
  REQUIRE(far.has_value());
  REQUIRE(far.value().polygons.polygons.size() == 2);
  CHECK(area_of(far.value().polygons) == doctest::Approx(9.0));
}

// One simple convex polygon in, one polygon out: both orientations of the
// leftmost-vertex edge fan are covered (the pentagon has a lower-half edge, and
// the regular 100-gon has one on each side of the +x axis).
TEST_CASE("boolean2: single convex polygon unions to itself") {
  const Polygon2d pentagon =
      polygon_of({{0, 0}, {3, -1}, {4, 1}, {2, 3}, {-1, 1}});
  const Polygon2d pent_ops[1] = {pentagon};
  const auto pent = boolean2(pent_ops, BooleanOp::union_);
  REQUIRE(pent.has_value());
  REQUIRE(pent.value().polygons.polygons.size() == 1);
  CHECK(pent.value().polygons.polygons[0].ring_count() == 1);
  CHECK(polygon_area2(pent.value().polygons.polygons[0]) == doctest::Approx(11.0));

  std::vector<Point2d> circle;
  for (int i = 0; i < 100; ++i) {
    const double angle = 2.0 * 3.14159265358979323846 * i / 100.0;
    circle.push_back({std::cos(angle), std::sin(angle)});
  }
  const Polygon2d polygon = polygon_of(circle);
  const Polygon2d polygon_ops[1] = {polygon};
  const auto round = boolean2(polygon_ops, BooleanOp::union_);
  REQUIRE(round.has_value());
  REQUIRE(round.value().polygons.polygons.size() == 1);
  CHECK(round.value().polygons.polygons[0].ring_count() == 1);
  CHECK(polygon_area2(round.value().polygons.polygons[0]) ==
        doctest::Approx(3.139525976));
}


// Ring areas of 1e200 squares are 1e400, above the double range. Per-loop
// scaling must still recover the orientation and nesting of every ring. The
// island inside the hole makes the test sensitive to an unscaled area: NaN
// orientation happens to read as a hole, so only the polygon count shows that
// the island was lost.
TEST_CASE("boolean2: per-loop scaling survives ring areas above the double range") {
  Polygon2d poly;
  append_ring(poly, {{0.0, 0.0}, {1e200, 0.0}, {1e200, 1e200}, {0.0, 1e200}});  // CCW
  append_ring(poly, {{0.5e200, 0.5e200},
                     {0.5e200, 0.9e200},
                     {0.9e200, 0.9e200},
                     {0.9e200, 0.5e200}});  // CW hole
  append_ring(poly, {{0.55e200, 0.55e200},
                     {0.85e200, 0.55e200},
                     {0.85e200, 0.85e200},
                     {0.55e200, 0.85e200}});  // CCW island inside the hole
  const Polygon2d ops[1] = {poly};

  const auto out = boolean2(ops, BooleanOp::union_);
  REQUIRE(out.has_value());
  REQUIRE(out.value().polygons.polygons.size() == 2);
  const Polygon2d* outer = nullptr;
  const Polygon2d* island = nullptr;
  for (const Polygon2d& candidate : out.value().polygons.polygons) {
    if (candidate.ring_count() == 2) outer = &candidate;
    if (candidate.ring_count() == 1) island = &candidate;
  }
  REQUIRE(outer != nullptr);
  REQUIRE(island != nullptr);
  CHECK(outer->ring(0).size() == 4);
  CHECK(outer->ring(1).size() == 4);
  CHECK(ring_orientation2(outer->ring(0)) == 1);
  CHECK(ring_orientation2(outer->ring(1)) == -1);
  CHECK(island->ring(0).size() == 4);
  CHECK(ring_orientation2(island->ring(0)) == 1);
}

// The horizontal segment's endpoint difference overflows double. Before the
// split fix the T junction at (0, 0) was not recorded (max_residual = +inf) and
// the intersection came out empty instead of the unit square.
TEST_CASE("boolean2: overflowing endpoint difference keeps the T junction topology") {
  const Polygon2d huge = rect(-1e308, 0.0, 1e308, 1e308);
  const Polygon2d band = rect(0.0, -1.0, 1.0, 1.0);
  const Polygon2d ops[2] = {huge, band};

  const auto inter = boolean2(ops, BooleanOp::intersection);
  REQUIRE(inter.has_value());
  REQUIRE(inter.value().polygons.polygons.size() == 1);
  const Polygon2d& unit = inter.value().polygons.polygons[0];
  CHECK(unit.ring_count() == 1);
  CHECK(unit.vertices.size() == 4);
  CHECK(inter.value().report.max_residual == 0.0);
  CHECK(polygon_area2(unit) == doctest::Approx(1.0));

  const auto uni = boolean2(ops, BooleanOp::union_);
  REQUIRE(uni.has_value());
  REQUIRE(uni.value().polygons.polygons.size() == 1);
  CHECK(uni.value().polygons.polygons[0].ring_count() == 1);
  CHECK(uni.value().polygons.polygons[0].vertices.size() == 8);

  const auto diff = boolean2(ops, BooleanOp::difference);
  REQUIRE(diff.has_value());
  REQUIRE(diff.value().polygons.polygons.size() == 1);
  CHECK(diff.value().polygons.polygons[0].ring_count() == 1);
  CHECK(diff.value().polygons.polygons[0].vertices.size() == 8);
}
