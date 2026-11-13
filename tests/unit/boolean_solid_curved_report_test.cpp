#include <doctest/doctest.h>

#include <geom/brep/boolean_solid_curved.hpp>
#include <geom/brep/primitives.hpp>

using namespace geom::v1;

TEST_CASE("curved boolean: a failed operation still reports the stage counters") {
  const Result<Brep> a = make_box(Point3d{0, 0, 0}, Point3d{1, 1, 1}, 1e-9);
  const Result<Brep> b = make_box(Point3d{10, 0, 0}, Point3d{11, 1, 1}, 1e-9);
  REQUIRE(a.has_value());
  REQUIRE(b.has_value());
  const SolidId a_solid = a.value().solid_ids().front();
  const SolidId b_solid = b.value().solid_ids().front();
  CurvedBooleanReport report{};
  report.face_pairs = 999;  // must be overwritten by the operation
  const Result<CurvedBooleanResult> result =
      boolean_solid_curved(a.value(), a_solid, b.value(), b_solid, SolidBooleanOp::intersection,
                           CurvedBooleanOptions{}, &report);
  REQUIRE_FALSE(result.has_value());
  CHECK(result.error().code == Errc::numeric_failure);
  CHECK(report.face_pairs == 0);
  CHECK(report.cuts == 0);
  CHECK(report.kept == 0);
}

TEST_CASE("curved boolean: a failure after pairing copies the accumulated counters") {
  const Result<Brep> a = make_box(Point3d{0, 0, 0}, Point3d{1, 1, 1}, 1e-9);
  const Result<Brep> b = make_box(Point3d{3, 0, 0}, Point3d{4, 1, 1}, 1e-9);
  REQUIRE(a.has_value());
  REQUIRE(b.has_value());
  const SolidId a_solid = a.value().solid_ids().front();
  const SolidId b_solid = b.value().solid_ids().front();
  CurvedBooleanReport report{};
  const Result<CurvedBooleanResult> result =
      boolean_solid_curved(a.value(), a_solid, b.value(), b_solid, SolidBooleanOp::union_,
                           CurvedBooleanOptions{}, &report);
  REQUIRE_FALSE(result.has_value());
  CHECK(result.error().code == Errc::unsupported);
  CHECK(report.fragments == 12u);
  CHECK(report.kept == 12u);
}

TEST_CASE("curved boolean: success leaves the failure report untouched") {
  const Result<Brep> a = make_box(Point3d{0, 0, 0}, Point3d{2, 2, 2}, 1e-9);
  const Result<Brep> b = make_box(Point3d{1, 1, 1}, Point3d{3, 3, 3}, 1e-9);
  REQUIRE(a.has_value());
  REQUIRE(b.has_value());
  const SolidId a_solid = a.value().solid_ids().front();
  const SolidId b_solid = b.value().solid_ids().front();
  CurvedBooleanReport report{};
  report.face_pairs = 999;
  report.cuts = 999;
  report.fragments = 999;
  report.kept = 999;
  report.dropped = 999;
  report.flipped = 999;
  report.max_cut_residual = 1.0;
  const Result<CurvedBooleanResult> result =
      boolean_solid_curved(a.value(), a_solid, b.value(), b_solid, SolidBooleanOp::union_,
                           CurvedBooleanOptions{}, &report);
  REQUIRE(result.has_value());
  CHECK((report.face_pairs == 0u && report.cuts == 0u && report.fragments == 0u &&
         report.kept == 0u && report.dropped == 0u && report.flipped == 0u &&
         report.max_cut_residual == 0.0));
}
