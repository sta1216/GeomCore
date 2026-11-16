#include <doctest/doctest.h>

#include <cmath>
#include <numbers>
#include <optional>
#include <string_view>

#include "brep_boolean_curved_fixture.hpp"

#include <geom/brep/boolean_solid_curved.hpp>
#include <geom/brep/primitives.hpp>
#include <geom/brep/validate.hpp>
#include <geom/core/error.hpp>
#include <geom/surfaces/analytic.hpp>

namespace {

constexpr double kPi = std::numbers::pi_v<double>;

std::string_view message_of(const geom::Error& error) { return error.message; }

}  // namespace

// QA-17: the whole-sphere boolean must not depend on a chart-aligned placement.
// The task's original defect appeared when the centre line was 15 degrees off an
// axis; the corpus pins that angle and the 30/10-degree grid pins the quadrants,
// but 37 degrees is not on any of those grids and is the acceptance number for
// the general re-chart. The lens of two unit spheres at centre distance 1.5 is
// 11*pi/96 at every direction; the metric is the cap-sum `sphere_face_volume`,
// not the tessellated volume, so a seam-only face cannot hide in an approximation.
TEST_CASE("QA-17: a sphere lens at an arbitrary placement is rebuilt exactly") {
  const double tolerance = 1e-9;
  const double angle = 37.0 * kPi / 180.0;
  const geom::Result<geom::Brep> south =
      geom::make_sphere(geom::Point3d{0.0, 0.0, 0.0}, 1.0, tolerance);
  const geom::Result<geom::Brep> north = geom::make_sphere(
      geom::Point3d{1.5 * std::cos(angle), 1.5 * std::sin(angle), 0.0}, 1.0, tolerance);
  REQUIRE_MESSAGE(south.has_value(), message_of(south.error()));
  REQUIRE_MESSAGE(north.has_value(), message_of(north.error()));

  const geom::Result<geom::CurvedBooleanResult> lens = geom::boolean_solid_curved(
      south.value(), south.value().solid_ids().front(), north.value(),
      north.value().solid_ids().front(), geom::SolidBooleanOp::intersection,
      geom::CurvedBooleanOptions{tolerance});
  REQUIRE_MESSAGE(lens.has_value(), message_of(lens.error()));
  CHECK(geom::validate_topology(lens.value().brep).ok());
  REQUIRE(lens.value().brep.solid_shells(lens.value().brep.solid_ids().front()).size() == 1);

  const std::optional<double> volume = geom::test::sphere_face_volume(lens.value().brep);
  REQUIRE_MESSAGE(volume.has_value(), "the lens has a face that is not a spherical cap");
  CHECK(*volume == doctest::Approx(11.0 * kPi / 96.0).epsilon(1e-6));

  // The rebuild is deterministic, not merely close: the two runs must agree bit
  // for bit, so a successful run that re-charts about a differently rounded axis
  // cannot pass this test.
  const geom::Result<geom::CurvedBooleanResult> again = geom::boolean_solid_curved(
      south.value(), south.value().solid_ids().front(), north.value(),
      north.value().solid_ids().front(), geom::SolidBooleanOp::intersection,
      geom::CurvedBooleanOptions{tolerance});
  REQUIRE_MESSAGE(again.has_value(), message_of(again.error()));
  const std::optional<double> repeated = geom::test::sphere_face_volume(again.value().brep);
  REQUIRE_MESSAGE(repeated.has_value(), "the repeated lens has a face that is not a spherical cap");
  CHECK(*repeated == *volume);
}
