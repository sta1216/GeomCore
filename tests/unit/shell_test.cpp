#include <doctest/doctest.h>

#include <array>
#include <cmath>
#include <cstdint>
#include <limits>
#include <span>
#include <string>
#include <string_view>
#include <variant>

#include <geom/brep/mass_properties.hpp>
#include <geom/brep/primitives.hpp>
#include <geom/brep/shell.hpp>
#include <geom/brep/validate.hpp>
#include <geom/core/error.hpp>

using namespace geom;

namespace {

Brep box_solid(Point3d lo, Point3d hi) {
  const Result<Brep> box = make_box(lo, hi, 1e-9);
  REQUIRE(box.has_value());
  return box.value();
}

double volume_of(const Brep& brep) {
  const Result<MassProperties> mass = mass_properties(brep);
  REQUIRE(mass.has_value());
  return mass.value().volume;
}

}  // namespace

TEST_CASE("closed shell leaves a uniform cavity") {
  const Brep box = box_solid(Point3d{0.0, 0.0, 0.0}, Point3d{2.0, 3.0, 4.0});
  ShellOptions options;
  options.thickness = 0.25;
  const std::array<FaceId, 0> open{};
  const Result<ShellResult> result = shell_solid(box, open, options);
  REQUIRE(result.has_value());
  CHECK(result.value().cavity_shells == 1u);
  CHECK(result.value().brep.solid_count() == 1u);
  const auto shells = result.value().brep.solid_shells(result.value().brep.solid_ids().front());
  CHECK(shells.size() == 2u);
  CHECK(validate_topology(result.value().brep).ok());
  // 24 - (1.5 * 2.5 * 3.5).
  CHECK(std::abs(volume_of(result.value().brep) - 10.875) <= 1e-9);
}

TEST_CASE("shell rejects a thickness that leaves no interior") {
  const Brep box = box_solid(Point3d{0.0, 0.0, 0.0}, Point3d{2.0, 3.0, 4.0});
  ShellOptions options;
  options.thickness = 1.0;
  const std::array<FaceId, 0> open{};
  const Result<ShellResult> result = shell_solid(box, open, options);
  REQUIRE_FALSE(result.has_value());
  CHECK(result.error().code == Errc::degenerate);
  CHECK(std::string_view(result.error().message).find("leaves no interior") !=
        std::string::npos);
}
