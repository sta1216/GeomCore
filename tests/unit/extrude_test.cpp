#include <doctest/doctest.h>

#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <map>
#include <utility>
#include <vector>

#include <geom/brep/extrude.hpp>
#include <geom/brep/mass_properties.hpp>
#include <geom/brep/tessellate.hpp>
#include <geom/brep/topology.hpp>
#include <geom/brep/validate.hpp>
#include <geom/cgeom/polygon2.hpp>
#include <geom/core/result.hpp>
#include <geom/primitives/point.hpp>
#include <geom/surfaces/plane.hpp>

namespace {

using namespace geom;

Polygon2d square_profile() {
  Polygon2d profile;
  profile.vertices = {Point2d{0.0, 0.0}, Point2d{2.0, 0.0}, Point2d{2.0, 2.0},
                      Point2d{0.0, 2.0}};
  profile.ring_starts = {0u, 4u};
  return profile;
}

double volume_of(const Brep& brep) {
  const Result<MassProperties> props = mass_properties(brep);
  REQUIRE(props.has_value());
  return props.value().volume;
}

std::size_t boundary_edges(const TriMesh& mesh) {
  std::map<std::pair<std::uint32_t, std::uint32_t>, int> uses;
  for (const std::array<std::uint32_t, 3>& triangle : mesh.triangles) {
    for (int i = 0; i < 3; ++i) {
      std::uint32_t a = triangle[static_cast<std::size_t>(i)];
      std::uint32_t b = triangle[static_cast<std::size_t>((i + 1) % 3)];
      if (a > b) std::swap(a, b);
      ++uses[{a, b}];
    }
  }
  std::size_t count = 0;
  for (const auto& entry : uses) {
    if (entry.second == 1) ++count;
  }
  return count;
}

}  // namespace

TEST_CASE("extrude a square profile") {
  const Polygon2d square = square_profile();
  const Plane3 plane{geom::Vec3d{0.0, 0.0, 1.0}, 0.0};
  const Result<Brep> solid = extrude_profile(plane, square, 3.0, 1e-9);
  REQUIRE(solid.has_value());
  CHECK(solid.value().solid_count() == 1u);
  CHECK(solid.value().shell_count() == 1u);
  CHECK(solid.value().vertex_count() == 8u);
  CHECK(solid.value().edge_count() == 12u);
  CHECK(solid.value().face_count() == 6u);
  CHECK(volume_of(solid.value()) == doctest::Approx(12.0).epsilon(1e-12));
  REQUIRE(validate_topology(solid.value()).ok());

  const Result<TriMesh> mesh = tessellate(solid.value(), TessellationOptions{});
  REQUIRE(mesh.has_value());
  CHECK(boundary_edges(mesh.value()) == 0u);
  CHECK(!mesh.value().triangles.empty());
}

TEST_CASE("extrude rejects invalid options and height") {
  const Polygon2d square = square_profile();
  const Plane3 plane{geom::Vec3d{0.0, 0.0, 1.0}, 0.0};

  CHECK(extrude_profile(plane, square, 1.0, 0.0).error().code ==
        Errc::tolerance_violation);
  CHECK(extrude_profile(plane, square, 1.0, -1e-9).error().code ==
        Errc::tolerance_violation);
  CHECK(extrude_profile(plane, square, 1.0, std::numeric_limits<double>::quiet_NaN())
            .error().code == Errc::tolerance_violation);
  CHECK(extrude_profile(plane, square, 0.0, 1e-9).error().code == Errc::degenerate);
  CHECK(extrude_profile(plane, square, -1.0, 1e-9).error().code == Errc::degenerate);
  CHECK(extrude_profile(plane, square, std::numeric_limits<double>::infinity(), 1e-9)
            .error().code == Errc::non_finite);
}

TEST_CASE("extrude rejects a structurally invalid profile") {
  const Plane3 plane{geom::Vec3d{0.0, 0.0, 1.0}, 0.0};
  const Polygon2d empty;
  CHECK(extrude_profile(plane, empty, 1.0, 1e-9).error().code ==
        Errc::degenerate);

  Polygon2d duplicate = square_profile();
  duplicate.vertices.insert(duplicate.vertices.begin(), Point2d{0.0, 0.0});
  duplicate.ring_starts = {0u, 5u};
  CHECK(extrude_profile(plane, duplicate, 1.0, 1e-9).error().code ==
        Errc::degenerate);
}
