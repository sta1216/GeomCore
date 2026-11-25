#include <doctest/doctest.h>

#include <array>
#include <cmath>
#include <cstdint>
#include <limits>
#include <span>
#include <string>
#include <vector>

#include <geom/brep/boolean_solid.hpp>
#include <geom/brep/chamfer.hpp>
#include <geom/brep/mass_properties.hpp>
#include <geom/brep/primitives.hpp>
#include <geom/brep/validate.hpp>
#include <geom/core/error.hpp>

using namespace geom;

namespace {

Brep make_box_solid(Point3d lo, Point3d hi) {
  const Result<Brep> box = make_box(lo, hi, 1e-9);
  REQUIRE(box.has_value());
  return box.value();
}

double volume_of(const Brep& brep) {
  const Result<MassProperties> mass = mass_properties(brep);
  REQUIRE(mass.has_value());
  return mass.value().volume;
}

// Finds the straight edge whose endpoints are the two given corners.
EdgeId edge_between(const Brep& brep, Point3d a, Point3d b) {
  for (const EdgeId id : brep.edge_ids()) {
    const Point3d p0 = brep.vertex(brep.edge(id).v0).point;
    const Point3d p1 = brep.vertex(brep.edge(id).v1).point;
    const auto same = [](Point3d p, Point3d q) {
      return std::abs(p[0] - q[0]) + std::abs(p[1] - q[1]) +
                 std::abs(p[2] - q[2]) <
             1e-12;
    };
    if ((same(p0, a) && same(p1, b)) || (same(p0, b) && same(p1, a))) return id;
  }
  REQUIRE_MESSAGE(false, "edge not found");
  return {};
}

}  // namespace

TEST_CASE("chamfer one box edge removes a triangular prism") {
  const Brep box = make_box_solid(Point3d{0.0, 0.0, 0.0}, Point3d{2.0, 3.0, 4.0});
  const EdgeId edge =
      edge_between(box, Point3d{0.0, 0.0, 0.0}, Point3d{0.0, 0.0, 4.0});
  const std::array<EdgeId, 1> edges{edge};

  ChamferOptions options;
  options.distance = 0.25;
  const Result<ChamferResult> result = chamfer_edges(box, edges, options);
  REQUIRE(result.has_value());
  CHECK(result.value().chamfered_edges == 1u);
  CHECK(result.value().brep.solid_count() == 1u);
  CHECK(result.value().brep.shell_count() == 1u);
  CHECK(validate_topology(result.value().brep).ok());
  // Removed cross-section is d^2/2 = 0.03125 over the 4-long edge.
  CHECK(std::abs(volume_of(result.value().brep) - 23.875) <= 1e-9);
  // The input document is untouched.
  CHECK(box.face_count() == 6u);
  CHECK(std::abs(volume_of(box) - 24.0) <= 1e-12);
}

TEST_CASE("chamfer rejects a non-positive distance") {
  const Brep box = make_box_solid(Point3d{0.0, 0.0, 0.0}, Point3d{1.0, 1.0, 1.0});
  const std::array<EdgeId, 1> edges{box.edge_ids().front()};
  ChamferOptions options;
  options.distance = 0.0;
  const Result<ChamferResult> result = chamfer_edges(box, edges, options);
  REQUIRE_FALSE(result.has_value());
  CHECK(result.error().code == Errc::degenerate);
}

// Guard: the shared planar-convex helper Task 2/Task 3 build on must describe a
// box as six outward planes with unit normals, without touching the input.
TEST_CASE("chamfer helper describes a box as six planes") {
  const Brep box = make_box_solid(Point3d{0.0, 0.0, 0.0}, Point3d{2.0, 3.0, 4.0});
  const Result<detail::BlendPlanarSolid> view =
      detail::blend_planar_convex_solid(box, 1e-9, "chamfer_edges");
  REQUIRE(view.has_value());
  CHECK(view.value().solid == box.solid_ids().front());
  CHECK(view.value().faces.size() == 6u);
  CHECK(view.value().scale == 4.0);
  for (const Vec3d outward : view.value().outward) {
    CHECK(std::abs(norm2(outward) - 1.0) <= 1e-12);
  }
}
