#include <geom/brep/mesh_inside.hpp>
#include <geom/brep/primitives.hpp>

#include <doctest/doctest.h>

#include <array>
#include <cmath>
#include <cstdint>
#include <vector>

using geom::v1::point_in_tessellation;
using geom::v1::Point3d;
using geom::v1::TriMesh;

namespace {

constexpr double kTolerance = 1e-9;

TriMesh cube_mesh() {
  TriMesh mesh;
  mesh.positions = {
      Point3d{0.0, 0.0, 0.0}, Point3d{1.0, 0.0, 0.0}, Point3d{1.0, 1.0, 0.0},
      Point3d{0.0, 1.0, 0.0}, Point3d{0.0, 0.0, 1.0}, Point3d{1.0, 0.0, 1.0},
      Point3d{1.0, 1.0, 1.0}, Point3d{0.0, 1.0, 1.0},
  };
  mesh.triangles = {
      {0u, 2u, 1u}, {0u, 3u, 2u},  // z = 0
      {4u, 5u, 6u}, {4u, 6u, 7u},  // z = 1
      {0u, 1u, 5u}, {0u, 5u, 4u},  // y = 0
      {3u, 7u, 6u}, {3u, 6u, 2u},  // y = 1, diagonal 3-6
      {0u, 4u, 7u}, {0u, 7u, 3u},  // x = 0
      {1u, 2u, 6u}, {1u, 6u, 5u},  // x = 1, diagonal 1-6
  };
  mesh.triangle_face.assign(mesh.triangles.size(), 0u);
  return mesh;
}

}  // namespace

TEST_CASE("point_in_tessellation classifies a tessellated cube") {
  const auto box = geom::make_box(Point3d{0.0, 0.0, 0.0}, Point3d{1.0, 1.0, 1.0}, kTolerance);
  REQUIRE(box.has_value());
  geom::TessellationOptions options;
  options.deflection = 1e-3;
  const auto mesh = geom::tessellate(box.value(), options);
  REQUIRE(mesh.has_value());

  const auto center = point_in_tessellation(mesh.value(), Point3d{0.5, 0.5, 0.5}, kTolerance);
  REQUIRE(center.has_value());
  CHECK(center.value());

  const auto outside = point_in_tessellation(mesh.value(), Point3d{2.0, 0.0, 0.0}, kTolerance);
  REQUIRE(outside.has_value());
  CHECK_FALSE(outside.value());

  const auto surface =
      point_in_tessellation(mesh.value(), Point3d{0.5, 0.5, 1.0 + 1e-12}, kTolerance);
  REQUIRE_FALSE(surface.has_value());
  CHECK(surface.error().code == geom::v1::Errc::numeric_failure);
}

TEST_CASE("point_in_tessellation falls back when a ray grazes a shared diagonal") {
  const TriMesh mesh = cube_mesh();

  // The +x ray from this point hits the x = 1 diagonal at (1, 0.3, 0.3), so the
  // classifier must fall through to another direction and still answer inside.
  const auto inside = point_in_tessellation(mesh, Point3d{0.4, 0.3, 0.3}, kTolerance);
  REQUIRE(inside.has_value());
  CHECK(inside.value());

  // An outside point whose +x ray enters through the x = 0 diagonal is also
  // ambiguous on that ray and must fall through to another direction.
  const auto outside = point_in_tessellation(mesh, Point3d{-0.5, 0.3, 0.3}, kTolerance);
  REQUIRE(outside.has_value());
  CHECK_FALSE(outside.value());
}

TEST_CASE("mesh_inside_triangle_distance2 survives a degenerate triangle") {
  using geom::v1::detail::mesh_inside_detail::mesh_inside_triangle_distance2;
  const Point3d a{0.0, 0.0, 0.0};
  const Point3d b{0.0, 0.0, 0.0};
  const Point3d c{1.0, 0.0, 0.0};
  const double d2 = mesh_inside_triangle_distance2(Point3d{0.5, 1.0, 0.0}, a, b, c);
  CHECK(std::isfinite(d2));
  CHECK(d2 == doctest::Approx(1.0));
}

TEST_CASE("point_in_tessellation rejects degenerate inputs") {
  const TriMesh mesh = cube_mesh();

  const auto nan_point =
      point_in_tessellation(mesh, Point3d{std::nan(""), 0.5, 0.5}, kTolerance);
  REQUIRE_FALSE(nan_point.has_value());
  CHECK(nan_point.error().code == geom::v1::Errc::non_finite);

  const auto bad_tolerance = point_in_tessellation(mesh, Point3d{0.5, 0.5, 0.5}, 0.0);
  REQUIRE_FALSE(bad_tolerance.has_value());
  CHECK(bad_tolerance.error().code == geom::v1::Errc::tolerance_violation);

  const TriMesh empty;
  const auto empty_mesh = point_in_tessellation(empty, Point3d{0.5, 0.5, 0.5}, kTolerance);
  REQUIRE_FALSE(empty_mesh.has_value());
  CHECK(empty_mesh.error().code == geom::v1::Errc::degenerate);
}
