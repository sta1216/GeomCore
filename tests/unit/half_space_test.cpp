#include <doctest/doctest.h>

#include <array>
#include <cmath>
#include <cstddef>
#include <limits>
#include <span>
#include <string_view>
#include <vector>

#include <geom/brep/detail/half_space.hpp>

using namespace geom;

namespace {

std::vector<detail::HalfSpacePlane> unit_box_planes() {
  return {
      {{1.0, 0.0, 0.0}, 2.0},
      {{-1.0, 0.0, 0.0}, 0.0},
      {{0.0, 1.0, 0.0}, 3.0},
      {{0.0, -1.0, 0.0}, 0.0},
      {{0.0, 0.0, 1.0}, 4.0},
      {{0.0, 0.0, -1.0}, 0.0},
  };
}

// Six half-spaces of the cube [0, extent]^3.
std::vector<detail::HalfSpacePlane> cube_planes(double extent) {
  return {
      {{1.0, 0.0, 0.0}, extent},
      {{-1.0, 0.0, 0.0}, 0.0},
      {{0.0, 1.0, 0.0}, extent},
      {{0.0, -1.0, 0.0}, 0.0},
      {{0.0, 0.0, 1.0}, extent},
      {{0.0, 0.0, -1.0}, 0.0},
  };
}

}  // namespace

TEST_CASE("half-space vertices of a box are its eight corners") {
  const auto planes = unit_box_planes();
  const Result<std::vector<Point3d>> vertices =
      detail::half_space_vertices(planes, 1e-9);
  REQUIRE(vertices.has_value());
  CHECK(vertices.value().size() == 8u);
  for (const Point3d& p : vertices.value()) {
    CHECK((p[0] == 0.0 || p[0] == 2.0));
    CHECK((p[1] == 0.0 || p[1] == 3.0));
    CHECK((p[2] == 0.0 || p[2] == 4.0));
  }
}

TEST_CASE("half-space vertex order is fixed for a fixed plane order") {
  const auto planes = unit_box_planes();
  const Result<std::vector<Point3d>> first =
      detail::half_space_vertices(planes, 1e-9);
  REQUIRE(first.has_value());
  const Result<std::vector<Point3d>> second =
      detail::half_space_vertices(planes, 1e-9);
  REQUIRE(second.has_value());
  REQUIRE(second.value().size() == first.value().size());
  for (std::size_t i = 0; i < first.value().size(); ++i) {
    CHECK(first.value()[i][0] == second.value()[i][0]);
    CHECK(first.value()[i][1] == second.value()[i][1]);
    CHECK(first.value()[i][2] == second.value()[i][2]);
  }
}

TEST_CASE("half-space vertices ignore a redundant plane") {
  auto planes = unit_box_planes();
  planes.push_back({{0.0, 0.0, 1.0}, 10.0});
  const Result<std::vector<Point3d>> vertices =
      detail::half_space_vertices(planes, 1e-9);
  REQUIRE(vertices.has_value());
  CHECK(vertices.value().size() == 8u);
}

TEST_CASE("half-space vertices ignore a far redundant plane") {
  auto planes = unit_box_planes();
  planes.push_back({{0.0, 0.0, 1.0}, 1e12});
  const Result<std::vector<Point3d>> vertices =
      detail::half_space_vertices(planes, 1e-9);
  REQUIRE(vertices.has_value());
  CHECK(vertices.value().size() == 8u);
}

TEST_CASE("half-space vertices deduplicate coincident planes") {
  auto planes = unit_box_planes();
  planes.push_back({{1.0, 0.0, 0.0}, 2.0});
  const Result<std::vector<Point3d>> vertices =
      detail::half_space_vertices(planes, 1e-9);
  REQUIRE(vertices.has_value());
  CHECK(vertices.value().size() == 8u);
}

TEST_CASE("half-space vertices of a huge cube are its eight corners") {
  const auto planes = cube_planes(1e150);
  const Result<std::vector<Point3d>> vertices =
      detail::half_space_vertices(planes, 1e-9);
  REQUIRE(vertices.has_value());
  CHECK(vertices.value().size() == 8u);
}

TEST_CASE("half-space vertices of a tiny cube are its eight corners") {
  const auto planes = cube_planes(1e-150);
  const Result<std::vector<Point3d>> vertices =
      detail::half_space_vertices(planes, 1e-155);
  REQUIRE(vertices.has_value());
  CHECK(vertices.value().size() == 8u);
}

TEST_CASE("half-space vertices skip near-parallel plane triples") {
  // Both extra planes pass through the corner (2, 2, 2) and lean by ~1e-9, so
  // their triples with the box planes are nearly singular. The cell is still
  // the untouched cube; an ill-conditioned solve must not add a vertex.
  auto planes = cube_planes(2.0);
  planes.push_back({{1.0, 1e-9, 1e-9}, 2.0 + 4e-9});
  planes.push_back({{1.0, 2e-9, 3e-9}, 2.0 + 10e-9});
  const Result<std::vector<Point3d>> vertices =
      detail::half_space_vertices(planes, 1e-9);
  REQUIRE(vertices.has_value());
  CHECK(vertices.value().size() == 8u);
}

TEST_CASE("half-space vertices report an empty four-plane cell") {
  const std::array<detail::HalfSpacePlane, 4> planes{{
      {{1.0, 0.0, 0.0}, 0.0},
      {{0.0, 1.0, 0.0}, 0.0},
      {{0.0, 0.0, 1.0}, 0.0},
      {{-1.0, -1.0, -1.0}, -1.0},
  }};
  const Result<std::vector<Point3d>> vertices =
      detail::half_space_vertices(planes, 1e-9);
  REQUIRE(vertices.has_value());
  CHECK(vertices.value().empty());
}

TEST_CASE("half-space vertices report an empty intersection") {
  const std::array<detail::HalfSpacePlane, 2> planes{{
      {{1.0, 0.0, 0.0}, 0.0},
      {{-1.0, 0.0, 0.0}, -1.0},
  }};
  const Result<std::vector<Point3d>> vertices =
      detail::half_space_vertices(planes, 1e-9);
  REQUIRE(vertices.has_value());
  CHECK(vertices.value().empty());
}

TEST_CASE("half-space vertices reject invalid tolerances") {
  const auto planes = unit_box_planes();
  const std::array<double, 4> tolerances{
      0.0, -1.0, std::numeric_limits<double>::quiet_NaN(),
      std::numeric_limits<double>::infinity()};
  for (const double tolerance : tolerances) {
    const Result<std::vector<Point3d>> vertices =
        detail::half_space_vertices(planes, tolerance);
    REQUIRE_FALSE(vertices.has_value());
    CHECK(vertices.error().code == Errc::tolerance_violation);
  }
}

TEST_CASE("half-space vertices reject a non-finite plane offset") {
  const std::array<detail::HalfSpacePlane, 1> planes{{
      {{1.0, 0.0, 0.0}, std::nan("")},
  }};
  const Result<std::vector<Point3d>> vertices =
      detail::half_space_vertices(planes, 1e-9);
  REQUIRE_FALSE(vertices.has_value());
  CHECK(vertices.error().code == Errc::non_finite);
}

TEST_CASE("half-space vertices reject a zero normal") {
  const std::array<detail::HalfSpacePlane, 1> planes{{
      {{0.0, 0.0, 0.0}, 1.0},
  }};
  const Result<std::vector<Point3d>> vertices =
      detail::half_space_vertices(planes, 1e-9);
  REQUIRE_FALSE(vertices.has_value());
  CHECK(vertices.error().code == Errc::degenerate);
}
