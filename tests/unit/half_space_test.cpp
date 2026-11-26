#include <doctest/doctest.h>

#include <array>
#include <cmath>
#include <cstddef>
#include <limits>
#include <span>
#include <string_view>
#include <vector>

#include <geom/brep/assemble.hpp>
#include <geom/brep/detail/half_space.hpp>
#include <geom/brep/mass_properties.hpp>
#include <geom/brep/validate.hpp>

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

// Four half-spaces of the tetrahedron with corners (2,0,0), (0,2,0), (0,0,2)
// and the origin. `factor` scales every equation, which must leave the cell and
// the emitted faces unchanged.
std::array<detail::HalfSpacePlane, 4> tetrahedron_planes(double factor = 1.0) {
  return {{
      {{-factor, 0.0, 0.0}, 0.0},
      {{0.0, -factor, 0.0}, 0.0},
      {{0.0, 0.0, -factor}, 0.0},
      {{factor, factor, factor}, 2.0 * factor},
  }};
}

// Twice the signed area vector of an open ring, dotted with the face normal.
// It is positive exactly when the ring winds counter-clockwise about the
// outward normal.
double face_ring_area_dot_normal(const AssembledFace& face) {
  Vec3d sum{};
  for (std::size_t i = 0; i < face.outer.size(); ++i) {
    const Point3d& a = face.outer[i];
    const Point3d& b = face.outer[(i + 1u) % face.outer.size()];
    sum = sum + cross(Vec3d{a[0], a[1], a[2]}, Vec3d{b[0], b[1], b[2]});
  }
  return 0.5 * dot(sum, face.normal);
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

TEST_CASE("half-space faces of a box are six outward quads") {
  const auto planes = unit_box_planes();
  const Result<std::vector<Point3d>> vertices =
      detail::half_space_vertices(planes, 1e-9);
  REQUIRE(vertices.has_value());
  const Result<std::vector<AssembledFace>> faces =
      detail::half_space_faces(planes, vertices.value(), 1e-9);
  REQUIRE(faces.has_value());
  CHECK(faces.value().size() == 6u);
  for (const AssembledFace& face : faces.value()) {
    CHECK(face.outer.size() == 4u);
    CHECK(face.holes.empty());
    CHECK(face_ring_area_dot_normal(face) > 0.0);
  }
  const Result<Brep> solid =
      assemble_solid(faces.value(), 1e-9);
  REQUIRE(solid.has_value());
  const Result<MassProperties> mass = mass_properties(solid.value());
  REQUIRE(mass.has_value());
  CHECK(std::abs(mass.value().volume - 24.0) <= 1e-9);
}

TEST_CASE("half-space faces of a tetrahedron are four outward triangles") {
  const auto planes = tetrahedron_planes();
  const Result<std::vector<Point3d>> vertices =
      detail::half_space_vertices(planes, 1e-9);
  REQUIRE(vertices.has_value());
  REQUIRE(vertices.value().size() == 4u);
  const Result<std::vector<AssembledFace>> faces =
      detail::half_space_faces(planes, vertices.value(), 1e-9);
  REQUIRE(faces.has_value());
  REQUIRE(faces.value().size() == 4u);
  for (const AssembledFace& face : faces.value()) {
    CHECK(face.outer.size() == 3u);
    CHECK(face.holes.empty());
    CHECK(face_ring_area_dot_normal(face) > 0.0);
  }
  const Result<Brep> solid = assemble_solid(faces.value(), 1e-9);
  REQUIRE(solid.has_value());
  const Result<MassProperties> mass = mass_properties(solid.value());
  REQUIRE(mass.has_value());
  CHECK(std::abs(mass.value().volume - 4.0 / 3.0) <= 1e-9);
}

TEST_CASE("half-space faces ignore a common plane scale") {
  const auto planes = tetrahedron_planes(3.0);
  const Result<std::vector<Point3d>> vertices =
      detail::half_space_vertices(planes, 1e-9);
  REQUIRE(vertices.has_value());
  REQUIRE(vertices.value().size() == 4u);
  const Result<std::vector<AssembledFace>> faces =
      detail::half_space_faces(planes, vertices.value(), 1e-9);
  REQUIRE(faces.has_value());
  REQUIRE(faces.value().size() == 4u);
  for (const AssembledFace& face : faces.value()) {
    CHECK(face.outer.size() == 3u);
    CHECK(face_ring_area_dot_normal(face) > 0.0);
  }
  const Result<Brep> solid = assemble_solid(faces.value(), 1e-9);
  REQUIRE(solid.has_value());
  const Result<MassProperties> mass = mass_properties(solid.value());
  REQUIRE(mass.has_value());
  CHECK(std::abs(mass.value().volume - 4.0 / 3.0) <= 1e-9);
}

TEST_CASE("half-space faces skip a redundant far plane") {
  auto planes = unit_box_planes();
  planes.push_back({{0.0, 0.0, 1.0}, 1e12});
  const Result<std::vector<Point3d>> vertices =
      detail::half_space_vertices(planes, 1e-9);
  REQUIRE(vertices.has_value());
  const Result<std::vector<AssembledFace>> faces =
      detail::half_space_faces(planes, vertices.value(), 1e-9);
  REQUIRE(faces.has_value());
  CHECK(faces.value().size() == 6u);
  for (const AssembledFace& face : faces.value()) {
    CHECK(face.outer.size() == 4u);
    CHECK(face_ring_area_dot_normal(face) > 0.0);
  }
}

TEST_CASE("half-space faces reject an empty vertex set") {
  const auto planes = unit_box_planes();
  const std::vector<Point3d> vertices;
  const Result<std::vector<AssembledFace>> faces =
      detail::half_space_faces(planes, vertices, 1e-9);
  REQUIRE_FALSE(faces.has_value());
  CHECK(faces.error().code == Errc::degenerate);
  CHECK(std::string_view(faces.error().message) ==
        "half_space: fewer than three planes bound a cell");
}

TEST_CASE("half-space faces reject invalid tolerances") {
  const auto planes = unit_box_planes();
  const Result<std::vector<Point3d>> vertices =
      detail::half_space_vertices(planes, 1e-9);
  REQUIRE(vertices.has_value());
  const std::array<double, 4> tolerances{
      0.0, -1.0, std::numeric_limits<double>::quiet_NaN(),
      std::numeric_limits<double>::infinity()};
  for (const double tolerance : tolerances) {
    const Result<std::vector<AssembledFace>> faces =
        detail::half_space_faces(planes, vertices.value(), tolerance);
    REQUIRE_FALSE(faces.has_value());
    CHECK(faces.error().code == Errc::tolerance_violation);
  }
}

TEST_CASE("half-space faces reject invalid planes") {
  const std::vector<Point3d> vertices{
      Point3d{0.0, 0.0, 0.0}, Point3d{1.0, 0.0, 0.0},
      Point3d{0.0, 1.0, 0.0}, Point3d{0.0, 0.0, 1.0}};

  SUBCASE("zero normal") {
    const std::array<detail::HalfSpacePlane, 1> planes{{
        {{0.0, 0.0, 0.0}, 1.0},
    }};
    const Result<std::vector<AssembledFace>> faces =
        detail::half_space_faces(planes, vertices, 1e-9);
    REQUIRE_FALSE(faces.has_value());
    CHECK(faces.error().code == Errc::degenerate);
    CHECK(std::string_view(faces.error().message) ==
          "half_space: plane normal is zero or non-finite");
  }

  SUBCASE("non-finite offset") {
    const std::array<detail::HalfSpacePlane, 1> planes{{
        {{1.0, 0.0, 0.0}, std::nan("")},
    }};
    const Result<std::vector<AssembledFace>> faces =
        detail::half_space_faces(planes, vertices, 1e-9);
    REQUIRE_FALSE(faces.has_value());
    CHECK(faces.error().code == Errc::non_finite);
    CHECK(std::string_view(faces.error().message) ==
          "half_space: plane offset is not finite");
  }
}
