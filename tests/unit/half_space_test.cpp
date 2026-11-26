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

constexpr double kPi = 3.141592653589793238462643383279502884;

// The axis-aligned box [0, extent]^3 rotated by Rz(30 deg) * Rx(20 deg), then
// translated. The rotation is orthonormal, so the volume stays extent^3.
struct RotatedBox {
  std::vector<detail::HalfSpacePlane> planes;
  std::vector<Point3d> vertices;
  double volume = 0.0;
};

RotatedBox rotated_box(double extent, const Vec3d& translation) {
  const double a = 30.0 * kPi / 180.0;
  const double b = 20.0 * kPi / 180.0;
  const double ca = std::cos(a);
  const double sa = std::sin(a);
  const double cb = std::cos(b);
  const double sb = std::sin(b);
  const Vec3d axes[3] = {
      Vec3d{ca, sa, 0.0},
      Vec3d{-sa * cb, ca * cb, sb},
      Vec3d{sa * sb, -ca * sb, cb},
  };

  RotatedBox box;
  box.planes.reserve(6u);
  box.vertices.reserve(8u);
  for (const Vec3d& axis : axes) {
    for (const double sign : {1.0, -1.0}) {
      const Vec3d normal = axis * sign;
      const double base = sign > 0.0 ? extent : 0.0;
      box.planes.push_back({normal, base + dot(normal, translation)});
    }
  }
  for (int mask = 0; mask < 8; ++mask) {
    Vec3d corner = translation;
    for (int axis = 0; axis < 3; ++axis) {
      if ((mask & (1 << axis)) != 0) corner = corner + axes[axis] * extent;
    }
    box.vertices.push_back(Point3d{corner[0], corner[1], corner[2]});
  }
  box.volume = extent * extent * extent;
  return box;
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

TEST_CASE("half-space faces of a rotated box at planetary scale") {
  // The corners sit at ~1e12, so evaluating a plane equation rounds by ~1e-4.
  // That is far above the 1e-9 tolerance and only the per-plane relative slack
  // (1e-12 * scale) keeps every corner incident; a slack of plain `tolerance`
  // would drop the whole cell.
  const RotatedBox box = rotated_box(1e12, Vec3d{5e12, -3e12, 7e12});
  const Result<std::vector<Point3d>> vertices =
      detail::half_space_vertices(box.planes, 1e-9);
  REQUIRE(vertices.has_value());
  REQUIRE(vertices.value().size() == 8u);
  const Result<std::vector<AssembledFace>> faces =
      detail::half_space_faces(box.planes, vertices.value(), 1e-9);
  REQUIRE(faces.has_value());
  REQUIRE(faces.value().size() == 6u);
  for (const AssembledFace& face : faces.value()) {
    CHECK(face.outer.size() == 4u);
    CHECK(face.holes.empty());
    CHECK(face_ring_area_dot_normal(face) > 0.0);
  }
  // assemble_solid takes an absolute tolerance, and at 1e12 the ring planarity
  // residual is ~1e-4, so it needs 1e-3 rather than the 1e-9 used for the cell
  // size of the other boxes. The geometry is unchanged: distinct corners are
  // 1e12 apart.
  const Result<Brep> solid = assemble_solid(faces.value(), 1e-3);
  REQUIRE(solid.has_value());
  const Result<MassProperties> mass = mass_properties(solid.value());
  REQUIRE(mass.has_value());
  CHECK(std::abs(mass.value().volume - box.volume) <= 1e-9 * box.volume);
}

TEST_CASE("half-space faces ignore a far plane on a unit cell") {
  // A global slack of 1e-12 * max|offset| would be 1.0 here, because of the far
  // plane's 1e12 offset, and would then capture every vertex within one unit of
  // the unit cube's planes.
  auto planes = cube_planes(1.0);
  planes.push_back({{0.0, 0.0, 1.0}, 1e12});
  const Result<std::vector<Point3d>> vertices =
      detail::half_space_vertices(planes, 1e-9);
  REQUIRE(vertices.has_value());
  REQUIRE(vertices.value().size() == 8u);
  const Result<std::vector<AssembledFace>> faces =
      detail::half_space_faces(planes, vertices.value(), 1e-9);
  REQUIRE(faces.has_value());
  REQUIRE(faces.value().size() == 6u);
  for (const AssembledFace& face : faces.value()) {
    CHECK(face.outer.size() == 4u);
    CHECK(face_ring_area_dot_normal(face) > 0.0);
  }
  const Result<Brep> solid = assemble_solid(faces.value(), 1e-9);
  REQUIRE(solid.has_value());
  const Result<MassProperties> mass = mass_properties(solid.value());
  REQUIRE(mass.has_value());
  CHECK(std::abs(mass.value().volume - 1.0) <= 1e-9);
}

TEST_CASE("half-space faces skip a plane that supports only an edge") {
  // x - y = 2 touches the box [0,2] x [0,3] x [0,4] only along the edge x = 2,
  // y = 0, so it has exactly two incident vertices and must not become a face.
  auto planes = unit_box_planes();
  planes.push_back({{1.0, -1.0, 0.0}, 2.0});
  const Result<std::vector<Point3d>> vertices =
      detail::half_space_vertices(planes, 1e-9);
  REQUIRE(vertices.has_value());
  const Result<std::vector<AssembledFace>> faces =
      detail::half_space_faces(planes, vertices.value(), 1e-9);
  REQUIRE(faces.has_value());
  REQUIRE(faces.value().size() == 6u);
  for (const AssembledFace& face : faces.value()) {
    CHECK(face.outer.size() == 4u);
    CHECK(face_ring_area_dot_normal(face) > 0.0);
  }
  const Result<Brep> solid = assemble_solid(faces.value(), 1e-9);
  REQUIRE(solid.has_value());
  const Result<MassProperties> mass = mass_properties(solid.value());
  REQUIRE(mass.has_value());
  CHECK(std::abs(mass.value().volume - 24.0) <= 1e-9);
}

TEST_CASE("half-space faces keep one face per plane for duplicate planes") {
  // Coincident planes are not deduplicated: the duplicate supports the same
  // four vertices and is emitted as a second, coincident face. Callers must
  // deduplicate coincident planes before calling; two coincident faces cannot
  // be assembled into one solid, so only the ring shape is checked here.
  auto planes = unit_box_planes();
  planes.push_back({{1.0, 0.0, 0.0}, 2.0});
  const Result<std::vector<Point3d>> vertices =
      detail::half_space_vertices(planes, 1e-9);
  REQUIRE(vertices.has_value());
  const Result<std::vector<AssembledFace>> faces =
      detail::half_space_faces(planes, vertices.value(), 1e-9);
  REQUIRE(faces.has_value());
  REQUIRE(faces.value().size() == 7u);
  for (const AssembledFace& face : faces.value()) {
    CHECK(face.outer.size() == 4u);
    CHECK(face_ring_area_dot_normal(face) > 0.0);
  }
}

TEST_CASE("half-space solid of a right tetrahedron has the analytic volume") {
  const double s = std::sqrt(3.0);
  // Scaling every equation of the tetrahedron below by s keeps the same
  // half-spaces; using s here (rather than an unused variable) keeps the test
  // warning-clean under /W4 /WX.
  const std::array<detail::HalfSpacePlane, 4> planes{{
      {{-s, 0.0, 0.0}, 0.0},
      {{0.0, -s, 0.0}, 0.0},
      {{0.0, 0.0, -s}, 0.0},
      {{s, s, s}, 2.0 * s},
  }};
  const Result<Brep> solid = detail::half_space_solid(planes, 1e-9);
  REQUIRE(solid.has_value());
  CHECK(validate_topology(solid.value()).ok());
  const Result<MassProperties> mass = mass_properties(solid.value());
  REQUIRE(mass.has_value());
  CHECK(std::abs(mass.value().volume - 4.0 / 3.0) <= 1e-9);
}

TEST_CASE("half-space solid rejects an empty cell") {
  const std::array<detail::HalfSpacePlane, 2> planes{{
      {{1.0, 0.0, 0.0}, 0.0},
      {{-1.0, 0.0, 0.0}, -1.0},
  }};
  const Result<Brep> solid = detail::half_space_solid(planes, 1e-9);
  REQUIRE_FALSE(solid.has_value());
  CHECK(solid.error().code == Errc::degenerate);
  CHECK(std::string_view(solid.error().message) ==
        "half_space: planes have no common point");
}

TEST_CASE("half-space solid of a box is a closed six-face box") {
  const auto planes = unit_box_planes();
  const Result<Brep> solid = detail::half_space_solid(planes, 1e-9);
  REQUIRE(solid.has_value());
  CHECK(validate_topology(solid.value()).ok());
  CHECK(solid.value().shell_count() == 1u);
  CHECK(solid.value().face_count() == 6u);
  CHECK(solid.value().edge_count() == 12u);
  CHECK(solid.value().vertex_count() == 8u);
  const Result<MassProperties> mass = mass_properties(solid.value());
  REQUIRE(mass.has_value());
  CHECK(std::abs(mass.value().volume - 24.0) <= 1e-9);
}

TEST_CASE("half-space solid is unchanged by a common plane scale") {
  const auto planes = tetrahedron_planes(3.0);
  const Result<Brep> solid = detail::half_space_solid(planes, 1e-9);
  REQUIRE(solid.has_value());
  CHECK(validate_topology(solid.value()).ok());
  const Result<MassProperties> mass = mass_properties(solid.value());
  REQUIRE(mass.has_value());
  CHECK(std::abs(mass.value().volume - 4.0 / 3.0) <= 1e-9);
}

TEST_CASE("half-space solid is deterministic for a fixed plane order") {
  const auto planes = unit_box_planes();
  const Result<Brep> first = detail::half_space_solid(planes, 1e-9);
  const Result<Brep> second = detail::half_space_solid(planes, 1e-9);
  REQUIRE(first.has_value());
  REQUIRE(second.has_value());
  CHECK(first.value().shell_count() == second.value().shell_count());
  CHECK(first.value().face_count() == second.value().face_count());
  CHECK(first.value().edge_count() == second.value().edge_count());
  CHECK(first.value().vertex_count() == second.value().vertex_count());
  const Result<MassProperties> first_mass = mass_properties(first.value());
  const Result<MassProperties> second_mass = mass_properties(second.value());
  REQUIRE(first_mass.has_value());
  REQUIRE(second_mass.has_value());
  CHECK(first_mass.value().volume == second_mass.value().volume);
}

TEST_CASE("half-space solid propagates an assembly failure") {
  // The coincident duplicate plane is not removed by half_space_faces, so the
  // assembler sees a non-manifold edge. Its error must reach the caller
  // unchanged, so pin the exact message as well as the code: a rewrite such as
  // "half_space: non-manifold cell" must not pass.
  auto planes = unit_box_planes();
  planes.push_back({{1.0, 0.0, 0.0}, 2.0});
  const Result<Brep> solid = detail::half_space_solid(planes, 1e-9);
  REQUIRE_FALSE(solid.has_value());
  CHECK(solid.error().code == Errc::invalid_topology);
  CHECK(std::string_view(solid.error().message) ==
        "assemble_solid: non-manifold edge");
}
