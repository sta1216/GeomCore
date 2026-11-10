#include <array>
#include <cmath>
#include <cstdint>
#include <limits>
#include <vector>

#include <doctest/doctest.h>
#include <geom/brep/detail/ssi_patch.hpp>
#include <geom/primitives/box.hpp>
#include <geom/primitives/point.hpp>
#include <geom/primitives/vec.hpp>
#include <geom/surfaces/analytic.hpp>
#include <geom/surfaces/nurbs_surface.hpp>
#include <geom/surfaces/plane.hpp>

using namespace geom::v1;

TEST_CASE("ssi: a planar patch tessellates into four triangles") {
  const Plane3 plane{{0.0, 0.0, 1.0}, 0.0};
  const Result<SsiPatch> patch = make_ssi_patch(SurfaceGeometry{plane}, {-1.0, 1.0}, {-1.0, 1.0});
  REQUIRE(patch.has_value());
  const Result<SsiPatchMesh> mesh = tessellate_ssi_patch(patch.value(), 1e-9, 12);
  REQUIRE(mesh.has_value());
  CHECK(mesh.value().triangles.size() == 4);
  CHECK(mesh.value().points.size() == 5);
  for (const Point3d& p : mesh.value().points) CHECK(std::fabs(p[2]) < 1e-15);
}

TEST_CASE("ssi: a sphere patch obeys the chord tolerance") {
  const Sphere3 sphere{Point3d{0, 0, 0}, Vec3d{0, 0, 1}, Vec3d{1, 0, 0}, 1.0};
  const Result<SsiPatch> patch = make_ssi_patch(SurfaceGeometry{sphere}, {0.0, 1.0}, {-0.5, 0.5});
  REQUIRE(patch.has_value());
  const Result<SsiPatchMesh> mesh = tessellate_ssi_patch(patch.value(), 1e-3, 12);
  REQUIRE(mesh.has_value());
  CHECK(mesh.value().triangles.size() > 8);
  for (const Point3d& p : mesh.value().points) {
    CHECK(std::fabs(std::sqrt(p[0] * p[0] + p[1] * p[1] + p[2] * p[2]) - 1.0) < 1e-12);
  }
  for (const auto& tri : mesh.value().triangles) {
    const Point3d& a = mesh.value().points[tri[0]];
    const Point3d& b = mesh.value().points[tri[1]];
    const Point3d& c = mesh.value().points[tri[2]];
    const Vec3d normal = cross(b - a, c - a);
    CHECK(norm2(normal) > 0.0);
  }
}

TEST_CASE("ssi: patch bounds contain every tessellated point") {
  const Torus3 torus{Point3d{0, 0, 0}, Vec3d{0, 0, 1}, Vec3d{1, 0, 0}, 3.0, 1.0};
  const Result<SsiPatch> patch = make_ssi_patch(SurfaceGeometry{torus}, {0.0, 2.0}, {0.0, 2.0});
  REQUIRE(patch.has_value());
  const Result<SsiPatchMesh> mesh = tessellate_ssi_patch(patch.value(), 1e-2, 12);
  REQUIRE(mesh.has_value());
  const Box3d bounds = patch_bounds(patch.value(), 1e-2);
  CHECK_FALSE(bounds.is_empty());
  for (const Point3d& p : mesh.value().points) {
    for (int i = 0; i < 3; ++i) {
      CHECK(p[i] >= bounds.lo[i] - 1e-12);
      CHECK(p[i] <= bounds.hi[i] + 1e-12);
    }
  }
}

TEST_CASE("ssi: tessellation rejects invalid tolerances") {
  const Plane3 plane{};
  const Result<SsiPatch> patch = make_ssi_patch(SurfaceGeometry{plane}, {-1.0, 1.0}, {-1.0, 1.0});
  REQUIRE(patch.has_value());
  CHECK_FALSE(tessellate_ssi_patch(patch.value(), 0.0, 12).has_value());
  CHECK_FALSE(tessellate_ssi_patch(patch.value(), -1.0, 12).has_value());
  CHECK_FALSE(tessellate_ssi_patch(patch.value(), 1e-9, 0).has_value());
}

TEST_CASE("ssi: tessellation rejects non-finite tolerance, bad depth and bad domains") {
  const Plane3 plane{};
  const Result<SsiPatch> patch = make_ssi_patch(SurfaceGeometry{plane}, {-1.0, 1.0}, {-1.0, 1.0});
  REQUIRE(patch.has_value());

  const Result<SsiPatchMesh> infinite =
      tessellate_ssi_patch(patch.value(), std::numeric_limits<double>::infinity(), 12);
  REQUIRE_FALSE(infinite.has_value());
  CHECK(infinite.error().code == Errc::tolerance_violation);

  const Result<SsiPatchMesh> negative_depth = tessellate_ssi_patch(patch.value(), 1e-9, -1);
  REQUIRE_FALSE(negative_depth.has_value());
  CHECK(negative_depth.error().code == Errc::index_out_of_range);

  SsiPatch invalid_surface;
  invalid_surface.surface =
      SurfaceGeometry{Cylinder3{Point3d{0, 0, 0}, Vec3d{0, 0, 1}, Vec3d{1, 0, 0}, 0.0}};
  invalid_surface.u = {0.0, 1.0};
  invalid_surface.v = {0.0, 1.0};
  const Result<SsiPatchMesh> bad_surface = tessellate_ssi_patch(invalid_surface, 1e-9, 12);
  REQUIRE_FALSE(bad_surface.has_value());
  CHECK(bad_surface.error().code == Errc::degenerate);

  SsiPatch inverted = patch.value();
  inverted.u = {1.0, -1.0};
  const Result<SsiPatchMesh> bad_domain = tessellate_ssi_patch(inverted, 1e-9, 12);
  REQUIRE_FALSE(bad_domain.has_value());
  CHECK(bad_domain.error().code == Errc::degenerate);
}

TEST_CASE("ssi: tessellation is deterministic") {
  const Sphere3 sphere{Point3d{0, 0, 0}, Vec3d{0, 0, 1}, Vec3d{1, 0, 0}, 1.0};
  const Result<SsiPatch> patch =
      make_ssi_patch(SurfaceGeometry{sphere}, {0.0, 1.5}, {-0.4, 0.7});
  REQUIRE(patch.has_value());
  const Result<SsiPatchMesh> first = tessellate_ssi_patch(patch.value(), 1e-2, 12);
  const Result<SsiPatchMesh> second = tessellate_ssi_patch(patch.value(), 1e-2, 12);
  REQUIRE(first.has_value());
  REQUIRE(second.has_value());
  REQUIRE(first.value().points.size() == second.value().points.size());
  REQUIRE(first.value().normals.size() == second.value().normals.size());
  REQUIRE(first.value().uv.size() == second.value().uv.size());
  REQUIRE(first.value().triangles.size() == second.value().triangles.size());
  for (std::size_t i = 0; i < first.value().points.size(); ++i) {
    for (int axis = 0; axis < 3; ++axis) {
      CHECK(first.value().points[i][axis] == second.value().points[i][axis]);
      CHECK(first.value().normals[i][axis] == second.value().normals[i][axis]);
    }
    CHECK(first.value().uv[i][0] == second.value().uv[i][0]);
    CHECK(first.value().uv[i][1] == second.value().uv[i][1]);
  }
  for (std::size_t i = 0; i < first.value().triangles.size(); ++i) {
    for (int k = 0; k < 3; ++k) {
      CHECK(first.value().triangles[i][k] == second.value().triangles[i][k]);
    }
  }
}

TEST_CASE("ssi: patch normals are unit and follow the sphere") {
  const Sphere3 sphere{Point3d{0, 0, 0}, Vec3d{0, 0, 1}, Vec3d{1, 0, 0}, 1.0};
  const Result<SsiPatch> patch = make_ssi_patch(SurfaceGeometry{sphere}, {0.0, 1.0}, {-0.5, 0.5});
  REQUIRE(patch.has_value());
  const Result<SsiPatchMesh> mesh = tessellate_ssi_patch(patch.value(), 1e-2, 12);
  REQUIRE(mesh.has_value());
  REQUIRE(mesh.value().normals.size() == mesh.value().points.size());
  REQUIRE(mesh.value().uv.size() == mesh.value().points.size());
  for (std::size_t i = 0; i < mesh.value().points.size(); ++i) {
    const Vec3d& normal = mesh.value().normals[i];
    CHECK(std::fabs(std::hypot(normal[0], normal[1], normal[2]) - 1.0) < 1e-12);
    for (int axis = 0; axis < 3; ++axis) {
      CHECK(normal[axis] == doctest::Approx(mesh.value().points[i][axis]).epsilon(1e-12));
    }
  }
}

TEST_CASE("ssi: a NURBS patch tessellates from its derivatives") {
  const NurbsSurface surface(
      1, 1, {0.0, 0.0, 1.0, 1.0}, {0.0, 0.0, 1.0, 1.0},
      {Point3d{0, 0, 0}, Point3d{0, 1, 0}, Point3d{1, 0, 0}, Point3d{1, 1, 0}},
      {1.0, 1.0, 1.0, 1.0});
  REQUIRE(surface.is_valid());
  const Result<SsiPatch> patch =
      make_ssi_patch(SurfaceGeometry{surface}, {0.0, 1.0}, {0.0, 1.0});
  REQUIRE(patch.has_value());
  const Result<SsiPatchMesh> mesh = tessellate_ssi_patch(patch.value(), 1e-9, 12);
  REQUIRE(mesh.has_value());
  CHECK(mesh.value().triangles.size() == 4);
  CHECK(mesh.value().points.size() == 5);
  for (const Point3d& p : mesh.value().points) CHECK(std::fabs(p[2]) < 1e-15);
  for (const Vec3d& normal : mesh.value().normals) {
    CHECK(std::fabs(std::fabs(normal[2]) - 1.0) < 1e-12);
  }
}

TEST_CASE("ssi: patch_bounds is empty for invalid inputs") {
  const Plane3 plane{};
  const Result<SsiPatch> patch = make_ssi_patch(SurfaceGeometry{plane}, {-1.0, 1.0}, {-1.0, 1.0});
  REQUIRE(patch.has_value());
  CHECK(patch_bounds(patch.value(), 0.0).is_empty());
  CHECK(patch_bounds(patch.value(), -1.0).is_empty());
  SsiPatch inverted = patch.value();
  inverted.u = {1.0, -1.0};
  CHECK(patch_bounds(inverted, 1e-3).is_empty());
}
