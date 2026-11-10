#include <cmath>
#include <vector>

#include <doctest/doctest.h>
#include <geom/brep/ssi.hpp>

using namespace geom::v1;

namespace {
NurbsSurface planar_nurbs() {
  return NurbsSurface(1, 1, {0.0, 0.0, 1.0, 1.0}, {0.0, 0.0, 1.0, 1.0},
                      {Point3d{0, 0, 0}, Point3d{0, 1, 0}, Point3d{1, 0, 0}, Point3d{1, 1, 0}},
                      {1.0, 1.0, 1.0, 1.0});
}

NurbsSurface wavy_nurbs() {
  std::vector<Point3d> control;
  for (int i = 0; i < 4; ++i) {
    for (int j = 0; j < 4; ++j) {
      const double x = static_cast<double>(i) / 3.0;
      const double y = static_cast<double>(j) / 3.0;
      control.push_back(Point3d{x, y, std::sin(3.0 * x) * std::sin(3.0 * y)});
    }
  }
  return NurbsSurface(3, 3, {0, 0, 0, 0, 1, 1, 1, 1}, {0, 0, 0, 0, 1, 1, 1, 1}, control,
                      std::vector<double>(16, 1.0));
}
}  // namespace

TEST_CASE("ssi: tangent and coincident inputs return diagnostics instead of curves") {
  const Sphere3 sphere{Point3d{0, 0, 0}, Vec3d{0, 0, 1}, Vec3d{1, 0, 0}, 1.0};
  const Result<SsiPatch> sphere_patch = make_ssi_patch(
      SurfaceGeometry{sphere}, {0, 2 * 3.14159265358979323846},
      {-0.5 * 3.14159265358979323846, 0.5 * 3.14159265358979323846});
  const Result<SsiPatch> tangent_plane = make_ssi_patch(
      SurfaceGeometry{Plane3{{0, 0, 1}, -1.0}}, {-2, 2}, {-2, 2});
  REQUIRE(sphere_patch.has_value());
  REQUIRE(tangent_plane.has_value());
  const Result<SsiResult> tangent =
      ssi(sphere_patch.value(), tangent_plane.value(), SsiOptions{1e-9});
  REQUIRE(tangent.has_value());
  CHECK(tangent.value().curves.empty());
  REQUIRE(tangent.value().issues.size() == 1);
  CHECK(tangent.value().issues[0].code == SsiDiagnostic::tangent_contact);

  const Result<SsiPatch> plane_a = make_ssi_patch(
      SurfaceGeometry{Plane3{{0, 0, 1}, 2.0}}, {-2, 2}, {-2, 2});
  const Result<SsiPatch> plane_b = make_ssi_patch(
      SurfaceGeometry{Plane3{{0, 0, -1}, -2.0}}, {-2, 2}, {-2, 2});
  REQUIRE(plane_a.has_value());
  REQUIRE(plane_b.has_value());
  const Result<SsiResult> coincident =
      ssi(plane_a.value(), plane_b.value(), SsiOptions{1e-9});
  REQUIRE(coincident.has_value());
  CHECK(coincident.value().curves.empty());
  REQUIRE(coincident.value().issues.size() == 1);
  CHECK(coincident.value().issues[0].code == SsiDiagnostic::coincident);
}

TEST_CASE("ssi: a cylinder patch intersects a plane in two lines") {
  const Cylinder3 cylinder{Point3d{0, 0, 0}, Vec3d{0, 0, 1}, Vec3d{1, 0, 0}, 1.0};
  const Result<SsiPatch> a = make_ssi_patch(
      SurfaceGeometry{cylinder}, {0, 2 * 3.14159265358979323846}, {0, 1});
  const Result<SsiPatch> b = make_ssi_patch(
      SurfaceGeometry{Plane3{{1, 0, 0}, 0.5}}, {-2, 2}, {-2, 2});
  REQUIRE(a.has_value());
  REQUIRE(b.has_value());
  const Result<SsiResult> result = ssi(a.value(), b.value(), SsiOptions{1e-9});
  REQUIRE(result.has_value());
  CHECK(result.value().curves.size() == 2);
  for (const SsiCurve& curve : result.value().curves) {
    for (const SsiSegment& s : curve.segments) {
      CHECK(s.residual <= 1e-9);
      CHECK(std::fabs(s.a[0] - 0.5) < 1e-8);
      CHECK(std::fabs(s.b[0] - 0.5) < 1e-8);
    }
  }
}

TEST_CASE("ssi: a cone patch intersects a plane in a circle") {
  const Cone3 cone{Point3d{0, 0, 0}, Vec3d{0, 0, 1}, Vec3d{1, 0, 0}, 1.0, 0.5};
  const Result<SsiPatch> a = make_ssi_patch(
      SurfaceGeometry{cone}, {0, 2 * 3.14159265358979323846}, {0, 1});
  const Result<SsiPatch> b = make_ssi_patch(
      SurfaceGeometry{Plane3{{0, 0, 1}, 0.5}}, {-2, 2}, {-2, 2});
  REQUIRE(a.has_value());
  REQUIRE(b.has_value());
  const Result<SsiResult> result = ssi(a.value(), b.value(), SsiOptions{1e-9});
  REQUIRE(result.has_value());
  REQUIRE(result.value().curves.size() == 1);
  CHECK(result.value().curves[0].closed);
  const double expected_radius = 1.0 - 0.5 * std::tan(0.5);
  for (const SsiSegment& s : result.value().curves[0].segments) {
    CHECK(s.residual <= 1e-9);
    CHECK(std::fabs(s.a[2] - 0.5) < 1e-8);
    CHECK(std::fabs(std::hypot(s.a[0], s.a[1]) - expected_radius) < 1e-8);
  }
}

TEST_CASE("ssi: parallel offset planes return no curve and no issue") {
  const Result<SsiPatch> a = make_ssi_patch(
      SurfaceGeometry{Plane3{{0, 0, 1}, 0.0}}, {-1, 1}, {-1, 1});
  const Result<SsiPatch> b = make_ssi_patch(
      SurfaceGeometry{Plane3{{0, 0, 1}, 1.0}}, {-1, 1}, {-1, 1});
  REQUIRE(a.has_value());
  REQUIRE(b.has_value());
  const Result<SsiResult> result = ssi(a.value(), b.value(), SsiOptions{1e-9});
  REQUIRE(result.has_value());
  CHECK(result.value().curves.empty());
  CHECK(result.value().issues.empty());
}

TEST_CASE("ssi: a planar NURBS patch intersects a plane in one line") {
  const NurbsSurface patch = planar_nurbs();
  REQUIRE(patch.is_valid());
  const Result<SsiPatch> a = make_ssi_patch(SurfaceGeometry{patch}, {0, 1}, {0, 1});
  const Result<SsiPatch> b = make_ssi_patch(SurfaceGeometry{Plane3{{1, 0, 0}, 0.5}}, {-2, 2}, {-2, 2});
  REQUIRE(a.has_value());
  REQUIRE(b.has_value());
  const Result<SsiResult> result = ssi(a.value(), b.value(), SsiOptions{1e-9});
  REQUIRE(result.has_value());
  REQUIRE(result.value().curves.size() == 1);
  CHECK(result.value().curves[0].segments.size() >= 1);
  for (const SsiSegment& s : result.value().curves[0].segments) {
    CHECK(std::fabs(s.a[0] - 0.5) < 1e-8);
    CHECK(std::fabs(s.b[0] - 0.5) < 1e-8);
  }
}

TEST_CASE("ssi: a wavy NURBS patch intersects a plane within tolerance") {
  const NurbsSurface patch = wavy_nurbs();
  REQUIRE(patch.is_valid());
  const Result<SsiPatch> a = make_ssi_patch(SurfaceGeometry{patch}, {0, 1}, {0, 1});
  const Result<SsiPatch> b = make_ssi_patch(SurfaceGeometry{Plane3{{0, 0, 1}, 0.25}}, {-2, 2}, {-2, 2});
  REQUIRE(a.has_value());
  REQUIRE(b.has_value());
  const Result<SsiResult> result = ssi(a.value(), b.value(), SsiOptions{1e-8});
  REQUIRE(result.has_value());
  CHECK_FALSE(result.value().curves.empty());
  for (const SsiCurve& curve : result.value().curves) {
    for (const SsiSegment& s : curve.segments) {
      CHECK(s.residual <= 1e-8);
      CHECK(std::fabs(s.a[2] - 0.25) < 1e-8);
      CHECK(std::fabs(s.b[2] - 0.25) < 1e-8);
    }
  }
}
