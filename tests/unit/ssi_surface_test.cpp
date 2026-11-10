#include <cmath>

#include <doctest/doctest.h>
#include <geom/brep/detail/ssi_patch.hpp>
#include <geom/brep/detail/ssi_surface.hpp>
#include <geom/brep/plane_frame.hpp>

using namespace geom::v1;

TEST_CASE("ssi: plane patch evaluation uses the plane frame") {
  const Plane3 plane{{0.0, 0.0, 1.0}, 2.0};
  const Result<PlaneFrame> frame = make_plane_frame(plane);
  REQUIRE(frame.has_value());
  const SurfaceGeometry surface = plane;
  const Result<SurfacePoint> eval = surface_eval(surface, 2.0, -3.0);
  REQUIRE(eval.has_value());
  const Point3d expected = frame.value().origin + frame.value().u * 2.0 + frame.value().v * (-3.0);
  for (int i = 0; i < 3; ++i) CHECK(eval.value().p[i] == doctest::Approx(expected[i]).epsilon(1e-15));
  CHECK(dot(eval.value().du, frame.value().u) == doctest::Approx(1.0));
  CHECK(dot(eval.value().dv, frame.value().v) == doctest::Approx(1.0));
  CHECK(std::fabs(eval.value().du[2]) < 1e-15);
}

TEST_CASE("ssi: analytic surfaces expose derivatives and periodicity") {
  const Cylinder3 cylinder{Point3d{0, 0, 0}, Vec3d{0, 0, 1}, Vec3d{1, 0, 0}, 2.0};
  REQUIRE(cylinder.is_valid());
  const SurfaceGeometry surface = cylinder;
  CHECK(surface_periodic(surface, 0));
  CHECK_FALSE(surface_periodic(surface, 1));
  const Result<SurfacePoint> eval = surface_eval(surface, 0.3, 1.7);
  REQUIRE(eval.has_value());
  CHECK(std::fabs(std::hypot(eval.value().p[0], eval.value().p[1]) - 2.0) < 1e-15);
  CHECK(eval.value().p[2] == doctest::Approx(1.7));
  CHECK(std::fabs(dot(eval.value().du, eval.value().dv)) < 1e-15);

  const Sphere3 sphere{Point3d{1, 2, 3}, Vec3d{0, 0, 1}, Vec3d{1, 0, 0}, 2.0};
  REQUIRE(sphere.is_valid());
  CHECK(surface_periodic(SurfaceGeometry{sphere}, 0));
  CHECK_FALSE(surface_periodic(SurfaceGeometry{sphere}, 1));

  const Torus3 torus{Point3d{0, 0, 0}, Vec3d{0, 0, 1}, Vec3d{1, 0, 0}, 3.0, 1.0};
  REQUIRE(torus.is_valid());
  CHECK(surface_periodic(SurfaceGeometry{torus}, 0));
  CHECK(surface_periodic(SurfaceGeometry{torus}, 1));
}

TEST_CASE("ssi: wrap_surface_param folds periodic parameters into [-pi, pi)") {
  const Cylinder3 cylinder{Point3d{0, 0, 0}, Vec3d{0, 0, 1}, Vec3d{1, 0, 0}, 1.0};
  const SurfaceGeometry surface = cylinder;
  const Result<double> wrapped = wrap_surface_param(surface, 0, 6.5);
  REQUIRE(wrapped.has_value());
  CHECK(wrapped.value() == doctest::Approx(6.5 - 2.0 * 3.14159265358979323846).epsilon(1e-15));
  const Result<double> negative = wrap_surface_param(surface, 0, -0.5);
  REQUIRE(negative.has_value());
  CHECK(negative.value() == doctest::Approx(-0.5));
  const Result<double> plane_wrap = wrap_surface_param(SurfaceGeometry{Plane3{}}, 0, 6.5);
  REQUIRE(plane_wrap.has_value());
  CHECK(plane_wrap.value() == doctest::Approx(6.5));
}

TEST_CASE("ssi: patch validation rejects inverted, non-finite and oversized domains") {
  const Sphere3 sphere{Point3d{0, 0, 0}, Vec3d{0, 0, 1}, Vec3d{1, 0, 0}, 1.0};
  const SurfaceGeometry surface = sphere;
  CHECK(make_ssi_patch(surface, {0.0, 2.0}, {-1.0, 1.0}).has_value());
  const Result<SsiPatch> inverted = make_ssi_patch(surface, {2.0, 0.0}, {-1.0, 1.0});
  REQUIRE_FALSE(inverted.has_value());
  CHECK(inverted.error().code == Errc::degenerate);
  const Result<SsiPatch> nan = make_ssi_patch(surface, {0.0, std::nan("")}, {-1.0, 1.0});
  REQUIRE_FALSE(nan.has_value());
  CHECK(nan.error().code == Errc::non_finite);
  const Result<SsiPatch> wide = make_ssi_patch(surface, {0.0, 7.0}, {-1.0, 1.0});
  REQUIRE_FALSE(wide.has_value());
  CHECK(wide.error().code == Errc::unsupported);
}

TEST_CASE("ssi: patch_contains wraps periodic coordinates") {
  const Cylinder3 cylinder{Point3d{0, 0, 0}, Vec3d{0, 0, 1}, Vec3d{1, 0, 0}, 1.0};
  const Result<SsiPatch> patch = make_ssi_patch(SurfaceGeometry{cylinder}, {-0.2, 0.2}, {0.0, 1.0});
  REQUIRE(patch.has_value());
  CHECK(patch_contains(patch.value(), 0.1, 0.5));
  CHECK(patch_contains(patch.value(), 6.283185307179586 - 0.1, 0.5));
  CHECK_FALSE(patch_contains(patch.value(), 0.3, 0.5));
  CHECK_FALSE(patch_contains(patch.value(), 0.1, 1.2));
}

TEST_CASE("ssi: surface_normal stays unit at analytic poles") {
  const Sphere3 sphere{Point3d{0, 0, 0}, Vec3d{0, 0, 1}, Vec3d{1, 0, 0}, 2.0};
  const Result<Vec3d> south = surface_normal(
      SurfaceGeometry{sphere}, 0.0, -0.5 * 3.14159265358979323846);
  REQUIRE(south.has_value());
  CHECK(south.value()[2] == doctest::Approx(-1.0).epsilon(1e-12));
  CHECK(std::fabs(std::hypot(south.value()[0], south.value()[1], south.value()[2]) - 1.0) < 1e-12);
}

TEST_CASE("ssi: NURBS surface evaluation matches the surface directly") {
  const NurbsSurface surface(1, 1, {0.0, 0.0, 1.0, 1.0}, {0.0, 0.0, 1.0, 1.0},
                             {Point3d{0, 0, 0}, Point3d{0, 1, 0}, Point3d{1, 0, 0}, Point3d{1, 1, 0}},
                             {1.0, 1.0, 1.0, 1.0});
  REQUIRE(surface.is_valid());
  const Result<SurfacePoint> eval = surface_eval(SurfaceGeometry{surface}, 0.25, 0.75);
  REQUIRE(eval.has_value());
  const Point3d direct = surface.eval(0.25, 0.75);
  for (int i = 0; i < 3; ++i) CHECK(eval.value().p[i] == doctest::Approx(direct[i]).epsilon(1e-15));
}
