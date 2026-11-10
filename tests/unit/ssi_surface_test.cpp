#include <cmath>
#include <limits>
#include <numbers>

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

TEST_CASE("ssi: plane evaluation propagates plane frame errors") {
  const SurfaceGeometry nan_plane = Plane3{Vec3d{std::nan(""), 0.0, 1.0}, 0.0};
  const Result<SurfacePoint> nan_eval = surface_eval(nan_plane, 0.0, 0.0);
  REQUIRE_FALSE(nan_eval.has_value());
  CHECK(nan_eval.error().code == Errc::non_finite);
  const Result<Vec3d> nan_normal = surface_normal(nan_plane, 0.0, 0.0);
  REQUIRE_FALSE(nan_normal.has_value());
  CHECK(nan_normal.error().code == Errc::non_finite);

  const SurfaceGeometry zero_plane = Plane3{Vec3d{0.0, 0.0, 0.0}, 1.0};
  const Result<SurfacePoint> zero_eval = surface_eval(zero_plane, 0.0, 0.0);
  REQUIRE_FALSE(zero_eval.has_value());
  CHECK(zero_eval.error().code == Errc::degenerate);
  const Result<Vec3d> zero_normal = surface_normal(zero_plane, 0.0, 0.0);
  REQUIRE_FALSE(zero_normal.has_value());
  CHECK(zero_normal.error().code == Errc::degenerate);
}

TEST_CASE("ssi: plane normal matches its frame") {
  const Plane3 plane{Vec3d{1.0, 2.0, 3.0}, 4.0};
  const Result<PlaneFrame> frame = make_plane_frame(plane);
  REQUIRE(frame.has_value());
  const Result<Vec3d> normal = surface_normal(SurfaceGeometry{plane}, 1.0, -2.0);
  REQUIRE(normal.has_value());
  const Vec3d expected = cross(frame.value().u, frame.value().v);
  for (int i = 0; i < 3; ++i) {
    CHECK(normal.value()[i] == doctest::Approx(expected[i]).epsilon(1e-15));
  }
  CHECK(std::fabs(std::hypot(normal.value()[0], normal.value()[1], normal.value()[2]) - 1.0) <
        1e-15);
}

TEST_CASE("ssi: NURBS normal is unit and rejects a degenerate net") {
  const NurbsSurface plane_net(
      1, 1, {0.0, 0.0, 1.0, 1.0}, {0.0, 0.0, 1.0, 1.0},
      {Point3d{0, 0, 0}, Point3d{0, 1, 0}, Point3d{1, 0, 0}, Point3d{1, 1, 0}},
      {1.0, 1.0, 1.0, 1.0});
  REQUIRE(plane_net.is_valid());
  const Result<Vec3d> normal = surface_normal(SurfaceGeometry{plane_net}, 0.25, 0.75);
  REQUIRE(normal.has_value());
  CHECK(std::fabs(normal.value()[2]) == doctest::Approx(1.0).epsilon(1e-12));
  CHECK(std::fabs(normal.value()[0]) < 1e-12);
  CHECK(std::fabs(normal.value()[1]) < 1e-12);
  CHECK(std::fabs(std::hypot(normal.value()[0], normal.value()[1], normal.value()[2]) - 1.0) <
        1e-12);

  const NurbsSurface line_net(
      1, 1, {0.0, 0.0, 1.0, 1.0}, {0.0, 0.0, 1.0, 1.0},
      {Point3d{0, 0, 0}, Point3d{1, 0, 0}, Point3d{2, 0, 0}, Point3d{3, 0, 0}},
      {1.0, 1.0, 1.0, 1.0});
  REQUIRE(line_net.is_valid());
  const Result<Vec3d> degenerate = surface_normal(SurfaceGeometry{line_net}, 0.25, 0.75);
  REQUIRE_FALSE(degenerate.has_value());
  CHECK(degenerate.error().code == Errc::numeric_failure);
}

TEST_CASE("ssi: cone and invalid surfaces are handled") {
  const Cone3 cone{Point3d{0, 0, 0}, Vec3d{0, 0, 1}, Vec3d{1, 0, 0}, 2.0,
                   std::numbers::pi_v<double> / 6.0};
  REQUIRE(cone.is_valid());
  const SurfaceGeometry cone_surface = cone;
  CHECK(surface_periodic(cone_surface, 0));
  CHECK_FALSE(surface_periodic(cone_surface, 1));
  const Result<SurfacePoint> cone_eval = surface_eval(cone_surface, 0.3, 0.5);
  REQUIRE(cone_eval.has_value());
  CHECK(surface_point_finite(cone_eval.value()));
  const Result<Vec3d> cone_normal = surface_normal(cone_surface, 0.3, 0.5);
  REQUIRE(cone_normal.has_value());
  CHECK(std::fabs(std::hypot(cone_normal.value()[0], cone_normal.value()[1],
                            cone_normal.value()[2]) -
                  1.0) < 1e-15);

  const Cylinder3 zero_radius{Point3d{0, 0, 0}, Vec3d{0, 0, 1}, Vec3d{1, 0, 0}, 0.0};
  REQUIRE_FALSE(zero_radius.is_valid());
  const SurfaceGeometry invalid = zero_radius;
  const Result<SurfacePoint> bad_eval = surface_eval(invalid, 0.0, 0.0);
  REQUIRE_FALSE(bad_eval.has_value());
  CHECK(bad_eval.error().code == Errc::degenerate);
  const Result<Vec3d> bad_normal = surface_normal(invalid, 0.0, 0.0);
  REQUIRE_FALSE(bad_normal.has_value());
  CHECK(bad_normal.error().code == Errc::degenerate);
  const Result<SsiPatch> bad_patch = make_ssi_patch(invalid, {0.0, 1.0}, {0.0, 1.0});
  REQUIRE_FALSE(bad_patch.has_value());
  CHECK(bad_patch.error().code == Errc::degenerate);
}

TEST_CASE("ssi: surface queries reject non-finite parameters") {
  const SurfaceGeometry cylinder =
      Cylinder3{Point3d{0, 0, 0}, Vec3d{0, 0, 1}, Vec3d{1, 0, 0}, 1.0};
  const Result<SurfacePoint> nan_u = surface_eval(cylinder, std::nan(""), 0.0);
  REQUIRE_FALSE(nan_u.has_value());
  CHECK(nan_u.error().code == Errc::non_finite);
  const Result<SurfacePoint> inf_v =
      surface_eval(cylinder, 0.0, std::numeric_limits<double>::infinity());
  REQUIRE_FALSE(inf_v.has_value());
  CHECK(inf_v.error().code == Errc::non_finite);
  const Result<Vec3d> nan_normal = surface_normal(cylinder, 0.0, std::nan(""));
  REQUIRE_FALSE(nan_normal.has_value());
  CHECK(nan_normal.error().code == Errc::non_finite);
  const Result<Vec3d> inf_normal =
      surface_normal(cylinder, std::numeric_limits<double>::infinity(), 0.0);
  REQUIRE_FALSE(inf_normal.has_value());
  CHECK(inf_normal.error().code == Errc::non_finite);
}

TEST_CASE("ssi: wrap maps pi exactly and rejects bad axes") {
  const SurfaceGeometry cylinder =
      Cylinder3{Point3d{0, 0, 0}, Vec3d{0, 0, 1}, Vec3d{1, 0, 0}, 1.0};
  const Result<double> seam = wrap_surface_param(cylinder, 0, std::numbers::pi_v<double>);
  REQUIRE(seam.has_value());
  CHECK(seam.value() == -std::numbers::pi_v<double>);
  const Result<double> nan = wrap_surface_param(cylinder, 0, std::nan(""));
  REQUIRE_FALSE(nan.has_value());
  CHECK(nan.error().code == Errc::non_finite);
  const Result<double> bad_axis = wrap_surface_param(cylinder, 2, 0.5);
  REQUIRE_FALSE(bad_axis.has_value());
  CHECK(bad_axis.error().code == Errc::index_out_of_range);
  CHECK_FALSE(surface_periodic(cylinder, 2));
  CHECK_FALSE(surface_periodic(cylinder, -1));
}

TEST_CASE("ssi: patch spans and closed boundaries") {
  const SurfaceGeometry cylinder =
      Cylinder3{Point3d{0, 0, 0}, Vec3d{0, 0, 1}, Vec3d{1, 0, 0}, 1.0};
  const double two_pi = 2.0 * std::numbers::pi_v<double>;
  const Result<SsiPatch> full = make_ssi_patch(cylinder, {0.0, two_pi}, {0.0, 1.0});
  REQUIRE(full.has_value());
  CHECK(patch_contains(full.value(), two_pi, 0.5));
  CHECK(patch_contains(full.value(), -0.5 * two_pi, 0.5));
  CHECK_FALSE(patch_contains(full.value(), std::nan(""), 0.5));
  CHECK_FALSE(patch_contains(full.value(), 0.5, std::nan("")));

  const Result<SsiPatch> closed = make_ssi_patch(cylinder, {-0.2, 0.2}, {0.0, 1.0});
  REQUIRE(closed.has_value());
  CHECK(patch_contains(closed.value(), -0.2, 0.0));
  CHECK(patch_contains(closed.value(), 0.2, 1.0));
  CHECK(patch_contains(closed.value(), 0.0, 0.5));
  CHECK_FALSE(patch_contains(closed.value(), 0.2000001, 0.5));
  CHECK_FALSE(patch_contains(closed.value(), std::nan(""), 0.5));
  CHECK_FALSE(patch_contains(closed.value(), 0.0, std::nan("")));

  const SurfaceGeometry torus =
      Torus3{Point3d{0, 0, 0}, Vec3d{0, 0, 1}, Vec3d{1, 0, 0}, 3.0, 1.0};
  const Result<SsiPatch> wide_v = make_ssi_patch(torus, {0.0, 1.0}, {0.0, two_pi + 1e-12});
  REQUIRE_FALSE(wide_v.has_value());
  CHECK(wide_v.error().code == Errc::unsupported);
}

TEST_CASE("ssi: seam-crossing patch contains both sides") {
  const SurfaceGeometry cylinder =
      Cylinder3{Point3d{0, 0, 0}, Vec3d{0, 0, 1}, Vec3d{1, 0, 0}, 1.0};
  const Result<SsiPatch> patch = make_ssi_patch(
      cylinder, {3.0, 2.0 * std::numbers::pi_v<double> + 0.5}, {0.0, 1.0});
  REQUIRE(patch.has_value());
  CHECK(patch_contains(patch.value(), 3.5, 0.5));
  CHECK(patch_contains(patch.value(), 0.2, 0.5));
  // -1 rad is the same angle as 5.283... rad, which lies in [3, 2*pi) of the
  // seam-crossing interval; the gap is (0.5, 3.0).
  CHECK(patch_contains(patch.value(), -1.0, 0.5));
  CHECK_FALSE(patch_contains(patch.value(), 2.0, 0.5));
}

TEST_CASE("ssi: NURBS patch outside the surface domain is rejected") {
  const NurbsSurface surface(
      1, 1, {0.0, 0.0, 1.0, 1.0}, {0.0, 0.0, 1.0, 1.0},
      {Point3d{0, 0, 0}, Point3d{0, 1, 0}, Point3d{1, 0, 0}, Point3d{1, 1, 0}},
      {1.0, 1.0, 1.0, 1.0});
  REQUIRE(surface.is_valid());
  const SurfaceGeometry geometry = surface;
  CHECK(make_ssi_patch(geometry, {0.0, 1.0}, {0.0, 1.0}).has_value());
  const Result<SsiPatch> outside_u = make_ssi_patch(geometry, {-0.5, 0.5}, {0.0, 1.0});
  REQUIRE_FALSE(outside_u.has_value());
  CHECK(outside_u.error().code == Errc::unsupported);
  const Result<SsiPatch> outside_v = make_ssi_patch(geometry, {0.0, 1.0}, {0.5, 1.5});
  REQUIRE_FALSE(outside_v.has_value());
  CHECK(outside_v.error().code == Errc::unsupported);
}

TEST_CASE("ssi: surface_point_finite rejects non-finite components") {
  SurfacePoint point;
  point.p = Point3d{0.0, std::nan(""), 0.0};
  CHECK_FALSE(surface_point_finite(point));
  point.p = Point3d{0.0, 0.0, 0.0};
  point.du = Vec3d{std::numeric_limits<double>::infinity(), 0.0, 0.0};
  CHECK_FALSE(surface_point_finite(point));
  point.du = Vec3d{};
  point.dv = Vec3d{0.0, -std::numeric_limits<double>::infinity(), 0.0};
  CHECK_FALSE(surface_point_finite(point));
}

TEST_CASE("ssi: NURBS normal rejects an overflowing cross product") {
  const NurbsSurface huge(
      1, 1, {0.0, 0.0, 1.0, 1.0}, {0.0, 0.0, 1.0, 1.0},
      {Point3d{-1e300, -1e300, 0.0}, Point3d{1e300, -1e300, 0.0},
       Point3d{-1e300, 1e300, 0.0}, Point3d{1e300, 1e300, 0.0}},
      {1.0, 1.0, 1.0, 1.0});
  REQUIRE(huge.is_valid());
  const Result<Vec3d> normal = surface_normal(SurfaceGeometry{huge}, 0.5, 0.5);
  REQUIRE_FALSE(normal.has_value());
  CHECK(normal.error().code == Errc::non_finite);
}

TEST_CASE("ssi: NURBS normal rejects an invalid surface") {
  const NurbsSurface broken(
      1, 1, {0.0, 0.0, 1.0, 1.0}, {0.0, 0.0, 1.0, 1.0},
      {Point3d{0, 0, 0}, Point3d{0, 1, 0}, Point3d{1, 0, 0}}, {1.0, 1.0, 1.0});
  REQUIRE_FALSE(broken.is_valid());
  const Result<Vec3d> normal = surface_normal(SurfaceGeometry{broken}, 0.25, 0.75);
  REQUIRE_FALSE(normal.has_value());
  CHECK(normal.error().code == Errc::degenerate);
}

TEST_CASE("ssi: patch validation propagates plane frame errors") {
  const SurfaceGeometry nan_plane = Plane3{Vec3d{std::nan(""), 0.0, 1.0}, 0.0};
  const Result<SsiPatch> nan_patch = make_ssi_patch(nan_plane, {0.0, 1.0}, {0.0, 1.0});
  REQUIRE_FALSE(nan_patch.has_value());
  CHECK(nan_patch.error().code == Errc::non_finite);

  const SurfaceGeometry zero_plane = Plane3{Vec3d{0.0, 0.0, 0.0}, 1.0};
  const Result<SsiPatch> zero_patch = make_ssi_patch(zero_plane, {0.0, 1.0}, {0.0, 1.0});
  REQUIRE_FALSE(zero_patch.has_value());
  CHECK(zero_patch.error().code == Errc::degenerate);
}
