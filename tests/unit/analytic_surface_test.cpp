#include <doctest/doctest.h>

#include <cmath>
#include <limits>

#include <geom/surfaces/analytic.hpp>

TEST_CASE("cylinder evaluation, derivatives and outward normal") {
  const geom::Cylinder3 cyl{{0, 0, 0}, {0, 0, 1}, {1, 0, 0}, 2.0};
  CHECK(cyl.is_valid());
  const geom::Point3d p = cyl.eval(0.0, 3.0);
  CHECK(p[0] == doctest::Approx(2.0));
  CHECK(p[1] == doctest::Approx(0.0));
  CHECK(p[2] == doctest::Approx(3.0));
  CHECK(cyl.eval(3.14159265358979323846 / 2, 0.0)[1] == doctest::Approx(2.0).epsilon(1e-12));
  const geom::Vec3d d = cyl.du(0.0, 0.0);
  CHECK(d[0] == doctest::Approx(0.0));
  CHECK(d[1] == doctest::Approx(2.0));
  CHECK(cyl.dv(0.0, 0.0)[2] == doctest::Approx(1.0));
  const geom::Vec3d n0 = cyl.normal(0.0, 0.0);
  CHECK(n0[0] == doctest::Approx(1.0).epsilon(1e-12));
  CHECK(geom::dot(cyl.normal(1.2, 5.0), cyl.normal(1.2, 5.0)) == doctest::Approx(1.0).epsilon(1e-12));
  // Mid-longitude: the normal must rotate with u; a constant normal fails here.
  const geom::Vec3d n_mid = cyl.normal(3.14159265358979323846 / 2, 0.5);
  CHECK(n_mid[0] == doctest::Approx(0.0).epsilon(1e-12));
  CHECK(n_mid[1] == doctest::Approx(1.0).epsilon(1e-12));
  CHECK(n_mid[2] == doctest::Approx(0.0).epsilon(1e-12));
}

TEST_CASE("sphere normal is unit and points outward at poles") {
  const geom::Sphere3 s{{0, 0, 0}, {0, 0, 1}, {1, 0, 0}, 3.0};
  CHECK(s.is_valid());
  const geom::Point3d q = s.eval(0.0, 0.0);
  CHECK(q[0] == doctest::Approx(3.0));
  CHECK(q[1] == doctest::Approx(0.0));
  const geom::Vec3d np = s.normal(0.0, 3.14159265358979323846 / 2);
  CHECK(np[2] == doctest::Approx(1.0).epsilon(1e-12));
  CHECK(geom::dot(np, np) == doctest::Approx(1.0).epsilon(1e-12));
  // Pole normal is the directional limit: ±axis plus a cos(pi/2)*radial(u)
  // residual of about one ulp, compared at 1e-15 rather than with exact ==.
  CHECK(np[0] == doctest::Approx(0.0).epsilon(1e-15));
  CHECK(np[1] == doctest::Approx(0.0).epsilon(1e-15));
  // Mid-latitude, mid-longitude: the normal rotates with both u and v.
  constexpr double kPi = 3.14159265358979323846;
  const geom::Vec3d nm = s.normal(kPi / 2, kPi / 4);
  CHECK(nm[0] == doctest::Approx(0.0).epsilon(1e-12));
  CHECK(nm[1] == doctest::Approx(std::cos(kPi / 4)).epsilon(1e-12));
  CHECK(nm[2] == doctest::Approx(std::sin(kPi / 4)).epsilon(1e-12));
}

TEST_CASE("torus geometry is rejected when r >= R") {
  geom::Torus3 t{};
  t.major_radius = 1.0;
  t.minor_radius = 1.0;
  CHECK_FALSE(t.is_valid());
  t.minor_radius = 0.5;
  CHECK(t.is_valid());
}

TEST_CASE("cone derivative follows the slant") {
  const geom::Cone3 cone{{0, 0, 0}, {0, 0, 1}, {1, 0, 0}, 1.0, 0.25};
  CHECK(cone.is_valid());
  CHECK(cone.eval(0.0, 2.0)[2] == doctest::Approx(2.0).epsilon(1e-12));
  CHECK(cone.eval(0.0, 2.0)[0] ==
        doctest::Approx(1.0 - 2.0 * std::tan(0.25)).epsilon(1e-12));
  CHECK(cone.dv(0.0, 2.0)[0] == doctest::Approx(-std::tan(0.25)).epsilon(1e-12));
  CHECK(cone.du(0.0, 2.0)[1] == doctest::Approx(1.0 - 2.0 * std::tan(0.25)).epsilon(1e-12));
}

TEST_CASE("analytic surfaces reject non-finite fields and out-of-range parameters") {
  using Limit = std::numeric_limits<double>;

  const geom::Cylinder3 cylinder0{};
  REQUIRE(cylinder0.is_valid());
  geom::Cylinder3 cylinder = cylinder0;
  cylinder.origin = geom::Point3d{Limit::quiet_NaN(), 0, 0};
  CHECK_FALSE(cylinder.is_valid());
  cylinder = cylinder0;
  cylinder.axis = {0, Limit::infinity(), 0};
  CHECK_FALSE(cylinder.is_valid());
  cylinder = cylinder0;
  cylinder.x_axis = {Limit::quiet_NaN(), 0, 0};
  CHECK_FALSE(cylinder.is_valid());
  cylinder = cylinder0;
  cylinder.radius = 0.0;
  CHECK_FALSE(cylinder.is_valid());
  cylinder = cylinder0;
  cylinder.radius = -1.0;
  CHECK_FALSE(cylinder.is_valid());
  cylinder = cylinder0;
  cylinder.radius = Limit::infinity();
  CHECK_FALSE(cylinder.is_valid());

  const geom::Cone3 cone0{};
  REQUIRE(cone0.is_valid());
  geom::Cone3 cone = cone0;
  cone.origin = {0, 0, Limit::infinity()};
  CHECK_FALSE(cone.is_valid());
  cone = cone0;
  cone.axis = {0, 0, 0};
  CHECK_FALSE(cone.is_valid());
  cone = cone0;
  cone.x_axis = {0, 2, 0};
  CHECK_FALSE(cone.is_valid());
  cone = cone0;
  cone.radius = 0.0;
  CHECK_FALSE(cone.is_valid());
  cone = cone0;
  cone.half_angle = 0.0;
  CHECK_FALSE(cone.is_valid());
  cone = cone0;
  cone.half_angle = 3.14159265358979323846 / 2;
  CHECK_FALSE(cone.is_valid());
  cone = cone0;
  cone.half_angle = -0.25;
  CHECK_FALSE(cone.is_valid());
  cone = cone0;
  cone.half_angle = Limit::quiet_NaN();
  CHECK_FALSE(cone.is_valid());

  const geom::Sphere3 sphere0{};
  REQUIRE(sphere0.is_valid());
  geom::Sphere3 sphere = sphere0;
  sphere.center = {0, Limit::quiet_NaN(), 0};
  CHECK_FALSE(sphere.is_valid());
  sphere = sphere0;
  sphere.axis = {1, 1, 0};  // finite but not unit length
  CHECK_FALSE(sphere.is_valid());
  sphere = sphere0;
  sphere.x_axis = {0, 0, 1};  // parallel to axis
  CHECK_FALSE(sphere.is_valid());
  sphere = sphere0;
  sphere.radius = 0.0;
  CHECK_FALSE(sphere.is_valid());

  const geom::Torus3 torus0{};
  REQUIRE(torus0.is_valid());
  geom::Torus3 torus = torus0;
  torus.center = {0, 0, Limit::infinity()};
  CHECK_FALSE(torus.is_valid());
  torus = torus0;
  torus.minor_radius = 0.0;
  CHECK_FALSE(torus.is_valid());
  torus = torus0;
  torus.minor_radius = -0.5;
  CHECK_FALSE(torus.is_valid());
  torus = torus0;
  torus.major_radius = Limit::quiet_NaN();
  CHECK_FALSE(torus.is_valid());
}

TEST_CASE("analytic surface frames use the unit-norm tolerance") {
  const geom::Cylinder3 unit{};
  REQUIRE(unit.is_valid());
  // 5e-13 from unit length: valid under | |x| - 1 | <= 1e-12, but rejected by a
  // squared-norm check, because then |dot(x, x) - 1| is 1.0001e-12.
  geom::Cylinder3 c = unit;
  c.x_axis = {1.0 + 5e-13, 0, 0};
  CHECK(c.is_valid());
  c = unit;
  c.x_axis = {1.0 + 3e-12, 0, 0};
  CHECK_FALSE(c.is_valid());
  c = unit;
  c.axis = {0, 0, 1.0 - 3e-12};
  CHECK_FALSE(c.is_valid());
  // Unit x_axis tilted off the normal plane by 5e-13 is orthogonal enough;
  // 3e-12 is not.
  c = unit;
  c.x_axis = {1, 0, 5e-13};
  CHECK(c.is_valid());
  c = unit;
  c.x_axis = {1, 0, 3e-12};
  CHECK_FALSE(c.is_valid());
}

TEST_CASE("sphere derivatives and south pole normal") {
  constexpr double kPi = 3.14159265358979323846;
  const geom::Sphere3 s{{0, 0, 0}, {0, 0, 1}, {1, 0, 0}, 2.0};
  CHECK(s.is_valid());
  const geom::Vec3d d = s.du(0.0, 0.0);
  CHECK(d[0] == doctest::Approx(0.0));
  CHECK(d[1] == doctest::Approx(2.0).epsilon(1e-12));
  CHECK(d[2] == doctest::Approx(0.0));
  const geom::Vec3d e = s.dv(0.0, 0.0);
  CHECK(e[0] == doctest::Approx(0.0));
  CHECK(e[2] == doctest::Approx(2.0).epsilon(1e-12));
  const geom::Vec3d ns = s.normal(0.7, -kPi / 2);
  CHECK(ns[2] == doctest::Approx(-1.0).epsilon(1e-12));
  CHECK(geom::dot(ns, ns) == doctest::Approx(1.0).epsilon(1e-12));
  CHECK(ns[0] == doctest::Approx(0.0).epsilon(1e-15));
  CHECK(ns[1] == doctest::Approx(0.0).epsilon(1e-15));
  CHECK(s.eval(0.0, -kPi / 2)[2] == doctest::Approx(-2.0).epsilon(1e-12));
}

TEST_CASE("torus derivatives, evaluation and normals") {
  constexpr double kPi = 3.14159265358979323846;
  const geom::Torus3 t{{0, 0, 0}, {0, 0, 1}, {1, 0, 0}, 2.0, 0.5};
  CHECK(t.is_valid());
  const geom::Point3d outer = t.eval(0.0, 0.0);
  CHECK(outer[0] == doctest::Approx(2.5).epsilon(1e-12));
  CHECK(outer[1] == doctest::Approx(0.0));
  CHECK(outer[2] == doctest::Approx(0.0));
  const geom::Point3d top = t.eval(0.0, kPi / 2);
  CHECK(top[0] == doctest::Approx(2.0).epsilon(1e-12));
  CHECK(top[2] == doctest::Approx(0.5).epsilon(1e-12));
  const geom::Vec3d du = t.du(0.0, 0.0);
  CHECK(du[1] == doctest::Approx(2.5).epsilon(1e-12));
  CHECK(du[2] == doctest::Approx(0.0));
  const geom::Vec3d dv = t.dv(0.0, 0.0);
  CHECK(dv[0] == doctest::Approx(0.0));
  CHECK(dv[2] == doctest::Approx(0.5).epsilon(1e-12));
  const geom::Vec3d n0 = t.normal(0.0, 0.0);
  CHECK(n0[0] == doctest::Approx(1.0).epsilon(1e-12));
  const geom::Vec3d nt = t.normal(0.0, kPi / 2);
  CHECK(nt[2] == doctest::Approx(1.0).epsilon(1e-12));
  CHECK(geom::dot(t.normal(1.1, 0.3), t.normal(1.1, 0.3)) == doctest::Approx(1.0).epsilon(1e-12));
  // Mid-longitude: the normal must rotate with u; a constant normal fails here.
  const geom::Vec3d n_mid = t.normal(kPi / 2, 0.0);
  CHECK(n_mid[0] == doctest::Approx(0.0).epsilon(1e-12));
  CHECK(n_mid[1] == doctest::Approx(1.0).epsilon(1e-12));
  CHECK(n_mid[2] == doctest::Approx(0.0).epsilon(1e-12));
}

// A zero v-origin must be indistinguishable from the surface without the field:
// the closed forms are written out here exactly as the old implementation
// evaluated them, so an accidental extra rounding step in eval/du/dv/normal
// fails a coordinate comparison.
TEST_CASE("torus with a zero v-origin is the written closed form, coordinate for coordinate") {
  constexpr double kPi = 3.14159265358979323846;
  const geom::Point3d center{0.5, -1.5, 2.0};
  const geom::Torus3 t{center, {0.0, 1.0, 0.0}, {1.0, 0.0, 0.0}, 2.0, 0.5};
  CHECK(t.v_origin == 0.0);
  CHECK(geom::Torus3{}.v_origin == 0.0);
  CHECK(geom::Torus3{}.is_valid());

  for (const double u : {0.0, 0.3, kPi / 2.0, 2.1, -1.7}) {
    for (const double v : {0.0, 0.4, kPi / 2.0, -0.9, 2.6}) {
      const double c = std::cos(v);
      const double s = std::sin(v);
      const double w = 2.0 + 0.5 * c;
      // axis = +y, x_axis = +x, so y_axis = cross(axis, x_axis) = -z and
      // radial(u) = cos(u) * x_axis - sin(u) * z.
      const geom::Point3d expected{center[0] + std::cos(u) * w, center[1] + 0.5 * s,
                                   center[2] - std::sin(u) * w};
      const geom::Point3d actual = t.eval(u, v);
      CHECK(actual[0] == expected[0]);
      CHECK(actual[1] == expected[1]);
      CHECK(actual[2] == expected[2]);

      const geom::Vec3d expected_du{-std::sin(u) * w, 0.0, -std::cos(u) * w};
      const geom::Vec3d actual_du = t.du(u, v);
      CHECK(actual_du[0] == expected_du[0]);
      CHECK(actual_du[1] == expected_du[1]);
      CHECK(actual_du[2] == expected_du[2]);

      const geom::Vec3d expected_dv{-std::cos(u) * (0.5 * s), 0.5 * c, std::sin(u) * (0.5 * s)};
      const geom::Vec3d actual_dv = t.dv(u, v);
      CHECK(actual_dv[0] == expected_dv[0]);
      CHECK(actual_dv[1] == expected_dv[1]);
      CHECK(actual_dv[2] == expected_dv[2]);

      const geom::Vec3d expected_normal{std::cos(u) * c, s, -std::sin(u) * c};
      const geom::Vec3d actual_normal = t.normal(u, v);
      CHECK(actual_normal[0] == expected_normal[0]);
      CHECK(actual_normal[1] == expected_normal[1]);
      CHECK(actual_normal[2] == expected_normal[2]);
    }
  }
}

// eval(u, v) on the shifted chart is the base surface at v + v_origin, with the
// same identity for the partials and the normal. The comparisons are exact:
// shifting the parameter by an exact addition must not add a rounding step.
TEST_CASE("torus v-origin shifts the parameter chart") {
  constexpr double kPi = 3.14159265358979323846;
  constexpr double kOrigin = 0.7;
  const geom::Torus3 base{{0.5, -1.5, 2.0}, {0.0, 1.0, 0.0}, {1.0, 0.0, 0.0}, 2.0, 0.5};
  geom::Torus3 shifted = base;
  shifted.v_origin = kOrigin;
  CHECK(shifted.is_valid());

  const double inf = std::numeric_limits<double>::infinity();
  const double nan = std::numeric_limits<double>::quiet_NaN();
  CHECK_FALSE((geom::Torus3{base.center, base.axis, base.x_axis, 2.0, 0.5, inf}.is_valid()));
  CHECK_FALSE((geom::Torus3{base.center, base.axis, base.x_axis, 2.0, 0.5, nan}.is_valid()));

  for (const double u : {0.0, 0.35, 1.2, kPi, 4.9, -2.2}) {
    for (const double v : {0.0, 0.1, kPi / 2.0, 3.0, 5.5, -1.4}) {
      const geom::Point3d moved = shifted.eval(u, v);
      const geom::Point3d reference = base.eval(u, v + kOrigin);
      CHECK(moved[0] == reference[0]);
      CHECK(moved[1] == reference[1]);
      CHECK(moved[2] == reference[2]);

      const geom::Vec3d moved_du = shifted.du(u, v);
      const geom::Vec3d reference_du = base.du(u, v + kOrigin);
      CHECK(moved_du[0] == reference_du[0]);
      CHECK(moved_du[1] == reference_du[1]);
      CHECK(moved_du[2] == reference_du[2]);

      const geom::Vec3d moved_dv = shifted.dv(u, v);
      const geom::Vec3d reference_dv = base.dv(u, v + kOrigin);
      CHECK(moved_dv[0] == reference_dv[0]);
      CHECK(moved_dv[1] == reference_dv[1]);
      CHECK(moved_dv[2] == reference_dv[2]);

      const geom::Vec3d moved_normal = shifted.normal(u, v);
      const geom::Vec3d reference_normal = base.normal(u, v + kOrigin);
      CHECK(moved_normal[0] == reference_normal[0]);
      CHECK(moved_normal[1] == reference_normal[1]);
      CHECK(moved_normal[2] == reference_normal[2]);
    }
  }
}

TEST_CASE("cone normal is unit and perpendicular to the slant") {
  const geom::Cone3 cone{{0, 0, 0}, {0, 0, 1}, {1, 0, 0}, 1.0, 0.25};
  CHECK(cone.is_valid());
  const geom::Vec3d n = cone.normal(0.0, 0.0);
  CHECK(n[0] == doctest::Approx(std::cos(0.25)).epsilon(1e-12));
  CHECK(n[2] == doctest::Approx(std::sin(0.25)).epsilon(1e-12));
  CHECK(geom::dot(n, cone.dv(0.0, 2.0)) == doctest::Approx(0.0).epsilon(1e-12));
  CHECK(geom::dot(n, cone.du(0.0, 2.0)) == doctest::Approx(0.0).epsilon(1e-12));
  CHECK(geom::dot(n, n) == doctest::Approx(1.0).epsilon(1e-12));
  // Mid-longitude: the normal must rotate with u; a constant normal fails here.
  const geom::Vec3d n_mid = cone.normal(3.14159265358979323846 / 2, 0.0);
  CHECK(n_mid[0] == doctest::Approx(0.0).epsilon(1e-12));
  CHECK(n_mid[1] == doctest::Approx(std::cos(0.25)).epsilon(1e-12));
  CHECK(n_mid[2] == doctest::Approx(std::sin(0.25)).epsilon(1e-12));
}
