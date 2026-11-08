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
