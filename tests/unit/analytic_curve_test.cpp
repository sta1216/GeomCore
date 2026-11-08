#include <doctest/doctest.h>

#include <limits>

#include <geom/curves/analytic_curve.hpp>

TEST_CASE("line curve evaluates and bounds") {
  const geom::LineCurve3 line{{1, 2, 3}, {0, 1, 0}};
  CHECK(line.is_valid());
  const geom::Point3d p = line.eval(4.0);
  CHECK(p[0] == doctest::Approx(1.0));
  CHECK(p[1] == doctest::Approx(6.0));
  CHECK(p[2] == doctest::Approx(3.0));
  CHECK(line.derivative(4.0)[1] == doctest::Approx(1.0));
  const geom::Box3d b = line.bounds(0.0, 4.0);
  CHECK(b.lo[0] == doctest::Approx(1.0));
  CHECK(b.lo[1] == doctest::Approx(2.0));
  CHECK(b.hi[0] == doctest::Approx(1.0));
  CHECK(b.hi[1] == doctest::Approx(6.0));
}

TEST_CASE("circle curve quarter arc bounds are exact") {
  const geom::CircleCurve3 c{{0, 0, 0}, {1, 0, 0}, {0, 1, 0}, 2.0};
  CHECK(c.is_valid());
  const geom::Point3d p0 = c.eval(0.0);
  CHECK(p0[0] == doctest::Approx(2.0));
  CHECK(p0[1] == doctest::Approx(0.0));
  const geom::Box3d b = c.bounds(0.0, 3.14159265358979323846 / 2);
  CHECK(b.lo[0] == doctest::Approx(0.0).epsilon(1e-12));
  CHECK(b.hi[0] == doctest::Approx(2.0).epsilon(1e-12));
  CHECK(b.lo[1] == doctest::Approx(0.0).epsilon(1e-12));
  CHECK(b.hi[1] == doctest::Approx(2.0).epsilon(1e-12));
}

TEST_CASE("circle curve full turn bounds are exact") {
  constexpr double kPi = 3.14159265358979323846;
  const geom::CircleCurve3 c{{0, 0, 0}, {1, 0, 0}, {0, 1, 0}, 2.0};
  CHECK(c.is_valid());
  // Every extremum of a full turn is interior to the arc: an implementation
  // that only boxes eval(a0) and eval(a1) returns lo[1] == hi[1] == 0 here.
  const geom::Box3d full = c.bounds(0.0, 2.0 * kPi);
  CHECK(full.lo[0] == doctest::Approx(-2.0).epsilon(1e-12));
  CHECK(full.hi[0] == doctest::Approx(2.0).epsilon(1e-12));
  CHECK(full.lo[1] == doctest::Approx(-2.0).epsilon(1e-12));
  CHECK(full.hi[1] == doctest::Approx(2.0).epsilon(1e-12));
  // Same arc shifted by pi, so both stationary angles sit at the same distance
  // from the endpoints.
  const geom::Box3d shifted = c.bounds(-kPi, kPi);
  CHECK(shifted.lo[0] == doctest::Approx(-2.0).epsilon(1e-12));
  CHECK(shifted.hi[0] == doctest::Approx(2.0).epsilon(1e-12));
  CHECK(shifted.lo[1] == doctest::Approx(-2.0).epsilon(1e-12));
  CHECK(shifted.hi[1] == doctest::Approx(2.0).epsilon(1e-12));
}

TEST_CASE("circle curve axis tolerance is measured on the axis norm") {
  const geom::CircleCurve3 unit{};
  REQUIRE(unit.is_valid());
  // 5e-13 from unit length: valid under | |x| - 1 | <= 1e-12, but rejected by a
  // squared-norm check, because then |dot(x, x) - 1| is 1.0001e-12.
  geom::CircleCurve3 c = unit;
  c.x_axis = {1.0 + 5e-13, 0, 0};
  CHECK(c.is_valid());
  c = unit;
  c.y_axis = {0, 1.0 + 5e-13, 0};
  CHECK(c.is_valid());
  c = unit;
  c.x_axis = {1.0 - 1e-12, 0, 0};
  CHECK(c.is_valid());
  c = unit;
  c.x_axis = {1.0 + 3e-12, 0, 0};
  CHECK_FALSE(c.is_valid());
  c = unit;
  c.y_axis = {0, 1.0 - 3e-12, 0};
  CHECK_FALSE(c.is_valid());
}

TEST_CASE("analytic curves reject non-finite and unnormalized data") {
  geom::LineCurve3 line{{0, 0, 0}, {0, 2, 0}};
  CHECK_FALSE(line.is_valid());
  line.direction = {0, 0, 0};
  CHECK_FALSE(line.is_valid());
  line.direction = {std::numeric_limits<double>::quiet_NaN(), 0, 0};
  CHECK_FALSE(line.is_valid());
  geom::CircleCurve3 c{};
  c.radius = std::numeric_limits<double>::infinity();
  CHECK_FALSE(c.is_valid());
  c = geom::CircleCurve3{};
  c.y_axis = {1, 0, 0};  // parallel to x_axis
  CHECK_FALSE(c.is_valid());
}

TEST_CASE("analytic curves reject non-finite fields, non-unit axes and radius <= 0") {
  using Limit = std::numeric_limits<double>;
  const geom::LineCurve3 line0{};
  REQUIRE(line0.is_valid());
  geom::LineCurve3 line = line0;
  line.origin = geom::Point3d{Limit::quiet_NaN(), 0, 0};
  CHECK_FALSE(line.is_valid());
  line = line0;
  line.origin = geom::Point3d{0, Limit::infinity(), 0};
  CHECK_FALSE(line.is_valid());
  line = line0;
  line.origin = geom::Point3d{0, 0, Limit::quiet_NaN()};
  CHECK_FALSE(line.is_valid());
  line = line0;
  line.direction = {0, Limit::infinity(), 0};
  CHECK_FALSE(line.is_valid());
  line = line0;
  line.direction = {0, 0, Limit::quiet_NaN()};
  CHECK_FALSE(line.is_valid());

  const geom::CircleCurve3 circle0{};
  REQUIRE(circle0.is_valid());
  geom::CircleCurve3 c = circle0;
  c.center = geom::Point3d{Limit::quiet_NaN(), 0, 0};
  CHECK_FALSE(c.is_valid());
  c = circle0;
  c.center = geom::Point3d{0, Limit::infinity(), 0};
  CHECK_FALSE(c.is_valid());
  c = circle0;
  c.center = geom::Point3d{0, 0, Limit::quiet_NaN()};
  CHECK_FALSE(c.is_valid());
  c = circle0;
  c.x_axis = {Limit::quiet_NaN(), 0, 0};
  CHECK_FALSE(c.is_valid());
  c = circle0;
  c.x_axis = {0, 2, 0};  // finite but not unit length
  CHECK_FALSE(c.is_valid());
  c = circle0;
  c.y_axis = {0, Limit::infinity(), 0};
  CHECK_FALSE(c.is_valid());
  c = circle0;
  c.y_axis = {0, 2, 0};  // finite but not unit length
  CHECK_FALSE(c.is_valid());
  c = circle0;
  c.radius = 0.0;
  CHECK_FALSE(c.is_valid());
  c = circle0;
  c.radius = -1.0;
  CHECK_FALSE(c.is_valid());
}
