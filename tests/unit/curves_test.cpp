#include <doctest/doctest.h>

#include <cmath>
#include <cstddef>
#include <vector>

#include <geom/curves/bezier.hpp>
#include <geom/curves/nurbs_curve.hpp>

using geom::v1::BezierCurve;
using geom::v1::NurbsCurve;
using geom::v1::Point2d;
using geom::v1::Point3d;

TEST_CASE("quadratic Bezier evaluates to the midpoint formula") {
  const BezierCurve<double, 2> c({Point2d{0, 0}, Point2d{1, 2}, Point2d{2, 0}});
  const Point2d p = c.eval(0.5);
  CHECK(p[0] == doctest::Approx(1.0));
  CHECK(p[1] == doctest::Approx(1.0));
  CHECK(c.eval(0.0)[0] == doctest::Approx(0.0));
  CHECK(c.eval(1.0)[0] == doctest::Approx(2.0));
}

TEST_CASE("Bezier derivative matches central differences") {
  const BezierCurve<double, 3> c({Point3d{0, 0, 0}, Point3d{1, 2, 0}, Point3d{2, 2, 1},
                                  Point3d{3, 0, 1}});
  const double h = 1e-5;
  const double t = 0.37;
  const auto fd = (c.eval(t + h) - c.eval(t - h)) * (1.0 / (2.0 * h));
  const auto exact = c.derivative(t);
  CHECK(exact[0] == doctest::Approx(fd[0]).epsilon(1e-6));
  CHECK(exact[1] == doctest::Approx(fd[1]).epsilon(1e-6));
  CHECK(exact[2] == doctest::Approx(fd[2]).epsilon(1e-6));
}

TEST_CASE("rational quadratic NURBS represents a quarter circle") {
  const double w = std::sqrt(0.5);
  const NurbsCurve c(2, {0.0, 0.0, 0.0, 1.0, 1.0, 1.0},
                     {Point3d{1, 0, 0}, Point3d{1, 1, 0}, Point3d{0, 1, 0}},
                     {1.0, w, 1.0});
  REQUIRE(c.is_valid());
  const auto [u0, u1] = c.domain();
  CHECK(u0 == doctest::Approx(0.0));
  CHECK(u1 == doctest::Approx(1.0));

  const Point3d p = c.eval(0.5);
  CHECK(p[0] == doctest::Approx(std::sqrt(0.5)).epsilon(1e-14));
  CHECK(p[1] == doctest::Approx(std::sqrt(0.5)).epsilon(1e-14));

  const double h = 1e-6;
  const auto fd = (c.eval(0.5 + h) - c.eval(0.5 - h)) * (1.0 / (2.0 * h));
  const auto exact = c.derivative(0.5);
  CHECK(exact[0] == doctest::Approx(fd[0]).epsilon(1e-5));
  CHECK(exact[1] == doctest::Approx(fd[1]).epsilon(1e-5));
}

TEST_CASE("NURBS with unit weights equals a B-spline") {
  const NurbsCurve c(2, {0.0, 0.0, 0.0, 1.0, 2.0, 2.0, 2.0},
                     {Point3d{0, 0, 0}, Point3d{1, 1, 0}, Point3d{2, 1, 0}, Point3d{3, 0, 0}},
                     {1.0, 1.0, 1.0, 1.0});
  REQUIRE(c.is_valid());
  CHECK(c.degree() == 2);
  CHECK(c.control_count() == 4);
  const Point3d p = c.eval(1.0);
  CHECK(p[0] == doctest::Approx(1.5).epsilon(1e-14));
}

TEST_CASE("invalid NURBS data is rejected") {
  const NurbsCurve bad(2, {0.0, 0.0, 1.0, 1.0}, {Point3d{0, 0, 0}, Point3d{1, 0, 0}},
                       {1.0, 1.0});
  CHECK_FALSE(bad.is_valid());
}

TEST_CASE("Bezier degenerate inputs and endpoint derivatives") {
  const BezierCurve<double, 2> empty({});
  CHECK(empty.size() == 0);
  CHECK(empty.eval(0.25).v[0] == 0.0);
  CHECK(empty.derivative(0.25).v[1] == 0.0);

  const BezierCurve<double, 2> one({Point2d{3.0, -4.0}});
  CHECK(one.size() == 1);
  CHECK(one.eval(0.75).v[0] == 3.0);
  CHECK(one.eval(0.75).v[1] == -4.0);
  CHECK(one.derivative(0.5).v[0] == 0.0);

  // For a quadratic, C'(0) = 2 (P1 - P0) and C'(1) = 2 (P2 - P1).
  const BezierCurve<double, 2> quad({Point2d{0.0, 0.0}, Point2d{1.0, 2.0}, Point2d{2.0, 0.0}});
  CHECK(quad.size() == 3);
  const auto d0 = quad.derivative(0.0);
  CHECK(d0.v[0] == doctest::Approx(2.0));
  CHECK(d0.v[1] == doctest::Approx(4.0));
  const auto d1 = quad.derivative(1.0);
  CHECK(d1.v[0] == doctest::Approx(2.0));
  CHECK(d1.v[1] == doctest::Approx(-4.0));
  const auto dm = quad.derivative(0.5);
  CHECK(dm.v[0] == doctest::Approx(2.0));
  CHECK(dm.v[1] == doctest::Approx(0.0));
}

TEST_CASE("NURBS derivative at the ends of a rational quarter circle") {
  const double w = std::sqrt(0.5);
  const NurbsCurve c(2, {0.0, 0.0, 0.0, 1.0, 1.0, 1.0},
                     {Point3d{1, 0, 0}, Point3d{1, 1, 0}, Point3d{0, 1, 0}}, {1.0, w, 1.0});
  REQUIRE(c.is_valid());
  // C'(0) = (0, 2w) and C'(1) = (-2w, 0) with w = sqrt(1/2).
  const auto d0 = c.derivative(0.0);
  CHECK(d0.v[0] == doctest::Approx(0.0).epsilon(1e-12));
  CHECK(d0.v[1] == doctest::Approx(std::sqrt(2.0)).epsilon(1e-12));
  const auto d1 = c.derivative(1.0);
  CHECK(d1.v[0] == doctest::Approx(-std::sqrt(2.0)).epsilon(1e-12));
  CHECK(d1.v[1] == doctest::Approx(0.0).epsilon(1e-12));
}

TEST_CASE("NURBS derivative divides by the weight at the parameter") {
  // Non-uniform weights: scaling derivative control points by their own weights
  // (instead of differentiating the homogeneous curve and dividing by w(u))
  // gives different values for this curve.
  const NurbsCurve c(2, {0.0, 0.0, 0.0, 1.0, 1.0, 1.0},
                     {Point3d{0, 0, 0}, Point3d{1, 2, 0}, Point3d{3, 0, 1}},
                     {1.0, 2.0, 0.5});
  REQUIRE(c.is_valid());
  for (double u : {0.2, 0.35, 0.5, 0.8}) {
    const double h = 1e-6;
    const auto fd = (c.eval(u + h) - c.eval(u - h)) * (1.0 / (2.0 * h));
    const auto exact = c.derivative(u);
    CHECK(exact[0] == doctest::Approx(fd[0]).epsilon(1e-5));
    CHECK(exact[1] == doctest::Approx(fd[1]).epsilon(1e-5));
    CHECK(exact[2] == doctest::Approx(fd[2]).epsilon(1e-5));
  }
}

TEST_CASE("NURBS clamps queries to its domain") {
  const NurbsCurve c(2, {0.0, 0.0, 0.0, 1.0, 2.0, 2.0, 2.0},
                     {Point3d{0, 0, 0}, Point3d{1, 1, 0}, Point3d{2, 1, 0}, Point3d{3, 0, 0}},
                     {1.0, 1.0, 1.0, 1.0});
  REQUIRE(c.is_valid());
  const auto [u0, u1] = c.domain();
  CHECK(u0 == 0.0);
  CHECK(u1 == 2.0);
  CHECK(c.eval(-1.0).v[0] == c.eval(u0).v[0]);
  CHECK(c.eval(3.0).v[0] == c.eval(u1).v[0]);
  CHECK(c.eval(3.0).v[2] == c.eval(u1).v[2]);
  CHECK(c.derivative(-1.0).v[0] == c.derivative(u0).v[0]);
  CHECK(c.derivative(99.0).v[1] == c.derivative(u1).v[1]);
}

TEST_CASE("NURBS validity rejects degenerate and unsafe data") {
  const std::vector<Point3d> pts{Point3d{0, 0, 0}, Point3d{1, 0, 0}, Point3d{2, 0, 0}};
  const std::vector<double> w{1.0, 1.0, 1.0};

  const NurbsCurve degree_zero(0, {0.0, 0.0, 0.0, 0.0}, pts, w);
  CHECK_FALSE(degree_zero.is_valid());
  CHECK(degree_zero.domain().first == 0.0);
  CHECK(degree_zero.domain().second == 0.0);

  const NurbsCurve no_control(2, {0.0, 0.0, 0.0, 1.0, 1.0, 1.0}, {}, {});
  CHECK_FALSE(no_control.is_valid());

  const NurbsCurve too_few(2, {0.0, 0.0, 0.0, 1.0, 1.0, 1.0},
                           {Point3d{0, 0, 0}, Point3d{1, 0, 0}}, {1.0, 1.0});
  CHECK_FALSE(too_few.is_valid());

  const NurbsCurve non_monotone(2, {0.0, 0.0, 0.0, 1.0, 0.5, 1.0}, pts, w);
  CHECK_FALSE(non_monotone.is_valid());

  const NurbsCurve zero_weight(2, {0.0, 0.0, 0.0, 1.0, 1.0, 1.0}, pts, {1.0, 0.0, 1.0});
  CHECK_FALSE(zero_weight.is_valid());

  const NurbsCurve negative_weight(2, {0.0, 0.0, 0.0, 1.0, 1.0, 1.0}, pts, {1.0, -2.0, 1.0});
  CHECK_FALSE(negative_weight.is_valid());

  // Degree above the de Boor scratch bound would index past its fixed array.
  std::vector<Point3d> pts17(18, Point3d{0.0, 0.0, 0.0});
  std::vector<double> knots17(36);
  for (std::size_t i = 0; i < knots17.size(); ++i) {
    knots17[i] = static_cast<double>(i) / static_cast<double>(knots17.size() - 1);
  }
  const NurbsCurve too_high(17, knots17, pts17, std::vector<double>(18, 1.0));
  CHECK_FALSE(too_high.is_valid());
  CHECK(too_high.control_count() == 18);

  // The largest supported degree still evaluates.
  std::vector<Point3d> pts16(17, Point3d{0.0, 0.0, 0.0});
  std::vector<double> knots16(34);
  for (std::size_t i = 0; i < knots16.size(); ++i) {
    knots16[i] = static_cast<double>(i) / static_cast<double>(knots16.size() - 1);
  }
  const NurbsCurve max_degree(16, knots16, pts16, std::vector<double>(17, 1.0));
  CHECK(max_degree.is_valid());
  CHECK(max_degree.eval(0.5).v[0] == 0.0);
}
