#include <doctest/doctest.h>

#include <array>
#include <cmath>
#include <cstddef>
#include <limits>
#include <span>

#include <geom/brep/detail/curve_polyline.hpp>
#include <geom/core/error.hpp>
#include <geom/curves/nurbs_curve.hpp>
#include <geom/primitives/point.hpp>

using namespace geom;

namespace {

// Control-point interpolation is exact in exact arithmetic, so a tight relative
// bound is enough to absorb the rounding of constructing knots as i / (n - 1).
bool same_point(const Point3d& a, const Point3d& b) {
  for (int axis = 0; axis < 3; ++axis) {
    if (!std::isfinite(a[axis]) || !std::isfinite(b[axis])) return false;
    if (std::abs(a[axis] - b[axis]) > 1e-12 * (1.0 + std::abs(b[axis]))) return false;
  }
  return true;
}

Point3d midpoint(const Point3d& a, const Point3d& b) {
  return Point3d{0.5 * (a[0] + b[0]), 0.5 * (a[1] + b[1]), 0.5 * (a[2] + b[2])};
}

}  // namespace

TEST_CASE("a polyline curve interpolates every vertex") {
  const std::array<Point3d, 3> points{Point3d{0.0, 0.0, 0.0}, Point3d{1.0, 0.0, 0.0},
                                      Point3d{1.0, 1.0, 0.0}};
  const auto result = make_polyline_curve(points);
  REQUIRE(result.has_value());
  const NurbsCurve& curve = result.value();
  CHECK(curve.is_valid());
  CHECK(curve.degree() == 1);
  CHECK(curve.control_count() == points.size());
  CHECK(curve.domain().first == 0.0);
  CHECK(curve.domain().second == 1.0);
  for (std::size_t i = 0; i < points.size(); ++i) {
    const double t = static_cast<double>(i) / static_cast<double>(points.size() - 1);
    CHECK_MESSAGE(same_point(curve.eval(t), points[i]), "control point " << i);
  }
  const Point3d p = curve.eval(0.5);
  CHECK(p[0] == doctest::Approx(1.0));
  CHECK(p[1] == doctest::Approx(0.0));
  const Point3d q = curve.eval(1.0);
  CHECK(q[1] == doctest::Approx(1.0));
}

TEST_CASE("a two-point polyline curve is the straight segment between them") {
  const std::array<Point3d, 2> points{Point3d{0.0, 0.0, 0.0}, Point3d{4.0, 2.0, 6.0}};
  const auto result = make_polyline_curve(points);
  REQUIRE(result.has_value());
  const NurbsCurve& curve = result.value();
  CHECK(curve.degree() == 1);
  CHECK(curve.control_count() == 2u);
  CHECK(same_point(curve.eval(0.0), points[0]));
  CHECK(same_point(curve.eval(0.25), Point3d{1.0, 0.5, 1.5}));
  CHECK(same_point(curve.eval(0.5), Point3d{2.0, 1.0, 3.0}));
  CHECK(same_point(curve.eval(1.0), points[1]));
}

TEST_CASE("a polyline curve spaces vertices uniformly in parameter") {
  const std::array<Point3d, 5> points{Point3d{0.0, 0.0, 0.0}, Point3d{2.0, 0.0, 0.0},
                                      Point3d{2.0, 3.0, 0.0}, Point3d{5.0, 3.0, 0.0},
                                      Point3d{5.0, 3.0, 4.0}};
  const auto result = make_polyline_curve(points);
  REQUIRE(result.has_value());
  const NurbsCurve& curve = result.value();
  CHECK(curve.control_count() == points.size());
  for (std::size_t i = 0; i < points.size(); ++i) {
    const double t = static_cast<double>(i) / static_cast<double>(points.size() - 1);
    CHECK_MESSAGE(same_point(curve.eval(t), points[i]), "control point " << i);
  }
  // Each segment spans 1 / (n - 1) of the domain, so its midpoint is an exact
  // binary fraction and must land halfway between two consecutive vertices.
  for (std::size_t i = 0; i + 1 < points.size(); ++i) {
    const double t = (static_cast<double>(i) + 0.5) / static_cast<double>(points.size() - 1);
    const Point3d expected = midpoint(points[i], points[i + 1]);
    CHECK_MESSAGE(same_point(curve.eval(t), expected), "segment midpoint " << i);
  }
}

TEST_CASE("a polyline curve needs at least two points") {
  const std::span<const Point3d> empty;
  const auto from_empty = make_polyline_curve(empty);
  REQUIRE_FALSE(from_empty.has_value());
  CHECK(from_empty.error().code == Errc::degenerate);

  const std::array<Point3d, 1> single{Point3d{1.0, 2.0, 3.0}};
  const auto from_single = make_polyline_curve(single);
  REQUIRE_FALSE(from_single.has_value());
  CHECK(from_single.error().code == Errc::degenerate);
}

TEST_CASE("a polyline curve rejects non-finite points") {
  const Point3d origin{0.0, 0.0, 0.0};
  const std::array<Point3d, 3> bad_points{
      Point3d{std::numeric_limits<double>::quiet_NaN(), 0.0, 0.0},
      Point3d{0.0, std::numeric_limits<double>::infinity(), 0.0},
      Point3d{0.0, 0.0, -std::numeric_limits<double>::infinity()}};
  for (const Point3d& bad : bad_points) {
    const std::array<Point3d, 2> first{bad, origin};
    const auto from_first = make_polyline_curve(first);
    REQUIRE_FALSE(from_first.has_value());
    CHECK(from_first.error().code == Errc::non_finite);

    const std::array<Point3d, 2> last{origin, bad};
    const auto from_last = make_polyline_curve(last);
    REQUIRE_FALSE(from_last.has_value());
    CHECK(from_last.error().code == Errc::non_finite);

    const std::array<Point3d, 3> middle{origin, bad, Point3d{1.0, 1.0, 1.0}};
    const auto from_middle = make_polyline_curve(middle);
    REQUIRE_FALSE(from_middle.has_value());
    CHECK(from_middle.error().code == Errc::non_finite);
  }
}
