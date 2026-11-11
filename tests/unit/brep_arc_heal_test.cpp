#include <geom/brep/detail/arc_heal.hpp>
#include <geom/brep/trim_face.hpp>

#include <doctest/doctest.h>

#include <cmath>
#include <vector>

using geom::v1::CircleCurve3;
using geom::v1::CurvedArcSpan;
using geom::v1::CurvedFragment;
using geom::v1::Point2d;
using geom::v1::Point3d;
using geom::v1::Vec3d;
using geom::v1::detail::collapse_circular_chains;

namespace {

constexpr double kPi = 3.14159265358979323846;

CircleCurve3 unit_circle_z0() {
  CircleCurve3 circle;
  circle.center = Point3d{0.0, 0.0, 0.0};
  circle.x_axis = Vec3d{1.0, 0.0, 0.0};
  circle.y_axis = Vec3d{0.0, 1.0, 0.0};
  circle.radius = 2.0;
  return circle;
}

std::vector<Point3d> circle_points(std::size_t count, double radius, double first, double last) {
  std::vector<Point3d> points;
  points.reserve(count);
  for (std::size_t i = 0; i < count; ++i) {
    const double t = first + (last - first) * static_cast<double>(i) /
                                 static_cast<double>(count - 1u);
    points.push_back(Point3d{radius * std::cos(t), radius * std::sin(t), 0.0});
  }
  return points;
}

bool same_point(const Point3d& a, const Point3d& b) {
  return a[0] == b[0] && a[1] == b[1] && a[2] == b[2];
}

std::vector<Point2d> dummy_uv(std::size_t count) {
  std::vector<Point2d> uv;
  uv.reserve(count);
  for (std::size_t i = 0; i < count; ++i) {
    uv.push_back(Point2d{static_cast<double>(i), 0.0});
  }
  return uv;
}

}  // namespace

TEST_CASE("arc_heal collapses a full circle into two half-turn spans") {
  CurvedFragment fragment;
  fragment.loops.push_back(circle_points(16u, 2.0, 0.0, 2.0 * kPi - 2.0 * kPi / 16.0));
  fragment.loops_uv.push_back(dummy_uv(16u));
  fragment.arc_spans.resize(1u);

  REQUIRE(collapse_circular_chains(fragment, unit_circle_z0(), 1e-9).value());

  REQUIRE(fragment.loops[0].size() == 2u);
  REQUIRE(fragment.arc_spans[0].size() == 2u);
  REQUIRE(fragment.arc_spans[0][0].has_value());
  REQUIRE(fragment.arc_spans[0][1].has_value());
  const double first = fragment.arc_spans[0][0]->t1 - fragment.arc_spans[0][0]->t0;
  const double second = fragment.arc_spans[0][1]->t1 - fragment.arc_spans[0][1]->t0;
  CHECK(std::abs(std::abs(first) - kPi) <= 1e-9);
  CHECK(std::abs(std::abs(second) - kPi) <= 1e-9);
  CHECK(first * second > 0.0);  // same traversal direction
  CHECK(fragment.loops_uv[0].size() == 2u);
}

TEST_CASE("arc_heal leaves the straight radii of a quarter disk alone") {
  CurvedFragment fragment;
  std::vector<Point3d> points = circle_points(7u, 2.0, 0.0, 0.5 * kPi);
  points.push_back(Point3d{0.0, 0.0, 0.0});  // centre
  fragment.loops.push_back(points);
  fragment.loops_uv.push_back(dummy_uv(points.size()));
  fragment.arc_spans.resize(1u);

  REQUIRE(collapse_circular_chains(fragment, unit_circle_z0(), 1e-9).value());

  REQUIRE(fragment.loops[0].size() == 3u);
  REQUIRE(fragment.arc_spans[0].size() == 3u);
  // The retained vertices are the centre, the start and the end of the arc.
  CHECK(same_point(fragment.loops[0][0], Point3d{0.0, 0.0, 0.0}));
  CHECK(same_point(fragment.loops[0][1], points[0]));
  CHECK(same_point(fragment.loops[0][2], points[6]));
  std::size_t engaged = 0u;
  for (const auto& span : fragment.arc_spans[0]) {
    if (span.has_value()) ++engaged;
  }
  CHECK(engaged == 1u);
  const auto& span = fragment.arc_spans[0][1];
  REQUIRE(span.has_value());
  CHECK(std::abs(std::abs(span->t1 - span->t0) - 0.5 * kPi) <= 1e-9);
}

TEST_CASE("arc_heal leaves a non-circular ring unchanged") {
  CurvedFragment fragment;
  fragment.loops.push_back({Point3d{0.0, 0.0, 0.0}, Point3d{1.0, 0.0, 0.0},
                            Point3d{1.0, 1.0, 0.0}, Point3d{0.0, 1.0, 0.0}});
  fragment.loops_uv.push_back(dummy_uv(4u));
  fragment.arc_spans.resize(1u);

  CHECK_FALSE(collapse_circular_chains(fragment, unit_circle_z0(), 1e-9).value());
  CHECK(fragment.loops[0].size() == 4u);
}

TEST_CASE("arc_heal rejects an invalid circle and mismatched uv rings") {
  CurvedFragment fragment;
  fragment.loops.push_back(circle_points(8u, 2.0, 0.0, kPi));
  fragment.loops_uv.push_back(dummy_uv(7u));

  CircleCurve3 bad;
  bad.radius = -1.0;
  const auto invalid = collapse_circular_chains(fragment, bad, 1e-9);
  REQUIRE_FALSE(invalid.has_value());
  CHECK(invalid.error().code == geom::v1::Errc::degenerate);

  const auto mismatch = collapse_circular_chains(fragment, unit_circle_z0(), 1e-9);
  REQUIRE_FALSE(mismatch.has_value());
  CHECK(mismatch.error().code == geom::v1::Errc::invalid_topology);
}
