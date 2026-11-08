#include <doctest/doctest.h>

#include <cmath>
#include <cstddef>
#include <limits>
#include <variant>
#include <vector>

#include <geom/brep/primitives.hpp>
#include <geom/brep/validate.hpp>
#include <geom/core/error.hpp>

namespace {

geom::Vec3d as_vec(geom::Point3d p) noexcept { return geom::Vec3d{p[0], p[1], p[2]}; }

geom::Vec3d unit_of(geom::Vec3d v) {
  const double h = std::hypot(v[0], v[1], v[2]);
  return geom::Vec3d{v[0] / h, v[1] / h, v[2] / h};
}

// Corners of the loop in traversal order: the start vertex of every coedge.
std::vector<geom::Point3d> loop_points(const geom::Brep& b, geom::LoopId loop) {
  std::vector<geom::Point3d> points;
  const geom::CoedgeId start = b.loop(loop).start;
  geom::CoedgeId current = start;
  for (std::size_t step = 0; step <= b.coedge_count(); ++step) {
    const geom::Coedge& c = b.coedge(current);
    const geom::Edge& e = b.edge(c.edge);
    points.push_back(b.vertex(c.reversed ? e.v1 : e.v0).point);
    current = c.next;
    if (current == start || current.is_null()) break;
  }
  return points;
}

}  // namespace

TEST_CASE("box has the expected topology") {
  const auto r = geom::make_box({0, 0, 0}, {2, 3, 4}, 0.0);
  REQUIRE(r.has_value());
  const geom::Brep& b = r.value();
  CHECK(b.vertex_count() == 8);
  CHECK(b.edge_count() == 12);
  CHECK(b.face_count() == 6);
  CHECK(b.coedge_count() == 24);
  CHECK(b.solid_count() == 1);
  CHECK(geom::validate_topology(b).ok());
  CHECK(geom::validate_tolerances(b).ok());
  const geom::SolidId s = b.solid_ids().front();
  CHECK(geom::euler_characteristic(b, s) == 2);
  CHECK(b.shell(b.solid_shells(s).front()).closed);
}

TEST_CASE("box face normals point outward") {
  const auto r = geom::make_box({-1, -2, -3}, {1, 2, 3}, 0.0);
  REQUIRE(r.has_value());
  const geom::Brep& b = r.value();
  for (const geom::Face& f : b.faces()) {
    const auto& plane = std::get<geom::Plane3>(b.geometry().surfaces[f.surface]);
    const double h = std::hypot(plane.normal[0], plane.normal[1], plane.normal[2]);
    const geom::Vec3d unit{plane.normal[0] / h, plane.normal[1] / h,
                           plane.normal[2] / h};
    const geom::CoedgeId c0 = b.loop(f.first_loop).start;
    const geom::Point3d sample = b.vertex(b.edge(b.coedge(c0).edge).v0).point;
    CHECK(geom::dot(unit, geom::Vec3d{sample[0], sample[1], sample[2]}) > 0.0);
  }
}

TEST_CASE("degenerate and non-finite boxes are rejected") {
  CHECK_FALSE(geom::make_box({0, 0, 0}, {0, 1, 1}, 0.0).has_value());
  CHECK_FALSE(geom::make_box({0, 0, 0}, {1, 1, 1}, -1.0).has_value());
  const double nan = std::numeric_limits<double>::quiet_NaN();
  CHECK_FALSE(geom::make_box({0, 0, 0}, {1, nan, 1}, 0.0).has_value());
}

TEST_CASE("box errors carry specific codes") {
  const double nan = std::numeric_limits<double>::quiet_NaN();
  const double inf = std::numeric_limits<double>::infinity();

  const auto equal_corner = geom::make_box({0, 0, 0}, {0, 1, 1}, 0.0);
  REQUIRE_FALSE(equal_corner.has_value());
  CHECK(equal_corner.error().code == geom::Errc::degenerate);

  const auto inverted = geom::make_box({1, 0, 0}, {0, 1, 1}, 0.0);
  REQUIRE_FALSE(inverted.has_value());
  CHECK(inverted.error().code == geom::Errc::degenerate);

  const auto negative_tol = geom::make_box({0, 0, 0}, {1, 1, 1}, -1.0);
  REQUIRE_FALSE(negative_tol.has_value());
  CHECK(negative_tol.error().code == geom::Errc::tolerance_violation);

  const auto nan_hi = geom::make_box({0, 0, 0}, {1, nan, 1}, 0.0);
  REQUIRE_FALSE(nan_hi.has_value());
  CHECK(nan_hi.error().code == geom::Errc::non_finite);

  const auto inf_lo = geom::make_box({0, inf, 0}, {1, 1, 1}, 0.0);
  REQUIRE_FALSE(inf_lo.has_value());
  CHECK(inf_lo.error().code == geom::Errc::non_finite);

  // A non-finite tolerance has no meaningful sign; rejecting it is enough.
  CHECK_FALSE(geom::make_box({0, 0, 0}, {1, 1, 1}, nan).has_value());
}

TEST_CASE("box vertices, edges, faces and the solid carry the tolerance") {
  const double tol = 1e-3;
  const auto r = geom::make_box({0, 0, 0}, {2, 3, 4}, tol);
  REQUIRE(r.has_value());
  const geom::Brep& b = r.value();
  for (const geom::VertexId id : b.vertex_ids()) CHECK(b.vertex(id).tol == tol);
  for (const geom::EdgeId id : b.edge_ids()) CHECK(b.edge(id).tol == tol);
  for (const geom::FaceId id : b.face_ids()) CHECK(b.face(id).tol == tol);
  const geom::SolidId s = b.solid_ids().front();
  CHECK(b.solid(s).tol == tol);
  CHECK(geom::validate_tolerances(b).ok());
}

TEST_CASE("box faces lie on their planes and wind around the outward normal") {
  const auto r = geom::make_box({-1, -2, -3}, {1, 2, 3}, 0.0);
  REQUIRE(r.has_value());
  const geom::Brep& b = r.value();
  REQUIRE(b.face_count() == 6);
  for (const geom::FaceId id : b.face_ids()) {
    const geom::Face& f = b.face(id);
    const auto& plane = std::get<geom::Plane3>(b.geometry().surfaces[f.surface]);
    const geom::Vec3d normal = unit_of(plane.normal);
    const std::vector<geom::Point3d> points = loop_points(b, f.first_loop);
    REQUIRE(points.size() == 4);

    geom::Vec3d area{0.0, 0.0, 0.0};
    for (std::size_t i = 0; i < points.size(); ++i) {
      CHECK(std::fabs(plane.signed_distance(points[i])) <= 1e-12);
      area = area + geom::cross(as_vec(points[i]), as_vec(points[(i + 1) % points.size()]));
    }
    // Right-hand-rule (Newell) normal of the traversal must agree with the
    // stored outward normal, so the loop is counter-clockwise when viewed from
    // outside the box.
    CHECK(geom::dot(unit_of(area), normal) == doctest::Approx(1.0).epsilon(1e-12));
  }
}

TEST_CASE("box edges are unit-speed lines between their vertices") {
  const auto r = geom::make_box({0, 0, 0}, {2, 3, 4}, 0.0);
  REQUIRE(r.has_value());
  const geom::Brep& b = r.value();
  REQUIRE(b.edge_count() == 12);
  for (const geom::EdgeId id : b.edge_ids()) {
    const geom::Edge& e = b.edge(id);
    CHECK(e.t0 == 0.0);
    REQUIRE(e.t1 > 0.0);
    const geom::Point3d p0 = b.vertex(e.v0).point;
    const geom::Point3d p1 = b.vertex(e.v1).point;
    const bool distinct = p0[0] != p1[0] || p0[1] != p1[1] || p0[2] != p1[2];
    CHECK(distinct);

    const geom::Vec3d delta = p1 - p0;
    const double length = std::hypot(delta[0], delta[1], delta[2]);
    CHECK(e.t1 == doctest::Approx(length).epsilon(1e-12));

    const auto& line = std::get<geom::LineCurve3>(b.geometry().curves[e.curve]);
    CHECK(geom::dot(line.direction, line.direction) == doctest::Approx(1.0).epsilon(1e-12));
    const geom::Point3d start = line.eval(e.t0);
    const geom::Point3d end = line.eval(e.t1);
    for (int i = 0; i < 3; ++i) {
      CHECK(start[i] == p0[i]);
      CHECK(end[i] == doctest::Approx(p1[i]).epsilon(1e-12));
    }
  }
}
