#include <doctest/doctest.h>

#include <cmath>
#include <cstddef>
#include <limits>
#include <numbers>
#include <variant>
#include <vector>

#include <geom/brep/primitives.hpp>
#include <geom/brep/validate.hpp>
#include <geom/core/error.hpp>
#include <geom/surfaces/analytic.hpp>

namespace {

constexpr double kTwoPi = 2.0 * std::numbers::pi_v<double>;

geom::Vec3d as_vec(geom::Point3d p) noexcept { return geom::Vec3d{p[0], p[1], p[2]}; }

// Coedges of `loop` in traversal order, starting at Loop::start.
std::vector<geom::Coedge> loop_uses(const geom::Brep& b, geom::LoopId loop) {
  std::vector<geom::Coedge> uses;
  const geom::CoedgeId start = b.loop(loop).start;
  geom::CoedgeId current = start;
  for (std::size_t step = 0; step <= b.coedge_count(); ++step) {
    const geom::Coedge& coedge = b.coedge(current);
    uses.push_back(coedge);
    current = coedge.next;
    if (current == start || current.is_null()) break;
  }
  return uses;
}

struct TorusSeams {
  geom::EdgeId u{};
  geom::EdgeId v{};
};

// The u seam is the outer equator (radius major + minor), the v seam is the
// tube circle (radius minor).
TorusSeams torus_seams(const geom::Brep& b, double major_radius, double minor_radius) {
  TorusSeams seams;
  for (const geom::EdgeId id : b.edge_ids()) {
    const auto& circle = std::get<geom::CircleCurve3>(b.geometry().curves[b.edge(id).curve]);
    if (circle.radius == major_radius + minor_radius) {
      seams.u = id;
    } else if (circle.radius == minor_radius) {
      seams.v = id;
    }
  }
  return seams;
}

}  // namespace

TEST_CASE("torus has Euler characteristic zero") {
  const auto r = geom::make_torus({0, 0, 0}, {0, 0, 1}, 3.0, 1.0, 0.0);
  REQUIRE(r.has_value());
  const geom::Brep& b = r.value();
  CHECK(b.vertex_count() == 1);   // 两个接缝的交点
  CHECK(b.edge_count() == 2);     // u 接缝、v 接缝（都是闭合边）
  CHECK(b.face_count() == 1);
  CHECK(b.coedge_count() == 4);
  CHECK(geom::validate_topology(b).ok());
  CHECK(geom::euler_characteristic(b, b.solid_ids().front()) == 0);
  CHECK(b.shell(b.solid_shells(b.solid_ids().front()).front()).closed);
}

TEST_CASE("torus rejects r >= R and non-unit axis") {
  CHECK_FALSE(geom::make_torus({0, 0, 0}, {0, 0, 1}, 1.0, 1.0, 0.0).has_value());
  CHECK_FALSE(geom::make_torus({0, 0, 0}, {0, 0, 2}, 3.0, 1.0, 0.0).has_value());
  CHECK_FALSE(geom::make_torus({0, 0, 0}, {0, 0, 1}, 0.0, 0.1, 0.0).has_value());
}

TEST_CASE("torus seams are full circles meeting at one vertex") {
  const geom::Point3d center{1.0, 2.0, 3.0};
  const geom::Vec3d axis{0.0, 0.0, 1.0};
  const double major = 3.0;
  const double minor = 1.0;
  const auto r = geom::make_torus(center, axis, major, minor, 0.0);
  REQUIRE(r.has_value());
  const geom::Brep& b = r.value();
  REQUIRE(b.vertex_count() == 1);
  REQUIRE(b.edge_count() == 2);
  const geom::VertexId vertex = b.vertex_ids().front();
  const geom::Point3d seam_point = b.vertex(vertex).point;

  const TorusSeams seams = torus_seams(b, major, minor);
  REQUIRE_FALSE(seams.u.is_null());
  REQUIRE_FALSE(seams.v.is_null());

  const geom::Edge& u_edge = b.edge(seams.u);
  const auto& u_circle = std::get<geom::CircleCurve3>(b.geometry().curves[u_edge.curve]);
  CHECK(u_edge.v0 == vertex);
  CHECK(u_edge.v1 == vertex);
  CHECK(u_edge.t0 == 0.0);
  CHECK(u_edge.t1 == doctest::Approx(kTwoPi).epsilon(1e-12));
  CHECK(u_circle.is_valid());
  for (int i = 0; i < 3; ++i) CHECK(u_circle.center[i] == center[i]);
  CHECK(u_circle.radius == major + minor);
  const geom::Vec3d x_axis = u_circle.x_axis;
  const geom::Vec3d y_axis = u_circle.y_axis;
  CHECK(std::hypot(x_axis[0], x_axis[1], x_axis[2]) == doctest::Approx(1.0).epsilon(1e-12));
  CHECK(geom::dot(x_axis, axis) == doctest::Approx(0.0).epsilon(1e-12));
  CHECK(geom::dot(y_axis, axis) == doctest::Approx(0.0).epsilon(1e-12));
  CHECK(geom::dot(x_axis, y_axis) == doctest::Approx(0.0).epsilon(1e-12));
  const geom::Vec3d frame_y = geom::cross(axis, x_axis);
  for (int i = 0; i < 3; ++i) {
    CHECK(frame_y[i] == doctest::Approx(y_axis[i]).epsilon(1e-12));
  }
  for (int i = 0; i < 3; ++i) {
    CHECK(u_circle.eval(0.0)[i] == doctest::Approx(seam_point[i]).epsilon(1e-12));
  }

  const geom::Edge& v_edge = b.edge(seams.v);
  const auto& v_circle = std::get<geom::CircleCurve3>(b.geometry().curves[v_edge.curve]);
  CHECK(v_edge.v0 == vertex);
  CHECK(v_edge.v1 == vertex);
  CHECK(v_edge.t0 == 0.0);
  CHECK(v_edge.t1 == doctest::Approx(kTwoPi).epsilon(1e-12));
  CHECK(v_circle.is_valid());
  const geom::Point3d v_center{center[0] + major * x_axis[0], center[1] + major * x_axis[1],
                               center[2] + major * x_axis[2]};
  for (int i = 0; i < 3; ++i) {
    CHECK(v_circle.center[i] == doctest::Approx(v_center[i]).epsilon(1e-12));
    CHECK(v_circle.x_axis[i] == doctest::Approx(x_axis[i]).epsilon(1e-12));
    CHECK(v_circle.y_axis[i] == doctest::Approx(axis[i]).epsilon(1e-12));
  }
  CHECK(v_circle.radius == minor);
  for (int i = 0; i < 3; ++i) {
    CHECK(v_circle.eval(0.0)[i] == doctest::Approx(seam_point[i]).epsilon(1e-12));
  }
}

TEST_CASE("torus loop uses both seams in both directions") {
  const auto r = geom::make_torus({0, 0, 0}, {0, 0, 1}, 3.0, 1.0, 0.0);
  REQUIRE(r.has_value());
  const geom::Brep& b = r.value();
  REQUIRE(b.face_count() == 1);
  REQUIRE(b.coedge_count() == 4);
  const geom::FaceId face = b.face_ids().front();
  const geom::LoopId loop = b.face(face).first_loop;
  CHECK(b.loop(loop).outer);

  const std::vector<geom::Coedge> uses = loop_uses(b, loop);
  REQUIRE(uses.size() == 4);
  CHECK(b.face(face).reversed == false);
  // The loop starts at the forward use of the u seam (outer equator, R + r).
  const geom::Coedge& first = b.coedge(b.loop(loop).start);
  CHECK(first.edge == uses[0].edge);
  CHECK(first.reversed == false);
  const auto* first_circle =
      std::get_if<geom::CircleCurve3>(&b.geometry().curves[b.edge(first.edge).curve]);
  REQUIRE(first_circle != nullptr);
  CHECK(first_circle->radius == 4.0);
  CHECK(uses[0].reversed == false);
  CHECK(uses[1].reversed == false);
  CHECK(uses[2].reversed == true);
  CHECK(uses[3].reversed == true);
  CHECK(uses[0].edge == uses[2].edge);
  CHECK(uses[1].edge == uses[3].edge);
  CHECK(uses[0].edge != uses[1].edge);
  for (const geom::Coedge& use : uses) {
    CHECK(use.face == face);
    CHECK(use.loop == loop);
  }

  std::vector<int> forward(b.edges().size(), 0);
  std::vector<int> reversed(b.edges().size(), 0);
  for (const geom::CoedgeId id : b.coedge_ids()) {
    const geom::Coedge& coedge = b.coedge(id);
    const std::size_t index = coedge.edge.index();
    if (coedge.reversed) {
      ++reversed[index];
    } else {
      ++forward[index];
    }
  }
  for (const geom::EdgeId id : b.edge_ids()) {
    CHECK(forward[id.index()] == 1);
    CHECK(reversed[id.index()] == 1);
  }
}

TEST_CASE("torus surface matches the seam circles and faces outward") {
  const geom::Point3d center{1.0, 2.0, 3.0};
  const geom::Vec3d axis{0.0, 0.0, 1.0};
  const double major = 3.0;
  const double minor = 1.0;
  const auto r = geom::make_torus(center, axis, major, minor, 0.0);
  REQUIRE(r.has_value());
  const geom::Brep& b = r.value();
  REQUIRE(b.face_count() == 1);
  const geom::Face& face = b.face(b.face_ids().front());
  const auto& torus = std::get<geom::Torus3>(b.geometry().surfaces[face.surface]);
  CHECK(torus.is_valid());
  CHECK(torus.major_radius == major);
  CHECK(torus.minor_radius == minor);
  for (int i = 0; i < 3; ++i) {
    CHECK(torus.center[i] == center[i]);
    CHECK(torus.axis[i] == axis[i]);
  }

  const TorusSeams seams = torus_seams(b, major, minor);
  REQUIRE_FALSE(seams.u.is_null());
  REQUIRE_FALSE(seams.v.is_null());
  const auto& u_circle =
      std::get<geom::CircleCurve3>(b.geometry().curves[b.edge(seams.u).curve]);
  const auto& v_circle =
      std::get<geom::CircleCurve3>(b.geometry().curves[b.edge(seams.v).curve]);

  const double us[3] = {0.0, 0.7, 2.5};
  for (const double u : us) {
    for (int i = 0; i < 3; ++i) {
      CHECK(torus.eval(u, 0.0)[i] == doctest::Approx(u_circle.eval(u)[i]).epsilon(1e-12));
    }
    CHECK(geom::dot(as_vec(torus.eval(u, 0.0)) - as_vec(center), axis) ==
          doctest::Approx(0.0).epsilon(1e-12));
  }
  const double vs[3] = {0.0, 1.1, 2.4};
  for (const double v : vs) {
    for (int i = 0; i < 3; ++i) {
      CHECK(torus.eval(0.0, v)[i] == doctest::Approx(v_circle.eval(v)[i]).epsilon(1e-12));
    }
  }

  // The stored normal is the outward one: +x on the outer equator, +axis on the
  // top, -x on the inner equator and -axis on the bottom.
  const geom::Vec3d x_axis = u_circle.x_axis;
  const geom::Vec3d outer = torus.normal(0.0, 0.0);
  const geom::Vec3d top = torus.normal(0.0, std::numbers::pi_v<double> / 2.0);
  const geom::Vec3d inner = torus.normal(0.0, std::numbers::pi_v<double>);
  const geom::Vec3d bottom = torus.normal(0.0, -std::numbers::pi_v<double> / 2.0);
  for (int i = 0; i < 3; ++i) {
    CHECK(outer[i] == doctest::Approx(x_axis[i]).epsilon(1e-12));
    CHECK(top[i] == doctest::Approx(axis[i]).epsilon(1e-12));
    CHECK(inner[i] == doctest::Approx(-x_axis[i]).epsilon(1e-12));
    CHECK(bottom[i] == doctest::Approx(-axis[i]).epsilon(1e-12));
  }
}

TEST_CASE("torus carries the tolerance") {
  const double tol = 1e-3;
  const auto r = geom::make_torus({0, 0, 0}, {0, 0, 1}, 3.0, 1.0, tol);
  REQUIRE(r.has_value());
  const geom::Brep& b = r.value();
  for (const geom::VertexId id : b.vertex_ids()) CHECK(b.vertex(id).tol == tol);
  for (const geom::EdgeId id : b.edge_ids()) CHECK(b.edge(id).tol == tol);
  for (const geom::FaceId id : b.face_ids()) CHECK(b.face(id).tol == tol);
  CHECK(b.solid(b.solid_ids().front()).tol == tol);
  CHECK(geom::validate_tolerances(b).ok());
  CHECK(b.geometry().pcurves.empty());
  for (const geom::CoedgeId id : b.coedge_ids()) {
    CHECK(b.coedge(id).pcurve == geom::kNoGeometry);
  }
}

TEST_CASE("torus errors carry specific codes") {
  const double nan = std::numeric_limits<double>::quiet_NaN();
  const double inf = std::numeric_limits<double>::infinity();

  const auto zero_major = geom::make_torus({0, 0, 0}, {0, 0, 1}, 0.0, 0.1, 0.0);
  REQUIRE_FALSE(zero_major.has_value());
  CHECK(zero_major.error().code == geom::Errc::degenerate);

  const auto negative_major = geom::make_torus({0, 0, 0}, {0, 0, 1}, -1.0, 0.1, 0.0);
  REQUIRE_FALSE(negative_major.has_value());
  CHECK(negative_major.error().code == geom::Errc::degenerate);

  const auto zero_minor = geom::make_torus({0, 0, 0}, {0, 0, 1}, 3.0, 0.0, 0.0);
  REQUIRE_FALSE(zero_minor.has_value());
  CHECK(zero_minor.error().code == geom::Errc::degenerate);

  const auto negative_minor = geom::make_torus({0, 0, 0}, {0, 0, 1}, 3.0, -0.5, 0.0);
  REQUIRE_FALSE(negative_minor.has_value());
  CHECK(negative_minor.error().code == geom::Errc::degenerate);

  const auto equal_radii = geom::make_torus({0, 0, 0}, {0, 0, 1}, 2.0, 2.0, 0.0);
  REQUIRE_FALSE(equal_radii.has_value());
  CHECK(equal_radii.error().code == geom::Errc::degenerate);

  const auto spindle = geom::make_torus({0, 0, 0}, {0, 0, 1}, 1.0, 2.0, 0.0);
  REQUIRE_FALSE(spindle.has_value());
  CHECK(spindle.error().code == geom::Errc::degenerate);

  const auto zero_axis = geom::make_torus({0, 0, 0}, {0, 0, 0}, 3.0, 1.0, 0.0);
  REQUIRE_FALSE(zero_axis.has_value());
  CHECK(zero_axis.error().code == geom::Errc::degenerate);

  const auto long_axis = geom::make_torus({0, 0, 0}, {0, 0, 2}, 3.0, 1.0, 0.0);
  REQUIRE_FALSE(long_axis.has_value());
  CHECK(long_axis.error().code == geom::Errc::degenerate);

  // R and r are finite and r < R, but the outer seam radius R + r overflows.
  const auto overflow = geom::make_torus({0, 0, 0}, {0, 0, 1}, 1.5e308, 1e308, 0.0);
  REQUIRE_FALSE(overflow.has_value());
  CHECK(overflow.error().code == geom::Errc::degenerate);

  const auto nan_center = geom::make_torus({nan, 0, 0}, {0, 0, 1}, 3.0, 1.0, 0.0);
  REQUIRE_FALSE(nan_center.has_value());
  CHECK(nan_center.error().code == geom::Errc::non_finite);

  const auto inf_axis = geom::make_torus({0, 0, 0}, {0, inf, 0}, 3.0, 1.0, 0.0);
  REQUIRE_FALSE(inf_axis.has_value());
  CHECK(inf_axis.error().code == geom::Errc::non_finite);

  const auto inf_major = geom::make_torus({0, 0, 0}, {0, 0, 1}, inf, 1.0, 0.0);
  REQUIRE_FALSE(inf_major.has_value());
  CHECK(inf_major.error().code == geom::Errc::non_finite);

  const auto inf_minor = geom::make_torus({0, 0, 0}, {0, 0, 1}, 3.0, inf, 0.0);
  REQUIRE_FALSE(inf_minor.has_value());
  CHECK(inf_minor.error().code == geom::Errc::non_finite);

  const auto negative_tol = geom::make_torus({0, 0, 0}, {0, 0, 1}, 3.0, 1.0, -1.0);
  REQUIRE_FALSE(negative_tol.has_value());
  CHECK(negative_tol.error().code == geom::Errc::tolerance_violation);

  const auto nan_tol = geom::make_torus({0, 0, 0}, {0, 0, 1}, 3.0, 1.0, nan);
  REQUIRE_FALSE(nan_tol.has_value());
  CHECK(nan_tol.error().code == geom::Errc::tolerance_violation);
}

TEST_CASE("torus accepts an axis within the analytic unit tolerance") {
  const auto accepted = geom::make_torus({0, 0, 0}, {0, 0, 1.0 + 5e-13}, 3.0, 1.0, 0.0);
  REQUIRE(accepted.has_value());
  CHECK(geom::validate_topology(accepted.value()).ok());
  CHECK(geom::euler_characteristic(accepted.value(),
                                   accepted.value().solid_ids().front()) == 0);

  const auto rejected = geom::make_torus({0, 0, 0}, {0, 0, 1.0 + 3e-12}, 3.0, 1.0, 0.0);
  REQUIRE_FALSE(rejected.has_value());
  CHECK(rejected.error().code == geom::Errc::degenerate);
}

TEST_CASE("torus frame adapts to a tilted axis") {
  const geom::Point3d center{1.0, 2.0, 3.0};
  const double diagonal = 1.0 / std::sqrt(3.0);
  const geom::Vec3d axis{diagonal, diagonal, diagonal};
  const double major = 3.0;
  const double minor = 1.0;
  const auto r = geom::make_torus(center, axis, major, minor, 0.0);
  REQUIRE(r.has_value());
  const geom::Brep& b = r.value();
  CHECK(geom::validate_topology(b).ok());
  CHECK(geom::euler_characteristic(b, b.solid_ids().front()) == 0);

  const auto& torus = std::get<geom::Torus3>(
      b.geometry().surfaces[b.face(b.face_ids().front()).surface]);
  const geom::Vec3d x = torus.x_axis;
  CHECK(std::hypot(x[0], x[1], x[2]) == doctest::Approx(1.0).epsilon(1e-12));
  CHECK(geom::dot(x, axis) == doctest::Approx(0.0).epsilon(1e-12));
  for (int i = 0; i < 3; ++i) CHECK(torus.axis[i] == axis[i]);

  const TorusSeams seams = torus_seams(b, major, minor);
  REQUIRE_FALSE(seams.u.is_null());
  REQUIRE_FALSE(seams.v.is_null());
  const auto& u_circle =
      std::get<geom::CircleCurve3>(b.geometry().curves[b.edge(seams.u).curve]);
  const auto& v_circle =
      std::get<geom::CircleCurve3>(b.geometry().curves[b.edge(seams.v).curve]);
  const geom::Vec3d frame_y = geom::cross(axis, x);
  for (int i = 0; i < 3; ++i) {
    CHECK(u_circle.x_axis[i] == doctest::Approx(x[i]).epsilon(1e-12));
    CHECK(u_circle.y_axis[i] == doctest::Approx(frame_y[i]).epsilon(1e-12));
    CHECK(v_circle.y_axis[i] == doctest::Approx(axis[i]).epsilon(1e-12));
  }
  const geom::Point3d expected{center[0] + (major + minor) * x[0],
                               center[1] + (major + minor) * x[1],
                               center[2] + (major + minor) * x[2]};
  const geom::Point3d seam_point = b.vertex(b.vertex_ids().front()).point;
  for (int i = 0; i < 3; ++i) {
    CHECK(u_circle.eval(0.0)[i] == doctest::Approx(expected[i]).epsilon(1e-12));
    CHECK(v_circle.eval(0.0)[i] == doctest::Approx(seam_point[i]).epsilon(1e-12));
  }
}
