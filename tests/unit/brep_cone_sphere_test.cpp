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
#include <geom/curves/analytic_curve.hpp>
#include <geom/surfaces/analytic.hpp>
#include <geom/surfaces/plane.hpp>

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

// Samples the loop's traversal in order. Closed curves are sampled at several
// parameters per coedge so the Newell sum sees the actual polygon.
std::vector<geom::Point3d> loop_polyline(const geom::Brep& b, geom::LoopId loop,
                                         int samples_per_edge = 32) {
  std::vector<geom::Point3d> points;
  for (const geom::Coedge& use : loop_uses(b, loop)) {
    const geom::Edge& edge = b.edge(use.edge);
    if (edge.curve == geom::kNoGeometry) {
      points.push_back(b.vertex(use.reversed ? edge.v1 : edge.v0).point);
      continue;
    }
    const geom::CurveGeometry& curve = b.geometry().curves[edge.curve];
    for (int i = 0; i < samples_per_edge; ++i) {
      const double s = static_cast<double>(i) / samples_per_edge;
      const double t = use.reversed ? edge.t1 - s * (edge.t1 - edge.t0)
                                    : edge.t0 + s * (edge.t1 - edge.t0);
      points.push_back(std::visit([t](const auto& c) { return c.eval(t); }, curve));
    }
  }
  return points;
}

// Newell normal of a closed polygon: it points along the polygon's
// counter-clockwise normal, so a planar loop wound CCW about its outward face
// normal gives a positive dot product with that normal.
geom::Vec3d newell_normal(const std::vector<geom::Point3d>& points) {
  double x = 0.0;
  double y = 0.0;
  double z = 0.0;
  for (std::size_t i = 0; i < points.size(); ++i) {
    const geom::Point3d& a = points[i];
    const geom::Point3d& c = points[(i + 1) % points.size()];
    x += a[1] * c[2] - a[2] * c[1];
    y += a[2] * c[0] - a[0] * c[2];
    z += a[0] * c[1] - a[1] * c[0];
  }
  return geom::Vec3d{x, y, z};
}

}  // namespace

TEST_CASE("cone topology shares the base circle and the seam") {
  const auto r = geom::make_cone({0, 0, 0}, {0, 0, 1}, 2.0, 3.0, 0.0);
  REQUIRE(r.has_value());
  const geom::Brep& b = r.value();
  CHECK(b.vertex_count() == 2);   // 底接缝点、锥顶
  CHECK(b.edge_count() == 2);     // 底圆、接缝
  CHECK(b.face_count() == 2);     // 底盖、侧面
  CHECK(b.coedge_count() == 4);
  CHECK(geom::validate_topology(b).ok());
  CHECK(geom::euler_characteristic(b, b.solid_ids().front()) == 2);
}

TEST_CASE("sphere is a single face with a two-sided seam") {
  const auto r = geom::make_sphere({0, 0, 0}, 3.0, 0.0);
  REQUIRE(r.has_value());
  const geom::Brep& b = r.value();
  CHECK(b.vertex_count() == 2);   // 南北极
  CHECK(b.edge_count() == 1);     // 子午线接缝
  CHECK(b.face_count() == 1);
  CHECK(b.coedge_count() == 2);
  CHECK(geom::validate_topology(b).ok());
  CHECK(geom::euler_characteristic(b, b.solid_ids().front()) == 2);
}

TEST_CASE("cone and sphere reject degenerate parameters") {
  CHECK_FALSE(geom::make_cone({0, 0, 0}, {0, 0, 1}, 0.0, 1.0, 0.0).has_value());
  CHECK_FALSE(geom::make_cone({0, 0, 0}, {0, 0, 1}, 1.0, 0.0, 0.0).has_value());
  CHECK_FALSE(geom::make_cone({0, 0, 0}, {0, 0, 2}, 1.0, 1.0, 0.0).has_value());
  CHECK_FALSE(geom::make_sphere({0, 0, 0}, 0.0, 0.0).has_value());
  CHECK_FALSE(geom::make_sphere({0, 0, 0}, -1.0, 0.0).has_value());
}

TEST_CASE("cone surface and seam follow the derived half angle") {
  const geom::Point3d base{1.0, 2.0, 3.0};
  const geom::Vec3d axis{0.0, 0.0, 1.0};
  const double radius = 2.0;
  const double height = 3.0;
  const auto r = geom::make_cone(base, axis, radius, height, 0.0);
  REQUIRE(r.has_value());
  const geom::Brep& b = r.value();
  REQUIRE(b.face_count() == 2);
  REQUIRE(b.edge_count() == 2);
  REQUIRE(b.vertex_count() == 2);

  const geom::Cone3* cone = nullptr;
  for (const geom::FaceId id : b.face_ids()) {
    const geom::SurfaceGeometry& surface = b.geometry().surfaces[b.face(id).surface];
    if (const auto* candidate = std::get_if<geom::Cone3>(&surface)) cone = candidate;
  }
  REQUIRE(cone != nullptr);
  CHECK(cone->is_valid());
  CHECK(cone->half_angle == doctest::Approx(std::atan2(radius, height)).epsilon(1e-15));
  CHECK(cone->radius == radius);
  for (int i = 0; i < 3; ++i) {
    CHECK(cone->origin[i] == base[i]);
    CHECK(cone->axis[i] == axis[i]);
  }
  const geom::Vec3d x_axis = cone->x_axis;
  CHECK(std::hypot(x_axis[0], x_axis[1], x_axis[2]) == doctest::Approx(1.0).epsilon(1e-12));
  CHECK(geom::dot(x_axis, axis) == doctest::Approx(0.0).epsilon(1e-12));

  geom::EdgeId base_circle{};
  geom::EdgeId seam{};
  for (const geom::EdgeId id : b.edge_ids()) {
    const geom::Edge& edge = b.edge(id);
    if (edge.v0 == edge.v1) {
      base_circle = id;
    } else {
      seam = id;
    }
  }
  REQUIRE_FALSE(base_circle.is_null());
  REQUIRE_FALSE(seam.is_null());
  CHECK(base_circle != seam);

  const geom::Edge& circle_edge = b.edge(base_circle);
  const auto& circle = std::get<geom::CircleCurve3>(b.geometry().curves[circle_edge.curve]);
  CHECK(circle.is_valid());
  CHECK(circle_edge.t0 == 0.0);
  CHECK(circle_edge.t1 == doctest::Approx(kTwoPi).epsilon(1e-12));
  CHECK(circle.radius == radius);
  for (int i = 0; i < 3; ++i) CHECK(circle.center[i] == base[i]);
  for (int i = 0; i < 3; ++i) CHECK(circle.x_axis[i] == x_axis[i]);
  const geom::Vec3d frame_y = geom::cross(axis, x_axis);
  for (int i = 0; i < 3; ++i) {
    CHECK(circle.y_axis[i] == doctest::Approx(frame_y[i]).epsilon(1e-12));
  }

  const geom::Edge& seam_edge = b.edge(seam);
  const auto& seam_line = std::get<geom::LineCurve3>(b.geometry().curves[seam_edge.curve]);
  const geom::Point3d apex{base[0] + height * axis[0], base[1] + height * axis[1],
                           base[2] + height * axis[2]};
  const geom::Point3d p0{base[0] + radius * x_axis[0], base[1] + radius * x_axis[1],
                         base[2] + radius * x_axis[2]};
  const double slant = std::hypot(radius, height);
  CHECK(seam_edge.v0 == circle_edge.v0);
  CHECK(seam_edge.t0 == 0.0);
  CHECK(seam_edge.t1 == doctest::Approx(slant).epsilon(1e-12));
  const geom::Point3d v0 = b.vertex(seam_edge.v0).point;
  const geom::Point3d v1 = b.vertex(seam_edge.v1).point;
  for (int i = 0; i < 3; ++i) {
    CHECK(v0[i] == doctest::Approx(p0[i]).epsilon(1e-12));
    CHECK(v1[i] == doctest::Approx(apex[i]).epsilon(1e-12));
    CHECK(seam_line.origin[i] == doctest::Approx(p0[i]).epsilon(1e-12));
  }
  // The seam is the slant generatrix of the Cone3 frame.
  for (int i = 0; i < 3; ++i) {
    const double direction = (height * axis[i] - radius * x_axis[i]) / slant;
    CHECK(seam_line.direction[i] == doctest::Approx(direction).epsilon(1e-12));
  }
  for (int i = 0; i < 3; ++i) {
    CHECK(seam_line.eval(seam_edge.t1)[i] == doctest::Approx(apex[i]).epsilon(1e-12));
    CHECK(cone->eval(0.0, 0.0)[i] == doctest::Approx(p0[i]).epsilon(1e-12));
    CHECK(cone->eval(0.0, height)[i] == doctest::Approx(apex[i]).epsilon(1e-12));
  }
  // Arc length along the seam grows as v / cos(half_angle).
  const double cos_half = std::cos(cone->half_angle);
  for (const double v : {0.0, 0.5 * height, height}) {
    const geom::Point3d surface_point = cone->eval(0.0, v);
    const geom::Point3d seam_point = seam_line.eval(v / cos_half);
    for (int i = 0; i < 3; ++i) {
      CHECK(surface_point[i] == doctest::Approx(seam_point[i]).epsilon(1e-12));
    }
  }
}

TEST_CASE("cone cap and lateral loops bound an outward solid") {
  const geom::Vec3d axis{0.0, 0.0, 1.0};
  const auto r = geom::make_cone({0, 0, 0}, axis, 2.0, 3.0, 0.0);
  REQUIRE(r.has_value());
  const geom::Brep& b = r.value();

  geom::FaceId cap{};
  geom::FaceId lateral{};
  const geom::Plane3* cap_plane = nullptr;
  const geom::Cone3* cone = nullptr;
  for (const geom::FaceId id : b.face_ids()) {
    const geom::SurfaceGeometry& surface = b.geometry().surfaces[b.face(id).surface];
    if (const auto* plane = std::get_if<geom::Plane3>(&surface)) {
      cap = id;
      cap_plane = plane;
    } else if (const auto* candidate = std::get_if<geom::Cone3>(&surface)) {
      lateral = id;
      cone = candidate;
    }
  }
  REQUIRE_FALSE(cap.is_null());
  REQUIRE_FALSE(lateral.is_null());
  REQUIRE(cap_plane != nullptr);
  REQUIRE(cone != nullptr);
  CHECK(geom::dot(cap_plane->normal, axis) == doctest::Approx(-1.0).epsilon(1e-12));
  CHECK(cap_plane->signed_distance({0.0, 0.0, 0.0}) == doctest::Approx(0.0).epsilon(1e-12));
  // The solid interior is on the negative side of the base cap.
  CHECK(cap_plane->signed_distance({0.0, 0.0, 1.5}) < 0.0);

  const std::vector<geom::Coedge> cap_uses = loop_uses(b, b.face(cap).first_loop);
  REQUIRE(cap_uses.size() == 1);
  CHECK(cap_uses[0].reversed == true);
  const std::vector<geom::Coedge> lateral_uses = loop_uses(b, b.face(lateral).first_loop);
  REQUIRE(lateral_uses.size() == 3);
  CHECK(lateral_uses[0].edge == cap_uses[0].edge);
  CHECK(lateral_uses[0].reversed == false);
  CHECK(lateral_uses[1].reversed == false);
  CHECK(lateral_uses[2].reversed == true);
  CHECK(lateral_uses[1].edge == lateral_uses[2].edge);  // the seam, up then down
  CHECK(lateral_uses[1].edge != lateral_uses[0].edge);

  std::vector<int> forward(b.edges().size(), 0);
  std::vector<int> reversed(b.edges().size(), 0);
  for (const geom::CoedgeId id : b.coedge_ids()) {
    const geom::Coedge& coedge = b.coedge(id);
    if (coedge.reversed) {
      ++reversed[coedge.edge.index()];
    } else {
      ++forward[coedge.edge.index()];
    }
  }
  for (const geom::EdgeId id : b.edge_ids()) {
    CHECK(forward[id.index()] == 1);
    CHECK(reversed[id.index()] == 1);
  }

  // The cap loop winds counter-clockwise about the outward cap normal.
  const geom::Vec3d cap_newell = newell_normal(loop_polyline(b, b.face(cap).first_loop));
  CHECK(geom::dot(cap_newell, cap_plane->normal) > 0.0);

  // The lateral parameterization normal agrees with cross(du, dv) and points
  // away from the axis and toward the apex.
  const geom::Vec3d du = cone->du(0.0, 1.0);
  const geom::Vec3d dv = cone->dv(0.0, 1.0);
  const geom::Vec3d normal = cone->normal(0.0, 1.0);
  CHECK(std::hypot(normal[0], normal[1], normal[2]) == doctest::Approx(1.0).epsilon(1e-12));
  CHECK(geom::dot(normal, geom::cross(du, dv)) > 0.0);
  CHECK(normal[0] > 0.0);
  CHECK(normal[2] > 0.0);
}

TEST_CASE("sphere seam is a half great circle") {
  const geom::Point3d center{1.0, 2.0, 3.0};
  const double radius = 3.0;
  const auto r = geom::make_sphere(center, radius, 0.0);
  REQUIRE(r.has_value());
  const geom::Brep& b = r.value();
  REQUIRE(b.vertex_count() == 2);
  REQUIRE(b.edge_count() == 1);
  REQUIRE(b.face_count() == 1);
  REQUIRE(b.coedge_count() == 2);

  const geom::FaceId face = b.face_ids().front();
  const geom::SurfaceGeometry& surface = b.geometry().surfaces[b.face(face).surface];
  const auto* sphere = std::get_if<geom::Sphere3>(&surface);
  REQUIRE(sphere != nullptr);
  CHECK(sphere->is_valid());
  CHECK(sphere->radius == radius);
  for (int i = 0; i < 3; ++i) CHECK(sphere->center[i] == center[i]);
  CHECK(sphere->axis[2] == 1.0);
  CHECK(sphere->x_axis[0] == 1.0);

  const geom::EdgeId seam = b.edge_ids().front();
  const geom::Edge& edge = b.edge(seam);
  const auto& circle = std::get<geom::CircleCurve3>(b.geometry().curves[edge.curve]);
  CHECK(circle.is_valid());
  CHECK(circle.radius == radius);
  for (int i = 0; i < 3; ++i) CHECK(circle.center[i] == center[i]);
  CHECK(circle.x_axis[0] == 1.0);
  CHECK(circle.y_axis[2] == 1.0);
  const double half_pi = std::numbers::pi_v<double> / 2.0;
  CHECK(edge.t0 == doctest::Approx(-half_pi).epsilon(1e-15));
  CHECK(edge.t1 == doctest::Approx(half_pi).epsilon(1e-15));

  const geom::Point3d south{center[0], center[1], center[2] - radius};
  const geom::Point3d north{center[0], center[1], center[2] + radius};
  const geom::Point3d v0 = b.vertex(edge.v0).point;
  const geom::Point3d v1 = b.vertex(edge.v1).point;
  for (int i = 0; i < 3; ++i) {
    CHECK(v0[i] == doctest::Approx(south[i]).epsilon(1e-12));
    CHECK(v1[i] == doctest::Approx(north[i]).epsilon(1e-12));
  }
  CHECK(edge.v0 != edge.v1);
  const geom::Point3d start = circle.eval(edge.t0);
  const geom::Point3d end = circle.eval(edge.t1);
  for (int i = 0; i < 3; ++i) {
    CHECK(start[i] == doctest::Approx(v0[i]).epsilon(1e-12));
    CHECK(end[i] == doctest::Approx(v1[i]).epsilon(1e-12));
  }

  const std::vector<geom::Coedge> uses = loop_uses(b, b.face(face).first_loop);
  REQUIRE(uses.size() == 2);
  CHECK(uses[0].edge == seam);
  CHECK(uses[1].edge == seam);
  CHECK(uses[0].reversed == false);
  CHECK(uses[1].reversed == true);

  // Sphere3::eval(0, v) is exactly the seam curve.
  for (const double t : {-half_pi, -0.5, 0.0, 0.5, half_pi}) {
    const geom::Point3d p = sphere->eval(0.0, t);
    const geom::Point3d q = circle.eval(t);
    for (int i = 0; i < 3; ++i) {
      CHECK(p[i] == doctest::Approx(q[i]).epsilon(1e-12));
    }
  }

  // Outward unit normals, including the pole limits.
  const geom::Vec3d equator = sphere->normal(0.0, 0.0);
  for (int i = 0; i < 3; ++i) {
    CHECK(equator[i] == doctest::Approx(sphere->x_axis[i]).epsilon(1e-12));
  }
  const geom::Vec3d north_normal = sphere->normal(0.0, half_pi);
  const geom::Vec3d south_normal = sphere->normal(0.0, -half_pi);
  CHECK(std::hypot(north_normal[0], north_normal[1], north_normal[2]) ==
        doctest::Approx(1.0).epsilon(1e-12));
  for (int i = 0; i < 2; ++i) {
    CHECK(north_normal[i] == doctest::Approx(0.0).epsilon(1e-15));
    CHECK(south_normal[i] == doctest::Approx(0.0).epsilon(1e-15));
  }
  CHECK(north_normal[2] == doctest::Approx(1.0).epsilon(1e-12));
  CHECK(south_normal[2] == doctest::Approx(-1.0).epsilon(1e-12));

  const geom::Point3d point = sphere->eval(1.1, 0.3);
  const geom::Vec3d outward = as_vec(point) - as_vec(center);
  CHECK(geom::dot(sphere->normal(1.1, 0.3), outward) > 0.0);
}

TEST_CASE("cone and sphere reject non-finite, invalid tolerance and collapsing geometry") {
  const double nan = std::numeric_limits<double>::quiet_NaN();
  const double inf = std::numeric_limits<double>::infinity();

  const auto nan_base = geom::make_cone({nan, 0, 0}, {0, 0, 1}, 1.0, 1.0, 0.0);
  REQUIRE_FALSE(nan_base.has_value());
  CHECK(nan_base.error().code == geom::Errc::non_finite);
  const auto inf_axis = geom::make_cone({0, 0, 0}, {0, inf, 0}, 1.0, 1.0, 0.0);
  REQUIRE_FALSE(inf_axis.has_value());
  CHECK(inf_axis.error().code == geom::Errc::non_finite);
  const auto nan_radius = geom::make_cone({0, 0, 0}, {0, 0, 1}, nan, 1.0, 0.0);
  REQUIRE_FALSE(nan_radius.has_value());
  CHECK(nan_radius.error().code == geom::Errc::non_finite);
  const auto inf_height = geom::make_cone({0, 0, 0}, {0, 0, 1}, 1.0, inf, 0.0);
  REQUIRE_FALSE(inf_height.has_value());
  CHECK(inf_height.error().code == geom::Errc::non_finite);

  const auto negative_tol = geom::make_cone({0, 0, 0}, {0, 0, 1}, 1.0, 1.0, -1.0);
  REQUIRE_FALSE(negative_tol.has_value());
  CHECK(negative_tol.error().code == geom::Errc::tolerance_violation);
  const auto nan_tol = geom::make_cone({0, 0, 0}, {0, 0, 1}, 1.0, 1.0, nan);
  REQUIRE_FALSE(nan_tol.has_value());
  CHECK(nan_tol.error().code == geom::Errc::tolerance_violation);

  const auto zero_axis = geom::make_cone({0, 0, 0}, {0, 0, 0}, 1.0, 1.0, 0.0);
  REQUIRE_FALSE(zero_axis.has_value());
  CHECK(zero_axis.error().code == geom::Errc::degenerate);
  const auto long_axis = geom::make_cone({0, 0, 0}, {0, 0, 2}, 1.0, 1.0, 0.0);
  REQUIRE_FALSE(long_axis.has_value());
  CHECK(long_axis.error().code == geom::Errc::degenerate);
  const auto zero_radius = geom::make_cone({0, 0, 0}, {0, 0, 1}, 0.0, 1.0, 0.0);
  REQUIRE_FALSE(zero_radius.has_value());
  CHECK(zero_radius.error().code == geom::Errc::degenerate);
  const auto negative_radius = geom::make_cone({0, 0, 0}, {0, 0, 1}, -1.0, 1.0, 0.0);
  REQUIRE_FALSE(negative_radius.has_value());
  CHECK(negative_radius.error().code == geom::Errc::degenerate);
  const auto zero_height = geom::make_cone({0, 0, 0}, {0, 0, 1}, 1.0, 0.0, 0.0);
  REQUIRE_FALSE(zero_height.has_value());
  CHECK(zero_height.error().code == geom::Errc::degenerate);
  const auto negative_height = geom::make_cone({0, 0, 0}, {0, 0, 1}, 1.0, -1.0, 0.0);
  REQUIRE_FALSE(negative_height.has_value());
  CHECK(negative_height.error().code == geom::Errc::degenerate);

  // atan2 underflows to zero, so the derived half angle is not representable.
  const auto flat_cone = geom::make_cone({0, 0, 0}, {0, 0, 1}, 1e-300, 1e308, 0.0);
  REQUIRE_FALSE(flat_cone.has_value());
  CHECK(flat_cone.error().code == geom::Errc::degenerate);
  // The seam point base + radius * x overflows.
  const auto wide_seam = geom::make_cone({1e308, 0, 0}, {0, 0, 1}, 1e308, 1.0, 0.0);
  REQUIRE_FALSE(wide_seam.has_value());
  CHECK(wide_seam.error().code == geom::Errc::degenerate);
  // The apex base + height * axis overflows.
  const auto tall_apex = geom::make_cone({1e308, 0, 0}, {1, 0, 0}, 1.0, 1e308, 0.0);
  REQUIRE_FALSE(tall_apex.has_value());
  CHECK(tall_apex.error().code == geom::Errc::degenerate);
  // Both seam endpoints round to the same point.
  const auto zero_seam = geom::make_cone({1e308, 0, 1e308}, {0, 0, 1}, 1e-300, 1e-300, 0.0);
  REQUIRE_FALSE(zero_seam.has_value());
  CHECK(zero_seam.error().code == geom::Errc::degenerate);

  const auto nan_center = geom::make_sphere({nan, 0, 0}, 1.0, 0.0);
  REQUIRE_FALSE(nan_center.has_value());
  CHECK(nan_center.error().code == geom::Errc::non_finite);
  const auto inf_sphere_radius = geom::make_sphere({0, 0, 0}, inf, 0.0);
  REQUIRE_FALSE(inf_sphere_radius.has_value());
  CHECK(inf_sphere_radius.error().code == geom::Errc::non_finite);
  const auto negative_sphere_tol = geom::make_sphere({0, 0, 0}, 1.0, -1.0);
  REQUIRE_FALSE(negative_sphere_tol.has_value());
  CHECK(negative_sphere_tol.error().code == geom::Errc::tolerance_violation);
  const auto nan_sphere_tol = geom::make_sphere({0, 0, 0}, 1.0, nan);
  REQUIRE_FALSE(nan_sphere_tol.has_value());
  CHECK(nan_sphere_tol.error().code == geom::Errc::tolerance_violation);
  const auto zero_sphere_radius = geom::make_sphere({0, 0, 0}, 0.0, 0.0);
  REQUIRE_FALSE(zero_sphere_radius.has_value());
  CHECK(zero_sphere_radius.error().code == geom::Errc::degenerate);
  const auto negative_sphere_radius = geom::make_sphere({0, 0, 0}, -1.0, 0.0);
  REQUIRE_FALSE(negative_sphere_radius.has_value());
  CHECK(negative_sphere_radius.error().code == geom::Errc::degenerate);
  // Both poles round to the same point.
  const auto flat_sphere = geom::make_sphere({0, 0, 1e308}, 1e-300, 0.0);
  REQUIRE_FALSE(flat_sphere.has_value());
  CHECK(flat_sphere.error().code == geom::Errc::degenerate);
  // The north pole center + radius * axis overflows.
  const auto tall_sphere = geom::make_sphere({0, 0, 1e307}, 1.7e308, 0.0);
  REQUIRE_FALSE(tall_sphere.has_value());
  CHECK(tall_sphere.error().code == geom::Errc::degenerate);
}

TEST_CASE("cone and sphere carry the tolerance") {
  const double tol = 1e-3;
  const auto check = [](const geom::Brep& b, double expected) {
    for (const geom::VertexId id : b.vertex_ids()) CHECK(b.vertex(id).tol == expected);
    for (const geom::EdgeId id : b.edge_ids()) CHECK(b.edge(id).tol == expected);
    for (const geom::FaceId id : b.face_ids()) CHECK(b.face(id).tol == expected);
    CHECK(b.solid(b.solid_ids().front()).tol == expected);
    CHECK(geom::validate_tolerances(b).ok());
    CHECK(b.geometry().pcurves.empty());
    for (const geom::CoedgeId id : b.coedge_ids()) {
      CHECK(b.coedge(id).pcurve == geom::kNoGeometry);
    }
  };

  const auto cone = geom::make_cone({1.0, 2.0, 3.0}, {0.0, 0.0, 1.0}, 2.0, 3.0, tol);
  REQUIRE(cone.has_value());
  check(cone.value(), tol);

  const auto sphere = geom::make_sphere({1.0, 2.0, 3.0}, 3.0, tol);
  REQUIRE(sphere.has_value());
  check(sphere.value(), tol);
}
