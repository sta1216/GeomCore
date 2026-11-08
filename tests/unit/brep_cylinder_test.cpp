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

struct CylinderParts {
  geom::FaceId lateral{};
  geom::FaceId bottom_cap{};
  geom::FaceId top_cap{};
  geom::EdgeId bottom_circle{};
  geom::EdgeId top_circle{};
  geom::EdgeId seam{};
};

// Faces of a freshly built cylinder and the lateral loop's three distinct
// edges; the lateral loop is bottom circle -> seam -> top circle -> seam.
CylinderParts cylinder_parts(const geom::Brep& b, const geom::Vec3d& axis) {
  CylinderParts parts;
  for (const geom::FaceId id : b.face_ids()) {
    const geom::Face& face = b.face(id);
    const auto& surface = b.geometry().surfaces[face.surface];
    if (std::holds_alternative<geom::Cylinder3>(surface)) {
      parts.lateral = id;
      continue;
    }
    const auto& plane = std::get<geom::Plane3>(surface);
    if (geom::dot(plane.normal, axis) < 0.0) {
      parts.bottom_cap = id;
    } else {
      parts.top_cap = id;
    }
  }
  if (parts.lateral.is_null()) return parts;
  const std::vector<geom::Coedge> uses = loop_uses(b, b.face(parts.lateral).first_loop);
  if (uses.size() == 4) {
    parts.bottom_circle = uses[0].edge;
    parts.seam = uses[1].edge;
    parts.top_circle = uses[2].edge;
  }
  return parts;
}

// Samples the loop's traversal in order. Closed curves are sampled at several
// parameters per coedge so the Newell sum sees the actual polygon.
std::vector<geom::Point3d> loop_polyline(const geom::Brep& b, geom::LoopId loop,
                                         int samples_per_edge = 16) {
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

TEST_CASE("cylinder topology uses a seam edge twice") {
  const auto r = geom::make_cylinder({0, 0, 0}, {0, 0, 1}, 2.0, 5.0, 0.0);
  REQUIRE(r.has_value());
  const geom::Brep& b = r.value();
  CHECK(b.vertex_count() == 2);   // 底/顶接缝点
  CHECK(b.edge_count() == 3);     // 底圆、顶圆、接缝
  CHECK(b.face_count() == 3);     // 底盖、顶盖、侧面
  CHECK(b.coedge_count() == 6);
  CHECK(geom::validate_topology(b).ok());
  CHECK(geom::euler_characteristic(b, b.solid_ids().front()) == 2);
  CHECK(b.shell(b.solid_shells(b.solid_ids().front()).front()).closed);
}

TEST_CASE("cylinder faces have consistent outward normals") {
  const auto r = geom::make_cylinder({0, 0, 0}, {0, 0, 1}, 2.0, 5.0, 0.0);
  REQUIRE(r.has_value());
  const geom::Brep& b = r.value();
  for (const geom::Face& f : b.faces()) {
    if (!std::holds_alternative<geom::Cylinder3>(b.geometry().surfaces[f.surface])) continue;
    const auto sample = b.vertex(b.edge(b.coedge(b.loop(f.first_loop).start).edge).v0).point;
    const auto& cyl = std::get<geom::Cylinder3>(b.geometry().surfaces[f.surface]);
    const double u = std::atan2(sample[1], sample[0]);
    CHECK(geom::dot(cyl.normal(u, sample[2]), geom::Vec3d{sample[0], sample[1], 0}) > 0.0);
  }
}

TEST_CASE("cylinder rejects bad radius, height and axis") {
  CHECK_FALSE(geom::make_cylinder({0, 0, 0}, {0, 0, 1}, 0.0, 1.0, 0.0).has_value());
  CHECK_FALSE(geom::make_cylinder({0, 0, 0}, {0, 0, 1}, 1.0, -1.0, 0.0).has_value());
  CHECK_FALSE(geom::make_cylinder({0, 0, 0}, {0, 0, 2}, 1.0, 1.0, 0.0).has_value());
}

TEST_CASE("cylinder circles and seam match the frame") {
  const geom::Point3d base{1.0, 2.0, 3.0};
  const geom::Vec3d axis{0.0, 0.0, 1.0};
  const geom::Point3d top{base[0], base[1], base[2] + 5.0};
  const auto r = geom::make_cylinder(base, axis, 2.0, 5.0, 0.0);
  REQUIRE(r.has_value());
  const geom::Brep& b = r.value();
  REQUIRE(b.vertex_count() == 2);
  REQUIRE(b.edge_count() == 3);

  geom::VertexId bottom{};
  geom::VertexId upper{};
  for (const geom::VertexId id : b.vertex_ids()) {
    if (b.vertex(id).point[2] == base[2]) {
      bottom = id;
    } else {
      upper = id;
    }
  }
  REQUIRE_FALSE(bottom.is_null());
  REQUIRE_FALSE(upper.is_null());
  const geom::Point3d p0 = b.vertex(bottom).point;
  const geom::Point3d p1 = b.vertex(upper).point;
  // The seam points sit on the radial +x_axis of the frame.
  const geom::Vec3d radial = as_vec(p0) - as_vec(base);
  CHECK(std::hypot(radial[0], radial[1], radial[2]) == doctest::Approx(2.0).epsilon(1e-12));
  CHECK(p1[0] == p0[0]);
  CHECK(p1[1] == p0[1]);
  CHECK(p1[2] == doctest::Approx(base[2] + 5.0).epsilon(1e-12));

  geom::EdgeId bottom_circle{};
  geom::EdgeId top_circle{};
  geom::EdgeId seam{};
  for (const geom::EdgeId id : b.edge_ids()) {
    const geom::Edge& edge = b.edge(id);
    if (edge.v0 != edge.v1) {
      seam = id;
      continue;
    }
    const auto& circle = std::get<geom::CircleCurve3>(b.geometry().curves[edge.curve]);
    if (circle.center[2] == base[2]) {
      bottom_circle = id;
    } else {
      top_circle = id;
    }
  }
  REQUIRE_FALSE(bottom_circle.is_null());
  REQUIRE_FALSE(top_circle.is_null());
  REQUIRE_FALSE(seam.is_null());

  const geom::Edge& bottom_edge = b.edge(bottom_circle);
  const auto& bottom_curve = std::get<geom::CircleCurve3>(b.geometry().curves[bottom_edge.curve]);
  CHECK(bottom_edge.v0 == bottom);
  CHECK(bottom_edge.v1 == bottom);
  CHECK(bottom_edge.t0 == 0.0);
  CHECK(bottom_edge.t1 == doctest::Approx(kTwoPi).epsilon(1e-12));
  CHECK(bottom_curve.is_valid());
  for (int i = 0; i < 3; ++i) CHECK(bottom_curve.center[i] == base[i]);
  CHECK(bottom_curve.radius == 2.0);
  const geom::Vec3d x_axis = bottom_curve.x_axis;
  CHECK(geom::dot(x_axis, axis) == doctest::Approx(0.0).epsilon(1e-12));
  CHECK(geom::dot(x_axis, x_axis) == doctest::Approx(1.0).epsilon(1e-12));
  const geom::Vec3d y_axis = bottom_curve.y_axis;
  CHECK(geom::dot(y_axis, y_axis) == doctest::Approx(1.0).epsilon(1e-12));
  CHECK(geom::dot(y_axis, axis) == doctest::Approx(0.0).epsilon(1e-12));
  const geom::Vec3d frame_y = geom::cross(axis, x_axis);
  for (int i = 0; i < 3; ++i) CHECK(y_axis[i] == doctest::Approx(frame_y[i]).epsilon(1e-12));
  // eval(0) is the seam point.
  for (int i = 0; i < 3; ++i) {
    CHECK(bottom_curve.eval(0.0)[i] == doctest::Approx(p0[i]).epsilon(1e-12));
  }

  const geom::Edge& top_edge = b.edge(top_circle);
  const auto& top_curve = std::get<geom::CircleCurve3>(b.geometry().curves[top_edge.curve]);
  CHECK(top_edge.v0 == upper);
  CHECK(top_edge.v1 == upper);
  CHECK(top_edge.t1 == doctest::Approx(kTwoPi).epsilon(1e-12));
  for (int i = 0; i < 3; ++i) CHECK(top_curve.center[i] == top[i]);
  CHECK(top_curve.radius == 2.0);
  for (int i = 0; i < 3; ++i) CHECK(top_curve.x_axis[i] == bottom_curve.x_axis[i]);
  for (int i = 0; i < 3; ++i) CHECK(top_curve.y_axis[i] == bottom_curve.y_axis[i]);

  const geom::Edge& seam_edge = b.edge(seam);
  const auto& seam_line = std::get<geom::LineCurve3>(b.geometry().curves[seam_edge.curve]);
  CHECK(seam_edge.v0 == bottom);
  CHECK(seam_edge.v1 == upper);
  CHECK(seam_edge.t0 == 0.0);
  CHECK(seam_edge.t1 == 5.0);
  for (int i = 0; i < 3; ++i) CHECK(seam_line.origin[i] == p0[i]);
  for (int i = 0; i < 3; ++i) CHECK(seam_line.direction[i] == axis[i]);
  // The line's far end is the top vertex.
  const geom::Point3d seam_top = b.vertex(upper).point;
  for (int i = 0; i < 3; ++i) {
    CHECK(seam_line.eval(seam_edge.t1)[i] == doctest::Approx(seam_top[i]).epsilon(1e-12));
  }
}

TEST_CASE("cylinder loops use every edge twice in opposite directions") {
  const auto r = geom::make_cylinder({0, 0, 0}, {0, 0, 1}, 2.0, 5.0, 0.0);
  REQUIRE(r.has_value());
  const geom::Brep& b = r.value();
  REQUIRE(b.face_count() == 3);
  REQUIRE(b.edge_count() == 3);
  REQUIRE(b.coedge_count() == 6);

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

  geom::FaceId lateral{};
  geom::FaceId bottom_cap{};
  geom::FaceId top_cap{};
  for (const geom::FaceId id : b.face_ids()) {
    const geom::Face& face = b.face(id);
    CHECK(b.loop(face.first_loop).outer);
    const auto& surface = b.geometry().surfaces[face.surface];
    if (std::holds_alternative<geom::Cylinder3>(surface)) {
      lateral = id;
    } else {
      const auto& plane = std::get<geom::Plane3>(surface);
      if (plane.normal[2] < 0.0) {
        bottom_cap = id;
      } else {
        top_cap = id;
      }
    }
  }
  REQUIRE_FALSE(lateral.is_null());
  REQUIRE_FALSE(bottom_cap.is_null());
  REQUIRE_FALSE(top_cap.is_null());

  const std::vector<geom::Coedge> lateral_uses =
      loop_uses(b, b.face(lateral).first_loop);
  REQUIRE(lateral_uses.size() == 4);
  CHECK(lateral_uses[0].reversed == false);
  CHECK(lateral_uses[1].reversed == false);
  CHECK(lateral_uses[2].reversed == true);
  CHECK(lateral_uses[3].reversed == true);
  CHECK(lateral_uses[1].edge == lateral_uses[3].edge);  // the seam, up then down
  CHECK(lateral_uses[0].edge != lateral_uses[2].edge);  // bottom vs top circle
  const auto& bottom_curve = std::get<geom::CircleCurve3>(
      b.geometry().curves[b.edge(lateral_uses[0].edge).curve]);
  const auto& top_curve = std::get<geom::CircleCurve3>(
      b.geometry().curves[b.edge(lateral_uses[2].edge).curve]);
  CHECK(bottom_curve.center[2] == 0.0);
  CHECK(top_curve.center[2] == 5.0);

  // The caps use their circle once: the bottom reversed (its outward normal is
  // -axis), the top forward.
  const std::vector<geom::Coedge> bottom_uses =
      loop_uses(b, b.face(bottom_cap).first_loop);
  REQUIRE(bottom_uses.size() == 1);
  CHECK(bottom_uses[0].reversed == true);
  CHECK(bottom_uses[0].edge == lateral_uses[0].edge);
  const std::vector<geom::Coedge> top_uses = loop_uses(b, b.face(top_cap).first_loop);
  REQUIRE(top_uses.size() == 1);
  CHECK(top_uses[0].reversed == false);
  CHECK(top_uses[0].edge == lateral_uses[2].edge);
}

TEST_CASE("cylinder caps are planar and point out of the solid") {
  const geom::Point3d base{1.0, 2.0, 3.0};
  const geom::Vec3d axis{0.0, 0.0, 1.0};
  const geom::Point3d interior{base[0], base[1], base[2] + 2.5};
  const auto r = geom::make_cylinder(base, axis, 2.0, 5.0, 0.0);
  REQUIRE(r.has_value());
  const geom::Brep& b = r.value();

  int caps = 0;
  for (const geom::FaceId id : b.face_ids()) {
    const geom::Face& face = b.face(id);
    const auto& surface = b.geometry().surfaces[face.surface];
    if (!std::holds_alternative<geom::Plane3>(surface)) continue;
    ++caps;
    const auto& plane = std::get<geom::Plane3>(surface);
    const std::vector<geom::Coedge> uses = loop_uses(b, face.first_loop);
    REQUIRE(uses.size() == 1);
    const geom::Edge& circle_edge = b.edge(uses[0].edge);
    const auto& circle = std::get<geom::CircleCurve3>(b.geometry().curves[circle_edge.curve]);
    CHECK(plane.signed_distance(circle.center) == doctest::Approx(0.0).epsilon(1e-12));
    CHECK(plane.signed_distance(circle.eval(0.0)) == doctest::Approx(0.0).epsilon(1e-12));
    // The interior of the solid lies on the negative side of every cap plane
    // and the stored normal points out of the solid.
    CHECK(plane.signed_distance(interior) < 0.0);
    if (std::get<geom::Plane3>(surface).normal[2] < 0.0) {
      CHECK(plane.signed_distance(base) == doctest::Approx(0.0).epsilon(1e-12));
    } else {
      CHECK(geom::dot(plane.normal, axis) > 0.0);
    }
  }
  CHECK(caps == 2);
}

TEST_CASE("cylinder carries the tolerance") {
  const double tol = 1e-3;
  const auto r = geom::make_cylinder({0, 0, 0}, {0, 0, 1}, 2.0, 5.0, tol);
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

TEST_CASE("cylinder errors carry specific codes") {
  const double nan = std::numeric_limits<double>::quiet_NaN();
  const double inf = std::numeric_limits<double>::infinity();

  const auto zero_radius = geom::make_cylinder({0, 0, 0}, {0, 0, 1}, 0.0, 1.0, 0.0);
  REQUIRE_FALSE(zero_radius.has_value());
  CHECK(zero_radius.error().code == geom::Errc::degenerate);

  const auto negative_radius = geom::make_cylinder({0, 0, 0}, {0, 0, 1}, -1.0, 1.0, 0.0);
  REQUIRE_FALSE(negative_radius.has_value());
  CHECK(negative_radius.error().code == geom::Errc::degenerate);

  const auto zero_height = geom::make_cylinder({0, 0, 0}, {0, 0, 1}, 1.0, 0.0, 0.0);
  REQUIRE_FALSE(zero_height.has_value());
  CHECK(zero_height.error().code == geom::Errc::degenerate);

  const auto negative_height = geom::make_cylinder({0, 0, 0}, {0, 0, 1}, 1.0, -1.0, 0.0);
  REQUIRE_FALSE(negative_height.has_value());
  CHECK(negative_height.error().code == geom::Errc::degenerate);

  const auto zero_axis = geom::make_cylinder({0, 0, 0}, {0, 0, 0}, 1.0, 1.0, 0.0);
  REQUIRE_FALSE(zero_axis.has_value());
  CHECK(zero_axis.error().code == geom::Errc::degenerate);

  const auto long_axis = geom::make_cylinder({0, 0, 0}, {0, 0, 2}, 1.0, 1.0, 0.0);
  REQUIRE_FALSE(long_axis.has_value());
  CHECK(long_axis.error().code == geom::Errc::degenerate);

  // A seam whose endpoints round to the same point is degenerate even though
  // base and height are finite and height is positive.
  const auto zero_seam = geom::make_cylinder({0, 0, 1e308}, {0, 0, 1}, 1.0, 1e-300, 0.0);
  REQUIRE_FALSE(zero_seam.has_value());
  CHECK(zero_seam.error().code == geom::Errc::degenerate);

  const auto nan_base = geom::make_cylinder({nan, 0, 0}, {0, 0, 1}, 1.0, 1.0, 0.0);
  REQUIRE_FALSE(nan_base.has_value());
  CHECK(nan_base.error().code == geom::Errc::non_finite);

  const auto inf_axis = geom::make_cylinder({0, 0, 0}, {0, inf, 0}, 1.0, 1.0, 0.0);
  REQUIRE_FALSE(inf_axis.has_value());
  CHECK(inf_axis.error().code == geom::Errc::non_finite);

  const auto inf_radius = geom::make_cylinder({0, 0, 0}, {0, 0, 1}, inf, 1.0, 0.0);
  REQUIRE_FALSE(inf_radius.has_value());
  CHECK(inf_radius.error().code == geom::Errc::non_finite);

  const auto inf_height = geom::make_cylinder({0, 0, 0}, {0, 0, 1}, 1.0, inf, 0.0);
  REQUIRE_FALSE(inf_height.has_value());
  CHECK(inf_height.error().code == geom::Errc::non_finite);

  const auto negative_tol = geom::make_cylinder({0, 0, 0}, {0, 0, 1}, 1.0, 1.0, -1.0);
  REQUIRE_FALSE(negative_tol.has_value());
  CHECK(negative_tol.error().code == geom::Errc::tolerance_violation);

  CHECK_FALSE(geom::make_cylinder({0, 0, 0}, {0, 0, 1}, 1.0, 1.0, nan).has_value());
}

TEST_CASE("cylinder frame adapts to tilted axes") {
  const geom::Point3d base{1.0, 2.0, 3.0};
  const double radius = 2.0;
  const double height = 5.0;
  const double diagonal = 1.0 / std::sqrt(3.0);
  const geom::Vec3d axes[2] = {{1.0, 0.0, 0.0}, {diagonal, diagonal, diagonal}};
  for (const geom::Vec3d& axis : axes) {
    INFO("axis = (", axis[0], ", ", axis[1], ", ", axis[2], ")");
    const auto r = geom::make_cylinder(base, axis, radius, height, 0.0);
    REQUIRE(r.has_value());
    const geom::Brep& b = r.value();
    const CylinderParts parts = cylinder_parts(b, axis);
    REQUIRE_FALSE(parts.lateral.is_null());
    REQUIRE_FALSE(parts.bottom_cap.is_null());
    REQUIRE_FALSE(parts.top_cap.is_null());
    REQUIRE_FALSE(parts.bottom_circle.is_null());
    REQUIRE_FALSE(parts.top_circle.is_null());
    REQUIRE_FALSE(parts.seam.is_null());
    CHECK(geom::validate_topology(b).ok());

    const auto& cyl = std::get<geom::Cylinder3>(
        b.geometry().surfaces[b.face(parts.lateral).surface]);
    CHECK(cyl.radius == radius);
    CHECK(geom::dot(cyl.axis, axis) == doctest::Approx(1.0).epsilon(1e-12));
    const geom::Vec3d x = cyl.x_axis;
    CHECK(std::hypot(x[0], x[1], x[2]) == doctest::Approx(1.0).epsilon(1e-12));
    CHECK(geom::dot(x, axis) == doctest::Approx(0.0).epsilon(1e-12));
    // When the axis is x-aligned the candidate rule must pick (0,1,0), which is
    // the branch a non-x-aligned axis does not exercise.
    if (axis[0] == 1.0) {
      CHECK(geom::dot(x, geom::Vec3d{0.0, 1.0, 0.0}) == doctest::Approx(1.0).epsilon(1e-12));
    }
    const geom::Edge& bottom_edge = b.edge(parts.bottom_circle);
    const auto& bottom_circle =
        std::get<geom::CircleCurve3>(b.geometry().curves[bottom_edge.curve]);
    const geom::Vec3d y_axis = bottom_circle.y_axis;
    const geom::Vec3d frame_y = geom::cross(axis, x);
    for (int i = 0; i < 3; ++i) {
      CHECK(y_axis[i] == doctest::Approx(frame_y[i]).epsilon(1e-12));
    }
    // The bottom circle starts at base + radius * x.
    const geom::Point3d expected_p0{base[0] + radius * x[0], base[1] + radius * x[1],
                                    base[2] + radius * x[2]};
    for (int i = 0; i < 3; ++i) {
      CHECK(bottom_circle.eval(0.0)[i] == doctest::Approx(expected_p0[i]).epsilon(1e-12));
    }

    // The cap planes carry the outward unit normals -axis and +axis and contain
    // base and top respectively.
    const auto& bottom_plane = std::get<geom::Plane3>(
        b.geometry().surfaces[b.face(parts.bottom_cap).surface]);
    const auto& top_plane =
        std::get<geom::Plane3>(b.geometry().surfaces[b.face(parts.top_cap).surface]);
    const geom::Point3d top{base[0] + height * axis[0], base[1] + height * axis[1],
                            base[2] + height * axis[2]};
    CHECK(std::hypot(bottom_plane.normal[0], bottom_plane.normal[1], bottom_plane.normal[2]) ==
          doctest::Approx(1.0).epsilon(1e-12));
    CHECK(std::hypot(top_plane.normal[0], top_plane.normal[1], top_plane.normal[2]) ==
          doctest::Approx(1.0).epsilon(1e-12));
    CHECK(geom::dot(bottom_plane.normal, axis) == doctest::Approx(-1.0).epsilon(1e-12));
    CHECK(geom::dot(top_plane.normal, axis) == doctest::Approx(1.0).epsilon(1e-12));
    CHECK(bottom_plane.signed_distance(base) == doctest::Approx(0.0).epsilon(1e-12));
    CHECK(top_plane.signed_distance(top) == doctest::Approx(0.0).epsilon(1e-12));

    // The seam line reaches the top vertex exactly.
    const geom::Edge& seam_edge = b.edge(parts.seam);
    const auto& seam_line = std::get<geom::LineCurve3>(b.geometry().curves[seam_edge.curve]);
    const geom::Point3d top_point = b.vertex(seam_edge.v1).point;
    for (int i = 0; i < 3; ++i) {
      CHECK(seam_line.eval(height)[i] == doctest::Approx(top_point[i]).epsilon(1e-12));
    }

    // Both caps wind counter-clockwise about their outward normal.
    for (const geom::FaceId cap : {parts.bottom_cap, parts.top_cap}) {
      const geom::Face& face = b.face(cap);
      const auto& plane = std::get<geom::Plane3>(b.geometry().surfaces[face.surface]);
      const geom::Vec3d newell = newell_normal(loop_polyline(b, face.first_loop));
      CHECK(geom::dot(newell, plane.normal) > 0.0);
    }
  }
}
