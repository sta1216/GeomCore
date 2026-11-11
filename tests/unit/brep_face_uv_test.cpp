#include <doctest/doctest.h>

#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <numbers>
#include <optional>
#include <string_view>
#include <variant>
#include <vector>

#include <geom/brep/face_uv.hpp>
#include <geom/brep/plane_frame.hpp>
#include <geom/brep/primitives.hpp>
#include <geom/brep/solid_editor.hpp>
#include <geom/brep/topology.hpp>
#include <geom/core/error.hpp>
#include <geom/curves/analytic_curve.hpp>
#include <geom/primitives/point.hpp>
#include <geom/surfaces/analytic.hpp>
#include <geom/surfaces/nurbs_surface.hpp>
#include <geom/surfaces/plane.hpp>

#include "brep_curved_fixture.hpp"

namespace {

constexpr double kPi = std::numbers::pi_v<double>;
constexpr double kTwoPi = 2.0 * kPi;
constexpr double kTolerance = 1e-9;

template <class Surface>
geom::FaceId face_with_surface(const geom::Brep& brep) {
  for (const geom::FaceId id : brep.face_ids()) {
    const geom::Face& face = brep.face(id);
    if (face.surface == geom::kNoGeometry) continue;
    if (std::holds_alternative<Surface>(brep.geometry().surfaces[face.surface])) return id;
  }
  return geom::FaceId{};
}

// The cap whose plane faces `axis` from the given side.
geom::FaceId cap_facing(const geom::Brep& brep, const geom::Vec3d& axis, bool bottom) {
  for (const geom::FaceId id : brep.face_ids()) {
    const geom::Face& face = brep.face(id);
    if (face.surface == geom::kNoGeometry) continue;
    const auto* plane = std::get_if<geom::Plane3>(&brep.geometry().surfaces[face.surface]);
    if (plane == nullptr) continue;
    const double side = geom::dot(plane->normal, axis);
    if (bottom ? side < 0.0 : side > 0.0) return id;
  }
  return geom::FaceId{};
}

// Unit square in the z = 0 plane with a stored pcurve per edge. The loop runs
// A -> B -> C -> D when `reversed_loop` is false and A -> D -> C -> B when it is
// true, so the reversed variant exercises Coedge::reversed on every coedge.
struct PcurveQuad {
  geom::Brep brep;
  geom::FaceId face{};
  geom::PlaneFrame frame{};
};

PcurveQuad make_pcurve_quad(bool reversed_loop, bool non_finite_pcurve = false) {
  geom::SolidEditor editor;
  const geom::Plane3 plane{geom::Vec3d{0.0, 0.0, 1.0}, 0.0};
  const geom::PlaneFrame frame = geom::make_plane_frame(plane).value();

  const geom::Point3d corners[4] = {{0.0, 0.0, 0.0}, {1.0, 0.0, 0.0}, {1.0, 1.0, 0.0},
                                    {0.0, 1.0, 0.0}};
  geom::VertexId vertices[4];
  for (int i = 0; i < 4; ++i) vertices[i] = editor.add_vertex(corners[i], 0.0).value();

  const std::uint32_t surface = editor.add_surface(plane);

  // One line curve and one edge per side, each stored from its lower-index
  // corner to the next one around the square.
  geom::EdgeId edges[4];
  for (int i = 0; i < 4; ++i) {
    const geom::Point3d& from = corners[i];
    const geom::Point3d& to = corners[(i + 1) % 4];
    const geom::Vec3d direction{to[0] - from[0], to[1] - from[1], to[2] - from[2]};
    const std::uint32_t curve = editor.add_curve(geom::LineCurve3{from, direction});
    edges[i] = editor.add_edge(vertices[i], vertices[(i + 1) % 4], curve, 0.0, 1.0, 0.0).value();
  }

  // Pcurves follow the edges from v0 to v1, never the coedge direction.
  std::uint32_t pcurves[4];
  for (int i = 0; i < 4; ++i) {
    const geom::Point2d from = geom::plane_frame_to_2d(frame, corners[i]);
    const geom::Point2d to = geom::plane_frame_to_2d(frame, corners[(i + 1) % 4]);
    if (non_finite_pcurve && i == 1) {
      const double nan = std::numeric_limits<double>::quiet_NaN();
      pcurves[i] = editor.add_pcurve(geom::Pcurve2d{from, geom::Point2d{nan, to[1]}});
    } else {
      pcurves[i] = editor.add_pcurve(geom::Pcurve2d{from, to});
    }
  }

  const geom::FaceId face = editor.add_face(surface, false, 0.0).value();
  std::array<geom::CoedgeSpec, 4> coedges{};
  if (reversed_loop) {
    // A -> D -> C -> B, i.e. the same square traversed backwards.
    coedges[0] = geom::CoedgeSpec{edges[3], true, pcurves[3], 0.0};
    coedges[1] = geom::CoedgeSpec{edges[2], true, pcurves[2], 0.0};
    coedges[2] = geom::CoedgeSpec{edges[1], true, pcurves[1], 0.0};
    coedges[3] = geom::CoedgeSpec{edges[0], true, pcurves[0], 0.0};
  } else {
    coedges[0] = geom::CoedgeSpec{edges[0], false, pcurves[0], 0.0};
    coedges[1] = geom::CoedgeSpec{edges[1], false, pcurves[1], 0.0};
    coedges[2] = geom::CoedgeSpec{edges[2], false, pcurves[2], 0.0};
    coedges[3] = geom::CoedgeSpec{edges[3], false, pcurves[3], 0.0};
  }
  editor.add_loop(face, coedges, true).value();
  return PcurveQuad{editor.brep(), face, frame};
}

}  // namespace

TEST_CASE("box faces expose four uv corners") {
  const auto box = geom::make_box(geom::Point3d{-1.0, -1.0, -1.0}, geom::Point3d{1.0, 1.0, 1.0},
                                  kTolerance);
  REQUIRE(box.has_value());
  const geom::Brep& brep = box.value();
  REQUIRE(brep.face_count() == 6u);
  const geom::FaceId first_face = brep.face_ids().front();
  const auto loops = geom::face_loops_uv(brep, first_face, kTolerance);
  REQUIRE(loops.has_value());
  REQUIRE(loops.value().size() == 1u);
  CHECK(loops.value()[0].outer);
  CHECK(loops.value()[0].points.size() == 4u);

  for (const geom::FaceId id : brep.face_ids()) {
    const auto face_loops = geom::face_loops_uv(brep, id, kTolerance);
    REQUIRE(face_loops.has_value());
    REQUIRE(face_loops.value().size() == 1u);
    CHECK(face_loops.value()[0].outer);
    CHECK(face_loops.value()[0].points.size() == 4u);
  }
}

TEST_CASE("box uv loops map back to the box corners") {
  const auto box = geom::make_box(geom::Point3d{-1.0, -1.0, -1.0}, geom::Point3d{1.0, 1.0, 1.0},
                                  kTolerance);
  REQUIRE(box.has_value());
  const geom::Brep& brep = box.value();
  for (const geom::FaceId id : brep.face_ids()) {
    const geom::Face& face = brep.face(id);
    const auto& surface = brep.geometry().surfaces[face.surface];
    const auto* plane = std::get_if<geom::Plane3>(&surface);
    REQUIRE(plane != nullptr);
    const geom::PlaneFrame frame = geom::make_plane_frame(*plane).value();

    const auto loops = geom::face_loops_uv(brep, id, kTolerance);
    REQUIRE(loops.has_value());
    REQUIRE(loops.value().size() == 1u);
    REQUIRE(loops.value()[0].points.size() == 4u);

    std::array<geom::Point3d, 4> corners{};
    for (std::size_t i = 0; i < loops.value()[0].points.size(); ++i) {
      const geom::Point3d point = geom::plane_frame_to_3d(frame, loops.value()[0].points[i]);
      for (int axis = 0; axis < 3; ++axis) {
        CHECK(std::fabs(point[axis]) == doctest::Approx(1.0));
      }
      corners[i] = point;
    }
    for (std::size_t i = 0; i < corners.size(); ++i) {
      for (std::size_t j = i + 1u; j < corners.size(); ++j) {
        CHECK(geom::norm2(corners[i] - corners[j]) > 1.0);
      }
    }
  }
}

TEST_CASE("cylinder lateral loop unwraps one full turn of u") {
  const geom::Vec3d axis{0.0, 0.0, 1.0};
  const auto cylinder = geom::make_cylinder(geom::Point3d{0.0, 0.0, 0.0}, axis, 2.0, 3.0,
                                            kTolerance);
  REQUIRE(cylinder.has_value());
  const geom::Brep& brep = cylinder.value();
  const geom::FaceId lateral = face_with_surface<geom::Cylinder3>(brep);
  REQUIRE_FALSE(lateral.is_null());

  const auto loops = geom::face_loops_uv(brep, lateral, kTolerance);
  REQUIRE(loops.has_value());
  REQUIRE(loops.value().size() == 1u);
  const std::vector<geom::Point2d>& points = loops.value()[0].points;
  REQUIRE(points.size() == 66u);

  double min_u = points[0][0];
  double max_u = points[0][0];
  double min_v = points[0][1];
  double max_v = points[0][1];
  for (const geom::Point2d& point : points) {
    min_u = std::min(min_u, point[0]);
    max_u = std::max(max_u, point[0]);
    min_v = std::min(min_v, point[1]);
    max_v = std::max(max_v, point[1]);
  }
  CHECK(max_u - min_u == doctest::Approx(kTwoPi));
  CHECK(max_v - min_v == doctest::Approx(3.0));
  for (std::size_t i = 0; i + 1u < points.size(); ++i) {
    CHECK(std::fabs(points[i + 1u][0] - points[i][0]) <= kPi);
  }

  // Every parameter pair evaluates back onto the cylinder wall.
  const auto& surface =
      std::get<geom::Cylinder3>(brep.geometry().surfaces[brep.face(lateral).surface]);
  for (const geom::Point2d& point : points) {
    const geom::Point3d point3d = surface.eval(point[0], point[1]);
    const geom::Vec3d offset = point3d - surface.origin;
    const double height = geom::dot(offset, surface.axis);
    const geom::Vec3d radial = offset - surface.axis * height;
    CHECK(geom::norm2(radial) == doctest::Approx(2.0).epsilon(1e-12));
    CHECK(height >= -kTolerance);
    CHECK(height <= 3.0 + kTolerance);
  }
}

TEST_CASE("a circular cap edge yields one closed uv ring") {
  const geom::Vec3d axis{0.0, 0.0, 1.0};
  const auto cylinder = geom::make_cylinder(geom::Point3d{0.0, 0.0, 0.0}, axis, 2.0, 3.0,
                                            kTolerance);
  REQUIRE(cylinder.has_value());
  const geom::Brep& brep = cylinder.value();
  const geom::FaceId bottom = cap_facing(brep, axis, true);
  REQUIRE_FALSE(bottom.is_null());

  const geom::Face& face = brep.face(bottom);
  const auto& plane = std::get<geom::Plane3>(brep.geometry().surfaces[face.surface]);
  const geom::PlaneFrame frame = geom::make_plane_frame(plane).value();

  const auto loops = geom::face_loops_uv(brep, bottom, kTolerance);
  REQUIRE(loops.has_value());
  REQUIRE(loops.value().size() == 1u);
  CHECK(loops.value()[0].outer);
  const std::vector<geom::Point2d>& points = loops.value()[0].points;
  REQUIRE(points.size() == 32u);
  const bool same_u = points.front()[0] == points.back()[0];
  const bool same_v = points.front()[1] == points.back()[1];
  const bool repeated_closing_point = same_u && same_v;
  CHECK_FALSE(repeated_closing_point);

  for (const geom::Point2d& point : points) {
    const geom::Point3d point3d = geom::plane_frame_to_3d(frame, point);
    const geom::Vec3d offset = point3d - geom::Point3d{0.0, 0.0, 0.0};
    const double height = geom::dot(offset, axis);
    const geom::Vec3d radial = offset - axis * height;
    CHECK(std::fabs(height) <= kTolerance);
    CHECK(geom::norm2(radial) == doctest::Approx(2.0).epsilon(1e-12));
  }
}

TEST_CASE("pcurve direction follows the coedge") {
  const geom::Point3d corners[4] = {{0.0, 0.0, 0.0}, {1.0, 0.0, 0.0}, {1.0, 1.0, 0.0},
                                    {0.0, 1.0, 0.0}};
  const PcurveQuad forward = make_pcurve_quad(false);
  const auto forward_loops = geom::face_loops_uv(forward.brep, forward.face, kTolerance);
  REQUIRE(forward_loops.has_value());
  REQUIRE(forward_loops.value().size() == 1u);
  REQUIRE(forward_loops.value()[0].points.size() == 4u);
  for (std::size_t i = 0; i < 4u; ++i) {
    const geom::Point2d expected = geom::plane_frame_to_2d(forward.frame, corners[i]);
    CHECK(forward_loops.value()[0].points[i][0] == doctest::Approx(expected[0]));
    CHECK(forward_loops.value()[0].points[i][1] == doctest::Approx(expected[1]));
  }

  const PcurveQuad reversed = make_pcurve_quad(true);
  const auto reversed_loops = geom::face_loops_uv(reversed.brep, reversed.face, kTolerance);
  REQUIRE(reversed_loops.has_value());
  REQUIRE(reversed_loops.value().size() == 1u);
  REQUIRE(reversed_loops.value()[0].points.size() == 4u);
  const std::size_t order[4] = {0u, 3u, 2u, 1u};
  for (std::size_t i = 0; i < 4u; ++i) {
    const geom::Point2d expected = geom::plane_frame_to_2d(reversed.frame, corners[order[i]]);
    CHECK(reversed_loops.value()[0].points[i][0] == doctest::Approx(expected[0]));
    CHECK(reversed_loops.value()[0].points[i][1] == doctest::Approx(expected[1]));
  }
}

TEST_CASE("face_loops_uv rejects a dead face and a bad tolerance") {
  const auto box = geom::make_box(geom::Point3d{-1.0, -1.0, -1.0}, geom::Point3d{1.0, 1.0, 1.0},
                                  kTolerance);
  REQUIRE(box.has_value());
  const geom::Brep& brep = box.value();

  const auto dead = geom::face_loops_uv(brep, geom::FaceId{}, kTolerance);
  REQUIRE_FALSE(dead.has_value());
  CHECK(dead.error().code == geom::Errc::invalid_topology);
  CHECK(std::string_view(dead.error().message).starts_with("face_loops_uv:"));

  const geom::FaceId first_face = brep.face_ids().front();
  const double nan = std::numeric_limits<double>::quiet_NaN();
  const double inf = std::numeric_limits<double>::infinity();
  for (const double bad : {-1.0, nan, inf}) {
    const auto result = geom::face_loops_uv(brep, first_face, bad);
    REQUIRE_FALSE(result.has_value());
    CHECK(result.error().code == geom::Errc::tolerance_violation);
    CHECK(std::string_view(result.error().message).starts_with("face_loops_uv:"));
  }
}

TEST_CASE("a face without a loop is rejected") {
  geom::SolidEditor editor;
  const std::uint32_t surface = editor.add_surface(geom::Plane3{geom::Vec3d{0.0, 0.0, 1.0}, 0.0});
  const geom::FaceId face = editor.add_face(surface, false, 0.0).value();

  const auto loops = geom::face_loops_uv(editor.brep(), face, kTolerance);
  REQUIRE_FALSE(loops.has_value());
  CHECK(loops.error().code == geom::Errc::invalid_topology);
  CHECK(std::string_view(loops.error().message).starts_with("face_loops_uv:"));
}

TEST_CASE("a NURBS surface face is unsupported") {
  geom::SolidEditor editor;
  const std::uint32_t surface = editor.add_surface(
      geom::NurbsSurface(1, 1, {0.0, 0.0, 1.0, 1.0}, {0.0, 0.0, 1.0, 1.0},
                         {{0.0, 0.0, 0.0}, {1.0, 0.0, 0.0}, {0.0, 1.0, 0.0}, {1.0, 1.0, 0.0}},
                         {1.0, 1.0, 1.0, 1.0}));
  const geom::VertexId a = editor.add_vertex(geom::Point3d{0.0, 0.0, 0.0}, 0.0).value();
  const geom::VertexId b = editor.add_vertex(geom::Point3d{1.0, 0.0, 0.0}, 0.0).value();
  const std::uint32_t curve = editor.add_curve(
      geom::LineCurve3{geom::Point3d{0.0, 0.0, 0.0}, geom::Vec3d{1.0, 0.0, 0.0}});
  const geom::EdgeId ab = editor.add_edge(a, b, curve, 0.0, 1.0, 0.0).value();
  const geom::FaceId face = editor.add_face(surface, false, 0.0).value();
  const std::array<geom::CoedgeSpec, 1> coedges{geom::CoedgeSpec{ab, false, geom::kNoGeometry,
                                                                 0.0}};
  editor.add_loop(face, coedges, true).value();

  const auto loops = geom::face_loops_uv(editor.brep(), face, kTolerance);
  REQUIRE_FALSE(loops.has_value());
  CHECK(loops.error().code == geom::Errc::unsupported);
}

TEST_CASE("a non-finite pcurve point is rejected") {
  const PcurveQuad quad = make_pcurve_quad(false, true);
  const auto loops = geom::face_loops_uv(quad.brep, quad.face, kTolerance);
  REQUIRE_FALSE(loops.has_value());
  CHECK(loops.error().code == geom::Errc::non_finite);
  CHECK(std::string_view(loops.error().message).starts_with("face_loops_uv:"));
}

TEST_CASE("an open loop is rejected as not closed") {
  geom::SolidEditor editor;
  const std::uint32_t surface = editor.add_surface(geom::Plane3{geom::Vec3d{0.0, 0.0, 1.0}, 0.0});
  const geom::Point3d a_point{0.0, 0.0, 0.0};
  const geom::Point3d b_point{1.0, 0.0, 0.0};
  const std::uint32_t curve =
      editor.add_curve(geom::LineCurve3{a_point, geom::Vec3d{1.0, 0.0, 0.0}});
  const geom::VertexId a = editor.add_vertex(a_point, 0.0).value();
  const geom::VertexId b = editor.add_vertex(b_point, 0.0).value();
  const geom::EdgeId ab = editor.add_edge(a, b, curve, 0.0, 1.0, 0.0).value();
  const geom::FaceId face = editor.add_face(surface, false, 0.0).value();
  const std::array<geom::CoedgeSpec, 1> coedges{geom::CoedgeSpec{ab, false, geom::kNoGeometry,
                                                                 0.0}};
  editor.add_loop(face, coedges, true).value();

  const auto loops = geom::face_loops_uv(editor.brep(), face, kTolerance);
  REQUIRE_FALSE(loops.has_value());
  CHECK(loops.error().code == geom::Errc::numeric_failure);
  CHECK(std::string_view(loops.error().message).starts_with("face_loops_uv: loop is not closed"));
}

TEST_CASE("a loop that collapses to two points is rejected") {
  geom::SolidEditor editor;
  const std::uint32_t surface = editor.add_surface(geom::Plane3{geom::Vec3d{0.0, 0.0, 1.0}, 0.0});
  const geom::Point3d a_point{0.0, 0.0, 0.0};
  const geom::Point3d b_point{1.0, 0.0, 0.0};
  const std::uint32_t curve =
      editor.add_curve(geom::LineCurve3{a_point, geom::Vec3d{1.0, 0.0, 0.0}});
  const geom::VertexId a = editor.add_vertex(a_point, 0.0).value();
  const geom::VertexId b = editor.add_vertex(b_point, 0.0).value();
  const geom::EdgeId ab = editor.add_edge(a, b, curve, 0.0, 1.0, 0.0).value();
  const geom::FaceId face = editor.add_face(surface, false, 0.0).value();
  // The edge is used once in each direction: the ring is A -> B -> A, which
  // closes but holds only two distinct points.
  const std::array<geom::CoedgeSpec, 2> coedges{
      geom::CoedgeSpec{ab, false, geom::kNoGeometry, 0.0},
      geom::CoedgeSpec{ab, true, geom::kNoGeometry, 0.0}};
  editor.add_loop(face, coedges, true).value();

  const auto loops = geom::face_loops_uv(editor.brep(), face, kTolerance);
  REQUIRE_FALSE(loops.has_value());
  CHECK(loops.error().code == geom::Errc::invalid_topology);
  CHECK(std::string_view(loops.error().message).starts_with("face_loops_uv:"));
}

TEST_CASE("a broken coedge chain is rejected") {
  geom::SolidEditor editor;
  const std::uint32_t surface = editor.add_surface(geom::Plane3{geom::Vec3d{0.0, 0.0, 1.0}, 0.0});
  const std::uint32_t curve =
      editor.add_curve(geom::LineCurve3{geom::Point3d{0.0, 0.0, 0.0}, geom::Vec3d{1.0, 0.0, 0.0}});
  const geom::VertexId a = editor.add_vertex(geom::Point3d{0.0, 0.0, 0.0}, 0.0).value();
  const geom::VertexId b = editor.add_vertex(geom::Point3d{1.0, 0.0, 0.0}, 0.0).value();
  const geom::VertexId c = editor.add_vertex(geom::Point3d{10.0, 0.0, 0.0}, 0.0).value();
  const geom::VertexId d = editor.add_vertex(geom::Point3d{11.0, 0.0, 0.0}, 0.0).value();
  const geom::EdgeId ab = editor.add_edge(a, b, curve, 0.0, 1.0, 0.0).value();
  const geom::EdgeId cd = editor.add_edge(c, d, curve, 0.0, 1.0, 0.0).value();
  const geom::FaceId face = editor.add_face(surface, false, 0.0).value();
  const std::array<geom::CoedgeSpec, 2> coedges{
      geom::CoedgeSpec{ab, false, geom::kNoGeometry, 0.0},
      geom::CoedgeSpec{cd, false, geom::kNoGeometry, 0.0}};
  editor.add_loop(face, coedges, true).value();

  const auto loops = geom::face_loops_uv(editor.brep(), face, kTolerance);
  REQUIRE_FALSE(loops.has_value());
  CHECK(loops.error().code == geom::Errc::invalid_topology);
  CHECK(std::string_view(loops.error().message).starts_with("face_loops_uv:"));
}

TEST_CASE("a sphere seam loop stays finite at the poles") {
  const auto sphere = geom::make_sphere(geom::Point3d{0.0, 0.0, 0.0}, 2.0, kTolerance);
  REQUIRE(sphere.has_value());
  const geom::Brep& brep = sphere.value();
  const geom::FaceId face = face_with_surface<geom::Sphere3>(brep);
  REQUIRE_FALSE(face.is_null());

  const auto loops = geom::face_loops_uv(brep, face, kTolerance);
  REQUIRE(loops.has_value());
  REQUIRE(loops.value().size() == 1u);
  const std::vector<geom::Point2d>& points = loops.value()[0].points;
  REQUIRE(points.size() >= 4u);

  const auto& surface =
      std::get<geom::Sphere3>(brep.geometry().surfaces[brep.face(face).surface]);
  double min_v = points[0][1];
  double max_v = points[0][1];
  for (const geom::Point2d& point : points) {
    // The projection keeps the previous u at a pole instead of returning NaN.
    CHECK(std::isfinite(point[0]));
    CHECK(std::isfinite(point[1]));
    CHECK(std::fabs(point[0]) <= kPi + kTolerance);
    CHECK(point[1] >= -0.5 * kPi - kTolerance);
    CHECK(point[1] <= 0.5 * kPi + kTolerance);
    min_v = std::min(min_v, point[1]);
    max_v = std::max(max_v, point[1]);

    const geom::Point3d point3d = surface.eval(point[0], point[1]);
    CHECK(geom::norm2(point3d - surface.center) == doctest::Approx(2.0).epsilon(1e-12));
  }
  // The seam meridian spans both poles.
  CHECK(max_v - min_v == doctest::Approx(kPi).epsilon(1e-12));
}

TEST_CASE("a torus seam loop covers both parameter periods") {
  const double major = 2.0;
  const double minor = 0.5;
  const auto torus = geom::make_torus(geom::Point3d{0.0, 0.0, 0.0}, geom::Vec3d{0.0, 0.0, 1.0},
                                      major, minor, kTolerance);
  REQUIRE(torus.has_value());
  const geom::Brep& brep = torus.value();
  const geom::FaceId face = face_with_surface<geom::Torus3>(brep);
  REQUIRE_FALSE(face.is_null());

  const auto loops = geom::face_loops_uv(brep, face, kTolerance);
  REQUIRE(loops.has_value());
  REQUIRE(loops.value().size() == 1u);
  const std::vector<geom::Point2d>& points = loops.value()[0].points;
  REQUIRE(points.size() >= 4u);

  double min_u = points[0][0];
  double max_u = points[0][0];
  double min_v = points[0][1];
  double max_v = points[0][1];
  for (std::size_t i = 0; i < points.size(); ++i) {
    min_u = std::min(min_u, points[i][0]);
    max_u = std::max(max_u, points[i][0]);
    min_v = std::min(min_v, points[i][1]);
    max_v = std::max(max_v, points[i][1]);
    if (i + 1u < points.size()) {
      CHECK(std::fabs(points[i + 1u][0] - points[i][0]) <= kPi);
      CHECK(std::fabs(points[i + 1u][1] - points[i][1]) <= kPi);
    }
  }
  // Both seam circles are traversed once, so each period appears exactly once.
  CHECK(max_u - min_u == doctest::Approx(kTwoPi));
  CHECK(max_v - min_v == doctest::Approx(kTwoPi));

  const auto& surface = std::get<geom::Torus3>(brep.geometry().surfaces[brep.face(face).surface]);
  const double outer_radius = major + minor;
  for (const geom::Point2d& point : points) {
    const geom::Point3d point3d = surface.eval(point[0], point[1]);
    const geom::Vec3d offset = point3d - surface.center;
    const double axial = geom::dot(offset, surface.axis);
    const double x = geom::dot(offset, surface.x_axis);
    const double radius = geom::norm2(offset - surface.axis * axial);
    // Each parameter pair must reconstruct a point of one of the two seam
    // circles the primitive stored: the outer equator (v = 0) or the tube
    // circle through center + major * x (u = 0).
    const double u_seam_error = std::max(std::fabs(radius - outer_radius), std::fabs(axial));
    const double v_seam_error = std::fabs(std::hypot(x - major, axial) - minor);
    CHECK(std::min(u_seam_error, v_seam_error) <= kTolerance);
  }
}
TEST_CASE("a cone lateral loop that reaches the apex is unsupported") {
  const auto cone = geom::make_cone(geom::Point3d{0.0, 0.0, 0.0}, geom::Vec3d{0.0, 0.0, 1.0}, 1.0,
                                    2.0, kTolerance);
  REQUIRE(cone.has_value());
  const geom::Brep& brep = cone.value();
  const geom::FaceId lateral = face_with_surface<geom::Cone3>(brep);
  REQUIRE_FALSE(lateral.is_null());

  // The lateral loop runs base circle -> apex -> base circle, and the apex has
  // no azimuth, so the uv preimage of the loop does not exist.
  const auto loops = geom::face_loops_uv(brep, lateral, kTolerance);
  REQUIRE_FALSE(loops.has_value());
  CHECK(loops.error().code == geom::Errc::unsupported);
  CHECK(std::string_view(loops.error().message).starts_with("face_loops_uv: cone projection"));

  // The flat cap of the same cone still exposes a planar ring.
  const geom::FaceId cap = cap_facing(brep, geom::Vec3d{0.0, 0.0, 1.0}, true);
  REQUIRE_FALSE(cap.is_null());
  const auto cap_loops = geom::face_loops_uv(brep, cap, kTolerance);
  REQUIRE(cap_loops.has_value());
  REQUIRE(cap_loops.value().size() == 1u);
  CHECK(cap_loops.value()[0].outer);
  CHECK(cap_loops.value()[0].points.size() == 32u);
}

TEST_CASE("an apex-free cone face exposes a partial u ring") {
  const double half_angle = std::atan2(1.0, 2.0);
  const double base_radius = 1.0;
  const double tan_half_angle = std::tan(half_angle);  // 0.5
  const double v_param = 1.0;
  const double rho = base_radius - v_param * tan_half_angle;  // 0.5
  REQUIRE(tan_half_angle == doctest::Approx(0.5));

  geom::SolidEditor editor;
  const std::uint32_t surface = editor.add_surface(
      geom::Cone3{geom::Point3d{0.0, 0.0, 0.0}, geom::Vec3d{0.0, 0.0, 1.0},
                  geom::Vec3d{1.0, 0.0, 0.0}, base_radius, half_angle});
  const std::uint32_t curve =
      editor.add_curve(geom::CircleCurve3{geom::Point3d{0.0, 0.0, v_param},
                                          geom::Vec3d{1.0, 0.0, 0.0},
                                          geom::Vec3d{0.0, 1.0, 0.0}, rho});
  const geom::VertexId seam = editor.add_vertex(geom::Point3d{rho, 0.0, v_param}, 0.0).value();
  const geom::EdgeId edge = editor.add_edge(seam, seam, curve, 0.0, kTwoPi, 0.0).value();
  const geom::FaceId face = editor.add_face(surface, false, 0.0).value();
  const std::array<geom::CoedgeSpec, 1> coedges{
      geom::CoedgeSpec{edge, false, geom::kNoGeometry, 0.0}};
  editor.add_loop(face, coedges, true).value();

  const geom::Brep& brep = editor.brep();
  const auto loops = geom::face_loops_uv(brep, face, kTolerance);
  REQUIRE(loops.has_value());
  REQUIRE(loops.value().size() == 1u);
  CHECK(loops.value()[0].outer);
  const std::vector<geom::Point2d>& points = loops.value()[0].points;
  // 33 samples of the closed circle, minus the repeated closing point.
  REQUIRE(points.size() == 32u);

  const auto& cone = std::get<geom::Cone3>(brep.geometry().surfaces[brep.face(face).surface]);
  const double step = kTwoPi / 32.0;
  double min_u = points[0][0];
  double max_u = points[0][0];
  for (std::size_t i = 0; i < points.size(); ++i) {
    const geom::Point2d& point = points[i];
    min_u = std::min(min_u, point[0]);
    max_u = std::max(max_u, point[0]);
    // The circle is a full turn, so every parameter pair sits at v = 1.
    CHECK(point[1] == doctest::Approx(v_param).epsilon(1e-12));
    // The returned u must follow the circle parameter t = i*step, not just be
    // self-consistent: this pins the azimuth direction and origin. A sign flip
    // (u = -atan2) or a constant offset fails here because the ring is sampled
    // at fixed t = i*step.
    double delta_u = point[0] - static_cast<double>(i) * step;
    while (delta_u > kPi) delta_u -= kTwoPi;
    while (delta_u < -kPi) delta_u += kTwoPi;
    CHECK(std::fabs(delta_u) <= 1e-12);

    const geom::Point3d reconstructed = cone.eval(point[0], point[1]);
    const geom::Point3d expected{rho * std::cos(point[0]), rho * std::sin(point[0]), v_param};
    CHECK(geom::norm2(reconstructed - expected) <= 1e-12);
    const double implicit = std::hypot(reconstructed[0], reconstructed[1]) -
                            (base_radius - tan_half_angle * reconstructed[2]);
    CHECK(std::fabs(implicit) <= 1e-12);
  }
  // The ring stops one sample short of a full turn: 2*pi - 2*pi/32.
  CHECK(max_u - min_u == doctest::Approx(kTwoPi - step).epsilon(1e-12));
}

TEST_CASE("two loops are reported outer first and keep their winding") {
  geom::SolidEditor editor;
  const std::uint32_t surface =
      editor.add_surface(geom::Plane3{geom::Vec3d{0.0, 0.0, 1.0}, 0.0});
  const geom::FaceId face = editor.add_face(surface, false, 0.0).value();

  const auto add_ring = [&](const std::array<geom::Point3d, 4>& corners, bool outer) {
    std::array<geom::VertexId, 4> vertices{};
    std::array<geom::CoedgeSpec, 4> coedges{};
    for (std::size_t i = 0; i < corners.size(); ++i) {
      vertices[i] = editor.add_vertex(corners[i], 0.0).value();
    }
    for (std::size_t i = 0; i < corners.size(); ++i) {
      const geom::Point3d& from = corners[i];
      const geom::Point3d& to = corners[(i + 1u) % corners.size()];
      const std::uint32_t curve = editor.add_curve(
          geom::LineCurve3{from, geom::Vec3d{to[0] - from[0], to[1] - from[1], to[2] - from[2]}});
      const geom::EdgeId ring_edge =
          editor.add_edge(vertices[i], vertices[(i + 1u) % corners.size()], curve, 0.0, 1.0, 0.0)
              .value();
      coedges[i] = geom::CoedgeSpec{ring_edge, false, geom::kNoGeometry, 0.0};
    }
    REQUIRE(editor.add_loop(face, coedges, outer).has_value());
  };

  // The inner ring is stored first and runs clockwise, the outer ring second and
  // counter-clockwise, so only the outer-first ordering can put them in this
  // order.
  const std::array<geom::Point3d, 4> inner_corners{
      geom::Point3d{0.75, 0.25, 0.0}, geom::Point3d{0.25, 0.25, 0.0},
      geom::Point3d{0.25, 0.75, 0.0}, geom::Point3d{0.75, 0.75, 0.0}};
  const std::array<geom::Point3d, 4> outer_corners{
      geom::Point3d{0.0, 0.0, 0.0}, geom::Point3d{1.0, 0.0, 0.0},
      geom::Point3d{1.0, 1.0, 0.0}, geom::Point3d{0.0, 1.0, 0.0}};
  add_ring(inner_corners, false);
  add_ring(outer_corners, true);

  const geom::Brep& brep = editor.brep();
  const auto loops = geom::face_loops_uv(brep, face, kTolerance);
  REQUIRE(loops.has_value());
  REQUIRE(loops.value().size() == 2u);
  CHECK(loops.value()[0].outer);
  CHECK_FALSE(loops.value()[1].outer);

  const auto signed_area = [](const std::vector<geom::Point2d>& ring) {
    double sum = 0.0;
    for (std::size_t i = 0; i < ring.size(); ++i) {
      const geom::Point2d& a = ring[i];
      const geom::Point2d& b = ring[(i + 1u) % ring.size()];
      sum += a[0] * b[1] - b[0] * a[1];
    }
    return 0.5 * sum;
  };
  const double outer_area = signed_area(loops.value()[0].points);
  const double inner_area = signed_area(loops.value()[1].points);
  CHECK(std::fabs(outer_area) == doctest::Approx(1.0));
  CHECK(std::fabs(inner_area) == doctest::Approx(0.25));
  // Opposite winding: a CCW island inside a CW hole (or the other way round).
  CHECK(outer_area * inner_area < 0.0);
}

TEST_CASE("cylinder cap and lateral loops carry circle provenance") {
  const geom::Vec3d axis{0.0, 0.0, 1.0};
  const auto cylinder =
      geom::make_cylinder(geom::Point3d{0.0, 0.0, 0.0}, axis, 2.0, 3.0, kTolerance);
  REQUIRE(cylinder.has_value());
  const geom::Brep& brep = cylinder.value();

  // Both planar caps are bounded by a rim circle, so every cap sample must carry
  // the rim circle and the azimuth of its 3D point around it.
  for (const bool bottom : {true, false}) {
    const geom::FaceId cap = cap_facing(brep, axis, bottom);
    REQUIRE_FALSE(cap.is_null());
    const geom::Face& face = brep.face(cap);
    const auto& plane = std::get<geom::Plane3>(brep.geometry().surfaces[face.surface]);
    const geom::PlaneFrame frame = geom::make_plane_frame(plane).value();
    const geom::EdgeId border = brep.coedge(brep.loop(face.first_loop).start).edge;
    REQUIRE_FALSE(border.is_null());
    const auto& edge_circle = std::get<geom::CircleCurve3>(
        brep.geometry().curves[brep.edge(border).curve]);

    const auto loops = geom::face_loops_uv(brep, cap, kTolerance);
    REQUIRE(loops.has_value());
    REQUIRE(loops.value().size() == 1u);
    const geom::FaceUvLoop& uv = loops.value()[0];
    REQUIRE(uv.points.size() == 32u);
    // The provenance is per segment: entry i describes points[i] -> points[i+1].
    REQUIRE(uv.arcs.size() == uv.points.size());

    int direction = 0;
    for (std::size_t i = 0; i < uv.points.size(); ++i) {
      REQUIRE(uv.arcs[i].has_value());
      const geom::FaceUvArcRef& ref = *uv.arcs[i];
      CHECK(geom::norm2(ref.circle.center - edge_circle.center) <= kTolerance);
      CHECK(std::fabs(ref.circle.radius - edge_circle.radius) <= kTolerance);
      CHECK(geom::norm2(ref.circle.x_axis - edge_circle.x_axis) <= kTolerance);
      CHECK(geom::norm2(ref.circle.y_axis - edge_circle.y_axis) <= kTolerance);
      CHECK(ref.t0 >= 0.0);
      CHECK(ref.t0 < kTwoPi);

      const std::size_t next = (i + 1u) % uv.points.size();
      const geom::Point3d start3d = geom::plane_frame_to_3d(frame, uv.points[i]);
      const geom::Point3d end3d = geom::plane_frame_to_3d(frame, uv.points[next]);
      CHECK(geom::norm2(ref.circle.eval(ref.t0) - start3d) <= 1e-12);
      CHECK(geom::norm2(ref.circle.eval(ref.t1) - end3d) <= 1e-12);

      const double sweep = ref.t1 - ref.t0;
      CHECK(std::fabs(std::fabs(sweep) - kTwoPi / 32.0) <= 1e-12);
      const int sign = sweep > 0.0 ? 1 : -1;
      if (direction == 0) direction = sign;
      CHECK(sign == direction);
    }
  }

  const geom::FaceId lateral = face_with_surface<geom::Cylinder3>(brep);
  REQUIRE_FALSE(lateral.is_null());
  const auto loops = geom::face_loops_uv(brep, lateral, kTolerance);
  REQUIRE(loops.has_value());
  REQUIRE(loops.value().size() == 1u);
  const geom::FaceUvLoop& uv = loops.value()[0];
  REQUIRE(uv.points.size() == 66u);
  // Entry i describes points[i] -> points[i+1]; the two seam segments are lines.
  REQUIRE(uv.arcs.size() == uv.points.size());

  const auto& wall =
      std::get<geom::Cylinder3>(brep.geometry().surfaces[brep.face(lateral).surface]);
  std::size_t engaged = 0;
  int bottom_direction = 0;
  int top_direction = 0;
  for (std::size_t i = 0; i < uv.points.size(); ++i) {
    const geom::Point3d point3d = wall.eval(uv.points[i][0], uv.points[i][1]);
    // The seam contributes only its two endpoints, and both are shared with a
    // rim, so every lateral sample lies on one of the two rim circles.
    const bool on_rim = std::fabs(point3d[2]) <= kTolerance ||
                        std::fabs(point3d[2] - 3.0) <= kTolerance;
    CHECK(on_rim);
    if (i == 32u || i == 65u) {
      // Both vertical seam segments are straight lines, so they must stay
      // disengaged even though their endpoints sit on a rim circle.
      CHECK_FALSE(uv.arcs[i].has_value());
      continue;
    }
    REQUIRE(uv.arcs[i].has_value());
    ++engaged;

    const geom::FaceUvArcRef& ref = *uv.arcs[i];
    const std::size_t next = (i + 1u) % uv.points.size();
    const geom::Point3d end3d = wall.eval(uv.points[next][0], uv.points[next][1]);
    CHECK(std::fabs(ref.circle.radius - 2.0) <= 1e-12);
    CHECK(std::fabs(ref.circle.center[2] - point3d[2]) <= 1e-12);
    CHECK(ref.t0 >= 0.0);
    CHECK(ref.t0 < kTwoPi);
    CHECK(geom::norm2(ref.circle.eval(ref.t0) - point3d) <= 1e-12);
    CHECK(geom::norm2(ref.circle.eval(ref.t1) - end3d) <= 1e-12);

    const double sweep = ref.t1 - ref.t0;
    CHECK(std::fabs(std::fabs(sweep) - kTwoPi / 32.0) <= 1e-12);
    int& direction = i < 32u ? bottom_direction : top_direction;
    const int sign = sweep > 0.0 ? 1 : -1;
    if (direction == 0) direction = sign;
    CHECK(sign == direction);
  }
  CHECK(engaged == 64u);
  CHECK(bottom_direction == 1);
  CHECK(top_direction == -1);

  // The seam joint is carried by the neighbouring segment, not duplicated on a
  // point, so the top rim starts at t0 = 0 and runs backwards from there.
  CHECK(std::fabs(uv.points[33][1] - 3.0) <= 1e-12);
  REQUIRE(uv.arcs[33].has_value());
  CHECK(std::fabs(uv.arcs[33]->t0) <= 1e-12);
  CHECK(uv.arcs[33]->t1 < 0.0);
}

TEST_CASE("box faces carry no circle provenance") {
  const auto box =
      geom::make_box(geom::Point3d{-1.0, -1.0, -1.0}, geom::Point3d{1.0, 1.0, 1.0}, kTolerance);
  REQUIRE(box.has_value());
  const geom::Brep& brep = box.value();
  for (const geom::FaceId id : brep.face_ids()) {
    const auto loops = geom::face_loops_uv(brep, id, kTolerance);
    REQUIRE(loops.has_value());
    for (const geom::FaceUvLoop& loop : loops.value()) {
      REQUIRE(loop.arcs.size() == loop.points.size());
      for (const auto& ref : loop.arcs) CHECK_FALSE(ref.has_value());
    }
  }
}

TEST_CASE("a half disk reports normalized circle angles") {
  geom::SolidEditor editor;
  const std::uint32_t surface =
      editor.add_surface(geom::Plane3{geom::Vec3d{0.0, 0.0, 1.0}, 0.0});
  const std::uint32_t curve =
      editor.add_curve(geom::CircleCurve3{geom::Point3d{0.0, 0.0, 0.0},
                                          geom::Vec3d{1.0, 0.0, 0.0},
                                          geom::Vec3d{0.0, 1.0, 0.0}, 1.0});
  const geom::VertexId low = editor.add_vertex(geom::Point3d{0.0, -1.0, 0.0}, 0.0).value();
  const geom::VertexId high = editor.add_vertex(geom::Point3d{0.0, 1.0, 0.0}, 0.0).value();
  // The edge starts at angle -pi/2, so its parameter is negative while the
  // stored provenance angle must be normalized to [0, 2*pi).
  const geom::EdgeId arc =
      editor.add_edge(low, high, curve, -kPi / 2.0, kPi / 2.0, 0.0).value();
  const std::uint32_t chord = editor.add_curve(
      geom::LineCurve3{geom::Point3d{0.0, 1.0, 0.0}, geom::Vec3d{0.0, -2.0, 0.0}});
  const geom::EdgeId line = editor.add_edge(high, low, chord, 0.0, 1.0, 0.0).value();
  const geom::FaceId face = editor.add_face(surface, false, 0.0).value();
  const std::array<geom::CoedgeSpec, 2> coedges{geom::CoedgeSpec{arc, false, geom::kNoGeometry, 0.0},
                                                geom::CoedgeSpec{line, false, geom::kNoGeometry,
                                                                 0.0}};
  REQUIRE(editor.add_loop(face, coedges, true).has_value());

  const geom::Brep& brep = editor.brep();
  const geom::Face& stored = brep.face(face);
  const auto& plane = std::get<geom::Plane3>(brep.geometry().surfaces[stored.surface]);
  const geom::PlaneFrame frame = geom::make_plane_frame(plane).value();
  const auto loops = geom::face_loops_uv(brep, face, kTolerance);
  REQUIRE(loops.has_value());
  REQUIRE(loops.value().size() == 1u);
  const geom::FaceUvLoop& uv = loops.value()[0];
  // 33 arc samples; the chord only adds its far endpoint, which the closing pop
  // removes again.
  REQUIRE(uv.points.size() == 33u);
  REQUIRE(uv.arcs.size() == uv.points.size());
  for (std::size_t i = 0; i + 1u < uv.points.size(); ++i) {
    REQUIRE(uv.arcs[i].has_value());
    const geom::FaceUvArcRef& ref = *uv.arcs[i];
    CHECK(ref.t0 >= 0.0);
    CHECK(ref.t0 < kTwoPi);
    const geom::Point3d start3d = geom::plane_frame_to_3d(frame, uv.points[i]);
    const geom::Point3d end3d = geom::plane_frame_to_3d(frame, uv.points[i + 1u]);
    CHECK(geom::norm2(ref.circle.eval(ref.t0) - start3d) <= 1e-12);
    CHECK(geom::norm2(ref.circle.eval(ref.t1) - end3d) <= 1e-12);
    CHECK(std::fabs(std::fabs(ref.t1 - ref.t0) - kPi / 32.0) <= 1e-12);
  }
  // The closing chord is a straight diameter, so it must stay disengaged.
  CHECK_FALSE(uv.arcs[32].has_value());
  // The first sample sits at angle -pi/2 and must come back as 3*pi/2.
  CHECK(uv.arcs.front()->t0 > kPi);
}

TEST_CASE("a chordal quarter cylinder carries no circle provenance") {
  const geom::Brep brep = geom::test::make_quarter_cylinder(1.0, 2.0, kTolerance);
  const geom::FaceId lateral = geom::test::quarter_cylinder_lateral(brep);
  REQUIRE_FALSE(lateral.is_null());

  const auto loops = geom::face_loops_uv(brep, lateral, kTolerance);
  REQUIRE(loops.has_value());
  REQUIRE(loops.value().size() == 1u);
  const geom::FaceUvLoop& uv = loops.value()[0];
  REQUIRE(uv.points.size() == 34u);
  REQUIRE(uv.arcs.size() == uv.points.size());
  for (const auto& ref : uv.arcs) CHECK_FALSE(ref.has_value());

  // The fixture approximates both quarter arcs with 16 line chords, so no coedge
  // can engage. The brief's "engaged refs only on the two arcs" needs a fixture
  // whose arc edges are real CircleCurve3 curves.
  const geom::Loop& boundary = brep.loop(brep.face(lateral).first_loop);
  geom::CoedgeId current = boundary.start;
  std::size_t coedges = 0;
  do {
    const geom::Coedge& coedge = brep.coedge(current);
    const geom::Edge& edge = brep.edge(coedge.edge);
    CHECK(std::holds_alternative<geom::LineCurve3>(brep.geometry().curves[edge.curve]));
    current = coedge.next;
    ++coedges;
  } while (current != boundary.start && coedges < 100u);
  CHECK(coedges == 34u);
}

TEST_CASE("an arc corner between two circles keeps the earlier provenance") {
  geom::SolidEditor editor;
  const std::uint32_t surface =
      editor.add_surface(geom::Plane3{geom::Vec3d{0.0, 0.0, 1.0}, 0.0});
  const double half = std::sqrt(0.5);
  // Two arcs meet at (0, 1, 0) but run on different circles: a corner between
  // two arcs, like the seam of a torus.
  const std::uint32_t first_curve =
      editor.add_curve(geom::CircleCurve3{geom::Point3d{0.0, 0.0, 0.0},
                                          geom::Vec3d{1.0, 0.0, 0.0},
                                          geom::Vec3d{0.0, 1.0, 0.0}, 1.0});
  const std::uint32_t second_curve =
      editor.add_curve(geom::CircleCurve3{geom::Point3d{0.5, 0.5, 0.0},
                                          geom::Vec3d{1.0, 0.0, 0.0},
                                          geom::Vec3d{0.0, 1.0, 0.0}, half});
  const geom::VertexId right = editor.add_vertex(geom::Point3d{1.0, 0.0, 0.0}, 0.0).value();
  const geom::VertexId top = editor.add_vertex(geom::Point3d{0.0, 1.0, 0.0}, 0.0).value();
  const geom::EdgeId first =
      editor.add_edge(right, top, first_curve, 0.0, kPi / 2.0, 0.0).value();
  const geom::EdgeId second = editor.add_edge(top, right, second_curve, 3.0 * kPi / 4.0,
                                              7.0 * kPi / 4.0, 0.0)
                                  .value();
  const geom::FaceId face = editor.add_face(surface, false, 0.0).value();
  const std::array<geom::CoedgeSpec, 2> coedges{
      geom::CoedgeSpec{first, false, geom::kNoGeometry, 0.0},
      geom::CoedgeSpec{second, false, geom::kNoGeometry, 0.0}};
  REQUIRE(editor.add_loop(face, coedges, true).has_value());

  const auto loops = geom::face_loops_uv(editor.brep(), face, kTolerance);
  REQUIRE(loops.has_value());
  REQUIRE(loops.value().size() == 1u);
  const geom::FaceUvLoop& uv = loops.value()[0];
  // 33 samples from each arc, minus the shared junction and the closing pop.
  REQUIRE(uv.points.size() == 64u);
  REQUIRE(uv.arcs.size() == uv.points.size());

  // The shared point survives from the first arc; the segment that leaves it
  // belongs to the second arc and must carry the second circle.
  for (const auto& ref : uv.arcs) REQUIRE(ref.has_value());
  CHECK(std::fabs(uv.arcs[31]->circle.radius - 1.0) <= 1e-12);
  CHECK(geom::norm2(uv.arcs[31]->circle.center - geom::Point3d{0.0, 0.0, 0.0}) <= 1e-12);
  CHECK(std::fabs(uv.arcs[32]->circle.radius - half) <= 1e-12);
  CHECK(geom::norm2(uv.arcs[32]->circle.center - geom::Point3d{0.5, 0.5, 0.0}) <= 1e-12);
}

TEST_CASE("a stale pcurve off its circle edge is rejected") {
  geom::SolidEditor editor;
  const std::uint32_t surface =
      editor.add_surface(geom::Plane3{geom::Vec3d{0.0, 0.0, 1.0}, 0.0});
  const std::uint32_t curve =
      editor.add_curve(geom::CircleCurve3{geom::Point3d{0.0, 0.0, 0.0},
                                          geom::Vec3d{1.0, 0.0, 0.0},
                                          geom::Vec3d{0.0, 1.0, 0.0}, 1.0});
  // The pcurve maps to 3D points at radius sqrt(2) instead of 1, like a stale
  // parameterization left behind by an earlier edit.
  const std::uint32_t stale =
      editor.add_pcurve(geom::Pcurve2d{geom::Point2d{-1.0, -1.0}, geom::Point2d{1.0, 1.0}});
  const geom::VertexId start = editor.add_vertex(geom::Point3d{1.0, 0.0, 0.0}, 0.0).value();
  const geom::EdgeId edge = editor.add_edge(start, start, curve, 0.0, kTwoPi, 0.0).value();
  const geom::FaceId face = editor.add_face(surface, false, 0.0).value();
  const std::array<geom::CoedgeSpec, 1> coedges{
      geom::CoedgeSpec{edge, false, stale, 0.0}};
  REQUIRE(editor.add_loop(face, coedges, true).has_value());

  const auto loops = geom::face_loops_uv(editor.brep(), face, kTolerance);
  REQUIRE_FALSE(loops.has_value());
  CHECK(loops.error().code == geom::Errc::numeric_failure);
  CHECK(std::string_view(loops.error().message) ==
        "face_loops_uv: point is off its circle");
}

TEST_CASE("a two point semicircle pcurve has no resolvable arc direction") {
  // The band around a half turn is 1e-6 wide: an exact half turn and an arc
  // 1e-6 short of one both have two nearly equal candidates, so the direction
  // is rejected rather than guessed.
  for (const double end_angle : {kPi, kPi - 1e-6}) {
    geom::SolidEditor editor;
    const std::uint32_t surface =
        editor.add_surface(geom::Plane3{geom::Vec3d{0.0, 0.0, 1.0}, 0.0});
    const std::uint32_t curve =
        editor.add_curve(geom::CircleCurve3{geom::Point3d{0.0, 0.0, 0.0},
                                            geom::Vec3d{1.0, 0.0, 0.0},
                                            geom::Vec3d{0.0, 1.0, 0.0}, 1.0});
    const std::uint32_t flat = editor.add_pcurve(geom::Pcurve2d{
        geom::Point2d{1.0, 0.0},
        geom::Point2d{std::cos(end_angle), std::sin(end_angle)}});
    const geom::VertexId start = editor.add_vertex(geom::Point3d{1.0, 0.0, 0.0}, 0.0).value();
    const geom::VertexId end =
        editor
            .add_vertex(geom::Point3d{std::cos(end_angle), std::sin(end_angle), 0.0}, 0.0)
            .value();
    const geom::EdgeId edge = editor.add_edge(start, end, curve, 0.0, end_angle, 0.0).value();
    const geom::FaceId face = editor.add_face(surface, false, 0.0).value();
    const std::array<geom::CoedgeSpec, 1> coedges{
        geom::CoedgeSpec{edge, false, flat, 0.0}};
    REQUIRE(editor.add_loop(face, coedges, true).has_value());

    const auto loops = geom::face_loops_uv(editor.brep(), face, kTolerance);
    REQUIRE_FALSE(loops.has_value());
    CHECK(loops.error().code == geom::Errc::numeric_failure);
    CHECK(std::string_view(loops.error().message) ==
          "face_loops_uv: circular arc direction is ambiguous");
    CHECK(loops.error().context_a == 0u);
    CHECK(loops.error().context_b == 2u);
  }
}
