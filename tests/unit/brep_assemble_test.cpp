#include <doctest/doctest.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <limits>
#include <span>
#include <string>
#include <string_view>
#include <variant>
#include <vector>

#include <geom/brep/assemble.hpp>
#include <geom/brep/mass_properties.hpp>
#include <geom/brep/validate.hpp>
#include <geom/core/error.hpp>

using namespace geom;

namespace {

AssembledFace quad(Point3d p0, Point3d p1, Point3d p2, Point3d p3,
                   Vec3d n) {
  return AssembledFace{{p0, p1, p2, p3}, {}, n};
}

std::vector<AssembledFace> box_faces(Point3d lo, Point3d hi) {
  const double x0 = lo[0], y0 = lo[1], z0 = lo[2];
  const double x1 = hi[0], y1 = hi[1], z1 = hi[2];
  return {
      quad({x0, y0, z1}, {x1, y0, z1}, {x1, y1, z1}, {x0, y1, z1}, {0, 0, 1}),
      quad({x0, y0, z0}, {x0, y1, z0}, {x1, y1, z0}, {x1, y0, z0}, {0, 0, -1}),
      quad({x0, y1, z0}, {x0, y1, z1}, {x1, y1, z1}, {x1, y1, z0}, {0, 1, 0}),
      quad({x0, y0, z0}, {x1, y0, z0}, {x1, y0, z1}, {x0, y0, z1}, {0, -1, 0}),
      quad({x1, y0, z0}, {x1, y1, z0}, {x1, y1, z1}, {x1, y0, z1}, {1, 0, 0}),
      quad({x0, y0, z0}, {x0, y0, z1}, {x0, y1, z1}, {x0, y1, z0}, {-1, 0, 0}),
  };
}

bool validate_solid(const Brep& brep) {
  return validate_topology(brep).ok() && validate_tolerances(brep).ok();
}

// Pcurves are stored in edge direction: p0 projects edge.v0 and p1 projects
// edge.v1 into the coedge's face frame, regardless of coedge traversal.
double pcurve_endpoint_error(const Brep& brep) {
  double error = 0.0;
  for (const CoedgeId coedge_id : brep.coedge_ids()) {
    const Coedge& coedge = brep.coedge(coedge_id);
    REQUIRE(coedge.pcurve != kNoGeometry);
    const Edge& edge = brep.edge(coedge.edge);
    const Face& face = brep.face(coedge.face);
    REQUIRE(face.surface != kNoGeometry);
    const auto* plane =
        std::get_if<Plane3>(&brep.geometry().surfaces[face.surface]);
    REQUIRE(plane != nullptr);
    const Result<PlaneFrame> frame = make_plane_frame(*plane);
    REQUIRE(frame.has_value());
    const Point3d edge_start = brep.vertex(edge.v0).point;
    const Point3d edge_end = brep.vertex(edge.v1).point;
    const Point2d expected_start =
        plane_frame_to_2d(frame.value(), edge_start);
    const Point2d expected_end = plane_frame_to_2d(frame.value(), edge_end);
    const Pcurve2d& pcurve = brep.geometry().pcurves[coedge.pcurve];
    error = std::max(
        {error,
         std::hypot(pcurve.points.front()[0] - expected_start[0],
                    pcurve.points.front()[1] - expected_start[1]),
         std::hypot(pcurve.points.back()[0] - expected_end[0],
                    pcurve.points.back()[1] - expected_end[1])});
  }
  return error;
}

Vec3d side_normal(Point3d a, Point3d b) {
  const Vec3d along{b[0] - a[0], b[1] - a[1], 0.0};
  const Vec3d n = cross(along, Vec3d{0.0, 0.0, 1.0});
  const double len = std::hypot(n[0], n[1], n[2]);
  return Vec3d{n[0] / len, n[1] / len, n[2] / len};
}

AssembledFace side_quad(Point3d a, Point3d b, double z0, double z1) {
  return AssembledFace{
      {{a[0], a[1], z0}, {b[0], b[1], z0}, {b[0], b[1], z1}, {a[0], a[1], z1}},
      {},
      side_normal(a, b)};
}

// A 4x4x1 box with a 1x1 through hole. The top face is the keyhole ring:
// the through hole is connected to the outside by a zero-width slit. The
// bottom face is the ordinary square-with-hole ring. The slit is a zero-depth
// notch in the top face, so there are no slit side walls and the volume is 15.
std::vector<AssembledFace> keyhole_box_faces(double z0, double z1) {
  const Point3d top[12] = {
      {0.0, 0.0, z1}, {4.0, 0.0, z1}, {4.0, 1.5, z1}, {2.0, 1.5, z1},
      {2.0, 1.0, z1}, {1.0, 1.0, z1}, {1.0, 2.0, z1}, {2.0, 2.0, z1},
      {2.0, 1.5, z1}, {4.0, 1.5, z1}, {4.0, 4.0, z1}, {0.0, 4.0, z1}};
  const Point3d bottom_outer[5] = {
      {0.0, 0.0, z0}, {0.0, 4.0, z0}, {4.0, 4.0, z0}, {4.0, 1.5, z0},
      {4.0, 0.0, z0}};
  const Point3d bottom_hole[5] = {
      {2.0, 1.5, z0}, {2.0, 2.0, z0}, {1.0, 2.0, z0}, {1.0, 1.0, z0},
      {2.0, 1.0, z0}};
  const Point3d outer_ring[5] = {
      {0.0, 0.0, z1}, {4.0, 0.0, z1}, {4.0, 1.5, z1}, {4.0, 4.0, z1},
      {0.0, 4.0, z1}};
  const Point3d hole_ring[5] = {
      {2.0, 1.5, z1}, {2.0, 1.0, z1}, {1.0, 1.0, z1}, {1.0, 2.0, z1},
      {2.0, 2.0, z1}};

  std::vector<AssembledFace> faces;
  faces.push_back(
      AssembledFace{std::vector<Point3d>(top, top + 12), {}, {0.0, 0.0, 1.0}});
  faces.push_back(AssembledFace{
      std::vector<Point3d>(bottom_outer, bottom_outer + 5),
      {std::vector<Point3d>(bottom_hole, bottom_hole + 5)},
      {0.0, 0.0, -1.0}});
  for (int i = 0; i < 5; ++i) {
    faces.push_back(side_quad(outer_ring[i], outer_ring[(i + 1) % 5], z0, z1));
    faces.push_back(side_quad(hole_ring[i], hole_ring[(i + 1) % 5], z0, z1));
  }
  return faces;
}

// A 4x4x1 box with a 2x2 through hole, volume (16 - 4) * 1 = 12.
std::vector<AssembledFace> square_hole_box_faces(double z0, double z1) {
  const Point3d top_outer[4] = {
      {0.0, 0.0, z1}, {4.0, 0.0, z1}, {4.0, 4.0, z1}, {0.0, 4.0, z1}};
  const Point3d top_hole[4] = {
      {1.0, 1.0, z1}, {1.0, 3.0, z1}, {3.0, 3.0, z1}, {3.0, 1.0, z1}};
  const Point3d bottom_outer[4] = {
      {0.0, 0.0, z0}, {0.0, 4.0, z0}, {4.0, 4.0, z0}, {4.0, 0.0, z0}};
  const Point3d bottom_hole[4] = {
      {1.0, 1.0, z0}, {3.0, 1.0, z0}, {3.0, 3.0, z0}, {1.0, 3.0, z0}};

  std::vector<AssembledFace> faces;
  faces.push_back(AssembledFace{std::vector<Point3d>(top_outer, top_outer + 4),
                                {std::vector<Point3d>(top_hole, top_hole + 4)},
                                {0.0, 0.0, 1.0}});
  faces.push_back(AssembledFace{
      std::vector<Point3d>(bottom_outer, bottom_outer + 4),
      {std::vector<Point3d>(bottom_hole, bottom_hole + 4)},
      {0.0, 0.0, -1.0}});
  for (int i = 0; i < 4; ++i) {
    const int j = (i + 1) % 4;
    faces.push_back(side_quad(top_outer[i], top_outer[j], z0, z1));
    faces.push_back(side_quad(top_hole[i], top_hole[j], z0, z1));
  }
  return faces;
}

}  // namespace

TEST_CASE("assemble: six quads rebuild a closed box") {
  const auto faces = box_faces({0, 0, 0}, {2, 3, 4});
  const auto brep = assemble_solid(faces, 1e-9);
  REQUIRE(brep.has_value());
  CHECK(brep.value().vertex_count() == 8);
  CHECK(brep.value().edge_count() == 12);
  CHECK(brep.value().face_count() == 6);
  CHECK(brep.value().shell_count() == 1);
  CHECK(brep.value().solid_count() == 1);
  CHECK(mass_properties(brep.value()).value().volume ==
        doctest::Approx(24.0).epsilon(1e-12));
  CHECK(validate_solid(brep.value()));
}

TEST_CASE("assemble: a reversed outer ring is reoriented") {
  auto faces = box_faces({0, 0, 0}, {2, 3, 4});
  std::reverse(faces[0].outer.begin(), faces[0].outer.end());
  const auto brep = assemble_solid(faces, 1e-9);
  REQUIRE(brep.has_value());
  CHECK(mass_properties(brep.value()).value().volume ==
        doctest::Approx(24.0).epsilon(1e-12));
  CHECK(validate_solid(brep.value()));
}

TEST_CASE("assemble: a reversed hole ring is reoriented") {
  auto faces = square_hole_box_faces(0.0, 1.0);
  std::reverse(faces[0].holes[0].begin(), faces[0].holes[0].end());
  const auto brep = assemble_solid(faces, 1e-9);
  REQUIRE(brep.has_value());
  CHECK(mass_properties(brep.value()).value().volume ==
        doctest::Approx(12.0).epsilon(1e-12));
  CHECK(validate_solid(brep.value()));
}

TEST_CASE("assemble: a missing face is a diagnosed not_closed failure") {
  auto faces = box_faces({0, 0, 0}, {1, 1, 1});
  faces.pop_back();
  const auto brep = assemble_solid(faces, 1e-9);
  REQUIRE_FALSE(brep.has_value());
  CHECK(brep.error().code == Errc::not_closed);
  CHECK(brep.error().context_a == 3u);
  CHECK(std::string_view(brep.error().message) ==
        "assemble_solid: edge is used once");
}

TEST_CASE("assemble: welding joins independently computed intersection points") {
  // Two quads share an edge whose endpoints differ by 1e-12: build four
  // vertices and two triangles forming a square; the result must weld to four
  // vertices and fail closed (a single quad pair is an open sheet).
  const Vec3d n{0, 0, 1};
  const Point3d a{0, 0, 0}, b{1, 0, 0};
  const Point3d c{1 + 1e-12, 1, 0}, d{0, 1, 0};
  const std::vector<AssembledFace> faces{
      AssembledFace{{a, b, c}, {}, n}, AssembledFace{{a, c, d}, {}, n}};
  const auto brep = assemble_solid(faces, 1e-9);
  REQUIRE_FALSE(brep.has_value());
  CHECK(brep.error().code == Errc::not_closed);
}

TEST_CASE("assemble: keyhole notch assembles into a validated solid") {
  const auto faces = keyhole_box_faces(0.0, 1.0);
  const auto brep = assemble_solid(faces, 1e-9);
  REQUIRE(brep.has_value());
  CHECK(brep.value().vertex_count() == 20);
  CHECK(brep.value().edge_count() == 31);
  CHECK(brep.value().face_count() == 12);
  CHECK(brep.value().shell_count() == 1);
  CHECK(brep.value().solid_count() == 1);
  CHECK(brep.value().geometry().surfaces.size() == 10u);
  CHECK(mass_properties(brep.value()).value().volume ==
        doctest::Approx(15.0).epsilon(1e-12));
  CHECK(validate_solid(brep.value()));
}

TEST_CASE("assemble: square face with a rectangular hole assembles") {
  const auto faces = square_hole_box_faces(0.0, 1.0);
  const auto brep = assemble_solid(faces, 1e-9);
  REQUIRE(brep.has_value());
  CHECK(brep.value().vertex_count() == 16);
  CHECK(brep.value().edge_count() == 24);
  CHECK(brep.value().face_count() == 10);
  CHECK(brep.value().shell_count() == 1);
  CHECK(brep.value().solid_count() == 1);
  CHECK(mass_properties(brep.value()).value().volume ==
        doctest::Approx(12.0).epsilon(1e-12));
  CHECK(validate_solid(brep.value()));
}

TEST_CASE("assemble: an inward cavity becomes its own shell") {
  auto faces = box_faces({0, 0, 0}, {4, 4, 4});
  auto inner = box_faces({1, 1, 1}, {3, 3, 3});
  for (AssembledFace& face : inner) {
    std::reverse(face.outer.begin(), face.outer.end());
    face.normal = Vec3d{-face.normal[0], -face.normal[1], -face.normal[2]};
  }
  faces.insert(faces.end(), inner.begin(), inner.end());
  const auto brep = assemble_solid(faces, 1e-9);
  REQUIRE(brep.has_value());
  CHECK(brep.value().vertex_count() == 16);
  CHECK(brep.value().edge_count() == 24);
  CHECK(brep.value().face_count() == 12);
  // The outer boundary and the cavity wall are separate closed shells.
  CHECK(brep.value().shell_count() == 2);
  CHECK(brep.value().solid_count() == 1);
  CHECK(mass_properties(brep.value()).value().volume ==
        doctest::Approx(56.0).epsilon(1e-12));
  CHECK(validate_solid(brep.value()));
}

TEST_CASE("assemble: three faces at one edge are non-manifold") {
  const Point3d a{0, 0, 0}, b{1, 0, 0};
  const std::vector<AssembledFace> faces{
      AssembledFace{{a, b, {1, 1, 0}, {0, 1, 0}}, {}, {0, 0, 1}},
      AssembledFace{{a, b, {1, 0, 1}, {0, 0, 1}}, {}, {0, -1, 0}},
      AssembledFace{{a, b, {1, -1, 0}, {0, -1, 0}}, {}, {0, 0, -1}}};
  const auto brep = assemble_solid(faces, 1e-9);
  REQUIRE_FALSE(brep.has_value());
  CHECK(brep.error().code == Errc::invalid_topology);
  CHECK(std::string_view(brep.error().message) ==
        "assemble_solid: non-manifold edge");
}

TEST_CASE("assemble: two faces in the same direction at one edge are non-manifold") {
  const Point3d a{0, 0, 0}, b{1, 0, 0};
  const std::vector<AssembledFace> faces{
      AssembledFace{{a, b, {1, 1, 0}, {0, 1, 0}}, {}, {0, 0, 1}},
      AssembledFace{{a, b, {1, 0, 1}, {0, 0, 1}}, {}, {0, -1, 0}}};
  const auto brep = assemble_solid(faces, 1e-9);
  REQUIRE_FALSE(brep.has_value());
  CHECK(brep.error().code == Errc::invalid_topology);
  CHECK(std::string_view(brep.error().message) ==
        "assemble_solid: non-manifold edge");
}

TEST_CASE("assemble: invalid tolerances are rejected") {
  const auto faces = box_faces({0, 0, 0}, {1, 1, 1});
  for (const double tol : {0.0, -1e-9,
                           std::numeric_limits<double>::quiet_NaN(),
                           std::numeric_limits<double>::infinity()}) {
    const auto brep = assemble_solid(faces, tol);
    REQUIRE_FALSE(brep.has_value());
    CHECK(brep.error().code == Errc::tolerance_violation);
  }
}

TEST_CASE("assemble: non-finite vertices are rejected") {
  auto faces = box_faces({0, 0, 0}, {1, 1, 1});
  faces[0].outer[0][0] = std::numeric_limits<double>::quiet_NaN();
  const auto brep = assemble_solid(faces, 1e-9);
  REQUIRE_FALSE(brep.has_value());
  CHECK(brep.error().code == Errc::non_finite);
}

TEST_CASE("assemble: an empty face span is rejected") {
  const std::span<const AssembledFace> faces;
  const auto brep = assemble_solid(faces, 1e-9);
  REQUIRE_FALSE(brep.has_value());
  CHECK(brep.error().code == Errc::invalid_topology);
}

TEST_CASE("assemble: a ring with fewer than three points is degenerate") {
  const Vec3d n{0, 0, 1};
  const std::vector<AssembledFace> faces{
      AssembledFace{{Point3d{0, 0, 0}, Point3d{1, 0, 0}}, {}, n}};
  const auto brep = assemble_solid(faces, 1e-9);
  REQUIRE_FALSE(brep.has_value());
  CHECK(brep.error().code == Errc::degenerate);
}

TEST_CASE("assemble: duplicate consecutive ring points weld away") {
  auto faces = box_faces({0, 0, 0}, {1, 1, 1});
  faces[0].outer.insert(faces[0].outer.begin() + 1, faces[0].outer[0]);
  const auto brep = assemble_solid(faces, 1e-9);
  REQUIRE(brep.has_value());
  CHECK(brep.value().vertex_count() == 8);
  CHECK(brep.value().edge_count() == 12);
  CHECK(brep.value().face_count() == 6);
  CHECK(brep.value().shell_count() == 1);
  CHECK(brep.value().solid_count() == 1);
  CHECK(mass_properties(brep.value()).value().volume ==
        doctest::Approx(1.0).epsilon(1e-12));
  CHECK(validate_solid(brep.value()));
}

TEST_CASE("assemble: welder merges points exactly one tolerance apart") {
  detail::assemble_welder welder(1e-9);
  const auto first = welder.weld(Point3d{0.5e-9, 0.0, 0.0});
  const auto second = welder.weld(Point3d{-0.5e-9, 0.0, 0.0});
  REQUIRE(first.has_value());
  REQUIRE(second.has_value());
  CHECK(first.value() == second.value());
  CHECK(welder.points().size() == 1);
}

TEST_CASE("assemble: a non-planar ring is rejected") {
  auto faces = box_faces({0, 0, 0}, {1, 1, 1});
  // Move the (0,0,1) corner out of the top plane in all three faces that
  // share it. The side faces remain planar; the top face does not.
  faces[0].outer[0] = Point3d{0.0, 0.0, 1.1};
  faces[3].outer[3] = Point3d{0.0, 0.0, 1.1};
  faces[5].outer[1] = Point3d{0.0, 0.0, 1.1};
  const auto brep = assemble_solid(faces, 1e-9);
  REQUIRE_FALSE(brep.has_value());
  CHECK(brep.error().code == Errc::tolerance_violation);
  CHECK(std::string_view(brep.error().message).find("planar") !=
        std::string_view::npos);
}

TEST_CASE("assemble: pcurves follow edge direction and survive split_edge") {
  const auto box = assemble_solid(box_faces({0, 0, 0}, {1, 1, 1}), 1e-9);
  REQUIRE(box.has_value());
  CHECK(pcurve_endpoint_error(box.value()) <= 1e-12);

  const auto keyhole = assemble_solid(keyhole_box_faces(0.0, 1.0), 1e-9);
  REQUIRE(keyhole.has_value());
  CHECK(pcurve_endpoint_error(keyhole.value()) <= 1e-12);

  SolidEditor editor(box.value());
  const auto middle = editor.split_edge(box.value().edge_ids().front(), 0.5);
  REQUIRE(middle.has_value());
  const auto split = editor.commit();
  REQUIRE(split.has_value());
  CHECK(pcurve_endpoint_error(split.value()) <= 1e-12);
  CHECK(validate_solid(split.value()));
}
