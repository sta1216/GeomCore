#include <array>
#include <cmath>
#include <cstdint>
#include <limits>
#include <string>
#include <vector>

#include <doctest/doctest.h>
#include <geom/brep/inside.hpp>
#include <geom/brep/primitives.hpp>

using namespace geom;

namespace {

// Every Result access goes through this helper so a failed construction reports
// its Errc instead of dereferencing an empty result.
template <class T>
T take(const Result<T>& result) {
  REQUIRE(result.has_value());
  return result.value();
}

Brep make_box_fixture(Point3d lo, Point3d hi) {
  const Result<Brep> made = make_box(lo, hi, 1e-12);
  REQUIRE(made.has_value());
  return made.value();
}

// 10 x 10 x 2 plate with a centred 2 x 2 through-hole: ten faces, sixteen
// vertices, twenty-four edges, chi = 0. Copied from the local fixture in
// mass_properties_test.cpp so this file stays self-contained.
Brep make_plate_with_hole_fixture() {
  const double x0 = 0.0, x1 = 10.0, y0 = 0.0, y1 = 10.0, z0 = 0.0, z1 = 2.0;
  const double hx0 = 4.0, hx1 = 6.0, hy0 = 4.0, hy1 = 6.0;
  const Point3d points[16] = {
      {x0, y0, z0}, {x1, y0, z0}, {x1, y1, z0}, {x0, y1, z0},
      {x0, y0, z1}, {x1, y0, z1}, {x1, y1, z1}, {x0, y1, z1},
      {hx0, hy0, z0}, {hx1, hy0, z0}, {hx1, hy1, z0}, {hx0, hy1, z0},
      {hx0, hy0, z1}, {hx1, hy0, z1}, {hx1, hy1, z1}, {hx0, hy1, z1},
  };
  SolidEditor editor;
  std::array<VertexId, 16> vertices{};
  for (int i = 0; i < 16; ++i) vertices[i] = take(editor.add_vertex(points[i], 0.0));

  struct EdgeSpec {
    int a;
    int b;
  };
  const EdgeSpec edge_specs[24] = {
      {0, 1},  {1, 2},  {2, 3},  {3, 0},    // 0..3 outer bottom
      {4, 5},  {5, 6},  {6, 7},  {7, 4},    // 4..7 outer top
      {0, 4},  {1, 5},  {2, 6},  {3, 7},    // 8..11 outer verticals
      {8, 9},  {9, 10}, {10, 11}, {11, 8},  // 12..15 inner bottom
      {12, 13}, {13, 14}, {14, 15}, {15, 12},  // 16..19 inner top
      {8, 12}, {9, 13}, {10, 14}, {11, 15},    // 20..23 inner verticals
  };
  std::array<EdgeId, 24> edges{};
  for (int i = 0; i < 24; ++i) {
    const Point3d from = points[edge_specs[i].a];
    const Point3d to = points[edge_specs[i].b];
    const Vec3d delta = to - from;
    const double length = std::hypot(delta[0], delta[1], delta[2]);
    const Vec3d direction{delta[0] / length, delta[1] / length, delta[2] / length};
    const std::uint32_t curve = editor.add_curve(LineCurve3{from, direction});
    edges[i] = take(editor.add_edge(vertices[edge_specs[i].a], vertices[edge_specs[i].b],
                                    curve, 0.0, length, 0.0));
  }

  const auto make_face = [&](Vec3d normal, int corner_index,
                             const std::vector<CoedgeSpec>& outer,
                             const std::vector<CoedgeSpec>& hole) {
    const Point3d corner = points[corner_index];
    const std::uint32_t surface = editor.add_surface(
        Plane3{normal, dot(normal, Vec3d{corner[0], corner[1], corner[2]})});
    const FaceId face = take(editor.add_face(surface, false, 0.0));
    take(editor.add_loop(face, outer, true));
    if (!hole.empty()) take(editor.add_loop(face, hole, false));
    return face;
  };
  const auto quad = [&](int a, int b, int c, int d) {
    const int indices[4] = {a, b, c, d};
    std::vector<CoedgeSpec> specs;
    for (int i = 0; i < 4; ++i) {
      const int from = indices[i];
      const int to = indices[(i + 1) % 4];
      for (int e = 0; e < 24; ++e) {
        if (edge_specs[e].a == from && edge_specs[e].b == to) {
          specs.push_back(CoedgeSpec{edges[e], false});
          break;
        }
        if (edge_specs[e].b == from && edge_specs[e].a == to) {
          specs.push_back(CoedgeSpec{edges[e], true});
          break;
        }
      }
    }
    return specs;
  };

  std::array<FaceId, 10> faces{};
  faces[0] = make_face({0.0, 0.0, -1.0}, 0, quad(0, 3, 2, 1), quad(8, 9, 10, 11));
  faces[1] = make_face({0.0, 0.0, 1.0}, 4, quad(4, 5, 6, 7), quad(12, 15, 14, 13));
  faces[2] = make_face({0.0, -1.0, 0.0}, 0, quad(0, 1, 5, 4), {});
  faces[3] = make_face({0.0, 1.0, 0.0}, 2, quad(2, 3, 7, 6), {});
  faces[4] = make_face({-1.0, 0.0, 0.0}, 0, quad(3, 0, 4, 7), {});
  faces[5] = make_face({1.0, 0.0, 0.0}, 1, quad(1, 2, 6, 5), {});
  faces[6] = make_face({1.0, 0.0, 0.0}, 8, quad(8, 11, 15, 12), {});
  faces[7] = make_face({-1.0, 0.0, 0.0}, 9, quad(9, 13, 14, 10), {});
  faces[8] = make_face({0.0, 1.0, 0.0}, 8, quad(9, 8, 12, 13), {});
  faces[9] = make_face({0.0, -1.0, 0.0}, 11, quad(11, 10, 14, 15), {});

  const ShellId shell = take(editor.add_shell(faces, true));
  const std::array<ShellId, 1> shells{shell};
  take(editor.add_solid(shells, 0.0));
  const Result<Brep> made = editor.commit();
  REQUIRE(made.has_value());
  return made.value();
}

void check_inside(const Brep& brep, SolidId solid, const Point3d& p, double tol) {
  const Result<bool> inside = point_in_solid(brep, solid, p, tol);
  REQUIRE(inside.has_value());
  CHECK(inside.value());
}

void check_outside(const Brep& brep, SolidId solid, const Point3d& p, double tol) {
  const Result<bool> inside = point_in_solid(brep, solid, p, tol);
  REQUIRE(inside.has_value());
  CHECK_FALSE(inside.value());
}

void check_error(const Result<bool>& result, Errc code) {
  REQUIRE_FALSE(result.has_value());
  CHECK(result.error().code == code);
}

}  // namespace

TEST_CASE("inside: a box classifies interior and exterior points") {
  const Brep box = make_box_fixture(Point3d{0.0, 0.0, 0.0}, Point3d{2.0, 2.0, 2.0});
  const SolidId solid = box.solid_ids().front();
  check_inside(box, solid, Point3d{1.0, 1.0, 1.0}, 1e-9);
  check_inside(box, solid, Point3d{0.1, 0.1, 0.1}, 1e-9);
  check_outside(box, solid, Point3d{2.5, 1.0, 1.0}, 1e-9);
  check_outside(box, solid, Point3d{-1.0, -1.0, -1.0}, 1e-9);
}

TEST_CASE("inside: a translated box classifies points at large coordinates") {
  const Brep box =
      make_box_fixture(Point3d{1e4, 1e4, 1e4}, Point3d{1e4 + 2.0, 1e4 + 2.0, 1e4 + 2.0});
  const SolidId solid = box.solid_ids().front();
  check_inside(box, solid, Point3d{1e4 + 1.0, 1e4 + 1.0, 1e4 + 1.0}, 1e-9);
  check_outside(box, solid, Point3d{1e4 + 2.0 + 1e-6, 1e4 + 1.0, 1e4 + 1.0}, 1e-9);
}

TEST_CASE("inside: boundary points belong to the closed solid") {
  const Brep box = make_box_fixture(Point3d{0.0, 0.0, 0.0}, Point3d{2.0, 2.0, 2.0});
  const SolidId solid = box.solid_ids().front();

  // Face interior, edge and corner points are all on the solid boundary.
  check_inside(box, solid, Point3d{1.0, 1.0, 0.0}, 1e-9);
  check_inside(box, solid, Point3d{0.0, 1.0, 1.0}, 1e-9);
  check_inside(box, solid, Point3d{0.0, 0.0, 0.0}, 1e-9);

  // The closed-solid convention includes a point within tol of a face even
  // when it lies just outside the solid.
  check_inside(box, solid, Point3d{1.0, 1.0, -1e-10}, 1e-9);
  check_outside(box, solid, Point3d{1.0, 1.0, -1e-6}, 1e-9);
}

TEST_CASE("inside: a grazing probe ray is discarded and retried") {
  const Brep box = make_box_fixture(Point3d{0.0, 0.0, 0.0}, Point3d{2.0, 2.0, 2.0});
  const SolidId solid = box.solid_ids().front();

  // The first probe direction is the box diagonal (1,1,1). From (1,1,0.5) it
  // meets the vertical edge x = 2, y = 2 at z = 1.5, so that face hit is
  // ambiguous and the direction is discarded. The answer is still interior,
  // which requires the retry direction to succeed.
  check_inside(box, solid, Point3d{1.0, 1.0, 0.5}, 1e-9);
}

TEST_CASE("inside: a plate with a through-hole respects the hole rim and void") {
  const Brep plate = make_plate_with_hole_fixture();
  const SolidId solid = plate.solid_ids().front();

  check_inside(plate, solid, Point3d{2.0, 5.0, 1.0}, 1e-9);
  check_outside(plate, solid, Point3d{5.3, 5.1, 1.0}, 1e-9);

  // Points on the hole rim belong to the closed solid.
  check_inside(plate, solid, Point3d{5.0, 4.0, 0.0}, 1e-9);
  check_inside(plate, solid, Point3d{5.0, 4.0, 2.0}, 1e-9);
}

TEST_CASE("inside: an ambiguous configuration fails instead of guessing") {
  const Brep box = make_box_fixture(Point3d{0.0, 0.0, 0.0}, Point3d{1.0, 1.0, 1.0});
  const SolidId solid = box.solid_ids().front();

  // With tol = 1e-3 every tetrahedral ray from the centre leaves within tol of
  // a box corner, so all four directions are ambiguous.
  const Result<bool> result = point_in_solid(box, solid, Point3d{0.5, 0.5, 0.5}, 1e-3);
  REQUIRE_FALSE(result.has_value());
  CHECK(result.error().code == Errc::numeric_failure);
  CHECK(std::string(result.error().message) ==
        "point_in_solid: every probe ray is degenerate");
}

TEST_CASE("inside: invalid inputs and unknown solids are rejected") {
  const Brep box = make_box_fixture(Point3d{0.0, 0.0, 0.0}, Point3d{2.0, 2.0, 2.0});
  const SolidId solid = box.solid_ids().front();
  const double nan = std::numeric_limits<double>::quiet_NaN();
  const double inf = std::numeric_limits<double>::infinity();

  check_error(point_in_solid(box, SolidId::null(), Point3d{1.0, 1.0, 1.0}, 1e-9),
              Errc::invalid_topology);

  // The handle is alive in `box` but unknown in the empty B-rep.
  const Brep empty;
  check_error(point_in_solid(empty, solid, Point3d{1.0, 1.0, 1.0}, 1e-9),
              Errc::invalid_topology);

  // A solid with no shells or no faces cannot be constructed: SolidEditor
  // rejects empty member lists, so only the empty-B-rep case is reachable.

  check_error(point_in_solid(box, solid, Point3d{nan, 1.0, 1.0}, 1e-9),
              Errc::non_finite);
  check_error(point_in_solid(box, solid, Point3d{1.0, 1.0, 1.0}, 0.0),
              Errc::degenerate);
  check_error(point_in_solid(box, solid, Point3d{1.0, 1.0, 1.0}, nan),
              Errc::degenerate);
  check_error(point_in_solid(box, solid, Point3d{1.0, 1.0, 1.0}, -1e-9),
              Errc::degenerate);
  check_error(point_in_solid(box, solid, Point3d{1.0, 1.0, 1.0}, inf),
              Errc::degenerate);
}
