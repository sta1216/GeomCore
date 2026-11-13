#include <doctest/doctest.h>

#include <array>
#include <cmath>
#include <cstdint>
#include <cstring>
#include <limits>
#include <numbers>
#include <string>
#include <utility>

#include <geom/brep/boolean_solid_curved.hpp>
#include <geom/brep/detail/curved_volume.hpp>
#include <geom/brep/mass_properties.hpp>
#include <geom/brep/primitives.hpp>
#include <geom/brep/solid_editor.hpp>
#include <geom/brep/split_solid.hpp>
#include <geom/core/error.hpp>

namespace {

struct SquareFace {
  geom::SolidEditor editor;
  geom::FaceId face;
  geom::LoopId loop;
  std::array<geom::VertexId, 4> vertices;
};

SquareFace make_square_face() {
  geom::SolidEditor ed;
  SquareFace s{std::move(ed), {}, {}, {}};
  const geom::Point3d pts[4] = {{0, 0, 0}, {1, 0, 0}, {1, 1, 0}, {0, 1, 0}};
  for (int i = 0; i < 4; ++i) s.vertices[i] = s.editor.add_vertex(pts[i], 0.0).value();
  const std::uint32_t plane = s.editor.add_surface(geom::Plane3{{0, 0, 1}, 0.0});
  const std::uint32_t line = s.editor.add_curve(
      geom::LineCurve3{{0, 0, 0}, {1, 0, 0}});
  std::array<geom::EdgeId, 4> edges{};
  for (int i = 0; i < 4; ++i) {
    edges[i] = s.editor.add_edge(s.vertices[i], s.vertices[(i + 1) % 4], line, 0.0, 1.0, 0.0)
                   .value();
  }
  s.face = s.editor.add_face(plane, false, 0.0).value();
  std::array<geom::CoedgeSpec, 4> specs{};
  for (int i = 0; i < 4; ++i) specs[i] = geom::CoedgeSpec{edges[i], false};
  s.loop = s.editor.add_loop(s.face, specs, true).value();
  auto shell = s.editor.add_shell(std::array<geom::FaceId, 1>{s.face}, false).value();
  s.editor.add_solid(std::array<geom::ShellId, 1>{shell}, 0.0).value();
  return s;
}

constexpr double kPi = std::numbers::pi_v<double>;

}  // namespace

TEST_CASE("box mass properties match closed forms") {
  const auto b = geom::make_box({0, 0, 0}, {2, 3, 4}, 0.0);
  REQUIRE(b.has_value());
  const auto m = geom::mass_properties(b.value());
  REQUIRE(m.has_value());
  const geom::MassProperties& mp = m.value();
  CHECK(mp.volume == doctest::Approx(24.0).epsilon(1e-9));
  CHECK(mp.centroid[0] == doctest::Approx(1.0).epsilon(1e-9));
  CHECK(mp.centroid[1] == doctest::Approx(1.5).epsilon(1e-9));
  CHECK(mp.centroid[2] == doctest::Approx(2.0).epsilon(1e-9));
  // I_xx = m/12 (b^2 + c^2) = 24/12 (9 + 16) = 50
  CHECK(mp.inertia[0][0] == doctest::Approx(50.0).epsilon(1e-6));
}

TEST_CASE("sphere and torus volumes match closed forms") {
  const auto s = geom::make_sphere({0, 0, 0}, 3.0, 0.0);
  REQUIRE(s.has_value());
  const auto sm = geom::mass_properties(s.value());
  REQUIRE(sm.has_value());
  CHECK(sm.value().volume ==
        doctest::Approx(36.0 * 3.14159265358979323846).epsilon(1e-6));
  CHECK(sm.value().centroid[0] == doctest::Approx(0.0).epsilon(1e-9));
  CHECK(sm.value().centroid[1] == doctest::Approx(0.0).epsilon(1e-9));
  CHECK(sm.value().centroid[2] == doctest::Approx(0.0).epsilon(1e-9));

  const auto t = geom::make_torus({0, 0, 0}, {0, 0, 1}, 3.0, 1.0, 0.0);
  REQUIRE(t.has_value());
  const auto tm = geom::mass_properties(t.value());
  REQUIRE(tm.has_value());
  const double expected = 2.0 * 3.14159265358979323846 * 3.14159265358979323846 * 3.0;
  CHECK(tm.value().volume == doctest::Approx(expected).epsilon(1e-6));
  CHECK(tm.value().centroid[0] == doctest::Approx(0.0).epsilon(1e-9));
  CHECK(tm.value().centroid[1] == doctest::Approx(0.0).epsilon(1e-9));
  CHECK(tm.value().centroid[2] == doctest::Approx(0.0).epsilon(1e-9));
}

// 需要 T4 的 SquareFace / make_square_face()（复制到本文件匿名命名空间）。
TEST_CASE("open shells have no mass properties") {
  SquareFace s = make_square_face();
  const auto b = s.editor.commit();
  REQUIRE(b.has_value());
  const auto m = geom::mass_properties(b.value());
  REQUIRE_FALSE(m.has_value());
  CHECK(m.error().code == geom::Errc::not_closed);
}

TEST_CASE("box inertia is diagonal about the centroid") {
  const auto b = geom::make_box({0, 0, 0}, {2, 3, 4}, 0.0);
  REQUIRE(b.has_value());
  const auto m = geom::mass_properties(b.value());
  REQUIRE(m.has_value());
  const geom::MassProperties& mp = m.value();
  CHECK(mp.inertia[0][0] == doctest::Approx(50.0).epsilon(1e-9));
  CHECK(mp.inertia[1][1] == doctest::Approx(40.0).epsilon(1e-9));
  CHECK(mp.inertia[2][2] == doctest::Approx(26.0).epsilon(1e-9));
  CHECK(std::fabs(mp.inertia[0][1]) <= 1e-9);
  CHECK(std::fabs(mp.inertia[0][2]) <= 1e-9);
  CHECK(std::fabs(mp.inertia[1][2]) <= 1e-9);
}

TEST_CASE("cylinder mass properties match closed forms") {
  const auto b = geom::make_cylinder({0, 0, 0}, {0, 0, 1}, 1.0, 2.0, 0.0);
  REQUIRE(b.has_value());
  const auto m = geom::mass_properties(b.value());
  REQUIRE(m.has_value());
  const geom::MassProperties& mp = m.value();
  const double mass = 2.0 * kPi;
  CHECK(mp.volume == doctest::Approx(mass).epsilon(1e-6));
  CHECK(mp.centroid[0] == doctest::Approx(0.0).epsilon(1e-9));
  CHECK(mp.centroid[1] == doctest::Approx(0.0).epsilon(1e-9));
  CHECK(mp.centroid[2] == doctest::Approx(1.0).epsilon(1e-9));
  CHECK(mp.inertia[2][2] == doctest::Approx(0.5 * mass).epsilon(1e-6));
  const double transverse = mass * (3.0 + 4.0) / 12.0;
  CHECK(mp.inertia[0][0] == doctest::Approx(transverse).epsilon(1e-6));
  CHECK(mp.inertia[1][1] == doctest::Approx(transverse).epsilon(1e-6));
}

TEST_CASE("cone mass properties match closed forms") {
  const auto b = geom::make_cone({0, 0, 0}, {0, 0, 1}, 1.0, 3.0, 0.0);
  REQUIRE(b.has_value());
  const auto m = geom::mass_properties(b.value());
  REQUIRE(m.has_value());
  const geom::MassProperties& mp = m.value();
  const double mass = kPi * 3.0 / 3.0;
  CHECK(mp.volume == doctest::Approx(kPi).epsilon(1e-6));
  CHECK(mp.centroid[0] == doctest::Approx(0.0).epsilon(1e-9));
  CHECK(mp.centroid[1] == doctest::Approx(0.0).epsilon(1e-9));
  CHECK(mp.centroid[2] == doctest::Approx(0.75).epsilon(1e-6));
  CHECK(mp.inertia[2][2] == doctest::Approx(0.3 * mass).epsilon(1e-6));
  const double transverse = 3.0 * mass * (4.0 + 9.0) / 80.0;
  CHECK(mp.inertia[0][0] == doctest::Approx(transverse).epsilon(1e-6));
  CHECK(mp.inertia[1][1] == doctest::Approx(transverse).epsilon(1e-6));
}

TEST_CASE("sphere and torus inertia match closed forms") {
  const auto s = geom::make_sphere({0, 0, 0}, 3.0, 0.0);
  REQUIRE(s.has_value());
  const auto sm = geom::mass_properties(s.value());
  REQUIRE(sm.has_value());
  const double sphere_mass = 36.0 * kPi;
  const double sphere_i = 0.4 * sphere_mass * 9.0;
  CHECK(sm.value().inertia[0][0] == doctest::Approx(sphere_i).epsilon(1e-6));
  CHECK(sm.value().inertia[1][1] == doctest::Approx(sphere_i).epsilon(1e-6));
  CHECK(sm.value().inertia[2][2] == doctest::Approx(sphere_i).epsilon(1e-6));
  CHECK(std::fabs(sm.value().inertia[0][1]) <= 1e-6 * sphere_i);

  const auto t = geom::make_torus({0, 0, 0}, {0, 0, 1}, 3.0, 1.0, 0.0);
  REQUIRE(t.has_value());
  const auto tm = geom::mass_properties(t.value());
  REQUIRE(tm.has_value());
  const double torus_mass = 2.0 * kPi * kPi * 3.0;
  const double axial = torus_mass * (9.0 + 0.75);
  const double transverse = torus_mass * (4.5 + 0.625);
  CHECK(tm.value().inertia[2][2] == doctest::Approx(axial).epsilon(1e-6));
  CHECK(tm.value().inertia[0][0] == doctest::Approx(transverse).epsilon(1e-6));
  CHECK(tm.value().inertia[1][1] == doctest::Approx(transverse).epsilon(1e-6));
}

TEST_CASE("mass properties reject invalid options") {
  const auto b = geom::make_box({0, 0, 0}, {2, 3, 4}, 0.0);
  REQUIRE(b.has_value());
  for (const double deflection : {-1.0, 0.0, std::numeric_limits<double>::quiet_NaN()}) {
    geom::MassPropertiesOptions options{};
    options.sampling_deflection = deflection;
    const auto m = geom::mass_properties(b.value(), options);
    REQUIRE_FALSE(m.has_value());
    CHECK(m.error().code == geom::Errc::tolerance_violation);
  }
  for (const int order : {0, -4}) {
    geom::MassPropertiesOptions options{};
    options.quadrature_order = order;
    const auto m = geom::mass_properties(b.value(), options);
    REQUIRE_FALSE(m.has_value());
    CHECK(m.error().code == geom::Errc::tolerance_violation);
  }
}

TEST_CASE("empty breps have no mass properties") {
  const geom::Brep empty;
  const auto m = geom::mass_properties(empty);
  REQUIRE_FALSE(m.has_value());
  CHECK(m.error().code == geom::Errc::numeric_failure);
}

namespace {

// A box whose six faces are all flagged reversed, so the accumulated volume is
// negative and mass_properties must reject it instead of returning |V|.
geom::Result<geom::Brep> make_reversed_box(geom::Point3d lo, geom::Point3d hi) {
  const geom::Point3d corners[8] = {
      {lo[0], lo[1], lo[2]}, {hi[0], lo[1], lo[2]}, {hi[0], hi[1], lo[2]}, {lo[0], hi[1], lo[2]},
      {lo[0], lo[1], hi[2]}, {hi[0], lo[1], hi[2]}, {hi[0], hi[1], hi[2]}, {lo[0], hi[1], hi[2]},
  };
  geom::SolidEditor editor;
  geom::VertexId vertices[8];
  for (int i = 0; i < 8; ++i) vertices[i] = editor.add_vertex(corners[i], 0.0).value();
  const int edge_v0[12] = {0, 1, 2, 3, 4, 5, 6, 7, 1, 2, 3, 0};
  const int edge_v1[12] = {1, 2, 3, 0, 5, 6, 7, 4, 5, 6, 7, 4};
  geom::EdgeId edges[12];
  for (int i = 0; i < 12; ++i) {
    const geom::Point3d from = corners[edge_v0[i]];
    const geom::Point3d to = corners[edge_v1[i]];
    const geom::Vec3d delta = to - from;
    const double length = std::hypot(delta[0], delta[1], delta[2]);
    const geom::Vec3d direction{delta[0] / length, delta[1] / length, delta[2] / length};
    const std::uint32_t curve = editor.add_curve(geom::LineCurve3{from, direction});
    edges[i] = editor.add_edge(vertices[edge_v0[i]], vertices[edge_v1[i]], curve, 0.0, length, 0.0)
                   .value();
  }
  struct LoopUse {
    int edge;
    bool reversed;
  };
  struct FaceSpec {
    geom::Vec3d normal;
    int corner;
    LoopUse loop[4];
  };
  const FaceSpec face_specs[6] = {
      {{0.0, 0.0, -1.0}, 0, {{3, true}, {2, true}, {1, true}, {0, true}}},
      {{0.0, 0.0, 1.0}, 4, {{4, false}, {5, false}, {6, false}, {7, false}}},
      {{0.0, -1.0, 0.0}, 0, {{0, false}, {8, false}, {4, true}, {11, true}}},
      {{0.0, 1.0, 0.0}, 2, {{2, false}, {10, false}, {6, true}, {9, true}}},
      {{-1.0, 0.0, 0.0}, 0, {{11, false}, {7, true}, {10, true}, {3, false}}},
      {{1.0, 0.0, 0.0}, 1, {{1, false}, {9, false}, {5, true}, {8, true}}},
  };
  std::array<geom::FaceId, 6> faces{};
  for (int f = 0; f < 6; ++f) {
    const FaceSpec& spec = face_specs[f];
    const geom::Point3d corner = corners[spec.corner];
    const std::uint32_t surface = editor.add_surface(geom::Plane3{
        spec.normal, geom::dot(spec.normal, geom::Vec3d{corner[0], corner[1], corner[2]})});
    const geom::FaceId face = editor.add_face(surface, true, 0.0).value();
    faces[static_cast<std::size_t>(f)] = face;
    std::array<geom::CoedgeSpec, 4> coedges{};
    for (int i = 0; i < 4; ++i) {
      coedges[static_cast<std::size_t>(i)] =
          geom::CoedgeSpec{edges[spec.loop[i].edge], spec.loop[i].reversed};
    }
    editor.add_loop(face, coedges, true).value();
  }
  const geom::ShellId shell = editor.add_shell(faces, true).value();
  const std::array<geom::ShellId, 1> shells{shell};
  editor.add_solid(shells, 0.0).value();
  return editor.commit();
}

// 10 x 10 x 2 plate with a centred 2 x 2 through-hole, built directly on the
// editor: ten faces, sixteen vertices, twenty-four edges, chi = 0.
geom::Result<geom::Brep> make_plate_with_hole() {
  const double x0 = 0.0, x1 = 10.0, y0 = 0.0, y1 = 10.0, z0 = 0.0, z1 = 2.0;
  const double hx0 = 4.0, hx1 = 6.0, hy0 = 4.0, hy1 = 6.0;
  const geom::Point3d points[16] = {
      {x0, y0, z0}, {x1, y0, z0}, {x1, y1, z0}, {x0, y1, z0},
      {x0, y0, z1}, {x1, y0, z1}, {x1, y1, z1}, {x0, y1, z1},
      {hx0, hy0, z0}, {hx1, hy0, z0}, {hx1, hy1, z0}, {hx0, hy1, z0},
      {hx0, hy0, z1}, {hx1, hy0, z1}, {hx1, hy1, z1}, {hx0, hy1, z1},
  };
  geom::SolidEditor editor;
  geom::VertexId vertices[16];
  for (int i = 0; i < 16; ++i) vertices[i] = editor.add_vertex(points[i], 0.0).value();
  struct EdgeSpec {
    int a;
    int b;
  };
  const EdgeSpec edge_specs[24] = {
      {0, 1}, {1, 2}, {2, 3}, {3, 0},          // 0..3 outer bottom
      {4, 5}, {5, 6}, {6, 7}, {7, 4},          // 4..7 outer top
      {0, 4}, {1, 5}, {2, 6}, {3, 7},          // 8..11 outer verticals
      {8, 9}, {9, 10}, {10, 11}, {11, 8},      // 12..15 inner bottom
      {12, 13}, {13, 14}, {14, 15}, {15, 12},  // 16..19 inner top
      {8, 12}, {9, 13}, {10, 14}, {11, 15},    // 20..23 inner verticals
  };
  geom::EdgeId edges[24];
  for (int i = 0; i < 24; ++i) {
    const geom::Point3d from = points[edge_specs[i].a];
    const geom::Point3d to = points[edge_specs[i].b];
    const geom::Vec3d delta = to - from;
    const double length = std::hypot(delta[0], delta[1], delta[2]);
    const geom::Vec3d direction{delta[0] / length, delta[1] / length, delta[2] / length};
    const std::uint32_t curve = editor.add_curve(geom::LineCurve3{from, direction});
    edges[i] = editor.add_edge(vertices[edge_specs[i].a], vertices[edge_specs[i].b], curve, 0.0,
                               length, 0.0)
                   .value();
  }
  const auto make_face = [&](geom::Vec3d normal, int corner_index,
                             const std::vector<geom::CoedgeSpec>& outer,
                             const std::vector<geom::CoedgeSpec>& hole) {
    const geom::Point3d corner = points[corner_index];
    const std::uint32_t surface = editor.add_surface(geom::Plane3{
        normal, geom::dot(normal, geom::Vec3d{corner[0], corner[1], corner[2]})});
    const geom::FaceId face = editor.add_face(surface, false, 0.0).value();
    editor.add_loop(face, outer, true).value();
    if (!hole.empty()) editor.add_loop(face, hole, false).value();
    return face;
  };
  const auto quad = [&](int a, int b, int c, int d) {
    const int indices[4] = {a, b, c, d};
    std::vector<geom::CoedgeSpec> specs;
    for (int i = 0; i < 4; ++i) {
      const int from = indices[i];
      const int to = indices[(i + 1) % 4];
      for (int e = 0; e < 24; ++e) {
        if (edge_specs[e].a == from && edge_specs[e].b == to) {
          specs.push_back(geom::CoedgeSpec{edges[e], false});
          break;
        }
        if (edge_specs[e].b == from && edge_specs[e].a == to) {
          specs.push_back(geom::CoedgeSpec{edges[e], true});
          break;
        }
      }
    }
    return specs;
  };
  std::array<geom::FaceId, 10> faces{};
  faces[0] = make_face({0.0, 0.0, -1.0}, 0, quad(0, 3, 2, 1), quad(8, 9, 10, 11));   // bottom
  faces[1] = make_face({0.0, 0.0, 1.0}, 4, quad(4, 5, 6, 7), quad(12, 15, 14, 13));  // top
  faces[2] = make_face({0.0, -1.0, 0.0}, 0, quad(0, 1, 5, 4), {});                   // front
  faces[3] = make_face({0.0, 1.0, 0.0}, 2, quad(2, 3, 7, 6), {});                    // back
  faces[4] = make_face({-1.0, 0.0, 0.0}, 0, quad(3, 0, 4, 7), {});                   // left
  faces[5] = make_face({1.0, 0.0, 0.0}, 1, quad(1, 2, 6, 5), {});                    // right
  faces[6] = make_face({1.0, 0.0, 0.0}, 8, quad(8, 11, 15, 12), {});   // hole wall +x
  faces[7] = make_face({-1.0, 0.0, 0.0}, 9, quad(9, 13, 14, 10), {});  // hole wall -x
  faces[8] = make_face({0.0, 1.0, 0.0}, 8, quad(9, 8, 12, 13), {});    // hole wall +y
  faces[9] = make_face({0.0, -1.0, 0.0}, 11, quad(11, 10, 14, 15), {});  // hole wall -y
  const geom::ShellId shell = editor.add_shell(faces, true).value();
  const std::array<geom::ShellId, 1> shells{shell};
  editor.add_solid(shells, 0.0).value();
  return editor.commit();
}

}  // namespace

TEST_CASE("mass properties are translation invariant") {
  const double cylinder_volume = 2.0 * kPi;
  const double cylinder_transverse = cylinder_volume * 7.0 / 12.0;
  for (const geom::Point3d base : {geom::Point3d{1e3, 1e3, 0.0}, geom::Point3d{1e4, 1e4, 0.0}}) {
    const auto b = geom::make_cylinder(base, {0, 0, 1}, 1.0, 2.0, 0.0);
    REQUIRE(b.has_value());
    const auto m = geom::mass_properties(b.value());
    REQUIRE(m.has_value());
    CHECK(std::fabs(m.value().volume - cylinder_volume) / cylinder_volume <= 1e-6);
    CHECK(std::fabs(m.value().inertia[0][0] - cylinder_transverse) / cylinder_transverse <= 1e-6);
    CHECK(std::fabs(m.value().inertia[1][1] - cylinder_transverse) / cylinder_transverse <= 1e-6);
    CHECK(m.value().centroid[0] == doctest::Approx(base[0]).epsilon(1e-9));
    CHECK(m.value().centroid[1] == doctest::Approx(base[1]).epsilon(1e-9));
    CHECK(m.value().centroid[2] == doctest::Approx(1.0).epsilon(1e-9));
  }

  const auto s = geom::make_sphere({1e4, 1e4, 1e4}, 3.0, 0.0);
  REQUIRE(s.has_value());
  const auto sm = geom::mass_properties(s.value());
  REQUIRE(sm.has_value());
  const double sphere_volume = 36.0 * kPi;
  CHECK(std::fabs(sm.value().volume - sphere_volume) / sphere_volume <= 1e-6);
  const double sphere_inertia = 0.4 * sphere_volume * 9.0;
  CHECK(std::fabs(sm.value().inertia[0][0] - sphere_inertia) / sphere_inertia <= 1e-6);
  CHECK(sm.value().centroid[2] == doctest::Approx(1e4).epsilon(1e-9));

  const auto box = geom::make_box({1e6, 1e6, 1e6}, {1e6 + 2.0, 1e6 + 3.0, 1e6 + 4.0}, 0.0);
  REQUIRE(box.has_value());
  const auto bm = geom::mass_properties(box.value());
  REQUIRE(bm.has_value());
  CHECK(std::fabs(bm.value().volume - 24.0) / 24.0 <= 1e-6);
  const double expected[3] = {50.0, 40.0, 26.0};
  for (int axis = 0; axis < 3; ++axis) {
    CHECK(std::fabs(bm.value().inertia[axis][axis] - expected[axis]) / expected[axis] <= 1e-6);
  }
}

TEST_CASE("plate with a through hole matches closed forms") {
  const auto b = make_plate_with_hole();
  REQUIRE(b.has_value());
  CHECK(geom::validate_topology(b.value()).ok());
  const auto m = geom::mass_properties(b.value());
  REQUIRE(m.has_value());
  CHECK(m.value().volume == doctest::Approx(192.0).epsilon(1e-9));
  CHECK(m.value().centroid[0] == doctest::Approx(5.0).epsilon(1e-9));
  CHECK(m.value().centroid[1] == doctest::Approx(5.0).epsilon(1e-9));
  CHECK(m.value().centroid[2] == doctest::Approx(1.0).epsilon(1e-9));
  const double inertia[3] = {1728.0, 1728.0, 3328.0};
  for (int axis = 0; axis < 3; ++axis) {
    CHECK(m.value().inertia[axis][axis] == doctest::Approx(inertia[axis]).epsilon(1e-6));
  }
  CHECK(std::fabs(m.value().inertia[0][1]) <= 1e-9);
  CHECK(std::fabs(m.value().inertia[0][2]) <= 1e-9);
  CHECK(std::fabs(m.value().inertia[1][2]) <= 1e-9);
}

TEST_CASE("all-reversed box has no mass properties") {
  const auto b = make_reversed_box({0, 0, 0}, {2, 3, 4});
  REQUIRE(b.has_value());
  CHECK(geom::validate_topology(b.value()).ok());
  const auto m = geom::mass_properties(b.value());
  REQUIRE_FALSE(m.has_value());
  CHECK(m.error().code == geom::Errc::numeric_failure);
}

namespace {

geom::Brep windowed_cylinder_brep() {
  const geom::Brep cylinder =
      geom::make_cylinder(geom::Point3d{0.0, 0.0, 0.0}, geom::Vec3d{0.0, 0.0, 1.0}, 1.0, 2.0,
                          1e-9)
          .value();
  const geom::Brep box =
      geom::make_box(geom::Point3d{0.5, -0.5, 0.5}, geom::Point3d{2.0, 0.5, 1.5}, 1e-9).value();
  return geom::boolean_solid_curved(cylinder, cylinder.solid_ids().front(), box,
                                    box.solid_ids().front(), geom::SolidBooleanOp::difference,
                                    geom::CurvedBooleanOptions{1e-9})
      .value()
      .brep;
}

geom::Brep hole_window_brep() {
  const geom::Brep cylinder =
      geom::make_cylinder(geom::Point3d{0.0, 0.0, 0.0}, geom::Vec3d{0.0, 0.0, 1.0}, 1.0, 2.0,
                          1e-9)
          .value();
  const geom::Brep box =
      geom::make_box(geom::Point3d{-0.5, 0.5, 0.5}, geom::Point3d{0.5, 2.0, 1.5}, 1e-9).value();
  return geom::boolean_solid_curved(cylinder, cylinder.solid_ids().front(), box,
                                    box.solid_ids().front(), geom::SolidBooleanOp::difference,
                                    geom::CurvedBooleanOptions{1e-9})
      .value()
      .brep;
}

geom::Brep slanted_trim_brep() {
  const geom::Brep cylinder =
      geom::make_cylinder(geom::Point3d{0.0, 0.0, 0.0}, geom::Vec3d{0.0, 0.0, 1.0}, 1.0, 2.0,
                          1e-9)
          .value();
  const geom::Vec3d normal{0.3, 0.0, 1.0};
  const double length = std::sqrt(geom::dot(normal, normal));
  const geom::Plane3 plane{geom::Vec3d{normal[0] / length, normal[1] / length, normal[2] / length},
                           1.0 / length};
  return geom::split_solid_by_plane(cylinder, cylinder.solid_ids().front(), plane, 1e-9)
      .value()
      .above;
}

geom::Brep scaled_window_brep() {
  // r=10, h=5 with a rectilinear window that crosses the seam (x >= 5, |y| <= 5)
  // between z=1 and z=4.
  const geom::Brep cylinder =
      geom::make_cylinder(geom::Point3d{0.0, 0.0, 0.0}, geom::Vec3d{0.0, 0.0, 1.0}, 10.0, 5.0,
                          1e-9)
          .value();
  const geom::Brep box =
      geom::make_box(geom::Point3d{5.0, -5.0, 1.0}, geom::Point3d{20.0, 5.0, 4.0}, 1e-9).value();
  return geom::boolean_solid_curved(cylinder, cylinder.solid_ids().front(), box,
                                    box.solid_ids().front(), geom::SolidBooleanOp::difference,
                                    geom::CurvedBooleanOptions{1e-9})
      .value()
      .brep;
}

std::array<std::uint64_t, 10> mass_bits(const geom::MassProperties& mass) {
  const double values[10] = {mass.volume,
                             mass.centroid[0],
                             mass.centroid[1],
                             mass.centroid[2],
                             mass.inertia[0][0],
                             mass.inertia[1][1],
                             mass.inertia[2][2],
                             mass.inertia[0][1],
                             mass.inertia[0][2],
                             mass.inertia[1][2]};
  std::array<std::uint64_t, 10> bits{};
  for (std::size_t i = 0; i < 10; ++i) {
    std::memcpy(&bits[i], &values[i], sizeof(double));
  }
  return bits;
}

}  // namespace

TEST_CASE("mass properties: canonical paths stay bit-identical") {
  // Pins the exact rectangular Gauss-Legendre output of the untrimmed paths so
  // the trimmed-face work cannot perturb them. The cylinder coefficients come
  // from sin/cos, whose last bits differ between standard libraries; there the
  // exact bits are an MSVC pin and other platforms check the analytic values.
  {
    const geom::Brep cylinder =
        geom::make_cylinder(geom::Point3d{0.0, 0.0, 0.0}, geom::Vec3d{0.0, 0.0, 1.0}, 1.0, 2.0,
                            1e-9)
            .value();
    const auto mass = geom::mass_properties(cylinder);
    REQUIRE(mass.has_value());
#if defined(_MSC_VER)
    const std::array<std::uint64_t, 10> expected = {
        0x401921fb541433a0ull, 0xbdde8ac000000000ull, 0xbcc100b76a374858ull, 0x3ff0000000000006ull,
        0x400d524fe1f92f62ull, 0x400d524fe1ef9732ull, 0x400921fb5401033eull, 0x3ca091760dcabfd0ull,
        0xbcf89b18e670b3aaull, 0xbc70342f650004a5ull};
    CHECK(mass_bits(mass.value()) == expected);
#else
    const double pi = std::numbers::pi_v<double>;
    CHECK(mass.value().volume == doctest::Approx(2.0 * pi).epsilon(1e-8));
    CHECK(std::fabs(mass.value().centroid[0]) < 1e-8);
    CHECK(std::fabs(mass.value().centroid[1]) < 1e-8);
    CHECK(mass.value().centroid[2] == doctest::Approx(1.0).epsilon(1e-8));
    CHECK(mass.value().inertia[0][0] == doctest::Approx(7.0 * pi / 6.0).epsilon(1e-8));
    CHECK(mass.value().inertia[1][1] == doctest::Approx(7.0 * pi / 6.0).epsilon(1e-8));
    CHECK(mass.value().inertia[2][2] == doctest::Approx(pi).epsilon(1e-8));
    CHECK(std::fabs(mass.value().inertia[0][1]) < 1e-8);
    CHECK(std::fabs(mass.value().inertia[0][2]) < 1e-8);
    CHECK(std::fabs(mass.value().inertia[1][2]) < 1e-8);
#endif
  }
  {
    const geom::Brep box =
        geom::make_box(geom::Point3d{0.0, 0.0, 0.0}, geom::Point3d{1.0, 2.0, 3.0}, 1e-9).value();
    const auto mass = geom::mass_properties(box);
    REQUIRE(mass.has_value());
    const std::array<std::uint64_t, 10> expected = {
        0x4018000000000000ull, 0x3fe0000000000000ull, 0x3ff0000000000000ull, 0x3ff8000000000000ull,
        0x401a000000000000ull, 0x4013ffffffffffffull, 0x4004000000000000ull, 0x0000000000000000ull,
        0x0000000000000000ull, 0x0000000000000000ull};
    CHECK(mass_bits(mass.value()) == expected);
  }
}

TEST_CASE("mass properties: a windowed cylinder matches its reference moments") {
  const geom::Brep brep = windowed_cylinder_brep();
  const geom::Result<geom::MassProperties> mass =
      geom::mass_properties(brep, geom::MassPropertiesOptions{1e-9, 24, 1e-5});
  REQUIRE(mass.has_value());
  // Cylinder r=1, h=2 minus the prism {x >= 0.5, |y| <= 0.5, 0.5 <= z <= 1.5}.
  // Closed-form reference (pocket integral by Simpson quadrature); all tensor
  // components are about the centroid. These are within-1e-7 sanity checks
  // against the analytic values; the exact floored output is pinned separately.
  CHECK(mass.value().volume == doctest::Approx(5.8265738297177867).epsilon(1e-7));
  CHECK(mass.value().centroid[0] == doctest::Approx(-0.057209149501289074).epsilon(1e-9));
  CHECK(mass.value().centroid[1] == doctest::Approx(0.0).epsilon(1e-7));
  CHECK(mass.value().centroid[2] == doctest::Approx(1.0).epsilon(1e-7));
  CHECK(mass.value().inertia[0][0] == doctest::Approx(3.5920340332365619).epsilon(1e-7));
  // Parallel-axis check: the full cylinder's centroidal Iyy is 3.6651914, the
  // pocket removes 0.2896637 about its own centroid, and the volume shift
  // contributes V*Cx^2 = 0.0190705; 3.6651914 - 0.2896637 - 0.0190705 =
  // 3.3564579658. Ixx is unaffected (Cy = 0), which is why only Iyy disagreed
  // with the original brief value 3.5729643167320004 (that was the pre-shift,
  // origin-plane value).
  CHECK(mass.value().inertia[1][1] == doctest::Approx(3.3564579658203764).epsilon(1e-7));
  CHECK(mass.value().inertia[2][2] == doctest::Approx(2.8358037071808475).epsilon(1e-7));
  CHECK(mass.value().inertia[0][1] == doctest::Approx(0.0).epsilon(1e-6));
  CHECK(mass.value().inertia[0][2] == doctest::Approx(0.0).epsilon(1e-6));
  CHECK(mass.value().inertia[1][2] == doctest::Approx(0.0).epsilon(1e-6));
  // The tensor must be exactly symmetric.
  CHECK(mass.value().inertia[0][1] == mass.value().inertia[1][0]);
  CHECK(mass.value().inertia[0][2] == mass.value().inertia[2][0]);
  CHECK(mass.value().inertia[1][2] == mass.value().inertia[2][1]);
}

TEST_CASE("mass properties: trimmed output is pinned at the quadrature floor") {
  // The floor is load-bearing: with curved_sampling_deflection below 3.125e-4
  // the interior step stays at 5e-2 rad, so the output is a fixed, reproducible
  // value. These checks are tight enough to fail if the std::max floor is
  // removed (the unfloored 1e-5 step changes the centroid and inertia; the
  // volume is subdivision-invariant).
  const geom::Brep brep = windowed_cylinder_brep();
  const auto mass = geom::mass_properties(brep, geom::MassPropertiesOptions{1e-9, 24, 1e-5});
  REQUIRE(mass.has_value());
  CHECK(mass.value().volume == doctest::Approx(5.826573827129446).epsilon(1e-12));
  CHECK(mass.value().centroid[0] == doctest::Approx(-0.057209149585365).epsilon(1e-12));
  CHECK(mass.value().centroid[1] == doctest::Approx(0.0).epsilon(1e-12));
  CHECK(mass.value().centroid[2] == doctest::Approx(1.0).epsilon(1e-12));
  CHECK(mass.value().inertia[0][0] == doctest::Approx(3.592034030767950).epsilon(1e-12));
  CHECK(mass.value().inertia[1][1] == doctest::Approx(3.356457963297069).epsilon(1e-12));
  CHECK(mass.value().inertia[2][2] == doctest::Approx(2.835803705474689).epsilon(1e-12));
  CHECK(mass.value().inertia[0][1] == doctest::Approx(0.0).epsilon(1e-12));
  CHECK(mass.value().inertia[0][2] == doctest::Approx(0.0).epsilon(1e-12));
  CHECK(mass.value().inertia[1][2] == doctest::Approx(0.0).epsilon(1e-12));
}

TEST_CASE("mass properties: trimmed volume agrees with the volume integrator") {
  const geom::Brep brep = windowed_cylinder_brep();
  const geom::Result<double> reference =
      geom::detail::curved_solid_volume(brep, brep.solid_ids().front(), 1e-9);
  REQUIRE(reference.has_value());
  const auto mass = geom::mass_properties(brep, geom::MassPropertiesOptions{1e-9, 24, 1e-5});
  REQUIRE(mass.has_value());
  CHECK(mass.value().volume == doctest::Approx(reference.value()).epsilon(1e-7));
}

TEST_CASE("mass properties: a scaled trimmed face stays accurate at the default") {
  // r=10, h=5 cylinder minus the seam-crossing window {x >= 5, |y| <= 5,
  // 1 <= z <= 4}: the trim boundary is still rectilinear, so this is the large
  // counterpart of the r=1 fixture. The exact volume is the closed form
  // pi*100*5 - 3*(integral of 2*min(sqrt(100-x^2), 5) for x in [5,10]).
  const geom::Brep brep = scaled_window_brep();
  const geom::Result<double> reference =
      geom::detail::curved_solid_volume(brep, brep.solid_ids().front(), 1e-9);
  REQUIRE(reference.has_value());
  CHECK(reference.value() == doctest::Approx(1433.8128835477412).epsilon(1e-9));
  // 1e-7 is below the floor (which binds for curved_sampling_deflection <=
  // 3.125e-4), so it is accepted but produces exactly the same floored step as
  // 1e-5: the second row documents option acceptance, not a finer mesh.
  for (const double deflection : {1e-5, 1e-7}) {
    const auto mass =
        geom::mass_properties(brep, geom::MassPropertiesOptions{1e-9, 24, deflection});
    REQUIRE(mass.has_value());
    const double relative =
        std::abs(mass.value().volume - reference.value()) / reference.value();
    CHECK(relative < 1e-9);
  }
}

TEST_CASE("mass properties: trimmed sampling deflection stays at the quadrature floor") {
  const geom::Brep brep = windowed_cylinder_brep();
  const geom::Result<double> reference =
      geom::detail::curved_solid_volume(brep, brep.solid_ids().front(), 1e-9);
  REQUIRE(reference.has_value());
  // The uv edges of this rectilinear trimmed domain are exact parameter lines,
  // so the deflection only controls the interior quadrature density: 1e-5 sits at
  // the floor (5e-2 rad); 1e-3 is above the floor (step 8.9e-2 rad) and both
  // stay within 2.6e-9 absolute of the reference.
  const auto fine = geom::mass_properties(brep, geom::MassPropertiesOptions{1e-9, 24, 1e-3});
  const auto finer = geom::mass_properties(brep, geom::MassPropertiesOptions{1e-9, 24, 1e-5});
  REQUIRE(fine.has_value());
  REQUIRE(finer.has_value());
  CHECK(std::abs(fine.value().volume - reference.value()) <= 1e-8 * reference.value());
  CHECK(std::abs(finer.value().volume - reference.value()) <= 1e-8 * reference.value());
  // A deliberately coarse deflection does coarsen the interior quadrature; only
  // there does the expected error ordering become observable.
  const auto coarse = geom::mass_properties(brep, geom::MassPropertiesOptions{1e-9, 24, 1e-1});
  REQUIRE(coarse.has_value());
  const double coarse_error = std::abs(coarse.value().volume - reference.value());
  const double fine_error = std::abs(fine.value().volume - reference.value());
  CHECK(fine_error < coarse_error);
}

TEST_CASE("mass properties: trimmed dents and holes report their limitation") {
  {
    // The window at +y (u in [45, 135] deg) does not cross the seam, so it
    // closes in uv and the trimmed face has an outer ring plus a hole loop; the
    // parametric integration does not support holes yet (deferred to M3c-2).
    const geom::Brep brep = hole_window_brep();
    const auto mass = geom::mass_properties(brep, geom::MassPropertiesOptions{1e-9, 24, 1e-5});
    REQUIRE(!mass.has_value());
    CHECK(mass.error().code == geom::Errc::unsupported);
    CHECK(std::string(mass.error().message) ==
          "mass_properties: trimmed curved face has a hole");
  }
  {
    // An oblique cut leaves an ellipse in uv, which is not rectilinear.
    const geom::Brep brep = slanted_trim_brep();
    const auto mass = geom::mass_properties(brep, geom::MassPropertiesOptions{1e-9, 24, 1e-5});
    REQUIRE(!mass.has_value());
    CHECK(mass.error().code == geom::Errc::unsupported);
    CHECK(std::string(mass.error().message) ==
          "mass_properties: trimmed curved face boundary is not rectilinear");
  }
}

TEST_CASE("mass properties: trimmed sampling deflection is validated") {
  // A plain full cylinder has no trimmed face, so the option check cannot be
  // masked by a face_uv_projection error.
  const geom::Brep brep =
      geom::make_cylinder(geom::Point3d{0.0, 0.0, 0.0}, geom::Vec3d{0.0, 0.0, 1.0}, 1.0, 2.0,
                          1e-9)
          .value();
  for (const double deflection : {0.0, -1.0, std::numeric_limits<double>::quiet_NaN()}) {
    const auto mass = geom::mass_properties(brep, geom::MassPropertiesOptions{1e-9, 24, deflection});
    REQUIRE(!mass.has_value());
    CHECK(mass.error().code == geom::Errc::tolerance_violation);
    CHECK(std::string(mass.error().message) == "mass_properties: options are out of range");
  }
}
