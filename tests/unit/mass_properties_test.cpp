#include <doctest/doctest.h>

#include <array>
#include <cmath>
#include <cstdint>
#include <limits>
#include <numbers>
#include <utility>

#include <geom/brep/mass_properties.hpp>
#include <geom/brep/primitives.hpp>
#include <geom/brep/solid_editor.hpp>
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
