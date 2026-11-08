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
