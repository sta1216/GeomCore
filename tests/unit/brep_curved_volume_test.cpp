#include <doctest/doctest.h>

#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <numbers>
#include <string_view>
#include <variant>
#include <vector>

#include <geom/brep/detail/curved_volume.hpp>
#include <geom/brep/face_uv.hpp>
#include <geom/brep/primitives.hpp>
#include <geom/brep/solid_editor.hpp>
#include <geom/brep/topology.hpp>
#include <geom/core/error.hpp>
#include <geom/surfaces/nurbs_surface.hpp>
#include <geom/surfaces/plane.hpp>

#include "brep_curved_fixture.hpp"

namespace {

constexpr double kPi = std::numbers::pi_v<double>;
constexpr double kTwoPi = 2.0 * kPi;
constexpr double kTolerance = 1e-9;

double relative_error(double actual, double expected) {
  return std::abs(actual - expected) / std::abs(expected);
}

void check_error(const geom::Result<double>& result, geom::Errc code, std::string_view message) {
  REQUIRE_FALSE(result.has_value());
  CHECK(result.error().code == code);
  CHECK(std::string_view(result.error().message) == message);
}

geom::FaceId top_plane_face(const geom::Brep& brep) {
  for (const geom::FaceId id : brep.face_ids()) {
    const geom::Face& face = brep.face(id);
    if (face.surface == geom::kNoGeometry) continue;
    const auto& surface = brep.geometry().surfaces[face.surface];
    if (!std::holds_alternative<geom::Plane3>(surface)) continue;
    const geom::Plane3& plane = std::get<geom::Plane3>(surface);
    if (plane.normal[2] > 0.5) return id;
  }
  return geom::FaceId{};
}

// A one-face, topologically closed shell whose only face is a bilinear NURBS
// patch. curved_solid_volume must reject the face before it attempts to
// integrate it.
geom::Brep make_nurbs_face_solid() {
  geom::SolidEditor editor;
  const std::uint32_t surface = editor.add_surface(geom::NurbsSurface(
      1, 1, {0.0, 0.0, 1.0, 1.0}, {0.0, 0.0, 1.0, 1.0},
      {{0.0, 0.0, 0.0}, {1.0, 0.0, 0.0}, {0.0, 1.0, 0.0}, {1.0, 1.0, 0.0}},
      {1.0, 1.0, 1.0, 1.0}));

  const std::array<geom::Point3d, 3> points{
      geom::Point3d{0.0, 0.0, 0.0},
      geom::Point3d{1.0, 0.0, 0.0},
      geom::Point3d{0.0, 1.0, 0.0},
  };
  std::array<geom::VertexId, 3> vertices{};
  for (std::size_t i = 0; i < vertices.size(); ++i) {
    vertices[i] = editor.add_vertex(points[i], 0.0).value();
  }
  std::array<geom::EdgeId, 3> edges{};
  for (std::size_t i = 0; i < edges.size(); ++i) {
    edges[i] = editor
                   .add_edge(vertices[i], vertices[(i + 1) % vertices.size()],
                             geom::kNoGeometry, 0.0, 1.0, 0.0)
                   .value();
  }

  const geom::FaceId face = editor.add_face(surface, false, 0.0).value();
  std::array<geom::CoedgeSpec, 3> coedges{};
  for (std::size_t i = 0; i < coedges.size(); ++i) {
    coedges[i] = geom::CoedgeSpec{edges[i], false, geom::kNoGeometry, 0.0};
  }
  editor.add_loop(face, coedges, true).value();

  const geom::ShellId shell = editor.add_shell(std::array{face}, true).value();
  editor.add_solid(std::array{shell}, 0.0).value();
  return editor.brep();
}

}  // namespace

TEST_CASE("curved_solid_volume: box") {
  const geom::Brep brep = geom::make_box({0.0, 0.0, 0.0}, {2.0, 2.0, 2.0}, kTolerance).value();
  const geom::Result<double> volume =
      geom::detail::curved_solid_volume(brep, brep.solid_ids().front(), kTolerance);
  REQUIRE(volume.has_value());
  CHECK(relative_error(volume.value(), 8.0) <= 1e-12);
}

TEST_CASE("curved_solid_volume: cylinder") {
  constexpr double radius = 2.0;
  constexpr double height = 3.0;
  const geom::Brep brep =
      geom::make_cylinder({0.0, 0.0, 0.0}, {0.0, 0.0, 1.0}, radius, height, kTolerance).value();
  const geom::FaceId cap = top_plane_face(brep);
  REQUIRE_FALSE(cap.is_null());

  const geom::Result<std::vector<geom::FaceUvLoop>> loops =
      geom::face_loops_uv(brep, cap, kTolerance);
  REQUIRE(loops.has_value());
  REQUIRE(loops.value().size() == 1);
  REQUIRE(loops.value().front().outer);
  const std::size_t segments = loops.value().front().points.size();
  CHECK(segments == 32);

  const geom::Result<double> volume =
      geom::detail::curved_solid_volume(brep, brep.solid_ids().front(), kTolerance);
  REQUIRE(volume.has_value());

  // The caps are the 32-gon used by face_loops_uv, but the lateral face is
  // integrated on the exact cylinder. The prescribed quadrature therefore
  // returns h * (2*pi*r^2 + A) / 3 with A the inscribed 32-gon area.
  const double n = static_cast<double>(segments);
  const double cap_area = 0.5 * n * radius * radius * std::sin(kTwoPi / n);
  const double expected = height * (kTwoPi * radius * radius + cap_area) / 3.0;
  CHECK(relative_error(volume.value(), expected) <= 1e-9);

  // 32 chords make the sampled solid about 0.21% smaller than the true
  // cylinder, so the brief's 1e-3 analytic bound is unattainable here.
  CHECK(relative_error(volume.value(), kPi * radius * radius * height) <= 3e-3);
}

TEST_CASE("curved_solid_volume: quarter cylinder") {
  constexpr double radius = 1.0;
  constexpr double height = 1.0;
  // The fixture commits a closed shell but no solid (and must not be edited),
  // so attach the shell to a solid here.
  geom::SolidEditor editor(geom::test::make_quarter_cylinder(radius, height, kTolerance));
  REQUIRE(editor.brep().shell_count() == 1);
  const geom::ShellId shell = editor.brep().shell_ids().front();
  editor.add_solid(std::array{shell}, kTolerance).value();
  const geom::Brep brep = editor.commit().value();
  const geom::FaceId lateral = geom::test::quarter_cylinder_lateral(brep);
  REQUIRE_FALSE(lateral.is_null());

  const geom::Result<double> volume =
      geom::detail::curved_solid_volume(brep, brep.solid_ids().front(), kTolerance);
  REQUIRE(volume.has_value());

  // The fixture has 16 chords per quarter arc. The lateral quarter stays a
  // smooth cylinder patch, so the result is the half-cylinder lateral term
  // plus the inscribed 16-gon top-cap term.
  constexpr double n = 16.0;
  const double cap_area = 0.5 * n * radius * radius * std::sin((kPi / 2.0) / n);
  const double expected = height * (kPi * radius * radius / 2.0 + cap_area) / 3.0;
  CHECK(relative_error(volume.value(), expected) <= 1e-9);
  CHECK(relative_error(volume.value(), kPi * radius * radius * height / 4.0) <= 2e-3);
}

TEST_CASE("curved_solid_volume: error paths") {
  const geom::Brep brep = geom::make_box({0.0, 0.0, 0.0}, {2.0, 2.0, 2.0}, kTolerance).value();
  const geom::SolidId solid = brep.solid_ids().front();

  check_error(geom::detail::curved_solid_volume(brep, geom::SolidId{}, kTolerance),
              geom::Errc::index_out_of_range,
              "curved_solid_volume: solid is not alive");

  const geom::Brep empty;
  check_error(geom::detail::curved_solid_volume(empty, geom::SolidId::make(0, 1), kTolerance),
              geom::Errc::index_out_of_range,
              "curved_solid_volume: solid is not alive");

  check_error(geom::detail::curved_solid_volume(brep, solid, 0.0),
              geom::Errc::tolerance_violation,
              "curved_solid_volume: tolerance must be finite and positive");
  check_error(geom::detail::curved_solid_volume(brep, solid,
                                                std::numeric_limits<double>::quiet_NaN()),
              geom::Errc::tolerance_violation,
              "curved_solid_volume: tolerance must be finite and positive");
}

TEST_CASE("curved_solid_volume: rejects NURBS faces") {
  const geom::Brep brep = make_nurbs_face_solid();
  REQUIRE(brep.solid_count() == 1);
  check_error(geom::detail::curved_solid_volume(brep, brep.solid_ids().front(), kTolerance),
              geom::Errc::unsupported,
              "curved_solid_volume: NURBS surfaces are not supported");
}
