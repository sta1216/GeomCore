#include <doctest/doctest.h>

#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <iomanip>
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
#include <geom/curves/analytic_curve.hpp>
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

// make_box() always builds a positively oriented solid and SolidEditor has no
// setter for Face::reversed, so this repeats the make_box construction with
// the top face marked reversed. Topology and geometry are otherwise identical
// to make_box({0,0,0}, {2,2,2}, kTolerance).
geom::Brep make_reversed_top_box() {
  constexpr double side = 2.0;
  const geom::Point3d corners[8] = {
      {0.0, 0.0, 0.0}, {side, 0.0, 0.0}, {side, side, 0.0}, {0.0, side, 0.0},
      {0.0, 0.0, side}, {side, 0.0, side}, {side, side, side}, {0.0, side, side},
  };

  geom::SolidEditor editor;
  geom::VertexId vertices[8];
  for (int i = 0; i < 8; ++i) {
    vertices[i] = editor.add_vertex(corners[i], kTolerance).value();
  }

  struct EdgeSpec {
    int v0;
    int v1;
  };
  const EdgeSpec edge_specs[12] = {
      {0, 1}, {1, 2}, {2, 3}, {3, 0},  // bottom ring
      {4, 5}, {5, 6}, {6, 7}, {7, 4},  // top ring
      {1, 5}, {2, 6}, {3, 7}, {0, 4},  // verticals
  };
  geom::EdgeId edges[12];
  for (int i = 0; i < 12; ++i) {
    const geom::Point3d from = corners[edge_specs[i].v0];
    const geom::Point3d to = corners[edge_specs[i].v1];
    const geom::Vec3d delta = to - from;
    const double length = std::hypot(delta[0], delta[1], delta[2]);
    const geom::Vec3d direction{delta[0] / length, delta[1] / length, delta[2] / length};
    const std::uint32_t curve = editor.add_curve(geom::LineCurve3{from, direction});
    edges[i] = editor
                   .add_edge(vertices[edge_specs[i].v0], vertices[edge_specs[i].v1], curve,
                             0.0, length, kTolerance)
                   .value();
  }

  struct CoedgeUse {
    int edge;
    bool reversed;
  };
  struct FaceSpec {
    geom::Vec3d normal;
    int corner;
    CoedgeUse loop[4];
  };
  // Loop order is counter-clockwise seen from outside along `normal`.
  const FaceSpec face_specs[6] = {
      {{0.0, 0.0, -1.0}, 0, {{3, true}, {2, true}, {1, true}, {0, true}}},    // bottom
      {{0.0, 0.0, 1.0}, 4, {{4, false}, {5, false}, {6, false}, {7, false}}},  // top
      {{0.0, -1.0, 0.0}, 0, {{0, false}, {8, false}, {4, true}, {11, true}}},  // front
      {{0.0, 1.0, 0.0}, 2, {{2, false}, {10, false}, {6, true}, {9, true}}},   // back
      {{-1.0, 0.0, 0.0}, 0, {{11, false}, {7, true}, {10, true}, {3, false}}}, // left
      {{1.0, 0.0, 0.0}, 1, {{1, false}, {9, false}, {5, true}, {8, true}}},    // right
  };

  std::array<geom::FaceId, 6> faces{};
  for (int f = 0; f < 6; ++f) {
    const FaceSpec& spec = face_specs[f];
    const geom::Point3d corner = corners[spec.corner];
    const std::uint32_t surface = editor.add_surface(geom::Plane3{
        spec.normal,
        geom::dot(spec.normal, geom::Vec3d{corner[0], corner[1], corner[2]})});
    // Only the top face is marked reversed.
    faces[static_cast<std::size_t>(f)] =
        editor.add_face(surface, f == 1, kTolerance).value();

    std::array<geom::CoedgeSpec, 4> coedges{};
    for (int i = 0; i < 4; ++i) {
      coedges[static_cast<std::size_t>(i)] =
          geom::CoedgeSpec{edges[spec.loop[i].edge], spec.loop[i].reversed};
    }
    editor.add_loop(faces[static_cast<std::size_t>(f)], coedges, true).value();
  }

  const geom::ShellId shell = editor.add_shell(faces, true).value();
  editor.add_solid(std::array{shell}, kTolerance).value();
  return editor.commit().value();
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

  // The caps carry exact circle provenance, so the arc-aware planar area uses
  // pi*r^2 rather than the 32-gon; the lateral face is integrated on the exact
  // cylinder. The result must therefore equal the analytic volume.
  CHECK(relative_error(volume.value(), kPi * radius * radius * height) <= 1e-9);
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

TEST_CASE("curved_solid_volume: translated cylinder") {
  constexpr double radius = 2.0;
  constexpr double height = 3.0;
  const geom::Brep origin_brep =
      geom::make_cylinder({0.0, 0.0, 0.0}, {0.0, 0.0, 1.0}, radius, height, kTolerance).value();
  const geom::Brep translated_brep =
      geom::make_cylinder({1000.0, 0.0, 0.0}, {0.0, 0.0, 1.0}, radius, height, kTolerance).value();

  const geom::Result<double> origin_volume =
      geom::detail::curved_solid_volume(origin_brep, origin_brep.solid_ids().front(), kTolerance);
  const geom::Result<double> translated_volume = geom::detail::curved_solid_volume(
      translated_brep, translated_brep.solid_ids().front(), kTolerance);
  REQUIRE(origin_volume.has_value());
  REQUIRE(translated_volume.has_value());

  MESSAGE("origin cylinder volume: " << std::setprecision(17) << origin_volume.value());
  MESSAGE("translated cylinder volume: " << std::setprecision(17) << translated_volume.value());

  // Quadrature over the raw uv triangles used to depend on where the solid sat.
  CHECK(relative_error(translated_volume.value(), origin_volume.value()) <= 1e-9);
  CHECK(relative_error(translated_volume.value(), kPi * radius * radius * height) <= 1e-9);
}

TEST_CASE("curved_solid_volume: box far from the origin") {
  constexpr double base = 1e9;
  const geom::Brep brep =
      geom::make_box({base, base, base}, {base + 2.0, base + 2.0, base + 2.0}, kTolerance)
          .value();
  const geom::Result<double> volume =
      geom::detail::curved_solid_volume(brep, brep.solid_ids().front(), kTolerance);
  REQUIRE(volume.has_value());
  MESSAGE("box at 1e9 volume: " << std::setprecision(17) << volume.value());

  // The absolute plane offsets are 1e9; the uv polygon shift plus the reference
  // point of the divergence integrand keep the two large opposite face terms
  // from cancelling the 8 away.
  CHECK(relative_error(volume.value(), 8.0) <= 1e-9);
}

TEST_CASE("curved_solid_volume: torus") {
  constexpr double major_radius = 2.0;
  constexpr double minor_radius = 0.5;
  const geom::Brep brep =
      geom::make_torus({0.0, 0.0, 0.0}, {0.0, 0.0, 1.0}, major_radius, minor_radius, kTolerance)
          .value();
  const geom::Result<double> volume =
      geom::detail::curved_solid_volume(brep, brep.solid_ids().front(), kTolerance);
  REQUIRE(volume.has_value());
  MESSAGE("torus volume: " << std::setprecision(17) << volume.value());

  const double expected =
      2.0 * kPi * kPi * major_radius * minor_radius * minor_radius;
  CHECK(relative_error(volume.value(), expected) <= 1e-8);
}

TEST_CASE("curved_solid_volume: reversed top face") {
  const geom::Brep brep = make_reversed_top_box();
  const geom::Result<double> volume =
      geom::detail::curved_solid_volume(brep, brep.solid_ids().front(), kTolerance);
  REQUIRE(volume.has_value());
  MESSAGE("reversed top box volume: " << std::setprecision(17) << volume.value());

  // The six outward faces give 8; reversing the top face (which contributes
  // 2 * 4 / 3 = 8/3) subtracts it twice.
  CHECK(relative_error(volume.value(), 8.0 / 3.0) <= 1e-9);
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
