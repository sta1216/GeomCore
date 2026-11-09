#include <array>
#include <cstdint>
#include <string_view>
#include <utility>

#include <doctest/doctest.h>
#include <geom/brep/boolean_solid.hpp>
#include <geom/brep/face_query.hpp>
#include <geom/brep/mass_properties.hpp>
#include <geom/brep/primitives.hpp>
#include <geom/brep/solid_editor.hpp>
#include <geom/brep/validate.hpp>
#include <geom/core/error.hpp>

using namespace geom;

namespace {

template <class T>
T take(const Result<T>& result) {
  REQUIRE(result.has_value());
  return result.value();
}

// Adds one closed box shell to `editor` and returns it. The six faces get their
// own Plane3 surface with the stored normal pointing out of the box, and each
// loop winds counter-clockwise about that normal; `inward` negates the normal
// and reverses every loop so the shell encloses a negative signed volume, which
// is how a cavity shell is represented.
ShellId add_box_shell(SolidEditor& editor, Point3d lo, Point3d hi, bool inward) {
  const Point3d corners[8] = {
      {lo[0], lo[1], lo[2]}, {hi[0], lo[1], lo[2]}, {hi[0], hi[1], lo[2]},
      {lo[0], hi[1], lo[2]}, {lo[0], lo[1], hi[2]}, {hi[0], lo[1], hi[2]},
      {hi[0], hi[1], hi[2]}, {lo[0], hi[1], hi[2]},
  };
  std::array<VertexId, 8> vertices{};
  for (int i = 0; i < 8; ++i) {
    vertices[i] = editor.add_vertex(corners[i], 0.0).value();
  }
  const int edge_specs[12][2] = {{0, 1},  {1, 2},  {2, 3},  {3, 0},
                                 {4, 5},  {5, 6},  {6, 7},  {7, 4},
                                 {1, 5},  {2, 6},  {3, 7},  {0, 4}};
  std::array<EdgeId, 12> edges{};
  for (int i = 0; i < 12; ++i) {
    edges[i] = editor
                   .add_edge(vertices[edge_specs[i][0]], vertices[edge_specs[i][1]],
                             kNoGeometry, 0.0, 1.0, 0.0)
                   .value();
  }

  const Vec3d normals[6] = {{0, 0, -1}, {0, 0, 1}, {0, -1, 0},
                            {0, 1, 0},  {-1, 0, 0}, {1, 0, 0}};
  const int corner_index[6] = {0, 4, 0, 2, 0, 1};
  const int face_edges[6][4] = {{3, 2, 1, 0},   {4, 5, 6, 7},  {0, 8, 4, 11},
                                {2, 10, 6, 9},  {11, 7, 10, 3}, {1, 9, 5, 8}};
  const bool face_reversed[6][4] = {{true, true, true, true},
                                    {false, false, false, false},
                                    {false, false, true, true},
                                    {false, false, true, true},
                                    {false, true, true, false},
                                    {false, false, true, true}};
  std::array<FaceId, 6> faces{};
  for (int f = 0; f < 6; ++f) {
    Vec3d normal = normals[f];
    if (inward) normal = Vec3d{-normal[0], -normal[1], -normal[2]};
    const Point3d& corner = corners[corner_index[f]];
    const double offset = dot(normal, Vec3d{corner[0], corner[1], corner[2]});
    const std::uint32_t surface = editor.add_surface(Plane3{normal, offset});
    faces[f] = editor.add_face(surface, false, 0.0).value();
    std::array<CoedgeSpec, 4> specs{};
    for (int j = 0; j < 4; ++j) {
      const int source = inward ? 3 - j : j;
      specs[j] = CoedgeSpec{edges[face_edges[f][source]],
                            inward ? !face_reversed[f][source]
                                   : face_reversed[f][source]};
    }
    editor.add_loop(faces[f], specs, true).value();
  }
  return editor.add_shell(faces, true).value();
}

}  // namespace

TEST_CASE("a fully interior difference opens a cavity") {
  const Brep outer = take(make_box({0, 0, 0}, {10, 10, 10}, 1e-12));
  const Brep inner = take(make_box({3, 3, 3}, {7, 7, 7}, 1e-12));
  const auto result = boolean_solid(outer, inner, SolidBooleanOp::difference);
  REQUIRE(result.has_value());
  const Brep& brep = result.value().brep;

  CHECK(brep.shell_count() == 2u);
  CHECK(brep.solid_count() == 1u);
  CHECK(validate_topology(brep).ok());
  // The outer shell is emitted first and the cavity shell after it.
  const auto shells = brep.shell_ids();
  REQUIRE(shells.size() == 2u);
  const Result<bool> outer_contains_corner =
      point_in_shell(brep, shells[0], Point3d{1, 1, 1}, 1e-9);
  REQUIRE(outer_contains_corner.has_value());
  CHECK(outer_contains_corner.value());
  const Result<bool> cavity_contains_corner =
      point_in_shell(brep, shells[1], Point3d{1, 1, 1}, 1e-9);
  REQUIRE(cavity_contains_corner.has_value());
  CHECK_FALSE(cavity_contains_corner.value());
  CHECK(mass_properties(brep).value().volume ==
        doctest::Approx(936.0).epsilon(1e-9));
  CHECK(euler_characteristic(brep, brep.solid_ids().front()) == 4);
}

TEST_CASE("a through-hole keeps one shell") {
  const Brep plate = take(make_box({0, 0, 0}, {10, 10, 2}, 1e-12));
  const Brep punch = take(make_box({4, 4, -1}, {6, 6, 3}, 1e-12));
  const auto result = boolean_solid(plate, punch, SolidBooleanOp::difference);
  REQUIRE(result.has_value());
  const Brep& brep = result.value().brep;

  CHECK(brep.shell_count() == 1u);
  CHECK(validate_topology(brep).ok());
  CHECK(mass_properties(brep).value().volume ==
        doctest::Approx(192.0).epsilon(1e-9));
  CHECK(euler_characteristic(brep, brep.solid_ids().front()) == 0);
}

TEST_CASE("a solid with two outer shells is rejected") {
  SolidEditor editor;
  const ShellId first = add_box_shell(editor, {0, 0, 0}, {1, 1, 1}, false);
  const ShellId second = add_box_shell(editor, {2, 2, 2}, {3, 3, 3}, false);
  REQUIRE(editor.add_solid(std::array<ShellId, 2>{first, second}, 0.0).has_value());

  const ValidationReport report = validate_topology(editor.brep());
  REQUIRE_FALSE(report.ok());
  CHECK(report.issues.front().code == Errc::invalid_topology);
  CHECK(std::string_view(report.issues.front().message) ==
        "solid has more than one outer shell");
}

TEST_CASE("a cavity outside its outer shell is rejected") {
  SolidEditor editor;
  const ShellId outer = add_box_shell(editor, {0, 0, 0}, {2, 2, 2}, false);
  const ShellId cavity = add_box_shell(editor, {3, 3, 3}, {4, 4, 4}, true);
  REQUIRE(editor.add_solid(std::array<ShellId, 2>{outer, cavity}, 0.0).has_value());

  const ValidationReport report = validate_topology(editor.brep());
  REQUIRE_FALSE(report.ok());
  CHECK(report.issues.front().code == Errc::invalid_topology);
  CHECK(std::string_view(report.issues.front().message) ==
        "cavity shell is not contained in the outer shell");
}

TEST_CASE("overlapping cavity shells are rejected") {
  SolidEditor editor;
  const ShellId outer = add_box_shell(editor, {0, 0, 0}, {10, 10, 10}, false);
  const ShellId first = add_box_shell(editor, {1, 1, 1}, {6, 6, 6}, true);
  const ShellId second = add_box_shell(editor, {5, 5, 5}, {9, 9, 9}, true);
  REQUIRE(editor.add_solid(std::array<ShellId, 3>{outer, first, second}, 0.0)
              .has_value());

  const ValidationReport report = validate_topology(editor.brep());
  REQUIRE_FALSE(report.ok());
  CHECK(report.issues.front().code == Errc::invalid_topology);
  CHECK(std::string_view(report.issues.front().message) ==
        "cavity shells overlap");
}

TEST_CASE("a contained cavity shell is accepted") {
  SolidEditor editor;
  const ShellId outer = add_box_shell(editor, {0, 0, 0}, {4, 4, 4}, false);
  const ShellId cavity = add_box_shell(editor, {1, 1, 1}, {3, 3, 3}, true);
  REQUIRE(editor.add_solid(std::array<ShellId, 2>{outer, cavity}, 0.0).has_value());

  const Brep& brep = editor.brep();
  CHECK(validate_topology(brep).ok());
  CHECK(validate_tolerances(brep).ok());
  CHECK(mass_properties(brep).value().volume ==
        doctest::Approx(56.0).epsilon(1e-12));
  CHECK(euler_characteristic(brep, brep.solid_ids().front()) == 4);
}

TEST_CASE("crossing cavity shells are rejected") {
  // The two cavities intersect, but every sample point and vertex of one lies
  // outside the other, so only the bounds test catches the crossing.
  SolidEditor editor;
  const ShellId outer = add_box_shell(editor, {0, 0, 0}, {100, 100, 100}, false);
  const ShellId first =
      add_box_shell(editor, {10, 10, 10}, {20, 11, 11}, true);
  const ShellId second =
      add_box_shell(editor, {12, 0.5, 10.2}, {12.5, 30, 10.5}, true);
  REQUIRE(editor.add_solid(std::array<ShellId, 3>{outer, first, second}, 0.0)
              .has_value());

  const ValidationReport report = validate_topology(editor.brep());
  REQUIRE_FALSE(report.ok());
  CHECK(report.issues.front().code == Errc::invalid_topology);
  CHECK(std::string_view(report.issues.front().message) ==
        "cavity shells overlap");
}

TEST_CASE("a zero-volume cavity shell is rejected") {
  SolidEditor editor;
  const ShellId outer = add_box_shell(editor, {0, 0, 0}, {4, 4, 4}, false);
  const ShellId degenerate = add_box_shell(editor, {2, 2, 2}, {2, 2, 2}, true);
  REQUIRE(editor.add_solid(std::array<ShellId, 2>{outer, degenerate}, 0.0)
              .has_value());

  const ValidationReport report = validate_topology(editor.brep());
  REQUIRE_FALSE(report.ok());
  CHECK(report.issues.front().code == Errc::invalid_topology);
  CHECK(std::string_view(report.issues.front().message) ==
        "shell has zero volume");
}

TEST_CASE("a contained cavity is accepted at zero tolerance") {
  // Non-round cavity coordinates and a zero solid tolerance exercise the probe
  // slack: without it the point queries can report an ambiguous ray hit.
  SolidEditor editor;
  const ShellId outer = add_box_shell(editor, {0, 0, 0}, {100, 100, 100}, false);
  const ShellId cavity =
      add_box_shell(editor, {10, 10, 10}, {20, 11, 11}, true);
  REQUIRE(editor.add_solid(std::array<ShellId, 2>{outer, cavity}, 0.0).has_value());

  const ValidationReport report = validate_topology(editor.brep());
  CHECK(report.ok());
  CHECK(mass_properties(editor.brep()).value().volume ==
        doctest::Approx(1000000.0 - 10.0).epsilon(1e-15));
}
