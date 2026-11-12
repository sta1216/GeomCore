#include <algorithm>
#include <array>
#include <cstddef>
#include <cstdint>
#include <string_view>
#include <vector>

#include <doctest/doctest.h>
#include <geom/brep/heal.hpp>
#include <geom/brep/primitives.hpp>
#include <geom/brep/solid_editor.hpp>

using namespace geom;

namespace {

// Adds four vertices and four edges of an axis-aligned square in the z = 0
// plane, in counter-clockwise order, and returns the edges in that order.
std::array<EdgeId, 4> add_quad_edges(SolidEditor& editor, double x, double y, double side) {
  const auto v00 = editor.add_vertex({x, y, 0.0}, 0).value();
  const auto v10 = editor.add_vertex({x + side, y, 0.0}, 0).value();
  const auto v11 = editor.add_vertex({x + side, y + side, 0.0}, 0).value();
  const auto v01 = editor.add_vertex({x, y + side, 0.0}, 0).value();
  return {editor.add_edge(v00, v10, kNoGeometry, 0, 1, 0).value(),
          editor.add_edge(v10, v11, kNoGeometry, 0, 1, 0).value(),
          editor.add_edge(v11, v01, kNoGeometry, 0, 1, 0).value(),
          editor.add_edge(v01, v00, kNoGeometry, 0, 1, 0).value()};
}

// Adds one face with a single outer loop, one open shell and one solid over
// `edges`, and returns the solid.
SolidId add_open_sheet(SolidEditor& editor, std::uint32_t surface,
                       const std::array<EdgeId, 4>& edges) {
  const auto face = editor.add_face(surface, false, 0).value();
  const CoedgeSpec loop[] = {
      {edges[0], false}, {edges[1], false}, {edges[2], false}, {edges[3], false}};
  REQUIRE(editor.add_loop(face, loop, true).has_value());
  const auto shell = editor.add_shell(std::array<FaceId, 1>{face}, false).value();
  return editor.add_solid(std::array<ShellId, 1>{shell}, 0).value();
}

}  // namespace

TEST_CASE("heal: boundary edges of a closed box are empty") {
  const Brep box = make_box({0, 0, 0}, {1, 1, 1}, 0.0).value();
  const SolidId solid = box.solid_ids().front();
  CHECK(detail::heal_boundary_edges(box, solid).empty());
}

TEST_CASE("heal: an open quad sheet exposes four boundary edges") {
  SolidEditor editor;
  const auto s = editor.add_surface(Plane3{{0.0, 0.0, 1.0}, 0.0});
  const auto v00 = editor.add_vertex({0, 0, 0}, 0).value();
  const auto v10 = editor.add_vertex({1, 0, 0}, 0).value();
  const auto v11 = editor.add_vertex({1, 1, 0}, 0).value();
  const auto v01 = editor.add_vertex({0, 1, 0}, 0).value();
  const auto e0 = editor.add_edge(v00, v10, kNoGeometry, 0, 1, 0).value();
  const auto e1 = editor.add_edge(v10, v11, kNoGeometry, 0, 1, 0).value();
  const auto e2 = editor.add_edge(v11, v01, kNoGeometry, 0, 1, 0).value();
  const auto e3 = editor.add_edge(v01, v00, kNoGeometry, 0, 1, 0).value();
  const auto face = editor.add_face(s, false, 0).value();
  const CoedgeSpec loop[] = {{e0, false}, {e1, false}, {e2, false}, {e3, false}};
  REQUIRE(editor.add_loop(face, loop, true).has_value());
  const std::array<FaceId, 1> faces{face};
  const auto shell = editor.add_shell(faces, false).value();
  const std::array<ShellId, 1> shells{shell};
  REQUIRE(editor.add_solid(shells, 0).has_value());
  const Brep brep = editor.commit().value();

  const std::vector<EdgeId> boundary =
      detail::heal_boundary_edges(brep, brep.solid_ids().front());
  REQUIRE(boundary.size() == 4u);
  // The result is sorted ascending by edge index.
  for (std::size_t i = 1; i < boundary.size(); ++i) {
    CHECK(boundary[i - 1].index() < boundary[i].index());
  }
}

TEST_CASE("heal: boundary discovery walks every loop of a face") {
  SolidEditor editor;
  const auto s = editor.add_surface(Plane3{{0.0, 0.0, 1.0}, 0.0});
  const std::array<EdgeId, 4> outer_edges = add_quad_edges(editor, 0.0, 0.0, 1.0);
  const std::array<EdgeId, 4> hole_edges = add_quad_edges(editor, 0.25, 0.25, 0.5);
  const auto face = editor.add_face(s, false, 0).value();
  const CoedgeSpec outer_loop[] = {{outer_edges[0], false},
                                   {outer_edges[1], false},
                                   {outer_edges[2], false},
                                   {outer_edges[3], false}};
  const CoedgeSpec hole_loop[] = {{hole_edges[0], false},
                                  {hole_edges[1], false},
                                  {hole_edges[2], false},
                                  {hole_edges[3], false}};
  REQUIRE(editor.add_loop(face, outer_loop, true).has_value());
  REQUIRE(editor.add_loop(face, hole_loop, false).has_value());
  const auto shell = editor.add_shell(std::array<FaceId, 1>{face}, false).value();
  REQUIRE(editor.add_solid(std::array<ShellId, 1>{shell}, 0).has_value());
  const Brep brep = editor.commit().value();

  const std::vector<EdgeId> boundary =
      detail::heal_boundary_edges(brep, brep.solid_ids().front());
  std::vector<std::uint32_t> expected;
  for (const EdgeId id : outer_edges) expected.push_back(id.index());
  for (const EdgeId id : hole_edges) expected.push_back(id.index());
  std::sort(expected.begin(), expected.end());
  REQUIRE(boundary.size() == expected.size());
  for (std::size_t i = 0; i < expected.size(); ++i) {
    CHECK(boundary[i].index() == expected[i]);
  }
}

TEST_CASE("heal: boundary discovery is restricted to the target solid") {
  SolidEditor editor;
  const auto s = editor.add_surface(Plane3{{0.0, 0.0, 1.0}, 0.0});
  const std::array<EdgeId, 4> first_edges = add_quad_edges(editor, 0.0, 0.0, 1.0);
  const std::array<EdgeId, 4> second_edges = add_quad_edges(editor, 5.0, 0.0, 1.0);
  const SolidId first = add_open_sheet(editor, s, first_edges);
  const SolidId second = add_open_sheet(editor, s, second_edges);
  const Brep brep = editor.commit().value();
  REQUIRE(first != second);

  const std::vector<EdgeId> first_boundary = detail::heal_boundary_edges(brep, first);
  REQUIRE(first_boundary.size() == 4u);
  for (std::size_t i = 0; i < first_boundary.size(); ++i) {
    CHECK(first_boundary[i] == first_edges[i]);
  }

  // The other solid still has its own four boundary edges.
  const std::vector<EdgeId> second_boundary = detail::heal_boundary_edges(brep, second);
  REQUIRE(second_boundary.size() == 4u);
  for (std::size_t i = 0; i < second_boundary.size(); ++i) {
    CHECK(second_boundary[i] == second_edges[i]);
  }
}

TEST_CASE("heal: an unknown solid has no boundary edges") {
  const Brep box = make_box({0, 0, 0}, {1, 1, 1}, 0.0).value();
  CHECK(detail::heal_boundary_edges(box, SolidId{}).empty());
}

TEST_CASE("heal: a stale solid handle has no boundary edges") {
  const Brep box = make_box({0, 0, 0}, {1, 1, 1}, 0.0).value();
  const SolidId solid = box.solid_ids().front();
  const std::uint8_t bumped = static_cast<std::uint8_t>((solid.generation() + 1u) & 0xFFu);
  const SolidId stale = SolidId::make(solid.index(), bumped);
  REQUIRE(!box.solid_alive(stale));
  CHECK(detail::heal_boundary_edges(box, stale).empty());
}

TEST_CASE("heal: issue kinds map to stable names") {
  CHECK(std::string_view(to_string(HealIssueKind::none)) == "none");
  CHECK(std::string_view(to_string(HealIssueKind::gap_too_large)) == "gap_too_large");
  CHECK(std::string_view(to_string(HealIssueKind::ambiguous_gap)) == "ambiguous_gap");
  CHECK(std::string_view(to_string(HealIssueKind::inconsistent_orientation)) ==
        "inconsistent_orientation");
  CHECK(std::string_view(to_string(HealIssueKind::invalid_result)) == "invalid_result");
  CHECK(std::string_view(to_string(static_cast<HealIssueKind>(200))) == "unknown");
}
