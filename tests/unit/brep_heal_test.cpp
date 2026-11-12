#include <array>
#include <cstddef>
#include <string_view>
#include <vector>

#include <doctest/doctest.h>
#include <geom/brep/heal.hpp>
#include <geom/brep/primitives.hpp>
#include <geom/brep/solid_editor.hpp>

using namespace geom;

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

TEST_CASE("heal: an unknown solid has no boundary edges") {
  const Brep box = make_box({0, 0, 0}, {1, 1, 1}, 0.0).value();
  CHECK(detail::heal_boundary_edges(box, SolidId{}).empty());
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
