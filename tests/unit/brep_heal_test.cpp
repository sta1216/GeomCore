#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <iomanip>
#include <sstream>
#include <string>
#include <string_view>
#include <variant>
#include <vector>

#include <doctest/doctest.h>
#include <geom/brep/detail/curved_volume.hpp>
#include <geom/brep/heal.hpp>
#include <geom/brep/primitives.hpp>
#include <geom/brep/solid_editor.hpp>
#include <geom/brep/validate.hpp>

#include "brep_heal_fixture.hpp"

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

// A stable rendering of everything heal_brep must leave alone: pool counts,
// every live entity's stored fields, the vertex coordinates and the geometry
// store. The source Brep is passed by const reference, so comparing this
// before and after a heal pins that the input is never mutated.
std::string input_fingerprint(const Brep& brep) {
  std::ostringstream out;
  out << std::setprecision(17);
  out << "counts " << brep.vertex_count() << ' ' << brep.edge_count() << ' '
      << brep.coedge_count() << ' ' << brep.loop_count() << ' ' << brep.face_count()
      << ' ' << brep.shell_count() << ' ' << brep.solid_count() << '\n';
  for (const VertexId id : brep.vertex_ids()) {
    const Vertex& vertex = brep.vertex(id);
    out << "v " << id.index() << '.' << static_cast<unsigned>(id.generation()) << ' '
        << vertex.point[0] << ' ' << vertex.point[1] << ' ' << vertex.point[2] << ' '
        << vertex.tol << '\n';
  }
  for (const EdgeId id : brep.edge_ids()) {
    const Edge& edge = brep.edge(id);
    out << "e " << id.index() << ' ' << edge.curve << ' ' << edge.t0 << ' ' << edge.t1
        << ' ' << edge.v0.index() << ' ' << edge.v1.index() << ' ' << edge.tol << '\n';
  }
  for (const CoedgeId id : brep.coedge_ids()) {
    const Coedge& coedge = brep.coedge(id);
    out << "c " << id.index() << ' ' << coedge.edge.index() << ' '
        << coedge.face.index() << ' ' << coedge.loop.index() << ' '
        << coedge.next.index() << ' ' << coedge.prev.index() << ' ' << coedge.pcurve
        << ' ' << coedge.reversed << ' ' << coedge.tol << '\n';
  }
  for (const LoopId id : brep.loop_ids()) {
    const Loop& loop = brep.loop(id);
    out << "l " << id.index() << ' ' << loop.face.index() << ' ' << loop.start.index()
        << ' ' << loop.next.index() << ' ' << loop.outer << '\n';
  }
  for (const FaceId id : brep.face_ids()) {
    const Face& face = brep.face(id);
    out << "f " << id.index() << ' ' << face.surface << ' ' << face.first_loop.index()
        << ' ' << face.reversed << ' ' << face.tol << '\n';
  }
  for (const ShellId id : brep.shell_ids()) {
    const Shell& shell = brep.shell(id);
    out << "s " << id.index() << ' ' << shell.first_face << ' ' << shell.face_count
        << ' ' << shell.closed << '\n';
  }
  for (const SolidId id : brep.solid_ids()) {
    const SolidEntity& solid = brep.solid(id);
    out << "d " << id.index() << ' ' << solid.first_shell << ' '
        << solid.shell_count << ' ' << solid.tol << '\n';
  }
  const GeometryStore& geometry = brep.geometry();
  out << "geometry " << geometry.curves.size() << ' ' << geometry.surfaces.size() << ' '
      << geometry.pcurves.size() << '\n';
  for (const SurfaceGeometry& surface : geometry.surfaces) {
    out << "surface " << surface.index();
    if (const auto* plane = std::get_if<Plane3>(&surface)) {
      out << ' ' << plane->normal[0] << ' ' << plane->normal[1] << ' '
          << plane->normal[2] << ' ' << plane->offset;
    }
    out << '\n';
  }
  for (const Pcurve2d& pcurve : geometry.pcurves) {
    out << "pcurve";
    for (const Point2d& point : pcurve.points) {
      out << ' ' << point[0] << ' ' << point[1];
    }
    out << '\n';
  }
  return out.str();
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

TEST_CASE("heal fixture: a torn sheet is a valid open solid") {
  for (const double gap : {1e-6, 1e-3}) {
    const test::TornSheet sheet = test::make_torn_sheet(gap);
    CHECK(sheet.brep.solid_alive(sheet.solid));
    SolidEditor editor(sheet.brep);
    CHECK(editor.commit().has_value());
    CHECK(detail::heal_boundary_edges(sheet.brep, sheet.solid).size() == 8u);
  }
}

TEST_CASE("heal fixture: the auxiliary sheets commit with eight boundary edges") {
  const test::TornSheet split = test::make_torn_split_sheet(1e-6);
  const test::TornSheet partial = test::make_partial_torn_sheet(1e-6, 1e-3);
  const test::TornSheet reversed = test::make_reversed_keep_sheet(1e-6);
  for (const test::TornSheet* sheet : {&split, &partial, &reversed}) {
    SolidEditor editor(sheet->brep);
    CHECK(editor.commit().has_value());
    CHECK(detail::heal_boundary_edges(sheet->brep, sheet->solid).size() == 8u);
  }
}

TEST_CASE("heal: a 1e-6 seam produces four mutual pairs") {
  const test::TornSheet sheet = test::make_torn_sheet(1e-6);
  const std::vector<EdgeId> boundary =
      detail::heal_boundary_edges(sheet.brep, sheet.solid);
  REQUIRE(boundary.size() == 8u);

  const auto plan_result = detail::heal_plan(sheet.brep, boundary, 1e-6);
  REQUIRE(plan_result.has_value());
  const detail::HealPlan& plan = plan_result.value();
  CHECK(plan.issues.empty());
  REQUIRE(plan.pairs.size() == 4u);
  for (const detail::HealEdgePair& pair : plan.pairs) {
    CHECK(pair.keep.index() < pair.drop.index());
    CHECK(pair.gap == doctest::Approx(1e-6).epsilon(1e-6));
    CHECK_FALSE(pair.same_direction);
    for (const detail::HealVertexWeld& weld : pair.welds) {
      CHECK(weld.keep.index() < weld.drop.index());
    }
  }
}

TEST_CASE("heal: a gap past max_gap is measured and reported") {
  const test::TornSheet sheet = test::make_torn_sheet(1e-3);
  const std::vector<EdgeId> boundary =
      detail::heal_boundary_edges(sheet.brep, sheet.solid);
  REQUIRE(boundary.size() == 8u);

  const auto plan_result = detail::heal_plan(sheet.brep, boundary, 1e-6);
  REQUIRE(plan_result.has_value());
  const detail::HealPlan& plan = plan_result.value();
  CHECK(plan.pairs.empty());
  REQUIRE(plan.issues.size() == 1u);
  const HealIssue& issue = plan.issues.front();
  CHECK(issue.kind == HealIssueKind::gap_too_large);
  CHECK(issue.gap == doctest::Approx(1e-3).epsilon(1e-6));
  REQUIRE(issue.entity_a != kNoEntity);
  REQUIRE(issue.entity_b != kNoEntity);
  CHECK(issue.entity_a < issue.entity_b);
}

TEST_CASE("heal: the same traversal direction is not sewn") {
  const test::TornSheet sheet = test::make_torn_sheet(1e-6, true);
  const std::vector<EdgeId> boundary =
      detail::heal_boundary_edges(sheet.brep, sheet.solid);
  REQUIRE(boundary.size() == 8u);

  const auto plan_result = detail::heal_plan(sheet.brep, boundary, 1e-6);
  REQUIRE(plan_result.has_value());
  const detail::HealPlan& plan = plan_result.value();
  CHECK(plan.pairs.empty());
  REQUIRE(plan.issues.size() == 4u);
  for (const HealIssue& issue : plan.issues) {
    CHECK(issue.kind == HealIssueKind::inconsistent_orientation);
    CHECK(issue.entity_a != kNoEntity);
    CHECK(issue.entity_b != kNoEntity);
  }
}

TEST_CASE("heal: ambiguous rims report one issue per middle edge") {
  const test::TornSheet sheet = test::make_ambiguous_rims(1e-6);
  const std::vector<EdgeId> boundary =
      detail::heal_boundary_edges(sheet.brep, sheet.solid);
  REQUIRE(boundary.size() == 12u);

  const auto plan_result = detail::heal_plan(sheet.brep, boundary, 1.5e-6);
  REQUIRE(plan_result.has_value());
  const detail::HealPlan& plan = plan_result.value();
  CHECK(plan.pairs.empty());
  REQUIRE(plan.issues.size() == 4u);
  std::vector<std::uint32_t> ambiguous_edges;
  for (const HealIssue& issue : plan.issues) {
    CHECK(issue.kind == HealIssueKind::ambiguous_gap);
    ambiguous_edges.push_back(issue.entity_a);
  }
  CHECK(ambiguous_edges == std::vector<std::uint32_t>{0, 1, 2, 3});
}

TEST_CASE("heal: a torn sheet becomes a closed shell") {
  const test::TornSheet sheet = test::make_torn_sheet(1e-6);
  const auto healed = heal_brep(sheet.brep, sheet.solid, {});
  REQUIRE(healed.has_value());
  const HealResult& result = healed.value();

  CHECK(result.report.changed);
  CHECK(result.report.boundary_edges == 8u);
  CHECK(result.report.welded_vertices == 4u);
  CHECK(result.report.sewn_edges == 4u);
  CHECK(result.report.issues.empty());
  CHECK(result.report.max_closed_gap == doctest::Approx(1e-6).epsilon(1e-6));
  CHECK(result.report.max_tolerance == doctest::Approx(1e-6).epsilon(1e-6));
  CHECK(detail::heal_boundary_edges(result.brep, sheet.solid).empty());
  CHECK(result.brep.vertex_count() == 4u);
  CHECK(result.brep.edge_count() == 4u);
  CHECK(result.brep.face_count() == 2u);
  const auto shells = result.brep.solid_shells(sheet.solid);
  REQUIRE(shells.size() == 1u);
  CHECK(result.brep.shell(shells.front()).closed);
  CHECK(validate_topology(result.brep).ok());
}

TEST_CASE("heal: a cracked box becomes the unit box") {
  const test::CrackedBox cracked = test::make_cracked_box(1e-6);
  const std::string fingerprint = input_fingerprint(cracked.brep);

  REQUIRE(validate_topology(cracked.brep).ok());
  CHECK(cracked.brep.vertex_count() == 16u);
  CHECK(cracked.brep.edge_count() == 24u);
  CHECK(cracked.brep.coedge_count() == 40u);
  CHECK(cracked.brep.face_count() == 10u);
  const auto cracked_shells = cracked.brep.solid_shells(cracked.solid);
  REQUIRE(cracked_shells.size() == 1u);
  CHECK_FALSE(cracked.brep.shell(cracked_shells.front()).closed);
  CHECK(detail::heal_boundary_edges(cracked.brep, cracked.solid).size() == 8u);
  CHECK(euler_characteristic(cracked.brep, cracked.solid) == 2);
  // The fixture must be committable on its own: the healer is not asked to
  // repair tolerances that were already invalid.
  CHECK(SolidEditor(cracked.brep).commit().has_value());

  // The fixture clamps the upper rim so the seam never exceeds 1e-6; the seam
  // is just below 1e-6 because 0.5 + 1e-6 rounds up and is stepped down one ulp.
  const auto healed =
      heal_brep(cracked.brep, cracked.solid,
                HealOptions{.tolerance = 1e-9, .max_gap = 1e-6, .max_tolerance = 1e-6});
  REQUIRE(healed.has_value());
  const HealResult& result = healed.value();
  CHECK(result.report.changed);
  CHECK(result.report.boundary_edges == 8u);
  CHECK(result.report.welded_vertices == 4u);
  CHECK(result.report.sewn_edges == 4u);
  CHECK(result.report.issues.empty());
  CHECK(result.report.max_closed_gap == doctest::Approx(1e-6).epsilon(1e-6));
  CHECK(result.report.max_tolerance == doctest::Approx(1e-6).epsilon(1e-6));

  CHECK(result.brep.vertex_count() == 12u);
  CHECK(result.brep.edge_count() == 20u);
  CHECK(result.brep.coedge_count() == 40u);
  CHECK(result.brep.face_count() == 10u);
  CHECK(detail::heal_boundary_edges(result.brep, cracked.solid).empty());
  CHECK(euler_characteristic(result.brep, cracked.solid) == 2);
  CHECK(validate_topology(result.brep).ok());
  CHECK(validate_tolerances(result.brep).ok());

  const auto shells = result.brep.solid_shells(cracked.solid);
  REQUIRE(shells.size() == 1u);
  CHECK(result.brep.shell(shells.front()).closed);

  // The rims are sewn at z = 0.5, so the two halves enclose the unit box.
  const auto volume = detail::curved_solid_volume(result.brep, cracked.solid, 1e-9);
  REQUIRE(volume.has_value());
  CHECK(volume.value() == doctest::Approx(1.0).epsilon(1e-9));

  // The healed tolerance chain must satisfy the nesting rules the validator
  // enforces: vertex <= edge <= face <= solid, with the four welded rim
  // vertices lifted to at least the closed seam.
  std::size_t welded_vertices = 0;
  for (const VertexId id : result.brep.vertex_ids()) {
    const Vertex& vertex = result.brep.vertex(id);
    if (vertex.tol > 0.0) {
      ++welded_vertices;
      CHECK(std::abs(vertex.tol - result.report.max_closed_gap) <=
            1e-12 * result.report.max_closed_gap);
    }
  }
  CHECK(welded_vertices == 4u);
  for (const EdgeId id : result.brep.edge_ids()) {
    const Edge& edge = result.brep.edge(id);
    CHECK(edge.tol >= result.brep.vertex(edge.v0).tol);
    CHECK(edge.tol >= result.brep.vertex(edge.v1).tol);
  }
  for (const CoedgeId id : result.brep.coedge_ids()) {
    const Coedge& coedge = result.brep.coedge(id);
    const Face& face = result.brep.face(coedge.face);
    CHECK(face.tol >= result.brep.edge(coedge.edge).tol);
    CHECK(face.tol >= coedge.tol);
  }
  const SolidEntity& solid = result.brep.solid(cracked.solid);
  for (const ShellId shell : shells) {
    for (const FaceId face : result.brep.shell_faces(shell)) {
      CHECK(solid.tol >= result.brep.face(face).tol);
    }
  }

  CHECK(input_fingerprint(cracked.brep) == fingerprint);
}

TEST_CASE("heal: a closed box is returned unchanged") {
  const Brep box = make_box({0, 0, 0}, {1, 1, 1}, 0).value();
  const SolidId solid = box.solid_ids().front();
  const auto healed = heal_brep(box, solid, {});
  REQUIRE(healed.has_value());
  const HealResult& result = healed.value();
  CHECK_FALSE(result.report.changed);
  CHECK(result.report.boundary_edges == 0u);
  CHECK(result.report.welded_vertices == 0u);
  CHECK(result.report.sewn_edges == 0u);
  CHECK(result.report.issues.empty());
  CHECK(result.brep.vertex_count() == box.vertex_count());
  CHECK(result.brep.edge_count() == box.edge_count());
}

TEST_CASE("heal: a seam past max_gap is rejected") {
  const test::CrackedBox cracked = test::make_cracked_box(1.5e-6);
  const std::string fingerprint = input_fingerprint(cracked.brep);

  const auto healed =
      heal_brep(cracked.brep, cracked.solid,
                HealOptions{.tolerance = 1e-9, .max_gap = 1e-6, .max_tolerance = 1e-6});
  REQUIRE(healed.has_value());
  const HealResult& result = healed.value();
  CHECK_FALSE(result.report.changed);
  CHECK(result.report.welded_vertices == 0u);
  CHECK(result.report.sewn_edges == 0u);
  REQUIRE(result.report.issues.size() == 1u);
  CHECK(result.report.issues.front().kind == HealIssueKind::gap_too_large);
  CHECK(result.report.issues.front().gap == doctest::Approx(1.5e-6).epsilon(1e-6));
  CHECK(input_fingerprint(result.brep) == fingerprint);
}

TEST_CASE("heal: a zero-length boundary edge is invalid topology") {
  SolidEditor editor;
  const auto first = editor.add_vertex({0, 0, 0}, 0).value();
  const auto second = editor.add_vertex({0, 0, 0}, 0).value();
  const auto edge = editor.add_edge(first, second, kNoGeometry, 0, 1, 0).value();
  const std::array<EdgeId, 1> boundary{edge};

  const auto plan = detail::heal_plan(editor.brep(), boundary, 1e-6);
  REQUIRE_FALSE(plan.has_value());
  CHECK(plan.error().code == Errc::invalid_topology);
  CHECK(std::string_view(plan.error().message) ==
        "heal_plan: boundary edge has zero length");
}

TEST_CASE("heal: a keep coedge whose raw reversed flag is set still pairs") {
  const test::TornSheet sheet = test::make_reversed_keep_sheet(1e-6);
  const std::vector<EdgeId> boundary =
      detail::heal_boundary_edges(sheet.brep, sheet.solid);
  REQUIRE(boundary.size() == 8u);

  const auto plan_result = detail::heal_plan(sheet.brep, boundary, 1e-6);
  REQUIRE(plan_result.has_value());
  const detail::HealPlan& plan = plan_result.value();
  CHECK(plan.issues.empty());
  REQUIRE(plan.pairs.size() == 4u);
  for (const detail::HealEdgePair& pair : plan.pairs) {
    CoedgeId use{};
    for (const CoedgeId id : sheet.brep.coedge_ids()) {
      if (sheet.brep.coedge(id).edge == pair.keep) {
        use = id;
        break;
      }
    }
    REQUIRE(!use.is_null());
    CHECK(sheet.brep.coedge(use).reversed);
  }

  const auto healed = heal_brep(sheet.brep, sheet.solid, {});
  REQUIRE(healed.has_value());
  CHECK(healed.value().report.sewn_edges == 4u);
  CHECK(detail::heal_boundary_edges(healed.value().brep, sheet.solid).empty());
}

TEST_CASE("heal: a sheet with only non-candidate boundary edges is rejected") {
  SolidEditor editor;
  const auto surface = editor.add_surface(Plane3{{0, 0, 1}, 0});
  const std::array<EdgeId, 4> edges = add_quad_edges(editor, 0.0, 0.0, 1.0);
  const SolidId solid = add_open_sheet(editor, surface, edges);
  const Brep brep = editor.commit().value();

  const auto healed = heal_brep(brep, solid, {});
  REQUIRE(healed.has_value());
  const HealResult& result = healed.value();
  CHECK_FALSE(result.report.changed);
  REQUIRE(result.report.issues.size() == 1u);
  CHECK(result.report.issues.front().kind == HealIssueKind::invalid_result);
  CHECK(std::string_view(result.report.issues.front().message) ==
        "heal_brep: no mutual pair");
}

TEST_CASE("heal: pcurves are repaired on interior edges at welded vertices") {
  const test::TornSheet sheet = test::make_torn_split_sheet(1e-6);
  const auto healed = heal_brep(sheet.brep, sheet.solid, {});
  REQUIRE(healed.has_value());
  const HealResult& result = healed.value();
  CHECK(result.report.changed);
  CHECK(result.report.welded_vertices == 4u);
  CHECK(result.report.sewn_edges == 4u);
  CHECK(detail::heal_boundary_edges(result.brep, sheet.solid).empty());

  const auto coedges = result.brep.coedge_ids();
  REQUIRE(!coedges.empty());
  for (const CoedgeId id : coedges) {
    CHECK(result.brep.coedge(id).pcurve != kNoGeometry);
  }
}

TEST_CASE("heal: a partially sewn sheet is rejected without a partial stitch") {
  const test::TornSheet sheet = test::make_partial_torn_sheet(1e-6, 1e-3);
  const auto healed = heal_brep(sheet.brep, sheet.solid, {});
  REQUIRE(healed.has_value());
  const HealResult& result = healed.value();
  CHECK_FALSE(result.report.changed);
  CHECK(result.report.boundary_edges == 8u);
  CHECK(result.report.welded_vertices == 0u);
  CHECK(result.report.sewn_edges == 0u);
  CHECK(result.report.max_closed_gap == 0.0);
  CHECK(result.report.max_tolerance == 0.0);
  CHECK(result.brep.vertex_count() == sheet.brep.vertex_count());
  CHECK(result.brep.edge_count() == sheet.brep.edge_count());
  REQUIRE(result.report.issues.size() == 2u);
  CHECK(result.report.issues[0].kind == HealIssueKind::invalid_result);
  CHECK(std::string_view(result.report.issues[0].message) ==
        "heal_brep: a shell is still open after stitching");
  CHECK(result.report.issues[1].kind == HealIssueKind::gap_too_large);
}

TEST_CASE("heal: the tolerance cap is checked before welding") {
  const test::TornSheet sheet =
      test::make_torn_sheet_with_tolerance(2e-6, false, 2e-6);
  HealOptions options;
  options.max_gap = 1e-5;
  options.max_tolerance = 1e-6;

  const auto healed = heal_brep(sheet.brep, sheet.solid, options);
  REQUIRE(healed.has_value());
  const HealResult& result = healed.value();
  CHECK_FALSE(result.report.changed);
  CHECK(result.report.welded_vertices == 0u);
  CHECK(result.report.sewn_edges == 0u);
  REQUIRE(result.report.issues.size() == 1u);
  CHECK(result.report.issues.front().kind == HealIssueKind::invalid_result);
  CHECK(std::string_view(result.report.issues.front().message) ==
        "heal_brep: tolerance cap exceeded");
}

TEST_CASE("heal: a one-sided nearest preference is not sewn") {
  const test::AsymmetricPreference fixture =
      test::make_asymmetric_preference_sheet(1e-6);
  const std::array<EdgeId, 5> boundary{fixture.three[0], fixture.three[1],
                                       fixture.three[2], fixture.two[0],
                                       fixture.two[1]};

  const auto plan_result = detail::heal_plan(fixture.brep, boundary, 1.0);
  REQUIRE(plan_result.has_value());
  const detail::HealPlan& plan = plan_result.value();
  CHECK(plan.issues.empty());
  REQUIRE(plan.pairs.size() == 2u);
  for (const detail::HealEdgePair& pair : plan.pairs) {
    CHECK(pair.keep.index() != fixture.three[2].index());
    CHECK(pair.drop.index() != fixture.three[2].index());
  }
}

TEST_CASE("heal: different faces sharing a vertex are not candidates") {
  const test::SharedVertexSheets fixture = test::make_shared_vertex_sheets();
  const std::array<EdgeId, 2> boundary{fixture.first, fixture.second};

  const auto plan_result = detail::heal_plan(fixture.brep, boundary, 1.0);
  REQUIRE(plan_result.has_value());
  CHECK(plan_result.value().pairs.empty());
  CHECK(plan_result.value().issues.empty());
}

TEST_CASE("heal: weld chains collapse to one live keep") {
  SolidEditor editor;
  const auto v0 = editor.add_vertex({0.0, 0.0, 0.0}, 0).value();
  const auto v1 = editor.add_vertex({1.0, 0.0, 0.0}, 0).value();
  const auto v2 = editor.add_vertex({2.0, 0.0, 0.0}, 0).value();
  const auto v3 = editor.add_vertex({3.0, 0.0, 0.0}, 0).value();
  const auto dead = editor.add_vertex({4.0, 0.0, 0.0}, 0).value();
  REQUIRE(editor.merge_vertices(v0, dead).has_value());
  const Brep base = editor.brep();

  std::vector<detail::HealEdgePair> pairs(4);
  pairs[0].welds = {detail::HealVertexWeld{v0, v1},
                    detail::HealVertexWeld{v0, v1}};
  pairs[1].welds = {detail::HealVertexWeld{v1, v2},
                    detail::HealVertexWeld{v1, v2}};
  pairs[2].welds = {detail::HealVertexWeld{v2, v3},
                    detail::HealVertexWeld{v2, v3}};
  pairs[3].welds = {detail::HealVertexWeld{v0, dead},
                    detail::HealVertexWeld{v0, dead}};

  const auto deduped = detail::heal_dedup_welds(base, pairs);
  REQUIRE(deduped.has_value());
  CHECK(deduped.value().size() == 3u);

  SolidEditor apply_editor(base);
  const auto applied = detail::heal_apply_welds(apply_editor, deduped.value());
  REQUIRE(applied.has_value());
  CHECK(applied.value().size() == 3u);
  CHECK(detail::heal_live_weld_count(apply_editor.brep(), applied.value()) == 1u);
  CHECK(apply_editor.brep().vertex_alive(v0));
  CHECK_FALSE(apply_editor.brep().vertex_alive(v1));
  CHECK_FALSE(apply_editor.brep().vertex_alive(v2));
  CHECK_FALSE(apply_editor.brep().vertex_alive(v3));
  CHECK_FALSE(apply_editor.brep().vertex_alive(dead));
}

TEST_CASE("heal: a curved face is rejected before planning") {
  const Brep cylinder =
      make_cylinder({0.0, 0.0, 0.0}, {0.0, 0.0, 1.0}, 1.0, 1.0, 0.0).value();
  const auto healed = heal_brep(cylinder, cylinder.solid_ids().front(), {});
  REQUIRE_FALSE(healed.has_value());
  CHECK(healed.error().code == Errc::unsupported);
  CHECK(std::string_view(healed.error().message) ==
        "heal_brep: curved geometry is not supported yet");
}

TEST_CASE("heal: unknown solid handles are index errors") {
  const test::TornSheet sheet = test::make_torn_sheet(1e-6);
  const SolidId stale = SolidId::make(
      sheet.solid.index(),
      static_cast<std::uint8_t>(sheet.solid.generation() + 1u));
  const SolidId out_of_range = SolidId::make(
      static_cast<std::uint32_t>(sheet.brep.solids().size() + 7u), 0);
  for (const SolidId solid : {SolidId{}, stale, out_of_range}) {
    const auto healed = heal_brep(sheet.brep, solid, {});
    REQUIRE_FALSE(healed.has_value());
    CHECK(healed.error().code == Errc::index_out_of_range);
    CHECK(std::string_view(healed.error().message) ==
          "heal_brep: unknown solid");
  }
}

TEST_CASE("heal: a boundary edge shared with another solid is rejected") {
  SolidEditor editor;
  const auto surface = editor.add_surface(Plane3{{0.0, 0.0, 1.0}, 0.0});
  const auto v0 = editor.add_vertex({0.0, 0.0, 0.0}, 0).value();
  const auto v1 = editor.add_vertex({1.0, 0.0, 0.0}, 0).value();
  const auto v2 = editor.add_vertex({0.0, 1.0, 0.0}, 0).value();
  const auto v3 = editor.add_vertex({-1.0, 0.0, 0.0}, 0).value();

  const auto shared = editor.add_edge(v0, v1, kNoGeometry, 0, 1, 0).value();
  const auto a1 = editor.add_edge(v1, v2, kNoGeometry, 0, 1, 0).value();
  const auto a2 = editor.add_edge(v2, v0, kNoGeometry, 0, 1, 0).value();
  const auto a_face = editor.add_face(surface, false, 0).value();
  const CoedgeSpec a_loop[] = {{shared, false}, {a1, false}, {a2, false}};
  REQUIRE(editor.add_loop(a_face, a_loop, true).has_value());
  const auto a_shell =
      editor.add_shell(std::array<FaceId, 1>{a_face}, false).value();
  const SolidId solid_a =
      editor.add_solid(std::array<ShellId, 1>{a_shell}, 0).value();

  const auto b1 = editor.add_edge(v0, v3, kNoGeometry, 0, 1, 0).value();
  const auto b2 = editor.add_edge(v3, v1, kNoGeometry, 0, 1, 0).value();
  const auto b_face = editor.add_face(surface, false, 0).value();
  const CoedgeSpec b_loop[] = {{shared, true}, {b1, false}, {b2, false}};
  REQUIRE(editor.add_loop(b_face, b_loop, true).has_value());
  const auto b_shell =
      editor.add_shell(std::array<FaceId, 1>{b_face}, false).value();
  REQUIRE(editor.add_solid(std::array<ShellId, 1>{b_shell}, 0).has_value());
  const Brep brep = editor.commit().value();

  const auto healed = heal_brep(brep, solid_a, {});
  REQUIRE_FALSE(healed.has_value());
  CHECK(healed.error().code == Errc::invalid_topology);
  CHECK(std::string_view(healed.error().message) ==
        "heal_brep: boundary edge is shared with another solid");
}
