#include <doctest/doctest.h>

#include <array>
#include <cstdint>
#include <limits>
#include <span>
#include <string_view>
#include <utility>

#include <geom/brep/solid_editor.hpp>
#include <geom/core/error.hpp>

#include "m4b_heal_fixture.hpp"

namespace {

// Walks a loop's next-cycle and returns its coedge count. The guard keeps a
// corrupt chain from turning a failed assertion into an infinite loop.
std::uint32_t loop_coedge_count(const geom::Brep& brep, geom::LoopId loop) {
  const std::uint32_t guard = brep.coedge_count() + 1u;
  std::uint32_t count = 0;
  geom::CoedgeId current = brep.loop(loop).start;
  do {
    current = brep.coedge(current).next;
    ++count;
  } while (current != brep.loop(loop).start && count < guard);
  return count;
}

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

}  // namespace

TEST_CASE("editor commits a single-face open shell") {
  SquareFace s = make_square_face();
  const auto r = s.editor.commit();
  REQUIRE(r.has_value());
  const geom::Brep& b = r.value();
  CHECK(b.vertex_count() == 4);
  CHECK(b.edge_count() == 4);
  CHECK(b.coedge_count() == 4);
  CHECK(b.face_count() == 1);
  CHECK(b.solid_count() == 1);
  CHECK(b.shell(b.solid_shells(b.solid_ids().front()).front()).closed == false);
}

TEST_CASE("loop coedges form a closed cycle") {
  SquareFace s = make_square_face();
  const auto r = s.editor.commit();
  REQUIRE(r.has_value());
  const geom::Brep& b = r.value();
  const geom::Loop& loop = b.loop(s.loop);
  geom::CoedgeId id = loop.start;
  int steps = 0;
  do {
    const geom::Coedge& c = b.coedge(id);
    CHECK(c.loop == s.loop);
    CHECK(c.face == s.face);
    id = c.next;
    ++steps;
  } while (id != loop.start && steps <= 4);
  CHECK(steps == 4);
}

TEST_CASE("commit rejects a loop with no coedges") {
  geom::SolidEditor ed;
  const std::uint32_t plane = ed.add_surface(geom::Plane3{{0, 0, 1}, 0.0});
  const geom::FaceId f = ed.add_face(plane, false, 0.0).value();
  const auto r = ed.add_loop(f, std::span<const geom::CoedgeSpec>{}, true);
  CHECK_FALSE(r.has_value());
  CHECK(r.error().code == geom::Errc::invalid_topology);
  CHECK(ed.brep().face_count() == 1);  // 失败的 add_loop 没有破坏已有实体
  CHECK(ed.brep().loop_count() == 0);
}

TEST_CASE("rollback discards staged entities") {
  geom::SolidEditor ed;
  CHECK(ed.add_vertex({1, 2, 3}, 0.0).has_value());
  CHECK(ed.brep().vertex_count() == 1);
  ed.rollback();
  CHECK(ed.brep().vertex_count() == 0);
}

TEST_CASE("rollback restores the construction base") {
  SquareFace source = make_square_face();
  const auto committed = source.editor.commit();
  REQUIRE(committed.has_value());
  geom::SolidEditor ed(committed.value());
  const std::size_t base_curves = ed.brep().geometry().curves.size();
  const std::size_t base_surfaces = ed.brep().geometry().surfaces.size();
  const std::size_t base_pcurves = ed.brep().geometry().pcurves.size();
  REQUIRE(base_curves == 1u);
  REQUIRE(base_surfaces == 1u);
  REQUIRE(base_pcurves == 0u);
  REQUIRE(ed.add_vertex({2, 2, 2}, 0.0).has_value());
  const std::uint32_t curve = ed.add_curve(geom::LineCurve3{{0, 0, 0}, {0, 1, 0}});
  const std::uint32_t surface = ed.add_surface(geom::Plane3{{1, 0, 0}, 0.0});
  const std::uint32_t pcurve = ed.add_pcurve(geom::Pcurve2d{{0, 0}, {1, 1}});
  CHECK(curve == static_cast<std::uint32_t>(base_curves));
  CHECK(surface == static_cast<std::uint32_t>(base_surfaces));
  CHECK(pcurve == static_cast<std::uint32_t>(base_pcurves));
  CHECK(ed.brep().vertex_count() == 5);
  CHECK(ed.brep().geometry().curves.size() == base_curves + 1u);
  CHECK(ed.brep().geometry().surfaces.size() == base_surfaces + 1u);
  CHECK(ed.brep().geometry().pcurves.size() == base_pcurves + 1u);
  ed.rollback();
  CHECK(ed.brep().vertex_count() == 4);
  CHECK(ed.brep().face_count() == 1);
  CHECK(ed.brep().geometry().curves.size() == base_curves);
  CHECK(ed.brep().geometry().surfaces.size() == base_surfaces);
  CHECK(ed.brep().geometry().pcurves.size() == base_pcurves);
}

TEST_CASE("editor rejects non-finite vertex positions") {
  geom::SolidEditor ed;
  const auto r = ed.add_vertex({std::numeric_limits<double>::infinity(), 0.0, 0.0}, 0.0);
  CHECK_FALSE(r.has_value());
  CHECK(r.error().code == geom::Errc::non_finite);
  CHECK(ed.brep().vertex_count() == 0);
}

TEST_CASE("editor rejects a negative vertex tolerance") {
  geom::SolidEditor ed;
  const auto r = ed.add_vertex({0, 0, 0}, -1e-6);
  CHECK_FALSE(r.has_value());
  CHECK(r.error().code == geom::Errc::tolerance_violation);
  CHECK(ed.brep().vertex_count() == 0);
}

TEST_CASE("editor rejects an edge with a dead vertex") {
  geom::SolidEditor ed;
  const geom::VertexId v0 = ed.add_vertex({0, 0, 0}, 0.0).value();
  const geom::VertexId stale = geom::VertexId::make(0, 2);
  const auto r = ed.add_edge(v0, stale, geom::kNoGeometry, 0.0, 1.0, 0.0);
  CHECK_FALSE(r.has_value());
  CHECK(r.error().code == geom::Errc::index_out_of_range);
  CHECK(ed.brep().edge_count() == 0);
}

TEST_CASE("editor rejects an empty edge parameter interval") {
  geom::SolidEditor ed;
  const geom::VertexId v0 = ed.add_vertex({0, 0, 0}, 0.0).value();
  const geom::VertexId v1 = ed.add_vertex({1, 0, 0}, 0.0).value();
  const auto r = ed.add_edge(v0, v1, geom::kNoGeometry, 1.0, 1.0, 0.0);
  CHECK_FALSE(r.has_value());
  CHECK(r.error().code == geom::Errc::degenerate);
  CHECK(ed.brep().edge_count() == 0);
}

TEST_CASE("editor rejects an out-of-range curve index") {
  geom::SolidEditor ed;
  const geom::VertexId v0 = ed.add_vertex({0, 0, 0}, 0.0).value();
  const geom::VertexId v1 = ed.add_vertex({1, 0, 0}, 0.0).value();
  const auto r = ed.add_edge(v0, v1, 0u, 0.0, 1.0, 0.0);
  CHECK_FALSE(r.has_value());
  CHECK(r.error().code == geom::Errc::index_out_of_range);
  CHECK(ed.brep().edge_count() == 0);
}

TEST_CASE("editor rejects an out-of-range surface index") {
  geom::SolidEditor ed;
  const auto r = ed.add_face(0u, false, 0.0);
  CHECK_FALSE(r.has_value());
  CHECK(r.error().code == geom::Errc::index_out_of_range);
  CHECK(ed.brep().face_count() == 0);
}

TEST_CASE("editor rejects a shell with no faces") {
  geom::SolidEditor ed;
  const auto r = ed.add_shell(std::span<const geom::FaceId>{}, false);
  CHECK_FALSE(r.has_value());
  CHECK(r.error().code == geom::Errc::invalid_topology);
  CHECK(ed.brep().shell_count() == 0);
}

TEST_CASE("editor rejects a solid with no shells") {
  geom::SolidEditor ed;
  const auto r = ed.add_solid(std::span<const geom::ShellId>{}, 0.0);
  CHECK_FALSE(r.has_value());
  CHECK(r.error().code == geom::Errc::invalid_topology);
  CHECK(ed.brep().solid_count() == 0);
}

TEST_CASE("editor accepts a coedge pcurve that exists") {
  geom::SolidEditor ed;
  const geom::VertexId v0 = ed.add_vertex({0, 0, 0}, 0.0).value();
  const geom::VertexId v1 = ed.add_vertex({1, 0, 0}, 0.0).value();
  const geom::EdgeId edge = ed.add_edge(v0, v1, geom::kNoGeometry, 0.0, 1.0, 0.0).value();
  const std::uint32_t plane = ed.add_surface(geom::Plane3{{0, 0, 1}, 0.0});
  const geom::FaceId face = ed.add_face(plane, false, 0.0).value();
  const std::uint32_t pcurve = ed.add_pcurve(geom::Pcurve2d{{0, 0}, {1, 0}});
  CHECK(pcurve == 0u);
  const std::array<geom::CoedgeSpec, 1> specs{geom::CoedgeSpec{edge, false, pcurve, 0.0}};
  const auto loop = ed.add_loop(face, specs, true);
  REQUIRE(loop.has_value());
  const auto committed = ed.commit();
  REQUIRE(committed.has_value());
  const geom::Brep& b = committed.value();
  CHECK(b.geometry().pcurves.size() == 1u);
  CHECK(b.coedge(b.loop(loop.value()).start).pcurve == pcurve);
}

TEST_CASE("editor rejects an out-of-range pcurve index") {
  geom::SolidEditor ed;
  const geom::VertexId v0 = ed.add_vertex({0, 0, 0}, 0.0).value();
  const geom::VertexId v1 = ed.add_vertex({1, 0, 0}, 0.0).value();
  const geom::EdgeId edge = ed.add_edge(v0, v1, geom::kNoGeometry, 0.0, 1.0, 0.0).value();
  const std::uint32_t plane = ed.add_surface(geom::Plane3{{0, 0, 1}, 0.0});
  const geom::FaceId face = ed.add_face(plane, false, 0.0).value();
  const std::array<geom::CoedgeSpec, 1> specs{geom::CoedgeSpec{edge, false, 0u, 0.0}};
  const auto r = ed.add_loop(face, specs, true);
  CHECK_FALSE(r.has_value());
  CHECK(r.error().code == geom::Errc::index_out_of_range);
  CHECK(ed.brep().loop_count() == 0);
  CHECK(ed.brep().coedge_count() == 0);
}

TEST_CASE("editor rejects shells and solids with dead members") {
  SquareFace s = make_square_face();
  geom::SolidEditor& ed = s.editor;
  // The populated editor owns face 0 and shell 0, so make(0, 2) points at a live
  // slot with the wrong generation: the rejection must come from the generation
  // check, not from an index-versus-size check.
  REQUIRE(ed.brep().face_alive(geom::FaceId::make(0, 1)));
  REQUIRE(ed.brep().shell_alive(geom::ShellId::make(0, 1)));
  const std::size_t base_shells = ed.brep().shell_count();
  const std::size_t base_solids = ed.brep().solid_count();
  REQUIRE(base_shells == 1u);
  REQUIRE(base_solids == 1u);

  const auto shell = ed.add_shell(std::array<geom::FaceId, 1>{geom::FaceId::make(0, 2)}, false);
  CHECK_FALSE(shell.has_value());
  CHECK(shell.error().code == geom::Errc::index_out_of_range);
  CHECK(ed.brep().shell_count() == base_shells);

  const auto solid = ed.add_solid(std::array<geom::ShellId, 1>{geom::ShellId::make(0, 2)}, 0.0);
  CHECK_FALSE(solid.has_value());
  CHECK(solid.error().code == geom::Errc::index_out_of_range);
  CHECK(ed.brep().solid_count() == base_solids);
}

TEST_CASE("empty editor commits an empty brep") {
  const geom::SolidEditor ed;
  const auto r = ed.commit();
  REQUIRE(r.has_value());
  CHECK(r.value().vertex_count() == 0u);
  CHECK(r.value().edge_count() == 0u);
  CHECK(r.value().loop_count() == 0u);
  CHECK(r.value().solid_count() == 0u);
  CHECK(r.value().geometry().curves.empty());
}

TEST_CASE("commit rejects a face with two outer loops") {
  SquareFace s = make_square_face();
  const geom::Point3d pts[4] = {{0, 0, 0}, {1, 0, 0}, {1, 1, 0}, {0, 1, 0}};
  std::array<geom::VertexId, 4> vertices{};
  for (int i = 0; i < 4; ++i) vertices[i] = s.editor.add_vertex(pts[i], 0.0).value();
  std::array<geom::EdgeId, 4> edges{};
  for (int i = 0; i < 4; ++i) {
    edges[i] = s.editor
                   .add_edge(vertices[i], vertices[(i + 1) % 4], geom::kNoGeometry, 0.0, 1.0,
                             0.0)
                   .value();
  }
  std::array<geom::CoedgeSpec, 4> specs{};
  for (int i = 0; i < 4; ++i) specs[i] = geom::CoedgeSpec{edges[i], false};
  REQUIRE(s.editor.add_loop(s.face, specs, true).has_value());
  const auto r = s.editor.commit();
  CHECK_FALSE(r.has_value());
  CHECK(r.error().code == geom::Errc::invalid_topology);
  CHECK(s.editor.brep().loop_count() == 2);
}

TEST_CASE("solid editor: dissolve_edge of a pinhole rim arc keeps the loop") {
  geom::SolidEditor ed;
  const geom::test::m4b::ConePinhole cone =
      geom::test::m4b::add_cone_pinhole(ed, 0.0, 0.25, true);
  REQUIRE(ed.brep().edge_count() == 5u);
  REQUIRE(ed.brep().coedge_count() == 8u);
  REQUIRE(loop_coedge_count(ed.brep(), cone.lateral_loop) == 6u);
  REQUIRE(loop_coedge_count(ed.brep(), cone.base_loop) == 2u);
  const geom::CoedgeId old_start = ed.brep().loop(cone.lateral_loop).start;
  REQUIRE(ed.brep().coedge(old_start).edge == cone.base0);

  const auto r = ed.dissolve_edge(cone.top_a);
  REQUIRE(r.has_value());

  CHECK_FALSE(ed.brep().edge_alive(cone.top_a));
  CHECK(ed.brep().vertex_alive(cone.t0));
  CHECK(ed.brep().vertex_alive(cone.t1));
  CHECK(ed.brep().edge_alive(cone.top_b));
  CHECK(ed.brep().edge_alive(cone.seam));
  CHECK(loop_coedge_count(ed.brep(), cone.lateral_loop) == 5u);
  CHECK(loop_coedge_count(ed.brep(), cone.base_loop) == 2u);
  CHECK(ed.brep().loop(cone.lateral_loop).start == old_start);
  CHECK(ed.brep().edge_count() == 4u);
  CHECK(ed.brep().coedge_count() == 7u);
  CHECK(ed.brep().loop_count() == 2u);
  CHECK(ed.brep().shell_count() == 1u);
  CHECK(ed.brep().solid_count() == 1u);
  CHECK(ed.commit().has_value());
}

TEST_CASE("solid editor: dissolve_edge of a twice-used edge relinks one loop") {
  geom::SolidEditor ed;
  const geom::test::m4b::ConePinch cone = geom::test::m4b::add_cone_pinch(ed, 0.0);
  REQUIRE(ed.brep().edge_count() == 3u);
  REQUIRE(ed.brep().coedge_count() == 6u);
  REQUIRE(loop_coedge_count(ed.brep(), cone.lateral_loop) == 4u);
  REQUIRE(loop_coedge_count(ed.brep(), cone.base_loop) == 2u);

  const auto r = ed.dissolve_edge(cone.seam);
  REQUIRE(r.has_value());

  CHECK_FALSE(ed.brep().edge_alive(cone.seam));
  CHECK(ed.brep().vertex_alive(cone.apex));
  CHECK(loop_coedge_count(ed.brep(), cone.lateral_loop) == 2u);
  CHECK(loop_coedge_count(ed.brep(), cone.base_loop) == 2u);
  CHECK(ed.brep().edge_count() == 2u);
  CHECK(ed.brep().coedge_count() == 4u);
  CHECK(ed.brep().loop_count() == 2u);
  CHECK(ed.brep().shell_count() == 1u);
  CHECK(ed.brep().solid_count() == 1u);
  const auto committed = ed.commit();
  REQUIRE(committed.has_value());
  CHECK(geom::validate_topology(committed.value()).ok());
}

TEST_CASE("solid editor: dissolve_edge shared by two loops relocates both starts") {
  geom::SolidEditor ed;
  const geom::test::m4b::ConePinhole cone =
      geom::test::m4b::add_cone_pinhole(ed, 0.0, 0.25, true);
  const geom::CoedgeId lateral_start = ed.brep().loop(cone.lateral_loop).start;
  const geom::CoedgeId base_start = ed.brep().loop(cone.base_loop).start;
  REQUIRE(ed.brep().coedge(lateral_start).edge == cone.base0);
  REQUIRE(ed.brep().coedge(base_start).edge == cone.base0);

  const auto r = ed.dissolve_edge(cone.base0);
  REQUIRE(r.has_value());

  CHECK_FALSE(ed.brep().edge_alive(cone.base0));
  CHECK(ed.brep().vertex_alive(cone.b0));
  CHECK(ed.brep().vertex_alive(cone.b1));
  CHECK(loop_coedge_count(ed.brep(), cone.lateral_loop) == 5u);
  CHECK(loop_coedge_count(ed.brep(), cone.base_loop) == 1u);
  const geom::CoedgeId new_lateral_start = ed.brep().loop(cone.lateral_loop).start;
  const geom::CoedgeId new_base_start = ed.brep().loop(cone.base_loop).start;
  CHECK(new_lateral_start != lateral_start);
  CHECK(new_base_start != base_start);
  CHECK(ed.brep().coedge(new_lateral_start).edge == cone.base1);
  CHECK(ed.brep().coedge(new_base_start).edge == cone.base1);
  CHECK(ed.brep().edge_count() == 4u);
  CHECK(ed.brep().coedge_count() == 6u);
  const auto committed = ed.commit();
  REQUIRE(committed.has_value());
  CHECK(geom::validate_topology(committed.value()).ok());
}

TEST_CASE("solid editor: dissolve_edge refuses to empty a loop") {
  geom::SolidEditor ed;
  const geom::test::m4b::ConePinch cone = geom::test::m4b::add_cone_pinch(ed, 0.0);
  REQUIRE(ed.dissolve_edge(cone.base0).has_value());
  REQUIRE(loop_coedge_count(ed.brep(), cone.base_loop) == 1u);
  REQUIRE(loop_coedge_count(ed.brep(), cone.lateral_loop) == 3u);

  const auto r = ed.dissolve_edge(cone.base1);
  REQUIRE_FALSE(r.has_value());
  CHECK(r.error().code == geom::Errc::unsupported);
  CHECK(std::string_view(r.error().message) == "solid_editor: dissolve would empty a loop");
  CHECK(r.error().context_a == cone.base1.index());
  // The refused call staged nothing.
  CHECK(ed.brep().edge_alive(cone.base1));
  CHECK(loop_coedge_count(ed.brep(), cone.base_loop) == 1u);
  CHECK(loop_coedge_count(ed.brep(), cone.lateral_loop) == 3u);
  CHECK(ed.brep().edge_count() == 2u);
  CHECK(ed.brep().coedge_count() == 4u);
  CHECK(ed.commit().has_value());
}

TEST_CASE("solid editor: dissolve_edge rejects a dead edge") {
  geom::SolidEditor ed;
  const geom::test::m4b::ConePinch cone = geom::test::m4b::add_cone_pinch(ed, 0.0);
  REQUIRE(ed.dissolve_edge(cone.seam).has_value());
  REQUIRE(ed.brep().coedge_count() == 4u);
  REQUIRE(ed.brep().edge_count() == 2u);

  const auto r = ed.dissolve_edge(cone.seam);
  REQUIRE_FALSE(r.has_value());
  CHECK(r.error().code == geom::Errc::index_out_of_range);
  CHECK(std::string_view(r.error().message) == "solid_editor: edge is not alive");
  CHECK(r.error().context_a == cone.seam.index());
  CHECK(ed.brep().coedge_count() == 4u);
  CHECK(ed.brep().edge_count() == 2u);
}

TEST_CASE("solid editor: dissolve_edge keeps a closed manifold valid") {
  geom::SolidEditor ed;
  const geom::test::m4b::ConePinch cone = geom::test::m4b::add_cone_pinch(ed, 0.0);
  const auto closed = ed.commit();
  REQUIRE(closed.has_value());
  CHECK(geom::validate_topology(closed.value()).ok());
  REQUIRE(geom::euler_characteristic(closed.value(), cone.solid) == 2);

  REQUIRE(ed.dissolve_edge(cone.seam).has_value());
  const auto committed = ed.commit();
  REQUIRE(committed.has_value());
  const geom::Brep& brep = committed.value();
  CHECK(geom::validate_topology(brep).ok());
  CHECK(geom::euler_characteristic(brep, cone.solid) == 2);
  CHECK(brep.shell(cone.shell).closed);
  CHECK(brep.vertex_alive(cone.apex));
  CHECK(brep.edge_count() == 2u);
  CHECK(brep.coedge_count() == 4u);
  CHECK(brep.face_count() == 2u);
  CHECK(brep.loop_count() == 2u);
  CHECK(brep.solid_count() == 1u);
}
