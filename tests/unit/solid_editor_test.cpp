#include <doctest/doctest.h>

#include <array>
#include <cstdint>
#include <limits>
#include <span>
#include <utility>

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
  REQUIRE(ed.add_vertex({2, 2, 2}, 0.0).has_value());
  CHECK(ed.brep().vertex_count() == 5);
  ed.rollback();
  CHECK(ed.brep().vertex_count() == 4);
  CHECK(ed.brep().face_count() == 1);
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
