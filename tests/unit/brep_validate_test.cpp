#include <doctest/doctest.h>

#include <array>
#include <cstdint>
#include <span>
#include <utility>

#include <geom/brep/solid_editor.hpp>
#include <geom/brep/validate.hpp>
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

struct Box {
  geom::SolidEditor editor;
  geom::SolidId solid;
  geom::ShellId shell;
};

// Closed box built only from topology (no curve geometry): 8 vertices, 12 edges,
// 6 quad faces, one closed shell. Every edge is used once per direction. With
// `split_shells` the six faces are split into two closed three-face shells.
Box make_box(bool split_shells = false) {
  geom::SolidEditor ed;
  const geom::Point3d pts[8] = {{0, 0, 0}, {1, 0, 0}, {1, 1, 0}, {0, 1, 0},
                                {0, 0, 1}, {1, 0, 1}, {1, 1, 1}, {0, 1, 1}};
  std::array<geom::VertexId, 8> vertices{};
  for (int i = 0; i < 8; ++i) vertices[i] = ed.add_vertex(pts[i], 0.0).value();
  const int pairs[12][2] = {{0, 1}, {1, 2}, {2, 3}, {3, 0}, {4, 5}, {5, 6},
                            {6, 7}, {7, 4}, {0, 4}, {1, 5}, {2, 6}, {3, 7}};
  std::array<geom::EdgeId, 12> edges{};
  for (int i = 0; i < 12; ++i) {
    edges[i] = ed.add_edge(vertices[pairs[i][0]], vertices[pairs[i][1]],
                           geom::kNoGeometry, 0.0, 1.0, 0.0)
                   .value();
  }
  const std::uint32_t plane = ed.add_surface(geom::Plane3{{0, 0, 1}, 0.0});
  const int face_edges[6][4] = {{0, 1, 2, 3},  {7, 6, 5, 4},  {8, 4, 9, 0},
                                {9, 5, 10, 1}, {10, 6, 11, 2}, {11, 7, 8, 3}};
  const bool face_reversed[6][4] = {{false, false, false, false},
                                    {true, true, true, true},
                                    {false, false, true, true},
                                    {false, false, true, true},
                                    {false, false, true, true},
                                    {false, false, true, true}};
  std::array<geom::FaceId, 6> faces{};
  for (int i = 0; i < 6; ++i) {
    faces[i] = ed.add_face(plane, false, 0.0).value();
    std::array<geom::CoedgeSpec, 4> specs{};
    for (int j = 0; j < 4; ++j) {
      specs[j] = geom::CoedgeSpec{edges[face_edges[i][j]], face_reversed[i][j]};
    }
    ed.add_loop(faces[i], specs, true).value();
  }
  Box box{std::move(ed), {}, {}};
  if (split_shells) {
    const geom::ShellId first = box.editor
                                    .add_shell(std::array<geom::FaceId, 3>{
                                                   faces[0], faces[1], faces[2]},
                                               true)
                                    .value();
    const geom::ShellId second = box.editor
                                     .add_shell(std::array<geom::FaceId, 3>{
                                                    faces[3], faces[4], faces[5]},
                                                true)
                                     .value();
    box.shell = first;
    box.solid =
        box.editor.add_solid(std::array<geom::ShellId, 2>{first, second}, 0.0).value();
  } else {
    box.shell = box.editor.add_shell(faces, true).value();
    box.solid = box.editor.add_solid(std::array<geom::ShellId, 1>{box.shell}, 0.0).value();
  }
  return box;
}

}  // namespace

TEST_CASE("square face passes topology validation") {
  SquareFace s = make_square_face();
  const auto r = s.editor.commit();
  REQUIRE(r.has_value());
  CHECK(geom::validate_topology(r.value()).ok());
  CHECK(geom::euler_characteristic(r.value(), r.value().solid_ids().front()) == 1);
}

TEST_CASE("empty brep passes topology validation") {
  const geom::SolidEditor ed;
  const auto r = ed.commit();
  REQUIRE(r.has_value());
  CHECK(geom::validate_topology(r.value()).ok());
}

TEST_CASE("closed box passes topology validation and has Euler characteristic two") {
  Box box = make_box();
  const auto r = box.editor.commit();
  REQUIRE(r.has_value());
  CHECK(geom::validate_topology(r.value()).ok());
  CHECK(geom::euler_characteristic(r.value(), box.solid) == 2);
}

TEST_CASE("dangling edge reference is rejected at commit") {
  SquareFace s = make_square_face();
  const geom::EdgeId e0 = s.editor.brep().coedge(s.editor.brep().loop(s.loop).start).edge;
  REQUIRE(s.editor.erase_edge(e0).has_value());
  CHECK(s.editor.brep().edge_count() == 3);
  const auto r = s.editor.commit();
  REQUIRE_FALSE(r.has_value());
  CHECK(r.error().code == geom::Errc::invalid_topology);
  CHECK(r.error().context_a == e0.index());
}

TEST_CASE("face without loops is rejected") {
  geom::SolidEditor ed;
  ed.add_face(geom::kNoGeometry, false, 0.0).value();
  const auto r = ed.commit();
  REQUIRE_FALSE(r.has_value());
  CHECK(r.error().code == geom::Errc::invalid_topology);
}

TEST_CASE("erasing a referenced vertex leaves dangling edges and is rejected") {
  SquareFace s = make_square_face();
  const geom::VertexId v0 = s.vertices[0];
  REQUIRE(s.editor.erase_vertex(v0).has_value());
  CHECK_FALSE(s.editor.brep().vertex_alive(v0));
  const auto r = s.editor.commit();
  REQUIRE_FALSE(r.has_value());
  CHECK(r.error().code == geom::Errc::invalid_topology);
}

TEST_CASE("erase rejects dead and out-of-range ids without changing the editor") {
  SquareFace s = make_square_face();
  const geom::EdgeId live_edge =
      s.editor.brep().coedge(s.editor.brep().loop(s.loop).start).edge;
  const geom::VertexId live_vertex = s.vertices[0];
  const std::uint32_t edges_before = s.editor.brep().edge_count();
  const std::uint32_t vertices_before = s.editor.brep().vertex_count();
  REQUIRE(edges_before == 4);
  REQUIRE(vertices_before == 4);

  const geom::EdgeId dead_edge = geom::EdgeId::make(
      live_edge.index(), static_cast<std::uint8_t>(live_edge.generation() + 1u));
  const auto dead = s.editor.erase_edge(dead_edge);
  CHECK_FALSE(dead.has_value());
  CHECK(dead.error().code == geom::Errc::index_out_of_range);

  const auto out_of_range = s.editor.erase_edge(geom::EdgeId::make(1000u, 1u));
  CHECK_FALSE(out_of_range.has_value());
  CHECK(out_of_range.error().code == geom::Errc::index_out_of_range);

  const geom::VertexId dead_vertex = geom::VertexId::make(
      live_vertex.index(), static_cast<std::uint8_t>(live_vertex.generation() + 1u));
  const auto dead_v = s.editor.erase_vertex(dead_vertex);
  CHECK_FALSE(dead_v.has_value());
  CHECK(dead_v.error().code == geom::Errc::index_out_of_range);

  const auto out_of_range_v = s.editor.erase_vertex(geom::VertexId::make(1000u, 1u));
  CHECK_FALSE(out_of_range_v.has_value());
  CHECK(out_of_range_v.error().code == geom::Errc::index_out_of_range);

  CHECK(s.editor.brep().edge_count() == edges_before);
  CHECK(s.editor.brep().vertex_count() == vertices_before);
  CHECK(s.editor.brep().edge_alive(live_edge));
  CHECK(s.editor.brep().vertex_alive(live_vertex));
  CHECK(s.editor.commit().has_value());
}

TEST_CASE("edge referenced by three coedges is rejected") {
  SquareFace s = make_square_face();
  const geom::EdgeId e0 = s.editor.brep().coedge(s.editor.brep().loop(s.loop).start).edge;
  const std::array<geom::CoedgeSpec, 3> specs{
      geom::CoedgeSpec{e0, false}, geom::CoedgeSpec{e0, true},
      geom::CoedgeSpec{e0, false}};
  const auto r = s.editor.add_loop(s.face, specs, false);
  CHECK_FALSE(r.has_value());
  CHECK(r.error().code == geom::Errc::invalid_topology);
}

TEST_CASE("edge used twice in the same direction is rejected") {
  SquareFace s = make_square_face();
  const geom::Loop loop = s.editor.brep().loop(s.loop);
  std::array<geom::CoedgeSpec, 4> specs{};
  geom::CoedgeId current = loop.start;
  for (int i = 0; i < 4; ++i) {
    const geom::Coedge& coedge = s.editor.brep().coedge(current);
    specs[i] = geom::CoedgeSpec{coedge.edge, false};
    current = coedge.next;
  }
  const geom::FaceId second = s.editor.add_face(geom::kNoGeometry, false, 0.0).value();
  REQUIRE(s.editor.add_loop(second, specs, true).has_value());
  const auto r = s.editor.commit();
  REQUIRE_FALSE(r.has_value());
  CHECK(r.error().code == geom::Errc::invalid_topology);
}

TEST_CASE("edge that no coedge uses is rejected") {
  SquareFace s = make_square_face();
  const geom::VertexId v = s.editor.add_vertex({2, 2, 2}, 0.0).value();
  const geom::VertexId w = s.editor.add_vertex({3, 3, 3}, 0.0).value();
  REQUIRE(s.editor.add_edge(v, w, geom::kNoGeometry, 0.0, 1.0, 0.0).has_value());
  const auto r = s.editor.commit();
  REQUIRE_FALSE(r.has_value());
  CHECK(r.error().code == geom::Errc::invalid_topology);
  CHECK(r.error().context_a == 4);  // the unused edge's index; no coedge shares it
}

TEST_CASE("closed shell with a boundary edge reports not_closed") {
  geom::SolidEditor ed;
  const geom::Point3d pts[4] = {{0, 0, 0}, {1, 0, 0}, {1, 1, 0}, {0, 1, 0}};
  std::array<geom::VertexId, 4> vertices{};
  for (int i = 0; i < 4; ++i) vertices[i] = ed.add_vertex(pts[i], 0.0).value();
  std::array<geom::EdgeId, 4> edges{};
  for (int i = 0; i < 4; ++i) {
    edges[i] = ed.add_edge(vertices[i], vertices[(i + 1) % 4], geom::kNoGeometry, 0.0, 1.0,
                           0.0)
                   .value();
  }
  const geom::FaceId face = ed.add_face(geom::kNoGeometry, false, 0.0).value();
  std::array<geom::CoedgeSpec, 4> specs{};
  for (int i = 0; i < 4; ++i) specs[i] = geom::CoedgeSpec{edges[i], false};
  ed.add_loop(face, specs, true).value();
  const geom::ShellId shell = ed.add_shell(std::array<geom::FaceId, 1>{face}, true).value();
  ed.add_solid(std::array<geom::ShellId, 1>{shell}, 0.0).value();
  const auto r = ed.commit();
  REQUIRE_FALSE(r.has_value());
  CHECK(r.error().code == geom::Errc::not_closed);
}

TEST_CASE("two open sheets in closed shells are not closed shells") {
  SquareFace s = make_square_face();
  const geom::Loop loop = s.editor.brep().loop(s.loop);
  std::array<geom::CoedgeSpec, 4> specs{};
  geom::CoedgeId current = loop.start;
  for (int i = 0; i < 4; ++i) {
    const geom::Coedge& coedge = s.editor.brep().coedge(current);
    specs[i] = geom::CoedgeSpec{coedge.edge, true};
    current = coedge.next;
  }
  const geom::FaceId back = s.editor.add_face(geom::kNoGeometry, false, 0.0).value();
  REQUIRE(s.editor.add_loop(back, specs, true).has_value());
  const geom::ShellId front_shell =
      s.editor.add_shell(std::array<geom::FaceId, 1>{s.face}, true).value();
  const geom::ShellId back_shell =
      s.editor.add_shell(std::array<geom::FaceId, 1>{back}, true).value();
  REQUIRE(s.editor.add_solid(std::array<geom::ShellId, 2>{front_shell, back_shell}, 0.0)
              .has_value());
  const auto r = s.editor.commit();
  REQUIRE_FALSE(r.has_value());
  CHECK(r.error().code == geom::Errc::not_closed);
}

TEST_CASE("box split into two closed three-face shells is not closed") {
  Box box = make_box(true);
  const auto r = box.editor.commit();
  REQUIRE_FALSE(r.has_value());
  CHECK(r.error().code == geom::Errc::not_closed);
}

TEST_CASE("shell listing a face twice is rejected") {
  SquareFace s = make_square_face();
  REQUIRE(s.editor.add_shell(std::array<geom::FaceId, 2>{s.face, s.face}, true).has_value());
  const auto r = s.editor.commit();
  REQUIRE_FALSE(r.has_value());
  CHECK(r.error().code == geom::Errc::invalid_topology);
}

TEST_CASE("solid listing a shell twice is rejected") {
  Box box = make_box();
  REQUIRE(box.editor.add_solid(std::array<geom::ShellId, 2>{box.shell, box.shell}, 0.0)
              .has_value());
  const auto r = box.editor.commit();
  REQUIRE_FALSE(r.has_value());
  CHECK(r.error().code == geom::Errc::invalid_topology);
}
