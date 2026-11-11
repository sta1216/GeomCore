#include <doctest/doctest.h>

#include <array>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <utility>
#include <variant>

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
  const std::uint32_t line = s.editor.add_curve(geom::LineCurve3{{0, 0, 0}, {1, 0, 0}});
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

// One edge, one closed single-coedge loop on one face: the smallest editor state that can carry a
// pcurve and be committed.
struct SingleEdgeFace {
  geom::SolidEditor editor;
  geom::EdgeId edge;
  geom::FaceId face;
  std::array<geom::VertexId, 2> vertices;
};

SingleEdgeFace make_single_edge_face(double tol) {
  geom::SolidEditor ed;
  SingleEdgeFace s{std::move(ed), {}, {}, {}};
  s.vertices[0] = s.editor.add_vertex({0, 0, 0}, 0.0).value();
  s.vertices[1] = s.editor.add_vertex({1, 0, 0}, 0.0).value();
  const std::uint32_t line = s.editor.add_curve(geom::LineCurve3{{0, 0, 0}, {1, 0, 0}});
  s.edge = s.editor.add_edge(s.vertices[0], s.vertices[1], line, 0.0, 1.0, tol).value();
  const std::uint32_t plane = s.editor.add_surface(geom::Plane3{{0, 0, 1}, 0.0});
  s.face = s.editor.add_face(plane, false, tol).value();
  return s;
}

// Entity and geometry-store sizes, used to prove rejected operations do not
// mutate staging.
struct EditorCounts {
  std::uint32_t vertices = 0;
  std::uint32_t edges = 0;
  std::uint32_t coedges = 0;
  std::size_t curves = 0;
  std::size_t surfaces = 0;
  std::size_t pcurves = 0;
};

EditorCounts snapshot(const geom::SolidEditor& editor) {
  const geom::Brep& b = editor.brep();
  return EditorCounts{b.vertex_count(), b.edge_count(), b.coedge_count(),
                      b.geometry().curves.size(), b.geometry().surfaces.size(),
                      b.geometry().pcurves.size()};
}

void check_unchanged(const geom::SolidEditor& editor, const EditorCounts& expected) {
  const EditorCounts actual = snapshot(editor);
  CHECK(actual.vertices == expected.vertices);
  CHECK(actual.edges == expected.edges);
  CHECK(actual.coedges == expected.coedges);
  CHECK(actual.curves == expected.curves);
  CHECK(actual.surfaces == expected.surfaces);
  CHECK(actual.pcurves == expected.pcurves);
}

}  // namespace

TEST_CASE("split edge inserts a vertex and keeps the loop closed") {
  SquareFace s = make_square_face();
  const geom::Brep before = s.editor.commit().value();
  const geom::EdgeId e0 = before.coedge(before.loop(s.loop).start).edge;
  REQUIRE(s.editor.split_edge(e0, 0.5).has_value());
  const auto after = s.editor.commit();
  REQUIRE(after.has_value());
  const geom::Brep& b = after.value();
  CHECK(b.vertex_count() == before.vertex_count() + 1);
  CHECK(b.edge_count() == before.edge_count() + 1);
  CHECK(b.coedge_count() == before.coedge_count() + 1);
  CHECK(geom::validate_topology(b).ok());
  CHECK(geom::euler_characteristic(b, b.solid_ids().front()) == 1);
}

TEST_CASE("join vertex restores the original counts") {
  SquareFace s = make_square_face();
  const geom::EdgeId e0 = s.editor.brep().coedge(s.editor.brep().loop(s.loop).start).edge;
  const geom::VertexId mid = s.editor.split_edge(e0, 0.5).value();
  CHECK(s.editor.join_vertex(mid).has_value());
  const auto r = s.editor.commit();
  REQUIRE(r.has_value());
  CHECK(r.value().vertex_count() == 4);
  CHECK(r.value().edge_count() == 4);
  CHECK(geom::validate_topology(r.value()).ok());
}

TEST_CASE("split at an endpoint is rejected") {
  SquareFace s = make_square_face();
  const geom::EdgeId e0 = s.editor.brep().coedge(s.editor.brep().loop(s.loop).start).edge;
  CHECK_FALSE(s.editor.split_edge(e0, 0.0).has_value());
  CHECK_FALSE(s.editor.split_edge(e0, 1.0).has_value());
  CHECK(s.editor.split_edge(e0, 0.25).has_value());
}

TEST_CASE("join rejects a vertex with three incident edges") {
  SquareFace s = make_square_face();
  const geom::VertexId extra =
      s.editor.add_vertex({0.5, -1.0, 0.0}, 0.0).value();
  const std::uint32_t line =
      s.editor.add_curve(geom::LineCurve3{{1, 0, 0}, {0, 1, 0}});
  REQUIRE(s.editor.add_edge(s.vertices[1], extra, line, 0.0, 1.0, 0.0).has_value());
  const auto joined = s.editor.join_vertex(s.vertices[1]);
  REQUIRE_FALSE(joined.has_value());
  CHECK(joined.error().code == geom::Errc::invalid_topology);
}

TEST_CASE("split edge duplicates every coedge use") {
  geom::SolidEditor ed;
  const geom::VertexId v0 = ed.add_vertex({0, 0, 0}, 0.0).value();
  const geom::VertexId v1 = ed.add_vertex({1, 0, 0}, 0.0).value();
  const std::uint32_t line = ed.add_curve(geom::LineCurve3{{0, 0, 0}, {1, 0, 0}});
  const geom::EdgeId edge = ed.add_edge(v0, v1, line, 0.0, 1.0, 0.0).value();
  const std::uint32_t plane = ed.add_surface(geom::Plane3{{0, 0, 1}, 0.0});
  const geom::FaceId face = ed.add_face(plane, false, 0.0).value();
  // A lens loop visits the same edge once forward and once reversed.
  const std::array<geom::CoedgeSpec, 2> specs{geom::CoedgeSpec{edge, false},
                                              geom::CoedgeSpec{edge, true}};
  const geom::LoopId loop = ed.add_loop(face, specs, true).value();
  REQUIRE(ed.split_edge(edge, 0.5).has_value());
  const auto r = ed.commit();
  REQUIRE(r.has_value());
  const geom::Brep& b = r.value();
  CHECK(b.vertex_count() == 3);
  CHECK(b.edge_count() == 2);
  CHECK(b.coedge_count() == 4);
  CHECK(geom::validate_topology(b).ok());
  const geom::CoedgeId first = b.loop(loop).start;
  const geom::EdgeId half0 = b.coedge(first).edge;
  const geom::EdgeId half1 = b.coedge(b.coedge(first).next).edge;
  CHECK(half0 != half1);
  std::size_t forward0 = 0;
  std::size_t reversed0 = 0;
  std::size_t forward1 = 0;
  std::size_t reversed1 = 0;
  for (const geom::CoedgeId id : b.coedge_ids()) {
    const geom::Coedge& c = b.coedge(id);
    if (c.edge == half0) (c.reversed ? reversed0 : forward0) += 1;
    if (c.edge == half1) (c.reversed ? reversed1 : forward1) += 1;
  }
  CHECK(forward0 == 1);
  CHECK(reversed0 == 1);
  CHECK(forward1 == 1);
  CHECK(reversed1 == 1);
}

TEST_CASE("split edge subdivides the curve parameter interval") {
  SquareFace s = make_square_face();
  const geom::EdgeId e0 = s.editor.brep().coedge(s.editor.brep().loop(s.loop).start).edge;
  const geom::Edge old = s.editor.brep().edge(e0);
  REQUIRE(s.editor.split_edge(e0, 0.5).has_value());
  const auto r = s.editor.commit();
  REQUIRE(r.has_value());
  const geom::Brep& b = r.value();
  const geom::Coedge& c0 = b.coedge(b.loop(s.loop).start);
  const geom::Coedge& c1 = b.coedge(c0.next);
  const geom::Edge& first = b.edge(c0.edge);
  const geom::Edge& second = b.edge(c1.edge);
  CHECK(first.curve == old.curve);
  CHECK(second.curve == old.curve);
  CHECK(first.t0 == old.t0);
  CHECK(first.t1 == 0.5);
  CHECK(second.t0 == 0.5);
  CHECK(second.t1 == old.t1);
  CHECK(first.v0 == old.v0);
  CHECK(first.v1 == second.v0);
  CHECK(second.v1 == old.v1);
  const geom::VertexId mid = first.v1;
  const geom::Point3d expected =
      std::get<geom::LineCurve3>(b.geometry().curves[first.curve]).eval(0.5);
  const geom::Point3d& actual = b.vertex(mid).point;
  CHECK(actual[0] == expected[0]);
  CHECK(actual[1] == expected[1]);
  CHECK(actual[2] == expected[2]);
}

TEST_CASE("split edge splits a forward pcurve at the parameter") {
  SingleEdgeFace s = make_single_edge_face(0.0);
  const std::uint32_t pcurve = s.editor.add_pcurve(geom::Pcurve2d{{0, 0}, {1, 1}});
  const std::array<geom::CoedgeSpec, 1> specs{geom::CoedgeSpec{s.edge, false, pcurve, 0.0}};
  const geom::LoopId loop = s.editor.add_loop(s.face, specs, true).value();
  REQUIRE(s.editor.split_edge(s.edge, 0.5).has_value());
  const auto r = s.editor.commit();
  REQUIRE(r.has_value());
  const geom::Brep& b = r.value();
  CHECK(b.geometry().pcurves.size() == 3u);
  const geom::Coedge& head = b.coedge(b.loop(loop).start);
  const geom::Coedge& tail = b.coedge(head.next);
  CHECK_FALSE(head.reversed);
  CHECK_FALSE(tail.reversed);
  const geom::Pcurve2d& first = b.geometry().pcurves[head.pcurve];
  const geom::Pcurve2d& second = b.geometry().pcurves[tail.pcurve];
  CHECK(first.points.front()[0] == 0.0);
  CHECK(first.points.front()[1] == 0.0);
  CHECK(first.points.back()[0] == 0.5);
  CHECK(first.points.back()[1] == 0.5);
  CHECK(second.points.front()[0] == 0.5);
  CHECK(second.points.front()[1] == 0.5);
  CHECK(second.points.back()[0] == 1.0);
  CHECK(second.points.back()[1] == 1.0);
}

TEST_CASE("split edge splits a reversed pcurve at the parameter") {
  SingleEdgeFace s = make_single_edge_face(0.0);
  const std::uint32_t pcurve = s.editor.add_pcurve(geom::Pcurve2d{{0, 0}, {1, 1}});
  const std::array<geom::CoedgeSpec, 1> specs{geom::CoedgeSpec{s.edge, true, pcurve, 0.0}};
  const geom::LoopId loop = s.editor.add_loop(s.face, specs, true).value();
  REQUIRE(s.editor.split_edge(s.edge, 0.25).has_value());
  const auto r = s.editor.commit();
  REQUIRE(r.has_value());
  const geom::Brep& b = r.value();
  CHECK(b.geometry().pcurves.size() == 3u);
  // Reversed use: the traversal starts on the high sub-edge, so the first coedge in the loop is
  // the e1 half, still marked reversed, with the pcurve segment that covers [0.25, 1].
  const geom::Coedge& head = b.coedge(b.loop(loop).start);
  const geom::Coedge& tail = b.coedge(head.next);
  CHECK(head.reversed);
  CHECK(tail.reversed);
  const geom::Edge& first_edge = b.edge(head.edge);
  const geom::Edge& second_edge = b.edge(tail.edge);
  CHECK(first_edge.t0 == 0.25);
  CHECK(first_edge.t1 == 1.0);
  CHECK(second_edge.t0 == 0.0);
  CHECK(second_edge.t1 == 0.25);
  const geom::Pcurve2d& high = b.geometry().pcurves[head.pcurve];
  const geom::Pcurve2d& low = b.geometry().pcurves[tail.pcurve];
  CHECK(low.points.front()[0] == 0.0);
  CHECK(low.points.front()[1] == 0.0);
  CHECK(low.points.back()[0] == 0.25);
  CHECK(low.points.back()[1] == 0.25);
  CHECK(high.points.front()[0] == 0.25);
  CHECK(high.points.front()[1] == 0.25);
  CHECK(high.points.back()[0] == 1.0);
  CHECK(high.points.back()[1] == 1.0);
}

TEST_CASE("split edge uses the interval-relative pcurve fraction") {
  geom::SolidEditor ed;
  const geom::VertexId v0 = ed.add_vertex({0, 0, 0}, 0.0).value();
  const geom::VertexId v1 = ed.add_vertex({1, 0, 0}, 0.0).value();
  const std::uint32_t line = ed.add_curve(geom::LineCurve3{{0, 0, 0}, {1, 0, 0}});
  const geom::EdgeId edge = ed.add_edge(v0, v1, line, 0.4, 0.9, 0.0).value();
  const std::uint32_t plane = ed.add_surface(geom::Plane3{{0, 0, 1}, 0.0});
  const geom::FaceId face = ed.add_face(plane, false, 0.0).value();
  // A non-45-degree pcurve so a wrong parameter (absolute vs interval-relative)
  // shows up in both coordinates.
  const std::uint32_t pcurve = ed.add_pcurve(geom::Pcurve2d{{0, 0}, {2, 1}});
  const std::array<geom::CoedgeSpec, 1> specs{geom::CoedgeSpec{edge, false, pcurve, 0.0}};
  const geom::LoopId loop = ed.add_loop(face, specs, true).value();
  REQUIRE(ed.split_edge(edge, 0.65).has_value());
  const auto r = ed.commit();
  REQUIRE(r.has_value());
  const geom::Brep& b = r.value();
  CHECK(b.geometry().pcurves.size() == 3u);
  const geom::Coedge& head = b.coedge(b.loop(loop).start);
  const geom::Coedge& tail = b.coedge(head.next);
  const geom::Pcurve2d& low = b.geometry().pcurves[head.pcurve];
  const geom::Pcurve2d& high = b.geometry().pcurves[tail.pcurve];
  // f = (0.65 - 0.4) / (0.9 - 0.4) = 0.5, so the split point is (1, 0.5).
  CHECK(low.points.front()[0] == 0.0);
  CHECK(low.points.front()[1] == 0.0);
  CHECK(low.points.back()[0] == doctest::Approx(1.0).epsilon(1e-12));
  CHECK(low.points.back()[1] == doctest::Approx(0.5).epsilon(1e-12));
  CHECK(high.points.front()[0] == doctest::Approx(1.0).epsilon(1e-12));
  CHECK(high.points.front()[1] == doctest::Approx(0.5).epsilon(1e-12));
  CHECK(high.points.back()[0] == 2.0);
  CHECK(high.points.back()[1] == 1.0);
}

TEST_CASE("join rebuilds the pcurve over the merged interval") {
  SingleEdgeFace s = make_single_edge_face(0.0);
  const std::uint32_t pcurve = s.editor.add_pcurve(geom::Pcurve2d{{0, 0}, {2, 1}});
  const std::array<geom::CoedgeSpec, 1> specs{geom::CoedgeSpec{s.edge, false, pcurve, 0.0}};
  const geom::LoopId loop = s.editor.add_loop(s.face, specs, true).value();
  const geom::VertexId mid = s.editor.split_edge(s.edge, 0.5).value();
  REQUIRE(s.editor.join_vertex(mid).has_value());
  const auto r = s.editor.commit();
  REQUIRE(r.has_value());
  const geom::Brep& b = r.value();
  const geom::Coedge& survivor = b.coedge(b.loop(loop).start);
  CHECK(survivor.pcurve != geom::kNoGeometry);
  const geom::Pcurve2d& merged = b.geometry().pcurves[survivor.pcurve];
  CHECK(merged.points.front()[0] == 0.0);
  CHECK(merged.points.front()[1] == 0.0);
  CHECK(merged.points.back()[0] == 2.0);
  CHECK(merged.points.back()[1] == 1.0);
  CHECK(geom::validate_topology(b).ok());
}

TEST_CASE("join rebuilds the pcurve for a reversed coedge") {
  SingleEdgeFace s = make_single_edge_face(0.0);
  const std::uint32_t pcurve = s.editor.add_pcurve(geom::Pcurve2d{{0, 0}, {2, 1}});
  const std::array<geom::CoedgeSpec, 1> specs{geom::CoedgeSpec{s.edge, true, pcurve, 0.0}};
  const geom::LoopId loop = s.editor.add_loop(s.face, specs, true).value();
  const geom::VertexId mid = s.editor.split_edge(s.edge, 0.5).value();
  REQUIRE(s.editor.join_vertex(mid).has_value());
  const auto r = s.editor.commit();
  REQUIRE(r.has_value());
  const geom::Brep& b = r.value();
  const geom::Coedge& survivor = b.coedge(b.loop(loop).start);
  CHECK(survivor.reversed);
  CHECK(survivor.pcurve != geom::kNoGeometry);
  // The pcurve always runs along the edge from t0 to t1, independent of the
  // coedge direction, so the merged segment is the original one.
  const geom::Pcurve2d& merged = b.geometry().pcurves[survivor.pcurve];
  CHECK(merged.points.front()[0] == 0.0);
  CHECK(merged.points.front()[1] == 0.0);
  CHECK(merged.points.back()[0] == 2.0);
  CHECK(merged.points.back()[1] == 1.0);
}

TEST_CASE("join drops the pcurve when either half has none") {
  geom::SolidEditor ed;
  const geom::VertexId v0 = ed.add_vertex({0, 0, 0}, 0.0).value();
  const geom::VertexId v1 = ed.add_vertex({0.5, 0, 0}, 0.0).value();
  const geom::VertexId v2 = ed.add_vertex({1, 0, 0}, 0.0).value();
  const std::uint32_t line = ed.add_curve(geom::LineCurve3{{0, 0, 0}, {1, 0, 0}});
  const geom::EdgeId first = ed.add_edge(v0, v1, line, 0.0, 0.5, 0.0).value();
  const geom::EdgeId second = ed.add_edge(v1, v2, line, 0.5, 1.0, 0.0).value();
  const std::uint32_t plane = ed.add_surface(geom::Plane3{{0, 0, 1}, 0.0});
  const geom::FaceId face = ed.add_face(plane, false, 0.0).value();
  const std::uint32_t pcurve = ed.add_pcurve(geom::Pcurve2d{{0, 0}, {0.5, 0.5}});
  const std::array<geom::CoedgeSpec, 2> specs{geom::CoedgeSpec{first, false, pcurve, 0.0},
                                              geom::CoedgeSpec{second, false}};
  const geom::LoopId loop = ed.add_loop(face, specs, true).value();
  REQUIRE(ed.join_vertex(v1).has_value());
  const auto r = ed.commit();
  REQUIRE(r.has_value());
  const geom::Brep& b = r.value();
  const geom::Coedge& survivor = b.coedge(b.loop(loop).start);
  CHECK(survivor.pcurve == geom::kNoGeometry);
}

TEST_CASE("split edge gives the new vertex and edges the edge tolerance") {
  SingleEdgeFace s = make_single_edge_face(1e-3);
  const std::array<geom::CoedgeSpec, 1> specs{geom::CoedgeSpec{s.edge, false}};
  const geom::LoopId loop = s.editor.add_loop(s.face, specs, true).value();
  const geom::VertexId mid = s.editor.split_edge(s.edge, 0.5).value();
  const auto r = s.editor.commit();
  REQUIRE(r.has_value());
  const geom::Brep& b = r.value();
  CHECK(b.vertex(mid).tol == 1e-3);
  const geom::Coedge& head = b.coedge(b.loop(loop).start);
  const geom::Coedge& tail = b.coedge(head.next);
  CHECK(b.edge(head.edge).tol == 1e-3);
  CHECK(b.edge(tail.edge).tol == 1e-3);
}

TEST_CASE("split rejects an invalid edge or parameter without touching staging") {
  SquareFace s = make_square_face();
  geom::SolidEditor& ed = s.editor;
  const geom::EdgeId e0 = ed.brep().coedge(ed.brep().loop(s.loop).start).edge;
  const std::uint32_t vertices = ed.brep().vertex_count();
  const std::uint32_t edges = ed.brep().edge_count();
  const std::uint32_t coedges = ed.brep().coedge_count();
  const double nan = std::numeric_limits<double>::quiet_NaN();
  const double inf = std::numeric_limits<double>::infinity();
  CHECK(ed.split_edge(geom::EdgeId::make(0, 2), 0.5).error().code ==
        geom::Errc::index_out_of_range);
  CHECK(ed.split_edge(e0, nan).error().code == geom::Errc::degenerate);
  CHECK(ed.split_edge(e0, inf).error().code == geom::Errc::degenerate);
  CHECK(ed.split_edge(e0, 0.0).error().code == geom::Errc::degenerate);
  CHECK(ed.split_edge(e0, 1.0).error().code == geom::Errc::degenerate);
  CHECK(ed.brep().vertex_count() == vertices);
  CHECK(ed.brep().edge_count() == edges);
  CHECK(ed.brep().coedge_count() == coedges);
  CHECK(ed.commit().has_value());
}

TEST_CASE("join rejects a dead vertex and a vertex with one incident edge") {
  geom::SolidEditor ed;
  const geom::VertexId v0 = ed.add_vertex({0, 0, 0}, 0.0).value();
  const geom::VertexId v1 = ed.add_vertex({1, 0, 0}, 0.0).value();
  const std::uint32_t line = ed.add_curve(geom::LineCurve3{{0, 0, 0}, {1, 0, 0}});
  REQUIRE(ed.add_edge(v0, v1, line, 0.0, 1.0, 0.0).has_value());
  const EditorCounts before = snapshot(ed);
  CHECK(ed.join_vertex(geom::VertexId::make(0, 2)).error().code ==
        geom::Errc::index_out_of_range);
  CHECK(ed.join_vertex(v0).error().code == geom::Errc::invalid_topology);
  CHECK(ed.join_vertex(v1).error().code == geom::Errc::invalid_topology);
  check_unchanged(ed, before);
}

TEST_CASE("join rejects adjacent edges on different curves") {
  geom::SolidEditor ed;
  const geom::VertexId v0 = ed.add_vertex({0, 0, 0}, 0.0).value();
  const geom::VertexId v1 = ed.add_vertex({1, 0, 0}, 0.0).value();
  const geom::VertexId v2 = ed.add_vertex({1, 1, 0}, 0.0).value();
  const std::uint32_t line_x = ed.add_curve(geom::LineCurve3{{0, 0, 0}, {1, 0, 0}});
  const std::uint32_t line_y = ed.add_curve(geom::LineCurve3{{0, 0, 0}, {0, 1, 0}});
  // The intervals are adjacent (0..1 then 1..2) and the coedges meet
  // head-to-tail, so only the curve-equality check can reject this join.
  const geom::EdgeId edge_a = ed.add_edge(v0, v1, line_x, 0.0, 1.0, 0.0).value();
  const geom::EdgeId edge_b = ed.add_edge(v1, v2, line_y, 1.0, 2.0, 0.0).value();
  const geom::EdgeId edge_c = ed.add_edge(v2, v0, geom::kNoGeometry, 0.0, 1.0, 0.0).value();
  const std::uint32_t plane = ed.add_surface(geom::Plane3{{0, 0, 1}, 0.0});
  const geom::FaceId face = ed.add_face(plane, false, 0.0).value();
  const std::array<geom::CoedgeSpec, 3> specs{geom::CoedgeSpec{edge_a, false},
                                              geom::CoedgeSpec{edge_b, false},
                                              geom::CoedgeSpec{edge_c, false}};
  REQUIRE(ed.add_loop(face, specs, true).has_value());
  const EditorCounts before = snapshot(ed);
  CHECK(ed.join_vertex(v1).error().code == geom::Errc::degenerate);
  check_unchanged(ed, before);
}

TEST_CASE("join rejects edges whose parameter intervals are not adjacent") {
  geom::SolidEditor ed;
  const geom::VertexId v0 = ed.add_vertex({0, 0, 0}, 0.0).value();
  const geom::VertexId v1 = ed.add_vertex({1, 0, 0}, 0.0).value();
  const geom::VertexId v2 = ed.add_vertex({2, 0, 0}, 0.0).value();
  const std::uint32_t line = ed.add_curve(geom::LineCurve3{{0, 0, 0}, {1, 0, 0}});
  const geom::EdgeId edge_a = ed.add_edge(v0, v1, line, 0.0, 1.0, 0.0).value();
  const geom::EdgeId edge_b = ed.add_edge(v1, v2, line, 2.0, 3.0, 0.0).value();
  const geom::EdgeId edge_c = ed.add_edge(v2, v0, geom::kNoGeometry, 0.0, 1.0, 0.0).value();
  const std::uint32_t plane = ed.add_surface(geom::Plane3{{0, 0, 1}, 0.0});
  const geom::FaceId face = ed.add_face(plane, false, 0.0).value();
  const std::array<geom::CoedgeSpec, 3> specs{geom::CoedgeSpec{edge_a, false},
                                              geom::CoedgeSpec{edge_b, false},
                                              geom::CoedgeSpec{edge_c, false}};
  REQUIRE(ed.add_loop(face, specs, true).has_value());
  const EditorCounts before = snapshot(ed);
  CHECK(ed.join_vertex(v1).error().code == geom::Errc::degenerate);
  check_unchanged(ed, before);
}

TEST_CASE("join merges the adjacent parameter intervals") {
  SquareFace s = make_square_face();
  const geom::EdgeId e0 = s.editor.brep().coedge(s.editor.brep().loop(s.loop).start).edge;
  const geom::VertexId mid = s.editor.split_edge(e0, 0.25).value();
  REQUIRE(s.editor.join_vertex(mid).has_value());
  const auto r = s.editor.commit();
  REQUIRE(r.has_value());
  const geom::Brep& b = r.value();
  const geom::Coedge& survivor = b.coedge(b.loop(s.loop).start);
  const geom::Edge& merged = b.edge(survivor.edge);
  CHECK(merged.t0 == 0.0);
  CHECK(merged.t1 == 1.0);
  CHECK(merged.curve == 0u);
  CHECK(merged.v0 == s.vertices[0]);
  CHECK(merged.v1 == s.vertices[1]);
  CHECK(b.vertex_count() == 4);
  CHECK(b.edge_count() == 4);
}

TEST_CASE("join rejects a vertex with four coedge uses") {
  geom::SolidEditor ed;
  const geom::VertexId v0 = ed.add_vertex({0, 0, 0}, 0.0).value();
  const geom::VertexId v1 = ed.add_vertex({1, 0, 0}, 0.0).value();
  const std::uint32_t line = ed.add_curve(geom::LineCurve3{{0, 0, 0}, {1, 0, 0}});
  const geom::EdgeId edge = ed.add_edge(v0, v1, line, 0.0, 1.0, 0.0).value();
  const std::uint32_t plane = ed.add_surface(geom::Plane3{{0, 0, 1}, 0.0});
  const geom::FaceId face = ed.add_face(plane, false, 0.0).value();
  const std::array<geom::CoedgeSpec, 2> specs{geom::CoedgeSpec{edge, false},
                                              geom::CoedgeSpec{edge, true}};
  REQUIRE(ed.add_loop(face, specs, true).has_value());
  const geom::VertexId mid = ed.split_edge(edge, 0.5).value();
  const EditorCounts before = snapshot(ed);
  // Two incident edges, but four coedge uses: outside the KEV contract of this task.
  CHECK(ed.join_vertex(mid).error().code == geom::Errc::invalid_topology);
  check_unchanged(ed, before);
}

TEST_CASE("join failure with three incident edges leaves staging unchanged") {
  SquareFace s = make_square_face();
  const geom::VertexId extra = s.editor.add_vertex({0.5, -1.0, 0.0}, 0.0).value();
  const std::uint32_t line = s.editor.add_curve(geom::LineCurve3{{1, 0, 0}, {0, 1, 0}});
  REQUIRE(s.editor.add_edge(s.vertices[1], extra, line, 0.0, 1.0, 0.0).has_value());
  const EditorCounts before = snapshot(s.editor);
  CHECK(s.editor.join_vertex(s.vertices[1]).error().code == geom::Errc::invalid_topology);
  check_unchanged(s.editor, before);
}
