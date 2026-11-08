#include <doctest/doctest.h>

#include <array>
#include <cstdint>
#include <limits>
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
  CHECK(first.t1 == doctest::Approx(0.5));
  CHECK(second.t0 == doctest::Approx(0.5));
  CHECK(second.t1 == old.t1);
  CHECK(first.v0 == old.v0);
  CHECK(first.v1 == second.v0);
  CHECK(second.v1 == old.v1);
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
  CHECK(first.p0[0] == doctest::Approx(0.0));
  CHECK(first.p0[1] == doctest::Approx(0.0));
  CHECK(first.p1[0] == doctest::Approx(0.5));
  CHECK(first.p1[1] == doctest::Approx(0.5));
  CHECK(second.p0[0] == doctest::Approx(0.5));
  CHECK(second.p0[1] == doctest::Approx(0.5));
  CHECK(second.p1[0] == doctest::Approx(1.0));
  CHECK(second.p1[1] == doctest::Approx(1.0));
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
  CHECK(first_edge.t0 == doctest::Approx(0.25));
  CHECK(first_edge.t1 == doctest::Approx(1.0));
  CHECK(second_edge.t0 == doctest::Approx(0.0));
  CHECK(second_edge.t1 == doctest::Approx(0.25));
  const geom::Pcurve2d& high = b.geometry().pcurves[head.pcurve];
  const geom::Pcurve2d& low = b.geometry().pcurves[tail.pcurve];
  CHECK(low.p0[0] == doctest::Approx(0.0));
  CHECK(low.p0[1] == doctest::Approx(0.0));
  CHECK(low.p1[0] == doctest::Approx(0.25));
  CHECK(low.p1[1] == doctest::Approx(0.25));
  CHECK(high.p0[0] == doctest::Approx(0.25));
  CHECK(high.p0[1] == doctest::Approx(0.25));
  CHECK(high.p1[0] == doctest::Approx(1.0));
  CHECK(high.p1[1] == doctest::Approx(1.0));
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
  CHECK(ed.join_vertex(geom::VertexId::make(0, 2)).error().code ==
        geom::Errc::index_out_of_range);
  CHECK(ed.join_vertex(v0).error().code == geom::Errc::invalid_topology);
  CHECK(ed.join_vertex(v1).error().code == geom::Errc::invalid_topology);
}

TEST_CASE("join rejects adjacent edges on different curves") {
  geom::SolidEditor ed;
  const geom::VertexId v0 = ed.add_vertex({0, 0, 0}, 0.0).value();
  const geom::VertexId v1 = ed.add_vertex({1, 0, 0}, 0.0).value();
  const geom::VertexId v2 = ed.add_vertex({1, 1, 0}, 0.0).value();
  const std::uint32_t line_x = ed.add_curve(geom::LineCurve3{{0, 0, 0}, {1, 0, 0}});
  const std::uint32_t line_y = ed.add_curve(geom::LineCurve3{{0, 0, 0}, {0, 1, 0}});
  const geom::EdgeId edge_a = ed.add_edge(v0, v1, line_x, 0.0, 1.0, 0.0).value();
  const geom::EdgeId edge_b = ed.add_edge(v1, v2, line_y, 0.0, 1.0, 0.0).value();
  const geom::EdgeId edge_c = ed.add_edge(v2, v0, geom::kNoGeometry, 0.0, 1.0, 0.0).value();
  const std::uint32_t plane = ed.add_surface(geom::Plane3{{0, 0, 1}, 0.0});
  const geom::FaceId face = ed.add_face(plane, false, 0.0).value();
  const std::array<geom::CoedgeSpec, 3> specs{geom::CoedgeSpec{edge_a, false},
                                              geom::CoedgeSpec{edge_b, false},
                                              geom::CoedgeSpec{edge_c, false}};
  REQUIRE(ed.add_loop(face, specs, true).has_value());
  CHECK(ed.join_vertex(v1).error().code == geom::Errc::degenerate);
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
  CHECK(ed.join_vertex(v1).error().code == geom::Errc::degenerate);
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
  // Two incident edges, but four coedge uses: outside the KEV contract of this task.
  CHECK(ed.join_vertex(mid).error().code == geom::Errc::invalid_topology);
}
