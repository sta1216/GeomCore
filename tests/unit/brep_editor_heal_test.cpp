#include <doctest/doctest.h>

#include <array>
#include <cstddef>
#include <cstdint>
#include <string_view>
#include <vector>

#include <geom/brep/solid_editor.hpp>
#include <geom/brep/topology.hpp>
#include <geom/core/error.hpp>
#include <geom/curves/analytic_curve.hpp>
#include <geom/surfaces/plane.hpp>

namespace {

using namespace geom;

// Four boundary edges of an axis-aligned unit quad on z = 0, running
// counter-clockwise. The face and loop are created but no shell: the callers
// decide how the face joins a shell, which keeps the sewing fixtures small.
struct QuadFace {
  std::array<VertexId, 4> vertices;
  std::array<EdgeId, 4> edges;
  FaceId face{};
  LoopId loop{};
};

QuadFace add_quad_face(SolidEditor& editor, std::uint32_t surface, double x, double y,
                       double side, double vertex_tol, double edge_tol, double face_tol,
                       std::uint32_t curve = kNoGeometry, double t1 = 1.0) {
  QuadFace quad;
  quad.vertices = {editor.add_vertex(Point3d{x, y, 0.0}, vertex_tol).value(),
                   editor.add_vertex(Point3d{x + side, y, 0.0}, vertex_tol).value(),
                   editor.add_vertex(Point3d{x + side, y + side, 0.0}, vertex_tol).value(),
                   editor.add_vertex(Point3d{x, y + side, 0.0}, vertex_tol).value()};
  for (std::size_t i = 0; i < quad.edges.size(); ++i) {
    quad.edges[i] = editor
                        .add_edge(quad.vertices[i], quad.vertices[(i + 1u) % 4u], curve, 0.0, t1,
                                  edge_tol)
                        .value();
  }
  quad.face = editor.add_face(surface, false, face_tol).value();
  std::vector<CoedgeSpec> loop;
  for (const EdgeId edge : quad.edges) loop.push_back(CoedgeSpec{edge});
  quad.loop = editor.add_loop(quad.face, loop, true).value();
  return quad;
}

TEST_CASE("heal editor: merging a duplicate vertex retargets its edges") {
  SolidEditor editor;
  const std::uint32_t surface = editor.add_surface(Plane3{});
  const VertexId a = editor.add_vertex(Point3d{0.0, 0.0, 0.0}, 0.1).value();
  const VertexId b = editor.add_vertex(Point3d{1.0, 0.0, 0.0}, 0.1).value();
  const VertexId c = editor.add_vertex(Point3d{1.0, 1.0, 0.0}, 0.1).value();
  const VertexId d = editor.add_vertex(Point3d{0.0, 1.0, 0.0}, 0.1).value();
  const VertexId duplicate = editor.add_vertex(Point3d{0.5, 0.0, 0.0}, 0.2).value();
  const EdgeId e0 = editor.add_edge(a, b, kNoGeometry, 0.0, 1.0, 0.5).value();
  const EdgeId e1 = editor.add_edge(b, c, kNoGeometry, 0.0, 1.0, 0.5).value();
  const EdgeId e2 = editor.add_edge(c, d, kNoGeometry, 0.0, 1.0, 0.5).value();
  const EdgeId e3 = editor.add_edge(d, duplicate, kNoGeometry, 0.0, 1.0, 0.5).value();
  const FaceId face = editor.add_face(surface, false, 0.5).value();
  std::vector<CoedgeSpec> loop{CoedgeSpec{e0}, CoedgeSpec{e1}, CoedgeSpec{e2}, CoedgeSpec{e3}};
  REQUIRE(editor.add_loop(face, loop, true).has_value());
  const ShellId shell =
      editor.add_shell(std::array<FaceId, 1>{face}, false).value();
  REQUIRE(editor.add_solid(std::array<ShellId, 1>{shell}, 0.5).has_value());

  REQUIRE(editor.merge_vertices(a, duplicate).has_value());

  const Brep& staged = editor.brep();
  CHECK(staged.vertex(a).tol == 0.5);  // the weld distance dominates both tols
  CHECK_FALSE(staged.vertex_alive(duplicate));
  CHECK(staged.edge(e3).v1 == a);
  CHECK(staged.edge(e3).v0 == d);
  REQUIRE(editor.commit().has_value());
}

TEST_CASE("heal editor: a vertex merge rejects itself and stale handles") {
  SolidEditor editor;
  const VertexId vertex = editor.add_vertex(Point3d{}, 0.0).value();
  CHECK(editor.merge_vertices(vertex, vertex).error().code == Errc::invalid_topology);
  CHECK(editor.merge_vertices(vertex, VertexId{}).error().code == Errc::invalid_topology);
  CHECK(editor.merge_vertices(VertexId{}, vertex).error().code == Errc::invalid_topology);
  CHECK(editor.brep().vertex_alive(vertex));
}

TEST_CASE("heal editor: a vertex merge cannot create a self-loop edge") {
  SolidEditor editor;
  const VertexId a = editor.add_vertex(Point3d{0.0, 0.0, 0.0}, 0.0).value();
  const VertexId b = editor.add_vertex(Point3d{1.0, 0.0, 0.0}, 0.0).value();
  const EdgeId edge = editor.add_edge(a, b, kNoGeometry, 0.0, 1.0, 0.0).value();

  CHECK(editor.merge_vertices(a, b).error().code == Errc::unsupported);
  CHECK(editor.brep().vertex_alive(b));
  CHECK(editor.brep().edge(edge).v0 == a);
  CHECK(editor.brep().edge(edge).v1 == b);
}

TEST_CASE("heal editor: merging welded boundary edges preserves direction and tolerance") {
  SolidEditor editor;
  const std::uint32_t surface = editor.add_surface(Plane3{});
  const std::uint32_t curve = editor.add_curve(LineCurve3{});
  const QuadFace left =
      add_quad_face(editor, surface, 0.0, 0.0, 1.0, 0.001, 0.003, 0.003);
  const QuadFace right =
      add_quad_face(editor, surface, 1.0, 0.0, 1.0, 0.001, 0.009, 0.009, curve, 2.0);
  const ShellId left_shell =
      editor.add_shell(std::array<FaceId, 1>{left.face}, false).value();
  const ShellId right_shell =
      editor.add_shell(std::array<FaceId, 1>{right.face}, false).value();
  REQUIRE(editor
              .add_solid(std::array<ShellId, 1>{left_shell}, 0.009)
              .has_value());
  REQUIRE(editor
              .add_solid(std::array<ShellId, 1>{right_shell}, 0.009)
              .has_value());

  // left.edges[1] runs (1,0) -> (1,1); right.edges[3] runs (1,1) -> (1,0), so
  // the shared segment is traversed in opposite directions.
  CHECK(editor.merge_edges(left.edges[1], right.edges[3]).error().code ==
        Errc::invalid_topology);
  REQUIRE(editor.merge_vertices(left.vertices[1], right.vertices[0]).has_value());
  REQUIRE(editor.merge_vertices(left.vertices[2], right.vertices[3]).has_value());
  REQUIRE(editor.merge_edges(left.edges[1], right.edges[3]).has_value());

  const Brep& staged = editor.brep();
  CHECK_FALSE(staged.edge_alive(right.edges[3]));
  CHECK(staged.edge(left.edges[1]).tol == 0.003);
  CHECK(staged.edge(left.edges[1]).curve == curve);
  CHECK(staged.edge(left.edges[1]).t0 == 0.0);
  CHECK(staged.edge(left.edges[1]).t1 == 2.0);

  int forward = 0;
  int reversed = 0;
  for (const CoedgeId coedge : staged.coedge_ids()) {
    if (staged.coedge(coedge).edge != left.edges[1]) continue;
    if (staged.coedge(coedge).reversed) {
      ++reversed;
    } else {
      ++forward;
    }
  }
  CHECK(forward == 1);
  CHECK(reversed == 1);
  REQUIRE(editor.commit().has_value());
}

TEST_CASE("heal editor: closing a shell requires every edge used twice") {
  SolidEditor editor;
  const std::uint32_t surface = editor.add_surface(Plane3{});
  const QuadFace sheet = add_quad_face(editor, surface, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0);
  const ShellId open_shell =
      editor.add_shell(std::array<FaceId, 1>{sheet.face}, false).value();
  const SolidId open_solid =
      editor.add_solid(std::array<ShellId, 1>{open_shell}, 0.0).value();
  CHECK(editor.brep().solid(open_solid).tol == 0.0);

  CHECK(editor.set_shell_closed(open_shell, true).error().code == Errc::not_closed);
  CHECK_FALSE(editor.brep().shell(open_shell).closed);
  CHECK(editor.commit().has_value());

  // Two faces using the same four edges in opposite directions form a closed
  // surface topologically: every edge is used exactly twice.
  SolidEditor doubled;
  const std::uint32_t doubled_surface = doubled.add_surface(Plane3{});
  const QuadFace front =
      add_quad_face(doubled, doubled_surface, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0);
  const FaceId back = doubled.add_face(doubled_surface, false, 0.0).value();
  std::vector<CoedgeSpec> back_loop{CoedgeSpec{front.edges[3], true},
                                    CoedgeSpec{front.edges[2], true},
                                    CoedgeSpec{front.edges[1], true},
                                    CoedgeSpec{front.edges[0], true}};
  REQUIRE(doubled.add_loop(back, back_loop, true).has_value());
  const ShellId closed_shell =
      doubled.add_shell(std::array<FaceId, 2>{front.face, back}, false).value();
  REQUIRE(doubled.add_solid(std::array<ShellId, 1>{closed_shell}, 0.0).has_value());

  REQUIRE(doubled.set_shell_closed(closed_shell, true).has_value());
  CHECK(doubled.brep().shell(closed_shell).closed);
  REQUIRE(doubled.commit().has_value());

  const ShellId stale =
      ShellId::make(open_shell.index(),
                    static_cast<std::uint8_t>(open_shell.generation() + 1u));
  CHECK(editor.set_shell_closed(stale, true).error().code == Errc::invalid_topology);
}

TEST_CASE("heal editor: pcurve assignment round-trips and validates") {
  SolidEditor editor;
  const std::uint32_t surface = editor.add_surface(Plane3{});
  const QuadFace quad = add_quad_face(editor, surface, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0);
  const CoedgeId coedge = editor.brep().coedge_ids().front();
  const std::uint32_t pcurve =
      editor.add_pcurve(Pcurve2d{{Point2d{0.0, 0.0}, Point2d{1.0, 1.0}}});

  REQUIRE(editor.set_coedge_pcurve(coedge, pcurve).has_value());
  CHECK(editor.brep().coedge(coedge).pcurve == pcurve);
  CHECK(editor.brep().coedge_count() == 4u);
  CHECK(editor.brep().loop(quad.loop).start == coedge);

  const std::uint32_t out_of_range =
      static_cast<std::uint32_t>(editor.brep().geometry().pcurves.size());
  CHECK(editor.set_coedge_pcurve(coedge, out_of_range).error().code ==
        Errc::index_out_of_range);
  CHECK(editor.brep().coedge(coedge).pcurve == pcurve);

  REQUIRE(editor.set_coedge_pcurve(coedge, kNoGeometry).has_value());
  CHECK(editor.brep().coedge(coedge).pcurve == kNoGeometry);

  const CoedgeId stale =
      CoedgeId::make(coedge.index(), static_cast<std::uint8_t>(coedge.generation() + 1u));
  CHECK(editor.set_coedge_pcurve(stale, pcurve).error().code == Errc::invalid_topology);
}

TEST_CASE("heal editor: edge geometry round-trips and validates") {
  SolidEditor editor;
  const std::uint32_t surface = editor.add_surface(Plane3{});
  const QuadFace quad = add_quad_face(editor, surface, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0);
  const EdgeId edge = quad.edges[0];
  const std::uint32_t curve = editor.add_curve(LineCurve3{});

  REQUIRE(editor.set_edge_geometry(edge, curve, 0.0, 2.0).has_value());
  CHECK(editor.brep().edge(edge).curve == curve);
  CHECK(editor.brep().edge(edge).t0 == 0.0);
  CHECK(editor.brep().edge(edge).t1 == 2.0);

  const std::uint32_t out_of_range =
      static_cast<std::uint32_t>(editor.brep().geometry().curves.size());
  CHECK(editor.set_edge_geometry(edge, out_of_range, 0.0, 1.0).error().code ==
        Errc::index_out_of_range);
  CHECK(editor.brep().edge(edge).curve == curve);

  CHECK(editor.set_edge_geometry(edge, kNoGeometry, 1.0, 1.0).error().code ==
        Errc::degenerate);
  CHECK(editor.brep().edge(edge).curve == curve);
  REQUIRE(editor.set_edge_geometry(edge, kNoGeometry, 0.0, 1.0).has_value());
  CHECK(editor.brep().edge(edge).curve == kNoGeometry);

  const EdgeId stale =
      EdgeId::make(edge.index(), static_cast<std::uint8_t>(edge.generation() + 1u));
  CHECK(editor.set_edge_geometry(stale, curve, 0.0, 1.0).error().code ==
        Errc::invalid_topology);
}

TEST_CASE("heal editor: raising tolerances propagates through the tolerance chain") {
  SolidEditor editor;
  const std::uint32_t surface = editor.add_surface(Plane3{});
  const QuadFace quad = add_quad_face(editor, surface, 0.0, 0.0, 1.0, 0.02, 0.0, 0.0);
  const ShellId shell = editor.add_shell(std::array<FaceId, 1>{quad.face}, false).value();
  const SolidId solid = editor.add_solid(std::array<ShellId, 1>{shell}, 0.0).value();

  const Result<double> raised = editor.raise_tolerances(0.02);
  REQUIRE(raised.has_value());
  CHECK(raised.value() == 0.02);
  const Brep& staged = editor.brep();
  CHECK(staged.edge(quad.edges[0]).tol == 0.02);
  CHECK(staged.face(quad.face).tol == 0.02);
  CHECK(staged.solid(solid).tol == 0.02);
  CHECK(editor.raise_tolerances(0.02).value() == 0.0);
  REQUIRE(editor.commit().has_value());
}

TEST_CASE("heal editor: a tolerance cap below a needed value changes nothing") {
  SolidEditor editor;
  const std::uint32_t surface = editor.add_surface(Plane3{});
  const QuadFace quad = add_quad_face(editor, surface, 0.0, 0.0, 1.0, 0.02, 0.0, 0.0);
  const ShellId shell = editor.add_shell(std::array<FaceId, 1>{quad.face}, false).value();
  const SolidId solid = editor.add_solid(std::array<ShellId, 1>{shell}, 0.0).value();

  const Result<double> rejected = editor.raise_tolerances(0.01);
  CHECK(rejected.error().code == Errc::tolerance_violation);
  CHECK(std::string_view(rejected.error().message) == "heal_brep: tolerance cap exceeded");
  CHECK(editor.brep().edge(quad.edges[0]).tol == 0.0);
  CHECK(editor.brep().face(quad.face).tol == 0.0);
  CHECK(editor.brep().solid(solid).tol == 0.0);
}

TEST_CASE("heal editor: the tolerance cap is enforced at the edge level") {
  SolidEditor editor;
  const std::uint32_t surface = editor.add_surface(Plane3{});
  // Only the edges must rise (to the vertex tolerance); the face and solid are
  // already above the needed value, so a skipped edge-level cap check would let
  // this pass silently.
  const QuadFace quad = add_quad_face(editor, surface, 0.0, 0.0, 1.0, 0.025, 0.0, 0.03);
  const ShellId shell = editor.add_shell(std::array<FaceId, 1>{quad.face}, false).value();
  CHECK(editor.add_solid(std::array<ShellId, 1>{shell}, 0.03).has_value());

  const Result<double> rejected = editor.raise_tolerances(0.02);
  CHECK(rejected.error().code == Errc::tolerance_violation);
  CHECK(editor.brep().edge(quad.edges[0]).tol == 0.0);
}

}  // namespace
