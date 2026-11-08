#include <doctest/doctest.h>

#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <string_view>

#include <geom/brep/solid_editor.hpp>
#include <geom/brep/validate.hpp>

namespace {

geom::SolidEditor make_square_editor(double vertex_tol, double edge_tol,
                                     double face_tol, double solid_tol) {
  geom::SolidEditor ed;
  const geom::Point3d pts[4] = {{0, 0, 0}, {1, 0, 0}, {1, 1, 0}, {0, 1, 0}};
  std::array<geom::VertexId, 4> v{};
  for (int i = 0; i < 4; ++i) v[i] = ed.add_vertex(pts[i], vertex_tol).value();
  const std::uint32_t plane = ed.add_surface(geom::Plane3{{0, 0, 1}, 0.0});
  const std::uint32_t line = ed.add_curve(geom::LineCurve3{{0, 0, 0}, {1, 0, 0}});
  std::array<geom::EdgeId, 4> e{};
  for (int i = 0; i < 4; ++i) {
    e[i] = ed.add_edge(v[i], v[(i + 1) % 4], line, 0.0, 1.0, edge_tol).value();
  }
  const geom::FaceId f = ed.add_face(plane, false, face_tol).value();
  std::array<geom::CoedgeSpec, 4> specs{};
  for (int i = 0; i < 4; ++i) specs[i] = geom::CoedgeSpec{e[i], false};
  ed.add_loop(f, specs, true).value();
  const auto sh = ed.add_shell(std::array<geom::FaceId, 1>{f}, false).value();
  ed.add_solid(std::array<geom::ShellId, 1>{sh}, solid_tol).value();
  return ed;
}

// A single square face with per-edge and per-coedge tolerances and no shell or
// solid. Used to attribute rule-3 findings to a non-zero edge/coedge index: with
// the first use valid and a later one violating, entity attribution is
// distinguishable from the index-0 coincidence of the shared helper.
geom::SolidEditor make_square_with_tolerances(const std::array<double, 4>& edge_tols,
                                              double face_tol,
                                              const std::array<double, 4>& coedge_tols) {
  geom::SolidEditor ed;
  const geom::Point3d pts[4] = {{0, 0, 0}, {1, 0, 0}, {1, 1, 0}, {0, 1, 0}};
  std::array<geom::VertexId, 4> v{};
  for (int i = 0; i < 4; ++i) v[i] = ed.add_vertex(pts[i], 0.0).value();
  const std::uint32_t plane = ed.add_surface(geom::Plane3{{0, 0, 1}, 0.0});
  const std::uint32_t line = ed.add_curve(geom::LineCurve3{{0, 0, 0}, {1, 0, 0}});
  std::array<geom::EdgeId, 4> e{};
  for (int i = 0; i < 4; ++i) {
    e[i] = ed.add_edge(v[i], v[(i + 1) % 4], line, 0.0, 1.0, edge_tols[i]).value();
  }
  const geom::FaceId f = ed.add_face(plane, false, face_tol).value();
  std::array<geom::CoedgeSpec, 4> specs{};
  for (int i = 0; i < 4; ++i) {
    specs[i] = geom::CoedgeSpec{e[i], false, geom::kNoGeometry, coedge_tols[i]};
  }
  ed.add_loop(f, specs, true).value();
  return ed;
}

// Live entity counts of the staging area, for asserting that a rejected add_*
// call did not mutate it.
struct EntityCounts {
  std::size_t vertices;
  std::size_t edges;
  std::size_t coedges;
  std::size_t loops;
  std::size_t faces;
  std::size_t shells;
  std::size_t solids;
};

EntityCounts counts_of(const geom::SolidEditor& ed) {
  const geom::Brep& brep = ed.brep();
  return EntityCounts{brep.vertex_count(), brep.edge_count(), brep.coedge_count(),
                      brep.loop_count(),   brep.face_count(), brep.shell_count(),
                      brep.solid_count()};
}

bool same_counts(const EntityCounts& lhs, const EntityCounts& rhs) {
  return lhs.vertices == rhs.vertices && lhs.edges == rhs.edges &&
         lhs.coedges == rhs.coedges && lhs.loops == rhs.loops && lhs.faces == rhs.faces &&
         lhs.shells == rhs.shells && lhs.solids == rhs.solids;
}

}  // namespace

TEST_CASE("edge tolerance must dominate its vertices") {
  auto ed = make_square_editor(/*vertex=*/0.1, /*edge=*/0.05, /*face=*/0.2,
                               /*solid=*/0.2);
  const auto r = ed.commit();
  REQUIRE_FALSE(r.has_value());
  CHECK(r.error().code == geom::Errc::tolerance_violation);
}

TEST_CASE("solid tolerance must dominate its faces") {
  auto ed = make_square_editor(0.0, 0.0, /*face=*/1.0, /*solid=*/0.1);
  const auto r = ed.commit();
  REQUIRE_FALSE(r.has_value());
  CHECK(r.error().code == geom::Errc::tolerance_violation);
}

TEST_CASE("tolerance growth beyond ten times is a violation") {
  auto ed = make_square_editor(0.0, /*edge=*/0.1, /*face=*/2.0, /*solid=*/2.0);
  const auto r = ed.commit();
  REQUIRE_FALSE(r.has_value());
  CHECK(r.error().code == geom::Errc::tolerance_violation);
  CHECK(std::string_view(r.error().message).find("growth") != std::string_view::npos);
}

TEST_CASE("valid tolerance chain is accepted") {
  auto ed = make_square_editor(0.01, 0.02, 0.03, 0.03);
  const auto r = ed.commit();
  REQUIRE(r.has_value());
  CHECK(geom::validate_tolerances(r.value()).ok());
}

TEST_CASE("negative tolerance is rejected at construction") {
  auto ed = make_square_editor(0.0, 0.0, 0.0, 0.0);
  CHECK_FALSE(ed.add_vertex({0, 0, 0}, -1.0).has_value());
  CHECK_FALSE(ed.add_face(0, false, -0.5).has_value());
}

TEST_CASE("tolerance exactly at the growth factor is accepted") {
  // edge.tol == kToleranceGrowthFactor * vertex.tol must not trip the growth
  // alarm: the check is a strict greater-than.
  auto ed = make_square_editor(/*vertex=*/0.01, /*edge=*/0.1, /*face=*/0.1,
                               /*solid=*/0.1);
  const auto r = ed.commit();
  REQUIRE(r.has_value());
  CHECK(geom::validate_tolerances(r.value()).ok());
}

TEST_CASE("tolerance just above the growth factor is rejected") {
  const double parent = std::nextafter(1.0, 2.0);
  auto ed = make_square_editor(/*vertex=*/0.1, /*edge=*/parent, /*face=*/parent,
                               /*solid=*/parent);
  const auto r = ed.commit();
  REQUIRE_FALSE(r.has_value());
  CHECK(r.error().code == geom::Errc::tolerance_violation);
  CHECK(std::string_view(r.error().message).find("growth") != std::string_view::npos);
}

TEST_CASE("face tolerance below an edge tolerance is a violation") {
  auto ed = make_square_with_tolerances(
      /*edge=*/{0.05, 0.2, 0.2, 0.2}, /*face=*/0.1, /*coedge=*/{0.0, 0.0, 0.0, 0.0});
  const auto r = ed.commit();
  REQUIRE_FALSE(r.has_value());
  CHECK(r.error().code == geom::Errc::tolerance_violation);
  // The first edge is fine (0.05 <= 0.1); the second edge is the first one whose
  // tolerance exceeds its face. The finding is attributed to that edge.
  CHECK(r.error().context_a == 1u);
}

TEST_CASE("coedge tolerance above its face tolerance is a violation") {
  auto ed = make_square_with_tolerances(
      /*edge=*/{0.02, 0.02, 0.02, 0.02}, /*face=*/0.1,
      /*coedge=*/{0.0, 0.2, 0.0, 0.0});
  const auto r = ed.commit();
  REQUIRE_FALSE(r.has_value());
  CHECK(r.error().code == geom::Errc::tolerance_violation);
  // All face/edge tolerances are valid here; the second coedge is the first one
  // whose tolerance exceeds its face. The finding is attributed to that coedge.
  CHECK(r.error().context_a == 1u);
}

TEST_CASE("edge growth beyond ten times its vertex is a violation") {
  auto ed = make_square_editor(/*vertex=*/0.01, /*edge=*/0.2, /*face=*/0.2,
                               /*solid=*/0.2);
  const auto r = ed.commit();
  REQUIRE_FALSE(r.has_value());
  CHECK(r.error().code == geom::Errc::tolerance_violation);
  CHECK(std::string_view(r.error().message).find("growth") != std::string_view::npos);
}

TEST_CASE("solid growth beyond ten times a face is a violation") {
  auto ed = make_square_editor(/*vertex=*/0.01, /*edge=*/0.01, /*face=*/0.01,
                               /*solid=*/0.2);
  const auto r = ed.commit();
  REQUIRE_FALSE(r.has_value());
  CHECK(r.error().code == geom::Errc::tolerance_violation);
  CHECK(std::string_view(r.error().message).find("growth") != std::string_view::npos);
}

TEST_CASE("non-finite or negative tolerances are rejected at construction") {
  auto ed = make_square_editor(0.0, 0.0, 0.0, 0.0);
  const double nan = std::numeric_limits<double>::quiet_NaN();
  const geom::Brep& brep = ed.brep();
  const auto vertices = brep.vertex_ids();
  const auto edges = brep.edge_ids();
  const auto faces = brep.face_ids();
  const auto shells = brep.shell_ids();
  const EntityCounts before = counts_of(ed);

  const auto rejected_without_mutation = [&ed, &before](const auto& result) {
    CHECK_FALSE(result.has_value());
    CHECK(result.error().code == geom::Errc::tolerance_violation);
    CHECK(same_counts(counts_of(ed), before));
  };

  rejected_without_mutation(ed.add_edge(vertices[0], vertices[1], 0u, 0.0, 1.0, nan));
  rejected_without_mutation(ed.add_edge(vertices[0], vertices[1], 0u, 0.0, 1.0, -0.1));
  rejected_without_mutation(ed.add_face(0u, false, nan));
  rejected_without_mutation(ed.add_face(0u, false, -0.5));
  const std::array<geom::CoedgeSpec, 1> nan_coedge{
      geom::CoedgeSpec{edges[0], false, geom::kNoGeometry, nan}};
  rejected_without_mutation(ed.add_loop(faces[0], nan_coedge, false));
  const std::array<geom::CoedgeSpec, 1> negative_coedge{
      geom::CoedgeSpec{edges[0], false, geom::kNoGeometry, -0.2}};
  rejected_without_mutation(ed.add_loop(faces[0], negative_coedge, false));
  rejected_without_mutation(ed.add_solid(std::array<geom::ShellId, 1>{shells[0]}, nan));
  rejected_without_mutation(ed.add_solid(std::array<geom::ShellId, 1>{shells[0]}, -0.2));
}

TEST_CASE("commit reports the attributed entity for a tolerance-only failure") {
  auto ed = make_square_editor(/*vertex=*/0.1, /*edge=*/0.05, /*face=*/0.2,
                               /*solid=*/0.2);
  const geom::VertexId first_vertex = ed.brep().vertex_ids().front();
  const auto r = ed.commit();
  REQUIRE_FALSE(r.has_value());
  CHECK(r.error().code == geom::Errc::tolerance_violation);
  // The edge's tolerance is below both vertex tolerances; the first offending
  // vertex (v0) owns the finding.
  CHECK(r.error().context_a == first_vertex.index());
}

TEST_CASE("topology errors are reported before tolerance errors") {
  auto ed = make_square_editor(/*vertex=*/0.1, /*edge=*/0.05, /*face=*/0.2,
                               /*solid=*/0.2);
  const geom::EdgeId first_edge = ed.brep().edge_ids().front();
  REQUIRE(ed.erase_edge(first_edge).has_value());
  const auto r = ed.commit();
  REQUIRE_FALSE(r.has_value());
  // Both a dangling edge reference and an edge/vertex tolerance violation are
  // staged; validate_topology() runs first, so invalid_topology wins.
  CHECK(r.error().code == geom::Errc::invalid_topology);
}
