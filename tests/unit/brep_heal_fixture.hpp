#pragma once

// Test fixtures for the healing pipeline. Each fixture builds a small open
// solid whose boundary edges are meant to be sewn by detail::heal_plan and
// healed by heal_brep.

#include <array>
#include <cmath>
#include <cstdint>
#include <span>
#include <stdexcept>
#include <vector>

#include <geom/brep/face_uv.hpp>
#include <geom/brep/solid_editor.hpp>
#include <geom/brep/topology.hpp>
#include <geom/core/result.hpp>
#include <geom/surfaces/plane.hpp>

namespace geom::test {

struct TornSheet {
  Brep brep;
  SolidId solid{};
};

namespace fixture_detail {

template <typename T>
inline void require_ok(const Result<T>& result) {
  if (!result.has_value()) throw std::runtime_error(result.error().message);
}

struct Sheet {
  std::array<VertexId, 4> vertices{};
  std::array<EdgeId, 4> edges{};
  FaceId face{};
};

// Adds a unit square at z whose vertices are always stored in the order
// (0,0), (1,0), (1,1), (0,1). `reverse_edges` stores each edge from its second
// vertex to its first; `reverse_loop` makes the coedge traverse against that
// stored direction, so the loop always walks the square counter-clockwise.
inline Sheet add_sheet(SolidEditor& editor, std::uint32_t surface, double z,
                       bool reverse_edges, bool reverse_loop, double tol) {
  Sheet sheet;
  sheet.vertices = {editor.add_vertex({0.0, 0.0, z}, 0).value(),
                    editor.add_vertex({1.0, 0.0, z}, 0).value(),
                    editor.add_vertex({1.0, 1.0, z}, 0).value(),
                    editor.add_vertex({0.0, 1.0, z}, 0).value()};
  const auto& v = sheet.vertices;
  if (reverse_edges) {
    sheet.edges = {editor.add_edge(v[1], v[0], kNoGeometry, 0, 1, tol).value(),
                   editor.add_edge(v[2], v[1], kNoGeometry, 0, 1, tol).value(),
                   editor.add_edge(v[3], v[2], kNoGeometry, 0, 1, tol).value(),
                   editor.add_edge(v[0], v[3], kNoGeometry, 0, 1, tol).value()};
  } else {
    sheet.edges = {editor.add_edge(v[0], v[1], kNoGeometry, 0, 1, tol).value(),
                   editor.add_edge(v[1], v[2], kNoGeometry, 0, 1, tol).value(),
                   editor.add_edge(v[2], v[3], kNoGeometry, 0, 1, tol).value(),
                   editor.add_edge(v[3], v[0], kNoGeometry, 0, 1, tol).value()};
  }
  sheet.face = editor.add_face(surface, false, tol).value();
  const CoedgeSpec loop[] = {{sheet.edges[0], reverse_loop},
                             {sheet.edges[1], reverse_loop},
                             {sheet.edges[2], reverse_loop},
                             {sheet.edges[3], reverse_loop}};
  require_ok(editor.add_loop(sheet.face, loop, true));
  return sheet;
}

inline TornSheet finish(SolidEditor& editor, std::span<const FaceId> faces,
                        double tol) {
  const auto shell = editor.add_shell(faces, false).value();
  const std::array<ShellId, 1> shells{shell};
  const SolidId solid = editor.add_solid(shells, tol).value();
  return {editor.commit().value(), solid};
}

}  // namespace fixture_detail

// Two 1x1 quads, one at z = 0 and one at z = gap. The quads do not share any
// vertex. By default the second loop is reversed so the two faces traverse
// their common boundary in opposite directions; `same_direction` turns that
// off. Both faces live in one open shell of one solid.
inline TornSheet make_torn_sheet_with_tolerance(double gap, bool same_direction,
                                                double tol) {
  SolidEditor editor;
  const auto below = editor.add_surface(Plane3{{0.0, 0.0, 1.0}, 0.0});
  const auto above = editor.add_surface(Plane3{{0.0, 0.0, 1.0}, gap});
  const auto a = fixture_detail::add_sheet(editor, below, 0.0, false, false, tol);
  const auto b =
      fixture_detail::add_sheet(editor, above, gap, false, !same_direction, tol);
  const std::array<FaceId, 2> faces{a.face, b.face};
  return fixture_detail::finish(editor, faces, tol);
}

inline TornSheet make_torn_sheet(double gap, bool same_direction = false) {
  return make_torn_sheet_with_tolerance(gap, same_direction, 0.0);
}

// Three parallel quads at z = 0, +g and -g. Every edge of the middle quad is
// equidistant from the matching edges of the two outer quads.
inline TornSheet make_ambiguous_rims(double g) {
  SolidEditor editor;
  const auto mid_surface = editor.add_surface(Plane3{{0.0, 0.0, 1.0}, 0.0});
  const auto high_surface = editor.add_surface(Plane3{{0.0, 0.0, 1.0}, g});
  const auto low_surface = editor.add_surface(Plane3{{0.0, 0.0, 1.0}, -g});
  const auto mid = fixture_detail::add_sheet(editor, mid_surface, 0.0, false, false, 0);
  const auto high = fixture_detail::add_sheet(editor, high_surface, g, false, true, 0);
  const auto low = fixture_detail::add_sheet(editor, low_surface, -g, false, false, 0);
  const std::array<FaceId, 3> faces{mid.face, high.face, low.face};
  return fixture_detail::finish(editor, faces, 0.0);
}

// A quad at z = 0 sewn to a square at z = gap that is split into two clockwise
// triangles. The diagonal is an interior edge: it is never sewn, but both of
// its endpoints are welded, so its coedges still need repaired pcurves.
inline TornSheet make_torn_split_sheet(double gap) {
  SolidEditor editor;
  const auto below = editor.add_surface(Plane3{{0.0, 0.0, 1.0}, 0.0});
  const auto above = editor.add_surface(Plane3{{0.0, 0.0, 1.0}, gap});
  const auto a = fixture_detail::add_sheet(editor, below, 0.0, false, false, 0);

  const auto b00 = editor.add_vertex({0.0, 0.0, gap}, 0).value();
  const auto b10 = editor.add_vertex({1.0, 0.0, gap}, 0).value();
  const auto b11 = editor.add_vertex({1.0, 1.0, gap}, 0).value();
  const auto b01 = editor.add_vertex({0.0, 1.0, gap}, 0).value();
  const auto bottom = editor.add_edge(b10, b00, kNoGeometry, 0, 1, 0).value();
  const auto right = editor.add_edge(b11, b10, kNoGeometry, 0, 1, 0).value();
  const auto top = editor.add_edge(b01, b11, kNoGeometry, 0, 1, 0).value();
  const auto left = editor.add_edge(b00, b01, kNoGeometry, 0, 1, 0).value();
  const auto diagonal = editor.add_edge(b00, b11, kNoGeometry, 0, 1, 0).value();

  const auto first = editor.add_face(above, false, 0).value();
  const CoedgeSpec first_loop[] = {{left, false}, {top, false}, {diagonal, true}};
  fixture_detail::require_ok(editor.add_loop(first, first_loop, true));
  const auto second = editor.add_face(above, false, 0).value();
  const CoedgeSpec second_loop[] = {{diagonal, false}, {right, false}, {bottom, false}};
  fixture_detail::require_ok(editor.add_loop(second, second_loop, true));

  const std::array<FaceId, 3> faces{a.face, first, second};
  return fixture_detail::finish(editor, faces, 0.0);
}

// A quad at z = 0 and a quad at z = gap whose (1,0) corner is shifted by
// `offset` in x. Only the two edges away from the shifted corner pair; the
// other two stay farther apart than a 1e-6 max_gap.
inline TornSheet make_partial_torn_sheet(double gap, double offset) {
  SolidEditor editor;
  const auto below = editor.add_surface(Plane3{{0.0, 0.0, 1.0}, 0.0});
  const auto above = editor.add_surface(Plane3{{0.0, 0.0, 1.0}, gap});
  const auto a = fixture_detail::add_sheet(editor, below, 0.0, false, false, 0);

  const auto b00 = editor.add_vertex({0.0, 0.0, gap}, 0).value();
  const auto b10 = editor.add_vertex({1.0 + offset, 0.0, gap}, 0).value();
  const auto b11 = editor.add_vertex({1.0, 1.0, gap}, 0).value();
  const auto b01 = editor.add_vertex({0.0, 1.0, gap}, 0).value();
  const auto b0 = editor.add_edge(b00, b10, kNoGeometry, 0, 1, 0).value();
  const auto b1 = editor.add_edge(b10, b11, kNoGeometry, 0, 1, 0).value();
  const auto b2 = editor.add_edge(b11, b01, kNoGeometry, 0, 1, 0).value();
  const auto b3 = editor.add_edge(b01, b00, kNoGeometry, 0, 1, 0).value();
  const auto face = editor.add_face(above, false, 0).value();
  const CoedgeSpec loop[] = {
      {b0, true}, {b1, true}, {b2, true}, {b3, true}};
  fixture_detail::require_ok(editor.add_loop(face, loop, true));

  const std::array<FaceId, 2> faces{a.face, face};
  return fixture_detail::finish(editor, faces, 0.0);
}

// A quad at z = 0 whose edges and coedges are both stored reversed, sewn to a
// quad at z = gap whose coedges are reversed. Both keep coedges therefore have
// their raw `reversed` flag set, but the geometric traversals are opposite.
inline TornSheet make_reversed_keep_sheet(double gap) {
  SolidEditor editor;
  const auto below = editor.add_surface(Plane3{{0.0, 0.0, 1.0}, 0.0});
  const auto above = editor.add_surface(Plane3{{0.0, 0.0, 1.0}, gap});
  const auto a = fixture_detail::add_sheet(editor, below, 0.0, true, true, 0);
  const auto b = fixture_detail::add_sheet(editor, above, gap, false, true, 0);
  const std::array<FaceId, 2> faces{a.face, b.face};
  return fixture_detail::finish(editor, faces, 0.0);
}

// A 3-edge triangle and a 2-edge part of a second triangle. The nearest
// candidates form A0-B0 and A1-B1; A2's nearest is B1 but B1 prefers A1, so
// A2 is a one-sided preference and must not be sewn.
struct AsymmetricPreference {
  Brep brep;
  SolidId solid{};
  std::array<EdgeId, 3> three{};
  std::array<EdgeId, 2> two{};
};

inline AsymmetricPreference make_asymmetric_preference_sheet(double gap) {
  SolidEditor editor;
  const auto below = editor.add_surface(Plane3{{0.0, 0.0, 1.0}, 0.0});
  const auto above = editor.add_surface(Plane3{{0.0, 0.0, 1.0}, gap});
  const auto a0 = editor.add_vertex({0.229, -0.565, 0.0}, 0).value();
  const auto a1 = editor.add_vertex({-0.854, -0.197, 0.0}, 0).value();
  const auto a2 = editor.add_vertex({-0.169, 0.788, 0.0}, 0).value();
  const auto ae0 = editor.add_edge(a0, a1, kNoGeometry, 0, 1, 0).value();
  const auto ae1 = editor.add_edge(a1, a2, kNoGeometry, 0, 1, 0).value();
  const auto ae2 = editor.add_edge(a2, a0, kNoGeometry, 0, 1, 0).value();
  const auto a_face = editor.add_face(below, false, 0).value();
  const CoedgeSpec a_loop[] = {{ae0, false}, {ae1, false}, {ae2, false}};
  fixture_detail::require_ok(editor.add_loop(a_face, a_loop, true));

  const auto b0 = editor.add_vertex({-0.957, -0.301, gap}, 0).value();
  const auto b1 = editor.add_vertex({-0.474, -0.852, gap}, 0).value();
  const auto b2 = editor.add_vertex({-0.463, 0.742, gap}, 0).value();
  const auto be0 = editor.add_edge(b0, b1, kNoGeometry, 0, 1, 0).value();
  const auto be1 = editor.add_edge(b1, b2, kNoGeometry, 0, 1, 0).value();
  const auto be2 = editor.add_edge(b2, b0, kNoGeometry, 0, 1, 0).value();
  const auto b_face = editor.add_face(above, false, 0).value();
  const CoedgeSpec b_loop[] = {
      {be0, false}, {be1, true}, {be2, true}};
  fixture_detail::require_ok(editor.add_loop(b_face, b_loop, true));

  const std::array<FaceId, 2> faces{a_face, b_face};
  const TornSheet sheet = fixture_detail::finish(editor, faces, 0.0);
  return {sheet.brep, sheet.solid, {ae0, ae1, ae2}, {be0, be1}};
}

// Two triangles that meet only at the shared vertex s. Each returned edge is a
// boundary edge of a different face and both are incident to s, so they are
// not a candidate pair.
struct SharedVertexSheets {
  Brep brep;
  SolidId solid{};
  EdgeId first{};
  EdgeId second{};
};

inline SharedVertexSheets make_shared_vertex_sheets() {
  SolidEditor editor;
  const auto plane = editor.add_surface(Plane3{{0.0, 0.0, 1.0}, 0.0});
  const auto s = editor.add_vertex({0.0, 0.0, 0.0}, 0).value();
  const auto a1 = editor.add_vertex({1.0, 0.0, 0.0}, 0).value();
  const auto a2 = editor.add_vertex({0.0, 1.0, 0.0}, 0).value();
  const auto ae0 = editor.add_edge(s, a1, kNoGeometry, 0, 1, 0).value();
  const auto ae1 = editor.add_edge(a1, a2, kNoGeometry, 0, 1, 0).value();
  const auto ae2 = editor.add_edge(a2, s, kNoGeometry, 0, 1, 0).value();
  const auto a_face = editor.add_face(plane, false, 0).value();
  const CoedgeSpec a_loop[] = {{ae0, false}, {ae1, false}, {ae2, false}};
  fixture_detail::require_ok(editor.add_loop(a_face, a_loop, true));

  const auto b1 = editor.add_vertex({0.8, 0.2, 0.0}, 0).value();
  const auto b2 = editor.add_vertex({0.2, 0.8, 0.0}, 0).value();
  const auto be0 = editor.add_edge(s, b1, kNoGeometry, 0, 1, 0).value();
  const auto be1 = editor.add_edge(b1, b2, kNoGeometry, 0, 1, 0).value();
  const auto be2 = editor.add_edge(b2, s, kNoGeometry, 0, 1, 0).value();
  const auto b_face = editor.add_face(plane, false, 0).value();
  const CoedgeSpec b_loop[] = {{be0, false}, {be1, false}, {be2, false}};
  fixture_detail::require_ok(editor.add_loop(b_face, b_loop, true));

  const std::array<FaceId, 2> faces{a_face, b_face};
  const TornSheet sheet = fixture_detail::finish(editor, faces, 0.0);
  return {sheet.brep, sheet.solid, ae0, be0};
}

// A unit box cut at z = 0.5. The lower half ends at z = 0.5 and the upper half
// starts at the smallest representable z at or below z = 0.5 + gap, so the two
// four-edge rims are boundary loops of one open shell and the seam never exceeds
// gap; the eight rim edges are the only boundary edges. Volume before healing is
// 1 minus the seam.
struct CrackedBox {
  Brep brep;
  SolidId solid{};
};

namespace fixture_detail {

// One coedge of a fixture face, in the order the loop walks it.
struct CoedgeUse {
  EdgeId edge{};
  bool reversed = false;
};

// Adds `plane` as a face with one outer loop over `uses` and a two-point
// pcurve on every coedge. A pcurve runs along its edge from t0 to t1, so the
// stored points are the projections of the edge's start and end vertices; a
// coedge with `reversed` reads them backwards.
inline FaceId add_pcurved_face(SolidEditor& editor, const Plane3& plane, double tol,
                               std::span<const CoedgeUse> uses) {
  const std::uint32_t surface = editor.add_surface(plane);
  const FaceId face = editor.add_face(surface, false, tol).value();
  using Projection = detail::FaceUvProjection;
  const Result<Projection> projection =
      detail::face_uv_projection(editor.brep(), editor.brep().face(face), tol);
  require_ok(projection);

  std::vector<CoedgeSpec> loop;
  loop.reserve(uses.size());
  for (const CoedgeUse& use : uses) {
    const Edge& edge = editor.brep().edge(use.edge);
    const auto project = [&](VertexId id) {
      const Result<Point2d> uv = detail::face_uv_project(
          projection.value(), editor.brep().vertex(id).point, 0.0, false);
      require_ok(uv);
      return uv.value();
    };
    const std::uint32_t pcurve =
        editor.add_pcurve(Pcurve2d{{project(edge.v0), project(edge.v1)}});
    loop.push_back(CoedgeSpec{use.edge, use.reversed, pcurve, 0.0});
  }
  require_ok(editor.add_loop(face, loop, true));
  return face;
}

}  // namespace fixture_detail

// Builds the cracked box described above. `chain_tol` is written to every edge,
// every face and the solid; every vertex keeps tolerance 0. Every face is a
// plane with an outward normal and a loop that is counter-clockwise seen from
// outside, which is what detail::curved_solid_volume reads.
inline CrackedBox make_cracked_box(double gap, double chain_tol = 0.0) {
  SolidEditor editor;

  // Bottom corners, lower rim, upper rim, top corners.
  const VertexId b00 = editor.add_vertex({0.0, 0.0, 0.0}, 0).value();
  const VertexId b10 = editor.add_vertex({1.0, 0.0, 0.0}, 0).value();
  const VertexId b11 = editor.add_vertex({1.0, 1.0, 0.0}, 0).value();
  const VertexId b01 = editor.add_vertex({0.0, 1.0, 0.0}, 0).value();
  const VertexId l00 = editor.add_vertex({0.0, 0.0, 0.5}, 0).value();
  const VertexId l10 = editor.add_vertex({1.0, 0.0, 0.5}, 0).value();
  const VertexId l11 = editor.add_vertex({1.0, 1.0, 0.5}, 0).value();
  const VertexId l01 = editor.add_vertex({0.0, 1.0, 0.5}, 0).value();
  // 0.5 + gap rounds up for many gaps; step down one ulp so the seam never
  // exceeds the requested gap and the documented max_gap contract holds.
  double upper_z = 0.5 + gap;
  if (upper_z - 0.5 > gap) upper_z = std::nextafter(upper_z, 0.0);
  const VertexId u00 = editor.add_vertex({0.0, 0.0, upper_z}, 0).value();
  const VertexId u10 = editor.add_vertex({1.0, 0.0, upper_z}, 0).value();
  const VertexId u11 = editor.add_vertex({1.0, 1.0, upper_z}, 0).value();
  const VertexId u01 = editor.add_vertex({0.0, 1.0, upper_z}, 0).value();
  const VertexId t00 = editor.add_vertex({0.0, 0.0, 1.0}, 0).value();
  const VertexId t10 = editor.add_vertex({1.0, 0.0, 1.0}, 0).value();
  const VertexId t11 = editor.add_vertex({1.0, 1.0, 1.0}, 0).value();
  const VertexId t01 = editor.add_vertex({0.0, 1.0, 1.0}, 0).value();

  // Shared edges are created once and used by both adjacent faces.
  const VertexId b_ring[] = {b00, b10, b11, b01};
  const VertexId l_ring[] = {l00, l10, l11, l01};
  const VertexId u_ring[] = {u00, u10, u11, u01};
  const VertexId t_ring[] = {t00, t10, t11, t01};

  const EdgeId b0 = editor.add_edge(b_ring[0], b_ring[1], kNoGeometry, 0, 1, chain_tol).value();
  const EdgeId b1 = editor.add_edge(b_ring[1], b_ring[2], kNoGeometry, 0, 1, chain_tol).value();
  const EdgeId b2 = editor.add_edge(b_ring[2], b_ring[3], kNoGeometry, 0, 1, chain_tol).value();
  const EdgeId b3 = editor.add_edge(b_ring[3], b_ring[0], kNoGeometry, 0, 1, chain_tol).value();
  const EdgeId t0 = editor.add_edge(t_ring[0], t_ring[1], kNoGeometry, 0, 1, chain_tol).value();
  const EdgeId t1 = editor.add_edge(t_ring[1], t_ring[2], kNoGeometry, 0, 1, chain_tol).value();
  const EdgeId t2 = editor.add_edge(t_ring[2], t_ring[3], kNoGeometry, 0, 1, chain_tol).value();
  const EdgeId t3 = editor.add_edge(t_ring[3], t_ring[0], kNoGeometry, 0, 1, chain_tol).value();
  const EdgeId lv0 = editor.add_edge(b_ring[0], l_ring[0], kNoGeometry, 0, 1, chain_tol).value();
  const EdgeId lv1 = editor.add_edge(b_ring[1], l_ring[1], kNoGeometry, 0, 1, chain_tol).value();
  const EdgeId lv2 = editor.add_edge(b_ring[2], l_ring[2], kNoGeometry, 0, 1, chain_tol).value();
  const EdgeId lv3 = editor.add_edge(b_ring[3], l_ring[3], kNoGeometry, 0, 1, chain_tol).value();
  const EdgeId uv0 = editor.add_edge(u_ring[0], t_ring[0], kNoGeometry, 0, 1, chain_tol).value();
  const EdgeId uv1 = editor.add_edge(u_ring[1], t_ring[1], kNoGeometry, 0, 1, chain_tol).value();
  const EdgeId uv2 = editor.add_edge(u_ring[2], t_ring[2], kNoGeometry, 0, 1, chain_tol).value();
  const EdgeId uv3 = editor.add_edge(u_ring[3], t_ring[3], kNoGeometry, 0, 1, chain_tol).value();
  const EdgeId lr0 = editor.add_edge(l_ring[0], l_ring[1], kNoGeometry, 0, 1, chain_tol).value();
  const EdgeId lr1 = editor.add_edge(l_ring[1], l_ring[2], kNoGeometry, 0, 1, chain_tol).value();
  const EdgeId lr2 = editor.add_edge(l_ring[2], l_ring[3], kNoGeometry, 0, 1, chain_tol).value();
  const EdgeId lr3 = editor.add_edge(l_ring[3], l_ring[0], kNoGeometry, 0, 1, chain_tol).value();
  const EdgeId ur0 = editor.add_edge(u_ring[0], u_ring[1], kNoGeometry, 0, 1, chain_tol).value();
  const EdgeId ur1 = editor.add_edge(u_ring[1], u_ring[2], kNoGeometry, 0, 1, chain_tol).value();
  const EdgeId ur2 = editor.add_edge(u_ring[2], u_ring[3], kNoGeometry, 0, 1, chain_tol).value();
  const EdgeId ur3 = editor.add_edge(u_ring[3], u_ring[0], kNoGeometry, 0, 1, chain_tol).value();

  // Each loop below is counter-clockwise seen from outside the box. The rim
  // edges (lr*, ur*) are used once each and are therefore the boundary.
  const std::array<fixture_detail::CoedgeUse, 4> bottom_loop{{
      {b3, true}, {b2, true}, {b1, true}, {b0, true}}};
  const std::array<fixture_detail::CoedgeUse, 4> top_loop{{
      {t0, false}, {t1, false}, {t2, false}, {t3, false}}};
  const std::array<fixture_detail::CoedgeUse, 4> lower_south_loop{{
      {b0, false}, {lv1, false}, {lr0, true}, {lv0, true}}};
  const std::array<fixture_detail::CoedgeUse, 4> lower_east_loop{{
      {b1, false}, {lv2, false}, {lr1, true}, {lv1, true}}};
  const std::array<fixture_detail::CoedgeUse, 4> lower_north_loop{{
      {b2, false}, {lv3, false}, {lr2, true}, {lv2, true}}};
  const std::array<fixture_detail::CoedgeUse, 4> lower_west_loop{{
      {b3, false}, {lv0, false}, {lr3, true}, {lv3, true}}};
  const std::array<fixture_detail::CoedgeUse, 4> upper_south_loop{{
      {ur0, false}, {uv1, false}, {t0, true}, {uv0, true}}};
  const std::array<fixture_detail::CoedgeUse, 4> upper_east_loop{{
      {ur1, false}, {uv2, false}, {t1, true}, {uv1, true}}};
  const std::array<fixture_detail::CoedgeUse, 4> upper_north_loop{{
      {ur2, false}, {uv3, false}, {t2, true}, {uv2, true}}};
  const std::array<fixture_detail::CoedgeUse, 4> upper_west_loop{{
      {ur3, false}, {uv0, false}, {t3, true}, {uv3, true}}};

  const FaceId bottom = fixture_detail::add_pcurved_face(
      editor, Plane3{{0.0, 0.0, -1.0}, 0.0}, chain_tol, bottom_loop);
  const FaceId top = fixture_detail::add_pcurved_face(
      editor, Plane3{{0.0, 0.0, 1.0}, 1.0}, chain_tol, top_loop);
  const FaceId lower_south = fixture_detail::add_pcurved_face(
      editor, Plane3{{0.0, -1.0, 0.0}, 0.0}, chain_tol, lower_south_loop);
  const FaceId lower_east = fixture_detail::add_pcurved_face(
      editor, Plane3{{1.0, 0.0, 0.0}, 1.0}, chain_tol, lower_east_loop);
  const FaceId lower_north = fixture_detail::add_pcurved_face(
      editor, Plane3{{0.0, 1.0, 0.0}, 1.0}, chain_tol, lower_north_loop);
  const FaceId lower_west = fixture_detail::add_pcurved_face(
      editor, Plane3{{-1.0, 0.0, 0.0}, 0.0}, chain_tol, lower_west_loop);
  const FaceId upper_south = fixture_detail::add_pcurved_face(
      editor, Plane3{{0.0, -1.0, 0.0}, 0.0}, chain_tol, upper_south_loop);
  const FaceId upper_east = fixture_detail::add_pcurved_face(
      editor, Plane3{{1.0, 0.0, 0.0}, 1.0}, chain_tol, upper_east_loop);
  const FaceId upper_north = fixture_detail::add_pcurved_face(
      editor, Plane3{{0.0, 1.0, 0.0}, 1.0}, chain_tol, upper_north_loop);
  const FaceId upper_west = fixture_detail::add_pcurved_face(
      editor, Plane3{{-1.0, 0.0, 0.0}, 0.0}, chain_tol, upper_west_loop);

  // Both components live in one open shell, so the shell has Euler
  // characteristic 2 before healing.
  const std::array<FaceId, 10> faces{bottom, top,      lower_south, lower_east,
                                     lower_north, lower_west, upper_south, upper_east,
                                     upper_north, upper_west};
  const ShellId shell = editor.add_shell(faces, false).value();
  const std::array<ShellId, 1> shells{shell};
  const SolidId solid = editor.add_solid(shells, chain_tol).value();
  return {editor.commit().value(), solid};
}

}  // namespace geom::test
