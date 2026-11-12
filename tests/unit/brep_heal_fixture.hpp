#pragma once

// Test fixtures for the healing pipeline. Each fixture builds a small open
// solid whose boundary edges are meant to be sewn by detail::heal_plan and
// healed by heal_brep.

#include <array>
#include <cstdint>
#include <stdexcept>

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

}  // namespace geom::test
