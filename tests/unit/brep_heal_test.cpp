#include <algorithm>
#include <array>
#include <charconv>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <iomanip>
#include <limits>
#include <set>
#include <sstream>
#include <string>
#include <string_view>
#include <variant>
#include <vector>

#include <doctest/doctest.h>
#include <geom/brep/detail/curved_volume.hpp>
#include <geom/curves/analytic_curve.hpp>
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

// Appends a value in a form that round-trips exactly, so fingerprints compare
// bit-for-bit rather than at a decimal display precision.
template <typename T>
void append_integer(std::string& out, T value) {
  char buffer[64];
  const auto result = std::to_chars(buffer, buffer + sizeof buffer, value);
  out.append(buffer, result.ptr);
}

void append_double(std::string& out, double value) {
  char buffer[64];
  const auto result = std::to_chars(buffer, buffer + sizeof buffer, value,
                                    std::chars_format::hex);
  out.append(buffer, result.ptr);
}

template <typename Id>
void append_handle(std::string& out, Id id) {
  append_integer(out, id.index());
  out.push_back('.');
  if (id.is_null()) {
    out.push_back('-');
  } else {
    append_integer(out, static_cast<unsigned>(id.generation()));
  }
}

// A stable, exact rendering of the whole Brep: every table slot in index order
// (including dead slots), alive flags and generations, all geometry/tolerance
// fields, loop coedge lists, shell face lists, solid shell lists, and the
// append-only geometry store. Comparing this before and after pins that a
// rejected heal never mutates its input and that an idempotent heal changes
// nothing.
std::string heal_fingerprint(const Brep& brep) {
  std::string out;
  out.reserve(1u << 16);

  const auto generations = [](const auto& ids, std::size_t slots) {
    std::vector<int> result(slots, -1);
    for (const auto id : ids) {
      result[id.index()] = static_cast<int>(id.generation());
    }
    return result;
  };

  const std::vector<int> vertex_generations =
      generations(brep.vertex_ids(), brep.vertices().size());
  const std::vector<int> edge_generations =
      generations(brep.edge_ids(), brep.edges().size());
  const std::vector<int> coedge_generations =
      generations(brep.coedge_ids(), brep.coedges().size());
  const std::vector<int> loop_generations =
      generations(brep.loop_ids(), brep.loops().size());
  const std::vector<int> face_generations =
      generations(brep.face_ids(), brep.faces().size());
  const std::vector<int> shell_generations =
      generations(brep.shell_ids(), brep.shells().size());
  const std::vector<int> solid_generations =
      generations(brep.solid_ids(), brep.solids().size());

  out += "counts ";
  append_integer(out, brep.vertex_count());
  out.push_back(' ');
  append_integer(out, brep.edge_count());
  out.push_back(' ');
  append_integer(out, brep.coedge_count());
  out.push_back(' ');
  append_integer(out, brep.loop_count());
  out.push_back(' ');
  append_integer(out, brep.face_count());
  out.push_back(' ');
  append_integer(out, brep.shell_count());
  out.push_back(' ');
  append_integer(out, brep.solid_count());
  out.push_back('\n');

  const auto& vertices = brep.vertices();
  for (std::size_t index = 0; index < vertices.size(); ++index) {
    const Vertex& vertex = vertices[index];
    out += "v ";
    append_integer(out, index);
    out.push_back(' ');
    append_integer(out, vertex_generations[index]);
    out.push_back(' ');
    out.push_back(vertex_generations[index] >= 0 ? '1' : '0');
    out.push_back(' ');
    append_double(out, vertex.point[0]);
    out.push_back(' ');
    append_double(out, vertex.point[1]);
    out.push_back(' ');
    append_double(out, vertex.point[2]);
    out.push_back(' ');
    append_double(out, vertex.tol);
    out.push_back('\n');
  }

  const auto& edges = brep.edges();
  for (std::size_t index = 0; index < edges.size(); ++index) {
    const Edge& edge = edges[index];
    out += "e ";
    append_integer(out, index);
    out.push_back(' ');
    append_integer(out, edge_generations[index]);
    out.push_back(' ');
    out.push_back(edge_generations[index] >= 0 ? '1' : '0');
    out.push_back(' ');
    append_handle(out, edge.v0);
    out.push_back(' ');
    append_handle(out, edge.v1);
    out.push_back(' ');
    append_integer(out, edge.curve);
    out.push_back(' ');
    append_double(out, edge.t0);
    out.push_back(' ');
    append_double(out, edge.t1);
    out.push_back(' ');
    append_double(out, edge.tol);
    out.push_back('\n');
  }

  const auto& coedges = brep.coedges();
  for (std::size_t index = 0; index < coedges.size(); ++index) {
    const Coedge& coedge = coedges[index];
    out += "c ";
    append_integer(out, index);
    out.push_back(' ');
    append_integer(out, coedge_generations[index]);
    out.push_back(' ');
    out.push_back(coedge_generations[index] >= 0 ? '1' : '0');
    out.push_back(' ');
    append_handle(out, coedge.edge);
    out.push_back(' ');
    append_handle(out, coedge.face);
    out.push_back(' ');
    append_handle(out, coedge.loop);
    out.push_back(' ');
    append_handle(out, coedge.next);
    out.push_back(' ');
    append_handle(out, coedge.prev);
    out.push_back(' ');
    out.push_back(coedge.reversed ? '1' : '0');
    out.push_back(' ');
    append_integer(out, coedge.pcurve);
    out.push_back(' ');
    append_double(out, coedge.tol);
    out.push_back('\n');
  }

  const auto& loops = brep.loops();
  for (std::size_t index = 0; index < loops.size(); ++index) {
    const Loop& loop = loops[index];
    const bool alive = loop_generations[index] >= 0;
    out += "l ";
    append_integer(out, index);
    out.push_back(' ');
    append_integer(out, loop_generations[index]);
    out.push_back(' ');
    out.push_back(alive ? '1' : '0');
    out.push_back(' ');
    append_handle(out, loop.face);
    out.push_back(' ');
    append_handle(out, loop.start);
    out.push_back(' ');
    append_handle(out, loop.next);
    out.push_back(' ');
    out.push_back(loop.outer ? '1' : '0');
    if (alive && !loop.start.is_null() && brep.coedge_alive(loop.start)) {
      CoedgeId current = loop.start;
      std::size_t seen = 0;
      out.push_back(' ');
      do {
        out += "use ";
        append_handle(out, current);
        out.push_back(' ');
        if (!brep.coedge_alive(brep.coedge(current).next)) break;
        current = brep.coedge(current).next;
        ++seen;
      } while (current != loop.start && seen <= brep.coedge_count());
    }
    out.push_back('\n');
  }

  const auto& faces = brep.faces();
  for (std::size_t index = 0; index < faces.size(); ++index) {
    const Face& face = faces[index];
    out += "f ";
    append_integer(out, index);
    out.push_back(' ');
    append_integer(out, face_generations[index]);
    out.push_back(' ');
    out.push_back(face_generations[index] >= 0 ? '1' : '0');
    out.push_back(' ');
    append_integer(out, face.surface);
    out.push_back(' ');
    append_handle(out, face.first_loop);
    out.push_back(' ');
    out.push_back(face.reversed ? '1' : '0');
    out.push_back(' ');
    append_double(out, face.tol);
    out.push_back('\n');
  }

  const auto& shells = brep.shells();
  for (std::size_t index = 0; index < shells.size(); ++index) {
    const Shell& shell = shells[index];
    const bool alive = shell_generations[index] >= 0;
    out += "s ";
    append_integer(out, index);
    out.push_back(' ');
    append_integer(out, shell_generations[index]);
    out.push_back(' ');
    out.push_back(alive ? '1' : '0');
    out.push_back(' ');
    append_integer(out, shell.first_face);
    out.push_back(' ');
    append_integer(out, shell.face_count);
    out.push_back(' ');
    out.push_back(shell.closed ? '1' : '0');
    if (alive) {
      const ShellId id =
          ShellId::make(static_cast<std::uint32_t>(index),
                        static_cast<std::uint8_t>(shell_generations[index]));
      for (const FaceId face : brep.shell_faces(id)) {
        out.push_back(' ');
        append_handle(out, face);
      }
    }
    out.push_back('\n');
  }

  const auto& solids = brep.solids();
  for (std::size_t index = 0; index < solids.size(); ++index) {
    const SolidEntity& solid = solids[index];
    const bool alive = solid_generations[index] >= 0;
    out += "d ";
    append_integer(out, index);
    out.push_back(' ');
    append_integer(out, solid_generations[index]);
    out.push_back(' ');
    out.push_back(alive ? '1' : '0');
    out.push_back(' ');
    append_integer(out, solid.first_shell);
    out.push_back(' ');
    append_integer(out, solid.shell_count);
    out.push_back(' ');
    append_double(out, solid.tol);
    if (alive) {
      const SolidId id =
          SolidId::make(static_cast<std::uint32_t>(index),
                        static_cast<std::uint8_t>(solid_generations[index]));
      for (const ShellId shell : brep.solid_shells(id)) {
        out.push_back(' ');
        append_handle(out, shell);
      }
    }
    out.push_back('\n');
  }

  const GeometryStore& geometry = brep.geometry();
  out += "geometry ";
  append_integer(out, geometry.curves.size());
  out.push_back(' ');
  append_integer(out, geometry.surfaces.size());
  out.push_back(' ');
  append_integer(out, geometry.pcurves.size());
  out.push_back('\n');
  for (const Pcurve2d& pcurve : geometry.pcurves) {
    out += "pcurve";
    for (const Point2d& point : pcurve.points) {
      out.push_back(' ');
      append_double(out, point[0]);
      out.push_back(' ');
      append_double(out, point[1]);
    }
    out.push_back('\n');
  }
  return out;
}

template <typename Id>
Id find_live(const std::vector<Id>& ids, std::uint32_t index) {
  for (const Id id : ids) {
    if (id.index() == index) return id;
  }
  return Id{};
}

// Fingerprints only the entities reachable from one solid. Used to prove that
// healing a target solid leaves a second solid in the same Brep untouched.
std::string solid_fingerprint(const Brep& brep, SolidId solid) {
  std::set<std::uint32_t> shell_ids;
  std::set<std::uint32_t> face_ids;
  std::set<std::uint32_t> loop_ids;
  std::set<std::uint32_t> coedge_ids;
  std::set<std::uint32_t> edge_ids;
  std::set<std::uint32_t> vertex_ids;

  for (const ShellId shell : brep.solid_shells(solid)) {
    shell_ids.insert(shell.index());
    for (const FaceId face : brep.shell_faces(shell)) {
      face_ids.insert(face.index());
      LoopId loop_id = brep.face(face).first_loop;
      std::size_t loops_seen = 0;
      while (brep.loop_alive(loop_id) && loops_seen <= brep.loop_count()) {
        const Loop& loop = brep.loop(loop_id);
        loop_ids.insert(loop_id.index());
        if (brep.coedge_alive(loop.start)) {
          CoedgeId current = loop.start;
          std::size_t coedges_seen = 0;
          while (brep.coedge_alive(current) &&
                 coedges_seen <= brep.coedge_count()) {
            const Coedge& coedge = brep.coedge(current);
            coedge_ids.insert(current.index());
            if (brep.edge_alive(coedge.edge)) {
              edge_ids.insert(coedge.edge.index());
              const Edge& edge = brep.edge(coedge.edge);
              if (brep.vertex_alive(edge.v0)) vertex_ids.insert(edge.v0.index());
              if (brep.vertex_alive(edge.v1)) vertex_ids.insert(edge.v1.index());
            }
            if (coedge.next == loop.start) break;
            current = coedge.next;
            ++coedges_seen;
          }
        }
        loop_id = loop.next;
        ++loops_seen;
      }
    }
  }

  std::string out;
  out += "solid ";
  append_handle(out, solid);
  out.push_back(' ');
  append_double(out, brep.solid(solid).tol);
  out.push_back('\n');

  const std::vector<ShellId> live_shells = brep.shell_ids();
  for (const std::uint32_t index : shell_ids) {
    const ShellId id = find_live(live_shells, index);
    if (id.is_null()) continue;
    const Shell& shell = brep.shell(id);
    out += "shell ";
    append_handle(out, id);
    out.push_back(' ');
    out.push_back(shell.closed ? '1' : '0');
    for (const FaceId face : brep.shell_faces(id)) {
      out.push_back(' ');
      append_handle(out, face);
    }
    out.push_back('\n');
  }

  const std::vector<FaceId> live_faces = brep.face_ids();
  for (const std::uint32_t index : face_ids) {
    const FaceId id = find_live(live_faces, index);
    if (id.is_null()) continue;
    const Face& face = brep.face(id);
    out += "face ";
    append_handle(out, id);
    out.push_back(' ');
    append_integer(out, face.surface);
    out.push_back(' ');
    append_handle(out, face.first_loop);
    out.push_back(' ');
    out.push_back(face.reversed ? '1' : '0');
    out.push_back(' ');
    append_double(out, face.tol);
    out.push_back('\n');
  }

  const std::vector<LoopId> live_loops = brep.loop_ids();
  for (const std::uint32_t index : loop_ids) {
    const LoopId id = find_live(live_loops, index);
    if (id.is_null()) continue;
    const Loop& loop = brep.loop(id);
    out += "loop ";
    append_handle(out, id);
    out.push_back(' ');
    append_handle(out, loop.face);
    out.push_back(' ');
    out.push_back(loop.outer ? '1' : '0');
    if (brep.coedge_alive(loop.start)) {
      CoedgeId current = loop.start;
      std::size_t seen = 0;
      do {
        out.push_back(' ');
        append_handle(out, current);
        current = brep.coedge(current).next;
        ++seen;
      } while (current != loop.start && brep.coedge_alive(current) &&
               seen <= brep.coedge_count());
    }
    out.push_back('\n');
  }

  const std::vector<CoedgeId> live_coedges = brep.coedge_ids();
  for (const std::uint32_t index : coedge_ids) {
    const CoedgeId id = find_live(live_coedges, index);
    if (id.is_null()) continue;
    const Coedge& coedge = brep.coedge(id);
    out += "coedge ";
    append_handle(out, id);
    out.push_back(' ');
    append_handle(out, coedge.edge);
    out.push_back(' ');
    append_handle(out, coedge.face);
    out.push_back(' ');
    out.push_back(coedge.reversed ? '1' : '0');
    out.push_back(' ');
    append_integer(out, coedge.pcurve);
    out.push_back(' ');
    append_double(out, coedge.tol);
    out.push_back('\n');
  }

  const std::vector<EdgeId> live_edges = brep.edge_ids();
  for (const std::uint32_t index : edge_ids) {
    const EdgeId id = find_live(live_edges, index);
    if (id.is_null()) continue;
    const Edge& edge = brep.edge(id);
    out += "edge ";
    append_handle(out, id);
    out.push_back(' ');
    append_handle(out, edge.v0);
    out.push_back(' ');
    append_handle(out, edge.v1);
    out.push_back(' ');
    append_integer(out, edge.curve);
    out.push_back(' ');
    append_double(out, edge.t0);
    out.push_back(' ');
    append_double(out, edge.t1);
    out.push_back(' ');
    append_double(out, edge.tol);
    out.push_back('\n');
  }

  const std::vector<VertexId> live_vertices = brep.vertex_ids();
  for (const std::uint32_t index : vertex_ids) {
    const VertexId id = find_live(live_vertices, index);
    if (id.is_null()) continue;
    const Vertex& vertex = brep.vertex(id);
    out += "vertex ";
    append_handle(out, id);
    out.push_back(' ');
    append_double(out, vertex.point[0]);
    out.push_back(' ');
    append_double(out, vertex.point[1]);
    out.push_back(' ');
    append_double(out, vertex.point[2]);
    out.push_back(' ');
    append_double(out, vertex.tol);
    out.push_back('\n');
  }
  return out;
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
  const std::string fingerprint = heal_fingerprint(cracked.brep);

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

  CHECK(heal_fingerprint(cracked.brep) == fingerprint);
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
  const std::string fingerprint = heal_fingerprint(cracked.brep);

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
  CHECK(heal_fingerprint(result.brep) == fingerprint);
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
  const std::array<EdgeId, 6> boundary{
      fixture.three[0], fixture.three[1], fixture.three[2],
      fixture.b_edges[0], fixture.b_edges[1], fixture.b_edges[2]};

  const auto plan_result = detail::heal_plan(fixture.brep, boundary, 1.0);
  REQUIRE(plan_result.has_value());
  const detail::HealPlan& plan = plan_result.value();
  CHECK(plan.issues.empty());
  REQUIRE(plan.pairs.size() == 2u);
  // A2's nearest candidate is B1 and their welded traversal directions agree,
  // but B1 prefers A1. Only mutual preferences may be sewn.
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

TEST_CASE("heal: QA-12 rejects a 1e-3 cracked box at its measured gap") {
  const test::CrackedBox cracked = test::make_cracked_box(1e-3);
  const std::string fingerprint = heal_fingerprint(cracked.brep);
  const std::vector<EdgeId> boundary =
      detail::heal_boundary_edges(cracked.brep, cracked.solid);
  REQUIRE(boundary.size() == 8u);

  const auto healed =
      heal_brep(cracked.brep, cracked.solid,
                HealOptions{.tolerance = 1e-9, .max_gap = 1e-6,
                            .max_tolerance = 1e-6});
  REQUIRE(healed.has_value());
  const HealResult& result = healed.value();
  CHECK_FALSE(result.report.changed);
  CHECK(heal_fingerprint(result.brep) == fingerprint);
  REQUIRE(result.report.issues.size() == 1u);
  const HealIssue& issue = result.report.issues.front();
  CHECK(issue.kind == HealIssueKind::gap_too_large);
  CHECK(issue.gap == doctest::Approx(1e-3).epsilon(1e-6));
  CHECK(issue.entity_a != kNoEntity);
  CHECK(issue.entity_b > issue.entity_a);

  const auto in_boundary = [&](std::uint32_t index) {
    return std::any_of(boundary.begin(), boundary.end(),
                       [&](const EdgeId edge) { return edge.index() == index; });
  };
  CHECK(in_boundary(issue.entity_a));
  CHECK(in_boundary(issue.entity_b));

  EdgeId edge_a{};
  EdgeId edge_b{};
  for (const EdgeId edge : boundary) {
    if (edge.index() == issue.entity_a) edge_a = edge;
    if (edge.index() == issue.entity_b) edge_b = edge;
  }
  REQUIRE(!edge_a.is_null());
  REQUIRE(!edge_b.is_null());
  const Edge& a = cracked.brep.edge(edge_a);
  const Edge& b = cracked.brep.edge(edge_b);
  const Point3d& a0 = cracked.brep.vertex(a.v0).point;
  const Point3d& a1 = cracked.brep.vertex(a.v1).point;
  const Point3d& b0 = cracked.brep.vertex(b.v0).point;
  const Point3d& b1 = cracked.brep.vertex(b.v1).point;
  const double d0 = norm2(a0 - b0);
  const double d1 = norm2(a1 - b1);
  const Point3d expected =
      d0 <= d1 ? Point3d{(a0[0] + b0[0]) * 0.5, (a0[1] + b0[1]) * 0.5,
                         (a0[2] + b0[2]) * 0.5}
               : Point3d{(a1[0] + b1[0]) * 0.5, (a1[1] + b1[1]) * 0.5,
                         (a1[2] + b1[2]) * 0.5};
  CHECK(norm2(issue.location - expected) <= 1e-6);
}

TEST_CASE("heal: healing a closed result again is idempotent") {
  const test::CrackedBox cracked = test::make_cracked_box(1e-6);
  const auto first = heal_brep(cracked.brep, cracked.solid, {});
  REQUIRE(first.has_value());
  const std::string fingerprint = heal_fingerprint(first.value().brep);

  const auto second = heal_brep(first.value().brep, cracked.solid, {});
  REQUIRE(second.has_value());
  CHECK_FALSE(second.value().report.changed);
  CHECK(second.value().report.issues.empty());
  CHECK(second.value().report.welded_vertices == 0u);
  CHECK(second.value().report.sewn_edges == 0u);
  // Nothing is raised on a healed input, so the field stays zero even though
  // the chain keeps its 1e-6 tolerances.
  CHECK(second.value().report.max_tolerance == 0.0);
  CHECK(heal_fingerprint(second.value().brep) == fingerprint);
}

TEST_CASE("heal: an ambiguous rim is reported instead of sewn") {
  const test::TornSheet rims = test::make_ambiguous_rims(0.9e-6);
  const std::string fingerprint = heal_fingerprint(rims.brep);
  HealOptions options;
  options.max_gap = 1e-6;

  const auto healed = heal_brep(rims.brep, rims.solid, options);
  REQUIRE(healed.has_value());
  const HealResult& result = healed.value();
  CHECK_FALSE(result.report.changed);
  CHECK(result.report.sewn_edges == 0u);
  CHECK(heal_fingerprint(result.brep) == fingerprint);
  REQUIRE(result.report.issues.size() == 4u);
  for (const HealIssue& issue : result.report.issues) {
    CHECK(issue.kind == HealIssueKind::ambiguous_gap);
  }
}

TEST_CASE("heal: opposite traversal orientation is reported") {
  const test::TornSheet sheet = test::make_torn_sheet(1e-6, true);
  const std::string fingerprint = heal_fingerprint(sheet.brep);

  const auto healed = heal_brep(sheet.brep, sheet.solid, {});
  REQUIRE(healed.has_value());
  const HealResult& result = healed.value();
  CHECK_FALSE(result.report.changed);
  CHECK(result.report.sewn_edges == 0u);
  CHECK(heal_fingerprint(result.brep) == fingerprint);
  REQUIRE(!result.report.issues.empty());
  for (const HealIssue& issue : result.report.issues) {
    CHECK(issue.kind == HealIssueKind::inconsistent_orientation);
  }
}

TEST_CASE("heal: a tolerance cap rejects a 1e-6 seam") {
  const test::CrackedBox cracked = test::make_cracked_box(1e-6);
  const std::string fingerprint = heal_fingerprint(cracked.brep);
  HealOptions options;
  options.max_tolerance = 1e-7;

  const auto healed = heal_brep(cracked.brep, cracked.solid, options);
  REQUIRE(healed.has_value());
  const HealResult& result = healed.value();
  CHECK_FALSE(result.report.changed);
  CHECK(result.report.welded_vertices == 0u);
  CHECK(result.report.sewn_edges == 0u);
  CHECK(heal_fingerprint(result.brep) == fingerprint);
  REQUIRE(result.report.issues.size() == 1u);
  CHECK(result.report.issues.front().kind == HealIssueKind::invalid_result);
  CHECK(std::string_view(result.report.issues.front().message) ==
        "heal_brep: tolerance cap exceeded");
}

TEST_CASE("heal: the cap pre-check also rejects an already-tall chain") {
  const test::CrackedBox cracked = test::make_cracked_box(1e-6, 1e-6);
  const std::string fingerprint = heal_fingerprint(cracked.brep);
  HealOptions options;
  options.max_gap = 1e-6;
  options.max_tolerance = 1e-7;

  const auto healed = heal_brep(cracked.brep, cracked.solid, options);
  REQUIRE(healed.has_value());
  const HealResult& result = healed.value();
  CHECK_FALSE(result.report.changed);
  CHECK(heal_fingerprint(result.brep) == fingerprint);
  REQUIRE(result.report.issues.size() == 1u);
  CHECK(result.report.issues.front().kind == HealIssueKind::invalid_result);
  CHECK(std::string_view(result.report.issues.front().message) ==
        "heal_brep: tolerance cap exceeded");
  for (const VertexId id : result.brep.vertex_ids()) {
    CHECK(result.brep.vertex(id).tol <= 1e-7);
  }
}

TEST_CASE("heal: invalid options are tolerance violations") {
  const test::TornSheet sheet = test::make_torn_sheet(1e-6);
  const auto options = [](double tolerance, double max_gap, double max_tolerance) {
    HealOptions value;
    value.tolerance = tolerance;
    value.max_gap = max_gap;
    value.max_tolerance = max_tolerance;
    return value;
  };
  const double nan = std::numeric_limits<double>::quiet_NaN();
  const double inf = std::numeric_limits<double>::infinity();
  const std::vector<std::pair<HealOptions, std::string_view>> cases{
      {options(-1.0, 1e-6, 1e-6),
       "heal_brep: tolerance must be finite and non-negative"},
      {options(0.0, 0.0, 1e-6),
       "heal_brep: max_gap must be finite and positive"},
      {options(0.0, -1.0, 1e-6),
       "heal_brep: max_gap must be finite and positive"},
      {options(0.0, nan, 1e-6),
       "heal_brep: max_gap must be finite and positive"},
      {options(0.0, inf, 1e-6),
       "heal_brep: max_gap must be finite and positive"},
      {options(0.0, 1e-6, 0.0),
       "heal_brep: max_tolerance must be finite and positive"},
      {options(0.0, 1e-6, -1.0),
       "heal_brep: max_tolerance must be finite and positive"},
      {options(0.0, 1e-6, nan),
       "heal_brep: max_tolerance must be finite and positive"},
      {options(0.0, 1e-6, inf),
       "heal_brep: max_tolerance must be finite and positive"}};
  for (const auto& [bad, message] : cases) {
    const auto healed = heal_brep(sheet.brep, sheet.solid, bad);
    REQUIRE_FALSE(healed.has_value());
    CHECK(healed.error().code == Errc::tolerance_violation);
    CHECK(std::string_view(healed.error().message) == message);
  }

  const SolidId stale = SolidId::make(
      sheet.solid.index(),
      static_cast<std::uint8_t>(sheet.solid.generation() + 1u));
  for (const SolidId solid : {SolidId{}, stale}) {
    const auto healed = heal_brep(sheet.brep, solid, {});
    REQUIRE_FALSE(healed.has_value());
    CHECK(healed.error().code == Errc::index_out_of_range);
    CHECK(std::string_view(healed.error().message) == "heal_brep: unknown solid");
  }
}

TEST_CASE("heal: a closed box is unchanged and has no counters") {
  const Brep box = make_box({0, 0, 0}, {1, 1, 1}, 0).value();
  const SolidId solid = box.solid_ids().front();
  const std::string fingerprint = heal_fingerprint(box);

  const auto healed = heal_brep(box, solid, {});
  REQUIRE(healed.has_value());
  const HealResult& result = healed.value();
  CHECK_FALSE(result.report.changed);
  CHECK(result.report.boundary_edges == 0u);
  CHECK(result.report.welded_vertices == 0u);
  CHECK(result.report.sewn_edges == 0u);
  CHECK(result.report.max_closed_gap == 0.0);
  CHECK(result.report.max_tolerance == 0.0);
  CHECK(result.report.issues.empty());
  CHECK(heal_fingerprint(result.brep) == fingerprint);
}

TEST_CASE("heal: a curved operand is rejected and left untouched") {
  const Brep cylinder =
      make_cylinder({0.0, 0.0, 0.0}, {0.0, 0.0, 1.0}, 1.0, 1.0, 0.0).value();
  const std::string fingerprint = heal_fingerprint(cylinder);

  const auto healed = heal_brep(cylinder, cylinder.solid_ids().front(), {});
  REQUIRE_FALSE(healed.has_value());
  CHECK(healed.error().code == Errc::unsupported);
  CHECK(std::string_view(healed.error().message) ==
        "heal_brep: curved geometry is not supported yet");
  CHECK(heal_fingerprint(cylinder) == fingerprint);
}

TEST_CASE("heal: a curved incident edge is rejected even with planar faces") {
  const test::CrackedBox cracked = test::make_cracked_box(1e-6);
  std::vector<EdgeId> boundary =
      detail::heal_boundary_edges(cracked.brep, cracked.solid);
  REQUIRE(!boundary.empty());

  SolidEditor editor(cracked.brep);
  const std::uint32_t curve = editor.add_curve(CircleCurve3{});
  REQUIRE(editor.set_edge_geometry(boundary.front(), curve, 0.0, 1.0).has_value());
  const Brep brep = editor.commit().value();

  const auto healed = heal_brep(brep, cracked.solid, {});
  REQUIRE_FALSE(healed.has_value());
  CHECK(healed.error().code == Errc::unsupported);
  CHECK(std::string_view(healed.error().message) ==
        "heal_brep: curved geometry is not supported yet");
}

TEST_CASE("heal: a second solid is not touched") {
  const test::TornSheet sheet = test::make_torn_sheet(1e-6);
  SolidEditor editor(sheet.brep);
  const std::uint32_t surface = editor.add_surface(Plane3{{0, 0, 1}, 0});
  const std::array<EdgeId, 4> edges = add_quad_edges(editor, 5.0, 0.0, 1.0);
  const SolidId other = add_open_sheet(editor, surface, edges);
  const Brep brep = editor.commit().value();
  const std::string other_fingerprint = solid_fingerprint(brep, other);

  CHECK(detail::heal_boundary_edges(brep, sheet.solid).size() == 8u);
  CHECK(detail::heal_boundary_edges(brep, other).size() == 4u);
  const auto healed = heal_brep(brep, sheet.solid, {});
  REQUIRE(healed.has_value());
  CHECK(healed.value().report.changed);
  CHECK(detail::heal_boundary_edges(healed.value().brep, sheet.solid).empty());
  CHECK(detail::heal_boundary_edges(healed.value().brep, other).size() == 4u);
  CHECK(solid_fingerprint(healed.value().brep, other) == other_fingerprint);
}

TEST_CASE("heal: a triangle with no mutual pair is rejected") {
  SolidEditor editor;
  const std::uint32_t surface = editor.add_surface(Plane3{{0, 0, 1}, 0});
  const auto v0 = editor.add_vertex({0, 0, 0}, 0).value();
  const auto v1 = editor.add_vertex({1, 0, 0}, 0).value();
  const auto v2 = editor.add_vertex({0, 1, 0}, 0).value();
  const auto e0 = editor.add_edge(v0, v1, kNoGeometry, 0, 1, 0).value();
  const auto e1 = editor.add_edge(v1, v2, kNoGeometry, 0, 1, 0).value();
  const auto e2 = editor.add_edge(v2, v0, kNoGeometry, 0, 1, 0).value();
  const auto face = editor.add_face(surface, false, 0).value();
  const CoedgeSpec loop[] = {{e0, false}, {e1, false}, {e2, false}};
  REQUIRE(editor.add_loop(face, loop, true).has_value());
  const auto shell = editor.add_shell(std::array<FaceId, 1>{face}, false).value();
  const SolidId solid =
      editor.add_solid(std::array<ShellId, 1>{shell}, 0).value();
  const Brep brep = editor.commit().value();
  const std::string fingerprint = heal_fingerprint(brep);

  const auto healed = heal_brep(brep, solid, {});
  REQUIRE(healed.has_value());
  const HealResult& result = healed.value();
  CHECK_FALSE(result.report.changed);
  CHECK(heal_fingerprint(result.brep) == fingerprint);
  REQUIRE(result.report.issues.size() == 1u);
  CHECK(result.report.issues.front().kind == HealIssueKind::invalid_result);
  CHECK(std::string_view(result.report.issues.front().message) ==
        "heal_brep: no mutual pair");
}

TEST_CASE("heal: the direction predicate follows weld-mapped starts") {
  SolidEditor editor;
  const std::uint32_t surface = editor.add_surface(Plane3{{0, 0, 1}, 0});
  const auto k0 = editor.add_vertex({0, 0, 0}, 0).value();
  const auto k1 = editor.add_vertex({1, 0, 0}, 0).value();
  const auto d0 = editor.add_vertex({0, 1, 0}, 0).value();
  const auto d1 = editor.add_vertex({1, 1, 0}, 0).value();
  const auto keep = editor.add_edge(k0, k1, kNoGeometry, 0, 1, 0).value();
  const auto drop = editor.add_edge(d0, d1, kNoGeometry, 0, 1, 0).value();
  const auto face = editor.add_face(surface, false, 0).value();
  const CoedgeSpec loop[] = {{keep, false}, {drop, true}};
  REQUIRE(editor.add_loop(face, loop, true).has_value());
  const Brep brep = editor.brep();

  CoedgeId keep_use{};
  CoedgeId drop_use{};
  for (const CoedgeId id : brep.coedge_ids()) {
    if (brep.coedge(id).edge == keep) keep_use = id;
    if (brep.coedge(id).edge == drop) drop_use = id;
  }
  REQUIRE(!keep_use.is_null());
  REQUIRE(!drop_use.is_null());
  CHECK(brep.coedge(keep_use).reversed != brep.coedge(drop_use).reversed);

  detail::HealEdgePair pair;
  pair.keep = keep;
  pair.drop = drop;
  pair.welds = {detail::HealVertexWeld{k0, d1},
                detail::HealVertexWeld{k1, d0}};
  CHECK(detail::heal_detail::directions_agree(brep, pair, keep_use, drop_use));
}

TEST_CASE("heal: duplicate weld mappings are inconsistent") {
  SolidEditor editor;
  const auto v0 = editor.add_vertex({0, 0, 0}, 0).value();
  const auto v1 = editor.add_vertex({1, 0, 0}, 0).value();
  const auto v2 = editor.add_vertex({2, 0, 0}, 0).value();
  const Brep brep = editor.brep();

  const auto expect_inconsistent = [&](const std::vector<detail::HealEdgePair>& pairs) {
    const auto deduped = detail::heal_dedup_welds(brep, pairs);
    REQUIRE_FALSE(deduped.has_value());
    CHECK(deduped.error().code == Errc::invalid_topology);
    CHECK(std::string_view(deduped.error().message) ==
          "heal_brep: weld map is inconsistent");
  };

  std::vector<detail::HealEdgePair> drop_to_two_keeps(2);
  drop_to_two_keeps[0].welds = {detail::HealVertexWeld{v0, v2},
                                detail::HealVertexWeld{v0, v2}};
  drop_to_two_keeps[1].welds = {detail::HealVertexWeld{v1, v2},
                                detail::HealVertexWeld{v1, v2}};
  expect_inconsistent(drop_to_two_keeps);

  std::vector<detail::HealEdgePair> keep_to_two_drops(2);
  keep_to_two_drops[0].welds = {detail::HealVertexWeld{v0, v1},
                                detail::HealVertexWeld{v0, v1}};
  keep_to_two_drops[1].welds = {detail::HealVertexWeld{v0, v2},
                                detail::HealVertexWeld{v0, v2}};
  expect_inconsistent(keep_to_two_drops);

  const VertexId out_of_range = VertexId::make(
      static_cast<std::uint32_t>(brep.vertices().size() + 7u), 0);
  std::vector<detail::HealEdgePair> bad_handle(1);
  bad_handle[0].welds = {detail::HealVertexWeld{out_of_range, v0},
                         detail::HealVertexWeld{out_of_range, v0}};
  expect_inconsistent(bad_handle);
}
