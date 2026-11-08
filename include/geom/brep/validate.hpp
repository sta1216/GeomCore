#pragma once

#include <cassert>
#include <cstddef>
#include <cstdint>
#include <vector>

#include <geom/brep/topology.hpp>
#include <geom/core/error.hpp>

namespace geom {
inline namespace v1 {

// One finding of validate_topology(). `entity` is the index() of the handle the
// finding is attributed to; 0xFFFFFFFFu means the finding is global (no single
// owning entity). `message` is a static string, never owned.
struct ValidationIssue {
  Errc code = Errc::none;
  std::uint32_t entity = 0xFFFFFFFFu;
  const char* message = "";
};

struct ValidationReport {
  std::vector<ValidationIssue> issues;
  bool ok() const noexcept { return issues.empty(); }
};

namespace detail {
// validate_topology() stops collecting after this many findings: a single
// broken reference usually cascades into many derived findings, and callers only
// need the first one (SolidEditor::commit() reports it).
inline constexpr std::size_t kMaxValidationIssues = 8;
}  // namespace detail

// Checks the M2 topology invariants of `brep`:
//
//  1. every live entity references only live handles (Coedge edge/face/loop/
//     next/prev, Loop face/start/next, Face first_loop, Edge v0/v1) and existing
//     geometry-store entries;
//  2. every loop's coedge chain is a closed next-cycle with a consistent prev
//     link, and each coedge agrees with its loop's face and loop id;
//  3. every face has a loop, owns at most one outer loop and its loop chain is
//     acyclic;
//  4. every edge is used by one or two coedges, and two uses run in opposite
//     directions (a cylinder seam is exactly this two-direction use);
//  5. every closed==true shell has only interior edges: each edge its faces use
//     is used twice inside that shell; closed==false shells may have boundary
//     edges (an open sheet). Shell contents are also validated to list no face
//     twice, and a solid may not list a shell twice.
//
// A null Loop::next is the legitimate end of a face's loop chain, and a null
// Face::first_loop denotes a face with no loops, so neither is a dangling
// reference. The checks collect at most detail::kMaxValidationIssues findings;
// the report is ordered by check and then by entity index.
//
// The function allocates (handle enumeration and per-edge use counts), so it is
// not noexcept; every other failure mode is reported as an issue instead of a
// thrown exception.
inline ValidationReport validate_topology(const Brep& brep) {
  ValidationReport report;
  const auto fail = [&report](Errc code, std::uint32_t entity, const char* message) {
    if (report.issues.size() >= detail::kMaxValidationIssues) return;
    report.issues.push_back(ValidationIssue{code, entity, message});
  };
  const auto full = [&report] {
    return report.issues.size() >= detail::kMaxValidationIssues;
  };

  // 1. References from live entities must be live. Edges first so a dangling
  // vertex is attributed to the edge that still points at it.
  for (const EdgeId id : brep.edge_ids()) {
    const Edge& edge = brep.edge(id);
    if (!brep.vertex_alive(edge.v0)) {
      fail(Errc::invalid_topology, id.index(), "edge start vertex is not alive");
    }
    if (!brep.vertex_alive(edge.v1)) {
      fail(Errc::invalid_topology, id.index(), "edge end vertex is not alive");
    }
    if (edge.curve != kNoGeometry &&
        static_cast<std::size_t>(edge.curve) >= brep.geometry().curves.size()) {
      fail(Errc::invalid_topology, id.index(), "edge curve index is out of range");
    }
  }
  for (const CoedgeId id : brep.coedge_ids()) {
    const Coedge& coedge = brep.coedge(id);
    if (!brep.edge_alive(coedge.edge)) {
      fail(Errc::invalid_topology, id.index(), "coedge edge is not alive");
    }
    if (!brep.face_alive(coedge.face)) {
      fail(Errc::invalid_topology, id.index(), "coedge face is not alive");
    }
    if (!brep.loop_alive(coedge.loop)) {
      fail(Errc::invalid_topology, id.index(), "coedge loop is not alive");
    }
    if (!brep.coedge_alive(coedge.next)) {
      fail(Errc::invalid_topology, id.index(), "coedge next is not alive");
    }
    if (!brep.coedge_alive(coedge.prev)) {
      fail(Errc::invalid_topology, id.index(), "coedge prev is not alive");
    }
    if (coedge.pcurve != kNoGeometry &&
        static_cast<std::size_t>(coedge.pcurve) >= brep.geometry().pcurves.size()) {
      fail(Errc::invalid_topology, id.index(), "coedge pcurve index is out of range");
    }
  }
  for (const LoopId id : brep.loop_ids()) {
    const Loop& loop = brep.loop(id);
    if (!brep.face_alive(loop.face)) {
      fail(Errc::invalid_topology, id.index(), "loop face is not alive");
    }
    if (!brep.coedge_alive(loop.start)) {
      fail(Errc::invalid_topology, id.index(), "loop start coedge is not alive");
    }
    if (!loop.next.is_null() && !brep.loop_alive(loop.next)) {
      fail(Errc::invalid_topology, id.index(), "loop next is not alive");
    }
  }
  for (const FaceId id : brep.face_ids()) {
    const Face& face = brep.face(id);
    if (face.surface != kNoGeometry &&
        static_cast<std::size_t>(face.surface) >= brep.geometry().surfaces.size()) {
      fail(Errc::invalid_topology, id.index(), "face surface index is out of range");
    }
    if (!face.first_loop.is_null() && !brep.loop_alive(face.first_loop)) {
      fail(Errc::invalid_topology, id.index(), "face first loop is not alive");
    }
  }
  for (const ShellId id : brep.shell_ids()) {
    const Shell& shell = brep.shell(id);
    const auto faces = brep.shell_faces(id);
    if (faces.size() != shell.face_count) {
      fail(Errc::invalid_topology, id.index(), "shell face range is out of bounds");
      continue;
    }
    bool duplicate_face = false;
    for (std::size_t i = 0; i < faces.size() && !duplicate_face; ++i) {
      for (std::size_t j = 0; j < i; ++j) {
        if (faces[i] == faces[j]) {
          duplicate_face = true;
          break;
        }
      }
    }
    if (duplicate_face) {
      fail(Errc::invalid_topology, id.index(), "shell lists a face more than once");
      continue;
    }
    for (const FaceId face : faces) {
      if (!brep.face_alive(face)) {
        fail(Errc::invalid_topology, id.index(), "shell face is not alive");
      }
    }
  }
  for (const SolidId id : brep.solid_ids()) {
    const SolidEntity& solid = brep.solid(id);
    const auto shells = brep.solid_shells(id);
    if (shells.size() != solid.shell_count) {
      fail(Errc::invalid_topology, id.index(), "solid shell range is out of bounds");
      continue;
    }
    bool duplicate_shell = false;
    for (std::size_t i = 0; i < shells.size() && !duplicate_shell; ++i) {
      for (std::size_t j = 0; j < i; ++j) {
        if (shells[i] == shells[j]) {
          duplicate_shell = true;
          break;
        }
      }
    }
    if (duplicate_shell) {
      fail(Errc::invalid_topology, id.index(), "solid lists a shell more than once");
      continue;
    }
    for (const ShellId shell : shells) {
      if (!brep.shell_alive(shell)) {
        fail(Errc::invalid_topology, id.index(), "solid shell is not alive");
      }
    }
  }
  if (full()) return report;

  // 2. Every loop's coedge chain is a closed cycle and agrees on loop and face.
  const std::size_t coedge_limit = brep.coedge_count();
  for (const LoopId loop_id : brep.loop_ids()) {
    const Loop& loop = brep.loop(loop_id);
    if (!brep.face_alive(loop.face) || !brep.coedge_alive(loop.start)) continue;
    CoedgeId current = loop.start;
    std::size_t steps = 0;
    for (;;) {
      const Coedge& coedge = brep.coedge(current);
      if (coedge.loop != loop_id || coedge.face != loop.face) {
        fail(Errc::invalid_topology, loop_id.index(),
             "loop coedge disagrees on loop or face");
        break;
      }
      const CoedgeId next = coedge.next;
      if (!brep.coedge_alive(next)) {
        fail(Errc::invalid_topology, loop_id.index(), "loop coedge next is not alive");
        break;
      }
      if (brep.coedge(next).prev != current) {
        fail(Errc::invalid_topology, loop_id.index(), "loop next/prev links disagree");
        break;
      }
      const CoedgeId prev = coedge.prev;
      if (!brep.coedge_alive(prev) || brep.coedge(prev).next != current) {
        fail(Errc::invalid_topology, loop_id.index(), "loop prev link disagrees");
        break;
      }
      if (next == loop.start) break;
      current = next;
      if (++steps > coedge_limit) {
        fail(Errc::invalid_topology, loop_id.index(),
             "loop coedge cycle does not return to its start");
        break;
      }
    }
  }
  if (full()) return report;

  // 3. Every face has at least one loop, at most one outer loop, and an acyclic
  // loop chain. A loopless face would otherwise pass check 5 and contribute an
  // Euler term without any boundary.
  const std::size_t loop_limit = brep.loop_count();
  for (const FaceId face_id : brep.face_ids()) {
    const Face& face = brep.face(face_id);
    if (face.first_loop.is_null()) {
      fail(Errc::invalid_topology, face_id.index(), "face has no loops");
      continue;
    }
    std::size_t outer_loops = 0;
    std::size_t steps = 0;
    LoopId current = face.first_loop;
    while (!current.is_null()) {
      if (!brep.loop_alive(current)) break;  // reported above
      const Loop& loop = brep.loop(current);
      if (loop.face != face_id) {
        fail(Errc::invalid_topology, current.index(), "face loop chain leaves the face");
        break;
      }
      if (loop.outer && ++outer_loops > 1) {
        fail(Errc::invalid_topology, face_id.index(),
             "face has more than one outer loop");
        break;
      }
      current = loop.next;
      if (++steps > loop_limit) {
        fail(Errc::invalid_topology, face_id.index(), "face loop chain contains a cycle");
        break;
      }
    }
  }
  if (full()) return report;

  // 4. Each edge is used by one or two coedges; two uses must run in opposite
  // directions. Counts are indexed by edge slot, so a released slot is skipped.
  std::vector<std::uint32_t> forward(brep.edges().size(), 0u);
  std::vector<std::uint32_t> reversed(brep.edges().size(), 0u);
  for (const CoedgeId id : brep.coedge_ids()) {
    const Coedge& coedge = brep.coedge(id);
    if (!brep.edge_alive(coedge.edge)) continue;  // reported above
    const std::size_t index = coedge.edge.index();
    if (index >= forward.size()) continue;
    if (coedge.reversed) {
      ++reversed[index];
    } else {
      ++forward[index];
    }
  }
  for (const EdgeId id : brep.edge_ids()) {
    const std::size_t index = id.index();
    if (index >= forward.size()) continue;
    const std::uint32_t total = forward[index] + reversed[index];
    if (total == 0u) {
      fail(Errc::invalid_topology, id.index(), "edge is not used by any coedge");
    } else if (total > 2u) {
      fail(Errc::invalid_topology, id.index(), "edge is used by more than two coedges");
    } else if (total == 2u && (forward[index] != 1u || reversed[index] != 1u)) {
      fail(Errc::invalid_topology, id.index(),
           "both coedges on an edge run in the same direction");
    }
  }
  if (full()) return report;

  // 5. A closed shell may not have a boundary edge: every edge used by one of
  // its faces must be used exactly twice inside that same shell. The count is
  // per shell, not global: two closed shells each contributing one coedge to a
  // shared edge are two open sheets, not one closed surface. Open shells
  // deliberately allow boundary edges.
  std::vector<std::uint32_t> shell_uses(brep.edges().size(), 0u);
  std::vector<std::uint32_t> touched_edges;
  for (const ShellId shell_id : brep.shell_ids()) {
    const Shell& shell = brep.shell(shell_id);
    if (!shell.closed) continue;
    const auto faces = brep.shell_faces(shell_id);
    if (faces.size() != shell.face_count) continue;  // reported above
    for (const std::uint32_t index : touched_edges) shell_uses[index] = 0u;
    touched_edges.clear();
    for (const FaceId face_id : faces) {
      if (!brep.face_alive(face_id)) continue;  // reported above
      LoopId loop_id = brep.face(face_id).first_loop;
      std::size_t loop_steps = 0;
      while (!loop_id.is_null() && brep.loop_alive(loop_id) && loop_steps <= loop_limit) {
        const Loop& loop = brep.loop(loop_id);
        if (brep.coedge_alive(loop.start)) {
          CoedgeId current = loop.start;
          std::size_t steps = 0;
          for (;;) {
            if (!brep.coedge_alive(current)) break;
            const Coedge& coedge = brep.coedge(current);
            if (brep.edge_alive(coedge.edge)) {
              const std::size_t index = coedge.edge.index();
              if (index < shell_uses.size()) {
                if (shell_uses[index]++ == 0u) {
                  touched_edges.push_back(static_cast<std::uint32_t>(index));
                }
              }
            }
            const CoedgeId next = coedge.next;
            if (!brep.coedge_alive(next) || next == loop.start) break;
            current = next;
            if (++steps > coedge_limit) break;
          }
        }
        loop_id = loop.next;
        ++loop_steps;
      }
    }
    for (const std::uint32_t index : touched_edges) {
      if (shell_uses[index] != 2u) {
        fail(Errc::not_closed, index, "closed shell has an edge without two coedges");
      }
    }
    if (full()) return report;
  }

  return report;
}

// Euler characteristic of `solid`: V - E + F - (L - F), where V/E/F are the
// deduplicated vertices/edges/faces reached through solid -> shell -> face ->
// loop -> coedge and L is the number of loops that bound those faces (so L - F
// counts inner boundaries). The value is intentionally not asserted here: a
// torus is 0, a single open face is 1 and a closed box is 2; each primitive's
// tests assert its own value.
//
// Precondition: solid_alive(solid). The function allocates its de-duplication
// flags, so it is not noexcept; a malformed topology yields whatever the walk
// could collect (validate_topology reports the actual defect).
inline std::int64_t euler_characteristic(const Brep& brep, SolidId solid) {
  assert(brep.solid_alive(solid));
  if (!brep.solid_alive(solid)) return 0;

  std::vector<bool> seen_vertices(brep.vertices().size(), false);
  std::vector<bool> seen_edges(brep.edges().size(), false);
  std::vector<bool> seen_faces(brep.faces().size(), false);
  std::int64_t vertex_total = 0;
  std::int64_t edge_total = 0;
  std::int64_t face_total = 0;
  std::int64_t loop_total = 0;

  const std::size_t coedge_limit = brep.coedge_count();
  const auto visit_loop = [&](LoopId loop_id) {
    const Loop& loop = brep.loop(loop_id);
    ++loop_total;
    if (!brep.coedge_alive(loop.start)) return;
    CoedgeId current = loop.start;
    std::size_t steps = 0;
    for (;;) {
      if (!brep.coedge_alive(current)) return;
      const Coedge& coedge = brep.coedge(current);
      if (brep.edge_alive(coedge.edge)) {
        const std::size_t edge_index = coedge.edge.index();
        if (!seen_edges[edge_index]) {
          seen_edges[edge_index] = true;
          ++edge_total;
        }
        const Edge& edge = brep.edge(coedge.edge);
        if (brep.vertex_alive(edge.v0) && !seen_vertices[edge.v0.index()]) {
          seen_vertices[edge.v0.index()] = true;
          ++vertex_total;
        }
        if (brep.vertex_alive(edge.v1) && !seen_vertices[edge.v1.index()]) {
          seen_vertices[edge.v1.index()] = true;
          ++vertex_total;
        }
      }
      const CoedgeId next = coedge.next;
      if (!brep.coedge_alive(next) || next == loop.start) return;
      current = next;
      if (++steps > coedge_limit) return;
    }
  };

  const std::size_t loop_limit = brep.loop_count();
  for (const ShellId shell_id : brep.solid_shells(solid)) {
    if (!brep.shell_alive(shell_id)) continue;
    for (const FaceId face_id : brep.shell_faces(shell_id)) {
      if (!brep.face_alive(face_id)) continue;
      // A face shared by two shells (unusual but constructible) is counted once,
      // including its loops; only the first visit walks them.
      if (seen_faces[face_id.index()]) continue;
      seen_faces[face_id.index()] = true;
      ++face_total;
      LoopId loop_id = brep.face(face_id).first_loop;
      std::size_t chain = 0;
      while (!loop_id.is_null() && brep.loop_alive(loop_id) && chain <= loop_limit) {
        visit_loop(loop_id);
        loop_id = brep.loop(loop_id).next;
        ++chain;
      }
    }
  }

  return vertex_total - edge_total + face_total - (loop_total - face_total);
}

}  // namespace v1
}  // namespace geom
