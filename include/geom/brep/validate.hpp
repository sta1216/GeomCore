#pragma once

#include <cassert>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <span>
#include <utility>
#include <vector>

#include <geom/brep/face_query.hpp>
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

// Signed volume enclosed by one shell, recovered from the direction its face
// loops run. An outward shell encloses positive volume; a cavity wall's loops
// run inward and enclose negative volume. `computable` is false when the shell
// contains curved surface geometry: the planar polygon walk cannot infer the
// orientation of a curved face.
struct validate_shell_volume {
  bool computable = false;
  double volume = 0.0;
};

// Point overview of one shell: every loop as a point ring in traversal order,
// the distinct vertex points, and their average. `ok` is false when a reference
// is dangling or a loop is not a closed ring, so callers never infer orientation
// from a partial walk.
struct validate_shell_walk {
  bool ok = false;
  std::vector<std::vector<Point3d>> rings;
  std::vector<Point3d> vertices;
  Point3d sample{};
};

inline bool validate_shell_planar(const Brep& brep, const Face& face) {
  if (face.surface == kNoGeometry) return true;
  if (static_cast<std::size_t>(face.surface) >= brep.geometry().surfaces.size()) {
    return false;
  }
  return std::holds_alternative<Plane3>(brep.geometry().surfaces[face.surface]);
}

inline validate_shell_walk validate_walk_shell(const Brep& brep,
                                               std::span<const FaceId> shell_faces) {
  validate_shell_walk walk;
  const std::size_t coedge_limit = brep.coedge_count();
  const std::size_t loop_limit = brep.loop_count();
  std::vector<bool> seen(brep.vertices().size(), false);
  for (const FaceId face_id : shell_faces) {
    if (!brep.face_alive(face_id)) return walk;
    const Face& face = brep.face(face_id);
    if (!validate_shell_planar(brep, face)) return walk;
    LoopId loop_id = face.first_loop;
    std::size_t loop_steps = 0;
    while (!loop_id.is_null()) {
      if (!brep.loop_alive(loop_id) || loop_steps++ > loop_limit) return walk;
      const Loop& loop = brep.loop(loop_id);
      if (!brep.coedge_alive(loop.start)) return walk;
      std::vector<Point3d> ring;
      CoedgeId current = loop.start;
      std::size_t steps = 0;
      for (;;) {
        if (!brep.coedge_alive(current)) return walk;
        const Coedge& coedge = brep.coedge(current);
        if (!brep.edge_alive(coedge.edge)) return walk;
        const Edge& edge = brep.edge(coedge.edge);
        const VertexId vertex = coedge.reversed ? edge.v1 : edge.v0;
        if (!brep.vertex_alive(vertex)) return walk;
        if (vertex.index() >= seen.size()) return walk;
        ring.push_back(brep.vertex(vertex).point);
        if (!seen[vertex.index()]) {
          seen[vertex.index()] = true;
          walk.vertices.push_back(brep.vertex(vertex).point);
        }
        const CoedgeId next = coedge.next;
        if (!brep.coedge_alive(next)) return walk;
        if (next == loop.start) break;
        current = next;
        if (++steps > coedge_limit) return walk;
      }
      if (ring.size() < 3u) return walk;
      walk.rings.push_back(std::move(ring));
      loop_id = loop.next;
    }
  }
  if (walk.rings.empty() || walk.vertices.empty()) return walk;
  for (const Point3d& point : walk.vertices) {
    walk.sample[0] += point[0];
    walk.sample[1] += point[1];
    walk.sample[2] += point[2];
  }
  const double count = static_cast<double>(walk.vertices.size());
  walk.sample[0] /= count;
  walk.sample[1] /= count;
  walk.sample[2] /= count;
  walk.ok = std::isfinite(walk.sample[0]) && std::isfinite(walk.sample[1]) &&
            std::isfinite(walk.sample[2]);
  return walk;
}

inline validate_shell_volume validate_signed_shell_volume(
    const Brep& brep, std::span<const FaceId> shell_faces) {
  validate_shell_volume result;
  const validate_shell_walk walk = validate_walk_shell(brep, shell_faces);
  if (!walk.ok) return result;
  const Point3d& reference = walk.rings.front().front();
  double volume = 0.0;
  for (const std::vector<Point3d>& ring : walk.rings) {
    const Vec3d origin = ring.front() - reference;
    for (std::size_t i = 1; i + 1u < ring.size(); ++i) {
      const Vec3d first = ring[i] - reference;
      const Vec3d second = ring[i + 1u] - reference;
      volume += dot(origin, cross(first, second)) / 6.0;
    }
  }
  if (!std::isfinite(volume)) return result;
  result.computable = true;
  result.volume = volume;
  return result;
}

// True when every probe point of `walk` (its sample and each vertex) lies
// strictly inside `shell`; a point on the shell surface counts as outside.
inline Result<bool> validate_walk_inside(const Brep& brep, ShellId shell,
                                         const validate_shell_walk& walk, double tol) {
  const Result<bool> sample = point_in_shell(brep, shell, walk.sample, tol);
  if (!sample) return sample;
  if (!sample.value()) return false;
  for (const Point3d& vertex : walk.vertices) {
    const Result<bool> inside = point_in_shell(brep, shell, vertex, tol);
    if (!inside) return inside;
    if (!inside.value()) return false;
  }
  return true;
}

// True when any probe point of `walk` lies strictly inside `shell`. Two cavities
// overlap when either one's shell reaches into the other.
inline Result<bool> validate_walk_touches(const Brep& brep, ShellId shell,
                                          const validate_shell_walk& walk, double tol) {
  const Result<bool> sample = point_in_shell(brep, shell, walk.sample, tol);
  if (!sample) return sample;
  if (sample.value()) return true;
  for (const Point3d& vertex : walk.vertices) {
    const Result<bool> inside = point_in_shell(brep, shell, vertex, tol);
    if (!inside) return inside;
    if (inside.value()) return true;
  }
  return false;
}
}  // namespace detail

// A positive parent tolerance more than this factor above a positive child
// tolerance is suspicious growth rather than a legitimate refinement (§4.6).
inline constexpr double kToleranceGrowthFactor = 10.0;

// Checks the M2 topology invariants of `brep`:
//
//  1. every live entity references only live handles (Coedge edge/face/loop/
//     next/prev, Loop face/start/next, Face first_loop, Edge v0/v1) and existing
//     geometry-store entries;
//  2. every loop's coedge chain is a closed next-cycle with a consistent prev
//     link, and each coedge agrees with its loop's face and loop id;
//  3. every face has a loop, owns exactly one outer loop and its loop chain is
//     acyclic;
//  4. every edge is used by one or two coedges, and two uses run in opposite
//     directions (a cylinder seam is exactly this two-direction use);
//  5. every closed==true shell has only interior edges: each edge its faces use
//     is used twice inside that shell; closed==false shells may have boundary
//     edges (an open sheet). Shell contents are also validated to list no face
//     twice, and a solid may not list a shell twice.
//  6. every solid has at least one shell. A solid with several shells is a
//     solid with cavities: all shells must be closed, exactly one must wind
//     outward (a positive signed volume), and every other shell must wind
//     inward (a negative signed volume). Each cavity must lie strictly inside
//     the outer shell and outside the other cavities. A single-shell solid
//     keeps the M2 behavior: open sheets and curved shells are accepted without
//     an orientation test. Shell orientation is recovered from the loop
//     winding, so a multi-shell solid whose faces are curved cannot be
//     classified and is reported as Errc::unsupported.
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
    if (shells.empty()) {
      fail(Errc::invalid_topology, id.index(), "solid has no shells");
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
    if (current.is_null() && outer_loops == 0) {
      fail(Errc::invalid_topology, face_id.index(), "face has no outer loop");
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

  // 6. A multi-shell solid is a solid with cavities: exactly one shell winds
  // outward and every other shell is a cavity strictly inside it. A single-shell
  // solid keeps the M2 behavior, so open sheets and curved shells are accepted
  // without an orientation test.
  for (const SolidId id : brep.solid_ids()) {
    const SolidEntity& solid = brep.solid(id);
    const auto shells = brep.solid_shells(id);
    if (shells.size() != solid.shell_count || shells.size() < 2u) continue;

    bool closed = true;
    for (const ShellId shell : shells) {
      if (!brep.shell(shell).closed) {
        fail(Errc::invalid_topology, id.index(),
             "solid with multiple shells has an open shell");
        closed = false;
        break;
      }
    }
    if (!closed) continue;

    std::vector<detail::validate_shell_volume> volumes;
    volumes.reserve(shells.size());
    bool orientation_known = true;
    for (const ShellId shell : shells) {
      const detail::validate_shell_volume volume =
          detail::validate_signed_shell_volume(brep, brep.shell_faces(shell));
      if (!volume.computable) {
        orientation_known = false;
        break;
      }
      volumes.push_back(volume);
    }
    if (!orientation_known) {
      fail(Errc::unsupported, id.index(),
           "multi-shell shell orientation is not computable");
      continue;
    }

    std::size_t outer_count = 0;
    std::size_t outer = 0;
    for (std::size_t i = 0; i < volumes.size(); ++i) {
      if (volumes[i].volume > 0.0) {
        ++outer_count;
        outer = i;
      }
    }
    if (outer_count == 0u) {
      fail(Errc::invalid_topology, id.index(), "solid has no outer shell");
      continue;
    }
    if (outer_count > 1u) {
      fail(Errc::invalid_topology, id.index(), "solid has more than one outer shell");
      continue;
    }

    double tol = solid.tol;
    if (!std::isfinite(tol) || tol < 0.0) tol = 0.0;
    std::vector<ShellId> cavity_shells;
    std::vector<detail::validate_shell_walk> cavity_walks;
    bool valid = true;
    for (std::size_t i = 0; i < shells.size() && valid; ++i) {
      if (i == outer) continue;
      detail::validate_shell_walk walk =
          detail::validate_walk_shell(brep, brep.shell_faces(shells[i]));
      if (!walk.ok) {
        fail(Errc::invalid_topology, id.index(), "cavity shell cannot be walked");
        valid = false;
        break;
      }
      const Result<bool> inside =
          detail::validate_walk_inside(brep, shells[outer], walk, tol);
      if (!inside) {
        fail(inside.error().code, id.index(), inside.error().message);
        valid = false;
        break;
      }
      if (!inside.value()) {
        fail(Errc::invalid_topology, id.index(),
             "cavity shell is not contained in the outer shell");
        valid = false;
        break;
      }
      cavity_shells.push_back(shells[i]);
      cavity_walks.push_back(std::move(walk));
    }
    if (!valid) continue;

    for (std::size_t i = 0; i < cavity_shells.size() && valid; ++i) {
      for (std::size_t j = i + 1u; j < cavity_shells.size() && valid; ++j) {
        const Result<bool> first =
            detail::validate_walk_touches(brep, cavity_shells[j], cavity_walks[i], tol);
        if (!first) {
          fail(first.error().code, id.index(), first.error().message);
          valid = false;
          break;
        }
        if (first.value()) {
          fail(Errc::invalid_topology, id.index(), "cavity shells overlap");
          valid = false;
          break;
        }
        const Result<bool> second =
            detail::validate_walk_touches(brep, cavity_shells[i], cavity_walks[j], tol);
        if (!second) {
          fail(second.error().code, id.index(), second.error().message);
          valid = false;
          break;
        }
        if (second.value()) {
          fail(Errc::invalid_topology, id.index(), "cavity shells overlap");
          valid = false;
          break;
        }
      }
    }
  }
  if (full()) return report;

  return report;
}

// Checks the M2 tolerance chain (§4.6):
//
//  1. every tolerance is finite and non-negative;
//  2. an edge's tolerance is at least the tolerance of either of its vertices;
//  3. for every coedge, its face's tolerance is at least the edge's tolerance,
//     and the coedge's own tolerance is at most the face's tolerance (a pcurve
//     tolerance is a face-local quantity, so it may not exceed its face);
//  4. a solid's tolerance is at least every face tolerance in its shells;
//  5. growth: a positive parent tolerance greater than kToleranceGrowthFactor
//     times a positive child tolerance is suspicious, with message
//     "tolerance growth exceeds 10x".
//
// Each finding is attributed to the finer (child) entity of the failing
// relation: the vertex for rule 2, the edge for the face/edge half of rule 3,
// the coedge for the coedge/face half, the face for rule 4, and the child of the
// growing pair for rule 5. Dangling references are left to validate_topology()
// and skipped here instead of asserted. Like validate_topology(), this function
// allocates and collects at most detail::kMaxValidationIssues findings, ordered
// by rule and then by entity index.
inline ValidationReport validate_tolerances(const Brep& brep) {
  ValidationReport report;
  const auto fail = [&report](Errc code, std::uint32_t entity, const char* message) {
    if (report.issues.size() >= detail::kMaxValidationIssues) return;
    report.issues.push_back(ValidationIssue{code, entity, message});
  };
  const auto full = [&report] {
    return report.issues.size() >= detail::kMaxValidationIssues;
  };

  // 1. Tolerances must be finite and non-negative.
  for (const VertexId id : brep.vertex_ids()) {
    const double tol = brep.vertex(id).tol;
    if (!std::isfinite(tol) || tol < 0.0) {
      fail(Errc::tolerance_violation, id.index(), "tolerance must be finite and non-negative");
    }
  }
  for (const EdgeId id : brep.edge_ids()) {
    const double tol = brep.edge(id).tol;
    if (!std::isfinite(tol) || tol < 0.0) {
      fail(Errc::tolerance_violation, id.index(), "tolerance must be finite and non-negative");
    }
  }
  for (const CoedgeId id : brep.coedge_ids()) {
    const double tol = brep.coedge(id).tol;
    if (!std::isfinite(tol) || tol < 0.0) {
      fail(Errc::tolerance_violation, id.index(), "tolerance must be finite and non-negative");
    }
  }
  for (const FaceId id : brep.face_ids()) {
    const double tol = brep.face(id).tol;
    if (!std::isfinite(tol) || tol < 0.0) {
      fail(Errc::tolerance_violation, id.index(), "tolerance must be finite and non-negative");
    }
  }
  for (const SolidId id : brep.solid_ids()) {
    const double tol = brep.solid(id).tol;
    if (!std::isfinite(tol) || tol < 0.0) {
      fail(Errc::tolerance_violation, id.index(), "tolerance must be finite and non-negative");
    }
  }
  if (full()) return report;

  // 2. Edge tolerances dominate their bounding vertices. One finding per edge,
  // attributed to the vertex with the largest offending tolerance.
  for (const EdgeId id : brep.edge_ids()) {
    const Edge& edge = brep.edge(id);
    VertexId offender{};
    bool has_offender = false;
    const auto consider = [&](VertexId vertex_id) {
      if (!brep.vertex_alive(vertex_id)) return;
      const double vertex_tol = brep.vertex(vertex_id).tol;
      if (vertex_tol > edge.tol && (!has_offender || vertex_tol > brep.vertex(offender).tol)) {
        offender = vertex_id;
        has_offender = true;
      }
    };
    consider(edge.v0);
    consider(edge.v1);
    if (has_offender) {
      fail(Errc::tolerance_violation, offender.index(),
           "vertex tolerance exceeds the edge tolerance");
    }
  }
  if (full()) return report;

  // 3. Face tolerances dominate the tolerances of the edges they use, and a
  // coedge tolerance stays inside its face. Each failing coedge use is reported:
  // the same edge may be used by two faces with different tolerances.
  for (const CoedgeId id : brep.coedge_ids()) {
    const Coedge& coedge = brep.coedge(id);
    if (!brep.edge_alive(coedge.edge) || !brep.face_alive(coedge.face)) continue;
    const double edge_tol = brep.edge(coedge.edge).tol;
    const double face_tol = brep.face(coedge.face).tol;
    if (face_tol < edge_tol) {
      fail(Errc::tolerance_violation, coedge.edge.index(),
           "edge tolerance exceeds the face tolerance");
    }
    if (coedge.tol > face_tol) {
      fail(Errc::tolerance_violation, id.index(), "coedge tolerance exceeds the face tolerance");
    }
  }
  if (full()) return report;

  // 4. Solid tolerances dominate every face tolerance in their shells. One
  // finding per solid, attributed to the face with the largest offending
  // tolerance.
  for (const SolidId id : brep.solid_ids()) {
    const SolidEntity& solid = brep.solid(id);
    const auto shells = brep.solid_shells(id);
    if (shells.size() != solid.shell_count) continue;  // reported by rule-check 1
    FaceId offender{};
    double worst = -1.0;
    bool has_offender = false;
    for (const ShellId shell_id : shells) {
      if (!brep.shell_alive(shell_id)) continue;
      for (const FaceId face_id : brep.shell_faces(shell_id)) {
        if (!brep.face_alive(face_id)) continue;
        const double face_tol = brep.face(face_id).tol;
        if (!has_offender || face_tol > worst) {
          worst = face_tol;
          offender = face_id;
          has_offender = true;
        }
      }
    }
    if (has_offender && worst > solid.tol) {
      fail(Errc::tolerance_violation, offender.index(),
           "face tolerance exceeds the solid tolerance");
    }
  }
  if (full()) return report;

  // 5. Growth. A parent that needs a tolerance an order of magnitude larger
  // than its child suggests the child was modeled at the wrong scale.
  for (const EdgeId id : brep.edge_ids()) {
    const Edge& edge = brep.edge(id);
    if (!(edge.tol > 0.0)) continue;
    const VertexId ends[2] = {edge.v0, edge.v1};
    for (const VertexId vertex_id : ends) {
      if (!brep.vertex_alive(vertex_id)) continue;
      const double vertex_tol = brep.vertex(vertex_id).tol;
      if (vertex_tol > 0.0 && edge.tol > kToleranceGrowthFactor * vertex_tol) {
        fail(Errc::tolerance_violation, vertex_id.index(), "tolerance growth exceeds 10x");
      }
    }
  }
  for (const CoedgeId id : brep.coedge_ids()) {
    const Coedge& coedge = brep.coedge(id);
    if (!brep.edge_alive(coedge.edge) || !brep.face_alive(coedge.face)) continue;
    const double edge_tol = brep.edge(coedge.edge).tol;
    const double face_tol = brep.face(coedge.face).tol;
    if (edge_tol > 0.0 && face_tol > 0.0 && face_tol > kToleranceGrowthFactor * edge_tol) {
      fail(Errc::tolerance_violation, coedge.edge.index(), "tolerance growth exceeds 10x");
    }
  }
  for (const SolidId id : brep.solid_ids()) {
    const SolidEntity& solid = brep.solid(id);
    if (!(solid.tol > 0.0)) continue;
    for (const ShellId shell_id : brep.solid_shells(id)) {
      if (!brep.shell_alive(shell_id)) continue;
      for (const FaceId face_id : brep.shell_faces(shell_id)) {
        if (!brep.face_alive(face_id)) continue;
        const double face_tol = brep.face(face_id).tol;
        if (face_tol > 0.0 && solid.tol > kToleranceGrowthFactor * face_tol) {
          fail(Errc::tolerance_violation, face_id.index(), "tolerance growth exceeds 10x");
        }
      }
    }
  }

  return report;
}

// Euler characteristic of `solid`: V - E + F - (L - F), where V/E/F are the
// deduplicated vertices/edges/faces reached through solid -> shell -> face ->
// loop -> coedge and L is the number of loops that bound those faces (so L - F
// counts inner boundaries). For a closed solid the value reduces to
// 2 * (S - H), where S is the number of shells and H the number of handles: a
// box is 2, a through-hole or torus is 0, a solid with one fully enclosed cavity
// has two shells and is 4, and a single open face is 1. The value is
// intentionally not asserted here: each primitive's tests assert its own value.
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
