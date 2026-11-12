#pragma once

// Healing surface for B-rep solids.
//
// `heal_brep` closes small gaps along the boundary of a solid: coincident
// vertices are welded, boundary edges of neighbouring faces are sewn, and the
// remaining open edges are reported. Repairing a pcurve projects the edge's two
// endpoint vertices onto the surface at `HealOptions::tolerance` resolution;
// gaps wider than `HealOptions::max_gap` are left alone and reported as issues
// so a caller can decide whether the tolerance is appropriate.
//
// `HealReport` summarises what a run changed and how close it came to the
// requested tolerances. `detail::heal_boundary_edges` exposes the boundary
// discovery used by the first phase so it can be tested in isolation.

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <span>
#include <utility>
#include <vector>

#include <geom/brep/face_uv.hpp>
#include <geom/brep/pcurve2d.hpp>
#include <geom/brep/solid_editor.hpp>
#include <geom/brep/topology.hpp>
#include <geom/brep/validate.hpp>
#include <geom/core/result.hpp>
#include <geom/primitives/vec.hpp>

namespace geom {
inline namespace v1 {

// Why a healing attempt could not repair one spot of the model.
enum class HealIssueKind : std::uint8_t {
  none,
  gap_too_large,
  ambiguous_gap,
  inconsistent_orientation,
  invalid_result,
};

inline const char* to_string(HealIssueKind kind) noexcept {
  switch (kind) {
    case HealIssueKind::none:
      return "none";
    case HealIssueKind::gap_too_large:
      return "gap_too_large";
    case HealIssueKind::ambiguous_gap:
      return "ambiguous_gap";
    case HealIssueKind::inconsistent_orientation:
      return "inconsistent_orientation";
    case HealIssueKind::invalid_result:
      return "invalid_result";
    default:
      break;
  }
  return "unknown";
}

// Sentinel stored in HealIssue when the issue does not name an entity.
inline constexpr std::uint32_t kNoEntity = 0xFFFFFFFFu;

struct HealOptions {
  // Tolerance used when projecting an edge's endpoints onto their surface to
  // repair a pcurve.
  double tolerance = 1e-9;
  // Gaps wider than this are reported and left open.
  double max_gap = 1e-6;
  // Largest tolerance a welded vertex or sewn edge may end up with.
  double max_tolerance = 1e-6;
};

struct HealIssue {
  HealIssueKind kind = HealIssueKind::none;
  Point3d location{};
  double gap = 0.0;
  std::uint32_t entity_a = kNoEntity;
  std::uint32_t entity_b = kNoEntity;
  const char* message = "";
};

struct HealReport {
  std::uint32_t boundary_edges = 0;
  std::uint32_t welded_vertices = 0;
  std::uint32_t sewn_edges = 0;
  double max_closed_gap = 0.0;
  double max_tolerance = 0.0;
  bool changed = false;
  std::vector<HealIssue> issues;
};

struct HealResult {
  Brep brep;
  HealReport report;
};

// Returns a healed copy of `brep` and a report of what changed. A solid that
// already closes is returned unchanged with `report.changed == false`; issues
// that cannot be repaired within `options` are reported instead of failing the
// run. `heal_brep` fails only when the input is unusable (for example an
// unknown solid handle).
Result<HealResult> heal_brep(const Brep& brep, SolidId solid, const HealOptions& options);

namespace detail {

struct HealVertexWeld {
  VertexId keep{};  // lower index
  VertexId drop{};  // higher index
};

struct HealEdgePair {
  EdgeId keep{};      // lower index
  EdgeId drop{};      // higher index
  std::array<HealVertexWeld, 2> welds{};  // endpoint correspondence, keep < drop
  bool same_direction = false;  // coedge traversal directions agree
  double gap = 0.0;             // max endpoint distance under the chosen mapping
};

struct HealPlan {
  std::vector<HealEdgePair> pairs;
  std::vector<HealIssue> issues;
};

// Builds the deduplicated vertex-weld map for a plan. Entries whose keep or
// drop handle is no longer alive are skipped. A drop mapped to two different
// keeps (or a keep to two different drops) is inconsistent.
inline Result<std::vector<HealVertexWeld>> heal_dedup_welds(
    const Brep& brep, std::span<const HealEdgePair> pairs) {
  const std::size_t vertex_slots = brep.vertices().size();
  std::vector<std::uint32_t> drop_to_keep(vertex_slots, kNoEntity);
  std::vector<std::uint32_t> keep_to_drop(vertex_slots, kNoEntity);
  std::vector<HealVertexWeld> welds;
  for (const HealEdgePair& pair : pairs) {
    for (const HealVertexWeld& weld : pair.welds) {
      const std::uint32_t keep_index = weld.keep.index();
      const std::uint32_t drop_index = weld.drop.index();
      if (keep_index >= vertex_slots || drop_index >= vertex_slots) {
        return Error{Errc::invalid_topology,
                     "heal_brep: weld map is inconsistent", 0, 0};
      }
      if (!brep.vertex_alive(weld.keep) || !brep.vertex_alive(weld.drop)) {
        continue;
      }
      if (drop_to_keep[drop_index] != kNoEntity &&
          drop_to_keep[drop_index] != keep_index) {
        return Error{Errc::invalid_topology,
                     "heal_brep: weld map is inconsistent", 0, 0};
      }
      if (keep_to_drop[keep_index] != kNoEntity &&
          keep_to_drop[keep_index] != drop_index) {
        return Error{Errc::invalid_topology,
                     "heal_brep: weld map is inconsistent", 0, 0};
      }
      if (drop_to_keep[drop_index] == kNoEntity) {
        drop_to_keep[drop_index] = keep_index;
        keep_to_drop[keep_index] = drop_index;
        welds.push_back(weld);
      }
    }
  }
  return welds;
}

// Applies welds in descending drop order so a chain v0<-v1<-v2 collapses to a
// single live keep. Handles that became dead earlier in the chain are skipped.
inline Result<std::vector<HealVertexWeld>> heal_apply_welds(
    SolidEditor& editor, std::span<const HealVertexWeld> welds) {
  std::vector<HealVertexWeld> ordered(welds.begin(), welds.end());
  std::sort(ordered.begin(), ordered.end(),
            [](const HealVertexWeld& a, const HealVertexWeld& b) {
              return a.drop.index() > b.drop.index();
            });
  std::vector<HealVertexWeld> applied;
  applied.reserve(ordered.size());
  for (const HealVertexWeld& weld : ordered) {
    if (!editor.brep().vertex_alive(weld.keep) ||
        !editor.brep().vertex_alive(weld.drop)) {
      continue;
    }
    const Result<void> merged = editor.merge_vertices(weld.keep, weld.drop);
    if (!merged) return merged.error();
    applied.push_back(weld);
  }
  return applied;
}

inline std::uint32_t heal_live_weld_count(
    const Brep& brep, std::span<const HealVertexWeld> welds) {
  std::vector<char> seen(brep.vertices().size(), 0);
  std::uint32_t count = 0;
  for (const HealVertexWeld& weld : welds) {
    if (!brep.vertex_alive(weld.keep) || weld.keep.index() >= seen.size()) {
      continue;
    }
    if (seen[weld.keep.index()] == 0) {
      seen[weld.keep.index()] = 1;
      ++count;
    }
  }
  return count;
}

namespace heal_detail {

inline double point_distance(const Point3d& a, const Point3d& b) {
  return norm2(b - a);
}

inline Point3d midpoint(const Point3d& a, const Point3d& b) {
  return Point3d{(a[0] + b[0]) * 0.5, (a[1] + b[1]) * 0.5,
                 (a[2] + b[2]) * 0.5};
}

inline HealVertexWeld normalize_weld(VertexId a, VertexId b) {
  if (a.index() < b.index()) return HealVertexWeld{a, b};
  return HealVertexWeld{b, a};
}

// Maps a vertex of the kept edge to its counterpart on the dropped edge.
inline VertexId pair_counterpart(const HealEdgePair& pair, VertexId vertex,
                                 bool keep_to_drop) {
  for (const HealVertexWeld& weld : pair.welds) {
    if (keep_to_drop && weld.keep == vertex) return weld.drop;
    if (!keep_to_drop && weld.drop == vertex) return weld.keep;
  }
  return vertex;
}

// Resolves the vertex that survives after this pair's welds are applied.
inline VertexId pair_survivor(const HealEdgePair& pair, VertexId vertex) {
  for (int step = 0; step < 3; ++step) {
    VertexId next = vertex;
    for (const HealVertexWeld& weld : pair.welds) {
      if (vertex == weld.drop) next = weld.keep;
    }
    if (next == vertex) break;
    vertex = next;
  }
  return vertex;
}

// The direction predicate used by SolidEditor::merge_edges, lifted through the
// endpoint welds. Two uses are compatible when their post-weld traversal starts
// at the same kept vertex (and will therefore end up with opposite flags).
inline bool directions_agree(const Brep& brep, const HealEdgePair& pair,
                             CoedgeId keep_coedge, CoedgeId drop_coedge) {
  const Edge& keep_edge = brep.edge(pair.keep);
  const Edge& drop_edge = brep.edge(pair.drop);
  const Coedge& keep_use = brep.coedge(keep_coedge);
  const Coedge& drop_use = brep.coedge(drop_coedge);
  const VertexId keep_start = keep_use.reversed ? keep_edge.v1 : keep_edge.v0;
  const VertexId drop_start = drop_use.reversed ? drop_edge.v1 : drop_edge.v0;
  const VertexId mapped_keep =
      pair_survivor(pair, pair_counterpart(pair, keep_start, true));
  const VertexId mapped_drop =
      pair_survivor(pair, pair_counterpart(pair, drop_start, false));
  return mapped_keep == mapped_drop;
}

inline Point3d edge_midpoint(const Brep& brep, EdgeId edge_id) {
  const Edge& edge = brep.edge(edge_id);
  return midpoint(brep.vertex(edge.v0).point, brep.vertex(edge.v1).point);
}

struct Candidate {
  HealEdgePair pair;
  std::size_t keep_position = 0;
  std::size_t drop_position = 0;
  double sum = 0.0;
};

inline bool candidate_less(const Candidate& a, const Candidate& b) {
  if (a.pair.gap != b.pair.gap) return a.pair.gap < b.pair.gap;
  if (a.sum != b.sum) return a.sum < b.sum;
  if (a.pair.keep.index() != b.pair.keep.index()) {
    return a.pair.keep.index() < b.pair.keep.index();
  }
  return a.pair.drop.index() < b.pair.drop.index();
}

inline bool gap_issue_less(const Candidate& a, const Candidate& b) {
  if (a.pair.gap != b.pair.gap) return a.pair.gap < b.pair.gap;
  if (a.pair.keep.index() != b.pair.keep.index()) {
    return a.pair.keep.index() < b.pair.keep.index();
  }
  return a.pair.drop.index() < b.pair.drop.index();
}

// Every coedge used by a face of `shell`, in face/loop order.
inline void collect_shell_coedges(const Brep& brep, ShellId shell,
                                  std::vector<CoedgeId>& out) {
  for (const FaceId face_id : brep.shell_faces(shell)) {
    if (!brep.face_alive(face_id)) continue;
    LoopId loop_id = brep.face(face_id).first_loop;
    std::size_t loops_seen = 0;
    while (brep.loop_alive(loop_id) && loops_seen <= brep.loop_count()) {
      const Loop& loop = brep.loop(loop_id);
      if (brep.coedge_alive(loop.start)) {
        CoedgeId current = loop.start;
        std::size_t coedges_seen = 0;
        while (brep.coedge_alive(current) && coedges_seen <= brep.coedge_count()) {
          out.push_back(current);
          if (brep.coedge(current).next == loop.start) break;
          current = brep.coedge(current).next;
          ++coedges_seen;
        }
      }
      loop_id = loop.next;
      ++loops_seen;
    }
  }
}

inline bool shell_is_closed_by_uses(const Brep& brep, ShellId shell) {
  std::vector<CoedgeId> coedges;
  collect_shell_coedges(brep, shell, coedges);
  std::vector<std::uint32_t> uses(brep.edges().size(), 0u);
  for (const CoedgeId coedge_id : coedges) {
    const Coedge& coedge = brep.coedge(coedge_id);
    if (brep.edge_alive(coedge.edge) && coedge.edge.index() < uses.size()) {
      ++uses[coedge.edge.index()];
    }
  }
  for (const CoedgeId coedge_id : coedges) {
    const Coedge& coedge = brep.coedge(coedge_id);
    if (!brep.edge_alive(coedge.edge) || coedge.edge.index() >= uses.size() ||
        uses[coedge.edge.index()] != 2u) {
      return false;
    }
  }
  return true;
}

inline std::size_t selected_weld(const Brep& brep, const HealEdgePair& pair) {
  const double d0 = point_distance(brep.vertex(pair.welds[0].keep).point,
                                   brep.vertex(pair.welds[0].drop).point);
  const double d1 = point_distance(brep.vertex(pair.welds[1].keep).point,
                                   brep.vertex(pair.welds[1].drop).point);
  return d0 <= d1 ? 0u : 1u;
}

}  // namespace heal_detail

inline Result<HealPlan> heal_plan(const Brep& brep,
                                  std::span<const EdgeId> boundary,
                                  double max_gap) {
  HealPlan plan;
  if (boundary.empty()) return plan;

  // When heal_brep calls this, it has already rejected boundary edges shared
  // with another solid, so every coedge of a boundary edge belongs to the
  // target solid. The first use is therefore a sound direction reference.
  std::vector<CoedgeId> first_coedge(brep.edges().size(), CoedgeId{});
  for (const CoedgeId coedge_id : brep.coedge_ids()) {
    const Coedge& coedge = brep.coedge(coedge_id);
    if (coedge.edge.index() < first_coedge.size() &&
        first_coedge[coedge.edge.index()].is_null()) {
      first_coedge[coedge.edge.index()] = coedge_id;
    }
  }

  std::vector<std::uint32_t> position(brep.edges().size(), kNoEntity);
  for (std::size_t index = 0; index < boundary.size(); ++index) {
    const EdgeId edge_id = boundary[index];
    if (!brep.edge_alive(edge_id)) {
      return Error{Errc::invalid_topology,
                   "heal_plan: boundary edge is not alive", 0, 0};
    }
    if (position[edge_id.index()] != kNoEntity) continue;
    position[edge_id.index()] = static_cast<std::uint32_t>(index);
    const Edge& edge = brep.edge(edge_id);
    if (!brep.vertex_alive(edge.v0) || !brep.vertex_alive(edge.v1)) {
      return Error{Errc::invalid_topology,
                   "heal_plan: boundary edge has a dead vertex", 0, 0};
    }
    if (heal_detail::point_distance(brep.vertex(edge.v0).point,
                                    brep.vertex(edge.v1).point) == 0.0) {
      return Error{Errc::invalid_topology,
                   "heal_plan: boundary edge has zero length", 0, 0};
    }
    if (first_coedge[edge_id.index()].is_null()) {
      return Error{Errc::invalid_topology,
                   "heal_plan: boundary edge has no coedge", 0, 0};
    }
  }

  std::vector<heal_detail::Candidate> candidates;
  std::vector<std::vector<std::size_t>> touching(boundary.size());
  for (std::size_t i = 0; i < boundary.size(); ++i) {
    const EdgeId first_id = boundary[i];
    if (position[first_id.index()] != i) continue;
    const Edge& first = brep.edge(first_id);
    const FaceId first_face = brep.coedge(first_coedge[first_id.index()]).face;
    for (std::size_t j = i + 1; j < boundary.size(); ++j) {
      const EdgeId second_id = boundary[j];
      if (position[second_id.index()] != j) continue;
      if (first_id == second_id) continue;
      const Edge& second = brep.edge(second_id);

      if (first.v0 == second.v0 || first.v0 == second.v1 ||
          first.v1 == second.v0 || first.v1 == second.v1) {
        continue;
      }
      const FaceId second_face =
          brep.coedge(first_coedge[second_id.index()]).face;
      if (first_face == second_face) continue;

      EdgeId keep = first_id;
      EdgeId drop = second_id;
      std::size_t keep_position = i;
      std::size_t drop_position = j;
      if (keep.index() > drop.index()) {
        std::swap(keep, drop);
        std::swap(keep_position, drop_position);
      }
      const Edge& keep_edge = brep.edge(keep);
      const Edge& drop_edge = brep.edge(drop);
      const Point3d& k0 = brep.vertex(keep_edge.v0).point;
      const Point3d& k1 = brep.vertex(keep_edge.v1).point;
      const Point3d& d0 = brep.vertex(drop_edge.v0).point;
      const Point3d& d1 = brep.vertex(drop_edge.v1).point;

      const double d00 = heal_detail::point_distance(k0, d0);
      const double d11 = heal_detail::point_distance(k1, d1);
      const double d01 = heal_detail::point_distance(k0, d1);
      const double d10 = heal_detail::point_distance(k1, d0);
      const double same_gap = std::max(d00, d11);
      const double flipped_gap = std::max(d01, d10);
      const double same_sum = d00 + d11;
      const double flipped_sum = d01 + d10;
      const bool same_mapping =
          same_gap < flipped_gap ||
          (same_gap == flipped_gap && same_sum <= flipped_sum);

      heal_detail::Candidate candidate;
      candidate.pair.keep = keep;
      candidate.pair.drop = drop;
      candidate.keep_position = keep_position;
      candidate.drop_position = drop_position;
      if (same_mapping) {
        candidate.pair.welds = {
            heal_detail::normalize_weld(keep_edge.v0, drop_edge.v0),
            heal_detail::normalize_weld(keep_edge.v1, drop_edge.v1)};
        candidate.pair.gap = same_gap;
        candidate.sum = same_sum;
      } else {
        candidate.pair.welds = {
            heal_detail::normalize_weld(keep_edge.v0, drop_edge.v1),
            heal_detail::normalize_weld(keep_edge.v1, drop_edge.v0)};
        candidate.pair.gap = flipped_gap;
        candidate.sum = flipped_sum;
      }
      candidate.pair.same_direction = heal_detail::directions_agree(
          brep, candidate.pair, first_coedge[keep.index()],
          first_coedge[drop.index()]);

      const std::size_t candidate_index = candidates.size();
      candidates.push_back(candidate);
      touching[keep_position].push_back(candidate_index);
      touching[drop_position].push_back(candidate_index);
    }
  }

  std::vector<std::vector<std::size_t>> accepted(boundary.size());
  for (std::size_t candidate_index = 0; candidate_index < candidates.size();
       ++candidate_index) {
    const heal_detail::Candidate& candidate = candidates[candidate_index];
    if (candidate.pair.gap <= max_gap) {
      accepted[candidate.keep_position].push_back(candidate_index);
      accepted[candidate.drop_position].push_back(candidate_index);
    }
  }

  constexpr double kAmbiguityGap = 1e-12;
  std::vector<char> ambiguous(boundary.size(), 0);
  for (std::size_t index = 0; index < boundary.size(); ++index) {
    if (accepted[index].empty()) continue;
    std::sort(accepted[index].begin(), accepted[index].end(),
              [&](std::size_t a, std::size_t b) {
                return heal_detail::candidate_less(candidates[a], candidates[b]);
              });
    if (accepted[index].size() < 2u) continue;
    const double best_gap = candidates[accepted[index][0]].pair.gap;
    const double next_gap = candidates[accepted[index][1]].pair.gap;
    if (std::fabs(best_gap - next_gap) > kAmbiguityGap) continue;

    ambiguous[index] = 1;
    const EdgeId edge_id = boundary[index];
    HealIssue issue;
    issue.kind = HealIssueKind::ambiguous_gap;
    issue.location = heal_detail::edge_midpoint(brep, edge_id);
    issue.gap = best_gap;
    issue.entity_a = edge_id.index();
    plan.issues.push_back(issue);
  }

  const std::size_t no_candidate = std::numeric_limits<std::size_t>::max();
  std::vector<std::size_t> preferred(boundary.size(), no_candidate);
  for (std::size_t index = 0; index < boundary.size(); ++index) {
    if (ambiguous[index] != 0 || accepted[index].empty()) continue;
    preferred[index] = accepted[index][0];
  }

  for (std::size_t candidate_index = 0; candidate_index < candidates.size();
       ++candidate_index) {
    const heal_detail::Candidate& candidate = candidates[candidate_index];
    if (ambiguous[candidate.keep_position] != 0 ||
        ambiguous[candidate.drop_position] != 0) {
      continue;
    }
    if (preferred[candidate.keep_position] != candidate_index ||
        preferred[candidate.drop_position] != candidate_index) {
      continue;
    }
    if (candidate.pair.same_direction) {
      HealIssue issue;
      issue.kind = HealIssueKind::inconsistent_orientation;
      issue.location = heal_detail::edge_midpoint(brep, candidate.pair.keep);
      issue.gap = candidate.pair.gap;
      issue.entity_a = candidate.pair.keep.index();
      issue.entity_b = candidate.pair.drop.index();
      plan.issues.push_back(issue);
      continue;
    }
    plan.pairs.push_back(candidate.pair);
  }

  std::size_t nearest_index = no_candidate;
  for (std::size_t index = 0; index < boundary.size(); ++index) {
    if (ambiguous[index] != 0 || !accepted[index].empty()) continue;
    for (const std::size_t candidate_index : touching[index]) {
      if (nearest_index == no_candidate ||
          heal_detail::gap_issue_less(candidates[candidate_index],
                                      candidates[nearest_index])) {
        nearest_index = candidate_index;
      }
    }
  }
  if (nearest_index != no_candidate) {
    const heal_detail::Candidate& candidate = candidates[nearest_index];
    if (candidate.pair.gap > max_gap) {
      const std::size_t weld_index =
          heal_detail::selected_weld(brep, candidate.pair);
      const HealVertexWeld& weld = candidate.pair.welds[weld_index];
      HealIssue issue;
      issue.kind = HealIssueKind::gap_too_large;
      issue.location = heal_detail::midpoint(
          brep.vertex(weld.keep).point, brep.vertex(weld.drop).point);
      issue.gap = candidate.pair.gap;
      issue.entity_a = candidate.pair.keep.index();
      issue.entity_b = candidate.pair.drop.index();
      plan.issues.push_back(issue);
    }
  }

  return plan;
}

// Edges of `solid` that are used by exactly one coedge of the solid's faces, in
// ascending edge-index order. An empty result means the solid is closed (or the
// handle is unknown). Coedges of faces outside the solid are ignored.
inline std::vector<EdgeId> heal_boundary_edges(const Brep& brep, SolidId solid) {
  if (!brep.solid_alive(solid)) return {};
  std::vector<std::uint32_t> uses(brep.edges().size(), 0u);
  for (const ShellId shell : brep.solid_shells(solid)) {
    if (!brep.shell_alive(shell)) continue;
    for (const FaceId face_id : brep.shell_faces(shell)) {
      if (!brep.face_alive(face_id)) continue;
      LoopId loop_id = brep.face(face_id).first_loop;
      std::size_t loops_seen = 0;
      while (brep.loop_alive(loop_id) && loops_seen <= brep.loop_count()) {
        const Loop& loop = brep.loop(loop_id);
        if (brep.coedge_alive(loop.start)) {
          CoedgeId current = loop.start;
          std::size_t coedges_seen = 0;
          while (brep.coedge_alive(current) && coedges_seen <= brep.coedge_count()) {
            const Coedge& coedge = brep.coedge(current);
            if (brep.edge_alive(coedge.edge) && coedge.edge.index() < uses.size()) {
              ++uses[coedge.edge.index()];
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
  std::vector<EdgeId> boundary;
  for (const EdgeId id : brep.edge_ids()) {
    if (id.index() < uses.size() && uses[id.index()] == 1u) boundary.push_back(id);
  }
  return boundary;
}

}  // namespace detail

inline Result<HealResult> heal_brep(const Brep& brep, SolidId solid,
                                    const HealOptions& options) {
  if (!std::isfinite(options.tolerance) || options.tolerance < 0.0) {
    return Error{Errc::tolerance_violation,
                 "heal_brep: tolerance must be finite and non-negative", 0, 0};
  }
  if (!std::isfinite(options.max_gap) || !(options.max_gap > 0.0)) {
    return Error{Errc::tolerance_violation,
                 "heal_brep: max_gap must be finite and positive", 0, 0};
  }
  if (!std::isfinite(options.max_tolerance) || !(options.max_tolerance > 0.0)) {
    return Error{Errc::tolerance_violation,
                 "heal_brep: max_tolerance must be finite and positive", 0, 0};
  }
  if (!brep.solid_alive(solid)) {
    return Error{Errc::index_out_of_range, "heal_brep: unknown solid", 0, 0};
  }

  const ValidationReport topology = validate_topology(brep);
  if (!topology.ok()) {
    const ValidationIssue& issue = topology.issues.front();
    return Error{issue.code, issue.message, issue.entity, 0};
  }

  const std::vector<EdgeId> boundary = detail::heal_boundary_edges(brep, solid);

  // A boundary edge is only safe to stitch when both of its uses belong to the
  // target solid. If the other use is on another solid, the edge is shared and
  // heal_plan would otherwise pair the target's face with a face it must not
  // modify.
  std::vector<char> solid_faces(brep.faces().size(), 0);
  for (const ShellId shell : brep.solid_shells(solid)) {
    if (!brep.shell_alive(shell)) continue;
    for (const FaceId face_id : brep.shell_faces(shell)) {
      if (brep.face_alive(face_id) && face_id.index() < solid_faces.size()) {
        solid_faces[face_id.index()] = 1;
      }
    }
  }
  std::vector<std::vector<CoedgeId>> edge_coedges(brep.edges().size());
  for (const CoedgeId coedge_id : brep.coedge_ids()) {
    const Coedge& coedge = brep.coedge(coedge_id);
    if (brep.edge_alive(coedge.edge) &&
        coedge.edge.index() < edge_coedges.size()) {
      edge_coedges[coedge.edge.index()].push_back(coedge_id);
    }
  }
  for (const EdgeId edge_id : boundary) {
    for (const CoedgeId coedge_id : edge_coedges[edge_id.index()]) {
      const FaceId face_id = brep.coedge(coedge_id).face;
      if (!brep.face_alive(face_id) || face_id.index() >= solid_faces.size() ||
          solid_faces[face_id.index()] == 0) {
        return Error{Errc::invalid_topology,
                     "heal_brep: boundary edge is shared with another solid", 0,
                     0};
      }
    }
  }

  const auto& surfaces = brep.geometry().surfaces;
  const auto& curves = brep.geometry().curves;
  for (const ShellId shell : brep.solid_shells(solid)) {
    if (!brep.shell_alive(shell)) continue;
    for (const FaceId face_id : brep.shell_faces(shell)) {
      if (!brep.face_alive(face_id)) continue;
      const Face& face = brep.face(face_id);
      if (face.surface == kNoGeometry || face.surface >= surfaces.size() ||
          std::get_if<Plane3>(&surfaces[face.surface]) == nullptr) {
        return Error{Errc::unsupported,
                     "heal_brep: curved geometry is not supported yet", 0, 0};
      }
      LoopId loop_id = face.first_loop;
      std::size_t loops_seen = 0;
      while (brep.loop_alive(loop_id) && loops_seen <= brep.loop_count()) {
        const Loop& loop = brep.loop(loop_id);
        if (brep.coedge_alive(loop.start)) {
          CoedgeId current = loop.start;
          std::size_t coedges_seen = 0;
          while (brep.coedge_alive(current) &&
                 coedges_seen <= brep.coedge_count()) {
            const Coedge& coedge = brep.coedge(current);
            if (brep.edge_alive(coedge.edge)) {
              const Edge& edge = brep.edge(coedge.edge);
              if (edge.curve != kNoGeometry &&
                  (edge.curve >= curves.size() ||
                   std::get_if<LineCurve3>(&curves[edge.curve]) == nullptr)) {
                return Error{Errc::unsupported,
                             "heal_brep: curved geometry is not supported yet",
                             0, 0};
              }
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

  HealReport report;
  report.boundary_edges = static_cast<std::uint32_t>(boundary.size());
  const std::int64_t euler_before = euler_characteristic(brep, solid);

  Result<detail::HealPlan> plan_result =
      detail::heal_plan(brep, boundary, options.max_gap);
  if (!plan_result) return plan_result.error();
  detail::HealPlan plan = std::move(plan_result).value();

  auto reject = [&](const char* message, bool carry_plan) -> Result<HealResult> {
    HealReport rejected;
    rejected.boundary_edges = report.boundary_edges;
    HealIssue issue;
    issue.kind = HealIssueKind::invalid_result;
    issue.message = message;
    rejected.issues.push_back(issue);
    if (carry_plan) {
      rejected.issues.insert(rejected.issues.end(), plan.issues.begin(),
                             plan.issues.end());
    }
    return HealResult{brep, std::move(rejected)};
  };

  if (plan.pairs.empty()) {
    report.issues = plan.issues;
    if (!boundary.empty() && report.issues.empty()) {
      HealIssue issue;
      issue.kind = HealIssueKind::invalid_result;
      issue.message = "heal_brep: no mutual pair";
      report.issues.push_back(issue);
    }
    return HealResult{brep, std::move(report)};
  }

  double largest_weld_tolerance = 0.0;
  for (const detail::HealEdgePair& pair : plan.pairs) {
    for (const detail::HealVertexWeld& weld : pair.welds) {
      const double weld_tolerance =
          std::max({brep.vertex(weld.keep).tol, brep.vertex(weld.drop).tol,
                    pair.gap});
      largest_weld_tolerance = std::max(largest_weld_tolerance, weld_tolerance);
      if (weld_tolerance > options.max_tolerance) {
        return reject("heal_brep: tolerance cap exceeded", false);
      }
    }
  }

  const Result<std::vector<detail::HealVertexWeld>> welds_result =
      detail::heal_dedup_welds(brep, plan.pairs);
  if (!welds_result) return welds_result.error();

  SolidEditor editor(brep);
  const Result<std::vector<detail::HealVertexWeld>> applied_result =
      detail::heal_apply_welds(editor, welds_result.value());
  if (!applied_result) return reject(applied_result.error().message, true);
  const std::vector<detail::HealVertexWeld>& applied_welds =
      applied_result.value();

  std::uint32_t sewn_edges = 0;
  for (const detail::HealEdgePair& pair : plan.pairs) {
    if (!editor.brep().edge_alive(pair.keep) ||
        !editor.brep().edge_alive(pair.drop)) {
      return reject("heal_brep: an edge disappeared during stitching", true);
    }
    const Result<void> merged = editor.merge_edges(pair.keep, pair.drop);
    if (!merged) return reject(merged.error().message, true);
    ++sewn_edges;
  }

  std::vector<char> welded_vertices(editor.brep().vertices().size(), 0);
  for (const detail::HealVertexWeld& weld : applied_welds) {
    if (editor.brep().vertex_alive(weld.keep) &&
        weld.keep.index() < welded_vertices.size()) {
      welded_vertices[weld.keep.index()] = 1;
    }
  }
  std::vector<char> incident(editor.brep().edges().size(), 0);
  for (const EdgeId edge_id : editor.brep().edge_ids()) {
    const Edge& edge = editor.brep().edge(edge_id);
    if (welded_vertices[edge.v0.index()] ||
        welded_vertices[edge.v1.index()]) {
      incident[edge_id.index()] = 1;
    }
  }
  for (const CoedgeId coedge_id : editor.brep().coedge_ids()) {
    const Coedge& coedge = editor.brep().coedge(coedge_id);
    if (coedge.edge.index() >= incident.size() ||
        incident[coedge.edge.index()] == 0) {
      continue;
    }
    const Face& face = editor.brep().face(coedge.face);
    const Result<detail::FaceUvProjection> projection =
        detail::face_uv_projection(editor.brep(), face, options.tolerance);
    if (!projection) return reject(projection.error().message, true);
    const EdgeId edge_id = coedge.edge;
    const Edge& edge = editor.brep().edge(edge_id);
    const Result<Point2d> uv0 = detail::face_uv_project(
        projection.value(), editor.brep().vertex(edge.v0).point, 0.0, false);
    if (!uv0) return reject(uv0.error().message, true);
    const Result<Point2d> uv1 = detail::face_uv_project(
        projection.value(), editor.brep().vertex(edge.v1).point, 0.0, false);
    if (!uv1) return reject(uv1.error().message, true);
    const std::uint32_t pcurve =
        editor.add_pcurve(Pcurve2d{uv0.value(), uv1.value()});
    const Result<void> set_pcurve =
        editor.set_coedge_pcurve(coedge_id, pcurve);
    if (!set_pcurve) return reject(set_pcurve.error().message, true);
    const Result<void> set_geometry =
        editor.set_edge_geometry(edge_id, kNoGeometry, 0.0, 1.0);
    if (!set_geometry) return reject(set_geometry.error().message, true);
  }

  for (const ShellId shell : editor.brep().solid_shells(solid)) {
    if (!editor.brep().shell_alive(shell)) continue;
    if (detail::heal_detail::shell_is_closed_by_uses(editor.brep(), shell)) {
      const Result<void> closed = editor.set_shell_closed(shell, true);
      if (!closed) return reject(closed.error().message, true);
    }
  }

  const Result<double> raised_result =
      editor.raise_tolerances(options.max_tolerance);
  if (!raised_result) return reject(raised_result.error().message, true);
  const double raised = raised_result.value();

  if (!detail::heal_boundary_edges(editor.brep(), solid).empty()) {
    return reject("heal_brep: a shell is still open after stitching", true);
  }

  if (euler_characteristic(editor.brep(), solid) != euler_before) {
    return reject("heal_brep: Euler characteristic changed", false);
  }

  Result<Brep> committed = editor.commit();
  if (!committed) return reject(committed.error().message, true);

  report.welded_vertices = detail::heal_live_weld_count(editor.brep(), applied_welds);
  report.sewn_edges = sewn_edges;
  report.changed = report.welded_vertices + report.sewn_edges > 0;
  double max_closed_gap = 0.0;
  for (const detail::HealEdgePair& pair : plan.pairs) {
    max_closed_gap = std::max(max_closed_gap, pair.gap);
  }
  report.max_closed_gap = max_closed_gap;
  report.max_tolerance = std::max(raised, largest_weld_tolerance);
  report.issues = std::move(plan.issues);
  return HealResult{std::move(committed).value(), std::move(report)};
}

}  // namespace v1
}  // namespace geom
