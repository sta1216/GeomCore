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
  // Best candidate pair for every boundary edge, including edges whose pairing
  // stayed ambiguous. Geometric validation (coincident carriers, on-surface
  // rims) uses this so a corrupt input is rejected even when the tie-break
  // cannot pick a seam.
  std::vector<HealEdgePair> considered;
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

// Whether the face using `coedge` carries a Plane3 surface. Missing geometry is
// reported as planar so the tie-break leaves the pre-existing ambiguity
// behavior untouched when a carrier is unknown.
inline bool adjacent_face_is_planar(const Brep& brep, CoedgeId coedge_id) {
  if (!brep.coedge_alive(coedge_id)) return true;
  const FaceId face_id = brep.coedge(coedge_id).face;
  if (!brep.face_alive(face_id)) return true;
  const Face& face = brep.face(face_id);
  const auto& surfaces = brep.geometry().surfaces;
  if (face.surface == kNoGeometry || face.surface >= surfaces.size()) {
    return true;
  }
  return std::holds_alternative<Plane3>(surfaces[face.surface]);
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
    const double tie_gap = kAmbiguityGap * std::max(1.0, best_gap);
    if (std::fabs(best_gap - next_gap) > tie_gap) continue;

    // Two gap edges on a central-symmetric analytic carrier (the two rims of
    // one cylinder, two latitudes of one sphere) have identical endpoint
    // distances, so the distance mapping alone cannot tell the real seam from
    // the central-symmetric decoy. Only the real seam welds the two coedge
    // starts to each other, which is what directions_agree reports. Resolve
    // such a tie by orientation only when the boundary edge's adjacent face is
    // not planar: planar rims must keep reporting ambiguous_gap even when one
    // candidate happens to be orientation-consistent.
    if (!heal_detail::adjacent_face_is_planar(
            brep, first_coedge[boundary[index].index()])) {
      std::size_t tied_opposite = 0;
      std::size_t tied_opposite_index = std::numeric_limits<std::size_t>::max();
      for (const std::size_t candidate_index : accepted[index]) {
        if (std::fabs(candidates[candidate_index].pair.gap - best_gap) >
            tie_gap) {
          continue;
        }
        if (!candidates[candidate_index].pair.same_direction) {
          ++tied_opposite;
          tied_opposite_index = candidate_index;
        }
      }
      if (tied_opposite == 1u) {
        accepted[index].assign(1u, tied_opposite_index);
        continue;
      }
    }

    ambiguous[index] = 1;
    const EdgeId edge_id = boundary[index];
    HealIssue issue;
    issue.kind = HealIssueKind::ambiguous_gap;
    issue.location = heal_detail::edge_midpoint(brep, edge_id);
    issue.gap = best_gap;
    issue.entity_a = edge_id.index();
    plan.issues.push_back(issue);
  }

  // Record the best candidate for every boundary edge, ambiguous or not, so the
  // geometric guards can reject a bad carrier before the ambiguity is reported.
  for (std::size_t index = 0; index < boundary.size(); ++index) {
    if (accepted[index].empty()) continue;
    const HealEdgePair& best = candidates[accepted[index][0]].pair;
    bool seen = false;
    for (const HealEdgePair& known : plan.considered) {
      if (known.keep == best.keep && known.drop == best.drop) {
        seen = true;
        break;
      }
    }
    if (!seen) plan.considered.push_back(best);
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

// Largest radius/centre/origin magnitude of a surface. Multiplied by 1e-12 it
// supplies the geometric floor for the coincidence and on-surface guards.
inline double heal_surface_scale(const SurfaceGeometry& surface) {
  double scale = 0.0;
  const auto include_point = [&scale](const Point3d& point) {
    scale = std::max(scale, std::hypot(point[0], point[1], point[2]));
  };
  if (const auto* plane = std::get_if<Plane3>(&surface)) {
    scale = std::max(scale,
                     std::hypot(plane->normal[0], plane->normal[1],
                                plane->normal[2]));
    scale = std::max(scale, std::fabs(plane->offset));
  } else if (const auto* cylinder = std::get_if<Cylinder3>(&surface)) {
    include_point(cylinder->origin);
    scale = std::max(scale, cylinder->radius);
  } else if (const auto* sphere = std::get_if<Sphere3>(&surface)) {
    include_point(sphere->center);
    scale = std::max(scale, sphere->radius);
  } else if (const auto* cone = std::get_if<Cone3>(&surface)) {
    include_point(cone->origin);
    scale = std::max(scale, cone->radius);
  } else if (const auto* torus = std::get_if<Torus3>(&surface)) {
    include_point(torus->center);
    scale = std::max(scale, torus->major_radius);
    scale = std::max(scale, torus->minor_radius);
  }
  return scale;
}

// True when `a` and `b` describe the same carrier surface to `tolerance`.
// Planes are always compatible: planar healing never depended on their exact
// placement, and the planar tests rely on that. Mixed kinds and NURBS are not
// coincident.
inline bool heal_surfaces_coincide(const SurfaceGeometry& a,
                                   const SurfaceGeometry& b,
                                   double tolerance) {
  const bool a_plane = std::holds_alternative<Plane3>(a);
  const bool b_plane = std::holds_alternative<Plane3>(b);
  if (a_plane && b_plane) return true;

  if (const auto* cylinder_a = std::get_if<Cylinder3>(&a)) {
    const auto* cylinder_b = std::get_if<Cylinder3>(&b);
    if (cylinder_b == nullptr) return false;
    if (std::fabs(cylinder_a->radius - cylinder_b->radius) > tolerance) {
      return false;
    }
    if (norm2(cross(cylinder_a->axis, cylinder_b->axis)) > tolerance) {
      return false;
    }
    const Vec3d delta = cylinder_a->origin - cylinder_b->origin;
    const Vec3d offset =
        delta - cylinder_b->axis * dot(delta, cylinder_b->axis);
    return norm2(offset) <= tolerance;
  }

  if (const auto* sphere_a = std::get_if<Sphere3>(&a)) {
    const auto* sphere_b = std::get_if<Sphere3>(&b);
    if (sphere_b == nullptr) return false;
    return std::fabs(sphere_a->radius - sphere_b->radius) <= tolerance &&
           norm2(sphere_a->center - sphere_b->center) <= tolerance;
  }

  if (const auto* cone_a = std::get_if<Cone3>(&a)) {
    const auto* cone_b = std::get_if<Cone3>(&b);
    if (cone_b == nullptr) return false;
    if (std::fabs(cone_a->radius - cone_b->radius) > tolerance ||
        std::fabs(cone_a->half_angle - cone_b->half_angle) > tolerance ||
        norm2(cross(cone_a->axis, cone_b->axis)) > tolerance) {
      return false;
    }
    const double tan_a = std::tan(cone_a->half_angle);
    const double tan_b = std::tan(cone_b->half_angle);
    if (!(tan_a > 0.0) || !(tan_b > 0.0)) return false;
    const Point3d apex_a =
        cone_a->origin + cone_a->axis * (cone_a->radius / tan_a);
    const Point3d apex_b =
        cone_b->origin + cone_b->axis * (cone_b->radius / tan_b);
    return norm2(apex_a - apex_b) <= tolerance;
  }

  if (const auto* torus_a = std::get_if<Torus3>(&a)) {
    const auto* torus_b = std::get_if<Torus3>(&b);
    if (torus_b == nullptr) return false;
    return norm2(cross(torus_a->axis, torus_b->axis)) <= tolerance &&
           norm2(torus_a->center - torus_b->center) <= tolerance &&
           std::fabs(torus_a->major_radius - torus_b->major_radius) <=
               tolerance &&
           std::fabs(torus_a->minor_radius - torus_b->minor_radius) <=
               tolerance;
  }

  return false;
}

// Distance from `point` to the implicit surface of `surface`. An unsupported
// kind returns infinity; heal_brep rejects those kinds before calling this.
inline double heal_surface_residual(const SurfaceGeometry& surface,
                                    const Point3d& point) {
  if (const auto* plane = std::get_if<Plane3>(&surface)) {
    return std::fabs(plane->signed_distance(point));
  }
  if (const auto* cylinder = std::get_if<Cylinder3>(&surface)) {
    const Vec3d delta = point - cylinder->origin;
    const Vec3d radial = delta - cylinder->axis * dot(delta, cylinder->axis);
    return std::fabs(norm2(radial) - cylinder->radius);
  }
  if (const auto* sphere = std::get_if<Sphere3>(&surface)) {
    return std::fabs(norm2(point - sphere->center) - sphere->radius);
  }
  if (const auto* cone = std::get_if<Cone3>(&surface)) {
    const double tan_half = std::tan(cone->half_angle);
    const Vec3d delta = point - cone->origin;
    const double axial = dot(delta, cone->axis);
    const Vec3d radial = delta - cone->axis * axial;
    const double expected = cone->radius - axial * tan_half;
    return std::fabs(norm2(radial) - expected) /
           std::sqrt(1.0 + tan_half * tan_half);
  }
  if (const auto* torus = std::get_if<Torus3>(&surface)) {
    const Vec3d delta = point - torus->center;
    const double axial = dot(delta, torus->axis);
    const Vec3d radial = delta - torus->axis * axial;
    const double ring = norm2(radial) - torus->major_radius;
    return std::fabs(std::hypot(ring, axial) - torus->minor_radius);
  }
  return std::numeric_limits<double>::infinity();
}

// Distance from `point` to a line or circle curve; infinity for a kind the
// caller does not understand (NURBS is rejected by the heal_brep guards).
inline double heal_curve_residual(const CurveGeometry& curve,
                                  const Point3d& point) {
  if (const auto* line = std::get_if<LineCurve3>(&curve)) {
    const Vec3d offset = point - line->origin;
    return norm2(offset - line->direction * dot(offset, line->direction));
  }
  if (const auto* circle = std::get_if<CircleCurve3>(&curve)) {
    const Vec3d normal = cross(circle->x_axis, circle->y_axis);
    const Vec3d delta = point - circle->center;
    const double axial = dot(delta, normal);
    const Vec3d radial = delta - normal * axial;
    return std::hypot(norm2(radial) - circle->radius, axial);
  }
  return std::numeric_limits<double>::infinity();
}

// Samples the edge curve at t0, the parameter midpoint and t1 and requires every
// sample to sit on `surface` within `tolerance`. Edges without geometry, and
// curve kinds without an on-surface test here, pass through.
inline bool heal_edge_lies_on_surface(const Brep& brep, const Edge& edge,
                                      const SurfaceGeometry& surface,
                                      double tolerance) {
  if (edge.curve == kNoGeometry ||
      edge.curve >= brep.geometry().curves.size()) {
    return true;
  }
  const CurveGeometry& curve = brep.geometry().curves[edge.curve];
  for (int sample = 0; sample < 3; ++sample) {
    const double fraction = 0.5 * static_cast<double>(sample);
    const double t = edge.t0 + (edge.t1 - edge.t0) * fraction;
    Point3d point;
    if (const auto* line = std::get_if<LineCurve3>(&curve)) {
      point = line->eval(t);
    } else if (const auto* circle = std::get_if<CircleCurve3>(&curve)) {
      point = circle->eval(t);
    } else {
      return true;
    }
    if (!(heal_surface_residual(surface, point) <= tolerance)) return false;
  }
  return true;
}

// Samples the 3D path of `coedge` in traversal order, from the coedge's start
// vertex to its end vertex. The endpoint samples come from the vertices, so a
// vertex that was just welded projects at its new position. Line and circle
// carriers are supported; a circle segment spans at most a tenth of a turn and
// the sample count is capped so a closed rim stays bounded.
inline Result<std::vector<Point3d>> heal_coedge_path(const Brep& brep,
                                                     const Coedge& coedge) {
  if (!brep.edge_alive(coedge.edge)) {
    return Error{Errc::invalid_topology, "heal_brep: coedge edge is not alive",
                 coedge.edge.index(), 0u};
  }
  const Edge& edge = brep.edge(coedge.edge);
  const Point3d start = brep.vertex(coedge.reversed ? edge.v1 : edge.v0).point;
  const Point3d end = brep.vertex(coedge.reversed ? edge.v0 : edge.v1).point;
  std::vector<Point3d> path{start, end};
  if (edge.curve == kNoGeometry) return path;
  if (edge.curve >= brep.geometry().curves.size()) {
    return Error{Errc::invalid_topology,
                 "heal_brep: edge curve index is out of range", edge.curve, 0u};
  }
  const CurveGeometry& curve = brep.geometry().curves[edge.curve];
  if (std::get_if<LineCurve3>(&curve) != nullptr) return path;
  const auto* circle = std::get_if<CircleCurve3>(&curve);
  if (circle == nullptr) {
    return Error{Errc::unsupported, "heal_brep: curve type is not supported yet",
                 edge.curve, 0u};
  }

  double span = coedge.reversed ? edge.t0 - edge.t1 : edge.t1 - edge.t0;
  if (edge.v0 == edge.v1) {
    span = span < 0.0 ? -2.0 * kFaceUvPi : 2.0 * kFaceUvPi;
  } else {
    while (span > 2.0 * kFaceUvPi) span -= 2.0 * kFaceUvPi;
    while (span <= -2.0 * kFaceUvPi) span += 2.0 * kFaceUvPi;
  }
  if (span == 0.0) return path;

  std::size_t count =
      static_cast<std::size_t>(std::ceil(std::fabs(span) / 0.1)) + 1u;
  count = std::max<std::size_t>(count, 3u);
  count = std::min<std::size_t>(count, 257u);
  const double first = coedge.reversed ? edge.t1 : edge.t0;
  std::vector<Point3d> samples;
  samples.reserve(count);
  samples.push_back(start);
  for (std::size_t i = 1; i + 1 < count; ++i) {
    const double fraction =
        static_cast<double>(i) / static_cast<double>(count - 1);
    samples.push_back(circle->eval(first + span * fraction));
  }
  samples.push_back(end);
  return samples;
}

// True when `point` is the singular pole of a sphere or the apex of a cone,
// where the azimuth is undefined and the projection has to borrow a neighbour's.
inline bool heal_is_singular_sample(const SurfaceGeometry& surface,
                                    const Point3d& point) {
  if (const auto* sphere = std::get_if<Sphere3>(&surface)) {
    const double sin_v =
        std::clamp(dot(point - sphere->center, sphere->axis) / sphere->radius,
                   -1.0, 1.0);
    return std::fabs(std::fabs(sin_v) - 1.0) <= 1e-12;
  }
  if (const auto* cone = std::get_if<Cone3>(&surface)) {
    const double v_apex = cone->radius / std::tan(cone->half_angle);
    const Point3d apex = cone->origin + cone->axis * v_apex;
    const double scale = std::max({1.0, std::fabs(v_apex), cone->radius});
    return norm2(point - apex) <= 1e-9 * scale;
  }
  return false;
}

// Adds a carrier that traverses `curve` backwards over the same interval:
// eval_reversed(t) = eval(t0 + t1 - t). `merge_edges` copies a dropped edge's
// carrier into a curve-less keep edge without swapping t0/t1, so a flipped
// welded correspondence would otherwise leave the keep edge running against its
// own direction.
inline Result<std::uint32_t> heal_reverse_curve(SolidEditor& editor,
                                                EdgeId edge_id) {
  const Edge& edge = editor.brep().edge(edge_id);
  if (edge.curve == kNoGeometry ||
      edge.curve >= editor.brep().geometry().curves.size()) {
    return Error{Errc::invalid_topology,
                 "heal_brep: edge curve index is out of range", edge_id.index(),
                 0u};
  }
  const CurveGeometry& curve = editor.brep().geometry().curves[edge.curve];
  const double sum = edge.t0 + edge.t1;
  CurveGeometry reversed;
  if (const auto* line = std::get_if<LineCurve3>(&curve)) {
    reversed = LineCurve3{line->origin + line->direction * sum,
                          line->direction * -1.0};
  } else if (const auto* circle = std::get_if<CircleCurve3>(&curve)) {
    const double c = std::cos(sum);
    const double s = std::sin(sum);
    reversed = CircleCurve3{circle->center,
                            circle->x_axis * c + circle->y_axis * s,
                            circle->x_axis * s - circle->y_axis * c,
                            circle->radius};
  } else {
    return Error{Errc::unsupported, "heal_brep: curve type is not supported yet",
                 edge.curve, 0u};
  }
  return editor.add_curve(std::move(reversed));
}

// Rebuilds every pcurve of `face` by walking its loops in coedge order and
// projecting one polyline per coedge with the previous parameter as the unwrap
// hint. The edge carriers are never touched. All loop runs are projected and
// simplified before the first pcurve is written, but SolidEditor has no
// transaction: a write failure leaves the editor partially updated and the
// caller must propagate the error without committing.
inline Result<void> heal_rebuild_face_pcurves(SolidEditor& editor, FaceId face_id,
                                              double tolerance) {
  const Brep& brep = editor.brep();
  if (!brep.face_alive(face_id)) {
    return Error{Errc::invalid_topology, "heal_brep: face is not alive",
                 face_id.index(), 0u};
  }
  const Result<FaceUvProjection> projection_result =
      face_uv_projection(brep, brep.face(face_id), tolerance);
  if (!projection_result) return projection_result.error();
  const FaceUvProjection& projection = projection_result.value();

  struct Run {
    CoedgeId coedge;
    bool curved = false;
    std::vector<Point2d> points;
    // The 3D samples in the same (stored) order as `points`, used to recognise
    // the sphere poles and the cone apex after the run has been oriented.
    std::vector<Point3d> samples;
  };
  std::vector<Run> runs;
  double extent = 1.0;

  const std::size_t coedge_limit = brep.coedge_count() + 1u;
  const std::size_t loop_limit = brep.loop_count() + 1u;
  std::size_t loop_steps = 0u;
  for (LoopId loop_id = brep.face(face_id).first_loop; !loop_id.is_null();
       loop_id = brep.loop(loop_id).next) {
    if (!brep.loop_alive(loop_id)) {
      return Error{Errc::invalid_topology,
                   "heal_brep: face uses a loop that is not alive",
                   loop_id.index(), 0u};
    }
    const Loop& loop = brep.loop(loop_id);
    if (loop.face != face_id) {
      return Error{Errc::invalid_topology,
                   "heal_brep: loop belongs to a different face",
                   loop_id.index(), 0u};
    }
    if (++loop_steps > loop_limit) {
      return Error{Errc::invalid_topology,
                   "heal_brep: face loop chain does not terminate",
                   loop_id.index(), 0u};
    }

    double previous_u = 0.0;
    double previous_v = 0.0;
    Point3d previous_sample{};
    bool has_previous = false;
    bool closed = false;
    const std::size_t loop_first_run = runs.size();
    CoedgeId coedge_id = loop.start;
    for (std::size_t step = 0; step <= coedge_limit; ++step) {
      if (!brep.coedge_alive(coedge_id)) {
        return Error{Errc::invalid_topology,
                     "heal_brep: loop uses a coedge that is not alive",
                     coedge_id.index(), 0u};
      }
      const Coedge& coedge = brep.coedge(coedge_id);
      if (coedge.face != face_id) {
        return Error{Errc::invalid_topology,
                     "heal_brep: coedge belongs to a different face",
                     coedge_id.index(), 0u};
      }
      const Result<std::vector<Point3d>> path_result =
          heal_coedge_path(brep, coedge);
      if (!path_result) return path_result.error();
      const std::vector<Point3d>& path = path_result.value();

      Run run;
      run.coedge = coedge_id;
      const Edge& edge = brep.edge(coedge.edge);
      run.curved = edge.curve != kNoGeometry &&
                   edge.curve < brep.geometry().curves.size() &&
                   std::get_if<CircleCurve3>(&brep.geometry().curves[edge.curve]) !=
                       nullptr;

      run.points.reserve(path.size());
      for (std::size_t i = 0; i < path.size(); ++i) {
        Result<Point2d> projected =
            face_uv_project(projection, path[i], previous_u, has_previous);
        if (!projected && has_previous) {
          // A pole or apex is singular at the exact sample; step a hair back
          // along the path and keep the previous parameter.
          const Point3d nudged =
              previous_sample + (path[i] - previous_sample) * (1.0 - 1e-9);
          projected = face_uv_project(projection, nudged, previous_u, true);
        }
        if (!projected) return projected.error();
        double u = projected.value()[0];
        double v = projected.value()[1];
        if (has_previous) {
          if (projection.periodic_u()) u = face_uv_unwrap(u, previous_u);
          if (projection.periodic_v()) v = face_uv_unwrap(v, previous_v);
        }
        previous_u = u;
        previous_v = v;
        previous_sample = path[i];
        has_previous = true;
        run.points.push_back(Point2d{u, v});
        run.samples.push_back(path[i]);
        extent = std::max({extent, std::fabs(u), std::fabs(v)});
      }
      if (coedge.reversed) {
        std::reverse(run.points.begin(), run.points.end());
        std::reverse(run.samples.begin(), run.samples.end());
      }
      runs.push_back(std::move(run));

      if (coedge.next == loop.start) {
        closed = true;
        break;
      }
      coedge_id = coedge.next;
    }
    if (!closed) {
      return Error{Errc::invalid_topology,
                   "heal_brep: loop coedge chain does not terminate",
                   loop.start.index(), 0u};
    }
    // A seam that retraces itself across the periodic cut ends the loop one
    // full turn away from its start, which leaves the stored uv polygon cut
    // open (a cone apex fan then integrates to zero). Shift the closing
    // constant-u run onto the loop's first column so the chart closes.
    if (projection.periodic_u() && runs.size() > loop_first_run) {
      const std::vector<Point2d>& first_points = runs[loop_first_run].points;
      const std::vector<Point2d>& last_points = runs.back().points;
      if (!first_points.empty() && !last_points.empty() &&
          std::fabs(last_points.front()[0] - last_points.back()[0]) <= 1e-9) {
        const double delta = last_points.back()[0] - first_points.front()[0];
        if (std::fabs(delta) > kFaceUvPi) {
          const double turns = std::round(delta / (2.0 * kFaceUvPi));
          for (Point2d& point : runs.back().points) {
            point[0] -= turns * 2.0 * kFaceUvPi;
          }
        }
      }
    }
  }
  if (runs.empty()) {
    return Error{Errc::invalid_topology, "heal_brep: face has no loops",
                 face_id.index(), 0u};
  }

  // At a sphere pole or a cone apex the azimuth is undefined, so the projection
  // returns an arbitrary value (usually the raw atan2 of the rounded sample). A
  // meridian or seam that ends at the pole would then leave the chart column of
  // its own carrier. Borrow the azimuth from the adjacent sample in stored
  // order: the previous point when there is one, otherwise the next point (the
  // same value a probe a hair inside the path would project). The latitude stays
  // exactly at the pole/apex.
  const Face& face = brep.face(face_id);
  if (face.surface != kNoGeometry &&
      face.surface < brep.geometry().surfaces.size()) {
    const SurfaceGeometry& surface = brep.geometry().surfaces[face.surface];
    for (Run& run : runs) {
      for (std::size_t i = 0; i < run.points.size(); ++i) {
        if (!heal_is_singular_sample(surface, run.samples[i])) continue;
        double u = run.points[i][0];
        if (i > 0 && !heal_is_singular_sample(surface, run.samples[i - 1])) {
          u = run.points[i - 1][0];
        } else if (i + 1 < run.points.size() &&
                   !heal_is_singular_sample(surface, run.samples[i + 1])) {
          u = run.points[i + 1][0];
        }
        run.points[i][0] = u;
      }
    }
  }

  const double cross_tolerance = 1e-12 * extent * extent;
  for (Run& run : runs) {
    if (run.points.size() < 3u) continue;
    std::vector<Point2d> simplified;
    simplified.reserve(run.points.size());
    simplified.push_back(run.points.front());
    for (std::size_t i = 1; i + 1 < run.points.size(); ++i) {
      const Point2d& a = simplified.back();
      const Point2d& b = run.points[i];
      const Point2d& c = run.points[i + 1];
      const double cross = (b[0] - a[0]) * (c[1] - a[1]) -
                           (b[1] - a[1]) * (c[0] - a[0]);
      const double between =
          (b[0] - a[0]) * (c[0] - b[0]) + (b[1] - a[1]) * (c[1] - b[1]);
      if (std::fabs(cross) <= cross_tolerance && between >= 0.0) continue;
      simplified.push_back(b);
    }
    simplified.push_back(run.points.back());
    // A circle keeps one interior sample so its arc direction stays defined and
    // the polyline can still be read as an arc.
    if (run.curved && simplified.size() < 3u) {
      simplified.insert(simplified.begin() + 1,
                        run.points[run.points.size() / 2u]);
    }
    run.points = std::move(simplified);
  }

  for (const Run& run : runs) {
    if (run.points.size() < 2u) {
      return Error{Errc::degenerate, "heal_brep: pcurve has fewer than two points",
                   run.coedge.index(), 0u};
    }
    Pcurve2d pcurve;
    pcurve.points = run.points;
    const std::uint32_t index = editor.add_pcurve(std::move(pcurve));
    const Result<void> stored = editor.set_coedge_pcurve(run.coedge, index);
    if (!stored) return stored.error();
  }
  return {};
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
  // target solid. A shared edge would let heal_plan read the foreign coedge for
  // the traversal direction, and solid-local closure counts would then describe
  // the wrong shell.
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
      if (face.surface == kNoGeometry || face.surface >= surfaces.size()) {
        return Error{Errc::invalid_topology, "heal_brep: face has no surface",
                     face_id.index(), 0};
      }
      const SurfaceGeometry& surface = surfaces[face.surface];
      if (!std::holds_alternative<Plane3>(surface) &&
          !std::holds_alternative<Cylinder3>(surface) &&
          !std::holds_alternative<Cone3>(surface) &&
          !std::holds_alternative<Sphere3>(surface) &&
          !std::holds_alternative<Torus3>(surface)) {
        return Error{Errc::unsupported,
                     "heal_brep: surface type is not supported yet",
                     face_id.index(), 0};
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
                   (!std::holds_alternative<LineCurve3>(curves[edge.curve]) &&
                    !std::holds_alternative<CircleCurve3>(curves[edge.curve])))) {
                return Error{Errc::unsupported,
                             "heal_brep: curve type is not supported yet", 0, 0};
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

  // The geometric guards below run on the immutable input after the plan is
  // known and before any SolidEditor mutation, so a rejection leaves the brep
  // exactly as it arrived.
  double surface_scale = 0.0;
  for (const SurfaceGeometry& surface : surfaces) {
    surface_scale = std::max(surface_scale, detail::heal_surface_scale(surface));
  }
  const double check_tolerance =
      std::max(options.tolerance, 1e-12 * surface_scale);

  // The first coedge use of a boundary edge belongs to the target solid (the
  // shared-solid check above guarantees it), so it names the face whose surface
  // the edge is expected to lie on.
  const auto edge_surface = [&](EdgeId edge_id) -> const SurfaceGeometry* {
    for (const CoedgeId coedge_id : edge_coedges[edge_id.index()]) {
      const FaceId face_id = brep.coedge(coedge_id).face;
      if (!brep.face_alive(face_id)) continue;
      const Face& face = brep.face(face_id);
      if (face.surface != kNoGeometry && face.surface < surfaces.size()) {
        return &surfaces[face.surface];
      }
    }
    return nullptr;
  };

  for (const detail::HealEdgePair& pair : plan.considered) {
    const SurfaceGeometry* keep_surface = edge_surface(pair.keep);
    const SurfaceGeometry* drop_surface = edge_surface(pair.drop);
    if (keep_surface == nullptr || drop_surface == nullptr) continue;
    if (!detail::heal_surfaces_coincide(*keep_surface, *drop_surface,
                                        check_tolerance)) {
      return Error{Errc::unsupported,
                   "heal_brep: gap faces do not share the same surface",
                   pair.keep.index(), pair.drop.index()};
    }
  }

  for (const detail::HealEdgePair& pair : plan.considered) {
    const SurfaceGeometry* keep_surface = edge_surface(pair.keep);
    const SurfaceGeometry* drop_surface = edge_surface(pair.drop);
    const Edge& keep_edge = brep.edge(pair.keep);
    const Edge& drop_edge = brep.edge(pair.drop);
    const bool keep_off_surface =
        keep_surface != nullptr &&
        !detail::heal_edge_lies_on_surface(brep, keep_edge, *keep_surface,
                                           check_tolerance);
    const bool drop_off_surface =
        drop_surface != nullptr &&
        !detail::heal_edge_lies_on_surface(brep, drop_edge, *drop_surface,
                                           check_tolerance);
    if (keep_off_surface || drop_off_surface) {
      return Error{Errc::unsupported,
                   "heal_brep: gap edge is not on its face surface",
                   pair.keep.index(), pair.drop.index()};
    }
  }

  // The kept edge's curve survives merge_edges, so after welding it must still
  // pass through both of the edge's final keep vertices. Resolve chains to the
  // single surviving vertex before measuring.
  const Result<std::vector<detail::HealVertexWeld>> welds_result =
      detail::heal_dedup_welds(brep, plan.pairs);
  if (!welds_result) return welds_result.error();

  std::vector<std::uint32_t> final_keep(brep.vertices().size());
  for (std::size_t index = 0; index < final_keep.size(); ++index) {
    final_keep[index] = static_cast<std::uint32_t>(index);
  }
  for (const detail::HealVertexWeld& weld : welds_result.value()) {
    if (weld.keep.index() < final_keep.size() &&
        weld.drop.index() < final_keep.size()) {
      final_keep[weld.drop.index()] = weld.keep.index();
    }
  }
  const auto root_of = [&final_keep](std::uint32_t index) {
    std::uint32_t root = index;
    while (final_keep[root] != root) {
      final_keep[root] = final_keep[final_keep[root]];
      root = final_keep[root];
    }
    return root;
  };

  for (const detail::HealEdgePair& pair : plan.pairs) {
    const Edge& keep_edge = brep.edge(pair.keep);
    if (keep_edge.curve == kNoGeometry ||
        keep_edge.curve >= curves.size()) {
      continue;
    }
    const CurveGeometry& curve = curves[keep_edge.curve];
    for (const VertexId endpoint : {keep_edge.v0, keep_edge.v1}) {
      const Point3d mapped = brep.vertices()[root_of(endpoint.index())].point;
      if (!(detail::heal_curve_residual(curve, mapped) <= check_tolerance)) {
        return Error{Errc::unsupported,
                     "heal_brep: kept edge curve does not match its welded "
                     "vertices",
                     pair.keep.index(), pair.drop.index()};
      }
    }
  }

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
    const Edge keep_before = editor.brep().edge(pair.keep);
    const Edge drop_before = editor.brep().edge(pair.drop);
    const bool inherits_curve = keep_before.curve == kNoGeometry &&
                                drop_before.curve != kNoGeometry;
    // heal_apply_welds already rewired the dropped edge endpoints, so the
    // dropped edge's v0 now names the surviving keep vertex at its traversal
    // start. When that is the keep edge's v1 the welded correspondence is
    // flipped: merge_edges still copies the dropped carrier, which then runs
    // against the keep edge's own direction.
    const bool flipped = keep_before.v0 != drop_before.v0;
    const Result<void> merged = editor.merge_edges(pair.keep, pair.drop);
    if (!merged) return reject(merged.error().message, true);
    ++sewn_edges;
    if (inherits_curve && flipped) {
      const Result<std::uint32_t> reversed =
          detail::heal_reverse_curve(editor, pair.keep);
      if (!reversed) return reject(reversed.error().message, true);
      const Edge& inherited = editor.brep().edge(pair.keep);
      const Result<void> restored = editor.set_edge_geometry(
          pair.keep, reversed.value(), inherited.t0, inherited.t1);
      if (!restored) return reject(restored.error().message, true);
    }
  }

  std::vector<char> welded_vertices(editor.brep().vertices().size(), 0);
  for (const detail::HealVertexWeld& weld : applied_welds) {
    if (editor.brep().vertex_alive(weld.keep) &&
        weld.keep.index() < welded_vertices.size()) {
      welded_vertices[weld.keep.index()] = 1;
    }
  }
  // merge_edges can hand a curve-less keep edge the dropped edge's carrier. That
  // copy follows the dropped edge's parameterization and direction, which the
  // pre-merge guard cannot see. Re-check every curve that now touches a welded
  // vertex: it must still pass through the surviving endpoints. The Brep here is
  // already mutated, so a failure returns an error instead of silently keeping a
  // curve that does not match its vertices.
  for (const EdgeId edge_id : editor.brep().edge_ids()) {
    const Edge& edge = editor.brep().edge(edge_id);
    if (edge.curve == kNoGeometry ||
        edge.curve >= editor.brep().geometry().curves.size()) {
      continue;
    }
    if (!welded_vertices[edge.v0.index()] &&
        !welded_vertices[edge.v1.index()]) {
      continue;
    }
    const CurveGeometry& curve = editor.brep().geometry().curves[edge.curve];
    for (const VertexId endpoint : {edge.v0, edge.v1}) {
      const Point3d mapped =
          editor.brep().vertices()[root_of(endpoint.index())].point;
      if (!(detail::heal_curve_residual(curve, mapped) <= check_tolerance)) {
        return Error{Errc::unsupported,
                     "heal_brep: kept edge curve does not match its welded "
                     "vertices",
                     edge_id.index(), 0u};
      }
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
  std::vector<char> dirty_faces(editor.brep().faces().size(), 0);
  for (const CoedgeId coedge_id : editor.brep().coedge_ids()) {
    const Coedge& coedge = editor.brep().coedge(coedge_id);
    if (coedge.edge.index() >= incident.size() ||
        incident[coedge.edge.index()] == 0) {
      continue;
    }
    if (coedge.face.index() < dirty_faces.size()) {
      dirty_faces[coedge.face.index()] = 1;
    }
  }
  for (const FaceId face_id : editor.brep().face_ids()) {
    if (face_id.index() >= dirty_faces.size() ||
        dirty_faces[face_id.index()] == 0) {
      continue;
    }
    const Result<void> rebuilt = detail::heal_rebuild_face_pcurves(
        editor, face_id, options.tolerance);
    if (!rebuilt) return reject(rebuilt.error().message, true);
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
