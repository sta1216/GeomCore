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

#include <cstddef>
#include <cstdint>
#include <vector>

#include <geom/brep/topology.hpp>
#include <geom/core/result.hpp>

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

}  // namespace v1
}  // namespace geom
