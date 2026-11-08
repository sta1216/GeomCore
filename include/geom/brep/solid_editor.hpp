#pragma once

#include <cmath>
#include <cstddef>
#include <cstdint>
#include <span>
#include <utility>
#include <vector>

#include <geom/brep/topology.hpp>
#include <geom/brep/validate.hpp>
#include <geom/core/error.hpp>
#include <geom/core/handle.hpp>
#include <geom/core/result.hpp>

namespace geom {
inline namespace v1 {

// One coedge of a new loop: the edge the coedge uses, whether it runs against
// that edge's orientation, the optional pcurve store index and the coincidence
// tolerance. An omitted pcurve is kNoGeometry, an omitted tolerance is 0.
struct CoedgeSpec {
  EdgeId edge{};
  bool reversed = false;
  std::uint32_t pcurve = kNoGeometry;
  double tol = 0.0;
};

// Transactional editor for a Brep.
//
// Topology changes go through add_vertex/add_edge/add_face/add_loop/add_shell/
// add_solid, which mutate a private staging copy. Every method validates its
// arguments before touching the staging copy and reports failure through
// Result, so a failed call leaves the staged topology exactly as it was.
// rollback() discards the staging changes and restores the base handed to the
// constructor; commit() validates the staging copy and returns it as a fresh
// Brep on success, leaving the editor untouched on failure.
//
// SolidEditor is Brep's only friend, so it is the sole path that can allocate
// entity slots. Curve/surface/pcurve geometry is appended to the staging
// geometry store and addressed by the 32-bit index the add_* method returns
// (kNoGeometry means "no geometry"); passing an index that does not exist yet
// is rejected as Errc::index_out_of_range.
//
// Allocation failure (std::bad_alloc from vector growth or a member-pool resize)
// is not caught: it propagates out of these methods instead of being mapped to
// an Errc. This is the repo-wide container-allocation policy and the editor adds
// no recovery of its own; every other failure is reported as a Result error.
class SolidEditor {
 public:
  SolidEditor() = default;
  explicit SolidEditor(const Brep& base) : base_(base), staged_(base) {}

  // Staged view: what commit() would publish if it validated successfully.
  const Brep& brep() const noexcept { return staged_; }

  // Store-index append helpers. Geometry is addressed by position in the
  // staging store and can never be released, so the returned index is stable
  // until rollback()/destruction. A store cannot exceed kNoGeometry entries
  // (the index simply does not fit the 32-bit fields that reference it), which
  // is far beyond what memory can hold.
  std::uint32_t add_curve(CurveGeometry curve) {
    staged_.geometry_.curves.push_back(std::move(curve));
    return static_cast<std::uint32_t>(staged_.geometry_.curves.size() - 1u);
  }
  std::uint32_t add_surface(SurfaceGeometry surface) {
    staged_.geometry_.surfaces.push_back(std::move(surface));
    return static_cast<std::uint32_t>(staged_.geometry_.surfaces.size() - 1u);
  }
  std::uint32_t add_pcurve(Pcurve2d pcurve) {
    staged_.geometry_.pcurves.push_back(pcurve);
    return static_cast<std::uint32_t>(staged_.geometry_.pcurves.size() - 1u);
  }

  Result<VertexId> add_vertex(Point3d p, double tol) {
    if (!(std::isfinite(p[0]) && std::isfinite(p[1]) && std::isfinite(p[2]))) {
      return Error{Errc::non_finite, "solid_editor: vertex position is not finite", 0u, 0u};
    }
    if (tol < 0.0 || !std::isfinite(tol)) {
      return Error{Errc::tolerance_violation,
                   "solid_editor: vertex tolerance must be finite and non-negative", 0u, 0u};
    }
    const Handle id = staged_.vertex_alloc_.alloc();
    if (id.is_null()) return capacity_error();
    Vertex vertex;
    vertex.point = p;
    vertex.tol = tol;
    place(staged_.vertices_, id, vertex);
    return id;
  }

  Result<EdgeId> add_edge(VertexId v0, VertexId v1, std::uint32_t curve, double t0, double t1,
                          double tol) {
    if (!staged_.vertex_alive(v0) || !staged_.vertex_alive(v1)) {
      return Error{Errc::index_out_of_range, "solid_editor: edge vertex is not alive",
                   v0.index(), v1.index()};
    }
    if (tol < 0.0 || !std::isfinite(tol)) {
      return Error{Errc::tolerance_violation,
                   "solid_editor: edge tolerance must be finite and non-negative", 0u, 0u};
    }
    if (!(std::isfinite(t0) && std::isfinite(t1)) || !(t0 < t1)) {
      return Error{Errc::degenerate,
                   "solid_editor: edge parameter interval is empty or non-finite", 0u, 0u};
    }
    if (curve != kNoGeometry &&
        static_cast<std::size_t>(curve) >= staged_.geometry_.curves.size()) {
      return Error{Errc::index_out_of_range, "solid_editor: edge curve index out of range",
                   curve, 0u};
    }
    const Handle id = staged_.edge_alloc_.alloc();
    if (id.is_null()) return capacity_error();
    Edge edge;
    edge.curve = curve;
    edge.t0 = t0;
    edge.t1 = t1;
    edge.v0 = v0;
    edge.v1 = v1;
    edge.tol = tol;
    place(staged_.edges_, id, edge);
    return id;
  }

  Result<FaceId> add_face(std::uint32_t surface, bool reversed, double tol) {
    if (surface != kNoGeometry &&
        static_cast<std::size_t>(surface) >= staged_.geometry_.surfaces.size()) {
      return Error{Errc::index_out_of_range, "solid_editor: face surface index out of range",
                   surface, 0u};
    }
    if (tol < 0.0 || !std::isfinite(tol)) {
      return Error{Errc::tolerance_violation,
                   "solid_editor: face tolerance must be finite and non-negative", 0u, 0u};
    }
    const Handle id = staged_.face_alloc_.alloc();
    if (id.is_null()) return capacity_error();
    Face face;
    face.surface = surface;
    face.reversed = reversed;
    face.tol = tol;
    place(staged_.faces_, id, face);
    return id;
  }

  // Creates a loop from coedges.size() fresh Coedges in the given order:
  // coedge i -> coedge i + 1 -> ... -> coedge 0. Each edge may be used at most
  // twice in one loop, so an edge can carry the two opposite-oriented uses a
  // manifold boundary needs; a third use would make the cycle non-manifold.
  Result<LoopId> add_loop(FaceId face, std::span<const CoedgeSpec> coedges, bool outer) {
    if (!staged_.face_alive(face)) {
      return Error{Errc::index_out_of_range, "solid_editor: loop face is not alive",
                   face.index(), 0u};
    }
    if (coedges.empty()) {
      return Error{Errc::invalid_topology, "solid_editor: loop needs at least one coedge", 0u,
                   0u};
    }
    for (const CoedgeSpec& spec : coedges) {
      if (!staged_.edge_alive(spec.edge)) {
        return Error{Errc::index_out_of_range, "solid_editor: loop coedge edge is not alive",
                     spec.edge.index(), 0u};
      }
      if (spec.tol < 0.0 || !std::isfinite(spec.tol)) {
        return Error{Errc::tolerance_violation,
                     "solid_editor: coedge tolerance must be finite and non-negative",
                     spec.edge.index(), 0u};
      }
      if (spec.pcurve != kNoGeometry &&
          static_cast<std::size_t>(spec.pcurve) >= staged_.geometry_.pcurves.size()) {
        return Error{Errc::index_out_of_range,
                     "solid_editor: loop coedge pcurve index out of range", spec.pcurve, 0u};
      }
    }
    // Occurrence counting is quadratic in the loop size. That is deliberate:
    // loops are short, and this keeps the validation allocation-free so the
    // failure paths cannot throw.
    for (std::size_t i = 0; i < coedges.size(); ++i) {
      std::size_t occurrences = 0;
      for (std::size_t j = 0; j < coedges.size(); ++j) {
        if (coedges[j].edge == coedges[i].edge) ++occurrences;
      }
      if (occurrences > 2) {
        return Error{Errc::invalid_topology, "solid_editor: edge used more than twice in a loop",
                     coedges[i].edge.index(), 0u};
      }
    }
    if (coedges.size() > Handle::kNullIndex ||
        !can_alloc(staged_.coedge_alloc_, static_cast<std::uint32_t>(coedges.size())) ||
        !can_alloc(staged_.loop_alloc_, 1u)) {
      return capacity_error();
    }

    // Locate the tail of the face's loop chain before creating anything: if the
    // chain is corrupt (a link to a dead or foreign loop, or a cycle) the call
    // fails without leaving a half-built loop behind. A null tail means the
    // face has no loops yet.
    LoopId tail{};
    {
      const std::size_t limit = staged_.loop_alloc_.live_count();
      std::size_t steps = 0;
      LoopId current = staged_.faces_[static_cast<std::size_t>(face.index())].first_loop;
      while (!current.is_null()) {
        if (!staged_.loop_alive(current) ||
            staged_.loops_[static_cast<std::size_t>(current.index())].face != face) {
          return Error{Errc::invalid_topology, "solid_editor: face loop chain is corrupt",
                       current.index(), 0u};
        }
        tail = current;
        current = staged_.loops_[static_cast<std::size_t>(current.index())].next;
        if (++steps > limit) {
          return Error{Errc::invalid_topology, "solid_editor: face loop chain is corrupt",
                       current.index(), 0u};
        }
      }
    }

    const std::uint32_t count = static_cast<std::uint32_t>(coedges.size());
    std::vector<CoedgeId> ids;
    ids.reserve(count);
    for (std::uint32_t i = 0; i < count; ++i) ids.push_back(staged_.coedge_alloc_.alloc());

    const Handle loop_id = staged_.loop_alloc_.alloc();
    for (std::uint32_t i = 0; i < count; ++i) {
      Coedge coedge;
      coedge.edge = coedges[i].edge;
      coedge.face = face;
      coedge.loop = loop_id;
      coedge.reversed = coedges[i].reversed;
      coedge.pcurve = coedges[i].pcurve;
      coedge.next = ids[(i + 1u) % count];
      coedge.prev = ids[(i + count - 1u) % count];
      coedge.tol = coedges[i].tol;
      place(staged_.coedges_, ids[i], coedge);
    }
    Loop loop;
    loop.face = face;
    loop.start = ids.front();
    loop.outer = outer;
    place(staged_.loops_, loop_id, loop);

    // The chain was validated above, so appending here cannot fail.
    if (tail.is_null()) {
      staged_.faces_[static_cast<std::size_t>(face.index())].first_loop = loop_id;
    } else {
      staged_.loops_[static_cast<std::size_t>(tail.index())].next = loop_id;
    }
    return loop_id;
  }

  Result<ShellId> add_shell(std::span<const FaceId> faces, bool closed) {
    if (faces.empty()) {
      return Error{Errc::invalid_topology, "solid_editor: shell needs at least one face", 0u, 0u};
    }
    for (const FaceId face : faces) {
      if (!staged_.face_alive(face)) {
        return Error{Errc::index_out_of_range, "solid_editor: shell face is not alive",
                     face.index(), 0u};
      }
    }
    const std::size_t first = staged_.shell_faces_.size();
    if (first > kNoGeometry || faces.size() > static_cast<std::size_t>(kNoGeometry) - first ||
        !can_alloc(staged_.shell_alloc_, 1u)) {
      return capacity_error();
    }
    const Handle id = staged_.shell_alloc_.alloc();
    Shell shell;
    shell.first_face = static_cast<std::uint32_t>(first);
    shell.face_count = static_cast<std::uint32_t>(faces.size());
    shell.closed = closed;
    staged_.shell_faces_.insert(staged_.shell_faces_.end(), faces.begin(), faces.end());
    place(staged_.shells_, id, shell);
    return id;
  }

  Result<SolidId> add_solid(std::span<const ShellId> shells, double tol) {
    if (shells.empty()) {
      return Error{Errc::invalid_topology, "solid_editor: solid needs at least one shell", 0u, 0u};
    }
    for (const ShellId shell : shells) {
      if (!staged_.shell_alive(shell)) {
        return Error{Errc::index_out_of_range, "solid_editor: solid shell is not alive",
                     shell.index(), 0u};
      }
    }
    if (tol < 0.0 || !std::isfinite(tol)) {
      return Error{Errc::tolerance_violation,
                   "solid_editor: solid tolerance must be finite and non-negative", 0u, 0u};
    }
    const std::size_t first = staged_.solid_shells_.size();
    if (first > kNoGeometry || shells.size() > static_cast<std::size_t>(kNoGeometry) - first ||
        !can_alloc(staged_.solid_alloc_, 1u)) {
      return capacity_error();
    }
    const Handle id = staged_.solid_alloc_.alloc();
    SolidEntity solid;
    solid.first_shell = static_cast<std::uint32_t>(first);
    solid.shell_count = static_cast<std::uint32_t>(shells.size());
    solid.tol = tol;
    staged_.solid_shells_.insert(staged_.solid_shells_.end(), shells.begin(), shells.end());
    place(staged_.solids_, id, solid);
    return id;
  }

  // Releases the edge's handle without repairing references to it. This is the
  // deliberate low-level corruption path: the staged topology now contains
  // dangling handles, and commit() rejects it until the affected loops are
  // rebuilt. A dead or out-of-range id is rejected with index_out_of_range and
  // leaves the editor unchanged.
  Result<void> erase_edge(EdgeId id) {
    if (!staged_.edge_alive(id)) {
      return Error{Errc::index_out_of_range, "solid_editor: edge is not alive", id.index(), 0u};
    }
    staged_.edge_alloc_.release(id);
    return {};
  }

  // Counterpart of erase_edge() for vertices; edges that still bound the vertex
  // are left dangling and commit() reports them. Same id rules.
  Result<void> erase_vertex(VertexId id) {
    if (!staged_.vertex_alive(id)) {
      return Error{Errc::index_out_of_range, "solid_editor: vertex is not alive", id.index(),
                   0u};
    }
    staged_.vertex_alloc_.release(id);
    return {};
  }

  // Validates the staging copy and returns it as a fresh Brep. Topology and
  // tolerance validation are delegated to validate_topology() and
  // validate_tolerances() so commit() and a direct validator call cannot
  // disagree; the first issue's code, entity and message become the Result
  // error's code, context_a and message.
  Result<Brep> commit() const {
    Brep candidate = staged_;
    ValidationReport report = validate_topology(candidate);
    if (report.ok()) report = validate_tolerances(candidate);
    if (!report.ok()) {
      const ValidationIssue& issue = report.issues.front();
      return Error{issue.code, issue.message, issue.entity, 0u};
    }
    return candidate;
  }

  // Drops every staged change and restores the base passed to the constructor
  // (an empty Brep for a default-constructed editor).
  void rollback() { staged_ = base_; }

 private:
  static Error capacity_error() noexcept {
    return Error{Errc::unsupported, "solid_editor: entity capacity exhausted", 0u, 0u};
  }

  static bool can_alloc(const HandleAllocator& allocator, std::uint32_t count) noexcept {
    const std::uint32_t live = allocator.live_count();
    return count <= Handle::kNullIndex - live;
  }

  // Stores `entity` at the slot the allocator handed out. A fresh allocation
  // always lands at the end of the array; a slot reused after a release lands
  // inside it and is overwritten. Growing first keeps index() == slot index,
  // which is the invariant every Brep accessor relies on.
  template <class Entity>
  static void place(std::vector<Entity>& slots, Handle id, const Entity& entity) {
    const std::size_t index = static_cast<std::size_t>(id.index());
    if (index >= slots.size()) slots.resize(index + 1u);
    slots[index] = entity;
  }

  // Base state for rollback(), and the mutable staging copy. Keeping the base
  // by value costs one extra copy per editor but makes rollback trivial and
  // exact; the alternative (replaying the inverse edits) would be fragile.
  Brep base_;
  Brep staged_;
};

}  // namespace v1
}  // namespace geom
