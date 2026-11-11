#pragma once

#include <cassert>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <span>
#include <utility>
#include <variant>
#include <vector>

#include <geom/brep/pcurve2d.hpp>
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
    staged_.geometry_.pcurves.push_back(std::move(pcurve));
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

  // MEV: splits `edge` at the interior curve parameter `t` by inserting a new
  // vertex there. The curve store is never subdivided, so the two new edges
  // keep the original curve index and only differ in their parameter intervals
  // ([t0, t] and [t, t1]); both inherit the edge tolerance, and the new vertex
  // inherits it too. Every coedge that uses the edge is replaced in its loop by
  // two coedges in traversal order with the same orientation. Every pcurve that
  // references the edge is split at the fraction of its chord length that
  // corresponds to the curve-parameter fraction (pcurves run along the edge from
  // t0 to t1, not along the coedge's traversal direction); both halves share the
  // computed junction point exactly, so join_vertex() can rebuild the polyline.
  // A pcurve that cannot be split aborts the call with its error and nothing is
  // staged. The old edge and its old coedges are released; the new vertex is
  // returned. Nothing is staged unless every check passes.
  Result<VertexId> split_edge(EdgeId edge, double t) {
    if (!staged_.edge_alive(edge)) {
      return Error{Errc::index_out_of_range, "solid_editor: split edge is not alive",
                   edge.index(), 0u};
    }
    const Edge split = staged_.edge(edge);
    if (!std::isfinite(t) || !(split.t0 < t) || !(t < split.t1)) {
      return Error{Errc::degenerate,
                   "solid_editor: split parameter is not inside the edge interval", 0u, 0u};
    }
    if (split.curve == kNoGeometry) {
      return Error{Errc::degenerate, "solid_editor: split edge has no curve", 0u, 0u};
    }
    if (static_cast<std::size_t>(split.curve) >= staged_.geometry_.curves.size()) {
      return Error{Errc::index_out_of_range, "solid_editor: split edge curve index out of range",
                   split.curve, 0u};
    }
    const Point3d split_point = evaluate_curve(staged_.geometry_.curves[split.curve], t);
    if (!point_is_finite(split_point)) {
      return Error{Errc::non_finite, "solid_editor: split point is not finite", 0u, 0u};
    }

    std::vector<CoedgeId> uses;
    std::vector<Coedge> old_coedges;
    std::uint32_t pcurves_needed = 0;
    for (std::size_t i = 0; i < staged_.coedges_.size(); ++i) {
      const Handle id = staged_.coedge_alloc_.handle_at(static_cast<std::uint32_t>(i));
      if (id.is_null() || staged_.coedges_[i].edge != edge) continue;
      uses.push_back(id);
      old_coedges.push_back(staged_.coedges_[i]);
      if (staged_.coedges_[i].pcurve != kNoGeometry) ++pcurves_needed;
    }
    const std::uint32_t use_count = static_cast<std::uint32_t>(uses.size());

    // Split every referenced pcurve before any handle is allocated: appending to
    // the store can reallocate it, so the halves are computed from copies first,
    // and a malformed pcurve must leave the staged topology unchanged.
    std::vector<Pcurve2d> low_parts(use_count);
    std::vector<Pcurve2d> high_parts(use_count);
    const double fraction = (t - split.t0) / (split.t1 - split.t0);
    for (std::uint32_t i = 0; i < use_count; ++i) {
      if (old_coedges[i].pcurve == kNoGeometry) continue;
      const Pcurve2d source = staged_.geometry_.pcurves[old_coedges[i].pcurve];
      Result<std::pair<Pcurve2d, Pcurve2d>> halves =
          detail::pcurve_split_halves(source, fraction);
      if (!halves.has_value()) return halves.error();
      low_parts[i] = std::move(halves.value().first);
      high_parts[i] = std::move(halves.value().second);
    }

    // One vertex + two edges are added, and each use trades one coedge for two,
    // so the coedge pool grows by use_count. pcurves grow by two per use that
    // references one. Check everything before any handle is allocated.
    const std::uint32_t coedge_live = staged_.coedge_alloc_.live_count();
    const std::size_t pcurve_live = staged_.geometry_.pcurves.size();
    if (!can_alloc(staged_.vertex_alloc_, 1u) || !can_alloc(staged_.edge_alloc_, 2u) ||
        use_count > Handle::kNullIndex - coedge_live ||
        pcurve_live > static_cast<std::size_t>(kNoGeometry) -
                          2u * static_cast<std::size_t>(pcurves_needed)) {
      return capacity_error();
    }

    const Handle vertex_id = staged_.vertex_alloc_.alloc();
    const Handle first_edge = staged_.edge_alloc_.alloc();
    const Handle second_edge = staged_.edge_alloc_.alloc();
    if (vertex_id.is_null() || first_edge.is_null() || second_edge.is_null()) {
      if (!first_edge.is_null()) staged_.edge_alloc_.release(first_edge);
      if (!second_edge.is_null()) staged_.edge_alloc_.release(second_edge);
      if (!vertex_id.is_null()) staged_.vertex_alloc_.release(vertex_id);
      return capacity_error();
    }
    Vertex inserted;
    inserted.point = split_point;
    inserted.tol = split.tol;
    place(staged_.vertices_, vertex_id, inserted);

    Edge first_half;
    first_half.curve = split.curve;
    first_half.t0 = split.t0;
    first_half.t1 = t;
    first_half.v0 = split.v0;
    first_half.v1 = vertex_id;
    first_half.tol = split.tol;
    place(staged_.edges_, first_edge, first_half);

    Edge second_half = first_half;
    second_half.t0 = t;
    second_half.t1 = split.t1;
    second_half.v0 = vertex_id;
    second_half.v1 = split.v1;
    place(staged_.edges_, second_edge, second_half);

    // Append the pre-split halves before releasing the old coedges. A use whose
    // pcurve was absent has an empty part and keeps kNoGeometry.
    std::vector<std::uint32_t> low_pcurves(use_count, kNoGeometry);
    std::vector<std::uint32_t> high_pcurves(use_count, kNoGeometry);
    for (std::uint32_t i = 0; i < use_count; ++i) {
      if (low_parts[i].points.empty()) continue;
      low_pcurves[i] = add_pcurve(std::move(low_parts[i]));
      high_pcurves[i] = add_pcurve(std::move(high_parts[i]));
    }

    for (const CoedgeId id : uses) staged_.coedge_alloc_.release(id);

    // The allocator bumps a generation on release, so a reused slot can never
    // alias an old handle; pair_of() below can safely match old neighbours.
    std::vector<std::pair<Handle, Handle>> pairs;
    pairs.reserve(use_count);
    for (std::uint32_t i = 0; i < use_count; ++i) {
      const Coedge& old = old_coedges[i];
      const Handle head_id = staged_.coedge_alloc_.alloc();
      const Handle tail_id = staged_.coedge_alloc_.alloc();
      // The capacity pre-check above reserved two coedges per use, so alloc()
      // cannot fail here; asserting documents that and keeps the loop atomic.
      assert(!head_id.is_null() && !tail_id.is_null());
      Coedge head;
      head.face = old.face;
      head.loop = old.loop;
      head.reversed = old.reversed;
      head.tol = old.tol;
      head.next = tail_id;
      if (old.reversed) {
        head.edge = second_edge;
        head.pcurve = high_pcurves[i];
      } else {
        head.edge = first_edge;
        head.pcurve = low_pcurves[i];
      }
      Coedge tail = head;
      tail.next = Handle{};
      tail.prev = head_id;
      if (old.reversed) {
        tail.edge = first_edge;
        tail.pcurve = low_pcurves[i];
      } else {
        tail.edge = second_edge;
        tail.pcurve = high_pcurves[i];
      }
      place(staged_.coedges_, head_id, head);
      place(staged_.coedges_, tail_id, tail);
      pairs.emplace_back(head_id, tail_id);
    }

    const auto pair_of = [&uses, &pairs](Handle id) -> const std::pair<Handle, Handle>* {
      for (std::size_t i = 0; i < uses.size(); ++i) {
        if (uses[i] == id) return &pairs[i];
      }
      return nullptr;
    };
    for (std::uint32_t i = 0; i < use_count; ++i) {
      const Coedge& old = old_coedges[i];
      const Handle head_id = pairs[i].first;
      const Handle tail_id = pairs[i].second;
      const std::pair<Handle, Handle>* pred_pair = pair_of(old.prev);
      if (pred_pair != nullptr) {
        staged_.coedges_[head_id.index()].prev = pred_pair->second;
      } else {
        staged_.coedges_[head_id.index()].prev = old.prev;
        staged_.coedges_[old.prev.index()].next = head_id;
      }
      const std::pair<Handle, Handle>* succ_pair = pair_of(old.next);
      if (succ_pair != nullptr) {
        staged_.coedges_[tail_id.index()].next = succ_pair->first;
      } else {
        staged_.coedges_[tail_id.index()].next = old.next;
        staged_.coedges_[old.next.index()].prev = tail_id;
      }
      Loop& loop = staged_.loops_[old.loop.index()];
      if (loop.start == uses[i]) loop.start = head_id;
    }
    staged_.edge_alloc_.release(edge);
    return vertex_id;
  }

  // KEV: joins the two edges meeting at `vertex` into one edge. The vertex must
  // be alive, be an endpoint of exactly two edges, and be used by exactly two
  // coedges (one per edge) that run through it head-to-tail inside one loop.
  // The edges must share one curve, their parameter intervals must meet exactly
  // at the vertex, and the curve tangents there must point the same way (dot >
  // 0). The merged edge keeps the lower interval's start and curve, extends to
  // the upper interval's end, and takes the larger of the two tolerances; the
  // two coedges collapse into one, and the vertex and the other edge are
  // released. The survivor's pcurve is rebuilt over the merged interval by
  // concatenating the lower half's points with the upper half's points, with the
  // shared junction stored once (pcurves run along the edge from t0 to t1, as in
  // split_edge); if either half has no pcurve the merged edge has none and the
  // two segments' continuity is not re-validated. Nothing is staged unless every
  // check passes.
  Result<bool> join_vertex(VertexId vertex) {
    if (!staged_.vertex_alive(vertex)) {
      return Error{Errc::index_out_of_range, "solid_editor: join vertex is not alive",
                   vertex.index(), 0u};
    }

    std::vector<EdgeId> incident;
    for (std::size_t i = 0; i < staged_.edges_.size(); ++i) {
      const Handle id = staged_.edge_alloc_.handle_at(static_cast<std::uint32_t>(i));
      if (id.is_null()) continue;
      const Edge& edge = staged_.edges_[i];
      if (edge.v0 == vertex || edge.v1 == vertex) incident.push_back(id);
    }
    if (incident.size() != 2u) {
      return Error{Errc::invalid_topology,
                   "solid_editor: join vertex does not have exactly two incident edges",
                   vertex.index(), 0u};
    }

    std::vector<CoedgeId> uses;
    for (std::size_t i = 0; i < staged_.coedges_.size(); ++i) {
      const Handle id = staged_.coedge_alloc_.handle_at(static_cast<std::uint32_t>(i));
      if (id.is_null()) continue;
      const Handle used_edge = staged_.coedges_[i].edge;
      if (used_edge == incident[0] || used_edge == incident[1]) uses.push_back(id);
    }
    if (uses.size() != 2u || staged_.coedge(uses[0]).edge == staged_.coedge(uses[1]).edge) {
      return Error{Errc::invalid_topology,
                   "solid_editor: join vertex does not have two incident coedges",
                   vertex.index(), 0u};
    }

    const Coedge coedge_a = staged_.coedge(uses[0]);
    const Coedge coedge_b = staged_.coedge(uses[1]);
    const Edge edge_a = staged_.edge(coedge_a.edge);
    const Edge edge_b = staged_.edge(coedge_b.edge);
    const auto ends_at = [vertex](const Coedge& coedge, const Edge& edge) {
      return (coedge.reversed ? edge.v0 : edge.v1) == vertex;
    };
    if (ends_at(coedge_a, edge_a) == ends_at(coedge_b, edge_b)) {
      return Error{Errc::invalid_topology,
                   "solid_editor: join coedges do not run in opposite directions at the vertex",
                   vertex.index(), 0u};
    }
    const CoedgeId end_id = ends_at(coedge_a, edge_a) ? uses[0] : uses[1];
    const CoedgeId start_id = ends_at(coedge_a, edge_a) ? uses[1] : uses[0];
    const Coedge end_coedge = staged_.coedge(end_id);
    const Coedge start_coedge = staged_.coedge(start_id);
    const Edge end_edge = staged_.edge(end_coedge.edge);
    const Edge start_edge = staged_.edge(start_coedge.edge);
    if (end_coedge.loop != start_coedge.loop || end_coedge.next != start_id) {
      return Error{Errc::invalid_topology,
                   "solid_editor: join coedges are not adjacent in one loop", vertex.index(),
                   0u};
    }
    if (end_edge.curve == kNoGeometry || start_edge.curve != end_edge.curve) {
      return Error{Errc::degenerate, "solid_editor: join edges do not share one curve",
                   start_edge.curve, end_edge.curve};
    }

    // Parameter order: the lower interval must end exactly where the upper
    // interval starts. Orientations must agree, otherwise the loop doubles back
    // at the vertex and there is no single merged interval.
    Edge merged;
    if (!end_coedge.reversed && !start_coedge.reversed) {
      if (end_edge.t1 != start_edge.t0) {
        return Error{Errc::degenerate,
                     "solid_editor: join edge intervals do not meet at the vertex",
                     vertex.index(), 0u};
      }
      merged = end_edge;
      merged.t1 = start_edge.t1;
      merged.v1 = start_edge.v1;
    } else if (end_coedge.reversed && start_coedge.reversed) {
      if (start_edge.t1 != end_edge.t0) {
        return Error{Errc::degenerate,
                     "solid_editor: join edge intervals do not meet at the vertex",
                     vertex.index(), 0u};
      }
      merged = start_edge;
      merged.t1 = end_edge.t1;
      merged.v1 = end_edge.v1;
    } else {
      return Error{Errc::degenerate,
                   "solid_editor: join coedges have inconsistent orientation", vertex.index(),
                   0u};
    }
    merged.tol = end_edge.tol > start_edge.tol ? end_edge.tol : start_edge.tol;

    if (static_cast<std::size_t>(merged.curve) >= staged_.geometry_.curves.size()) {
      return Error{Errc::index_out_of_range, "solid_editor: join curve index out of range",
                   merged.curve, 0u};
    }
    const double end_parameter = end_coedge.reversed ? end_edge.t0 : end_edge.t1;
    const double start_parameter = start_coedge.reversed ? start_edge.t1 : start_edge.t0;
    const Vec3d end_tangent =
        curve_derivative(staged_.geometry_.curves[merged.curve], end_parameter);
    const Vec3d start_tangent =
        curve_derivative(staged_.geometry_.curves[merged.curve], start_parameter);
    if (!vec_is_finite(end_tangent) || !vec_is_finite(start_tangent) ||
        !(dot(end_tangent, start_tangent) > 0.0)) {
      return Error{Errc::degenerate, "solid_editor: join edges are not tangent at the vertex",
                   vertex.index(), 0u};
    }

    const bool end_is_lower = !end_coedge.reversed;
    const EdgeId kept_edge = end_is_lower ? end_coedge.edge : start_coedge.edge;
    const EdgeId dropped_edge = end_is_lower ? start_coedge.edge : end_coedge.edge;

    // Rebuild the survivor's pcurve over the merged interval. pcurves run along
    // the edge from t0 to t1 (the convention split_edge uses), so the merged
    // polyline is the lower half followed by the upper half with the exactly
    // shared junction stored once. The two source segments' continuity is not
    // re-validated. If either half has no pcurve the merged edge has none either.
    const std::uint32_t lower_pcurve = end_is_lower ? end_coedge.pcurve : start_coedge.pcurve;
    const std::uint32_t upper_pcurve = end_is_lower ? start_coedge.pcurve : end_coedge.pcurve;
    std::uint32_t merged_pcurve = kNoGeometry;
    if (lower_pcurve != kNoGeometry && upper_pcurve != kNoGeometry) {
      if (staged_.geometry_.pcurves.size() >= static_cast<std::size_t>(kNoGeometry)) {
        return capacity_error();
      }
      const Pcurve2d lower = staged_.geometry_.pcurves[lower_pcurve];
      const Pcurve2d upper = staged_.geometry_.pcurves[upper_pcurve];
      merged_pcurve = add_pcurve(pcurve_merge(lower, upper));
    }

    staged_.edges_[kept_edge.index()] = merged;

    Coedge survivor = end_coedge;
    survivor.edge = kept_edge;
    survivor.pcurve = merged_pcurve;
    survivor.next = start_coedge.next;
    if (survivor.prev == start_id) survivor.prev = end_id;
    staged_.coedges_[end_id.index()] = survivor;
    staged_.coedges_[start_coedge.next.index()].prev = end_id;
    Loop& loop = staged_.loops_[end_coedge.loop.index()];
    if (loop.start == start_id) loop.start = end_id;
    staged_.coedge_alloc_.release(start_id);
    staged_.edge_alloc_.release(dropped_edge);
    staged_.vertex_alloc_.release(vertex);
    return true;
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

  static bool point_is_finite(const Point3d& point) noexcept {
    return std::isfinite(point[0]) && std::isfinite(point[1]) && std::isfinite(point[2]);
  }

  static bool vec_is_finite(const Vec3d& vector) noexcept {
    return std::isfinite(vector[0]) && std::isfinite(vector[1]) && std::isfinite(vector[2]);
  }

  static Point3d evaluate_curve(const CurveGeometry& curve, double t) {
    return std::visit([t](const auto& geometry) { return geometry.eval(t); }, curve);
  }

  static Vec3d curve_derivative(const CurveGeometry& curve, double t) {
    return std::visit([t](const auto& geometry) { return geometry.derivative(t); }, curve);
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
