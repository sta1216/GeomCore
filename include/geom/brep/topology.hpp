#pragma once

#include <cassert>
#include <cstddef>
#include <cstdint>
#include <span>
#include <variant>
#include <vector>

#include <geom/core/handle.hpp>
#include <geom/curves/analytic_curve.hpp>
#include <geom/curves/nurbs_curve.hpp>
#include <geom/primitives/point.hpp>
#include <geom/surfaces/analytic.hpp>
#include <geom/surfaces/nurbs_surface.hpp>
#include <geom/surfaces/plane.hpp>

namespace geom {
inline namespace v1 {

// Flat, index-addressed B-rep storage.
//
// Every entity type owns one HandleAllocator and one std::vector<Entity>; the
// vector is indexed by handle.index(), so a slot and the handle that owns it
// carry the same index. A vector only grows when a fresh index is allocated: a
// released slot keeps its old value in place and is handed out again by a later
// allocation, which overwrites it. The *_alive queries therefore check both the
// array bounds and the allocator generation, and every accessor that takes a
// handle is a precondition violation (debug assert) for a dead handle.
//
// Geometry is not stored in the entities. CurveGeometry, SurfaceGeometry and
// Pcurve2d objects live in a shared GeometryStore and an entity carries the
// 32-bit store index (kNoGeometry when it has none), so an edge with a line
// curve and an edge with a NURBS curve have identical topology layout. Shell
// membership (faces of a shell, shells of a solid) goes through the global
// shell_faces_ / solid_shells_ pools; a shell/solid records the first entry and
// the count, which makes the member lists contiguous spans instead of
// per-entity containers.

inline constexpr std::uint32_t kNoGeometry = 0xFFFFFFFFu;

using VertexId = Handle;
using EdgeId = Handle;
using CoedgeId = Handle;
using LoopId = Handle;
using FaceId = Handle;
using ShellId = Handle;
using SolidId = Handle;

struct Vertex {
  Point3d point{};
  double tol = 0.0;
};

// Edge of the curve store; t0/t1 are the curve parameter interval, v0/v1 the
// bounding vertices.
struct Edge {
  std::uint32_t curve = kNoGeometry;
  double t0 = 0.0;
  double t1 = 1.0;
  VertexId v0{};
  VertexId v1{};
  double tol = 0.0;
};

// One oriented use of an edge on a face. next/prev form the loop cycle in the
// direction the face traverses it; reversed records whether that direction runs
// against the edge's own orientation.
struct Coedge {
  EdgeId edge{};
  FaceId face{};
  LoopId loop{};
  bool reversed = false;
  std::uint32_t pcurve = kNoGeometry;
  CoedgeId next{};
  CoedgeId prev{};
  double tol = 0.0;
};

// A face boundary. outer marks the single outer boundary of the face, next
// links the loops of the same face.
struct Loop {
  FaceId face{};
  CoedgeId start{};
  LoopId next{};
  bool outer = false;
};

struct Face {
  std::uint32_t surface = kNoGeometry;
  LoopId first_loop{};
  bool reversed = false;
  double tol = 0.0;
};

// A contiguous run of shell_faces_ entries: [first_face, first_face + face_count).
struct Shell {
  std::uint32_t first_face = 0;
  std::uint32_t face_count = 0;
  bool closed = false;
};

// A contiguous run of solid_shells_ entries: [first_shell, first_shell + shell_count).
struct SolidEntity {
  std::uint32_t first_shell = 0;
  std::uint32_t shell_count = 0;
  double tol = 0.0;
};

using CurveGeometry = std::variant<LineCurve3, CircleCurve3, NurbsCurve>;
using SurfaceGeometry = std::variant<Plane3, Cylinder3, Cone3, Sphere3, Torus3, NurbsSurface>;

// Straight pcurve segment in the surface parameter domain.
struct Pcurve2d {
  Point2d p0{};
  Point2d p1{};
};

struct GeometryStore {
  std::vector<CurveGeometry> curves;
  std::vector<SurfaceGeometry> surfaces;
  std::vector<Pcurve2d> pcurves;
};

// Read-only view of a B-rep. All topology is created through SolidEditor, which
// is the only friend that can allocate slots; a committed Brep is immutable and
// is copied by value (allocators, entity arrays, pools and geometry all travel
// together), so returning a copy never aliases the editor's staging state.
//
// The *_count queries return the number of live entities: released slots are
// excluded even though they remain in the flat arrays until reused. The
// vertices()/edges()/... arrays are indexed by handle index and can therefore be
// longer than the corresponding count; iterate them together with *_alive or
// follow handles from the ownership chains instead.
class Brep {
 public:
  Brep() = default;

  bool vertex_alive(VertexId id) const noexcept {
    return vertex_alloc_.alive(id) && id.index() < vertices_.size();
  }
  // Precondition: vertex_alive(id).
  const Vertex& vertex(VertexId id) const noexcept {
    assert(vertex_alive(id));
    return vertices_[static_cast<std::size_t>(id.index())];
  }
  std::uint32_t vertex_count() const noexcept { return vertex_alloc_.live_count(); }

  bool edge_alive(EdgeId id) const noexcept {
    return edge_alloc_.alive(id) && id.index() < edges_.size();
  }
  // Precondition: edge_alive(id).
  const Edge& edge(EdgeId id) const noexcept {
    assert(edge_alive(id));
    return edges_[static_cast<std::size_t>(id.index())];
  }
  std::uint32_t edge_count() const noexcept { return edge_alloc_.live_count(); }

  bool coedge_alive(CoedgeId id) const noexcept {
    return coedge_alloc_.alive(id) && id.index() < coedges_.size();
  }
  // Precondition: coedge_alive(id).
  const Coedge& coedge(CoedgeId id) const noexcept {
    assert(coedge_alive(id));
    return coedges_[static_cast<std::size_t>(id.index())];
  }
  std::uint32_t coedge_count() const noexcept { return coedge_alloc_.live_count(); }

  bool loop_alive(LoopId id) const noexcept {
    return loop_alloc_.alive(id) && id.index() < loops_.size();
  }
  // Precondition: loop_alive(id).
  const Loop& loop(LoopId id) const noexcept {
    assert(loop_alive(id));
    return loops_[static_cast<std::size_t>(id.index())];
  }
  std::uint32_t loop_count() const noexcept { return loop_alloc_.live_count(); }

  bool face_alive(FaceId id) const noexcept {
    return face_alloc_.alive(id) && id.index() < faces_.size();
  }
  // Precondition: face_alive(id).
  const Face& face(FaceId id) const noexcept {
    assert(face_alive(id));
    return faces_[static_cast<std::size_t>(id.index())];
  }
  std::uint32_t face_count() const noexcept { return face_alloc_.live_count(); }

  bool shell_alive(ShellId id) const noexcept {
    return shell_alloc_.alive(id) && id.index() < shells_.size();
  }
  // Precondition: shell_alive(id).
  const Shell& shell(ShellId id) const noexcept {
    assert(shell_alive(id));
    return shells_[static_cast<std::size_t>(id.index())];
  }
  std::uint32_t shell_count() const noexcept { return shell_alloc_.live_count(); }

  bool solid_alive(SolidId id) const noexcept {
    return solid_alloc_.alive(id) && id.index() < solids_.size();
  }
  // Precondition: solid_alive(id).
  const SolidEntity& solid(SolidId id) const noexcept {
    assert(solid_alive(id));
    return solids_[static_cast<std::size_t>(id.index())];
  }
  std::uint32_t solid_count() const noexcept { return solid_alloc_.live_count(); }

  const std::vector<Vertex>& vertices() const noexcept { return vertices_; }
  const std::vector<Edge>& edges() const noexcept { return edges_; }
  const std::vector<Coedge>& coedges() const noexcept { return coedges_; }
  const std::vector<Loop>& loops() const noexcept { return loops_; }
  const std::vector<Face>& faces() const noexcept { return faces_; }
  const std::vector<Shell>& shells() const noexcept { return shells_; }
  const std::vector<SolidEntity>& solids() const noexcept { return solids_; }

  // Live solid ids in handle-index (construction) order. Allocates, so the
  // noexcept marker terminates on allocation failure rather than throwing
  // across the boundary; callers that cannot tolerate that should not call it
  // under memory pressure.
  std::vector<SolidId> solid_ids() const noexcept {
    std::vector<SolidId> ids;
    ids.reserve(solid_alloc_.live_count());
    for (std::uint32_t i = 0; i < static_cast<std::uint32_t>(solids_.size()); ++i) {
      const Handle id = solid_alloc_.handle_at(i);
      if (!id.is_null()) ids.push_back(id);
    }
    return ids;
  }

  const GeometryStore& geometry() const noexcept { return geometry_; }

  // Faces of `id`, as recorded in the global shell_faces_ pool. Precondition:
  // shell_alive(id); a stale shell whose pool range no longer fits yields an
  // empty span instead of reading out of bounds.
  std::span<const FaceId> shell_faces(ShellId id) const noexcept {
    assert(shell_alive(id));
    if (!shell_alive(id)) return {};
    const Shell& s = shells_[static_cast<std::size_t>(id.index())];
    if (s.first_face > shell_faces_.size()) return {};
    if (s.face_count > shell_faces_.size() - s.first_face) return {};
    return std::span<const FaceId>{shell_faces_}.subspan(s.first_face, s.face_count);
  }

  // Shells of `id`, as recorded in the global solid_shells_ pool. Precondition:
  // solid_alive(id), with the same tolerance for a stale pool range as
  // shell_faces().
  std::span<const ShellId> solid_shells(SolidId id) const noexcept {
    assert(solid_alive(id));
    if (!solid_alive(id)) return {};
    const SolidEntity& s = solids_[static_cast<std::size_t>(id.index())];
    if (s.first_shell > solid_shells_.size()) return {};
    if (s.shell_count > solid_shells_.size() - s.first_shell) return {};
    return std::span<const ShellId>{solid_shells_}.subspan(s.first_shell, s.shell_count);
  }

 private:
  friend class SolidEditor;

  HandleAllocator vertex_alloc_;
  HandleAllocator edge_alloc_;
  HandleAllocator coedge_alloc_;
  HandleAllocator loop_alloc_;
  HandleAllocator face_alloc_;
  HandleAllocator shell_alloc_;
  HandleAllocator solid_alloc_;

  std::vector<Vertex> vertices_;
  std::vector<Edge> edges_;
  std::vector<Coedge> coedges_;
  std::vector<Loop> loops_;
  std::vector<Face> faces_;
  std::vector<Shell> shells_;
  std::vector<SolidEntity> solids_;

  // Contiguous member pools addressed by Shell::first_face/face_count and
  // SolidEntity::first_shell/shell_count.
  std::vector<FaceId> shell_faces_;
  std::vector<ShellId> solid_shells_;

  GeometryStore geometry_;
};

}  // namespace v1
}  // namespace geom
