#pragma once

// Boolean operations between analytic B-rep solids whose boundary pieces meet
// in trimmed faces.
//
// The driver is the curved counterpart of `boolean_solid` (which is planar
// only). Operands must be built from planes and cylinders: every other analytic
// surface kind is reported as `unsupported` before any work is done (cone, sphere
// and torus operands are deferred to a later milestone, NURBS patches have no
// analytic parameterization here).
//
// The driver pairs face bounding boxes, asks the surface/surface intersector for
// a cut curve per pair, clips each cut to the face it will split, splits every
// face of both operands along the cuts and keeps the fragments that the
// requested membership rule puts on the result boundary. Cut curves are chordal
// polylines; whenever a pair's intersection is an analytic circle (an
// axis-parallel cylinder meeting a plane) the chords are healed back into exact
// `CircleCurve3` arc spans, exactly as `split_solid_by_plane` does, so a trimmed
// cap keeps its exact area instead of an inscribed polygon. The healing helpers
// are shared with `split_solid.hpp` (they live in `detail` there and are
// included, not modified). A cut that crosses the periodic seam of a cylindrical
// face is placed so that its ends land on the seam edges (u = 0 and u = 2 pi),
// which is what lets a closed rim cut the full lateral face. Cuts are clipped to
// both faces of the pair, and where two cuts of the same face meet at a boundary
// edge the shared corner is recomputed as the exact circle/plane crossing (a
// uv-scaled junction tolerance decides when two ends meet, an endpoint that stays
// 1e-2 inside the face is still rejected). Healing only rewrites runs whose
// samples lie on the circle: a straight clipped chord between two on-circle
// vertices stays a chord, and a run may wrap the end of a ring.
//
// Classification is two-sided, mirroring `detail::solid_boolean_classify_face`
// in `boolean_solid.hpp`: a fragment is kept when the membership rule differs on
// its two sides, and its normal is flipped when the result interior lies on the
// side its own surface normal points to. Sampling membership once per fragment is
// enough here because the other operand's mesh must not contain the
// representative point (that case is reported as a numeric failure). The
// classification mesh comes from `tessellate` when that accepts the operand and
// from the uv triangulation otherwise. The uv path conforms the triangulation to
// every boundary sample (the triangulator drops collinear vertices, and a
// cylinder rim is a straight line in uv), so faces that share a boundary agree on
// their subdivision. In both cases the mesh is welded and checked to be
// watertight before use, and a mesh that is not closed is reported as a numeric
// failure instead of being classified against.
//
// A flipped fragment is emitted with the flipped normal baked into its geometry
// rather than through a face flag: a planar piece is re-expressed on the opposite
// plane, its rings are reversed and its parameters reprojected, so every output
// face is counter-clockwise about its surface normal with `reversed == false`,
// exactly like `boolean_solid`'s output. That keeps two faces sharing a cut
// edge walking it in opposite directions, and keeps `curved_solid_volume`'s
// planar area positively oriented.
//
// Deferred: coplanar and tangent face contact, cone, sphere and torus operands,
// NURBS operands, flipped fragments on curved surfaces (no reparameterization
// here), results that fall apart into several connected components, and
// results that touch themselves along an edge (for example the
// `symmetric_difference` of two boxes that share a face corner edge). All of
// them are reported as `unsupported` instead of producing a wrong solid.

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <numeric>
#include <span>
#include <string_view>
#include <unordered_map>
#include <utility>
#include <variant>
#include <vector>

#include <geom/brep/assemble_curved.hpp>
#include <geom/brep/boolean_solid.hpp>
#include <geom/brep/detail/ssi_patch.hpp>
#include <geom/brep/detail/ssi_solver.hpp>
#include <geom/brep/detail/ssi_surface.hpp>
#include <geom/brep/face_uv.hpp>
#include <geom/brep/mesh_inside.hpp>
#include <geom/brep/plane_frame.hpp>
#include <geom/brep/split_solid.hpp>
#include <geom/brep/ssi.hpp>
#include <geom/brep/tessellate.hpp>
#include <geom/brep/topology.hpp>
#include <geom/brep/trim_face.hpp>
#include <geom/brep/validate.hpp>
#include <geom/cgeom/polygon2.hpp>
#include <geom/cgeom/triangulate2.hpp>
#include <geom/core/error.hpp>
#include <geom/core/result.hpp>
#include <geom/curves/analytic_curve.hpp>
#include <geom/surfaces/analytic.hpp>
#include <geom/surfaces/nurbs_surface.hpp>
#include <geom/surfaces/plane.hpp>

namespace geom {
inline namespace v1 {

struct CurvedBooleanOptions {
  double tolerance = 1e-9;
};

struct CurvedBooleanReport {
  std::uint32_t face_pairs = 0;
  std::uint32_t cuts = 0;
  std::uint32_t fragments = 0;
  std::uint32_t kept = 0;
  std::uint32_t dropped = 0;
  std::uint32_t flipped = 0;
  double max_cut_residual = 0.0;
};

struct CurvedBooleanResult {
  Brep brep;
  CurvedBooleanReport report;
};

namespace detail {

inline constexpr double kCurvedBooleanPi = 3.14159265358979323846;
inline constexpr double kCurvedBooleanTwoPi = 2.0 * kCurvedBooleanPi;

inline bool curved_boolean_valid_tolerance(double tolerance) noexcept {
  return std::isfinite(tolerance) && tolerance > 0.0;
}

// Axis-aligned 3D box used to reject face pairs that cannot meet. It is built
// from a 5x5 sample of the face's parameter domain, so it is a conservative
// underestimate of the face only up to one grid cell of a curved surface.
struct CurvedBooleanBounds {
  std::array<double, 3> lo{};
  std::array<double, 3> hi{};
};

inline bool curved_boolean_bounds_overlap(const CurvedBooleanBounds& a,
                                          const CurvedBooleanBounds& b,
                                          double tolerance) noexcept {
  for (int axis = 0; axis < 3; ++axis) {
    const std::size_t index = static_cast<std::size_t>(axis);
    if (a.lo[index] > b.hi[index] + tolerance) return false;
    if (b.lo[index] > a.hi[index] + tolerance) return false;
  }
  return true;
}

// The surface kinds the driver can cut. Cone operands have an apex that
// produces tangent loops the splitter cannot represent yet; sphere and torus
// operands have uv domains that fold (periodic u with poles, both axes periodic)
// and their classification mesh is deferred to M3c; NURBS patches have no
// analytic parameterization here.
struct CurvedBooleanSurfaceGate {
  bool supported = false;
  const char* message = nullptr;
};

inline CurvedBooleanSurfaceGate curved_boolean_surface_gate(
    const SurfaceGeometry& surface) noexcept {
  if (std::holds_alternative<Plane3>(surface) || std::holds_alternative<Cylinder3>(surface)) {
    return CurvedBooleanSurfaceGate{true, nullptr};
  }
  if (std::holds_alternative<NurbsSurface>(surface)) {
    return CurvedBooleanSurfaceGate{
        false, "boolean_solid_curved: NURBS operands are not supported"};
  }
  if (std::holds_alternative<Sphere3>(surface)) {
    return CurvedBooleanSurfaceGate{
        false, "boolean_solid_curved: sphere operands are deferred to M3c"};
  }
  if (std::holds_alternative<Torus3>(surface)) {
    return CurvedBooleanSurfaceGate{
        false, "boolean_solid_curved: torus operands are deferred to M3c"};
  }
  return CurvedBooleanSurfaceGate{
      false,
      "boolean_solid_curved: cone operands are not supported (apex-touching loops are "
      "deferred)"};
}

// One face of an operand in its surface parameter domain, with the domain, the
// face's 3D bounds and the projection needed to lift uv points back to 3D.
struct CurvedBooleanFaceInfo {
  FaceId id{};
  const SurfaceGeometry* surface = nullptr;
  FaceUvProjection projection{};
  std::vector<FaceUvLoop> loops;
  std::array<double, 2> u{};
  std::array<double, 2> v{};
  CurvedBooleanBounds bounds{};
};

inline Result<CurvedBooleanFaceInfo> curved_boolean_face_info(const Brep& brep, FaceId id,
                                                              double tolerance) {
  Result<std::vector<FaceUvLoop>> loops = face_loops_uv(brep, id, tolerance);
  if (!loops.has_value()) return loops.error();
  Result<FaceUvProjection> projection =
      face_uv_projection(brep, brep.face(id), tolerance);
  if (!projection.has_value()) return projection.error();

  CurvedBooleanFaceInfo info;
  info.id = id;
  info.surface = projection.value().surface;
  info.projection = projection.value();
  info.loops = std::move(loops.value());
  if (info.surface == nullptr) {
    return Error{Errc::invalid_topology, "boolean_solid_curved: face has no surface geometry",
                 id.index(), 0u};
  }

  double u_lo = std::numeric_limits<double>::infinity();
  double u_hi = -std::numeric_limits<double>::infinity();
  double v_lo = std::numeric_limits<double>::infinity();
  double v_hi = -std::numeric_limits<double>::infinity();
  for (const FaceUvLoop& loop : info.loops) {
    for (const Point2d& point : loop.points) {
      u_lo = std::min(u_lo, point[0]);
      u_hi = std::max(u_hi, point[0]);
      v_lo = std::min(v_lo, point[1]);
      v_hi = std::max(v_hi, point[1]);
    }
  }
  if (!(u_hi >= u_lo) || !(v_hi >= v_lo)) {
    return Error{Errc::degenerate, "boolean_solid_curved: face has an empty parameter domain",
                 id.index(), 0u};
  }
  // A periodic axis must not be sampled across more than one period: the
  // projection folds the seam, so a wider patch would trace it twice.
  if (info.projection.periodic_u() && u_hi - u_lo > kCurvedBooleanTwoPi) {
    u_hi = u_lo + kCurvedBooleanTwoPi;
  }
  info.u = {u_lo, u_hi};
  info.v = {v_lo, v_hi};

  bool have_bounds = false;
  for (int i = 0; i < 5; ++i) {
    for (int j = 0; j < 5; ++j) {
      const double fu = u_lo + (u_hi - u_lo) * (static_cast<double>(i) / 4.0);
      const double fv = v_lo + (v_hi - v_lo) * (static_cast<double>(j) / 4.0);
      const Result<SurfacePoint> sample = surface_eval(*info.surface, fu, fv);
      if (!sample.has_value()) return sample.error();
      const Point3d& point = sample.value().p;
      if (!have_bounds) {
        info.bounds.lo = {point[0], point[1], point[2]};
        info.bounds.hi = {point[0], point[1], point[2]};
        have_bounds = true;
        continue;
      }
      for (int axis = 0; axis < 3; ++axis) {
        const std::size_t slot = static_cast<std::size_t>(axis);
        info.bounds.lo[slot] = std::min(info.bounds.lo[slot], point[axis]);
        info.bounds.hi[slot] = std::max(info.bounds.hi[slot], point[axis]);
      }
    }
  }
  return info;
}

// Chord height a sampled boundary polyline may deviate from its exact circle.
constexpr double kCurvedBooleanMeshSag = 1e-3;

inline std::uint64_t curved_boolean_edge_key(std::uint32_t a, std::uint32_t b) noexcept {
  const std::uint32_t lo = a < b ? a : b;
  const std::uint32_t hi = a < b ? b : a;
  return (static_cast<std::uint64_t>(hi) << 32) | static_cast<std::uint64_t>(lo);
}

// Number of chords a circular boundary segment needs to stay within
// kCurvedBooleanMeshSag of its arc. It depends only on the circle and the sweep,
// so two faces that share the arc choose the same count.
inline int curved_boolean_arc_chords(const CircleCurve3& circle, double t0, double t1) noexcept {
  const double sweep = std::fabs(t1 - t0);
  if (!std::isfinite(sweep) || !(sweep > 0.0)) return 1;
  double step = 0.5 * kCurvedBooleanPi;
  if (std::isfinite(circle.radius) && circle.radius > kCurvedBooleanMeshSag) {
    const double cosine = 1.0 - kCurvedBooleanMeshSag / circle.radius;
    if (cosine > -1.0 && cosine < 1.0) step = 2.0 * std::acos(cosine);
  }
  if (!std::isfinite(step) || !(step > 1e-6)) step = 1e-6;
  const double chords = std::ceil(sweep / step);
  if (!(chords >= 1.0)) return 1;
  if (chords > 4096.0) return 4096;
  return static_cast<int>(chords);
}

// One face loop as a uv ring with every circular segment sampled on its own
// circle: the samples are exact curve points, and a face reached through a
// reversed coedge walks the same parameter range backwards, so both uses of an
// arc produce the same point set.
inline Result<std::vector<Point2d>> curved_boolean_mesh_ring(const CurvedBooleanFaceInfo& info,
                                                             const FaceUvLoop& loop) {
  const std::size_t count = loop.points.size();
  if (count < 3) {
    return Error{Errc::invalid_topology,
                 "boolean_solid_curved: a face loop has fewer than three points", info.id.index(),
                 0u};
  }
  std::vector<Point2d> ring;
  ring.reserve(count * 4);
  double previous_u = 0.0;
  bool has_previous_u = false;
  for (std::size_t i = 0; i < count; ++i) {
    const std::optional<FaceUvArcRef>* arc = i < loop.arcs.size() ? &loop.arcs[i] : nullptr;
    if (arc == nullptr || !arc->has_value()) {
      const Point2d uv = loop.points[i];
      ring.push_back(uv);
      previous_u = uv[0];
      has_previous_u = true;
      continue;
    }
    const FaceUvArcRef& span = arc->value();
    const int chords = curved_boolean_arc_chords(span.circle, span.t0, span.t1);
    for (int k = 0; k < chords; ++k) {
      const double t =
          span.t0 + (span.t1 - span.t0) * (static_cast<double>(k) / static_cast<double>(chords));
      const Result<Point2d> projected =
          face_uv_project(info.projection, span.circle.eval(t), previous_u, has_previous_u);
      if (!projected.has_value()) return projected.error();
      ring.push_back(projected.value());
      previous_u = projected.value()[0];
      has_previous_u = true;
    }
  }
  return ring;
}

// triangulate_polygon2 drops exactly collinear ring vertices, and a curved
// face's boundary samples are collinear in uv (a cylinder rim is a straight line
// in uv). A face that shares the boundary keeps every sample, so this mesh must
// keep them too: every triangle edge that spans a run of dropped ring vertices is
// split into a fan over the run. The fan only subdivides the triangle - the area
// and the winding are unchanged - and it makes each boundary edge identical to
// the neighbour's.
inline void curved_boolean_conform_boundary(
    const std::vector<Point2d>& ring_points, const std::vector<std::uint32_t>& ring_starts,
    std::vector<std::array<std::uint32_t, 3>>& triangles) {
  if (ring_starts.size() < 2u || ring_points.empty()) return;
  const std::size_t ring_count = ring_starts.size() - 1u;
  std::vector<bool> used(ring_points.size(), false);
  bool dropped = false;
  for (const std::array<std::uint32_t, 3>& triangle : triangles) {
    for (const std::uint32_t index : triangle) {
      if (index < used.size() && !used[index]) used[index] = true;
    }
  }
  for (const bool value : used) {
    if (!value) {
      dropped = true;
      break;
    }
  }
  if (!dropped) return;
  const auto ring_of = [&](std::uint32_t index) -> std::size_t {
    for (std::size_t r = 0; r < ring_count; ++r) {
      if (index >= ring_starts[r] && index < ring_starts[r + 1u]) return r;
    }
    return ring_count;
  };
  const auto same_point = [](const Point2d& a, const Point2d& b) noexcept {
    return a[0] == b[0] && a[1] == b[1];
  };
  std::vector<std::array<std::uint32_t, 3>> pending = triangles;
  std::vector<std::array<std::uint32_t, 3>> result;
  result.reserve(triangles.size());
  while (!pending.empty()) {
    const std::array<std::uint32_t, 3> triangle = pending.back();
    pending.pop_back();
    bool split = false;
    for (int k = 0; k < 3 && !split; ++k) {
      const std::uint32_t from = triangle[static_cast<std::size_t>(k)];
      const std::uint32_t to = triangle[static_cast<std::size_t>((k + 1) % 3)];
      const std::uint32_t apex = triangle[static_cast<std::size_t>((k + 2) % 3)];
      const std::size_t ring = ring_of(from);
      if (ring == ring_count || ring_of(to) != ring) continue;
      const std::uint32_t first = ring_starts[ring];
      const std::size_t size = ring_starts[ring + 1u] - first;
      const auto chain = [&](std::uint32_t start, std::uint32_t end) {
        std::vector<std::uint32_t> points;
        std::size_t i = start - first;
        const std::size_t last = end - first;
        while (true) {
          i = (i + 1u) % size;
          if (i == last) break;
          const std::uint32_t index = first + static_cast<std::uint32_t>(i);
          if (used[index]) {
            points.clear();
            return points;
          }
          points.push_back(index);
        }
        return points;
      };
      std::vector<std::uint32_t> run = chain(from, to);
      std::vector<std::uint32_t> reversed;
      if (run.empty()) {
        reversed = chain(to, from);
        run.swap(reversed);
        if (!run.empty()) std::reverse(run.begin(), run.end());
      }
      if (run.empty()) continue;
      std::vector<std::uint32_t> chain_indices;
      chain_indices.reserve(run.size() + 2u);
      chain_indices.push_back(from);
      for (const std::uint32_t index : run) {
        if (!same_point(ring_points[chain_indices.back()], ring_points[index])) {
          chain_indices.push_back(index);
        }
      }
      if (!same_point(ring_points[chain_indices.back()], ring_points[to])) {
        chain_indices.push_back(to);
      }
      if (chain_indices.size() < 3u) continue;
      for (std::size_t i = 0; i + 1u < chain_indices.size(); ++i) {
        pending.push_back({chain_indices[i], chain_indices[i + 1u], apex});
      }
      split = true;
    }
    if (!split) result.push_back(triangle);
  }
  triangles = std::move(result);
}

// Appends one face's triangles. Both planar and curved faces are triangulated
// from the projected boundary loops only: the mesh is welded and checked for
// closure afterwards, so it must not add interior vertices that a neighbouring
// face would not share.
inline Result<void> curved_boolean_mesh_face(const CurvedBooleanFaceInfo& info,
                                             std::vector<Point3d>& positions,
                                             std::vector<std::array<std::uint32_t, 3>>& indices) {
  std::vector<Point2d> ring_points;
  std::vector<std::uint32_t> ring_starts;
  for (const FaceUvLoop& loop : info.loops) {
    Result<std::vector<Point2d>> ring = curved_boolean_mesh_ring(info, loop);
    if (!ring.has_value()) return ring.error();
    if (ring.value().size() < 3) continue;
    ring_starts.push_back(static_cast<std::uint32_t>(ring_points.size()));
    ring_points.insert(ring_points.end(), ring.value().begin(), ring.value().end());
  }
  if (ring_starts.empty() || ring_starts.front() != 0u) {
    return Error{Errc::invalid_topology, "boolean_solid_curved: a face has no usable outer loop",
                 info.id.index(), 0u};
  }
  const std::uint32_t ring_count = static_cast<std::uint32_t>(ring_points.size());
  ring_starts.push_back(ring_count);

  Polygon2d polygon;
  polygon.vertices = ring_points;
  polygon.ring_starts = ring_starts;
  const Triangulation2 triangulation = triangulate_polygon2(polygon);
  if (!triangulation.ok || triangulation.triangles.empty()) {
    return Error{Errc::numeric_failure,
                 "boolean_solid_curved: a face outline could not be triangulated", info.id.index(),
                 0u};
  }

  const std::uint32_t base = static_cast<std::uint32_t>(positions.size());
  for (const Point2d& uv : ring_points) positions.push_back(face_uv_eval(info.projection, uv));
  std::vector<std::array<std::uint32_t, 3>> triangles = triangulation.triangles;
  curved_boolean_conform_boundary(ring_points, ring_starts, triangles);
  for (const std::array<std::uint32_t, 3>& triangle : triangles) {
    indices.push_back({base + triangle[0], base + triangle[1], base + triangle[2]});
  }
  return {};
}

// Edge-use census of a triangle mesh: a closed (watertight) mesh uses every
// undirected edge exactly twice.
struct CurvedBooleanMeshClosure {
  std::size_t edges = 0;
  std::size_t used_once = 0;
  std::size_t nonmanifold = 0;

  bool closed() const noexcept {
    return edges > 0 && used_once == 0 && nonmanifold == 0;
  }
};

inline CurvedBooleanMeshClosure curved_boolean_mesh_closure(const TriMesh& mesh) noexcept {
  CurvedBooleanMeshClosure closure;
  std::vector<std::uint64_t> keys;
  keys.reserve(mesh.triangles.size() * 3u);
  for (const std::array<std::uint32_t, 3>& triangle : mesh.triangles) {
    const std::uint32_t a = triangle[0];
    const std::uint32_t b = triangle[1];
    const std::uint32_t c = triangle[2];
    if (a >= mesh.positions.size() || b >= mesh.positions.size() || c >= mesh.positions.size()) {
      ++closure.nonmanifold;
      continue;
    }
    if (a == b || b == c || c == a) continue;
    keys.push_back(curved_boolean_edge_key(a, b));
    keys.push_back(curved_boolean_edge_key(b, c));
    keys.push_back(curved_boolean_edge_key(c, a));
  }
  std::sort(keys.begin(), keys.end());
  for (std::size_t i = 0; i < keys.size();) {
    std::size_t next = i + 1;
    while (next < keys.size() && keys[next] == keys[i]) ++next;
    const std::size_t uses = next - i;
    ++closure.edges;
    if (uses == 1) {
      ++closure.used_once;
    } else if (uses > 2) {
      ++closure.nonmanifold;
    }
    i = next;
  }
  return closure;
}

// Merges mesh vertices closer than `weld` into one and drops the triangles that
// collapse onto an edge. Positions keep their first-seen order so the result is
// deterministic.
inline void curved_boolean_mesh_weld(TriMesh& mesh, double weld) {
  if (mesh.positions.empty()) return;
  const double cell = weld > 0.0 ? weld : 1e-12;
  struct CellHash {
    std::size_t operator()(const std::array<long long, 3>& cell_value) const noexcept {
      std::uint64_t hash = 1469598103934665603ull;
      for (int axis = 0; axis < 3; ++axis) {
        const std::uint64_t bits =
            static_cast<std::uint64_t>(cell_value[static_cast<std::size_t>(axis)]);
        hash = (hash ^ bits) * 1099511628211ull;
      }
      return static_cast<std::size_t>(hash);
    }
  };
  std::unordered_map<std::array<long long, 3>, std::vector<std::uint32_t>, CellHash> cells;
  std::vector<std::uint32_t> remap(mesh.positions.size(), 0);
  std::vector<Point3d> kept;
  kept.reserve(mesh.positions.size());
  for (std::size_t i = 0; i < mesh.positions.size(); ++i) {
    const Point3d& point = mesh.positions[i];
    std::array<long long, 3> home{};
    for (int axis = 0; axis < 3; ++axis) {
      home[static_cast<std::size_t>(axis)] =
          static_cast<long long>(std::floor(point[axis] / cell));
    }
    std::uint32_t match = std::numeric_limits<std::uint32_t>::max();
    for (int dx = -1; dx <= 1 && match == std::numeric_limits<std::uint32_t>::max(); ++dx) {
      for (int dy = -1; dy <= 1 && match == std::numeric_limits<std::uint32_t>::max(); ++dy) {
        for (int dz = -1; dz <= 1 && match == std::numeric_limits<std::uint32_t>::max(); ++dz) {
          const std::array<long long, 3> key{home[0] + dx, home[1] + dy, home[2] + dz};
          const auto found = cells.find(key);
          if (found == cells.end()) continue;
          for (const std::uint32_t candidate : found->second) {
            const Point3d& other = kept[candidate];
            const double d0 = point[0] - other[0];
            const double d1 = point[1] - other[1];
            const double d2 = point[2] - other[2];
            if (d0 * d0 + d1 * d1 + d2 * d2 <= weld * weld) {
              match = candidate;
              break;
            }
          }
        }
      }
    }
    if (match != std::numeric_limits<std::uint32_t>::max()) {
      remap[i] = match;
      continue;
    }
    const std::uint32_t index = static_cast<std::uint32_t>(kept.size());
    kept.push_back(point);
    cells[home].push_back(index);
    remap[i] = index;
  }

  std::vector<std::array<std::uint32_t, 3>> triangles;
  triangles.reserve(mesh.triangles.size());
  for (const std::array<std::uint32_t, 3>& triangle : mesh.triangles) {
    if (triangle[0] >= remap.size() || triangle[1] >= remap.size() ||
        triangle[2] >= remap.size()) {
      continue;
    }
    const std::uint32_t a = remap[triangle[0]];
    const std::uint32_t b = remap[triangle[1]];
    const std::uint32_t c = remap[triangle[2]];
    if (a == b || b == c || c == a) continue;
    triangles.push_back({a, b, c});
  }

  mesh.positions = std::move(kept);
  mesh.triangles = std::move(triangles);
  mesh.triangle_face.resize(mesh.triangles.size(), 0u);
}

// Triangulates an operand's faces into one closed mesh for membership tests.
// Full primitives go through `tessellate`, which already welds its vertex pool
// and samples cylinder caps and bands consistently; trimmed operands (open rim
// arcs) are rejected there, so the driver samples the face parameter domains
// itself and welds the result. Either way the mesh must be watertight before it
// is used: classification against an open mesh is meaningless.
inline Result<TriMesh> curved_boolean_solid_mesh(const Brep& brep,
                                                 const std::vector<CurvedBooleanFaceInfo>& infos,
                                                 double tolerance) {
  double scale = 0.0;
  for (const CurvedBooleanFaceInfo& info : infos) {
    for (std::size_t axis = 0; axis < 3; ++axis) {
      scale = std::max(scale, info.bounds.hi[axis] - info.bounds.lo[axis]);
    }
  }
  if (!(scale > 0.0)) scale = 1.0;
  const double weld = std::max(64.0 * tolerance, 1e-9 * scale);

  const TessellationOptions tessellation_options;
  Result<TriMesh> tessellated = tessellate(brep, tessellation_options);
  if (tessellated.has_value() && !tessellated.value().triangles.empty()) {
    TriMesh mesh = std::move(tessellated.value());
    curved_boolean_mesh_weld(mesh, weld);
    if (!mesh.triangles.empty() && curved_boolean_mesh_closure(mesh).closed()) {
      return mesh;
    }
  }

  TriMesh mesh;
  for (const CurvedBooleanFaceInfo& info : infos) {
    const Result<void> status = curved_boolean_mesh_face(info, mesh.positions, mesh.triangles);
    if (!status.has_value()) return status.error();
    mesh.triangle_face.resize(mesh.triangles.size(), info.id.index());
  }
  if (mesh.triangles.empty()) {
    return Error{Errc::invalid_topology,
                 "boolean_solid_curved: an operand has no triangulable face", 0u, 0u};
  }
  curved_boolean_mesh_weld(mesh, weld);
  if (!curved_boolean_mesh_closure(mesh).closed()) {
    return Error{Errc::numeric_failure,
                 "boolean_solid_curved: classification mesh is not watertight", 0u, 0u};
  }
  return mesh;
}

// Wraps a periodic parameter into [base, base + 2 pi).
inline double curved_boolean_wrap_periodic(double value, double base) noexcept {
  const double offset = value - base;
  const double wrapped =
      offset - kCurvedBooleanTwoPi * std::floor(offset / kCurvedBooleanTwoPi);
  return base + wrapped;
}

// Unwraps one parameter axis of a cut chain so consecutive samples differ by
// less than half a period; a non-periodic axis is only clamped to the face.
inline void curved_boolean_unwrap_axis(std::vector<Point2d>& points, int axis, bool periodic,
                                       const std::array<double, 2>& range) {
  if (points.empty()) return;
  if (!periodic) {
    for (Point2d& point : points) {
      point[axis] = std::min(std::max(point[axis], range[0]), range[1]);
    }
    return;
  }
  for (std::size_t i = 1; i < points.size(); ++i) {
    const double delta = points[i][axis] - points[i - 1][axis];
    const double wrapped =
        delta - kCurvedBooleanTwoPi * std::floor((delta + kCurvedBooleanPi) / kCurvedBooleanTwoPi);
    points[i][axis] = points[i - 1][axis] + wrapped;
  }
}

// Splits a chain whose u values are wrapped into [base, base + 2 pi) where it
// steps across the periodic seam, inserting the seam point itself so both new
// ends land exactly on the seam edge they cross.
inline std::vector<std::vector<Point2d>> curved_boolean_split_at_seam(
    const std::vector<Point2d>& points, double base, double tolerance) {
  const double high = base + kCurvedBooleanTwoPi;
  std::vector<std::vector<Point2d>> chains;
  if (points.size() < 2) return chains;
  std::vector<Point2d> current{points.front()};
  const auto finish = [&chains, &current, tolerance]() {
    if (current.size() < 2) {
      current.clear();
      return;
    }
    std::vector<Point2d> chain;
    chain.reserve(current.size());
    for (const Point2d& point : current) {
      if (!chain.empty()) {
        const double du = point[0] - chain.back()[0];
        const double dv = point[1] - chain.back()[1];
        if (du * du + dv * dv <= tolerance * tolerance) continue;
      }
      chain.push_back(point);
    }
    if (chain.size() >= 2) chains.push_back(std::move(chain));
    current.clear();
  };
  for (std::size_t i = 1; i < points.size(); ++i) {
    const Point2d& previous = points[i - 1];
    const Point2d& point = points[i];
    const double delta = point[0] - previous[0];
    if (std::abs(delta) <= kCurvedBooleanPi) {
      current.push_back(point);
      continue;
    }
    // The step wraps across the seam: `rising` marks a jump from the high seam
    // edge down to the low one.
    const bool rising = delta < 0.0;
    const double low_part = previous[0] - base;
    const double high_part = high - point[0];
    const double span = rising ? (high - previous[0]) + (point[0] - base) : low_part + high_part;
    double t = 0.5;
    if (span > 0.0) t = rising ? (high - previous[0]) / span : low_part / span;
    const double crossing_u = rising ? high : base;
    const Point2d crossing{crossing_u, previous[1] + t * (point[1] - previous[1])};
    current.push_back(crossing);
    finish();
    current.push_back(Point2d{rising ? base : high, crossing[1]});
    current.push_back(point);
  }
  finish();
  return chains;
}

// Turns one intersector curve into the cut chains of one side. `segment.ua` /
// `ub` are the A-side parameters at the segment's two endpoints, `va` / `vb` the
// B-side ones. On a face with a periodic u a chain that winds once around the
// face is re-expressed as a single chain whose ends land exactly on the two seam
// edges, and a chain that steps across the seam is split there; both forms give
// the splitter ends that lie on the seam boundary instead of in the middle of
// the face.
inline std::vector<std::vector<Point2d>> curved_boolean_cut_chains(
    const SsiCurve& curve, bool a_side, const CurvedBooleanFaceInfo& face, double tolerance) {
  std::vector<Point2d> points;
  points.reserve(curve.segments.size() + 1);
  auto push = [&points](const std::array<double, 2>& params) {
    points.push_back(Point2d{params[0], params[1]});
  };
  for (const SsiSegment& segment : curve.segments) {
    push(a_side ? segment.ua : segment.va);
  }
  if (!curve.segments.empty()) {
    push(a_side ? curve.segments.back().ub : curve.segments.back().vb);
  }

  curved_boolean_unwrap_axis(points, 0, face.projection.periodic_u(), face.u);
  curved_boolean_unwrap_axis(points, 1, face.projection.periodic_v(), face.v);

  // Drop the points that only duplicate their predecessor: a segment shorter
  // than the tolerance carries no cutting information.
  std::vector<Point2d> unique;
  unique.reserve(points.size());
  for (const Point2d& point : points) {
    if (!unique.empty()) {
      const double du = point[0] - unique.back()[0];
      const double dv = point[1] - unique.back()[1];
      if (du * du + dv * dv <= tolerance * tolerance) continue;
    }
    unique.push_back(point);
  }
  if (unique.size() < 2) return {};

  if (!face.projection.periodic_u()) return {std::move(unique)};

  const double base = face.u[0];
  const double high = base + kCurvedBooleanTwoPi;
  if (curve.closed) {
    const double winding = unique.back()[0] - unique.front()[0];
    if (std::fabs(winding) > kCurvedBooleanPi) {
      const bool forward = winding > 0.0;
      bool monotone = true;
      for (std::size_t i = 1; i < unique.size(); ++i) {
        const double delta = unique[i][0] - unique[i - 1][0];
        if (forward ? delta < -tolerance : delta > tolerance) {
          monotone = false;
          break;
        }
      }
      const double span = std::fabs(winding);
      if (monotone && span <= kCurvedBooleanTwoPi + 64.0 * tolerance) {
        if (!forward) std::reverse(unique.begin(), unique.end());
        // The intersector answers in its own period ([-pi, pi] for a cylinder),
        // which is offset from the face's seam by an arbitrary number of turns.
        // Wrapping every sample into the face's period puts the single seam
        // crossing back where it belongs; the chain is then rotated to start on
        // the low seam edge and end on the high one, so both ends are boundary
        // points and the cut closes over the seam instead of through the middle
        // of the face.
        const double start = unique.front()[0];
        const double seam =
            base + std::ceil((start - base) / kCurvedBooleanTwoPi - 1e-12) * kCurvedBooleanTwoPi;
        double crossing_v = unique.front()[1];
        std::size_t after = unique.size();
        for (std::size_t i = 1; i < unique.size(); ++i) {
          const double u0 = unique[i - 1][0];
          const double u1 = unique[i][0];
          if (u1 < seam) continue;
          const double length = u1 - u0;
          const double t = length > 0.0 ? (seam - u0) / length : 0.0;
          crossing_v = unique[i - 1][1] + t * (unique[i][1] - unique[i - 1][1]);
          after = i;
          break;
        }
        if (after < unique.size()) {
          // Re-order the lap so it starts at the seam: the samples after the
          // crossing come first, the samples before it follow shifted by one
          // turn.
          std::vector<Point2d> chain;
          chain.reserve(unique.size() + 2);
          chain.push_back(Point2d{base, crossing_v});
          for (std::size_t step = 0; step < unique.size(); ++step) {
            const std::size_t i = (after + step) % unique.size();
            double u = unique[i][0] + (i < after ? kCurvedBooleanTwoPi : 0.0);
            u = base + (u - seam);
            if (u >= high - tolerance) continue;
            if (u <= base + tolerance) continue;
            const Point2d candidate{u, unique[i][1]};
            const double du = candidate[0] - chain.back()[0];
            const double dv = candidate[1] - chain.back()[1];
            if (du * du + dv * dv <= tolerance * tolerance) continue;
            chain.push_back(candidate);
          }
          chain.push_back(Point2d{high, crossing_v});
          return {std::move(chain)};
        }
      }
    }
  }
  for (Point2d& point : unique) point[0] = curved_boolean_wrap_periodic(point[0], base);
  if (curve.closed) {
    bool crosses = false;
    for (std::size_t i = 1; i < unique.size(); ++i) {
      if (std::abs(unique[i][0] - unique[i - 1][0]) > kCurvedBooleanPi) {
        crosses = true;
        break;
      }
    }
    if (!crosses) return {std::move(unique)};
  }
  return curved_boolean_split_at_seam(unique, base, tolerance);
}

// A closed cut that crosses a periodic seam becomes a single closed rim on the
// periodic face, running from the seam crossing at u = base to the same point
// at u = base + 2 pi. The other face carries the same rim without a seam, so
// its chain has to be re-based onto that crossing point: both faces then share
// one closed edge used twice, instead of two self-loops that never meet.
inline bool curved_boolean_rebase_closed_chain(std::vector<Point2d>& chain, const Point3d& start,
                                               const FaceUvProjection& projection,
                                               double tolerance) {
  if (chain.size() < 5) return false;
  const std::size_t count = chain.size() - 1;
  std::vector<Point3d> points;
  points.reserve(chain.size());
  for (const Point2d& uv : chain) {
    const Result<Point3d> point = face_uv_eval(projection, uv);
    if (!point.has_value()) return false;
    points.push_back(point.value());
  }
  double best = std::numeric_limits<double>::infinity();
  std::size_t segment = count;
  double fraction = 0.0;
  for (std::size_t i = 0; i < count; ++i) {
    const Point3d& begin = points[i];
    const Point3d& end = points[(i + 1) % count];
    const Vec3d span = end - begin;
    const double length2 = dot(span, span);
    double t = 0.0;
    if (length2 > 0.0) {
      t = dot(start - begin, span) / length2;
      t = std::min(1.0, std::max(0.0, t));
    }
    const Vec3d offset = (begin + span * t) - start;
    const double distance = std::sqrt(dot(offset, offset));
    if (distance < best) {
      best = distance;
      segment = i;
      fraction = t;
    }
  }
  if (segment >= count) return false;
  // The seam point sits on the rim arc between two samples, so it is only
  // within a chord's sag of the closed polyline: require it to be closer than
  // half the longest chord, which keeps a stray point from being inserted.
  double chord = 0.0;
  for (std::size_t i = 0; i < count; ++i) {
    const Vec3d span = points[(i + 1) % count] - points[i];
    chord = std::max(chord, std::sqrt(dot(span, span)));
  }
  if (!(best <= 0.5 * chord)) {
    return false;
  }
  const Result<Point2d> seam = face_uv_project(projection, start, chain[segment][0], true);
  if (!seam.has_value()) {
    return false;
  }
  const Result<Point3d> seam_point = face_uv_eval(projection, seam.value());
  if (!seam_point.has_value()) return false;
  const Vec3d seam_offset = seam_point.value() - start;
  const double seam_distance = std::sqrt(dot(seam_offset, seam_offset));
  const double start_scale = 1.0 + std::sqrt(start[0] * start[0] + start[1] * start[1] +
                                            start[2] * start[2]);
  if (!(seam_distance <= 1e-6 * start_scale)) return false;
  const std::size_t last_step = fraction <= 1e-12 ? count - 1 : count;
  std::vector<Point2d> rebased;
  rebased.reserve(last_step + 2);
  rebased.push_back(seam.value());
  for (std::size_t step = 1; step <= last_step; ++step) {
    Point2d point = chain[(segment + step) % count];
    if (projection.periodic_u()) point[0] = face_uv_unwrap(point[0], rebased.back()[0]);
    rebased.push_back(point);
  }
  rebased.push_back(seam.value());
  std::vector<Point2d> cleaned;
  cleaned.reserve(rebased.size());
  for (const Point2d& point : rebased) {
    if (!cleaned.empty()) {
      const double du = point[0] - cleaned.back()[0];
      const double dv = point[1] - cleaned.back()[1];
      if (du * du + dv * dv <= tolerance * tolerance) continue;
    }
    cleaned.push_back(point);
  }
  if (cleaned.size() < 4) return false;
  const double du = cleaned.back()[0] - cleaned.front()[0];
  const double dv = cleaned.back()[1] - cleaned.front()[1];
  if (du * du + dv * dv > tolerance * tolerance) cleaned.push_back(cleaned.front());
  chain = std::move(cleaned);
  return true;
}

// Finds the point where the 3D segment p0 -> p1 leaves the trimmed region of a
// face: a boundary arc contributes its exact circle crossing, a straight
// boundary edge a segment-to-segment contact. When the cut lies on a known
// circle, a straight boundary edge is intersected with that circle exactly: the
// chord between two samples lands inside the circle by the sagitta, which is far
// above the assembly weld and would leave the two faces' rims subdivided
// differently. The crossing nearest the middle of the segment wins, so a chain
// vertex that grazes the boundary is treated as a single crossing.
inline Result<Point3d> curved_boolean_boundary_crossing(const CurvedBooleanFaceInfo& face,
                                                        const Point3d& p0, const Point3d& p1,
                                                        const CircleCurve3* curve,
                                                        const Point3d* prefer,
                                                        double tolerance) {
  Point3d best;
  double best_score = std::numeric_limits<double>::infinity();
  bool found = false;
  const auto consider = [&](const Point3d& point) {
    const Vec3d span = p1 - p0;
    const double length2 = dot(span, span);
    double score = 0.5;
    if (prefer != nullptr) {
      const Vec3d gap = point - *prefer;
      const double distance = std::sqrt(dot(gap, gap));
      if (distance > 1e-4 * (1.0 + std::sqrt(length2))) return;
      score = distance;
    } else {
      if (length2 > 0.0) score = std::abs(dot(point - p0, span) / length2 - 0.5);
    }
    if (score < best_score) {
      best_score = score;
      best = point;
      found = true;
    }
  };
  for (const FaceUvLoop& loop : face.loops) {
    const std::size_t count = loop.points.size();
    for (std::size_t i = 0; i < count; ++i) {
      const std::size_t next = (i + 1) % count;
      if (i < loop.arcs.size() && loop.arcs[i].has_value()) {
        const FaceUvArcRef& ref = *loop.arcs[i];
        for (const TrimFaceCircleHit& hit : trim_face_circle_hits(ref.circle, p0, p1, tolerance)) {
          double offset = hit.t - ref.t0;
          while (offset > kCurvedBooleanPi) offset -= kCurvedBooleanTwoPi;
          while (offset < -kCurvedBooleanPi) offset += kCurvedBooleanTwoPi;
          const double lower = std::min(0.0, ref.t1 - ref.t0);
          const double upper = std::max(0.0, ref.t1 - ref.t0);
          if (offset < lower - tolerance || offset > upper + tolerance) continue;
          consider(hit.point);
        }
        continue;
      }
      const Result<Point3d> begin = face_uv_eval(face.projection, loop.points[i]);
      const Result<Point3d> end = face_uv_eval(face.projection, loop.points[next]);
      if (!begin.has_value() || !end.has_value()) continue;
      if (curve != nullptr) {
        const Vec3d normal = [&]() {
          const Vec3d raw = cross(curve->x_axis, curve->y_axis);
          const double length = std::sqrt(dot(raw, raw));
          if (!(length > 0.0)) return raw;
          return Vec3d{raw[0] / length, raw[1] / length, raw[2] / length};
        }();
        const auto in_plane = [&](const Point3d& point) {
          const double offset = dot(point - curve->center, normal);
          return Point3d{point[0] - normal[0] * offset, point[1] - normal[1] * offset,
                         point[2] - normal[2] * offset};
        };
        const Point3d plane_begin = in_plane(begin.value());
        const Point3d plane_end = in_plane(end.value());
        const Vec3d edge = plane_end - plane_begin;
        const Vec3d offset = plane_begin - curve->center;
        const double a = dot(edge, edge);
        if (a > 1e-30) {
          const double b = 2.0 * dot(offset, edge);
          const double c = dot(offset, offset) - curve->radius * curve->radius;
          const double discriminant = b * b - 4.0 * a * c;
          if (discriminant >= 0.0) {
            const double root = std::sqrt(discriminant);
            for (const double t : {(-b - root) / (2.0 * a), (-b + root) / (2.0 * a)}) {
              if (t < -1e-12 || t > 1.0 + 1e-12) continue;
              consider(Point3d{plane_begin[0] + edge[0] * t, plane_begin[1] + edge[1] * t,
                               plane_begin[2] + edge[2] * t});
            }
          }
        }
        continue;
      }
      const Vec3d cut = p1 - p0;
      const Vec3d edge = end.value() - begin.value();
      const Vec3d offset = p0 - begin.value();
      const double cut2 = dot(cut, cut);
      const double edge2 = dot(edge, edge);
      const double mixed = dot(cut, edge);
      const double denom = cut2 * edge2 - mixed * mixed;
      double cut_t = 0.0;
      double edge_t = 0.0;
      if (denom > 1e-24) {
        cut_t = (mixed * dot(edge, offset) - edge2 * dot(cut, offset)) / denom;
        edge_t = (cut2 * dot(edge, offset) - mixed * dot(cut, offset)) / denom;
      }
      cut_t = std::min(1.0, std::max(0.0, cut_t));
      edge_t = std::min(1.0, std::max(0.0, edge_t));
      const Point3d cut_point = p0 + cut * cut_t;
      const Point3d edge_point = begin.value() + edge * edge_t;
      const Vec3d gap = cut_point - edge_point;
      if (std::sqrt(dot(gap, gap)) <= tolerance) {
        consider(Point3d{(cut_point[0] + edge_point[0]) * 0.5, (cut_point[1] + edge_point[1]) * 0.5,
                         (cut_point[2] + edge_point[2]) * 0.5});
      }
    }
  }
  if (!found) {
    return Error{Errc::numeric_failure,
                 "boolean_solid_curved: a cut leaves an operand face without crossing its boundary",
                 0u, 0u};
  }
  return best;
}

// Keeps the part of a cut that lies inside both trimmed regions. The cut is the
// segment intersection of the two surfaces, so it is only a real edge of the
// result where it runs over both faces; the parameterizations are clipped at
// the exact boundary crossing so both faces meet in the same three-dimensional
// vertex.
inline std::vector<std::vector<Point2d>> curved_boolean_clip_pair(
    const CurvedBooleanFaceInfo& face, const std::vector<Point2d>& chain,
    const CurvedBooleanFaceInfo& other, const CircleCurve3* curve, double tolerance) {
  const std::size_t count = chain.size();
  if (count < 2) return {};
  std::vector<std::vector<Point2d>> loops;
  loops.reserve(face.loops.size());
  for (const FaceUvLoop& loop : face.loops) loops.push_back(loop.points);
  std::vector<std::vector<Point2d>> other_loops;
  other_loops.reserve(other.loops.size());
  for (const FaceUvLoop& loop : other.loops) other_loops.push_back(loop.points);

  std::vector<Point3d> points;
  points.reserve(count);
  for (const Point2d& uv : chain) {
    const Result<Point3d> point = face_uv_eval(face.projection, uv);
    if (!point.has_value()) return {chain};
    points.push_back(point.value());
  }
  std::vector<bool> inside(count, false);
  bool all_inside = true;
  double other_u_hint = chain.front()[0];
  for (std::size_t i = 0; i < count; ++i) {
    const bool own = trim_face_uv_inside(loops, chain[i], tolerance);
    bool foreign = false;
    const Result<Point2d> project = face_uv_project(other.projection, points[i], other_u_hint, true);
    if (project.has_value()) {
      other_u_hint = project.value()[0];
      Point2d foreign_uv = project.value();
      if (other.projection.periodic_u()) {
        foreign_uv[0] = curved_boolean_wrap_periodic(foreign_uv[0], other.u[0]);
      }
      foreign = trim_face_uv_inside(other_loops, foreign_uv, tolerance);
    }
    inside[i] = own && foreign;
    if (!inside[i]) all_inside = false;
  }
  std::vector<std::vector<Point2d>> pieces;
  if (all_inside) {
    pieces.push_back(chain);
  } else {
  std::vector<Point2d> clipped;
  const auto flush = [&]() {
    if (clipped.size() >= 2) pieces.push_back(std::move(clipped));
    clipped.clear();
  };
  const auto append = [&](const Point2d& point) {
    if (!clipped.empty()) {
      const double du = point[0] - clipped.back()[0];
      const double dv = point[1] - clipped.back()[1];
      if (du * du + dv * dv <= tolerance * tolerance) return;
    }
    clipped.push_back(point);
  };
  const auto crossing = [&](std::size_t from, std::size_t to) {
    const Result<Point3d> hit =
        curved_boolean_boundary_crossing(face, points[from], points[to], curve, nullptr, tolerance);
    if (hit.has_value()) {
      const Result<Point2d> uv = face_uv_project(face.projection, hit.value(), chain[from][0], true);
      if (uv.has_value()) return uv.value();
    }
    const Result<Point3d> other_hit =
        curved_boolean_boundary_crossing(other, points[from], points[to], curve, nullptr, tolerance);
    if (other_hit.has_value()) {
      const Result<Point2d> uv = face_uv_project(face.projection, other_hit.value(), chain[from][0], true);
      if (uv.has_value()) return uv.value();
    }
    return Point2d{(chain[from][0] + chain[to][0]) * 0.5, (chain[from][1] + chain[to][1]) * 0.5};
  };
  for (std::size_t i = 0; i + 1 < count; ++i) {
    const bool a = inside[i];
    const bool b = inside[i + 1];
    if (a && b) {
      append(chain[i]);
      append(chain[i + 1]);
      continue;
    }
    if (!a && !b) {
      flush();
      continue;
    }
    const Point2d edge = crossing(i, i + 1);
    if (a) {
      append(chain[i]);
      append(edge);
      flush();
    } else {
      flush();
      append(edge);
      append(chain[i + 1]);
    }
  }
  flush();
  }
  // The intersector traces within the surfaces' parameter rectangles, so an end
  // of the clipped chain can sit a few micro-units inside the true boundary
  // crossing. When the cut lies on a known circle, both faces snap their open
  // ends onto the exact circle-boundary crossing: otherwise one face heals its
  // rim into an arc while the other keeps a chord, and the assembler cannot pair
  // the two subdivisions.
  const double diagonal_u = face.u[1] - face.u[0];
  const double diagonal_v = face.v[1] - face.v[0];
  const double junction_tol =
      std::max(64.0 * tolerance, 1e-5 * std::sqrt(diagonal_u * diagonal_u + diagonal_v * diagonal_v));
  for (std::vector<Point2d>& piece : pieces) {
    if (piece.size() < 2) continue;
    for (const int end : {0, 1}) {
      const std::size_t index = end == 0 ? 0 : piece.size() - 1;
      const std::size_t neighbor = end == 0 ? 1 : piece.size() - 2;
      const Result<Point3d> tip = face_uv_eval(face.projection, piece[index]);
      const Result<Point3d> next = face_uv_eval(face.projection, piece[neighbor]);
      if (!tip.has_value() || !next.has_value()) continue;
      Result<Point3d> hit = curved_boolean_boundary_crossing(face, next.value(), tip.value(),
                                                             curve, &tip.value(), tolerance);
      if (!hit.has_value()) {
        hit = curved_boolean_boundary_crossing(other, next.value(), tip.value(), curve,
                                               &tip.value(), tolerance);
      }
      if (!hit.has_value()) continue;
      const Vec3d gap = hit.value() - tip.value();
      if (!(std::sqrt(dot(gap, gap)) <= junction_tol)) continue;
      const Result<Point2d> projected =
          face_uv_project(face.projection, hit.value(), piece[index][0], true);
      if (projected.has_value()) piece[index] = projected.value();
    }
  }
  // The intersector and `face_uv_project` report a periodic coordinate in
  // (-pi, pi], while the face's own parameter window is [u_lo, u_lo + 2*pi]. A
  // chain that reaches the splitter with u below u_lo reads as "outside the
  // face" there, so the splitter trims its end and tears the ring it should
  // close. Wrapping every sample into the window and splitting at the seam keeps
  // each piece inside the window with its ends exactly on the seam edges.
  if (face.projection.periodic_u()) {
    const double base = face.u[0];
    std::vector<std::vector<Point2d>> windowed;
    windowed.reserve(pieces.size());
    for (const std::vector<Point2d>& piece : pieces) {
      if (piece.size() < 2) continue;
      std::vector<Point2d> wrapped = piece;
      for (Point2d& uv : wrapped) uv[0] = curved_boolean_wrap_periodic(uv[0], base);
      for (std::vector<Point2d>& run : curved_boolean_split_at_seam(wrapped, base, tolerance)) {
        windowed.push_back(std::move(run));
      }
    }
    pieces.swap(windowed);
  }
  return pieces;
}

// Appends the cut when it spans at least the tolerance, and reports whether it
// did.
inline bool curved_boolean_push_cut(std::vector<FaceCut>& cuts, std::vector<Point2d> points,
                                    double tolerance) {
  if (points.size() < 2) return false;
  bool spans = false;
  for (std::size_t i = 1; i < points.size(); ++i) {
    const double du = points[i][0] - points[i - 1][0];
    const double dv = points[i][1] - points[i - 1][1];
    if (du * du + dv * dv > tolerance * tolerance) {
      spans = true;
      break;
    }
  }
  if (!spans) return false;
  cuts.push_back(FaceCut{std::move(points)});
  return true;
}

// A representative point of a fragment, in the face's parameter domain.
//
// The centroid of the largest triangulation triangle is used rather than the
// first ear: a near-degenerate first ear (three nearly collinear ring vertices,
// common where a healed cut arc runs into a fragment) has its centroid on the
// fragment's border, which would classify the fragment by the surface it lies on
// instead of by its body. The whole fragment is triangulated, holes included, so
// a ring-shaped fragment is represented by a point in its band rather than by a
// point of the hole. `split_solid_by_plane` documents the same failure mode and
// makes the same choice.
inline Point2d curved_boolean_representative(const std::vector<std::vector<Point2d>>& loops) {
  Polygon2d polygon;
  polygon.ring_starts.push_back(0u);
  for (const std::vector<Point2d>& loop : loops) {
    polygon.vertices.insert(polygon.vertices.end(), loop.begin(), loop.end());
    polygon.ring_starts.push_back(static_cast<std::uint32_t>(polygon.vertices.size()));
  }
  const std::vector<Point2d>& ring = polygon.vertices;
  const Triangulation2 triangulation = triangulate_polygon2(polygon);
  Point2d representative{0.0, 0.0};
  bool have = false;
  if (triangulation.ok && !triangulation.triangles.empty()) {
    double best_area = -1.0;
    for (const std::array<std::uint32_t, 3>& triangle : triangulation.triangles) {
      const Point2d& a = ring[triangle[0]];
      const Point2d& b = ring[triangle[1]];
      const Point2d& c = ring[triangle[2]];
      const double area = std::abs((b[0] - a[0]) * (c[1] - a[1]) - (b[1] - a[1]) * (c[0] - a[0]));
      if (area <= best_area) continue;
      best_area = area;
      representative = Point2d{0.0, 0.0};
      for (const std::uint32_t index : triangle) {
        representative[0] += ring[index][0];
        representative[1] += ring[index][1];
      }
      representative[0] *= (1.0 / 3.0);
      representative[1] *= (1.0 / 3.0);
      have = true;
    }
  }
  if (!have && !loops.empty()) {
    for (const Point2d& point : loops.front()) {
      representative[0] += point[0];
      representative[1] += point[1];
    }
    if (!loops.front().empty()) {
      const double inverse = 1.0 / static_cast<double>(loops.front().size());
      representative[0] *= inverse;
      representative[1] *= inverse;
    }
  }
  return representative;
}

// Signed area of a loop's parameter polygon. Split fragments are produced in the
// face's (u, v) frame, so a positive value means the ring runs counter-clockwise
// about the surface normal, which is the winding every consumer of
// `loops_uv` (the assembler's pcurves and `curved_solid_volume`'s planar area)
// expects.
inline double curved_boolean_uv_area(const std::vector<Point2d>& ring) noexcept {
  const std::size_t count = ring.size();
  double twice = 0.0;
  for (std::size_t i = 0; i < count; ++i) {
    const std::size_t next = (i + 1) % count;
    twice += ring[i][0] * ring[next][1] - ring[next][0] * ring[i][1];
  }
  return 0.5 * twice;
}

// Reverses every ring of a fragment in place, keeping points, parameters and arc
// spans in step.
inline void curved_boolean_reverse_rings(std::vector<std::vector<Point3d>>& loops,
                                         std::vector<std::vector<Point2d>>& loops_uv,
                                         std::vector<std::vector<std::optional<CurvedArcSpan>>>&
                                             arc_spans) {
  const std::size_t loop_count = loops.size();
  if (loops_uv.size() != loop_count) loops_uv.resize(loop_count);
  if (arc_spans.size() != loop_count) arc_spans.resize(loop_count);
  for (std::size_t loop = 0; loop < loop_count; ++loop) {
    split_solid_reverse_loop(loops[loop], loops_uv[loop], arc_spans[loop]);
  }
}

// The exact circle in which a face pair meets, when one side is a plane and the
// other a cylinder or a sphere with a parallel axis. Both faces of the pair
// receive the same circle object, which is what lets the trimmed fragments and
// the healed arcs agree edge for edge.
inline bool curved_boolean_pair_circle(const SurfaceGeometry& a_surface,
                                       const PlaneFrame& a_frame,
                                       const SurfaceGeometry& b_surface,
                                       const PlaneFrame& b_frame,
                                       CircleCurve3& circle) noexcept {
  const Plane3* plane = std::get_if<Plane3>(&a_surface);
  const SurfaceGeometry* other = &b_surface;
  const PlaneFrame* frame = &a_frame;
  if (plane == nullptr) {
    plane = std::get_if<Plane3>(&b_surface);
    other = &a_surface;
    frame = &b_frame;
  }
  if (plane == nullptr) return false;
  const double length = std::sqrt(plane->normal[0] * plane->normal[0] +
                                  plane->normal[1] * plane->normal[1] +
                                  plane->normal[2] * plane->normal[2]);
  if (!(length > 0.0) || !std::isfinite(length)) return false;
  const Vec3d unit{plane->normal[0] / length, plane->normal[1] / length,
                   plane->normal[2] / length};
  return split_solid_face_circle(*other, unit, plane->offset / length, *frame, circle);
}

// Number of connected components of the kept faces, glued by the 3D points their
// loops share. `assemble_curved_solid` wraps every input face in one shell, so a
// result that falls apart has to be rejected before it is assembled.
inline std::uint32_t curved_boolean_component_count(
    const std::vector<CurvedAssembledFace>& faces, double weld) {
  if (faces.empty()) return 0;
  std::vector<Point3d> points;
  std::vector<std::vector<std::size_t>> face_points(faces.size());
  for (std::size_t f = 0; f < faces.size(); ++f) {
    for (const std::vector<Point3d>& loop : faces[f].loops) {
      for (const Point3d& point : loop) {
        std::size_t slot = points.size();
        for (std::size_t i = 0; i < points.size(); ++i) {
          const Vec3d delta{point[0] - points[i][0], point[1] - points[i][1],
                            point[2] - points[i][2]};
          if (norm2(delta) <= weld) {
            slot = i;
            break;
          }
        }
        if (slot == points.size()) points.push_back(point);
        face_points[f].push_back(slot);
      }
    }
  }

  std::vector<std::size_t> parent(faces.size());
  std::iota(parent.begin(), parent.end(), std::size_t{0});
  auto find = [&parent](std::size_t start) {
    std::size_t root = start;
    while (parent[root] != root) root = parent[root];
    while (parent[start] != root) {
      const std::size_t next = parent[start];
      parent[start] = root;
      start = next;
    }
    return root;
  };
  std::vector<std::size_t> owner(points.size(), faces.size());
  for (std::size_t f = 0; f < faces.size(); ++f) {
    for (const std::size_t slot : face_points[f]) {
      if (owner[slot] == faces.size()) {
        owner[slot] = f;
        continue;
      }
      const std::size_t left = find(f);
      const std::size_t right = find(owner[slot]);
      if (left != right) parent[left] = right;
    }
  }
  std::uint32_t components = 0;
  for (std::size_t f = 0; f < faces.size(); ++f) {
    if (find(f) == f) ++components;
  }
  return components;
}

// Validates one operand and returns its faces. `solid` must be the document's
// only solid, and that solid must own exactly one shell whose faces carry
// supported surfaces.
inline Result<std::vector<FaceId>> curved_boolean_operand_faces(const Brep& brep, SolidId solid) {
  const std::vector<SolidId> solids = brep.solid_ids();
  bool alive = false;
  for (const SolidId candidate : solids) {
    if (candidate == solid) alive = true;
  }
  if (!alive || solids.size() != 1) {
    return Error{Errc::invalid_topology,
                 "boolean_solid_curved: operand must have exactly one solid", solid.index(), 0u};
  }
  const std::span<const ShellId> shells = brep.solid_shells(solid);
  if (shells.size() != 1) {
    return Error{Errc::unsupported, "boolean_solid_curved: operand must have exactly one shell",
                 solid.index(), 0u};
  }
  const ValidationReport topology = validate_topology(brep);
  if (!topology.ok()) {
    const ValidationIssue& issue = topology.issues.front();
    return Error{issue.code, issue.message, issue.entity, 0u};
  }

  std::vector<FaceId> faces;
  for (const FaceId id : brep.shell_faces(shells.front())) {
    const Face& entity = brep.face(id);
    if (entity.surface == kNoGeometry) {
      return Error{Errc::invalid_topology, "boolean_solid_curved: face has no surface geometry",
                   id.index(), 0u};
    }
    const SurfaceGeometry& surface = brep.geometry().surfaces[entity.surface];
    const CurvedBooleanSurfaceGate gate = curved_boolean_surface_gate(surface);
    if (!gate.supported) return Error{Errc::unsupported, gate.message, id.index(), 0u};
    faces.push_back(id);
  }
  return faces;
}

}  // namespace detail

// Boolean combination of two analytic B-rep solids.
//
// `a` / `b` are single-solid documents, `a_solid` / `b_solid` their solids, and
// `op` the regularized operation. The result is a new document with one solid;
// counters in the report are a pure function of the inputs.
inline Result<CurvedBooleanResult> boolean_solid_curved(const Brep& a, SolidId a_solid,
                                                        const Brep& b, SolidId b_solid,
                                                        SolidBooleanOp op,
                                                        const CurvedBooleanOptions& options = {}) {
  if (!detail::curved_boolean_valid_tolerance(options.tolerance)) {
    return Error{Errc::tolerance_violation,
                 "boolean_solid_curved: tolerance must be positive and finite", 0u, 0u};
  }
  const double tolerance = options.tolerance;

  const Result<std::vector<FaceId>> faces_a = detail::curved_boolean_operand_faces(a, a_solid);
  if (!faces_a.has_value()) return faces_a.error();
  const Result<std::vector<FaceId>> faces_b = detail::curved_boolean_operand_faces(b, b_solid);
  if (!faces_b.has_value()) return faces_b.error();

  std::vector<detail::CurvedBooleanFaceInfo> infos_a;
  infos_a.reserve(faces_a.value().size());
  for (const FaceId id : faces_a.value()) {
    const Result<detail::CurvedBooleanFaceInfo> info =
        detail::curved_boolean_face_info(a, id, tolerance);
    if (!info.has_value()) return info.error();
    infos_a.push_back(info.value());
  }
  std::vector<detail::CurvedBooleanFaceInfo> infos_b;
  infos_b.reserve(faces_b.value().size());
  for (const FaceId id : faces_b.value()) {
    const Result<detail::CurvedBooleanFaceInfo> info =
        detail::curved_boolean_face_info(b, id, tolerance);
    if (!info.has_value()) return info.error();
    infos_b.push_back(info.value());
  }

  // A curved fragment whose result interior lies on the side its own surface
  // normal points to must be emitted with a flipped normal, and this driver has
  // no reparameterization for a curved surface (see the harvest lambda below).
  // Decide that from the operation and the operand bounds before any cut work:
  // such a fragment is discarded anyway, and a cut that cannot be traced would
  // otherwise mask the real limitation with a numeric failure.
  {
    const auto operand_bounds =
        [](const std::vector<detail::CurvedBooleanFaceInfo>& infos) {
          detail::CurvedBooleanBounds bounds{};
          bool valid = false;
          for (const detail::CurvedBooleanFaceInfo& info : infos) {
            for (std::size_t axis = 0; axis < 3; ++axis) {
              if (!valid || info.bounds.lo[axis] < bounds.lo[axis]) {
                bounds.lo[axis] = info.bounds.lo[axis];
              }
              if (!valid || info.bounds.hi[axis] > bounds.hi[axis]) {
                bounds.hi[axis] = info.bounds.hi[axis];
              }
            }
            valid = true;
          }
          return bounds;
        };
    const detail::CurvedBooleanBounds bounds_a = operand_bounds(infos_a);
    const detail::CurvedBooleanBounds bounds_b = operand_bounds(infos_b);
    if (detail::curved_boolean_bounds_overlap(bounds_a, bounds_b, tolerance)) {
      const auto flipped_face =
          [op](const Brep& operand,
               const std::vector<detail::CurvedBooleanFaceInfo>& infos,
               bool from_a) -> std::optional<std::uint32_t> {
        for (const detail::CurvedBooleanFaceInfo& info : infos) {
          if (std::holds_alternative<Plane3>(*info.surface)) continue;
          const bool reversed = operand.face(info.id).reversed;
          for (const bool other_in : {false, true}) {
            const bool keep_plus =
                from_a ? detail::solid_boolean_member(false, other_in, op)
                       : detail::solid_boolean_member(other_in, false, op);
            const bool keep_minus =
                from_a ? detail::solid_boolean_member(true, other_in, op)
                       : detail::solid_boolean_member(other_in, true, op);
            if (keep_plus == keep_minus) continue;
            if (reversed == keep_plus) continue;
            return static_cast<std::uint32_t>(info.id.index());
          }
        }
        return std::nullopt;
      };
      std::optional<std::uint32_t> flipped = flipped_face(a, infos_a, true);
      if (!flipped.has_value()) flipped = flipped_face(b, infos_b, false);
      if (flipped.has_value()) {
        return Error{Errc::unsupported,
                     "boolean_solid_curved: flipped curved fragments are not supported yet",
                     flipped.value(), 0u};
      }
    }
  }

  // Membership is decided against a triangulated copy of the other operand.
  const Result<TriMesh> mesh_a = detail::curved_boolean_solid_mesh(a, infos_a, tolerance);
  if (!mesh_a.has_value()) return mesh_a.error();
  const Result<TriMesh> mesh_b = detail::curved_boolean_solid_mesh(b, infos_b, tolerance);
  if (!mesh_b.has_value()) return mesh_b.error();

  double scale = 0.0;
  for (const detail::CurvedBooleanFaceInfo& info : infos_a) {
    for (std::size_t axis = 0; axis < 3; ++axis) {
      scale = std::max(scale, info.bounds.hi[axis] - info.bounds.lo[axis]);
    }
  }
  for (const detail::CurvedBooleanFaceInfo& info : infos_b) {
    for (std::size_t axis = 0; axis < 3; ++axis) {
      scale = std::max(scale, info.bounds.hi[axis] - info.bounds.lo[axis]);
    }
  }
  if (!(scale > 0.0)) scale = 1.0;

  CurvedBooleanReport report;
  SsiOptions ssi_options;
  ssi_options.tolerance = tolerance;

  std::vector<std::vector<FaceCut>> cuts_a(infos_a.size());
  std::vector<std::vector<FaceCut>> cuts_b(infos_b.size());
  std::vector<std::vector<CircleCurve3>> circles_a(infos_a.size());
  std::vector<std::vector<CircleCurve3>> circles_b(infos_b.size());

  for (std::size_t ia = 0; ia < infos_a.size(); ++ia) {
    for (std::size_t ib = 0; ib < infos_b.size(); ++ib) {
      const detail::CurvedBooleanFaceInfo& face_a = infos_a[ia];
      const detail::CurvedBooleanFaceInfo& face_b = infos_b[ib];
      if (!detail::curved_boolean_bounds_overlap(face_a.bounds, face_b.bounds, tolerance)) {
        continue;
      }
      ++report.face_pairs;

      SsiPatch patch_a;
      patch_a.surface = *face_a.surface;
      patch_a.u = face_a.u;
      patch_a.v = face_a.v;
      SsiPatch patch_b;
      patch_b.surface = *face_b.surface;
      patch_b.u = face_b.u;
      patch_b.v = face_b.v;

      const Result<SsiResult> intersected = ssi(patch_a, patch_b, ssi_options);
      if (!intersected.has_value()) {
        return Error{Errc::numeric_failure, "boolean_solid_curved: cut curve did not converge",
                     face_a.id.index(), face_b.id.index()};
      }
      if (!intersected.value().issues.empty()) {
        bool skip_pair = true;
        for (const SsiIssue& issue : intersected.value().issues) {
          switch (issue.code) {
            case SsiDiagnostic::coincident:
              return Error{Errc::unsupported,
                           "boolean_solid_curved: coplanar face contact is not supported yet",
                           face_a.id.index(), face_b.id.index()};
            case SsiDiagnostic::tangent_contact:
              return Error{Errc::unsupported,
                           "boolean_solid_curved: tangent face contact is not supported yet",
                           face_a.id.index(), face_b.id.index()};
            case SsiDiagnostic::no_convergence:
              return Error{Errc::numeric_failure,
                           "boolean_solid_curved: cut curve did not converge",
                           face_a.id.index(), face_b.id.index()};
            case SsiDiagnostic::degenerate_patch:
              break;
            case SsiDiagnostic::none:
              skip_pair = false;
              break;
          }
        }
        if (skip_pair) continue;
      }
      report.max_cut_residual = std::max(report.max_cut_residual, intersected.value().max_residual);

      CircleCurve3 circle;
      bool has_circle = false;
      if (detail::curved_boolean_pair_circle(*face_a.surface, face_a.projection.frame,
                                             *face_b.surface, face_b.projection.frame, circle)) {
        has_circle = true;
        circles_a[ia].push_back(circle);
        circles_b[ib].push_back(circle);
      }

      for (const SsiCurve& curve : intersected.value().curves) {
        std::vector<std::vector<Point2d>> chains_a =
            detail::curved_boolean_cut_chains(curve, true, face_a, tolerance);
        std::vector<std::vector<Point2d>> chains_b =
            detail::curved_boolean_cut_chains(curve, false, face_b, tolerance);
        // Both faces must name the rim by the same vertex sequence. A closed
        // cut on the periodic side starts at the seam crossing, so the other
        // side's copy is re-based onto that crossing point.
        if (curve.closed && chains_a.size() == 1 && chains_b.size() == 1 &&
            face_a.projection.periodic_u() != face_b.projection.periodic_u()) {
          const bool periodic_a = face_a.projection.periodic_u();
          const auto& seam_face = periodic_a ? face_a : face_b;
          const auto& plain_face = periodic_a ? face_b : face_a;
          std::vector<Point2d>& other_chain = periodic_a ? chains_b.front() : chains_a.front();
          const Result<Point3d> rim_start = detail::face_uv_eval(
              seam_face.projection, periodic_a ? chains_a.front().front() : chains_b.front().front());
          if (rim_start.has_value()) {
            (void)detail::curved_boolean_rebase_closed_chain(other_chain, rim_start.value(),
                                                             plain_face.projection, tolerance);
          }
        }
        for (std::vector<Point2d>& chain : chains_a) {
          for (std::vector<Point2d>& piece : detail::curved_boolean_clip_pair(
                   face_a, chain, face_b, has_circle ? &circle : nullptr, tolerance)) {
            if (detail::curved_boolean_push_cut(cuts_a[ia], std::move(piece), tolerance)) {
              ++report.cuts;
            }
          }
        }
        for (std::vector<Point2d>& chain : chains_b) {
          for (std::vector<Point2d>& piece : detail::curved_boolean_clip_pair(
                   face_b, chain, face_a, has_circle ? &circle : nullptr, tolerance)) {
            if (detail::curved_boolean_push_cut(cuts_b[ib], std::move(piece), tolerance)) {
              ++report.cuts;
            }
          }
        }
      }
    }
  }

  struct CurvedBooleanCandidate {
    CurvedAssembledFace face;
    bool from_a = true;
    FaceId source{};
  };
  std::vector<CurvedBooleanCandidate> kept;

  // Splits every face, heals its circular cut chains and keeps the fragments the
  // membership rule puts on the result boundary.
  auto harvest = [&](const Brep& brep, const std::vector<detail::CurvedBooleanFaceInfo>& infos,
                     const std::vector<std::vector<FaceCut>>& cuts,
                     const std::vector<std::vector<CircleCurve3>>& circles, const TriMesh& other,
                     bool from_a) -> Result<void> {
    for (std::size_t index = 0; index < infos.size(); ++index) {
      const detail::CurvedBooleanFaceInfo& info = infos[index];
      const Result<std::vector<CurvedFragment>> split =
          split_face_curved(brep, info.id, cuts[index], tolerance);
      if (!split.has_value()) return split.error();
      for (CurvedFragment fragment : split.value()) {
        if (fragment.loops.empty()) continue;
        ++report.fragments;
        const Result<void> healed =
            detail::split_solid_heal_fragment(fragment, circles[index], tolerance);
        if (!healed.has_value()) return healed.error();

        const Point2d representative =
            detail::curved_boolean_representative(fragment.loops_uv);
        const Point3d point = detail::face_uv_eval(info.projection, representative);
        const Result<bool> inside = point_in_tessellation(other, point, tolerance);
        if (!inside.has_value()) {
          if (inside.error().code == Errc::numeric_failure) {
            return Error{Errc::numeric_failure,
                         "boolean_solid_curved: fragment representative lies on the other solid",
                         info.id.index(), 0u};
          }
          return inside.error();
        }
        const bool other_in = inside.value();
        const bool keep_plus =
            from_a ? detail::solid_boolean_member(false, other_in, op)
                   : detail::solid_boolean_member(other_in, false, op);
        const bool keep_minus =
            from_a ? detail::solid_boolean_member(true, other_in, op)
                   : detail::solid_boolean_member(other_in, true, op);
        if (keep_plus == keep_minus) {
          ++report.dropped;
          continue;
        }

        CurvedAssembledFace assembled;
        assembled.surface = *info.surface;
        assembled.loops = std::move(fragment.loops);
        assembled.loops_uv = std::move(fragment.loops_uv);
        assembled.arc_spans = std::move(fragment.arc_spans);
        assembled.reversed = false;
        if (assembled.arc_spans.size() != assembled.loops.size()) {
          assembled.arc_spans.resize(assembled.loops.size());
        }
        if (assembled.loops_uv.size() != assembled.loops.size()) {
          assembled.loops_uv.resize(assembled.loops.size());
        }
        // Ring winding invariant: counter-clockwise about the surface normal in
        // the face's (u, v) frame. `split_face_curved` inherits the winding of
        // the source face's loops, which is only guaranteed for an unreversed
        // face, so a reversed fragment is normalized here rather than handed to
        // the assembler and to the volume integral as-is.
        if (!assembled.loops_uv.empty() &&
            detail::curved_boolean_uv_area(assembled.loops_uv.front()) < 0.0) {
          detail::curved_boolean_reverse_rings(assembled.loops, assembled.loops_uv,
                                               assembled.arc_spans);
        }
        // A fragment whose result interior lies on the side its own surface
        // normal points to needs that normal flipped, exactly as
        // `detail::solid_boolean_classify_face` negates the plane it emits. The
        // `reversed` flag alone is not enough: two faces sharing a cut edge must
        // walk it in opposite directions, so the rings are reversed too, and a
        // planar face is re-expressed with the flipped normal baked into its
        // plane (parameters reprojected into the new frame) so that `reversed`
        // stays false and `curved_solid_volume` never sees a negatively oriented
        // planar ring. Curved surfaces have no such reparameterization here, so
        // they are deferred instead of being emitted with a winding the volume
        // integral would misread.
        const bool flip = brep.face(info.id).reversed != keep_plus;
        if (flip) {
          const Plane3* plane = std::get_if<Plane3>(&assembled.surface);
          if (plane == nullptr) {
            return Error{Errc::unsupported,
                         "boolean_solid_curved: flipped curved fragments are not supported yet",
                         info.id.index(), 0u};
          }
          detail::curved_boolean_reverse_rings(assembled.loops, assembled.loops_uv,
                                               assembled.arc_spans);
          const Plane3 flipped{
              Vec3d{-plane->normal[0], -plane->normal[1], -plane->normal[2]}, -plane->offset};
          assembled.surface = SurfaceGeometry{flipped};
          const Result<PlaneFrame> frame = make_plane_frame(flipped);
          if (!frame.has_value()) return frame.error();
          detail::FaceUvProjection projection;
          projection.surface = &assembled.surface;
          projection.frame = frame.value();
          projection.tolerance = tolerance;
          std::vector<std::vector<Point2d>> reprojected(assembled.loops.size());
          for (std::size_t loop = 0; loop < assembled.loops.size(); ++loop) {
            reprojected[loop].reserve(assembled.loops[loop].size());
            double previous_u = 0.0;
            bool has_previous = false;
            for (const Point3d& surface_point : assembled.loops[loop]) {
              const Result<Point2d> uv =
                  detail::face_uv_project(projection, surface_point, previous_u, has_previous);
              if (!uv.has_value()) return uv.error();
              reprojected[loop].push_back(uv.value());
              previous_u = uv.value()[0];
              has_previous = true;
            }
          }
          assembled.loops_uv = std::move(reprojected);
        }
        detail::split_solid_canonicalize_face(assembled, tolerance);
        kept.push_back(CurvedBooleanCandidate{std::move(assembled), from_a, info.id});
        if (flip) ++report.flipped;
      }
    }
    return Result<void>{};
  };

  const Result<void> harvested_a =
      harvest(a, infos_a, cuts_a, circles_a, mesh_b.value(), true);
  if (!harvested_a.has_value()) return harvested_a.error();
  const Result<void> harvested_b =
      harvest(b, infos_b, cuts_b, circles_b, mesh_a.value(), false);
  if (!harvested_b.has_value()) return harvested_b.error();

  if (kept.empty()) {
    return Error{Errc::numeric_failure, "boolean_solid_curved: result is empty", 0u, 0u};
  }
  std::stable_sort(kept.begin(), kept.end(),
                   [](const CurvedBooleanCandidate& left, const CurvedBooleanCandidate& right) {
                     if (left.from_a != right.from_a) return left.from_a;
                     return left.source.index() < right.source.index();
                   });

  report.kept = static_cast<std::uint32_t>(kept.size());
  std::vector<CurvedAssembledFace> assembled_faces;
  assembled_faces.reserve(kept.size());
  for (CurvedBooleanCandidate& candidate : kept) {
    assembled_faces.push_back(std::move(candidate.face));
  }

  const double weld = std::max(64.0 * tolerance, 1e-9 * scale);
  const std::size_t components = detail::curved_boolean_component_count(assembled_faces, weld);
  if (components > 1) {
    return Error{Errc::unsupported,
                 "boolean_solid_curved: result is not a single connected component", 0u, 0u};
  }

  Result<Brep> assembled = assemble_curved_solid(assembled_faces, tolerance);
  if (!assembled.has_value()) {
    // A result that touches itself along an edge (for example the symmetric
    // difference of two boxes that share a corner edge) is non-manifold: the
    // assembler reports the shared edge as used once because it cannot decide
    // how many face uses it needs. That is a documented limit of this driver,
    // not a bug to work around, so report it as unsupported.
    if (assembled.error().code == Errc::invalid_topology &&
        assembled.error().message != nullptr &&
        std::string_view{assembled.error().message} ==
            "assemble_curved_solid: edge is used once") {
      return Error{Errc::unsupported,
                   "boolean_solid_curved: result touches itself along an edge", 0u, 0u};
    }
    return assembled.error();
  }
  if (!validate_topology(assembled.value()).ok()) {
    return Error{Errc::numeric_failure,
                 "boolean_solid_curved: assembled result failed validation", 0u, 0u};
  }

  CurvedBooleanResult result;
  result.brep = std::move(assembled.value());
  result.report = report;
  return result;
}

}  // namespace v1
}  // namespace geom

