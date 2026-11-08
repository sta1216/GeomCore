#pragma once

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <span>
#include <utility>
#include <variant>
#include <vector>

#include <geom/brep/topology.hpp>
#include <geom/cgeom/polygon2.hpp>
#include <geom/core/result.hpp>
#include <geom/cgeom/triangulate2.hpp>
#include <geom/curves/analytic_curve.hpp>
#include <geom/primitives/point.hpp>
#include <geom/primitives/vec.hpp>
#include <geom/surfaces/plane.hpp>

namespace geom {
inline namespace v1 {

// Triangulated view of a B-rep. `positions` is the welded vertex pool, every
// `triangles` entry indexes it, and `triangle_face[i]` names the face that
// produced `triangles[i]`. The weld guarantees that two coedges of the same edge
// emit the same position indices, so a mesh built from a closed solid has no
// cracks along shared edges.
struct TriMesh {
  std::vector<Point3d> positions;
  std::vector<std::array<std::uint32_t, 3>> triangles;
  std::vector<std::uint32_t> triangle_face;
};

struct TessellationOptions {
  // Maximum distance between a chord and the true curve; a lower value emits more
  // circle edge segments.
  double deflection = 1e-3;
  // Reserved for curved faces, which are not tessellated yet; it must still be a
  // valid positive number so callers can pass one option set everywhere.
  double angle_tolerance = 0.35;
};

namespace detail {

inline bool tessellate_point_finite(const Point3d& point) noexcept {
  return std::isfinite(point[0]) && std::isfinite(point[1]) && std::isfinite(point[2]);
}

inline bool tessellate_vec_finite(const Vec3d& value) noexcept {
  return std::isfinite(value[0]) && std::isfinite(value[1]) && std::isfinite(value[2]);
}

// Samples of one edge in edge-parameter order. The cache is keyed by the edge slot
// so a coedge and its twin on the neighbouring face share the exact same
// positions instead of sampling the curve twice.
struct TessellateEdgeSamples {
  std::vector<std::uint32_t> positions;
  bool sampled = false;
};

struct TessellateState {
  TriMesh mesh;
  std::vector<TessellateEdgeSamples> edges;
  std::vector<std::uint32_t> vertices;
  std::vector<char> vertex_sampled;
};

inline Result<std::uint32_t> tessellate_append_point(TessellateState& state,
                                                     const Point3d& point) {
  if (!tessellate_point_finite(point)) {
    return Error{Errc::non_finite, "tessellate: sample position is not finite", 0u, 0u};
  }
  state.mesh.positions.push_back(point);
  return static_cast<std::uint32_t>(state.mesh.positions.size() - 1u);
}

// Position index of a topological vertex. All edges incident to the vertex reuse
// this one index, which is what welds endpoint samples across faces.
inline Result<std::uint32_t> tessellate_vertex_position(TessellateState& state, const Brep& brep,
                                                        VertexId id) {
  if (id.index() >= state.vertices.size()) {
    return Error{Errc::index_out_of_range, "tessellate: vertex index is out of range", id.index(),
                 0u};
  }
  if (state.vertex_sampled[id.index()] == 0) {
    const Result<std::uint32_t> index = tessellate_append_point(state, brep.vertex(id).point);
    if (!index.has_value()) return index.error();
    state.vertices[id.index()] = index.value();
    state.vertex_sampled[id.index()] = 1;
  }
  return state.vertices[id.index()];
}

// Number of chord segments needed to keep the sagitta of a circle edge within
// `deflection`. The sagitta of a step `theta` is radius * (1 - cos(theta / 2)), so
// the largest admissible step is 2 * acos(1 - deflection / radius). The clamp
// keeps acos defined when deflection >= 2 * radius, where the half turn is the
// natural limit.
inline Result<std::uint32_t> tessellate_circle_segments(const CircleCurve3& circle, double t0,
                                                        double t1, double deflection) {
  if (!tessellate_point_finite(circle.center) || !tessellate_vec_finite(circle.x_axis) ||
      !tessellate_vec_finite(circle.y_axis) || !std::isfinite(circle.radius)) {
    return Error{Errc::non_finite, "tessellate: circle edge geometry is not finite", 0u, 0u};
  }
  if (!(circle.radius > 0.0)) {
    return Error{Errc::degenerate, "tessellate: circle edge radius must be positive", 0u, 0u};
  }
  const double span = t1 - t0;
  if (!std::isfinite(span) || !(span > 0.0)) {
    return Error{Errc::degenerate, "tessellate: circle edge interval must be positive", 0u, 0u};
  }
  const double cosine = 1.0 - deflection / circle.radius;
  const double clamped = cosine < -1.0 ? -1.0 : cosine;
  const double step = 2.0 * std::acos(clamped);
  double count = std::ceil(span / step);
  if (!(count >= 2.0)) count = 2.0;
  constexpr double kMaximumCount =
      static_cast<double>(std::numeric_limits<std::uint32_t>::max());
  if (!(count <= kMaximumCount)) {
    return Error{Errc::numeric_failure, "tessellate: circle segmentation count is too large", 0u,
                 0u};
  }
  return static_cast<std::uint32_t>(count);
}

// Samples an edge in edge-parameter order. A line or an edge without curve
// geometry contributes its two bounding vertex positions; a circle contributes
// its endpoint plus the interior chord samples. A closed circle edge (both
// bounding vertices are the same) stops at the last interior sample so the seam
// vertex is not repeated.
inline Result<std::vector<Point3d>> tessellate_sample_edge(const Brep& brep, EdgeId id,
                                                           double deflection) {
  const Edge& edge = brep.edge(id);
  const Point3d start = brep.vertex(edge.v0).point;
  const Point3d end = brep.vertex(edge.v1).point;
  const bool closed = edge.v0 == edge.v1;

  std::vector<Point3d> samples;
  if (edge.curve == kNoGeometry || edge.curve >= brep.geometry().curves.size()) {
    samples.push_back(start);
    samples.push_back(end);
    return samples;
  }

  const CurveGeometry& curve = brep.geometry().curves[edge.curve];
  if (std::holds_alternative<LineCurve3>(curve)) {
    samples.push_back(start);
    samples.push_back(end);
    return samples;
  }
  const CircleCurve3* circle = std::get_if<CircleCurve3>(&curve);
  if (circle == nullptr) {
    return Error{Errc::unsupported, "tessellate: only line and circle edges are supported",
                 id.index(), 0u};
  }

  const Result<std::uint32_t> segments =
      tessellate_circle_segments(*circle, edge.t0, edge.t1, deflection);
  if (!segments.has_value()) return segments.error();
  const std::uint32_t count = segments.value();
  const double step = (edge.t1 - edge.t0) / static_cast<double>(count);
  if (!std::isfinite(step)) {
    return Error{Errc::degenerate, "tessellate: circle edge step is not finite", id.index(), 0u};
  }
  samples.reserve(closed ? count : count + 1u);
  samples.push_back(start);
  for (std::uint32_t k = 1; k < count; ++k) {
    samples.push_back(circle->eval(edge.t0 + step * static_cast<double>(k)));
  }
  if (!closed) samples.push_back(end);
  return samples;
}

// Position indices of one edge in edge-parameter order, sampled at most once per
// edge. Endpoints are welded by vertex, interior samples by (edge slot, index).
inline Result<std::vector<std::uint32_t>> tessellate_edge_positions(TessellateState& state,
                                                                    const Brep& brep, EdgeId id,
                                                                    double deflection) {
  if (id.index() >= state.edges.size()) {
    return Error{Errc::index_out_of_range, "tessellate: edge index is out of range", id.index(),
                 0u};
  }
  TessellateEdgeSamples& slot = state.edges[id.index()];
  if (!slot.sampled) {
    const Result<std::vector<Point3d>> samples = tessellate_sample_edge(brep, id, deflection);
    if (!samples.has_value()) return samples.error();
    const std::vector<Point3d>& points = samples.value();
    const Edge& edge = brep.edge(id);
    const bool closed = edge.v0 == edge.v1;

    const Result<std::uint32_t> first = tessellate_vertex_position(state, brep, edge.v0);
    if (!first.has_value()) return first.error();
    slot.positions.push_back(first.value());
    for (std::size_t k = 1; k + 1 < points.size(); ++k) {
      const Result<std::uint32_t> index = tessellate_append_point(state, points[k]);
      if (!index.has_value()) return index.error();
      slot.positions.push_back(index.value());
    }
    if (closed) {
      if (!points.empty()) {
        const Result<std::uint32_t> index = tessellate_append_point(state, points.back());
        if (!index.has_value()) return index.error();
        slot.positions.push_back(index.value());
      }
    } else {
      const Result<std::uint32_t> last = tessellate_vertex_position(state, brep, edge.v1);
      if (!last.has_value()) return last.error();
      slot.positions.push_back(last.value());
    }
    slot.sampled = true;
  }
  return slot.positions;
}

// One loop of a planar face after projection into the face frame.
struct TessellateRing {
  std::vector<Point2d> points;
  std::vector<std::uint32_t> positions;
  bool outer = false;
};

inline Result<void> tessellate_planar_face(TessellateState& state, const Brep& brep, FaceId id,
                                           const TessellationOptions& options) {
  const Face& face = brep.face(id);
  if (face.surface >= brep.geometry().surfaces.size()) {
    return Error{Errc::unsupported, "tessellate: only planar faces are supported", id.index(), 0u};
  }
  const Plane3* plane = std::get_if<Plane3>(&brep.geometry().surfaces[face.surface]);
  if (plane == nullptr) {
    return Error{Errc::unsupported, "tessellate: only planar faces are supported", id.index(), 0u};
  }
  if (!std::isfinite(plane->offset)) {
    return Error{Errc::non_finite, "tessellate: plane offset is not finite", id.index(), 0u};
  }
  if (!tessellate_vec_finite(plane->normal)) {
    return Error{Errc::non_finite, "tessellate: plane normal is not finite", id.index(), 0u};
  }
  const double normal_length = std::hypot(plane->normal[0], plane->normal[1], plane->normal[2]);
  if (!(normal_length > 0.0)) {
    return Error{Errc::degenerate, "tessellate: plane normal is zero", id.index(), 0u};
  }
  const Vec3d unit{plane->normal[0] / normal_length, plane->normal[1] / normal_length,
                   plane->normal[2] / normal_length};

  // Orthonormal frame with cross(u_axis, v_axis) == unit: a ring that is CCW in
  // (u, v) has world normal +unit. The final swap mirrors that direction for a
  // face that is stored reversed.
  const Vec3d reference =
      std::fabs(unit[0]) < 0.9 ? Vec3d{1.0, 0.0, 0.0} : Vec3d{0.0, 1.0, 0.0};
  const Vec3d raw_u = cross(reference, unit);
  const double u_length = std::hypot(raw_u[0], raw_u[1], raw_u[2]);
  if (!(u_length > 0.0)) {
    return Error{Errc::degenerate, "tessellate: plane frame is degenerate", id.index(), 0u};
  }
  const Vec3d u_axis{raw_u[0] / u_length, raw_u[1] / u_length, raw_u[2] / u_length};
  const Vec3d v_axis = cross(unit, u_axis);
  const double offset = plane->offset / normal_length;
  const Point3d origin{unit[0] * offset, unit[1] * offset, unit[2] * offset};

  std::vector<TessellateRing> rings;
  const std::size_t loop_limit = brep.loops().size();
  const std::size_t coedge_limit = brep.coedges().size();
  LoopId loop_id = face.first_loop;
  std::size_t loop_count = 0;
  while (!loop_id.is_null()) {
    if (loop_count++ > loop_limit) {
      return Error{Errc::invalid_topology, "tessellate: loop chain is cyclic", id.index(), 0u};
    }
    const Loop& loop = brep.loop(loop_id);

    TessellateRing ring;
    ring.outer = loop.outer;
    CoedgeId coedge_id = loop.start;
    std::size_t coedge_count = 0;
    if (coedge_id.is_null()) {
      return Error{Errc::degenerate, "tessellate: loop has no coedges", id.index(), 0u};
    }
    do {
      if (coedge_count++ > coedge_limit) {
        return Error{Errc::invalid_topology, "tessellate: coedge cycle is cyclic", id.index(), 0u};
      }
      const Coedge& coedge = brep.coedge(coedge_id);
      const Result<std::vector<std::uint32_t>> samples =
          tessellate_edge_positions(state, brep, coedge.edge, options.deflection);
      if (!samples.has_value()) return samples.error();
      const std::vector<std::uint32_t>& order = samples.value();
      if (coedge.reversed) {
        for (std::size_t k = order.size(); k > 0; --k) {
          const std::uint32_t index = order[k - 1];
          if (ring.positions.empty() || ring.positions.back() != index) {
            ring.positions.push_back(index);
          }
        }
      } else {
        for (const std::uint32_t index : order) {
          if (ring.positions.empty() || ring.positions.back() != index) {
            ring.positions.push_back(index);
          }
        }
      }
      coedge_id = coedge.next;
    } while (!coedge_id.is_null() && coedge_id != loop.start);
    if (coedge_id.is_null()) {
      return Error{Errc::invalid_topology, "tessellate: loop is not closed", id.index(), 0u};
    }
    if (!ring.positions.empty() && ring.positions.front() == ring.positions.back()) {
      ring.positions.pop_back();
    }
    if (ring.positions.size() < 3) {
      return Error{Errc::degenerate, "tessellate: loop has fewer than three samples",
                   id.index(), 0u};
    }
    rings.push_back(std::move(ring));
    loop_id = loop.next;
  }

  const TessellateRing* outer_ring = nullptr;
  for (const TessellateRing& ring : rings) {
    if (!ring.outer) continue;
    if (outer_ring != nullptr) {
      return Error{Errc::invalid_topology, "tessellate: face has more than one outer loop",
                   id.index(), 0u};
    }
    outer_ring = &ring;
  }
  if (outer_ring == nullptr) {
    return Error{Errc::degenerate, "tessellate: face has no outer loop", id.index(), 0u};
  }

  Polygon2d polygon;
  polygon.ring_starts.push_back(0u);
  std::vector<std::uint32_t> polygon_positions;
  const auto append_ring = [&](const TessellateRing& ring, bool want_outer) -> Errc {
    std::vector<Point2d> points;
    std::vector<std::uint32_t> positions = ring.positions;
    points.reserve(positions.size());
    for (const std::uint32_t index : positions) {
      const Vec3d relative = state.mesh.positions[index] - origin;
      const Point2d projected{dot(relative, u_axis), dot(relative, v_axis)};
      if (!(std::isfinite(projected[0]) && std::isfinite(projected[1]))) {
        return Errc::non_finite;
      }
      points.push_back(projected);
    }
    const int orientation = ring_orientation2(points);
    if (orientation == 0) return Errc::degenerate;
    if ((orientation > 0) != want_outer) {
      std::reverse(points.begin(), points.end());
      std::reverse(positions.begin(), positions.end());
    }
    for (std::size_t k = 0; k < points.size(); ++k) {
      polygon.vertices.push_back(points[k]);
      polygon_positions.push_back(positions[k]);
    }
    polygon.ring_starts.push_back(static_cast<std::uint32_t>(polygon.vertices.size()));
    return Errc::none;
  };

  Errc status = append_ring(*outer_ring, true);
  if (status != Errc::none) {
    return Error{status, "tessellate: outer loop cannot be projected", id.index(), 0u};
  }
  for (const TessellateRing& ring : rings) {
    if (&ring == outer_ring) continue;
    status = append_ring(ring, false);
    if (status != Errc::none) {
      return Error{status, "tessellate: hole loop cannot be projected", id.index(), 0u};
    }
  }

  const double area = polygon_area2(polygon);
  if (!std::isfinite(area) || !(area > 0.0)) {
    return Error{Errc::degenerate, "tessellate: planar face has no area", id.index(), 0u};
  }

  const Triangulation2 triangulation = triangulate_polygon2(polygon);
  if (!triangulation.ok) {
    return Error{Errc::degenerate, "tessellate: planar face triangulation failed", id.index(), 0u};
  }
  for (const std::array<std::uint32_t, 3>& triangle : triangulation.triangles) {
    std::array<std::uint32_t, 3> mapped{polygon_positions[triangle[0]],
                                        polygon_positions[triangle[1]],
                                        polygon_positions[triangle[2]]};
    if (face.reversed) std::swap(mapped[1], mapped[2]);
    state.mesh.triangles.push_back(mapped);
    state.mesh.triangle_face.push_back(id.index());
  }
  return {};
}

}  // namespace detail

// Triangulates the planar faces of `brep`. Every edge is sampled once in
// edge-parameter order and the samples are shared by both coedges, so the mesh
// welds along shared edges. Faces on any other surface type are rejected with
// Errc::unsupported until curved-face tessellation lands.
inline Result<TriMesh> tessellate(const Brep& brep, const TessellationOptions& options) {
  if (!std::isfinite(options.deflection) || !(options.deflection > 0.0)) {
    return Error{Errc::tolerance_violation, "tessellate: deflection must be finite and positive",
                 0u, 0u};
  }
  if (!std::isfinite(options.angle_tolerance) || !(options.angle_tolerance > 0.0)) {
    return Error{Errc::tolerance_violation,
                 "tessellate: angle tolerance must be finite and positive", 0u, 0u};
  }

  const std::vector<FaceId> faces = brep.face_ids();
  for (const FaceId id : faces) {
    const Face& face = brep.face(id);
    if (face.surface >= brep.geometry().surfaces.size() ||
        !std::holds_alternative<Plane3>(brep.geometry().surfaces[face.surface])) {
      return Error{Errc::unsupported, "tessellate: only planar faces are supported", id.index(),
                   0u};
    }
  }

  detail::TessellateState state;
  state.edges.resize(brep.edges().size());
  state.vertices.resize(brep.vertices().size());
  state.vertex_sampled.assign(brep.vertices().size(), 0);
  for (const FaceId id : faces) {
    const Result<void> status = detail::tessellate_planar_face(state, brep, id, options);
    if (!status.has_value()) return status.error();
  }
  return std::move(state.mesh);
}

}  // namespace v1
}  // namespace geom
