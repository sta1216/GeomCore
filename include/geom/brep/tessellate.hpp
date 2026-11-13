#pragma once

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <optional>
#include <span>
#include <utility>
#include <variant>
#include <vector>

#include <geom/brep/face_uv.hpp>
#include <geom/brep/topology.hpp>
#include <geom/cgeom/polygon2.hpp>
#include <geom/cgeom/trimmed_domain.hpp>
#include <geom/core/result.hpp>
#include <geom/cgeom/triangulate2.hpp>
#include <geom/curves/analytic_curve.hpp>
#include <geom/primitives/point.hpp>
#include <geom/primitives/vec.hpp>
#include <geom/surfaces/analytic.hpp>
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

constexpr double kPi = 3.14159265358979323846;

inline bool tessellate_point_finite(const Point3d& point) noexcept {
  return std::isfinite(point[0]) && std::isfinite(point[1]) && std::isfinite(point[2]);
}

inline bool tessellate_vec_finite(const Vec3d& value) noexcept {
  return std::isfinite(value[0]) && std::isfinite(value[1]) && std::isfinite(value[2]);
}

// Endpoint agreement between a sampled curve and the edge's vertices. The mesh
// is welded through the vertex pool, so a curve that does not meet its vertices
// cannot contribute a consistent boundary ring.
inline bool tessellate_point_close(const Point3d& a, const Point3d& b) noexcept {
  const double scale = std::max({1.0, std::fabs(a[0]), std::fabs(a[1]), std::fabs(a[2]),
                                 std::fabs(b[0]), std::fabs(b[1]), std::fabs(b[2])});
  return std::fabs(a[0] - b[0]) <= 1e-9 * scale && std::fabs(a[1] - b[1]) <= 1e-9 * scale &&
         std::fabs(a[2] - b[2]) <= 1e-9 * scale;
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
// keeps acos defined when deflection >= 2 * radius, where a half turn already
// meets the tolerance; a closed edge still needs three samples so the resulting
// disk remains a valid triangle, while an open arc only needs its two segments.
inline Result<std::uint32_t> tessellate_circle_segments(const CircleCurve3& circle, double t0,
                                                        double t1, double deflection,
                                                        bool closed) {
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
  const double minimum = closed ? 3.0 : 2.0;
  if (!(count >= minimum)) count = minimum;
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
    const NurbsCurve* nurbs = std::get_if<NurbsCurve>(&curve);
    if (nurbs == nullptr || nurbs->degree() != 1 || !nurbs->is_valid() ||
        !std::isfinite(edge.t0) || !std::isfinite(edge.t1) || !(edge.t1 > edge.t0)) {
      return Error{Errc::unsupported, "tessellate: only line and circle edges are supported",
                   id.index(), 0u};
    }
    // The curved boolean stores every split as a degree-1 polyline, so sampling
    // at the uniform control parameters reproduces the control polygon exactly.
    // A curve whose endpoints disagree with the edge vertices cannot be welded
    // and keeps the historical unsupported result.
    if (!tessellate_point_close(nurbs->eval(edge.t0), start) ||
        !tessellate_point_close(nurbs->eval(edge.t1), end)) {
      return Error{Errc::unsupported, "tessellate: only line and circle edges are supported",
                   id.index(), 0u};
    }
    const std::size_t count = std::max<std::size_t>(2, nurbs->control_count());
    const std::size_t sample_count = closed ? count - 1 : count;
    samples.reserve(sample_count);
    for (std::size_t k = 0; k < sample_count; ++k) {
      if (k == 0) {
        samples.push_back(start);
        continue;
      }
      if (!closed && k + 1 == count) {
        samples.push_back(end);
        continue;
      }
      const double fraction = static_cast<double>(k) / static_cast<double>(count - 1);
      samples.push_back(nurbs->eval(edge.t0 + (edge.t1 - edge.t0) * fraction));
    }
    return samples;
  }

  const Result<std::uint32_t> segments =
      tessellate_circle_segments(*circle, edge.t0, edge.t1, deflection, closed);
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

// Point layout of the concatenated boundary rings of a face.
struct TessellateTrimmedRingLayout {
  std::vector<std::uint32_t> ring_of;
  std::vector<std::uint32_t> ring_begin;
  std::vector<std::uint32_t> ring_count;
};

// Returns the ring points strictly between `from` and `to` when every one of
// them is exactly collinear with the chord and strictly inside it. `std::nullopt`
// means the ring path is not a straight run. The ear clipper removes exactly
// collinear uv vertices, so this recovers the boundary samples it merged.
inline std::optional<std::vector<std::uint32_t>> tessellate_trimmed_run(
    const std::vector<Point2d>& points, std::uint32_t begin, std::uint32_t count,
    std::uint32_t from, std::uint32_t to) {
  const Point2d& a = points[from];
  const Point2d& b = points[to];
  const double dx = b[0] - a[0];
  const double dy = b[1] - a[1];
  const double length2 = dx * dx + dy * dy;
  if (!(length2 > 1e-30)) return std::nullopt;
  std::vector<std::uint32_t> run;
  std::uint32_t offset = (from - begin + 1u) % count;
  for (std::uint32_t step = 0; step + 1u < count; ++step) {
    if (begin + offset == to) return run;
    const std::uint32_t index = begin + offset;
    const Point2d& point = points[index];
    if (sign_orient2d(a[0], a[1], b[0], b[1], point[0], point[1]) != 0) {
      return std::nullopt;
    }
    const double t = ((point[0] - a[0]) * dx + (point[1] - a[1]) * dy) / length2;
    if (!(t > 0.0 && t < 1.0)) return std::nullopt;
    run.push_back(index);
    offset = (offset + 1u) % count;
  }
  return std::nullopt;
}

// The ear clipper merges exactly collinear uv boundary runs into one chord. The
// neighbouring faces keep every sample of a shared trim curve, so each chord that
// spans a straight boundary run is fanned back into the individual sub-segments.
inline Result<std::vector<std::array<std::uint32_t, 3>>> tessellate_expand_trimmed_runs(
    const std::vector<Point2d>& points, const TessellateTrimmedRingLayout& layout,
    const std::vector<std::array<std::uint32_t, 3>>& triangles, FaceId id) {
  std::vector<std::array<std::uint32_t, 3>> expanded;
  std::vector<std::array<std::uint32_t, 3>> pending = triangles;
  while (!pending.empty()) {
    const std::array<std::uint32_t, 3> triangle = pending.back();
    pending.pop_back();
    bool split = false;
    for (std::uint32_t edge = 0; edge < 3u && !split; ++edge) {
      const std::uint32_t a = triangle[edge];
      const std::uint32_t b = triangle[(edge + 1u) % 3u];
      const std::uint32_t c = triangle[(edge + 2u) % 3u];
      if (layout.ring_of[a] != layout.ring_of[b]) continue;
      const std::uint32_t ring = layout.ring_of[a];
      const std::uint32_t begin = layout.ring_begin[ring];
      const std::uint32_t ring_count = layout.ring_count[ring];
      std::optional<std::vector<std::uint32_t>> run =
          tessellate_trimmed_run(points, begin, ring_count, a, b);
      if (run.has_value() && run->empty()) run.reset();
      if (!run.has_value()) {
        run = tessellate_trimmed_run(points, begin, ring_count, b, a);
      }
      if (!run.has_value() || run->empty()) continue;
      const Point2d& pa = points[a];
      const Point2d& pb = points[b];
      const double dx = pb[0] - pa[0];
      const double dy = pb[1] - pa[1];
      std::sort(run->begin(), run->end(), [&](std::uint32_t lhs, std::uint32_t rhs) {
        const Point2d& pl = points[lhs];
        const Point2d& pr = points[rhs];
        return ((pl[0] - pa[0]) * dx + (pl[1] - pa[1]) * dy) <
               ((pr[0] - pa[0]) * dx + (pr[1] - pa[1]) * dy);
      });
      std::uint32_t previous = a;
      for (const std::uint32_t point : *run) {
        // The run points lie on the chord, so every fan triangle keeps the
        // winding of the triangle it replaces.
        if (sign_orient2d(points[previous][0], points[previous][1], points[point][0],
                          points[point][1], points[c][0], points[c][1]) <= 0) {
          return Error{Errc::degenerate, "tessellate: boundary fan is not convex", id.index(), 0u};
        }
        pending.push_back({previous, point, c});
        previous = point;
      }
      pending.push_back({previous, b, c});
      split = true;
    }
    if (!split) expanded.push_back(triangle);
  }
  return expanded;
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
  TessellateTrimmedRingLayout layout;
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
    layout.ring_begin.push_back(static_cast<std::uint32_t>(polygon.vertices.size()));
    layout.ring_count.push_back(static_cast<std::uint32_t>(points.size()));
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
  layout.ring_of.resize(polygon.vertices.size());
  for (std::uint32_t ring = 0; ring < layout.ring_begin.size(); ++ring) {
    for (std::uint32_t sample = 0; sample < layout.ring_count[ring]; ++sample) {
      layout.ring_of[layout.ring_begin[ring] + sample] = ring;
    }
  }
  const Result<std::vector<std::array<std::uint32_t, 3>>> expanded = tessellate_expand_trimmed_runs(
      polygon.vertices, layout, triangulation.triangles, id);
  if (!expanded.has_value()) return expanded.error();
  for (const std::array<std::uint32_t, 3>& triangle : expanded.value()) {
    std::array<std::uint32_t, 3> mapped{polygon_positions[triangle[0]],
                                        polygon_positions[triangle[1]],
                                        polygon_positions[triangle[2]]};
    if (face.reversed) std::swap(mapped[1], mapped[2]);
    state.mesh.triangles.push_back(mapped);
    state.mesh.triangle_face.push_back(id.index());
  }
  return {};
}

// One boundary edge of a face with the welded position indices of its samples in
// edge-parameter order. Curved faces derive their grids from these samples, so a
// seam column or a cap ring is the same sequence of vertices the adjacent face
// sees and the mesh stays watertight.
struct TessellateBoundaryEdge {
  EdgeId edge{};
  std::vector<std::uint32_t> positions;
};

// Collects every distinct boundary edge of a face once, in loop order.
inline Result<std::vector<TessellateBoundaryEdge>> tessellate_boundary_edges(
    TessellateState& state, const Brep& brep, FaceId id, double deflection) {
  const Face& face = brep.face(id);
  std::vector<TessellateBoundaryEdge> edges;
  const std::size_t loop_limit = brep.loops().size();
  const std::size_t coedge_limit = brep.coedges().size();
  LoopId loop_id = face.first_loop;
  std::size_t loop_count = 0;
  if (loop_id.is_null()) {
    return Error{Errc::degenerate, "tessellate: face has no loop", id.index(), 0u};
  }
  while (!loop_id.is_null()) {
    if (loop_count++ > loop_limit) {
      return Error{Errc::invalid_topology, "tessellate: loop chain is cyclic", id.index(), 0u};
    }
    const Loop& loop = brep.loop(loop_id);
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
      bool known = false;
      for (const TessellateBoundaryEdge& edge : edges) {
        if (edge.edge == coedge.edge) {
          known = true;
          break;
        }
      }
      if (!known) {
        const Result<std::vector<std::uint32_t>> samples =
            tessellate_edge_positions(state, brep, coedge.edge, deflection);
        if (!samples.has_value()) return samples.error();
        TessellateBoundaryEdge boundary;
        boundary.edge = coedge.edge;
        boundary.positions = samples.value();
        edges.push_back(std::move(boundary));
      }
      coedge_id = coedge.next;
    } while (!coedge_id.is_null() && coedge_id != loop.start);
    if (coedge_id.is_null()) {
      return Error{Errc::invalid_topology, "tessellate: loop is not closed", id.index(), 0u};
    }
    loop_id = loop.next;
  }
  return edges;
}

inline const CircleCurve3* tessellate_circle_curve(const Brep& brep, EdgeId id) {
  const Edge& edge = brep.edge(id);
  if (edge.curve == kNoGeometry || edge.curve >= brep.geometry().curves.size()) return nullptr;
  return std::get_if<CircleCurve3>(&brep.geometry().curves[edge.curve]);
}

inline void tessellate_push_triangle(TessellateState& state, FaceId id, bool reversed,
                                     std::uint32_t a, std::uint32_t b, std::uint32_t c) {
  std::array<std::uint32_t, 3> triangle{a, b, c};
  if (reversed) std::swap(triangle[1], triangle[2]);
  state.mesh.triangles.push_back(triangle);
  state.mesh.triangle_face.push_back(id.index());
}

// Cylinder lateral face: the two cap circles supply the sample count and the
// welded rings, so the caps and the lateral band share their boundary vertices.
// The lateral grid is the same edge in edge-parameter order, and the seam column
// is the seam line's own samples rather than new positions.
inline Result<void> tessellate_cylinder_face(TessellateState& state, const Brep& brep, FaceId id,
                                             const Cylinder3& cylinder, double deflection) {
  if (!tessellate_point_finite(cylinder.origin) || !tessellate_vec_finite(cylinder.axis) ||
      !std::isfinite(cylinder.radius)) {
    return Error{Errc::non_finite, "tessellate: cylinder surface is not finite", id.index(), 0u};
  }
  if (!(cylinder.radius > 0.0)) {
    return Error{Errc::degenerate, "tessellate: cylinder radius must be positive", id.index(), 0u};
  }
  const Result<std::vector<TessellateBoundaryEdge>> collected =
      tessellate_boundary_edges(state, brep, id, deflection);
  if (!collected.has_value()) return collected.error();
  const TessellateBoundaryEdge* first_circle = nullptr;
  const TessellateBoundaryEdge* second_circle = nullptr;
  for (const TessellateBoundaryEdge& boundary : collected.value()) {
    if (tessellate_circle_curve(brep, boundary.edge) == nullptr) continue;
    const Edge& edge = brep.edge(boundary.edge);
    if (!(edge.v0 == edge.v1)) continue;
    if (first_circle == nullptr) {
      first_circle = &boundary;
    } else if (second_circle == nullptr) {
      second_circle = &boundary;
    } else {
      return Error{Errc::invalid_topology, "tessellate: cylinder face has too many circles",
                   id.index(), 0u};
    }
  }
  if (first_circle == nullptr || second_circle == nullptr) {
    return Error{Errc::invalid_topology, "tessellate: cylinder face is missing a cap circle",
                 id.index(), 0u};
  }
  // The full-band grid only covers a face whose boundary is the two closed cap
  // circles plus straight seam lines. A trimmed band carries the window rim as
  // open arcs or polyline trim edges and must not be tessellated as if the
  // window were still there.
  for (const TessellateBoundaryEdge& boundary : collected.value()) {
    if (&boundary == first_circle || &boundary == second_circle) continue;
    const Edge& edge = brep.edge(boundary.edge);
    const bool line = edge.curve < brep.geometry().curves.size() &&
                      std::holds_alternative<LineCurve3>(brep.geometry().curves[edge.curve]);
    if (!line) {
      return Error{Errc::invalid_topology, "tessellate: cylinder face is not a full band",
                   id.index(), 0u};
    }
  }
  const auto axial = [&](const TessellateBoundaryEdge& boundary) {
    const Point3d& point = state.mesh.positions[boundary.positions.front()];
    return dot(point - cylinder.origin, cylinder.axis);
  };
  const TessellateBoundaryEdge* base = first_circle;
  const TessellateBoundaryEdge* top = second_circle;
  if (axial(*base) > axial(*top)) std::swap(base, top);
  const std::size_t count = base->positions.size();
  if (count < 3 || top->positions.size() != count) {
    return Error{Errc::degenerate, "tessellate: cylinder cap circles disagree", id.index(), 0u};
  }
  const bool reversed = brep.face(id).reversed;
  for (std::size_t k = 0; k < count; ++k) {
    const std::size_t next = (k + 1) % count;
    tessellate_push_triangle(state, id, reversed, base->positions[k], base->positions[next],
                             top->positions[next]);
    tessellate_push_triangle(state, id, reversed, base->positions[k], top->positions[next],
                             top->positions[k]);
  }
  return {};
}

// Cone lateral face: the base circle is used unchanged and every base sample is
// joined to the single apex vertex, which also makes the apex a fan pole rather
// than a ring copy.
inline Result<void> tessellate_cone_face(TessellateState& state, const Brep& brep, FaceId id,
                                         const Cone3& cone, double deflection) {
  if (!tessellate_point_finite(cone.origin) || !tessellate_vec_finite(cone.axis) ||
      !std::isfinite(cone.radius) || !std::isfinite(cone.half_angle)) {
    return Error{Errc::non_finite, "tessellate: cone surface is not finite", id.index(), 0u};
  }
  if (!(cone.radius > 0.0) || !(cone.half_angle > 0.0) || !(cone.half_angle < 0.5 * kPi)) {
    return Error{Errc::degenerate, "tessellate: cone surface is degenerate", id.index(), 0u};
  }
  const Result<std::vector<TessellateBoundaryEdge>> collected =
      tessellate_boundary_edges(state, brep, id, deflection);
  if (!collected.has_value()) return collected.error();
  const TessellateBoundaryEdge* base = nullptr;
  const TessellateBoundaryEdge* seam = nullptr;
  for (const TessellateBoundaryEdge& boundary : collected.value()) {
    const Edge& edge = brep.edge(boundary.edge);
    if (tessellate_circle_curve(brep, boundary.edge) != nullptr && edge.v0 == edge.v1) {
      if (base != nullptr) {
        return Error{Errc::invalid_topology, "tessellate: cone face has too many circles",
                     id.index(), 0u};
      }
      base = &boundary;
    } else {
      if (seam != nullptr) {
        return Error{Errc::invalid_topology, "tessellate: cone face has too many seam edges",
                     id.index(), 0u};
      }
      seam = &boundary;
    }
  }
  if (base == nullptr || seam == nullptr || seam->positions.size() != 2) {
    return Error{Errc::invalid_topology, "tessellate: cone face is missing its seam", id.index(),
                 0u};
  }
  if (base->positions.size() < 3) {
    return Error{Errc::degenerate, "tessellate: cone base circle is too coarse", id.index(), 0u};
  }
  const std::uint32_t first_index = seam->positions.front();
  const std::uint32_t last_index = seam->positions.back();
  const Point3d& first_point = state.mesh.positions[first_index];
  const Point3d& last_point = state.mesh.positions[last_index];
  const double first_axial = dot(first_point - cone.origin, cone.axis);
  const double last_axial = dot(last_point - cone.origin, cone.axis);
  const std::uint32_t apex = last_axial > first_axial ? last_index : first_index;
  const double dx = state.mesh.positions[apex][0] - state.mesh.positions[base->positions[0]][0];
  const double dy = state.mesh.positions[apex][1] - state.mesh.positions[base->positions[0]][1];
  const double dz = state.mesh.positions[apex][2] - state.mesh.positions[base->positions[0]][2];
  if (!(dx * dx + dy * dy + dz * dz > 0.0)) {
    return Error{Errc::degenerate, "tessellate: cone apex coincides with its base ring",
                 id.index(), 0u};
  }
  const bool reversed = brep.face(id).reversed;
  const std::size_t count = base->positions.size();
  for (std::size_t k = 0; k < count; ++k) {
    const std::size_t next = (k + 1) % count;
    tessellate_push_triangle(state, id, reversed, base->positions[k], base->positions[next], apex);
  }
  return {};
}

// Sphere face: the meridian seam supplies the v grid; every u column is a full
// turn at the same v, and both poles collapse to the seam endpoints so the mesh
// has a single south and a single north vertex.
inline Result<void> tessellate_sphere_face(TessellateState& state, const Brep& brep, FaceId id,
                                           const Sphere3& sphere, double deflection) {
  if (!tessellate_point_finite(sphere.center) || !tessellate_vec_finite(sphere.axis) ||
      !std::isfinite(sphere.radius)) {
    return Error{Errc::non_finite, "tessellate: sphere surface is not finite", id.index(), 0u};
  }
  if (!(sphere.radius > 0.0)) {
    return Error{Errc::degenerate, "tessellate: sphere radius must be positive", id.index(), 0u};
  }
  const Result<std::vector<TessellateBoundaryEdge>> collected =
      tessellate_boundary_edges(state, brep, id, deflection);
  if (!collected.has_value()) return collected.error();
  const TessellateBoundaryEdge* seam = nullptr;
  for (const TessellateBoundaryEdge& boundary : collected.value()) {
    const Edge& edge = brep.edge(boundary.edge);
    if (tessellate_circle_curve(brep, boundary.edge) != nullptr && !(edge.v0 == edge.v1)) {
      if (seam != nullptr) {
        return Error{Errc::invalid_topology, "tessellate: sphere face has too many seam arcs",
                     id.index(), 0u};
      }
      seam = &boundary;
    }
  }
  if (seam == nullptr || seam->positions.size() < 3) {
    return Error{Errc::invalid_topology, "tessellate: sphere face is missing its meridian",
                 id.index(), 0u};
  }
  const Edge& seam_edge = brep.edge(seam->edge);
  const std::size_t segments = seam->positions.size() - 1;
  const double span = seam_edge.t1 - seam_edge.t0;
  if (!std::isfinite(span) || !(span > 0.0)) {
    return Error{Errc::degenerate, "tessellate: sphere meridian interval is invalid", id.index(),
                 0u};
  }
  const double step_v = span / static_cast<double>(segments);
  double count = std::ceil(2.0 * kPi / step_v);
  if (!(count >= 3.0)) count = 3.0;
  const double maximum = static_cast<double>(std::numeric_limits<std::uint32_t>::max());
  if (!(count <= maximum)) {
    return Error{Errc::numeric_failure, "tessellate: sphere u count is too large", id.index(), 0u};
  }
  const std::uint32_t n_u = static_cast<std::uint32_t>(count);
  const double step_u = 2.0 * kPi / static_cast<double>(n_u);

  std::vector<std::vector<std::uint32_t>> rows(segments + 1);
  rows.front().push_back(seam->positions.front());
  rows.back().push_back(seam->positions.back());
  for (std::size_t j = 1; j < segments; ++j) {
    rows[j].push_back(seam->positions[j]);
    const double v = seam_edge.t0 + step_v * static_cast<double>(j);
    for (std::uint32_t k = 1; k < n_u; ++k) {
      const Result<std::uint32_t> index =
          tessellate_append_point(state, sphere.eval(step_u * static_cast<double>(k), v));
      if (!index.has_value()) return index.error();
      rows[j].push_back(index.value());
    }
  }

  const bool reversed = brep.face(id).reversed;
  for (std::size_t j = 0; j < segments; ++j) {
    const std::vector<std::uint32_t>& a = rows[j];
    const std::vector<std::uint32_t>& b = rows[j + 1];
    const bool a_pole = a.size() == 1;
    const bool b_pole = b.size() == 1;
    for (std::uint32_t k = 0; k < n_u; ++k) {
      const std::uint32_t next = (k + 1) % n_u;
      if (a_pole) {
        tessellate_push_triangle(state, id, reversed, a[0], b[next], b[k]);
      } else if (b_pole) {
        tessellate_push_triangle(state, id, reversed, a[k], a[next], b[0]);
      } else {
        tessellate_push_triangle(state, id, reversed, a[k], a[next], b[next]);
        tessellate_push_triangle(state, id, reversed, a[k], b[next], b[k]);
      }
    }
  }
  return {};
}

// Torus face: both circles are closed seams. The u seam is the ring whose
// samples are nearly axial-offset free, the v seam is the minor circle; both
// grids wrap in their own direction and reuse the seam samples at k == 0 / j == 0.
inline Result<void> tessellate_torus_face(TessellateState& state, const Brep& brep, FaceId id,
                                          const Torus3& torus, double deflection) {
  if (!tessellate_point_finite(torus.center) || !tessellate_vec_finite(torus.axis) ||
      !std::isfinite(torus.major_radius) || !std::isfinite(torus.minor_radius)) {
    return Error{Errc::non_finite, "tessellate: torus surface is not finite", id.index(), 0u};
  }
  if (!(torus.minor_radius > 0.0) || !(torus.major_radius > torus.minor_radius)) {
    return Error{Errc::degenerate, "tessellate: torus radii are degenerate", id.index(), 0u};
  }
  const Result<std::vector<TessellateBoundaryEdge>> collected =
      tessellate_boundary_edges(state, brep, id, deflection);
  if (!collected.has_value()) return collected.error();
  const TessellateBoundaryEdge* first_circle = nullptr;
  const TessellateBoundaryEdge* second_circle = nullptr;
  for (const TessellateBoundaryEdge& boundary : collected.value()) {
    if (tessellate_circle_curve(brep, boundary.edge) == nullptr) continue;
    const Edge& edge = brep.edge(boundary.edge);
    if (!(edge.v0 == edge.v1)) continue;
    if (first_circle == nullptr) {
      first_circle = &boundary;
    } else if (second_circle == nullptr) {
      second_circle = &boundary;
    } else {
      return Error{Errc::invalid_topology, "tessellate: torus face has too many circles",
                   id.index(), 0u};
    }
  }
  if (first_circle == nullptr || second_circle == nullptr) {
    return Error{Errc::invalid_topology, "tessellate: torus face is missing a seam", id.index(),
                 0u};
  }
  const auto axial_extent = [&](const TessellateBoundaryEdge& boundary) {
    double maximum = 0.0;
    for (const std::uint32_t index : boundary.positions) {
      const double offset = std::fabs(dot(state.mesh.positions[index] - torus.center, torus.axis));
      if (!(offset <= maximum)) maximum = offset;
    }
    return maximum;
  };
  const TessellateBoundaryEdge* u_seam = first_circle;
  const TessellateBoundaryEdge* v_seam = second_circle;
  if (axial_extent(*u_seam) > axial_extent(*v_seam)) std::swap(u_seam, v_seam);
  const std::size_t n_u = u_seam->positions.size();
  const std::size_t n_v = v_seam->positions.size();
  if (n_u < 3 || n_v < 3) {
    return Error{Errc::degenerate, "tessellate: torus seam circles are too coarse", id.index(), 0u};
  }
  const Edge& u_edge = brep.edge(u_seam->edge);
  const Edge& v_edge = brep.edge(v_seam->edge);
  const double step_u = (u_edge.t1 - u_edge.t0) / static_cast<double>(n_u);
  const double step_v = (v_edge.t1 - v_edge.t0) / static_cast<double>(n_v);
  if (!std::isfinite(step_u) || !std::isfinite(step_v) || !(step_u > 0.0) || !(step_v > 0.0)) {
    return Error{Errc::degenerate, "tessellate: torus seam interval is invalid", id.index(), 0u};
  }

  std::vector<std::vector<std::uint32_t>> rows(n_v);
  rows.front() = u_seam->positions;
  for (std::size_t j = 1; j < n_v; ++j) {
    rows[j].push_back(v_seam->positions[j]);
    const double v = v_edge.t0 + step_v * static_cast<double>(j);
    for (std::size_t k = 1; k < n_u; ++k) {
      const double u = u_edge.t0 + step_u * static_cast<double>(k);
      const Result<std::uint32_t> index = tessellate_append_point(state, torus.eval(u, v));
      if (!index.has_value()) return index.error();
      rows[j].push_back(index.value());
    }
  }

  const bool reversed = brep.face(id).reversed;
  for (std::size_t j = 0; j < n_v; ++j) {
    const std::vector<std::uint32_t>& a = rows[j];
    const std::vector<std::uint32_t>& b = rows[(j + 1) % n_v];
    for (std::size_t k = 0; k < n_u; ++k) {
      const std::size_t next = (k + 1) % n_u;
      tessellate_push_triangle(state, id, reversed, a[k], a[next], b[next]);
      tessellate_push_triangle(state, id, reversed, a[k], b[next], b[k]);
    }
  }
  return {};
}

// One boundary ring of a trimmed analytic face: the uv points and the welded
// position index of every point. The positions are exactly the edge samples the
// neighbouring faces see, so the trimmed triangles reuse the shared weld.
struct TessellateTrimmedRing {
  std::vector<Point2d> points;
  std::vector<std::uint32_t> positions;
  bool outer = false;
};

// Upper bound on the triangles a single tessellation may produce. The trimmed uv
// path checks it after triangulating the (already bounded) boundary rings, before
// appending anything to the mesh.
constexpr std::size_t kTrimmedMaxTriangles = 2000000u;

// Triangulates one trimmed polygon (outer ring first, then holes) and maps the
// uv triangle indices back to the shared position pool.
inline Result<std::vector<std::array<std::uint32_t, 3>>>
tessellate_trimmed_polygon_triangles(std::vector<TessellateTrimmedRing> rings, FaceId id) {
  std::vector<TrimmedBoundary> boundaries;
  TessellateTrimmedRingLayout layout;
  std::vector<std::uint32_t> positions;
  boundaries.reserve(rings.size());
  for (TessellateTrimmedRing& ring : rings) {
    TrimmedBoundary boundary;
    boundary.outer = ring.outer;
    boundary.points = std::move(ring.points);
    boundary.arcs.resize(boundary.points.size());
    layout.ring_begin.push_back(static_cast<std::uint32_t>(positions.size()));
    layout.ring_count.push_back(static_cast<std::uint32_t>(boundary.points.size()));
    positions.insert(positions.end(), ring.positions.begin(), ring.positions.end());
    boundaries.push_back(std::move(boundary));
  }
  layout.ring_of.resize(positions.size());
  for (std::uint32_t ring = 0; ring < layout.ring_begin.size(); ++ring) {
    for (std::uint32_t sample = 0; sample < layout.ring_count[ring]; ++sample) {
      layout.ring_of[layout.ring_begin[ring] + sample] = ring;
    }
  }
  const Result<TrimmedUvTriangulation> triangulation = triangulate_trimmed_domain(boundaries);
  if (!triangulation.has_value()) return triangulation.error();
  const Result<std::vector<std::array<std::uint32_t, 3>>> expanded =
      tessellate_expand_trimmed_runs(triangulation.value().points, layout,
                                     triangulation.value().triangles, id);
  if (!expanded.has_value()) return expanded.error();
  std::vector<std::array<std::uint32_t, 3>> result;
  result.reserve(expanded.value().size());
  for (const std::array<std::uint32_t, 3>& triangle : expanded.value()) {
    result.push_back({positions[triangle[0]], positions[triangle[1]], positions[triangle[2]]});
  }
  return result;
}

// A periodic trimmed ring repeats its seam vertices: the loop climbs the seam at
// one end of the unwrapped u range and descends it at the other end. Triangulating
// the cut rectangle in one piece lets the ear clipper pick diagonals from both
// copies to the same far vertex, so after welding that 3D edge has four incident
// triangles. Splitting the ring at a cut placed half a period away from the seam
// keeps every piece on one side of the weld.
inline std::optional<std::array<TessellateTrimmedRing, 2>> tessellate_split_periodic_ring(
    const TessellateTrimmedRing& ring) {
  const std::size_t count = ring.points.size();
  if (count < 6) return std::nullopt;
  std::vector<std::uint32_t> sorted = ring.positions;
  std::sort(sorted.begin(), sorted.end());
  bool wrapped = false;
  for (std::size_t i = 1; i < sorted.size(); ++i) {
    if (sorted[i] == sorted[i - 1]) {
      wrapped = true;
      break;
    }
  }
  if (!wrapped) return std::nullopt;
  double u_min = ring.points[0][0];
  double u_max = u_min;
  double v_min = ring.points[0][1];
  double v_max = v_min;
  for (const Point2d& point : ring.points) {
    u_min = std::min(u_min, point[0]);
    u_max = std::max(u_max, point[0]);
    v_min = std::min(v_min, point[1]);
    v_max = std::max(v_max, point[1]);
  }
  if (!(u_max - u_min > 0.0) || u_max - u_min > 2.0 * kPi + 1e-9) return std::nullopt;
  const double target = u_min + 0.5 * (u_max - u_min);
  const double band = 0.25 * (v_max - v_min);
  std::size_t low = count;
  std::size_t high = count;
  double low_gap = 0.0;
  double high_gap = 0.0;
  for (std::size_t i = 0; i < count; ++i) {
    const Point2d& point = ring.points[i];
    const double gap = std::abs(point[0] - target);
    if (point[1] <= v_min + band && (low == count || gap < low_gap)) {
      low = i;
      low_gap = gap;
    }
    if (point[1] >= v_max - band && (high == count || gap < high_gap)) {
      high = i;
      high_gap = gap;
    }
  }
  if (low == count || high == count || low == high) return std::nullopt;
  std::array<TessellateTrimmedRing, 2> pieces;
  pieces[0].outer = true;
  pieces[1].outer = true;
  for (std::size_t piece = 0; piece < 2; ++piece) {
    const std::size_t begin = piece == 0 ? low : high;
    const std::size_t end = piece == 0 ? high : low;
    std::size_t index = begin;
    do {
      pieces[piece].points.push_back(ring.points[index]);
      pieces[piece].positions.push_back(ring.positions[index]);
      if (index == end) break;
      index = (index + 1) % count;
    } while (index != begin);
  }
  if (pieces[0].points.size() < 3 || pieces[1].points.size() < 3) return std::nullopt;
  return pieces;
}

// Triangulates a trimmed analytic face from its boundary samples. Every loop of
// the face is projected through `projection`, unwrapped along the ring so the
// seam stays connected, normalised to outer CCW / holes CW, and handed to
// `triangulate_trimmed_domain`. The uv triangle vertices map back to the shared
// position pool, so a trim curve welds to the adjacent faces.
inline Result<void> tessellate_trimmed_face(TessellateState& state, const Brep& brep, FaceId id,
                                            const FaceUvProjection& projection,
                                            double deflection) {
  const Face& face = brep.face(id);
  if (face.first_loop.is_null()) {
    return Error{Errc::degenerate, "tessellate: face has no loop", id.index(), 0u};
  }

  std::vector<TessellateTrimmedRing> rings;
  const std::size_t loop_limit = brep.loops().size();
  const std::size_t coedge_limit = brep.coedges().size();
  LoopId loop_id = face.first_loop;
  std::size_t loop_count = 0;
  while (!loop_id.is_null()) {
    if (loop_count++ > loop_limit) {
      return Error{Errc::invalid_topology, "tessellate: loop chain is cyclic", id.index(), 0u};
    }
    const Loop& loop = brep.loop(loop_id);
    TessellateTrimmedRing ring;
    ring.outer = loop.outer;
    double previous_u = 0.0;
    bool has_previous_u = false;
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
          tessellate_edge_positions(state, brep, coedge.edge, deflection);
      if (!samples.has_value()) return samples.error();
      const std::vector<std::uint32_t>& order = samples.value();
      for (std::size_t k = 0; k < order.size(); ++k) {
        const std::uint32_t index = coedge.reversed ? order[order.size() - 1 - k] : order[k];
        if (!ring.positions.empty() && ring.positions.back() == index) continue;
        const Result<Point2d> uv =
            face_uv_project(projection, state.mesh.positions[index], previous_u, has_previous_u);
        if (!uv.has_value()) return uv.error();
        Point2d value = uv.value();
        if (projection.periodic_u() && has_previous_u) {
          value[0] = face_uv_unwrap(value[0], previous_u);
        }
        previous_u = value[0];
        has_previous_u = true;
        ring.points.push_back(value);
        ring.positions.push_back(index);
      }
      coedge_id = coedge.next;
    } while (!coedge_id.is_null() && coedge_id != loop.start);
    if (coedge_id.is_null()) {
      return Error{Errc::invalid_topology, "tessellate: loop is not closed", id.index(), 0u};
    }
    if (!ring.points.empty() && ring.positions.front() == ring.positions.back()) {
      ring.points.pop_back();
      ring.positions.pop_back();
    }
    if (ring.points.size() < 3) {
      return Error{Errc::degenerate, "tessellate: trimmed loop has fewer than three samples",
                   id.index(), 0u};
    }
    const int winding = ring_orientation2(ring.points);
    if (winding == 0) {
      return Error{Errc::degenerate, "tessellate: trimmed loop has no area", id.index(), 0u};
    }
    if ((winding > 0) != ring.outer) {
      std::reverse(ring.points.begin(), ring.points.end());
      std::reverse(ring.positions.begin(), ring.positions.end());
    }
    rings.push_back(std::move(ring));
    loop_id = loop.next;
  }

  std::size_t outer_count = 0;
  for (const TessellateTrimmedRing& ring : rings) {
    if (ring.outer) ++outer_count;
  }
  if (outer_count != 1) {
    return Error{Errc::invalid_topology, "tessellate: trimmed face needs one outer loop",
                 id.index(), static_cast<std::uint32_t>(outer_count)};
  }

  // The outer ring is triangulated first, then the holes in ring order. A
  // periodic outer ring is split away from the seam so the welded seam stays
  // manifold; rings with holes keep the single-polygon path.
  std::vector<TessellateTrimmedRing> outer;
  std::vector<TessellateTrimmedRing> holes;
  for (TessellateTrimmedRing& ring : rings) {
    if (ring.outer) {
      outer.push_back(std::move(ring));
    } else {
      holes.push_back(std::move(ring));
    }
  }

  std::vector<std::array<std::uint32_t, 3>> triangles;
  const std::optional<std::array<TessellateTrimmedRing, 2>> pieces =
      projection.periodic_u() && holes.empty()
          ? tessellate_split_periodic_ring(outer.front())
          : std::nullopt;
  if (pieces.has_value()) {
    for (const TessellateTrimmedRing& piece : pieces.value()) {
      const Result<std::vector<std::array<std::uint32_t, 3>>> part =
          tessellate_trimmed_polygon_triangles(std::vector<TessellateTrimmedRing>{piece}, id);
      if (!part.has_value()) return part.error();
      triangles.insert(triangles.end(), part.value().begin(), part.value().end());
    }
  } else {
    std::vector<TessellateTrimmedRing> polygon;
    polygon.push_back(std::move(outer.front()));
    for (TessellateTrimmedRing& hole : holes) polygon.push_back(std::move(hole));
    const Result<std::vector<std::array<std::uint32_t, 3>>> part =
        tessellate_trimmed_polygon_triangles(std::move(polygon), id);
    if (!part.has_value()) return part.error();
    triangles = std::move(part.value());
  }
  if (triangles.size() > kTrimmedMaxTriangles ||
      state.mesh.triangles.size() > kTrimmedMaxTriangles - triangles.size()) {
    return Error{Errc::numeric_failure,
                 "tessellate: trimmed face exceeds the triangle limit", id.index(),
                 static_cast<std::uint32_t>(kTrimmedMaxTriangles)};
  }

  std::size_t kept = 0;
  for (const std::array<std::uint32_t, 3>& triangle : triangles) {
    const std::uint32_t a = triangle[0];
    const std::uint32_t b = triangle[1];
    const std::uint32_t c = triangle[2];
    const Point3d& pa = state.mesh.positions[a];
    const Point3d& pb = state.mesh.positions[b];
    const Point3d& pc = state.mesh.positions[c];
    const Vec3d normal = cross(pb - pa, pc - pa);
    // A triangle that collapses in 3D (a pole ring, a seam duplicate) has no
    // normal and would be a zero-area triangle in the mesh.
    if (!(dot(normal, normal) > 0.0)) continue;
    tessellate_push_triangle(state, id, face.reversed, a, b, c);
    ++kept;
  }
  if (kept == 0) {
    return Error{Errc::degenerate, "tessellate: trimmed face has no triangle", id.index(), 0u};
  }
  return {};
}

// A face whose boundary contains a degree-1 NURBS edge is a trimmed face
// produced by the curved boolean, which rewrites every split as a polyline.
// Faces bounded only by lines and circles stay on the analytic path: a failure
// there is a shape the analytic builders do not support, and the generic uv
// triangulation must not mask it.
inline bool tessellate_face_has_polyline_edge(const Brep& brep, FaceId id) {
  const Face& face = brep.face(id);
  const std::size_t loop_limit = brep.loops().size();
  const std::size_t coedge_limit = brep.coedges().size();
  LoopId loop_id = face.first_loop;
  std::size_t loop_count = 0;
  while (!loop_id.is_null()) {
    if (loop_count++ > loop_limit) return false;
    const Loop& loop = brep.loop(loop_id);
    CoedgeId coedge_id = loop.start;
    std::size_t coedge_count = 0;
    if (coedge_id.is_null()) {
      loop_id = loop.next;
      continue;
    }
    do {
      if (coedge_count++ > coedge_limit) return false;
      const Coedge& coedge = brep.coedge(coedge_id);
      const Edge& edge = brep.edge(coedge.edge);
      if (edge.curve < brep.geometry().curves.size() &&
          std::holds_alternative<NurbsCurve>(brep.geometry().curves[edge.curve])) {
        return true;
      }
      coedge_id = coedge.next;
    } while (!coedge_id.is_null() && coedge_id != loop.start);
    loop_id = loop.next;
  }
  return false;
}

// Dispatches one face to its analytic grid. A face that the analytic builder
// cannot represent as a full primitive (a trimmed band, a cap with a window) is
// triangulated from its boundary samples in the face uv domain instead. Nurbs
// faces remain unsupported.
inline Result<void> sample_curved_face(TessellateState& state, const Brep& brep, FaceId id,
                                           const TessellationOptions& options) {
  const Face& face = brep.face(id);
  if (face.surface >= brep.geometry().surfaces.size()) {
    return Error{Errc::unsupported, "tessellate: face has no surface geometry", id.index(), 0u};
  }
  const SurfaceGeometry& surface = brep.geometry().surfaces[face.surface];
  // The analytic grid runs on a scratch copy so a rejected attempt leaves no
  // partial samples in the shared pool.
  TessellateState attempt = state;
  const auto canonical = [&](TessellateState& target) -> Result<void> {
    if (const Cylinder3* cylinder = std::get_if<Cylinder3>(&surface)) {
      return tessellate_cylinder_face(target, brep, id, *cylinder, options.deflection);
    }
    if (const Cone3* cone = std::get_if<Cone3>(&surface)) {
      return tessellate_cone_face(target, brep, id, *cone, options.deflection);
    }
    if (const Sphere3* sphere = std::get_if<Sphere3>(&surface)) {
      return tessellate_sphere_face(target, brep, id, *sphere, options.deflection);
    }
    if (const Torus3* torus = std::get_if<Torus3>(&surface)) {
      return tessellate_torus_face(target, brep, id, *torus, options.deflection);
    }
    return Error{Errc::unsupported, "tessellate: unsupported surface type", id.index(), 0u};
  };
  const Result<void> canonical_status = canonical(attempt);
  if (canonical_status.has_value()) {
    state = std::move(attempt);
    return canonical_status;
  }
  // The generic uv triangulation only serves trim faces whose boundary carries
  // polyline trim edges. Faces with pure line/circle boundaries keep the analytic
  // error so an unsupported shape is not silently remeshed.
  if (!tessellate_face_has_polyline_edge(brep, id)) return canonical_status.error();
  const Result<FaceUvProjection> projection = face_uv_projection(brep, face, options.deflection);
  if (!projection.has_value()) return canonical_status.error();
  const Result<void> trimmed =
      tessellate_trimmed_face(state, brep, id, projection.value(), options.deflection);
  if (!trimmed.has_value()) return canonical_status.error();
  return trimmed;
}

inline bool tessellate_supported_surface(const SurfaceGeometry& surface) {
  return std::holds_alternative<Plane3>(surface) ||
         std::holds_alternative<Cylinder3>(surface) || std::holds_alternative<Cone3>(surface) ||
         std::holds_alternative<Sphere3>(surface) || std::holds_alternative<Torus3>(surface);
}

}  // namespace detail

// Triangulates the faces of `brep`. Planar faces are triangulated from their
// projected rings; analytic cylinder, cone, sphere and torus faces are sampled on
// a grid whose seam and cap columns reuse the welded edge samples, so a closed
// solid has no cracks along shared edges. Every edge is sampled once in
// edge-parameter order and the samples are shared by both coedges.
//
// `angle_tolerance` does not resample the curved grids: the edge chord rule fixes
// every seam sample, and re-deriving interior columns from a normal angle would
// break the bit-identical weld to those seams. It is validated here and reserved
// for a future adaptive pass.
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
        !detail::tessellate_supported_surface(brep.geometry().surfaces[face.surface])) {
      return Error{Errc::unsupported, "tessellate: unsupported surface type", id.index(), 0u};
    }
  }

  detail::TessellateState state;
  state.edges.resize(brep.edges().size());
  state.vertices.resize(brep.vertices().size());
  state.vertex_sampled.assign(brep.vertices().size(), 0);
  for (const FaceId id : faces) {
    const Face& face = brep.face(id);
    const bool planar =
        std::holds_alternative<Plane3>(brep.geometry().surfaces[face.surface]);
    const Result<void> status = planar
                                    ? detail::tessellate_planar_face(state, brep, id, options)
                                    : detail::sample_curved_face(state, brep, id, options);
    if (!status.has_value()) return status.error();
  }
  return std::move(state.mesh);
}

}  // namespace v1
}  // namespace geom
