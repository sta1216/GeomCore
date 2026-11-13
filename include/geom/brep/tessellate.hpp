#pragma once

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <map>
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

// A degree-1 NURBS edge whose control polygon starts and ends exactly on the edge
// vertices is the polyline the curved boolean writes for a split: it can be
// sampled and welded. Degree >= 2 trim edges are a documented limitation: they are
// unsupported. Any other NURBS edge keeps the historical unsupported result, so the
// trimmed-path gate and the sampler agree on what is usable.
inline bool tessellate_polyline_edge_supported(const Brep& brep, const Edge& edge) {
  if (edge.curve == kNoGeometry || edge.curve >= brep.geometry().curves.size()) return false;
  const CurveGeometry& curve = brep.geometry().curves[edge.curve];
  const NurbsCurve* nurbs = std::get_if<NurbsCurve>(&curve);
  if (nurbs == nullptr || nurbs->degree() != 1 || !nurbs->is_valid() ||
      !std::isfinite(edge.t0) || !std::isfinite(edge.t1) || !(edge.t1 > edge.t0)) {
    return false;
  }
  return tessellate_point_close(nurbs->eval(edge.t0), brep.vertex(edge.v0).point) &&
         tessellate_point_close(nurbs->eval(edge.t1), brep.vertex(edge.v1).point);
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
    // The curved boolean stores every split as a degree-1 polyline, so sampling
    // at the uniform control parameters reproduces the control polygon exactly.
    // A curve whose endpoints disagree with the edge vertices cannot be welded
    // and keeps the historical unsupported result.
    if (!tessellate_polyline_edge_supported(brep, edge)) {
      return Error{Errc::unsupported, "tessellate: only line and circle edges are supported",
                   id.index(), 0u};
    }
    const NurbsCurve* nurbs = std::get_if<NurbsCurve>(&curve);
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
  std::vector<std::array<std::uint32_t, 3>> pending;
  pending.reserve(triangles.size());
  for (auto it = triangles.rbegin(); it != triangles.rend(); ++it) pending.push_back(*it);
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
      std::vector<std::array<std::uint32_t, 3>> fan;
      fan.reserve(run->size() + 1u);
      std::uint32_t previous = a;
      for (const std::uint32_t point : *run) {
        // The run points lie on the chord, so every fan triangle keeps the
        // winding of the triangle it replaces.
        if (sign_orient2d(points[previous][0], points[previous][1], points[point][0],
                          points[point][1], points[c][0], points[c][1]) <= 0) {
          return Error{Errc::degenerate, "tessellate: boundary fan is not convex", id.index(), 0u};
        }
        fan.push_back({previous, point, c});
        previous = point;
      }
      fan.push_back({previous, b, c});
      // The worklist is a stack, so push the fan back to front to keep the
      // original triangle order (the canonical output must not be reordered).
      for (auto it = fan.rbegin(); it != fan.rend(); ++it) pending.push_back(*it);
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

// One boundary-aligned chain of the rectilinear strip decomposition, sorted by u.
struct TessellateStripChain {
  double v = 0.0;
  std::vector<double> u;
};

// Triangulates a rectilinear trimmed uv domain by stacking boundary-aligned
// strips instead of ear clipping the whole polygon. A domain is rectilinear when
// every boundary segment is axis aligned, which is the case for an analytic face
// cut by boxes: the cap circles and the window rims are lines of constant v, the
// trim lines are lines of constant u. Ear clipping such a polygon keeps only its
// corners and spans the band with long chords, which inflates the 3D area and
// folds the mesh; the stripes instead never span more than one row, so triangles
// follow the surface and every boundary sub-segment stays present on both sides.
// Returns `false` when the ring is not rectilinear (the caller keeps the generic
// polygon path) and an error result when the plan itself is invalid.
inline Result<bool> tessellate_trimmed_strip_triangles(
    TessellateState& state, const FaceUvProjection& projection, const TessellateTrimmedRing& ring,
    FaceId id, std::vector<std::array<std::uint32_t, 3>>& triangles) {
  const std::size_t count = ring.points.size();
  if (count < 4) return false;
  constexpr double kStripEps = 1e-12;

  struct Horizontal {
    double v;
    double u0;
    double u1;
  };
  struct Vertical {
    double u;
    double v0;
    double v1;
  };
  std::vector<double> xs;
  std::vector<double> ys;
  std::vector<Horizontal> horizontals;
  std::vector<Vertical> verticals;
  xs.reserve(count);
  ys.reserve(count);
  for (std::size_t i = 0; i < count; ++i) {
    const Point2d& a = ring.points[i];
    const Point2d& b = ring.points[(i + 1) % count];
    xs.push_back(a[0]);
    ys.push_back(a[1]);
    const bool flat_v = std::fabs(a[1] - b[1]) <= kStripEps;
    const bool flat_u = std::fabs(a[0] - b[0]) <= kStripEps;
    if (!flat_v && !flat_u) return false;
    if (flat_v && !flat_u) {
      horizontals.push_back(Horizontal{a[1], std::min(a[0], b[0]), std::max(a[0], b[0])});
    } else if (flat_u && !flat_v) {
      verticals.push_back(Vertical{a[0], std::min(a[1], b[1]), std::max(a[1], b[1])});
    }
  }
  auto sort_unique = [](std::vector<double>& values) {
    std::sort(values.begin(), values.end());
    std::vector<double> unique;
    unique.reserve(values.size());
    for (const double value : values) {
      if (unique.empty() || value - unique.back() > kStripEps) unique.push_back(value);
    }
    values = std::move(unique);
  };
  sort_unique(xs);
  sort_unique(ys);
  std::vector<double> levels;
  levels.reserve(horizontals.size());
  for (const Horizontal& segment : horizontals) levels.push_back(segment.v);
  sort_unique(levels);
  if (xs.size() < 2 || levels.size() < 2) return false;

  const auto on_boundary = [&](double v, double u) {
    for (const Horizontal& segment : horizontals) {
      if (std::fabs(segment.v - v) <= kStripEps && u >= segment.u0 - kStripEps &&
          u <= segment.u1 + kStripEps) {
        return true;
      }
    }
    for (const Vertical& segment : verticals) {
      if (std::fabs(segment.u - u) <= kStripEps && v >= segment.v0 - kStripEps &&
          v <= segment.v1 + kStripEps) {
        return true;
      }
    }
    return false;
  };

  const double v_low = levels.front();
  const double v_high = levels.back();
  const auto is_cap = [&](double v) {
    return std::fabs(v - v_low) <= kStripEps || std::fabs(v - v_high) <= kStripEps;
  };

  // Builds the sorted chain of row `v` inside [u0, u1]. Cap rows (the circles the
  // face shares with its cap neighbours) stay sample-only: every vertex of that
  // edge must be one of the neighbour's samples. Interior rows also carry the
  // domain's strip columns so the bands above and below a trim row share vertices.
  const auto build_chain = [&](double v, double u0,
                               double u1) -> std::optional<TessellateStripChain> {
    TessellateStripChain chain;
    chain.v = v;
    for (std::size_t i = 0; i < count; ++i) {
      const double pu = ring.points[i][0];
      const double pv = ring.points[i][1];
      if (std::fabs(pv - v) <= kStripEps && pu >= u0 - kStripEps && pu <= u1 + kStripEps) {
        chain.u.push_back(pu);
      }
    }
    if (!is_cap(v)) {
      for (const double u : xs) {
        if (u > u0 + kStripEps && u < u1 - kStripEps && !on_boundary(v, u)) {
          chain.u.push_back(u);
        }
      }
    }
    std::sort(chain.u.begin(), chain.u.end());
    chain.u.erase(std::unique(chain.u.begin(), chain.u.end(),
                              [](double a, double b) { return b - a <= kStripEps; }),
                  chain.u.end());
    if (chain.u.size() < 2) return std::nullopt;
    if (std::fabs(chain.u.front() - u0) > 1e-9 || std::fabs(chain.u.back() - u1) > 1e-9) {
      return std::nullopt;
    }
    return chain;
  };

  struct StripPair {
    TessellateStripChain bottom;
    TessellateStripChain top;
    std::vector<std::uint32_t> left_wall;
    std::vector<std::uint32_t> right_wall;
  };
  std::vector<StripPair> plan;
  std::size_t triangle_count = 0;

  // Samples of one trim wall, ordered from the bottom row to the top row. A band
  // wall shared with a planar neighbour may be split differently on each face, so
  // every sample of the wall is a vertex of this face as well.
  const auto wall_samples = [&](double u, double va, double vb) {
    std::vector<std::pair<double, std::uint32_t>> ordered;
    for (std::size_t i = 0; i < count; ++i) {
      const Point2d& point = ring.points[i];
      if (std::fabs(point[0] - u) <= kStripEps && point[1] >= va - kStripEps &&
          point[1] <= vb + kStripEps) {
        ordered.emplace_back(point[1], ring.positions[i]);
      }
    }
    std::sort(ordered.begin(), ordered.end(),
              [](const auto& a, const auto& b) { return a.first < b.first; });
    std::vector<std::uint32_t> samples;
    for (const auto& entry : ordered) {
      if (samples.empty() || samples.back() != entry.second) samples.push_back(entry.second);
    }
    return samples;
  };

  // A band joins two rows and is closed with a zipper instead of a grid: the rows
  // carry their own samples and their chains may disagree. The trim walls keep
  // their own samples too; the emission fans the wall edge of the corner triangle
  // so no wall sample is lost.
  const auto add_zipper = [&](double va, double vb, double u0, double u1) -> bool {
    const std::optional<TessellateStripChain> bottom = build_chain(va, u0, u1);
    const std::optional<TessellateStripChain> top = build_chain(vb, u0, u1);
    if (!bottom.has_value() || !top.has_value()) return false;
    triangle_count += bottom.value().u.size() + top.value().u.size() - 2u;
    std::vector<std::uint32_t> left_wall = wall_samples(u0, va, vb);
    std::vector<std::uint32_t> right_wall = wall_samples(u1, va, vb);
    if (left_wall.size() > 2u) triangle_count += left_wall.size() - 2u;
    if (right_wall.size() > 2u) triangle_count += right_wall.size() - 2u;
    plan.push_back(StripPair{bottom.value(), top.value(), std::move(left_wall),
                             std::move(right_wall)});
    return true;
  };

  for (std::size_t band = 0; band + 1 < levels.size(); ++band) {
    const double va = levels[band];
    const double vb = levels[band + 1];
    const double vm = 0.5 * (va + vb);
    std::size_t run_begin = xs.size();
    for (std::size_t i = 0; i + 1 < xs.size(); ++i) {
      const double um = 0.5 * (xs[i] + xs[i + 1]);
      std::size_t crossings = 0;
      for (const Horizontal& segment : horizontals) {
        if (segment.v < vm && segment.u0 < um && um < segment.u1) ++crossings;
      }
      const bool inside = (crossings % 2u) == 1u;
      if (inside) {
        if (run_begin == xs.size()) run_begin = i;
        continue;
      }
      if (run_begin != xs.size()) {
        const double u0 = xs[run_begin];
        const double u1 = xs[i];
        const bool ok = add_zipper(va, vb, u0, u1);
        if (!ok) return false;
        run_begin = xs.size();
      }
    }
    if (run_begin != xs.size()) {
      const double u0 = xs[run_begin];
      const double u1 = xs.back();
      const bool ok = add_zipper(va, vb, u0, u1);
      if (!ok) return false;
    }
  }
  if (plan.empty()) return false;
  if (triangle_count > kTrimmedMaxTriangles ||
      state.mesh.triangles.size() > kTrimmedMaxTriangles - triangle_count) {
    return Error{Errc::numeric_failure, "tessellate: trimmed face exceeds the triangle limit",
                 id.index(), static_cast<std::uint32_t>(kTrimmedMaxTriangles)};
  }

  std::map<std::pair<double, double>, std::uint32_t> vertices;
  for (std::size_t i = 0; i < count; ++i) {
    vertices.insert_or_assign(std::make_pair(ring.points[i][1], ring.points[i][0]),
                              ring.positions[i]);
  }
  const auto vertex_at = [&](double v, double u) -> Result<std::uint32_t> {
    const auto found = vertices.find(std::make_pair(v, u));
    if (found != vertices.end()) return found->second;
    const Point2d uv{u, v};
    const Point3d point = face_uv_eval(projection, uv);
    const Result<std::uint32_t> index = tessellate_append_point(state, point);
    if (!index.has_value()) return index.error();
    vertices.emplace(std::make_pair(v, u), index.value());
    return index.value();
  };

  triangles.clear();
  triangles.reserve(triangle_count);
  // Replaces the corner triangle of a band whose wall carries more samples than
  // its two corner vertices by a fan over those wall samples: the wall is a shared
  // boundary edge, so every sample of it must be a vertex of this face too.
  const auto fan_wall = [&](std::size_t first,
                            const std::vector<std::uint32_t>& wall) -> Result<void> {
    if (wall.size() <= 2) return {};
    for (std::size_t t = first; t < triangles.size(); ++t) {
      std::array<std::uint32_t, 3>& tri = triangles[t];
      const bool forward =
          (tri[0] == wall.front() && tri[1] == wall.back()) ||
          (tri[1] == wall.front() && tri[2] == wall.back()) ||
          (tri[2] == wall.front() && tri[0] == wall.back());
      const bool backward =
          (tri[0] == wall.back() && tri[1] == wall.front()) ||
          (tri[1] == wall.back() && tri[2] == wall.front()) ||
          (tri[2] == wall.back() && tri[0] == wall.front());
      if (!forward && !backward) continue;
      std::uint32_t corner = tri[0];
      for (int k = 0; k < 3; ++k) {
        if (tri[k] != wall.front() && tri[k] != wall.back()) corner = tri[k];
      }
      tri = forward ? std::array<std::uint32_t, 3>{corner, wall[0], wall[1]}
                    : std::array<std::uint32_t, 3>{corner, wall[1], wall[0]};
      for (std::size_t i = 1; i + 1 < wall.size(); ++i) {
        triangles.push_back(forward ? std::array<std::uint32_t, 3>{corner, wall[i], wall[i + 1]}
                                    : std::array<std::uint32_t, 3>{corner, wall[i + 1], wall[i]});
      }
      return {};
    }
    return Error{Errc::degenerate, "tessellate: trimmed band wall is not closed", id.index(), 0u};
  };
  for (const StripPair& pair : plan) {
    const std::size_t nb = pair.bottom.u.size();
    const std::size_t nt = pair.top.u.size();
    std::vector<std::uint32_t> bottom(nb);
    std::vector<std::uint32_t> top(nt);
    for (std::size_t i = 0; i < nb; ++i) {
      const Result<std::uint32_t> index = vertex_at(pair.bottom.v, pair.bottom.u[i]);
      if (!index.has_value()) return index.error();
      bottom[i] = index.value();
    }
    for (std::size_t i = 0; i < nt; ++i) {
      const Result<std::uint32_t> index = vertex_at(pair.top.v, pair.top.u[i]);
      if (!index.has_value()) return index.error();
      top[i] = index.value();
    }
    const std::size_t first = triangles.size();
    std::size_t i = 0;
    std::size_t j = 0;
    while (i + 1 < nb || j + 1 < nt) {
      if (j + 1 >= nt || (i + 1 < nb && pair.bottom.u[i + 1] <= pair.top.u[j + 1])) {
        triangles.push_back({bottom[i], bottom[i + 1], top[j]});
        ++i;
      } else {
        triangles.push_back({top[j], bottom[i], top[j + 1]});
        ++j;
      }
    }
    const Result<void> left = fan_wall(first, pair.left_wall);
    if (!left.has_value()) return left.error();
    const Result<void> right = fan_wall(first, pair.right_wall);
    if (!right.has_value()) return right.error();
  }
  return true;
}

// Splits interior edges of a trimmed-face triangulation that are longer than the
// chord bound implied by `deflection`. Boundary edges (pairs of ring samples) are
// never split: they carry the vertices a neighbouring face sampled and must weld
// bit-identically. The refinement only adds interior vertices, so the mesh stays
// watertight while every triangle follows the surface closely.
inline Result<void> tessellate_refine_trimmed_edges(
    TessellateState& state, const FaceUvProjection& projection,
    std::vector<std::array<std::uint32_t, 3>>& triangles,
    std::vector<std::uint64_t> constrained, double deflection, FaceId id) {
  // An edge is refined only when the surface point above its chord midpoint is
  // farther than this from the chord: that is the exact deviation bound the
  // trimmed mesh promises, so refining to a length bound would over-split the
  // thin fan triangles that the boundary expansion creates.
  const double sag_limit = 0.9 * deflection;
  if (!(sag_limit > 0.0) || triangles.empty()) return {};
  std::sort(constrained.begin(), constrained.end());
  const auto key = [](std::uint32_t a, std::uint32_t b) -> std::uint64_t {
    if (a > b) std::swap(a, b);
    return (static_cast<std::uint64_t>(a) << 32) | static_cast<std::uint64_t>(b);
  };
  const auto edge_length = [&](std::uint32_t a, std::uint32_t b) {
    const Point3d& pa = state.mesh.positions[a];
    const Point3d& pb = state.mesh.positions[b];
    double square = 0.0;
    for (int k = 0; k < 3; ++k) {
      const double d = pa[k] - pb[k];
      square += d * d;
    }
    return std::sqrt(square);
  };
  const auto is_boundary = [&](std::uint64_t edge) {
    return std::binary_search(constrained.begin(), constrained.end(), edge);
  };
  std::vector<bool> alive(triangles.size(), true);
  std::map<std::uint64_t, std::vector<std::uint32_t>> adjacency;
  for (std::size_t i = 0; i < triangles.size(); ++i) {
    for (int k = 0; k < 3; ++k) {
      adjacency[key(triangles[i][k], triangles[i][(k + 1) % 3])].push_back(
          static_cast<std::uint32_t>(i));
    }
  }
  std::size_t splits = 0;
  const std::size_t split_budget = 8 * triangles.size() + 64;
  std::vector<std::uint64_t> queue;
  for (const auto& entry : adjacency) {
    if (is_boundary(entry.first)) continue;
    queue.push_back(entry.first);
  }
  while (!queue.empty()) {
    // Refinement is bounded: a domain whose base triangulation is chordal cannot
    // be healed by splitting, so stop after a fixed budget of splits. The caller
    // then rejects the mesh through its deflection check instead of running to
    // the global triangle limit.
    if (splits >= split_budget) break;
    const std::uint64_t edge = queue.back();
    queue.pop_back();
    const auto found = adjacency.find(edge);
    if (found == adjacency.end()) continue;
    if (is_boundary(edge)) continue;
    const std::uint32_t a = static_cast<std::uint32_t>(edge >> 32);
    const std::uint32_t b = static_cast<std::uint32_t>(edge & 0xffffffffu);
    // The split point is the chord midpoint projected onto the surface. The
    // midpoint parameter is estimated from the two endpoint parameters so a
    // periodic seam cannot send the projection to the far side of the face;
    // the estimate is clamped to the segment's parameter box.
    const Point3d& pa = state.mesh.positions[a];
    const Point3d& pb = state.mesh.positions[b];
    const Result<Point2d> first = face_uv_project(projection, pa, 0.0, false);
    if (!first.has_value()) continue;
    const Point2d ua = first.value();
    const Result<Point2d> second = face_uv_project(projection, pb, ua[0], true);
    if (!second.has_value()) continue;
    Point2d ub = second.value();
    if (projection.periodic_u()) ub[0] = face_uv_unwrap(ub[0], ua[0]);
    if (projection.periodic_v()) ub[1] = face_uv_unwrap(ub[1], ua[1]);
    Point2d middle_uv{0.5 * (ua[0] + ub[0]), 0.5 * (ua[1] + ub[1])};
    middle_uv[0] = std::min(std::max(middle_uv[0], std::min(ua[0], ub[0])),
                            std::max(ua[0], ub[0]));
    middle_uv[1] = std::min(std::max(middle_uv[1], std::min(ua[1], ub[1])),
                            std::max(ua[1], ub[1]));
    const Point3d chord_middle{0.5 * (pa[0] + pb[0]), 0.5 * (pa[1] + pb[1]),
                               0.5 * (pa[2] + pb[2])};
    const Point3d surface_middle = face_uv_eval(projection, middle_uv);
    double sag_square = 0.0;
    for (int k = 0; k < 3; ++k) {
      const double d = chord_middle[k] - surface_middle[k];
      sag_square += d * d;
    }
    if (sag_square <= sag_limit * sag_limit) continue;
    const Result<std::uint32_t> appended = tessellate_append_point(state, surface_middle);
    if (!appended.has_value()) return appended.error();
    const std::uint32_t mid = appended.value();
    if (!(edge_length(a, mid) < 0.999 * edge_length(a, b)) ||
        !(edge_length(mid, b) < 0.999 * edge_length(a, b))) {
      continue;
    }
    ++splits;
    const std::vector<std::uint32_t> adjacent = found->second;
    adjacency.erase(found);
    for (const std::uint32_t t : adjacent) {
      if (!alive[t]) continue;
      alive[t] = false;
      const std::array<std::uint32_t, 3> triangle = triangles[t];
      std::uint32_t c = triangle[0];
      for (int k = 0; k < 3; ++k) {
        if (triangle[k] != a && triangle[k] != b) c = triangle[k];
      }
      for (int k = 0; k < 3; ++k) {
        const std::uint64_t other = key(triangle[k], triangle[(k + 1) % 3]);
        if (other == edge) continue;
        const auto it = adjacency.find(other);
        if (it == adjacency.end()) continue;
        std::vector<std::uint32_t>& list = it->second;
        list.erase(std::remove(list.begin(), list.end(), t), list.end());
        if (list.empty()) adjacency.erase(it);
      }
      // Keep the winding: the split edge is traversed a->b or b->a in the triangle.
      bool forward = false;
      for (int k = 0; k < 3; ++k) {
        if (triangle[k] == a && triangle[(k + 1) % 3] == b) forward = true;
      }
      const std::array<std::array<std::uint32_t, 3>, 2> added =
          forward ? std::array<std::array<std::uint32_t, 3>, 2>{{{a, mid, c}, {mid, b, c}}}
                  : std::array<std::array<std::uint32_t, 3>, 2>{{{b, mid, c}, {mid, a, c}}};
      for (const std::array<std::uint32_t, 3>& piece : added) {
        const std::uint32_t index = static_cast<std::uint32_t>(triangles.size());
        triangles.push_back(piece);
        alive.push_back(true);
        for (int k = 0; k < 3; ++k) {
          const std::uint64_t other = key(piece[k], piece[(k + 1) % 3]);
          adjacency[other].push_back(index);
          if (!is_boundary(other)) queue.push_back(other);
        }
      }
      if (triangles.size() > kTrimmedMaxTriangles) {
        return Error{Errc::numeric_failure, "tessellate: trimmed face exceeds the triangle limit",
                     id.index(), static_cast<std::uint32_t>(kTrimmedMaxTriangles)};
      }
    }
  }
  std::vector<std::array<std::uint32_t, 3>> compact;
  compact.reserve(triangles.size());
  for (std::size_t i = 0; i < triangles.size(); ++i) {
    if (alive[i]) compact.push_back(triangles[i]);
  }
  triangles = std::move(compact);
  return {};
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

  // A cheap upper bound on the triangles the generic path can emit is checked
  // before that path allocates its triangulation: a simple polygon of n vertices
  // produces at most n - 2 triangles, plus at most two per bridged hole, and the
  // fan expansion never adds more triangles than the ring has samples.
  std::size_t triangle_bound = 2u * rings.size();
  for (const TessellateTrimmedRing& ring : rings) {
    triangle_bound += ring.points.size();
  }
  if (triangle_bound > kTrimmedMaxTriangles ||
      state.mesh.triangles.size() > kTrimmedMaxTriangles - triangle_bound) {
    return Error{Errc::numeric_failure,
                 "tessellate: trimmed face exceeds the triangle limit", id.index(),
                 static_cast<std::uint32_t>(kTrimmedMaxTriangles)};
  }

  std::vector<std::array<std::uint32_t, 3>> triangles;
  bool stripped = false;
  if (holes.empty()) {
    const Result<bool> strip = tessellate_trimmed_strip_triangles(
        state, projection, outer.front(), id, triangles);
    if (!strip.has_value()) return strip.error();
    stripped = strip.value();
  }
  if (!stripped) {
    triangles.clear();
    const std::optional<std::array<TessellateTrimmedRing, 2>> pieces =
        projection.periodic_u()
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
  }
  // The sample vertices are fixed by the neighbouring faces; only interior edges
  // of the triangulation may be split, so the trimmed mesh follows the surface
  // without ever splitting a shared boundary edge.
  std::vector<std::uint64_t> constrained;
  const auto add_ring_uvs = [&](const TessellateTrimmedRing& ring) {
    const std::size_t ring_size = ring.positions.size();
    constrained.reserve(constrained.size() + ring_size);
    for (std::size_t i = 0; i < ring_size; ++i) {
      std::uint32_t a = ring.positions[i];
      std::uint32_t b = ring.positions[(i + 1) % ring_size];
      if (a > b) std::swap(a, b);
      constrained.push_back((static_cast<std::uint64_t>(a) << 32) | static_cast<std::uint64_t>(b));
    }
  };
  for (const TessellateTrimmedRing& ring : outer) add_ring_uvs(ring);
  for (const TessellateTrimmedRing& ring : holes) add_ring_uvs(ring);
  const Result<void> refined = tessellate_refine_trimmed_edges(
      state, projection, triangles, std::move(constrained), deflection, id);
  if (!refined.has_value()) return refined.error();
  if (triangles.size() > kTrimmedMaxTriangles ||
      state.mesh.triangles.size() > kTrimmedMaxTriangles - triangles.size()) {
    return Error{Errc::numeric_failure,
                 "tessellate: trimmed face exceeds the triangle limit", id.index(),
                 static_cast<std::uint32_t>(kTrimmedMaxTriangles)};
  }
  // Every kept triangle must follow the surface: the surface point above the
  // triangle centroid has to stay within the requested deflection of the
  // centroid. A chordal triangulation of a domain that the strip decomposition
  // cannot cover is rejected here instead of silently breaking the contract.
  for (const std::array<std::uint32_t, 3>& triangle : triangles) {
    const Point3d& p0 = state.mesh.positions[triangle[0]];
    const Point3d& p1 = state.mesh.positions[triangle[1]];
    const Point3d& p2 = state.mesh.positions[triangle[2]];
    const Result<Point2d> uv0 = face_uv_project(projection, p0, 0.0, false);
    if (!uv0.has_value()) continue;
    const Result<Point2d> uv1 = face_uv_project(projection, p1, uv0.value()[0], true);
    if (!uv1.has_value()) continue;
    const Result<Point2d> uv2 = face_uv_project(projection, p2, uv0.value()[0], true);
    if (!uv2.has_value()) continue;
    Point2d second = uv1.value();
    Point2d third = uv2.value();
    if (projection.periodic_u()) {
      second[0] = face_uv_unwrap(second[0], uv0.value()[0]);
      third[0] = face_uv_unwrap(third[0], uv0.value()[0]);
    }
    if (projection.periodic_v()) {
      second[1] = face_uv_unwrap(second[1], uv0.value()[1]);
      third[1] = face_uv_unwrap(third[1], uv0.value()[1]);
    }
    const Point2d centroid_uv{(uv0.value()[0] + second[0] + third[0]) / 3.0,
                              (uv0.value()[1] + second[1] + third[1]) / 3.0};
    const Point3d surface = face_uv_eval(projection, centroid_uv);
    const Point3d centroid{(p0[0] + p1[0] + p2[0]) / 3.0, (p0[1] + p1[1] + p2[1]) / 3.0,
                           (p0[2] + p1[2] + p2[2]) / 3.0};
    const double dx = surface[0] - centroid[0];
    const double dy = surface[1] - centroid[1];
    const double dz = surface[2] - centroid[2];
    if (std::sqrt(dx * dx + dy * dy + dz * dz) > deflection) {
      return Error{Errc::unsupported,
                   "tessellate: trimmed face does not follow the surface", id.index(), 0u};
    }
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
      if (tessellate_polyline_edge_supported(brep, edge)) return true;
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
  // error so an unsupported shape is not silently remeshed. A polyline-bounded
  // face that the trimmed path also rejects surfaces the trimmed-path error: the
  // analytic error would describe a primitive the face is not.
  if (!tessellate_face_has_polyline_edge(brep, id)) return canonical_status.error();
  const Result<FaceUvProjection> projection = face_uv_projection(brep, face, options.deflection);
  if (!projection.has_value()) return projection.error();
  const Result<void> trimmed =
      tessellate_trimmed_face(state, brep, id, projection.value(), options.deflection);
  if (!trimmed.has_value()) return trimmed.error();
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
// A trimmed analytic face (a curved-boolean face whose boundary carries degree-1
// polyline trim edges) first tries its analytic builder; when that fails it is
// triangulated from its boundary uv rings and refined in uv until the chord
// deviation is under `deflection`.
//
// `angle_tolerance` does not resample the curved grids or the trimmed bands: the
// edge chord rule fixes every boundary sample, and re-deriving interior columns
// from a normal angle would break the bit-identical weld to those seams. It is
// validated here and reserved for a future adaptive pass.
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



