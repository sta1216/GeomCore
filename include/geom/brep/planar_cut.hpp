#pragma once

// Cuts one planar B-rep face into 2D fragments along 3D segment cuts.
//
// The face boundary and the cuts are projected into the face frame. A cut
// endpoint within `tol` of a boundary vertex snaps to that vertex; otherwise an
// endpoint within `tol` of a boundary edge is projected onto the edge and a new
// boundary vertex is inserted there. Cut endpoints within `tol` of each other
// are welded to one point. The resulting loops (source 0) and cuts (source
// 1..n) are handed to Arrangement2.
//
// Every bounded arrangement face owns exactly one boundary cycle, so a nested
// cycle is a separate face. Fragments are related by a containment tree: a
// fragment's outer ring is its own cycle (CCW) and its holes are the cycles of
// its immediate children (CW). A fragment's sample is the centroid of a
// triangle of the fragment; it is kept only when the corresponding 3D point is
// strictly inside the original face. Fragments are returned by descending area,
// then by the lexicographic first outer-ring point, so the result is stable.
//
// Errors:
//   - invalid_topology for a missing face, a face without a usable boundary,
//     or a cut that dangles (an open endpoint strictly inside the face).
//   - non_finite for a non-finite tolerance or cut coordinate.
//   - degenerate for a negative tolerance, a cut that projects to a point, or a
//     face with no interior sample.
//   - numeric_failure when split_segments cannot represent an intersection.
//   - unsupported for a face that is not planar (propagated from face_query).

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <optional>
#include <span>
#include <vector>

#include <geom/brep/face_intersect.hpp>
#include <geom/brep/face_query.hpp>
#include <geom/brep/plane_frame.hpp>
#include <geom/brep/topology.hpp>
#include <geom/cgeom/arrangement2.hpp>
#include <geom/cgeom/polygon2.hpp>
#include <geom/cgeom/triangulate2.hpp>
#include <geom/core/error.hpp>
#include <geom/core/result.hpp>
#include <geom/primitives/point.hpp>
#include <geom/primitives/vec.hpp>
#include <geom/surfaces/plane.hpp>

namespace geom {
inline namespace v1 {

// One fragment of a planar face: an outer ring plus zero or more hole rings,
// all in the face frame. The rings are closed implicitly (the first point is
// not repeated). `sample` is a point strictly inside the fragment.
struct FaceFragment2d {
  std::vector<Point2d> outer;
  std::vector<std::vector<Point2d>> holes;
  Point2d sample{};
};

namespace detail {

inline bool planar_cut_finite2(const Point2d& p) noexcept {
  return std::isfinite(p[0]) && std::isfinite(p[1]);
}

inline bool planar_cut_finite3(const Point3d& p) noexcept {
  return std::isfinite(p[0]) && std::isfinite(p[1]) && std::isfinite(p[2]);
}

inline bool planar_cut_lex_less(const Point2d& a, const Point2d& b) noexcept {
  if (a[0] != b[0]) return a[0] < b[0];
  return a[1] < b[1];
}

inline bool planar_cut_same_xy(const Point2d& a, const Point2d& b) noexcept {
  return a[0] == b[0] && a[1] == b[1];
}

inline Point2d planar_cut_project_on_edge(const Point2d& p, const Point2d& a,
                                          const Point2d& b) noexcept {
  const Vec2d edge = b - a;
  const double length2 = dot(edge, edge);
  if (!(length2 > 0.0)) return a;
  double t = dot(p - a, edge) / length2;
  if (!(t > 0.0)) {
    t = 0.0;
  } else if (t >= 1.0) {
    t = 1.0;
  }
  return Point2d{a[0] + edge[0] * t, a[1] + edge[1] * t};
}

// Reverses the ring when needed and rotates it so the lexicographically
// smallest vertex comes first. The rotation makes the ring start independent of
// the arrangement's internal vertex ids.
inline void planar_cut_normalize_ring(std::vector<Point2d>& ring, bool want_ccw) {
  if (ring.size() < 3) return;
  const int orientation = ring_orientation2(ring);
  if (orientation != 0 && ((orientation > 0) != want_ccw)) {
    std::reverse(ring.begin(), ring.end());
  }
  std::size_t first = 0;
  for (std::size_t i = 1; i < ring.size(); ++i) {
    if (planar_cut_lex_less(ring[i], ring[first])) first = i;
  }
  if (first != 0) {
    std::rotate(ring.begin(), ring.begin() + static_cast<std::ptrdiff_t>(first),
                ring.end());
  }
}

inline double planar_cut_fragment_area(const FaceFragment2d& fragment) noexcept {
  double area = std::fabs(signed_area2(fragment.outer));
  for (const std::vector<Point2d>& hole : fragment.holes) {
    area -= std::fabs(signed_area2(hole));
  }
  return area;
}

// Where one cut endpoint landed relative to the projected face boundary.
struct PlanarCutSnap {
  Point2d point{};
  bool boundary = false;
  bool vertex = false;
  std::size_t ring = 0;
  std::size_t edge = 0;
};

inline PlanarCutSnap planar_cut_find_snap(
    const std::vector<std::vector<Point2d>>& loops, const Point2d& p, double tol) {
  PlanarCutSnap snap;
  snap.point = p;

  // A vertex within tol wins over an edge, even when the edge is closer: two
  // cuts ending at the same corner must land on the same exact coordinate.
  bool found_vertex = false;
  double best_vertex = 0.0;
  for (std::size_t r = 0; r < loops.size(); ++r) {
    for (const Point2d& vertex : loops[r]) {
      const double distance = norm2(vertex - p);
      if (!(distance <= tol)) continue;
      if (!found_vertex || distance < best_vertex) {
        found_vertex = true;
        best_vertex = distance;
        snap.point = vertex;
        snap.boundary = true;
        snap.vertex = true;
        snap.ring = r;
        snap.edge = 0;
      }
    }
  }
  if (found_vertex) return snap;

  bool found_edge = false;
  double best_edge = 0.0;
  for (std::size_t r = 0; r < loops.size(); ++r) {
    const std::vector<Point2d>& ring = loops[r];
    for (std::size_t i = 0; i < ring.size(); ++i) {
      const Point2d& a = ring[i];
      const Point2d& b = ring[(i + 1) % ring.size()];
      const double distance = face_query_point_segment_distance(p, a, b);
      if (!(distance <= tol)) continue;
      if (!found_edge || distance < best_edge) {
        found_edge = true;
        best_edge = distance;
        snap.point = planar_cut_project_on_edge(p, a, b);
        snap.boundary = true;
        snap.vertex = false;
        snap.ring = r;
        snap.edge = i;
      }
    }
  }
  return snap;
}

// Prefers a boundary snap, then a vertex snap, then the lexicographically
// smallest coordinate. This picks one exact weld point per snapped group.
inline bool planar_cut_canonical_better(const PlanarCutSnap& candidate,
                                        const PlanarCutSnap& current) noexcept {
  if (candidate.boundary != current.boundary) return candidate.boundary;
  if (candidate.vertex != current.vertex) return candidate.vertex;
  return planar_cut_lex_less(candidate.point, current.point);
}

inline std::size_t planar_cut_find_root(std::vector<std::size_t>& parent,
                                        std::size_t index) noexcept {
  while (parent[index] != index) {
    parent[index] = parent[parent[index]];
    index = parent[index];
  }
  return index;
}

// One new boundary vertex, inserted into the ring edge it snapped to.
struct PlanarCutInsert {
  std::size_t ring = 0;
  std::size_t edge = 0;
  Point2d point{};
};

// Inserts every snap point into its ring edge. The ring is rebuilt once, so all
// edge indices still refer to the original ring.
inline void planar_cut_apply_inserts(std::vector<std::vector<Point2d>>& rings,
                                     std::vector<PlanarCutInsert> inserts) {
  if (inserts.empty()) return;
  std::sort(inserts.begin(), inserts.end(),
            [](const PlanarCutInsert& lhs, const PlanarCutInsert& rhs) noexcept {
              if (lhs.ring != rhs.ring) return lhs.ring < rhs.ring;
              if (lhs.edge != rhs.edge) return lhs.edge < rhs.edge;
              return planar_cut_lex_less(lhs.point, rhs.point);
            });
  inserts.erase(std::unique(inserts.begin(), inserts.end(),
                            [](const PlanarCutInsert& lhs,
                               const PlanarCutInsert& rhs) noexcept {
                              return lhs.ring == rhs.ring && lhs.edge == rhs.edge &&
                                     planar_cut_same_xy(lhs.point, rhs.point);
                            }),
                inserts.end());

  std::size_t index = 0;
  while (index < inserts.size()) {
    const std::size_t ring_index = inserts[index].ring;
    std::size_t group_end = index;
    while (group_end < inserts.size() && inserts[group_end].ring == ring_index) {
      ++group_end;
    }
    std::vector<Point2d>& ring = rings[ring_index];
    const std::size_t vertex_count = ring.size();
    std::vector<std::vector<std::pair<double, Point2d>>> extra(vertex_count);
    for (std::size_t i = index; i < group_end; ++i) {
      const std::size_t edge = inserts[i].edge;
      const Point2d& a = ring[edge];
      const Point2d& b = ring[(edge + 1) % vertex_count];
      const Vec2d edge_vector = b - a;
      const double length2 = dot(edge_vector, edge_vector);
      const double t =
          length2 > 0.0 ? dot(inserts[i].point - a, edge_vector) / length2 : 0.0;
      extra[edge].emplace_back(t, inserts[i].point);
    }
    std::size_t extra_total = 0;
    for (std::vector<std::pair<double, Point2d>>& edge_extra : extra) {
      std::sort(edge_extra.begin(), edge_extra.end(),
                [](const std::pair<double, Point2d>& lhs,
                   const std::pair<double, Point2d>& rhs) noexcept {
                  if (lhs.first != rhs.first) return lhs.first < rhs.first;
                  return planar_cut_lex_less(lhs.second, rhs.second);
                });
      extra_total += edge_extra.size();
    }
    std::vector<Point2d> rebuilt;
    rebuilt.reserve(vertex_count + extra_total);
    for (std::size_t i = 0; i < vertex_count; ++i) {
      rebuilt.push_back(ring[i]);
      for (const std::pair<double, Point2d>& entry : extra[i]) {
        rebuilt.push_back(entry.second);
      }
    }
    ring.swap(rebuilt);
    index = group_end;
  }
}

// Centroid of the first triangle of the fragment (trying every triangle), or of
// a triangle of the outer ring alone when the full polygon fails to
// triangulate. The returned point is strictly inside the fragment: it is inside
// the outer ring and outside every hole.
inline std::optional<Point2d> planar_cut_fragment_sample(
    const std::vector<Point2d>& outer,
    const std::vector<std::vector<Point2d>>& holes) {
  const auto inside = [&outer, &holes](const Point2d& p) noexcept {
    if (!face_query_point_in_ring(outer, p)) return false;
    for (const std::vector<Point2d>& hole : holes) {
      if (face_query_point_in_ring(hole, p)) return false;
    }
    return true;
  };
  const auto centroid = [](const Point2d& a, const Point2d& b,
                           const Point2d& c) noexcept {
    return Point2d{(a[0] + b[0] + c[0]) / 3.0, (a[1] + b[1] + c[1]) / 3.0};
  };

  std::vector<std::vector<Point2d>> rings;
  rings.reserve(1 + holes.size());
  rings.push_back(outer);
  rings.insert(rings.end(), holes.begin(), holes.end());
  const Polygon2d polygon = detail::face_intersect_polygon(rings);
  const Triangulation2 triangulation = triangulate_polygon2(polygon);
  if (triangulation.ok) {
    for (const std::array<std::uint32_t, 3>& triangle : triangulation.triangles) {
      const Point2d candidate =
          centroid(polygon.vertices[triangle[0]], polygon.vertices[triangle[1]],
                   polygon.vertices[triangle[2]]);
      if (inside(candidate)) return candidate;
    }
  }

  // A ring with a keyhole slit may fail the full triangulation. Triangulate the
  // outer ring alone and still reject candidates that fall inside a hole.
  Polygon2d outer_polygon;
  outer_polygon.vertices = outer;
  outer_polygon.ring_starts.push_back(0u);
  outer_polygon.ring_starts.push_back(
      static_cast<std::uint32_t>(outer.size()));
  const Triangulation2 outer_triangulation = triangulate_polygon2(outer_polygon);
  if (outer_triangulation.ok) {
    for (const std::array<std::uint32_t, 3>& triangle :
         outer_triangulation.triangles) {
      const Point2d candidate = centroid(outer_polygon.vertices[triangle[0]],
                                         outer_polygon.vertices[triangle[1]],
                                         outer_polygon.vertices[triangle[2]]);
      if (inside(candidate)) return candidate;
    }
  }
  return std::nullopt;
}

// Normalizes the rings and computes their sample. Returns nullopt when no
// strictly interior sample can be found.
inline std::optional<FaceFragment2d> planar_cut_make_fragment(
    const std::vector<Point2d>& outer,
    const std::vector<std::vector<Point2d>>& holes) {
  FaceFragment2d fragment;
  fragment.outer = outer;
  planar_cut_normalize_ring(fragment.outer, true);
  fragment.holes.reserve(holes.size());
  for (const std::vector<Point2d>& hole : holes) {
    std::vector<Point2d> normalized = hole;
    planar_cut_normalize_ring(normalized, false);
    fragment.holes.push_back(std::move(normalized));
  }
  std::optional<Point2d> sample =
      planar_cut_fragment_sample(fragment.outer, fragment.holes);
  if (!sample.has_value()) return std::nullopt;
  fragment.sample = sample.value();
  return fragment;
}

struct PlanarCutFace {
  std::vector<Point2d> ring;
  double area = 0.0;
  Point2d rep{};
  bool has_rep = false;
  std::size_t parent = static_cast<std::size_t>(-1);
  std::vector<std::size_t> children;
};

// Walks the single boundary cycle owned by an arrangement face.
inline bool planar_cut_walk_cycle(const Arrangement2& arrangement,
                                  std::uint32_t face,
                                  std::vector<Point2d>& ring) {
  const std::span<const Arrangement2::HalfEdge> half_edges =
      arrangement.half_edges();
  const std::span<const Point2d> vertices = arrangement.vertices();
  const std::uint32_t start = arrangement.faces()[face].edge;
  std::uint32_t half_edge = start;
  std::size_t guard = 0;
  do {
    ring.push_back(vertices[half_edges[half_edge].origin]);
    half_edge = half_edges[half_edge].next;
    if (++guard > half_edges.size()) return false;
  } while (half_edge != start);
  return true;
}

}  // namespace detail

// Splits `face` into fragments along `segments`; see the file comment for the
// snapping, containment and error rules.
inline Result<std::vector<FaceFragment2d>> cut_face_by_segments(
    const Brep& brep, FaceId face, std::span<const Segment3> segments,
    double tol) {
  if (face.is_null() || !brep.face_alive(face)) {
    return Error{Errc::invalid_topology, "planar_cut: face is not alive",
                 face.index(), 0u};
  }
  if (!std::isfinite(tol)) {
    return Error{Errc::non_finite, "planar_cut: tolerance is not finite",
                 face.index(), 0u};
  }
  if (tol < 0.0) {
    return Error{Errc::degenerate, "planar_cut: tolerance is negative",
                 face.index(), 0u};
  }

  const Result<Plane3> plane = face_plane(brep, face);
  if (!plane) return plane.error();
  const Result<PlaneFrame> frame_result = make_plane_frame(plane.value());
  if (!frame_result) return frame_result.error();
  const PlaneFrame& frame = frame_result.value();

  const Result<std::vector<std::vector<Point2d>>> loops_result =
      face_loops_2d(brep, face, frame);
  if (!loops_result) return loops_result.error();
  const std::vector<std::vector<Point2d>>& loops = loops_result.value();
  if (loops.empty()) {
    return Error{Errc::invalid_topology, "planar_cut: face has no boundary loops",
                 face.index(), 0u};
  }
  for (const std::vector<Point2d>& ring : loops) {
    for (const Point2d& point : ring) {
      if (!detail::planar_cut_finite2(point)) {
        return Error{Errc::non_finite,
                     "planar_cut: boundary point is not finite", face.index(),
                     0u};
      }
    }
  }

  // No cuts: the whole face is one fragment, holes included.
  if (segments.empty()) {
    std::vector<std::vector<Point2d>> holes(loops.begin() + 1, loops.end());
    std::optional<FaceFragment2d> fragment =
        detail::planar_cut_make_fragment(loops.front(), holes);
    if (!fragment.has_value()) {
      return Error{Errc::degenerate, "planar_cut: face has no interior sample",
                   face.index(), 0u};
    }
    std::vector<FaceFragment2d> fragments;
    fragments.push_back(std::move(fragment.value()));
    return fragments;
  }

  struct Endpoint {
    Point2d projected{};
    detail::PlanarCutSnap snap{};
    Point2d final{};
  };

  std::vector<Endpoint> endpoints;
  endpoints.reserve(segments.size() * 2);
  for (const Segment3& segment : segments) {
    if (!detail::planar_cut_finite3(segment.a) ||
        !detail::planar_cut_finite3(segment.b)) {
      return Error{Errc::non_finite,
                   "planar_cut: cut point is not finite", face.index(), 0u};
    }
    const Point2d a = plane_frame_to_2d(frame, segment.a);
    const Point2d b = plane_frame_to_2d(frame, segment.b);
    if (!detail::planar_cut_finite2(a) || !detail::planar_cut_finite2(b)) {
      return Error{Errc::non_finite,
                   "planar_cut: projected cut point is not finite",
                   face.index(), 0u};
    }
    if (detail::planar_cut_same_xy(a, b)) {
      return Error{Errc::degenerate, "planar_cut: cut projects to a point",
                   face.index(), 0u};
    }
    endpoints.push_back(Endpoint{a, detail::planar_cut_find_snap(loops, a, tol)});
    endpoints.push_back(Endpoint{b, detail::planar_cut_find_snap(loops, b, tol)});
  }

  // Weld cut endpoints that coincide within tol, then give every member of a
  // welded group the same exact coordinate.
  const std::size_t endpoint_count = endpoints.size();
  std::vector<std::size_t> parent(endpoint_count);
  for (std::size_t i = 0; i < endpoint_count; ++i) parent[i] = i;
  for (std::size_t i = 0; i < endpoint_count; ++i) {
    for (std::size_t j = i + 1; j < endpoint_count; ++j) {
      if (norm2(endpoints[i].projected - endpoints[j].projected) <= tol) {
        const std::size_t root_i = detail::planar_cut_find_root(parent, i);
        const std::size_t root_j = detail::planar_cut_find_root(parent, j);
        if (root_i != root_j) parent[root_j] = root_i;
      }
    }
  }
  std::vector<detail::PlanarCutSnap> canonical(endpoint_count);
  std::vector<char> has_canonical(endpoint_count, 0);
  for (std::size_t i = 0; i < endpoint_count; ++i) {
    const std::size_t root = detail::planar_cut_find_root(parent, i);
    if (has_canonical[root] == 0 ||
        detail::planar_cut_canonical_better(endpoints[i].snap, canonical[root])) {
      canonical[root] = endpoints[i].snap;
      has_canonical[root] = 1;
    }
  }
  for (std::size_t i = 0; i < endpoint_count; ++i) {
    endpoints[i].final =
        canonical[detail::planar_cut_find_root(parent, i)].point;
  }

  // Insert the welded edge snaps into the boundary rings. Vertex snaps already
  // name an existing ring point and need no insertion.
  std::vector<std::vector<Point2d>> rings = loops;
  std::vector<detail::PlanarCutInsert> inserts;
  std::vector<char> inserted(endpoint_count, 0);
  for (std::size_t i = 0; i < endpoint_count; ++i) {
    const std::size_t root = detail::planar_cut_find_root(parent, i);
    if (inserted[root] != 0) continue;
    inserted[root] = 1;
    const detail::PlanarCutSnap& snap = canonical[root];
    if (!snap.boundary || snap.vertex) continue;
    bool already_present = false;
    for (const Point2d& vertex : rings[snap.ring]) {
      if (detail::planar_cut_same_xy(vertex, snap.point)) {
        already_present = true;
        break;
      }
    }
    if (!already_present) {
      inserts.push_back({snap.ring, snap.edge, snap.point});
    }
  }
  detail::planar_cut_apply_inserts(rings, std::move(inserts));

  // Boundary rings are source 0; each cut is its own source so the arrangement
  // keeps them distinct.
  std::vector<ArrangementSegment> arrangement_segments;
  std::size_t segment_total = 0;
  for (const std::vector<Point2d>& ring : rings) segment_total += ring.size();
  arrangement_segments.reserve(segment_total + segments.size());
  for (const std::vector<Point2d>& ring : rings) {
    for (std::size_t i = 0; i < ring.size(); ++i) {
      const Point2d& a = ring[i];
      const Point2d& b = ring[(i + 1) % ring.size()];
      if (detail::planar_cut_same_xy(a, b)) continue;
      arrangement_segments.push_back(ArrangementSegment{a, b, 0u});
    }
  }
  for (std::size_t i = 0; i < segments.size(); ++i) {
    const Point2d& a = endpoints[i * 2].final;
    const Point2d& b = endpoints[i * 2 + 1].final;
    if (detail::planar_cut_same_xy(a, b)) {
      return Error{Errc::degenerate, "planar_cut: cut collapses to a point",
                   face.index(), static_cast<std::uint32_t>(i)};
    }
    arrangement_segments.push_back(
        ArrangementSegment{a, b, static_cast<std::uint32_t>(i + 1)});
  }

  const detail::SplitResult split = detail::split_segments(arrangement_segments);
  if (!std::isfinite(split.max_residual)) {
    return Error{Errc::numeric_failure,
                 "planar_cut: segment intersection is not representable",
                 face.index(), 0u};
  }
  const Arrangement2 arrangement = Arrangement2::build(arrangement_segments);

  // A cut endpoint strictly inside the face leaves a degree-1 vertex: the cut
  // does not close, and using it would only produce a zero-width ring spike.
  // Endpoints outside the face (a cut that misses, or one that overshoots) are
  // harmless and are filtered out later by the sample test.
  std::vector<std::uint32_t> degree(arrangement.vertex_count(), 0u);
  for (const Arrangement2::HalfEdge& half_edge : arrangement.half_edges()) {
    ++degree[half_edge.origin];
  }
  for (std::uint32_t vertex = 0; vertex < arrangement.vertex_count(); ++vertex) {
    if (degree[vertex] != 1u) continue;
    const Point3d point =
        plane_frame_to_3d(frame, arrangement.vertices()[vertex]);
    const Result<bool> inside = point_in_face(brep, face, point, tol);
    if (!inside) return inside.error();
    if (inside.value()) {
      return Error{Errc::invalid_topology,
                   "planar_cut: a cut dangles inside the face", face.index(),
                   vertex};
    }
  }

  // Collect every bounded arrangement face. A nested cycle is its own face.
  std::vector<detail::PlanarCutFace> faces;
  for (std::uint32_t index = 0; index < arrangement.face_count(); ++index) {
    if (!arrangement.faces()[index].bounded) continue;
    detail::PlanarCutFace record;
    if (!detail::planar_cut_walk_cycle(arrangement, index, record.ring)) {
      return Error{Errc::invalid_topology,
                   "planar_cut: arrangement face cycle is not closed",
                   face.index(), index};
    }
    record.area = std::fabs(signed_area2(record.ring));
    if (!(record.area > 0.0)) continue;

    // A representative point strictly inside the (hole-free) cycle decides
    // containment. The first triangle's centroid is used because the cycle is
    // a simple ring.
    Polygon2d polygon;
    polygon.vertices = record.ring;
    polygon.ring_starts.push_back(0u);
    polygon.ring_starts.push_back(
        static_cast<std::uint32_t>(record.ring.size()));
    const Triangulation2 triangulation = triangulate_polygon2(polygon);
    if (triangulation.ok && !triangulation.triangles.empty()) {
      const std::array<std::uint32_t, 3>& triangle =
          triangulation.triangles.front();
      record.rep = Point2d{
          (polygon.vertices[triangle[0]][0] + polygon.vertices[triangle[1]][0] +
           polygon.vertices[triangle[2]][0]) /
              3.0,
          (polygon.vertices[triangle[0]][1] + polygon.vertices[triangle[1]][1] +
           polygon.vertices[triangle[2]][1]) /
              3.0};
      record.has_rep = true;
    }
    faces.push_back(std::move(record));
  }

  // Parent of a face is the smallest-area larger face whose cycle contains its
  // representative point. Immediate children become the fragment's holes.
  for (std::size_t index = 0; index < faces.size(); ++index) {
    if (!faces[index].has_rep) continue;
    std::size_t best = static_cast<std::size_t>(-1);
    for (std::size_t candidate = 0; candidate < faces.size(); ++candidate) {
      if (candidate == index) continue;
      if (!(faces[candidate].area > faces[index].area)) continue;
      if (!detail::face_query_point_in_ring(faces[candidate].ring,
                                            faces[index].rep)) {
        continue;
      }
      if (best == static_cast<std::size_t>(-1) ||
          faces[candidate].area < faces[best].area) {
        best = candidate;
      }
    }
    faces[index].parent = best;
  }
  for (std::size_t index = 0; index < faces.size(); ++index) {
    if (faces[index].parent != static_cast<std::size_t>(-1)) {
      faces[faces[index].parent].children.push_back(index);
    }
  }

  std::vector<FaceFragment2d> fragments;
  for (const detail::PlanarCutFace& record : faces) {
    if (!record.has_rep) continue;
    std::vector<std::vector<Point2d>> holes;
    holes.reserve(record.children.size());
    for (const std::size_t child : record.children) {
      holes.push_back(faces[child].ring);
    }
    std::optional<FaceFragment2d> fragment =
        detail::planar_cut_make_fragment(record.ring, holes);
    if (!fragment.has_value()) continue;
    const Point3d sample3 = plane_frame_to_3d(frame, fragment.value().sample);
    const Result<bool> inside = point_in_face(brep, face, sample3, tol);
    if (!inside) return inside.error();
    if (!inside.value()) continue;
    fragments.push_back(std::move(fragment.value()));
  }

  std::sort(fragments.begin(), fragments.end(),
            [](const FaceFragment2d& lhs, const FaceFragment2d& rhs) noexcept {
              const double area_lhs = detail::planar_cut_fragment_area(lhs);
              const double area_rhs = detail::planar_cut_fragment_area(rhs);
              if (area_lhs != area_rhs) return area_lhs > area_rhs;
              return detail::planar_cut_lex_less(lhs.outer.front(),
                                                 rhs.outer.front());
            });
  return fragments;
}

}  // namespace v1
}  // namespace geom
