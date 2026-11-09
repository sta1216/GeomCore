#pragma once

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <limits>
#include <map>
#include <numbers>
#include <numeric>
#include <span>
#include <unordered_map>
#include <utility>
#include <vector>

#include <geom/brep/plane_frame.hpp>
#include <geom/brep/solid_editor.hpp>
#include <geom/cgeom/polygon2.hpp>
#include <geom/cgeom/triangulate2.hpp>
#include <geom/core/error.hpp>
#include <geom/core/result.hpp>

namespace geom {
inline namespace v1 {

// One planar fragment of a solid boundary. `outer` and every ring in `holes`
// are open rings: the last point is not repeated. `normal` is the outward
// normal of the fragment and is used to orient the face; an outer ring must
// wind counter-clockwise about `normal`, and hole rings must wind clockwise.
// Every ring point must lie within `tol` of the plane through the first outer
// point with normal `normal`; a non-planar ring is rejected with
// `Errc::tolerance_violation`. The assembler welds coincident endpoints, builds
// the B-rep topology, groups the fragments into closed components, and returns
// one solid per outward-oriented component.
//
// Closed-solid convention: a component whose faces use every undirected edge
// exactly twice, once in each direction, encloses a positive signed volume and
// is an outer shell. A closed component whose normals point inward encloses a
// negative signed volume and is a cavity; it must be contained by exactly one
// outer component. The current topology validator defers multi-shell solids,
// so an outer component and its cavities are placed in one closed shell.
struct AssembledFace {
  std::vector<Point3d> outer;
  std::vector<std::vector<Point3d>> holes;
  Vec3d normal;
};

namespace detail {

inline bool assemble_point_finite(Point3d p) noexcept {
  return std::isfinite(p[0]) && std::isfinite(p[1]) && std::isfinite(p[2]);
}

inline bool assemble_vec_finite(Vec3d v) noexcept {
  return std::isfinite(v[0]) && std::isfinite(v[1]) && std::isfinite(v[2]);
}

// Unit vector for a finite, nonzero input. The largest component is scaled out
// first so the intermediate length stays in [1, sqrt(3)] and cannot overflow.
inline Vec3d assemble_unit(Vec3d v) noexcept {
  const double scale =
      std::max({std::fabs(v[0]), std::fabs(v[1]), std::fabs(v[2])});
  if (!(scale > 0.0) || !std::isfinite(scale)) return Vec3d{};
  const double x = v[0] / scale;
  const double y = v[1] / scale;
  const double z = v[2] / scale;
  const double length = std::sqrt(x * x + y * y + z * z);
  return Vec3d{x / length, y / length, z / length};
}

struct assemble_cell {
  std::int64_t x = 0;
  std::int64_t y = 0;
  std::int64_t z = 0;

  bool operator==(const assemble_cell& other) const noexcept {
    return x == other.x && y == other.y && z == other.z;
  }
};

struct assemble_cell_hash {
  std::size_t operator()(const assemble_cell& cell) const noexcept {
    const std::size_t x = std::hash<std::int64_t>{}(cell.x);
    const std::size_t y = std::hash<std::int64_t>{}(cell.y);
    const std::size_t z = std::hash<std::int64_t>{}(cell.z);
    return x ^ (y << 1) ^ (z << 2);
  }
};

// Welds points within `tol` of an earlier point. Coordinates are binned with
// floor(x / tol), so any two points at distance <= tol land in the same cell or
// in adjacent cells, and the 27 neighbouring cells are searched. The cell
// coordinate is clamped to +/-9e18 before the conversion: finite inputs can
// still make x / tol overflow to infinity (a tiny tol and a huge coordinate),
// and converting an out-of-range floating value to long long is undefined
// behaviour. Clamping keeps the conversion defined; the exact-distance test
// below still decides whether the point welds, so clamped points are not merged
// unless they are actually within `tol`.
class assemble_welder {
 public:
  explicit assemble_welder(double tol) noexcept : tol_(tol) {}

  Result<std::uint32_t> weld(Point3d point) {
    if (!assemble_point_finite(point)) {
      return Error{Errc::non_finite, "assemble_solid: vertex is not finite"};
    }
    const assemble_cell cell = cell_of(point);
    for (std::int64_t dx = -1; dx <= 1; ++dx) {
      for (std::int64_t dy = -1; dy <= 1; ++dy) {
        for (std::int64_t dz = -1; dz <= 1; ++dz) {
          const auto found =
              cells_.find(assemble_cell{cell.x + dx, cell.y + dy, cell.z + dz});
          if (found == cells_.end()) continue;
          for (const std::uint32_t index : found->second) {
            const Point3d other = points_[index];
            const double distance = std::hypot(point[0] - other[0],
                                               point[1] - other[1],
                                               point[2] - other[2]);
            if (distance <= tol_) return index;
          }
        }
      }
    }
    const std::uint32_t index = static_cast<std::uint32_t>(points_.size());
    points_.push_back(point);
    cells_[cell].push_back(index);
    return index;
  }

  const std::vector<Point3d>& points() const noexcept { return points_; }

  double tolerance() const noexcept { return tol_; }

 private:
  assemble_cell cell_of(Point3d point) const noexcept {
    constexpr double kCellLimit = 9.0e18;
    const auto snap = [this](double coordinate) noexcept -> std::int64_t {
      double scaled = coordinate / tol_;
      if (scaled > kCellLimit) scaled = kCellLimit;
      if (scaled < -kCellLimit) scaled = -kCellLimit;
      return static_cast<std::int64_t>(std::floor(scaled));
    };
    return assemble_cell{snap(point[0]), snap(point[1]), snap(point[2])};
  }

  double tol_;
  std::vector<Point3d> points_;
  std::unordered_map<assemble_cell, std::vector<std::uint32_t>,
                     assemble_cell_hash>
      cells_;
};

inline Result<std::vector<std::uint32_t>> assemble_weld_ring(
    const std::vector<Point3d>& ring, assemble_welder& welder) {
  std::vector<std::uint32_t> welded;
  welded.reserve(ring.size());
  for (const Point3d& point : ring) {
    const Result<std::uint32_t> index = welder.weld(point);
    if (!index) return index.error();
    if (welded.empty() || welded.back() != index.value()) {
      welded.push_back(index.value());
    }
  }
  while (welded.size() > 1u && welded.front() == welded.back()) {
    welded.pop_back();
  }
  if (welded.size() < 3u) {
    return Error{Errc::degenerate,
                 "assemble_solid: ring has fewer than three distinct points"};
  }
  return welded;
}

struct assemble_face_data {
  Vec3d normal{};
  double offset = 0.0;
  PlaneFrame frame{};
  std::vector<std::vector<std::uint32_t>> rings;
  std::vector<std::array<std::uint32_t, 3>> triangles;
};

inline Result<assemble_face_data> assemble_make_face(const AssembledFace& face,
                                                     assemble_welder& welder) {
  if (!assemble_vec_finite(face.normal)) {
    return Error{Errc::non_finite,
                 "assemble_solid: face normal is not finite"};
  }
  const double normal_length =
      std::hypot(face.normal[0], face.normal[1], face.normal[2]);
  if (!(normal_length > 0.0)) {
    return Error{Errc::degenerate, "assemble_solid: face normal is zero"};
  }
  const Vec3d normal = assemble_unit(face.normal);

  std::vector<std::vector<std::uint32_t>> rings;
  Result<std::vector<std::uint32_t>> outer =
      assemble_weld_ring(face.outer, welder);
  if (!outer) return outer.error();
  rings.push_back(std::move(outer.value()));
  for (const std::vector<Point3d>& hole : face.holes) {
    Result<std::vector<std::uint32_t>> welded =
        assemble_weld_ring(hole, welder);
    if (!welded) return welded.error();
    rings.push_back(std::move(welded.value()));
  }

  const Point3d anchor = welder.points()[rings.front().front()];
  const double offset =
      normal[0] * anchor[0] + normal[1] * anchor[1] +
      normal[2] * anchor[2];
  if (!std::isfinite(offset)) {
    return Error{Errc::non_finite,
                 "assemble_solid: face offset is not finite"};
  }
  // Every ring must lie on the plane selected by the supplied normal. Without
  // this check, a non-planar input is silently projected onto the plane and can
  // produce a plausible but wrong volume.
  const double tolerance = welder.tolerance();
  for (const std::vector<std::uint32_t>& ring : rings) {
    for (const std::uint32_t index : ring) {
      const Point3d point = welder.points()[index];
      const double distance =
          normal[0] * (point[0] - anchor[0]) +
          normal[1] * (point[1] - anchor[1]) +
          normal[2] * (point[2] - anchor[2]);
      if (!std::isfinite(distance) || std::fabs(distance) > tolerance) {
        return Error{
            Errc::tolerance_violation,
            "assemble_solid: face ring is not planar within tolerance"};
      }
    }
  }
  const Result<PlaneFrame> frame_result =
      make_plane_frame(Plane3{normal, offset});
  if (!frame_result) return frame_result.error();
  const PlaneFrame frame = frame_result.value();

  // Make the loop winding agree with the supplied outward normal: outer rings
  // are counter-clockwise and holes clockwise in the face frame. This keeps
  // the signed-volume and mass-property integration consistent even when a
  // caller supplies a ring in the opposite order.
  for (std::size_t ring_index = 0; ring_index < rings.size(); ++ring_index) {
    std::vector<Point2d> projected;
    projected.reserve(rings[ring_index].size());
    for (const std::uint32_t index : rings[ring_index]) {
      projected.push_back(plane_frame_to_2d(frame, welder.points()[index]));
    }
    const double area = signed_area2(std::span<const Point2d>(projected));
    if (!std::isfinite(area) || area == 0.0) {
      return Error{Errc::degenerate, "assemble_solid: ring has zero area"};
    }
    const bool want_positive = ring_index == 0u;
    if ((area > 0.0) != want_positive) {
      std::reverse(rings[ring_index].begin(), rings[ring_index].end());
    }
  }

  Polygon2d polygon;
  polygon.ring_starts.push_back(0u);
  std::vector<std::uint32_t> flat_indices;
  for (const std::vector<std::uint32_t>& ring : rings) {
    for (const std::uint32_t index : ring) {
      polygon.vertices.push_back(
          plane_frame_to_2d(frame, welder.points()[index]));
      flat_indices.push_back(index);
    }
    polygon.ring_starts.push_back(
        static_cast<std::uint32_t>(flat_indices.size()));
  }
  if (validate_polygon2(polygon) != Errc::none) {
    return Error{Errc::degenerate, "assemble_solid: face polygon is invalid"};
  }
  const double face_area = polygon_area2(polygon);
  if (!std::isfinite(face_area) || !(face_area > 0.0)) {
    return Error{Errc::degenerate, "assemble_solid: face has zero area"};
  }
  const Triangulation2 triangulation = triangulate_polygon2(polygon);
  if (!triangulation.ok || triangulation.triangles.empty()) {
    return Error{Errc::degenerate,
                 "assemble_solid: face triangulation failed"};
  }

  assemble_face_data data;
  data.normal = normal;
  data.offset = offset;
  data.frame = frame;
  data.rings = std::move(rings);
  data.triangles.reserve(triangulation.triangles.size());
  for (const std::array<std::uint32_t, 3>& triangle :
       triangulation.triangles) {
    data.triangles.push_back(
        {flat_indices[triangle[0]], flat_indices[triangle[1]],
         flat_indices[triangle[2]]});
  }
  return data;
}

struct assemble_edge_use {
  std::uint32_t face = 0;
  bool forward = false;
};

struct assemble_edge_record {
  std::uint32_t a = 0;
  std::uint32_t b = 0;
  std::vector<assemble_edge_use> uses;
};

class assemble_disjoint_set {
 public:
  explicit assemble_disjoint_set(std::size_t size)
      : parent_(size), rank_(size, 0u) {
    std::iota(parent_.begin(), parent_.end(), 0u);
  }

  std::uint32_t find(std::uint32_t value) {
    std::uint32_t root = value;
    while (parent_[root] != root) root = parent_[root];
    while (parent_[value] != value) {
      const std::uint32_t next = parent_[value];
      parent_[value] = root;
      value = next;
    }
    return root;
  }

  void unite(std::uint32_t first, std::uint32_t second) {
    first = find(first);
    second = find(second);
    if (first == second) return;
    if (rank_[first] < rank_[second]) std::swap(first, second);
    parent_[second] = first;
    if (rank_[first] == rank_[second]) ++rank_[first];
  }

 private:
  std::vector<std::uint32_t> parent_;
  std::vector<std::uint32_t> rank_;
};

struct assemble_component {
  std::vector<std::uint32_t> faces;
  double volume = 0.0;
  Point3d sample{};
};

inline Result<double> assemble_component_volume(
    const std::vector<std::uint32_t>& component_faces,
    const std::vector<assemble_face_data>& face_data,
    const std::vector<Point3d>& points) {
  const Point3d reference =
      points[face_data[component_faces.front()].rings.front().front()];
  double volume = 0.0;
  for (const std::uint32_t face : component_faces) {
    for (const std::array<std::uint32_t, 3>& triangle :
         face_data[face].triangles) {
      const Vec3d a = points[triangle[0]] - reference;
      const Vec3d b = points[triangle[1]] - reference;
      const Vec3d c = points[triangle[2]] - reference;
      volume += dot(a, cross(b, c)) / 6.0;
    }
  }
  if (!std::isfinite(volume)) {
    return Error{Errc::numeric_failure,
                 "assemble_solid: component volume is not finite"};
  }
  return volume;
}

inline Point3d assemble_component_sample(
    const std::vector<std::uint32_t>& component_faces,
    const std::vector<assemble_face_data>& face_data,
    const std::vector<Point3d>& points) {
  std::vector<std::uint32_t> indices;
  for (const std::uint32_t face : component_faces) {
    for (const std::vector<std::uint32_t>& ring : face_data[face].rings) {
      indices.insert(indices.end(), ring.begin(), ring.end());
    }
  }
  std::sort(indices.begin(), indices.end());
  indices.erase(std::unique(indices.begin(), indices.end()), indices.end());
  Point3d sample{};
  if (indices.empty()) return sample;
  for (const std::uint32_t index : indices) {
    sample[0] += points[index][0];
    sample[1] += points[index][1];
    sample[2] += points[index][2];
  }
  const double count = static_cast<double>(indices.size());
  sample[0] /= count;
  sample[1] /= count;
  sample[2] /= count;
  return sample;
}

// Winding number of the closed triangle set around `sample`, using the signed
// solid angle of each triangle. A full turn (4*pi) means the point is inside.
inline bool assemble_point_inside(
    Point3d sample, const std::vector<std::uint32_t>& component_faces,
    const std::vector<assemble_face_data>& face_data,
    const std::vector<Point3d>& points) {
  double angle = 0.0;
  for (const std::uint32_t face : component_faces) {
    for (const std::array<std::uint32_t, 3>& triangle :
         face_data[face].triangles) {
      const Vec3d a = points[triangle[0]] - sample;
      const Vec3d b = points[triangle[1]] - sample;
      const Vec3d c = points[triangle[2]] - sample;
      const double length_a = std::hypot(a[0], a[1], a[2]);
      const double length_b = std::hypot(b[0], b[1], b[2]);
      const double length_c = std::hypot(c[0], c[1], c[2]);
      const double numerator = dot(a, cross(b, c));
      const double denominator =
          length_a * length_b * length_c + dot(a, b) * length_c +
          dot(b, c) * length_a + dot(c, a) * length_b;
      angle += 2.0 * std::atan2(numerator, denominator);
    }
  }
  return std::isfinite(angle) && std::fabs(angle) > 2.0 * std::numbers::pi;
}

inline bool assemble_same_plane(const Vec3d& first, double first_offset,
                                const Vec3d& second, double second_offset,
                                double tol) noexcept {
  return std::fabs(first[0] - second[0]) <= tol &&
         std::fabs(first[1] - second[1]) <= tol &&
         std::fabs(first[2] - second[2]) <= tol &&
         std::fabs(first_offset - second_offset) <= tol;
}

}  // namespace detail

inline Result<Brep> assemble_solid(std::span<const AssembledFace> faces,
                                   double tol) {
  if (!std::isfinite(tol) || !(tol > 0.0)) {
    return Error{Errc::tolerance_violation,
                 "assemble_solid: tolerance must be finite and positive"};
  }
  if (faces.empty()) {
    return Error{Errc::invalid_topology, "assemble_solid: no faces"};
  }

  detail::assemble_welder welder(tol);
  std::vector<detail::assemble_face_data> face_data;
  face_data.reserve(faces.size());
  for (const AssembledFace& face : faces) {
    Result<detail::assemble_face_data> made =
        detail::assemble_make_face(face, welder);
    if (!made) return made.error();
    face_data.push_back(std::move(made.value()));
  }

  std::map<std::pair<std::uint32_t, std::uint32_t>, std::uint32_t> edge_lookup;
  std::vector<detail::assemble_edge_record> edge_records;
  for (std::uint32_t face = 0; face < face_data.size(); ++face) {
    for (const std::vector<std::uint32_t>& ring : face_data[face].rings) {
      for (std::size_t i = 0; i < ring.size(); ++i) {
        const std::uint32_t from = ring[i];
        const std::uint32_t to = ring[(i + 1u) % ring.size()];
        if (from == to) {
          return Error{Errc::degenerate,
                       "assemble_solid: ring has a zero-length edge"};
        }
        const std::uint32_t a = std::min(from, to);
        const std::uint32_t b = std::max(from, to);
        const auto key = std::make_pair(a, b);
        auto found = edge_lookup.find(key);
        std::uint32_t edge_index = 0;
        if (found == edge_lookup.end()) {
          edge_index = static_cast<std::uint32_t>(edge_records.size());
          edge_lookup.emplace(key, edge_index);
          edge_records.push_back(detail::assemble_edge_record{a, b, {}});
        } else {
          edge_index = found->second;
        }
        edge_records[edge_index].uses.push_back(
            detail::assemble_edge_use{face, from == a});
      }
    }
  }

  // Report non-manifold edges before boundary edges so a three-face fan is
  // diagnosed as non-manifold rather than as an unrelated open boundary.
  for (std::size_t i = 0; i < edge_records.size(); ++i) {
    const detail::assemble_edge_record& edge = edge_records[i];
    if (edge.uses.size() > 2u ||
        (edge.uses.size() == 2u &&
         edge.uses[0].forward == edge.uses[1].forward)) {
      return Error{Errc::invalid_topology,
                   "assemble_solid: non-manifold edge",
                   static_cast<std::uint32_t>(i), 0u};
    }
  }
  for (std::size_t i = 0; i < edge_records.size(); ++i) {
    if (edge_records[i].uses.size() == 1u) {
      return Error{Errc::not_closed, "assemble_solid: edge is used once",
                   static_cast<std::uint32_t>(i), 0u};
    }
  }

  SolidEditor editor;
  std::vector<VertexId> vertex_ids(welder.points().size());
  for (std::size_t i = 0; i < welder.points().size(); ++i) {
    const Result<VertexId> vertex = editor.add_vertex(welder.points()[i], tol);
    if (!vertex) return vertex.error();
    vertex_ids[i] = vertex.value();
  }

  std::vector<EdgeId> edge_ids(edge_records.size());
  for (std::size_t i = 0; i < edge_records.size(); ++i) {
    const detail::assemble_edge_record& edge = edge_records[i];
    const Point3d first = welder.points()[edge.a];
    const Point3d second = welder.points()[edge.b];
    const Vec3d delta = second - first;
    const double scale =
        std::max({std::fabs(delta[0]), std::fabs(delta[1]),
                  std::fabs(delta[2])});
    if (!(scale > 0.0) || !std::isfinite(scale)) {
      return Error{Errc::degenerate, "assemble_solid: edge has zero length"};
    }
    const Vec3d scaled{delta[0] / scale, delta[1] / scale, delta[2] / scale};
    const double unit_length =
        std::sqrt(dot(scaled, scaled));
    const double length = scale * unit_length;
    if (!(length > 0.0) || !std::isfinite(length)) {
      return Error{Errc::degenerate, "assemble_solid: edge has zero length"};
    }
    const Vec3d direction{scaled[0] / unit_length, scaled[1] / unit_length,
                          scaled[2] / unit_length};
    const std::uint32_t curve =
        editor.add_curve(CurveGeometry{LineCurve3{first, direction}});
    const Result<EdgeId> made =
        editor.add_edge(vertex_ids[edge.a], vertex_ids[edge.b], curve, 0.0,
                        length, tol);
    if (!made) return made.error();
    edge_ids[i] = made.value();
  }

  std::vector<std::pair<Vec3d, double>> planes;
  std::vector<FaceId> face_ids(face_data.size());
  for (std::size_t face = 0; face < face_data.size(); ++face) {
    const detail::assemble_face_data& data = face_data[face];
    std::uint32_t surface = kNoGeometry;
    for (std::size_t plane = 0; plane < planes.size(); ++plane) {
      if (detail::assemble_same_plane(planes[plane].first,
                                      planes[plane].second, data.normal,
                                      data.offset, tol)) {
        surface = static_cast<std::uint32_t>(plane);
        break;
      }
    }
    if (surface == kNoGeometry) {
      surface = editor.add_surface(SurfaceGeometry{Plane3{data.normal,
                                                          data.offset}});
      planes.push_back({data.normal, data.offset});
    }
    const Result<FaceId> made = editor.add_face(surface, false, tol);
    if (!made) return made.error();
    face_ids[face] = made.value();
  }

  for (std::size_t face = 0; face < face_data.size(); ++face) {
    const detail::assemble_face_data& data = face_data[face];
    for (std::size_t ring_index = 0; ring_index < data.rings.size();
         ++ring_index) {
      const std::vector<std::uint32_t>& ring = data.rings[ring_index];
      std::vector<CoedgeSpec> coedges;
      coedges.reserve(ring.size());
      for (std::size_t i = 0; i < ring.size(); ++i) {
        const std::uint32_t from = ring[i];
        const std::uint32_t to = ring[(i + 1u) % ring.size()];
        const std::uint32_t a = std::min(from, to);
        const std::uint32_t b = std::max(from, to);
        const auto found = edge_lookup.find(std::make_pair(a, b));
        const std::uint32_t edge_index = found->second;
        // Pcurves follow the edge from v0 to v1, not the coedge traversal; the
        // coedge's reversed flag already records traversal direction. This is
        // the convention split_edge relies on.
        const Point2d p0 =
            plane_frame_to_2d(data.frame, welder.points()[a]);
        const Point2d p1 =
            plane_frame_to_2d(data.frame, welder.points()[b]);
        const std::uint32_t pcurve = editor.add_pcurve(Pcurve2d{p0, p1});
        coedges.push_back(
            CoedgeSpec{edge_ids[edge_index], from != a, pcurve, tol});
      }
      const Result<LoopId> loop =
          editor.add_loop(face_ids[face], coedges, ring_index == 0u);
      if (!loop) return loop.error();
    }
  }

  detail::assemble_disjoint_set components(face_data.size());
  for (const detail::assemble_edge_record& edge : edge_records) {
    for (std::size_t i = 1; i < edge.uses.size(); ++i) {
      components.unite(edge.uses.front().face, edge.uses[i].face);
    }
  }
  std::map<std::uint32_t, std::vector<std::uint32_t>> grouped;
  for (std::uint32_t face = 0; face < face_data.size(); ++face) {
    grouped[components.find(face)].push_back(face);
  }

  std::vector<detail::assemble_component> closed_components;
  closed_components.reserve(grouped.size());
  for (auto& entry : grouped) {
    detail::assemble_component component;
    component.faces = std::move(entry.second);
    const Result<double> volume = detail::assemble_component_volume(
        component.faces, face_data, welder.points());
    if (!volume) return volume.error();
    component.volume = volume.value();
    component.sample = detail::assemble_component_sample(
        component.faces, face_data, welder.points());
    closed_components.push_back(std::move(component));
  }

  std::vector<std::uint32_t> outer_components;
  std::vector<std::uint32_t> cavity_components;
  for (std::uint32_t i = 0; i < closed_components.size(); ++i) {
    if (!std::isfinite(closed_components[i].volume)) {
      return Error{Errc::numeric_failure,
                   "assemble_solid: component volume is not finite"};
    }
    if (closed_components[i].volume > 0.0) {
      outer_components.push_back(i);
    } else if (closed_components[i].volume < 0.0) {
      cavity_components.push_back(i);
    } else {
      return Error{Errc::degenerate,
                   "assemble_solid: component has zero volume"};
    }
  }
  if (outer_components.empty()) {
    return Error{Errc::invalid_topology,
                 "assemble_solid: no outward-oriented component"};
  }

  std::vector<std::vector<std::uint32_t>> shell_faces(outer_components.size());
  for (std::size_t i = 0; i < outer_components.size(); ++i) {
    shell_faces[i] = closed_components[outer_components[i]].faces;
  }
  for (const std::uint32_t cavity : cavity_components) {
    std::vector<std::size_t> containers;
    for (std::size_t i = 0; i < outer_components.size(); ++i) {
      if (detail::assemble_point_inside(
              closed_components[cavity].sample,
              closed_components[outer_components[i]].faces, face_data,
              welder.points())) {
        containers.push_back(i);
      }
    }
    if (containers.empty()) {
      return Error{Errc::invalid_topology,
                   "assemble_solid: cavity has no containing outer shell"};
    }
    if (containers.size() > 1u) {
      return Error{
          Errc::invalid_topology,
          "assemble_solid: cavity is contained by multiple outer shells"};
    }
    std::vector<std::uint32_t>& destination = shell_faces[containers.front()];
    const std::vector<std::uint32_t>& cavity_faces =
        closed_components[cavity].faces;
    destination.insert(destination.end(), cavity_faces.begin(),
                       cavity_faces.end());
  }

  for (const std::vector<std::uint32_t>& component_faces : shell_faces) {
    std::vector<FaceId> faces_of_shell;
    faces_of_shell.reserve(component_faces.size());
    for (const std::uint32_t face : component_faces) {
      faces_of_shell.push_back(face_ids[face]);
    }
    const Result<ShellId> shell = editor.add_shell(faces_of_shell, true);
    if (!shell) return shell.error();
    const std::array<ShellId, 1> shells{shell.value()};
    const Result<SolidId> solid = editor.add_solid(shells, tol);
    if (!solid) return solid.error();
  }
  return editor.commit();
}

}  // namespace v1
}  // namespace geom
