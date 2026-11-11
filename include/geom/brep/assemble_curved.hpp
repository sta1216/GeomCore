#pragma once

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <optional>
#include <span>
#include <unordered_map>
#include <utility>
#include <variant>
#include <vector>

#include <geom/brep/detail/curve_polyline.hpp>
#include <geom/brep/face_uv.hpp>
#include <geom/brep/solid_editor.hpp>
#include <geom/brep/topology.hpp>
#include <geom/brep/trim_face.hpp>
#include <geom/brep/validate.hpp>
#include <geom/core/error.hpp>
#include <geom/core/result.hpp>
#include <geom/curves/analytic_curve.hpp>
#include <geom/primitives/point.hpp>
#include <geom/primitives/vec.hpp>
#include <geom/surfaces/analytic.hpp>
#include <geom/surfaces/plane.hpp>

namespace geom {
inline namespace v1 {

// One trimmed face handed to the curved assembler. Every loop lists only its
// corner points (the ring is implicitly closed); `loops_uv` holds the same
// rings in the surface's parameter domain and `arc_spans` is parallel to the
// segment point[i] -> point[(i + 1) % size]. A segment is a circular arc
// exactly when its span entry is engaged; all other segments are straight.
//
// The pieces of a loop share endpoints only by position: the assembler welds
// them within the tolerance and never trusts the caller's indexing.
struct CurvedAssembledFace {
  SurfaceGeometry surface;
  std::vector<std::vector<Point3d>> loops;      // outer first
  std::vector<std::vector<Point2d>> loops_uv;   // same sizes as loops
  std::vector<std::vector<std::optional<CurvedArcSpan>>> arc_spans;  // same shape; may be empty
};

namespace detail {
namespace assemble_curved_detail {

// ---------------------------------------------------------------------------
// Small geometric helpers.
// ---------------------------------------------------------------------------

inline double assemble_curved_length(const Vec3d& v) noexcept { return std::sqrt(dot(v, v)); }

inline double assemble_curved_distance(const Point3d& a, const Point3d& b) noexcept {
  return assemble_curved_length(a - b);
}

inline Vec3d assemble_curved_unit(Vec3d v) noexcept {
  const double scale = std::max({std::abs(v[0]), std::abs(v[1]), std::abs(v[2])});
  if (!(scale > 0.0) || !std::isfinite(scale)) return Vec3d{0.0, 0.0, 0.0};
  v = v * (1.0 / scale);
  const double length = assemble_curved_length(v);
  if (!(length > 0.0)) return Vec3d{0.0, 0.0, 0.0};
  return v * (1.0 / length);
}

inline double assemble_curved_wrapped_delta(double delta) noexcept {
  const double turn = 2.0 * kFaceUvPi;
  delta = std::fmod(delta, turn);
  if (delta <= -kFaceUvPi) {
    delta += turn;
  } else if (delta > kFaceUvPi) {
    delta -= turn;
  }
  return delta;
}

inline bool assemble_curved_circle_normal(const CircleCurve3& circle, Vec3d& normal) noexcept {
  normal = assemble_curved_unit(cross(circle.x_axis, circle.y_axis));
  return assemble_curved_length(normal) > 0.0;
}

inline bool assemble_curved_same_circle(const CircleCurve3& lhs, const CircleCurve3& rhs,
                                        double tolerance) noexcept {
  if (!(lhs.radius > 0.0) || !(rhs.radius > 0.0)) return false;
  if (!(std::abs(lhs.radius - rhs.radius) <= tolerance)) return false;
  const double slack = std::max(tolerance, 1e-12 * std::max({lhs.radius, rhs.radius, 1.0}));
  if (!(assemble_curved_distance(lhs.center, rhs.center) <= slack)) return false;
  Vec3d lhs_normal{0.0, 0.0, 0.0};
  Vec3d rhs_normal{0.0, 0.0, 0.0};
  if (!assemble_curved_circle_normal(lhs, lhs_normal)) return false;
  if (!assemble_curved_circle_normal(rhs, rhs_normal)) return false;
  return assemble_curved_length(cross(lhs_normal, rhs_normal)) <= tolerance;
}

inline bool assemble_curved_same_surface(const SurfaceGeometry& lhs, const SurfaceGeometry& rhs,
                                         double tolerance) noexcept {
  if (lhs.index() != rhs.index()) return false;
  if (const Plane3* a = std::get_if<Plane3>(&lhs)) {
    const Plane3& b = std::get<Plane3>(rhs);
    const Vec3d an = assemble_curved_unit(a->normal);
    const Vec3d bn = assemble_curved_unit(b.normal);
    const double an_length = assemble_curved_length(a->normal);
    const double bn_length = assemble_curved_length(b.normal);
    if (!(assemble_curved_length(an) > 0.0) || !(assemble_curved_length(bn) > 0.0)) return false;
    const double alignment = dot(an, bn);
    if (std::abs(alignment) < 1.0 - tolerance) return false;
    const double lhs_offset = a->offset / an_length;
    const double rhs_offset = b.offset / bn_length;
    const double delta = alignment >= 0.0 ? lhs_offset - rhs_offset : lhs_offset + rhs_offset;
    return std::abs(delta) <= tolerance;
  }
  if (const Cylinder3* a = std::get_if<Cylinder3>(&lhs)) {
    const Cylinder3& b = std::get<Cylinder3>(rhs);
    return std::abs(a->radius - b.radius) <= tolerance &&
           assemble_curved_length(a->axis - b.axis) <= tolerance &&
           assemble_curved_distance(a->origin, b.origin) <= tolerance &&
           assemble_curved_length(a->x_axis - b.x_axis) <= tolerance;
  }
  if (const Cone3* a = std::get_if<Cone3>(&lhs)) {
    const Cone3& b = std::get<Cone3>(rhs);
    return std::abs(a->radius - b.radius) <= tolerance &&
           std::abs(a->half_angle - b.half_angle) <= tolerance &&
           assemble_curved_length(a->axis - b.axis) <= tolerance &&
           assemble_curved_distance(a->origin, b.origin) <= tolerance &&
           assemble_curved_length(a->x_axis - b.x_axis) <= tolerance;
  }
  if (const Sphere3* a = std::get_if<Sphere3>(&lhs)) {
    const Sphere3& b = std::get<Sphere3>(rhs);
    return std::abs(a->radius - b.radius) <= tolerance &&
           assemble_curved_distance(a->center, b.center) <= tolerance &&
           assemble_curved_length(a->axis - b.axis) <= tolerance &&
           assemble_curved_length(a->x_axis - b.x_axis) <= tolerance;
  }
  if (const Torus3* a = std::get_if<Torus3>(&lhs)) {
    const Torus3& b = std::get<Torus3>(rhs);
    return std::abs(a->major_radius - b.major_radius) <= tolerance &&
           std::abs(a->minor_radius - b.minor_radius) <= tolerance &&
           assemble_curved_distance(a->center, b.center) <= tolerance &&
           assemble_curved_length(a->axis - b.axis) <= tolerance &&
           assemble_curved_length(a->x_axis - b.x_axis) <= tolerance;
  }
  return false;
}

// ---------------------------------------------------------------------------
// Vertex welding. Endpoints are the only vertices; equal positions within the
// tolerance must reach the same welded id, so the grid is searched in the
// 27 neighbouring cells as well as the point's own cell.
// ---------------------------------------------------------------------------

struct assemble_curved_cell {
  std::int64_t x = 0;
  std::int64_t y = 0;
  std::int64_t z = 0;

  bool operator==(const assemble_curved_cell& other) const noexcept {
    return x == other.x && y == other.y && z == other.z;
  }
};

struct assemble_curved_cell_hash {
  std::size_t operator()(const assemble_curved_cell& cell) const noexcept {
    std::uint64_t hash = 1469598103934665603ull;
    hash = (hash ^ static_cast<std::uint64_t>(cell.x)) * 1099511628211ull;
    hash = (hash ^ static_cast<std::uint64_t>(cell.y)) * 1099511628211ull;
    hash = (hash ^ static_cast<std::uint64_t>(cell.z)) * 1099511628211ull;
    return static_cast<std::size_t>(hash ^ (hash >> 29));
  }
};

inline std::int64_t assemble_curved_cell_coordinate(double value, double size) noexcept {
  double cell = std::floor(value / size);
  if (!(cell < 9.0e18)) cell = 9.0e18;
  if (!(cell > -9.0e18)) cell = -9.0e18;
  return static_cast<std::int64_t>(cell);
}

class assemble_curved_welder {
 public:
  explicit assemble_curved_welder(double tolerance) noexcept : tolerance_(tolerance) {}

  std::uint32_t weld(const Point3d& point) {
    const assemble_curved_cell base = cell_of(point);
    std::uint32_t best = kNoVertex;
    double best_distance = std::numeric_limits<double>::infinity();
    for (std::int64_t dx = -1; dx <= 1; ++dx) {
      for (std::int64_t dy = -1; dy <= 1; ++dy) {
        for (std::int64_t dz = -1; dz <= 1; ++dz) {
          const assemble_curved_cell cell{base.x + dx, base.y + dy, base.z + dz};
          const auto found = grid_.find(cell);
          if (found == grid_.end()) continue;
          for (const std::uint32_t index : found->second) {
            const double distance = assemble_curved_distance(points_[index], point);
            if (distance <= tolerance_ && distance < best_distance) {
              best_distance = distance;
              best = index;
            }
          }
        }
      }
    }
    if (best != kNoVertex) return best;
    const std::uint32_t index = static_cast<std::uint32_t>(points_.size());
    points_.push_back(point);
    grid_[cell_of(point)].push_back(index);
    return index;
  }

  const std::vector<Point3d>& points() const noexcept { return points_; }

 private:
  static constexpr std::uint32_t kNoVertex = 0xFFFFFFFFu;

  assemble_curved_cell cell_of(const Point3d& point) const noexcept {
    const double size = 2.0 * tolerance_;
    return assemble_curved_cell{
        assemble_curved_cell_coordinate(point[0], size),
        assemble_curved_cell_coordinate(point[1], size),
        assemble_curved_cell_coordinate(point[2], size)};
  }

  double tolerance_ = 0.0;
  std::vector<Point3d> points_;
  std::unordered_map<assemble_curved_cell, std::vector<std::uint32_t>, assemble_curved_cell_hash>
      grid_;
};

// ---------------------------------------------------------------------------
// Loop splitting: maximal same-circle arc runs and single straight segments.
// ---------------------------------------------------------------------------

struct assemble_curved_piece {
  bool arc = false;
  std::size_t start = 0;  // index of the first point of the piece in its loop
  std::size_t count = 1;  // number of loop segments in the piece
  double sweep = 0.0;     // signed sweep of an arc piece
  CircleCurve3 circle;
  std::uint32_t v0 = 0;  // welded id of the loop-order start point
  std::uint32_t v1 = 0;  // welded id of the loop-order end point
  std::size_t use = 0;   // index into the collected uses
};

inline Result<assemble_curved_piece> assemble_curved_make_run(
    const std::vector<Point3d>& points, const std::vector<std::optional<CurvedArcSpan>>& spans,
    std::size_t start, std::size_t length, double tolerance) {
  const std::size_t size = points.size();
  const CircleCurve3& circle = spans[start].value().circle;
  if (!(circle.radius > 0.0)) {
    return Error{Errc::invalid_topology, "assemble_curved_solid: arc run is inconsistent", 0u, 0u};
  }
  const double angle_tolerance =
      std::min(kFaceUvPi, std::max(1e-12, 2.0 * tolerance / circle.radius));
  double total = 0.0;
  int sign = 0;
  for (std::size_t offset = 0; offset < length; ++offset) {
    const std::size_t index = (start + offset) % size;
    const std::size_t next = (index + 1) % size;
    const double width = spans[index].value().t1 - spans[index].value().t0;
    if (!(width > 0.0) || width > 2.0 * kFaceUvPi + angle_tolerance) {
      return Error{Errc::invalid_topology, "assemble_curved_solid: arc run is inconsistent", 0u,
                   0u};
    }
    if (width >= 2.0 * kFaceUvPi - angle_tolerance) {
      // A single segment covering a whole turn: its endpoints coincide, so the
      // direction is not recoverable from them. Restricting this to a closed
      // piece keeps the assembler from inventing a direction for an open arc.
      if (length != 1) {
        return Error{Errc::invalid_topology, "assemble_curved_solid: arc run is inconsistent", 0u,
                     0u};
      }
      assemble_curved_piece piece;
      piece.arc = true;
      piece.start = start;
      piece.count = 1;
      piece.sweep = 2.0 * kFaceUvPi;
      piece.circle = circle;
      return piece;
    }
    const double begin_angle = face_uv_circle_angle(circle, points[index]);
    const double end_angle = face_uv_circle_angle(circle, points[next]);
    const double delta = assemble_curved_wrapped_delta(end_angle - begin_angle);
    int segment_sign = 0;
    if (std::abs(delta - width) <= angle_tolerance) {
      segment_sign = 1;
    } else if (std::abs(delta + width) <= angle_tolerance) {
      segment_sign = -1;
    } else {
      return Error{Errc::invalid_topology, "assemble_curved_solid: arc run is inconsistent", 0u,
                   0u};
    }
    if (sign != 0 && sign != segment_sign) {
      return Error{Errc::invalid_topology, "assemble_curved_solid: arc run is inconsistent", 0u,
                   0u};
    }
    sign = segment_sign;
    total += segment_sign > 0 ? width : -width;
  }
  assemble_curved_piece piece;
  piece.arc = true;
  piece.start = start;
  piece.count = length;
  piece.sweep = total;
  piece.circle = circle;
  return piece;
}

inline Result<std::vector<assemble_curved_piece>> assemble_curved_split_loop(
    const std::vector<Point3d>& points, const std::vector<std::optional<CurvedArcSpan>>& spans,
    double tolerance) {
  const std::size_t size = points.size();
  std::vector<bool> arc(size, false);
  if (!spans.empty()) {
    for (std::size_t i = 0; i < size; ++i) arc[i] = spans[i].has_value();
  }
  const auto continues = [&](std::size_t previous, std::size_t current) {
    if (!arc[previous] || !arc[current]) return false;
    return assemble_curved_same_circle(spans[previous].value().circle,
                                       spans[current].value().circle, tolerance);
  };

  std::size_t start = size;
  for (std::size_t i = 0; i < size; ++i) {
    const std::size_t previous = (i + size - 1) % size;
    if (!continues(previous, i)) {
      start = i;
      break;
    }
  }
  if (start == size) start = 0;  // every segment continues the previous one

  std::vector<assemble_curved_piece> pieces;
  std::size_t processed = 0;
  while (processed < size) {
    const std::size_t index = (start + processed) % size;
    if (!arc[index]) {
      assemble_curved_piece piece;
      piece.start = index;
      piece.count = 1;
      pieces.push_back(piece);
      ++processed;
      continue;
    }
    std::size_t length = 1;
    while (length < size - processed) {
      const std::size_t next = (start + processed + length) % size;
      const std::size_t previous = (next + size - 1) % size;
      if (!continues(previous, next)) break;
      ++length;
    }
    Result<assemble_curved_piece> run =
        assemble_curved_make_run(points, spans, index, length, tolerance);
    if (!run) return run.error();
    pieces.push_back(run.value());
    processed += length;
  }
  return pieces;
}

// ---------------------------------------------------------------------------
// Grouping by undirected edge key.
// ---------------------------------------------------------------------------

struct assemble_curved_use {
  std::size_t face = 0;
  std::size_t loop = 0;
  std::size_t piece = 0;
  std::uint32_t v0 = 0;
  std::uint32_t v1 = 0;
  bool arc = false;
  double sweep = 0.0;
  std::size_t edge = 0;
  bool reversed = false;
};

struct assemble_curved_edge_slot {
  std::uint32_t v0 = 0;
  std::uint32_t v1 = 0;
  std::vector<std::size_t> uses;
  bool arc = false;
  double sweep = 0.0;
  CircleCurve3 circle;
  double t0 = 0.0;
  double t1 = 1.0;
  std::uint32_t curve = kNoGeometry;
  EdgeId id{};
};

inline std::uint64_t assemble_curved_edge_key(std::uint32_t v0, std::uint32_t v1) noexcept {
  const std::uint32_t low = std::min(v0, v1);
  const std::uint32_t high = std::max(v0, v1);
  return (static_cast<std::uint64_t>(low) << 32) | static_cast<std::uint64_t>(high);
}

// Projects a polyline of 3D samples into the face's parameter space, keeping
// the periodic branch continuous with the fragment's uv ring.
inline Result<std::uint32_t> assemble_curved_add_pcurve(
    SolidEditor& editor, const FaceUvProjection& projection,
    const std::vector<Point3d>& samples, const Point2d& seed) {
  std::vector<Point2d> points;
  points.reserve(samples.size());
  double previous_u = seed[0];
  for (const Point3d& sample : samples) {
    const Result<Point2d> projected = face_uv_project(projection, sample, previous_u, true);
    if (!projected) return projected.error();
    Point2d uv = projected.value();
    if (projection.periodic_u()) uv[0] = face_uv_unwrap(uv[0], previous_u);
    previous_u = uv[0];
    points.push_back(uv);
  }
  Pcurve2d pcurve;
  pcurve.points = std::move(points);
  return editor.add_pcurve(std::move(pcurve));
}

inline Result<FaceUvProjection> assemble_curved_projection(const SurfaceGeometry& surface,
                                                           double tolerance) {
  if (std::holds_alternative<NurbsSurface>(surface)) {
    return Error{Errc::unsupported, "assemble_curved_solid: NURBS surfaces are not supported", 0u,
                 0u};
  }
  FaceUvProjection projection;
  projection.surface = &surface;
  projection.tolerance = tolerance;
  if (const Plane3* plane = std::get_if<Plane3>(&surface)) {
    const Result<PlaneFrame> frame = make_plane_frame(*plane);
    if (!frame) {
      return Error{Errc::invalid_topology,
                   "assemble_curved_solid: plane surface is degenerate", 0u, 0u};
    }
    projection.frame = frame.value();
  }
  return projection;
}

// ---------------------------------------------------------------------------
// The assembler itself.
// ---------------------------------------------------------------------------

inline Result<Brep> assemble_curved_impl(std::span<const CurvedAssembledFace> faces,
                                         double tolerance) {
  if (!std::isfinite(tolerance) || !(tolerance > 0.0)) {
    return Error{Errc::tolerance_violation,
                 "assemble_curved_solid: tolerance must be finite and positive", 0u, 0u};
  }
  if (faces.empty()) {
    return Error{Errc::invalid_topology, "assemble_curved_solid: at least one face is required",
                 0u, 0u};
  }

  // Step 1: validate the fragments before touching any geometry.
  for (std::size_t f = 0; f < faces.size(); ++f) {
    const CurvedAssembledFace& face = faces[f];
    const std::uint32_t face_index = static_cast<std::uint32_t>(f);
    if (face.loops.size() != face.loops_uv.size()) {
      return Error{Errc::invalid_topology,
                   "assemble_curved_solid: loop and uv loop counts differ", face_index, 0u};
    }
    if (!face.arc_spans.empty() && face.arc_spans.size() != face.loops.size()) {
      return Error{Errc::invalid_topology,
                   "assemble_curved_solid: arc provenance does not match the loops", face_index,
                   0u};
    }
    for (std::size_t l = 0; l < face.loops.size(); ++l) {
      const std::uint32_t loop_index = static_cast<std::uint32_t>(l);
      if (face.loops[l].size() < 3u) {
        return Error{Errc::invalid_topology,
                     "assemble_curved_solid: loop has fewer than three points", face_index,
                     loop_index};
      }
      if (face.loops_uv[l].size() != face.loops[l].size()) {
        return Error{Errc::invalid_topology,
                     "assemble_curved_solid: loop and uv loop sizes differ", face_index,
                     loop_index};
      }
      if (!face.arc_spans.empty() && face.arc_spans[l].size() != face.loops[l].size()) {
        return Error{Errc::invalid_topology,
                     "assemble_curved_solid: arc provenance does not match the loops", face_index,
                     loop_index};
      }
      for (const Point3d& point : face.loops[l]) {
        if (!face_uv_point3d_finite(point)) {
          return Error{Errc::non_finite, "assemble_curved_solid: loop point is not finite",
                       face_index, loop_index};
        }
      }
      for (const Point2d& uv : face.loops_uv[l]) {
        if (!face_uv_point2d_finite(uv)) {
          return Error{Errc::non_finite, "assemble_curved_solid: loop uv point is not finite",
                       face_index, loop_index};
        }
      }
    }
  }

  // Step 2: split every loop into arc runs and straight segments, and step 3:
  // weld the piece endpoints. Only those endpoints become vertices.
  assemble_curved_welder welder(tolerance);
  std::vector<std::vector<std::vector<assemble_curved_piece>>> layouts(faces.size());
  for (std::size_t f = 0; f < faces.size(); ++f) {
    const CurvedAssembledFace& face = faces[f];
    layouts[f].resize(face.loops.size());
    for (std::size_t l = 0; l < face.loops.size(); ++l) {
      const std::vector<std::optional<CurvedArcSpan>> empty;
      const std::vector<std::optional<CurvedArcSpan>>& spans =
          face.arc_spans.empty() ? empty : face.arc_spans[l];
      Result<std::vector<assemble_curved_piece>> pieces =
          assemble_curved_split_loop(face.loops[l], spans, tolerance);
      if (!pieces) return pieces.error();
      layouts[f][l] = std::move(pieces.value());
      const std::size_t size = face.loops[l].size();
      for (assemble_curved_piece& piece : layouts[f][l]) {
        piece.v0 = welder.weld(face.loops[l][piece.start]);
        piece.v1 = welder.weld(face.loops[l][(piece.start + piece.count) % size]);
      }
    }
  }
  if (welder.points().size() > 0xFFFFFFFFull) {
    return Error{Errc::invalid_topology, "assemble_curved_solid: too many vertices", 0u, 0u};
  }

  // Step 4: group the uses by undirected key and check that every edge is
  // shared by exactly two compatible pieces.
  std::vector<assemble_curved_use> uses;
  std::vector<assemble_curved_edge_slot> edges;
  std::unordered_map<std::uint64_t, std::size_t> edge_of_key;
  for (std::size_t f = 0; f < faces.size(); ++f) {
    for (std::size_t l = 0; l < layouts[f].size(); ++l) {
      for (std::size_t p = 0; p < layouts[f][l].size(); ++p) {
        assemble_curved_piece& piece = layouts[f][l][p];
        assemble_curved_use use;
        use.face = f;
        use.loop = l;
        use.piece = p;
        use.v0 = piece.v0;
        use.v1 = piece.v1;
        use.arc = piece.arc;
        use.sweep = piece.sweep;
        const std::uint64_t key = assemble_curved_edge_key(piece.v0, piece.v1);
        const auto found = edge_of_key.find(key);
        if (found == edge_of_key.end()) {
          use.edge = edges.size();
          assemble_curved_edge_slot slot;
          slot.v0 = piece.v0;
          slot.v1 = piece.v1;
          slot.arc = piece.arc;
          slot.sweep = piece.sweep;
          slot.circle = piece.circle;
          edges.push_back(slot);
          edge_of_key.emplace(key, use.edge);
        } else {
          use.edge = found->second;
        }
        piece.use = uses.size();
        uses.push_back(use);
        edges[use.edge].uses.push_back(uses.size() - 1);
      }
    }
  }

  for (const assemble_curved_edge_slot& slot : edges) {
    if (slot.uses.size() != 2u) {
      return Error{Errc::invalid_topology, "assemble_curved_solid: edge is used once", 0u, 0u};
    }
  }
  for (const assemble_curved_edge_slot& slot : edges) {
    for (std::size_t k = 1; k < slot.uses.size(); ++k) {
      assemble_curved_use& use = uses[slot.uses[k]];
      const assemble_curved_piece& piece = layouts[use.face][use.loop][use.piece];
      if (use.arc != slot.arc) {
        return Error{Errc::invalid_topology,
                     "assemble_curved_solid: shared edge curve mismatch", 0u, 0u};
      }
      if (slot.arc) {
        if (!assemble_curved_same_circle(piece.circle, slot.circle, tolerance)) {
          return Error{Errc::invalid_topology,
                       "assemble_curved_solid: shared edge curve mismatch", 0u, 0u};
        }
        const double sweep_tolerance =
            std::max(1e-12, 4.0 * tolerance / std::max(slot.circle.radius, 1e-12));
        if (!(std::abs(std::abs(piece.sweep) - std::abs(slot.sweep)) <= sweep_tolerance)) {
          return Error{Errc::invalid_topology,
                       "assemble_curved_solid: shared edge curve mismatch", 0u, 0u};
        }
        use.reversed = (piece.sweep < 0.0) != (slot.sweep < 0.0);
      } else {
        use.reversed = use.v0 != slot.v0;
      }
    }
  }

  // Step 5: shared vertices, edges and surfaces.
  SolidEditor editor;
  std::vector<VertexId> vertex_ids;
  vertex_ids.reserve(welder.points().size());
  for (const Point3d& point : welder.points()) {
    const Result<VertexId> vertex = editor.add_vertex(point, tolerance);
    if (!vertex) return vertex.error();
    vertex_ids.push_back(vertex.value());
  }

  std::vector<std::uint32_t> surface_ids(faces.size(), kNoGeometry);
  std::vector<const SurfaceGeometry*> unique_surfaces;
  std::vector<std::uint32_t> unique_surface_ids;
  for (std::size_t f = 0; f < faces.size(); ++f) {
    std::uint32_t id = kNoGeometry;
    for (std::size_t k = 0; k < unique_surfaces.size(); ++k) {
      if (assemble_curved_same_surface(faces[f].surface, *unique_surfaces[k], tolerance)) {
        id = unique_surface_ids[k];
        break;
      }
    }
    if (id == kNoGeometry) {
      id = editor.add_surface(faces[f].surface);
      unique_surfaces.push_back(&faces[f].surface);
      unique_surface_ids.push_back(id);
    }
    surface_ids[f] = id;
  }

  for (assemble_curved_edge_slot& slot : edges) {
    if (slot.arc) {
      Vec3d normal{0.0, 0.0, 0.0};
      if (!assemble_curved_circle_normal(slot.circle, normal)) {
        return Error{Errc::invalid_topology, "assemble_curved_solid: arc run is inconsistent", 0u,
                     0u};
      }
      const Vec3d start_radial =
          assemble_curved_unit(welder.points()[slot.v0] - slot.circle.center);
      if (!(assemble_curved_length(start_radial) > 0.0)) {
        return Error{Errc::invalid_topology, "assemble_curved_solid: arc run is inconsistent", 0u,
                     0u};
      }
      Vec3d y_axis = assemble_curved_unit(cross(normal, start_radial));
      if (!(assemble_curved_length(y_axis) > 0.0)) {
        return Error{Errc::invalid_topology, "assemble_curved_solid: arc run is inconsistent", 0u,
                     0u};
      }
      // The edge starts at the run's start point (angle zero) and parameterises
      // towards the run's sweep: a negative run flips the handedness instead of
      // moving the start angle.
      if (slot.sweep < 0.0) y_axis = y_axis * -1.0;
      double t1 = std::abs(slot.sweep);
      if (slot.v0 == slot.v1 && t1 >= 2.0 * kFaceUvPi - 1e-12) t1 = 2.0 * kFaceUvPi;
      slot.circle = CircleCurve3{slot.circle.center, start_radial, y_axis, slot.circle.radius};
      slot.t0 = 0.0;
      slot.t1 = t1;
      slot.curve = editor.add_curve(CurveGeometry{slot.circle});
    } else {
      const std::array<Point3d, 2> span_points{welder.points()[slot.v0],
                                               welder.points()[slot.v1]};
      const Result<NurbsCurve> curve = make_polyline_curve(span_points);
      if (!curve) return curve.error();
      slot.curve = editor.add_curve(CurveGeometry{curve.value()});
      slot.t0 = 0.0;
      slot.t1 = 1.0;
    }
    const Result<EdgeId> edge =
        editor.add_edge(vertex_ids[slot.v0], vertex_ids[slot.v1], slot.curve, slot.t0, slot.t1,
                        tolerance);
    if (!edge) return edge.error();
    slot.id = edge.value();
  }

  std::vector<FaceId> face_ids;
  face_ids.reserve(faces.size());
  for (std::size_t f = 0; f < faces.size(); ++f) {
    const Result<FaceId> face = editor.add_face(surface_ids[f], false, tolerance);
    if (!face) return face.error();
    face_ids.push_back(face.value());
  }

  // Step 6: one dense pcurve per coedge, always in the edge direction.
  for (std::size_t f = 0; f < faces.size(); ++f) {
    const CurvedAssembledFace& face = faces[f];
    const Result<FaceUvProjection> projection =
        assemble_curved_projection(face.surface, tolerance);
    if (!projection) return projection.error();
    for (std::size_t l = 0; l < face.loops.size(); ++l) {
      const std::size_t size = face.loops[l].size();
      std::vector<CoedgeSpec> coedges;
      coedges.reserve(layouts[f][l].size());
      for (const assemble_curved_piece& piece : layouts[f][l]) {
        const assemble_curved_use& use = uses[piece.use];
        const assemble_curved_edge_slot& slot = edges[use.edge];
        std::vector<Point3d> samples;
        if (slot.arc) {
          const double span = slot.t1 - slot.t0;
          const double ratio = std::min(std::min(tolerance, 1e-6) / slot.circle.radius, 1.0);
          const double step = 2.0 * std::acos(1.0 - ratio);
          double count = std::ceil(span / step);
          if (!std::isfinite(count) || count < 32.0) count = 32.0;
          if (count > 4096.0) count = 4096.0;
          const std::size_t segments = static_cast<std::size_t>(count);
          samples.reserve(segments + 1);
          for (std::size_t k = 0; k <= segments; ++k) {
            const double fraction =
                static_cast<double>(k) / static_cast<double>(segments);
            samples.push_back(slot.circle.eval(slot.t0 + span * fraction));
          }
        } else {
          samples.push_back(welder.points()[slot.v0]);
          samples.push_back(welder.points()[slot.v1]);
        }
        const std::size_t seed_index =
            use.reversed ? (piece.start + piece.count) % size : piece.start;
        const Result<std::uint32_t> pcurve = assemble_curved_add_pcurve(
            editor, projection.value(), samples, face.loops_uv[l][seed_index]);
        if (!pcurve) return pcurve.error();
        coedges.push_back(CoedgeSpec{slot.id, use.reversed, pcurve.value(), tolerance});
      }
      const Result<LoopId> loop = editor.add_loop(face_ids[f], coedges, l == 0);
      if (!loop) return loop.error();
    }
  }

  // Step 7 and 8: one closed shell, one solid, then commit.
  const Result<ShellId> shell = editor.add_shell(face_ids, true);
  if (!shell) return shell.error();
  const std::array<ShellId, 1> shells{shell.value()};
  const Result<SolidId> solid = editor.add_solid(shells, tolerance);
  if (!solid) return solid.error();
  Result<Brep> committed = editor.commit();
  if (!committed) return committed.error();

  // Step 9: surface the first topology finding.
  const ValidationReport report = validate_topology(committed.value());
  if (!report.ok()) {
    const ValidationIssue& issue = report.issues.front();
    return Error{Errc::invalid_topology, issue.message, issue.entity, 0u};
  }
  return std::move(committed).value();
}

}  // namespace assemble_curved_detail
}  // namespace detail

inline Result<Brep> assemble_curved_solid(std::span<const CurvedAssembledFace> faces,
                                          double tolerance) {
  return detail::assemble_curved_detail::assemble_curved_impl(faces, tolerance);
}

}  // namespace v1
}  // namespace geom
