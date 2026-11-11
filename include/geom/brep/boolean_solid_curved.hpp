#pragma once

// Boolean operations between analytic B-rep solids whose boundary pieces meet
// in trimmed faces.
//
// The driver is the curved counterpart of `boolean_solid` (which is planar
// only). It pairs face bounding boxes, asks the surface/surface intersector for
// a cut curve per pair, splits every face of both operands along the cuts and
// keeps the fragments that the requested membership rule puts on the result
// boundary. Cut curves are chordal polylines; whenever a pair's intersection is
// an analytic circle (a cylinder or sphere meeting a plane) the chords are
// healed back into exact `CircleCurve3` arc spans, exactly as
// `split_solid_by_plane` does, so a trimmed cap keeps its exact area instead of
// an inscribed polygon. The healing helpers are shared with `split_solid.hpp`
// (they live in `detail` there and are included, not modified).
//
// Classification is two-sided, mirroring `detail::solid_boolean_classify_face`
// in `boolean_solid.hpp`: a fragment is kept when the membership rule differs on
// its two sides, and its normal is flipped when the result interior lies on the
// side its own surface normal points to. Sampling membership once per fragment is
// enough here because the other operand's mesh must not contain the
// representative point (that case is reported as a numeric failure).
//
// A flipped fragment is emitted with the flipped normal baked into its geometry
// rather than through a face flag: a planar piece is re-expressed on the opposite
// plane, its rings are reversed and its parameters reprojected, so every output
// face is counter-clockwise about its surface normal with `reversed == false`,
// exactly like `boolean_solid`'s output. That keeps two faces sharing a cut
// edge walking it in opposite directions, and keeps `curved_solid_volume`'s
// planar area positively oriented.
//
// Deferred: coplanar and tangent face contact, NURBS and cone operands, flipped
// fragments on curved surfaces (no reparameterization here), and results that
// fall apart into several connected components. All of them are reported as
// `unsupported` instead of producing a wrong solid.

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <numeric>
#include <span>
#include <string_view>
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

// The surface kinds the driver can cut. NURBS patches have no analytic
// parameterization here and a cone's apex produces tangent loops that the
// splitter cannot yet represent, so both are rejected up front.
struct CurvedBooleanSurfaceGate {
  bool supported = false;
  const char* message = nullptr;
};

inline CurvedBooleanSurfaceGate curved_boolean_surface_gate(
    const SurfaceGeometry& surface) noexcept {
  if (std::holds_alternative<Plane3>(surface) || std::holds_alternative<Cylinder3>(surface) ||
      std::holds_alternative<Sphere3>(surface) || std::holds_alternative<Torus3>(surface)) {
    return CurvedBooleanSurfaceGate{true, nullptr};
  }
  if (std::holds_alternative<NurbsSurface>(surface)) {
    return CurvedBooleanSurfaceGate{
        false, "boolean_solid_curved: NURBS operands are not supported"};
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
// Longest 3D edge of a curved face's interior mesh, as a fraction of the face
// extent. Interior refinement never splits a boundary edge, so two faces that
// share an edge keep sharing their whole boundary polyline and the mesh of a
// closed solid stays closed.
constexpr double kCurvedBooleanMeshEdge = 0.05;
constexpr int kCurvedBooleanMeshDepth = 4;

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

// Appends one face's triangles. Planar faces are triangulated from the projected
// loops directly; a curved face is additionally refined in its parameter domain
// so an interior chord cannot cut a large corner off the surface.
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

  if (info.projection.is_plane()) {
    const std::uint32_t base = static_cast<std::uint32_t>(positions.size());
    for (const Point2d& uv : ring_points) positions.push_back(face_uv_eval(info.projection, uv));
    for (const std::array<std::uint32_t, 3>& triangle : triangulation.triangles) {
      indices.push_back(
          {base + triangle[0], base + triangle[1], base + triangle[2]});
    }
    return {};
  }

  std::vector<Point2d> uv_points = ring_points;
  std::vector<Point3d> uv_positions;
  uv_positions.reserve(ring_count * 8u);
  for (const Point2d& uv : ring_points) uv_positions.push_back(face_uv_eval(info.projection, uv));

  // Every edge between two consecutive ring vertices is a boundary edge that a
  // neighbouring face reproduces exactly.
  std::vector<std::uint64_t> boundary;
  boundary.reserve(ring_count);
  for (std::size_t r = 0; r + 1 < ring_starts.size(); ++r) {
    const std::uint32_t begin = ring_starts[r];
    const std::uint32_t end = ring_starts[r + 1];
    for (std::uint32_t i = begin; i < end; ++i) {
      const std::uint32_t j = i + 1 < end ? i + 1 : begin;
      boundary.push_back(curved_boolean_edge_key(i, j));
    }
  }
  std::sort(boundary.begin(), boundary.end());

  double extent = 0.0;
  {
    std::array<double, 3> lo{};
    std::array<double, 3> hi{};
    for (int axis = 0; axis < 3; ++axis) {
      lo[static_cast<std::size_t>(axis)] = uv_positions.front()[axis];
      hi[static_cast<std::size_t>(axis)] = uv_positions.front()[axis];
    }
    for (const Point3d& point : uv_positions) {
      for (int axis = 0; axis < 3; ++axis) {
        const std::size_t slot = static_cast<std::size_t>(axis);
        lo[slot] = std::min(lo[slot], point[axis]);
        hi[slot] = std::max(hi[slot], point[axis]);
      }
    }
    for (std::size_t axis = 0; axis < 3; ++axis) extent = std::max(extent, hi[axis] - lo[axis]);
  }
  const double target = std::max(kCurvedBooleanMeshEdge * extent, 4.0 * kCurvedBooleanMeshSag);

  const auto is_boundary = [&](std::uint32_t a, std::uint32_t b) noexcept {
    if (a >= ring_count || b >= ring_count) return false;
    return std::binary_search(boundary.begin(), boundary.end(), curved_boolean_edge_key(a, b));
  };
  const auto midpoint = [&](std::uint32_t x, std::uint32_t y) {
    const Point2d uv{0.5 * (uv_points[x][0] + uv_points[y][0]),
                     0.5 * (uv_points[x][1] + uv_points[y][1])};
    uv_points.push_back(uv);
    uv_positions.push_back(face_uv_eval(info.projection, uv));
    return static_cast<std::uint32_t>(uv_points.size() - 1);
  };
  std::vector<std::array<std::uint32_t, 3>> refined;
  refined.reserve(triangulation.triangles.size() * 8u);
  const auto refine = [&](auto&& self, std::uint32_t a, std::uint32_t b, std::uint32_t c,
                          int depth) -> void {
    bool ab = is_boundary(a, b);
    bool bc = is_boundary(b, c);
    bool ca = is_boundary(c, a);
    const int interior = (!ab ? 1 : 0) + (!bc ? 1 : 0) + (!ca ? 1 : 0);
    double longest = 0.0;
    longest = std::max(longest, face_uv_distance(uv_positions[a], uv_positions[b]));
    longest = std::max(longest, face_uv_distance(uv_positions[b], uv_positions[c]));
    longest = std::max(longest, face_uv_distance(uv_positions[c], uv_positions[a]));
    if (interior == 0 || depth <= 0 || !(longest > target)) {
      refined.push_back({a, b, c});
      return;
    }
    // Rotate the corners so the interior edges match a canonical pattern.
    if (interior == 1) {
      if (!bc) {
        std::swap(a, b);
        std::swap(b, c);
        ab = is_boundary(a, b);
        bc = is_boundary(b, c);
        ca = is_boundary(c, a);
      } else if (!ca) {
        std::swap(b, c);
        std::swap(a, b);
        ab = is_boundary(a, b);
        bc = is_boundary(b, c);
        ca = is_boundary(c, a);
      }
      const std::uint32_t m_ab = midpoint(a, b);
      self(self, a, m_ab, c, depth - 1);
      self(self, m_ab, b, c, depth - 1);
      return;
    }
    if (interior == 2) {
      if (ab && bc) {
        std::swap(a, b);
        std::swap(b, c);
      } else if (ab && ca) {
        std::swap(b, c);
        std::swap(a, b);
      }
      // Boundary edge ca; interior edges ab and bc.
      const std::uint32_t m_ab = midpoint(a, b);
      const std::uint32_t m_bc = midpoint(b, c);
      self(self, a, m_ab, m_bc, depth - 1);
      self(self, a, m_bc, c, depth - 1);
      self(self, m_ab, b, m_bc, depth - 1);
      return;
    }
    const std::uint32_t m_ab = midpoint(a, b);
    const std::uint32_t m_bc = midpoint(b, c);
    const std::uint32_t m_ca = midpoint(c, a);
    self(self, a, m_ab, m_ca, depth - 1);
    self(self, m_ab, b, m_bc, depth - 1);
    self(self, m_ca, m_bc, c, depth - 1);
    self(self, m_ab, m_bc, m_ca, depth - 1);
  };
  for (const std::array<std::uint32_t, 3>& triangle : triangulation.triangles) {
    refine(refine, triangle[0], triangle[1], triangle[2], kCurvedBooleanMeshDepth);
  }

  const std::uint32_t base = static_cast<std::uint32_t>(positions.size());
  positions.insert(positions.end(), uv_positions.begin(), uv_positions.end());
  for (const std::array<std::uint32_t, 3>& triangle : refined) {
    indices.push_back({base + triangle[0], base + triangle[1], base + triangle[2]});
  }
  return {};
}

// Triangulates an operand's faces into one mesh for membership tests. The B-rep
// tessellator cannot sample a partial analytic cylinder (it requires closed cap
// circles on the lateral face), so the driver samples the parameter domains
// itself.
inline Result<TriMesh> curved_boolean_solid_mesh(const std::vector<CurvedBooleanFaceInfo>& infos) {
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
  return mesh;
}

// Unwraps one parameter axis of a cut chain so consecutive parameters differ by
// less than half a period, then slides the chain by whole turns next to the
// face's own domain.
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
  const double middle = 0.5 * (points.front()[axis] + points.back()[axis]);
  const double center = 0.5 * (range[0] + range[1]);
  const double turns = std::floor((center - middle) / kCurvedBooleanTwoPi + 0.5);
  for (Point2d& point : points) {
    point[axis] += turns * kCurvedBooleanTwoPi;
    point[axis] = std::min(std::max(point[axis], range[0]), range[1]);
  }
}

// Turns one intersector curve into the cut polyline of one side. `segment.ua` /
// `ub` are the A-side parameters at the segment's two endpoints, `va` / `vb` the
// B-side ones.
inline std::vector<Point2d> curved_boolean_cut_params(const SsiCurve& curve, bool a_side,
                                                      const CurvedBooleanFaceInfo& face,
                                                      double tolerance) {
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
  return unique;
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
// instead of by its body. `split_solid_by_plane` documents the same failure mode
// and makes the same choice.
inline Point2d curved_boolean_representative(const std::vector<Point2d>& ring) {
  Polygon2d polygon;
  polygon.vertices = ring;
  polygon.ring_starts = {0u, static_cast<std::uint32_t>(ring.size())};
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
  if (!have) {
    for (const Point2d& point : ring) {
      representative[0] += point[0];
      representative[1] += point[1];
    }
    if (!ring.empty()) {
      const double inverse = 1.0 / static_cast<double>(ring.size());
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

  // Membership is decided against a triangulated copy of the other operand.
  const Result<TriMesh> mesh_a = detail::curved_boolean_solid_mesh(infos_a);
  if (!mesh_a.has_value()) return mesh_a.error();
  const Result<TriMesh> mesh_b = detail::curved_boolean_solid_mesh(infos_b);
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
      if (detail::curved_boolean_pair_circle(*face_a.surface, face_a.projection.frame,
                                             *face_b.surface, face_b.projection.frame, circle)) {
        circles_a[ia].push_back(circle);
        circles_b[ib].push_back(circle);
      }

      for (const SsiCurve& curve : intersected.value().curves) {
        std::vector<Point2d> uv_a =
            detail::curved_boolean_cut_params(curve, true, face_a, tolerance);
        if (detail::curved_boolean_push_cut(cuts_a[ia], std::move(uv_a), tolerance)) {
          ++report.cuts;
        }
        std::vector<Point2d> uv_b =
            detail::curved_boolean_cut_params(curve, false, face_b, tolerance);
        if (detail::curved_boolean_push_cut(cuts_b[ib], std::move(uv_b), tolerance)) {
          ++report.cuts;
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
            detail::curved_boolean_representative(fragment.loops_uv.front());
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
  if (detail::curved_boolean_component_count(assembled_faces, weld) > 1) {
    return Error{Errc::unsupported,
                 "boolean_solid_curved: result is not a single connected component", 0u, 0u};
  }

  Result<Brep> assembled = assemble_curved_solid(assembled_faces, tolerance);
  if (!assembled.has_value()) return assembled.error();
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
