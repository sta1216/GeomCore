#pragma once

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <optional>
#include <variant>
#include <vector>

#include <geom/brep/plane_frame.hpp>
#include <geom/brep/pcurve2d.hpp>
#include <geom/brep/topology.hpp>
#include <geom/core/error.hpp>
#include <geom/core/result.hpp>
#include <geom/curves/analytic_curve.hpp>
#include <geom/curves/nurbs_curve.hpp>
#include <geom/primitives/point.hpp>
#include <geom/primitives/vec.hpp>
#include <geom/surfaces/analytic.hpp>
#include <geom/surfaces/plane.hpp>

namespace geom {
inline namespace v1 {

// One boundary loop of a face, lifted into the parameter domain of the face's
// surface. The points are ordered the way the loop is traversed; the loop is
// implicitly closed, so the first point is not repeated at the end.
//
// `outer` mirrors Loop::outer. face_loops_uv() reports the outer loop first.
//
// `arcs` is parallel to the loop's segments: entry i describes the segment from
// points[i] to points[(i + 1) % points.size()], and the last entry is the
// closing segment. An entry is engaged exactly when its segment is a true
// circular arc from a coedge whose edge curve is a CircleCurve3; straight
// segments (cut lines, diameters, line pcurves) stay disengaged even when both
// endpoints lie on some circle.
//
// `t0` is the start azimuth in [0, 2*pi) and `t1` is `t0` plus the signed sweep
// of the traversal, so `|t1 - t0| < pi` and `t1` may fall outside [0, 2*pi).
struct FaceUvArcRef {
  CircleCurve3 circle;
  double t0 = 0.0;
  double t1 = 0.0;
};

struct FaceUvLoop {
  std::vector<Point2d> points;
  bool outer = false;
  std::vector<std::optional<FaceUvArcRef>> arcs;
};

namespace detail {

inline constexpr double kFaceUvPi = 3.14159265358979323846;

// Circle and NURBS edges are sampled with this many segments, so they
// contribute 33 points before junction deduplication.
inline constexpr std::size_t kFaceUvSegments = 32u;

inline bool face_uv_finite(double value) noexcept { return std::isfinite(value); }

inline bool face_uv_point2d_finite(const Point2d& point) noexcept {
  return face_uv_finite(point[0]) && face_uv_finite(point[1]);
}

inline bool face_uv_point3d_finite(const Point3d& point) noexcept {
  return face_uv_finite(point[0]) && face_uv_finite(point[1]) && face_uv_finite(point[2]);
}

inline double face_uv_distance(const Point3d& a, const Point3d& b) noexcept {
  return std::hypot(a[0] - b[0], a[1] - b[1], a[2] - b[2]);
}

// Slides `value` by whole turns so that it lies within half a turn of
// `reference`. Both arguments must be finite.
inline double face_uv_unwrap(double value, double reference) noexcept {
  const double turn = 2.0 * kFaceUvPi;
  double result = value;
  while (result - reference > kFaceUvPi) result -= turn;
  while (reference - result > kFaceUvPi) result += turn;
  return result;
}

// The surface a face is parameterized on, plus the frame used to project onto a
// plane. A cylinder, cone, sphere or torus is periodic in u; a torus is also
// periodic in v.
struct FaceUvProjection {
  const SurfaceGeometry* surface = nullptr;
  PlaneFrame frame{};
  double tolerance = 0.0;

  bool is_plane() const noexcept {
    return surface != nullptr && std::holds_alternative<Plane3>(*surface);
  }

  bool periodic_u() const noexcept { return surface != nullptr && !is_plane(); }

  bool periodic_v() const noexcept {
    return surface != nullptr && std::holds_alternative<Torus3>(*surface);
  }
};

inline Result<Point2d> face_uv_project_plane(const PlaneFrame& frame, const Point3d& point) {
  const Point2d uv = plane_frame_to_2d(frame, point);
  if (!face_uv_point2d_finite(uv)) {
    return Error{Errc::non_finite, "face_loops_uv: plane projection is not finite", 0u, 0u};
  }
  return uv;
}

inline Result<Point2d> face_uv_project_cylinder(const Cylinder3& surface, const Point3d& point) {
  const Vec3d offset = point - surface.origin;
  const double v = dot(offset, surface.axis);
  const Vec3d y_axis = cross(surface.axis, surface.x_axis);
  const Vec3d radial = offset - surface.axis * v;
  const double x = dot(radial, surface.x_axis);
  const double y = dot(radial, y_axis);
  // Points on the axis have no azimuth: rho is zero and the preimage is not
  // recoverable.
  if (!(std::hypot(x, y) > 0.0)) {
    return Error{Errc::invalid_topology,
                 "face_loops_uv: cylinder projection is singular on the axis", 0u, 0u};
  }
  const double u = std::atan2(y, x);
  if (!face_uv_finite(u) || !face_uv_finite(v)) {
    return Error{Errc::non_finite, "face_loops_uv: cylinder inverse is not finite", 0u, 0u};
  }
  return Point2d{u, v};
}

inline Result<Point2d> face_uv_project_cone(const Cone3& surface, const Point3d& point) {
  const Vec3d offset = point - surface.origin;
  const double v = dot(offset, surface.axis);
  const double rho = surface.radius - v * std::tan(surface.half_angle);
  // At and beyond the apex the radius vanishes (or turns negative on the other
  // nappe), so the azimuth is undefined. The threshold is relative to the base
  // radius so that a rounding-level positive rho at the apex does not turn into
  // an arbitrary azimuth.
  const double min_rho = surface.radius * 1e-12;
  if (!(rho > min_rho)) {
    return Error{Errc::unsupported,
                 "face_loops_uv: cone projection is undefined at or beyond the apex", 0u, 0u};
  }
  const Vec3d y_axis = cross(surface.axis, surface.x_axis);
  const Vec3d radial = offset - surface.axis * v;
  const double x = dot(radial, surface.x_axis);
  const double y = dot(radial, y_axis);
  if (!(std::hypot(x, y) > 0.0)) {
    return Error{Errc::invalid_topology,
                 "face_loops_uv: cone projection is singular on the axis", 0u, 0u};
  }
  const double u = std::atan2(y, x);
  if (!face_uv_finite(u) || !face_uv_finite(v)) {
    return Error{Errc::non_finite, "face_loops_uv: cone inverse is not finite", 0u, 0u};
  }
  return Point2d{u, v};
}

// `previous_u` and `has_previous_u` carry the last azimuth of the same loop:
// at the poles the radial component vanishes, and the previous azimuth is kept
// so that the unwrapped ring does not jump there.
inline Result<Point2d> face_uv_project_sphere(const Sphere3& surface, const Point3d& point,
                                              double previous_u, bool has_previous_u) {
  const Vec3d offset = point - surface.center;
  const double sin_v = std::clamp(dot(offset, surface.axis) / surface.radius, -1.0, 1.0);
  const double v = std::asin(sin_v);
  const Vec3d y_axis = cross(surface.axis, surface.x_axis);
  const Vec3d radial = offset - surface.axis * (surface.radius * sin_v);
  const double x = dot(radial, surface.x_axis);
  const double y = dot(radial, y_axis);
  double u = 0.0;
  if (std::hypot(x, y) > surface.radius * 1e-12) {
    u = std::atan2(y, x);
  } else if (has_previous_u) {
    u = previous_u;
  }
  if (!face_uv_finite(u) || !face_uv_finite(v)) {
    return Error{Errc::non_finite, "face_loops_uv: sphere inverse is not finite", 0u, 0u};
  }
  return Point2d{u, v};
}

inline Result<Point2d> face_uv_project_torus(const Torus3& surface, const Point3d& point) {
  const Vec3d offset = point - surface.center;
  const double axial = dot(offset, surface.axis);
  const Vec3d y_axis = cross(surface.axis, surface.x_axis);
  const double x = dot(offset, surface.x_axis);
  const double y = dot(offset, y_axis);
  if (!(std::hypot(x, y) > 0.0)) {
    return Error{Errc::invalid_topology,
                 "face_loops_uv: torus projection is singular on the axis", 0u, 0u};
  }
  const double u = std::atan2(y, x);
  // q = major_radius + minor_radius * cos(v) is the distance from the axis,
  // which is non-negative because the tube radius is smaller than the major
  // radius; axial = minor_radius * sin(v).
  const double q = std::hypot(x, y);
  const double v = std::atan2(axial, q - surface.major_radius);
  if (!face_uv_finite(u) || !face_uv_finite(v)) {
    return Error{Errc::non_finite, "face_loops_uv: torus inverse is not finite", 0u, 0u};
  }
  return Point2d{u, v};
}

inline Result<Point2d> face_uv_project(const FaceUvProjection& projection, const Point3d& point,
                                       double previous_u, bool has_previous_u) {
  if (!face_uv_point3d_finite(point)) {
    return Error{Errc::non_finite, "face_loops_uv: sample point is not finite", 0u, 0u};
  }
  if (projection.is_plane()) return face_uv_project_plane(projection.frame, point);
  const SurfaceGeometry& surface = *projection.surface;
  if (const auto* cylinder = std::get_if<Cylinder3>(&surface)) {
    return face_uv_project_cylinder(*cylinder, point);
  }
  if (const auto* cone = std::get_if<Cone3>(&surface)) return face_uv_project_cone(*cone, point);
  if (const auto* sphere = std::get_if<Sphere3>(&surface)) {
    return face_uv_project_sphere(*sphere, point, previous_u, has_previous_u);
  }
  if (const auto* torus = std::get_if<Torus3>(&surface)) {
    return face_uv_project_torus(*torus, point);
  }
  return Error{Errc::unsupported, "face_loops_uv: surface type is not supported", 0u, 0u};
}

// Evaluates the parameter pair back on the surface. The caller only passes
// parameters that came out of face_uv_project(), so the surface is analytic.
inline Point3d face_uv_eval(const FaceUvProjection& projection, const Point2d& uv) noexcept {
  if (projection.is_plane()) return plane_frame_to_3d(projection.frame, uv);
  const SurfaceGeometry& surface = *projection.surface;
  if (const auto* cylinder = std::get_if<Cylinder3>(&surface)) {
    return cylinder->eval(uv[0], uv[1]);
  }
  if (const auto* cone = std::get_if<Cone3>(&surface)) return cone->eval(uv[0], uv[1]);
  if (const auto* sphere = std::get_if<Sphere3>(&surface)) return sphere->eval(uv[0], uv[1]);
  if (const auto* torus = std::get_if<Torus3>(&surface)) return torus->eval(uv[0], uv[1]);
  return Point3d{};
}

inline Result<FaceUvProjection> face_uv_projection(const Brep& brep, const Face& face,
                                                   double tolerance) {
  FaceUvProjection projection;
  projection.tolerance = tolerance;
  const GeometryStore& geometry = brep.geometry();
  if (face.surface == kNoGeometry ||
      static_cast<std::size_t>(face.surface) >= geometry.surfaces.size()) {
    return Error{Errc::invalid_topology, "face_loops_uv: face has no surface", 0u, 0u};
  }
  projection.surface = &geometry.surfaces[face.surface];
  if (const auto* plane = std::get_if<Plane3>(projection.surface)) {
    const Result<PlaneFrame> frame = make_plane_frame(*plane);
    if (!frame.has_value()) {
      return Error{Errc::invalid_topology, "face_loops_uv: plane surface is degenerate", 0u, 0u};
    }
    projection.frame = frame.value();
    return projection;
  }
  if (std::holds_alternative<NurbsSurface>(*projection.surface)) {
    return Error{Errc::unsupported, "face_loops_uv: NURBS surfaces are not supported", 0u, 0u};
  }
  return projection;
}

// Samples an edge in its own t0 -> t1 direction. A line needs only its two
// endpoints; circles and NURBS curves are sampled uniformly with
// kFaceUvSegments segments.
inline Result<std::vector<Point3d>> face_uv_edge_samples(const Brep& brep, const Edge& edge) {
  const GeometryStore& geometry = brep.geometry();
  if (edge.curve == kNoGeometry ||
      static_cast<std::size_t>(edge.curve) >= geometry.curves.size()) {
    return Error{Errc::invalid_topology, "face_loops_uv: edge has no curve", 0u, 0u};
  }
  if (!face_uv_finite(edge.t0) || !face_uv_finite(edge.t1)) {
    return Error{Errc::non_finite, "face_loops_uv: edge parameter is not finite", 0u, 0u};
  }
  const CurveGeometry& curve = geometry.curves[edge.curve];
  std::vector<Point3d> samples;
  if (const auto* line = std::get_if<LineCurve3>(&curve)) {
    samples.reserve(2u);
    samples.push_back(line->eval(edge.t0));
    samples.push_back(line->eval(edge.t1));
  } else if (std::get_if<CircleCurve3>(&curve) != nullptr ||
             std::get_if<NurbsCurve>(&curve) != nullptr) {
    samples.reserve(kFaceUvSegments + 1u);
    for (std::size_t i = 0; i <= kFaceUvSegments; ++i) {
      const double fraction =
          static_cast<double>(i) / static_cast<double>(kFaceUvSegments);
      const double t = edge.t0 + (edge.t1 - edge.t0) * fraction;
      if (const auto* circle = std::get_if<CircleCurve3>(&curve)) {
        samples.push_back(circle->eval(t));
      } else {
        samples.push_back(std::get<NurbsCurve>(curve).eval(t));
      }
    }
  } else {
    return Error{Errc::unsupported, "face_loops_uv: curve type is not supported", 0u, 0u};
  }
  for (const Point3d& point : samples) {
    if (!face_uv_point3d_finite(point)) {
      return Error{Errc::non_finite, "face_loops_uv: edge sample is not finite", 0u, 0u};
    }
  }
  return samples;
}

// Azimuth of `point` around `circle`, normalized to [0, 2*pi). The angle is
// measured from x_axis towards y_axis and is derived from the 3D position, not
// from an edge parameter: pcurves may interpolate between the edge endpoints and
// several edges may share the same circle with different parameter ranges.
inline double face_uv_circle_angle(const CircleCurve3& circle, const Point3d& point) {
  const Vec3d offset = point - circle.center;
  double angle = std::atan2(dot(offset, circle.y_axis), dot(offset, circle.x_axis));
  if (angle < 0.0) angle += 2.0 * kFaceUvPi;
  // atan2 returns -0.0 for a point on the +x axis approached from below, and
  // adding a full turn to that lands exactly on 2*pi; keep the range half-open.
  if (angle >= 2.0 * kFaceUvPi) angle = 0.0;
  return angle;
}

// The circle a coedge's edge lies on, or nullptr when the edge is a line, a
// NURBS curve, or has no usable curve index.
inline const CircleCurve3* face_uv_coedge_circle(const Brep& brep, const Coedge& coedge) {
  if (!brep.edge_alive(coedge.edge)) return nullptr;
  const GeometryStore& geometry = brep.geometry();
  const Edge& edge = brep.edge(coedge.edge);
  if (edge.curve == kNoGeometry ||
      static_cast<std::size_t>(edge.curve) >= geometry.curves.size()) {
    return nullptr;
  }
  return std::get_if<CircleCurve3>(&geometry.curves[edge.curve]);
}

// The parameter-domain polyline of one oriented coedge, together with the 3D
// position of every point. `from_pcurve` records whether the points are the
// stored pcurve (exact parameter values) or projections of 3D samples.
struct FaceUvRun {
  std::vector<Point2d> points;
  std::vector<Point3d> positions;
  bool from_pcurve = false;
  // segments[i] describes points[i] -> points[i + 1].
  std::vector<std::optional<FaceUvArcRef>> segments;
};

// Fills `run.segments` for a coedge whose edge is a circle. Each segment is
// engaged with the signed sweep of the traversal; a straight edge leaves every
// segment disengaged. The direction comes from the chord against the tangent at
// the middle azimuth, so a concave arc (travelled clockwise about its circle)
// keeps its negative sweep.
inline Result<void> face_uv_fill_arcs(const CircleCurve3* circle, FaceUvRun& run,
                                      double tolerance) {
  run.segments.assign(run.points.size() > 0u ? run.points.size() - 1u : 0u, std::nullopt);
  if (circle == nullptr) return Result<void>{};
  const double slack = std::max(tolerance, 1e-12 * std::max(1.0, circle->radius));
  for (std::size_t i = 0; i + 1u < run.positions.size(); ++i) {
    const Point3d& start_point = run.positions[i];
    const Point3d& end_point = run.positions[i + 1u];
    // A stale pcurve can leave the circle its edge lies on, which would make the
    // provenance angle meaningless. The slack also covers the rounding of the 3D
    // position that a pcurve evaluates to.
    const double start_radial = face_uv_distance(start_point, circle->center);
    const double end_radial = face_uv_distance(end_point, circle->center);
    if (std::fabs(start_radial - circle->radius) > slack ||
        std::fabs(end_radial - circle->radius) > slack) {
      return Error{Errc::numeric_failure, "face_loops_uv: point is off its circle", 0u, 0u};
    }
    const double start = face_uv_circle_angle(*circle, start_point);
    const double end = face_uv_circle_angle(*circle, end_point);
    double sweep = end - start;
    if (sweep > kFaceUvPi) {
      sweep -= 2.0 * kFaceUvPi;
    } else if (sweep < -kFaceUvPi) {
      sweep += 2.0 * kFaceUvPi;
    }
    // Half a turn (or more) between two consecutive points has two equally short
    // arcs, so without an intermediate point the direction is unknowable.
    if (std::fabs(sweep) >= kFaceUvPi - 1e-9) {
      return Error{Errc::numeric_failure,
                   "face_loops_uv: circular arc direction is ambiguous", 0u, 0u};
    }
    const double middle = start + 0.5 * sweep;
    const Vec3d chord = end_point - start_point;
    const Vec3d tangent =
        circle->x_axis * (-std::sin(middle)) + circle->y_axis * std::cos(middle);
    // A short arc's chord is parallel to the tangent at its middle azimuth, so a
    // vanishing projection means the two directions are indistinguishable (a
    // half turn, for instance). The wrapped difference above already carries the
    // direction, so it is not flipped here.
    const double alignment = dot(chord, tangent);
    if (std::fabs(alignment) <= 1e-9 * norm2(chord)) {
      return Error{Errc::numeric_failure,
                   "face_loops_uv: circular arc direction is ambiguous", 0u, 0u};
    }
    run.segments[i] = FaceUvArcRef{*circle, start, start + sweep};
  }
  return Result<void>{};
}

inline Result<FaceUvRun> face_uv_coedge_run(const Brep& brep, const Coedge& coedge,
                                            const FaceUvProjection& projection, double tolerance,
                                            double& previous_u, bool& has_previous_u) {
  const GeometryStore& geometry = brep.geometry();
  FaceUvRun run;
  if (coedge.pcurve != kNoGeometry) {
    if (static_cast<std::size_t>(coedge.pcurve) >= geometry.pcurves.size()) {
      return Error{Errc::invalid_topology, "face_loops_uv: pcurve index is out of range", 0u, 0u};
    }
    const std::vector<Point2d>& stored = geometry.pcurves[coedge.pcurve].points;
    const std::size_t count = stored.size();
    if (count < 2u) {
      return Error{Errc::invalid_topology, "face_loops_uv: pcurve has fewer than two points", 0u,
                   0u};
    }
    run.points.reserve(count);
    run.positions.reserve(count);
    // The pcurve always runs along the edge from t0 to t1, so a reversed
    // coedge traverses it backwards.
    for (std::size_t i = 0; i < count; ++i) {
      const std::size_t index = coedge.reversed ? (count - 1u - i) : i;
      const Point2d uv = stored[index];
      if (!face_uv_point2d_finite(uv)) {
        return Error{Errc::non_finite, "face_loops_uv: pcurve point is not finite", 0u, 0u};
      }
      const Point3d position = face_uv_eval(projection, uv);
      if (!face_uv_point3d_finite(position)) {
        return Error{Errc::non_finite,
                     "face_loops_uv: pcurve point evaluates to a non-finite position", 0u, 0u};
      }
      run.points.push_back(uv);
      run.positions.push_back(position);
    }
    const Result<void> segments = face_uv_fill_arcs(face_uv_coedge_circle(brep, coedge), run, tolerance);
    if (!segments.has_value()) return segments.error();
    previous_u = run.points.back()[0];
    has_previous_u = true;
    run.from_pcurve = true;
    return run;
  }

  if (!brep.edge_alive(coedge.edge)) {
    return Error{Errc::invalid_topology, "face_loops_uv: coedge edge is not alive", 0u, 0u};
  }
  Result<std::vector<Point3d>> samples = face_uv_edge_samples(brep, brep.edge(coedge.edge));
  if (!samples.has_value()) return samples.error();
  std::vector<Point3d>& points = samples.value();
  if (coedge.reversed) std::reverse(points.begin(), points.end());
  run.points.reserve(points.size());
  run.positions.reserve(points.size());
  for (const Point3d& point : points) {
    const Result<Point2d> uv = face_uv_project(projection, point, previous_u, has_previous_u);
    if (!uv.has_value()) return uv.error();
    previous_u = uv.value()[0];
    has_previous_u = true;
    run.points.push_back(uv.value());
    run.positions.push_back(point);
  }
  const Result<void> segments = face_uv_fill_arcs(face_uv_coedge_circle(brep, coedge), run, tolerance);
  if (!segments.has_value()) return segments.error();
  return run;
}

inline Result<FaceUvLoop> face_uv_read_loop(const Brep& brep, FaceId face_id, LoopId loop_id,
                                            const Loop& loop, const FaceUvProjection& projection,
                                            double tolerance) {
  FaceUvLoop result;
  result.outer = loop.outer;

  const std::size_t coedge_limit = static_cast<std::size_t>(brep.coedge_count()) + 1u;
  double previous_u = 0.0;
  bool has_previous_u = false;
  bool chain_closed = false;
  bool started = false;
  Point3d first_position{};
  Point3d last_position{};
  double extent = 0.0;
  bool previous_from_pcurve = false;

  CoedgeId current = loop.start;
  std::size_t steps = 0u;
  while (!current.is_null()) {
    if (!brep.coedge_alive(current)) {
      return Error{Errc::invalid_topology, "face_loops_uv: loop uses a coedge that is not alive",
                   current.index(), 0u};
    }
    const Coedge& coedge = brep.coedge(current);
    if (coedge.loop != loop_id || coedge.face != face_id) {
      return Error{Errc::invalid_topology,
                   "face_loops_uv: coedge belongs to a different loop or face",
                   current.index(), 0u};
    }
    const Result<FaceUvRun> run_result = face_uv_coedge_run(brep, coedge, projection, tolerance,
                                                            previous_u, has_previous_u);
    if (!run_result.has_value()) return run_result.error();
    const FaceUvRun& run = run_result.value();

    for (std::size_t i = 0; i < run.points.size(); ++i) {
      Point2d uv = run.points[i];
      if (!face_uv_point2d_finite(uv) || !face_uv_point3d_finite(run.positions[i])) {
        return Error{Errc::non_finite, "face_loops_uv: loop point is not finite", 0u, 0u};
      }
      if (started) {
        const Point2d& previous = result.points.back();
        if (projection.periodic_u()) uv[0] = face_uv_unwrap(uv[0], previous[0]);
        if (projection.periodic_v()) uv[1] = face_uv_unwrap(uv[1], previous[1]);
        if (i == 0u) {
          // One point per junction: the coedges meet either at exactly the same
          // parameter value (stored pcurves) or at the same 3D vertex. The
          // segment leaving the junction is appended with the next point, so
          // dropping the duplicate point keeps the segment refs aligned.
          const bool uv_junction = previous_from_pcurve && run.from_pcurve &&
                                   uv[0] == previous[0] && uv[1] == previous[1];
          const bool position_junction =
              face_uv_distance(run.positions[i], last_position) <= tolerance;
          if (uv_junction || position_junction) continue;
          return Error{Errc::invalid_topology, "face_loops_uv: coedge chain is not continuous",
                       current.index(), 0u};
        }
      }
      result.points.push_back(uv);
      // The segment ending at this point starts at the previous loop point.
      if (i > 0u) result.arcs.push_back(run.segments[i - 1u]);
      if (!started) {
        first_position = run.positions[i];
        started = true;
      }
      last_position = run.positions[i];
      const double reach = face_uv_distance(run.positions[i], first_position);
      if (reach > extent) extent = reach;
      previous_from_pcurve = run.from_pcurve;
    }

    if (coedge.next == loop.start) {
      chain_closed = true;
      break;
    }
    if (++steps >= coedge_limit) {
      return Error{Errc::invalid_topology, "face_loops_uv: loop coedge cycle does not close",
                   current.index(), 0u};
    }
    current = coedge.next;
  }
  if (!chain_closed) {
    return Error{Errc::invalid_topology,
                 "face_loops_uv: loop coedge cycle does not return to its start", loop_id.index(),
                 0u};
  }
  if (result.points.empty()) {
    return Error{Errc::invalid_topology, "face_loops_uv: loop has no points", loop_id.index(), 0u};
  }
  // Closure is judged in 3D: parameter periods are angles, not lengths, so a
  // uv comparison against `tolerance` would be meaningless on a cylinder and
  // would reject the torus seam, where u wraps by one full turn.
  const double closure = std::max(tolerance, 1e-12 * extent);
  if (face_uv_distance(last_position, first_position) > closure) {
    return Error{Errc::numeric_failure, "face_loops_uv: loop is not closed", loop_id.index(), 0u};
  }
  // The last point repeats the first one; dropping it leaves the closing
  // segment in `arcs`, which is the segment from the last point back to the
  // first.
  result.points.pop_back();
  if (result.arcs.size() != result.points.size()) {
    return Error{Errc::invalid_topology, "face_loops_uv: arc provenance size mismatch",
                 loop_id.index(), 0u};
  }
  if (result.points.size() < 3u) {
    return Error{Errc::invalid_topology, "face_loops_uv: loop has fewer than three points",
                 loop_id.index(), 0u};
  }
  return result;
}

}  // namespace detail

// Reads every boundary loop of `face` in the parameter domain of the face's
// surface. Consecutive parameters are unwrapped to differ by less than half a
// turn, and the closing point is not repeated. A ring that wraps a full period
// (for example the base circle of a cone) is legal: its last-to-first step is
// then 2*pi rather than a small one. The outer loop is reported first.
//
// A coedge uses its stored pcurve when it has one and otherwise samples the
// edge in 3D and projects the samples with the surface inverse. Parametric
// surfaces (cylinder, cone, sphere, torus) are supported, planes use their own
// frame, and NURBS surfaces are rejected as unsupported.
//
// Errors: a negative or non-finite tolerance is Errc::tolerance_violation, a
// dead face or a face without loops is Errc::invalid_topology, a surface or
// curve type that cannot be parameterized here is Errc::unsupported, non-finite
// pcurve points or samples are Errc::non_finite, and a chain that does not close
// is Errc::numeric_failure. Every message begins with "face_loops_uv:".
inline Result<std::vector<FaceUvLoop>> face_loops_uv(const Brep& brep, FaceId face_id,
                                                     double tolerance) {
  if (!detail::face_uv_finite(tolerance) || !(tolerance >= 0.0)) {
    return Error{Errc::tolerance_violation,
                 "face_loops_uv: tolerance must be finite and non-negative", 0u, 0u};
  }
  if (!brep.face_alive(face_id)) {
    return Error{Errc::invalid_topology, "face_loops_uv: face is not alive", face_id.index(), 0u};
  }
  const Result<detail::FaceUvProjection> projection =
      detail::face_uv_projection(brep, brep.face(face_id), tolerance);
  if (!projection.has_value()) return projection.error();

  std::vector<FaceUvLoop> loops;
  const std::size_t loop_limit = static_cast<std::size_t>(brep.loop_count()) + 1u;
  LoopId current = brep.face(face_id).first_loop;
  std::size_t steps = 0u;
  while (!current.is_null()) {
    if (!brep.loop_alive(current)) {
      return Error{Errc::invalid_topology, "face_loops_uv: face uses a loop that is not alive",
                   current.index(), 0u};
    }
    const Loop& loop = brep.loop(current);
    if (loop.face != face_id) {
      return Error{Errc::invalid_topology, "face_loops_uv: loop belongs to a different face",
                   current.index(), 0u};
    }
    Result<FaceUvLoop> loop_result =
        detail::face_uv_read_loop(brep, face_id, current, loop, projection.value(), tolerance);
    if (!loop_result.has_value()) return loop_result.error();
    loops.push_back(std::move(loop_result.value()));
    if (++steps >= loop_limit) {
      return Error{Errc::invalid_topology, "face_loops_uv: loop chain does not terminate",
                   current.index(), 0u};
    }
    current = loop.next;
  }
  if (loops.empty()) {
    return Error{Errc::invalid_topology, "face_loops_uv: face has no loops", face_id.index(), 0u};
  }

  // Outer loop first, remaining loops in their stored order.
  std::vector<FaceUvLoop> ordered;
  ordered.reserve(loops.size());
  for (const FaceUvLoop& loop : loops) {
    if (loop.outer) ordered.push_back(loop);
  }
  for (const FaceUvLoop& loop : loops) {
    if (!loop.outer) ordered.push_back(loop);
  }
  return ordered;
}

}  // namespace v1
}  // namespace geom
