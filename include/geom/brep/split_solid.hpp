#pragma once

// M3b-2 Task 6: split a solid by a plane into two trimmed halves.
//
// The driver trims every face against the cut plane with the curved
// surface-surface intersector, heals the circular cut chains into exact
// `CircleCurve3` arc spans, stitches the cut polylines into planar cap loops and
// assembles one closed shell per side. Both sides are reported in a
// `SplitSolidResult`: `below` is the negative side of the plane normal and
// `above` the positive side.
//
// The exact circles used for healing are the intersections of the cut plane with
// the analytic surfaces of the solid: a cylinder whose axis is parallel to the
// plane normal and a sphere that meets the plane in a real circle. Both the
// trimmed lateral fragments and the two cap faces are healed against the same
// circle objects, so their shared cut edges are bit-identical arcs.

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <numeric>
#include <optional>
#include <span>
#include <utility>
#include <vector>

#include <geom/brep/assemble_curved.hpp>
#include <geom/brep/detail/arc_heal.hpp>
#include <geom/brep/detail/ssi_patch.hpp>
#include <geom/brep/detail/ssi_surface.hpp>
#include <geom/brep/face_uv.hpp>
#include <geom/brep/plane_frame.hpp>
#include <geom/brep/ssi.hpp>
#include <geom/brep/topology.hpp>
#include <geom/brep/trim_face.hpp>
#include <geom/cgeom/polygon2.hpp>
#include <geom/cgeom/triangulate2.hpp>
#include <geom/core/error.hpp>
#include <geom/core/result.hpp>
#include <geom/curves/analytic_curve.hpp>
#include <geom/surfaces/analytic.hpp>
#include <geom/surfaces/plane.hpp>

namespace geom {
inline namespace v1 {

struct SplitSolidResult {
  Brep below;  // negative side of the plane normal
  Brep above;  // positive side of the plane normal
};

namespace detail {

inline constexpr double kSplitSolidPi = 3.141592653589793238462643383279502884;
inline constexpr double kSplitSolidTwoPi = 2.0 * kSplitSolidPi;

inline Vec3d split_solid_unit(const Vec3d& v) noexcept {
  const double length = norm2(v);
  if (!(length > 0.0) || !std::isfinite(length)) return Vec3d{0.0, 0.0, 0.0};
  return v * (1.0 / length);
}

inline double split_solid_wrapped_delta(double delta) noexcept {
  delta = std::fmod(delta, kSplitSolidTwoPi);
  if (delta <= -kSplitSolidPi) delta += kSplitSolidTwoPi;
  if (delta > kSplitSolidPi) delta -= kSplitSolidTwoPi;
  return delta;
}

inline Vec3d split_solid_circle_normal(const CircleCurve3& circle) noexcept {
  return cross(circle.x_axis, circle.y_axis);
}

inline double split_solid_circle_angle(const CircleCurve3& circle,
                                       const Point3d& point) noexcept {
  const Vec3d offset = point - circle.center;
  return std::atan2(dot(offset, circle.y_axis), dot(offset, circle.x_axis));
}

inline bool split_solid_near(const Point3d& a, const Point3d& b, double tolerance) noexcept {
  return norm2(a - b) <= tolerance;
}

inline Vec3d split_solid_as_vec(const Point3d& point) noexcept {
  return Vec3d{point[0], point[1], point[2]};
}

inline bool split_solid_same2(const Point2d& a, const Point2d& b, double tolerance) noexcept {
  return std::abs(a[0] - b[0]) <= tolerance && std::abs(a[1] - b[1]) <= tolerance;
}

inline bool split_solid_on_circle(const CircleCurve3& circle, const Point3d& point,
                                  double tolerance) noexcept {
  const double slack = std::max(10.0 * tolerance, 1e-9 * circle.radius);
  const Vec3d normal = split_solid_circle_normal(circle);
  const Vec3d offset = point - circle.center;
  const double axial = dot(offset, normal);
  if (std::abs(axial) > slack) return false;
  const Vec3d in_plane = offset - normal * axial;
  return std::abs(norm2(in_plane) - circle.radius) <= slack;
}

// Signed area of a closed uv ring: positive when the ring is counter-clockwise.
inline double split_solid_ring_area(const std::vector<Point2d>& ring) noexcept {
  double twice_area = 0.0;
  const std::size_t count = ring.size();
  for (std::size_t i = 0; i < count; ++i) {
    const Point2d& a = ring[i];
    const Point2d& b = ring[(i + 1) % count];
    twice_area += a[0] * b[1] - b[0] * a[1];
  }
  return 0.5 * twice_area;
}

// Even-odd ray cast; the point is assumed not to lie exactly on the ring.
inline bool split_solid_point_in_ring(const std::vector<Point2d>& ring,
                                      const Point2d& point) noexcept {
  bool inside = false;
  const std::size_t count = ring.size();
  for (std::size_t i = 0, j = count - 1; i < count; j = i++) {
    const Point2d& a = ring[i];
    const Point2d& b = ring[j];
    if (((a[1] > point[1]) != (b[1] > point[1])) &&
        (point[0] < (b[0] - a[0]) * (point[1] - a[1]) / (b[1] - a[1]) + a[0])) {
      inside = !inside;
    }
  }
  return inside;
}

// Reverses a loop: point order, and the per-segment spans so span i still
// describes the segment from point i to point (i + 1) % size.
inline std::vector<std::optional<CurvedArcSpan>> split_solid_reversed_spans(
    const std::vector<std::optional<CurvedArcSpan>>& spans) {
  const std::size_t count = spans.size();
  std::vector<std::optional<CurvedArcSpan>> out(count);
  for (std::size_t i = 0; i < count; ++i) {
    out[i] = spans[(2 * count + count - 2 - i) % count];
  }
  return out;
}

inline void split_solid_reverse_loop(std::vector<Point3d>& points, std::vector<Point2d>& uv,
                                     std::vector<std::optional<CurvedArcSpan>>& spans) {
  const std::size_t count = points.size();
  std::reverse(points.begin(), points.end());
  if (uv.size() == count) std::reverse(uv.begin(), uv.end());
  if (spans.size() == count) spans = split_solid_reversed_spans(spans);
}

// Reverses a closed loop while keeping its first point first: [p0, p_{n-1}, ...,
// p1]. Cap loops must keep their start vertex because a full-circle cut rim is
// keyed by its single welded vertex: if the cap chain and the matching trimmed
// fragment started at different welded vertices the assembler would see the arc
// edge used only once. New segment i runs p_{n-i} -> p_{n-1-i}, which is old
// segment (n - 1 - i) mod n.
inline void split_solid_reverse_loop_keep_start(
    std::vector<Point3d>& points, std::vector<Point2d>& uv,
    std::vector<std::optional<CurvedArcSpan>>& spans) {
  const std::size_t count = points.size();
  if (count < 2) return;
  std::vector<Point3d> rotated_points;
  rotated_points.reserve(count);
  rotated_points.push_back(points.front());
  for (std::size_t i = count; i-- > 1;) rotated_points.push_back(points[i]);
  points = std::move(rotated_points);
  if (uv.size() == count) {
    std::vector<Point2d> rotated_uv;
    rotated_uv.reserve(count);
    rotated_uv.push_back(uv.front());
    for (std::size_t i = count; i-- > 1;) rotated_uv.push_back(uv[i]);
    uv = std::move(rotated_uv);
  }
  if (spans.size() == count) {
    std::vector<std::optional<CurvedArcSpan>> rotated_spans(count);
    const std::size_t n = count;
    for (std::size_t i = 0; i < count; ++i) rotated_spans[i] = spans[(n - 1 - i) % n];
    spans = std::move(rotated_spans);
  }
}

// A single loop in surface parameter space plus its healed arc spans.
struct SplitSolidLoop {
  std::vector<Point3d> points;
  std::vector<Point2d> uv;
  std::vector<std::optional<CurvedArcSpan>> spans;
};

// Engages spans on every on-circle run of one loop.
//
// Runs shorter than half a turn are left unengaged for
// `collapse_circular_chains`, which merges them into a single exact span. Runs
// of half a turn or more are engaged segment by segment, because
// `assemble_curved_solid` recovers a span's traversal direction from the point
// order with a wrapped delta in (-pi, pi]: a span wider than pi can never match
// that test, and two spans that split a full turn cannot both be shorter than
// half a turn. Per-segment spans are exactly the representation the assembler
// receives from `face_loops_uv` for the full-circle rims of `make_cylinder`.
inline Result<void> split_solid_heal_loop(SplitSolidLoop& loop, const CircleCurve3& circle,
                                          double tolerance) {
  const std::size_t count = loop.points.size();
  if (count < 3) return Result<void>{};
  if (loop.uv.size() != count) return Result<void>{};
  if (!circle.is_valid()) return Result<void>{};
  if (loop.spans.size() != count) loop.spans.assign(count, std::nullopt);
  const double angle_tolerance =
      std::min(kSplitSolidPi, std::max(1e-12, 2.0 * tolerance / circle.radius));

  std::vector<char> on(count, 0);
  bool all_on = true;
  for (std::size_t i = 0; i < count; ++i) {
    on[i] = split_solid_on_circle(circle, loop.points[i], tolerance) ? 1 : 0;
    if (!on[i]) all_on = false;
  }

  std::vector<char> usable(count, 0);
  std::size_t usable_count = 0;
  for (std::size_t i = 0; i < count; ++i) {
    const std::size_t next = (i + 1) % count;
    if (loop.spans[i].has_value() || !on[i] || !on[next]) continue;
    usable[i] = 1;
    ++usable_count;
  }

  auto heal_run = [&](std::size_t run_start, std::size_t length) {
    double sweep = 0.0;
    for (std::size_t k = 0; k < length; ++k) {
      const std::size_t a = (run_start + k) % count;
      const std::size_t b = (a + 1) % count;
      sweep += split_solid_wrapped_delta(split_solid_circle_angle(circle, loop.points[b]) -
                                         split_solid_circle_angle(circle, loop.points[a]));
    }
    if (std::abs(sweep) < kSplitSolidPi - angle_tolerance) return;
    for (std::size_t k = 0; k < length; ++k) {
      const std::size_t a = (run_start + k) % count;
      const std::size_t b = (a + 1) % count;
      const double t0 = split_solid_circle_angle(circle, loop.points[a]);
      const double delta = split_solid_wrapped_delta(
          split_solid_circle_angle(circle, loop.points[b]) - t0);
      // A half turn cannot record a direction; leave that segment alone.
      if (!(std::abs(delta) < kSplitSolidPi)) continue;
      const double t1 = t0 + delta;
      loop.spans[a] = CurvedArcSpan{circle, std::min(t0, t1), std::max(t0, t1)};
    }
  };

  if (usable_count > 0) {
    std::size_t start = count;
    for (std::size_t i = 0; i < count; ++i) {
      if (usable[i] && !usable[(i + count - 1) % count]) {
        start = i;
        break;
      }
    }
    if (start == count) {
      heal_run(0, count);
    } else {
      std::size_t index = start;
      std::size_t walked = 0;
      while (walked < count) {
        while (walked < count && !usable[index]) {
          index = (index + 1) % count;
          ++walked;
        }
        if (walked >= count) break;
        const std::size_t run_start = index;
        std::size_t length = 0;
        while (walked < count && usable[index]) {
          index = (index + 1) % count;
          ++length;
          ++walked;
        }
        if (length >= 2) heal_run(run_start, length);
      }
    }
  }

  bool unengaged_on = false;
  for (std::size_t i = 0; i < count; ++i) {
    if (loop.spans[i].has_value() || !on[i] || !on[(i + 1) % count]) continue;
    unengaged_on = true;
    break;
  }
  // A loop that lies entirely on the circle has already been handled above:
  // handing it to `collapse_circular_chains` would rewrite it into two
  // half-turn spans and drop any span the pre-pass just engaged.
  if (unengaged_on && !all_on) {
    CurvedFragment fragment;
    fragment.loops.push_back(loop.points);
    fragment.loops_uv.push_back(loop.uv);
    fragment.arc_spans.push_back(loop.spans);
    fragment.reversed = false;
    const Result<bool> healed = collapse_circular_chains(fragment, circle, tolerance);
    if (!healed.has_value()) return healed.error();
    if (healed.value() && fragment.loops.size() == 1 && fragment.loops_uv.size() == 1 &&
        fragment.arc_spans.size() == 1) {
      const std::size_t healed_count = fragment.loops[0].size();
      if (fragment.loops_uv[0].size() == healed_count &&
          fragment.arc_spans[0].size() == healed_count) {
        loop.points = std::move(fragment.loops[0]);
        loop.uv = std::move(fragment.loops_uv[0]);
        loop.spans = std::move(fragment.arc_spans[0]);
      }
    }
  }

  for (std::optional<CurvedArcSpan>& span : loop.spans) {
    if (span.has_value() && span->t1 < span->t0) std::swap(span->t0, span->t1);
  }
  return Result<void>{};
}

inline Result<void> split_solid_heal_fragment(CurvedFragment& fragment,
                                              const std::vector<CircleCurve3>& circles,
                                              double tolerance) {
  if (circles.empty()) return Result<void>{};
  if (fragment.loops_uv.size() != fragment.loops.size()) return Result<void>{};
  if (fragment.arc_spans.size() != fragment.loops.size()) {
    fragment.arc_spans.resize(fragment.loops.size());
  }
  for (std::size_t l = 0; l < fragment.loops.size(); ++l) {
    SplitSolidLoop loop;
    loop.points = fragment.loops[l];
    loop.uv = fragment.loops_uv[l];
    const std::size_t count = loop.points.size();
    if (loop.uv.size() != count) continue;
    if (fragment.arc_spans[l].size() == count) {
      loop.spans = fragment.arc_spans[l];
    } else {
      loop.spans.assign(count, std::nullopt);
    }
    for (const CircleCurve3& circle : circles) {
      const Result<void> healed = split_solid_heal_loop(loop, circle, tolerance);
      if (!healed.has_value()) return healed.error();
    }
    fragment.loops[l] = std::move(loop.points);
    fragment.loops_uv[l] = std::move(loop.uv);
    fragment.arc_spans[l] = std::move(loop.spans);
  }
  return Result<void>{};
}

// A loop whose every segment is an arc span on the same circle is assembled as
// a single edge: its one welded vertex becomes the edge's start key. The two
// uses of such an edge (for example a cylinder rim shared by a cap and the
// lateral face) must therefore start at the SAME 3D point. A trimmed lateral
// run always starts at the face's u = 0 column, i.e. at angle 0 of the shared
// edge circle, so every all-arc loop is rotated to begin at its smallest
// non-negative angle around the circle. Both uses read the circle from the same
// edge geometry (or from the same constructed cut circle), so the choice is
// consistent. Rotation cannot split a run because the whole loop is one run.
inline void split_solid_canonicalize_all_arc_loop(std::vector<Point3d>& points,
                                                  std::vector<Point2d>& uv,
                                                  std::vector<std::optional<CurvedArcSpan>>& spans,
                                                  double tolerance) {
  const std::size_t count = points.size();
  if (count < 3 || spans.size() != count || uv.size() != count) return;
  if (!spans.front().has_value()) return;
  const CircleCurve3 circle = spans.front()->circle;
  if (!circle.is_valid()) return;
  const Vec3d axis = split_solid_circle_normal(circle);
  const double slack = std::max(10.0 * tolerance, 1e-9 * circle.radius);
  for (const std::optional<CurvedArcSpan>& span : spans) {
    if (!span.has_value()) return;
    const Vec3d other_axis = split_solid_circle_normal(span->circle);
    if (norm2(other_axis - axis) > 1e-9 && norm2(other_axis + axis) > 1e-9) return;
    if (std::abs(span->circle.radius - circle.radius) > slack) return;
    if (norm2(span->circle.center - circle.center) > slack) return;
  }
  for (const Point3d& point : points) {
    if (!split_solid_on_circle(circle, point, tolerance)) return;
  }
  std::size_t start = 0;
  double best_angle = std::numeric_limits<double>::infinity();
  for (std::size_t i = 0; i < count; ++i) {
    double angle = split_solid_circle_angle(circle, points[i]);
    angle -= kSplitSolidTwoPi * std::floor(angle / kSplitSolidTwoPi);
    // A point sampled at t = 2*pi evaluates with a tiny negative sine, which
    // wraps to a full turn; treat that as the zero angle it represents.
    if (angle > kSplitSolidTwoPi - 1e-9) angle = 0.0;
    if (angle < best_angle - 1e-12) {
      best_angle = angle;
      start = i;
    }
  }
  if (start == 0) return;
  std::vector<Point3d> rotated_points(count);
  std::vector<Point2d> rotated_uv(count);
  std::vector<std::optional<CurvedArcSpan>> rotated_spans(count);
  for (std::size_t i = 0; i < count; ++i) {
    const std::size_t source = (start + i) % count;
    rotated_points[i] = points[source];
    rotated_uv[i] = uv[source];
    rotated_spans[i] = spans[source];
  }
  points = std::move(rotated_points);
  uv = std::move(rotated_uv);
  spans = std::move(rotated_spans);
}

inline void split_solid_canonicalize_face(CurvedAssembledFace& face, double tolerance) {
  for (std::size_t l = 0; l < face.loops.size(); ++l) {
    if (l >= face.arc_spans.size() || l >= face.loops_uv.size()) continue;
    split_solid_canonicalize_all_arc_loop(face.loops[l], face.loops_uv[l], face.arc_spans[l],
                                          tolerance);
  }
}

// The exact circle in which `surface` meets the plane n . x = offset, when that
// intersection is a circle known analytically. The circle axis is +n for every
// face so the trimmed fragments and both caps share one circle object.
inline bool split_solid_face_circle(const SurfaceGeometry& surface, const Vec3d& normal,
                                    double offset, const PlaneFrame& frame,
                                    CircleCurve3& circle) noexcept {
  if (const auto* cylinder = std::get_if<Cylinder3>(&surface)) {
    const double alignment = dot(normal, cylinder->axis);
    if (!(std::abs(alignment) >= 1.0 - 1e-12)) return false;
    if (!(cylinder->radius > 0.0)) return false;
    const double t = (offset - dot(normal, split_solid_as_vec(cylinder->origin))) / alignment;
    if (!std::isfinite(t)) return false;
    Point3d center = cylinder->origin + cylinder->axis * t;
    center = center - normal * (dot(normal, split_solid_as_vec(center)) - offset);
    Vec3d x_axis = split_solid_unit(cylinder->x_axis - normal * dot(cylinder->x_axis, normal));
    if (!(norm2(x_axis) > 0.5)) x_axis = frame.u;
    const Vec3d y_axis = split_solid_unit(cross(normal, x_axis));
    circle = CircleCurve3{center, x_axis, y_axis, cylinder->radius};
    return circle.is_valid();
  }
  if (const auto* sphere = std::get_if<Sphere3>(&surface)) {
    if (!(sphere->radius > 0.0)) return false;
    const double distance = dot(normal, split_solid_as_vec(sphere->center)) - offset;
    if (!(std::abs(distance) < sphere->radius)) return false;
    const double radius2 = sphere->radius * sphere->radius - distance * distance;
    if (!(radius2 > 0.0)) return false;
    const Point3d center = sphere->center - normal * distance;
    Vec3d x_axis = split_solid_unit(sphere->x_axis - normal * dot(sphere->x_axis, normal));
    if (!(norm2(x_axis) > 0.5)) x_axis = frame.u;
    const Vec3d y_axis = split_solid_unit(cross(normal, x_axis));
    circle = CircleCurve3{center, x_axis, y_axis, std::sqrt(radius2)};
    return circle.is_valid();
  }
  return false;
}

// Liang-Barsky clipping of one uv polyline to an axis-aligned rectangle. The
// pieces that stay inside are chained; a segment that leaves the rectangle
// flushes the current chain.
inline std::vector<std::vector<Point2d>> split_solid_clip_polyline(
    const std::vector<Point2d>& polyline, double u0, double u1, double v0, double v1,
    double weld) {
  std::vector<std::vector<Point2d>> chains;
  if (polyline.size() < 2) return chains;
  std::vector<Point2d> current;
  for (std::size_t i = 0; i + 1 < polyline.size(); ++i) {
    const Point2d a = polyline[i];
    const Point2d b = polyline[i + 1];
    const double du = b[0] - a[0];
    const double dv = b[1] - a[1];
    double t0 = 0.0;
    double t1 = 1.0;
    bool accepted = true;
    auto clip = [&t0, &t1](double p, double q) {
      if (p == 0.0) return q >= 0.0;
      const double r = q / p;
      if (p < 0.0) {
        if (r > t1) return false;
        if (r > t0) t0 = r;
      } else {
        if (r < t0) return false;
        if (r < t1) t1 = r;
      }
      return true;
    };
    accepted = clip(-du, a[0] - u0) && clip(du, u1 - a[0]) && clip(-dv, a[1] - v0) &&
               clip(dv, v1 - a[1]);
    if (!accepted || t1 < t0) {
      if (current.size() >= 2) chains.push_back(std::move(current));
      current.clear();
      continue;
    }
    const Point2d start{a[0] + t0 * du, a[1] + t0 * dv};
    const Point2d end{a[0] + t1 * du, a[1] + t1 * dv};
    if (current.empty() || !split_solid_same2(current.back(), start, weld)) {
      if (current.size() >= 2) chains.push_back(std::move(current));
      current.clear();
      current.push_back(start);
    }
    if (!split_solid_same2(start, end, weld)) current.push_back(end);
  }
  if (current.size() >= 2) chains.push_back(std::move(current));
  return chains;
}

// Solves one free parameter so a point already pinned to a face boundary
// column/row also lies exactly on the cut plane.
inline void split_solid_solve_free(const SurfaceGeometry& surface, const Vec3d& normal,
                                   double offset, Point2d& point, bool fixed_u) {
  for (int iteration = 0; iteration < 4; ++iteration) {
    const Result<SurfacePoint> evaluated = surface_eval(surface, point[0], point[1]);
    if (!evaluated.has_value()) return;
    const double gap = dot(normal, split_solid_as_vec(evaluated.value().p)) - offset;
    if (std::abs(gap) <= 1e-15 * (1.0 + std::abs(offset))) return;
    const double slope =
        fixed_u ? dot(normal, evaluated.value().dv) : dot(normal, evaluated.value().du);
    if (!(std::abs(slope) > 0.0)) return;
    if (fixed_u) {
      point[1] -= gap / slope;
    } else {
      point[0] -= gap / slope;
    }
  }
}

// The intersector can end a whole tracing step short of a patch boundary when
// its last seed falls outside the patch domain. Pull such an endpoint onto the
// face uv rectangle along the last chord, then re-project it onto the cut plane
// so the cut reaches the face boundary exactly.
inline void split_solid_extend_cut_end(std::vector<Point2d>& polyline, bool front,
                                       const SurfaceGeometry& surface, const Vec3d& normal,
                                       double offset, double u0, double u1, double v0, double v1,
                                       double snap, double limit) {
  if (polyline.size() < 2) return;
  const std::size_t end = front ? 0 : polyline.size() - 1;
  const std::size_t inner = front ? 1 : polyline.size() - 2;
  const Point2d point = polyline[end];
  const Point2d& other = polyline[inner];
  const double near_u = std::min(std::abs(point[0] - u0), std::abs(point[0] - u1));
  const double near_v = std::min(std::abs(point[1] - v0), std::abs(point[1] - v1));
  if (std::min(near_u, near_v) <= snap) return;
  if (point[0] < u0 || point[0] > u1 || point[1] < v0 || point[1] > v1) return;
  const double du = point[0] - other[0];
  const double dv = point[1] - other[1];
  const double length = std::sqrt(du * du + dv * dv);
  const double scale = std::max({1.0, std::abs(u1 - u0), std::abs(v1 - v0)});
  const double epsilon = 1e-15 * scale;
  if (!(length > epsilon)) return;
  const double ux = du / length;
  const double vx = dv / length;
  double best_t = limit;
  int best_kind = -1;
  double best_u = 0.0;
  double best_v = 0.0;
  auto consider = [&](double t, double u, double v, int kind) {
    if (!(t > 0.0) || t > best_t) return;
    if (u < u0 - epsilon || u > u1 + epsilon || v < v0 - epsilon || v > v1 + epsilon) return;
    best_t = t;
    best_kind = kind;
    best_u = u;
    best_v = v;
  };
  if (std::abs(ux) > 1e-12) {
    consider((u0 - point[0]) / ux, u0, point[1] + (u0 - point[0]) / ux * vx, 0);
    consider((u1 - point[0]) / ux, u1, point[1] + (u1 - point[0]) / ux * vx, 0);
  }
  if (std::abs(vx) > 1e-12) {
    consider((v0 - point[1]) / vx, point[0] + (v0 - point[1]) / vx * ux, v0, 1);
    consider((v1 - point[1]) / vx, point[0] + (v1 - point[1]) / vx * ux, v1, 1);
  }
  if (best_kind < 0) return;
  Point2d extended{best_u, best_v};
  split_solid_solve_free(surface, normal, offset, extended, best_kind == 0);
  polyline[end] = extended;
}

struct SplitSolidFaceInfo {
  FaceId face{};
  const SurfaceGeometry* surface = nullptr;
  FaceUvProjection projection{};
  std::vector<FaceUvLoop> loops;
  double distance_min = 0.0;
  double distance_max = 0.0;
};

struct SplitSolidCutChain {
  std::vector<Point3d> points;
  std::optional<CircleCurve3> circle;
};

struct SplitSolidCutLoop {
  std::vector<Point3d> points;
  std::vector<std::optional<CurvedArcSpan>> spans;
  std::vector<CircleCurve3> circles;
};

inline void split_solid_append_circles(std::vector<CircleCurve3>& circles,
                                       const CircleCurve3& circle) {
  for (const CircleCurve3& existing : circles) {
    if (split_solid_near(existing.center, circle.center, 1e-9 * std::max(1.0, circle.radius)) &&
        std::abs(existing.radius - circle.radius) <= 1e-9 * std::max(1.0, circle.radius)) {
      return;
    }
  }
  circles.push_back(circle);
}

// Drops chains that trace the same intersection curve twice. The tracer can
// report two coincident curves for one face when several seeds land on it; the
// duplicate cut would otherwise make the cap stitch fail.
inline void split_solid_dedupe_chains(std::vector<SplitSolidCutChain>& chains, double weld) {
  std::vector<SplitSolidCutChain> unique;
  unique.reserve(chains.size());
  for (SplitSolidCutChain& chain : chains) {
    bool duplicate = false;
    for (const SplitSolidCutChain& kept : unique) {
      if (kept.points.size() != chain.points.size()) continue;
      bool same = true;
      bool mirrored = true;
      for (std::size_t i = 0; i < kept.points.size(); ++i) {
        if (!split_solid_near(kept.points[i], chain.points[i], weld)) same = false;
        if (!split_solid_near(kept.points[i], chain.points[chain.points.size() - 1 - i], weld)) {
          mirrored = false;
        }
      }
      if (same || mirrored) {
        duplicate = true;
        break;
      }
    }
    if (!duplicate) unique.push_back(std::move(chain));
  }
  chains = std::move(unique);
}

// Drops consecutive points that collapse onto one another. Extending a cut end
// onto a face boundary can leave a sub-tolerance backward kink against the
// neighbouring trace point; the arrangement would keep the resulting sliver as
// an extra fragment sitting exactly on the cut plane.
inline void split_solid_dedupe_polyline(std::vector<Point2d>& points, double weld) {
  std::vector<Point2d> unique;
  unique.reserve(points.size());
  for (const Point2d& point : points) {
    if (!unique.empty() && std::hypot(point[0] - unique.back()[0], point[1] - unique.back()[1]) <=
                               weld) {
      continue;
    }
    unique.push_back(point);
  }
  points = std::move(unique);
}

// Stitches the cut polylines of all faces into closed loops in the cut plane.
inline Result<std::vector<SplitSolidCutLoop>> split_solid_stitch_cuts(
    std::vector<SplitSolidCutChain>& chains, double weld) {
  std::vector<SplitSolidCutLoop> loops;
  std::vector<char> used(chains.size(), 0);
  for (std::size_t i = 0; i < chains.size(); ++i) {
    if (used[i]) continue;
    std::vector<Point3d> points = chains[i].points;
    std::vector<CircleCurve3> circles;
    if (chains[i].circle.has_value()) circles.push_back(*chains[i].circle);
    used[i] = 1;
    if (points.size() < 3) continue;
    if (split_solid_near(points.front(), points.back(), weld)) {
      points.pop_back();
      if (points.size() < 3) continue;
      SplitSolidCutLoop loop;
      loop.points = std::move(points);
      loop.spans.assign(loop.points.size(), std::nullopt);
      loop.circles = std::move(circles);
      loops.push_back(std::move(loop));
      continue;
    }
    bool closed = false;
    while (!closed) {
      bool extended = false;
      for (std::size_t j = 0; j < chains.size(); ++j) {
        if (used[j]) continue;
        const std::vector<Point3d>& other = chains[j].points;
        if (other.size() < 2) {
          used[j] = 1;
          continue;
        }
        if (split_solid_near(other.front(), points.back(), weld)) {
          points.insert(points.end(), other.begin() + 1, other.end());
          if (chains[j].circle.has_value()) {
            split_solid_append_circles(circles, *chains[j].circle);
          }
          used[j] = 1;
          extended = true;
          break;
        }
        if (split_solid_near(other.back(), points.back(), weld)) {
          points.insert(points.end(), other.rbegin() + 1, other.rend());
          if (chains[j].circle.has_value()) {
            split_solid_append_circles(circles, *chains[j].circle);
          }
          used[j] = 1;
          extended = true;
          break;
        }
        if (split_solid_near(other.front(), points.front(), weld)) {
          std::vector<Point3d> merged(other.rbegin(), other.rend());
          merged.insert(merged.end(), points.begin() + 1, points.end());
          points = std::move(merged);
          if (chains[j].circle.has_value()) {
            split_solid_append_circles(circles, *chains[j].circle);
          }
          used[j] = 1;
          extended = true;
          break;
        }
        if (split_solid_near(other.back(), points.front(), weld)) {
          std::vector<Point3d> merged(other.begin(), other.end());
          merged.insert(merged.end(), points.begin() + 1, points.end());
          points = std::move(merged);
          if (chains[j].circle.has_value()) {
            split_solid_append_circles(circles, *chains[j].circle);
          }
          used[j] = 1;
          extended = true;
          break;
        }
      }
      if (!extended) {
        return Error{Errc::numeric_failure, "split_solid_by_plane: cut loops do not close",
                     static_cast<std::uint32_t>(i), 0u};
      }
      closed = split_solid_near(points.front(), points.back(), weld);
    }
    if (points.size() < 4) {
      return Error{Errc::numeric_failure, "split_solid_by_plane: cut loops do not close",
                   static_cast<std::uint32_t>(i), 0u};
    }
    points.pop_back();
    SplitSolidCutLoop loop;
    loop.points = std::move(points);
    loop.spans.assign(loop.points.size(), std::nullopt);
    loop.circles = std::move(circles);
    loops.push_back(std::move(loop));
  }
  return loops;
}

// Heals and orients one cap loop. `want_ccw` is true for the outer ring of a cap
// face and false for its holes.
inline Result<void> split_solid_finish_cap_loop(std::vector<Point3d>& points,
                                                std::vector<Point2d>& uv,
                                                std::vector<std::optional<CurvedArcSpan>>& spans,
                                                const std::vector<CircleCurve3>& circles,
                                                double tolerance, bool want_ccw) {
  if (points.size() < 3 || uv.size() != points.size() || spans.size() != points.size()) {
    return Error{Errc::invalid_topology, "split_solid_by_plane: malformed cap loop", 0u, 0u};
  }
  SplitSolidLoop loop;
  loop.points = std::move(points);
  loop.uv = std::move(uv);
  loop.spans = std::move(spans);
  for (const CircleCurve3& circle : circles) {
    const Result<void> healed = split_solid_heal_loop(loop, circle, tolerance);
    if (!healed.has_value()) return healed.error();
  }
  const double area = split_solid_ring_area(loop.uv);
  const bool is_ccw = area > 0.0;
  if (is_ccw != want_ccw) {
    split_solid_reverse_loop_keep_start(loop.points, loop.uv, loop.spans);
  }
  points = std::move(loop.points);
  uv = std::move(loop.uv);
  spans = std::move(loop.spans);
  return Result<void>{};
}

// Builds one planar cap face per top-level cut loop. Loops nested directly in a
// cap face become its holes; a loop nested in a hole starts a new cap face.
inline Result<std::vector<CurvedAssembledFace>> split_solid_build_caps(
    const std::vector<SplitSolidCutLoop>& loops, const Vec3d& normal, double offset,
    double tolerance) {
  std::vector<CurvedAssembledFace> faces;
  if (loops.empty()) return faces;
  const Result<PlaneFrame> frame_result = make_plane_frame(Plane3{normal, offset});
  if (!frame_result.has_value()) return frame_result.error();
  const PlaneFrame frame = frame_result.value();
  const std::size_t count = loops.size();
  std::vector<std::vector<Point2d>> uv(count);
  std::vector<double> area(count, 0.0);
  for (std::size_t i = 0; i < count; ++i) {
    uv[i].reserve(loops[i].points.size());
    for (const Point3d& point : loops[i].points) {
      uv[i].push_back(plane_frame_to_2d(frame, point));
    }
    area[i] = split_solid_ring_area(uv[i]);
  }
  std::vector<std::size_t> order(count);
  std::iota(order.begin(), order.end(), std::size_t{0});
  std::sort(order.begin(), order.end(), [&area](std::size_t a, std::size_t b) {
    return std::abs(area[a]) > std::abs(area[b]);
  });
  std::vector<int> parent(count, -1);
  for (std::size_t a = 0; a < count; ++a) {
    const std::size_t i = order[a];
    if (loops[i].points.size() < 3) continue;
    for (std::size_t b = a + 1; b < count; ++b) {
      const std::size_t j = order[b];
      if (!(std::abs(area[j]) > std::abs(area[i]))) continue;
      if (loops[j].points.size() < 3) continue;
      if (split_solid_point_in_ring(uv[j], uv[i][0])) parent[i] = static_cast<int>(j);
    }
  }
  std::vector<int> depth(count, 0);
  for (std::size_t i = 0; i < count; ++i) {
    int current = parent[i];
    int level = 0;
    while (current >= 0 && level < 64) {
      ++level;
      current = parent[static_cast<std::size_t>(current)];
    }
    depth[i] = level;
  }

  for (std::size_t i = 0; i < count; ++i) {
    if (depth[i] % 2 != 0) continue;
    CurvedAssembledFace face;
    face.surface = SurfaceGeometry{Plane3{normal, offset}};
    face.reversed = false;

    std::vector<std::size_t> members;
    members.push_back(i);
    for (std::size_t j = 0; j < count; ++j) {
      if (parent[j] == static_cast<int>(i)) members.push_back(j);
    }
    std::vector<CircleCurve3> circles;
    for (std::size_t member : members) {
      for (const CircleCurve3& circle : loops[member].circles) {
        split_solid_append_circles(circles, circle);
      }
    }
    for (std::size_t k = 0; k < members.size(); ++k) {
      const std::size_t member = members[k];
      std::vector<Point3d> points = loops[member].points;
      std::vector<Point2d> loop_uv = uv[member];
      std::vector<std::optional<CurvedArcSpan>> spans = loops[member].spans;
      if (spans.size() != points.size()) spans.assign(points.size(), std::nullopt);
      const Result<void> finished =
          split_solid_finish_cap_loop(points, loop_uv, spans, circles, tolerance, k == 0);
      if (!finished.has_value()) return finished.error();
      face.loops.push_back(std::move(points));
      face.loops_uv.push_back(std::move(loop_uv));
      face.arc_spans.push_back(std::move(spans));
    }
    faces.push_back(std::move(face));
  }
  return faces;
}

// Lifts one intersector curve into the face parameter domain: unwraps the
// periodic parameter, re-projects every point onto the cut plane with Newton
// steps on the analytic partials, and places the chain inside the face patch.
inline Result<std::vector<Point2d>> split_solid_cut_uv(const SsiCurve& curve,
                                                       const SurfaceGeometry& surface,
                                                       const FaceUvProjection& projection,
                                                       const std::vector<FaceUvLoop>& loops,
                                                       const Vec3d& normal, double offset,
                                                       double scale) {
  if (curve.segments.empty()) {
    return Error{Errc::numeric_failure, "split_solid_by_plane: empty intersection curve", 0u, 0u};
  }
  std::vector<Point2d> polyline;
  polyline.reserve(curve.segments.size() + 1);
  for (const SsiSegment& segment : curve.segments) {
    polyline.push_back(Point2d{segment.ua[0], segment.ua[1]});
  }
  const SsiSegment& last = curve.segments.back();
  polyline.push_back(Point2d{last.ub[0], last.ub[1]});
  if (projection.periodic_u()) {
    for (std::size_t i = 1; i < polyline.size(); ++i) {
      polyline[i][0] = polyline[i - 1][0] +
                       split_solid_wrapped_delta(polyline[i][0] - polyline[i - 1][0]);
    }
  }

  double face_u_min = loops.front().points.front()[0];
  double face_u_max = face_u_min;
  double face_v_min = loops.front().points.front()[1];
  double face_v_max = face_v_min;
  for (const FaceUvLoop& loop : loops) {
    for (const Point2d& point : loop.points) {
      face_u_min = std::min(face_u_min, point[0]);
      face_u_max = std::max(face_u_max, point[0]);
      face_v_min = std::min(face_v_min, point[1]);
      face_v_max = std::max(face_v_max, point[1]);
    }
  }

  if (curve.closed) {
    if (projection.periodic_u()) {
      double min_u = polyline.front()[0];
      for (const Point2d& point : polyline) min_u = std::min(min_u, point[0]);
      const double shift = face_u_min - min_u;
      for (Point2d& point : polyline) point[0] += shift;
    } else if (!split_solid_same2(polyline.front(), polyline.back(), 1e-12 * scale)) {
      polyline.push_back(polyline.front());
    }
  } else if (projection.periodic_u()) {
    double min_u = polyline.front()[0];
    double max_u = polyline.front()[0];
    for (const Point2d& point : polyline) {
      min_u = std::min(min_u, point[0]);
      max_u = std::max(max_u, point[0]);
    }
    double best_shift = 0.0;
    double best_overlap = -std::numeric_limits<double>::infinity();
    for (int turn = -1; turn <= 1; ++turn) {
      const double shift = static_cast<double>(turn) * kSplitSolidTwoPi;
      const double overlap = std::min(max_u + shift, face_u_max) -
                             std::max(min_u + shift, face_u_min);
      if (overlap > best_overlap) {
        best_overlap = overlap;
        best_shift = shift;
      }
    }
    for (Point2d& point : polyline) point[0] += best_shift;
  }

  // Place every sample exactly on the cut plane: the intersector stops at the
  // patch boundary within its own tolerance, while the assembler compares
  // circles at max(tolerance, 1e-9 * radius).
  const double plane_tolerance = 1e-15 * scale;
  for (Point2d& point : polyline) {
    for (int iteration = 0; iteration < 3; ++iteration) {
      const Result<SurfacePoint> evaluated = surface_eval(surface, point[0], point[1]);
      if (!evaluated.has_value()) break;
      const double gap = dot(normal, split_solid_as_vec(evaluated.value().p)) - offset;
      if (std::abs(gap) <= plane_tolerance) break;
      const double gu = dot(normal, evaluated.value().du);
      const double gv = dot(normal, evaluated.value().dv);
      const double denominator = gu * gu + gv * gv;
      if (!(denominator > 0.0)) break;
      point[0] -= gap * gu / denominator;
      point[1] -= gap * gv / denominator;
    }
  }

  const double weld = 1e-12 * scale;
  std::vector<Point2d> deduped;
  deduped.reserve(polyline.size());
  for (const Point2d& point : polyline) {
    if (!deduped.empty() && split_solid_same2(deduped.back(), point, weld)) continue;
    deduped.push_back(point);
  }
  if (deduped.size() < 2) {
    return Error{Errc::numeric_failure, "split_solid_by_plane: empty intersection curve", 0u, 0u};
  }
  return deduped;
}

}  // namespace detail

inline Result<SplitSolidResult> split_solid_by_plane(const Brep& brep, SolidId solid,
                                                     const Plane3& plane, double tolerance) {
  if (!std::isfinite(tolerance) || !(tolerance > 0.0)) {
    return Error{Errc::tolerance_violation,
                 "split_solid_by_plane: tolerance must be finite and positive"};
  }
  if (!brep.solid_alive(solid)) {
    return Error{Errc::index_out_of_range, "split_solid_by_plane: solid is not alive"};
  }
  const double raw_length = norm2(plane.normal);
  if (!std::isfinite(raw_length) || !(raw_length > 0.0)) {
    return Error{Errc::degenerate, "split_solid_by_plane: plane normal is zero"};
  }
  const Vec3d normal = plane.normal * (1.0 / raw_length);
  const double offset = plane.offset / raw_length;
  if (!std::isfinite(offset)) {
    return Error{Errc::non_finite, "split_solid_by_plane: plane is not finite"};
  }

  std::vector<FaceId> faces;
  for (ShellId shell : brep.solid_shells(solid)) {
    for (FaceId face : brep.shell_faces(shell)) faces.push_back(face);
  }
  if (faces.empty()) {
    return Error{Errc::invalid_topology, "split_solid_by_plane: solid has no faces"};
  }

  double scale = 1.0;
  for (VertexId vertex : brep.vertex_ids()) {
    const Point3d& point = brep.vertex(vertex).point;
    for (int component = 0; component < 3; ++component) {
      scale = std::max(scale, std::abs(point[component]));
    }
  }

  const Result<PlaneFrame> plane_frame_result = make_plane_frame(Plane3{normal, offset});
  if (!plane_frame_result.has_value()) return plane_frame_result.error();
  const PlaneFrame plane_frame = plane_frame_result.value();

  std::vector<detail::SplitSolidFaceInfo> infos;
  infos.reserve(faces.size());

  // Cut-plane patch bounds come from the face boundary points, not the B-rep
  // vertices: a full circular rim is a single closed edge with one vertex, so
  // the vertex set does not span the solid.
  double patch_u_min = std::numeric_limits<double>::infinity();
  double patch_u_max = -std::numeric_limits<double>::infinity();
  double patch_v_min = std::numeric_limits<double>::infinity();
  double patch_v_max = -std::numeric_limits<double>::infinity();

  for (FaceId face : faces) {
    detail::SplitSolidFaceInfo info;
    info.face = face;
    const Result<std::vector<FaceUvLoop>> loops_result = face_loops_uv(brep, face, tolerance);
    if (!loops_result.has_value()) return loops_result.error();
    info.loops = std::move(loops_result.value());
    if (info.loops.empty()) {
      return Error{Errc::invalid_topology, "split_solid_by_plane: face has no boundary loops",
                   face.index(), 0u};
    }
    const Result<detail::FaceUvProjection> projection_result =
        detail::face_uv_projection(brep, brep.face(face), tolerance);
    if (!projection_result.has_value()) return projection_result.error();
    info.projection = projection_result.value();
    info.surface = info.projection.surface;
    bool first = true;
    for (const FaceUvLoop& loop : info.loops) {
      for (const Point2d& point : loop.points) {
        const Point3d point3d = detail::face_uv_eval(info.projection, point);
        const double distance = dot(normal, detail::split_solid_as_vec(point3d)) - offset;
        for (int component = 0; component < 3; ++component) {
          scale = std::max(scale, std::abs(point3d[component]));
        }
        const Point2d plane_uv = plane_frame_to_2d(plane_frame, point3d);
        patch_u_min = std::min(patch_u_min, plane_uv[0]);
        patch_u_max = std::max(patch_u_max, plane_uv[0]);
        patch_v_min = std::min(patch_v_min, plane_uv[1]);
        patch_v_max = std::max(patch_v_max, plane_uv[1]);
        if (first) {
          info.distance_min = distance;
          info.distance_max = distance;
          first = false;
        } else {
          info.distance_min = std::min(info.distance_min, distance);
          info.distance_max = std::max(info.distance_max, distance);
        }
      }
    }
    if (!first && std::abs(info.distance_min) <= tolerance &&
        std::abs(info.distance_max) <= tolerance) {
      return Error{Errc::unsupported,
                   "split_solid_by_plane: plane is tangent or coincident with a face",
                   face.index(), 0u};
    }
    infos.push_back(std::move(info));
  }

  const double expand = std::max(tolerance, 1e-12 * scale);
  // The cut curves must end on a face boundary, never on the plane patch
  // border: a cut that runs along (or touches) the plane patch border makes the
  // tracer leave the patch and lose whole steps. Give the plane patch a
  // generous geometric margin so the face patch is the limiting domain, and so
  // a straight prediction that cuts the corner of a tangency to the border
  // cannot overshoot it.
  const double plane_margin =
      std::max(expand, 0.25 * std::max(patch_u_max - patch_u_min, patch_v_max - patch_v_min));
  SsiPatch plane_patch;
  plane_patch.surface = SurfaceGeometry{Plane3{normal, offset}};
  plane_patch.u = {patch_u_min - plane_margin, patch_u_max + plane_margin};
  plane_patch.v = {patch_v_min - plane_margin, patch_v_max + plane_margin};

  SsiOptions options;
  options.tolerance = tolerance;

  std::vector<detail::SplitSolidCutChain> chains;
  std::vector<CurvedAssembledFace> below_faces;
  std::vector<CurvedAssembledFace> above_faces;

  auto add_face = [&](std::vector<CurvedAssembledFace>& target, const SurfaceGeometry& surface,
                      const CurvedFragment& fragment, bool reversed) {
    CurvedAssembledFace face;
    face.surface = surface;
    face.loops = fragment.loops;
    face.loops_uv = fragment.loops_uv;
    face.arc_spans = fragment.arc_spans;
    face.reversed = reversed;
    detail::split_solid_canonicalize_face(face, tolerance);
    target.push_back(std::move(face));
  };

  auto classify = [&](const detail::SplitSolidFaceInfo& info,
                      const CurvedFragment& fragment) -> int {
    const std::vector<Point2d>& ring = fragment.loops_uv.front();
    Polygon2d polygon;
    polygon.vertices = ring;
    polygon.ring_starts = {0u, static_cast<std::uint32_t>(ring.size())};
    const Triangulation2 triangulation = triangulate_polygon2(polygon);
    Point2d representative = ring.front();
    if (triangulation.ok && !triangulation.triangles.empty()) {
      // Use the largest triangle: a near-degenerate first ear (three nearly
      // collinear ring vertices) has a centroid on the cut itself, which would
      // classify a fragment by the plane it lies on instead of by its body.
      double best_area = -1.0;
      for (const std::array<std::uint32_t, 3>& triangle : triangulation.triangles) {
        const Point2d& a = ring[triangle[0]];
        const Point2d& b = ring[triangle[1]];
        const Point2d& c = ring[triangle[2]];
        const double area =
            std::abs((b[0] - a[0]) * (c[1] - a[1]) - (b[1] - a[1]) * (c[0] - a[0]));
        if (area <= best_area) continue;
        best_area = area;
        representative = Point2d{0.0, 0.0};
        for (std::uint32_t index : triangle) {
          representative[0] += ring[index][0];
          representative[1] += ring[index][1];
        }
        representative[0] *= (1.0 / 3.0);
        representative[1] *= (1.0 / 3.0);
      }
    }
    const Point3d point = detail::face_uv_eval(info.projection, representative);
    const double signed_distance = dot(normal, detail::split_solid_as_vec(point)) - offset;
    return signed_distance >= 0.0 ? 1 : -1;
  };

  for (const detail::SplitSolidFaceInfo& info : infos) {
    const bool needs_cut =
        info.distance_min <= tolerance && info.distance_max >= -tolerance;
    double face_u_min = std::numeric_limits<double>::infinity();
    double face_u_max = -std::numeric_limits<double>::infinity();
    double face_v_min = std::numeric_limits<double>::infinity();
    double face_v_max = -std::numeric_limits<double>::infinity();
    for (const FaceUvLoop& loop : info.loops) {
      for (const Point2d& point : loop.points) {
        face_u_min = std::min(face_u_min, point[0]);
        face_u_max = std::max(face_u_max, point[0]);
        face_v_min = std::min(face_v_min, point[1]);
        face_v_max = std::max(face_v_max, point[1]);
      }
    }
    if (info.projection.periodic_u() && face_u_max - face_u_min > detail::kSplitSolidTwoPi) {
      face_u_max = face_u_min + detail::kSplitSolidTwoPi;
    }
    const double face_extent = std::max(face_u_max - face_u_min, face_v_max - face_v_min);
    std::vector<CircleCurve3> circles;
    std::vector<FaceCut> cuts;
    if (needs_cut) {
      CircleCurve3 circle;
      if (detail::split_solid_face_circle(*info.surface, normal, offset, plane_frame, circle)) {
        circles.push_back(circle);
      }
      // Clamp a periodic patch to exactly one period. Keep any padding inside
      // the period: a patch wider than the period would fold the seam twice.
      double u_pad = tolerance;
      if (info.projection.periodic_u()) {
        const double width = face_u_max - face_u_min;
        u_pad = std::min(tolerance,
                         std::max(0.0, 0.5 * (detail::kSplitSolidTwoPi - width)));
      }
      SsiPatch face_patch;
      face_patch.surface = *info.surface;
      face_patch.u = {face_u_min - u_pad, face_u_max + u_pad};
      face_patch.v = {face_v_min - tolerance, face_v_max + tolerance};
      if (face_patch.u[1] - face_patch.u[0] <= 0.0 ||
          face_patch.v[1] - face_patch.v[0] <= 0.0) {
        return Error{Errc::degenerate, "split_solid_by_plane: face patch is degenerate",
                     info.face.index(), 0u};
      }

      const Result<SsiResult> intersected = ssi(face_patch, plane_patch, options);
      if (!intersected.has_value()) return intersected.error();
      if (!intersected.value().issues.empty()) {
        const SsiIssue& issue = intersected.value().issues.front();
        if (issue.code == SsiDiagnostic::tangent_contact ||
            issue.code == SsiDiagnostic::coincident) {
          return Error{Errc::unsupported,
                       "split_solid_by_plane: plane is tangent or coincident with a face",
                       info.face.index(), 0u};
        }
        return Error{Errc::numeric_failure, "split_solid_by_plane: intersection failed",
                     info.face.index(), 0u};
      }
      const std::vector<FaceUvLoop>& loops = info.loops;
      const double snap_radius = std::max(tolerance, 1e-12 * face_extent);
      const double extend_limit = std::max(0.05 * face_extent, 8.0 * tolerance);
      for (const SsiCurve& curve : intersected.value().curves) {
        const Result<std::vector<Point2d>> polyline = detail::split_solid_cut_uv(
            curve, *info.surface, info.projection, loops, normal, offset, scale);
        if (!polyline.has_value()) return polyline.error();
        std::vector<Point2d> lifted = polyline.value();
        detail::split_solid_extend_cut_end(lifted, true, *info.surface, normal, offset,
                                           face_u_min, face_u_max, face_v_min, face_v_max,
                                           snap_radius, extend_limit);
        detail::split_solid_extend_cut_end(lifted, false, *info.surface, normal, offset,
                                           face_u_min, face_u_max, face_v_min, face_v_max,
                                           snap_radius, extend_limit);
        detail::split_solid_dedupe_polyline(lifted, std::max(tolerance, 1e-12 * scale));
        std::vector<std::vector<Point2d>> clipped = detail::split_solid_clip_polyline(
            lifted, face_u_min, face_u_max, face_v_min, face_v_max,
            std::max(1e-12 * scale, tolerance));
        const std::vector<std::vector<Point2d>> loop_points = [&info] {
          std::vector<std::vector<Point2d>> points;
          points.reserve(info.loops.size());
          for (const FaceUvLoop& loop : info.loops) points.push_back(loop.points);
          return points;
        }();
        for (std::vector<Point2d>& piece : clipped) {
          if (piece.size() < 2) continue;
          const detail::TrimFaceSnap front =
              detail::trim_face_find_snap(loop_points, piece.front(), snap_radius);
          if (front.boundary) piece.front() = front.point;
          const detail::TrimFaceSnap back =
              detail::trim_face_find_snap(loop_points, piece.back(), snap_radius);
          if (back.boundary) piece.back() = back.point;
          cuts.push_back(FaceCut{piece});
        }
      }
    }

    const Result<std::vector<CurvedFragment>> split_result =
        split_face_curved(brep, info.face, cuts, tolerance);
    if (!split_result.has_value()) return split_result.error();
    for (CurvedFragment fragment : split_result.value()) {
      if (fragment.loops.empty()) continue;
      const int side = classify(info, fragment);
      const Result<void> healed =
          detail::split_solid_heal_fragment(fragment, circles, tolerance);
      if (!healed.has_value()) return healed.error();
      add_face(side >= 0 ? above_faces : below_faces, *info.surface, fragment,
               brep.face(info.face).reversed);
    }
    // One chain per cut polyline per face, independent of how many fragments
    // the cut produced.
    if (needs_cut) {
      for (const FaceCut& cut : cuts) {
        detail::SplitSolidCutChain chain;
        chain.points.reserve(cut.uv.size());
        for (const Point2d& point : cut.uv) {
          chain.points.push_back(detail::face_uv_eval(info.projection, point));
        }
        if (!circles.empty()) chain.circle = circles.front();
        chains.push_back(std::move(chain));
      }
    }
  }

  if (below_faces.empty() || above_faces.empty()) {
    return Error{Errc::degenerate, "split_solid_by_plane: plane does not split the solid"};
  }

  const double weld = std::max(64.0 * tolerance, 1e-9 * scale);
  detail::split_solid_dedupe_chains(chains, weld);
  Result<std::vector<detail::SplitSolidCutLoop>> stitched =
      detail::split_solid_stitch_cuts(chains, weld);
  if (!stitched.has_value()) return stitched.error();

  const Result<std::vector<CurvedAssembledFace>> below_caps =
      detail::split_solid_build_caps(stitched.value(), normal, offset, tolerance);
  if (!below_caps.has_value()) return below_caps.error();
  const Result<std::vector<CurvedAssembledFace>> above_caps = detail::split_solid_build_caps(
      stitched.value(), normal * -1.0, -offset, tolerance);
  if (!above_caps.has_value()) return above_caps.error();
  for (CurvedAssembledFace cap : below_caps.value()) {
    detail::split_solid_canonicalize_face(cap, tolerance);
    below_faces.push_back(std::move(cap));
  }
  for (CurvedAssembledFace cap : above_caps.value()) {
    detail::split_solid_canonicalize_face(cap, tolerance);
    above_faces.push_back(std::move(cap));
  }

  const Result<Brep> below = assemble_curved_solid(below_faces, tolerance);
  if (!below.has_value()) return below.error();
  const Result<Brep> above = assemble_curved_solid(above_faces, tolerance);
  if (!above.has_value()) return above.error();

  SplitSolidResult result;
  result.below = std::move(below.value());
  result.above = std::move(above.value());
  return result;
}

}  // namespace v1
}  // namespace geom
