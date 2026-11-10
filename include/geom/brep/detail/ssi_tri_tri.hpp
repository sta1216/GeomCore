#pragma once

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <geom/primitives/point.hpp>
#include <geom/primitives/vec.hpp>

namespace geom {
inline namespace v1 {

// Classification of the intersection of two closed triangles.
//   none     - disjoint (strictly separated by one of the two planes)
//   segment  - the intersection is a non-degenerate segment (a and b set)
//   coplanar - the triangles lie in the same plane; no segment is produced
//   touching - a single contact point, or a degenerate/non-finite input
enum class TriTriKind : std::uint8_t { none, segment, coplanar, touching };

struct TriTriResult {
  TriTriKind kind = TriTriKind::none;
  Point3d a{};
  Point3d b{};
};

// True when all three vertices are finite.
inline bool tri_finite(Point3d p0, Point3d p1, Point3d p2) noexcept {
  return std::isfinite(p0[0]) && std::isfinite(p0[1]) && std::isfinite(p0[2]) &&
         std::isfinite(p1[0]) && std::isfinite(p1[1]) && std::isfinite(p1[2]) &&
         std::isfinite(p2[0]) && std::isfinite(p2[1]) && std::isfinite(p2[2]);
}

namespace detail {

inline bool ssi_tri_vec_finite(Vec3d v) noexcept {
  return std::isfinite(v[0]) && std::isfinite(v[1]) && std::isfinite(v[2]);
}

inline bool ssi_tri_less(Point3d a, Point3d b) noexcept {
  if (a[0] != b[0]) {
    return a[0] < b[0];
  }
  if (a[1] != b[1]) {
    return a[1] < b[1];
  }
  return a[2] < b[2];
}

// A triangle is degenerate when an edge is no longer than the tolerance or its
// area is no larger than tolerance^2. The normal is passed in so callers can
// validate its finiteness first.
inline bool ssi_tri_degenerate(Point3d p0, Point3d p1, Point3d p2, Vec3d normal,
                               double tolerance) noexcept {
  if (norm2(p1 - p0) <= tolerance || norm2(p2 - p1) <= tolerance || norm2(p0 - p2) <= tolerance) {
    return true;
  }
  const double area = 0.5 * norm2(normal);
  return !(area > tolerance * tolerance);
}

// Signed distance from p to the plane through origin with the given (not
// necessarily unit) normal.
inline double ssi_tri_signed_distance(Point3d p, Point3d origin, Vec3d normal) noexcept {
  return dot(p - origin, normal) / norm2(normal);
}

inline int ssi_tri_dominant_axis(Vec3d normal) noexcept {
  const double ax = std::fabs(normal[0]);
  const double ay = std::fabs(normal[1]);
  const double az = std::fabs(normal[2]);
  if (ax >= ay && ax >= az) {
    return 0;
  }
  if (ay >= az) {
    return 1;
  }
  return 2;
}

inline void ssi_tri_project(Point3d p, int axis, double& u, double& v) noexcept {
  if (axis == 0) {
    u = p[1];
    v = p[2];
  } else if (axis == 1) {
    u = p[0];
    v = p[2];
  } else {
    u = p[0];
    v = p[1];
  }
}

// Point-in-triangle test in the projection that drops the triangle normal's
// dominant axis. Points on the boundary count as inside; tolerance is compared
// against the true (unprojected) edge distance.
inline bool ssi_tri_point_in_triangle(Point3d p, Point3d t0, Point3d t1, Point3d t2, int axis,
                                      double tolerance) noexcept {
  double px = 0.0;
  double py = 0.0;
  double ax = 0.0;
  double ay = 0.0;
  double bx = 0.0;
  double by = 0.0;
  double cx = 0.0;
  double cy = 0.0;
  ssi_tri_project(p, axis, px, py);
  ssi_tri_project(t0, axis, ax, ay);
  ssi_tri_project(t1, axis, bx, by);
  ssi_tri_project(t2, axis, cx, cy);

  const double area2 = (bx - ax) * (cy - ay) - (by - ay) * (cx - ax);
  if (area2 == 0.0) {
    return false;
  }
  const double sign = area2 > 0.0 ? 1.0 : -1.0;

  const auto inside_edge = [&](double ex0, double ey0, double ex1, double ey1) {
    const double ex = ex1 - ex0;
    const double ey = ey1 - ey0;
    const double edge_length = std::hypot(ex, ey);
    const double cross = ex * (py - ey0) - ey * (px - ex0);
    return sign * cross >= -tolerance * edge_length;
  };

  return inside_edge(ax, ay, bx, by) && inside_edge(bx, by, cx, cy) && inside_edge(cx, cy, ax, ay);
}

inline void ssi_tri_add_point(std::array<Point3d, 8>& points, std::size_t& count, Point3d p,
                              double tolerance) noexcept {
  for (std::size_t i = 0; i < count; ++i) {
    if (norm2(p - points[i]) <= tolerance) {
      return;
    }
  }
  if (count < points.size()) {
    points[count] = p;
    ++count;
  }
}

// Adds the part of edge (a, b) lying on the plane, using the signed distances
// da and db of its endpoints. Endpoints on the plane are added directly.
inline void ssi_tri_collect_edge(Point3d a, Point3d b, double da, double db, double tolerance,
                                 std::array<Point3d, 8>& points, std::size_t& count) noexcept {
  const bool a_on = std::fabs(da) <= tolerance;
  const bool b_on = std::fabs(db) <= tolerance;
  if (a_on) {
    ssi_tri_add_point(points, count, a, tolerance);
  }
  if (b_on) {
    ssi_tri_add_point(points, count, b, tolerance);
  }
  if (!a_on && !b_on && ((da > 0.0) != (db > 0.0))) {
    const double t = da / (da - db);
    const Point3d p{a[0] + (b[0] - a[0]) * t, a[1] + (b[1] - a[1]) * t, a[2] + (b[2] - a[2]) * t};
    ssi_tri_add_point(points, count, p, tolerance);
  }
}

inline TriTriResult ssi_tri_result_from_points(const std::array<Point3d, 8>& points,
                                               std::size_t count) noexcept {
  TriTriResult result;
  if (count == 0) {
    return result;
  }
  if (count == 1) {
    result.kind = TriTriKind::touching;
    result.a = points[0];
    result.b = points[0];
    return result;
  }

  std::size_t ia = 0;
  std::size_t ib = 1;
  double farthest = -1.0;
  for (std::size_t i = 0; i < count; ++i) {
    for (std::size_t j = i + 1; j < count; ++j) {
      const double d = norm2(points[j] - points[i]);
      if (d > farthest) {
        farthest = d;
        ia = i;
        ib = j;
      }
    }
  }
  if (!(farthest > 0.0)) {
    result.kind = TriTriKind::touching;
    result.a = points[0];
    result.b = points[0];
    return result;
  }

  result.kind = TriTriKind::segment;
  if (ssi_tri_less(points[ia], points[ib])) {
    result.a = points[ia];
    result.b = points[ib];
  } else {
    result.a = points[ib];
    result.b = points[ia];
  }
  return result;
}

}  // namespace detail

// Intersects the closed triangles (p0, p1, p2) and (q0, q1, q2).
//
// The tolerance is the only geometric authority; there is no global epsilon.
// It must be finite and strictly positive: zero, negative, NaN and infinity all
// classify as touching. The six vertices must be finite, and the two triangles
// should have coordinates of comparable scale (within roughly 1e80 of each
// other); beyond that the smaller triangle's normal underflows in the
// normalized frame and the call conservatively classifies as touching.
// The call is normalized into a power-of-two frame before any geometric
// predicate is evaluated, so ordinary finite coordinates cannot overflow or
// underflow the cross/dot products. A tolerance that underflows to zero in
// that frame is allowed, and only exact coplanarity then counts.
//
// Degenerate or non-finite input classifies as touching. Coplanar triangles
// always classify as coplanar, even when they are disjoint. Segment endpoints
// are ordered lexicographically by (x, y, z) for determinism; for a touching
// result a and b both hold the contact point.
inline TriTriResult tri_tri_intersect(Point3d p0, Point3d p1, Point3d p2, Point3d q0, Point3d q1,
                                      Point3d q2, double tolerance) {
  TriTriResult result;
  if (!tri_finite(p0, p1, p2) || !tri_finite(q0, q1, q2)) {
    result.kind = TriTriKind::touching;
    return result;
  }
  if (!(tolerance > 0.0) || !std::isfinite(tolerance)) {
    result.kind = TriTriKind::touching;
    return result;
  }

  // Choose a power of two that puts the largest coordinate in [1, 2). The
  // scale stays finite for every finite input because exponent - 1 <= 1023.
  const std::array<Point3d, 6> input{p0, p1, p2, q0, q1, q2};
  double max_abs = 0.0;
  for (const Point3d& p : input) {
    max_abs = std::max({max_abs, std::fabs(p[0]), std::fabs(p[1]), std::fabs(p[2])});
  }
  if (max_abs == 0.0) {
    result.kind = TriTriKind::touching;
    return result;
  }
  int exponent = 0;
  std::frexp(max_abs, &exponent);
  const double scale = std::ldexp(1.0, exponent - 1);
  const double tol = tolerance / scale;
  const auto rescale = [scale](Point3d p) {
    return Point3d{p[0] / scale, p[1] / scale, p[2] / scale};
  };
  p0 = rescale(p0);
  p1 = rescale(p1);
  p2 = rescale(p2);
  q0 = rescale(q0);
  q1 = rescale(q1);
  q2 = rescale(q2);

  const Vec3d np = cross(p1 - p0, p2 - p0);
  const Vec3d nq = cross(q1 - q0, q2 - q0);
  if (!detail::ssi_tri_vec_finite(np) || !detail::ssi_tri_vec_finite(nq)) {
    result.kind = TriTriKind::touching;
    return result;
  }
  if (detail::ssi_tri_degenerate(p0, p1, p2, np, tol) ||
      detail::ssi_tri_degenerate(q0, q1, q2, nq, tol)) {
    result.kind = TriTriKind::touching;
    return result;
  }

  double dq[3] = {};
  double dp[3] = {};
  bool all_q_near = true;
  bool all_p_near = true;
  for (int i = 0; i < 3; ++i) {
    const Point3d pi = i == 0 ? p0 : (i == 1 ? p1 : p2);
    const Point3d qi = i == 0 ? q0 : (i == 1 ? q1 : q2);
    dq[i] = detail::ssi_tri_signed_distance(pi, q0, nq);
    dp[i] = detail::ssi_tri_signed_distance(qi, p0, np);
    all_q_near = all_q_near && std::fabs(dq[i]) <= tol;
    all_p_near = all_p_near && std::fabs(dp[i]) <= tol;
  }
  if (all_q_near || all_p_near) {
    result.kind = TriTriKind::coplanar;
    return result;
  }

  const bool p_above = dq[0] > tol && dq[1] > tol && dq[2] > tol;
  const bool p_below = dq[0] < -tol && dq[1] < -tol && dq[2] < -tol;
  const bool q_above = dp[0] > tol && dp[1] > tol && dp[2] > tol;
  const bool q_below = dp[0] < -tol && dp[1] < -tol && dp[2] < -tol;
  if (p_above || p_below || q_above || q_below) {
    return result;
  }

  const int q_axis = detail::ssi_tri_dominant_axis(nq);
  std::array<Point3d, 8> p_points{};
  std::size_t p_count = 0;
  detail::ssi_tri_collect_edge(p0, p1, dq[0], dq[1], tol, p_points, p_count);
  detail::ssi_tri_collect_edge(p1, p2, dq[1], dq[2], tol, p_points, p_count);
  detail::ssi_tri_collect_edge(p2, p0, dq[2], dq[0], tol, p_points, p_count);

  const int p_axis = detail::ssi_tri_dominant_axis(np);
  std::array<Point3d, 8> q_points{};
  std::size_t q_count = 0;
  detail::ssi_tri_collect_edge(q0, q1, dp[0], dp[1], tol, q_points, q_count);
  detail::ssi_tri_collect_edge(q1, q2, dp[1], dp[2], tol, q_points, q_count);
  detail::ssi_tri_collect_edge(q2, q0, dp[2], dp[0], tol, q_points, q_count);

  std::array<Point3d, 8> inside{};
  std::size_t inside_count = 0;
  for (std::size_t i = 0; i < p_count; ++i) {
    if (detail::ssi_tri_point_in_triangle(p_points[i], q0, q1, q2, q_axis, tol)) {
      detail::ssi_tri_add_point(inside, inside_count, p_points[i], tol);
    }
  }
  for (std::size_t i = 0; i < q_count; ++i) {
    if (detail::ssi_tri_point_in_triangle(q_points[i], p0, p1, p2, p_axis, tol)) {
      detail::ssi_tri_add_point(inside, inside_count, q_points[i], tol);
    }
  }

  result = detail::ssi_tri_result_from_points(inside, inside_count);
  if (result.kind == TriTriKind::segment || result.kind == TriTriKind::touching) {
    for (int i = 0; i < 3; ++i) {
      result.a[i] *= scale;
      result.b[i] *= scale;
    }
  }
  return result;
}

}  // namespace v1
}  // namespace geom
