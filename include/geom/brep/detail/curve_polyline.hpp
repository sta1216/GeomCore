#pragma once

#include <cmath>
#include <cstddef>
#include <span>
#include <utility>
#include <vector>

#include <geom/core/error.hpp>
#include <geom/core/result.hpp>
#include <geom/curves/nurbs_curve.hpp>
#include <geom/primitives/point.hpp>

namespace geom {
inline namespace v1 {

// Builds a degree-1 (polyline) NURBS curve through the given points.
//
// The parameter domain is [0, 1] with clamped, uniform knots, so the curve
// passes through points[i] exactly at parameter i / (n - 1) for a span of n
// points. All weights are 1, which makes every segment the straight chord
// between two consecutive input points. Fewer than two points is degenerate;
// a point with any non-finite coordinate is rejected as non_finite.
inline Result<NurbsCurve> make_polyline_curve(std::span<const Point3d> points) {
  if (points.size() < 2) {
    return Error{Errc::degenerate, "make_polyline_curve: needs at least two points", 0u, 0u};
  }
  for (const Point3d& point : points) {
    if (!std::isfinite(point[0]) || !std::isfinite(point[1]) || !std::isfinite(point[2])) {
      return Error{Errc::non_finite, "make_polyline_curve: point is not finite", 0u, 0u};
    }
  }

  const std::size_t count = points.size();
  const double last = static_cast<double>(count - 1);
  // Degree 1 needs count + 2 knots: a double knot at each end and one interior
  // knot per segment boundary.
  std::vector<double> knots(count + 2, 0.0);
  for (std::size_t i = 1; i + 1 < count; ++i) {
    knots[i + 1] = static_cast<double>(i) / last;
  }
  knots[count] = 1.0;
  knots[count + 1] = 1.0;

  std::vector<Point3d> control(points.begin(), points.end());
  std::vector<double> weights(count, 1.0);
  return NurbsCurve(1, std::move(knots), std::move(control), std::move(weights));
}

}  // namespace v1
}  // namespace geom
