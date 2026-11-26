#pragma once

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>
#include <span>
#include <vector>

#include <geom/core/error.hpp>
#include <geom/core/result.hpp>
#include <geom/primitives/point.hpp>
#include <geom/primitives/vec.hpp>

namespace geom {
inline namespace v1 {
namespace detail {

// One half-space dot(normal, x) <= offset. `normal` need not be unit length:
// the enumeration normalizes it and divides `offset` by the same length, so
// callers may pass plane equations in any convenient scaling.
struct HalfSpacePlane {
  Vec3d normal{};
  double offset = 0.0;
};

// Enumerates the vertices of the convex cell that is the intersection of the
// given half-spaces. Every candidate is the solution of three plane equations,
// kept only when it satisfies all the other planes within `tolerance`, and
// deduplicated. The result is double-precision by design (the open-shell
// operation is allowed to model approximately), so no exact predicates are
// used here. An empty intersection is an empty vector, not an error.
inline Result<std::vector<Point3d>> half_space_vertices(
    std::span<const HalfSpacePlane> planes, double tolerance) {
  if (!std::isfinite(tolerance) || tolerance <= 0.0) {
    return Error{Errc::tolerance_violation,
                 "half_space: tolerance must be finite and positive"};
  }

  const std::size_t count = planes.size();
  std::vector<Vec3d> normals(count);
  std::vector<double> offsets(count);
  for (std::size_t i = 0; i < count; ++i) {
    const Vec3d n = planes[i].normal;
    const double length = std::sqrt(dot(n, n));
    if (!std::isfinite(length) || length <= 0.0) {
      return Error{Errc::degenerate,
                   "half_space: plane normal is zero or non-finite"};
    }
    if (!std::isfinite(planes[i].offset)) {
      return Error{Errc::non_finite, "half_space: plane offset is not finite"};
    }
    normals[i] = Vec3d{n[0] / length, n[1] / length, n[2] / length};
    offsets[i] = planes[i].offset / length;
  }

  std::vector<Point3d> accepted;
  for (std::size_t i = 0; i < count; ++i) {
    for (std::size_t j = i + 1; j < count; ++j) {
      for (std::size_t k = j + 1; k < count; ++k) {
        // Cramer's rule: p solves n_i . p = d_i for the three planes, with the
        // determinant written in mixed-product form.
        const Vec3d c0 = cross(normals[j], normals[k]);
        const Vec3d c1 = cross(normals[k], normals[i]);
        const Vec3d c2 = cross(normals[i], normals[j]);
        const double det = dot(normals[i], c0);
        // Exactly parallel planes have no single intersection point. Skip
        // rather than fail, since other triples may describe the same vertex.
        if (!std::isfinite(det) || det == 0.0) continue;
        const Vec3d p{
            (offsets[i] * c0[0] + offsets[j] * c1[0] + offsets[k] * c2[0]) / det,
            (offsets[i] * c0[1] + offsets[j] * c1[1] + offsets[k] * c2[1]) / det,
            (offsets[i] * c0[2] + offsets[j] * c1[2] + offsets[k] * c2[2]) / det};
        if (!std::isfinite(p[0]) || !std::isfinite(p[1]) ||
            !std::isfinite(p[2])) {
          continue;
        }

        // Relative slack for this candidate, with no absolute floor: it scales
        // with the triple's offsets and the solution's own magnitude so that
        // tiny and huge models keep the same relative behaviour.
        const double point_scale =
            std::max({std::fabs(offsets[i]), std::fabs(offsets[j]),
                      std::fabs(offsets[k]), std::fabs(p[0]),
                      std::fabs(p[1]), std::fabs(p[2])});
        const double slack = std::max(tolerance, 1e-12 * point_scale);

        // Conditioning check: if the intersection is so ill-conditioned that
        // the rounding error bound of the Cramer solve (about
        // eps * point_scale / |det|) exceeds the slack, the recovered point
        // says nothing reliable. Skip the triple instead of admitting a
        // spurious vertex.
        if (32.0 * std::numeric_limits<double>::epsilon() * point_scale /
                std::fabs(det) >
            slack) {
          continue;
        }

        bool inside = true;
        for (std::size_t l = 0; l < count && inside; ++l) {
          inside = dot(normals[l], p) <= offsets[l] + slack;
        }
        if (!inside) continue;

        // Duplicates arise when several plane triples meet at the same vertex.
        // Keep the first representative so the output order is deterministic
        // for a fixed plane order.
        bool duplicate = false;
        for (const Point3d& q : accepted) {
          duplicate = duplicate ||
                      (std::fabs(p[0] - q[0]) <= slack &&
                       std::fabs(p[1] - q[1]) <= slack &&
                       std::fabs(p[2] - q[2]) <= slack);
        }
        if (!duplicate) accepted.push_back(Point3d{p[0], p[1], p[2]});
      }
    }
  }
  return accepted;
}

}  // namespace detail
}  // namespace v1
}  // namespace geom
