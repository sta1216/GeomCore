#pragma once

#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <span>
#include <vector>

#include <geom/numeric/predicates2.hpp>
#include <geom/primitives/point.hpp>

namespace geom {
inline namespace v1 {

// Counter-clockwise convex hull of 2D points, returned as indices into `pts`.
// The first index is not repeated, collinear points keep only the two extreme
// ones, exact duplicate coordinates are collapsed, and degenerate inputs return
// 0, 1 or 2 indices without a repeated point.
//
// Ordering rule: points are sorted lexicographically by (x, then y); when two
// coordinates are exactly equal the smaller input index comes first, so sorting
// is deterministic and duplicate coordinates keep their smallest input index.
// Because of that rule an all-collinear input yields its leftmost endpoint first
// and its rightmost endpoint second, independent of the input order.
//
// Collinearity is decided by sign_orient2d, never by a floating-point cross
// product, so near-collinear configurations are classified exactly.
//
// Preconditions (documented, checked by assert in debug builds):
//   - every coordinate is finite. A NaN coordinate violates the strict weak
//     ordering the sort comparator relies on, and sign_orient2d would classify
//     the affected orientation as 0.
//   - pts.size() <= 2^32 - 1, because indices are std::uint32_t. Larger inputs
//     are unsupported and would silently wrap their indices; the assert makes
//     that a debug-build failure instead of a wrong hull.
// Both are caller errors: this function has no Result return and does not try
// to repair its input, so neither check changes the answer in release builds.
inline std::vector<std::uint32_t> convex_hull2d(std::span<const Point2d> pts) {
  assert(pts.size() <= (std::size_t{1} << 32) - 1);
#ifndef NDEBUG
  for (const Point2d& p : pts) {
    assert(std::isfinite(p[0]) && std::isfinite(p[1]));
  }
#endif

  std::vector<std::uint32_t> idx(pts.size());
  for (std::uint32_t i = 0; i < pts.size(); ++i) idx[i] = i;
  std::sort(idx.begin(), idx.end(), [&](std::uint32_t a, std::uint32_t b) {
    if (pts[a][0] != pts[b][0]) return pts[a][0] < pts[b][0];
    if (pts[a][1] != pts[b][1]) return pts[a][1] < pts[b][1];
    return a < b;
  });
  // Remove exact duplicates. The sort placed equal coordinates next to each
  // other in input-index order, so unique() keeps the smallest input index.
  idx.erase(std::unique(idx.begin(), idx.end(), [&](std::uint32_t a, std::uint32_t b) {
              return pts[a][0] == pts[b][0] && pts[a][1] == pts[b][1];
            }),
            idx.end());

  if (idx.size() <= 2) return idx;

  const auto cross = [&](std::uint32_t o, std::uint32_t a, std::uint32_t b) {
    return sign_orient2d(pts[o][0], pts[o][1], pts[a][0], pts[a][1], pts[b][0], pts[b][1]);
  };

  // Monotone chain: one index buffer reused as the hull stack, no recursion.
  std::vector<std::uint32_t> hull(2 * idx.size());
  std::size_t k = 0;
  for (std::size_t i = 0; i < idx.size(); ++i) {
    while (k >= 2 && cross(hull[k - 2], hull[k - 1], idx[i]) <= 0) --k;
    hull[k++] = idx[i];
  }
  const std::size_t lower = k + 1;
  for (std::size_t i = idx.size() - 1; i > 0; --i) {
    while (k >= lower && cross(hull[k - 2], hull[k - 1], idx[i - 1]) <= 0) --k;
    hull[k++] = idx[i - 1];
  }
  hull.resize(k - 1);
  return hull;
}

}  // namespace v1
}  // namespace geom
