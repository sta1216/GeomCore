#pragma once

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <limits>
#include <span>
#include <utility>
#include <vector>

#include <geom/primitives/box.hpp>
#include <geom/primitives/point.hpp>
#include <geom/primitives/vec.hpp>

namespace geom {
inline namespace v1 {

// Minimal triangle and ray types for the M1 index. They are deliberately not
// more general: the BVH only needs a primitive with a bounding box and a
// closest-point query.
struct Triangle3d {
  Point3d a, b, c;
};

struct Ray3d {
  Point3d origin;
  Vec3d dir;
};

inline Box3d triangle_box(const Triangle3d& t) noexcept {
  Box3d b = Box3d::empty();
  b.expand(t.a);
  b.expand(t.b);
  b.expand(t.c);
  return b;
}

// Closest point on a triangle (Ericson, Real-Time Collision Detection 5.1.5),
// squared so callers can compare distances without a square root. Uses
// floating point, not the exact predicates: it is a distance query, not a
// degenerate-input decision.
inline double squared_distance(Point3d p, const Triangle3d& t) noexcept {
  const Vec3d ab = t.b - t.a, ac = t.c - t.a, ap = p - t.a;
  const double d1 = dot(ab, ap), d2 = dot(ac, ap);
  if (d1 <= 0.0 && d2 <= 0.0) return dot(ap, ap);
  const Vec3d bp = p - t.b;
  const double d3 = dot(ab, bp), d4 = dot(ac, bp);
  if (d3 >= 0.0 && d4 <= d3) return dot(bp, bp);
  const double vc = d1 * d4 - d3 * d2;
  if (vc <= 0.0 && d1 >= 0.0 && d3 <= 0.0) {
    const double v = d1 / (d1 - d3);
    const Vec3d q = ab * v - ap;
    return dot(q, q);
  }
  const Vec3d cp = p - t.c;
  const double d5 = dot(ab, cp), d6 = dot(ac, cp);
  if (d6 >= 0.0 && d5 <= d6) return dot(cp, cp);
  const double vb = d5 * d2 - d1 * d6;
  if (vb <= 0.0 && d2 >= 0.0 && d6 <= 0.0) {
    const double w = d2 / (d2 - d6);
    const Vec3d q = ac * w - ap;
    return dot(q, q);
  }
  const double va = d3 * d6 - d5 * d4;
  if (va <= 0.0 && (d4 - d3) >= 0.0 && (d5 - d6) >= 0.0) {
    const double w = (d4 - d3) / ((d4 - d3) + (d5 - d6));
    const Vec3d q = (t.c - t.b) * w - bp;
    return dot(q, q);
  }
  const double denom = 1.0 / (va + vb + vc);
  const double v = vb * denom;
  const double w = vc * denom;
  const Vec3d q = ab * v + ac * w - ap;
  return dot(q, q);
}

// Bounding volume hierarchy over triangles, built with a 12-bin surface area
// heuristic and queried iteratively (no recursion on the query path).
//
// Memory layout. All primitives live in one contiguous vector in build order;
// `order_` is a permutation of primitive indices and every leaf covers a
// contiguous range `[first, first + count)` of that permutation. Nodes live in
// one contiguous vector; a node is a leaf when `count != 0`, otherwise it is an
// interior node with `left`/`right` child indices (int32, never negative for a
// valid node). A depth guard bounds the tree at 64 levels, so the explicit
// query stacks never exceed 65 entries.
//
// Build: 12 bins over the centroid bounds of the current range, cost
// `traversal + (left_area / area) * nleft + (right_area / area) * nright`
// (the traversal constant cancels between candidates but is kept so the value
// is comparable to a classic SAH cost). If the best bin split is not strictly
// cheaper than a median split by centroid, or would leave one side empty, the
// median split is used instead; it sorts the range by (centroid axis, index),
// so the tree is deterministic. Leaves hold at most 4 primitives unless the
// depth guard forces a larger one.
//
// `bounds()` is the root node's box, therefore exactly the union of all
// primitive boxes. An empty tree has no node at all and returns an empty box.
//
// Ray contract: `dir` is expected to be finite and non-zero; t is measured in
// units of |dir|, so a non-normalised direction is fine. Hits with t <= 0 are
// behind the origin and are ignored; hits must satisfy t < tmax. Equal t ties
// resolve to the smallest primitive index.
class Bvh {
 public:
  struct Hit {
    std::uint32_t prim = 0;
    double t = 0.0;
  };

  void build(std::span<const Triangle3d> triangles) {
    prims_.clear();
    order_.clear();
    nodes_.clear();
    prims_.reserve(triangles.size());
    order_.reserve(triangles.size());
    for (std::size_t i = 0; i < triangles.size(); ++i) {
      Prim p;
      p.tri = triangles[i];
      p.box = triangle_box(triangles[i]);
      p.centroid = Point3d{(p.box.lo[0] + p.box.hi[0]) * 0.5,
                           (p.box.lo[1] + p.box.hi[1]) * 0.5,
                           (p.box.lo[2] + p.box.hi[2]) * 0.5};
      prims_.push_back(p);
      order_.push_back(static_cast<std::uint32_t>(i));
    }
    if (prims_.empty()) return;
    nodes_.reserve(2 * prims_.size());
    build_node(0, static_cast<std::uint32_t>(prims_.size()), 0);
  }

  const Box3d& bounds() const noexcept {
    static const Box3d kEmpty = Box3d::empty();
    return nodes_.empty() ? kEmpty : nodes_[0].box;
  }

  std::uint32_t size() const noexcept { return static_cast<std::uint32_t>(prims_.size()); }

  bool raycast(const Ray3d& ray, double tmax, Hit& hit) const {
    if (nodes_.empty()) return false;
    const Vec3d inv{inverse_dir(ray.dir[0]), inverse_dir(ray.dir[1]), inverse_dir(ray.dir[2])};
    std::array<std::int32_t, 128> stack{};
    int sp = 0;
    stack[sp++] = 0;
    bool found = false;
    double best = tmax;
    std::uint32_t best_prim = 0;
    while (sp > 0) {
      const Node& n = nodes_[static_cast<std::size_t>(stack[--sp])];
      if (!slab(n.box, ray, inv, best)) continue;
      if (n.count > 0) {
        for (std::uint32_t i = n.first; i < n.first + n.count; ++i) {
          const std::uint32_t prim = order_[i];
          double t = 0.0;
          if (!intersect_tri(ray, prims_[prim].tri, t)) continue;
          if (t < best || (t == best && found && prim < best_prim)) {
            best = t;
            best_prim = prim;
            found = true;
          }
        }
      } else {
        stack[sp++] = n.left;
        stack[sp++] = n.right;
      }
    }
    if (found) {
      hit.prim = best_prim;
      hit.t = best;
    }
    return found;
  }

  // Returns the index of the nearest primitive and writes its squared distance;
  // equal distances resolve to the smallest primitive index. An empty tree
  // returns 0xFFFFFFFF and an infinite distance.
  std::uint32_t nearest(const Point3d& p, double& best_d2) const {
    if (nodes_.empty()) {
      best_d2 = std::numeric_limits<double>::infinity();
      return 0xFFFFFFFFu;
    }
    std::array<std::int32_t, 128> stack{};
    int sp = 0;
    stack[sp++] = 0;
    best_d2 = std::numeric_limits<double>::infinity();
    std::uint32_t best = 0xFFFFFFFFu;
    while (sp > 0) {
      const Node& n = nodes_[static_cast<std::size_t>(stack[--sp])];
      // Equality is not pruned: a different primitive at exactly best_d2 may
      // still win the index tie-break.
      if (box_distance2(n.box, p) > best_d2) continue;
      if (n.count > 0) {
        for (std::uint32_t i = n.first; i < n.first + n.count; ++i) {
          const std::uint32_t prim = order_[i];
          const double d2 = squared_distance(p, prims_[prim].tri);
          if (d2 < best_d2 || (d2 == best_d2 && prim < best)) {
            best_d2 = d2;
            best = prim;
          }
        }
      } else {
        stack[sp++] = n.left;
        stack[sp++] = n.right;
      }
    }
    return best;
  }

 private:
  struct Prim {
    Triangle3d tri;
    Box3d box;
    Point3d centroid;
  };

  struct Node {
    Box3d box;
    std::uint32_t first = 0;   // leaf: first index into order_
    std::uint32_t count = 0;   // leaf: number of primitives; 0 => interior
    std::int32_t left = -1;    // interior child index
    std::int32_t right = -1;   // interior child index
  };

  static constexpr int kBins = 12;
  static constexpr std::uint32_t kMaxLeaf = 4;
  static constexpr int kMaxDepth = 64;
  static constexpr double kTraversalCost = 1.0;

  std::int32_t build_node(std::uint32_t begin, std::uint32_t end, int depth) {
    const auto node_index = static_cast<std::int32_t>(nodes_.size());
    nodes_.push_back(Node{});
    Box3d box = Box3d::empty();
    for (std::uint32_t i = begin; i < end; ++i) box.expand(prims_[order_[i]].box);
    nodes_[static_cast<std::size_t>(node_index)].box = box;

    const std::uint32_t count = end - begin;
    if (count <= kMaxLeaf || depth >= kMaxDepth) {
      nodes_[static_cast<std::size_t>(node_index)].first = begin;
      nodes_[static_cast<std::size_t>(node_index)].count = count;
      return node_index;
    }

    Box3d centroid_box = Box3d::empty();
    for (std::uint32_t i = begin; i < end; ++i) centroid_box.expand(prims_[order_[i]].centroid);
    const Vec3d extent = centroid_box.diagonal();
    int axis = 0;
    if (extent[1] > extent[axis]) axis = 1;
    if (extent[2] > extent[axis]) axis = 2;
    const double lo = centroid_box.lo[axis];
    const double hi = centroid_box.hi[axis];
    const double area = surface_area(box);

    // Candidate 1: 12-bin SAH split over the chosen centroid axis.
    int best_bin = -1;
    double sah_cost = std::numeric_limits<double>::infinity();
    if (hi > lo && area > 0.0) {
      std::array<int, kBins> counts{};
      std::array<Box3d, kBins> bins{};
      for (auto& b : bins) b = Box3d::empty();
      const double scale = static_cast<double>(kBins) / (hi - lo);
      for (std::uint32_t i = begin; i < end; ++i) {
        const int b = bin_of(prims_[order_[i]].centroid[axis], lo, scale);
        ++counts[static_cast<std::size_t>(b)];
        bins[static_cast<std::size_t>(b)].expand(prims_[order_[i]].box);
      }
      std::array<double, kBins> left_area{};
      std::array<double, kBins> right_area{};
      std::array<int, kBins> left_count{};
      std::array<int, kBins> right_count{};
      Box3d acc = Box3d::empty();
      int nc = 0;
      for (int b = 0; b < kBins; ++b) {
        acc.expand(bins[static_cast<std::size_t>(b)]);
        nc += counts[static_cast<std::size_t>(b)];
        left_count[static_cast<std::size_t>(b)] = nc;
        left_area[static_cast<std::size_t>(b)] = surface_area(acc);
      }
      acc = Box3d::empty();
      nc = 0;
      for (int b = kBins - 1; b >= 0; --b) {
        acc.expand(bins[static_cast<std::size_t>(b)]);
        nc += counts[static_cast<std::size_t>(b)];
        right_count[static_cast<std::size_t>(b)] = nc;
        right_area[static_cast<std::size_t>(b)] = surface_area(acc);
      }
      for (int b = 0; b + 1 < kBins; ++b) {
        const int lc = left_count[static_cast<std::size_t>(b)];
        const int rc = right_count[static_cast<std::size_t>(b + 1)];
        if (lc == 0 || rc == 0) continue;
        const double cost = kTraversalCost + (left_area[static_cast<std::size_t>(b)] / area) * lc +
                            (right_area[static_cast<std::size_t>(b + 1)] / area) * rc;
        if (cost < sah_cost) {
          sah_cost = cost;
          best_bin = b;
        }
      }
    }

    // Candidate 2: median split by centroid, deterministic via (axis, index).
    std::sort(order_.begin() + begin, order_.begin() + end,
              [&](std::uint32_t x, std::uint32_t y) {
                const double cx = prims_[x].centroid[axis];
                const double cy = prims_[y].centroid[axis];
                if (cx != cy) return cx < cy;
                return x < y;
              });
    const std::uint32_t median_mid = begin + count / 2;
    double median_cost = std::numeric_limits<double>::infinity();
    if (area > 0.0) {
      Box3d lb = Box3d::empty();
      Box3d rb = Box3d::empty();
      for (std::uint32_t i = begin; i < median_mid; ++i) lb.expand(prims_[order_[i]].box);
      for (std::uint32_t i = median_mid; i < end; ++i) rb.expand(prims_[order_[i]].box);
      median_cost = kTraversalCost + (surface_area(lb) / area) * (median_mid - begin) +
                    (surface_area(rb) / area) * (end - median_mid);
    }

    std::uint32_t mid = median_mid;
    if (best_bin >= 0 && sah_cost < median_cost) {
      const double scale = static_cast<double>(kBins) / (hi - lo);
      auto it = std::partition(order_.begin() + begin, order_.begin() + end,
                               [&](std::uint32_t i) {
                                 return bin_of(prims_[i].centroid[axis], lo, scale) <= best_bin;
                               });
      mid = static_cast<std::uint32_t>(it - order_.begin());
    }
    if (mid == begin || mid == end) mid = begin + count / 2;

    const std::int32_t left = build_node(begin, mid, depth + 1);
    const std::int32_t right = build_node(mid, end, depth + 1);
    nodes_[static_cast<std::size_t>(node_index)].left = left;
    nodes_[static_cast<std::size_t>(node_index)].right = right;
    return node_index;
  }

  // 1 / d, mapping a zero component to infinity. A non-finite inverse tells
  // the slab test that the ray is parallel to that axis.
  static double inverse_dir(double d) noexcept {
    if (d == 0.0) return std::numeric_limits<double>::infinity();
    return 1.0 / d;
  }

  static int bin_of(double c, double lo, double scale) noexcept {
    int b = static_cast<int>((c - lo) * scale);
    if (b < 0) b = 0;
    if (b >= kBins) b = kBins - 1;
    return b;
  }

  static double surface_area(const Box3d& b) noexcept {
    if (b.is_empty()) return 0.0;
    const Vec3d d = b.diagonal();
    return 2.0 * (d[0] * d[1] + d[1] * d[2] + d[2] * d[0]);
  }

  static double box_distance2(const Box3d& b, Point3d p) noexcept {
    double sum = 0.0;
    for (int i = 0; i < 3; ++i) {
      double d = 0.0;
      if (p[i] < b.lo[i]) d = b.lo[i] - p[i];
      if (p[i] > b.hi[i]) d = p[i] - b.hi[i];
      sum += d * d;
    }
    return sum;
  }

  // Slab test. A zero direction component makes 1/dir infinite; that case is
  // handled explicitly so a ray starting exactly on a slab plane cannot form
  // 0 * inf = NaN and poison the interval comparisons.
  static bool slab(const Box3d& b, const Ray3d& r, const Vec3d& inv, double tmax) noexcept {
    double tmin = 0.0;
    double t_limit = tmax;
    for (int i = 0; i < 3; ++i) {
      if (!std::isfinite(inv[i])) {
        if (r.origin[i] < b.lo[i] || r.origin[i] > b.hi[i]) return false;
        continue;
      }
      double t0 = (b.lo[i] - r.origin[i]) * inv[i];
      double t1 = (b.hi[i] - r.origin[i]) * inv[i];
      if (t0 > t1) std::swap(t0, t1);
      tmin = std::max(tmin, t0);
      t_limit = std::min(t_limit, t1);
      if (tmin > t_limit) return false;
    }
    return true;
  }

  // Möller-Trumbore. The 1e-14 determinant cutoff is a deliberate
  // floating-point guard (parallel rays and sliver triangles), not an exact
  // predicate; the BVH is an index, not a decision structure.
  static bool intersect_tri(const Ray3d& r, const Triangle3d& t, double& tout) noexcept {
    const Vec3d e1 = t.b - t.a;
    const Vec3d e2 = t.c - t.a;
    const Vec3d pv = cross(r.dir, e2);
    const double det = dot(e1, pv);
    if (std::abs(det) < 1e-14) return false;
    const double inv = 1.0 / det;
    const Vec3d tv{r.origin[0] - t.a[0], r.origin[1] - t.a[1], r.origin[2] - t.a[2]};
    const double u = dot(tv, pv) * inv;
    if (u < 0.0 || u > 1.0) return false;
    const Vec3d qv = cross(tv, e1);
    const double v = dot(r.dir, qv) * inv;
    if (v < 0.0 || u + v > 1.0) return false;
    const double tt = dot(e2, qv) * inv;
    if (tt <= 0.0) return false;
    tout = tt;
    return true;
  }

  std::vector<Prim> prims_;
  std::vector<std::uint32_t> order_;
  std::vector<Node> nodes_;
};

}  // namespace v1
}  // namespace geom
