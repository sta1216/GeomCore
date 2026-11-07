#pragma once

#include <algorithm>
#include <array>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <span>
#include <vector>

#include <geom/primitives/box.hpp>
#include <geom/primitives/point.hpp>
#include <geom/primitives/vec.hpp>

namespace geom {
inline namespace v1 {

// Static KD tree over 3D points, used for nearest-neighbour and box queries.
//
// Storage. `points_` is a copy of the build input, `order_` is a permutation of
// point indices and every leaf covers a contiguous range
// `[first, first + count)` of that permutation. Nodes live in one contiguous
// vector; a node is a leaf when `count > 0`, otherwise it is an interior node
// with `left`/`right` child indices (int32, never negative for a valid node).
//
// Build. Each node's box is the union of the points in its range. A range of at
// most 8 points, or a depth of 64, becomes a leaf. Otherwise the split axis is
// the longest side of the node box (ties resolve x, then y, then z) and the
// range is split at its median with `nth_element` over the key
// `(coordinate[axis], index)`, which makes the tree deterministic for duplicate
// coordinates. Rebuilding clears the previous tree.
//
// Preconditions: coordinates must be finite. NaN coordinates make the median
// comparator non-transitive and can leave `nearest` returning the 0xFFFFFFFF
// sentinel with an infinite distance. Indices are `uint32_t`, so `build`
// supports at most `UINT32_MAX` points; a larger input is silently truncated.
//
// Queries use an explicit stack (the depth guard bounds it at 65 entries) and
// prune by the squared distance from the query point to a node box. Ties in
// `nearest` and `nearest_k` resolve to the smallest point index: a node whose
// box is exactly at the current k-th distance is still visited, because a
// smaller index may win the tie there. Coordinates are expected to be finite;
// squared distances (dx*dx + dy*dy + dz*dz) must stay finite, so coordinates
// around 1e200 overflow and NaN coordinates are outside the supported range.
//
// Empty-tree contract: `size()` is 0, `bounds()` is an empty box,
// `nearest` returns 0xFFFFFFFF with `best_d2 = +inf`, `nearest_k` writes
// nothing, and `range` returns an empty vector.
class KdTree {
 public:
  // Replaces the tree with one over `points`. The input is copied, so the span
  // may be released after this call.
  void build(std::span<const Point3d> points) {
    points_.assign(points.begin(), points.end());
    order_.clear();
    nodes_.clear();
    if (points_.empty()) return;
    order_.reserve(points_.size());
    for (std::size_t i = 0; i < points_.size(); ++i) {
      order_.push_back(static_cast<std::uint32_t>(i));
    }
    nodes_.reserve(2 * points_.size());
    build_node(0, static_cast<std::uint32_t>(points_.size()), 0);
  }

  std::uint32_t size() const noexcept { return static_cast<std::uint32_t>(points_.size()); }

  // Union of all point boxes; an empty box when the tree is empty.
  const Box3d& bounds() const noexcept {
    static const Box3d kEmpty = Box3d::empty();
    return nodes_.empty() ? kEmpty : nodes_[0].box;
  }

  // Returns the index of the nearest point and writes its squared distance;
  // equal distances resolve to the smallest index. An empty tree returns
  // 0xFFFFFFFF and writes +inf.
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
      // Equality is not pruned: another point at exactly best_d2 may still win
      // the index tie-break.
      if (box_distance2(n.box, p) > best_d2) continue;
      if (n.count > 0) {
        for (std::uint32_t i = n.first; i < n.first + n.count; ++i) {
          const std::uint32_t index = order_[i];
          const double d2 = point_distance2(p, points_[index]);
          if (d2 < best_d2 || (d2 == best_d2 && index < best)) {
            best_d2 = d2;
            best = index;
          }
        }
      } else {
        stack[sp++] = n.left;
        stack[sp++] = n.right;
      }
    }
    return best;
  }

  // Writes at most min(k, out.size()) indices, ordered by (distance, index)
  // ascending, and returns the number written. An empty tree or k == 0 writes
  // nothing. The k nearest set is correct even when several points share a
  // distance, because pruning stops only when the box distance is strictly
  // greater than the current worst entry.
  std::uint32_t nearest_k(const Point3d& p, std::uint32_t k,
                          std::span<std::uint32_t> out) const {
    const auto limit =
        static_cast<std::uint32_t>(std::min<std::size_t>(k, out.size()));
    if (limit == 0 || nodes_.empty()) return 0;
    std::vector<Entry> heap;
    heap.reserve(limit);
    std::array<std::int32_t, 128> stack{};
    int sp = 0;
    stack[sp++] = 0;
    while (sp > 0) {
      const Node& n = nodes_[static_cast<std::size_t>(stack[--sp])];
      if (heap.size() >= limit && box_distance2(n.box, p) > heap.front().d2) continue;
      if (n.count > 0) {
        for (std::uint32_t i = n.first; i < n.first + n.count; ++i) {
          const Entry candidate{point_distance2(p, points_[order_[i]]), order_[i]};
          if (heap.size() < limit) {
            heap.push_back(candidate);
            std::push_heap(heap.begin(), heap.end(), better);
          } else if (better(candidate, heap.front())) {
            std::pop_heap(heap.begin(), heap.end(), better);
            heap.back() = candidate;
            std::push_heap(heap.begin(), heap.end(), better);
          }
        }
      } else {
        stack[sp++] = n.left;
        stack[sp++] = n.right;
      }
    }
    std::sort(heap.begin(), heap.end(), better);
    for (std::size_t i = 0; i < heap.size(); ++i) out[i] = heap[i].index;
    return static_cast<std::uint32_t>(heap.size());
  }

  // Indices of the points inside `box`, boundary included, in ascending order.
  // An empty tree or an empty box returns an empty vector. The cost is the
  // number of visited nodes plus the output size, with a final O(m log m) sort
  // that puts the m collected indices in ascending order.
  std::vector<std::uint32_t> range(const Box3d& box) const {
    std::vector<std::uint32_t> hits;
    if (nodes_.empty() || box.is_empty()) return hits;
    std::array<std::int32_t, 128> stack{};
    int sp = 0;
    stack[sp++] = 0;
    while (sp > 0) {
      const Node& n = nodes_[static_cast<std::size_t>(stack[--sp])];
      if (!intersects(n.box, box)) continue;
      if (n.count > 0) {
        for (std::uint32_t i = n.first; i < n.first + n.count; ++i) {
          const std::uint32_t index = order_[i];
          if (box.contains(points_[index])) hits.push_back(index);
        }
      } else {
        stack[sp++] = n.left;
        stack[sp++] = n.right;
      }
    }
    std::sort(hits.begin(), hits.end());
    return hits;
  }

 private:
  struct Node {
    Box3d box;
    std::int32_t left = -1;    // interior child index
    std::int32_t right = -1;   // interior child index
    std::uint32_t first = 0;   // leaf: first index into order_
    std::uint32_t count = 0;   // leaf: number of points; 0 => interior
  };

  struct Entry {
    double d2;
    std::uint32_t index;
  };

  static constexpr std::uint32_t kLeafSize = 8;
  static constexpr int kMaxDepth = 64;

  // True when `a` is closer than `b`, with the index resolving exact ties.
  // Used both as the max-heap comparator (which keeps the worst entry on top)
  // and as the final ascending sort.
  static bool better(const Entry& a, const Entry& b) noexcept {
    if (a.d2 != b.d2) return a.d2 < b.d2;
    return a.index < b.index;
  }

  static double point_distance2(const Point3d& a, const Point3d& b) noexcept {
    const double dx = a[0] - b[0];
    const double dy = a[1] - b[1];
    const double dz = a[2] - b[2];
    return dx * dx + dy * dy + dz * dz;
  }

  static double box_distance2(const Box3d& b, const Point3d& p) noexcept {
    double sum = 0.0;
    for (int i = 0; i < 3; ++i) {
      double d = 0.0;
      if (p[i] < b.lo[i]) d = b.lo[i] - p[i];
      if (p[i] > b.hi[i]) d = p[i] - b.hi[i];
      sum += d * d;
    }
    return sum;
  }

  static bool intersects(const Box3d& a, const Box3d& b) noexcept {
    for (int i = 0; i < 3; ++i) {
      if (a.lo[i] > b.hi[i] || a.hi[i] < b.lo[i]) return false;
    }
    return true;
  }

  std::int32_t build_node(std::uint32_t begin, std::uint32_t end, int depth) {
    const auto node_index = static_cast<std::int32_t>(nodes_.size());
    nodes_.push_back(Node{});
    Box3d box = Box3d::empty();
    for (std::uint32_t i = begin; i < end; ++i) box.expand(points_[order_[i]]);
    {
      Node& node = nodes_[static_cast<std::size_t>(node_index)];
      node.box = box;
      const std::uint32_t count = end - begin;
      if (count <= kLeafSize || depth >= kMaxDepth) {
        node.first = begin;
        node.count = count;
        return node_index;
      }
    }

    const Vec3d extent = box.diagonal();
    int axis = 0;
    if (extent[1] > extent[0]) axis = 1;
    if (extent[2] > extent[axis]) axis = 2;
    const std::uint32_t mid = begin + (end - begin) / 2;
    std::nth_element(order_.begin() + begin, order_.begin() + mid, order_.begin() + end,
                     [&](std::uint32_t x, std::uint32_t y) {
                       const double cx = points_[x][axis];
                       const double cy = points_[y][axis];
                       if (cx != cy) return cx < cy;
                       return x < y;
                     });

    const std::int32_t left = build_node(begin, mid, depth + 1);
    const std::int32_t right = build_node(mid, end, depth + 1);
    Node& node = nodes_[static_cast<std::size_t>(node_index)];
    node.left = left;
    node.right = right;
    return node_index;
  }

  std::vector<Point3d> points_;
  std::vector<std::uint32_t> order_;
  std::vector<Node> nodes_;
};

}  // namespace v1
}  // namespace geom
