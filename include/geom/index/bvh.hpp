#pragma once

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <future>
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

// Squared distance from a point to the closest point of a segment. A
// zero-length segment degenerates to the point itself.
inline double squared_distance_segment(Point3d p, Point3d a, Point3d b) noexcept {
  const Vec3d ab = b - a;
  const double len2 = dot(ab, ab);
  if (len2 == 0.0) {
    const Vec3d d = p - a;
    return dot(d, d);
  }
  const Vec3d ap = p - a;
  double t = dot(ap, ab) / len2;
  t = std::clamp(t, 0.0, 1.0);
  const Vec3d q = a + ab * t - p;
  return dot(q, q);
}

// Squared distance from a point to the closest point on a triangle. Uses
// floating point, not the exact predicates: it is a distance query, not a
// degenerate-input decision.
//
// A degenerate triangle (repeated vertex, or exactly collinear) has no plane
// to project onto, so the result is the minimum over its three edges, which is
// zero-length-safe. A non-degenerate triangle projects the point onto its
// plane and uses the plane distance when the projection lies inside the
// triangle (barycentric u, v >= 0 and u + v <= 1); otherwise the closest point
// lies on the boundary and the minimum edge distance is used. The result is
// always finite for finite inputs and finite squared intermediates.
inline double squared_distance(Point3d p, const Triangle3d& t) noexcept {
  const Vec3d e1 = t.b - t.a;
  const Vec3d e2 = t.c - t.a;
  const Vec3d n = cross(e1, e2);
  // Exactly degenerate: no plane. Edge distances also cover a == b or a == c.
  const auto edge_min = [&p, &t]() noexcept {
    return std::min({squared_distance_segment(p, t.a, t.b), squared_distance_segment(p, t.b, t.c),
                     squared_distance_segment(p, t.c, t.a)});
  };
  if (n[0] == 0.0 && n[1] == 0.0 && n[2] == 0.0) return edge_min();
  const Vec3d ap = p - t.a;
  const double nn = dot(n, n);
  // A non-degenerate triangle has nn > 0, but the squared cross can underflow
  // or overflow for extreme coordinates; fall back to the edges then.
  if (!(nn > 0.0) || !std::isfinite(nn)) return edge_min();
  const double s = dot(ap, n) / nn;
  const Point3d q{p[0] - n[0] * s, p[1] - n[1] * s, p[2] - n[2] * s};
  const Vec3d aq = q - t.a;
  const double d00 = dot(e1, e1), d01 = dot(e1, e2), d11 = dot(e2, e2);
  const double den = d00 * d11 - d01 * d01;
  if (!(den > 0.0) || !std::isfinite(den)) return edge_min();
  const double d20 = dot(aq, e1), d21 = dot(aq, e2);
  const double u = (d11 * d20 - d01 * d21) / den;
  const double v = (d00 * d21 - d01 * d20) / den;
  if (u >= 0.0 && v >= 0.0 && u + v <= 1.0) {
    const Vec3d d = p - q;
    return dot(d, d);
  }
  return edge_min();
}

// Bounding volume hierarchy over triangles, built with a 12-bin surface area
// heuristic and queried iteratively (no recursion on the query path).
//
// Memory layout. The primitives are stored structure-of-arrays (triangles and
// boxes in separate contiguous vectors; a primitive's centroid is the midpoint of
// its box, so it needs no array of its own); `order_` is a permutation of
// primitive indices and every leaf covers a contiguous range
// `[first, first + count)` of that permutation. Nodes live in one contiguous
// vector; a node is a leaf when `count != 0`, otherwise it is an interior node
// with `left`/`right` child indices (int32, never negative for a valid node). A
// depth guard bounds the tree at 64 levels, so the explicit query stacks never
// exceed 65 entries.
//
// Build: 12 bins over the centroid bounds of the current range, cost
// `traversal + (left_area / area) * nleft + (right_area / area) * nright`
// (the traversal constant cancels between candidates but is kept so the value
// is comparable to a classic SAH cost). If the best bin split is not strictly
// cheaper than a median split by centroid, or would leave one side empty, the
// median split is used instead; it keeps the `count / 2` smallest (centroid
// axis, index) keys, so the tree is deterministic. Both candidates reuse the
// one binned scan: the prefix and suffix unions of the per-bin boxes are the
// exact child bounds of every candidate, and because bins are ordered by
// centroid the median split is a prefix of the bin order, so only the single bin
// straddling the median needs a selection pass. A node's box and centroid box
// are handed down by its parent (the root unions its primitives once), so a
// primitive is examined once per level. Leaves hold at most 4 primitives unless
// the depth guard forces a larger one.
//
// Parallel build: above a primitive threshold the top levels are built on the
// calling thread and the frontier subtrees at a fixed depth are built by
// independent tasks, each into its own node array and over its own slice of
// `order_`. The arrays are then appended in frontier order. Both the frontier
// depth and the assembly order are fixed, so the node array and the permutation
// are byte-identical whatever the thread count or completion order.
//
// `bounds()` is the root node's box, therefore exactly the union of all
// primitive boxes. An empty tree has no node at all and returns an empty box.
//
// Ray contract: `dir` is expected to be finite and non-zero; t is measured in
// units of |dir|, so a non-normalised direction is fine. The accepted interval
// is (0, tmax]: hits with t <= 0 are behind the origin and are ignored, and a
// hit exactly at t == tmax is accepted. Equal t ties resolve to the smallest
// primitive index.
class Bvh {
 public:
  struct Hit {
    std::uint32_t prim = 0;
    double t = 0.0;
  };

  void build(std::span<const Triangle3d> triangles) { build_impl(triangles, true); }

  const Box3d& bounds() const noexcept {
    static const Box3d kEmpty = Box3d::empty();
    return nodes_.empty() ? kEmpty : nodes_[0].box;
  }

  std::uint32_t size() const noexcept { return static_cast<std::uint32_t>(tris_.size()); }

  bool raycast(const Ray3d& ray, double tmax, Hit& hit) const {
    if (nodes_.empty()) return false;
    const Vec3d inv{inverse_dir(ray.dir[0]), inverse_dir(ray.dir[1]), inverse_dir(ray.dir[2])};
    std::array<std::int32_t, 128> stack;
    int sp = 0;
    stack[static_cast<std::size_t>(sp++)] = 0;
    bool found = false;
    double best = tmax;
    std::uint32_t best_prim = 0;
    while (sp > 0) {
      const std::size_t level = static_cast<std::size_t>(--sp);
      const Node& n = nodes_[static_cast<std::size_t>(stack[level])];
      if (!slab(n.box, ray, inv, best)) continue;
      if (n.count > 0) {
        for (std::uint32_t i = n.first; i < n.first + n.count; ++i) {
          const std::uint32_t prim = order_[i];
          double t = 0.0;
          if (!intersect_tri(ray, tris_[prim], t)) continue;
          // The interval is (0, tmax]: a hit exactly at tmax is included, so the
          // slab pruning limit above also keeps equality.
          const bool better = !found || t < best || (t == best && prim < best_prim);
          if (t <= best && better) {
            best = t;
            best_prim = prim;
            found = true;
          }
        }
      } else {
        stack[static_cast<std::size_t>(sp++)] = n.left;
        stack[static_cast<std::size_t>(sp++)] = n.right;
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
  // returns 0xFFFFFFFF and an infinite distance. The squared intermediates
  // (dx*dx + dy*dy + dz*dz) must stay finite: with coordinates around 1e200 the
  // sum overflows to infinity and any distance involving a NaN coordinate is
  // NaN, which is outside the supported coordinate range.
  std::uint32_t nearest(const Point3d& p, double& best_d2) const {
    if (nodes_.empty()) {
      best_d2 = std::numeric_limits<double>::infinity();
      return 0xFFFFFFFFu;
    }
    std::array<std::int32_t, 128> stack;
    int sp = 0;
    stack[static_cast<std::size_t>(sp++)] = 0;
    best_d2 = std::numeric_limits<double>::infinity();
    std::uint32_t best = 0xFFFFFFFFu;
    while (sp > 0) {
      const std::size_t level = static_cast<std::size_t>(--sp);
      const Node& n = nodes_[static_cast<std::size_t>(stack[level])];
      // Equality is not pruned: a different primitive at exactly best_d2 may
      // still win the index tie-break.
      if (box_distance2(n.box, p) > best_d2) continue;
      if (n.count > 0) {
        for (std::uint32_t i = n.first; i < n.first + n.count; ++i) {
          const std::uint32_t prim = order_[i];
          const double d2 = squared_distance(p, tris_[prim]);
          if (d2 < best_d2 || (d2 == best_d2 && prim < best)) {
            best_d2 = d2;
            best = prim;
          }
        }
      } else {
        stack[static_cast<std::size_t>(sp++)] = n.left;
        stack[static_cast<std::size_t>(sp++)] = n.right;
      }
    }
    return best;
  }

 private:
  // Test-only accessor for the determinism digest; it is not part of the
  // public API and does not expose the node layout beyond the tests.
  friend struct BvhTestAccess;

  // `allow_parallel` is a test seam: the determinism test builds one input both
  // ways to check that the two modes agree on the tree shape and on query
  // results. Mode state is reset here, so a build that throws part-way through
  // the parallel phase cannot leak into the next one.
  void build_impl(std::span<const Triangle3d> triangles, bool allow_parallel) {
    pending_.clear();
    cut_depth_ = -1;
    tris_.assign(triangles.begin(), triangles.end());
    boxes_.resize(triangles.size());
    order_.clear();
    nodes_.clear();
    order_.reserve(triangles.size());
    Box3d root_box = Box3d::empty();
    Box3d root_centroid_box = Box3d::empty();
    for (std::size_t i = 0; i < triangles.size(); ++i) {
      const Box3d box = triangle_box(triangles[i]);
      boxes_[i] = box;
      root_box.expand(box);
      root_centroid_box.expand(center_of(box));
      order_.push_back(static_cast<std::uint32_t>(i));
    }
    if (tris_.empty()) return;
    scratch_.resize(tris_.size());
    bin_index_.resize(tris_.size());
    nodes_.reserve(2 * tris_.size());
    const auto n = static_cast<std::uint32_t>(tris_.size());
    if (allow_parallel && n >= kParallelMinPrims) {
      build_parallel(n, root_box, root_centroid_box);
    } else {
      build_node(nodes_, 0, n, 0, root_box, root_centroid_box);
    }
  }

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
  // Subtrees at this depth are built independently. The depth is fixed, so the
  // resulting node array does not depend on the number of hardware threads or on
  // the order in which the subtrees finish.
  static constexpr int kParallelDepth = 3;
  static constexpr std::uint32_t kParallelMinPrims = 1u << 15;

  // A subtree handed to the parallel phase: its range, the boxes its parent
  // already computed, and the index of its (already emitted) root node.
  struct Pending {
    std::uint32_t begin;
    std::uint32_t end;
    int depth;
    Box3d box;
    Box3d centroid_box;
    std::int32_t node;
  };

  // Builds the subtree over order_[begin, end) into `out`, returning the node's
  // index in that vector. `box` and `centroid_box` are the exact unions of the
  // range, computed by the caller from the data it already touched for its own
  // split. When `cut_depth_` is set, a node at that depth is emitted without its
  // children and recorded in `pending_` instead.
  std::int32_t build_node(std::vector<Node>& out, std::uint32_t begin, std::uint32_t end, int depth,
                          const Box3d& box, const Box3d& centroid_box) {
    const auto node_index = static_cast<std::int32_t>(out.size());
    out.push_back(Node{});
    out.back().box = box;

    const std::uint32_t count = end - begin;
    if (count <= kMaxLeaf || depth >= kMaxDepth) {
      out[static_cast<std::size_t>(node_index)].first = begin;
      out[static_cast<std::size_t>(node_index)].count = count;
      return node_index;
    }
    if (depth == cut_depth_) {
      pending_.push_back(Pending{begin, end, depth, box, centroid_box, node_index});
      return node_index;
    }

    const Vec3d extent = centroid_box.diagonal();
    int axis = 0;
    if (extent[1] > extent[axis]) axis = 1;
    if (extent[2] > extent[axis]) axis = 2;
    const double lo = centroid_box.lo[axis];
    const double hi = centroid_box.hi[axis];
    const double area = surface_area(box);

    // Candidate 1: 12-bin SAH split over the chosen centroid axis. The same
    // scan also records each primitive's bin, so grouping by bin is a counting
    // sort instead of a second predicate pass, and the per-bin boxes become the
    // exact child bounds.
    const std::uint32_t median_mid = begin + count / 2;
    std::uint32_t mid = median_mid;
    Box3d left_box = Box3d::empty();
    Box3d left_cbox = Box3d::empty();
    Box3d right_box = Box3d::empty();
    Box3d right_cbox = Box3d::empty();
    if (hi > lo && area > 0.0) {
      std::array<Box3d, kBins> bin_boxes;
      std::array<Box3d, kBins> bin_cboxes;
      for (auto& b : bin_boxes) b = Box3d::empty();
      for (auto& b : bin_cboxes) b = Box3d::empty();
      std::array<int, kBins> counts{};
      std::array<int, kBins> prefix_count{};
      std::array<int, kBins> suffix_count{};
      std::array<double, kBins> prefix_area{};
      std::array<double, kBins> suffix_area{};
      // The prefix and suffix unions are kept as boxes, not just as areas, so a
      // candidate's child bounds are a copy instead of a second union of bins.
      std::array<Box3d, kBins> prefix_box;
      std::array<Box3d, kBins> prefix_cbox;
      std::array<Box3d, kBins> suffix_box;
      std::array<Box3d, kBins> suffix_cbox;
      const double scale = static_cast<double>(kBins) / (hi - lo);
      for (std::uint32_t i = begin; i < end; ++i) {
        const std::uint32_t prim = order_[i];
        const Box3d& prim_box = boxes_[prim];
        const Point3d centroid = center_of(prim_box);
        const int b = bin_of(centroid[axis], lo, scale);
        ++counts[static_cast<std::size_t>(b)];
        bin_boxes[static_cast<std::size_t>(b)].expand(prim_box);
        bin_cboxes[static_cast<std::size_t>(b)].expand(centroid);
        bin_index_[i] = static_cast<std::uint8_t>(b);
      }
      Box3d acc = Box3d::empty();
      Box3d acc_c = Box3d::empty();
      int nc = 0;
      for (int b = 0; b < kBins; ++b) {
        const auto index = static_cast<std::size_t>(b);
        acc.expand(bin_boxes[index]);
        acc_c.expand(bin_cboxes[index]);
        nc += counts[index];
        prefix_count[index] = nc;
        prefix_box[index] = acc;
        prefix_cbox[index] = acc_c;
        prefix_area[index] = surface_area(acc);
      }
      acc = Box3d::empty();
      acc_c = Box3d::empty();
      nc = 0;
      for (int b = kBins - 1; b >= 0; --b) {
        const auto index = static_cast<std::size_t>(b);
        acc.expand(bin_boxes[index]);
        acc_c.expand(bin_cboxes[index]);
        nc += counts[index];
        suffix_count[index] = nc;
        suffix_box[index] = acc;
        suffix_cbox[index] = acc_c;
        suffix_area[index] = surface_area(acc);
      }
      int best_bin = -1;
      double sah_cost = std::numeric_limits<double>::infinity();
      for (int b = 0; b + 1 < kBins; ++b) {
        const int lc = prefix_count[static_cast<std::size_t>(b)];
        const int rc = suffix_count[static_cast<std::size_t>(b + 1)];
        if (lc == 0 || rc == 0) continue;
        const double cost =
            kTraversalCost + (prefix_area[static_cast<std::size_t>(b)] / area) * lc +
            (suffix_area[static_cast<std::size_t>(b + 1)] / area) * rc;
        if (cost < sah_cost) {
          sah_cost = cost;
          best_bin = b;
        }
      }

      // Candidate 2's split point: the median keeps the count / 2 smallest
      // (centroid, index) keys. Bins are monotone in the centroid, so that set
      // is every bin below the straddling bin plus its head.
      int median_bin = 0;
      while (prefix_count[static_cast<std::size_t>(median_bin)] < static_cast<int>(count / 2)) {
        ++median_bin;
      }
      const std::uint32_t block_begin =
          begin + (median_bin == 0
                       ? 0u
                       : static_cast<std::uint32_t>(
                             prefix_count[static_cast<std::size_t>(median_bin) - 1]));
      const std::uint32_t block_end =
          begin + static_cast<std::uint32_t>(prefix_count[static_cast<std::size_t>(median_bin)]);

      // Group the range by bin (counting sort through `scratch_`). This puts
      // the straddling bin in one contiguous block, so the median selection
      // only touches that block, and it makes the SAH split a plain index.
      std::array<std::uint32_t, kBins> cursor{};
      std::uint32_t pos = begin;
      for (int b = 0; b < kBins; ++b) {
        cursor[static_cast<std::size_t>(b)] = pos;
        pos += static_cast<std::uint32_t>(counts[static_cast<std::size_t>(b)]);
      }
      for (std::uint32_t i = begin; i < end; ++i) {
        scratch_[cursor[bin_index_[i]]++] = order_[i];
      }
      std::copy_n(scratch_.begin() + begin, count, order_.begin() + begin);

      // Candidate 2: median split by centroid, deterministic via (axis, index).
      // Only the straddling bin needs selecting: everything before it is
      // already smaller and everything after it is already larger.
      if (median_mid < block_end) {
        std::nth_element(order_.begin() + block_begin, order_.begin() + median_mid,
                         order_.begin() + block_end, [&](std::uint32_t x, std::uint32_t y) {
                           const double cx = center_of(boxes_[x])[axis];
                           const double cy = center_of(boxes_[y])[axis];
                           if (cx != cy) return cx < cy;
                           return x < y;
                         });
      }

      if (median_bin > 0) {
        left_box = prefix_box[static_cast<std::size_t>(median_bin) - 1];
        left_cbox = prefix_cbox[static_cast<std::size_t>(median_bin) - 1];
      }
      if (median_bin + 1 < kBins) {
        right_box = suffix_box[static_cast<std::size_t>(median_bin) + 1];
        right_cbox = suffix_cbox[static_cast<std::size_t>(median_bin) + 1];
      }
      for (std::uint32_t i = block_begin; i < block_end; ++i) {
        const std::uint32_t prim = order_[i];
        if (i < median_mid) {
          left_box.expand(boxes_[prim]);
          left_cbox.expand(center_of(boxes_[prim]));
        } else {
          right_box.expand(boxes_[prim]);
          right_cbox.expand(center_of(boxes_[prim]));
        }
      }

      const double median_cost =
          kTraversalCost + (surface_area(left_box) / area) * (median_mid - begin) +
          (surface_area(right_box) / area) * (end - median_mid);
      if (best_bin >= 0 && sah_cost < median_cost) {
        mid = begin + static_cast<std::uint32_t>(prefix_count[static_cast<std::size_t>(best_bin)]);
        left_box = prefix_box[static_cast<std::size_t>(best_bin)];
        left_cbox = prefix_cbox[static_cast<std::size_t>(best_bin)];
        right_box = suffix_box[static_cast<std::size_t>(best_bin) + 1];
        right_cbox = suffix_cbox[static_cast<std::size_t>(best_bin) + 1];
      }
    } else {
      // No usable axis range, or a zero-area parent: keep the whole-range
      // median selection and derive both children by scanning.
      if (median_mid < end) {
        std::nth_element(order_.begin() + begin, order_.begin() + median_mid, order_.begin() + end,
                         [&](std::uint32_t x, std::uint32_t y) {
                           const double cx = center_of(boxes_[x])[axis];
                           const double cy = center_of(boxes_[y])[axis];
                           if (cx != cy) return cx < cy;
                           return x < y;
                         });
      }
      for (std::uint32_t i = begin; i < median_mid; ++i) {
        left_box.expand(boxes_[order_[i]]);
        left_cbox.expand(center_of(boxes_[order_[i]]));
      }
      for (std::uint32_t i = median_mid; i < end; ++i) {
        right_box.expand(boxes_[order_[i]]);
        right_cbox.expand(center_of(boxes_[order_[i]]));
      }
    }
    if (mid == begin || mid == end) {
      // The candidates above never leave a side empty; this keeps a malformed
      // cost (e.g. a non-finite area) from producing one anyway.
      mid = median_mid;
      left_box = Box3d::empty();
      left_cbox = Box3d::empty();
      right_box = Box3d::empty();
      right_cbox = Box3d::empty();
      for (std::uint32_t i = begin; i < mid; ++i) {
        left_box.expand(boxes_[order_[i]]);
        left_cbox.expand(center_of(boxes_[order_[i]]));
      }
      for (std::uint32_t i = mid; i < end; ++i) {
        right_box.expand(boxes_[order_[i]]);
        right_cbox.expand(center_of(boxes_[order_[i]]));
      }
    }

    const std::int32_t left = build_node(out, begin, mid, depth + 1, left_box, left_cbox);
    const std::int32_t right = build_node(out, mid, end, depth + 1, right_box, right_cbox);
    out[static_cast<std::size_t>(node_index)].left = left;
    out[static_cast<std::size_t>(node_index)].right = right;
    return node_index;
  }

  // Builds the top kParallelDepth levels on this thread, then the frontier
  // subtrees on separate tasks. Each subtree writes its own node array and its
  // own slice of `order_`, so the work is independent; the arrays are then
  // appended in frontier order. The bytes produced are therefore the same for
  // any thread count and any completion order.
  void build_parallel(std::uint32_t n, const Box3d& root_box, const Box3d& root_centroid_box) {
    pending_.clear();
    cut_depth_ = kParallelDepth;
    build_node(nodes_, 0, n, 0, root_box, root_centroid_box);
    cut_depth_ = -1;
    if (pending_.empty()) {
      // The top phase already finished the whole tree.
      return;
    }
    // A single frontier subtree still goes through the append path below: its
    // root node was emitted by the top phase and still needs its children.
    std::vector<std::vector<Node>> subtrees(pending_.size());
    try {
      std::vector<std::future<void>> futures;
      futures.reserve(pending_.size());
      for (std::size_t t = 0; t < pending_.size(); ++t) {
        const Pending task = pending_[t];
        futures.push_back(std::async(std::launch::async, [this, &subtrees, t, task] {
          build_node(subtrees[t], task.begin, task.end, task.depth, task.box, task.centroid_box);
        }));
      }
      for (auto& future : futures) future.get();
    } catch (...) {
      // A thread or an allocation failed: rebuild on this thread. The tree does
      // not depend on how the subtrees were scheduled, so the result is the same.
      subtrees.clear();
      nodes_.clear();
      pending_.clear();
      build_node(nodes_, 0, n, 0, root_box, root_centroid_box);
      return;
    }
    for (std::size_t t = 0; t < pending_.size(); ++t) {
      const std::vector<Node>& subtree = subtrees[t];
      const auto base = static_cast<std::int32_t>(nodes_.size());
      // The subtree's root is the node already emitted by the top phase, so its
      // children are packed from `base` on: local index j >= 1 maps to base+j-1.
      const std::int32_t shift = base - 1;
      for (std::size_t j = 1; j < subtree.size(); ++j) {
        Node node = subtree[j];
        if (node.left >= 0) node.left += shift;
        if (node.right >= 0) node.right += shift;
        nodes_.push_back(node);
      }
      const Node& root = subtree.front();
      nodes_[static_cast<std::size_t>(pending_[t].node)].left = root.left + shift;
      nodes_[static_cast<std::size_t>(pending_[t].node)].right = root.right + shift;
    }
    pending_.clear();
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

  // A primitive's centroid is the midpoint of its box, so it never needs its
  // own array: the build derives it from the only field it has to read anyway.
  static Point3d center_of(const Box3d& b) noexcept {
    return Point3d{(b.lo[0] + b.hi[0]) * 0.5, (b.lo[1] + b.hi[1]) * 0.5,
                   (b.lo[2] + b.hi[2]) * 0.5};
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

  // Möller-Trumbore. The determinant guard is relative:
  // |det| <= 1e-14 * |e1| * |e2| * |dir| rejects only genuinely near-parallel
  // rays and sliver triangles, independent of the coordinate scale of the
  // triangle. It is a floating-point guard, not an exact predicate; the product
  // is computed in double, so it saturates once |e1| * |e2| * |dir| exceeds
  // DBL_MAX (~1e308). Past that point the guard rejects every ray rather than
  // misclassifying one, and the sliver measure itself has already overflowed.
  static bool intersect_tri(const Ray3d& r, const Triangle3d& t, double& tout) noexcept {
    const Vec3d e1 = t.b - t.a;
    const Vec3d e2 = t.c - t.a;
    const Vec3d pv = cross(r.dir, e2);
    const double det = dot(e1, pv);
    const double scale =
        std::sqrt(dot(e1, e1)) * std::sqrt(dot(e2, e2)) * std::sqrt(dot(r.dir, r.dir));
    if (!(std::abs(det) > 1e-14 * scale)) return false;
    const Vec3d tv{r.origin[0] - t.a[0], r.origin[1] - t.a[1], r.origin[2] - t.a[2]};
    // Divide by det instead of multiplying by 1/det: the reciprocal is
    // rounded once, so num * (1/det) can land one ulp off even when num == det
    // (t = 1). Division rounds the quotient directly and returns exactly 1.
    const double u = dot(tv, pv) / det;
    if (u < 0.0 || u > 1.0) return false;
    const Vec3d qv = cross(tv, e1);
    const double v = dot(r.dir, qv) / det;
    if (v < 0.0 || u + v > 1.0) return false;
    const double tt = dot(e2, qv) / det;
    if (tt <= 0.0) return false;
    tout = tt;
    return true;
  }

  // Structure-of-arrays primitive storage: the build touches only boxes and
  // centroids, and the binning scan reads the centroid axis alone, so keeping
  // each field in its own contiguous array cuts the random-access footprint of
  // a level by roughly the size of the triangle data.
  std::vector<Triangle3d> tris_;
  std::vector<Box3d> boxes_;
  std::vector<std::uint32_t> order_;
  std::vector<Node> nodes_;
  // Build scratch, sized to the primitive count: the binned range grouped into
  // `scratch_`, and the bin of each order_ slot for the counting sort.
  std::vector<std::uint32_t> scratch_;
  std::vector<std::uint8_t> bin_index_;
  // Frontier subtrees recorded by the top phase of a parallel build.
  std::vector<Pending> pending_;
  int cut_depth_ = -1;
};

}  // namespace v1
}  // namespace geom
