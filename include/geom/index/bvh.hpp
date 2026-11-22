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

#if defined(_MSC_VER) || defined(__SSE__)
#include <xmmintrin.h>
#endif

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
// its box, so it needs no array of its own). The build permutes the primitives
// into leaf order: `order_[i]` is the original index of the primitive stored at
// slot i, and every leaf covers a contiguous range of slots, so a leaf's
// primitive tests walk contiguous triangle storage. The tree is widened after
// the build: a wide node stores `kFanout` child boxes and child references in
// one record, so a visited traversal level costs one 224-byte node record instead of one
// line per binary node. A child reference is either a leaf range (tagged payload)
// or the index of another wide node. The binary tree produced by the SAH build
// is converted depth-first, so the wide array is also deterministic. A depth
// guard bounds the binary tree at 64 levels, hence the query stack bound below.
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
// `bounds()` is the exact union of all primitive boxes. An empty tree has no
// node at all and returns an empty box.
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
    return tris_.empty() ? kEmpty : root_box_;
  }

  std::uint32_t size() const noexcept { return static_cast<std::uint32_t>(tris_.size()); }

  bool raycast(const Ray3d& ray, double tmax, Hit& hit) const {
    if (tris_.empty() || root_payload_ == kInvalid) return false;
    const Vec3d inv{inverse_dir(ray.dir[0]), inverse_dir(ray.dir[1]), inverse_dir(ray.dir[2])};
    // |dir| is shared by every primitive test of this ray; it used to be
    // recomputed (three square roots) inside each test.
    const double dir_len = std::sqrt(dot(ray.dir, ray.dir));
    std::array<StackEntry, kStackSize> stack;
    int sp = 0;
    bool found = false;
    double best = tmax;
    std::uint32_t best_prim = 0;
    // The root box is tested once; every other child box arrives with its
    // interior node's record, so a visited level costs one 224-byte node record.
    double root_t = 0.0;
    if (!slab(root_box_, ray, inv, best * (1.0 + kPruneSlack), root_t)) return false;
    stack[static_cast<std::size_t>(sp++)] = StackEntry{root_t, root_payload_};
    while (sp > 0) {
      const StackEntry entry = stack[static_cast<std::size_t>(--sp)];
      // A child box entered behind the current best (allowing for the rounding
      // slack described at kPruneSlack) can only contain hits at or beyond it;
      // equality is kept so a tie with a smaller index still wins.
      if (entry.t > best * (1.0 + kPruneSlack)) continue;
      const std::uint64_t payload = entry.payload;
      if (is_leaf(payload)) {
        const std::uint32_t first = leaf_first(payload);
        const std::uint32_t count = leaf_count(payload);
        for (std::uint32_t i = first; i < first + count; ++i) {
          double t = 0.0;
          if (!intersect_tri(ray, tris_[i], dir_len, t)) continue;
          // The interval is (0, tmax]: a hit exactly at tmax is included, so the
          // slab pruning limit above also keeps equality.
          if (t > best) continue;
          const std::uint32_t prim = order_[i];
          const bool better = !found || t < best || (t == best && prim < best_prim);
          if (better) {
            best = t;
            best_prim = prim;
            found = true;
          }
        }
      } else {
        const Node& n = nodes_[static_cast<std::size_t>(payload)];
        for (int j = 0; j < kFanout; ++j) {
          const std::uint64_t child = n.child[j];
          if (child == kInvalid) break;
          double t = 0.0;
          if (!slab(n.lo[j], n.hi[j], ray, inv, best * (1.0 + kPruneSlack), t)) continue;
          if (is_leaf(child)) {
            // A leaf is small: testing it here avoids a stack round trip and
            // lets its hits tighten `best` before the remaining slots of this
            // node are tested. The result is identical to visiting the leaf
            // later, because every primitive is still compared with the same
            // tie-break.
            const std::uint32_t first = leaf_first(child);
            const std::uint32_t count = leaf_count(child);
            for (std::uint32_t i = first; i < first + count; ++i) {
              double tt = 0.0;
              if (!intersect_tri(ray, tris_[i], dir_len, tt)) continue;
              if (tt > best) continue;
              const std::uint32_t prim = order_[i];
              const bool better = !found || tt < best || (tt == best && prim < best_prim);
              if (better) {
                best = tt;
                best_prim = prim;
                found = true;
              }
            }
          } else {
            prefetch(&nodes_[static_cast<std::size_t>(child)]);
            stack[static_cast<std::size_t>(sp++)] = StackEntry{t, child};
          }
        }
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
    if (tris_.empty() || root_payload_ == kInvalid) {
      best_d2 = std::numeric_limits<double>::infinity();
      return 0xFFFFFFFFu;
    }
    std::array<StackEntry, kStackSize> stack;
    int sp = 0;
    best_d2 = std::numeric_limits<double>::infinity();
    std::uint32_t best = 0xFFFFFFFFu;
    stack[static_cast<std::size_t>(sp++)] =
        StackEntry{box_distance2(root_box_, p), root_payload_};
    while (sp > 0) {
      const StackEntry entry = stack[static_cast<std::size_t>(--sp)];
      // Equality is not pruned: a different primitive at exactly best_d2 may
      // still win the index tie-break. The slack covers the same box-distance
      // rounding as kPruneSlack does for raycast.
      if (entry.t > best_d2 * (1.0 + kPruneSlack)) continue;
      const std::uint64_t payload = entry.payload;
      if (is_leaf(payload)) {
        const std::uint32_t first = leaf_first(payload);
        const std::uint32_t count = leaf_count(payload);
        for (std::uint32_t i = first; i < first + count; ++i) {
          const double d2 = squared_distance(p, tris_[i]);
          if (d2 > best_d2) continue;
          const std::uint32_t prim = order_[i];
          if (d2 < best_d2 || prim < best) {
            best_d2 = d2;
            best = prim;
          }
        }
      } else {
        const Node& n = nodes_[static_cast<std::size_t>(payload)];
        for (int j = 0; j < kFanout; ++j) {
          const std::uint64_t child = n.child[j];
          if (child == kInvalid) break;
          const double d2 = box_distance2(n.lo[j], n.hi[j], p);
          if (d2 > best_d2 * (1.0 + kPruneSlack)) continue;
          stack[static_cast<std::size_t>(sp++)] = StackEntry{d2, child};
        }
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
    tris_.clear();
    root_box_ = Box3d::empty();
    root_payload_ = kInvalid;
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
    if (triangles.empty()) return;
    scratch_.resize(triangles.size());
    bin_index_.resize(triangles.size());
    const auto n = static_cast<std::uint32_t>(triangles.size());
    // The SAH tree is built into a local binary array and widened into `nodes_`
    // once complete; that keeps the query-node type out of the build entirely.
    std::vector<BuildNode> bin;
    bin.reserve(2 * triangles.size());
    if (allow_parallel && n >= kParallelMinPrims) {
      build_parallel(bin, n, root_box, root_centroid_box);
    } else {
      build_node(bin, 0, n, 0, root_box, root_centroid_box);
    }
    root_box_ = root_box;
    // The binary tree above is the exact SAH tree (and the thing the parallel
    // build packs identically). Widen it once into the k-ary traversal layout;
    // this also releases the binary array before the leaf-ordered primitives are
    // materialised.
    widen_nodes(bin);
    // The build only ever permutes `order_`, so the leaf ranges are known now.
    // Materialise the triangles directly in that leaf order: a leaf's primitive
    // tests then walk a contiguous run of the triangle array instead of
    // gathering through `order_`, and `order_[i]` still reports the original
    // index of the triangle that sits at leaf slot i.
    tris_.resize(order_.size());
    for (std::size_t i = 0; i < order_.size(); ++i) tris_[i] = triangles[order_[i]];
  }

  // How many binary levels one traversal node collapses. The constants live
  // before the node type because its array bounds need them.
  static constexpr int kFanoutBits = 2;
  static constexpr int kFanout = 1 << kFanoutBits;
  static constexpr int kBins = 12;
  static constexpr std::uint32_t kMaxLeaf = 4;
  static constexpr int kMaxDepth = 64;
  static constexpr double kTraversalCost = 1.0;
  static constexpr int kParallelDepth = 3;
  static constexpr std::uint32_t kParallelMinPrims = 1u << 15;
  static constexpr std::uint64_t kLeafTag = 1ull << 63;
  static constexpr std::uint64_t kInvalid = ~0ull;
  // A wide node collapses kFanoutBits levels of the binary tree, so its query
  // depth is at most ceil(kMaxDepth / kFanoutBits).
  static constexpr int kStackSize = 1 + (kMaxDepth / kFanoutBits + 1) * (kFanout - 1);

  // A box entry parameter and a primitive distance are computed by different
  // floating-point expressions, so an entry can round a few ulps above a hit
  // that lies exactly on the box's own boundary; coplanar geometry does that for
  // every ray. Pruning therefore compares against the current best widened by
  // this relative slack. Without it the traversal order alone decides whether
  // such a leaf is visited, and the reported closest hit changes by ulps. The
  // slack is far below any distance a caller can distinguish, and it only keeps
  // boxes that were within a few thousand ulps of the best hit anyway.
  static constexpr double kPruneSlack = 1e-12;

  // Binary build node: the SAH tree produced by build_node. It is converted to
  // the wide traversal layout once the build is complete, so this type is never
  // touched by a query.
  struct BuildNode {
    Box3d box;
    std::uint32_t first = 0;   // leaf: first index into order_
    std::uint32_t count = 0;   // leaf: number of primitives; 0 => interior
    std::int32_t left = -1;    // interior child index
    std::int32_t right = -1;   // interior child index
  };

  // Traversal node: kFanout child boxes and child references in one 224-byte
  // record (three and a half 64-byte cache lines).
  // Each child box is the exact binary child's box, stored in double. The entry
  // parameter is compared against the primitive distances computed from the same
  // ray; a narrowed (float) bound can round the entry time above a hit and prune
  // the leaf that holds it, so the bounds must stay exact. `child[j]` is either a
  // leaf payload (kLeafTag | first | count << 32) or the index of another wide
  // node; unused slots are kInvalid and the valid ones are a prefix of the array.
  struct Node {
    double lo[kFanout][3];
    double hi[kFanout][3];
    std::uint64_t child[kFanout];
  };

  // Query stack entry: a node payload plus the entry parameter (raycast) or the
  // box distance lower bound (nearest) that prunes it.
  struct StackEntry {
    double t;
    std::uint64_t payload;
  };

  // Issues a best-effort hardware prefetch for a node that is about to be
  // visited. It has no effect on results.

  static void prefetch(const void* p) noexcept {
#if defined(_MSC_VER) || defined(__SSE__)
    _mm_prefetch(static_cast<const char*>(p), _MM_HINT_T0);
#else
    (void)p;
#endif
  }

  static std::uint64_t leaf_payload(std::uint32_t first, std::uint32_t count) noexcept {
    return kLeafTag | static_cast<std::uint64_t>(first) |
           (static_cast<std::uint64_t>(count) << 32);
  }
  static bool is_leaf(std::uint64_t payload) noexcept { return (payload & kLeafTag) != 0; }
  static std::uint32_t leaf_first(std::uint64_t payload) noexcept {
    return static_cast<std::uint32_t>(payload);
  }
  static std::uint32_t leaf_count(std::uint64_t payload) noexcept {
    // The tag lives in bit 63, which sits inside the high word; mask it off
    // before extracting the count.
    return static_cast<std::uint32_t>((payload >> 32) & 0x7FFFFFFFull);
  }

  // One slot of a wide node while the binary tree is being collapsed.
  struct SlotDesc {
    Box3d box;
    std::int32_t node = -1;  // interior: binary node to expand
    std::uint32_t first = 0;
    std::uint32_t count = 0;
  };

  // Collects the descendants at relative binary depth kFanoutBits that become
  // this wide node's children: a leaf ends the descent early, an interior node
  // exactly at the target depth becomes a child wide node.
  void collect_slots(const std::vector<BuildNode>& bin, std::int32_t bi, int rel,
                     SlotDesc* slots, int& used) const {
    const BuildNode& node = bin[static_cast<std::size_t>(bi)];
    if (node.count > 0 || rel == kFanoutBits) {
      SlotDesc& slot = slots[used++];
      slot.box = node.box;
      slot.first = node.first;
      slot.count = node.count;
      slot.node = node.count > 0 ? -1 : bi;
      return;
    }
    collect_slots(bin, node.left, rel + 1, slots, used);
    collect_slots(bin, node.right, rel + 1, slots, used);
  }

  // Emits the wide node for binary interior node `bi` into `out` and returns its
  // index. Children are emitted depth-first, so the array is a pre-order layout
  // independent of thread count.
  std::uint32_t convert_node(const std::vector<BuildNode>& bin, std::vector<Node>& out,
                             std::int32_t bi) {
    std::array<SlotDesc, kFanout> slots;
    int used = 0;
    const BuildNode& node = bin[static_cast<std::size_t>(bi)];
    collect_slots(bin, node.left, 1, slots.data(), used);
    collect_slots(bin, node.right, 1, slots.data(), used);
    const auto index = static_cast<std::uint32_t>(out.size());
    out.push_back(Node{});
    Node wide{};
    for (int j = 0; j < kFanout; ++j) {
      if (j >= used) {
        wide.child[j] = kInvalid;
        continue;
      }
      const SlotDesc& slot = slots[static_cast<std::size_t>(j)];
      for (int a = 0; a < 3; ++a) {
        wide.lo[j][a] = slot.box.lo[a];
        wide.hi[j][a] = slot.box.hi[a];
      }
      if (slot.node < 0) wide.child[j] = leaf_payload(slot.first, slot.count);
    }
    // Expand the interior children last: `out` may reallocate while they are
    // emitted, so this node is written by index once the recursion returns.
    for (int j = 0; j < used; ++j) {
      const SlotDesc& slot = slots[static_cast<std::size_t>(j)];
      if (slot.node >= 0) wide.child[j] = convert_node(bin, out, slot.node);
    }
    out[index] = wide;
    return index;
  }

  // Replaces the binary node array with the wide traversal layout. A single
  // leaf tree becomes an empty node array with a root leaf payload.
  void widen_nodes(const std::vector<BuildNode>& bin) {
    if (bin.empty()) {
      root_payload_ = kInvalid;
      return;
    }
    const BuildNode& root = bin[0];
    if (root.count > 0) {
      root_payload_ = leaf_payload(root.first, root.count);
      return;
    }
    std::vector<Node> wide;
    wide.reserve(bin.size() / 2 + 1);
    const std::uint32_t root_index = convert_node(bin, wide, 0);
    root_payload_ = root_index;
    nodes_.swap(wide);
  }

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
  std::int32_t build_node(std::vector<BuildNode>& out, std::uint32_t begin, std::uint32_t end,
                          int depth, const Box3d& box, const Box3d& centroid_box) {
    const auto node_index = static_cast<std::int32_t>(out.size());
    out.push_back(BuildNode{});
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
  void build_parallel(std::vector<BuildNode>& out, std::uint32_t n, const Box3d& root_box,
                      const Box3d& root_centroid_box) {
    pending_.clear();
    cut_depth_ = kParallelDepth;
    build_node(out, 0, n, 0, root_box, root_centroid_box);
    cut_depth_ = -1;
    if (pending_.empty()) {
      // The top phase already finished the whole tree.
      return;
    }
    // A single frontier subtree still goes through the append path below: its
    // root node was emitted by the top phase and still needs its children.
    std::vector<std::vector<BuildNode>> subtrees(pending_.size());
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
      out.clear();
      pending_.clear();
      build_node(out, 0, n, 0, root_box, root_centroid_box);
      return;
    }
    for (std::size_t t = 0; t < pending_.size(); ++t) {
      const std::vector<BuildNode>& subtree = subtrees[t];
      const auto base = static_cast<std::int32_t>(out.size());
      // The subtree's root is the node already emitted by the top phase, so its
      // children are packed from `base` on: local index j >= 1 maps to base+j-1.
      const std::int32_t shift = base - 1;
      for (std::size_t j = 1; j < subtree.size(); ++j) {
        BuildNode node = subtree[j];
        if (node.left >= 0) node.left += shift;
        if (node.right >= 0) node.right += shift;
        out.push_back(node);
      }
      const BuildNode& root = subtree.front();
      out[static_cast<std::size_t>(pending_[t].node)].left = root.left + shift;
      out[static_cast<std::size_t>(pending_[t].node)].right = root.right + shift;
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

  static double box_distance2(const double* lo, const double* hi, Point3d p) noexcept {
    double sum = 0.0;
    for (int i = 0; i < 3; ++i) {
      const double l = lo[i];
      const double h = hi[i];
      double d = 0.0;
      if (p[i] < l) d = l - p[i];
      if (p[i] > h) d = p[i] - h;
      sum += d * d;
    }
    return sum;
  }

  // Slab test. A zero direction component makes 1/dir infinite; that case is
  // handled explicitly so a ray starting exactly on a slab plane cannot form
  // 0 * inf = NaN and poison the interval comparisons. On a hit, `t_out` is the
  // entry parameter (0 when the origin is already inside the box); a child box
  // entered behind the current best is only pruned by the caller's limit.
  static bool slab(const Box3d& b, const Ray3d& r, const Vec3d& inv, double tmax,
                   double& t_out) noexcept {
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
    t_out = tmin;
    return true;
  }

  // Same test for the child boxes stored inline in a wide node.
  static bool slab(const double* lo, const double* hi, const Ray3d& r, const Vec3d& inv,
                   double tmax, double& t_out) noexcept {
    double tmin = 0.0;
    double t_limit = tmax;
    for (int i = 0; i < 3; ++i) {
      const double l = lo[i];
      const double h = hi[i];
      if (!std::isfinite(inv[i])) {
        if (r.origin[i] < l || r.origin[i] > h) return false;
        continue;
      }
      double t0 = (l - r.origin[i]) * inv[i];
      double t1 = (h - r.origin[i]) * inv[i];
      if (t0 > t1) std::swap(t0, t1);
      tmin = std::max(tmin, t0);
      t_limit = std::min(t_limit, t1);
      if (tmin > t_limit) return false;
    }
    t_out = tmin;
    return true;
  }

  // Möller-Trumbore. The determinant guard is relative:
  // |det| <= 1e-14 * |e1| * |e2| * |dir| rejects only genuinely near-parallel
  // rays and sliver triangles, independent of the coordinate scale of the
  // triangle. It is a floating-point guard, not an exact predicate; the product
  // is computed in double, so it saturates once |e1| * |e2| * |dir| exceeds
  // DBL_MAX (~1e308). Past that point the guard rejects every ray rather than
  // misclassifying one, and the sliver measure itself has already overflowed.
  static bool intersect_tri(const Ray3d& r, const Triangle3d& t, double dir_len,
                            double& tout) noexcept {
    const Vec3d e1 = t.b - t.a;
    const Vec3d e2 = t.c - t.a;
    const Vec3d pv = cross(r.dir, e2);
    const double det = dot(e1, pv);
    const double scale =
        std::sqrt(dot(e1, e1)) * std::sqrt(dot(e2, e2)) * dir_len;
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
  // Wide traversal nodes (the binary build array is converted and released by
  // `widen_nodes` at the end of each build).
  std::vector<Node> nodes_;
  // Exact union of every primitive box; the query's first slab test.
  Box3d root_box_;
  // Leaf payload when the whole tree is a single leaf, otherwise the index of
  // the root wide node. kInvalid for an empty tree.
  std::uint64_t root_payload_ = kInvalid;
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
