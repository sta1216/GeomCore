#pragma once

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <span>
#include <utility>
#include <vector>

#include <geom/numeric/predicates2.hpp>
#include <geom/primitives/point.hpp>

namespace geom {
inline namespace v1 {

// One input edge of a planar arrangement. `source` identifies the ring or
// operand the edge came from and is copied unchanged to every split piece.
struct ArrangementSegment {
  Point2d a{};
  Point2d b{};
  std::uint32_t source = 0;
};

namespace detail {

// A basic (unsplit) piece produced by split_segments, oriented like its input
// segment.
struct SplitSegment {
  Point2d a;
  Point2d b;
  std::uint32_t source;
};

struct SplitResult {
  std::vector<SplitSegment> segments;
  // Original input-segment pairs that touch, overlap or cross, each once, i < j.
  std::vector<std::pair<std::uint32_t, std::uint32_t>> input_pairs;
  // Maximum perpendicular distance of a computed intersection point to the two
  // supporting lines, in coordinate units.
  double max_residual = 0.0;
};

// A split parameter on one input segment. `p` is the contact coordinate that
// belongs to `t`; endpoints are always emitted from the original coordinates.
struct SplitParam {
  long double t;
  Point2d p;
};

// Input segment plus its original index and cached bounds, used by the sweep.
struct SplitInput {
  ArrangementSegment seg;
  std::uint32_t original;
  double min_x;
  double max_x;
  double min_y;
  double max_y;
};

inline bool split_same_xy(const Point2d& a, const Point2d& b) noexcept {
  return a[0] == b[0] && a[1] == b[1];
}

// Inclusive bounding-box test. For a point already known to be collinear with
// ab (exact sign_orient2d == 0) this is exactly closed-segment containment.
inline bool split_in_box(const Point2d& a, const Point2d& b, const Point2d& p) noexcept {
  return p[0] >= std::min(a[0], b[0]) && p[0] <= std::max(a[0], b[0]) &&
         p[1] >= std::min(a[1], b[1]) && p[1] <= std::max(a[1], b[1]);
}

// Parameter of p on ab, clamped to [0, 1]. p is assumed to lie on the line.
inline long double split_param_on(const Point2d& p, const Point2d& a,
                                  const Point2d& b) noexcept {
  const long double dx = static_cast<long double>(b[0]) - static_cast<long double>(a[0]);
  const long double dy = static_cast<long double>(b[1]) - static_cast<long double>(a[1]);
  const long double len2 = dx * dx + dy * dy;
  if (!(len2 > 0.0L)) return 0.0L;
  const long double px = static_cast<long double>(p[0]) - static_cast<long double>(a[0]);
  const long double py = static_cast<long double>(p[1]) - static_cast<long double>(a[1]);
  const long double t = (px * dx + py * dy) / len2;
  if (!(t > 0.0L)) return 0.0L;
  if (t >= 1.0L) return 1.0L;
  return t;
}

inline void split_add_param(std::vector<SplitParam>& params, long double t,
                            const Point2d& p) {
  if (!(t > 0.0L)) {
    t = 0.0L;
  } else if (t >= 1.0L) {
    t = 1.0L;
  }
  params.push_back(SplitParam{t, p});
}

// |cross(b - a, p - a)| / |b - a|, evaluated in long double.
inline long double split_perp_distance(const Point2d& p, const Point2d& a,
                                       const Point2d& b) noexcept {
  const long double dx = static_cast<long double>(b[0]) - static_cast<long double>(a[0]);
  const long double dy = static_cast<long double>(b[1]) - static_cast<long double>(a[1]);
  const long double len = std::sqrt(dx * dx + dy * dy);
  if (!(len > 0.0L)) return 0.0L;
  const long double px = static_cast<long double>(p[0]) - static_cast<long double>(a[0]);
  const long double py = static_cast<long double>(p[1]) - static_cast<long double>(a[1]);
  return std::fabs(dx * py - dy * px) / len;
}

inline void split_update_residual(long double& max_residual, const Point2d& p,
                                  const Point2d& a, const Point2d& b,
                                  const Point2d& c, const Point2d& d) noexcept {
  const long double to_ab = split_perp_distance(p, a, b);
  const long double to_cd = split_perp_distance(p, c, d);
  const long double worst = std::max(to_ab, to_cd);
  if (worst > max_residual) max_residual = worst;
}

// Classifies one candidate pair with the four exact orientation signs and adds
// the resulting split parameters. Returns true when the pair touches, overlaps
// or crosses.
inline bool split_record_pair(const SplitInput& sa, std::vector<SplitParam>& pa,
                              const SplitInput& sb, std::vector<SplitParam>& pb,
                              long double& max_residual) {
  const Point2d& a = sa.seg.a;
  const Point2d& b = sa.seg.b;
  const Point2d& c = sb.seg.a;
  const Point2d& d = sb.seg.b;

  const int o1 = sign_orient2d(a[0], a[1], b[0], b[1], c[0], c[1]);
  const int o2 = sign_orient2d(a[0], a[1], b[0], b[1], d[0], d[1]);
  const int o3 = sign_orient2d(c[0], c[1], d[0], d[1], a[0], a[1]);
  const int o4 = sign_orient2d(c[0], c[1], d[0], d[1], b[0], b[1]);

  // Collinear: all four signs vanish. The pair is relevant only if the closed
  // bounding boxes overlap; a zero-length overlap is a point touch, a positive
  // one is a collinear overlap.
  if (o1 == 0 && o2 == 0 && o3 == 0 && o4 == 0) {
    if (sa.max_x < sb.min_x || sb.max_x < sa.min_x || sa.max_y < sb.min_y ||
        sb.max_y < sa.min_y) {
      return false;
    }
    if (split_in_box(a, b, c)) {
      const long double t = split_param_on(c, a, b);
      split_add_param(pa, t, c);
      split_update_residual(max_residual, c, a, b, c, d);
    }
    if (split_in_box(a, b, d)) {
      const long double t = split_param_on(d, a, b);
      split_add_param(pa, t, d);
      split_update_residual(max_residual, d, a, b, c, d);
    }
    if (split_in_box(c, d, a)) {
      const long double u = split_param_on(a, c, d);
      split_add_param(pb, u, a);
      split_update_residual(max_residual, a, a, b, c, d);
    }
    if (split_in_box(c, d, b)) {
      const long double u = split_param_on(b, c, d);
      split_add_param(pb, u, b);
      split_update_residual(max_residual, b, a, b, c, d);
    }
    return true;
  }

  // Proper crossing: each segment's endpoints lie strictly on opposite sides
  // of the other segment.
  if (((o1 > 0 && o2 < 0) || (o1 < 0 && o2 > 0)) &&
      ((o3 > 0 && o4 < 0) || (o3 < 0 && o4 > 0))) {
    const long double ax = static_cast<long double>(a[0]);
    const long double ay = static_cast<long double>(a[1]);
    const long double rx = static_cast<long double>(b[0]) - ax;
    const long double ry = static_cast<long double>(b[1]) - ay;
    const long double cx = static_cast<long double>(c[0]);
    const long double cy = static_cast<long double>(c[1]);
    const long double sx = static_cast<long double>(d[0]) - cx;
    const long double sy = static_cast<long double>(d[1]) - cy;
    const long double denom = rx * sy - ry * sx;
    if (denom == 0.0L) return false;
    long double t = ((cx - ax) * sy - (cy - ay) * sx) / denom;
    long double u = ((cx - ax) * ry - (cy - ay) * rx) / denom;
    if (!(t > 0.0L)) {
      t = 0.0L;
    } else if (t >= 1.0L) {
      t = 1.0L;
    }
    if (!(u > 0.0L)) {
      u = 0.0L;
    } else if (u >= 1.0L) {
      u = 1.0L;
    }

    // Reuse the exact coordinate of an existing endpoint whenever the
    // intersection degenerates onto one; otherwise take a + t*(b - a).
    Point2d contact;
    if (t <= 0.0L) {
      contact = a;
    } else if (t >= 1.0L) {
      contact = b;
    } else if (u <= 0.0L) {
      contact = c;
    } else if (u >= 1.0L) {
      contact = d;
    } else {
      contact = Point2d{static_cast<double>(ax + t * rx), static_cast<double>(ay + t * ry)};
    }
    split_add_param(pa, t, contact);
    split_add_param(pb, u, contact);
    split_update_residual(max_residual, contact, a, b, c, d);
    return true;
  }

  // Endpoint touch: an endpoint of one segment lies on the other. The signs are
  // exact, so the only remaining work is the closed bounding-box check.
  bool touched = false;
  if (o1 == 0 && split_in_box(a, b, c)) {
    touched = true;
    split_add_param(pa, split_param_on(c, a, b), c);
    split_update_residual(max_residual, c, a, b, c, d);
  }
  if (o2 == 0 && split_in_box(a, b, d)) {
    touched = true;
    split_add_param(pa, split_param_on(d, a, b), d);
    split_update_residual(max_residual, d, a, b, c, d);
  }
  if (o3 == 0 && split_in_box(c, d, a)) {
    touched = true;
    split_add_param(pb, split_param_on(a, c, d), a);
    split_update_residual(max_residual, a, a, b, c, d);
  }
  if (o4 == 0 && split_in_box(c, d, b)) {
    touched = true;
    split_add_param(pb, split_param_on(b, c, d), b);
    split_update_residual(max_residual, b, a, b, c, d);
  }
  return touched;
}

// Splits every non-degenerate input segment at all intersection points with the
// other segments.
//
// Contract:
//   - zero-length inputs (a == b) are dropped;
//   - the output keeps the input order and each piece keeps its input's
//     direction and `source`;
//   - endpoint touches, collinear overlaps and proper crossings all split both
//     segments and are recorded once in `input_pairs` as (i, j), i < j, sorted;
//   - `max_residual` is the largest perpendicular distance of a computed
//     intersection point to the two supporting lines.
//
// Candidate pairs come from a min-x sweep with an active set; the exact
// classification uses sign_orient2d, so disjoint and exact endpoint cases never
// depend on a tolerance. Intersection parameters are evaluated in long double
// and clamped to [0, 1]. The function throws no exceptions of its own
// (standard allocation failures may still propagate).
inline detail::SplitResult split_segments(std::span<const ArrangementSegment> segments) {
  detail::SplitResult result;

  std::vector<detail::SplitInput> inputs;
  inputs.reserve(segments.size());
  for (std::size_t i = 0; i < segments.size(); ++i) {
    const ArrangementSegment& s = segments[i];
    if (detail::split_same_xy(s.a, s.b)) continue;
    detail::SplitInput in;
    in.seg = s;
    in.original = static_cast<std::uint32_t>(i);
    in.min_x = std::min(s.a[0], s.b[0]);
    in.max_x = std::max(s.a[0], s.b[0]);
    in.min_y = std::min(s.a[1], s.b[1]);
    in.max_y = std::max(s.a[1], s.b[1]);
    inputs.push_back(in);
  }

  std::vector<std::vector<detail::SplitParam>> params(inputs.size());
  for (std::size_t i = 0; i < inputs.size(); ++i) {
    params[i].push_back(detail::SplitParam{0.0L, inputs[i].seg.a});
    params[i].push_back(detail::SplitParam{1.0L, inputs[i].seg.b});
  }

  // Broad phase: sweep by min x, keep an active set, prune by max x and by the
  // y interval. Every surviving pair has overlapping closed bounding boxes.
  std::vector<std::size_t> order(inputs.size());
  for (std::size_t i = 0; i < order.size(); ++i) order[i] = i;
  std::sort(order.begin(), order.end(), [&inputs](std::size_t lhs, std::size_t rhs) {
    if (inputs[lhs].min_x != inputs[rhs].min_x) return inputs[lhs].min_x < inputs[rhs].min_x;
    return inputs[lhs].original < inputs[rhs].original;
  });

  std::vector<std::size_t> active;
  std::vector<std::pair<std::uint32_t, std::uint32_t>> pairs;
  long double max_residual = 0.0L;

  for (const std::size_t current : order) {
    const detail::SplitInput& cur = inputs[current];
    active.erase(std::remove_if(active.begin(), active.end(),
                                [&cur, &inputs](std::size_t other) {
                                  return inputs[other].max_x < cur.min_x;
                                }),
                 active.end());
    for (const std::size_t other : active) {
      const detail::SplitInput& prev = inputs[other];
      if (prev.max_y < cur.min_y || cur.max_y < prev.min_y) continue;
      if (detail::split_record_pair(prev, params[other], cur, params[current],
                                    max_residual)) {
        const std::uint32_t lo = std::min(prev.original, cur.original);
        const std::uint32_t hi = std::max(prev.original, cur.original);
        pairs.push_back(std::make_pair(lo, hi));
      }
    }
    active.push_back(current);
  }

  std::sort(pairs.begin(), pairs.end());
  pairs.erase(std::unique(pairs.begin(), pairs.end()), pairs.end());
  result.input_pairs = std::move(pairs);

  for (std::size_t i = 0; i < inputs.size(); ++i) {
    std::vector<detail::SplitParam>& ps = params[i];
    std::sort(ps.begin(), ps.end(), [](const detail::SplitParam& lhs,
                                       const detail::SplitParam& rhs) {
      if (lhs.t != rhs.t) return lhs.t < rhs.t;
      if (lhs.p[0] != rhs.p[0]) return lhs.p[0] < rhs.p[0];
      return lhs.p[1] < rhs.p[1];
    });
    ps.erase(std::unique(ps.begin(), ps.end(),
                         [](const detail::SplitParam& lhs, const detail::SplitParam& rhs) {
                           return lhs.t == rhs.t;
                         }),
             ps.end());

    const ArrangementSegment& source = inputs[i].seg;
    for (std::size_t k = 0; k + 1 < ps.size(); ++k) {
      const Point2d p0 = ps[k].t <= 0.0L
                             ? source.a
                             : (ps[k].t >= 1.0L ? source.b : ps[k].p);
      const Point2d p1 = ps[k + 1].t <= 0.0L
                             ? source.a
                             : (ps[k + 1].t >= 1.0L ? source.b : ps[k + 1].p);
      if (detail::split_same_xy(p0, p1)) continue;
      result.segments.push_back(detail::SplitSegment{p0, p1, source.source});
    }
  }

  result.max_residual = static_cast<double>(max_residual);
  return result;
}

}  // namespace detail

}  // namespace v1
}  // namespace geom
