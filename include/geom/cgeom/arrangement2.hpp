#pragma once

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>
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
// The direction is normalized by its largest component so that the dot products
// neither overflow nor underflow at extreme (but finite) coordinate scales.
inline long double split_param_on(const Point2d& p, const Point2d& a,
                                  const Point2d& b) noexcept {
  long double dx = static_cast<long double>(b[0]) - static_cast<long double>(a[0]);
  long double dy = static_cast<long double>(b[1]) - static_cast<long double>(a[1]);
  const long double dm = std::max(std::fabs(dx), std::fabs(dy));
  if (!(dm > 0.0L)) return 0.0L;
  dx /= dm;
  dy /= dm;
  const long double px =
      (static_cast<long double>(p[0]) - static_cast<long double>(a[0])) / dm;
  const long double py =
      (static_cast<long double>(p[1]) - static_cast<long double>(a[1])) / dm;
  const long double t = (px * dx + py * dy) / (dx * dx + dy * dy);
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

// Exact orientation determinant det = cross(b - a, c - a), evaluated as an
// Expansion. The caller passes a power-of-two-scaled frame so the products stay
// in range; the sign of this determinant is exact in that frame.
inline Expansion split_orient_det(double ax, double ay, double bx, double by,
                                  double cx, double cy) {
  const Expansion acx = Expansion::difference(ax, cx);
  const Expansion acy = Expansion::difference(ay, cy);
  const Expansion bcx = Expansion::difference(bx, cx);
  const Expansion bcy = Expansion::difference(by, cy);

  Expansion det = Expansion::mul(acx, bcy);
  det.add_scaled(Expansion::mul(acy, bcx), -1.0);
  return det;
}

// Sum an Expansion in long double, largest component first. On MSVC long double
// and double are the same type, but the ordering still avoids needless loss.
inline long double split_exp_to_long_double(const Expansion& e) noexcept {
  long double v = 0.0L;
  for (std::size_t i = e.size(); i > 0; --i) {
    v += static_cast<long double>(e[i - 1]);
  }
  return v;
}

// Clamp an intersection parameter. A non-finite value must not drop the pair
// (the pair was already classified exactly); +-inf clamps to the nearer end.
// NaN carries no direction information, so it is split in the middle.
inline long double split_clamp_parameter(long double t) noexcept {
  if (std::isnan(t)) return 0.5L;
  if (t <= 0.0L) return 0.0L;
  if (t >= 1.0L) return 1.0L;
  return t;
}

// |cross(b - a, p - a)| / |b - a|, evaluated in long double with both the
// direction and the offset normalized. This keeps 1e200-scale inputs finite and
// 1e-200-scale inputs nonzero. A non-finite result is saturated to infinity.
inline long double split_perp_distance(const Point2d& p, const Point2d& a,
                                       const Point2d& b) noexcept {
  const long double dx = static_cast<long double>(b[0]) - static_cast<long double>(a[0]);
  const long double dy = static_cast<long double>(b[1]) - static_cast<long double>(a[1]);
  const long double dm = std::max(std::fabs(dx), std::fabs(dy));
  if (!(dm > 0.0L) || !std::isfinite(dm)) {
    return std::numeric_limits<long double>::infinity();
  }
  const long double ux = dx / dm;
  const long double uy = dy / dm;

  const long double px = static_cast<long double>(p[0]) - static_cast<long double>(a[0]);
  const long double py = static_cast<long double>(p[1]) - static_cast<long double>(a[1]);
  const long double pm = std::max(std::fabs(px), std::fabs(py));
  if (pm == 0.0L) return 0.0L;
  if (!std::isfinite(pm)) {
    return std::numeric_limits<long double>::infinity();
  }
  const long double qx = px / pm;
  const long double qy = py / pm;

  const long double cross = ux * qy - uy * qx;
  const long double len = std::hypot(ux, uy);
  long double distance = std::fabs(cross) / len * pm;
  if (!std::isfinite(distance)) {
    return std::numeric_limits<long double>::infinity();
  }
  return distance;
}

inline void split_update_residual(long double& max_residual, const Point2d& p,
                                  const Point2d& a, const Point2d& b,
                                  const Point2d& c, const Point2d& d) noexcept {
  const long double to_ab = split_perp_distance(p, a, b);
  const long double to_cd = split_perp_distance(p, c, d);
  if (!std::isfinite(to_ab) || !std::isfinite(to_cd)) {
    max_residual = std::numeric_limits<long double>::infinity();
    return;
  }
  const long double worst = to_ab > to_cd ? to_ab : to_cd;
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
  // of the other segment. The classification above is exact, so this branch
  // must never be vetoed by a floating-point denominator test. Instead the
  // parameters are ratios of exact orientation determinants evaluated in a
  // per-pair power-of-two frame: scaling by 2^k is exact, preserves signs, and
  // keeps both the near-parallel (2^60 products cancelling to 2) and the
  // 1e+-200 uniform-scale cases representable.
  if (((o1 > 0 && o2 < 0) || (o1 < 0 && o2 > 0)) &&
      ((o3 > 0 && o4 < 0) || (o3 < 0 && o4 > 0))) {
    const double max_coord = std::max({std::fabs(a[0]), std::fabs(a[1]), std::fabs(b[0]),
                                       std::fabs(b[1]), std::fabs(c[0]), std::fabs(c[1]),
                                       std::fabs(d[0]), std::fabs(d[1])});
    const int shift =
        (max_coord > 0.0 && std::isfinite(max_coord)) ? -std::ilogb(max_coord) : 0;
    const double sax = std::scalbn(a[0], shift);
    const double say = std::scalbn(a[1], shift);
    const double sbx = std::scalbn(b[0], shift);
    const double sby = std::scalbn(b[1], shift);
    const double scx = std::scalbn(c[0], shift);
    const double scy = std::scalbn(c[1], shift);
    const double sdx = std::scalbn(d[0], shift);
    const double sdy = std::scalbn(d[1], shift);

    const Expansion det1 = split_orient_det(sax, say, sbx, sby, scx, scy);
    const Expansion det2 = split_orient_det(sax, say, sbx, sby, sdx, sdy);
    const Expansion det3 = split_orient_det(scx, scy, sdx, sdy, sax, say);
    const Expansion det4 = split_orient_det(scx, scy, sdx, sdy, sbx, sby);

    // t = O3 / (O3 - O4), u = O1 / (O1 - O2). For a proper crossing the two
    // terms in each denominator have opposite signs, so the denominator is a
    // sum of magnitudes and the exact difference is evaluated as an Expansion
    // before the conversion.
    Expansion denom_t = det3;
    denom_t.add_scaled(det4, -1.0);
    Expansion denom_u = det1;
    denom_u.add_scaled(det2, -1.0);

    const long double o1_ld = split_exp_to_long_double(det1);
    const long double o3_ld = split_exp_to_long_double(det3);
    const long double denom_t_ld = split_exp_to_long_double(denom_t);
    const long double denom_u_ld = split_exp_to_long_double(denom_u);

    long double t = 0.5L;
    long double u = 0.5L;
    if (denom_t_ld != 0.0L && o3_ld != 0.0L) t = o3_ld / denom_t_ld;
    if (denom_u_ld != 0.0L && o1_ld != 0.0L) u = o1_ld / denom_u_ld;
    t = split_clamp_parameter(t);
    u = split_clamp_parameter(u);

    // Reuse the exact coordinate of an existing endpoint whenever the
    // intersection degenerates onto one; otherwise evaluate the convex
    // combination (1 - t) * a + t * b, which cannot overflow when the endpoint
    // coordinates are finite even if b - a would.
    const long double ax = static_cast<long double>(a[0]);
    const long double ay = static_cast<long double>(a[1]);
    const long double bx = static_cast<long double>(b[0]);
    const long double by = static_cast<long double>(b[1]);
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
      const long double omt = 1.0L - t;
      contact = Point2d{static_cast<double>(omt * ax + t * bx),
                        static_cast<double>(omt * ay + t * by)};
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
// depend on a tolerance. Proper-crossing parameters are exact determinant
// ratios evaluated in a power-of-two-scaled frame, then clamped to [0, 1].
// The function throws no exceptions of its own (standard allocation failures
// may still propagate).
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
