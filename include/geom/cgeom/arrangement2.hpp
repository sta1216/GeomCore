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

// Largest supported exponent span between finite input coordinates. The exact
// determinant ratios are evaluated in the unscaled frame first and in one
// power-of-two scaled frame second; within kMaxSpan those two frames cover the
// representable ratios. Inputs whose finite coordinates span more than this are
// outside the documented contract.
inline constexpr int kMaxSpan = 1314;

inline bool split_same_xy(const Point2d& a, const Point2d& b) noexcept {
  return a[0] == b[0] && a[1] == b[1];
}

// Inclusive bounding-box test. For a point already known to be collinear with
// ab (exact sign_orient2d == 0) this is exactly closed-segment containment.
inline bool split_in_box(const Point2d& a, const Point2d& b, const Point2d& p) noexcept {
  return p[0] >= std::min(a[0], b[0]) && p[0] <= std::max(a[0], b[0]) &&
         p[1] >= std::min(a[1], b[1]) && p[1] <= std::max(a[1], b[1]);
}

// Power-of-two shift that brings the largest coordinate of a, b, p into
// [0.5, 1). Used when a coordinate difference overflows double (MSVC long
// double is double): dividing inf by inf would otherwise produce NaN.
inline int split_scale_shift(const Point2d& a, const Point2d& b,
                             const Point2d& p) noexcept {
  double m = 0.0;
  for (const double c : {a[0], a[1], b[0], b[1], p[0], p[1]}) {
    const double mag = std::fabs(c);
    if (mag > m) m = mag;
  }
  if (!(m > 0.0) || !std::isfinite(m)) return 0;
  return -std::ilogb(m);
}

// Parameter of p on ab, clamped to [0, 1]. p is assumed to lie on the line.
// The direction is normalized by its largest component so that the dot products
// neither overflow nor underflow at extreme (but finite) coordinate scales.
// When b - a overflows (for example -1e308 to 1e308), every involved coordinate
// is first rescaled by a power of two; t is invariant under that scaling.
inline long double split_param_on(const Point2d& p, const Point2d& a,
                                  const Point2d& b) noexcept {
  double ax = a[0];
  double ay = a[1];
  double bx = b[0];
  double by = b[1];
  double px_in = p[0];
  double py_in = p[1];
  long double dx = static_cast<long double>(bx) - static_cast<long double>(ax);
  long double dy = static_cast<long double>(by) - static_cast<long double>(ay);
  long double dm = std::max(std::fabs(dx), std::fabs(dy));
  if (!std::isfinite(dm)) {
    const int shift = split_scale_shift(a, b, p);
    ax = std::scalbn(ax, shift);
    ay = std::scalbn(ay, shift);
    bx = std::scalbn(bx, shift);
    by = std::scalbn(by, shift);
    px_in = std::scalbn(px_in, shift);
    py_in = std::scalbn(py_in, shift);
    dx = static_cast<long double>(bx) - static_cast<long double>(ax);
    dy = static_cast<long double>(by) - static_cast<long double>(ay);
    dm = std::max(std::fabs(dx), std::fabs(dy));
  }
  if (!(dm > 0.0L) || !std::isfinite(dm)) return 0.0L;
  dx /= dm;
  dy /= dm;
  const long double px = (static_cast<long double>(px_in) - static_cast<long double>(ax)) / dm;
  const long double py = (static_cast<long double>(py_in) - static_cast<long double>(ay)) / dm;
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

// Orders two distinct contacts that share a (rounded) parameter. The dominant
// coordinate of the segment is monotone from a to b, so comparing it with the
// segment direction recovers the spatial order that the parameter lost: on a
// 1e308-long segment two contacts one unit apart can round to the same t.
inline bool split_contact_before(const Point2d& a, const Point2d& b,
                                 const Point2d& lhs, const Point2d& rhs) noexcept {
  const double dx = b[0] - a[0];
  const double dy = b[1] - a[1];
  const bool use_x = std::fabs(dx) >= std::fabs(dy);
  const double direction = use_x ? dx : dy;
  const double l = use_x ? lhs[0] : lhs[1];
  const double r = use_x ? rhs[0] : rhs[1];
  return direction >= 0.0 ? l < r : l > r;
}

// Exact orientation determinant det = cross(b - a, c - a), evaluated as an
// Expansion. The same routine is used both on the original coordinates and on
// a power-of-two-scaled frame; scaling preserves the exact sign.
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

// True when every component of e is finite. An unscaled product can overflow to
// inf (or produce NaN through inf - inf), in which case that frame must be
// discarded and the scaled retry used instead.
inline bool split_exp_finite(const Expansion& e) noexcept {
  for (std::size_t i = 0; i < e.size(); ++i) {
    if (!std::isfinite(static_cast<long double>(e[i]))) return false;
  }
  return true;
}

// Computes t = O3/(O3 - O4) and u = O1/(O1 - O2) from the four exact
// orientation determinants. Returns false when the frame overflowed or a
// determinant that must be nonzero for a proper crossing underflowed to zero;
// the caller then retries in the scaled frame.
inline bool split_exact_params(const Expansion& det1, const Expansion& det2,
                               const Expansion& det3, const Expansion& det4,
                               long double& t, long double& u) {
  if (!split_exp_finite(det1) || !split_exp_finite(det2) || !split_exp_finite(det3) ||
      !split_exp_finite(det4)) {
    return false;
  }
  Expansion denom_t = det3;
  denom_t.add_scaled(det4, -1.0);
  Expansion denom_u = det1;
  denom_u.add_scaled(det2, -1.0);
  if (!split_exp_finite(denom_t) || !split_exp_finite(denom_u)) return false;

  const long double o1 = split_exp_to_long_double(det1);
  const long double o2 = split_exp_to_long_double(det2);
  const long double o3 = split_exp_to_long_double(det3);
  const long double o4 = split_exp_to_long_double(det4);
  const long double denom_t_ld = split_exp_to_long_double(denom_t);
  const long double denom_u_ld = split_exp_to_long_double(denom_u);
  if (!std::isfinite(o1) || !std::isfinite(o2) || !std::isfinite(o3) ||
      !std::isfinite(o4) || !std::isfinite(denom_t_ld) || !std::isfinite(denom_u_ld)) {
    return false;
  }
  // For a proper crossing all four signs are strict, so a zero here means the
  // frame lost information (underflow), not a degenerate configuration.
  if (o1 == 0.0L || o2 == 0.0L || o3 == 0.0L || o4 == 0.0L ||
      denom_t_ld == 0.0L || denom_u_ld == 0.0L) {
    return false;
  }
  const long double tt = o3 / denom_t_ld;
  const long double uu = o1 / denom_u_ld;
  if (!std::isfinite(tt) || !std::isfinite(uu)) return false;
  // Materialize the parameters as double so the decision does not depend on the
  // width of long double (80-bit on Linux, 64-bit on MSVC/macOS).
  const double t_double = static_cast<double>(tt);
  const double u_double = static_cast<double>(uu);
  if (!std::isfinite(t_double) || !std::isfinite(u_double)) return false;
  // The signs are strict for a proper crossing, so both ratios are strictly
  // inside (0, 1). A ratio that rounds to an endpoint (or underflows to zero in
  // double) is not representable: emitting it would add the contact at t = 0 or
  // t = 1, where the t-only dedupe drops it and leaves the other segment's
  // vertex on this segment's interior without a split. Report the frame as
  // unusable instead; the caller then retries scaled, and finally leaves both
  // inputs unsplit with max_residual = +infinity rather than silently losing
  // the vertex.
  if (t_double <= 0.0 || t_double >= 1.0 || u_double <= 0.0 || u_double >= 1.0) return false;
  t = t_double;
  u = u_double;
  return true;
}

// Clamp a finite intersection parameter. Non-finite values are rejected by
// split_exact_params before reaching this point, so NaN here is a defensive
// fallback to the start endpoint, never a silent 0.5 midpoint.
inline long double split_clamp_parameter(long double t) noexcept {
  if (std::isnan(t)) return 0.0L;
  if (t <= 0.0L) return 0.0L;
  if (t >= 1.0L) return 1.0L;
  return t;
}

// Point at parameter t on ab. Endpoints are returned exactly; the interior is
// evaluated as a convex combination, which stays finite for finite endpoints
// even when b - a would overflow.
inline Point2d split_point_at(const Point2d& a, const Point2d& b, long double t) noexcept {
  if (t <= 0.0L) return a;
  if (t >= 1.0L) return b;
  const long double omt = 1.0L - t;
  return Point2d{static_cast<double>(omt * static_cast<long double>(a[0]) +
                                     t * static_cast<long double>(b[0])),
                 static_cast<double>(omt * static_cast<long double>(a[1]) +
                                     t * static_cast<long double>(b[1]))};
}

// |cross(b - a, p - a)| / |b - a|, evaluated in long double with both the
// direction and the offset normalized. This keeps 1e200-scale inputs finite and
// 1e-200-scale inputs nonzero. A non-finite result is saturated to infinity.
// If b - a overflows double, the frame is rescaled by a power of two and the
// distance is converted back to the original units (divided by the scale).
inline long double split_perp_distance(const Point2d& p, const Point2d& a,
                                       const Point2d& b) noexcept {
  double ax = a[0];
  double ay = a[1];
  double bx = b[0];
  double by = b[1];
  double px_in = p[0];
  double py_in = p[1];
  long double dx = static_cast<long double>(bx) - static_cast<long double>(ax);
  long double dy = static_cast<long double>(by) - static_cast<long double>(ay);
  long double dm = std::max(std::fabs(dx), std::fabs(dy));
  int shift = 0;
  if (!std::isfinite(dm)) {
    shift = split_scale_shift(a, b, p);
    ax = std::scalbn(ax, shift);
    ay = std::scalbn(ay, shift);
    bx = std::scalbn(bx, shift);
    by = std::scalbn(by, shift);
    px_in = std::scalbn(px_in, shift);
    py_in = std::scalbn(py_in, shift);
    dx = static_cast<long double>(bx) - static_cast<long double>(ax);
    dy = static_cast<long double>(by) - static_cast<long double>(ay);
    dm = std::max(std::fabs(dx), std::fabs(dy));
  }
  if (!(dm > 0.0L) || !std::isfinite(dm)) {
    return std::numeric_limits<long double>::infinity();
  }
  const long double ux = dx / dm;
  const long double uy = dy / dm;

  const long double px = static_cast<long double>(px_in) - static_cast<long double>(ax);
  const long double py = static_cast<long double>(py_in) - static_cast<long double>(ay);
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
  if (shift != 0) distance = std::ldexp(distance, -shift);
  if (!std::isfinite(distance)) {
    return std::numeric_limits<long double>::infinity();
  }
  return distance;
}

// Worst perpendicular distance from p to the two supporting lines, in the same
// normalized units as split_perp_distance. A non-finite distance is returned as
// infinity so the caller can mark the split point as not representable.
inline long double split_contact_error(const Point2d& p, const Point2d& a,
                                       const Point2d& b, const Point2d& c,
                                       const Point2d& d) noexcept {
  const long double to_ab = split_perp_distance(p, a, b);
  const long double to_cd = split_perp_distance(p, c, d);
  if (!std::isfinite(to_ab) || !std::isfinite(to_cd)) {
    return std::numeric_limits<long double>::infinity();
  }
  return to_ab > to_cd ? to_ab : to_cd;
}

inline void split_update_residual(long double& max_residual, const Point2d& p,
                                  const Point2d& a, const Point2d& b,
                                  const Point2d& c, const Point2d& d) noexcept {
  const long double worst = split_contact_error(p, a, b, c, d);
  if (!std::isfinite(worst)) {
    max_residual = std::numeric_limits<long double>::infinity();
    return;
  }
  if (worst > max_residual) max_residual = worst;
}

// Cheap, non-exact rejection for the broad-phase candidate stream. When both
// endpoints of one segment lie clearly on the same side of the other segment's
// line, the two segments cannot touch. The threshold is deliberately far above
// the rounding error of the double determinants (about 1e-16 relative), so a
// pair is only rejected when the double evaluation is unambiguous; exact
// collinearity, endpoint touches, overflow to infinity and underflow to zero all
// fall through to the exact path. This removes the near-parallel disjoint pairs
// that dominate a dense polygon-ring candidate set without weakening exactness.
inline bool split_clearly_disjoint(const Point2d& a, const Point2d& b,
                                   const Point2d& c, const Point2d& d) noexcept {
  const double abx = b[0] - a[0];
  const double aby = b[1] - a[1];
  const double acx = c[0] - a[0];
  const double acy = c[1] - a[1];
  const double adx = d[0] - a[0];
  const double ady = d[1] - a[1];
  const double d1 = abx * acy - aby * acx;
  const double d2 = abx * ady - aby * adx;
  const double t1 = 1e-12 * (std::fabs(abx * acy) + std::fabs(aby * acx));
  const double t2 = 1e-12 * (std::fabs(abx * ady) + std::fabs(aby * adx));
  if ((d1 > t1 && d2 > t2) || (d1 < -t1 && d2 < -t2)) return true;

  const double cdx = d[0] - c[0];
  const double cdy = d[1] - c[1];
  const double cax = a[0] - c[0];
  const double cay = a[1] - c[1];
  const double cbx = b[0] - c[0];
  const double cby = b[1] - c[1];
  const double d3 = cdx * cay - cdy * cax;
  const double d4 = cdx * cby - cdy * cbx;
  const double t3 = 1e-12 * (std::fabs(cdx * cay) + std::fabs(cdy * cax));
  const double t4 = 1e-12 * (std::fabs(cdx * cby) + std::fabs(cdy * cbx));
  return (d3 > t3 && d4 > t4) || (d3 < -t3 && d4 < -t4);
}

// Double-double arithmetic for the proper-crossing parameters. Four exact
// determinant ratios dominate the cost of a dense arrangement (the Expansion
// path allocates for every product), while the parameters only need enough
// precision for the contact point and the residual check. Two doubles carry
// about 106 bits, so a plain double-double evaluation of the orientation
// determinants is far more accurate than the long double division that follows;
// the exact signs still come from sign_orient2d. Overflow, underflow and any
// cancellation fall back to the exact Expansion path.
struct SplitDd {
  double hi = 0.0;
  double lo = 0.0;
};

inline SplitDd split_dd_add(SplitDd a, SplitDd b) noexcept {
  const double s = a.hi + b.hi;
  const double v = s - a.hi;
  double e = (a.hi - (s - v)) + (b.hi - v);
  e += a.lo + b.lo;
  const double t = s + e;
  e -= (t - s);
  return SplitDd{t, e};
}

inline SplitDd split_dd_sub(SplitDd a, SplitDd b) noexcept {
  return split_dd_add(a, SplitDd{-b.hi, -b.lo});
}

inline SplitDd split_dd_mul(SplitDd a, SplitDd b) noexcept {
  const double p = a.hi * b.hi;
  double e = std::fma(a.hi, b.hi, -p);
  e += a.hi * b.lo + a.lo * b.hi;
  const double t = p + e;
  e -= (t - p);
  return SplitDd{t, e};
}

inline SplitDd split_dd_diff(double x, double y) noexcept {
  double s = 0.0;
  double e = 0.0;
  two_diff(x, y, s, e);
  return SplitDd{s, e};
}

inline SplitDd split_det_dd(double ax, double ay, double bx, double by,
                            double cx, double cy) noexcept {
  const SplitDd abx = split_dd_diff(bx, ax);
  const SplitDd aby = split_dd_diff(by, ay);
  const SplitDd acx = split_dd_diff(cx, ax);
  const SplitDd acy = split_dd_diff(cy, ay);
  return split_dd_sub(split_dd_mul(abx, acy), split_dd_mul(aby, acx));
}

inline bool split_dd_finite(SplitDd v) noexcept {
  return std::isfinite(v.hi) && std::isfinite(v.lo);
}

inline bool split_dd_params(double ax, double ay, double bx, double by, double cx,
                            double cy, double dx, double dy, long double& t,
                            long double& u) noexcept {
  const SplitDd det1 = split_det_dd(ax, ay, bx, by, cx, cy);
  const SplitDd det2 = split_det_dd(ax, ay, bx, by, dx, dy);
  const SplitDd det3 = split_det_dd(cx, cy, dx, dy, ax, ay);
  const SplitDd det4 = split_det_dd(cx, cy, dx, dy, bx, by);
  if (!split_dd_finite(det1) || !split_dd_finite(det2) || !split_dd_finite(det3) ||
      !split_dd_finite(det4)) {
    return false;
  }
  const long double o1 = static_cast<long double>(det1.hi) + det1.lo;
  const long double o2 = static_cast<long double>(det2.hi) + det2.lo;
  const long double o3 = static_cast<long double>(det3.hi) + det3.lo;
  const long double o4 = static_cast<long double>(det4.hi) + det4.lo;
  // Subnormal determinants lose enough precision that the ratio is no longer
  // trustworthy; the exact path handles those rare extreme-span pairs.
  constexpr long double kMinDdMagnitude = 1e-290L;
  if (std::fabs(o1) < kMinDdMagnitude || std::fabs(o2) < kMinDdMagnitude ||
      std::fabs(o3) < kMinDdMagnitude || std::fabs(o4) < kMinDdMagnitude) {
    return false;
  }
  // For a proper crossing the two endpoints of each segment lie on opposite
  // sides, so the denominators are sums of two same-sign magnitudes: the
  // subtraction below cannot cancel.
  const long double denom_t = o3 - o4;
  const long double denom_u = o1 - o2;
  if (!(denom_t > 0.0L) && !(denom_t < 0.0L)) return false;
  if (!(denom_u > 0.0L) && !(denom_u < 0.0L)) return false;
  if (o1 == 0.0L || o2 == 0.0L || o3 == 0.0L || o4 == 0.0L) return false;

  const long double tt = o3 / denom_t;
  const long double uu = o1 / denom_u;
  if (!std::isfinite(tt) || !std::isfinite(uu)) return false;
  // Match the exact path: materialize as double and reject ratios that round to
  // an endpoint, which would silently lose the split point.
  const double t_double = static_cast<double>(tt);
  const double u_double = static_cast<double>(uu);
  if (!std::isfinite(t_double) || !std::isfinite(u_double)) return false;
  if (t_double <= 0.0 || t_double >= 1.0 || u_double <= 0.0 || u_double >= 1.0) {
    return false;
  }
  t = t_double;
  u = u_double;
  return true;
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

  if (split_clearly_disjoint(a, b, c, d)) return false;

  // Fast path: two non-collinear segments that share an exact endpoint meet only
  // at that point, which is already a vertex of both segments. Adjacent ring
  // edges share an endpoint by construction and dominate the candidate count, so
  // recording the pair here avoids four exact orientation tests per pair. A
  // collinear pair still needs the general path below because the overlap can
  // extend past the shared endpoint.
  const Point2d* shared = nullptr;
  const Point2d* far_a = nullptr;
  const Point2d* far_b = nullptr;
  if (a[0] == c[0] && a[1] == c[1]) {
    shared = &a;
    far_a = &b;
    far_b = &d;
  } else if (a[0] == d[0] && a[1] == d[1]) {
    shared = &a;
    far_a = &b;
    far_b = &c;
  } else if (b[0] == c[0] && b[1] == c[1]) {
    shared = &b;
    far_a = &a;
    far_b = &d;
  } else if (b[0] == d[0] && b[1] == d[1]) {
    shared = &b;
    far_a = &a;
    far_b = &c;
  }
  if (shared != nullptr &&
      sign_orient2d((*shared)[0], (*shared)[1], (*far_a)[0], (*far_a)[1],
                    (*far_b)[0], (*far_b)[1]) != 0) {
    return true;
  }

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
  // must never be vetoed by a floating-point denominator test. Parameters are
  // ratios of exact orientation determinants, evaluated in the unscaled frame
  // first; a mixed-span pair (for example 2e200 and 1e-160) is exact there even
  // though a single globally scaled frame would collapse the tiny coordinate to
  // zero. Only when the unscaled products overflow/underflow is the 2^k-scaled
  // frame retried. If neither frame yields a ratio strictly inside (0, 1), the
  // split point is not representable: the pair is still recorded and both input
  // segments are left unsplit, while max_residual becomes +infinity so the
  // failure is visible instead of silently bending the geometry.
  if (((o1 > 0 && o2 < 0) || (o1 < 0 && o2 > 0)) &&
      ((o3 > 0 && o4 < 0) || (o3 < 0 && o4 > 0))) {
    long double t = 0.0L;
    long double u = 0.0L;
    bool have_params = split_dd_params(a[0], a[1], b[0], b[1], c[0], c[1], d[0], d[1], t, u);
    if (!have_params) {
      const Expansion det1 =
          split_orient_det(a[0], a[1], b[0], b[1], c[0], c[1]);
      const Expansion det2 =
          split_orient_det(a[0], a[1], b[0], b[1], d[0], d[1]);
      const Expansion det3 =
          split_orient_det(c[0], c[1], d[0], d[1], a[0], a[1]);
      const Expansion det4 =
          split_orient_det(c[0], c[1], d[0], d[1], b[0], b[1]);

      have_params = split_exact_params(det1, det2, det3, det4, t, u);
      if (!have_params) {
        const double max_coord = std::max(
            {std::fabs(a[0]), std::fabs(a[1]), std::fabs(b[0]), std::fabs(b[1]),
             std::fabs(c[0]), std::fabs(c[1]), std::fabs(d[0]), std::fabs(d[1])});
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

        const Expansion sdet1 = split_orient_det(sax, say, sbx, sby, scx, scy);
        const Expansion sdet2 = split_orient_det(sax, say, sbx, sby, sdx, sdy);
        const Expansion sdet3 = split_orient_det(scx, scy, sdx, sdy, sax, say);
        const Expansion sdet4 = split_orient_det(scx, scy, sdx, sdy, sbx, sby);
        have_params = split_exact_params(sdet1, sdet2, sdet3, sdet4, t, u);
      }
    }

    if (!have_params) {
      // No representable split point: keep both input segments unsplit (the
      // caller still records the pair) and make the failure visible.
      max_residual = std::numeric_limits<long double>::infinity();
      return true;
    }
    t = split_clamp_parameter(t);
    u = split_clamp_parameter(u);

    // Evaluate the contact from both parameterizations and keep the candidate
    // closer to both supporting lines. One of them is exact whenever the
    // intersection is representable (for example the 2e200 x 1e-160 case is
    // exact on the short segment but rounded on the long one).
    const Point2d contact_a = split_point_at(a, b, t);
    const Point2d contact_b = split_point_at(c, d, u);
    const long double error_a = split_contact_error(contact_a, a, b, c, d);
    const long double error_b = split_contact_error(contact_b, a, b, c, d);
    const Point2d contact = error_b < error_a ? contact_b : contact_a;

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
//   - endpoint touches, collinear overlaps and proper crossings split both
//     segments and are recorded once in `input_pairs` as (i, j), i < j, sorted;
//     a crossing whose parameter ratio is not representable leaves both
//     segments unsplit (it is still recorded) and sets `max_residual` to +inf;
//   - `max_residual` is the largest perpendicular distance of a computed
//     intersection point to the two supporting lines; it is +infinity when a
//     crossing is not representable (see the supported domain below).
//
// Candidate pairs come from a min-x sweep with an active set; a cheap
// double-precision same-side test rejects the near-parallel disjoint pairs that
// dominate a dense ring before the exact classification runs. The exact
// classification uses sign_orient2d, so disjoint and exact endpoint cases never
// depend on a tolerance. Pairs that share an exact endpoint are recorded without
// further work unless the two segments are collinear. Proper-crossing parameters
// are determinant ratios: they are evaluated in double-double first (about 106
// bits, no allocation), then in exact Expansions in the unscaled frame, and
// finally in a 2^k-scaled frame. A ratio that rounds to an endpoint is rejected
// rather than emitted. The supported domain is finite coordinates whose input exponent
// span is at most detail::kMaxSpan (1314). Within it, differences and products
// are formed exactly as Expansions, so mixed spans such as 2e200 and 1e-160 are
// handled. Parameters are materialized as double regardless of the width of
// long double, so a ratio such as 1e-360 that underflows to zero is treated the
// same on every platform: the pair is recorded but both inputs stay unsplit and
// `max_residual` is set to +infinity ("split point not representable") rather
// than silently emitting a wrong split. Inputs beyond the span are outside the
// contract.
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
    const ArrangementSegment& source = inputs[i].seg;
    std::sort(ps.begin(), ps.end(), [&source](const detail::SplitParam& lhs,
                                              const detail::SplitParam& rhs) {
      if (lhs.t != rhs.t) return lhs.t < rhs.t;
      // Contacts that round to the same parameter are ordered along the
      // segment. On a 1e308-long segment two contacts 1 apart share t exactly,
      // so t alone cannot order them; the dominant coordinate is monotone.
      return detail::split_contact_before(source.a, source.b, lhs.p, rhs.p);
    });
    ps.erase(std::unique(ps.begin(), ps.end(),
                         [](const detail::SplitParam& lhs, const detail::SplitParam& rhs) {
                           return lhs.p[0] == rhs.p[0] && lhs.p[1] == rhs.p[1];
                         }),
             ps.end());

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

// |ma| * 2^ea < |mb| * 2^eb. A zero-area face (mantissa 0) has no
// scale and sorts below every nonzero face.
inline bool area_abs_less(double ma, int ea, double mb, int eb) noexcept {
  const bool zero_a = (ma == 0.0);
  const bool zero_b = (mb == 0.0);
  if (zero_a || zero_b) return zero_a && !zero_b;
  if (ea != eb) return ea < eb;
  return std::fabs(ma) < std::fabs(mb);
}

// Signed comparison of ma*2^ea and mb*2^eb: -1, 0 or 1.
inline int area_compare(double ma, int ea, double mb, int eb) noexcept {
  if (ma == 0.0 || mb == 0.0) {
    if (ma == 0.0 && mb == 0.0) return 0;
    if (ma == 0.0) return (mb < 0.0) ? 1 : -1;
    return (ma < 0.0) ? -1 : 1;
  }
  const bool negative_a = ma < 0.0;
  const bool negative_b = mb < 0.0;
  if (negative_a != negative_b) return negative_a ? -1 : 1;
  int magnitude = 0;
  if (ea != eb) {
    magnitude = ea < eb ? -1 : 1;
  } else {
    const double aa = std::fabs(ma);
    const double ab = std::fabs(mb);
    magnitude = aa < ab ? -1 : (aa > ab ? 1 : 0);
  }
  return negative_a ? -magnitude : magnitude;
}

}  // namespace detail

// A planar arrangement (DCEL) of closed input chains with per-source winding.
//
// Build contract:
//   - The input segments are split at every intersection first
//     (`detail::split_segments`), so the resulting pieces meet only at shared
//     endpoints; those endpoints are merged into vertices by exact coordinate
//     equality.
//   - An undirected edge aggregates the signed contribution of every source that
//     traverses it: +1 when a source traverses it from the lower to the higher
//     vertex id, -1 in the opposite direction. Duplicate and reversed pieces
//     therefore accumulate instead of overwriting each other.
//   - Half-edges are linked with the standard face-on-left convention:
//     `next(h)` is the clockwise predecessor of `twin(h)` in the angular order
//     of the outgoing half-edges at the target vertex of `h`.
//   - Every closed boundary walk is one `Face`. Its signed shoelace area is
//     positive for a bounded face and negative for the unbounded side. A face
//     with several boundary rings (a hole, or a component nested inside another
//     component) contributes one `Face` record per ring; all of those records
//     carry the same winding. Point location only considers bounded faces, but
//     the boolean layer classifies every face record by its per-source winding:
//     for a hole the region lies on the unbounded side, so the kept record is
//     the unbounded one. Boundary emission is therefore seeded from every
//     accepted record, bounded or not.
//   - Per-source winding follows w[left] = w[right] + s * delta_low_to_high, with
//     s = +1 for the low->high half-edge. The unbounded side of a component
//     starts at zero; a component nested inside another one inherits the winding
//     of the face that contains it. `winding(face, source)` returns that one
//     source's value; `Face::winding` is the sum over all sources.
//   - Inputs must be closed chains (the segments of one source must join end to
//     end); otherwise the winding numbers are only defined relative to the
//     initialized outer faces.
//   - Face areas are computed after an exact power-of-two normalization of the
//     coordinates (e = ilogb of the largest |coordinate|), so faces are
//     classified and ordered correctly from subnormal magnitudes up to about
//     1e308, provided every face's area stays representable after that scaling
//     (feature-size dynamic range below roughly 1e154). Ordering a 1e308 face
//     against a 1e-200 nested face would need per-face exponents; that mixed
//     scale is outside this contract. locate() uses the original coordinates
//     and exact predicates.
class Arrangement2 {
 public:
  struct HalfEdge {
    std::uint32_t origin = 0;
    std::uint32_t twin = 0;
    std::uint32_t next = 0;
    std::uint32_t face = 0;
  };

  struct Face {
    std::uint32_t edge = 0;
    std::int32_t winding = 0;
    bool bounded = false;
  };

  // Sentinel returned by locate() when the point is outside every bounded face
  // or lies exactly on an edge.
  static constexpr std::uint32_t kNoFace = 0xFFFFFFFFu;

  // Builds the arrangement from closed input chains. Intersections are split
  // before the DCEL is assembled, so no two pieces cross in their interiors.
  static Arrangement2 build(std::span<const ArrangementSegment> segments) {
    Arrangement2 arrangement;
    const detail::SplitResult split = detail::split_segments(segments);
    arrangement.input_pairs_ = split.input_pairs;
    arrangement.max_residual_ = split.max_residual;

    if (!segments.empty()) {
      std::uint32_t max_source = 0;
      for (const ArrangementSegment& segment : segments) {
        if (segment.source > max_source) max_source = segment.source;
      }
      arrangement.source_count_ = max_source + 1;
    }

    // 1. Merge the split endpoints into vertices by exact coordinate equality.
    // The endpoints below are exact doubles produced by the split stage, so
    // bitwise equality is the right identity and no epsilon is involved.
    std::vector<Point2d> points;
    points.reserve(split.segments.size() * 2);
    for (const detail::SplitSegment& piece : split.segments) {
      points.push_back(piece.a);
      points.push_back(piece.b);
    }
    const auto point_less = [](const Point2d& lhs, const Point2d& rhs) noexcept {
      if (lhs[0] != rhs[0]) return lhs[0] < rhs[0];
      return lhs[1] < rhs[1];
    };
    std::sort(points.begin(), points.end(), point_less);
    points.erase(std::unique(points.begin(), points.end(),
                             [](const Point2d& lhs, const Point2d& rhs) noexcept {
                               return lhs[0] == rhs[0] && lhs[1] == rhs[1];
                             }),
                 points.end());
    arrangement.vertices_ = std::move(points);

    const auto vertex_id = [&arrangement, &point_less](const Point2d& p) noexcept {
      const auto it = std::lower_bound(arrangement.vertices_.begin(),
                                       arrangement.vertices_.end(), p, point_less);
      return static_cast<std::uint32_t>(it - arrangement.vertices_.begin());
    };

    // 2. Aggregate undirected edges, keeping the signed low->high delta of every
    // source. A source may traverse the same edge twice (shared or duplicated
    // pieces); the deltas add up so the accumulated winding stays exact.
    struct EdgePiece {
      std::uint32_t low = 0;
      std::uint32_t high = 0;
      std::uint32_t source = 0;
      std::int32_t delta = 0;
    };
    std::vector<EdgePiece> pieces;
    pieces.reserve(split.segments.size());
    for (const detail::SplitSegment& piece : split.segments) {
      const std::uint32_t u = vertex_id(piece.a);
      const std::uint32_t v = vertex_id(piece.b);
      if (u == v) continue;
      const bool forward = u < v;
      EdgePiece edge_piece;
      edge_piece.low = forward ? u : v;
      edge_piece.high = forward ? v : u;
      edge_piece.source = piece.source;
      edge_piece.delta = forward ? 1 : -1;
      pieces.push_back(edge_piece);
    }
    std::sort(pieces.begin(), pieces.end(),
              [](const EdgePiece& lhs, const EdgePiece& rhs) noexcept {
                if (lhs.low != rhs.low) return lhs.low < rhs.low;
                if (lhs.high != rhs.high) return lhs.high < rhs.high;
                return lhs.source < rhs.source;
              });

    struct EdgeRecord {
      std::uint32_t low = 0;
      std::uint32_t high = 0;
      std::vector<std::pair<std::uint32_t, std::int32_t>> deltas;
    };
    std::vector<EdgeRecord> edges;
    for (std::size_t i = 0; i < pieces.size();) {
      const std::uint32_t low = pieces[i].low;
      const std::uint32_t high = pieces[i].high;
      EdgeRecord record;
      record.low = low;
      record.high = high;
      while (i < pieces.size() && pieces[i].low == low && pieces[i].high == high) {
        const std::uint32_t source = pieces[i].source;
        std::int32_t delta = 0;
        while (i < pieces.size() && pieces[i].low == low &&
               pieces[i].high == high && pieces[i].source == source) {
          delta += pieces[i].delta;
          ++i;
        }
        if (delta != 0) record.deltas.emplace_back(source, delta);
      }
      edges.push_back(std::move(record));
    }

    std::vector<std::uint32_t> sources;
    sources.reserve(segments.size());
    for (const ArrangementSegment& segment : segments) {
      sources.push_back(segment.source);
    }
    std::sort(sources.begin(), sources.end());
    sources.erase(std::unique(sources.begin(), sources.end()), sources.end());
    arrangement.sources_ = sources;
    const std::size_t source_total = arrangement.sources_.size();

    for (EdgeRecord& record : edges) {
      for (std::pair<std::uint32_t, std::int32_t>& delta : record.deltas) {
        const auto it = std::lower_bound(arrangement.sources_.begin(),
                                         arrangement.sources_.end(), delta.first);
        delta.first = static_cast<std::uint32_t>(it - arrangement.sources_.begin());
      }
    }

    // 3. One half-edge pair per undirected edge, with the outgoing half-edges
    // of each vertex sorted counter-clockwise by angle.
    const std::size_t edge_total = edges.size();
    arrangement.half_edges_.resize(edge_total * 2);
    for (std::size_t e = 0; e < edge_total; ++e) {
      HalfEdge& low = arrangement.half_edges_[e * 2];
      HalfEdge& high = arrangement.half_edges_[e * 2 + 1];
      low.origin = edges[e].low;
      low.twin = static_cast<std::uint32_t>(e * 2 + 1);
      high.origin = edges[e].high;
      high.twin = static_cast<std::uint32_t>(e * 2);
    }

    const std::uint32_t half_edge_total =
        static_cast<std::uint32_t>(arrangement.half_edges_.size());
    std::vector<std::vector<std::uint32_t>> outgoing(arrangement.vertices_.size());
    for (std::uint32_t h = 0; h < half_edge_total; ++h) {
      outgoing[arrangement.half_edges_[h].origin].push_back(h);
    }

    // Upper half-plane group first (positive x axis included), then the exact
    // orientation predicate. The predicate is exact, so the order is a strict
    // weak ordering without any epsilon; equal directions fall back to the
    // target vertex id and the half-edge id.
    const auto half_plane = [](double dx, double dy) noexcept {
      return (dy > 0.0 || (dy == 0.0 && dx > 0.0)) ? 0 : 1;
    };
    for (std::size_t v = 0; v < arrangement.vertices_.size(); ++v) {
      const Point2d origin = arrangement.vertices_[v];
      std::vector<std::uint32_t>& ring = outgoing[v];
      std::sort(ring.begin(), ring.end(), [&](std::uint32_t lhs, std::uint32_t rhs) {
        const std::uint32_t lhs_target =
            arrangement.half_edges_[arrangement.half_edges_[lhs].twin].origin;
        const std::uint32_t rhs_target =
            arrangement.half_edges_[arrangement.half_edges_[rhs].twin].origin;
        const Point2d& a = arrangement.vertices_[lhs_target];
        const Point2d& b = arrangement.vertices_[rhs_target];
        const int lhs_half = half_plane(a[0] - origin[0], a[1] - origin[1]);
        const int rhs_half = half_plane(b[0] - origin[0], b[1] - origin[1]);
        if (lhs_half != rhs_half) return lhs_half < rhs_half;
        const int orientation =
            sign_orient2d(origin[0], origin[1], a[0], a[1], b[0], b[1]);
        if (orientation != 0) return orientation > 0;
        if (lhs_target != rhs_target) return lhs_target < rhs_target;
        return lhs < rhs;
      });
    }

    std::vector<std::uint32_t> position(half_edge_total);
    for (std::size_t v = 0; v < outgoing.size(); ++v) {
      for (std::size_t i = 0; i < outgoing[v].size(); ++i) {
        position[outgoing[v][i]] = static_cast<std::uint32_t>(i);
      }
    }

    // 4. next(h): at the target vertex of h, the clockwise predecessor of
    // twin(h) in the outgoing angular order (standard face-on-left DCEL).
    for (std::uint32_t h = 0; h < half_edge_total; ++h) {
      const std::uint32_t twin = arrangement.half_edges_[h].twin;
      const std::uint32_t target = arrangement.half_edges_[twin].origin;
      const std::vector<std::uint32_t>& ring = outgoing[target];
      const std::uint32_t count = static_cast<std::uint32_t>(ring.size());
      const std::uint32_t index = position[twin];
      arrangement.half_edges_[h].next = ring[(index + count - 1) % count];
    }

    // 5. Walk next() to assign one Face per closed boundary cycle and to
    // accumulate the signed shoelace area. Each cycle is scaled by its own
    // exact power of two (e = ilogb(max |coordinate| of that cycle)): long
    // double is double on MSVC, so a raw shoelace would overflow to
    // inf - inf = NaN around 1e200 and underflow to 0 around 1e-200, and one
    // global frame collapses a unit-scale ring next to a 1e308 ring to zero
    // area. The normalized area is stored as mantissa * 2^exponent so
    // `mantissa > 0` gives the face orientation and the |area| ordering used by
    // locate() and by the nesting seed stays exact across mixed scales. Point
    // location keeps operating on the original coordinates.
    std::vector<std::uint32_t> face_of(half_edge_total, kNoFace);
    std::vector<std::uint32_t> cycle;
    cycle.reserve(half_edge_total);
    for (std::uint32_t start = 0; start < half_edge_total; ++start) {
      if (face_of[start] != kNoFace) continue;
      const std::uint32_t face_id =
          static_cast<std::uint32_t>(arrangement.faces_.size());
      cycle.clear();
      double max_magnitude = 0.0;
      std::uint32_t h = start;
      std::uint32_t guard = 0;
      do {
        face_of[h] = face_id;
        arrangement.half_edges_[h].face = face_id;
        const Point2d& p = arrangement.vertices_[arrangement.half_edges_[h].origin];
        const double x = std::fabs(p[0]);
        const double y = std::fabs(p[1]);
        if (x > max_magnitude) max_magnitude = x;
        if (y > max_magnitude) max_magnitude = y;
        cycle.push_back(h);
        h = arrangement.half_edges_[h].next;
        ++guard;
      } while (h != start && guard <= half_edge_total);

      int area_exponent = 0;
      if (max_magnitude > 0.0 && std::isfinite(max_magnitude)) {
        area_exponent = std::ilogb(max_magnitude);
      }
      long double twice_area = 0.0L;
      for (const std::uint32_t edge : cycle) {
        const Point2d& a = arrangement.vertices_[arrangement.half_edges_[edge].origin];
        const Point2d& b = arrangement.vertices_[
            arrangement.half_edges_[arrangement.half_edges_[edge].twin].origin];
        const double ax = std::ldexp(a[0], -area_exponent);
        const double ay = std::ldexp(a[1], -area_exponent);
        const double bx = std::ldexp(b[0], -area_exponent);
        const double by = std::ldexp(b[1], -area_exponent);
        twice_area += static_cast<long double>(ax) * static_cast<long double>(by) -
                      static_cast<long double>(bx) * static_cast<long double>(ay);
      }
      const long double area = twice_area * 0.5L;
      double mantissa = 0.0;
      int mantissa_exponent = 0;
      if (area != 0.0L) {
        mantissa = static_cast<double>(std::frexp(area, &mantissa_exponent));
      }
      arrangement.face_area_mantissa_.push_back(mantissa);
      arrangement.face_area_exponent_.push_back(
          area != 0.0L ? mantissa_exponent + 2 * area_exponent : 0);
      arrangement.faces_.push_back(Face{start, 0, mantissa > 0.0});
    }

    // 6. Propagate the per-source winding across every edge. Each connected
    // component seeds the winding of its unbounded face: zero when it is the
    // outermost component covering that spot, otherwise the winding of the
    // face that contains it. Components are processed from the largest outer
    // boundary to the smallest so a container is always known first.
    const std::uint32_t face_total =
        static_cast<std::uint32_t>(arrangement.faces_.size());
    arrangement.windings_.assign(
        static_cast<std::size_t>(face_total) * source_total, 0);
    if (face_total != 0 && source_total != 0) {
      std::vector<std::int32_t>& windings = arrangement.windings_;

      std::vector<std::uint32_t> parent(face_total);
      for (std::uint32_t f = 0; f < face_total; ++f) parent[f] = f;
      const auto find_root = [&parent](std::uint32_t f) noexcept {
        while (parent[f] != f) {
          parent[f] = parent[parent[f]];
          f = parent[f];
        }
        return f;
      };
      for (std::size_t e = 0; e < edge_total; ++e) {
        const std::uint32_t left = face_of[e * 2];
        const std::uint32_t right = face_of[e * 2 + 1];
        const std::uint32_t left_root = find_root(left);
        const std::uint32_t right_root = find_root(right);
        if (left_root == right_root) continue;
        if (left_root < right_root) {
          parent[right_root] = left_root;
        } else {
          parent[left_root] = right_root;
        }
      }

      std::vector<std::vector<std::uint32_t>> adjacency(face_total);
      for (std::size_t e = 0; e < edge_total; ++e) {
        const std::uint32_t left = face_of[e * 2];
        const std::uint32_t right = face_of[e * 2 + 1];
        if (left == right) continue;
        adjacency[left].push_back(static_cast<std::uint32_t>(e));
        adjacency[right].push_back(static_cast<std::uint32_t>(e));
      }

      std::vector<std::vector<std::uint32_t>> members(face_total);
      for (std::uint32_t f = 0; f < face_total; ++f) {
        members[find_root(f)].push_back(f);
      }

      std::vector<std::uint32_t> seed_face(face_total, kNoFace);
      std::vector<std::uint32_t> seed_vertex(face_total, kNoFace);
      std::vector<std::uint32_t> component_order;
      for (std::uint32_t root = 0; root < face_total; ++root) {
        if (members[root].empty()) continue;
        std::uint32_t min_vertex = kNoFace;
        for (const std::uint32_t f : members[root]) {
          const std::uint32_t start = arrangement.faces_[f].edge;
          std::uint32_t h = start;
          std::uint32_t guard = 0;
          do {
            const std::uint32_t origin = arrangement.half_edges_[h].origin;
            if (min_vertex == kNoFace || origin < min_vertex) min_vertex = origin;
            h = arrangement.half_edges_[h].next;
            ++guard;
          } while (h != start && guard <= half_edge_total);
        }
        // The seed is the unbounded side of this component: the incident face
        // record with the most negative signed area. Each closed ring has one
        // positive (bounded) and one negative (unbounded) cycle, and at the
        // lexicographically smallest vertex both cycles touch one of the
        // outgoing edges. Selecting by signed area is robust where the previous
        // angular rule was not: `outgoing` is ordered [0..pi, pi..2pi], so for
        // a leftmost vertex whose edges straddle the +x axis `.back()` is the
        // lower-half-plane edge and its left side is the interior. Restricting
        // candidates to this component also keeps touching components apart.
        std::uint32_t seed = kNoFace;
        for (const std::uint32_t h : outgoing[min_vertex]) {
          const std::uint32_t sides[2] = {face_of[h],
                                          face_of[arrangement.half_edges_[h].twin]};
          for (const std::uint32_t candidate : sides) {
            if (find_root(candidate) != root) continue;
            if (seed == kNoFace) {
              seed = candidate;
              continue;
            }
            const int order = detail::area_compare(
                arrangement.face_area_mantissa_[candidate],
                arrangement.face_area_exponent_[candidate],
                arrangement.face_area_mantissa_[seed],
                arrangement.face_area_exponent_[seed]);
            if (order < 0 || (order == 0 && candidate < seed)) {
              seed = candidate;
            }
          }
        }
        seed_face[root] = seed;
        seed_vertex[root] = min_vertex;
        component_order.push_back(root);
      }

      std::sort(component_order.begin(), component_order.end(),
                [&](std::uint32_t lhs, std::uint32_t rhs) noexcept {
                  const std::uint32_t lhs_seed = seed_face[lhs];
                  const std::uint32_t rhs_seed = seed_face[rhs];
                  const int order = detail::area_compare(
                      std::fabs(arrangement.face_area_mantissa_[lhs_seed]),
                      arrangement.face_area_exponent_[lhs_seed],
                      std::fabs(arrangement.face_area_mantissa_[rhs_seed]),
                      arrangement.face_area_exponent_[rhs_seed]);
                  if (order != 0) return order > 0;
                  return lhs < rhs;
                });

      std::vector<char> assigned(face_total, 0);
      std::vector<std::int32_t> seed_values(source_total, 0);
      std::vector<std::uint32_t> queue;
      for (const std::uint32_t root : component_order) {
        const std::uint32_t seed = seed_face[root];
        const Point2d representative = arrangement.vertices_[seed_vertex[root]];
        double best_mantissa = 0.0;
        int best_exponent = 0;
        bool have_best = false;
        for (std::size_t s = 0; s < source_total; ++s) seed_values[s] = 0;
        for (std::uint32_t f = 0; f < face_total; ++f) {
          if (assigned[f] == 0) continue;
          if (!(arrangement.face_area_mantissa_[f] > 0.0)) continue;
          if (have_best &&
              !detail::area_abs_less(arrangement.face_area_mantissa_[f],
                                     arrangement.face_area_exponent_[f],
                                     best_mantissa, best_exponent)) {
            continue;
          }
          if (!arrangement.point_in_face(f, representative)) continue;
          have_best = true;
          best_mantissa = arrangement.face_area_mantissa_[f];
          best_exponent = arrangement.face_area_exponent_[f];
          for (std::size_t s = 0; s < source_total; ++s) {
            seed_values[s] = windings[static_cast<std::size_t>(f) * source_total + s];
          }
        }
        for (std::size_t s = 0; s < source_total; ++s) {
          windings[static_cast<std::size_t>(seed) * source_total + s] = seed_values[s];
        }
        assigned[seed] = 1;
        queue.clear();
        queue.push_back(seed);
        while (!queue.empty()) {
          const std::uint32_t face = queue.back();
          queue.pop_back();
          for (const std::uint32_t e : adjacency[face]) {
            const std::uint32_t left = face_of[e * 2];
            const std::uint32_t right = face_of[e * 2 + 1];
            const std::uint32_t other = left == face ? right : left;
            if (assigned[other] != 0) continue;
            const std::int32_t* base =
                &windings[static_cast<std::size_t>(face) * source_total];
            std::int32_t* dest =
                &windings[static_cast<std::size_t>(other) * source_total];
            for (std::size_t s = 0; s < source_total; ++s) dest[s] = base[s];
            // w[left] = w[right] + delta_low_to_high, with `left` and `right`
            // taken from the low->high half-edge of the edge.
            for (const std::pair<std::uint32_t, std::int32_t>& delta : edges[e].deltas) {
              if (left == face) {
                dest[delta.first] -= delta.second;
              } else {
                dest[delta.first] += delta.second;
              }
            }
            assigned[other] = 1;
            queue.push_back(other);
          }
        }
      }

      for (std::uint32_t f = 0; f < face_total; ++f) {
        std::int64_t sum = 0;
        for (std::size_t s = 0; s < source_total; ++s) {
          sum += windings[static_cast<std::size_t>(f) * source_total + s];
        }
        arrangement.faces_[f].winding = static_cast<std::int32_t>(sum);
      }
    }

    return arrangement;
  }

  std::uint32_t vertex_count() const noexcept {
    return static_cast<std::uint32_t>(vertices_.size());
  }
  std::uint32_t edge_count() const noexcept {
    return static_cast<std::uint32_t>(half_edges_.size() / 2);
  }
  std::uint32_t face_count() const noexcept {
    return static_cast<std::uint32_t>(faces_.size());
  }
  std::uint32_t source_count() const noexcept { return source_count_; }

  std::span<const Point2d> vertices() const noexcept { return vertices_; }
  std::span<const HalfEdge> half_edges() const noexcept { return half_edges_; }
  std::span<const Face> faces() const noexcept { return faces_; }

  // Winding number of `face` contributed by `source`. The two-argument call is
  // always per source, which is what a boolean fill rule needs; it never depends
  // on how many other sources cover the same face. Face::winding holds the
  // aggregate sum over all sources. A kNoFace face id, any face id beyond
  // faces().size(), or a source that never appears in the input yields 0.
  std::int32_t winding(std::uint32_t face, std::uint32_t source) const noexcept {
    if (face >= faces_.size()) return 0;
    const auto it = std::lower_bound(sources_.begin(), sources_.end(), source);
    if (it == sources_.end() || *it != source) return 0;
    const std::size_t index = static_cast<std::size_t>(it - sources_.begin());
    return windings_[static_cast<std::size_t>(face) * sources_.size() + index];
  }

  // Alias for winding(); both return the per-source contribution. Kept so the
  // boolean layer can spell out that it wants one operand at a time.
  std::int32_t winding_of_source(std::uint32_t face, std::uint32_t source) const noexcept {
    return winding(face, source);
  }

  // Face containing `p`, or kNoFace when `p` is on an edge or outside every
  // bounded face. Among the bounded faces that contain the point the innermost
  // one (smallest |area|) wins.
  std::uint32_t locate(const Point2d& p) const noexcept {
    const std::size_t edge_total = half_edges_.size() / 2;
    for (std::size_t e = 0; e < edge_total; ++e) {
      const Point2d& a = vertices_[half_edges_[e * 2].origin];
      const Point2d& b = vertices_[half_edges_[e * 2 + 1].origin];
      if (sign_orient2d(a[0], a[1], b[0], b[1], p[0], p[1]) == 0 &&
          detail::split_in_box(a, b, p)) {
        return kNoFace;
      }
    }

    std::uint32_t best = kNoFace;
    double best_mantissa = 0.0;
    int best_exponent = 0;
    const std::uint32_t face_total = static_cast<std::uint32_t>(faces_.size());
    for (std::uint32_t f = 0; f < face_total; ++f) {
      if (!faces_[f].bounded) continue;
      if (best != kNoFace &&
          !detail::area_abs_less(face_area_mantissa_[f], face_area_exponent_[f],
                                 best_mantissa, best_exponent)) {
        continue;
      }
      if (!point_in_face(f, p)) continue;
      best_mantissa = face_area_mantissa_[f];
      best_exponent = face_area_exponent_[f];
      best = f;
    }
    return best;
  }

  std::span<const std::pair<std::uint32_t, std::uint32_t>> input_pairs() const noexcept {
    return input_pairs_;
  }
  double max_residual() const noexcept { return max_residual_; }

 private:
  // Exact right-ray crossing test for the even-odd rule. A vertex at the ray's
  // height belongs to the edge that leaves upward, which is the usual half-open
  // convention and keeps shared vertices from being counted twice.
  static bool crosses_right(const Point2d& a, const Point2d& b,
                            const Point2d& p) noexcept {
    const bool a_above = a[1] > p[1];
    const bool b_above = b[1] > p[1];
    if (a_above == b_above) return false;
    const int orientation = sign_orient2d(a[0], a[1], b[0], b[1], p[0], p[1]);
    return b[1] > a[1] ? orientation > 0 : orientation < 0;
  }

  bool point_in_face(std::uint32_t face, const Point2d& p) const noexcept {
    bool inside = false;
    const std::uint32_t start = faces_[face].edge;
    const std::uint32_t limit = static_cast<std::uint32_t>(half_edges_.size());
    std::uint32_t h = start;
    std::uint32_t guard = 0;
    do {
      const std::uint32_t target = half_edges_[half_edges_[h].twin].origin;
      if (crosses_right(vertices_[half_edges_[h].origin], vertices_[target], p)) {
        inside = !inside;
      }
      h = half_edges_[h].next;
      ++guard;
    } while (h != start && guard <= limit);
    return inside;
  }

  std::vector<Point2d> vertices_;
  std::vector<HalfEdge> half_edges_;
  std::vector<Face> faces_;
  // Signed face-cycle area as mantissa * 2^exponent after the exact per-cycle
  // power-of-two normalization in build(); only the sign and the |area|
  // ordering are meaningful.
  std::vector<double> face_area_mantissa_;
  std::vector<int> face_area_exponent_;
  std::vector<std::pair<std::uint32_t, std::uint32_t>> input_pairs_;
  // Distinct input sources, sorted; `windings_` is a faces x sources table.
  std::vector<std::uint32_t> sources_;
  std::vector<std::int32_t> windings_;
  std::uint32_t source_count_ = 0;
  double max_residual_ = 0.0;
};

}  // namespace v1
}  // namespace geom
