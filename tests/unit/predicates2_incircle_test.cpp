#include <doctest/doctest.h>

#include <algorithm>
#include <array>
#include <cstdint>
#include <cstdio>
#include <cstring>
#include <limits>
#include <random>

#if defined(_MSC_VER)
#include <intrin.h>
#endif

#include <geom/numeric/predicates2.hpp>

using geom::v1::sign_incircle;

namespace {

// Coordinates are small integers, so the lifted determinant fits in int64.
int incircle_oracle(std::int64_t ax, std::int64_t ay, std::int64_t bx, std::int64_t by,
                    std::int64_t cx, std::int64_t cy, std::int64_t dx, std::int64_t dy) {
  const std::int64_t adx = ax - dx, ady = ay - dy;
  const std::int64_t bdx = bx - dx, bdy = by - dy;
  const std::int64_t cdx = cx - dx, cdy = cy - dy;
  const std::int64_t alift = adx * adx + ady * ady;
  const std::int64_t blift = bdx * bdx + bdy * bdy;
  const std::int64_t clift = cdx * cdx + cdy * cdy;
  const std::int64_t det =
      alift * (bdx * cdy - cdx * bdy) - blift * (adx * cdy - cdx * ady) +
      clift * (adx * bdy - bdx * ady);
  return det > 0 ? 1 : (det < 0 ? -1 : 0);
}

}  // namespace

TEST_CASE("incircle sign convention") {
  // CCW triangle (0,0), (1,0), (0,1): center (0.5,0.5), point (0.5,0.5) is inside.
  CHECK(sign_incircle(0, 0, 1, 0, 0, 1, 0.5, 0.5) > 0);
  CHECK(sign_incircle(0, 0, 1, 0, 0, 1, 5.0, 5.0) < 0);
}

TEST_CASE("incircle agrees with the integer oracle") {
  std::mt19937 rng(7654321u);
  std::uniform_int_distribution<int> dist(-512, 511);
  for (int i = 0; i < 100000; ++i) {
    const int ax = dist(rng), ay = dist(rng);
    const int bx = dist(rng), by = dist(rng);
    const int cx = dist(rng), cy = dist(rng);
    const int dx = dist(rng), dy = dist(rng);
    REQUIRE(sign_incircle(ax, ay, bx, by, cx, cy, dx, dy) ==
            incircle_oracle(ax, ay, bx, by, cx, cy, dx, dy));
  }
}

TEST_CASE("cocircular points are exactly degenerate") {
  CHECK(sign_incircle(1, 0, 0, 1, -1, 0, 0, -1) == 0);
  CHECK(sign_incircle(0, 0, 2, 0, 0, 2, 2, 2) == 0);
}

// ---------------------------------------------------------------------------
// Cases appended after the brief's tests: non-finite input, degenerate input,
// subnormal and huge coordinates, mixed magnitude, and the exact path just off
// the circle. Expected signs are exact values computed in rational arithmetic
// over the actual double inputs.
// ---------------------------------------------------------------------------

TEST_CASE("incircle rejects non-finite coordinates in every position") {
  const double nan = std::numeric_limits<double>::quiet_NaN();
  const double inf = std::numeric_limits<double>::infinity();
  CHECK(sign_incircle(nan, 0, 1, 0, 0, 1, 0.5, 0.5) == 0);
  CHECK(sign_incircle(0, nan, 1, 0, 0, 1, 0.5, 0.5) == 0);
  CHECK(sign_incircle(0, 0, nan, 0, 0, 1, 0.5, 0.5) == 0);
  CHECK(sign_incircle(0, 0, 1, nan, 0, 1, 0.5, 0.5) == 0);
  CHECK(sign_incircle(0, 0, 1, 0, nan, 1, 0.5, 0.5) == 0);
  CHECK(sign_incircle(0, 0, 1, 0, 0, nan, 0.5, 0.5) == 0);
  CHECK(sign_incircle(0, 0, 1, 0, 0, 1, nan, 0.5) == 0);
  CHECK(sign_incircle(0, 0, 1, 0, 0, 1, 0.5, nan) == 0);
  CHECK(sign_incircle(inf, 0, 1, 0, 0, 1, 0.5, 0.5) == 0);
  CHECK(sign_incircle(0, 0, 1, 0, 0, -inf, 0.5, 0.5) == 0);
  CHECK(sign_incircle(0, 0, 1, 0, 0, 1, 0.5, inf) == 0);
}

TEST_CASE("incircle degenerate configurations are exactly zero") {
  CHECK(sign_incircle(0, 0, 0, 0, 0, 0, 0, 0) == 0);
  CHECK(sign_incircle(1, 1, 1, 1, 2, 3, 0.5, 0.5) == 0);  // repeated vertex
  CHECK(sign_incircle(0, 0, 1, 1, 2, 2, 3, 3) == 0);      // collinear
  CHECK(sign_incircle(0, 0, 4, 0, 0, 4, 4, 4) == 0);      // cocircular
  // Clockwise triangle: an inside point is negative.
  CHECK(sign_incircle(0, 0, 0, 1, 1, 0, 0.5, 0.5) < 0);
}

TEST_CASE("incircle handles subnormal and huge coordinates") {
  const double s = std::ldexp(1.0, -1070);  // subnormal, exactly representable
  CHECK(sign_incircle(0, 0, s, 0, 0, s, 0.25 * s, 0.25 * s) > 0);
  CHECK(sign_incircle(0, 0, 0, s, s, 0, 0.25 * s, 0.25 * s) < 0);
  CHECK(sign_incircle(0, 0, s, 0, 0, s, s, s) == 0);

  const double h = std::ldexp(1.0, 600);  // lifts overflow without rescaling
  CHECK(sign_incircle(0, 0, h, 0, 0, h, h / 2, h / 2) > 0);
  CHECK(sign_incircle(0, 0, h, 0, 0, h, 1.5 * h, 1.5 * h) < 0);
  CHECK(sign_incircle(0, 0, h, 0, 0, h, h, h) == 0);
}

TEST_CASE("incircle resolves mixed-magnitude coordinates") {
  // The determinant spans 800 decades, but the lifted products stay finite, so
  // the unscaled filter decides without any rescaling.
  const double big = std::ldexp(1.0, 400);
  const double small = std::ldexp(1.0, -400);
  CHECK(sign_incircle(0, 0, big, 0, 0, small, big / 2, small / 2) > 0);
  CHECK(sign_incircle(0, 0, big, 0, 0, small, 2.0 * big, 2.0 * small) < 0);
  CHECK(sign_incircle(0, 0, big, 0, 0, small, big, small) == 0);

  // The 1e200 lifts overflow, so the predicate rescales by a power of two; the
  // 1e100 coordinates survive the rescale and the sign is still recovered.
  CHECK(sign_incircle(0, 0, 1e200, 0, 0, 1e100, 5e199, 5e99) > 0);
  CHECK(sign_incircle(0, 0, 1e200, 0, 0, 1e100, 2e200, 2e100) < 0);
  CHECK(sign_incircle(0, 0, 1e200, 0, 0, 1e100, 1e200, 1e100) == 0);

  // 1e200 over 1e-200 spans about 2^1330 binary exponents. The exact determinants
  // here are about 2^1325 and 2^1330, both far above DBL_MAX, so neither evaluation
  // frame can represent the result and the predicate returns a conservative zero;
  // the 1e-200 coordinates flushing under the maximal headroom rescale is a second,
  // independent reason.
  CHECK(sign_incircle(0, 0, 1e200, 0, 0, 1e-200, 1e199, 1e-201) == 0);
  CHECK(sign_incircle(0, 0, 1e200, 0, 0, 1e-200, 2e200, 2e-200) == 0);
}

TEST_CASE("incircle exact path decides near-cocircular points") {
  // One ulp outside and inside the circle through (0,0), (4,0), (0,4), whose
  // centre is (2,2). The differences are at the level of the filter error bound.
  const double outside = std::nextafter(4.0, std::numeric_limits<double>::infinity());
  const double inside = std::nextafter(4.0, 0.0);
  CHECK(sign_incircle(0, 0, 4, 0, 0, 4, 4, outside) < 0);
  CHECK(sign_incircle(0, 0, 4, 0, 0, 4, 4, inside) > 0);

  // Cancellation: the products are about 2^105 and round to the same double,
  // while the exact lifted determinant is negative.
  const double k = 4503599627370496.0;  // 2^52
  const double m = 9007199254740994.0;  // 2^53 + 2
  CHECK(sign_incircle(0, 0, k, m, k + 1.0, m + 2.0, 0.5, 0.5) < 0);
}

// ---------------------------------------------------------------------------
// Cases appended after the fix round 1 review: at 2^-270 the fourth-power
// determinant lands on the smallest subnormal double while the permanent is
// itself subnormal, so kIncircleErrBoundA * permanent underflows to zero. The
// filter must not trust such a determinant; the expected signs are the exact
// integer determinants of the same coordinates divided by 2^-270.
// ---------------------------------------------------------------------------

TEST_CASE("incircle rescales when the filter error bound is subnormal") {
  const double s = std::ldexp(1.0, -270);
  CHECK(sign_incircle(-26 * s, 17 * s, 3 * s, -4 * s, -13 * s, 20 * s, -5 * s, -15 * s) < 0);
  CHECK(sign_incircle(13 * s, -13 * s, 14 * s, -13 * s, 31 * s, -18 * s, 2 * s, -15 * s) > 0);
  CHECK(sign_incircle(-2 * s, -17 * s, 30 * s, -6 * s, 1 * s, -27 * s, 27 * s, -3 * s) < 0);
  CHECK(sign_incircle(36 * s, 3 * s, 38 * s, 5 * s, 9 * s, 0.0, -2 * s, 28 * s) < 0);
  CHECK(sign_incircle(-8 * s, 16 * s, 0.0, -12 * s, -5 * s, 14 * s, -38 * s, -5 * s) > 0);
  // Exactly cocircular at the same scale: the rounded determinant is meaningless,
  // the answer must still be 0.
  CHECK(sign_incircle(0.0, 0.0, s, 0.0, 0.0, s, s, s) == 0);
}

TEST_CASE("incircle at tiny scale matches the exact integer oracle") {
  const double s = std::ldexp(1.0, -270);
  std::mt19937 rng(20261107u);
  std::uniform_int_distribution<int> dist(-40, 40);
  int nonzero = 0;
  for (int i = 0; i < 3000; ++i) {
    const int ax = dist(rng), ay = dist(rng);
    const int bx = dist(rng), by = dist(rng);
    const int cx = dist(rng), cy = dist(rng);
    const int dx = dist(rng), dy = dist(rng);
    const int expected = incircle_oracle(ax, ay, bx, by, cx, cy, dx, dy);
    const int got = sign_incircle(ax * s, ay * s, bx * s, by * s, cx * s, cy * s, dx * s, dy * s);
    REQUIRE(got == expected);
    if (expected != 0) ++nonzero;
  }
  CHECK(nonzero > 1500);
}

// ---------------------------------------------------------------------------
// Exact sign oracle for the randomized ulp suite below.
//
// Every finite double is an integer multiple of 2^-1074, so every coordinate and
// every lifted-determinant term is an integer. The oracle widens each coordinate
// into a fixed-point integer (2^64 limbs, scale 2^-1074) and evaluates the
// determinant with exact integer arithmetic, so no rounding can occur.
namespace {

struct Wide {
  int sign = 0;
  int used = 0;
  std::array<std::uint64_t, 200> limb{};
};

std::uint64_t wide_mul64(std::uint64_t a, std::uint64_t b, std::uint64_t* hi) {
#if defined(_MSC_VER)
  return _umul128(a, b, hi);
#else
  const unsigned __int128 product = static_cast<unsigned __int128>(a) * b;
  *hi = static_cast<std::uint64_t>(product >> 64);
  return static_cast<std::uint64_t>(product);
#endif
}

Wide wide_from_double(double v) {
  Wide r;
  if (v == 0.0) {
    return r;
  }
  std::uint64_t bits = 0;
  std::memcpy(&bits, &v, sizeof bits);
  r.sign = (bits >> 63) != 0 ? -1 : 1;
  const int exponent = static_cast<int>((bits >> 52) & 0x7FFu);
  const std::uint64_t fraction = bits & ((std::uint64_t{1} << 52) - 1);
  const std::uint64_t mantissa = exponent == 0 ? fraction : (fraction | (std::uint64_t{1} << 52));
  const int shift = exponent == 0 ? 0 : exponent - 1;
  const int word = shift / 64;
  const int bit = shift % 64;
  r.limb[static_cast<std::size_t>(word)] |= mantissa << bit;
  if (bit != 0) {
    r.limb[static_cast<std::size_t>(word) + 1] |= mantissa >> (64 - bit);
  }
  r.used = word + (bit != 0 ? 2 : 1);
  while (r.used > 0 && r.limb[static_cast<std::size_t>(r.used - 1)] == 0) {
    --r.used;
  }
  if (r.used == 0) {
    r.sign = 0;
  }
  return r;
}

int wide_compare_magnitude(const Wide& a, const Wide& b) {
  if (a.used != b.used) {
    return a.used < b.used ? -1 : 1;
  }
  for (int i = a.used - 1; i >= 0; --i) {
    const std::uint64_t av = a.limb[static_cast<std::size_t>(i)];
    const std::uint64_t bv = b.limb[static_cast<std::size_t>(i)];
    if (av != bv) {
      return av < bv ? -1 : 1;
    }
  }
  return 0;
}

Wide wide_add_magnitude(const Wide& a, const Wide& b) {
  Wide r;
  const int count = a.used > b.used ? a.used : b.used;
  std::uint64_t carry = 0;
  for (int i = 0; i < count; ++i) {
    const std::uint64_t av = i < a.used ? a.limb[static_cast<std::size_t>(i)] : 0;
    const std::uint64_t bv = i < b.used ? b.limb[static_cast<std::size_t>(i)] : 0;
    const std::uint64_t sum = av + bv;
    const std::uint64_t carry1 = sum < av ? 1 : 0;
    const std::uint64_t sum2 = sum + carry;
    const std::uint64_t carry2 = sum2 < sum ? 1 : 0;
    r.limb[static_cast<std::size_t>(i)] = sum2;
    carry = carry1 + carry2;
  }
  r.used = count;
  if (carry != 0) {
    r.limb[static_cast<std::size_t>(count)] = carry;
    r.used = count + 1;
  }
  return r;
}

// Requires |a| >= |b|.
Wide wide_sub_magnitude(const Wide& a, const Wide& b) {
  Wide r;
  std::uint64_t borrow = 0;
  for (int i = 0; i < a.used; ++i) {
    const std::uint64_t av = a.limb[static_cast<std::size_t>(i)];
    const std::uint64_t bv = i < b.used ? b.limb[static_cast<std::size_t>(i)] : 0;
    const std::uint64_t diff = av - bv;
    const std::uint64_t borrow1 = av < bv ? 1 : 0;
    const std::uint64_t diff2 = diff - borrow;
    const std::uint64_t borrow2 = diff < borrow ? 1 : 0;
    r.limb[static_cast<std::size_t>(i)] = diff2;
    borrow = borrow1 + borrow2;
  }
  r.used = a.used;
  while (r.used > 0 && r.limb[static_cast<std::size_t>(r.used - 1)] == 0) {
    --r.used;
  }
  return r;
}

Wide wide_add(const Wide& a, const Wide& b) {
  if (a.sign == 0) {
    return b;
  }
  if (b.sign == 0) {
    return a;
  }
  if (a.sign == b.sign) {
    Wide r = wide_add_magnitude(a, b);
    r.sign = a.sign;
    return r;
  }
  const int cmp = wide_compare_magnitude(a, b);
  if (cmp == 0) {
    return Wide{};
  }
  Wide r = cmp > 0 ? wide_sub_magnitude(a, b) : wide_sub_magnitude(b, a);
  r.sign = cmp > 0 ? a.sign : b.sign;
  return r;
}

Wide wide_neg(Wide a) {
  a.sign = -a.sign;
  return a;
}

Wide wide_sub(const Wide& a, const Wide& b) {
  return wide_add(a, wide_neg(b));
}

Wide wide_mul(const Wide& a, const Wide& b) {
  Wide r;
  if (a.sign == 0 || b.sign == 0) {
    return r;
  }
  for (int i = 0; i < a.used; ++i) {
    std::uint64_t carry = 0;
    for (int j = 0; j < b.used; ++j) {
      std::uint64_t hi = 0;
      const std::uint64_t lo =
          wide_mul64(a.limb[static_cast<std::size_t>(i)], b.limb[static_cast<std::size_t>(j)], &hi);
      const std::size_t at = static_cast<std::size_t>(i + j);
      const std::uint64_t acc = r.limb[at];
      const std::uint64_t sum = acc + lo;
      const std::uint64_t carry1 = sum < acc ? 1 : 0;
      const std::uint64_t sum2 = sum + carry;
      const std::uint64_t carry2 = sum2 < sum ? 1 : 0;
      r.limb[at] = sum2;
      carry = hi + carry1 + carry2;
    }
    int k = i + b.used;
    while (carry != 0) {
      const std::uint64_t sum = r.limb[static_cast<std::size_t>(k)] + carry;
      carry = sum < carry ? 1 : 0;
      r.limb[static_cast<std::size_t>(k)] = sum;
      ++k;
    }
  }
  r.used = a.used + b.used;
  while (r.used > 0 && r.limb[static_cast<std::size_t>(r.used - 1)] == 0) {
    --r.used;
  }
  r.sign = a.sign * b.sign;
  return r;
}

int wide_bit_length(const Wide& a) {
  if (a.used == 0) {
    return 0;
  }
  std::uint64_t top = a.limb[static_cast<std::size_t>(a.used - 1)];
  int bits = 0;
  while (top != 0) {
    ++bits;
    top >>= 1;
  }
  return (a.used - 1) * 64 + bits;
}

Wide exact_incircle_determinant(const double p[8]) {
  const Wide ax = wide_from_double(p[0]);
  const Wide ay = wide_from_double(p[1]);
  const Wide bx = wide_from_double(p[2]);
  const Wide by = wide_from_double(p[3]);
  const Wide cx = wide_from_double(p[4]);
  const Wide cy = wide_from_double(p[5]);
  const Wide dx = wide_from_double(p[6]);
  const Wide dy = wide_from_double(p[7]);

  const Wide adx = wide_sub(ax, dx);
  const Wide ady = wide_sub(ay, dy);
  const Wide bdx = wide_sub(bx, dx);
  const Wide bdy = wide_sub(by, dy);
  const Wide cdx = wide_sub(cx, dx);
  const Wide cdy = wide_sub(cy, dy);

  const Wide alift = wide_add(wide_mul(adx, adx), wide_mul(ady, ady));
  const Wide blift = wide_add(wide_mul(bdx, bdx), wide_mul(bdy, bdy));
  const Wide clift = wide_add(wide_mul(cdx, cdx), wide_mul(cdy, cdy));

  const Wide cross_bc = wide_sub(wide_mul(bdx, cdy), wide_mul(cdx, bdy));
  const Wide cross_ac = wide_sub(wide_mul(adx, cdy), wide_mul(cdx, ady));
  const Wide cross_ab = wide_sub(wide_mul(adx, bdy), wide_mul(bdx, ady));

  return wide_add(wide_sub(wide_mul(alift, cross_bc), wide_mul(blift, cross_ac)),
                  wide_mul(clift, cross_ab));
}

}  // namespace

TEST_CASE("incircle recovers subnormal perturbations of a circle point") {
  const double sub = std::ldexp(1.0, -1074);  // smallest positive double, 5e-324
  // Circle through (0,0), (4,0), (0,4). The vertex (4,0) is moved into the
  // subnormal range; the exact determinant is about 2^-1068, so the fast filter
  // cannot resolve it and the exact path must run on a scale that keeps the
  // perturbation representable.
  CHECK(sign_incircle(0, 0, 4, 0, 0, 4, 4, sub) > 0);
  CHECK(sign_incircle(0, 0, 4, 0, 0, 4, 4, -sub) < 0);
  CHECK(sign_incircle(0, 0, 4, 0, 0, 4, sub, 4) > 0);
  CHECK(sign_incircle(0, 0, 4, 0, 0, 4, -sub, 4) < 0);
  // A diagonal move by (+sub, -sub) is a tangential perturbation: the exact
  // determinant is about 2^-2143, below the smallest double even after the
  // maximal headroom rescale, so the predicate stays conservative.
  CHECK(sign_incircle(0, 0, 4, 0, 0, 4, sub, -sub) == 0);
  CHECK(sign_incircle(0, 0, 4, 0, 0, 4, 4, 2.0 * sub) > 0);
}

TEST_CASE("incircle falls back to the distributed product for deep cancellation") {
  // Two wild-exponent configurations from the randomized suite. In the scaled
  // frame the grouped form computes cross_ab first and that intermediate cancels
  // below 2^-1074, losing the determinant; distributing the lift into the product
  // (clift*adx)*bdy keeps the same value as three chained exact products and
  // recovers it. The unscaled frame cannot help here: its lift products overflow
  // (the Expansion fills with inf/inf = NaN, and sign() would read NaN as -1).
  // Signs come from an exact rational determinant.
  const double p[8] = {0x1.de8a08ef0603dp+162, 0x1.54caa94adeae7p-181, 0x1.23dc67ad237c2p+108,
                       0x1.cb90cceb66074p-81,  0x1.ae03803c9de2cp+152, 0x1.393c3471a2641p+971,
                       0x1.670d258a5599fp+358,  0x1.d67baed6c1e4ep+155};
  CHECK(sign_incircle(p[0], p[1], p[2], p[3], p[4], p[5], p[6], p[7]) < 0);

  const double q[8] = {0x1.07f677ed8fef2p-373, 0x1.3d42b170871dbp-282, 0x1.f383a3ef05caap-272,
                       0x1.f53195bf09cb0p+608, 0x1.6a2f3e15ad717p-562, 0x1.02bf204c2dc22p-423,
                       0x1.22e43ff675baep-675, 0x1.19d381e7dcbd8p-281};
  CHECK(sign_incircle(q[0], q[1], q[2], q[3], q[4], q[5], q[6], q[7]) > 0);
}

TEST_CASE("incircle ulp-double suite matches the exact integer oracle") {
  std::mt19937 rng(20261107u);
  std::uniform_int_distribution<int> exponent_dist(-1074, 1020);
  std::uniform_real_distribution<double> mantissa_dist(1.0, 2.0);
  const double sub = std::ldexp(1.0, -1074);

  int decided = 0;
  int zeros = 0;
  int zeros_by_flush = 0;
  int zeros_below_representation = 0;
  int zeros_unresolved = 0;
  int wrong_nonzero = 0;
  int oracle_zeros = 0;
  int max_span_bits = 0;
  int max_shift = 0;

  for (int i = 0; i < 3000; ++i) {
    double p[8] = {};
    if (i % 3 == 0) {
      // Wildly mixed magnitudes: each coordinate gets its own exponent.
      for (int k = 0; k < 8; ++k) {
        p[k] = std::ldexp(mantissa_dist(rng), exponent_dist(rng));
      }
    } else if (i % 3 == 1) {
      // The reviewer's shape: a circle point whose second coordinate is subnormal.
      const int k = 2 + (i % 4);
      const double s = std::ldexp(1.0, exponent_dist(rng));
      const double delta = sub * static_cast<double>(1 + (i % 3));
      p[0] = 0.0;
      p[1] = 0.0;
      p[2] = k * s;
      p[3] = 0.0;
      p[4] = 0.0;
      p[5] = k * s;
      p[6] = k * s;
      p[7] = (i % 2 == 0) ? delta : -delta;
    } else {
      // One ulp inside or outside a circle point at a random overall scale.
      const int k = 2 + (i % 4);
      const double s = std::ldexp(1.0, exponent_dist(rng));
      const double v = k * s;
      p[0] = 0.0;
      p[1] = 0.0;
      p[2] = v;
      p[3] = 0.0;
      p[4] = 0.0;
      p[5] = v;
      p[6] = v;
      p[7] = (i % 2 == 0) ? std::nextafter(v, 0.0)
                          : std::nextafter(v, std::numeric_limits<double>::infinity());
    }

    double maximum = 0.0;
    double minimum = std::numeric_limits<double>::infinity();
    for (int k = 0; k < 8; ++k) {
      const double a = std::abs(p[k]);
      maximum = a > maximum ? a : maximum;
      if (a != 0.0 && a < minimum) {
        minimum = a;
      }
    }
    int max_exp = 0;
    int min_exp = 0;
    std::frexp(maximum, &max_exp);
    std::frexp(minimum, &min_exp);
    max_span_bits = std::max(max_span_bits, max_exp - min_exp);

    const Wide det = exact_incircle_determinant(p);
    if (det.sign == 0) {
      ++oracle_zeros;
      continue;
    }
    const int got = sign_incircle(p[0], p[1], p[2], p[3], p[4], p[5], p[6], p[7]);
    if (got == 0) {
      ++zeros;
      int shift = 240 - max_exp;
      shift = shift > 1023 ? 1023 : shift;
      shift = shift < -1074 ? -1074 : shift;
      max_shift = std::max(max_shift, shift);
      const double s = std::ldexp(1.0, shift);
      bool flush = false;
      for (int k = 0; k < 8; ++k) {
        if (p[k] * s == 0.0 && p[k] != 0.0) {
          flush = true;
        }
      }
      if (flush) {
        ++zeros_by_flush;
      } else if (wide_bit_length(det) - 1 < 3222 - 4 * shift + 32) {
        ++zeros_below_representation;
      } else {
        // The determinant is representable and no coordinate flushed, yet both
        // the scaled and the unscaled exact frame returned zero. The dual-frame
        // rule exists to prevent this, so any case here is a real gap.
        ++zeros_unresolved;
      }
    } else if (got != det.sign) {
      ++wrong_nonzero;
    } else {
      ++decided;
    }
  }

  INFO("decided=" << decided << " zeros=" << zeros << " flush=" << zeros_by_flush
                  << " below=" << zeros_below_representation << " unresolved=" << zeros_unresolved
                  << " wrong=" << wrong_nonzero << " oracle_zeros=" << oracle_zeros
                  << " max_span_bits=" << max_span_bits << " max_shift=" << max_shift);
  CHECK(wrong_nonzero == 0);
  CHECK(zeros_unresolved == 0);
  CHECK(decided + zeros + oracle_zeros == 3000);
  // Every case that survives the rescale without flushing a coordinate must be
  // decided: the dual-frame rule leaves only flush zeros. In this suite that is
  // 1068 of the 3000 cases (span up to 2084 bits), so 1932 are decided.
  CHECK(zeros == zeros_by_flush + zeros_below_representation);
}
