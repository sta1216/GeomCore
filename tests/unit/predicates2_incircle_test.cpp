#include <doctest/doctest.h>

#include <cstdint>
#include <limits>
#include <random>

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

  // 400 decades do not survive a single global rescale: the 1e-200 coordinates
  // flush to zero. The contract prefers a conservative 0 over a guessed sign.
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
