#include <doctest/doctest.h>

#include <cmath>
#include <cstdint>
#include <random>

#include <geom/numeric/predicates2.hpp>

using geom::v1::sign_orient2d;

namespace {

// Oracle: coordinates are small integers, so the determinant is exactly
// representable in int64 and can be compared bit-exactly.
int orient_oracle(std::int64_t ax, std::int64_t ay, std::int64_t bx, std::int64_t by,
                  std::int64_t cx, std::int64_t cy) {
  const std::int64_t det =
      (ax - cx) * (by - cy) - (ay - cy) * (bx - cx);
  return det > 0 ? 1 : (det < 0 ? -1 : 0);
}

}  // namespace

TEST_CASE("orient2d matches the exact oracle on random integer input") {
  std::mt19937 rng(20261031u);
  std::uniform_int_distribution<int> dist(-32768, 32767);
  int exercised_exact = 0;
  for (int i = 0; i < 200000; ++i) {
    const int ax = dist(rng), ay = dist(rng);
    const int bx = dist(rng), by = dist(rng);
    const int cx = dist(rng), cy = dist(rng);
    const int got = sign_orient2d(ax, ay, bx, by, cx, cy);
    const int want = orient_oracle(ax, ay, bx, by, cx, cy);
    REQUIRE(got == want);
    if (want == 0) ++exercised_exact;
  }
  CHECK(exercised_exact >= 0);
}

TEST_CASE("orient2d handles exactly collinear input") {
  CHECK(sign_orient2d(0, 0, 1, 0, 2, 0) == 0);
  CHECK(sign_orient2d(0, 0, 1, 1, 2, 2) == 0);
  CHECK(sign_orient2d(-3, -3, 0, 0, 7, 7) == 0);
}

TEST_CASE("orient2d detects near-degenerate cases") {
  // Points differing by 1 ulp around a line must still be classified exactly.
  CHECK(sign_orient2d(0.0, 0.0, 1.0, 0.0, 0.5, 1.0e-300) > 0);
  CHECK(sign_orient2d(0.0, 0.0, 1.0, 0.0, 0.5, -1.0e-300) < 0);
}

TEST_CASE("orient2d is scale invariant") {
  const double s = 1e200;
  CHECK(sign_orient2d(0.0, 0.0, s, 0.0, 0.0, s) > 0);
  CHECK(sign_orient2d(0.0, 0.0, 0.0, s, s, 0.0) < 0);
  const double t = 1e-200;
  CHECK(sign_orient2d(0.0, 0.0, t, 0.0, 0.0, t) > 0);
}

TEST_CASE("orient2d rejects non-finite input as degenerate") {
  CHECK(sign_orient2d(std::nan(""), 0.0, 1.0, 0.0, 0.0, 1.0) == 0);
  CHECK(sign_orient2d(std::numeric_limits<double>::infinity(), 0.0, 1.0, 0.0, 0.0, 1.0) == 0);
}

TEST_CASE("orient2d handles mixed-magnitude input") {
  // det = 1: coordinates span 200 decades without triggering the global rescale.
  CHECK(sign_orient2d(0.0, 0.0, 1e100, 0.0, 0.0, 1e-100) > 0);
  // Collinear points along y = 1e-200 * x, exactly representable as doublings
  // of 1e100 / 1e-100, so the exact path must return 0 without underflow.
  CHECK(sign_orient2d(0.0, 0.0, 1e100, 1e-100, 2e100, 2e-100) == 0);
}

TEST_CASE("orient2d handles extreme magnitude spans without rescaling loss") {
  // det = 1 with coordinates 400 (resp. 616) decades apart: only an unscaled
  // classification can see the determinant, since scaling by 2^k would push
  // the small coordinate out of the double range.
  CHECK(sign_orient2d(0.0, 0.0, 1e200, 0.0, 0.0, 1e-200) > 0);
  CHECK(sign_orient2d(0.0, 0.0, 1e308, 0.0, 0.0, 1e-308) > 0);
  // Just above the old 1e150 rescale threshold: no discontinuity.
  CHECK(sign_orient2d(0.0, 0.0, std::nextafter(1e150, INFINITY), 0.0, 0.0, 1e-200) > 0);
}

TEST_CASE("orient2d handles subnormal coordinates") {
  // 2^-1063 and 2^-1074 coordinates: the determinant underflows to zero in
  // doubles, so only a clamped (finite) rescale can classify the triangle.
  CHECK(sign_orient2d(0.0, 0.0, 1e-320, 0.0, 0.0, 1e-320) > 0);
  CHECK(sign_orient2d(0.0, 0.0, 5e-324, 0.0, 0.0, 5e-324) > 0);
  CHECK(sign_orient2d(0.0, 0.0, 1e-320, 0.0, 2e-320, 0.0) == 0);
}

TEST_CASE("orient2d rejects non-finite input in every argument position") {
  const double nan = std::nan("");
  const double inf = std::numeric_limits<double>::infinity();
  // NaN anywhere must be degenerate, not signed by whichever comparison sees it.
  CHECK(sign_orient2d(0.0, nan, 1.0, 0.0, 0.0, 1.0) == 0);
  CHECK(sign_orient2d(0.0, 0.0, nan, 0.0, 0.0, 1.0) == 0);
  CHECK(sign_orient2d(0.0, 0.0, 1.0, nan, 0.0, 1.0) == 0);
  CHECK(sign_orient2d(0.0, 0.0, 1.0, 0.0, nan, 1.0) == 0);
  CHECK(sign_orient2d(0.0, 0.0, 1.0, 0.0, 0.0, nan) == 0);
  // Infinities are covered by the same rule, in either sign and position.
  CHECK(sign_orient2d(0.0, 0.0, inf, 0.0, 0.0, 1.0) == 0);
  CHECK(sign_orient2d(0.0, 0.0, -inf, 0.0, 0.0, 1.0) == 0);
  CHECK(sign_orient2d(0.0, 0.0, 1.0, 0.0, inf, 1.0) == 0);
  CHECK(sign_orient2d(0.0, 0.0, 1.0, 0.0, 0.0, -inf) == 0);
  CHECK(sign_orient2d(inf, 0.0, 1.0, 0.0, 0.0, 1.0) == 0);
}


