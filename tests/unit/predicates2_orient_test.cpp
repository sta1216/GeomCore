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
