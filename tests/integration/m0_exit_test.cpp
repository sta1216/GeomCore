#include <doctest/doctest.h>

#include <geom/numeric/predicates2.hpp>

#include <cmath>
#include <cstdint>
#include <random>

TEST_CASE("M0 exit: orient2d is consistent on a triangulated grid") {
  int positives = 0, negatives = 0;
  for (int y = 0; y < 8; ++y) {
    for (int x = 0; x < 8; ++x) {
      const int s = geom::v1::sign_orient2d(0, 0, 1, 0, x, y);
      if (s > 0) ++positives;
      if (s < 0) ++negatives;
    }
  }
  // y == 0 is collinear with the base edge; y > 0 yields 7 * 8 positive cases.
  CHECK(positives == 56);
  CHECK(negatives == 0);
}

namespace {

// Exact reference for integer coordinates: the orientation determinant of
// (ax, ay), (bx, by), (cx, cy) is an int64 value in this test's range, so its
// sign is computed without any floating point at all.
std::int64_t det64(std::int64_t ax, std::int64_t ay, std::int64_t bx, std::int64_t by,
                   std::int64_t cx, std::int64_t cy) {
  return (bx - ax) * (cy - ay) - (by - ay) * (cx - ax);
}

int sign64(std::int64_t det) { return det > 0 ? 1 : (det < 0 ? -1 : 0); }

std::int64_t random_coord(std::mt19937& rng) {
  return static_cast<std::int64_t>(rng() % 8193u) - 4096;
}

}  // namespace

TEST_CASE("M0 exit: every grid point reports its own orientation sign") {
  // Triangle (0,0), (1,0), (x,y): positive above the base edge, exactly zero on
  // it, negative below it. Asserting the sign per point (rather than totals)
  // makes a predicate that returns a constant sign fail immediately.
  for (int y = -4; y < 8; ++y) {
    for (int x = -8; x < 8; ++x) {
      const int expected = (y > 0) - (y < 0);
      const int actual = geom::v1::sign_orient2d(0, 0, 1, 0, x, y);
      CHECK_MESSAGE(actual == expected, "grid point (", x, ", ", y, ")");
    }
  }
}

TEST_CASE("M0 exit: pseudorandom integer coordinates match the exact int64 determinant") {
  // One million cases with a fixed seed; every 997th case is collinear by
  // construction so that the degenerate sign is exercised as well.
  std::mt19937 rng(20261107u);
  constexpr int kCases = 1000000;
  int mismatches = 0;
  int first_mismatch = -1;
  int zeros = 0;
  int negatives = 0;
  for (int i = 0; i < kCases; ++i) {
    const std::int64_t ax = random_coord(rng);
    const std::int64_t ay = random_coord(rng);
    std::int64_t bx = 0, by = 0, cx = 0, cy = 0;
    if (i % 997 == 0) {
      std::int64_t dx = random_coord(rng);
      const std::int64_t dy = random_coord(rng);
      if (dx == 0 && dy == 0) dx = 1;
      bx = ax + dx;
      by = ay + dy;
      cx = ax + 2 * dx;
      cy = ay + 2 * dy;
    } else {
      bx = random_coord(rng);
      by = random_coord(rng);
      cx = random_coord(rng);
      cy = random_coord(rng);
    }
    const int expected = sign64(det64(ax, ay, bx, by, cx, cy));
    const int actual =
        geom::v1::sign_orient2d(static_cast<double>(ax), static_cast<double>(ay),
                                static_cast<double>(bx), static_cast<double>(by),
                                static_cast<double>(cx), static_cast<double>(cy));
    if (actual != expected) {
      if (first_mismatch < 0) first_mismatch = i;
      ++mismatches;
    }
    if (expected == 0) ++zeros;
    if (expected < 0) ++negatives;
  }
  CHECK_MESSAGE(mismatches == 0, "first mismatch at case ", first_mismatch);
  CHECK(zeros > 0);
  CHECK(negatives > 0);
}

TEST_CASE("M0 exit: degenerate, near-degenerate and wide-span cases keep their sign") {
  using geom::v1::sign_orient2d;

  // Exactly degenerate configurations.
  CHECK(sign_orient2d(0.0, 0.0, 4.0, 4.0, 2.0, 2.0) == 0);
  CHECK(sign_orient2d(-3.0, -3.0, 5.0, 5.0, 1.0, 1.0) == 0);
  CHECK(sign_orient2d(0.0, 0.0, 0.0, 0.0, 7.0, 9.0) == 0);
  CHECK(sign_orient2d(0.0, 0.0, 4.0, 4.0, 4.0, 4.0) == 0);

  // One ulp off the line y = x: the determinant is 2^-49 and -2^-50, far below
  // what the fast filter can resolve, so the exact path decides.
  const double above = std::nextafter(2.0, 4.0);
  const double below = std::nextafter(2.0, 0.0);
  CHECK(sign_orient2d(0.0, 0.0, 4.0, 4.0, 2.0, above) > 0);
  CHECK(sign_orient2d(0.0, 0.0, 4.0, 4.0, 2.0, below) < 0);

  // Cancellation: the products are about 2^105 and round to the same double,
  // while the exact determinant is -2. A naive double determinant reports 0.
  const double k = 4503599627370496.0;  // 2^52
  const double m = 9007199254740994.0;  // 2^53 + 2
  CHECK(sign_orient2d(0.0, 0.0, k, m, k + 1.0, m + 2.0) < 0);

  // Wide span: the determinant is 1 although the coordinates differ by 200
  // decades, and exact collinearity is still reported as 0.
  CHECK(sign_orient2d(0.0, 0.0, 1e100, 0.0, 0.0, 1e-100) > 0);
  CHECK(sign_orient2d(0.0, 0.0, 1e200, 0.0, 0.0, 1e-200) > 0);
  CHECK(sign_orient2d(0.0, 0.0, 1e308, 0.0, 0.0, 1e-308) > 0);
  CHECK(sign_orient2d(0.0, 0.0, 1e200, 1e200, 1e-200, 1e-200) == 0);

  // Subnormal coordinates: both unscaled products underflow to zero, so the
  // scaled path has to recover the sign; exact collinearity stays 0.
  CHECK(sign_orient2d(0.0, 0.0, 1e-320, 0.0, 0.0, 1e-320) > 0);
  CHECK(sign_orient2d(0.0, 0.0, 0.0, 1e-320, 1e-320, 0.0) < 0);
  CHECK(sign_orient2d(0.0, 0.0, 1e-320, 0.0, 2e-320, 0.0) == 0);
}
