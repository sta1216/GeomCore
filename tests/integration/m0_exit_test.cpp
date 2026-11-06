#include <doctest/doctest.h>

#include <geom/numeric/predicates2.hpp>

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
