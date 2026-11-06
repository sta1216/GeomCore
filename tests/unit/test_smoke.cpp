#include <doctest/doctest.h>

TEST_CASE("toolchain is C++20 or newer") {
  static_assert(__cplusplus >= 202002L);
  CHECK(sizeof(double) == 8);
}
