#include <doctest/doctest.h>

#include <cmath>

#include <geom/numeric/expansion.hpp>

using geom::v1::detail::Expansion;
using geom::v1::detail::two_product;
using geom::v1::detail::two_sum;

TEST_CASE("two_sum keeps the residual of 2^53 + 1") {
  double s = 0.0, e = 0.0;
  two_sum(9007199254740992.0, 1.0, s, e);
  CHECK(s == 9007199254740992.0);
  CHECK(e == 1.0);
}

TEST_CASE("two_product keeps the residual of (2^27+1)*(2^27-1)") {
  const double a = 134217729.0;  // 2^27 + 1
  const double b = 134217727.0;  // 2^27 - 1
  double p = 0.0, e = 0.0;
  two_product(a, b, p, e);
  CHECK(p == 18014398509481984.0);  // 2^54
  CHECK(e == -1.0);                 // exact product is 2^54 - 1
}

TEST_CASE("Expansion adds small values to large ones without loss") {
  Expansion e(1.0e16);
  e.add_term(1.0);
  REQUIRE(e.size() >= 2);
  CHECK(e.approx() == 1.0e16);
  CHECK(e[0] == 1.0);
  CHECK(e[e.size() - 1] == 1.0e16);
  CHECK(e.sign() > 0);
}

TEST_CASE("Expansion subtraction is exact") {
  const Expansion d = Expansion::difference(1.0, 1.0e16);
  CHECK(d.size() >= 2);
  CHECK(d.approx() == -1.0e16);
  CHECK(d.sign() < 0);
}

TEST_CASE("Expansion multiplication keeps the exact product") {
  const Expansion a = Expansion::difference(1.0e8, -1.0);  // 1e8 + 1
  const Expansion b = Expansion::difference(1.0e8, 1.0);   // 1e8 - 1
  const Expansion p = Expansion::mul(a, b);
  CHECK(p.size() >= 2);
  CHECK(p.approx() == 1.0e16);  // (1e16 - 1) rounds to 1e16 in double
  CHECK(p.sign() > 0);
  const Expansion zero = Expansion::difference(p.approx(), p.approx());
  CHECK(zero.sign() == 0);
}

TEST_CASE("Expansion sign scans the most significant component") {
  Expansion e = Expansion::difference(1.0e16, 1.0);
  CHECK(e.sign() > 0);
  e = Expansion::difference(1.0, 1.0e16);
  CHECK(e.sign() < 0);
  CHECK(Expansion(0.0).sign() == 0);
}
