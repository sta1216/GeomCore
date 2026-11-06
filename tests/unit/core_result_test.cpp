#include <doctest/doctest.h>

#include <string_view>

#include <geom/core/config.hpp>
#include <geom/core/error.hpp>
#include <geom/core/result.hpp>

using geom::v1::Errc;
using geom::v1::Error;
using geom::v1::Result;

TEST_CASE("version constants exist") {
  static_assert(geom::v1::version_major == 0);
  CHECK(std::string_view(geom::v1::version_string) == "0.1.0");
}

TEST_CASE("Errc has stable string names") {
  CHECK(std::string_view(geom::v1::to_string(Errc::degenerate)) == "degenerate");
  CHECK(std::string_view(geom::v1::to_string(Errc::numeric_failure)) == "numeric_failure");
}

TEST_CASE("Result carries either a value or an error") {
  Result<int> ok{42};
  REQUIRE(ok.has_value());
  CHECK(ok.value() == 42);

  Result<int> bad{Error{Errc::degenerate, "degenerate triangle", 7u, 0u}};
  REQUIRE_FALSE(bad.has_value());
  CHECK(bad.error().code == Errc::degenerate);
  CHECK(bad.error().context_a == 7u);
}

TEST_CASE("Result<void> models success and failure") {
  Result<void> ok{};
  CHECK(ok.has_value());

  Result<void> bad{Error{Errc::cancelled, "cancelled", 0u, 0u}};
  REQUIRE_FALSE(bad.has_value());
  CHECK(bad.error().code == Errc::cancelled);
}
