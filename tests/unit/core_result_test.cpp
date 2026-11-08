#include <doctest/doctest.h>

#include <stdexcept>
#include <string_view>
#include <utility>

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
  CHECK(std::string_view(geom::v1::to_string(Errc::none)) == "none");
  CHECK(std::string_view(geom::v1::to_string(Errc::non_finite)) == "non_finite");
  CHECK(std::string_view(geom::v1::to_string(Errc::index_out_of_range)) ==
        "index_out_of_range");
  CHECK(std::string_view(geom::v1::to_string(Errc::io_error)) == "io_error");
  CHECK(std::string_view(geom::v1::to_string(Errc::unsupported)) == "unsupported");
  CHECK(std::string_view(geom::v1::to_string(Errc::cancelled)) == "cancelled");
  CHECK(std::string_view(geom::v1::to_string(Errc::self_intersecting)) ==
        "self_intersecting");
  CHECK(std::string_view(geom::v1::to_string(static_cast<Errc>(0x7f))) == "unknown");
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

TEST_CASE("Result value access throws on an error state in every ref category") {
  Result<int> bad{Error{Errc::degenerate, "degenerate triangle", 0u, 0u}};
  CHECK_THROWS_AS(bad.value(), std::runtime_error);
  const Result<int>& const_bad = bad;
  CHECK_THROWS_AS(const_bad.value(), std::runtime_error);
  CHECK_THROWS_AS(std::move(bad).value(), std::runtime_error);

  Result<int> no_message{Error{Errc::numeric_failure, nullptr, 0u, 0u}};
  CHECK_THROWS_WITH_AS(no_message.value(), "geometry error", std::runtime_error);

  Result<int> ok{7};
  CHECK(static_cast<bool>(ok));
  CHECK_FALSE(static_cast<bool>(no_message));
  const Result<int>& const_ok = ok;
  CHECK(const_ok.value() == 7);
}

TEST_CASE("Result<void> models success and failure") {
  Result<void> ok{};
  CHECK(ok.has_value());

  Result<void> bad{Error{Errc::cancelled, "cancelled", 0u, 0u}};
  REQUIRE_FALSE(bad.has_value());
  CHECK(bad.error().code == Errc::cancelled);
}
