#pragma once

#include <cstdlib>
#include <optional>
#include <stdexcept>
#include <utility>
#include <variant>

#include <geom/core/error.hpp>

namespace geom {
inline namespace v1 {

namespace detail {

[[noreturn]] inline void result_contract_violation(const Error& e) {
#if defined(__cpp_exceptions)
  throw std::runtime_error(e.message != nullptr ? e.message : "geometry error");
#else
  (void)e;
  std::abort();
#endif
}

}  // namespace detail

template <class T>
class Result {
 public:
  Result(T value) : v_(std::in_place_index<1>, std::move(value)) {}
  Result(Error error) : v_(std::in_place_index<0>, error) {}

  bool has_value() const noexcept { return v_.index() == 1; }
  explicit operator bool() const noexcept { return has_value(); }

  T& value() & {
    if (!has_value()) detail::result_contract_violation(std::get<0>(v_));
    return std::get<1>(v_);
  }
  const T& value() const& {
    if (!has_value()) detail::result_contract_violation(std::get<0>(v_));
    return std::get<1>(v_);
  }
  T&& value() && {
    if (!has_value()) detail::result_contract_violation(std::get<0>(v_));
    return std::move(std::get<1>(v_));
  }

  Error& error() & { return std::get<0>(v_); }
  const Error& error() const& { return std::get<0>(v_); }

 private:
  std::variant<Error, T> v_;
};

template <>
class Result<void> {
 public:
  Result() = default;
  Result(Error error) : e_(error) {}

  bool has_value() const noexcept { return !e_.has_value(); }
  explicit operator bool() const noexcept { return has_value(); }

  Error& error() & { return *e_; }
  const Error& error() const& { return *e_; }

 private:
  std::optional<Error> e_;
};

}  // namespace v1
}  // namespace geom
