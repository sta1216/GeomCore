#pragma once

#include <cmath>
#include <limits>
#include <type_traits>

namespace geom {
inline namespace v1 {

template <class T>
struct scalar_traits {
  static_assert(std::is_floating_point_v<T>, "only IEEE floating point scalars are supported");
  using value_type = T;
  static constexpr T epsilon() noexcept { return std::numeric_limits<T>::epsilon(); }
  static constexpr T infinity() noexcept { return std::numeric_limits<T>::infinity(); }
  static bool is_finite(T v) noexcept { return std::isfinite(v); }
};

}  // namespace v1
}  // namespace geom
