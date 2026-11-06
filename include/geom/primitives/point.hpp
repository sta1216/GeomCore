#pragma once

#include <geom/primitives/vec.hpp>

namespace geom {
inline namespace v1 {

template <class T, int D>
struct Point {
  T v[D]{};

  constexpr T& operator[](int i) noexcept { return v[i]; }
  constexpr const T& operator[](int i) const noexcept { return v[i]; }
};

template <class T, int D>
constexpr Point<T, D> operator+(Point<T, D> p, Vec<T, D> v) noexcept {
  Point<T, D> r;
  for (int i = 0; i < D; ++i) r[i] = p[i] + v[i];
  return r;
}

template <class T, int D>
constexpr Point<T, D> operator-(Point<T, D> p, Vec<T, D> v) noexcept {
  Point<T, D> r;
  for (int i = 0; i < D; ++i) r[i] = p[i] - v[i];
  return r;
}

template <class T, int D>
constexpr Vec<T, D> operator-(Point<T, D> a, Point<T, D> b) noexcept {
  Vec<T, D> r;
  for (int i = 0; i < D; ++i) r[i] = a[i] - b[i];
  return r;
}

using Point2d = Point<double, 2>;
using Point3d = Point<double, 3>;

}  // namespace v1
}  // namespace geom
