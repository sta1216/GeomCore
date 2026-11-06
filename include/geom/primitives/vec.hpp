#pragma once

#include <cmath>

namespace geom {
inline namespace v1 {

template <class T, int D>
struct Vec {
  T v[D]{};

  constexpr T& operator[](int i) noexcept { return v[i]; }
  constexpr const T& operator[](int i) const noexcept { return v[i]; }
};

template <class T, int D>
constexpr Vec<T, D> operator+(Vec<T, D> a, Vec<T, D> b) noexcept {
  Vec<T, D> r;
  for (int i = 0; i < D; ++i) r[i] = a[i] + b[i];
  return r;
}

template <class T, int D>
constexpr Vec<T, D> operator-(Vec<T, D> a, Vec<T, D> b) noexcept {
  Vec<T, D> r;
  for (int i = 0; i < D; ++i) r[i] = a[i] - b[i];
  return r;
}

template <class T, int D>
constexpr Vec<T, D> operator*(Vec<T, D> a, T s) noexcept {
  Vec<T, D> r;
  for (int i = 0; i < D; ++i) r[i] = a[i] * s;
  return r;
}

template <class T, int D>
constexpr T dot(Vec<T, D> a, Vec<T, D> b) noexcept {
  T s = T{0};
  for (int i = 0; i < D; ++i) s += a[i] * b[i];
  return s;
}

template <class T>
constexpr Vec<T, 3> cross(Vec<T, 3> a, Vec<T, 3> b) noexcept {
  return Vec<T, 3>{a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2],
                   a[0] * b[1] - a[1] * b[0]};
}

template <class T, int D>
inline T norm2(Vec<T, D> a) noexcept {
  return std::sqrt(dot(a, a));
}

using Vec2d = Vec<double, 2>;
using Vec3d = Vec<double, 3>;

}  // namespace v1
}  // namespace geom
