#pragma once

#include <geom/primitives/point.hpp>

namespace geom {
inline namespace v1 {

template <class T, int D>
struct Segment {
  Point<T, D> a;
  Point<T, D> b;

  constexpr Point<T, D> closest_point(Point<T, D> p) const noexcept {
    const Vec<T, D> ab = b - a;
    const T denom = dot(ab, ab);
    if (denom == T{0}) return a;
    T t = dot(p - a, ab) / denom;
    if (t < T{0}) t = T{0};
    if (t > T{1}) t = T{1};
    return a + ab * t;
  }
};

template <class T, int D>
constexpr T squared_distance(Point<T, D> p, Segment<T, D> s) noexcept {
  const Vec<T, D> d = p - s.closest_point(p);
  return dot(d, d);
}

using Segment2d = Segment<double, 2>;
using Segment3d = Segment<double, 3>;

}  // namespace v1
}  // namespace geom
