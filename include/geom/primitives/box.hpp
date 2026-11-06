#pragma once

#include <cmath>
#include <limits>
#include <type_traits>

#include <geom/primitives/point.hpp>

namespace geom {
inline namespace v1 {

template <class T, int D>
struct Box {
  Point<T, D> lo;
  Point<T, D> hi;

  static constexpr Box empty() noexcept {
    Box b;
    for (int i = 0; i < D; ++i) {
      b.lo[i] = std::numeric_limits<T>::infinity();
      b.hi[i] = -std::numeric_limits<T>::infinity();
    }
    return b;
  }

  constexpr bool is_empty() const noexcept {
    for (int i = 0; i < D; ++i) {
      if (lo[i] > hi[i]) return true;
    }
    return false;
  }

  constexpr void expand(Point<T, D> p) noexcept {
    for (int i = 0; i < D; ++i) {
      if (p[i] < lo[i]) lo[i] = p[i];
      if (p[i] > hi[i]) hi[i] = p[i];
    }
  }

  constexpr void expand(Box other) noexcept {
    for (int i = 0; i < D; ++i) {
      if (other.lo[i] < lo[i]) lo[i] = other.lo[i];
      if (other.hi[i] > hi[i]) hi[i] = other.hi[i];
    }
  }

  // A point is contained only if it is finite and every coordinate compares
  // inside the box in both directions. NaN fails the comparisons; an explicit
  // finite check also rejects +/-inf, which would otherwise be "inside" a box
  // with infinite bounds.
  constexpr bool contains(Point<T, D> p) const noexcept {
    for (int i = 0; i < D; ++i) {
      if constexpr (std::is_floating_point_v<T>) {
        if (!std::isfinite(p[i])) return false;
      }
      if (!(p[i] >= lo[i] && p[i] <= hi[i])) return false;
    }
    return true;
  }

  constexpr Box merged(Box other) const noexcept {
    Box r = *this;
    r.expand(other);
    return r;
  }

  constexpr Vec<T, D> diagonal() const noexcept {
    Vec<T, D> d;
    for (int i = 0; i < D; ++i) d[i] = hi[i] - lo[i];
    return d;
  }
};

using Box2d = Box<double, 2>;
using Box3d = Box<double, 3>;

}  // namespace v1
}  // namespace geom
