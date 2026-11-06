#pragma once

#include <utility>
#include <vector>

#include <geom/primitives/point.hpp>
#include <geom/primitives/vec.hpp>

namespace geom {
inline namespace v1 {

template <class T, int D>
class BezierCurve {
 public:
  explicit BezierCurve(std::vector<Point<T, D>> control) : c_(std::move(control)) {}

  std::size_t size() const noexcept { return c_.size(); }

  Point<T, D> eval(T t) const noexcept {
    if (c_.empty()) return Point<T, D>{};
    std::vector<Point<T, D>> d = c_;
    for (std::size_t k = 1; k < d.size(); ++k) {
      for (std::size_t i = 0; i + k < d.size(); ++i) {
        Point<T, D> q;
        for (int c = 0; c < D; ++c) q[c] = (T{1} - t) * d[i][c] + t * d[i + 1][c];
        d[i] = q;
      }
    }
    return d[0];
  }

  Vec<T, D> derivative(T t) const noexcept {
    if (c_.size() < 2) return Vec<T, D>{};
    const T p = static_cast<T>(c_.size() - 1);
    std::vector<Vec<T, D>> d(c_.size() - 1);
    for (std::size_t i = 0; i + 1 < c_.size(); ++i) d[i] = (c_[i + 1] - c_[i]) * p;
    for (std::size_t k = 1; k < d.size(); ++k) {
      for (std::size_t i = 0; i + k < d.size(); ++i) {
        Vec<T, D> q;
        for (int c = 0; c < D; ++c) q[c] = (T{1} - t) * d[i][c] + t * d[i + 1][c];
        d[i] = q;
      }
    }
    return d[0];
  }

 private:
  std::vector<Point<T, D>> c_;
};

}  // namespace v1
}  // namespace geom
