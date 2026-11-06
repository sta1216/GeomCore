#pragma once

#include <array>
#include <cstddef>
#include <utility>
#include <vector>

#include <geom/primitives/point.hpp>
#include <geom/primitives/vec.hpp>

namespace geom {
inline namespace v1 {
namespace detail {

constexpr int kMaxCurveDegree = 16;

inline int find_span(const std::vector<double>& U, int p, int n, double u) noexcept {
  if (u >= U[static_cast<std::size_t>(n)]) return n;
  int low = p;
  int high = n + 1;
  int mid = (low + high) / 2;
  while (u < U[static_cast<std::size_t>(mid)] || u >= U[static_cast<std::size_t>(mid + 1)]) {
    if (u < U[static_cast<std::size_t>(mid)]) {
      high = mid;
    } else {
      low = mid;
    }
    mid = (low + high) / 2;
  }
  return mid;
}

// de Boor for a 4D homogeneous B-spline. pts.size() == U.size() - p - 1.
// Requires p <= kMaxCurveDegree and U.size() == pts.size() + p + 1.
inline std::array<double, 4> eval_curve(const std::vector<double>& U, int p,
                                        const std::vector<std::array<double, 4>>& pts,
                                        double u) noexcept {
  const int n = static_cast<int>(pts.size()) - 1;
  const int k = find_span(U, p, n, u);
  std::array<std::array<double, 4>, kMaxCurveDegree + 1> d{};
  for (int j = 0; j <= p; ++j) {
    d[static_cast<std::size_t>(j)] = pts[static_cast<std::size_t>(k - p + j)];
  }
  for (int r = 1; r <= p; ++r) {
    for (int j = p; j >= r; --j) {
      const int i0 = k - p + j;
      const double denom =
          U[static_cast<std::size_t>(k + j + 1 - r)] - U[static_cast<std::size_t>(i0)];
      const double a = denom == 0.0 ? 0.0 : (u - U[static_cast<std::size_t>(i0)]) / denom;
      for (int c = 0; c < 4; ++c) {
        d[static_cast<std::size_t>(j)][static_cast<std::size_t>(c)] =
            (1.0 - a) * d[static_cast<std::size_t>(j - 1)][static_cast<std::size_t>(c)] +
            a * d[static_cast<std::size_t>(j)][static_cast<std::size_t>(c)];
      }
    }
  }
  return d[static_cast<std::size_t>(p)];
}

}  // namespace detail

class NurbsCurve {
 public:
  NurbsCurve(int degree, std::vector<double> knots, std::vector<Point3d> control,
             std::vector<double> weights)
      : p_(degree), U_(std::move(knots)), P_(std::move(control)), W_(std::move(weights)) {}

  bool is_valid() const noexcept {
    if (p_ < 1 || P_.empty()) return false;
    if (p_ > detail::kMaxCurveDegree) return false;
    if (P_.size() <= static_cast<std::size_t>(p_)) return false;
    if (W_.size() != P_.size()) return false;
    if (U_.size() != P_.size() + static_cast<std::size_t>(p_) + 1) return false;
    for (std::size_t i = 1; i < U_.size(); ++i) {
      if (!(U_[i] >= U_[i - 1])) return false;
    }
    for (double w : W_) {
      if (!(w > 0.0)) return false;
    }
    return true;
  }

  int degree() const noexcept { return p_; }
  std::size_t control_count() const noexcept { return P_.size(); }

  std::pair<double, double> domain() const noexcept {
    if (!is_valid()) return {0.0, 0.0};
    return {U_[static_cast<std::size_t>(p_)],
            U_[U_.size() - static_cast<std::size_t>(p_) - 1]};
  }

  // Both queries require is_valid(); the parameter is clamped to domain().
  Point3d eval(double u) const noexcept {
    const auto h = eval_h(clamp(u));
    return Point3d{h[0] / h[3], h[1] / h[3], h[2] / h[3]};
  }

  Vec3d derivative(double u) const noexcept {
    const auto h = eval_h(clamp(u));
    const auto dh = deriv_h(clamp(u));
    const double w = h[3];
    const Vec3d c{h[0] / w, h[1] / w, h[2] / w};
    return Vec3d{(dh[0] - c[0] * dh[3]) / w, (dh[1] - c[1] * dh[3]) / w,
                 (dh[2] - c[2] * dh[3]) / w};
  }

 private:
  double clamp(double u) const noexcept {
    const auto [a, b] = domain();
    return u < a ? a : (u > b ? b : u);
  }

  std::vector<std::array<double, 4>> homogeneous() const {
    std::vector<std::array<double, 4>> H(P_.size());
    for (std::size_t i = 0; i < P_.size(); ++i) {
      H[i] = {W_[i] * P_[i][0], W_[i] * P_[i][1], W_[i] * P_[i][2], W_[i]};
    }
    return H;
  }

  std::array<double, 4> eval_h(double u) const noexcept {
    const auto H = homogeneous();
    return detail::eval_curve(U_, p_, H, u);
  }

  // Derivative of the homogeneous curve: degree p-1 control points
  // p (H[i+1] - H[i]) / (U[i+p+1] - U[i+1]) with the first and last knot removed.
  std::array<double, 4> deriv_h(double u) const noexcept {
    const auto H = homogeneous();
    std::vector<double> U1(U_.begin() + 1, U_.end() - 1);
    std::vector<std::array<double, 4>> D(H.size() - 1);
    for (std::size_t i = 0; i + 1 < H.size(); ++i) {
      const double denom = U_[i + static_cast<std::size_t>(p_) + 1] - U_[i + 1];
      const double f = denom == 0.0 ? 0.0 : static_cast<double>(p_) / denom;
      for (int c = 0; c < 4; ++c) {
        D[i][static_cast<std::size_t>(c)] =
            f * (H[i + 1][static_cast<std::size_t>(c)] - H[i][static_cast<std::size_t>(c)]);
      }
    }
    return detail::eval_curve(U1, p_ - 1, D, u);
  }

  int p_ = 0;
  std::vector<double> U_;
  std::vector<Point3d> P_;
  std::vector<double> W_;
};

}  // namespace v1
}  // namespace geom
