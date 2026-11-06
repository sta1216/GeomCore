#pragma once

#include <array>
#include <cmath>
#include <cstddef>
#include <utility>
#include <vector>

#include <geom/curves/nurbs_curve.hpp>
#include <geom/primitives/point.hpp>
#include <geom/primitives/vec.hpp>

namespace geom {
inline namespace v1 {

// Tensor-product NURBS surface. Control points are indexed as i * nv + j
// (i along u, j along v). All evaluation routes through the 4D homogeneous
// de Boor helpers in <geom/curves/nurbs_curve.hpp>.
//
// Validity rules mirror NurbsCurve and apply to both directions: finite
// knots, finite control coordinates, finite positive weights, a nondegenerate
// domain, knot multiplicity <= degree + 1, and degree in 1..kMaxCurveDegree.
//
// Precondition: eval(), du() and dv() require is_valid(). Calling them on an
// invalid surface (degree above kMaxCurveDegree, mismatched control/weight
// sizes, non-finite data, degenerate domain, ...) is undefined behaviour, not
// a reported error. Parameters outside the domain are clamped to it.
// Evaluation allocates scratch vectors, so the noexcept markers are
// conservative for the M0 scope.
class NurbsSurface {
 public:
  NurbsSurface(int pu, int pv, std::vector<double> ku, std::vector<double> kv,
               std::vector<Point3d> control, std::vector<double> weights)
      : pu_(pu), pv_(pv), ku_(std::move(ku)), kv_(std::move(kv)),
        P_(std::move(control)), W_(std::move(weights)) {}

  bool is_valid() const noexcept {
    if (pu_ < 1 || pv_ < 1) return false;
    if (pu_ > detail::kMaxCurveDegree || pv_ > detail::kMaxCurveDegree) return false;
    if (ku_.size() < static_cast<std::size_t>(pu_) + 2) return false;
    if (kv_.size() < static_cast<std::size_t>(pv_) + 2) return false;
    const std::size_t nu_ = nu();
    const std::size_t nv_ = nv();
    if (nu_ < 2 || nv_ < 2) return false;
    if (P_.size() != nu_ * nv_ || W_.size() != P_.size()) return false;
    if (!valid_knots(ku_, pu_) || !valid_knots(kv_, pv_)) return false;
    if (!(ku_[static_cast<std::size_t>(pu_)] <
          ku_[ku_.size() - static_cast<std::size_t>(pu_) - 1])) {
      return false;
    }
    if (!(kv_[static_cast<std::size_t>(pv_)] <
          kv_[kv_.size() - static_cast<std::size_t>(pv_) - 1])) {
      return false;
    }
    for (const Point3d& p : P_) {
      if (!(std::isfinite(p[0]) && std::isfinite(p[1]) && std::isfinite(p[2]))) return false;
    }
    for (double w : W_) {
      if (!(std::isfinite(w) && w > 0.0)) return false;
    }
    return true;
  }

  std::size_t nu() const noexcept {
    return ku_.size() - static_cast<std::size_t>(pu_) - 1;
  }
  std::size_t nv() const noexcept {
    return kv_.size() - static_cast<std::size_t>(pv_) - 1;
  }

  // Domain queries return {0, 0} when is_valid() is false.
  std::pair<double, double> domain_u() const noexcept {
    if (!is_valid()) return {0.0, 0.0};
    return {ku_[static_cast<std::size_t>(pu_)],
            ku_[ku_.size() - static_cast<std::size_t>(pu_) - 1]};
  }
  std::pair<double, double> domain_v() const noexcept {
    if (!is_valid()) return {0.0, 0.0};
    return {kv_[static_cast<std::size_t>(pv_)],
            kv_[kv_.size() - static_cast<std::size_t>(pv_) - 1]};
  }

  Point3d eval(double u, double v) const noexcept {
    const auto h = eval_h(clamp_u(u), clamp_v(v));
    return Point3d{h[0] / h[3], h[1] / h[3], h[2] / h[3]};
  }

  Vec3d du(double u, double v) const noexcept {
    const double uu = clamp_u(u), vv = clamp_v(v);
    return project(eval_h(uu, vv), deriv_h(uu, vv, true));
  }

  Vec3d dv(double u, double v) const noexcept {
    const double uu = clamp_u(u), vv = clamp_v(v);
    return project(eval_h(uu, vv), deriv_h(uu, vv, false));
  }

 private:
  // Finite knots, nondecreasing order, and no knot repeated more than p + 1
  // times. Finiteness is checked before ordering so infinities cannot slip in.
  static bool valid_knots(const std::vector<double>& k, int p) noexcept {
    for (double x : k) {
      if (!std::isfinite(x)) return false;
    }
    int run = 1;
    for (std::size_t i = 1; i < k.size(); ++i) {
      if (!(k[i] >= k[i - 1])) return false;
      if (k[i] == k[i - 1]) {
        ++run;
        if (run > p + 1) return false;
      } else {
        run = 1;
      }
    }
    return true;
  }

  double clamp_u(double u) const noexcept {
    const auto [a, b] = domain_u();
    return u < a ? a : (u > b ? b : u);
  }
  double clamp_v(double v) const noexcept {
    const auto [a, b] = domain_v();
    return v < a ? a : (v > b ? b : v);
  }

  std::array<double, 4> hom(std::size_t i, std::size_t j) const noexcept {
    const std::size_t idx = i * nv() + j;
    return {W_[idx] * P_[idx][0], W_[idx] * P_[idx][1], W_[idx] * P_[idx][2], W_[idx]};
  }

  std::vector<std::array<double, 4>> rows(double v) const {
    std::vector<std::array<double, 4>> R(nu());
    std::vector<std::array<double, 4>> row(nv());
    for (std::size_t i = 0; i < nu(); ++i) {
      for (std::size_t j = 0; j < nv(); ++j) row[j] = hom(i, j);
      R[i] = detail::eval_curve(kv_, pv_, row, v);
    }
    return R;
  }

  std::vector<std::array<double, 4>> cols(double u) const {
    std::vector<std::array<double, 4>> C(nv());
    std::vector<std::array<double, 4>> col(nu());
    for (std::size_t j = 0; j < nv(); ++j) {
      for (std::size_t i = 0; i < nu(); ++i) col[i] = hom(i, j);
      C[j] = detail::eval_curve(ku_, pu_, col, u);
    }
    return C;
  }

  std::array<double, 4> eval_h(double u, double v) const noexcept {
    const auto R = rows(v);
    return detail::eval_curve(ku_, pu_, R, u);
  }

  std::array<double, 4> deriv_h(double u, double v, bool along_u) const noexcept {
    if (along_u) {
      const auto R = rows(v);
      std::vector<std::array<double, 4>> D(R.size() - 1);
      for (std::size_t i = 0; i + 1 < R.size(); ++i) {
        const double denom = ku_[i + static_cast<std::size_t>(pu_) + 1] - ku_[i + 1];
        const double f = denom == 0.0 ? 0.0 : static_cast<double>(pu_) / denom;
        for (int c = 0; c < 4; ++c) {
          D[i][static_cast<std::size_t>(c)] =
              f * (R[i + 1][static_cast<std::size_t>(c)] - R[i][static_cast<std::size_t>(c)]);
        }
      }
      std::vector<double> U1(ku_.begin() + 1, ku_.end() - 1);
      return detail::eval_curve(U1, pu_ - 1, D, u);
    }
    const auto C = cols(u);
    std::vector<std::array<double, 4>> D(C.size() - 1);
    for (std::size_t j = 0; j + 1 < C.size(); ++j) {
      const double denom = kv_[j + static_cast<std::size_t>(pv_) + 1] - kv_[j + 1];
      const double f = denom == 0.0 ? 0.0 : static_cast<double>(pv_) / denom;
      for (int c = 0; c < 4; ++c) {
        D[j][static_cast<std::size_t>(c)] =
            f * (C[j + 1][static_cast<std::size_t>(c)] - C[j][static_cast<std::size_t>(c)]);
      }
    }
    std::vector<double> V1(kv_.begin() + 1, kv_.end() - 1);
    return detail::eval_curve(V1, pv_ - 1, D, v);
  }

  static Vec3d project(const std::array<double, 4>& h,
                       const std::array<double, 4>& dh) noexcept {
    const double w = h[3];
    const Vec3d c{h[0] / w, h[1] / w, h[2] / w};
    return Vec3d{(dh[0] - c[0] * dh[3]) / w, (dh[1] - c[1] * dh[3]) / w,
                 (dh[2] - c[2] * dh[3]) / w};
  }

  int pu_ = 0;
  int pv_ = 0;
  std::vector<double> ku_;
  std::vector<double> kv_;
  std::vector<Point3d> P_;
  std::vector<double> W_;
};

}  // namespace v1
}  // namespace geom
