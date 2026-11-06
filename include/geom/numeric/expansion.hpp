#pragma once

#include <cmath>
#include <utility>
#include <vector>

namespace geom {
inline namespace v1 {
namespace detail {

// Error-free transformations (Dekker / Knuth). For finite doubles the returned
// pair satisfies s + e == a op b exactly.
inline void two_sum(double a, double b, double& s, double& e) noexcept {
  s = a + b;
  const double bv = s - a;
  const double av = s - bv;
  e = (a - av) + (b - bv);
}

inline void two_diff(double a, double b, double& s, double& e) noexcept {
  s = a - b;
  const double bv = a - s;
  const double av = s + bv;
  e = (a - av) + (bv - b);
}

inline void two_product(double a, double b, double& p, double& e) noexcept {
  p = a * b;
  e = std::fma(a, b, -p);
}

// Exact sum of doubles stored as a non-overlapping expansion ordered by
// increasing magnitude. Only the four arithmetic primitives above are used,
// so every operation is exact for finite inputs.
class Expansion {
 public:
  Expansion() = default;
  explicit Expansion(double a) { add_term(a); }

  static Expansion difference(double a, double b) noexcept {
    double s = 0.0, e = 0.0;
    two_diff(a, b, s, e);
    Expansion r;
    r.add_term(e);
    r.add_term(s);
    return r;
  }

  void add_term(double b) {
    if (b == 0.0) return;
    std::vector<double> h;
    h.reserve(v_.size() + 1);
    double q = b;
    for (double e : v_) {
      double s = 0.0, t = 0.0;
      two_sum(q, e, s, t);
      h.push_back(t);
      q = s;
    }
    if (q != 0.0) h.push_back(q);
    v_.swap(h);
  }

  void add(const Expansion& e) { for (double x : e.v_) add_term(x); }

  void add_scaled(const Expansion& e, double s) {
    for (double x : e.v_) {
      double p = 0.0, err = 0.0;
      two_product(x, s, p, err);
      add_term(err);
      add_term(p);
    }
  }

  static Expansion mul(const Expansion& a, const Expansion& b) {
    Expansion r;
    for (double y : b.v_) r.add_scaled(a, y);
    return r;
  }

  int sign() const noexcept {
    for (std::size_t i = v_.size(); i > 0; --i) {
      const double x = v_[i - 1];
      if (x != 0.0) return x > 0.0 ? 1 : -1;
    }
    return 0;
  }

  double approx() const noexcept {
    double s = 0.0;
    for (double x : v_) s += x;
    return s;
  }

  std::size_t size() const noexcept { return v_.size(); }
  double operator[](std::size_t i) const noexcept { return v_[i]; }
  bool empty() const noexcept { return v_.empty(); }

 private:
  std::vector<double> v_;
};

}  // namespace detail
}  // namespace v1
}  // namespace geom
