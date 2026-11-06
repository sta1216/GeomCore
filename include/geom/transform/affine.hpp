#pragma once

#include <cmath>

#include <geom/core/result.hpp>
#include <geom/primitives/point.hpp>
#include <geom/primitives/vec.hpp>

namespace geom {
inline namespace v1 {

// 3x4 affine transform: y = M * x + t, stored row-major in m[3][4].
struct Affine {
  double m[3][4] = {{1.0, 0.0, 0.0, 0.0}, {0.0, 1.0, 0.0, 0.0}, {0.0, 0.0, 1.0, 0.0}};

  static constexpr Affine identity() noexcept { return Affine{}; }

  static constexpr Affine translation(Vec3d t) noexcept {
    Affine a{};
    a.m[0][3] = t[0];
    a.m[1][3] = t[1];
    a.m[2][3] = t[2];
    return a;
  }

  constexpr Point3d transform(Point3d p) const noexcept {
    Point3d r;
    for (int i = 0; i < 3; ++i) {
      r[i] = m[i][0] * p[0] + m[i][1] * p[1] + m[i][2] * p[2] + m[i][3];
    }
    return r;
  }

  constexpr Vec3d transform_vector(Vec3d v) const noexcept {
    Vec3d r;
    for (int i = 0; i < 3; ++i) {
      r[i] = m[i][0] * v[0] + m[i][1] * v[1] + m[i][2] * v[2];
    }
    return r;
  }

  // Returns the inverse transform, or an error:
  //  - Errc::non_finite if any of the 12 stored entries is not finite;
  //  - Errc::degenerate if the determinant is zero, non-finite, or too small
  //    for its reciprocal (1/det overflows), or if any inverse entry is
  //    non-finite.
  // Success implies every entry of the returned matrix is finite.
  //
  // The determinant test is conservative: products that underflow to 0 or
  // overflow to +/-inf (for example diag(1e-200,1e-200,1e-200) or
  // diag(1e300,1e300,1e300)) are reported degenerate even when the exact
  // inverse would be representable. No condition-number estimate is made, so a
  // near-singular but finite matrix may produce a finite yet inaccurate
  // inverse; checking accuracy is the caller's responsibility.
  Result<Affine> inverse() const noexcept {
    for (int row = 0; row < 3; ++row) {
      for (int col = 0; col < 4; ++col) {
        if (!std::isfinite(m[row][col])) {
          return Error{Errc::non_finite, "non-finite affine transform entry", 0u, 0u};
        }
      }
    }

    const double a = m[0][0], b = m[0][1], c = m[0][2];
    const double d = m[1][0], e = m[1][1], f = m[1][2];
    const double g = m[2][0], h = m[2][1], i = m[2][2];

    const double A = e * i - f * h;
    const double B = -(d * i - f * g);
    const double C = d * h - e * g;
    const double det = a * A + b * B + c * C;
    if (det == 0.0 || !std::isfinite(det)) {
      return Error{Errc::degenerate, "singular affine transform", 0u, 0u};
    }
    const double inv_det = 1.0 / det;
    if (!std::isfinite(inv_det)) {
      return Error{Errc::degenerate, "affine determinant too small for inverse", 0u, 0u};
    }

    Affine r{};
    r.m[0][0] = A * inv_det;
    r.m[0][1] = -(b * i - c * h) * inv_det;
    r.m[0][2] = (b * f - c * e) * inv_det;
    r.m[1][0] = B * inv_det;
    r.m[1][1] = (a * i - c * g) * inv_det;
    r.m[1][2] = -(a * f - c * d) * inv_det;
    r.m[2][0] = C * inv_det;
    r.m[2][1] = -(a * h - b * g) * inv_det;
    r.m[2][2] = (a * e - b * d) * inv_det;

    for (int row = 0; row < 3; ++row) {
      r.m[row][3] = -(r.m[row][0] * m[0][3] + r.m[row][1] * m[1][3] + r.m[row][2] * m[2][3]);
    }

    for (int row = 0; row < 3; ++row) {
      for (int col = 0; col < 4; ++col) {
        if (!std::isfinite(r.m[row][col])) {
          return Error{Errc::degenerate, "non-finite affine inverse entry", 0u, 0u};
        }
      }
    }
    return r;
  }
};

}  // namespace v1
}  // namespace geom
