#pragma once

#include <algorithm>
#include <cmath>
#include <limits>

#include <geom/numeric/expansion.hpp>

namespace geom {
inline namespace v1 {
namespace detail {

constexpr double kEpsilon = std::numeric_limits<double>::epsilon() * 0.5;  // 2^-53
constexpr double kCcwErrBoundA = (3.0 + 16.0 * kEpsilon) * kEpsilon;

// 精确符号：det = (ax-cx)*(by-cy) - (ay-cy)*(bx-cx)，用 Expansion 计算。
//
// 精确性前提：本路径只在没有任何中间量上溢/下溢时才是精确的。Expansion 的
// 每一项都必须是有限 double；若某个乘积（或展开后的分量）落到次正规区之外，
// 精确性即失效。行列式量级低于约 2^-1074 或因乘积下溢而丢失全部有效位时，
// 结果可能被报告为 0（保守的“退化”判定），而不是抛错。快速过滤路径则保持
// 其通常的误差界（kCcwErrBoundA * detsum），只在过滤无法判定时才进入这里。
inline int orient2d_exact(double ax, double ay, double bx, double by, double cx, double cy) noexcept {
  const Expansion acx = Expansion::difference(ax, cx);
  const Expansion acy = Expansion::difference(ay, cy);
  const Expansion bcx = Expansion::difference(bx, cx);
  const Expansion bcy = Expansion::difference(by, cy);

  Expansion det = Expansion::mul(acx, bcy);
  det.add_scaled(Expansion::mul(acy, bcx), -1.0);
  return det.sign();
}

}  // namespace detail

// 符号函数：返回 -1/0/1。调用方保证有限输入；非有限输入返回 0。
//
// 精确性前提（与 detail::orient2d_exact 相同）：仅当所有中间量都不上溢/下溢时
// 精确。量级低于约 2^-1074 或乘积下溢的行列式可能被报告为 0；不会静默给出
// 错误的非零符号。过滤路径保持通常误差界，仅在无法判定时才回退到精确路径，
// 因此“接近退化”的输入也按精确符号分类。
inline int sign_orient2d(double ax, double ay, double bx, double by, double cx,
                         double cy) noexcept {
  const double m = std::max({std::abs(ax), std::abs(ay), std::abs(bx), std::abs(by),
                             std::abs(cx), std::abs(cy)});
  if (!std::isfinite(m)) return 0;
  if (m == 0.0) return 0;

  // 保持尺度在安全区间内：按 2 的幂缩放不改变行列式符号（齐次 3 次）。
  constexpr double kMaxSafe = 1e150;
  constexpr double kMinSafe = 1e-150;
  double scale = 1.0;
  if (m > kMaxSafe || m < kMinSafe) {
    int exp2 = 0;
    std::frexp(m, &exp2);
    scale = std::ldexp(1.0, -exp2);
  }
  ax *= scale; ay *= scale; bx *= scale; by *= scale; cx *= scale; cy *= scale;

  const double detleft = (ax - cx) * (by - cy);
  const double detright = (ay - cy) * (bx - cx);
  const double det = detleft - detright;
  const double detsum = std::abs(detleft) + std::abs(detright);
  const double errbound = detail::kCcwErrBoundA * detsum;

  if (det > errbound) return 1;
  if (-det > errbound) return -1;
  return detail::orient2d_exact(ax, ay, bx, by, cx, cy);
}

}  // namespace v1
}  // namespace geom
