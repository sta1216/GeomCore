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
// 精确性前提：本路径只在没有任何中间量上溢/下溢时才是精确的。Expansion 的每一项
// 都必须是有限 double；若某个乘积（或展开后的分量）落到次正规区之外，精确性即失效。
// 行列式量级低于约 2^-1074，或因乘积下溢而丢失全部有效位时，结果可能被报告为 0
// （保守的“退化”判定），而不是抛错。快速过滤路径则保持其通常的误差界
// （kCcwErrBoundA * detsum），只在过滤无法判定时才进入这里。
inline int orient2d_exact(double ax, double ay, double bx, double by, double cx, double cy) noexcept {
  const Expansion acx = Expansion::difference(ax, cx);
  const Expansion acy = Expansion::difference(ay, cy);
  const Expansion bcx = Expansion::difference(bx, cx);
  const Expansion bcy = Expansion::difference(by, cy);

  Expansion det = Expansion::mul(acx, bcy);
  det.add_scaled(Expansion::mul(acy, bcx), -1.0);
  return det.sign();
}

// 在给定尺度上的分类结果。rescale 表示 double 过滤本身已不可用：乘积上溢，或两个
// 乘积都下溢为 0；此时只能换个尺度重新分类。
enum class OrientDecision { pos, neg, exact, rescale };

// 浮点过滤：det > errbound 时符号可信（pos/neg），否则 exact 交给精确路径；
// 过滤本身失效（上溢/全下溢）时返回 rescale，表示“当前尺度不可用”。
inline OrientDecision classify_orient2d(double ax, double ay, double bx, double by,
                                        double cx, double cy) noexcept {
  const double detleft = (ax - cx) * (by - cy);
  const double detright = (ay - cy) * (bx - cx);
  const double det = detleft - detright;
  const double detsum = std::abs(detleft) + std::abs(detright);
  if (!std::isfinite(det) || !std::isfinite(detsum)) return OrientDecision::rescale;
  if (detsum == 0.0) return OrientDecision::rescale;  // 两个乘积都下溢了
  const double errbound = kCcwErrBoundA * detsum;
  if (det > errbound) return OrientDecision::pos;
  if (-det > errbound) return OrientDecision::neg;
  return OrientDecision::exact;
}

}  // namespace detail

// 符号函数：返回 -1/0/1。六个坐标中任意一个非有限（NaN 或 ±inf）都返回 0，与它
// 出现在第几个参数无关。
//
// 精确性前提：对有限输入，只要缩放后的行列式不低于约 2^-1074（即展开后的每个分量
// 仍是有限 double），结果就是精确的。本函数不会返回错误的非零符号：当下溢使行列式
// 无法表示时，最坏情况是保守地返回 0。
//
// 先按原始尺度分类，仅当 double 过滤上溢、或两个乘积都下溢为 0 时才按 2 的幂重缩放。
// 未缩放优先是必要的：坐标量级跨度超过 300 个十进位的输入，一旦缩放就会把小坐标
// 推下溢，而未缩放时过滤依然能给出正确的符号。
inline int sign_orient2d(double ax, double ay, double bx, double by, double cx,
                         double cy) noexcept {
  // 逐个坐标检查：std::max({...}) 对 NaN 只保留首个元素，NaN 出现在后面的位置时
  // 会被静默跳过，随后 Expansion 会把 NaN 当成负项符号化。
  if (!(std::isfinite(ax) && std::isfinite(ay) && std::isfinite(bx) && std::isfinite(by) &&
        std::isfinite(cx) && std::isfinite(cy))) {
    return 0;
  }
  const double m = std::max({std::abs(ax), std::abs(ay), std::abs(bx), std::abs(by),
                             std::abs(cx), std::abs(cy)});
  if (m == 0.0) return 0;

  switch (detail::classify_orient2d(ax, ay, bx, by, cx, cy)) {
    case detail::OrientDecision::pos:
      return 1;
    case detail::OrientDecision::neg:
      return -1;
    case detail::OrientDecision::exact:
      return detail::orient2d_exact(ax, ay, bx, by, cx, cy);
    case detail::OrientDecision::rescale:
      break;
  }

  // 行列式关于坐标是 2 次齐次的，所以按 2 的幂缩放不改变符号。
  // 夹紧 shift 是下溢安全的关键：次正规输入（exp2 <= -1024）若直接取 -exp2，
  // ldexp(1.0, -exp2) 会得到 +inf，坐标随即变成 NaN/0，符号就全错了。
  int exp2 = 0;
  std::frexp(m, &exp2);
  int shift = -exp2;
  if (shift > 1023) shift = 1023;  // 保证 s 有限，ldexp 绝不上溢
  const double s = std::ldexp(1.0, shift);

  ax *= s; ay *= s; bx *= s; by *= s; cx *= s; cy *= s;

  switch (detail::classify_orient2d(ax, ay, bx, by, cx, cy)) {
    case detail::OrientDecision::pos:
      return 1;
    case detail::OrientDecision::neg:
      return -1;
    case detail::OrientDecision::exact:
    case detail::OrientDecision::rescale:
      return detail::orient2d_exact(ax, ay, bx, by, cx, cy);
  }
  return 0;
}

}  // namespace v1
}  // namespace geom
