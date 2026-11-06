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

namespace detail {

constexpr double kIncircleErrBoundA = (10.0 + 96.0 * kEpsilon) * kEpsilon;

// 精确符号：d 相对过 a,b,c 的圆的位置，a,b,c 逆时针时圆内为正。
//
// 精确性前提：与 orient2d_exact 相同，只有所有中间量都保持有限、既不上溢也不下溢时才
// 精确。行列式是 4 次齐次的，可表示范围比 orient2d 窄得多：|坐标| 到 1e75 量级时升降项
// 已经接近 1e300。若展开后的分量落到次正规区之外，结果可能被保守地报告为 0。
// 精确行列式本身（返回 Expansion），符号由调用者取。分开暴露是为了让 sign_incircle
// 能在缩放帧丢掉行列式时检查未缩放帧的展开是否仍然全部有限，再决定是否采信它。
inline Expansion incircle_det(double ax, double ay, double bx, double by, double cx, double cy,
                              double dx, double dy) noexcept {
  const Expansion adx = Expansion::difference(ax, dx);
  const Expansion ady = Expansion::difference(ay, dy);
  const Expansion bdx = Expansion::difference(bx, dx);
  const Expansion bdy = Expansion::difference(by, dy);
  const Expansion cdx = Expansion::difference(cx, dx);
  const Expansion cdy = Expansion::difference(cy, dy);

  Expansion alift = Expansion::mul(adx, adx);
  alift.add(Expansion::mul(ady, ady));
  Expansion blift = Expansion::mul(bdx, bdx);
  blift.add(Expansion::mul(bdy, bdy));
  Expansion clift = Expansion::mul(cdx, cdx);
  clift.add(Expansion::mul(cdy, cdy));

  Expansion cross_bc = Expansion::mul(bdx, cdy);
  cross_bc.add_scaled(Expansion::mul(cdx, bdy), -1.0);
  Expansion cross_ac = Expansion::mul(adx, cdy);
  cross_ac.add_scaled(Expansion::mul(cdx, ady), -1.0);
  Expansion cross_ab = Expansion::mul(adx, bdy);
  cross_ab.add_scaled(Expansion::mul(bdx, ady), -1.0);

  Expansion det = Expansion::mul(alift, cross_bc);
  det.add_scaled(Expansion::mul(blift, cross_ac), -1.0);
  det.add_scaled(Expansion::mul(clift, cross_ab), 1.0);
  return det;
}

inline int incircle_exact(double ax, double ay, double bx, double by, double cx, double cy,
                          double dx, double dy) noexcept {
  return incircle_det(ax, ay, bx, by, cx, cy, dx, dy).sign();
}

// 同一行列式的“分布式”展开：把 alift*cross_bc 这类项拆成 6 个链式三项乘积
// （例如 clift*adx*bdy），而不是先算 cross 再乘 lift。两种形式在精确算术下相等，
// 但乘法顺序不同：分组形式在 cross 内部就先相消，深处相消时这个中间量可能整个
// 掉到 2^-1074 以下；分布式形式把每个乘积单独保留下来，能救回这些符号。
// 只在分组形式返回 0 时作为第二层回退使用。
inline Expansion incircle_det_distributed(double ax, double ay, double bx, double by, double cx,
                                          double cy, double dx, double dy) noexcept {
  const Expansion adx = Expansion::difference(ax, dx);
  const Expansion ady = Expansion::difference(ay, dy);
  const Expansion bdx = Expansion::difference(bx, dx);
  const Expansion bdy = Expansion::difference(by, dy);
  const Expansion cdx = Expansion::difference(cx, dx);
  const Expansion cdy = Expansion::difference(cy, dy);

  Expansion alift = Expansion::mul(adx, adx);
  alift.add(Expansion::mul(ady, ady));
  Expansion blift = Expansion::mul(bdx, bdx);
  blift.add(Expansion::mul(bdy, bdy));
  Expansion clift = Expansion::mul(cdx, cdx);
  clift.add(Expansion::mul(cdy, cdy));

  Expansion det = Expansion::mul(Expansion::mul(alift, bdx), cdy);
  det.add_scaled(Expansion::mul(Expansion::mul(alift, cdx), bdy), -1.0);
  det.add_scaled(Expansion::mul(Expansion::mul(blift, adx), cdy), -1.0);
  det.add_scaled(Expansion::mul(Expansion::mul(blift, cdx), ady), 1.0);
  det.add_scaled(Expansion::mul(Expansion::mul(clift, adx), bdy), 1.0);
  det.add_scaled(Expansion::mul(Expansion::mul(clift, bdx), ady), -1.0);
  return det;
}

// 在给定尺度上的分类结果，语义与 OrientDecision 相同：rescale 表示 double 过滤本身
// 已经不可用——升降项上溢、全部下溢为 0，或误差界下溢到非正规区。
enum class IncircleDecision { pos, neg, exact, rescale };

inline IncircleDecision classify_incircle(double ax, double ay, double bx, double by, double cx,
                                          double cy, double dx, double dy) noexcept {
  const double adx = ax - dx, ady = ay - dy;
  const double bdx = bx - dx, bdy = by - dy;
  const double cdx = cx - dx, cdy = cy - dy;
  const double alift = adx * adx + ady * ady;
  const double blift = bdx * bdx + bdy * bdy;
  const double clift = cdx * cdx + cdy * cdy;
  const double det = alift * (bdx * cdy - cdx * bdy) - blift * (adx * cdy - cdx * ady) +
                     clift * (adx * bdy - bdx * ady);
  const double permanent =
      (std::abs(bdx * cdy) + std::abs(cdx * bdy)) * alift +
      (std::abs(adx * cdy) + std::abs(cdx * ady)) * blift +
      (std::abs(adx * bdy) + std::abs(bdx * ady)) * clift;
  if (!std::isfinite(det) || !std::isfinite(permanent)) return IncircleDecision::rescale;
  if (permanent == 0.0) return IncircleDecision::rescale;  // 全部下溢为 0
  const double errbound = kIncircleErrBoundA * permanent;
  // 误差界必须是正规 double。permanent 落在次正规区时 errbound 会下溢为 0，而此时的 det
  // 本身也已是次正规数，舍入误差与它同量级，任何非零比较都可能给出错误符号。这种情形
  // 必须换一个尺度，不能当作精确判定。
  if (!(errbound >= std::numeric_limits<double>::min())) return IncircleDecision::rescale;
  if (det > errbound) return IncircleDecision::pos;
  if (-det > errbound) return IncircleDecision::neg;
  return IncircleDecision::exact;
}

}  // namespace detail

// incircle 符号：返回 -1/0/1。正号表示 d 在以 a,b,c 为顶点的圆内（a,b,c 逆时针）。
// 八个坐标中任意一个非有限（NaN 或 ±inf）都返回 0，与它出现在第几个参数无关。
//
// 数值前提：过滤把 det 的舍入误差限制在 kIncircleErrBoundA * permanent 之内，前提是
// permanent 是正规 double，从而误差界本身也是正规数。行列式关于坐标是 4 次齐次的，
// 可表示范围比 orient2d 窄：一次全局缩放最多把 max|坐标| 抬到 2^240（见下方 rescale），
// 只要 max|坐标| / min|非零坐标| 不超过约 2^1314（约 396 个十进制数量级），
// 缩放后每个非零坐标都仍可表示且在缩放帧里可精确求值，符号就是精确的。
//
// 返回 0 的三种情况，都只表示“未知/退化”，绝不是错误符号：
//   (a) 行列式本身超出 double 范围（> DBL_MAX 或 < 最小次正规数），两个求值帧都
//       无法表示它；
//   (b) 缩放把某个非零坐标冲成 0，缩放后的构型已经不是输入构型；这只可能发生在
//       max|坐标| > 2^240 且 min|非零坐标| < 2^(exp2(max) - 1315) 时；
//   (c) 坐标的全局动态范围超过上面那段上限；
// 在支持的范围内，有限输入总是得到精确符号；范围之外返回 0 而不是猜测的符号，
// 因此本函数在任何情况下都不会返回错误的非零符号。
//
// 先按原始尺度分类，仅当产物上溢、全部下溢为 0 或误差界不是正规数时才按 2 的幂重缩放。
// 未缩放优先是必要的：1e200 与 1e100 这类混合量级一旦直接缩放，小坐标就会下溢，未缩放时
// 过滤仍能给出符号。
//
// 缩放帧的展开式会在深处相消时丢掉行列式：某个分量乘积在主项相消之前就下溢到
// 2^-1074 以下。回退顺序是“分组展开 → 分布式展开 → 未缩放帧”：分布式把 lift 乘进
// 三项乘积，不让 cross 先相消；未缩放帧只在它的展开全部有限时才采信（溢出的展开
// 会含 inf/NaN，而 sign() 会把 NaN 读成负号）。任一层给出的非零符号都是精确符号。
inline int sign_incircle(double ax, double ay, double bx, double by, double cx, double cy,
                         double dx, double dy) noexcept {
  if (!(std::isfinite(ax) && std::isfinite(ay) && std::isfinite(bx) && std::isfinite(by) &&
        std::isfinite(cx) && std::isfinite(cy) && std::isfinite(dx) && std::isfinite(dy))) {
    return 0;
  }
  const double m = std::max({std::abs(ax), std::abs(ay), std::abs(bx), std::abs(by),
                             std::abs(cx), std::abs(cy), std::abs(dx), std::abs(dy)});
  if (m == 0.0) return 0;

  switch (detail::classify_incircle(ax, ay, bx, by, cx, cy, dx, dy)) {
    case detail::IncircleDecision::pos:
      return 1;
    case detail::IncircleDecision::neg:
      return -1;
    case detail::IncircleDecision::exact:
      return detail::incircle_exact(ax, ay, bx, by, cx, cy, dx, dy);
    case detail::IncircleDecision::rescale:
      break;
  }

  // 行列式关于坐标是 4 次齐次的，按统一的 2 的幂缩放不改变符号。取 shift = 240 - exp2
  // 而不是 -exp2：缩放后 max|坐标| <= 2^240，差值 <= 2^241，lift 与 cross <= 2^483，
  // 每个 4 次项 <= 2^966，permanent <= 2^968，全部有限；同时 shift 为正时小坐标被放大，
  // 而不是像“把最大坐标归一到 1”那样被冲掉（后者会把次正规坐标直接压成 0）。
  // 夹紧到 [-1074, 1023]：次正规输入的 exp2 约 -1073，shift 可到 1313，ldexp 会得到 +inf。
  int exp2 = 0;
  std::frexp(m, &exp2);
  int shift = 240 - exp2;
  if (shift > 1023) shift = 1023;
  if (shift < -1074) shift = -1074;
  const double s = std::ldexp(1.0, shift);

  const double sax = ax * s, say = ay * s;
  const double sbx = bx * s, sby = by * s;
  const double scx = cx * s, scy = cy * s;
  const double sdx = dx * s, sdy = dy * s;

  // 缩放把某个非零坐标冲成 0 时，缩放后的构型已经不是输入构型，任何非零结论都不可信；
  // 返回保守的 0，而不是让被裁剪过的构型决定符号。
  if ((sax == 0.0 && ax != 0.0) || (say == 0.0 && ay != 0.0) || (sbx == 0.0 && bx != 0.0) ||
      (sby == 0.0 && by != 0.0) || (scx == 0.0 && cx != 0.0) || (scy == 0.0 && cy != 0.0) ||
      (sdx == 0.0 && dx != 0.0) || (sdy == 0.0 && dy != 0.0)) {
    return 0;
  }

  switch (detail::classify_incircle(sax, say, sbx, sby, scx, scy, sdx, sdy)) {
    case detail::IncircleDecision::pos:
      return 1;
    case detail::IncircleDecision::neg:
      return -1;
    case detail::IncircleDecision::exact:
    case detail::IncircleDecision::rescale: {
      // 缩放帧的精确计算只在中间量不下溢时精确：某个分量乘积可能在主导项相消
      // 之前掉到 2^-1074 以下，使 Expansion 丢掉整个行列式。三层回退，任一层的
      // 非零符号都是精确符号：
      //   1. 缩放帧的分组展开（最快）；
      //   2. 缩放帧的分布式展开——cross 项不再先相消，单个三项乘积能保住符号；
      //   3. 未缩放帧（检查全部有限，避免把溢出的 NaN 当成负号）。
      // 三层都为 0 才返回保守的 0。
      const detail::Expansion scaled =
          detail::incircle_det(sax, say, sbx, sby, scx, scy, sdx, sdy);
      const int scaled_sign = scaled.sign();
      if (scaled_sign != 0) return scaled_sign;

      const detail::Expansion distributed =
          detail::incircle_det_distributed(sax, say, sbx, sby, scx, scy, sdx, sdy);
      const int distributed_sign = distributed.sign();
      if (distributed_sign != 0) return distributed_sign;

      const detail::Expansion unscaled = detail::incircle_det(ax, ay, bx, by, cx, cy, dx, dy);
      return unscaled.all_finite() ? unscaled.sign() : 0;
    }
  }
  return 0;
}

}  // namespace v1
}  // namespace geom
