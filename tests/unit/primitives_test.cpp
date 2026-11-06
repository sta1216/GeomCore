#include <doctest/doctest.h>

#include <cmath>
#include <limits>

#include <geom/primitives/box.hpp>
#include <geom/primitives/point.hpp>
#include <geom/primitives/segment.hpp>
#include <geom/primitives/vec.hpp>
#include <geom/transform/affine.hpp>

using geom::v1::Affine;
using geom::v1::Box3d;
using geom::v1::Point3d;
using geom::v1::Segment3d;
using geom::v1::Vec3d;

TEST_CASE("point and vector affine algebra") {
  const Point3d p{1.0, 2.0, 3.0};
  const Vec3d v{0.5, -0.5, 1.0};
  const Point3d q = p + v;
  CHECK(q.v[0] == 1.5);
  CHECK(q.v[1] == 1.5);
  CHECK(q.v[2] == 4.0);

  const Vec3d d = q - p;
  CHECK(d.v[0] == v.v[0]);
  CHECK(geom::v1::dot(v, v) == doctest::Approx(1.5));
  CHECK(geom::v1::dot(geom::v1::cross(Vec3d{1, 0, 0}, Vec3d{0, 1, 0}), Vec3d{0, 0, 1}) ==
        doctest::Approx(1.0));
}

TEST_CASE("box is empty until expanded") {
  Box3d b = Box3d::empty();
  CHECK(b.is_empty());
  b.expand(Point3d{1.0, 2.0, 3.0});
  CHECK_FALSE(b.is_empty());
  CHECK(b.contains(Point3d{1.0, 2.0, 3.0}));
  b.expand(Point3d{-1.0, 0.0, 5.0});
  CHECK(b.lo.v[0] == -1.0);
  CHECK(b.hi.v[2] == 5.0);
}

TEST_CASE("closest point on a segment clamps to the endpoints") {
  const Segment3d s{Point3d{0, 0, 0}, Point3d{1, 0, 0}};
  CHECK(s.closest_point(Point3d{0.25, 1.0, 0.0}).v[0] == doctest::Approx(0.25));
  CHECK(s.closest_point(Point3d{-3.0, 0.0, 0.0}).v[0] == doctest::Approx(0.0));
  CHECK(s.closest_point(Point3d{3.0, 0.0, 0.0}).v[0] == doctest::Approx(1.0));
  CHECK(geom::v1::squared_distance(Point3d{0.5, 2.0, 0.0}, s) == doctest::Approx(4.0));
}

TEST_CASE("affine transforms round-trip through the inverse") {
  Affine a = Affine::identity();
  a = Affine::translation(Vec3d{1.0, 2.0, 3.0});
  const Point3d p{10.0, 20.0, 30.0};
  const Point3d q = a.transform(p);
  CHECK(q.v[0] == 11.0);

  const auto inv = a.inverse();
  REQUIRE(inv.has_value());
  const Point3d back = inv.value().transform(q);
  CHECK(back.v[0] == doctest::Approx(p.v[0]));
  CHECK(back.v[1] == doctest::Approx(p.v[1]));
  CHECK(back.v[2] == doctest::Approx(p.v[2]));
}

TEST_CASE("vector algebra and norm") {
  const Vec3d a{1.0, 2.0, 3.0};
  const Vec3d b{4.0, 5.0, 6.0};
  const Vec3d sum = a + b;
  CHECK(sum.v[0] == 5.0);
  CHECK(sum.v[2] == 9.0);
  const Vec3d diff = b - a;
  CHECK(diff.v[1] == 3.0);
  const Vec3d scaled = a * 2.0;
  CHECK(scaled.v[2] == 6.0);
  CHECK(geom::v1::norm2(geom::v1::Vec2d{3.0, 4.0}) == doctest::Approx(5.0));
  CHECK(geom::v1::dot(geom::v1::Vec2d{1.0, 2.0}, geom::v1::Vec2d{3.0, 4.0}) ==
        doctest::Approx(11.0));
  // cross is right-handed and anticommutative.
  CHECK(geom::v1::cross(Vec3d{1.0, 0.0, 0.0}, Vec3d{0.0, 1.0, 0.0}).v[2] == 1.0);
  CHECK(geom::v1::cross(Vec3d{0.0, 1.0, 0.0}, Vec3d{1.0, 0.0, 0.0}).v[2] == -1.0);
}

TEST_CASE("box merge and containment") {
  Box3d a = Box3d::empty();
  a.expand(Point3d{0.0, 0.0, 0.0});
  a.expand(Point3d{2.0, 4.0, 6.0});
  Box3d b = Box3d::empty();
  b.expand(Point3d{-1.0, 1.0, 1.0});
  b.expand(Point3d{1.0, 2.0, 2.0});

  CHECK_FALSE(Box3d::empty().contains(Point3d{0.0, 0.0, 0.0}));

  const Box3d m = a.merged(b);
  CHECK(m.lo.v[0] == -1.0);
  CHECK(m.hi.v[2] == 6.0);
  CHECK(m.diagonal().v[0] == doctest::Approx(3.0));
  CHECK(m.diagonal().v[2] == doctest::Approx(6.0));
  CHECK(m.contains(Point3d{0.5, 1.5, 2.0}));
  CHECK_FALSE(m.contains(Point3d{0.5, 5.0, 2.0}));

  // NaN fails both comparisons and must not be reported as contained.
  const double nan = std::numeric_limits<double>::quiet_NaN();
  const double inf = std::numeric_limits<double>::infinity();
  CHECK_FALSE(m.contains(Point3d{nan, 1.5, 2.0}));
  CHECK_FALSE(m.contains(Point3d{0.5, nan, 2.0}));
  CHECK_FALSE(m.contains(Point3d{0.5, 1.5, nan}));
  CHECK_FALSE(m.contains(Point3d{inf, 1.5, 2.0}));
  CHECK_FALSE(m.contains(Point3d{-inf, 1.5, 2.0}));
}

TEST_CASE("degenerate segment collapses to its point") {
  const Segment3d s{Point3d{1.0, 2.0, 3.0}, Point3d{1.0, 2.0, 3.0}};
  const Point3d c = s.closest_point(Point3d{4.0, 6.0, 3.0});
  CHECK(c.v[0] == 1.0);
  CHECK(c.v[1] == 2.0);
  CHECK(geom::v1::squared_distance(Point3d{4.0, 6.0, 3.0}, s) == doctest::Approx(25.0));
}

TEST_CASE("affine failure and vector paths") {
  // transform_vector ignores the translation column.
  const Affine t = Affine::translation(Vec3d{1.0, 2.0, 3.0});
  const Vec3d v = t.transform_vector(Vec3d{1.0, 1.0, 1.0});
  CHECK(v.v[0] == 1.0);
  CHECK(v.v[1] == 1.0);
  CHECK(v.v[2] == 1.0);
  // identity leaves points unchanged.
  const Point3d p{4.0, 5.0, 6.0};
  const Point3d q = Affine::identity().transform(p);
  CHECK(q.v[0] == 4.0);
  CHECK(q.v[2] == 6.0);

  // Singular linear part -> Err::degenerate, not a NaN transform.
  Affine singular{};
  for (int i = 0; i < 3; ++i) {
    for (int j = 0; j < 3; ++j) singular.m[i][j] = 0.0;
  }
  const auto inv = singular.inverse();
  CHECK_FALSE(inv.has_value());
  CHECK(inv.error().code == geom::v1::Errc::degenerate);
  CHECK(inv.error().message != nullptr);

  // A determinant that overflows to infinity is rejected the same way.
  Affine overflow{};
  overflow.m[0][0] = 1e308;
  overflow.m[1][1] = 1e308;
  overflow.m[2][2] = 1e308;
  const auto inv2 = overflow.inverse();
  CHECK_FALSE(inv2.has_value());
  CHECK(inv2.error().code == geom::v1::Errc::degenerate);
}

TEST_CASE("affine inverse rejects non-finite entries in every slot") {
  using geom::v1::Errc;
  const double nan_v = std::numeric_limits<double>::quiet_NaN();
  const double inf_v = std::numeric_limits<double>::infinity();

  auto expect_non_finite = [](const Affine& a) {
    const auto inv = a.inverse();
    CHECK_FALSE(inv.has_value());
    if (inv.has_value()) return;
    CHECK(inv.error().code == Errc::non_finite);
  };

  // Translation column slots.
  for (int row = 0; row < 3; ++row) {
    Affine a{};
    a.m[row][3] = inf_v;
    expect_non_finite(a);
    Affine b{};
    b.m[row][3] = nan_v;
    expect_non_finite(b);
  }
  // Diagonal slots.
  for (int row = 0; row < 3; ++row) {
    Affine a{};
    a.m[row][row] = -inf_v;
    expect_non_finite(a);
  }
  // Off-diagonal (rotation) slots.
  Affine rotation{};
  rotation.m[0][1] = nan_v;
  expect_non_finite(rotation);
  Affine rotation2{};
  rotation2.m[2][0] = inf_v;
  expect_non_finite(rotation2);
}

TEST_CASE("affine inverse finite-output contract") {
  using geom::v1::Errc;

  // A subnormal determinant makes 1/det overflow: rejected, never a NaN matrix.
  Affine subnormal{};
  subnormal.m[0][0] = 1e-320;
  const auto inv_sub = subnormal.inverse();
  REQUIRE_FALSE(inv_sub.has_value());
  CHECK(inv_sub.error().code == Errc::degenerate);

  // Determinant products that underflow to 0 or overflow to inf are conservative
  // rejects even though the exact inverse would be representable.
  Affine underflow{};
  underflow.m[0][0] = 1e-200;
  underflow.m[1][1] = 1e-200;
  underflow.m[2][2] = 1e-200;
  const auto inv_under = underflow.inverse();
  REQUIRE_FALSE(inv_under.has_value());
  CHECK(inv_under.error().code == Errc::degenerate);

  Affine overflow{};
  overflow.m[0][0] = 1e300;
  overflow.m[1][1] = 1e300;
  overflow.m[2][2] = 1e300;
  const auto inv_over = overflow.inverse();
  REQUIRE_FALSE(inv_over.has_value());
  CHECK(inv_over.error().code == Errc::degenerate);

  // A finite determinant can still yield an unrepresentable cofactor ratio.
  Affine mixed{};
  mixed.m[0][0] = 1e-300;
  mixed.m[0][1] = 1e300;
  mixed.m[1][1] = 1.0;
  mixed.m[2][0] = 1.0;
  mixed.m[2][2] = 1.0;
  const auto inv_mixed = mixed.inverse();
  REQUIRE_FALSE(inv_mixed.has_value());
  CHECK(inv_mixed.error().code == Errc::degenerate);

  // One extreme-but-finite scale factor keeps unit determinant: still invertible.
  Affine small_scale{};
  small_scale.m[0][0] = 1e-200;
  const auto inv_small = small_scale.inverse();
  REQUIRE(inv_small.has_value());
  CHECK(std::isfinite(inv_small.value().m[0][0]));
  CHECK(inv_small.value().m[1][1] == doctest::Approx(1.0));

  Affine big_scale{};
  big_scale.m[0][0] = 1e300;
  const auto inv_big = big_scale.inverse();
  REQUIRE(inv_big.has_value());
  CHECK(std::isfinite(inv_big.value().m[0][0]));
  CHECK(inv_big.value().m[0][0] == doctest::Approx(1e-300));
}

TEST_CASE("affine rotation and scale inverse round-trips") {
  Affine a{};
  a.m[0][0] = 0.0;
  a.m[0][1] = -2.0;
  a.m[0][3] = 5.0;
  a.m[1][0] = 2.0;
  a.m[1][1] = 0.0;
  a.m[1][3] = -3.0;
  a.m[2][2] = 1.0;
  a.m[2][3] = 2.0;

  const auto inv = a.inverse();
  REQUIRE(inv.has_value());
  const Point3d p{1.0, 2.0, 3.0};
  const Point3d q = a.transform(p);
  const Point3d back = inv.value().transform(q);
  CHECK(back.v[0] == doctest::Approx(p.v[0]).epsilon(1e-12));
  CHECK(back.v[1] == doctest::Approx(p.v[1]).epsilon(1e-12));
  CHECK(back.v[2] == doctest::Approx(p.v[2]).epsilon(1e-12));

  const Vec3d v = a.transform_vector(Vec3d{1.0, 0.0, 0.0});
  CHECK(v.v[0] == 0.0);
  CHECK(v.v[1] == doctest::Approx(2.0));
}