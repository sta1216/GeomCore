#include <doctest/doctest.h>

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
