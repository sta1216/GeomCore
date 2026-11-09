#include <doctest/doctest.h>

#include <limits>

#include <geom/brep/plane_frame.hpp>

using namespace geom;

TEST_CASE("plane_frame: orthonormal right-handed frame on a tilted plane") {
  const Vec3d n{1.0, 2.0, 2.0};  // |n| = 3
  const Plane3 plane{n, 7.0};
  const auto frame = make_plane_frame(plane);
  REQUIRE(frame.has_value());
  const PlaneFrame& f = frame.value();
  CHECK(dot(f.normal, f.normal) == doctest::Approx(1.0).epsilon(1e-15));
  CHECK(dot(f.u, f.u) == doctest::Approx(1.0).epsilon(1e-15));
  CHECK(dot(f.v, f.v) == doctest::Approx(1.0).epsilon(1e-15));
  CHECK(dot(f.u, f.v) == doctest::Approx(0.0).epsilon(1e-15));
  CHECK(dot(cross(f.u, f.v), f.normal) == doctest::Approx(1.0).epsilon(1e-15));
  CHECK(f.origin[0] + 2.0 * f.origin[1] + 2.0 * f.origin[2] ==
        doctest::Approx(7.0).epsilon(1e-15));
}

TEST_CASE("plane_frame: 3d to 2d round trip is exact on the plane") {
  const Plane3 plane{Vec3d{0.0, 1.0, 0.0}, -2.0};  // y = -2
  const auto frame = make_plane_frame(plane);
  REQUIRE(frame.has_value());
  for (const Point3d p : {Point3d{1.0, -2.0, 0.5}, Point3d{-3.25, -2.0, 4.0},
                          Point3d{0.0, -2.0, 0.0}}) {
    const Point3d back = plane_frame_to_3d(frame.value(), plane_frame_to_2d(frame.value(), p));
    CHECK(back[0] == doctest::Approx(p[0]).epsilon(1e-14));
    CHECK(back[1] == doctest::Approx(p[1]).epsilon(1e-14));
    CHECK(back[2] == doctest::Approx(p[2]).epsilon(1e-14));
  }
}

TEST_CASE("plane_frame: degenerate and non-finite planes are rejected") {
  const auto zero = make_plane_frame(Plane3{Vec3d{0.0, 0.0, 0.0}, 0.0});
  REQUIRE_FALSE(zero.has_value());
  CHECK(zero.error().code == Errc::degenerate);
  const auto nan = make_plane_frame(Plane3{Vec3d{0.0, 0.0, std::numeric_limits<double>::quiet_NaN()}, 0.0});
  REQUIRE_FALSE(nan.has_value());
  CHECK(nan.error().code == Errc::non_finite);
}

TEST_CASE("plane_frame: the basis only depends on the normal direction") {
  const Vec3d n{3.0, 4.0, 0.0};
  const auto a = make_plane_frame(Plane3{n, 1.0});
  const auto b = make_plane_frame(Plane3{n, 9.0});
  REQUIRE(a.has_value());
  REQUIRE(b.has_value());
  CHECK(dot(a.value().u, b.value().u) == doctest::Approx(1.0).epsilon(1e-15));
  CHECK(dot(a.value().v, b.value().v) == doctest::Approx(1.0).epsilon(1e-15));
  CHECK(plane_frame_plane(a.value()).offset == doctest::Approx(0.2).epsilon(1e-15));
  CHECK(plane_frame_plane(b.value()).offset == doctest::Approx(1.8).epsilon(1e-15));
}
