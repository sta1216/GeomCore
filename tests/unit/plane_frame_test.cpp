#include <doctest/doctest.h>

#include <cmath>
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

TEST_CASE("plane_frame: huge finite normal components normalize safely") {
  const auto frame =
      make_plane_frame(Plane3{Vec3d{1.04e308, 1.04e308, 1.04e308}, 0.0});
  REQUIRE(frame.has_value());
  const PlaneFrame& f = frame.value();
  for (int i = 0; i < 3; ++i) {
    CHECK(std::isfinite(f.origin[i]));
    CHECK(std::isfinite(f.u[i]));
    CHECK(std::isfinite(f.v[i]));
    CHECK(std::isfinite(f.normal[i]));
  }
  CHECK(std::sqrt(dot(f.normal, f.normal)) == doctest::Approx(1.0).epsilon(1e-15));
  CHECK(dot(f.u, f.u) == doctest::Approx(1.0).epsilon(1e-15));
  CHECK(dot(f.v, f.v) == doctest::Approx(1.0).epsilon(1e-15));
  CHECK(dot(f.u, f.v) == doctest::Approx(0.0).epsilon(1e-15));
  const Vec3d w = cross(f.u, f.v);
  for (int i = 0; i < 3; ++i) {
    CHECK(w[i] == doctest::Approx(f.normal[i]).epsilon(1e-15));
  }
}

TEST_CASE("plane_frame: denormal normal components normalize safely") {
  const double d = std::numeric_limits<double>::denorm_min();
  const auto frame = make_plane_frame(Plane3{Vec3d{d, d, d}, 0.0});
  REQUIRE(frame.has_value());
  const PlaneFrame& f = frame.value();
  CHECK(dot(f.normal, f.normal) == doctest::Approx(1.0).epsilon(1e-12));
  CHECK(dot(f.u, f.u) == doctest::Approx(1.0).epsilon(1e-12));
  CHECK(dot(f.v, f.v) == doctest::Approx(1.0).epsilon(1e-12));
  CHECK(dot(f.u, f.v) == doctest::Approx(0.0).epsilon(1e-12));
  const Vec3d w = cross(f.u, f.v);
  for (int i = 0; i < 3; ++i) {
    CHECK(w[i] == doctest::Approx(f.normal[i]).epsilon(1e-12));
  }
}

TEST_CASE("plane_frame: huge offsets stay finite with a representable distance") {
  // x = 1e200 expressed with a unit normal.
  const auto unit = make_plane_frame(Plane3{Vec3d{1.0, 0.0, 0.0}, 1e200});
  REQUIRE(unit.has_value());
  CHECK(std::isfinite(unit.value().origin[0]));
  CHECK(unit.value().origin[0] == doctest::Approx(1e200).epsilon(1e-12));
  CHECK(unit.value().origin[1] == doctest::Approx(0.0).epsilon(1e-12));
  CHECK(unit.value().origin[2] == doctest::Approx(0.0).epsilon(1e-12));
  CHECK(plane_frame_plane(unit.value()).offset == doctest::Approx(1e200).epsilon(1e-12));

  // The same geometric plane with a tiny normal: offset / |normal| = 1e200.
  const auto scaled = make_plane_frame(Plane3{Vec3d{1e-300, 0.0, 0.0}, 1e-100});
  REQUIRE(scaled.has_value());
  CHECK(std::isfinite(scaled.value().origin[0]));
  CHECK(scaled.value().origin[0] == doctest::Approx(1e200).epsilon(1e-12));
  CHECK(plane_frame_plane(scaled.value()).offset == doctest::Approx(1e200).epsilon(1e-12));
}

TEST_CASE("plane_frame: non-finite plane offset is rejected") {
  const auto nan = make_plane_frame(
      Plane3{Vec3d{1.0, 0.0, 0.0}, std::numeric_limits<double>::quiet_NaN()});
  REQUIRE_FALSE(nan.has_value());
  CHECK(nan.error().code == Errc::non_finite);
}

TEST_CASE("plane_frame: off-plane point projects orthogonally") {
  const Plane3 plane{Vec3d{0.0, 0.0, 3.0}, 2.4};
  const auto frame = make_plane_frame(plane);
  REQUIRE(frame.has_value());
  const Point3d p{1.0, 2.0, 5.0};
  const Point3d back = plane_frame_to_3d(frame.value(), plane_frame_to_2d(frame.value(), p));
  const Vec3d n = frame.value().normal;
  const double offset = plane_frame_plane(frame.value()).offset;
  const double distance = dot(n, Vec3d{p[0], p[1], p[2]}) - offset;
  CHECK(back[0] == doctest::Approx(p[0] - distance * n[0]).epsilon(1e-12));
  CHECK(back[1] == doctest::Approx(p[1] - distance * n[1]).epsilon(1e-12));
  CHECK(back[2] == doctest::Approx(p[2] - distance * n[2]).epsilon(1e-12));
}
