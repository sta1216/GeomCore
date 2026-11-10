#include <doctest/doctest.h>

#include <array>
#include <cmath>
#include <geom/brep/detail/ssi_solver.hpp>
#include <limits>

using namespace geom::v1;

TEST_CASE("ssi: Newton refines a seed on a plane-plane intersection") {
  const Result<SsiPatch> a =
      make_ssi_patch(SurfaceGeometry{Plane3{{0, 0, 1}, 0}}, {-2, 2}, {-2, 2});
  const Result<SsiPatch> b =
      make_ssi_patch(SurfaceGeometry{Plane3{{1, 0, 0}, 0}}, {-2, 2}, {-2, 2});
  REQUIRE(a.has_value());
  REQUIRE(b.has_value());
  const Result<SsiPoint> p = refine_ssi_point(a.value(), b.value(), {0.3, 0.4, 0.2, 0.9}, 1e-10);
  REQUIRE(p.has_value());
  CHECK(p.value().residual <= 1e-10);
  CHECK(std::fabs(p.value().p[0]) < 1e-9);
  CHECK(std::fabs(p.value().p[2]) < 1e-9);
  const Result<SurfacePoint> on_a =
      surface_eval(SurfaceGeometry{Plane3{{0, 0, 1}, 0}}, p.value().ua[0], p.value().ua[1]);
  const Result<SurfacePoint> on_b =
      surface_eval(SurfaceGeometry{Plane3{{1, 0, 0}, 0}}, p.value().va[0], p.value().va[1]);
  REQUIRE(on_a.has_value());
  REQUIRE(on_b.has_value());
  for (int i = 0; i < 3; ++i) {
    CHECK(std::fabs(on_a.value().p[i] - p.value().p[i]) < 1e-9);
    CHECK(std::fabs(on_b.value().p[i] - p.value().p[i]) < 1e-9);
  }
}

TEST_CASE("ssi: Newton converges on a sphere-plane circle") {
  const Sphere3 sphere{Point3d{0, 0, 0}, Vec3d{0, 0, 1}, Vec3d{1, 0, 0}, 1.0};
  const Result<SsiPatch> a =
      make_ssi_patch(SurfaceGeometry{sphere}, {0, 2 * 3.141592653589793}, {-1.5, 1.5});
  const Result<SsiPatch> b =
      make_ssi_patch(SurfaceGeometry{Plane3{{0, 0, 1}, 0.25}}, {-2, 2}, {-2, 2});
  REQUIRE(a.has_value());
  REQUIRE(b.has_value());
  const Result<SsiPoint> p = refine_ssi_point(a.value(), b.value(), {0.1, 0.2, 0.3, 0.4}, 1e-12);
  REQUIRE(p.has_value());
  CHECK(p.value().residual <= 1e-12);
  CHECK(std::fabs(p.value().p[2] - 0.25) < 1e-9);
  const Point3d& point = p.value().p;
  const double radius = std::sqrt(point[0] * point[0] + point[1] * point[1] + point[2] * point[2]);
  CHECK(std::fabs(radius - 1.0) < 1e-9);
}

TEST_CASE("ssi: divergent seeds are reported as numeric_failure") {
  const Result<SsiPatch> a =
      make_ssi_patch(SurfaceGeometry{Plane3{{0, 0, 1}, 0}}, {-1, 1}, {-1, 1});
  const Result<SsiPatch> b =
      make_ssi_patch(SurfaceGeometry{Plane3{{0, 0, 1}, 5}}, {-1, 1}, {-1, 1});
  REQUIRE(a.has_value());
  REQUIRE(b.has_value());
  const Result<SsiPoint> p = refine_ssi_point(a.value(), b.value(), {0.1, 0.1, 0.1, 0.1}, 1e-12, 4);
  REQUIRE_FALSE(p.has_value());
  CHECK(p.error().code == Errc::numeric_failure);
}

TEST_CASE("ssi: invalid tolerances are rejected before any math") {
  const Result<SsiPatch> a =
      make_ssi_patch(SurfaceGeometry{Plane3{{0, 0, 1}, 0}}, {-2, 2}, {-2, 2});
  const Result<SsiPatch> b =
      make_ssi_patch(SurfaceGeometry{Plane3{{1, 0, 0}, 0}}, {-2, 2}, {-2, 2});
  REQUIRE(a.has_value());
  REQUIRE(b.has_value());
  const std::array<double, 6> invalid{0.0,
                                      -0.0,
                                      -1e-9,
                                      std::numeric_limits<double>::quiet_NaN(),
                                      std::numeric_limits<double>::infinity(),
                                      -std::numeric_limits<double>::infinity()};
  for (const double tolerance : invalid) {
    const Result<SsiPoint> p =
        refine_ssi_point(a.value(), b.value(), {0.3, 0.4, 0.2, 0.9}, tolerance);
    REQUIRE_FALSE(p.has_value());
    CHECK(p.error().code == Errc::tolerance_violation);
  }
  // The tolerance is the first contract checked: a bad seed cannot mask it.
  const Result<SsiPoint> ordered = refine_ssi_point(
      a.value(), b.value(), {std::numeric_limits<double>::quiet_NaN(), 0.0, 0.0, 0.0}, 0.0);
  REQUIRE_FALSE(ordered.has_value());
  CHECK(ordered.error().code == Errc::tolerance_violation);
}

TEST_CASE("ssi: non-finite seeds are rejected") {
  const Result<SsiPatch> a =
      make_ssi_patch(SurfaceGeometry{Plane3{{0, 0, 1}, 0}}, {-2, 2}, {-2, 2});
  const Result<SsiPatch> b =
      make_ssi_patch(SurfaceGeometry{Plane3{{1, 0, 0}, 0}}, {-2, 2}, {-2, 2});
  REQUIRE(a.has_value());
  REQUIRE(b.has_value());
  for (int index = 0; index < 4; ++index) {
    std::array<double, 4> seed{0.3, 0.4, 0.2, 0.9};
    seed[static_cast<std::size_t>(index)] = std::numeric_limits<double>::quiet_NaN();
    const Result<SsiPoint> p = refine_ssi_point(a.value(), b.value(), seed, 1e-10);
    REQUIRE_FALSE(p.has_value());
    CHECK(p.error().code == Errc::non_finite);
    seed[static_cast<std::size_t>(index)] = std::numeric_limits<double>::infinity();
    const Result<SsiPoint> q = refine_ssi_point(a.value(), b.value(), seed, 1e-10);
    REQUIRE_FALSE(q.has_value());
    CHECK(q.error().code == Errc::non_finite);
  }
}

TEST_CASE("ssi: a non-finite surface evaluation is rejected") {
  // Valid cylinder whose evaluation overflows: origin + radial * radius is inf.
  const Cylinder3 huge{Point3d{1e308, 0, 0}, Vec3d{0, 0, 1}, Vec3d{1, 0, 0}, 1e308};
  REQUIRE(huge.is_valid());
  const Result<SsiPatch> a = make_ssi_patch(SurfaceGeometry{huge}, {-1, 1}, {-1, 1});
  const Result<SsiPatch> b =
      make_ssi_patch(SurfaceGeometry{Plane3{{0, 0, 1}, 0}}, {-2, 2}, {-2, 2});
  REQUIRE(a.has_value());
  REQUIRE(b.has_value());
  const Result<SsiPoint> p = refine_ssi_point(a.value(), b.value(), {0.0, 0.0, 0.0, 0.0}, 1e-12);
  REQUIRE_FALSE(p.has_value());
  CHECK(p.error().code == Errc::non_finite);
}

TEST_CASE("ssi: a numerically singular 3x3 system is rejected") {
  // The two planes meet at an angle of 1e-14: the four tangential columns are
  // nearly rank two, so every 3x3 minor is below the relative pivot threshold.
  const Result<SsiPatch> a =
      make_ssi_patch(SurfaceGeometry{Plane3{{0, 0, 1}, 0}}, {-2, 2}, {-2, 2});
  const Result<SsiPatch> b =
      make_ssi_patch(SurfaceGeometry{Plane3{{-1e-14, 0, 1}, 0}}, {-2, 2}, {-2, 2});
  REQUIRE(a.has_value());
  REQUIRE(b.has_value());
  const Result<SsiPoint> p = refine_ssi_point(a.value(), b.value(), {0.3, 0.4, 0.2, 0.9}, 1e-12);
  REQUIRE_FALSE(p.has_value());
  CHECK(p.error().code == Errc::numeric_failure);
}

TEST_CASE("ssi: non-positive iteration budgets are rejected") {
  const Result<SsiPatch> a =
      make_ssi_patch(SurfaceGeometry{Plane3{{0, 0, 1}, 0}}, {-2, 2}, {-2, 2});
  const Result<SsiPatch> b =
      make_ssi_patch(SurfaceGeometry{Plane3{{1, 0, 0}, 0}}, {-2, 2}, {-2, 2});
  REQUIRE(a.has_value());
  REQUIRE(b.has_value());
  for (const int budget : {0, -1, -24}) {
    const Result<SsiPoint> p =
        refine_ssi_point(a.value(), b.value(), {0.3, 0.4, 0.2, 0.9}, 1e-10, budget);
    REQUIRE_FALSE(p.has_value());
    CHECK(p.error().code == Errc::numeric_failure);
  }
}

TEST_CASE("ssi: a seed outside either patch is rejected without clamping") {
  const Result<SsiPatch> a =
      make_ssi_patch(SurfaceGeometry{Plane3{{0, 0, 1}, 0}}, {-1, 1}, {-1, 1});
  const Result<SsiPatch> b =
      make_ssi_patch(SurfaceGeometry{Plane3{{1, 0, 0}, 0}}, {-1, 1}, {-1, 1});
  REQUIRE(a.has_value());
  REQUIRE(b.has_value());
  // Outside A in v: without the seed check the first step would move v into the
  // patch and converge, so a success here would prove the seed was clamped.
  const Result<SsiPoint> outside_a =
      refine_ssi_point(a.value(), b.value(), {0.2, 1.5, 0.2, 0.9}, 1e-10);
  REQUIRE_FALSE(outside_a.has_value());
  CHECK(outside_a.error().code == Errc::numeric_failure);
  // Outside B in t, again on a component the solve would move back inside.
  const Result<SsiPoint> outside_b =
      refine_ssi_point(a.value(), b.value(), {0.2, 0.2, 0.2, 1.5}, 1e-10);
  REQUIRE_FALSE(outside_b.has_value());
  CHECK(outside_b.error().code == Errc::numeric_failure);
  // The closed domain is inclusive: a boundary seed is accepted.
  const Result<SsiPoint> boundary =
      refine_ssi_point(a.value(), b.value(), {1.0, 0.2, 0.2, 0.9}, 1e-10);
  REQUIRE(boundary.has_value());
  CHECK(boundary.value().residual <= 1e-10);
  CHECK(patch_contains(a.value(), boundary.value().ua[0], boundary.value().ua[1]));
  CHECK(patch_contains(b.value(), boundary.value().va[0], boundary.value().va[1]));
}

TEST_CASE("ssi: refinement is deterministic and wraps a periodic parameter") {
  const Sphere3 sphere{Point3d{0, 0, 0}, Vec3d{0, 0, 1}, Vec3d{1, 0, 0}, 1.0};
  const Result<SsiPatch> a =
      make_ssi_patch(SurfaceGeometry{sphere}, {0, 2 * 3.141592653589793}, {-1.5, 1.5});
  const Result<SsiPatch> b =
      make_ssi_patch(SurfaceGeometry{Plane3{{0, 0, 1}, 0.25}}, {-2, 2}, {-2, 2});
  REQUIRE(a.has_value());
  REQUIRE(b.has_value());
  // u = 7.0 is the same angle as 7.0 - 2*pi, but outside [-pi, pi): the stored
  // parameter must be folded, not carried as 7.0.
  const Result<SsiPoint> first =
      refine_ssi_point(a.value(), b.value(), {7.0, 0.2, 0.3, 0.4}, 1e-12);
  REQUIRE(first.has_value());
  CHECK(first.value().residual <= 1e-12);
  CHECK(std::fabs(first.value().ua[0]) <= 3.141592653589793);
  CHECK(std::fabs(first.value().ua[1]) <= 1.5);
  CHECK(patch_contains(a.value(), first.value().ua[0], first.value().ua[1]));
  const Result<SurfacePoint> on_a =
      surface_eval(a.value().surface, first.value().ua[0], first.value().ua[1]);
  REQUIRE(on_a.has_value());
  for (int i = 0; i < 3; ++i) {
    CHECK(std::fabs(on_a.value().p[i] - first.value().p[i]) < 1e-12);
  }
  CHECK(std::fabs(ssi_residual(a.value(), b.value(), first.value()) - first.value().residual) <=
        1e-15);
  const Result<SsiPoint> second =
      refine_ssi_point(a.value(), b.value(), {7.0, 0.2, 0.3, 0.4}, 1e-12);
  REQUIRE(second.has_value());
  for (int i = 0; i < 2; ++i) {
    CHECK(first.value().ua[i] == second.value().ua[i]);
    CHECK(first.value().va[i] == second.value().va[i]);
  }
  for (int i = 0; i < 3; ++i) {
    CHECK(first.value().p[i] == second.value().p[i]);
  }
  CHECK(first.value().residual == second.value().residual);
}

TEST_CASE("ssi: both torus parameters are wrapped into the stored domain") {
  const Torus3 torus{Point3d{0, 0, 0}, Vec3d{0, 0, 1}, Vec3d{1, 0, 0}, 2.0, 0.5};
  const Result<SsiPatch> a =
      make_ssi_patch(SurfaceGeometry{torus}, {-3.141592653589793, 3.141592653589793},
                     {-3.141592653589793, 3.141592653589793});
  const Result<SsiPatch> b =
      make_ssi_patch(SurfaceGeometry{Plane3{{0, 0, 1}, 0}}, {-3, 3}, {-3, 3});
  REQUIRE(a.has_value());
  REQUIRE(b.has_value());
  // Both torus axes start one turn away from their stored representatives.
  const Result<SsiPoint> p = refine_ssi_point(a.value(), b.value(), {7.0, 6.5, 0.2, 0.9}, 1e-10);
  REQUIRE(p.has_value());
  CHECK(p.value().residual <= 1e-10);
  CHECK(std::fabs(p.value().ua[0]) <= 3.141592653589793);
  CHECK(std::fabs(p.value().ua[1]) <= 3.141592653589793);
  CHECK(std::fabs(p.value().p[2]) < 1e-9);
  CHECK(patch_contains(a.value(), p.value().ua[0], p.value().ua[1]));
  const Result<SurfacePoint> on_a =
      surface_eval(a.value().surface, p.value().ua[0], p.value().ua[1]);
  REQUIRE(on_a.has_value());
  for (int i = 0; i < 3; ++i) {
    CHECK(std::fabs(on_a.value().p[i] - p.value().p[i]) < 1e-10);
  }
}

TEST_CASE("ssi: residual measures the gap and both point-to-surface distances") {
  const Result<SsiPatch> a =
      make_ssi_patch(SurfaceGeometry{Plane3{{0, 0, 1}, 0}}, {-2, 2}, {-2, 2});
  const Result<SsiPatch> b =
      make_ssi_patch(SurfaceGeometry{Plane3{{1, 0, 0}, 0}}, {-2, 2}, {-2, 2});
  REQUIRE(a.has_value());
  REQUIRE(b.has_value());
  // The stored p is the exact intersection point, but the parameters are not:
  // a residual that trusts p, or that returns a constant, is caught here.
  SsiPoint off;
  off.ua = {0.5, 0.5};
  off.va = {0.5, 0.5};
  off.p = {0.0, 0.3, 0.0};
  const double off_residual = ssi_residual(a.value(), b.value(), off);
  CHECK(std::isfinite(off_residual));
  CHECK(std::fabs(off_residual - std::sqrt(1.5)) < 1e-12);
  CHECK(off_residual > 1.0);
  // Genuine intersection parameters have zero residual.
  SsiPoint hit;
  hit.ua = {-0.3, 0.0};
  hit.va = {0.0, 0.3};
  hit.p = {0.0, 0.3, 0.0};
  CHECK(ssi_residual(a.value(), b.value(), hit) <= 1e-12);
}
