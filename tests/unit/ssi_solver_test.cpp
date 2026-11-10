#include <doctest/doctest.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <geom/brep/detail/ssi_solver.hpp>
#include <limits>
#include <string_view>
#include <vector>

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
  // The tolerance also outranks a non-positive iteration budget.
  const Result<SsiPoint> over_budget =
      refine_ssi_point(a.value(), b.value(), {0.3, 0.4, 0.2, 0.9}, 0.0, 0);
  REQUIRE_FALSE(over_budget.has_value());
  CHECK(over_budget.error().code == Errc::tolerance_violation);
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
  // The budget outranks a non-finite seed: the contract order is tolerance,
  // budget, seed, patches.
  const Result<SsiPoint> ordered = refine_ssi_point(
      a.value(), b.value(), {std::numeric_limits<double>::quiet_NaN(), 0.0, 0.0, 0.0}, 1e-10, 0);
  REQUIRE_FALSE(ordered.has_value());
  CHECK(ordered.error().code == Errc::numeric_failure);
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

TEST_CASE("ssi: residual is the parameter gap and ignores the stored point") {
  const Result<SsiPatch> a =
      make_ssi_patch(SurfaceGeometry{Plane3{{0, 0, 1}, 0}}, {-2, 2}, {-2, 2});
  const Result<SsiPatch> b =
      make_ssi_patch(SurfaceGeometry{Plane3{{1, 0, 0}, 0}}, {-2, 2}, {-2, 2});
  REQUIRE(a.has_value());
  REQUIRE(b.has_value());
  // The stored p is the exact intersection point, but the parameters are not:
  // a residual that trusts p, or that returns a constant, is caught here. The
  // residual is only the parameter gap: it certifies neither the correspondence
  // of the two parameter pairs nor their patch membership.
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

TEST_CASE("ssi: residual is infinite when a surface evaluation fails") {
  // R and r are both finite and 0 < r < R, but R + r overflows in eval(), so
  // this is a genuine evaluation failure rather than a rejected patch.
  const Torus3 huge{Point3d{0, 0, 0}, Vec3d{0, 0, 1}, Vec3d{1, 0, 0}, 1.5e308, 1.4e308};
  const Result<SsiPatch> a =
      make_ssi_patch(SurfaceGeometry{huge}, {0, 2 * 3.141592653589793}, {0, 2 * 3.141592653589793});
  REQUIRE(a.has_value());
  SsiPoint point;
  point.ua = {0.0, 0.0};
  point.va = {0.0, 0.0};
  const double residual = ssi_residual(a.value(), a.value(), point);
  CHECK(std::isinf(residual));
  CHECK(residual > 0.0);
}

TEST_CASE("ssi: the Newton step does not depend on the model scale") {
  // The Jacobian columns are divided by their largest entry so the 3x3 minors
  // stay O(1); the solved step is in those scaled units and must be converted
  // back. Without that conversion the step is `scale` times too large and the
  // line search cannot recover it, so a large or tiny model stalls.
  auto gap_between = [](const Result<SurfacePoint>& first, const Result<SurfacePoint>& second) {
    double squared = 0.0;
    for (int i = 0; i < 3; ++i) {
      const double difference = first.value().p[i] - second.value().p[i];
      squared += difference * difference;
    }
    return std::sqrt(squared);
  };
  for (const double radius : {1e-3, 1e3}) {
    const double tolerance = 1e-12 * radius;
    const Sphere3 first{Point3d{0, 0, 0}, Vec3d{0, 0, 1}, Vec3d{1, 0, 0}, radius};
    const Sphere3 second{Point3d{radius, 0, 0}, Vec3d{0, 0, 1}, Vec3d{1, 0, 0}, radius};
    const Result<SsiPatch> a =
        make_ssi_patch(SurfaceGeometry{first}, {0, 2 * 3.141592653589793}, {-1.5, 1.5});
    const Result<SsiPatch> b =
        make_ssi_patch(SurfaceGeometry{second}, {0, 2 * 3.141592653589793}, {-1.5, 1.5});
    REQUIRE(a.has_value());
    REQUIRE(b.has_value());
    const Result<SsiPoint> p =
        refine_ssi_point(a.value(), b.value(), {1.2472, 0.1, 1.9444, -0.1}, tolerance);
    REQUIRE(p.has_value());
    CHECK(p.value().residual <= tolerance);
    const Result<SurfacePoint> on_a =
        surface_eval(a.value().surface, p.value().ua[0], p.value().ua[1]);
    const Result<SurfacePoint> on_b =
        surface_eval(b.value().surface, p.value().va[0], p.value().va[1]);
    REQUIRE(on_a.has_value());
    REQUIRE(on_b.has_value());
    CHECK(gap_between(on_a, on_b) <= tolerance);
    CHECK(patch_contains(a.value(), p.value().ua[0], p.value().ua[1]));
    CHECK(patch_contains(b.value(), p.value().va[0], p.value().va[1]));
  }
  // The same factor scale for a curved/planar pair at a large size.
  const double radius = 1e3;
  const Sphere3 sphere{Point3d{0, 0, 0}, Vec3d{0, 0, 1}, Vec3d{1, 0, 0}, radius};
  const Result<SsiPatch> a =
      make_ssi_patch(SurfaceGeometry{sphere}, {0, 2 * 3.141592653589793}, {-1.5, 1.5});
  const Result<SsiPatch> b = make_ssi_patch(SurfaceGeometry{Plane3{{0, 0, 1}, 0.25 * radius}},
                                            {-2 * radius, 2 * radius}, {-2 * radius, 2 * radius});
  REQUIRE(a.has_value());
  REQUIRE(b.has_value());
  const Result<SsiPoint> p =
      refine_ssi_point(a.value(), b.value(), {0.1, 0.2, -97.83, 975.2}, 1e-12 * radius);
  REQUIRE(p.has_value());
  CHECK(p.value().residual <= 1e-12 * radius);
  const Result<SurfacePoint> on_a =
      surface_eval(a.value().surface, p.value().ua[0], p.value().ua[1]);
  const Result<SurfacePoint> on_b =
      surface_eval(b.value().surface, p.value().va[0], p.value().va[1]);
  REQUIRE(on_a.has_value());
  REQUIRE(on_b.has_value());
  CHECK(gap_between(on_a, on_b) <= 1e-12 * radius);
  CHECK(patch_contains(a.value(), p.value().ua[0], p.value().ua[1]));
  CHECK(patch_contains(b.value(), p.value().va[0], p.value().va[1]));
}

TEST_CASE("ssi: a solution across the periodic seam is stored folded") {
  // The raw seed parameter is just past +pi and the intersection point next to
  // it stays there, so the solution's raw u is in (pi, pi + 0.1) and the stored
  // representative must be folded into (-pi, -3.0). A solver that omits the
  // fold returns the raw value and fails the sign, the bound, and the patch
  // containment check.
  const Sphere3 sphere{Point3d{0, 0, 0}, Vec3d{0, 0, 1}, Vec3d{1, 0, 0}, 1.0};
  const Result<SsiPatch> a =
      make_ssi_patch(SurfaceGeometry{sphere}, {0, 2 * 3.141592653589793}, {-1.5, 1.5});
  const Result<SsiPatch> b =
      make_ssi_patch(SurfaceGeometry{Plane3{{1, 0, 0}, -0.5}}, {-2, 2}, {-2, 2});
  REQUIRE(a.has_value());
  REQUIRE(b.has_value());
  const Result<SsiPoint> p = refine_ssi_point(a.value(), b.value(), {3.2, 1.0, -0.8, -0.3}, 1e-10);
  REQUIRE(p.has_value());
  CAPTURE(p.value().ua[0]);
  CAPTURE(p.value().ua[1]);
  CHECK(p.value().residual <= 1e-10);
  CHECK(p.value().ua[0] < -3.0);
  CHECK(p.value().ua[0] >= -3.141592653589793);
  CHECK(patch_contains(a.value(), p.value().ua[0], p.value().ua[1]));
  const Result<SurfacePoint> on_a =
      surface_eval(a.value().surface, p.value().ua[0], p.value().ua[1]);
  const Result<SurfacePoint> on_b =
      surface_eval(b.value().surface, p.value().va[0], p.value().va[1]);
  REQUIRE(on_a.has_value());
  REQUIRE(on_b.has_value());
  for (int i = 0; i < 3; ++i) {
    CHECK(std::fabs(on_a.value().p[i] - on_b.value().p[i]) <= 1e-10);
  }
}

TEST_CASE("ssi: a seed that is already a solution returns without a step") {
  // The spheres touch the plane at one point and the seed is that point exactly,
  // so the convergence test at the top of the loop must succeed before any
  // Newton step (tangency is otherwise a numeric_failure).
  const Sphere3 sphere{Point3d{0, 0, 1}, Vec3d{0, 0, 1}, Vec3d{1, 0, 0}, 1.0};
  const Result<SsiPatch> a =
      make_ssi_patch(SurfaceGeometry{sphere}, {0, 2 * 3.141592653589793}, {-2, 2});
  const Result<SsiPatch> b =
      make_ssi_patch(SurfaceGeometry{Plane3{{0, 0, 1}, 0}}, {-2, 2}, {-2, 2});
  REQUIRE(a.has_value());
  REQUIRE(b.has_value());
  const Result<SsiPoint> p =
      refine_ssi_point(a.value(), b.value(), {0.5, -1.5707963267948966, 0.0, 0.0}, 1e-12);
  REQUIRE(p.has_value());
  CHECK(p.value().residual <= 1e-12);
  CHECK(p.value().ua[0] == 0.5);
  CHECK(p.value().ua[1] == -1.5707963267948966);
  CHECK(p.value().va[0] == 0.0);
  CHECK(p.value().va[1] == 0.0);
  CHECK(std::fabs(p.value().p[2]) <= 1e-15);
}

TEST_CASE("ssi: a displaced tangency seed does not converge") {
  const Sphere3 sphere{Point3d{0, 0, 1}, Vec3d{0, 0, 1}, Vec3d{1, 0, 0}, 1.0};
  const Result<SsiPatch> a =
      make_ssi_patch(SurfaceGeometry{sphere}, {0, 2 * 3.141592653589793}, {-1.5, 1.5});
  const Result<SsiPatch> b =
      make_ssi_patch(SurfaceGeometry{Plane3{{0, 0, 1}, 0}}, {-2, 2}, {-2, 2});
  REQUIRE(a.has_value());
  REQUIRE(b.has_value());
  const Result<SsiPoint> p =
      refine_ssi_point(a.value(), b.value(), {0.5, -1.4707963267948966, 0.1, 0.1}, 1e-12);
  REQUIRE_FALSE(p.has_value());
  CHECK(p.error().code == Errc::numeric_failure);
}

namespace {

constexpr double kTracePi = 3.141592653589793;

SsiPoint traced_endpoint(const SsiSegment& segment, bool second) {
  SsiPoint point;
  point.ua = second ? segment.ub : segment.ua;
  point.va = second ? segment.vb : segment.va;
  point.p = second ? segment.b : segment.a;
  return point;
}

void check_segment_certified(const SsiPatch& a, const SsiPatch& b, const SsiSegment& segment,
                             double tolerance) {
  const SsiPoint first = traced_endpoint(segment, false);
  const SsiPoint second = traced_endpoint(segment, true);
  CHECK(patch_contains(a, first.ua[0], first.ua[1]));
  CHECK(patch_contains(a, second.ua[0], second.ua[1]));
  CHECK(patch_contains(b, first.va[0], first.va[1]));
  CHECK(patch_contains(b, second.va[0], second.va[1]));
  CHECK(segment.residual >= 0.0);
  CHECK(segment.residual <= tolerance);
  // The stored residual is gap-only: the maximum of the two endpoint residuals,
  // not a distance that could stand in for patch membership.
  CHECK(segment.residual == std::max(ssi_residual(a, b, first), ssi_residual(a, b, second)));
}

}  // namespace

TEST_CASE("ssi: trace closes the sphere-plane circle") {
  const Sphere3 sphere{Point3d{0, 0, 0}, Vec3d{0, 0, 1}, Vec3d{1, 0, 0}, 1.0};
  const Result<SsiPatch> a =
      make_ssi_patch(SurfaceGeometry{sphere}, {0, 2 * kTracePi}, {-1.5, 1.5});
  const Result<SsiPatch> b =
      make_ssi_patch(SurfaceGeometry{Plane3{{0, 0, 1}, 0}}, {-2, 2}, {-2, 2});
  REQUIRE(a.has_value());
  REQUIRE(b.has_value());
  const Result<SsiPoint> seed = refine_ssi_point(a.value(), b.value(), {0.1, 0.2, 0.3, 0.4}, 1e-12);
  REQUIRE(seed.has_value());
  const Result<std::vector<SsiCurve>> traced =
      trace_ssi_curves(a.value(), b.value(), {seed.value()}, 1e-9);
  REQUIRE(traced.has_value());
  REQUIRE(traced.value().size() == 1);

  const SsiCurve& curve = traced.value().front();
  CHECK(curve.closed);
  CHECK(curve.segments.size() >= 8);
  for (std::size_t i = 0; i < curve.segments.size(); ++i) {
    const SsiSegment& segment = curve.segments[i];
    check_segment_certified(a.value(), b.value(), segment, 1e-9);
    // Every endpoint is on the unit circle in the plane z = 0. A vertex that
    // was emitted from a raw prediction (or a stale seed) shows up here.
    const double radius_a = std::sqrt(segment.a[0] * segment.a[0] + segment.a[1] * segment.a[1] +
                                      segment.a[2] * segment.a[2]);
    const double radius_b = std::sqrt(segment.b[0] * segment.b[0] + segment.b[1] * segment.b[1] +
                                      segment.b[2] * segment.b[2]);
    CHECK(std::fabs(radius_a - 1.0) < 1e-8);
    CHECK(std::fabs(radius_b - 1.0) < 1e-8);
    CHECK(std::fabs(segment.a[2]) < 1e-8);
    CHECK(std::fabs(segment.b[2]) < 1e-8);
    // Every emitted vertex is itself a refined point: re-refining its stored
    // parameters reproduces it, so no unrefined prediction leaked into a
    // segment.
    const SsiPoint first = traced_endpoint(segment, false);
    const Result<SsiPoint> again = refine_ssi_point(
        a.value(), b.value(), {first.ua[0], first.ua[1], first.va[0], first.va[1]}, 1e-9);
    REQUIRE(again.has_value());
    for (int axis = 0; axis < 3; ++axis) CHECK(again.value().p[axis] == first.p[axis]);
    // Adjacent segments share an endpoint bitwise, and the closing segment ends
    // at the start vertex without duplicating it as an extra zero-length
    // segment.
    const SsiSegment& next = curve.segments[(i + 1) % curve.segments.size()];
    for (int axis = 0; axis < 3; ++axis) CHECK(segment.b[axis] == next.a[axis]);
    CHECK(segment.ub == next.ua);
    CHECK(segment.vb == next.va);
    // A duplicated closing vertex would show up as a zero-length segment.
    CHECK(norm2(segment.b - segment.a) > 0.0);
  }
}

TEST_CASE("ssi: trace stops at a bounded patch boundary") {
  const Sphere3 sphere{Point3d{0, 0, 0}, Vec3d{0, 0, 1}, Vec3d{1, 0, 0}, 1.0};
  const Result<SsiPatch> a = make_ssi_patch(SurfaceGeometry{sphere}, {-0.5, 0.5}, {-0.5, 0.5});
  const Result<SsiPatch> b =
      make_ssi_patch(SurfaceGeometry{Plane3{{0, 0, 1}, 0}}, {-2, 2}, {-2, 2});
  REQUIRE(a.has_value());
  REQUIRE(b.has_value());
  const Result<SsiPoint> seed = refine_ssi_point(a.value(), b.value(), {0.1, 0.1, 0.2, 0.2}, 1e-12);
  REQUIRE(seed.has_value());
  const Result<std::vector<SsiCurve>> traced =
      trace_ssi_curves(a.value(), b.value(), {seed.value()}, 1e-9);
  REQUIRE(traced.has_value());
  REQUIRE(traced.value().size() == 1);

  const SsiCurve& curve = traced.value().front();
  CHECK_FALSE(curve.closed);
  CHECK(curve.segments.size() >= 2);
  for (const SsiSegment& segment : curve.segments) {
    check_segment_certified(a.value(), b.value(), segment, 1e-9);
  }
  // The curve must end where the patch does, not at a predicted point past the
  // boundary: both open ends are refined vertices inside both patches.
  const SsiSegment& first = curve.segments.front();
  const SsiSegment& last = curve.segments.back();
  CHECK(patch_contains(a.value(), first.ua[0], first.ua[1]));
  CHECK(patch_contains(b.value(), first.va[0], first.va[1]));
  CHECK(patch_contains(a.value(), last.ub[0], last.ub[1]));
  CHECK(patch_contains(b.value(), last.vb[0], last.vb[1]));
  CHECK(std::fabs(ssi_residual(a.value(), b.value(), traced_endpoint(last, true))) <= 1e-9);
  // The arc is bounded by the narrow A patch (u in [-0.5, 0.5], v = 0), so its
  // u coordinates must not run past the domain.
  for (const SsiSegment& segment : curve.segments) {
    CHECK(segment.ua[0] >= -0.5);
    CHECK(segment.ua[0] <= 0.5);
    CHECK(segment.ub[0] >= -0.5);
    CHECK(segment.ub[0] <= 0.5);
  }
}

TEST_CASE("ssi: traced segments are contiguous") {
  const Result<SsiPatch> a =
      make_ssi_patch(SurfaceGeometry{Plane3{{0, 0, 1}, 0}}, {-2, 2}, {-2, 2});
  const Result<SsiPatch> b =
      make_ssi_patch(SurfaceGeometry{Plane3{{1, 0, 0}, 0}}, {-2, 2}, {-2, 2});
  REQUIRE(a.has_value());
  REQUIRE(b.has_value());
  const Result<SsiPoint> seed = refine_ssi_point(a.value(), b.value(), {0.1, 0.1, 0.1, 0.1}, 1e-12);
  REQUIRE(seed.has_value());
  const Result<std::vector<SsiCurve>> traced =
      trace_ssi_curves(a.value(), b.value(), {seed.value()}, 1e-9, 64);
  REQUIRE(traced.has_value());
  REQUIRE(traced.value().size() == 1);

  const SsiCurve& curve = traced.value().front();
  CHECK_FALSE(curve.closed);
  CHECK(curve.segments.size() >= 2);
  for (std::size_t i = 0; i + 1 < curve.segments.size(); ++i) {
    for (int axis = 0; axis < 3; ++axis) {
      CHECK(std::fabs(curve.segments[i].b[axis] - curve.segments[i + 1].a[axis]) < 1e-12);
      CHECK(curve.segments[i].b[axis] == curve.segments[i + 1].a[axis]);
    }
    CHECK(curve.segments[i].ub == curve.segments[i + 1].ua);
    CHECK(curve.segments[i].vb == curve.segments[i + 1].va);
  }
  for (const SsiSegment& segment : curve.segments) {
    check_segment_certified(a.value(), b.value(), segment, 1e-9);
  }
}

TEST_CASE("ssi: coincident seeds on one circle trace a single closed curve") {
  const Sphere3 sphere{Point3d{0, 0, 0}, Vec3d{0, 0, 1}, Vec3d{1, 0, 0}, 1.0};
  const Result<SsiPatch> a =
      make_ssi_patch(SurfaceGeometry{sphere}, {0, 2 * kTracePi}, {-1.5, 1.5});
  const Result<SsiPatch> b =
      make_ssi_patch(SurfaceGeometry{Plane3{{0, 0, 1}, 0}}, {-2, 2}, {-2, 2});
  REQUIRE(a.has_value());
  REQUIRE(b.has_value());
  const Result<SsiPoint> seed = refine_ssi_point(a.value(), b.value(), {0.1, 0.2, 0.3, 0.4}, 1e-12);
  const Result<SsiPoint> folded =
      refine_ssi_point(a.value(), b.value(), {0.1 + 2 * kTracePi, 0.2, 0.3, 0.4}, 1e-12);
  REQUIRE(seed.has_value());
  REQUIRE(folded.has_value());
  // The two seeds name the same point after the periodic fold, so only one
  // curve is traced. The surviving seed is chosen by the deterministic seed
  // ordering, so swapping the inputs cannot change the output.
  const Result<std::vector<SsiCurve>> traced =
      trace_ssi_curves(a.value(), b.value(), {seed.value(), folded.value()}, 1e-9);
  const Result<std::vector<SsiCurve>> swapped =
      trace_ssi_curves(a.value(), b.value(), {folded.value(), seed.value()}, 1e-9);
  REQUIRE(traced.has_value());
  REQUIRE(swapped.has_value());
  REQUIRE(traced.value().size() == 1);
  REQUIRE(swapped.value().size() == 1);
  CHECK(traced.value().front().closed);
  REQUIRE(traced.value().front().segments.size() == swapped.value().front().segments.size());
  for (std::size_t i = 0; i < traced.value().front().segments.size(); ++i) {
    const SsiSegment& one = traced.value().front().segments[i];
    const SsiSegment& two = swapped.value().front().segments[i];
    CHECK(one.ua == two.ua);
    CHECK(one.ub == two.ub);
    CHECK(one.va == two.va);
    CHECK(one.vb == two.vb);
    CHECK(one.residual == two.residual);
    for (int axis = 0; axis < 3; ++axis) {
      CHECK(one.a[axis] == two.a[axis]);
      CHECK(one.b[axis] == two.b[axis]);
    }
  }
}

TEST_CASE("ssi: two disjoint intersection segments trace two curves") {
  const Cylinder3 cylinder{Point3d{0, 0, 0}, Vec3d{0, 0, 1}, Vec3d{1, 0, 0}, 1.0};
  const Result<SsiPatch> a = make_ssi_patch(SurfaceGeometry{cylinder}, {0, 2 * kTracePi}, {-1, 1});
  const Result<SsiPatch> b =
      make_ssi_patch(SurfaceGeometry{Plane3{{1, 0, 0}, 0.5}}, {-2, 2}, {-2, 2});
  REQUIRE(a.has_value());
  REQUIRE(b.has_value());
  const Result<SsiPoint> right =
      refine_ssi_point(a.value(), b.value(), {kTracePi / 3.0, 0.0, 0.0, 0.8660254037844386}, 1e-12);
  const Result<SsiPoint> left = refine_ssi_point(
      a.value(), b.value(), {-kTracePi / 3.0, 0.0, 0.0, -0.8660254037844386}, 1e-12);
  REQUIRE(right.has_value());
  REQUIRE(left.has_value());
  const Result<std::vector<SsiCurve>> traced =
      trace_ssi_curves(a.value(), b.value(), {right.value(), left.value()}, 1e-9);
  REQUIRE(traced.has_value());
  REQUIRE(traced.value().size() == 2);
  for (const SsiCurve& curve : traced.value()) {
    CHECK_FALSE(curve.closed);
    CHECK(curve.segments.size() >= 2);
    for (const SsiSegment& segment : curve.segments) {
      check_segment_certified(a.value(), b.value(), segment, 1e-9);
    }
    // Each curve is one vertical line: the cylinder angle is constant.
    const double u = curve.segments.front().ua[0];
    CHECK(std::fabs(std::fabs(u) - kTracePi / 3.0) < 1e-9);
    for (const SsiSegment& segment : curve.segments) {
      CHECK(std::fabs(segment.ua[0] - u) < 1e-12);
    }
  }
  // Deterministic output: the curves are ordered by their first vertex.
  const SsiCurve& first_curve = traced.value()[0];
  const SsiCurve& second_curve = traced.value()[1];
  CHECK(first_curve.segments.front().ua < second_curve.segments.front().ua);
}

TEST_CASE("ssi: an unrefinable seed is skipped without corrupting the trace") {
  const Result<SsiPatch> a =
      make_ssi_patch(SurfaceGeometry{Plane3{{0, 0, 1}, 0}}, {-2, 2}, {-2, 2});
  const Result<SsiPatch> b =
      make_ssi_patch(SurfaceGeometry{Plane3{{1, 0, 0}, 0}}, {-2, 2}, {-2, 2});
  REQUIRE(a.has_value());
  REQUIRE(b.has_value());
  const Result<SsiPoint> good = refine_ssi_point(a.value(), b.value(), {0.1, 0.1, 0.1, 0.1}, 1e-12);
  REQUIRE(good.has_value());
  SsiPoint bad;
  bad.ua = {5.0, 5.0};
  bad.va = {5.0, 5.0};
  for (const std::vector<SsiPoint>& seeds :
       {std::vector<SsiPoint>{bad, good.value()}, std::vector<SsiPoint>{good.value(), bad}}) {
    const Result<std::vector<SsiCurve>> traced =
        trace_ssi_curves(a.value(), b.value(), seeds, 1e-9, 64);
    REQUIRE(traced.has_value());
    REQUIRE(traced.value().size() == 1);
    CHECK_FALSE(traced.value().front().closed);
    CHECK(traced.value().front().segments.size() >= 2);
  }
}

TEST_CASE("ssi: a closed trace exhausts a three-step budget") {
  const Sphere3 sphere{Point3d{0, 0, 0}, Vec3d{0, 0, 1}, Vec3d{1, 0, 0}, 1.0};
  const Result<SsiPatch> a =
      make_ssi_patch(SurfaceGeometry{sphere}, {0, 2 * kTracePi}, {-1.5, 1.5});
  const Result<SsiPatch> b =
      make_ssi_patch(SurfaceGeometry{Plane3{{0, 0, 1}, 0}}, {-2, 2}, {-2, 2});
  REQUIRE(a.has_value());
  REQUIRE(b.has_value());
  const Result<SsiPoint> seed = refine_ssi_point(a.value(), b.value(), {0.1, 0.2, 0.3, 0.4}, 1e-12);
  REQUIRE(seed.has_value());
  const Result<std::vector<SsiCurve>> traced =
      trace_ssi_curves(a.value(), b.value(), {seed.value()}, 1e-9, 3);
  REQUIRE_FALSE(traced.has_value());
  CHECK(traced.error().code == Errc::numeric_failure);
  CHECK(std::string_view(traced.error().message) == "trace_ssi_curves: step budget exhausted");
}

TEST_CASE("ssi: trace validates tolerance and step budget") {
  const Result<SsiPatch> a =
      make_ssi_patch(SurfaceGeometry{Plane3{{0, 0, 1}, 0}}, {-2, 2}, {-2, 2});
  const Result<SsiPatch> b =
      make_ssi_patch(SurfaceGeometry{Plane3{{1, 0, 0}, 0}}, {-2, 2}, {-2, 2});
  REQUIRE(a.has_value());
  REQUIRE(b.has_value());
  const Result<SsiPoint> seed = refine_ssi_point(a.value(), b.value(), {0.1, 0.1, 0.1, 0.1}, 1e-12);
  REQUIRE(seed.has_value());
  for (const double tolerance :
       {0.0, -0.0, -1e-9, std::numeric_limits<double>::quiet_NaN(),
        std::numeric_limits<double>::infinity(), -std::numeric_limits<double>::infinity()}) {
    const Result<std::vector<SsiCurve>> traced =
        trace_ssi_curves(a.value(), b.value(), {seed.value()}, tolerance);
    REQUIRE_FALSE(traced.has_value());
    CHECK(traced.error().code == Errc::tolerance_violation);
  }
  // The tolerance outranks a bad budget, matching refine_ssi_point.
  const Result<std::vector<SsiCurve>> ordered =
      trace_ssi_curves(a.value(), b.value(), {seed.value()}, 0.0, 0);
  REQUIRE_FALSE(ordered.has_value());
  CHECK(ordered.error().code == Errc::tolerance_violation);
  for (const int budget : {0, -1, -4096}) {
    const Result<std::vector<SsiCurve>> traced =
        trace_ssi_curves(a.value(), b.value(), {seed.value()}, 1e-9, budget);
    REQUIRE_FALSE(traced.has_value());
    CHECK(traced.error().code == Errc::numeric_failure);
  }
}

TEST_CASE("ssi: trace output is independent of seed order") {
  const Cylinder3 cylinder{Point3d{0, 0, 0}, Vec3d{0, 0, 1}, Vec3d{1, 0, 0}, 1.0};
  const Result<SsiPatch> a = make_ssi_patch(SurfaceGeometry{cylinder}, {0, 2 * kTracePi}, {-1, 1});
  const Result<SsiPatch> b =
      make_ssi_patch(SurfaceGeometry{Plane3{{1, 0, 0}, 0.5}}, {-2, 2}, {-2, 2});
  REQUIRE(a.has_value());
  REQUIRE(b.has_value());
  const Result<SsiPoint> right =
      refine_ssi_point(a.value(), b.value(), {kTracePi / 3.0, 0.0, 0.0, 0.8660254037844386}, 1e-12);
  const Result<SsiPoint> left = refine_ssi_point(
      a.value(), b.value(), {-kTracePi / 3.0, 0.0, 0.0, -0.8660254037844386}, 1e-12);
  REQUIRE(right.has_value());
  REQUIRE(left.has_value());
  const Result<std::vector<SsiCurve>> forward =
      trace_ssi_curves(a.value(), b.value(), {right.value(), left.value()}, 1e-9);
  const Result<std::vector<SsiCurve>> reversed =
      trace_ssi_curves(a.value(), b.value(), {left.value(), right.value()}, 1e-9);
  REQUIRE(forward.has_value());
  REQUIRE(reversed.has_value());
  REQUIRE(forward.value().size() == reversed.value().size());
  for (std::size_t curve_index = 0; curve_index < forward.value().size(); ++curve_index) {
    const SsiCurve& first = forward.value()[curve_index];
    const SsiCurve& second = reversed.value()[curve_index];
    CHECK(first.closed == second.closed);
    REQUIRE(first.segments.size() == second.segments.size());
    for (std::size_t i = 0; i < first.segments.size(); ++i) {
      const SsiSegment& one = first.segments[i];
      const SsiSegment& two = second.segments[i];
      for (int axis = 0; axis < 3; ++axis) {
        CHECK(one.a[axis] == two.a[axis]);
        CHECK(one.b[axis] == two.b[axis]);
      }
      CHECK(one.ua == two.ua);
      CHECK(one.ub == two.ub);
      CHECK(one.va == two.va);
      CHECK(one.vb == two.vb);
      CHECK(one.residual == two.residual);
    }
  }
}
