#include <cmath>
#include <cstring>
#include <vector>

#include <doctest/doctest.h>
#include <geom/brep/ssi.hpp>

using namespace geom::v1;

namespace {
NurbsSurface planar_nurbs() {
  return NurbsSurface(1, 1, {0.0, 0.0, 1.0, 1.0}, {0.0, 0.0, 1.0, 1.0},
                      {Point3d{0, 0, 0}, Point3d{0, 1, 0}, Point3d{1, 0, 0}, Point3d{1, 1, 0}},
                      {1.0, 1.0, 1.0, 1.0});
}

NurbsSurface wavy_nurbs() {
  std::vector<Point3d> control;
  for (int i = 0; i < 4; ++i) {
    for (int j = 0; j < 4; ++j) {
      const double x = static_cast<double>(i) / 3.0;
      const double y = static_cast<double>(j) / 3.0;
      control.push_back(Point3d{x, y, std::sin(3.0 * x) * std::sin(3.0 * y)});
    }
  }
  return NurbsSurface(3, 3, {0, 0, 0, 0, 1, 1, 1, 1}, {0, 0, 0, 0, 1, 1, 1, 1}, control,
                      std::vector<double>(16, 1.0));
}
}  // namespace

TEST_CASE("ssi: tangent and coincident inputs return diagnostics instead of curves") {
  const Sphere3 sphere{Point3d{0, 0, 0}, Vec3d{0, 0, 1}, Vec3d{1, 0, 0}, 1.0};
  const Result<SsiPatch> sphere_patch = make_ssi_patch(
      SurfaceGeometry{sphere}, {0, 2 * 3.14159265358979323846},
      {-0.5 * 3.14159265358979323846, 0.5 * 3.14159265358979323846});
  const Result<SsiPatch> tangent_plane = make_ssi_patch(
      SurfaceGeometry{Plane3{{0, 0, 1}, -1.0}}, {-2, 2}, {-2, 2});
  REQUIRE(sphere_patch.has_value());
  REQUIRE(tangent_plane.has_value());
  const Result<SsiResult> tangent =
      ssi(sphere_patch.value(), tangent_plane.value(), SsiOptions{1e-9});
  REQUIRE(tangent.has_value());
  CHECK(tangent.value().curves.empty());
  REQUIRE(tangent.value().issues.size() == 1);
  CHECK(tangent.value().issues[0].code == SsiDiagnostic::tangent_contact);

  const Result<SsiPatch> plane_a = make_ssi_patch(
      SurfaceGeometry{Plane3{{0, 0, 1}, 2.0}}, {-2, 2}, {-2, 2});
  const Result<SsiPatch> plane_b = make_ssi_patch(
      SurfaceGeometry{Plane3{{0, 0, -1}, -2.0}}, {-2, 2}, {-2, 2});
  REQUIRE(plane_a.has_value());
  REQUIRE(plane_b.has_value());
  const Result<SsiResult> coincident =
      ssi(plane_a.value(), plane_b.value(), SsiOptions{1e-9});
  REQUIRE(coincident.has_value());
  CHECK(coincident.value().curves.empty());
  REQUIRE(coincident.value().issues.size() == 1);
  CHECK(coincident.value().issues[0].code == SsiDiagnostic::coincident);
}

TEST_CASE("ssi: a cylinder patch intersects a plane in two lines") {
  const Cylinder3 cylinder{Point3d{0, 0, 0}, Vec3d{0, 0, 1}, Vec3d{1, 0, 0}, 1.0};
  const Result<SsiPatch> a = make_ssi_patch(
      SurfaceGeometry{cylinder}, {0, 2 * 3.14159265358979323846}, {0, 1});
  const Result<SsiPatch> b = make_ssi_patch(
      SurfaceGeometry{Plane3{{1, 0, 0}, 0.5}}, {-2, 2}, {-2, 2});
  REQUIRE(a.has_value());
  REQUIRE(b.has_value());
  const Result<SsiResult> result = ssi(a.value(), b.value(), SsiOptions{1e-9});
  REQUIRE(result.has_value());
  CHECK(result.value().curves.size() == 2);
  for (const SsiCurve& curve : result.value().curves) {
    for (const SsiSegment& s : curve.segments) {
      CHECK(s.residual <= 1e-9);
      CHECK(std::fabs(s.a[0] - 0.5) < 1e-8);
      CHECK(std::fabs(s.b[0] - 0.5) < 1e-8);
    }
  }
}

TEST_CASE("ssi: a cone patch intersects a plane in a circle") {
  const Cone3 cone{Point3d{0, 0, 0}, Vec3d{0, 0, 1}, Vec3d{1, 0, 0}, 1.0, 0.5};
  const Result<SsiPatch> a = make_ssi_patch(
      SurfaceGeometry{cone}, {0, 2 * 3.14159265358979323846}, {0, 1});
  const Result<SsiPatch> b = make_ssi_patch(
      SurfaceGeometry{Plane3{{0, 0, 1}, 0.5}}, {-2, 2}, {-2, 2});
  REQUIRE(a.has_value());
  REQUIRE(b.has_value());
  const Result<SsiResult> result = ssi(a.value(), b.value(), SsiOptions{1e-9});
  REQUIRE(result.has_value());
  REQUIRE(result.value().curves.size() == 1);
  CHECK(result.value().curves[0].closed);
  const double expected_radius = 1.0 - 0.5 * std::tan(0.5);
  for (const SsiSegment& s : result.value().curves[0].segments) {
    CHECK(s.residual <= 1e-9);
    CHECK(std::fabs(s.a[2] - 0.5) < 1e-8);
    CHECK(std::fabs(std::hypot(s.a[0], s.a[1]) - expected_radius) < 1e-8);
  }
}

TEST_CASE("ssi: parallel offset planes return no curve and no issue") {
  const Result<SsiPatch> a = make_ssi_patch(
      SurfaceGeometry{Plane3{{0, 0, 1}, 0.0}}, {-1, 1}, {-1, 1});
  const Result<SsiPatch> b = make_ssi_patch(
      SurfaceGeometry{Plane3{{0, 0, 1}, 1.0}}, {-1, 1}, {-1, 1});
  REQUIRE(a.has_value());
  REQUIRE(b.has_value());
  const Result<SsiResult> result = ssi(a.value(), b.value(), SsiOptions{1e-9});
  REQUIRE(result.has_value());
  CHECK(result.value().curves.empty());
  CHECK(result.value().issues.empty());
}

TEST_CASE("ssi: a planar NURBS patch intersects a plane in one line") {
  const NurbsSurface patch = planar_nurbs();
  REQUIRE(patch.is_valid());
  const Result<SsiPatch> a = make_ssi_patch(SurfaceGeometry{patch}, {0, 1}, {0, 1});
  const Result<SsiPatch> b = make_ssi_patch(SurfaceGeometry{Plane3{{1, 0, 0}, 0.5}}, {-2, 2}, {-2, 2});
  REQUIRE(a.has_value());
  REQUIRE(b.has_value());
  const Result<SsiResult> result = ssi(a.value(), b.value(), SsiOptions{1e-9});
  REQUIRE(result.has_value());
  REQUIRE(result.value().curves.size() == 1);
  CHECK(result.value().curves[0].segments.size() >= 1);
  for (const SsiSegment& s : result.value().curves[0].segments) {
    CHECK(std::fabs(s.a[0] - 0.5) < 1e-8);
    CHECK(std::fabs(s.b[0] - 0.5) < 1e-8);
  }
}

TEST_CASE("ssi: a full sphere cut by a plane closes across the u seam") {
  const Sphere3 sphere{Point3d{0, 0, 0}, Vec3d{0, 0, 1}, Vec3d{1, 0, 0}, 1.0};
  const Result<SsiPatch> a = make_ssi_patch(SurfaceGeometry{sphere}, {0, 2 * 3.141592653589793}, {-1.5707963267948966, 1.5707963267948966});
  const Result<SsiPatch> b = make_ssi_patch(SurfaceGeometry{Plane3{{0, 0, 1}, 0.2}}, {-2, 2}, {-2, 2});
  REQUIRE(a.has_value());
  REQUIRE(b.has_value());
  const Result<SsiResult> result = ssi(a.value(), b.value(), SsiOptions{1e-9});
  REQUIRE(result.has_value());
  REQUIRE(result.value().curves.size() == 1);
  CHECK(result.value().curves[0].closed);
  const auto& segments = result.value().curves[0].segments;
  REQUIRE(segments.size() >= 8);
  for (const SsiSegment& s : segments) {
    CHECK(std::fabs(s.a[2] - 0.2) < 1e-8);
    CHECK(std::fabs(std::sqrt(s.a[0] * s.a[0] + s.a[1] * s.a[1] + s.a[2] * s.a[2]) - 1.0) < 1e-8);
  }
}

TEST_CASE("ssi: ssi_param_distance wraps periodic axes") {
  const Cylinder3 cylinder{Point3d{0, 0, 0}, Vec3d{0, 0, 1}, Vec3d{1, 0, 0}, 1.0};
  const Result<SsiPatch> patch = make_ssi_patch(SurfaceGeometry{cylinder}, {0, 6.28}, {0, 1});
  REQUIRE(patch.has_value());
  CHECK(ssi_param_distance(patch.value(), {0.01, 0.5}, {6.273185307179586, 0.5}) ==
        doctest::Approx(0.02).epsilon(1e-12));
}

TEST_CASE("ssi: oblique sphere cuts yield one closed curve for every plane extent") {
  const Sphere3 sphere{Point3d{0, 0, 0}, Vec3d{0, 0, 1}, Vec3d{1, 0, 0}, 1.0};
  const Result<SsiPatch> a = make_ssi_patch(SurfaceGeometry{sphere}, {0, 2 * 3.141592653589793}, {-1.5707963267948966, 1.5707963267948966});
  REQUIRE(a.has_value());
  const double inv_root2 = 1.0 / std::sqrt(2.0);
  for (const double offset : {0.5, 0.75}) {
    for (const double extent : {1.5, 2.0, 4.0, 6.0}) {
      INFO("offset=" << offset << " extent=" << extent);
      const Result<SsiPatch> b = make_ssi_patch(
          SurfaceGeometry{Plane3{{inv_root2, 0.0, inv_root2}, offset}}, {-extent, extent}, {-extent, extent});
      REQUIRE(b.has_value());
      const Result<SsiResult> result = ssi(a.value(), b.value(), SsiOptions{1e-9});
      REQUIRE(result.has_value());
      REQUIRE(result.value().curves.size() == 1);
      CHECK(result.value().curves[0].closed);
      for (const SsiSegment& s : result.value().curves[0].segments) {
        CHECK(std::fabs(std::sqrt(s.a[0] * s.a[0] + s.a[1] * s.a[1] + s.a[2] * s.a[2]) - 1.0) < 1e-8);
        CHECK(std::fabs(inv_root2 * (s.a[0] + s.a[2]) - offset) < 1e-8);
        CHECK(s.residual <= 1e-9);
      }
    }
  }
}

TEST_CASE("ssi: a boundary-clipped arc stays open") {
  // The plane cuts the sphere in a full latitude circle; trimming the sphere
  // patch just short of one full turn clips that circle into a single arc whose
  // ends sit on the patch u boundary. The arc is almost closed in space, so a
  // heuristic that closes a trace by comparing its end gap to its overall size
  // would wrongly report it closed (and invent a chord through the patch).
  const Sphere3 sphere{Point3d{0, 0, 0}, Vec3d{0, 0, 1}, Vec3d{1, 0, 0}, 1.0};
  const Result<SsiPatch> a =
      make_ssi_patch(SurfaceGeometry{sphere}, {0, 2 * 3.141592653589793 - 0.4},
                     {-1.5707963267948966, 1.5707963267948966});
  const Result<SsiPatch> b =
      make_ssi_patch(SurfaceGeometry{Plane3{{0, 0, 1}, 0.2}}, {-2, 2}, {-2, 2});
  REQUIRE(a.has_value());
  REQUIRE(b.has_value());
  const Result<SsiResult> result = ssi(a.value(), b.value(), SsiOptions{1e-9});
  REQUIRE(result.has_value());
  REQUIRE(result.value().curves.size() == 1);
  const SsiCurve& curve = result.value().curves[0];
  CHECK_FALSE(curve.closed);
  REQUIRE(curve.segments.size() >= 2);
  // Both ends are trimmed by the patch u boundary, u = 0 or u = 2*pi - 0.4
  // (stored folded into [-pi, pi)).
  const auto on_boundary = [](double u) {
    return std::fabs(u) < 1e-6 || std::fabs(u + 0.4) < 1e-6;
  };
  CHECK(on_boundary(curve.segments.front().ua[0]));
  CHECK(on_boundary(curve.segments.back().ub[0]));
  // The two trimmed ends sit at opposite sides of the removed wedge.
  CHECK(std::fabs(std::remainder(curve.segments.front().ua[0] - curve.segments.back().ub[0],
                                 2 * 3.141592653589793)) > 0.3);
}

TEST_CASE("ssi: a cut next to the sphere pole closes") {
  // The intersection circle at z=0.995 has radius 0.0999. Near the pole a small
  // 3D gap maps to a much larger parameter gap (|du| shrinks with cos(v)), so
  // the closure window has to be wide enough in parameter space to accept the
  // closing step instead of halving it away.
  const Sphere3 sphere{Point3d{0, 0, 0}, Vec3d{0, 0, 1}, Vec3d{1, 0, 0}, 1.0};
  const Result<SsiPatch> a =
      make_ssi_patch(SurfaceGeometry{sphere}, {0, 2 * 3.141592653589793},
                     {-1.5707963267948966, 1.5707963267948966});
  const Result<SsiPatch> b =
      make_ssi_patch(SurfaceGeometry{Plane3{{0, 0, 1}, 0.995}}, {-2, 2}, {-2, 2});
  REQUIRE(a.has_value());
  REQUIRE(b.has_value());
  const Result<SsiResult> result = ssi(a.value(), b.value(), SsiOptions{1e-9});
  REQUIRE(result.has_value());
  REQUIRE(result.value().curves.size() == 1);
  CHECK(result.value().curves[0].closed);
}

TEST_CASE("ssi: torus-plane transverse cuts are argument-order independent") {
  const Torus3 torus{Point3d{0, 0, 0}, Vec3d{0, 0, 1}, Vec3d{1, 0, 0}, 2.0, 0.5};
  const Result<SsiPatch> a = make_ssi_patch(SurfaceGeometry{torus}, {0, 2 * 3.141592653589793},
                                            {0, 2 * 3.141592653589793});
  REQUIRE(a.has_value());
  for (const double z : {0.1, 0.2, 0.25, 0.3, 0.4}) {
    for (const int depth : {6, 8, 12}) {
      INFO("z=" << z << " depth=" << depth);
      const Result<SsiPatch> b = make_ssi_patch(
          SurfaceGeometry{Plane3{{0, 0, 1}, z}}, {-4, 4}, {-4, 4});
      REQUIRE(b.has_value());
      for (const bool swapped : {false, true}) {
        const Result<SsiResult> result = swapped ? ssi(b.value(), a.value(), SsiOptions{1e-9, depth})
                                                 : ssi(a.value(), b.value(), SsiOptions{1e-9, depth});
        REQUIRE(result.has_value());
        REQUIRE(result.value().curves.size() == 2);
        for (const SsiCurve& curve : result.value().curves) {
          CHECK(curve.closed);
          for (const SsiSegment& s : curve.segments) {
            const double radial = std::sqrt(s.a[0] * s.a[0] + s.a[1] * s.a[1]);
            CHECK(std::fabs((radial - 2.0) * (radial - 2.0) + s.a[2] * s.a[2] - 0.25) < 1e-7);
            CHECK(std::fabs(s.a[2] - z) < 1e-8);
            CHECK(s.residual <= 1e-9);
          }
        }
      }
    }
  }
}

TEST_CASE("ssi: a wavy NURBS patch intersects a plane within tolerance") {
  const NurbsSurface patch = wavy_nurbs();
  REQUIRE(patch.is_valid());
  const Result<SsiPatch> a = make_ssi_patch(SurfaceGeometry{patch}, {0, 1}, {0, 1});
  const Result<SsiPatch> b = make_ssi_patch(SurfaceGeometry{Plane3{{0, 0, 1}, 0.25}}, {-2, 2}, {-2, 2});
  REQUIRE(a.has_value());
  REQUIRE(b.has_value());
  const Result<SsiResult> result = ssi(a.value(), b.value(), SsiOptions{1e-8});
  REQUIRE(result.has_value());
  CHECK_FALSE(result.value().curves.empty());
  for (const SsiCurve& curve : result.value().curves) {
    for (const SsiSegment& s : curve.segments) {
      CHECK(s.residual <= 1e-8);
      CHECK(std::fabs(s.a[2] - 0.25) < 1e-8);
      CHECK(std::fabs(s.b[2] - 0.25) < 1e-8);
    }
  }
}

TEST_CASE("ssi: seed components merge periodic parameter branches") {
  const Torus3 torus{Point3d{0, 0, 0}, Vec3d{0, 0, 1}, Vec3d{1, 0, 0}, 2.0, 0.5};
  const double pi = 3.14159265358979323846;
  const Result<SsiPatch> patch =
      make_ssi_patch(SurfaceGeometry{torus}, {0, 2 * pi}, {0, 2 * pi});
  REQUIRE(patch.has_value());
  SsiPoint seam_a;
  seam_a.ua = {-pi + 1e-9, 0.3};
  seam_a.va = {0.4, 0.5};
  seam_a.p = Point3d{2.1, 0.0, 0.4};
  SsiPoint seam_b = seam_a;
  seam_b.ua = {pi - 1e-9, 0.3};
  REQUIRE(ssi_param_distance(patch.value(), seam_a.ua, seam_b.ua) < 1e-6);
  const std::vector<SsiPoint> representatives = detail::ssi_driver_curve_seeds(
      patch.value(), patch.value(), {seam_a, seam_b}, {}, 1e-6);
  CHECK(representatives.size() == 1);
}

TEST_CASE("ssi: a segment refined at only one end keeps its refined seed") {
  const Sphere3 sphere{Point3d{0, 0, 0}, Vec3d{0, 0, 1}, Vec3d{1, 0, 0}, 1.0};
  const Result<SsiPatch> patch = make_ssi_patch(
      SurfaceGeometry{sphere}, {0, 2 * 3.14159265358979323846}, {-1.5, 1.5});
  REQUIRE(patch.has_value());
  SsiPoint refined;
  refined.ua = {0.2, 0.1};
  refined.va = {0.3, 0.4};
  refined.p = Point3d{0.9, 0.1, 0.2};
  SsiPoint unreferenced = refined;
  unreferenced.ua = {1.0, 0.9};
  unreferenced.va = {-0.4, 0.7};
  unreferenced.p = Point3d{-0.9, 0.2, -0.1};
  const std::vector<std::array<std::size_t, 2>> segments = {
      {0, detail::kSsiDriverNoSeed}};
  const std::vector<SsiPoint> representatives = detail::ssi_driver_curve_seeds(
      patch.value(), patch.value(), {refined, unreferenced}, segments, 1e-6);
  REQUIRE(representatives.size() == 2);
  bool kept_refined = false;
  for (const SsiPoint& candidate : representatives) {
    if (candidate.ua == refined.ua && candidate.va == refined.va) kept_refined = true;
  }
  CHECK(kept_refined);
}

TEST_CASE("ssi: an unreferenced seed is its own component") {
  const Sphere3 sphere{Point3d{0, 0, 0}, Vec3d{0, 0, 1}, Vec3d{1, 0, 0}, 1.0};
  const Result<SsiPatch> patch = make_ssi_patch(
      SurfaceGeometry{sphere}, {0, 2 * 3.14159265358979323846}, {-1.5, 1.5});
  REQUIRE(patch.has_value());
  SsiPoint one;
  one.ua = {0.1, 0.2};
  one.va = {0.3, 0.4};
  one.p = Point3d{0.8, 0.1, 0.2};
  SsiPoint two = one;
  two.ua = {-1.0, -0.5};
  two.va = {-0.3, 0.1};
  two.p = Point3d{-0.2, 0.7, 0.3};
  const std::vector<SsiPoint> representatives = detail::ssi_driver_curve_seeds(
      patch.value(), patch.value(), {one, two}, {}, 1e-6);
  CHECK(representatives.size() == 2);
}

TEST_CASE("ssi: max_depth below the clustering floor is rejected") {
  // The seed clustering assumes at least 64 cells per patch side; coarser
  // meshes produce duplicate and argument-order-dependent curves, so the
  // option is rejected before any tessellation or predicate work.
  const Torus3 torus{Point3d{0, 0, 0}, Vec3d{0, 0, 1}, Vec3d{1, 0, 0}, 2.0, 0.5};
  const Result<SsiPatch> a = make_ssi_patch(SurfaceGeometry{torus}, {0, 2 * 3.141592653589793},
                                            {0, 2 * 3.141592653589793});
  const Result<SsiPatch> b =
      make_ssi_patch(SurfaceGeometry{Plane3{{0, 0, 1}, 0.1}}, {-4, 4}, {-4, 4});
  REQUIRE(a.has_value());
  REQUIRE(b.has_value());
  for (const int depth : {1, 3, 5}) {
    INFO("depth=" << depth);
    const Result<SsiResult> rejected = ssi(a.value(), b.value(), SsiOptions{1e-9, depth});
    REQUIRE_FALSE(rejected.has_value());
    CHECK(rejected.error().code == Errc::unsupported);
    CHECK(std::strcmp(rejected.error().message, "ssi: max_depth must be at least 6") == 0);
  }
  const Result<SsiResult> accepted = ssi(a.value(), b.value(), SsiOptions{1e-9, 6});
  REQUIRE(accepted.has_value());
  CHECK(accepted.value().curves.size() == 2);
}

TEST_CASE("ssi: a sphere cut by the equator plane closes once") {
  // At depth 4 the polar triangles make the refined seeds split a single
  // equator circle into seven coincident closed curves; this is a sanity
  // regression for the depth-6 floor (the RED case is depth 4 on the old HEAD).
  const Sphere3 sphere{Point3d{0, 0, 0}, Vec3d{0, 0, 1}, Vec3d{1, 0, 0}, 1.0};
  const Result<SsiPatch> a = make_ssi_patch(SurfaceGeometry{sphere}, {0, 2 * 3.141592653589793},
                                            {-1.5707963267948966, 1.5707963267948966});
  const Result<SsiPatch> b =
      make_ssi_patch(SurfaceGeometry{Plane3{{0, 0, 1}, 0.0}}, {-2, 2}, {-2, 2});
  REQUIRE(a.has_value());
  REQUIRE(b.has_value());
  const Result<SsiResult> result = ssi(a.value(), b.value(), SsiOptions{1e-9, 6});
  REQUIRE(result.has_value());
  CHECK(result.value().curves.size() == 1);
  CHECK(result.value().curves.front().closed);
}

TEST_CASE("ssi: near-crest torus cuts keep two branches apart") {
  // The two intersection circles at z=0.45..0.499 are only 0.063 apart near the
  // crest, while a bridged duplicate seed pair on one circle is up to 0.0485
  // apart (92.5% of the 0.0525 union radius). The radius must cover the
  // duplicates without bridging the branches.
  const Torus3 torus{Point3d{0, 0, 0}, Vec3d{0, 0, 1}, Vec3d{1, 0, 0}, 2.0, 0.5};
  const Result<SsiPatch> a = make_ssi_patch(SurfaceGeometry{torus}, {0, 2 * 3.141592653589793},
                                            {0, 2 * 3.141592653589793});
  REQUIRE(a.has_value());
  for (const double z : {0.45, 0.48, 0.49, 0.499}) {
    for (const int depth : {6, 12}) {
      INFO("z=" << z << " depth=" << depth);
      const Result<SsiPatch> b = make_ssi_patch(
          SurfaceGeometry{Plane3{{0, 0, 1}, z}}, {-4, 4}, {-4, 4});
      REQUIRE(b.has_value());
      for (const bool swapped : {false, true}) {
        const Result<SsiResult> result = swapped ? ssi(b.value(), a.value(), SsiOptions{1e-9, depth})
                                                 : ssi(a.value(), b.value(), SsiOptions{1e-9, depth});
        REQUIRE(result.has_value());
        const auto& curves = result.value().curves;
        REQUIRE(curves.size() == 2);
        CHECK(curves[0].closed);
        CHECK(curves[1].closed);
        double radius[2] = {0.0, 0.0};
        for (std::size_t i = 0; i < curves.size(); ++i) {
          double sum = 0.0;
          for (const auto& segment : curves[i].segments) {
            sum += std::hypot(segment.a[0], segment.a[1]);
            sum += std::hypot(segment.b[0], segment.b[1]);
          }
          radius[i] = sum / (2.0 * static_cast<double>(curves[i].segments.size()));
        }
        const double offset = std::sqrt(0.25 - z * z);
        const double lo = radius[0] < radius[1] ? radius[0] : radius[1];
        const double hi = radius[0] < radius[1] ? radius[1] : radius[0];
        CHECK(std::fabs(lo - (2.0 - offset)) < 1e-3);
        CHECK(std::fabs(hi - (2.0 + offset)) < 1e-3);
      }
    }
  }
  // z=0.4995 puts the branches 0.0447 apart, below the 0.0525 union radius, so
  // the union necessarily bridges them; the measured onset is z=0.49932
  // (separation < radius). That near-crest band is a documented M3b-1 limit
  // rather than a regression; z=0.499 has a 1.205x margin and is pinned.
}

TEST_CASE("ssi: a wide plane patch does not merge the torus branches") {
  // The two torus-plane circles at z=0.4 have radii 1.7 and 2.3, so they are
  // about 0.6 apart. A seed radius scaled by the plane patch's bounding box
  // grows with the plane extent and merges them; the union must instead be a
  // 3D distance bounded by the tessellation error.
  const Torus3 torus{Point3d{0, 0, 0}, Vec3d{0, 0, 1}, Vec3d{1, 0, 0}, 2.0, 0.5};
  const Result<SsiPatch> a = make_ssi_patch(SurfaceGeometry{torus}, {0, 2 * 3.141592653589793},
                                            {0, 2 * 3.141592653589793});
  REQUIRE(a.has_value());
  for (const double extent : {4.0, 12.0}) {
    INFO("extent=" << extent);
    const Result<SsiPatch> b = make_ssi_patch(
        SurfaceGeometry{Plane3{{0, 0, 1}, 0.4}}, {-extent, extent}, {-extent, extent});
    REQUIRE(b.has_value());
    const Result<SsiResult> result = ssi(a.value(), b.value(), SsiOptions{1e-9});
    REQUIRE(result.has_value());
    REQUIRE(result.value().curves.size() == 2);
    for (const SsiCurve& curve : result.value().curves) {
      CHECK(curve.closed);
      double length = 0.0;
      for (const SsiSegment& s : curve.segments) length += norm2(s.b - s.a);
      CHECK(length >= 2 * 3.141592653589793 * 1.69);
      CHECK(length <= 2 * 3.141592653589793 * 2.31);
    }
  }
}

TEST_CASE("ssi: a trace that cannot close or reach a boundary is rejected") {
  // The intersection circle leaves the v band at v=1.5028, but the sphere patch
  // ends at v=1.55, so the circle is entirely inside and must close. The walk
  // instead keeps circling; the lap guard has to stop it instead of returning a
  // seven-lap arc.
  const Sphere3 sphere{Point3d{0, 0, 0}, Vec3d{0, 0, 1}, Vec3d{1, 0, 0}, 1.0};
  const Result<SsiPatch> a = make_ssi_patch(SurfaceGeometry{sphere}, {0, 2 * 3.141592653589793},
                                            {-1.55, 1.55});
  const double inv_root2 = 1.0 / std::sqrt(2.0);
  const Result<SsiPatch> b = make_ssi_patch(
      SurfaceGeometry{Plane3{{inv_root2, 0.0, inv_root2}, 0.75}}, {-2, 2}, {-2, 2});
  REQUIRE(a.has_value());
  REQUIRE(b.has_value());
  const Result<SsiResult> result = ssi(a.value(), b.value(), SsiOptions{1e-9});
  if (!result.has_value()) {
    CHECK(result.error().code == Errc::numeric_failure);
    CHECK(std::strcmp(result.error().message,
                      "trace_ssi_curves: curve does not close or reach a boundary") == 0);
  } else {
    REQUIRE(result.value().curves.size() == 1);
    CHECK_FALSE(result.value().curves[0].closed);
    double length = 0.0;
    for (const SsiSegment& s : result.value().curves[0].segments) length += norm2(s.b - s.a);
    CHECK(length <= 6.0);
  }
}

TEST_CASE("ssi: a boundary-clipped arc ends on the v boundary") {
  // The circle leaves the v band at v=1.5028, so both trimmed ends must sit on
  // |v|=1.5 within the boundary tolerance rather than a fraction of a radian
  // short of it.
  const Sphere3 sphere{Point3d{0, 0, 0}, Vec3d{0, 0, 1}, Vec3d{1, 0, 0}, 1.0};
  const Result<SsiPatch> a = make_ssi_patch(SurfaceGeometry{sphere}, {0, 2 * 3.141592653589793},
                                            {-1.5, 1.5});
  const double inv_root2 = 1.0 / std::sqrt(2.0);
  const Result<SsiPatch> b = make_ssi_patch(
      SurfaceGeometry{Plane3{{inv_root2, 0.0, inv_root2}, 0.75}}, {-2, 2}, {-2, 2});
  REQUIRE(a.has_value());
  REQUIRE(b.has_value());
  const Result<SsiResult> result = ssi(a.value(), b.value(), SsiOptions{1e-9});
  REQUIRE(result.has_value());
  REQUIRE(result.value().curves.size() == 1);
  const SsiCurve& curve = result.value().curves[0];
  CHECK_FALSE(curve.closed);
  REQUIRE(curve.segments.size() >= 2);
  CHECK(std::fabs(std::fabs(curve.segments.front().ua[1]) - 1.5) <= 1e-6);
  CHECK(std::fabs(std::fabs(curve.segments.back().ub[1]) - 1.5) <= 1e-6);
}

TEST_CASE("ssi: crossing branches keep separate seed components") {
  // A widened seed radius can reach a different branch that passes close by.
  // Two seeds at the same location and parameters but leaving along different
  // directions must not be merged; the tangent guard keeps them apart.
  const Sphere3 sphere{Point3d{0, 0, 0}, Vec3d{0, 0, 1}, Vec3d{1, 0, 0}, 1.0};
  const Result<SsiPatch> patch = make_ssi_patch(
      SurfaceGeometry{sphere}, {0, 2 * 3.14159265358979323846},
      {-1.5707963267948966, 1.5707963267948966});
  REQUIRE(patch.has_value());
  SsiPoint first;
  first.ua = {0.2, 0.1};
  first.va = {0.3, 0.4};
  first.p = Point3d{0.9, 0.1, 0.2};
  SsiPoint second = first;
  SsiPoint third = first;
  third.p = Point3d{1.2, 0.1, 0.2};
  SsiPoint fourth = first;
  fourth.p = Point3d{0.9, 0.1, 0.5};
  const std::vector<std::array<std::size_t, 2>> segments = {{0, 2}, {1, 3}};
  const std::vector<SsiPoint> representatives = detail::ssi_driver_curve_seeds(
      patch.value(), patch.value(), {first, second, third, fourth}, segments, 1e-6);
  CHECK(representatives.size() == 2);
}
