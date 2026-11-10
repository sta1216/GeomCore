#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstring>
#include <limits>
#include <vector>

#include <doctest/doctest.h>
#include <geom/brep/ssi.hpp>

using namespace geom::v1;

namespace {
// Counts the triangle-pair intersection segments exactly as ssi() extracts them,
// so the statistics tests can compare accepted_segments + rejected_segments
// against an independently measured extraction count.
std::uint64_t extracted_segments(const SsiPatch& a, const SsiPatch& b, double tolerance,
                                 int max_depth) {
  const Box3d bounds_a = patch_bounds(a, tolerance);
  const Box3d bounds_b = patch_bounds(b, tolerance);
  Box3d bounds = bounds_a;
  bounds.expand(bounds_b);
  const double chord_tolerance = std::max(1e-3 * norm2(bounds.diagonal()), 64.0 * tolerance);
  const Result<SsiPatchMesh> mesh_a = tessellate_ssi_patch(a, chord_tolerance, max_depth);
  const Result<SsiPatchMesh> mesh_b = tessellate_ssi_patch(b, chord_tolerance, max_depth);
  if (!mesh_a.has_value() || !mesh_b.has_value()) return 0;
  const SsiPatchMesh& patch_mesh_a = mesh_a.value();
  const SsiPatchMesh& patch_mesh_b = mesh_b.value();
  std::vector<Box3d> boxes_a;
  boxes_a.reserve(patch_mesh_a.triangles.size());
  for (const std::array<std::uint32_t, 3>& triangle : patch_mesh_a.triangles) {
    boxes_a.push_back(detail::ssi_driver_triangle_bounds(patch_mesh_a, triangle));
  }
  std::vector<Box3d> boxes_b;
  boxes_b.reserve(patch_mesh_b.triangles.size());
  for (const std::array<std::uint32_t, 3>& triangle : patch_mesh_b.triangles) {
    boxes_b.push_back(detail::ssi_driver_triangle_bounds(patch_mesh_b, triangle));
  }
  std::uint64_t count = 0;
  for (std::size_t index_a = 0; index_a < patch_mesh_a.triangles.size(); ++index_a) {
    const std::array<std::uint32_t, 3>& triangle_a = patch_mesh_a.triangles[index_a];
    for (std::size_t index_b = 0; index_b < patch_mesh_b.triangles.size(); ++index_b) {
      if (!detail::ssi_driver_bounds_overlap(boxes_a[index_a], boxes_b[index_b])) continue;
      const std::array<std::uint32_t, 3>& triangle_b = patch_mesh_b.triangles[index_b];
      const TriTriResult hit = tri_tri_intersect(
          patch_mesh_a.points[triangle_a[0]], patch_mesh_a.points[triangle_a[1]],
          patch_mesh_a.points[triangle_a[2]], patch_mesh_b.points[triangle_b[0]],
          patch_mesh_b.points[triangle_b[1]], patch_mesh_b.points[triangle_b[2]], tolerance);
      if (hit.kind == TriTriKind::segment) ++count;
    }
  }
  return count;
}

// Mirrors the two-ulp floor computed by ssi() from the patch bounds.
double coordinate_floor(const SsiPatch& a, const SsiPatch& b) {
  const Box3d bounds_a = patch_bounds(a, 1e-9);
  const Box3d bounds_b = patch_bounds(b, 1e-9);
  double max_abs = 0.0;
  for (const Box3d* box : {&bounds_a, &bounds_b}) {
    for (int axis = 0; axis < 3; ++axis) {
      max_abs = std::max(max_abs, std::fabs(box->lo[axis]));
      max_abs = std::max(max_abs, std::fabs(box->hi[axis]));
    }
  }
  return 2.0 * std::numeric_limits<double>::epsilon() * std::max(1.0, max_abs);
}

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

TEST_CASE("ssi: validates options and patches") {
  const Result<SsiPatch> a = make_ssi_patch(SurfaceGeometry{Plane3{}}, {-1, 1}, {-1, 1});
  REQUIRE(a.has_value());
  const auto expect_error = [&a](const SsiOptions& options, Errc code, const char* message) {
    const Result<SsiResult> result = ssi(a.value(), a.value(), options);
    REQUIRE_FALSE(result.has_value());
    CHECK(result.error().code == code);
    CHECK(std::strcmp(result.error().message, message) == 0);
  };
  expect_error(SsiOptions{0.0}, Errc::numeric_failure,
               "ssi: tolerance must be positive and finite");
  expect_error(SsiOptions{-1.0}, Errc::numeric_failure,
               "ssi: tolerance must be positive and finite");
  expect_error(SsiOptions{std::nan("")}, Errc::numeric_failure,
               "ssi: tolerance must be positive and finite");
  expect_error(SsiOptions{std::numeric_limits<double>::infinity()}, Errc::numeric_failure,
               "ssi: tolerance must be positive and finite");
  SsiOptions bad_depth{1e-9};
  bad_depth.max_depth = 5;
  expect_error(bad_depth, Errc::unsupported, "ssi: max_depth must be at least 6");
  SsiOptions no_segments{1e-9};
  no_segments.max_segments = 0;
  expect_error(no_segments, Errc::numeric_failure,
               "ssi: max_segments must be positive");
  // Tolerance is validated before the other options, so a request that breaks
  // all three reports the tolerance first.
  SsiOptions broken{0.0};
  broken.max_depth = 0;
  broken.max_segments = 0;
  expect_error(broken, Errc::numeric_failure,
               "ssi: tolerance must be positive and finite");
  // Depth is validated before the segment budget.
  SsiOptions broken_depth{1e-9};
  broken_depth.max_depth = 0;
  broken_depth.max_segments = 0;
  expect_error(broken_depth, Errc::unsupported, "ssi: max_depth must be at least 6");
}

TEST_CASE("ssi: identical inputs produce identical output") {
  const Sphere3 sphere{Point3d{0, 0, 0}, Vec3d{0, 0, 1}, Vec3d{1, 0, 0}, 1.0};
  const Result<SsiPatch> a = make_ssi_patch(
      SurfaceGeometry{sphere}, {0, 2 * 3.141592653589793}, {-1.5, 1.5});
  const Result<SsiPatch> b = make_ssi_patch(
      SurfaceGeometry{Plane3{{0, 0, 1}, 0.1}}, {-2, 2}, {-2, 2});
  REQUIRE(a.has_value());
  REQUIRE(b.has_value());
  const Result<SsiResult> r1 = ssi(a.value(), b.value(), SsiOptions{1e-9});
  const Result<SsiResult> r2 = ssi(a.value(), b.value(), SsiOptions{1e-9});
  REQUIRE(r1.has_value());
  REQUIRE(r2.has_value());
  REQUIRE(r1.value().curves.size() == r2.value().curves.size());
  for (std::size_t i = 0; i < r1.value().curves.size(); ++i) {
    CHECK(r1.value().curves[i].closed == r2.value().curves[i].closed);
    const auto& s1 = r1.value().curves[i].segments;
    const auto& s2 = r2.value().curves[i].segments;
    REQUIRE(s1.size() == s2.size());
    for (std::size_t k = 0; k < s1.size(); ++k) {
      for (int c = 0; c < 3; ++c) {
        CHECK(s1[k].a[c] == s2[k].a[c]);
        CHECK(s1[k].b[c] == s2[k].b[c]);
      }
      CHECK(s1[k].ua == s2[k].ua);
      CHECK(s1[k].ub == s2[k].ub);
      CHECK(s1[k].va == s2[k].va);
      CHECK(s1[k].vb == s2[k].vb);
      CHECK(s1[k].residual == s2[k].residual);
    }
  }
  REQUIRE(r1.value().issues.size() == r2.value().issues.size());
  for (std::size_t i = 0; i < r1.value().issues.size(); ++i) {
    CHECK(r1.value().issues[i].code == r2.value().issues[i].code);
    for (int p = 0; p < 4; ++p) {
      CHECK(r1.value().issues[i].params[p] == r2.value().issues[i].params[p]);
    }
  }
  CHECK(r1.value().candidate_pairs == r2.value().candidate_pairs);
  CHECK(r1.value().accepted_segments == r2.value().accepted_segments);
  CHECK(r1.value().rejected_segments == r2.value().rejected_segments);
  CHECK(r1.value().max_residual == r2.value().max_residual);
}

TEST_CASE("ssi: a segment budget is rejected instead of truncating") {
  const Cylinder3 cylinder{Point3d{0, 0, 0}, Vec3d{0, 0, 1}, Vec3d{1, 0, 0}, 1.0};
  const Result<SsiPatch> a = make_ssi_patch(
      SurfaceGeometry{cylinder}, {0, 2 * 3.141592653589793}, {0, 1});
  const Result<SsiPatch> b = make_ssi_patch(
      SurfaceGeometry{Plane3{{1, 0, 0}, 0.5}}, {-2, 2}, {-2, 2});
  REQUIRE(a.has_value());
  REQUIRE(b.has_value());
  const Result<SsiResult> reference = ssi(a.value(), b.value(), SsiOptions{1e-9});
  REQUIRE(reference.has_value());
  REQUIRE(reference.value().curves.size() == 2);
  const std::uint64_t accepted = reference.value().accepted_segments;
  REQUIRE(accepted > 1);

  SsiOptions exact{1e-9};
  exact.max_segments = accepted;
  const Result<SsiResult> exact_result = ssi(a.value(), b.value(), exact);
  REQUIRE(exact_result.has_value());
  CHECK(exact_result.value().accepted_segments == accepted);
  CHECK(exact_result.value().curves.size() == 2);

  SsiOptions tight{1e-9};
  tight.max_segments = accepted - 1;
  const Result<SsiResult> rejected = ssi(a.value(), b.value(), tight);
  REQUIRE_FALSE(rejected.has_value());
  CHECK(rejected.error().code == Errc::numeric_failure);
  CHECK(std::strcmp(rejected.error().message, "ssi: segment budget exhausted") == 0);
}

TEST_CASE("ssi: a collapsed patch is diagnosed as degenerate") {
  const Sphere3 sphere{Point3d{0, 0, 0}, Vec3d{0, 0, 1}, Vec3d{1, 0, 0}, 1.0};
  const Result<SsiPatch> tiny = make_ssi_patch(
      SurfaceGeometry{sphere}, {0, 1e-300}, {0, 1e-300});
  const Result<SsiPatch> plane = make_ssi_patch(
      SurfaceGeometry{Plane3{{0, 0, 1}, 1.0}}, {-2, 2}, {-2, 2});
  REQUIRE(tiny.has_value());
  REQUIRE(plane.has_value());
  const Result<SsiResult> result = ssi(tiny.value(), plane.value(), SsiOptions{1e-9});
  REQUIRE(result.has_value());
  CHECK(result.value().curves.empty());
  REQUIRE(result.value().issues.size() == 1);
  CHECK(result.value().issues[0].code == SsiDiagnostic::degenerate_patch);
  CHECK(result.value().candidate_pairs == 0);
  CHECK(result.value().accepted_segments == 0);
  CHECK(result.value().rejected_segments == 0);
  CHECK(result.value().max_residual == 0.0);
}

TEST_CASE("ssi: an unreachable tolerance is reported as no_convergence") {
  const Sphere3 sphere{Point3d{0, 0, 0}, Vec3d{0, 0, 1}, Vec3d{1, 0, 0}, 1.0};
  const Result<SsiPatch> a = make_ssi_patch(
      SurfaceGeometry{sphere}, {0, 2 * 3.14159265358979323846}, {-1.5, 1.5});
  const Result<SsiPatch> b = make_ssi_patch(
      SurfaceGeometry{Plane3{{0, 0, 1}, 0.1}}, {-2, 2}, {-2, 2});
  REQUIRE(a.has_value());
  REQUIRE(b.has_value());
  const Result<SsiResult> result = ssi(a.value(), b.value(), SsiOptions{1e-300});
  REQUIRE(result.has_value());
  CHECK(result.value().curves.empty());
  REQUIRE(result.value().issues.size() == 1);
  CHECK(result.value().issues[0].code == SsiDiagnostic::no_convergence);
  CHECK(result.value().candidate_pairs == 0);
  CHECK(result.value().accepted_segments == 0);
  CHECK(result.value().rejected_segments == 0);
  CHECK(result.value().max_residual == 0.0);
}

TEST_CASE("ssi: large coordinates trace at an achievable tolerance") {
  const Sphere3 sphere{Point3d{0, 0, 0}, Vec3d{0, 0, 1}, Vec3d{1, 0, 0}, 1.0e6};
  const Result<SsiPatch> a = make_ssi_patch(
      SurfaceGeometry{sphere}, {0, 2 * 3.14159265358979323846}, {-1.5, 1.5});
  const Result<SsiPatch> b = make_ssi_patch(
      SurfaceGeometry{Plane3{{0, 0, 1}, 5.0e5}}, {-2.0e6, 2.0e6}, {-2.0e6, 2.0e6});
  REQUIRE(a.has_value());
  REQUIRE(b.has_value());
  const Result<SsiResult> result = ssi(a.value(), b.value(), SsiOptions{1e-9});
  REQUIRE(result.has_value());
  CHECK(result.value().issues.empty());
  REQUIRE(result.value().curves.size() == 1);
  CHECK(result.value().curves[0].closed);
  CHECK(result.value().max_residual <= 1e-9);
}

TEST_CASE("ssi: a tight but achievable tolerance still traces") {
  const Sphere3 sphere{Point3d{0, 0, 0}, Vec3d{0, 0, 1}, Vec3d{1, 0, 0}, 1.0};
  const Result<SsiPatch> a = make_ssi_patch(
      SurfaceGeometry{sphere}, {0, 2 * 3.14159265358979323846}, {-1.5, 1.5});
  const Result<SsiPatch> b = make_ssi_patch(
      SurfaceGeometry{Plane3{{0, 0, 1}, 0.1}}, {-2, 2}, {-2, 2});
  REQUIRE(a.has_value());
  REQUIRE(b.has_value());
  const Result<SsiResult> result = ssi(a.value(), b.value(), SsiOptions{1e-14});
  REQUIRE(result.has_value());
  CHECK_FALSE(result.value().curves.empty());
  CHECK(result.value().max_residual <= 1e-14);
}

TEST_CASE("ssi: the tolerance floor is two ulp of the coordinate magnitude") {
  const Sphere3 sphere{Point3d{0, 0, 0}, Vec3d{0, 0, 1}, Vec3d{1, 0, 0}, 1.0};
  const Result<SsiPatch> a = make_ssi_patch(
      SurfaceGeometry{sphere}, {0, 2 * 3.14159265358979323846}, {-1.5, 1.5});
  const Result<SsiPatch> b = make_ssi_patch(
      SurfaceGeometry{Plane3{{0, 0, 1}, 0.1}}, {-2, 2}, {-2, 2});
  REQUIRE(a.has_value());
  REQUIRE(b.has_value());
  const double floor = coordinate_floor(a.value(), b.value());
  const double below = std::nextafter(floor, -std::numeric_limits<double>::infinity());
  const Result<SsiResult> gated = ssi(a.value(), b.value(), SsiOptions{below});
  REQUIRE(gated.has_value());
  CHECK(gated.value().curves.empty());
  REQUIRE(gated.value().issues.size() == 1);
  CHECK(gated.value().issues[0].code == SsiDiagnostic::no_convergence);
  const Result<SsiResult> at_floor = ssi(a.value(), b.value(), SsiOptions{floor});
  REQUIRE(at_floor.has_value());
  REQUIRE(at_floor.value().curves.size() == 1);
  CHECK(at_floor.value().curves[0].closed);
  for (const SsiIssue& issue : at_floor.value().issues) {
    CHECK(issue.code != SsiDiagnostic::no_convergence);
  }
}

TEST_CASE("ssi: rejected segments report the saturated extracted-vs-accepted surplus") {
  // Perpendicular planes: every extracted segment belongs to the returned line.
  const Result<SsiPatch> plane_z = make_ssi_patch(
      SurfaceGeometry{Plane3{{0, 0, 1}, 0}}, {-2, 2}, {-2, 2});
  const Result<SsiPatch> plane_x = make_ssi_patch(
      SurfaceGeometry{Plane3{{1, 0, 0}, 0}}, {-2, 2}, {-2, 2});
  REQUIRE(plane_z.has_value());
  REQUIRE(plane_x.has_value());
  const Result<SsiResult> planes = ssi(plane_z.value(), plane_x.value(), SsiOptions{1e-9});
  REQUIRE(planes.has_value());
  REQUIRE(planes.value().curves.size() == 1);
  const std::uint64_t plane_extracted =
      extracted_segments(plane_z.value(), plane_x.value(), 1e-9, 12);
  CHECK(plane_extracted == 2);
  CHECK(planes.value().accepted_segments == 21);
  // A coarse planar tessellation yields fewer extracted segments than trace
  // steps (accepted can exceed extracted), so the difference saturates at zero.
  CHECK(planes.value().rejected_segments == 0);

  // Sphere/plane: accepted 31 of the 96 extracted segments; the saturated
  // surplus is 65 (trace segmentation is finer than seed segmentation).
  const Sphere3 sphere{Point3d{0, 0, 0}, Vec3d{0, 0, 1}, Vec3d{1, 0, 0}, 1.0};
  const Result<SsiPatch> sphere_patch = make_ssi_patch(
      SurfaceGeometry{sphere}, {0, 2 * 3.14159265358979323846}, {-1.5, 1.5});
  const Result<SsiPatch> plane_z01 = make_ssi_patch(
      SurfaceGeometry{Plane3{{0, 0, 1}, 0.1}}, {-2, 2}, {-2, 2});
  REQUIRE(sphere_patch.has_value());
  REQUIRE(plane_z01.has_value());
  const Result<SsiResult> sphere_result =
      ssi(sphere_patch.value(), plane_z01.value(), SsiOptions{1e-9});
  REQUIRE(sphere_result.has_value());
  REQUIRE(sphere_result.value().curves.size() == 1);
  const std::uint64_t sphere_extracted =
      extracted_segments(sphere_patch.value(), plane_z01.value(), 1e-9, 12);
  CHECK(sphere_extracted == 96);
  CHECK(sphere_result.value().accepted_segments == 31);
  CHECK(sphere_result.value().rejected_segments == 65);
  CHECK(sphere_result.value().accepted_segments + sphere_result.value().rejected_segments ==
        sphere_extracted);
  double worst_residual = 0.0;
  for (const SsiSegment& segment : sphere_result.value().curves[0].segments) {
    worst_residual = std::max(worst_residual, segment.residual);
  }
  CHECK(sphere_result.value().max_residual == worst_residual);
  CHECK(sphere_result.value().max_residual > 0.0);
  CHECK(sphere_result.value().max_residual <= 1e-9);

  // Cylinder/plane: two open lines; the saturated surplus is reported.
  const Cylinder3 cylinder{Point3d{0, 0, 0}, Vec3d{0, 0, 1}, Vec3d{1, 0, 0}, 1.0};
  const Result<SsiPatch> cylinder_patch = make_ssi_patch(
      SurfaceGeometry{cylinder}, {0, 2 * 3.141592653589793}, {0, 1});
  REQUIRE(cylinder_patch.has_value());
  const Result<SsiResult> cylinder_result =
      ssi(cylinder_patch.value(), plane_x.value(), SsiOptions{1e-9});
  REQUIRE(cylinder_result.has_value());
  REQUIRE(cylinder_result.value().curves.size() == 2);
  const std::uint64_t cylinder_extracted =
      extracted_segments(cylinder_patch.value(), plane_x.value(), 1e-9, 12);
  CHECK(cylinder_result.value().accepted_segments > 0);
  CHECK(cylinder_result.value().accepted_segments + cylinder_result.value().rejected_segments ==
        cylinder_extracted);
  CHECK(cylinder_result.value().rejected_segments < cylinder_extracted);
}

TEST_CASE("ssi: curves are ordered by their first point") {
  const Torus3 torus{Point3d{0, 0, 0}, Vec3d{0, 0, 1}, Vec3d{1, 0, 0}, 2.0, 0.5};
  const Result<SsiPatch> torus_patch = make_ssi_patch(
      SurfaceGeometry{torus}, {0, 2 * 3.141592653589793}, {0, 2 * 3.141592653589793});
  const Result<SsiPatch> torus_plane = make_ssi_patch(
      SurfaceGeometry{Plane3{{0, 0, 1}, 0.4}}, {-4, 4}, {-4, 4});
  REQUIRE(torus_patch.has_value());
  REQUIRE(torus_plane.has_value());
  const Result<SsiResult> torus_result = ssi(torus_patch.value(), torus_plane.value(), SsiOptions{1e-9});
  REQUIRE(torus_result.has_value());
  REQUIRE(torus_result.value().curves.size() == 2);
  const Point3d& first = torus_result.value().curves[0].segments.front().a;
  const Point3d& second = torus_result.value().curves[1].segments.front().a;
  CHECK(first[0] < second[0]);
  CHECK(std::fabs(first[0] + 2.3) < 1e-6);
  CHECK(std::fabs(second[0] + 1.7) < 1e-6);

  const Cylinder3 cylinder{Point3d{0, 0, 0}, Vec3d{0, 0, 1}, Vec3d{1, 0, 0}, 1.0};
  const Result<SsiPatch> cylinder_patch = make_ssi_patch(
      SurfaceGeometry{cylinder}, {0, 2 * 3.141592653589793}, {0, 1});
  const Result<SsiPatch> axial_plane = make_ssi_patch(
      SurfaceGeometry{Plane3{{0, 1, 0}, 0.0}}, {-2, 2}, {-2, 2});
  REQUIRE(cylinder_patch.has_value());
  REQUIRE(axial_plane.has_value());
  const Result<SsiResult> cylinder_result =
      ssi(cylinder_patch.value(), axial_plane.value(), SsiOptions{1e-9});
  REQUIRE(cylinder_result.has_value());
  REQUIRE(cylinder_result.value().curves.size() == 2);
  const Point3d& lower = cylinder_result.value().curves[0].segments.front().a;
  const Point3d& upper = cylinder_result.value().curves[1].segments.front().a;
  CHECK(lower[0] < upper[0]);
  CHECK(std::fabs(lower[0] + 1.0) < 1e-9);
  CHECK(std::fabs(upper[0] - 1.0) < 1e-9);

  // The same geometry with the cylinder frame rotated by 180 degrees gives the
  // two lines opposite parameter order (u = +2.09 for y = -0.866 and u = -2.09
  // for y = +0.866). The output order follows the geometry, not the parameters.
  const Result<SsiPatch> flipped_patch = make_ssi_patch(
      SurfaceGeometry{Cylinder3{Point3d{0, 0, 0}, Vec3d{0, 0, 1}, Vec3d{-1, 0, 0}, 1.0}},
      {0, 2 * 3.141592653589793}, {0, 1});
  const Result<SsiPatch> offset_plane = make_ssi_patch(
      SurfaceGeometry{Plane3{{1, 0, 0}, 0.5}}, {-2, 2}, {-2, 2});
  REQUIRE(flipped_patch.has_value());
  REQUIRE(offset_plane.has_value());
  const Result<SsiResult> flipped_result =
      ssi(flipped_patch.value(), offset_plane.value(), SsiOptions{1e-9});
  REQUIRE(flipped_result.has_value());
  REQUIRE(flipped_result.value().curves.size() == 2);
  const Point3d& flipped_lower = flipped_result.value().curves[0].segments.front().a;
  const Point3d& flipped_upper = flipped_result.value().curves[1].segments.front().a;
  CHECK(flipped_lower[1] < flipped_upper[1]);
  CHECK(std::fabs(flipped_lower[1] + 0.8660254037844386) < 1e-9);
  CHECK(std::fabs(flipped_upper[1] - 0.8660254037844386) < 1e-9);
}

TEST_CASE("ssi: a cut below the tessellation resolution is diagnosed") {
  const Sphere3 sphere{Point3d{0, 0, 0}, Vec3d{0, 0, 1}, Vec3d{1, 0, 0}, 1.0};
  const Result<SsiPatch> a = make_ssi_patch(
      SurfaceGeometry{sphere}, {0, 2 * 3.14159265358979323846}, {-0.5 * 3.14159265358979323846, 0.5 * 3.14159265358979323846});
  const Result<SsiPatch> b = make_ssi_patch(
      SurfaceGeometry{Plane3{{0, 0, 1}, 1.0 - 1e-10}}, {-2, 2}, {-2, 2});
  REQUIRE(a.has_value());
  REQUIRE(b.has_value());
  const Result<SsiResult> result = ssi(a.value(), b.value(), SsiOptions{1e-9});
  REQUIRE(result.has_value());
  CHECK(result.value().curves.empty());
  REQUIRE(result.value().issues.size() == 1);
  CHECK(result.value().issues[0].code == SsiDiagnostic::no_convergence);
  CHECK(result.value().candidate_pairs > 0);
}

TEST_CASE("ssi: coincident NURBS patches are diagnosed, not silently empty") {
  const NurbsSurface patch = planar_nurbs();
  const Result<SsiPatch> a = make_ssi_patch(SurfaceGeometry{patch}, {0, 1}, {0, 1});
  const Result<SsiPatch> b = make_ssi_patch(SurfaceGeometry{patch}, {0, 1}, {0, 1});
  REQUIRE(a.has_value());
  REQUIRE(b.has_value());
  const Result<SsiResult> result = ssi(a.value(), b.value(), SsiOptions{1e-9});
  REQUIRE(result.has_value());
  CHECK(result.value().curves.empty());
  REQUIRE(result.value().issues.size() == 1);
  CHECK(result.value().issues[0].code == SsiDiagnostic::no_convergence);
}
