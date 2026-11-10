#include <doctest/doctest.h>

#include <array>
#include <cmath>
#include <limits>
#include <string_view>
#include <type_traits>

#include <geom/brep/ssi.hpp>

using namespace geom::v1;

namespace {

constexpr double kPi = 3.14159265358979323846;
constexpr double kTwoPi = 2.0 * kPi;
constexpr double kHalfPi = 0.5 * kPi;
constexpr double kTolerance = 1e-9;

const Sphere3 kUnitSphere{Point3d{0, 0, 0}, Vec3d{0, 0, 1}, Vec3d{1, 0, 0}, 1.0};

SurfaceGeometry sphere_surface(const Sphere3& sphere) { return SurfaceGeometry{sphere}; }

Result<SsiPatch> full_sphere_patch(const Sphere3& sphere) {
  return make_ssi_patch(sphere_surface(sphere), {0.0, kTwoPi}, {-kHalfPi, kHalfPi});
}

Result<SsiPatch> square_plane_patch(const Vec3d& normal, double offset, double extent) {
  return make_ssi_patch(SurfaceGeometry{Plane3{normal, offset}}, {-extent, extent},
                        {-extent, extent});
}

Result<SsiPatch> full_cylinder_patch(const Cylinder3& cylinder) {
  return make_ssi_patch(SurfaceGeometry{cylinder}, {0.0, kTwoPi}, {0.0, 1.0});
}

Result<SsiPatch> full_cone_patch(const Cone3& cone) {
  return make_ssi_patch(SurfaceGeometry{cone}, {0.0, kTwoPi}, {0.0, 1.0});
}

Result<SsiPatch> full_torus_patch(const Torus3& torus) {
  return make_ssi_patch(SurfaceGeometry{torus}, {0.0, kTwoPi}, {0.0, kTwoPi});
}

}  // namespace

TEST_CASE("ssi: diagnostics header exposes the documented defaults") {
  using Entry = Result<SsiResult> (*)(const SsiPatch&, const SsiPatch&, const SsiOptions&);
  static_assert(std::is_same_v<decltype(&ssi), Entry>);

  const SsiOptions options{};
  CHECK(options.tolerance == 1e-9);
  CHECK(options.max_depth == 12);
  CHECK(options.max_segments == 100000);

  const SsiIssue issue{};
  CHECK(issue.code == SsiDiagnostic::none);
  CHECK(issue.params == std::array<double, 4>{});
  CHECK(std::string_view{issue.message}.empty());

  const SsiResult result{};
  CHECK(result.curves.empty());
  CHECK(result.issues.empty());
  CHECK(result.candidate_pairs == 0);
  CHECK(result.accepted_segments == 0);
  CHECK(result.rejected_segments == 0);
  CHECK(result.max_residual == 0.0);
  CHECK(SsiDiagnostic::tangent_contact != SsiDiagnostic::coincident);
}

TEST_CASE("ssi: a plane tangent to a sphere at its pole is tangent, not coincident") {
  // South pole is exactly (0,0,-1); the plane z=-1 touches there with parallel
  // normals. The 5x5 tangent grid includes both (0,0,-1) samples.
  const Result<SsiPatch> sphere = full_sphere_patch(kUnitSphere);
  const Result<SsiPatch> plane = square_plane_patch({0, 0, 1}, -1.0, 2.0);
  REQUIRE(sphere.has_value());
  REQUIRE(plane.has_value());
  CHECK(patches_tangent(sphere.value(), plane.value(), kTolerance));
  CHECK_FALSE(surfaces_coincident(sphere.value(), plane.value(), kTolerance));
}

TEST_CASE("ssi: identical planes are coincident, not tangent") {
  const Result<SsiPatch> a = square_plane_patch({0, 0, 1}, 2.0, 2.0);
  const Result<SsiPatch> b = square_plane_patch({0, 0, 3}, 6.0, 2.0);
  REQUIRE(a.has_value());
  REQUIRE(b.has_value());
  CHECK(surfaces_coincident(a.value(), b.value(), kTolerance));
  CHECK_FALSE(patches_tangent(a.value(), b.value(), kTolerance));
}

TEST_CASE("ssi: coincident planes with opposite normal orientation are coincident") {
  // Both planes are z=2, but b stores the opposite normal, so its unit offset
  // has the opposite sign too.
  const Result<SsiPatch> a = square_plane_patch({0, 0, 1}, 2.0, 2.0);
  const Result<SsiPatch> b = square_plane_patch({0, 0, -3}, -6.0, 2.0);
  REQUIRE(a.has_value());
  REQUIRE(b.has_value());
  CHECK(surfaces_coincident(a.value(), b.value(), kTolerance));
  CHECK_FALSE(patches_tangent(a.value(), b.value(), kTolerance));
}

TEST_CASE("ssi: coaxial equal-radius cylinders are coincident") {
  const Cylinder3 a_surface{Point3d{0, 0, 0}, Vec3d{0, 0, 1}, Vec3d{1, 0, 0}, 1.0};
  const Cylinder3 b_surface{Point3d{0, 0, 1}, Vec3d{0, 0, 1}, Vec3d{1, 0, 0}, 1.0};
  const Result<SsiPatch> a = full_cylinder_patch(a_surface);
  const Result<SsiPatch> b = full_cylinder_patch(b_surface);
  REQUIRE(a.has_value());
  REQUIRE(b.has_value());
  // The origins differ along the axis: a cylinder is one carrier, not a
  // particular origin. The carrier parameters must still agree.
  CHECK(surfaces_coincident(a.value(), b.value(), kTolerance));
  CHECK_FALSE(patches_tangent(a.value(), b.value(), kTolerance));
}

TEST_CASE("ssi: parallel offset planes are neither tangent nor coincident") {
  const Result<SsiPatch> a = square_plane_patch({0, 0, 1}, 0.0, 1.0);
  const Result<SsiPatch> b = square_plane_patch({0, 0, 1}, 1.0, 1.0);
  REQUIRE(a.has_value());
  REQUIRE(b.has_value());
  CHECK_FALSE(surfaces_coincident(a.value(), b.value(), kTolerance));
  CHECK_FALSE(patches_tangent(a.value(), b.value(), kTolerance));
}

TEST_CASE("ssi: crossing planes are neither tangent nor coincident") {
  const Result<SsiPatch> a = square_plane_patch({0, 0, 1}, 0.0, 2.0);
  const Result<SsiPatch> b = square_plane_patch({1, 0, 0}, 0.0, 2.0);
  REQUIRE(a.has_value());
  REQUIRE(b.has_value());
  // The grids share the origin, so the contact test passes; only the normal
  // test can reject this transversal intersection.
  CHECK_FALSE(surfaces_coincident(a.value(), b.value(), kTolerance));
  CHECK_FALSE(patches_tangent(a.value(), b.value(), kTolerance));
}

TEST_CASE("ssi: a sphere-plane circle is a transversal intersection, not a tangency") {
  const Result<SsiPatch> sphere = full_sphere_patch(kUnitSphere);
  const Result<SsiPatch> plane = square_plane_patch({0, 0, 1}, 0.25, 2.0);
  REQUIRE(sphere.has_value());
  REQUIRE(plane.has_value());
  CHECK_FALSE(surfaces_coincident(sphere.value(), plane.value(), kTolerance));
  CHECK_FALSE(patches_tangent(sphere.value(), plane.value(), kTolerance));
}

TEST_CASE("ssi: externally tangent spheres are tangent, not coincident") {
  const Sphere3 a_surface{Point3d{0, 0, 0}, Vec3d{0, 0, 1}, Vec3d{1, 0, 0}, 1.0};
  const Sphere3 b_surface{Point3d{2, 0, 0}, Vec3d{0, 0, 1}, Vec3d{1, 0, 0}, 1.0};
  const Result<SsiPatch> a = full_sphere_patch(a_surface);
  const Result<SsiPatch> b = full_sphere_patch(b_surface);
  REQUIRE(a.has_value());
  REQUIRE(b.has_value());
  CHECK(patches_tangent(a.value(), b.value(), kTolerance));
  CHECK_FALSE(surfaces_coincident(a.value(), b.value(), kTolerance));
}

TEST_CASE("ssi: internally tangent spheres are tangent, not coincident") {
  const Sphere3 a_surface{Point3d{0, 0, 0}, Vec3d{0, 0, 1}, Vec3d{1, 0, 0}, 2.0};
  const Sphere3 b_surface{Point3d{1, 0, 0}, Vec3d{0, 0, 1}, Vec3d{1, 0, 0}, 1.0};
  const Result<SsiPatch> a = full_sphere_patch(a_surface);
  const Result<SsiPatch> b = full_sphere_patch(b_surface);
  REQUIRE(a.has_value());
  REQUIRE(b.has_value());
  CHECK(patches_tangent(a.value(), b.value(), kTolerance));
  CHECK_FALSE(surfaces_coincident(a.value(), b.value(), kTolerance));
}

TEST_CASE("ssi: separated spheres are neither tangent nor coincident") {
  const Sphere3 a_surface{Point3d{0, 0, 0}, Vec3d{0, 0, 1}, Vec3d{1, 0, 0}, 1.0};
  const Sphere3 b_surface{Point3d{3, 0, 0}, Vec3d{0, 0, 1}, Vec3d{1, 0, 0}, 1.0};
  const Result<SsiPatch> a = full_sphere_patch(a_surface);
  const Result<SsiPatch> b = full_sphere_patch(b_surface);
  REQUIRE(a.has_value());
  REQUIRE(b.has_value());
  CHECK_FALSE(patches_tangent(a.value(), b.value(), kTolerance));
  CHECK_FALSE(surfaces_coincident(a.value(), b.value(), kTolerance));
}

TEST_CASE("ssi: interpenetrating spheres are neither tangent nor coincident") {
  const Sphere3 a_surface{Point3d{0, 0, 0}, Vec3d{0, 0, 1}, Vec3d{1, 0, 0}, 1.0};
  const Sphere3 b_surface{Point3d{0.5, 0, 0}, Vec3d{0, 0, 1}, Vec3d{1, 0, 0}, 1.0};
  const Result<SsiPatch> a = full_sphere_patch(a_surface);
  const Result<SsiPatch> b = full_sphere_patch(b_surface);
  REQUIRE(a.has_value());
  REQUIRE(b.has_value());
  CHECK_FALSE(patches_tangent(a.value(), b.value(), kTolerance));
  CHECK_FALSE(surfaces_coincident(a.value(), b.value(), kTolerance));
}

TEST_CASE("ssi: a plane patch that misses the tangent point is not tangent") {
  // The plane z=-1 is tangent at the sphere's south pole (0,0,-1), but this
  // patch only covers y in [-3,-1]; the nearest sample is (0,-1,-1), which is
  // sqrt(2)-1 away from the sphere. A single-point sample at the patch origin
  // would not catch this; the 5x5 grid does.
  const Result<SsiPatch> sphere = full_sphere_patch(kUnitSphere);
  const Result<SsiPatch> plane =
      make_ssi_patch(SurfaceGeometry{Plane3{{0, 0, 1}, -1.0}}, {1.0, 3.0}, {-1.0, 1.0});
  REQUIRE(sphere.has_value());
  REQUIRE(plane.has_value());
  CHECK_FALSE(patches_tangent(sphere.value(), plane.value(), kTolerance));
}

TEST_CASE("ssi: equal spheres are coincident, not tangent") {
  const Result<SsiPatch> a = full_sphere_patch(kUnitSphere);
  const Result<SsiPatch> b = full_sphere_patch(kUnitSphere);
  REQUIRE(a.has_value());
  REQUIRE(b.has_value());
  CHECK(surfaces_coincident(a.value(), b.value(), kTolerance));
  CHECK_FALSE(patches_tangent(a.value(), b.value(), kTolerance));
}

TEST_CASE("ssi: coaxially reparameterized equal cones are coincident") {
  // Both describe the same double cone with its apex at (0,0,1): the cone
  // parameters origin/radius trade off along the axis. The apex, axis and
  // half angle are the carrier invariants.
  const Cone3 a_surface{Point3d{0, 0, 0}, Vec3d{0, 0, 1}, Vec3d{1, 0, 0}, 1.0, 0.25 * kPi};
  const Cone3 b_surface{Point3d{0, 0, -0.5}, Vec3d{0, 0, 1}, Vec3d{1, 0, 0}, 1.5, 0.25 * kPi};
  const Result<SsiPatch> a = full_cone_patch(a_surface);
  const Result<SsiPatch> b = full_cone_patch(b_surface);
  REQUIRE(a.has_value());
  REQUIRE(b.has_value());
  CHECK(surfaces_coincident(a.value(), b.value(), kTolerance));
  CHECK_FALSE(patches_tangent(a.value(), b.value(), kTolerance));
}

TEST_CASE("ssi: equal tori are coincident, not tangent") {
  const Torus3 torus{Point3d{0, 0, 0}, Vec3d{0, 0, 1}, Vec3d{1, 0, 0}, 2.0, 0.5};
  const Result<SsiPatch> a = full_torus_patch(torus);
  const Result<SsiPatch> b = full_torus_patch(torus);
  REQUIRE(a.has_value());
  REQUIRE(b.has_value());
  CHECK(surfaces_coincident(a.value(), b.value(), kTolerance));
  CHECK_FALSE(patches_tangent(a.value(), b.value(), kTolerance));
}

TEST_CASE("ssi: equal carriers with a different defining parameter are not coincident") {
  const Result<SsiPatch> unit_sphere_patch = full_sphere_patch(kUnitSphere);
  const Sphere3 big_sphere{Point3d{0, 0, 0}, Vec3d{0, 0, 1}, Vec3d{1, 0, 0}, 1.001};
  const Result<SsiPatch> big_sphere_patch = full_sphere_patch(big_sphere);

  const Cylinder3 thin_cylinder{Point3d{0, 0, 0}, Vec3d{0, 0, 1}, Vec3d{1, 0, 0}, 1.0};
  const Cylinder3 fat_cylinder{Point3d{0, 0, 0}, Vec3d{0, 0, 1}, Vec3d{1, 0, 0}, 1.001};
  const Result<SsiPatch> thin_cylinder_patch = full_cylinder_patch(thin_cylinder);
  const Result<SsiPatch> fat_cylinder_patch = full_cylinder_patch(fat_cylinder);

  const Cone3 steep_cone{Point3d{0, 0, 0}, Vec3d{0, 0, 1}, Vec3d{1, 0, 0}, 1.0, 0.25 * kPi};
  const Cone3 shallow_cone{Point3d{0, 0, 0}, Vec3d{0, 0, 1}, Vec3d{1, 0, 0}, 1.0, kPi / 6.0};
  const Result<SsiPatch> steep_cone_patch = full_cone_patch(steep_cone);
  const Result<SsiPatch> shallow_cone_patch = full_cone_patch(shallow_cone);

  const Torus3 thick_torus{Point3d{0, 0, 0}, Vec3d{0, 0, 1}, Vec3d{1, 0, 0}, 2.0, 0.5};
  const Torus3 thin_torus{Point3d{0, 0, 0}, Vec3d{0, 0, 1}, Vec3d{1, 0, 0}, 2.0, 0.45};
  const Result<SsiPatch> thick_torus_patch = full_torus_patch(thick_torus);
  const Result<SsiPatch> thin_torus_patch = full_torus_patch(thin_torus);

  REQUIRE(unit_sphere_patch.has_value());
  REQUIRE(big_sphere_patch.has_value());
  REQUIRE(thin_cylinder_patch.has_value());
  REQUIRE(fat_cylinder_patch.has_value());
  REQUIRE(steep_cone_patch.has_value());
  REQUIRE(shallow_cone_patch.has_value());
  REQUIRE(thick_torus_patch.has_value());
  REQUIRE(thin_torus_patch.has_value());
  CHECK_FALSE(surfaces_coincident(unit_sphere_patch.value(), big_sphere_patch.value(), kTolerance));
  CHECK_FALSE(
      surfaces_coincident(thin_cylinder_patch.value(), fat_cylinder_patch.value(), kTolerance));
  CHECK_FALSE(
      surfaces_coincident(steep_cone_patch.value(), shallow_cone_patch.value(), kTolerance));
  CHECK_FALSE(
      surfaces_coincident(thick_torus_patch.value(), thin_torus_patch.value(), kTolerance));
}

TEST_CASE("ssi: non-finite and non-positive tolerances are rejected") {
  const Result<SsiPatch> sphere = full_sphere_patch(kUnitSphere);
  const Result<SsiPatch> plane = square_plane_patch({0, 0, 1}, -1.0, 2.0);
  REQUIRE(sphere.has_value());
  REQUIRE(plane.has_value());
  const double nan = std::numeric_limits<double>::quiet_NaN();
  const double inf = std::numeric_limits<double>::infinity();
  for (const double tolerance : {nan, inf, -inf, 0.0, -1e-9}) {
    CHECK_FALSE(surfaces_coincident(sphere.value(), plane.value(), tolerance));
    CHECK_FALSE(patches_tangent(sphere.value(), plane.value(), tolerance));
  }
}

TEST_CASE("ssi: patches with non-finite or inverted domains are rejected") {
  const Result<SsiPatch> sphere = full_sphere_patch(kUnitSphere);
  const Result<SsiPatch> plane = square_plane_patch({0, 0, 1}, -1.0, 2.0);
  REQUIRE(sphere.has_value());
  REQUIRE(plane.has_value());

  const double nan = std::numeric_limits<double>::quiet_NaN();
  const double inf = std::numeric_limits<double>::infinity();

  SsiPatch nan_bound = sphere.value();
  nan_bound.u[1] = nan;
  SsiPatch infinite_bound = plane.value();
  infinite_bound.v[0] = -inf;
  SsiPatch inverted = sphere.value();
  inverted.u = {1.0, -1.0};
  SsiPatch invalid_surface{
      SurfaceGeometry{Sphere3{Point3d{0, 0, 0}, Vec3d{0, 0, 1}, Vec3d{1, 0, 0}, nan}},
      {0.0, 1.0},
      {-1.0, 1.0}};

  // The baseline pair is a genuine tangency, so every rejection below is due to
  // the malformed patch rather than a universally false predicate.
  CHECK(patches_tangent(sphere.value(), plane.value(), kTolerance));
  for (const SsiPatch& bad : {nan_bound, infinite_bound, inverted, invalid_surface}) {
    CHECK_FALSE(surfaces_coincident(bad, plane.value(), kTolerance));
    CHECK_FALSE(surfaces_coincident(sphere.value(), bad, kTolerance));
    CHECK_FALSE(patches_tangent(bad, plane.value(), kTolerance));
    CHECK_FALSE(patches_tangent(sphere.value(), bad, kTolerance));
  }
}

TEST_CASE("ssi: a patch whose evaluation overflows is rejected") {
  const Cylinder3 huge{Point3d{1e308, 0, 0}, Vec3d{0, 0, 1}, Vec3d{1, 0, 0}, 1e308};
  REQUIRE(huge.is_valid());
  const Result<SsiPatch> patch = make_ssi_patch(SurfaceGeometry{huge}, {-1.0, 1.0}, {-1.0, 1.0});
  REQUIRE(patch.has_value());
  CHECK_FALSE(surfaces_coincident(patch.value(), patch.value(), kTolerance));
  CHECK_FALSE(patches_tangent(patch.value(), patch.value(), kTolerance));
}

TEST_CASE("ssi: coincidence of NURBS surfaces is documented as unsupported") {
  const NurbsSurface surface(1, 1, {0.0, 0.0, 1.0, 1.0}, {0.0, 0.0, 1.0, 1.0},
                             {Point3d{0, 0, 0}, Point3d{0, 1, 0}, Point3d{1, 0, 0},
                              Point3d{1, 1, 0}},
                             {1.0, 1.0, 1.0, 1.0});
  REQUIRE(surface.is_valid());
  const Result<SsiPatch> a = make_ssi_patch(SurfaceGeometry{surface}, {0.0, 1.0}, {0.0, 1.0});
  REQUIRE(a.has_value());
  // Identical NURBS patches are deliberately not classified: sampling could not
  // distinguish coincidence from tangency, so both predicates report false.
  CHECK_FALSE(surfaces_coincident(a.value(), a.value(), kTolerance));
  CHECK_FALSE(patches_tangent(a.value(), a.value(), kTolerance));
}
