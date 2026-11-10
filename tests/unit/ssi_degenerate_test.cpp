#include <doctest/doctest.h>

#include <array>
#include <cmath>
#include <cstddef>
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

Vec3d rotate_about(const Vec3d& value, const Vec3d& axis, double angle) {
  const double cosine = std::cos(angle);
  const double sine = std::sin(angle);
  return value * cosine + cross(axis, value) * sine +
         axis * (dot(axis, value) * (1.0 - cosine));
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

TEST_CASE("ssi: near-parallel crossing planes are not tangent") {
  // Both planes pass through the origin, so the gap at the contact is exactly
  // zero. The 1e-7 rad crossing is far above the strict 1e-10 normal guard and
  // the global sign test sees both sides across the crossing line.
  for (const double angle : {1e-7, 1e-6}) {
    CAPTURE(angle);
    const Result<SsiPatch> a = square_plane_patch({0, 0, 1}, 0.0, 2.0);
    const Result<SsiPatch> b =
        square_plane_patch({std::sin(angle), 0.0, std::cos(angle)}, 0.0, 2.0);
    REQUIRE(a.has_value());
    REQUIRE(b.has_value());
    CHECK_FALSE(surfaces_coincident(a.value(), b.value(), kTolerance));
    CHECK_FALSE(patches_tangent(a.value(), b.value(), kTolerance));
  }
}

TEST_CASE("ssi: a plane cutting a sphere just below its pole is not tangent") {
  // z = 1 - 0.5*tolerance cuts a circle of radius sqrt(tolerance) ~ 3.2e-5. The
  // minimum gap is inside the band and the normals are parallel, so only the
  // one-sided test can reject it.
  const Result<SsiPatch> sphere = full_sphere_patch(kUnitSphere);
  const Result<SsiPatch> plane =
      square_plane_patch({0, 0, 1}, 1.0 - 0.5 * kTolerance, 2.0);
  REQUIRE(sphere.has_value());
  REQUIRE(plane.has_value());
  CHECK_FALSE(patches_tangent(sphere.value(), plane.value(), kTolerance));
}

TEST_CASE("ssi: parallel planes half a tolerance apart are coincident, not tangent") {
  const Result<SsiPatch> a = square_plane_patch({0, 0, 1}, 0.0, 2.0);
  const Result<SsiPatch> b = square_plane_patch({0, 0, 1}, 0.5 * kTolerance, 2.0);
  REQUIRE(a.has_value());
  REQUIRE(b.has_value());
  // The separation lies inside the coincidence band, so the carriers count as
  // the same plane; tangency is reserved for contact between distinct carriers.
  CHECK(surfaces_coincident(a.value(), b.value(), kTolerance));
  CHECK_FALSE(patches_tangent(a.value(), b.value(), kTolerance));
}

TEST_CASE("ssi: a sphere tangent to a plane at its equator with a rotated frame is tangent") {
  const double angle = 0.3;
  // Rotating the frame about the contact normal x moves the poles and the seam,
  // not the sphere or the contact point (1, 0, 0).
  const Sphere3 sphere{Point3d{0, 0, 0}, Vec3d{0.0, std::cos(angle), std::sin(angle)},
                       Vec3d{0.0, -std::sin(angle), std::cos(angle)}, 1.0};
  REQUIRE(sphere.is_valid());
  const Result<SsiPatch> sphere_patch = full_sphere_patch(sphere);
  const Result<SsiPatch> plane = square_plane_patch({1, 0, 0}, 1.0, 2.0);
  REQUIRE(sphere_patch.has_value());
  REQUIRE(plane.has_value());
  CHECK(patches_tangent(sphere_patch.value(), plane.value(), kTolerance));
  CHECK_FALSE(surfaces_coincident(sphere_patch.value(), plane.value(), kTolerance));
}

TEST_CASE("ssi: a cylinder tangent to a plane along a line with a rotated frame is tangent") {
  const double angle = 0.3;
  const Cylinder3 cylinder{Point3d{0, 0, 0}, Vec3d{0, 0, 1},
                           Vec3d{std::cos(angle), std::sin(angle), 0.0}, 1.0};
  REQUIRE(cylinder.is_valid());
  const Result<SsiPatch> cylinder_patch = full_cylinder_patch(cylinder);
  const Result<SsiPatch> plane = square_plane_patch({1, 0, 0}, 1.0, 2.0);
  REQUIRE(cylinder_patch.has_value());
  REQUIRE(plane.has_value());
  CHECK(patches_tangent(cylinder_patch.value(), plane.value(), kTolerance));
  CHECK_FALSE(surfaces_coincident(cylinder_patch.value(), plane.value(), kTolerance));
}

TEST_CASE("ssi: a torus tangent to a plane away from the plane sample grid is tangent") {
  // The contact circle has radius R = 2, off the plane patch's 9x9 sample grid,
  // so refinement has to walk the pair onto the contact.
  const Torus3 torus{Point3d{0, 0, 0}, Vec3d{0, 0, 1}, Vec3d{1, 0, 0}, 2.0, 0.5};
  const Result<SsiPatch> torus_patch = full_torus_patch(torus);
  const Result<SsiPatch> plane = square_plane_patch({0, 0, 1}, 0.5, 3.0);
  REQUIRE(torus_patch.has_value());
  REQUIRE(plane.has_value());
  CHECK(patches_tangent(torus_patch.value(), plane.value(), kTolerance));
  CHECK_FALSE(surfaces_coincident(torus_patch.value(), plane.value(), kTolerance));
}

TEST_CASE("ssi: tangent spheres with rotated frames are tangent") {
  const double angle = 0.3;
  const Vec3d axis{0.0, std::cos(angle), std::sin(angle)};
  const Vec3d x_axis{0.0, -std::sin(angle), std::cos(angle)};

  const Result<SsiPatch> external_a =
      full_sphere_patch(Sphere3{Point3d{0, 0, 0}, axis, x_axis, 1.0});
  const Result<SsiPatch> external_b =
      full_sphere_patch(Sphere3{Point3d{2, 0, 0}, axis, x_axis, 1.0});
  REQUIRE(external_a.has_value());
  REQUIRE(external_b.has_value());
  CHECK(patches_tangent(external_a.value(), external_b.value(), kTolerance));

  const Result<SsiPatch> internal_a =
      full_sphere_patch(Sphere3{Point3d{0, 0, 0}, axis, x_axis, 2.0});
  const Result<SsiPatch> internal_b =
      full_sphere_patch(Sphere3{Point3d{1, 0, 0}, axis, x_axis, 1.0});
  REQUIRE(internal_a.has_value());
  REQUIRE(internal_b.has_value());
  CHECK(patches_tangent(internal_a.value(), internal_b.value(), kTolerance));
}

TEST_CASE("ssi: coincidence ignores the seam parameter of symmetric carriers") {
  const double angle = 0.3;
  const Vec3d tilted_x{std::cos(angle), std::sin(angle), 0.0};

  // A sphere's axis only names its poles: rotating the axis and the seam
  // together leaves the same sphere.
  const Result<SsiPatch> sphere_a =
      full_sphere_patch(Sphere3{Point3d{0, 0, 0}, Vec3d{0, 0, 1}, Vec3d{1, 0, 0}, 1.0});
  const Result<SsiPatch> sphere_b = full_sphere_patch(
      Sphere3{Point3d{0, 0, 0}, Vec3d{0.0, std::cos(angle), std::sin(angle)},
              Vec3d{0.0, -std::sin(angle), std::cos(angle)}, 1.0});
  REQUIRE(sphere_a.has_value());
  REQUIRE(sphere_b.has_value());
  CHECK(surfaces_coincident(sphere_a.value(), sphere_b.value(), kTolerance));

  const Result<SsiPatch> cylinder_a = full_cylinder_patch(
      Cylinder3{Point3d{0, 0, 0}, Vec3d{0, 0, 1}, Vec3d{1, 0, 0}, 1.0});
  const Result<SsiPatch> cylinder_b = full_cylinder_patch(
      Cylinder3{Point3d{0, 0, 0}, Vec3d{0, 0, 1}, tilted_x, 1.0});
  REQUIRE(cylinder_a.has_value());
  REQUIRE(cylinder_b.has_value());
  CHECK(surfaces_coincident(cylinder_a.value(), cylinder_b.value(), kTolerance));
  CHECK_FALSE(patches_tangent(cylinder_a.value(), cylinder_b.value(), kTolerance));

  const Result<SsiPatch> cone_a = full_cone_patch(
      Cone3{Point3d{0, 0, 0}, Vec3d{0, 0, 1}, Vec3d{1, 0, 0}, 1.0, 0.25 * kPi});
  const Result<SsiPatch> cone_b = full_cone_patch(
      Cone3{Point3d{0, 0, 0}, Vec3d{0, 0, 1}, tilted_x, 1.0, 0.25 * kPi});
  REQUIRE(cone_a.has_value());
  REQUIRE(cone_b.has_value());
  CHECK(surfaces_coincident(cone_a.value(), cone_b.value(), kTolerance));
  CHECK_FALSE(patches_tangent(cone_a.value(), cone_b.value(), kTolerance));

  const Result<SsiPatch> torus_a = full_torus_patch(
      Torus3{Point3d{0, 0, 0}, Vec3d{0, 0, 1}, Vec3d{1, 0, 0}, 2.0, 0.5});
  const Result<SsiPatch> torus_b = full_torus_patch(
      Torus3{Point3d{0, 0, 0}, Vec3d{0, 0, 1}, tilted_x, 2.0, 0.5});
  REQUIRE(torus_a.has_value());
  REQUIRE(torus_b.has_value());
  CHECK(surfaces_coincident(torus_a.value(), torus_b.value(), kTolerance));
  CHECK_FALSE(patches_tangent(torus_a.value(), torus_b.value(), kTolerance));
}

TEST_CASE("ssi: coincidence samples between the 5x5 nodes and rejects a mid-cell escape") {
  // B is A rotated by a tiny angle about a tangent direction at u = pi/4. The
  // deviation 2*epsilon*|sin v * cos(u - u0)| keeps every 5x5 node inside the
  // band while several 9x9 mid-cell samples leave it, so a node-only test would
  // wrongly report coincidence.
  const Torus3 torus{Point3d{0, 0, 0}, Vec3d{0, 0, 1}, Vec3d{1, 0, 0}, 2.0, 0.5};
  const double epsilon = 6e-10;
  const double spin = 0.25 * kPi;
  const Vec3d tilt{-std::sin(spin), std::cos(spin), 0.0};
  const Torus3 tilted{Point3d{0, 0, 0}, rotate_about({0, 0, 1}, tilt, epsilon),
                      rotate_about({1, 0, 0}, tilt, epsilon), 2.0, 0.5};
  REQUIRE(tilted.is_valid());

  const Result<SsiPatch> a =
      make_ssi_patch(SurfaceGeometry{torus}, {0.0, kTwoPi}, {0.1, 2.0});
  const Result<SsiPatch> b = full_torus_patch(tilted);
  REQUIRE(a.has_value());
  REQUIRE(b.has_value());

  const auto distance = [](const Torus3& carrier, const Point3d& point) {
    const Vec3d delta = point - carrier.center;
    const double axial = dot(delta, carrier.axis);
    const double radial = norm2(delta - carrier.axis * axial);
    return std::fabs(std::hypot(radial - carrier.major_radius, axial) - carrier.minor_radius);
  };

  double worst_node = 0.0;
  for (int iu = 0; iu < 5; ++iu) {
    for (int iv = 0; iv < 5; ++iv) {
      const double u = (kTwoPi / 4.0) * iu;
      const double v = 0.1 + (1.9 / 4.0) * iv;
      const Result<SurfacePoint> evaluated = surface_eval(a.value().surface, u, v);
      REQUIRE(evaluated.has_value());
      worst_node = std::max(worst_node, distance(tilted, evaluated.value().p));
    }
  }
  CHECK(worst_node <= kTolerance);

  std::size_t outside = 0;
  for (int iu = 0; iu < 9; ++iu) {
    for (int iv = 0; iv < 9; ++iv) {
      const double u = (kTwoPi / 8.0) * iu;
      const double v = 0.1 + (1.9 / 8.0) * iv;
      const Result<SurfacePoint> evaluated = surface_eval(a.value().surface, u, v);
      REQUIRE(evaluated.has_value());
      if (distance(tilted, evaluated.value().p) > kTolerance) {
        ++outside;
      }
    }
  }
  CHECK(outside >= 1);
  CHECK_FALSE(surfaces_coincident(a.value(), b.value(), kTolerance));
}

TEST_CASE("ssi: coincidence rejects a domain sample outside the band") {
  // Two near-parallel planes whose 9x9 grid has exactly one sample outside the
  // band: the (u=-2, v=2) corner, where v - u = 4. A coincidence test that
  // tolerates a stray sample would report the pair as coincident.
  const double slope = 0.13e-9;
  const Vec3d normal{slope, slope, std::sqrt(1.0 - 2.0 * slope * slope)};
  const Result<SsiPatch> a = square_plane_patch({0, 0, 1}, 0.0, 2.0);
  const Result<SsiPatch> b = square_plane_patch(normal, -0.5 * kTolerance, 2.0);
  REQUIRE(a.has_value());
  REQUIRE(b.has_value());

  std::size_t outside = 0;
  for (int iu = 0; iu < 9; ++iu) {
    for (int iv = 0; iv < 9; ++iv) {
      const double u = -2.0 + (4.0 / 8.0) * iu;
      const double v = -2.0 + (4.0 / 8.0) * iv;
      const Point3d point{v, -u, 0.0};
      if (std::fabs(dot(normal, Vec3d{point[0], point[1], point[2]}) + 0.5 * kTolerance) >
          kTolerance) {
        ++outside;
      }
    }
  }
  CHECK(outside == 1);
  CHECK_FALSE(surfaces_coincident(a.value(), b.value(), kTolerance));
}

TEST_CASE("ssi: crossing planes at a nanosized angle are not tangent") {
  // The normal mismatch (1e-9 rad) is above the strict 1e-10 normal guard, and
  // the patches straddle the crossing line, so the global sign test also sees
  // both sides. The off-origin variant crosses near z = 1 with unequal extents.
  const double angle = 1e-9;
  const Result<SsiPatch> origin_a = square_plane_patch({0, 0, 1}, 0.0, 2.0);
  const Result<SsiPatch> origin_b =
      square_plane_patch({std::sin(angle), 0.0, std::cos(angle)}, 0.0, 3.0);
  REQUIRE(origin_a.has_value());
  REQUIRE(origin_b.has_value());
  CHECK_FALSE(patches_tangent(origin_a.value(), origin_b.value(), kTolerance));
  CHECK_FALSE(patches_tangent(origin_b.value(), origin_a.value(), kTolerance));

  const double off_angle = 2.0 * angle;
  const Result<SsiPatch> offset_a = square_plane_patch({0, 0, 1}, 1.0, 2.0);
  const Result<SsiPatch> offset_b =
      square_plane_patch({std::sin(off_angle), 0.0, std::cos(off_angle)}, 1.0, 3.0);
  REQUIRE(offset_a.has_value());
  REQUIRE(offset_b.has_value());
  CHECK_FALSE(patches_tangent(offset_a.value(), offset_b.value(), kTolerance));
}

TEST_CASE("ssi: a sphere cut within the tolerance band is not tangent") {
  // The plane sits a fraction of the tolerance below the pole, so the gap is
  // inside the band and the normals are parallel. The cap above the plane is
  // only a fraction of the tolerance deep, which the tolerance band alone would
  // treat as neutral; the geometric band sees it while the far side of the
  // sphere is negative, so the global sign test rejects every cut.
  const Result<SsiPatch> sphere = full_sphere_patch(kUnitSphere);
  REQUIRE(sphere.has_value());
  for (const double fraction : {0.1, 0.2, 0.25, 0.3, 0.5, 0.9, 1.5}) {
    CAPTURE(fraction);
    const Result<SsiPatch> plane =
        square_plane_patch({0, 0, 1}, 1.0 - fraction * kTolerance, 2.0);
    REQUIRE(plane.has_value());
    CHECK_FALSE(patches_tangent(sphere.value(), plane.value(), kTolerance));
  }
  // Exactly at the pole the pair touches from one side and is tangent.
  const Result<SsiPatch> touching = square_plane_patch({0, 0, 1}, 1.0, 2.0);
  REQUIRE(touching.has_value());
  CHECK(patches_tangent(sphere.value(), touching.value(), kTolerance));
}

TEST_CASE("ssi: cone coincidence scales the half-angle comparison by the carrier size") {
  // The half angle is an angle: a 1e-10 rad difference is negligible on a
  // millimetre-scale cone and invisible on a kilometre-scale one, so the
  // comparison is scaled by the carrier radius (translation invariant).
  const double tiny = 1e-3;
  const Result<SsiPatch> tiny_a = make_ssi_patch(
      SurfaceGeometry{
          Cone3{Point3d{0, 0, 0}, Vec3d{0, 0, 1}, Vec3d{1, 0, 0}, tiny, 0.25 * kPi}},
      {0.0, kTwoPi}, {0.0, tiny});
  const Result<SsiPatch> tiny_b = make_ssi_patch(
      SurfaceGeometry{
          Cone3{Point3d{0, 0, 0}, Vec3d{0, 0, 1}, Vec3d{1, 0, 0}, tiny,
                0.25 * kPi + 1e-10}},
      {0.0, kTwoPi}, {0.0, tiny});
  REQUIRE(tiny_a.has_value());
  REQUIRE(tiny_b.has_value());
  CHECK(surfaces_coincident(tiny_a.value(), tiny_b.value(), 1e-12));

  const Result<SsiPatch> different = make_ssi_patch(
      SurfaceGeometry{
          Cone3{Point3d{0, 0, 0}, Vec3d{0, 0, 1}, Vec3d{1, 0, 0}, tiny,
                0.25 * kPi + 1e-6}},
      {0.0, kTwoPi}, {0.0, tiny});
  REQUIRE(different.has_value());
  CHECK_FALSE(surfaces_coincident(tiny_a.value(), different.value(), 1e-12));

  const double large = 1e3;
  const Result<SsiPatch> large_a = make_ssi_patch(
      SurfaceGeometry{
          Cone3{Point3d{0, 0, 0}, Vec3d{0, 0, 1}, Vec3d{1, 0, 0}, large, 0.25 * kPi}},
      {0.0, kTwoPi}, {0.0, large});
  const Result<SsiPatch> large_b = make_ssi_patch(
      SurfaceGeometry{
          Cone3{Point3d{0, 0, 0}, Vec3d{0, 0, 1}, Vec3d{1, 0, 0}, large,
                0.25 * kPi + 1e-10}},
      {0.0, kTwoPi}, {0.0, large});
  REQUIRE(large_a.has_value());
  REQUIRE(large_b.has_value());
  CHECK(surfaces_coincident(large_a.value(), large_b.value(), 1e-6));

  // The half-angle scale is the local characteristic radius, not the absolute
  // apex position: translating the same pair must not change the verdict. The
  // translated checks use tolerance 1e-9 because coordinates near 1e6 quantize
  // to ~1.2e-10 in double (ulp(1e6) = 2.3e-10), so no sampled coincidence
  // below that scale is representable; the angular comparison itself is
  // unaffected. With the previous apex-position scaling, 1e-10 rad * 1e6 =
  // 1e-4 would fail these checks.
  for (const double translation : {1e3, 1e6}) {
    CAPTURE(translation);
    const Result<SsiPatch> moved_a = make_ssi_patch(
        SurfaceGeometry{Cone3{Point3d{translation, 0, 0}, Vec3d{0, 0, 1},
                              Vec3d{1, 0, 0}, tiny, 0.25 * kPi}},
        {0.0, kTwoPi}, {0.0, tiny});
    const Result<SsiPatch> moved_b = make_ssi_patch(
        SurfaceGeometry{Cone3{Point3d{translation, 0, 0}, Vec3d{0, 0, 1},
                              Vec3d{1, 0, 0}, tiny, 0.25 * kPi + 1e-10}},
        {0.0, kTwoPi}, {0.0, tiny});
    REQUIRE(moved_a.has_value());
    REQUIRE(moved_b.has_value());
    CHECK(surfaces_coincident(moved_a.value(), moved_b.value(), 1e-9));
  }
  const Result<SsiPatch> moved_different = make_ssi_patch(
      SurfaceGeometry{Cone3{Point3d{1e6, 0, 0}, Vec3d{0, 0, 1}, Vec3d{1, 0, 0}, tiny,
                            0.25 * kPi + 1e-6}},
      {0.0, kTwoPi}, {0.0, tiny});
  REQUIRE(moved_different.has_value());
  const Result<SsiPatch> moved_reference = make_ssi_patch(
      SurfaceGeometry{Cone3{Point3d{1e6, 0, 0}, Vec3d{0, 0, 1}, Vec3d{1, 0, 0}, tiny,
                            0.25 * kPi}},
      {0.0, kTwoPi}, {0.0, tiny});
  REQUIRE(moved_reference.has_value());
  CHECK_FALSE(surfaces_coincident(moved_reference.value(), moved_different.value(), 1e-9));
}

TEST_CASE("ssi: the global sign test covers a finite cylinder patch") {
  const Result<SsiPatch> cylinder = make_ssi_patch(
      SurfaceGeometry{
          Cylinder3{Point3d{0, 0, 0}, Vec3d{0, 0, 1}, Vec3d{1, 0, 0}, 1.0}},
      {0.0, kTwoPi}, {-0.5, 0.5});
  REQUIRE(cylinder.has_value());
  const Result<SsiPatch> tangent_plane = square_plane_patch({1, 0, 0}, 1.0, 2.0);
  REQUIRE(tangent_plane.has_value());
  CHECK(patches_tangent(cylinder.value(), tangent_plane.value(), kTolerance));

  // Half a tolerance inside the cylinder the carriers cut along a strip: the
  // gap is zero, but the samples of the cylinder span both sides of the plane.
  const Result<SsiPatch> cutting_plane =
      square_plane_patch({1, 0, 0}, 1.0 - 0.5 * kTolerance, 2.0);
  REQUIRE(cutting_plane.has_value());
  CHECK_FALSE(patches_tangent(cylinder.value(), cutting_plane.value(), kTolerance));
}

TEST_CASE("ssi: a shallow crossing inside the coincidence band is still not tangent") {
  // The patches are small enough that a 1e-7 rad crossing moves their samples
  // by only ~1e-11, and the equal offsets make the carriers the same plane
  // within the tolerance: the pair is coincident, so tangency is rejected at
  // that gate. It is still a crossing, not a contact.
  const double angle = 1e-7;
  const double extent = 1e-4;
  const Result<SsiPatch> a = square_plane_patch({0, 0, 1}, 0.0, extent);
  const Result<SsiPatch> b =
      square_plane_patch({std::sin(angle), 0.0, std::cos(angle)}, 0.0, extent);
  REQUIRE(a.has_value());
  REQUIRE(b.has_value());
  CHECK(surfaces_coincident(a.value(), b.value(), kTolerance));
  CHECK_FALSE(patches_tangent(a.value(), b.value(), kTolerance));
}

TEST_CASE("ssi: a boundary-crossing plane pair is not tangent") {
  // The two patches share the x = 0 boundary and cross along it: plane B is
  // z = tan(1e-7) * x, plane A is z = 0. They are not coincident (the separation
  // grows to 2e-7 over the domain), and the contact at the shared boundary has
  // normals 1e-7 rad apart. Without the strict normal guard the alignment pass
  // accepts that mismatch, the gap at the boundary is zero, and the one-sided
  // test sees only non-negative distances on both patches, so this fixture is
  // what pins the guard.
  const double angle = 1e-7;
  const Result<SsiPatch> a = make_ssi_patch(
      SurfaceGeometry{Plane3{{0, 0, 1}, 0.0}}, {-2.0, 2.0}, {0.0, 2.0});
  const Result<SsiPatch> b = make_ssi_patch(
      SurfaceGeometry{Plane3{{-std::sin(angle), 0.0, std::cos(angle)}, 0.0}},
      {0.0, 2.0}, {-2.0, 2.0});
  REQUIRE(a.has_value());
  REQUIRE(b.has_value());
  CHECK_FALSE(surfaces_coincident(a.value(), b.value(), kTolerance));
  CHECK_FALSE(patches_tangent(a.value(), b.value(), kTolerance));
  CHECK_FALSE(patches_tangent(b.value(), a.value(), kTolerance));
}

TEST_CASE("ssi: a cut below the documented detection floor is tangent contact") {
  // The geometric band max(1e-12 * (1 + max|p|), 0.05 * tolerance) is the
  // detection floor: a cut shallower than that leaves every sample of the cut
  // carrier inside the band, so it is reported as tangent_contact by contract.
  // The nominal threshold is 0.05 * tolerance; exactly 0.05 is not tested
  // because the final comparison is ulp-fragile.
  const Result<SsiPatch> sphere = full_sphere_patch(kUnitSphere);
  REQUIRE(sphere.has_value());
  const Result<SsiPatch> inside_floor =
      square_plane_patch({0, 0, 1}, 1.0 - 0.04 * kTolerance, 2.0);
  REQUIRE(inside_floor.has_value());
  CHECK(patches_tangent(sphere.value(), inside_floor.value(), kTolerance));

  for (const double fraction : {0.06, 0.1, 0.25, 0.5}) {
    CAPTURE(fraction);
    const Result<SsiPatch> cutting =
        square_plane_patch({0, 0, 1}, 1.0 - fraction * kTolerance, 2.0);
    REQUIRE(cutting.has_value());
    CHECK_FALSE(patches_tangent(sphere.value(), cutting.value(), kTolerance));
  }

  const Result<SsiPatch> pole = square_plane_patch({0, 0, 1}, 1.0, 2.0);
  REQUIRE(pole.has_value());
  CHECK(patches_tangent(sphere.value(), pole.value(), kTolerance));
}

TEST_CASE("ssi: QA-9 degenerate corpus is either stitched or diagnosed") {
  int cases = 0;
  int curves = 0;
  int diagnosed = 0;
  int tangent = 0;
  int coincident = 0;
  const auto check_result = [&](const SsiResult& result) {
    ++cases;
    if (result.curves.empty()) {
      REQUIRE_FALSE(result.issues.empty());
      for (const SsiIssue& issue : result.issues) {
        CHECK(issue.code != SsiDiagnostic::none);
      }
      if (result.issues[0].code == SsiDiagnostic::tangent_contact) ++tangent;
      if (result.issues[0].code == SsiDiagnostic::coincident) ++coincident;
      ++diagnosed;
      return;
    }
    ++curves;
    for (const SsiCurve& curve : result.curves) {
      for (const SsiSegment& s : curve.segments) CHECK(s.residual <= 1e-8);
    }
  };

  // Family 1: sphere tangent to a plane at its south pole (50 radii).
  for (int i = 0; i < 50; ++i) {
    const double r = 0.5 + 0.01 * static_cast<double>(i);
    const Result<SsiPatch> a = make_ssi_patch(
        SurfaceGeometry{Sphere3{Point3d{0, 0, 0}, Vec3d{0, 0, 1}, Vec3d{1, 0, 0}, r}},
        {0, 2 * kPi}, {-0.5 * kPi, 0.5 * kPi});
    const Result<SsiPatch> b = make_ssi_patch(
        SurfaceGeometry{Plane3{{0, 0, 1}, -r}}, {-3, 3}, {-3, 3});
    REQUIRE(a.has_value());
    REQUIRE(b.has_value());
    const Result<SsiResult> result = ssi(a.value(), b.value(), SsiOptions{1e-8});
    REQUIRE(result.has_value());
    check_result(result.value());
  }

  // Family 2: coincident planes written with different normals/scales (50).
  for (int i = 1; i <= 50; ++i) {
    const double d = 0.1 * static_cast<double>(i);
    const Result<SsiPatch> a = make_ssi_patch(
        SurfaceGeometry{Plane3{{0, 0, 1}, d}}, {-2, 2}, {-2, 2});
    const Result<SsiPatch> b = make_ssi_patch(
        SurfaceGeometry{Plane3{{0, 0, -3}, -3 * d}}, {-2, 2}, {-2, 2});
    REQUIRE(a.has_value());
    REQUIRE(b.has_value());
    const Result<SsiResult> result = ssi(a.value(), b.value(), SsiOptions{1e-8});
    REQUIRE(result.has_value());
    check_result(result.value());
  }

  // Family 3: coaxial equal-radius cylinders, offset along the axis (50).
  for (int i = 0; i < 50; ++i) {
    const Cylinder3 a_surface{Point3d{0, 0, 0}, Vec3d{0, 0, 1}, Vec3d{1, 0, 0}, 1.0};
    const Cylinder3 b_surface{Point3d{0, 0, 0.01 * static_cast<double>(i)},
                              Vec3d{0, 0, 1}, Vec3d{1, 0, 0}, 1.0};
    const Result<SsiPatch> a = make_ssi_patch(
        SurfaceGeometry{a_surface}, {0, 2 * kPi}, {0, 1});
    const Result<SsiPatch> b = make_ssi_patch(
        SurfaceGeometry{b_surface}, {0, 2 * kPi}, {0, 1});
    REQUIRE(a.has_value());
    REQUIRE(b.has_value());
    const Result<SsiResult> result = ssi(a.value(), b.value(), SsiOptions{1e-8});
    REQUIRE(result.has_value());
    check_result(result.value());
  }

  // Family 4: torus tangent to the plane z = minor_radius (50 major radii).
  for (int i = 0; i < 50; ++i) {
    const double major = 2.0 + 0.01 * static_cast<double>(i);
    const Torus3 torus{Point3d{0, 0, 0}, Vec3d{0, 0, 1}, Vec3d{1, 0, 0}, major, 0.5};
    const Result<SsiPatch> a = make_ssi_patch(
        SurfaceGeometry{torus}, {0, 2 * kPi}, {0, 2 * kPi});
    const Result<SsiPatch> b = make_ssi_patch(
        SurfaceGeometry{Plane3{{0, 0, 1}, 0.5}}, {-4, 4}, {-4, 4});
    REQUIRE(a.has_value());
    REQUIRE(b.has_value());
    const Result<SsiResult> result = ssi(a.value(), b.value(), SsiOptions{1e-8});
    REQUIRE(result.has_value());
    check_result(result.value());
  }

  CHECK(curves == 0);
  CHECK(cases == 200);
  CHECK(diagnosed == 200);
  CHECK(tangent == 100);
  CHECK(coincident == 100);
}
