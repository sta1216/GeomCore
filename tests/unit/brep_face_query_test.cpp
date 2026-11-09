#include <cmath>
#include <cstddef>
#include <limits>
#include <variant>
#include <vector>

#include <doctest/doctest.h>
#include <geom/brep/face_query.hpp>
#include <geom/brep/primitives.hpp>

using namespace geom;

namespace {
// Finds the face whose plane has the given outward unit normal and offset.
// Copy this helper into every test file that needs to name a face.
FaceId plane_face(const Brep& brep, Vec3d want_normal, double want_offset) {
  for (const FaceId id : brep.face_ids()) {
    const auto& variant = brep.geometry().surfaces[brep.face(id).surface];
    const auto* plane = std::get_if<Plane3>(&variant);
    if (plane == nullptr) continue;
    const double len = std::sqrt(dot(plane->normal, plane->normal));
    const Vec3d unit{plane->normal[0] / len, plane->normal[1] / len,
                     plane->normal[2] / len};
    if (dot(unit, want_normal) > 0.999 &&
        std::fabs(plane->offset / len - want_offset) <= 1e-9) {
      return id;
    }
  }
  return FaceId::null();
}
}

TEST_CASE("face_query: a box top face contains its centre, not its rim") {
  const Brep box = make_box(Point3d{0.0, 0.0, 0.0}, Point3d{2.0, 2.0, 2.0}, 1e-12).value();
  const FaceId top = plane_face(box, Vec3d{0.0, 0.0, 1.0}, 2.0);
  REQUIRE_FALSE(top.is_null());
  REQUIRE(point_in_face(box, top, Point3d{1.0, 1.0, 2.0}, 1e-9).has_value());
  CHECK(point_in_face(box, top, Point3d{1.0, 1.0, 2.0}, 1e-9).value());
  CHECK_FALSE(point_in_face(box, top, Point3d{0.0, 0.0, 2.0}, 1e-9).value());
  CHECK_FALSE(point_in_face(box, top, Point3d{3.0, 1.0, 2.0}, 1e-9).value());
  CHECK(point_on_face_boundary(box, top, Point3d{0.0, 1.0, 2.0}, 1e-9).value());
}

TEST_CASE("face_query: a shell contains an interior point and not an exterior one") {
  const Brep box = make_box(Point3d{0.0, 0.0, 0.0}, Point3d{2.0, 2.0, 2.0}, 1e-12).value();
  const ShellId shell = box.shell_ids().front();
  CHECK(point_in_shell(box, shell, Point3d{1.0, 1.0, 1.0}, 1e-9).value());
  CHECK_FALSE(point_in_shell(box, shell, Point3d{5.0, 1.0, 1.0}, 1e-9).value());
}

TEST_CASE("face_query: box top face projects to one four-point ring of area four") {
  const Brep box = make_box(Point3d{0.0, 0.0, 0.0}, Point3d{2.0, 2.0, 2.0}, 1e-12).value();
  const FaceId top = plane_face(box, Vec3d{0.0, 0.0, 1.0}, 2.0);
  REQUIRE_FALSE(top.is_null());
  const Result<Plane3> plane = face_plane(box, top);
  REQUIRE(plane.has_value());
  const Result<PlaneFrame> frame = make_plane_frame(plane.value());
  REQUIRE(frame.has_value());

  const Result<std::vector<std::vector<Point2d>>> loops =
      face_loops_2d(box, top, frame.value());
  REQUIRE(loops.has_value());
  REQUIRE(loops.value().size() == std::size_t{1});
  const std::vector<Point2d>& ring = loops.value().front();
  REQUIRE(ring.size() == std::size_t{4});

  double twice_area = 0.0;
  for (std::size_t i = 0; i < ring.size(); ++i) {
    const Point2d& a = ring[i];
    const Point2d& b = ring[(i + 1) % ring.size()];
    twice_area += a[0] * b[1] - b[0] * a[1];
    // to_3d must land back on the face plane, z = 2 for the box top face.
    const Point3d back = plane_frame_to_3d(frame.value(), a);
    CHECK(back[2] == doctest::Approx(2.0).epsilon(1e-12));
  }
  CHECK(std::fabs(0.5 * twice_area) == doctest::Approx(4.0).epsilon(1e-12));
}

TEST_CASE("face_query: box top plane is the unit plus-z plane at offset two") {
  const Brep box = make_box(Point3d{0.0, 0.0, 0.0}, Point3d{2.0, 2.0, 2.0}, 1e-12).value();
  const FaceId top = plane_face(box, Vec3d{0.0, 0.0, 1.0}, 2.0);
  REQUIRE_FALSE(top.is_null());
  const Result<Plane3> plane = face_plane(box, top);
  REQUIRE(plane.has_value());
  const double len = std::sqrt(dot(plane.value().normal, plane.value().normal));
  CHECK(len == doctest::Approx(1.0).epsilon(1e-12));
  CHECK(plane.value().normal[0] == doctest::Approx(0.0).epsilon(1e-12));
  CHECK(plane.value().normal[1] == doctest::Approx(0.0).epsilon(1e-12));
  CHECK(plane.value().normal[2] == doctest::Approx(1.0).epsilon(1e-12));
  CHECK(plane.value().offset / len == doctest::Approx(2.0).epsilon(1e-12));

  // A cylinder's lateral surface is a Cylinder3 variant, not a Plane3.
  const Brep cylinder =
      make_cylinder(Point3d{0.0, 0.0, 0.0}, Vec3d{0.0, 0.0, 1.0}, 1.0, 2.0, 1e-12).value();
  bool saw_curved = false;
  for (const FaceId id : cylinder.face_ids()) {
    const auto& surface = cylinder.geometry().surfaces[cylinder.face(id).surface];
    if (std::get_if<Cylinder3>(&surface) == nullptr) continue;
    const Result<Plane3> curved = face_plane(cylinder, id);
    REQUIRE_FALSE(curved.has_value());
    CHECK(curved.error().code == Errc::unsupported);
    saw_curved = true;
  }
  CHECK(saw_curved);
}

TEST_CASE("face_query: boundary uses point-segment distance within the tolerance") {
  const Brep box = make_box(Point3d{0.0, 0.0, 0.0}, Point3d{2.0, 2.0, 2.0}, 1e-12).value();
  const FaceId top = plane_face(box, Vec3d{0.0, 0.0, 1.0}, 2.0);
  REQUIRE_FALSE(top.is_null());

  const Result<bool> corner =
      point_on_face_boundary(box, top, Point3d{0.0, 0.0, 2.0}, 1e-9);
  REQUIRE(corner.has_value());
  CHECK(corner.value());

  const Result<bool> just_outside =
      point_on_face_boundary(box, top, Point3d{2.0 + 1e-6, 1.0, 2.0}, 1e-9);
  REQUIRE(just_outside.has_value());
  CHECK_FALSE(just_outside.value());

  // One micrometre off the plane is far outside a 1e-9 tolerance.
  const Result<bool> off_plane =
      point_on_face_boundary(box, top, Point3d{1.0, 1.0, 2.0 + 1e-6}, 1e-9);
  REQUIRE(off_plane.has_value());
  CHECK_FALSE(off_plane.value());
}

TEST_CASE("face_query: a shell excludes a point just outside a face") {
  const Brep box = make_box(Point3d{0.0, 0.0, 0.0}, Point3d{2.0, 2.0, 2.0}, 1e-12).value();
  const ShellId shell = box.shell_ids().front();
  const Result<bool> outside =
      point_in_shell(box, shell, Point3d{1.0, 1.0, 2.0 + 1e-6}, 1e-9);
  REQUIRE(outside.has_value());
  CHECK_FALSE(outside.value());
}

TEST_CASE("face_query: the shell boundary is reported as outside") {
  const Brep box = make_box(Point3d{0.0, 0.0, 0.0}, Point3d{2.0, 2.0, 2.0}, 1e-12).value();
  const ShellId shell = box.shell_ids().front();
  // A point exactly on the top face is on the shell boundary. The parity
  // convention here skips plane hits with t <= tol, so the first ray direction
  // (positive z) leaves the solid without a crossing and reports outside.
  const Result<bool> boundary = point_in_shell(box, shell, Point3d{1.0, 1.0, 2.0}, 1e-9);
  REQUIRE(boundary.has_value());
  CHECK_FALSE(boundary.value());
}

TEST_CASE("face_query: a shell with a curved face is unsupported") {
  const Brep cylinder =
      make_cylinder(Point3d{0.0, 0.0, 0.0}, Vec3d{0.0, 0.0, 1.0}, 1.0, 2.0, 1e-12).value();
  const ShellId shell = cylinder.shell_ids().front();
  const Result<bool> inside =
      point_in_shell(cylinder, shell, Point3d{0.0, 0.0, 1.0}, 1e-9);
  REQUIRE_FALSE(inside.has_value());
  CHECK(inside.error().code == Errc::unsupported);
}

TEST_CASE("face_query: invalid tolerances and non-finite points are rejected") {
  const Brep box = make_box(Point3d{0.0, 0.0, 0.0}, Point3d{2.0, 2.0, 2.0}, 1e-12).value();
  const FaceId top = plane_face(box, Vec3d{0.0, 0.0, 1.0}, 2.0);
  REQUIRE_FALSE(top.is_null());
  const double nan = std::numeric_limits<double>::quiet_NaN();

  const Result<bool> nan_tol = point_in_face(box, top, Point3d{1.0, 1.0, 2.0}, nan);
  REQUIRE_FALSE(nan_tol.has_value());
  CHECK(nan_tol.error().code == Errc::tolerance_violation);

  const Result<bool> negative_tol = point_in_face(box, top, Point3d{1.0, 1.0, 2.0}, -1e-9);
  REQUIRE_FALSE(negative_tol.has_value());
  CHECK(negative_tol.error().code == Errc::tolerance_violation);

  const Result<bool> nan_point = point_in_face(box, top, Point3d{nan, 1.0, 2.0}, 1e-9);
  REQUIRE_FALSE(nan_point.has_value());
  CHECK(nan_point.error().code == Errc::non_finite);
}

TEST_CASE("face_query: null handles are index_out_of_range") {
  const Brep box = make_box(Point3d{0.0, 0.0, 0.0}, Point3d{2.0, 2.0, 2.0}, 1e-12).value();

  const Result<Plane3> null_plane = face_plane(box, FaceId::null());
  REQUIRE_FALSE(null_plane.has_value());
  CHECK(null_plane.error().code == Errc::index_out_of_range);

  const Result<bool> null_face =
      point_in_face(box, FaceId::null(), Point3d{1.0, 1.0, 2.0}, 1e-9);
  REQUIRE_FALSE(null_face.has_value());
  CHECK(null_face.error().code == Errc::index_out_of_range);

  const Result<bool> null_shell =
      point_in_shell(box, ShellId::null(), Point3d{1.0, 1.0, 1.0}, 1e-9);
  REQUIRE_FALSE(null_shell.has_value());
  CHECK(null_shell.error().code == Errc::index_out_of_range);
}
