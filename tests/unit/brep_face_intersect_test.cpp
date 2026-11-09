#include <cmath>
#include <limits>
#include <variant>
#include <vector>

#include <doctest/doctest.h>
#include <geom/brep/face_intersect.hpp>
#include <geom/brep/primitives.hpp>

using namespace geom;

namespace {
// Finds the face whose plane has the given outward normal direction and whose
// distance from the origin is want_offset. The plan names that distance with a
// positive sign even for an inward-facing normal (for example the bottom face
// of a box), so compare magnitudes.
FaceId plane_face(const Brep& brep, Vec3d want_normal, double want_offset) {
  for (const FaceId id : brep.face_ids()) {
    const auto& variant = brep.geometry().surfaces[brep.face(id).surface];
    const auto* plane = std::get_if<Plane3>(&variant);
    if (plane == nullptr) continue;
    const double len = std::sqrt(dot(plane->normal, plane->normal));
    const Vec3d unit{plane->normal[0] / len, plane->normal[1] / len,
                     plane->normal[2] / len};
    if (dot(unit, want_normal) > 0.999 &&
        std::fabs(std::fabs(plane->offset / len) - want_offset) <= 1e-9) {
      return id;
    }
  }
  return FaceId::null();
}

// First face whose surface is not a plane; a cylinder's lateral face.
FaceId curved_face(const Brep& brep) {
  for (const FaceId id : brep.face_ids()) {
    const auto& variant = brep.geometry().surfaces[brep.face(id).surface];
    if (std::get_if<Plane3>(&variant) == nullptr) return id;
  }
  return FaceId::null();
}

double segment_length(const Segment3& segment) {
  const Vec3d delta = segment.b - segment.a;
  return std::hypot(delta[0], delta[1], delta[2]);
}

bool point_on_plane(const Point3d& point, const Plane3& plane, double tol) {
  const double len = std::hypot(plane.normal[0], plane.normal[1], plane.normal[2]);
  const Vec3d unit{plane.normal[0] / len, plane.normal[1] / len,
                   plane.normal[2] / len};
  return std::fabs(dot(unit, Vec3d{point[0], point[1], point[2]}) - plane.offset / len) <=
         tol;
}
}  // namespace

TEST_CASE("face_intersect: perpendicular boxes cut a top face along one segment") {
  const Brep a = make_box(Point3d{0, 0, 0}, Point3d{2, 2, 2}, 1e-12).value();
  const Brep b = make_box(Point3d{1, 1, 0}, Point3d{3, 3, 3}, 1e-12).value();
  const FaceId top = plane_face(a, Vec3d{0.0, 0.0, 1.0}, 2.0);
  const FaceId side = plane_face(b, Vec3d{-1.0, 0.0, 0.0}, 1.0);
  REQUIRE_FALSE(top.is_null());
  REQUIRE_FALSE(side.is_null());
  const auto cut = intersect_planar_faces(a, top, b, side, 1e-12);
  REQUIRE(cut.has_value());
  CHECK_FALSE(cut.value().coplanar);
  REQUIRE(cut.value().segments.size() == 1);
  CHECK(cut.value().segments[0].a[0] == doctest::Approx(1.0).epsilon(1e-12));
  CHECK(cut.value().segments[0].a[2] == doctest::Approx(2.0).epsilon(1e-12));
  CHECK(cut.value().segments[0].b[1] == doctest::Approx(2.0).epsilon(1e-12));

  // The visible cut is shorter than either face and lies on both planes.
  CHECK(segment_length(cut.value().segments[0]) == doctest::Approx(1.0).epsilon(1e-12));
  const auto plane_a = face_plane(a, top);
  const auto plane_b = face_plane(b, side);
  REQUIRE(plane_a.has_value());
  REQUIRE(plane_b.has_value());
  CHECK(point_on_plane(cut.value().segments[0].a, plane_a.value(), 1e-12));
  CHECK(point_on_plane(cut.value().segments[0].a, plane_b.value(), 1e-12));
  CHECK(point_on_plane(cut.value().segments[0].b, plane_a.value(), 1e-12));
  CHECK(point_on_plane(cut.value().segments[0].b, plane_b.value(), 1e-12));
}

TEST_CASE("face_intersect: coplanar faces return the overlap boundary") {
  const Brep a = make_box(Point3d{0, 0, 0}, Point3d{2, 2, 2}, 1e-12).value();
  const Brep b = make_box(Point3d{0, 0, 2}, Point3d{1, 1, 4}, 1e-12).value();
  const FaceId top_a = plane_face(a, Vec3d{0.0, 0.0, 1.0}, 2.0);
  const FaceId bottom_b = plane_face(b, Vec3d{0.0, 0.0, -1.0}, 2.0);
  REQUIRE_FALSE(top_a.is_null());
  REQUIRE_FALSE(bottom_b.is_null());
  const auto cut = intersect_planar_faces(a, top_a, b, bottom_b, 1e-12);
  REQUIRE(cut.has_value());
  CHECK(cut.value().coplanar);
  CHECK(cut.value().segments.size() == 4);
}

TEST_CASE("face_intersect: disjoint faces produce no segments") {
  const Brep a = make_box(Point3d{0, 0, 0}, Point3d{1, 1, 1}, 1e-12).value();
  const Brep b = make_box(Point3d{5, 5, 5}, Point3d{6, 6, 6}, 1e-12).value();
  const auto cut = intersect_planar_faces(
      a, plane_face(a, Vec3d{0, 0, 1}, 1.0), b, plane_face(b, Vec3d{0, 0, -1}, 5.0), 1e-12);
  REQUIRE(cut.has_value());
  CHECK(cut.value().segments.empty());
  CHECK_FALSE(cut.value().coplanar);
}

TEST_CASE("face_intersect: parallel faces on different planes produce no cut") {
  const Brep a = make_box(Point3d{0, 0, 0}, Point3d{2, 2, 2}, 1e-12).value();
  const Brep b = make_box(Point3d{0, 0, 0}, Point3d{2, 2, 1}, 1e-12).value();
  const FaceId top = plane_face(a, Vec3d{0, 0, 1}, 2.0);
  const FaceId bottom = plane_face(b, Vec3d{0, 0, -1}, 0.0);
  REQUIRE_FALSE(top.is_null());
  REQUIRE_FALSE(bottom.is_null());
  const auto cut = intersect_planar_faces(a, top, b, bottom, 1e-12);
  REQUIRE(cut.has_value());
  CHECK_FALSE(cut.value().coplanar);
  CHECK(cut.value().segments.empty());
}

TEST_CASE("face_intersect: coplanar disjoint faces produce no cut") {
  const Brep a = make_box(Point3d{0, 0, 0}, Point3d{2, 2, 2}, 1e-12).value();
  const Brep b = make_box(Point3d{5, 5, 2}, Point3d{7, 7, 4}, 1e-12).value();
  const FaceId top = plane_face(a, Vec3d{0, 0, 1}, 2.0);
  const FaceId bottom = plane_face(b, Vec3d{0, 0, -1}, 2.0);
  REQUIRE_FALSE(top.is_null());
  REQUIRE_FALSE(bottom.is_null());
  const auto cut = intersect_planar_faces(a, top, b, bottom, 1e-12);
  REQUIRE(cut.has_value());
  CHECK(cut.value().coplanar);
  CHECK(cut.value().segments.empty());
}

TEST_CASE("face_intersect: coplanar overlap can be smaller than either face") {
  const Brep a = make_box(Point3d{0, 0, 0}, Point3d{4, 4, 2}, 1e-12).value();
  const Brep b = make_box(Point3d{1, 1, 2}, Point3d{3, 3, 4}, 1e-12).value();
  const FaceId top = plane_face(a, Vec3d{0, 0, 1}, 2.0);
  const FaceId bottom = plane_face(b, Vec3d{0, 0, -1}, 2.0);
  REQUIRE_FALSE(top.is_null());
  REQUIRE_FALSE(bottom.is_null());
  const auto cut = intersect_planar_faces(a, top, b, bottom, 1e-12);
  REQUIRE(cut.has_value());
  CHECK(cut.value().coplanar);
  REQUIRE(cut.value().segments.size() == 4);
  double total = 0.0;
  for (const Segment3& segment : cut.value().segments) total += segment_length(segment);
  CHECK(total == doctest::Approx(8.0).epsilon(1e-12));
}

TEST_CASE("face_intersect: a tangential rim touch produces no cut") {
  const Brep a = make_box(Point3d{0, 0, 0}, Point3d{2, 2, 2}, 1e-12).value();
  const Brep b = make_box(Point3d{2, 0, 0}, Point3d{3, 2, 3}, 1e-12).value();
  const FaceId top = plane_face(a, Vec3d{0, 0, 1}, 2.0);
  const FaceId side = plane_face(b, Vec3d{-1, 0, 0}, 2.0);
  REQUIRE_FALSE(top.is_null());
  REQUIRE_FALSE(side.is_null());
  const auto cut = intersect_planar_faces(a, top, b, side, 1e-12);
  REQUIRE(cut.has_value());
  CHECK_FALSE(cut.value().coplanar);
  CHECK(cut.value().segments.empty());
}

TEST_CASE("face_intersect: curved faces and invalid tolerances are rejected") {
  const Brep box = make_box(Point3d{0, 0, 0}, Point3d{2, 2, 2}, 1e-12).value();
  const Brep cylinder =
      make_cylinder(Point3d{0, 0, 0}, Vec3d{0, 0, 1}, 1.0, 2.0, 1e-12).value();
  const FaceId top = plane_face(box, Vec3d{0, 0, 1}, 2.0);
  const FaceId bottom = plane_face(box, Vec3d{0, 0, -1}, 0.0);
  const FaceId curved = curved_face(cylinder);
  REQUIRE_FALSE(top.is_null());
  REQUIRE_FALSE(bottom.is_null());
  REQUIRE_FALSE(curved.is_null());

  const auto unsupported = intersect_planar_faces(cylinder, curved, box, top, 1e-12);
  REQUIRE_FALSE(unsupported.has_value());
  CHECK(unsupported.error().code == Errc::unsupported);

  const auto negative = intersect_planar_faces(box, top, box, bottom, -1e-9);
  REQUIRE_FALSE(negative.has_value());
  CHECK(negative.error().code == Errc::tolerance_violation);

  const auto nan = intersect_planar_faces(
      box, top, box, bottom, std::numeric_limits<double>::quiet_NaN());
  REQUIRE_FALSE(nan.has_value());
  CHECK(nan.error().code == Errc::tolerance_violation);

  const auto null_face = intersect_planar_faces(box, FaceId::null(), box, bottom, 1e-12);
  REQUIRE_FALSE(null_face.has_value());
  CHECK(null_face.error().code == Errc::index_out_of_range);
}
