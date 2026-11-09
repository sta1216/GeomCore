#include <array>
#include <cmath>
#include <limits>
#include <variant>
#include <vector>

#include <doctest/doctest.h>
#include <geom/brep/face_intersect.hpp>
#include <geom/brep/primitives.hpp>
#include <geom/brep/solid_editor.hpp>

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

// Builds a one-face open B-rep from four coplanar corners. The corner order
// must be cyclic and its winding must agree with `normal`.
Brep make_planar_face(const std::array<Point3d, 4>& corners, Vec3d normal,
                      double offset) {
  SolidEditor editor;
  const std::uint32_t surface = editor.add_surface(Plane3{normal, offset});
  std::array<VertexId, 4> vertices{};
  for (std::size_t i = 0; i < corners.size(); ++i) {
    vertices[i] = editor.add_vertex(corners[i], 0.0).value();
  }
  std::array<EdgeId, 4> edges{};
  for (std::size_t i = 0; i < corners.size(); ++i) {
    edges[i] = editor
                   .add_edge(vertices[i], vertices[(i + 1) % corners.size()],
                             kNoGeometry, 0.0, 1.0, 0.0)
                   .value();
  }
  const FaceId face = editor.add_face(surface, false, 0.0).value();
  std::array<CoedgeSpec, 4> specs{};
  for (std::size_t i = 0; i < specs.size(); ++i) specs[i].edge = edges[i];
  editor.add_loop(face, specs, true).value();
  const std::array<FaceId, 1> faces{face};
  editor.add_shell(faces, false).value();
  return editor.commit().value();
}

FaceId only_face(const Brep& brep) {
  const std::vector<FaceId> ids = brep.face_ids();
  return ids.empty() ? FaceId::null() : ids.front();
}

void check_point(const Point3d& actual, const Point3d& expected) {
  CHECK(std::fabs(actual[0] - expected[0]) <= 1e-12);
  CHECK(std::fabs(actual[1] - expected[1]) <= 1e-12);
  CHECK(std::fabs(actual[2] - expected[2]) <= 1e-12);
}

// point_in_face is strict, so the exact ends of a chord that runs into a rim
// are checked on the planes and just inside the rim instead.
void check_chord_in_both(const Brep& a, FaceId fa, const Brep& b, FaceId fb,
                         const Segment3& segment, double tol) {
  const auto plane_a = face_plane(a, fa);
  const auto plane_b = face_plane(b, fb);
  REQUIRE(plane_a.has_value());
  REQUIRE(plane_b.has_value());
  CHECK(point_on_plane(segment.a, plane_a.value(), tol));
  CHECK(point_on_plane(segment.a, plane_b.value(), tol));
  CHECK(point_on_plane(segment.b, plane_a.value(), tol));
  CHECK(point_on_plane(segment.b, plane_b.value(), tol));

  const Vec3d delta = segment.b - segment.a;
  for (const double f : {1e-6, 0.5, 1.0 - 1e-6}) {
    const Point3d point = segment.a + delta * f;
    const auto inside_a = point_in_face(a, fa, point, tol);
    const auto inside_b = point_in_face(b, fb, point, tol);
    REQUIRE(inside_a.has_value());
    REQUIRE(inside_b.has_value());
    CHECK(inside_a.value());
    CHECK(inside_b.value());
  }
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

TEST_CASE("face_intersect: tilted plane keeps the full chord in the unit line parameter") {
  const Brep a = make_box(Point3d{-1, -1, -1}, Point3d{1, 1, 0}, 1e-12).value();
  const Brep b = make_planar_face(
      std::array<Point3d, 4>{Point3d{-1, -1, -1}, Point3d{1, -1, 1},
                             Point3d{1, 1, 1}, Point3d{-1, 1, -1}},
      Vec3d{-1.0, 0.0, 1.0}, 0.0);
  const FaceId top = plane_face(a, Vec3d{0, 0, 1}, 0.0);
  const FaceId side = only_face(b);
  REQUIRE_FALSE(top.is_null());
  REQUIRE_FALSE(side.is_null());

  const auto cut = intersect_planar_faces(a, top, b, side, 1e-9);
  REQUIRE(cut.has_value());
  REQUIRE(cut.value().segments.size() == 1);
  const Segment3& segment = cut.value().segments[0];
  check_point(segment.a, Point3d{0.0, -1.0, 0.0});
  check_point(segment.b, Point3d{0.0, 1.0, 0.0});
  CHECK(segment_length(segment) == doctest::Approx(2.0).epsilon(1e-12));
  check_chord_in_both(a, top, b, side, segment, 1e-9);
}

TEST_CASE("face_intersect: translated tilted plane keeps the chord near the faces") {
  const Brep a = make_box(Point3d{9, 9, -1}, Point3d{11, 11, 0}, 1e-12).value();
  const Brep b = make_planar_face(
      std::array<Point3d, 4>{Point3d{9, 9, -1}, Point3d{11, 9, 1},
                             Point3d{11, 11, 1}, Point3d{9, 11, -1}},
      Vec3d{-1.0, 0.0, 1.0}, -10.0);
  const FaceId top = plane_face(a, Vec3d{0, 0, 1}, 0.0);
  const FaceId side = only_face(b);
  REQUIRE_FALSE(top.is_null());
  REQUIRE_FALSE(side.is_null());

  const auto cut = intersect_planar_faces(a, top, b, side, 1e-9);
  REQUIRE(cut.has_value());
  REQUIRE(cut.value().segments.size() == 1);
  const Segment3& segment = cut.value().segments[0];
  check_point(segment.a, Point3d{10.0, 9.0, 0.0});
  check_point(segment.b, Point3d{10.0, 11.0, 0.0});
  CHECK(segment_length(segment) == doctest::Approx(2.0).epsilon(1e-12));
  check_chord_in_both(a, top, b, side, segment, 1e-9);
}

TEST_CASE("face_intersect: shallow angle chord follows the analytic length") {
  const double theta = 0.1;
  const double sine = std::sin(theta);
  const double cosine = std::cos(theta);
  const Vec3d normal_b{-cosine * sine, sine * sine, cosine};
  const Vec3d line_direction{sine, cosine, 0.0};
  const Vec3d across = cross(normal_b, line_direction);
  const double half_length = 15.0;
  const double half_width = 5.0;
  const auto on_line = [&](double t) {
    return Point3d{t * sine, t * cosine, 0.0};
  };
  const Brep a = make_box(Point3d{-1, -20, -1}, Point3d{1, 20, 0}, 1e-12).value();
  const Brep b = make_planar_face(
      std::array<Point3d, 4>{
          on_line(-half_length) + across * half_width,
          on_line(-half_length) - across * half_width,
          on_line(half_length) - across * half_width,
          on_line(half_length) + across * half_width},
      normal_b, 0.0);
  const FaceId top = plane_face(a, Vec3d{0, 0, 1}, 0.0);
  const FaceId side = only_face(b);
  REQUIRE_FALSE(top.is_null());
  REQUIRE_FALSE(side.is_null());

  const auto cut = intersect_planar_faces(a, top, b, side, 1e-9);
  REQUIRE(cut.has_value());
  REQUIRE(cut.value().segments.size() == 1);
  const Segment3& segment = cut.value().segments[0];
  const double cotangent = cosine / sine;
  check_point(segment.a, Point3d{-1.0, -cotangent, 0.0});
  check_point(segment.b, Point3d{1.0, cotangent, 0.0});
  CHECK(segment_length(segment) == doctest::Approx(2.0 / sine).epsilon(1e-12));
  check_chord_in_both(a, top, b, side, segment, 1e-9);
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
  REQUIRE(cut.value().segments.size() == 4);

  const std::array<Point3d, 4> corners{Point3d{0, 0, 2}, Point3d{1, 0, 2},
                                       Point3d{1, 1, 2}, Point3d{0, 1, 2}};
  const auto is_corner = [&](const Point3d& point) {
    for (const Point3d& corner : corners) {
      if (std::fabs(point[0] - corner[0]) <= 1e-12 &&
          std::fabs(point[1] - corner[1]) <= 1e-12 &&
          std::fabs(point[2] - corner[2]) <= 1e-12) {
        return true;
      }
    }
    return false;
  };
  const auto plane_a = face_plane(a, top_a);
  const auto plane_b = face_plane(b, bottom_b);
  REQUIRE(plane_a.has_value());
  REQUIRE(plane_b.has_value());
  for (const Point3d& corner : corners) {
    CHECK(point_on_plane(corner, plane_a.value(), 1e-9));
    CHECK(point_on_plane(corner, plane_b.value(), 1e-9));
    const auto inside_a = point_in_face(a, top_a, corner, 1e-9);
    const auto rim_a = point_on_face_boundary(a, top_a, corner, 1e-9);
    const auto inside_b = point_in_face(b, bottom_b, corner, 1e-9);
    const auto rim_b = point_on_face_boundary(b, bottom_b, corner, 1e-9);
    REQUIRE(inside_a.has_value());
    REQUIRE(rim_a.has_value());
    REQUIRE(inside_b.has_value());
    REQUIRE(rim_b.has_value());
    const bool on_a = inside_a.value() || rim_a.value();
    const bool on_b = inside_b.value() || rim_b.value();
    CHECK(on_a);
    CHECK(on_b);
  }
  for (const Segment3& segment : cut.value().segments) {
    CHECK(segment_length(segment) == doctest::Approx(1.0).epsilon(1e-12));
    CHECK(is_corner(segment.a));
    CHECK(is_corner(segment.b));
  }
}

TEST_CASE("face_intersect: a crossing through a face corner is one segment") {
  const Brep a = make_box(Point3d{0, 0, 0}, Point3d{2, 2, 2}, 1e-12).value();
  const Brep b = make_planar_face(
      std::array<Point3d, 4>{Point3d{0, 0, 1}, Point3d{2, 2, 1},
                             Point3d{2, 2, 3}, Point3d{0, 0, 3}},
      Vec3d{1.0, -1.0, 0.0}, 0.0);
  const FaceId top = plane_face(a, Vec3d{0, 0, 1}, 2.0);
  const FaceId side = only_face(b);
  REQUIRE_FALSE(top.is_null());
  REQUIRE_FALSE(side.is_null());

  const auto cut = intersect_planar_faces(a, top, b, side, 1e-9);
  REQUIRE(cut.has_value());
  REQUIRE(cut.value().segments.size() == 1);
  const Segment3& segment = cut.value().segments[0];
  check_point(segment.a, Point3d{0.0, 0.0, 2.0});
  check_point(segment.b, Point3d{2.0, 2.0, 2.0});
  CHECK(segment_length(segment) == doctest::Approx(2.0 * std::sqrt(2.0)).epsilon(1e-12));
  check_chord_in_both(a, top, b, side, segment, 1e-9);
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
