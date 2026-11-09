#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <span>
#include <variant>
#include <vector>

#include <doctest/doctest.h>
#include <geom/brep/planar_cut.hpp>
#include <geom/brep/primitives.hpp>
#include <geom/brep/solid_editor.hpp>
#include <geom/cgeom/polygon2.hpp>
#include <geom/cgeom/triangulate2.hpp>
#include <geom/numeric/predicates2.hpp>

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

double shoelace(std::span<const Point2d> ring) {
  double sum = 0.0;
  for (std::size_t i = 0; i < ring.size(); ++i) {
    const Point2d& a = ring[i];
    const Point2d& b = ring[(i + 1) % ring.size()];
    sum += a[0] * b[1] - b[0] * a[1];
  }
  return 0.5 * sum;
}

double fragment_area(const FaceFragment2d& f) {
  double area = std::fabs(shoelace(f.outer));
  for (const auto& hole : f.holes) area -= std::fabs(shoelace(hole));
  return area;
}

using Segments = std::vector<Segment3>;

Polygon2d fragment_polygon(const FaceFragment2d& f) {
  Polygon2d poly;
  poly.ring_starts.push_back(0u);
  poly.vertices.insert(poly.vertices.end(), f.outer.begin(), f.outer.end());
  poly.ring_starts.push_back(static_cast<std::uint32_t>(poly.vertices.size()));
  for (const auto& hole : f.holes) {
    poly.vertices.insert(poly.vertices.end(), hole.begin(), hole.end());
    poly.ring_starts.push_back(static_cast<std::uint32_t>(poly.vertices.size()));
  }
  return poly;
}

bool in_triangle(const Point2d& a, const Point2d& b, const Point2d& c,
                 const Point2d& p) {
  return sign_orient2d(a[0], a[1], b[0], b[1], p[0], p[1]) >= 0 &&
         sign_orient2d(b[0], b[1], c[0], c[1], p[0], p[1]) >= 0 &&
         sign_orient2d(c[0], c[1], a[0], a[1], p[0], p[1]) >= 0;
}

// The sample must be inside the original face and inside the fragment itself.
// The fragment is checked by triangulating its rings; the sample is the
// centroid of one of those triangles, so it must lie in a closed triangle.
void check_samples_inside(const Brep& brep, FaceId id,
                          const std::vector<FaceFragment2d>& fragments, double tol) {
  const auto plane = face_plane(brep, id);
  CHECK(plane.has_value());
  if (!plane.has_value()) return;
  const auto frame = make_plane_frame(plane.value());
  CHECK(frame.has_value());
  if (!frame.has_value()) return;

  for (const FaceFragment2d& f : fragments) {
    const Point3d sample3 = plane_frame_to_3d(frame.value(), f.sample);
    const auto inside_face = point_in_face(brep, id, sample3, tol);
    CHECK(inside_face.has_value());
    if (inside_face.has_value()) CHECK(inside_face.value());

    const Polygon2d poly = fragment_polygon(f);
    CHECK(validate_polygon2(poly) == Errc::none);
    const Triangulation2 tri = triangulate_polygon2(poly);
    CHECK(tri.ok);
    if (!tri.ok) continue;
    bool found = false;
    for (const auto& t : tri.triangles) {
      if (in_triangle(poly.vertices[t[0]], poly.vertices[t[1]], poly.vertices[t[2]],
                      f.sample)) {
        found = true;
        break;
      }
    }
    CHECK(found);
  }
}
// Builds a single planar face from an outer ring and hole rings. The outer ring
// must wind counter-clockwise around normal and each hole ring clockwise, both
// in world space. Used by the keyhole test, which a box cannot express.
Brep make_planar_face(std::span<const Point3d> outer,
                      const std::vector<std::vector<Point3d>>& holes,
                      Vec3d normal, double offset) {
  SolidEditor editor;
  const std::uint32_t surface = editor.add_surface(Plane3{normal, offset});
  const FaceId face = editor.add_face(surface, false, 0.0).value();
  const auto add_ring = [&editor, face](std::span<const Point3d> ring,
                                        bool outer_ring) {
    std::vector<VertexId> vertices;
    vertices.reserve(ring.size());
    for (const Point3d& p : ring) {
      vertices.push_back(editor.add_vertex(p, 0.0).value());
    }
    std::vector<EdgeId> edges;
    edges.reserve(ring.size());
    for (std::size_t i = 0; i < ring.size(); ++i) {
      edges.push_back(editor
                          .add_edge(vertices[i], vertices[(i + 1) % ring.size()],
                                    kNoGeometry, 0.0, 1.0, 0.0)
                          .value());
    }
    std::vector<CoedgeSpec> specs(ring.size());
    for (std::size_t i = 0; i < specs.size(); ++i) specs[i].edge = edges[i];
    editor.add_loop(face, specs, outer_ring).value();
  };
  add_ring(outer, true);
  for (const std::vector<Point3d>& hole : holes) add_ring(hole, false);
  const std::array<FaceId, 1> faces{face};
  editor.add_shell(faces, false).value();
  return editor.commit().value();
}
}  // namespace

TEST_CASE("planar_cut: one transverse cut splits a square into two fragments") {
  const Brep box = make_box(Point3d{0, 0, 0}, Point3d{2, 2, 2}, 1e-12).value();
  const FaceId top = plane_face(box, Vec3d{0, 0, 1}, 2.0);
  const Segments cuts{{Point3d{1, 0, 2}, Point3d{1, 2, 2}}};
  const auto fragments = cut_face_by_segments(box, top, cuts, 1e-12);
  REQUIRE(fragments.has_value());
  REQUIRE(fragments.value().size() == 2);
  double total = 0.0;
  for (const FaceFragment2d& f : fragments.value()) {
    CHECK(f.holes.empty());
    total += fragment_area(f);
  }
  CHECK(total == doctest::Approx(4.0).epsilon(1e-12));
}

TEST_CASE("planar_cut: no cuts returns the whole face") {
  const Brep box = make_box(Point3d{0, 0, 0}, Point3d{2, 2, 2}, 1e-12).value();
  const FaceId top = plane_face(box, Vec3d{0, 0, 1}, 2.0);
  const Segments none;
  const auto fragments = cut_face_by_segments(box, top, none, 1e-12);
  REQUIRE(fragments.has_value());
  REQUIRE(fragments.value().size() == 1);
  CHECK(fragments.value()[0].outer.size() == 4);
  CHECK(fragment_area(fragments.value()[0]) == doctest::Approx(4.0).epsilon(1e-12));
}

TEST_CASE("planar_cut: a cut that misses the face returns the whole face") {
  const Brep box = make_box(Point3d{0, 0, 0}, Point3d{2, 2, 2}, 1e-12).value();
  const FaceId top = plane_face(box, Vec3d{0, 0, 1}, 2.0);
  const Segments outside{{Point3d{5, -1, 2}, Point3d{5, 3, 2}}};
  const auto fragments = cut_face_by_segments(box, top, outside, 1e-12);
  REQUIRE(fragments.has_value());
  REQUIRE(fragments.value().size() == 1);
  CHECK(fragment_area(fragments.value()[0]) == doctest::Approx(4.0).epsilon(1e-12));
}

TEST_CASE("planar_cut: crossing cuts make four fragments whose areas sum to four") {
  const Brep box = make_box(Point3d{0, 0, 0}, Point3d{2, 2, 2}, 1e-12).value();
  const FaceId top = plane_face(box, Vec3d{0, 0, 1}, 2.0);
  const Segments cuts{{Point3d{1, 0, 2}, Point3d{1, 2, 2}},
                      {Point3d{0, 1, 2}, Point3d{2, 1, 2}}};
  const auto fragments = cut_face_by_segments(box, top, cuts, 1e-12);
  REQUIRE(fragments.has_value());
  REQUIRE(fragments.value().size() == 4);
  double total = 0.0;
  for (const FaceFragment2d& f : fragments.value()) total += fragment_area(f);
  CHECK(total == doctest::Approx(4.0).epsilon(1e-12));
}

TEST_CASE("planar_cut: a closed cut loop makes an island and a holed remainder") {
  const Brep box = make_box(Point3d{0, 0, 0}, Point3d{2, 2, 2}, 1e-12).value();
  const FaceId top = plane_face(box, Vec3d{0, 0, 1}, 2.0);
  const Segments loop{{Point3d{1.0, 1.0, 2}, Point3d{1.5, 1.0, 2}},
                      {Point3d{1.5, 1.0, 2}, Point3d{1.5, 1.5, 2}},
                      {Point3d{1.5, 1.5, 2}, Point3d{1.0, 1.5, 2}},
                      {Point3d{1.0, 1.5, 2}, Point3d{1.0, 1.0, 2}}};
  const auto fragments = cut_face_by_segments(box, top, loop, 1e-12);
  REQUIRE(fragments.has_value());
  REQUIRE(fragments.value().size() == 2);
  double total = 0.0;
  double island = 0.0;
  std::size_t holed = 0;
  for (const FaceFragment2d& f : fragments.value()) {
    total += fragment_area(f);
    if (f.holes.empty()) {
      island += fragment_area(f);
    } else {
      ++holed;
      CHECK(f.holes.size() == 1);
    }
  }
  CHECK(holed == 1);
  CHECK(island == doctest::Approx(0.25).epsilon(1e-12));
  CHECK(total == doctest::Approx(4.0).epsilon(1e-12));
}

TEST_CASE("planar_cut: a dangling cut is rejected") {
  const Brep box = make_box(Point3d{0, 0, 0}, Point3d{2, 2, 2}, 1e-12).value();
  const FaceId top = plane_face(box, Vec3d{0, 0, 1}, 2.0);
  const Segments dangling{{Point3d{1, 0, 2}, Point3d{1, 1, 2}}};
  const auto fragments = cut_face_by_segments(box, top, dangling, 1e-12);
  REQUIRE_FALSE(fragments.has_value());
  CHECK(fragments.error().code == Errc::invalid_topology);
}

TEST_CASE("planar_cut: fragment samples lie strictly inside the fragment") {
  const Brep box = make_box(Point3d{0, 0, 0}, Point3d{2, 2, 2}, 1e-12).value();
  const FaceId top = plane_face(box, Vec3d{0, 0, 1}, 2.0);
  REQUIRE_FALSE(top.is_null());

  const Segments crossing{{Point3d{1, 0, 2}, Point3d{1, 2, 2}},
                          {Point3d{0, 1, 2}, Point3d{2, 1, 2}}};
  const auto square = cut_face_by_segments(box, top, crossing, 1e-12);
  REQUIRE(square.has_value());
  REQUIRE(square.value().size() == 4);
  check_samples_inside(box, top, square.value(), 1e-12);

  const Segments loop{{Point3d{1.0, 1.0, 2}, Point3d{1.5, 1.0, 2}},
                      {Point3d{1.5, 1.0, 2}, Point3d{1.5, 1.5, 2}},
                      {Point3d{1.5, 1.5, 2}, Point3d{1.0, 1.5, 2}},
                      {Point3d{1.0, 1.5, 2}, Point3d{1.0, 1.0, 2}}};
  const auto island = cut_face_by_segments(box, top, loop, 1e-12);
  REQUIRE(island.has_value());
  REQUIRE(island.value().size() == 2);
  check_samples_inside(box, top, island.value(), 1e-12);
}

TEST_CASE("planar_cut: nested cut loops are separate holes of nested owners") {
  const Brep box = make_box(Point3d{0, 0, 0}, Point3d{2, 2, 2}, 1e-12).value();
  const FaceId top = plane_face(box, Vec3d{0, 0, 1}, 2.0);
  const Segments loops{
      {Point3d{0.5, 0.5, 2}, Point3d{1.5, 0.5, 2}},
      {Point3d{1.5, 0.5, 2}, Point3d{1.5, 1.5, 2}},
      {Point3d{1.5, 1.5, 2}, Point3d{0.5, 1.5, 2}},
      {Point3d{0.5, 1.5, 2}, Point3d{0.5, 0.5, 2}},
      {Point3d{0.75, 0.75, 2}, Point3d{1.25, 0.75, 2}},
      {Point3d{1.25, 0.75, 2}, Point3d{1.25, 1.25, 2}},
      {Point3d{1.25, 1.25, 2}, Point3d{0.75, 1.25, 2}},
      {Point3d{0.75, 1.25, 2}, Point3d{0.75, 0.75, 2}}};
  const auto fragments = cut_face_by_segments(box, top, loops, 1e-12);
  REQUIRE(fragments.has_value());
  const auto& f = fragments.value();
  REQUIRE(f.size() == 3);
  CHECK(fragment_area(f[0]) == doctest::Approx(3.0).epsilon(1e-12));
  CHECK(fragment_area(f[1]) == doctest::Approx(0.75).epsilon(1e-12));
  CHECK(fragment_area(f[2]) == doctest::Approx(0.25).epsilon(1e-12));
  CHECK(fragment_area(f[0]) > fragment_area(f[1]));
  std::size_t holed = 0;
  double total = 0.0;
  for (const FaceFragment2d& frag : f) {
    total += fragment_area(frag);
    if (!frag.holes.empty()) {
      ++holed;
      CHECK(frag.holes.size() == 1);
    }
  }
  CHECK(holed == 2);
  CHECK(total == doctest::Approx(4.0).epsilon(1e-12));
  REQUIRE(f[0].holes.size() == 1);
  CHECK(std::fabs(shoelace(f[0].holes[0])) ==
        doctest::Approx(1.0).epsilon(1e-12));
  REQUIRE(f[1].holes.size() == 1);
  CHECK(std::fabs(shoelace(f[1].holes[0])) ==
        doctest::Approx(0.25).epsilon(1e-12));
  CHECK(f[2].holes.empty());
}

TEST_CASE("planar_cut: cut endpoints near the boundary snap onto it") {
  const Brep box = make_box(Point3d{0, 0, 0}, Point3d{2, 2, 2}, 1e-12).value();
  const FaceId top = plane_face(box, Vec3d{0, 0, 1}, 2.0);
  const double tol = 1e-6;
  const std::array<Point3d, 2> starts{Point3d{1.0, 1e-15, 2.0},
                                      Point3d{1.0, -1e-13, 2.0}};
  for (const Point3d& start : starts) {
    const Segments cut{{start, Point3d{1.0, 2.0 - 1e-15, 2.0}}};
    const auto fragments = cut_face_by_segments(box, top, cut, tol);
    REQUIRE(fragments.has_value());
    REQUIRE(fragments.value().size() == 2);
    CHECK(fragment_area(fragments.value()[0]) ==
          doctest::Approx(2.0).epsilon(1e-12));
    CHECK(fragment_area(fragments.value()[1]) ==
          doctest::Approx(2.0).epsilon(1e-12));
    for (const FaceFragment2d& frag : fragments.value()) {
      CHECK(frag.holes.empty());
    }
  }
}

TEST_CASE("planar_cut: a cut loop outside the face is filtered out") {
  const Brep box = make_box(Point3d{0, 0, 0}, Point3d{2, 2, 2}, 1e-12).value();
  const FaceId top = plane_face(box, Vec3d{0, 0, 1}, 2.0);
  const Segments loop{{Point3d{3.0, 3.0, 2}, Point3d{3.5, 3.0, 2}},
                      {Point3d{3.5, 3.0, 2}, Point3d{3.5, 3.5, 2}},
                      {Point3d{3.5, 3.5, 2}, Point3d{3.0, 3.5, 2}},
                      {Point3d{3.0, 3.5, 2}, Point3d{3.0, 3.0, 2}}};
  const auto fragments = cut_face_by_segments(box, top, loop, 1e-12);
  REQUIRE(fragments.has_value());
  REQUIRE(fragments.value().size() == 1);
  CHECK(fragment_area(fragments.value()[0]) ==
        doctest::Approx(4.0).epsilon(1e-12));
  CHECK(fragments.value()[0].holes.empty());
}

TEST_CASE("planar_cut: a keyhole cut yields one slit fragment") {
  const std::array<Point3d, 4> outer{Point3d{0, 0, 2}, Point3d{4, 0, 2},
                                     Point3d{4, 4, 2}, Point3d{0, 4, 2}};
  const std::array<Point3d, 4> hole{Point3d{1, 1, 2}, Point3d{1, 2, 2},
                                    Point3d{2, 2, 2}, Point3d{2, 1, 2}};
  const std::vector<std::vector<Point3d>> holes{
      std::vector<Point3d>(hole.begin(), hole.end())};
  const Brep face = make_planar_face(outer, holes, Vec3d{0, 0, 1}, 2.0);
  const FaceId id = face.face_ids().front();
  const Segments cut{{Point3d{4, 1.5, 2}, Point3d{2, 1.5, 2}}};
  const auto fragments = cut_face_by_segments(face, id, cut, 1e-12);
  REQUIRE(fragments.has_value());
  REQUIRE(fragments.value().size() == 1);
  const FaceFragment2d& f = fragments.value()[0];
  CHECK(f.outer.size() == 12);
  CHECK(f.holes.empty());
  CHECK(fragment_area(f) == doctest::Approx(15.0).epsilon(1e-12));
  check_samples_inside(face, id, fragments.value(), 1e-12);
}

TEST_CASE("planar_cut: fragments are ordered by area then lexicographically") {
  const Brep box = make_box(Point3d{0, 0, 0}, Point3d{2, 2, 2}, 1e-12).value();
  const FaceId top = plane_face(box, Vec3d{0, 0, 1}, 2.0);
  const Segments loops{
      {Point3d{0.5, 0.5, 2}, Point3d{1.0, 0.5, 2}},
      {Point3d{1.0, 0.5, 2}, Point3d{1.0, 1.0, 2}},
      {Point3d{1.0, 1.0, 2}, Point3d{0.5, 1.0, 2}},
      {Point3d{0.5, 1.0, 2}, Point3d{0.5, 0.5, 2}},
      {Point3d{0.5, 1.2, 2}, Point3d{1.0, 1.2, 2}},
      {Point3d{1.0, 1.2, 2}, Point3d{1.0, 1.7, 2}},
      {Point3d{1.0, 1.7, 2}, Point3d{0.5, 1.7, 2}},
      {Point3d{0.5, 1.7, 2}, Point3d{0.5, 1.2, 2}}};
  const auto fragments = cut_face_by_segments(box, top, loops, 1e-12);
  REQUIRE(fragments.has_value());
  const auto& f = fragments.value();
  REQUIRE(f.size() == 3);
  CHECK(fragment_area(f[0]) == doctest::Approx(3.5).epsilon(1e-12));
  CHECK(fragment_area(f[0]) >= fragment_area(f[1]));
  CHECK(fragment_area(f[0]) >= fragment_area(f[2]));
  CHECK(fragment_area(f[1]) == doctest::Approx(0.25).epsilon(1e-12));
  CHECK(fragment_area(f[2]) == doctest::Approx(0.25).epsilon(1e-12));
  CHECK(f[0].holes.size() == 2);
  // Equal-area islands fall back to the lexicographic order of their first
  // point. In the face frame of the +z top face, the island at world y 1.2..1.7
  // maps to smaller frame x than the island at world y 0.5..1.0.
  const Point2d a = f[1].outer.front();
  const Point2d b = f[2].outer.front();
  CHECK(a[0] < b[0]);
}
