#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <initializer_list>
#include <limits>
#include <numbers>
#include <string_view>
#include <variant>
#include <vector>

#include <doctest/doctest.h>
#include <geom/brep/face_query.hpp>
#include <geom/brep/face_uv.hpp>
#include <geom/brep/plane_frame.hpp>
#include <geom/brep/primitives.hpp>
#include <geom/brep/solid_editor.hpp>
#include <geom/brep/topology.hpp>
#include <geom/brep/trim_face.hpp>
#include <geom/cgeom/polygon2.hpp>
#include <geom/core/error.hpp>
#include <geom/curves/analytic_curve.hpp>
#include <geom/primitives/point.hpp>
#include <geom/primitives/vec.hpp>
#include <geom/surfaces/analytic.hpp>
#include <geom/surfaces/plane.hpp>

using namespace geom;

namespace {

constexpr double kTol = 1e-12;
constexpr double kPi = std::numbers::pi_v<double>;

template <class Surface>
FaceId face_with_surface(const Brep& brep) {
  for (const FaceId id : brep.face_ids()) {
    const Face& face = brep.face(id);
    if (face.surface == kNoGeometry) continue;
    if (std::holds_alternative<Surface>(brep.geometry().surfaces[face.surface])) {
      return id;
    }
  }
  return FaceId{};
}

// The face whose plane has the given outward normal and offset magnitude.
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
  return FaceId{};
}

double ring_area(std::span<const Point2d> ring) {
  return std::fabs(signed_area2(ring));
}

double fragment_area(const CurvedFragment& fragment) {
  REQUIRE(!fragment.loops_uv.empty());
  double area = ring_area(fragment.loops_uv.front());
  for (std::size_t i = 1; i < fragment.loops_uv.size(); ++i) {
    area -= ring_area(fragment.loops_uv[i]);
  }
  return area;
}

// Converts world points on a planar face into the face's uv frame.
FaceCut uv_cut(const Brep& brep, FaceId face,
               std::initializer_list<Point3d> points) {
  const Result<Plane3> plane = face_plane(brep, face);
  REQUIRE(plane.has_value());
  const Result<PlaneFrame> frame = make_plane_frame(plane.value());
  REQUIRE(frame.has_value());
  FaceCut cut;
  cut.uv.reserve(points.size());
  for (const Point3d& p : points) {
    cut.uv.push_back(plane_frame_to_2d(frame.value(), p));
  }
  return cut;
}

// Builds a single planar face from an outer ring and hole rings. The outer ring
// must wind counter-clockwise around normal and each hole ring clockwise, both
// in world space.
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
      const Point3d& from = ring[i];
      const Point3d& to = ring[(i + 1) % ring.size()];
      const std::uint32_t curve = editor.add_curve(
          LineCurve3{from, Vec3d{to[0] - from[0], to[1] - from[1],
                                 to[2] - from[2]}});
      edges.push_back(editor
                          .add_edge(vertices[i], vertices[(i + 1) % ring.size()],
                                    curve, 0.0, 1.0, 0.0)
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

// Checks that the uv rings are the cycle rings (outer counter-clockwise, holes
// clockwise) and that every 3D vertex is the planar lift of its uv vertex.
void check_plane_fragment(const Brep& brep, FaceId face,
                          const CurvedFragment& fragment, double extent) {
  CHECK(fragment.reversed == false);
  REQUIRE(!fragment.loops_uv.empty());
  CHECK(fragment.loops.size() == fragment.loops_uv.size());
  for (std::size_t r = 0; r < fragment.loops_uv.size(); ++r) {
    CHECK(fragment.loops[r].size() == fragment.loops_uv[r].size());
    if (r == 0) {
      CHECK(signed_area2(fragment.loops_uv[r]) > 0.0);
    } else {
      CHECK(signed_area2(fragment.loops_uv[r]) < 0.0);
    }
  }

  const Result<Plane3> plane = face_plane(brep, face);
  REQUIRE(plane.has_value());
  const Result<PlaneFrame> frame = make_plane_frame(plane.value());
  REQUIRE(frame.has_value());
  const double bound = 1e-12 * extent;
  for (std::size_t r = 0; r < fragment.loops_uv.size(); ++r) {
    for (std::size_t i = 0; i < fragment.loops_uv[r].size(); ++i) {
      const Point3d expected =
          plane_frame_to_3d(frame.value(), fragment.loops_uv[r][i]);
      CHECK(norm2(fragment.loops[r][i] - expected) <= bound * bound);
    }
  }
}

}  // namespace

TEST_CASE("trim_face: one chord splits a plane face into two fragments") {
  const Brep box = make_box(Point3d{0, 0, 0}, Point3d{2, 2, 2}, kTol).value();
  const FaceId top = plane_face(box, Vec3d{0, 0, 1}, 2.0);
  REQUIRE_FALSE(top.is_null());

  const std::vector<FaceCut> cuts{
      uv_cut(box, top, {Point3d{1, 0, 2}, Point3d{1, 2, 2}})};
  const Result<std::vector<CurvedFragment>> fragments =
      split_face_curved(box, top, cuts, kTol);
  REQUIRE(fragments.has_value());
  REQUIRE(fragments.value().size() == 2);

  double total = 0.0;
  for (const CurvedFragment& fragment : fragments.value()) {
    CHECK(fragment.loops_uv.size() == 1);
    CHECK(fragment_area(fragment) == doctest::Approx(2.0).epsilon(1e-12));
    total += fragment_area(fragment);
    check_plane_fragment(box, top, fragment, 2.0);
  }
  CHECK(total == doctest::Approx(4.0).epsilon(1e-12));

  // A second identical call must produce the same rings.
  const Result<std::vector<CurvedFragment>> again =
      split_face_curved(box, top, cuts, kTol);
  REQUIRE(again.has_value());
  REQUIRE(again.value().size() == fragments.value().size());
  for (std::size_t f = 0; f < fragments.value().size(); ++f) {
    REQUIRE(again.value()[f].loops_uv.size() ==
            fragments.value()[f].loops_uv.size());
    for (std::size_t r = 0; r < fragments.value()[f].loops_uv.size(); ++r) {
      REQUIRE(again.value()[f].loops_uv[r].size() ==
              fragments.value()[f].loops_uv[r].size());
      for (std::size_t i = 0; i < fragments.value()[f].loops_uv[r].size(); ++i) {
        CHECK(again.value()[f].loops_uv[r][i][0] ==
              fragments.value()[f].loops_uv[r][i][0]);
        CHECK(again.value()[f].loops_uv[r][i][1] ==
              fragments.value()[f].loops_uv[r][i][1]);
      }
    }
  }
}

TEST_CASE("trim_face: no cuts returns the whole plane face") {
  const Brep box = make_box(Point3d{0, 0, 0}, Point3d{2, 2, 2}, kTol).value();
  const FaceId top = plane_face(box, Vec3d{0, 0, 1}, 2.0);
  const std::vector<FaceCut> none;
  const Result<std::vector<CurvedFragment>> fragments =
      split_face_curved(box, top, none, kTol);
  REQUIRE(fragments.has_value());
  REQUIRE(fragments.value().size() == 1);
  CHECK(fragments.value()[0].loops_uv.size() == 1);
  CHECK(fragments.value()[0].loops_uv[0].size() == 4);
  CHECK(fragment_area(fragments.value()[0]) ==
        doctest::Approx(4.0).epsilon(1e-12));
  check_plane_fragment(box, top, fragments.value()[0], 2.0);
}

TEST_CASE("trim_face: a full cylinder lateral face is split by an axial cut") {
  const Brep cylinder =
      make_cylinder(Point3d{0, 0, 0}, Vec3d{0, 0, 1}, 2.0, 3.0, kTol).value();
  const FaceId lateral = face_with_surface<Cylinder3>(cylinder);
  REQUIRE_FALSE(lateral.is_null());

  const Result<std::vector<FaceUvLoop>> loops =
      face_loops_uv(cylinder, lateral, kTol);
  REQUIRE(loops.has_value());
  REQUIRE(loops.value().size() == 1);
  double min_u = std::numeric_limits<double>::infinity();
  double max_u = -std::numeric_limits<double>::infinity();
  double min_v = std::numeric_limits<double>::infinity();
  double max_v = -std::numeric_limits<double>::infinity();
  for (const Point2d& p : loops.value()[0].points) {
    min_u = std::min(min_u, p[0]);
    max_u = std::max(max_u, p[0]);
    min_v = std::min(min_v, p[1]);
    max_v = std::max(max_v, p[1]);
  }
  CHECK(max_u - min_u == doctest::Approx(2.0 * kPi).epsilon(1e-12));

  const double u_mid = 0.5 * (min_u + max_u);
  FaceCut cut;
  cut.uv = {{u_mid, min_v}, {u_mid, max_v}};
  const std::vector<FaceCut> cuts{cut};
  const Result<std::vector<CurvedFragment>> fragments =
      split_face_curved(cylinder, lateral, cuts, kTol);
  REQUIRE(fragments.has_value());
  REQUIRE(fragments.value().size() == 2);

  const auto& surface =
      std::get<Cylinder3>(cylinder.geometry().surfaces[cylinder.face(lateral).surface]);
  double total = 0.0;
  for (const CurvedFragment& fragment : fragments.value()) {
    CHECK(fragment.reversed == false);
    CHECK(fragment.loops_uv.size() == 1);
    CHECK(fragment.loops.size() == fragment.loops_uv.size());
    CHECK(ring_area(fragment.loops_uv[0]) ==
          doctest::Approx(kPi * 3.0).epsilon(1e-12));
    total += fragment_area(fragment);
    REQUIRE(fragment.loops[0].size() == fragment.loops_uv[0].size());
    for (std::size_t i = 0; i < fragment.loops_uv[0].size(); ++i) {
      const Point2d& uv = fragment.loops_uv[0][i];
      const Point3d& point = fragment.loops[0][i];
      CHECK(std::fabs(std::hypot(point[0], point[1]) - 2.0) <= 1e-12 * 3.0);
      CHECK(point[2] >= -1e-12 * 3.0);
      CHECK(point[2] <= 3.0 + 1e-12 * 3.0);
      const Point3d expected = surface.eval(uv[0], uv[1]);
      CHECK(norm2(point - expected) <= (1e-12 * 3.0) * (1e-12 * 3.0));
    }
  }
  CHECK(total == doctest::Approx(2.0 * kPi * 3.0).epsilon(1e-12));
}

TEST_CASE("trim_face: a cut from the outer boundary to a hole opens the hole") {
  const std::array<Point3d, 4> outer{Point3d{0, 0, 0}, Point3d{4, 0, 0},
                                     Point3d{4, 4, 0}, Point3d{0, 4, 0}};
  const std::array<Point3d, 4> hole{Point3d{1, 1, 0}, Point3d{1, 2, 0},
                                    Point3d{2, 2, 0}, Point3d{2, 1, 0}};
  const std::vector<std::vector<Point3d>> holes{
      std::vector<Point3d>(hole.begin(), hole.end())};
  const Brep face = make_planar_face(outer, holes, Vec3d{0, 0, 1}, 0.0);
  const FaceId id = face.face_ids().front();

  const std::vector<FaceCut> cuts{
      uv_cut(face, id, {Point3d{4, 1.5, 0}, Point3d{2, 1.5, 0}})};
  const Result<std::vector<CurvedFragment>> fragments =
      split_face_curved(face, id, cuts, kTol);
  REQUIRE(fragments.has_value());
  REQUIRE(fragments.value().size() == 1);
  const CurvedFragment& fragment = fragments.value()[0];
  CHECK(fragment.loops_uv.size() == 1);
  CHECK(fragment.loops_uv[0].size() == 12);
  CHECK(fragment_area(fragment) == doctest::Approx(15.0).epsilon(1e-12));
  check_plane_fragment(face, id, fragment, 4.0);
}

TEST_CASE("trim_face: a closed cut loop makes an interior fragment") {
  const Brep box = make_box(Point3d{0, 0, 0}, Point3d{2, 2, 2}, kTol).value();
  const FaceId top = plane_face(box, Vec3d{0, 0, 1}, 2.0);
  REQUIRE_FALSE(top.is_null());

  const std::vector<FaceCut> cuts{uv_cut(
      box, top,
      {Point3d{0.5, 0.5, 2}, Point3d{1.5, 0.5, 2}, Point3d{1.5, 1.5, 2},
       Point3d{0.5, 1.5, 2}, Point3d{0.5, 0.5, 2}})};
  const Result<std::vector<CurvedFragment>> fragments =
      split_face_curved(box, top, cuts, kTol);
  REQUIRE(fragments.has_value());
  REQUIRE(fragments.value().size() == 2);

  double total = 0.0;
  double island = -1.0;
  double holed = -1.0;
  for (const CurvedFragment& fragment : fragments.value()) {
    total += fragment_area(fragment);
    if (fragment.loops_uv.size() == 1) {
      island = fragment_area(fragment);
    } else {
      REQUIRE(fragment.loops_uv.size() == 2);
      holed = fragment_area(fragment);
      CHECK(ring_area(fragment.loops_uv[1]) ==
            doctest::Approx(1.0).epsilon(1e-12));
    }
    check_plane_fragment(box, top, fragment, 2.0);
  }
  CHECK(island == doctest::Approx(1.0).epsilon(1e-12));
  CHECK(holed == doctest::Approx(3.0).epsilon(1e-12));
  CHECK(total == doctest::Approx(4.0).epsilon(1e-12));
}

TEST_CASE("trim_face: a closed cut inside a hole stays void") {
  const std::array<Point3d, 4> outer{Point3d{0, 0, 0}, Point3d{4, 0, 0},
                                     Point3d{4, 4, 0}, Point3d{0, 4, 0}};
  const std::array<Point3d, 4> hole{Point3d{1, 1, 0}, Point3d{1, 2, 0},
                                    Point3d{2, 2, 0}, Point3d{2, 1, 0}};
  const std::vector<std::vector<Point3d>> holes{
      std::vector<Point3d>(hole.begin(), hole.end())};
  const Brep face = make_planar_face(outer, holes, Vec3d{0, 0, 1}, 0.0);
  const FaceId id = face.face_ids().front();

  const std::vector<FaceCut> cuts{uv_cut(
      face, id,
      {Point3d{1.25, 1.25, 0}, Point3d{1.75, 1.25, 0}, Point3d{1.75, 1.75, 0},
       Point3d{1.25, 1.75, 0}, Point3d{1.25, 1.25, 0}})};
  const Result<std::vector<CurvedFragment>> fragments =
      split_face_curved(face, id, cuts, kTol);
  REQUIRE(fragments.has_value());
  // The hole is void, so a loop carved inside it encloses no material.
  REQUIRE(fragments.value().size() == 1);
  CHECK(fragments.value()[0].loops_uv.size() == 2);
  CHECK(fragment_area(fragments.value()[0]) ==
        doctest::Approx(15.0).epsilon(1e-12));
  check_plane_fragment(face, id, fragments.value()[0], 4.0);
}

TEST_CASE("trim_face: a grazing cut drops a below-tolerance sliver") {
  const Brep box = make_box(Point3d{0, 0, 0}, Point3d{2, 2, 2}, kTol).value();
  const FaceId top = plane_face(box, Vec3d{0, 0, 1}, 2.0);
  REQUIRE_FALSE(top.is_null());

  const double tolerance = 1e-3;
  const std::vector<FaceCut> cuts{uv_cut(
      box, top, {Point3d{0, 2.0 - 1e-7, 2}, Point3d{2, 2.0 - 1e-7, 2}})};
  const Result<std::vector<CurvedFragment>> fragments =
      split_face_curved(box, top, cuts, tolerance);
  REQUIRE(fragments.has_value());
  REQUIRE(fragments.value().size() == 1);
  CHECK(fragments.value()[0].loops_uv.size() == 1);
  CHECK(fragment_area(fragments.value()[0]) ==
        doctest::Approx(4.0).epsilon(1e-6));
}

TEST_CASE("trim_face: an open cut endpoint off the boundary is rejected") {
  const Brep box = make_box(Point3d{0, 0, 0}, Point3d{2, 2, 2}, kTol).value();
  const FaceId top = plane_face(box, Vec3d{0, 0, 1}, 2.0);
  REQUIRE_FALSE(top.is_null());

  const std::vector<FaceCut> cuts{
      uv_cut(box, top, {Point3d{0.5, 0.6, 2}, Point3d{1.2, 1.1, 2}})};
  const Result<std::vector<CurvedFragment>> fragments =
      split_face_curved(box, top, cuts, kTol);
  REQUIRE_FALSE(fragments.has_value());
  CHECK(fragments.error().code == Errc::numeric_failure);
  CHECK(std::string_view(fragments.error().message)
            .starts_with("split_face_curved:"));
}

TEST_CASE("trim_face: a cone lateral face propagates unsupported") {
  const Brep cone =
      make_cone(Point3d{0, 0, 0}, Vec3d{0, 0, 1}, 2.0, 3.0, kTol).value();
  const FaceId lateral = face_with_surface<Cone3>(cone);
  REQUIRE_FALSE(lateral.is_null());

  const std::vector<FaceCut> none;
  const Result<std::vector<CurvedFragment>> fragments =
      split_face_curved(cone, lateral, none, kTol);
  REQUIRE_FALSE(fragments.has_value());
  CHECK(fragments.error().code == Errc::unsupported);
  CHECK(std::string_view(fragments.error().message).starts_with("face_loops_uv:"));
}
