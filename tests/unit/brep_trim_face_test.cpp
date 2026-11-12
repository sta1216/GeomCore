#include <algorithm>
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
#include <geom/brep/tessellate.hpp>
#include <geom/brep/topology.hpp>
#include <geom/brep/trim_face.hpp>
#include <geom/brep/validate.hpp>
#include <geom/cgeom/polygon2.hpp>
#include <geom/core/error.hpp>
#include <geom/curves/analytic_curve.hpp>
#include <geom/primitives/point.hpp>
#include <geom/primitives/vec.hpp>
#include <geom/surfaces/analytic.hpp>
#include <geom/surfaces/plane.hpp>

#include "brep_curved_fixture.hpp"

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

// The plane face whose uv window is exactly extent x extent: the face of a box
// built from the origin to (extent, extent, extent) that the junction fixtures
// below cut in a 10 x 10 parameter plane.
FaceId box_plane_face(const Brep& brep, double extent) {
  for (const FaceId id : brep.face_ids()) {
    const auto& variant = brep.geometry().surfaces[brep.face(id).surface];
    if (!std::holds_alternative<Plane3>(variant)) continue;
    const Result<std::vector<FaceUvLoop>> uv = face_loops_uv(brep, id, 1e-9);
    if (!uv.has_value() || uv.value().empty()) continue;
    double umin = 1e300;
    double umax = -1e300;
    double vmin = 1e300;
    double vmax = -1e300;
    for (const Point2d& p : uv.value().front().points) {
      umin = std::min(umin, p[0]);
      umax = std::max(umax, p[0]);
      vmin = std::min(vmin, p[1]);
      vmax = std::max(vmax, p[1]);
    }
    if (std::fabs(umax - umin - extent) <= 1e-9 &&
        std::fabs(vmax - vmin - extent) <= 1e-9) {
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
  // 2e-7 from the edge: close enough to snap, far enough that the sliver still
  // triangulates, so the area floor (not a triangulation failure) drops it.
  const std::vector<FaceCut> cuts{uv_cut(
      box, top, {Point3d{0, 2.0 - 2e-7, 2}, Point3d{2, 2.0 - 2e-7, 2}})};
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

TEST_CASE("trim_face: two collinear cuts weld at their shared endpoint") {
  const Brep box = make_box(Point3d{0, 0, 0}, Point3d{2, 2, 2}, kTol).value();
  const FaceId top = plane_face(box, Vec3d{0, 0, 1}, 2.0);
  REQUIRE_FALSE(top.is_null());

  const std::vector<FaceCut> cuts{
      uv_cut(box, top, {Point3d{0, 0.5, 2}, Point3d{1, 0.5, 2}}),
      uv_cut(box, top, {Point3d{1, 0.5, 2}, Point3d{2, 0.5, 2}})};
  const Result<std::vector<CurvedFragment>> fragments =
      split_face_curved(box, top, cuts, kTol);
  REQUIRE(fragments.has_value());
  REQUIRE(fragments.value().size() == 2);
  double total = 0.0;
  std::vector<double> areas;
  for (const CurvedFragment& fragment : fragments.value()) {
    CHECK(fragment.loops_uv.size() == 1);
    areas.push_back(fragment_area(fragment));
    total += fragment_area(fragment);
    check_plane_fragment(box, top, fragment, 2.0);
  }
  // The chord at y = 0.5 cuts the 2 x 2 square into 3 and 1, not 2 and 2 as the
  // review message stated.
  CHECK(total == doctest::Approx(4.0).epsilon(1e-12));
  std::sort(areas.begin(), areas.end());
  CHECK(areas[0] == doctest::Approx(1.0).epsilon(1e-12));
  CHECK(areas[1] == doctest::Approx(3.0).epsilon(1e-12));
}

TEST_CASE("trim_face: a cut T-junction onto another cut splits both regions") {
  const Brep box = make_box(Point3d{0, 0, 0}, Point3d{2, 2, 2}, kTol).value();
  const FaceId top = plane_face(box, Vec3d{0, 0, 1}, 2.0);
  REQUIRE_FALSE(top.is_null());

  // The second cut stops 1e-9 short of the first one, so the endpoint is only
  // accepted because it T-junctions into the host cut and is projected into it.
  const double tolerance = 1e-6;
  const std::vector<FaceCut> cuts{
      uv_cut(box, top, {Point3d{0, 0.25, 2}, Point3d{2, 0.25, 2}}),
      uv_cut(box, top,
             {Point3d{1, 0.25 + 1e-9, 2}, Point3d{1, 2.0, 2}})};
  const Result<std::vector<CurvedFragment>> fragments =
      split_face_curved(box, top, cuts, tolerance);
  REQUIRE(fragments.has_value());
  REQUIRE(fragments.value().size() == 3);
  double total = 0.0;
  std::vector<double> areas;
  for (const CurvedFragment& fragment : fragments.value()) {
    CHECK(fragment.loops_uv.size() == 1);
    areas.push_back(fragment_area(fragment));
    total += fragment_area(fragment);
    check_plane_fragment(box, top, fragment, 2.0);
  }
  CHECK(total == doctest::Approx(4.0).epsilon(1e-12));
  std::sort(areas.begin(), areas.end());
  CHECK(areas[0] == doctest::Approx(0.5).epsilon(1e-12));
  CHECK(areas[1] == doctest::Approx(1.75).epsilon(1e-12));
  CHECK(areas[2] == doctest::Approx(1.75).epsilon(1e-12));
}

TEST_CASE("trim_face: a tiny closed island below tolerance is dropped") {
  const Brep box = make_box(Point3d{0, 0, 0}, Point3d{2, 2, 2}, kTol).value();
  const FaceId top = plane_face(box, Vec3d{0, 0, 1}, 2.0);
  REQUIRE_FALSE(top.is_null());

  const double tolerance = 1e-3;
  const double half = 5e-5;
  const std::vector<FaceCut> cuts{uv_cut(
      box, top,
      {Point3d{1 - half, 1 - half, 2}, Point3d{1 + half, 1 - half, 2},
       Point3d{1 + half, 1 + half, 2}, Point3d{1 - half, 1 + half, 2},
       Point3d{1 - half, 1 - half, 2}})};
  const Result<std::vector<CurvedFragment>> fragments =
      split_face_curved(box, top, cuts, tolerance);
  REQUIRE(fragments.has_value());
  // The 1e-8 island falls under the tolerance^2 area floor: only the punched
  // remainder is emitted, with the island's ring as its hole.
  REQUIRE(fragments.value().size() == 1);
  CHECK(fragments.value()[0].loops_uv.size() == 2);
  CHECK(fragment_area(fragments.value()[0]) ==
        doctest::Approx(4.0).epsilon(1e-6));
}

TEST_CASE("trim_face: a thin frame below the area floor yields no fragments") {
  const std::array<Point3d, 4> outer{Point3d{0, 0, 0}, Point3d{2, 0, 0},
                                     Point3d{2, 2, 0}, Point3d{0, 2, 0}};
  const std::array<Point3d, 4> hole{Point3d{0.025, 0.025, 0},
                                    Point3d{0.025, 1.975, 0},
                                    Point3d{1.975, 1.975, 0},
                                    Point3d{1.975, 0.025, 0}};
  const std::vector<std::vector<Point3d>> holes{
      std::vector<Point3d>(hole.begin(), hole.end())};
  const Brep face = make_planar_face(outer, holes, Vec3d{0, 0, 1}, 0.0);
  const FaceId id = face.face_ids().front();

  const std::vector<FaceCut> none;
  const Result<std::vector<CurvedFragment>> fragments =
      split_face_curved(face, id, none, 0.5);
  REQUIRE(fragments.has_value());
  // Net area 4 - 1.95^2 = 0.1975 is below tolerance^2 = 0.25.
  CHECK(fragments.value().empty());
}

TEST_CASE("trim_face: an interior sub-tolerance step keeps the cut connected") {
  const Brep box = make_box(Point3d{0, 0, 0}, Point3d{2, 2, 2}, kTol).value();
  const FaceId top = plane_face(box, Vec3d{0, 0, 1}, 2.0);
  REQUIRE_FALSE(top.is_null());

  // The middle step is 1e-9 long, far below the 1e-6 tolerance. Dropping that
  // segment would leave two disconnected pieces and the cut would not split.
  const std::vector<FaceCut> cuts{uv_cut(
      box, top,
      {Point3d{0, 0.5, 2}, Point3d{0.5, 0.5, 2}, Point3d{0.5 + 1e-9, 0.5, 2},
       Point3d{2, 0.5, 2}})};
  const Result<std::vector<CurvedFragment>> fragments =
      split_face_curved(box, top, cuts, 1e-6);
  REQUIRE(fragments.has_value());
  REQUIRE(fragments.value().size() == 2);
  const double first = fragment_area(fragments.value()[0]);
  const double second = fragment_area(fragments.value()[1]);
  CHECK(first + second == doctest::Approx(4.0).epsilon(1e-12));
  CHECK(std::min(first, second) == doctest::Approx(1.0).epsilon(1e-12));
  CHECK(std::max(first, second) == doctest::Approx(3.0).epsilon(1e-12));
}

TEST_CASE("trim_face: a sub-tolerance first step of an open cut is collapsed") {
  const Brep box = make_box(Point3d{0, 0, 0}, Point3d{2, 2, 2}, kTol).value();
  const FaceId top = plane_face(box, Vec3d{0, 0, 1}, 2.0);
  REQUIRE_FALSE(top.is_null());

  // The first step is 1.4e-9 long, below the 1e-6 tolerance. Collapsing it
  // leaves the diagonal from the corner to (2,1), which splits the face into
  // areas 1 and 3; keeping it would shift the split by ~5e-10 in area.
  const std::vector<FaceCut> cuts{uv_cut(
      box, top,
      {Point3d{0, 0, 2}, Point3d{1e-9, 1e-9, 2}, Point3d{2, 1, 2}})};
  const Result<std::vector<CurvedFragment>> fragments =
      split_face_curved(box, top, cuts, 1e-6);
  REQUIRE(fragments.has_value());
  REQUIRE(fragments.value().size() == 2);
  const double first = fragment_area(fragments.value()[0]);
  const double second = fragment_area(fragments.value()[1]);
  CHECK(first + second == doctest::Approx(4.0).epsilon(1e-12));
  CHECK(std::min(first, second) == doctest::Approx(1.0).epsilon(1e-12));
  CHECK(std::max(first, second) == doctest::Approx(3.0).epsilon(1e-12));
}

TEST_CASE("trim_face: a fragment's holes are ordered by descending area") {
  const std::array<Point3d, 4> outer{Point3d{0, 0, 0}, Point3d{4, 0, 0},
                                     Point3d{4, 4, 0}, Point3d{0, 4, 0}};
  // Small hole (area 0.25) at (3,3) and large hole (area 1) at (0.25,0.25).
  // The arrangement discovers the small hole's cycle first, so the area sort is
  // what puts the large hole first in the emitted fragment.
  const std::array<Point3d, 4> small{Point3d{3, 3, 0}, Point3d{3, 3.5, 0},
                                     Point3d{3.5, 3.5, 0}, Point3d{3.5, 3, 0}};
  const std::array<Point3d, 4> large{Point3d{0.25, 0.25, 0},
                                     Point3d{0.25, 1.25, 0},
                                     Point3d{1.25, 1.25, 0},
                                     Point3d{1.25, 0.25, 0}};
  const std::vector<std::vector<Point3d>> holes{
      std::vector<Point3d>(small.begin(), small.end()),
      std::vector<Point3d>(large.begin(), large.end())};
  const Brep face = make_planar_face(outer, holes, Vec3d{0, 0, 1}, 0.0);
  const FaceId id = face.face_ids().front();

  const std::vector<FaceCut> none;
  const Result<std::vector<CurvedFragment>> fragments =
      split_face_curved(face, id, none, kTol);
  REQUIRE(fragments.has_value());
  REQUIRE(fragments.value().size() == 1);
  const CurvedFragment& fragment = fragments.value()[0];
  REQUIRE(fragment.loops_uv.size() == 3);
  CHECK(ring_area(fragment.loops_uv[1]) ==
        doctest::Approx(1.0).epsilon(1e-12));
  CHECK(ring_area(fragment.loops_uv[2]) ==
        doctest::Approx(0.25).epsilon(1e-12));
  CHECK(fragment_area(fragment) == doctest::Approx(14.75).epsilon(1e-12));
}

TEST_CASE("trim_face: invalid input is rejected with a diagnostic") {
  const Brep box = make_box(Point3d{0, 0, 0}, Point3d{2, 2, 2}, kTol).value();
  const FaceId top = plane_face(box, Vec3d{0, 0, 1}, 2.0);
  REQUIRE_FALSE(top.is_null());

  const std::vector<FaceCut> none;
  const auto check_error = [](const Result<std::vector<CurvedFragment>>& result,
                              Errc code) {
    REQUIRE_FALSE(result.has_value());
    CHECK(result.error().code == code);
    CHECK(std::string_view(result.error().message)
              .starts_with("split_face_curved:"));
  };

  check_error(split_face_curved(box, FaceId{}, none, kTol),
              Errc::invalid_topology);
  check_error(split_face_curved(box, top, none,
                                std::numeric_limits<double>::quiet_NaN()),
              Errc::non_finite);
  check_error(split_face_curved(box, top, none, -1e-6),
              Errc::tolerance_violation);

  FaceCut one_point;
  one_point.uv = {Point2d{0.5, 0.5}};
  const std::vector<FaceCut> short_cut{one_point};
  check_error(split_face_curved(box, top, short_cut, kTol),
              Errc::invalid_topology);

  FaceCut non_finite;
  non_finite.uv = {Point2d{0.5, 0.5},
                   Point2d{std::numeric_limits<double>::infinity(), 0.5}};
  const std::vector<FaceCut> bad_point{non_finite};
  check_error(split_face_curved(box, top, bad_point, kTol), Errc::non_finite);

  FaceCut collapsed;
  collapsed.uv = {Point2d{0.5, 0.5}, Point2d{0.5, 0.5}};
  const std::vector<FaceCut> collapsed_cut{collapsed};
  check_error(split_face_curved(box, top, collapsed_cut, kTol),
              Errc::degenerate);
}

TEST_CASE("trim_face: the quarter cylinder fixture is a clean solid") {
  const double radius = 1.0;
  const double height = 2.0;
  const Brep solid = test::make_quarter_cylinder(radius, height, kTol);
  const ValidationReport report = validate_topology(solid);
  CHECK_MESSAGE(report.ok(), "quarter cylinder topology is invalid");
  CHECK(solid.face_ids().size() == 5);

  // Self-check: every vertex is either on the axis or exactly on the cylinder
  // radius, and on one of the two cap planes.
  for (const VertexId vertex : solid.vertex_ids()) {
    const Point3d& point = solid.vertex(vertex).point;
    const double radial = std::hypot(point[0], point[1]);
    CHECK((radial == 0.0 || std::fabs(radial - radius) <= 1e-12 * radius));
    CHECK((std::fabs(point[2]) <= 1e-12 * height ||
           std::fabs(point[2] - height) <= 1e-12 * height));
  }

  const FaceId lateral = test::quarter_cylinder_lateral(solid);
  REQUIRE_FALSE(lateral.is_null());

  const std::vector<FaceCut> none;
  const Result<std::vector<CurvedFragment>> fragments =
      split_face_curved(solid, lateral, none, 1e-9);
  REQUIRE(fragments.has_value());
  REQUIRE(fragments.value().size() == 1);
  CHECK(fragments.value()[0].loops_uv.size() == 1);
  CHECK(fragment_area(fragments.value()[0]) ==
        doctest::Approx(0.5 * kPi * height).epsilon(1e-12));

  // geom::tessellate() only stitches full-circle Cylinder3 caps today, so the
  // partial lateral face is not tessellable; check the chordal polyhedron
  // volume of the fixture instead. Sixteen chords undershoot the true volume
  // by 1.6e-3, far more than 1e-6 (the brief's tolerance would need thousands
  // of chords).
  const double chord_area =
      8.0 * radius * radius * std::sin(kPi / 32.0);
  const double volume = chord_area * height;
  CHECK(volume == doctest::Approx(0.25 * kPi * radius * radius * height)
                       .epsilon(2e-3));
  const Result<TriMesh> mesh = tessellate(solid, TessellationOptions{});
  CHECK_FALSE(mesh.has_value());
}

namespace {

double point_distance(const Point3d& lhs, const Point3d& rhs) {
  return std::hypot(std::hypot(lhs[0] - rhs[0], lhs[1] - rhs[1]),
                    lhs[2] - rhs[2]);
}

std::size_t engaged_span_count(const CurvedFragment& fragment) {
  std::size_t count = 0;
  for (const std::vector<std::optional<CurvedArcSpan>>& loop :
       fragment.arc_spans) {
    for (const std::optional<CurvedArcSpan>& span : loop) {
      if (span) ++count;
    }
  }
  return count;
}

// Checks the alignment contract of every reported span and that each span
// lies on its circle between the two ring vertices it describes. The midpoint
// check pins *which* half of the circle the span covers: reporting the
// complementary arc (same endpoints, opposite direction) must fail.
void check_arc_spans(const Brep& brep, const FaceId face,
                     const CurvedFragment& fragment) {
  REQUIRE(fragment.arc_spans.size() == fragment.loops.size());
  REQUIRE(fragment.arc_spans.size() == fragment.loops_uv.size());
  const Result<detail::FaceUvProjection> projection =
      detail::face_uv_projection(brep, brep.face(face), 1e-12);
  REQUIRE(projection.has_value());
  for (std::size_t l = 0; l < fragment.loops.size(); ++l) {
    const std::size_t n = fragment.loops[l].size();
    REQUIRE(fragment.arc_spans[l].size() == n);
    REQUIRE(fragment.loops_uv[l].size() == n);
    for (std::size_t i = 0; i < n; ++i) {
      const std::optional<CurvedArcSpan>& span = fragment.arc_spans[l][i];
      if (!span) continue;
      const double radius = span->circle.radius;
      CHECK(span->t0 >= 0.0);
      CHECK(span->t1 > span->t0);
      CHECK(span->t1 <= span->t0 + 2.0 * kPi + 1e-12);

      const Point3d& a = fragment.loops[l][i];
      const Point3d& b = fragment.loops[l][(i + 1) % n];
      const Point2d& uv_a = fragment.loops_uv[l][i];
      const Point2d& uv_b = fragment.loops_uv[l][(i + 1) % n];
      const Point3d start = span->circle.eval(span->t0);
      const Point3d end = span->circle.eval(span->t1);
      const double eps = 1e-12 * radius;
      // The ring may follow the arc in either direction; the stored span is the
      // geometric interval between its two vertices.
      const bool forward = point_distance(start, a) <= eps &&
                           point_distance(end, b) <= eps;
      const bool backward = point_distance(start, b) <= eps &&
                            point_distance(end, a) <= eps;
      CHECK((forward || backward));

      const double middle_t = 0.5 * (span->t0 + span->t1);
      const Point3d middle = span->circle.eval(middle_t);
      CHECK(std::fabs(point_distance(middle, span->circle.center) - radius) <=
            eps);

      // The parameter midpoint of the ring segment must be the parameter
      // midpoint of the reported span.
      const Point2d middle_uv{0.5 * (uv_a[0] + uv_b[0]),
                              0.5 * (uv_a[1] + uv_b[1])};
      const Point3d middle_point =
          detail::face_uv_eval(projection.value(), middle_uv);
      const Vec3d radial = middle_point - span->circle.center;
      const double azimuth =
          std::atan2(dot(radial, span->circle.y_axis),
                     dot(radial, span->circle.x_axis));
      double offset = azimuth - middle_t;
      offset -= 2.0 * kPi * std::round(offset / (2.0 * kPi));
      CHECK(std::fabs(offset) <= 1e-9);

      // The width must be the traversed arc, not its complement.
      const Vec3d start_radial = start - span->circle.center;
      const double half_width =
          std::atan2(norm2(cross(start_radial, radial)), dot(start_radial, radial));
      CHECK(std::fabs(0.5 * (span->t1 - span->t0) - half_width) <= 1e-9);
    }
  }
}

// Lengths of the maximal runs of consecutive engaged segments in loop l.
std::vector<double> span_run_lengths(const CurvedFragment& fragment,
                                     std::size_t l) {
  std::vector<double> lengths;
  const std::vector<std::optional<CurvedArcSpan>>& spans = fragment.arc_spans[l];
  const std::size_t n = spans.size();
  if (n == 0) return lengths;
  std::size_t first = n;
  for (std::size_t i = 0; i < n; ++i) {
    if (spans[i] && !spans[(i + n - 1) % n]) {
      first = i;
      break;
    }
  }
  if (first == n) {
    double total = 0.0;
    for (const std::optional<CurvedArcSpan>& span : spans) {
      total += span->t1 - span->t0;
    }
    lengths.push_back(total);
    return lengths;
  }
  double total = 0.0;
  for (std::size_t step = 0; step < n; ++step) {
    const std::size_t index = (first + step) % n;
    if (spans[index]) {
      total += spans[index]->t1 - spans[index]->t0;
    } else if (total > 0.0) {
      lengths.push_back(total);
      total = 0.0;
    }
  }
  if (total > 0.0) lengths.push_back(total);
  return lengths;
}

// Total reported arc length on one circle across every fragment.
double span_length_on_circle(const std::vector<CurvedFragment>& fragments,
                             const Point3d& center, double radius) {
  double total = 0.0;
  for (const CurvedFragment& fragment : fragments) {
    for (const std::vector<std::optional<CurvedArcSpan>>& loop :
         fragment.arc_spans) {
      for (const std::optional<CurvedArcSpan>& span : loop) {
        if (!span) continue;
        const double scale = std::max(1.0, radius);
        if (std::fabs(span->circle.radius - radius) > 1e-12 * scale) continue;
        if (point_distance(span->circle.center, center) > 1e-12 * scale) {
          continue;
        }
        total += span->t1 - span->t0;
      }
    }
  }
  return total;
}

// A single half-cylinder lateral patch (u in [0, pi]) bounded by two exact
// semicircle rims. Each rim is one CircleCurve3 edge carrying a 17-point pcurve,
// so the patch has four vertices and each rim is sampled by 16 arc segments.
Brep make_half_cylinder_patch(double radius, double height, double tol) {
  SolidEditor editor;
  const std::uint32_t surface = editor.add_surface(
      Cylinder3{Point3d{0, 0, 0}, Vec3d{0, 0, 1}, Vec3d{1, 0, 0}, radius});
  const CircleCurve3 bottom_circle{Point3d{0, 0, 0}, Vec3d{1, 0, 0},
                                   Vec3d{0, 1, 0}, radius};
  const CircleCurve3 top_circle{Point3d{0, 0, height}, Vec3d{1, 0, 0},
                                Vec3d{0, 1, 0}, radius};
  const Point3d b0 = bottom_circle.eval(0.0);
  const Point3d b1 = bottom_circle.eval(kPi);
  const Point3d t0 = top_circle.eval(0.0);
  const Point3d t1 = top_circle.eval(kPi);
  const VertexId vb0 = editor.add_vertex(b0, 0.0).value();
  const VertexId vb1 = editor.add_vertex(b1, 0.0).value();
  const VertexId vt0 = editor.add_vertex(t0, 0.0).value();
  const VertexId vt1 = editor.add_vertex(t1, 0.0).value();
  const EdgeId bottom_edge =
      editor
          .add_edge(vb0, vb1,
                    editor.add_curve(CircleCurve3{bottom_circle}), 0.0, kPi, 0.0)
          .value();
  const EdgeId top_edge =
      editor
          .add_edge(vt0, vt1, editor.add_curve(CircleCurve3{top_circle}), 0.0,
                    kPi, 0.0)
          .value();
  const EdgeId rule_zero =
      editor
          .add_edge(vb0, vt0,
                    editor.add_curve(LineCurve3{b0, Vec3d{0.0, 0.0, height}}),
                    0.0, 1.0, 0.0)
          .value();
  const EdgeId rule_pi =
      editor
          .add_edge(vb1, vt1,
                    editor.add_curve(LineCurve3{b1, Vec3d{0.0, 0.0, height}}),
                    0.0, 1.0, 0.0)
          .value();

  // A two-point pcurve would leave each half-turn segment's direction ambiguous,
  // so the rims carry 17 samples and 16 arc segments each.
  std::vector<Point2d> bottom_uv;
  std::vector<Point2d> top_uv;
  for (std::size_t i = 0; i <= 16u; ++i) {
    const double angle = kPi * static_cast<double>(i) / 16.0;
    bottom_uv.push_back(Point2d{angle, 0.0});
    top_uv.push_back(Point2d{angle, height});
  }
  const std::uint32_t bottom_pcurve = editor.add_pcurve(Pcurve2d{bottom_uv});
  const std::uint32_t top_pcurve = editor.add_pcurve(Pcurve2d{top_uv});
  const std::uint32_t rule_zero_pcurve =
      editor.add_pcurve(Pcurve2d{{Point2d{0.0, 0.0}, Point2d{0.0, height}}});
  const std::uint32_t rule_pi_pcurve =
      editor.add_pcurve(Pcurve2d{{Point2d{kPi, 0.0}, Point2d{kPi, height}}});

  const FaceId lateral = editor.add_face(surface, false, tol).value();
  std::vector<CoedgeSpec> specs;
  const auto push = [&specs](EdgeId edge, bool reversed, std::uint32_t pcurve) {
    CoedgeSpec spec;
    spec.edge = edge;
    spec.reversed = reversed;
    spec.pcurve = pcurve;
    spec.tol = 0.0;
    specs.push_back(spec);
  };
  // u = 0 -> u = pi along the bottom rim, up the pi rule, back along the top
  // rim, and down the zero rule.
  push(bottom_edge, false, bottom_pcurve);
  push(rule_pi, false, rule_pi_pcurve);
  push(top_edge, true, top_pcurve);
  push(rule_zero, true, rule_zero_pcurve);
  editor.add_loop(lateral, specs, true).value();
  editor.add_shell(std::array<FaceId, 1>{lateral}, false).value();
  return editor.commit().value();
}

}  // namespace

TEST_CASE("trim_face: an axial cut splits the cylinder wall into two arc runs") {
  const double radius = 2.0;
  const double height = 3.0;
  const Brep cylinder =
      make_cylinder(Point3d{0, 0, 0}, Vec3d{0, 0, 1}, radius, height, kTol)
          .value();
  const FaceId lateral = face_with_surface<Cylinder3>(cylinder);
  REQUIRE_FALSE(lateral.is_null());

  FaceCut cut;
  cut.uv = {Point2d{kPi, 0.0}, Point2d{kPi, height}};
  const std::vector<FaceCut> cuts{cut};
  const Result<std::vector<CurvedFragment>> result =
      split_face_curved(cylinder, lateral, cuts, 1e-9);
  REQUIRE(result.has_value());
  REQUIRE(result.value().size() == 2);

  for (const CurvedFragment& fragment : result.value()) {
    check_arc_spans(cylinder, lateral, fragment);
    CHECK(fragment_area(fragment) ==
          doctest::Approx(kPi * height).epsilon(1e-12));
    // Each half follows half of the bottom rim and half of the top rim.
    const std::vector<double> lengths = span_run_lengths(fragment, 0);
    REQUIRE(lengths.size() == 2);
    CHECK(lengths[0] == doctest::Approx(kPi).epsilon(1e-12));
    CHECK(lengths[1] == doctest::Approx(kPi).epsilon(1e-12));
  }

  // The two fragments together cover each rim exactly once.
  CHECK(span_length_on_circle(result.value(), Point3d{0, 0, 0}, radius) ==
        doctest::Approx(2.0 * kPi).epsilon(1e-12));
  CHECK(span_length_on_circle(result.value(), Point3d{0, 0, height}, radius) ==
        doctest::Approx(2.0 * kPi).epsilon(1e-12));
}

TEST_CASE("trim_face: a diameter cut splits both cylinder caps along their arcs") {
  const double radius = 2.0;
  const double height = 3.0;
  const Brep cylinder =
      make_cylinder(Point3d{0, 0, 0}, Vec3d{0, 0, 1}, radius, height, kTol)
          .value();
  const FaceId bottom = plane_face(cylinder, Vec3d{0, 0, -1}, 0.0);
  const FaceId top = plane_face(cylinder, Vec3d{0, 0, 1}, height);
  REQUIRE_FALSE(bottom.is_null());
  REQUIRE_FALSE(top.is_null());

  const auto check_cap = [&](const FaceId face, const Point3d& center) {
    const CircleCurve3 circle{center, Vec3d{1, 0, 0}, Vec3d{0, 1, 0}, radius};
    const std::vector<FaceCut> cuts{uv_cut(
        cylinder, face, {circle.eval(kPi / 2.0), circle.eval(3.0 * kPi / 2.0)})};
    const Result<std::vector<CurvedFragment>> result =
        split_face_curved(cylinder, face, cuts, 1e-9);
    REQUIRE(result.has_value());
    REQUIRE(result.value().size() == 2);
    for (const CurvedFragment& fragment : result.value()) {
      check_arc_spans(cylinder, face, fragment);
      const std::vector<double> lengths = span_run_lengths(fragment, 0);
      REQUIRE(lengths.size() == 1);
      CHECK(lengths[0] == doctest::Approx(kPi).epsilon(1e-12));
    }
    CHECK(span_length_on_circle(result.value(), center, radius) ==
          doctest::Approx(2.0 * kPi).epsilon(1e-12));
  };

  check_cap(bottom, Point3d{0, 0, 0});
  check_cap(top, Point3d{0, 0, height});
}

TEST_CASE("trim_face: a mid-height cut leaves one full-circle arc run per half") {
  const double radius = 2.0;
  const double height = 3.0;
  const Brep cylinder =
      make_cylinder(Point3d{0, 0, 0}, Vec3d{0, 0, 1}, radius, height, kTol)
          .value();
  const FaceId lateral = face_with_surface<Cylinder3>(cylinder);
  REQUIRE_FALSE(lateral.is_null());

  FaceCut cut;
  cut.uv = {Point2d{0.0, height / 2.0}, Point2d{2.0 * kPi, height / 2.0}};
  const std::vector<FaceCut> cuts{cut};
  const Result<std::vector<CurvedFragment>> result =
      split_face_curved(cylinder, lateral, cuts, 1e-9);
  REQUIRE(result.has_value());
  REQUIRE(result.value().size() == 2);

  double bottom = 0.0;
  double top = 0.0;
  for (const CurvedFragment& fragment : result.value()) {
    check_arc_spans(cylinder, lateral, fragment);
    const std::vector<double> lengths = span_run_lengths(fragment, 0);
    REQUIRE(lengths.size() == 1);
    CHECK(lengths[0] == doctest::Approx(2.0 * kPi).epsilon(1e-12));
    if (std::fabs(fragment.loops[0][0][2]) < 1e-12) {
      bottom += lengths[0];
    } else {
      top += lengths[0];
    }
  }
  CHECK(bottom == doctest::Approx(2.0 * kPi).epsilon(1e-12));
  CHECK(top == doctest::Approx(2.0 * kPi).epsilon(1e-12));
}

TEST_CASE("trim_face: a polygonal face reports no arc spans") {
  const Brep box = make_box(Point3d{0, 0, 0}, Point3d{2, 2, 2}, kTol).value();
  const FaceId top = plane_face(box, Vec3d{0, 0, 1}, 2.0);
  REQUIRE_FALSE(top.is_null());

  const std::vector<FaceCut> none;
  const Result<std::vector<CurvedFragment>> whole =
      split_face_curved(box, top, none, kTol);
  REQUIRE(whole.has_value());
  REQUIRE(whole.value().size() == 1);
  check_arc_spans(box, top, whole.value()[0]);
  CHECK(engaged_span_count(whole.value()[0]) == 0);

  const std::vector<FaceCut> cuts{
      uv_cut(box, top, {Point3d{0, 0.5, 2}, Point3d{2, 1.5, 2}})};
  const Result<std::vector<CurvedFragment>> split =
      split_face_curved(box, top, cuts, kTol);
  REQUIRE(split.has_value());
  REQUIRE(split.value().size() == 2);
  for (const CurvedFragment& fragment : split.value()) {
    check_arc_spans(box, top, fragment);
    CHECK(engaged_span_count(fragment) == 0);
  }
}

TEST_CASE("trim_face: a cut ending on an arc splits it at the exact junction") {
  const double radius = 2.0;
  const double height = 3.0;
  const Brep cylinder =
      make_cylinder(Point3d{0, 0, 0}, Vec3d{0, 0, 1}, radius, height, kTol)
          .value();
  const FaceId lateral = face_with_surface<Cylinder3>(cylinder);
  REQUIRE_FALSE(lateral.is_null());

  FaceCut cut;
  cut.uv = {Point2d{0.0, height / 2.0}, Point2d{0.3, 0.0}};
  const std::vector<FaceCut> cuts{cut};
  const Result<std::vector<CurvedFragment>> result =
      split_face_curved(cylinder, lateral, cuts, 1e-9);
  REQUIRE(result.has_value());
  REQUIRE(result.value().size() == 2);

  std::vector<double> lengths;
  for (const CurvedFragment& fragment : result.value()) {
    check_arc_spans(cylinder, lateral, fragment);
    const std::vector<double> run = span_run_lengths(fragment, 0);
    lengths.insert(lengths.end(), run.begin(), run.end());
  }
  std::sort(lengths.begin(), lengths.end());
  REQUIRE(lengths.size() == 3);
  CHECK(lengths[0] == doctest::Approx(0.3).epsilon(1e-12));
  CHECK(lengths[1] == doctest::Approx(2.0 * kPi - 0.3).epsilon(1e-12));
  CHECK(lengths[2] == doctest::Approx(2.0 * kPi).epsilon(1e-12));
  CHECK(span_length_on_circle(result.value(), Point3d{0, 0, 0}, radius) ==
        doctest::Approx(2.0 * kPi).epsilon(1e-12));
  CHECK(span_length_on_circle(result.value(), Point3d{0, 0, height}, radius) ==
        doctest::Approx(2.0 * kPi).epsilon(1e-12));
}

TEST_CASE("trim_face: a chord cuts a cap into a major and a minor arc run") {
  const double radius = 2.0;
  const double height = 3.0;
  const Brep cylinder =
      make_cylinder(Point3d{0, 0, 0}, Vec3d{0, 0, 1}, radius, height, kTol)
          .value();
  const FaceId bottom = plane_face(cylinder, Vec3d{0, 0, -1}, 0.0);
  REQUIRE_FALSE(bottom.is_null());

  const CircleCurve3 circle{Point3d{0, 0, 0}, Vec3d{1, 0, 0}, Vec3d{0, 1, 0},
                            radius};
  const std::vector<FaceCut> cuts{uv_cut(
      cylinder, bottom, {circle.eval(0.4), circle.eval(1.1)})};
  const Result<std::vector<CurvedFragment>> result =
      split_face_curved(cylinder, bottom, cuts, 1e-9);
  REQUIRE(result.has_value());
  REQUIRE(result.value().size() == 2);

  std::vector<double> lengths;
  for (const CurvedFragment& fragment : result.value()) {
    check_arc_spans(cylinder, bottom, fragment);
    const std::vector<double> run = span_run_lengths(fragment, 0);
    REQUIRE(run.size() == 1);
    lengths.push_back(run[0]);
  }
  std::sort(lengths.begin(), lengths.end());
  REQUIRE(lengths.size() == 2);
  CHECK(lengths[0] == doctest::Approx(0.7).epsilon(1e-12));
  CHECK(lengths[1] == doctest::Approx(2.0 * kPi - 0.7).epsilon(1e-12));
  CHECK(span_length_on_circle(result.value(), circle.center, radius) ==
        doctest::Approx(2.0 * kPi).epsilon(1e-12));
}

TEST_CASE("trim_face: the chordal quarter cylinder reports no arc spans") {
  const double radius = 1.0;
  const double height = 2.0;
  const Brep solid = test::make_quarter_cylinder(radius, height, kTol);
  const FaceId lateral = test::quarter_cylinder_lateral(solid);
  REQUIRE_FALSE(lateral.is_null());

  const std::vector<FaceCut> none;
  const Result<std::vector<CurvedFragment>> whole =
      split_face_curved(solid, lateral, none, 1e-9);
  REQUIRE(whole.has_value());
  REQUIRE(whole.value().size() == 1);
  check_arc_spans(solid, lateral, whole.value()[0]);
  CHECK(engaged_span_count(whole.value()[0]) == 0);

  FaceCut cut;
  cut.uv = {Point2d{kPi / 4.0, 0.0}, Point2d{kPi / 4.0, height}};
  const std::vector<FaceCut> cuts{cut};
  const Result<std::vector<CurvedFragment>> split =
      split_face_curved(solid, lateral, cuts, 1e-9);
  REQUIRE(split.has_value());
  REQUIRE(split.value().size() == 2);
  for (const CurvedFragment& fragment : split.value()) {
    check_arc_spans(solid, lateral, fragment);
    CHECK(engaged_span_count(fragment) == 0);
  }
}

TEST_CASE("trim_face: a cut tangent to a cap circle does not split it") {
  const double radius = 2.0;
  const double height = 3.0;
  const Brep cylinder =
      make_cylinder(Point3d{0, 0, 0}, Vec3d{0, 0, 1}, radius, height, kTol)
          .value();
  const FaceId bottom = plane_face(cylinder, Vec3d{0, 0, -1}, 0.0);
  REQUIRE_FALSE(bottom.is_null());

  const std::vector<FaceCut> none;
  const Result<std::vector<CurvedFragment>> whole =
      split_face_curved(cylinder, bottom, none, 1e-9);
  REQUIRE(whole.has_value());
  REQUIRE(whole.value().size() == 1);
  const double whole_area = fragment_area(whole.value()[0]);
  const std::size_t whole_points = whole.value()[0].loops_uv[0].size();

  // A closed square whose left edge touches the cap circle at one point.
  const Point2d center_uv = uv_cut(cylinder, bottom, {Point3d{0, 0, 0}}).uv[0];
  FaceCut square;
  square.uv = {Point2d{center_uv[0] - radius, center_uv[1] - radius - 1.0},
               Point2d{center_uv[0] + radius + 1.0, center_uv[1] - radius - 1.0},
               Point2d{center_uv[0] + radius + 1.0, center_uv[1] + radius + 1.0},
               Point2d{center_uv[0] - radius, center_uv[1] + radius + 1.0},
               Point2d{center_uv[0] - radius, center_uv[1] - radius - 1.0}};
  const std::vector<FaceCut> cuts{square};
  const Result<std::vector<CurvedFragment>> result =
      split_face_curved(cylinder, bottom, cuts, 1e-9);
  REQUIRE(result.has_value());
  REQUIRE(result.value().size() == 1);
  check_arc_spans(cylinder, bottom, result.value()[0]);
  CHECK(result.value()[0].loops_uv[0].size() == whole_points);
  CHECK(fragment_area(result.value()[0]) ==
        doctest::Approx(whole_area).epsilon(1e-12));
  const std::vector<double> lengths = span_run_lengths(result.value()[0], 0);
  REQUIRE(lengths.size() == 1);
  CHECK(lengths[0] == doctest::Approx(2.0 * kPi).epsilon(1e-12));
  CHECK(engaged_span_count(result.value()[0]) == whole_points);
}

TEST_CASE("trim_face: a half-cylinder patch reports its exact half-turn spans") {
  const double radius = 2.0;
  const double height = 3.0;
  const Brep patch = make_half_cylinder_patch(radius, height, kTol);
  const FaceId lateral = face_with_surface<Cylinder3>(patch);
  REQUIRE_FALSE(lateral.is_null());

  // Uncut, each rim is a half turn sampled by 16 arc segments, so 32 spans are
  // engaged and the two rims cover pi each. A per-point provenance scheme would
  // instead report one ref per vertex and no way to tell the rims apart from the
  // two straight rules.
  const std::vector<FaceCut> none;
  const Result<std::vector<CurvedFragment>> whole =
      split_face_curved(patch, lateral, none, 1e-9);
  REQUIRE(whole.has_value());
  REQUIRE(whole.value().size() == 1);
  const CurvedFragment& fragment = whole.value()[0];
  check_arc_spans(patch, lateral, fragment);
  REQUIRE(fragment.loops_uv[0].size() == 34);
  REQUIRE(engaged_span_count(fragment) == 32);
  double whole_bottom = 0.0;
  double whole_top = 0.0;
  for (const std::optional<CurvedArcSpan>& span : fragment.arc_spans[0]) {
    if (!span) continue;
    CHECK(span->t1 - span->t0 == doctest::Approx(kPi / 16.0).epsilon(1e-12));
    if (span->circle.center[2] == 0.0) {
      whole_bottom += span->t1 - span->t0;
    } else {
      whole_top += span->t1 - span->t0;
    }
  }
  CHECK(whole_bottom == doctest::Approx(kPi).epsilon(1e-12));
  CHECK(whole_top == doctest::Approx(kPi).epsilon(1e-12));

  // Cutting at u = pi/2 splits each half turn at its midpoint segment boundary.
  FaceCut cut;
  cut.uv = {Point2d{0.5 * kPi, 0.0}, Point2d{0.5 * kPi, height}};
  const std::vector<FaceCut> cuts{cut};
  const Result<std::vector<CurvedFragment>> split =
      split_face_curved(patch, lateral, cuts, 1e-9);
  REQUIRE(split.has_value());
  REQUIRE(split.value().size() == 2);
  std::size_t engaged = 0;
  double bottom = 0.0;
  double top = 0.0;
  for (const CurvedFragment& piece : split.value()) {
    check_arc_spans(patch, lateral, piece);
    engaged += engaged_span_count(piece);
    for (const std::optional<CurvedArcSpan>& span : piece.arc_spans[0]) {
      if (!span) continue;
      CHECK(span->t1 - span->t0 ==
            doctest::Approx(kPi / 16.0).epsilon(1e-12));
      if (span->circle.center[2] == 0.0) {
        bottom += span->t1 - span->t0;
      } else {
        top += span->t1 - span->t0;
      }
    }
  }
  // Sixteen half-turn spans per piece, and the two pieces cover each semicircle
  // rim exactly once (a semicircle rim has parameter length pi, so the two rims
  // together contribute 2*pi).
  CHECK(engaged == 32);
  CHECK(bottom == doctest::Approx(kPi).epsilon(1e-12));
  CHECK(top == doctest::Approx(kPi).epsilon(1e-12));
}

TEST_CASE("trim_face: a tangency between boundary vertices keeps the cap whole") {
  const double radius = 2.0;
  const double height = 3.0;
  const Brep cylinder =
      make_cylinder(Point3d{0, 0, 0}, Vec3d{0, 0, 1}, radius, height, kTol)
          .value();
  const FaceId bottom = plane_face(cylinder, Vec3d{0, 0, -1}, 0.0);
  REQUIRE_FALSE(bottom.is_null());

  const std::vector<FaceCut> none;
  const Result<std::vector<CurvedFragment>> whole =
      split_face_curved(cylinder, bottom, none, 1e-9);
  REQUIRE(whole.has_value());
  REQUIRE(whole.value().size() == 1);
  const double whole_area = fragment_area(whole.value()[0]);
  const std::size_t whole_points = whole.value()[0].loops_uv[0].size();

  // A square whose inner edge is tangent to the cap circle at uv angle 0.5,
  // between boundary vertices, sinking 1e-10 below the radius. That is a real
  // (very shallow) chord of the analytic circle but a tangency within the
  // 1e-9 tolerance, so it must be ignored.
  const Point2d center_uv = uv_cut(cylinder, bottom, {Point3d{0, 0, 0}}).uv[0];
  const double angle = 0.5;
  const Point2d radial{std::cos(angle), std::sin(angle)};
  const Point2d tangent{-std::sin(angle), std::cos(angle)};
  const double penetration = 1e-10;
  const auto corner = [&](double outward, double side) {
    return Point2d{center_uv[0] + outward * radial[0] + side * tangent[0],
                   center_uv[1] + outward * radial[1] + side * tangent[1]};
  };
  FaceCut square;
  square.uv = {corner(radius - penetration, -1.0),
               corner(radius - penetration, 1.0),
               corner(radius - penetration + 2.0, 1.0),
               corner(radius - penetration + 2.0, -1.0),
               corner(radius - penetration, -1.0)};
  const std::vector<FaceCut> cuts{square};
  const Result<std::vector<CurvedFragment>> result =
      split_face_curved(cylinder, bottom, cuts, 1e-9);
  REQUIRE(result.has_value());
  REQUIRE(result.value().size() == 1);
  check_arc_spans(cylinder, bottom, result.value()[0]);
  CHECK(result.value()[0].loops_uv[0].size() == whole_points);
  CHECK(engaged_span_count(result.value()[0]) == whole_points);
  CHECK(fragment_area(result.value()[0]) ==
        doctest::Approx(whole_area).epsilon(1e-12));
}

TEST_CASE("trim_face: circle hits accept a cut that penetrates past the tolerance") {
  const geom::CircleCurve3 circle{Point3d{0, 0, 0}, Vec3d{1, 0, 0}, Vec3d{0, 1, 0}, 2.0};
  // A chord 1e-7 inside the radius is a real penetration at tol 1e-9.
  const auto inner = geom::detail::trim_face_circle_hits(
      circle, Point3d{-1.5, 2.0 - 1e-7, 0.0}, Point3d{1.5, 2.0 - 1e-7, 0.0}, 1e-9);
  REQUIRE(inner.size() == 2);
  for (const geom::detail::TrimFaceCircleHit& hit : inner) {
    CHECK(std::fabs(geom::norm2(hit.point - circle.center) - circle.radius) <= 1e-12);
    CHECK(hit.s > 0.0);
    CHECK(hit.s < 1.0);
    CHECK(hit.point[1] == doctest::Approx(2.0 - 1e-7).epsilon(1e-12));
  }
  // 1e-10 inside is within the tolerance: a graze, not a cut.
  CHECK(geom::detail::trim_face_circle_hits(
            circle, Point3d{-1.5, 2.0 - 1e-10, 0.0}, Point3d{1.5, 2.0 - 1e-10, 0.0},
            1e-9)
            .empty());
  // A chord outside the circle does not hit at all.
  CHECK(geom::detail::trim_face_circle_hits(
            circle, Point3d{-1.5, 2.0 + 1e-7, 0.0}, Point3d{1.5, 2.0 + 1e-7, 0.0},
            1e-9)
            .empty());
}

TEST_CASE("trim_face: two cuts weld at one exact interior corner") {
  constexpr double kSplitTol = 1e-9;
  const Brep box =
      make_box(Point3d{0, 0, 0}, Point3d{10, 10, 10}, kSplitTol).value();
  const FaceId face = box_plane_face(box, 10.0);
  REQUIRE_FALSE(face.is_null());

  // Two cuts end 5e-5 apart with different slopes. Welding onto either raw end
  // would drag the other cut off its line, so the corner is the intersection of
  // the two end segments.
  const double d = 5e-5;
  std::vector<FaceCut> cuts(2);
  cuts[0].uv = {Point2d{2.0, 0.0}, Point2d{5.0, 4.0}};
  cuts[1].uv = {Point2d{8.0, 0.0}, Point2d{5.0 + d, 4.0 + d}};
  const Result<std::vector<CurvedFragment>> result =
      split_face_curved(box, face, cuts, kSplitTol);
  REQUIRE(result.has_value());
  REQUIRE(result.value().size() == 2);

  std::vector<double> areas;
  areas.reserve(2);
  for (const CurvedFragment& fragment : result.value()) {
    areas.push_back(fragment_area(fragment));
  }
  std::sort(areas.begin(), areas.end());
  CHECK(areas[0] == doctest::Approx(12.000175000).epsilon(1e-9));
  CHECK(areas[1] == doctest::Approx(87.999825000).epsilon(1e-9));
  CHECK(areas[0] + areas[1] == doctest::Approx(100.0).epsilon(1e-12));

  // Both fragments carry the same corner, and it is the exact intersection of
  // the two end segments.
  const Point2d expected{5.000043750, 4.000058333};
  const double junction_tol =
      std::max(64.0 * kSplitTol, 1e-5 * std::hypot(10.0, 10.0));
  std::vector<Point2d> corners;
  for (const CurvedFragment& fragment : result.value()) {
    Point2d nearest{};
    double best = std::numeric_limits<double>::max();
    for (const Point2d& p : fragment.loops_uv.front()) {
      const double distance = std::hypot(p[0] - 5.0, p[1] - 4.0);
      if (distance < best) {
        best = distance;
        nearest = p;
      }
    }
    REQUIRE(best <= 1e-3);
    CHECK(std::hypot(nearest[0] - expected[0], nearest[1] - expected[1]) <=
          1e-9);
    // The welded corner stays a uv displacement inside the junction tolerance
    // (a bound on the raw end-segment parameter would allow 4x more here).
    CHECK(best <= 2.0 * junction_tol);
    corners.push_back(nearest);
  }
  CHECK(geom::detail::trim_face_same_xy(corners[0], corners[1]));
}

TEST_CASE("trim_face: a near-parallel digon collapses without a corner") {
  constexpr double kSplitTol = 1e-9;
  const Brep box =
      make_box(Point3d{0, 0, 0}, Point3d{10, 10, 10}, kSplitTol).value();
  const FaceId face = box_plane_face(box, 10.0);
  REQUIRE_FALSE(face.is_null());

  // The two end segments are nearly parallel: their lines meet far outside the
  // face, so the corner recompute must fall back instead of producing a point.
  const double d = 5e-5;
  std::vector<FaceCut> cuts(2);
  cuts[0].uv = {Point2d{3.0, 3.0}, Point2d{7.0, 3.0}};
  cuts[1].uv = {Point2d{7.0 + d, 3.0 + d}, Point2d{3.0 - d, 3.0 - d}};
  const Result<std::vector<CurvedFragment>> result =
      split_face_curved(box, face, cuts, kSplitTol);
  REQUIRE(result.has_value());
  REQUIRE(result.value().size() == 1);
  const double area = fragment_area(result.value().front());
  CHECK(std::isfinite(area));
  CHECK(area == doctest::Approx(100.0).epsilon(1e-9));
}

TEST_CASE("trim_face: cut ends on a boundary edge weld onto one vertex") {
  constexpr double kSplitTol = 1e-9;
  const Brep box =
      make_box(Point3d{0, 0, 0}, Point3d{10, 10, 10}, kSplitTol).value();
  const FaceId face = box_plane_face(box, 10.0);
  REQUIRE_FALSE(face.is_null());

  // A closed triangle network whose two cut ends land 5e-5 apart on the same
  // straight boundary edge: both ends must weld onto the one boundary point.
  std::vector<FaceCut> cuts(3);
  cuts[0].uv = {Point2d{2.0, 2.0}, Point2d{5.0, 0.0}};
  cuts[1].uv = {Point2d{5.00005, 0.0}, Point2d{8.0, 2.0}};
  cuts[2].uv = {Point2d{8.0, 2.0}, Point2d{2.0, 2.0}};
  const Result<std::vector<CurvedFragment>> result =
      split_face_curved(box, face, cuts, kSplitTol);
  REQUIRE(result.has_value());
  REQUIRE(result.value().size() == 2);

  std::vector<double> areas;
  areas.reserve(2);
  for (const CurvedFragment& fragment : result.value()) {
    areas.push_back(fragment_area(fragment));
  }
  std::sort(areas.begin(), areas.end());
  CHECK(areas[0] == doctest::Approx(6.0).epsilon(1e-9));
  CHECK(areas[1] == doctest::Approx(94.0).epsilon(1e-9));

  // The outer ring passes through the single welded point; no vertex is left in
  // the 5e-5 gap between the two raw ends.
  std::vector<Point2d> welded;
  for (const CurvedFragment& fragment : result.value()) {
    const std::vector<Point2d>& ring = fragment.loops_uv.front();
    for (const Point2d& p : ring) {
      const double distance = std::hypot(p[0] - 5.0, p[1]);
      if (distance <= 1e-12) {
        bool seen = false;
        for (const Point2d& q : welded) {
          if (geom::detail::trim_face_same_xy(q, p)) seen = true;
        }
        if (!seen) welded.push_back(p);
      }
      const bool in_gap = distance > 1e-12 && distance <= 1e-6;
      CHECK_FALSE(in_gap);
    }
  }
  REQUIRE(welded.size() == 1);
  CHECK(std::hypot(welded[0][0] - 5.0, welded[0][1]) <= 1e-12);
}

TEST_CASE("trim_face: an endpoint just inside a rim arc snaps onto the circle") {
  constexpr double kSplitTol = 1e-9;
  const Brep cylinder =
      make_cylinder(Point3d{0, 0, 0}, Vec3d{0, 0, 1}, 1.0, 2.0, kSplitTol)
          .value();
  const FaceId lateral = face_with_surface<Cylinder3>(cylinder);
  REQUIRE_FALSE(lateral.is_null());

  // The endpoint stops 5e-5 short of the rim, midway between two rim samples
  // (u = 1.5 rad). It must be projected radially onto the rim circle, not
  // dropped onto the chord, so the cut reaches the arc exactly.
  FaceCut cut;
  cut.uv = {Point2d{0.0, 1.0}, Point2d{1.5, 5e-5}};
  const std::vector<FaceCut> cuts{cut};
  const Result<std::vector<CurvedFragment>> result =
      split_face_curved(cylinder, lateral, cuts, kSplitTol);
  REQUIRE(result.has_value());
  REQUIRE(result.value().size() == 2);

  std::vector<double> areas;
  areas.reserve(2);
  for (const CurvedFragment& fragment : result.value()) {
    areas.push_back(fragment_area(fragment));
  }
  std::sort(areas.begin(), areas.end());
  CHECK(areas[0] == doctest::Approx(0.75).epsilon(1e-9));
  CHECK(areas[1] == doctest::Approx(4.0 * kPi - 0.75).epsilon(1e-9));

  // The welded vertex is the exact circle point at u = 1.5; the raw endpoint
  // 5e-5 inside is gone.
  bool on_rim = false;
  for (const CurvedFragment& fragment : result.value()) {
    for (const Point2d& p : fragment.loops_uv.front()) {
      const double to_rim = std::hypot(p[0] - 1.5, p[1]);
      if (to_rim <= 1e-12) on_rim = true;
      const bool in_gap = to_rim > 1e-12 && to_rim <= 1e-6;
      CHECK_FALSE(in_gap);
    }
  }
  CHECK(on_rim);
}
