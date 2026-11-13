#include <doctest/doctest.h>

#include <geom/brep/boolean_solid.hpp>
#include <geom/brep/boolean_solid_curved.hpp>
#include <geom/brep/face_uv.hpp>
#include <geom/brep/primitives.hpp>
#include <geom/cgeom/polygon2.hpp>
#include <geom/cgeom/trimmed_domain.hpp>
#include <geom/core/result.hpp>

#include <array>
#include <cmath>
#include <numbers>
#include <optional>
#include <vector>

namespace {

constexpr double kPi = std::numbers::pi_v<double>;

geom::Brep cut_cylinder() {
  const geom::Brep cylinder =
      geom::make_cylinder(geom::Point3d{0.0, 0.0, 0.0}, geom::Vec3d{0.0, 0.0, 1.0}, 1.0, 2.0,
                          1e-9)
          .value();
  const geom::Brep box =
      geom::make_box(geom::Point3d{0.5, -0.5, 0.5}, geom::Point3d{2.0, 0.5, 1.5}, 1e-9).value();
  const geom::Result<geom::CurvedBooleanResult> result = geom::boolean_solid_curved(
      cylinder, cylinder.solid_ids().front(), box, box.solid_ids().front(),
      geom::SolidBooleanOp::difference, geom::CurvedBooleanOptions{1e-9});
  REQUIRE(result.has_value());
  return result.value().brep;
}

std::vector<geom::FaceId> cylinder_faces(const geom::Brep& brep) {
  std::vector<geom::FaceId> faces;
  for (const geom::FaceId id : brep.face_ids()) {
    const geom::Face& face = brep.face(id);
    if (face.surface < brep.geometry().surfaces.size() &&
        std::holds_alternative<geom::Cylinder3>(brep.geometry().surfaces[face.surface])) {
      faces.push_back(id);
    }
  }
  return faces;
}

double polygon_area(const std::vector<geom::TrimmedBoundary>& boundaries) {
  double total = 0.0;
  for (const geom::TrimmedBoundary& ring : boundaries) {
    total += geom::signed_area2(ring.points);
  }
  return total;
}

double triangle_area(const geom::TrimmedUvTriangulation& triangulation) {
  double total = 0.0;
  for (const std::array<std::uint32_t, 3>& triangle : triangulation.triangles) {
    const geom::Point2d& a = triangulation.points[triangle[0]];
    const geom::Point2d& b = triangulation.points[triangle[1]];
    const geom::Point2d& c = triangulation.points[triangle[2]];
    total += 0.5 * ((b[0] - a[0]) * (c[1] - a[1]) - (b[1] - a[1]) * (c[0] - a[0]));
  }
  return total;
}

// Planar plate 4 x 4 x 1 with a 1 x 1 through-hole at its centre. Both z faces
// carry an outer loop and a hole loop, so they exercise boundary ordering and
// hole winding that the curved window fixture cannot reach.
geom::Brep plate_with_hole() {
  const geom::Brep plate =
      geom::make_box(geom::Point3d{0.0, 0.0, 0.0}, geom::Point3d{4.0, 4.0, 1.0}, 1e-9).value();
  const geom::Brep cutter =
      geom::make_box(geom::Point3d{1.5, 1.5, -1.0}, geom::Point3d{2.5, 2.5, 2.0}, 1e-9).value();
  const geom::Result<geom::SolidBooleanResult> result = geom::boolean_solid(
      plate, cutter, geom::SolidBooleanOp::difference, geom::SolidBooleanOptions{1e-9});
  REQUIRE(result.has_value());
  return result.value().brep;
}

geom::detail::FaceUvProjection uv_projection(const geom::Brep& brep, geom::FaceId face) {
  const geom::Result<geom::detail::FaceUvProjection> projection =
      geom::detail::face_uv_projection(brep, brep.face(face), 1e-9);
  REQUIRE(projection.has_value());
  return projection.value();
}

// Drives the extraction/orientation path directly: no triangulation in
// between, so the winding normalization in trimmed_build_ring is what is
// actually asserted.
geom::Result<geom::TrimmedBoundary> build_ring_direct(
    const geom::FaceUvLoop& loop, const geom::detail::FaceUvProjection& projection,
    geom::FaceId face) {
  const std::vector<geom::detail::TrimmedRingRun> runs =
      geom::detail::trimmed_ring_runs(loop.arcs);
  const geom::Result<std::size_t> planned = geom::detail::trimmed_plan_ring(runs, 1e-2, face);
  if (!planned) return planned.error();
  return geom::detail::trimmed_build_ring(loop, projection, 1e-2, runs, planned.value(),
                                          loop.outer, face);
}

}  // namespace

TEST_CASE("trimmed domain: a side window leaves one notched outer ring") {
  const geom::Brep brep = cut_cylinder();
  const std::vector<geom::FaceId> faces = cylinder_faces(brep);
  REQUIRE(!faces.empty());

  std::size_t holes = 0;
  double lateral_area = 0.0;
  for (const geom::FaceId id : faces) {
    const geom::Result<geom::detail::FaceUvProjection> projection =
        geom::detail::face_uv_projection(brep, brep.face(id), 1e-9);
    REQUIRE(projection.has_value());
    const geom::Result<std::vector<geom::TrimmedBoundary>> boundaries =
        geom::trimmed_face_boundaries(brep, id, projection.value(), 1e-2);
    REQUIRE(boundaries.has_value());
    REQUIRE(!boundaries.value().empty());
    std::size_t outer = 0;
    for (const geom::TrimmedBoundary& ring : boundaries.value()) {
      CHECK(ring.points.size() >= 3);
      CHECK(ring.arcs.size() == ring.points.size());
      CHECK(geom::signed_area2(ring.points) != 0.0);
      for (std::size_t i = 0; i < ring.points.size(); ++i) {
        const geom::Point2d& a = ring.points[i];
        const geom::Point2d& b = ring.points[(i + 1) % ring.points.size()];
        CHECK(std::abs(b[0] - a[0]) < kPi);  // unwrapped: no jump at the seam
      }
      if (ring.outer) {
        ++outer;
        CHECK(geom::signed_area2(ring.points) > 0.0);
      } else {
        ++holes;
        CHECK(geom::signed_area2(ring.points) < 0.0);
      }
    }
    CHECK(outer == 1);
    lateral_area += polygon_area(boundaries.value());
  }
  CHECK(holes == 0);
  // The seam-crossing window is a notch in the single outer ring: lateral band
  // 4*pi minus the window. y = +-0.5 cuts the unit circle at theta = +-pi/6, so
  // the window spans pi/3 radians and height 1.
  CHECK(lateral_area == doctest::Approx(2.0 * kPi * 1.0 * 2.0 - kPi / 3.0).epsilon(1e-9));
}

TEST_CASE("trimmed domain: deflection adds boundary samples monotonically") {
  const geom::Brep brep = cut_cylinder();
  const geom::FaceId id = cylinder_faces(brep).front();
  const geom::detail::FaceUvProjection projection =
      geom::detail::face_uv_projection(brep, brep.face(id), 1e-9).value();

  const auto count = [&](double deflection) {
    const auto boundaries = geom::trimmed_face_boundaries(brep, id, projection, deflection);
    REQUIRE(boundaries.has_value());
    std::size_t total = 0;
    for (const geom::TrimmedBoundary& ring : boundaries.value()) total += ring.points.size();
    return total;
  };
  CHECK(count(1e-3) > count(1e-2));
  CHECK(count(1e-4) >= count(1e-3));
}

TEST_CASE("trimmed domain: triangulation covers the ring area exactly") {
  const geom::Brep brep = cut_cylinder();
  const geom::FaceId id = cylinder_faces(brep).front();
  const geom::detail::FaceUvProjection projection =
      geom::detail::face_uv_projection(brep, brep.face(id), 1e-9).value();
  const auto boundaries = geom::trimmed_face_boundaries(brep, id, projection, 5e-3);
  REQUIRE(boundaries.has_value());

  const auto triangulation = geom::triangulate_trimmed_domain(boundaries.value());
  REQUIRE(triangulation.has_value());
  CHECK(triangulation.value().points.size() ==
        [&] {
          std::size_t total = 0;
          for (const auto& ring : boundaries.value()) total += ring.points.size();
          return total;
        }());
  REQUIRE(!triangulation.value().triangles.empty());

  double area = 0.0;
  for (const std::array<std::uint32_t, 3>& triangle : triangulation.value().triangles) {
    const geom::Point2d& a = triangulation.value().points[triangle[0]];
    const geom::Point2d& b = triangulation.value().points[triangle[1]];
    const geom::Point2d& c = triangulation.value().points[triangle[2]];
    const double twice = (b[0] - a[0]) * (c[1] - a[1]) - (b[1] - a[1]) * (c[0] - a[0]);
    CHECK(twice > 0.0);
    area += 0.5 * twice;
  }
  CHECK(area == doctest::Approx(polygon_area(boundaries.value())).epsilon(1e-12));
}

TEST_CASE("trimmed domain: triangulation subtracts a hole ring") {
  geom::TrimmedBoundary outer;
  outer.outer = true;
  outer.points = {{0.0, 0.0}, {4.0, 0.0}, {4.0, 4.0}, {0.0, 4.0}};  // CCW
  outer.arcs.resize(outer.points.size());

  geom::TrimmedBoundary hole;
  hole.outer = false;
  hole.points = {{1.0, 1.0}, {1.0, 3.0}, {3.0, 3.0}, {3.0, 1.0}};  // CW
  hole.arcs.resize(hole.points.size());

  const std::vector<geom::TrimmedBoundary> rings{outer, hole};
  const auto triangulation = geom::triangulate_trimmed_domain(rings);
  REQUIRE(triangulation.has_value());
  CHECK(triangulation.value().points.size() == 8);
  double area = 0.0;
  for (const std::array<std::uint32_t, 3>& triangle : triangulation.value().triangles) {
    const geom::Point2d& a = triangulation.value().points[triangle[0]];
    const geom::Point2d& b = triangulation.value().points[triangle[1]];
    const geom::Point2d& c = triangulation.value().points[triangle[2]];
    area += 0.5 * ((b[0] - a[0]) * (c[1] - a[1]) - (b[1] - a[1]) * (c[0] - a[0]));
  }
  CHECK(area == doctest::Approx(12.0).epsilon(1e-12));  // 4x4 minus the 2x2 hole
}

TEST_CASE("trimmed domain: duplicate or non-finite rings are rejected") {
  const geom::Brep brep = cut_cylinder();
  const geom::FaceId id = cylinder_faces(brep).front();
  const geom::detail::FaceUvProjection projection =
      geom::detail::face_uv_projection(brep, brep.face(id), 1e-9).value();
  auto boundaries = geom::trimmed_face_boundaries(brep, id, projection, 1e-2).value();

  std::vector<geom::TrimmedBoundary> duplicated = boundaries;
  duplicated.push_back(boundaries.front());
  CHECK(!geom::triangulate_trimmed_domain(duplicated).has_value());

  std::vector<geom::TrimmedBoundary> non_finite = boundaries;
  non_finite.front().points.front() = geom::Point2d{std::nan(""), 0.0};
  CHECK(!geom::triangulate_trimmed_domain(non_finite).has_value());

  const auto too_fine = geom::trimmed_face_boundaries(brep, id, projection, 1e-12);
  REQUIRE(!too_fine.has_value());
  CHECK(too_fine.error().code == geom::Errc::numeric_failure);
}

TEST_CASE("trimmed domain: extraction normalizes ring winding directly") {
  const geom::Brep brep =
      geom::make_box(geom::Point3d{0.0, 0.0, 0.0}, geom::Point3d{4.0, 4.0, 4.0}, 1e-9).value();
  const geom::FaceId face = brep.face_ids().front();
  const geom::detail::FaceUvProjection projection = uv_projection(brep, face);

  // Clockwise in the face's own uv frame: the outer ring must come back CCW.
  geom::FaceUvLoop clockwise_outer;
  clockwise_outer.outer = true;
  clockwise_outer.points = {{0.0, 0.0}, {0.0, 1.0}, {1.0, 1.0}, {1.0, 0.0}};
  clockwise_outer.arcs.resize(clockwise_outer.points.size());
  const geom::Result<geom::TrimmedBoundary> outer_ring =
      build_ring_direct(clockwise_outer, projection, face);
  REQUIRE(outer_ring.has_value());
  CHECK(outer_ring.value().points.size() == 4);
  CHECK(geom::ring_orientation2(outer_ring.value().points) > 0);
  CHECK(geom::signed_area2(outer_ring.value().points) > 0.0);

  // Counter-clockwise hole input must come back clockwise.
  geom::FaceUvLoop counter_clockwise_hole;
  counter_clockwise_hole.outer = false;
  counter_clockwise_hole.points = {{0.0, 0.0}, {1.0, 0.0}, {1.0, 1.0}, {0.0, 1.0}};
  counter_clockwise_hole.arcs.resize(counter_clockwise_hole.points.size());
  const geom::Result<geom::TrimmedBoundary> hole_ring =
      build_ring_direct(counter_clockwise_hole, projection, face);
  REQUIRE(hole_ring.has_value());
  CHECK(hole_ring.value().points.size() == 4);
  CHECK(geom::ring_orientation2(hole_ring.value().points) < 0);
  CHECK(geom::signed_area2(hole_ring.value().points) < 0.0);
}

TEST_CASE("trimmed domain: extraction rejects zero-area and collapsing rings") {
  const geom::Brep brep =
      geom::make_box(geom::Point3d{0.0, 0.0, 0.0}, geom::Point3d{4.0, 4.0, 4.0}, 1e-9).value();
  const geom::FaceId face = brep.face_ids().front();
  const geom::detail::FaceUvProjection projection = uv_projection(brep, face);

  geom::FaceUvLoop collinear;
  collinear.outer = true;
  collinear.points = {{0.0, 0.0}, {1.0, 1.0}, {2.0, 2.0}};
  collinear.arcs.resize(collinear.points.size());
  const geom::Result<geom::TrimmedBoundary> zero_area =
      build_ring_direct(collinear, projection, face);
  REQUIRE(!zero_area.has_value());
  CHECK(zero_area.error().code == geom::Errc::degenerate);

  // Three identical vertices collapse in the duplicate pass and leave fewer
  // than three points.
  geom::FaceUvLoop collapsing;
  collapsing.outer = true;
  collapsing.points = {{0.0, 0.0}, {0.0, 0.0}, {0.0, 0.0}};
  collapsing.arcs.resize(collapsing.points.size());
  const geom::Result<geom::TrimmedBoundary> collapsed =
      build_ring_direct(collapsing, projection, face);
  REQUIRE(!collapsed.has_value());
  CHECK(collapsed.error().code == geom::Errc::invalid_topology);
}

TEST_CASE("trimmed domain: full-period ring closes across the seam modulo 2*pi") {
  const geom::Brep brep =
      geom::make_cylinder(geom::Point3d{0.0, 0.0, 0.0}, geom::Vec3d{0.0, 0.0, 1.0}, 1.0, 2.0, 1e-9)
          .value();
  geom::FaceId lateral = brep.face_ids().front();
  for (const geom::FaceId id : brep.face_ids()) {
    const geom::Face& face = brep.face(id);
    if (std::holds_alternative<geom::Cylinder3>(brep.geometry().surfaces[face.surface])) {
      lateral = id;
      break;
    }
  }
  const geom::detail::FaceUvProjection projection = uv_projection(brep, lateral);

  // Synthetic full-period loop on the barrel that starts at u = 0.5 (off the
  // seam); its closing vertex is the front vertex one period later, so only the
  // modulo-2*pi closure path can accept it.
  geom::FaceUvLoop loop;
  loop.outer = true;
  loop.points = {{0.5, 0.0}, {2.5, 0.5}, {4.5, 0.5}, {6.5, 0.5}};
  loop.arcs.resize(loop.points.size());
  const geom::Result<geom::TrimmedBoundary> ring = build_ring_direct(loop, projection, lateral);
  REQUIRE(ring.has_value());
  CHECK(ring.value().points.size() == 4);
  CHECK(geom::ring_orientation2(ring.value().points) > 0);
  for (std::size_t i = 0; i < ring.value().points.size(); ++i) {
    const geom::Point2d& a = ring.value().points[i];
    const geom::Point2d& b = ring.value().points[(i + 1) % ring.value().points.size()];
    CHECK(std::abs(a[0] - b[0]) + std::abs(a[1] - b[1]) > 1e-12);  // no zero-length segment
  }
}

TEST_CASE("trimmed domain: opposite arc sweeps do not merge into one run") {
  const geom::CircleCurve3 circle{geom::Point3d{0.0, 0.0, 0.0}, geom::Vec3d{1.0, 0.0, 0.0},
                                  geom::Vec3d{0.0, 1.0, 0.0}, 1.0};

  std::vector<std::optional<geom::FaceUvArcRef>> same_direction(2);
  same_direction[0] = geom::FaceUvArcRef{circle, 0.0, 0.5};
  same_direction[1] = geom::FaceUvArcRef{circle, 0.5, 1.0};
  CHECK(geom::detail::trimmed_ring_runs(same_direction).size() == 1);

  std::vector<std::optional<geom::FaceUvArcRef>> opposite(2);
  opposite[0] = geom::FaceUvArcRef{circle, 0.0, 0.5};
  opposite[1] = geom::FaceUvArcRef{circle, 0.5, 0.2};  // retraces the first span
  CHECK(geom::detail::trimmed_ring_runs(opposite).size() == 2);
}

TEST_CASE("trimmed domain: triangulation rejects duplicate holes and a missing outer") {
  geom::TrimmedBoundary outer;
  outer.outer = true;
  outer.points = {{0.0, 0.0}, {4.0, 0.0}, {4.0, 4.0}, {0.0, 4.0}};
  outer.arcs.resize(outer.points.size());

  geom::TrimmedBoundary hole;
  hole.outer = false;
  hole.points = {{1.0, 1.0}, {1.0, 3.0}, {3.0, 3.0}, {3.0, 1.0}};
  hole.arcs.resize(hole.points.size());

  // Two identical holes must be rejected by the duplicate-ring check, not by
  // the outer-count check.
  const std::vector<geom::TrimmedBoundary> duplicate_holes{outer, hole, hole};
  const auto duplicated = geom::triangulate_trimmed_domain(duplicate_holes);
  REQUIRE(!duplicated.has_value());
  CHECK(duplicated.error().code == geom::Errc::degenerate);

  // A hole with no outer ring is an invalid domain.
  const std::vector<geom::TrimmedBoundary> hole_only{hole};
  const auto missing_outer = geom::triangulate_trimmed_domain(hole_only);
  REQUIRE(!missing_outer.has_value());
  CHECK(missing_outer.error().code == geom::Errc::invalid_topology);
}

TEST_CASE("trimmed domain: planar plate extracts outer and hole rings in order") {
  const geom::Brep brep = plate_with_hole();
  std::size_t plate_faces = 0;
  for (const geom::FaceId id : brep.face_ids()) {
    const geom::Result<std::vector<geom::FaceUvLoop>> loops = geom::face_loops_uv(brep, id, 1e-9);
    REQUIRE(loops.has_value());
    if (loops.value().size() != 2) continue;
    ++plate_faces;

    const geom::detail::FaceUvProjection projection =
        geom::detail::face_uv_projection(brep, brep.face(id), 1e-9).value();
    const auto boundaries = geom::trimmed_face_boundaries(brep, id, projection, 1e-2);
    REQUIRE(boundaries.has_value());
    REQUIRE(boundaries.value().size() == 2);
    CHECK(boundaries.value()[0].outer);
    CHECK(!boundaries.value()[1].outer);
    CHECK(boundaries.value()[0].points.size() == 4);
    CHECK(boundaries.value()[1].points.size() == 4);
    CHECK(geom::ring_orientation2(boundaries.value()[0].points) > 0);
    CHECK(geom::ring_orientation2(boundaries.value()[1].points) < 0);
    CHECK(boundaries.value()[0].arcs.size() == boundaries.value()[0].points.size());
    CHECK(boundaries.value()[1].arcs.size() == boundaries.value()[1].points.size());

    const auto triangulation = geom::triangulate_trimmed_domain(boundaries.value());
    REQUIRE(triangulation.has_value());
    CHECK(triangulation.value().points.size() == 8);
    CHECK(triangle_area(triangulation.value()) == doctest::Approx(15.0).epsilon(1e-12));
  }
  CHECK(plate_faces == 2);
}
