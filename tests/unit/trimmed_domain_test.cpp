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

TEST_CASE("trimmed domain: reversed outer ring still triangulates counter-clockwise") {
  geom::TrimmedBoundary outer;
  outer.outer = true;
  outer.points = {{0.0, 0.0}, {0.0, 4.0}, {4.0, 4.0}, {4.0, 0.0}};  // CW on purpose
  outer.arcs.resize(outer.points.size());

  const std::vector<geom::TrimmedBoundary> rings{outer};
  const auto triangulation = geom::triangulate_trimmed_domain(rings);
  REQUIRE(triangulation.has_value());
  CHECK(triangulation.value().points.size() == 4);
  REQUIRE(!triangulation.value().triangles.empty());
  for (const std::array<std::uint32_t, 3>& triangle : triangulation.value().triangles) {
    const geom::Point2d& a = triangulation.value().points[triangle[0]];
    const geom::Point2d& b = triangulation.value().points[triangle[1]];
    const geom::Point2d& c = triangulation.value().points[triangle[2]];
    CHECK(geom::sign_orient2d(a[0], a[1], b[0], b[1], c[0], c[1]) > 0);
  }
  CHECK(triangle_area(triangulation.value()) == doctest::Approx(16.0).epsilon(1e-12));
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
