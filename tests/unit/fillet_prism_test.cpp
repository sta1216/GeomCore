#include <doctest/doctest.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <map>
#include <numbers>
#include <span>
#include <utility>
#include <vector>

#include <geom/brep/assemble.hpp>
#include <geom/brep/fillet.hpp>
#include <geom/brep/detail/curved_volume.hpp>
#include <geom/brep/tessellate.hpp>
#include <geom/brep/topology.hpp>
#include <geom/brep/validate.hpp>
#include <geom/cgeom/polygon2.hpp>
#include <geom/core/error.hpp>
#include <geom/core/result.hpp>
#include <geom/curves/analytic_curve.hpp>
#include <geom/primitives/point.hpp>
#include <geom/primitives/vec.hpp>

namespace {

constexpr double kTolerance = 1e-9;
constexpr double kPi = std::numbers::pi_v<double>;

std::vector<geom::Point2d> triangle() {
  return {geom::Point2d{0.0, 0.0}, geom::Point2d{2.0, 0.0},
          geom::Point2d{1.0, std::sqrt(3.0)}};
}

// CCW polygon in the z = [z0, z1] extrusion. Bottom is clockwise about -Z;
// sides are outward-facing quads. The topology is built by assemble_solid.
geom::Result<geom::Brep> make_prism(std::span<const geom::Point2d> polygon, double z0, double z1,
                                   double tolerance) {
  std::vector<geom::AssembledFace> faces;
  std::vector<geom::Point3d> bottom;
  std::vector<geom::Point3d> top;
  for (const geom::Point2d& p : polygon) {
    bottom.push_back(geom::Point3d{p[0], p[1], z0});
    top.push_back(geom::Point3d{p[0], p[1], z1});
  }

  geom::AssembledFace bottom_face;
  bottom_face.outer.assign(bottom.rbegin(), bottom.rend());
  bottom_face.normal = geom::Vec3d{0.0, 0.0, -1.0};
  faces.push_back(std::move(bottom_face));

  geom::AssembledFace top_face;
  top_face.outer = top;
  top_face.normal = geom::Vec3d{0.0, 0.0, 1.0};
  faces.push_back(std::move(top_face));

  for (std::size_t i = 0; i < polygon.size(); ++i) {
    const std::size_t j = (i + 1) % polygon.size();
    const geom::Vec3d edge{polygon[j][0] - polygon[i][0], polygon[j][1] - polygon[i][1], 0.0};
    const geom::Vec3d normal = geom::detail::assemble_unit(
        geom::Vec3d{edge[1], -edge[0], 0.0});
    geom::AssembledFace side;
    side.normal = normal;
    side.outer = {bottom[i], bottom[j], top[j], top[i]};
    faces.push_back(std::move(side));
  }
  return geom::assemble_solid(faces, tolerance);
}

std::vector<geom::EdgeId> vertical_edges(const geom::Brep& brep) {
  std::vector<geom::EdgeId> result;
  for (const geom::EdgeId edge_id : brep.edge_ids()) {
    const geom::Edge& edge = brep.edge(edge_id);
    const geom::Point3d& a = brep.vertex(edge.v0).point;
    const geom::Point3d& b = brep.vertex(edge.v1).point;
    if (std::hypot(b[0] - a[0], b[1] - a[1]) <= 1e-9 && std::abs(b[2] - a[2]) > 1e-9) {
      result.push_back(edge_id);
    }
  }
  return result;
}

geom::EdgeId vertical_edge_at(const geom::Brep& brep, geom::Point2d xy) {
  for (const geom::EdgeId edge_id : vertical_edges(brep)) {
    const geom::Edge& edge = brep.edge(edge_id);
    const geom::Point3d& a = brep.vertex(edge.v0).point;
    if (std::hypot(a[0] - xy[0], a[1] - xy[1]) <= 1e-9) return edge_id;
  }
  return geom::EdgeId{};
}

double corner_area(double radius, double theta) {
  return radius * radius *
         (1.0 / std::tan(theta / 2.0) - (kPi - theta) / 2.0);
}

double volume_of(const geom::Brep& brep) {
  return geom::detail::curved_solid_volume(brep, brep.solid_ids().front(), 1e-9).value();
}

}  // namespace

TEST_CASE("fillet_edges rounds one vertical edge of a triangular prism") {
  const geom::Result<geom::Brep> prism = make_prism(triangle(), 0.0, 3.0, kTolerance);
  REQUIRE(prism.has_value());
  const geom::Brep& solid = prism.value();
  const std::vector<geom::EdgeId> edges{vertical_edge_at(solid, geom::Point2d{0.0, 0.0})};
  REQUIRE(edges.front() != geom::EdgeId{});

  const geom::Result<geom::FilletResult> result =
      geom::fillet_edges(solid, edges, geom::FilletOptions{0.2, kTolerance});
  REQUIRE(result.has_value());
  const geom::Brep& rounded = result.value().brep;
  CHECK(result.value().filleted_edges == 1u);
  CHECK(geom::validate_topology(rounded).ok());
  CHECK(rounded.solid_ids().size() == 1u);
  CHECK(rounded.shell_count() == 1u);

  const double expected =
      3.0 * std::sqrt(3.0) - 3.0 * corner_area(0.2, kPi / 3.0);
  CHECK(volume_of(rounded) == doctest::Approx(expected).epsilon(1e-9));
  CHECK(volume_of(rounded) == doctest::Approx(5.113970031941959).epsilon(1e-9));
  CHECK(rounded.vertex_count() == 8u);
  CHECK(rounded.edge_count() == 12u);
  CHECK(rounded.face_count() == 6u);
}

std::vector<geom::Point2d> regular_polygon(int n, double radius) {
  std::vector<geom::Point2d> polygon;
  for (int i = 0; i < n; ++i) {
    const double angle = 2.0 * kPi * static_cast<double>(i) / static_cast<double>(n);
    polygon.push_back(geom::Point2d{radius * std::cos(angle), radius * std::sin(angle)});
  }
  return polygon;
}

TEST_CASE("fillet_edges rounds every vertical edge of a triangular prism") {
  const geom::Result<geom::Brep> prism = make_prism(triangle(), 0.0, 3.0, kTolerance);
  REQUIRE(prism.has_value());
  const std::vector<geom::EdgeId> edges = vertical_edges(prism.value());
  REQUIRE(edges.size() == 3u);
  const geom::Result<geom::FilletResult> result =
      geom::fillet_edges(prism.value(), edges, geom::FilletOptions{0.2, kTolerance});
  REQUIRE(result.has_value());
  CHECK(result.value().filleted_edges == 3u);
  CHECK(volume_of(result.value().brep) ==
        doctest::Approx(3.0 * std::sqrt(3.0) - 9.0 * corner_area(0.2, kPi / 3.0))
            .epsilon(1e-9));
  CHECK(result.value().brep.vertex_count() == 12u);
  CHECK(result.value().brep.edge_count() == 18u);
  CHECK(result.value().brep.face_count() == 8u);
}

TEST_CASE("fillet_edges rounds an adjacent pair on a hexagonal prism") {
  const std::vector<geom::Point2d> polygon = regular_polygon(6, 1.0);
  const geom::Result<geom::Brep> prism = make_prism(polygon, 0.0, 2.0, kTolerance);
  REQUIRE(prism.has_value());
  std::vector<geom::EdgeId> edges = vertical_edges(prism.value());
  REQUIRE(edges.size() == 6u);
  const std::vector<geom::EdgeId> selected{edges[0], edges[1]};
  const geom::Result<geom::FilletResult> result =
      geom::fillet_edges(prism.value(), selected, geom::FilletOptions{0.2, kTolerance});
  REQUIRE(result.has_value());
  const double area = std::abs(geom::signed_area2(polygon));
  const double expected =
      area * 2.0 - 2.0 * 2.0 * corner_area(0.2, 2.0 * kPi / 3.0);
  CHECK(volume_of(result.value().brep) == doctest::Approx(expected).epsilon(1e-9));
  CHECK(result.value().brep.vertex_count() == 16u);
  CHECK(result.value().brep.edge_count() == 24u);
  CHECK(result.value().brep.face_count() == 10u);
}

TEST_CASE("fillet_edges deduplicates and reorders vertical edge selections") {
  const geom::Result<geom::Brep> prism = make_prism(triangle(), 0.0, 3.0, kTolerance);
  REQUIRE(prism.has_value());
  const std::vector<geom::EdgeId> edges = vertical_edges(prism.value());
  REQUIRE(edges.size() == 3u);

  const std::vector<geom::EdgeId> duplicated{edges[2], edges[0], edges[2]};
  const std::vector<geom::EdgeId> unique{edges[0], edges[2]};
  const geom::Result<geom::FilletResult> from_duplicated =
      geom::fillet_edges(prism.value(), duplicated, geom::FilletOptions{0.2, kTolerance});
  const geom::Result<geom::FilletResult> from_unique =
      geom::fillet_edges(prism.value(), unique, geom::FilletOptions{0.2, kTolerance});
  REQUIRE(from_duplicated.has_value());
  REQUIRE(from_unique.has_value());

  CHECK(from_duplicated.value().filleted_edges == 2u);
  CHECK(from_unique.value().filleted_edges == 2u);
  CHECK(volume_of(from_duplicated.value().brep) ==
        doctest::Approx(volume_of(from_unique.value().brep)).epsilon(1e-12));
  const double expected =
      3.0 * std::sqrt(3.0) - 6.0 * corner_area(0.2, kPi / 3.0);
  CHECK(volume_of(from_unique.value().brep) == doctest::Approx(expected).epsilon(1e-9));
  CHECK(from_duplicated.value().brep.vertex_count() == 10u);
  CHECK(from_duplicated.value().brep.edge_count() == 15u);
  CHECK(from_duplicated.value().brep.face_count() == 7u);
  CHECK(geom::validate_topology(from_duplicated.value().brep).ok());
}
