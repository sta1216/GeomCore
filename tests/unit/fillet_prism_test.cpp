#include <doctest/doctest.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <map>
#include <numbers>
#include <span>
#include <string>
#include <utility>
#include <vector>

#include <geom/brep/assemble.hpp>
#include <geom/brep/fillet.hpp>
#include <geom/brep/primitives.hpp>
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

// Same extrusion as make_prism, but in the orthonormal frame
// (u, cross(axis, u), axis) rooted at origin.
geom::Result<geom::Brep> make_prism_in_frame(std::span<const geom::Point2d> polygon,
                                             geom::Point3d origin, geom::Vec3d axis,
                                             geom::Vec3d u, double height, double tolerance) {
  axis = geom::detail::assemble_unit(axis);
  u = geom::detail::assemble_unit(u);
  const geom::Vec3d v = geom::detail::assemble_unit(cross(axis, u));
  const auto place = [&](const geom::Point2d& point, double lift) {
    return geom::Point3d{origin[0] + u[0] * point[0] + v[0] * point[1] + axis[0] * lift,
                         origin[1] + u[1] * point[0] + v[1] * point[1] + axis[1] * lift,
                         origin[2] + u[2] * point[0] + v[2] * point[1] + axis[2] * lift};
  };

  std::vector<geom::Point3d> bottom;
  std::vector<geom::Point3d> top;
  for (const geom::Point2d& point : polygon) {
    bottom.push_back(place(point, 0.0));
    top.push_back(place(point, height));
  }

  std::vector<geom::AssembledFace> faces;
  geom::AssembledFace bottom_face;
  bottom_face.outer.assign(bottom.rbegin(), bottom.rend());
  bottom_face.normal = axis * -1.0;
  faces.push_back(std::move(bottom_face));

  geom::AssembledFace top_face;
  top_face.outer = top;
  top_face.normal = axis;
  faces.push_back(std::move(top_face));

  for (std::size_t i = 0; i < polygon.size(); ++i) {
    const std::size_t j = (i + 1) % polygon.size();
    const geom::Vec3d edge = bottom[j] - bottom[i];
    geom::AssembledFace side;
    side.normal = geom::detail::assemble_unit(cross(edge, axis));
    side.outer = {bottom[i], bottom[j], top[j], top[i]};
    faces.push_back(std::move(side));
  }
  return geom::assemble_solid(faces, tolerance);
}

// Edges of the prism that run parallel to its axis (the filletable ones in any
// frame; vertical_edges assumes the axis is +Z).
std::vector<geom::EdgeId> axis_edges(const geom::Brep& brep, const geom::Vec3d& axis) {
  std::vector<geom::EdgeId> result;
  for (const geom::EdgeId edge_id : brep.edge_ids()) {
    const geom::Edge& edge = brep.edge(edge_id);
    const geom::Vec3d direction = geom::detail::assemble_unit(
        brep.vertex(edge.v1).point - brep.vertex(edge.v0).point);
    if (std::abs(geom::dot(direction, axis)) > 1.0 - 1e-9) result.push_back(edge_id);
  }
  return result;
}

bool same_positions(const std::vector<geom::Point3d>& a,
                    const std::vector<geom::Point3d>& b) {
  if (a.size() != b.size()) return false;
  for (std::size_t i = 0; i < a.size(); ++i) {
    if (a[i][0] != b[i][0] || a[i][1] != b[i][1] || a[i][2] != b[i][2]) return false;
  }
  return true;
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
  CHECK(geom::validate_topology(result.value().brep).ok());
  CHECK(result.value().brep.solid_ids().size() == 1u);
  CHECK(result.value().brep.shell_count() == 1u);
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
  CHECK(geom::validate_topology(result.value().brep).ok());
  CHECK(result.value().brep.solid_ids().size() == 1u);
  CHECK(result.value().brep.shell_count() == 1u);
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

TEST_CASE("fillet_edges enforces the prism face-width guard") {
  const geom::Result<geom::Brep> prism = make_prism(triangle(), 0.0, 3.0, kTolerance);
  REQUIRE(prism.has_value());
  const geom::Brep& solid = prism.value();
  const std::vector<geom::EdgeId> edges{
      vertical_edge_at(solid, geom::Point2d{0.0, 0.0}),
      vertical_edge_at(solid, geom::Point2d{2.0, 0.0})};
  REQUIRE(edges[0] != geom::EdgeId{});
  REQUIRE(edges[1] != geom::EdgeId{});

  // Each base corner keeps a tangent length t = radius * sqrt(3), so the two
  // fillets leave a straight strip of 2 - 2 * t = 5e-7 on the base side. A
  // tolerance at or above the strip rejects the radius; one below it accepts.
  const double radius = (2.0 - 5e-7) / (2.0 * std::sqrt(3.0));

  const geom::Result<geom::FilletResult> rejected =
      geom::fillet_edges(solid, edges, geom::FilletOptions{radius, 1e-6});
  REQUIRE(!rejected.has_value());
  CHECK(rejected.error().code == geom::Errc::unsupported);
  CHECK(std::string(rejected.error().message) ==
        "fillet_edges: radius leaves a prism face narrower than the tolerance");

  const geom::Result<geom::FilletResult> accepted =
      geom::fillet_edges(solid, edges, geom::FilletOptions{radius, 1e-7});
  REQUIRE(accepted.has_value());
  CHECK(accepted.value().filleted_edges == 2u);
  CHECK(geom::validate_topology(accepted.value().brep).ok());
  CHECK(accepted.value().brep.solid_ids().size() == 1u);
  CHECK(accepted.value().brep.shell_count() == 1u);
  CHECK(volume_of(accepted.value().brep) ==
        doctest::Approx(3.82644659481533).epsilon(1e-9));
  CHECK(accepted.value().brep.vertex_count() == 10u);
  CHECK(accepted.value().brep.edge_count() == 15u);
  CHECK(accepted.value().brep.face_count() == 7u);
}

TEST_CASE("fillet_edges rejects a non-prism with the combined message") {
  const geom::Result<geom::Brep> cylinder = geom::make_cylinder(
      geom::Point3d{0.0, 0.0, 0.0}, geom::Vec3d{0.0, 0.0, 1.0}, 1.0, 2.0, kTolerance);
  REQUIRE(cylinder.has_value());
  const std::vector<geom::EdgeId> edges{cylinder.value().edge_ids().front()};
  const geom::Result<geom::FilletResult> result =
      geom::fillet_edges(cylinder.value(), edges, geom::FilletOptions{0.2, kTolerance});
  REQUIRE(!result.has_value());
  CHECK(result.error().code == geom::Errc::unsupported);
  CHECK(std::string(result.error().message) ==
        "fillet_edges: only a box or a convex prism is supported");
}

TEST_CASE("fillet_edges rejects a horizontal prism edge") {
  const geom::Result<geom::Brep> prism = make_prism(triangle(), 0.0, 3.0, kTolerance);
  REQUIRE(prism.has_value());
  geom::EdgeId horizontal{};
  for (const geom::EdgeId edge_id : prism.value().edge_ids()) {
    const geom::Edge& edge = prism.value().edge(edge_id);
    if (std::abs(prism.value().vertex(edge.v0).point[2] - 3.0) <= kTolerance &&
        std::abs(prism.value().vertex(edge.v1).point[2] - 3.0) <= kTolerance) {
      horizontal = edge_id;
      break;
    }
  }
  REQUIRE(horizontal != geom::EdgeId{});
  const std::vector<geom::EdgeId> edges{horizontal};
  const geom::Result<geom::FilletResult> result =
      geom::fillet_edges(prism.value(), edges, geom::FilletOptions{0.2, kTolerance});
  REQUIRE(!result.has_value());
  CHECK(result.error().code == geom::Errc::unsupported);
  CHECK(std::string(result.error().message) ==
        "fillet_edges: only vertical prism edges are supported");
}

TEST_CASE("fillet_edges rejects a radius that leaves a prism face too narrow") {
  const geom::Result<geom::Brep> prism = make_prism(triangle(), 0.0, 3.0, kTolerance);
  REQUIRE(prism.has_value());
  const std::vector<geom::EdgeId> edges = vertical_edges(prism.value());
  REQUIRE(edges.size() == 3u);
  const geom::Result<geom::FilletResult> result =
      geom::fillet_edges(prism.value(), edges, geom::FilletOptions{1.0, kTolerance});
  REQUIRE(!result.has_value());
  CHECK(result.error().code == geom::Errc::unsupported);
  CHECK(std::string(result.error().message) ==
        "fillet_edges: radius leaves a prism face narrower than the tolerance");
}

TEST_CASE("fillet_edges rejects a radius above the prism cross-section") {
  const geom::Result<geom::Brep> prism = make_prism(triangle(), 0.0, 3.0, kTolerance);
  REQUIRE(prism.has_value());
  const std::vector<geom::EdgeId> edges = vertical_edges(prism.value());
  REQUIRE(edges.size() == 3u);
  const geom::Result<geom::FilletResult> result =
      geom::fillet_edges(prism.value(), edges, geom::FilletOptions{1.2, kTolerance});
  REQUIRE(!result.has_value());
  CHECK(result.error().code == geom::Errc::unsupported);
  CHECK(std::string(result.error().message) ==
        "fillet_edges: radius exceeds the prism cross-section");
}

TEST_CASE("fillet_edges rounds one prism corner beyond the polygon inradius") {
  const geom::Result<geom::Brep> prism = make_prism(triangle(), 0.0, 3.0, kTolerance);
  REQUIRE(prism.has_value());
  const std::array<geom::EdgeId, 1> edges{
      vertical_edge_at(prism.value(), geom::Point2d{0.0, 0.0})};
  const geom::Result<geom::FilletResult> result =
      geom::fillet_edges(prism.value(), edges, geom::FilletOptions{0.6, kTolerance});
  REQUIRE(result.has_value());
  CHECK(result.value().filleted_edges == 1u);
  CHECK(volume_of(result.value().brep) ==
        doctest::Approx(4.45651090582474).epsilon(1e-9));
  CHECK(result.value().brep.vertex_count() == 8u);
  CHECK(result.value().brep.edge_count() == 12u);
  CHECK(result.value().brep.face_count() == 6u);
  CHECK(geom::validate_topology(result.value().brep).ok());
  CHECK(result.value().brep.solid_ids().size() == 1u);
  CHECK(result.value().brep.shell_count() == 1u);
}

TEST_CASE("fillet_edges rejects a dead edge and an empty selection") {
  const geom::Result<geom::Brep> prism = make_prism(triangle(), 0.0, 3.0, kTolerance);
  REQUIRE(prism.has_value());

  const geom::Result<geom::FilletResult> empty = geom::fillet_edges(
      prism.value(), std::span<const geom::EdgeId>{}, geom::FilletOptions{0.2, kTolerance});
  REQUIRE(!empty.has_value());
  CHECK(empty.error().code == geom::Errc::unsupported);
  CHECK(std::string(empty.error().message) == "fillet_edges: at least one edge is required");

  const std::array<geom::EdgeId, 1> dead{geom::EdgeId::make(999u, 0u)};
  const geom::Result<geom::FilletResult> dead_result =
      geom::fillet_edges(prism.value(), dead, geom::FilletOptions{0.2, kTolerance});
  REQUIRE(!dead_result.has_value());
  CHECK(dead_result.error().code == geom::Errc::index_out_of_range);
  CHECK(std::string(dead_result.error().message) ==
        "fillet_edges: edge does not belong to the solid");
}

TEST_CASE("fillet_edges rounds a prism in a tilted frame") {
  const geom::Vec3d axis = geom::detail::assemble_unit(geom::Vec3d{1.0, 1.0, 1.0});
  const geom::Vec3d u = geom::detail::assemble_unit(cross(axis, geom::Vec3d{0.0, 0.0, 1.0}));
  const geom::Result<geom::Brep> prism = make_prism_in_frame(
      triangle(), geom::Point3d{5.0, -2.0, 7.0}, axis, u, 3.0, kTolerance);
  REQUIRE(prism.has_value());
  const std::vector<geom::EdgeId> edges = axis_edges(prism.value(), axis);
  REQUIRE(edges.size() == 3u);
  const geom::Result<geom::FilletResult> result =
      geom::fillet_edges(prism.value(), edges, geom::FilletOptions{0.2, kTolerance});
  REQUIRE(result.has_value());
  CHECK(volume_of(result.value().brep) ==
        doctest::Approx(4.949605250412612).epsilon(1e-9));
  CHECK(result.value().brep.vertex_count() == 12u);
  CHECK(result.value().brep.edge_count() == 18u);
  CHECK(result.value().brep.face_count() == 8u);
  CHECK(geom::validate_topology(result.value().brep).ok());
  CHECK(result.value().brep.solid_ids().size() == 1u);
  CHECK(result.value().brep.shell_count() == 1u);
}

TEST_CASE("fillet_edges leaves its input unchanged and is deterministic") {
  const geom::Result<geom::Brep> prism = make_prism(triangle(), 0.0, 3.0, kTolerance);
  REQUIRE(prism.has_value());
  const geom::Brep& solid = prism.value();
  const std::vector<geom::EdgeId> edges = vertical_edges(solid);
  REQUIRE(edges.size() == 3u);

  const std::size_t vertex_count = solid.vertex_count();
  const std::size_t edge_count = solid.edge_count();
  const std::size_t face_count = solid.face_count();
  const double volume = volume_of(solid);

  const geom::Result<geom::FilletResult> first =
      geom::fillet_edges(solid, edges, geom::FilletOptions{0.2, kTolerance});
  const geom::Result<geom::FilletResult> second =
      geom::fillet_edges(solid, edges, geom::FilletOptions{0.2, kTolerance});
  REQUIRE(first.has_value());
  REQUIRE(second.has_value());

  CHECK(solid.vertex_count() == vertex_count);
  CHECK(solid.edge_count() == edge_count);
  CHECK(solid.face_count() == face_count);
  CHECK(volume_of(solid) == volume);

  CHECK(first.value().brep.vertex_count() == second.value().brep.vertex_count());
  CHECK(first.value().brep.edge_count() == second.value().brep.edge_count());
  CHECK(first.value().brep.face_count() == second.value().brep.face_count());
  CHECK(volume_of(first.value().brep) == volume_of(second.value().brep));

  const geom::Result<geom::TriMesh> first_mesh =
      geom::tessellate(first.value().brep, geom::TessellationOptions{});
  const geom::Result<geom::TriMesh> second_mesh =
      geom::tessellate(second.value().brep, geom::TessellationOptions{});
  REQUIRE(first_mesh.has_value());
  REQUIRE(second_mesh.has_value());
  CHECK(same_positions(first_mesh.value().positions, second_mesh.value().positions));
  CHECK(first_mesh.value().triangles == second_mesh.value().triangles);

  const geom::Result<geom::FilletResult> rejected =
      geom::fillet_edges(solid, edges, geom::FilletOptions{1.0, kTolerance});
  REQUIRE(!rejected.has_value());
  CHECK(solid.vertex_count() == vertex_count);
  CHECK(solid.edge_count() == edge_count);
  CHECK(solid.face_count() == face_count);
  CHECK(volume_of(solid) == volume);
}
