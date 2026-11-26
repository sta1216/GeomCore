#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <map>
#include <numbers>
#include <span>
#include <string_view>
#include <utility>
#include <variant>
#include <vector>

#include <doctest/doctest.h>

#include <geom/brep/detail/curved_volume.hpp>
#include <geom/brep/fillet.hpp>
#include <geom/brep/primitives.hpp>
#include <geom/brep/tessellate.hpp>
#include <geom/brep/topology.hpp>
#include <geom/brep/validate.hpp>
#include <geom/core/error.hpp>
#include <geom/curves/analytic_curve.hpp>
#include <geom/surfaces/analytic.hpp>
#include <geom/surfaces/plane.hpp>

namespace {

using geom::Brep;
using geom::EdgeId;
using geom::Errc;
using geom::Point3d;
using geom::Vec3d;

Brep box_2x3x4() {
  return geom::make_box(Point3d{0.0, 0.0, 0.0}, Point3d{2.0, 3.0, 4.0}, 1e-9).value();
}

double volume_of(const Brep& brep) {
  return geom::detail::curved_solid_volume(brep, brep.solid_ids().front(), 1e-9).value();
}

// Exact coordinate equality: two runs of the same deterministic code must
// produce bit-identical vertices.
bool exactly_equal(const Point3d& a, const Point3d& b) {
  return a[0] == b[0] && a[1] == b[1] && a[2] == b[2];
}

EdgeId vertical_edge_at(const Brep& brep, double x, double y) {
  for (const EdgeId id : brep.edge_ids()) {
    const geom::Edge& edge = brep.edge(id);
    const Point3d a = brep.vertex(edge.v0).point;
    const Point3d b = brep.vertex(edge.v1).point;
    if (std::abs(a[0] - x) <= 1e-12 && std::abs(a[1] - y) <= 1e-12 &&
        std::abs(b[0] - x) <= 1e-12 && std::abs(b[1] - y) <= 1e-12 &&
        std::abs(a[2] - 0.0) <= 1e-12 && std::abs(b[2] - 4.0) <= 1e-12) {
      return id;
    }
  }
  return EdgeId{};
}

// Horizontal family edge at the given (y, z): one endpoint at x=0, the other at
// x=2.
EdgeId horizontal_edge_at(const Brep& brep, double z, double y) {
  for (const EdgeId id : brep.edge_ids()) {
    const geom::Edge& edge = brep.edge(id);
    const Point3d a = brep.vertex(edge.v0).point;
    const Point3d b = brep.vertex(edge.v1).point;
    const bool spans_x = (std::abs(a[0] - 0.0) <= 1e-12 && std::abs(b[0] - 2.0) <= 1e-12) ||
                         (std::abs(a[0] - 2.0) <= 1e-12 && std::abs(b[0] - 0.0) <= 1e-12);
    if (spans_x && std::abs(a[1] - y) <= 1e-12 && std::abs(b[1] - y) <= 1e-12 &&
        std::abs(a[2] - z) <= 1e-12 && std::abs(b[2] - z) <= 1e-12) {
      return id;
    }
  }
  return EdgeId{};
}

// Counts the undirected triangle edges that occur exactly once. A closed
// triangle mesh has none.
std::size_t boundary_edges(const geom::TriMesh& mesh) {
  std::map<std::pair<std::uint32_t, std::uint32_t>, int> uses;
  for (const auto& triangle : mesh.triangles) {
    for (int i = 0; i < 3; ++i) {
      std::uint32_t a = triangle[static_cast<std::size_t>(i)];
      std::uint32_t b = triangle[(i + 1) % 3];
      if (a > b) std::swap(a, b);
      ++uses[{a, b}];
    }
  }
  std::size_t boundary = 0;
  for (const auto& [key, count] : uses) {
    (void)key;
    if (count != 2) ++boundary;
  }
  return boundary;
}

std::vector<EdgeId> family_edges(const Brep& brep) {
  std::vector<EdgeId> result;
  for (const std::array<double, 2>& corner : {std::array{0.0, 0.0}, std::array{2.0, 0.0},
                                              std::array{2.0, 3.0}, std::array{0.0, 3.0}}) {
    result.push_back(vertical_edge_at(brep, corner[0], corner[1]));
  }
  return result;
}

// Edge spanning the whole box along `axis` at the given coordinates on the two
// cross axes (passed in ascending axis order).
EdgeId parallel_edge_at(const Brep& brep, int axis, double first, double second) {
  const int axis_u = (axis + 1) % 3;
  const int axis_v = (axis + 2) % 3;
  for (const EdgeId id : brep.edge_ids()) {
    const geom::Edge& edge = brep.edge(id);
    const Point3d a = brep.vertex(edge.v0).point;
    const Point3d b = brep.vertex(edge.v1).point;
    if (std::abs(a[axis_u] - first) <= 1e-12 && std::abs(b[axis_u] - first) <= 1e-12 &&
        std::abs(a[axis_v] - second) <= 1e-12 && std::abs(b[axis_v] - second) <= 1e-12 &&
        std::abs(b[axis] - a[axis]) > 0.0) {
      return id;
    }
  }
  return EdgeId{};
}

std::vector<EdgeId> family_edges_on(const Brep& brep, int axis) {
  const int axis_u = (axis + 1) % 3;
  const int axis_v = (axis + 2) % 3;
  const double extent[3] = {2.0, 3.0, 4.0};
  std::vector<EdgeId> result;
  for (const double first : {0.0, extent[axis_u]}) {
    for (const double second : {0.0, extent[axis_v]}) {
      result.push_back(parallel_edge_at(brep, axis, first, second));
    }
  }
  return result;
}

}  // namespace

TEST_CASE("fillet_edges rejects options and non-box shapes") {
  const Brep box = box_2x3x4();
  const std::span<const EdgeId> empty;
  CHECK(geom::fillet_edges(box, empty, {}).error().code == Errc::unsupported);

  const EdgeId edge = vertical_edge_at(box, 2.0, 3.0);
  const std::span<const EdgeId> one(&edge, 1);
  CHECK(geom::fillet_edges(box, one, {0.0, 1e-9}).error().code == Errc::degenerate);
  CHECK(geom::fillet_edges(box, one,
                           {std::numeric_limits<double>::quiet_NaN(), 1e-9})
            .error()
            .code == Errc::non_finite);
  CHECK(geom::fillet_edges(box, one, {0.25, 0.0}).error().code ==
        Errc::tolerance_violation);

  const Brep cylinder = geom::make_cylinder(Point3d{0.0, 0.0, 0.0}, Vec3d{0.0, 0.0, 1.0},
                                            1.0, 2.0, 1e-9)
                            .value();
  const std::array<EdgeId, 1> dead{EdgeId::make(999u, 0u)};
  CHECK(geom::fillet_edges(cylinder, dead, {0.25, 1e-9}).error().code == Errc::unsupported);
}

TEST_CASE("fillet_edges rounds one vertical box edge exactly") {
  const Brep box = box_2x3x4();
  const EdgeId edge = vertical_edge_at(box, 2.0, 3.0);
  REQUIRE(edge.index() != 0xFFFFFFFFu);

  const auto result = geom::fillet_edges(box, std::span<const EdgeId>(&edge, 1), {0.25, 1e-9});
  REQUIRE(result.has_value());
  const Brep& filleted = result.value().brep;
  CHECK(result.value().filleted_edges == 1u);
  CHECK(geom::validate_topology(filleted).ok());
  CHECK(filleted.vertex_count() == 10u);
  CHECK(filleted.edge_count() == 15u);
  CHECK(filleted.face_count() == 7u);
  CHECK(std::abs(volume_of(filleted) - 23.94634954084936) <= 1e-9);

  std::size_t cylinders = 0;
  std::size_t arcs = 0;
  for (const geom::FaceId face : filleted.face_ids()) {
    const auto& surface =
        filleted.geometry().surfaces[filleted.face(face).surface];
    if (std::holds_alternative<geom::Cylinder3>(surface)) {
      ++cylinders;
      CHECK(std::abs(std::get<geom::Cylinder3>(surface).radius - 0.25) <= 1e-12);
    }
  }
  for (const EdgeId id : filleted.edge_ids()) {
    const auto& curve = filleted.geometry().curves[filleted.edge(id).curve];
    if (std::holds_alternative<geom::CircleCurve3>(curve)) {
      ++arcs;
      CHECK(std::abs(std::get<geom::CircleCurve3>(curve).radius - 0.25) <= 1e-12);
    }
  }
  CHECK(cylinders == 1u);
  CHECK(arcs == 2u);
}

TEST_CASE("fillet_edges rounds every non-empty subset of every box family") {
  const Brep box = box_2x3x4();
  const double pi = std::numbers::pi_v<double>;
  const double extent[3] = {2.0, 3.0, 4.0};

  for (int axis = 0; axis < 3; ++axis) {
    const std::vector<EdgeId> all = family_edges_on(box, axis);
    REQUIRE(all.size() == 4u);
    for (const EdgeId id : all) REQUIRE(id.index() != 0xFFFFFFFFu);

    for (unsigned mask = 1; mask < 16u; ++mask) {
      std::vector<EdgeId> selected;
      for (int bit = 0; bit < 4; ++bit) {
        if ((mask & (1u << bit)) != 0u) selected.push_back(all[static_cast<std::size_t>(bit)]);
      }
      const auto k = static_cast<double>(selected.size());
      const double per_edge = (1.0 - pi / 4.0) * 0.25 * 0.25 * extent[axis];
      const auto result = geom::fillet_edges(box, selected, {0.25, 1e-9});
      REQUIRE(result.has_value());
      const Brep& filleted = result.value().brep;
      CHECK(result.value().filleted_edges == selected.size());
      CHECK(filleted.vertex_count() == 8u + 2u * selected.size());
      CHECK(filleted.edge_count() == 12u + 3u * selected.size());
      CHECK(filleted.face_count() == 6u + selected.size());
      CHECK(std::abs(volume_of(filleted) - (24.0 - k * per_edge)) <= 1e-9);
      CHECK(geom::validate_topology(filleted).ok());

      std::size_t cylinders = 0;
      for (const geom::FaceId face : filleted.face_ids()) {
        const auto& surface = filleted.geometry().surfaces[filleted.face(face).surface];
        if (!std::holds_alternative<geom::Cylinder3>(surface)) continue;
        const geom::Cylinder3& cylinder = std::get<geom::Cylinder3>(surface);
        ++cylinders;
        CHECK(std::abs(cylinder.radius - 0.25) <= 1e-12);
        CHECK(std::abs(cylinder.origin[axis] - 0.0) <= 1e-12);
        for (int other = 0; other < 3; ++other) {
          if (other == axis) continue;
          const double coordinate = cylinder.origin[other];
          const bool at_low = std::abs(coordinate - 0.25) <= 1e-12;
          const bool at_high = std::abs(coordinate - (extent[other] - 0.25)) <= 1e-12;
          CHECK((at_low || at_high));
        }
      }
      CHECK(cylinders == selected.size());
    }
  }
}

TEST_CASE("fillet_edges rounds the four vertical prefixes") {
  const Brep box = box_2x3x4();
  const std::vector<EdgeId> all = family_edges(box);
  const double pi = std::numbers::pi_v<double>;
  const double per_edge = (1.0 - pi / 4.0) * 0.25 * 0.25 * 4.0;

  for (std::size_t k = 1; k <= 4; ++k) {
    const std::vector<EdgeId> selected(all.begin(), all.begin() + static_cast<std::ptrdiff_t>(k));
    const auto result = geom::fillet_edges(box, selected, {0.25, 1e-9});
    REQUIRE(result.has_value());
    const Brep& filleted = result.value().brep;
    CHECK(result.value().filleted_edges == k);
    CHECK(filleted.vertex_count() == 8u + 2u * k);
    CHECK(filleted.edge_count() == 12u + 3u * k);
    CHECK(filleted.face_count() == 6u + k);
    CHECK(std::abs(volume_of(filleted) - (24.0 - static_cast<double>(k) * per_edge)) <= 1e-9);
    CHECK(geom::validate_topology(filleted).ok());

    if (k == 4) {
      // Per-corner identity: four cylinder patches, each axis along +z from the
      // family-axis low end, and four distinct cross-section corners inset by
      // the radius from a lo/hi extreme.
      std::vector<Point3d> origins;
      std::size_t cylinders = 0;
      for (const geom::FaceId face : filleted.face_ids()) {
        const auto& surface =
            filleted.geometry().surfaces[filleted.face(face).surface];
        if (!std::holds_alternative<geom::Cylinder3>(surface)) continue;
        const geom::Cylinder3& cylinder = std::get<geom::Cylinder3>(surface);
        ++cylinders;
        CHECK(std::abs(cylinder.radius - 0.25) <= 1e-12);
        CHECK(std::abs(cylinder.axis[0]) <= 1e-12);
        CHECK(std::abs(cylinder.axis[1]) <= 1e-12);
        CHECK(std::abs(cylinder.axis[2] - 1.0) <= 1e-12);
        origins.push_back(cylinder.origin);
      }
      CHECK(cylinders == 4u);
      REQUIRE(origins.size() == 4u);
      for (const double x : {0.25, 1.75}) {
        for (const double y : {0.25, 2.75}) {
          bool found = false;
          for (const Point3d& origin : origins) {
            if (std::abs(origin[0] - x) <= 1e-12 && std::abs(origin[1] - y) <= 1e-12 &&
                std::abs(origin[2] - 0.0) <= 1e-12) {
              found = true;
            }
          }
          CHECK(found);
        }
      }
    }
  }
}

TEST_CASE("fillet_edges deduplicates and ignores selection order") {
  const Brep box = box_2x3x4();
  const std::vector<EdgeId> all = family_edges(box);
  std::vector<EdgeId> duplicated{all[0], all[1], all[0], all[2], all[1]};
  const auto forward = geom::fillet_edges(box, duplicated, {0.25, 1e-9});
  REQUIRE(forward.has_value());
  CHECK(forward.value().filleted_edges == 3u);

  std::vector<EdgeId> reversed(all.begin(), all.begin() + 3);
  std::reverse(reversed.begin(), reversed.end());
  const auto backward = geom::fillet_edges(box, reversed, {0.25, 1e-9});
  REQUIRE(backward.has_value());
  CHECK(backward.value().brep.face_count() == forward.value().brep.face_count());
  CHECK(std::abs(volume_of(backward.value().brep) - volume_of(forward.value().brep)) <= 1e-9);
}

TEST_CASE("fillet_edges rejects a radius at the section limit") {
  const Brep box = box_2x3x4();
  const auto edge = vertical_edge_at(box, 2.0, 3.0);
  const auto too_large = geom::fillet_edges(box, std::span<const EdgeId>(&edge, 1), {1.0, 1e-9});
  REQUIRE_FALSE(too_large.has_value());
  CHECK(too_large.error().code == Errc::unsupported);

  const auto just_below =
      geom::fillet_edges(box, std::span<const EdgeId>(&edge, 1), {0.999, 1e-9});
  REQUIRE(just_below.has_value());
  CHECK(geom::validate_topology(just_below.value().brep).ok());

  // A radius closer than the tolerance to the limit leaves a 2e-10 side face,
  // which the assembler cannot represent; reject it with a stable public error
  // instead of leaking the assembler message.
  const std::vector<EdgeId> all = family_edges(box);
  const auto sliver = geom::fillet_edges(box, all, {0.9999999999, 1e-9});
  REQUIRE_FALSE(sliver.has_value());
  CHECK(sliver.error().code == Errc::unsupported);
  CHECK(std::string_view(sliver.error().message).find("narrower than the tolerance") !=
        std::string_view::npos);
}

TEST_CASE("fillet_edges rejects a mixed family and a dead edge") {
  const Brep box = box_2x3x4();
  const EdgeId vertical = vertical_edge_at(box, 2.0, 3.0);
  const EdgeId horizontal = horizontal_edge_at(box, 4.0, 3.0);
  const std::array<EdgeId, 2> mixed{vertical, horizontal};
  const auto mixed_result = geom::fillet_edges(box, mixed, {0.25, 1e-9});
  REQUIRE_FALSE(mixed_result.has_value());
  CHECK(mixed_result.error().code == Errc::unsupported);
  CHECK(std::string_view(mixed_result.error().message).find("parallel") !=
        std::string_view::npos);

  const std::array<EdgeId, 1> dead{EdgeId::make(999u, 0u)};
  const auto dead_result = geom::fillet_edges(box, dead, {0.25, 1e-9});
  REQUIRE_FALSE(dead_result.has_value());
  CHECK(dead_result.error().code == Errc::index_out_of_range);
}

TEST_CASE("fillet_edges is watertight and deterministic") {
  const Brep box = box_2x3x4();
  const std::vector<EdgeId> all = family_edges(box);
  const std::size_t box_vertices = box.vertex_count();
  const std::size_t box_edges = box.edge_count();
  const std::size_t box_faces = box.face_count();
  const double box_volume = volume_of(box);
  const auto check_box_unchanged = [&] {
    CHECK(box.vertex_count() == box_vertices);
    CHECK(box.edge_count() == box_edges);
    CHECK(box.face_count() == box_faces);
    CHECK(std::abs(volume_of(box) - box_volume) <= 1e-12);
  };

  const auto first = geom::fillet_edges(box, all, {0.25, 1e-9});
  REQUIRE(first.has_value());
  check_box_unchanged();

  const auto second = geom::fillet_edges(box, all, {0.25, 1e-9});
  REQUIRE(second.has_value());
  check_box_unchanged();

  // A failing call on the same box must not mutate it either.
  const auto rejected = geom::fillet_edges(box, all, {1.0, 1e-9});
  REQUIRE_FALSE(rejected.has_value());
  check_box_unchanged();

  const auto tessellated =
      geom::tessellate(first.value().brep, geom::TessellationOptions{});
  REQUIRE(tessellated.has_value());
  CHECK(boundary_edges(tessellated.value()) == 0u);

  CHECK(first.value().brep.vertex_count() == second.value().brep.vertex_count());
  CHECK(first.value().brep.edge_count() == second.value().brep.edge_count());
  CHECK(first.value().brep.face_count() == second.value().brep.face_count());
  CHECK(std::abs(volume_of(first.value().brep) - volume_of(second.value().brep)) <= 1e-12);

  // Determinism beyond counts: the vertex slots and their coordinates must
  // match in vertex_ids() order.
  const std::vector<geom::VertexId> first_vertices = first.value().brep.vertex_ids();
  const std::vector<geom::VertexId> second_vertices = second.value().brep.vertex_ids();
  REQUIRE(first_vertices.size() == second_vertices.size());
  for (std::size_t index = 0; index < first_vertices.size(); ++index) {
    CHECK(exactly_equal(first.value().brep.vertex(first_vertices[index]).point,
                        second.value().brep.vertex(second_vertices[index]).point));
  }
}
