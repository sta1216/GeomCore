#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <map>
#include <numbers>
#include <span>
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

std::vector<EdgeId> family_edges(const Brep& brep) {
  std::vector<EdgeId> result;
  for (const std::array<double, 2>& corner : {std::array{0.0, 0.0}, std::array{2.0, 0.0},
                                              std::array{2.0, 3.0}, std::array{0.0, 3.0}}) {
    result.push_back(vertical_edge_at(brep, corner[0], corner[1]));
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

TEST_CASE("fillet_edges rounds every subset of one box family") {
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
