#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <map>
#include <span>
#include <string_view>
#include <utility>
#include <variant>
#include <vector>

#include <doctest/doctest.h>

#include <geom/brep/detail/curved_volume.hpp>
#include <geom/brep/fillet.hpp>
#include <geom/brep/mass_properties.hpp>
#include <geom/brep/primitives.hpp>
#include <geom/brep/tessellate.hpp>
#include <geom/brep/topology.hpp>
#include <geom/brep/validate.hpp>
#include <geom/curves/analytic_curve.hpp>
#include <geom/surfaces/analytic.hpp>

namespace {
using geom::Brep;
using geom::EdgeId;
using geom::Point3d;

Brep box_2x3x4() {
  return geom::make_box(Point3d{0.0, 0.0, 0.0}, Point3d{2.0, 3.0, 4.0}, 1e-9).value();
}

std::vector<EdgeId> all_edges(const Brep& brep) {
  const auto ids = brep.edge_ids();
  return std::vector<EdgeId>(ids.begin(), ids.end());
}

double volume_of(const Brep& brep) {
  return geom::detail::curved_solid_volume(brep, brep.solid_ids().front(), 1e-9).value();
}
}  // namespace

TEST_CASE("rounded box rounds all twelve edges with exact corner patches") {
  const Brep box = box_2x3x4();
  const std::vector<EdgeId> edges = all_edges(box);
  REQUIRE(edges.size() == 12u);

  const auto result = geom::fillet_edges(box, edges, geom::FilletOptions{0.25, 1e-9});
  REQUIRE(result.has_value());
  const Brep& rounded = result.value().brep;
  CHECK(result.value().filleted_edges == 12u);
  CHECK(rounded.solid_count() == 1u);
  CHECK(rounded.shell_count() == 1u);
  CHECK(rounded.vertex_count() == 24u);
  CHECK(rounded.edge_count() == 48u);
  CHECK(rounded.face_count() == 26u);
  CHECK(geom::validate_topology(rounded).ok());
  CHECK(volume_of(rounded) == doctest::Approx(23.53807140332).epsilon(1e-12));

  std::uint32_t planes = 0;
  std::uint32_t cylinders = 0;
  std::uint32_t spheres = 0;
  for (const geom::FaceId id : rounded.face_ids()) {
    const geom::SurfaceGeometry& surface =
        rounded.geometry().surfaces[rounded.face(id).surface];
    if (std::holds_alternative<geom::Plane3>(surface)) ++planes;
    if (std::holds_alternative<geom::Cylinder3>(surface)) ++cylinders;
    if (std::holds_alternative<geom::Sphere3>(surface)) ++spheres;
  }
  CHECK(planes == 6u);
  CHECK(cylinders == 12u);
  CHECK(spheres == 8u);
}

TEST_CASE("rounded box tessellates watertight") {
  const Brep box = box_2x3x4();
  const auto rounded = geom::fillet_edges(box, all_edges(box), geom::FilletOptions{0.25, 1e-9});
  REQUIRE(rounded.has_value());
  const auto mesh = geom::tessellate(rounded.value().brep, geom::TessellationOptions{});
  REQUIRE(mesh.has_value());
  CHECK(mesh.value().triangles.size() > 0u);
  std::map<std::pair<std::uint32_t, std::uint32_t>, std::uint32_t> uses;
  for (const auto& triangle : mesh.value().triangles) {
    for (int e = 0; e < 3; ++e) {
      std::uint32_t a = triangle[static_cast<std::size_t>(e)];
      std::uint32_t b = triangle[static_cast<std::size_t>((e + 1) % 3)];
      if (a > b) std::swap(a, b);
      ++uses[{a, b}];
    }
  }
  std::uint32_t boundary = 0;
  for (const auto& [edge, count] : uses) {
    if (count == 1) ++boundary;
  }
  CHECK(boundary == 0u);
  // No degenerate or duplicate triangles: a patch tessellator that fans from a
  // pole can emit collapsed rows if the pole column is handled twice.
  for (const auto& triangle : mesh.value().triangles) {
    const geom::Point3d& a = mesh.value().positions[triangle[0]];
    const geom::Point3d& b = mesh.value().positions[triangle[1]];
    const geom::Point3d& c = mesh.value().positions[triangle[2]];
    const geom::Vec3d n = geom::cross(b - a, c - a);
    CHECK(geom::dot(n, n) > 0.0);
  }
}

TEST_CASE("rounded box has exact bounds and corner patches") {
  const Brep box = box_2x3x4();
  const auto rounded = geom::fillet_edges(box, all_edges(box), geom::FilletOptions{0.25, 1e-9});
  REQUIRE(rounded.has_value());
  const Brep& solid = rounded.value().brep;
  double lo[3] = {1e300, 1e300, 1e300};
  double hi[3] = {-1e300, -1e300, -1e300};
  for (const geom::VertexId id : solid.vertex_ids()) {
    const geom::Point3d& p = solid.vertex(id).point;
    for (int axis = 0; axis < 3; ++axis) {
      lo[axis] = std::min(lo[axis], p[axis]);
      hi[axis] = std::max(hi[axis], p[axis]);
    }
  }
  CHECK(lo[0] == 0.0);
  CHECK(lo[1] == 0.0);
  CHECK(lo[2] == 0.0);
  CHECK(hi[0] == 2.0);
  CHECK(hi[1] == 3.0);
  CHECK(hi[2] == 4.0);
  bool saw_hi_plane = false;
  bool saw_corner_sphere = false;
  for (const geom::SurfaceGeometry& surface : solid.geometry().surfaces) {
    if (const geom::Plane3* plane = std::get_if<geom::Plane3>(&surface)) {
      if (plane->normal[0] == 1.0 && plane->offset == 2.0) saw_hi_plane = true;
    }
    if (const geom::Sphere3* sphere = std::get_if<geom::Sphere3>(&surface)) {
      if (sphere->center[0] == 1.75 && sphere->center[1] == 2.75 && sphere->center[2] == 3.75) {
        saw_corner_sphere = true;
      }
    }
  }
  CHECK(saw_hi_plane);
  CHECK(saw_corner_sphere);
}

TEST_CASE("rounded box mass properties match the Steiner volume") {
  const Brep box = box_2x3x4();
  const auto rounded = geom::fillet_edges(box, all_edges(box), geom::FilletOptions{0.25, 1e-9});
  REQUIRE(rounded.has_value());
  const auto mass = geom::mass_properties(rounded.value().brep);
  REQUIRE(mass.has_value());
  CHECK(mass.value().volume == doctest::Approx(23.53807140332).epsilon(1e-9));
  CHECK(mass.value().centroid[0] == doctest::Approx(1.0).epsilon(1e-9));
  CHECK(mass.value().centroid[1] == doctest::Approx(1.5).epsilon(1e-9));
  CHECK(mass.value().centroid[2] == doctest::Approx(2.0).epsilon(1e-9));
}
