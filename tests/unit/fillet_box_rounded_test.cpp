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
