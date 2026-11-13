#include <doctest/doctest.h>

#include <geom/brep/boolean_solid_curved.hpp>
#include <geom/brep/detail/curved_volume.hpp>
#include <geom/brep/primitives.hpp>
#include <geom/brep/tessellate.hpp>
#include <geom/brep/validate.hpp>

#include <cmath>
#include <map>
#include <numbers>
#include <utility>

namespace {

geom::Brep cut_cylinder() {
  // Same fixture as the trimmed-domain tests: a side window that crosses the seam.
  const geom::Brep cylinder =
      geom::make_cylinder(geom::Point3d{0.0, 0.0, 0.0}, geom::Vec3d{0.0, 0.0, 1.0}, 1.0, 2.0,
                          1e-9)
          .value();
  const geom::Brep box =
      geom::make_box(geom::Point3d{0.5, -0.5, 0.5}, geom::Point3d{2.0, 0.5, 1.5}, 1e-9).value();
  return geom::boolean_solid_curved(cylinder, cylinder.solid_ids().front(), box,
                                    box.solid_ids().front(), geom::SolidBooleanOp::difference,
                                    geom::CurvedBooleanOptions{1e-9})
      .value()
      .brep;
}

void check_watertight(const geom::TriMesh& mesh) {
  std::map<std::pair<std::uint32_t, std::uint32_t>, int> uses;
  for (const auto& triangle : mesh.triangles) {
    for (int k = 0; k < 3; ++k) {
      const std::uint32_t a = triangle[static_cast<std::size_t>(k)];
      const std::uint32_t b = triangle[static_cast<std::size_t>((k + 1) % 3)];
      uses[{std::min(a, b), std::max(a, b)}] += 1;
    }
  }
  for (const auto& [edge, count] : uses) CHECK(count == 2);
}

double mesh_volume(const geom::TriMesh& mesh) {
  double volume = 0.0;
  for (const auto& triangle : mesh.triangles) {
    const geom::Point3d& p = mesh.positions[triangle[0]];
    const geom::Point3d& q = mesh.positions[triangle[1]];
    const geom::Point3d& r = mesh.positions[triangle[2]];
    const geom::Vec3d a{p[0], p[1], p[2]};
    const geom::Vec3d b{q[0], q[1], q[2]};
    const geom::Vec3d c{r[0], r[1], r[2]};
    volume += geom::dot(a, geom::cross(b, c)) / 6.0;
  }
  return volume;
}

}  // namespace

TEST_CASE("trimmed tessellation: windowed cylinder is watertight and keeps its faces") {
  const geom::Brep brep = cut_cylinder();
  REQUIRE(geom::validate_topology(brep).ok());
  const auto mesh = geom::tessellate(brep, geom::TessellationOptions{1e-3, 0.35});
  REQUIRE(mesh.has_value());
  check_watertight(mesh.value());
  CHECK(mesh.value().triangle_face.size() == mesh.value().triangles.size());
  for (const std::uint32_t face : mesh.value().triangle_face) {
    CHECK(face < brep.face_ids().size());
  }
}

TEST_CASE("trimmed tessellation: every vertex lies on its analytic surface") {
  const geom::Brep brep = cut_cylinder();
  const auto mesh = geom::tessellate(brep, geom::TessellationOptions{2e-3, 0.35});
  REQUIRE(mesh.has_value());
  for (std::size_t i = 0; i < mesh.value().positions.size(); ++i) {
    const geom::Point3d& p = mesh.value().positions[i];
    const double radial = std::hypot(p[0], p[1]);
    const bool on_side = std::abs(radial - 1.0) <= 2e-3;
    const bool on_cap = std::abs(p[2]) <= 2e-3 || std::abs(p[2] - 2.0) <= 2e-3;
    const bool on_window = std::abs(p[0] - 0.5) <= 2e-3 || std::abs(p[0] + 0.5) <= 2e-3 || std::abs(p[1] - 0.5) <= 2e-3 || std::abs(p[1] + 0.5) <= 2e-3 ||
                           std::abs(p[2] - 0.5) <= 2e-3 || std::abs(p[2] - 1.5) <= 2e-3;
    CHECK((on_side || on_cap || on_window));
  }
}

TEST_CASE("trimmed tessellation: canonical primitives are unchanged") {
  const geom::Brep cylinder =
      geom::make_cylinder(geom::Point3d{0.0, 0.0, 0.0}, geom::Vec3d{0.0, 0.0, 1.0}, 1.0, 2.0,
                          1e-9)
          .value();
  // A 1e-2 deflection leaves a 23-gon cap with a 1.24% chord deficit, so the
  // brief's 1e-3 volume tolerance cannot pass even on the unchanged canonical
  // path; 1e-4 keeps the asserted property and tolerance.
  const auto mesh = geom::tessellate(cylinder, geom::TessellationOptions{1e-4, 0.35});
  REQUIRE(mesh.has_value());
  check_watertight(mesh.value());
  CHECK(mesh_volume(mesh.value()) == doctest::Approx(2.0 * std::numbers::pi_v<double>).epsilon(1e-3));
}
