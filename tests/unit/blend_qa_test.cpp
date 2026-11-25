#include <doctest/doctest.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <map>
#include <span>
#include <utility>
#include <vector>

#include <geom/brep/chamfer.hpp>
#include <geom/brep/mass_properties.hpp>
#include <geom/brep/primitives.hpp>
#include <geom/brep/shell.hpp>
#include <geom/brep/tessellate.hpp>
#include <geom/brep/validate.hpp>
#include <geom/core/error.hpp>

using namespace geom;

namespace {

Brep box_solid(Point3d lo, Point3d hi) {
  const Result<Brep> box = make_box(lo, hi, 1e-9);
  REQUIRE(box.has_value());
  return box.value();
}

// Returns the number of triangle edges used exactly once.
std::size_t boundary_edges(const TriMesh& mesh) {
  std::map<std::pair<std::uint32_t, std::uint32_t>, int> uses;
  for (const auto& tri : mesh.triangles) {
    for (int i = 0; i < 3; ++i) {
      std::uint32_t a = tri[static_cast<std::size_t>(i)];
      std::uint32_t b = tri[static_cast<std::size_t>((i + 1) % 3)];
      if (a > b) std::swap(a, b);
      ++uses[{a, b}];
    }
  }
  std::size_t count = 0;
  for (const auto& [edge, n] : uses) {
    (void)edge;
    if (n == 1) ++count;
  }
  return count;
}

}  // namespace

TEST_CASE("chamfer and shell outputs tessellate watertight") {
  const Brep box = box_solid(Point3d{0.0, 0.0, 0.0}, Point3d{2.0, 3.0, 4.0});
  ChamferOptions chamfer;
  chamfer.distance = 0.25;
  std::vector<EdgeId> edges = box.edge_ids();
  const Result<ChamferResult> beveled = chamfer_edges(box, edges, chamfer);
  REQUIRE(beveled.has_value());
  const Result<TriMesh> bevel_mesh = tessellate(beveled.value().brep, TessellationOptions{});
  REQUIRE(bevel_mesh.has_value());
  CHECK_FALSE(bevel_mesh.value().triangles.empty());
  CHECK(boundary_edges(bevel_mesh.value()) == 0u);

  ShellOptions shell;
  shell.thickness = 0.25;
  const std::array<FaceId, 0> open{};
  const Result<ShellResult> hollow = shell_solid(box, open, shell);
  REQUIRE(hollow.has_value());
  const Result<TriMesh> shell_mesh = tessellate(hollow.value().brep, TessellationOptions{});
  REQUIRE(shell_mesh.has_value());
  CHECK(boundary_edges(shell_mesh.value()) == 0u);
}

TEST_CASE("blend outputs are deterministic and leave inputs unchanged") {
  const Brep box = box_solid(Point3d{0.0, 0.0, 0.0}, Point3d{2.0, 2.0, 2.0});
  ChamferOptions chamfer;
  chamfer.distance = 0.25;
  const std::vector<EdgeId> edges = box.edge_ids();
  const Result<ChamferResult> first = chamfer_edges(box, edges, chamfer);
  const Result<ChamferResult> second = chamfer_edges(box, edges, chamfer);
  REQUIRE(first.has_value());
  REQUIRE(second.has_value());
  CHECK(first.value().brep.vertex_count() == second.value().brep.vertex_count());
  CHECK(first.value().brep.edge_count() == second.value().brep.edge_count());
  CHECK(first.value().brep.face_count() == second.value().brep.face_count());
  const Result<MassProperties> a = mass_properties(first.value().brep);
  const Result<MassProperties> b = mass_properties(second.value().brep);
  REQUIRE(a.has_value());
  REQUIRE(b.has_value());
  CHECK(a.value().volume == b.value().volume);

  const Result<TriMesh> mesh_a = tessellate(first.value().brep, TessellationOptions{});
  const Result<TriMesh> mesh_b = tessellate(second.value().brep, TessellationOptions{});
  REQUIRE(mesh_a.has_value());
  REQUIRE(mesh_b.has_value());
  CHECK(mesh_a.value().positions.size() == mesh_b.value().positions.size());
  CHECK(mesh_a.value().triangles == mesh_b.value().triangles);

  CHECK(box.face_count() == 6u);
  const Result<MassProperties> original = mass_properties(box);
  REQUIRE(original.has_value());
  CHECK(std::abs(original.value().volume - 8.0) <= 1e-12);
}
