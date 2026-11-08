#include <doctest/doctest.h>

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <map>
#include <set>
#include <utility>
#include <vector>

#include <geom/brep/primitives.hpp>
#include <geom/brep/tessellate.hpp>
#include <geom/primitives/point.hpp>

namespace {

using geom::Point3d;

// Every undirected mesh edge must be used exactly twice by a closed solid.
void check_closed_manifold(const geom::TriMesh& mesh) {
  std::map<std::pair<std::uint32_t, std::uint32_t>, int> use;
  for (const auto& t : mesh.triangles) {
    const geom::Vec3d a = mesh.positions[t[1]] - mesh.positions[t[0]];
    const geom::Vec3d b = mesh.positions[t[2]] - mesh.positions[t[0]];
    const geom::Vec3d normal = geom::cross(a, b);
    CHECK(geom::dot(normal, normal) > 0.0);
    for (int i = 0; i < 3; ++i) {
      auto p = t[i];
      auto q = t[(i + 1) % 3];
      if (p > q) std::swap(p, q);
      ++use[{p, q}];
    }
  }
  REQUIRE(!use.empty());
  for (const auto& [edge, count] : use) CHECK(count == 2);
}

std::int64_t euler_characteristic(const geom::TriMesh& mesh) {
  std::set<std::pair<std::uint32_t, std::uint32_t>> edges;
  for (const auto& t : mesh.triangles) {
    for (int i = 0; i < 3; ++i) {
      auto a = t[i];
      auto b = t[(i + 1) % 3];
      if (a > b) std::swap(a, b);
      edges.insert({a, b});
    }
  }
  return static_cast<std::int64_t>(mesh.positions.size()) -
         static_cast<std::int64_t>(edges.size()) +
         static_cast<std::int64_t>(mesh.triangles.size());
}

std::size_t count_position(const geom::TriMesh& mesh, const Point3d& point) {
  std::size_t count = 0;
  for (const Point3d& p : mesh.positions) {
    if (p[0] == point[0] && p[1] == point[1] && p[2] == point[2]) ++count;
  }
  return count;
}

}  // namespace

TEST_CASE("cylinder mesh is watertight across the seam and caps") {
  const auto b = geom::make_cylinder({0, 0, 0}, {0, 0, 1}, 2.0, 5.0, 0.0);
  REQUIRE(b.has_value());
  const auto result = geom::tessellate(b.value(), {0.02, 0.35});
  REQUIRE(result.has_value());
  const geom::TriMesh& mesh = result.value();
  // Closed 2-manifold: every undirected edge is used exactly twice, and no
  // triangle has zero area.
  std::map<std::pair<std::uint32_t, std::uint32_t>, int> use;
  for (const auto& t : mesh.triangles) {
    const geom::Point3d a = mesh.positions[t[0]];
    const geom::Point3d c = mesh.positions[t[1]];
    const geom::Point3d d = mesh.positions[t[2]];
    const geom::Vec3d n = geom::cross(c - a, d - a);
    CHECK(geom::dot(n, n) > 0.0);
    for (int i = 0; i < 3; ++i) {
      auto p = t[i];
      auto q = t[(i + 1) % 3];
      if (p > q) std::swap(p, q);
      ++use[{p, q}];
    }
  }
  for (const auto& [e, count] : use) CHECK(count == 2);
}

TEST_CASE("sphere and torus meshes have Euler characteristic 2 and 0") {
  const auto chi = [](const geom::TriMesh& m) {
    std::set<std::pair<std::uint32_t, std::uint32_t>> edges;
    for (const auto& t : m.triangles) {
      for (int i = 0; i < 3; ++i) {
        auto a = t[i];
        auto b = t[(i + 1) % 3];
        if (a > b) std::swap(a, b);
        edges.insert({a, b});
      }
    }
    return static_cast<std::int64_t>(m.positions.size()) -
           static_cast<std::int64_t>(edges.size()) +
           static_cast<std::int64_t>(m.triangles.size());
  };
  const auto s = geom::make_sphere({0, 0, 0}, 2.0, 0.0);
  REQUIRE(s.has_value());
  const auto sm = geom::tessellate(s.value(), {0.02, 0.35});
  REQUIRE(sm.has_value());
  const geom::TriMesh& sphere_mesh = sm.value();
  CHECK(sphere_mesh.positions.size() > 100);
  CHECK(chi(sphere_mesh) == 2);
  const auto t = geom::make_torus({0, 0, 0}, {0, 0, 1}, 3.0, 1.0, 0.0);
  REQUIRE(t.has_value());
  const auto tm = geom::tessellate(t.value(), {0.02, 0.35});
  REQUIRE(tm.has_value());
  CHECK(chi(tm.value()) == 0);
}

TEST_CASE("chord height stays within deflection on a circle") {
  const auto b = geom::make_cylinder({0, 0, 0}, {0, 0, 1}, 2.0, 1.0, 0.0);
  REQUIRE(b.has_value());
  const auto result = geom::tessellate(b.value(), {0.01, 0.35});
  REQUIRE(result.has_value());
  const geom::TriMesh& mesh = result.value();
  std::vector<geom::Point3d> ring;
  for (const geom::Point3d& p : mesh.positions) {
    if (std::abs(p[2]) < 1e-12) ring.push_back(p);
  }
  REQUIRE(ring.size() > 8);
  std::sort(ring.begin(), ring.end(), [](const geom::Point3d& a, const geom::Point3d& c) {
    return std::atan2(a[1], a[0]) < std::atan2(c[1], c[0]);
  });
  for (std::size_t i = 0; i < ring.size(); ++i) {
    const geom::Point3d& a = ring[i];
    const geom::Point3d& c = ring[(i + 1) % ring.size()];
    const double mx = 0.5 * (a[0] + c[0]);
    const double my = 0.5 * (a[1] + c[1]);
    CHECK(2.0 - std::hypot(mx, my) <= 0.01 + 1e-12);  // 矢高
  }
}

TEST_CASE("cone mesh is watertight with a single apex") {
  const auto b = geom::make_cone({0, 0, 0}, {0, 0, 1}, 2.0, 5.0, 0.0);
  REQUIRE(b.has_value());
  const auto result = geom::tessellate(b.value(), {0.02, 0.35});
  REQUIRE(result.has_value());
  const geom::TriMesh& mesh = result.value();
  check_closed_manifold(mesh);
  CHECK(euler_characteristic(mesh) == 2);
  // The apex is one position, not one copy per base sample.
  CHECK(count_position(mesh, {0, 0, 5}) == 1);
  for (const auto& t : mesh.triangles) {
    const geom::Vec3d a = mesh.positions[t[1]] - mesh.positions[t[0]];
    const geom::Vec3d c = mesh.positions[t[2]] - mesh.positions[t[0]];
    const geom::Vec3d n = geom::cross(a, c);
    CHECK(geom::dot(n, n) > 0.0);
  }
}

TEST_CASE("cylinder caps weld to the lateral ring") {
  const auto b = geom::make_cylinder({0, 0, 0}, {0, 0, 1}, 2.0, 5.0, 0.0);
  REQUIRE(b.has_value());
  const auto result = geom::tessellate(b.value(), {0.05, 0.35});
  REQUIRE(result.has_value());
  const geom::TriMesh& mesh = result.value();
  std::size_t bottom = 0;
  std::size_t top = 0;
  for (const geom::Point3d& p : mesh.positions) {
    if (std::abs(p[2]) < 1e-12) ++bottom;
    if (std::abs(p[2] - 5.0) < 1e-12) ++top;
  }
  // Each cap boundary and the lateral ring must be the same welded samples.
  CHECK(bottom > 8);
  CHECK(bottom == top);
  CHECK(mesh.positions.size() == 2 * bottom);
}

TEST_CASE("sphere poles are single welded vertices") {
  const auto b = geom::make_sphere({0, 0, 0}, 2.0, 0.0);
  REQUIRE(b.has_value());
  const auto result = geom::tessellate(b.value(), {0.05, 0.35});
  REQUIRE(result.has_value());
  const geom::TriMesh& mesh = result.value();
  check_closed_manifold(mesh);
  CHECK(count_position(mesh, {0, 0, -2}) == 1);
  CHECK(count_position(mesh, {0, 0, 2}) == 1);
  CHECK(euler_characteristic(mesh) == 2);
}

TEST_CASE("curved tessellation is deterministic") {
  const auto b = geom::make_torus({0, 0, 0}, {0, 0, 1}, 3.0, 1.0, 0.0);
  REQUIRE(b.has_value());
  const auto first = geom::tessellate(b.value(), {0.05, 0.35});
  const auto second = geom::tessellate(b.value(), {0.05, 0.35});
  REQUIRE(first.has_value());
  REQUIRE(second.has_value());
  const geom::TriMesh& a = first.value();
  const geom::TriMesh& b_mesh = second.value();
  check_closed_manifold(a);
  CHECK(euler_characteristic(a) == 0);
  REQUIRE(a.positions.size() == b_mesh.positions.size());
  for (std::size_t i = 0; i < a.positions.size(); ++i) {
    CHECK(a.positions[i][0] == b_mesh.positions[i][0]);
    CHECK(a.positions[i][1] == b_mesh.positions[i][1]);
    CHECK(a.positions[i][2] == b_mesh.positions[i][2]);
  }
  CHECK(a.triangles == b_mesh.triangles);
  CHECK(a.triangle_face == b_mesh.triangle_face);
}
