#include <doctest/doctest.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <map>
#include <set>
#include <utility>
#include <variant>
#include <vector>

#include <geom/brep/primitives.hpp>
#include <geom/brep/tessellate.hpp>
#include <geom/primitives/point.hpp>

namespace {

using geom::Point3d;

// Reversed curved faces: the analytic primitives create faces with
// reversed=false, and SolidEditor has no operation that flips an existing
// face, so a curved reversed face cannot be authored here without editing the
// Brep internals. The planar "reversed closed-circle disk" case covers the
// reversed-face triangle swap; curved reversal waits for an API that can
// construct it. The same applies to the `angle_tolerance` option, which is
// validated only (see tessellate.hpp).

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

geom::Vec3d triangle_normal(const geom::TriMesh& mesh, std::size_t index) {
  const auto& t = mesh.triangles[index];
  return geom::cross(mesh.positions[t[1]] - mesh.positions[t[0]],
                     mesh.positions[t[2]] - mesh.positions[t[0]]);
}

geom::Vec3d triangle_centroid(const geom::TriMesh& mesh, std::size_t index) {
  const auto& t = mesh.triangles[index];
  const geom::Point3d& a = mesh.positions[t[0]];
  const geom::Point3d& b = mesh.positions[t[1]];
  const geom::Point3d& c = mesh.positions[t[2]];
  return geom::Vec3d{(a[0] + b[0] + c[0]) / 3.0, (a[1] + b[1] + c[1]) / 3.0,
                     (a[2] + b[2] + c[2]) / 3.0};
}

// Outward reference direction: each triangle normal must agree with the
// surface's outward direction at that triangle, up to a scale-relative
// rounding tolerance (never a sign flip).
template <typename Outward>
void check_triangles_point_outward(const geom::TriMesh& mesh, Outward outward) {
  REQUIRE(!mesh.triangles.empty());
  for (std::size_t i = 0; i < mesh.triangles.size(); ++i) {
    const geom::Vec3d n = triangle_normal(mesh, i);
    const geom::Vec3d ref = outward(i);
    const double scale = geom::norm2(n) * geom::norm2(ref);
    REQUIRE(scale > 0.0);
    CHECK(geom::dot(n, ref) > -1e-12 * scale);
  }
}

// Face roles for a cylinder: 1 = bottom cap, 2 = lateral, 3 = top cap.
std::vector<int> cylinder_face_roles(const geom::Brep& brep, const geom::Vec3d& axis) {
  std::vector<int> roles(brep.faces().size(), 0);
  for (const geom::FaceId id : brep.face_ids()) {
    const geom::Face& face = brep.face(id);
    if (face.surface == geom::kNoGeometry) continue;
    const geom::SurfaceGeometry& geometry = brep.geometry().surfaces[face.surface];
    if (std::get_if<geom::Cylinder3>(&geometry) != nullptr) {
      roles[id.index()] = 2;
    } else if (const auto* plane = std::get_if<geom::Plane3>(&geometry)) {
      roles[id.index()] = geom::dot(plane->normal, axis) > 0.0 ? 3 : 1;
    }
  }
  return roles;
}

// Face roles for a cone: 1 = base cap, 2 = lateral.
std::vector<int> cone_face_roles(const geom::Brep& brep) {
  std::vector<int> roles(brep.faces().size(), 0);
  for (const geom::FaceId id : brep.face_ids()) {
    const geom::Face& face = brep.face(id);
    if (face.surface == geom::kNoGeometry) continue;
    const geom::SurfaceGeometry& geometry = brep.geometry().surfaces[face.surface];
    roles[id.index()] = std::get_if<geom::Cone3>(&geometry) != nullptr ? 2 : 1;
  }
  return roles;
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
    CHECK(2.0 - std::hypot(mx, my) <= 0.01 + 1e-12);  // 鐭㈤珮
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

TEST_CASE("cylinder triangles point outward and carry their face ids") {
  const geom::Point3d base{0, 0, 0};
  const geom::Vec3d axis{0, 0, 1};
  const auto b = geom::make_cylinder(base, axis, 2.0, 5.0, 0.0);
  REQUIRE(b.has_value());
  const geom::Brep& brep = b.value();
  const auto result = geom::tessellate(brep, {0.05, 0.35});
  REQUIRE(result.has_value());
  const geom::TriMesh& mesh = result.value();
  check_closed_manifold(mesh);
  CHECK(euler_characteristic(mesh) == 2);
  REQUIRE(mesh.triangle_face.size() == mesh.triangles.size());
  const std::vector<int> roles = cylinder_face_roles(brep, axis);
  std::size_t counts[4] = {0, 0, 0, 0};
  for (std::size_t i = 0; i < mesh.triangles.size(); ++i) {
    const int role = roles[mesh.triangle_face[i]];
    REQUIRE(role != 0);
    ++counts[role];
    const geom::Vec3d c = triangle_centroid(mesh, i);
    // The id -> role mapping must match the geometry: caps are the z = 0 and
    // z = 5 rings, the lateral triangles sit strictly between them.
    if (role == 1) CHECK(std::abs(c[2]) < 1e-12);
    if (role == 3) CHECK(std::abs(c[2] - 5.0) < 1e-12);
    if (role == 2) {
      CHECK(c[2] > 0.0);
      CHECK(c[2] < 5.0);
    }
    const geom::Vec3d n = triangle_normal(mesh, i);
    geom::Vec3d outward{};
    if (role == 1) {
      outward = axis * -1.0;
    } else if (role == 3) {
      outward = axis;
    } else {
      const geom::Point3d point{c[0], c[1], c[2]};
      outward = (point - base) - axis * geom::dot(point - base, axis);
    }
    const double scale = geom::norm2(n) * geom::norm2(outward);
    REQUIRE(scale > 0.0);
    CHECK(geom::dot(n, outward) > -1e-12 * scale);
  }
  CHECK(counts[1] > 0);
  CHECK(counts[2] > 0);
  CHECK(counts[3] > 0);
}

TEST_CASE("cone triangles point outward along the slant") {
  const geom::Point3d origin{0, 0, 0};
  const geom::Vec3d axis{0, 0, 1};
  const auto b = geom::make_cone(origin, axis, 2.0, 5.0, 0.0);
  REQUIRE(b.has_value());
  const geom::Brep& brep = b.value();
  const auto result = geom::tessellate(brep, {0.05, 0.35});
  REQUIRE(result.has_value());
  const geom::TriMesh& mesh = result.value();
  check_closed_manifold(mesh);
  CHECK(euler_characteristic(mesh) == 2);
  const double half_angle = std::atan2(2.0, 5.0);
  const double cos_a = std::cos(half_angle);
  const double sin_a = std::sin(half_angle);
  const std::vector<int> roles = cone_face_roles(brep);
  std::size_t lateral = 0;
  std::size_t cap = 0;
  for (std::size_t i = 0; i < mesh.triangles.size(); ++i) {
    const int role = roles[mesh.triangle_face[i]];
    REQUIRE(role != 0);
    const geom::Vec3d n = triangle_normal(mesh, i);
    const geom::Vec3d c = triangle_centroid(mesh, i);
    geom::Vec3d outward{};
    if (role == 1) {
      ++cap;
      outward = axis * -1.0;
    } else {
      ++lateral;
      const geom::Point3d point{c[0], c[1], c[2]};
      const geom::Vec3d radial = (point - origin) - axis * geom::dot(point - origin, axis);
      outward = radial * cos_a + axis * sin_a;
    }
    const double scale = geom::norm2(n) * geom::norm2(outward);
    REQUIRE(scale > 0.0);
    CHECK(geom::dot(n, outward) > -1e-12 * scale);
  }
  CHECK(cap > 0);
  CHECK(lateral > 0);
}

TEST_CASE("sphere samples lie on the sphere and point outward") {
  const geom::Point3d center{0, 0, 0};
  const auto b = geom::make_sphere(center, 2.0, 0.0);
  REQUIRE(b.has_value());
  const auto result = geom::tessellate(b.value(), {0.05, 0.35});
  REQUIRE(result.has_value());
  const geom::TriMesh& mesh = result.value();
  check_closed_manifold(mesh);
  CHECK(euler_characteristic(mesh) == 2);
  for (const geom::Point3d& p : mesh.positions) {
    const double radius = geom::norm2(p - center);
    // Required bound: every sample lies within the deflection of the sphere.
    CHECK(std::abs(radius - 2.0) <= 0.05);
    // The samples are exact surface evaluations, so the residual is far below
    // the chord bound; this catches drift to a wrong radius.
    CHECK(std::abs(radius - 2.0) <= 1e-12);
  }
  check_triangles_point_outward(mesh, [&](std::size_t i) {
    const geom::Vec3d c = triangle_centroid(mesh, i);
    return geom::Point3d{c[0], c[1], c[2]} - center;
  });
}

TEST_CASE("torus samples satisfy the implicit equation and point outward") {
  const geom::Point3d center{0, 0, 0};
  const geom::Vec3d axis{0, 0, 1};
  const double major = 3.0;
  const double minor = 1.0;
  const auto b = geom::make_torus(center, axis, major, minor, 0.0);
  REQUIRE(b.has_value());
  const auto result = geom::tessellate(b.value(), {0.05, 0.35});
  REQUIRE(result.has_value());
  const geom::TriMesh& mesh = result.value();
  check_closed_manifold(mesh);
  CHECK(euler_characteristic(mesh) == 0);
  for (const geom::Point3d& p : mesh.positions) {
    const geom::Vec3d d = p - center;
    const double xy = d[0] * d[0] + d[1] * d[1];
    const double rr = xy + d[2] * d[2];
    const double lhs = rr + major * major - minor * minor;
    const double f = lhs * lhs - 4.0 * major * major * xy;
    const double scale = rr * rr + 4.0 * major * major * xy;
    // f = (|p|^2 + R^2 - r^2)^2 - 4 R^2 (x^2 + y^2) is quartic; bound it
    // relative to its own magnitude. Samples are exact evaluations, so the
    // residual is near machine precision.
    CHECK(std::abs(f) <= 1e-9 * (1.0 + scale));
  }
  check_triangles_point_outward(mesh, [&](std::size_t i) {
    const geom::Vec3d c = triangle_centroid(mesh, i);
    const geom::Point3d point{c[0], c[1], c[2]};
    const geom::Vec3d radial = (point - center) - axis * geom::dot(point - center, axis);
    const double r = geom::norm2(radial);
    REQUIRE(r > 0.0);
    const geom::Point3d tube_center = center + radial * (major / r);
    return point - tube_center;
  });
}
