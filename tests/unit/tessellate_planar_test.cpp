#include <doctest/doctest.h>

#include <array>
#include <cmath>
#include <cstdint>
#include <limits>
#include <map>
#include <utility>
#include <vector>

#include <geom/brep/primitives.hpp>
#include <geom/brep/tessellate.hpp>
#include <geom/core/error.hpp>
#include <geom/primitives/point.hpp>

namespace {

using geom::CoedgeSpec;
using geom::EdgeId;
using geom::FaceId;
using geom::Point3d;
using geom::VertexId;

// T4's SquareFace helper, copied into this translation unit (each test target
// compiles independently, so anonymous-namespace symbols cannot be shared).
struct SquareFace {
  geom::SolidEditor editor;
  geom::FaceId face;
  geom::LoopId loop;
  std::array<geom::VertexId, 4> vertices;
};

// Same topology as T4's helper; the surface is a parameter so the validation
// cases can build a square with a non-finite or zero plane normal.
SquareFace make_square_face_with_plane(const geom::Plane3& plane) {
  geom::SolidEditor ed;
  SquareFace s{std::move(ed), {}, {}, {}};
  const geom::Point3d pts[4] = {{0, 0, 0}, {1, 0, 0}, {1, 1, 0}, {0, 1, 0}};
  for (int i = 0; i < 4; ++i) s.vertices[i] = s.editor.add_vertex(pts[i], 0.0).value();
  const std::uint32_t surface = s.editor.add_surface(plane);
  const std::uint32_t line = s.editor.add_curve(geom::LineCurve3{{0, 0, 0}, {1, 0, 0}});
  std::array<geom::EdgeId, 4> edges{};
  for (int i = 0; i < 4; ++i) {
    edges[i] = s.editor.add_edge(s.vertices[i], s.vertices[(i + 1) % 4], line, 0.0, 1.0, 0.0)
                   .value();
  }
  s.face = s.editor.add_face(surface, false, 0.0).value();
  std::array<geom::CoedgeSpec, 4> specs{};
  for (int i = 0; i < 4; ++i) specs[i] = geom::CoedgeSpec{edges[i], false};
  s.loop = s.editor.add_loop(s.face, specs, true).value();
  auto shell = s.editor.add_shell(std::array<geom::FaceId, 1>{s.face}, false).value();
  s.editor.add_solid(std::array<geom::ShellId, 1>{shell}, 0.0).value();
  return s;
}

SquareFace make_square_face() {
  return make_square_face_with_plane(geom::Plane3{{0, 0, 1}, 0.0});
}

// Two unit squares sharing the edge (1,0,0)-(1,1,0), both wound CCW seen from
// +z. The shared edge is one Edge used forward by the left face and reversed by
// the right face; the tessellation must weld both uses to the same samples.
struct TwoSquares {
  geom::SolidEditor editor;
  std::array<FaceId, 2> faces;
  std::array<VertexId, 6> vertices;
};

TwoSquares make_two_squares() {
  geom::SolidEditor ed;
  TwoSquares s{std::move(ed), {}, {}};
  const geom::Point3d pts[6] = {{0, 0, 0}, {1, 0, 0}, {1, 1, 0},
                                {0, 1, 0}, {2, 0, 0}, {2, 1, 0}};
  for (int i = 0; i < 6; ++i) s.vertices[i] = s.editor.add_vertex(pts[i], 0.0).value();
  const std::uint32_t plane = s.editor.add_surface(geom::Plane3{{0, 0, 1}, 0.0});
  const int pairs[7][2] = {{0, 1}, {1, 2}, {2, 3}, {3, 0},
                           {1, 4}, {4, 5}, {5, 2}};
  std::array<EdgeId, 7> edges{};
  for (int i = 0; i < 7; ++i) {
    edges[i] = s.editor
                   .add_edge(s.vertices[pairs[i][0]], s.vertices[pairs[i][1]],
                             geom::kNoGeometry, 0.0, 1.0, 0.0)
                   .value();
  }
  const EdgeId shared_edge = edges[1];
  s.faces[0] = s.editor.add_face(plane, false, 0.0).value();
  s.faces[1] = s.editor.add_face(plane, false, 0.0).value();
  const std::array<CoedgeSpec, 4> left{CoedgeSpec{edges[0], false},
                                       CoedgeSpec{shared_edge, false},
                                       CoedgeSpec{edges[2], false},
                                       CoedgeSpec{edges[3], false}};
  const std::array<CoedgeSpec, 4> right{CoedgeSpec{edges[4], false},
                                        CoedgeSpec{edges[5], false},
                                        CoedgeSpec{edges[6], false},
                                        CoedgeSpec{shared_edge, true}};
  s.editor.add_loop(s.faces[0], left, true).value();
  s.editor.add_loop(s.faces[1], right, true).value();
  auto shell = s.editor.add_shell(std::array<geom::FaceId, 2>{s.faces[0], s.faces[1]}, false)
                   .value();
  s.editor.add_solid(std::array<geom::ShellId, 1>{shell}, 0.0).value();
  return s;
}

// Square annulus: outer 2x2 square, inner 1x1 hole, one planar face with an
// outer loop and a hole loop (both built CCW; the tessellator normalizes the
// hole ring to clockwise).
struct HoleFace {
  geom::SolidEditor editor;
  FaceId face;
  std::array<VertexId, 8> vertices;
};

HoleFace make_hole_face() {
  geom::SolidEditor ed;
  HoleFace h{std::move(ed), {}, {}};
  const geom::Point3d pts[8] = {{0, 0, 0}, {2, 0, 0}, {2, 2, 0}, {0, 2, 0},
                                {0.5, 0.5, 0}, {1.5, 0.5, 0}, {1.5, 1.5, 0}, {0.5, 1.5, 0}};
  for (int i = 0; i < 8; ++i) h.vertices[i] = h.editor.add_vertex(pts[i], 0.0).value();
  const std::uint32_t plane = h.editor.add_surface(geom::Plane3{{0, 0, 1}, 0.0});
  const int pairs[8][2] = {{0, 1}, {1, 2}, {2, 3}, {3, 0},
                           {4, 5}, {5, 6}, {6, 7}, {7, 4}};
  std::array<geom::EdgeId, 8> edges{};
  for (int i = 0; i < 8; ++i) {
    edges[i] = h.editor
                   .add_edge(h.vertices[pairs[i][0]], h.vertices[pairs[i][1]],
                             geom::kNoGeometry, 0.0, 1.0, 0.0)
                   .value();
  }
  h.face = h.editor.add_face(plane, false, 0.0).value();
  std::array<CoedgeSpec, 4> outer{};
  std::array<CoedgeSpec, 4> hole{};
  for (int i = 0; i < 4; ++i) {
    outer[i] = CoedgeSpec{edges[i], false};
    hole[i] = CoedgeSpec{edges[4 + i], false};
  }
  h.editor.add_loop(h.face, outer, true).value();
  h.editor.add_loop(h.face, hole, false).value();
  auto shell = h.editor.add_shell(std::array<geom::FaceId, 1>{h.face}, false).value();
  h.editor.add_solid(std::array<geom::ShellId, 1>{shell}, 0.0).value();
  return h;
}

// Index of the mesh position equal to `p`, or mesh.positions.size() when absent.
std::uint32_t find_position(const geom::TriMesh& mesh, const Point3d& p) {
  for (std::uint32_t i = 0; i < mesh.positions.size(); ++i) {
    if (mesh.positions[i][0] == p[0] && mesh.positions[i][1] == p[1] &&
        mesh.positions[i][2] == p[2]) {
      return i;
    }
  }
  return static_cast<std::uint32_t>(mesh.positions.size());
}

}  // namespace

TEST_CASE("single square face tessellates into two triangles") {
  SquareFace s = make_square_face();
  const auto b = s.editor.commit();
  REQUIRE(b.has_value());
  const auto result = geom::tessellate(b.value(), {});
  REQUIRE(result.has_value());
  const geom::TriMesh& mesh = result.value();
  CHECK(mesh.positions.size() == 4);
  CHECK(mesh.triangles.size() == 2);
  CHECK(mesh.triangle_face.size() == 2);
}

TEST_CASE("box mesh is closed and edge-manifold") {
  const auto b = geom::make_box({0, 0, 0}, {2, 3, 4}, 0.0);
  REQUIRE(b.has_value());
  const auto result = geom::tessellate(b.value(), {});
  REQUIRE(result.has_value());
  const geom::TriMesh& mesh = result.value();
  CHECK(mesh.positions.size() == 8);
  CHECK(mesh.triangles.size() == 12);
  std::map<std::pair<std::uint32_t, std::uint32_t>, int> edge_use;
  for (const auto& t : mesh.triangles) {
    for (int i = 0; i < 3; ++i) {
      auto a = t[i];
      auto c = t[(i + 1) % 3];
      if (a > c) std::swap(a, c);
      ++edge_use[{a, c}];
    }
  }
  CHECK(edge_use.size() == 18);
  for (const auto& [edge, count] : edge_use) CHECK(count == 2);
}

TEST_CASE("tessellation rejects non-positive deflection") {
  const auto b = geom::make_box({0, 0, 0}, {1, 1, 1}, 0.0);
  REQUIRE(b.has_value());
  CHECK_FALSE(geom::tessellate(b.value(), {0.0, 0.35}).has_value());
  CHECK_FALSE(geom::tessellate(b.value(), {-1.0, 0.35}).has_value());
}

TEST_CASE("adjacent faces share welded edge samples") {
  TwoSquares s = make_two_squares();
  const auto b = s.editor.commit();
  REQUIRE(b.has_value());
  const auto result = geom::tessellate(b.value(), {});
  REQUIRE(result.has_value());
  const geom::TriMesh& mesh = result.value();
  CHECK(mesh.positions.size() == 6);
  CHECK(mesh.triangles.size() == 4);
  // The shared edge endpoints appear exactly once, so both faces index the
  // same welded samples instead of duplicating them.
  CHECK(find_position(mesh, {1, 0, 0}) < mesh.positions.size());
  CHECK(find_position(mesh, {1, 1, 0}) < mesh.positions.size());
  std::size_t duplicated = 0;
  for (std::uint32_t i = 0; i < mesh.positions.size(); ++i) {
    for (std::uint32_t j = i + 1; j < mesh.positions.size(); ++j) {
      if (mesh.positions[i][0] == mesh.positions[j][0] &&
          mesh.positions[i][1] == mesh.positions[j][1] &&
          mesh.positions[i][2] == mesh.positions[j][2]) {
        ++duplicated;
      }
    }
  }
  CHECK(duplicated == 0);
}

TEST_CASE("square with a hole tessellates the annulus") {
  HoleFace h = make_hole_face();
  const auto b = h.editor.commit();
  REQUIRE(b.has_value());
  const auto result = geom::tessellate(b.value(), {});
  REQUIRE(result.has_value());
  const geom::TriMesh& mesh = result.value();
  CHECK(mesh.positions.size() == 8);
  REQUIRE(mesh.triangles.size() == 8);
  CHECK(mesh.triangle_face.size() == 8);
  for (std::uint32_t face : mesh.triangle_face) CHECK(face == h.face.index());
  double area = 0.0;
  for (const auto& t : mesh.triangles) {
    const geom::Vec3d ab = mesh.positions[t[1]] - mesh.positions[t[0]];
    const geom::Vec3d ac = mesh.positions[t[2]] - mesh.positions[t[0]];
    const geom::Vec3d n = geom::cross(ab, ac);
    area += 0.5 * std::hypot(n[0], n[1], n[2]);
  }
  CHECK(area == doctest::Approx(3.0).epsilon(1e-12));
  // This face is open along the hole: the hole ring is a mesh boundary, so every
  // hole edge belongs to exactly one triangle. Only the bridges and the diagonals
  // inside the annulus are shared by two triangles.
  const Point3d hole_pts[4] = {{0.5, 0.5, 0}, {1.5, 0.5, 0}, {1.5, 1.5, 0}, {0.5, 1.5, 0}};
  std::uint32_t hole_index[4] = {0, 0, 0, 0};
  for (int i = 0; i < 4; ++i) {
    hole_index[i] = find_position(mesh, hole_pts[i]);
    REQUIRE(hole_index[i] < mesh.positions.size());
  }
  std::map<std::pair<std::uint32_t, std::uint32_t>, int> edge_use;
  for (const auto& t : mesh.triangles) {
    for (int i = 0; i < 3; ++i) {
      auto a = t[i];
      auto c = t[(i + 1) % 3];
      if (a > c) std::swap(a, c);
      ++edge_use[{a, c}];
    }
  }
  for (int i = 0; i < 4; ++i) {
    auto a = hole_index[i];
    auto c = hole_index[(i + 1) % 4];
    if (a > c) std::swap(a, c);
    CHECK(edge_use[{a, c}] == 1);
  }
}

TEST_CASE("triangles follow the face normal") {
  const auto b = geom::make_box({0, 0, 0}, {1, 1, 1}, 0.0);
  REQUIRE(b.has_value());
  const auto result = geom::tessellate(b.value(), {});
  REQUIRE(result.has_value());
  const geom::TriMesh& mesh = result.value();
  const Point3d center{0.5, 0.5, 0.5};
  for (const auto& t : mesh.triangles) {
    const geom::Vec3d ab = mesh.positions[t[1]] - mesh.positions[t[0]];
    const geom::Vec3d ac = mesh.positions[t[2]] - mesh.positions[t[0]];
    const geom::Vec3d n = geom::cross(ab, ac);
    REQUIRE(std::hypot(n[0], n[1], n[2]) > 0.0);
    const Point3d centroid{
        (mesh.positions[t[0]][0] + mesh.positions[t[1]][0] + mesh.positions[t[2]][0]) / 3.0,
        (mesh.positions[t[0]][1] + mesh.positions[t[1]][1] + mesh.positions[t[2]][1]) / 3.0,
        (mesh.positions[t[0]][2] + mesh.positions[t[1]][2] + mesh.positions[t[2]][2]) / 3.0};
    const geom::Vec3d outward = centroid - center;
    CHECK(geom::dot(n, outward) > 0.0);
  }
}

TEST_CASE("tessellation rejects invalid tolerance options") {
  const auto b = geom::make_box({0, 0, 0}, {1, 1, 1}, 0.0);
  REQUIRE(b.has_value());
  const auto check_rejected = [&](double deflection, double angle) {
    const auto r = geom::tessellate(b.value(), {deflection, angle});
    REQUIRE_FALSE(r.has_value());
    CHECK(r.error().code == geom::Errc::tolerance_violation);
  };
  check_rejected(0.0, 0.35);
  check_rejected(-1.0, 0.35);
  check_rejected(std::numeric_limits<double>::quiet_NaN(), 0.35);
  check_rejected(1e-3, 0.0);
  check_rejected(1e-3, -0.5);
  check_rejected(1e-3, std::numeric_limits<double>::quiet_NaN());
  check_rejected(1e-3, std::numeric_limits<double>::infinity());
}

TEST_CASE("tessellation rejects non-finite and degenerate planes") {
  {
    SquareFace s = make_square_face_with_plane(
        geom::Plane3{{std::numeric_limits<double>::quiet_NaN(), 0.0, 1.0}, 0.0});
    const auto b = s.editor.commit();
    REQUIRE(b.has_value());
    const auto r = geom::tessellate(b.value(), {});
    REQUIRE_FALSE(r.has_value());
    CHECK(r.error().code == geom::Errc::non_finite);
  }
  {
    SquareFace s = make_square_face_with_plane(geom::Plane3{{0, 0, 0}, 0.0});
    const auto b = s.editor.commit();
    REQUIRE(b.has_value());
    const auto r = geom::tessellate(b.value(), {});
    REQUIRE_FALSE(r.has_value());
    CHECK(r.error().code == geom::Errc::degenerate);
  }
}

TEST_CASE("non-planar faces are not tessellated yet") {
  const auto cylinder = geom::make_cylinder({0, 0, 0}, {0, 0, 1}, 1.0, 2.0, 0.0);
  REQUIRE(cylinder.has_value());
  const auto r = geom::tessellate(cylinder.value(), {});
  REQUIRE_FALSE(r.has_value());
  CHECK(r.error().code == geom::Errc::unsupported);
}
