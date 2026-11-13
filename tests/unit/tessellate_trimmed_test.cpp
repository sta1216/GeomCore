#include <doctest/doctest.h>

#include <geom/brep/boolean_solid_curved.hpp>
#include <geom/brep/detail/curved_volume.hpp>
#include <geom/brep/face_uv.hpp>
#include <geom/brep/primitives.hpp>
#include <geom/brep/solid_editor.hpp>
#include <geom/brep/split_solid.hpp>
#include <geom/brep/tessellate.hpp>
#include <geom/brep/validate.hpp>

#include <cmath>
#include <cstdint>
#include <cstring>
#include <map>
#include <numbers>
#include <string>
#include <utility>
#include <variant>

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

std::uint64_t hash_mesh(const geom::TriMesh& mesh) {
  std::uint64_t hash = 1469598103934665603ull;
  const auto mix = [&](const void* data, std::size_t size) {
    const auto* bytes = static_cast<const unsigned char*>(data);
    for (std::size_t i = 0; i < size; ++i) {
      hash ^= bytes[i];
      hash *= 1099511628211ull;
    }
  };
  for (const geom::Point3d& point : mesh.positions) {
    const double values[3] = {point[0], point[1], point[2]};
    mix(values, sizeof(values));
  }
  for (const auto& triangle : mesh.triangles) mix(triangle.data(), sizeof(triangle));
  for (const std::uint32_t face : mesh.triangle_face) mix(&face, sizeof(face));
  return hash;
}

geom::FaceId lateral_face(const geom::Brep& brep) {
  for (const geom::FaceId id : brep.face_ids()) {
    if (std::holds_alternative<geom::Cylinder3>(
            brep.geometry().surfaces[brep.face(id).surface])) {
      return id;
    }
  }
  return geom::FaceId::null();
}

double half_cylinder_exact_volume() {
  // The x >= 0.5 part of the unit cylinder of height 2: 2 * (pi/2 - 0.956611477...).
  return 1.228369698608801;
}

geom::Brep half_cylinder() {
  const geom::Brep cylinder =
      geom::make_cylinder(geom::Point3d{0.0, 0.0, 0.0}, geom::Vec3d{0.0, 0.0, 1.0}, 1.0, 2.0,
                          1e-9)
          .value();
  const geom::Brep box =
      geom::make_box(geom::Point3d{-2.0, -2.0, -1.0}, geom::Point3d{0.5, 2.0, 3.0}, 1e-9).value();
  return geom::boolean_solid_curved(cylinder, cylinder.solid_ids().front(), box,
                                    box.solid_ids().front(), geom::SolidBooleanOp::difference,
                                    geom::CurvedBooleanOptions{1e-9})
      .value()
      .brep;
}

// A cylinder cut by a slanted plane: the trim curve is an ellipse that is not
// rectilinear in the cylinder uv chart, so the strip decomposition cannot cover it.
geom::Brep slanted_cylinder(bool upper) {
  const geom::Brep cylinder =
      geom::make_cylinder(geom::Point3d{0.0, 0.0, 0.0}, geom::Vec3d{0.0, 0.0, 1.0}, 1.0, 2.0,
                          1e-9)
          .value();
  const geom::Vec3d normal{0.3, 0.0, 1.0};
  const double length = std::sqrt(geom::dot(normal, normal));
  const geom::Plane3 plane{geom::Vec3d{normal[0] / length, normal[1] / length, normal[2] / length},
                           1.0 / length};
  const auto split =
      geom::split_solid_by_plane(cylinder, cylinder.solid_ids().front(), plane, 1e-9);
  REQUIRE(split.has_value());
  return upper ? split.value().above : split.value().below;
}

double max_wall_offset(const geom::TriMesh& mesh, int axis, double value) {
  double worst = 0.0;
  for (const geom::Point3d& point : mesh.positions) {
    const double offset = std::abs(point[axis] - value);
    if (offset > 1e-6) continue;
    if (offset > worst) worst = offset;
  }
  return worst;
}

// A window fully inside the lateral face (does not touch the seam), produced by a
// box that pierces the +y wall: the lateral face then carries an inner hole loop.
geom::Brep hole_window_cylinder() {
  const geom::Brep cylinder =
      geom::make_cylinder(geom::Point3d{0.0, 0.0, 0.0}, geom::Vec3d{0.0, 0.0, 1.0}, 1.0, 2.0,
                          1e-9)
          .value();
  const geom::Brep box =
      geom::make_box(geom::Point3d{-0.5, 0.5, 0.5}, geom::Point3d{0.5, 2.0, 1.5}, 1e-9).value();
  return geom::boolean_solid_curved(cylinder, cylinder.solid_ids().front(), box,
                                    box.solid_ids().front(), geom::SolidBooleanOp::difference,
                                    geom::CurvedBooleanOptions{1e-9})
      .value()
      .brep;
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

TEST_CASE("trimmed tessellation: canonical primitive output hashes are stable") {
  // Hashes pin the byte-exact output of the canonical builders: the trimmed
  // fallback must never take over a face the analytic builder accepts.
  const geom::Brep box =
      geom::make_box(geom::Point3d{-1.0, -1.0, -1.0}, geom::Point3d{1.0, 1.0, 1.0}, 1e-9).value();
  const geom::Brep cylinder =
      geom::make_cylinder(geom::Point3d{0.0, 0.0, 0.0}, geom::Vec3d{0.0, 0.0, 1.0}, 1.0, 2.0,
                          1e-9)
          .value();
  const geom::Brep torus =
      geom::make_torus(geom::Point3d{0.0, 0.0, 0.0}, geom::Vec3d{0.0, 0.0, 1.0}, 2.0, 0.5, 1e-9)
          .value();
  const geom::TessellationOptions options{1e-2, 0.35};
  const auto box_mesh = geom::tessellate(box, options);
  const auto cylinder_mesh = geom::tessellate(cylinder, options);
  const auto torus_mesh = geom::tessellate(torus, options);
  REQUIRE(box_mesh.has_value());
  REQUIRE(cylinder_mesh.has_value());
  REQUIRE(torus_mesh.has_value());
  CHECK(box_mesh.value().triangles.size() == 12u);
  CHECK(box_mesh.value().positions.size() == 8u);
  CHECK(hash_mesh(box_mesh.value()) == 0x2abc6bb21fd6d800ull);
  CHECK(cylinder_mesh.value().triangles.size() == 88u);
  CHECK(cylinder_mesh.value().positions.size() == 46u);
  CHECK(hash_mesh(cylinder_mesh.value()) == 0x687a72048f8f4296ull);
  CHECK(torus_mesh.value().triangles.size() == 1152u);
  CHECK(torus_mesh.value().positions.size() == 576u);
  CHECK(hash_mesh(torus_mesh.value()) == 0xc808f6ee4d9c72dbull);
}

TEST_CASE("trimmed tessellation: windowed cylinder area is exact and outward") {
  const geom::Brep brep = cut_cylinder();
  const geom::FaceId lateral = lateral_face(brep);
  REQUIRE(!lateral.is_null());
  const auto projection =
      geom::detail::face_uv_projection(brep, brep.face(lateral), 1e-9);
  REQUIRE(projection.has_value());
  const double pi = std::numbers::pi_v<double>;
  const double ring_area = 4.0 * pi - pi / 3.0;
  const auto exact = geom::detail::curved_solid_volume(brep, brep.solid_ids().front(), 1e-9);
  REQUIRE(exact.has_value());
  CHECK(exact.value() == doctest::Approx(5.826573829689068).epsilon(1e-9));
  for (const double deflection : {1e-2, 1e-3}) {
    const auto mesh = geom::tessellate(brep, geom::TessellationOptions{deflection, 0.35});
    REQUIRE(mesh.has_value());
    check_watertight(mesh.value());
    double uv_sum = 0.0;
    double area_3d_sum = 0.0;
    for (std::size_t t = 0; t < mesh.value().triangles.size(); ++t) {
      if (mesh.value().triangle_face[t] != lateral.index()) continue;
      const auto& triangle = mesh.value().triangles[t];
      geom::Point2d uv[3];
      for (int k = 0; k < 3; ++k) {
        const auto projected = geom::detail::face_uv_project(
            projection.value(),
            mesh.value().positions[triangle[static_cast<std::size_t>(k)]], 0.0, false);
        REQUIRE(projected.has_value());
        uv[k] = projected.value();
      }
      for (int k = 1; k < 3; ++k) {
        while (uv[k][0] - uv[0][0] > pi) uv[k][0] -= 2.0 * pi;
        while (uv[0][0] - uv[k][0] > pi) uv[k][0] += 2.0 * pi;
      }
      const double twice = (uv[1][0] - uv[0][0]) * (uv[2][1] - uv[0][1]) -
                           (uv[1][1] - uv[0][1]) * (uv[2][0] - uv[0][0]);
      CHECK(twice > 0.0);
      uv_sum += 0.5 * twice;
      for (int k = 0; k < 3; ++k) {
        for (int j = k + 1; j < 3; ++j) {
          CHECK(std::abs(uv[k][0] - uv[j][0]) < pi / 2.0);
        }
      }
      const geom::Point3d& p = mesh.value().positions[triangle[0]];
      const geom::Point3d& q = mesh.value().positions[triangle[1]];
      const geom::Point3d& r = mesh.value().positions[triangle[2]];
      const geom::Vec3d normal =
          geom::cross(geom::Vec3d{q[0] - p[0], q[1] - p[1], q[2] - p[2]},
                      geom::Vec3d{r[0] - p[0], r[1] - p[1], r[2] - p[2]});
      const double area_3d = 0.5 * std::sqrt(geom::dot(normal, normal));
      area_3d_sum += area_3d;
      const geom::Vec3d centroid{(p[0] + q[0] + r[0]) / 3.0, (p[1] + q[1] + r[1]) / 3.0, 0.0};
      CHECK(geom::dot(normal, centroid) > 0.0);
      for (int k = 0; k < 3; ++k) {
        const geom::Point3d& a = mesh.value().positions[triangle[static_cast<std::size_t>(k)]];
        const geom::Point3d& b =
            mesh.value().positions[triangle[static_cast<std::size_t>((k + 1) % 3)]];
        CHECK(std::abs(std::hypot(0.5 * (a[0] + b[0]), 0.5 * (a[1] + b[1])) - 1.0) <=
              deflection);
      }
    }
    CHECK(uv_sum == doctest::Approx(ring_area).epsilon(1e-9));
    // The inscribed mesh may not gain area: only the total is comparable, since a
    // chord triangle can exceed its UV triangle near the window corners.
    CHECK(area_3d_sum <= uv_sum + 1e-12);
    const double volume = mesh_volume(mesh.value());
    CHECK(std::abs(volume - exact.value()) <= 2.0 * deflection * exact.value());
  }
}

TEST_CASE("trimmed tessellation: a half cylinder sliver stays watertight") {
  const geom::Brep brep = half_cylinder();
  REQUIRE(geom::validate_topology(brep).ok());
  const auto exact = geom::detail::curved_solid_volume(brep, brep.solid_ids().front(), 1e-9);
  REQUIRE(exact.has_value());
  CHECK(exact.value() == doctest::Approx(half_cylinder_exact_volume()).epsilon(1e-6));
  for (const double deflection : {1e-2, 1e-3}) {
    const auto mesh = geom::tessellate(brep, geom::TessellationOptions{deflection, 0.35});
    REQUIRE(mesh.has_value());
    check_watertight(mesh.value());
    CHECK(std::abs(mesh_volume(mesh.value()) - exact.value()) <=
          2.0 * deflection * exact.value());
  }
}
TEST_CASE("trimmed tessellation: a slanted trim reports the trimmed-path error") {
  // The ellipse boundary is not rectilinear in the uv chart, so this domain is
  // rejected by the surface-following check. The error must describe the trimmed
  // path, never the analytic primitive.
  for (const bool upper : {true, false}) {
    const geom::Brep brep = slanted_cylinder(upper);
    REQUIRE(geom::validate_topology(brep).ok());
    for (const double deflection : {1e-2, 1e-3}) {
      const auto mesh = geom::tessellate(brep, geom::TessellationOptions{deflection, 0.35});
      REQUIRE(!mesh.has_value());
      CHECK(mesh.error().code == geom::Errc::unsupported);
      CHECK(std::string(mesh.error().message) ==
            "tessellate: trimmed face does not follow the surface");
      CHECK(std::string(mesh.error().message) != "tessellate: cylinder face is missing a cap circle");
    }
  }
}

TEST_CASE("trimmed tessellation: trim samples sit exactly on the wall planes") {
  const geom::Brep windowed = cut_cylinder();
  const auto window_mesh = geom::tessellate(windowed, geom::TessellationOptions{1e-3, 0.35});
  REQUIRE(window_mesh.has_value());
  CHECK(max_wall_offset(window_mesh.value(), 1, 0.5) <= 1e-12);
  CHECK(max_wall_offset(window_mesh.value(), 1, -0.5) <= 1e-12);
  CHECK(max_wall_offset(window_mesh.value(), 2, 0.5) <= 1e-12);
  CHECK(max_wall_offset(window_mesh.value(), 2, 1.5) <= 1e-12);

  const geom::Brep sliver = half_cylinder();
  const auto sliver_mesh = geom::tessellate(sliver, geom::TessellationOptions{1e-3, 0.35});
  REQUIRE(sliver_mesh.has_value());
  CHECK(max_wall_offset(sliver_mesh.value(), 0, 0.5) <= 1e-12);
}

TEST_CASE("trimmed tessellation: the sag refinement subdivides the sliver band") {
  // Refinement is live: grinding the band from 1e-2 to 1e-3 must add triangles
  // beyond what the coarser sampling alone produces. Disabling the refinement
  // leaves the count unchanged and fails this check.
  const geom::Brep brep = half_cylinder();
  const auto coarse = geom::tessellate(brep, geom::TessellationOptions{1e-2, 0.35});
  const auto fine = geom::tessellate(brep, geom::TessellationOptions{1e-3, 0.35});
  REQUIRE(coarse.has_value());
  REQUIRE(fine.has_value());
  CHECK(fine.value().triangles.size() >=
        static_cast<std::size_t>(1.5 * static_cast<double>(coarse.value().triangles.size())));
}
TEST_CASE("trimmed tessellation: a hole loop is rejected with the trimmed-path error") {
  // The window sits fully inside the lateral face, so the face carries two loops.
  // The strip path covers rings without holes, and the generic path's chords
  // exceed the deviation bound, so the domain is rejected rather than remeshed
  // with a chordal patch. This pins the honest error for the 2-loop fixture.
  const geom::Brep brep = hole_window_cylinder();
  REQUIRE(geom::validate_topology(brep).ok());
  const geom::FaceId lateral = lateral_face(brep);
  REQUIRE(!lateral.is_null());
  const auto loops = geom::face_loops_uv(brep, lateral, 1e-9);
  REQUIRE(loops.has_value());
  REQUIRE(loops.value().size() == 2);
  for (const double deflection : {1e-2, 1e-3}) {
    const auto mesh = geom::tessellate(brep, geom::TessellationOptions{deflection, 0.35});
    REQUIRE(!mesh.has_value());
    CHECK(mesh.error().code == geom::Errc::unsupported);
    CHECK(std::string(mesh.error().message) ==
          "tessellate: trimmed face does not follow the surface");
  }
}
TEST_CASE("trimmed tessellation: a band-shaped domain refines with smaller deflection") {
  const geom::Brep brep = cut_cylinder();
  const auto coarse = geom::tessellate(brep, geom::TessellationOptions{5e-2, 0.35});
  const auto fine = geom::tessellate(brep, geom::TessellationOptions{2e-3, 0.35});
  REQUIRE(coarse.has_value());
  REQUIRE(fine.has_value());
  CHECK(fine.value().triangles.size() > coarse.value().triangles.size());
  check_watertight(coarse.value());
  check_watertight(fine.value());
}

TEST_CASE("trimmed tessellation: a band-shaped domain chord height meets deflection") {
  const geom::Brep brep = cut_cylinder();
  const double deflection = 2e-3;
  const auto mesh = geom::tessellate(brep, geom::TessellationOptions{deflection, 0.35});
  REQUIRE(mesh.has_value());
  const geom::FaceId lateral = lateral_face(brep);
  REQUIRE(!lateral.is_null());
  double worst = 0.0;
  for (std::size_t t = 0; t < mesh.value().triangles.size(); ++t) {
    if (mesh.value().triangle_face[t] != lateral.index()) continue;
    const geom::Point3d& a = mesh.value().positions[mesh.value().triangles[t][0]];
    const geom::Point3d& b = mesh.value().positions[mesh.value().triangles[t][1]];
    const geom::Point3d& c = mesh.value().positions[mesh.value().triangles[t][2]];
    const geom::Point3d points[3] = {a, b, c};
    for (int e = 0; e < 3; ++e) {
      const geom::Point3d& p = points[e];
      const geom::Point3d& q = points[(e + 1) % 3];
      const double mid[2] = {(p[0] + q[0]) / 2.0, (p[1] + q[1]) / 2.0};
      worst = std::max(worst, std::abs(std::hypot(mid[0], mid[1]) - 1.0));
    }
  }
  CHECK(worst <= deflection * 1.05);
}

TEST_CASE("trimmed tessellation: analytic failure does not fall back without a polyline edge") {
  // A curved face bounded only by line edges: the analytic cylinder builder
  // rejects it, and the dispatcher must surface that analytic error instead of
  // silently trying the uv path (which would fail with a different message).
  geom::SolidEditor editor;
  const std::uint32_t surface = editor.add_surface(
      geom::Cylinder3{geom::Point3d{0.0, 0.0, 0.0}, geom::Vec3d{0.0, 0.0, 1.0},
                      geom::Vec3d{1.0, 0.0, 0.0}, 1.0});
  const geom::Point3d corners[3] = {geom::Point3d{1.0, 0.0, 0.0}, geom::Point3d{0.0, 1.0, 0.0},
                                    geom::Point3d{0.0, 0.0, 1.0}};
  std::vector<geom::VertexId> vertices;
  for (const geom::Point3d& corner : corners) {
    vertices.push_back(editor.add_vertex(corner, 0.0).value());
  }
  std::vector<geom::CoedgeSpec> specs;
  for (std::size_t i = 0; i < 3; ++i) {
    const geom::Point3d from = corners[i];
    const geom::Point3d to = corners[(i + 1) % 3];
    const geom::Vec3d delta = to - from;
    const double length = std::hypot(delta[0], delta[1], delta[2]);
    const std::uint32_t curve = editor.add_curve(
        geom::LineCurve3{from, geom::Vec3d{delta[0] / length, delta[1] / length, delta[2] / length}});
    const geom::EdgeId edge =
        editor.add_edge(vertices[i], vertices[(i + 1) % 3], curve, 0.0, length, 0.0).value();
    specs.push_back(geom::CoedgeSpec{edge, false});
  }
  const geom::FaceId face = editor.add_face(surface, false, 0.0).value();
  REQUIRE(editor.add_loop(face, specs, true).has_value());
  const auto mesh = geom::tessellate(editor.brep(), {});
  REQUIRE(!mesh.has_value());
  CHECK(mesh.error().code == geom::Errc::invalid_topology);
  CHECK(std::string(mesh.error().message) == "tessellate: cylinder face is missing a cap circle");
}

TEST_CASE("QA-13: trimmed tessellation is watertight, finite and deterministic") {
  const geom::Brep brep = cut_cylinder();
  const geom::TessellationOptions options{2e-3, 0.35};
  const auto first = geom::tessellate(brep, options);
  const auto second = geom::tessellate(brep, options);
  REQUIRE(first.has_value());
  REQUIRE(second.has_value());
  const geom::TriMesh& a = first.value();
  const geom::TriMesh& b = second.value();
  check_watertight(a);
  CHECK(a.triangles == b.triangles);
  CHECK(a.triangle_face == b.triangle_face);
  CHECK(a.positions.size() == b.positions.size());
  for (const auto& triangle : a.triangles) {
    for (const std::uint32_t index : triangle) {
      CHECK(index < a.positions.size());
    }
  }
  for (std::size_t i = 0; i < a.positions.size(); ++i) {
    for (int axis = 0; axis < 3; ++axis) {
      CHECK(std::isfinite(a.positions[i][static_cast<std::size_t>(axis)]));
      double left = a.positions[i][static_cast<std::size_t>(axis)];
      double right = b.positions[i][static_cast<std::size_t>(axis)];
      CHECK(std::memcmp(&left, &right, sizeof(double)) == 0);
    }
  }
}
