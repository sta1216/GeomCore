#include <doctest/doctest.h>

#include <array>
#include <cmath>
#include <cstdint>
#include <limits>
#include <optional>
#include <span>
#include <string>
#include <string_view>
#include <vector>

#include <geom/brep/boolean_solid.hpp>
#include <geom/brep/chamfer.hpp>
#include <geom/brep/mass_properties.hpp>
#include <geom/brep/primitives.hpp>
#include <geom/brep/solid_editor.hpp>
#include <geom/brep/validate.hpp>
#include <geom/core/error.hpp>

using namespace geom;

namespace {

Brep make_box_solid(Point3d lo, Point3d hi) {
  const Result<Brep> box = make_box(lo, hi, 1e-9);
  REQUIRE(box.has_value());
  return box.value();
}

double volume_of(const Brep& brep) {
  const Result<MassProperties> mass = mass_properties(brep);
  REQUIRE(mass.has_value());
  return mass.value().volume;
}

// Finds the straight edge whose endpoints are the two given corners.
EdgeId edge_between(const Brep& brep, Point3d a, Point3d b) {
  for (const EdgeId id : brep.edge_ids()) {
    const Point3d p0 = brep.vertex(brep.edge(id).v0).point;
    const Point3d p1 = brep.vertex(brep.edge(id).v1).point;
    const auto same = [](Point3d p, Point3d q) {
      return std::abs(p[0] - q[0]) + std::abs(p[1] - q[1]) +
                 std::abs(p[2] - q[2]) <
             1e-12;
    };
    if ((same(p0, a) && same(p1, b)) || (same(p0, b) && same(p1, a))) return id;
  }
  REQUIRE_MESSAGE(false, "edge not found");
  return {};
}

}  // namespace

TEST_CASE("chamfer one box edge removes a triangular prism") {
  const Brep box = make_box_solid(Point3d{0.0, 0.0, 0.0}, Point3d{2.0, 3.0, 4.0});
  const EdgeId edge =
      edge_between(box, Point3d{0.0, 0.0, 0.0}, Point3d{0.0, 0.0, 4.0});
  const std::array<EdgeId, 1> edges{edge};

  ChamferOptions options;
  options.distance = 0.25;
  const Result<ChamferResult> result = chamfer_edges(box, edges, options);
  REQUIRE(result.has_value());
  CHECK(result.value().chamfered_edges == 1u);
  CHECK(result.value().brep.solid_count() == 1u);
  CHECK(result.value().brep.shell_count() == 1u);
  CHECK(validate_topology(result.value().brep).ok());
  // Removed cross-section is d^2/2 = 0.03125 over the 4-long edge.
  CHECK(std::abs(volume_of(result.value().brep) - 23.875) <= 1e-9);
  // The input document is untouched.
  CHECK(box.face_count() == 6u);
  CHECK(std::abs(volume_of(box) - 24.0) <= 1e-12);
}

TEST_CASE("chamfer rejects a non-positive distance") {
  const Brep box = make_box_solid(Point3d{0.0, 0.0, 0.0}, Point3d{1.0, 1.0, 1.0});
  const std::array<EdgeId, 1> edges{box.edge_ids().front()};
  ChamferOptions options;
  options.distance = 0.0;
  const Result<ChamferResult> result = chamfer_edges(box, edges, options);
  REQUIRE_FALSE(result.has_value());
  CHECK(result.error().code == Errc::degenerate);
  CHECK(std::string_view(result.error().message).find("distance must be positive") !=
        std::string_view::npos);
}

// Guard: the shared planar-convex helper Task 2/Task 3 build on must describe a
// box as six outward planes with unit normals, without touching the input.
TEST_CASE("chamfer helper describes a box as six planes") {
  const Brep box = make_box_solid(Point3d{0.0, 0.0, 0.0}, Point3d{2.0, 3.0, 4.0});
  const Result<detail::BlendPlanarSolid> view =
      detail::blend_planar_convex_solid(box, 1e-9, "chamfer_edges");
  REQUIRE(view.has_value());
  CHECK(view.value().solid == box.solid_ids().front());
  CHECK(view.value().faces.size() == 6u);
  CHECK(view.value().scale == 4.0);
  for (const Vec3d outward : view.value().outward) {
    CHECK(std::abs(norm2(outward) - 1.0) <= 1e-12);
  }
}

TEST_CASE("chamfer two edges sharing a vertex subtracts their overlap once") {
  const Brep box = make_box_solid(Point3d{0.0, 0.0, 0.0}, Point3d{2.0, 2.0, 2.0});
  const std::array<EdgeId, 2> edges{
      edge_between(box, Point3d{0.0, 0.0, 0.0}, Point3d{0.0, 0.0, 2.0}),
      edge_between(box, Point3d{0.0, 0.0, 0.0}, Point3d{0.0, 2.0, 0.0})};
  ChamferOptions options;
  options.distance = 0.25;
  const Result<ChamferResult> result = chamfer_edges(box, edges, options);
  REQUIRE(result.has_value());
  CHECK(result.value().chamfered_edges == 2u);
  CHECK(validate_topology(result.value().brep).ok());
  // Two prisms of 0.0625 each, their d^3/3 corner overlap removed once.
  CHECK(std::abs(volume_of(result.value().brep) - 7.880208333333333) <= 1e-9);
}

TEST_CASE("chamfer all twelve cube edges leaves the analytic volume") {
  const Brep box = make_box_solid(Point3d{0.0, 0.0, 0.0}, Point3d{2.0, 2.0, 2.0});
  std::vector<EdgeId> edges = box.edge_ids();
  REQUIRE(edges.size() == 12u);
  ChamferOptions options;
  options.distance = 0.25;
  const Result<ChamferResult> result = chamfer_edges(box, edges, options);
  REQUIRE(result.has_value());
  CHECK(result.value().chamfered_edges == 12u);
  CHECK(validate_topology(result.value().brep).ok());
  // Inclusion-exclusion over the twelve prisms: a^3 - 6*a*d^2 + 6*d^3.
  CHECK(std::abs(volume_of(result.value().brep) - 7.34375) <= 1e-9);
}

TEST_CASE("chamfer counts a repeated edge once") {
  const Brep box = make_box_solid(Point3d{0.0, 0.0, 0.0}, Point3d{1.0, 1.0, 1.0});
  const EdgeId edge = edge_between(box, Point3d{0.0, 0.0, 0.0}, Point3d{0.0, 0.0, 1.0});
  const std::array<EdgeId, 2> edges{edge, edge};
  ChamferOptions options;
  options.distance = 0.1;
  const Result<ChamferResult> result = chamfer_edges(box, edges, options);
  REQUIRE(result.has_value());
  CHECK(result.value().chamfered_edges == 1u);
}

TEST_CASE("chamfer rejects a curved operand") {
  const Result<Brep> cylinder = make_cylinder(Point3d{0.0, 0.0, 0.0},
                                              Vec3d{0.0, 0.0, 1.0}, 1.0, 2.0, 1e-9);
  REQUIRE(cylinder.has_value());
  const std::array<EdgeId, 1> edges{cylinder.value().edge_ids().front()};
  ChamferOptions options;
  options.distance = 0.1;
  const Result<ChamferResult> result = chamfer_edges(cylinder.value(), edges, options);
  REQUIRE_FALSE(result.has_value());
  CHECK(result.error().code == Errc::unsupported);
  CHECK(std::string_view(result.error().message).find("all faces must be planar") !=
        std::string_view::npos);
}

TEST_CASE("chamfer rejects a non-convex solid") {
  const Brep a = make_box_solid(Point3d{0.0, 0.0, 0.0}, Point3d{2.0, 1.0, 1.0});
  const Brep b = make_box_solid(Point3d{0.0, 0.5, 0.0}, Point3d{1.0, 1.5, 1.0});
  const Result<SolidBooleanResult> formed =
      boolean_solid(a, b, SolidBooleanOp::union_);
  REQUIRE(formed.has_value());
  const std::array<EdgeId, 1> edges{formed.value().brep.edge_ids().front()};
  // Guard the fixture: the rejection below must come from the union's shape,
  // not from a stale handle, so the chosen edge must be alive.
  CHECK(formed.value().brep.edge_alive(edges[0]));
  ChamferOptions options;
  options.distance = 0.1;
  const Result<ChamferResult> result = chamfer_edges(formed.value().brep, edges, options);
  REQUIRE_FALSE(result.has_value());
  CHECK(result.error().code == Errc::unsupported);
  CHECK(std::string_view(result.error().message).find("solid must be convex") !=
        std::string_view::npos);
}

TEST_CASE("chamfer rejects a distance at the face width") {
  const Brep box = make_box_solid(Point3d{0.0, 0.0, 0.0}, Point3d{2.0, 2.0, 2.0});
  const std::array<EdgeId, 1> edges{
      edge_between(box, Point3d{0.0, 0.0, 0.0}, Point3d{0.0, 0.0, 2.0})};
  ChamferOptions options;
  options.distance = 2.0;
  const Result<ChamferResult> result = chamfer_edges(box, edges, options);
  REQUIRE_FALSE(result.has_value());
  CHECK(result.error().code == Errc::unsupported);
  CHECK(std::string_view(result.error().message).find("exceeds the adjacent face width") !=
        std::string_view::npos);
}

TEST_CASE("chamfer rejects an empty edge list and a dead edge") {
  const Brep box = make_box_solid(Point3d{0.0, 0.0, 0.0}, Point3d{1.0, 1.0, 1.0});
  ChamferOptions options;
  options.distance = 0.1;
  const std::array<EdgeId, 0> none{};
  const Result<ChamferResult> empty = chamfer_edges(box, none, options);
  REQUIRE_FALSE(empty.has_value());
  CHECK(empty.error().code == Errc::unsupported);
  CHECK(std::string_view(empty.error().message).find("at least one edge is required") !=
        std::string_view::npos);

  const std::array<EdgeId, 1> dead{EdgeId{}};
  const Result<ChamferResult> bad = chamfer_edges(box, dead, options);
  REQUIRE_FALSE(bad.has_value());
  CHECK(bad.error().code == Errc::index_out_of_range);
  CHECK(std::string_view(bad.error().message).find("edge does not belong to the solid") !=
        std::string_view::npos);
}

TEST_CASE("chamfer rejects an option error before touching the solid") {
  const Brep box = make_box_solid(Point3d{0.0, 0.0, 0.0}, Point3d{1.0, 1.0, 1.0});
  const std::array<EdgeId, 1> edges{box.edge_ids().front()};
  ChamferOptions options;
  options.distance = std::numeric_limits<double>::quiet_NaN();
  const Result<ChamferResult> non_finite = chamfer_edges(box, edges, options);
  REQUIRE_FALSE(non_finite.has_value());
  CHECK(non_finite.error().code == Errc::non_finite);
  CHECK(std::string_view(non_finite.error().message).find("distance must be finite") !=
        std::string_view::npos);

  options.distance = 0.1;
  options.tolerance = 0.0;
  const Result<ChamferResult> bad_tolerance = chamfer_edges(box, edges, options);
  REQUIRE_FALSE(bad_tolerance.has_value());
  CHECK(bad_tolerance.error().code == Errc::tolerance_violation);
  CHECK(std::string_view(bad_tolerance.error().message).find("tolerance must be positive") !=
        std::string_view::npos);
}

namespace {

// Builds a one-solid, one-shell open shell through the public editor API: a
// single square face whose shell is flagged open. This reaches the shape
// stage's closed-shell rejection without private access.
Result<Brep> make_open_shell_solid() {
  SolidEditor editor;
  const Point3d points[4] = {
      {0.0, 0.0, 0.0}, {1.0, 0.0, 0.0}, {1.0, 1.0, 0.0}, {0.0, 1.0, 0.0}};
  std::array<VertexId, 4> vertices{};
  for (int i = 0; i < 4; ++i) vertices[i] = editor.add_vertex(points[i], 0.0).value();
  const std::uint32_t plane = editor.add_surface(Plane3{{0.0, 0.0, 1.0}, 0.0});
  const std::uint32_t line =
      editor.add_curve(LineCurve3{{0.0, 0.0, 0.0}, {1.0, 0.0, 0.0}});
  std::array<EdgeId, 4> edges{};
  for (int i = 0; i < 4; ++i) {
    edges[i] =
        editor.add_edge(vertices[i], vertices[(i + 1) % 4], line, 0.0, 1.0, 0.0).value();
  }
  const FaceId face = editor.add_face(plane, false, 0.0).value();
  std::array<CoedgeSpec, 4> specs{};
  for (int i = 0; i < 4; ++i) specs[i] = CoedgeSpec{edges[i], false};
  editor.add_loop(face, specs, true).value();
  const ShellId shell = editor.add_shell(std::array<FaceId, 1>{face}, false).value();
  editor.add_solid(std::array<ShellId, 1>{shell}, 0.0).value();
  return editor.commit();
}

}  // namespace

TEST_CASE("chamfer rejects an empty solid") {
  ChamferOptions options;
  options.distance = 0.1;
  const std::array<EdgeId, 1> edges{EdgeId{}};
  const Result<ChamferResult> shape = chamfer_edges(Brep{}, edges, options);
  REQUIRE_FALSE(shape.has_value());
  CHECK(shape.error().code == Errc::unsupported);
  CHECK(std::string_view(shape.error().message)
            .find("exactly one solid with one shell is required") != std::string_view::npos);

  // The option stage runs before the shape stage: the same empty operand with
  // a non-finite distance reports the option error, not the missing shell.
  options.distance = std::numeric_limits<double>::quiet_NaN();
  const Result<ChamferResult> ordered = chamfer_edges(Brep{}, edges, options);
  REQUIRE_FALSE(ordered.has_value());
  CHECK(ordered.error().code == Errc::non_finite);
  CHECK(std::string_view(ordered.error().message).find("distance must be finite") !=
        std::string_view::npos);
}

TEST_CASE("chamfer rejects an open shell") {
  const Result<Brep> open = make_open_shell_solid();
  REQUIRE(open.has_value());
  CHECK(open.value().solid_count() == 1u);
  CHECK(open.value().shell_count() == 1u);
  const std::array<EdgeId, 1> edges{open.value().edge_ids().front()};
  ChamferOptions options;
  options.distance = 0.1;
  const Result<ChamferResult> result = chamfer_edges(open.value(), edges, options);
  REQUIRE_FALSE(result.has_value());
  CHECK(result.error().code == Errc::invalid_topology);
  CHECK(std::string_view(result.error().message).find("solid shell must be closed") !=
        std::string_view::npos);
}

// The union of two overlapping boxes is one valid convex solid whose coincident
// inner faces are not merged, so its seam edges border two coplanar faces. At
// distance 0.01 this fixture has 8 such edges out of 28.
TEST_CASE("chamfer rejects coplanar adjacent faces") {
  const Brep a = make_box_solid(Point3d{0.0, 0.0, 0.0}, Point3d{1.0, 1.0, 1.0});
  const Brep b = make_box_solid(Point3d{0.5, 0.0, 0.0}, Point3d{2.0, 1.0, 1.0});
  const Result<SolidBooleanResult> formed = boolean_solid(a, b, SolidBooleanOp::union_);
  REQUIRE(formed.has_value());
  const Brep& united = formed.value().brep;
  REQUIRE(validate_topology(united).ok());

  std::optional<EdgeId> coplanar;
  std::uint32_t coplanar_count = 0;
  for (const EdgeId id : united.edge_ids()) {
    ChamferOptions options;
    options.distance = 0.01;
    const Result<ChamferResult> one =
        chamfer_edges(united, std::array<EdgeId, 1>{id}, options);
    if (one.has_value()) continue;
    if (std::string_view(one.error().message).find("adjacent faces are coplanar") ==
        std::string_view::npos) {
      continue;
    }
    ++coplanar_count;
    if (!coplanar.has_value()) coplanar = id;
  }
  REQUIRE(coplanar.has_value());
  CHECK(coplanar_count >= 1u);

  ChamferOptions options;
  options.distance = 0.01;
  const Result<ChamferResult> result =
      chamfer_edges(united, std::array<EdgeId, 1>{coplanar.value()}, options);
  REQUIRE_FALSE(result.has_value());
  CHECK(result.error().code == Errc::unsupported);
  CHECK(std::string_view(result.error().message).find("adjacent faces are coplanar") !=
        std::string_view::npos);
}

TEST_CASE("chamfer rejects a non-finite and a negative tolerance") {
  const Brep box = make_box_solid(Point3d{0.0, 0.0, 0.0}, Point3d{1.0, 1.0, 1.0});
  const std::array<EdgeId, 1> edges{box.edge_ids().front()};
  ChamferOptions options;
  options.distance = 0.1;
  options.tolerance = std::numeric_limits<double>::quiet_NaN();
  const Result<ChamferResult> non_finite = chamfer_edges(box, edges, options);
  REQUIRE_FALSE(non_finite.has_value());
  CHECK(non_finite.error().code == Errc::non_finite);
  CHECK(std::string_view(non_finite.error().message).find("tolerance must be finite") !=
        std::string_view::npos);

  options.tolerance = -1e-9;
  const Result<ChamferResult> negative = chamfer_edges(box, edges, options);
  REQUIRE_FALSE(negative.has_value());
  CHECK(negative.error().code == Errc::tolerance_violation);
  CHECK(std::string_view(negative.error().message).find("tolerance must be positive") !=
        std::string_view::npos);
}
