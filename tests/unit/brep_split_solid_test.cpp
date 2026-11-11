#include <doctest/doctest.h>

#include <cmath>
#include <cstddef>
#include <cstdint>
#include <numbers>
#include <optional>
#include <string_view>
#include <variant>
#include <vector>

#include <geom/brep/detail/curved_volume.hpp>
#include <geom/brep/face_uv.hpp>
#include <geom/brep/primitives.hpp>
#include <geom/brep/solid_editor.hpp>
#include <geom/brep/split_solid.hpp>
#include <geom/brep/validate.hpp>
#include <geom/core/error.hpp>
#include <geom/curves/analytic_curve.hpp>
#include <geom/surfaces/analytic.hpp>
#include <geom/surfaces/plane.hpp>

namespace {

constexpr double kPi = std::numbers::pi_v<double>;

// Exact-arc fixture for the split tests. A quarter cylinder of radius
// `radius` and height `height` whose bottom rim lies at z = 0 and whose
// lateral surface is the cylinder (origin, +z, +x x_axis). Both rims are true
// CircleCurve3 arcs, so the volume of either half of a z = height / 2 cut is
// exactly pi * radius^2 * height / 8; the chordal `geom::test` fixture cannot
// represent that (its top cap is a 16-gon).
geom::Brep make_exact_quarter_cylinder(double radius, double height, double tol) {
  constexpr int kSegments = 8;
  const double delta = 0.5 * kPi / static_cast<double>(kSegments);

  geom::SolidEditor editor;

  std::vector<geom::VertexId> bottom;
  std::vector<geom::VertexId> top;
  bottom.reserve(kSegments + 1);
  top.reserve(kSegments + 1);
  for (int i = 0; i <= kSegments; ++i) {
    const double angle = static_cast<double>(i) * delta;
    const double x = radius * std::cos(angle);
    const double y = radius * std::sin(angle);
    bottom.push_back(editor.add_vertex(geom::Point3d{x, y, 0.0}, tol).value());
    top.push_back(editor.add_vertex(geom::Point3d{x, y, height}, tol).value());
  }
  const geom::VertexId base = editor.add_vertex(geom::Point3d{0.0, 0.0, 0.0}, tol).value();
  const geom::VertexId apex =
      editor.add_vertex(geom::Point3d{0.0, 0.0, height}, tol).value();

  const geom::CircleCurve3 bottom_circle{geom::Point3d{0.0, 0.0, 0.0},
                                         geom::Vec3d{1.0, 0.0, 0.0},
                                         geom::Vec3d{0.0, 1.0, 0.0}, radius};
  const geom::CircleCurve3 top_circle{geom::Point3d{0.0, 0.0, height},
                                      geom::Vec3d{1.0, 0.0, 0.0},
                                      geom::Vec3d{0.0, 1.0, 0.0}, radius};
  const std::uint32_t bottom_curve = editor.add_curve(bottom_circle);
  const std::uint32_t top_curve = editor.add_curve(top_circle);

  std::vector<geom::EdgeId> bottom_arcs;
  std::vector<geom::EdgeId> top_arcs;
  bottom_arcs.reserve(kSegments);
  top_arcs.reserve(kSegments);
  for (int i = 0; i < kSegments; ++i) {
    const double t0 = static_cast<double>(i) * delta;
    const double t1 = static_cast<double>(i + 1) * delta;
    bottom_arcs.push_back(
        editor.add_edge(bottom[i], bottom[i + 1], bottom_curve, t0, t1, tol).value());
    top_arcs.push_back(editor.add_edge(top[i], top[i + 1], top_curve, t0, t1, tol).value());
  }

  auto line = [&](geom::VertexId v0, geom::VertexId v1, geom::Point3d origin,
                  geom::Vec3d direction, double length) {
    const std::uint32_t curve = editor.add_curve(geom::LineCurve3{origin, direction});
    return editor.add_edge(v0, v1, curve, 0.0, length, tol).value();
  };

  const geom::Point3d seam0_start{radius, 0.0, 0.0};
  const geom::Point3d seam1_start{0.0, radius, 0.0};
  const geom::EdgeId seam0 = line(bottom[0], top[0], seam0_start, geom::Vec3d{0.0, 0.0, 1.0},
                                  height);
  const geom::EdgeId seam1 = line(bottom[kSegments], top[kSegments], seam1_start,
                                  geom::Vec3d{0.0, 0.0, 1.0}, height);
  const geom::EdgeId bottom_x =
      line(base, bottom[0], geom::Point3d{0.0, 0.0, 0.0}, geom::Vec3d{1.0, 0.0, 0.0}, radius);
  const geom::EdgeId bottom_y = line(base, bottom[kSegments], geom::Point3d{0.0, 0.0, 0.0},
                                     geom::Vec3d{0.0, 1.0, 0.0}, radius);
  const geom::EdgeId top_x =
      line(apex, top[0], geom::Point3d{0.0, 0.0, height}, geom::Vec3d{1.0, 0.0, 0.0}, radius);
  const geom::EdgeId top_y = line(apex, top[kSegments], geom::Point3d{0.0, 0.0, height},
                                  geom::Vec3d{0.0, 1.0, 0.0}, radius);
  const geom::EdgeId axis =
      line(base, apex, geom::Point3d{0.0, 0.0, 0.0}, geom::Vec3d{0.0, 0.0, 1.0}, height);

  using geom::CoedgeSpec;
  using geom::FaceId;

  const std::uint32_t lateral_surface = editor.add_surface(
      geom::Cylinder3{geom::Point3d{0.0, 0.0, 0.0}, geom::Vec3d{0.0, 0.0, 1.0},
                      geom::Vec3d{1.0, 0.0, 0.0}, radius});
  const std::uint32_t bottom_surface =
      editor.add_surface(geom::Plane3{geom::Vec3d{0.0, 0.0, -1.0}, 0.0});
  const std::uint32_t top_surface =
      editor.add_surface(geom::Plane3{geom::Vec3d{0.0, 0.0, 1.0}, height});
  const std::uint32_t side_y_surface =
      editor.add_surface(geom::Plane3{geom::Vec3d{0.0, -1.0, 0.0}, 0.0});
  const std::uint32_t side_x_surface =
      editor.add_surface(geom::Plane3{geom::Vec3d{-1.0, 0.0, 0.0}, 0.0});

  const FaceId lateral = editor.add_face(lateral_surface, false, tol).value();
  const FaceId bottom_cap = editor.add_face(bottom_surface, false, tol).value();
  const FaceId top_cap = editor.add_face(top_surface, false, tol).value();
  const FaceId side_y = editor.add_face(side_y_surface, false, tol).value();
  const FaceId side_x = editor.add_face(side_x_surface, false, tol).value();

  std::vector<CoedgeSpec> lateral_loop;
  for (const geom::EdgeId edge : bottom_arcs) lateral_loop.push_back(CoedgeSpec{edge});
  lateral_loop.push_back(CoedgeSpec{seam1});
  for (std::size_t i = top_arcs.size(); i-- > 0;) {
    lateral_loop.push_back(CoedgeSpec{top_arcs[i], true});
  }
  lateral_loop.push_back(CoedgeSpec{seam0, true});
  CHECK(editor.add_loop(lateral, lateral_loop, true).has_value());

  std::vector<CoedgeSpec> bottom_loop;
  bottom_loop.push_back(CoedgeSpec{bottom_y});
  for (std::size_t i = bottom_arcs.size(); i-- > 0;) {
    bottom_loop.push_back(CoedgeSpec{bottom_arcs[i], true});
  }
  bottom_loop.push_back(CoedgeSpec{bottom_x, true});
  CHECK(editor.add_loop(bottom_cap, bottom_loop, true).has_value());

  std::vector<CoedgeSpec> top_loop;
  for (const geom::EdgeId edge : top_arcs) top_loop.push_back(CoedgeSpec{edge});
  top_loop.push_back(CoedgeSpec{top_y, true});
  top_loop.push_back(CoedgeSpec{top_x});
  CHECK(editor.add_loop(top_cap, top_loop, true).has_value());

  std::vector<CoedgeSpec> side_y_loop;
  side_y_loop.push_back(CoedgeSpec{bottom_x});
  side_y_loop.push_back(CoedgeSpec{seam0});
  side_y_loop.push_back(CoedgeSpec{top_x, true});
  side_y_loop.push_back(CoedgeSpec{axis, true});
  CHECK(editor.add_loop(side_y, side_y_loop, true).has_value());

  std::vector<CoedgeSpec> side_x_loop;
  side_x_loop.push_back(CoedgeSpec{axis});
  side_x_loop.push_back(CoedgeSpec{top_y});
  side_x_loop.push_back(CoedgeSpec{seam1, true});
  side_x_loop.push_back(CoedgeSpec{bottom_y, true});
  CHECK(editor.add_loop(side_x, side_x_loop, true).has_value());

  const std::vector<FaceId> faces{lateral, bottom_cap, top_cap, side_y, side_x};
  const geom::ShellId shell = editor.add_shell(faces, true).value();
  const std::vector<geom::ShellId> shells{shell};
  CHECK(editor.add_solid(shells, tol).has_value());
  return editor.commit().value();
}

double volume_of(const geom::Brep& brep) {
  const std::vector<geom::SolidId> solids = brep.solid_ids();
  REQUIRE(solids.size() == 1);
  const geom::Result<double> volume = geom::detail::curved_solid_volume(brep, solids.front(), 1e-9);
  REQUIRE(volume.has_value());
  return volume.value();
}

// Number of coedges that reference a CircleCurve3 with the given axis height
// and radius. A cut circle healed into an exact arc edge appears exactly twice
// (once on the trimmed lateral face, once on the cap).
std::size_t cut_circle_uses(const geom::Brep& brep, double height, double radius) {
  std::size_t uses = 0;
  for (const geom::CoedgeId id : brep.coedge_ids()) {
    const geom::Edge& edge = brep.edge(brep.coedge(id).edge);
    if (edge.curve == geom::kNoGeometry) continue;
    const geom::CurveGeometry& curve = brep.geometry().curves[edge.curve];
    if (!std::holds_alternative<geom::CircleCurve3>(curve)) continue;
    const geom::CircleCurve3& circle = std::get<geom::CircleCurve3>(curve);
    if (std::abs(circle.center[2] - height) <= 1e-9 &&
        std::abs(circle.radius - radius) <= 1e-9) {
      ++uses;
    }
  }
  return uses;
}

// No loop may carry a full-turn arc span: the assembler rejects those, so a
// solid that assembles at all proves the split, but this asserts the span
// contract directly as well.
void check_no_full_turn_spans(const geom::Brep& brep) {
  for (const geom::FaceId face : brep.face_ids()) {
    const auto loops = geom::face_loops_uv(brep, face, 1e-9);
    REQUIRE(loops.has_value());
    for (const geom::FaceUvLoop& loop : loops.value()) {
      for (const std::optional<geom::FaceUvArcRef>& arc : loop.arcs) {
        if (!arc.has_value()) continue;
        CHECK(std::abs(arc->t1 - arc->t0) < kPi);
      }
    }
  }
}

}  // namespace

TEST_CASE("brep: split_solid_by_plane cuts a quarter cylinder at z = 1") {
  const double tolerance = 1e-9;
  const geom::Brep solid = make_exact_quarter_cylinder(1.0, 2.0, tolerance);
  REQUIRE(geom::validate_topology(solid).ok());

  const std::vector<geom::SolidId> solids = solid.solid_ids();
  REQUIRE(solids.size() == 1);
  const geom::Result<geom::SplitSolidResult> split =
      geom::split_solid_by_plane(solid, solids.front(),
                                 geom::Plane3{geom::Vec3d{0.0, 0.0, 1.0}, 1.0}, tolerance);
  REQUIRE(split.has_value());

  const geom::Brep& below = split.value().below;
  const geom::Brep& above = split.value().above;
  CHECK(geom::validate_topology(below).ok());
  CHECK(geom::validate_topology(above).ok());

  const double expected = kPi * 1.0 * 1.0 * 2.0 / 8.0;
  CHECK(std::abs(volume_of(below) - expected) <= 1e-9);
  CHECK(std::abs(volume_of(above) - expected) <= 1e-9);

  // The cut boundary is an exact circle shared by the trimmed lateral face and
  // the cap; a chordal approximation would show up as polyline edges instead.
  CHECK(cut_circle_uses(below, 1.0, 1.0) == 2);
  CHECK(cut_circle_uses(above, 1.0, 1.0) == 2);
}

TEST_CASE("brep: split_solid_by_plane cuts a box at z = 0.5") {
  const double tolerance = 1e-9;
  const geom::Brep solid =
      geom::make_box(geom::Point3d{0.0, 0.0, 0.0}, geom::Point3d{4.0, 2.0, 1.0}, tolerance).value();
  const std::vector<geom::SolidId> solids = solid.solid_ids();
  REQUIRE(solids.size() == 1);

  const geom::Result<geom::SplitSolidResult> split = geom::split_solid_by_plane(
      solid, solids.front(), geom::Plane3{geom::Vec3d{0.0, 0.0, 1.0}, 0.5}, tolerance);
  REQUIRE(split.has_value());
  CHECK(geom::validate_topology(split.value().below).ok());
  CHECK(geom::validate_topology(split.value().above).ok());
  CHECK(std::abs(volume_of(split.value().below) - 4.0) <= 1e-12);
  CHECK(std::abs(volume_of(split.value().above) - 4.0) <= 1e-12);
}

TEST_CASE("brep: split_solid_by_plane cuts a full cylinder at half height") {
  const double tolerance = 1e-9;
  const geom::Brep solid =
      geom::make_cylinder(geom::Point3d{0.0, 0.0, 0.0}, geom::Vec3d{0.0, 0.0, 1.0}, 1.0, 2.0,
                          tolerance)
          .value();
  const std::vector<geom::SolidId> solids = solid.solid_ids();
  REQUIRE(solids.size() == 1);

  const geom::Result<geom::SplitSolidResult> split = geom::split_solid_by_plane(
      solid, solids.front(), geom::Plane3{geom::Vec3d{0.0, 0.0, 1.0}, 1.0}, tolerance);
  REQUIRE(split.has_value());
  CHECK(geom::validate_topology(split.value().below).ok());
  CHECK(geom::validate_topology(split.value().above).ok());

  const double expected = kPi * 1.0 * 1.0 * 2.0 / 2.0;
  CHECK(std::abs(volume_of(split.value().below) - expected) <= 1e-9);
  CHECK(std::abs(volume_of(split.value().above) - expected) <= 1e-9);

  // Both halves expose the full cut circle as an exact arc edge.
  CHECK(cut_circle_uses(split.value().below, 1.0, 1.0) == 2);
  CHECK(cut_circle_uses(split.value().above, 1.0, 1.0) == 2);

  // A single 2*pi span is direction-ambiguous; the full cut circle must be
  // carried as several shorter spans.
  check_no_full_turn_spans(split.value().below);
  check_no_full_turn_spans(split.value().above);
}

TEST_CASE("brep: split_solid_by_plane rejects a tangent plane") {
  const double tolerance = 1e-9;
  const geom::Brep solid = make_exact_quarter_cylinder(1.0, 2.0, tolerance);
  const std::vector<geom::SolidId> solids = solid.solid_ids();
  REQUIRE(solids.size() == 1);

  const geom::Result<geom::SplitSolidResult> split = geom::split_solid_by_plane(
      solid, solids.front(), geom::Plane3{geom::Vec3d{0.0, 0.0, 1.0}, 2.0}, tolerance);
  REQUIRE(!split.has_value());
  CHECK(split.error().code == geom::Errc::unsupported);
  CHECK(std::string_view(split.error().message) ==
        "split_solid_by_plane: plane is tangent or coincident with a face");
}
