#include <doctest/doctest.h>

#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <numbers>
#include <optional>
#include <string_view>
#include <utility>
#include <variant>
#include <vector>

#include <geom/brep/boolean_solid.hpp>
#include <geom/brep/boolean_solid_curved.hpp>
#include <geom/brep/detail/curved_volume.hpp>
#include <geom/brep/face_uv.hpp>
#include <geom/brep/primitives.hpp>
#include <geom/brep/solid_editor.hpp>
#include <geom/brep/split_solid.hpp>
#include <geom/brep/validate.hpp>
#include <geom/core/error.hpp>
#include <geom/curves/analytic_curve.hpp>
#include <geom/surfaces/analytic.hpp>
#include <geom/surfaces/nurbs_surface.hpp>
#include <geom/surfaces/plane.hpp>

namespace {

constexpr double kPi = std::numbers::pi_v<double>;

// Exact-arc fixture for the curved boolean tests. A quarter cylinder of radius
// `radius` and height `height` whose bottom rim lies at z = 0 and whose lateral
// surface is the cylinder (origin, +z, +x x_axis). Both rims are true
// CircleCurve3 arcs, so the volume of the lower half of a z = height / 2 cut is
// exactly pi * radius^2 * height / 8.
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

// A closed box whose bottom face carries a planar degree-1 NURBS patch instead
// of a Plane3. The topology matches make_box, so the operand validates and only
// the surface-kind check of the driver can reject it.
geom::Brep make_nurbs_box(geom::Point3d lo, geom::Point3d hi, double tol) {
  const double x0 = lo[0];
  const double y0 = lo[1];
  const double z0 = lo[2];
  const double x1 = hi[0];
  const double y1 = hi[1];
  const double z1 = hi[2];
  const geom::Point3d corners[8] = {
      {x0, y0, z0}, {x1, y0, z0}, {x1, y1, z0}, {x0, y1, z0},
      {x0, y0, z1}, {x1, y0, z1}, {x1, y1, z1}, {x0, y1, z1},
  };

  geom::SolidEditor editor;
  geom::VertexId vertices[8];
  for (int i = 0; i < 8; ++i) {
    vertices[i] = editor.add_vertex(corners[i], tol).value();
  }

  struct EdgeSpec {
    int v0;
    int v1;
  };
  const EdgeSpec edge_specs[12] = {
      {0, 1}, {1, 2}, {2, 3}, {3, 0},
      {4, 5}, {5, 6}, {6, 7}, {7, 4},
      {1, 5}, {2, 6}, {3, 7}, {0, 4},
  };
  geom::EdgeId edges[12];
  for (int i = 0; i < 12; ++i) {
    const geom::Point3d from = corners[edge_specs[i].v0];
    const geom::Point3d to = corners[edge_specs[i].v1];
    const geom::Vec3d delta = to - from;
    const double length = std::hypot(delta[0], delta[1], delta[2]);
    const geom::Vec3d direction{delta[0] / length, delta[1] / length, delta[2] / length};
    const std::uint32_t curve = editor.add_curve(geom::LineCurve3{from, direction});
    edges[i] = editor.add_edge(vertices[edge_specs[i].v0], vertices[edge_specs[i].v1], curve,
                               0.0, length, tol)
                   .value();
  }

  struct CoedgeUse {
    int edge;
    bool reversed;
  };
  struct FaceSpec {
    geom::Vec3d normal;
    int corner;
    CoedgeUse loop[4];
  };
  const FaceSpec face_specs[6] = {
      {{0.0, 0.0, -1.0}, 0, {{3, true}, {2, true}, {1, true}, {0, true}}},
      {{0.0, 0.0, 1.0}, 4, {{4, false}, {5, false}, {6, false}, {7, false}}},
      {{0.0, -1.0, 0.0}, 0, {{0, false}, {8, false}, {4, true}, {11, true}}},
      {{0.0, 1.0, 0.0}, 2, {{2, false}, {10, false}, {6, true}, {9, true}}},
      {{-1.0, 0.0, 0.0}, 0, {{11, false}, {7, true}, {10, true}, {3, false}}},
      {{1.0, 0.0, 0.0}, 1, {{1, false}, {9, false}, {5, true}, {8, true}}},
  };

  const geom::NurbsSurface patch(
      1, 1, {0.0, 0.0, 1.0, 1.0}, {0.0, 0.0, 1.0, 1.0},
      {geom::Point3d{x0, y0, z0}, geom::Point3d{x0, y1, z0}, geom::Point3d{x1, y0, z0},
       geom::Point3d{x1, y1, z0}},
      {1.0, 1.0, 1.0, 1.0});

  std::array<geom::FaceId, 6> faces{};
  for (int f = 0; f < 6; ++f) {
    const FaceSpec& spec = face_specs[f];
    const geom::Point3d corner = corners[spec.corner];
    const std::uint32_t surface =
        f == 0 ? editor.add_surface(patch)
               : editor.add_surface(geom::Plane3{
                     spec.normal, geom::dot(spec.normal, geom::Vec3d{corner[0], corner[1],
                                                                     corner[2]})});
    faces[static_cast<std::size_t>(f)] = editor.add_face(surface, false, tol).value();
    std::array<geom::CoedgeSpec, 4> coedges{};
    for (int i = 0; i < 4; ++i) {
      coedges[static_cast<std::size_t>(i)] =
          geom::CoedgeSpec{edges[spec.loop[i].edge], spec.loop[i].reversed};
    }
    CHECK(editor.add_loop(faces[static_cast<std::size_t>(f)], coedges, true).has_value());
  }

  const geom::ShellId shell = editor.add_shell(faces, true).value();
  const std::array<geom::ShellId, 1> shells{shell};
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

geom::SolidId only_solid(const geom::Brep& brep) {
  const std::vector<geom::SolidId> solids = brep.solid_ids();
  REQUIRE(solids.size() == 1);
  return solids.front();
}

std::string_view message_of(const geom::Error& error) { return std::string_view(error.message); }

}  // namespace

TEST_CASE("brep: boolean_solid_curved subtracts a box from a quarter cylinder") {
  const double tolerance = 1e-9;
  const geom::Brep quarter = make_exact_quarter_cylinder(1.0, 2.0, tolerance);
  // The box covers z >= 1 and has no face on x = 0 or y = 0, so the quarter
  // cylinder's two radial faces are cut instead of lying coplanar with a box
  // face.
  const geom::Brep box =
      geom::make_box(geom::Point3d{-1.0, -1.0, 1.0}, geom::Point3d{10.0, 10.0, 10.0},
                     tolerance)
          .value();
  REQUIRE(geom::validate_topology(quarter).ok());
  REQUIRE(geom::validate_topology(box).ok());

  // Reference: the same cut taken by the plane splitter.
  const geom::Result<geom::SplitSolidResult> split = geom::split_solid_by_plane(
      quarter, only_solid(quarter), geom::Plane3{geom::Vec3d{0.0, 0.0, 1.0}, 1.0}, tolerance);
  REQUIRE(split.has_value());
  const geom::Brep& reference = split.value().below;
  REQUIRE(geom::validate_topology(reference).ok());
  const double expected = volume_of(reference);

  const geom::Result<geom::CurvedBooleanResult> result = geom::boolean_solid_curved(
      quarter, only_solid(quarter), box, only_solid(box), geom::SolidBooleanOp::difference,
      geom::CurvedBooleanOptions{tolerance});
  REQUIRE_MESSAGE(result.has_value(), "boolean_solid_curved failed: " << message_of(result.error()));
  const geom::Brep& brep = result.value().brep;
  CHECK(geom::validate_topology(brep).ok());

  const double actual = volume_of(brep);
  CHECK(std::abs(actual - expected) <= 1e-6 * expected);
  // The cut circle is analytic: a chordal cap would leave the trimmed box face
  // as a polygon and the volume far below the reference.
  CHECK(std::abs(actual - kPi / 4.0) <= 1e-6);
  CHECK(brep.face_ids().size() == reference.face_ids().size());

  const geom::CurvedBooleanReport& report = result.value().report;
  CHECK(report.fragments == report.kept + report.dropped);
  CHECK(report.kept == static_cast<std::uint32_t>(brep.face_ids().size()));
  CHECK(report.max_cut_residual <= tolerance);

  // The counters are a pure function of the inputs.
  const geom::Result<geom::CurvedBooleanResult> again = geom::boolean_solid_curved(
      quarter, only_solid(quarter), box, only_solid(box), geom::SolidBooleanOp::difference,
      geom::CurvedBooleanOptions{tolerance});
  REQUIRE(again.has_value());
  CHECK(again.value().report.face_pairs == report.face_pairs);
  CHECK(again.value().report.cuts == report.cuts);
  CHECK(again.value().report.fragments == report.fragments);
  CHECK(again.value().report.kept == report.kept);
  CHECK(again.value().report.dropped == report.dropped);
  CHECK(again.value().report.flipped == report.flipped);
  CHECK(again.value().report.max_cut_residual == report.max_cut_residual);
}

TEST_CASE("brep: boolean_solid_curved unions two boxes like boolean_solid") {
  const double tolerance = 1e-9;
  const geom::Brep a =
      geom::make_box(geom::Point3d{0.0, 0.0, 0.0}, geom::Point3d{2.0, 2.0, 2.0}, tolerance)
          .value();
  const geom::Brep b =
      geom::make_box(geom::Point3d{1.0, 1.0, 1.0}, geom::Point3d{3.0, 3.0, 3.0}, tolerance)
          .value();

  const geom::Result<geom::SolidBooleanResult> reference =
      geom::boolean_solid(a, b, geom::SolidBooleanOp::union_, geom::SolidBooleanOptions{tolerance});
  REQUIRE(reference.has_value());
  REQUIRE(geom::validate_topology(reference.value().brep).ok());

  const geom::Result<geom::CurvedBooleanResult> result = geom::boolean_solid_curved(
      a, only_solid(a), b, only_solid(b), geom::SolidBooleanOp::union_,
      geom::CurvedBooleanOptions{tolerance});
  REQUIRE(result.has_value());
  CHECK(geom::validate_topology(result.value().brep).ok());

  const double expected = volume_of(reference.value().brep);
  const double actual = volume_of(result.value().brep);
  CHECK(std::abs(actual - expected) <= 1e-9);
  CHECK(std::abs(actual - 15.0) <= 1e-9);
  CHECK(result.value().report.kept ==
        static_cast<std::uint32_t>(result.value().brep.face_ids().size()));
}

TEST_CASE("brep: boolean_solid_curved rejects NURBS operands") {
  const double tolerance = 1e-9;
  const geom::Brep nurbs =
      make_nurbs_box(geom::Point3d{0.0, 0.0, 0.0}, geom::Point3d{1.0, 1.0, 1.0}, tolerance);
  REQUIRE(geom::validate_topology(nurbs).ok());
  const geom::Brep box =
      geom::make_box(geom::Point3d{0.5, 0.5, 0.5}, geom::Point3d{2.0, 2.0, 2.0}, tolerance)
          .value();

  const geom::Result<geom::CurvedBooleanResult> result = geom::boolean_solid_curved(
      nurbs, only_solid(nurbs), box, only_solid(box), geom::SolidBooleanOp::union_,
      geom::CurvedBooleanOptions{tolerance});
  REQUIRE(!result.has_value());
  CHECK(result.error().code == geom::Errc::unsupported);
  CHECK(message_of(result.error()) == "boolean_solid_curved: NURBS operands are not supported");
}

TEST_CASE("brep: boolean_solid_curved rejects cone operands") {
  const double tolerance = 1e-9;
  const geom::Brep cone =
      geom::make_cone(geom::Point3d{0.0, 0.0, 0.0}, geom::Vec3d{0.0, 0.0, 1.0}, 1.0, 1.0,
                      tolerance)
          .value();
  REQUIRE(geom::validate_topology(cone).ok());
  const geom::Brep box =
      geom::make_box(geom::Point3d{-2.0, -2.0, 0.0}, geom::Point3d{2.0, 2.0, 2.0}, tolerance)
          .value();

  const geom::Result<geom::CurvedBooleanResult> result = geom::boolean_solid_curved(
      cone, only_solid(cone), box, only_solid(box), geom::SolidBooleanOp::union_,
      geom::CurvedBooleanOptions{tolerance});
  REQUIRE(!result.has_value());
  CHECK(result.error().code == geom::Errc::unsupported);
  CHECK(message_of(result.error()) ==
        "boolean_solid_curved: cone operands are not supported (apex-touching loops are deferred)");
}

TEST_CASE("brep: boolean_solid_curved rejects coplanar face contact") {
  const double tolerance = 1e-9;
  const geom::Brep a =
      geom::make_box(geom::Point3d{0.0, 0.0, 0.0}, geom::Point3d{1.0, 1.0, 1.0}, tolerance)
          .value();
  // b's left face lies in a's right face plane and covers it.
  const geom::Brep b =
      geom::make_box(geom::Point3d{1.0, -1.0, -1.0}, geom::Point3d{2.0, 2.0, 2.0}, tolerance)
          .value();
  REQUIRE(geom::validate_topology(a).ok());
  REQUIRE(geom::validate_topology(b).ok());

  const geom::Result<geom::CurvedBooleanResult> result = geom::boolean_solid_curved(
      a, only_solid(a), b, only_solid(b), geom::SolidBooleanOp::union_,
      geom::CurvedBooleanOptions{tolerance});
  REQUIRE(!result.has_value());
  CHECK(result.error().code == geom::Errc::unsupported);
  CHECK(message_of(result.error()) ==
        "boolean_solid_curved: coplanar face contact is not supported yet");
}

TEST_CASE("brep: boolean_solid_curved rejects a disconnected result") {
  const double tolerance = 1e-9;
  const geom::Brep a =
      geom::make_box(geom::Point3d{0.0, 0.0, 0.0}, geom::Point3d{1.0, 1.0, 1.0}, tolerance)
          .value();
  const geom::Brep b =
      geom::make_box(geom::Point3d{3.0, 0.0, 0.0}, geom::Point3d{4.0, 1.0, 1.0}, tolerance)
          .value();
  REQUIRE(geom::validate_topology(a).ok());
  REQUIRE(geom::validate_topology(b).ok());

  const geom::Result<geom::CurvedBooleanResult> result = geom::boolean_solid_curved(
      a, only_solid(a), b, only_solid(b), geom::SolidBooleanOp::union_,
      geom::CurvedBooleanOptions{tolerance});
  REQUIRE(!result.has_value());
  CHECK(result.error().code == geom::Errc::unsupported);
  CHECK(message_of(result.error()) ==
        "boolean_solid_curved: result is not a single connected component");
}

TEST_CASE("brep: boolean_solid_curved heals the cut rim into exact arcs") {
  const double tolerance = 1e-9;
  const geom::Brep quarter = make_exact_quarter_cylinder(1.0, 2.0, tolerance);
  const geom::Result<geom::Brep> box_result =
      geom::make_box(geom::Point3d{-1.0, -1.0, 1.0}, geom::Point3d{10.0, 10.0, 10.0}, tolerance);
  REQUIRE(box_result.has_value());
  const geom::Brep& box = box_result.value();
  const geom::Result<geom::CurvedBooleanResult> result = geom::boolean_solid_curved(
      quarter, only_solid(quarter), box, only_solid(box), geom::SolidBooleanOp::difference,
      geom::CurvedBooleanOptions{tolerance});
  REQUIRE_MESSAGE(result.has_value(),
                  "boolean_solid_curved failed: " << message_of(result.error()));
  // The trimmed box cap at z = 1 is new boundary: its rim against the lateral
  // face keeps the exact area the volume check above relies on only when the
  // chordal cut is healed back onto the cylinder's circle. A healed cap is
  // emitted with the outward normal baked into its plane, so it is the only +z
  // facing plane of the result and its rings carry engaged arc references.
  bool healed_cap = false;
  for (const geom::FaceId id : result.value().brep.face_ids()) {
    const geom::Face& face = result.value().brep.face(id);
    const geom::SurfaceGeometry& surface = result.value().brep.geometry().surfaces[face.surface];
    const geom::Plane3* plane = std::get_if<geom::Plane3>(&surface);
    if (plane == nullptr || !(plane->normal[2] > 0.5)) continue;
    CHECK_FALSE(face.reversed);
    const geom::Result<std::vector<geom::FaceUvLoop>> loops =
        geom::face_loops_uv(result.value().brep, id, tolerance);
    REQUIRE(loops.has_value());
    std::size_t arcs = 0;
    for (const geom::FaceUvLoop& loop : loops.value()) {
      for (const std::optional<geom::FaceUvArcRef>& arc : loop.arcs) {
        if (arc.has_value()) ++arcs;
      }
    }
    CHECK(arcs > 0);
    healed_cap = true;
  }
  CHECK(healed_cap);
}