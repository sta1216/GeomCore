#include <doctest/doctest.h>

#include <algorithm>
#include <array>
#include <bit>
#include <limits>
#include <chrono>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <numbers>
#include <optional>
#include <random>
#include <string>
#include <string_view>
#include <utility>
#include <variant>
#include <vector>

#include "brep_boolean_curved_fixture.hpp"

#include <geom/brep/boolean_solid.hpp>
#include <geom/brep/boolean_solid_curved.hpp>
#include <geom/brep/detail/curved_volume.hpp>
#include <geom/brep/face_uv.hpp>
#include <geom/brep/heal.hpp>
#include <geom/brep/primitives.hpp>
#include <geom/brep/solid_editor.hpp>
#include <geom/brep/split_solid.hpp>
#include <geom/brep/tessellate.hpp>
#include <geom/brep/validate.hpp>
#include <geom/core/error.hpp>
#include <geom/curves/analytic_curve.hpp>
#include <geom/surfaces/analytic.hpp>
#include <geom/surfaces/nurbs_surface.hpp>
#include <geom/surfaces/plane.hpp>

namespace {

constexpr double kPi = std::numbers::pi_v<double>;

using geom::test::sphere_face_cap_volume;
using geom::test::sphere_face_volume;

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

// Structural fingerprint of an operand. The sphere/sphere tests compare one
// before and after a rejected boolean so a diagnostic cannot quietly re-chart or
// otherwise mutate its inputs.
struct BrepFingerprint {
  double volume = 0.0;
  std::size_t vertices = 0;
  std::size_t edges = 0;
  std::size_t faces = 0;
  std::size_t shells = 0;
  std::size_t solids = 0;
  bool operator==(const BrepFingerprint& other) const = default;
};

BrepFingerprint fingerprint_of(const geom::Brep& brep) {
  BrepFingerprint result;
  result.volume = volume_of(brep);
  result.vertices = brep.vertex_ids().size();
  result.edges = brep.edge_ids().size();
  result.faces = brep.face_ids().size();
  result.solids = brep.solid_ids().size();
  result.shells = brep.solid_shells(brep.solid_ids().front()).size();
  return result;
}

// Full analytic cylinder r = 1, h = 2 with its bottom cap at z = 0. Unlike the
// quarter fixture above its lateral face carries the full periodic seam, so its
// cuts have to land on the seam edges.
geom::Brep make_test_cylinder(double tolerance) {
  const geom::Result<geom::Brep> result =
      geom::make_cylinder(geom::Point3d{0.0, 0.0, 0.0}, geom::Vec3d{0.0, 0.0, 1.0}, 1.0, 2.0,
                          tolerance);
  REQUIRE(result.has_value());
  return result.value();
}

// QA corpus helpers. A positive case must produce a validated single solid whose
// volume matches the analytic value, and a repeated call must reproduce that
// volume bit for bit. A rejection must reproduce the same code and message.
void check_qa_volume(const char* name, const geom::Brep& a, const geom::Brep& b,
                     geom::SolidBooleanOp op, double expected) {
  INFO(name);
  const geom::Result<geom::CurvedBooleanResult> first = geom::boolean_solid_curved(
      a, only_solid(a), b, only_solid(b), op, geom::CurvedBooleanOptions{1e-9});
  REQUIRE_MESSAGE(first.has_value(), name << ": " << message_of(first.error()));
  CHECK_MESSAGE(geom::validate_topology(first.value().brep).ok(), name << ": invalid topology");
  const double volume = volume_of(first.value().brep);
  CHECK_MESSAGE(std::abs(volume - expected) <= 1e-6 * std::abs(expected),
                name << ": volume " << volume << " expected " << expected);

  const geom::Result<geom::CurvedBooleanResult> again = geom::boolean_solid_curved(
      a, only_solid(a), b, only_solid(b), op, geom::CurvedBooleanOptions{1e-9});
  REQUIRE_MESSAGE(again.has_value(), name << ": the repeated call failed");
  CHECK_MESSAGE(volume_of(again.value().brep) == volume, name << ": volume is not deterministic");
}

void check_qa_rejection(const char* name, const geom::Brep& a, const geom::Brep& b,
                        geom::SolidBooleanOp op, std::string_view message) {
  INFO(name);
  const geom::Result<geom::CurvedBooleanResult> first = geom::boolean_solid_curved(
      a, only_solid(a), b, only_solid(b), op, geom::CurvedBooleanOptions{1e-9});
  REQUIRE_MESSAGE(!first.has_value(), name << ": expected a diagnostic");
  CHECK_MESSAGE(first.error().code == geom::Errc::unsupported, name << ": unexpected error code");
  CHECK_MESSAGE(!message_of(first.error()).empty(), name << ": empty message");
  CHECK_MESSAGE(message_of(first.error()) == message,
                name << ": message " << message_of(first.error()));

  const geom::Result<geom::CurvedBooleanResult> again = geom::boolean_solid_curved(
      a, only_solid(a), b, only_solid(b), op, geom::CurvedBooleanOptions{1e-9});
  REQUIRE_MESSAGE(!again.has_value(), name << ": the repeated call succeeded");
  CHECK(again.error().code == first.error().code);
  CHECK(message_of(again.error()) == message_of(first.error()));
}

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
        "boolean_solid_curved: cone operands need closed-cut and apex topology; booleans are "
        "deferred to M3c-3");
}

TEST_CASE("curved boolean: cone operands are deferred with a diagnostic") {
  const double tolerance = 1e-9;
  const geom::Result<geom::Brep> cone =
      geom::make_cone(geom::Point3d{0.0, 0.0, 0.0}, geom::Vec3d{0.0, 0.0, 1.0}, 1.0, 2.0,
                      tolerance);
  const geom::Result<geom::Brep> box =
      geom::make_box(geom::Point3d{-2.0, -2.0, -1.0}, geom::Point3d{2.0, 2.0, 2.0}, tolerance);
  REQUIRE(cone.has_value());
  REQUIRE(box.has_value());
  for (const geom::SolidBooleanOp op : {geom::SolidBooleanOp::union_,
                                        geom::SolidBooleanOp::intersection,
                                        geom::SolidBooleanOp::difference,
                                        geom::SolidBooleanOp::symmetric_difference}) {
    geom::CurvedBooleanReport report{};
    report.face_pairs = 999;
    const geom::Result<geom::CurvedBooleanResult> out = geom::boolean_solid_curved(
        cone.value(), only_solid(cone.value()), box.value(), only_solid(box.value()), op,
        geom::CurvedBooleanOptions{tolerance}, &report);
    REQUIRE_FALSE(out.has_value());
    CHECK(out.error().code == geom::Errc::unsupported);
    const std::string message(out.error().message);
    CHECK(message.find("deferred to M3c-3") != std::string::npos);
    CHECK(report.face_pairs == 0);
    CHECK(report.fragments == 0);
    CHECK(report.kept == 0);
  }
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

TEST_CASE("curved boolean: clipping keeps a sphere latitude cut") {
  const double tolerance = 1e-9;
  const geom::Result<geom::Brep> sphere =
      geom::make_sphere(geom::Point3d{0.0, 0.0, 0.0}, 1.0, tolerance);
  const geom::Result<geom::Brep> box = geom::make_box(
      geom::Point3d{-2.0, -2.0, 0.5}, geom::Point3d{2.0, 2.0, 2.0}, tolerance);
  REQUIRE(sphere.has_value());
  REQUIRE(box.has_value());
  const auto sphere_info = geom::detail::curved_boolean_face_info(
      sphere.value(), sphere.value().face_ids().front(), tolerance);
  const auto box_info = geom::detail::curved_boolean_face_info(
      box.value(), box.value().face_ids().front(), tolerance);
  REQUIRE(sphere_info.has_value());
  REQUIRE(box_info.has_value());

  std::vector<geom::Point2d> chain;
  constexpr std::size_t kSamples = 32;
  for (std::size_t i = 0; i <= kSamples; ++i) {
    const double u = 2.0 * 3.14159265358979323846 *
                     static_cast<double>(i) / static_cast<double>(kSamples);
    chain.push_back(geom::Point2d{u, std::asin(0.5)});
  }
  const std::vector<std::vector<geom::Point2d>> pieces =
      geom::detail::curved_boolean_clip_pair(sphere_info.value(), chain,
                                             box_info.value(), nullptr, tolerance);
  REQUIRE(pieces.size() == 1);
  CHECK(pieces.front().size() == chain.size());
}

TEST_CASE("curved boolean: a sphere cap keeps one exact pole point") {
  const double tolerance = 1e-9;
  const geom::Result<geom::Brep> sphere =
      geom::make_sphere(geom::Point3d{0.0, 0.0, 0.0}, 1.0, tolerance);
  REQUIRE(sphere.has_value());
  const geom::FaceId face = sphere.value().face_ids().front();

  std::vector<geom::Point2d> ring;
  constexpr std::size_t kSamples = 64;
  const double v = std::asin(0.5);
  for (std::size_t i = 0; i <= kSamples; ++i) {
    const double u = 2.0 * 3.14159265358979323846 *
                     static_cast<double>(i) / static_cast<double>(kSamples);
    ring.push_back(geom::Point2d{u, v});
  }
  const std::vector<geom::FaceCut> cuts{geom::FaceCut{ring}};
  const geom::Result<std::vector<geom::CurvedFragment>> fragments =
      geom::split_face_curved(sphere.value(), face, cuts, tolerance);
  REQUIRE(fragments.has_value());
  REQUIRE(fragments.value().size() == 2);

  const auto uv_area = [](const std::vector<geom::Point2d>& loop) {
    double twice = 0.0;
    for (std::size_t i = 0; i < loop.size(); ++i) {
      const geom::Point2d& a = loop[i];
      const geom::Point2d& b = loop[(i + 1) % loop.size()];
      twice += a[0] * b[1] - b[0] * a[1];
    }
    return std::fabs(twice) * 0.5;
  };
  bool found_cap = false;
  for (const geom::CurvedFragment& fragment : fragments.value()) {
    REQUIRE(fragment.loops.size() == fragment.loops_uv.size());
    for (std::size_t l = 0; l < fragment.loops_uv.size(); ++l) {
      CHECK(fragment.loops[l].size() == fragment.loops_uv[l].size());
      std::size_t poles = 0;
      for (const geom::Point2d& uv : fragment.loops_uv[l]) {
        if (uv[1] == 0.5 * 3.14159265358979323846) ++poles;
      }
      if (poles == 1) {
        found_cap = true;
        CHECK(uv_area(fragment.loops_uv[l]) > 0.1);
      }
    }
  }
  CHECK(found_cap);
}

TEST_CASE("curved boolean: a sphere band joins its two cut runs") {
  const double tolerance = 1e-9;
  const geom::Result<geom::Brep> sphere =
      geom::make_sphere(geom::Point3d{0.0, 0.0, 0.0}, 1.0, tolerance);
  REQUIRE(sphere.has_value());
  const geom::FaceId face = sphere.value().face_ids().front();

  constexpr std::size_t kSamples = 64;
  const double v = std::asin(0.8);
  std::vector<geom::FaceCut> cuts(2);
  for (std::size_t i = 0; i <= kSamples; ++i) {
    const double u = 2.0 * 3.14159265358979323846 *
                     static_cast<double>(i) / static_cast<double>(kSamples);
    cuts[0].uv.push_back(geom::Point2d{u, -v});
    cuts[1].uv.push_back(geom::Point2d{u, v});
  }
  const geom::Result<std::vector<geom::CurvedFragment>> fragments =
      geom::split_face_curved(sphere.value(), face, cuts, tolerance);
  REQUIRE(fragments.has_value());
  REQUIRE(fragments.value().size() == 3);

  const auto uv_area = [](const std::vector<geom::Point2d>& loop) {
    double twice = 0.0;
    for (std::size_t i = 0; i < loop.size(); ++i) {
      const geom::Point2d& a = loop[i];
      const geom::Point2d& b = loop[(i + 1) % loop.size()];
      twice += a[0] * b[1] - b[0] * a[1];
    }
    return std::fabs(twice) * 0.5;
  };
  bool found_band = false;
  for (const geom::CurvedFragment& fragment : fragments.value()) {
    if (fragment.loops_uv.size() != 1) continue;
    double v_min = 1e300;
    double v_max = -1e300;
    for (const geom::Point2d& uv : fragment.loops_uv.front()) {
      v_min = std::min(v_min, uv[1]);
      v_max = std::max(v_max, uv[1]);
    }
    if (v_min < -0.5 && v_max > 0.5) {
      found_band = true;
      CHECK(uv_area(fragment.loops_uv.front()) ==
            doctest::Approx(2.0 * 3.14159265358979323846 * 2.0 * v)
                .epsilon(1e-9));
    }
  }
  CHECK(found_band);
}

TEST_CASE("curved boolean: sphere cap intersection") {
  const double tolerance = 1e-9;
  const geom::Result<geom::Brep> sphere =
      geom::make_sphere(geom::Point3d{0.0, 0.0, 0.0}, 1.0, tolerance);
  const geom::Result<geom::Brep> box = geom::make_box(
      geom::Point3d{-2.0, -2.0, 0.5}, geom::Point3d{2.0, 2.0, 2.0}, tolerance);
  REQUIRE(sphere.has_value());
  REQUIRE(box.has_value());
  const geom::SolidId a = only_solid(sphere.value());
  const geom::SolidId b = only_solid(box.value());
  geom::CurvedBooleanReport report;
  const geom::Result<geom::CurvedBooleanResult> result = geom::boolean_solid_curved(
      sphere.value(), a, box.value(), b, geom::SolidBooleanOp::intersection,
      geom::CurvedBooleanOptions{tolerance}, &report);
  REQUIRE_MESSAGE(result.has_value(), "boolean failed: " << message_of(result.error()));
  CHECK(geom::validate_topology(result.value().brep).ok());
  // `report` is the failure report and stays zeroed on success; the returned
  // report carries the counters.
  CHECK(report.cuts == 0);
  CHECK(result.value().report.cuts >= 1);
  CHECK(result.value().report.fragments >= 2);
  CHECK(std::fabs(volume_of(result.value().brep) - 0.6544984694978736) <= 1e-6);
}

TEST_CASE("curved boolean: sphere complement cap") {
  const double tolerance = 1e-9;
  const geom::Result<geom::Brep> sphere =
      geom::make_sphere(geom::Point3d{0.0, 0.0, 0.0}, 1.0, tolerance);
  const geom::Result<geom::Brep> box = geom::make_box(
      geom::Point3d{-2.0, -2.0, 0.5}, geom::Point3d{2.0, 2.0, 2.0}, tolerance);
  REQUIRE(sphere.has_value());
  REQUIRE(box.has_value());
  const geom::Result<geom::CurvedBooleanResult> result = geom::boolean_solid_curved(
      sphere.value(), only_solid(sphere.value()), box.value(), only_solid(box.value()),
      geom::SolidBooleanOp::difference, geom::CurvedBooleanOptions{tolerance});
  REQUIRE_MESSAGE(result.has_value(), "boolean failed: " << message_of(result.error()));
  CHECK(geom::validate_topology(result.value().brep).ok());
  CHECK(std::fabs(volume_of(result.value().brep) - 3.534291735288517) <= 1e-6);
}

TEST_CASE("curved boolean: sphere slab band") {
  const double tolerance = 1e-9;
  const geom::Result<geom::Brep> sphere =
      geom::make_sphere(geom::Point3d{0.0, 0.0, 0.0}, 1.0, tolerance);
  const geom::Result<geom::Brep> slab = geom::make_box(
      geom::Point3d{-2.0, -2.0, -0.8}, geom::Point3d{2.0, 2.0, 0.8}, tolerance);
  REQUIRE(sphere.has_value());
  REQUIRE(slab.has_value());
  const geom::Result<geom::CurvedBooleanResult> result = geom::boolean_solid_curved(
      sphere.value(), only_solid(sphere.value()), slab.value(), only_solid(slab.value()),
      geom::SolidBooleanOp::intersection, geom::CurvedBooleanOptions{tolerance});
  REQUIRE_MESSAGE(result.has_value(), "boolean failed: " << message_of(result.error()));
  CHECK(geom::validate_topology(result.value().brep).ok());
  CHECK(result.value().brep.face_ids().size() == 3);
  CHECK(result.value().brep.edge_ids().size() == 3);
  CHECK(std::fabs(volume_of(result.value().brep) - 3.9542179533184135) <= 1e-6);
}

// QA-15: the three M3c-3a end-to-end sphere cases, each run twice. Every run
// must be a validated single-shell solid whose boundary closes, and the
// repeated run must reproduce the volume bit for bit.
TEST_CASE("QA-15: sphere cap booleans are closed and deterministic") {
  const double tolerance = 1e-9;
  const geom::Brep sphere =
      geom::make_sphere(geom::Point3d{0.0, 0.0, 0.0}, 1.0, tolerance).value();
  const geom::Brep cap_box =
      geom::make_box(geom::Point3d{-2.0, -2.0, 0.5}, geom::Point3d{2.0, 2.0, 2.0}, tolerance)
          .value();
  const geom::Brep slab_box =
      geom::make_box(geom::Point3d{-2.0, -2.0, -0.8}, geom::Point3d{2.0, 2.0, 0.8}, tolerance)
          .value();

  struct Case {
    const char* name;
    const geom::Brep* b;
    geom::SolidBooleanOp op;
    double expected;
    // The slab keeps one sphere band whose 8194-point uv outline joins its two
    // boundary circles across the seam; triangulate_polygon2 rejects that
    // outline, so the classification-mesh route is unavailable for that face.
    // Cap and complement close their classification mesh; all three close by
    // edge uses (every edge is used by exactly two coedges).
    bool mesh_closable;
  };
  const std::array<Case, 3> cases{{
      {"cap", &cap_box, geom::SolidBooleanOp::intersection, 0.6544984694978736, true},
      {"complement", &cap_box, geom::SolidBooleanOp::difference, 3.534291735288517, true},
      {"slab", &slab_box, geom::SolidBooleanOp::intersection, 3.9542179533184135, false},
  }};

  // `tessellate` rejects trimmed boolean faces, so close the check the way the
  // driver builds its classification mesh: project every face, weld and census.
  const auto check_mesh_closed = [&tolerance](const char* name, const geom::Brep& brep) {
    std::vector<geom::detail::CurvedBooleanFaceInfo> infos;
    for (const geom::FaceId id : brep.face_ids()) {
      const geom::Result<geom::detail::CurvedBooleanFaceInfo> info =
          geom::detail::curved_boolean_face_info(brep, id, tolerance);
      REQUIRE_MESSAGE(info.has_value(),
                      name << ": face info failed: " << message_of(info.error()));
      infos.push_back(info.value());
    }
    const geom::Result<geom::TriMesh> mesh =
        geom::detail::curved_boolean_solid_mesh(brep, infos, tolerance);
    REQUIRE_MESSAGE(mesh.has_value(),
                    name << ": classification mesh failed: " << message_of(mesh.error()));
    const geom::detail::CurvedBooleanMeshClosure closure =
        geom::detail::curved_boolean_mesh_closure(mesh.value());
    CHECK_MESSAGE(closure.closed(), name << ": mesh is not closed");
    CHECK_MESSAGE(closure.used_once == 0, name << ": mesh has boundary edges");
    CHECK_MESSAGE(closure.nonmanifold == 0, name << ": mesh is non-manifold");
  };

  for (const Case& entry : cases) {
    INFO(std::string(entry.name));
    std::optional<double> first_volume;
    for (int run = 0; run < 2; ++run) {
      const geom::Result<geom::CurvedBooleanResult> result = geom::boolean_solid_curved(
          sphere, only_solid(sphere), *entry.b, only_solid(*entry.b), entry.op,
          geom::CurvedBooleanOptions{tolerance});
      REQUIRE_MESSAGE(result.has_value(), entry.name << ": " << message_of(result.error()));
      const geom::Brep& brep = result.value().brep;
      CHECK(geom::validate_topology(brep).ok());
      const auto shells = brep.solid_shells(only_solid(brep));
      REQUIRE(shells.size() == 1);
      const double volume = volume_of(brep);
      CHECK(std::abs(volume - entry.expected) <= 1e-6 * std::abs(entry.expected));
      if (entry.mesh_closable) {
        check_mesh_closed(entry.name, brep);
      } else {
        CHECK_MESSAGE(
            geom::detail::heal_detail::shell_is_closed_by_uses(brep, shells.front()),
            entry.name << ": an edge is not used by exactly two coedges");
      }
      if (run == 0) {
        first_volume = volume;
      } else {
        CHECK(volume == *first_volume);
      }
    }
  }
}

// Task 1 of M3c-3b: the analytic sphere/sphere pair circle. Two unit spheres
// whose centres are 1.5 apart cross in the plane x = 0.75 with radius
// sqrt(1 - 0.75^2) = sqrt(7)/4, perpendicular to the centre line.
TEST_CASE("curved boolean: sphere/sphere has one analytic circle") {
  const geom::PlaneFrame no_frame{};
  const geom::Sphere3 sphere_a{geom::Point3d{0.0, 0.0, 0.0}, geom::Vec3d{0.0, 0.0, 1.0},
                               geom::Vec3d{1.0, 0.0, 0.0}, 1.0};
  const geom::Sphere3 sphere_b{geom::Point3d{1.5, 0.0, 0.0}, geom::Vec3d{0.0, 0.0, 1.0},
                               geom::Vec3d{1.0, 0.0, 0.0}, 1.0};
  const std::vector<geom::CircleCurve3> circles = geom::detail::curved_boolean_pair_circles(
      geom::SurfaceGeometry{sphere_a}, no_frame, geom::SurfaceGeometry{sphere_b}, no_frame,
      1e-9);
  REQUIRE(circles.size() == 1);
  CHECK(circles[0].center[0] == doctest::Approx(0.75));
  CHECK(circles[0].center[1] == doctest::Approx(0.0));
  CHECK(circles[0].center[2] == doctest::Approx(0.0));
  CHECK(circles[0].radius == doctest::Approx(std::sqrt(7.0) / 4.0));
  CHECK(circles[0].is_valid());

  // Tangent, contained and separate boundaries produce no cutting circle. The
  // containment case uses radii 1.0 and 0.2 so that a centre distance of 0.5
  // really is `d + rb <= ra`; with two unit spheres 0.5 apart is a lens, not
  // containment.
  const auto sphere_at = [](geom::Point3d center, double radius) {
    return geom::SurfaceGeometry{geom::Sphere3{center, geom::Vec3d{0.0, 0.0, 1.0},
                                               geom::Vec3d{1.0, 0.0, 0.0}, radius}};
  };
  const geom::SurfaceGeometry origin = sphere_at(geom::Point3d{0.0, 0.0, 0.0}, 1.0);
  CHECK(geom::detail::curved_boolean_pair_circles(
            origin, no_frame, sphere_at(geom::Point3d{2.0, 0.0, 0.0}, 1.0), no_frame, 1e-9)
            .empty());
  CHECK(geom::detail::curved_boolean_pair_circles(
            origin, no_frame, sphere_at(geom::Point3d{0.5, 0.0, 0.0}, 0.2), no_frame, 1e-9)
            .empty());
  CHECK(geom::detail::curved_boolean_pair_circles(
            origin, no_frame, sphere_at(geom::Point3d{2.5, 0.0, 0.0}, 1.0), no_frame, 1e-9)
            .empty());
  // Near external tangency: the d < ra + rb guard must fire before the radius
  // is computed.
  CHECK(geom::detail::curved_boolean_pair_circles(
            origin, no_frame, sphere_at(geom::Point3d{2.0 - 1e-10, 0.0, 0.0}, 1.0), no_frame,
            1e-9)
            .empty());
  // Containment boundary exactly d == |ra - rb| + tolerance: the guard is strict.
  CHECK(geom::detail::curved_boolean_pair_circles(
            origin, no_frame, sphere_at(geom::Point3d{0.75, 0.0, 0.0}, 0.5), no_frame, 0.25)
            .empty());
  // rho^2 <= tolerance^2 with a tolerance on the scale of the radii.
  const geom::SurfaceGeometry centered_101 = sphere_at(geom::Point3d{0.0, 0.0, 0.0}, 1.01);
  CHECK(geom::detail::curved_boolean_pair_circles(
            centered_101, no_frame, sphere_at(geom::Point3d{1.02 - 1e-7, 0.0, 0.0}, 1.01),
            no_frame, 1.0)
            .empty());
  // Unequal radii: the (ra^2 - rb^2) term moves the circle plane to x = 0.98.
  const std::vector<geom::CircleCurve3> unequal = geom::detail::curved_boolean_pair_circles(
      origin, no_frame, sphere_at(geom::Point3d{1.0, 0.0, 0.0}, 0.2), no_frame, 1e-9);
  REQUIRE(unequal.size() == 1);
  CHECK(unequal[0].center[0] == doctest::Approx(0.98));
  CHECK(unequal[0].radius == doctest::Approx(std::sqrt(0.0396)));
  // Frame handedness for n = (1, 0, 0): candidate is {0, 1, 0}, so the frame is
  // x_axis = {0, 0, 1}, y_axis = {0, -1, 0}.
  CHECK(circles[0].x_axis[2] == doctest::Approx(1.0));
  CHECK(circles[0].y_axis[1] == doctest::Approx(-1.0));
}

TEST_CASE("curved boolean: torus pairs need a plane cut") {
  const double tolerance = 1e-9;
  const geom::Result<geom::Brep> torus_a =
      geom::make_torus(geom::Point3d{0.0, 0.0, 0.0}, geom::Vec3d{0.0, 0.0, 1.0}, 1.0, 0.5,
                       tolerance);
  const geom::Result<geom::Brep> torus_b =
      geom::make_torus(geom::Point3d{0.3, 0.0, 0.0}, geom::Vec3d{0.0, 0.0, 1.0}, 1.0, 0.5,
                       tolerance);
  const geom::Result<geom::Brep> cylinder =
      geom::make_cylinder(geom::Point3d{0.0, 0.0, -2.0}, geom::Vec3d{0.0, 0.0, 1.0}, 1.2, 4.0,
                         tolerance);
  REQUIRE(torus_a.has_value());
  REQUIRE(torus_b.has_value());
  REQUIRE(cylinder.has_value());
  const std::string expected =
      "boolean_solid_curved: torus operands only support plane cuts; torus/torus and "
      "torus/curved pairs are deferred to M3c-3";
  for (const geom::SolidBooleanOp op : {geom::SolidBooleanOp::union_,
                                        geom::SolidBooleanOp::intersection,
                                        geom::SolidBooleanOp::difference,
                                        geom::SolidBooleanOp::symmetric_difference}) {
    for (const geom::Brep* other : {&torus_b.value(), &cylinder.value()}) {
      geom::CurvedBooleanReport report{};
      report.face_pairs = 999;
      const geom::Result<geom::CurvedBooleanResult> out = geom::boolean_solid_curved(
          torus_a.value(), only_solid(torus_a.value()), *other, only_solid(*other), op,
          geom::CurvedBooleanOptions{tolerance}, &report);
      REQUIRE_FALSE(out.has_value());
      CHECK(out.error().code == geom::Errc::unsupported);
      CHECK(std::string(out.error().message) == expected);
      CHECK(report.face_pairs == 0);
    }
  }
}

TEST_CASE("curved boolean: torus and half space") {
  const double tolerance = 1e-9;
  const geom::Result<geom::Brep> torus_result =
      geom::make_torus(geom::Point3d{0.0, 0.0, 0.0}, geom::Vec3d{0.0, 0.0, 1.0}, 1.0, 0.5,
                       tolerance);
  const geom::Result<geom::Brep> box_result =
      geom::make_box(geom::Point3d{-2.0, -2.0, 0.0}, geom::Point3d{2.0, 2.0, 2.0}, tolerance);
  REQUIRE(torus_result.has_value());
  REQUIRE(box_result.has_value());
  const geom::Brep& torus = torus_result.value();
  const geom::Brep& box = box_result.value();

  const geom::Result<geom::CurvedBooleanResult> half = geom::boolean_solid_curved(
      torus, only_solid(torus), box, only_solid(box), geom::SolidBooleanOp::intersection,
      geom::CurvedBooleanOptions{tolerance});
  REQUIRE_MESSAGE(half.has_value(),
                  (half.has_value() ? std::string_view{} : message_of(half.error())));
  CHECK(half.value().report.max_cut_residual <= 1e-3);
  CHECK(geom::validate_topology(half.value().brep).ok());
  CHECK(volume_of(half.value().brep) ==
        doctest::Approx(2.4674011002723395).epsilon(1e-6));

  const geom::Result<geom::CurvedBooleanResult> other_half = geom::boolean_solid_curved(
      torus, only_solid(torus), box, only_solid(box), geom::SolidBooleanOp::difference,
      geom::CurvedBooleanOptions{tolerance});
  REQUIRE_MESSAGE(other_half.has_value(),
                  (other_half.has_value() ? std::string_view{}
                                          : message_of(other_half.error())));
  CHECK(geom::validate_topology(other_half.value().brep).ok());
  CHECK(volume_of(other_half.value().brep) ==
        doctest::Approx(2.4674011002723395).epsilon(1e-6));
}

TEST_CASE("curved boolean: torus cut across the v seam") {
  const double tolerance = 1e-9;
  const geom::Result<geom::Brep> torus_result =
      geom::make_torus(geom::Point3d{0.0, 0.0, 0.0}, geom::Vec3d{0.0, 0.0, 1.0}, 1.0, 0.5,
                       tolerance);
  const geom::Result<geom::Brep> box_result =
      geom::make_box(geom::Point3d{0.0, -2.0, -2.0}, geom::Point3d{2.0, 2.0, 2.0}, tolerance);
  REQUIRE(torus_result.has_value());
  REQUIRE(box_result.has_value());
  const geom::Brep& torus = torus_result.value();
  const geom::Brep& box = box_result.value();

  const geom::Result<geom::CurvedBooleanResult> half = geom::boolean_solid_curved(
      torus, only_solid(torus), box, only_solid(box), geom::SolidBooleanOp::intersection,
      geom::CurvedBooleanOptions{tolerance});
  REQUIRE_MESSAGE(half.has_value(),
                  (half.has_value() ? std::string_view{} : message_of(half.error())));
  CHECK(geom::validate_topology(half.value().brep).ok());
  CHECK(volume_of(half.value().brep) ==
        doctest::Approx(2.4674011002723395).epsilon(1e-6));
  CHECK(half.value().report.max_cut_residual <= 1e-3);
}

TEST_CASE("curved boolean: torus minus half space") {
  const double tolerance = 1e-9;
  const geom::Result<geom::Brep> torus_result =
      geom::make_torus(geom::Point3d{0.0, 0.0, 0.0}, geom::Vec3d{0.0, 0.0, 1.0}, 1.0, 0.5,
                       tolerance);
  const geom::Result<geom::Brep> box_result =
      geom::make_box(geom::Point3d{-2.0, -2.0, 0.0}, geom::Point3d{2.0, 2.0, 2.0}, tolerance);
  REQUIRE(torus_result.has_value());
  REQUIRE(box_result.has_value());
  const geom::Brep& torus = torus_result.value();
  const geom::Brep& box = box_result.value();

  const geom::Result<geom::CurvedBooleanResult> half = geom::boolean_solid_curved(
      torus, only_solid(torus), box, only_solid(box), geom::SolidBooleanOp::difference,
      geom::CurvedBooleanOptions{tolerance});
  REQUIRE_MESSAGE(half.has_value(),
                  (half.has_value() ? std::string_view{} : message_of(half.error())));
  CHECK(geom::validate_topology(half.value().brep).ok());
}

TEST_CASE("curved boolean: a box that contains a torus keeps it whole") {
  const double tolerance = 1e-9;
  const geom::Result<geom::Brep> torus_result =
      geom::make_torus(geom::Point3d{0.0, 0.0, 0.0}, geom::Vec3d{0.0, 0.0, 1.0}, 1.0, 0.5,
                       tolerance);
  const geom::Result<geom::Brep> box_result =
      geom::make_box(geom::Point3d{-2.0, -2.0, -2.0}, geom::Point3d{2.0, 2.0, 2.0}, tolerance);
  REQUIRE(torus_result.has_value());
  REQUIRE(box_result.has_value());
  const geom::Brep& torus = torus_result.value();
  const geom::Brep& box = box_result.value();

  const geom::Result<geom::CurvedBooleanResult> whole = geom::boolean_solid_curved(
      torus, only_solid(torus), box, only_solid(box), geom::SolidBooleanOp::intersection,
      geom::CurvedBooleanOptions{tolerance});
  REQUIRE_MESSAGE(whole.has_value(),
                  (whole.has_value() ? std::string_view{} : message_of(whole.error())));
  CHECK(geom::validate_topology(whole.value().brep).ok());
  CHECK(volume_of(whole.value().brep) ==
        doctest::Approx(4.934802200544679).epsilon(1e-6));
}

// A closed slab whose bottom face is the tilted plane z = 0.3 x. The shear keeps
// every face planar and the shell closed, so the driver can reach its torus cut
// plane check instead of failing on the operand itself. Only the bottom face
// overlaps the torus, which makes the rejection deterministic.
geom::Brep make_tilted_slab(double tolerance) {
  const double slope = 0.3;
  const geom::Point3d corners[8] = {
      {-4.0, -4.0, 0.0 + slope * -4.0}, {4.0, -4.0, 0.0 + slope * 4.0},
      {4.0, 4.0, 0.0 + slope * 4.0},    {-4.0, 4.0, 0.0 + slope * -4.0},
      {-4.0, -4.0, 4.0 + slope * -4.0}, {4.0, -4.0, 4.0 + slope * 4.0},
      {4.0, 4.0, 4.0 + slope * 4.0},    {-4.0, 4.0, 4.0 + slope * -4.0},
  };
  geom::SolidEditor editor;
  geom::VertexId vertices[8];
  for (int i = 0; i < 8; ++i) {
    vertices[i] = editor.add_vertex(corners[i], tolerance).value();
  }
  const int edge_specs[12][2] = {{0, 1}, {1, 2}, {2, 3},  {3, 0},  {4, 5}, {5, 6},
                                 {6, 7}, {7, 4}, {1, 5},  {2, 6},  {3, 7}, {0, 4}};
  geom::EdgeId edges[12];
  for (int i = 0; i < 12; ++i) {
    const geom::Point3d from = corners[edge_specs[i][0]];
    const geom::Point3d to = corners[edge_specs[i][1]];
    const geom::Vec3d delta = to - from;
    const double length = std::hypot(delta[0], delta[1], delta[2]);
    const geom::Vec3d direction{delta[0] / length, delta[1] / length, delta[2] / length};
    const std::uint32_t curve = editor.add_curve(geom::LineCurve3{from, direction});
    edges[i] = editor.add_edge(vertices[edge_specs[i][0]], vertices[edge_specs[i][1]], curve,
                               0.0, length, tolerance)
                   .value();
  }
  const double side = std::sqrt(0.09 + 1.0);
  const geom::Vec3d bottom_normal{0.3 / side, 0.0, -1.0 / side};
  const geom::Vec3d top_normal{-0.3 / side, 0.0, 1.0 / side};
  struct FaceSpec {
    geom::Vec3d normal;
    int corner;
    int loop[4];
    bool reversed[4];
  };
  const FaceSpec face_specs[6] = {
      {bottom_normal, 0, {3, 2, 1, 0}, {true, true, true, true}},
      {top_normal, 4, {4, 5, 6, 7}, {false, false, false, false}},
      {{0.0, -1.0, 0.0}, 0, {0, 8, 4, 11}, {false, false, true, true}},
      {{0.0, 1.0, 0.0}, 2, {2, 10, 6, 9}, {false, false, true, true}},
      {{-1.0, 0.0, 0.0}, 0, {11, 7, 10, 3}, {false, true, true, false}},
      {{1.0, 0.0, 0.0}, 1, {1, 9, 5, 8}, {false, false, true, true}},
  };
  std::array<geom::FaceId, 6> faces{};
  int face_index = 0;
  for (const FaceSpec& spec : face_specs) {
    const geom::Point3d corner = corners[spec.corner];
    const std::uint32_t surface = editor.add_surface(geom::Plane3{
        spec.normal,
        geom::dot(spec.normal, geom::Vec3d{corner[0], corner[1], corner[2]})});
    const geom::FaceId face = editor.add_face(surface, false, tolerance).value();
    std::array<geom::CoedgeSpec, 4> coedges{};
    for (int i = 0; i < 4; ++i) {
      coedges[static_cast<std::size_t>(i)] = geom::CoedgeSpec{
          edges[spec.loop[i]], spec.reversed[i]};
    }
    editor.add_loop(face, coedges, true).value();
    faces[static_cast<std::size_t>(face_index++)] = face;
  }
  const geom::ShellId shell = editor.add_shell(faces, true).value();
  const std::array<geom::ShellId, 1> shells{shell};
  editor.add_solid(shells, tolerance).value();
  return editor.commit().value();
}

TEST_CASE("curved boolean: torus cut by an offset parallel plane is rejected") {
  const double tolerance = 1e-9;
  const geom::Result<geom::Brep> torus =
      geom::make_torus(geom::Point3d{0.0, 0.0, 0.0}, geom::Vec3d{0.0, 0.0, 1.0}, 1.0, 0.5,
                       tolerance);
  const geom::Result<geom::Brep> box =
      geom::make_box(geom::Point3d{0.2, -2.0, -2.0}, geom::Point3d{2.0, 2.0, 2.0}, tolerance);
  REQUIRE(torus.has_value());
  REQUIRE(box.has_value());
  geom::CurvedBooleanReport report;
  const geom::Result<geom::CurvedBooleanResult> result = geom::boolean_solid_curved(
      torus.value(), only_solid(torus.value()), box.value(), only_solid(box.value()),
      geom::SolidBooleanOp::intersection, geom::CurvedBooleanOptions{tolerance}, &report);
  REQUIRE_FALSE(result.has_value());
  CHECK(result.error().code == geom::Errc::unsupported);
  CHECK(message_of(result.error()) ==
        "boolean_solid_curved: torus cut plane must be perpendicular to or contain the torus "
        "axis");
  CHECK(report.face_pairs == 0);
}

TEST_CASE("curved boolean: torus cut by a tilted plane is rejected") {
  const double tolerance = 1e-9;
  const geom::Result<geom::Brep> torus =
      geom::make_torus(geom::Point3d{0.0, 0.0, 0.0}, geom::Vec3d{0.0, 0.0, 1.0}, 1.0, 0.5,
                       tolerance);
  REQUIRE(torus.has_value());
  const geom::Brep slab = make_tilted_slab(tolerance);
  REQUIRE(geom::validate_topology(slab).ok());
  geom::CurvedBooleanReport report;
  const geom::Result<geom::CurvedBooleanResult> result = geom::boolean_solid_curved(
      torus.value(), only_solid(torus.value()), slab, only_solid(slab),
      geom::SolidBooleanOp::intersection, geom::CurvedBooleanOptions{tolerance}, &report);
  REQUIRE_FALSE(result.has_value());
  CHECK(result.error().code == geom::Errc::unsupported);
  CHECK(message_of(result.error()) ==
        "boolean_solid_curved: torus cut plane must be perpendicular to or contain the torus "
        "axis");
  CHECK(report.face_pairs == 0);
}

// A closed slab whose bottom face is the plane through the origin with the given
// normal, extending four units along that normal. Only the bottom face reaches
// the torus bounds (|x|,|y|,|z| <= 1.5), so one torus/plane pair decides the
// cut-plane guard. The corner ring and coedge pattern follow make_tilted_slab;
// each face plane is derived from the corner sequence its loop traverses.
geom::Brep make_plane_slab(const geom::Vec3d& normal, double tolerance) {
  const double length = std::sqrt(geom::dot(normal, normal));
  REQUIRE(length > 0.0);
  const geom::Vec3d n{normal[0] / length, normal[1] / length, normal[2] / length};
  const geom::Vec3d seed =
      std::fabs(n[0]) < 0.9 ? geom::Vec3d{1.0, 0.0, 0.0} : geom::Vec3d{0.0, 1.0, 0.0};
  const geom::Vec3d raw_u = geom::cross(seed, n);
  const double u_length = std::sqrt(geom::dot(raw_u, raw_u));
  REQUIRE(u_length > 0.0);
  const geom::Vec3d u{raw_u[0] / u_length, raw_u[1] / u_length, raw_u[2] / u_length};
  const geom::Vec3d v = geom::cross(n, u);
  const auto corner = [&](double su, double sv, double sn) {
    return geom::Point3d{su * u[0] + sv * v[0] + sn * n[0],
                         su * u[1] + sv * v[1] + sn * n[1],
                         su * u[2] + sv * v[2] + sn * n[2]};
  };
  const geom::Point3d corners[8] = {
      corner(-4.0, -4.0, 0.0), corner(4.0, -4.0, 0.0), corner(4.0, 4.0, 0.0),
      corner(-4.0, 4.0, 0.0),  corner(-4.0, -4.0, 4.0), corner(4.0, -4.0, 4.0),
      corner(4.0, 4.0, 4.0),   corner(-4.0, 4.0, 4.0),
  };
  geom::SolidEditor editor;
  geom::VertexId vertices[8];
  for (int i = 0; i < 8; ++i) {
    vertices[i] = editor.add_vertex(corners[i], tolerance).value();
  }
  const int edge_specs[12][2] = {{0, 1}, {1, 2}, {2, 3},  {3, 0},  {4, 5}, {5, 6},
                                 {6, 7}, {7, 4}, {1, 5},  {2, 6},  {3, 7}, {0, 4}};
  geom::EdgeId edges[12];
  for (int i = 0; i < 12; ++i) {
    const geom::Point3d from = corners[edge_specs[i][0]];
    const geom::Point3d to = corners[edge_specs[i][1]];
    const geom::Vec3d delta = to - from;
    const double edge_length = std::hypot(delta[0], delta[1], delta[2]);
    const geom::Vec3d direction{delta[0] / edge_length, delta[1] / edge_length,
                                delta[2] / edge_length};
    const std::uint32_t curve = editor.add_curve(geom::LineCurve3{from, direction});
    edges[i] = editor.add_edge(vertices[edge_specs[i][0]], vertices[edge_specs[i][1]], curve,
                               0.0, edge_length, tolerance)
                   .value();
  }
  // Edges traversed by each loop, in loop order, with the coedge direction flags.
  const int loop_edges[6][4] = {{3, 2, 1, 0}, {4, 5, 6, 7}, {0, 8, 4, 11},
                                {2, 10, 6, 9}, {11, 7, 10, 3}, {1, 9, 5, 8}};
  const bool loop_reversed[6][4] = {{true, true, true, true},   {false, false, false, false},
                                    {false, false, true, true}, {false, false, true, true},
                                    {false, true, true, false}, {false, false, true, true}};
  // Corners the same loops visit, used to derive each face plane exactly.
  const int loop_corners[6][4] = {{0, 3, 2, 1}, {4, 5, 6, 7}, {0, 1, 5, 4},
                                  {2, 3, 7, 6}, {0, 4, 7, 3}, {1, 2, 6, 5}};
  std::array<geom::FaceId, 6> faces{};
  for (int f = 0; f < 6; ++f) {
    const geom::Point3d a = corners[loop_corners[f][0]];
    const geom::Point3d b = corners[loop_corners[f][1]];
    const geom::Point3d c = corners[loop_corners[f][2]];
    const geom::Vec3d ab{b[0] - a[0], b[1] - a[1], b[2] - a[2]};
    const geom::Vec3d ac{c[0] - a[0], c[1] - a[1], c[2] - a[2]};
    const geom::Vec3d raw = geom::cross(ab, ac);
    const double normal_length = std::sqrt(geom::dot(raw, raw));
    REQUIRE(normal_length > 0.0);
    const geom::Vec3d face_normal{raw[0] / normal_length, raw[1] / normal_length,
                                  raw[2] / normal_length};
    const std::uint32_t surface = editor.add_surface(geom::Plane3{
        face_normal, geom::dot(face_normal, geom::Vec3d{a[0], a[1], a[2]})});
    const geom::FaceId face = editor.add_face(surface, false, tolerance).value();
    std::array<geom::CoedgeSpec, 4> coedges{};
    for (int i = 0; i < 4; ++i) {
      coedges[static_cast<std::size_t>(i)] =
          geom::CoedgeSpec{edges[loop_edges[f][i]], loop_reversed[f][i]};
    }
    editor.add_loop(face, coedges, true).value();
    faces[static_cast<std::size_t>(f)] = face;
  }
  const geom::ShellId shell = editor.add_shell(faces, true).value();
  const std::array<geom::ShellId, 1> shells{shell};
  editor.add_solid(shells, tolerance).value();
  return editor.commit().value();
}

TEST_CASE("curved boolean: torus cut by a near-axis plane is rejected quickly") {
  const double tolerance = 1e-9;
  // The cutting plane normal is (1, 0, 1e-10): it contains the torus axis only up
  // to rounding. Before the 1e-12 guard this pair was accepted and the exact cut
  // path ran for a very long time; the rejection must be immediate.
  const geom::Result<geom::Brep> torus =
      geom::make_torus(geom::Point3d{0.0, 0.0, 0.0}, geom::Vec3d{0.0, 0.0, 1.0}, 1.0, 0.5,
                       tolerance);
  REQUIRE(torus.has_value());
  const geom::Brep slab = make_plane_slab(geom::Vec3d{1.0, 0.0, 1e-10}, tolerance);
  REQUIRE(geom::validate_topology(slab).ok());
  geom::CurvedBooleanReport report;
  const std::chrono::steady_clock::time_point start = std::chrono::steady_clock::now();
  const geom::Result<geom::CurvedBooleanResult> result = geom::boolean_solid_curved(
      torus.value(), only_solid(torus.value()), slab, only_solid(slab),
      geom::SolidBooleanOp::intersection, geom::CurvedBooleanOptions{tolerance}, &report);
  const double elapsed =
      std::chrono::duration<double>(std::chrono::steady_clock::now() - start).count();
  REQUIRE_FALSE(result.has_value());
  CHECK(result.error().code == geom::Errc::unsupported);
  CHECK(message_of(result.error()) ==
        "boolean_solid_curved: torus cut plane must be perpendicular to or contain the torus "
        "axis");
  CHECK(report.face_pairs == 0);
  CHECK(elapsed < 5.0);
}

TEST_CASE("curved boolean: torus cut by a near-perpendicular plane is rejected") {
  const double tolerance = 1e-9;
  // The axis is 5e-10 short of the exact box face normal. The exact cut path
  // cannot represent that loop (trim_face reports "segment intersection is not
  // representable"), so the guard must reject it instead of letting it through.
  const double axial = 1.0 - 5e-10;
  const double radial = std::sqrt(1.0 - axial * axial);
  const geom::Result<geom::Brep> torus =
      geom::make_torus(geom::Point3d{0.0, 0.0, 0.0}, geom::Vec3d{radial, 0.0, axial}, 1.0, 0.5,
                       tolerance);
  const geom::Result<geom::Brep> box =
      geom::make_box(geom::Point3d{-2.0, -2.0, 0.0}, geom::Point3d{2.0, 2.0, 2.0}, tolerance);
  REQUIRE(torus.has_value());
  REQUIRE(box.has_value());
  geom::CurvedBooleanReport report;
  const geom::Result<geom::CurvedBooleanResult> result = geom::boolean_solid_curved(
      torus.value(), only_solid(torus.value()), box.value(), only_solid(box.value()),
      geom::SolidBooleanOp::intersection, geom::CurvedBooleanOptions{tolerance}, &report);
  REQUIRE_FALSE(result.has_value());
  CHECK(result.error().code == geom::Errc::unsupported);
  CHECK(message_of(result.error()) ==
        "boolean_solid_curved: torus cut plane must be perpendicular to or contain the torus "
        "axis");
  CHECK(report.face_pairs == 0);
}

TEST_CASE("curved boolean: torus cut by a 5e-13 short axis is rejected quickly") {
  const double tolerance = 1e-9;
  // This is the accepted side of the old 1e-12 band: 5e-13 short of perpendicular.
  // No snapping exists, so the exact cut path still dies in trim_face; the guard
  // must reject the pair immediately (this case took the full suite down at
  // 8cc43af with "segment intersection is not representable").
  const double axial = 1.0 - 5e-13;
  const double radial = std::sqrt(1.0 - axial * axial);
  const geom::Result<geom::Brep> torus =
      geom::make_torus(geom::Point3d{0.0, 0.0, 0.0}, geom::Vec3d{radial, 0.0, axial}, 1.0, 0.5,
                       tolerance);
  const geom::Result<geom::Brep> box =
      geom::make_box(geom::Point3d{-2.0, -2.0, 0.0}, geom::Point3d{2.0, 2.0, 2.0}, tolerance);
  REQUIRE(torus.has_value());
  REQUIRE(box.has_value());
  geom::CurvedBooleanReport report;
  const std::chrono::steady_clock::time_point start = std::chrono::steady_clock::now();
  const geom::Result<geom::CurvedBooleanResult> result = geom::boolean_solid_curved(
      torus.value(), only_solid(torus.value()), box.value(), only_solid(box.value()),
      geom::SolidBooleanOp::intersection, geom::CurvedBooleanOptions{tolerance}, &report);
  const double elapsed =
      std::chrono::duration<double>(std::chrono::steady_clock::now() - start).count();
  REQUIRE_FALSE(result.has_value());
  CHECK(result.error().code == geom::Errc::unsupported);
  CHECK(message_of(result.error()) ==
        "boolean_solid_curved: torus cut plane must be perpendicular to or contain the torus "
        "axis");
  CHECK(report.face_pairs == 0);
  CHECK(elapsed < 5.0);
}

TEST_CASE("curved boolean: torus cut by a 1e-12 short axis is rejected") {
  const double tolerance = 1e-9;
  // Boundary of the old 1e-12 band: 1e-12 short of perpendicular was accepted
  // by the 8cc43af guard. The normalized-axis 8-ulp band rejects it as well.
  const double axial = 1.0 - 1e-12;
  const double radial = std::sqrt(1.0 - axial * axial);
  const geom::Result<geom::Brep> torus =
      geom::make_torus(geom::Point3d{0.0, 0.0, 0.0}, geom::Vec3d{radial, 0.0, axial}, 1.0, 0.5,
                       tolerance);
  const geom::Result<geom::Brep> box =
      geom::make_box(geom::Point3d{-2.0, -2.0, 0.0}, geom::Point3d{2.0, 2.0, 2.0}, tolerance);
  REQUIRE(torus.has_value());
  REQUIRE(box.has_value());
  geom::CurvedBooleanReport report;
  const geom::Result<geom::CurvedBooleanResult> result = geom::boolean_solid_curved(
      torus.value(), only_solid(torus.value()), box.value(), only_solid(box.value()),
      geom::SolidBooleanOp::intersection, geom::CurvedBooleanOptions{tolerance}, &report);
  REQUIRE_FALSE(result.has_value());
  CHECK(result.error().code == geom::Errc::unsupported);
  CHECK(message_of(result.error()) ==
        "boolean_solid_curved: torus cut plane must be perpendicular to or contain the torus "
        "axis");
  CHECK(report.face_pairs == 0);
}

TEST_CASE("curved boolean: rounded unit axis still cuts the torus in half") {
  const double tolerance = 1e-9;
  // make_torus stores the axis unnormalized and only requires it within
  // kAnalyticSurfaceUnitTolerance (1e-12), and even a unit axis can give a
  // normalized dot of 1 +/- 1 ulp. This theta (hypot exactly 1.0) is the sampled
  // case where the dot rounds below 1.0, so exact equality rejected a valid cut.
  const double angle = 0.00013823007675795089;
  const geom::Vec3d axis{std::cos(angle), std::sin(angle), 0.0};
  REQUIRE(std::hypot(axis[0], axis[1]) == 1.0);
  const geom::Result<geom::Brep> torus =
      geom::make_torus(geom::Point3d{0.0, 0.0, 0.0}, axis, 1.0, 0.5, tolerance);
  REQUIRE(torus.has_value());
  // The slab's bottom face is the plane through the origin with the same normal,
  // so it cuts the torus into two equal halves.
  const geom::Brep slab = make_plane_slab(axis, tolerance);
  REQUIRE(geom::validate_topology(slab).ok());
  geom::CurvedBooleanReport report;
  const geom::Result<geom::CurvedBooleanResult> half = geom::boolean_solid_curved(
      torus.value(), only_solid(torus.value()), slab, only_solid(slab),
      geom::SolidBooleanOp::intersection, geom::CurvedBooleanOptions{tolerance}, &report);
  REQUIRE_MESSAGE(half.has_value(),
                  (half.has_value() ? std::string_view{} : message_of(half.error())));
  CHECK(geom::validate_topology(half.value().brep).ok());
  CHECK(volume_of(half.value().brep) ==
        doctest::Approx(2.4674011002723395).epsilon(1e-6));
  CHECK(half.value().report.face_pairs >= 1);
}

TEST_CASE("curved boolean: a torus axis with a tiny norm excess still cuts") {
  const double tolerance = 1e-9;
  // make_torus accepts an axis only within 1e-12 of unit length but keeps it
  // unnormalized, so the alignment dot carries the norm excess. The guard must
  // normalize before comparing (5e-13 is far outside the 8 ulp band).
  const geom::Result<geom::Brep> torus = geom::make_torus(
      geom::Point3d{0.0, 0.0, 0.0}, geom::Vec3d{0.0, 0.0, 1.0 + 5e-13}, 1.0, 0.5, tolerance);
  REQUIRE(torus.has_value());
  const geom::Result<geom::Brep> box =
      geom::make_box(geom::Point3d{-2.0, -2.0, 0.0}, geom::Point3d{2.0, 2.0, 2.0}, tolerance);
  REQUIRE(box.has_value());
  geom::CurvedBooleanReport report;
  const geom::Result<geom::CurvedBooleanResult> half = geom::boolean_solid_curved(
      torus.value(), only_solid(torus.value()), box.value(), only_solid(box.value()),
      geom::SolidBooleanOp::intersection, geom::CurvedBooleanOptions{tolerance}, &report);
  REQUIRE_MESSAGE(half.has_value(),
                  (half.has_value() ? std::string_view{} : message_of(half.error())));
  CHECK(geom::validate_topology(half.value().brep).ok());
  CHECK(volume_of(half.value().brep) ==
        doctest::Approx(2.4674011002723395).epsilon(1e-6));
  CHECK(half.value().report.face_pairs >= 1);
}

TEST_CASE("brep: boolean_solid_curved cuts a full cylinder with a box above its equator") {
  const double tolerance = 1e-9;
  const geom::Brep cylinder = make_test_cylinder(tolerance);
  REQUIRE(geom::validate_topology(cylinder).ok());
  // z >= 1 leaves the lower half (pi) inside the box and nothing else: the box
  // volume is 4 * 4 * 9 = 144.
  const geom::Brep box =
      geom::make_box(geom::Point3d{-2.0, -2.0, 1.0}, geom::Point3d{2.0, 2.0, 10.0}, tolerance)
          .value();
  REQUIRE(geom::validate_topology(box).ok());

  struct Case {
    geom::SolidBooleanOp op;
    double volume;
    std::size_t faces;
  };
  const Case cases[] = {
      {geom::SolidBooleanOp::union_, 144.0 + kPi, 8},
      {geom::SolidBooleanOp::intersection, kPi, 3},
      {geom::SolidBooleanOp::difference, kPi, 3},
  };
  for (const Case& item : cases) {
    const geom::Result<geom::CurvedBooleanResult> result = geom::boolean_solid_curved(
        cylinder, only_solid(cylinder), box, only_solid(box), item.op,
        geom::CurvedBooleanOptions{tolerance});
    REQUIRE_MESSAGE(result.has_value(), "boolean_solid_curved failed: "
                                            << message_of(result.error()));
    const geom::Brep& brep = result.value().brep;
    CHECK(geom::validate_topology(brep).ok());
    CHECK(std::abs(volume_of(brep) - item.volume) <= 1e-6);
    CHECK(brep.face_ids().size() == item.faces);
    CHECK(result.value().report.kept ==
          static_cast<std::uint32_t>(brep.face_ids().size()));
  }
}

TEST_CASE("brep: boolean_solid_curved cuts a full cylinder with a vertical slab") {
  const double tolerance = 1e-9;
  const geom::Brep cylinder = make_test_cylinder(tolerance);
  const geom::Brep slab =
      geom::make_box(geom::Point3d{0.5, -2.0, -1.0}, geom::Point3d{2.0, 2.0, 10.0}, tolerance)
          .value();
  REQUIRE(geom::validate_topology(cylinder).ok());
  REQUIRE(geom::validate_topology(slab).ok());

  // The slab keeps x >= 0.5. One cap's area outside the slab is the circular
  // segment (pi / 3 - sqrt(3) / 4); the cylinder is 2 tall.
  const double segment = kPi / 3.0 - std::numbers::sqrt3_v<double> / 4.0;
  const double intersection = 2.0 * segment;
  const double difference = 2.0 * kPi - intersection;
  const double slab_volume = 1.5 * 4.0 * 11.0;

  const geom::Result<geom::CurvedBooleanResult> intersect = geom::boolean_solid_curved(
      cylinder, only_solid(cylinder), slab, only_solid(slab), geom::SolidBooleanOp::intersection,
      geom::CurvedBooleanOptions{tolerance});
  REQUIRE_MESSAGE(intersect.has_value(),
                  "boolean_solid_curved failed: " << message_of(intersect.error()));
  CHECK(geom::validate_topology(intersect.value().brep).ok());
  CHECK(std::abs(volume_of(intersect.value().brep) - intersection) <= 1e-6);
  CHECK(intersect.value().brep.face_ids().size() == 5);

  const geom::Result<geom::CurvedBooleanResult> take_away = geom::boolean_solid_curved(
      cylinder, only_solid(cylinder), slab, only_solid(slab), geom::SolidBooleanOp::difference,
      geom::CurvedBooleanOptions{tolerance});
  REQUIRE_MESSAGE(take_away.has_value(),
                  "boolean_solid_curved failed: " << message_of(take_away.error()));
  CHECK(geom::validate_topology(take_away.value().brep).ok());
  CHECK(std::abs(volume_of(take_away.value().brep) - difference) <= 1e-6);
  CHECK(take_away.value().brep.face_ids().size() == 4);

  const geom::Result<geom::CurvedBooleanResult> united = geom::boolean_solid_curved(
      cylinder, only_solid(cylinder), slab, only_solid(slab), geom::SolidBooleanOp::union_,
      geom::CurvedBooleanOptions{tolerance});
  REQUIRE_MESSAGE(united.has_value(),
                  "boolean_solid_curved failed: " << message_of(united.error()));
  CHECK(geom::validate_topology(united.value().brep).ok());
  CHECK(std::abs(volume_of(united.value().brep) - (slab_volume + difference)) <= 1e-6);
  CHECK(united.value().brep.face_ids().size() == 9);
}

TEST_CASE("brep: boolean_solid_curved defers a symmetric difference that touches itself") {
  const double tolerance = 1e-9;
  const geom::Brep a =
      geom::make_box(geom::Point3d{0.0, 0.0, 0.0}, geom::Point3d{2.0, 2.0, 2.0}, tolerance)
          .value();
  const geom::Brep b =
      geom::make_box(geom::Point3d{1.0, 1.0, 1.0}, geom::Point3d{3.0, 3.0, 3.0}, tolerance)
          .value();
  REQUIRE(geom::validate_topology(a).ok());
  REQUIRE(geom::validate_topology(b).ok());

  // Both remnants keep the shared cut faces, so the result is non-manifold along
  // their common edges: the assembler rejects it with invalid_topology, which the
  // driver reports as deferred instead of leaking the internal message.
  const geom::Result<geom::CurvedBooleanResult> result = geom::boolean_solid_curved(
      a, only_solid(a), b, only_solid(b), geom::SolidBooleanOp::symmetric_difference,
      geom::CurvedBooleanOptions{tolerance});
  REQUIRE(!result.has_value());
  CHECK(result.error().code == geom::Errc::unsupported);
  CHECK(message_of(result.error()) ==
        "boolean_solid_curved: result touches itself along an edge");
}

TEST_CASE("brep: boolean_solid_curved classification mesh is watertight") {
  const double tolerance = 1e-9;
  const geom::Brep quarter = make_exact_quarter_cylinder(1.0, 2.0, tolerance);
  REQUIRE(geom::validate_topology(quarter).ok());
  // tessellate only accepts trimmed faces that are full primitives, so the
  // quarter fixture with its open rim arcs exercises the uv fallback mesh.
  const geom::Result<geom::TriMesh> tessellated =
      geom::tessellate(quarter, geom::TessellationOptions{});
  CHECK(!tessellated.has_value());

  const auto mesh_of = [&tolerance](const geom::Brep& brep) {
    std::vector<geom::detail::CurvedBooleanFaceInfo> infos;
    for (const geom::FaceId id : brep.face_ids()) {
      const geom::Result<geom::detail::CurvedBooleanFaceInfo> info =
          geom::detail::curved_boolean_face_info(brep, id, tolerance);
      REQUIRE(info.has_value());
      infos.push_back(info.value());
    }
    const geom::Result<geom::TriMesh> mesh =
        geom::detail::curved_boolean_solid_mesh(brep, infos, tolerance);
    REQUIRE_MESSAGE(mesh.has_value(), "classification mesh failed: "
                                          << message_of(mesh.error()));
    return mesh.value();
  };

  const geom::detail::CurvedBooleanMeshClosure fallback =
      geom::detail::curved_boolean_mesh_closure(mesh_of(quarter));
  CHECK(fallback.edges > 0);
  CHECK(fallback.used_once == 0);
  CHECK(fallback.nonmanifold == 0);
  CHECK(fallback.closed());

  // A full primitive goes through tessellate instead; its mesh must be closed
  // too, otherwise the driver would fall back to a mesh it cannot trust.
  const geom::detail::CurvedBooleanMeshClosure full =
      geom::detail::curved_boolean_mesh_closure(mesh_of(make_test_cylinder(tolerance)));
  CHECK(full.closed());
}

TEST_CASE("brep: boolean_solid_curved cuts an offset cap box from a full cylinder") {
  const double tolerance = 1e-9;
  const geom::Brep cylinder = make_test_cylinder(tolerance);
  // The box removes the slab x in [0.5, 2], y in [-0.5, 0.5], z in [0.5, 1.5]
  // from above the z = 0.5 plane. The u = +-60 degree vertical lines (x = 0.5)
  // lie outside |y| <= 0.5 and are not cuts; the slab sides meet the lateral
  // face along the vertical lines u = +-30 degrees (y = +-0.5), so the cut
  // corners are the exact circle points at u = +-30 degrees, z = 0.5.
  const geom::Brep box =
      geom::make_box(geom::Point3d{0.5, -0.5, 0.5}, geom::Point3d{2.0, 0.5, 1.5}, tolerance)
          .value();
  REQUIRE(geom::validate_topology(cylinder).ok());
  REQUIRE(geom::validate_topology(box).ok());

  // Shared area of one cap plane and the unit disk {x >= 0.5, |y| <= 0.5}:
  //   int_{0.5}^{sqrt(3)/2} 1 dx + int_{sqrt(3)/2}^{1} 2 sqrt(1 - x^2) dx
  // = sqrt(3)/4 - 1/2 + pi/6. The overlap is one unit tall and its volume is
  // the same number.
  const double intersection = kPi / 6.0 + std::numbers::sqrt3_v<double> / 4.0 - 0.5;
  const double difference = 2.0 * kPi - intersection;
  const double box_volume = 1.5 * 1.0 * 1.0;

  struct Case {
    geom::SolidBooleanOp op;
    double volume;
    std::size_t faces;
  };
  const Case cases[] = {
      {geom::SolidBooleanOp::union_, box_volume + difference, 8},
      {geom::SolidBooleanOp::intersection, intersection, 7},
      {geom::SolidBooleanOp::difference, difference, 8},
  };
  for (const Case& item : cases) {
    const geom::Result<geom::CurvedBooleanResult> result = geom::boolean_solid_curved(
        cylinder, only_solid(cylinder), box, only_solid(box), item.op,
        geom::CurvedBooleanOptions{tolerance});
    REQUIRE_MESSAGE(result.has_value(), "boolean_solid_curved failed: "
                                            << message_of(result.error()));
    const geom::Brep& brep = result.value().brep;
    CHECK(geom::validate_topology(brep).ok());
    CHECK(std::abs(volume_of(brep) - item.volume) <= 1e-6);
    CHECK(brep.face_ids().size() == item.faces);
  }
}

TEST_CASE("brep: boolean_solid_curved cuts a box across both cylinder caps") {
  const double tolerance = 1e-9;
  const geom::Brep cylinder = make_test_cylinder(tolerance);
  // x in [-0.2, 0.4] reaches through both caps; the box is 1.5 tall from
  // z = 1.5 to 3, so the overlap is the 0.5 high slab.
  const geom::Brep box =
      geom::make_box(geom::Point3d{-0.2, -2.0, 1.5}, geom::Point3d{0.4, 2.0, 3.0}, tolerance)
          .value();
  REQUIRE(geom::validate_topology(cylinder).ok());
  REQUIRE(geom::validate_topology(box).ok());

  // Area of the disk strip between x = -0.2 and x = 0.4 is
  //   int_{-0.2}^{0.4} 2 sqrt(1 - x^2) dx
  // = [x sqrt(1 - x^2) + asin(x)]_{-0.2}^{0.4}, and the 0.5 height cancels the
  // two, so the intersection volume equals
  //   (0.4 sqrt(0.84) + asin(0.4) + 0.2 sqrt(0.96) + asin(0.2)) / 2.
  const auto cap = [](double x) {
    return (x * std::sqrt(1.0 - x * x) + std::asin(x)) / 2.0;
  };
  const double intersection = cap(0.4) - cap(-0.2);
  const double difference = 2.0 * kPi - intersection;
  const double box_volume = 0.6 * 4.0 * 1.5;

  struct Case {
    geom::SolidBooleanOp op;
    double volume;
    std::size_t faces;
  };
  const Case cases[] = {
      {geom::SolidBooleanOp::union_, box_volume + difference, 11},
      {geom::SolidBooleanOp::intersection, intersection, 6},
      {geom::SolidBooleanOp::difference, difference, 7},
  };
  for (const Case& item : cases) {
    const geom::Result<geom::CurvedBooleanResult> result = geom::boolean_solid_curved(
        cylinder, only_solid(cylinder), box, only_solid(box), item.op,
        geom::CurvedBooleanOptions{tolerance});
    REQUIRE_MESSAGE(result.has_value(), "boolean_solid_curved failed: "
                                            << message_of(result.error()));
    const geom::Brep& brep = result.value().brep;
    CHECK(geom::validate_topology(brep).ok());
    CHECK(std::abs(volume_of(brep) - item.volume) <= 1e-6);
    CHECK(brep.face_ids().size() == item.faces);
  }
}

TEST_CASE("brep: boolean_solid_curved defers a box minus a curved cap operand") {
  const double tolerance = 1e-9;
  const geom::Brep quarter = make_exact_quarter_cylinder(1.0, 2.0, tolerance);
  const geom::Brep box =
      geom::make_box(geom::Point3d{-1.0, -1.0, 1.0}, geom::Point3d{10.0, 10.0, 10.0}, tolerance)
          .value();
  REQUIRE(geom::validate_topology(quarter).ok());
  REQUIRE(geom::validate_topology(box).ok());

  // The kept cavity fragments of the quarter would have to be emitted with the
  // lateral face flipped. That is deferred, and the driver has to say so before
  // it starts cutting instead of failing somewhere inside the split.
  const geom::Result<geom::CurvedBooleanResult> result = geom::boolean_solid_curved(
      box, only_solid(box), quarter, only_solid(quarter), geom::SolidBooleanOp::difference,
      geom::CurvedBooleanOptions{tolerance});
  REQUIRE(!result.has_value());
  CHECK(result.error().code == geom::Errc::unsupported);
  CHECK(message_of(result.error()) ==
        "boolean_solid_curved: flipped curved fragments are not supported yet");
}

TEST_CASE("brep: split_face_curved rejects a cut that ends inside the face") {
  const double tolerance = 1e-9;
  const geom::Brep cylinder = make_test_cylinder(tolerance);
  geom::FaceId lateral{};
  bool found = false;
  for (const geom::FaceId id : cylinder.face_ids()) {
    const geom::SurfaceGeometry& surface =
        cylinder.geometry().surfaces[cylinder.face(id).surface];
    if (std::holds_alternative<geom::Cylinder3>(surface)) {
      lateral = id;
      found = true;
    }
  }
  REQUIRE(found);

  // The lower end stops 1e-2 short of the v = 0 rim, far beyond the junction
  // tolerance (about 6.6e-5 in uv); the reconciliation pass must not absorb
  // that gap.
  std::vector<geom::FaceCut> cuts(1);
  cuts[0].uv = {geom::Point2d{1.0, 0.01}, geom::Point2d{1.0, 1.0}};
  const geom::Result<std::vector<geom::CurvedFragment>> result =
      geom::split_face_curved(cylinder, lateral, cuts, tolerance);
  REQUIRE(!result.has_value());
  CHECK(result.error().code == geom::Errc::numeric_failure);
  CHECK(message_of(result.error()) ==
        "split_face_curved: an open cut endpoint is not on the face boundary");
}

TEST_CASE("brep: collapse_circular_chains heals a rim run beside the cut corner") {
  const double tolerance = 1e-9;
  const geom::CircleCurve3 circle{geom::Point3d{0.0, 0.0, 0.5}, geom::Vec3d{1.0, 0.0, 0.0},
                                  geom::Vec3d{0.0, 1.0, 0.0}, 1.0};
  const auto at_degrees = [&circle](double degrees) {
    const double radians = degrees * kPi / 180.0;
    return circle.center + circle.x_axis * (circle.radius * std::cos(radians)) +
           circle.y_axis * (circle.radius * std::sin(radians));
  };

  // The failing rim of the offset-cap fixture: a 30 degree arc sampled every
  // 1.875 degrees, ending exactly at the cut corner (cos(-30), sin(-30), 0.5)
  // where the box plane y = -0.5 crosses the circle, then the straight clipped
  // chords back to the start of the arc.
  geom::CurvedFragment fragment;
  std::vector<geom::Point3d>& loop = fragment.loops.emplace_back();
  for (int step = 0; step <= 16; ++step) {
    loop.push_back(at_degrees(-60.0 + 1.875 * static_cast<double>(step)));
  }
  loop.push_back(geom::Point3d{0.5, -0.5, 0.5});
  fragment.loops_uv.emplace_back(loop.size(), geom::Point2d{0.0, 0.0});
  CHECK(std::abs(loop[16][0] - 0.8660254037844386) <= 1e-12);
  // Copy the corner: the collapse replaces the whole ring vector.
  const geom::Point3d corner = loop[16];

  const geom::Result<bool> collapsed =
      geom::detail::collapse_circular_chains(fragment, circle, tolerance);
  REQUIRE_MESSAGE(collapsed.has_value(), "collapse failed: " << message_of(collapsed.error()));
  CHECK(collapsed.value());

  // The dense samples collapse into one span whose sweep advances in a single
  // direction (no backtracking), and the exact corner stays a ring vertex.
  REQUIRE(fragment.loops[0].size() == 3);
  REQUIRE(fragment.arc_spans[0].size() == fragment.loops[0].size());
  std::size_t engaged = fragment.arc_spans[0].size();
  for (std::size_t i = 0; i < fragment.arc_spans[0].size(); ++i) {
    if (!fragment.arc_spans[0][i].has_value()) continue;
    CHECK(engaged == fragment.arc_spans[0].size());
    engaged = i;
  }
  REQUIRE(engaged < fragment.arc_spans[0].size());
  const double sweep =
      fragment.arc_spans[0][engaged]->t1 - fragment.arc_spans[0][engaged]->t0;
  CHECK(std::abs(sweep - kPi / 6.0) <= 1e-9);
  bool keeps_corner = false;
  for (const geom::Point3d& point : fragment.loops[0]) {
    if (std::abs(point[0] - corner[0]) <= 1e-15 && std::abs(point[1] - corner[1]) <= 1e-15) {
      keeps_corner = true;
    }
  }
  CHECK(keeps_corner);
}
// The QA corpus: six positive cases with independently derived volumes and five
// diagnosed rejections. The positive volumes are
//   I = pi/6 + sqrt(3)/4 - 1/2                  (offset-cap footprint, unit tall)
//   offset box [0.5, 2] x [-0.5, 0.5] x [0.5, 1.5] on the unit cylinder:
//     intersection = I, union = 1.5 + 2*pi - I, difference = 2*pi - I
//   seam box [-2, 2]^2 x [1, 10]: union = 144 + pi
//   slab [0.5, 2] x [-2, 2] x [-1, 10]: intersection = 2*(pi/3 - sqrt(3)/4)
//   box [-0.2, 0.4] x [-2, 2] x [1.5, 3]: intersection = 0.587720000938..., i.e.
//   half of the disk strip x in [-0.2, 0.4].
TEST_CASE("brep: boolean_solid_curved QA corpus") {
  const double tolerance = 1e-9;
  const double i = kPi / 6.0 + std::sqrt(3.0) / 4.0 - 0.5;

  const geom::Brep cylinder = make_test_cylinder(tolerance);
  const geom::Brep offset_box =
      geom::make_box(geom::Point3d{0.5, -0.5, 0.5}, geom::Point3d{2.0, 0.5, 1.5}, tolerance)
          .value();
  const geom::Brep seam_box =
      geom::make_box(geom::Point3d{-2.0, -2.0, 1.0}, geom::Point3d{2.0, 2.0, 10.0}, tolerance)
          .value();
  const geom::Brep slab_box =
      geom::make_box(geom::Point3d{0.5, -2.0, -1.0}, geom::Point3d{2.0, 2.0, 10.0}, tolerance)
          .value();
  const geom::Brep caps_box =
      geom::make_box(geom::Point3d{-0.2, -2.0, 1.5}, geom::Point3d{0.4, 2.0, 3.0}, tolerance)
          .value();
  // Sphere ∩ box z >= 0.5: the cap of height 0.5 above the latitude z = 0.5.
  const geom::Brep sphere =
      geom::make_sphere(geom::Point3d{0.0, 0.0, 0.0}, 1.0, tolerance).value();
  const geom::Brep sphere_cap_box =
      geom::make_box(geom::Point3d{-2.0, -2.0, 0.5}, geom::Point3d{2.0, 2.0, 2.0}, tolerance)
          .value();

  struct PositiveCase {
    const char* name;
    const geom::Brep* a;
    const geom::Brep* b;
    geom::SolidBooleanOp op;
    double expected;
  };
  const std::array<PositiveCase, 7> positives{{
      {"offset cap union", &cylinder, &offset_box, geom::SolidBooleanOp::union_,
       1.5 + 2.0 * kPi - i},
      {"offset cap intersection", &cylinder, &offset_box, geom::SolidBooleanOp::intersection, i},
      {"offset cap difference", &cylinder, &offset_box, geom::SolidBooleanOp::difference,
       2.0 * kPi - i},
      {"seam box union", &cylinder, &seam_box, geom::SolidBooleanOp::union_, 144.0 + kPi},
      {"slab intersection", &cylinder, &slab_box, geom::SolidBooleanOp::intersection,
       2.0 * (kPi / 3.0 - std::sqrt(3.0) / 4.0)},
      {"box through both caps intersection", &cylinder, &caps_box,
       geom::SolidBooleanOp::intersection, 0.587720000938},
      {"sphere cap intersection", &sphere, &sphere_cap_box,
       geom::SolidBooleanOp::intersection, 0.6544984694978736},
  }};
  for (const PositiveCase& entry : positives) {
    check_qa_volume(entry.name, *entry.a, *entry.b, entry.op, entry.expected);
  }

  // Diagnosed rejections. The expected messages are the ones the driver
  // documents:
  //   cone     "boolean_solid_curved: cone operands need closed-cut and apex
  //             topology; booleans are deferred to M3c-3"
  //   NURBS    "boolean_solid_curved: NURBS operands are not supported"
  //   tangent  "boolean_solid_curved: tangent face contact is not supported yet"
  //   coplanar "boolean_solid_curved: coplanar face contact is not supported yet"
  const geom::Brep cone =
      geom::make_cone(geom::Point3d{0.0, 0.0, 0.0}, geom::Vec3d{0.0, 0.0, 1.0}, 1.0, 1.0,
                      tolerance)
          .value();
  const geom::Brep nurbs = make_nurbs_box(geom::Point3d{0.0, 0.0, 0.0},
                                          geom::Point3d{1.0, 1.0, 1.0}, tolerance);
  const geom::Brep wide_box =
      geom::make_box(geom::Point3d{-5.0, -5.0, -5.0}, geom::Point3d{5.0, 5.0, 5.0}, tolerance)
          .value();
  // The box faces x = +-1 and y = +-1 are tangent to the unit cylinder.
  const geom::Brep tangent_box =
      geom::make_box(geom::Point3d{-1.0, -1.0, -1.0}, geom::Point3d{1.0, 1.0, 1.0}, tolerance)
          .value();
  // The box bottom face lies in the top cap plane (z = 2) and overlaps it.
  const geom::Brep resting_box =
      geom::make_box(geom::Point3d{-0.5, -0.5, 2.0}, geom::Point3d{0.5, 0.5, 3.0}, tolerance)
          .value();

  struct RejectionCase {
    const char* name;
    const geom::Brep* a;
    const geom::Brep* b;
    geom::SolidBooleanOp op;
    const char* message;
  };
  const std::array<RejectionCase, 4> rejections{{
      {"cone box", &cone, &wide_box, geom::SolidBooleanOp::union_,
       "boolean_solid_curved: cone operands need closed-cut and apex topology; booleans are "
       "deferred to M3c-3"},
      {"NURBS box", &nurbs, &wide_box, geom::SolidBooleanOp::union_,
       "boolean_solid_curved: NURBS operands are not supported"},
      {"tangent box", &tangent_box, &cylinder, geom::SolidBooleanOp::union_,
       "boolean_solid_curved: tangent face contact is not supported yet"},
      {"coplanar cap", &cylinder, &resting_box, geom::SolidBooleanOp::union_,
       "boolean_solid_curved: coplanar face contact is not supported yet"},
  }};
  for (const RejectionCase& entry : rejections) {
    check_qa_rejection(entry.name, *entry.a, *entry.b, entry.op, entry.message);
  }
}


TEST_CASE("curved boolean: sphere/sphere is not a splitter problem") {
  const double tolerance = 1e-9;
  const geom::CurvedBooleanOptions options{tolerance};
  const geom::Brep unit_a =
      geom::make_sphere(geom::Point3d{0.0, 0.0, 0.0}, 1.0, tolerance).value();
  const geom::Brep unit_b =
      geom::make_sphere(geom::Point3d{1.5, 0.0, 0.0}, 1.0, tolerance).value();
  const BrepFingerprint unit_a_before = fingerprint_of(unit_a);
  const BrepFingerprint unit_b_before = fingerprint_of(unit_b);

  // Two equal spherical caps of height h = r - d / 2 = 0.25 meet at the rim.
  const geom::Result<geom::CurvedBooleanResult> lens = geom::boolean_solid_curved(
      unit_a, only_solid(unit_a), unit_b, only_solid(unit_b),
      geom::SolidBooleanOp::intersection, options);
  REQUIRE_MESSAGE(lens.has_value(), message_of(lens.error()));
  CHECK(geom::validate_topology(lens.value().brep).ok());
  CHECK(lens.value().brep.face_ids().size() == 2);
  CHECK(lens.value().brep.edge_ids().size() == 3);
  CHECK(lens.value().brep.vertex_ids().size() == 3);

  // Both caps are charted on the shared frame: each face is the polar cap above
  // the rim latitude v = pi / 2 - acos(x_rim) and it reaches the pole at pi / 2.
  const double rim_latitude = 0.5 * kPi - std::acos(0.75);
  for (const geom::FaceId face_id : lens.value().brep.face_ids()) {
    const geom::Result<std::vector<geom::FaceUvLoop>> loops =
        geom::face_loops_uv(lens.value().brep, face_id, tolerance);
    REQUIRE_MESSAGE(loops.has_value(), message_of(loops.error()));
    REQUIRE(loops.value().size() == 1);
    const std::vector<geom::Point2d>& ring = loops.value().front().points;
    REQUIRE(ring.size() >= 3);
    double min_v = 1e30;
    double max_v = -1e30;
    std::size_t poles = 0;
    for (const geom::Point2d& point : ring) {
      min_v = std::min(min_v, point[1]);
      max_v = std::max(max_v, point[1]);
      if (std::abs(point[1] - 0.5 * kPi) <= 1e-9) ++poles;
    }
    CHECK(poles == 1);
    CHECK(std::abs(max_v - 0.5 * kPi) <= 1e-9);
    CHECK(std::abs(min_v - rim_latitude) <= 1e-9);
  }

  // The output carries the re-charted frame, not the operands' canonical +z
  // charts: A's cap is on +x, B's on -x, and both share one x-axis meridian.
  std::vector<const geom::Sphere3*> cap_surfaces;
  for (const geom::FaceId face_id : lens.value().brep.face_ids()) {
    const geom::Face& face = lens.value().brep.face(face_id);
    const geom::SurfaceGeometry& surface = lens.value().brep.geometry().surfaces[face.surface];
    const geom::Sphere3* sphere = std::get_if<geom::Sphere3>(&surface);
    REQUIRE_MESSAGE(sphere != nullptr, "a lens face is not on a sphere");
    cap_surfaces.push_back(sphere);
  }
  REQUIRE(cap_surfaces.size() == 2);
  const geom::Vec3d plus_x{1.0, 0.0, 0.0};
  const geom::Vec3d minus_x{-1.0, 0.0, 0.0};
  const bool first_is_plus =
      std::abs(geom::dot(cap_surfaces[0]->axis, plus_x) - 1.0) <= 1e-12;
  const bool first_is_minus =
      std::abs(geom::dot(cap_surfaces[0]->axis, minus_x) - 1.0) <= 1e-12;
  CHECK(first_is_plus != first_is_minus);
  CHECK(std::abs(geom::dot(cap_surfaces[0]->x_axis, cap_surfaces[1]->x_axis) - 1.0) <= 1e-12);
  CHECK(std::abs(geom::dot(cap_surfaces[1]->axis, first_is_plus ? minus_x : plus_x) - 1.0) <=
        1e-12);

  // The two caps are stitched on one rim edge used in opposite senses; each
  // seam edge is a doubled use on a single face.
  std::optional<geom::EdgeId> rim_edge;
  for (const geom::EdgeId edge_id : lens.value().brep.edge_ids()) {
    std::vector<const geom::Coedge*> uses;
    for (const geom::CoedgeId coedge_id : lens.value().brep.coedge_ids()) {
      const geom::Coedge& coedge = lens.value().brep.coedge(coedge_id);
      if (coedge.edge == edge_id) uses.push_back(&coedge);
    }
    if (uses.size() == 2 && uses[0]->face != uses[1]->face) {
      REQUIRE(!rim_edge.has_value());
      rim_edge = edge_id;
      CHECK(uses[0]->reversed != uses[1]->reversed);
    }
  }
  REQUIRE(rim_edge.has_value());
  CHECK(fingerprint_of(unit_b) == unit_b_before);

  // The sphere check needs whole-sphere operands: a sphere intersected with a
  // box gives a two-face solid whose spherical face is only a cap, so the pair
  // is diagnosed instead of reaching the splitter.
  const geom::Brep cap_box =
      geom::make_box(geom::Point3d{-2.0, -2.0, 0.5}, geom::Point3d{2.0, 2.0, 2.0}, tolerance)
          .value();
  const geom::Result<geom::CurvedBooleanResult> cap_result = geom::boolean_solid_curved(
      unit_a, only_solid(unit_a), cap_box, only_solid(cap_box),
      geom::SolidBooleanOp::intersection, options);
  REQUIRE_MESSAGE(cap_result.has_value(), message_of(cap_result.error()));
  const geom::Brep cap = cap_result.value().brep;
  REQUIRE(cap.face_ids().size() == 2);
  const geom::Brep cap_other =
      geom::make_sphere(geom::Point3d{0.0, 0.0, 1.5}, 1.0, tolerance).value();

  const BrepFingerprint cap_before = fingerprint_of(cap);
  const BrepFingerprint cap_other_before = fingerprint_of(cap_other);
  const geom::Result<geom::CurvedBooleanResult> not_whole = geom::boolean_solid_curved(
      cap, only_solid(cap), cap_other, only_solid(cap_other),
      geom::SolidBooleanOp::intersection, options);
  REQUIRE_MESSAGE(!not_whole.has_value(), "a cap solid and a whole sphere must be rejected");
  CHECK(not_whole.error().code == geom::Errc::unsupported);
  CHECK(message_of(not_whole.error()) ==
        "boolean_solid_curved: sphere/sphere booleans need whole-sphere operands");
  CHECK(fingerprint_of(cap) == cap_before);
  CHECK(fingerprint_of(cap_other) == cap_other_before);

  // External tangency (unit/unit at d = 2), containment (d + 0.2 <= 1) and
  // separation (unit/unit at d = 2.5) have no circle of intersection.
  struct CrossingCase {
    const char* name;
    geom::Point3d center;
    double radius;
  };
  const std::array<CrossingCase, 3> crossing{{
      {"external tangency", geom::Point3d{2.0, 0.0, 0.0}, 1.0},
      {"containment", geom::Point3d{0.5, 0.0, 0.0}, 0.2},
      {"separation", geom::Point3d{2.5, 0.0, 0.0}, 1.0},
  }};
  for (const CrossingCase& entry : crossing) {
    INFO(entry.name);
    const geom::Brep other = geom::make_sphere(entry.center, entry.radius, tolerance).value();
    const BrepFingerprint other_before = fingerprint_of(other);
    const geom::Result<geom::CurvedBooleanResult> rejected = geom::boolean_solid_curved(
        unit_a, only_solid(unit_a), other, only_solid(other),
        geom::SolidBooleanOp::intersection, options);
    REQUIRE_MESSAGE(!rejected.has_value(), "expected the crossing-spheres diagnostic");
    CHECK(rejected.error().code == geom::Errc::unsupported);
    CHECK(message_of(rejected.error()) ==
          "boolean_solid_curved: sphere/sphere booleans need two crossing spheres");
    CHECK(fingerprint_of(other) == other_before);
  }

  // The re-chart uses the unit centre line directly, so a +y centre line has to
  // take the same path with unequal radii instead of only the canonical +x one.
  const geom::Brep y_a =
      geom::make_sphere(geom::Point3d{0.0, 0.0, 0.0}, 1.2, tolerance).value();
  const geom::Brep y_b =
      geom::make_sphere(geom::Point3d{0.0, 1.4, 0.0}, 0.8, tolerance).value();
  const geom::Result<geom::CurvedBooleanResult> y_lens = geom::boolean_solid_curved(
      y_a, only_solid(y_a), y_b, only_solid(y_b), geom::SolidBooleanOp::intersection, options);
  REQUIRE_MESSAGE(y_lens.has_value(), message_of(y_lens.error()));
  CHECK(geom::validate_topology(y_lens.value().brep).ok());
  CHECK(y_lens.value().brep.face_ids().size() == 2);
  CHECK(y_lens.value().brep.edge_ids().size() == 3);
  CHECK(y_lens.value().brep.vertex_ids().size() == 3);

  CHECK(fingerprint_of(unit_a) == unit_a_before);
}

namespace {
// Defined with the M3c-3b lens corpus later in this file; the off-axis case
// below reuses the same rim check so a rotated lens cannot close by a
// different rule.
void check_one_shared_rim(const char* name, const geom::Brep& brep);
}  // namespace

TEST_CASE("curved boolean: an off-axis sphere lens is built") {
  const double tolerance = 1e-9;
  const double angle = 15.0 * 3.14159265358979323846 / 180.0;
  const geom::Result<geom::Brep> south =
      geom::make_sphere(geom::Point3d{0.0, 0.0, 0.0}, 1.0, tolerance);
  const geom::Result<geom::Brep> north = geom::make_sphere(
      geom::Point3d{1.5 * std::cos(angle), 1.5 * std::sin(angle), 0.0}, 1.0,
      tolerance);
  REQUIRE(south.has_value());
  REQUIRE(north.has_value());

  const geom::Result<geom::CurvedBooleanResult> lens = geom::boolean_solid_curved(
      south.value(), south.value().solid_ids().front(), north.value(),
      north.value().solid_ids().front(), geom::SolidBooleanOp::intersection,
      geom::CurvedBooleanOptions{tolerance});
  REQUIRE_MESSAGE(lens.has_value(), message_of(lens.error()));
  CHECK(geom::validate_topology(lens.value().brep).ok());
  check_one_shared_rim("off-axis lens", lens.value().brep);
  const std::optional<double> volume = sphere_face_volume(lens.value().brep);
  REQUIRE(volume.has_value());
  CHECK(*volume == doctest::Approx(11.0 * 3.14159265358979323846 / 96.0)
                      .epsilon(1e-6));

  // The union is the complementary two-cap solid, so the same rim must also
  // stitch the outside of the lens: 8/3*pi - 11*pi/96 = 245*pi/96.
  const geom::Result<geom::CurvedBooleanResult> united = geom::boolean_solid_curved(
      south.value(), south.value().solid_ids().front(), north.value(),
      north.value().solid_ids().front(), geom::SolidBooleanOp::union_,
      geom::CurvedBooleanOptions{tolerance});
  REQUIRE_MESSAGE(united.has_value(), message_of(united.error()));
  CHECK(geom::validate_topology(united.value().brep).ok());
  check_one_shared_rim("off-axis union", united.value().brep);
  CHECK(*sphere_face_volume(united.value().brep) ==
        doctest::Approx(245.0 * 3.14159265358979323846 / 96.0).epsilon(1e-6));
}

TEST_CASE("curved boolean: sphere/sphere re-charts opposite meridians") {
  const double tolerance = 1e-9;
  const geom::CurvedBooleanOptions options{tolerance};
  // Both operands are already charted on the centre line (+x / -x), but their
  // seams lie on opposite meridians (+y versus -y). The recursion fixed point
  // has to include x_axis agreement, or the lens is stitched at two meridians.
  const geom::Brep a = geom::make_sphere(geom::Point3d{0.0, 0.0, 0.0},
                                         geom::Vec3d{1.0, 0.0, 0.0},
                                         geom::Vec3d{0.0, 1.0, 0.0}, 1.0, tolerance)
                           .value();
  const geom::Brep b = geom::make_sphere(geom::Point3d{1.5, 0.0, 0.0},
                                         geom::Vec3d{-1.0, 0.0, 0.0},
                                         geom::Vec3d{0.0, -1.0, 0.0}, 1.0, tolerance)
                           .value();
  const geom::Result<geom::CurvedBooleanResult> lens = geom::boolean_solid_curved(
      a, only_solid(a), b, only_solid(b), geom::SolidBooleanOp::intersection, options);
  REQUIRE_MESSAGE(lens.has_value(), message_of(lens.error()));
  CHECK(geom::validate_topology(lens.value().brep).ok());
  CHECK(lens.value().brep.face_ids().size() == 2);
  CHECK(lens.value().brep.edge_ids().size() == 3);
  CHECK(lens.value().brep.vertex_ids().size() == 3);
}

TEST_CASE("curved boolean: difference subtracts a sphere cap") {
  const double tolerance = 1e-9;
  const geom::CurvedBooleanOptions options{tolerance};
  const geom::Brep unit_a =
      geom::make_sphere(geom::Point3d{0.0, 0.0, 0.0}, 1.0, tolerance).value();
  const geom::Brep unit_b =
      geom::make_sphere(geom::Point3d{1.5, 0.0, 0.0}, 1.0, tolerance).value();

  // Equal unit spheres at d = 1.5 meet in caps of height 0.25: the lens is
  // 2 * pi * 0.25^2 * (3 - 0.25) / 3 = 11*pi/96, so either difference is
  // 4*pi/3 - 11*pi/96 = 39*pi/32.
  const double expected = 39.0 * std::numbers::pi / 32.0;
  struct DifferenceCase {
    const char* name;
    const geom::Brep* a;
    const geom::Brep* b;
    // x component of the chart axis of the subtracted (reversed) cap.
    double reversed_axis_x;
  };
  const std::array<DifferenceCase, 2> cases{{
      {"A - B", &unit_a, &unit_b, -1.0},
      {"B - A", &unit_b, &unit_a, 1.0},
  }};
  for (const DifferenceCase& entry : cases) {
    INFO(entry.name);
    const geom::Result<geom::CurvedBooleanResult> result = geom::boolean_solid_curved(
        *entry.a, only_solid(*entry.a), *entry.b, only_solid(*entry.b),
        geom::SolidBooleanOp::difference, options);
    REQUIRE_MESSAGE(result.has_value(), entry.name << ": " << message_of(result.error()));
    const geom::Brep& brep = result.value().brep;
    CHECK(geom::validate_topology(brep).ok());
    CHECK(brep.face_ids().size() == 2);
    CHECK(brep.edge_ids().size() == 3);
    CHECK(brep.vertex_ids().size() == 3);
    CHECK(brep.solid_shells(only_solid(brep)).size() == 1);

    const std::optional<double> volume = sphere_face_volume(brep);
    REQUIRE_MESSAGE(volume.has_value(), entry.name << ": a face is not a polar cap");
    CHECK(*volume == doctest::Approx(expected).epsilon(1e-9));

    // Only the subtracted B cap is flipped; the A complement stays outward.
    // `curved_solid_volume` takes each face's sign from the `reversed` flag,
    // not from the ring winding, so exactly one face may carry it, and it must
    // be the smaller cap (cap height 0.25: 11*pi/192 = 0.179987...).
    std::vector<geom::FaceId> reversed_faces;
    std::optional<double> reversed_cap_volume;
    for (const geom::FaceId face_id : brep.face_ids()) {
      const geom::Face& face = brep.face(face_id);
      const geom::SurfaceGeometry& surface = brep.geometry().surfaces[face.surface];
      if (!std::holds_alternative<geom::Sphere3>(surface)) continue;
      const geom::Sphere3& sphere = std::get<geom::Sphere3>(surface);
      // Re-charting puts the complement on the first operand's axis and the
      // subtracted cap on the opposite one; B - A swaps which chart is which.
      if (face.reversed) {
        CHECK(sphere.axis[0] == doctest::Approx(entry.reversed_axis_x).epsilon(1e-12));
      } else {
        CHECK(sphere.axis[0] == doctest::Approx(-entry.reversed_axis_x).epsilon(1e-12));
      }
      if (!face.reversed) continue;
      reversed_faces.push_back(face_id);
      reversed_cap_volume = sphere_face_cap_volume(brep, face_id);
    }
    REQUIRE(reversed_faces.size() == 1);
    REQUIRE_MESSAGE(reversed_cap_volume.has_value(), entry.name << ": reversed face is not a cap");
    const double smaller_cap = 11.0 * std::numbers::pi / 192.0;
    CHECK(std::abs(*reversed_cap_volume) == doctest::Approx(smaller_cap).epsilon(1e-9));

    // The A complement cap and the reversed B cap meet on one rim edge used by
    // both faces; each of the two seam edges is used twice by one face. The two
    // faces walk the rim in opposite directions, so their coedges disagree on
    // whether they run against the edge's own orientation.
    std::size_t rim_edges = 0;
    for (const geom::EdgeId edge_id : brep.edge_ids()) {
      std::vector<geom::FaceId> users;
      std::vector<bool> senses;
      std::size_t uses = 0;
      for (const geom::CoedgeId coedge_id : brep.coedge_ids()) {
        const geom::Coedge& coedge = brep.coedge(coedge_id);
        if (coedge.edge != edge_id) continue;
        ++uses;
        senses.push_back(coedge.reversed);
        if (std::find(users.begin(), users.end(), coedge.face) == users.end()) {
          users.push_back(coedge.face);
        }
      }
      CHECK(uses == 2);
      if (users.size() == 2) {
        ++rim_edges;
        REQUIRE(senses.size() == 2);
        CHECK(senses[0] != senses[1]);
      }
    }
    CHECK(rim_edges == 1);
  }

  // A spherical band has a pole at neither v end, so the helper must refuse it
  // instead of reporting zero volume for it.
  const geom::Brep slab =
      geom::make_box(geom::Point3d{-2.0, -2.0, -0.8}, geom::Point3d{2.0, 2.0, 0.8}, tolerance)
          .value();
  const geom::Result<geom::CurvedBooleanResult> band = geom::boolean_solid_curved(
      unit_a, only_solid(unit_a), slab, only_solid(slab), geom::SolidBooleanOp::intersection,
      options);
  REQUIRE_MESSAGE(band.has_value(), message_of(band.error()));
  CHECK(!sphere_face_volume(band.value().brep).has_value());
}

// QA-16: the M3c-3b sphere/sphere lens corpus. Three operations on unit-sphere
// lenses (intersection 11*pi/96, union 245*pi/96, difference 39*pi/32) and one
// unequal-radii lens on a non-zero origin. Each row must be a validated
// two-face single-shell solid that closes by coedge uses and carries exactly one
// shared rim edge; a repeated call must reproduce the chart volume bits and the
// topology counts.
//
// The lens rim at tolerance 1e-9 carries about 4098 points and
// `detail::curved_solid_volume` stalls on it in ear clipping (Task 4), so the
// rows measure the test-local `sphere_face_volume` chart sum instead. That helper
// is cross-checked once against `curved_solid_volume` on the M3c-3a cap fixture,
// whose single analytic rim circle the integrator does handle.
namespace {

struct LensTopology {
  std::size_t solids = 0;
  std::size_t shells = 0;
  std::size_t faces = 0;
  std::size_t edges = 0;
  std::size_t coedges = 0;
  std::size_t vertices = 0;
  bool operator==(const LensTopology&) const = default;
};

LensTopology lens_topology(const geom::Brep& brep) {
  LensTopology topology;
  topology.solids = brep.solid_ids().size();
  topology.faces = brep.face_ids().size();
  topology.edges = brep.edge_ids().size();
  topology.coedges = brep.coedge_ids().size();
  topology.vertices = brep.vertex_ids().size();
  if (topology.solids == 1) {
    topology.shells = brep.solid_shells(brep.solid_ids().front()).size();
  }
  return topology;
}

// Exactly one edge may be used by two different faces, and those two coedges
// must walk it in opposite directions; every other edge is used twice by one
// face. For these lenses the shared edge is the rim circle.
void check_one_shared_rim(const char* name, const geom::Brep& brep) {
  std::size_t rim_edges = 0;
  std::size_t edge_index = 0;
  for (const geom::EdgeId edge_id : brep.edge_ids()) {
    const std::size_t current = edge_index++;
    std::vector<geom::FaceId> users;
    std::vector<bool> senses;
    for (const geom::CoedgeId coedge_id : brep.coedge_ids()) {
      const geom::Coedge& coedge = brep.coedge(coedge_id);
      if (coedge.edge != edge_id) continue;
      senses.push_back(coedge.reversed);
      if (std::find(users.begin(), users.end(), coedge.face) == users.end()) {
        users.push_back(coedge.face);
      }
    }
    CHECK_MESSAGE(senses.size() == 2, name << ": edge " << current << " is not used twice");
    if (users.size() != 2) continue;
    ++rim_edges;
    REQUIRE(senses.size() == 2);
    CHECK_MESSAGE(senses[0] != senses[1], name << ": the rim coedges agree on orientation");
  }
  CHECK_MESSAGE(rim_edges == 1, name << ": expected one shared rim edge, saw " << rim_edges);
}

// One QA-16 lens run. `volume` is the signed chart cap sum from
// `sphere_face_volume`; its bits are kept so a repeat can be compared exactly.
struct LensRun {
  LensTopology topology;
  double volume = 0.0;
  std::uint64_t volume_bits = 0;
  std::size_t cuts = 0;
  std::size_t kept = 0;
};

LensRun run_lens_row(const char* name, const geom::Brep& a, const geom::Brep& b,
                     geom::SolidBooleanOp op, double expected, double relative_tolerance) {
  INFO(name);
  const geom::Result<geom::CurvedBooleanResult> result = geom::boolean_solid_curved(
      a, only_solid(a), b, only_solid(b), op, geom::CurvedBooleanOptions{1e-9});
  REQUIRE_MESSAGE(result.has_value(), name << ": " << message_of(result.error()));
  const geom::Brep& brep = result.value().brep;
  LensRun run;
  run.topology = lens_topology(brep);
  CHECK_MESSAGE(run.topology.solids == 1, name << ": expected one solid, saw " << run.topology.solids);
  CHECK_MESSAGE(run.topology.shells == 1, name << ": expected one shell, saw " << run.topology.shells);
  CHECK_MESSAGE(run.topology.faces == 2, name << ": expected two faces, saw " << run.topology.faces);
  CHECK_MESSAGE(geom::validate_topology(brep).ok(), name << ": invalid topology");
  CHECK_MESSAGE(geom::detail::heal_detail::shell_is_closed_by_uses(
                    brep, brep.solid_shells(only_solid(brep)).front()),
                name << ": an edge is not used by exactly two coedges");
  check_one_shared_rim(name, brep);

  const std::optional<double> chart_volume = sphere_face_volume(brep);
  REQUIRE_MESSAGE(chart_volume.has_value(), name << ": a face is not a polar cap");
  run.volume = *chart_volume;
  run.volume_bits = std::bit_cast<std::uint64_t>(run.volume);
  const double residual = std::abs(run.volume - expected);
  CHECK_MESSAGE(residual <= relative_tolerance * std::abs(expected),
                name << ": chart volume " << run.volume << " expected " << expected << " residual "
                     << residual);
  run.cuts = result.value().report.cuts;
  run.kept = result.value().report.kept;
  return run;
}

// Runs the row twice and requires the two runs to agree exactly.
void check_lens_row(const char* name, const geom::Brep& a, const geom::Brep& b,
                    geom::SolidBooleanOp op, double expected, double relative_tolerance) {
  INFO(name);
  const LensRun first = run_lens_row(name, a, b, op, expected, relative_tolerance);
  const LensRun second = run_lens_row(name, a, b, op, expected, relative_tolerance);
  CHECK_MESSAGE(first.volume_bits == second.volume_bits, name << ": chart volume bits differ");
  CHECK_MESSAGE(first.topology == second.topology, name << ": topology counts differ");
  CHECK_MESSAGE(first.cuts == second.cuts, name << ": cut counts differ");
  CHECK_MESSAGE(first.kept == second.kept, name << ": kept fragment counts differ");
}

}  // namespace

TEST_CASE("QA-16: sphere lens booleans are closed and deterministic") {
  constexpr double tolerance = 1e-9;

  // Helper cross-check on the M3c-3a cap fixture: the library volume must
  // reproduce the analytic cap (h = 0.5) and the chart helper must agree with it
  // to 1e-12. The fixture is the two-face cap solid whose two faces share the one
  // analytic rim circle.
  const geom::Brep unit_sphere =
      geom::make_sphere(geom::Point3d{0.0, 0.0, 0.0}, 1.0, tolerance).value();
  const geom::Brep cap_box =
      geom::make_box(geom::Point3d{-2.0, -2.0, 0.5}, geom::Point3d{2.0, 2.0, 2.0}, tolerance)
          .value();
  {
    const geom::Result<geom::CurvedBooleanResult> cap = geom::boolean_solid_curved(
        unit_sphere, only_solid(unit_sphere), cap_box, only_solid(cap_box),
        geom::SolidBooleanOp::intersection, geom::CurvedBooleanOptions{tolerance});
    REQUIRE_MESSAGE(cap.has_value(), "cap fixture: " << message_of(cap.error()));
    const geom::Brep& cap_brep = cap.value().brep;
    const geom::Result<double> library =
        geom::detail::curved_solid_volume(cap_brep, only_solid(cap_brep), tolerance);
    REQUIRE_MESSAGE(library.has_value(), "cap fixture: " << message_of(library.error()));
    constexpr double kCapReference = 0.6544984694978736;
    CHECK_MESSAGE(std::abs(library.value() - kCapReference) <= 1e-6 * kCapReference,
                  "cap fixture: library volume " << library.value() << " expected "
                                                << kCapReference);
    const std::optional<double> helper = sphere_face_volume(cap_brep);
    REQUIRE_MESSAGE(helper.has_value(), "cap fixture: a face is not a polar cap");
    CHECK_MESSAGE(std::abs(*helper - library.value()) <= 1e-12 * std::abs(library.value()),
                  "cap fixture: helper " << *helper << " vs library " << library.value());

    // Pin the fixture: the spherical hat and its planar cap share one rim edge.
    CHECK_MESSAGE(cap_brep.face_ids().size() == 2,
                  "cap fixture: expected two faces, saw " << cap_brep.face_ids().size());
    check_one_shared_rim("cap fixture", cap_brep);
  }

  // Equal unit spheres 1.5 apart: the pair circle has radius sqrt(7)/4 and caps
  // of height 0.25. The lens is 2 * pi * 0.25^2 * (3 - 0.25) / 3 = 11*pi/96, the
  // union is two complements (4*pi/3 - 11*pi/192 each), and either difference is
  // 4*pi/3 - 11*pi/96 = 39*pi/32. Only the intersection match is exact in the
  // chart formula, so it is held to 1e-9; the union and difference are held to
  // 1e-6 like the rest of the curved boolean corpus.
  const geom::Brep lens_a =
      geom::make_sphere(geom::Point3d{0.0, 0.0, 0.0}, 1.0, tolerance).value();
  const geom::Brep lens_b =
      geom::make_sphere(geom::Point3d{1.5, 0.0, 0.0}, 1.0, tolerance).value();
  struct Row {
    const char* name;
    geom::SolidBooleanOp op;
    double expected;
    double relative_tolerance;
  };
  const std::array<Row, 3> rows{{
      {"intersection", geom::SolidBooleanOp::intersection, 11.0 * kPi / 96.0, 1e-9},
      {"union", geom::SolidBooleanOp::union_, 245.0 * kPi / 96.0, 1e-6},
      {"difference", geom::SolidBooleanOp::difference, 39.0 * kPi / 32.0, 1e-6},
  }};
  for (const Row& row : rows) {
    check_lens_row(row.name, lens_a, lens_b, row.op, row.expected, row.relative_tolerance);
  }

  // Generality: the unequal pair from the M3c-3b alignment-band fixture
  // (ra = 1.2, rb = 0.8, d = 1.4) but placed along +x on a non-zero origin. The
  // reference is the cap formula, never the implementation: each cap height
  // follows from the centre distance, and the lens is their sum.
  const double radius_a = 1.2;
  const double radius_b = 0.8;
  const double distance = 1.4;
  const double height_a =
      radius_a -
      (distance * distance + radius_a * radius_a - radius_b * radius_b) / (2.0 * distance);
  const double height_b =
      radius_b -
      (distance * distance + radius_b * radius_b - radius_a * radius_a) / (2.0 * distance);
  const auto cap_volume = [](double radius, double height) {
    return kPi * height * height * (3.0 * radius - height) / 3.0;
  };
  // ha = 3/14, hb = 27/70: 0.21428571428571441 / 0.38571428571428584, caps
  // 0.16280410908909224 / 0.31382037635553112, lens 0.47662448544462332.
  const double lens_reference = cap_volume(radius_a, height_a) + cap_volume(radius_b, height_b);
  const geom::Point3d center_a{0.3, -0.2, 0.1};
  const geom::Point3d center_b{center_a[0] + distance, center_a[1], center_a[2]};
  const geom::Brep general_a = geom::make_sphere(center_a, radius_a, tolerance).value();
  const geom::Brep general_b = geom::make_sphere(center_b, radius_b, tolerance).value();
  check_lens_row("general", general_a, general_b, geom::SolidBooleanOp::intersection,
                 lens_reference, 1e-9);
}

// M3c-3c-1 general-placement acceptance corpus. The M3c-3b re-chart aligns each
// whole sphere's chart with the centre line, so a lens whose centre line points
// anywhere in the z = 0 plane must satisfy the same placement, volume and
// topology assertions as the +x fixture it is a rotation of. The operation count
// stays inside the 33-operation Debug budget set by the plan: 12 + 6 + 10 + 5.
namespace {

const char* boolean_op_name(geom::SolidBooleanOp op) {
  switch (op) {
    case geom::SolidBooleanOp::intersection:
      return "intersection";
    case geom::SolidBooleanOp::union_:
      return "union";
    case geom::SolidBooleanOp::difference:
      return "difference";
  }
  return "unknown";
}

// Two unit spheres at centre distance d meet in two spherical caps of height
// h = 1 - d / 2: the lens is 2 * pi * h^2 * (3 - h) / 3, the union is two full
// spheres minus the lens, and subtracting the second operand leaves one full
// sphere minus the lens. The cap-sum reference is analytic, never the
// implementation.
double rotated_lens_volume(double distance, geom::SolidBooleanOp op) {
  const double height = 1.0 - 0.5 * distance;
  const double lens = 2.0 * kPi * height * height * (3.0 - height) / 3.0;
  switch (op) {
    case geom::SolidBooleanOp::intersection:
      return lens;
    case geom::SolidBooleanOp::union_:
      return 8.0 * kPi / 3.0 - lens;
    case geom::SolidBooleanOp::difference:
      return 4.0 * kPi / 3.0 - lens;
  }
  return 0.0;
}

// One rotated placement: the first operand is the unit sphere at the origin and
// the second is the unit sphere centred at
// `distance * (cos(degrees), sin(degrees), 0)`, so `degrees` is the direction of
// the centre line. Operands stay ordered south then north, which makes
// `difference` the south complement (39*pi/32 at distance 1.5).
void check_rotated_lens(double degrees, double distance, geom::SolidBooleanOp op,
                        double relative_tolerance) {
  INFO("degrees " << degrees << " distance " << distance << " op " << boolean_op_name(op));
  const double angle = degrees * kPi / 180.0;
  const geom::Result<geom::Brep> south =
      geom::make_sphere(geom::Point3d{0.0, 0.0, 0.0}, 1.0, 1e-9);
  REQUIRE(south.has_value());
  const geom::Result<geom::Brep> north = geom::make_sphere(
      geom::Point3d{distance * std::cos(angle), distance * std::sin(angle), 0.0}, 1.0, 1e-9);
  REQUIRE(north.has_value());

  const geom::Result<geom::CurvedBooleanResult> result = geom::boolean_solid_curved(
      south.value(), south.value().solid_ids().front(), north.value(),
      north.value().solid_ids().front(), op, geom::CurvedBooleanOptions{1e-9});
  REQUIRE_MESSAGE(result.has_value(), "degrees " << degrees << " distance " << distance << " op "
                                                << boolean_op_name(op) << ": "
                                                << message_of(result.error()));
  CHECK(geom::validate_topology(result.value().brep).ok());
  check_one_shared_rim(boolean_op_name(op), result.value().brep);

  const std::optional<double> volume = sphere_face_volume(result.value().brep);
  REQUIRE_MESSAGE(volume.has_value(), "a face is not a polar cap");
  const double expected = rotated_lens_volume(distance, op);
  CHECK_MESSAGE(*volume == doctest::Approx(expected).epsilon(relative_tolerance),
                "volume " << *volume << " expected " << expected);
}

}  // namespace

TEST_CASE("curved boolean: sphere lenses accept rotated placements") {
  // Intersection at every 30 degrees covers all four quadrants at the 1.5
  // distance. The known-hard 15-degree placement stays pinned by the existing
  // off-axis case above, not by this 30/10-degree grid.
  for (int degrees = 0; degrees < 360; degrees += 30) {
    check_rotated_lens(degrees, 1.5, geom::SolidBooleanOp::intersection, 1e-9);
  }

  // Union and difference take the other two classification paths through the
  // same cut, so they are pinned at the axis, the 45-degree diagonal and the +y
  // axis.
  for (const int degrees : {0, 45, 90}) {
    check_rotated_lens(degrees, 1.5, geom::SolidBooleanOp::union_, 1e-6);
    check_rotated_lens(degrees, 1.5, geom::SolidBooleanOp::difference, 1e-6);
  }

  // Fine sweep over the first quadrant at the fixed 1.5 distance.
  for (int degrees = 0; degrees <= 90; degrees += 10) {
    check_rotated_lens(degrees, 1.5, geom::SolidBooleanOp::intersection, 1e-9);
  }

  // Five deterministic draws over the whole placement circle and the whole
  // crossing distance band (1.05..1.95 for unit spheres). The seed pins the
  // sequence so a failure names the exact placement above it.
  std::mt19937 generator{20261116u};
  std::uniform_real_distribution<double> degree_draw{0.0, 360.0};
  std::uniform_real_distribution<double> distance_draw{1.05, 1.95};
  for (int draw = 0; draw < 5; ++draw) {
    const double degrees = degree_draw(generator);
    const double distance = distance_draw(generator);
    INFO("random draw " << draw);
    check_rotated_lens(degrees, distance, geom::SolidBooleanOp::intersection, 1e-9);
  }
}
