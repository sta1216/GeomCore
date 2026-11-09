#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <span>
#include <string_view>
#include <utility>
#include <vector>

#include <doctest/doctest.h>
#include <geom/brep/boolean_solid.hpp>
#include <geom/brep/face_query.hpp>
#include <geom/brep/primitives.hpp>
#include <geom/brep/solid_editor.hpp>
#include <geom/brep/topology.hpp>
#include <geom/brep/validate.hpp>
#include <geom/core/error.hpp>

using namespace geom;

namespace {

constexpr std::string_view kCurvedMessage =
    "boolean_solid: curved faces are deferred to M3b";
constexpr std::string_view kMultiShellMessage =
    "boolean_solid: multi-shell operands are deferred to M3c";
constexpr std::string_view kExactlyOneMessage =
    "boolean_solid: expected exactly one solid per operand";
constexpr std::string_view kToleranceMessage =
    "boolean_solid: tolerance must be positive and finite";
constexpr std::string_view kEmptyResultMessage =
    "boolean_solid: result is empty";
constexpr std::string_view kZeroVolumeMessage =
    "boolean_solid: operand has zero volume";

template <class T>
T take(const Result<T>& result) {
  REQUIRE(result.has_value());
  return result.value();
}

Brep box(Point3d lo, Point3d hi) {
  return take(make_box(lo, hi, 1e-12));
}

const char* op_name(SolidBooleanOp op) {
  switch (op) {
    case SolidBooleanOp::union_:
      return "union";
    case SolidBooleanOp::intersection:
      return "intersection";
    case SolidBooleanOp::difference:
      return "difference";
    case SolidBooleanOp::symmetric_difference:
      return "symmetric_difference";
  }
  return "unknown";
}

// A deep copy of the topology counts and vertex positions is enough to prove a
// failed boolean left both operands exactly as they were. The Brep API is
// read-only, so this also catches accidental edits made through a future cache
// or editor path.
struct BrepSnapshot {
  std::uint32_t vertices = 0;
  std::uint32_t edges = 0;
  std::uint32_t coedges = 0;
  std::uint32_t loops = 0;
  std::uint32_t faces = 0;
  std::uint32_t shells = 0;
  std::uint32_t solids = 0;
  std::size_t curves = 0;
  std::size_t surfaces = 0;
  std::size_t pcurves = 0;
  std::vector<Point3d> points;
};

BrepSnapshot snapshot(const Brep& brep) {
  BrepSnapshot result;
  result.vertices = brep.vertex_count();
  result.edges = brep.edge_count();
  result.coedges = brep.coedge_count();
  result.loops = brep.loop_count();
  result.faces = brep.face_count();
  result.shells = brep.shell_count();
  result.solids = brep.solid_count();
  result.curves = brep.geometry().curves.size();
  result.surfaces = brep.geometry().surfaces.size();
  result.pcurves = brep.geometry().pcurves.size();
  for (const VertexId id : brep.vertex_ids()) {
    result.points.push_back(brep.vertex(id).point);
  }
  return result;
}

void check_unchanged(const Brep& brep, const BrepSnapshot& before) {
  CHECK(brep.vertex_count() == before.vertices);
  CHECK(brep.edge_count() == before.edges);
  CHECK(brep.coedge_count() == before.coedges);
  CHECK(brep.loop_count() == before.loops);
  CHECK(brep.face_count() == before.faces);
  CHECK(brep.shell_count() == before.shells);
  CHECK(brep.solid_count() == before.solids);
  CHECK(brep.geometry().curves.size() == before.curves);
  CHECK(brep.geometry().surfaces.size() == before.surfaces);
  CHECK(brep.geometry().pcurves.size() == before.pcurves);
  const std::vector<VertexId> vertices = brep.vertex_ids();
  REQUIRE(vertices.size() == before.points.size());
  for (std::size_t i = 0; i < vertices.size(); ++i) {
    const Point3d& point = brep.vertex(vertices[i]).point;
    for (int axis = 0; axis < 3; ++axis) {
      CHECK(point[axis] == before.points[i][axis]);
    }
  }
}

void check_error(const Result<SolidBooleanResult>& result, Errc code,
                 std::string_view message, std::uint32_t context_a = 0,
                 std::uint32_t context_b = 0) {
  REQUIRE_FALSE(result.has_value());
  CHECK(result.error().code == code);
  CHECK(std::string_view(result.error().message) == message);
  CHECK(result.error().context_a == context_a);
  CHECK(result.error().context_b == context_b);
}

Brep make_flat_box(Point3d lo, Point3d hi) {
  SolidEditor editor;
  const ShellId shell = [&editor, lo, hi]() {
    const Point3d corners[8] = {
        {lo[0], lo[1], lo[2]}, {hi[0], lo[1], lo[2]}, {hi[0], hi[1], lo[2]},
        {lo[0], hi[1], lo[2]}, {lo[0], lo[1], hi[2]}, {hi[0], lo[1], hi[2]},
        {hi[0], hi[1], hi[2]}, {lo[0], hi[1], hi[2]},
    };
    std::array<VertexId, 8> vertices{};
    for (int i = 0; i < 8; ++i) {
      vertices[i] = take(editor.add_vertex(corners[i], 0.0));
    }
    const int edge_specs[12][2] = {{0, 1}, {1, 2},  {2, 3},  {3, 0},
                                   {4, 5}, {5, 6},  {6, 7},  {7, 4},
                                   {1, 5}, {2, 6},  {3, 7},  {0, 4}};
    std::array<EdgeId, 12> edges{};
    for (int i = 0; i < 12; ++i) {
      edges[i] = take(editor.add_edge(vertices[edge_specs[i][0]],
                                      vertices[edge_specs[i][1]], kNoGeometry,
                                      0.0, 1.0, 0.0));
    }
    const Vec3d normals[6] = {{0, 0, -1}, {0, 0, 1}, {0, -1, 0},
                              {0, 1, 0},  {-1, 0, 0}, {1, 0, 0}};
    const int corner_index[6] = {0, 4, 0, 2, 0, 1};
    const int face_edges[6][4] = {{3, 2, 1, 0},  {4, 5, 6, 7},
                                  {0, 8, 4, 11}, {2, 10, 6, 9},
                                  {11, 7, 10, 3}, {1, 9, 5, 8}};
    const bool face_reversed[6][4] = {{true, true, true, true},
                                      {false, false, false, false},
                                      {false, false, true, true},
                                      {false, false, true, true},
                                      {false, true, true, false},
                                      {false, false, true, true}};
    std::array<FaceId, 6> faces{};
    for (int face = 0; face < 6; ++face) {
      Vec3d normal = normals[face];
      const Point3d& corner = corners[corner_index[face]];
      const double offset =
          dot(normal, Vec3d{corner[0], corner[1], corner[2]});
      const std::uint32_t surface =
          editor.add_surface(Plane3{normal, offset});
      faces[face] = take(editor.add_face(surface, false, 0.0));
      std::array<CoedgeSpec, 4> specs{};
      for (int j = 0; j < 4; ++j) {
        specs[j] = CoedgeSpec{edges[face_edges[face][j]],
                              face_reversed[face][j]};
      }
      take(editor.add_loop(faces[face], specs, true));
    }
    return take(editor.add_shell(faces, true));
  }();
  take(editor.add_solid(std::array<ShellId, 1>{shell}, 1e-12));
  return take(editor.commit());
}

// A topology-only multi-shell solid with one genuinely curved surface. The
// validator refuses this shape (multi-shell orientation is not computable), so
// the staged view is used instead of commit(). boolean_solid must report the
// multi-shell shape error before it ever inspects the surfaces, which is the
// precedence this fixture exists to pin.
Brep make_multi_shell_with_curved_face() {
  SolidEditor editor;
  const std::uint32_t curved = editor.add_surface(
      Cylinder3{Point3d{0, 0, 0}, Vec3d{0, 0, 1}, Vec3d{1, 0, 0}, 1.0});
  const FaceId curved_face = take(editor.add_face(curved, false, 0.0));
  const FaceId planar_face = take(editor.add_face(kNoGeometry, false, 0.0));
  const std::array<FaceId, 1> curved_faces{curved_face};
  const std::array<FaceId, 1> planar_faces{planar_face};
  const ShellId first = take(editor.add_shell(curved_faces, true));
  const ShellId second = take(editor.add_shell(planar_faces, true));
  take(editor.add_solid(std::array<ShellId, 2>{first, second}, 0.0));
  return editor.brep();
}

std::vector<FaceId> non_planar_faces(const Brep& brep) {
  std::vector<FaceId> faces;
  for (const FaceId id : brep.face_ids()) {
    if (!face_plane(brep, id).has_value()) faces.push_back(id);
  }
  return faces;
}

}  // namespace

TEST_CASE("boolean_solid errors: curved operands are diagnosed on either side") {
  const Brep unit = box({0, 0, 0}, {1, 1, 1});
  const Brep cylinder =
      take(make_cylinder(Point3d{0, 0, 0}, Vec3d{0, 0, 1}, 0.5, 1.0, 1e-12));
  const std::vector<FaceId> curved = non_planar_faces(cylinder);
  REQUIRE(!curved.empty());
  const std::uint32_t curved_face = curved.front().index();
  for (const SolidBooleanOp op :
       {SolidBooleanOp::union_, SolidBooleanOp::intersection,
        SolidBooleanOp::difference, SolidBooleanOp::symmetric_difference}) {
    for (int side = 0; side < 2; ++side) {
      CAPTURE(op_name(op));
      CAPTURE(side);
      const Brep& a = side == 0 ? unit : cylinder;
      const Brep& b = side == 0 ? cylinder : unit;
      const BrepSnapshot before_a = snapshot(a);
      const BrepSnapshot before_b = snapshot(b);
      const auto result = boolean_solid(a, b, op);
      check_error(result, Errc::unsupported, kCurvedMessage,
                  side == 1 ? curved_face : 0u,
                  side == 0 ? curved_face : 0u);
      check_unchanged(a, before_a);
      check_unchanged(b, before_b);
    }
  }
}

TEST_CASE("boolean_solid errors: multi-shell operands are diagnosed on either side") {
  const Brep split = make_multi_shell_with_curved_face();
  REQUIRE(split.solid_count() == 1);
  const std::vector<SolidId> solids = split.solid_ids();
  REQUIRE(solids.size() == 1);
  REQUIRE(split.solid_shells(solids.front()).size() == 2);
  const std::uint32_t solid = solids.front().index();
  const Brep operand = box({0, 0, 0}, {1, 1, 1});
  for (int side = 0; side < 2; ++side) {
    CAPTURE(side);
    const Brep& a = side == 0 ? split : operand;
    const Brep& b = side == 0 ? operand : split;
    const BrepSnapshot before_a = snapshot(a);
    const BrepSnapshot before_b = snapshot(b);
    const auto result = boolean_solid(a, b, SolidBooleanOp::union_);
    check_error(result, Errc::unsupported, kMultiShellMessage,
                side == 0 ? solid : 0u, side == 1 ? solid : 0u);
    check_unchanged(a, before_a);
    check_unchanged(b, before_b);
  }
}

TEST_CASE("boolean_solid errors: an empty Brep is an invalid operand document") {
  // An empty Brep (no solids) is a valid empty document, not a malformed one.
  // boolean_solid still requires exactly one solid per operand, so it reports
  // the operand-shape error instead of treating the empty input as an empty
  // result.
  const Brep empty;
  REQUIRE(validate_topology(empty).ok());
  REQUIRE(empty.solid_count() == 0);
  const Brep operand = box({0, 0, 0}, {1, 1, 1});
  for (int side = 0; side < 2; ++side) {
    CAPTURE(side);
    const Brep& a = side == 0 ? empty : operand;
    const Brep& b = side == 0 ? operand : empty;
    const BrepSnapshot before_a = snapshot(a);
    const BrepSnapshot before_b = snapshot(b);
    const auto result = boolean_solid(a, b, SolidBooleanOp::union_);
    check_error(result, Errc::unsupported, kExactlyOneMessage);
    check_unchanged(a, before_a);
    check_unchanged(b, before_b);
  }
}

TEST_CASE("boolean_solid errors: subtracting an identical box is empty") {
  const Brep operand = box({0, 0, 0}, {1, 1, 1});
  const BrepSnapshot before = snapshot(operand);
  const auto result =
      boolean_solid(operand, operand, SolidBooleanOp::difference);
  check_error(result, Errc::degenerate, kEmptyResultMessage);
  check_unchanged(operand, before);
}

TEST_CASE("boolean_solid errors: edge-touching boxes are not silently empty") {
  const Brep first = box({0, 0, 0}, {1, 1, 1});
  const Brep second = box({1, 1, 0}, {2, 2, 1});
  for (int side = 0; side < 2; ++side) {
    CAPTURE(side);
    const Brep& a = side == 0 ? first : second;
    const Brep& b = side == 0 ? second : first;
    const BrepSnapshot before_a = snapshot(a);
    const BrepSnapshot before_b = snapshot(b);
    const auto result = boolean_solid(a, b, SolidBooleanOp::union_);
    REQUIRE_FALSE(result.has_value());
    const Error& error = result.error();
    CHECK((error.code == Errc::not_closed ||
           error.code == Errc::invalid_topology));
    CHECK_FALSE(std::string_view(error.message).empty());
    check_unchanged(a, before_a);
    check_unchanged(b, before_b);
  }
}

TEST_CASE("boolean_solid errors: every invalid tolerance uses tolerance_violation") {
  const Brep a = box({0, 0, 0}, {1, 1, 1});
  const Brep b = box({0.5, 0.5, 0.5}, {1.5, 1.5, 1.5});
  const BrepSnapshot before_a = snapshot(a);
  const BrepSnapshot before_b = snapshot(b);
  const std::array<double, 5> invalid = {
      0.0,
      -1e-9,
      std::numeric_limits<double>::quiet_NaN(),
      std::numeric_limits<double>::infinity(),
      -std::numeric_limits<double>::infinity(),
  };
  for (const double tolerance : invalid) {
    CAPTURE(tolerance);
    SolidBooleanOptions options;
    options.tolerance = tolerance;
    const auto result =
        boolean_solid(a, b, SolidBooleanOp::union_, options);
    check_error(result, Errc::tolerance_violation, kToleranceMessage);
  }
  check_unchanged(a, before_a);
  check_unchanged(b, before_b);
}

TEST_CASE("boolean_solid errors: zero-height operands are diagnosed") {
  const Brep flat = make_flat_box({0, 0, 0}, {1, 1, 0});
  REQUIRE(validate_topology(flat).ok());
  REQUIRE(validate_tolerances(flat).ok());
  const std::vector<SolidId> solids = flat.solid_ids();
  REQUIRE(solids.size() == 1);
  const std::span<const ShellId> shells = flat.solid_shells(solids.front());
  REQUIRE(shells.size() == 1);
  const std::uint32_t shell = shells.front().index();
  const Brep operand = box({0, 0, 0}, {1, 1, 1});
  for (int side = 0; side < 2; ++side) {
    CAPTURE(side);
    const Brep& a = side == 0 ? flat : operand;
    const Brep& b = side == 0 ? operand : flat;
    const BrepSnapshot before_a = snapshot(a);
    const BrepSnapshot before_b = snapshot(b);
    const auto result = boolean_solid(a, b, SolidBooleanOp::union_);
    check_error(result, Errc::degenerate, kZeroVolumeMessage,
                side == 0 ? shell : 0u, side == 1 ? shell : 0u);
    check_unchanged(a, before_a);
    check_unchanged(b, before_b);
  }
}

TEST_CASE("boolean_solid errors: precedence is options, shape, then surfaces") {
  const Brep unit = box({0, 0, 0}, {1, 1, 1});
  const Brep cylinder =
      take(make_cylinder(Point3d{0, 0, 0}, Vec3d{0, 0, 1}, 0.5, 1.0, 1e-12));
  const Brep multi = make_multi_shell_with_curved_face();
  const std::vector<SolidId> multi_solids = multi.solid_ids();
  REQUIRE(multi_solids.size() == 1);
  const std::uint32_t multi_solid = multi_solids.front().index();

  SUBCASE("an invalid tolerance is reported before a curved operand") {
    SolidBooleanOptions options;
    options.tolerance = 0.0;
    for (int side = 0; side < 2; ++side) {
      CAPTURE(side);
      const Brep& a = side == 0 ? unit : cylinder;
      const Brep& b = side == 0 ? cylinder : unit;
      const auto result =
          boolean_solid(a, b, SolidBooleanOp::union_, options);
      check_error(result, Errc::tolerance_violation, kToleranceMessage);
    }
  }

  SUBCASE("the multi-shell shape is reported before a curved face") {
    for (int side = 0; side < 2; ++side) {
      CAPTURE(side);
      const Brep& a = side == 0 ? multi : unit;
      const Brep& b = side == 0 ? unit : multi;
      const auto result = boolean_solid(a, b, SolidBooleanOp::union_);
      check_error(result, Errc::unsupported, kMultiShellMessage,
                  side == 0 ? multi_solid : 0u,
                  side == 1 ? multi_solid : 0u);
    }
  }

  SUBCASE("an invalid tolerance is reported before a multi-shell shape") {
    SolidBooleanOptions options;
    options.tolerance = -1.0;
    const auto result =
        boolean_solid(multi, unit, SolidBooleanOp::union_, options);
    check_error(result, Errc::tolerance_violation, kToleranceMessage);
  }
}

TEST_CASE("boolean_solid: overlapping union report is deterministic") {
  const Brep a = box({0, 0, 0}, {2, 2, 2});
  const Brep b = box({1, 1, 1}, {3, 3, 3});
  const auto first = boolean_solid(a, b, SolidBooleanOp::union_);
  const auto second = boolean_solid(a, b, SolidBooleanOp::union_);
  REQUIRE(first.has_value());
  REQUIRE(second.has_value());
  const SolidBooleanReport& one = first.value().report;
  const SolidBooleanReport& two = second.value().report;
  CHECK(one.face_pairs == two.face_pairs);
  CHECK(one.cut_segments == two.cut_segments);
  CHECK(one.fragments == two.fragments);
  CHECK(one.kept == two.kept);
  CHECK(one.dropped == two.dropped);
  CHECK(one.flipped == two.flipped);
  CHECK(one.max_cut_residual == two.max_cut_residual);
  CHECK(one.face_pairs > 0);
  CHECK(one.cut_segments > 0);
  CHECK(one.fragments > 0);
  CHECK(one.kept > 0);
  CHECK(one.dropped > 0);
  CHECK(std::isfinite(one.max_cut_residual));
  CHECK(one.max_cut_residual <= SolidBooleanOptions{}.tolerance);
}
