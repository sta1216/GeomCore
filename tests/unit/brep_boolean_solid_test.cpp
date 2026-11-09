#include <array>
#include <cmath>
#include <cstdint>
#include <limits>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

#include <doctest/doctest.h>
#include <geom/brep/boolean_solid.hpp>
#include <geom/brep/face_query.hpp>
#include <geom/brep/inside.hpp>
#include <geom/brep/mass_properties.hpp>
#include <geom/brep/plane_frame.hpp>
#include <geom/brep/primitives.hpp>
#include <geom/brep/solid_editor.hpp>
#include <geom/brep/validate.hpp>
#include <geom/core/error.hpp>

using namespace geom;

namespace {

template <class T>
T take(const Result<T>& result) {
  REQUIRE(result.has_value());
  return result.value();
}

Brep box(Point3d lo, Point3d hi) {
  return take(make_box(lo, hi, 1e-12));
}

double volume_of(const Brep& brep) {
  return take(mass_properties(brep)).volume;
}

void check_valid(const Brep& brep) {
  CHECK(validate_topology(brep).ok());
  CHECK(validate_tolerances(brep).ok());
}

// A topology-only box with the six faces split across two closed shells. This
// is the multi-shell operand shape that the current validator defers. Only the
// topology is built: the single shared plane is never used for point tests.
Brep make_split_shell_box() {
  SolidEditor ed;
  const Point3d pts[8] = {{0, 0, 0}, {1, 0, 0}, {1, 1, 0}, {0, 1, 0},
                          {0, 0, 1}, {1, 0, 1}, {1, 1, 1}, {0, 1, 1}};
  std::array<VertexId, 8> vertices{};
  for (int i = 0; i < 8; ++i) vertices[i] = take(ed.add_vertex(pts[i], 0.0));
  const int pairs[12][2] = {{0, 1}, {1, 2}, {2, 3}, {3, 0}, {4, 5}, {5, 6},
                            {6, 7}, {7, 4}, {0, 4}, {1, 5}, {2, 6}, {3, 7}};
  std::array<EdgeId, 12> edges{};
  for (int i = 0; i < 12; ++i) {
    edges[i] = take(ed.add_edge(vertices[pairs[i][0]], vertices[pairs[i][1]],
                                kNoGeometry, 0.0, 1.0, 0.0));
  }
  const std::uint32_t plane = ed.add_surface(Plane3{{0, 0, 1}, 0.0});
  const int face_edges[6][4] = {{0, 1, 2, 3},  {7, 6, 5, 4},  {8, 4, 9, 0},
                                {9, 5, 10, 1}, {10, 6, 11, 2}, {11, 7, 8, 3}};
  const bool face_reversed[6][4] = {{false, false, false, false},
                                    {true, true, true, true},
                                    {false, false, true, true},
                                    {false, false, true, true},
                                    {false, false, true, true},
                                    {false, false, true, true}};
  std::array<FaceId, 6> faces{};
  for (int i = 0; i < 6; ++i) {
    faces[i] = take(ed.add_face(plane, false, 0.0));
    std::array<CoedgeSpec, 4> specs{};
    for (int j = 0; j < 4; ++j) {
      specs[j] = CoedgeSpec{edges[face_edges[i][j]], face_reversed[i][j]};
    }
    take(ed.add_loop(faces[i], specs, true));
  }
  const ShellId first =
      take(ed.add_shell(std::array<FaceId, 3>{faces[0], faces[1], faces[2]},
                        true));
  const ShellId second =
      take(ed.add_shell(std::array<FaceId, 3>{faces[3], faces[4], faces[5]},
                        true));
  take(ed.add_solid(std::array<ShellId, 2>{first, second}, 0.0));
  return ed.brep();
}

void check_case(const Brep& a, const Brep& b, SolidBooleanOp op,
                double expected_volume, std::uint32_t expected_solids = 1) {
  const auto result = boolean_solid(a, b, op);
  REQUIRE(result.has_value());
  const Brep& out = result.value().brep;
  check_valid(out);
  CHECK(out.solid_count() == expected_solids);
  CHECK(volume_of(out) == doctest::Approx(expected_volume).epsilon(1e-9));
}

// The symmetric difference of partially overlapping boxes removes a cube that
// touches the outer boundary, which makes the regularized result edge-connected
// and therefore non-manifold. The contract is a diagnosed error, not a
// malformed Brep.
void check_symdiff_diagnosed(const Brep& a, const Brep& b) {
  const auto result = boolean_solid(a, b, SolidBooleanOp::symmetric_difference);
  REQUIRE_FALSE(result.has_value());
  CHECK(result.error().code == Errc::invalid_topology);
}

}  // namespace

TEST_CASE("boolean_solid: overlapping boxes, four operations") {
  const Brep a = make_box(Point3d{0, 0, 0}, Point3d{2, 2, 2}, 1e-12).value();
  const Brep b = make_box(Point3d{1, 1, 1}, Point3d{3, 3, 3}, 1e-12).value();
  struct Case {
    SolidBooleanOp op;
    double volume;
  };
  for (const Case c : {Case{SolidBooleanOp::union_, 15.0},
                       Case{SolidBooleanOp::intersection, 1.0},
                       Case{SolidBooleanOp::difference, 7.0}}) {
    const auto result = boolean_solid(a, b, c.op);
    REQUIRE(result.has_value());
    const Brep& out = result.value().brep;
    CHECK(validate_topology(out).ok());
    CHECK(validate_tolerances(out).ok());
    CHECK(mass_properties(out).value().volume ==
          doctest::Approx(c.volume).epsilon(1e-9));
    CHECK(out.solid_count() == 1);
  }
  // Pinned contract: the symmetric difference (volume 14) is edge-connected.
  // The removed overlap cube is interior to the union, but its walls meet the
  // union boundary along edges, so four faces would share those edges. A single
  // manifold B-rep cannot represent that; assemble_solid rejects it and
  // boolean_solid propagates the diagnosis instead of returning a malformed
  // Brep. Task 9 owns richer diagnostics for edge-touching results.
  const auto symdiff = boolean_solid(a, b, SolidBooleanOp::symmetric_difference);
  REQUIRE_FALSE(symdiff.has_value());
  CHECK(symdiff.error().code == Errc::invalid_topology);
  CHECK(std::string_view(symdiff.error().message) ==
        "assemble_solid: non-manifold edge");
}

TEST_CASE("boolean_solid: stacked boxes merge their shared face away") {
  const Brep a = make_box(Point3d{0, 0, 0}, Point3d{1, 1, 1}, 1e-12).value();
  const Brep b = make_box(Point3d{0, 0, 1}, Point3d{1, 1, 2}, 1e-12).value();
  const auto result = boolean_solid(a, b, SolidBooleanOp::union_);
  REQUIRE(result.has_value());
  CHECK(mass_properties(result.value().brep).value().volume ==
        doctest::Approx(2.0).epsilon(1e-9));
  CHECK(result.value().brep.shell_count() == 1);
  CHECK(result.value().brep.solid_count() == 1);
  check_valid(result.value().brep);
}

TEST_CASE("boolean_solid: stacked difference keeps only the lower box") {
  const Brep a = box({0, 0, 0}, {1, 1, 1});
  const Brep b = box({0, 0, 1}, {1, 1, 2});
  const auto result = boolean_solid(a, b, SolidBooleanOp::difference);
  REQUIRE(result.has_value());
  const Brep& out = result.value().brep;
  CHECK(out.solid_count() == 1);
  CHECK(out.shell_count() == 1);
  CHECK(out.face_count() == 6);
  CHECK(volume_of(out) == doctest::Approx(1.0).epsilon(1e-12));
  check_valid(out);
}

TEST_CASE("boolean_solid: partial overlaps follow the volume formula") {
  struct Boxes {
    Point3d a_lo;
    Point3d a_hi;
    Point3d b_lo;
    Point3d b_hi;
    double overlap;
  };
  const Boxes cases[3] = {
      {{0, 0, 0}, {2, 2, 2}, {0.5, 0.5, 0.5}, {1.5, 2.5, 1.5}, 1.5},
      {{0, 0, 0}, {1, 1, 1}, {0.25, 0.25, -0.5}, {1.75, 1.25, 0.75},
       0.421875},
      {{-0.5, 0, 0}, {1.5, 3, 2}, {0.25, 0.75, 0.5}, {2.75, 1.75, 3.5},
       1.875},
  };
  for (const Boxes& c : cases) {
    const Brep a = box(c.a_lo, c.a_hi);
    const Brep b = box(c.b_lo, c.b_hi);
    const double va = volume_of(a);
    const double vb = volume_of(b);
    const double vi = c.overlap;
    check_case(a, b, SolidBooleanOp::intersection, vi);
    check_case(a, b, SolidBooleanOp::union_, va + vb - vi);
    check_case(a, b, SolidBooleanOp::difference, va - vi);
    check_symdiff_diagnosed(a, b);
  }
}

TEST_CASE("boolean_solid: side-by-side boxes share only a face") {
  const Brep a = box({0, 0, 0}, {1, 1, 1});
  const Brep b = box({1, 0, 0}, {2, 1, 1});
  {
    const auto result = boolean_solid(a, b, SolidBooleanOp::union_);
    REQUIRE(result.has_value());
    check_valid(result.value().brep);
    CHECK(result.value().brep.solid_count() == 1);
    CHECK(volume_of(result.value().brep) ==
          doctest::Approx(2.0).epsilon(1e-9));
  }
  {
    // The intersection is the shared face: a zero-volume result. The documented
    // contract is a diagnosed empty result rather than a degenerate solid.
    const auto result = boolean_solid(a, b, SolidBooleanOp::intersection);
    REQUIRE_FALSE(result.has_value());
    CHECK(result.error().code == Errc::degenerate);
    CHECK(std::string_view(result.error().message) ==
          "boolean_solid: result is empty");
  }
  {
    const auto result = boolean_solid(a, b, SolidBooleanOp::difference);
    REQUIRE(result.has_value());
    check_valid(result.value().brep);
    CHECK(result.value().brep.solid_count() == 1);
    CHECK(volume_of(result.value().brep) ==
          doctest::Approx(1.0).epsilon(1e-9));
  }
}

TEST_CASE("boolean_solid: report counters describe the work") {
  const Brep a = box({0, 0, 0}, {2, 2, 2});
  const Brep b = box({1, 1, 1}, {3, 3, 3});
  const auto result = boolean_solid(a, b, SolidBooleanOp::difference);
  REQUIRE(result.has_value());
  const SolidBooleanReport& report = result.value().report;
  CHECK(report.face_pairs > 0);
  CHECK(report.cut_segments > 0);
  CHECK(report.fragments > 0);
  CHECK(report.kept > 0);
  CHECK(report.dropped > 0);
  CHECK(report.flipped > 0);
  CHECK(std::isfinite(report.max_cut_residual));
  CHECK(report.max_cut_residual <= 1e-9);
}

TEST_CASE("boolean_solid: disjoint boxes stay two solids") {
  const Brep a = box({0, 0, 0}, {1, 1, 1});
  const Brep b = box({2, 0, 0}, {3, 1, 1});
  const auto result = boolean_solid(a, b, SolidBooleanOp::union_);
  REQUIRE(result.has_value());
  const Brep& out = result.value().brep;
  check_valid(out);
  // Multi-solid output is supported (assemble_solid emits one solid per outward
  // component); the two touching-free boxes therefore stay separate.
  CHECK(out.solid_count() == 2);
  CHECK(volume_of(out) == doctest::Approx(2.0).epsilon(1e-12));
  CHECK(result.value().report.face_pairs == 0);
  CHECK(result.value().report.cut_segments == 0);
  CHECK(result.value().report.kept == 12);
  CHECK(result.value().report.dropped == 0);
}

TEST_CASE("boolean_solid errors: curved operands are unsupported") {
  const Brep box_brep = box({0, 0, 0}, {1, 1, 1});
  const Brep cylinder =
      take(make_cylinder(Point3d{0, 0, 0}, Vec3d{0, 0, 1}, 0.5, 1.0, 1e-12));
  const auto result = boolean_solid(box_brep, cylinder, SolidBooleanOp::union_);
  REQUIRE_FALSE(result.has_value());
  CHECK(result.error().code == Errc::unsupported);
}

TEST_CASE("boolean_solid errors: multi-shell operands are unsupported") {
  const Brep split = make_split_shell_box();
  REQUIRE_FALSE(validate_topology(split).ok());
  CHECK(validate_topology(split).issues.front().code == Errc::unsupported);
  const Brep operand = box({0, 0, 0}, {1, 1, 1});
  const auto result = boolean_solid(split, operand, SolidBooleanOp::union_);
  REQUIRE_FALSE(result.has_value());
  CHECK(result.error().code == Errc::unsupported);
}

TEST_CASE("boolean_solid errors: tolerances must be positive and finite") {
  const Brep a = box({0, 0, 0}, {1, 1, 1});
  const Brep b = box({0.5, 0.5, 0.5}, {1.5, 1.5, 1.5});
  for (const double tol : {0.0, -1e-9}) {
    SolidBooleanOptions options;
    options.tolerance = tol;
    const auto result = boolean_solid(a, b, SolidBooleanOp::union_, options);
    REQUIRE_FALSE(result.has_value());
    CHECK(result.error().code == Errc::degenerate);
  }
  for (const double tol : {std::numeric_limits<double>::quiet_NaN(),
                           std::numeric_limits<double>::infinity()}) {
    SolidBooleanOptions options;
    options.tolerance = tol;
    const auto result = boolean_solid(a, b, SolidBooleanOp::union_, options);
    REQUIRE_FALSE(result.has_value());
    CHECK(result.error().code == Errc::non_finite);
  }
}

TEST_CASE("boolean_solid errors: invalid topology is rejected") {
  SolidEditor ed;
  take(ed.add_face(kNoGeometry, false, 0.0));
  const Brep invalid = ed.brep();
  REQUIRE_FALSE(validate_topology(invalid).ok());
  const Brep operand = box({0, 0, 0}, {1, 1, 1});
  const auto result = boolean_solid(invalid, operand, SolidBooleanOp::union_);
  REQUIRE_FALSE(result.has_value());
  CHECK(result.error().code == Errc::invalid_topology);
}

TEST_CASE("boolean_solid: intersection face normals point outward") {
  const Brep a = box({0, 0, 0}, {2, 2, 2});
  const Brep b = box({1, 1, 1}, {3, 3, 3});
  const auto result = boolean_solid(a, b, SolidBooleanOp::intersection);
  REQUIRE(result.has_value());
  const Brep& out = result.value().brep;
  REQUIRE(out.solid_count() == 1);
  const SolidId solid = out.solid_ids().front();
  CHECK(volume_of(out) > 0.0);
  for (const FaceId face : out.face_ids()) {
    const Plane3 plane = take(face_plane(out, face));
    const PlaneFrame frame = take(make_plane_frame(plane));
    const std::vector<std::vector<Point2d>> loops =
        take(face_loops_2d(out, face, frame));
    REQUIRE(!loops.empty());
    Point3d centre{0, 0, 0};
    for (const Point2d& p : loops.front()) {
      const Point3d p3 = plane_frame_to_3d(frame, p);
      centre = Point3d{centre[0] + p3[0], centre[1] + p3[1],
                       centre[2] + p3[2]};
    }
    const double count = static_cast<double>(loops.front().size());
    centre = Point3d{centre[0] / count, centre[1] / count, centre[2] / count};
    const Vec3d& n = frame.normal;
    const Point3d outside{centre[0] + n[0] * 1e-3, centre[1] + n[1] * 1e-3,
                          centre[2] + n[2] * 1e-3};
    const Point3d inside{centre[0] - n[0] * 1e-3, centre[1] - n[1] * 1e-3,
                         centre[2] - n[2] * 1e-3};
    CHECK_FALSE(take(point_in_solid(out, solid, outside, 1e-9)));
    CHECK(take(point_in_solid(out, solid, inside, 1e-9)));
  }
}
