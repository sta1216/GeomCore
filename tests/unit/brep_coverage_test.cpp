// Negative-path and boundary coverage for the M2 B-rep modules.
//
// Every case here drives a public-API guard that the happy-path suites never
// reach: primitive validation returns, SolidEditor add_*/erase_*/split/join
// failures, validator findings on deliberately malformed topologies, and the
// tessellate/mass_properties unsupported-domain errors. Assertions check the
// exact Errc, the exact message and, where the API exposes it, that the editor
// state is unchanged afterwards.

#include <doctest/doctest.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <initializer_list>
#include <limits>
#include <map>
#include <numbers>
#include <span>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

#include <geom/brep/mass_properties.hpp>
#include <geom/brep/primitives.hpp>
#include <geom/brep/solid_editor.hpp>
#include <geom/brep/tessellate.hpp>
#include <geom/brep/validate.hpp>
#include <geom/core/error.hpp>
#include <geom/curves/analytic_curve.hpp>
#include <geom/curves/nurbs_curve.hpp>
#include <geom/surfaces/analytic.hpp>
#include <geom/surfaces/nurbs_surface.hpp>
#include <geom/surfaces/plane.hpp>

namespace {

constexpr double kNaN = std::numeric_limits<double>::quiet_NaN();
constexpr double kInf = std::numeric_limits<double>::infinity();
constexpr double kPi = std::numbers::pi_v<double>;
constexpr double kTwoPi = 2.0 * kPi;

geom::Point3d point_with_nan(int axis) {
  geom::Point3d p{0.0, 0.0, 0.0};
  p[axis] = kNaN;
  return p;
}

geom::Point3d point_with_inf(int axis) {
  geom::Point3d p{0.0, 0.0, 0.0};
  p[axis] = kInf;
  return p;
}

geom::Vec3d vec_with_nan(int axis) {
  geom::Vec3d v{0.0, 0.0, 0.0};
  v[axis] = kNaN;
  return v;
}

geom::Vec3d vec_with_inf(int axis) {
  geom::Vec3d v{0.0, 0.0, 0.0};
  v[axis] = kInf;
  return v;
}

template <typename T>
void expect_error(const geom::Result<T>& result, geom::Errc code, std::string_view message) {
  INFO("expected error: " << message);
  REQUIRE_FALSE(result.has_value());
  CHECK(result.error().code == code);
  CHECK(std::string_view(result.error().message) == message);
}

bool has_issue(const geom::ValidationReport& report, geom::Errc code, std::string_view message) {
  for (const geom::ValidationIssue& issue : report.issues) {
    if (issue.code == code && std::string_view(issue.message) == message) return true;
  }
  return false;
}

std::string_view first_message(const geom::ValidationReport& report) {
  REQUIRE_FALSE(report.issues.empty());
  return report.issues.front().message;
}

// --- shared topology builders -------------------------------------------------

struct Square {
  geom::SolidEditor editor;
  std::array<geom::VertexId, 4> vertices{};
  std::array<geom::EdgeId, 4> edges{};
  geom::FaceId face{};
  geom::LoopId loop{};
  geom::ShellId shell{};
  geom::SolidId solid{};
  std::uint32_t surface = geom::kNoGeometry;
};

struct SquareOptions {
  geom::Plane3 plane{{0.0, 0.0, 1.0}, 0.0};
  std::array<double, 4> vertex_tol{0.0, 0.0, 0.0, 0.0};
  double edge_tol = 0.0;
  double coedge_tol = 0.0;
  double face_tol = 0.0;
  double solid_tol = 0.0;
  bool outer = true;
  bool closed_shell = false;
  bool with_loop = true;
  bool with_shell = true;
  bool with_solid = true;
};

// Unit square in the z = 0 plane, counter-clockwise seen from +z. Each edge has
// its own unit line curve so the sampled geometry agrees with the vertices.
Square make_square(const SquareOptions& options = {}) {
  Square s;
  const geom::Point3d pts[4] = {{0.0, 0.0, 0.0}, {1.0, 0.0, 0.0}, {1.0, 1.0, 0.0}, {0.0, 1.0, 0.0}};
  for (int i = 0; i < 4; ++i) {
    s.vertices[i] =
        s.editor.add_vertex(pts[i], options.vertex_tol[static_cast<std::size_t>(i)]).value();
  }
  s.surface = s.editor.add_surface(options.plane);
  for (int i = 0; i < 4; ++i) {
    const geom::Point3d from = pts[i];
    const std::uint32_t curve = s.editor.add_curve(geom::LineCurve3{from, {1.0, 0.0, 0.0}});
    s.edges[i] =
        s.editor.add_edge(s.vertices[i], s.vertices[(i + 1) % 4], curve, 0.0, 1.0, options.edge_tol)
            .value();
  }
  s.face = s.editor.add_face(s.surface, false, options.face_tol).value();
  if (options.with_loop) {
    std::array<geom::CoedgeSpec, 4> specs{};
    for (int i = 0; i < 4; ++i) {
      specs[static_cast<std::size_t>(i)] =
          geom::CoedgeSpec{s.edges[static_cast<std::size_t>(i)], false, geom::kNoGeometry,
                           options.coedge_tol};
    }
    s.loop = s.editor.add_loop(s.face, specs, options.outer).value();
  }
  if (options.with_shell) {
    s.shell = s.editor
                  .add_shell(std::array<geom::FaceId, 1>{s.face}, options.closed_shell)
                  .value();
    if (options.with_solid) {
      s.solid = s.editor.add_solid(std::array<geom::ShellId, 1>{s.shell}, options.solid_tol).value();
    }
  }
  return s;
}

// Two vertices joined by one edge, used by self-cycles.
struct TwoEdgeFace {
  geom::SolidEditor editor;
  geom::FaceId face{};
  geom::LoopId loop{};
};

// A face whose single loop walks one edge forward and back, so the ring has only
// two distinct samples.
TwoEdgeFace make_two_edge_face() {
  TwoEdgeFace result;
  const geom::Point3d a{0.0, 0.0, 0.0};
  const geom::Point3d b{1.0, 0.0, 0.0};
  const geom::VertexId v0 = result.editor.add_vertex(a, 0.0).value();
  const geom::VertexId v1 = result.editor.add_vertex(b, 0.0).value();
  const std::uint32_t curve = result.editor.add_curve(geom::LineCurve3{a, {1.0, 0.0, 0.0}});
  const geom::EdgeId edge = result.editor.add_edge(v0, v1, curve, 0.0, 1.0, 0.0).value();
  result.face = result.editor.add_face(result.editor.add_surface(geom::Plane3{{0.0, 0.0, 1.0}, 0.0}),
                                       false, 0.0)
                    .value();
  const std::array<geom::CoedgeSpec, 2> specs{geom::CoedgeSpec{edge, false},
                                              geom::CoedgeSpec{edge, true}};
  result.loop = result.editor.add_loop(result.face, specs, true).value();
  return result;
}

struct Cube {
  geom::SolidEditor editor;
  std::array<geom::VertexId, 8> vertices{};
  std::vector<geom::EdgeId> edges;
  std::array<geom::FaceId, 6> faces{};
  geom::ShellId shell{};
  geom::SolidId solid{};
};

struct CubeOptions {
  bool reversed = false;
  bool closed_shell = true;
  std::array<double, 6> face_tol{0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
  double solid_tol = 0.0;
  double vertex_tol = 0.0;
  double edge_tol = 0.0;
};

// Unit cube with outward-facing planar quads, one closed shell and one solid.
Cube make_cube(const CubeOptions& options = {}) {
  Cube cube;
  const geom::Point3d pts[8] = {{0.0, 0.0, 0.0}, {1.0, 0.0, 0.0}, {1.0, 1.0, 0.0}, {0.0, 1.0, 0.0},
                                {0.0, 0.0, 1.0}, {1.0, 0.0, 1.0}, {1.0, 1.0, 1.0}, {0.0, 1.0, 1.0}};
  const int quads[6][4] = {{0, 3, 2, 1}, {4, 5, 6, 7}, {0, 1, 5, 4},
                           {1, 2, 6, 5}, {2, 3, 7, 6}, {3, 0, 4, 7}};
  const geom::Vec3d normals[6] = {{0.0, 0.0, -1.0}, {0.0, 0.0, 1.0}, {0.0, -1.0, 0.0},
                                  {1.0, 0.0, 0.0},  {0.0, 1.0, 0.0},  {-1.0, 0.0, 0.0}};
  const int reference[6] = {0, 4, 0, 1, 2, 3};
  for (int i = 0; i < 8; ++i) cube.vertices[static_cast<std::size_t>(i)] =
      cube.editor.add_vertex(pts[i], options.vertex_tol).value();
  std::map<std::pair<int, int>, geom::EdgeId> shared;
  for (int f = 0; f < 6; ++f) {
    const geom::Point3d corner = pts[reference[f]];
    const std::uint32_t surface = cube.editor.add_surface(
        geom::Plane3{normals[f], geom::dot(normals[f], geom::Vec3d{corner[0], corner[1], corner[2]})});
    cube.faces[static_cast<std::size_t>(f)] =
        cube.editor.add_face(surface, options.reversed, options.face_tol[static_cast<std::size_t>(f)])
            .value();
    std::vector<geom::CoedgeSpec> specs;
    specs.reserve(4);
    for (int i = 0; i < 4; ++i) {
      const int a = quads[f][i];
      const int b = quads[f][(i + 1) % 4];
      const std::pair<int, int> key{std::min(a, b), std::max(a, b)};
      auto found = shared.find(key);
      if (found == shared.end()) {
        const geom::Point3d from = pts[key.first];
        const geom::Point3d to = pts[key.second];
        const geom::Vec3d delta = to - from;
        const double length = std::hypot(delta[0], delta[1], delta[2]);
        const std::uint32_t curve = cube.editor.add_curve(
            geom::LineCurve3{from, {delta[0] / length, delta[1] / length, delta[2] / length}});
        const geom::EdgeId edge =
            cube.editor
                .add_edge(cube.vertices[static_cast<std::size_t>(key.first)],
                          cube.vertices[static_cast<std::size_t>(key.second)], curve, 0.0, length,
                          options.edge_tol)
                .value();
        found = shared.emplace(key, edge).first;
        cube.edges.push_back(edge);
      }
      specs.push_back(geom::CoedgeSpec{found->second, a != key.first, geom::kNoGeometry, 0.0});
    }
    cube.editor.add_loop(cube.faces[static_cast<std::size_t>(f)], specs, true).value();
  }
  cube.shell = cube.editor
                   .add_shell(cube.faces, options.closed_shell)
                   .value();
  cube.solid = cube.editor.add_solid(std::array<geom::ShellId, 1>{cube.shell}, options.solid_tol).value();
  return cube;
}

geom::LoopId add_loop_of(geom::SolidEditor& editor, geom::FaceId face,
                         const std::vector<geom::CoedgeSpec>& specs, bool outer = true) {
  return editor.add_loop(face, std::span<const geom::CoedgeSpec>(specs), outer).value();
}

geom::EdgeId add_line_edge(geom::SolidEditor& editor, geom::Point3d from, geom::Point3d to,
                           double t0 = 0.0, double t1 = 1.0) {
  const geom::Vec3d delta = to - from;
  const double length = std::hypot(delta[0], delta[1], delta[2]);
  const std::uint32_t curve = editor.add_curve(
      geom::LineCurve3{from, {delta[0] / length, delta[1] / length, delta[2] / length}});
  const geom::VertexId v0 = editor.add_vertex(from, 0.0).value();
  const geom::VertexId v1 = editor.add_vertex(to, 0.0).value();
  return editor.add_edge(v0, v1, curve, t0, t1, 0.0).value();
}

}  // namespace

// =============================================================================
// primitives.hpp
// =============================================================================

TEST_CASE("primitive finiteness helpers short-circuit on each coordinate") {
  using geom::detail::primitive_point_equal;
  using geom::detail::primitive_point_finite;
  using geom::detail::primitive_vec_finite;
  CHECK(primitive_point_finite({0.0, 0.0, 0.0}));
  CHECK_FALSE(primitive_point_finite({kNaN, 0.0, 0.0}));
  CHECK_FALSE(primitive_point_finite({0.0, kInf, 0.0}));
  CHECK_FALSE(primitive_point_finite({0.0, 0.0, kNaN}));
  CHECK(primitive_vec_finite({0.0, 0.0, 0.0}));
  CHECK_FALSE(primitive_vec_finite({kNaN, 0.0, 0.0}));
  CHECK_FALSE(primitive_vec_finite({0.0, kInf, 0.0}));
  CHECK_FALSE(primitive_vec_finite({0.0, 0.0, kNaN}));
  CHECK(primitive_point_equal({1.0, 2.0, 3.0}, {1.0, 2.0, 3.0}));
  CHECK_FALSE(primitive_point_equal({0.0, 2.0, 3.0}, {1.0, 2.0, 3.0}));
  CHECK_FALSE(primitive_point_equal({1.0, 0.0, 3.0}, {1.0, 2.0, 3.0}));
  CHECK_FALSE(primitive_point_equal({1.0, 2.0, 0.0}, {1.0, 2.0, 3.0}));
}

TEST_CASE("make_box rejects non-finite corners, tolerances and degeneracy") {
  for (int axis = 0; axis < 3; ++axis) {
    expect_error(geom::make_box(point_with_nan(axis), {1.0, 1.0, 1.0}, 0.0), geom::Errc::non_finite,
                 "make_box: box corner is not finite");
    expect_error(geom::make_box({0.0, 0.0, 0.0}, point_with_inf(axis), 0.0),
                 geom::Errc::non_finite, "make_box: box corner is not finite");
  }
  expect_error(geom::make_box({0.0, 0.0, 0.0}, {1.0, 1.0, 1.0}, kNaN),
               geom::Errc::tolerance_violation,
               "make_box: tolerance must be finite and non-negative");
  expect_error(geom::make_box({0.0, 0.0, 0.0}, {1.0, 1.0, 1.0}, kInf),
               geom::Errc::tolerance_violation,
               "make_box: tolerance must be finite and non-negative");
  expect_error(geom::make_box({0.0, 0.0, 0.0}, {1.0, 1.0, 1.0}, -1.0),
               geom::Errc::tolerance_violation,
               "make_box: tolerance must be finite and non-negative");
  for (int axis = 0; axis < 3; ++axis) {
    geom::Point3d lo{0.0, 0.0, 0.0};
    geom::Point3d hi{1.0, 1.0, 1.0};
    lo[axis] = 1.0;
    expect_error(geom::make_box(lo, hi, 0.0), geom::Errc::degenerate,
                 "make_box: box is degenerate or inverted");
    lo[axis] = 2.0;
    expect_error(geom::make_box(lo, hi, 0.0), geom::Errc::degenerate,
                 "make_box: box is degenerate or inverted");
  }
  expect_error(geom::make_box({-1e308, -1e308, -1e308}, {1e308, 1e308, 1e308}, 0.0),
               geom::Errc::degenerate, "make_box: box edge length is not representable");
}

TEST_CASE("make_cylinder rejects every invalid input class") {
  const geom::Vec3d z{0.0, 0.0, 1.0};
  for (int axis = 0; axis < 3; ++axis) {
    expect_error(geom::make_cylinder(point_with_nan(axis), z, 1.0, 1.0, 0.0),
                 geom::Errc::non_finite, "make_cylinder: input is not finite");
    expect_error(geom::make_cylinder({0.0, 0.0, 0.0}, vec_with_nan(axis), 1.0, 1.0, 0.0),
                 geom::Errc::non_finite, "make_cylinder: input is not finite");
    expect_error(geom::make_cylinder({0.0, 0.0, 0.0}, vec_with_inf(axis), 1.0, 1.0, 0.0),
                 geom::Errc::non_finite, "make_cylinder: input is not finite");
  }
  expect_error(geom::make_cylinder({0.0, 0.0, 0.0}, z, kNaN, 1.0, 0.0), geom::Errc::non_finite,
               "make_cylinder: input is not finite");
  expect_error(geom::make_cylinder({0.0, 0.0, 0.0}, z, 1.0, kNaN, 0.0), geom::Errc::non_finite,
               "make_cylinder: input is not finite");
  expect_error(geom::make_cylinder({0.0, 0.0, 0.0}, z, 1.0, 1.0, kNaN),
               geom::Errc::tolerance_violation,
               "make_cylinder: tolerance must be finite and non-negative");
  expect_error(geom::make_cylinder({0.0, 0.0, 0.0}, z, 1.0, 1.0, -1.0),
               geom::Errc::tolerance_violation,
               "make_cylinder: tolerance must be finite and non-negative");
  expect_error(geom::make_cylinder({0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}, 1.0, 1.0, 0.0),
               geom::Errc::degenerate, "make_cylinder: axis must be unit length");
  expect_error(geom::make_cylinder({0.0, 0.0, 0.0}, {0.0, 0.0, 2.0}, 1.0, 1.0, 0.0),
               geom::Errc::degenerate, "make_cylinder: axis must be unit length");
  expect_error(geom::make_cylinder({0.0, 0.0, 0.0}, {1.0001, 0.0, 0.0}, 1.0, 1.0, 0.0),
               geom::Errc::degenerate, "make_cylinder: axis must be unit length");
  expect_error(geom::make_cylinder({0.0, 0.0, 0.0}, z, 0.0, 1.0, 0.0), geom::Errc::degenerate,
               "make_cylinder: radius must be positive");
  expect_error(geom::make_cylinder({0.0, 0.0, 0.0}, z, -1.0, 1.0, 0.0), geom::Errc::degenerate,
               "make_cylinder: radius must be positive");
  expect_error(geom::make_cylinder({0.0, 0.0, 0.0}, z, 1.0, 0.0, 0.0), geom::Errc::degenerate,
               "make_cylinder: height must be positive");
  expect_error(geom::make_cylinder({0.0, 0.0, 0.0}, z, 1.0, -1.0, 0.0), geom::Errc::degenerate,
               "make_cylinder: height must be positive");
  expect_error(geom::make_cylinder({1e308, 0.0, 0.0}, {1.0, 0.0, 0.0}, 1.0, 1e308, 0.0),
               geom::Errc::degenerate, "make_cylinder: top center is not representable");
  expect_error(geom::make_cylinder({1e308, 0.0, 0.0}, {1.0, 0.0, 0.0}, 1.0, 1.0, 0.0),
               geom::Errc::degenerate, "make_cylinder: seam has zero length");
  expect_error(geom::make_cylinder({1e308, 0.0, 0.0}, z, 1e308, 1.0, 0.0), geom::Errc::degenerate,
               "make_cylinder: seam point is not representable");
}

TEST_CASE("make_cone rejects every invalid input class") {
  const geom::Vec3d z{0.0, 0.0, 1.0};
  for (int axis = 0; axis < 3; ++axis) {
    expect_error(geom::make_cone(point_with_nan(axis), z, 1.0, 1.0, 0.0), geom::Errc::non_finite,
                 "make_cone: input is not finite");
    expect_error(geom::make_cone({0.0, 0.0, 0.0}, vec_with_nan(axis), 1.0, 1.0, 0.0),
                 geom::Errc::non_finite, "make_cone: input is not finite");
  }
  expect_error(geom::make_cone({0.0, 0.0, 0.0}, z, kNaN, 1.0, 0.0), geom::Errc::non_finite,
               "make_cone: input is not finite");
  expect_error(geom::make_cone({0.0, 0.0, 0.0}, z, 1.0, kInf, 0.0), geom::Errc::non_finite,
               "make_cone: input is not finite");
  expect_error(geom::make_cone({0.0, 0.0, 0.0}, z, 1.0, 1.0, -1.0), geom::Errc::tolerance_violation,
               "make_cone: tolerance must be finite and non-negative");
  expect_error(geom::make_cone({0.0, 0.0, 0.0}, z, 1.0, 1.0, kNaN), geom::Errc::tolerance_violation,
               "make_cone: tolerance must be finite and non-negative");
  expect_error(geom::make_cone({0.0, 0.0, 0.0}, {0.0, 1.0001, 0.0}, 1.0, 1.0, 0.0),
               geom::Errc::degenerate, "make_cone: axis must be unit length");
  expect_error(geom::make_cone({0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}, 1.0, 1.0, 0.0),
               geom::Errc::degenerate, "make_cone: axis must be unit length");
  expect_error(geom::make_cone({0.0, 0.0, 0.0}, z, 0.0, 1.0, 0.0), geom::Errc::degenerate,
               "make_cone: base radius must be positive");
  expect_error(geom::make_cone({0.0, 0.0, 0.0}, z, -1.0, 1.0, 0.0), geom::Errc::degenerate,
               "make_cone: base radius must be positive");
  expect_error(geom::make_cone({0.0, 0.0, 0.0}, z, 1.0, 0.0, 0.0), geom::Errc::degenerate,
               "make_cone: height must be positive");
  expect_error(geom::make_cone({0.0, 0.0, 0.0}, z, 1.0, -1.0, 0.0), geom::Errc::degenerate,
               "make_cone: height must be positive");
  expect_error(geom::make_cone({0.0, 0.0, 0.0}, z, 1e-300, 1e300, 0.0), geom::Errc::degenerate,
               "make_cone: half angle is not representable");
  expect_error(geom::make_cone({0.0, 0.0, 0.0}, z, 1e300, 1e-300, 0.0), geom::Errc::degenerate,
               "make_cone: half angle must be below pi/2");
  expect_error(geom::make_cone({1e308, 0.0, 0.0}, {1.0, 0.0, 0.0}, 1.0, 1e308, 0.0),
               geom::Errc::degenerate, "make_cone: apex is not representable");
  expect_error(geom::make_cone({1e308, 0.0, 0.0}, z, 9e307, 9e307, 0.0), geom::Errc::degenerate,
               "make_cone: seam point is not representable");
  const double sqrt_half = std::sqrt(0.5);
  expect_error(geom::make_cone({0.0, 0.0, 0.0}, {sqrt_half, sqrt_half, 0.0}, 1.5e308, 1.5e308, 0.0),
               geom::Errc::degenerate, "make_cone: seam is not representable");
}

TEST_CASE("make_sphere rejects every invalid input class") {
  for (int axis = 0; axis < 3; ++axis) {
    expect_error(geom::make_sphere(point_with_nan(axis), 1.0, 0.0), geom::Errc::non_finite,
                 "make_sphere: input is not finite");
    expect_error(geom::make_sphere({0.0, 0.0, 0.0}, axis == 0 ? kNaN : kInf, 0.0),
                 geom::Errc::non_finite, "make_sphere: input is not finite");
  }
  expect_error(geom::make_sphere({0.0, 0.0, 0.0}, 1.0, kNaN), geom::Errc::tolerance_violation,
               "make_sphere: tolerance must be finite and non-negative");
  expect_error(geom::make_sphere({0.0, 0.0, 0.0}, 1.0, -1.0), geom::Errc::tolerance_violation,
               "make_sphere: tolerance must be finite and non-negative");
  expect_error(geom::make_sphere({0.0, 0.0, 0.0}, 0.0, 0.0), geom::Errc::degenerate,
               "make_sphere: radius must be positive");
  expect_error(geom::make_sphere({0.0, 0.0, 0.0}, -1.0, 0.0), geom::Errc::degenerate,
               "make_sphere: radius must be positive");
  expect_error(geom::make_sphere({1e308, 1e308, 1e308}, 1e308, 0.0), geom::Errc::degenerate,
               "make_sphere: pole is not representable");
  expect_error(geom::make_sphere({0.0, 0.0, 1e308}, 1.0, 0.0), geom::Errc::degenerate,
               "make_sphere: seam has zero length");
}

TEST_CASE("make_torus rejects every invalid input class") {
  const geom::Vec3d z{0.0, 0.0, 1.0};
  for (int axis = 0; axis < 3; ++axis) {
    expect_error(geom::make_torus(point_with_nan(axis), z, 2.0, 0.5, 0.0), geom::Errc::non_finite,
                 "make_torus: input is not finite");
    expect_error(geom::make_torus({0.0, 0.0, 0.0}, vec_with_nan(axis), 2.0, 0.5, 0.0),
                 geom::Errc::non_finite, "make_torus: input is not finite");
  }
  expect_error(geom::make_torus({0.0, 0.0, 0.0}, z, kNaN, 0.5, 0.0), geom::Errc::non_finite,
               "make_torus: input is not finite");
  expect_error(geom::make_torus({0.0, 0.0, 0.0}, z, 2.0, kNaN, 0.0), geom::Errc::non_finite,
               "make_torus: input is not finite");
  expect_error(geom::make_torus({0.0, 0.0, 0.0}, z, 2.0, 0.5, -1.0), geom::Errc::tolerance_violation,
               "make_torus: tolerance must be finite and non-negative");
  expect_error(geom::make_torus({0.0, 0.0, 0.0}, z, 2.0, 0.5, kNaN), geom::Errc::tolerance_violation,
               "make_torus: tolerance must be finite and non-negative");
  expect_error(geom::make_torus({0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}, 2.0, 0.5, 0.0),
               geom::Errc::degenerate, "make_torus: axis must be unit length");
  expect_error(geom::make_torus({0.0, 0.0, 0.0}, {0.0, 0.0, 3.0}, 2.0, 0.5, 0.0),
               geom::Errc::degenerate, "make_torus: axis must be unit length");
  expect_error(geom::make_torus({0.0, 0.0, 0.0}, z, 0.0, 0.5, 0.0), geom::Errc::degenerate,
               "make_torus: major radius must be positive");
  expect_error(geom::make_torus({0.0, 0.0, 0.0}, z, -2.0, 0.5, 0.0), geom::Errc::degenerate,
               "make_torus: major radius must be positive");
  expect_error(geom::make_torus({0.0, 0.0, 0.0}, z, 2.0, 0.0, 0.0), geom::Errc::degenerate,
               "make_torus: minor radius must be positive");
  expect_error(geom::make_torus({0.0, 0.0, 0.0}, z, 2.0, -0.5, 0.0), geom::Errc::degenerate,
               "make_torus: minor radius must be positive");
  expect_error(geom::make_torus({0.0, 0.0, 0.0}, z, 2.0, 2.0, 0.0), geom::Errc::degenerate,
               "make_torus: minor radius must be smaller than the major radius");
  expect_error(geom::make_torus({0.0, 0.0, 0.0}, z, 2.0, 3.0, 0.0), geom::Errc::degenerate,
               "make_torus: minor radius must be smaller than the major radius");
  expect_error(geom::make_torus({0.0, 0.0, 0.0}, z, 1.5e308, 1e308, 0.0), geom::Errc::degenerate,
               "make_torus: seam circle radius is not representable");
  expect_error(geom::make_torus({1e308, 0.0, 0.0}, z, 6e307, 3e307, 0.0), geom::Errc::degenerate,
               "make_torus: seam point is not representable");
}

// =============================================================================
// solid_editor.hpp
// =============================================================================

TEST_CASE("SolidEditor rejects invalid vertices and edge specifications") {
  geom::SolidEditor editor;
  expect_error(editor.add_vertex({kNaN, 0.0, 0.0}, 0.0), geom::Errc::non_finite,
               "solid_editor: vertex position is not finite");
  expect_error(editor.add_vertex({0.0, kInf, 0.0}, 0.0), geom::Errc::non_finite,
               "solid_editor: vertex position is not finite");
  expect_error(editor.add_vertex({0.0, 0.0, 0.0}, kNaN), geom::Errc::tolerance_violation,
               "solid_editor: vertex tolerance must be finite and non-negative");
  expect_error(editor.add_vertex({0.0, 0.0, 0.0}, -1.0), geom::Errc::tolerance_violation,
               "solid_editor: vertex tolerance must be finite and non-negative");
  CHECK(editor.brep().vertex_count() == 0);

  const geom::VertexId v0 = editor.add_vertex({0.0, 0.0, 0.0}, 0.0).value();
  const geom::VertexId v1 = editor.add_vertex({1.0, 0.0, 0.0}, 0.0).value();
  const geom::VertexId dead{};
  expect_error(editor.add_edge(dead, v1, geom::kNoGeometry, 0.0, 1.0, 0.0),
               geom::Errc::index_out_of_range, "solid_editor: edge vertex is not alive");
  expect_error(editor.add_edge(v0, dead, geom::kNoGeometry, 0.0, 1.0, 0.0),
               geom::Errc::index_out_of_range, "solid_editor: edge vertex is not alive");
  expect_error(editor.add_edge(v0, v1, geom::kNoGeometry, 0.0, 1.0, kNaN),
               geom::Errc::tolerance_violation,
               "solid_editor: edge tolerance must be finite and non-negative");
  expect_error(editor.add_edge(v0, v1, geom::kNoGeometry, 0.0, 1.0, -1.0),
               geom::Errc::tolerance_violation,
               "solid_editor: edge tolerance must be finite and non-negative");
  expect_error(editor.add_edge(v0, v1, geom::kNoGeometry, 1.0, 1.0, 0.0), geom::Errc::degenerate,
               "solid_editor: edge parameter interval is empty or non-finite");
  expect_error(editor.add_edge(v0, v1, geom::kNoGeometry, 1.0, 0.0, 0.0), geom::Errc::degenerate,
               "solid_editor: edge parameter interval is empty or non-finite");
  expect_error(editor.add_edge(v0, v1, geom::kNoGeometry, kNaN, 1.0, 0.0), geom::Errc::degenerate,
               "solid_editor: edge parameter interval is empty or non-finite");
  expect_error(editor.add_edge(v0, v1, geom::kNoGeometry, 0.0, kInf, 0.0), geom::Errc::degenerate,
               "solid_editor: edge parameter interval is empty or non-finite");
  expect_error(editor.add_edge(v0, v1, 5u, 0.0, 1.0, 0.0), geom::Errc::index_out_of_range,
               "solid_editor: edge curve index out of range");
  CHECK(editor.brep().edge_count() == 0);
}

TEST_CASE("SolidEditor rejects invalid faces, loops, shells and solids") {
  geom::SolidEditor editor;
  const geom::VertexId v0 = editor.add_vertex({0.0, 0.0, 0.0}, 0.0).value();
  const geom::VertexId v1 = editor.add_vertex({1.0, 0.0, 0.0}, 0.0).value();
  const std::uint32_t plane = editor.add_surface(geom::Plane3{{0.0, 0.0, 1.0}, 0.0});
  const geom::EdgeId edge = add_line_edge(editor, {0.0, 0.0, 0.0}, {1.0, 0.0, 0.0});

  expect_error(editor.add_face(9u, false, 0.0), geom::Errc::index_out_of_range,
               "solid_editor: face surface index out of range");
  expect_error(editor.add_face(geom::kNoGeometry, false, kNaN), geom::Errc::tolerance_violation,
               "solid_editor: face tolerance must be finite and non-negative");
  const geom::FaceId face = editor.add_face(plane, false, 0.0).value();

  expect_error(editor.add_loop(geom::FaceId{}, std::vector<geom::CoedgeSpec>{}, true),
               geom::Errc::index_out_of_range, "solid_editor: loop face is not alive");
  expect_error(editor.add_loop(face, std::vector<geom::CoedgeSpec>{}, true),
               geom::Errc::invalid_topology, "solid_editor: loop needs at least one coedge");
  expect_error(editor.add_loop(face,
                               std::vector<geom::CoedgeSpec>{geom::CoedgeSpec{geom::EdgeId{}, false}},
                               true),
               geom::Errc::index_out_of_range, "solid_editor: loop coedge edge is not alive");
  expect_error(editor.add_loop(face,
                               std::vector<geom::CoedgeSpec>{geom::CoedgeSpec{
                                   edge, false, geom::kNoGeometry, kNaN}},
                               true),
               geom::Errc::tolerance_violation,
               "solid_editor: coedge tolerance must be finite and non-negative");
  expect_error(editor.add_loop(face,
                               std::vector<geom::CoedgeSpec>{geom::CoedgeSpec{edge, false, 4u, 0.0}},
                               true),
               geom::Errc::index_out_of_range,
               "solid_editor: loop coedge pcurve index out of range");
  expect_error(
      editor.add_loop(face,
                      std::vector<geom::CoedgeSpec>{geom::CoedgeSpec{edge, false},
                                                    geom::CoedgeSpec{edge, false},
                                                    geom::CoedgeSpec{edge, true}},
                      true),
      geom::Errc::invalid_topology, "solid_editor: edge used more than twice in a loop");
  CHECK(editor.brep().loop_count() == 0);
  // Two uses of one edge are legal (the seam pattern).
  const geom::LoopId loop = editor.add_loop(
      face, std::array<geom::CoedgeSpec, 2>{geom::CoedgeSpec{edge, false}, geom::CoedgeSpec{edge, true}},
      true)
                                .value();
  CHECK(editor.brep().loop_count() == 1);

  expect_error(editor.add_shell(std::span<const geom::FaceId>(), false), geom::Errc::invalid_topology,
               "solid_editor: shell needs at least one face");
  expect_error(editor.add_shell(std::array<geom::FaceId, 1>{geom::FaceId{}}, false),
               geom::Errc::index_out_of_range, "solid_editor: shell face is not alive");
  const geom::ShellId shell = editor.add_shell(std::array<geom::FaceId, 1>{face}, false).value();

  expect_error(editor.add_solid(std::span<const geom::ShellId>(), 0.0), geom::Errc::invalid_topology,
               "solid_editor: solid needs at least one shell");
  expect_error(editor.add_solid(std::array<geom::ShellId, 1>{geom::ShellId{}}, 0.0),
               geom::Errc::index_out_of_range, "solid_editor: solid shell is not alive");
  expect_error(editor.add_solid(std::array<geom::ShellId, 1>{shell}, kNaN),
               geom::Errc::tolerance_violation,
               "solid_editor: solid tolerance must be finite and non-negative");
  CHECK(editor.brep().solid_count() == 0);
  CHECK(loop.is_null() == false);
  CHECK(v0 != v1);
  CHECK(editor.add_solid(std::array<geom::ShellId, 1>{shell}, 0.0).has_value());
}

TEST_CASE("SolidEditor erase and split_edge validate their arguments") {
  Square square = make_square();
  geom::SolidEditor& editor = square.editor;
  const std::size_t edge_count = editor.brep().edge_count();
  const std::size_t vertex_count = editor.brep().vertex_count();

  expect_error(editor.erase_edge(geom::EdgeId{}), geom::Errc::index_out_of_range,
               "solid_editor: edge is not alive");
  expect_error(editor.erase_vertex(geom::VertexId{}), geom::Errc::index_out_of_range,
               "solid_editor: vertex is not alive");

  expect_error(editor.split_edge(geom::EdgeId{}, 0.5), geom::Errc::index_out_of_range,
               "solid_editor: split edge is not alive");
  expect_error(editor.split_edge(square.edges[0], 0.0), geom::Errc::degenerate,
               "solid_editor: split parameter is not inside the edge interval");
  expect_error(editor.split_edge(square.edges[0], 1.0), geom::Errc::degenerate,
               "solid_editor: split parameter is not inside the edge interval");
  expect_error(editor.split_edge(square.edges[0], kNaN), geom::Errc::degenerate,
               "solid_editor: split parameter is not inside the edge interval");
  const geom::EdgeId no_curve = editor.add_edge(square.vertices[0], square.vertices[1],
                                                geom::kNoGeometry, 0.0, 1.0, 0.0)
                                   .value();
  expect_error(editor.split_edge(no_curve, 0.5), geom::Errc::degenerate,
               "solid_editor: split edge has no curve");
  const std::uint32_t huge = editor.add_curve(
      geom::LineCurve3{{1e308, 0.0, 0.0}, {1e308, 0.0, 0.0}});
  const geom::EdgeId exploding = editor.add_edge(square.vertices[0], square.vertices[1], huge, 0.0,
                                                 2.0, 0.0)
                                     .value();
  expect_error(editor.split_edge(exploding, 1.5), geom::Errc::non_finite,
               "solid_editor: split point is not finite");
  CHECK(editor.brep().edge_count() == edge_count + 2);
  CHECK(editor.brep().vertex_count() == vertex_count);

  // A successful split adds one vertex and one edge, and rewrites the loop.
  const geom::VertexId middle = editor.split_edge(square.edges[0], 0.5).value();
  CHECK(editor.brep().vertex_count() == vertex_count + 1);
  CHECK(editor.brep().edge_count() == edge_count + 3);
  CHECK(middle.index() < editor.brep().vertices().size());
  CHECK(editor.brep().edge_ids().size() == editor.brep().edge_count());

  // Erasing a live entity succeeds and is reported as dead afterwards.
  CHECK(editor.erase_edge(no_curve).has_value());
  expect_error(editor.erase_edge(no_curve), geom::Errc::index_out_of_range,
               "solid_editor: edge is not alive");
  CHECK(editor.erase_vertex(square.vertices[3]).has_value());
  expect_error(editor.erase_vertex(square.vertices[3]), geom::Errc::index_out_of_range,
               "solid_editor: vertex is not alive");
}

TEST_CASE("SolidEditor split_edge splits pcurves along the edge") {
  geom::SolidEditor editor;
  const geom::Point3d pts[4] = {{0.0, 0.0, 0.0}, {1.0, 0.0, 0.0}, {1.0, 1.0, 0.0}, {0.0, 1.0, 0.0}};
  std::array<geom::VertexId, 4> vertices{};
  for (int i = 0; i < 4; ++i) vertices[static_cast<std::size_t>(i)] = editor.add_vertex(pts[i], 0.0).value();
  const std::uint32_t plane = editor.add_surface(geom::Plane3{{0.0, 0.0, 1.0}, 0.0});
  std::array<geom::EdgeId, 4> edges{};
  for (int i = 0; i < 4; ++i) {
    const geom::Point3d from = pts[i];
    const std::uint32_t curve = editor.add_curve(geom::LineCurve3{from, {1.0, 0.0, 0.0}});
    edges[static_cast<std::size_t>(i)] =
        editor.add_edge(vertices[static_cast<std::size_t>(i)],
                        vertices[static_cast<std::size_t>((i + 1) % 4)], curve, 0.0, 1.0, 0.0)
            .value();
  }
  const std::uint32_t pcurve = editor.add_pcurve(geom::Pcurve2d{{0.0, 0.0}, {1.0, 1.0}});
  const geom::FaceId face = editor.add_face(plane, false, 0.0).value();
  std::array<geom::CoedgeSpec, 4> specs{};
  for (int i = 0; i < 4; ++i) {
    specs[static_cast<std::size_t>(i)] =
        geom::CoedgeSpec{edges[static_cast<std::size_t>(i)], false, pcurve, 0.0};
  }
  editor.add_loop(face, specs, true).value();

  CHECK(editor.split_edge(edges[0], 0.5).has_value());
  CHECK(editor.brep().geometry().pcurves.size() >= 2u);
}

TEST_CASE("SolidEditor join_vertex reports every rejected configuration") {
  // Dead vertex.
  {
    geom::SolidEditor editor;
    expect_error(editor.join_vertex(geom::VertexId{}), geom::Errc::index_out_of_range,
                 "solid_editor: join vertex is not alive");
  }
  // Zero and three incident edges.
  {
    geom::SolidEditor editor;
    const geom::VertexId lone = editor.add_vertex({0.0, 0.0, 0.0}, 0.0).value();
    expect_error(editor.join_vertex(lone), geom::Errc::invalid_topology,
                 "solid_editor: join vertex does not have exactly two incident edges");
    const geom::VertexId a = editor.add_vertex({1.0, 0.0, 0.0}, 0.0).value();
    const geom::VertexId b = editor.add_vertex({-1.0, 0.0, 0.0}, 0.0).value();
    const geom::VertexId c = editor.add_vertex({0.0, 1.0, 0.0}, 0.0).value();
    add_line_edge(editor, {0.0, 0.0, 0.0}, {1.0, 0.0, 0.0});
    add_line_edge(editor, {0.0, 0.0, 0.0}, {-1.0, 0.0, 0.0});
    add_line_edge(editor, {0.0, 0.0, 0.0}, {0.0, 1.0, 0.0});
    CHECK(editor.brep().edge_count() == 3);
    CHECK(a != b);
    CHECK(b != c);
    expect_error(editor.join_vertex(lone), geom::Errc::invalid_topology,
                 "solid_editor: join vertex does not have exactly two incident edges");
  }
  // Two incident edges, but only one is used by a coedge.
  {
    geom::SolidEditor editor;
    const geom::VertexId a = editor.add_vertex({0.0, 0.0, 0.0}, 0.0).value();
    const geom::VertexId v = editor.add_vertex({1.0, 0.0, 0.0}, 0.0).value();
    const geom::VertexId b = editor.add_vertex({2.0, 0.0, 0.0}, 0.0).value();
    const std::uint32_t plane = editor.add_surface(geom::Plane3{{0.0, 0.0, 1.0}, 0.0});
    const std::uint32_t curve = editor.add_curve(geom::LineCurve3{{0.0, 0.0, 0.0}, {1.0, 0.0, 0.0}});
    const geom::EdgeId left = editor.add_edge(a, v, curve, 0.0, 1.0, 0.0).value();
    editor.add_edge(v, b, curve, 1.0, 2.0, 0.0).value();
    const geom::FaceId face = editor.add_face(plane, false, 0.0).value();
    add_loop_of(editor, face, {geom::CoedgeSpec{left, false}});
    expect_error(editor.join_vertex(v), geom::Errc::invalid_topology,
                 "solid_editor: join vertex does not have two incident coedges");
  }
  // Both coedge uses are the same edge.
  {
    geom::SolidEditor editor;
    const geom::VertexId a = editor.add_vertex({0.0, 0.0, 0.0}, 0.0).value();
    const geom::VertexId v = editor.add_vertex({1.0, 0.0, 0.0}, 0.0).value();
    const geom::VertexId b = editor.add_vertex({2.0, 0.0, 0.0}, 0.0).value();
    const std::uint32_t plane = editor.add_surface(geom::Plane3{{0.0, 0.0, 1.0}, 0.0});
    const std::uint32_t curve = editor.add_curve(geom::LineCurve3{{0.0, 0.0, 0.0}, {1.0, 0.0, 0.0}});
    const geom::EdgeId left = editor.add_edge(a, v, curve, 0.0, 1.0, 0.0).value();
    editor.add_edge(v, b, curve, 1.0, 2.0, 0.0).value();
    const geom::FaceId face = editor.add_face(plane, false, 0.0).value();
    add_loop_of(editor, face, {geom::CoedgeSpec{left, false}, geom::CoedgeSpec{left, true}});
    expect_error(editor.join_vertex(v), geom::Errc::invalid_topology,
                 "solid_editor: join vertex does not have two incident coedges");
  }
  // Both coedges run into the vertex.
  {
    geom::SolidEditor editor;
    const geom::Point3d pa{0.0, 0.0, 0.0};
    const geom::Point3d pv{1.0, 0.0, 0.0};
    const geom::Point3d pb{2.0, 0.0, 0.0};
    const geom::VertexId a = editor.add_vertex(pa, 0.0).value();
    const geom::VertexId v = editor.add_vertex(pv, 0.0).value();
    const geom::VertexId b = editor.add_vertex(pb, 0.0).value();
    const std::uint32_t plane = editor.add_surface(geom::Plane3{{0.0, 0.0, 1.0}, 0.0});
    const std::uint32_t curve = editor.add_curve(geom::LineCurve3{pa, {1.0, 0.0, 0.0}});
    const geom::EdgeId left = editor.add_edge(a, v, curve, 0.0, 1.0, 0.0).value();
    const geom::EdgeId right = editor.add_edge(v, b, curve, 1.0, 2.0, 0.0).value();
    const geom::FaceId face = editor.add_face(plane, false, 0.0).value();
    // right: v0 == v, so traversing it reversed arrives at v as well.
    add_loop_of(editor, face, {geom::CoedgeSpec{left, false}, geom::CoedgeSpec{right, true}});
    expect_error(editor.join_vertex(v), geom::Errc::invalid_topology,
                 "solid_editor: join coedges do not run in opposite directions at the vertex");
  }
  // Two uses in different loops.
  {
    geom::SolidEditor editor;
    const geom::Point3d pa{0.0, 0.0, 0.0};
    const geom::Point3d pv{1.0, 0.0, 0.0};
    const geom::Point3d pb{2.0, 0.0, 0.0};
    const geom::VertexId a = editor.add_vertex(pa, 0.0).value();
    const geom::VertexId v = editor.add_vertex(pv, 0.0).value();
    const geom::VertexId b = editor.add_vertex(pb, 0.0).value();
    const std::uint32_t plane = editor.add_surface(geom::Plane3{{0.0, 0.0, 1.0}, 0.0});
    const std::uint32_t curve = editor.add_curve(geom::LineCurve3{pa, {1.0, 0.0, 0.0}});
    const geom::EdgeId left = editor.add_edge(a, v, curve, 0.0, 1.0, 0.0).value();
    const geom::EdgeId right = editor.add_edge(v, b, curve, 1.0, 2.0, 0.0).value();
    const geom::FaceId face_a = editor.add_face(plane, false, 0.0).value();
    const geom::FaceId face_b = editor.add_face(plane, false, 0.0).value();
    add_loop_of(editor, face_a, {geom::CoedgeSpec{left, false}});
    add_loop_of(editor, face_b, {geom::CoedgeSpec{right, false}});
    expect_error(editor.join_vertex(v), geom::Errc::invalid_topology,
                 "solid_editor: join coedges are not adjacent in one loop");
  }
  // Edges do not share one curve.
  {
    geom::SolidEditor editor;
    const geom::Point3d pa{0.0, 0.0, 0.0};
    const geom::Point3d pv{1.0, 0.0, 0.0};
    const geom::Point3d pb{2.0, 0.0, 0.0};
    const geom::VertexId a = editor.add_vertex(pa, 0.0).value();
    const geom::VertexId v = editor.add_vertex(pv, 0.0).value();
    const geom::VertexId b = editor.add_vertex(pb, 0.0).value();
    const std::uint32_t plane = editor.add_surface(geom::Plane3{{0.0, 0.0, 1.0}, 0.0});
    const geom::VertexId unused = a;
    const geom::EdgeId left =
        editor.add_edge(a, v, geom::kNoGeometry, 0.0, 1.0, 0.0).value();
    const geom::EdgeId right = editor.add_edge(v, b, geom::kNoGeometry, 1.0, 2.0, 0.0).value();
    const geom::FaceId face = editor.add_face(plane, false, 0.0).value();
    add_loop_of(editor, face, {geom::CoedgeSpec{left, false}, geom::CoedgeSpec{right, false}});
    CHECK(unused == a);
    expect_error(editor.join_vertex(v), geom::Errc::degenerate,
                 "solid_editor: join edges do not share one curve");
  }
  // Shared curve but intervals that do not meet.
  {
    geom::SolidEditor editor;
    const geom::Point3d pa{0.0, 0.0, 0.0};
    const geom::Point3d pv{1.0, 0.0, 0.0};
    const geom::Point3d pb{2.0, 0.0, 0.0};
    const geom::VertexId a = editor.add_vertex(pa, 0.0).value();
    const geom::VertexId v = editor.add_vertex(pv, 0.0).value();
    const geom::VertexId b = editor.add_vertex(pb, 0.0).value();
    const std::uint32_t plane = editor.add_surface(geom::Plane3{{0.0, 0.0, 1.0}, 0.0});
    const std::uint32_t curve = editor.add_curve(geom::LineCurve3{pa, {1.0, 0.0, 0.0}});
    const geom::EdgeId left = editor.add_edge(a, v, curve, 0.0, 0.4, 0.0).value();
    const geom::EdgeId right = editor.add_edge(v, b, curve, 0.5, 1.0, 0.0).value();
    const geom::FaceId face = editor.add_face(plane, false, 0.0).value();
    add_loop_of(editor, face, {geom::CoedgeSpec{left, false}, geom::CoedgeSpec{right, false}});
    expect_error(editor.join_vertex(v), geom::Errc::degenerate,
                 "solid_editor: join edge intervals do not meet at the vertex");
  }
  // Mixed orientations.
  {
    geom::SolidEditor editor;
    const geom::Point3d pv{0.0, 0.0, 0.0};
    const geom::Point3d pa{1.0, 0.0, 0.0};
    const geom::Point3d pb{2.0, 0.0, 0.0};
    const geom::VertexId v = editor.add_vertex(pv, 0.0).value();
    const geom::VertexId a = editor.add_vertex(pa, 0.0).value();
    const geom::VertexId b = editor.add_vertex(pb, 0.0).value();
    const std::uint32_t plane = editor.add_surface(geom::Plane3{{0.0, 0.0, 1.0}, 0.0});
    const std::uint32_t curve = editor.add_curve(geom::LineCurve3{pv, {1.0, 0.0, 0.0}});
    const geom::EdgeId left = editor.add_edge(v, a, curve, 0.0, 1.0, 0.0).value();
    const geom::EdgeId right = editor.add_edge(v, b, curve, 1.0, 2.0, 0.0).value();
    const geom::FaceId face = editor.add_face(plane, false, 0.0).value();
    add_loop_of(editor, face, {geom::CoedgeSpec{left, true}, geom::CoedgeSpec{right, false}});
    expect_error(editor.join_vertex(v), geom::Errc::degenerate,
                 "solid_editor: join coedges have inconsistent orientation");
  }
  // Zero-length tangent: a line curve with a zero direction has no tangent.
  {
    geom::SolidEditor editor;
    const geom::Point3d pa{0.0, 0.0, 0.0};
    const geom::Point3d pv{1.0, 0.0, 0.0};
    const geom::Point3d pb{2.0, 0.0, 0.0};
    const geom::VertexId a = editor.add_vertex(pa, 0.0).value();
    const geom::VertexId v = editor.add_vertex(pv, 0.0).value();
    const geom::VertexId b = editor.add_vertex(pb, 0.0).value();
    const std::uint32_t plane = editor.add_surface(geom::Plane3{{0.0, 0.0, 1.0}, 0.0});
    const std::uint32_t curve = editor.add_curve(geom::LineCurve3{pa, {0.0, 0.0, 0.0}});
    const geom::EdgeId left = editor.add_edge(a, v, curve, 0.0, 0.5, 0.0).value();
    const geom::EdgeId right = editor.add_edge(v, b, curve, 0.5, 1.0, 0.0).value();
    const geom::FaceId face = editor.add_face(plane, false, 0.0).value();
    add_loop_of(editor, face, {geom::CoedgeSpec{left, false}, geom::CoedgeSpec{right, false}});
    expect_error(editor.join_vertex(v), geom::Errc::degenerate,
                 "solid_editor: join edges are not tangent at the vertex");
  }
}

TEST_CASE("SolidEditor split then join restores the original topology") {
  Square square = make_square();
  const std::size_t vertices = square.editor.brep().vertex_count();
  const std::size_t edges = square.editor.brep().edge_count();
  const geom::VertexId middle = square.editor.split_edge(square.edges[0], 0.5).value();
  CHECK(square.editor.brep().vertex_count() == vertices + 1);
  CHECK(square.editor.brep().edge_count() == edges + 1);
  CHECK(square.editor.join_vertex(middle).has_value());
  CHECK(square.editor.brep().vertex_count() == vertices);
  CHECK(square.editor.brep().edge_count() == edges);
  CHECK(geom::validate_topology(square.editor.brep()).ok());
}

// =============================================================================
// validate.hpp
// =============================================================================

TEST_CASE("validate_topology accepts a well-formed cube and counts Euler") {
  Cube cube = make_cube();
  const geom::Brep& brep = cube.editor.brep();
  CHECK(geom::validate_topology(brep).ok());
  CHECK(geom::validate_tolerances(brep).ok());
  CHECK(geom::euler_characteristic(brep, cube.solid) == 2);
  CHECK(cube.editor.commit().has_value());
}

TEST_CASE("validate_topology reports dangling edges and vertices") {
  {
    Cube cube = make_cube();
    geom::SolidEditor editor(cube.editor.brep());
    REQUIRE(editor.erase_edge(cube.edges[0]).has_value());
    const geom::ValidationReport report = geom::validate_topology(editor.brep());
    REQUIRE_FALSE(report.ok());
    CHECK(report.issues.front().code == geom::Errc::invalid_topology);
    CHECK(first_message(report) == "coedge edge is not alive");
    CHECK(report.issues.size() == 2u);
    CHECK(editor.commit().error().code == geom::Errc::invalid_topology);
    CHECK(std::string_view(editor.commit().error().message) == "coedge edge is not alive");
  }
  {
    Cube cube = make_cube();
    geom::SolidEditor editor(cube.editor.brep());
    REQUIRE(editor.erase_vertex(cube.vertices[0]).has_value());
    const geom::ValidationReport report = geom::validate_topology(editor.brep());
    CHECK(has_issue(report, geom::Errc::invalid_topology, "edge start vertex is not alive"));
  }
  {
    Cube cube = make_cube();
    geom::SolidEditor editor(cube.editor.brep());
    REQUIRE(editor.erase_vertex(cube.vertices[7]).has_value());
    const geom::ValidationReport report = geom::validate_topology(editor.brep());
    CHECK(has_issue(report, geom::Errc::invalid_topology, "edge end vertex is not alive"));
  }
}

TEST_CASE("validate_topology reports face loop rule violations") {
  {
    geom::SolidEditor editor;
    const std::uint32_t plane = editor.add_surface(geom::Plane3{{0.0, 0.0, 1.0}, 0.0});
    const geom::FaceId face = editor.add_face(plane, false, 0.0).value();
    editor.add_shell(std::array<geom::FaceId, 1>{face}, false).value();
    const geom::ValidationReport report = geom::validate_topology(editor.brep());
    CHECK(first_message(report) == "face has no loops");
    CHECK(report.issues.front().code == geom::Errc::invalid_topology);
    CHECK(report.issues.front().entity == face.index());
  }
  {
    SquareOptions options;
    options.outer = false;
    Square square = make_square(options);
    const geom::ValidationReport report = geom::validate_topology(square.editor.brep());
    CHECK(first_message(report) == "face has no outer loop");
  }
  {
    Square square = make_square();
    std::vector<geom::CoedgeSpec> second;
    for (const geom::EdgeId edge : square.edges) second.push_back(geom::CoedgeSpec{edge, true});
    add_loop_of(square.editor, square.face, second, true);
    const geom::ValidationReport report = geom::validate_topology(square.editor.brep());
    CHECK(first_message(report) == "face has more than one outer loop");
  }
}

TEST_CASE("validate_topology reports edge use rule violations") {
  {
    Square square = make_square();
    add_line_edge(square.editor, {0.0, 0.0, 0.0}, {1.0, 1.0, 0.0});
    const geom::ValidationReport report = geom::validate_topology(square.editor.brep());
    CHECK(first_message(report) == "edge is not used by any coedge");
  }
  {
    geom::SolidEditor editor;
    const std::uint32_t plane = editor.add_surface(geom::Plane3{{0.0, 0.0, 1.0}, 0.0});
    const geom::EdgeId edge = add_line_edge(editor, {0.0, 0.0, 0.0}, {1.0, 0.0, 0.0});
    for (int i = 0; i < 3; ++i) {
      const geom::FaceId face = editor.add_face(plane, false, 0.0).value();
      add_loop_of(editor, face, {geom::CoedgeSpec{edge, false}});
    }
    const geom::ValidationReport report = geom::validate_topology(editor.brep());
    CHECK(first_message(report) == "edge is used by more than two coedges");
    CHECK(report.issues.front().entity == edge.index());
  }
  {
    geom::SolidEditor editor;
    const std::uint32_t plane = editor.add_surface(geom::Plane3{{0.0, 0.0, 1.0}, 0.0});
    const geom::EdgeId edge = add_line_edge(editor, {0.0, 0.0, 0.0}, {1.0, 0.0, 0.0});
    for (int i = 0; i < 2; ++i) {
      const geom::FaceId face = editor.add_face(plane, false, 0.0).value();
      add_loop_of(editor, face, {geom::CoedgeSpec{edge, false}});
    }
    const geom::ValidationReport report = geom::validate_topology(editor.brep());
    CHECK(first_message(report) == "both coedges on an edge run in the same direction");
  }
}

TEST_CASE("validate_topology reports shell and solid container violations") {
  {
    Square square = make_square();
    square.editor.add_shell(std::array<geom::FaceId, 2>{square.face, square.face}, false).value();
    const geom::ValidationReport report = geom::validate_topology(square.editor.brep());
    CHECK(first_message(report) == "shell lists a face more than once");
  }
  {
    SquareOptions options;
    options.closed_shell = true;
    Square square = make_square(options);
    const geom::ValidationReport report = geom::validate_topology(square.editor.brep());
    REQUIRE_FALSE(report.ok());
    CHECK(report.issues.front().code == geom::Errc::not_closed);
    CHECK(first_message(report) == "closed shell has an edge without two coedges");
    CHECK(report.issues.size() == 4u);
  }
  {
    Square square = make_square();
    const geom::ShellId shell =
        square.editor.add_shell(std::array<geom::FaceId, 1>{square.face}, false).value();
    square.editor.add_solid(std::array<geom::ShellId, 2>{shell, shell}, 0.0).value();
    const geom::ValidationReport report = geom::validate_topology(square.editor.brep());
    CHECK(first_message(report) == "solid lists a shell more than once");
  }
  {
    Square square = make_square();
    const geom::ShellId first =
        square.editor.add_shell(std::array<geom::FaceId, 1>{square.face}, false).value();
    const geom::ShellId second =
        square.editor.add_shell(std::array<geom::FaceId, 1>{square.face}, false).value();
    square.editor.add_solid(std::array<geom::ShellId, 2>{first, second}, 0.0).value();
    const geom::ValidationReport report = geom::validate_topology(square.editor.brep());
    CHECK(first_message(report) == "multi-shell solids are deferred to M3");
    CHECK(report.issues.front().code == geom::Errc::unsupported);
    CHECK(square.editor.commit().error().code == geom::Errc::unsupported);
  }
}

TEST_CASE("validate_topology stops after the issue cap") {
  Square square = make_square();
  for (int i = 0; i < 10; ++i) {
    const double offset = static_cast<double>(i);
    add_line_edge(square.editor, {0.0, 0.0, 0.0}, {offset + 2.0, 0.0, 0.0});
  }
  const geom::ValidationReport report = geom::validate_topology(square.editor.brep());
  CHECK(report.issues.size() == geom::detail::kMaxValidationIssues);
  for (const geom::ValidationIssue& issue : report.issues) {
    CHECK(std::string_view(issue.message) == "edge is not used by any coedge");
  }
  CHECK(geom::euler_characteristic(square.editor.brep(), square.solid) == 1);
}

TEST_CASE("validate_tolerances enforces the parent/child hierarchy") {
  // Rule 2: the finer vertex is reported, and the larger of two offenders wins.
  {
    SquareOptions options;
    options.vertex_tol = {2.0, 1.0, 0.0, 0.0};
    Square square = make_square(options);
    const geom::ValidationReport report = geom::validate_tolerances(square.editor.brep());
    REQUIRE_FALSE(report.ok());
    CHECK(first_message(report) == "vertex tolerance exceeds the edge tolerance");
    CHECK(report.issues.front().entity == square.vertices[0].index());
    CHECK(report.issues.size() == 3u);
  }
  {
    SquareOptions options;
    options.vertex_tol = {1.0, 2.0, 0.0, 0.0};
    Square square = make_square(options);
    const geom::ValidationReport report = geom::validate_tolerances(square.editor.brep());
    REQUIRE_FALSE(report.ok());
    CHECK(first_message(report) == "vertex tolerance exceeds the edge tolerance");
    CHECK(report.issues.front().entity == square.vertices[1].index());
  }
  // Rule 3a: an edge tolerance may not exceed its face tolerance.
  {
    SquareOptions options;
    options.edge_tol = 0.5;
    Square square = make_square(options);
    const geom::ValidationReport report = geom::validate_tolerances(square.editor.brep());
    REQUIRE_FALSE(report.ok());
    CHECK(first_message(report) == "edge tolerance exceeds the face tolerance");
    CHECK(report.issues.front().entity == square.edges[0].index());
  }
  // Rule 3b: a coedge tolerance may not exceed its face tolerance.
  {
    SquareOptions options;
    options.coedge_tol = 1.0;
    Square square = make_square(options);
    const geom::ValidationReport report = geom::validate_tolerances(square.editor.brep());
    REQUIRE_FALSE(report.ok());
    CHECK(first_message(report) == "coedge tolerance exceeds the face tolerance");
  }
  // Rule 4: a face tolerance may not exceed the solid tolerance.
  {
    SquareOptions options;
    options.face_tol = 1.0;
    Square square = make_square(options);
    const geom::ValidationReport report = geom::validate_tolerances(square.editor.brep());
    REQUIRE_FALSE(report.ok());
    CHECK(first_message(report) == "face tolerance exceeds the solid tolerance");
    CHECK(report.issues.front().entity == square.face.index());
  }
  // Rule 5a: edge tolerance growing more than 10x over a vertex tolerance.
  {
    SquareOptions options;
    options.vertex_tol = {1e-6, 1e-6, 1e-6, 1e-6};
    options.edge_tol = 1e-3;
    options.face_tol = 1e-3;
    options.solid_tol = 1e-3;
    Square square = make_square(options);
    const geom::ValidationReport report = geom::validate_tolerances(square.editor.brep());
    REQUIRE_FALSE(report.ok());
    CHECK(first_message(report) == "tolerance growth exceeds 10x");
    CHECK(report.issues.front().entity == square.vertices[0].index());
  }
  // Rule 5b: face tolerance growing more than 10x over an edge tolerance.
  {
    SquareOptions options;
    options.edge_tol = 1e-6;
    options.face_tol = 1e-3;
    options.solid_tol = 1e-3;
    Square square = make_square(options);
    const geom::ValidationReport report = geom::validate_tolerances(square.editor.brep());
    REQUIRE_FALSE(report.ok());
    CHECK(first_message(report) == "tolerance growth exceeds 10x");
  }
  // Rule 5c: solid tolerance growing more than 10x over a face tolerance.
  {
    SquareOptions options;
    options.face_tol = 1e-6;
    options.solid_tol = 1e-3;
    Square square = make_square(options);
    const geom::ValidationReport report = geom::validate_tolerances(square.editor.brep());
    REQUIRE_FALSE(report.ok());
    CHECK(first_message(report) == "tolerance growth exceeds 10x");
    CHECK(report.issues.front().entity == square.face.index());
  }
  // A cube with mixed face tolerances reports the worst offending face.
  {
    CubeOptions options;
    options.face_tol = {0.3, 0.1, 0.2, 0.05, 0.4, 0.05};
    options.solid_tol = 0.02;
    Cube cube = make_cube(options);
    const geom::ValidationReport report = geom::validate_tolerances(cube.editor.brep());
    REQUIRE_FALSE(report.ok());
    CHECK(first_message(report) == "face tolerance exceeds the solid tolerance");
    CHECK(report.issues.front().entity == cube.faces[4].index());
  }
}

// =============================================================================
// mass_properties.hpp
// =============================================================================

namespace {

struct PlanarFace {
  geom::SolidEditor editor;
  std::vector<geom::VertexId> vertices;
  std::vector<geom::EdgeId> edges;
  geom::FaceId face{};
  geom::LoopId loop{};
  std::uint32_t surface = geom::kNoGeometry;
};

struct PlanarOptions {
  geom::Plane3 plane{{0.0, 0.0, 1.0}, 0.0};
  bool outer = true;
  bool with_loop = true;
  bool use_curves = true;
  bool closed_shell = false;
  bool with_shell = false;
};

// Arbitrary planar polygon; edges either carry line geometry matching the
// vertices or no geometry at all (the tessellator then uses the vertex points).
PlanarFace make_planar_face(const std::vector<geom::Point3d>& pts, const PlanarOptions& options = {}) {
  PlanarFace result;
  result.surface = result.editor.add_surface(options.plane);
  for (const geom::Point3d& point : pts) {
    result.vertices.push_back(result.editor.add_vertex(point, 0.0).value());
  }
  std::vector<geom::CoedgeSpec> specs;
  specs.reserve(pts.size());
  for (std::size_t i = 0; i < pts.size(); ++i) {
    const geom::Point3d from = pts[i];
    const geom::Point3d to = pts[(i + 1) % pts.size()];
    std::uint32_t curve = geom::kNoGeometry;
    double t1 = 1.0;
    if (options.use_curves) {
      const geom::Vec3d delta = to - from;
      const double length = std::hypot(delta[0], delta[1], delta[2]);
      curve = result.editor.add_curve(
          geom::LineCurve3{from, {delta[0] / length, delta[1] / length, delta[2] / length}});
      t1 = length;
    }
    result.edges.push_back(
        result.editor
            .add_edge(result.vertices[i], result.vertices[(i + 1) % pts.size()], curve, 0.0, t1, 0.0)
            .value());
    specs.push_back(geom::CoedgeSpec{result.edges.back(), false});
  }
  result.face = result.editor.add_face(result.surface, false, 0.0).value();
  if (options.with_loop) {
    result.loop = result.editor.add_loop(result.face, specs, options.outer).value();
  }
  if (options.with_shell) {
    result.editor
        .add_shell(std::array<geom::FaceId, 1>{result.face}, options.closed_shell)
        .value();
  }
  return result;
}

void wrap_in_closed_solid(PlanarFace& face) {
  const geom::ShellId shell =
      face.editor.add_shell(std::array<geom::FaceId, 1>{face.face}, true).value();
  face.editor.add_solid(std::array<geom::ShellId, 1>{shell}, 0.0).value();
}

// Planar disk bounded by one closed circle edge.
PlanarFace make_circle_face(const geom::CircleCurve3& circle, double t0, double t1) {
  PlanarFace result;
  result.surface = result.editor.add_surface(geom::Plane3{{0.0, 0.0, 1.0}, 0.0});
  const std::uint32_t curve = result.editor.add_curve(circle);
  const geom::VertexId start = result.editor.add_vertex(circle.eval(t0), 0.0).value();
  const geom::VertexId end =
      result.editor.add_vertex(circle.eval(t1), 0.0).value();
  const geom::EdgeId edge =
      result.editor.add_edge(start, end, curve, t0, t1, 0.0).value();
  result.face = result.editor.add_face(result.surface, false, 0.0).value();
  result.loop =
      result.editor.add_loop(result.face, std::array<geom::CoedgeSpec, 1>{geom::CoedgeSpec{edge, false}},
                             true)
          .value();
  result.edges.push_back(edge);
  result.vertices.push_back(start);
  return result;
}

}  // namespace

TEST_CASE("mass_properties validates options and empty state") {
  Cube cube = make_cube();
  const geom::Brep& brep = cube.editor.brep();
  expect_error(geom::mass_properties(brep, geom::MassPropertiesOptions{0.0, 24}),
               geom::Errc::tolerance_violation, "mass_properties: options are out of range");
  expect_error(geom::mass_properties(brep, geom::MassPropertiesOptions{-1.0, 24}),
               geom::Errc::tolerance_violation, "mass_properties: options are out of range");
  expect_error(geom::mass_properties(brep, geom::MassPropertiesOptions{kNaN, 24}),
               geom::Errc::tolerance_violation, "mass_properties: options are out of range");
  expect_error(geom::mass_properties(brep, geom::MassPropertiesOptions{1e-9, 0}),
               geom::Errc::tolerance_violation, "mass_properties: options are out of range");
  expect_error(geom::mass_properties(brep, geom::MassPropertiesOptions{1e-9, -3}),
               geom::Errc::tolerance_violation, "mass_properties: options are out of range");

  const geom::Brep empty;
  expect_error(geom::mass_properties(empty), geom::Errc::numeric_failure,
               "mass_properties: brep has no shells");

  // A shell without any vertex cannot establish a reference point.
  {
    geom::SolidEditor editor;
    const std::uint32_t plane = editor.add_surface(geom::Plane3{{0.0, 0.0, 1.0}, 0.0});
    const geom::FaceId face = editor.add_face(plane, false, 0.0).value();
    const geom::ShellId shell = editor.add_shell(std::array<geom::FaceId, 1>{face}, true).value();
    editor.add_solid(std::array<geom::ShellId, 1>{shell}, 0.0).value();
    expect_error(geom::mass_properties(editor.brep()), geom::Errc::numeric_failure,
                 "mass_properties: brep has no vertices");
  }
  // A finite vertex set whose bounding-box centre overflows.
  {
    const std::vector<geom::Point3d> pts{{1e308, 0.0, 0.0},
                                         {1.7e308, 0.0, 0.0},
                                         {1.7e308, 1.0, 0.0},
                                         {1e308, 1.0, 0.0}};
    PlanarOptions options;
    options.use_curves = false;
    PlanarFace plan = make_planar_face(pts, options);
    wrap_in_closed_solid(plan);
    expect_error(geom::mass_properties(plan.editor.brep()), geom::Errc::numeric_failure,
                 "mass_properties: vertex bounds are not finite");
  }
}

TEST_CASE("mass_properties rejects open shells and malformed planar faces") {
  {
    Square square = make_square();
    const auto result = geom::mass_properties(square.editor.brep());
    expect_error(result, geom::Errc::not_closed, "mass_properties: shell is not closed");
    CHECK(result.error().context_a == square.shell.index());
  }
  {
    PlanarOptions options;
    options.with_loop = false;
    PlanarFace plan = make_planar_face({{0.0, 0.0, 0.0}, {1.0, 0.0, 0.0}, {1.0, 1.0, 0.0}}, options);
    wrap_in_closed_solid(plan);
    expect_error(geom::mass_properties(plan.editor.brep()), geom::Errc::invalid_topology,
                 "mass_properties: face has no loop");
  }
  {
    PlanarOptions options;
    options.plane = geom::Plane3{{0.0, 0.0, 0.0}, 0.0};
    PlanarFace plan = make_planar_face({{0.0, 0.0, 0.0}, {1.0, 0.0, 0.0}, {1.0, 1.0, 0.0}}, options);
    wrap_in_closed_solid(plan);
    expect_error(geom::mass_properties(plan.editor.brep()), geom::Errc::numeric_failure,
                 "mass_properties: face plane normal is not usable");
  }
  {
    PlanarOptions options;
    options.plane = geom::Plane3{{kNaN, 0.0, 1.0}, 0.0};
    PlanarFace plan = make_planar_face({{0.0, 0.0, 0.0}, {1.0, 0.0, 0.0}, {1.0, 1.0, 0.0}}, options);
    wrap_in_closed_solid(plan);
    expect_error(geom::mass_properties(plan.editor.brep()), geom::Errc::numeric_failure,
                 "mass_properties: face plane normal is not usable");
  }
  {
    // Face without surface geometry.
    geom::SolidEditor editor;
    const geom::VertexId a = editor.add_vertex({0.0, 0.0, 0.0}, 0.0).value();
    const geom::VertexId b = editor.add_vertex({1.0, 0.0, 0.0}, 0.0).value();
    const geom::VertexId c = editor.add_vertex({0.0, 1.0, 0.0}, 0.0).value();
    const geom::EdgeId ab = editor.add_edge(a, b, geom::kNoGeometry, 0.0, 1.0, 0.0).value();
    const geom::EdgeId bc = editor.add_edge(b, c, geom::kNoGeometry, 0.0, 1.0, 0.0).value();
    const geom::EdgeId ca = editor.add_edge(c, a, geom::kNoGeometry, 0.0, 1.0, 0.0).value();
    const geom::FaceId face = editor.add_face(geom::kNoGeometry, false, 0.0).value();
    add_loop_of(editor, face,
                {geom::CoedgeSpec{ab, false}, geom::CoedgeSpec{bc, false},
                 geom::CoedgeSpec{ca, false}});
    const geom::ShellId shell = editor.add_shell(std::array<geom::FaceId, 1>{face}, true).value();
    editor.add_solid(std::array<geom::ShellId, 1>{shell}, 0.0).value();
    expect_error(geom::mass_properties(editor.brep()), geom::Errc::index_out_of_range,
                 "mass_properties: face surface index out of range");
  }
  {
    // Face with an unsupported (NURBS) surface.
    geom::SolidEditor editor;
    const std::uint32_t surface = editor.add_surface(
        geom::NurbsSurface(1, 1, {0.0, 0.0, 1.0, 1.0}, {0.0, 0.0, 1.0, 1.0},
                           {{0.0, 0.0, 0.0}, {1.0, 0.0, 0.0}, {0.0, 1.0, 0.0}, {1.0, 1.0, 0.0}},
                           {1.0, 1.0, 1.0, 1.0}));
    const geom::VertexId a = editor.add_vertex({0.0, 0.0, 0.0}, 0.0).value();
    const geom::VertexId b = editor.add_vertex({1.0, 0.0, 0.0}, 0.0).value();
    const geom::VertexId c = editor.add_vertex({0.0, 1.0, 0.0}, 0.0).value();
    const geom::EdgeId ab = editor.add_edge(a, b, geom::kNoGeometry, 0.0, 1.0, 0.0).value();
    const geom::EdgeId bc = editor.add_edge(b, c, geom::kNoGeometry, 0.0, 1.0, 0.0).value();
    const geom::EdgeId ca = editor.add_edge(c, a, geom::kNoGeometry, 0.0, 1.0, 0.0).value();
    const geom::FaceId face = editor.add_face(surface, false, 0.0).value();
    add_loop_of(editor, face,
                {geom::CoedgeSpec{ab, false}, geom::CoedgeSpec{bc, false},
                 geom::CoedgeSpec{ca, false}});
    const geom::ShellId shell = editor.add_shell(std::array<geom::FaceId, 1>{face}, true).value();
    editor.add_solid(std::array<geom::ShellId, 1>{shell}, 0.0).value();
    expect_error(geom::mass_properties(editor.brep()), geom::Errc::unsupported,
                 "mass_properties: unsupported surface type");
  }
}

TEST_CASE("mass_properties integrates a closed cube and rejects inverted volume") {
  Cube cube = make_cube();
  const auto result = geom::mass_properties(cube.editor.brep());
  REQUIRE(result.has_value());
  const geom::MassProperties& mass = result.value();
  CHECK(mass.volume == doctest::Approx(1.0).epsilon(1e-12));
  CHECK(mass.centroid[0] == doctest::Approx(0.5).epsilon(1e-12));
  CHECK(mass.centroid[1] == doctest::Approx(0.5).epsilon(1e-12));
  CHECK(mass.centroid[2] == doctest::Approx(0.5).epsilon(1e-12));
  CHECK(mass.inertia[0][0] == doctest::Approx(1.0 / 6.0).epsilon(1e-9));
  CHECK(mass.inertia[1][1] == doctest::Approx(1.0 / 6.0).epsilon(1e-9));
  CHECK(mass.inertia[2][2] == doctest::Approx(1.0 / 6.0).epsilon(1e-9));
  CHECK(mass.inertia[0][1] == doctest::Approx(0.0).epsilon(1e-12));

  CubeOptions options;
  options.reversed = true;
  Cube inverted = make_cube(options);
  expect_error(geom::mass_properties(inverted.editor.brep()), geom::Errc::numeric_failure,
               "mass_properties: volume is not positive");
}

// =============================================================================
// tessellate.hpp
// =============================================================================

TEST_CASE("tessellate validates its options") {
  Square square = make_square();
  const geom::Brep& brep = square.editor.brep();
  for (const double deflection : {0.0, -1.0, kNaN, kInf}) {
    expect_error(geom::tessellate(brep, geom::TessellationOptions{deflection, 0.35}),
                 geom::Errc::tolerance_violation,
                 "tessellate: deflection must be finite and positive");
  }
  for (const double angle : {0.0, -1.0, kNaN, kInf}) {
    expect_error(geom::tessellate(brep, geom::TessellationOptions{1e-3, angle}),
                 geom::Errc::tolerance_violation,
                 "tessellate: angle tolerance must be finite and positive");
  }
}

TEST_CASE("tessellate rejects unsupported and empty input") {
  const geom::Brep empty;
  const auto mesh = geom::tessellate(empty, {});
  REQUIRE(mesh.has_value());
  CHECK(mesh.value().positions.empty());
  CHECK(mesh.value().triangles.empty());

  {
    geom::SolidEditor editor;
    const geom::FaceId face = editor.add_face(geom::kNoGeometry, false, 0.0).value();
    expect_error(geom::tessellate(editor.brep(), {}), geom::Errc::unsupported,
                 "tessellate: unsupported surface type");
    CHECK(face.index() == 0u);
  }
  {
    geom::SolidEditor editor;
    editor.add_surface(geom::NurbsSurface(1, 1, {0.0, 0.0, 1.0, 1.0}, {0.0, 0.0, 1.0, 1.0},
                                          {{0.0, 0.0, 0.0},
                                           {1.0, 0.0, 0.0},
                                           {0.0, 1.0, 0.0},
                                           {1.0, 1.0, 0.0}},
                                          {1.0, 1.0, 1.0, 1.0}));
    editor.add_face(0u, false, 0.0).value();
    expect_error(geom::tessellate(editor.brep(), {}), geom::Errc::unsupported,
                 "tessellate: unsupported surface type");
  }
}

TEST_CASE("tessellate rejects malformed planar faces") {
  const std::vector<geom::Point3d> triangle{{0.0, 0.0, 0.0}, {1.0, 0.0, 0.0}, {0.0, 1.0, 0.0}};
  {
    PlanarOptions options;
    options.plane = geom::Plane3{{0.0, 0.0, 1.0}, kNaN};
    PlanarFace plan = make_planar_face(triangle, options);
    expect_error(geom::tessellate(plan.editor.brep(), {}), geom::Errc::non_finite,
                 "tessellate: plane offset is not finite");
  }
  {
    PlanarOptions options;
    options.plane = geom::Plane3{{0.0, kNaN, 1.0}, 0.0};
    PlanarFace plan = make_planar_face(triangle, options);
    expect_error(geom::tessellate(plan.editor.brep(), {}), geom::Errc::non_finite,
                 "tessellate: plane normal is not finite");
  }
  {
    PlanarOptions options;
    options.plane = geom::Plane3{{0.0, 0.0, 0.0}, 0.0};
    PlanarFace plan = make_planar_face(triangle, options);
    expect_error(geom::tessellate(plan.editor.brep(), {}), geom::Errc::degenerate,
                 "tessellate: plane normal is zero");
  }
  {
    PlanarOptions options;
    options.with_loop = false;
    PlanarFace plan = make_planar_face(triangle, options);
    expect_error(geom::tessellate(plan.editor.brep(), {}), geom::Errc::degenerate,
                 "tessellate: face has no outer loop");
  }
  {
    PlanarOptions options;
    options.outer = false;
    PlanarFace plan = make_planar_face(triangle, options);
    expect_error(geom::tessellate(plan.editor.brep(), {}), geom::Errc::degenerate,
                 "tessellate: face has no outer loop");
  }
  {
    PlanarFace plan = make_planar_face(triangle);
    std::vector<geom::CoedgeSpec> second;
    for (const geom::EdgeId edge : plan.edges) second.push_back(geom::CoedgeSpec{edge, true});
    add_loop_of(plan.editor, plan.face, second, true);
    expect_error(geom::tessellate(plan.editor.brep(), {}), geom::Errc::invalid_topology,
                 "tessellate: face has more than one outer loop");
  }
  {
    PlanarFace plan = make_planar_face(
        {{0.0, 0.0, 0.0}, {1.0, 0.0, 0.0}, {2.0, 0.0, 0.0}, {3.0, 0.0, 0.0}});
    expect_error(geom::tessellate(plan.editor.brep(), {}), geom::Errc::degenerate,
                 "tessellate: outer loop cannot be projected");
  }
  {
    PlanarFace plan = make_planar_face(triangle);
    // Degenerate (collinear) hole ring.
    const geom::VertexId a = plan.editor.add_vertex({0.5, 0.0, 0.0}, 0.0).value();
    const geom::VertexId b = plan.editor.add_vertex({0.75, 0.0, 0.0}, 0.0).value();
    const geom::VertexId c = plan.editor.add_vertex({1.0, 0.0, 0.0}, 0.0).value();
    const geom::EdgeId ab = plan.editor.add_edge(a, b, geom::kNoGeometry, 0.0, 1.0, 0.0).value();
    const geom::EdgeId bc = plan.editor.add_edge(b, c, geom::kNoGeometry, 0.0, 1.0, 0.0).value();
    const geom::EdgeId ca = plan.editor.add_edge(c, a, geom::kNoGeometry, 0.0, 1.0, 0.0).value();
    add_loop_of(plan.editor, plan.face,
                {geom::CoedgeSpec{ab, false}, geom::CoedgeSpec{bc, false},
                 geom::CoedgeSpec{ca, false}},
                false);
    expect_error(geom::tessellate(plan.editor.brep(), {}), geom::Errc::degenerate,
                 "tessellate: hole loop cannot be projected");
  }
  {
    TwoEdgeFace two = make_two_edge_face();
    expect_error(geom::tessellate(two.editor.brep(), {}), geom::Errc::degenerate,
                 "tessellate: loop has fewer than three samples");
  }
}

TEST_CASE("tessellate rejects bad circle and non-analytic edge geometry") {
  const geom::CircleCurve3 circle{{0.0, 0.0, 0.0}, {1.0, 0.0, 0.0}, {0.0, 1.0, 0.0}, 1.0};
  {
    geom::CircleCurve3 broken = circle;
    broken.radius = 0.0;
    PlanarFace plan = make_circle_face(broken, 0.0, kTwoPi);
    expect_error(geom::tessellate(plan.editor.brep(), {}), geom::Errc::degenerate,
                 "tessellate: circle edge radius must be positive");
  }
  {
    geom::CircleCurve3 broken = circle;
    broken.center[1] = kNaN;
    // The vertices must stay finite so that the circle sampling is what fails.
    geom::SolidEditor editor;
    const std::uint32_t curve = editor.add_curve(broken);
    const geom::VertexId v0 = editor.add_vertex({0.0, 0.0, 0.0}, 0.0).value();
    const geom::VertexId v1 = editor.add_vertex({1.0, 0.0, 0.0}, 0.0).value();
    const geom::EdgeId edge = editor.add_edge(v0, v1, curve, 0.0, kTwoPi, 0.0).value();
    const geom::FaceId face = editor
                                  .add_face(editor.add_surface(geom::Plane3{{0.0, 0.0, 1.0}, 0.0}),
                                            false, 0.0)
                                  .value();
    add_loop_of(editor, face, {geom::CoedgeSpec{edge, false}});
    expect_error(geom::tessellate(editor.brep(), {}), geom::Errc::non_finite,
                 "tessellate: circle edge geometry is not finite");
  }
  {
    // A finite edge parameter interval whose span overflows.
    PlanarFace plan = make_circle_face(circle, -1e308, 1e308);
    expect_error(geom::tessellate(plan.editor.brep(), {}), geom::Errc::degenerate,
                 "tessellate: circle edge interval must be positive");
  }
  {
    geom::CircleCurve3 huge = circle;
    huge.radius = 1e308;
    PlanarFace plan = make_circle_face(huge, 0.0, kTwoPi);
    expect_error(geom::tessellate(plan.editor.brep(), geom::TessellationOptions{1e-9, 0.35}),
                 geom::Errc::numeric_failure,
                 "tessellate: circle segmentation count is too large");
  }
  {
    // A NURBS edge inside a planar loop.
    PlanarFace plan = make_planar_face(
        {{0.0, 0.0, 0.0}, {1.0, 0.0, 0.0}, {1.0, 1.0, 0.0}, {0.0, 1.0, 0.0}});
    const std::uint32_t curve = plan.editor.add_curve(
        geom::NurbsCurve(1, {0.0, 0.0, 1.0, 1.0}, {{0.0, 0.0, 0.0}, {1.0, 0.0, 0.0}}, {1.0, 1.0}));
    const geom::VertexId a = plan.editor.add_vertex({0.0, 0.0, 0.0}, 0.0).value();
    const geom::VertexId b = plan.editor.add_vertex({1.0, 1.0, 0.0}, 0.0).value();
    const geom::EdgeId edge = plan.editor.add_edge(a, b, curve, 0.0, 1.0, 0.0).value();
    const geom::FaceId face = plan.editor.add_face(
        plan.editor.add_surface(geom::Plane3{{0.0, 0.0, 1.0}, 0.0}), false, 0.0)
                                  .value();
    add_loop_of(plan.editor, face, {geom::CoedgeSpec{edge, false}});
    expect_error(geom::tessellate(plan.editor.brep(), {}), geom::Errc::unsupported,
                 "tessellate: only line and circle edges are supported");
  }
}

TEST_CASE("tessellate produces planar and holed meshes") {
  Square square = make_square();
  const auto mesh = geom::tessellate(square.editor.brep(), {});
  REQUIRE(mesh.has_value());
  CHECK(mesh.value().positions.size() == 4u);
  CHECK(mesh.value().triangles.size() == 2u);
  for (const std::uint32_t index : mesh.value().triangle_face) CHECK(index == 0u);

  PlanarFace plan = make_planar_face(
      {{0.0, 0.0, 0.0}, {4.0, 0.0, 0.0}, {4.0, 4.0, 0.0}, {0.0, 4.0, 0.0}});
  const geom::VertexId a = plan.editor.add_vertex({1.0, 1.0, 0.0}, 0.0).value();
  const geom::VertexId b = plan.editor.add_vertex({3.0, 1.0, 0.0}, 0.0).value();
  const geom::VertexId c = plan.editor.add_vertex({3.0, 3.0, 0.0}, 0.0).value();
  const geom::VertexId d = plan.editor.add_vertex({1.0, 3.0, 0.0}, 0.0).value();
  const geom::EdgeId ab = plan.editor.add_edge(a, b, geom::kNoGeometry, 0.0, 1.0, 0.0).value();
  const geom::EdgeId bc = plan.editor.add_edge(b, c, geom::kNoGeometry, 0.0, 1.0, 0.0).value();
  const geom::EdgeId cd = plan.editor.add_edge(c, d, geom::kNoGeometry, 0.0, 1.0, 0.0).value();
  const geom::EdgeId da = plan.editor.add_edge(d, a, geom::kNoGeometry, 0.0, 1.0, 0.0).value();
  add_loop_of(plan.editor, plan.face,
              {geom::CoedgeSpec{ab, false}, geom::CoedgeSpec{bc, false},
               geom::CoedgeSpec{cd, false}, geom::CoedgeSpec{da, false}},
              false);
  const auto holed = geom::tessellate(plan.editor.brep(), {});
  REQUIRE(holed.has_value());
  CHECK(holed.value().positions.size() == 8u);
  CHECK(holed.value().triangles.size() == 8u);
}

// =============================================================================
// Curved faces: the surface-level guards in tessellate and mass_properties fire
// before a face needs a valid curved boundary, so a bare triangular loop is
// enough to reach every "unsupported / degenerate curved face" outcome.
// =============================================================================

namespace {

PlanarFace make_curved_face(const geom::SurfaceGeometry& surface, bool with_loop = true) {
  PlanarFace result;
  result.surface = result.editor.add_surface(surface);
  const std::vector<geom::Point3d> pts{
      {0.0, 0.0, 0.0}, {1.0, 0.0, 0.0}, {0.0, 1.0, 0.0}};
  for (const geom::Point3d& point : pts) {
    result.vertices.push_back(result.editor.add_vertex(point, 0.0).value());
  }
  std::vector<geom::CoedgeSpec> specs;
  for (std::size_t i = 0; i < pts.size(); ++i) {
    result.edges.push_back(result.editor
                               .add_edge(result.vertices[i], result.vertices[(i + 1) % pts.size()],
                                         geom::kNoGeometry, 0.0, 1.0, 0.0)
                               .value());
    specs.push_back(geom::CoedgeSpec{result.edges.back(), false});
  }
  result.face = result.editor.add_face(result.surface, false, 0.0).value();
  if (with_loop) result.loop = result.editor.add_loop(result.face, specs, true).value();
  return result;
}

geom::SurfaceGeometry cylinder_surface(double radius) {
  return geom::Cylinder3{{0.0, 0.0, 0.0}, {0.0, 0.0, 1.0}, {1.0, 0.0, 0.0}, radius};
}
geom::SurfaceGeometry cone_surface(double radius) {
  return geom::Cone3{{0.0, 0.0, 0.0}, {0.0, 0.0, 1.0}, {1.0, 0.0, 0.0}, radius, kPi / 4.0};
}
geom::SurfaceGeometry sphere_surface(double radius) {
  return geom::Sphere3{{0.0, 0.0, 0.0}, {0.0, 0.0, 1.0}, {1.0, 0.0, 0.0}, radius};
}
geom::SurfaceGeometry torus_surface(double major, double minor) {
  return geom::Torus3{{0.0, 0.0, 0.0}, {0.0, 0.0, 1.0}, {1.0, 0.0, 0.0}, major, minor};
}

}  // namespace

TEST_CASE("tessellate rejects degenerate curved surfaces") {
  {
    PlanarFace face = make_curved_face(cylinder_surface(0.0));
    expect_error(geom::tessellate(face.editor.brep(), {}), geom::Errc::degenerate,
                 "tessellate: cylinder radius must be positive");
  }
  {
    geom::SurfaceGeometry surface = cylinder_surface(1.0);
    std::get<geom::Cylinder3>(surface).origin[0] = kNaN;
    PlanarFace face = make_curved_face(surface);
    expect_error(geom::tessellate(face.editor.brep(), {}), geom::Errc::non_finite,
                 "tessellate: cylinder surface is not finite");
  }
  {
    PlanarFace face = make_curved_face(cone_surface(0.0));
    expect_error(geom::tessellate(face.editor.brep(), {}), geom::Errc::degenerate,
                 "tessellate: cone surface is degenerate");
  }
  {
    PlanarFace face = make_curved_face(cone_surface(-1.0));
    expect_error(geom::tessellate(face.editor.brep(), {}), geom::Errc::degenerate,
                 "tessellate: cone surface is degenerate");
  }
  {
    geom::SurfaceGeometry surface = cone_surface(1.0);
    std::get<geom::Cone3>(surface).half_angle = kInf;
    PlanarFace face = make_curved_face(surface);
    expect_error(geom::tessellate(face.editor.brep(), {}), geom::Errc::non_finite,
                 "tessellate: cone surface is not finite");
  }
  {
    PlanarFace face = make_curved_face(sphere_surface(0.0));
    expect_error(geom::tessellate(face.editor.brep(), {}), geom::Errc::degenerate,
                 "tessellate: sphere radius must be positive");
  }
  {
    geom::SurfaceGeometry surface = sphere_surface(1.0);
    std::get<geom::Sphere3>(surface).center[2] = kNaN;
    PlanarFace face = make_curved_face(surface);
    expect_error(geom::tessellate(face.editor.brep(), {}), geom::Errc::non_finite,
                 "tessellate: sphere surface is not finite");
  }
  {
    PlanarFace face = make_curved_face(torus_surface(1.0, 1.0));
    expect_error(geom::tessellate(face.editor.brep(), {}), geom::Errc::degenerate,
                 "tessellate: torus radii are degenerate");
  }
  {
    geom::SurfaceGeometry surface = torus_surface(2.0, 1.0);
    std::get<geom::Torus3>(surface).axis[1] = kInf;
    PlanarFace face = make_curved_face(surface);
    expect_error(geom::tessellate(face.editor.brep(), {}), geom::Errc::non_finite,
                 "tessellate: torus surface is not finite");
  }
}

TEST_CASE("tessellate rejects flat loops on curved surfaces") {
  {
    PlanarFace face = make_curved_face(cylinder_surface(1.0));
    expect_error(geom::tessellate(face.editor.brep(), {}), geom::Errc::invalid_topology,
                 "tessellate: cylinder face is missing a cap circle");
  }
  {
    PlanarFace face = make_curved_face(cone_surface(1.0));
    expect_error(geom::tessellate(face.editor.brep(), {}), geom::Errc::invalid_topology,
                 "tessellate: cone face has too many seam edges");
  }
  {
    PlanarFace face = make_curved_face(sphere_surface(1.0));
    expect_error(geom::tessellate(face.editor.brep(), {}), geom::Errc::invalid_topology,
                 "tessellate: sphere face is missing its meridian");
  }
  {
    PlanarFace face = make_curved_face(torus_surface(2.0, 1.0));
    expect_error(geom::tessellate(face.editor.brep(), {}), geom::Errc::invalid_topology,
                 "tessellate: torus face is missing a seam");
  }
}

TEST_CASE("mass_properties rejects curved faces that are not full rectangles") {
  {
    PlanarFace face = make_curved_face(cylinder_surface(1.0));
    wrap_in_closed_solid(face);
    expect_error(geom::mass_properties(face.editor.brep()), geom::Errc::unsupported,
                 "mass_properties: curved face is not an untrimmed full rectangle");
  }
  {
    PlanarFace face = make_curved_face(cone_surface(1.0));
    wrap_in_closed_solid(face);
    expect_error(geom::mass_properties(face.editor.brep()), geom::Errc::unsupported,
                 "mass_properties: curved face is not an untrimmed full rectangle");
  }
  {
    PlanarFace face = make_curved_face(sphere_surface(1.0));
    wrap_in_closed_solid(face);
    expect_error(geom::mass_properties(face.editor.brep()), geom::Errc::unsupported,
                 "mass_properties: curved face is not an untrimmed full rectangle");
  }
  {
    PlanarFace face = make_curved_face(torus_surface(2.0, 1.0));
    wrap_in_closed_solid(face);
    expect_error(geom::mass_properties(face.editor.brep()), geom::Errc::unsupported,
                 "mass_properties: curved face is not an untrimmed full rectangle");
  }
}

// =============================================================================
// topology.hpp handle queries and SolidEditor::commit issue propagation
// =============================================================================

TEST_CASE("topology handle queries distinguish dead, live and null entities") {
  Square square = make_square();
  const geom::EdgeId edge = square.edges.front();
  const geom::VertexId vertex = square.vertices.front();
  const geom::Brep& brep = square.editor.brep();
  CHECK(brep.edge_alive(edge));
  CHECK(brep.vertex_alive(vertex));
  CHECK(brep.face_alive(square.face));
  CHECK(brep.loop_alive(square.loop));
  CHECK(brep.shell_alive(square.shell));

  REQUIRE(square.editor.erase_edge(edge).has_value());
  CHECK_FALSE(square.editor.brep().edge_alive(edge));
  CHECK(square.editor.brep().vertex_alive(vertex));
  REQUIRE(square.editor.erase_vertex(vertex).has_value());
  CHECK_FALSE(square.editor.brep().vertex_alive(vertex));

  CHECK_FALSE(square.editor.brep().vertex_alive(geom::VertexId{}));
  CHECK_FALSE(square.editor.brep().edge_alive(geom::EdgeId{}));
  CHECK_FALSE(square.editor.brep().coedge_alive(geom::CoedgeId{}));
  CHECK_FALSE(square.editor.brep().loop_alive(geom::LoopId{}));
  CHECK_FALSE(square.editor.brep().face_alive(geom::FaceId{}));
  CHECK_FALSE(square.editor.brep().shell_alive(geom::ShellId{}));
  CHECK_FALSE(square.editor.brep().solid_alive(geom::SolidId{}));
  CHECK(square.editor.brep().shell_faces(square.shell).size() == 1u);
}

TEST_CASE("solid_editor commit propagates the first validation issue") {
  {
    Cube cube = make_cube();
    const auto committed = cube.editor.commit();
    REQUIRE(committed.has_value());
    CHECK(committed.value().face_count() == 6u);
  }
  {
    CubeOptions options;
    options.face_tol[0] = 1e-3;
    Cube broken = make_cube(options);
    expect_error(broken.editor.commit(), geom::Errc::tolerance_violation,
                 "face tolerance exceeds the solid tolerance");
  }
  {
    Square square = make_square();
    REQUIRE(square.editor.erase_edge(square.edges.front()).has_value());
    expect_error(square.editor.commit(), geom::Errc::invalid_topology,
                 "coedge edge is not alive");
  }
}


TEST_CASE("curved faces without any loop are rejected") {
  for (const geom::SurfaceGeometry& surface :
       {cylinder_surface(1.0), cone_surface(1.0), sphere_surface(1.0), torus_surface(2.0, 1.0)}) {
    PlanarFace face = make_curved_face(surface, false);
    expect_error(geom::tessellate(face.editor.brep(), {}), geom::Errc::degenerate,
                 "tessellate: face has no loop");
  }
  for (const geom::SurfaceGeometry& surface :
       {cylinder_surface(1.0), cone_surface(1.0), sphere_surface(1.0), torus_surface(2.0, 1.0)}) {
    PlanarFace face = make_curved_face(surface, false);
    wrap_in_closed_solid(face);
    expect_error(geom::mass_properties(face.editor.brep()), geom::Errc::invalid_topology,
                 "mass_properties: face has no loop");
  }
}

TEST_CASE("tessellate rejects a hole larger than its outer loop") {
  PlanarFace plan = make_planar_face({{0.0, 0.0, 0.0}, {1.0, 0.0, 0.0}, {1.0, 1.0, 0.0}, {0.0, 1.0, 0.0}});
  const geom::Point3d big[4] = {{-2.0, -2.0, 0.0}, {3.0, -2.0, 0.0}, {3.0, 3.0, 0.0}, {-2.0, 3.0, 0.0}};
  std::array<geom::VertexId, 4> corners{};
  for (int i = 0; i < 4; ++i) corners[static_cast<std::size_t>(i)] = plan.editor.add_vertex(big[i], 0.0).value();
  const std::uint32_t curve = plan.editor.add_curve(geom::LineCurve3{{0.0, 0.0, 0.0}, {1.0, 0.0, 0.0}});
  std::vector<geom::CoedgeSpec> specs;
  for (int i = 0; i < 4; ++i) {
    const geom::EdgeId edge =
        plan.editor.add_edge(corners[static_cast<std::size_t>(i)],
                             corners[static_cast<std::size_t>((i + 1) % 4)], curve, 0.0, 1.0, 0.0).value();
    specs.push_back(geom::CoedgeSpec{edge, false});
  }
  add_loop_of(plan.editor, plan.face, specs, false);
  expect_error(geom::tessellate(plan.editor.brep(), {}), geom::Errc::degenerate,
               "tessellate: planar face has no area");
}
