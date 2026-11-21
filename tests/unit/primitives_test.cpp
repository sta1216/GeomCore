#include <doctest/doctest.h>

#include <cmath>
#include <cstdint>
#include <limits>
#include <numbers>
#include <string>
#include <variant>
#include <vector>

#include <geom/brep/detail/curved_volume.hpp>
#include <geom/brep/face_uv.hpp>
#include <geom/brep/primitives.hpp>
#include <geom/brep/tessellate.hpp>
#include <geom/brep/validate.hpp>
#include <geom/core/error.hpp>
#include <geom/core/result.hpp>
#include <geom/primitives/box.hpp>
#include <geom/primitives/point.hpp>
#include <geom/primitives/segment.hpp>
#include <geom/primitives/vec.hpp>
#include <geom/surfaces/analytic.hpp>
#include <geom/transform/affine.hpp>

using geom::v1::Affine;
using geom::v1::Box3d;
using geom::v1::Point3d;
using geom::v1::Segment3d;

namespace {

// contains() must stay usable in constant expressions (C++20): no isfinite().
constexpr Box3d kUnitBox{Point3d{0.0, 0.0, 0.0}, Point3d{1.0, 1.0, 1.0}};
static_assert(kUnitBox.contains(Point3d{0.5, 0.5, 0.5}));
static_assert(!kUnitBox.contains(Point3d{2.0, 0.5, 0.5}));
static_assert(!kUnitBox.contains(Point3d{std::numeric_limits<double>::infinity(), 0.5, 0.5}));
static_assert(!kUnitBox.contains(Point3d{std::numeric_limits<double>::quiet_NaN(), 0.5, 0.5}));

}  // namespace
using geom::v1::Vec3d;

TEST_CASE("point and vector affine algebra") {
  const Point3d p{1.0, 2.0, 3.0};
  const Vec3d v{0.5, -0.5, 1.0};
  const Point3d q = p + v;
  CHECK(q.v[0] == 1.5);
  CHECK(q.v[1] == 1.5);
  CHECK(q.v[2] == 4.0);

  const Vec3d d = q - p;
  CHECK(d.v[0] == v.v[0]);
  CHECK(geom::v1::dot(v, v) == doctest::Approx(1.5));
  CHECK(geom::v1::dot(geom::v1::cross(Vec3d{1, 0, 0}, Vec3d{0, 1, 0}), Vec3d{0, 0, 1}) ==
        doctest::Approx(1.0));
}

TEST_CASE("box is empty until expanded") {
  Box3d b = Box3d::empty();
  CHECK(b.is_empty());
  b.expand(Point3d{1.0, 2.0, 3.0});
  CHECK_FALSE(b.is_empty());
  CHECK(b.contains(Point3d{1.0, 2.0, 3.0}));
  b.expand(Point3d{-1.0, 0.0, 5.0});
  CHECK(b.lo.v[0] == -1.0);
  CHECK(b.hi.v[2] == 5.0);
}

TEST_CASE("closest point on a segment clamps to the endpoints") {
  const Segment3d s{Point3d{0, 0, 0}, Point3d{1, 0, 0}};
  CHECK(s.closest_point(Point3d{0.25, 1.0, 0.0}).v[0] == doctest::Approx(0.25));
  CHECK(s.closest_point(Point3d{-3.0, 0.0, 0.0}).v[0] == doctest::Approx(0.0));
  CHECK(s.closest_point(Point3d{3.0, 0.0, 0.0}).v[0] == doctest::Approx(1.0));
  CHECK(geom::v1::squared_distance(Point3d{0.5, 2.0, 0.0}, s) == doctest::Approx(4.0));
}

TEST_CASE("affine transforms round-trip through the inverse") {
  Affine a = Affine::identity();
  a = Affine::translation(Vec3d{1.0, 2.0, 3.0});
  const Point3d p{10.0, 20.0, 30.0};
  const Point3d q = a.transform(p);
  CHECK(q.v[0] == 11.0);

  const auto inv = a.inverse();
  REQUIRE(inv.has_value());
  const Point3d back = inv.value().transform(q);
  CHECK(back.v[0] == doctest::Approx(p.v[0]));
  CHECK(back.v[1] == doctest::Approx(p.v[1]));
  CHECK(back.v[2] == doctest::Approx(p.v[2]));
}

TEST_CASE("vector algebra and norm") {
  const Vec3d a{1.0, 2.0, 3.0};
  const Vec3d b{4.0, 5.0, 6.0};
  const Vec3d sum = a + b;
  CHECK(sum.v[0] == 5.0);
  CHECK(sum.v[2] == 9.0);
  const Vec3d diff = b - a;
  CHECK(diff.v[1] == 3.0);
  const Vec3d scaled = a * 2.0;
  CHECK(scaled.v[2] == 6.0);
  CHECK(geom::v1::norm2(geom::v1::Vec2d{3.0, 4.0}) == doctest::Approx(5.0));
  CHECK(geom::v1::dot(geom::v1::Vec2d{1.0, 2.0}, geom::v1::Vec2d{3.0, 4.0}) ==
        doctest::Approx(11.0));
  // cross is right-handed and anticommutative.
  CHECK(geom::v1::cross(Vec3d{1.0, 0.0, 0.0}, Vec3d{0.0, 1.0, 0.0}).v[2] == 1.0);
  CHECK(geom::v1::cross(Vec3d{0.0, 1.0, 0.0}, Vec3d{1.0, 0.0, 0.0}).v[2] == -1.0);
}

TEST_CASE("box merge and containment") {
  Box3d a = Box3d::empty();
  a.expand(Point3d{0.0, 0.0, 0.0});
  a.expand(Point3d{2.0, 4.0, 6.0});
  Box3d b = Box3d::empty();
  b.expand(Point3d{-1.0, 1.0, 1.0});
  b.expand(Point3d{1.0, 2.0, 2.0});

  CHECK_FALSE(Box3d::empty().contains(Point3d{0.0, 0.0, 0.0}));

  const Box3d m = a.merged(b);
  CHECK(m.lo.v[0] == -1.0);
  CHECK(m.hi.v[2] == 6.0);
  CHECK(m.diagonal().v[0] == doctest::Approx(3.0));
  CHECK(m.diagonal().v[2] == doctest::Approx(6.0));
  CHECK(m.contains(Point3d{0.5, 1.5, 2.0}));
  CHECK_FALSE(m.contains(Point3d{0.5, 5.0, 2.0}));

  // NaN fails both comparisons and must not be reported as contained.
  const double nan = std::numeric_limits<double>::quiet_NaN();
  const double inf = std::numeric_limits<double>::infinity();
  CHECK_FALSE(m.contains(Point3d{nan, 1.5, 2.0}));
  CHECK_FALSE(m.contains(Point3d{0.5, nan, 2.0}));
  CHECK_FALSE(m.contains(Point3d{0.5, 1.5, nan}));
  CHECK_FALSE(m.contains(Point3d{inf, 1.5, 2.0}));
  CHECK_FALSE(m.contains(Point3d{-inf, 1.5, 2.0}));

  // A box with infinite bounds is a valid "universe" box, but infinite points
  // are still not contained.
  const Box3d unbounded{Point3d{-inf, -inf, -inf}, Point3d{inf, inf, inf}};
  CHECK(unbounded.contains(Point3d{1e300, -1e300, 0.0}));
  CHECK_FALSE(unbounded.contains(Point3d{inf, 0.0, 0.0}));
  CHECK_FALSE(unbounded.contains(Point3d{-inf, 0.0, 0.0}));
  CHECK_FALSE(unbounded.contains(Point3d{nan, 0.0, 0.0}));
  CHECK_FALSE(Box3d::empty().contains(Point3d{inf, 0.0, 0.0}));
}

TEST_CASE("degenerate segment collapses to its point") {
  const Segment3d s{Point3d{1.0, 2.0, 3.0}, Point3d{1.0, 2.0, 3.0}};
  const Point3d c = s.closest_point(Point3d{4.0, 6.0, 3.0});
  CHECK(c.v[0] == 1.0);
  CHECK(c.v[1] == 2.0);
  CHECK(geom::v1::squared_distance(Point3d{4.0, 6.0, 3.0}, s) == doctest::Approx(25.0));
}

TEST_CASE("affine failure and vector paths") {
  // transform_vector ignores the translation column.
  const Affine t = Affine::translation(Vec3d{1.0, 2.0, 3.0});
  const Vec3d v = t.transform_vector(Vec3d{1.0, 1.0, 1.0});
  CHECK(v.v[0] == 1.0);
  CHECK(v.v[1] == 1.0);
  CHECK(v.v[2] == 1.0);
  // identity leaves points unchanged.
  const Point3d p{4.0, 5.0, 6.0};
  const Point3d q = Affine::identity().transform(p);
  CHECK(q.v[0] == 4.0);
  CHECK(q.v[2] == 6.0);

  // Singular linear part -> Err::degenerate, not a NaN transform.
  Affine singular{};
  for (int i = 0; i < 3; ++i) {
    for (int j = 0; j < 3; ++j) singular.m[i][j] = 0.0;
  }
  const auto inv = singular.inverse();
  CHECK_FALSE(inv.has_value());
  CHECK(inv.error().code == geom::v1::Errc::degenerate);
  CHECK(inv.error().message != nullptr);

  // A determinant that overflows to infinity is rejected the same way.
  Affine overflow{};
  overflow.m[0][0] = 1e308;
  overflow.m[1][1] = 1e308;
  overflow.m[2][2] = 1e308;
  const auto inv2 = overflow.inverse();
  CHECK_FALSE(inv2.has_value());
  CHECK(inv2.error().code == geom::v1::Errc::degenerate);
}

TEST_CASE("affine inverse rejects non-finite entries in every slot") {
  using geom::v1::Errc;
  const double nan_v = std::numeric_limits<double>::quiet_NaN();
  const double inf_v = std::numeric_limits<double>::infinity();

  auto expect_non_finite = [](const Affine& a) {
    const auto inv = a.inverse();
    CHECK_FALSE(inv.has_value());
    if (inv.has_value()) return;
    CHECK(inv.error().code == Errc::non_finite);
  };

  // Translation column slots.
  for (int row = 0; row < 3; ++row) {
    Affine a{};
    a.m[row][3] = inf_v;
    expect_non_finite(a);
    Affine b{};
    b.m[row][3] = nan_v;
    expect_non_finite(b);
  }
  // Diagonal slots.
  for (int row = 0; row < 3; ++row) {
    Affine a{};
    a.m[row][row] = -inf_v;
    expect_non_finite(a);
  }
  // Off-diagonal (rotation) slots.
  Affine rotation{};
  rotation.m[0][1] = nan_v;
  expect_non_finite(rotation);
  Affine rotation2{};
  rotation2.m[2][0] = inf_v;
  expect_non_finite(rotation2);
}

TEST_CASE("affine inverse finite-output contract") {
  using geom::v1::Errc;

  // A subnormal determinant makes 1/det overflow: rejected, never a NaN matrix.
  Affine subnormal{};
  subnormal.m[0][0] = 1e-320;
  const auto inv_sub = subnormal.inverse();
  REQUIRE_FALSE(inv_sub.has_value());
  CHECK(inv_sub.error().code == Errc::degenerate);

  // Determinant products that underflow to 0 or overflow to inf are conservative
  // rejects even though the exact inverse would be representable.
  Affine underflow{};
  underflow.m[0][0] = 1e-200;
  underflow.m[1][1] = 1e-200;
  underflow.m[2][2] = 1e-200;
  const auto inv_under = underflow.inverse();
  REQUIRE_FALSE(inv_under.has_value());
  CHECK(inv_under.error().code == Errc::degenerate);

  Affine overflow{};
  overflow.m[0][0] = 1e300;
  overflow.m[1][1] = 1e300;
  overflow.m[2][2] = 1e300;
  const auto inv_over = overflow.inverse();
  REQUIRE_FALSE(inv_over.has_value());
  CHECK(inv_over.error().code == Errc::degenerate);

  // A finite determinant can still yield an unrepresentable cofactor ratio.
  Affine mixed{};
  mixed.m[0][0] = 1e-300;
  mixed.m[0][1] = 1e300;
  mixed.m[1][1] = 1.0;
  mixed.m[2][0] = 1.0;
  mixed.m[2][2] = 1.0;
  const auto inv_mixed = mixed.inverse();
  REQUIRE_FALSE(inv_mixed.has_value());
  CHECK(inv_mixed.error().code == Errc::degenerate);

  // One extreme-but-finite scale factor keeps unit determinant: still invertible.
  Affine small_scale{};
  small_scale.m[0][0] = 1e-200;
  const auto inv_small = small_scale.inverse();
  REQUIRE(inv_small.has_value());
  CHECK(std::isfinite(inv_small.value().m[0][0]));
  CHECK(inv_small.value().m[1][1] == doctest::Approx(1.0));

  Affine big_scale{};
  big_scale.m[0][0] = 1e300;
  const auto inv_big = big_scale.inverse();
  REQUIRE(inv_big.has_value());
  CHECK(std::isfinite(inv_big.value().m[0][0]));
  CHECK(inv_big.value().m[0][0] == doctest::Approx(1e-300));
}

TEST_CASE("affine rotation and scale inverse round-trips") {
  Affine a{};
  a.m[0][0] = 0.0;
  a.m[0][1] = -2.0;
  a.m[0][3] = 5.0;
  a.m[1][0] = 2.0;
  a.m[1][1] = 0.0;
  a.m[1][3] = -3.0;
  a.m[2][2] = 1.0;
  a.m[2][3] = 2.0;

  const auto inv = a.inverse();
  REQUIRE(inv.has_value());
  const Point3d p{1.0, 2.0, 3.0};
  const Point3d q = a.transform(p);
  const Point3d back = inv.value().transform(q);
  CHECK(back.v[0] == doctest::Approx(p.v[0]).epsilon(1e-12));
  CHECK(back.v[1] == doctest::Approx(p.v[1]).epsilon(1e-12));
  CHECK(back.v[2] == doctest::Approx(p.v[2]).epsilon(1e-12));

  const Vec3d v = a.transform_vector(Vec3d{1.0, 0.0, 0.0});
  CHECK(v.v[0] == 0.0);
  CHECK(v.v[1] == doctest::Approx(2.0));
}

namespace {

double primitive_relative_error(double actual, double expected) {
  const double scale = std::fabs(expected) > 1.0 ? std::fabs(expected) : 1.0;
  return std::fabs(actual - expected) / scale;
}

}  // namespace

// Task 2 of M3c-3b: the sphere primitive must accept an arbitrary orthonormal
// chart so that a boolean can re-chart an operand in the shared frame. The
// sphere below is re-charted with axis = +x and x_axis = +y; the face surface
// must carry exactly that frame (the seam circle and the poles follow from it).
TEST_CASE("make_sphere accepts an arbitrary chart") {
  const auto sphere = geom::make_sphere(geom::Point3d{0.0, 0.0, 0.0}, geom::Vec3d{1.0, 0.0, 0.0},
                                        geom::Vec3d{0.0, 1.0, 0.0}, 1.0, 1e-9);
  REQUIRE(sphere.has_value());
  const geom::Face& face = sphere.value().face(sphere.value().face_ids().front());
  const geom::Sphere3* surface =
      std::get_if<geom::Sphere3>(&sphere.value().geometry().surfaces[face.surface]);
  REQUIRE(surface != nullptr);
  CHECK(surface->axis[0] == doctest::Approx(1.0));
  CHECK(surface->axis[1] == doctest::Approx(0.0));
  CHECK(surface->axis[2] == doctest::Approx(0.0));
  CHECK(surface->x_axis[0] == doctest::Approx(0.0));
  CHECK(surface->x_axis[1] == doctest::Approx(1.0));
  CHECK(surface->x_axis[2] == doctest::Approx(0.0));

  // The seam edge must be built from the same chart as the surface; a seam
  // left on the canonical frame would silently tear the chart.
  const geom::Loop& loop = sphere.value().loop(face.first_loop);
  const geom::Coedge& coedge = sphere.value().coedge(loop.start);
  const geom::Edge& edge = sphere.value().edge(coedge.edge);
  const auto* seam = std::get_if<geom::CircleCurve3>(
      &sphere.value().geometry().curves[edge.curve]);
  REQUIRE(seam != nullptr);
  CHECK(seam->x_axis[0] == doctest::Approx(0.0));
  CHECK(seam->x_axis[1] == doctest::Approx(1.0));
  CHECK(seam->x_axis[2] == doctest::Approx(0.0));
  CHECK(seam->y_axis[0] == doctest::Approx(1.0));
  CHECK(seam->y_axis[1] == doctest::Approx(0.0));
  CHECK(seam->y_axis[2] == doctest::Approx(0.0));
}

// Regression guard for the canonical overload: it must delegate to the frame
// overload without changing topology, chart, poles or volume.
TEST_CASE("make_sphere canonical chart is unchanged") {
  const auto sphere = geom::make_sphere(geom::Point3d{0.0, 0.0, 0.0}, 1.0, 1e-9);
  REQUIRE(sphere.has_value());
  const geom::Brep& brep = sphere.value();

  CHECK(brep.vertex_count() == 2);
  CHECK(brep.edge_count() == 1);
  CHECK(brep.face_count() == 1);
  CHECK(brep.shell_count() == 1);
  CHECK(brep.solid_count() == 1);

  const geom::Face& face = brep.face(brep.face_ids().front());
  const geom::Sphere3* surface =
      std::get_if<geom::Sphere3>(&brep.geometry().surfaces[face.surface]);
  REQUIRE(surface != nullptr);
  CHECK(surface->axis[0] == 0.0);
  CHECK(surface->axis[1] == 0.0);
  CHECK(surface->axis[2] == 1.0);
  CHECK(surface->x_axis[0] == 1.0);
  CHECK(surface->x_axis[1] == 0.0);
  CHECK(surface->x_axis[2] == 0.0);

  // The two vertices are the poles on the axis: (0, 0, -1) and (0, 0, +1).
  const std::vector<geom::VertexId> vertices = brep.vertex_ids();
  REQUIRE(vertices.size() == 2);
  bool saw_south = false;
  bool saw_north = false;
  for (const geom::VertexId id : vertices) {
    const Point3d& p = brep.vertex(id).point;
    CHECK(p[0] == 0.0);
    CHECK(p[1] == 0.0);
    if (p[2] == -1.0) saw_south = true;
    if (p[2] == 1.0) saw_north = true;
  }
  CHECK(saw_south);
  CHECK(saw_north);

  const geom::Result<double> volume =
      geom::detail::curved_solid_volume(brep, brep.solid_ids().front(), 1e-9);
  REQUIRE(volume.has_value());
  const double expected = 4.0 * std::numbers::pi_v<double> / 3.0;
  CHECK(primitive_relative_error(volume.value(), expected) <= 1e-12);
}

namespace {

// Pin both the error class and the exact message of a rejected sphere. The
// canonical overload delegates to the frame overload, so this contract is the
// old one: messages and codes must not drift when a frame is rejected.
void check_sphere_error(const geom::Result<geom::Brep>& result, geom::Errc code,
                        const char* message) {
  REQUIRE_FALSE(result.has_value());
  CHECK(result.error().code == code);
  CHECK(std::string(result.error().message) == std::string(message));
}

}  // namespace

// Frame-aware sphere: every rejection path, with the exact code and message.
TEST_CASE("make_sphere rejects an invalid chart") {
  const geom::Point3d origin{0.0, 0.0, 0.0};
  const geom::Vec3d z{0.0, 0.0, 1.0};
  const geom::Vec3d x{1.0, 0.0, 0.0};
  const double inf = std::numeric_limits<double>::infinity();
  const double nan = std::numeric_limits<double>::quiet_NaN();

  // Non-finite center, axis, x_axis and radius all take the first branch.
  // The 1e-12 frame bound is inclusive: 1e-10 and 2e-12 must be rejected so a
  // widened bound cannot pass an invalid Sphere3 past the library invariant.
  check_sphere_error(geom::make_sphere(origin, geom::Vec3d{1.0 + 1e-10, 0.0, 0.0}, x, 1.0, 1e-9),
                     geom::Errc::degenerate, "make_sphere: axis must be unit length");
  check_sphere_error(geom::make_sphere(origin, z, geom::Vec3d{1.0 + 1e-10, 0.0, 0.0}, 1.0, 1e-9),
                     geom::Errc::degenerate,
                     "make_sphere: x-axis must be perpendicular to the axis");
  check_sphere_error(geom::make_sphere(origin, z, geom::Vec3d{1.0, 0.0, 2e-12}, 1.0, 1e-9),
                     geom::Errc::degenerate,
                     "make_sphere: x-axis must be perpendicular to the axis");
  check_sphere_error(geom::make_sphere(geom::Point3d{inf, 0.0, 0.0}, z, x, 1.0, 1e-9),
                     geom::Errc::non_finite, "make_sphere: input is not finite");
  check_sphere_error(geom::make_sphere(geom::Point3d{nan, 0.0, 0.0}, z, x, 1.0, 1e-9),
                     geom::Errc::non_finite, "make_sphere: input is not finite");
  check_sphere_error(geom::make_sphere(origin, geom::Vec3d{0.0, nan, 0.0}, x, 1.0, 1e-9),
                     geom::Errc::non_finite, "make_sphere: input is not finite");
  check_sphere_error(geom::make_sphere(origin, z, geom::Vec3d{inf, 0.0, 0.0}, 1.0, 1e-9),
                     geom::Errc::non_finite, "make_sphere: input is not finite");
  check_sphere_error(geom::make_sphere(origin, z, x, nan, 1e-9), geom::Errc::non_finite,
                     "make_sphere: input is not finite");

  // A NaN tolerance is a tolerance violation, not a non-finite input: the
  // finite branch above covers center/axis/x_axis/radius only.
  check_sphere_error(geom::make_sphere(origin, z, x, 1.0, nan),
                     geom::Errc::tolerance_violation,
                     "make_sphere: tolerance must be finite and non-negative");
  check_sphere_error(geom::make_sphere(origin, z, x, 1.0, -1.0),
                     geom::Errc::tolerance_violation,
                     "make_sphere: tolerance must be finite and non-negative");
  check_sphere_error(geom::make_sphere(origin, z, x, 0.0, 1e-9), geom::Errc::degenerate,
                     "make_sphere: radius must be positive");
  check_sphere_error(geom::make_sphere(origin, z, x, -1.0, 1e-9), geom::Errc::degenerate,
                     "make_sphere: radius must be positive");
  check_sphere_error(geom::make_sphere(origin, geom::Vec3d{0.0, 0.0, 1.001}, x, 1.0, 1e-9),
                     geom::Errc::degenerate, "make_sphere: axis must be unit length");
  check_sphere_error(geom::make_sphere(origin, geom::Vec3d{0.0, 0.0, 2.0}, x, 1.0, 1e-9),
                     geom::Errc::degenerate, "make_sphere: axis must be unit length");
  check_sphere_error(geom::make_sphere(origin, geom::Vec3d{0.0, 0.0, 0.0}, x, 1.0, 1e-9),
                     geom::Errc::degenerate, "make_sphere: axis must be unit length");
  check_sphere_error(geom::make_sphere(origin, z, geom::Vec3d{1.001, 0.0, 0.0}, 1.0, 1e-9),
                     geom::Errc::degenerate,
                     "make_sphere: x-axis must be perpendicular to the axis");
  check_sphere_error(geom::make_sphere(origin, geom::Vec3d{1.0, 0.0, 0.0},
                                       geom::Vec3d{1.0, 0.0, 0.0}, 1.0, 1e-9),
                     geom::Errc::degenerate,
                     "make_sphere: x-axis must be perpendicular to the axis");

  // The canonical overload delegates, so it must reject with the same contract.
  check_sphere_error(geom::make_sphere(geom::Point3d{inf, 0.0, 0.0}, 1.0, 1e-9),
                     geom::Errc::non_finite, "make_sphere: input is not finite");
}

namespace {

double point_distance(const geom::Point3d& a, const geom::Point3d& b) {
  return std::hypot(a[0] - b[0], a[1] - b[1], a[2] - b[2]);
}

// FNV-1a over the mesh arrays, the same hash the tessellation stability test
// pins for the canonical primitive builders.
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

}  // namespace

// A shifted chart moves the u seam onto the surface curve with old v = v_origin,
// so the face's chart seam never crosses a cut that avoids that circle. The u
// seam circle must lie on the shifted surface and the v seam must still be the
// tube circle at u = 0.
TEST_CASE("make_torus with a v-origin re-seats the u seam on the shifted chart seam") {
  constexpr double kPi = 3.14159265358979323846;
  constexpr double kTwoPi = 2.0 * kPi;
  constexpr double kOrigin = 0.7;
  constexpr double kTolerance = 1e-9;
  constexpr double kMajor = 1.0;
  constexpr double kMinor = 0.5;
  const geom::Point3d center{0.25, -0.5, 1.5};
  const geom::Vec3d axis{0.0, 0.0, 1.0};

  const auto result = geom::make_torus(center, axis, kMajor, kMinor, kTolerance, kOrigin);
  REQUIRE(result.has_value());
  const geom::Brep& brep = result.value();
  CHECK(brep.vertex_count() == 1);
  CHECK(brep.edge_count() == 2);
  CHECK(brep.face_count() == 1);
  CHECK(brep.coedge_count() == 4);
  CHECK(geom::validate_topology(brep).ok());

  const geom::Face& face = brep.face(brep.face_ids().front());
  const geom::Torus3* torus = std::get_if<geom::Torus3>(&brep.geometry().surfaces[face.surface]);
  REQUIRE(torus != nullptr);
  CHECK(torus->v_origin == kOrigin);
  CHECK(torus->is_valid());

  const double cos_origin = std::cos(kOrigin);
  const double sin_origin = std::sin(kOrigin);
  const double seam_radius = kMajor + kMinor * cos_origin;
  const geom::Point3d seam_center{center[0], center[1], center[2] + kMinor * sin_origin};

  // The single vertex is the chart origin: old (u, v) = (0, v_origin).
  const geom::VertexId vertex_id = brep.vertex_ids().front();
  const geom::Point3d vertex = brep.vertex(vertex_id).point;
  CHECK(vertex[0] == doctest::Approx(center[0] + seam_radius).epsilon(1e-12));
  CHECK(vertex[1] == doctest::Approx(center[1]).epsilon(1e-12));
  CHECK(vertex[2] == doctest::Approx(seam_center[2]).epsilon(1e-12));

  // Every vertex is on the shifted surface: project it, evaluate back.
  const auto uv = geom::detail::face_uv_project_torus(*torus, vertex);
  REQUIRE(uv.has_value());
  CHECK(point_distance(torus->eval(uv.value()[0], uv.value()[1]), vertex) <= 1e-12);

  geom::EdgeId u_edge{};
  geom::EdgeId v_edge{};
  for (const geom::EdgeId id : brep.edge_ids()) {
    const auto& circle = std::get<geom::CircleCurve3>(brep.geometry().curves[brep.edge(id).curve]);
    if (std::fabs(circle.radius - kMinor) < 1e-12) {
      v_edge = id;
    } else {
      u_edge = id;
    }
  }
  REQUIRE_FALSE(u_edge.is_null());
  REQUIRE_FALSE(v_edge.is_null());

  const geom::Edge& u = brep.edge(u_edge);
  const auto& u_circle = std::get<geom::CircleCurve3>(brep.geometry().curves[u.curve]);
  CHECK(u.v0 == vertex_id);
  CHECK(u.v1 == vertex_id);
  CHECK(u.t0 == 0.0);
  CHECK(u.t1 == doctest::Approx(kTwoPi).epsilon(1e-12));
  CHECK(u_circle.radius == doctest::Approx(seam_radius).epsilon(1e-12));
  for (int i = 0; i < 3; ++i) {
    CHECK(u_circle.center[i] == doctest::Approx(seam_center[i]).epsilon(1e-12));
    CHECK(u_circle.x_axis[i] == doctest::Approx(i == 0 ? 1.0 : 0.0).epsilon(1e-12));
    CHECK(u_circle.y_axis[i] == doctest::Approx(i == 1 ? 1.0 : 0.0).epsilon(1e-12));
  }
  for (const double t : {0.0, 0.4, kPi / 2.0, 2.5, 4.0}) {
    const geom::Point3d on_circle = u_circle.eval(t);
    CHECK(point_distance(on_circle, torus->eval(t, 0.0)) <= 1e-12);
    const auto seam_uv = geom::detail::face_uv_project_torus(*torus, on_circle);
    REQUIRE(seam_uv.has_value());
    // The u seam is the chart curve v = 0.
    CHECK(std::fabs(seam_uv.value()[1]) <= 1e-12);
  }

  const geom::Edge& v = brep.edge(v_edge);
  const auto& v_circle = std::get<geom::CircleCurve3>(brep.geometry().curves[v.curve]);
  CHECK(v.v0 == vertex_id);
  CHECK(v.v1 == vertex_id);
  // The tube circle stays at u = 0; the interval starts at the origin angle so
  // that t0 maps to the shared vertex.
  CHECK(v.t0 == kOrigin);
  CHECK(v.t1 == doctest::Approx(kOrigin + kTwoPi).epsilon(1e-12));
  CHECK(v_circle.radius == kMinor);
  const geom::Point3d tube_center{center[0] + kMajor, center[1], center[2]};
  for (int i = 0; i < 3; ++i) {
    CHECK(v_circle.center[i] == doctest::Approx(tube_center[i]).epsilon(1e-12));
    CHECK(v_circle.x_axis[i] == doctest::Approx(i == 0 ? 1.0 : 0.0).epsilon(1e-12));
    CHECK(v_circle.y_axis[i] == doctest::Approx(axis[i]).epsilon(1e-12));
  }
  CHECK(point_distance(u_circle.eval(u.t0), v_circle.eval(v.t0)) <= 1e-12);
  for (const double w : {0.0, 0.3, 1.1, 1.9}) {
    const geom::Point3d on_circle = v_circle.eval(v.t0 + w);
    CHECK(point_distance(on_circle, torus->eval(0.0, w)) <= 1e-12);
    const auto tube_uv = geom::detail::face_uv_project_torus(*torus, on_circle);
    REQUIRE(tube_uv.has_value());
    CHECK(std::fabs(tube_uv.value()[0]) <= 1e-12);
    CHECK(tube_uv.value()[1] == doctest::Approx(w).epsilon(1e-12));
  }

  // The projection subtracts v_origin without wrapping: the far side of the
  // chart lands below zero instead of being folded back into [0, 2*pi).
  const geom::Point3d past = torus->eval(0.4, 3.6);
  const auto past_uv = geom::detail::face_uv_project_torus(*torus, past);
  REQUIRE(past_uv.has_value());
  CHECK(past_uv.value()[0] == doctest::Approx(0.4).epsilon(1e-12));
  CHECK(past_uv.value()[1] == doctest::Approx(3.6 - kTwoPi).epsilon(1e-12));

  // The enclosed volume is origin-independent: pi^2 / 2 for R = 1, r = 1/2.
  const geom::Result<double> volume =
      geom::detail::curved_solid_volume(brep, brep.solid_ids().front(), kTolerance);
  REQUIRE(volume.has_value());
  CHECK(std::fabs(volume.value() - 4.934802200544679) <= 1e-9);

  const double inf = std::numeric_limits<double>::infinity();
  const auto infinite_origin = geom::make_torus(center, axis, kMajor, kMinor, kTolerance, inf);
  REQUIRE_FALSE(infinite_origin.has_value());
  CHECK(infinite_origin.error().code == geom::Errc::non_finite);
  CHECK(std::string(infinite_origin.error().message) == "make_torus: input is not finite");
  CHECK_FALSE(geom::make_torus(center, axis, kMajor, kMinor, kTolerance,
                               std::numeric_limits<double>::quiet_NaN())
                  .has_value());
}

// An explicit zero v-origin must leave the canonical primitive untouched: same
// topology, same two full-turn seam circles, byte-identical tessellation.
TEST_CASE("make_torus with a zero v-origin keeps the canonical output") {
  constexpr double kPi = 3.14159265358979323846;
  constexpr double kTwoPi = 2.0 * kPi;
  const geom::Point3d center{0.0, 0.0, 0.0};
  const geom::Vec3d axis{0.0, 0.0, 1.0};
  const auto implicit_origin = geom::make_torus(center, axis, 2.0, 0.5, 1e-9);
  const auto explicit_origin = geom::make_torus(center, axis, 2.0, 0.5, 1e-9, 0.0);
  REQUIRE(implicit_origin.has_value());
  REQUIRE(explicit_origin.has_value());
  const geom::Brep& brep = explicit_origin.value();
  CHECK(brep.vertex_count() == 1);
  CHECK(brep.edge_count() == 2);
  CHECK(brep.face_count() == 1);
  CHECK(brep.coedge_count() == 4);
  CHECK(brep.shell_count() == 1);
  CHECK(brep.solid_count() == 1);

  const geom::Face& face = brep.face(brep.face_ids().front());
  const geom::Torus3* torus = std::get_if<geom::Torus3>(&brep.geometry().surfaces[face.surface]);
  REQUIRE(torus != nullptr);
  CHECK(torus->v_origin == 0.0);

  bool saw_outer = false;
  bool saw_tube = false;
  for (const geom::EdgeId id : brep.edge_ids()) {
    const geom::Edge& edge = brep.edge(id);
    const auto& circle = std::get<geom::CircleCurve3>(brep.geometry().curves[edge.curve]);
    CHECK(edge.t0 == 0.0);
    CHECK(edge.t1 == doctest::Approx(kTwoPi).epsilon(1e-12));
    if (circle.radius == 2.5) saw_outer = true;
    if (circle.radius == 0.5) saw_tube = true;
  }
  CHECK(saw_outer);
  CHECK(saw_tube);

  const geom::TessellationOptions options{1e-2, 0.35};
  const auto implicit_mesh = geom::tessellate(implicit_origin.value(), options);
  const auto explicit_mesh = geom::tessellate(explicit_origin.value(), options);
  REQUIRE(implicit_mesh.has_value());
  REQUIRE(explicit_mesh.has_value());
  CHECK(explicit_mesh.value().triangles.size() == 1152u);
  CHECK(explicit_mesh.value().positions.size() == 576u);
  CHECK(hash_mesh(explicit_mesh.value()) == hash_mesh(implicit_mesh.value()));
#if defined(_MSC_VER)
  CHECK(hash_mesh(explicit_mesh.value()) == 0xc808f6ee4d9c72dbull);
#endif
}
