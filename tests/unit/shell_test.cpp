#include <doctest/doctest.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <limits>
#include <span>
#include <string>
#include <string_view>
#include <variant>
#include <vector>

#include <geom/brep/assemble.hpp>
#include <geom/brep/mass_properties.hpp>
#include <geom/brep/primitives.hpp>
#include <geom/brep/shell.hpp>
#include <geom/brep/validate.hpp>
#include <geom/core/error.hpp>

using namespace geom;

namespace {

Brep box_solid(Point3d lo, Point3d hi) {
  const Result<Brep> box = make_box(lo, hi, 1e-9);
  REQUIRE(box.has_value());
  return box.value();
}

double volume_of(const Brep& brep) {
  const Result<MassProperties> mass = mass_properties(brep);
  REQUIRE(mass.has_value());
  return mass.value().volume;
}

}  // namespace

TEST_CASE("closed shell leaves a uniform cavity") {
  const Brep box = box_solid(Point3d{0.0, 0.0, 0.0}, Point3d{2.0, 3.0, 4.0});
  ShellOptions options;
  options.thickness = 0.25;
  const std::array<FaceId, 0> open{};
  const Result<ShellResult> result = shell_solid(box, open, options);
  REQUIRE(result.has_value());
  CHECK(result.value().cavity_shells == 1u);
  CHECK(result.value().brep.solid_count() == 1u);
  const auto shells = result.value().brep.solid_shells(result.value().brep.solid_ids().front());
  CHECK(shells.size() == 2u);
  CHECK(validate_topology(result.value().brep).ok());
  // 24 - (1.5 * 2.5 * 3.5).
  CHECK(std::abs(volume_of(result.value().brep) - 10.875) <= 1e-9);
}

TEST_CASE("shell rejects a thickness that leaves no interior") {
  const Brep box = box_solid(Point3d{0.0, 0.0, 0.0}, Point3d{2.0, 3.0, 4.0});
  ShellOptions options;
  options.thickness = 1.0;
  const std::array<FaceId, 0> open{};
  const Result<ShellResult> result = shell_solid(box, open, options);
  REQUIRE_FALSE(result.has_value());
  CHECK(result.error().code == Errc::degenerate);
  CHECK(std::string_view(result.error().message).find("leaves no interior") !=
        std::string::npos);
}

TEST_CASE("shell refuses an open face") {
  const Brep box = box_solid(Point3d{0.0, 0.0, 0.0}, Point3d{2.0, 3.0, 4.0});
  FaceId top{};
  for (const FaceId face : box.face_ids()) {
    const auto& plane = std::get<Plane3>(box.geometry().surfaces[box.face(face).surface]);
    const Vec3d normal = box.face(face).reversed ? plane.normal * -1.0 : plane.normal;
    if (normal[2] > 0.5) top = face;
  }
  REQUIRE_FALSE(top.is_null());
  ShellOptions options;
  options.thickness = 0.25;
  const std::array<FaceId, 1> open{top};
  const Result<ShellResult> result = shell_solid(box, open, options);
  REQUIRE_FALSE(result.has_value());
  CHECK(result.error().code == Errc::unsupported);
  CHECK(std::string_view(result.error().message) ==
        "shell_solid: open faces are not supported in M4c-1");
}

TEST_CASE("closed axis-aligned shell is independent of face order") {
  const auto p = [](double x, double y, double z) { return Point3d{x, y, z}; };
  const std::array<AssembledFace, 6> faces{{
      {{p(0, 0, 0), p(0, 0, 4), p(0, 3, 4), p(0, 3, 0)}, {}, {-1, 0, 0}},
      {{p(2, 0, 0), p(2, 3, 0), p(2, 3, 4), p(2, 0, 4)}, {}, {1, 0, 0}},
      {{p(0, 0, 0), p(2, 0, 0), p(2, 0, 4), p(0, 0, 4)}, {}, {0, -1, 0}},
      {{p(0, 3, 0), p(0, 3, 4), p(2, 3, 4), p(2, 3, 0)}, {}, {0, 1, 0}},
      {{p(0, 0, 0), p(0, 3, 0), p(2, 3, 0), p(2, 0, 0)}, {}, {0, 0, -1}},
      {{p(0, 0, 4), p(2, 0, 4), p(2, 3, 4), p(0, 3, 4)}, {}, {0, 0, 1}},
  }};
  std::array<std::size_t, 6> order{0, 1, 2, 3, 4, 5};
  ShellOptions options;
  options.thickness = 0.25;
  const std::array<FaceId, 0> closed{};
  // All 720 orders in Release; a fixed 24-order prefix keeps the unoptimized
  // Debug suite bounded while still exercising several permutations.
#ifdef NDEBUG
  constexpr std::size_t target = 720;
#else
  constexpr std::size_t target = 24;
#endif
  std::size_t checked = 0;
  std::string_view first_failure;
  do {
    std::vector<AssembledFace> input;
    input.reserve(faces.size());
    for (const std::size_t index : order) input.push_back(faces[index]);
    const Result<Brep> solid = assemble_solid(input, 1e-9);
    if (!solid.has_value()) {
      first_failure = solid.error().message;
      break;
    }
    const Result<ShellResult> result = shell_solid(solid.value(), closed, options);
    if (!result.has_value()) {
      first_failure = result.error().message;
      break;
    }
    if (std::abs(volume_of(result.value().brep) - 10.875) > 1e-9) {
      first_failure = "wrong volume";
      break;
    }
    ++checked;
  } while (checked < target && std::next_permutation(order.begin(), order.end()));
  INFO("permutations checked: " << checked << " of " << target
                                << "; first failure: " << first_failure);
  CHECK(checked == target);
}

TEST_CASE("shell refuses a wall below the boolean resolution") {
  const Brep box = box_solid(Point3d{0.0, 0.0, 0.0}, Point3d{2.0, 3.0, 4.0});
  ShellOptions options;
  options.thickness = 1e-8;
  const std::array<FaceId, 0> open{};
  const Result<ShellResult> result = shell_solid(box, open, options);
  REQUIRE_FALSE(result.has_value());
  CHECK(result.error().code == Errc::degenerate);
  CHECK(std::string_view(result.error().message).find("no measurable interior") !=
        std::string::npos);
}

TEST_CASE("shell rejects a cavity thinner than the tolerance") {
  const Brep box = box_solid(Point3d{0.0, 0.0, 0.0}, Point3d{2.0, 3.0, 4.0});
  ShellOptions options;
  options.thickness = 0.999999999;
  const std::array<FaceId, 0> open{};
  const Result<ShellResult> result = shell_solid(box, open, options);
  REQUIRE_FALSE(result.has_value());
  CHECK(result.error().code == Errc::degenerate);
  CHECK(std::string_view(result.error().message).find("no measurable interior") !=
        std::string::npos);
}

TEST_CASE("shell rejects opening every face") {
  const Brep box = box_solid(Point3d{0.0, 0.0, 0.0}, Point3d{2.0, 3.0, 4.0});
  const auto faces = box.face_ids();
  ShellOptions options;
  options.thickness = 0.25;
  const Result<ShellResult> result = shell_solid(box, faces, options);
  REQUIRE_FALSE(result.has_value());
  CHECK(result.error().code == Errc::unsupported);
  CHECK(std::string_view(result.error().message) ==
        "shell_solid: open faces are not supported in M4c-1");
}

TEST_CASE("shell refuses an oblique face") {
  const double s = std::sqrt(2.0);
  const auto p = [](double x, double y, double z) { return Point3d{x, y, z}; };
  const AssembledFace bottom{{p(0, 0, 0), p(0, 1.5, 0), p(1.5, 0, 0)}, {}, {0, 0, -1}};
  const AssembledFace top{{p(0, 0, 4), p(1.5, 0, 4), p(0, 1.5, 4)}, {}, {0, 0, 1}};
  const AssembledFace front{{p(0, 0, 0), p(1.5, 0, 0), p(1.5, 0, 4), p(0, 0, 4)}, {}, {0, -1, 0}};
  const AssembledFace diagonal{
      {p(1.5, 0, 0), p(0, 1.5, 0), p(0, 1.5, 4), p(1.5, 0, 4)}, {}, {1 / s, 1 / s, 0}};
  const AssembledFace left{{p(0, 1.5, 0), p(0, 0, 0), p(0, 0, 4), p(0, 1.5, 4)}, {}, {-1, 0, 0}};
  const std::vector<AssembledFace> faces{bottom, top, front, diagonal, left};
  const Result<Brep> prism = assemble_solid(faces, 1e-9);
  REQUIRE(prism.has_value());

  ShellOptions options;
  options.thickness = 0.25;
  const std::array<FaceId, 0> closed{};
  const Result<ShellResult> result = shell_solid(prism.value(), closed, options);
  REQUIRE_FALSE(result.has_value());
  CHECK(result.error().code == Errc::unsupported);
  CHECK(std::string_view(result.error().message) ==
        "shell_solid: only axis-aligned faces are supported");
}

TEST_CASE("shell half-space classifier sees a bounded interior") {
  const std::vector<Plane3> box{
      {{-1, 0, 0}, -0.25}, {{1, 0, 0}, 1.75}, {{0, -1, 0}, -0.25},
      {{0, 1, 0}, 1.75},   {{0, 0, -1}, -0.25}, {{0, 0, 1}, 1.75}};
  CHECK(detail::shell_half_spaces_have_interior(box, 1e-9));

  std::vector<Plane3> thin = box;
  thin[1].offset = 0.250000002;  // 2e-9 wide, so the half-width equals the tolerance
  CHECK_FALSE(detail::shell_half_spaces_have_interior(thin, 1e-9));
}

TEST_CASE("shell refuses any oblique face in any face order") {
  const double a = 30.0 * 3.14159265358979323846 / 180.0;
  const double b = 20.0 * 3.14159265358979323846 / 180.0;
  const auto rot = [a, b](const Vec3d& p) {
    const double x = std::cos(a) * p[0] - std::sin(a) * p[1];
    const double y = std::sin(a) * p[0] + std::cos(a) * p[1];
    const double y2 = std::cos(b) * y - std::sin(b) * p[2];
    const double z2 = std::sin(b) * y + std::cos(b) * p[2];
    return Vec3d{x, y2, z2};
  };
  const Point3d center{1.0, 1.0, 1.0};
  std::array<Point3d, 8> corners{};
  std::size_t c = 0;
  for (int sx = -1; sx <= 1; sx += 2) {
    for (int sy = -1; sy <= 1; sy += 2) {
      for (int sz = -1; sz <= 1; sz += 2) {
        const Vec3d offset = rot(Vec3d{static_cast<double>(sx), static_cast<double>(sy),
                                       static_cast<double>(sz)});
        corners[c++] = Point3d{center[0] + offset[0], center[1] + offset[1], center[2] + offset[2]};
      }
    }
  }
  std::vector<AssembledFace> faces;
  for (int axis = 0; axis < 3; ++axis) {
    for (int sign = -1; sign <= 1; sign += 2) {
      Vec3d n{0, 0, 0};
      n[axis] = static_cast<double>(sign);
      const Vec3d normal = rot(n);
      Vec3d basis{0, 0, 0};
      basis[(axis + 1) % 3] = 1.0;
      const Vec3d u = rot(basis);
      const Vec3d v = cross(normal, u);
      std::vector<Point3d> points;
      for (const Point3d& p : corners) {
        const Vec3d d = p - center;
        if (dot(d, normal) > 0.0) points.push_back(p);
      }
      REQUIRE(points.size() == 4u);
      std::sort(points.begin(), points.end(), [&](const Point3d& l, const Point3d& r) {
        const Vec3d dl = l - center;
        const Vec3d dr = r - center;
        return std::atan2(dot(dl, v), dot(dl, u)) < std::atan2(dot(dr, v), dot(dr, u));
      });
      faces.push_back(AssembledFace{points, {}, normal});
    }
  }
  ShellOptions options;
  options.thickness = 0.25;
  const auto check = [&options](std::vector<AssembledFace> input, const char* order) {
    INFO("face order: " << order);
    const Result<Brep> solid = assemble_solid(input, 1e-9);
    REQUIRE(solid.has_value());
    REQUIRE(validate_topology(solid.value()).ok());
    const Result<ShellResult> result = shell_solid(solid.value(), {}, options);
    REQUIRE_FALSE(result.has_value());
    CHECK(result.error().code == Errc::unsupported);
    CHECK(std::string_view(result.error().message) ==
          "shell_solid: only axis-aligned faces are supported");
  };
  check(faces, "canonical");
  std::reverse(faces.begin(), faces.end());
  check(faces, "reversed");
}


TEST_CASE("shell rejects bad options and dead faces") {
  const Brep box = box_solid(Point3d{0.0, 0.0, 0.0}, Point3d{2.0, 2.0, 2.0});
  const std::array<FaceId, 0> open{};
  ShellOptions options;
  options.thickness = 0.0;
  CHECK(shell_solid(box, open, options).error().code == Errc::degenerate);

  options.thickness = 0.25;
  options.tolerance = std::numeric_limits<double>::infinity();
  CHECK(shell_solid(box, open, options).error().code == Errc::non_finite);

  const std::array<FaceId, 1> dead{FaceId{}};
  options.tolerance = 1e-9;
  CHECK(shell_solid(box, dead, options).error().code == Errc::index_out_of_range);
}

TEST_CASE("a duplicate open face request is still refused") {
  const Brep box = box_solid(Point3d{0.0, 0.0, 0.0}, Point3d{2.0, 3.0, 4.0});
  FaceId top{};
  for (const FaceId face : box.face_ids()) {
    const auto& plane = std::get<Plane3>(box.geometry().surfaces[box.face(face).surface]);
    const Vec3d normal = box.face(face).reversed ? plane.normal * -1.0 : plane.normal;
    if (normal[2] > 0.5) top = face;
  }
  REQUIRE_FALSE(top.is_null());
  ShellOptions options;
  options.thickness = 0.25;
  const std::array<FaceId, 2> open{top, top};
  const Result<ShellResult> result = shell_solid(box, open, options);
  REQUIRE_FALSE(result.has_value());
  CHECK(result.error().code == Errc::unsupported);
  CHECK(std::string_view(result.error().message) ==
        "shell_solid: open faces are not supported in M4c-1");
}

TEST_CASE("shell keeps a representable seed below one coordinate ulp") {
  const Brep box = box_solid(Point3d{0.0, 0.0, 0.0}, Point3d{2.0, 3.0, 4.0});
  ShellOptions options;
  options.thickness = 1e-300;
  options.tolerance = 1e-9;
  const std::array<FaceId, 0> open{};
  const Result<ShellResult> result = shell_solid(box, open, options);
  REQUIRE_FALSE(result.has_value());
  CHECK(result.error().code == Errc::degenerate);
  CHECK(std::string_view(result.error().message) ==
        "shell_solid: wall thickness leaves no measurable interior at this tolerance");
}
