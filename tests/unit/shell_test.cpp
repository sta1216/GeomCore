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

FaceId top_face(const Brep& brep) {
  for (const FaceId face : brep.face_ids()) {
    const auto& plane =
        std::get<Plane3>(brep.geometry().surfaces[brep.face(face).surface]);
    const Vec3d normal = brep.face(face).reversed ? plane.normal * -1.0 : plane.normal;
    if (normal[2] > 0.5) return face;
  }
  FAIL("box has no +Z face");
  return FaceId{};
}

FaceId bottom_face(const Brep& brep) {
  for (const FaceId face : brep.face_ids()) {
    const auto& plane =
        std::get<Plane3>(brep.geometry().surfaces[brep.face(face).surface]);
    const Vec3d normal = brep.face(face).reversed ? plane.normal * -1.0 : plane.normal;
    if (normal[2] < -0.5) return face;
  }
  FAIL("box has no -Z face");
  return FaceId{};
}

FaceId right_face(const Brep& brep) {
  for (const FaceId face : brep.face_ids()) {
    const auto& plane =
        std::get<Plane3>(brep.geometry().surfaces[brep.face(face).surface]);
    const Vec3d normal = brep.face(face).reversed ? plane.normal * -1.0 : plane.normal;
    if (normal[0] > 0.5) return face;
  }
  FAIL("box has no +X face");
  return FaceId{};
}

Brep right_tetrahedron(double a) {
  const auto p = [](double x, double y, double z) { return Point3d{x, y, z}; };
  const std::array<AssembledFace, 4> faces{{
      {{p(0, 0, 0), p(0, 0, a), p(0, a, 0)}, {}, {-1, 0, 0}},
      {{p(0, 0, 0), p(a, 0, 0), p(0, 0, a)}, {}, {0, -1, 0}},
      {{p(0, a, 0), p(a, 0, 0), p(0, 0, 0)}, {}, {0, 0, -1}},
      {{p(a, 0, 0), p(0, a, 0), p(0, 0, a)}, {}, {1 / std::sqrt(3.0), 1 / std::sqrt(3.0),
                                                  1 / std::sqrt(3.0)}},
  }};
  const Result<Brep> tetra = assemble_solid(faces, 1e-9);
  REQUIRE(tetra.has_value());
  return tetra.value();
}

FaceId slant_face(const Brep& brep) {
  for (const FaceId face : brep.face_ids()) {
    const auto& plane =
        std::get<Plane3>(brep.geometry().surfaces[brep.face(face).surface]);
    const Vec3d normal = brep.face(face).reversed ? plane.normal * -1.0 : plane.normal;
    if (normal[0] > 0.5 && normal[1] > 0.5 && normal[2] > 0.5) return face;
  }
  FAIL("tetrahedron has no slant face");
  return FaceId{};
}

// A 2 x 3 x 4 box rotated about +Z, assembled from transformed corners with
// transformed outward normals. At a coarse tolerance this is exactly the input
// that an axis-aligned test based on tolerance would misclassify as a box.
Brep rotated_box(double angle) {
  const double c = std::cos(angle);
  const double s = std::sin(angle);
  const auto p = [c, s](double x, double y, double z) {
    return Point3d{c * x - s * y, s * x + c * y, z};
  };
  const auto n = [c, s](double x, double y, double z) {
    return Vec3d{c * x - s * y, s * x + c * y, z};
  };
  const std::array<AssembledFace, 6> faces{{
      {{p(0, 0, 0), p(0, 0, 4), p(0, 3, 4), p(0, 3, 0)}, {}, n(-1, 0, 0)},
      {{p(2, 0, 0), p(2, 3, 0), p(2, 3, 4), p(2, 0, 4)}, {}, n(1, 0, 0)},
      {{p(0, 0, 0), p(2, 0, 0), p(2, 0, 4), p(0, 0, 4)}, {}, n(0, -1, 0)},
      {{p(0, 3, 0), p(0, 3, 4), p(2, 3, 4), p(2, 3, 0)}, {}, n(0, 1, 0)},
      {{p(0, 0, 0), p(0, 3, 0), p(2, 3, 0), p(2, 0, 0)}, {}, n(0, 0, -1)},
      {{p(0, 0, 4), p(2, 0, 4), p(2, 3, 4), p(0, 3, 4)}, {}, n(0, 0, 1)},
  }};
  const Result<Brep> box = assemble_solid(faces, 1e-9);
  REQUIRE(box.has_value());
  return box.value();
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

TEST_CASE("shell of a rotated box keeps its wall at coarse tolerance") {
  const Brep box = rotated_box(0.04);
  const std::array<FaceId, 0> open{};
  ShellOptions options;
  options.thickness = 0.25;
  options.tolerance = 1e-3;
  const Result<ShellResult> result = shell_solid(box, open, options);
  REQUIRE(result.has_value());
  CHECK(result.value().cavity_shells == 1u);
  CHECK(result.value().brep.shell_count() == 2u);
  CHECK(validate_topology(result.value().brep).ok());
  // A rigid rotation does not change the wall volume; only a genuine
  // half-space cavity can recover it at this coarse tolerance.
  CHECK(std::abs(volume_of(result.value().brep) - 10.875) <= 1e-6);
}

TEST_CASE("shell rejects a wall at or beyond the half-thickness") {
  const Brep box = box_solid(Point3d{0.0, 0.0, 0.0}, Point3d{2.0, 3.0, 4.0});
  ShellOptions options;
  options.thickness = 1.0;
  options.tolerance = 1e-9;
  const std::array<FaceId, 0> open{};
  const Result<ShellResult> result = shell_solid(box, open, options);
  REQUIRE_FALSE(result.has_value());
  CHECK(result.error().code == Errc::degenerate);
  CHECK(std::string_view(result.error().message) ==
        "shell_solid: wall thickness leaves no interior");
}

TEST_CASE("shell rejects an oblique wall that leaves no interior") {
  const Brep tetra = right_tetrahedron(2.0);
  ShellOptions options;
  options.thickness = 0.7;
  options.tolerance = 1e-9;
  const std::array<FaceId, 0> open{};
  const Result<ShellResult> result = shell_solid(tetra, open, options);
  REQUIRE_FALSE(result.has_value());
  CHECK(result.error().code == Errc::degenerate);
  CHECK(std::string_view(result.error().message) ==
        "shell_solid: wall thickness leaves no measurable interior at this tolerance");
}

TEST_CASE("shell opens a box through one face") {
  const Brep box = box_solid(Point3d{0.0, 0.0, 0.0}, Point3d{2.0, 3.0, 4.0});
  const std::array<FaceId, 1> open{top_face(box)};
  ShellOptions options;
  options.thickness = 0.25;
  options.tolerance = 1e-9;
  const Result<ShellResult> result = shell_solid(box, open, options);
  REQUIRE(result.has_value());
  CHECK(result.value().cavity_shells == 0u);
  CHECK(result.value().brep.shell_count() == 1u);
  CHECK(validate_topology(result.value().brep).ok());
  const Result<MassProperties> mass = mass_properties(result.value().brep);
  REQUIRE(mass.has_value());
  // 24 - (1.5 * 2.5 * 3.75).
  CHECK(std::abs(mass.value().volume - 9.9375) <= 1e-9);
}

TEST_CASE("shell opens a box through two opposite faces") {
  const Brep box = box_solid(Point3d{0.0, 0.0, 0.0}, Point3d{2.0, 3.0, 4.0});
  const std::array<FaceId, 2> open{top_face(box), bottom_face(box)};
  ShellOptions options;
  options.thickness = 0.25;
  options.tolerance = 1e-9;
  const Result<ShellResult> result = shell_solid(box, open, options);
  REQUIRE(result.has_value());
  CHECK(result.value().cavity_shells == 0u);
  CHECK(result.value().brep.shell_count() == 1u);
  CHECK(validate_topology(result.value().brep).ok());
  // 24 - (1.5 * 2.5 * 4.0): the pocket runs through both open ends.
  CHECK(std::abs(volume_of(result.value().brep) - 9.0) <= 1e-9);
}

TEST_CASE("shell opens a box through one side face") {
  const Brep box = box_solid(Point3d{0.0, 0.0, 0.0}, Point3d{2.0, 3.0, 4.0});
  const std::array<FaceId, 1> open{right_face(box)};
  ShellOptions options;
  options.thickness = 0.25;
  options.tolerance = 1e-9;
  const Result<ShellResult> result = shell_solid(box, open, options);
  REQUIRE(result.has_value());
  CHECK(result.value().cavity_shells == 0u);
  CHECK(result.value().brep.shell_count() == 1u);
  CHECK(validate_topology(result.value().brep).ok());
  // 24 - (1.75 * 2.5 * 3.5): opening is on +X, so it must not depend on Z.
  CHECK(std::abs(volume_of(result.value().brep) - 8.6875) <= 1e-9);
}

TEST_CASE("shell of an oblique tetrahedron matches the inset volume") {
  const Brep tetra = right_tetrahedron(2.0);
  const std::array<FaceId, 0> open{};
  ShellOptions options;
  options.thickness = 0.1;
  options.tolerance = 1e-9;
  const Result<ShellResult> result = shell_solid(tetra, open, options);
  REQUIRE(result.has_value());
  CHECK(result.value().cavity_shells == 1u);
  CHECK(result.value().brep.shell_count() == 2u);
  CHECK(validate_topology(result.value().brep).ok());
  const double outer = 4.0 / 3.0;
  const double s = 2.0 - options.thickness * (3.0 + std::sqrt(3.0));
  const double inner = s * s * s / 6.0;
  CHECK(std::abs(volume_of(result.value().brep) - (outer - inner)) <= 1e-9);
}

TEST_CASE("shell rejects an open tetrahedron with no removable interior") {
  const Brep tetra = right_tetrahedron(2.0);
  const char* const no_interior =
      "shell_solid: wall thickness leaves no measurable interior at this tolerance";

  // The relaxed slant plane lies wholly outside the solid once the three legs
  // are inset by 0.7, so the pocket never enters the tetrahedron.
  ShellOptions options;
  options.thickness = 0.7;
  options.tolerance = 1e-9;
  const std::array<FaceId, 1> slant{slant_face(tetra)};
  const Result<ShellResult> slant_result = shell_solid(tetra, slant, options);
  REQUIRE_FALSE(slant_result.has_value());
  CHECK(slant_result.error().code == Errc::degenerate);
  CHECK(std::string_view(slant_result.error().message) == no_interior);

  // The open leg cavity is non-empty but sits entirely below z = 0, outside
  // the tetrahedron, so an unchanged result must be rejected rather than
  // returned as a successful open shell.
  options.thickness = 0.54;
  const std::array<FaceId, 1> leg{bottom_face(tetra)};
  const Result<ShellResult> leg_result = shell_solid(tetra, leg, options);
  REQUIRE_FALSE(leg_result.has_value());
  CHECK(leg_result.error().code == Errc::degenerate);
  CHECK(std::string_view(leg_result.error().message) == no_interior);
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
    if (!validate_topology(result.value().brep).ok()) {
      first_failure = "invalid output topology";
      break;
    }
    if (std::abs(volume_of(result.value().brep) - 10.875) > 1e-9) {
      first_failure = "wrong volume";
      break;
    }
    const std::array<FaceId, 1> open{top_face(solid.value())};
    const Result<ShellResult> opened = shell_solid(solid.value(), open, options);
    if (!opened.has_value()) {
      first_failure = opened.error().message;
      break;
    }
    if (opened.value().cavity_shells != 0u) {
      first_failure = "open shell retained a cavity shell";
      break;
    }
    if (opened.value().brep.shell_count() != 1u) {
      first_failure = "open shell has more than one shell";
      break;
    }
    if (!validate_topology(opened.value().brep).ok()) {
      first_failure = "invalid open output topology";
      break;
    }
    if (std::abs(volume_of(opened.value().brep) - 9.9375) > 1e-9) {
      first_failure = "wrong open volume";
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
  CHECK(std::string_view(result.error().message) ==
        "shell_solid: wall thickness leaves no measurable interior at this tolerance");
}

TEST_CASE("shell accepts a thin but removable open wall") {
  const Brep box = box_solid(Point3d{0.0, 0.0, 0.0}, Point3d{2.0, 3.0, 4.0});
  const std::array<FaceId, 1> open{top_face(box)};
  ShellOptions options;
  options.thickness = 1e-6;
  options.tolerance = 1e-9;
  const Result<ShellResult> result = shell_solid(box, open, options);
  REQUIRE(result.has_value());
  CHECK(result.value().cavity_shells == 0u);
  CHECK(result.value().brep.shell_count() == 1u);
  CHECK(validate_topology(result.value().brep).ok());
  // The inset walls are one micrometre thick: above the boolean resolution,
  // so the cavity must survive instead of being rejected as unmeasurable.
  const double expected = 24.0 - (2.0 - 2e-6) * (3.0 - 2e-6) * (4.0 - 1e-6);
  CHECK(std::abs(volume_of(result.value().brep) - expected) <= 1e-9);
}

TEST_CASE("shell keeps a small open shell above its own volume scale") {
  // The removed-material guard must be relative to the model: 2x3x4 scaled by
  // 1e-5 has total volume 2.4e-14, and the open pocket removes 1.40625e-14.
  // An absolute 1e-12 floor would reject this valid shell.
  const Brep box = box_solid(Point3d{0.0, 0.0, 0.0},
                             Point3d{2e-5, 3e-5, 4e-5});
  const std::array<FaceId, 1> open{top_face(box)};
  ShellOptions options;
  options.thickness = 2.5e-6;
  options.tolerance = 1e-14;
  const Result<ShellResult> result = shell_solid(box, open, options);
  REQUIRE(result.has_value());
  CHECK(result.value().cavity_shells == 0u);
  CHECK(result.value().brep.shell_count() == 1u);
  CHECK(validate_topology(result.value().brep).ok());
  const double expected = 9.9375e-15;
  CHECK(std::abs(volume_of(result.value().brep) - expected) <= 1e-9 * expected);
}

TEST_CASE("shell rejects opening every face") {
  const Brep box = box_solid(Point3d{0.0, 0.0, 0.0}, Point3d{2.0, 3.0, 4.0});
  const auto ids = box.face_ids();
  std::vector<FaceId> faces(ids.begin(), ids.end());
  faces.push_back(top_face(box));  // Repeating a face must not evade all-open.
  ShellOptions options;
  options.thickness = 0.25;
  const Result<ShellResult> result = shell_solid(box, faces, options);
  REQUIRE_FALSE(result.has_value());
  CHECK(result.error().code == Errc::unsupported);
  CHECK(std::string_view(result.error().message) ==
        "shell_solid: at least one face must stay closed");
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

TEST_CASE("shell accepts a duplicate open face once") {
  const Brep box = box_solid(Point3d{0.0, 0.0, 0.0}, Point3d{2.0, 3.0, 4.0});
  const FaceId top = top_face(box);
  const std::array<FaceId, 2> open{top, top};
  ShellOptions options;
  options.thickness = 0.25;
  options.tolerance = 1e-9;
  const Result<ShellResult> result = shell_solid(box, open, options);
  REQUIRE(result.has_value());
  CHECK(result.value().cavity_shells == 0u);
  CHECK(result.value().brep.shell_count() == 1u);
  CHECK(validate_topology(result.value().brep).ok());
  CHECK(std::abs(volume_of(result.value().brep) - 9.9375) <= 1e-9);
}

TEST_CASE("shell rejects a wall below one coordinate ulp") {
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
