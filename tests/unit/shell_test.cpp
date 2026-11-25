#include <doctest/doctest.h>

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

TEST_CASE("open-top shell opens one pocket") {
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
  REQUIRE(result.has_value());
  CHECK(result.value().cavity_shells == 0u);
  CHECK(result.value().brep.solid_count() == 1u);
  CHECK(result.value().brep.solid_shells(result.value().brep.solid_ids().front()).size() == 1u);
  CHECK(validate_topology(result.value().brep).ok());
  // 24 - (1.5 * 2.5 * 3.75).
  CHECK(std::abs(volume_of(result.value().brep) - 9.9375) <= 1e-9);
}

TEST_CASE("opposite open faces make a shell through pocket") {
  const Brep box = box_solid(Point3d{0.0, 0.0, 0.0}, Point3d{2.0, 3.0, 4.0});
  std::array<FaceId, 2> open{};
  std::size_t count = 0;
  for (const FaceId face : box.face_ids()) {
    const auto& plane = std::get<Plane3>(box.geometry().surfaces[box.face(face).surface]);
    const Vec3d normal = box.face(face).reversed ? plane.normal * -1.0 : plane.normal;
    if (std::abs(normal[2]) > 0.5) open[count++] = face;
  }
  REQUIRE(count == 2u);
  ShellOptions options;
  options.thickness = 0.25;
  const Result<ShellResult> result = shell_solid(box, open, options);
  REQUIRE(result.has_value());
  CHECK(result.value().cavity_shells == 0u);
  CHECK(result.value().brep.solid_shells(result.value().brep.solid_ids().front()).size() == 1u);
  CHECK(std::abs(volume_of(result.value().brep) - 9.0) <= 1e-9);
}

TEST_CASE("shell rejects a wall below the boolean resolution") {
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
  CHECK(std::string_view(result.error().message).find("at least one face must stay closed") !=
        std::string::npos);
}

TEST_CASE("shell cuts axis-aligned faces before an oblique face") {
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
  const Result<ShellResult> shelled = shell_solid(prism.value(), closed, options);
  REQUIRE(shelled.has_value());
  CHECK(shelled.value().cavity_shells == 1u);
  CHECK(std::abs(volume_of(shelled.value().brep) - 3.76868686707646) <= 1e-9);

  FaceId top_face{};
  for (const FaceId face : prism.value().face_ids()) {
    const auto& plane =
        std::get<Plane3>(prism.value().geometry().surfaces[prism.value().face(face).surface]);
    const double nz = prism.value().face(face).reversed ? -plane.normal[2] : plane.normal[2];
    if (nz > 0.5) top_face = face;
  }
  REQUIRE_FALSE(top_face.is_null());
  const std::array<FaceId, 1> open{top_face};
  const Result<ShellResult> pocket = shell_solid(prism.value(), open, options);
  REQUIRE(pocket.has_value());
  CHECK(pocket.value().cavity_shells == 0u);
  CHECK(std::abs(volume_of(pocket.value().brep) - 3.71645021472478) <= 1e-9);
}
