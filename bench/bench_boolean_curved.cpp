// Task 9 bench: `boolean_solid_curved` throughput on an analytic cylinder/box
// pair.
//
// The operands are a unit cylinder (r = 1, h = 2, bottom cap at z = 0) and the
// offset box [0.5, 2] x [-0.5, 0.5] x [0.5, 1.5]. The three operations have
// closed-form volumes: with I = pi/6 + sqrt(3)/4 - 1/2 (the disk cap intersection
// area, the box is one unit tall),
//
//   intersection = I
//   union        = 1.5 + 2*pi - I   (box volume 1.5, cylinder volume 2*pi)
//   difference   = 2*pi - I
//
// Each operation is built once through a warm-up call, then timed over five
// runs; every run is verified with `validate_topology` and
// `detail::curved_solid_volume` against the closed form at 1e-9 relative. The
// timed region is one whole `boolean_solid_curved` call, i.e. both classification
// meshes, the surface-surface intersections, the cuts and the assembly.
//
// The classification mesh is built with the default `TessellationOptions`
// (deflection 1e-3). That is about 71 chord segments per full unit circle, which
// is recorded below instead of the plan's original "5000 lateral samples": the
// driver exposes no deflection knob, so claiming a finer sampling would be
// false. The target is the median: 2 seconds per operation.
//
// Exit codes: 0 all volumes/topologies verified (the time target may still be
// missed, which prints FAIL), 1 a result failed validation or the volume check,
// 2 setup or driver error. Progress goes to stderr; per-operation numbers and the
// machine/compiler line go to stdout.
#include <geom/brep/boolean_solid_curved.hpp>
#include <geom/brep/detail/curved_volume.hpp>
#include <geom/brep/primitives.hpp>
#include <geom/brep/tessellate.hpp>
#include <geom/brep/validate.hpp>

#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstdio>
#include <numbers>
#include <vector>

using namespace geom::v1;

namespace {

constexpr double kPi = std::numbers::pi_v<double>;
constexpr double kTolerance = 1e-9;
constexpr double kTargetMs = 2000.0;
constexpr int kRuns = 5;

double ms_since(std::chrono::steady_clock::time_point start) {
  return std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - start)
      .count();
}

void print_error(const Error& error) {
  std::printf("error=%s: %s\n", to_string(error.code), error.message);
  std::fflush(stdout);
}

std::size_t face_count(const Brep& brep, SolidId solid) {
  std::size_t faces = 0;
  for (const ShellId shell : brep.solid_shells(solid)) {
    faces += brep.shell_faces(shell).size();
  }
  return faces;
}

const char* compiler_name() {
#if defined(_MSC_VER)
  return "msvc";
#elif defined(__clang__)
  return "clang";
#elif defined(__GNUC__)
  return "gcc";
#else
  return "unknown";
#endif
}

const char* os_name() {
#if defined(_WIN32)
  return "windows";
#elif defined(__linux__)
  return "linux";
#elif defined(__APPLE__)
  return "macos";
#else
  return "unknown";
#endif
}

struct OpCase {
  const char* name;
  SolidBooleanOp op;
  double expected;
};

}  // namespace

int main() {
  const Result<Brep> cylinder =
      make_cylinder(Point3d{0.0, 0.0, 0.0}, Vec3d{0.0, 0.0, 1.0}, 1.0, 2.0, kTolerance);
  if (!cylinder.has_value()) {
    print_error(cylinder.error());
    return 2;
  }
  const Result<Brep> box = make_box(Point3d{0.5, -0.5, 0.5}, Point3d{2.0, 0.5, 1.5}, kTolerance);
  if (!box.has_value()) {
    print_error(box.error());
    return 2;
  }

  const Brep& a = cylinder.value();
  const Brep& b = box.value();
  const std::vector<SolidId> a_solids = a.solid_ids();
  const std::vector<SolidId> b_solids = b.solid_ids();
  if (a_solids.size() != 1 || b_solids.size() != 1) {
    std::printf("error=operands must each hold exactly one solid\n");
    return 2;
  }
  const SolidId a_solid = a_solids.front();
  const SolidId b_solid = b_solids.front();

  const TessellationOptions tessellation_options;
  const CircleCurve3 unit_circle{Point3d{0.0, 0.0, 0.0}, Vec3d{1.0, 0.0, 0.0},
                                 Vec3d{0.0, 1.0, 0.0}, 1.0};
  const Result<std::uint32_t> circle_segments = detail::tessellate_circle_segments(
      unit_circle, 0.0, 2.0 * kPi, tessellation_options.deflection, true);
  if (!circle_segments.has_value()) {
    print_error(circle_segments.error());
    return 2;
  }
  const Result<TriMesh> mesh = tessellate(a, tessellation_options);
  if (!mesh.has_value()) {
    print_error(mesh.error());
    return 2;
  }

  const double offset_segment = kPi / 6.0 + std::sqrt(3.0) / 4.0 - 0.5;
  const std::array<OpCase, 3> cases{{
      {"union", SolidBooleanOp::union_, 1.5 + 2.0 * kPi - offset_segment},
      {"intersection", SolidBooleanOp::intersection, offset_segment},
      {"difference", SolidBooleanOp::difference, 2.0 * kPi - offset_segment},
  }};

  std::fprintf(stderr, "[bench] curved boolean cylinder r=1 h=2 vs box offset tol=%g runs=%d\n",
               kTolerance, kRuns);
  std::fflush(stderr);

  bool verified = true;
  bool target_met = true;
  std::size_t last_faces = 0;
  for (const OpCase& entry : cases) {
    // One untimed call warms the allocator and the tessellation caches; its
    // result is discarded and it is not part of the measured sample count.
    const Result<CurvedBooleanResult> warmup = boolean_solid_curved(
        a, a_solid, b, b_solid, entry.op, CurvedBooleanOptions{kTolerance});
    if (!warmup.has_value()) {
      print_error(warmup.error());
      return 2;
    }

    std::vector<double> samples;
    samples.reserve(kRuns);
    for (int run = 0; run < kRuns; ++run) {
      const auto start = std::chrono::steady_clock::now();
      const Result<CurvedBooleanResult> result = boolean_solid_curved(
          a, a_solid, b, b_solid, entry.op, CurvedBooleanOptions{kTolerance});
      samples.push_back(ms_since(start));
      if (!result.has_value()) {
        print_error(result.error());
        return 2;
      }
      const Brep& brep = result.value().brep;
      if (!validate_topology(brep).ok()) {
        std::printf("error=%s: result failed validation\n", entry.name);
        verified = false;
        continue;
      }
      const std::vector<SolidId> solids = brep.solid_ids();
      if (solids.size() != 1) {
        std::printf("error=%s: result is not a single solid\n", entry.name);
        verified = false;
        continue;
      }
      const Result<double> volume = detail::curved_solid_volume(brep, solids.front(), kTolerance);
      if (!volume.has_value()) {
        print_error(volume.error());
        verified = false;
        continue;
      }
      const double error = std::abs(volume.value() - entry.expected);
      if (!(error <= 1e-9 * std::abs(entry.expected))) {
        std::printf("error=%s: volume %.12f expected %.12f residual %.3e\n", entry.name,
                    volume.value(), entry.expected, error);
        verified = false;
      }
      last_faces = face_count(brep, solids.front());
    }

    std::vector<double> sorted = samples;
    std::sort(sorted.begin(), sorted.end());
    const double min_ms = sorted.front();
    const double max_ms = sorted.back();
    const double median_ms = sorted[sorted.size() / 2];
    const bool passed = median_ms <= kTargetMs;
    target_met = target_met && passed;
    std::printf(
        "m3b2-curved-boolean op=%s runs=%d min_ms=%.3f median_ms=%.3f max_ms=%.3f "
        "target_ms=%.0f faces=%zu %s\n",
        entry.name, kRuns, min_ms, median_ms, max_ms, kTargetMs, last_faces,
        passed ? "PASS" : "FAIL");
    std::fflush(stdout);
  }

  std::printf(
      "m3b2-curved-boolean circle_segments=%u deflection=%.0e tessellation_vertices=%zu "
      "tessellation_triangles=%zu\n",
      circle_segments.value(), tessellation_options.deflection, mesh.value().positions.size(),
      mesh.value().triangles.size());
  std::printf("m3b2-curved-boolean compiler=%s os=%s pointer_bits=%zu target_met=%s verified=%s\n",
              compiler_name(), os_name(), sizeof(void*) * 8u, target_met ? "yes" : "no",
              verified ? "yes" : "no");
  std::fflush(stdout);

  return verified ? 0 : 1;
}
