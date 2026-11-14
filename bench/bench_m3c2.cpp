// Task 8 bench: `boolean_solid_curved` on the M3c-2 torus corpus.
//
// Three cases over one torus (R = 1, r = 0.5, axis +z):
//   - intersection with the half space z >= 0 (axis-perpendicular cut at the
//     seam-only band, v = pi),
//   - intersection with the half space x >= 0 (plane through the axis, cut
//     crosses the v seam),
//   - intersection with a containing box (no cut).
//
// Each case runs one untimed warm-up call, then `kRuns` timed calls. Two
// intervals are measured per iteration: the whole `boolean_solid_curved` call
// and the `detail::curved_solid_volume` call on its result, because the torus
// results are trimmed analytic faces whose volume integration dominates the
// wall clock. Every iteration verifies the volume against the analytic value
// (1e-6 relative); the first and last results are also checked with
// `validate_topology` and for closed shells. The benchmark is informational:
// there is no time gate, so it only prints PASS/FAIL for verification.
//
// Exit codes: 0 all volumes/topologies verified, 1 a result failed validation
// or the volume check, 2 setup or driver error. Progress goes to stderr; the
// per-case numbers and the machine/compiler line go to stdout.
#include <geom/brep/boolean_solid_curved.hpp>
#include <geom/brep/detail/curved_volume.hpp>
#include <geom/brep/primitives.hpp>
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

constexpr double kTolerance = 1e-9;
constexpr double kVolumeRelative = 1e-6;
constexpr int kRuns = 100;
// (2 * pi^2 * R * r^2) / 2 and 2 * pi^2 * R * r^2 for R = 1, r = 0.5.
constexpr double kHalfVolume = 2.4674011002723395;
constexpr double kWholeVolume = 4.934802200544679;

double ms_since(std::chrono::steady_clock::time_point start) {
  return std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - start)
      .count();
}

void print_error(const Error& error) {
  std::printf("error=%s: %s\n", to_string(error.code), error.message);
  std::fflush(stdout);
}

double percentile(const std::vector<double>& sorted, double fraction) {
  if (sorted.empty()) return 0.0;
  const double position = fraction * static_cast<double>(sorted.size());
  std::size_t index = static_cast<std::size_t>(std::ceil(position));
  if (index == 0u) index = 1u;
  if (index > sorted.size()) index = sorted.size();
  return sorted[index - 1u];
}

void print_times(const char* label, std::vector<double> samples) {
  std::sort(samples.begin(), samples.end());
  std::printf("%s_min_ms=%.3f %s_median_ms=%.3f %s_p95_ms=%.3f %s_max_ms=%.3f", label,
              samples.front(), label, percentile(samples, 0.5), label, percentile(samples, 0.95),
              label, samples.back());
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

struct Case {
  const char* name;
  Point3d lower;
  Point3d upper;
  double expected;
};

bool verify_topology(const char* name, const Brep& brep) {
  if (!validate_topology(brep).ok()) {
    std::printf("error=%s: result failed validation\n", name);
    return false;
  }
  const std::vector<SolidId> solids = brep.solid_ids();
  if (solids.size() != 1) {
    std::printf("error=%s: result is not a single solid\n", name);
    return false;
  }
  for (const ShellId shell : brep.solid_shells(solids.front())) {
    if (!brep.shell(shell).closed) {
      std::printf("error=%s: result has an open shell\n", name);
      return false;
    }
  }
  return true;
}

}  // namespace

int main() {
  const Result<Brep> torus_result =
      make_torus(Point3d{0.0, 0.0, 0.0}, Vec3d{0.0, 0.0, 1.0}, 1.0, 0.5, kTolerance);
  if (!torus_result.has_value()) {
    print_error(torus_result.error());
    return 2;
  }
  const Brep& torus = torus_result.value();
  const std::vector<SolidId> torus_solids = torus.solid_ids();
  if (torus_solids.size() != 1) {
    std::printf("error=torus must hold exactly one solid\n");
    return 2;
  }

  const std::array<Case, 3> cases{{
      {"torus-z-half", {-2.0, -2.0, 0.0}, {2.0, 2.0, 2.0}, kHalfVolume},
      {"torus-x-half", {0.0, -2.0, -2.0}, {2.0, 2.0, 2.0}, kHalfVolume},
      {"torus-in-box", {-2.0, -2.0, -2.0}, {2.0, 2.0, 2.0}, kWholeVolume},
  }};

  std::fprintf(stderr, "[bench] m3c2 torus R=1 r=0.5 tol=%g runs=%d\n", kTolerance, kRuns);
  std::fflush(stderr);

  bool verified = true;
  for (const Case& entry : cases) {
    const Result<Brep> box_result =
        make_box(entry.lower, entry.upper, kTolerance);
    if (!box_result.has_value()) {
      print_error(box_result.error());
      return 2;
    }
    const Brep& box = box_result.value();
    const std::vector<SolidId> box_solids = box.solid_ids();
    if (box_solids.size() != 1) {
      std::printf("error=%s: box must hold exactly one solid\n", entry.name);
      return 2;
    }
    const SolidId torus_solid = torus_solids.front();
    const SolidId box_solid = box_solids.front();

    // One untimed call warms the allocator and the analytic caches; its result
    // is verified but not part of the measured sample count.
    const Result<CurvedBooleanResult> warmup = boolean_solid_curved(
        torus, torus_solid, box, box_solid, SolidBooleanOp::intersection,
        CurvedBooleanOptions{kTolerance});
    if (!warmup.has_value()) {
      print_error(warmup.error());
      return 2;
    }
    if (!verify_topology(entry.name, warmup.value().brep)) return 1;

    std::vector<double> boolean_samples;
    std::vector<double> volume_samples;
    boolean_samples.reserve(kRuns);
    volume_samples.reserve(kRuns);
    std::size_t faces = 0;
    double measured_volume = 0.0;
    double max_cut_residual = 0.0;
    std::size_t cuts = 0;
    std::size_t kept = 0;
    for (int run = 0; run < kRuns; ++run) {
      const auto boolean_start = std::chrono::steady_clock::now();
      const Result<CurvedBooleanResult> result = boolean_solid_curved(
          torus, torus_solid, box, box_solid, SolidBooleanOp::intersection,
          CurvedBooleanOptions{kTolerance});
      boolean_samples.push_back(ms_since(boolean_start));
      if (!result.has_value()) {
        print_error(result.error());
        return 2;
      }
      const Brep& brep = result.value().brep;
      const std::vector<SolidId> solids = brep.solid_ids();
      if (solids.size() != 1) {
        std::printf("error=%s: result is not a single solid\n", entry.name);
        verified = false;
        continue;
      }
      const auto volume_start = std::chrono::steady_clock::now();
      const Result<double> volume =
          detail::curved_solid_volume(brep, solids.front(), kTolerance);
      volume_samples.push_back(ms_since(volume_start));
      if (!volume.has_value()) {
        print_error(volume.error());
        verified = false;
        continue;
      }
      const double error = std::abs(volume.value() - entry.expected);
      if (!(error <= kVolumeRelative * std::abs(entry.expected))) {
        std::printf("error=%s: volume %.12f expected %.12f residual %.3e\n", entry.name,
                    volume.value(), entry.expected, error);
        verified = false;
      }
      if (run == 0 || run == kRuns - 1) {
        if (!verify_topology(entry.name, brep)) verified = false;
      }
      faces = 0;
      for (const ShellId shell : brep.solid_shells(solids.front())) {
        faces += brep.shell_faces(shell).size();
      }
      measured_volume = volume.value();
      max_cut_residual = result.value().report.max_cut_residual;
      cuts = result.value().report.cuts;
      kept = result.value().report.kept;
    }

    std::printf("m3c2-curved-boolean case=%s runs=%d ", entry.name, kRuns);
    print_times("boolean", boolean_samples);
    std::printf(" ");
    print_times("volume", volume_samples);
    std::printf(
        " faces=%zu expected_volume=%.15f measured_volume=%.15f cuts=%zu kept=%zu "
        "max_cut_residual=%.6g verified=%s\n",
        faces, entry.expected, measured_volume, cuts, kept, max_cut_residual,
        verified ? "yes" : "no");
    std::fflush(stdout);
  }

  std::printf(
      "m3c2-curved-boolean compiler=%s os=%s pointer_bits=%zu runs=%d verified=%s\n",
      compiler_name(), os_name(), sizeof(void*) * 8u, kRuns, verified ? "yes" : "no");
  std::fflush(stdout);

  return verified ? 0 : 1;
}
