// Task 6 bench: `boolean_solid_curved` on the M3c-3a sphere cap corpus.
//
// Three cases over one unit sphere and axis-aligned boxes:
//   - intersection with box z >= 0.5 (spherical cap, h = 0.5),
//   - difference with the same box (complement of the cap),
//   - intersection with the slab -0.8 <= z <= 0.8 (band between two caps).
//
// Each case runs one untimed warm-up call, then `kRuns` timed calls. Two
// intervals are measured per iteration: the whole `boolean_solid_curved` call
// and the `detail::curved_solid_volume` call on its result. Every iteration
// verifies the volume against the analytic value (1e-6 relative) and against
// the first iteration bit for bit; the first and last results are also checked
// with `validate_topology` and for closed shells. A 64-bit FNV-1a checksum over
// every measured volume bit pattern is printed so a mismatching repeat run is
// visible in the output as well as in the exit code. The benchmark is
// informational: there is no time gate.
//
// Exit codes: 0 all volumes/topologies verified, 2 setup, driver or
// verification error (volume, topology or checksum mismatch). Progress goes to
// stderr; the per-case numbers and the machine/compiler line go to stdout.
#include <geom/brep/boolean_solid_curved.hpp>
#include <geom/brep/detail/curved_volume.hpp>
#include <geom/brep/primitives.hpp>
#include <geom/brep/validate.hpp>

#include <algorithm>
#include <array>
#include <bit>
#include <chrono>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstdio>
#include <vector>

using namespace geom::v1;

namespace {

constexpr double kTolerance = 1e-9;
constexpr double kVolumeRelative = 1e-6;
constexpr int kRuns = 100;
// pi * h^2 * (3R - h) / 3 for R = 1: h = 0.5, h = 0.2 (twice for the band).
constexpr double kCapVolume = 0.6544984694978736;
constexpr double kComplementVolume = 3.534291735288517;
constexpr double kSlabVolume = 3.9542179533184135;
constexpr std::uint64_t kFnvOffset = 1469598103934665603ull;
constexpr std::uint64_t kFnvPrime = 1099511628211ull;

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

std::uint64_t checksum_volume(std::uint64_t checksum, double volume) {
  const std::uint64_t bits = std::bit_cast<std::uint64_t>(volume);
  for (int shift = 0; shift < 64; shift += 8) {
    checksum = (checksum ^ ((bits >> shift) & 0xffu)) * kFnvPrime;
  }
  return checksum;
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
  SolidBooleanOp op;
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
  const Result<Brep> sphere_result =
      make_sphere(Point3d{0.0, 0.0, 0.0}, 1.0, kTolerance);
  if (!sphere_result.has_value()) {
    print_error(sphere_result.error());
    return 2;
  }
  const Brep& sphere = sphere_result.value();
  const std::vector<SolidId> sphere_solids = sphere.solid_ids();
  if (sphere_solids.size() != 1) {
    std::printf("error=sphere must hold exactly one solid\n");
    return 2;
  }

  const std::array<Case, 3> cases{{
      {"cap", {-2.0, -2.0, 0.5}, {2.0, 2.0, 2.0}, SolidBooleanOp::intersection, kCapVolume},
      {"complement",
       {-2.0, -2.0, 0.5},
       {2.0, 2.0, 2.0},
       SolidBooleanOp::difference,
       kComplementVolume},
      {"slab", {-2.0, -2.0, -0.8}, {2.0, 2.0, 0.8}, SolidBooleanOp::intersection, kSlabVolume},
  }};

  std::fprintf(stderr, "[bench] m3c3 sphere R=1 tol=%g runs=%d\n", kTolerance, kRuns);
  std::fflush(stderr);

  bool verified = true;
  std::uint64_t corpus_checksum = kFnvOffset;
  for (const Case& entry : cases) {
    const Result<Brep> box_result = make_box(entry.lower, entry.upper, kTolerance);
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
    const SolidId sphere_solid = sphere_solids.front();
    const SolidId box_solid = box_solids.front();

    // One untimed call warms the allocator and the analytic caches; its result
    // is verified but not part of the measured sample count.
    const Result<CurvedBooleanResult> warmup = boolean_solid_curved(
        sphere, sphere_solid, box, box_solid, entry.op, CurvedBooleanOptions{kTolerance});
    if (!warmup.has_value()) {
      print_error(warmup.error());
      return 2;
    }
    if (!verify_topology(entry.name, warmup.value().brep)) return 2;

    std::vector<double> boolean_samples;
    std::vector<double> volume_samples;
    boolean_samples.reserve(kRuns);
    volume_samples.reserve(kRuns);
    std::size_t faces = 0;
    double measured_volume = 0.0;
    double max_cut_residual = 0.0;
    std::size_t cuts = 0;
    std::size_t kept = 0;
    std::uint64_t baseline_bits = 0;
    std::size_t volume_mismatches = 0;
    std::uint64_t checksum = kFnvOffset;
    for (int run = 0; run < kRuns; ++run) {
      const auto boolean_start = std::chrono::steady_clock::now();
      const Result<CurvedBooleanResult> result = boolean_solid_curved(
          sphere, sphere_solid, box, box_solid, entry.op, CurvedBooleanOptions{kTolerance});
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
      const double value = volume.value();
      checksum = checksum_volume(checksum, value);
      corpus_checksum = checksum_volume(corpus_checksum, value);
      if (run == 0) {
        baseline_bits = std::bit_cast<std::uint64_t>(value);
      } else if (std::bit_cast<std::uint64_t>(value) != baseline_bits) {
        ++volume_mismatches;
      }
      const double error = std::abs(value - entry.expected);
      if (!(error <= kVolumeRelative * std::abs(entry.expected))) {
        std::printf("error=%s: volume %.12f expected %.12f residual %.3e\n", entry.name, value,
                    entry.expected, error);
        verified = false;
      }
      if (run == 0 || run == kRuns - 1) {
        if (!verify_topology(entry.name, brep)) verified = false;
      }
      faces = 0;
      for (const ShellId shell : brep.solid_shells(solids.front())) {
        faces += brep.shell_faces(shell).size();
      }
      measured_volume = value;
      max_cut_residual = result.value().report.max_cut_residual;
      cuts = result.value().report.cuts;
      kept = result.value().report.kept;
    }
    if (volume_mismatches != 0u) {
      std::printf("error=%s: %zu of %d volumes are not bit-identical\n", entry.name,
                  volume_mismatches, kRuns);
      verified = false;
    }

    std::printf("m3c3-sphere-cap case=%s runs=%d ", entry.name, kRuns);
    print_times("boolean", boolean_samples);
    std::printf(" ");
    print_times("volume", volume_samples);
    std::printf(
        " faces=%zu expected_volume=%.15f measured_volume=%.15f cuts=%zu kept=%zu "
        "max_cut_residual=%.6g checksum=0x%016llx verified=%s\n",
        faces, entry.expected, measured_volume, cuts, kept, max_cut_residual,
        static_cast<unsigned long long>(checksum), verified ? "yes" : "no");
    std::fflush(stdout);
  }

  std::printf(
      "m3c3-sphere-cap compiler=%s os=%s pointer_bits=%zu runs=%d checksum=0x%016llx "
      "verified=%s\n",
      compiler_name(), os_name(), sizeof(void*) * 8u, kRuns,
      static_cast<unsigned long long>(corpus_checksum), verified ? "yes" : "no");
  std::fflush(stdout);

  return verified ? 0 : 2;
}
