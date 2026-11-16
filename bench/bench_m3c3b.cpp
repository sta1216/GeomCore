// Task 5 bench: `boolean_solid_curved` on the M3c-3b sphere/sphere lens corpus.
//
// Three operations over two unit spheres whose centres are 1.5 apart:
//   - intersection: the lens, 11*pi/96,
//   - union: the two complements, 245*pi/96,
//   - difference: one complement minus the small cap, 39*pi/32.
//
// Each operation runs one untimed warm-up call, then `kRuns` timed calls. Two
// intervals are measured per iteration: the whole `boolean_solid_curved` call and
// the signed chart cap volume of its result. The chart sum is used instead of
// `detail::curved_solid_volume` because the lens rim at 1e-9 carries about 4098
// points and that integrator stalls on it in ear clipping; the chart sum is
// cross-checked once against `curved_solid_volume` on the M3c-3a cap fixture by
// QA-16 (`tests/unit/brep_boolean_curved_test.cpp`), not here. Every iteration
// verifies the chart volume against the closed form (1e-9 relative for the
// intersection, 1e-6 for the union and difference) and against the first
// iteration bit for bit; the first and last results are also checked with
// `validate_topology`, one solid, one shell, two faces, closed shells and
// `heal_detail::shell_is_closed_by_uses`. A 64-bit FNV-1a checksum over every
// measured chart volume bit pattern is printed so a mismatching repeat run is
// visible in the output as well as in the exit code. The benchmark is
// informational: there is no time gate.
//
// Exit codes: 0 all volumes/topologies verified, 2 setup, driver or verification
// error (volume, topology or checksum mismatch). Progress goes to stderr; the
// per-case numbers and the machine/compiler line go to stdout.
#include <geom/brep/boolean_solid_curved.hpp>
#include <geom/brep/face_uv.hpp>
#include <geom/brep/heal.hpp>
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
#include <numbers>
#include <optional>
#include <span>
#include <vector>

using namespace geom::v1;

namespace {

constexpr double kTolerance = 1e-9;
constexpr double kPi = std::numbers::pi_v<double>;
constexpr int kRuns = 20;
// Equal unit spheres 1.5 apart: caps of height 0.25.
constexpr double kLensVolume = 11.0 * kPi / 96.0;
constexpr double kUnionVolume = 245.0 * kPi / 96.0;
constexpr double kDifferenceVolume = 39.0 * kPi / 32.0;
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
  SolidBooleanOp op;
  double expected;
  double relative_tolerance;
};

// Signed chart volume of a solid built only from polar sphere caps. Mirrors the
// QA-16 helper so the benchmark shares its measurement without touching the
// stalled trim-loop integrator. A face without a pole at either v end is not a
// polar cap and returns nullopt.
std::optional<double> chart_volume(const Brep& brep) {
  double volume = 0.0;
  bool any_sphere = false;
  for (const FaceId face_id : brep.face_ids()) {
    const Face& face = brep.face(face_id);
    const auto* sphere = std::get_if<Sphere3>(&brep.geometry().surfaces[face.surface]);
    if (sphere == nullptr) continue;
    const Result<std::vector<FaceUvLoop>> loops = face_loops_uv(brep, face_id, kTolerance);
    if (!loops.has_value() || loops.value().size() != 1) return std::nullopt;
    const std::vector<Point2d>& ring = loops.value().front().points;
    if (ring.size() < 3) return std::nullopt;

    double twice_area = 0.0;
    double v_lo = 1e30;
    double v_hi = -1e30;
    for (std::size_t i = 0; i < ring.size(); ++i) {
      const Point2d& point = ring[i];
      const Point2d& next = ring[(i + 1) % ring.size()];
      twice_area += point[0] * next[1] - next[0] * point[1];
      v_lo = std::min(v_lo, point[1]);
      v_hi = std::max(v_hi, point[1]);
    }
    if (twice_area == 0.0) return std::nullopt;

    constexpr double kHalfPi = kPi / 2.0;
    double height = 0.0;
    if (std::abs(v_lo + kHalfPi) <= kTolerance) {
      height = sphere->radius * (1.0 + std::sin(v_hi));
    } else if (std::abs(v_hi - kHalfPi) <= kTolerance) {
      height = sphere->radius * (1.0 - std::sin(v_lo));
    } else {
      return std::nullopt;
    }
    const double sign = twice_area > 0.0 ? 1.0 : -1.0;
    const double radius = sphere->radius;
    volume += sign * kPi * height * height * (3.0 * radius - height) / 3.0;
    any_sphere = true;
  }
  if (!any_sphere) return std::nullopt;
  return volume;
}

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
  const std::span<const ShellId> shells = brep.solid_shells(solids.front());
  if (shells.size() != 1) {
    std::printf("error=%s: result is not a single shell\n", name);
    return false;
  }
  if (!brep.shell(shells.front()).closed) {
    std::printf("error=%s: result has an open shell\n", name);
    return false;
  }
  if (!geom::detail::heal_detail::shell_is_closed_by_uses(brep, shells.front())) {
    std::printf("error=%s: an edge is not used by exactly two coedges\n", name);
    return false;
  }
  std::size_t faces = 0;
  for (const ShellId shell : shells) {
    faces += brep.shell_faces(shell).size();
  }
  if (faces != 2) {
    std::printf("error=%s: result has %zu faces, expected 2\n", name, faces);
    return false;
  }
  return true;
}

}  // namespace

int main() {
  const Result<Brep> a_result = make_sphere(Point3d{0.0, 0.0, 0.0}, 1.0, kTolerance);
  const Result<Brep> b_result = make_sphere(Point3d{1.5, 0.0, 0.0}, 1.0, kTolerance);
  if (!a_result.has_value()) {
    print_error(a_result.error());
    return 2;
  }
  if (!b_result.has_value()) {
    print_error(b_result.error());
    return 2;
  }
  const Brep& a = a_result.value();
  const Brep& b = b_result.value();
  const std::vector<SolidId> a_solids = a.solid_ids();
  const std::vector<SolidId> b_solids = b.solid_ids();
  if (a_solids.size() != 1 || b_solids.size() != 1) {
    std::printf("error=operands must each hold exactly one solid\n");
    return 2;
  }

  const std::array<Case, 3> cases{{
      {"intersection", SolidBooleanOp::intersection, kLensVolume, 1e-9},
      {"union", SolidBooleanOp::union_, kUnionVolume, 1e-6},
      {"difference", SolidBooleanOp::difference, kDifferenceVolume, 1e-6},
  }};

  std::fprintf(stderr, "[bench] m3c3b sphere lens R=1 d=1.5 tol=%g runs=%d\n", kTolerance, kRuns);
  std::fflush(stderr);

  bool verified = true;
  std::uint64_t corpus_checksum = kFnvOffset;
  for (const Case& entry : cases) {
    // One untimed call warms the allocator and the analytic caches; its result
    // is verified but not part of the measured sample count.
    const Result<CurvedBooleanResult> warmup =
        boolean_solid_curved(a, a_solids.front(), b, b_solids.front(), entry.op,
                             CurvedBooleanOptions{kTolerance});
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
      const Result<CurvedBooleanResult> result =
          boolean_solid_curved(a, a_solids.front(), b, b_solids.front(), entry.op,
                               CurvedBooleanOptions{kTolerance});
      boolean_samples.push_back(ms_since(boolean_start));
      if (!result.has_value()) {
        print_error(result.error());
        return 2;
      }
      const Brep& brep = result.value().brep;
      const auto volume_start = std::chrono::steady_clock::now();
      const std::optional<double> volume = chart_volume(brep);
      volume_samples.push_back(ms_since(volume_start));
      if (!volume.has_value()) {
        std::printf("error=%s: a face is not a polar cap\n", entry.name);
        verified = false;
        continue;
      }
      const double value = *volume;
      checksum = checksum_volume(checksum, value);
      corpus_checksum = checksum_volume(corpus_checksum, value);
      if (run == 0) {
        baseline_bits = std::bit_cast<std::uint64_t>(value);
      } else if (std::bit_cast<std::uint64_t>(value) != baseline_bits) {
        ++volume_mismatches;
      }
      const double error = std::abs(value - entry.expected);
      if (!(error <= entry.relative_tolerance * std::abs(entry.expected))) {
        std::printf("error=%s: volume %.12f expected %.12f residual %.3e\n", entry.name, value,
                    entry.expected, error);
        verified = false;
      }
      if (run == 0 || run == kRuns - 1) {
        if (!verify_topology(entry.name, brep)) verified = false;
      }
      faces = 0;
      for (const ShellId shell : brep.solid_shells(brep.solid_ids().front())) {
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

    // The residual is printed relative to the analytic value; a value that is
    // bit-identical across runs but wrong still fails `verified`.
    const double residual = std::abs(measured_volume - entry.expected);
    std::printf("m3c3b-sphere-lens case=%s runs=%d ", entry.name, kRuns);
    print_times("boolean", boolean_samples);
    std::printf(" ");
    print_times("volume", volume_samples);
    std::printf(
        " faces=%zu expected_volume=%.15f measured_volume=%.15f residual=%.3e cuts=%zu kept=%zu "
        "max_cut_residual=%.6g checksum=0x%016llx verified=%s\n",
        faces, entry.expected, measured_volume, residual, cuts, kept, max_cut_residual,
        static_cast<unsigned long long>(checksum), verified ? "yes" : "no");
    std::fflush(stdout);
  }

  std::printf(
      "m3c3b-sphere-lens compiler=%s os=%s pointer_bits=%zu runs=%d checksum=0x%016llx "
      "verified=%s\n",
      compiler_name(), os_name(), sizeof(void*) * 8u, kRuns,
      static_cast<unsigned long long>(corpus_checksum), verified ? "yes" : "no");
  std::fflush(stdout);

  return verified ? 0 : 2;
}
