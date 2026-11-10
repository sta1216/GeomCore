// QA-9: sphere-plane surface-surface intersection throughput.
//
// Builds a unit sphere patch (u in [0, 2*pi], v in [-1.5, 1.5]) and the plane
// z = 0.1 over [-2, 2]^2, both well conditioned and crossing in one closed
// circle, then times exactly `ssi()` for each of the requested calls (default
// 1000). Patch construction, result checks and statistics are outside the timed
// region, so each sample is one whole intersection request. Every timed call
// must return at least one curve; a call that does not marks the run unusable
// (FAIL, exit 1) rather than a timing result. Setup failures and ssi() errors
// exit 2.
//
// The gate is the median call: 100 microseconds per pair (§13.1 of the
// architecture). min/median/max are all printed so a missed gate can be read
// without rerunning, and the target is not relaxed when it is missed.
//
// Exit codes: 0 gate passed, 1 gate missed or a result was unusable, 2 setup or
// SSI error. Progress goes to stderr and the machine-readable QA9 line plus the
// curve/segment/residual summary go to stdout.
#include <geom/brep/ssi.hpp>

#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <numbers>
#include <string_view>
#include <vector>

using namespace geom::v1;

namespace {

constexpr double kPi = std::numbers::pi;
constexpr double kTolerance = 1e-9;
constexpr double kTargetUs = 100.0;
constexpr int kDefaultCalls = 1000;
constexpr int kMinCalls = 1;
constexpr int kMaxCalls = 1000000;

double us_since(std::chrono::steady_clock::time_point start) {
  return std::chrono::duration<double, std::micro>(
             std::chrono::steady_clock::now() - start)
      .count();
}

void print_error(const Error& error) {
  std::printf("error=%s: %s\n", to_string(error.code), error.message);
  std::fflush(stdout);
}

bool parse_calls(std::string_view text, int* out) {
  if (text.empty()) return false;
  int value = 0;
  for (const char digit : text) {
    if (digit < '0' || digit > '9') return false;
    value = value * 10 + (digit - '0');
    if (value > kMaxCalls) return false;
  }
  if (value < kMinCalls) return false;
  *out = value;
  return true;
}

}  // namespace

int main(int argc, char** argv) {
  if (argc > 2) {
    std::printf("error=usage: geom_bench_ssi [calls]\n");
    return 2;
  }

  int calls = kDefaultCalls;
  if (argc > 1 && !parse_calls(argv[1], &calls)) {
    std::printf("error=usage: calls must be an integer in [%d, %d]\n", kMinCalls,
                kMaxCalls);
    return 2;
  }

  const Result<SsiPatch> sphere = make_ssi_patch(
      SurfaceGeometry{Sphere3{Point3d{0, 0, 0}, Vec3d{0, 0, 1}, Vec3d{1, 0, 0}, 1.0}},
      {0.0, 2.0 * kPi}, {-1.5, 1.5});
  if (!sphere.has_value()) {
    print_error(sphere.error());
    return 2;
  }
  // z = 0.1 crosses the full sphere latitude circle; the domain is larger than
  // the circle radius sqrt(1 - 0.1^2) so every point of the curve is interior.
  const Result<SsiPatch> plane = make_ssi_patch(
      SurfaceGeometry{Plane3{{0, 0, 1}, 0.1}}, {-2.0, 2.0}, {-2.0, 2.0});
  if (!plane.has_value()) {
    print_error(plane.error());
    return 2;
  }

  const SsiOptions options{kTolerance, 12, 100000};

  // One untimed call warms the allocation paths and validates the fixture; its
  // result is discarded and it is not part of the measured sample count.
  const Result<SsiResult> warmup = ssi(sphere.value(), plane.value(), options);
  if (!warmup.has_value()) {
    print_error(warmup.error());
    return 2;
  }
  if (warmup.value().curves.empty()) {
    std::printf("error=ssi: warmup produced no intersection curve\n");
    return 2;
  }

  std::vector<double> samples(static_cast<std::size_t>(calls));
  std::uint64_t curve_count = 0;
  std::uint64_t segment_count = 0;
  std::uint64_t candidate_pairs = 0;
  std::uint64_t rejected_segments = 0;
  double max_residual = 0.0;
  bool unusable = false;
  for (int i = 0; i < calls; ++i) {
    const auto start = std::chrono::steady_clock::now();
    const Result<SsiResult> result = ssi(sphere.value(), plane.value(), options);
    samples[static_cast<std::size_t>(i)] = us_since(start);
    if (!result.has_value()) {
      print_error(result.error());
      return 2;
    }
    if (result.value().curves.empty()) {
      std::printf("error=ssi: call %d produced no intersection curve\n", i);
      unusable = true;
      continue;
    }
    curve_count += result.value().curves.size();
    segment_count += result.value().accepted_segments;
    candidate_pairs += result.value().candidate_pairs;
    rejected_segments += result.value().rejected_segments;
    max_residual = std::max(max_residual, result.value().max_residual);
  }

  std::vector<double> sorted = samples;
  std::sort(sorted.begin(), sorted.end());
  const double min_us = sorted.front();
  const double max_us = sorted.back();
  const double median_us = sorted.size() % 2 == 1
                               ? sorted[sorted.size() / 2]
                               : 0.5 * (sorted[sorted.size() / 2 - 1] +
                                        sorted[sorted.size() / 2]);
  const bool passed = !unusable && median_us <= kTargetUs;

  std::fprintf(stderr, "[bench] sphere-plane calls=%d tol=%g\n", calls, kTolerance);
  std::fflush(stderr);
  std::printf(
      "QA9 calls=%d min_us=%.3f median_us=%.3f max_us=%.3f target_us=%.0f %s\n",
      calls, min_us, median_us, max_us, kTargetUs, passed ? "PASS" : "FAIL");
  std::printf(
      "QA9 curves=%llu segments=%llu candidate_pairs=%llu rejected_segments=%llu "
      "max_residual=%.3e\n",
      static_cast<unsigned long long>(curve_count),
      static_cast<unsigned long long>(segment_count),
      static_cast<unsigned long long>(candidate_pairs),
      static_cast<unsigned long long>(rejected_segments), max_residual);
  std::fflush(stdout);

  return passed ? 0 : 1;
}
