// M4c-3c benchmark: round all twelve edges of a 100 mm cube with one radius.
//
// The timed region is one `fillet_edges` call per sample: recognition of the
// twelve-edge selection plus the 26-fragment construction (6 inset planes, 12
// quarter cylinders, 8 spherical octants). No tessellation is timed or
// validated here; watertightness is a unit-test claim, not a benchmark claim.

#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstdio>
#include <vector>

#include <geom/brep/detail/curved_volume.hpp>
#include <geom/brep/fillet.hpp>
#include <geom/brep/primitives.hpp>
#include <geom/brep/topology.hpp>
#include <geom/brep/validate.hpp>
#include <geom/core/error.hpp>
#include <geom/primitives/point.hpp>

namespace {

using Clock = std::chrono::steady_clock;

double seconds_since(Clock::time_point start) {
  return std::chrono::duration<double>(Clock::now() - start).count();
}

constexpr double kSize = 100.0;
constexpr double kRadius = 1.0;
constexpr double kTolerance = 1e-9;
constexpr std::size_t kSamples = 5;

// Steiner volume with A = B = C = 98:
//   ABC + 2(AB + BC + CA)r + pi(A + B + C)r^2 + 4pi r^3 / 3
//   = 941192 + 57624 + 923.628240155... + 4.18879... = 999743.8170303602.
constexpr double kExpectedVolume = 999743.8170303602;

struct Stats {
  double min_ms = 0.0;
  double median_ms = 0.0;
  double max_ms = 0.0;
};

Stats summarize(std::array<double, kSamples> samples_ms) {
  std::sort(samples_ms.begin(), samples_ms.end());
  return Stats{samples_ms.front(), samples_ms[kSamples / 2], samples_ms.back()};
}

}  // namespace

int main() {
  const geom::Result<geom::Brep> box = geom::make_box(
      geom::Point3d{0.0, 0.0, 0.0}, geom::Point3d{kSize, kSize, kSize}, kTolerance);
  if (!box.has_value()) {
    std::fprintf(stderr, "make_box(100 cube) failed: %s\n",
                 geom::to_string(box.error().code));
    return 2;
  }

  const std::vector<geom::EdgeId> edges = box.value().edge_ids();
  if (edges.size() != 12u) {
    std::fprintf(stderr, "expected 12 box edges, found %zu\n", edges.size());
    return 2;
  }

  const geom::FilletOptions options{kRadius, kTolerance};
  std::array<double, kSamples> samples{};
  double volume = 0.0;

  for (std::size_t i = 0; i < kSamples; ++i) {
    const Clock::time_point start = Clock::now();
    const geom::Result<geom::FilletResult> result =
        geom::fillet_edges(box.value(), edges, options);
    samples[i] = seconds_since(start) * 1000.0;

    if (!result.has_value()) {
      std::fprintf(stderr, "sample %zu: fillet_edges failed: %s: %s\n", i,
                   geom::to_string(result.error().code), result.error().message);
      return 2;
    }
    const geom::FilletResult& fillet = result.value();
    if (fillet.filleted_edges != 12u) {
      std::fprintf(stderr, "sample %zu: filleted %u of %zu edges\n", i,
                   fillet.filleted_edges, edges.size());
      return 2;
    }
    const geom::Brep& rounded = fillet.brep;
    if (rounded.solid_ids().size() != 1u || rounded.shell_count() != 1u) {
      std::fprintf(stderr,
                   "sample %zu: fillet result has %zu solid(s) and %u shell(s), expected 1/1\n",
                   i, rounded.solid_ids().size(), rounded.shell_count());
      return 2;
    }
    if (rounded.vertex_count() != 24u || rounded.edge_count() != 48u ||
        rounded.face_count() != 26u) {
      std::fprintf(stderr, "sample %zu: fillet topology counts V/E/F = %u/%u/%u, expected 24/48/26\n",
                   i, rounded.vertex_count(), rounded.edge_count(), rounded.face_count());
      return 2;
    }
    const geom::ValidationReport report = geom::validate_topology(rounded);
    if (!report.ok()) {
      std::fprintf(stderr, "sample %zu: validate_topology(rounded box) reported %zu issue(s)\n",
                   i, report.issues.size());
      return 2;
    }
    const geom::Result<double> computed =
        geom::detail::curved_solid_volume(rounded, rounded.solid_ids().front(), kTolerance);
    if (!computed.has_value()) {
      std::fprintf(stderr, "sample %zu: curved_solid_volume(rounded box) failed: %s\n", i,
                   geom::to_string(computed.error().code));
      return 2;
    }
    volume = computed.value();
    if (std::abs(volume - kExpectedVolume) > 1e-9 * std::abs(kExpectedVolume)) {
      std::fprintf(stderr, "sample %zu: volume checksum mismatch: %.15f (expected %.15f)\n", i,
                   volume, kExpectedVolume);
      return 2;
    }
  }

  std::printf("rounded box volume: %.15f (expected %.15f)\n", volume, kExpectedVolume);
  const Stats stats = summarize(samples);
  std::printf("fillet 12 box edges r=%.1f: min %.3f ms, median %.3f ms, max %.3f ms "
              "(%zu runs)\n",
              kRadius, stats.min_ms, stats.median_ms, stats.max_ms, samples.size());
  return 0;
}
