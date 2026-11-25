#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <vector>

#include <geom/brep/chamfer.hpp>
#include <geom/brep/mass_properties.hpp>
#include <geom/brep/primitives.hpp>
#include <geom/brep/shell.hpp>
#include <geom/core/error.hpp>

namespace {

using Clock = std::chrono::steady_clock;

double seconds_since(Clock::time_point start) {
  return std::chrono::duration<double>(Clock::now() - start).count();
}

struct Stats {
  double min_ms = 0.0;
  double median_ms = 0.0;
  double max_ms = 0.0;
};

Stats summarize(std::array<double, 5> samples_ms) {
  std::sort(samples_ms.begin(), samples_ms.end());
  return Stats{samples_ms.front(), samples_ms[2], samples_ms.back()};
}

void print_stats(const char* label, const Stats& stats, std::size_t runs) {
  std::printf("%s: min %.3f ms, median %.3f ms, max %.3f ms (%zu runs)\n", label,
              stats.min_ms, stats.median_ms, stats.max_ms, runs);
}

bool close_to(double actual, double expected) {
  return std::abs(actual - expected) <= 1e-6 * std::abs(expected);
}

}  // namespace

int main() {
  const geom::Result<geom::Brep> box =
      geom::make_box(geom::Point3d{0.0, 0.0, 0.0}, geom::Point3d{100.0, 100.0, 100.0}, 1e-9);
  if (!box.has_value()) {
    std::fprintf(stderr, "make_box failed: %s\n", geom::to_string(box.error().code));
    return 2;
  }

  const std::vector<geom::EdgeId> edges = box.value().edge_ids();
  const geom::ChamferOptions chamfer_options{1.0, 1e-9};

  std::array<double, 5> chamfer_samples{};
  double chamfer_volume = 0.0;
  for (std::size_t i = 0; i < chamfer_samples.size(); ++i) {
    const Clock::time_point start = Clock::now();
    const geom::Result<geom::ChamferResult> result =
        geom::chamfer_edges(box.value(), edges, chamfer_options);
    chamfer_samples[i] = seconds_since(start) * 1000.0;
    if (!result.has_value()) {
      std::fprintf(stderr, "chamfer_edges failed: %s: %s\n",
                   geom::to_string(result.error().code), result.error().message);
      return 2;
    }
    const geom::Result<geom::MassProperties> props = geom::mass_properties(result.value().brep);
    if (!props.has_value()) {
      std::fprintf(stderr, "mass_properties(chamfer) failed: %s\n",
                   geom::to_string(props.error().code));
      return 2;
    }
    chamfer_volume = props.value().volume;
    if (i == 0 && result.value().chamfered_edges != edges.size()) {
      std::fprintf(stderr, "chamfered %u of %zu edges\n", result.value().chamfered_edges,
                   edges.size());
      return 2;
    }
  }

  const geom::ShellOptions shell_options{1.0, 1e-9};
  const std::vector<geom::FaceId> open_faces;
  std::array<double, 5> shell_samples{};
  double shell_volume = 0.0;
  for (std::size_t i = 0; i < shell_samples.size(); ++i) {
    const Clock::time_point start = Clock::now();
    const geom::Result<geom::ShellResult> result =
        geom::shell_solid(box.value(), open_faces, shell_options);
    shell_samples[i] = seconds_since(start) * 1000.0;
    if (!result.has_value()) {
      std::fprintf(stderr, "shell_solid failed: %s: %s\n",
                   geom::to_string(result.error().code), result.error().message);
      return 2;
    }
    const geom::Result<geom::MassProperties> props = geom::mass_properties(result.value().brep);
    if (!props.has_value()) {
      std::fprintf(stderr, "mass_properties(shell) failed: %s\n",
                   geom::to_string(props.error().code));
      return 2;
    }
    shell_volume = props.value().volume;
    if (result.value().cavity_shells != 1) {
      std::fprintf(stderr, "expected 1 cavity shell, got %u\n", result.value().cavity_shells);
      return 2;
    }
  }

  std::printf("chamfer volume: %.12f (expected 999406)\n", chamfer_volume);
  std::printf("shell volume:   %.12f (expected 58808)\n", shell_volume);
  print_stats("chamfer all 12 edges", summarize(chamfer_samples), chamfer_samples.size());
  print_stats("closed shell", summarize(shell_samples), shell_samples.size());

  if (!close_to(chamfer_volume, 999406.0) || !close_to(shell_volume, 58808.0)) {
    std::fprintf(stderr, "volume checksum mismatch\n");
    return 3;
  }
  return 0;
}
