#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <cstddef>
#include <cstdio>
#include <vector>

#include <geom/brep/detail/curved_volume.hpp>
#include <geom/brep/fillet.hpp>
#include <geom/brep/primitives.hpp>
#include <geom/brep/validate.hpp>
#include <geom/core/error.hpp>

namespace {

using Clock = std::chrono::steady_clock;

double seconds_since(Clock::time_point start) {
  return std::chrono::duration<double>(Clock::now() - start).count();
}

// 100^3 minus four quarter-cylinder edge cutouts of radius 1 and length 100:
// 1,000,000 - 4 * (1 - pi/4) * 1^2 * 100 = 999,914.159265359.
constexpr double kExpectedVolume = 999914.159265359;

bool close_to(double actual, double expected) {
  return std::abs(actual - expected) <= 1e-6 * std::abs(expected);
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

// The four edges of the box that run parallel to +Z: both endpoints share x and
// y and differ in z. `make_box` emits no other edge with that property.
std::vector<geom::EdgeId> vertical_edges(const geom::Brep& brep) {
  std::vector<geom::EdgeId> edges;
  for (const geom::EdgeId id : brep.edge_ids()) {
    const geom::Edge& edge = brep.edge(id);
    const geom::Point3d& start = brep.vertex(edge.v0).point;
    const geom::Point3d& end = brep.vertex(edge.v1).point;
    const bool same_x = std::abs(end[0] - start[0]) <= 1e-12;
    const bool same_y = std::abs(end[1] - start[1]) <= 1e-12;
    const bool spans_z = std::abs(end[2] - start[2]) > 1e-12;
    if (same_x && same_y && spans_z) edges.push_back(id);
  }
  return edges;
}

}  // namespace

int main() {
  const geom::Result<geom::Brep> box =
      geom::make_box(geom::Point3d{0.0, 0.0, 0.0}, geom::Point3d{100.0, 100.0, 100.0}, 1e-9);
  if (!box.has_value()) {
    std::fprintf(stderr, "make_box failed: %s\n", geom::to_string(box.error().code));
    return 2;
  }

  const std::vector<geom::EdgeId> edges = vertical_edges(box.value());
  if (edges.size() != 4) {
    std::fprintf(stderr, "expected 4 vertical edges, found %zu\n", edges.size());
    return 2;
  }
  const geom::FilletOptions options{1.0, 1e-9};

  std::array<double, 5> samples{};
  double volume = 0.0;
  for (std::size_t i = 0; i < samples.size(); ++i) {
    const Clock::time_point start = Clock::now();
    const geom::Result<geom::FilletResult> result =
        geom::fillet_edges(box.value(), edges, options);
    samples[i] = seconds_since(start) * 1000.0;
    if (!result.has_value()) {
      std::fprintf(stderr, "fillet_edges failed: %s: %s\n",
                   geom::to_string(result.error().code), result.error().message);
      return 2;
    }
    if (result.value().filleted_edges != 4) {
      std::fprintf(stderr, "filleted %u of %zu edges\n", result.value().filleted_edges,
                   edges.size());
      return 2;
    }
    const geom::Brep& filleted = result.value().brep;
    if (filleted.solid_ids().empty()) {
      std::fprintf(stderr, "fillet result has no solid\n");
      return 2;
    }
    if (filleted.vertex_count() != 16 || filleted.edge_count() != 24 ||
        filleted.face_count() != 10) {
      std::fprintf(stderr, "fillet topology counts V/E/F = %u/%u/%u, expected 16/24/10\n",
                   filleted.vertex_count(), filleted.edge_count(), filleted.face_count());
      return 2;
    }
    const geom::ValidationReport report = geom::validate_topology(filleted);
    if (!report.ok()) {
      std::fprintf(stderr, "validate_topology(fillet) reported %zu issue(s)\n",
                   report.issues.size());
      return 2;
    }
    const geom::Result<double> computed =
        geom::detail::curved_solid_volume(filleted, filleted.solid_ids().front(), 1e-9);
    if (!computed.has_value()) {
      std::fprintf(stderr, "curved_solid_volume(fillet) failed: %s\n",
                   geom::to_string(computed.error().code));
      return 2;
    }
    volume = computed.value();
    if (!close_to(volume, kExpectedVolume)) {
      std::fprintf(stderr, "volume checksum mismatch: %.15f (expected %.15f)\n", volume,
                   kExpectedVolume);
      return 2;
    }
  }

  std::printf("four-edge fillet volume: %.15f (expected %.15f)\n", volume, kExpectedVolume);
  print_stats("fillet four vertical edges r=1", summarize(samples), samples.size());
  return 0;
}
