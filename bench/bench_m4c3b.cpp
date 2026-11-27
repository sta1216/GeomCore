#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstdio>
#include <numbers>
#include <span>
#include <vector>

#include <geom/brep/assemble.hpp>
#include <geom/brep/detail/curved_volume.hpp>
#include <geom/brep/fillet.hpp>
#include <geom/brep/topology.hpp>
#include <geom/brep/validate.hpp>
#include <geom/core/error.hpp>
#include <geom/primitives/point.hpp>
#include <geom/primitives/vec.hpp>

namespace {

using Clock = std::chrono::steady_clock;

double seconds_since(Clock::time_point start) {
  return std::chrono::duration<double>(Clock::now() - start).count();
}

constexpr double kTolerance = 1e-9;
constexpr double kRadius = 0.05;
constexpr double kPi = std::numbers::pi_v<double>;
constexpr int kSides = 64;

// Area removed from the cross-section by one filleted corner of interior angle
// theta: the tangent triangle minus the circular sector.
double corner_area(double radius, double theta) {
  return radius * radius *
         (1.0 / std::tan(theta / 2.0) - (kPi - theta) / 2.0);
}

// Same fixture as tests/unit/fillet_prism_test.cpp::make_prism: a CCW polygon
// in z = [0, height], clockwise bottom about -Z, outward-facing side quads.
geom::Result<geom::Brep> make_prism(const std::span<const geom::Point2d> polygon,
                                    double height, double tolerance) {
  std::vector<geom::AssembledFace> faces;
  std::vector<geom::Point3d> bottom;
  std::vector<geom::Point3d> top;
  bottom.reserve(polygon.size());
  top.reserve(polygon.size());
  for (const geom::Point2d& p : polygon) {
    bottom.push_back(geom::Point3d{p[0], p[1], 0.0});
    top.push_back(geom::Point3d{p[0], p[1], height});
  }

  geom::AssembledFace bottom_face;
  bottom_face.outer.assign(bottom.rbegin(), bottom.rend());
  bottom_face.normal = geom::Vec3d{0.0, 0.0, -1.0};
  faces.push_back(std::move(bottom_face));

  geom::AssembledFace top_face;
  top_face.outer = top;
  top_face.normal = geom::Vec3d{0.0, 0.0, 1.0};
  faces.push_back(std::move(top_face));

  for (std::size_t i = 0; i < polygon.size(); ++i) {
    const std::size_t j = (i + 1) % polygon.size();
    const geom::Vec3d edge{polygon[j][0] - polygon[i][0], polygon[j][1] - polygon[i][1], 0.0};
    geom::AssembledFace side;
    side.normal = geom::detail::assemble_unit(geom::Vec3d{edge[1], -edge[0], 0.0});
    side.outer = {bottom[i], bottom[j], top[j], top[i]};
    faces.push_back(std::move(side));
  }
  return geom::assemble_solid(faces, tolerance);
}

// A vertical prism edge joins two ring vertices with the same x/y and different
// z. No other prism edge has that property.
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

struct Stats {
  double min_ms = 0.0;
  double median_ms = 0.0;
  double max_ms = 0.0;
};

Stats summarize(std::array<double, 5> samples_ms) {
  std::sort(samples_ms.begin(), samples_ms.end());
  return Stats{samples_ms.front(), samples_ms[2], samples_ms.back()};
}

}  // namespace

int main() {
  constexpr double kHeight = 10.0;
  std::array<geom::Point2d, kSides> polygon{};
  for (int i = 0; i < kSides; ++i) {
    const double angle = 2.0 * kPi * static_cast<double>(i) / static_cast<double>(kSides);
    polygon[static_cast<std::size_t>(i)] =
        geom::Point2d{2.0 * std::cos(angle), 2.0 * std::sin(angle)};
  }

  const geom::Result<geom::Brep> prism = make_prism(polygon, kHeight, kTolerance);
  if (!prism.has_value()) {
    std::fprintf(stderr, "assemble_solid(prism) failed: %s\n",
                 geom::to_string(prism.error().code));
    return 2;
  }

  const std::vector<geom::EdgeId> edges = vertical_edges(prism.value());
  if (edges.size() != static_cast<std::size_t>(kSides)) {
    std::fprintf(stderr, "expected %d vertical edges, found %zu\n", kSides, edges.size());
    return 2;
  }
  const geom::FilletOptions options{kRadius, kTolerance};

  const double theta = kPi - 2.0 * kPi / static_cast<double>(kSides);
  double area = 0.0;
  for (int i = 0; i < kSides; ++i) {
    const geom::Point2d& a = polygon[static_cast<std::size_t>(i)];
    const geom::Point2d& b = polygon[static_cast<std::size_t>((i + 1) % kSides)];
    area += a[0] * b[1] - b[0] * a[1];
  }
  area = 0.5 * std::abs(area);
  const double expected_volume =
      area * kHeight -
      static_cast<double>(kSides) * kHeight * corner_area(kRadius, theta);

  std::array<double, 5> samples{};
  double volume = 0.0;
  for (std::size_t i = 0; i < samples.size(); ++i) {
    const Clock::time_point start = Clock::now();
    const geom::Result<geom::FilletResult> result =
        geom::fillet_edges(prism.value(), edges, options);
    samples[i] = seconds_since(start) * 1000.0;
    if (!result.has_value()) {
      std::fprintf(stderr, "fillet_edges failed: %s: %s\n",
                   geom::to_string(result.error().code), result.error().message);
      return 2;
    }
    if (result.value().filleted_edges != static_cast<std::uint32_t>(kSides)) {
      std::fprintf(stderr, "filleted %u of %zu edges\n", result.value().filleted_edges,
                   edges.size());
      return 2;
    }
    const geom::Brep& filleted = result.value().brep;
    if (filleted.solid_ids().size() != 1u || filleted.shell_count() != 1u) {
      std::fprintf(stderr, "fillet result has %zu solid(s) and %zu shell(s), expected 1/1\n",
                   filleted.solid_ids().size(),
                   static_cast<std::size_t>(filleted.shell_count()));
      return 2;
    }
    if (filleted.vertex_count() != 256u || filleted.edge_count() != 384u ||
        filleted.face_count() != 130u) {
      std::fprintf(stderr, "fillet topology counts V/E/F = %u/%u/%u, expected 256/384/130\n",
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
        geom::detail::curved_solid_volume(filleted, filleted.solid_ids().front(), kTolerance);
    if (!computed.has_value()) {
      std::fprintf(stderr, "curved_solid_volume(fillet) failed: %s\n",
                   geom::to_string(computed.error().code));
      return 2;
    }
    volume = computed.value();
    if (std::abs(volume - expected_volume) > 1e-9 * std::abs(expected_volume)) {
      std::fprintf(stderr, "volume checksum mismatch: %.15f (expected %.15f)\n", volume,
                   expected_volume);
      return 2;
    }
  }

  std::printf("64-gon prism fillet volume: %.15f (expected %.15f)\n", volume, expected_volume);
  const Stats stats = summarize(samples);
  std::printf("fillet 64 vertical edges r=0.05: min %.3f ms, median %.3f ms, max %.3f ms "
              "(%zu runs)\n",
              stats.min_ms, stats.median_ms, stats.max_ms, samples.size());
  return 0;
}
