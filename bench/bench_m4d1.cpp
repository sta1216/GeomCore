// M4d-1 profile extrusion benchmark.
//
// Two workloads exercise `extrude_profile` on the shapes the milestone accepts:
// a 1024-gon prism (the plan's headline case) and a 512-gon ring with a 256-gon
// hole (the shape used to measure the certificate-triangulation reuse). Each
// workload is checked, not just timed: topology counts, `validate_topology`, the
// analytic volume and a non-empty tessellation must all hold, otherwise the
// benchmark exits 2.
//
// A third timing reports `triangulate_polygon2` on the holed profile, the one
// remaining triangulation inside `extrude_profile`. Together the two numbers
// reproduce the reuse claim: the extrusion costs about one certificate, not
// three triangulations.

#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <numbers>
#include <vector>

#include <geom/brep/extrude.hpp>
#include <geom/brep/mass_properties.hpp>
#include <geom/brep/tessellate.hpp>
#include <geom/brep/topology.hpp>
#include <geom/brep/validate.hpp>
#include <geom/cgeom/polygon2.hpp>
#include <geom/cgeom/triangulate2.hpp>
#include <geom/surfaces/plane.hpp>

namespace {

using namespace geom;

constexpr int kSamples = 5;
constexpr double kTolerance = 1e-12;
constexpr double kPi = std::numbers::pi_v<double>;

using Clock = std::chrono::steady_clock;

double seconds_since(Clock::time_point start) {
  return std::chrono::duration<double>(Clock::now() - start).count();
}

// Vertices of a counter-clockwise regular ring.
std::vector<Point2d> regular_ring(int sides, double radius) {
  std::vector<Point2d> ring;
  ring.reserve(static_cast<std::size_t>(sides));
  for (int i = 0; i < sides; ++i) {
    const double angle = 2.0 * kPi * static_cast<double>(i) /
                         static_cast<double>(sides);
    ring.push_back(Point2d{radius * std::cos(angle), radius * std::sin(angle)});
  }
  return ring;
}

Polygon2d regular_polygon(int sides, double radius) {
  Polygon2d profile;
  profile.vertices = regular_ring(sides, radius);
  profile.ring_starts = {0u, static_cast<std::uint32_t>(sides)};
  return profile;
}

// Outer ring followed by one hole ring, both wound counter-clockwise; the
// extrusion normalizes the hole to clockwise itself.
Polygon2d ring_with_hole(int outer_sides, double outer_radius, int hole_sides,
                         double hole_radius) {
  Polygon2d profile;
  profile.vertices = regular_ring(outer_sides, outer_radius);
  const std::vector<Point2d> hole = regular_ring(hole_sides, hole_radius);
  profile.vertices.insert(profile.vertices.end(), hole.begin(), hole.end());
  profile.ring_starts = {0u,
                         static_cast<std::uint32_t>(outer_sides),
                         static_cast<std::uint32_t>(outer_sides + hole_sides)};
  return profile;
}

// Absolute area of the first ring minus the absolute area of every hole ring,
// independent of the stored winding.
double profile_area(const Polygon2d& profile) {
  double area = 0.0;
  for (std::size_t r = 0; r + 1u < profile.ring_starts.size(); ++r) {
    const auto begin = profile.vertices.cbegin() + profile.ring_starts[r];
    const auto end = profile.vertices.cbegin() + profile.ring_starts[r + 1u];
    double twice = 0.0;
    for (auto it = begin; it != end; ++it) {
      const auto next = (it + 1 == end) ? begin : it + 1;
      twice += (*it)[0] * (*next)[1] - (*next)[0] * (*it)[1];
    }
    const double ring_area = 0.5 * std::abs(twice);
    area += r == 0u ? ring_area : -ring_area;
  }
  return area;
}

struct Stats {
  double min_ms = 0.0;
  double median_ms = 0.0;
  double max_ms = 0.0;
};

Stats summarize(std::array<double, kSamples> samples_ms) {
  std::sort(samples_ms.begin(), samples_ms.end());
  return Stats{samples_ms.front(), samples_ms[kSamples / 2], samples_ms.back()};
}

struct Workload {
  const char* name;
  Polygon2d profile;
  double height;
  double expected_volume;
  std::uint32_t vertices;
  std::uint32_t edges;
  std::uint32_t faces;
};

// Extrudes a workload kSamples times and checks every run. Returns 0 on
// success, 2 on any mismatch.
int run_workload(const Workload& workload, const Plane3& plane) {
  std::array<double, kSamples> samples{};
  double volume = 0.0;
  for (std::size_t i = 0; i < samples.size(); ++i) {
    const Clock::time_point start = Clock::now();
    const Result<Brep> solid = extrude_profile(plane, workload.profile,
                                               workload.height, kTolerance);
    samples[i] = seconds_since(start) * 1000.0;
    if (!solid.has_value()) {
      std::fprintf(stderr, "%s: extrude failed: %s\n", workload.name,
                   solid.error().message);
      return 2;
    }
    const Brep& brep = solid.value();
    if (brep.vertex_count() != workload.vertices ||
        brep.edge_count() != workload.edges ||
        brep.face_count() != workload.faces) {
      std::fprintf(stderr, "%s: topology V/E/F = %u/%u/%u, expected %u/%u/%u\n",
                   workload.name, brep.vertex_count(), brep.edge_count(),
                   brep.face_count(), workload.vertices, workload.edges,
                   workload.faces);
      return 2;
    }
    const ValidationReport report = validate_topology(brep);
    if (!report.ok()) {
      std::fprintf(stderr, "%s: validate_topology reported %zu issue(s)\n",
                   workload.name, report.issues.size());
      return 2;
    }
    const Result<MassProperties> properties = mass_properties(brep);
    if (!properties.has_value()) {
      std::fprintf(stderr, "%s: mass_properties failed\n", workload.name);
      return 2;
    }
    volume = properties.value().volume;
    if (std::abs(volume - workload.expected_volume) >
        1e-9 * workload.expected_volume) {
      std::fprintf(stderr,
                   "%s: volume %.17g, expected %.17g (relative error %.3g)\n",
                   workload.name, volume, workload.expected_volume,
                   std::abs(volume - workload.expected_volume) /
                       workload.expected_volume);
      return 2;
    }
    const Result<TriMesh> mesh = tessellate(brep, TessellationOptions{});
    if (!mesh.has_value() || mesh.value().triangles.empty()) {
      std::fprintf(stderr, "%s: tessellation failed\n", workload.name);
      return 2;
    }
  }

  const Stats stats = summarize(samples);
  std::printf("%s: min %.3f ms, median %.3f ms, max %.3f ms (%d runs)\n",
              workload.name, stats.min_ms, stats.median_ms, stats.max_ms,
              kSamples);
  std::printf("%s: volume %.15f (expected %.15f)\n", workload.name, volume,
              workload.expected_volume);
  return 0;
}

}  // namespace

int main() {
  constexpr int kSides = 1024;
  constexpr double kRadius = 10.0;
  constexpr double kHeight = 1.0;

  const Polygon2d disc = regular_polygon(kSides, kRadius);
  const Plane3 plane{Vec3d{0.0, 0.0, 1.0}, 0.0};
  const double disc_area =
      0.5 * static_cast<double>(kSides) * kRadius * kRadius *
      std::sin(2.0 * kPi / static_cast<double>(kSides));

  const Workload disc_workload{"extrude 1024-gon h=1", disc, kHeight,
                               disc_area * kHeight,
                               static_cast<std::uint32_t>(2 * kSides),
                               static_cast<std::uint32_t>(3 * kSides),
                               static_cast<std::uint32_t>(kSides + 2)};
  if (const int status = run_workload(disc_workload, plane); status != 0) {
    return status;
  }

  constexpr int kOuterSides = 512;
  constexpr int kHoleSides = 256;
  constexpr double kHoleRadius = 5.0;
  const Polygon2d holed =
      ring_with_hole(kOuterSides, kRadius, kHoleSides, kHoleRadius);
  const Workload hole_workload{
      "extrude 512-gon + 256-gon hole h=1", holed, kHeight,
      profile_area(holed) * kHeight,
      static_cast<std::uint32_t>(2 * (kOuterSides + kHoleSides)),
      static_cast<std::uint32_t>(3 * (kOuterSides + kHoleSides)),
      static_cast<std::uint32_t>(kOuterSides + kHoleSides + 2)};
  if (const int status = run_workload(hole_workload, plane); status != 0) {
    return status;
  }

  // The certificate: the only triangulation `extrude_profile` still performs.
  std::array<double, kSamples> certificate_samples{};
  std::size_t certificate_triangles = 0u;
  for (std::size_t i = 0; i < certificate_samples.size(); ++i) {
    const Clock::time_point start = Clock::now();
    const Triangulation2 triangulation = triangulate_polygon2(holed);
    certificate_samples[i] = seconds_since(start) * 1000.0;
    if (!triangulation.ok || triangulation.triangles.empty()) {
      std::fprintf(stderr, "certificate: triangulate_polygon2 failed\n");
      return 2;
    }
    certificate_triangles = triangulation.triangles.size();
  }
  const Stats certificate = summarize(certificate_samples);
  std::printf("certificate 512-gon + 256-gon hole: min %.3f ms, median %.3f "
              "ms, max %.3f ms (%d runs, %zu triangles)\n",
              certificate.min_ms, certificate.median_ms, certificate.max_ms,
              kSamples, certificate_triangles);
  return 0;
}
