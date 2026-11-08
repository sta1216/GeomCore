// QA-10: tessellation throughput on 16,666 unit boxes (100k faces nominal).
// Builds the boxes once, tessellates the whole set three times and reports the
// fastest run. The mesh is also checked once for watertightness and outward
// winding; a miss is reported but never fails the process (benchmarks report).
#include <geom/brep/primitives.hpp>
#include <geom/brep/tessellate.hpp>

#include <array>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <map>
#include <utility>
#include <vector>

using namespace geom::v1;

namespace {

constexpr std::uint32_t kBoxes = 16666;
constexpr double kTargetMtri = 1.0;

double ms_since(std::chrono::steady_clock::time_point t0) {
  return std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - t0).count();
}

Vec3d box_center(std::uint32_t index) {
  const std::uint32_t group = index / 1000;
  const std::uint32_t within = index % 1000;
  const double x = 1.5 * static_cast<double>(within % 10);
  const double y = 1.5 * static_cast<double>((within / 10) % 10);
  const double z = 1.5 * static_cast<double>(within / 100);
  return Vec3d{x + 20.0 * static_cast<double>(group), y, z};
}

}  // namespace

int main() {
  std::vector<Brep> boxes;
  boxes.reserve(kBoxes);
  for (std::uint32_t i = 0; i < kBoxes; ++i) {
    const Vec3d center = box_center(i);
    const Result<Brep> box =
        make_box({center[0] - 0.5, center[1] - 0.5, center[2] - 0.5},
                 {center[0] + 0.5, center[1] + 0.5, center[2] + 0.5}, 0.0);
    if (!box.has_value()) {
      std::printf("QA10 build error code=%d message=%s\n", static_cast<int>(box.error().code),
                  box.error().message);
      return 0;
    }
    boxes.push_back(box.value());
  }

  const TessellationOptions options{0.01, 0.35};
  double best_ms = 1e300;
  std::size_t triangles = 0;
  bool watertight = true;
  bool outward = true;
  for (int run = 0; run < 3; ++run) {
    std::size_t run_triangles = 0;
    const auto t0 = std::chrono::steady_clock::now();
    for (std::size_t i = 0; i < boxes.size(); ++i) {
      const Result<TriMesh> mesh = tessellate(boxes[i], options);
      if (!mesh.has_value()) {
        std::printf("QA10 tessellation error\n");
        return 0;
      }
      run_triangles += mesh.value().triangles.size();
      if (run == 0) {
        const TriMesh& m = mesh.value();
        if (m.triangles.size() != 12 || m.positions.size() != 8) watertight = false;
        std::map<std::pair<std::uint32_t, std::uint32_t>, int> use;
        for (const auto& t : m.triangles) {
          for (int k = 0; k < 3; ++k) {
            std::uint32_t a = t[k];
            std::uint32_t b = t[(k + 1) % 3];
            if (a > b) std::swap(a, b);
            ++use[{a, b}];
          }
          const Vec3d ab = m.positions[t[1]] - m.positions[t[0]];
          const Vec3d ac = m.positions[t[2]] - m.positions[t[0]];
          const Vec3d normal = cross(ab, ac);
          if (!(dot(normal, normal) > 0.0)) outward = false;
          const Point3d& pa = m.positions[t[0]];
          const Point3d& pb = m.positions[t[1]];
          const Point3d& pc = m.positions[t[2]];
          const Vec3d centroid{(pa[0] + pb[0] + pc[0]) / 3.0, (pa[1] + pb[1] + pc[1]) / 3.0,
                               (pa[2] + pb[2] + pc[2]) / 3.0};
          const Vec3d outward_dir = centroid - (box_center(static_cast<std::uint32_t>(i)));
          if (!(dot(normal, outward_dir) > 0.0)) outward = false;
        }
        for (const auto& [edge, count] : use) {
          if (count != 2) watertight = false;
        }
      }
    }
    const double ms = ms_since(t0);
    if (ms < best_ms) best_ms = ms;
    triangles = run_triangles;
  }

  const double mtri = static_cast<double>(triangles) / (best_ms / 1000.0) / 1e6;
  std::printf("QA10 boxes=%u faces=%u tris=%zu ms=%.2f Mtri/s=%.2f\n", kBoxes,
              kBoxes * 6u, triangles, best_ms, mtri);
  std::printf("QA10 watertight=%s outward=%s target=%.1f Mtri/s %s\n", watertight ? "yes" : "no",
              outward ? "yes" : "no", kTargetMtri, mtri >= kTargetMtri ? "PASS" : "FAIL");
  return 0;
}
