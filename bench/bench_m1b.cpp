#include <geom/cgeom/boolean2.hpp>
#include <geom/cgeom/polygon2.hpp>
#include <geom/index/bvh.hpp>

#include <chrono>
#include <cmath>
#include <cstdio>
#include <random>
#include <vector>

using namespace geom::v1;

namespace {

double ms_since(std::chrono::steady_clock::time_point t0) {
  return std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - t0)
      .count();
}

Polygon2d star(std::uint32_t n, double r0, double r1) {
  Polygon2d p;
  p.vertices.resize(n);
  for (std::uint32_t i = 0; i < n; ++i) {
    const double a = 2.0 * 3.14159265358979323846 * i / n;
    const double r = (i % 2 == 0) ? r0 : r1;
    p.vertices[i] = Point2d{std::cos(a) * r, std::sin(a) * r};
  }
  p.ring_starts = {0, n};
  return p;
}

}  // namespace

int main() {
  // QA-1: two 100k-vertex polygons, area identity |A| + |B| = |A u B| + |A n B|.
  Polygon2d a = star(50000, 100.0, 99.0);
  Polygon2d b = a;
  for (auto& q : b.vertices) q[0] += 0.35;
  std::vector<Polygon2d> ops{a, b};

  const auto t0 = std::chrono::steady_clock::now();
  auto uni = boolean2(ops, BooleanOp::union_);
  const double union_ms = ms_since(t0);
  auto inter = boolean2(ops, BooleanOp::intersection);
  if (!uni || !inter) {
    std::printf("QA-1 failed: boolean error\n");
    return 2;
  }
  const double area_a = polygon_area2(a);
  const double area_b = polygon_area2(b);
  const double lhs = area_a + area_b;
  const double rhs = polygon_set_area2(uni.value().polygons) +
                     polygon_set_area2(inter.value().polygons);
  const double rel = std::fabs(lhs - rhs) / std::fabs(lhs);
  std::printf("QA-1 union %.2f ms, area rel err %.3e (limit 1e-12, 80 ms)\n", union_ms,
              rel);

  // QA-4: 1M triangles, build and 8M rays.
  std::vector<Triangle3d> tris;
  tris.reserve(1000000);
  for (int y = 0; y < 500; ++y) {
    for (int x = 0; x < 1000; ++x) {
      const double fx = x, fy = y;
      tris.push_back(Triangle3d{Point3d{fx, fy, 0}, Point3d{fx + 1, fy, 0},
                                Point3d{fx, fy + 1, 0}});
      tris.push_back(Triangle3d{Point3d{fx + 1, fy, 0}, Point3d{fx + 1, fy + 1, 0},
                                Point3d{fx, fy + 1, 0}});
    }
  }
  Bvh bvh;
  const auto tb = std::chrono::steady_clock::now();
  bvh.build(tris);
  const double build_ms = ms_since(tb);

  std::mt19937 rng(20261107);
  std::uniform_real_distribution<double> uni_x(0.0, 1000.0);
  std::uniform_real_distribution<double> uni_y(0.0, 500.0);
  constexpr std::uint32_t kRays = 8000000;
  std::uint64_t hits = 0;
  const auto tr = std::chrono::steady_clock::now();
  for (std::uint32_t i = 0; i < kRays; ++i) {
    const Ray3d ray{Point3d{uni_x(rng), uni_y(rng), -1}, Vec3d{0, 0, 1}};
    Bvh::Hit hit;
    if (bvh.raycast(ray, 10.0, hit)) ++hits;
  }
  const double ray_s = ms_since(tr);
  const double mrays = kRays / (ray_s / 1000.0) / 1e6;
  std::printf("QA-4 build %.1f ms (limit 350), %.2f Mray/s (limit 8), hits %llu\n",
              build_ms, mrays, static_cast<unsigned long long>(hits));

  bool ok = rel <= 1e-12 && union_ms <= 80.0 && build_ms <= 350.0 && mrays >= 8.0;
#ifndef NDEBUG
  ok = true;  // Debug builds report only.
#endif
  return ok ? 0 : 2;
}
