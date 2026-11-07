#include <doctest/doctest.h>

#include <geom/index/bvh.hpp>

#include <array>
#include <cmath>
#include <cstdint>
#include <thread>
#include <vector>

using geom::v1::Bvh;
using geom::v1::Point3d;
using geom::v1::Ray3d;
using geom::v1::Triangle3d;

TEST_CASE("qa2: 4 threads x 1e6 rays agree with serial") {
  std::vector<Triangle3d> tris;
  for (int y = 0; y < 100; ++y) {
    for (int x = 0; x < 100; ++x) {
      const double fx = static_cast<double>(x), fy = static_cast<double>(y);
      tris.push_back(Triangle3d{Point3d{fx, fy, 0}, Point3d{fx + 0.9, fy, 0},
                                Point3d{fx, fy + 0.9, 0}});
    }
  }
  Bvh bvh;
  bvh.build(tris);

  constexpr std::uint32_t kRays = 1000000;
  auto hash_of = [&](std::uint32_t seed) {
    std::uint64_t h = 1469598103934665603ull;
    for (std::uint32_t i = 0; i < kRays; ++i) {
      seed = seed * 1664525u + 1013904223u;
      const double ox = (seed % 10000) * 0.01;
      seed = seed * 1664525u + 1013904223u;
      const double oy = (seed % 10000) * 0.01;
      const Ray3d ray{Point3d{ox, oy, -1}, geom::v1::Vec3d{0, 0, 1}};
      Bvh::Hit hit;
      const bool ok = bvh.raycast(ray, 100.0, hit);
      h ^= ok ? (static_cast<std::uint64_t>(hit.prim) + 1) : 0;
      h *= 1099511628211ull;
    }
    return h;
  };

  std::array<std::uint64_t, 4> parallel{};
  std::array<std::thread, 4> workers;
  for (std::uint32_t t = 0; t < 4; ++t) {
    workers[t] = std::thread([&, t] { parallel[t] = hash_of(1000 + t); });
  }
  for (auto& w : workers) w.join();
  for (std::uint32_t t = 0; t < 4; ++t) CHECK(parallel[t] == hash_of(1000 + t));
}
