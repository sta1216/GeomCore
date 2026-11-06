#include <doctest/doctest.h>

#include <cmath>
#include <cstdint>
#include <limits>
#include <random>
#include <vector>

#include <geom/index/bvh.hpp>

using geom::v1::Bvh;
using geom::v1::Point3d;
using geom::v1::Ray3d;
using geom::v1::Triangle3d;
using geom::v1::Vec3d;

namespace {

bool ray_tri(const Ray3d& r, const Triangle3d& t, double& tout) {
  const Vec3d e1 = t.b - t.a;
  const Vec3d e2 = t.c - t.a;
  const Vec3d pv = geom::v1::cross(r.dir, e2);
  const double det = geom::v1::dot(e1, pv);
  const double scale = std::sqrt(geom::v1::dot(e1, e1)) *
                       std::sqrt(geom::v1::dot(e2, e2)) *
                       std::sqrt(geom::v1::dot(r.dir, r.dir));
  if (!(std::abs(det) > 1e-14 * scale)) return false;
  const Vec3d tv{r.origin[0] - t.a[0], r.origin[1] - t.a[1], r.origin[2] - t.a[2]};
  const double u = geom::v1::dot(tv, pv) / det;
  if (u < 0.0 || u > 1.0) return false;
  const Vec3d qv = geom::v1::cross(tv, e1);
  const double v = geom::v1::dot(r.dir, qv) / det;
  if (v < 0.0 || u + v > 1.0) return false;
  const double tt = geom::v1::dot(e2, qv) / det;
  if (tt <= 0.0) return false;
  tout = tt;
  return true;
}

}  // namespace

TEST_CASE("BVH raycast matches brute force") {
  std::mt19937 rng(42u);
  std::uniform_real_distribution<double> coord(-1.0, 1.0);
  std::vector<Triangle3d> tris;
  for (int i = 0; i < 200; ++i) {
    tris.push_back(Triangle3d{Point3d{coord(rng), coord(rng), coord(rng)},
                              Point3d{coord(rng), coord(rng), coord(rng)},
                              Point3d{coord(rng), coord(rng), coord(rng)}});
  }
  Bvh bvh;
  bvh.build(tris);
  REQUIRE_FALSE(bvh.bounds().is_empty());

  for (int r = 0; r < 500; ++r) {
    Ray3d ray{Point3d{coord(rng) * 3.0, coord(rng) * 3.0, 3.0},
              Vec3d{coord(rng) * 0.1, coord(rng) * 0.1, -1.0}};
    double best_brute = std::numeric_limits<double>::infinity();
    for (const auto& t : tris) {
      double tt = 0.0;
      if (ray_tri(ray, t, tt) && tt < best_brute) best_brute = tt;
    }
    Bvh::Hit hit;
    const bool got = bvh.raycast(ray, 1e30, hit);
    CHECK(got == std::isfinite(best_brute));
    if (got && std::isfinite(best_brute)) {
      CHECK(hit.t == doctest::Approx(best_brute).epsilon(1e-12));
    }
  }
}

TEST_CASE("BVH nearest matches brute force") {
  std::mt19937 rng(7u);
  std::uniform_real_distribution<double> coord(-1.0, 1.0);
  std::vector<Triangle3d> tris;
  for (int i = 0; i < 150; ++i) {
    tris.push_back(Triangle3d{Point3d{coord(rng), coord(rng), coord(rng)},
                              Point3d{coord(rng), coord(rng), coord(rng)},
                              Point3d{coord(rng), coord(rng), coord(rng)}});
  }
  Bvh bvh;
  bvh.build(tris);
  for (int q = 0; q < 200; ++q) {
    const Point3d p{coord(rng) * 2.0, coord(rng) * 2.0, coord(rng) * 2.0};
    double best = std::numeric_limits<double>::infinity();
    for (const auto& t : tris) {
      const auto d = geom::v1::squared_distance(p, t);
      best = std::min(best, d);
    }
    double d2 = std::numeric_limits<double>::infinity();
    const std::uint32_t idx = bvh.nearest(p, d2);
    CHECK(idx != 0xFFFFFFFFu);
    CHECK(d2 == doctest::Approx(best).epsilon(1e-10));
  }
}

// ---------------------------------------------------------------------------
// Appended cases: independent brute-force oracles, contract edges and ties.
// ---------------------------------------------------------------------------

namespace {

struct BruteHit {
  bool hit = false;
  std::uint32_t prim = 0;
  double t = 0.0;
};

// Nearest hit with the same tie rule as the BVH: smallest primitive index.
BruteHit brute_raycast(const std::vector<Triangle3d>& tris, const Ray3d& ray, double tmax) {
  BruteHit best;
  for (std::uint32_t i = 0; i < tris.size(); ++i) {
    double t = 0.0;
    if (!ray_tri(ray, tris[i], t) || t > tmax) continue;
    if (!best.hit || t < best.t || (t == best.t && i < best.prim)) {
      best.hit = true;
      best.prim = i;
      best.t = t;
    }
  }
  return best;
}

// Nearest primitive with the same tie rule: smallest primitive index.
std::uint32_t brute_nearest(const std::vector<Triangle3d>& tris, const Point3d& p, double& d2) {
  d2 = std::numeric_limits<double>::infinity();
  std::uint32_t best = 0xFFFFFFFFu;
  for (std::uint32_t i = 0; i < tris.size(); ++i) {
    const double d = geom::v1::squared_distance(p, tris[i]);
    if (d < d2 || (d == d2 && i < best)) {
      d2 = d;
      best = i;
    }
  }
  return best;
}

std::vector<Triangle3d> random_triangles(std::mt19937& rng, int count) {
  std::uniform_real_distribution<double> coord(-1.0, 1.0);
  std::vector<Triangle3d> tris;
  tris.reserve(static_cast<std::size_t>(count));
  for (int i = 0; i < count; ++i) {
    tris.push_back(Triangle3d{Point3d{coord(rng), coord(rng), coord(rng)},
                              Point3d{coord(rng), coord(rng), coord(rng)},
                              Point3d{coord(rng), coord(rng), coord(rng)}});
  }
  return tris;
}

bool same_box(const geom::v1::Box3d& lhs, const geom::v1::Box3d& rhs) {
  for (int i = 0; i < 3; ++i) {
    if (lhs.lo[i] != rhs.lo[i] || lhs.hi[i] != rhs.hi[i]) return false;
  }
  return true;
}

}  // namespace

TEST_CASE("BVH raycast matches an independent brute force exactly") {
  std::mt19937 rng(20261107u);
  const auto tris = random_triangles(rng, 300);
  Bvh bvh;
  bvh.build(tris);
  std::uniform_real_distribution<double> coord(-1.0, 1.0);
  int hits = 0;
  for (int r = 0; r < 1000; ++r) {
    const Point3d origin{coord(rng) * 4.0, coord(rng) * 4.0, 4.0};
    const Point3d target{coord(rng) * 1.5, coord(rng) * 1.5, coord(rng) * 1.5};
    const Ray3d ray{origin, target - origin};
    const BruteHit ref = brute_raycast(tris, ray, 1e30);
    Bvh::Hit hit;
    const bool got = bvh.raycast(ray, 1e30, hit);
    CHECK(got == ref.hit);
    if (got && ref.hit) {
      // Bit-identical t: both sides use the same arithmetic, and ties resolve
      // to the smallest primitive index on both sides.
      CHECK(hit.t == ref.t);
      CHECK(hit.prim == ref.prim);
      ++hits;
    }
  }
  CHECK(hits > 400);  // measured 539 of 1000 for this fixed seed
}

TEST_CASE("BVH nearest matches an independent brute force exactly") {
  std::mt19937 rng(20261108u);
  const auto tris = random_triangles(rng, 250);
  Bvh bvh;
  bvh.build(tris);
  std::uniform_real_distribution<double> coord(-1.0, 1.0);
  for (int q = 0; q < 500; ++q) {
    const Point3d p{coord(rng) * 3.0, coord(rng) * 3.0, coord(rng) * 3.0};
    double ref_d2 = 0.0;
    const std::uint32_t ref = brute_nearest(tris, p, ref_d2);
    double d2 = 0.0;
    const std::uint32_t idx = bvh.nearest(p, d2);
    CHECK(idx == ref);
    CHECK(d2 == ref_d2);
  }
}

TEST_CASE("BVH returns the root box as the union of the primitive boxes") {
  std::mt19937 rng(20261109u);
  const auto tris = random_triangles(rng, 120);
  Bvh bvh;
  bvh.build(tris);
  geom::v1::Box3d expected = geom::v1::Box3d::empty();
  for (const auto& t : tris) expected.expand(geom::v1::triangle_box(t));
  CHECK_FALSE(bvh.bounds().is_empty());
  CHECK(same_box(bvh.bounds(), expected));
  CHECK(bvh.size() == tris.size());

  Bvh empty;
  empty.build({});
  CHECK(empty.bounds().is_empty());
  CHECK(empty.size() == 0u);
  Bvh::Hit hit;
  CHECK_FALSE(empty.raycast(Ray3d{Point3d{0, 0, 0}, Vec3d{0, 0, -1}}, 1e30, hit));
  double d2 = 123.0;
  CHECK(empty.nearest(Point3d{0, 0, 0}, d2) == 0xFFFFFFFFu);
  CHECK(std::isinf(d2));
}

TEST_CASE("BVH handles axis-parallel rays without NaN") {
  std::mt19937 rng(20261110u);
  const auto tris = random_triangles(rng, 60);
  Bvh bvh;
  bvh.build(tris);
  std::uniform_real_distribution<double> coord(-2.0, 2.0);
  const Vec3d dirs[6] = {{1, 0, 0},  {0, 1, 0},  {0, 0, 1},
                         {-1, 0, 0}, {0, -1, 0}, {0, 0, -1}};
  for (const Vec3d& dir : dirs) {
    for (int r = 0; r < 100; ++r) {
      const Ray3d ray{Point3d{coord(rng), coord(rng), coord(rng)}, dir};
      const BruteHit ref = brute_raycast(tris, ray, 1e30);
      Bvh::Hit hit;
      const bool got = bvh.raycast(ray, 1e30, hit);
      CHECK(got == ref.hit);
      if (got && ref.hit) {
        CHECK(std::isfinite(hit.t));
        CHECK(hit.t == ref.t);
        CHECK(hit.prim == ref.prim);
      }
    }
  }
}

TEST_CASE("BVH ray parallel to a slab plane at an exact boundary does not produce NaN") {
  // Origin coordinates equal a box plane along the two zero direction axes:
  // a naive (plane - origin) * (1 / 0) evaluates 0 * inf = NaN.
  const Triangle3d tri{Point3d{0, 0, 0}, Point3d{1, 0, 0}, Point3d{0, 1, 0}};
  std::vector<Triangle3d> tris{tri};
  Bvh bvh;
  bvh.build(tris);
  const Ray3d ray{Point3d{-1.0, 0.0, 0.0}, Vec3d{1.0, 0.0, 0.0}};
  Bvh::Hit hit;
  const bool got = bvh.raycast(ray, 1e30, hit);
  CHECK(got == brute_raycast(tris, ray, 1e30).hit);
  if (got) CHECK(std::isfinite(hit.t));
}

TEST_CASE("BVH ray origin inside the root box still finds forward hits") {
  const Triangle3d tri{Point3d{-2, -2, -1}, Point3d{2, -2, -1}, Point3d{0, 2, 2}};
  std::vector<Triangle3d> tris{tri};
  Bvh bvh;
  bvh.build(tris);
  REQUIRE_FALSE(bvh.bounds().is_empty());
  const Point3d inside{0.0, 0.0, 0.0};  // strictly inside the root box
  CHECK(bvh.bounds().contains(inside));
  const Ray3d ray{inside, Vec3d{0, 0, 1}};
  const BruteHit ref = brute_raycast(tris, ray, 1e30);
  REQUIRE(ref.hit);
  Bvh::Hit hit;
  CHECK(bvh.raycast(ray, 1e30, hit));
  CHECK(hit.prim == 0u);
  CHECK(hit.t == ref.t);
}

TEST_CASE("BVH handles zero-length boxes, coincident boxes and tie-breaks") {
  // Primitive 0 is a degenerate point triangle: a zero-length bounding box.
  const Triangle3d point{Point3d{0, 0, 0}, Point3d{0, 0, 0}, Point3d{0, 0, 0}};
  const Triangle3d real{Point3d{0, 0, 1}, Point3d{1, 0, 1}, Point3d{0, 1, 1}};
  std::vector<Triangle3d> tris{point, real};
  Bvh bvh;
  bvh.build(tris);
  CHECK_FALSE(bvh.bounds().is_empty());
  const Ray3d ray{Point3d{0.25, 0.25, 5.0}, Vec3d{0, 0, -1}};
  Bvh::Hit hit;
  CHECK(bvh.raycast(ray, 1e30, hit));
  CHECK(hit.prim == 1u);
  double d2 = 0.0;
  CHECK(bvh.nearest(Point3d{0.0, 0.0, 0.0}, d2) == 0u);
  CHECK(d2 == 0.0);

  // Three coincident triangles: every hit distance ties, so the smallest
  // primitive index must win for both queries.
  std::vector<Triangle3d> coins{real, real, real};
  Bvh tie_bvh;
  tie_bvh.build(coins);
  Bvh::Hit tie_hit;
  CHECK(tie_bvh.raycast(ray, 1e30, tie_hit));
  CHECK(tie_hit.prim == 0u);
  double tie_d2 = 0.0;
  CHECK(tie_bvh.nearest(Point3d{0.25, 0.25, 3.0}, tie_d2) == 0u);
}

TEST_CASE("BVH nearest breaks equal distances by the smallest index") {
  // Mirrored triangles give exactly equal squared distances from (0.2, 0.2, 0).
  const Triangle3d upper{Point3d{0, 0, 1}, Point3d{1, 0, 1}, Point3d{0, 1, 1}};
  const Triangle3d lower{Point3d{0, 0, -1}, Point3d{1, 0, -1}, Point3d{0, 1, -1}};
  std::vector<Triangle3d> tris{upper, lower};
  Bvh bvh;
  bvh.build(tris);
  const Point3d q{0.2, 0.2, 0.0};
  CHECK(geom::v1::squared_distance(q, upper) == geom::v1::squared_distance(q, lower));
  double d2 = 0.0;
  CHECK(bvh.nearest(q, d2) == 0u);
  CHECK(d2 == geom::v1::squared_distance(q, upper));

  // Reversed input order: the winner is the new index 0 (formerly lower).
  std::vector<Triangle3d> reversed{lower, upper};
  Bvh reversed_bvh;
  reversed_bvh.build(reversed);
  double rd2 = 0.0;
  CHECK(reversed_bvh.nearest(q, rd2) == 0u);
}

TEST_CASE("BVH single primitive is found by both queries") {
  const Triangle3d tri{Point3d{0, 0, 0}, Point3d{1, 0, 0}, Point3d{0, 1, 0}};
  std::vector<Triangle3d> tris{tri};
  Bvh bvh;
  bvh.build(tris);
  CHECK(bvh.size() == 1u);
  CHECK(same_box(bvh.bounds(), geom::v1::triangle_box(tri)));
  Bvh::Hit hit;
  CHECK(bvh.raycast(Ray3d{Point3d{0.25, 0.25, 1.0}, Vec3d{0, 0, -1}}, 1e30, hit));
  CHECK(hit.prim == 0u);
  CHECK(hit.t == doctest::Approx(1.0));
  double d2 = 0.0;
  CHECK(bvh.nearest(Point3d{0.25, 0.25, 0.0}, d2) == 0u);
  CHECK(d2 == 0.0);
}

TEST_CASE("BVH builds many coincident-centroid primitives without degenerating") {
  // Every primitive has the same centroid box, so the 12-bin SAH has no axis
  // extent and must fall back to the deterministic median split. The build
  // must terminate and the queries must still resolve to the smallest index.
  const Triangle3d tri{Point3d{0, 0, 1}, Point3d{1, 0, 1}, Point3d{0, 1, 1}};
  std::vector<Triangle3d> tris(500, tri);
  Bvh bvh;
  bvh.build(tris);
  CHECK(bvh.size() == 500u);
  CHECK_FALSE(bvh.bounds().is_empty());
  Bvh::Hit hit;
  CHECK(bvh.raycast(Ray3d{Point3d{0.25, 0.25, 5.0}, Vec3d{0, 0, -1}}, 1e30, hit));
  CHECK(hit.prim == 0u);
  double d2 = 0.0;
  CHECK(bvh.nearest(Point3d{0.25, 0.25, 2.0}, d2) == 0u);
  CHECK(d2 == 1.0);
}
// ---------------------------------------------------------------------------
// Fix round: degenerate distance, the inclusive tmax boundary, and tiny
// triangles. Expected values are hand-computed, independent of the brief's
// ray_tri helper.
// ---------------------------------------------------------------------------

TEST_CASE("BVH nearest handles degenerate triangles without NaN") {
  // a == b: no plane exists, so the distance is the minimum over the edges.
  const Triangle3d collapsed{Point3d{0, 0, 0}, Point3d{0, 0, 0}, Point3d{1, 0, 0}};
  CHECK(squared_distance(Point3d{0.5, 0, 0}, collapsed) == 0.0);
  CHECK(squared_distance(Point3d{0.5, 1, 0}, collapsed) == 1.0);
  CHECK(squared_distance(Point3d{-1, 0, 0}, collapsed) == 1.0);

  // Exactly collinear vertices: the point is 5 away from the segment on x.
  const Triangle3d collinear{Point3d{0, 0, 0}, Point3d{2, 0, 0}, Point3d{1, 0, 0}};
  CHECK(squared_distance(Point3d{1, 5, 0}, collinear) == 25.0);
  CHECK(squared_distance(Point3d{0.5, 0, 3}, collinear) == 9.0);

  // Fully collapsed triangle: a point.
  const Triangle3d point{Point3d{0, 0, 0}, Point3d{0, 0, 0}, Point3d{0, 0, 0}};
  CHECK(squared_distance(Point3d{1, 2, 2}, point) == 9.0);

  std::vector<Triangle3d> tris{collapsed};
  Bvh bvh;
  bvh.build(tris);
  double d2 = -1.0;
  CHECK(bvh.nearest(Point3d{0.5, 0, 0}, d2) == 0u);
  CHECK(d2 == 0.0);
  CHECK(bvh.nearest(Point3d{0.5, 1, 0}, d2) == 0u);
  CHECK(d2 == 1.0);
}

TEST_CASE("BVH raycast includes a hit exactly at tmax") {
  const Triangle3d tri{Point3d{0, 0, 0}, Point3d{1, 0, 0}, Point3d{0, 1, 0}};
  std::vector<Triangle3d> tris{tri};
  Bvh bvh;
  bvh.build(tris);
  const Ray3d ray{Point3d{0.25, 0.25, 1.0}, Vec3d{0, 0, -1}};

  Bvh::Hit hit;
  CHECK(bvh.raycast(ray, 1.0, hit));
  CHECK(hit.prim == 0u);
  CHECK(hit.t == 1.0);

  // One ulp below the hit distance is a miss, and so is any smaller budget.
  CHECK_FALSE(bvh.raycast(ray, std::nextafter(1.0, 0.0), hit));
  CHECK_FALSE(bvh.raycast(ray, 0.5, hit));

  // A budget just above the hit still reports the same t.
  CHECK(bvh.raycast(ray, std::nextafter(1.0, 2.0), hit));
  CHECK(hit.t == 1.0);
}

TEST_CASE("BVH raycast keeps an exact t = 1 hit exact") {
  // Regression: tt = num * (1/det) rounded num == det to 0x1.fffffffffffffp-1.
  // The six literals are exact hex floats; the hit point is constructed at
  // t = 1 by origin = p - dir.
  const Triangle3d tri{
      Point3d{-0x1.97516e455b738p+0, 0x1.6edbddb7820bap+0, -0x1.4b6b2cd729280p-3},
      Point3d{0x1.d5400339b7474p+0, 0x1.db1788b645f04p-1, 0x1.30a65e6874de0p-2},
      Point3d{-0x1.2ab1c98db2b68p-2, 0x1.7d3114c3d6edcp-1, -0x1.07ce2d175799ap+0}};
  std::vector<Triangle3d> tris{tri};
  Bvh bvh;
  bvh.build(tris);
  const Ray3d ray{
      Point3d{-0x1.cbc66c5d49575p-1, 0x1.f5244fde8fbc8p-3, 0x1.41a2b266eaa0ep-1},
      Vec3d{0x1.4634d534141f2p-1, 0x1.cca0934339b4ep-1, -0x1.9c7d031b3f67ep-1}};

  Bvh::Hit hit;
  REQUIRE(bvh.raycast(ray, 2.0, hit));
  CHECK(hit.prim == 0u);
  CHECK(hit.t == 1.0);
}

TEST_CASE("BVH raycast and nearest are scale-invariant for tiny triangles") {
  for (const double s : {1e-8, 1e-7, 1e-6}) {
    CAPTURE(s);
    const Triangle3d tri{Point3d{0, 0, 0}, Point3d{s, 0, 0}, Point3d{0, s, 0}};
    std::vector<Triangle3d> tris{tri};
    Bvh bvh;
    bvh.build(tris);

    const Ray3d ray{Point3d{0.25 * s, 0.25 * s, 1.0}, Vec3d{0, 0, -1}};
    Bvh::Hit hit;
    CHECK(bvh.raycast(ray, 2.0, hit));
    CHECK(hit.prim == 0u);
    CHECK(hit.t == doctest::Approx(1.0));

    double d2 = -1.0;
    CHECK(bvh.nearest(Point3d{0.25 * s, 0.25 * s, 1.0}, d2) == 0u);
    CHECK(d2 == 1.0);
    CHECK(squared_distance(Point3d{s, s, 1.0}, tri) == doctest::Approx(1.0));
  }
}