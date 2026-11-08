#include <doctest/doctest.h>

#include <cmath>
#include <string_view>
#include <vector>

#include <geom/cgeom/delaunay2.hpp>
#include <geom/cgeom/hull2d.hpp>
#include <geom/core/handle.hpp>
#include <geom/curves/nurbs_curve.hpp>
#include <geom/index/bvh.hpp>
#include <geom/numeric/expansion.hpp>
#include <geom/transform/affine.hpp>

using geom::v1::Affine;
using geom::v1::Bvh;
using geom::v1::convex_hull2d;
using geom::v1::delaunay2;
using geom::v1::detail::Expansion;
using geom::v1::Errc;
using geom::v1::Handle;
using geom::v1::HandleAllocator;
using geom::v1::NurbsCurve;
using geom::v1::Point2d;
using geom::v1::Point3d;
using geom::v1::Ray3d;
using geom::v1::Vec3d;

TEST_CASE("singular affine inverse reports degenerate") {
  Affine a{};
  for (int r = 0; r < 3; ++r) {
    for (int c = 0; c < 3; ++c) a.m[r][c] = 0.0;
  }
  const auto inv = a.inverse();
  REQUIRE_FALSE(inv.has_value());
  CHECK(inv.error().code == Errc::degenerate);
}

TEST_CASE("NURBS rejects zero weights, bad knots and size mismatches") {
  CHECK_FALSE(NurbsCurve(1, {0.0, 0.0, 1.0, 1.0},
                         {Point3d{0, 0, 0}, Point3d{1, 0, 0}}, {1.0, 0.0}).is_valid());
  CHECK_FALSE(NurbsCurve(1, {0.0, 1.0, 0.0, 1.0},
                         {Point3d{0, 0, 0}, Point3d{1, 0, 0}}, {1.0, 1.0}).is_valid());
  CHECK_FALSE(NurbsCurve(1, {0.0, 0.0, 1.0, 1.0},
                         {Point3d{0, 0, 0}, Point3d{1, 0, 0}}, {1.0}).is_valid());
}

TEST_CASE("NURBS clamps evaluation outside the domain") {
  const NurbsCurve c(1, {0.0, 0.0, 1.0, 1.0},
                     {Point3d{0, 0, 0}, Point3d{2, 0, 0}}, {1.0, 1.0});
  REQUIRE(c.is_valid());
  CHECK(c.eval(-5.0)[0] == doctest::Approx(0.0));
  CHECK(c.eval(5.0)[0] == doctest::Approx(2.0));
}

TEST_CASE("empty BVH answers queries safely") {
  Bvh b;
  Bvh::Hit hit;
  CHECK_FALSE(b.raycast(Ray3d{Point3d{0, 0, 0}, Vec3d{0, 0, 1}}, 1e9, hit));
  double d2 = -1.0;
  CHECK(b.nearest(Point3d{0, 0, 0}, d2) == 0xFFFFFFFFu);
  CHECK(std::isinf(d2));
  CHECK(b.bounds().is_empty());
}

TEST_CASE("hull collapses identical and two-point input") {
  const std::vector<Point2d> same{{1, 1}, {1, 1}, {1, 1}};
  CHECK(convex_hull2d(same).size() == 1);

  const std::vector<Point2d> two{{0, 0}, {2, 2}};
  CHECK(convex_hull2d(two).size() == 2);
}

TEST_CASE("delaunay of a single triangle yields one triangle") {
  const std::vector<Point2d> tri{{0, 0}, {1, 0}, {0, 1}};
  const auto tris = delaunay2(tri);
  REQUIRE(tris.size() == 1);
}

TEST_CASE("expansion handles scalar scaling and null handles") {
  Expansion e;
  e.add_scaled(Expansion(1.0), 2.0);
  CHECK(e.approx() == doctest::Approx(2.0));

  HandleAllocator alloc;
  alloc.release(Handle::null());
  CHECK_FALSE(alloc.alive(Handle::null()));
}

TEST_CASE("self_intersecting error code has a name") {
  CHECK(std::string_view(geom::v1::to_string(Errc::self_intersecting)) ==
        "self_intersecting");
}

TEST_CASE("brep error codes have names") {
  CHECK(std::string_view(geom::v1::to_string(Errc::invalid_topology)) == "invalid_topology");
  CHECK(std::string_view(geom::v1::to_string(Errc::tolerance_violation)) ==
        "tolerance_violation");
  CHECK(std::string_view(geom::v1::to_string(Errc::not_closed)) == "not_closed");
}

TEST_CASE("allocator counts and enumerates live handles") {
  HandleAllocator alloc;
  CHECK(alloc.live_count() == 0);
  CHECK(alloc.handle_at(0).is_null());

  const Handle a = alloc.alloc();
  const Handle b = alloc.alloc();
  CHECK(alloc.live_count() == 2);
  CHECK(alloc.handle_at(0) == a);
  CHECK(alloc.handle_at(1) == b);
  CHECK(alloc.handle_at(2).is_null());

  alloc.release(a);
  CHECK(alloc.live_count() == 1);
  CHECK(alloc.handle_at(0).is_null());
  CHECK(alloc.handle_at(1) == b);

  const Handle reused = alloc.alloc();
  CHECK(reused.index() == a.index());
  CHECK(reused.generation() == static_cast<std::uint8_t>(a.generation() + 1));
  CHECK(alloc.live_count() == 2);
  CHECK(alloc.handle_at(a.index()) == reused);
}

