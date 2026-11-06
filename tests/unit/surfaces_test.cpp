#include <doctest/doctest.h>

#include <cmath>
#include <cstddef>
#include <limits>
#include <vector>

#include <geom/surfaces/nurbs_surface.hpp>
#include <geom/surfaces/plane.hpp>

using geom::v1::NurbsSurface;
using geom::v1::Plane3;
using geom::v1::Point3d;

TEST_CASE("bilinear patch evaluates to the corner average") {
  const NurbsSurface s(1, 1, {0.0, 0.0, 1.0, 1.0}, {0.0, 0.0, 1.0, 1.0},
                       {Point3d{0, 0, 0}, Point3d{1, 0, 1}, Point3d{0, 1, 2},
                        Point3d{1, 1, 3}},
                       {1.0, 1.0, 1.0, 1.0});
  REQUIRE(s.is_valid());
  const Point3d p = s.eval(0.5, 0.5);
  CHECK(p[0] == doctest::Approx(0.5));
  CHECK(p[1] == doctest::Approx(0.5));
  CHECK(p[2] == doctest::Approx(1.5));
}

TEST_CASE("surface derivatives match central differences") {
  const NurbsSurface s(1, 1, {0.0, 0.0, 1.0, 1.0}, {0.0, 0.0, 1.0, 1.0},
                       {Point3d{0, 0, 0}, Point3d{1, 0, 1}, Point3d{0, 1, 2},
                        Point3d{1, 1, 3}},
                       {1.0, 1.0, 1.0, 1.0});
  const double h = 1e-6;
  const auto fd_u = (s.eval(0.5 + h, 0.3) - s.eval(0.5 - h, 0.3)) * (1.0 / (2.0 * h));
  const auto exact_u = s.du(0.5, 0.3);
  CHECK(exact_u[0] == doctest::Approx(fd_u[0]).epsilon(1e-6));
  CHECK(exact_u[1] == doctest::Approx(fd_u[1]).epsilon(1e-6));
  CHECK(exact_u[2] == doctest::Approx(fd_u[2]).epsilon(1e-6));

  const auto fd_v = (s.eval(0.2, 0.5 + h) - s.eval(0.2, 0.5 - h)) * (1.0 / (2.0 * h));
  const auto exact_v = s.dv(0.2, 0.5);
  CHECK(exact_v[0] == doctest::Approx(fd_v[0]).epsilon(1e-6));
  CHECK(exact_v[1] == doctest::Approx(fd_v[1]).epsilon(1e-6));
  CHECK(exact_v[2] == doctest::Approx(fd_v[2]).epsilon(1e-6));
}

TEST_CASE("plane distance and side are consistent") {
  const Plane3 plane{geom::v1::Vec3d{0, 0, 1}, 5.0};
  const Point3d p{1, 2, 7};
  CHECK(plane.signed_distance(p) == doctest::Approx(2.0));
  CHECK(plane.side(p) == 1);
  CHECK(plane.project(p)[2] == doctest::Approx(5.0));
  CHECK(plane.side(Point3d{1, 2, 5}) == 0);
}

TEST_CASE("bilinear surface has constant analytic derivatives") {
  // Control points are indexed i * nv + j (i along u), so this patch is
  // S(u,v) = (v, u, v + 2u): dS/du = (0, 1, 2) and dS/dv = (1, 0, 1).
  const NurbsSurface s(1, 1, {0.0, 0.0, 1.0, 1.0}, {0.0, 0.0, 1.0, 1.0},
                       {Point3d{0, 0, 0}, Point3d{1, 0, 1}, Point3d{0, 1, 2},
                        Point3d{1, 1, 3}},
                       {1.0, 1.0, 1.0, 1.0});
  REQUIRE(s.is_valid());
  const auto du = s.du(0.3, 0.7);
  CHECK(du[0] == doctest::Approx(0.0).epsilon(1e-13));
  CHECK(du[1] == doctest::Approx(1.0).epsilon(1e-13));
  CHECK(du[2] == doctest::Approx(2.0).epsilon(1e-13));
  const auto dv = s.dv(0.3, 0.7);
  CHECK(dv[0] == doctest::Approx(1.0).epsilon(1e-13));
  CHECK(dv[1] == doctest::Approx(0.0).epsilon(1e-13));
  CHECK(dv[2] == doctest::Approx(1.0).epsilon(1e-13));
}

TEST_CASE("rational surface derivatives divide by the evaluated weight") {
  // Non-uniform weights: differentiating weighted control points without the
  // quotient rule gives different values for this patch.
  const NurbsSurface s(1, 1, {0.0, 0.0, 1.0, 1.0}, {0.0, 0.0, 1.0, 1.0},
                       {Point3d{0, 0, 0}, Point3d{1, 0, 1}, Point3d{0, 1, 2},
                        Point3d{1, 1, 3}},
                       {1.0, 2.0, 0.5, 3.0});
  REQUIRE(s.is_valid());
  const double h = 1e-6;
  const auto fd_u = (s.eval(0.3 + h, 0.7) - s.eval(0.3 - h, 0.7)) * (1.0 / (2.0 * h));
  const auto du = s.du(0.3, 0.7);
  CHECK(du[0] == doctest::Approx(fd_u[0]).epsilon(1e-5));
  CHECK(du[1] == doctest::Approx(fd_u[1]).epsilon(1e-5));
  CHECK(du[2] == doctest::Approx(fd_u[2]).epsilon(1e-5));

  const auto fd_v = (s.eval(0.6, 0.2 + h) - s.eval(0.6, 0.2 - h)) * (1.0 / (2.0 * h));
  const auto dv = s.dv(0.6, 0.2);
  CHECK(dv[0] == doctest::Approx(fd_v[0]).epsilon(1e-5));
  CHECK(dv[1] == doctest::Approx(fd_v[1]).epsilon(1e-5));
  CHECK(dv[2] == doctest::Approx(fd_v[2]).epsilon(1e-5));
}

TEST_CASE("NURBS surface clamps parameters to its domain") {
  const NurbsSurface s(1, 1, {0.0, 0.0, 1.0, 1.0}, {0.0, 0.0, 1.0, 1.0},
                       {Point3d{0, 0, 0}, Point3d{1, 0, 1}, Point3d{0, 1, 2},
                        Point3d{1, 1, 3}},
                       {1.0, 1.0, 1.0, 1.0});
  REQUIRE(s.is_valid());
  const auto [u0, u1] = s.domain_u();
  const auto [v0, v1] = s.domain_v();
  CHECK(u0 == 0.0);
  CHECK(u1 == 1.0);
  CHECK(v0 == 0.0);
  CHECK(v1 == 1.0);
  CHECK(s.eval(-5.0, 7.0).v[0] == s.eval(u0, v1).v[0]);
  CHECK(s.eval(-5.0, 7.0).v[2] == s.eval(u0, v1).v[2]);
  CHECK(s.du(-5.0, 7.0).v[0] == s.du(u0, v1).v[0]);
  CHECK(s.dv(-5.0, 7.0).v[1] == s.dv(u0, v1).v[1]);
}

TEST_CASE("NURBS surface validity rejects non-finite and degenerate data") {
  const std::vector<Point3d> pts{Point3d{0, 0, 0}, Point3d{1, 0, 1}, Point3d{0, 1, 2},
                                 Point3d{1, 1, 3}};
  const std::vector<double> w{1.0, 1.0, 1.0, 1.0};
  const std::vector<double> ku{0.0, 0.0, 1.0, 1.0};
  const std::vector<double> kv{0.0, 0.0, 1.0, 1.0};
  const double inf = std::numeric_limits<double>::infinity();
  const double nan = std::numeric_limits<double>::quiet_NaN();

  CHECK_FALSE(NurbsSurface(1, 1, {0.0, 0.0, 1.0, inf}, kv, pts, w).is_valid());
  CHECK_FALSE(NurbsSurface(1, 1, ku, {0.0, 0.0, nan, 1.0}, pts, w).is_valid());
  CHECK_FALSE(NurbsSurface(1, 1, ku, kv, pts, {1.0, 1.0, inf, 1.0}).is_valid());
  CHECK_FALSE(NurbsSurface(1, 1, ku, kv, pts, {1.0, 1.0, nan, 1.0}).is_valid());
  CHECK_FALSE(NurbsSurface(1, 1, ku, kv,
                           {Point3d{0, 0, 0}, Point3d{inf, 0, 1}, Point3d{0, 1, 2},
                            Point3d{1, 1, 3}},
                           w)
                  .is_valid());
  CHECK_FALSE(NurbsSurface(1, 1, ku, kv,
                           {Point3d{0, 0, 0}, Point3d{nan, 0, 1}, Point3d{0, 1, 2},
                            Point3d{1, 1, 3}},
                           w)
                  .is_valid());

  // Degenerate parameter domains in either direction.
  const NurbsSurface flat_u(1, 1, {0.0, 0.0, 0.0, 0.0}, kv, pts, w);
  CHECK_FALSE(flat_u.is_valid());
  CHECK(flat_u.domain_u().first == 0.0);
  CHECK(flat_u.domain_u().second == 0.0);
  CHECK(flat_u.domain_v().first == 0.0);
  CHECK(flat_u.domain_v().second == 0.0);
  CHECK_FALSE(NurbsSurface(1, 1, ku, {0.0, 0.0, 0.0, 0.0}, pts, w).is_valid());

  // Interior multiplicity p + 2 in u exceeds the maximal p + 1.
  const std::vector<double> ku_multi{0.0, 0.0, 0.5, 0.5, 0.5, 1.0};
  CHECK_FALSE(NurbsSurface(1, 1, ku_multi, kv,
                           std::vector<Point3d>(8, Point3d{0.0, 0.0, 0.0}),
                           std::vector<double>(8, 1.0))
                  .is_valid());

  // Degrees above the de Boor scratch bound are rejected; degree 16 is accepted.
  std::vector<double> ku17(36), kv17(4);
  for (std::size_t i = 0; i < ku17.size(); ++i) {
    ku17[i] = static_cast<double>(i) / static_cast<double>(ku17.size() - 1);
  }
  kv17 = {0.0, 0.0, 1.0, 1.0};
  const NurbsSurface too_high(17, 1, ku17, kv17,
                              std::vector<Point3d>(36, Point3d{0.0, 0.0, 0.0}),
                              std::vector<double>(36, 1.0));
  CHECK_FALSE(too_high.is_valid());

  std::vector<double> ku16(34);
  for (std::size_t i = 0; i < ku16.size(); ++i) {
    ku16[i] = static_cast<double>(i) / static_cast<double>(ku16.size() - 1);
  }
  const NurbsSurface max_degree(16, 1, ku16, kv17,
                                std::vector<Point3d>(34, Point3d{0.0, 0.0, 0.0}),
                                std::vector<double>(34, 1.0));
  CHECK(max_degree.is_valid());
  CHECK(max_degree.eval(0.5, 0.5).v[0] == 0.0);
}

TEST_CASE("plane projection and side for tilted and on-plane points") {
  const Plane3 tilted{geom::v1::Vec3d{0, 1, 0}, -2.0};
  const Point3d above{3, 4, 5};
  CHECK(tilted.signed_distance(above) == doctest::Approx(6.0));
  CHECK(tilted.side(above) == 1);
  const Point3d below{3, -5, 5};
  CHECK(tilted.signed_distance(below) == doctest::Approx(-3.0));
  CHECK(tilted.side(below) == -1);
  const Point3d projected = tilted.project(above);
  CHECK(projected.v[0] == doctest::Approx(3.0));
  CHECK(projected.v[1] == doctest::Approx(-2.0));
  CHECK(projected.v[2] == doctest::Approx(5.0));
  CHECK(tilted.signed_distance(projected) == doctest::Approx(0.0).epsilon(1e-15));

  const Point3d on_plane{7, -2, 1};
  CHECK(tilted.signed_distance(on_plane) == doctest::Approx(0.0));
  CHECK(tilted.side(on_plane) == 0);
  CHECK(tilted.project(on_plane).v[1] == doctest::Approx(-2.0));
}

TEST_CASE("plane queries are normal-scale invariant") {
  // Same geometric plane as 2y = 4, expressed with three different normals.
  const Point3d p{0, 4, 0};
  const Plane3 scaled{geom::v1::Vec3d{0, 2, 0}, 4.0};
  CHECK(scaled.signed_distance(p) == doctest::Approx(2.0));
  CHECK(scaled.side(p) == 1);
  const Point3d projected = scaled.project(p);
  CHECK(projected.v[0] == doctest::Approx(0.0));
  CHECK(projected.v[1] == doctest::Approx(2.0));
  CHECK(projected.v[2] == doctest::Approx(0.0));
  CHECK(scaled.signed_distance(projected) == doctest::Approx(0.0).epsilon(1e-15));

  const Plane3 shrunk{geom::v1::Vec3d{0, 0.5, 0}, 1.0};
  CHECK(shrunk.signed_distance(p) == doctest::Approx(2.0));
  CHECK(shrunk.project(p).v[1] == doctest::Approx(2.0));
  CHECK(shrunk.side(p) == 1);

  // Flipping the normal flips the orientation but not the projected point.
  const Plane3 flipped{geom::v1::Vec3d{0, -2, 0}, 4.0};
  CHECK(flipped.signed_distance(p) == doctest::Approx(-6.0));
  CHECK(flipped.side(p) == -1);
  CHECK(flipped.project(p).v[1] == doctest::Approx(-2.0));
  CHECK(flipped.signed_distance(flipped.project(p)) == doctest::Approx(0.0).epsilon(1e-15));

  // On-plane query with a non-unit, non-axis-aligned normal: 3x = 6, i.e. x = 2.
  const Plane3 oblique{geom::v1::Vec3d{3, 0, 0}, 6.0};
  const Point3d on{2, 5, 7};
  CHECK(oblique.signed_distance(on) == doctest::Approx(0.0));
  CHECK(oblique.side(on) == 0);
  const Point3d kept = oblique.project(on);
  CHECK(kept.v[0] == doctest::Approx(2.0));
  CHECK(kept.v[1] == doctest::Approx(5.0));
  CHECK(kept.v[2] == doctest::Approx(7.0));
  CHECK(oblique.signed_distance(Point3d{5, 0, 0}) == doctest::Approx(3.0));
  CHECK(oblique.project(Point3d{5, 0, 0}).v[0] == doctest::Approx(2.0));
}
