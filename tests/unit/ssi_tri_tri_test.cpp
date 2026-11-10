#include <doctest/doctest.h>

#include <cmath>
#include <geom/brep/detail/ssi_tri_tri.hpp>
#include <limits>

using namespace geom::v1;

TEST_CASE("ssi: crossing triangles yield one segment") {
  // Deviation from the brief: its q1 (0,-1,1) puts Q's z=0 chord at y in
  // [-1,0], so the closed triangles meet only at (0,0,0). Moving q1 to
  // (0,1,-1) restores the intended crossing segment (0,0,0)-(0,1,0).
  const TriTriResult r =
      tri_tri_intersect(Point3d{-1, 0, 0}, Point3d{1, 0, 0}, Point3d{0, 1, 0}, Point3d{0, -1, -1},
                        Point3d{0, 1, -1}, Point3d{0, 1, 1}, 1e-12);
  REQUIRE(r.kind == TriTriKind::segment);
  CHECK(std::fabs(r.a[0]) < 1e-12);
  CHECK(std::fabs(r.b[0]) < 1e-12);
  CHECK(std::fabs(r.a[1] - r.b[1]) > 1e-6);
}

TEST_CASE("ssi: planes crossing outside the triangles is a single-point contact") {
  // The brief's original Q configuration: the planes cross, but the closed
  // triangles meet only at the shared vertex (0,0,0).
  const TriTriResult r =
      tri_tri_intersect(Point3d{-1, 0, 0}, Point3d{1, 0, 0}, Point3d{0, 1, 0}, Point3d{0, -1, -1},
                        Point3d{0, -1, 1}, Point3d{0, 1, 1}, 1e-12);
  REQUIRE(r.kind == TriTriKind::touching);
  CHECK(r.a[0] == 0.0);
  CHECK(r.a[1] == 0.0);
  CHECK(r.a[2] == 0.0);
  CHECK(r.b[0] == 0.0);
  CHECK(r.b[1] == 0.0);
  CHECK(r.b[2] == 0.0);
}

TEST_CASE("ssi: separated and vertex-touching triangles are classified") {
  const TriTriResult separated =
      tri_tri_intersect(Point3d{0, 0, 0}, Point3d{1, 0, 0}, Point3d{0, 1, 0}, Point3d{0, 0, 1},
                        Point3d{1, 0, 1}, Point3d{0, 1, 1}, 1e-12);
  CHECK(separated.kind == TriTriKind::none);
  const TriTriResult touch =
      tri_tri_intersect(Point3d{0, 0, 0}, Point3d{1, 0, 0}, Point3d{0, 1, 0}, Point3d{0, 0, 0},
                        Point3d{0, -1, 1}, Point3d{0, 1, 1}, 1e-12);
  CHECK(touch.kind == TriTriKind::touching);
}

TEST_CASE("ssi: coplanar triangles are reported, not silently intersected") {
  const TriTriResult r =
      tri_tri_intersect(Point3d{0, 0, 0}, Point3d{2, 0, 0}, Point3d{0, 2, 0}, Point3d{0.5, 0.5, 0},
                        Point3d{2.5, 0.5, 0}, Point3d{0.5, 2.5, 0}, 1e-12);
  CHECK(r.kind == TriTriKind::coplanar);
}

TEST_CASE("ssi: degenerate and non-finite triangles are touching (no segment)") {
  const TriTriResult collapsed =
      tri_tri_intersect(Point3d{0, 0, 0}, Point3d{0, 0, 0}, Point3d{0, 0, 0}, Point3d{-1, -1, 0},
                        Point3d{1, -1, 0}, Point3d{0, 1, 0}, 1e-12);
  CHECK(collapsed.kind == TriTriKind::touching);
  CHECK_FALSE(tri_finite(Point3d{0, 0, std::nan("")}, Point3d{1, 0, 0}, Point3d{0, 1, 0}));
}

TEST_CASE("ssi: segment endpoints are deterministic and lexicographically ordered") {
  const auto intersect = [] {
    return tri_tri_intersect(Point3d{-1, 0, 0}, Point3d{1, 0, 0}, Point3d{0, 1, 0},
                             Point3d{0, -1, -1}, Point3d{0, 1, -1}, Point3d{0, 1, 1}, 1e-12);
  };
  const TriTriResult first = intersect();
  const TriTriResult second = intersect();
  REQUIRE(first.kind == TriTriKind::segment);
  REQUIRE(second.kind == TriTriKind::segment);
  CHECK(first.a[0] == second.a[0]);
  CHECK(first.a[1] == second.a[1]);
  CHECK(first.a[2] == second.a[2]);
  CHECK(first.b[0] == second.b[0]);
  CHECK(first.b[1] == second.b[1]);
  CHECK(first.b[2] == second.b[2]);
  CHECK(first.a[0] == 0.0);
  CHECK(first.a[1] == 0.0);
  CHECK(first.a[2] == 0.0);
  CHECK(first.b[0] == 0.0);
  CHECK(first.b[1] == 1.0);
  CHECK(first.b[2] == 0.0);
  const bool ordered =
      first.a[0] < first.b[0] ||
      (first.a[0] == first.b[0] &&
       (first.a[1] < first.b[1] || (first.a[1] == first.b[1] && first.a[2] <= first.b[2])));
  CHECK(ordered);
}

TEST_CASE("ssi: disjoint coplanar triangles are still coplanar") {
  const TriTriResult r =
      tri_tri_intersect(Point3d{0, 0, 0}, Point3d{1, 0, 0}, Point3d{0, 1, 0}, Point3d{10, 10, 0},
                        Point3d{11, 10, 0}, Point3d{10, 11, 0}, 1e-12);
  CHECK(r.kind == TriTriKind::coplanar);
}

TEST_CASE("ssi: only finite positive tolerance is accepted") {
  const auto intersect = [](double tolerance) {
    return tri_tri_intersect(Point3d{-1, 0, 0}, Point3d{1, 0, 0}, Point3d{0, 1, 0},
                             Point3d{0, -1, -1}, Point3d{0, 1, -1}, Point3d{0, 1, 1}, tolerance);
  };
  CHECK(intersect(0.0).kind == TriTriKind::touching);
  CHECK(intersect(-0.0).kind == TriTriKind::touching);
  CHECK(intersect(-1e-12).kind == TriTriKind::touching);
  CHECK(intersect(std::nan("")).kind == TriTriKind::touching);
  CHECK(intersect(std::numeric_limits<double>::infinity()).kind == TriTriKind::touching);
}

TEST_CASE("ssi: a non-finite vertex makes the intersection touching") {
  const TriTriResult r = tri_tri_intersect(Point3d{std::numeric_limits<double>::infinity(), 0, 0},
                                           Point3d{1, 0, 0}, Point3d{0, 1, 0}, Point3d{0, 0, 0},
                                           Point3d{0, -1, 1}, Point3d{0, 1, 1}, 1e-12);
  CHECK(r.kind == TriTriKind::touching);
}

TEST_CASE("ssi: triangles sharing an edge yield that edge as a segment") {
  const TriTriResult r =
      tri_tri_intersect(Point3d{0, 0, 0}, Point3d{1, 0, 0}, Point3d{0, 1, 0}, Point3d{0, 0, 0},
                        Point3d{1, 0, 0}, Point3d{0, 1, 1}, 1e-12);
  REQUIRE(r.kind == TriTriKind::segment);
  CHECK(r.a[0] == 0.0);
  CHECK(r.a[1] == 0.0);
  CHECK(r.a[2] == 0.0);
  CHECK(r.b[0] == 1.0);
  CHECK(r.b[1] == 0.0);
  CHECK(r.b[2] == 0.0);
}

TEST_CASE("ssi: a collinear triangle is degenerate even when its edges are long") {
  const TriTriResult r =
      tri_tri_intersect(Point3d{0, 0, 0}, Point3d{0.5, 0, 0}, Point3d{1, 0, 0}, Point3d{-1, -1, 0},
                        Point3d{1, -1, 0}, Point3d{0, 1, 0}, 1e-12);
  CHECK(r.kind == TriTriKind::touching);
}

TEST_CASE("ssi: tolerance is the only authority for near-coplanar input") {
  const auto intersect = [](double tolerance) {
    return tri_tri_intersect(Point3d{0, 0, 0}, Point3d{1, 0, 0}, Point3d{0, 1, 0},
                             Point3d{0, 0, 1e-6}, Point3d{1, 0, 1e-6}, Point3d{0, 1, 1e-6},
                             tolerance);
  };
  CHECK(intersect(1e-9).kind == TriTriKind::none);
  CHECK(intersect(1e-3).kind == TriTriKind::coplanar);
}

TEST_CASE("ssi: detail point-in-triangle includes exactly the tolerance band") {
  const double tol = 1e-6;
  const Point3d t0{0, 0, 0};
  const Point3d t1{1, 0, 0};
  const Point3d t2{0, 1, 0};
  CHECK(detail::ssi_tri_point_in_triangle(Point3d{0.25, -0.5 * tol, 0}, t0, t1, t2, 2, tol));
  CHECK_FALSE(detail::ssi_tri_point_in_triangle(Point3d{0.25, -1.5 * tol, 0}, t0, t1, t2, 2, tol));
}

TEST_CASE("ssi: detail degeneracy accepts an area above tolerance squared") {
  const double tol = 1e-3;
  const Point3d p0{0, 0, 0};
  const Point3d p1{1, 0, 0};
  const Point3d p2{0.5, 2e-5, 0};
  const Vec3d normal = cross(p1 - p0, p2 - p0);
  // Area 1e-5 lies between tol^2 (1e-6) and tol (1e-3).
  CHECK_FALSE(detail::ssi_tri_degenerate(p0, p1, p2, normal, tol));
}

TEST_CASE("ssi: detail degeneracy rejects an edge shorter than tolerance") {
  const double tol = 1e-3;
  const Point3d p0{0, 0, 0};
  const Point3d p1{5e-4, 0, 0};
  const Point3d p2{0, 1, 0};
  const Vec3d normal = cross(p1 - p0, p2 - p0);
  // Edge 5e-4 is at most tol, while the area 2.5e-4 exceeds tol^2.
  CHECK(detail::ssi_tri_degenerate(p0, p1, p2, normal, tol));
}

TEST_CASE("ssi: coplanar is detected from either triangle's plane") {
  // Tiny P lies within tol of large tilted Q's plane; Q's vertices are far
  // from P's plane (z = 0).
  const TriTriResult tiny_near_q =
      tri_tri_intersect(Point3d{0, 0, 0}, Point3d{1e-3, 0, 0}, Point3d{0, 1e-3, 0},
                        Point3d{0, 0, 0}, Point3d{1000, 0, 1}, Point3d{0, 1000, 0}, 1e-6);
  CHECK(tiny_near_q.kind == TriTriKind::coplanar);

  // Mirror case: tiny Q lies within tol of large tilted P's plane.
  const TriTriResult tiny_near_p =
      tri_tri_intersect(Point3d{0, 0, 0}, Point3d{1000, 0, 1}, Point3d{0, 1000, 0},
                        Point3d{0, 0, 0}, Point3d{1e-3, 0, 0}, Point3d{0, 1e-3, 0}, 1e-6);
  CHECK(tiny_near_p.kind == TriTriKind::coplanar);
}

TEST_CASE("ssi: large coordinates still classify") {
  const double s = 1e200;
  const TriTriResult crossing =
      tri_tri_intersect(Point3d{-s, 0, 0}, Point3d{s, 0, 0}, Point3d{0, s, 0}, Point3d{0, -s, -s},
                        Point3d{0, s, -s}, Point3d{0, s, s}, 1e-12);
  REQUIRE(crossing.kind == TriTriKind::segment);
  CHECK(crossing.a[0] == 0.0);
  CHECK(crossing.a[1] == 0.0);
  CHECK(crossing.a[2] == 0.0);
  CHECK(crossing.b[0] == 0.0);
  CHECK(std::fabs(crossing.b[1] - s) <= 1e188);
  CHECK(std::fabs(crossing.b[2]) <= 1e188);

  const TriTriResult separated =
      tri_tri_intersect(Point3d{0, 0, 0}, Point3d{s, 0, 0}, Point3d{0, s, 0}, Point3d{0, 0, s},
                        Point3d{s, 0, s}, Point3d{0, s, s}, 1e-12);
  CHECK(separated.kind == TriTriKind::none);
}

TEST_CASE("ssi: small coordinates still classify when the tolerance is scaled") {
  const double s = 1e-100;
  const auto intersect = [s](double tolerance) {
    return tri_tri_intersect(Point3d{-s, 0, 0}, Point3d{s, 0, 0}, Point3d{0, s, 0},
                             Point3d{0, -s, -s}, Point3d{0, s, -s}, Point3d{0, s, s}, tolerance);
  };
  const TriTriResult crossing = intersect(1e-112);
  REQUIRE(crossing.kind == TriTriKind::segment);
  CHECK(std::fabs(crossing.a[1]) <= 1e-112);
  CHECK(std::fabs(crossing.a[2]) <= 1e-112);
  CHECK(std::fabs(crossing.b[1] - s) <= 1e-112);
  CHECK(std::fabs(crossing.b[2]) <= 1e-112);

  // At an absolute tolerance far above the geometry the triangle is
  // degenerate, so the result is conservatively touching.
  CHECK(intersect(1e-12).kind == TriTriKind::touching);
}
