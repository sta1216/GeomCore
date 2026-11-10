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

TEST_CASE("ssi: invalid tolerance is touching") {
  const TriTriResult nan_tol =
      tri_tri_intersect(Point3d{-1, 0, 0}, Point3d{1, 0, 0}, Point3d{0, 1, 0}, Point3d{0, -1, -1},
                        Point3d{0, 1, -1}, Point3d{0, 1, 1}, std::nan(""));
  CHECK(nan_tol.kind == TriTriKind::touching);
  const TriTriResult negative_tol =
      tri_tri_intersect(Point3d{-1, 0, 0}, Point3d{1, 0, 0}, Point3d{0, 1, 0}, Point3d{0, -1, -1},
                        Point3d{0, 1, -1}, Point3d{0, 1, 1}, -1e-12);
  CHECK(negative_tol.kind == TriTriKind::touching);
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
