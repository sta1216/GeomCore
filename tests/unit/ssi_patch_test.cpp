#include <doctest/doctest.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <geom/brep/detail/ssi_patch.hpp>
#include <geom/primitives/box.hpp>
#include <geom/primitives/point.hpp>
#include <geom/primitives/vec.hpp>
#include <geom/surfaces/analytic.hpp>
#include <geom/surfaces/nurbs_surface.hpp>
#include <geom/surfaces/plane.hpp>
#include <limits>
#include <set>
#include <utility>
#include <vector>

using namespace geom::v1;

namespace {

double UvSignedArea(const std::array<double, 2>& a, const std::array<double, 2>& b,
                    const std::array<double, 2>& c) {
  return 0.5 * ((b[0] - a[0]) * (c[1] - a[1]) - (c[0] - a[0]) * (b[1] - a[1]));
}

bool FindContainingTriangle(const SsiPatchMesh& mesh, double u, double v, std::size_t* found) {
  for (std::size_t t = 0; t < mesh.triangles.size(); ++t) {
    const std::array<std::uint32_t, 3>& indices = mesh.triangles[t];
    const std::array<double, 2>& a = mesh.uv[indices[0]];
    const std::array<double, 2>& b = mesh.uv[indices[1]];
    const std::array<double, 2>& c = mesh.uv[indices[2]];
    const double area = (b[0] - a[0]) * (c[1] - a[1]) - (c[0] - a[0]) * (b[1] - a[1]);
    if (!(std::fabs(area) > 0.0)) continue;
    const double w0 = ((b[0] - u) * (c[1] - v) - (c[0] - u) * (b[1] - v)) / area;
    const double w1 = ((c[0] - u) * (a[1] - v) - (a[0] - u) * (c[1] - v)) / area;
    const double w2 = 1.0 - w0 - w1;
    const double slack = -1e-12;
    if (w0 >= slack && w1 >= slack && w2 >= slack) {
      *found = t;
      return true;
    }
  }
  return false;
}

double PointTriangleDistance(const Point3d& p, const Point3d& a, const Point3d& b,
                             const Point3d& c) {
  const Vec3d ab = b - a;
  const Vec3d ac = c - a;
  const Vec3d ap = p - a;
  const double d1 = dot(ab, ap);
  const double d2 = dot(ac, ap);
  if (d1 <= 0.0 && d2 <= 0.0) return norm2(ap);

  const Vec3d bp = p - b;
  const double d3 = dot(ab, bp);
  const double d4 = dot(ac, bp);
  if (d3 >= 0.0 && d4 <= d3) return norm2(bp);

  const double vc = d1 * d4 - d3 * d2;
  if (vc <= 0.0 && d1 >= 0.0 && d3 <= 0.0) {
    const double t = d1 / (d1 - d3);
    return norm2(p - (a + ab * t));
  }

  const Vec3d cp = p - c;
  const double d5 = dot(ab, cp);
  const double d6 = dot(ac, cp);
  if (d6 >= 0.0 && d5 <= d6) return norm2(cp);

  const double vb = d5 * d2 - d1 * d6;
  if (vb <= 0.0 && d2 >= 0.0 && d6 <= 0.0) {
    const double t = d2 / (d2 - d6);
    return norm2(p - (a + ac * t));
  }

  const double va = d3 * d6 - d5 * d4;
  if (va <= 0.0 && (d4 - d3) >= 0.0 && (d5 - d6) >= 0.0) {
    const double t = (d4 - d3) / ((d4 - d3) + (d5 - d6));
    return norm2(p - (b + (c - b) * t));
  }

  const double denominator = 1.0 / (va + vb + vc);
  const double v = vb * denominator;
  const double w = vc * denominator;
  return norm2(p - (a + ab * v + ac * w));
}

// NURBS surface whose z control values are the outer product of
// {0, 1/6, -1/6, 0}, giving the bicubic saddle z = P(u) * P(v).
NurbsSurface MakeSaddle() {
  const std::array<double, 4> coefficients{0.0, 1.0 / 6.0, -1.0 / 6.0, 0.0};
  std::vector<Point3d> controls;
  std::vector<double> weights;
  for (int iu = 0; iu < 4; ++iu) {
    for (int iv = 0; iv < 4; ++iv) {
      controls.push_back(Point3d{
          static_cast<double>(iu) / 3.0, static_cast<double>(iv) / 3.0,
          coefficients[static_cast<std::size_t>(iu)] * coefficients[static_cast<std::size_t>(iv)]});
      weights.push_back(1.0);
    }
  }
  return NurbsSurface(3, 3, {0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0},
                      {0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0}, controls, weights);
}

Box3d RawGridBounds(const SsiPatch& patch) {
  Box3d bounds = Box3d::empty();
  for (int iu = 0; iu <= 16; ++iu) {
    for (int iv = 0; iv <= 16; ++iv) {
      const double u = patch.u[0] + (patch.u[1] - patch.u[0]) * (static_cast<double>(iu) / 16.0);
      const double v = patch.v[0] + (patch.v[1] - patch.v[0]) * (static_cast<double>(iv) / 16.0);
      const Result<SurfacePoint> sample = surface_eval(patch.surface, u, v);
      REQUIRE(sample.has_value());
      bounds.expand(sample.value().p);
    }
  }
  return bounds;
}

void CheckBoundsContain(const SsiPatch& patch) {
  const Box3d bounds = patch_bounds(patch, 1e-6);
  REQUIRE_FALSE(bounds.is_empty());
  const Box3d raw = RawGridBounds(patch);
  const double eps = 1e-12 * norm2(bounds.hi - bounds.lo);

  for (int axis = 0; axis < 3; ++axis) {
    CHECK(bounds.lo[axis] < raw.lo[axis]);
    CHECK(bounds.hi[axis] > raw.hi[axis]);
  }

  const int samples = 45;
  for (int iu = 0; iu < samples; ++iu) {
    for (int iv = 0; iv < samples; ++iv) {
      const double u = patch.u[0] + (patch.u[1] - patch.u[0]) * ((static_cast<double>(iu) + 0.5) /
                                                                 static_cast<double>(samples));
      const double v = patch.v[0] + (patch.v[1] - patch.v[0]) * ((static_cast<double>(iv) + 0.5) /
                                                                 static_cast<double>(samples));
      const Result<SurfacePoint> sample = surface_eval(patch.surface, u, v);
      REQUIRE(sample.has_value());
      for (int axis = 0; axis < 3; ++axis) {
        CHECK(sample.value().p[axis] >= bounds.lo[axis] - eps);
        CHECK(sample.value().p[axis] <= bounds.hi[axis] + eps);
      }
    }
  }
}

}  // namespace

TEST_CASE("ssi: a planar patch tessellates into four triangles") {
  const Plane3 plane{{0.0, 0.0, 1.0}, 0.0};
  const Result<SsiPatch> patch = make_ssi_patch(SurfaceGeometry{plane}, {-1.0, 1.0}, {-1.0, 1.0});
  REQUIRE(patch.has_value());
  const Result<SsiPatchMesh> mesh = tessellate_ssi_patch(patch.value(), 1e-9, 12);
  REQUIRE(mesh.has_value());
  CHECK(mesh.value().triangles.size() == 4);
  CHECK(mesh.value().points.size() == 5);
  for (const Point3d& p : mesh.value().points) CHECK(std::fabs(p[2]) < 1e-15);

  double uv_area = 0.0;
  for (const std::array<std::uint32_t, 3>& triangle : mesh.value().triangles) {
    const double area = UvSignedArea(mesh.value().uv[triangle[0]], mesh.value().uv[triangle[1]],
                                     mesh.value().uv[triangle[2]]);
    CHECK(area > 0.0);
    uv_area += area;
  }
  CHECK(uv_area == doctest::Approx(4.0).epsilon(1e-12));
  CHECK(mesh.value().max_chord_error == 0.0);
  CHECK(mesh.value().tolerance_met);
}

TEST_CASE("ssi: a sphere patch obeys the chord tolerance") {
  const Sphere3 sphere{Point3d{0, 0, 0}, Vec3d{0, 0, 1}, Vec3d{1, 0, 0}, 1.0};
  const Result<SsiPatch> patch = make_ssi_patch(SurfaceGeometry{sphere}, {0.0, 1.0}, {-0.5, 0.5});
  REQUIRE(patch.has_value());
  const Result<SsiPatchMesh> mesh = tessellate_ssi_patch(patch.value(), 1e-3, 12);
  REQUIRE(mesh.has_value());
  CHECK(mesh.value().triangles.size() > 8);
  for (const Point3d& p : mesh.value().points) {
    CHECK(std::fabs(std::sqrt(p[0] * p[0] + p[1] * p[1] + p[2] * p[2]) - 1.0) < 1e-12);
  }
  for (const auto& tri : mesh.value().triangles) {
    const Point3d& a = mesh.value().points[tri[0]];
    const Point3d& b = mesh.value().points[tri[1]];
    const Point3d& c = mesh.value().points[tri[2]];
    const Vec3d normal = cross(b - a, c - a);
    CHECK(norm2(normal) > 0.0);
  }

  CHECK(mesh.value().tolerance_met);
  CHECK(mesh.value().max_chord_error <= 1e-3);

  double measured = 0.0;
  const int grid = 8;
  for (int iu = 0; iu < grid; ++iu) {
    for (int iv = 0; iv < grid; ++iv) {
      const double u = (static_cast<double>(iu) + 0.5) / static_cast<double>(grid);
      const double v = -0.5 + (static_cast<double>(iv) + 0.5) / static_cast<double>(grid);
      std::size_t triangle = 0;
      REQUIRE(FindContainingTriangle(mesh.value(), u, v, &triangle));
      const std::array<std::uint32_t, 3>& indices = mesh.value().triangles[triangle];
      const Result<SurfacePoint> sample = surface_eval(sphere, u, v);
      REQUIRE(sample.has_value());
      measured = std::max(
          measured,
          PointTriangleDistance(sample.value().p, mesh.value().points[indices[0]],
                                mesh.value().points[indices[1]], mesh.value().points[indices[2]]));
    }
  }
  CHECK(measured <= 1e-3);
}

TEST_CASE("ssi: patch bounds contain every tessellated point") {
  const Torus3 torus{Point3d{0, 0, 0}, Vec3d{0, 0, 1}, Vec3d{1, 0, 0}, 3.0, 1.0};
  const Result<SsiPatch> patch = make_ssi_patch(SurfaceGeometry{torus}, {0.0, 2.0}, {0.0, 2.0});
  REQUIRE(patch.has_value());
  const Result<SsiPatchMesh> mesh = tessellate_ssi_patch(patch.value(), 1e-2, 12);
  REQUIRE(mesh.has_value());
  const Box3d bounds = patch_bounds(patch.value(), 1e-2);
  CHECK_FALSE(bounds.is_empty());
  for (const Point3d& p : mesh.value().points) {
    for (int i = 0; i < 3; ++i) {
      CHECK(p[i] >= bounds.lo[i] - 1e-12);
      CHECK(p[i] <= bounds.hi[i] + 1e-12);
    }
  }
}

TEST_CASE("ssi: tessellation rejects invalid tolerances") {
  const Plane3 plane{};
  const Result<SsiPatch> patch = make_ssi_patch(SurfaceGeometry{plane}, {-1.0, 1.0}, {-1.0, 1.0});
  REQUIRE(patch.has_value());
  CHECK_FALSE(tessellate_ssi_patch(patch.value(), 0.0, 12).has_value());
  CHECK_FALSE(tessellate_ssi_patch(patch.value(), -1.0, 12).has_value());
  CHECK_FALSE(tessellate_ssi_patch(patch.value(), 1e-9, 0).has_value());
}

TEST_CASE("ssi: tessellation rejects non-finite tolerance, bad depth and bad domains") {
  const Plane3 plane{};
  const Result<SsiPatch> patch = make_ssi_patch(SurfaceGeometry{plane}, {-1.0, 1.0}, {-1.0, 1.0});
  REQUIRE(patch.has_value());

  const Result<SsiPatchMesh> infinite =
      tessellate_ssi_patch(patch.value(), std::numeric_limits<double>::infinity(), 12);
  REQUIRE_FALSE(infinite.has_value());
  CHECK(infinite.error().code == Errc::tolerance_violation);

  const Result<SsiPatchMesh> negative_depth = tessellate_ssi_patch(patch.value(), 1e-9, -1);
  REQUIRE_FALSE(negative_depth.has_value());
  CHECK(negative_depth.error().code == Errc::index_out_of_range);

  SsiPatch invalid_surface;
  invalid_surface.surface =
      SurfaceGeometry{Cylinder3{Point3d{0, 0, 0}, Vec3d{0, 0, 1}, Vec3d{1, 0, 0}, 0.0}};
  invalid_surface.u = {0.0, 1.0};
  invalid_surface.v = {0.0, 1.0};
  const Result<SsiPatchMesh> bad_surface = tessellate_ssi_patch(invalid_surface, 1e-9, 12);
  REQUIRE_FALSE(bad_surface.has_value());
  CHECK(bad_surface.error().code == Errc::degenerate);

  SsiPatch inverted = patch.value();
  inverted.u = {1.0, -1.0};
  const Result<SsiPatchMesh> bad_domain = tessellate_ssi_patch(inverted, 1e-9, 12);
  REQUIRE_FALSE(bad_domain.has_value());
  CHECK(bad_domain.error().code == Errc::degenerate);
}

TEST_CASE("ssi: tessellation is deterministic") {
  const Sphere3 sphere{Point3d{0, 0, 0}, Vec3d{0, 0, 1}, Vec3d{1, 0, 0}, 1.0};
  const Result<SsiPatch> patch = make_ssi_patch(SurfaceGeometry{sphere}, {0.0, 1.5}, {-0.4, 0.7});
  REQUIRE(patch.has_value());
  const Result<SsiPatchMesh> first = tessellate_ssi_patch(patch.value(), 1e-2, 12);
  const Result<SsiPatchMesh> second = tessellate_ssi_patch(patch.value(), 1e-2, 12);
  REQUIRE(first.has_value());
  REQUIRE(second.has_value());
  CHECK(first.value().max_chord_error == second.value().max_chord_error);
  CHECK(first.value().tolerance_met == second.value().tolerance_met);
  REQUIRE(first.value().points.size() == second.value().points.size());
  REQUIRE(first.value().normals.size() == second.value().normals.size());
  REQUIRE(first.value().uv.size() == second.value().uv.size());
  REQUIRE(first.value().triangles.size() == second.value().triangles.size());
  for (std::size_t i = 0; i < first.value().points.size(); ++i) {
    for (int axis = 0; axis < 3; ++axis) {
      CHECK(first.value().points[i][axis] == second.value().points[i][axis]);
      CHECK(first.value().normals[i][axis] == second.value().normals[i][axis]);
    }
    CHECK(first.value().uv[i][0] == second.value().uv[i][0]);
    CHECK(first.value().uv[i][1] == second.value().uv[i][1]);
  }
  for (std::size_t i = 0; i < first.value().triangles.size(); ++i) {
    for (int k = 0; k < 3; ++k) {
      CHECK(first.value().triangles[i][k] == second.value().triangles[i][k]);
    }
  }
}

TEST_CASE("ssi: patch normals are unit and follow the sphere") {
  const Sphere3 sphere{Point3d{0, 0, 0}, Vec3d{0, 0, 1}, Vec3d{1, 0, 0}, 1.0};
  const Result<SsiPatch> patch = make_ssi_patch(SurfaceGeometry{sphere}, {0.0, 1.0}, {-0.5, 0.5});
  REQUIRE(patch.has_value());
  const Result<SsiPatchMesh> mesh = tessellate_ssi_patch(patch.value(), 1e-2, 12);
  REQUIRE(mesh.has_value());
  REQUIRE(mesh.value().normals.size() == mesh.value().points.size());
  REQUIRE(mesh.value().uv.size() == mesh.value().points.size());
  for (std::size_t i = 0; i < mesh.value().points.size(); ++i) {
    const Vec3d& normal = mesh.value().normals[i];
    CHECK(std::fabs(std::hypot(normal[0], normal[1], normal[2]) - 1.0) < 1e-12);
    for (int axis = 0; axis < 3; ++axis) {
      CHECK(normal[axis] == doctest::Approx(mesh.value().points[i][axis]).epsilon(1e-12));
    }
  }
}

TEST_CASE("ssi: a NURBS patch tessellates from its derivatives") {
  const NurbsSurface surface(
      1, 1, {0.0, 0.0, 1.0, 1.0}, {0.0, 0.0, 1.0, 1.0},
      {Point3d{0, 0, 0}, Point3d{0, 1, 0}, Point3d{1, 0, 0}, Point3d{1, 1, 0}},
      {1.0, 1.0, 1.0, 1.0});
  REQUIRE(surface.is_valid());
  const Result<SsiPatch> patch = make_ssi_patch(SurfaceGeometry{surface}, {0.0, 1.0}, {0.0, 1.0});
  REQUIRE(patch.has_value());
  const Result<SsiPatchMesh> mesh = tessellate_ssi_patch(patch.value(), 1e-9, 12);
  REQUIRE(mesh.has_value());
  CHECK(mesh.value().triangles.size() == 4);
  CHECK(mesh.value().points.size() == 5);
  for (const Point3d& p : mesh.value().points) CHECK(std::fabs(p[2]) < 1e-15);
  for (const Vec3d& normal : mesh.value().normals) {
    CHECK(std::fabs(std::fabs(normal[2]) - 1.0) < 1e-12);
  }
}

TEST_CASE("ssi: a saddle above tolerance at the depth cap is reported") {
  const NurbsSurface saddle = MakeSaddle();
  REQUIRE(saddle.is_valid());
  const Result<SsiPatch> patch = make_ssi_patch(SurfaceGeometry{saddle}, {0.0, 1.0}, {0.0, 1.0});
  REQUIRE(patch.has_value());
  const Result<SsiPatchMesh> mesh = tessellate_ssi_patch(patch.value(), 1e-9, 4);
  REQUIRE(mesh.has_value());
  CHECK_FALSE(mesh.value().tolerance_met);
  CHECK(mesh.value().max_chord_error > 1e-4);
  CHECK(mesh.value().triangles.size() > 4);
}

TEST_CASE("ssi: patch_bounds contains analytic and saddle patches") {
  const Plane3 plane{};
  const Sphere3 sphere{Point3d{0, 0, 0}, Vec3d{0, 0, 1}, Vec3d{1, 0, 0}, 1.0};
  const Cylinder3 cylinder{Point3d{0, 0, 0}, Vec3d{0, 0, 1}, Vec3d{1, 0, 0}, 1.0};
  const Torus3 torus{Point3d{0, 0, 0}, Vec3d{0, 0, 1}, Vec3d{1, 0, 0}, 3.0, 1.0};
  const NurbsSurface saddle = MakeSaddle();

  const Result<SsiPatch> plane_patch =
      make_ssi_patch(SurfaceGeometry{plane}, {-1.0, 1.0}, {-1.0, 1.0});
  const Result<SsiPatch> sphere_patch =
      make_ssi_patch(SurfaceGeometry{sphere}, {0.0, 1.0}, {-0.7, 0.7});
  const Result<SsiPatch> cylinder_patch =
      make_ssi_patch(SurfaceGeometry{cylinder}, {0.0, 2.0}, {-1.0, 1.0});
  const Result<SsiPatch> torus_patch =
      make_ssi_patch(SurfaceGeometry{torus}, {0.0, 2.0}, {0.0, 2.0});
  const Result<SsiPatch> saddle_patch =
      make_ssi_patch(SurfaceGeometry{saddle}, {0.0, 1.0}, {0.0, 1.0});
  REQUIRE(plane_patch.has_value());
  REQUIRE(sphere_patch.has_value());
  REQUIRE(cylinder_patch.has_value());
  REQUIRE(torus_patch.has_value());
  REQUIRE(saddle_patch.has_value());

  CheckBoundsContain(plane_patch.value());
  CheckBoundsContain(sphere_patch.value());
  CheckBoundsContain(cylinder_patch.value());
  CheckBoundsContain(torus_patch.value());
  CheckBoundsContain(saddle_patch.value());
}

TEST_CASE("ssi: a planar mesh deduplicates vertices") {
  const Plane3 plane{};
  const Result<SsiPatch> patch = make_ssi_patch(SurfaceGeometry{plane}, {-1.0, 1.0}, {-1.0, 1.0});
  REQUIRE(patch.has_value());
  const Result<SsiPatchMesh> mesh = tessellate_ssi_patch(patch.value(), 1e-9, 12);
  REQUIRE(mesh.has_value());
  CHECK(mesh.value().points.size() == 5);
  std::set<std::uint32_t> used;
  for (const std::array<std::uint32_t, 3>& triangle : mesh.value().triangles) {
    used.insert(triangle.begin(), triangle.end());
  }
  CHECK(used.size() == 5);

  const Sphere3 sphere{Point3d{0, 0, 0}, Vec3d{0, 0, 1}, Vec3d{1, 0, 0}, 1.0};
  const Result<SsiPatch> sphere_patch =
      make_ssi_patch(SurfaceGeometry{sphere}, {0.0, 1.0}, {-0.5, 0.5});
  REQUIRE(sphere_patch.has_value());
  const Result<SsiPatchMesh> sphere_mesh = tessellate_ssi_patch(sphere_patch.value(), 1e-2, 12);
  REQUIRE(sphere_mesh.has_value());
  std::set<std::pair<double, double>> unique_uv;
  for (const std::array<double, 2>& uv : sphere_mesh.value().uv) {
    unique_uv.insert({uv[0], uv[1]});
  }
  CHECK(unique_uv.size() == sphere_mesh.value().uv.size());
}

TEST_CASE("ssi: a patch reaching a sphere pole stays finite") {
  const Sphere3 sphere{Point3d{0, 0, 0}, Vec3d{0, 0, 1}, Vec3d{1, 0, 0}, 1.0};
  const double half_pi = std::acos(0.0);
  const Result<SsiPatch> patch =
      make_ssi_patch(SurfaceGeometry{sphere}, {0.0, 2.0}, {1.2, half_pi});
  REQUIRE(patch.has_value());
  const Result<SsiPatchMesh> mesh = tessellate_ssi_patch(patch.value(), 1e-2, 12);
  REQUIRE(mesh.has_value());
  for (const Point3d& p : mesh.value().points) {
    CHECK(std::isfinite(p[0]));
    CHECK(std::isfinite(p[1]));
    CHECK(std::isfinite(p[2]));
  }
  for (const Vec3d& normal : mesh.value().normals) {
    CHECK(std::fabs(std::hypot(normal[0], normal[1], normal[2]) - 1.0) < 1e-12);
  }
}

TEST_CASE("ssi: patch_bounds is empty for invalid inputs") {
  const Plane3 plane{};
  const Result<SsiPatch> patch = make_ssi_patch(SurfaceGeometry{plane}, {-1.0, 1.0}, {-1.0, 1.0});
  REQUIRE(patch.has_value());
  CHECK(patch_bounds(patch.value(), 0.0).is_empty());
  CHECK(patch_bounds(patch.value(), -1.0).is_empty());
  SsiPatch inverted = patch.value();
  inverted.u = {1.0, -1.0};
  CHECK(patch_bounds(inverted, 1e-3).is_empty());
}
