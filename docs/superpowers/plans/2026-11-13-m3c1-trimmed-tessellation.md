# M3c-1 Trimmed Curved Tessellation and Mass Properties Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use subagent-driven-development (recommended) or executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Tessellate trimmed analytic B-rep faces (holes, partial rims, seam crossings) into a watertight mesh within `deflection`, and compute volume, centroid and inertia for solids that contain those faces.

**Architecture:** A trimmed face contributes a set of UV rings (outer CCW, holes CW) obtained by projecting its coedge edge samples, so every boundary vertex is the same pool index the neighbouring face uses. The rings are triangulated by the existing exact `triangulate_polygon2`, lifted with `face_uv_eval`, and refined uniformly per face until the chord height meets `deflection`. `mass_properties` reuses the same ring/triangulation machinery with arc-aware boundary resampling and a seven-point Dunavant rule that accumulates the first and second moments. Canonical untrimmed primitives keep their existing grid builders, so the current weld and pole guarantees are untouched.

**Tech Stack:** C++20, CMake 4.4.4 + VS2022 (MSVC 19.44), doctest 2.4.11, `geom::v1`, exact predicates from `geom/cgeom`, `/W4 /WX /utf-8`.

**Spec:** `docs/architecture/geometry-library-architecture.md` §8.10 (tessellation contract), §8.11 (mass properties), §16.3; `docs/quality/testing-strategy.md`.

## Scope Split

M3c is too large for one plan. This plan is **M3c-1**:

- trimmed tessellation for the surfaces M3b-2 already produces (cylinder bands and caps, plus any analytic fragment);
- trimmed mass properties for those solids;
- QA-13 and documentation.

**M3c-2** (separate plan) covers sphere/torus/cone operands in `boolean_solid_curved` and their pole/apex tessellation. NURBS faces stay `unsupported` in both plans: no pcurve projection exists for a NURBS surface, so a NURBS trim contour cannot be formed yet.

## Global Constraints

- C++20, `#pragma once`, everything in `namespace geom { inline namespace v1 { ... } }`.
- No global epsilon and no tolerance widening to make a test pass; every rejection is a typed `Error{Errc, message, entity, context}`.
- Topology decisions use the exact predicates (`sign_orient2d`, `ring_orientation2`); floating point is only for sampling and diagnostics.
- Every new test file is added to `GEOM_TEST_SOURCES` in `tests/CMakeLists.txt`; Debug and Release must both pass, with zero warnings under `/W4 /WX /utf-8`.
- `tessellate()` keeps its contract: shared edge samples are bit-identical between coedges, triangles carry `triangle_face`, orientation follows `Face::reversed`.
- `mass_properties()` keeps rejecting open shells, non-positive accumulated volume and non-finite integrals; it must never silently integrate a wrong region.
- Each task is TDD: failing test first (`RED` evidence recorded in the report), minimal implementation, Green in both configurations, one commit.

## Review Focus

1. **Hole touching the outer ring** at a single point or along an edge: triangulation must not emit zero-area or overlapping triangles; if the domain is not a valid polygon, return `Errc::degenerate` with the face id instead of a partial mesh.
2. **Seam-crossing rim**: a loop whose u coordinate leaves `[0, 2π)` must stay one connected ring; no duplicated shadow ring and no T-junction at the seam.
3. **Trimmed face sharing its boundary with a planar cut face**: boundary indices must be shared exactly; a mesh that is watertight only after a tolerance weld is a failure.
4. **Deflection far below the feature size**: boundary resampling must cap the samples per ring (1,000,000) and return `Errc::numeric_failure` before allocating, never attempt a billion-point ring.
5. **Mass properties of a solid with a trimmed curved face**: the trimmed region must not be double counted with its planar neighbours; volume must agree with `detail::curved_solid_volume` and with the analytic value.

---

### Task 1: Trimmed UV domain extraction and triangulation

**Files:**
- Create: `include/geom/brep/trimmed_domain.hpp`
- Create: `tests/unit/trimmed_domain_test.cpp`
- Modify: `tests/CMakeLists.txt` (append `unit/trimmed_domain_test.cpp`)

**Interfaces:**
- Consumes: `face_loops_uv(const Brep&, FaceId, double)` → `Result<std::vector<FaceUvLoop>>`; `FaceUvLoop{points, outer, arcs}`; `FaceUvArcRef{CircleCurve3 circle; double t0; double t1;}`; `triangulate_polygon2(const Polygon2d&)` → `Triangulation2{triangles, ok}`; `ring_orientation2`; `signed_area2`.
- Produces:

```cpp
struct TrimmedBoundary {
  std::vector<Point2d> points;   // implicitly closed, no repeated last point
  std::vector<std::optional<FaceUvArcRef>> arcs;  // arcs[i]: points[i] -> points[i+1]
  bool outer = false;
};

struct TrimmedUvTriangulation {
  std::vector<Point2d> points;                      // concatenated ring vertices
  std::vector<std::array<std::uint32_t, 3>> triangles;
};

Result<std::vector<TrimmedBoundary>> trimmed_face_boundaries(const Brep& brep, FaceId face,
                                                             const FaceUvProjection& projection,
                                                             double deflection);

Result<TrimmedUvTriangulation> triangulate_trimmed_domain(
    std::span<const TrimmedBoundary> boundaries);
```

- [x] **Step 1: Write the failing tests**

`tests/unit/trimmed_domain_test.cpp`:

```cpp
#include <doctest/doctest.h>

#include <geom/brep/boolean_solid_curved.hpp>
#include <geom/brep/face_uv.hpp>
#include <geom/brep/primitives.hpp>
#include <geom/brep/trimmed_domain.hpp>
#include <geom/core/result.hpp>
#include <geom/cgeom/polygon2.hpp>

#include <cmath>
#include <numbers>
#include <vector>

namespace {

constexpr double kPi = std::numbers::pi_v<double>;

geom::Brep cut_cylinder() {
  const geom::Brep cylinder =
      geom::make_cylinder(geom::Point3d{0.0, 0.0, 0.0}, geom::Vec3d{0.0, 0.0, 1.0}, 1.0, 2.0,
                          1e-9)
          .value();
  const geom::Brep box =
      geom::make_box(geom::Point3d{0.5, -2.0, 0.5}, geom::Point3d{2.0, 2.0, 1.5}, 1e-9).value();
  const geom::Result<geom::CurvedBooleanResult> result = geom::boolean_solid_curved(
      cylinder, cylinder.solid_ids().front(), box, box.solid_ids().front(),
      geom::SolidBooleanOp::difference, geom::CurvedBooleanOptions{1e-9});
  REQUIRE(result.has_value());
  return result.value().brep;
}

std::vector<geom::FaceId> cylinder_faces(const geom::Brep& brep) {
  std::vector<geom::FaceId> faces;
  for (const geom::FaceId id : brep.face_ids()) {
    const geom::Face& face = brep.face(id);
    if (face.surface < brep.geometry().surfaces.size() &&
        std::holds_alternative<geom::Cylinder3>(brep.geometry().surfaces[face.surface])) {
      faces.push_back(id);
    }
  }
  return faces;
}

double polygon_area(const std::vector<geom::TrimmedBoundary>& boundaries) {
  double total = 0.0;
  for (const geom::TrimmedBoundary& ring : boundaries) {
    total += geom::signed_area2(ring.points);
  }
  return total;
}

}  // namespace

TEST_CASE("trimmed domain: a side window leaves one hole across the seam") {
  const geom::Brep brep = cut_cylinder();
  const std::vector<geom::FaceId> faces = cylinder_faces(brep);
  REQUIRE(!faces.empty());

  std::size_t holes = 0;
  double lateral_area = 0.0;
  for (const geom::FaceId id : faces) {
    const geom::Result<geom::FaceUvProjection> projection =
        geom::face_uv_projection(brep, brep.face(id), 1e-9);
    REQUIRE(projection.has_value());
    const geom::Result<std::vector<geom::TrimmedBoundary>> boundaries =
        geom::trimmed_face_boundaries(brep, id, projection.value(), 1e-2);
    REQUIRE(boundaries.has_value());
    REQUIRE(!boundaries.value().empty());
    std::size_t outer = 0;
    for (const geom::TrimmedBoundary& ring : boundaries.value()) {
      CHECK(ring.points.size() >= 3);
      CHECK(ring.arcs.size() == ring.points.size());
      CHECK(geom::signed_area2(ring.points) != 0.0);
      for (std::size_t i = 0; i < ring.points.size(); ++i) {
        const geom::Point2d& a = ring.points[i];
        const geom::Point2d& b = ring.points[(i + 1) % ring.points.size()];
        CHECK(std::abs(b[0] - a[0]) < kPi);  // unwrapped: no jump at the seam
      }
      if (ring.outer) {
        ++outer;
        CHECK(geom::signed_area2(ring.points) > 0.0);
      } else {
        ++holes;
        CHECK(geom::signed_area2(ring.points) < 0.0);
      }
    }
    CHECK(outer == 1);
    lateral_area += polygon_area(boundaries.value());
  }
  CHECK(holes == 1);
  // Lateral band 4*pi minus the window: x = 0.5 cuts the unit circle at
  // theta = +-pi/3, so the window spans 2*pi/3 radians and height 1.
  CHECK(lateral_area ==
        doctest::Approx(2.0 * kPi * 1.0 * 2.0 - 2.0 * kPi / 3.0).epsilon(2e-3));
}

TEST_CASE("trimmed domain: deflection adds boundary samples monotonically") {
  const geom::Brep brep = cut_cylinder();
  const geom::FaceId id = cylinder_faces(brep).front();
  const geom::FaceUvProjection projection =
      geom::face_uv_projection(brep, brep.face(id), 1e-9).value();

  const auto count = [&](double deflection) {
    const auto boundaries = geom::trimmed_face_boundaries(brep, id, projection, deflection);
    REQUIRE(boundaries.has_value());
    std::size_t total = 0;
    for (const geom::TrimmedBoundary& ring : boundaries.value()) total += ring.points.size();
    return total;
  };
  CHECK(count(1e-3) > count(1e-2));
  CHECK(count(1e-4) >= count(1e-3));
}

TEST_CASE("trimmed domain: triangulation covers the ring area exactly") {
  const geom::Brep brep = cut_cylinder();
  const geom::FaceId id = cylinder_faces(brep).front();
  const geom::FaceUvProjection projection =
      geom::face_uv_projection(brep, brep.face(id), 1e-9).value();
  const auto boundaries = geom::trimmed_face_boundaries(brep, id, projection, 5e-3);
  REQUIRE(boundaries.has_value());

  const auto triangulation = geom::triangulate_trimmed_domain(boundaries.value());
  REQUIRE(triangulation.has_value());
  CHECK(triangulation.value().points.size() ==
        [&] {
          std::size_t total = 0;
          for (const auto& ring : boundaries.value()) total += ring.points.size();
          return total;
        }());
  REQUIRE(!triangulation.value().triangles.empty());

  double area = 0.0;
  for (const std::array<std::uint32_t, 3>& triangle : triangulation.value().triangles) {
    const geom::Point2d& a = triangulation.value().points[triangle[0]];
    const geom::Point2d& b = triangulation.value().points[triangle[1]];
    const geom::Point2d& c = triangulation.value().points[triangle[2]];
    const double twice = (b[0] - a[0]) * (c[1] - a[1]) - (b[1] - a[1]) * (c[0] - a[0]);
    CHECK(twice > 0.0);
    area += 0.5 * twice;
  }
  CHECK(area == doctest::Approx(polygon_area(boundaries.value())).epsilon(1e-12));
}

TEST_CASE("trimmed domain: triangulation subtracts a hole ring") {
  geom::TrimmedBoundary outer;
  outer.outer = true;
  outer.points = {{0.0, 0.0}, {4.0, 0.0}, {4.0, 4.0}, {0.0, 4.0}};  // CCW
  outer.arcs.resize(outer.points.size());

  geom::TrimmedBoundary hole;
  hole.outer = false;
  hole.points = {{1.0, 1.0}, {1.0, 3.0}, {3.0, 3.0}, {3.0, 1.0}};  // CW
  hole.arcs.resize(hole.points.size());

  const std::vector<geom::TrimmedBoundary> rings{outer, hole};
  const auto triangulation = geom::triangulate_trimmed_domain(rings);
  REQUIRE(triangulation.has_value());
  CHECK(triangulation.value().points.size() == 8);
  double area = 0.0;
  for (const std::array<std::uint32_t, 3>& triangle : triangulation.value().triangles) {
    const geom::Point2d& a = triangulation.value().points[triangle[0]];
    const geom::Point2d& b = triangulation.value().points[triangle[1]];
    const geom::Point2d& c = triangulation.value().points[triangle[2]];
    area += 0.5 * ((b[0] - a[0]) * (c[1] - a[1]) - (b[1] - a[1]) * (c[0] - a[0]));
  }
  CHECK(area == doctest::Approx(12.0).epsilon(1e-12));  // 4x4 minus the 2x2 hole
}

TEST_CASE("trimmed domain: duplicate or non-finite rings are rejected") {
  const geom::Brep brep = cut_cylinder();
  const geom::FaceId id = cylinder_faces(brep).front();
  const geom::FaceUvProjection projection =
      geom::face_uv_projection(brep, brep.face(id), 1e-9).value();
  auto boundaries = geom::trimmed_face_boundaries(brep, id, projection, 1e-2).value();

  std::vector<geom::TrimmedBoundary> duplicated = boundaries;
  duplicated.push_back(boundaries.front());
  CHECK(!geom::triangulate_trimmed_domain(duplicated).has_value());

  std::vector<geom::TrimmedBoundary> non_finite = boundaries;
  non_finite.front().points.front() = geom::Point2d{std::nan(""), 0.0};
  CHECK(!geom::triangulate_trimmed_domain(non_finite).has_value());

  const auto too_fine = geom::trimmed_face_boundaries(brep, id, projection, 1e-12);
  REQUIRE(!too_fine.has_value());
  CHECK(too_fine.error().code == geom::Errc::numeric_failure);
}
```

- [x] **Step 2: Run the tests to verify they fail**

Run:
```powershell
cmake --build build --config Debug --target geom_tests
.\build\tests\Debug\geom_tests.exe --test-case="trimmed domain:*"
```
Expected: `trimmed_domain.hpp: No such file or directory` (C1083).

- [x] **Step 3: Implement the domain and triangulation**

Implementation notes (signature and tests leave the body open, these points are not):
- Call `face_loops_uv(brep, face, projection.tolerance)`; require one `outer` loop and that every other loop is a hole; otherwise `Errc::invalid_topology`.
- For each segment `i -> i+1`, if `arcs[i]` is engaged, resample it from `circle` with step `2 * acos(max(-1, min(1, 1 - deflection / radius)))`, append `face_uv_project(projection, eval(t), previous_u, true)`; append the segment end if the last sample is not within `1e-12` of it. Straight segments contribute their end point only. Compute the sample count before allocating: an arc that needs more than 1,000,000 segments, or a ring whose total exceeds 1,000,000 points, returns `Errc::numeric_failure` naming the face.
- Drop consecutive duplicates (within `1e-12`) and a closing duplicate; require at least 3 points per ring and a non-zero area; re-unwrap each ring so consecutive `u` differ by less than `pi` (add or subtract `2*pi` as needed, keeping the seam-crossing ring connected); normalize outer to CCW and holes to CW.
- `triangulate_trimmed_domain`: reject fewer than one outer ring, more than one outer ring, non-finite points, or duplicate rings (same first point and vertex count); build `Polygon2d` with `ring_starts` (outer first), translate by the first vertex before calling `triangulate_polygon2` only if needed for the predicate range (keep a mapping back); return `Errc::degenerate` when `!triangulation.ok`. Verify every returned triangle has positive orientation and non-zero area.

- [x] **Step 4: Run the tests to verify they pass**

Run the same command; expected: all four cases pass, plus the whole suite still green in Debug and Release.

- [x] **Step 5: Commit**

```powershell
git add include/geom/brep/trimmed_domain.hpp tests/unit/trimmed_domain_test.cpp tests/CMakeLists.txt
git commit -m "feat(m3c1): extract and triangulate trimmed face domains"
```

---

### Task 2: Trimmed face tessellation inside `tessellate()`

**Files:**
- Modify: `include/geom/brep/tessellate.hpp` (after `tessellate_planar_face`, before `sample_curved_face`)
- Create: `tests/unit/tessellate_trimmed_test.cpp`
- Modify: `tests/CMakeLists.txt`

**Interfaces:**
- Consumes: `TessellateState`, `tessellate_edge_positions(state, brep, EdgeId, deflection)`, `face_uv_projection`, `face_uv_project`, `face_uv_eval`, `trimmed_domain.hpp`, `triangulate_polygon2`.
- Produces (internal, used by Task 3 and by the fallback):

```cpp
namespace detail {
Result<void> tessellate_trimmed_face(TessellateState& state, const Brep& brep, FaceId id,
                                     const FaceUvProjection& projection, double deflection);
}
```

- [x] **Step 1: Write the failing test**

`tests/unit/tessellate_trimmed_test.cpp`:

```cpp
#include <doctest/doctest.h>

#include <geom/brep/boolean_solid_curved.hpp>
#include <geom/brep/detail/curved_volume.hpp>
#include <geom/brep/primitives.hpp>
#include <geom/brep/tessellate.hpp>
#include <geom/brep/validate.hpp>

#include <cmath>
#include <map>
#include <numbers>
#include <utility>

namespace {

geom::Brep cut_cylinder() {
  // Same fixture as the trimmed-domain tests: a side window that crosses the seam.
  const geom::Brep cylinder =
      geom::make_cylinder(geom::Point3d{0.0, 0.0, 0.0}, geom::Vec3d{0.0, 0.0, 1.0}, 1.0, 2.0,
                          1e-9)
          .value();
  const geom::Brep box =
      geom::make_box(geom::Point3d{0.5, -2.0, 0.5}, geom::Point3d{2.0, 2.0, 1.5}, 1e-9).value();
  return geom::boolean_solid_curved(cylinder, cylinder.solid_ids().front(), box,
                                    box.solid_ids().front(), geom::SolidBooleanOp::difference,
                                    geom::CurvedBooleanOptions{1e-9})
      .value()
      .brep;
}

void check_watertight(const geom::TriMesh& mesh) {
  std::map<std::pair<std::uint32_t, std::uint32_t>, int> uses;
  for (const auto& triangle : mesh.triangles) {
    for (int k = 0; k < 3; ++k) {
      const std::uint32_t a = triangle[static_cast<std::size_t>(k)];
      const std::uint32_t b = triangle[static_cast<std::size_t>((k + 1) % 3)];
      uses[{std::min(a, b), std::max(a, b)}] += 1;
    }
  }
  for (const auto& [edge, count] : uses) CHECK(count == 2);
}

double mesh_volume(const geom::TriMesh& mesh) {
  double volume = 0.0;
  for (const auto& triangle : mesh.triangles) {
    const geom::Point3d& p = mesh.positions[triangle[0]];
    const geom::Point3d& q = mesh.positions[triangle[1]];
    const geom::Point3d& r = mesh.positions[triangle[2]];
    const geom::Vec3d a{p[0], p[1], p[2]};
    const geom::Vec3d b{q[0], q[1], q[2]};
    const geom::Vec3d c{r[0], r[1], r[2]};
    volume += geom::dot(a, geom::cross(b, c)) / 6.0;
  }
  return volume;
}

}  // namespace

TEST_CASE("trimmed tessellation: windowed cylinder is watertight and keeps its faces") {
  const geom::Brep brep = cut_cylinder();
  REQUIRE(geom::validate_topology(brep).ok());
  const auto mesh = geom::tessellate(brep, geom::TessellationOptions{1e-3, 0.35});
  REQUIRE(mesh.has_value());
  check_watertight(mesh.value());
  CHECK(mesh.value().triangle_face.size() == mesh.value().triangles.size());
  for (const std::uint32_t face : mesh.value().triangle_face) {
    CHECK(face < brep.face_ids().size());
  }
}

TEST_CASE("trimmed tessellation: every vertex lies on its analytic surface") {
  const geom::Brep brep = cut_cylinder();
  const auto mesh = geom::tessellate(brep, geom::TessellationOptions{2e-3, 0.35});
  REQUIRE(mesh.has_value());
  for (std::size_t i = 0; i < mesh.value().positions.size(); ++i) {
    const geom::Point3d& p = mesh.value().positions[i];
    const double radial = std::hypot(p[0], p[1]);
    const bool on_side = std::abs(radial - 1.0) <= 2e-3;
    const bool on_cap = std::abs(p[2]) <= 2e-3 || std::abs(p[2] - 2.0) <= 2e-3;
    const bool on_window = std::abs(p[0] - 0.5) <= 2e-3 || std::abs(p[0] + 0.5) <= 2e-3 ||
                           std::abs(p[2] - 0.5) <= 2e-3 || std::abs(p[2] - 1.5) <= 2e-3;
    CHECK((on_side || on_cap || on_window));
  }
}

TEST_CASE("trimmed tessellation: canonical primitives are unchanged") {
  const geom::Brep cylinder =
      geom::make_cylinder(geom::Point3d{0.0, 0.0, 0.0}, geom::Vec3d{0.0, 0.0, 1.0}, 1.0, 2.0,
                          1e-9)
          .value();
  const auto mesh = geom::tessellate(cylinder, geom::TessellationOptions{1e-2, 0.35});
  REQUIRE(mesh.has_value());
  check_watertight(mesh.value());
  CHECK(mesh_volume(mesh.value()) == doctest::Approx(2.0 * std::numbers::pi_v<double>).epsilon(1e-3));
}
```

- [x] **Step 2: Run the tests to verify they fail**

Expected: the windowed case fails with the current `unsupported` error (the canonical builder rejects the trimmed band), and the canonical-primitives case still passes.

- [x] **Step 3: Implement the trimmed builder**

Implementation notes:
- Build one ring per face loop by walking coedges in loop order. For each coedge, call `tessellate_edge_positions(state, brep, coedge.edge, deflection)`, reverse the sample order when `coedge.reversed`, and project each position with `face_uv_project(projection, point, previous_u, has_previous_u)`. Append the sample positions (not freshly sampled points) to the ring and remember the pool index, so shared edges reuse the neighbour's indices.
- Classify the loop as outer/hole with `face_loops_uv` (same loop order) or, when that is not available, by nesting depth of the projected rings; require exactly one outer.
- Triangulate with `triangulate_polygon2` on the UV points (translate the first point to the origin for predicate range, add it back for evaluation).
- For each triangle, append the three pool indices; reverse the winding when `face.reversed`; drop triangles whose 3D area is zero; return `Errc::degenerate` if fewer than one triangle remains.
- Do not touch `sample_curved_face` yet; Task 3 wires the fallback.

- [x] **Step 4: Run the tests to verify they pass**

Same command; expected: all three cases pass.

- [x] **Step 5: Commit**

```powershell
git add include/geom/brep/tessellate.hpp tests/unit/tessellate_trimmed_test.cpp tests/CMakeLists.txt
git commit -m "feat(m3c1): tessellate trimmed curved faces from shared boundary samples"
```

---

### Task 3: Fallback dispatch and adaptive refinement

**Files:**
- Modify: `include/geom/brep/tessellate.hpp` (`sample_curved_face` and the per-face loop in `tessellate()`)
- Modify: `tests/unit/tessellate_trimmed_test.cpp`

**Interfaces:**
- Consumes: `tessellate_trimmed_face` from Task 2.
- Produces: no new public name; `tessellate()` now succeeds on trimmed analytic faces.

- [x] **Step 1: Add the failing refinement tests**

Append to `tests/unit/tessellate_trimmed_test.cpp`:

```cpp
TEST_CASE("trimmed tessellation: smaller deflection refines the trimmed band") {
  const geom::Brep brep = cut_cylinder();
  const auto coarse = geom::tessellate(brep, geom::TessellationOptions{5e-2, 0.35});
  const auto fine = geom::tessellate(brep, geom::TessellationOptions{2e-3, 0.35});
  REQUIRE(coarse.has_value());
  REQUIRE(fine.has_value());
  CHECK(fine.value().triangles.size() > coarse.value().triangles.size());
  check_watertight(coarse.value());
  check_watertight(fine.value());
}

TEST_CASE("trimmed tessellation: chord height of the windowed band meets deflection") {
  const geom::Brep brep = cut_cylinder();
  const double deflection = 2e-3;
  const auto mesh = geom::tessellate(brep, geom::TessellationOptions{deflection, 0.35});
  REQUIRE(mesh.has_value());
  double worst = 0.0;
  for (const auto& triangle : mesh.value().triangles) {
    const geom::Point3d& a = mesh.value().positions[triangle[0]];
    const geom::Point3d& b = mesh.value().positions[triangle[1]];
    const geom::Point3d& c = mesh.value().positions[triangle[2]];
    const geom::Point3d mid{(a[0] + b[0]) / 2.0, (a[1] + b[1]) / 2.0, (a[2] + b[2]) / 2.0};
    if (std::abs(std::hypot(a[0], a[1]) - 1.0) < 1e-9) {
      worst = std::max(worst, std::abs(std::hypot(mid[0], mid[1]) - 1.0));
    }
  }
  CHECK(worst <= deflection * 1.05);
}
```

- [x] **Step 2: Run to verify failure**

Expected: the windowed fixture still fails at the dispatch (or the refinement case fails because the coarse and fine meshes have the same triangle count).

- [x] **Step 3: Implement fallback and refinement**

Implementation notes:
- In `sample_curved_face`, first build into a scratch copy of the state (`TessellateState attempt = state;`). Run the existing analytic builder on the copy; on success assign `state = std::move(attempt)` and return. On failure, discard the copy, obtain the projection with `face_uv_projection`, call `tessellate_trimmed_face` on `state`, and if that fails return the original analytic error (deterministic message).
- Refinement: after lifting, subdivide each triangle uniformly at a per-face integer level. Choose the smallest level `L <= 10` whose midpoint deviation test passes for all triangles: midpoint of the three UV edges evaluated on the surface and compared with the chord midpoint, max distance `<= deflection`; if level 10 still fails, or `4^L * triangle_count` would exceed 2,000,000, return `Errc::numeric_failure` naming the face instead of allocating. Subdivide each triangle into `4^L` children with shared midpoint indices cached by ordered UV-edge key so no T-junction appears.
- Keep the `TessellateState.edges` cache authoritative: refinement only adds interior points.
- Update the doc comment of `tessellate()` to state that trimmed analytic faces use the UV triangulation path, and that `angle_tolerance` is still validated but not used for sampling.

- [x] **Step 4: Run Debug and Release suites**

```powershell
cmake --build build --config Debug --target geom_tests
.\build\tests\Debug\geom_tests.exe
cmake --build build-rel --config Release --target geom_tests
.\build-rel\tests\Release\geom_tests.exe
```
Expected: all cases pass, including the existing `tessellate_curved_test.cpp` cases unchanged.

- [x] **Step 5: Commit**

```powershell
git add include/geom/brep/tessellate.hpp tests/unit/tessellate_trimmed_test.cpp
git commit -m "feat(m3c1): refine trimmed meshes and fall back for trimmed analytic faces"
```

---

### Task 4: Mass properties for trimmed curved faces

**Files:**
- Modify: `include/geom/brep/mass_properties.hpp`
- Modify: `tests/unit/mass_properties_test.cpp`
- Modify: `include/geom/brep/trimmed_domain.hpp` if Task 1's resampling hook needs the deflection (no API change expected)

**Interfaces:**
- Consumes: `trimmed_face_boundaries` (with `projection.tolerance` and `options.curved_sampling_deflection`), `triangulate_trimmed_domain`, `surface_eval`.
- Produces:

```cpp
struct MassPropertiesOptions {
  double sampling_deflection = 1e-9;
  int quadrature_order = 24;
  double curved_sampling_deflection = 1e-5;  // boundary chord for trimmed curved faces
};
```

- [x] **Step 1: Write the failing tests**

Append to `tests/unit/mass_properties_test.cpp`:

```cpp
TEST_CASE("mass properties: a half-trimmed cylinder matches its analytic moments") {
  const geom::Brep cylinder =
      geom::make_cylinder(geom::Point3d{0.0, 0.0, 0.0}, geom::Vec3d{0.0, 0.0, 1.0}, 1.0, 2.0,
                          1e-9)
          .value();
  const geom::Brep box =
      geom::make_box(geom::Point3d{0.0, -2.0, 0.5}, geom::Point3d{2.0, 2.0, 1.5}, 1e-9).value();
  const geom::Brep brep =
      geom::boolean_solid_curved(cylinder, cylinder.solid_ids().front(), box,
                                 box.solid_ids().front(), geom::SolidBooleanOp::difference,
                                 geom::CurvedBooleanOptions{1e-9})
          .value()
          .brep;

  const geom::Result<geom::MassProperties> mass = geom::mass_properties(
      brep, geom::MassPropertiesOptions{1e-9, 24, 1e-5});
  REQUIRE(mass.has_value());
  constexpr double kPi = std::numbers::pi_v<double>;
  // Cylinder minus the half-slab {x >= 0, 0.5 <= z <= 1.5}: volume 3*pi/2,
  // centroid x = -4/(9*pi^2); inertia about the centroid.
  CHECK(mass.value().volume == doctest::Approx(1.5 * kPi).epsilon(1e-5));
  CHECK(mass.value().centroid[0] == doctest::Approx(-4.0 / (9.0 * kPi * kPi)).epsilon(1e-5));
  CHECK(mass.value().centroid[1] == doctest::Approx(0.0).epsilon(1e-9));
  CHECK(mass.value().centroid[2] == doctest::Approx(1.0).epsilon(1e-9));
  CHECK(mass.value().inertia[0][0] == doctest::Approx(kPi).epsilon(1e-4));
  CHECK(mass.value().inertia[1][1] ==
        doctest::Approx(kPi - 8.0 / (27.0 * kPi * kPi * kPi)).epsilon(1e-4));
  CHECK(mass.value().inertia[2][2] ==
        doctest::Approx(0.75 * kPi - 8.0 / (27.0 * kPi * kPi * kPi)).epsilon(1e-4));
}

TEST_CASE("mass properties: trimmed volume agrees with the volume integrator") {
  const geom::Brep cylinder =
      geom::make_cylinder(geom::Point3d{0.0, 0.0, 0.0}, geom::Vec3d{0.0, 0.0, 1.0}, 1.0, 2.0,
                          1e-9)
          .value();
  const geom::Brep box =
      geom::make_box(geom::Point3d{0.0, -2.0, 0.5}, geom::Point3d{2.0, 2.0, 1.5}, 1e-9).value();
  const geom::Brep brep =
      geom::boolean_solid_curved(cylinder, cylinder.solid_ids().front(), box,
                                 box.solid_ids().front(), geom::SolidBooleanOp::difference,
                                 geom::CurvedBooleanOptions{1e-9})
          .value()
          .brep;
  const geom::Result<double> reference = geom::detail::curved_solid_volume(
      brep, brep.solid_ids().front(), 1e-9);
  REQUIRE(reference.has_value());
  const auto mass = geom::mass_properties(brep, geom::MassPropertiesOptions{1e-9, 24, 1e-5});
  REQUIRE(mass.has_value());
  CHECK(mass.value().volume == doctest::Approx(reference.value()).epsilon(1e-5));
}
```

Both cases need `#include <geom/brep/boolean_solid_curved.hpp>`, `<geom/brep/primitives.hpp>` and `<geom/brep/detail/curved_volume.hpp>` at the top of `mass_properties_test.cpp` if they are not already there.

- [x] **Step 2: Run to verify failure**

Expected: `Errc::unsupported` "curved face is not an untrimmed full rectangle".

- [x] **Step 3: Implement trimmed integration**

Implementation notes:
- In `mass_add_face`, keep the existing untrimmed path. When it would return `mass_unsupported_rectangle`, build the face's trimmed boundaries with `trimmed_face_boundaries(brep, id, projection, options.curved_sampling_deflection)` and `triangulate_trimmed_domain`; if either fails return that error.
- For each UV triangle, evaluate the surface at the seven Dunavant points of `kCurvedVolumeQuadraturePoints`, compute `area = cross(du, dv)` and `point`, then call `mass_add_surface_element` with `weight * triangle_uv_area * rule_weight` exactly as `mass_integrate_parametric` does; flip the sign for `face.reversed`.
- Keep the planar path (exact polygon moments) untouched; only curved faces take the trimmed quadrature.
- If a face's rings are not representable (multiple outer loops, non-finite, degenerate), return `Errc::invalid_topology`/`degenerate` naming the face; do not fall back to the old unsupported message.
- Document in the header that trimmed curved faces are integrated by UV triangulation with `curved_sampling_deflection` (default 1e-5, relative accuracy ~1e-5) and that the untrimmed rectangle path remains exact.

- [x] **Step 4: Run the suites**

Expected: the two new cases pass and every existing mass-properties case is unchanged.

- [x] **Step 5: Commit**

```powershell
git add include/geom/brep/mass_properties.hpp tests/unit/mass_properties_test.cpp
git commit -m "feat(m3c1): integrate mass properties over trimmed curved faces"
```

---

### Task 5: QA-13, benchmark, ADR and architecture sync

**Files:**
- Modify: `tests/unit/tessellate_trimmed_test.cpp` (determinism and corpus case)
- Modify: `include/geom/brep/tessellate.hpp` (doc comment only, if not done)
- Create: `docs/benchmarks/m3c1.md`
- Create: `docs/adr/0027-trimmed-tessellation.md`
- Modify: `docs/architecture/geometry-library-architecture.md` §16.3 (M3c-1 status and limits)
- Modify: `docs/superpowers/plans/2026-11-13-m3c1-trimmed-tessellation.md` (check the M3c-1 exit boxes)

**Interfaces:** no code interface; the deliverable is evidence and documentation.

- [x] **Step 1: Add the QA-13 test**

```cpp
TEST_CASE("QA-13: trimmed tessellation is watertight, finite and deterministic") {
  const geom::Brep brep = cut_cylinder();
  const auto first = geom::tessellate(brep, geom::TessellationOptions{2e-3, 0.35});
  const auto second = geom::tessellate(brep, geom::TessellationOptions{2e-3, 0.35});
  REQUIRE(first.has_value());
  REQUIRE(second.has_value());
  check_watertight(first.value());
  CHECK(first.value().triangles == second.value().triangles);
  CHECK(first.value().positions.size() == second.value().positions.size());
  for (std::size_t i = 0; i < first.value().positions.size(); ++i) {
    CHECK(first.value().positions[i][0] == second.value().positions[i][0]);
    CHECK(first.value().positions[i][1] == second.value().positions[i][1]);
    CHECK(first.value().positions[i][2] == second.value().positions[i][2]);
    CHECK(geom::is_finite(first.value().positions[i][0]));
  }
}
```

Add the matching `tessellate` include and a `is_finite` helper as needed.

- [x] **Step 2: Run and record the numbers**

Run the Debug and Release suites and record: total cases/assertions, the trimmed-cylinder triangle count at `deflection = 2e-3`, the mesh-volume error against `3*pi/2`, and the trimmed mass-properties errors against the analytic values. Put them in `docs/benchmarks/m3c1.md` with the machine/compiler line, and in the ADR evidence section.

- [x] **Step 3: Write ADR-0027**

Record: trimmed analytic faces are meshed by UV triangulation of the projected ring set (not by extending the per-surface grids); boundary indices are shared with neighbouring faces; refinement is uniform per face; mass properties integrate the same triangulation with Dunavant quadrature at `curved_sampling_deflection = 1e-5` (relative accuracy ~1e-5); sphere/torus/cone operands and NURBS trim contours are M3c-2; canonical primitives keep the old grids.

- [x] **Step 4: Update §16.3 and the plan boxes**

Add an M3c-1 row: "trimmed cylinder tessellation + trimmed mass properties: PASS (QA-13, QA-11 extended); limits: M3c-2 for sphere/torus/cone operands and NURBS". Check the task boxes in this plan.

- [x] **Step 5: Commit**

```powershell
git add docs tests/unit/tessellate_trimmed_test.cpp
git commit -m "docs(m3c1): record trimmed tessellation evidence and limits"
```

---

## M3c-1 Exit Criteria

- [x] `tessellate()` returns a watertight, finite, per-face-grouped mesh for the seam-crossing window fixture and every canonical primitive test stays green.
- [x] The trimmed band meets the chord-height bound at `deflection = 2e-3`; a smaller deflection produces more triangles.
- [x] `mass_properties()` returns volume/centroid/inertia for the half-trimmed cylinder within the documented 1e-5 relative accuracy, matching both the analytic values and `detail::curved_solid_volume`.
- [x] QA-13 exists and is deterministic; the ADR and `docs/benchmarks/m3c1.md` record the measured evidence.
- [x] Sphere, torus, cone and NURBS remain explicit `unsupported`/deferred and are listed for M3c-2.

## Review Focus Notes for the Implementer

- Do not sample boundary curves twice with different parameterizations: the shared index path is what makes the mesh watertight.
- Do not merge vertices by tolerance at the end; that hides a T-junction instead of fixing it.
- When a ring crosses the periodic seam, keep the unwrapped u (it may exceed `2π`); translating it back would tear the ring.
- Refinement depth is per face, not per triangle, until you have a conforming red-green scheme; over-refining is acceptable, cracks are not.

## Self-Review

- Spec coverage: §8.10 watertightness (Tasks 2–3), deviation (Task 3), coverage/holes (Tasks 1–2); §8.11 divergence-theorem moments (Task 4); QA (Task 5). Sphere/torus/cone/NURBS are explicitly out of this plan and listed for M3c-2.
- No placeholders: every test and signature is complete; the two Task 4 cases build the same fixture explicitly so each case stays runnable on its own.
- Type consistency: `TrimmedBoundary`, `TrimmedUvTriangulation`, `tessellate_trimmed_face`, `MassPropertiesOptions.curved_sampling_deflection` are used with the same names and types in every task that mentions them.
- Review Focus: each of the five lines has a test — hole validity (Task 1 hole-ring and duplicate/non-finite cases plus Task 2 watertight), seam (Task 1 ring `|Δu| < pi` check on the seam-crossing window), shared boundary (Task 2 watertight and canonical case), sample cap (Task 1 `too_fine` numeric-failure case), double counting (Task 4 agreement with `curved_solid_volume`).
