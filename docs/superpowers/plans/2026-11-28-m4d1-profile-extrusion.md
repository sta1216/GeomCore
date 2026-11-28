# M4d-1 Planar Profile Extrusion Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use subagent-driven-development (recommended) or executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Add a B-rep `extrude_profile` operator that turns a planar polygon with holes into a closed, exact, watertight solid extruded along the sketch plane normal.

**Architecture:** The sketch plane and profile are validated, the profile rings are normalized to outer-CCW / holes-CW, and the solid is built directly from 2 (or 2 + n) planar caps and one planar quad per profile edge. `assemble_solid` welds the shared vertices and validates the topology; no boolean, no SSI, no faceting, no sampling. The cap faces carry their hole loops directly, so the result is exact for arbitrary polygonal profiles.

**Tech Stack:** C++20, existing `geom::v1` headers, doctest 2.4.11, CMake, MSVC `/W4 /WX /utf-8` and GCC/Clang warning parity.

**Spec:** `docs/architecture/geometry-library-architecture.md` §8.7 (B-rep modeling operators: extrude/revolve/loft/sweep) and §16.4 (milestone evidence). Also `docs/quality/testing-strategy.md`.

## Global Constraints

- C++20, no new third-party dependency; only existing headers under `include/geom/`.
- Build with `-DGEOM_WERROR=ON`; zero warnings on MSVC `/W4 /WX /utf-8` and on the existing CI compilers.
- No global epsilon. Geometric orientation decisions use the exact predicates already in `geom/numeric/`; tolerances are explicit arguments and never silent defaults.
- A failed call returns a stable `Errc` + message and must not modify the caller's profile or plane.
- Deterministic: the same inputs produce the same counts, volume, bounds and tessellation on every run.
- Topology changes are published only through `assemble_solid`; no partially built solid is returned.
- Every new test file is registered in `tests/CMakeLists.txt` under `GEOM_TEST_SOURCES`.
- One task = one commit, with RED evidence before the implementation.

## Scope

M4d-1 supports:

- one planar sketch plane (`Plane3`), mapped through `make_plane_frame`;
- one simple outer ring plus zero or more strictly interior, non-touching hole rings (`Polygon2d` semantics);
- a finite positive extrusion height along the plane's unit normal;
- exact planar caps with hole loops and exact planar side quads.

M4d-2 (later) will add oblique extrusion/shear, non-planar paths, revolve, loft and sweep with frame propagation; M4d-1 deliberately does not accept a direction vector or a negative height (flip the plane instead). Curved profiles are out of scope: a `Polygon2d` is polygonal, so every face and edge is straight.

## Review Focus

1. **Hole containment and orientation** — a hole outside, crossing or touching the outer ring must be rejected as `degenerate`, never assembled into a self-intersecting cap. Tested in Task 2 and Task 3.
2. **Self-intersecting profile** — a bowtie outer ring must be rejected, not triangulated into overlapping caps. Tested in Task 3.
3. **Tilted/non-origin plane** — the plane frame must map the profile without changing area or orientation; volume and outward normals must remain correct. Tested in Task 2.
4. **Repeated or collinear vertices** — consecutive duplicates must be rejected (they would create zero-length edges); a collinear run must either be accepted with a valid ring or rejected, never produce a non-manifold edge. Tested in Task 3.
5. **Scale robustness** — a profile at `1e6` and a height at `1e-6` must extrude with finite volume within tolerance; tolerance `1e-12` must not merge distinct ring vertices. Tested in Task 3.

---

### Task 1: Extrusion API and outer-ring solid

**Files:**
- Create: `include/geom/brep/extrude.hpp`
- Create: `tests/unit/extrude_test.cpp`
- Modify: `tests/CMakeLists.txt` (append `unit/extrude_test.cpp` to `GEOM_TEST_SOURCES`)

**Interfaces:**
- Consumes: `Plane3` (`geom/surfaces/plane.hpp`), `make_plane_frame`/`plane_frame_to_3d` (`geom/brep/plane_frame.hpp`), `Polygon2d`/`validate_polygon2`/`ring_orientation2`/`polygon_area2` (`geom/cgeom/polygon2.hpp`), `AssembledFace`/`assemble_solid` (`geom/brep/assemble.hpp`), `Brep` (`geom/brep/topology.hpp`), `Result`/`Errc`.
- Produces:

```cpp
Result<Brep> extrude_profile(const Plane3& plane,
                             const Polygon2d& profile,
                             double height,
                             double tolerance);
```

- Preconditions and errors, checked in this order:
  1. non-finite tolerance or `tolerance <= 0` -> `Errc::tolerance_violation` `"extrude_profile: tolerance must be finite and positive"`;
  2. invalid plane (non-finite normal/offset, zero normal) -> propagate `make_plane_frame`'s error unchanged;
  3. non-finite height -> `Errc::non_finite` `"extrude_profile: height is not finite"`; `height <= 0` -> `Errc::degenerate` `"extrude_profile: height must be positive"`;
  4. `validate_polygon2(profile)` -> propagate its `non_finite`/`degenerate`;
  5. `profile.ring_count() != 1` in Task 1 -> `Errc::unsupported` `"extrude_profile: hole rings are not supported yet"` (removed in Task 2);
  6. outer orientation `ring_orientation2(profile.ring(0)) != 1` in Task 1: if `0` -> `Errc::degenerate` `"extrude_profile: profile ring is degenerate"`, otherwise reverse the ring mentally (Task 1 may accept both orientations; the test uses CCW).
  7. consecutive duplicate points, including last==first -> `Errc::degenerate` `"extrude_profile: profile ring has a repeated vertex"`.
- Construction: bottom cap at `plane_frame_to_3d(frame, q)` with outward normal `-frame.normal`, top cap at `point + frame.normal * height` with outward normal `+frame.normal`, and one side quad `[a, b, b + n*h, a + n*h]` per consecutive ring pair. The side quad's outward normal is `normalize(cross(b - a, frame.normal))`. Bottom rings are the top rings reversed so that every cap's outer ring is CCW about its own outward normal and holes are CW.
- On assembly failure, return `Errc::invalid_topology` `"extrude_profile: solid could not be assembled"` (do not leak the assembler's internal message).

- [ ] **Step 1: Write the failing test**

Add `tests/unit/extrude_test.cpp` with:

```cpp
#include <doctest/doctest.h>

#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <map>
#include <utility>
#include <vector>

#include <geom/brep/extrude.hpp>
#include <geom/brep/mass_properties.hpp>
#include <geom/brep/tessellate.hpp>
#include <geom/brep/topology.hpp>
#include <geom/brep/validate.hpp>
#include <geom/cgeom/polygon2.hpp>
#include <geom/core/result.hpp>
#include <geom/primitives/point.hpp>
#include <geom/surfaces/plane.hpp>

namespace {

using namespace geom;

Polygon2d square_profile() {
  Polygon2d profile;
  profile.vertices = {Point2d{0.0, 0.0}, Point2d{2.0, 0.0}, Point2d{2.0, 2.0},
                      Point2d{0.0, 2.0}};
  profile.ring_starts = {0u, 4u};
  return profile;
}

double volume_of(const Brep& brep) {
  const Result<MassProperties> props = mass_properties(brep);
  REQUIRE(props.has_value());
  return props->volume;
}

std::size_t boundary_edges(const TriMesh& mesh) {
  std::map<std::pair<std::uint32_t, std::uint32_t>, int> uses;
  for (const std::array<std::uint32_t, 3>& triangle : mesh.triangles) {
    for (int i = 0; i < 3; ++i) {
      std::uint32_t a = triangle[static_cast<std::size_t>(i)];
      std::uint32_t b = triangle[static_cast<std::size_t>((i + 1) % 3)];
      if (a > b) std::swap(a, b);
      ++uses[{a, b}];
    }
  }
  std::size_t count = 0;
  for (const auto& entry : uses) {
    if (entry.second == 1) ++count;
  }
  return count;
}

}  // namespace

TEST_CASE("extrude a square profile") {
  const Polygon2d square = square_profile();
  const Plane3 plane{geom::Vec3d{0.0, 0.0, 1.0}, 0.0};
  const Result<Brep> solid = extrude_profile(plane, square, 3.0, 1e-9);
  REQUIRE(solid.has_value());
  CHECK(solid->solid_count() == 1u);
  CHECK(solid->shell_count() == 1u);
  CHECK(solid->vertex_count() == 8u);
  CHECK(solid->edge_count() == 12u);
  CHECK(solid->face_count() == 6u);
  CHECK(volume_of(solid.value()) == doctest::Approx(12.0).epsilon(1e-12));
  REQUIRE(validate_topology(solid.value()).ok());

  const Result<TriMesh> mesh = tessellate(solid.value(), TessellationOptions{});
  REQUIRE(mesh.has_value());
  CHECK(boundary_edges(mesh.value()) == 0u);
  CHECK(!mesh.value().triangles.empty());
}

TEST_CASE("extrude rejects invalid options and height") {
  const Polygon2d square = square_profile();
  const Plane3 plane{geom::Vec3d{0.0, 0.0, 1.0}, 0.0};

  CHECK(extrude_profile(plane, square, 1.0, 0.0).error().code ==
        Errc::tolerance_violation);
  CHECK(extrude_profile(plane, square, 1.0, -1e-9).error().code ==
        Errc::tolerance_violation);
  CHECK(extrude_profile(plane, square, 1.0, std::numeric_limits<double>::quiet_NaN())
            .error().code == Errc::tolerance_violation);
  CHECK(extrude_profile(plane, square, 0.0, 1e-9).error().code == Errc::degenerate);
  CHECK(extrude_profile(plane, square, -1.0, 1e-9).error().code == Errc::degenerate);
  CHECK(extrude_profile(plane, square, std::numeric_limits<double>::infinity(), 1e-9)
            .error().code == Errc::non_finite);
}

TEST_CASE("extrude rejects a structurally invalid profile") {
  const Plane3 plane{geom::Vec3d{0.0, 0.0, 1.0}, 0.0};
  const Polygon2d empty;
  CHECK(extrude_profile(plane, empty, 1.0, 1e-9).error().code ==
        Errc::degenerate);

  Polygon2d duplicate = square_profile();
  duplicate.vertices.insert(duplicate.vertices.begin(), Point2d{0.0, 0.0});
  duplicate.ring_starts = {0u, 5u};
  CHECK(extrude_profile(plane, duplicate, 1.0, 1e-9).error().code ==
        Errc::degenerate);
}
```

- [ ] **Step 2: Run the test to verify it fails**

Run (Windows, from the repository root):

```powershell
cmake -S . -B build -DGEOM_WERROR=ON
cmake --build build --config Debug --target geom_tests
.\build\tests\Debug\geom_tests.exe --test-case="*extrude*"
```

Expected: compile fails with `C1083: Cannot open include file: 'geom/brep/extrude.hpp'`.

- [ ] **Step 3: Implement `extrude_profile`**

Create `include/geom/brep/extrude.hpp` with the public signature above. Build a `std::vector<AssembledFace>`; map each ring point with `plane_frame_to_3d`; add the bottom cap first (outer ring reversed), then the top cap, then the side quads. Call `assemble_solid(faces, tolerance)` and map its failure to the `invalid_topology` message above. Include `<algorithm>`, `<cmath>`, `<span>`, `<vector>`, and the headers listed in Interfaces.

Keep the ring orientation check exact: use `ring_orientation2`; never `signed_area2` for the accept/reject decision. Only Task 1 rejects holes; Task 2 replaces that rejection.

- [ ] **Step 4: Run the test to verify it passes**

```powershell
cmake --build build --config Debug --target geom_tests
.\build\tests\Debug\geom_tests.exe --test-case="*extrude*"
```

Expected: 3 cases, 0 failures. Then repeat for `--config Release`.

- [ ] **Step 5: Commit**

```powershell
git add include/geom/brep/extrude.hpp tests/unit/extrude_test.cpp tests/CMakeLists.txt
git commit -m "feat(brep): extrude a planar outer profile"
```

---

### Task 2: Hole caps, arbitrary plane frames and orientation normalization

**Files:**
- Modify: `include/geom/brep/extrude.hpp`
- Modify: `tests/unit/extrude_test.cpp`

**Interfaces:**
- Consumes: Task 1's `extrude_profile`.
- Produces: the same public signature, now accepting `ring_count() >= 1`; outer ring normalized to CCW and every hole to CW before cap construction. `AssembledFace::holes` receives the hole rings; the bottom cap reverses every ring, the top cap keeps the normalized order. Side faces are built for outer and hole rings alike, with outward normal `normalize(cross(edge, frame.normal))`.

- [ ] **Step 1: Write the failing tests**

Append to `tests/unit/extrude_test.cpp`:

```cpp
TEST_CASE("extrude a profile with a hole") {
  Polygon2d profile;
  profile.vertices = {
      Point2d{0.0, 0.0}, Point2d{4.0, 0.0}, Point2d{4.0, 4.0}, Point2d{0.0, 4.0},
      Point2d{1.0, 1.0}, Point2d{1.0, 3.0}, Point2d{3.0, 3.0}, Point2d{3.0, 1.0}};
  profile.ring_starts = {0u, 4u, 8u};
  const Plane3 plane{geom::Vec3d{0.0, 0.0, 1.0}, 0.0};

  const Result<Brep> solid = extrude_profile(plane, profile, 2.0, 1e-9);
  REQUIRE(solid.has_value());
  CHECK(solid->vertex_count() == 16u);
  CHECK(solid->edge_count() == 24u);
  CHECK(solid->face_count() == 10u);
  CHECK(volume_of(solid.value()) == doctest::Approx(24.0).epsilon(1e-12));
  REQUIRE(validate_topology(solid.value()).ok());
}

TEST_CASE("extrude ignores the input ring orientations") {
  Polygon2d profile;
  profile.vertices = {
      Point2d{0.0, 4.0}, Point2d{4.0, 4.0}, Point2d{4.0, 0.0}, Point2d{0.0, 0.0},
      Point2d{1.0, 1.0}, Point2d{3.0, 1.0}, Point2d{3.0, 3.0}, Point2d{1.0, 3.0}};
  profile.ring_starts = {0u, 4u, 8u};
  const Plane3 plane{geom::Vec3d{0.0, 0.0, 1.0}, 0.0};

  const Result<Brep> solid = extrude_profile(plane, profile, 2.0, 1e-9);
  REQUIRE(solid.has_value());
  CHECK(volume_of(solid.value()) == doctest::Approx(24.0).epsilon(1e-12));
  REQUIRE(validate_topology(solid.value()).ok());
}

TEST_CASE("extrude an L-shaped profile") {
  Polygon2d profile;
  profile.vertices = {Point2d{0.0, 0.0}, Point2d{3.0, 0.0}, Point2d{3.0, 1.0},
                      Point2d{1.0, 1.0}, Point2d{1.0, 2.0}, Point2d{0.0, 2.0}};
  profile.ring_starts = {0u, 6u};
  const Plane3 plane{geom::Vec3d{0.0, 0.0, 1.0}, 0.0};

  const Result<Brep> solid = extrude_profile(plane, profile, 2.0, 1e-9);
  REQUIRE(solid.has_value());
  CHECK(solid->vertex_count() == 12u);
  CHECK(solid->edge_count() == 18u);
  CHECK(solid->face_count() == 8u);
  CHECK(volume_of(solid.value()) == doctest::Approx(8.0).epsilon(1e-12));
  REQUIRE(validate_topology(solid.value()).ok());
}

TEST_CASE("extrude on a tilted, offset plane") {
  Polygon2d profile;
  profile.vertices = {Point2d{0.0, 0.0}, Point2d{2.0, 0.0}, Point2d{2.0, 2.0},
                      Point2d{0.0, 2.0}};
  profile.ring_starts = {0u, 4u};
  const Plane3 plane{geom::Vec3d{1.0, 1.0, 1.0}, 0.5};

  const Result<Brep> solid = extrude_profile(plane, profile, 2.0, 1e-9);
  REQUIRE(solid.has_value());
  CHECK(volume_of(solid.value()) == doctest::Approx(8.0).epsilon(1e-12));
  REQUIRE(validate_topology(solid.value()).ok());
}
```

Add `<limits>` if the file does not already include it.

- [ ] **Step 2: Run the tests to verify they fail**

```powershell
cmake --build build --config Debug --target geom_tests
.\build\tests\Debug\geom_tests.exe --test-case="*extrude*"
```

Expected: the hole test fails at `REQUIRE(solid.has_value())` with
`unsupported` / `"extrude_profile: hole rings are not supported yet"`; the
L-shape and tilted-plane tests fail for the same missing implementation.

- [ ] **Step 3: Implement hole caps and normalization**

Replace the Task 1 single-ring restriction with a normalized ring list:

- for each ring, compute `ring_orientation2`; `0` -> `degenerate`;
- outer must end CCW (`== 1`), holes must end CW (`== -1`); reverse a ring's vertex order when needed;
- reject consecutive duplicates in every ring, including last==first;
- build `AssembledFace` caps: `outer` = outer ring, `holes` = hole rings; bottom cap uses the reverse of each normalized ring and normal `-n`, top cap uses the normalized rings and normal `+n`;
- build side quads for every ring (outer and holes), normal `normalize(cross(b - a, n))`.

Do not call `triangulate_polygon2` yet; Task 3 adds the simplicity check.

- [ ] **Step 4: Run the tests to verify they pass**

```powershell
cmake --build build --config Debug --target geom_tests
.\build\tests\Debug\geom_tests.exe --test-case="*extrude*"
```

Expected: 7 cases, 0 failures in Debug and Release.

- [ ] **Step 5: Commit**

```powershell
git add include/geom/brep/extrude.hpp tests/unit/extrude_test.cpp
git commit -m "feat(brep): extrude planar profiles with holes"
```

---

### Task 3: Simplicity certification, diagnostics and robustness

**Files:**
- Modify: `include/geom/brep/extrude.hpp`
- Modify: `tests/unit/extrude_test.cpp`

**Interfaces:**
- Consumes: `triangulate_polygon2(const Polygon2d&) -> Triangulation2` (`geom/cgeom/triangulate2.hpp`); Task 2's normalized rings.
- Produces: the same `extrude_profile` signature; a profile whose triangulation is not `ok` is rejected as `Errc::degenerate` `"extrude_profile: profile is not a simple polygon with holes"` before any face is built. All other errors are unchanged.

- [ ] **Step 1: Write the failing tests**

```cpp
TEST_CASE("extrude rejects a self-intersecting profile") {
  Polygon2d bowtie;
  bowtie.vertices = {Point2d{0.0, 0.0}, Point2d{2.0, 2.0}, Point2d{0.0, 2.0},
                     Point2d{2.0, 0.0}};
  bowtie.ring_starts = {0u, 4u};
  const Plane3 plane{geom::Vec3d{0.0, 0.0, 1.0}, 0.0};
  CHECK(extrude_profile(plane, bowtie, 1.0, 1e-9).error().code ==
        Errc::degenerate);
}

TEST_CASE("extrude rejects a hole outside the outer ring") {
  Polygon2d profile;
  profile.vertices = {Point2d{0.0, 0.0}, Point2d{2.0, 0.0}, Point2d{2.0, 2.0},
                      Point2d{0.0, 2.0}, Point2d{3.0, 3.0}, Point2d{4.0, 3.0},
                      Point2d{4.0, 4.0}, Point2d{3.0, 4.0}};
  profile.ring_starts = {0u, 4u, 8u};
  const Plane3 plane{geom::Vec3d{0.0, 0.0, 1.0}, 0.0};
  CHECK(extrude_profile(plane, profile, 1.0, 1e-9).error().code ==
        Errc::degenerate);
}

TEST_CASE("extrude is deterministic and leaves its inputs unchanged") {
  const Polygon2d square = square_profile();
  const Plane3 plane{geom::Vec3d{0.0, 0.0, 1.0}, 0.0};
  const std::vector<Point2d> before_vertices = square.vertices;
  const std::vector<std::uint32_t> before_starts = square.ring_starts;

  const Result<Brep> first = extrude_profile(plane, square, 2.5, 1e-9);
  const Result<Brep> second = extrude_profile(plane, square, 2.5, 1e-9);
  REQUIRE(first.has_value());
  REQUIRE(second.has_value());
  CHECK(first->vertex_count() == second->vertex_count());
  CHECK(first->edge_count() == second->edge_count());
  CHECK(first->face_count() == second->face_count());
  CHECK(volume_of(first.value()) == volume_of(second.value()));
  CHECK(square.vertices.size() == before_vertices.size());
  CHECK(square.ring_starts == before_starts);
  for (std::size_t i = 0; i < before_vertices.size(); ++i) {
    CHECK(square.vertices[i][0] == before_vertices[i][0]);
    CHECK(square.vertices[i][1] == before_vertices[i][1]);
  }
}

TEST_CASE("extrude keeps a small height on a large profile") {
  Polygon2d square;
  square.vertices = {Point2d{0.0, 0.0}, Point2d{1e6, 0.0}, Point2d{1e6, 1e6},
                     Point2d{0.0, 1e6}};
  square.ring_starts = {0u, 4u};
  const Plane3 plane{geom::Vec3d{0.0, 0.0, 1.0}, 0.0};
  const Result<Brep> solid = extrude_profile(plane, square, 1e-6, 1e-12);
  REQUIRE(solid.has_value());
  CHECK(volume_of(solid.value()) == doctest::Approx(1.0).epsilon(1e-9));
  REQUIRE(validate_topology(solid.value()).ok());
}
```

- [ ] **Step 2: Run the tests to verify they fail**

```powershell
cmake --build build --config Debug --target geom_tests
.\build\tests\Debug\geom_tests.exe --test-case="*extrude*"
```

Expected: the bowtie and outside-hole tests currently assemble a face list
without rejection (or fail with the assembler's `invalid_topology`), not the
required `degenerate`; the determinism/scale tests pass already. Only the two
new rejection tests must fail.

- [ ] **Step 3: Certify the profile and map assembly failures**

After normalization and before building faces, call
`triangulate_polygon2(profile)`; `!ok` -> `degenerate` with the message above.
Keep the existing consecutive-duplicate rejection and the exact orientation
checks. Do not use triangulation output for the caps: the caps remain the
polygon-with-holes rings; the triangulation is only a simplicity certificate.

- [ ] **Step 4: Run the tests to verify they pass**

```powershell
cmake --build build --config Debug --target geom_tests
.\build\tests\Debug\geom_tests.exe --test-case="*extrude*"
```

Expected: 11 cases, 0 failures in Debug and Release.

- [ ] **Step 5: Commit**

```powershell
git add include/geom/brep/extrude.hpp tests/unit/extrude_test.cpp
git commit -m "fix(brep): reject non-simple and out-of-place extrusion profiles"
```

---

### Task 4: Benchmark, ADR and milestone evidence

**Files:**
- Create: `bench/bench_m4d1.cpp`
- Modify: `CMakeLists.txt` (add `geom_bench_m4d1` under `GEOM_BUILD_BENCH`)
- Modify: `.github/workflows/ci.yml` (build/run/upload `bench-m4d1.txt`)
- Create: `docs/adr/0045-profile-extrusion.md`
- Create: `docs/benchmarks/m4d1-extrusion.md`
- Modify: `docs/architecture/geometry-library-architecture.md` (§8.7 and §16.4)
- Modify: this plan's Results table

**Interfaces:**
- Consumes: Task 1–3's `extrude_profile`.
- Produces: a Release benchmark `geom_bench_m4d1` that builds and tessellates a 1024-gon extrusion, checks the analytic volume and topology, prints `min/median/max`, and exits `2` on any mismatch.

- [ ] **Step 1: Add the benchmark**

`bench/bench_m4d1.cpp`:

```cpp
#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <numbers>
#include <vector>

#include <geom/brep/extrude.hpp>
#include <geom/brep/mass_properties.hpp>
#include <geom/brep/tessellate.hpp>
#include <geom/brep/topology.hpp>
#include <geom/brep/validate.hpp>
#include <geom/cgeom/polygon2.hpp>
#include <geom/surfaces/plane.hpp>

namespace {

using namespace geom;

Polygon2d regular_polygon(int sides, double radius) {
  Polygon2d profile;
  profile.ring_starts = {0u, static_cast<std::uint32_t>(sides)};
  profile.vertices.reserve(static_cast<std::size_t>(sides));
  for (int i = 0; i < sides; ++i) {
    const double angle = 2.0 * std::numbers::pi * static_cast<double>(i) /
                         static_cast<double>(sides);
    profile.vertices.push_back(
        Point2d{radius * std::cos(angle), radius * std::sin(angle)});
  }
  return profile;
}

}  // namespace

int main() {
  constexpr int kSides = 1024;
  constexpr double kRadius = 10.0;
  constexpr double kHeight = 1.0;
  const double area = 0.5 * kSides * kRadius * kRadius *
                      std::sin(2.0 * std::numbers::pi / kSides);
  const double expected = area * kHeight;

  const Polygon2d profile = regular_polygon(kSides, kRadius);
  const Plane3 plane{Vec3d{0.0, 0.0, 1.0}, 0.0};

  const std::array<double, 5> timings = [&] {
    std::array<double, 5> seconds{};
    for (std::size_t i = 0; i < seconds.size(); ++i) {
      const auto start = std::chrono::steady_clock::now();
      const Result<Brep> solid = extrude_profile(plane, profile, kHeight, 1e-12);
      if (!solid) {
        std::fprintf(stderr, "extrude failed: %s\n", solid.error().message);
        std::exit(2);
      }
      if (solid->vertex_count() != 2048u || solid->edge_count() != 3072u ||
          solid->face_count() != 1026u) {
        std::fprintf(stderr, "topology mismatch\n");
        std::exit(2);
      }
      if (!validate_topology(solid.value()).ok()) {
        std::fprintf(stderr, "validation failed\n");
        std::exit(2);
      }
      const Result<MassProperties> props = mass_properties(solid.value());
      if (!props || std::fabs(props->volume - expected) > 1e-9 * expected) {
        std::fprintf(stderr, "volume mismatch: %.17g expected %.17g\n",
                     props ? props->volume : 0.0, expected);
        std::exit(2);
      }
      const Result<TriMesh> mesh = tessellate(solid.value(), TessellationOptions{});
      if (!mesh || mesh->triangles.empty()) {
        std::fprintf(stderr, "tessellation failed\n");
        std::exit(2);
      }
      const auto stop = std::chrono::steady_clock::now();
      seconds[i] =
          std::chrono::duration<double, std::milli>(stop - start).count();
    }
    return seconds;
  }();

  std::array<double, 5> sorted = timings;
  std::sort(sorted.begin(), sorted.end());
  std::printf("extrude 1024-gon h=1: min %.3f ms median %.3f ms max %.3f ms\n",
              sorted.front(), sorted[2], sorted.back());
  std::printf("volume %.17g expected %.17g\n", expected, expected);
  return 0;
}
```

The volume printed is the analytic value; the per-sample check above is the
actual guard. `std::exit` requires `<cstdlib>` (add it).

- [ ] **Step 2: Wire the benchmark and CI**

In the root `CMakeLists.txt`, under `GEOM_BUILD_BENCH`, add an executable
`geom_bench_m4d1` from `bench/bench_m4d1.cpp` linked to `geom_core` with the
same warnings target as the other benchmarks. In `.github/workflows/ci.yml`,
add `geom_bench_m4d1` to the bench build targets, run it, and upload
`bench-m4d1.txt`; use `set -o pipefail` in the run step.

- [ ] **Step 3: Record the evidence**

Write ADR-0045 with: the exact construction (caps + quads), the orientation
normalization, the simplicity certificate, the scope limits (normal extrusion,
polygonal profiles only; oblique/shear/revolve/sweep deferred to M4d-2), and the
failure contract. Update architecture §8.7 with a short "M4d-1" paragraph and
§16.4 with the measured counts and benchmark. Create
`docs/benchmarks/m4d1-extrusion.md` with environment, workload, timings, volume
and reproduction command. Fill this plan's Results table with the exact measured
numbers.

- [ ] **Step 4: Run the full verification**

```powershell
cmake --build build --config Release --target geom_tests geom_bench_m4d1
.\build\tests\Release\geom_tests.exe
.\build\bench\Release\geom_bench_m4d1.exe
```

Expected: full Release suite green; benchmark exits `0` with the analytic volume
and finite timings. Debug `--test-case="*extrude*"` green. Record the full-suite
counts verbatim; do not claim a full Debug pass if it was not run.

- [ ] **Step 5: Commit**

```powershell
git add bench/bench_m4d1.cpp CMakeLists.txt .github/workflows/ci.yml \
        docs/adr/0045-profile-extrusion.md docs/benchmarks/m4d1-extrusion.md \
        docs/architecture/geometry-library-architecture.md \
        docs/superpowers/plans/2026-11-28-m4d1-profile-extrusion.md
git commit -m "docs(brep): record the profile extrusion contract and benchmark"
```

---

## M4d-1 Exit

- `extrude_profile` accepted for square, L-shape, hole, reversed-ring, tilted-plane and 1024-gon profiles.
- Rejections: bad tolerance, bad height, empty/malformed profile, repeated vertex, self-intersecting ring, hole outside outer ring, degenerate plane.
- Exact volume and topology checked for every accepted profile; zero boundary edges after tessellation; deterministic and input-preserving.
- ADR-0045 and `docs/benchmarks/m4d1-extrusion.md` recorded; architecture §8.7/§16.4 updated; CI bench wired.
- Full Release suite green; Debug extrusion lane green; full Debug explicitly PENDING if not run.

## Results

| Task | Commit | Tests | Notes |
| --- | --- | --- | --- |
| T1 | PENDING | PENDING | outer-ring square extrusion |
| T2 | PENDING | PENDING | holes, reversed input, L-shape, tilted plane |
| T3 | PENDING | PENDING | simplicity/diagnostics/robustness |
| T4 | PENDING | PENDING | benchmark, ADR-0045, docs, CI |

## Self-Check

- **Spec coverage:** §8.7 extrude is covered by Task 1–3; §16.4 evidence by Task 4. Revolve/loft/sweep are explicitly deferred to M4d-2, not silently omitted.
- **Step scan:** every step has one action and a checkable result; code is limited to tests and the benchmark (the algorithm is described, not transcribed).
- **Type consistency:** `extrude_profile(const Plane3&, const Polygon2d&, double, double) -> Result<Brep>` is used identically in all tasks; `AssembledFace` fields match `assemble.hpp`.
- **Review Focus:** hole containment (T2/T3), self-intersection (T3), tilted plane (T2), duplicate vertices (T1/T3), scale (T3) each have a test.
- **Proportion:** 4 tasks, ~11 tests; the plan is shorter than the code it will produce.
