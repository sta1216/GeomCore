# Rounded Box (All-Edge Fillet) Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use subagent-driven-development (recommended) or executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Let `fillet_edges` round all twelve edges of an axis-aligned box with one radius, producing six planar faces, twelve exact cylindrical patches and eight exact spherical corner octants.

**Architecture:** Detect the all-twelve-edges selection before the existing box-family selection, then build the twenty-six `CurvedAssembledFace` fragments directly: inset planar rectangles, quarter cylinders, and spherical corner patches whose three boundary arcs are great-circle quarter arcs. Assemble with `assemble_curved_solid` and validate topology. No boolean, no polyhedral approximation, no new tolerance.

**Tech Stack:** C++20, existing `geom::v1` B-rep, `CurvedAssembledFace`, `CircleCurve3`, `Cylinder3`, `Sphere3`, `Plane3`, `assemble_curved_solid`.

**Spec:** `docs/architecture/geometry-library-architecture.md` §8.9 (blend engine) and §16.4 (M4 evidence); this plan is M4c-3c of the fillet track after M4c-3a (prism/box family fillets) and M4c-3b (convex prism vertical edges).

## Global Constraints

- C++20; MSVC `/W4 /WX /utf-8`; GCC/Clang build with `-Wall -Wextra -Wconversion -Wshadow -Wpedantic -Werror` (CI).
- No new third-party dependency; no global epsilon; all tolerance decisions use the caller's `tolerance` or an explicit multiple of it.
- Only axis-aligned boxes, exactly all twelve box edges selected (duplicates allowed after deduplication), one radius.
- `radius` must be finite and strictly `radius < 0.5 * min(a, b, c)`; otherwise return `Errc::unsupported` with a stable message. No silent partial solid.
- Successful result: one solid, one closed shell, `filleted_edges == 12`, `validate_topology(brep).ok()`.
- Failure must not modify the input box. Deterministic: two calls produce bit-identical vertex positions.
- Existing box-family and prism fillet paths and messages must not change.
- Every task ends with a single commit; run Debug and Release filtered tests plus `ctest`.

## Review Focus

1. **Radius just below half the smallest dimension**: faces shrink toward zero; the operation must reject at the limit rather than assemble a degenerate cell. Test both `0.5*min - 1e-6` (accept) and `0.5*min` (reject).
2. **All eight corners**: signs vary per corner; a construction that only works at the low-low-low corner must fail the symmetry test.
3. **Spherical patch chart**: each octant is a polar cap around its corner direction; the loop must unwrap around the pole without a seam jump, and `face_loops_uv`/`mass_properties` must succeed.
4. **Tessellation and pcurves**: every face must tessellate into a watertight mesh; boundary edge count must be zero and triangles non-empty.
5. **Selection boundary**: 11 edges, a duplicated edge, or 12 edges of a non-box must not be mistaken for a rounded box.

## File Map

- Create `include/geom/brep/detail/fillet_box_rounded.hpp` — detection and construction of the twenty-six faces.
- Modify `include/geom/brep/fillet.hpp` — call the rounded path before `fillet_box_selection`.
- Create `tests/unit/fillet_box_rounded_test.cpp` — construction, symmetry, diagnostics.
- Modify `tests/CMakeLists.txt` — register the new test.
- Create `bench/bench_m4c3c.cpp` — 100 mm cube, r = 1 mm.
- Modify root `CMakeLists.txt` — `geom_bench_m4c3c` target.
- Modify `.github/workflows/ci.yml` — build/run/upload the new benchmark.
- Create `docs/adr/0044-rounded-box-spherical-corners.md`, `docs/benchmarks/m4c3c-rounded-box.md`.
- Modify `docs/architecture/geometry-library-architecture.md` §8.9/§16.4 and this plan's Results table.

---

### Task 1: Detect the all-edges selection and build the rounded box

**Files:**
- Create: `include/geom/brep/detail/fillet_box_rounded.hpp`
- Modify: `include/geom/brep/fillet.hpp`
- Create: `tests/unit/fillet_box_rounded_test.cpp`
- Modify: `tests/CMakeLists.txt`

**Interfaces:**
- Consumes: `detail::FilletBox{double lo[3]; double hi[3];}` from `detail/fillet_box.hpp`; `assemble_curved_solid(std::span<const CurvedAssembledFace>, double)`; `CurvedArcSpan{CircleCurve3 circle, double t0, double t1}`; `Sphere3{center, axis, x_axis, radius}`; `Cylinder3{origin, axis, x_axis, radius}`; `Plane3{normal, offset}`; `face_uv_project_sphere`.
- Produces:
  ```cpp
  namespace geom { inline namespace v1 { namespace detail {
  struct FilletBoxRounded { double lo[3]; double hi[3]; };
  // Returns true only when `edges` (after deduplication) are exactly the twelve
  // live edges of `box`, in any order.
  bool fillet_box_all_edges(const Brep& brep, std::span<const EdgeId> edges,
                            const FilletBox& box, double tolerance);
  // Builds the rounded box. Caller has already checked radius and tolerance.
  Result<Brep> fillet_box_rounded(const FilletBoxRounded& box, double radius,
                                  double tolerance);
  } } }
  ```
- `fillet_edges` gains, immediately after `fillet_box_bounds` succeeds:
  ```cpp
  if (detail::fillet_box_all_edges(brep, edges, box_bounds.value(), options.tolerance)) {
    const Result<Brep> rounded = detail::fillet_box_rounded(
        detail::FilletBoxRounded{{box_bounds.value().lo[0], box_bounds.value().lo[1],
                                  box_bounds.value().lo[2]},
                                 {box_bounds.value().hi[0], box_bounds.value().hi[1],
                                  box_bounds.value().hi[2]}},
        options.radius, options.tolerance);
    if (!rounded.has_value()) return rounded.error();
    return FilletResult{std::move(rounded.value()), 12u};
  }
  ```
  The check must happen before `fillet_box_selection`, whose shared-vertex rejection would otherwise fire on the twelve-edge set.

- [ ] **Step 1: Write the failing test**

```cpp
#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <map>
#include <span>
#include <string_view>
#include <utility>
#include <variant>
#include <vector>

#include <doctest/doctest.h>

#include <geom/brep/detail/curved_volume.hpp>
#include <geom/brep/fillet.hpp>
#include <geom/brep/primitives.hpp>
#include <geom/brep/tessellate.hpp>
#include <geom/brep/topology.hpp>
#include <geom/brep/validate.hpp>
#include <geom/curves/analytic_curve.hpp>
#include <geom/surfaces/analytic.hpp>

namespace {
using geom::Brep;
using geom::EdgeId;
using geom::Point3d;

Brep box_2x3x4() {
  return geom::make_box(Point3d{0.0, 0.0, 0.0}, Point3d{2.0, 3.0, 4.0}, 1e-9).value();
}

std::vector<EdgeId> all_edges(const Brep& brep) {
  return std::vector<EdgeId>(brep.edge_ids().begin(), brep.edge_ids().end());
}

double volume_of(const Brep& brep) {
  return geom::detail::curved_solid_volume(brep, brep.solid_ids().front(), 1e-9).value();
}
}  // namespace

TEST_CASE("rounded box rounds all twelve edges with exact corner patches") {
  const Brep box = box_2x3x4();
  const std::vector<EdgeId> edges = all_edges(box);
  REQUIRE(edges.size() == 12u);

  const auto result = geom::fillet_edges(box, edges, geom::FilletOptions{0.25, 1e-9});
  REQUIRE(result.has_value());
  const Brep& rounded = result.value().brep;
  CHECK(result.value().filleted_edges == 12u);
  CHECK(rounded.solid_count() == 1u);
  CHECK(rounded.shell_count() == 1u);
  CHECK(rounded.vertex_count() == 24u);
  CHECK(rounded.edge_count() == 48u);
  CHECK(rounded.face_count() == 26u);
  CHECK(geom::validate_topology(rounded).ok());
  CHECK(volume_of(rounded) == doctest::Approx(23.53807140332).epsilon(1e-12));

  std::uint32_t planes = 0;
  std::uint32_t cylinders = 0;
  std::uint32_t spheres = 0;
  for (const geom::FaceId id : rounded.face_ids()) {
    const geom::SurfaceGeometry& surface = rounded.surface(rounded.face(id).surface);
    if (std::holds_alternative<geom::Plane3>(surface)) ++planes;
    if (std::holds_alternative<geom::Cylinder3>(surface)) ++cylinders;
    if (std::holds_alternative<geom::Sphere3>(surface)) ++spheres;
  }
  CHECK(planes == 6u);
  CHECK(cylinders == 12u);
  CHECK(spheres == 8u);
}
```

- [ ] **Step 2: Run the test to verify it fails**

Run: `ctest --test-dir build-debug --output-on-failure -R geom_tests` or the built `geom_tests.exe --test-case="rounded box rounds all twelve edges*"`.
Expected: the call returns `Errc::unsupported` (the twelve-edge set reaches `fillet_box_selection`, which rejects shared vertices), or the test fails to compile because `fillet_box_rounded` is not implemented.

- [ ] **Step 3: Implement the construction**

In `detail/fillet_box_rounded.hpp`:

1. `fillet_box_all_edges`: deduplicate `edges` by `EdgeId`; require `unique.size() == 12`; require every id `brep.edge_alive(id)`; build the set of the box's twelve edges by matching each edge's endpoints to two box corners (all six coordinates equal `lo` or `hi` within `max(64*tolerance, 1e-12*scale)`), then require the unique set equals that set.

2. `fillet_box_rounded` builds local helper `make_face(surface, points, spans, uv)` like `fillet_box_detail::make_face`.

3. **Planar faces (6).** For each axis `i` and side `s` (lo/hi), the face is the rectangle on plane `x_i = lo_i` or `hi_i`, with the other two coordinates in `[lo+r, hi-r]`. Use outward normal `-e_i` for lo and `+e_i` for hi, offset `lo_i` / `hi_i` with the sign convention of `Plane3` (`normal*point == offset`). Order the four loop points counter-clockwise about the outward normal. UV may be the two remaining point coordinates in ascending axis order, as `fillet_box_detail::side_face` does for its axis-aligned planes.

4. **Cylinder faces (12).** For each axis `f` and the two cross-axis bits `(a_bit, b_bit)`, the axis line is at the inset coordinates `(a_bit ? hi_a - r : lo_a + r, b_bit ? hi_b - r : lo_b + r)`; the patch spans `f` from `lo_f + r` to `hi_f - r`, length `L = hi_f - lo_f - 2r`. Choose `x_axis = e_a` (or `-e_a` when the inset is at `lo_a`, matching the radial direction toward that face) and `y_axis = cross(e_f, x_axis)`; the quarter arc runs from the radial direction toward the a-face to the radial direction toward the b-face. Use the same point/span/uv layout as `fillet_box_detail::cylinder_face`: `evaluate(t0,0), evaluate(t1,0), evaluate(t1,L), evaluate(t0,L)`, with `CurvedArcSpan` on the first and third segments and uv `(t0,0),(t1,0),(t1,L),(t0,L)`.

5. **Sphere faces (8).** For a corner with bits `(bx,by,bz)`, the sphere center is `(bx ? hi_x-r : lo_x+r, ...)`; the corner direction `d = (bx ? +1 : -1, ...)` normalized is `Sphere3.axis`; choose any unit `x_axis` perpendicular to `d` (e.g. normalize of `cross(d, e_axis)` with the first axis where the cross is non-negligible). In local coordinates `(u,v,w)` measured from the corner into the box, the three boundary points are `A=(0,r,r)`, `B=(r,0,r)`, `C=(r,r,0)`; the outward loop is `C, B, A`. Each segment is a quarter arc centered at the sphere center: `C->B` in plane `v=r`, `B->A` in plane `w=r`, `A->C` in plane `u=r`; each `CurvedArcSpan` has `t0=0, t1=pi/2` and a `CircleCurve3` whose `x_axis` points to the segment start and whose `y_axis` is the quarter-turn direction (so evaluating `t in [0,pi/2]` follows the loop). Map local points with the corner frame (axis directions point into the box). Build `loops_uv` with `face_uv_project_sphere(sphere, point, previous_u, has_previous)`, carrying `previous_u` across the loop and **not wrapping** it; for a closed loop the last projected u may differ from the first by a full turn, which is allowed.

6. Call `assemble_curved_solid(faces, tolerance)`; return its `Result<Brep>` unchanged (the assembler already validates topology). If any face has fewer than three loop points or a non-finite coordinate, return `Errc::invalid_topology`/`non_finite` with a `fillet_edges: ...` prefix rather than letting the assembler report an anonymous error.

- [ ] **Step 4: Run the test to verify it passes**

Run the filtered test in Debug and Release. Expected: 1 case, 14 assertions, volume `23.53807140332` within `1e-12`.

- [ ] **Step 5: Commit**

```bash
git add include/geom/brep/detail/fillet_box_rounded.hpp include/geom/brep/fillet.hpp tests/unit/fillet_box_rounded_test.cpp tests/CMakeLists.txt
git commit -m "feat(brep): round all twelve edges of a box with spherical corners"
```

---

### Task 2: Symmetry, watertightness, determinism and input immutability

**Files:**
- Modify: `tests/unit/fillet_box_rounded_test.cpp`

**Interfaces:**
- Consumes: `fillet_edges`, `tessellate`, `validate_topology`, `curved_solid_volume`.
- Produces: no new API.

- [ ] **Step 1: Write the failing tests**

Add three test cases:

```cpp
TEST_CASE("rounded box is symmetric at all eight corners") {
  const Brep box = geom::make_box(Point3d{-2.0, -3.0, -4.0}, Point3d{0.0, 0.0, 0.0}, 1e-9).value();
  const auto result = geom::fillet_edges(box, all_edges(box), geom::FilletOptions{0.25, 1e-9});
  REQUIRE(result.has_value());
  CHECK(volume_of(result.value().brep) ==
        doctest::Approx(23.53807140332).epsilon(1e-12));
  CHECK(result.value().brep.vertex_count() == 24u);
  CHECK(result.value().brep.edge_count() == 48u);
  CHECK(result.value().brep.face_count() == 26u);
}

TEST_CASE("rounded box tessellates watertight") {
  const Brep box = box_2x3x4();
  const auto result = geom::fillet_edges(box, all_edges(box), geom::FilletOptions{0.25, 1e-9});
  REQUIRE(result.has_value());
  const auto mesh = geom::tessellate(result.value().brep, geom::TessellationOptions{});
  REQUIRE(mesh.has_value());
  CHECK(mesh.value().triangles.size() > 0u);
  std::map<std::pair<std::uint32_t, std::uint32_t>, std::uint32_t> uses;
  for (const auto& triangle : mesh.value().triangles) {
    for (int e = 0; e < 3; ++e) {
      std::uint32_t a = triangle[static_cast<std::size_t>(e)];
      std::uint32_t b = triangle[static_cast<std::size_t>((e + 1) % 3)];
      if (a > b) std::swap(a, b);
      ++uses[{a, b}];
    }
  }
  std::uint32_t boundary = 0;
  for (const auto& [edge, count] : uses) {
    if (count == 1) ++boundary;
  }
  CHECK(boundary == 0u);
}

TEST_CASE("rounded box is deterministic and leaves its input unchanged") {
  const Brep box = box_2x3x4();
  const auto first = geom::fillet_edges(box, all_edges(box), geom::FilletOptions{0.25, 1e-9});
  const auto second = geom::fillet_edges(box, all_edges(box), geom::FilletOptions{0.25, 1e-9});
  REQUIRE(first.has_value());
  REQUIRE(second.has_value());
  CHECK(first.value().brep.vertex_count() == second.value().brep.vertex_count());
  CHECK(volume_of(first.value().brep) == doctest::Approx(volume_of(second.value().brep)));
  for (const geom::VertexId id : first.value().brep.vertex_ids()) {
    CHECK(first.value().brep.vertex(id).point[0] == second.value().brep.vertex(id).point[0]);
    CHECK(first.value().brep.vertex(id).point[1] == second.value().brep.vertex(id).point[1]);
    CHECK(first.value().brep.vertex(id).point[2] == second.value().brep.vertex(id).point[2]);
  }
  CHECK(box.vertex_count() == 8u);
  CHECK(box.edge_count() == 12u);
  CHECK(box.face_count() == 6u);
  CHECK(volume_of(box) == doctest::Approx(24.0));
}
```

- [ ] **Step 2: Run the tests to verify they fail**

Expected: the symmetry test fails if the construction is hard-coded to the positive corner; the watertight/determinism tests may pass once Task 1 is correct — keep them as regression guards. Do not weaken Task 1.

- [ ] **Step 3: Fix the construction**

Generalize the corner frame to all eight sign combinations; do not special-case the low-low-low corner. If `face_uv_project_sphere` returns an error for a loop, carry `previous_u` across the whole loop (as specified in Task 1 Step 3.5) rather than resetting it per point.

- [ ] **Step 4: Run Debug and Release filtered tests**

Run `--test-case="rounded box*"` in both configurations. Expected: 4 cases, all pass.

- [ ] **Step 5: Commit**

```bash
git add tests/unit/fillet_box_rounded_test.cpp
git commit -m "test(brep): prove rounded-box symmetry, watertightness and determinism"
```

---

### Task 3: Diagnostics, radius boundary and selection validation

**Files:**
- Modify: `include/geom/brep/detail/fillet_box_rounded.hpp`
- Modify: `tests/unit/fillet_box_rounded_test.cpp`

**Interfaces:**
- Consumes: `fillet_edges`.
- Produces: stable errors:
  - `unsupported "fillet_edges: radius must be less than half the smallest box dimension"` when `radius >= 0.5*min(a,b,c)`.
  - `index_out_of_range "fillet_edges: edge does not belong to the solid"` for a dead edge in the twelve-edge selection.
  - existing `unsupported "fillet_edges: only an axis-aligned box is supported"` for a non-box.
  - existing option errors unchanged.

- [ ] **Step 1: Write the failing tests**

```cpp
TEST_CASE("rounded box rejects a radius at half the smallest dimension") {
  const Brep box = box_2x3x4();
  const auto edges = all_edges(box);
  const auto at_limit = geom::fillet_edges(box, edges, geom::FilletOptions{1.0, 1e-9});
  REQUIRE_FALSE(at_limit.has_value());
  CHECK(at_limit.error().code == geom::Errc::unsupported);
  CHECK(at_limit.error().message.find("half the smallest box dimension") != std::string_view::npos);
  const auto near_limit = geom::fillet_edges(box, edges, geom::FilletOptions{1.0 - 1e-6, 1e-9});
  CHECK(near_limit.has_value());
}

TEST_CASE("rounded box rejects an incomplete or dead selection") {
  const Brep box = box_2x3x4();
  auto edges = all_edges(box);
  edges.pop_back();
  const auto eleven = geom::fillet_edges(box, edges, geom::FilletOptions{0.25, 1e-9});
  CHECK_FALSE(eleven.has_value());
  edges = all_edges(box);
  edges.push_back(edges.front());
  const auto duplicated = geom::fillet_edges(box, edges, geom::FilletOptions{0.25, 1e-9});
  CHECK(duplicated.has_value());
  CHECK(duplicated.value().filleted_edges == 12u);
  edges = all_edges(box);
  edges.push_back(geom::EdgeId::make(999u, 0u));
  const auto dead = geom::fillet_edges(box, edges, geom::FilletOptions{0.25, 1e-9});
  REQUIRE_FALSE(dead.has_value());
  CHECK(dead.error().code == geom::Errc::index_out_of_range);
}
```

Add a non-box case using a prism/cylinder solid if available, or an L-shaped planar solid built with `SolidEditor`; assert `unsupported` and the box-only message.

- [ ] **Step 2: Run the tests to verify they fail**

Expected: `at_limit` currently assembles a degenerate solid or returns a different message; `near_limit` fails at radius `.999999999` because the nine-face construction cannot shrink; `dead` may be swallowed by `fillet_box_all_edges` returning false and the prism path returning a different message.

- [ ] **Step 3: Implement the guards**

In `fillet_box_rounded` (before building any face):
```cpp
const double smallest = std::min({box.hi[0] - box.lo[0], box.hi[1] - box.lo[1],
                                  box.hi[2] - box.lo[2]});
if (!(radius < 0.5 * smallest)) {
  return Error{Errc::unsupported,
               "fillet_edges: radius must be less than half the smallest box dimension", 0u, 0u};
}
```
In `fillet_box_all_edges`, a dead edge in an otherwise twelve-edge selection must not return `false` (which would fall through to the family path): return the `index_out_of_range` error from `fillet_edges` by checking aliveness before the rounded dispatch, or make the dispatch return `std::optional<Result<...>>`. Keep the public contract: dead edge always `index_out_of_range`, even when the selection is incomplete.

- [ ] **Step 4: Run Debug and Release filtered tests**

Expected: 6 cases, all pass; the near-limit case succeeds and the exact limit is rejected.

- [ ] **Step 5: Commit**

```bash
git add include/geom/brep/detail/fillet_box_rounded.hpp tests/unit/fillet_box_rounded_test.cpp
git commit -m "fix(brep): bound the rounded-box radius and reject dead edges"
```

---

### Task 4: Benchmark and evidence

**Files:**
- Create: `bench/bench_m4c3c.cpp`
- Modify: `CMakeLists.txt`
- Modify: `.github/workflows/ci.yml`
- Create: `docs/adr/0044-rounded-box-spherical-corners.md`
- Create: `docs/benchmarks/m4c3c-rounded-box.md`
- Modify: `docs/architecture/geometry-library-architecture.md` §8.9/§16.4
- Modify: this plan's Results table.

**Interfaces:**
- Consumes: `make_box`, `fillet_edges`, `curved_solid_volume`.
- Produces: `geom_bench_m4c3c` executable; benchmark artifact.

- [ ] **Step 1: Add the benchmark**

Workload: `100 x 100 x 100` box, radius `1.0`, tolerance `1e-9`, twelve edges; build once per sample, five samples. Verify each sample's volume against
`999743.8170303602` with relative tolerance `1e-9`, solid count 1, `filleted_edges == 12`, `validate_topology` and V/E/F `24/48/26`; exit 2 on mismatch. Print min/median/max milliseconds and the measured volume.

- [ ] **Step 2: Build and run**

Configure a Release bench build with `-DGEOM_BUILD_BENCH=ON`; run `geom_bench_m4c3c`; expect exit 0 and a volume within `1e-9` relative of the anchor.

- [ ] **Step 3: Document**

Write ADR-0044: decision (direct 26-face construction instead of boolean or faceted approximation), scope (axis-aligned box, all twelve edges, one radius), spherical corner geometry and volume formula, rejected alternatives, and failure contract. Add the benchmark table to `docs/benchmarks/m4c3c-rounded-box.md`, the M4c-3c paragraph to §8.9/§16.4, and fill the Results table with measured counts, timings, commit range, and the full-Release totals.

- [ ] **Step 4: Wire CI**

Add `geom_bench_m4c3c` to the build targets, a run step, and an artifact upload in the bench job, following `geom_bench_m4c3b`.

- [ ] **Step 5: Run the full verification**

Debug filtered `*rounded*` tests; Release full suite; `ctest` Debug and Release. Expected: all green, 0 warnings.

- [ ] **Step 6: Commit**

```bash
git add bench/bench_m4c3c.cpp CMakeLists.txt .github/workflows/ci.yml docs/adr/0044-rounded-box-spherical-corners.md docs/benchmarks/m4c3c-rounded-box.md docs/architecture/geometry-library-architecture.md docs/superpowers/plans/2026-11-27-m4c3c-rounded-box.md
git commit -m "docs(brep): record the rounded-box contract and benchmark"
```

---

## Acceptance

- `fillet_edges(box, all 12 edges, r)` returns a valid, closed, single-shell solid with 24 vertices, 48 edges, 26 faces (6 planes, 12 cylinders, 8 spheres).
- Volume matches the Steiner formula
  `ABC + 2(AB+BC+CA)r + pi(A+B+C)r^2 + 4pi r^3/3`, `A=a-2r`, `B=b-2r`, `C=c-2r`,
  within `1e-9` relative.
- Tessellation is watertight; the result is deterministic; the input is unchanged.
- `r >= 0.5*min(a,b,c)`, dead edges, non-boxes and partial selections fail with stable errors.
- Existing box-family and prism fillet tests are unchanged and pass.

## Results

| Item | Value |
|---|---|
| Bench (100 mm cube, r = 1 mm) | PENDING |
| Volume | 999743.8170303602 (analytic) |
| V / E / F | 24 / 48 / 26 |
| Debug filtered | PENDING |
| Release full suite | PENDING |
| Commit range | PENDING |
