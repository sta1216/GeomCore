# M3c-3c-4b3a: Trimmed Cone Tessellation Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use subagent-driven-development (recommended) or executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Make every trimmed cone frustum produced by `boolean_solid_curved` and `split_solid_by_plane` tessellate watertight by measuring chord sag on the actual cone surface instead of in parameter space.

**Architecture:** On a cone, the uv midpoint of a chord is not the chord midpoint in 3D because azimuth is angle while the chord is radius-weighted. The trimmed-face refiner and the final deviation check therefore measure a parameter-space error that halves only linearly and never converges. Both measurements become a nearest-surface distance for `Cone3` faces only; splitting stays at the parameter midpoint and every non-cone path is bitwise unchanged.

**Tech Stack:** C++20 header-only GeomCore, CMake 4.4.4 + VS2022 (MSVC 19.44), doctest 2.4.11, `/W4 /WX /utf-8`.

**Spec:** `docs/architecture/geometry-library-architecture.md` §8.10, `docs/adr/0035-cone-boolean.md`, probe report `.superpowers/sdd/m3c3c4b3-probe-report.md`.

## Global Constraints

- Predicates stay exact; no global epsilon. New tolerances are relative and documented at the use site.
- Only `Cone3` faces may take the geometric sag path; every other surface keeps the existing parameter metric bitwise. QA-13 cut-cylinder counts stay `triangles == 894`, `positions == 449` at `{deflection 2e-3, angle 0.35}`.
- Canonical primitive counts stay box 12/8, cylinder 88/46, torus 1152/576 (`tests/unit/tessellate_curved_test.cpp`).
- No signature change to `tessellate`, `TessellationOptions`, or `FaceUvProjection`.
- No apex fan, no seam rewrite, no uv-ring apex insertion in this slice. Apex-bearing cone faces keep using `tessellate_cone_face`; an apex ring reaching the trimmed path stays rejected by the rectilinear gate.
- MSBuild header dependency tracking is broken in this repository: after editing `tessellate.hpp`, touch a `.cpp` or run `/t:Rebuild` before trusting a test run.
- Debug and Release both build with `-DGEOM_WERROR=ON`; zero warnings.
- Every task ends with one commit; run the focused tests, then the full suite before committing.

## Review Focus

1. A cone frustum whose two rims have different sample counts (72 vs 51) must mesh watertight with zero boundary/nonmanifold edges.
2. An apex-bearing cone (intersection, above-cut small cone, whole cone) must keep its canonical single-apex fan and must not enter the trimmed path.
3. A non-cone trimmed face (cut cylinder, QA-13) must keep identical counts, positions and volume.
4. A translated cone that tessellates today (split/boolean at `(1,2,3)`) must remain unchanged; the coincident-seam-copy refusal must not move.
5. Repeated tessellation of the same `Brep` must be bitwise identical.

---

### Task 1: Cone-scoped geometric sag

**Files:**
- Modify: `include/geom/brep/tessellate.hpp` (refine predicate near `:1575`, final check near `:2003`)
- Test: `tests/unit/brep_boolean_curved_test.cpp` (after the QA-22 block, around `:554`)

**Interfaces:**
- Consumes: `FaceUvProjection{surface, periodic_u, periodic_v, bounds}`, `face_uv_project(projection, point, hint_u, along_tangent)`, `face_uv_eval(projection, uv)`, `Cone3`.
- Produces: `namespace detail` helper `std::optional<double> tessellate_nearest_sag_square(const FaceUvProjection& projection, const Point3d& point, double hint_u);` in `tessellate.hpp` — returns `nullopt` for non-cone faces or when the projection fails, otherwise the squared 3D distance from `point` to its nearest point on the cone.

- [ ] **Step 1: Write the failing tests**

Add a `TEST_CASE("tessellate: trimmed cone frusta are watertight")` to `tests/unit/brep_boolean_curved_test.cpp`. Use the existing `make_cone`/`make_box`/`boolean_solid_curved` helpers from the QA-22 case (cone `{0,0,0}`, `{0,0,1}`, radius 1, height 2, tol 1e-9; box `{-2,-2,1}..{2,2,4}`; low box `{...}..{2,2,3}`). For each result call:

```cpp
const geom::TessellationOptions options{1e-3, 0.35};
const geom::Result<geom::TriMesh> mesh = geom::tessellate(result.value().brep, options);
REQUIRE(mesh.has_value());
CHECK(mesh.value().positions.size() == 134u);   // difference; union 142, split 134
CHECK(mesh.value().triangles.size() == 264u);   // difference; union 280, split 264
// boundary_edges / nonmanifold_edges helpers: count undirected edges with 1 / >2 incident triangles.
CHECK(boundary_edge_count(mesh.value()) == 0u);
CHECK(nonmanifold_edge_count(mesh.value()) == 0u);
CHECK(mesh_volume(mesh.value()) == doctest::Approx(1.830287422577).epsilon(5e-3));
CHECK(geom::tessellate(result.value().brep, options).value().positions ==
      mesh.value().positions);  // bitwise repeat
```

Cases and expected values (measured at probe report §5):
- QA-22 difference: 134/264, volume `1.830287422577` (`7π/12`).
- QA-22 union: 142/280, volume `33.830287422577`.
- `split_solid_by_plane` at z=1, lower part: 134/264, volume `1.830287422577`.
- QA-22 intersection: 51/98, volume `0.261110903259`, exactly one apex vertex (the small cone keeps the analytic apex at z=2).
- split z=1 upper part (canonical apex): 51/98, volume `0.261110903259`, exactly one apex vertex.

Also assert `triangle_face.size() == triangles.size()` and every id `< brep.face_count()`.

- [ ] **Step 2: Run the tests to verify they fail**

Run: `.\build\tests\Debug\geom_tests.exe --test-case="tessellate: trimmed cone frusta are watertight"`
Expected: FAIL for difference/union/lower split with `Errc::numeric_failure`, message `tessellate: trimmed face refinement budget is exhausted`, context `(1, 0)`; intersection and upper split already pass.

- [ ] **Step 3: Implement the cone sag metric**

In `include/geom/brep/tessellate.hpp`, add the helper above `tessellate_refine_trimmed_edges`:

```cpp
inline std::optional<double> tessellate_nearest_sag_square(
    const FaceUvProjection& projection, const Point3d& point, double hint_u) {
  if (projection.surface == nullptr) return std::nullopt;
  if (!std::holds_alternative<Cone3>(*projection.surface)) return std::nullopt;
  const Result<Point2d> uv = face_uv_project(projection, point, hint_u, false);
  if (!uv.has_value()) return std::nullopt;
  const Point3d nearest = face_uv_eval(projection, uv.value());
  const double dx = point[0] - nearest[0];
  const double dy = point[1] - nearest[1];
  const double dz = point[2] - nearest[2];
  return dx * dx + dy * dy + dz * dz;
}
```

In the refiner (after `surface_middle` is computed and before the sag test), replace the sag with:

```cpp
if (const std::optional<double> geometric =
        tessellate_nearest_sag_square(projection, chord_middle, middle_uv[0])) {
  sag_square = *geometric;
}
```

where `chord_middle` is the existing `0.5 * (pa + pb)` used by the old sag. Keep the split point `surface_middle`.

In the final check loop, for the triangle centroid:

```cpp
double sag_square = /* existing value from centroid_uv */;
if (const std::optional<double> geometric =
        tessellate_nearest_sag_square(projection, centroid, centroid_uv[0])) {
  sag_square = *geometric;
}
```

Keep the error type/message unchanged. Do not touch the rectilinear gate, the strip builder, periodic split, or trimmed boundaries.

- [ ] **Step 4: Run the focused tests to verify they pass**

Run: `.\build\tests\Debug\geom_tests.exe --test-case="tessellate: trimmed cone frusta are watertight"`
Expected: PASS, 5 cases' assertions green.

- [ ] **Step 5: Run the full suites**

Run: `.\build\tests\Debug\geom_tests.exe` and `.\build-rel\tests\Release\geom_tests.exe`
Expected: all pass, zero warnings, 1010 cases (Task 1 adds one case; record the exact assertion totals in the report).

- [ ] **Step 6: Commit**

```bash
git add include/geom/brep/tessellate.hpp tests/unit/brep_boolean_curved_test.cpp
git commit -m "fix(m3c3c4b3a): measure cone chord sag on the surface"
```

---

### Task 2: Extended cone tessellation matrix

**Files:**
- Modify: `include/geom/brep/tessellate.hpp` (rename the Task 1 helper and document its conservative metric)
- Test: `tests/unit/brep_boolean_curved_test.cpp` (new `TEST_CASE`s after Task 1)

**Interfaces:**
- Consumes: Task 1 helper behavior; existing split/boolean builders and `mesh_volume`/edge-count helpers from Task 1.
- Produces: `std::optional<double> tessellate_cone_radial_sag_square(...)` — the renamed helper; its value is the squared distance from the point to the cone point at the same azimuth and axial height, an upper bound on the true nearest-surface distance.

- [ ] **Step 1: Write the matrix and helper tests**

Add `TEST_CASE("tessellate: split cone heights stay watertight")` and `TEST_CASE("tessellate: canonical cone faces keep their apex")` with:
- split z=0.5 lower: 140/276, volume `1.209417839998`.
- split z=1.5 lower: 124/244, volume `2.058752793525`.
- translated split z=4 (cone at `(1,2,3)`): lower 134/264, volume `1.830287422577`; upper 51/98, volume `0.261110903259`.
- translated boolean difference (same `(1,2,3)` cone): 134/264, volume `1.830287422577`.
- apex intersection z=0.5: 62/120, volume `0.882011366069`, exactly one apex vertex within 1e-12 of `(0,0,2)` (the analytic apex; the z=0.5 cap is an ear-clipped polygon).
- apex intersection z=1.5: 37/70, volume `0.032559033313`, one apex vertex.
- whole cone: 72/140, volume `2.091662475470`, one apex vertex.
Every case also asserts boundary 0, nonmanifold 0, `triangle_face` size/id bounds, and bitwise repeat. Compare volumes with `epsilon(1e-9)` against the 12-decimal literals below (stricter than the 5e-3 chordal bound).

Add `TEST_CASE("tessellate: the cone sag helper is cone-only")` that builds a plane projection and a cylinder projection (reuse projections from existing trimmed tests) and asserts `detail::tessellate_cone_radial_sag_square(projection, point, 0.0)` is `nullopt` for both.

- [ ] **Step 2: Rename and document the helper**

Rename `tessellate_nearest_sag_square` to `tessellate_cone_radial_sag_square` at its definition and both call sites. Add a comment: the projection keeps the point's azimuth and axial height, so the value is an upper bound on the true nearest-surface distance — conservative for both the split decision and the final check.

- [ ] **Step 3: Prove the RED is real**

Temporarily restore the pre-Task-1 header: `git checkout dcb5979 -- include/geom/brep/tessellate.hpp`, rebuild, run the new cases.
Expected: the split/translated/union cases fail with the refinement-budget error; the apex/canonical cases pass. Restore with `git checkout HEAD -- include/geom/brep/tessellate.hpp` and rebuild.

- [ ] **Step 4: Run the full suites**

Debug and Release, zero warnings; record exact counts.

- [ ] **Step 5: Commit**

```bash
git add tests/unit/brep_boolean_curved_test.cpp
git commit -m "test(m3c3c4b3a): pin the cone tessellation matrix"
```

---

### Task 3: Contract documentation

**Files:**
- Create: `docs/adr/0036-cone-tessellation-sag.md`
- Create: `docs/benchmarks/m3c3c4b3a.md`
- Modify: `docs/architecture/geometry-library-architecture.md` (§16.3)
- Modify: `docs/adr/0035-cone-boolean.md` (limitation: trimmed cone tessellation now supported for frusta; apex fan still canonical only)

- [x] **Step 1: Write the ADR**

Record: the parameter-vs-geometric sag discrepancy on cones (8.3x at the z=1 frustum, 1032 splits vs 13), the cone-only scoping that protects QA-13/canonical counts, why the strip zipper is sufficient for frusta, and why apex rings stay rejected until the fan slice. Include the measured table from Task 1/2.

- [x] **Step 2: Write the benchmark note**

`docs/benchmarks/m3c3c4b3a.md`: cases, V/T, mesh volumes, max absolute volume error `2.917e-3` (max relative `5.069e-3`), deflection/angle, determinism, platform.

- [x] **Step 3: Update architecture §16.3 and ADR-0035**

Add the 4b3a paragraph and narrow the ADR-0035 limitation sentence: frusta now tessellate; apex faces remain canonical-only.

- [x] **Step 4: Verify docs and commit**

Run `git diff --check`; confirm no production files changed.

```bash
git add docs/
git commit -m "docs(m3c3c4b3a): record the cone tessellation metric"
```

---

## Exit Criteria

- QA-22 difference/union and split z=1 both tessellate: boundary 0, nonmanifold 0, volume within `5e-3` of analytic, bitwise repeat.
- Cone frusta with 72/51 rim samples succeed; apex-bearing faces still produce exactly one apex vertex via the canonical path.
- Non-cone counts unchanged: QA-13 cut cylinder 894/449; canonical box/cylinder/torus counts unchanged; all non-cone face hashes/positions bitwise identical.
- Translated split/boolean that pass today are unchanged; the coincident-seam-copy refusal still fires at `y = -2`/`|y| >= 2.25`.
- Debug and Release full suites pass with zero warnings; ADR-0036, benchmark and §16.3 landed.

## Self-Review

- Spec coverage: §8.10 watertight/deviation contract → Task 1/2 invariants; ADR-0035 cone scope → Task 3.
- Each step names one action and its expected result; test values come from the probe report.
- Names are consistent: `tessellate_cone_radial_sag_square`, `mesh_volume`, `boundary_edge_count`, `nonmanifold_edge_count`.
- Review Focus: each line has a test in Task 1 or 2.
- Final counts: Debug and Release each 1,013 cases / 673,296 assertions, 0 failed; zero warnings with `GEOM_WERROR=ON`.
- Proportion: production change ~30 lines; tests and docs dominate, as expected for a metric fix.
