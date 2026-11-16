# M3c-3c-1 General Sphere Lenses Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use subagent-driven-development (recommended) or executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Make `boolean_solid_curved` support two whole spheres whose centres lie on any line, not only on a coordinate axis.

**Architecture:** The M3c-3b re-chart already aligns each sphere's `axis` with the centre line, which turns the intersection circle into a latitude of the chart. Two defects keep that path from working off-axis. First, an SSI sample projects onto the seam meridian a few ulp inside the parameter domain; the seam-pair closure accepts that drift, but the arrangement does not, so the cut stays a dangling chain and a zero-width face appears. Second, once the cut is attached, each fragment's rim is lifted through its own sphere chart and healed independently, so the two caps can carry different vertex or arc representations of the same physical rim and the assembler reports "edge is used once". The fix snaps near-seam ends onto the domain side and snaps cut-rim points onto the one analytic pair circle before healing, so both fragments share the same rim geometry.

**Tech Stack:** C++20, CMake 4.4.4 + MSVC 19.44, doctest 2.4.11, `/W4 /WX /utf-8`.

**Spec:** `docs/architecture/geometry-library-architecture.md` §8.3/§8.10, `docs/quality/testing-strategy.md`, ADR-0029, ADR-0030, and the M3c-3b probe report `.superpowers/sdd/m3c3b-probe-report.md`.

## Global Constraints

- Repository `D:\ai_geometry`; this plan branches from `f6e2814` (master after M3c-3b).
- No new dependency. Public API of `boolean_solid_curved` unchanged.
- Predicates decide topology; floating point only pre-filters. No global epsilon.
- Every error return sets `failure_report` through the existing `record_failure()` path.
- Tests are registered in `tests/CMakeLists.txt`; the four touched test files are already registered.
- Debug and Release both pass, `GEOM_WERROR=ON`, zero C++ warnings.
- MSBuild header dependency tracking is unreliable in this repository: after editing a header, `touch` one `.cpp` of the target or use `/t:Rebuild` before trusting a run.

## Scope and deviations

In scope: whole-sphere/sphere operands, crossing, any centre-line direction; union, intersection, difference; exact analytic pair circle; volumes through `sphere_face_volume` (the library `curved_solid_volume` stalls on 1000+ point rims and is not used here).

Out of scope, unchanged and documented: sphere faces that are not whole spheres, containment, tangency, separation, cone, torus, NURBS, boolean-result tessellation. Deferred to M3c-3c-2/3/4.

Approved deviation: the rim unification is done by projection onto the existing pair circle, not by extending the seam-only splitter to general closed cuts. The latter stays deferred.

## Review Focus

1. A cut endpoint 1 ulp inside the seam: the splitter must attach it to the side, not create a zero-width face.
2. A rim lifted through two different charts: both fragments must expose the same 3D rim vertices and the same arc spans, or the assembler must fail cleanly.
3. An off-axis placement whose centre line is within rounding of an axis: accepted and built, not rejected.
4. A non-crossing or non-whole-sphere operand: documented `unsupported`, input unchanged.
5. A rotated but otherwise identical pair: volume equals the axis-aligned value within 1e-6 relative.

---

### Task 1: Snap near-seam cut endpoints onto the domain side

**Files:**
- Modify: `include/geom/brep/trim_face.hpp:967-1027` (cut normalization in `trim_face_split_seam_only`)
- Test: `tests/unit/brep_trim_face_test.cpp`

**Interfaces:**
- Consumes: `split_face_curved(const Brep&, FaceId, std::span<const FaceCut>, double)`, `FaceCut{std::vector<Point2d> uv}`.
- Produces: no signature change; `trim_face_split_seam_only` snaps any cut coordinate within `1e-12 * max(1, |u_lo|, |u_hi|)` of `u_lo`/`u_hi` (and of `v_lo`/`v_hi` when `projection.periodic_v()`) onto that side before the arrangement.

- [ ] **Step 1: Write the failing test**

```cpp
TEST_CASE("curved boolean: a seam endpoint a few ulp inside is snapped") {
  const double tolerance = 1e-9;
  const geom::Result<geom::Brep> sphere =
      geom::make_sphere(geom::Point3d{0.0, 0.0, 0.0}, 1.0, tolerance);
  REQUIRE(sphere.has_value());
  const geom::FaceId face = sphere.value().face_ids().front();

  std::vector<geom::Point2d> ring;
  constexpr std::size_t kSamples = 25;
  const double v = std::asin(0.5);
  for (std::size_t i = 0; i <= kSamples; ++i) {
    const double u = 2.0 * 3.14159265358979323846 *
                     static_cast<double>(i) / static_cast<double>(kSamples);
    ring.push_back(geom::Point2d{u, v});
  }
  // The SSI projection of a seam point is not exactly on the seam meridian.
  ring.front()[0] = 2.9101132668452026e-18;
  const std::vector<geom::FaceCut> cuts{geom::FaceCut{ring}};
  const geom::Result<std::vector<geom::CurvedFragment>> fragments =
      geom::split_face_curved(sphere.value(), face, cuts, tolerance);
  REQUIRE(fragments.has_value());
  REQUIRE(fragments.value().size() == 2);
  for (const geom::CurvedFragment& fragment : fragments.value()) {
    REQUIRE(fragment.loops_uv.size() == 1);
    // The zero-width ring produced by a dangling chain has 2 * kSamples + 2
    // points; the two real pieces stay at kSamples + 1 (rim) plus one pole.
    CHECK(fragment.loops_uv.front().size() <= kSamples + 3);
    CHECK(geom::signed_area2(fragment.loops_uv.front()) != 0.0);
  }
}
```

- [ ] **Step 2: Run the test to verify it fails**

Run: `.\build\tests\Debug\geom_tests.exe --test-case="*a seam endpoint a few ulp inside*"`
Expected: FAIL — one fragment, its `loops_uv.front().size()` is 52, `signed_area2` is 0.

- [ ] **Step 3: Snap the endpoints**

In `trim_face_split_seam_only`, after the whole-chain u shift and before the bounds check, clamp each point:

```cpp
const double seam_eps =
    1e-12 * std::max({1.0, std::fabs(u_lo), std::fabs(u_hi)});
for (Point2d& point : points) {
  if (std::fabs(point[0] - u_lo) <= seam_eps) {
    point[0] = u_lo;
  } else if (std::fabs(point[0] - u_hi) <= seam_eps) {
    point[0] = u_hi;
  }
  if (periodic_v) {
    if (std::fabs(point[1] - v_lo) <= seam_eps) {
      point[1] = v_lo;
    } else if (std::fabs(point[1] - v_hi) <= seam_eps) {
      point[1] = v_hi;
    }
  }
}
```

Recompute `min_u`/`max_u` from the snapped points before the existing bounds check, because the original values were captured before snapping.

- [ ] **Step 4: Run the test to verify it passes**

Run: `.\build\tests\Debug\geom_tests.exe --test-case="*a seam endpoint a few ulp inside*"`
Expected: PASS, two fragments.

- [ ] **Step 5: Commit**

```bash
git add include/geom/brep/trim_face.hpp tests/unit/brep_trim_face_test.cpp
git commit -m "fix(m3c3c): snap near-seam cut ends onto the parameter side"
```

---

### Task 2: Remove the axis-aligned placement guard

**Files:**
- Modify: `include/geom/brep/boolean_solid_curved.hpp:1808-1826` (`curved_boolean_rechart_whole_spheres`)
- Test: `tests/unit/brep_boolean_curved_test.cpp`

**Interfaces:**
- Consumes: `curved_boolean_rechart_whole_spheres(const Brep&, SolidId, const Brep&, SolidId, double)`, `make_sphere(Point3d, Vec3d, Vec3d, double, double)`.
- Produces: the helper accepts any crossing placement and re-charts with `axis = unit(cb - ca)` (no dominant-axis snap). The fixed-point comparison still requires both axes and both `x_axis`s to agree within `8*epsilon`, so recursion depth stays two.

- [ ] **Step 1: Write the failing test**

```cpp
TEST_CASE("curved boolean: an off-axis sphere lens is built") {
  const double tolerance = 1e-9;
  const double angle = 15.0 * 3.14159265358979323846 / 180.0;
  const geom::Result<geom::Brep> south =
      geom::make_sphere(geom::Point3d{0.0, 0.0, 0.0}, 1.0, tolerance);
  const geom::Result<geom::Brep> north = geom::make_sphere(
      geom::Point3d{1.5 * std::cos(angle), 1.5 * std::sin(angle), 0.0}, 1.0,
      tolerance);
  REQUIRE(south.has_value());
  REQUIRE(north.has_value());

  const geom::Result<geom::CurvedBooleanResult> lens = geom::boolean_solid_curved(
      south.value(), south.value().solid_ids().front(), north.value(),
      north.value().solid_ids().front(), geom::SolidBooleanOp::intersection,
      geom::CurvedBooleanOptions{tolerance});
  REQUIRE_MESSAGE(lens.has_value(), message_of(lens.error()));
  const std::optional<double> volume = sphere_face_volume(lens.value().brep);
  REQUIRE(volume.has_value());
  CHECK(*volume == doctest::Approx(11.0 * 3.14159265358979323846 / 96.0)
                      .epsilon(1e-6));
}
```

- [ ] **Step 2: Run the test to verify it fails**

Run: `.\build\tests\Debug\geom_tests.exe --test-case="*an off-axis sphere lens is built*"`
Expected: FAIL with `boolean_solid_curved: sphere/sphere placement is not axis-aligned; general placements are deferred to M3c-3c`.

- [ ] **Step 3: Accept any crossing direction**

Delete the dominant-axis snap and the `kAxisAlignment` guard. Keep the crossing checks and the shared-frame choice:

```cpp
const Vec3d axis{delta[0] / distance, delta[1] / distance, delta[2] / distance};
const Vec3d negated_axis{-axis[0], -axis[1], -axis[2]};
```

The `candidate` selection already uses `|axis[0]| >= 0.9`, so it remains valid for any direction.

- [ ] **Step 4: Run the test to verify it fails for the next reason**

Run: `.\build\tests\Debug\geom_tests.exe --test-case="*an off-axis sphere lens is built*"`
Expected: FAIL with `boolean_solid_curved: result touches itself along an edge`. That is the rim-representation defect Task 3 owns; this failure is the acceptance signal for Step 3.

- [ ] **Step 5: Commit**

```bash
git add include/geom/brep/boolean_solid_curved.hpp tests/unit/brep_boolean_curved_test.cpp
git commit -m "feat(m3c3c): re-chart sphere pairs on any crossing direction"
```

---

### Task 3: Detect seam-only sphere charts by relative area

**Files:**
- Modify: `include/geom/brep/trim_face.hpp:1511-1524` (seam-only dispatch in `split_face_curved`)
- Modify: `tests/unit/brep_boolean_curved_test.cpp` (`sphere_face_cap_volume` pole threshold)

**Interfaces:**
- Consumes: `signed_area2(loops.front())`, `extent` (already computed at `trim_face.hpp:1493`), `detail::trim_face_seam_only_domain(const SurfaceGeometry&)`.
- Produces: a seam-only face is recognized when `|signed_area2(loop)| <= 1e-12 * extent * extent`, not only when it is exactly zero. The seam-only branch still requires `trim_face_seam_only_domain` to return a domain, so a genuinely trimmed face is unaffected.

- [ ] **Step 1: Keep the acceptance test and extend it**

The `curved boolean: an off-axis sphere lens is built` case (Task 2) already fails at its intersection assertion. Extend it with the union assertion from Task 3's original brief:

```cpp
  const geom::Result<geom::CurvedBooleanResult> united = geom::boolean_solid_curved(
      south.value(), south.value().solid_ids().front(), north.value(),
      north.value().solid_ids().front(), geom::SolidBooleanOp::union_,
      geom::CurvedBooleanOptions{tolerance});
  REQUIRE_MESSAGE(united.has_value(), message_of(united.error()));
  CHECK(*sphere_face_volume(united.value().brep) ==
        doctest::Approx(245.0 * 3.14159265358979323846 / 96.0).epsilon(1e-6));
```

- [ ] **Step 2: Run the test to verify it fails for the diagnosed reason**

Run: `.uild	ests\Debug\geom_tests.exe --test-case="*an off-axis sphere lens is built*"`
Expected: FAIL with `boolean_solid_curved: result touches itself along an edge`. The splitter produced no fragment for operand A: its 64-point seam loop has uv area `6.1629758220391547e-33`, not exactly zero, so the dispatch took the regular trim path and its area filters dropped every fragment.

- [ ] **Step 3: Use a relative area bound and widen the test helper's pole bound**

At `trim_face.hpp:1515`, replace the exact comparison:

```cpp
  // The seam loop is a slit: its uv area is a rounding artefact, not a region
  // area. A relative bound proportional to extent^2 recognizes it on charts
  // whose frame is not axis aligned; a genuinely trimmed boundary stays far
  // above it.
  if (std::fabs(signed_area2(loops.front())) <= 1e-12 * extent * extent) {
```

Keep the `seam_domain` gate and both returns unchanged.

In `tests/unit/brep_boolean_curved_test.cpp`, the helper `sphere_face_cap_volume` accepts a pole when `|v|` is within `1e-9` of `pi/2`. Off-axis charts recover the pole through `asin`, where the `|axis|^2 = 1 +/- 1e-16` rounding is amplified to about `2e-8`. Change that threshold to `1e-6` and add the reason as a comment; the cap height error introduced is `O(1e-12)`, far below the 1e-6 volume check.

Remove any rim-snapping helper added during the blocked attempt (`curved_boolean_snap_rim_to_circle`): the 15-degree probe shows both caps' rims are already identical, so the helper is inert and out of scope.

- [ ] **Step 4: Run the test to verify it passes**

Run: `.uild	ests\Debug\geom_tests.exe --test-case="*an off-axis sphere lens is built*"`
Expected: PASS; intersection `11π/96` and union `245π/96`.

Then run the full Debug suite and the full Release suite; both must be green.

- [ ] **Step 5: Commit**

```bash
git add include/geom/brep/trim_face.hpp tests/unit/brep_boolean_curved_test.cpp
git commit -m "fix(m3c3c): detect seam-only sphere charts by relative area"
```

### Task 4: General-placement acceptance corpus

**Files:**
- Test: `tests/unit/brep_boolean_curved_test.cpp`
- Modify: `docs/superpowers/plans/2026-11-13-m3c2-curved-booleans.md` only if it repeats the old axis-only claim (search for "not axis-aligned").

**Interfaces:**
- Consumes: the Task 2/3 driver, `sphere_face_volume`, `validate_topology`.
- Produces: no production API; the corpus pins placement, volume, topology, determinism and rejection behavior.

- [ ] **Step 1: Write the failing corpus test**

```cpp
TEST_CASE("curved boolean: sphere lenses accept rotated placements") {
  const double tolerance = 1e-9;
  const double pi = 3.14159265358979323846;
  const geom::Result<geom::Brep> south =
      geom::make_sphere(geom::Point3d{0.0, 0.0, 0.0}, 1.0, tolerance);
  REQUIRE(south.has_value());

  const auto make_pair = [&](double degrees, double distance) {
    const double angle = degrees * pi / 180.0;
    return geom::make_sphere(
        geom::Point3d{distance * std::cos(angle), distance * std::sin(angle),
                      0.0},
        1.0, tolerance);
  };

  for (int degrees = 0; degrees < 360; degrees += 15) {
    INFO("degrees " << degrees);
    const geom::Result<geom::Brep> north = make_pair(degrees, 1.5);
    REQUIRE(north.has_value());
    for (const geom::SolidBooleanOp op :
         {geom::SolidBooleanOp::intersection, geom::SolidBooleanOp::union_,
          geom::SolidBooleanOp::difference}) {
      const geom::Result<geom::CurvedBooleanResult> result =
          geom::boolean_solid_curved(
              south.value(), south.value().solid_ids().front(), north.value(),
              north.value().solid_ids().front(), op,
              geom::CurvedBooleanOptions{tolerance});
      REQUIRE_MESSAGE(result.has_value(),
                      "degrees " << degrees << ": " << message_of(result.error()));
      CHECK(geom::validate_topology(result.value().brep).ok());
      const std::optional<double> volume = sphere_face_volume(result.value().brep);
      REQUIRE(volume.has_value());
      const double expected =
          op == geom::SolidBooleanOp::intersection
              ? 11.0 * pi / 96.0
              : (op == geom::SolidBooleanOp::union_ ? 245.0 * pi / 96.0
                                                    : 39.0 * pi / 32.0);
      CHECK(*volume == doctest::Approx(expected).epsilon(1e-6));
    }
  }
}
```

- [ ] **Step 2: Run the corpus to verify it fails on the first non-axis angle**

Run: `.\build\tests\Debug\geom_tests.exe --test-case="*accept rotated placements*"`
Expected before Task 3: FAIL at degrees 15. After Task 3: PASS for every 15-degree step.

- [ ] **Step 3: Extend with a one-degree sweep and a fixed random set**

Add a second loop over `degrees = 0..359` that only runs `intersection` with `distance = 1.5`, and 200 deterministic `mt19937` draws (seed `20261116`) of `distance ∈ [1.05, 1.95]` and `degrees ∈ [0, 360)`, comparing to `11π/96` with `1e-3` relative (chord residual grows for large distances). Keep the runtime under 4 minutes on Debug; drop to 60 random draws if it does not.

- [ ] **Step 4: Run the extended corpus**

Run: `.\build\tests\Debug\geom_tests.exe --test-case="*accept rotated placements*"`
Expected: PASS; report the count in the task report.

- [ ] **Step 5: Commit**

```bash
git add tests/unit/brep_boolean_curved_test.cpp
git commit -m "test(m3c3c): cover rotated sphere lens placements"
```

---

### Task 5: QA-17, ADR-0031 and the status sections

**Files:**
- Create: `tests/unit/curved_boolean_qa17_test.cpp`
- Create: `docs/adr/0031-general-sphere-lenses.md`
- Create: `docs/benchmarks/m3c3c.md`
- Modify: `tests/CMakeLists.txt` (register the QA file)
- Modify: `docs/architecture/geometry-library-architecture.md` §16.3
- Modify: `docs/superpowers/plans/2026-11-13-m3c2-curved-booleans.md` if it still claims axis-only support

**Interfaces:**
- Consumes: everything above.
- Produces: QA-17 corpus and the documentation of the now-supported placements and remaining limits.

- [ ] **Step 1: Write the QA-17 test**

```cpp
TEST_CASE("QA-17: sphere lens at an arbitrary placement") {
  const double tolerance = 1e-9;
  const double pi = 3.14159265358979323846;
  const double angle = 37.0 * pi / 180.0;
  const geom::Result<geom::Brep> a =
      geom::make_sphere(geom::Point3d{0.0, 0.0, 0.0}, 1.0, tolerance);
  const geom::Result<geom::Brep> b = geom::make_sphere(
      geom::Point3d{1.5 * std::cos(angle), 1.5 * std::sin(angle), 0.0}, 1.0,
      tolerance);
  REQUIRE(a.has_value());
  REQUIRE(b.has_value());
  const auto lens = geom::boolean_solid_curved(
      a.value(), a.value().solid_ids().front(), b.value(),
      b.value().solid_ids().front(), geom::SolidBooleanOp::intersection,
      geom::CurvedBooleanOptions{tolerance});
  REQUIRE(lens.has_value());
  CHECK(geom::validate_topology(lens.value().brep).ok());
  REQUIRE(lens.value().brep.solid_shells(lens.value().brep.solid_ids().front()).size() == 1);
  CHECK(*sphere_face_volume(lens.value().brep) ==
        doctest::Approx(11.0 * pi / 96.0).epsilon(1e-6));
  // The result must be reproducible bit for bit.
  const auto again = geom::boolean_solid_curved(
      a.value(), a.value().solid_ids().front(), b.value(),
      b.value().solid_ids().front(), geom::SolidBooleanOp::intersection,
      geom::CurvedBooleanOptions{tolerance});
  REQUIRE(again.has_value());
  CHECK(*sphere_face_volume(again.value().brep) ==
        *sphere_face_volume(lens.value().brep));
}
```

The QA file needs a local copy of `sphere_face_volume` or a small shared test header; prefer moving the helper into `tests/unit/brep_boolean_curved_fixture.hpp` and including it from both files.

- [ ] **Step 2: Run it**

Run: `.\build\tests\Debug\geom_tests.exe --test-case="QA-17*"`
Expected: PASS.

- [ ] **Step 3: Write ADR-0031 and the benchmark page**

ADR-0031 records: accepted inputs (whole spheres, crossing, any centre line), the seam snapping rule (`1e-12` relative), the shared-rim rule (project onto the analytic pair circle within `max(64*tol, 1e-9*r)`), and remaining limits (non-crossing, non-whole-sphere, cone/torus/NURBS, tessellation of the result). `docs/benchmarks/m3c3c.md` records the 37-degree QA run and the 15-degree step cost.

- [ ] **Step 4: Update §16.3 and run the full suites**

Add an M3c-3c-1 paragraph with the QA-17 number and the new acceptance. Run:

```powershell
cmake --build build --config Debug --target geom_tests
.\build\tests\Debug\geom_tests.exe
cmake --build build-rel --config Release --target geom_tests
.\build-rel\tests\Release\geom_tests.exe
```

Expected: all pass, counts recorded in the task report, zero warnings.

- [ ] **Step 5: Commit**

```bash
git add tests/unit/curved_boolean_qa17_test.cpp tests/CMakeLists.txt docs/adr/0031-general-sphere-lenses.md docs/benchmarks/m3c3c.md docs/architecture/geometry-library-architecture.md
git commit -m "docs(m3c3c): record QA-17 and the general sphere lens limits"
```

---

## Exit criteria

- Tasks 1–5 committed; 15-degree, 1-degree and random rotated pairs all build valid single-shell lenses whose sphere-cap volume matches the analytic value.
- A seam endpoint 2.9e-18 inside the domain yields two fragments, never a zero-width ring.
- The axis-aligned guard is gone; non-crossing and non-whole-sphere inputs still return their documented `unsupported` messages with `face_pairs == 0`.
- Debug and Release suites pass with zero warnings; ADR-0031/QA-17/§16.3 updated.

## Self-review

- Spec coverage: every Review Focus line is pinned by a task test (Task 1 seam drift, Task 2 acceptance, Task 3 rim identity, Task 4 corpus, Task 5 QA).
- Type consistency: `curved_boolean_snap_rim_to_circle` is the only new production name; it is used only in `harvest`.
- No placeholders: all tests carry exact values; no step says "handle edge cases".
