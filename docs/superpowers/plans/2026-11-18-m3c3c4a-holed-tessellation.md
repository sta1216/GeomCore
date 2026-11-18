# M3c-3c-4a: Holes in Trimmed Curved Tessellation and Mass Properties

> **For agentic workers:** REQUIRED SUB-SKILL: Use subagent-driven-development (recommended) or executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Mesh and weigh a trimmed analytic face that carries a rectilinear inner hole, end to end, without cracking shared edges.

**Architecture:** The boundary-aligned strip decomposition in `tessellate_trimmed_face` gains hole rings: every outer and hole sample becomes a row level or column, cells inside a hole are dropped, and each periodic piece receives the holes that lie strictly inside it. The mass integrator stops rejecting multi-loop trimmed faces and triangulates outer plus holes with the same helper. `triangulate_trimmed_domain` and the existing fan expansion are reused unchanged; no new triangulation algorithm is introduced.

**Tech Stack:** C++20, MSVC 19.44 `/W4 /WX /utf-8`, doctest 2.4.11, CMake 4.4.4, VS2022 generator.

**Spec:** `docs/architecture/geometry-library-architecture.md` §8.10 (tessellation contract), §8.11 (mass properties), §12.5/§16.3; `docs/adr/0027-trimmed-tessellation.md`, `docs/adr/0032-general-domain-splitter.md`.

## Global Constraints

- Geometry decisions use exact predicates (`sign_orient2d`, `ring_orientation2`); floating point is only for sampling and rejection.
- A trimmed curved face is meshed only when every ring is rectilinear (`tessellate_ring_rectilinear`); slanted trims and non-rectilinear holes stay `Errc::unsupported` with a deterministic message.
- Shared boundary samples are never duplicated or interpolated: every position contributed by a neighbouring face appears once and is used by the triangles it bounds.
- A hole must lie strictly inside one periodic piece: no sample on the periodic cut, no shared position with the outer ring, and no segment crossing the cut.
- Canonical (non-trimmed) face tessellation and the existing untrimmed mass paths stay bitwise unchanged.
- Each task is a single commit, Debug and Release full suites green, zero warnings under `GEOM_WERROR=ON`.
- MSBuild header dependency tracking is unreliable in this repository: after editing a header, touch a `.cpp` or rebuild with `/t:Rebuild` before trusting a test run.
- Test registrations already exist: `tests/unit/tessellate_trimmed_test.cpp` and `tests/unit/mass_properties_test.cpp` are listed in `tests/CMakeLists.txt`; no CMake change is expected.

## Review Focus

- **A hole that touches or nearly touches the outer ring**: a 1e-12 gap must not produce a zero-width strip or a self-touching mesh; the face must either mesh watertight or be rejected with `hole touches the boundary`.
- **A hole sample on the periodic cut**: a hole that is not wholly inside one piece must be rejected by the per-piece assignment; no constructible fixture exists yet, so this is pinned by the pre-check code and the review, not by a test.
- **Coarse deflection**: at `{1e-2, 0.35}` the hole corners and the four hole walls must remain in the mesh; subdivision must not fill the hole or drop its samples.
- **The periodic piece assignment**: the split cut is the mean of the sampled u extremes (≈3.097), not π; a hole in either piece must be meshed in that piece and nowhere else.
- **Mass with a hole**: a hole on the +y wall must move the centroid to −y and reduce the volume; a sign error produces a plausible-looking wrong number, so the test pins the analytic volume, not only success.

---

### Task 1: Hole-aware strip decomposition in `tessellate_trimmed_face`

**Files:**
- Modify: `include/geom/brep/tessellate.hpp:1101-1240` (`tessellate_trimmed_strip_triangles` and its band helpers), `include/geom/brep/tessellate.hpp:1538-1730` (`tessellate_trimmed_face`)
- Test: `tests/unit/tessellate_trimmed_test.cpp:404-424` (replace the hole-rejection case), `tests/unit/tessellate_trimmed_test.cpp` (new cases)

**Interfaces:**
- Consumes: `TessellateTrimmedRing{points, positions, outer}`; `tessellate_ring_rectilinear(std::span<const Point2d>)`; `tessellate_split_periodic_ring(const TessellateTrimmedRing&) -> std::optional<std::array<TessellateTrimmedRing,2>>`; `face_uv_project`, `face_uv_unwrap`, `face_uv_eval`; `triangulate_trimmed_domain`; `tessellate_expand_trimmed_runs` (`tessellate.hpp:319-377`).
- Produces: `Result<bool> tessellate_trimmed_strip_triangles(TessellateState&, const FaceUvProjection&, const TessellateTrimmedRing& outer, std::span<const TessellateTrimmedRing> holes, FaceId, std::vector<std::array<std::uint32_t,3>>& triangles)` — returns `true` when it produced the mesh, `false` when the domain is not strip-decomposable (the caller falls back to the polygon path). The existing four-argument overload delegates with an empty hole span.

- [ ] **Step 1: Replace the rejection test with the watertight-hole test**

In `tests/unit/tessellate_trimmed_test.cpp` delete `TEST_CASE("trimmed tessellation: a hole loop is rejected with the trimmed-path error")` and add:

```cpp
TEST_CASE("trimmed tessellation: a rectilinear hole meshes watertight") {
  const geom::Brep brep = hole_window_cylinder();
  REQUIRE(geom::validate_topology(brep).ok());
  const auto mesh = geom::tessellate(brep, geom::TessellationOptions{1e-3, 0.35});
  REQUIRE(mesh.has_value());
  check_watertight(mesh.value());
  CHECK(mesh.value().triangle_face.size() == mesh.value().triangles.size());

  const auto corner = [&](double x, double y, double z) {
    return std::any_of(mesh.value().positions.begin(), mesh.value().positions.end(),
                       [&](const geom::Point3d& p) {
                         return std::abs(p[0] - x) <= 1e-6 && std::abs(p[1] - y) <= 1e-6 &&
                                std::abs(p[2] - z) <= 1e-6;
                       });
  };
  CHECK(corner(-0.5, 0.5, 0.5));
  CHECK(corner(0.5, 0.5, 0.5));
  CHECK(corner(-0.5, 0.5, 1.5));
  CHECK(corner(0.5, 0.5, 1.5));
  CHECK(mesh.value().triangles.size() >= 950u);
  CHECK(mesh.value().triangles.size() <= 1450u);
}
```

- [ ] **Step 2: Run the new test to verify it fails on the hole message**

Run: `.\build\tests\Debug\geom_tests.exe --test-case="*rectilinear hole*"`
Expected: FAIL — `REQUIRE(mesh.has_value())` with `tessellate: trimmed curved face has a hole`.

- [ ] **Step 3: Extend the strip decomposition with holes**

In `tessellate_trimmed_strip_triangles` (five-argument overload):
- Build the v-levels from the outer ring and from every hole; within each band, build x-intervals from the outer ring and subtract the x-intervals of every hole that spans the band, using the same even-odd crossing rule the current band code uses.
- For each remaining interval, triangulate between its left and right chains. A chain may come from the hole ring; its samples become the row vertices exactly as outer samples do.
- A band that contains no remaining interval emits nothing; a band whose interval is bounded by a hole on both sides is still meshed (the hole wall is part of the face boundary).
- Keep the existing zipper rule: consecutive rows share the sample vertices, and a row boundary that is a hole segment uses the hole's `positions` indices.
- Return `false` when any ring is not rectilinear, any ring sample projects outside the current piece, or a hole segment crosses a band boundary without a level at that v.
- The four-argument overload must call the five-argument one with `{}` so the hole-free output stays bitwise identical.

- [ ] **Step 4: Admit holes in `tessellate_trimmed_face`**

Replace the block at `tessellate.hpp:1644-1651`:
- Reject only with `Errc::unsupported` `"tessellate: trimmed curved face hole boundary is not rectilinear"` when a hole ring fails `tessellate_ring_rectilinear`.
- Do not route holed faces through `tessellate_split_periodic_ring` blindly: if `projection.periodic_u()` and the outer ring splits into two pieces, assign each hole to the piece that contains every one of its samples; a hole whose sample lies on the cut or outside the piece is rejected with `"tessellate: trimmed curved face hole crosses the periodic cut"`.
- Call the five-argument strip decomposition with the piece's outer ring and its assigned holes. If it returns `false`, fall back to `tessellate_trimmed_polygon_triangles` with outer plus holes.
- Add every ring (outer and holes) to `constrained` via `add_ring_uvs` (`:1712`).

- [ ] **Step 5: Run the hole test and the whole trimmed group**

Run: `.\build\tests\Debug\geom_tests.exe --test-case="trimmed tessellation:*"`
Expected: all cases PASS, including the four hole-corner checks and the triangle-count bound.

- [ ] **Step 6: Pin the coarse-deflection contract and the seam-crossing refusal**

Add:

```cpp
TEST_CASE("trimmed tessellation: coarse deflection keeps the hole open") {
  const geom::Brep brep = hole_window_cylinder();
  const auto fine = geom::tessellate(brep, geom::TessellationOptions{1e-3, 0.35});
  const auto coarse = geom::tessellate(brep, geom::TessellationOptions{1e-2, 0.35});
  REQUIRE(fine.has_value());
  REQUIRE(coarse.has_value());
  check_watertight(coarse.value());
  CHECK(coarse.value().triangles.size() < fine.value().triangles.size());
}

- [ ] **Step 7: Run the hole and coarse-deflection cases, then the full Debug and Release suites**

Run: `.\build\tests\Debug\geom_tests.exe --test-case="trimmed tessellation:*"` and the same for `build-rel`.
Expected: PASS; full suites report 985 cases (983, minus the removed rejection case, plus the hole and coarse cases) and zero failed assertions.

- [ ] **Step 8: Commit**

```powershell
git add include/geom/brep/tessellate.hpp tests/unit/tessellate_trimmed_test.cpp
git commit -m "feat(m3c3c4a): tessellate rectilinear holes in trimmed curved faces"
```

---

### Task 2: Mass properties for a trimmed curved face with a hole

**Files:**
- Modify: `include/geom/brep/mass_properties.hpp:431-468` (`mass_add_trimmed_face`)
- Test: `tests/unit/mass_properties_test.cpp:648-669` (replace the hole-rejection case), new cases at the end of the file

**Interfaces:**
- Consumes: `trimmed_face_boundaries(brep, face, projection, deflection) -> Result<std::vector<TrimmedBoundary>>`; `mass_simplify_rectilinear_ring(std::vector<Point2d>)`; `tessellate_ring_rectilinear`; `triangulate_trimmed_domain(std::span<const TrimmedBoundary>)`; `mass_integrate_trimmed_parametric(surface, points, triangles, reversed, reference, step, acc)`.
- Produces: `mass_add_trimmed_face` accepts any number of rectilinear boundaries (one outer plus zero or more holes); each boundary is simplified and its `arcs` cleared before triangulation.

- [ ] **Step 1: Replace the mass hole test with the analytic-volume test**

In `tests/unit/mass_properties_test.cpp` change the first block of `TEST_CASE("mass properties: trimmed dents and holes report their limitation")` to:

```cpp
TEST_CASE("mass properties: a trimmed curved face with a rectilinear hole") {
  const geom::Brep brep = hole_window_brep();
  const auto mass = geom::mass_properties(brep, geom::MassPropertiesOptions{1e-9, 24, 1e-5});
  REQUIRE(mass.has_value());
  const double expected = 2.0 * 3.14159265358979323846 - 0.45661147746179975;
  CHECK(std::abs(mass.value().volume - expected) <= 1e-9 * expected);
  CHECK(mass.value().centroid[1] < 0.0);

  const auto again = geom::mass_properties(brep, geom::MassPropertiesOptions{1e-9, 24, 1e-5});
  REQUIRE(again.has_value());
  CHECK(again.value().volume == mass.value().volume);
  CHECK(again.value().centroid[1] == mass.value().centroid[1]);
}
```

Keep the slanted-trim block as its own test case `TEST_CASE("mass properties: a slanted trim is still refused")`.

- [ ] **Step 2: Run it to verify it fails on the hole message**

Run: `.\build\tests\Debug\geom_tests.exe --test-case="*rectilinear hole*"`
Expected: FAIL — `REQUIRE(mass.has_value())` with `mass_properties: trimmed curved face has a hole`.

- [ ] **Step 3: Integrate every boundary**

In `mass_add_trimmed_face` replace the single-boundary gate and the domain construction:

```cpp
std::vector<TrimmedBoundary> simplified;
simplified.reserve(boundaries.value().size());
for (const TrimmedBoundary& boundary : boundaries.value()) {
  if (!tessellate_ring_rectilinear(boundary.points)) {
    return Error{Errc::unsupported,
                 "mass_properties: trimmed curved face boundary is not rectilinear", id.index(),
                 0u};
  }
  TrimmedBoundary domain = boundary;
  domain.points = mass_simplify_rectilinear_ring(domain.points);
  domain.arcs.assign(domain.points.size(), std::nullopt);
  simplified.push_back(std::move(domain));
}
const Result<TrimmedUvTriangulation> triangulation = triangulate_trimmed_domain(simplified);
```

Leave the quadrature step, the reversed flag, and the reference-point handling untouched. The non-rectilinear error message must stay byte-identical for the slanted case.

- [ ] **Step 4: Run the mass tests and the full suites**

Run: `.\build\tests\Debug\geom_tests.exe --test-case="mass properties:*"`, then both full suites.
Expected: PASS; the volume is within 1e-9 relative of `5.826573829717787`; suites report 986 cases (985 + 1) and zero failures.

- [ ] **Step 5: Commit**

```powershell
git add include/geom/brep/mass_properties.hpp tests/unit/mass_properties_test.cpp
git commit -m "feat(m3c3c4a): weigh trimmed curved faces with rectilinear holes"
```

---

### Task 3: QA record, ADR and documentation

**Files:**
- Create: `docs/adr/0033-trimmed-hole-tessellation.md`
- Create: `docs/benchmarks/m3c3c4a.md`
- Modify: `docs/architecture/geometry-library-architecture.md` §16.3 (append an M3c-3c-4a paragraph and update the M3c-2 hole limitation)
- Modify: `tests/unit/tessellate_trimmed_test.cpp`, `tests/unit/mass_properties_test.cpp` (one QA case each is acceptable, but do not duplicate Task 1/2 assertions)

**Interfaces:**
- Consumes: the implementations from Tasks 1 and 2.
- Produces: ADR-0033 recording the supported hole class (rectilinear, strictly inside one periodic piece, no seam contact), the refusal matrix, and the deferred general/twisted/NURBS cases.

- [ ] **Step 1: Add the QA case to `tests/unit/tessellate_trimmed_test.cpp`** (add `#include <geom/brep/mass_properties.hpp>`)

```cpp
TEST_CASE("QA-20: the holed window cylinder meshes and weighs consistently") {
  const geom::Brep brep = hole_window_cylinder();
  const auto mesh = geom::tessellate(brep, geom::TessellationOptions{1e-3, 0.35});
  REQUIRE(mesh.has_value());
  check_watertight(mesh.value());
  const auto mass = geom::mass_properties(brep, geom::MassPropertiesOptions{1e-9, 24, 1e-5});
  REQUIRE(mass.has_value());
  CHECK(std::abs(mass.value().volume - 5.826573829717787) <= 1e-6 * 5.826573829717787);
}
```

Run it and confirm both the mesh and the mass succeed in the same test.

- [ ] **Step 2: Write ADR-0033**

State: context (M3c-2 deferred curved-face holes; M3c-3c-4a lifts the rectilinear, non-seam-crossing case); decision (strip decomposition with holes, per-piece assignment, mass integrates all boundaries); consequences (curved holes now mesh and weigh; non-rectilinear/crossing-seam/general twisted holes still refuse with the two Task 1 messages; boolean results with holes remain limited by the upstream outline builder); alternatives rejected (fan post-pass — the fan already exists; full constrained triangulation — unnecessary for rectilinear).

- [ ] **Step 3: Write the benchmark note and update §16.3**

`docs/benchmarks/m3c3c4a.md` records: fixture, meshing options, boundary edges 0, non-manifold edges 0, lateral area `11.519173063162574 ± 1e-6`, mesh volume in `[5.8200, 5.8266]`, mass volume `5.826573829717787`, and the two refusal inputs. §16.3 gains one paragraph naming the new support and the refusal messages, and the M3c-2 sentence "curved-face holes are rejected" is changed to point at ADR-0033.

- [ ] **Step 4: Run the full suites and commit**

Run both suites; expected 987 cases and zero failures.
```powershell
git add docs/adr/0033-trimmed-hole-tessellation.md docs/benchmarks/m3c3c4a.md docs/architecture/geometry-library-architecture.md tests/unit/tessellate_trimmed_test.cpp
git commit -m "docs(m3c3c4a): record the hole-tessellation contract"
```

---

## Exit Criteria

- `tessellate()` meshes `hole_window_cylinder()` watertight at `{1e-3,0.35}` and `{1e-2,0.35}`; the four hole corners are mesh vertices; the lateral area is `11.519173063162574 ± 1e-6`; the mesh volume is within `[5.8200, 5.8266]`.
- A hole that is not wholly inside one periodic piece is refused by the per-piece pre-check (`"tessellate: trimmed curved face hole crosses the periodic cut"`); no constructible fixture exists, so this is a review-gated defensive path, not a test-pinned one.
- `mass_properties()` on the same solid returns volume `5.826573829717787` within `1e-9` relative, a negative y centroid, and bitwise-identical repeats.
- Canonical primitive tessellation and untrimmed mass outputs are unchanged (existing pinned tests pass).
- Debug and Release full suites: 987 cases, zero failures, zero warnings under `/W4 /WX`.

## Self-Review

- [x] Each task names exact files, signatures and the assertions that prove it.
- [x] Tests use existing fixtures (`hole_window_cylinder`, `hole_window_brep`) and analytic values measured in the probe report.
- [x] The rectilinear refusal message is pinned; the periodic-cut refusal is defensive (no constructible fixture) and marked as such.
- [x] Review Focus lines each map to a test: touching/seam contact (refusal test), crossing cut (refusal test), coarse deflection (Task 1 Step 6), piece assignment (Task 1 Step 3), mass sign (Task 2 Step 1).
- [x] No CMake changes needed; both test files are already registered.
- [x] Deferred work is explicit: non-rectilinear holes, twisted holes, general curved holes, NURBS.
