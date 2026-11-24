# M4b-2: pole and point-gap healing

## Goal

Close the two remaining hole classes that `heal_brep` cannot repair today:

1. **A pole gap** — two faces of a sphere of revolution (a band and a pole cap) separated by a gap so close to the pole that the latitude projection saturates. `heal_brep` currently fails with `face_loops_uv: point is off its circle`.
2. **A point gap** — a cone pinhole whose boundary rim collapses to a point (two vertices within `max_gap`, zero-length boundary arcs). Today this either hits `heal_plan: boundary edge has zero length`, or silently returns a success with an **open** shell.

At the end of this plan both classes heal to a validated, closed, single-shell solid with the analytic volume, and every unrepairable input returns a deterministic error that leaves the input untouched.

## Architecture

The change stays inside the existing healing pipeline, in this order:

```
validation -> closed-rim split (M4b-1) -> point-gap weld (new)
           -> boundary -> plan pairs -> carrier guards -> welds -> pcurves
           -> shell closed -> boundary check -> Euler -> commit
```

- `face_uv_project_sphere` becomes numerically stable at the pole (`atan2` instead of `asin`). This is a projection fix only; no pcurve-rebuild special case is added and the singular-azimuth borrow in `heal_rebuild_face_pcurves` stays.
- `SolidEditor::dissolve_edge` removes a zero-length boundary edge from its loop(s) while keeping its vertices, so the point-weld pass can turn a pinhole rim into a single welded vertex.
- `detail::heal_weld_point_gaps` clusters boundary vertices within `options.max_gap`, dissolves the zero-length boundary edges that join a cluster, merges the cluster, and then lets the ordinary pipeline continue. If it welds anything, the empty-plan early return must not report success while boundary edges remain.
- A malformed analytic surface (invalid frame) is rejected with `heal_brep: surface is invalid` instead of silently producing degenerate pcurves.

## Tech Stack

C++20; CMake 4.x with the existing presets; MSVC 19.44, GCC 15.2, Clang 21.1 in CI; doctest 2.4.11; no new dependencies.

## Spec

- `docs/architecture/geometry-library-architecture.md` §16.4 (healing contract) and §4.6 (tolerance chain).
- `docs/adr/0026-healing-stitch.md` (decisions to extend).
- `docs/quality/testing-strategy.md`: every public API gets positive / boundary / failure coverage; predicates stay exact; tolerance is only a decision input, never an exactness claim.
- M4b-1 contract: closed rims are split on a private working copy; rejections return the caller's input bit-identically.

## Global Constraints

- Never report `changed == true` for a solid that is still open: after the weld pass, the
  existing boundary/Euler/commit checks stay authoritative.
- All geometric decisions use exact predicates or the documented tolerance guards; no new
  global epsilon.
- Every rejection returns the original `Brep` and a full `HealReport`; no partial editor state
  is committed.
- One commit per task; TDD (red first, then green); Release fast lane at task level, full
  Debug + Release at the milestone review.
- Do not change public signatures of `heal_brep`, `HealOptions`, or the existing
  `HealReport` fields; new fields are additive.

## Scope and deviations (need approval)

1. **Pole fix is a projection fix only.** `face_uv_project_sphere` switches the latitude to
   `atan2(axial, hypot(x, y))`. The pcurve rebuild keeps its current singular-azimuth borrow;
   no pole-aware rebuild is introduced. This fixes rims closer than ~1.5e-8 to the pole.
2. **Point-gap welding is single-shell only.** The dissolve/weld pass runs on boundary edges
   of the target solid. Multi-shell point gaps stay unreachable because `validate` already
   rejects a solid with an open shell, and closing two shells at one point needs a different
   ownership model. Deferred.
3. **`dissolve_edge` is boundary-only.** It handles an edge used twice in one loop (a pinhole
   rim) or once in each of two loops. Interior-edge dissolve and T-junction welds are deferred.
4. **A pinhole with a small nonzero radius is closed only when its vertex pair is within
   `max_gap`.** With `r0 = 1e-6` and the default `max_gap = 1e-9` the two vertices are 2e-6
   apart, so the pass must not weld; the call must return a deterministic error instead of the
   current silent open-shell success.
5. **Invalid analytic frames become a rejection.** A surface whose `is_valid()` is false is
   reported as `heal_brep: surface is invalid` (context = surface id) before any planning.
   This only changes inputs that were already malformed; all primitive builders produce valid
   frames.

## Review Focus

1. Does the `atan2` latitude keep every existing sphere test bit-identical away from the pole,
   and is the singular-azimuth borrow still needed/reachable?
2. Is `dissolve_edge` loop surgery safe: `loop.start` relocation, last-coedge refusal, edges
   used once in each of two loops, pcurves, and commit invariants?
3. Does the point-weld pass run before `euler_before`, mark the right faces dirty, and never
   allow an open result to be reported as success?
4. Are the analytic volumes (upper hemisphere `2pi/3`, cone `2pi/3`) asserted with tolerances
   that the old code fails (real RED)?
5. Do all rejections return the input fingerprint unchanged, and are the new messages exact?
6. Does the wider suite stay green (`*sphere*`, `*heal*`, `*trim_face*`, canonical hashes)?

## File Map

- `include/geom/brep/face_uv.hpp` — `face_uv_project_sphere` latitude.
- `include/geom/brep/solid_editor.hpp` — `dissolve_edge`.
- `include/geom/brep/heal.hpp` — point-weld pass, plan/report fields, guard ordering, invalid
  surface check.
- `tests/unit/brep_face_uv_test.cpp` — projection round-trip near the pole.
- `tests/unit/brep_editor_test.cpp` — dissolve surgery.
- `tests/unit/brep_heal_curved_test.cpp` — pole gap, pinhole, point-gap failures.
- `tests/unit/m4b_heal_fixture.hpp` — shared pole-gap and cone-pinhole builders.
- `docs/adr/0026-healing-stitch.md` — addendum.
- `docs/benchmarks/m4b2-healing.md` — measured values.
- `docs/architecture/geometry-library-architecture.md` — §16.4 paragraph.
- `docs/superpowers/plans/2026-11-24-m4b2-pole-and-point-gaps.md` — this plan; results at the
  end of Task 5.

## Task 1: Stabilise the sphere latitude at the pole

**Files**

- `include/geom/brep/face_uv.hpp`
- `tests/unit/brep_face_uv_test.cpp`
- `tests/unit/brep_heal_curved_test.cpp`
- `tests/unit/m4b_heal_fixture.hpp` (create; pole builder only in this task)

**Interfaces**

- Consumes: `Sphere3{center, axis, x_axis, radius}`; `face_uv_project_sphere(const Sphere3&,
  const Point3d&, double previous_u, bool has_previous_u)`; `heal_brep`.
- Produces: unchanged signatures; the latitude of a point at `v = pi/2 - 1e-9` is no longer
  rounded to `pi/2`.

**Steps**

1. In the test fixture, port `make_sphere_pole_gap(delta, gap)` from the M4b-2 probe:
   upper hemisphere split into a band `v in [0, v_up]` and a pole cap `v in [v_cap, pi/2]`,
   with `v_up = pi/2 - delta` and `v_cap = v_up + gap`, both on the sphere chart
   `Sphere3{center, {0,0,1}, {1,0,0}, radius}`. The cap loop is
   `{capa,-} {capb,-} {seam_high,-} {seam_high,rev}`.
2. RED: `TEST_CASE("face_uv: a rim one nanometre below the pole keeps its latitude")` —
   build the point at `v = pi/2 - 1e-9` on the unit sphere, call
   `face_uv_project_sphere`, then `face_uv_eval`; assert the latitude is strictly below
   `pi/2` and the round-trip 3D error is `<= 1e-14`. On the current code the latitude is
   exactly `pi/2` and the error is `1e-9`.
3. RED: `TEST_CASE("heal: a sphere pole gap heals to the hemisphere volume")` — cases
   `(delta, gap) = (1e-8, 5e-9)`, `(1e-4, 5e-5)`, `(5e-8, 2.5e-8)`. Expected volume is the
   hemisphere `2*pi/3 = 2.0943951023931953`, relative `1e-12`; boundary 0, one shell, closed,
   deterministic. On the current code the first case fails with
   `face_loops_uv: point is off its circle` (radial error ~1e-9).
4. Implement the latitude:

   ```cpp
   const double axial = dot(offset, surface.axis);
   const Vec3d y_axis = cross(surface.axis, surface.x_axis);
   const Vec3d radial = offset - surface.axis * axial;
   const double x = dot(radial, surface.x_axis);
   const double y = dot(radial, y_axis);
   const double v = std::atan2(axial, std::hypot(x, y));
   ```

   Keep the azimuth threshold and the `previous_u` fallback unchanged.
5. Run the three cases; then `-tc="*face_uv*"` and `-tc="*sphere*"` to prove no regression.
   Commit `fix(face_uv): measure the sphere latitude away from the pole`.

## Task 2: Add `SolidEditor::dissolve_edge`

**Files**

- `include/geom/brep/solid_editor.hpp`
- `tests/unit/brep_editor_test.cpp`
- `tests/unit/m4b_heal_fixture.hpp` (pinhole builder)

**Interfaces**

- Produces:

  ```cpp
  // Removes `edge` from every loop that uses it, keeping its vertices alive.
  // Refuses an edge whose removal would empty a loop.
  Result<void> dissolve_edge(EdgeId edge);
  ```

  Errors: `index_out_of_range "solid_editor: edge is not alive"`;
  `unsupported "solid_editor: dissolve would empty a loop"` (context = edge id).

**Steps**

1. Port `add_cone_pinhole(editor, z_base, r0, split_rim)` and `add_cone_pinch(editor, z_base)`
   from the probe into the fixture header. The split pinhole has two zero-length top arcs
   `top_a`/`top_b` sharing both top vertices; the closed form has one `v0 == v1` arc.
2. RED: `TEST_CASE("solid editor: dissolving a twice-used pinhole rim keeps the loop")` —
   dissolve `top_a`; the lateral loop still has five coedges, `top_a` is dead, both top
   vertices and `top_b` are alive, commit succeeds and `edge_count()` drops by one.
   On the current code the call does not exist (compile RED).
3. RED: `TEST_CASE("solid editor: dissolving the last coedge of a loop is refused")` —
   `add_cone_pinch`; dissolve the base loop's two coedges; the second call returns
   `unsupported` with the exact message and leaves the editor usable (first dissolve still
   applied, commit still valid).
4. RED: `TEST_CASE("solid editor: dissolving an edge shared by two faces removes both uses")`
   — dissolve `base0` from the pinhole; the base face loop has one coedge, the lateral loop has
   five, and commit succeeds.
5. Implement: collect every live coedge that uses `edge`; group by loop; if any loop's live
   coedge count equals the number of uses being removed, return the empty-loop error before
   mutating. Otherwise, for each use, link `prev.next = next` and `next.prev = prev`; if
   `loop.start == use`, move it to `next`; release the coedge. Then release the edge.
6. Run the three cases and `-tc="*editor*"`; commit `feat(brep): add SolidEditor::dissolve_edge`.

## Task 3: Weld collapsed boundary points before pairing

**Files**

- `include/geom/brep/heal.hpp`
- `tests/unit/brep_heal_curved_test.cpp`
- `tests/unit/m4b_heal_fixture.hpp`

**Interfaces**

- Produces (detail):

  ```cpp
  struct HealPointWeld {
    std::vector<HealVertexWeld> welds;
    std::vector<EdgeId> dissolved;
  };
  Result<HealPointWeld> heal_weld_point_gaps(SolidEditor& editor,
                                             std::span<const EdgeId> boundary,
                                             double max_gap);
  ```

  `HealPlan` gains `std::vector<HealVertexWeld> point_welds;` and
  `std::vector<EdgeId> dissolved;`; `HealReport` gains
  `std::uint32_t collapsed_edges = 0;`.

**Steps**

1. RED: `TEST_CASE("heal: a cone pinhole closes at a welded point")` — split pinhole with
   `r0 = 0`, `max_gap = 1e-9`, `max_tolerance = 1e-9`. Expect
   `volume == 2*pi/3` relative `1e-12`, `boundary 0`, one shell, closed, chi 2,
   `welded_vertices == 2`, `collapsed_edges == 2`, `V/E/F = 3/3/2`, deterministic.
   On the current code this fails at `heal_plan: boundary edge has zero length`.
2. RED: `TEST_CASE("heal: a tiny pinhole rim welds only within max_gap")` — split pinhole
   `r0 = 1e-6`, `max_gap = 1e-5`, `max_tolerance = 1e-4`: volume `2*pi/3` relative `1e-9`,
   same topology and `collapsed_edges == 2`. Then the same fixture with `max_gap = 1e-9`
   must return `invalid_result "heal_brep: no mutual pair"`, `changed == false`, input
   fingerprint unchanged. On the current code the small-gap case silently reports success with
   an open shell.
3. Implement the pass after `heal_split_closed_rims` and before `report.boundary_edges` /
   `euler_before`:
   - union-find over the boundary edge endpoints, joining two vertices whose distance is
     `<= max_gap`; each cluster records its diameter;
   - a boundary edge is a dissolve candidate when both endpoints are in one cluster and the
     edge is used once (the boundary list is already unique), and its vertices are distinct;
   - dissolve all candidates (Task 2), then merge each cluster into its lowest live vertex with
     tolerance `max(existing tolerances, cluster diameter)`; propagate rejections unchanged;
   - record welds/dissolved in the plan and `collapsed_edges` in the report; add the
     `HealIssue` path only for `max_tolerance` violations.
4. Make the empty-plan guard honest: after the weld pass, if boundary edges remain and no pair
   was found, return `invalid_result "heal_brep: no mutual pair"` (do not early-return success
   when `boundary` is non-empty). If the weld pass closed the last boundary edge, skip pairing
   and continue to shell-closed/Euler/commit.
5. Run the two cases plus `-tc="*heal*"`; commit
   `feat(heal): weld collapsed boundary points before pairing`.

## Task 4: Reject malformed analytic surfaces

**Files**

- `include/geom/brep/heal.hpp`
- `tests/unit/brep_heal_curved_test.cpp`

**Interfaces**

- Error: `unsupported "heal_brep: surface is invalid"`, context = surface index.

**Steps**

1. RED: build the probe8-style malformed sphere (`Sphere3{center, radius}`, so
   `axis == x_axis`) on a valid B-rep topology; `heal_brep` must return that error and leave the
   input untouched. On the current code it proceeds and returns the off-circle failure.
2. Add the check in the surface whitelist loop, before closed-rim splitting:
   `Plane3/Cylinder3/Cone3/Sphere3/Torus3` each expose `is_valid()`; if false, return the error
   with `face.surface`.
3. Run the case and `-tc="*heal*"`; commit
   `fix(heal): reject analytic surfaces with an invalid frame`.

## Task 5: Consolidate fixtures, document the contract

**Files**

- `tests/unit/m4b_heal_fixture.hpp` (cleanup)
- `docs/adr/0026-healing-stitch.md`
- `docs/benchmarks/m4b2-healing.md`
- `docs/architecture/geometry-library-architecture.md`
- this plan

**Steps**

1. Move the builders used by Tasks 1–4 into `m4b_heal_fixture.hpp`; delete the local copies;
   keep the existing `m4a_curved_heal_fixture.hpp` untouched.
2. ADR-0026 addendum: pole latitude via `atan2`, `dissolve_edge` semantics, point-weld ordering
   and limits (single-shell only; `max_gap` gate; no interior dissolve), the honest
   no-mutual-pair rule, and the invalid-frame rejection.
3. Benchmark note: measured volumes/errors for the pole gaps (5e-9, 5e-5, 2.5e-8), pinhole
   `r0=0` and `r0=1e-6`, and the rejection cases; state the analytic anchors separately from the
   measured returns.
4. Architecture §16.4: one paragraph naming the two new classes and the deferred multi-shell
   point gap.
5. Full Debug and Release suites; update the count in every document; record results and
   commits in this plan's Results section.
6. Commit `docs(heal): record the pole and point-gap contract`.

## Exit Criteria

- Pole-gap cases `(1e-8, 5e-9)`, `(5e-8, 2.5e-8)`, `(1e-4, 5e-5)` heal to `2*pi/3`
  (relative `1e-12`), boundary 0, one closed shell, chi 2, deterministic.
- Split cone pinhole `r0 = 0` heals to `2*pi/3` (relative `1e-12`) with `V/E/F = 3/3/2`,
  `welded_vertices == 2`, `collapsed_edges == 2`.
- Split cone pinhole `r0 = 1e-6` with `max_gap = 1e-5` heals to `2*pi/3` (relative `1e-9`);
  with `max_gap = 1e-9` it returns `invalid_result "heal_brep: no mutual pair"`, `changed == false`,
  and the input fingerprint unchanged.
- Every rejection in the new tests returns the caller's `Brep` bit-identically.
- Invalid analytic frame returns `unsupported "heal_brep: surface is invalid"`.
- Existing sphere/trim/heal suites stay green; canonical mesh hashes unchanged.
- Full Debug and Release suites pass; counts and docs agree.

## Self-Review

- [ ] Each task starts from a failing test that the previous commit fails.
- [ ] No test asserts only counts; each has geometry (volume, boundary, shell, chi).
- [ ] `dissolve_edge` refuses the empty-loop case before mutating.
- [ ] The weld pass cannot produce a success while boundary edges remain.
- [ ] Rejections leave the input untouched.
- [ ] Pole fix does not change away-from-pole latitudes beyond round-off.
- [ ] Multi-shell point gap is explicitly deferred.
- [ ] Fast-lane commands are given; full suites run once at the end.
- [ ] Docs separate analytic anchors from measured returns.

## Results

2026-11-24, branch `sdd/m4b2`, base `3bad7e1`; docs land as `docs(heal): record the pole and
point-gap contract` (`9496510`), followed by the post-review fixes `fbc1ea9` and `c41c6e2`.

**Commits**

- `f7f277f` fix(face_uv): measure the sphere latitude away from the pole
- `68b0a43` test(face_uv): cover the remaining sphere pole-gap cases
- `dd1d76a` feat(solid_editor): dissolve an edge from its loops
- `8af9474` feat(heal): weld collapsed boundary points before pairing
- `feb257b` fix(heal): reject analytic surfaces with an invalid frame
- `faf0b70` test(heal): compare heal issue messages by value
- `9496510` docs(heal): record the pole and point-gap contract
- `fbc1ea9` test(brep): format failure messages without error() on success
- `c41c6e2` fix(brep): restrict dissolve_edge to a single loop

**Fixtures.** All builders used by Tasks 1–4 live in `tests/unit/m4b_heal_fixture.hpp`
(`make_sphere_pole_gap`, `make_malformed_sphere_pole_gap`, `add_cone_pinhole`,
`add_cone_pinch`); a grep over `tests/unit/*.cpp` and `*.hpp` found no local duplicates, so no
file had to move. `tests/unit/m4a_curved_heal_fixture.hpp` is untouched.

**Measured (Release, clean rebuild).** Pole `(1e-8, 5e-9)` and `(5e-8, 2.5e-8)` with default
options, `(1e-4, 5e-5)` with `{1e-9, 1e-4, 1e-4}`: boundary 4, welded 2 / sewn 2, chi 2, one
closed shell, volume `2*pi/3`; `max_closed_gap` = `5e-9` / `2.5e-8` /
`4.9999999994791671e-05`; volume deviations `1.19904e-14` / `4.88498e-15` / `1.95399e-14`.
Pinholes r0=0 (`{1e-9, 1e-9, 1e-9}`) and r0=1e-6 (`{1e-9, 1e-5, 1e-4}`): boundary 0, welded 2,
sewn 0, collapsed 2, V/E/F `3/3/2`, chi 2, volume deviations `1.33227e-15` / `4.44089e-15`.
r0=1e-6 with `{1e-9, 1e-9, 1e-4}`: `changed == false`, boundary 2, one `invalid_result
"heal_brep: no mutual pair"`, input fingerprint unchanged. r0=1e-6 with `{1e-9, 1e-5, 1e-9}`:
`changed == false`, boundary 2, welded/collapsed 0, one `invalid_result
"heal_brep: tolerance cap exceeded"`, input fingerprint unchanged (the cluster needs 2e-6 and the
caller capped the tolerance at 1e-9). Invalid frame: `unsupported
"heal_brep: surface is invalid"` with `context_a == 1` (surface store index), fingerprint
unchanged.

**Lanes (Release, clean rebuild).** `*heal*` 95/2,013; `*torus*` 59/10,754; merged
149/12,587 (936 skipped); `*uv*` 9/388; `*dissolve*` 7/85; `*pinhole*` 3/89; `*pole*` 17/1,389;
0 failed.

**Full suites.** Release (`--clean-first` rebuild after the header edits): 1086 cases /
680,749 assertions, 0 failed (1074/680,503 before this milestone, i.e. +12 cases/+246
assertions).
Debug full suite **PENDING — deferred to the milestone review**: sharded Debug runs covered
all but the 8 slow cases with 553,598 assertions and 0 failures over about 5 hours; the 8
cases are all `brep_boolean_curved_test.cpp` sphere/torus half-space booleans that take
>20 min each under `/Od`, so the lead stopped the run rather than block the milestone. Debug
lanes that did run: `*heal*` 95/2,013 and `*uv*` 9/388, both 0 failed. All Release lanes and
the Debug heal lane used the freshly rebuilt binaries; no incremental objects were trusted.

**Deviations.** The full Debug suite was not completed (see above; PENDING at the milestone).
The fixture-move step was a no-op because the builders
were already consolidated in M4b-2 Tasks 1–4. Final review found that the public
`dissolve_edge` accepted an edge shared by two loops and could commit a solid whose loops no
longer close (`validate_topology` does not check loop closure); `c41c6e2` now rejects that case
with `unsupported`, returns `invalid_topology` for an emptied loop, and the corresponding test
asserts rejection instead of the old relinking. The same commit adds the direct point-weld
`max_tolerance` rejection test (`{1e-9, 1e-5, 1e-9}`), so the cap is no longer covered only
positively. `fbc1ea9` fixes a Debug-only test bug (pointer comparison of `const char*` issue
messages) and routes boolean test messages through a guarded helper, so `--success` runs no
longer call `error()` on a success value.

The second review round found that a same-loop, single-use edge still committed a loop whose
consecutive coedges no longer met at a vertex (`face_loops_uv` reported "coedge chain is not
continuous", while `validate_topology` stayed green). `commit()` now re-checks vertex continuity
for exactly the loops a `dissolve_edge` touched, and rejects the candidate with
`invalid_topology "loop coedge chain is not vertex-continuous"`; the check is deliberately not
global because many older fixtures stage intentionally partial topology for validator tests and
the point-weld pass stages one rim arc before the other. The pinhole fixture is committed first,
then re-opened through `SolidEditor(const Brep&)`, so rollback after the refused commit returns
to the five-edge solid.
