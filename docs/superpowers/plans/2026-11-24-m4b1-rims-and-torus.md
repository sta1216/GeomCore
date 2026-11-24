# M4b-1: Closed-rim healing and torus gaps

## Goal
Extend `heal_brep` so open shells whose gap boundary is a **closed rim** (a full-turn circle with a single vertex) can be split, paired, welded and sewn, and prove the machinery on torus gaps. A cylinder split into two open tubes with a 1e-6 seam must heal to `pi*r^2*h`; torus bands/halves with 1e-6 gaps must heal to their analytic volumes.

## Architecture
Two production changes in `include/geom/brep/heal.hpp`, plus tests:

1. **Pre-split closed rims.** Before planning, copy the input Brep to a working copy and split every *boundary* edge with `v0 == v1` at the midpoint of its parameter span with `SolidEditor::split_edge`. Failures discard the copy and return the original input. Planning/mutation/rejection continue on the working copy; a rejection still returns the untouched input.
2. **Keep-edge/weld consistency.** `heal_plan` currently chooses the kept edge by lower edge index but normalizes each vertex weld by lower vertex index. After a rim split those disagree: the kept edge's carrier (e.g. circle at z=0) ends at the other rim's welded vertex (z=gap) and the kept-curve residual check fails. Build each pair's welds from the pair's keep/drop edge endpoints instead of `normalize_weld`, so the kept carrier always passes through the kept vertex.

No changes to `SolidEditor`, `face_uv`, pcurve reconstruction or assembly.

## Tech stack
C++20, CMake 4.4.4 + MSVC 19.44, doctest 2.4.11, `/W4 /WX /utf-8`, `-DGEOM_WERROR=ON`.

## Global constraints
- Rejections leave the caller's Brep bit-identical; only a successful result may carry the pre-split topology.
- Splitting is limited to boundary edges of the target solid whose `v0 == v1`; non-closed boundaries take the existing path unchanged.
- Existing M4a curved-gap tests, planar QA-12, and all negative classes must keep their values/messages.
- Planar and curved behavior outside closed rims must be bit-identical; curved pcurves are rebuilt per dirty face as in M4a.
- Fast lane: task-level Release targeted filters; full Debug+Release only at milestone review.
- MSBuild header tracking is unreliable: touch a TU or `/t:Rebuild` after header edits.

## Spec
Closed-rim fixture: a band shell with two circular boundary rims, `v0 == v1`, `t0 == 0`, `t1 == 2*pi`, radius `r`, at heights `z0` and `z1`, plus a cylindrical face per rim-connected strip. Two such bands separated by `gap` along z form the acceptance case.

After healing:
- boundary edges 0, `V = 6`, `E = 8`, `F = 4`, `chi = 2` for the two-band cylinder;
- `welded_vertices = 2`, `sewn_edges = 2`, `max_closed_gap ~= gap`;
- volume `pi*r^2*(h_total)` within 1e-12 relative; idempotent second heal; deterministic repeat.

Torus cases use open bands (already supported shapes) and closed-rim halves; see Task 4/5 values.

## Scope and deviations (need approval)
1. **Closed-rim splitting is automatic inside `heal_brep`**, not a caller pre-step; the split exists only in a working copy.
2. **Rim split parameter is the span midpoint** (for a full turn, opposite the seam vertex), not an arbitrary caller parameter.
3. **Apex/pole gaps** (boundary shrinking to a point; point-to-point weld) are **deferred to M4b-2**: the probe shows they need a distinct pole-pcurve fix and a zero-length edge weld path.
4. **Fillet/chamfer and shell** remain later M4 milestones.
5. Torus acceptance is limited to gaps between whole bands/halves; no partial-u trimming, no NURBS.

## Review Focus
1. A failed heal must return the original input, never the working copy with splits.
2. Splitting only closed *boundary* edges; no interior edges, no non-target solids; shared-solid check still precedes mutation.
3. `max_gap`/`max_tolerance`/tolerance raising unchanged; split vertices start with tolerance 0 and are raised through the existing chain.
4. Kept-carrier selection: the kept edge's curve must pass through the welded keep vertex for line and circle carriers; flipped mappings still reverse the inherited carrier.
5. Euler guard compares post-split working topology, not pre-split input.
6. Torus pcurves across the v seam: the face chart may be `[0,2pi]` while projection returns `(-pi,pi]`; the rebuilt pcurve must not jump.
7. Negative classes unchanged: radius mismatch, plane lid, off-surface, kept-curve mismatch, NURBS face/edge, ambiguous seams.

## File Map
- `include/geom/brep/heal.hpp` — pre-split pass; weld construction; rejection path returns input.
- `tests/unit/brep_heal_curved_fixture.hpp` — closed-rim cylinder bands, torus band/half fixtures.
- `tests/unit/brep_heal_curved_test.cpp` — acceptance and torus tests.
- `tests/unit/brep_heal_test.cpp` — planar regressions only if messages/counts are touched (expected none).
- `docs/adr/0026-healing-stitch.md`, `docs/benchmarks/m4b-healing.md`, `docs/architecture/geometry-library-architecture.md` §16.4 — Task 6.

## Tasks

### Task 1: split closed boundary rims in a working copy
**Files:** `include/geom/brep/heal.hpp`.
**Interfaces:**
- Consumes: `Brep`, `SolidId`, `HealOptions`, `detail::heal_boundary_edges`, `SolidEditor::split_edge`.
- Produces: internal `detail::heal_split_closed_rims(const Brep&, SolidId) -> Result<std::optional<Brep>>` (or an equivalent inline block).
**Steps:**
1. Add a failing test in `brep_heal_curved_test.cpp`: two-band cylinder with closed rims and gap 1e-6; today `unsupported "closed gap edge cannot be paired..."`.
2. Run it and capture the failure (`-tc="*closed rim*"`).
3. Copy input; for every boundary edge with `v0 == v1`, call `split_edge(edge, edge.t0 + 0.5*(edge.t1-edge.t0))`; commit. Recompute boundary/edge maps on the working copy; replace the closed-rim rejection with the split pass; rejection returns the original `input`.
4. Re-run the test (still fails on kept-curve, expected until Task 2) and the existing `*heal*` filters.
5. Commit `feat(heal): split closed gap rims before planning`.

### Task 2: consistent kept edge and welded vertices
**Files:** `include/geom/brep/heal.hpp`.
**Interfaces:** Consumes `HealEdgePair`/`HealVertexWeld`; produces corrected `pair.welds` built as `{keep_edge.vX, drop_edge.vY}`.
**Steps:**
1. Extend the Task 1 test to assert success (RED at the kept-curve message).
2. Run and confirm `"heal_brep: kept edge curve does not match its welded vertices"`.
3. Replace the two `normalize_weld(...)` calls in `heal_plan` with endpoint-directed welds; keep dedup/order semantics.
4. Re-run the closed-rim test and `-tc="*heal*"`.
5. Commit `fix(heal): keep the welded endpoint on the kept edge`.

### Task 3: cylinder closed-rim acceptance
**Files:** `tests/unit/m4a_curved_heal_fixture.hpp`, `tests/unit/brep_heal_curved_test.cpp`.
**Steps:**
1. Add fixture builders: `make_cylinder_closed_rim_band(radius, z0, z1)` (one shell, one cylinder face, two closed rim edges, 2 coedges, pcurves present) and `make_two_band_closed_rim_gap(radius, height, gap)` (two bands, z=[0,height], [height+gap, 2*height+gap]).
2. RED: `CHECK(heal_brep(...))` fails on the diagnostic; after T1/T2 GREEN.
3. Assert volume `pi*(2+gap)` (r=1), boundary 0, post-split working `V=6/E=8/F=4` and healed `V=4/E=6/F=4`, chi==2, welded 2/sewn 2, `max_closed_gap` ~1e-6, validate clean, second heal unchanged, repeat deterministic.
4. Run `-tc="*closed rim*"` and `-tc="*heal*"`.
5. Commit `test(heal): accept a cylinder split by a closed-rim seam`.

### Task 4: torus open-rim gaps
**Files:** fixture + test.
**Steps:**
1. Build torus band fixtures: (a) u-gap total span 3.0 with 1e-6 seam, (b) outer-branch v-gap band, (c) inner-branch band crossing the v seam.
2. RED on the current tree where they fail; record the first failure.
3. After T1/T2 assert volumes `2.3561944901923608`, `5.4334195273316128`, `0.8635058157199863` within 1e-9 relative, boundary 0, single shell, chi==2, deterministic.
4. Run `-tc="*torus*"` and `-tc="*heal*"`.
5. Commit `test(heal): accept torus band gaps`.

### Task 5: torus closed-rim halves and negatives
**Files:** fixture + test.
**Steps:**
1. Build two torus halves (cut at z=+-gap/2) with closed rims; heal to full torus volume `4.934802200544679` (R=1,r=0.5), boundary 0, chi==0 (not 2: V=4,E=6,F=2), one shell.
2. Add negatives: radius mismatch (gap_too_large), plane lid carrier mismatch, off-surface rim, NURBS face/edge, ambiguous seam; assert messages and `context_a/context_b`.
3. Run the torus filters; commit `test(heal): accept closed-rim torus halves with negatives`.

### Task 6: docs and results
**Files:** ADR-0026 addendum, `docs/benchmarks/m4b-healing.md`, architecture §16.4, plan results.
**Steps:**
1. Record pre-split behavior, kept/weld rule, acceptance values, deferrals (apex/pole to M4b-2, fillet/shell later).
2. Tick plan boxes, add Results (commits, counts).
3. Run fast-lane `-tc="*heal*"`; commit `docs(m4b): record the closed-rim and torus contract`.

## Exit criteria
- Closed-rim cylinder and torus cases heal to the stated volumes with zero boundary edges, `chi==2` for the cylinder/bands and `chi==0` for the closed torus halves, one shell, deterministic repeats.
- Failed heals return the original input bit-identically.
- Existing M4a/QA-12/negative tests unchanged in value and message.
- Task-level: Release targeted filters green with file-level counts; milestone: full Debug+Release suites green, then ff merge, push, CI green.

## Self-review
- [x] Each task has a RED test and a single commit.
- [x] No test asserts only counts; volumes/topology/messages are pinned.
- [x] Rejection paths are tested for input immutability.
- [x] Fast-lane commands are named.
- [x] Deferrals (apex/pole, fillet/chamfer, shell, NURBS, partial-u) are explicit.

## Results

| Task | Commit | Scope |
| --- | --- | --- |
| T1 | `3f72969` | split closed boundary rims on a private copy |
| T2 | `95bb910` | endpoint-directed welds + same-seam orientation normalization |
| T3 | `4936801` | cylinder closed-rim acceptance |
| T4 | `b0dd6b9` | torus band gaps (u, outer v, inner v seam) |
| T5 | `19773bc` | torus closed-rim halves + negatives |
| T6 | `docs(m4b): record the closed-rim and torus contract`（本提交） | ADR-0026 addendum, `docs/benchmarks/m4b-healing.md`, architecture §16.4, this plan, rejection lock |

Fast lane (Release `build-rel`, `geom_tests.exe`, `GEOM_WERROR=ON`, no warnings):

- `-tc="*heal*"`: 91 cases / 1,866 assertions, 0 failed.
- `-tc="*torus*"`: 59 cases / 10,754 assertions, 0 failed.
- `-tc="*heal*,*torus*"`: 145 cases / 12,440 assertions, 0 failed.

Pre-lock baseline at T5 was `*heal*` 90/1,847, `*torus*` 59/10,754, combined
144/12,421. The Task 6 lock adds exactly 1 case / 19 assertions, so the growth is
accounted for by the rejection test and no other test changed.

**Full Debug and Release suites were run at milestone close:** 1074 cases /
680,503 assertions each, 0 failed, 0 warnings, before ff merge/push and CI.

Corrections recorded during Task 6: torus halves are `chi == 0` (the brief and the
original Task 5 step said 2); Task 3's `V=6/E=8/F=4` is the post-split working copy
while the healed result is `V=4/E=6/F=4`; the fixture header is
`tests/unit/m4a_curved_heal_fixture.hpp` (not `brep_heal_curved_fixture.hpp`); and
on a post-split rejection the returned Brep is the caller's input while
`report.boundary_edges` is the post-split working count.
