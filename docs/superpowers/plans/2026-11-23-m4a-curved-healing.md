# M4a: Curved-Gap Healing Implementation Plan

## Goal

Extend `heal_brep` from planar, line-only shells to analytic curved shells: a
cylinder split into two bands, a sphere split into lunes, and a cone band ending
at its apex must heal across a sub-micron gap into a closed, valid solid whose
volume matches the input geometry, while mismatched carriers, off-surface rims,
unsupported geometry, and closed rims are rejected deterministically.

Scope: analytic surfaces `Plane3/Cylinder3/Cone3/Sphere3/Torus3` and edge curves
`LineCurve3/CircleCurve3`. NURBS surfaces/curves remain rejected with explicit
messages. Torus gap healing and apex/pole gap fixtures are deferred to M4b.

## Architecture

`heal_brep` stays a whole-solid transaction: plan matching (weld vertices, pair
edges), apply welds, merge edges, rebuild pcurves, close shells, raise
tolerances, validate. The existing planar path must remain bit-identical.

Three model changes:

1. Kind whitelists replace the `Plane3` / `LineCurve3` rejections.
2. `heal_plan` gains an orientation-aware tie-break for curved seams, where two
   candidates have identical gaps (two arcs of the same circle, or a sphere
   lune's two meridians). The tie-break is gated to non-planar pairs so planar
   ambiguity tests keep reporting `ambiguous_gap`.
3. Pcurves are rebuilt per face in loop order with periodic-u/v unwrapping and
   interior circle samples; the 3D curve is preserved (no
   `set_edge_geometry(kNoGeometry)`).

Matching must also prove the two faces carry the same analytic surface and that
each gap edge lies on its face's surface; otherwise welding silently absorbs a
carrier mismatch.

## Tech Stack

C++20, CMake 4.4.4, MSVC 19.44 (VS2022), doctest, PowerShell. No new
dependencies. `/W4 /WX /utf-8`.

## Spec

- `docs/architecture/geometry-library-architecture.md` §4.6 (tolerance
  invariants), §8.12 (healing), §16.4.
- `docs/adr/0026-healing-stitch.md` and its M4a addendum.
- `.superpowers/sdd/m4a-heal-probe-report.md`: measured prototype behaviour,
  fixture values, and dead ends.

## Global Constraints

- Task-level gate (fast lane, m30759/m30761): build and run only the touched
  translation units plus targeted doctest filters, Release first; run the full
  Debug and Release suites only at milestone review. Reviewers are reused within
  the milestone; they run targeted checks and trust implementer full-suite
  counts.
- Every task starts with a RED doctest on the current tree, then the minimal
  implementation, then targeted GREEN. One commit per task.
- Do not weaken the existing planar tests. The four planar ambiguity
  expectations (`make_ambiguous_rims`) must keep reporting `ambiguous_gap`.
- Exact predicates are used for topology decisions; doubles only measure gaps.
  No new global epsilon.
- Header changes invalidate dependents on this build: after editing `heal.hpp`
  delete `build/tests/geom_tests.dir/<Config>` or run `/t:Rebuild` before
  trusting a result (stale-COMDAT trap in the probe report).
- No silent success: every rejection returns the documented `Errc` and message
  with `changed=false` and the input unchanged.

## Scope and Deviations (need approval)

1. Minimum slice = Tasks 1-5. Task 6 is diagnostics/docs plus a cone-apex
   coverage case; torus gap healing is deferred to M4b because its v periodicity
   and two radii are unproven.
2. Tie-break gate is "at least one adjacent face surface is not `Plane3`", not
   "curve is `CircleCurve3`". The surface gate matches the real discriminator and
   keeps planar circular rims unchanged; a planar-rim behaviour change would need
   its own test and is out of scope.
3. Coincidence tolerance is `max(options.tolerance, 1e-12 * scale)` where
   `scale` is the largest of radius, centre distance, and origin magnitude. A
   one-ulp radius difference is accepted; a 2e-7 radius difference is rejected.
4. Closed gap rims (a single circle edge with `v0 == v1`, as `make_cylinder`
   emits) are not healed: they get the explicit
   `Errc::unsupported, "heal_brep: closed gap edge cannot be paired; split the
   rim first"` instead of `invalid_topology "heal_plan: boundary edge has zero
   length"`. Splitting inside `heal_brep` is prohibited; `collapse_circular_chains`
   is the M4b building block.
5. Apex/pole gap (a boundary shrunk to the apex) stays out of scope; the cone
   fixture ends the band exactly at the apex.

## Review Focus

1. **Planar parity**: box, cracked box, torn sheet, and the two
   `make_ambiguous_rims` expectations are unchanged, including messages.
2. **Carrier rejection**: cylinder r=1 vs r=1+2e-7, plane lid on a cylinder, and
   rims floating 2e-7 off an r=1 cylinder are rejected before mutation with
   carrier / edge-on-surface messages and context (keep edge, drop edge).
3. **Tie-break gating**: a curved seam resolves; a planar tie still reports
   `ambiguous_gap`; a tie whose candidates all share direction still reports
   ambiguity (do not turn ambiguity into a wrong weld).
4. **Pcurve continuity**: after a sphere heal, every lune's uv area is positive
   and its u band does not double across the seam; cylinder/cone charts close;
   the 3D circle curve survives.
5. **Acceptance values**: cylinder `pi*(2+gap)` rel 1e-12 and `2*pi` rel 1e-6;
   sphere `4*pi/3` rel 1e-12; cone `2*pi/3` rel 1e-12; `chi == 2`;
   report counts exactly as measured.

## File Map

- `include/geom/brep/heal.hpp` — guards, tie-break, carrier/edge checks,
  per-face pcurve rebuild, closed-rim diagnostic, report fields.
- `tests/unit/m4a_curved_heal_fixture.hpp` — new builders: two-band cylinder,
  shifted-chart cylinder, four-lune sphere, cone band ending at the apex,
  negatives (radius 2e-7 mismatch, off-surface rims, plane lid, closed rims,
  1e-3 mismatch).
- `tests/unit/brep_heal_curved_test.cpp` — new acceptance/negative tests.
- `tests/unit/brep_heal_test.cpp` — update the three old-message expectations.
- `tests/CMakeLists.txt` — register the new test file.
- `docs/adr/0026-healing-stitch.md`, `docs/benchmarks/m4a-healing.md`,
  architecture §16.4 — contract and evidence.
- `docs/superpowers/plans/2026-11-23-m4a-curved-healing.md` — this plan.

## Task 1: Analytic kind whitelist

**Files**: `include/geom/brep/heal.hpp` (guards at the current :663-704),
`tests/unit/m4a_curved_heal_fixture.hpp` (new), `tests/unit/brep_heal_curved_test.cpp`
(new), `tests/unit/brep_heal_test.cpp` (three expectations), `tests/CMakeLists.txt`.

**Interfaces**: no public API change.

**Step 1 — failing tests**:

```cpp
TEST_CASE("heal: a NURBS face is rejected as an unsupported surface kind") {
  // Build a closed box-like solid whose one face carries a flat NurbsSurface
  // (control net from brep_coverage_test.cpp:1350); the other faces are planes.
  const auto result = geom::v1::heal_brep(brep, solid, options);
  REQUIRE_FALSE(result);
  CHECK(result.error().code == geom::v1::Errc::unsupported);
  CHECK(std::string_view{result.error().message} ==
        "heal_brep: surface type is not supported yet");
  CHECK(result.error().context_a == nurbs_face.index());
}
```

Add the companion NURBS-curve test with
`"heal_brep: curve type is not supported yet"`, and a cylinder-gap test that
asserts the result is not the old `"curved geometry is not supported yet"`
message (at this task it may still report `ambiguous_gap`; assert only on the
absence of the old message).

**Step 2 — run**: `build-rel` Release `geom_tests --test-case="heal:*"`; expect
the new cases red (old guard message) and the three old expectations green.

**Step 3 — minimal implementation**:

- Faces: accept `Plane3`, `Cylinder3`, `Cone3`, `Sphere3`, `Torus3`; otherwise
  `Error{Errc::unsupported, "heal_brep: surface type is not supported yet",
  face.index(), 0}`.
- Edges: accept `kNoGeometry`, `LineCurve3`, `CircleCurve3`; otherwise
  `"heal_brep: curve type is not supported yet"`.
- Keep structural checks (alive references, shared-boundary-edge check).

**Step 4 — green**: targeted filters green; update the three
`brep_heal_test.cpp` expectations (`:1055`, `:1330`, `:1352`) to the new
messages; no planar behaviour change.

**Step 5 — commit**: `feat(heal): allow analytic surfaces and circular edges`.

## Task 2: Orientation-aware tie-break

**Files**: `include/geom/brep/heal.hpp` (`heal_plan` candidate selection),
`tests/unit/brep_heal_curved_test.cpp`, fixture.

**Interfaces**: none public. Keep `ambiguous_gap` for genuine ambiguity.

**Step 1 — failing tests**:

```cpp
TEST_CASE("heal: a cylinder radial gap resolves by seam orientation") {
  const auto [brep, solid] = make_two_band_cylinder(1.0, 1e-6);
  const auto healed = geom::v1::heal_brep(brep, solid, options);
  REQUIRE(healed);
  CHECK(healed.value().report.changed);
  CHECK(healed.value().report.boundary_edges == 4u);
  CHECK(healed.value().report.welded_vertices == 2u);
  CHECK(healed.value().report.sewn_edges == 2u);
}
```

Add the sphere lune version. On the current tree both report `ambiguous_gap`,
`changed = false`.

**Step 2 — run** the two cases; expect red.

**Step 3 — implementation**:

- When the best two candidates are within `kAmbiguityGap = 1e-12` of the best
  gap, and at least one adjacent face surface is not `Plane3`, pick the
  candidate whose `same_direction` flag is `false` (the true counterpart
  traverses the shared seam in the opposite direction). If all tied candidates
  share the same direction, keep the existing `ambiguous_gap` issue.
- Planar pairs keep the old behaviour exactly.

**Step 4 — green**: targeted tests green; run the full `heal:*` filter and
confirm 40/43-style planar cases (except Task 1's updated messages) still pass
and `make_ambiguous_rims` still yields four `ambiguous_gap` issues.

**Step 5 — commit**: `feat(heal): resolve curved seam ties by orientation`.

## Task 3: Carrier and edge-on-surface acceptance

**Files**: `include/geom/brep/heal.hpp` (helpers + checks before welding),
fixture negatives, tests.

**Interfaces** (detail):

```cpp
bool heal_surfaces_coincide(const SurfaceGeometry& a, const SurfaceGeometry& b,
                            double tol);
double heal_surface_residual(const SurfaceGeometry& s, const Point3d& p);
bool heal_edge_lies_on_surface(const Brep& brep, const Edge& edge,
                               const SurfaceGeometry& surface, double tol);
```

**Step 1 — failing tests** (exact messages/context):

| fixture | expected |
|---|---|
| upper band radius `1 + 2e-7` | `unsupported "heal_brep: gap faces do not share the same surface"`, context (keep edge, drop edge), input unchanged |
| plane lid on a cylinder wall | same carrier message |
| gap rims at `1 + 2e-7` on an `r = 1` cylinder | `unsupported "heal_brep: gap edge is not on its face surface"` |
| kept curve not through welded vertices | `unsupported "heal_brep: kept edge curve does not match its welded vertices"` |
| radius mismatch `1e-3` (> max_gap) | still a `gap_too_large` report, not an error |

**Step 2 — run**: red (some currently weld and fail later with
`face_loops_uv: point is off its circle`).

**Step 3 — implementation**:

- Compare surfaces by kind: planes always coincide (historical); cylinders need
  equal radius, parallel axes, perpendicular origin offset; spheres need centre
  distance and radius; cones need radius, half angle, parallel axes, apex
  distance; tori need centre, parallel axes, both radii; mixed kinds are false.
- Tolerance `max(options.tolerance, 1e-12 * scale)` where scale is the largest
  radius/centre/origin magnitude.
- Edge-on-surface: sample `t0`, midpoint, `t1` on the edge curve and test the
  implicit residual of the face surface; require each `<= tol`.
- Kept-curve check: both endpoints of the kept edge within `tol` of the curve
  (circle radius / line distance) after welding.
- Run these after plan matching, before any `SolidEditor` mutation; rejection
  returns the input Brep with `changed=false`.

**Step 4 — green**: four negatives exact; positives still proceed.

**Step 5 — commit**: `feat(heal): require coincident carriers and on-surface rims`.

## Task 4: Per-face loop-order pcurve rebuild

**Files**: `include/geom/brep/heal.hpp` (replace the per-edge rebuild),
fixture/test.

**Interfaces** (detail):

```cpp
Result<void> heal_rebuild_face_pcurves(SolidEditor& editor, FaceId face,
                                       double tolerance);
```

**Step 1 — failing tests**:

- After sphere heal, every lune's `face_loops_uv` area is positive; the south
  west lune's u band is `[-pi, 0]` (not `[0, 2pi]`); total volume is
  `4*pi/3` rel 1e-12.
- After cylinder/cone heal, every face projects through `face_loops_uv` and the
  cap charts have 64 points, the lateral chart 130 points.
- The 3D circle/line curve survives: `edge.curve != kNoGeometry` after healing.

**Step 2 — run**: red (`curved_solid_volume: face triangulation failed` for the
sphere; `face_loops_uv: circular arc direction is ambiguous` for a two-point
pcurve).

**Step 3 — implementation**:

- Mark a face dirty when it owns a coedge of an edge touching a welded vertex.
- For each dirty face: resolve `face_uv_projection` once; walk its loops in
  coedge order; carry `previous_u` across junctions; sample the edge curve
  (endpoints from vertices, 31 interior points for circles, 2 for lines);
  reverse samples for a reversed coedge; project with `face_uv_project` and
  unwrap u/v for periodic surfaces; build `Pcurve2d` in edge direction; store
  on the coedge.
- Do not call `set_edge_geometry`.
- Keep planar pcurves (2 points per line) rebuilt by the same path.

**Step 4 — green**: sphere/cylinder/cone charts valid; planar heal and QA-12
tests unchanged.

**Step 5 — commit**: `feat(heal): rebuild face pcurves in loop order`.

## Task 5: Curved acceptance, volume and report

**Files**: `tests/unit/brep_heal_curved_test.cpp`,
`include/geom/brep/heal.hpp` only if a report field needs adjustment.

**Step 1 — failing tests** (after Tasks 1-4 they should pass; keep them as the
acceptance gate):

| case | volume | tolerance |
|---|---|---|
| cylinder r=1, gap 1e-6 | `pi*(2+gap)` = 6.2831884487723997 | rel 1e-12 |
| same, `2*pi` | 6.283185307179586 | rel 1e-6 |
| sphere r=1, gap 1e-6 | 4.189... (four lunes) | rel 1e-12 |
| cone r=1, h=2, band to apex | 2.0943951023931864 | rel 1e-12 |

Also assert `validate_topology(...).ok()`, `euler_characteristic == 2`,
`boundary_edges == 4`, `welded_vertices == 2`, `sewn_edges == 2`,
`max_closed_gap` equals the measured gap, `max_tolerance <= 1e-6`,
`changed == true`, and a second run is idempotent.

**Step 2 — run**: green if Tasks 1-4 are correct; any failure is a real defect.

**Step 3 — fix only what fails**; do not relax assertions.

**Step 4 — full targeted regression**: `heal:*`, `curved volume:*`,
`tessellate:*`, `validate:*` in Release.

**Step 5 — commit**: `test(heal): accept analytic curved gaps`.

## Task 6: Cone apex coverage, closed-rim diagnostic, docs

**Files**: fixture, tests, `heal.hpp` (closed-rim message), ADR-0026,
`docs/benchmarks/m4a-healing.md`, architecture §16.4.

**Steps**:

1. Red: a closed rim edge (cylinder with closed rims) returns
   `unsupported "heal_brep: closed gap edge cannot be paired; split the rim
   first"` instead of `heal_plan: boundary edge has zero length`.
2. Implement the closed-edge check before the zero-length guard.
3. Add the cone apex band as a permanent acceptance case (already green);
   document that an apex gap is out of scope.
4. ADR addendum: whitelists, carrier predicate, orientation tie-break (gated to
   non-planar), per-face pcurve contract, closed-rim limitation, torus/apex
   deferral.
5. Benchmark note: volume/face counts for the three accept cases, no time gate.
6. Commit: `docs(heal): record the curved-gap contract`.

## Exit Criteria

- Cylinder, sphere, cone accept cases heal: closed, `chi == 2`, exact volumes,
  deterministic.
- Four negative classes reject before mutation with the exact messages above.
- Planar heal suite and QA-12 unchanged; no new global epsilon.
- Full Debug and Release suites green at milestone review, zero warnings.
- ADR/architecture/benchmark document the supported class and deferrals.

## Self-review

- [ ] Each task has a RED test that fails for the intended reason.
- [ ] No test asserts only counts where a value/message is available.
- [ ] Tie-break cannot turn a genuine ambiguity into a weld.
- [ ] Rejections leave the input bit-identical.
- [ ] Fast-lane commands are given; milestone full-suite run is explicit.
- [ ] Torus/apex/NURBS deferrals are stated, not implied.
