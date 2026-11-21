# M3c-3c-4b3c-1: Twisted Torus Cuts by Chart Re-Origin

## Goal

Let a whole torus be cut by a plane that crosses the torus seam twice at different
transverse coordinates (the `z = 0.3 x` fixture), by moving the torus parameter
origin into a cut-free band instead of teaching the splitter to stitch seam
fragments. The public `boolean_solid_curved` path then produces the two halves with
the analytic volume `pi^2 / 4 = 2.4674011002723395`.

This is the first slice of M3c-3c-4b3c. The cone axis-plane cut and the cone apex
ring rewrite are a separate plan (4b3c-2), because the cone failure is not a
twisted crossing at all: SSI returns no curves for a plane containing the axis.

## Architecture

The torus chart is `u in [0, 2*pi]`, `v in [0, 2*pi]`, both periodic. A cut that
leaves the chart at one `u` on the `v = 0` seam and re-enters at another `u`
cannot be represented as one rectangle ring without an off-surface chord. Moving
the `v` origin so the seam lies in a cut-free band turns the same geometry into
three ordinary rectangle fragments. This is the chart re-origin accepted in
ADR-0032 and already used for sphere/sphere lenses (ADR-0030).

The re-origin is a property of the torus surface (`Torus3::v_origin`), and the
primitive rebuilds its two seam circles at the shifted origin. The splitter, the
assembler, the clipper and the healer keep their current contracts: they receive a
face whose chart seam does not cross the cut.

## Tech Stack

- C++20, MSVC 19.44 (`/W4 /WX /utf-8`) and GCC/Clang in CI.
- CMake 4.4.4, doctest 2.4.11, `GEOM_WERROR=ON`.
- Existing exact predicates (`sign_orient2d`, `ring_orientation2`) and the
  existing `ssi` pipeline.

## Spec

- ADR-0032 (general domain splitter) and its declaration that a twisted cut needs
  a chart re-origin, not a different ring.
- ADR-0030 (sphere lens re-chart) as the precedent for a single re-chart recursion.
- `docs/benchmarks/m3c3c4b3c.md` (to be created in Task 4) records the fixture and
  the measured volumes.

## Global Constraints

- No global epsilon; tolerances are scaled to the geometry.
- `Torus3::v_origin = 0` must be bit-identical to today for `eval`, `du`, `dv`,
  `normal`, `face_uv_project_torus`, `make_torus` output and every existing test.
- Exact predicates decide orientation; floating point only measures gaps.
- No silent empty result: every new refusal returns `Errc::unsupported` or
  `Errc::numeric_failure` with a message that names the limitation.
- The standalone `split_face_curved` API keeps refusing twisted closed chains; the
  new support is reached through the boolean re-chart, which presents a
  non-twisted chart.
- Only a whole torus (one shell, one face, no holes) can be re-charted. A torus
  that is one face of a larger solid returns the existing explicit refusal.
- MSBuild header dependency tracking is unreliable in this repository: after a
  header edit, touch a `.cpp` or pass `/t:Rebuild` before trusting a test run.

## Scope and deviations (for approval)

1. **Torus/plane only.** The re-chart handles a whole torus against a closed
   all-planar operand (the tilted slab fixture), plus the axis-perpendicular and
   axis-containing plane cases that already work. Torus/torus, torus/cylinder and
   curved/curved cutters stay deferred.
2. **At most two seam crossings.** A cut with more than four seam crossings, or
   with a free gap narrower than the margin, is refused deterministically. The
   current message is reworded from "M3c-3c-4" to "M3c-3c-4b3c".
3. **Direct twisted calls still fail.** `split_face_curved` on a raw twisted chain
   keeps the twisted refusal; only the boolean driver re-charts.
4. **Oblique volume tolerance 5e-3.** The re-charted SSI chain is polygonal with
   a fixed chord budget (`ssi.hpp`), so the tilted fixture asserts 5e-3 relative
   and 1e-2 on the sum, while the axis-aligned suite keeps 1e-9.
5. **Cone is out of scope.** Axis-containing planes and the apex ring closure are
   4b3c-2.

## Review Focus

1. `Torus3::v_origin` semantics: `eval(u, v)` equals the old surface at
   `(u, v + v_origin)`; projection returns the old angle minus `v_origin` without
   wrapping; zero origin is bit-identical.
2. `make_torus` with `v_origin`: both seam circles lie on the shifted chart seam
   (`v = 0`), and the loop pcurves project back to the same `(u, v)` points.
3. The gap helper rejects ambiguous coverage (no free gap, gap below margin,
   more than four crossings) and never returns a point inside a cut band.
4. The re-chart runs at most once, before face pairing, and leaves every
   non-torus case bit-identical.
5. The tilted fixture's `intersection + difference` equals the torus volume within
   1e-9, and `union` equals the slab plus that half; no orientation flip.
6. Refusal messages and `CurvedBooleanReport` counters are explicit; no silent
   empty solid.

## File Map

- `include/geom/surfaces/analytic.hpp` — `Torus3::v_origin`, eval/du/dv/normal.
- `include/geom/brep/face_uv.hpp` — `face_uv_project_torus` offset.
- `include/geom/brep/primitives.hpp` — `make_torus` optional `v_origin`.
- `include/geom/brep/trim_face.hpp` — gap helper next to
  `trim_face_closed_cut_twisted`; message reword.
- `include/geom/brep/boolean_solid_curved.hpp` — torus operand classification,
  re-chart driver at the entry, guard relaxation, face-info rebuild, report.
- `tests/unit/analytic_surface_test.cpp` — torus shift identity.
- `tests/unit/primitives_test.cpp` — `make_torus` shifted seam/volume.
- `tests/unit/brep_trim_face_test.cpp` — gap helper and message.
- `tests/unit/brep_boolean_curved_test.cpp` — tilted torus matrix; flips the old
  rejection case.
- `docs/adr/0037-torus-twisted-cuts.md`, `docs/benchmarks/m3c3c4b3c.md`,
  `docs/architecture/geometry-library-architecture.md` §16.3,
  `docs/superpowers/plans/2026-11-21-m3c3c4b3c1-torus-twisted.md`.

## Task 1: torus v-origin

**Files:** `include/geom/surfaces/analytic.hpp`,
`include/geom/brep/face_uv.hpp`, `include/geom/brep/primitives.hpp`,
`tests/unit/analytic_surface_test.cpp`, `tests/unit/primitives_test.cpp`.

**Interfaces:**

```cpp
struct Torus3 {
  Point3d center{};
  Vec3d axis{0.0, 0.0, 1.0};
  Vec3d x_axis{1.0, 0.0, 0.0};
  double major_radius = 2.0;
  double minor_radius = 0.5;
  double v_origin = 0.0;   // chart v = 0 is the surface curve v = v_origin
  // eval/du/dv/normal evaluate at (v + v_origin); is_valid() requires a finite
  // v_origin.
};

inline Result<Brep> make_torus(Point3d center, Vec3d axis, double major_radius,
                               double minor_radius, double tol,
                               double v_origin = 0.0);
```

`face_uv_project_torus` returns `atan2(axial, q - R) - v_origin` (no wrapping;
`face_uv_read_loop` keeps its existing unwrap).

**Steps:**

1. Write the failing tests:
   - zero-origin bit identity for `eval`, `du`, `dv`, `normal`: a `Torus3` with
     `v_origin = 0` equals the existing closed-form expressions written out in
     the test, coordinate for coordinate;
   - shift identity: `shifted.eval(u, v) == base.eval(u, v + 0.7)` for a grid of
     `u, v`, same for `du`/`dv`, and `normal` equals the base normal;
   - `make_torus(..., tolerance, 0.7)`: `validate_topology` ok, every vertex is on
     the shifted surface within 1e-12, `curved_solid_volume` equals
     `4.934802200544679` within 1e-9 (volume is origin-independent);
   - `make_torus(..., tolerance, 0.0)` keeps the existing count and mesh hash.
2. Run the tests; expect compile failure on `v_origin` and on the six-argument
   `make_torus`.
3. Implement the field and the primitive: the `u_edge` seam circle sits at old
   `v = v_origin`, with `radius = R + r*cos(v_origin)` and centre
   `center + axis*(r*sin(v_origin))`; the `v_edge` tube circle stays at `u = 0`.
   Both pcurves come from the existing `add_loop(..., true)` projection.
4. Rebuild and run the focused tests, then the full Debug and Release suites.

**Acceptance:** `Torus3{}` and five-argument `make_torus` calls compile unchanged;
zero-origin output is bit-identical; shifted torus volume is unchanged.

## Task 2: twist-offset helper

**Files:** `include/geom/brep/trim_face.hpp`,
`tests/unit/brep_trim_face_test.cpp`.

**Interface (detail):**

```cpp
// Midpoint of the widest free interval of [0, period), given covered intervals.
// Returns nullopt when there is no free interval, when the widest free interval
// is shorter than 2*margin, or when more than max_intervals intervals are given.
std::optional<double> trim_face_periodic_free_midpoint(
    std::span<const std::array<double, 2>> covered, double period,
    double margin, std::size_t max_intervals = 4);
```

**Steps:**

1. Failing tests: covered `{[2.82082, 3.46237], [5.37949, 0.90369]}` with period
   `2*pi`, margin `1e-3` returns `1.862255` (first free band midpoint) or the
   other band midpoint `4.420925`; deterministic choice: the widest band, ties
   broken by the lower midpoint. No free gap (full coverage) -> nullopt; a gap of
   `1e-4` with margin `1e-3` -> nullopt; five intervals -> nullopt.
2. Implement with sorted, normalized endpoints; merge overlaps; compute the
   circular complement; pick the widest complement interval.
3. Reword the twisted refusal message to
   `"split_face_curved: closed cut crosses the seam at different transverse coordinates; twisted seam crossings are deferred to M3c-3c-4b3c"` and update the
   existing assertions.

**Acceptance:** helper returns a point at least `margin` away from every covered
interval; message test updated.

## Task 3: boolean re-chart for a twisted torus

**Files:** `include/geom/brep/boolean_solid_curved.hpp`,
`tests/unit/brep_boolean_curved_test.cpp`.

**Interfaces (detail):**

```cpp
struct CurvedBooleanRebasedTorus {
  Brep brep;
  SolidId solid = kNoEntity;
};

// Returns nullopt when the pair is not a whole-torus/plane pair that the current
// guards would reject as twisted. Returns an Error when a re-chart was attempted
// and failed.
Result<std::optional<CurvedBooleanRebasedTorus>> curved_boolean_rechart_twisted_torus(
    const Brep& torus_brep, SolidId torus_solid, const Brep& other,
    SolidId other_solid, double tolerance);
```

The driver runs after the sphere re-chart and before face pairing. It classifies
the torus operand as a whole torus (one shell, one face, a single `Torus3`, no
holes) and the other operand as all-planar. It collects the SSI cut chains with
the existing `ssi` + `curved_boolean_cut_chains` path, covers the `v` coordinates
(merged within `1e-9`), calls the gap helper, and when a delta exists rebuilds the
torus with `make_torus(center, axis, R, r, tol, v_origin + delta)`. It then reruns
the boolean once (single re-chart guard, mirroring
`curved_boolean_rechart_whole_spheres`). If the re-charted cut is still twisted,
return the existing refusal.

Guard changes: `curved_boolean_torus_plane_cut_supported` stays for the
perpendicular and axis-containing families; the tilted family is admitted only
through a successful re-chart. The refusal message for an unrecoverable tilt
remains the current torus guard message.

**Steps:**

1. Write failing tests in `tests/unit/brep_boolean_curved_test.cpp`:
   - replace the body of `curved boolean: torus cut by a tilted plane is rejected`
     with the success matrix below; keep an offset-parallel rejection case;
   - a miss/tangent plane still returns `unsupported` with the torus guard message.
2. Expected values (torus `R=1`, `r=0.5`, slab bottom `z = 0.3x`, `tol=1e-9`):

```cpp
constexpr double kHalfTorus = 2.4674011002723395;  // pi^2 / 4
constexpr double kSlab = 256.0;                    // 8 x 8 x 4 tilted slab
constexpr double kTorus = 4.934802200544679;       // pi^2 / 2
// intersection, difference: kHalfTorus (5e-3 relative)
// union: kSlab + kHalfTorus (5e-3 relative)
// intersection + difference == kTorus within 1e-2
```

   Each result: `validate_topology(...).ok()`, one shell, `closed`, two runs
   bit-identical, `report.face_pairs > 0`.
3. Implement the operand classification, cut collection, delta selection and the
   single recursive retry. Any failure leaves the previous guard result.
4. Run focused tests, then the full suites.

**Acceptance:** the tilted fixture succeeds with the analytic half volume; the
axis-aligned torus matrix (offset plane, union/intersection/difference) is
bit-identical; non-whole torus operands still refuse explicitly.

## Task 4: contract and evidence

**Files:** `docs/adr/0037-torus-twisted-cuts.md`,
`docs/benchmarks/m3c3c4b3c.md`,
`docs/architecture/geometry-library-architecture.md`,
`docs/superpowers/plans/2026-11-21-m3c3c4b3c1-torus-twisted.md`.

**Steps:**

1. ADR-0037: the seam-crossing problem, the re-origin decision, why stitching was
   rejected (the probe's 1.0-long off-surface chords), the single re-chart rule,
   the refusal class, and the deferred torus/torus and cone work. Note the splitter
   message change from 4b3 to 4b3c.
2. Benchmark note: fixture, tolerances, measured volumes for intersection,
   difference and union, and the fact that the direct splitter still refuses raw
   twisted chains.
3. §16.3: replace the "tilted torus planes deferred" sentence with the supported
   class and its limits; keep the M3c-3c-4b3c forward pointer for cone work.
4. Tick the plan checkboxes and record the two-configuration counts.

**Acceptance:** documents match the tests; no over-claim for torus/torus, tangent
cuts, multiple twists, or cone.

## Exit criteria

- Debug and Release: all tests pass, zero warnings under `GEOM_WERROR=ON`.
- Tilted torus intersection and difference are `2.4674011002723395`, union is
  `258.46740110027233`, each within 5e-3 relative; `intersection + difference`
  equals the torus volume within 1e-2.
- `v_origin = 0` is bit-identical across the surface, projection, primitive and
  all existing booleans.
- The offset-parallel and miss/tangent refusals are unchanged in code and message.
- The full branch is reviewed, merged fast-forward to `master`, pushed, and the
  six-job CI matrix is green.

## Self-review

- [ ] Every task has files, interfaces, five TDD steps and acceptance.
- [ ] No task depends on an interface introduced later.
- [ ] Tolerance choices are justified (5e-3/1e-2 for the oblique SSI chain, 1e-9 for
      exact invariants).
- [ ] Refusals are explicit; no test weakens an assertion to pass.
- [ ] Scope deviations are listed for approval.
- [ ] No placeholder or self-correction note remains in this plan.
