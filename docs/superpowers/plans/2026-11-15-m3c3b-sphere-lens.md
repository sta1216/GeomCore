# M3c-3b: Sphere/Sphere Lens Booleans — Implementation Plan

## Goal

Make `boolean_solid_curved` support two **whole-sphere** operands whose boundaries
cross: `union`, `intersection` and `difference` must produce a closed, validated
2-manifold whose volume matches the closed-form lens values, with no change to
the seam-only splitter.

## Architecture

Two spheres intersect in a circle. Re-charting each operand sphere so that its
`axis` points at the other centre turns that circle into a **latitude**, which the
already-shipped M3c-3a cap/pole path splits exactly. `make_sphere` gains an
axis/frame overload; a helper rebuilds whole-sphere operands in the shared chart;
the driver recurses once on the rebuilt copies. An analytic sphere/sphere circle
replaces the missing `pair_circles` branch, and the flipped-fragment harvest
accepts a reversed spherical cap so `difference` contributes negative volume.

## Tech Stack

- C++20, header-only library, existing `Brep`/`SolidEditor`/SSI/trim/assemble
  pipeline.
- CMake 4.4.4 + VS2022 (MSVC 19.44) on Windows, GCC/Clang on CI.
- doctest 2.4.11.

## Spec

- `docs/architecture/geometry-library-architecture.md` v2.1, §8.10 (tessellation
  contract), §8.12 (boolean contracts), §16.3 (milestone status).
- `docs/quality/testing-strategy.md` (five test kinds, differential corpus,
  invariants).
- Probe evidence: `.superpowers/sdd/m3c3b-probe-report.md` (not committed).

## Global Constraints

- No new third-party dependencies.
- `-DGEOM_WERROR=ON`: `/W4 /WX /utf-8` on MSVC, `-Wall -Wextra -Werror
  -Wshadow` on GCC/Clang. Every changed translation unit must build clean.
- No global epsilon; all comparisons are tolerance-scaled
  (`tolerance` passed by the caller, default 1e-9).
- Topological decisions use exact predicates where the code already does;
  volume/report comparisons use relative tolerances.
- Determinism: the same inputs must produce bitwise-identical volumes on
  repeated calls in one process.
- No silent degradation: every unsupported input returns `Errc::unsupported`
  with a message naming the limitation; no partial or invalid solid is returned.
- Public headers must stay self-contained and installed by the existing rules.
- Every task commits once, after Debug and Release suites pass.

## Scope and Deviations

1. **Whole-sphere operands only.** A sphere face that is one face of a larger
   solid is rejected with
   `"boolean_solid_curved: sphere/sphere booleans need whole-sphere operands"`.
   Re-charting a partial solid would move the seam of a face whose neighbours
   are already charted.
2. **Lens cases only.** The two spheres must cross:
   `|ra - rb| < d < ra + rb` (relative to `tolerance`). Containment, tangency and
   separation return
   `"boolean_solid_curved: sphere/sphere booleans need two crossing spheres"`.
   Containment is representable in principle (the result is one operand), but the
   difference case needs a spherical cavity and that is deferred with the rest of
   multi-shell curved work.
3. **Output chart differs from the input.** The output sphere faces carry the
   re-charted `axis`/`x_axis`. Geometry and volume are unchanged; face identity
   and pcurves are not preserved. This is documented in ADR-0030 and is the
   reason no test may compare face surface frames with the inputs.
4. **Split only here.** Cone apex topology (~250–400 LOC) and offset/tilted torus
   v-windows (~120–180 LOC) are separate plans (M3c-3c, M3c-3d). This plan does
   not touch the splitter, cone, torus, NURBS, or general closed-cut support.
5. **No performance gate.** QA-16 is a correctness corpus. The e2e runtime is
   recorded in `docs/benchmarks/m3c3b.md` for information, with no threshold.

## Review Focus

Reviewers should try to break these five claims:

1. **Re-chart transparency.** The rebuilt operand differs only in `axis`/`x_axis`;
   its `face_loops_uv`, classification mesh, healed rim and assembled shell are
   consistent, and the output solid is closed and validated.
2. **Whole-sphere precondition.** A sphere face inside a box/cylinder solid, a
   solid with two sphere faces, or a sphere with an extra face is rejected — not
   silently re-charted.
3. **Near and exact degeneracies.** `d = ra + rb` (external tangency),
   `d = |ra - rb|` (internal tangency), `d` within `tolerance` of either, and
   `d` underflow/overflow must return the documented error, never a thin lens
   with a wrong volume or a division by zero.
4. **Orientation.** `difference` must subtract the B cap: `A - B` and `B - A`
   must have different volumes, both negative-correct, with one shared rim edge
   used in opposite directions.
5. **Radii and placement generality.** Unequal radii and centres off every axis
   work; the shared `x_axis` choice is deterministic; the seam-crossing point is
   the same 3D point in both charts.

## File Map

- `include/geom/brep/primitives.hpp` — `make_sphere` axis/frame overload.
- `include/geom/brep/boolean_solid_curved.hpp` — sphere/sphere circle, rebase
  helper and driver hook, flipped spherical cap.
- `tests/unit/primitives_test.cpp` — T2 unit tests (already registered).
- `tests/unit/brep_boolean_curved_test.cpp` — T1/T3/T4/T5 tests (already
  registered).
- `docs/adr/0030-sphere-lens-booleans.md` — new ADR.
- `docs/benchmarks/m3c3b.md` — QA-16 evidence.
- `docs/architecture/geometry-library-architecture.md` §16.3 — status.

## Task 1: Analytic sphere/sphere intersection circle

**Files**
- Modify `include/geom/brep/boolean_solid_curved.hpp` (`curved_boolean_pair_circles`,
  currently `:1474-1546`).
- Modify `tests/unit/brep_boolean_curved_test.cpp`.

**Interfaces**
- Consumes: `Sphere3`, `PlaneFrame`, `CircleCurve3`,
  `detail::curved_boolean_pair_circles(const SurfaceGeometry&, const PlaneFrame&,
  const SurfaceGeometry&, const PlaneFrame&, double) noexcept`.
- Produces: same signature; a sphere/sphere pair now returns one `CircleCurve3`
  when the boundaries cross, and `{}` otherwise.

**Algorithm**

For spheres `A(center ca, radius ra)` and `B(center cb, radius rb)`, let
`n = unit(cb - ca)` and `d = |cb - ca|`:

```
t   = (d*d + ra*ra - rb*rb) / (2*d)      // distance from ca along n
rho = sqrt(max(0, ra*ra - t*t))          // circle radius
center = ca + n * t
```

Guards, all relative to `tolerance`: `d <= tolerance`; `d >= ra + rb - tolerance`;
`d <= |ra - rb| + tolerance`; `rho*rho <= tolerance*tolerance` → return `{}`.
Pick the circle frame deterministically: `x_axis = unit(cross(n, candidate))`
with `candidate = {1,0,0}` unless `|n[0]| >= 0.9`, then `{0,1,0}`;
`y_axis = cross(n, x_axis)`. Keep the existing torus branch untouched.

**Step 1 — failing test**

Add to `tests/unit/brep_boolean_curved_test.cpp`:

```cpp
TEST_CASE("curved boolean: sphere/sphere has one analytic circle") {
  geom::SurfaceGeometry a = geom::Sphere3{{0.0, 0.0, 0.0}, {0, 0, 1}, {1, 0, 0}, 1.0};
  geom::SurfaceGeometry b = geom::Sphere3{{1.5, 0.0, 0.0}, {0, 0, 1}, {1, 0, 0}, 1.0};
  const std::vector<geom::CircleCurve3> circles =
      geom::detail::curved_boolean_pair_circles(a, geom::PlaneFrame{}, b,
                                                geom::PlaneFrame{}, 1e-9);
  REQUIRE(circles.size() == 1);
  CHECK(circles[0].center[0] == doctest::Approx(0.75));
  CHECK(circles[0].radius == doctest::Approx(std::sqrt(7.0) / 4.0));
  CHECK(circles[0].is_valid());
}
```

Plus `{}` cases: `d = 2.0` (tangent), `d = 0.5` (one inside the other),
`d = 2.5` (separate).

**Step 2 — run and confirm failure**

`curved_boolean_pair_circles` returns `{}` for sphere/sphere today, so
`REQUIRE(circles.size() == 1)` fails. Record the failing assertion.

**Step 3 — implement**

Add the branch before the torus branch, after the plane branch's early return.
No allocation is needed beyond the returned vector.

**Step 4 — run tests**

Focused case and the full suite pass in Debug and Release.

**Step 5 — commit**

`feat(m3c3b): add the analytic sphere/sphere intersection circle`

## Task 2: Frame-aware sphere primitive

**Files**
- Modify `include/geom/brep/primitives.hpp` (`make_sphere`, `:481-536`).
- Modify `tests/unit/primitives_test.cpp`.

**Interfaces**
- Produces:
  `Result<Brep> make_sphere(Point3d center, Vec3d axis, Vec3d x_axis, double radius, double tol);`
- Consumes: existing `make_sphere(center, radius, tol)` behaviour.

**Algorithm**

The existing body already declares `const Vec3d axis{0,0,1}; const Vec3d
x_axis{1,0,0};` locally. Parameterise it:

- validate finite `center`, `axis`, `x_axis`, `radius`, `tol` → existing error
  codes/messages;
- `radius > 0` → `Errc::degenerate "make_sphere: radius must be positive"`;
- `axis` unit length within `detail::kAnalyticSurfaceUnitTolerance` →
  `"make_sphere: axis must be unit length"`;
- `x_axis` unit and `|dot(axis, x_axis)| <= kAnalyticSurfaceUnitTolerance` →
  `"make_sphere: x-axis must be perpendicular to the axis"`;
- keep the pole finite/zero-length seam checks, the `CircleCurve3{center,
  x_axis, axis, radius}` seam, the `Sphere3{center, axis, x_axis, radius}` face,
  the two-coedge loop and the shell/solid construction.

The old overload delegates:

```cpp
return make_sphere(center, Vec3d{0.0, 0.0, 1.0}, Vec3d{1.0, 0.0, 0.0}, radius, tol);
```

**Step 1 — failing test**

```cpp
TEST_CASE("make_sphere accepts an arbitrary chart") {
  const auto sphere = geom::make_sphere({0.0, 0.0, 0.0}, {1.0, 0.0, 0.0},
                                        {0.0, 1.0, 0.0}, 1.0, 1e-9);
  REQUIRE(sphere.has_value());
  const geom::Face& face = sphere.value().face(sphere.value().face_ids().front());
  const geom::Sphere3* surface = std::get_if<geom::Sphere3>(
      &sphere.value().geometry().surfaces[face.surface]);
  REQUIRE(surface != nullptr);
  CHECK(surface->axis[0] == doctest::Approx(1.0));
  CHECK(surface->x_axis[1] == doctest::Approx(1.0));
}
```

Compilation fails (no such overload) — record the compiler error as the RED
evidence.

**Step 2 — implement and run**

Debug/Release: new case plus the old canonical sphere case (volume 4π/3, 2
vertices, 1 edge, 1 face) passes.

**Step 3 — commit**

`feat(m3c3b): add a frame-aware sphere primitive`

## Task 3: Re-chart whole-sphere operands and recurse

**Files**
- Modify `include/geom/brep/boolean_solid_curved.hpp`.
- Modify `tests/unit/brep_boolean_curved_test.cpp`.

**Interfaces**

```cpp
namespace detail {

struct CurvedBooleanRebasedSpheres {
  Brep a;
  SolidId a_solid;
  Brep b;
  SolidId b_solid;
};

// Returns nullopt when the operands are not a sphere/sphere pair, or when they
// are already charted with the shared axis. Returns an error for a sphere/sphere
// pair whose operands are not whole spheres or whose boundaries do not cross.
Result<std::optional<CurvedBooleanRebasedSpheres>>
curved_boolean_rebase_sphere_pair(const Brep& a, SolidId a_solid, const Brep& b,
                                  SolidId b_solid, double tolerance);

}  // namespace detail
```

- Consumes: `Brep::face_ids/edge_ids/vertex_ids`, `solid_shells`, `shell_faces`,
  `make_sphere(center, axis, x_axis, radius, tol)` (Task 2),
  `curved_boolean_face_info(...).synthetic_domain`.
- Produces: the struct above.

**Algorithm**

1. Classify each operand as a whole sphere: exactly one solid, one shell, one
   face, one edge, two vertices, the face is `Sphere3`, and
   `curved_boolean_face_info(...).synthetic_domain` is true (the seam-only test).
   If neither operand is a sphere, return `nullopt`.
2. If exactly one is a sphere, return `nullopt` (sphere against other surfaces is
   handled by the existing code).
3. If both are spheres: if either is not a whole sphere, return
   `Errc::unsupported "boolean_solid_curved: sphere/sphere booleans need
   whole-sphere operands"`.
4. `n = unit(cb - ca)`, `d = |cb - ca|`. If
   `d <= tolerance || d >= ra + rb - tolerance || d <= |ra - rb| + tolerance`,
   return `Errc::unsupported "boolean_solid_curved: sphere/sphere booleans need
   two crossing spheres"`.
5. If `|dot(axis_a, n) - 1| <= 8*eps` and `|dot(axis_b, -n) - 1| <= 8*eps`,
   return `nullopt` (already re-charted; this terminates the recursion).
6. Choose one shared frame: `x_axis = unit(cross(n, candidate))` with the same
   candidate rule as Task 1. Rebuild
   `a2 = make_sphere(ca, n, x_axis, ra, tolerance)` and
   `b2 = make_sphere(cb, -n, x_axis, rb, tolerance)`. Propagate errors. Return
   the two Breps with `solid_ids().front()`.

**Driver hook** — after `infos_b` is built (`:1704`) and before the torus gate:

```cpp
if (const auto rebased = detail::curved_boolean_rebase_sphere_pair(
        a, a_solid, b, b_solid, tolerance);
    !rebased.has_value()) {
  record_failure();
  return rebased.error();
} else if (rebased.value().has_value()) {
  const detail::CurvedBooleanRebasedSpheres& pair = *rebased.value();
  return boolean_solid_curved(pair.a, pair.a_solid, pair.b, pair.b_solid, op,
                              options, failure_report);
}
```

The recursive call re-enters with aligned charts, so step 5 returns `nullopt` and
the normal pipeline runs.

**Step 1 — failing test**

```cpp
TEST_CASE("curved boolean: sphere/sphere is not a splitter problem") {
  const auto a = geom::make_sphere({0.0, 0.0, 0.0}, 1.0, 1e-9);
  const auto b = geom::make_sphere({1.5, 0.0, 0.0}, 1.0, 1e-9);
  REQUIRE(a.has_value());
  REQUIRE(b.has_value());
  const auto result = geom::boolean_solid_curved(
      a.value(), a.value().solid_ids().front(), b.value(),
      b.value().solid_ids().front(), geom::SolidBooleanOp::intersection);
  INFO((result.has_value() ? std::string{} : std::string(result.error().message)));
  CHECK(result.has_value());
}
```

Today it fails with
`"split_face_curved: open cuts on a seam-only face"`.

**Step 2 — implement, run**

Focused test passes; suite green.

**Step 3 — negative tests**

A sphere face extracted from a box (or a two-face solid) returns the
whole-sphere message; `d = 2.0`, `0.5`, `2.5` return the crossing message and
leave the input Breps unchanged (compare a volume/edge-count fingerprint).

**Step 4 — commit**

`feat(m3c3b): re-chart whole-sphere operands for lens booleans`

## Task 4: Reversed spherical cap for difference

**Files**
- Modify `include/geom/brep/boolean_solid_curved.hpp` (pre-check `:1746-1797`,
  harvest `:2028-2069`).
- Modify `tests/unit/brep_boolean_curved_test.cpp`.

**Interfaces**
- Produces: no API change. A `Sphere3` fragment whose result interior is on the
  normal side is emitted with reversed rings instead of an error.

**Algorithm**

- In the pre-check, skip `Sphere3` faces when looking for the "flipped curved
  fragment" so the difference still reaches the pipeline; keep rejecting every
  other curved type.
- In the harvest, when `flip` is true and the surface is `Sphere3`:
  - call `curved_boolean_reverse_rings(assembled.loops, assembled.loops_uv,
    assembled.arc_spans)`;
  - keep `assembled.surface` as the sphere;
  - keep `assembled.reversed = false` (the reversed ring winding already carries
    the sign; `curved_solid_volume` reads the face winding).

Do not run the plane reprojection block for a sphere.

**Step 1 — failing test**

```cpp
TEST_CASE("curved boolean: difference subtracts a sphere cap") {
  // A - B with d = 1.5, r = 1: 39*pi/32
  ...
  REQUIRE(result.has_value());
  const double volume =
      geom::detail::curved_solid_volume(result.value().brep, result.value().solid, 1e-9)
          .value();
  CHECK(volume == doctest::Approx(39.0 * std::numbers::pi / 32.0).epsilon(1e-9));
}
```

Today it fails at the pre-check with the flipped-fragment message.

**Step 2 — implement, run**

Both `A - B` and `B - A` pass with volumes 3.8288160465625605; both are
validated, single-shell, and share one rim edge.

**Step 3 — commit**

`feat(m3c3b): emit reversed spherical caps for difference`

## Task 5: End-to-end corpus, QA-16, ADR and status

**Files**
- Modify `tests/unit/brep_boolean_curved_test.cpp`.
- New `docs/adr/0030-sphere-lens-booleans.md`.
- New `docs/benchmarks/m3c3b.md`.
- Modify `docs/architecture/geometry-library-architecture.md` §16.3.
- Modify header comments in `boolean_solid_curved.hpp` (the "deferred" list).

**Interfaces**
- Consumes the whole pipeline:
  `boolean_solid_curved`, `validate_topology`, `detail::curved_solid_volume`,
  `shell_is_closed_by_uses`.

**Test corpus (QA-16)**

| op | expected volume | tolerance |
|---|---|---|
| intersection | `11π/96 = 0.35997415822383048` | 1e-9 relative |
| union | `245π/96 = 8.0176062513489515` | 1e-6 relative |
| difference | `39π/32 = 3.8288160465625605` | 1e-6 relative |

For each: exactly 2 faces, 1 shell, `validate_topology(brep).ok()`, closed by
coedge uses, one shared rim edge (the two faces' loops use the same edge in
opposite directions), and two consecutive calls produce the same volume bitwise.

Generality cases: unequal radii (`ra = 1.2, rb = 0.8`, `d = 1.4`), and centres off
all axes (`ca = {0.3,-0.2,0.1}`, `cb = ca + 1.5*unit({1,2,2})`). Compute the
reference with the same closed-form expressions, not with the implementation.

**Step 1 — failing tests**

All three rows fail today: intersection/union fail in the splitter, difference at
the flipped pre-check.

**Step 2 — implement, run**

Debug and Release: QA-16 plus the full suite pass.

**Step 3 — docs**

- ADR-0030: whole-sphere precondition, lens-only crossing precondition, chart
  change, orientation handling, deferred cases.
- `docs/benchmarks/m3c3b.md`: machine, per-op runtime in Release (informational,
  no threshold), volumes, cuts/kept, residual.
- §16.3: M3c-3b status plus the split list for M3c-3c/3d.

**Step 4 — commit**

`test(m3c3b): pin the sphere lens corpus and record QA-16`

## M3c-3b Exit Criteria

- `union`, `intersection`, `difference` of two crossing whole spheres return
  validated single-shell solids with the volumes above.
- `A - B` and `B - A` differ and both match their closed forms.
- Whole-sphere, containment, tangency and separation inputs return documented
  `unsupported` errors and leave inputs unchanged.
- Debug and Release: full suite green, zero warnings under `GEOM_WERROR=ON`.
- ADR-0030, QA-16 and §16.3 merged; CI six jobs green after the merge.

## Self-Check

- [ ] Every task has Files, Interfaces, and the five steps.
- [ ] Each task's RED is a concrete failing assertion or compile error.
- [ ] Values in the corpus are closed forms, not implementation outputs.
- [ ] The plan never claims cone, torus, NURBS or general closed-cut support.
- [ ] No placeholder or self-correction note remains.
