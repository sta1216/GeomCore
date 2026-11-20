# M3c-3c-4b3b: Canonical Rims for Split Cones

## Goal

Make `split_solid_by_plane` of a cone produce meshes that are watertight on every
platform, and restore the cone split rows to the tessellation matrix.

The split Brep is topologically correct: a probe of `cone z=1 below` shows three
faces, three edges and two vertices, with the bottom rim a single closed
`CircleCurve3` edge used by both the cap and the lateral band, `t=[0, 2pi]` on
both. The mesh nevertheless broke on glibc: 214 positions / 304 triangles / 122
boundary edges versus the pinned 134 / 264 / 0 on MSVC. Since both faces share one
edge and one parameter range, the defect was in the rim vertex *positions* produced
by the two independent meshing paths, not in the topology.

## Architecture

`tessellate()` meshes cap faces and the trimmed lateral band through different code
paths. Each path placed the shared rim vertices by its own evaluation (circle
evaluation versus parametric surface evaluation), so the two copies differed in the
last bit on glibc and exact-position closure failed.

The fix is a post-pass in `tessellate()` for solids that carry a `Cone3`:

1. For every `CircleCurve3` edge with `v0 == v1` and a full turn (normalize a
   reversed interval first, then `|(t1-t0) - 2pi| <= 1e-12`), compute the run's
   segment count `N = tessellate_circle_segments(circle, t0, t1, deflection, true)`
   and `step = (t1-t0)/N`.
2. Select a mesh position only when it lies in the circle's plane and on its
   radius: `|dot(rel, normal)| <= 1e-9 * max(1, r)` and
   `| |rel - axial*normal| - r | <= 1e-9 * max(1, r)`. Snap that position to the
   canonical sample `turns = round((atan2(dot(in_plane,y), dot(in_plane,x)) - t0)/step)`,
   wrapped modulo `N`, position `circle.eval(t0 + step*turns)`. The plane/radial
   split matters: a 45-degree cone puts its apex on the base-rim sphere.
3. Merge positions that are bitwise identical (exact `double` key), remap the
   triangles, and drop the zero-area triangles that merging produced (keeping
   `triangle_face` in sync).

The snap is confined to cone-carrying solids, so every non-cone mesh stays bitwise
identical. The merge is exact; no tolerance-based welding is introduced.

## Tech Stack

C++20, MSVC 19.44 (primary), GCC 15.2 / Clang on CI, doctest 2.4.11, CMake.

## Spec

- `docs/architecture/geometry-library-architecture.md` section 8.10 (meshing
  contract) and section 16.3.
- `docs/quality/testing-strategy.md`.
- ADR-0036 and the M3c-3c-4b3a matrix.

## Global Constraints

- Do not change `split_solid_by_plane` topology or the shared `[0, 2pi]` edge.
- Non-cone output must remain bitwise identical; QA-13, QA-15 and the canonical
  primitive hashes must not move.
- The snap must be deterministic and depend only on the circle, its range and the
  deflection, never on ring order or a previous run.
- Closure is still verified by exact position equality; the merge key is exact.
- TDD: the libm-specific failure is not reproducible on MSVC, so the new test pins
  the Linux numbers as well; the ubuntu CI job is the cross-platform gate.
- After header edits run `/t:Rebuild` or touch a `.cpp`; MSBuild header tracking is
  unreliable here.

## Review Focus

1. Snap formula and wrap (`round`, modulo `N`), including the `t0` phase.
2. Full-turn detection (`v0 == v1`, `2pi` within `1e-12`) and cone-only gating.
3. Merge correctness: remapped triangles, `triangle_face` length, no zero-area
   triangles left, boundary/once counts.
4. Non-cone output bitwise unchanged (hash pins).
5. Split rows close: positions/triangles/volume/closed/used-once/nonmanifold.
6. Docs: ADR-0036, benchmark and architecture section 16.3 no longer claim the
   limitation.

## File Map

- `include/geom/brep/tessellate.hpp` - snap/merge post-pass.
- `tests/unit/brep_boolean_curved_test.cpp` - five split rows and their invariants.
- `docs/adr/0036-cone-tessellation-sag.md`, `docs/benchmarks/m3c3c4b3a.md`,
  architecture section 16.3, this plan - remove the limitation.

## Task 1: Snap and merge cone rims

**Files**: `include/geom/brep/tessellate.hpp`, `tests/unit/brep_boolean_curved_test.cpp`

**Interfaces**: Produces
`detail::tessellate_snap_cone_rims(TriMesh& mesh, const Brep& brep, double deflection)`.
Consumes `tessellate_circle_segments`, `CircleCurve3::eval`, `TriMesh`.

**Steps**

1. Add a failing test `tessellate: split cone fragments stay watertight` with the
   five rows (below, above, z=0.5, z=1.5, translated), asserting positions,
   triangles, `closed`, `used_on_boundary == 0`, `nonmanifold == 0`, volume within
   `5e-3` and a second run bitwise identical.
2. Run it; on MSVC it already passes after the product change, so also record the
   glibc probe numbers in the report.
3. Implement the post-pass as described; call it at the end of `tessellate()`.
4. Debug and Release suites, zero warnings, commit.

**Test values**

| row | positions | triangles | volume |
| --- | --- | --- | --- |
| z=1 below | 134 | 264 | 1.830287422577 |
| z=1 above | 51 | 98 | 0.261110903259 |
| z=0.5 below | 140 | 276 | 1.209417839998 |
| z=1.5 below | 124 | 244 | 2.058752793525 |
| translated split | 134 | 264 | n/a |

## Task 2: Documentation

**Files**: ADR-0036, `docs/benchmarks/m3c3c4b3a.md`, architecture section 16.3.

Replace the "split deferred" text with the snap/merge contract, restore the split
rows to the matrix, and update the assertion counts to 1014 / 673,425.

## Exit criteria

- Both configurations: 1014 cases / 673,425 assertions, zero warnings.
- The five split rows close with `used_once == 0` and plausible volumes.
- Non-cone hashes unchanged.
- CI ubuntu and coverage jobs green (the original failure).
- ADR, benchmark and architecture describe the implemented snap, not the deferred
  limitation.

## Self-review

- [x] The task has a failing test and the test values come from the probe.
- [x] The post-pass is gated to cone-carrying solids.
- [x] Exact merge, no new epsilon.
- [x] Docs and counts match the implementation.
- [x] Plan file map matches the changed files.
