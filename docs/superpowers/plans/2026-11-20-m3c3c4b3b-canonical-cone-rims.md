# M3c-3c-4b3b: Canonical Rim Sampling for Split Cones

## Goal

Make `split_solid_by_plane` of a cone produce meshes that are watertight on every
platform, and restore the cone split rows to the tessellation matrix.

The split Brep is topologically correct: a probe of `cone z=1 below` shows three
faces, three edges and two vertices, with the bottom rim a single closed
`CircleCurve3` edge used by both the cap (`rev=0`) and the lateral band (`rev=1`),
`t=[0, 2π]` on both. The mesh nevertheless breaks on glibc: 214 positions / 304
triangles / 122 boundary edges versus the pinned 134 / 264 / 0 on MSVC. Because both
faces share one edge and one parameter range, the defect is in the rim vertex
*positions* produced by the two independent meshing paths, not in the topology.

## Architecture

`tessellate()` meshes the two cap faces and the trimmed lateral band through
different code paths. Each path places the shared rim vertices by its own
evaluation (parametric surface evaluation versus circle evaluation), so the two
rims differ in the last bit on glibc and the exact-position closure test fails.

The fix is to canonicalize every full-turn circular rim run: keep the run's
vertex count, but replace each position with `circle.eval(t0 + i·2π/N)` where `N`
is the run's segment count. Both faces then emit bitwise-identical rim positions
on every libm, with no change to triangle counts, ring order or topology.

## Tech Stack

C++20, MSVC 19.44 (primary), GCC 15.2 / Clang on CI, doctest 2.4.11, CMake.
Predicate-level decisions stay exact; this change only affects sampling positions.

## Spec

- `docs/architecture/geometry-library-architecture.md` §8.10 (meshing contract:
  watertight, no flipped triangles, deterministic) and §16.3.
- `docs/quality/testing-strategy.md` (invariants, differential tests).
- ADR-0036 (cone tessellation metric) and the M3c-3c-4b3a test matrix.

## Global Constraints

- Do not change `split_solid_by_plane` topology: the shared closed edge and its
  `[0, 2π]` range are already correct.
- Non-cone and non-circular output must remain bitwise identical; the existing
  hash pins (QA-13, QA-15, canonical primitives, trimmed cylinder) must not move.
- Full-turn snapping must be deterministic and independent of the ring's starting
  index: rotating a synthetic ring then snapping must give the same sequence.
- No tolerance-based welding of finished meshes; closure is still verified by exact
  position equality.
- TDD: each task starts with a failing test. The libm-specific failure cannot be
  reproduced on MSVC, so the unit test forces a phase shift (rotate the start by
  one segment) to prove start independence; the end-to-end guarantee is verified by
  the ubuntu CI job.
- After header edits run `/t:Rebuild` or touch a `.cpp`; MSBuild header tracking is
  unreliable in this repository.

## Review Focus

1. Snapping formula: `N` is the run's segment count; positions are
   `circle.eval(t0 + i·2π/N)`, `i = 0..N-1`; no angle accumulation.
2. Start independence: rotating the input ring by one segment yields the same
   snapped sequence.
3. Call sites: every full-turn circular run of a split cone face (cap and band),
   and no other surface.
4. Non-circular runs, partial arcs, seam columns and poles are untouched.
5. Existing hash pins unchanged; `split_solid` topology assertions unchanged.
6. Tests fail when either face bypasses the snapping.

## File Map

- `include/geom/brep/tessellate.hpp` — snapping helper and its use in the cap and
  trimmed-band ring construction.
- `include/geom/brep/split_solid.hpp` — only if a rim run needs its circle edge
  exposed; prefer no change.
- `tests/unit/brep_boolean_curved_test.cpp` — restore the cone split rows and add
  rim-identity assertions.
- `tests/unit/tessellate_trimmed_test.cpp` — helper-level start-independence test.
- `docs/adr/0036-cone-tessellation.md`, `docs/benchmarks/m3c3c4b3a.md`,
  architecture §16.3 — remove the "split deferred" limitation once green.

## Task 1: Full-turn rim snapping helper

**Files**

- `include/geom/brep/tessellate.hpp` (helper only)
- `tests/unit/tessellate_trimmed_test.cpp`

**Interfaces**

Produces:

```cpp
namespace geom::v1::detail {

// Replaces `segments` consecutive positions of a closed ring with the uniform
// samples of `circle`, starting at `t0`: position i becomes
// circle.eval(t0 + i * 2π / segments). The run covers a full turn, so index 0
// and index segments coincide; `points[start]` stays the shared seam vertex.
void tessellate_snap_full_circle_run(std::vector<Point3d>& points,
                                     std::size_t start, std::size_t segments,
                                     const CircleCurve3& circle, double t0);
}
```

The helper must not accumulate angles and must not depend on traversal order.

**Steps**

1. Add the test `tessellate: full-turn rim snapping is start independent`: build a
   ring of 41 points on radius 0.5 with arbitrary start angle, snap twice (once at
   the ring's own start, once rotated by one segment), and require the two snapped
   sequences to be identical after rotation. Require each snapped point to equal
   `circle.eval(t0 + i*2π/40)` bitwise.
2. Run the test; it fails to compile (the helper does not exist).
3. Implement the helper: validate `segments >= 3`, `circle.is_valid()`, finite
   `t0`; write `points[start + i] = circle.eval(t0 + 2π*i/segments)` for
   `i = 0..segments-1`.
4. Run the test, then the full suite. Expected counts unchanged (one new case and
   its assertions).
5. Commit: `feat(brep): add full-turn circle rim snapping`.

## Task 2: Use snapped rims for split cone meshes

**Files**

- `include/geom/brep/tessellate.hpp`
- `tests/unit/brep_boolean_curved_test.cpp`

**Interfaces**

Consumes: `tessellate_snap_full_circle_run`.

Produces: snapping applied in the ring construction used by both the planar cap
and the trimmed cone band, for runs whose points all lie on one closed full-turn
circle edge.

**Steps**

1. Add the split row back to `curved boolean: cone tessellation follows the
   surface` (z=1 below and above): counts 134/264, boundary 0, `closed()`, volume
   `1.830287422577` (below). Also require the cap rim positions and the band rim
   positions to be identical (sorted arrays, bitwise). On MSVC this passes; on
   ubuntu the assertion must turn green only after the fix (RED = CI run
   44715035487).
2. Run the focused test; record the MSVC pass and the CI reference.
3. Locate the boundary ring construction for the planar cap and for the trimmed
   cone band. When a run of ring points lies on the face's closed circle edge with
   `t1 - t0 == 2π` and its segment count matches the edge's canonical count, snap
   that run with the helper (keep the ring's seam vertex).
4. Rebuild and run the focused test, QA-22 and the trimmed tessellation suite.
5. Commit: `fix(brep): snap split cone rims onto one circle sampling`.

## Task 3: Restore the matrix and document

**Files**

- `tests/unit/brep_boolean_curved_test.cpp`
- `docs/adr/0036-cone-tessellation.md`
- `docs/benchmarks/m3c3c4b3a.md`
- `docs/architecture/geometry-library-architecture.md`
- `docs/superpowers/plans/2026-11-20-m3c3c4b3a-cone-tessellation.md`

**Steps**

1. Restore the full split matrix (z=0.5/1/1.5, translated split) with the
   canonical counts and volumes; update QA-22 accordingly.
2. Update ADR-0036, the benchmark and §16.3: remove the "split-solid cone meshes
   are deferred" sentence and state that all meshed faces of one split cone share
   canonical full-turn rim sampling.
3. Run Debug and Release full suites; confirm every other hash pin is unchanged.
4. Commit: `docs(brep): record canonical split-cone rims`.

## Exit Criteria

- The split matrix passes on windows-2022, ubuntu-24.04 and macos-15.
- Split cone meshes are closed with 0 boundary edges on every platform; cap and
  band rim positions are bitwise identical.
- Non-cone output hashes unchanged; QA-22 volumes within the existing bounds.
- Debug and Release full suites green with 0 warnings under `GEOM_WERROR=ON`.
- ADR-0036, benchmark, §16.3 and the plan reflect the final contract.

## Self-Review

- [ ] Every task has Files, Interfaces and a five-step TDD sequence.
- [ ] The helper test can fail on MSVC by rotating the start; the end-to-end RED is
      the ubuntu CI run.
- [ ] No placeholder or stale count remains.
- [ ] `git diff --check` clean under the repository's `core.whitespace` setting.
