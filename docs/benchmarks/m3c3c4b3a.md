# M3c-3c-4b3a: cone tessellation

Date: 2026-11-20. The matrix was measured at `75915ae` (Task 2 test commit `0a2556a`; Task 1 `3d0801a`). The full-suite counts were measured at `43218b4` after the split-solid rows left the matrix; the binaries are the Debug and Release builds of that commit.

## Scope

This note records the measured behavior of the cone chordal-sag change:

- trimmed cone frusta produced by `boolean_solid_curved` now tessellate watertight
  at `TessellationOptions{1e-3, 0.35}`; split-solid cone fragments are not claimed
  (see Limits);
- the split predicate and the final acceptance check measure the distance to the
  cone surface at the same azimuth and axial height, not the parameter-midpoint
  surface point (ADR-0036);
- apex-bearing cone fragments keep the canonical single-apex fan;
- every non-cone surface keeps its previous parameter metric and its previous
  counts.

The matrix is pinned by `tests/unit/brep_boolean_curved_test.cpp`:
`tessellate: trimmed cone frusta are watertight` (line 652, three QA-22 rows),
`tessellate: translated cone booleans stay watertight` (one translated boolean)
`tessellate: canonical cone faces keep their apex` (three apex rows) and `tessellate: the cone sag helper is cone-only` (line 785). The
shared census is `check_cone_tessellation` (line 439).

## Environment

- Windows 10 Enterprise 10.0.19045, build 19045, x64.
- Intel Core i7-8550U, 4 cores / 8 threads, 15.9 GB RAM.
- MSVC 19.44 (Visual Studio 17 2022, toolset 14.44.35207), CMake 4.4.4.
- Debug and Release, `GEOM_WERROR=ON`, `/W4 /WX /utf-8`.
- Single run per configuration on the host above; there is no timing gate.

## Fixture

The cone is `make_cone({0,0,0}, {0,0,1}, 1, 2, 1e-9)`: base radius `1`, height `2`,
so the radius at height `z` is `1 - z/2`. Boolean tolerance is `1e-9` and every
result is tessellated with `TessellationOptions{1e-3, 0.35}`. The `z = 1` split
leaves the upper small cone (radius `0.5`, height `1`, apex at `z = 2`) and the lower
frustum. The QA-22 box is `{-2,-2,1}..{2,2,4}`; the union uses the low box
`{-2,-2,1}..{2,2,3}`. The translated cone sits at `(1,2,3)`; its translated difference box is
`{-1,0,4}..{3,4,7}` at `z = 4`.

`V` and `T` are positions and triangles. Volume error is `mesh volume - analytic`;
the relative error is `|error| / analytic`.

## Matrix

| case | V | T | mesh volume | analytic | error | relative |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| QA-22 intersection | 51 | 98 | 0.261110903259 | `pi/12` = 0.261799387799 | -6.885e-4 | 2.630e-3 |
| QA-22 difference | 134 | 264 | 1.830287422577 | `7*pi/12` = 1.832595714594 | -2.308e-3 | 1.260e-3 |
| QA-22 union (low box) | 142 | 280 | 33.830287422577 | `32 + 7*pi/12` = 33.832595714594 | -2.308e-3 | 6.823e-5 |
| translated difference | 134 | 264 | 1.830287422577 | `7*pi/12` | -2.308e-3 | 1.260e-3 |
| apex intersection `z = 0.5` | 62 | 120 | 0.882011366069 | `0.28125*pi` = 0.883572933822 | -1.562e-3 | 1.767e-3 |
| apex intersection `z = 1.5` | 37 | 70 | 0.032559033313 | `0.0104166667*pi` = 0.032724923580 | -1.659e-4 | 5.069e-3 |
| whole cone | 72 | 140 | 2.091662475470 | `2*pi/3` = 2.094395102393 | -2.733e-3 | 1.305e-3 |

The maximum absolute volume error is `2.733e-3` at the whole cone. The maximum
relative error is `5.069e-3` at apex intersection `z = 1.5`, where the volume is only
`0.0327`; that row is still pinned to `1e-9` relative by its 12-decimal literal. The
original brief and plan text call `2.917e-3` the "maximum relative error"; the probe column is
absolute and the plan was corrected to "maximum absolute", so this note keeps both
numbers and labels the units. Task 1's
three-row QA-22 census uses a `5e-3` relative band; the other four rows use `1e-9`
relative.

## Watertightness and determinism

Every one of the seven rows reports:

- `closure.closed()` true, `used_once == 0` (no boundary edges) and
  `nonmanifold == 0`;
- `triangle_face.size() == triangles.size()`, with every face id below
  `face_count()`;
- a second `tessellate` call with bitwise-identical positions and identical
  `triangles` and `triangle_face` arrays;
- the volume above (within `5e-3` relative for the three QA-22 rows, `1e-9` relative
  for the other four).

The `z = 0.5` apex case has no vertex at `(0,0,0.5)`: its cap is an ear-clipped
polygon, so the cut plane is represented by the polygon boundary rather than a ring
of samples. That is the measured geometry, not a missing cap: the mesh is closed and
its volume is `0.882011366069` against `0.28125*pi`.

## Apex and unchanged surfaces

The apex-bearing rows have exactly one apex vertex and it is the cone's analytic apex:
`(0,0,2)` for the origin cases and `(1,2,5)` for the translated upper cone. The four
apex rows are the QA-22 intersection, the two apex intersections and the whole
cone. Every lower frustum and every difference/union has
zero apex vertices, as expected.

The canonical whole cone is 72/140 and the QA-22 intersection is 51/98; both keep
`tessellate_cone_face`'s single-apex fan. The QA-22 difference is 134/264.

Non-cone behavior is unchanged. QA-13's cut cylinder is still 894/449 at
`{2e-3, 0.35}`, the canonical box/cylinder/torus are still 12/8, 88/46 and 1152/576,
and `tessellate_cone_radial_sag_square` returns `nullopt` for a real plane face and a
real cylinder face (test line 801). A translated cone whose seam copies evaluate to
the same doubles is still the documented refusal, not an internal error.

## Focused and full runs

- Four-case matrix filter: 4 cases / 1,352 assertions, 0 failed.
- All cone tests with `*cone*`: 45 cases / 3,482 assertions, 0 failed.
- Full suites, both configurations: 1,013 cases / 671,966 assertions, 0 failed, zero
  warnings with `GEOM_WERROR=ON`. Debug and Release report the same numbers. Task 1's
  baseline was 1,010 / 671,694; Task 2 adds 3 cases and 272 assertions after the
  split-solid rows leave the matrix.

The timings are single-run wall clock and are recorded for scale only; nothing in the
test suite asserts them.

## Limits

- **The metric is a conservative upper bound, not the exact nearest distance.** The
  helper projects the point to the cone's azimuth and axial height and measures the
  distance to the surface point there. That point is not always the Euclidean nearest
  point, so the reported sag can exceed the true nearest-surface distance, and the
  refinement can split slightly more than a nearest-point metric would. It never
  under-reports: a mesh accepted by this check has true distance at most the reported
  value.
- **Apex faces are canonical only.** This slice does not rewrite the apex fan, the
  seam or the uv ring, and it does not insert an apex into a trimmed uv ring. An
  apex-bearing fragment is tessellated by `tessellate_cone_face`'s whole-ring fan;
  an apex ring that reached the trimmed path is not rectilinear and is rejected with
  `"tessellate: trimmed curved face boundary is not rectilinear"`.
- **The trimmed refiner is still bounded.** Each face has a `8*T + 64` split budget
  (`tessellate.hpp:1543`) and the trimmed mesh has a global cap of 2,000,000
  triangles (`tessellate.hpp:955`), so a pathological cone fragment still terminates
  with a refusal instead of running without bound.
- **No new tolerance.** The only tolerances are still `0.9 * deflection` for
  refinement and `deflection` for final acceptance. The sag helper itself is an exact
  double computation.
- **Split-solid cone fragments are not claimed.** `split_solid_by_plane` gives the
  two sides of the cut independent rim samples; MSVC happens to produce matching
  counts and the mesh closes, glibc does not (214 positions / 304 triangles / 122
  boundary edges against 134/264/0). The split-solid rows were removed from this
  matrix and the fix (shared edge samples) is M3c-3c-4b3b.
- **Not in this slice:** oblique or axis-containing cone planes, general
  curved/curved cone booleans, and re-centring translated cones whose seam copies
  collide. Those, together with the apex fan/seam/uv-ring rewrite, are M3c-3c-4b3
  work.
