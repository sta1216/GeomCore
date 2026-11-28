# M4d-1 profile extrusion benchmark

- Date: 2026-11-28
- Commit: `sdd/m4d1` (extrusion code under test `8df76cf..63e5cfe`: Task 1
  `ad11c24`, Task 2 `c223ef2` plus pinning `9ab7000`, Task 3 `c536cca`,
  certificate reuse `63e5cfe`; the benchmark and docs are added by the Task 4
  commit)
- Compiler: MSVC 19.44.35217, `/std:c++20 /O2`, `GEOM_WERROR=ON`
- Machine: Windows 10 Enterprise 22H2 (build 19045), Intel Core i7-8550U
  (4C/8T), Release
- Benchmark: `bench/bench_m4d1.cpp`, five timed runs, min/median/max

## Workload

Both profiles are swept along the `+z` normal by a height of `1` through
`extrude_profile(Plane3{Vec3d{0, 0, 1}, 0.0}, profile, 1.0, 1e-12)`. Each timed
sample builds a fresh result from the same immutable profile; the five samples
are timed, and the correctness checks run outside the timed region.

- **1024-gon.** A single counter-clockwise ring of 1024 vertices on radius `10`.
  The analytic area is `0.5 * 1024 * 10^2 * sin(2*pi/1024)`, so the volume is
  that area times height `1`:
  `314.157294036709118...`. Expected topology `V/E/F = 2048/3072/1026`
  (`V = 2N`, `E = 3N`, `F = N + 2`) and `chi = 2`.
- **512-gon with a 256-gon hole.** A counter-clockwise `512`-gon of radius `10`
  with a counter-clockwise `256`-gon hole of radius `5` at the origin; the
  profile area is the outer polygon area minus the hole area, and the expected
  volume is `235.619448841111847...`. Expected topology `V/E/F = 1536/2304/770`.
- **Certificate alone.** The second profile is also passed to
  `triangulate_polygon2` five times, to show what share of the extrusion is the
  simplicity certificate.

## Results

| Operation | min (ms) | median (ms) | max (ms) | Volume | Check |
| --- | ---: | ---: | ---: | ---: | --- |
| Extrude 1024-gon, height 1 | 286.401 | 297.212 | 322.401 | 314.157294036708493 | within 1e-9 relative of 314.157294036709118 |
| Extrude 512-gon + 256-gon hole, height 1 | 126.939 | 143.930 | 152.921 | 235.619448841109261 | within 1e-9 relative of 235.619448841111847 |
| Triangulate 512-gon + 256-gon hole (certificate) | 103.813 | 127.637 | 152.904 | 768 triangles | `ok == true` |

The table is the Task 4 evidence run, five timed iterations on an idle machine.
The measured volumes differ from the analytic values by about `2e-15` and
`1.1e-14` relative, well inside the required `1e-9` band. Every iteration checks,
outside the timed region: the result exists, `vertex_count`/`edge_count`/
`face_count` are `2048/3072/1026` and `1536/2304/770`, `validate_topology`
reports no issue, `mass_properties` gives the volume above, and `tessellate`
returns a non-empty triangle list.

The holed extrusion costs `126.9–152.9 ms` while the certificate alone costs
`103.8–152.9 ms`; the difference is the wall construction, the assembler and the
volume integration, and the caps are not triangulated a second time. That is the
point of the reuse: before `63e5cfe` the same extrusion paid three
triangulations per call, measured at `439226 us` per call on a 20-iteration
probe versus `176075 us` after (`2.49x`, Task 3 report).

Host noise is large in the tens of percent: repeat batches on the same idle host
spread from `286` to `322 ms` for the 1024-gon. The table is one representative
five-run batch rather than a best-of-N.

The binary exits 0 only when all five sample checks pass for each workload, and
2 with a clear message on an extrusion failure, wrong `V/E/F` counts, a
non-clean `validate_topology`, a failed `mass_properties` or a volume outside the
`1e-9` relative band. It is built only with `-DGEOM_BUILD_BENCH=ON` and is not
part of ctest.

## Test evidence

- Debug `*extrude*`: 14 cases / 99 assertions, 0 failed, exit 0
  (`.\build\tests\Debug\geom_tests.exe --test-case="*extrude*"`).
- Debug filtered extrusion lane `*extrude*,assemble:*,*triangulate2*`:
  49 cases / 10461 assertions, 0 failed, exit 0. Release rerun of the same lane:
  49 / 10461, 0 failed.
- Full Release (the gate): `.\build\tests\Release\geom_tests.exe` —
  1,208 cases / 685,744 assertions, 0 failed, exit 0, run on the tree of the
  fallback-pin commit (Task 4 and the pin commit add only the benchmark, tests
  and docs, no library code).
- Full Debug was **not** run and is never reported: the suite still contains the
  pre-existing >20-minute boolean cases documented in earlier milestones.

## Scope

Accepted: one planar profile with an optional set of non-touching holes wound in
any orientation, one positive finite height, one straight extrusion along the
plane normal, one solid, exact planar faces and exact volume. Rejected before any
face is built: non-finite inputs, a malformed layout, an undefined ring
orientation, repeated vertices, non-simple profiles (self-intersections, holes
outside the outer ring, nested holes) and assembly failures — all with the
stable codes and messages of the ADR-0045 table.

Deferred to M4d-2: oblique or sheared extrusion (a direction vector rather than
the plane normal), revolve, loft and sweep, variable height or taper, multiple
solids, and any repair or healing of a self-intersecting profile. The
simplicity certificate is an operational predicate, and ear clipping is O(n²),
which is why the 1024-gon costs about `300 ms`.

## Reproduce

```powershell
cmake -S . -B build-bench -DGEOM_BUILD_BENCH=ON -DGEOM_BUILD_TESTS=OFF -DGEOM_WERROR=ON
cmake --build build-bench --config Release --target geom_bench_m4d1
.\build-bench\Release\geom_bench_m4d1.exe
```

CI runs the same binary as `geom_bench_m4d1` in the bench job (the step uses
`set -o pipefail` so a failure is not masked by `tee`) and uploads
`bench-m4d1.txt` with the other benchmark artifacts.
