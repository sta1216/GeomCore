# M4c-3b convex prism fillet benchmark

- Date: 2026-11-27
- Commit: `sdd/m4c3b` (prism code under test `2f5a258..b77cd88`; benchmark and
  docs added in the Task 4 commit)
- Compiler: MSVC 19.44.35217, `/std:c++20 /O2`, `GEOM_WERROR=ON`
- Machine: Windows 10 Enterprise 22H2 (build 19045), Intel Core i7-8550U
  (4C/8T), Release
- Benchmark: `bench/bench_m4c3b.cpp`, five timed runs, min/median/max

## Workload

A regular 64-gon of radius 2 in `z = [0, 10]` is built with the public
`assemble_solid` fixture used by
`tests/unit/fillet_prism_test.cpp::make_prism`: CCW top ring, reversed bottom
ring, and outward-facing quad sides. All 64 vertical edges (equal endpoint
`x`/`y`, differing `z`) are selected and passed to `fillet_edges` with
`FilletOptions{0.05, 1e-9}`. The same immutable input is reused for all five
timed runs.

For a regular `n`-gon the interior angle is `theta = pi - 2*pi/n`; the removed
cross-section area per corner is
`r^2 * (cot(theta/2) - (pi - theta)/2)`. With `n = 64`, `r = 0.05`,
`h = 10`, the analytic volume is

`125.461876478546088`.

The expected topology is `V/E/F = 256/384/130` (`2n + 2k`, `3n + 3k`,
`n + 2 + k` with `n = k = 64`) and `chi = 2`.

## Results

| Operation | min (ms) | median (ms) | max (ms) | Volume | Check |
| --- | ---: | ---: | ---: | ---: | --- |
| Fillet 64 vertical edges, r=0.05 | 11.683 | 14.477 | 15.385 | 125.461876478546259 | within 1e-9 relative of 125.461876478546088 |

The table is the Task 4 evidence run, five timed iterations on an idle machine.
Every iteration checks, outside the timed region: the result exists,
`filleted_edges == 64`, exactly one solid and one shell, `vertex_count`,
`edge_count` and `face_count` are `256/384/130`, `validate_topology` reports no
issue, and the volume from `detail::curved_solid_volume` is inside the `1e-9`
relative band. The timed region covers `fillet_edges` only; checks and the
volume integration run after each sample is taken. Two repeat batches on the
same idle host measured `13.094 / 14.042 / 18.298 ms` and
`12.298 / 15.362 / 16.400 ms`; the few-millisecond spread is host noise, not a
code change, so the table is one representative five-run batch rather than a
best-of-N.

The benchmark does **not** tessellate. It therefore makes no watertightness
claim: watertightness is the tessellation boundary check exercised by the prism
unit tests (`fillet_edges rounds a prism in a tilted frame` and
`fillet_edges leaves its input unchanged and is deterministic`), not something
this benchmark measures. Topology validation here means
`validate_topology` accepts the B-rep, not that a triangle mesh is watertight.

The binary exits 0 only when all five sample checks pass, and 2 on assembly or
fillet failure, a wrong `filleted_edges` count, wrong solid/shell or `V/E/F`
counts, a non-clean `validate_topology`, a failed `curved_solid_volume`, or a
volume outside the `1e-9` relative band. It is built only with
`-DGEOM_BUILD_BENCH=ON` and is not part of ctest.

## Test evidence

- `*prism*`: 11 cases / 98 assertions in Debug and 11 / 98 in Release, 0
  failed. This filter also matches the pre-existing
  `chamfer one box edge removes a triangular prism` case and misses three prism
  fillet cases whose names do not contain "prism"; the prism file itself adds
  13 cases / 127 assertions.
- `*fillet*`: 21 cases / 1,014 assertions in Debug and 21 / 1,014 in Release,
  0 failed (8 box cases / 887 assertions plus the 13 prism cases / 127
  assertions).
- Full Release: 1,172 cases / 683,738 assertions, 0 failed, 956 s.
- Full Debug remains PENDING: the suite still contains the known >20-minute
  cases documented in earlier milestones, so the Debug `*prism*` lane is the
  prism evidence, not a full Debug pass.

## Scope

Accepted: a single-shell, closed, convex, all-planar prism whose filletable
edges are parallel to the prism axis; any non-empty subset of those vertical
edges; one uniform radius. At each selected corner the radius must satisfy
`r * cot(theta/2) < both incident side lengths`, and every side must keep a
positive strip after both trims. Each patch is an exact cylinder trimmed by
cap-plane arcs; no Boolean and no faceted approximation.

Deferred to M4c-3c: adjacent non-vertical edge blends, spherical vertex
patches, variable radius, non-convex polygons, multi-shell solids, curved side
surfaces and general operands.

## Reproduce

```powershell
cmake -S . -B build-bench -DGEOM_BUILD_BENCH=ON -DGEOM_BUILD_TESTS=OFF -DGEOM_WERROR=ON
cmake --build build-bench --config Release --target geom_bench_m4c3b
.\build-bench\Release\geom_bench_m4c3b.exe
```

CI runs the same binary as `geom_bench_m4c3b` in the bench job (the step uses
`set -o pipefail` so a failure is not masked by `tee`) and uploads
`bench-m4c3b.txt` with the other benchmark artifacts.
