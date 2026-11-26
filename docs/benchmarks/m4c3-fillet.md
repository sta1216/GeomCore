# M4c-3a box fillet benchmark

- Date: 2026-11-26
- Commit: `sdd/m4c3` (Task 4 tree; code under test `0fbe889..624bc44`; benchmark and docs amended in the Task 4 commit)
- Compiler: MSVC 19.44.35217, `/std:c++20 /O2`, `GEOM_WERROR=ON`
- Machine: Windows 10 Enterprise 22H2 (build 19045), Intel Core i7-8550U (4C/8T), Release
- Benchmark: `bench/bench_m4c3.cpp`, five timed runs, min/median/max

## Workload

`make_box({0,0,0}, {100,100,100}, 1e-9)`. The four vertical edges of the `+Z`
family (one at each corner of the `XY` cross-section) are selected by equal
endpoint `x`/`y` and differing `z`, and passed to `fillet_edges` with
`FilletOptions{1.0, 1e-9}`.

Analytic volume `100^3 - 4 * (1 - pi/4) * 1^2 * 100 = 999,914.159265359`; the
result has `V/E/F = 16/24/10` and `chi = 2`.

## Results

| Operation | min (ms) | median (ms) | max (ms) | Volume | Check |
| --- | ---: | ---: | ---: | ---: | --- |
| Fillet four vertical edges, r=1 | 7.587 | 8.247 | 12.766 | 999914.159265358 | within 1e-6 relative of 999914.159265359 |

The table is the Task 4 evidence run, five timed iterations on an idle machine.
Every iteration checks, outside the timed region: the result exists,
`filleted_edges == 4`, the result has at least one solid, `vertex_count`,
`edge_count` and `face_count` are `16/24/10`, `validate_topology` reports no
issue, and the volume is inside the `1e-6` relative band. The timed region
covers `fillet_edges` only; the checks and the volume integration
(`detail::curved_solid_volume`) run after each sample is taken. The benchmark
does not tessellate, so watertightness stays a unit-test claim, not a benchmark
check.

The binary exits 0 only when all five sample checks pass, and 2 on a missing
result, a `filleted_edges` count other than four, a missing solid, wrong
topology counts, a non-clean `validate_topology`, a failed
`curved_solid_volume`, or a volume outside the `1e-6` relative band. It is
built only with `-DGEOM_BUILD_BENCH=ON` and is not part of ctest.

## Test evidence

- `*fillet*`: 8 cases / 884 assertions in Debug and 8 / 884 in Release, 0
  failed (post-review Task 3 lane re-run on this tree).
- Debug lane `*fillet*,*blend*,*editor*`: 56 cases / 1,394 assertions, 0 failed.
- Full Release: 1,159 cases / 683,608 assertions, 0 failed (880 s).
- Full Debug remains PENDING: the suite still contains the known >20-minute
  cases documented in earlier milestones, so the Debug lane above is the fillet
  evidence, not a full Debug pass.

## Scope

Accepted: an axis-aligned box (single shell, closed, convex, all-planar) and
any non-empty subset of one parallel edge family (1–4 edges). The patch is an
exact quarter cylinder trimmed by circular arcs, `r < 0.5 * min(cross-section
extents)`, deterministic, and watertight when tessellated.

Deferred to M4c-3b: vertex blends, variable radius, concave edges, general
convex and curved operands, chamfer angles, offset and sweep. A faceted fillet
is explicitly not shipped as a substitute.

## Reproduce

```powershell
cmake -S . -B build-bench -DGEOM_BUILD_BENCH=ON -DGEOM_BUILD_TESTS=OFF -DGEOM_WERROR=ON
cmake --build build-bench --config Release --target geom_bench_m4c3
.\build-bench\Release\geom_bench_m4c3.exe
```

CI runs the same binary as `geom_bench_m4c3` in the bench job (the step uses
`set -o pipefail` so a failure is not masked by `tee`) and uploads
`bench-m4c3.txt` with the other benchmark artifacts.
