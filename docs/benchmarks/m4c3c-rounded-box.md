# M4c-3c rounded-box fillet benchmark

- Date: 2026-11-28
- Commit: `sdd/m4c3c` (rounded-box code under test `0b854c1..07cb196`: Task 1
  `2cb37df`, Task 2 `c536f32` plus fix `f5d40c9`, Task 3 `78c3726`, controller
  mass-gate fix `7f9df10`, controller test pinning `07cb196`; the benchmark and
  docs are added by the Task 4 commit)
- Compiler: MSVC 19.44.35217, `/std:c++20 /O2`, `GEOM_WERROR=ON`
- Machine: Windows 10 Enterprise 22H2 (build 19045), Intel Core i7-8550U
  (4C/8T), Release
- Benchmark: `bench/bench_m4c3c.cpp`, five timed runs, min/median/max

## Workload

The `100 x 100 x 100` box at the origin is built with
`make_box(Point3d{0, 0, 0}, Point3d{100, 100, 100}, 1e-9)`. All twelve edge ids
from `edge_ids()` are passed to `fillet_edges` with `FilletOptions{1.0, 1e-9}`.
The timed region is one `fillet_edges` call; it builds a fresh result per
sample, and the same immutable input B-rep is reused for all five timed runs.

With `A = B = C = 100 - 2*1 = 98` the Steiner volume of the rounded box is

`V = ABC + 2(AB + BC + CA) r + pi (A + B + C) r^2 + 4 pi r^3 / 3`

which for `r = 1` evaluates to

`941192 + 57624 + 923.628240155 + 4.188790205 = 999743.8170303602`.

The expected topology is `V/E/F = 24/48/26` (six planes, twelve quarter
cylinders, eight spherical octants; twenty-four straight tangency lines and
twenty-four great-circle arcs) and `chi = 2`.

## Results

| Operation | min (ms) | median (ms) | max (ms) | Volume | Check |
| --- | ---: | ---: | ---: | ---: | --- |
| Fillet all 12 box edges, r=1 | 29.784 | 31.382 | 33.426 | 999743.817030359060 | within 1e-9 relative of 999743.8170303602 |

The table is the Task 4 evidence run, five timed iterations on an idle machine.
The measured volume differs from the analytic value by about `1.2e-15`
relative, well inside the required `1e-9` band. Every iteration checks, outside
the timed region: the result exists, `filleted_edges == 12`, exactly one solid
and one shell, `vertex_count`/`edge_count`/`face_count` are `24/48/26`,
`validate_topology` reports no issue, and the volume from
`detail::curved_solid_volume` is inside the `1e-9` relative band. Five repeat
batches on the same idle host measured `30.359 / 36.664 / 47.081`,
`30.631 / 42.220 / 53.931`, `30.407 / 30.703 / 38.831`,
`29.368 / 33.514 / 36.774` and `28.709 / 31.024 / 38.051` ms; the
tens-of-percent spread is host noise, not a code change, so the table is one
representative five-run batch rather than a best-of-N.

The benchmark does **not** tessellate. It therefore makes no watertightness
claim: watertightness of the eight spherical octants is the tessellation
boundary check exercised by the rounded-box unit tests (`boundary_edges == 0`),
not something this benchmark measures. Topology validation here means
`validate_topology` accepts the B-rep, not that a triangle mesh is watertight.

The binary exits 0 only when all five sample checks pass, and 2 with a clear
message on assembly or fillet failure, a wrong `filleted_edges` count, wrong
solid/shell or `V/E/F` counts, a non-clean `validate_topology`, a failed
`curved_solid_volume`, or a volume outside the `1e-9` relative band. It is
built only with `-DGEOM_BUILD_BENCH=ON` and is not part of ctest.

## Test evidence

- Debug `rounded box*,sphere patch*`: 11 cases / 1,769 assertions, 0 failed,
  1,173 skipped (`.\build-debug\tests\Debug\geom_tests.exe
  --test-case="rounded box*,sphere patch*"`, exit 0). This lane covers the
  rounded-box suite and the trimmed-sphere tessellation cases the octants rely
  on.
- Full Release (the gate): `.\build-rel\tests\Release\geom_tests.exe` —
  1,184 cases / 685,512 assertions, 0 failed, 0 skipped, exit 0.
- ctest Release: 1/1 test passed (`geom_tests`), 865.52 s, exit 0.
- ctest Debug: **not a gate and not completed**. The run was stopped by the
  controller after 3,014 s because it was already inside the known
  >20-minute boolean cases; it was not rerun. The Debug
  `rounded box*,sphere patch*` lane above is the Debug evidence.
- Full Debug was **not** run: the suite still contains the pre-existing
  >20-minute boolean cases documented in earlier milestones, so this lane is
  not a full Debug pass and is never reported as one.

## Scope

Accepted: an axis-aligned box; exactly its twelve live edges selected in any
order, duplicate ids collapsed; one finite radius with `0 < r < 0.5 *
min(a, b, c)`. The twenty-six patches are exact — six planar inset rectangles,
twelve quarter cylinders, eight spherical octants — trimmed by exact lines and
great-circle arcs; no Boolean, no SSI, no faceted or sampled substitute.

Deferred to M4c-3d: partial corner selections (any subset that is not all
twelve edges keeps the M4c-3a/3b behaviour), mixed edge families, variable
radius, concave or non-box operands, and the general rolling-ball blend. A
radius just below `0.5 * min` is accepted: the inset strip is geometrically
thin but topologically valid, and only equality is rejected.

## Reproduce

```powershell
cmake -S . -B build-bench -DGEOM_BUILD_BENCH=ON -DGEOM_BUILD_TESTS=OFF -DGEOM_WERROR=ON
cmake --build build-bench --config Release --target geom_bench_m4c3c
.\build-bench\Release\geom_bench_m4c3c.exe
```

CI runs the same binary as `geom_bench_m4c3c` in the bench job (the step uses
`set -o pipefail` so a failure is not masked by `tee`) and uploads
`bench-m4c3c.txt` with the other benchmark artifacts.
