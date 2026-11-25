# M4c-1 planar chamfer and shell benchmarks

- Date: 2026-11-25
- Commit: `sdd/m4c1` (Task 7 evidence; see the plan Results for the merged hash)
- Compiler: MSVC 19.44.35217, `/std:c++20 /O2`, `GEOM_WERROR=ON`
- Machine: Windows 11, Intel Core i7-8550U (4C/8T), Release
- Benchmark: `bench/bench_m4c1.cpp`, five timed runs per operation, min/median/max

## Workload

`make_box({0,0,0}, {100,100,100}, 1e-9)`.

- **Chamfer all 12 edges**: `chamfer_edges` with `ChamferOptions{1.0, 1e-9}`.
  Analytic volume `a^3 - 6*a*d^2 + 6*d^3 = 1,000,000 - 600 + 6 = 999,406`.
- **Closed shell**: `shell_solid` with `ShellOptions{1.0, 1e-9}` and no open
  faces. Analytic volume `100^3 - 98^3 = 58,808`; one cavity shell.

## Results

| Operation | min (ms) | median (ms) | max (ms) | Volume | Check |
| --- | ---: | ---: | ---: | ---: | --- |
| Chamfer all 12 edges | 466.700 | 495.516 | 556.290 | 999405.999999998487 | within 1e-6 relative of 999406 |
| Closed shell | 1.105 | 1.158 | 1.596 | 58807.999999999709 | within 1e-6 relative of 58808 |

The benchmark exits 0 only when both checksums match, 2 on an operation failure,
and 3 on a checksum mismatch. The executable is built only with
`-DGEOM_BUILD_BENCH=ON` and is not part of ctest.

## Test evidence

- `*shell*`: 45 cases / 1,161 assertions in Release (the filter also matches the
  pre-existing `brep_shell_test.cpp` cases); Debug 45 / 465. New shell cases: 12.
- New QA cases: `chamfer and shell outputs tessellate watertight` (8 assertions)
  and `blend outputs are deterministic and leave inputs unchanged` (29
  assertions), 2 cases / 37 assertions in both configurations.
- Full Release: 1,115 cases / 681,668 assertions, 0 failed (786 s).
- Full Debug remains PENDING for the same eight long
  `brep_boolean_curved_test.cpp` half-space booleans documented in M4b-2; the
  Debug shell/QA filters were run separately.

## Scope

Accepted: a closed, convex, single-shell, all-planar solid, with straight
`LineCurve3` edges of exactly two incident faces for chamfer; shell requires the
faces to be axis-aligned (therefore a box) and no open face. Outputs are
validated, watertight when tessellated, deterministic, and the input is
unchanged.

Deferred to M4c-2: fillet (including variable radius), chamfer angle, open-face
shells, oblique faces, non-convex and curved operands, offset and sweep, and the
splitter/boolean face-order hardening. A faceted fillet is explicitly not
shipped as a substitute.

## Reproduce

```powershell
cmake -S . -B build-bench -DGEOM_BUILD_BENCH=ON -DGEOM_BUILD_TESTS=OFF -DGEOM_WERROR=ON
cmake --build build-bench --config Release --target geom_bench_m4c1
.\build-bench\Release\geom_bench_m4c1.exe
```
