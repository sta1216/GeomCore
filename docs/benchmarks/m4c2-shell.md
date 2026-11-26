# M4c-2 open-shell and oblique-shell benchmarks

- Date: 2026-11-26
- Commit: `sdd/m4c2` (Task 7 evidence; the measured tree is `4575ce0`)
- Compiler: MSVC 19.44.35217.0, `/std:c++20 /O2`, `GEOM_WERROR=ON`
- Machine: Windows 11, Intel Core i7-8550U (4C/8T), Release
- Benchmark: `bench/bench_m4c1.cpp`, five timed runs per workload, min/median/max

## Workloads

All three use the same binary. The 100^3 box is created with
`make_box({0,0,0}, {100,100,100}, 1e-9)`.

- **Closed shell**, `ShellOptions{1.0, 1e-9}`, no open faces: analytic volume
  `100^3 - 98^3 = 58,808`, one cavity shell.
- **Open top shell**, `ShellOptions{1.0, 1e-9}`, the `+Z` face is open: the
  cavity protrudes through `+Z`, so the removed block is `98 * 98 * 99 =
  950,796` and the shell volume is `1,000,000 - 950,796 = 49,204`; the result
  has `cavity_shells == 0` and one shell.
- **Oblique tetrahedron**, legs `a = 2`, `ShellOptions{0.1, 1e-9}`, no open
  faces: analytic volume
  `4/3 - (2 - 0.1 * (3 + sqrt(3)))^3 / 6 = 0.740147367097487`.

## Results

| Workload | min (ms) | median (ms) | max (ms) | Volume | Check |
| --- | ---: | ---: | ---: | ---: | --- |
| Closed shell, 100^3, t=1 | 1.097 | 1.154 | 1.337 | 58807.999999999709 | within 1e-6 relative of 58,808 |
| Open top shell, 100^3, t=1 | 1.372 | 1.454 | 1.845 | 49203.999999999942 | within 1e-6 relative of 49,204 |
| Oblique tetra, a=2, t=0.1 | 0.426 | 0.434 | 0.452 | 0.740147367097487 | within 1e-6 relative of the analytic value |

The table is the Task 7 evidence run, five timed iterations per workload on
an idle machine; the reported volume is the last iteration's value and the
checksum ran on every iteration.

The binary exits 0 only when every checksum matches; 2 is an operation
failure, 3 is an open-shell failure or the legacy chamfer/closed-shell
checksum, 4 is an open-top volume mismatch, 5 is a tetra operation failure and
6 is a tetra volume mismatch. The executable is built only with
`-DGEOM_BUILD_BENCH=ON` and is not part of ctest. The same binary still runs
the M4c-1 chamfer and closed-shell workloads; this record tabulates only the
M4c-2 shell cases and does not reuse M4c-1 timings.

## Reproduce

```powershell
cmake -S . -B build-bench -DGEOM_BUILD_BENCH=ON -DGEOM_BUILD_TESTS=OFF -DGEOM_WERROR=ON
cmake --build build-bench --config Release --target geom_bench_m4c1
.\build-bench\Release\geom_bench_m4c1.exe
```

CI runs the same binary and uploads `bench-m4c1.txt` with the other benchmark
artifacts.
