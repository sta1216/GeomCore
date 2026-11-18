# M3c-3c-4a benchmark note

## Scope

M3c-3c-4a lets a curved trimmed face with rectilinear hole rings mesh watertight and
weigh. It adds no benchmark executable and no timing gate: the acceptance oracles are
the exact lateral uv area, the controlled mesh volume band, the analytic mass volume,
the watertight edge count, and two pinned refusals. This note records only the numbers
measured for this milestone; `docs/benchmarks/m3c3c2b.md` covers the curved-boolean
side and no timing is claimed here.

## Environment

- MSVC 19.44 (Visual Studio 17 2022, toolset `14.44.35207`), Windows x64.
- CMake 4.4.4; `GEOM_WERROR=ON`, `/W4 /WX /utf-8`; zero warnings in Debug and Release.
- Mesh acceptance fixture options `TessellationOptions{1e-3, 0.35}` (fine) and
  `{1e-2, 0.35}` (coarse); mass options `MassPropertiesOptions{1e-9, 24, 1e-5}`.
- No timing gate applies; no benchmark was run for this task.

## Fixture

Cylinder `r = 1`, `h = 2` at the origin cut by the box
`{-0.5,0.5,0.5}..{0.5,2,1.5}`, tolerance `1e-9`: the window is centred on the `+y`
wall, does not cross the periodic cut, and produces 8 faces with a valid topology
(`validate_topology(...).ok()`). The open window removes the block
`0.45661147746179975`, so the exact surviving volume is
`2*pi - 0.45661147746179975 = 5.826573829717787`. The exact lateral uv area (the
cylinder side still covered by the face) is `4*pi - pi/3 = 11.519173063162574`.
The fixture lives in `tests/unit/tessellate_trimmed_test.cpp` (`hole_window_cylinder`).

## Mesh

| quantity | fine `{1e-3, 0.35}` | coarse `{1e-2, 0.35}` |
| --- | ---: | ---: |
| triangles | 1236 | 452 |
| positions | 620 | 228 |
| boundary edges | 0 | 0 |
| non-manifold edges | 0 | 0 |
| lateral uv area | `11.519173063162778` | `11.519173063162546` |
| lateral uv area deviation | `2.0e-13` (recording tolerance `1e-6`; test asserts `1e-12` relative) | `2.8e-14` (same) |
| mesh volume | `5.8238533084266857` | `5.8001219438259604` |
| volume against analytic | inside `[5.8200, 5.8266]` (acceptance band, no unit test) | below `[5.8200, 5.8266]` (not gated) |

The fine volume is inside the brief's band; the coarse volume is honestly below it, as
expected from the coarser chordal approximation (the coarse case only gates that the
hole stays open and the mesh stays watertight, not the volume band). The fine triangle
count is inside the `[950, 1450]` window asserted by the Task 1 case; both meshes reuse
the neighbouring faces' boundary vertices, which is why the uv area is preserved to
`1e-13`. The Task 1 uv-area checks assert `epsilon(1e-12)` relative; the `1e-6` figure
in the table is the recording tolerance from the brief.

## Mass properties

| quantity | measured | gate |
| --- | --- | --- |
| volume | `5.8265738271299785` | `4.44e-10` relative against `5.826573829717787`, gate `1e-9` |
| centroid | `(-2.7426505511130017e-11, -0.05720914955052997, 0.99999999999999789)` | moves to `-y`, away from the `+y` window |

The mass integrator simplifies the outer ring plus the hole ring and triangulates the
combined domain in one call; the measured volume is within `2.6e-9` absolute of the
analytic value.

## QA case

`QA-20: the holed window cylinder meshes and weighs consistently`
(`tests/unit/tessellate_trimmed_test.cpp`) meshes with `{1e-3, 0.35}`, checks
watertightness, then weighs the same `Brep` with `{1e-9, 24, 1e-5}` and asserts the
volume within `1e-6` relative of `5.826573829717787`. It adds no assertion that
Tasks 1 and 2 do not already pin; it only requires both operations to succeed on the
same input in one test.

- Focused run (`--test-case="QA-20*"`, Debug): 1 case / 1,857 assertions / 0 failed.
- Full suites: Debug 986 cases / 668,843 assertions / 0 failed, exit 0
  (`build\tests\Debug\geom_tests.exe`); Release 986 cases / 668,843 assertions /
  0 failed, exit 0 (`build-rel\tests\Release\geom_tests.exe`).
- The +1,857 assertions over the 985 / 666,986 baseline are dominated by
  `check_watertight`'s per-edge `CHECK` over the 1,236-triangle fine mesh.

## Refusals pinned

| input | result |
| --- | --- |
| oblique trim, plane normal `(0.3, 0, 1)` | tessellate `unsupported "tessellate: trimmed curved face boundary is not rectilinear"`; mass `unsupported "mass_properties: trimmed curved face boundary is not rectilinear"` |
| window straddling the artificial cut, box `{-2,-0.5,0.5}..{-0.5,0.5,1.5}` | tessellate `unsupported "tessellate: trimmed curved face hole crosses the periodic cut"`; mass succeeds with volume `5.8265738271295096` |

Both refusals are measured on constructible inputs. The second input is new evidence
from this milestone: the window is centred on the `-x` wall (`u = pi`) and its samples
straddle the artificial cut near `u = 3.097`, so the mesher refuses it per piece
while the mass path, which does not split the periodic ring, still weighs it to the
same analytic volume. It is recorded as an intended asymmetry, not as a test case; the
milestone's only new QA case is the strictly-inside window above.

Not measured / no fixture:

- The non-rectilinear **hole** refusal (`tessellate: trimmed curved face hole boundary
  is not rectilinear`, `tessellate.hpp:1747`) has no constructible fixture. A
  perpendicular cylindrical cutter (base `(0,0.5,1.0)`, axis `(0,1,0)`, radius `0.3`,
  height `1.5`) fails upstream in the boolean with
  `boolean_solid_curved: flipped curved fragments are not supported yet`, before a hole
  reaches the mesher. The refusal is pinned by code, not by this note.
- Stored holes on synthetic faces, twisted seam crossings, and NURBS trimmed contours
  are not exercised here; they keep the refusals recorded in ADR-0032 and ADR-0027.
- A second hole that unwraps outside its split piece is refused: cylinder minus two
  boxes (windows on +y at z 0.5..1.5 and on -y at z 0.5..1.5) returns
  `tessellate: trimmed curved face hole crosses the periodic cut`, while
  `mass_properties` succeeds at 5.3699623498711313. A same-piece two-hole face
  meshed in a probe (2098 triangles, no bad edges) but is neither tested nor
  claimed here; normalising each hole ring into its piece is a follow-up.
- No timings were taken; there is no timing gate and no benchmark binary.

## Full suites

- Debug: 986 cases / 668,843 assertions, 0 failed, exit 0.
- Release: 986 cases / 668,843 assertions, 0 failed, exit 0.
- Task 3 adds 1 case / 1,857 assertions over the 985 / 666,986 measured at
  `1b1eefa`.

## RED/GREEN notes

- Task 3 is test-plus-docs only and has no production change: `QA-20` passes
  immediately against `1b1eefa` (the implementation is `b120edd` + `3d5dc8d`). Its RED
  evidence is the Tasks 1 and 2 records above; the old blanket refusals
  (`tessellate: trimmed curved face has a hole`,
  `mass_properties: trimmed curved face has a hole`) are what the QA case would have
  hit before those commits.
- The counts above are from the committed tree at the Task 3 commit; no old binary was
  used.
