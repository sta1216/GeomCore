# M3c-3c-2b benchmark note

## Scope

M3c-3c-2b lets a horizontal offset torus plane split in all three boolean
operations and replaces the two zero-area rim rings of a periodic-v torus band
with one whole-cycle ring. It adds no benchmark executable and no timing gate: the
acceptance oracles are exact volumes, fragment/ring counts, `validate_topology`,
and deterministic repeats. Curved-boolean timing is unchanged and covered by
`docs/benchmarks/m3c3c.md`; this note records only the numbers measured for this
milestone.

## Environment

- MSVC 19.44 (Visual Studio 17 2022, toolset `14.44.35207`), Windows x64.
- `GEOM_WERROR=ON`, `/W4 /WX /utf-8`; zero warnings in Debug and Release.
- Boolean acceptance tolerance `1e-9`; splitter/ring tests `1e-12`.
- No timing gate applies; no benchmark was run for this task.

## Offset torus plane

Fixture: torus `R=1`, `r=0.5` at the origin, cut by the box
`{-2,-2,0.25}..{2,2,2}`; `CurvedBooleanOptions{1e-9}`; the volume comparison uses
relative `epsilon(1e-6)`.

| operation | faces | edges | vertices | shells | volume |
| --- | ---: | ---: | ---: | ---: | --- |
| intersection | 2 | 24 | 23 | 1 | `0.9647593052603948` |
| union | 9 | 122 | 116 | 1 | `31.970042895284287` |
| difference | 3 | 110 | 108 | 1 | `3.970042895284284` |

All three pass `validate_topology(...).ok()`, and a second intersection run
compares bitwise equal to the first (`volume_of(cut_again) == volume_of(cut)`).

## Ring contract

The torus band test cuts the full chart at `v = 0.6` and `v = 2.5` and expects
three fragments, one uv loop each, each ring longer than 100 points and starting
and ending on the `u = 0` seam column at different `v` (the closing column, not a
duplicated seam vertex). The ring areas equal `2*pi*height`:

- `3.7699111843077517` (`2*pi*0.6`);
- `11.938052083641214` (`2*pi*1.9`);
- `23.770454336408466` (`2*pi*(2*pi-2.5)`);

measured within `1.4e-14` of the exact values, against a relative test tolerance
of `1e-9`. The public offset split (`v = pi/6`, `v = 5*pi/6`) returns the same
three rings as the direct `trim_face_split_domain` call, with areas
`3.289868133696453`, `13.159472534785811`, `23.029076935899213`. A torus face
with no cuts is returned unchanged as one fragment of area `4*pi^2`.

## Refusals pinned

| input | result |
| --- | --- |
| stored inner loop on a sphere band | `unsupported "split_face_curved: stored holes are not supported yet"`, `context_b == 1u` |
| twisted u-seam crossing polygon `{(1,1),(6,1.1),(0.2,1.2),(3,1.15),(1,1)}` | `unsupported "split_face_curved: closed cut crosses the seam at different transverse coordinates; twisted seam crossings are deferred to M3c-3c-4"` |
| cut run on the parameter-domain side `u = 0` | `unsupported "split_face_curved: cut lies on the parameter domain boundary"` |
| interior island set on a periodic face | `unsupported "split_face_curved: closed cuts on a periodic face are not supported yet"` |
| trimmed `Torus3` patch | `unsupported "split_face_curved: trimmed torus faces are not supported yet"` |
| tilted or axis-parallel offset torus plane | `unsupported "boolean_solid_curved: torus cut plane must be perpendicular to or contain the torus axis"`, `report.face_pairs == 0` |

The tilted and axis-parallel offset refusals are pre-existing and were not
duplicated. The stored-hole public test, the no-cut torus test, and the
axis-parallel offset test also predate Task 4; Task 4 adds only the twisted
u-seam and domain-boundary pins.

## Full suites

- Debug: 983 cases / 660,301 assertions, 0 failed, exit 0
  (`build\tests\Debug\geom_tests.exe`).
- Release: 983 cases / 660,301 assertions, 0 failed, exit 0
  (`build-rel\tests\Release\geom_tests.exe`).
- Task 4 adds 2 cases / 11 assertions over the 981 / 660,290 measured at
  `6d77096`.

## RED/GREEN notes

- Task 1 RED (before `ed8caa0`): the torus latitude split emitted two 66-point
  zero-area rim rings instead of one ring per band; the test failed on
  `REQUIRE(fragment.loops_uv.size() == 1u)` (`2 == 1`) and
  `CHECK(fragment.loops_uv.front().size() > 100u)` (`66 > 100`). GREEN after the
  whole-cycle branch.
- Task 2 RED (before `8cb5b4c`): the public torus split returned one unchanged
  fragment instead of three (`REQUIRE(1 == 3)`) and accepted an interior island
  instead of refusing it. GREEN after the `trim_face_domain` dispatch.
- Task 3 RED (header at `ed8caa0`, before the dispatch): the offset intersection
  failed with `boolean_solid_curved: result touches itself along an edge` at
  `tests/unit/brep_boolean_curved_test.cpp:1148`. GREEN after `1738878`.
- Task 4 is test-plus-docs only and has no production change: the two new pins
  pass immediately against `6d77096`. Their RED evidence is the Tasks 1–3 record
  above; no old binary is used for the counts.
