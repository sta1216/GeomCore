# M3c-3c-4b3c-1: twisted torus cuts

Date: 2026-11-21. The volumes were measured with the development probe at `073f266`;
the guard fix `819a6c5` does not change this fixture's accepted path, and the final
tree's tests pin the same numbers at `5e-3`. The full-suite counts were measured on
the final tree with the Debug and Release builds of that tree.

## Scope

This note records the measured behavior of the twisted-torus re-chart:

- a whole torus cut by a plane that is neither perpendicular to the axis nor
  axis-containing now succeeds when the chart can be re-origined into a free `v`
  band (ADR-0037);
- the forward intersection, difference and union are checked against the analytic
  half-torus and slab volumes within the polygonal SSI budget;
- the reversed intersection and union are recorded as measured, while reversed
  difference remains refused by the pre-existing flipped-fragment limitation;
- the standalone `split_face_curved` still refuses a raw twisted chain.

The matrix is pinned by `curved boolean: torus cut by a tilted plane is re-charted`
(`tests/unit/brep_boolean_curved_test.cpp:1742`). The splitter refusal is pinned by
`tests/unit/brep_trim_face_test.cpp:3046`, and the free-band helper by the
`trim_face_periodic_free_midpoint` cases at `tests/unit/brep_trim_face_test.cpp:3052`.

## Environment

- Windows 10 Enterprise 10.0.19045, build 19045, x64.
- Intel Core i7-8550U, 4 cores / 8 threads, 15.9 GB RAM.
- MSVC 19.44 (Visual Studio 17 2022, toolset 14.44.35207), CMake 4.4.4.
- Debug and Release, `GEOM_WERROR=ON`, `/W4 /WX /utf-8`.
- Single run per configuration on the host above; there is no timing gate.

## Fixture

The torus is `make_torus({0,0,0}, {0,0,1}, 1, 0.5, 1e-9)`: major radius `1`, minor
radius `0.5`, canonical `v` origin before the re-chart. The cutter is a closed slab
with an 8x8
footprint in `x` and `y`, extending from `z = 0.3x` to `z = 0.3x + 4`; its bottom
face has normal `(0.3, 0, -1)/sqrt(1.09)` and contains the torus centre. Boolean
tolerance is `1e-9` for every case.

The plane passes through the torus centre, so by central symmetry the intersection
and difference are each `pi^2/4 = 2.4674011002723395`; the torus is
`pi^2/2 = 4.934802200544679`, and the slab is `256`. The union is
`256 + pi^2/4 = 258.46740110027235`.

Reversed cases put the slab first, so the torus contributes the second curved
operand. Volume error is `measured - analytic`; the relative error is
`|error| / analytic`.

## Matrix

| case | measured volume | analytic | error | relative |
| --- | ---: | ---: | ---: | ---: |
| forward intersection | 2.459259398081 | 2.4674011002723395 | -8.141702e-3 | -3.300e-3 |
| forward difference | 2.468547259432 | 2.4674011002723395 | +1.146159e-3 | +4.645e-4 |
| forward union | 258.468547259432 | 258.46740110027235 | +1.146159e-3 | +4.434e-6 |
| reversed intersection | 2.4622696393967654 | 2.4674011002723395 | -5.131461e-3 | -2.080e-3 |
| reversed union | 258.4815989254426 | 258.46740110027235 | +1.419783e-2 | +5.493e-5 |
| intersection + difference | 4.927806657513 | 4.934802200544679 | -6.995543e-3 | -1.418e-3 |

The intersection and difference are complementary halves: their sum is `6.996e-3`
below the torus volume, inside the `1e-2` sum check. The reversed difference is not
measured because it is refused (see Refusals). The reversed rows were measured on
the development probe's slab, whose torus side sits about `3e-6` relative away from
the committed `make_tilted_slab` fixture; the forward rows reproduce on that fixture
to better than `1e-9`, so the `5e-3` budget is unaffected.

## Tolerances

The new cases use `|V - V_analytic| <= 5e-3 * |V_analytic|` per volume and
`|(V_intersection + V_difference) - V_torus| <= 1e-2`. This is not a relaxation of
the exact suites: the perpendicular/axis-containing torus matrix and every
axis-aligned case keep their previous tolerances, `1e-9` where the cut circle is
analytic, and the topology checks remain `validate_topology`, one shell and
`shell_is_closed_by_uses`.

The reason is the fixed SSI chord budget. `ssi()` samples each patch with
`chord_tolerance = max(1e-3 * diagonal, 64 * tolerance)`
(`include/geom/brep/ssi.hpp:185`), so the oblique intersection is a polyline
inscribed in the true curve. The re-chart moves the seam but does not change that
sampling, and the resulting volume error is a chord-sag effect. The measured worst
half-volume relative error is `3.30e-3`, so `5e-3` is the measured envelope with
some headroom, not a target. The `1e-2` sum check is an absolute volume check on a
`4.9348` total, about `2.0e-3` relative.

## Band selection

The SSI chains cover the following `v` ranges (radians, printed to nine decimals):

| chain | `v` range |
| --- | --- |
| outer | `[0, 0.902491359]` |
| inner | `[2.820817100, 3.462368208]` |
| outer | `[5.379496164, 6.283185307]` |

Merging the `0`/`2*pi` copies leaves two free gaps: `(0.902491359, 2.820817100)`
and `(3.462368208, 5.379496164)`. The first is wider by about `1.2e-3` rad. The
helper returns its full-precision midpoint, which is the chosen
`v_origin = 1.8620149558509749`. A hand calculation from the nine-decimal printed
ends gives `1.8616542295`, `3.6e-4` away and still inside the same gap; the exact
chain ends decide the midpoint, so the test pins the engine value instead of a
recomputed one.

The chosen origin is stored on the rebuilt `Torus3` and on the rebuilt primitive's
seam edges; the boolean recursion then sees a canonical cut. `v_origin = 0` remains
bit-identical, and no existing axis-aligned case changes.

## Refusals and determinism

- Direct `split_face_curved` on the raw twisted chain still returns
  `unsupported "split_face_curved: closed cut crosses the seam at different
  transverse coordinates; twisted seam crossings are deferred to M3c-3c-4b3c"`
  (`include/geom/brep/trim_face.hpp:1393-1399`). The re-chart is in the boolean
  driver, not in the standalone splitter.
- A user-built torus with `v_origin = 4.0`, and a torus already re-charted to
  `1.8620149558509749`, are both refused by the public torus/plane guard with
  `"boolean_solid_curved: torus cut plane must be perpendicular to or contain the
  torus axis"` and `report.face_pairs == 0`.
- Near-axis and near-perpendicular planes (`1e-3` alignment bands, including normals
  tilted from the axis by `5e-10`, `5e-13` and `1e-12`) are refused by the same
  guard; a tilted plane that misses the torus is refused there as well when the
  operand pair overlaps (a fully disjoint pair reports an empty result instead).
- Reversed difference is refused by the pre-existing
  `"boolean_solid_curved: flipped curved fragments are not supported yet"`; the
  limitation predates the re-chart and is out of scope for 4b3c-1.
- Every accepted case is validated with `validate_topology`, has one shell and
  `shell_is_closed_by_uses`, and two runs produce bit-identical volume and counts.
  `report.face_pairs > 0` on the accepted cases.

## Focused and full runs

- Focused torus cases: 54 cases / 10,574 assertions, all passing.
- Full suites, Debug and Release: 1031 cases / 674,604 assertions, all passing,
  with `/W4 /WX /utf-8` and no warnings.
- The Task 2 baseline was 1023 cases / 674,283 assertions; Task 3 and its guard
  fixes add 8 cases and 321 assertions (the re-chart matrix, the miss and regression
  refusals, the reversed-order runs, and the near-aligned guard regression). Against
  the whole-milestone baseline (1014 cases / 673,425 assertions) the delta is
  17 cases / 1,179 assertions.

## Limits

- Only a whole torus (one shell, one torus face, no holes) against an all-planar
  operand is supported. Torus/torus, torus/cylinder and curved/curved pairs still
  return the plane-cuts-only message.
- The re-chart allowance is per pair and re-applies the `1e-3` near-aligned bands, so
  an oblique face that triggers a re-chart cannot widen the refusal of a
  near-perpendicular or near-axis plane in the same solid.
- The plane must be outside the `1e-3` near-aligned bands, must not be tangent, and
  must leave a free `v` band at least `2*margin` wide with at most four covered
  intervals. Cuts with more than two seam crossings are not supported.
- Only one re-chart is attempted; a nonzero `v_origin` cannot be re-charted.
- Reversed difference remains refused, and `split_face_curved` does not accept raw
  twisted chains.
- The cone axis-plane/apex-ring case is deferred to M3c-3c-4b3c-2 and is not
  measured here.
- These are polygonal SSI volumes. No performance claim, no exact-volume claim and
  no accuracy claim beyond the table are made.
