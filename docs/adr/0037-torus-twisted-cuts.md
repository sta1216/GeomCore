# ADR-0037: Twisted torus cuts and the chart re-origin

Status: accepted for M3c-3c-4b3c-1 (2026-11-21); extends the torus contract in
ADR-0032 and the tolerance accounting in ADR-0036.

## Context

`boolean_solid_curved` accepted a whole torus against planes only when the plane is
perpendicular to the torus axis or contains it
(`curved_boolean_torus_plane_cut_supported`,
`include/geom/brep/boolean_solid_curved.hpp:1393-1412`). Both families cut the torus
in circles that meet the chart seams without a twisted crossing, so the splitter has
nothing to rejoin. Every other tilt was refused before pairing with
`unsupported "boolean_solid_curved: torus cut plane must be perpendicular to or
contain the torus axis"` (`boolean_solid_curved.hpp:2183-2186`).

That left a real hole. A plane through the centre, `z = 0.3x`, bisects the torus by
central symmetry, so the two halves each have volume `pi^2/4`, yet the plane meets
the tube obliquely. On the canonical torus chart the cut is a closed chain that
crosses the `v = 0` seam (the outer equator, a circle around the axis) at two
different transverse `u` coordinates: for `R = 1`, `r = 0.5` and `z = 0.3x`, at
`u = pi/2` and `u = 3pi/2`. `split_face_curved` windows a closed chain at a periodic seam and
rejoins the pieces by translating the windowed copy along the seam. The rejoin is
area-preserving only when the two crossings keep the same transverse coordinate
(`include/geom/brep/trim_face.hpp:1379-1384`); a twisted crossing would need a shear
in the chart. The splitter detects the mismatch per axis
(`trim_face_closed_cut_twisted`, `trim_face.hpp:1234-1277`) and refuses before
arranging anything.

The obvious repair is to close the uv chain by joining the two crossings with a
straight segment. A probe did exactly that; the resulting boundary contained two
chords of length `1.0` whose endpoints are the two crossings but whose interiors are
off the torus surface. The
uv area is conserved and the arrangement accepts the loop, but the 3D boundary is no
longer a curve on the torus, so tessellation and mass properties would be built on a
false edge. Conserving uv area is a two-dimensional statement; it is not evidence
that the 3D stitch is valid.

The seam itself is not geometry. `v = 0` is an arbitrary outer equator (a circle
around the axis); changing it only reparameterizes the same surface. The cut occupies a bounded set of
`v` values and leaves free bands, so a different meridian can be chosen that the cut
never crosses. For the fixture the SSI chains cover `v` in
`[0, 0.902491359]`, `[2.820817100, 3.462368208]` and `[5.379496164, 6.283185307]`
(radians, printed to nine decimals), leaving gaps
`(0.902491359, 2.820817100)` and `(3.462368208, 5.379496164)`. With the chart origin
moved into a gap, the cut becomes an ordinary periodic cut and the existing splitter
applies unchanged.

## Decision

- **Re-origin the chart instead of stitching.** Add `double v_origin = 0.0` to
  `Torus3` (`include/geom/surfaces/analytic.hpp:217`). `eval`, `du`, `dv` and
  `normal` evaluate at `vv = v_origin == 0.0 ? v : v + v_origin`, so the zero-origin
  path is bit-identical. `face_uv_project` subtracts `surface.v_origin` instead of
  wrapping it (`include/geom/brep/face_uv.hpp:225-226`). `make_torus` takes the same
  optional origin and builds the `u` seam at the new meridian: centre
  `center + axis*r*sin(v_origin)`, radius `R + r*cos(v_origin)`, and starts the `v`
  seam interval at `v_origin` (`include/geom/brep/primitives.hpp:588-675`). No new
  surface type and no stored geometry; the same `Torus3` carries the chart.

- **Choose the midpoint of the widest free band.** The helper
  `curved_boolean_rechart_twisted_torus` (`boolean_solid_curved.hpp:1862-1980`) runs
  SSI on each torus/plane face pair, cuts the curves into chains, and collects each
  chain's `v` interval. It merges touching intervals and calls
  `trim_face_periodic_free_midpoint(merged, 2*pi, 1e-9*(1+margin), 4)`
  (`trim_face.hpp:1287-1353`). That helper returns `nullopt` when the circle is fully
  covered, the widest free gap is shorter than `2*margin`, or more than four
  intervals are supplied. The returned midpoint becomes `v_origin`; the torus is
  rebuilt and must contain exactly one solid. Eligibility is deliberately narrow:
  exactly one torus face in the operand, one loop, loop uv area equal to the chart
  area within `1e-9*chart_area`, `v_origin == 0`, every other operand face planar,
  and at least one overlapping plane with `alignment = |dot(normal, axis)|` satisfying
  `alignment > 1e-3 && 1 - alignment > 1e-3`
  (`kCurvedBooleanTorusRechartAlignmentSlack`, `boolean_solid_curved.hpp:1827`). This
  is what keeps the analytic perpendicular/axis-containing families on their exact
  path and keeps near-aligned planes out.

- **Re-chart once, on an internal pass.** The driver calls the helper with each
  operand first and, on success, recurses into `boolean_solid_curved` with the last
  parameter `allow_recharted_torus = true` (`boolean_solid_curved.hpp:2110-2147`).
  The parameter defaults to `false` (`:2013`) and only those two recursion sites
  pass `true`; the guard is `supported ||
  (allow_recharted_torus && plane_cut_recharted(...))`, and
  `plane_cut_recharted` re-applies the `1e-3` near-aligned bands, so a re-chart
  raised by one oblique face cannot admit a near-perpendicular or near-axis plane
  in the same pass (`:2190-2193`). A
  second re-chart is impossible by construction: the helper requires
  `v_origin == 0`, and the re-charted torus has a nonzero origin. Any failure after
  the helper's re-chart decision (a tracing error, a rebuilt torus with the wrong
  solid count) is an `Error`, not a silent fallback to the old chart, for example
  `Errc::invalid_topology` with
  `"boolean_solid_curved: re-charted torus must have exactly one solid"`. Residual
  risk: because the flag is a defaulted argument on a public function, a caller that
  deliberately passes `true` bypasses the exact-family guard without a re-chart for
  genuinely oblique planes. It is an opt-in escape hatch, not a supported mode, and
  it still cannot admit the near-aligned bands; the driver, the public tests and the
  documented API keep it false.

- **Refine the seam crossing, then rebase the partner chain.** The seam point is
  interpolated between SSI samples and can miss the plane by the chord sag (about
  `2.5e-3` on the fixture). `curved_boolean_refine_torus_seam`
  (`boolean_solid_curved.hpp:939-974`) keeps `u`, solves the plane residual as a
  function of `v` with Newton (at most 32 iterations, slope
  `dot(plane.normal, dv)/|normal|`), leaves an exact crossing untouched
  (`|offset| <= 1e-12*(1+|v|)`), requires the final residual to be within the boolean
  tolerance, and normalizes `v` into `[0, 2*pi)`. The driver then rebases the partner
  chain onto the refined rim point (`curved_boolean_rebase_closed_chain`) and drops
  interior partner samples within `2*tolerance` of the rim start
  (`boolean_solid_curved.hpp:2358-2409`). Without this the two rims are not the same
  polyline and assembly fails with "result touches itself along an edge".

- **Rejected: stitching the seam.** Joining the two crossings with a segment
  conserves the uv area, but the probe's closure produced two chords of length `1.0`
  that are off the torus surface between their endpoints. The arrangement and the
  area census would accept it, and every downstream 3D consumer would inherit a
  boundary that is not on the surface. The correct operation is to move the chart so
  no stitch is needed.

- **The refusal class is explicit.** The new path must not become a silent partial
  success:

  - Whole torus against an all-planar partner, but the plane is near-perpendicular
    (`||dot| - 1| <= 1e-3`), near-axis (`|dot| <= 1e-3`), missed, or leaves no usable
    free `v` band; or the torus is not whole (multiple faces, several loops, holes,
    loop area not equal to the chart area); or the operand chart already has a
    nonzero `v_origin`: `unsupported`
    `"boolean_solid_curved: torus cut plane must be perpendicular to or contain the
    torus axis"`, before `report.face_pairs` is incremented.
  - A tangent plane is refused by the generic contact diagnostic
    (`"boolean_solid_curved: tangent face contact is not supported yet"`) after
    pairing, so `report.face_pairs` may be `1`; some tangent variants instead reach
    the splitter and are refused as open cuts. Both are explicit errors with no
    solid, but they are not part of the zero-count guard class.
  - Torus against a non-plane operand: `unsupported`
    `"boolean_solid_curved: torus operands only support plane cuts; torus/torus and
    torus/curved pairs are deferred to M3c-3"`.
  - Re-chart postcondition failure: `Error(Errc::invalid_topology,
    "boolean_solid_curved: re-charted torus must have exactly one solid")`.
  - Direct `split_face_curved` on a raw twisted chain: `unsupported`
    `"split_face_curved: closed cut crosses the seam at different transverse
    coordinates; twisted seam crossings are deferred to M3c-3c-4b3c"`
    (`trim_face.hpp:1393-1399`). The `4b3c` suffix is new; the message previously
    ended at `M3c-3c-4`.
  - Reversed difference: the pre-existing `unsupported`
    `"boolean_solid_curved: flipped curved fragments are not supported yet"`
    (`boolean_solid_curved.hpp:2220-2250`) still fires; the re-chart does not change
    it.

- **The tolerance is the polygonal SSI budget.** `ssi()` samples every patch at
  `chord_tolerance = max(1e-3 * diagonal, 64 * tolerance)`
  (`include/geom/brep/ssi.hpp:185`). The oblique cut is therefore a polyline
  inscribed in the true intersection, and moving the seam does not reduce the sag.
  The new cases assert `|V - V_analytic| <= 5e-3 * |V_analytic|` for each half and
  the union, and `|(V_intersection + V_difference) - V_torus| <= 1e-2`. The `5e-3` is
  the measured envelope of the fixed budget (worst measured `3.30e-3`); it does not
  relax any exact case. The perpendicular/axis-containing suite and every other
  curved test keep their previous tolerances, `1e-9` where the cut is analytic.

## Evidence

Fixture: `make_torus({0,0,0}, {0,0,1}, 1, 0.5, 1e-9)` against a closed slab whose
bottom face is `z = 0.3x` (8x8 footprint, 4 high), boolean tolerance `1e-9`. The
plane passes through the torus centre, so the intersection and difference are each
`pi^2/4 = 2.4674011002723395` and the union is
`256 + pi^2/4 = 258.46740110027235`.

| case | measured | analytic | error | relative |
| --- | ---: | ---: | ---: | ---: |
| forward intersection | 2.459259398081 | 2.4674011002723395 | -8.141702e-3 | -3.300e-3 |
| forward difference | 2.468547259432 | 2.4674011002723395 | +1.146159e-3 | +4.645e-4 |
| forward union | 258.468547259432 | 258.46740110027235 | +1.146159e-3 | +4.434e-6 |
| reversed intersection | 2.4622696393967654 | 2.4674011002723395 | -5.131461e-3 | -2.080e-3 |
| reversed union | 258.4815989254426 | 258.46740110027235 | +1.419783e-2 | +5.493e-5 |
| intersection + difference | 4.927806657513 | 4.934802200544679 | -6.995543e-3 | -1.418e-3 |

The re-chart selected `v_origin = 1.8620149558509749`, the full-precision midpoint
of the wider gap `(0.902491359, 2.820817100)`. The two printed gaps differ by about
`1.2e-3` rad, and a hand calculation from the nine-decimal chain ends gives
`1.8616542295`, `3.6e-4` away but still inside the same gap; the exact ends decide
the midpoint, so the test pins the engine value.

Each case validates with `validate_topology`, one shell and
`shell_is_closed_by_uses`; two runs are bit-identical in volume and counts, and
`report.face_pairs > 0`. Full suites are 1031 cases / 674,604 assertions in both
Debug and Release with `/W4 /WX /utf-8` and no warnings; the focused torus tests are
54 cases / 10,574 assertions.

## Consequences

- A whole torus (one shell, one torus face, no holes) against an all-planar operand
  is supported when the plane is neither perpendicular nor axis-containing (outside
  the `1e-3` alignment bands), the cut leaves at most four covered `v` intervals and
  a free band at least `2*margin` wide, and the cut is reachable by one re-chart. The
  result meets the `5e-3`/`1e-2` polygonal budget above, not an exact-volume
  guarantee.
- `v_origin = 0` stays bit-identical across the surface, projection, primitive and
  every existing boolean. `v_origin` is not a general-purpose surface parameter: it
  is set only by the re-chart helper on a rebuilt whole torus. A user-built torus
  with a nonzero origin cannot be re-charted and is refused by the public guard.
- Raw splitter chains remain refused. `split_face_curved` is not taught to stitch;
  the boolean driver resolves the twist before calling it.
- Reversed difference is deferred by the flipped-fragment limitation, which
  pre-dates this change; reversed intersection and union work because their
  fragments keep the expected orientation.
- Deferred to M3c-3c-4b3c-2: torus/torus and torus/curved pairs, tangent planes,
  cuts with more than two seam crossings or more than four covered bands, multiple
  re-charts, and the cone axis-plane/apex-ring case. The cone work is the other half
  of 4b3c and is not covered by this ADR.
- The volumes are polygonal SSI volumes, not exact analytic ones. No performance or
  accuracy claim is made beyond the table; a change to the SSI chord budget moves
  these numbers and requires re-measuring this note.
