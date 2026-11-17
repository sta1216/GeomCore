# ADR-0032: General synthetic-domain splitter (periodic-v bands included)

Status: accepted for M3c-3c-2b (2026-11-18)

## Context

`split_face_curved` treated a curved face as *seam-only* when its stored uv loop had
zero area, then opened the surface onto a synthetic parameter rectangle. That worked
for sphere caps and bands but had three defects: a synthetic chart that already
carried a stored inner ring was split on its outer ring only (the hole was silently
dropped); whether a chart was synthetic was decided from the stored loop's area,
while a torus chart is the whole rectangle (`4*pi^2`), so the torus never took the
synthetic path; and torus bands that wrap the periodic-v seam had no representation
in the seam-only emitter: their two halves were emitted as zero-area slit rings or
as rings missing the v-seam rim, so assembly could not pair them.

M3c-3c-2a introduced `TrimFaceDomain` and the arrangement-based
`trim_face_split_domain`, routed sphere synthetic charts through it, and refused
stored holes, nested cuts, overlapping islands, and twisted seam crossings. It left
offset torus planes blocked: probes that shifted the v origin onto a cut latitude
and that emitted full-width v-side seam rings both failed the offset acceptance test
and regressed five centred-plane torus booleans (`touches itself along an edge`,
`arc run is inconsistent`, `not a single connected component`, plus one wrong error
code). The measured blocker was the two-rim band representation, not the
arrangement: harvested fragments carried healed spans (e.g. ring sizes 46/46, 45/49,
69/69), but the rims of adjacent torus bands and the planar annulus could not be
paired.

## Decision

- `split_face_curved` dispatches a `Torus3` face through `trim_face_domain` and
  `trim_face_split_domain`. The sphere keeps its historical seam-only wrapper (the
  zero-area / rounding-slit gate); its stored holes are still refused with
  `"split_face_curved: stored holes are not supported yet"` and
  `context_b = loops.size() - 1`.
- The torus router refuses an interior island set before the arrangement
  (`"split_face_curved: closed cuts on a periodic face are not supported yet"`) and
  a trimmed patch whose stored outer loop does not cover the full chart within
  `1e-9` relative area (`"split_face_curved: trimmed torus faces are not supported
  yet"`). A face with no cuts and no holes falls through to the ordinary path and is
  returned unchanged instead of building an empty synthetic arrangement.
- A periodic-v band carries exactly one whole-cycle ring. The ring starts at
  `(0, v0)`, runs along the lower rim to `(2*pi, v0)` with `+u`, turns up the right
  seam column to `(2*pi, v1)` with `+v`, returns along the upper rim to `(0, v1)`
  with `-u`, and closes down the left seam column with `-v`, ending at
  `(0, v0+delta)` on the same seam copy rather than at a duplicated vertex. The
  ring is CCW with area `2*pi*(v1-v0)`; no band is emitted as two zero-area rim
  halves and no full-width seam ring remains. Sphere bands keep their existing
  one-ring form.
- Dispatch changes where the public result comes from, not the ring contract: the
  public split of the offset torus equals the direct `trim_face_split_domain`
  result ring by ring.
- A horizontal offset plane is accepted by the torus plane guard: the cut plane
  must be perpendicular to the torus axis (`z = c`, any offset) or contain it.
  Tilted planes and axis-parallel offset planes are still refused with
  `"boolean_solid_curved: torus cut plane must be perpendicular to or contain the
  torus axis"` and `report.face_pairs == 0`.
- On the synthetic paths, stored holes, nested/mixed cuts, overlapping islands, and
  cuts on the parameter domain boundary are explicit errors, never silent losses:
  `"...stored holes are not supported yet"`, `"...nested closed cuts are not
  supported"`, `"...fragments do not partition the chart"`, and
  `"split_face_curved: cut lies on the parameter domain boundary"`. A non-synthetic
  face with stored holes still follows the ordinary clipping path and ignores them,
  exactly as before this milestone.
- Twisted seam crossings are still refused with `"split_face_curved: closed cut
  crosses the seam at different transverse coordinates; twisted seam crossings are
  deferred to M3c-3c-4"`. They need a chart re-origin, not merely a different ring;
  the offset acceptance test exercises only the same-transverse case.
- Tessellation of a trimmed periodic-v face is still unsupported. The measured gaps
  are `tessellate: torus face is missing a seam` (`tessellate.hpp:890`, the torus
  face's uv loop does not carry the periodic seam) and
  `tessellate: trimmed band wall is not closed` (`tessellate.hpp:1342`, the wall of
  a band with a periodic-v trim is not closed). This milestone adds no tessellation
  test and does not claim meshing support.

## Evidence

The acceptance fixture is a torus `R=1`, `r=0.5` at the origin cut by the box
`{-2,-2,0.25}..{2,2,2}`, tolerance `1e-9`, MSVC 19.44 (`/W4 /WX /utf-8`, zero
warnings):

| operation | faces | edges | vertices | volume (epsilon `1e-6`) |
| --- | ---: | ---: | ---: | --- |
| intersection | 2 | 24 | 23 | `0.9647593052603948` |
| union | 9 | 122 | 116 | `31.970042895284287` |
| difference | 3 | 110 | 108 | `3.970042895284284` |

All three result in one shell, pass `validate_topology`, and the volume of a
second intersection run matches the first bitwise. The torus band contract is pinned separately: two
latitude cuts at `v = 0.6` and `v = 2.5` yield three one-ring fragments with
positive areas `2*pi*height` (measured within `1.4e-14`), and the public offset
split at `v = pi/6` and `v = 5*pi/6` yields the same three rings and areas as the
direct splitter (`3.289868133696453`, `13.159472534785811`,
`23.029076935899213`, relative test tolerance `1e-9`). The no-cut torus face is
returned unchanged as one fragment of area `4*pi^2`.

Task 4 pins the twisted u-seam crossing and the domain-boundary cut and re-measures the suites: Debug and Release each report 983
cases / 660,301 assertions, zero warnings, with `GEOM_WERROR=ON`. The counts and
environment are recorded in `docs/benchmarks/m3c3c2b.md`.

## Consequences

- Horizontal offset torus planes now work in intersection, union, and difference
  with the measured volumes above; centred-plane torus booleans and the six guard
  refusals are unchanged.
- Only a canonical full-chart `Torus3` face is split. A trimmed `Torus3` patch is
  refused instead of reopened; a trimmed patch that carries no cuts and no holes is
  still returned unchanged (a pre-existing edge, unchanged here).
- A closed cut that touches a synthetic side is not classified as an island
  (`island = closed && no point on any side`), so the all-islands partition guard
  does not run: the cut is split and the complement ring is pinched at the contact
  point. A 32-gon of radius 0.2 tangent to `u = 0` at `(0, 0.5)` yields two
  fragments of areas `39.353559798267149` and `0.12485780609032211`. An interior
  island set on the periodic router is refused by the `islands_only` pre-check with
  `split_face_curved: closed cuts on a periodic face are not supported yet`; the
  general splitter refuses nested or mixed island cuts with
  `split_face_curved: nested closed cuts are not supported`.
- Sphere behaviour, including the seam-only wrapper and every refusal message, is
  byte-for-byte unchanged.
- Deferred to M3c-3c-4: multi-cycle holes, cone apex topology, general twisted
  crossings, and tessellation of trimmed periodic-v faces.
