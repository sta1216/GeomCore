# ADR-0035: Cone operands in the curved boolean

Status: accepted for M3c-3c-4b2 (2026-11-19)

## Context

M3c-3c-4b1 (ADR-0034) gave the cone a synthetic chart, an apex projection and open
cuts that end on domain sides. The boolean driver still refused cone operands, and
three downstream paths still had cone guards:

- `curved_boolean_surface_gate` returned `"boolean_solid_curved: cone operands need
  closed-cut and apex topology; booleans are deferred to M3c-3"`.
- `curved_solid_volume` returned `"curved_solid_volume: cone faces are not supported
  yet"` because a loop that collapses at the apex closes with a chord across the
  chart corner, so the apex fan was not integrated.
- `trimmed_face_boundaries` returned `"trimmed_face_boundaries: cone faces are not
  supported yet"` for the same collapsed side.
- `assemble_curved_solid` saw two uv copies of the apex in one ring (the cut end and
  the back-seam copy). Both weld to one 3D vertex, leaving a zero-length segment whose
  edge key is `v->v` and which is used once, so the two halves did not assemble.

M3c-3c-4b2 opens the boolean for one narrow class: a whole canonical cone cut by
planes. QA-22 is the end-to-end pin.

## Decision

### Supported class

A cone operand is accepted only when all of the following hold:

- it is a whole cone: a single shell whose lateral `Cone3` face has exactly one stored
  loop, that loop spans the full `u` period, and its `v` range reaches both the base
  row and the collapsed apex row;
- the other operand is entirely planar, so every cut of the cone is a plane cut;
- every cone/plane pair is either perpendicular to the cone axis (within `1e-9`) or
  provably disjoint (see the conservative predicate below).

No public signature changes. The operations are the existing union, intersection and
difference. A cone that is one face of a larger B-rep, a trimmed lateral face, or a
cap/band fragment is not a whole cone and stays refused.

### Whole-cone guard and analytic chart

`curved_boolean_face_info` (`include/geom/brep/boolean_solid_curved.hpp:296-325`)
keeps a cone only when `info.loops.size() == 1`, the `u` span covers the analytic
period, and the stored `v` endpoints match `0` and `radius/tan(half_angle)` within

```
slack = 1e-9 * max(1, 2*pi, |v_hi| + |v_lo|)
```

Otherwise it returns, with the face id in the error context,

`unsupported "boolean_solid_curved: cone operands must be whole cones; trimmed cone fragments are not supported"`

For an accepted cone the chart is replaced by the analytic rectangle
`[0, 2*pi] x [0, radius/tan(half_angle)]` and marked as a synthetic domain. The
base rim is a straight chord in the seam-only stored loop, while the cone's base cap
already carries exact arcs. Harvest therefore adds the operand cone's base circle

```
CircleCurve3{cone.origin, cone.x_axis, cross(cone.axis, cone.x_axis), cone.radius}
```

to every fragment's heal circles (`boolean_solid_curved.hpp:2152-2163`). Both sides
of the base rim then re-engage the same analytic arc run and the same subdivision,
instead of one side being chordal and the other exact.

### Analytic cone/plane cut circle

`split_solid_face_circle` (`include/geom/brep/split_solid.hpp:460-479`) adds a
`Cone3` branch. For a plane `n . x = offset` with `|dot(n, axis)| >= 1 - 1e-12`:

```
t   = (offset - dot(n, origin)) / dot(n, axis)
rho = radius - t * tan(half_angle)
```

If `rho > 0` and finite, the cut rim is the exact circle at `origin + axis*t`
projected onto the plane, with axis `+n`, and with `x_axis` taken from
`cone.x_axis` projected into the plane (`frame.u` as the degenerate fallback).
`y_axis = cross(n, x_axis)`. The same `CircleCurve3` object goes to both faces of
the pair, so both fragments and both caps agree on one circle and one arc
subdivision. This is what makes the QA-22 volumes exact to `1e-12`; a chordal rim
would miss by orders of magnitude. A plane at or beyond the apex (`rho <= 0`), a
non-perpendicular plane, and a degenerate normal or axis return `false`.

### Plane-only and perpendicular guards

After the face infos are built, a cone paired with any non-planar partner is refused
before any cut work (`boolean_solid_curved.hpp:1849-1854`):

`unsupported "boolean_solid_curved: cone operands only support plane cutters"`

This is tested on the operand faces, not only on overlapping pairs, so a disjoint
cone/torus or cone/cylinder pair is refused too. It covers cone/cone, cone/torus,
cone/cylinder and cone/sphere; curved/curved cutters have no analytic cone cut
circle, no apex fan and no chart contract in this slice.

Every cone/plane face pair is then checked with
`curved_boolean_cone_plane_cut_supported` (`boolean_solid_curved.hpp:1378-1397`):

- perpendicular within `kCurvedBooleanConePerpendicularSlack = 1e-9` is accepted;
- any other orientation is accepted only when the plane misses a conservative
  bounding sphere of the cone: centre `origin + axis*h/3`, radius
  `hypot(radius, 2*h/3)`, where `h = radius/tan(half_angle)`. The sphere contains
  the base disk, not just the lateral surface.

A plane that is disjoint from the cone but still reaches that sphere is refused.
This is deliberately conservative: for `radius = 1`, `height = 2`, the sphere has
radius `5/3`, so the parallel plane `x = 1.5` is refused even though it misses the
cone (whose maximum radius is 1). The orientation is harmless only when the plane
provably cuts nothing, which avoids producing half-turn rim arcs or an ellipse the
splitter cannot represent. The refusal is

`unsupported "boolean_solid_curved: cone cutters must be perpendicular to the cone axis; oblique and axis-containing planes are deferred to M3c-3c-4b3"`

with the cone face and the plane face in the error context.

### Apex fan volume

`curved_solid_volume` (`include/geom/brep/detail/curved_volume.hpp:334-340`,
`458-497`) now recognizes a `Cone3` and the apex row `v = radius/tan(half_angle)`
(exact equality, as the chart side is computed). When a ring point lies on that row,
the ring's closing chord cuts the chart corner between the two meridians meeting at
the apex. The fan inserts the two missing corners at the neighbouring azimuths,
keeping the apex-row `v`, before triangulating:

```
if (previous.v != point.v && |previous.u - point.u| > 1e-12)
  insert {previous.u, point.v}
insert point
if (next.v != point.v && |next.u - point.u| > 1e-12)
  insert {next.u, point.v}
```

Every inserted corner evaluates to the same 3D apex, so no surface is invented. The
existing zero-area triangle skip (`twice_area == 0`, lines 522-529) removes the
duplicate apex triangles, so the fan closes the chart exactly once and does not
double count the apex. The same expansion is used for a loop with one apex copy (a
cap) and a loop with the cut end plus the back-seam copy (a band/wedge). The checked
values are the whole cone `2*pi/3` and the cone/plane fragments `pi/12` (upper small
cone) and `7*pi/12` (lower frustum).

`trimmed_face_boundaries` needs no cone guard any more: `trimmed_collapsed_v`
recognizes the same apex row and `trimmed_expand_collapsed_side`
(`include/geom/brep/trimmed_domain.hpp:194-249`, invoked at line 385) inserts the
missing corners with the arc references preserved. `mass_properties` sends a whole
cone to `mass_cone_domain` as before and a trimmed cone to the trimmed parametric
path (`mass_properties.hpp:826-833`); the cone-face refusals are gone.

### Apex weld

The two parameter copies of the apex are still present in `loops_uv`; the weld does
not merge them there. `assemble_curved_solid` Step 3c
(`include/geom/brep/assemble_curved.hpp:611-626`) drops the zero-length assembled
pieces after the coincident 3D points have welded to one vertex:

```
erase piece when piece.v0 == piece.v1
  and (!piece.arc || |piece.sweep| < 2*pi - 1e-12)
```

A straight zero-length piece (the apex-to-apex segment) is dropped. A closed rim
that starts and ends at the same vertex carries a full turn and is kept. The
neighbouring coedges therefore meet at the welded apex with no `v->v` edge, while
the per-face uv loops keep both copies for the pcurves and the fan. The two halves of
a cone cut by a plane now assemble into one closed shell.

### Refusals

| input | result |
| --- | --- |
| cone/cone, cone/torus, cone/cylinder, cone/sphere (including a disjoint pair) | `unsupported "boolean_solid_curved: cone operands only support plane cutters"` |
| oblique or axis-containing plane that reaches the cone's bounding sphere | `unsupported "boolean_solid_curved: cone cutters must be perpendicular to the cone axis; oblique and axis-containing planes are deferred to M3c-3c-4b3"` |
| trimmed cone fragment (not one loop over the full chart) | `unsupported "boolean_solid_curved: cone operands must be whole cones; trimmed cone fragments are not supported"` |
| trimmed cone tessellation | deferred to M3c-3c-4b3; `tessellate_cone_face` still fans the whole base ring to one apex vertex, which does not match a fragment's apex row |

The first row covers every non-planar partner, so the cone gate message and the
plane-only guard stay distinct: the gate is about the operand's face kind, the
perpendicular guard is about the cut plane's orientation.

### Deferred

Trimmed cone fragments, cone/cone and cone/curved cutters, oblique and axis-containing
planes (which need half-turn rim arcs or an exact ellipse), trimmed cone tessellation,
and NURBS operands stay outside this slice. The whole-cone restriction is also why a
cone that is one face of a larger B-rep is refused rather than re-charted.

## Evidence

Fixture: `make_cone({0,0,0}, {0,0,1}, radius = 1, height = 2, tolerance = 1e-12)`,
boolean tolerance `1e-9`. The plane `z = 1` gives the upper small cone (base radius
`0.5`, height `1`) and the lower frustum; the box is `{-2,-2,1}..{2,2,4}`, and the
union uses the low box `{-2,-2,1}..{2,2,3}`.

| operation | measured volume | analytic |
| --- | ---: | ---: |
| intersection | `0.2617993877991494` | `pi/12` |
| difference | `1.8325957145940461` | `7*pi/12` |
| union | `33.83259571459405` | `32 + 7*pi/12` |

QA-22 (`tests/unit/brep_boolean_curved_test.cpp`) runs the three operations once each
through `check_cone_boolean`, which asserts:

- `validate_topology(result).ok()`;
- exactly one shell and `shell_is_closed_by_uses` (every edge used by exactly two
  coedges);
- `|volume - analytic| <= 1e-12` absolute;
- `volume > 0` (with the winding-derived volume this is the outward orientation);
- a second call reproduces the volume bit for bit.

The focused run reports 1 case / 54 assertions / 0 failed. Full suites: Debug and
Release each report 1008 cases / 670,604 assertions, 0 failed, zero warnings under
`/W4 /WX /utf-8` with `GEOM_WERROR=ON`. The 4b1 baseline was 1001 / 670,383, and
F1's whole-cone refusal test brought it to 1008 / 670,604 before this task; QA-22
renames the existing three-operation case, so it adds no case or assertion.

## Consequences

- `boolean_solid_curved` now accepts a whole cone cut by planes, with exact analytic
  rims and a positive, closed, deterministic result. Sphere, torus, cylinder and
  plane paths keep their existing results; no result or fixture in that corpus
  changed.
- The cone surface gate is open only through the guards above; a non-planar partner
  and a non-perpendicular plane keep their own stable messages, and the old
  `"cone operands need closed-cut and apex topology"` message no longer fires for a
  `Cone3`.
- `curved_solid_volume` and `trimmed_face_boundaries` no longer refuse cone faces;
  the apex fan closes the chart once and the expanded trimmed side keeps the apex
  corners. `mass_properties` follows the same split.
- The apex weld drops only zero-length straight pieces and keeps full-turn rims, so
  the seam-only chart's two uv apex copies still reach the pcurves while the
  assembled shell has one apex vertex.
- Trimmed cone tessellation is still the 4b3 item: the whole-ring fan of
  `tessellate_cone_face` is not the sampling of a split fragment. General
  curved/curved cone booleans and oblique/axis-containing planes remain deferred.
