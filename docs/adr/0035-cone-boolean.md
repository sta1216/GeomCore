# ADR-0035: Cone operands in the curved boolean

Status: accepted for M3c-3c-4b2 (2026-11-19); the trimmed-cone tessellation
limitation is narrowed by ADR-0036 (frusta now tessellate; apex-bearing fragments
keep the canonical fan). Axis-containing cone cuts are supported since ADR-0038
(M3c-3c-4b3c-2); oblique planes remain deferred to M3c-3c-4b3.

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
- every cone/plane pair is perpendicular to the cone axis (within `1e-9`), contains
  the cone axis (within `1e-9`, with the cone origin on the plane within the
  boolean tolerance), or is provably disjoint (see the conservative predicate below);
- the two parameter copies (`u = 0` and `u = 2*pi`) of a cut-circle seam point stay
  distinct in double precision. A placement where they round to the same point is
  refused deterministically (see the refusals table); re-centring the operands to
  recover it is deferred to M3c-3c-4b3.

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
`curved_boolean_cone_plane_cut_supported` (`boolean_solid_curved.hpp:1445-1466`):

- perpendicular within `kCurvedBooleanConePerpendicularSlack = 1e-9` is accepted;
- containing the cone axis (normal perpendicular to the axis within `1e-9`, cone
  origin on the plane within the boolean tolerance) is accepted (ADR-0038);
- any other orientation is accepted only when the plane misses a conservative
  bounding sphere of the cone: centre `origin + axis*h/3`, radius
  `hypot(radius, 2*h/3)`, where `h = radius/tan(half_angle)`. The sphere contains
  the base disk, not just the lateral surface.

A plane that is disjoint from the cone but still reaches that sphere is refused.
This is deliberately conservative: for `radius = 1`, `height = 2`, the sphere has
radius `5/3`, so the parallel plane `x = 1.5` is refused even though it misses the
cone (whose maximum radius is 1). The orientation is harmless only when the plane
provably cuts nothing, which avoids producing an ellipse the splitter cannot
represent (axis-containing cuts reuse the half-turn base rim, which the assembler
now supports). The refusal is

`unsupported "boolean_solid_curved: cone cutters must be perpendicular to or contain the cone axis; oblique planes are deferred to M3c-3c-4b3"`

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
| oblique plane that reaches the cone's bounding sphere | `unsupported "boolean_solid_curved: cone cutters must be perpendicular to or contain the cone axis; oblique planes are deferred to M3c-3c-4b3"` |
| trimmed cone fragment (not one loop over the full chart) | `unsupported "boolean_solid_curved: cone operands must be whole cones; trimmed cone fragments are not supported"` |
| cone placement where the two seam copies of a cut-circle point round to the same double (observed with base centre `y = -2`, or `|y| >= 2.25`, and radius `0.5`; `y = +2` still passes) | `unsupported "boolean_solid_curved: cone seam copies coincide at this placement; cone booleans at this offset are deferred to M3c-3c-4b3"` |
| apex-bearing trimmed ring that reaches the trimmed path | `unsupported "tessellate: trimmed curved face boundary is not rectilinear"`; frusta tessellate since M3c-3c-4b3a (ADR-0036), and apex-bearing fragments keep the canonical `tessellate_cone_face` fan |

The first row covers every non-planar partner. The cone gate message therefore
concerns the partner's face kind (for example `cone/cone` or `cone/torus`); it no
longer fires for a `Cone3` operand of the expected shape. The orientation guard (perpendicular,
axis-containing, or a harmless disjoint plane) is separate and concerns only the cut
plane.

The whole-cone gate is about boolean operands: a frustum produced by a supported cut
tessellates through the trimmed path (ADR-0036), while apex-bearing fragments keep the
canonical fan.

### Deferred

Trimmed cone fragments as boolean operands, cone/cone and cone/curved cutters, oblique
planes (ellipse sections), and
NURBS operands stay outside this slice. Frusta produced by a supported cut now
tessellate (ADR-0036), while apex-bearing fragments keep the canonical whole-ring fan.
The whole-cone restriction is also why a cone that is one face of a larger B-rep is
refused rather than re-charted.

A translated whole cone subtracted from a provably disjoint box at `|y| >= 10` still
fails inside the pre-existing base-rim heal (`assemble_curved_solid: shared edge curve
mismatch`). That is not a regression from this slice; it is tracked with the
re-centring work in M3c-3c-4b3.

## Evidence

Fixture: `make_cone({0,0,0}, {0,0,1}, radius = 1, height = 2, tolerance = 1e-9)`,
boolean tolerance `1e-9`. The plane `z = 1` gives the upper small cone (base radius
`0.5`, height `1`) and the lower frustum; the box is `{-2,-2,1}..{2,2,4}`, and the
union uses the low box `{-2,-2,1}..{2,2,3}`.

| operation | computed volume | analytic |
| --- | ---: | ---: |
| intersection | `0.26179938779915246` | `pi/12` |
| difference | `1.8325957145940359` | `7*pi/12` |
| union | `33.832595714594021` | `32 + 7*pi/12` |

QA-22 (`tests/unit/brep_boolean_curved_test.cpp`) runs the three operations once each
through `check_cone_boolean`, which asserts:

- `validate_topology(result).ok()`;
- exactly one shell and `shell_is_closed_by_uses` (every edge used by exactly two
  coedges);
- `|volume - analytic| <= 1e-12` absolute;
- `volume > 0` (with the winding-derived volume this is the outward orientation);
- a second call reproduces the volume bit for bit.

The focused run reports 1 case / 54 assertions / 0 failed. Full suites: Debug and
Release each report 1009 cases / 670,614 assertions, 0 failed, zero warnings under
`/W4 /WX /utf-8` with `GEOM_WERROR=ON`. The 4b1 baseline was 1001 / 670,383; Tasks 1-3 and the F1 whole-cone refusal test
brought it to 1008 / 670,604 before this task (F1 alone: +1 case, +16 assertions);
the translated-cone refusal adds one case and ten assertions for 1009 / 670,614.
QA-22 renames the existing three-operation case, so it adds no case or assertion.

## Consequences

- `boolean_solid_curved` now accepts a whole cone cut by planes, with exact analytic
  rims and a positive, closed, deterministic result. Sphere, torus, cylinder and
  plane paths keep their existing results. One intentional boundary detail changed:
  the trimmed boundary of a sphere boolean face now prepends the chart corner
  `(-2*pi, -pi/2)`, making the ring rectilinear; no volume or boolean result changed,
  and no test pins that corner yet.
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
- Trimmed cone frusta now tessellate (M3c-3c-4b3a, ADR-0036): the existing strip
  zipper joins the two independently sampled rims, and the sag metric now measures
  the surface at the same azimuth and axial height. Apex-bearing fragments still use
  the whole-ring fan of `tessellate_cone_face`; no supported path produces an apex
ring in the trimmed path. General curved/curved cone booleans and oblique
  planes remain deferred.
