# ADR-0034: Cone chart, apex projection and open cuts on domain sides

Status: accepted for M3c-3c-4b1 (2026-11-19)

## Context

M3c-3c-2b (ADR-0032) generalized the synthetic-domain splitter for torus faces and
kept the sphere on its historical seam-only wrapper. A cone lateral face was still
unusable even though its chart has the same periodic-`u` shape as a cylinder:

- `face_uv_project_cone` refused every point at or beyond the apex (`rho <= 1e-12 *
  radius`), so the stored slit loop could not be projected and every downstream
  consumer inherited the error.
- `trim_face_domain` knew the sphere and torus charts only; there was no cone chart.
- `split_face_curved` had no cone dispatch, so a canonical cone fell through to the
  generic clip path with a stored slit loop whose area is not zero.
- `trim_face_split_domain` paired open cut pieces seam-to-seam only. A meridian from
  the base rim to the apex ends on the base boundary and on the collapsed apex side,
  so it was refused with `"split_face_curved: open cuts on a seam-only face"`.

The canonical fixture is a cone at the origin with `radius = 1`, height `2`, axis
`+z` and tolerance `1e-12`; its lateral face is cut along `u = 0`. The stored slit
loop has area `pi * v_apex = 2*pi`, half of the chart rectangle's `4*pi`.

## Decision

- **Apex projection.** `face_uv_project_cone` (`include/geom/brep/face_uv.hpp:140`)
  accepts a point inside the relative rounding band `apex_tolerance = 1e-9 *
  max(1, v_apex)` when `v_apex > apex_tolerance`, and returns
  `{has_previous_u ? previous_u : 0, v_apex}`. The caller's `previous_u` is threaded
  through `face_uv_project` (`face_uv.hpp:241`), so the apex does not introduce a
  `u` jump into the ring. A point strictly below the apex still uses
  `rho = radius - v*tan(half_angle)`; a point at or beyond the apex outside the
  band still returns `"face_loops_uv: cone projection is undefined at or beyond the
  apex"`. No new absolute epsilon is introduced: the band scales with `v_apex`.
- **Cone chart.** `trim_face_domain` (`include/geom/brep/trim_domain.hpp:93`) maps a
  `Cone3` to `bounds = {0, 2*pi, 0, radius/tan(half_angle)}`, `periodic_u = true`,
  `periodic_v = false`, and `sides = {kSeamCopy, kSeamCopy, kRealBoundary,
  kCollapsed}`. `u = 0/2*pi` is the seam, `v = 0` is the base boundary and
  `v = v_apex` is collapsed.
- **Synthetic detection.** A cone face is synthetic when it stores exactly one loop,
  `||signed_area2(loop)| - 0.5 * chart_area| <= 1e-9 * chart_area`, every loop point
  lies within the parameter tolerance of `v_lo` or `v_hi`, and exactly one point lies
  on the collapsed `v_hi` row (`include/geom/brep/trim_face.hpp:2147-2167`). A
  synthetic face with cuts is routed to `trim_face_split_domain` with empty holes; a
  synthetic face with no cuts falls through unchanged. Any other cone face returns
  `"split_face_curved: trimmed cone faces are not supported yet"` with
  `(face.index(), 0)`; it is never reopened as the full chart.
- **Open cuts on domain sides.** After seam pairing fails, an open piece is accepted
  when both ends are on domain sides, neither end is on a periodic seam, and the two
  ends are on non-seam sides (`trim_face.hpp:1460-1513`). Each end is then snapped
  exactly to `u_lo`/`u_hi` or `v_lo`/`v_hi`, so the synthetic side segment is split
  at the cut endpoint and the fragment closes. A piece that leaves through a
  periodic seam and returns elsewhere returns
  `"split_face_curved: open cut crosses the seam"`; a piece that is not
  boundary-to-boundary, or that has a seam end on both sides, keeps
  `"split_face_curved: open cuts on a seam-only face"`.
- **Representation.** Fragments still carry `loops` and `loops_uv` of equal length.
  The collapsed apex side becomes apex-row points only: the narrow wedge carries the
  cut's apex end once, while the wide remainder carries the cut end and the
  back-seam copy at `u = 2*pi`. Every apex-row point evaluates to the same 3D apex
  `(0,0,2)` within `1e-12`; the plane is not welded or fanned across fragments.

## Refusal matrix

| input | result |
| --- | --- |
| cut chain ends strictly inside the chart | `unsupported "trim_face: cut chain ends inside the parameter domain"` (`trim_domain.hpp:431`) |
| cut chain lies on a domain side | `unsupported "split_face_curved: cut lies on the parameter domain boundary"` (`trim_face.hpp:1046`, `trim_face.hpp:1399`) |
| piece leaves through a seam and returns elsewhere | `unsupported "split_face_curved: open cut crosses the seam"` (`trim_face.hpp:1485`) |
| open piece is not boundary-to-boundary, or has a seam end | `unsupported "split_face_curved: open cuts on a seam-only face"` (`trim_face.hpp:1490`) |
| trimmed cone face (stored loop not the full half-chart) | `unsupported "split_face_curved: trimmed cone faces are not supported yet"`, context `(face.index(), 0)` (`trim_face.hpp:2172`) |
| cone volume | `unsupported "curved_solid_volume: cone faces are not supported yet"` (`include/geom/brep/detail/curved_volume.hpp:337`) |
| cone trimmed boundaries | `unsupported "trimmed_face_boundaries: cone faces are not supported yet"` (`include/geom/brep/trimmed_domain.hpp:348`) |
| cone boolean operands | `unsupported "boolean_solid_curved: cone operands need closed-cut and apex topology; booleans are deferred to M3c-3"` (`include/geom/brep/boolean_solid_curved.hpp:187`) |

## Evidence

Fixture: canonical cone `radius = 1`, height `2`, solid tolerance `1e-12`, split tolerance `1e-9`; MSVC 19.44
(`/W4 /WX /utf-8`, zero warnings). QA-21
(`tests/unit/brep_trim_face_test.cpp`) pins the public
`split_face_curved(cone, lateral, {{(0.7,0),(0.7,2)}}, 1e-9)`:

| quantity | measured |
| --- | --- |
| fragments | 2 |
| uv areas | `0.7` and `11.166370614359172` (`2*(2*pi - 0.7)`), test tolerance `1e-9` |
| apex-row point counts | 1 and 2 |
| apex evaluation | every `v == 2` point evaluates to `(0,0,2)` within `1e-12` |
| determinism | second run repeats every uv and 3D coordinate bitwise |
| focused run | 1 case / 32 assertions / 0 failed |

The Task 2 closed-cut case remains the companion pin: a horizontal cut at `v = 1`
gives a `pi` cap with one apex point and a `2*pi` band with none. Task 4 adds no
production behaviour; QA-21 intentionally repeats the Task 3 meridian acceptance as
the end-to-end pin.

The RED for this contract was measured on the parent commit `b9efed7`, using that
commit's `trim_face.hpp`/`trim_domain.hpp` over the current tree: the parent has
the cone chart and dispatch but not the open-side acceptance, so the same public
meridian call returns `unsupported "split_face_curved: open cuts on a seam-only
face"` with context `1/0`.

Full suites: Debug and Release each report 1001 cases / 670,383 assertions, 0 failed,
exit 0; QA-21 adds 1 case / 32 assertions over the 1000 / 670,351 baseline.

## Consequences

- A canonical cone lateral face now splits closed cuts and open cuts whose ends land
  on domain sides. The apex is one uv point per collapsed side; this does not create
  shared topology between fragments.
- **Assembly of the two halves is deferred.** The two fragments each carry an apex
  end, but assembling them fails with `"assemble_curved_solid: edge is used once"` /
  `"assemble_solid: edge is used once"`, and the collapsed side would produce a
  zero-length apex edge. No cone assembly, mass or volume path is enabled.
- **The boolean gate stays closed.** `boolean_solid_curved` still rejects cone
  operands with the message above, and no boolean operation changes result.
- **Trimmed cone faces are still refused.** Only the canonical full half-chart face
  is split; a trimmed cone face keeps the dispatch refusal. `trimmed_face_boundaries`
  and `curved_solid_volume` still refuse cone faces.
- **Apex tessellation and sampling stay mismatched.** The canonical cone mesh and
  mass properties are unchanged. `tessellate_cone_face` fans the whole base ring to
  the single apex vertex; a split fragment's apex row plus cut end does not match
  that sampling, so trimmed-cone tessellation and apex-fan volume remain deferred.
- General twisted seam crossings, other nappes, NURBS operands and trimmed cone
  meshing remain outside this milestone.
