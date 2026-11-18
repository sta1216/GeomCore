# ADR-0033: Mesh and mass properties for rectilinear trimmed-face holes

Status: accepted for M3c-3c-4a (2026-11-18)

## Context

M3c-1 tessellated a curved analytic face trimmed by a box. Its boundary is rectilinear
in the face's uv chart: cap circles and window rims are lines of constant `v`, box
walls are lines of constant `u`. The tessellator decomposed that rectilinear domain
into boundary-aligned strips and split a periodic outer ring into two pieces at an
artificial cut, so the periodic weld is never crossed by a diagonal. Mass properties
integrated the same trimmed face, but only when its boundary had a single loop.

M3c-2 deferred the hole case. When the cut window does not cross the seam, the trimmed
face carries an outer ring plus a hole ring: the window is open, yet the face's
boundary is two loops. The mesher refused it with
`tessellate: trimmed curved face has a hole`; the mass integrator refused it with
`mass_properties: trimmed curved face has a hole`. Both refusals were gated before any
geometry was examined, so the limitation applied even to a perfectly rectilinear
window whose volume has a closed form: the acceptance fixture is a cylinder `r = 1`,
`h = 2` cut by the box `{-0.5,0.5,0.5}..{0.5,2,1.5}`, i.e. `2*pi - 0.45661147746179975`.

M3c-3c-4a lifts the rectilinear, non-seam-crossing case. Two properties made that
tractable without a new triangulator:

- the strip grid is built from the boundary samples themselves, so a rectilinear hole
  only adds rows and columns to the same grid; and
- the even-odd crossing rule already returns several inside runs per band, so a hole
  splits a band into the intervals on either side of it.

The mesh has one extra constraint a plain uv triangulation does not have: on a
periodic face the two pieces meet on a shared seam wall and on the artificial cut
wall, and those wall vertices are shared position indices. A hole that crosses or
touches the artificial cut has no consistent piece assignment, so it must be refused
rather than clipped per piece.

## Decision

- **Strip decomposition with holes.** `tessellate_trimmed_strip_triangles_impl`
  (`include/geom/brep/tessellate.hpp:1101`) takes `rings = {outer, holes...}` and
  derives `xs`, `ys`, horizontals, verticals, chains, wall samples and the vertex map
  from every ring. A hole sample keeps the hole ring's weld position, so the zipper
  emits the same boundary vertices the neighbouring faces see. `tessellate_expand_trimmed_runs`
  and `triangulate_trimmed_domain` are reused unchanged; no new triangulation
  algorithm and no second fan pass are introduced.
- **Per-piece assignment.** For a periodic outer ring, `tessellate_trimmed_face`
  (`tessellate.hpp:1770`) splits the ring into its two pieces and assigns every hole
  to the piece whose bounds strictly contain every one of its samples. A hole sample
  on the artificial cut, or a hole not strictly inside one piece, returns
  `Errc::unsupported` with `tessellate: trimmed curved face hole crosses the periodic
  cut` (`tessellate.hpp:1836`). Hole row levels are passed to both pieces, so the
  shared seam and cut walls carry identical split vertices (no T-junction fan).
- **Rectilinearity is a precondition for every ring.** A non-rectilinear hole returns
  `tessellate: trimmed curved face hole boundary is not rectilinear`
  (`tessellate.hpp:1747`) before triangulation, independently of the deflection; a
  non-rectilinear outer ring keeps its existing message
  (`tessellate: trimmed curved face boundary is not rectilinear`). The check is
  `tessellate_ring_rectilinear` (`tessellate.hpp:1085`), shared with the mass path so
  the two rejections agree.
- **Mass integrates all boundaries.** `mass_add_trimmed_face`
  (`include/geom/brep/mass_properties.hpp`) drops the single-boundary gate, requires
  every boundary to be rectilinear, simplifies each one, and hands outer plus holes to
  `triangulate_trimmed_domain` in one call. A non-rectilinear boundary returns
  `mass_properties: trimmed curved face boundary is not rectilinear`
  (`mass_properties.hpp:447`). The mass chart has no artificial cut, so it integrates
  the outer ring and every hole directly.
- **Hole-free output is unchanged.** The four-argument tessellation overload delegates
  with an empty hole span, so its output stays bitwise identical; the hole-free
  canonical hash and count tests still pass.

## Evidence

Fixture: cylinder `r = 1`, `h = 2` at the origin cut by the box
`{-0.5,0.5,0.5}..{0.5,2,1.5}`, tolerance `1e-9`; the result has 8 faces and passes
`validate_topology`. MSVC 19.44 (`/W4 /WX /utf-8`, zero warnings). The exact lateral
uv area is `4*pi - pi/3 = 11.519173063162574`; the exact analytic volume is
`2*pi - 0.45661147746179975 = 5.826573829717787`.

| quantity | measured | gate |
| --- | --- | --- |
| fine mesh `{deflection 1e-3, 0.35}` | 1236 triangles, 620 positions | inside the `[950, 1450]` window |
| fine mesh boundary / non-manifold edges | 0 / 0 | watertight |
| fine mesh lateral uv area | `11.519173063162778` | within `1e-6` of exact (test asserts `epsilon(1e-12)` relative) |
| fine mesh volume | `5.8238533084266857` | inside `[5.8200, 5.8266]` acceptance band (no unit test asserts the mesh volume) |
| coarse mesh `{1e-2, 0.35}` | 452 triangles, 228 positions, boundary 0, non-manifold 0 | watertight |
| coarse mesh lateral uv area | `11.519173063162546` | within `1e-6` of exact (test asserts `epsilon(1e-12)` relative) |
| coarse mesh volume | `5.8001219438259604` | below `[5.8200, 5.8266]` (coarser chords; not a gate) |
| mass `{1e-9, 24, 1e-5}` volume | `5.8265738271299785` | `4.44e-10` relative to analytic, gate `1e-9` |
| mass centroid | `(-2.7e-11, -0.05720914955052997, 1.0)` | centroid moves away from the `+y` hole |

The mass volume and the coarse mesh volume are recorded as measured, not as ideal
values. The mass volume has the Task 2 `1e-9` relative gate; the fine mesh volume has
no unit assertion and is checked against the milestone acceptance band recorded here;
the uv area is asserted by the Task 1 cases at `epsilon(1e-12)` relative against the
exact `4*pi - pi/3`, well inside the `1e-6` recording tolerance. The QA case weighs the
same fixture in the same test as the mesh, with the looser `1e-6` relative volume bound
requested for the acceptance test (`tests/unit/tessellate_trimmed_test.cpp`, `QA-20`).

Refusals measured on constructible inputs:

| input | result |
| --- | --- |
| oblique trim, plane normal `(0.3, 0, 1)` | tessellate `unsupported "tessellate: trimmed curved face boundary is not rectilinear"`; mass `unsupported "mass_properties: trimmed curved face boundary is not rectilinear"` |
| window straddling the artificial cut, box `{-2,-0.5,0.5}..{-0.5,0.5,1.5}` | tessellate `unsupported "tessellate: trimmed curved face hole crosses the periodic cut"`; mass succeeds with volume `5.8265738271295096` |

The second row exposes an intended asymmetry: mass does not split the periodic ring, so
it can weigh a hole that the mesher refuses per piece; the same analytic volume is
returned. The non-rectilinear *hole* message has no constructible fixture yet. A
perpendicular cylindrical cutter (base `(0,0.5,1.0)`, axis `(0,1,0)`, radius `0.3`,
height `1.5`) is refused upstream by the boolean with
`boolean_solid_curved: flipped curved fragments are not supported yet`, before any hole
reaches the mesher. That refusal is therefore pinned by the pre-check code and review,
not by a test; the benchmark note marks it unmeasured.

## Consequences

- A curved face whose outer and hole rings are rectilinear, with the hole strictly
  inside the outer ring, now meshes watertight and weighs; on a periodic-`u` face the
  hole must lie strictly inside one split piece and no sample may touch the artificial
  cut. This is the supported class.
- Still refused, unchanged: non-rectilinear hole rings; holes crossing or touching the
  artificial cut; stored holes on synthetic faces
  (`split_face_curved: stored holes are not supported yet`, ADR-0032); general,
  twisted-seam, and NURBS-hole cases deferred to M3c-3c-4; and boolean results whose
  curved outline builder cannot produce a clean hole at all (the upstream
  near-coincident-point limits recorded in ADR-0027 remain).
- The mesh contract is stricter than the mass contract on a periodic chart: meshing
  refuses a hole straddling the artificial cut, while mass integrates it. That is
  intended for M3c-3c-4a; closing the gap would require a mass-side piece split or a
  chart re-origin, neither of which this milestone needs.
- Multiple holes are not claimed. The strip grid and the per-piece assignment are
  ring-count agnostic, and a same-piece two-hole face was observed to mesh, but only
  the one-hole case is tested. A hole whose sampled ring unwraps outside its split
  piece (for example a second hole on the far side of the chart) is refused with
  `tessellate: trimmed curved face hole crosses the periodic cut`; normalising each
  hole ring into its piece is a follow-up. The message is conservative: it also
  covers a hole that is geometrically inside one piece but was unwrapped outside
  it.
- Hole-free meshes and masses are bitwise unchanged.
- Deferred to M3c-3c-4: multi-cycle holes, cone apex topology, general twisted seam
  crossings, and meshing of trimmed periodic-`v` faces.

## Alternatives rejected

- **Fan post-pass.** `tessellate_trimmed_polygon_triangles` plus
  `tessellate_expand_trimmed_runs` already exist and fan boundary sub-segments, but
  applied to the whole holed rectangle the ear clipper spans the window with long
  chords: the 3D area inflates and the mesh folds, which is precisely why the strip
  decomposition exists. The fan expansion also only splits boundary runs after
  triangulation; it does not assign holes to periodic pieces.
- **Full constrained triangulation (CDT).** Unnecessary for this class: every segment
  of a rectilinear hole is already a row or column of the strip grid, so the grid
  contains the hole exactly. A CDT would add a new triangulator, new numeric edge
  cases, and a wider refusal surface for no gain.
- **Keep the M3c-2 refusal.** Rejected because the fixture has an exact analytic
  oracle, the strip grid covers rectilinear holes without a new algorithm, and both
  the mesh and the mass path were already one gate away from accepting it.
