# ADR-0044: Rounded box as direct spherical-corner fragments

- Status: accepted (M4c-3c)
- Date: 2026-11-28
- Related: architecture §8.9, §16.4; `docs/benchmarks/m4c3c-rounded-box.md`; plan
  `docs/superpowers/plans/2026-11-27-m4c3c-rounded-box.md`; ADR-0042 for the
  single box edge family, ADR-0043 for convex prism vertical edges and ADR-0041
  for chamfer/shell

## Context

§8.9 lists rolling-ball fillet next to chamfer and shell. ADR-0042 shipped the
first exact slice: 1–4 parallel edges of an axis-aligned box, each a quarter
cylinder between two cap-plane arcs. ADR-0043 shipped the vertical edges of a
convex prism. Both deliberately left the vertex blend out: where filleted edges
meet at a corner, the rolling-ball surface is a spherical octant, and neither
slice needed one.

Selecting all twelve edges of a box is the smallest closed case that cannot be
done with cylinders alone. Every edge is filleted, so all eight corners are
three-way junctions, and each junction must be both watertight against its three
cylinders and exactly the rolling-ball surface. The question was whether to
build that rounded box directly now, or to wait for the generic edge-chain/SSI
blend engine, or to cut the corners with booleans.

## Decision

M4c-3c extends `fillet_edges(const Brep&, span<const EdgeId>, const
FilletOptions&)` with one more exact path: when the deduplicated selection is
exactly the twelve live edges of an axis-aligned box, build the twenty-six
fragments directly and assemble them once.

- **Recognition.** After `detail::fillet_box_bounds` succeeds, the selection is
  deduplicated by `EdgeId` and must contain exactly twelve unique ids, all
  alive, whose endpoints match the box corners within
  `max(64 * tolerance, 1e-12 * scale)`. The matched set must equal the box's
  twelve edges, in any order; duplicates are allowed and collapse. This check
  runs before `fillet_box_selection`, whose shared-vertex rejection would
  otherwise fire on the twelve-edge set (all twelve edges touch at the corners).
- **Six inset planar rectangles.** Each box face stays on its original plane
  `x_i = lo_i` or `hi_i`, trimmed to the rectangle whose other two coordinates
  lie in `[lo + r, hi - r]`, wound counter-clockwise about the outward normal.
  No new plane is introduced and no tolerance is used: the inset is exactly
  `r`, taken from the caller's radius.
- **Twelve quarter cylinders.** Each box edge becomes an exact `Cylinder3` of
  radius `r`, axis parallel to that edge at the two inset coordinates, spanning
  the edge's inner interval `[lo_f + r, hi_f - r]`. Its cross-section arc runs
  from the tangency direction of one adjacent face to the other, so the patch is
  tangent to both planes along lines and is bounded at each end by a quarter
  `CircleCurve3` arc where it meets the corner sphere. The cylinder axis passes
  through both adjacent sphere centres, so the sphere/cylinder intersection is
  exactly a great circle and needs no SSI solve.
- **Eight spherical octants.** Each box corner becomes an exact `Sphere3` of
  radius `r` centred at the corner's inset point `(lo|hi)_i ± r`; its polar axis
  points from the centre to the box corner, and the patch is the octant facing
  that corner. The three boundary curves are great-circle quarter arcs lying in
  the planes through the sphere centre at the three inset coordinates; each arc
  is a `CircleCurve3` with `t0 = 0, t1 = pi/2`. In the corner frame with
  coordinates `(u, v, w)` measured from the centre toward the corner, the
  boundary points are `(0, r, r)`, `(r, 0, r)` and `(r, r, 0)`, and the outward
  loop visits them in that consistent order. The face's UV loop is the spherical
  chart projected with `face_uv_project_sphere`, carrying `previous_u` around
  the loop without wrapping it back by a full turn: the octant stays a single
  chart patch and `face_loops_uv`/`mass_properties` accept it.
- **Assembly.** The twenty-six `CurvedAssembledFace` fragments go through
  `assemble_curved_solid` once and are returned only after topology validation
  is clean. Every boundary is exact geometry: twenty-four straight degree-1
  polyline edges (each an exact line, represented as a degree-1 NURBS) and
  twenty-four `CircleCurve3` great-circle arcs. There is no Boolean, no SSI, no
  faceted cylinder or sphere, and no sampled fit. The result is
  `FilletResult{brep, 12u}`.
- **Topology is closed-form.** The twenty-four straight tangent lines (two per
  cylinder, one shared with each adjacent plane) and the twenty-four
  great-circle arcs (two per cylinder end, three per octant, each shared) give
  `V = 24`, `E = 48`, `F = 26` — six planes, twelve cylinders, eight spheres —
  `chi = 2`, one solid and one closed shell.

## Volume

Write `A = a - 2r`, `B = b - 2r`, `C = c - 2r` for a box `a x b x c`. The
rounded solid is the inset box `ABC` plus twelve quarter-cylinder edge shells
plus eight spherical octants; equivalently, start from `abc`, remove at each
edge the square prism minus its quarter cylinder, `r^2 L (1 - pi/4)` with each
of the four `A`, `B`, `C` lengths shared by four edges, and remove at each
corner the corner cube minus its octant, `r^3 (1 - pi/6)`. Both routes simplify
to the Steiner form

`V = ABC + 2(AB + BC + CA) r + pi (A + B + C) r^2 + 4 pi r^3 / 3`.

For the benchmark box `a = b = c = 100`, `r = 1` this is
`999743.8170303602`. The analytic formula is the volume authority; the
benchmark compares the integrated `detail::curved_solid_volume` against it, and
the unit tests pin the same formula for `2 x 3 x 4` at `r = 0.25`
(`23.53807140332`).

## Scope

- Axis-aligned box only; all twelve edges selected, in any order, with
  duplicates allowed after deduplication.
- One radius, uniform across the twelve edges, finite and strictly
  `0 < r < 0.5 * min(a, b, c)`.
- The input box is not modified; two calls produce bit-identical results.

## Rejection contract

Every input or geometry rejection happens before assembly and leaves the input unchanged:

- Empty selection: `unsupported "fillet_edges: at least one edge is required"`
  (checked before the options, so an empty selection with default options is not
  reported as a zero radius).
- Non-finite / non-positive options: the existing `non_finite`, `degenerate`
  and `tolerance_violation` messages are unchanged.
- Radius at or above half the smallest dimension, including exact equality:
  `unsupported "fillet_edges: radius must be less than half the smallest box
  dimension"`. Strict inequality is the contract: at equality the eight inset
  rectangles collapse to lines and the solid would be degenerate.
- A dead edge anywhere in the selection, complete or not:
  `index_out_of_range "fillet_edges: edge does not belong to the solid"`,
  raised before the rounded/box/prism dispatch so the diagnostic does not depend
  on which path would have matched.
- Non-box input: the existing shape recognition rejects it (`unsupported
  "fillet_edges: only a box or a convex prism is supported"`); a prism-shaped
  selection that is not a box never enters this path.
- Partial or mixed selections are **not** the rounded path. They fall through
  to the pre-existing family and prism paths unchanged, so eleven of twelve
  edges fails there with the stable shared-vertex / parallel / not-a-corner
  diagnostics rather than producing a partially rounded box. Dead edges are the
  exception that is pulled forward to `index_out_of_range`.
- An assembly failure that the guards did not predict:
  `invalid_topology "fillet_edges: fillet could not be assembled"`, never an
  anonymous assembler message. The rounded path has this single tail; the
  box-family and prism paths may additionally emit `"fillet_edges: fillet
  produced an invalid solid"` for their own validation step.

## Alternatives rejected

- **Boolean cuts** (`box - (corner wedge - sphere)`, or subtracting the eight
  corner cubes and adding back the octants): needs exact curved intersections,
  inherits splitter face-order sensitivity, and is numerically worse than the
  analytic fragments. ADR-0041/0042/0043 reject the same route for the same
  reasons.
- **Faceted cylinders or spheres** (a chain of planar strips, or a polyhedral
  corner): deviation is unbounded without a curvature bound, adds a sampling
  policy and a tolerance story, and contradicts the project rule that a faceted
  substitute must not be called a fillet.
- **Per-corner blends only** (fillet the edges, then patch the three-way
  junction with a separate blend surface): the junction is already exactly a
  spherical octant; a separate blend either leaves a gap, doubles the material,
  or needs the SSI machinery this case avoids. Treating each edge independently
  cannot produce a closed shell at a corner.
- **Generic sweep / edge-chain / SSI engine now.** That engine must cover
  non-parallel, concave and variable-radius edges, and it still owes curved
  intersection curves and general trimming. The all-edges box is exactly
  solvable today with plane, cylinder and sphere primitives, so shipping it does
  not block or pre-empt the general engine.
- **Sampled/fitted corner surface.** The surface is exactly a sphere octant, so
  a fit only adds error and a sampling policy.
- **Variable radius or a per-edge radius now.** Out of scope; the rolling-ball
  construction with a varying radius is a different problem (M4c-3d).

## Known limits

- **Near-limit radius is thin but valid.** Just below `0.5 * min`, the six
  inset rectangles shrink toward zero width and the corner octants nearly
  touch, but the solid remains a valid closed B-rep: the probes at
  `r = 0.5 * min - 1e-6` and even `-1e-9` assemble twenty-six faces, twenty-four
  vertices and pass `validate_topology`, with the Steiner volume matching the
  integrated volume. The guard rejects equality only, so this is accepted as a
  geometrically thin (not topologically degenerate) solid; no sliver face is
  silently dropped.
- **Partial corners, mixed families and variable radius** are M4c-3d, together
  with non-box operands and general rolling-ball blends. Selecting fewer than
  twelve edges is the existing ADR-0042/0043 behaviour, unchanged.
- **`mass_properties` gate.** The spherical octants are not full rectangles, so
  `mass_properties` admits a `Sphere3` face only when
  `detail::tessellate_face_is_sphere_patch` proves it is a single-loop patch
  with exactly three non-closed `CircleCurve3` arcs. The first Task-2 version
  exempted every `Sphere3` face, which broke the pre-existing rejection of an
  untrimmed curved face; the gate was narrowed to the patch predicate
  (`7f9df10`). Cylinder, cone and torus faces keep their original gate.
- **Tessellation is separate evidence.** The corner octants are tessellated by
  the dedicated sphere-patch tessellator; watertightness (`boundary_edges == 0`)
  is asserted by unit tests, not by this benchmark, which does not tessellate.
- **Benchmark scope.** The benchmark measures recognition plus the
  twenty-six-fragment construction and assembly for a `100` cube at `r = 1`; it
  is not a general performance claim and is not part of ctest.

## Consequences

- Constant-radius rounding of all twelve edges of a supported axis-aligned box
  is exact, deterministic and returned as a valid closed B-rep with the Steiner
  volume and `24/48/26` topology. The benchmark measures
  `29.784 / 31.382 / 33.426 ms` (min/median/max, one representative five-run
  batch) and volume `999743.817030359060` against the analytic
  `999743.8170303602` (relative error about `1.2e-15`).
- Unit tests pin construction, all eight sphere centres (one per corner sign),
  exact bounds and surface inventory, the `2 x 3 x 4 / r = 0.25` volume and
  centroid, watertight tessellation, radius rejection at and above `0.5 * min`,
  the near-limit success, duplicate and dead selections, exact messages, and
  rounded-path determinism and input immutability. The Debug
  `rounded box*, sphere patch*` lane is `11 cases / 1,862 assertions`.
- The path is additive: the box-family and prism fillet paths and their messages
  are unchanged, and a selection that is not exactly the twelve edges is still
  handled by them.
- Deferred to M4c-3d: partial corner selections, mixed edge families, variable
  radius, concave and non-box operands, and the general rolling-ball blend. No
  document claims general fillet support.
