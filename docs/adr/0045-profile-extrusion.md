# ADR-0045: Profile extrusion as planar caps plus quadrilateral walls

- Status: accepted (M4d-1)
- Date: 2026-11-28
- Related: architecture §8.7, §16.4; `docs/benchmarks/m4d1-extrusion.md`; plan
  `docs/superpowers/plans/2026-11-28-m4d1-profile-extrusion.md`; ADR-0041 for
  the assembler and half-space contract; `include/geom/cgeom/triangulate2.hpp`
  for the simplicity certificate

## Context

§8.7 lists `extrude` (straight direction) among the sweep operators next to
revolve, loft and sweep, and asks every sweep to detect self-intersection before
building. M4d-1 is the first narrow slice: a planar polygon with holes, swept
along the plane normal by a positive height. It has to answer four questions.

- **What is the solid?** Caps are planar polygons with holes and walls are
  planar quads; the caps must stay polygon faces, not triangle soups, because
  the rest of the library treats a planar face as one polygon-with-holes loop
  set. At the same time the volume integrator consumes triangles.
- **How is a non-simple profile rejected before any face exists?** The project
  already has `triangulate_polygon2`, whose covered-area invariant and ear
  clipping see self-intersections and misplaced holes; the question was whether
  that is enough to serve as the certificate, or whether a separate geometric
  simplicity test is needed.
- **Which winding is the contract?** Input rings arrive in either winding; the
  solid needs outer counter-clockwise and holes clockwise about the extrusion
  normal, stable enough that volume and side normals agree.
- **Exact volume without faceting.** B-rep volume comes from triangulating the
  caps internally, so the caps' triangulation is where the volume is decided.

## Decision

`extrude_profile(const Plane3& plane, const Polygon2d& profile, double height,
double tolerance) -> Result<Brep>` builds the closed B-rep directly, with no
Boolean, no SSI and no faceted substitute.

- **Construction.** The profile's rings are mapped through the plane frame
  (`make_plane_frame`), the top ring points are `base + frame.normal * height`,
  and the face list is one bottom cap, one top cap, then one quad per edge of
  every ring. The bottom cap reverses every normalized ring; the top cap keeps
  the normalized order. Each side quad is `{b_i, b_next, t_next, t_i}` with
  outward normal `cross(b_next - b_i, normal)`, normalized after scaling the
  largest component out first so the length stays in `[1, sqrt(3)]` for extreme
  coordinates. A zero length normal is left to the assembler as a degenerate
  face. Rings are passed to `assemble_solid` in profile order: outer first, then
  holes.
- **Orientation normalization.** The outer ring ends counter-clockwise and every
  hole clockwise, so the material lies to the left of each ring; an input ring
  wound the other way is reversed. A clockwise outer ring is therefore accepted
  and normalized, which is the documented input contract (ruling R3 of the
  plan). A ring whose orientation is undefined (`ring_orientation2 == 0`, e.g.
  fewer than three distinct points or a fully collinear ring) is rejected
  instead of being silently flipped.
- **Simplicity certificate.** Before any cap or wall is built the normalized
  rings are concatenated into one `Polygon2d` and passed to
  `triangulate_polygon2`. Success certifies that the profile is a simple
  polygon with holes: the pipeline validates rings, deduplicates, removes
  collinear vertices, normalizes orientation, bridges each hole into the outer
  ring, requires a positive expected area, ear-clips, and finally checks that the
  triangle areas add up to the expected area within `kTriangulateAreaRelTol =
  1e-9`. Failure returns `Errc::degenerate` with the fixed message
  `"extrude_profile: profile is not a simple polygon with holes"`.
- **The certificate is reused, not recomputed.** The caps carry it as
  `AssembledFace::triangulation`, the top cap verbatim and the bottom cap through
  `detail::extrude_reverse_triangulation(certificate, ring_sizes)`: the bottom
  cap is every ring reversed, so vertex `i` of a ring moves to `start + size - 1
  - (i - start)`; reading the top triangulation through that mapping preserves
  each triangle's cyclic order, and swapping the last two indices turns the
  resulting clockwise-about-`-normal` triangles into the counter-clockwise
  winding the assembler expects. The assembler **validates** any supplied
  triangulation (indices in range, counter-clockwise about the face normal,
  covered area within the tolerance) and re-triangulates when welding or winding
  normalization changed the vertex order it refers to. The caps themselves stay
  polygon-with-holes faces; the triangulation is a certificate and the volume
  integrator's input, never a face representation.
- **Closed-form topology and volume.** For `N` profile vertices in total:
  `V = 2N`, `E = 3N` (`N` bottom ring edges, `N` top ring edges, `N` verticals),
  `F = N + 2`, `chi = 2`. The volume is `(A_outer - sum A_holes) * height`,
  exact because the caps are planar and the walls are planar quads, so the
  internal triangulation is exact.

## Failure contract

Checks run in this order, before any B-rep is assembled; every rejection leaves
the input unchanged and returns a stable code and message. The assembler's
internal message is never exposed (the plan's ruling R9).

| Check | `Errc` | Message |
| --- | --- | --- |
| tolerance non-finite or not positive | `tolerance_violation` | `extrude_profile: tolerance must be finite and positive` |
| plane frame failure | the frame's own code | the frame's own message |
| height non-finite | `non_finite` | `extrude_profile: height is not finite` |
| height not positive | `degenerate` | `extrude_profile: height must be positive` |
| profile coordinate non-finite | `non_finite` | `extrude_profile: profile coordinates are not finite` |
| malformed ring layout | `degenerate` | `extrude_profile: profile is degenerate` |
| ring orientation undefined | `degenerate` | `extrude_profile: profile ring is degenerate` |
| repeated consecutive vertex, closing duplicate included | `degenerate` | `extrude_profile: profile ring has a repeated vertex` |
| not a simple polygon with holes | `degenerate` | `extrude_profile: profile is not a simple polygon with holes` |
| assembly failure | `invalid_topology` | `extrude_profile: solid could not be assembled` |

The hole checks are covered by the certificate: a hole outside the outer ring,
nested holes, and a self-intersecting ring are all rejected as "not a simple
polygon with holes". `tests/unit/triangulate2_test.cpp` pins the mechanism:
bridging a misplaced hole produces a zero-width slit that ear clipping cannot
consume (it requires strictly positive turns, `triangulate2.hpp:355`), and a hole
straddling the boundary fails the covered-area invariant.

## Alternatives rejected

- **Triangulated caps** (emit the certificate's triangles as faces): breaks the
  polygon-face representation the rest of the library relies on, multiplies the
  face count by the triangle count, and would make a planar cap look like a
  triangle soup to `mass_properties`, `tessellate` and future face queries. The
  triangles are used internally, not as faces.
- **An area-only or point-in-polygon gate:** neither sees self-intersection, and
  writing a second edge-crossing checker would create a second source of truth
  next to `triangulate2`. The certificate is the same predicate the caps need
  anyway, so the gate cannot be weaker than the code it guards.
- **A new geometric simplicity predicate now:** more machinery, a new tolerance
  story, and it would still have to agree with the triangulation used for the
  volume. Deferred until a case appears that the triangulation accepts but that
  is not simple.
- **Boolean union of one prism per ring, or of cap/wall half-spaces:** drags in
  the SSI and coplanar-face machinery for a shape whose topology is known in
  closed form, and risks slivers where cap and wall meet.
- **Faceting the profile or sampling the walls:** contradicts exactness; every
  face here is exactly planar.
- **Oblique extrusion (direction vector, shear) now:** caps stop being planar
  and walls stop being rectangles; that is M4d-2.
- **Re-triangulating the caps (the first Task 3 form):** three O(n²)
  triangulations per call (certificate plus both caps). Measured 439226 us per
  call on a 512+256 profile before the reuse and 176075 us after, 2.49x;
  rejected as waste.

## Known limits

- **The certificate is an operational predicate** ("the triangulation
  succeeds"), not a proof of simplicity. The area invariant alone is not a
  complete guard: for a self-touching ring, ear clipping emits two triangles
  that already cover the expected area and rejects only because it then gets
  stuck. A hypothetical non-simple profile that stays fully ear-clippable and
  area-consistent could pass; none of the probed shapes did, and the caps would
  still be triangulated by the same code, so the gate cannot be weaker than what
  it guards.
- **Holes must not touch each other or the outer ring.** Bridging needs a
  strictly positive-width slit; touching rings are rejected rather than repaired.
- **Ear clipping is O(n²).** A 1024-gon costs about 300 ms in Release; there is
  no monotone decomposition yet.
- **The reuse hint is advisory.** Any structural mismatch (welding merged a
  vertex, winding normalization reversed a ring) silently falls back to
  re-triangulating, so the result is correct but the optimization is skipped.
- **Scope.** Planar profile with optional holes, one straight extrusion along
  the plane normal, one positive height, one solid. No oblique/shear, no
  revolve/loft/sweep (M4d-2), no self-intersection repair, no multiple solids,
  no tapering.

## Consequences

- Extrusion of a simple planar profile with holes is exact, deterministic and
  returns a valid closed B-rep: `2N/3N/(N+2)` topology, `chi = 2`, volume
  `(A_outer - sum A_holes) * height`. For a 1024-gon of radius 10 and height 1
  the benchmark measures volume `314.157294036708493` against the analytic
  `314.157294036709118` (relative error about `2e-15`); for the 512-gon ring
  with a 256-gon hole, `235.619448841109261` against `235.619448841111847`
  (about `1.1e-14`).
- The benchmark (`bench/bench_m4d1.cpp`, Release, five runs) measures
  `286.401 / 297.212 / 322.401 ms` for the 1024-gon and
  `126.939 / 143.930 / 152.921 ms` for the holed profile; the certificate alone
  costs `103.813 / 127.637 / 152.904 ms` (768 triangles), so the extrusion is
  within about 10–15 % of one triangulation — the caps' reuse works and the
  remaining cost is the O(n²) certificate.
- Tests pin: square, L-shape, hole, reversed outer ring, tilted plane, holed
  profile on a tilted plane, determinism and input immutability, the 1e6 square
  with height 1e-6, every rejection message above, bowtie and outside-hole
  rejection, the assembler-side hint reuse and its three invalid forms, and the
  two structurally incompatible hints (reoriented ring, welded vertex) that must
  be dropped rather than validated.
  Debug and Release `*extrude*` are 14 cases / 99 assertions; the filtered
  extrusion lane (`*extrude*,assemble:*,*triangulate2*`) is 49 / 10461, and the
  full Release suite is 1208 cases / 685744 assertions, 0 failed.
- Deferred to M4d-2: oblique/shear extrusion, revolve, loft and sweep, plus any
  profile repair or self-intersection healing. No document claims general sweep
  support.
