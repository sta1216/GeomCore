# ADR-0043: Convex prism fillet as direct cylindrical fragments

- Status: accepted (M4c-3b)
- Date: 2026-11-27
- Related: architecture §8.9, §16.4; `docs/benchmarks/m4c3b-fillet.md`; plan
  `docs/superpowers/plans/2026-11-27-m4c3b-prism-fillet.md`; ADR-0042 for the
  box sibling and ADR-0041 for chamfer/shell

## Context

§8.9 lists rolling-ball fillet next to chamfer and shell. ADR-0042 shipped the
first exact slice: on an axis-aligned box, a selected edge joins two
perpendicular planes, so its fillet is a quarter cylinder bounded by cap-plane
arcs. It deliberately left the vertex-blend and general-operand work to
M4c-3b.

A convex prism is the next closed case. Its vertical edges are parallel, each
joins two planar side faces, and both caps are planes perpendicular to the
common axis. The rolling-ball surface along one vertical edge is therefore an
exact cylinder of the requested radius, and the surface ends where the cap
planes cut that cylinder in circular arcs. Selected corners do not interact in
three dimensions: their tangent trims only shorten the cap edges and the
vertical strips between them, which is a two-dimensional side-width check.

The question was whether to ship this slice now, with its own recognition and
radius guards, or to wait for the generic edge-chain/SSI blend engine.

## Decision

M4c-3b extends `fillet_edges(const Brep&, span<const EdgeId>, const
FilletOptions&)` to a single-shell, closed, convex, all-planar prism:

- **Recognition.** The solid must have exactly two cap faces with opposite
  outward normals and one common axis, matching cap rings, and a strictly
  convex polygon cross-section (checked with exact predicates, not a tolerance
  sweep). The candidate must also pass the existing planar-convex-solid
  recognition. All filletable edges are the ring-to-ring edges parallel to the
  prism axis.
- **Direct fragments, not a boolean.** Each selected vertical edge gets an
  exact `Cylinder3` patch of radius `r`, with its axis parallel to the prism
  axis, its center on the interior angle bisector, and its surface tangent to
  both adjacent side planes. The cylinder's ends are trimmed by
  `CircleCurve3` arcs in the two cap planes. The cap loops and the side faces
  are rebuilt once, trimmed at the tangent points, assembled with
  `assemble_curved_solid`, and returned only after `validate_topology` is
  clean. No Boolean, no sampled fit and no faceted ring are involved.
- **Scope.** Any non-empty subset of the vertical edges may be selected;
  duplicate ids collapse and the result is independent of selection order.
  The radius is uniform across the selection. Adjacent selected corners simply
  trim the shared cap edge from both ends. Non-convex polygons, curved side
  faces, multi-shell solids, vertical-edge subsets that do not exist, and
  non-vertical edges are out of scope.
- **Radius limit is per selected corner, not a global inradius cap.** At a
  corner with interior angle `theta`, the tangent length along each incident
  side is `t = r * cot(theta/2)`. The corner is accepted when `t` is strictly
  less than both incident side lengths. This is the exact condition for the
  cylinder to meet both sides inside the corner; it deliberately does not
  reject a radius merely because it exceeds the polygon's inradius. After all
  corners are trimmed, every side must retain a positive strip:
  `side_length - left_trim - right_trim > tolerance`; otherwise the case is
  rejected rather than producing a sliver face.
- **Volume and topology are analytic.** For a prism of cross-section area `A`
  and height `h`, with selected corners `S` and interior angles `theta_i`,

  `Vol = A*h - h * sum_{i in S} r^2 * (cot(theta_i/2) - (pi - theta_i)/2)`.

  The per-corner removed area is the tangent triangle minus the circular
  sector, which is exactly the material the cylinder removes. The counts are
  `V = 2n + 2k`, `E = 3n + 3k`, `F = n + 2 + k` for an `n`-gon with `k`
  selected vertical edges. The all-64-edges 64-gon fixture therefore has
  `256/384/130`, and `chi = 2`.
- **Failure contract.** Every rejection happens before assembly, with a stable
  `Errc` and message, and leaves the input unchanged:
  - empty selection: `unsupported "fillet_edges: at least one edge is
    required"` (checked before the options, because the default radius is zero);
  - non-finite radius or tolerance: `non_finite "fillet_edges: radius must be
    finite"` / `"fillet_edges: tolerance must be finite"`;
  - non-positive radius: `degenerate "fillet_edges: radius must be positive"`;
  - non-positive tolerance: `tolerance_violation "fillet_edges: tolerance must
    be positive"`;
  - a solid that is neither a supported box nor a convex prism: `unsupported
    "fillet_edges: only a box or a convex prism is supported"`;
  - a dead/foreign id on a recognised prism: `index_out_of_range
    "fillet_edges: edge does not belong to the solid"`;
  - a horizontal (non-vertical) prism edge: `unsupported "fillet_edges: only
    vertical prism edges are supported"`;
  - `r * cot(theta/2)` not strictly inside both incident sides: `unsupported
    "fillet_edges: radius exceeds the prism cross-section"`;
  - a post-trim side strip at or below the tolerance: `unsupported
    "fillet_edges: radius leaves a prism face narrower than the tolerance"`;
  - an assembly that cannot be closed or fails validation: `invalid_topology
    "fillet_edges: fillet could not be assembled"` or `"fillet_edges: fillet
    produced an invalid solid"` (defensive; the accepted cases above are
    expected to assemble).

The box path remains first and unchanged; the prism path is used only when box
recognition fails. This ADR records a narrow exact case, not general fillet
support.

## Alternatives rejected

- **Generic sweep / edge-chain blend engine now.** Edge chain → local patch →
  SSI → stitch would cover non-parallel and concave edges, but it needs curved
  intersection curves and trimming machinery that M4c-3c still owes. For a
  prism the arithmetic is already exact without it, so shipping the narrow case
  does not block or pre-empt the general engine.
- **Spherical vertex blend.** Where adjacent filleted edges meet at a vertex,
  a rolling-ball blend is a spherical patch. Prism vertical edges meet only at
  cap vertices, and their cap arcs already meet at the trimmed cap boundary;
  no spherical patch is needed. Blends for non-parallel or adjacent
  non-vertical edges remain M4c-3c.
- **Faceted cylinder approximation.** A chain of planar strips is not the
  requested surface, has unbounded deviation without a curvature bound, and
  contradicts the project rule that a faceted substitute must not be called a
  fillet (ADR-0041, ADR-0042).
- **Boolean round-trip** (`body - (corner wedge - cylinder)`): rejected for the
  same reasons as in ADR-0042 — it needs exact curved intersections, inherits
  splitter face-order sensitivity, and is numerically worse than the analytic
  patch.
- **Global inradius cap** (`r < inradius`): rejected as both too strict and not
  the real constraint. The cylinder is limited by the two incident side
  lengths through `r * cot(theta/2)`, which permits a single corner of a
  triangle to be rounded at `r = 0.6` even though that exceeds the polygon's
  inradius. The per-corner test is exact; the post-trim strip check keeps every
  remaining face non-degenerate.
- **Sampled/fitted surface instead of an exact cylinder.** The surface is
  exactly a cylinder, so a fit only adds error, a sampling policy and a
  tolerance story.

## Consequences

- Constant-radius fillet of any non-empty subset of a convex prism's vertical
  edges is exact, deterministic and returned as a valid closed B-rep. On the
  64-gon, radius 2, height 10 fixture with all 64 vertical edges at `r = 0.05`,
  the benchmark measures `11.683 / 14.477 / 15.385 ms` (min/median/max) and
  volume `125.461876478546259` against the analytic
  `125.461876478546088` (relative error about `1.4e-15`). Five timed runs reuse
  the same input; checks and volume integration run outside the timed region.
- Unit tests pin one-edge, all-edge, adjacent-pair, reordered/duplicated
  selection, the face-width guard, one corner beyond the polygon inradius, a
  tilted frame, dead/empty inputs, exact messages, and input immutability plus
  deterministic tessellation. The `*prism*` filter is 11 cases / 98 assertions
  in both configurations; the prism file itself adds 13 cases / 127 assertions
  (three of its cases do not contain the word "prism" in their names).
- The caps are arc-trimmed surfaces, so `mass_properties` does not accept them;
  volume evidence uses `detail::curved_solid_volume`.
- The benchmark checks topology, validation and the analytic volume. It does
  not tessellate, so it makes no watertightness claim; tessellation
  watertightness remains the unit-test boundary check.
- Deferred to M4c-3c: adjacent non-vertical edge blends, spherical vertex
  patches, variable radius, non-convex polygons, multi-shell solids, curved
  side surfaces and general operands. No document claims general fillet
  support.
