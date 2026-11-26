# ADR-0042: Constant-radius box fillet as a direct cylindrical assembly

- Status: accepted (M4c-3a)
- Date: 2026-11-26
- Related: architecture §8.9, §16.4; `docs/benchmarks/m4c3-fillet.md`; plan
  `docs/superpowers/plans/2026-11-26-m4c3-box-fillet.md`; ADR-0041 for the
  sibling chamfer and shell operations

## Context

§8.9 lists rolling-ball fillet next to chamfer and shell. M4c-1 and M4c-2
delivered chamfer and shell as exact half-space operations and explicitly
deferred fillet: the generic blend engine (edge chain → local patch → SSI →
stitch) needs curved intersection curves, and ADR-0041 rejects a faceted
approximation as a substitute.

One slice of that engine is analytically closed. On an axis-aligned box every
edge is straight and both adjacent faces are perpendicular planes, so the
rolling-ball surface along the edge is exactly a quarter cylinder of the
requested radius; where the edge ends, the cylinder is bounded by a circular
arc in the end plane. No intersection surface, no tolerance-driven fitting and
no boolean are involved. The question was whether to ship that narrow exact
case now or wait for the general engine.

## Decision

M4c-3a implements exactly that slice:

- `fillet_edges(const Brep&, span<const EdgeId>, const FilletOptions&)` on a
  single-shell, closed, convex, all-planar, axis-aligned box. It accepts any
  non-empty subset of one parallel edge family (1–4 edges); duplicate ids
  collapse to one edge and the result does not depend on the selection order.
- **Exact cylindrical assembly, not a boolean.** Each selected edge gets a
  quarter-cylinder `Cylinder3` patch of the requested radius. The two cap
  planes and the four side planes are rebuilt once and trimmed by
  `CircleCurve3` arcs at every selected corner. The faces are assembled with
  `assemble_curved_solid` and returned only after `validate_topology` is clean.
- The volume and topology are analytic: for `k` rounded edges of length `L`,
  `Vol = a*b*c - k * (1 - pi/4) * r^2 * L`, and the counts are `V = 8 + 2k`,
  `E = 12 + 3k`, `F = 6 + k`, `chi = 2`.
- Every unsupported input is refused before assembly with a stable `Errc` and
  message, and the input is unchanged on failure: empty span
  (`unsupported "fillet_edges: at least one edge is required"`); non-finite
  radius or tolerance (`non_finite "fillet_edges: radius must be finite"` /
  `"fillet_edges: tolerance must be finite"`); non-positive radius
  (`degenerate "fillet_edges: radius must be positive"`); non-positive
  tolerance (`tolerance_violation "fillet_edges: tolerance must be positive"`);
  a non-box solid
  (`unsupported "fillet_edges: only an axis-aligned box is supported"`); an
  edge that is not axis-aligned
  (`unsupported "fillet_edges: only axis-aligned box edges are supported"`);
  mixed families (`unsupported "fillet_edges: selected edges must be parallel"`);
  an edge that is not a box corner edge
  (`unsupported "fillet_edges: edge is not a box corner edge"`); a dead id
  (`index_out_of_range "fillet_edges: edge does not belong to the solid"`);
  two selected edges sharing a vertex
  (`unsupported "fillet_edges: filleted edges must not share a vertex"`); and
  `radius >= 0.5 * min(cross-section extents)`
  (`unsupported "fillet_edges: radius exceeds the box cross-section"`).
- **Empty span is checked before the options.** `FilletOptions{}` defaults to
  `radius = 0.0`, so validating the radius first would report `degenerate` for
  the empty+default call. "No edges selected" is an unsupported request, not a
  bad radius, so the empty check comes first and returns `unsupported`.

## Alternatives rejected

- **Faceted fillet ring**: approximate the quarter cylinder with a chain of
  planar strips. Rejected: the result is not the requested surface, the
  deviation is unbounded without a curvature bound, and this project already
  refuses facets as a fillet substitute: ADR-0041 defers fillet to M4c-3 and
  the M4c-2 architecture note repeats that multi-plane fitting must not stand
  in for it (chamfer and shell reject their own, different, approximations).
- **Multi-plane fit near the edge**: reconstruct a local surface from sampled
  points and fit planes. Rejected: the surface is exactly a cylinder, so the
  fit only adds error and a sampling policy; it is a numerically worse version
  of the direct patch.
- **Boolean round-trip** (`body - (corner box - cylinder)` or a swept tool):
  rejected for this milestone. It needs exact curved-curved intersection curves
  (the SSI work M4c-3 still owes), reintroduces the splitter's face-order
  sensitivity recorded in ADR-0041, and adds no accuracy over the analytic
  patch.
- **General SSI blend engine now**: edge chain → local patch → SSI → stitch for
  arbitrary convex and concave edges. Deferred: it is a large slice and the
  arithmetic for the box case is exact without it; shipping the narrow case
  keeps the contract honest instead of substituting facets.

## Consequences

- Constant-radius fillet of 1–4 parallel box edges is exact, deterministic and
  watertight when tessellated in the four-edge case tested here; the benchmark on a 100³ box with four vertical
  edges at `r = 1` measures 7.587 / 8.247 / 12.766 ms (min/median/max) and
  volume `999914.159265358`, within `1e-6` relative of `999914.159265359`.
- The tests pin the per-edge volume `23.94634954084936` and the four-edge volume
  `23.78539816339745` for a `2x3x4` box at `r = 0.25` within `1e-9`, the
  `V/E/F` counts for every `k`, watertight tessellation of the four-edge case,
  bit-identical repeated calls, an unchanged input on success and failure, the
  tested option/shape rejection codes, and the mixed-family rejection message
  (substring `parallel`).
- Deferred to M4c-3b: vertex blends where filleted edges meet, variable
  radius, concave edges, general convex and curved operands, chamfer angles,
  offset and sweep. No document claims general fillet support.
- Fillet patches are arc-trimmed surfaces, so `mass_properties` (which handles
  full cylinder parameter rectangles or polyline trims) does not accept them;
  volume evidence uses `detail::curved_solid_volume`.
