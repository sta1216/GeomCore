# ADR-0041: Planar chamfer and shell as exact half-space operations

- Status: accepted (M4c-1; extended by M4c-2)
- Date: 2026-11-25 (M4c-2 extension added 2026-11-26)
- Related: architecture §8.9, §16.4; `docs/benchmarks/m4c1-blends.md`,
  `docs/benchmarks/m4c2-shell.md`; plans
  `docs/superpowers/plans/2026-11-25-m4c1-planar-chamfer-shell.md`,
  `docs/superpowers/plans/2026-11-26-m4c2-open-shell.md`

## Context

M4a/M4b healed open rims on analytic faces, but the library had no blend
operations at all: §8.9 listed rolling-ball fillet, distance/angle chamfer and
shell without an implementation. A full blend engine (edge chain → local patch
→ SSI → stitch) is a large slice and its generic version needs curved
intersection surfaces.

Two operations of that table reduce exactly to half-space cuts on a convex,
planar solid:

- **Chamfer**: one plane through the setback points on the two faces adjacent to
  a straight edge, applied as a cut that keeps the material side. A chain of
  edges is a sequence of such cuts; adjacent cuts meet at a vertex facet
  naturally, without a special vertex blend.
- **Shell**: a uniform wall on a closed, axis-aligned convex solid. All faces
  are axis-aligned, so the solid is a box; its cavity is the same box inset by
  the wall thickness on every side, and the result is the regularized difference
  `A - cavity`.

The exact predicates from L1 decide every cut side; no new tolerance model and
no global epsilon are introduced.

## Decision

M4c-1 implements exactly the two operations above:

- `chamfer_edges(const Brep&, span<const EdgeId>, const ChamferOptions&)` on a
  single-shell, closed, convex, all-planar solid with straight `LineCurve3`
  edges. Failure returns a stable `Errc` and a message; the input is unchanged.
- `shell_solid(const Brep&, span<const FaceId>, const ShellOptions&)` on a
  single-shell, closed, convex, axis-aligned, all-planar solid, with no open
  face. The cavity is the inset box; a wall that leaves no positive interior is
  `Errc::degenerate`.
- Every unsupported shape is refused before any mutation: curved faces,
  non-convex bodies, open or multiple shells, oblique faces, open shell
  requests, non-`LineCurve3` edges, edges without exactly two incident faces
  and a chamfer that is wider than an adjacent face.
- **Fillet is not approximated with planar facets.** Variable-radius fillet,
  chamfer angle, non-convex and curved operands, open shells, offset and sweep
  remain M4c-2. A faceted "fillet" would look like a blend but fail the
  geometric contract, so it is explicitly rejected.

## M4c-2 extension (2026-11-26)

M4c-2 supersedes the two M4c-1 shell limits above:

- **Open faces are supported.** `shell_solid` accepts any non-empty proper
  subset of the input faces in `open_faces`. For an open face the cavity
  half-space is relaxed outward by the wall thickness, so the pocket breaks
  through that face and the result has `cavity_shells == 0` and one shell.
  Duplicate ids are accepted once; opening every face is refused with
  `Errc::unsupported` and the message
  `shell_solid: at least one face must stay closed`, because no wall would
  remain.
- **Oblique convex faces are supported.** The closed-shell domain is every
  closed, convex, single-shell, all-planar solid, not only an axis-aligned box.
  The cavity is the intersection of the inward-offset face half-spaces with a
  bounding seed box, and the result is the regularized difference `A - cavity`.
  The direct inset-box shortcut is retained only for the all-axis-aligned case.

The M4c-1 rejection messages `shell_solid: open faces are not supported in
M4c-1` and `shell_solid: only axis-aligned faces are supported` are therefore
superseded and no longer emitted. What remains rejected before mutation is:
curved or non-planar faces, non-convex bodies, open or multiple input shells,
non-`LineCurve3` edges, dead/foreign open-face ids, opening every face, and a
wall that leaves no measurable interior (`Errc::degenerate`).

**Fillet stays deferred.** Rolling-ball fillet, including variable radius,
is not part of M4c-2 and moves to M4c-3; a faceted approximation is still
explicitly rejected.

## Alternatives rejected

- **Polyhedral fillet approximation**: a chain of planar strips around a
  rolling-ball arc. Rejected: the result is not the requested surface and the
  deviation is unbounded without a curvature bound; it would be a half-finished
  feature.
- **Boolean tool solids for chamfer**: build a wedge and subtract it from the
  body. Rejected: constructing a valid wedge per edge needs the same setback
  planes, but the boolean adds stitching and order dependence without adding
  accuracy; the cut keeps `above` directly.
- **Shell by offsetting every face then sewing**: offsetting a face and
  extending it to intersections reintroduces the splitter's face-order
  sensitivity (observed while implementing the seed-box variant). For the
  axis-aligned slice the inset box is exact and order independent.
- **Open pockets via an outward cut**: shifting the open face's plane outward
  and keeping `below` is the identity for a convex body; a seed-box
  intersection produces the pocket but the splitter's cut-chain stitching is
  order sensitive (failed on 312/720 face orders for one box). Deferred in
  M4c-1 with the splitter/boolean face-order hardening. **Superseded by the
  M4c-2 extension above**: the order-independent half-space cavity plus the
  regularized difference now ships open shells without depending on splitter
  cut-chain order.

## Consequences

- Chamfer and shell are exact for the stated domain and reuse `SolidEditor`,
  `split_solid_by_plane`, `boolean_solid` and the existing validators.
- The M4c-1 shell domain was effectively boxes; M4c-2 widens it to every
  closed convex planar solid and to non-empty proper open-face subsets, as
  recorded in the extension above and in §8.9 and §16.4.
- The splitter's cut-chain clipping remains a documented weak point for oblique
  cuts and coplanar multi-edge unions; it is listed as an M4c-2 item, not hidden
  behind a silent fallback.
- Tests pin the analytic volumes, face-order independence, watertight
  tessellation, determinism and unchanged inputs; the benchmark records the
  measured cost.
