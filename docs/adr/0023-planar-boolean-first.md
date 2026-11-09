# ADR 0023: Planar-first B-rep boolean for M3a

- Status: Accepted
- Date: 2026-11-10
- Related: architecture §8.8, §16.2, §16.3; M3a plan Task 1–Task 11; ADR-0021, ADR-0022

## Context

M2 delivered planar and analytic-curved primitives but no boolean. M3 needs a boolean,
and it needs surface-surface intersection (SSI) to cut curved faces. SSI plus
curved-face splitting is the high-risk piece: surface-surface intersection,
pcurve-on-surface, and tolerance propagation across trimmed faces. A planar
cut-and-classify pipeline exercises the same topology, assembly and tolerance
machinery, but on a bounded problem where every intersection is a straight segment.

## Alternatives considered

1. **SSI-first.** Build SSI and curved-face splitting first, then the boolean on top.
   This is the correct end state, but it delays any usable boolean and couples M3 to
   its hardest numerical component.
2. **Planar cut-and-classify.** Support planar-face solids only: pair faces by AABB,
   intersect plane pairs into straight cut segments, cut faces into fragments,
   classify fragment samples against the other operand, and assemble a validated
   solid. Curved faces return `Errc::unsupported`.
3. **Mesh boolean.** Convert both operands to triangle meshes and run a mesh boolean.
   This is easy but loses B-rep identity, tolerance semantics, and exact topology.
4. **Exact-arithmetic arrangement in 3D.** Build a full 3D arrangement with exact
   arithmetic. This is out of scope for a first boolean.

## Decision

M3a implements option 2. `boolean_solid` supports planar-face solids only. It pairs
faces by AABB, intersects plane pairs into cut segments, cuts faces into fragments,
classifies fragment samples against the other operand, flips inward fragments, and
assembles a validated solid. Curved faces are rejected with `Errc::unsupported`;
every pcurve is a straight 2D segment. Multi-shell results (an outer shell plus cavity
shells) are supported, with conservative AABB cavity-disjointness.

## Consequences

- M3a cannot boolean cylinder/cone/sphere/torus solids, so QA-9 and curved-face
  booleans are deferred to M3b.
- Coplanar patches that touch a host face boundary are diagnosed
  (`Errc::not_closed`) rather than stitched with T-junction propagation.
- Symmetric difference of edge-touching L-shapes is diagnosed as non-manifold
  (`Errc::invalid_topology`).
- Cavity/shell intersection is conservative (AABB).
- The cut-and-classify pipeline, the tolerance contract, and the diagnostic
  vocabulary are reusable by M3b.

## Revisit trigger

M3b SSI work; or a user-visible need to boolean curved primitives or to propagate
T-junctions.
