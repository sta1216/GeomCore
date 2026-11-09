# ADR 0023: M3a 的平面优先 B-rep 布尔

- Status: Accepted
- Date: 2026-11-10
- Related: architecture §8.8, §16.2, §16.3; M3a plan Task 1–Task 11; ADR-0021, ADR-0022

## 背景

M2 delivered planar and analytic-curved primitives but no boolean. M3 needs a boolean,
and it needs surface-surface intersection (SSI) to cut curved faces. SSI plus
curved-face splitting is the high-risk piece: surface-surface intersection,
pcurve-on-surface, and tolerance propagation across trimmed faces. A planar
cut-and-classify pipeline exercises the same topology, assembly and tolerance
machinery, but on a bounded problem where every intersection is a straight segment.

## 选项

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

## 决策

M3a implements option 2. `boolean_solid` supports planar-face solids only. It pairs
faces by AABB, intersects plane pairs into cut segments, cuts faces into fragments,
classifies fragment samples against the other operand, flips inward fragments, and
assembles a validated solid. Curved faces are rejected with `Errc::unsupported`;
every pcurve is a straight 2D segment. Multi-shell results (an outer shell plus
cavity shells) are supported: `assemble_solid` groups the shells into outer and
cavity shells and uses solid-angle winding for cavity containment. `validate_topology`
— called by the tests and available to callers — adds a conservative
axis-aligned-bounds pre-check plus point/touch probes, not exact shell-shell
intersection. The Boolean entry point itself does not run `validate_topology` on
its result. A disjoint union can return a Brep with two solids (two shells,
volume summed). Operands are stricter: each must be exactly one solid with
exactly one shell. More than one solid returns `Errc::unsupported` with
`"boolean_solid: expected exactly one solid per operand"`; a multi-shell solid
returns `Errc::unsupported` with
`"boolean_solid: multi-shell operands are deferred to M3c"`.

## 后果

- M3a cannot boolean cylinder/cone/sphere/torus solids, so QA-9 and curved-face
  booleans are deferred to M3b.
- Only a partial coplanar patch that touches a host loop edge is diagnosed
  (`Errc::not_closed`) rather than stitched with T-junction propagation;
  full-face coplanar contact succeeds (stacked union volume 2, difference
  volume 1).
- Symmetric difference of edge-touching L-shapes is diagnosed as non-manifold
  (`Errc::invalid_topology`).
- Cavity containment is decided by `assemble_solid` with solid-angle winding;
  `validate_topology` adds the conservative AABB pre-check plus point/touch
  probes, not exact shell-shell intersection.
- A disjoint union can return a Brep with two solids (two shells, volume
  summed).
- Sampling offsets are max(64·tol, 1e-7·diagonal); solids thinner than that
  (for example a 1000×1000×2e-4 plate at tol 1e-9) fail loudly with
  `Errc::not_closed` (`"assemble_solid: edge is used once"`), not silently
  accepted.
- The cut-and-classify pipeline, the tolerance contract, and the diagnostic
  vocabulary are reusable by M3b.

## 重新评估的触发条件

M3b SSI work; or a user-visible need to boolean curved primitives or to propagate
T-junctions.
