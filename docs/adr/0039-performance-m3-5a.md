# ADR-0039: Query and meshing performance for M3.5a

Status: accepted for M3.5a (2026-11-23). Narrows the QA-4/QA-10 performance gaps; QA-1 and QA-9 remain open for the next slice.

## Context

The M1b/M2/M3 gates measured four gaps on this machine (Release, i7-8550U, 4 cores):

| gate | budget | baseline |
| --- | ---: | ---: |
| QA-4 BVH build (1M triangles) | 350 ms | 1,372 ms |
| QA-4 raycast | 8 Mray/s | 0.45 Mray/s |
| QA-10 box meshing | 1.0 Mtri/s | 0.22 Mtri/s (55.6 us/box) |
| QA-1 2D union | 80 ms | ~1.9 s |
| QA-9 SSI | 100 us | 9.6-19.5 ms |

A profile (`.superpowers/perf/m3-5-profile.md`) attributed BVH build time to SAH binning (24%), the median `nth_element` (35%), child-box rescans (23%), and primitive marshalling (4%). Ray traversal spent 78% of its time in node traversal, with 36.9 pop visits per ray. Box meshing allocated 288 times per box (8,688 B), and the planar path rebuilt rings and polygons per face.

## Decision

1. **BVH build**: one pass fills the 12 SAH bins (counting scatter), prefix/suffix box unions produce both child boxes without rescanning, and the median fallback only runs inside the straddling bin. Primitives are stored SoA (triangle, box, centroid). The top three levels are built as independent frontiers and concatenated in fixed order, so thread count changes speed but not the array contents or node numbering.
2. **BVH traversal**: four-way nodes with child references and exact **double** child boxes; primitives are permuted into leaf order so leaf tests read contiguous triangles; leaf tests are inlined and the hit index maps back through `order_`.
3. **Conservative pruning**: every box prune uses `best * (1 + 1e-12)` in the five prune comparisons. Leaf acceptance stays exact (`t > best`, `t == tmax`, tie to the smallest primitive index). The slack only widens the visited set; it exists because a box entry parameter and a Moller-Trumbore hit computed from the same ray can differ by an ulp, and float child boxes were measured to prune a leaf whose true hit was one ulp closer (4000 coplanar rays: 619 wrong primitives). Float boxes were therefore rejected. Tie-aware ordered traversal was measured 8-12% slower and is not used.
4. **Box meshing**: a canonical box (6 planar faces, one 4-coedge outer loop each, 12 line edges, 8 vertices, every edge used twice, strictly convex) emits the same 12 triangles directly from the same welded indices in the same order, preserving `triangle_face` and `reversed`. Any precondition failure clears partial state and falls back to the unchanged generic loop. Rings, projected points, polygons and layout buffers are reused across faces through `TessellatePlanarScratch`; `tessellate_edge_positions` returns a span instead of a per-coedge copy.
5. **No ray packets**: the measured single-ray latency floor with primitive work removed is 0.87-0.91 Mray/s, so the 8 Mray/s budget cannot be reached without packet/SIMD work. Packets are deferred; the 0.9 single-ray target is documented as a shortfall, not silently relaxed.

## Evidence

Commands and full tables are in `docs/benchmarks/m3-5.md`. Release, three clean process runs, median:

| gate | before | after | speedup | budget |
| --- | ---: | ---: | ---: | ---: |
| BVH build | 1,372 ms | 287-300 ms | ~4.6x | 350 ms PASS |
| raycast | 0.45 Mray/s | 0.66 Mray/s | 1.47x | 8 Mray/s FAIL (documented) |
| box meshing | 0.22 Mtri/s | 1.46-1.65 Mtri/s | 6.6-7.5x | 1.0 Mtri/s PASS |
| QA-1 union | ~1.9 s | ~0.9 s | unchanged | 80 ms FAIL (next slice) |
| QA-9 SSI | 9.6-19.5 ms | unchanged | - | 100 us FAIL (next slice) |

Correctness evidence: BVH build is bitwise deterministic on the same path; serial and parallel builds produce identical public query results (node numbering may differ); the coplanar two-plane regression (2048 triangles, exact repro ray, 300 rays, 200 nearest queries) matches brute force; box output is bitwise identical to the generic path for four boxes at two deflections, including reversed faces, and canonical hashes for box/cylinder/torus/sphere are unchanged; QA-10 results remain watertight and outward. Focused Release runs: BVH 23 cases / 8,873 assertions; tessellation 82 cases / 9,226 assertions.

## Consequences

- QA-4 build and QA-10 now meet their budgets. QA-4 raycast is 1.47x faster and remains below budget; the shortfall is recorded with its latency-floor evidence.
- The wide-node record is 224 bytes (three and a half cache lines), not the 64-byte line assumed when the task was written; that is why the raycast gain is smaller than the node-width estimate.
- `detail::tessellate_impl(..., allow_box_fast_path)` and `detail::tessellate_impl(..., false)` exist for parity testing only; the public `tessellate()` always enables the fast path.
- The fast path is intentionally narrow: any non-box topology, holed face, non-line edge, concavity or duplicate vertex falls back to the generic path. A malformed 8-vertex box with an isolated vertex is not explicitly rejected by the count checks; it is mirrored by the generic path today and remains a hardening gap.
- 8 Mray/s and the 0.9 single-ray target stay open. The next performance slice covers QA-1 (2D boolean) and QA-9 (SSI), not ray packets.
