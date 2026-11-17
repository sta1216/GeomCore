# M3c-3c-2a benchmark note

M3c-3c-2a generalizes the synthetic-domain splitter and adds explicit refusals
(stored holes, twisted seam crossings, nested/mixed closed cuts, overlapping
islands). It introduces no new algorithm path for the supported corpus: sphere
closed-cut splitting keeps the same fragment structure, and torus plane booleans
keep the M3c-2 contract.

Therefore this milestone adds no benchmark executable and no timing gate. The
existing measurements remain representative:

- M3c-1 trimmed-face tessellation and M3c-3a sphere caps are unchanged; see
  `docs/benchmarks/m3c1.md` and `docs/benchmarks/m3c3a.md`.
- QA-18a is a structural corpus (fragment and loop counts, determinism, explicit
  refusal messages), not a performance gate.
