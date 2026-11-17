# ADR-0032: General synthetic-domain splitter (partial: closed cuts on sphere charts)

Status: accepted for M3c-3c-2a (scope narrowed during implementation, 2026-11-18)

## Context

`split_face_curved` treated a curved face as *seam-only* when its stored uv loop had
zero area, then opened the surface onto a synthetic parameter rectangle. That worked
for sphere caps and bands but had three defects:

1. A synthetic chart that already carried a stored inner ring was split on its outer
   ring only; the hole was silently dropped.
2. Whether a chart was synthetic was decided from the stored loop's area. A torus
   chart is the whole rectangle (area `4*pi^2`), so the torus never took the
   synthetic path and the public entry could return zero fragments without an error.
3. Torus bands that wrap the periodic-v seam have no representation in the
   seam-only emitter: their two halves are emitted as zero-area slit rings or as
   rings missing the v-seam rim, so assembly cannot pair them.

## Decision

- Introduce a domain descriptor (`TrimFaceDomain`) and an arrangement-based splitter
  (`trim_face_split_domain`) for closed cuts, with explicit refusals:
  - `unsupported "split_face_curved: stored holes are not supported yet"` when the
    synthetic chart carries an inner loop;
  - `unsupported "...closed cut crosses the seam at different transverse
    coordinates; twisted seam crossings are deferred to M3c-3c-4"` for a twisted
    seam crossing;
  - `numeric_failure "...nested closed cuts are not supported"` for nested or mixed
    island/non-island cuts;
  - `numeric_failure "...fragments do not partition the chart"` when equal islands
    overlap.
- `split_face_curved` routes closed cuts on **sphere** synthetic charts through the
  general splitter. Stored holes on such charts are refused before the splitter
  runs; previously the inner ring was ignored.
- Torus planes keep the M3c-2 contract and guard. Offset/tilted torus plane
  booleans are deferred to M3c-3c-2b.

## Evidence

Probes for this milestone tried two ways to make offset horizontal torus planes work
without changing the fragment contract:

- shifting the v origin onto one cut latitude;
- emitting full-width v-side seam rings.

Neither passed `curved boolean: an offset torus plane splits both ways`; both
regressed five centered-plane torus booleans (`touches itself along an edge`,
`arc run is inconsistent`, `not a single connected component`, and a wrong error
code for the difference pin). The measured blocker is the assembler: harvested
fragments do carry healed spans (e.g. ring sizes 46/46, 45/49, 69/69), but the rims
of adjacent torus bands and the planar annulus cannot be paired, so one run is
assembled and an edge stays used once.

The working tree with those experiments is preserved only as
`.superpowers/sdd/m3c3c2a-splitter/task-3-blocked.patch` (not committed).

## Consequences

- No reachable public input silently loses geometry: holes, twisted/nested/mixed
  cuts and overlapping islands are explicit errors.
- Sphere closed-cut splitting is unchanged in behavior (same `1/2/1` band structure)
  and now refuses stored holes explicitly.
- M3c-3c-2b must define a two-rim band contract for periodic-v charts, pair rims in
  the assembler, and then enable offset/tilted torus planes and twisted crossings.
- Multi-cycle holes and the ear-clipping cost remain deferred (M3c-3c-4).
