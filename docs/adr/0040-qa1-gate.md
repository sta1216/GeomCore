# ADR-0040: QA-1's 80 ms union budget is unattainable; propose a revised gate

- Status: proposed (requires owner approval; the benchmark threshold is not changed by this ADR)
- Date: 2026-11-23
- Supersedes: the QA-1 acceptance number in the architecture's performance table, if approved
- Related: ADR-0039 (M3.5a performance), `docs/benchmarks/m3-5.md`

## Context

QA-1 measures the union of two 50,000-vertex stars, offset by 0.35, and has been
gated at 80 ms since M1b. The M3.5a slice deliberately left QA-1 untouched; this
slice measured it and tested the only cheap structural win suggested by the
profile.

The Release measurement is stable at **920-970 ms** on the reference machine
(median 923 ms over interleaved A/B runs; the same binary also produced 887-1,037
ms under load). The stage split for one union call (914 ms instrumented run) is:

| stage | ms | notes |
| --- | ---: | --- |
| `split_segments` | 534 | bbox enumeration 242; cheap disjoint reject 120; 44,462 proper crossings 77; shared-endpoint handling, setup, sorts, emit 83 |
| arrangement (DCEL) | 340 | edges 93, fan 73, faces 64, vertices 49, winding 38, rest 23 |
| `boolean2` post-build | 41 | classify, trace, assemble |
| total | 914 | 100,000 input segments |

Counts: 6,852,366 bbox candidates (68 per segment), 144,462 recorded pairs
(100,000 shared-endpoint plus 44,462 proper crossings), 188,924 output edges,
V/E/F = 144,462 / 188,924 / 44,464, one output polygon with 94,908 vertices,
area `31218.686993739982` (relative error 1.78e-13).

## Decision

Do not claim the 80 ms gate, and do not tune the benchmark to hide the gap.
Propose a revised QA-1 gate of **1,000 ms** for this fixture, with the honest
stage budget above, and keep 80 ms withdrawn as unattainable for this input.

Evidence that 80 ms is out of reach for this fixture and this contract:

- Pair enumeration alone costs 242 ms (6.85M candidate pairs before any exact
  test). No exact arrangement can be built without at least considering those
  pairs; a uniform-grid experiment produced 9.1-23.9M candidates and 2.0-4.2 s
  of deduplication, i.e. it is worse, not better.
- The realistic single-thread floor with a Bentley-Ottmann sweep and a flat DCEL
  is about 300 ms; 150-200 ms with eight threads. Both are above 80 ms.
- The output itself has 94,908 vertices and 188,924 edges; merely allocating and
  sorting that output is a measurable fraction of the budget.

The cheap inline disjoint check was implemented and measured A/B against the
baseline with interleaved runs: base 919/923/965 ms, inline 928/931/973 ms. The
MSVC Release build already inlines `split_record_pair`, so hoisting the check
does not recover the probe build's 100 ms; the change was reverted and is not
part of this slice. This is recorded so the next slice does not repeat it.

## Consequences

- The benchmark keeps failing at 80 ms until the owner approves a new number;
  the FAIL is documented rather than hidden.
- A revised 1,000 ms gate would be met by the current code with margin, and the
  remaining optimization path (flat DCEL, sweep, threading) is tracked as
  future work rather than silently bundled into a gate change.
- The QA-9 gate is closed separately by the exact analytic sphere-plane path
  (median 53.5 us, measured in `docs/benchmarks/m3-5.md`).

## Alternatives considered

- **Tune the fixture** (fewer vertices, no offset): rejected; QA-1 is meant to
  stress a large interacting union, and shrinking it would hide the cost.
- **Cache the union result across calls**: rejected; the benchmark times one
  call and caching would not represent the algorithm.
- **Declare 80 ms met by measuring only the split stage**: rejected; dishonest.
- **Adopt 400 ms now**: rejected as a gate because the current implementation
  does not meet it; 400 ms belongs in the optimization roadmap, not in the
  acceptance table.
