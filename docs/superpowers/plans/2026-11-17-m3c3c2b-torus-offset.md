# M3c-3c-2b: Periodic Torus Bands and Horizontal Offset Cuts

## Goal

Enable `boolean_solid_curved` on a torus cut by a horizontal plane at an arbitrary
offset (not only through the centre), by emitting each periodic-v torus band as one
whole-cycle ring and by routing synthetic torus faces through the general domain
splitter. Keep the torus plane guard, the sphere path, and the assembler unchanged.
Tilted planes, twisted seam crossings, stored holes, and trimmed periodic-v
tessellation stay explicitly unsupported and documented.

Acceptance (measured in the analyst probe, `m3c3c2b-probe-report.md`):

| op | faces | edges | vertices | shells | volume |
|---|---:|---:|---:|---:|---|
| intersection, torus ∩ box z≥0.25 | 2 | 24 | 23 | 1 | 0.9647593052603948 |
| union | 9 | 122 | 116 | 1 | 31.970042895284287 |
| difference | 3 | 110 | 108 | 1 | 3.970042895284284 |

All three `validate_topology(...).ok()`, each once with tolerance 1e-9.

## Architecture

The library represents a closed B-rep face as a list of closed 3D rings plus a
parallel list of UV rings and arc spans. For a periodic-v band (a torus chart strip
between two latitude cuts) the working representation is **one ring** that walks:

```
(0,v0) --lower rim +u--> (2π,v0) --u=2π column +v--> (2π,v1)
       --upper rim −u--> (0,v1) --u=0 column −v--> (0,v0+δ) --close--> (0,v0)
```

The ring is CCW in UV, has positive area `2π·(v1−v0)`, contains both rims and both
seam columns, and carries null arc spans on the columns. `assemble_curved` pairs
half-edge uses by undirected endpoint identity (`assemble_curved.hpp:571`), so no
assembler change is needed: the band shares its rim vertices with the neighbouring
planar annulus/cap.

The HEAD dispatch in `split_face_curved` decides "synthetic domain" by
`loop_area == 0.0 || sphere_slit`. That works for spheres (collapsed poles) but not
for a torus, whose canonical loop is the full rectangle (area `4π²`). The general
splitter also needs `trim_face_domain` classification of the four chart sides, so
`split_face_curved` must dispatch by domain, not by loop area.

## Tech Stack

C++20, CMake 4.4.4 + VS2022 (MSVC 19.44.35217), doctest 2.4.11, `/W4 /WX /utf-8`.
No new dependencies. Windows host (PowerShell); CI covers windows-2022,
ubuntu-24.04, macos-15, tsan, bench, coverage.

## Spec

- `docs/architecture/geometry-library-architecture.md` v2.1 §8.5–8.12 (curved
  booleans), §8.10 (tessellation contract), §16.3 (milestone status).
- `docs/adr/0032-general-domain-splitter.md` (domain descriptor, splitter scope,
  deferred items).
- `docs/adr/0028-curved-boolean-surfaces.md` (torus plane contract).
- `docs/quality/testing-strategy.md` (TDD, negative paths, exact messages).

## Global Constraints

- TDD per task: add the test, run it and capture the failing evidence, implement,
  then run Debug and Release full suites with `-DGEOM_WERROR=ON`.
- One commit per task, conventional message, no unrelated changes, no debug output.
- Never change `assemble_curved.hpp` in this milestone.
- Never remove or weaken `curved_boolean_torus_plane_cut_supported`; horizontal
  offset planes are already accepted by it (only tilted / axis-parallel offset
  planes are refused).
- Sphere, sphere caps, sphere bands, and sphere lenses must be byte-identical in
  behaviour; the sphere dispatch keeps `trim_face_split_seam_only`.
- Do not claim tessellation support for trimmed torus faces. `curved_solid_volume`
  and `validate_topology` are the acceptance oracles; meshing remains a known
  periodic-v gap.
- MSBuild in this repository does not track header dependencies reliably: after
  editing a header, `touch` a test source or use `/t:Rebuild`, otherwise a stale
  binary can make a test look green.
- Baseline at the start: 975 cases / 660,201 assertions, Debug and Release.

## Scope and Deviations (need approval)

1. **Horizontal offset only.** A torus cut plane must be perpendicular to the torus
   axis (`|n·axis| ≥ 1−8ε`) or contain it. Tilted planes and axis-parallel offset
   planes stay `unsupported` with the existing message. Tilt implies twisted seam
   crossings; those are deferred to M3c-3c-4.
2. **One ring per band, not two.** The blocked M3c-3c-2a Tree emitted two zero-area
   rim rings; that representation is abandoned. Sphere bands already use the
   one-ring contract; this milestone applies it to periodic-v torus bands.
3. **Tessellation stays deferred.** `tessellate()` of a trimmed torus face remains
   unsupported (pre-existing `tessellate: torus face is missing a seam` and
   `trimmed band wall is not closed` gaps). The milestone documents this and does
   not add a tessellation test.
4. **No multi-cycle holes, no cone.** Stored holes are refused with the existing
   message; nested/overlapping islands stay refused; cone apex topology is
   M3c-3c-2c/3c-4.
5. **Dispatch change is required.** Task 2 replaces the loop-area gate with a
   `trim_face_domain` dispatch; the analyst measured that the ring fix alone still
   yields `fragments=8 kept=1` and `touches itself along an edge` on the public
   path, while the dispatched path succeeds.

## Review Focus

- Task 1: ring order and directions (lower rim +u, upper rim −u, right column +v,
  left column −v), positive area `2π·height`, no zero-area ring, sphere untouched,
  all-islands complement gated by `!periodic_v`.
- Task 2: public `split_face_curved` equals the direct `trim_face_split_domain`
  result on the offset torus; a torus face with no cuts returns unchanged; stored
  holes are forwarded and refused with `context_b == 1`; sphere still uses the
  seam-only wrapper and keeps its rejection messages.
- Task 3: the three offset volumes/counts above, `validate_topology`, deterministic
  second run; centred torus regressions unchanged; all six guard refusals still
  `unsupported` with the exact message.
- Task 4: tilted / axis-parallel offset / twisted / domain-boundary / stored-hole
  refusals with exact messages; ADR wording does not promise tilted or
  tessellation support.

## File Map

| File | Task | Role |
|---|---|---|
| `include/geom/brep/trim_face.hpp` | 1, 2 | whole-cycle ring; domain dispatch |
| `tests/unit/brep_trim_face_test.cpp` | 1, 2, 4 | ring and dispatch tests, refusals |
| `tests/unit/brep_boolean_curved_test.cpp` | 3, 4 | offset acceptance, regressions, guard |
| `docs/adr/0032-general-domain-splitter.md` | 4 | contract and deferred items |
| `docs/benchmarks/m3c3c2b.md` | 4 | measured numbers and environment |
| `docs/architecture/geometry-library-architecture.md` | 4 | §16.3 status |
| `docs/superpowers/plans/2026-11-17-m3c3c2b-torus-offset.md` | 4 | this plan |

## Task 1: Emit one whole-cycle ring for a periodic torus band

**Files**
- `include/geom/brep/trim_face.hpp` (inside `trim_face_split_domain`, after the
  synthetic-edge `kept_edges == 0` guard, before the `trim_face_seam_runs` call;
  and the all-islands complement condition)
- `tests/unit/brep_trim_face_test.cpp` (`trim_face: torus latitudes split without
  the loop-area gate`)

**Interfaces**
- Consumes: `TrimFaceDomain{ bounds, periodic_u, periodic_v, sides,
  collapsed_v() }`, `Arrangement2::HalfEdge`, `vertices`, `synthetic`,
  `trim_face_normalize_ring`.
- Produces: no new public API. One `loops_uv` ring per torus fragment; the existing
  `CurvedFragment` layout is unchanged.

**Step 1: failing test**

Update the existing test at `tests/unit/brep_trim_face_test.cpp` (cuts at `v=0.6`
and `v=2.5`):

```cpp
REQUIRE(fragments.has_value());
REQUIRE(fragments.value().size() == 3u);
for (const geom::CurvedFragment& fragment : fragments.value()) {
  REQUIRE(fragment.loops_uv.size() == 1u);
  CHECK(fragment.loops_uv.front().size() > 100u);
}
CHECK(signed_area2(fragments.value().at(0).loops_uv.front()) ==
      doctest::Approx(3.7699111843077517).epsilon(1e-9));
CHECK(signed_area2(fragments.value().at(1).loops_uv.front()) ==
      doctest::Approx(11.938052083641214).epsilon(1e-9));
CHECK(signed_area2(fragments.value().at(2).loops_uv.front()) ==
      doctest::Approx(23.770454336408466).epsilon(1e-9));
for (const geom::CurvedFragment& fragment : fragments.value()) {
  const geom::Point2d front = fragment.loops_uv.front().front();
  const geom::Point2d back = fragment.loops_uv.front().back();
  CHECK(front.x == doctest::Approx(0.0).scale(1.0));
  CHECK(back.x == doctest::Approx(0.0).scale(1.0));
  CHECK(std::abs(front.y - back.y) > 1e-6);
}
```

**Step 2: run and confirm failure**

`.\build\tests\Debug\geom_tests.exe --test-case="trim_face: torus latitudes*"`
must fail: currently fragment 1 reports `loops_uv.size() == 2` and the areas are not
`2π·height`.

**Step 3: minimal implementation**

Add, immediately after the `kept_edges == 0` guard:

```cpp
bool built_whole_cycle = false;
if (domain.periodic_u && domain.periodic_v && !domain.collapsed_v()) {
  std::size_t first_kept = cycle.size();
  for (std::size_t i = 0; i < cycle.size(); ++i) {
    if (!synthetic[i]) { first_kept = i; break; }
  }
  if (first_kept == cycle.size()) first_kept = 0;
  std::vector<Point2d> whole;
  whole.reserve(cycle.size() + 1);
  bool seam = false;
  for (std::size_t k = 0; k < cycle.size(); ++k) {
    const std::size_t index = (first_kept + k) % cycle.size();
    const Arrangement2::HalfEdge& edge = half_edges[cycle[index]];
    whole.push_back(vertices[edge.origin]);
    seam = seam || synthetic[index];
  }
  whole.push_back(vertices[half_edges[cycle[first_kept]].origin]);
  trim_face_normalize_ring(whole, true);
  if (whole.size() >= 3) {
    rings.push_back(DomainRing{std::move(whole), seam});
    built_whole_cycle = true;
  }
}
if (!built_whole_cycle) {
  // existing trim_face_seam_runs / pole / seam-pair code, unchanged
}
```

Then gate the complement:

```cpp
if (all_islands && !domain.periodic_v) { /* existing keyhole complement */ }
```

`DomainRing { std::vector<Point2d> points; bool seam; }` is a local type declared
with the other ring vectors; reuse it for the existing branch to avoid duplicate
containers.

**Step 4: verify**

- The updated test passes with the three areas above (sum 39.478417604357432).
- `trim_face: a sphere band joins its two cut runs` still passes (`1` ring, area
  `2π·2v`), and all sphere/sphere-cap tests are unchanged.
- Debug and Release full suites: the case count stays 975; the assertion count
  grows by the new area/column/closure checks (about +12). Report the exact
  numbers; do not pin a guessed total in the code.
- `trim*` filter and `curved boolean:*` filter green.

**Step 5: commit**

`feat(trim_face): emit one seam-column ring for a periodic torus band`

## Task 2: Dispatch synthetic torus faces through the domain splitter

**Files**
- `include/geom/brep/trim_face.hpp` (`split_face_curved`, replacing the
  `loop_area == 0.0 || sphere_slit` block)
- `tests/unit/brep_trim_face_test.cpp` (public split of the offset torus face;
  stored-hole `context_b`)

**Interfaces**
- Consumes: `detail::trim_face_domain(brep, face, projection, tolerance)`,
  `detail::trim_face_seam_only_screen`, `detail::trim_face_split_seam_only`,
  `detail::trim_face_split_domain`.
- Produces: same public `split_face_curved` signature; synthetic torus faces with
  cuts now reach `trim_face_split_domain`.

**Step 1: failing test**

Add to `tests/unit/brep_trim_face_test.cpp`, using the same torus and cuts as the
Task 1 test but calling the public API:

```cpp
TEST_CASE("trim_face: the public split of an offset torus band is one ring") {
  const geom::Brep torus = geom::make_torus(
      geom::Point3d{0.0, 0.0, 0.0}, geom::Vec3d{0.0, 0.0, 1.0}, 1.0, 0.5,
      1e-12).value();
  const geom::FaceId tf = face_with_surface<geom::Torus3>(torus);
  std::vector<geom::Point2d> low;
  std::vector<geom::Point2d> high;
  for (int i = 0; i < 64; ++i) {
    const double u = 2.0 * kPi * i / 64.0;
    low.push_back(geom::Point2d{u, kPi / 6.0});
    high.push_back(geom::Point2d{u, 5.0 * kPi / 6.0});
  }
  low.push_back(low.front());
  high.push_back(high.front());
  const std::vector<geom::FaceCut> cuts = {geom::FaceCut{low},
                                           geom::FaceCut{high}};
  const auto fragments = geom::split_face_curved(torus, tf, cuts, 1e-12);
  REQUIRE(fragments.has_value());
  REQUIRE(fragments.value().size() == 3u);
  for (const geom::CurvedFragment& fragment : fragments.value()) {
    REQUIRE(fragment.loops_uv.size() == 1u);
  }
  CHECK(signed_area2(fragments.value().at(0).loops_uv.front()) ==
        doctest::Approx(3.289868133696453).epsilon(1e-9));
  CHECK(signed_area2(fragments.value().at(1).loops_uv.front()) ==
        doctest::Approx(13.159472534785811).epsilon(1e-9));
  CHECK(signed_area2(fragments.value().at(2).loops_uv.front()) ==
        doctest::Approx(23.029076935899213).epsilon(1e-9));
}
```

**Step 2: run and confirm failure**

With HEAD's gate the public path returns `fragments=8 kept=1` (or an
`assemble_curved_solid` self-touch when used inside a boolean); the test must fail
on `size() == 3`.

**Step 3: minimal implementation**

Replace the gate in `split_face_curved`:

```cpp
const Result<detail::TrimFaceDomain> split_domain =
    detail::trim_face_domain(brep, face, projection, tolerance);
if (split_domain.has_value()) {
  const detail::TrimFaceDomain& domain = split_domain.value();
  std::vector<std::vector<Point2d>> holes;
  if (loops.size() > 1) holes.assign(loops.begin() + 1, loops.end());
  if (holes.empty() && domain.collapsed_v()) {
    return detail::trim_face_split_seam_only(face, projection, cuts, domain.bounds);
  }
  if (cuts.empty() && holes.empty() && !domain.collapsed_v()) {
    // whole torus face: fall through to the ordinary path
  } else {
    if (domain.collapsed_v()) {
      const Result<void> screened = detail::trim_face_seam_only_screen(
          face, projection, cuts, domain.bounds);
      if (!screened.has_value()) return screened.error();
    }
    return detail::trim_face_split_domain(face, projection, cuts, domain, holes,
                                          tolerance);
  }
}
```

Notes: keep `trim_face_split_seam_only` for spheres so its historical messages and
1e-12 splitter stay; the no-cut torus fall-through is required because a synthetic
arrangement with no cuts has no material cycle.

**Step 4: verify**

- New test passes with the three areas; the public result equals the direct
  `trim_face_split_domain` result one-for-one (same ring counts and areas).
- A torus face with no cuts returns one unchanged fragment (add the assertion to the
  new test or a small second case).
- `trim face: a stored inner loop is refused publicly` gains
  `CHECK(error.context_b == 1u)` and still returns
  `unsupported "split_face_curved: stored holes are not supported yet"`.
- Sphere tests and `curved boolean: a sphere band joins its two cut runs` unchanged.
- Debug and Release full suites green; report the exact counts.

**Step 5: commit**

`feat(trim_face): dispatch synthetic torus faces to the domain splitter`

## Task 3: Accept horizontal offset torus planes in all three operations

**Files**
- `tests/unit/brep_boolean_curved_test.cpp` (new case before the offset-refusal
  case; keep the existing regressions)

**Interfaces**
- Consumes: `make_torus(center, axis, major, minor, tolerance)`,
  `make_box(lo, hi, tolerance)`, `boolean_solid_curved(brep, solid, brep, solid,
  op, CurvedBooleanOptions{tolerance})`, `validate_topology`, `volume_of`.
- Produces: no production API change; the boolean acceptance for horizontal offset.

**Step 1: failing test**

```cpp
TEST_CASE("curved boolean: an offset torus plane splits both ways") {
  const double tolerance = 1e-9;
  const geom::Result<geom::Brep> torus_result =
      geom::make_torus(geom::Point3d{0.0, 0.0, 0.0}, geom::Vec3d{0.0, 0.0, 1.0},
                       1.0, 0.5, tolerance);
  const geom::Result<geom::Brep> box_result =
      geom::make_box(geom::Point3d{-2.0, -2.0, 0.25},
                     geom::Point3d{2.0, 2.0, 2.0}, tolerance);
  REQUIRE(torus_result.has_value());
  REQUIRE(box_result.has_value());
  const geom::Brep& torus = torus_result.value();
  const geom::Brep& box = box_result.value();

  const auto run = [&](geom::SolidBooleanOp op) {
    return geom::boolean_solid_curved(
        torus, only_solid(torus), box, only_solid(box), op,
        geom::CurvedBooleanOptions{tolerance});
  };

  const auto cut = run(geom::SolidBooleanOp::intersection);
  REQUIRE_MESSAGE(cut.has_value(),
                  (cut.has_value() ? std::string_view{} : message_of(cut.error())));
  CHECK(geom::validate_topology(cut.value().brep).ok());
  CHECK(cut.value().brep.face_ids().size() == 2u);
  CHECK(cut.value().brep.edge_ids().size() == 24u);
  CHECK(cut.value().brep.vertex_ids().size() == 23u);
  CHECK(cut.value().brep.solid_shells(only_solid(cut.value().brep)).size() == 1u);
  CHECK(volume_of(cut.value().brep) ==
        doctest::Approx(0.9647593052603948).epsilon(1e-6));
  const auto cut_again = run(geom::SolidBooleanOp::intersection);
  REQUIRE(cut_again.has_value());
  CHECK(volume_of(cut_again.value().brep) == volume_of(cut.value().brep));

  const auto joined = run(geom::SolidBooleanOp::union_);
  REQUIRE_MESSAGE(joined.has_value(),
                  (joined.has_value() ? std::string_view{} : message_of(joined.error())));
  CHECK(geom::validate_topology(joined.value().brep).ok());
  CHECK(joined.value().brep.face_ids().size() == 9u);
  CHECK(joined.value().brep.edge_ids().size() == 122u);
  CHECK(joined.value().brep.vertex_ids().size() == 116u);
  CHECK(volume_of(joined.value().brep) ==
        doctest::Approx(31.970042895284287).epsilon(1e-6));

  const auto removed = run(geom::SolidBooleanOp::difference);
  REQUIRE_MESSAGE(removed.has_value(),
                  (removed.has_value() ? std::string_view{} : message_of(removed.error())));
  CHECK(geom::validate_topology(removed.value().brep).ok());
  CHECK(removed.value().brep.face_ids().size() == 3u);
  CHECK(removed.value().brep.edge_ids().size() == 110u);
  CHECK(removed.value().brep.vertex_ids().size() == 108u);
  CHECK(volume_of(removed.value().brep) ==
        doctest::Approx(3.970042895284284).epsilon(1e-6));
}
```

**Step 2: run and confirm failure**

The test must fail before Task 1/2 (typically `result touches itself along an edge`
for union, or wrong counts); after Task 1/2 it passes.

**Step 3: no production change required**

If the test fails for a legitimate reason, stop and report the exact failing
assertion and the splitter/assembler state rather than weakening the expected
values. Do not remove the torus plane guard, and do not change
`boolean_solid_curved.hpp`.

**Step 4: verify**

- The new case passes; run it twice and assert the same volume both times (add a
  second `run()` and compare bitwise if cheap).
- All centred torus tests and all six guard refusal tests unchanged.
- Debug and Release suites green; report exact counts.

**Step 5: commit**

`test(boolean): accept an offset torus-plane cut in all three operations`

## Task 4: Pin the refusals, the contract, and the non-goal

**Files**
- `tests/unit/brep_trim_face_test.cpp`
- `tests/unit/brep_boolean_curved_test.cpp`
- `docs/adr/0032-general-domain-splitter.md`
- `docs/benchmarks/m3c3c2b.md` (new)
- `docs/architecture/geometry-library-architecture.md` §16.3
- this plan's exit checklist

**Interfaces**
- Consumes: the Task 1–3 APIs and messages.
- Produces: refusal tests and documentation only; no production change.

**Step 1: failing tests**

The tilted-plane refusal already exists at
`tests/unit/brep_boolean_curved_test.cpp` (`curved boolean: torus cut by a tilted
plane is rejected`, using `make_plane_slab({1,0,0.2}, tolerance)`), and the
axis-parallel offset refusal exists at `curved boolean: torus cut by an offset
parallel plane is rejected`; both already assert `Errc::unsupported`, the exact
message, and `report.face_pairs == 0`. Do not duplicate them.

Add only the missing coverage:

- `trim face: a stored inner loop is refused publicly`: exact message and
  `context_b == 1u`.
- `trim_face: a closed cut crossing the u seam at two different v values is a
  twisted crossing`: polygon `{(1.0,1.0),(6.0,1.1),(0.2,1.2),(3.0,1.15),(1.0,1.0)}`
  with exact message `"split_face_curved: closed cut crosses the seam at different
  transverse coordinates; twisted seam crossings are deferred to M3c-3c-4"`.
- `trim_face: a cut on the parameter domain boundary is refused`: exact message
  `"split_face_curved: cut lies on the parameter domain boundary"`.
- If the axis-parallel offset case is not already covered, add one via
  `make_plane_slab({0,0,1}, tolerance)` translated off the torus centre; expect the
  same torus-plane guard message.
- `trim face: a torus face with no cuts is unchanged` (if not already covered in
  Task 2).

**Step 2: run and confirm**

All refusal cases fail (or are missing) before the assertions are added, then pass.

**Step 3: documentation**

Update `docs/adr/0032-general-domain-splitter.md`:

- one whole-cycle ring contract for periodic-v bands (with the sequence and
  direction rules);
- dispatch by `trim_face_domain` (sphere keeps the seam-only wrapper);
- horizontal offset planes accepted; tilted / axis-parallel offset / twisted /
  stored-hole / domain-boundary still refused;
- tessellation of trimmed periodic-v faces still unsupported, with the measured
  `tessellate: torus face is missing a seam` / `trimmed band wall is not closed`
  evidence;
- deferred: multi-cycle holes, cone apex, general twisted crossings.

Add `docs/benchmarks/m3c3c2b.md` with the measured volumes/counts, tolerance,
compiler, and the note that no timing gate applies. Update §16.3 status and the
assertion counts after the final full-suite run.

**Step 4: verify**

- All new refusals pass with exact messages; no silent success.
- Debug and Release full suites green; record exact counts in ADR/bench/§16.3.
- `git diff` contains no debug output; plan checklist updated.

**Step 5: commit**

`test(trim_face): pin the periodic-domain refusals and the M3c-3c-2b contract`

## Milestone Exit Criteria

1. Offset torus ∩/∪/− with box `z≥0.25` produce the measured volumes, counts,
   single shell, `validate_topology().ok()`, deterministic repeat.
2. All centred torus booleans and all six guard refusals unchanged.
3. Torus bands carry one ring each with positive `2π·height` area; sphere bands
   still one ring; no zero-area ring anywhere.
4. Public split of the offset torus equals the direct splitter result.
5. Stored holes, nested/mixed islands, twisted crossings, domain-boundary cuts,
   tilted and axis-parallel offset planes all refused with exact messages.
6. Debug/Release full suites green with `/W4 /WX /utf-8`, zero warnings; CI six
   jobs green.
7. ADR/bench/§16.3 state the horizontal-only and tessellation-deferred limits.

## Self-check

- [ ] Every task has a failing test and captured RED evidence.
- [ ] No assembler change; no guard change; no sphere behaviour change.
- [ ] Ring directions and area values are exact, not approximate tolerances beyond
      1e-9.
- [ ] Volumes/counts come from the analyst probe and are re-measured, not copied
      blindly.
- [ ] Test files registered (they already are); no new CMake entries needed.
- [ ] Header edits force a rebuild (touch or `/t:Rebuild`) before trusting results.
- [ ] No placeholders, TODOs, or self-correction notes remain in this plan.
- [ ] Documentation counts match the final full-suite output.
