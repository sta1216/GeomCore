# M3c-3c-2a — General Domain Splitter for Sphere and Torus

## Goal

Generalize the seam-only splitter into one synthetic-domain splitter for **closed cuts on sphere charts**, and remove the silent losses around it: a synthetic chart that already carries a stored inner ring is refused explicitly instead of dropping the hole, and the public tangent-loop rejection is pinned. Torus charts keep their current plane-only contract; their periodic-v bands need a two-rim fragment contract and rim pairing in the assembler, so offset/tilted torus plane booleans are deferred to M3c-3c-2b (ADR-0032). Cone apex and open cut chains stay out of scope.

## Architecture

The splitter stops deciding "is this a seam-only face?" from the area of `loops.front()`. Instead a domain decision function classifies each of the four uv rectangle sides as `kSeamCopy`, `kCollapsed`, or `kRealBoundary`, and the same decision is shared by `split_face_curved` and `boolean_solid_curved`. Cut chains are windowed onto the periodic domain once, then fed to the existing `split_segments` + `Arrangement2` machinery. Material faces are extracted from the arrangement; runs along seam copies are dropped, runs along collapsed sides collapse to the single physical point, runs along real boundaries stay as geometry. A face that already carries stored inner loops is refused with a documented diagnostic: the arrangement exposes one cycle per face, so hole preservation needs multi-cycle walking and is deferred to M3c-3c-4.

## Tech Stack

C++20; CMake 4.4.4 + VS2022 MSVC 19.44 (`/W4 /WX /utf-8`); doctest 2.4.11; existing `Arrangement2`, `split_segments`, `face_uv_projection`, `assemble_curved`. No new dependencies.

## Spec

- `docs/architecture/geometry-library-architecture.md` §8.2 (trimmed curved faces), §8.10 (tessellation contract), §16.3 (status).
- Probe report `D:\ai_geometry\.superpowers\sdd\m3c3c2-probe-report.md` §1–§6.
- Existing tests `tests/unit/brep_trim_face_test.cpp`, `tests/unit/brep_boolean_curved_test.cpp`.

## Global Constraints

- No global epsilon; all closure/seam tests use relative bounds (`1e-12 * max(1, |a|, |b|)`), matching `trim_face_same_xy` / `trim_face_seam_only_seam_pair`.
- Topology decisions only through exact predicates (`sign_orient2d`, `ring_orientation2`), never raw cross products.
- `split_face_curved` never returns an empty vector: no cut geometry consumed is `Errc::numeric_failure "split_face_curved: no fragments"`.
- Deterministic fragment order: seam rings first, then `|uv area|` descending, then point count, then ring key.
- Every task adds its tests to `tests/CMakeLists.txt` if it creates a new file; this plan only extends two existing test files, so no CMake change is expected.
- MSBuild header dependency tracking is unreliable in this repo: after touching a header, rebuild with `/t:Rebuild` or touch a `.cpp`.
- Debug and Release suites must both pass with zero warnings; final counts are recorded in the plan, ADR, benchmark and §16.3.

## Scope and Deviations (need approval)

1. **Closed cuts only.** Open chains (seam→pole, base rim→apex) and the cone apex are M3c-3c-2b. This milestone rejects an open chain with the existing message `"split_face_curved: open cuts on a seam-only face"` unless a task says otherwise.
2. **Torus operands keep their current plane-cut contract.** Torus/torus and torus/other-curved operands stay rejected with the M3c-2 message, and the "perpendicular to or containing the axis" guard stays. Offset/tilted torus plane booleans are deferred to M3c-3c-2b: their bands wrap the periodic-v seam and need a two-rim fragment contract plus rim pairing in the assembler. Both mechanisms tried here (shifting the v origin onto a cut; emitting full-width v-side rings) failed the offset acceptance and regressed five centered-plane torus booleans (ADR-0032).
3. **Stored holes are refused, not silently dropped.** A synthetic-domain face carrying `loops_uv` inner rings stops with `unsupported "split_face_curved: stored holes are not supported yet"`. Preserving them needs every boundary cycle of an arrangement face, which is deferred to M3c-3c-4 together with the ear-clip cost.
4. **Ear-clipping cost is not fixed here.** 4096-point rims still cost ~3.7 s per `triangulate_polygon2`; QA-18a asserts fragment and loop structure rather than paying that cost inside the suite. The cost fix is a separate follow-up.

## Review Focus

1. **Silent loss on synthetic charts.** A synthetic sphere chart that carries a stored inner ring used to drop the hole; the public `split_face_curved` must now return the exact `unsupported` message. The torus public path keeps its M3c-2 plane-cut contract in this milestone; the direct-API fragment oddity is documented, not fixed. Probe both the refusal and the unchanged torus guard.
2. **Seam identity.** For every seam ring, the repeated first point and the point before it must be the two parameter copies of one physical point within the relative bound.
3. **Orientation.** Ring direction is preserved from the arrangement traversal. The two sides of an interior island must come out with opposite traversal; normalizing both to CCW would make two faces claim the same material side. Check coedge direction after assembly.
4. **Area partition.** `sum(fragment_area) == whole-face area` within tolerance; UV and 3D areas for sphere caps; torus area `4π²Rr` for halves.
5. **No regression.** Existing sphere lens QA-16/17, sphere caps QA-15, sphere slit tests, planar trim tests, and torus QA rows must stay green. Torus gate message changes are the only expected test edits.

## File Map

| Path | Purpose |
| --- | --- |
| `include/geom/brep/trim_domain.hpp` (new) | `TrimSideKind`, `TrimFaceDomain`, `trim_face_domain`, `trim_face_window_chains`; seam windowing helpers moved out of `boolean_solid_curved.hpp` |
| `include/geom/brep/trim_face.hpp` | `trim_face_split_domain`; `trim_face_split_seam_only` becomes a wrapper; dispatch by domain; non-empty holes refused |
| `include/geom/brep/boolean_solid_curved.hpp` | include `trim_domain.hpp`, delete duplicated static helpers, torus gate relaxation |
| `tests/unit/brep_trim_face_test.cpp` | domain/windowing/splitter tests (hole refusal, not preservation) |
| `tests/unit/brep_boolean_curved_test.cpp` | torus offset/tilted plane e2e and gate-message updates |
| `tests/unit/curved_boolean_qa18_test.cpp` (new) | QA-18a sphere latitude-split corpus |
| `docs/adr/0032-general-domain-splitter.md`, `docs/benchmarks/m3c3c2a.md`, architecture §16.3 | decisions and measurements |

## Task 1 — Domain descriptor and cut-chain windowing

**Files:** create `include/geom/brep/trim_domain.hpp`; modify `include/geom/brep/boolean_solid_curved.hpp` (include it and delete `curved_boolean_split_at_seam`, `curved_boolean_window_chains`, `curved_boolean_cut_chains` duplicates), `include/geom/brep/trim_face.hpp` (include it); `tests/unit/brep_trim_face_test.cpp`.

**Interfaces (produces):**

```cpp
namespace geom::detail {

enum class TrimSideKind { kSeamCopy, kCollapsed, kRealBoundary };

struct TrimFaceDomain {
  std::array<double, 4> bounds{};        // u_lo, u_hi, v_lo, v_hi
  bool periodic_u = false;
  bool periodic_v = false;
  std::array<TrimSideKind, 4> sides{};   // u_lo, u_hi, v_lo, v_hi
  bool collapsed_v() const;              // any v side is kCollapsed
};

Result<TrimFaceDomain> trim_face_domain(
    const Brep& brep, FaceId face, const FaceUvProjection& projection,
    double tolerance);

// Closed chains stay closed; a chain crossing the seam is split into two
// chains ending exactly on the sides. An open chain whose endpoints are not on
// a side or a collapsed side is rejected.
Result<std::vector<std::vector<Point2d>>> trim_face_window_chains(
    std::span<const FaceCut> cuts, const TrimFaceDomain& domain,
    double tolerance);

}  // namespace geom::detail
```

Domain values: Sphere3 `{0, 2π, -π/2, π/2}`, `periodic_u = true`, sides `{kSeamCopy, kSeamCopy, kCollapsed, kCollapsed}`; Torus3 `{0, 2π, 0, 2π}`, both periodic, all four `kSeamCopy`. Anything else returns `Errc::unsupported "trim_face: no synthetic domain for this surface"`.

**Steps:**

1. Write the failing tests:

```cpp
TEST_CASE("trim_face: sphere and torus synthetic domains") {
  const geom::Brep sphere =
      geom::make_sphere(geom::Point3d{0.0, 0.0, 0.0}, 1.0, 1e-12).value();
  const geom::FaceId sf = face_with_surface<geom::Sphere3>(sphere);
  const geom::detail::FaceUvProjection sp =
      geom::detail::face_uv_projection(sphere, sphere.face(sf), 1e-12).value();
  const auto sd = geom::detail::trim_face_domain(sphere, sf, sp, 1e-12);
  REQUIRE(sd.has_value());
  CHECK(sd->bounds[0] == 0.0);
  CHECK(sd->bounds[1] == doctest::Approx(2.0 * kPi));
  CHECK(sd->bounds[2] == doctest::Approx(-kPi / 2.0));
  CHECK(sd->bounds[3] == doctest::Approx(kPi / 2.0));
  CHECK(sd->periodic_u);
  CHECK_FALSE(sd->periodic_v);
  CHECK(sd->sides[0] == geom::detail::TrimSideKind::kSeamCopy);
  CHECK(sd->sides[2] == geom::detail::TrimSideKind::kCollapsed);
  CHECK(sd->sides[3] == geom::detail::TrimSideKind::kCollapsed);

  const geom::Brep torus = geom::make_torus(
      geom::Point3d{0.0, 0.0, 0.0}, geom::Vec3d{0.0, 0.0, 1.0}, 1.0, 0.5,
      1e-12).value();
  const geom::FaceId tf = face_with_surface<geom::Torus3>(torus);
  const geom::detail::FaceUvProjection tp =
      geom::detail::face_uv_projection(torus, torus.face(tf), 1e-12).value();
  const auto td = geom::detail::trim_face_domain(torus, tf, tp, 1e-12);
  REQUIRE(td.has_value());
  CHECK(td->periodic_u);
  CHECK(td->periodic_v);
  for (int i = 0; i < 4; ++i) {
    CHECK(td->sides[i] == geom::detail::TrimSideKind::kSeamCopy);
  }
}

TEST_CASE("trim_face: windowing splits a chain at the periodic seam") {
  const geom::Brep sphere =
      geom::make_sphere(geom::Point3d{0.0, 0.0, 0.0}, 1.0, 1e-12).value();
  const geom::FaceId sf = face_with_surface<geom::Sphere3>(sphere);
  const geom::detail::FaceUvProjection sp =
      geom::detail::face_uv_projection(sphere, sphere.face(sf), 1e-12).value();
  const auto sd = geom::detail::trim_face_domain(sphere, sf, sp, 1e-12);
  REQUIRE(sd.has_value());
  const std::vector<geom::FaceCut> cuts = {
      geom::FaceCut{{geom::Point2d{0.2, 0.3}, geom::Point2d{3.0, 0.5},
                     geom::Point2d{6.2, 0.4}, geom::Point2d{0.3, 0.6}}}};
  const auto chains =
      geom::detail::trim_face_window_chains(cuts, *sd, 1e-12);
  REQUIRE(chains.has_value());
  REQUIRE(chains->size() == 2u);
  CHECK((*chains)[0].back()[0] == doctest::Approx(2.0 * kPi));
  CHECK((*chains)[1].front()[0] == 0.0);
  CHECK((*chains)[0].size() + (*chains)[1].size() >= 5u);
}

TEST_CASE("trim_face: an interior dangling cut endpoint is rejected") {
  const geom::Brep torus = geom::make_torus(
      geom::Point3d{0.0, 0.0, 0.0}, geom::Vec3d{0.0, 0.0, 1.0}, 1.0, 0.5,
      1e-12).value();
  const geom::FaceId tf = face_with_surface<geom::Torus3>(torus);
  const geom::detail::FaceUvProjection tp =
      geom::detail::face_uv_projection(torus, torus.face(tf), 1e-12).value();
  const auto td = geom::detail::trim_face_domain(torus, tf, tp, 1e-12);
  REQUIRE(td.has_value());
  const std::vector<geom::FaceCut> cuts = {
      geom::FaceCut{{geom::Point2d{1.0, 0.4}, geom::Point2d{2.0, 0.4}}}};
  const auto chains =
      geom::detail::trim_face_window_chains(cuts, *td, 1e-12);
  REQUIRE_FALSE(chains.has_value());
  CHECK(chains.error().code == geom::Errc::unsupported);
  CHECK(chains.error().message ==
        std::string_view("trim_face: cut chain ends inside the parameter "
                         "domain"));
}
```

2. Run: `.\build\tests\Debug\geom_tests.exe --test-case="trim_face:*"` — expect compile failure (missing header).
3. Implement the descriptor and windowing by moving the seam helpers from `boolean_solid_curved.hpp`; keep the old helper names as `inline` aliases only if a caller outside this plan uses them (grep first). Windowing preserves point order and inserts the seam point exactly at `u_lo`/`u_hi`.
3b. Seam splitting must be direction-aware: a chain dominated by decreasing u still has to split when it wraps from the low side to the high side, while a long step that stays inside the window must not be split. `trim_face_dominant_step_sign` derives the traversal direction from the short steps, and `trim_face_split_at_seam(..., dominant_sign)` only splits wraps against it. Closure and side tests use `max(tolerance, 1e-12 * scale)` so a documented tolerance is honoured.
4. Rebuild Debug/Release and run the cases below plus `curved boolean:*`.
5. Commit: `feat(m3c3c2a): describe synthetic domains and window cut chains`.

**Acceptance:** the descriptor, windowing, direction, tolerance and rejection cases pass; no existing test changes; counts increase by 7 cases / roughly 40 assertions. The committed cases pin: sphere/torus domain values, forward seam split, backward (decreasing) seam split, periodic-v split, closed chain preservation, side chain acceptance, unsupported surface, tolerance-widened closure, and the interior dangling-endpoint rejection.

## Task 2 — Arrangement splitter for closed cuts on a synthetic domain

**Files:** `include/geom/brep/trim_face.hpp`; `tests/unit/brep_trim_face_test.cpp`.

**Interfaces (produces):**

```cpp
namespace geom::detail {
Result<std::vector<CurvedFragment>> trim_face_split_domain(
    FaceId face, const FaceUvProjection& projection,
    std::span<const FaceCut> cuts, const TrimFaceDomain& domain,
    std::span<const std::vector<Point2d>> holes, double tolerance);
}
```

`holes` are stored inner uv rings; a non-empty span is refused with `unsupported "split_face_curved: stored holes are not supported yet"` — preserving them is M3c-3c-4.

`trim_face_split_seam_only` keeps its signature and delegates: build the sphere domain, window, split. Existing sphere-slit callers see the same fragments.

**Algorithm decisions:**

- Build the four synthetic sides from the domain; each side carries a `TrimSideKind` index so a cycle edge can be classified without geometric tests.
- Window every cut chain first (Task 1). Closed chains become one source each. A chain that crosses the seam is split into pieces ending on the seam sides; the pieces are joined only when every junction pairs the two seam copies at the **same** transverse coordinate (`v` for a `u` seam, `u` for a `v` seam), because that is the only case where translating the tail by one period leaves the shoelace area unchanged. A closed cut whose two seam crossings have different transverse coordinates is refused with `unsupported "split_face_curved: closed cut crosses the seam at different transverse coordinates; twisted seam crossings are deferred to M3c-3c-4"` instead of producing a chord or a wrong ring. After the fragments are emitted, verify that their absolute uv areas sum to the chart area within `1e-9 * chart_area`; otherwise return `numeric_failure "split_face_curved: fragments do not partition the chart"` (this is the overlapping-island guard).
- Run `split_segments` + `Arrangement2` exactly as today. For each bounded face with `winding(f, 0) != 0` (source 0 = synthetic sides), split its cycle into runs:
  - `kSeamCopy` run: drop it; the ring is closed through the seam pair when the run's neighbours are cut runs on opposite sides.
  - `kCollapsed` run: replace the run by the single collapse point (pole) at the side's v.
  - `kRealBoundary` run: keep the run as ring geometry.
- Preserve traversal direction; do not normalize a ring to CCW. The two faces of an interior island are visited separately and come out with opposite orientation.
- A non-empty `holes` span is refused before any arrangement work. Otherwise, if no fragment is produced, return `Errc::numeric_failure "split_face_curved: no fragments"` (never empty).
- Ring order and fragment order follow the existing `ring_key` rules.

**Steps:**

1. Failing tests:

```cpp
TEST_CASE("trim_face: a closed island splits into two fragments") {
  const geom::Brep sphere =
      geom::make_sphere(geom::Point3d{0.0, 0.0, 0.0}, 1.0, 1e-12).value();
  const geom::FaceId sf = face_with_surface<geom::Sphere3>(sphere);
  std::vector<geom::Point2d> ring;
  for (int i = 0; i < 33; ++i) {
    const double a = 2.0 * kPi * i / 33.0;
    ring.push_back(geom::Point2d{1.0 + 0.2 * std::cos(a),
                                 0.3 + 0.2 * std::sin(a)});
  }
  ring.push_back(ring.front());
  const geom::detail::FaceUvProjection sp =
      geom::detail::face_uv_projection(sphere, sphere.face(sf), 1e-12).value();
  const auto domain = geom::detail::trim_face_domain(sphere, sf, sp, 1e-12);
  REQUIRE(domain.has_value());
  const std::vector<geom::FaceCut> cuts = {geom::FaceCut{ring}};
  const auto fragments = geom::detail::trim_face_split_domain(
      sf, sp, cuts, domain.value(), std::span<const std::vector<geom::Point2d>>{},
      1e-12);
  REQUIRE(fragments.has_value());
  REQUIRE(fragments.value().size() == 2u);
  CHECK(fragments.value().at(0).loops_uv.size() == 1u);
  CHECK(fragments.value().at(1).loops_uv.size() == 1u);
  const double a0 = geom::signed_area2(fragments.value().at(0).loops_uv[0]);
  const double a1 = geom::signed_area2(fragments.value().at(1).loops_uv[0]);
  // The island is material on one side and the rest of the sphere on the
  // other: the two rings keep opposite traversal and partition the chart.
  CHECK(a0 * a1 < 0.0);
  CHECK(std::fabs(a0) + std::fabs(a1) ==
        doctest::Approx(2.0 * kPi * kPi).epsilon(1e-12));
  CHECK(std::min(std::fabs(a0), std::fabs(a1)) ==
        doctest::Approx(kPi * 0.04).epsilon(1e-9));
}

TEST_CASE("trim_face: torus latitudes split without the loop-area gate") {
  const geom::Brep torus = geom::make_torus(
      geom::Point3d{0.0, 0.0, 0.0}, geom::Vec3d{0.0, 0.0, 1.0}, 1.0, 0.5,
      1e-12).value();
  const geom::FaceId tf = face_with_surface<geom::Torus3>(torus);
  std::vector<geom::Point2d> low;
  std::vector<geom::Point2d> high;
  for (int i = 0; i < 64; ++i) {
    const double u = 2.0 * kPi * i / 64.0;
    low.push_back(geom::Point2d{u, 0.6});
    high.push_back(geom::Point2d{u, 2.5});
  }
  low.push_back(low.front());
  high.push_back(high.front());
  const geom::detail::FaceUvProjection tp =
      geom::detail::face_uv_projection(torus, torus.face(tf), 1e-12).value();
  const auto domain = geom::detail::trim_face_domain(torus, tf, tp, 1e-12);
  REQUIRE(domain.has_value());
  const std::vector<geom::FaceCut> cuts = {geom::FaceCut{low},
                                           geom::FaceCut{high}};
  const auto fragments = geom::detail::trim_face_split_domain(
      tf, tp, cuts, domain.value(), std::span<const std::vector<geom::Point2d>>{},
      1e-12);
  REQUIRE(fragments.has_value());
  REQUIRE(fragments.value().size() == 3u);
  CHECK(fragments.value().at(0).loops_uv.size() == 1u);
  CHECK(fragments.value().at(1).loops_uv.size() == 2u);
  CHECK(fragments.value().at(2).loops_uv.size() == 1u);
}
```

2. Run `--test-case="trim_face:*"`; expect the island to fail with `"closed cut does not reach the seam"` and the torus case to fail with 0 fragments (this is the silent-empty bug).
3. Implement `trim_face_split_domain`; keep the old seam-only tests green (sphere slit, caps, band join).
4. Rebuild and run `trim_face:*`, `curved boolean:*`, `tessellate: trimmed*`.
5. Commit: `feat(m3c3c2a): split synthetic domains with closed cuts`.

**Acceptance:** interior island → 2 fragments whose rings have opposite orientation and whose absolute areas sum to the chart area `2π²` (the smaller is the exact 33-gon area `0.12490582127787074`); torus two latitudes → 3 fragments with ring counts 1 / 2 / 1 (this test is structure-only: the old u-unwrap already accepted it, so it is not RED); a u-seam island with two different transverse crossings → the twisted-crossing `unsupported`; two overlapping islands → `numeric_failure "split_face_curved: fragments do not partition the chart"`; fragment order asserted explicitly; no empty vector; existing sphere slit and band tests unchanged.

## Task 3 — Explicit refusal of stored holes on synthetic charts

**Files:** `include/geom/brep/trim_face.hpp`, `tests/unit/brep_trim_face_test.cpp`.

**Context.** Routing every synthetic domain through `trim_face_split_domain` exposes periodic-v torus bands that today's assembler cannot pair. This milestone keeps the torus path and its plane guard unchanged, and removes the remaining silent hole loss.

**Changes:**

- In the seam-only branch of `split_face_curved`, when the face is a sphere slit and `loops.size() > 1`, return `unsupported "split_face_curved: stored holes are not supported yet"` before the splitter runs.
- Pin the public tangent-loop rejection so a later dispatch cannot hand it to the general splitter.
- Do not change `boolean_solid_curved`.

**Steps:**

1. Failing test: one-face sphere band, outer ring at v=0.3, stored inner ring at v=0.6 (SolidEditor, 64 `LineCurve3` segments per ring), `split_face_curved(band, face, {}, 1e-9)` must return exactly the message above. Today it succeeds and silently drops the inner ring.
2. Pin: tangent closed cut (centre (0.1,0.5), r=0.1, vertex (0,0.5)) through `split_face_curved` stays `unsupported "split_face_curved: closed cut crosses the seam twice"`.
3. Implement the refusal; run `trim*` and `curved boolean:*`; then Debug and Release full suites.
4. Commit `fix(m3c3c2a): refuse stored holes on synthetic charts`.

**Acceptance:** exact stored-hole message; tangent message unchanged; no other test edits; +1 case / +2 assertions.

## Task 4 — QA-18a, ADR-0032 and status

**Files:** `tests/unit/curved_boolean_qa18_test.cpp` (new), `docs/adr/0032-general-domain-splitter.md`, `docs/benchmarks/m3c3c2a.md`, architecture §16.3, this plan.

**Changes:**

- QA-18a (no new algorithm path, so no timing gate):
  1. a whole sphere with two latitude cuts through the public splitter → 3 fragments, each with exactly one uv loop (1/1/1);
  2. the sphere band with a stored inner loop → exact `unsupported` message;
  3. a twisted island and nested islands through the direct API → their exact messages;
  4. the sphere row twice → identical fragment and loop counts.
- `docs/benchmarks/m3c3c2a.md` records that no new benchmark was added; the milestone is infrastructure plus refusals, and existing trimmed/caps timings are unchanged.
- ADR-0032 records the splitter contract, the explicit refusals, and why torus offset/tilted planes are deferred (two-rim bands + assembler pairing; both attempted mechanisms failed and regressed five centered torus booleans). It lists the 3c-2b follow-ups: two-rim contract, rim pairing, offset/tilted planes, twisted crossings, multi-cycle holes.
- Architecture §16.3 gets one paragraph; plan checkboxes are updated.

**Acceptance:** four rows pass in both configurations; no claim that torus offset planes work; +1 case / about 12 assertions.

## Exit Criteria

- `split_face_curved` handles sphere and torus closed cuts through one domain decision; no silent empty vector.
- Torus plane booleans keep their pre-milestone behavior and messages (offset/tilted planes deferred to M3c-3c-2b); torus/torus and torus/curved stay explicitly rejected.
- Stored holes are refused with a documented diagnostic, not silently dropped.
- Debug and Release suites pass, zero warnings; M3c-3c-1/2/3a/3b QA rows unchanged.
- ADR-0032, `docs/benchmarks/m3c3c2a.md`, §16.3 updated.

## Self-Check

- [ ] Every task has Files, Interfaces, RED test, and a commit boundary.
- [ ] The deviations (closed cuts only, horizontal-plane torus, tilted/twisted crossings refused, stored holes refused, nested/mixed closed cuts refused, ear-clip deferred) are called out and will be confirmed with the user before execution.
- [ ] No placeholder or self-correction text.
- [ ] Tests use the existing test-local helpers (`face_with_surface`, `kPi`), `geom::detail::face_uv_projection`, `make_sphere`, `make_torus`, `make_box`, `validate_topology`, `signed_area2`.
- [ ] The plan does not claim cone/apex/open-chain support.

## Execution

Ask the user to approve the four deviations, then run subagent-driven: one implementer and one reviewer per task, ledger at `.superpowers/sdd/m3c3c2a-splitter/progress.md`.
