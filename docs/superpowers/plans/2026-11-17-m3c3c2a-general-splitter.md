# M3c-3c-2a — General Domain Splitter for Sphere and Torus

## Goal

Replace the seam-only special case in `split_face_curved` with one synthetic-domain splitter that works whenever a curved face's uv chart is a parameter rectangle with seam copies and collapsed sides. This milestone covers **closed cuts** on sphere and torus charts, preserves stored holes, fixes the torus silent-empty-success bug, and makes offset/tilted **plane** cuts of a torus produce valid boolean results. Cone apex and open cut chains are deliberately left for M3c-3c-2b.

## Architecture

The splitter stops deciding "is this a seam-only face?" from the area of `loops.front()`. Instead a domain decision function classifies each of the four uv rectangle sides as `kSeamCopy`, `kCollapsed`, or `kRealBoundary`, and the same decision is shared by `split_face_curved` and `boolean_solid_curved`. Cut chains are windowed onto the periodic domain once, then fed to the existing `split_segments` + `Arrangement2` machinery. Material faces are extracted from the arrangement; runs along seam copies are dropped, runs along collapsed sides collapse to the single physical point, runs along real boundaries stay as geometry. Stored inner loops are additional winding sources so a face that already has a hole keeps it.

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
2. **Torus operands remain plane-cut only.** Torus/torus and torus/other-curved operands stay rejected with the M3c-2 message. Offset/tilted planes become supported; the old "perpendicular to or containing the axis" guard is removed.
3. **Holes are preserved through stored-loop winding**, not through a general multi-face trimming engine: a synthetic-domain face may carry `loops_uv` inner rings; the splitter adds them as winding sources and drops fragments inside a hole.
4. **Ear-clipping cost is not fixed here.** 4096-point rims still cost ~3.7 s per `triangulate_polygon2`; QA-18a asserts structure plus analytic/partial volumes rather than paying that cost inside the suite. The cost fix is a separate follow-up.

## Review Focus

1. **Silent empty success.** `split_face_curved` on a torus with cuts `v=0.6` and `v=2.5` currently returns 0 fragments and no error through the public entry; the fix must return 3 fragments, and a genuinely consumed-everything input must return a documented error. Probe both.
2. **Seam identity.** For every seam ring, the repeated first point and the point before it must be the two parameter copies of one physical point within the relative bound.
3. **Orientation.** Ring direction is preserved from the arrangement traversal. The two sides of an interior island must come out with opposite traversal; normalizing both to CCW would make two faces claim the same material side. Check coedge direction after assembly.
4. **Area partition.** `sum(fragment_area) == whole-face area` within tolerance; UV and 3D areas for sphere caps; torus area `4π²Rr` for halves.
5. **No regression.** Existing sphere lens QA-16/17, sphere caps QA-15, sphere slit tests, planar trim tests, and torus QA rows must stay green. Torus gate message changes are the only expected test edits.

## File Map

| Path | Purpose |
| --- | --- |
| `include/geom/brep/trim_domain.hpp` (new) | `TrimSideKind`, `TrimFaceDomain`, `trim_face_domain`, `trim_face_window_chains`; seam windowing helpers moved out of `boolean_solid_curved.hpp` |
| `include/geom/brep/trim_face.hpp` | `trim_face_split_domain`; `trim_face_split_seam_only` becomes a wrapper; dispatch by domain, hole winding |
| `include/geom/brep/boolean_solid_curved.hpp` | include `trim_domain.hpp`, delete duplicated static helpers, torus gate relaxation |
| `tests/unit/brep_trim_face_test.cpp` | domain/windowing/splitter/hole tests |
| `tests/unit/brep_boolean_curved_test.cpp` | torus offset/tilted plane e2e and gate-message updates |
| `tests/unit/curved_boolean_qa18_test.cpp` (new) | QA-18a torus corpus |
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
4. Rebuild Debug/Release and run the three cases plus `curved boolean:*`.
5. Commit: `feat(m3c3c2a): describe synthetic domains and window cut chains`.

**Acceptance:** all three cases pass; no existing test changes; counts increase by 3 cases / 27 assertions.

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

`holes` are stored inner uv rings (the face's material excludes their interiors); Task 2 passes an empty span, Task 3 fills it from `face_loops_uv`.

`trim_face_split_seam_only` keeps its signature and delegates: build the sphere domain, window, split. Existing sphere-slit callers see the same fragments.

**Algorithm decisions:**

- Build the four synthetic sides from the domain; each side carries a `TrimSideKind` index so a cycle edge can be classified without geometric tests.
- Window every cut chain first (Task 1). Closed chains become one source each; a chain that crossed the seam was split into two open chains ending on the seam sides, and the two pieces are joined into one source through the seam column, preserving order.
- Run `split_segments` + `Arrangement2` exactly as today. For each bounded face with `winding(f, 0) != 0` (source 0 = synthetic sides), split its cycle into runs:
  - `kSeamCopy` run: drop it; the ring is closed through the seam pair when the run's neighbours are cut runs on opposite sides.
  - `kCollapsed` run: replace the run by the single collapse point (pole) at the side's v.
  - `kRealBoundary` run: keep the run as ring geometry.
- Preserve traversal direction; do not normalize a ring to CCW. The two faces of an interior island are visited separately and come out with opposite orientation.
- If a face has no cut edge and `holes` is empty, skip it. With holes and no cuts, emit one fragment whose ring is the hole boundary (traversed clockwise) — the seam sides are dropped. If no fragment is produced, return `Errc::numeric_failure "split_face_curved: no fragments"` (never empty).
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
      sf, sp, cuts, *domain, std::span<const std::vector<geom::Point2d>>{},
      1e-12);
  REQUIRE(fragments.has_value());
  REQUIRE(fragments->size() == 2u);
  CHECK(fragments->at(0).loops_uv.size() == 1u);
  CHECK(fragments->at(1).loops_uv.size() == 1u);
  const double a0 = geom::signed_area2(fragments->at(0).loops_uv[0]);
  const double a1 = geom::signed_area2(fragments->at(1).loops_uv[0]);
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
      tf, tp, cuts, *domain, std::span<const std::vector<geom::Point2d>>{},
      1e-12);
  REQUIRE(fragments.has_value());
  REQUIRE(fragments->size() == 3u);
  CHECK(fragments->at(0).loops_uv.size() == 1u);
  CHECK(fragments->at(1).loops_uv.size() == 2u);
  CHECK(fragments->at(2).loops_uv.size() == 1u);
}
```

2. Run `--test-case="trim_face:*"`; expect the island to fail with `"closed cut does not reach the seam"` and the torus case to fail with 0 fragments (this is the silent-empty bug).
3. Implement `trim_face_split_domain`; keep the old seam-only tests green (sphere slit, caps, band join).
4. Rebuild and run `trim_face:*`, `curved boolean:*`, `tessellate: trimmed*`.
5. Commit: `feat(m3c3c2a): split synthetic domains with closed cuts`.

**Acceptance:** island → 2 fragments whose rings have opposite orientation and whose absolute areas sum to the chart area `2π²` (the smaller is `0.04π`); torus → 3 fragments with ring counts 1 / 2 / 1; no empty vector; existing sphere slit and band tests unchanged.

## Task 3 — One domain decision, torus gate, stored holes

**Files:** `include/geom/brep/trim_face.hpp` (dispatch), `include/geom/brep/boolean_solid_curved.hpp` (gate + include), `tests/unit/brep_trim_face_test.cpp`, `tests/unit/brep_boolean_curved_test.cpp`.

**Changes:**

- `split_face_curved` first calls `trim_face_domain`; if it returns a domain, route through `trim_face_split_domain` regardless of `loops.front()` area. Ordinary planar/cylinder faces keep the current path. This removes the torus "loop area 4π² so never synthetic" bug.
- A synthetic-domain face with stored inner loops: each stored loop with `|signed_area2| > tolerance²` becomes an additional closed cut source with negative winding; fragments whose representative lies inside such a loop are dropped. Outer stored loops stay the face boundary.
- `boolean_solid_curved`: remove the torus "cut plane must be perpendicular to or contain the axis" guard (and its helper); keep torus/torus and torus/curved rejection unchanged. Torus plane cuts now flow through the general splitter.
- Replace the torus diagnostic tests that asserted `"torus cut plane must be..."` with acceptance tests for offset/tilted planes.

**Steps:**

1. Failing tests:

```cpp
TEST_CASE("curved boolean: an offset torus plane splits both ways") {
  const geom::Brep torus = geom::make_torus(
      geom::Point3d{0.0, 0.0, 0.0}, geom::Vec3d{0.0, 0.0, 1.0}, 1.0, 0.5,
      1e-9).value();
  const geom::Brep upper = geom::make_box(
      geom::Point3d{-2.0, -2.0, 0.25}, geom::Point3d{2.0, 2.0, 2.0},
      1e-9).value();
  // intersect torus with the half-space box -> the upper band
  const auto result = geom::boolean_solid_curved(
      torus, torus.solid_ids().front(), upper, upper.solid_ids().front(),
      geom::SolidBooleanOp::intersection);
  REQUIRE(result.has_value());
  CHECK(geom::validate_topology(result->brep).ok());
  const auto shells = result->brep.solid_shells(result->brep.solid_ids().front());
  CHECK(shells.size() == 1u);
}

TEST_CASE("curved boolean: torus difference stays deferred") {
  const geom::Brep torus = geom::make_torus(
      geom::Point3d{0.0, 0.0, 0.0}, geom::Vec3d{0.0, 0.0, 1.0}, 1.0, 0.5,
      1e-9).value();
  const geom::Brep upper = geom::make_box(
      geom::Point3d{-2.0, -2.0, 0.25}, geom::Point3d{2.0, 2.0, 2.0},
      1e-9).value();
  const auto result = geom::boolean_solid_curved(
      torus, torus.solid_ids().front(), upper, upper.solid_ids().front(),
      geom::SolidBooleanOp::difference);
  // The B-side cap of a torus would need a flipped torus fragment, which the
  // curved path still rejects; the failure must stay explicit.
  REQUIRE_FALSE(result.has_value());
  CHECK(result.error().code == geom::Errc::unsupported);
}

TEST_CASE("trim_face: a stored hole is not dropped") {
  const geom::Brep sphere =
      geom::make_sphere(geom::Point3d{0.0, 0.0, 0.0}, 1.0, 1e-12).value();
  const geom::FaceId sf = face_with_surface<geom::Sphere3>(sphere);
  const geom::detail::FaceUvProjection sp =
      geom::detail::face_uv_projection(sphere, sphere.face(sf), 1e-12).value();
  const auto domain =
      geom::detail::trim_face_domain(sphere, sf, sp, 1e-12);
  REQUIRE(domain.has_value());
  std::vector<geom::Point2d> hole;
  for (int i = 0; i < 33; ++i) {
    const double a = 2.0 * kPi * i / 33.0;
    hole.push_back(geom::Point2d{1.0 + 0.2 * std::cos(a),
                                 0.3 + 0.2 * std::sin(a)});
  }
  hole.push_back(hole.front());
  const std::vector<std::vector<geom::Point2d>> holes = {hole};
  const auto fragments = geom::detail::trim_face_split_domain(
      sf, sp, std::span<const geom::FaceCut>{}, *domain, holes, 1e-12);
  REQUIRE(fragments.has_value());
  REQUIRE(fragments->size() == 1u);
  CHECK(fragments->at(0).loops_uv.size() == 1u);
  // The material is the sphere minus the disk, so the only ring is the hole
  // boundary; it is traversed clockwise in uv (material outside, hole inner).
  const double area = geom::signed_area2(fragments->at(0).loops_uv[0]);
  CHECK(area == doctest::Approx(-kPi * 0.04).epsilon(1e-9));
}
```

The hole test calls the splitter-level API directly, so it needs no `SolidEditor` fixture. The production path feeds `holes` from `face_loops_uv` inner rings (Task 3).

2. Run; expect the torus case to fail at the old guard message and the hole case to keep the island (wrong).
3. Implement dispatch, hole sources, gate removal. Update the torus diagnostic tests (currently expecting `"torus cut plane must be perpendicular..."` and `"torus/torus..."`); keep torus/torus rejection tests.
4. Run `curved boolean:*`, `trim_face:*`, `*torus*`, `curved volume:*`; Debug/Release.
5. Commit: `feat(m3c3c2a): route torus planes through the general splitter`.

**Acceptance:** offset plane z=0.25 yields one valid shell for intersection and union; torus difference stays explicitly `unsupported` until flipped torus fragments exist (M3c-3c-2b); holes are preserved; torus/torus still rejected with the M3c-2 message.

## Task 4 — QA-18a, ADR, benchmark, status

**Files:** `tests/unit/curved_boolean_qa18_test.cpp` (new), `docs/adr/0032-general-domain-splitter.md` (new), `docs/benchmarks/m3c3c2a.md` (new), `docs/architecture/geometry-library-architecture.md` (§16.3), `cmake`/CI bench target if a new bench is added.

**QA-18a cases** (each: `validate_topology().ok()`, exactly one shell, closed by uses, V−E+F per shell, deterministic repeat; relative `1e-6`):

1. Torus `R=1, r=0.5` intersected with box `z >= 0` → half torus, volume `2.467401100272339`, surface area `9.869604401089358`;
2. Torus rotated 30° about x, intersected with the box `z >= 0` (a plane through the centre at 30° to the torus axis) → two fragments whose volumes are both in `(0, 4.934802200544679)` and sum to `4.934802200544679`; the halves are not asserted equal (a tilted plane is not a symmetry plane);
3. Torus intersected with box `z >= 0.25` → upper band; assert volume in `(0, 4.934802200544679)`, plus `V-E+F = 2`, not an analytic half;
4. Regression: sphere lens QA-17 row (37°) still exact.

Assertion counts are pinned after the first green run; the plan's expected totals are updated then.

**Steps:** write the RED tests (they fail on the current guard/messages), implement only the fixture wiring, run Debug and Release, record wall times, write ADR-0032 (decisions: descriptor, seam windowing, orientation preservation, holes via winding, torus gate scope, deferred open chains/cone), benchmark page, §16.3 entry, then commit `docs(m3c3c2a): record QA-18a and the general splitter limits`.

**Acceptance:** QA-18a green; ADR states the closed-cut-only and plane-only limits; benchmark records medians; §16.3 says M3c-3c-2a complete, 3c-2b/3/4 open.

## Exit Criteria

- `split_face_curved` handles sphere and torus closed cuts through one domain decision; no silent empty vector.
- Torus plane booleans (offset and tilted) produce valid single-shell results; torus/torus still explicitly rejected.
- Stored holes survive a synthetic-domain split.
- Debug and Release suites pass, zero warnings; M3c-3c-1/2/3a/3b QA rows unchanged.
- ADR-0032, `docs/benchmarks/m3c3c2a.md`, §16.3 updated.

## Self-Check

- [ ] Every task has Files, Interfaces, RED test, and a commit boundary.
- [ ] The deviations (closed cuts only, plane-only torus, hole-by-winding, ear-clip deferred) are called out and will be confirmed with the user before execution.
- [ ] No placeholder or self-correction text.
- [ ] Tests use the existing test-local helpers (`face_with_surface`, `kPi`), `geom::detail::face_uv_projection`, `make_sphere`, `make_torus`, `make_box`, `validate_topology`, `signed_area2`.
- [ ] The plan does not claim cone/apex/open-chain support.

## Execution

Ask the user to approve the four deviations, then run subagent-driven: one implementer and one reviewer per task, ledger at `.superpowers/sdd/m3c3c2a-splitter/progress.md`.
