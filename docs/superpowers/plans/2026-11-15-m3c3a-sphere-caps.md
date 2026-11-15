# M3c-3a: Sphere Caps Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use subagent-driven-development (recommended) or executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Let `boolean_solid_curved` cut `Sphere3` operands along closed latitude rings, producing validated cap, complement-cap and slab-band solids with exact-circle healing.

**Architecture:** The seam-only splitter already cuts a sphere latitude correctly, but the curved-boolean clip discards every cut because it clips against the sphere's zero-area seam loop instead of the synthesized parameter rectangle. The plan exposes that synthesized domain, clips against it, then represents pole-touching fragments with one exact pole uv point and joins the two latitude runs of a band into one ring. The sphere gate is then lifted and three end-to-end cases pin the volumes.

**Tech Stack:** C++20, CMake 4.4.4, VS2022 (MSVC 19.44) and GCC/Clang in CI, doctest 2.4.11, `/W4 /WX /utf-8`.

**Spec:** `docs/architecture/geometry-library-architecture.md` §8.12 (curved booleans), §16.3 (milestone status); probe evidence `.superpowers/sdd/m3c3-probe-report.md` (read-only analysis, 556 lines).

## Global Constraints

- Build both `Debug` and `Release` with `-DGEOM_WERROR=ON`; zero warnings on MSVC `/W4 /WX`.
- Every new test is registered in `tests/CMakeLists.txt` (`GEOM_TEST_SOURCES`) in the same commit that adds it.
- Exact predicates decide topology; floating point only samples and rejects. No global epsilon; tolerances are parameters.
- `CurvedFragment.loops_uv[l].size()` must equal `CurvedFragment.loops[l].size()` at all times (`assemble_curved.hpp:469,485`).
- Sphere pole points are exactly `v = -pi/2` or `v = +pi/2` (no projected approximations): `curved_volume.hpp:462-463` tests equality.
- Existing torus/cylinder/plane curved-boolean behaviour must stay byte-for-byte identical; gate every new branch on the synthetic-domain flag.
- Do not weaken or delete an existing assertion to make a new case pass; report the conflict instead.
- One commit per task, message `feat(m3c3a): ...` / `test(m3c3a): ...` / `docs(m3c3a): ...`.

## Scope and deviations (need approval)

This plan is **M3c-3a** of the M3c-3 milestone. It covers sphere operands with closed cuts only (latitude rings and the band between two of them). The following stay deferred and must be recorded in ADR-0029:

1. **sphere/sphere lens (d = 1.5)** — needs cross-chain joining of two seam contacts and an exact sphere/sphere circle; deferred to M3c-3b (probe §D).
2. **cone operands** — apex projection and closed/open cut topology; deferred to M3c-3b.
3. **offset/tilted torus cut planes** and exact torus intersection — deferred to M3c-3b.
4. **trimmed tessellation of boolean-result faces** (holes/non-rectilinear boundaries) — deferred to M3c-3c.
5. **NURBS operands** — deferred to M3c-3c.

The interior island cut (a closed ring that never touches the seam) remains rejected with the existing `"closed cut does not reach the seam"` message; the early rejections at `trim_face.hpp:1023-1058` are correct and must not be deleted.

## Review Focus

Inputs and failures a reasonable user would expect to work, most likely first:

1. A sphere cut by a plane above its centre (`sphere ∩ box z >= 0.5`): the cap must have volume `0.6544984694978736`, be a single closed shell and validate.
2. A sphere with its cap removed (`sphere - box`): volume `3.534291735288517`, and the south pole must appear as exactly one vertex, not a seam row.
3. A sphere sliced by a slab (`|z| <= 0.8`): the band is one fragment with one ring, volume `3.9542179533184135`, and the two latitude rings must join through the seam columns.
4. A sphere cut exactly at the equator (`z >= 0`, half sphere): volume `2π/3`; the cut ring is a full great circle.
5. A sphere face that is genuinely trimmed (a lune from `split_face_curved` with a non-synthetic loop): the synthetic rectangle must **not** be used, and clipping must keep the old behaviour.
6. A cut ring at `v = 0.5` that is sampled with an unwrapped seam end (`u = 2π`): the cap ring must be orientable and non-zero in the chart.

## File Map

- `include/geom/brep/boolean_solid_curved.hpp` — face-info synthetic flag/boundary (T1), region loops for clipping (T2), sphere gate (T5), doc comment on the deferral list.
- `include/geom/brep/trim_face.hpp` — pole point insertion and band run join (T3/T4), public doc comment update (T5).
- `tests/unit/brep_face_uv_test.cpp` — synthetic-domain white-box test (T1).
- `tests/unit/brep_boolean_curved_test.cpp` — clip test (T2), pole/band splitter tests (T3/T4), end-to-end sphere cases and QA table (T5), QA-15 (T6).
- `docs/adr/0029-sphere-cap-booleans.md` — new ADR (T6).
- `docs/benchmarks/m3c3a.md` — new QA-15 report (T6).
- `bench/bench_m3c3.cpp`, `CMakeLists.txt`, `.github/workflows/ci.yml` — new benchmark target and CI build (T6).
- `docs/architecture/geometry-library-architecture.md` §16.3 — M3c-3a status (T6).

---

### Task 1: Expose the synthesized parameter domain on the face info

**Files:**
- Modify: `include/geom/brep/boolean_solid_curved.hpp:209-217` (struct), `:264-272` (substitution)
- Test: `tests/unit/brep_face_uv_test.cpp`

**Interfaces:**
- Consumes: `detail::curved_boolean_face_info(const Brep&, FaceId, double) -> Result<CurvedBooleanFaceInfo>`.
- Produces: `CurvedBooleanFaceInfo::synthetic_domain` (`bool`) and `CurvedBooleanFaceInfo::synthetic_boundary` (`std::array<double,4>` = `{u_lo, v_lo, u_hi, v_hi}`); both default to `false` / zeros.

- [ ] **Step 1: Write the failing test**

Append to `tests/unit/brep_face_uv_test.cpp`, after the existing torus cases:

```cpp
TEST_CASE("curved boolean: the sphere face exposes a synthetic boundary") {
  const double tolerance = 1e-9;
  const geom::Result<geom::Brep> sphere =
      geom::make_sphere(geom::Point3d{0.0, 0.0, 0.0}, 1.0, tolerance);
  REQUIRE(sphere.has_value());
  const geom::Result<geom::detail::CurvedBooleanFaceInfo> sphere_info =
      geom::detail::curved_boolean_face_info(sphere.value(),
                                             sphere.value().face_ids().front(), tolerance);
  REQUIRE(sphere_info.has_value());
  CHECK(sphere_info.value().synthetic_domain);
  CHECK(sphere_info.value().u[0] == 0.0);
  CHECK(sphere_info.value().u[1] == doctest::Approx(2.0 * 3.14159265358979323846));
  CHECK(sphere_info.value().synthetic_boundary[0] == 0.0);
  CHECK(sphere_info.value().synthetic_boundary[1] ==
        doctest::Approx(-0.5 * 3.14159265358979323846));
  CHECK(sphere_info.value().synthetic_boundary[2] ==
        doctest::Approx(2.0 * 3.14159265358979323846));
  CHECK(sphere_info.value().synthetic_boundary[3] ==
        doctest::Approx(0.5 * 3.14159265358979323846));

  const geom::Result<geom::Brep> box = geom::make_box(
      geom::Point3d{-1.0, -1.0, -1.0}, geom::Point3d{1.0, 1.0, 1.0}, tolerance);
  REQUIRE(box.has_value());
  const geom::Result<geom::detail::CurvedBooleanFaceInfo> box_info =
      geom::detail::curved_boolean_face_info(box.value(), box.value().face_ids().front(),
                                             tolerance);
  REQUIRE(box_info.has_value());
  CHECK_FALSE(box_info.value().synthetic_domain);
}
```

- [ ] **Step 2: Run the test to verify it fails**

Run: `.\build\tests\Debug\geom_tests.exe --test-case="*synthetic boundary*"`
Expected: compile error `C2039: 'synthetic_domain' is not a member of 'geom::detail::CurvedBooleanFaceInfo'`.

- [ ] **Step 3: Add the fields**

In `CurvedBooleanFaceInfo` add, after `std::vector<FaceUvLoop> loops;`:

```cpp
  // True when `u`/`v` and `bounds` describe the synthesized full parameter
  // domain of a seam-only face (a sphere or torus whose stored loop is the
  // zero-area seam slit), not the stored loops. Clipping and classification
  // must then use `synthetic_boundary` as the region, because the stored loop
  // encloses no area.
  bool synthetic_domain = false;
  // {u_lo, v_lo, u_hi, v_hi} of that domain; zeros when the flag is false.
  std::array<double, 4> synthetic_boundary{};
```

Include `<array>` if it is not already included in the header.

In `curved_boolean_face_info`, inside the existing substitution branch at `:264-272`, set:

```cpp
    info.synthetic_domain = true;
    info.synthetic_boundary = {u_lo, v_lo, u_hi, v_hi};
```

immediately after the four bounds are replaced and before the 5x5 bounds loop (or immediately after it — the flag only describes the domain; do not set it when `analytic_domain` is nullopt or the collapsed-span predicate is false).

- [ ] **Step 4: Run the test to verify it passes**

Run: `.\build\tests\Debug\geom_tests.exe --test-case="*synthetic boundary*"`
Expected: `1 case, 8 assertions, all passed`.

- [ ] **Step 5: Run the full suite**

Run: `.\build\tests\Debug\geom_tests.exe`
Expected: `933 cases / 659233 assertions` (one new case, +8 assertions) with `0 failed`.

- [ ] **Step 6: Commit**

```bash
git add include/geom/brep/boolean_solid_curved.hpp tests/unit/brep_face_uv_test.cpp
git commit -m "feat(m3c3a): expose the synthesized seam-only domain"
```

---

### Task 2: Clip a seam-only face against its synthetic rectangle

**Files:**
- Modify: `include/geom/brep/boolean_solid_curved.hpp:1146-1156`
- Test: `tests/unit/brep_boolean_curved_test.cpp` (next to the sphere deferral case at `:556`)

**Interfaces:**
- Consumes: `CurvedBooleanFaceInfo::synthetic_domain` / `synthetic_boundary`; `detail::curved_boolean_clip_pair(const CurvedBooleanFaceInfo&, const std::vector<Point2d>&, const CurvedBooleanFaceInfo&, const CircleCurve3*, double) -> std::vector<std::vector<Point2d>>`.
- Produces: no public API change; own/foreign region loops come from `curved_boolean_region_loops(const CurvedBooleanFaceInfo&) -> std::vector<std::vector<Point2d>>`.

- [ ] **Step 1: Write the failing test**

Add to `tests/unit/brep_boolean_curved_test.cpp` right before the sphere deferral case:

```cpp
TEST_CASE("curved boolean: clipping keeps a sphere latitude cut") {
  const double tolerance = 1e-9;
  const geom::Result<geom::Brep> sphere =
      geom::make_sphere(geom::Point3d{0.0, 0.0, 0.0}, 1.0, tolerance);
  const geom::Result<geom::Brep> box = geom::make_box(
      geom::Point3d{-2.0, -2.0, 0.5}, geom::Point3d{2.0, 2.0, 2.0}, tolerance);
  REQUIRE(sphere.has_value());
  REQUIRE(box.has_value());
  const auto sphere_info = geom::detail::curved_boolean_face_info(
      sphere.value(), sphere.value().face_ids().front(), tolerance);
  const auto box_info = geom::detail::curved_boolean_face_info(
      box.value(), box.value().face_ids().front(), tolerance);
  REQUIRE(sphere_info.has_value());
  REQUIRE(box_info.has_value());

  std::vector<geom::Point2d> chain;
  constexpr std::size_t kSamples = 32;
  for (std::size_t i = 0; i <= kSamples; ++i) {
    const double u = 2.0 * 3.14159265358979323846 *
                     static_cast<double>(i) / static_cast<double>(kSamples);
    chain.push_back(geom::Point2d{u, std::asin(0.5)});
  }
  const std::vector<std::vector<geom::Point2d>> pieces =
      geom::detail::curved_boolean_clip_pair(sphere_info.value(), chain,
                                             box_info.value(), nullptr, tolerance);
  REQUIRE(pieces.size() == 1);
  CHECK(pieces.front().size() == chain.size());
}
```

- [ ] **Step 2: Run the test to verify it fails**

Run: `.\build\tests\Debug\geom_tests.exe --test-case="*clipping keeps a sphere*"`
Expected: `pieces.size() == 0`; `REQUIRE(pieces.size() == 1)` fails.

- [ ] **Step 3: Build region loops from the synthetic rectangle**

Add next to `curved_boolean_clip_pair`:

```cpp
// The region of a face in its own parameter chart. A seam-only face stores a
// zero-area seam slit, so membership must be tested against the synthesized
// rectangle instead (the caller's cut can only cross that rectangle).
inline std::vector<std::vector<Point2d>> curved_boolean_region_loops(
    const CurvedBooleanFaceInfo& info) {
  std::vector<std::vector<Point2d>> loops;
  if (info.synthetic_domain) {
    const double u_lo = info.synthetic_boundary[0];
    const double v_lo = info.synthetic_boundary[1];
    const double u_hi = info.synthetic_boundary[2];
    const double v_hi = info.synthetic_boundary[3];
    loops.push_back({Point2d{u_lo, v_lo}, Point2d{u_hi, v_lo}, Point2d{u_hi, v_hi},
                     Point2d{u_lo, v_hi}});
    return loops;
  }
  loops.reserve(info.loops.size());
  for (const FaceUvLoop& loop : info.loops) loops.push_back(loop.points);
  return loops;
}
```

Replace the two loop constructions at `:1151-1156` with:

```cpp
  const std::vector<std::vector<Point2d>> loops = curved_boolean_region_loops(face);
  const std::vector<std::vector<Point2d>> other_loops =
      curved_boolean_region_loops(other);
```

(`trim_face_uv_inside` accepts the same `std::vector<std::vector<Point2d>>` shape it already received.)

- [ ] **Step 4: Run the test to verify it passes**

Run: `.\build\tests\Debug\geom_tests.exe --test-case="*clipping keeps a sphere*"`
Expected: `1 case, 3 assertions, all passed`.

- [ ] **Step 5: Regression-test the non-synthetic curved faces**

Run: `.\build\tests\Debug\geom_tests.exe --test-case="*curved boolean*"`
Expected: all cases pass, including the torus/cylinder suite (`>= 58 cases`).

- [ ] **Step 6: Commit**

```bash
git add include/geom/brep/boolean_solid_curved.hpp tests/unit/brep_boolean_curved_test.cpp
git commit -m "feat(m3c3a): clip seam-only faces against the synthetic rectangle"
```

---

### Task 3: Represent a pole-touching fragment with one exact pole point

**Files:**
- Modify: `include/geom/brep/trim_face.hpp` (helper near `trim_face_seam_only_synthetic_edge` at `:882`; ring builder `:1149-1196`)
- Test: `tests/unit/brep_boolean_curved_test.cpp`

**Interfaces:**
- Consumes: `trim_face_seam_only_domain(const FaceUvProjection&) -> std::optional<std::array<double,4>>`.
- Produces: `bool trim_face_seam_only_pole_v(const std::array<double,4>& domain, double& pole_v) noexcept` — returns true and sets `pole_v` to `domain[2]` or `domain[3]` when the domain is the sphere chart (`v = ±pi/2`); no public API change.

- [ ] **Step 1: Write the failing test**

Add to `tests/unit/brep_boolean_curved_test.cpp`:

```cpp
TEST_CASE("curved boolean: a sphere cap keeps one exact pole point") {
  const double tolerance = 1e-9;
  const geom::Result<geom::Brep> sphere =
      geom::make_sphere(geom::Point3d{0.0, 0.0, 0.0}, 1.0, tolerance);
  REQUIRE(sphere.has_value());
  const geom::FaceId face = sphere.value().face_ids().front();

  std::vector<geom::Point2d> ring;
  constexpr std::size_t kSamples = 64;
  const double v = std::asin(0.5);
  for (std::size_t i = 0; i <= kSamples; ++i) {
    const double u = 2.0 * 3.14159265358979323846 *
                     static_cast<double>(i) / static_cast<double>(kSamples);
    ring.push_back(geom::Point2d{u, v});
  }
  const std::vector<geom::FaceCut> cuts{geom::FaceCut{ring}};
  const geom::Result<std::vector<geom::CurvedFragment>> fragments =
      geom::split_face_curved(sphere.value(), face, cuts, tolerance);
  REQUIRE(fragments.has_value());
  REQUIRE(fragments.value().size() == 2);

  const auto uv_area = [](const std::vector<geom::Point2d>& ring) {
    double twice = 0.0;
    for (std::size_t i = 0; i < ring.size(); ++i) {
      const geom::Point2d& a = ring[i];
      const geom::Point2d& b = ring[(i + 1) % ring.size()];
      twice += a[0] * b[1] - b[0] * a[1];
    }
    return std::fabs(twice) * 0.5;
  };
  bool found_cap = false;
  for (const geom::CurvedFragment& fragment : fragments.value()) {
    REQUIRE(fragment.loops.size() == fragment.loops_uv.size());
    for (std::size_t l = 0; l < fragment.loops_uv.size(); ++l) {
      CHECK(fragment.loops[l].size() == fragment.loops_uv[l].size());
      std::size_t poles = 0;
      for (const geom::Point2d& uv : fragment.loops_uv[l]) {
        if (uv[1] == 0.5 * 3.14159265358979323846) ++poles;
      }
      if (poles == 1) {
        found_cap = true;
        CHECK(uv_area(fragment.loops_uv[l]) > 0.1);
      }
    }
  }
  CHECK(found_cap);
}
```

- [ ] **Step 2: Run the test to verify it fails**

Run: `.\build\tests\Debug\geom_tests.exe --test-case="*one exact pole*"`
Expected: `found_cap` is false (the pole row is synthetic and never enters a ring) and the cap ring area is 0, so at least one CHECK fails.

- [ ] **Step 3: Insert the pole point while the ring is built**

Add the helper:

```cpp
// The sphere chart is the only seam-only chart whose v range is exactly
// [-pi/2, pi/2]; the cone's is [0, r/tan(theta)] and the torus's [0, 2pi].
inline bool trim_face_seam_only_pole_v(const std::array<double, 4>& domain,
                                       double& pole_v) noexcept {
  constexpr double kHalfPi = 1.57079632679489661923;
  if (domain[2] == -kHalfPi) {
    pole_v = domain[2];
    return true;
  }
  if (domain[3] == kHalfPi) {
    pole_v = domain[3];
    return true;
  }
  return false;
}
```

In `trim_face_split_seam_only`, while the cycle flags are built (`:1126-1137`), also record whether the cycle has a synthetic edge on the top side (`from[1] == domain[3] && to[1] == domain[3]`) or the bottom side (`from[1] == domain[2] && to[1] == domain[2]`). Declare before the loop:

```cpp
    double pole_v = 0.0;
    bool touches_pole = false;
```

and inside the `synthetic[i] = 1;` branch:

```cpp
        if (from[1] == domain[3] && to[1] == domain[3]) {
          pole_v = domain[3];
          touches_pole = true;
        } else if (from[1] == domain[2] && to[1] == domain[2]) {
          pole_v = domain[2];
          touches_pole = true;
        }
```

After the run polylines are built and before the closure block at `:1169`, when `touches_pole` is true, append the single pole point to the run and close the ring through it:

```cpp
        points.push_back(Point2d{points.back()[0], pole_v});
        points.push_back(points.front());
        if (points.size() < 3) continue;
        trim_face_normalize_ring(points, true);
        rings.push_back(SeamRing{std::move(points), false});
        continue;
```

This bypasses the seam-pair closure for that run, because the ring now closes through the real 3D meridian chord rim -> pole -> rim. `pole_v` comes from the domain, so it is exactly `±pi/2`.

- [ ] **Step 4: Run the test to verify it passes**

Run: `.\build\tests\Debug\geom_tests.exe --test-case="*one exact pole*"`
Expected: `2 fragments`; the cap ring contains exactly one `v == pi/2` point and its |uv area| is `2pi(pi/2 - asin 0.5)`.

- [ ] **Step 5: Run the seam-only regression cases**

Run: `.\build\tests\Debug\geom_tests.exe --test-case="*seam*"`
Expected: all pass (torus seam, closed-cut rejection, drift).

- [ ] **Step 6: Commit**

```bash
git add include/geom/brep/trim_face.hpp tests/unit/brep_boolean_curved_test.cpp
git commit -m "feat(m3c3a): keep one exact pole point on a cap fragment"
```

---

### Task 4: Join the two cut runs of a pole-free band

**Files:**
- Modify: `include/geom/brep/trim_face.hpp:1198-1216` (after rings are built, before the sort)
- Test: `tests/unit/brep_boolean_curved_test.cpp`

**Interfaces:**
- Consumes: the `rings`/`SeamRing` structure of Task 3; the ring's uv points (a local shoelace in the test).
- Produces: no public API change; a pole-free sphere fragment bounded by two latitude rings is emitted with one ring.

- [ ] **Step 1: Write the failing test**

Add:

```cpp
TEST_CASE("curved boolean: a sphere band joins its two cut runs") {
  const double tolerance = 1e-9;
  const geom::Result<geom::Brep> sphere =
      geom::make_sphere(geom::Point3d{0.0, 0.0, 0.0}, 1.0, tolerance);
  REQUIRE(sphere.has_value());
  const geom::FaceId face = sphere.value().face_ids().front();
  constexpr std::size_t kSamples = 64;
  const double v = std::asin(0.8);

  std::vector<geom::FaceCut> cuts;
  for (const double latitude : {-v, v}) {
    std::vector<geom::Point2d> ring;
    for (std::size_t i = 0; i <= kSamples; ++i) {
      const double u = 2.0 * 3.14159265358979323846 *
                       static_cast<double>(i) / static_cast<double>(kSamples);
      ring.push_back(geom::Point2d{u, latitude});
    }
    cuts.push_back(geom::FaceCut{ring});
  }
  const geom::Result<std::vector<geom::CurvedFragment>> fragments =
      geom::split_face_curved(sphere.value(), face, cuts, tolerance);
  REQUIRE(fragments.has_value());
  REQUIRE(fragments.value().size() == 3);

  const auto uv_area = [](const std::vector<geom::Point2d>& ring) {
    double twice = 0.0;
    for (std::size_t i = 0; i < ring.size(); ++i) {
      const geom::Point2d& a = ring[i];
      const geom::Point2d& b = ring[(i + 1) % ring.size()];
      twice += a[0] * b[1] - b[0] * a[1];
    }
    return std::fabs(twice) * 0.5;
  };
  bool found_band = false;
  for (const geom::CurvedFragment& fragment : fragments.value()) {
    if (fragment.loops_uv.size() != 1) continue;
    double v_min = 1e9;
    double v_max = -1e9;
    for (const geom::Point2d& uv : fragment.loops_uv.front()) {
      if (uv[1] < v_min) v_min = uv[1];
      if (uv[1] > v_max) v_max = uv[1];
    }
    if (v_min < -0.5 && v_max > 0.5) {
      found_band = true;
      CHECK(uv_area(fragment.loops_uv.front()) ==
            doctest::Approx(2.0 * 3.14159265358979323846 * 2.0 * v).epsilon(1e-9));
    }
  }
  CHECK(found_band);
}
```

- [ ] **Step 2: Run the test to verify it fails**

Run: `.\build\tests\Debug\geom_tests.exe --test-case="*band joins*"`
Expected: the band fragment has two zero-area rings, so `fragments.size() == 3` holds but `found_band` is false.

- [ ] **Step 3: Concatenate the rings when there is no pole contact**

After the `rings` vector is complete and before the sort at `:1208`, add:

```cpp
    // A pole-free sphere fragment is bounded by two latitude rings. In the
    // synthesized chart those rings are two horizontal lines; the u = 0 and
    // u = 2pi columns that connect them are synthetic and were dropped above.
    // Join the rings in cycle order through those columns so the fragment has
    // one ring with a real chart interior. A cap touches a pole instead and is
    // never joined.
    double pole_v = 0.0;
    if (rings.size() > 1 && trim_face_seam_only_pole_v(domain, pole_v)) {
      bool touches_pole = false;
      for (const SeamRing& ring : rings) {
        if (ring.seam) continue;
        for (const Point2d& point : ring.points) {
          if (point[1] == pole_v) {
            touches_pole = true;
            break;
          }
        }
        if (touches_pole) break;
      }
      if (!touches_pole) {
        std::vector<Point2d> merged;
        for (const SeamRing& ring : rings) {
          if (ring.seam) continue;
          for (std::size_t i = 0; i + 1 < ring.points.size(); ++i) {
            if (!merged.empty() && trim_face_same_xy(merged.back(), ring.points[i])) {
              continue;
            }
            merged.push_back(ring.points[i]);
          }
        }
        if (!merged.empty()) {
          merged.push_back(merged.front());
          trim_face_normalize_ring(merged, true);
          rings.clear();
          rings.push_back(SeamRing{std::move(merged), false});
        }
      }
    }
```

- [ ] **Step 4: Run the test to verify it passes**

Run: `.\build\tests\Debug\geom_tests.exe --test-case="*band joins*"`
Expected: `3 fragments`; one fragment has one ring with |uv area| `2pi * 2 * asin 0.8`.

- [ ] **Step 5: Run the full suite**

Run: `.\build\tests\Debug\geom_tests.exe`
Expected: `0 failed` (the new cases add 2 cases; existing torus/cylinder cases unchanged).

- [ ] **Step 6: Commit**

```bash
git add include/geom/brep/trim_face.hpp tests/unit/brep_boolean_curved_test.cpp
git commit -m "feat(m3c3a): join the latitude runs of a pole-free band"
```

---

### Task 5: Accept sphere operands and pin the three end-to-end volumes

**Files:**
- Modify: `include/geom/brep/boolean_solid_curved.hpp:172-181` (gate), `include/geom/brep/trim_face.hpp:1285-1307` (doc comment)
- Test: `tests/unit/brep_boolean_curved_test.cpp:556-578` (replace deferral case), QA table `:1489-1491`

**Interfaces:**
- Consumes: all prior tasks; `boolean_solid_curved(Brep, SolidId, Brep, SolidId, SolidBooleanOp, const CurvedBooleanOptions&, CurvedBooleanReport*)`.
- Produces: `Sphere3` accepted by the gate; cone still returns the M3c-3 diagnostic.

- [ ] **Step 1: Replace the deferral test with the cap cases**

Delete `TEST_CASE("curved boolean: sphere operands are rejected with the cap diagnostic")` and add:

```cpp
TEST_CASE("curved boolean: sphere cap intersection") {
  const double tolerance = 1e-9;
  const geom::Result<geom::Brep> sphere =
      geom::make_sphere(geom::Point3d{0.0, 0.0, 0.0}, 1.0, tolerance);
  const geom::Result<geom::Brep> box = geom::make_box(
      geom::Point3d{-2.0, -2.0, 0.5}, geom::Point3d{2.0, 2.0, 2.0}, tolerance);
  REQUIRE(sphere.has_value());
  REQUIRE(box.has_value());
  const geom::SolidId a = only_solid(sphere.value());
  const geom::SolidId b = only_solid(box.value());
  geom::CurvedBooleanReport report;
  const geom::Result<geom::CurvedBooleanResult> result = geom::boolean_solid_curved(
      sphere.value(), a, box.value(), b, geom::SolidBooleanOp::intersection,
      geom::CurvedBooleanOptions{tolerance}, &report);
  REQUIRE_MESSAGE(result.has_value(), "boolean failed: " << message_of(result.error()));
  CHECK(geom::validate_topology(result.value().brep).ok());
  CHECK(report.cuts >= 1);
  CHECK(report.fragments >= 2);
  CHECK(std::fabs(volume_of(result.value().brep) - 0.6544984694978736) <= 1e-6);
}

TEST_CASE("curved boolean: sphere complement cap") {
  const double tolerance = 1e-9;
  const geom::Result<geom::Brep> sphere =
      geom::make_sphere(geom::Point3d{0.0, 0.0, 0.0}, 1.0, tolerance);
  const geom::Result<geom::Brep> box = geom::make_box(
      geom::Point3d{-2.0, -2.0, 0.5}, geom::Point3d{2.0, 2.0, 2.0}, tolerance);
  REQUIRE(sphere.has_value());
  REQUIRE(box.has_value());
  const geom::Result<geom::CurvedBooleanResult> result = geom::boolean_solid_curved(
      sphere.value(), only_solid(sphere.value()), box.value(), only_solid(box.value()),
      geom::SolidBooleanOp::difference, geom::CurvedBooleanOptions{tolerance});
  REQUIRE_MESSAGE(result.has_value(), "boolean failed: " << message_of(result.error()));
  CHECK(geom::validate_topology(result.value().brep).ok());
  CHECK(std::fabs(volume_of(result.value().brep) - 3.534291735288517) <= 1e-6);
}

TEST_CASE("curved boolean: sphere slab band") {
  const double tolerance = 1e-9;
  const geom::Result<geom::Brep> sphere =
      geom::make_sphere(geom::Point3d{0.0, 0.0, 0.0}, 1.0, tolerance);
  const geom::Result<geom::Brep> slab = geom::make_box(
      geom::Point3d{-2.0, -2.0, -0.8}, geom::Point3d{2.0, 2.0, 0.8}, tolerance);
  REQUIRE(sphere.has_value());
  REQUIRE(slab.has_value());
  const geom::Result<geom::CurvedBooleanResult> result = geom::boolean_solid_curved(
      sphere.value(), only_solid(sphere.value()), slab.value(), only_solid(slab.value()),
      geom::SolidBooleanOp::intersection, geom::CurvedBooleanOptions{tolerance});
  REQUIRE_MESSAGE(result.has_value(), "boolean failed: " << message_of(result.error()));
  CHECK(geom::validate_topology(result.value().brep).ok());
  CHECK(result.value().brep.face_ids().size() == 3);
  CHECK(result.value().brep.edge_ids().size() == 3);
  CHECK(std::fabs(volume_of(result.value().brep) - 3.9542179533184135) <= 1e-6);
}
```

- [ ] **Step 2: Run the tests to verify they fail**

Run: `.\build\tests\Debug\geom_tests.exe --test-case="*sphere cap*" --test-case="*sphere complement*" --test-case="*sphere slab*"`
Expected: `Errc::unsupported` with `"sphere operands need cap topology; booleans are deferred to M3c-3"` (gate still closed).

- [ ] **Step 3: Lift the sphere gate**

In `curved_boolean_surface_gate`, replace the `Sphere3` rejection branch with:

```cpp
  if (std::get_if<Sphere3>(&surface) != nullptr) {
    return CurvedBooleanSurfaceGate{true, nullptr};
  }
```

Keep the cone branch and its message unchanged.

Update the `split_face_curved` doc comment at `trim_face.hpp:1285-1307`: a closed cut on a seam-only sphere face may be a latitude (one seam contact), a cap fragment touching a pole is emitted with exactly one pole point (`v = ±pi/2`), and a pole-free band is emitted as one ring joined through the seam columns. Remove the sentence claiming an interior closed cut is rejected; the island rejection stays.

- [ ] **Step 4: Run the three tests to verify they pass**

Run: `.\build\tests\Debug\geom_tests.exe --test-case="*sphere cap*" --test-case="*sphere complement*" --test-case="*sphere slab*"`
Expected: 3 cases pass; volumes within `1e-6`.

If `sphere complement cap` reports `"fragment representative lies on the other solid"`, do **not** weaken the volume assertion: record the failure with the report counters and stop (the representative fallback is out of this task's scope; it would be a plan defect to fix here).

- [ ] **Step 5: Update the QA corpus row**

In the QA corpus at `:1489-1491`, replace the sphere row that asserts the deferral message with one that runs `sphere ∩ box z >= 0.5` and checks volume `0.6544984694978736` with `<= 1e-6` relative error and `validate_topology(...).ok()`. Leave the cone deferral row and the torus row untouched.

- [ ] **Step 6: Run the full suite in both configurations**

Run:
```
.\build\tests\Debug\geom_tests.exe
.\build-rel\tests\Release\geom_tests.exe
```
Expected: both green. The case count is the baseline 933 plus the four earlier new cases (T1-T4) plus this task's net `+2` (one deferral case replaced by three): `939`. Record the exact counts in the task report; do not adjust an assertion to match a different number.

- [ ] **Step 7: Commit**

```bash
git add include/geom/brep/boolean_solid_curved.hpp include/geom/brep/trim_face.hpp tests/unit/brep_boolean_curved_test.cpp
git commit -m "feat(m3c3a): support sphere operands on closed cuts"
```

---

### Task 6: QA-15, ADR-0029, bench and status docs

**Files:**
- Create: `docs/adr/0029-sphere-cap-booleans.md`
- Create: `docs/benchmarks/m3c3a.md`
- Create: `bench/bench_m3c3.cpp`
- Modify: `CMakeLists.txt` (bench target), `.github/workflows/ci.yml:113` (bench job), `docs/architecture/geometry-library-architecture.md` §16.3
- Test: QA-15 cases in `tests/unit/brep_boolean_curved_test.cpp`

**Interfaces:**
- Consumes: `boolean_solid_curved`, `volume_of`, `check_qa_volume` (existing helpers).
- Produces: QA-15 corpus rows, ADR-0029, `geom_bench_m3c3` target.

- [ ] **Step 1: Write the QA-15 test rows**

Add a `TEST_CASE("QA-15: sphere cap booleans are closed and deterministic")` that runs the three cases from Task 5 twice each, asserts `check_qa_volume`-style agreement (relative `<= 1e-6`), `validate_topology(...).ok()`, single shell and closed mesh (`curved_boolean_mesh_closure(...).closed()`), and that the second run's volume equals the first exactly.

- [ ] **Step 2: Run it**

Run: `.\build\tests\Debug\geom_tests.exe --test-case="QA-15*"`
Expected: pass. If the mesh closure is not closed for one case, report which case rather than relaxing the closure check.

- [ ] **Step 3: Write ADR-0029 and the benchmark**

`docs/adr/0029-sphere-cap-booleans.md` records: the synthetic-domain exposure, clip against the rectangle, the single exact pole point, the band run join, the deferred items (lens, cone, offset torus, trimmed result tessellation, NURBS), and why island cuts stay rejected.

`bench/bench_m3c3.cpp` mirrors `bench/bench_m3c2.cpp` (100 iterations, median, volume column, checksum, exit 2 on checksum mismatch) with rows: cap intersection, complement cap, slab band; reference volumes as in Task 5. Register `geom_bench_m3c3` in `CMakeLists.txt` under `GEOM_BUILD_BENCH` and add it to the CI bench `--target` list.

- [ ] **Step 4: Run the bench**

Run (Release): `.\build-rel\bench\Release\geom_bench_m3c3.exe`
Expected: medians printed in ms, checksum matches, exit 0.

- [ ] **Step 5: Update status docs**

`docs/benchmarks/m3c3a.md` records the environment, the three volumes, timings, and the deferred list. Architecture §16.3 adds an M3c-3a row: completed items and the M3c-3b/3c deferred list. Do not claim whole M3c-3 complete.

- [ ] **Step 6: Commit**

```bash
git add docs/adr/0029-sphere-cap-booleans.md docs/benchmarks/m3c3a.md bench/bench_m3c3.cpp CMakeLists.txt .github/workflows/ci.yml docs/architecture/geometry-library-architecture.md tests/unit/brep_boolean_curved_test.cpp
git commit -m "docs(m3c3a): record QA-15, the ADR and the benchmark"
```

---

## M3c-3a exit criteria

- Debug and Release both green with zero warnings (`-DGEOM_WERROR=ON`).
- `sphere ∩ box z >= 0.5` = `0.6544984694978736`, `sphere - box` = `3.534291735288517`, `sphere ∩ slab |z| <= 0.8` = `3.9542179533184135`, each within `1e-6` relative, validated, single closed shell.
- Torus/cylinder/plane curved booleans unchanged (all prior tests pass without edits).
- Cone still returns `"boolean_solid_curved: cone operands need closed-cut and apex topology; booleans are deferred to M3c-3"`.
- QA-15 test, ADR-0029, bench target and §16.3 updated; CI bench job builds `geom_bench_m3c3`.

## Self-review

- **Spec coverage:** synthetic domain (T1/T2), pole representation (T3), band join (T4), gate + volumes (T5), QA/docs (T6). The deferred list is explicit in Scope and ADR-0029.
- **Step scan:** each task has a failing test, an exact failing output, a minimal implementation, a passing output, and a commit. No placeholder steps.
- **Type consistency:** `synthetic_domain`/`synthetic_boundary`, `curved_boolean_region_loops`, `trim_face_seam_only_pole_v`, and the three volume constants are used with the same names and types throughout.
- **Review Focus:** equator cut, genuine trimmed face, and unwrapped seam end are covered by T3/T5 assertions and the existing torus suite; the island cut keeps its rejection test from M3c-2.
- **Proportion:** five production tasks (~160 LOC) plus one docs/QA task; test code carries the exact values that the report verified.
