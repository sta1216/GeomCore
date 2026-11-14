# M3c-2 Curved Boolean Completion Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Make `boolean_solid_curved` accept torus operands, so the curved B-rep kernel can intersect, unite and subtract that analytic primitive; sphere and cone operands stay explicitly unsupported with the topology they need deferred to M3c-3.

**Architecture:** The driver already owns pairing, SSI cut extraction, classification, fragment healing and assembly. A cone face whose only real loop is the seam slit has a degenerate uv domain, so SSI never sees a patch; `split_face_curved` cannot cut a face whose outer boundary has zero uv area; `curved_solid_volume` cannot triangulate a face with poles. Torus needs the gate plus v-seam handling. Cone and sphere booleans hit further gaps (clip against the zero-area loop, cap ring pole detour, meridian folded into the cut arc span, and two seam contacts for sphere/sphere), so they are pinned as explicitly unsupported and deferred to M3c-3. The torus gate is opened and pinned by end-to-end volume tests; cone and sphere are pinned by failure diagnostics that name the missing topology.

**Tech Stack:** C++20, CMake, MSVC 19.44 / GCC / Clang, doctest 2.4.11, existing `geom::v1` headers only (no new dependency).

**Spec:** `docs/architecture/geometry-library-architecture.md` v2.1 §8.10 (tessellation contract), §8.11 (mass properties), §8.12 (boolean guarantees), §16.3 (M3 status); `docs/quality/testing-strategy.md`; probe report `.superpowers/sdd/m3c2-probe-report.md`.

---

## Global Constraints

- Public API stays inside `namespace geom::v1`; no new third-party dependency.
- Every new failure is an `Errc` with a message naming the operation and the offending entity; never silently return an empty result.
- Volumes in tests are analytic values; tolerances are relative and stated per test.
- Windows is MSVC 19.44 with `/W4 /WX /utf-8`; the Linux/macOS CI jobs build the same headers with GCC/Clang. All tests must pass on both Debug and Release.
- Every task follows TDD: add the failing test, show the failure, implement the minimum, show both configurations green, commit once.
- Do not change existing planar/cylindrical boolean results; the regression suites in `tests/unit/brep_boolean_curved_test.cpp` and `tests/unit/split_solid_test.cpp` must stay green.

## Scope and deviations (need review approval)

M3c-2 covers **torus operands of `boolean_solid_curved`**, plus the pole-aware sphere volume needed for measurement. Sphere and cone *booleans* are deferred to M3c-3 and pinned by diagnostics. The milestone is deliberately split from general trimmed-face tessellation:

1. **Closed cuts only on seam-only faces, anchored at one seam point.** A cut that ends in the middle of a face still returns `unsupported "split_face_curved: open cuts on a seam-only face"`; an interior closed ring (front and back points coincide) is rejected as `"closed cut crosses the seam twice"`. This is exactly why cone cuts cannot be delivered here: a cone/plane cut is either an interior closed circle or an open meridian chain ending on the base rim (see Task 7). Sphere booleans are out of scope (Task 5).
2. **Synthetic parameter boundary.** For a seam-only face the splitter synthesizes the missing rectangle of the parameter domain (`u ∈ [u0,u0+2π]`, `v ∈ [v_lo,v_hi]`), runs the existing arrangement, then strips every chain that came from the synthetic seam/pole edges and closes the remaining cut chain by 3D endpoint identity. The synthetic edges never appear in `CurvedFragment::loops`, so assembly sees only real cut edges.
3. **Torus cuts stay chordal.** There is no exact torus/plane or torus/cylinder healing (the intersection is quartic in general). Quality is reported by `CurvedBooleanReport::max_cut_residual`; the tests require it below `1e-3` for the chosen corpus.
4. **Cone booleans are deferred.** Every cone/plane cut needs either interior-closed-ring splitting or open-chain-to-real-boundary handling, plus apex welding and an apex volume fan. All of that is M3c-3 work; M3c-2 pins the gate message and zero pairing (Task 7).
4b. **Torus cut planes are restricted to axis-perpendicular or axis-containing planes, compared within a rounding band, with no geometric snapping.** `make_torus` stores the axis unnormalized (only within `kAnalyticSurfaceUnitTolerance = 1e-12` of unit), and even a unit axis can make the normalized dot round to `1 ± 1 ulp`, so exact equality rejected valid cuts. `curved_boolean_torus_plane_cut_supported` (and `curved_boolean_pair_circles`) normalize the torus axis first, then accept `| |n . axis| - 1 | <= 8 ulp` (perpendicular) or `|n . axis| <= 8 ulp` with `|n . center - offset| <= tolerance` (containing). Because the dot is stationary at 1, that band also admits genuine tilts up to about `6e-8` rad; when such a pair is admitted, the exact cut path fails loudly in `trim_face` with `"segment intersection is not representable"` and never produces wrong geometry (a plane `5e-13` short of perpendicular is rejected by the guard; a plane `1e-10` off the axis, if admitted, winds through a degenerate SSI loop), and M3c-2 does not snap, so every plane outside the band is rejected before pairing with `unsupported "boolean_solid_curved: torus cut plane must be perpendicular to or contain the torus axis"`. Offset-parallel (parallel but not through the axis) and tilted cuts are M3c-3.
5. **No trimmed tessellation of boolean results in this milestone.** QA-14 checks topology, closedness, volume and determinism; tessellating a result face with a non-rectilinear boundary is M3c-3 work (M3c-1 explicitly rejects those with `"... boundary is not rectilinear"`).

## Review Focus (feed these to every reviewer)

- Does the sphere/cone/torus uv-domain synthesis ever *widen* a legitimate trimmed domain? A face whose loop-derived area is non-zero must keep exactly its old range.
- Are the synthetic seam/pole chains really absent from every emitted `CurvedFragment` loop, and do the cut chains close with bitwise-identical 3D endpoints?
- Does the pole-aware volume integrate the same surface the tessellator would (no flat-triangle approximation of the sphere)?
- Do torus cuts that cross `v = 0 / 2π` produce loops with no `2π` jump?
- Are all error paths preserving the failure report and returning a message that names the stage?

## File Map

| File | Responsibility |
| --- | --- |
| `include/geom/brep/boolean_solid_curved.hpp` | failure report, uv-domain synthesis, gate widening, seam split axis, sphere-sphere exact circle |
| `include/geom/brep/trim_face.hpp` | seam-only boundary synthesis and synthetic-chain stripping |
| `include/geom/brep/detail/curved_volume.hpp` | pole-aware face triangulation |
| `include/geom/brep/split_solid.hpp` | sphere/sphere exact circle helper (plane construction) |
| `tests/unit/brep_boolean_curved_test.cpp` | end-to-end torus booleans and the sphere/cone diagnostics |
| `tests/unit/brep_trim_face_test.cpp` | seam-only fragment splitting |
| `tests/unit/brep_curved_volume_test.cpp` | pole-aware volume |
| `tests/unit/boolean_solid_curved_report_test.cpp` | failure report |
| `docs/adr/0028-curved-boolean-surfaces.md` | decisions and limits |
| `docs/benchmarks/m3c2.md` | QA-14 numbers |

---

### Task 1: Preserve the diagnostic report on failure

**Files:**
- Modify: `include/geom/brep/boolean_solid_curved.hpp:1430-1833`
- Create: `tests/unit/boolean_solid_curved_report_test.cpp`
- Modify: `tests/CMakeLists.txt` (register the new test)

**Interfaces:**
- Consumes: `CurvedBooleanReport` (`boolean_solid_curved.hpp:104`).
- Produces:

```cpp
inline Result<CurvedBooleanResult> boolean_solid_curved(
    const Brep& a, SolidId a_solid, const Brep& b, SolidId b_solid,
    SolidBooleanOp op, const CurvedBooleanOptions& options = {},
    CurvedBooleanReport* failure_report = nullptr);
```

Semantics: `failure_report`, when non-null, is zero-initialized on entry and filled with the counters accumulated at the point of every early `Error` return. On success it is left untouched (the returned `CurvedBooleanResult::report` is authoritative). Passing `nullptr` is exactly today's behavior.

- [ ] **Step 1: Write the failing test**

```cpp
// tests/unit/boolean_solid_curved_report_test.cpp
#include <doctest/doctest.h>

#include <geom/brep/boolean_solid_curved.hpp>
#include <geom/brep/primitives.hpp>

using namespace geom::v1;

TEST_CASE("curved boolean: a failed operation still reports the stage counters") {
  const Result<Brep> a = make_box(Point3d{0, 0, 0}, Point3d{1, 1, 1}, 1e-9);
  const Result<Brep> b = make_box(Point3d{10, 0, 0}, Point3d{11, 1, 1}, 1e-9);
  REQUIRE(a.has_value());
  REQUIRE(b.has_value());
  const SolidId solid_a = a.value().solid_ids().front();
  const SolidId solid_b = b.value().solid_ids().front();
  CurvedBooleanReport report{};
  report.face_pairs = 999;  // must be overwritten by the operation
  const Result<CurvedBooleanResult> result =
      boolean_solid_curved(a.value(), solid_a, b.value(), solid_b, SolidBooleanOp::intersection,
                           CurvedBooleanOptions{}, &report);
  REQUIRE_FALSE(result.has_value());
  CHECK(result.error().code == Errc::numeric_failure);
  CHECK(report.face_pairs == 0);
  CHECK(report.cuts == 0);
  CHECK(report.kept == 0);
}

TEST_CASE("curved boolean: a failure after pairing copies the accumulated counters") {
  const Result<Brep> a = make_box(Point3d{0, 0, 0}, Point3d{1, 1, 1}, 1e-9);
  const Result<Brep> b = make_box(Point3d{3, 0, 0}, Point3d{4, 1, 1}, 1e-9);
  REQUIRE(a.has_value());
  REQUIRE(b.has_value());
  CurvedBooleanReport report{};
  const Result<CurvedBooleanResult> result = boolean_solid_curved(
      a.value(), a.value().solid_ids().front(), b.value(), b.value().solid_ids().front(),
      SolidBooleanOp::union_, CurvedBooleanOptions{}, &report);
  REQUIRE_FALSE(result.has_value());
  CHECK(result.error().code == Errc::unsupported);
  CHECK(report.fragments == 12u);
  CHECK(report.kept == 12u);
}

TEST_CASE("curved boolean: success leaves the failure report untouched") {
  const Result<Brep> a = make_box(Point3d{0, 0, 0}, Point3d{2, 2, 2}, 1e-9);
  const Result<Brep> b = make_box(Point3d{1, 1, 1}, Point3d{3, 3, 3}, 1e-9);
  REQUIRE(a.has_value());
  REQUIRE(b.has_value());
  CurvedBooleanReport report{};
  const Result<CurvedBooleanResult> result = boolean_solid_curved(
      a.value(), a.value().solid_ids().front(), b.value(), b.value().solid_ids().front(),
      SolidBooleanOp::union_, CurvedBooleanOptions{}, &report);
  REQUIRE(result.has_value());
  CHECK(report.face_pairs == 0u);
  CHECK(report.cuts == 0u);
  CHECK(report.fragments == 0u);
  CHECK(report.kept == 0u);
  CHECK(report.dropped == 0u);
  CHECK(report.flipped == 0u);
  CHECK(report.max_cut_residual == 0.0);
}
```

- [ ] **Step 2: Run it and confirm it fails**

Run: `.\build\tests\Debug\geom_tests.exe --test-case="curved boolean: a failed operation still reports the stage counters"`
Expected: compile error — `boolean_solid_curved` takes 6 arguments.

- [ ] **Step 3: Implement**

Add the parameter, zero it on entry (`*failure_report = CurvedBooleanReport{};`), and replace every `return Error{...}` in the function body with a two-line form that copies the locally accumulated counters (the variables ending up in `report` today) into `*failure_report` before returning. Keep the counters in one `CurvedBooleanReport report{};` local used both for the success return and the failure copy; do not add a second set of counters.

- [ ] **Step 4: Verify**

Run Debug and Release: `.\build\tests\Debug\geom_tests.exe` and `.\build-rel\tests\Release\geom_tests.exe --test-case="*curved*"`.
Expected: all green; the existing curved boolean tests compile unchanged.

- [ ] **Step 5: Commit**

`git add -A; git commit -m "feat(m3c2): report curved boolean counters on failure"`

---

### Task 2: Synthesize the uv domain of a seam-only face

**Files:**
- Modify: `include/geom/brep/boolean_solid_curved.hpp:182-273`
- Test: `tests/unit/brep_face_uv_test.cpp`

**Interfaces:**
- Produces (in `geom::v1::detail`):

```cpp
// Full parameter rectangle of a surface whose face may have no real boundary
// loop (a sphere or torus whose only loop is its seam, or a cone whose lateral
// loop collapses at the apex). Returns {u_lo, u_hi, v_lo, v_hi}; nullopt for
// surfaces whose domain is not known analytically.
inline std::optional<std::array<double, 4>>
curved_boolean_surface_uv_domain(const SurfaceGeometry& surface) noexcept;
```

Values: `Sphere3` → `{0, 2π, -π/2, π/2}`; `Torus3` → `{0, 2π, 0, 2π}`; `Cone3` → `{0, 2π, 0, radius / tan(half_angle)}` (axial, i.e. the primitive height; `Cone3::eval` is axial, not ruling-length); everything else → `nullopt`.

- [ ] **Step 1: Write the failing test**

Add to `tests/unit/brep_face_uv_test.cpp` (include `<geom/surfaces/analytic.hpp>`):

```cpp
TEST_CASE("curved boolean: a sphere face advertises its full uv domain") {
  const Result<Brep> sphere = make_sphere(Point3d{0, 0, 0}, 1.0, 1e-9);
  REQUIRE(sphere.has_value());
  const std::vector<FaceId> faces = sphere.value().face_ids();
  REQUIRE(faces.size() == 1u);
  const Result<detail::CurvedBooleanFaceInfo> info =
      detail::curved_boolean_face_info(sphere.value(), faces.front(), 1e-9);
  REQUIRE(info.has_value());
  CHECK(info.value().u[0] == 0.0);
  CHECK(info.value().u[1] == doctest::Approx(2.0 * 3.14159265358979323846));
  CHECK(info.value().v[0] == doctest::Approx(-1.5707963267948966));
  CHECK(info.value().v[1] == doctest::Approx(1.5707963267948966));
  CHECK(info.value().bounds.lo[0] == doctest::Approx(-1.0));
  CHECK(info.value().bounds.hi[0] == doctest::Approx(1.0));
  CHECK(info.value().bounds.lo[1] == doctest::Approx(-1.0));
  CHECK(info.value().bounds.hi[1] == doctest::Approx(1.0));
  CHECK(info.value().bounds.lo[2] == doctest::Approx(-1.0));
  CHECK(info.value().bounds.hi[2] == doctest::Approx(1.0));
}

TEST_CASE("curved boolean: the torus advertises its full uv domain") {
  const auto domain = detail::curved_boolean_surface_uv_domain(
      Torus3{Point3d{0, 0, 0}, Vec3d{0, 0, 1}, Vec3d{1, 0, 0}, 2.0, 0.5});
  REQUIRE(domain.has_value());
  CHECK((*domain)[0] == 0.0);
  CHECK((*domain)[1] == doctest::Approx(2.0 * 3.14159265358979323846));
  CHECK((*domain)[2] == 0.0);
  CHECK((*domain)[3] == doctest::Approx(2.0 * 3.14159265358979323846));
}

TEST_CASE("curved boolean: the cone domain ends at the axial apex") {
  Cone3 cone;
  cone.origin = Point3d{0, 0, 0};
  cone.axis = Vec3d{0, 0, 1};
  cone.x_axis = Vec3d{1, 0, 0};
  cone.radius = 1.0;
  cone.half_angle = std::atan(0.5);  // height 2
  const auto domain = detail::curved_boolean_surface_uv_domain(cone);
  REQUIRE(domain.has_value());
  CHECK((*domain)[2] == 0.0);
  CHECK((*domain)[3] == doctest::Approx(2.0));
}
```

- [ ] **Step 2: Run it and confirm it fails**

Run: `.\build\tests\Debug\geom_tests.exe --test-case="curved boolean: a sphere face advertises its full uv domain"`.
Expected: `u[1]` is 0 (not 2π) and the bounds collapse to the meridian slice (`y=0`, `x` in `[0,1]`), so the test fails.

- [ ] **Step 3: Implement**

Add `curved_boolean_surface_uv_domain` above `curved_boolean_face_info`. In `curved_boolean_face_info`, after the loop-derived bounds and the periodic clamp, replace all four bounds only when `u_hi - u_lo <= 1e-12 * max(1, max(|u_lo|, |u_hi|))` and the analytic domain exists; recompute the 5×5 chart bounds with the existing `surface_eval` loop. Do not change plane/cylinder behavior.

Add a regression test for the guard: a trimmed torus patch (`R=2`, `r=0.5`, `u in [0, 0.1]`, `v in [0, 0.5]`, four `CircleCurve3` edges) must keep `u = [0, 0.1]`, `v = [0, 0.5]` and a narrow `bounds`; a version that substitutes unconditionally would report the full `4π²` rectangle.

- [ ] **Step 4: Verify**

Run Debug and Release full suites plus `--test-case="*curved*"`. Expected 898+ cases: the sphere/cone/torus/guard cases pass; `boolean_solid_curved` still rejects sphere operands (gate is Task 5).

- [ ] **Step 5: Commit**

`git add include/geom/brep/boolean_solid_curved.hpp tests/unit/brep_face_uv_test.cpp; git commit -m "feat(m3c2): synthesize the full uv domain of seam-only faces"`


---

### Task 3: Split a seam-only face by closed cuts

**Files:**
- Modify: `include/geom/brep/trim_face.hpp` (`split_face_curved`, around `:755-1750`)
- Test: `tests/unit/brep_trim_face_test.cpp`

**Interfaces:**
- Consumes: `detail::FaceUvProjection`, `face_uv_eval`, `Arrangement2`, `detail::split_segments`.
- Produces: unchanged `CurvedFragment` structure. New contract: for a face whose loop-derived uv area is zero, `split_face_curved` synthesizes the parameter rectangle and emits fragments whose loops contain **only** cut geometry; synthetic seam/pole chains never appear. Open cuts return `unsupported "split_face_curved: open cuts on a seam-only face"`; a cut that runs along a domain side (meridian or pole ring) returns `unsupported "split_face_curved: cut lies on the parameter domain boundary"`; a cycle set that yields no material fragment returns `numeric_failure "split_face_curved: no fragments"`. Closed cuts that do not meet the seam in exactly one identified point are rejected: `unsupported "split_face_curved: closed cut does not reach the seam"` / `"... closed cut crosses the seam twice"` (multi-chain joining is deferred; Task 5 reopens it if its SSI cuts cross the seam twice).

- [ ] **Step 1: Write the failing test**

```cpp
TEST_CASE("trim face: a sphere splits into caps and bands") {
  const Result<Brep> sphere = make_sphere(Point3d{0, 0, 0}, 1.0, 1e-9);
  REQUIRE(sphere.has_value());
  const FaceId face = sphere.value().face_ids().front();

  // z = 0.5, CCW seen from +z. Only the uv loop is consumed by the splitter.
  const double v = std::asin(0.5);
  std::vector<Point2d> uv;
  for (int i = 0; i <= 64; ++i) {
    const double t = 2.0 * 3.14159265358979323846 * i / 64.0;
    uv.emplace_back(t, v);  // includes u = 2*pi so both ends land on the seam copies
  }
  FaceCut cut;
  cut.uv = uv;
  const Result<std::vector<CurvedFragment>> fragments =
      split_face_curved(sphere.value(), face, std::span<const FaceCut>(&cut, 1), 1e-9);
  REQUIRE(fragments.has_value());
  REQUIRE(fragments.value().size() == 2u);
  // Both fragments are disks bounded by the same circle.
  CHECK(fragments.value()[0].loops.size() == 1u);
  CHECK(fragments.value()[1].loops.size() == 1u);
}
```

- [ ] **Step 2: Run it and confirm it fails**

Expected: the splitter emits degenerate sliver fragments (or an error), not 2 clean caps: the failing assertion is `fragments.value().size() == 2u`.

- [ ] **Step 3: Implement**

1. After the existing loop/bbox computation, if `std::fabs(signed_area2(loops.front())) == 0.0` and the surface is `Sphere3` or `Torus3`, mark the face seam-only. (Cone is Task 7.)
2. Get the full domain from `curved_boolean_surface_uv_domain`; synthesize four boundary chains in uv with a provenance tag `synthetic = true`: bottom `v = v_lo`, top `v = v_hi`, left `u = u_lo`, right `u = u_hi`. Sample each side at the same spacing used for cuts (at most `π/32`), and mark the four corners.
3. Build the `segments` vector exactly as today, but append the four synthetic chains (with a per-segment `synthetic` flag carried alongside).
4. Normalize every cut into the domain before use: walk the polyline and add whole `2π` multiples so consecutive `u` differ by at most `π` (the same unwrap the face projection uses), then shift the whole chain into `[u_lo, u_hi)`; shift `v` by whole `2π` only when the surface is periodic in `v` (torus). A cut with any `v` outside `[v_lo, v_hi]` by more than `1e-9 * max(1, |v|)` returns `unsupported "split_face_curved: cut leaves the parameter domain"`. If a cut has an endpoint that is **not** on the synthetic rectangle and is not closed, return `unsupported "split_face_curved: open cuts on a seam-only face"`.
5. Run `detail::split_segments` + `Arrangement2::build` as today. For every cycle array, rotate to a cut-run head (`synthetic[k] == false && synthetic[(k+n-1) % n] == true`) and walk exactly `n` edges modulo `n`; split into maximal non-synthetic runs, drop synthetic runs, and skip all-synthetic cycles. A run that does not start at a head is a bug, not a fragment.
6. Close a run only by one of: bitwise-equal 3D endpoints; the two u-seam copies (`u_lo`/`u_hi`) with `|dv| <= 1e-12 * max(1, |v|)` when periodic in `u`; the two v-seam copies with the analogous `du` bound when periodic in `v`. The seam-pair closure appends the first point (the two copies differ by rounding, not by geometry). Two different poles must never close each other. Otherwise return `numeric_failure "split_face_curved: synthetic boundary chain is open"`.
7. Reject unsupported closed cuts before any arrangement work: a closed cut with no seam contact -> `unsupported "... closed cut does not reach the seam"`; with two or more distinct seam contacts -> `unsupported "... closed cut crosses the seam twice"`; only a single identified seam point with matching `v` is accepted. Also reject side-coincident cuts: a cut with a maximal run of at least two consecutive points within `1e-12 * max(1, |coord|)` of one domain side (meridian ring, pole ring, or any segment along a side) returns `unsupported "split_face_curved: cut lies on the parameter domain boundary"`. Do not rely on geometric stripping to hide them: a real cut half that coincides with the synthetic side would be dropped.
8. Emit one `CurvedFragment` per material cycle with exactly those closed rings (seam-first outer ordering: a seam-spanning ring has uv area 0 because it is a chord; then descending |uv area|; the rest as holes). `arc_spans` stays empty: cuts are polylines, not analytic arcs. `reversed` stays false. If no material cycle yields a fragment, return `numeric_failure "split_face_curved: no fragments"` instead of an empty success.

- [ ] **Step 4: Verify**

Debug/Release: the new test yields exactly 2 fragments; add a second case with two latitude cuts (`v = ±asin(0.8)`) expecting 3 fragments with loop counts `{1, 2, 1}`. Compute each fragment's area by lifting its loop to 3D (`spherical_area`), not by uv shoelace: a latitude ring is a chord with uv area 0. One cut at `asin(0.5)` gives `{π, 3π}`; two cuts at `±asin(0.8)` give `{0.4π, 3.2π, 0.4π}`; both total `4π` within `1e-12`. Assert every loop closes in 3D within `2e-12` across the seam copy, and add the meridian/pole/open/interior/tangent/drift rejection-or-acceptance cases.

- [ ] **Step 5: Commit**

`git commit -m "feat(m3c2): split seam-only faces by closed cuts"`

---

### Task 4: Pole-aware volume of sphere faces

**Files:**
- Modify: `include/geom/brep/detail/curved_volume.hpp:161-380`
- Test: `tests/unit/brep_curved_volume_test.cpp`

**Interfaces:**
- Consumes: `curved_solid_volume(const Brep&, SolidId, double)`.
- Produces: same signature. New contract: faces on `Sphere3` with poles on their boundary (or a degenerate seam loop) integrate exactly; `curved_solid_volume` of a unit sphere is `4π/3` within `1e-12` relative.

- [ ] **Step 1: Write the failing test**

```cpp
TEST_CASE("curved volume: a sphere integrates through its poles") {
  const Result<Brep> sphere = make_sphere(Point3d{0, 0, 0}, 1.0, 1e-9);
  REQUIRE(sphere.has_value());
  const SolidId solid = sphere.value().solid_ids().front();
  const Result<double> volume = curved_solid_volume(sphere.value(), solid, 1e-9);
  REQUIRE(volume.has_value());
  CHECK(volume.value() == doctest::Approx(4.0 * 3.14159265358979323846 / 3.0).epsilon(1e-12));
}
```

- [ ] **Step 2: Run it and confirm it fails**

Expected: `Error{Errc::degenerate, "curved_solid_volume: face triangulation failed"}`.

- [ ] **Step 3: Implement**

In the face triangulation helper, when the face's loop-derived uv area is zero and the surface is `Sphere3`:
1. Build the full parameter rectangle from `curved_boolean_surface_uv_domain` (or an equivalent local helper — `detail` must not depend on `boolean_solid_curved.hpp`; duplicate the four constants with a comment).
2. Triangulate the rectangle into `(nu × nv)` quads with `nu = 2 * nv` and `nv` chosen so the sagitta is below `1e-4` (e.g. `nv = 64` for the unit sphere); collapse every `v = -π/2` vertex to the south pole and every `v = +π/2` vertex to the north pole, and weld `u = 0` to `u = 2π`.
3. For a trimmed sphere face (non-zero loop area), triangulate the uv loops as today, but when a loop contains a vertex at `|v| == π/2`, replace that vertex with the single pole point and close the fan there.

Do **not** approximate the surface: the existing quadrature evaluates `surface_eval` at the barycentric points, so the sphere remains exact. Only the triangle topology changes.

- [ ] **Step 4: Verify**

The unit sphere volume test passes at `1e-12`; the existing cylinder/torus/trimmed-face volume tests are unchanged; Debug and Release green.

- [ ] **Step 5: Commit**

`git commit -m "feat(m3c2): integrate sphere volumes through their poles"`

---

### Task 5: Keep sphere operands explicitly unsupported

**Why deferred:** a sphere face is one seam slit with zero uv area. Probe evidence (2026-11-14) shows four independent gaps beyond this milestone: clip drops every cut against the zero-area loop, a cap ring needs a pole detour, `split_solid_heal_fragment` then folds the detour meridian into the cut circle's arc span (`assemble_curved_solid: arc run is inconsistent`), and `sphere ∩ sphere` needs cross-chain joining of two seam contacts. Faking support would produce invalid solids, so this task only pins the diagnostic.

**Files:**
- Modify: `include/geom/brep/boolean_solid_curved.hpp` (gate message for `Sphere3`)
- Modify: `tests/unit/brep_boolean_curved_test.cpp`
- Modify: `docs/architecture/geometry-library-architecture.md` §16.3 (ADR-0028 is created in Task 8)

**Interfaces:**
- Consumes: `curved_boolean_surface_gate` (`:157`), `boolean_solid_curved`
- Produces: `Sphere3` gate error `"sphere operands need cap topology; booleans are deferred to M3c-3"`; no other behavior change. Torus/cone gates are opened by Tasks 6 and 7.

**Step 1: Write the failing test**

```cpp
TEST_CASE("curved boolean: sphere operands are rejected with the cap diagnostic") {
  const Result<Brep> sphere = make_sphere(Point3d{0, 0, 0}, 1.0, 1e-9);
  const Result<Brep> box = make_box(Point3d{-2, -2, 0.5}, Point3d{2, 2, 2}, 1e-9);
  REQUIRE(sphere.has_value());
  REQUIRE(box.has_value());
  const SolidId a = sphere.value().solid_ids().front();
  const SolidId b = box.value().solid_ids().front();
  CurvedBooleanReport report;
  for (SolidBooleanOp op : {SolidBooleanOp::union_, SolidBooleanOp::intersection,
                            SolidBooleanOp::difference}) {
    const Result<Brep> result =
        boolean_solid_curved(sphere.value(), a, box.value(), b, op, 1e-9, &report);
    REQUIRE_FALSE(result.has_value());
    CHECK(result.error().code == Errc::unsupported);
    CHECK(std::string_view(result.error().message).find("deferred to M3c-3") !=
          std::string_view::npos);
    CHECK(report.face_pairs == 0);
  }
}
```

**Step 2:** run it; it fails because the gate currently says `"sphere operands are deferred to M3c"`.

**Step 3:** replace the sphere gate message with the cap diagnostic. Keep the `Sphere3` pointer as null; do not touch torus/cone branches.

**Step 4:** Debug and Release suites green. Expected 918 cases; the new case adds 9 assertions.

**Step 5:** `git commit -m "test(m3c2): pin the sphere cap diagnostic"`

### Task 6: Enable torus operands and the v seam

**Probe findings (2026-11-14, R=2/r=0.5):** enabling the gate alone fails in three independent, small places: (1) `face_uv_project_torus` returns v in (-π,π] while the chart is [0,2π], and `curved_boolean_unwrap_axis` only unwraps deltas, so a cut at v=π emerges at -3.142 and is dropped; (2) the z=0 cut lies exactly on the v=0 seam edge of the parameter rectangle and produces 76 sliver fragments; (3) `assemble_curved_edge_key` collapses the two closed seam circles that share the point (2.5,0,0), giving one slot `uses=4`. Total ~50 LOC.

**Files:**
- Modify: `include/geom/brep/boolean_solid_curved.hpp` (gate ~:187, `curved_boolean_split_at_seam` :679, chart wrap helpers :697-712/:845/:1205, clip :1091)
- Modify: `include/geom/brep/assemble_curved.hpp` (`assemble_curved_edge_key` near :420)
- Test: `tests/unit/brep_boolean_curved_test.cpp`

**Interfaces:**
- Produces: `curved_boolean_split_at_seam(..., int axis)` where `axis == 0` splits on u and `axis == 1` on v; call sites pass `detail::FaceUvProjection::periodic_v() ? 1 : 0`. After unwrapping, chains are shifted by whole periods into `[0,2π]` on each periodic axis. Chains that lie entirely on a synthetic seam-only rectangle side are dropped before splitting. Torus/plane pairs outside the axis-perpendicular/axis-containing restriction return the documented `unsupported` message before pairing.

- [ ] **Step 1: Write the failing tests**

```cpp
TEST_CASE("curved boolean: torus and half space") {
  const Result<Brep> torus = make_torus(Point3d{0, 0, 0}, Vec3d{0, 0, 1}, 1.0, 0.5, 1e-9);
  const Result<Brep> box = make_box(Point3d{-2, -2, 0}, Point3d{2, 2, 2}, 1e-9);
  REQUIRE(torus.has_value());
  REQUIRE(box.has_value());
  const SolidId a = torus.value().solid_ids().front();
  const SolidId b = box.value().solid_ids().front();
  const Result<CurvedBooleanResult> half =
      boolean_solid_curved(torus.value(), a, box.value(), b,
                           SolidBooleanOp::intersection, 1e-9);
  REQUIRE(half.has_value());
  CHECK(half.value().report.max_cut_residual <= 1e-3);
  CHECK(validate_topology(half.value().brep).ok());
  CHECK(curved_solid_volume(half.value().brep, half.value().brep.solid_ids().front(), 1e-9).value() ==
        doctest::Approx(2.4674011002723395).epsilon(1e-6));
}

TEST_CASE("curved boolean: torus cut across the v seam") {
  const Result<Brep> torus = make_torus(Point3d{0, 0, 0}, Vec3d{0, 0, 1}, 1.0, 0.5, 1e-9);
  const Result<Brep> box = make_box(Point3d{0, -2, -2}, Point3d{2, 2, 2}, 1e-9);
  REQUIRE(torus.has_value());
  REQUIRE(box.has_value());
  const SolidId a = torus.value().solid_ids().front();
  const SolidId b = box.value().solid_ids().front();
  const Result<CurvedBooleanResult> half =
      boolean_solid_curved(torus.value(), a, box.value(), b,
                           SolidBooleanOp::intersection, 1e-9);
  REQUIRE(half.has_value());
  CHECK(validate_topology(half.value().brep).ok());
  CHECK(curved_solid_volume(half.value().brep, half.value().brep.solid_ids().front(), 1e-9).value() ==
        doctest::Approx(2.4674011002723395).epsilon(1e-6));
  CHECK(half.value().report.max_cut_residual <= 1e-3);
}

TEST_CASE("curved boolean: torus minus half space") {
  const Result<Brep> torus = make_torus(Point3d{0, 0, 0}, Vec3d{0, 0, 1}, 1.0, 0.5, 1e-9);
  const Result<Brep> box = make_box(Point3d{-2, -2, 0}, Point3d{2, 2, 2}, 1e-9);
  REQUIRE(torus.has_value());
  REQUIRE(box.has_value());
  const Result<CurvedBooleanResult> half =
      boolean_solid_curved(torus.value(), torus.value().solid_ids().front(), box.value(),
                           box.value().solid_ids().front(), SolidBooleanOp::difference, 1e-9);
  REQUIRE(half.has_value());
  CHECK(validate_topology(half.value().brep).ok());
}

TEST_CASE("curved boolean: a box that contains a torus keeps it whole") {
  const Result<Brep> torus = make_torus(Point3d{0, 0, 0}, Vec3d{0, 0, 1}, 1.0, 0.5, 1e-9);
  const Result<Brep> box = make_box(Point3d{-2, -2, -2}, Point3d{2, 2, 2}, 1e-9);
  REQUIRE(torus.has_value());
  REQUIRE(box.has_value());
  const Result<CurvedBooleanResult> whole =
      boolean_solid_curved(torus.value(), torus.value().solid_ids().front(), box.value(),
                           box.value().solid_ids().front(), SolidBooleanOp::intersection, 1e-9);
  REQUIRE(whole.has_value());
  CHECK(validate_topology(whole.value().brep).ok());
  CHECK(curved_solid_volume(whole.value().brep, whole.value().brep.solid_ids().front(), 1e-9).value() ==
        doctest::Approx(4.934802200544679).epsilon(1e-6));
}
```

- [ ] **Step 2: Run and confirm they fail**

Expected: gate error `"torus operands are deferred to M3c"` for every case.

- [ ] **Step 3: Implement**

1. Gate: `Torus3` → supported.
2. Chart wrap: add `curved_boolean_wrap_axis(chain, axis, period)` that shifts every point by whole periods into `[0,2π]` after unwrapping; call it for u always and for v when `periodic_v()` (both the cut chain and the foreign chain before clipping).
3. Seam-side chains: when a chain lies on a synthetic rectangle side (all points on that edge within `1e-12`), drop it before splitting; only cuts that cross the rectangle interior are real.
4. Thread `axis` through `curved_boolean_split_at_seam`; `curved_boolean_clip_pair` wraps `foreign_uv[axis]`.
5. `assemble_curved_edge_key`: add circle centre/normal/radius to the key so the two closed seam circles at (2.5,0,0) stay distinct.

- [ ] **Step 4: Verify**

All four torus cases pass with `max_cut_residual <= 1e-3`; `curved_solid_volume` of the untouched torus stays exact; all earlier suites unchanged.

- [ ] **Step 5: Commit**

`git commit -m "feat(m3c2): support torus operands and the v seam"`

---

### Task 7: Cone diagnostic (booleans deferred to M3c-3)

**Why deferred:** the seam-only splitter accepts only a closed cut that is anchored at exactly one seam point. A cone/plane cut is either (a) an interior closed circle when the plane is perpendicular to the axis below the apex — rejected because its front and back points coincide (`"closed cut crosses the seam twice"`), or (b) one or two open meridian chains ending on the base rim (plane containing the axis, possibly through the apex) — rejected as open cuts. Delivering cone booleans therefore needs the general cut splitter (closed rings not touching a seam, open chains ending on a real boundary) plus apex welding and an apex volume fan; `curved_solid_volume(cone)` also fails today at the apex projection. All of that moves to M3c-3.

**Files:**
- Modify: `include/geom/brep/boolean_solid_curved.hpp` (cone gate message only)
- Modify: `tests/unit/brep_boolean_curved_test.cpp`
- Modify: `docs/architecture/geometry-library-architecture.md` §16.3 (one paragraph)

**Interfaces:**
- Produces: cone gate returns `unsupported` with `"boolean_solid_curved: cone operands need closed-cut and apex topology; booleans are deferred to M3c-3"`; `report.face_pairs == 0`, `fragments == 0`, `kept == 0`. No projection, splitter or volume code changes.

- [ ] **Step 1: Write the failing test**

```cpp
TEST_CASE("curved boolean: cone operands are deferred with a diagnostic") {
  const Result<Brep> cone = make_cone(Point3d{0, 0, 0}, Vec3d{0, 0, 1}, 1.0, 2.0, 1e-9);
  const Result<Brep> box = make_box(Point3d{-2, -2, -1}, Point3d{2, 2, 2}, 1e-9);
  REQUIRE(cone.has_value());
  REQUIRE(box.has_value());
  for (const SolidBooleanOp op : {SolidBooleanOp::union_, SolidBooleanOp::intersection,
                                  SolidBooleanOp::difference}) {
    CurvedBooleanReport report{};
    report.face_pairs = 999;
    const Result<CurvedBooleanResult> out =
        boolean_solid_curved(cone.value(), cone.value().solid_ids().front(), box.value(),
                             box.value().solid_ids().front(), op, 1e-9, &report);
    REQUIRE_FALSE(out);
    CHECK(out.error().code == Errc::unsupported);
    const std::string message(out.error().message);
    CHECK(message.find("deferred to M3c-3") != std::string::npos);
    CHECK(report.face_pairs == 0);
    CHECK(report.fragments == 0);
    CHECK(report.kept == 0);
  }
}
```

- [ ] **Step 2: Run and confirm it fails**

Expected: the current gate message `"boolean_solid_curved: apex-touching loops are deferred"` lacks the M3c-3 wording and is reported after the box/cone faces are counted, so the message and `face_pairs == 0` assertions fail.

- [ ] **Step 3: Implement**

Replace the cone branch message with the exact diagnostic above; keep the cone pointer null and do not touch any other gate or code path.

- [ ] **Step 4: Verify**

Debug and Release full suites green (expect +1 case / +3 assertions per operation loop), zero warnings; the three torus `REQUIRE_MESSAGE` diagnostics and all existing tests unchanged.

- [ ] **Step 5: Commit**

`git commit -m "test(m3c2): pin the deferred cone diagnostic"`

---

### Task 8: QA-14, benchmark, ADR and status

**Files:**
- Create: `bench/bench_m3c2.cpp` (or extend `bench/bench_m3b2.cpp` if that target is already registered)
- Create: `docs/adr/0028-curved-boolean-surfaces.md`
- Create: `docs/benchmarks/m3c2.md`
- Modify: `docs/architecture/geometry-library-architecture.md` §16.3
- Modify: `tests/CMakeLists.txt` if a new test file is added

**QA-14 corpus** (each case: success, `validate_topology` clean, closed shell, volume within tolerance, two runs byte-identical, `report.max_cut_residual` within the stated bound):

| Case | Volume | Tolerance | Residual |
| --- | --- | --- | --- |
| torus ∩ half space `z ≥ 0` (`R=1`, `r=0.5`) | 2.4674011002723395 | 1e-6 rel | ≤ 1e-3 |
| torus ∩ half space `x ≥ 0` (crosses the v seam) | 2.4674011002723395 | 1e-6 rel | ≤ 1e-3 |
| torus inside a big box (no cut) | 4.934802200544679 | 1e-6 rel | n/a |

Cone and sphere rows are intentionally absent: both are pinned as `unsupported` diagnostics (Tasks 5 and 7), not measured.

- [ ] **Step 1: Add the QA test**

Create `tests/unit/curved_boolean_qa14_test.cpp` with one `TEST_CASE("QA-14: curved booleans of analytic primitives")` looping the table; require `report.cuts >= 1` for every cut row and `report.kept >= 1`; the no-cut row (torus inside a box) asserts `cuts == 0` and skips the residual bound. Run each operation twice and compare `curved_solid_volume` bitwise.

- [ ] **Step 2: Confirm it fails before the fixes are merged** (it should, on the torus gate)

- [ ] **Step 3: Add the benchmark**

Time 100 iterations of each case on a Release build; report median and p95 in `docs/benchmarks/m3c2.md`. This is informational only — no gate.

- [ ] **Step 4: Write ADR-0028**

Record: synthetic parameter boundary and why it is stripped; the seam-anchored closed-cuts-only limit; torus chordal cuts; why cone needs interior-closed-ring and open-chain splitting plus apex welding/fan; what remains for M3c-3 (general trimmed tessellation, those cut kinds, NURBS operands, exact torus/cone healing).

- [ ] **Step 5: Update §16.3**

Move torus operands from "deferred" to "supported with limits"; record sphere (cap topology) and cone (closed-cut/apex topology) as explicitly unsupported with their Task 5/7 diagnostics. Link ADR-0028 and the benchmark.

- [ ] **Step 6: Verify and commit**

`git commit -m "test(m3c2): add QA-14 curved boolean matrix"`

---

## M3c-2 exit criteria

- [ ] All tasks 1-8 complete; every commit reviewed and approved.
- [ ] Debug and Release: full suite green, zero warnings; 925+ cases.
- [ ] QA-14 torus table passes with the stated tolerances on MSVC; CI on ubuntu/macos passes with the documented relative tolerances (analytic values must not be bit-pinned across libms).
- [ ] Existing planar/cylinder curved boolean tests, `split_solid` tests and trimmed tessellation tests unchanged.
- [ ] ADR-0028 and §16.3 describe exactly the supported surface/operation matrix and the remaining limits.

## Self-check

- [ ] Every task names exact files, interfaces and a runnable failing test.
- [ ] No task depends on code from a later task (Task 5 depends on 1-4; Task 6 depends on 1 and 3; Task 7 depends on 1).
- [ ] No placeholders or TBDs.
- [ ] QA values are independent analytic formulas, not values read from the implementation.
- [ ] The plan states the deviations from the architecture (seam-anchored closed cuts, synthetic boundary, chordal torus, sphere/cone-boolean and general-tessellation deferral to M3c-3) explicitly.
