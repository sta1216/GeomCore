# M3c-2 Curved Boolean Completion Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Make `boolean_solid_curved` accept sphere, torus and cone operands, so the curved B-rep kernel can intersect, unite and subtract the analytic primitives it can already tessellate and measure.

**Architecture:** The driver already owns pairing, SSI cut extraction, classification, fragment healing and assembly. Three gaps block the new surface types: (1) a sphere/cone face whose only loop is the seam slit has a degenerate uv domain, so SSI never sees a patch; (2) `split_face_curved` cannot cut a face whose outer boundary has zero uv area; (3) `curved_solid_volume` cannot triangulate a face with poles. Torus needs only the gate plus v-seam handling; cone additionally needs an apex-tolerant projection. Each gap is closed in its own task, then the surface gates are opened and pinned by end-to-end volume tests.

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

M3c-2 covers **sphere, torus and cone operands of `boolean_solid_curved`**, plus the pole-aware volume needed to verify them. It is deliberately split from general trimmed-face tessellation:

1. **Closed cuts only on seam-only faces.** A sphere/cone face is split only when every cut is a closed ring (or a ring touching the face's own seam). A cut that ends in the middle of a face still returns `unsupported "split_face_curved: open cuts on a seam-only face"`. This covers all sphere/plane, sphere/sphere, sphere/cylinder, cone/plane and cone/cylinder cases in the QA corpus.
2. **Synthetic parameter boundary.** For a seam-only face the splitter synthesizes the missing rectangle of the parameter domain (`u ∈ [u0,u0+2π]`, `v ∈ [v_lo,v_hi]`), runs the existing arrangement, then strips every chain that came from the synthetic seam/pole edges and closes the remaining cut chain by 3D endpoint identity. The synthetic edges never appear in `CurvedFragment::loops`, so assembly sees only real cut edges.
3. **Torus cuts stay chordal.** There is no exact torus/plane or torus/cylinder healing (the intersection is quartic in general). Quality is reported by `CurvedBooleanReport::max_cut_residual`; the tests require it below `1e-3` for the chosen corpus.
4. **Cone support is plane/cylinder cuts**, including cuts through the apex; cuts tangent to the apex or to the base rim remain `unsupported`.
5. **No trimmed tessellation of boolean results in this milestone.** QA-14 checks topology, closedness, volume and determinism; tessellating a result face with a non-rectilinear boundary is M3c-3 work (M3c-1 explicitly rejects those with `"... boundary is not rectilinear"`).

## Review Focus (feed these to every reviewer)

- Does the sphere/cone uv-domain synthesis ever *widen* a legitimate trimmed domain? A face whose loop-derived area is non-zero must keep exactly its old range.
- Are the synthetic seam/pole chains really absent from every emitted `CurvedFragment` loop, and do the cut chains close with bitwise-identical 3D endpoints?
- Does the pole-aware volume integrate the same surface the tessellator would (no flat-triangle approximation of the sphere)?
- Do torus cuts that cross `v = 0 / 2π` produce loops with no `2π` jump?
- Is the apex a single vertex in the assembled Brep (no duplicate apex vertices, no zero-length edges)?
- Are all error paths preserving the failure report and returning a message that names the stage?

## File Map

| File | Responsibility |
| --- | --- |
| `include/geom/brep/boolean_solid_curved.hpp` | failure report, uv-domain synthesis, gate widening, seam split axis, sphere-sphere exact circle |
| `include/geom/brep/trim_face.hpp` | seam-only boundary synthesis and synthetic-chain stripping; apex collapse |
| `include/geom/brep/face_uv.hpp` | apex-tolerant cone projection |
| `include/geom/brep/detail/curved_volume.hpp` | pole-aware face triangulation |
| `include/geom/brep/split_solid.hpp` | sphere/sphere exact circle helper (plane construction) |
| `tests/unit/brep_boolean_curved_test.cpp` | end-to-end sphere/torus/cone booleans |
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
  const SolidId solid = sphere.value().solid_ids().front();
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
Expected: `u[1]` is 0 (not 2π) and the bounds are the xz slice `x=y=0`, so the test fails.

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
- Produces: unchanged `CurvedFragment` structure. New contract: for a face whose loop-derived uv area is zero, `split_face_curved` synthesizes the parameter rectangle and emits fragments whose loops contain **only** cut geometry; synthetic seam/pole chains never appear. Open cuts on such a face return `Errc::unsupported` with message `"split_face_curved: open cuts on a seam-only face"`.

- [ ] **Step 1: Write the failing test**

```cpp
TEST_CASE("trim face: a sphere splits into caps and bands") {
  const Result<Brep> sphere = make_sphere(Point3d{0, 0, 0}, 1.0, 1e-9);
  REQUIRE(sphere.has_value());
  const SolidId solid = sphere.value().solid_ids().front();
  const FaceId face = sphere.value().solid_faces(solid).front();

  // z = 0.5, CCW seen from +z. Only the uv loop is consumed by the splitter.
  const double v = std::asin(0.5);
  std::vector<Point2d> uv;
  for (int i = 0; i < 64; ++i) {
    const double t = 2.0 * 3.14159265358979323846 * i / 64.0;
    uv.emplace_back(t, v);
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

Expected: `fragments.value().size() == 0`.

- [ ] **Step 3: Implement**

1. After the existing loop/bbox computation, if `std::fabs(signed_area2(loops.front())) == 0.0` and the surface is `Sphere3` or `Torus3`, mark the face seam-only. (Cone is Task 7.)
2. Get the full domain from `curved_boolean_surface_uv_domain`; synthesize four boundary chains in uv with a provenance tag `synthetic = true`: bottom `v = v_lo`, top `v = v_hi`, left `u = u_lo`, right `u = u_hi`. Sample each side at the same spacing used for cuts (at most `π/32`), and mark the four corners.
3. Build the `segments` vector exactly as today, but append the four synthetic chains (with a per-segment `synthetic` flag carried alongside).
4. Normalize every cut into the domain before use: walk the polyline and add whole `2π` multiples so consecutive `u` differ by at most `π` (the same unwrap the face projection uses), then shift the whole chain into `[u_lo, u_hi)`; shift `v` by whole `2π` only when the surface is periodic in `v` (torus). A cut with any `v` outside `[v_lo, v_hi]` by more than `1e-9 * max(1, |v|)` returns `unsupported "split_face_curved: cut leaves the parameter domain"`. If a cut has an endpoint that is **not** on the synthetic rectangle and is not closed, return `unsupported "split_face_curved: open cuts on a seam-only face"`.
5. Run `detail::split_segments` + `Arrangement2::build` as today. When walking a cycle ring, split it into maximal runs of non-synthetic vertices; drop synthetic runs. Each remaining run is a chain of cut points:
   - if its two 3D endpoints are bitwise equal, close it into a ring;
   - else, if both endpoints project to the same synthetic seam side or to the same pole value (`v == v_lo` / `v == v_hi`), close it by appending the first point;
   - otherwise return `numeric_failure "split_face_curved: synthetic boundary chain is open"`.
6. Emit one `CurvedFragment` per material cycle with exactly those closed rings (outer first by descending |uv area|, the rest as holes). `arc_spans` stays empty: cuts are polylines, not analytic arcs. `reversed` stays false.

- [ ] **Step 4: Verify**

Debug/Release: the new test yields exactly 2 fragments; add a second case with two latitude cuts (`v = ±asin(0.8)`) expecting 3 fragments with loop counts `{1, 2, 1}` and total sphere area `4π` within `1e-12` when each loop is lifted and its spherical area integrated numerically (or assert the uv areas of the caps `0.4π`, band `2π`, cap `0.4π`).

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

### Task 5: Enable sphere operands end to end

**Files:**
- Modify: `include/geom/brep/boolean_solid_curved.hpp:152-200` (gate), `:1500-1800` (pair circle healing)
- Modify: `include/geom/brep/split_solid.hpp` (expose the plane/sphere circle helper if not already public in `detail`)
- Test: `tests/unit/brep_face_uv_test.cpp`

**Interfaces:**
- Consumes: Task 1-4.
- Produces: `curved_boolean_surface_gate` returns supported for `Sphere3`; `curved_boolean_pair_circle` gains a sphere/sphere branch:

```cpp
// center line n = (c_b - c_a)/d, d = |c_b - c_a|;
// plane n . x = (|c_b|^2 - |c_a|^2 + r_a^2 - r_b^2) / (2d)  [relative to origin]
```

then delegate to the existing sphere/plane circle construction. The assembly remap at `boolean_solid_curved.hpp:1919` is reworded to `"boolean_solid_curved: assembled result is an open shell"` (keep the underlying `assemble_curved_solid` error as `message` context by appending it).

- [ ] **Step 1: Write the failing tests**

```cpp
TEST_CASE("curved boolean: sphere and box") {
  const Result<Brep> sphere = make_sphere(Point3d{0, 0, 0}, 1.0, 1e-9);
  const Result<Brep> box = make_box(Point3d{-2, -2, 0.5}, Point3d{2, 2, 2.5},
                                    1e-9);
  REQUIRE(sphere.has_value());
  REQUIRE(box.has_value());
  const SolidId solid_a = sphere.value().solid_ids().front();
  const SolidId solid_b = box.value().solid_ids().front();
  const Result<CurvedBooleanResult> cut =
      boolean_solid_curved(sphere.value(), solid_a, box.value(), solid_b,
                           SolidBooleanOp::intersection);
  REQUIRE(cut.has_value());
  CHECK(cut.value().report.cuts >= 1u);
  CHECK(curved_solid_volume(cut.value().brep, cut.value().brep.solid_ids().front(), 1e-9).value() ==
        doctest::Approx(0.6544984694978736).epsilon(1e-9));
  CHECK(validate_topology(cut.value().brep).ok());
}
```

Include `<geom/brep/detail/curved_volume.hpp>` and `<geom/brep/validate.hpp>` at the top of the test file. Add sibling cases: `sphere - box` → `3.534291735288517`; `sphere ∩ sphere` at distance 1.5 → `0.3599741582`; `sphere ∩ cylinder r=0.6` → `3.954295`; and a regression `cylinder ∩ big box` unchanged (`4.523893`).

- [ ] **Step 2: Run and confirm they fail**

Expected: the current gate error `"sphere operands are deferred to M3c"` for the first three; the cylinder case may fail with `"touches itself"` until the caps are cut.

- [ ] **Step 3: Implement**

1. Gate: `Sphere3` → `{true, nullptr}`.
2. Exact circle: add the sphere/sphere branch described above.
3. Reword the open-shell remap; do not change `assemble_curved_solid`.
4. If `sphere ∩ cylinder` still leaves an open tube, inspect `kept`/`dropped` with the failure report: every cylinder cap must be split by the sphere's two cut circles. If the caps are dropped because their representative is misclassified, fix the representative selection (largest uv triangle whose centroid is not on the cut) rather than loosening assembly.

- [ ] **Step 4: Verify**

All five new cases pass on Debug and Release; existing curved/planar boolean suites unchanged; volume tolerance `1e-9` relative.

- [ ] **Step 5: Commit**

`git commit -m "feat(m3c2): support sphere operands in curved booleans"`

---

### Task 6: Enable torus operands and the v seam

**Files:**
- Modify: `include/geom/brep/boolean_solid_curved.hpp:187-190` (gate), `:679` (`curved_boolean_split_at_seam`), `:845`, `:1091`, `:1205`
- Test: `tests/unit/brep_face_uv_test.cpp`

**Interfaces:**
- Produces: `curved_boolean_split_at_seam(..., int axis)` where `axis == 0` splits on u and `axis == 1` on v; call sites pass `detail::FaceUvProjection::periodic_v() ? 1 : 0` (u for every periodic surface, v additionally for a torus). `curved_boolean_clip_pair` wraps the foreign uv coordinate selected by the same axis.

- [ ] **Step 1: Write the failing tests**

```cpp
TEST_CASE("curved boolean: torus and half space") {
  const Result<Brep> torus = make_torus(Point3d{0, 0, 0}, Vec3d{0, 0, 1}, 1.0, 0.5, 1e-9);
  const Result<Brep> box = make_box(Point3d{-2, -2, 0}, Point3d{2, 2, 2},
                                    1e-9);
  REQUIRE(torus.has_value());
  REQUIRE(box.has_value());
  const SolidId solid_a = torus.value().solid_ids().front();
  const SolidId solid_b = box.value().solid_ids().front();
  const Result<CurvedBooleanResult> half =
      boolean_solid_curved(torus.value(), solid_a, box.value(), solid_b,
                           SolidBooleanOp::intersection);
  REQUIRE(half.has_value());
  CHECK(half.value().report.max_cut_residual <= 1e-3);
  CHECK(curved_solid_volume(half.value().brep, half.value().brep.solid_ids().front(), 1e-9).value() ==
        doctest::Approx(2.4674011002723395).epsilon(1e-6));
  CHECK(validate_topology(half.value().brep).ok());
}
```

Add a tilted-plane case (`z = 0.3·x`, a plane crossing the v seam) that asserts every fragment loop has `|Δv| < π` between consecutive points.

- [ ] **Step 2: Run and confirm they fail**

Expected: gate error `"torus operands are deferred to M3c"`.

- [ ] **Step 3: Implement**

1. Gate: `Torus3` → supported.
2. Thread `axis` through `curved_boolean_split_at_seam`; at both call sites use the face's projection (`periodic_v()`).
3. In `curved_boolean_clip_pair`, wrap `foreign_uv[axis]` by `2π` when the projection is periodic in that axis.
4. In the cut unwrap (`:845` and `:1205`), unwrap whichever coordinate is periodic — u always, v when `periodic_v()`.

- [ ] **Step 4: Verify**

Both torus tests pass; the tilted cut has `max_cut_residual <= 1e-3`; all earlier suites unchanged.

- [ ] **Step 5: Commit**

`git commit -m "feat(m3c2): support torus operands and the v seam"`

---

### Task 7: Cone projection, apex splitting and gate

**Files:**
- Modify: `include/geom/brep/face_uv.hpp:132-160`
- Modify: `include/geom/brep/trim_face.hpp` (apex collapse in the seam-only splitter)
- Modify: `include/geom/brep/boolean_solid_curved.hpp:192-196` (gate)
- Test: `tests/unit/brep_boolean_curved_test.cpp`, `tests/unit/brep_face_uv_test.cpp`

**Interfaces:**
- Produces: `face_uv_project_cone` accepts a point within `tolerance` of the apex, returning `{previous_u (or 0), v_max}`; points beyond the apex remain `non_finite`/`degenerate` as documented. Each cone face's uv rectangle has a top edge `v = v_max` that collapses to one apex point; the splitter emits at most one apex vertex per fragment and never a zero-length edge.

- [ ] **Step 1: Write the failing tests**

```cpp
TEST_CASE("curved boolean: cone and box") {
  const Result<Brep> cone = make_cone(Point3d{0, 0, 0}, Vec3d{0, 0, 1}, 0.5, 1.0, 1e-9);
  const Result<Brep> box = make_box(Point3d{-2, -2, -0.5}, Point3d{2, 2, 1.5},
                                    1e-9);
  REQUIRE(cone.has_value());
  REQUIRE(box.has_value());
  const SolidId solid_a = cone.value().solid_ids().front();
  const SolidId solid_b = box.value().solid_ids().front();
  const Result<CurvedBooleanResult> kept =
      boolean_solid_curved(cone.value(), solid_a, box.value(), solid_b,
                           SolidBooleanOp::intersection);
  REQUIRE(kept.has_value());
  CHECK(curved_solid_volume(kept.value().brep, kept.value().brep.solid_ids().front(), 1e-9).value() ==
        doctest::Approx(0.2617993877991494).epsilon(1e-9));
  CHECK(validate_topology(kept.value().brep).ok());
}
```

Add a half-cone cut (`x = 0` plane through the apex) asserting volume `0.1308997` and exactly one apex vertex in the result (scan vertices, count those equal within `1e-12`).

- [ ] **Step 2: Run and confirm they fail**

Expected: `"face_loops_uv: cone projection is undefined at or beyond the apex"`.

- [ ] **Step 3: Implement**

1. `face_uv_project_cone`: if `v >= v_max` and the point is within `tolerance` of the apex, return `{previous_u, v_max}`; if `v > v_max` and far from the apex, keep the current error.
2. Seam-only synthesis (Task 3) includes `Cone3`; synthesize the rectangle `[0,2π] × [0,v_max]`. When stripping synthetic chains, the whole top edge is synthetic; collapse every remaining point with `v == v_max` to one 3D apex point before closing chains.
3. Cone gate → supported, message removed.
4. If the assembled result has duplicate apex vertices, weld them in `split_solid_heal_fragment`/`assemble_curved_solid` by 3D identity, not by parameter value.

- [ ] **Step 4: Verify**

Both cone cases pass; the full cone volume matches `πr²h/3` at `1e-9`; the apex appears once; existing cylinder/sphere/torus cases unchanged.

- [ ] **Step 5: Commit**

`git commit -m "feat(m3c2): support cone operands with a single apex vertex"`

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
| sphere ∩ box (`z ≥ 0.5`) | 0.6544984694978736 | 1e-9 rel | ≤ 1e-9 |
| sphere − box (same box) | 3.534291735288517 | 1e-9 rel | ≤ 1e-9 |
| sphere ∩ sphere (`d = 1.5`) | 0.3599741582 | 1e-9 rel | ≤ 1e-9 |
| sphere ∩ cylinder (`r = 0.6`, `z ∈ [-2,2]`) | 3.954295 | 1e-9 rel | ≤ 1e-9 |
| torus ∩ half space (`R=1`, `r=0.5`) | 2.4674011002723395 | 1e-6 rel | ≤ 1e-3 |
| cone ∩ box (whole cone) | 0.2617993877991494 | 1e-9 rel | ≤ 1e-9 |

- [ ] **Step 1: Add the QA test**

Create `tests/unit/curved_boolean_qa14_test.cpp` with one `TEST_CASE("QA-14: curved booleans of analytic primitives")` looping the table; require `report.cuts >= 1` for every cut case and `report.kept >= 1`; run each operation twice and compare `curved_solid_volume` bitwise.

- [ ] **Step 2: Confirm it fails before the fixes are merged** (it should, on the sphere gate)

- [ ] **Step 3: Add the benchmark**

Time 100 iterations of each case on a Release build; report median and p95 in `docs/benchmarks/m3c2.md`. This is informational only — no gate.

- [ ] **Step 4: Write ADR-0028**

Record: synthetic parameter boundary and why it is stripped; closed-cuts-only limit; torus chordal cuts; cone apex handling; what remains for M3c-3 (general trimmed tessellation, open cuts, NURBS operands, exact torus/cone healing).

- [ ] **Step 5: Update §16.3**

Move sphere/torus/cone operands from "deferred" to "supported with limits", link ADR-0028 and the benchmark.

- [ ] **Step 6: Verify and commit**

`git commit -m "test(m3c2): add QA-14 curved boolean matrix"`

---

## M3c-2 exit criteria

- [ ] All tasks 1-8 complete; every commit reviewed and approved.
- [ ] Debug and Release: full suite green, zero warnings; 900+ cases.
- [ ] QA-14 table passes with the stated tolerances on MSVC; CI on ubuntu/macos passes with the documented relative tolerances (analytic values must not be bit-pinned across libms).
- [ ] Existing planar/cylinder curved boolean tests, `split_solid` tests and trimmed tessellation tests unchanged.
- [ ] ADR-0028 and §16.3 describe exactly the supported surface/operation matrix and the remaining limits.

## Self-check

- [ ] Every task names exact files, interfaces and a runnable failing test.
- [ ] No task depends on code from a later task (Task 5 depends on 1-4; Task 6 and 7 depend on 1 and 3).
- [ ] No placeholders or TBDs.
- [ ] QA values are independent analytic formulas, not values read from the implementation.
- [ ] The plan states the deviations from the architecture (closed cuts, synthetic boundary, chordal torus, M3c-3 deferral) explicitly.
