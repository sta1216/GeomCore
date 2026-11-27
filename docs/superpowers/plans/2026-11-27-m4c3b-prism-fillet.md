# Convex Prism Vertical-Edge Fillet Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use subagent-driven-development (recommended) or executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Extend `fillet_edges` from axis-aligned boxes to convex prisms: round any subset of the vertical edges of a convex polygon extruded along one axis, with exact cylinders and planar caps bounded by circular arcs.

**Architecture:** A prism is recognised from the existing `Brep`: two parallel congruent convex cap faces, perpendicular side faces, and one vertical edge per polygon vertex. The fillet trims each selected vertex back by the tangent length `r * cot(theta/2)`, joins the two tangent points with a circular arc centred on the offset vertex, and builds one exact `Cylinder3` patch per selected edge. Caps, sides and cylinders are assembled directly with `assemble_curved_solid`; no Boolean, splitter, or tolerance surgery is involved.

**Tech Stack:** C++20, CMake 4.4.4, MSVC 19.44 (`/W4 /WX /utf-8`), doctest 2.4.11, existing GeomCore headers only.

**Spec:** `docs/architecture/geometry-library-architecture.md` v2.1 搂8.9 (blend engine, rolling-ball approximation, failure must not produce a partial solid), 搂16.4 (milestone evidence), `docs/quality/testing-strategy.md`, and ADR-0042 (M4c-3a box fillet). This plan is M4c-3b, the second slice of 搂8.9.

## Global Constraints

- C++20; no new third-party dependency.
- Warnings are errors: `-DGEOM_WERROR=ON`, `/W4 /WX /utf-8`; every new header compiles standalone.
- No global epsilon. Geometry decisions use exact predicates where a sign is required (convexity/orientation); distances use the caller's positive `tolerance`.
- One solid, one closed shell, planar faces only, convex cap polygon. Reject everything else with a stable `Errc` and message; never return a partial or silently wrong solid.
- Input `Brep` must be byte-for-byte unchanged on success and failure; the output is deterministic for the same input.
- Exact radius guard: tangent length `t = radius / tan(theta/2)`. The remaining straight length of every side is `L - t_i - t_{i+1}` and must exceed `tolerance`.
- Test-driven: every task writes the failing test first, shows the RED output, then implements.
- One commit per task; run Debug and Release; run only the filtered suites while iterating, then the full Release suite once per task or at the end of a task that touches shared headers.

## Review Focus

These are the inputs most likely to break a reasonable user, and each one has a test in the task that owns the code:

1. **Prism misdetection.** A tapered or oblique solid whose side faces are not parallel to the found axis, or a cap polygon that is not congruent, must be rejected as unsupported, not filleted into a wrong solid. The box path stays first and unchanged.
2. **Tangent overrun on a short side.** `t = r*cot(theta/2)` can exceed a short polygon side even when `r` is below the inradius; the guard must reject before assembly with a stable message.
3. **Horizontal or dead selected edge.** A valid but non-vertical prism edge must be `unsupported`, a dead `EdgeId` must be `index_out_of_range`, and the order/duplicates of the selected list must not change the result.
4. **Tilted prism and cap orientation.** A prism whose axis is not `+Z` must give the same volume/topology; top and bottom cap arcs must wind consistently about their own outward normals.
5. **Input immutability and determinism.** Volume, counts, topology and tessellation must be identical across two runs, and the input solid must be unchanged after success and after a rejected radius.

---

## File Map

- Create `include/geom/brep/detail/fillet_prism.hpp` 鈥?prism recognition, cap-ring construction, side/cylinder fragments, `fillet_prism_edges`.
- Modify `include/geom/brep/fillet.hpp` 鈥?dispatch box first, then prism; move the existing box body into a detail entry point only if needed by the dispatch.
- Create `tests/unit/fillet_prism_test.cpp` 鈥?prism fixture and all cases below.
- Modify `tests/CMakeLists.txt` 鈥?register the new test source.
- Modify `CMakeLists.txt` 鈥?add `geom_bench_m4c3b` in Task 4.
- Create `bench/bench_m4c3b.cpp` 鈥?Task 4.
- Create `docs/adr/0043-prism-fillets.md`, `docs/benchmarks/m4c3b-fillet.md` 鈥?Task 4.
- Modify `docs/architecture/geometry-library-architecture.md` 搂8.9/搂16.4 and this plan's Results 鈥?Task 4.
- Modify `.github/workflows/ci.yml` 鈥?build and run the new bench 鈥?Task 4.

---

### Task 1: Recognise a convex prism and round one vertical edge

**Files:**
- Create: `include/geom/brep/detail/fillet_prism.hpp`
- Modify: `include/geom/brep/fillet.hpp` (dispatch and, if needed, extract the box body)
- Create: `tests/unit/fillet_prism_test.cpp`
- Modify: `tests/CMakeLists.txt`

**Interfaces:**
- Consumes: `Brep`, `Face`, `Edge`, `EdgeId`, `face_loops_uv`, `face_uv_eval`, `make_plane_frame`, `plane_frame_to_2d`, `sign_orient2d`, `assemble_curved_solid`, `assemble_solid`, `CurvedAssembledFace`, `CurvedArcSpan`, `CircleCurve3`, `Cylinder3`, `Plane3`.
- Produces:
  - `struct FilletPrism { geom::Point3d bottom; geom::Vec3d axis; double height; std::vector<geom::Point3d> ring; std::vector<geom::EdgeId> vertical; };` in `geom::detail`.
  - `Result<FilletPrism> fillet_prism_bounds(const Brep&, double tolerance);`
  - `Result<Brep> fillet_prism_edges(const FilletPrism&, std::span<const std::uint32_t> vertices, double radius, double tolerance);` 鈥?`vertices[i]` indexes `prism.ring`/`prism.vertical`; duplicates are the caller's responsibility in this task (Task 2 normalises them).
  - `fillet_edges` keeps its exact current signature and messages for boxes. For a non-box it calls `fillet_prism_bounds`; if both recognisers fail it returns `Error{Errc::unsupported, "fillet_edges: only a box or a convex prism is supported"}`.

- [ ] **Step 1: Add the failing prism test**

Create `tests/unit/fillet_prism_test.cpp` with the fixture and one case. The fixture must build the prism from public geometry only, so a broken recogniser cannot make the test pass by construction.

```cpp
#include <doctest/doctest.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <map>
#include <numbers>
#include <span>
#include <utility>
#include <vector>

#include <geom/brep/assemble.hpp>
#include <geom/brep/fillet.hpp>
#include <geom/brep/detail/curved_volume.hpp>
#include <geom/brep/tessellate.hpp>
#include <geom/brep/topology.hpp>
#include <geom/brep/validate.hpp>
#include <geom/cgeom/polygon2.hpp>
#include <geom/core/error.hpp>
#include <geom/core/result.hpp>
#include <geom/curves/analytic_curve.hpp>
#include <geom/primitives/point.hpp>
#include <geom/primitives/vec.hpp>

namespace {

constexpr double kTolerance = 1e-9;
constexpr double kPi = std::numbers::pi_v<double>;

std::vector<geom::Point2d> triangle() {
  return {geom::Point2d{0.0, 0.0}, geom::Point2d{2.0, 0.0},
          geom::Point2d{1.0, std::sqrt(3.0)}};
}

// CCW polygon in the z = [z0, z1] extrusion. Bottom is clockwise about -Z;
// sides are outward-facing quads. The topology is built by assemble_solid.
geom::Result<geom::Brep> make_prism(std::span<const geom::Point2d> polygon, double z0, double z1,
                                   double tolerance) {
  std::vector<geom::AssembledFace> faces;
  std::vector<geom::Point3d> bottom;
  std::vector<geom::Point3d> top;
  for (const geom::Point2d& p : polygon) {
    bottom.push_back(geom::Point3d{p[0], p[1], z0});
    top.push_back(geom::Point3d{p[0], p[1], z1});
  }

  geom::AssembledFace bottom_face;
  bottom_face.outer.assign(bottom.rbegin(), bottom.rend());
  bottom_face.normal = geom::Vec3d{0.0, 0.0, -1.0};
  faces.push_back(std::move(bottom_face));

  geom::AssembledFace top_face;
  top_face.outer = top;
  top_face.normal = geom::Vec3d{0.0, 0.0, 1.0};
  faces.push_back(std::move(top_face));

  for (std::size_t i = 0; i < polygon.size(); ++i) {
    const std::size_t j = (i + 1) % polygon.size();
    const geom::Vec3d edge{polygon[j][0] - polygon[i][0], polygon[j][1] - polygon[i][1], 0.0};
    const geom::Vec3d normal = geom::detail::assemble_unit(
        geom::Vec3d{edge[1], -edge[0], 0.0});
    geom::AssembledFace side;
    side.normal = normal;
    side.outer = {bottom[i], bottom[j], top[j], top[i]};
    faces.push_back(std::move(side));
  }
  return geom::assemble_solid(faces, tolerance);
}

std::vector<geom::EdgeId> vertical_edges(const geom::Brep& brep) {
  std::vector<geom::EdgeId> result;
  for (const geom::EdgeId edge_id : brep.edge_ids()) {
    const geom::Edge& edge = brep.edge(edge_id);
    const geom::Point3d& a = brep.vertex(edge.v0).point;
    const geom::Point3d& b = brep.vertex(edge.v1).point;
    if (std::hypot(b[0] - a[0], b[1] - a[1]) <= 1e-9 && std::abs(b[2] - a[2]) > 1e-9) {
      result.push_back(edge_id);
    }
  }
  return result;
}

geom::EdgeId vertical_edge_at(const geom::Brep& brep, geom::Point2d xy) {
  for (const geom::EdgeId edge_id : vertical_edges(brep)) {
    const geom::Edge& edge = brep.edge(edge_id);
    const geom::Point3d& a = brep.vertex(edge.v0).point;
    if (std::hypot(a[0] - xy[0], a[1] - xy[1]) <= 1e-9) return edge_id;
  }
  return geom::EdgeId{};
}

double corner_area(double radius, double theta) {
  return radius * radius *
         (1.0 / std::tan(theta / 2.0) - (kPi - theta) / 2.0);
}

double volume_of(const geom::Brep& brep) {
  return geom::detail::curved_solid_volume(brep, brep.solid_ids().front(), 1e-9).value();
}

}  // namespace

TEST_CASE("fillet_edges rounds one vertical edge of a triangular prism") {
  const geom::Result<geom::Brep> prism = make_prism(triangle(), 0.0, 3.0, kTolerance);
  REQUIRE(prism.has_value());
  const geom::Brep& solid = prism.value();
  const std::vector<geom::EdgeId> edges{vertical_edge_at(solid, geom::Point2d{0.0, 0.0})};
  REQUIRE(edges.front() != geom::EdgeId{});

  const geom::Result<geom::FilletResult> result =
      geom::fillet_edges(solid, edges, geom::FilletOptions{0.2, kTolerance});
  REQUIRE(result.has_value());
  const geom::Brep& rounded = result.value().brep;
  CHECK(result.value().filleted_edges == 1u);
  CHECK(geom::validate_topology(rounded).ok());
  CHECK(rounded.solid_ids().size() == 1u);
  CHECK(rounded.shell_count() == 1u);

  const double expected =
      3.0 * std::sqrt(3.0) - 3.0 * corner_area(0.2, kPi / 3.0);
  CHECK(volume_of(rounded) == doctest::Approx(expected).epsilon(1e-9));
  CHECK(rounded.vertex_count() == 8u);
  CHECK(rounded.edge_count() == 12u);
  CHECK(rounded.face_count() == 6u);
}
```

Register the file in `tests/CMakeLists.txt` next to the other `unit/` sources, at the end of `GEOM_TEST_SOURCES`.

- [ ] **Step 2: Run the test and verify it fails**

Run: `cmake --build build --config Debug --target geom_tests` then
`build\tests\Debug\geom_tests.exe --test-case="*prism*"`.
Expected: compile error `C1083: Cannot open include file: 'geom/brep/detail/fillet_prism.hpp'` or the public dispatch returns `fillet_edges: only an axis-aligned box is supported`.

- [ ] **Step 3: Implement `fillet_prism_bounds` and `fillet_prism_edges`**

In `include/geom/brep/detail/fillet_prism.hpp`, implement the recognition:

- Require one solid, one closed shell, every face planar, one loop per face, `3 <= n <= 64`.
- Candidate axis: for each unique face normal `n` (normalised), count faces with `|dot(normal,n)| > 1 - 1e-9`; require exactly `2` and all other faces perpendicular (`|dot| <= 1e-9`). The two cap normals must be opposite. Try every face normal until one works; if none, return `unsupported "fillet_edges: only a box or a convex prism is supported"`.
- Bottom cap is the cap whose outward normal has `dot(normal, axis) < 0`; top cap is the other. Read each cap loop in 3D with `face_loops_uv` + `face_uv_eval`, require one loop each.
- Project onto an orthonormal basis `(u,v)` perpendicular to `axis`; reverse the ring when `ring_orientation2` is `-1`, then require it to be strictly convex and CCW (`sign_orient2d(prev, cur, next) > 0` for every triple; reject collinear or reflex).
- Cyclically align the top ring to the bottom ring by minimising the maximum distance; require every matched pair within `64*tolerance` and `top[i] - bottom[i] = axis*height` with `height > 0`.
- For each polygon vertex `i`, find exactly one edge whose endpoints match `bottom[i]` and `top[i]` within `tolerance`; store it in `vertical[i]`. Reject when missing or duplicated.
- `FilletPrism{bottom[0], axis, height, bottom_ring, vertical}`. `ring` is the CCW bottom ring.

For one selected vertex `i`:

- Compute the interior angle `theta` from the incoming and outgoing edges (must satisfy `0 < theta < pi`).
- `t = radius / tan(theta/2)`; tangent points `T_in = V - unit(V - prev)*t` and `T_out = V + unit(next - V)*t`.
- Solve the two inward-offset line equations for `A`: `n_in . p = n_in . V + radius`, `n_out . p = n_out . V + radius`, where each `n` is the unit inward normal of that edge (choose the sign whose dot with the polygon centroid is positive).
- Arc: `CircleCurve3{A, unit(T_in - A), cross(axis, unit(T_in - A)), radius}`; sweep `= pi - theta`; `CurvedArcSpan{circle, 0.0, sweep}` on the segment `T_in -> T_out`.
- Cap ring: walk vertices CCW; for a selected vertex push `T_in` and `T_out` with the arc span between them; otherwise push `V` with `std::nullopt`. Do the same on the top cap at `+axis*height`. The bottom cap loop is reversed (points and spans) so its loop winds CCW about `-axis`.
- Side face `i`: bottom points `T_out(i)` (or `V_i`) and `T_in(i+1)` (or `V_{i+1})`, top points `+ axis*height`; plane normal `= unit(cross(bottom_edge, axis))`; plane offset `= dot(normal, point)`; UV ring from `make_plane_frame`.
- Cylinder face: `Cylinder3{A, axis, unit(T_in - A), radius}`, points `{T_in(bottom), T_out(bottom), T_out(top), T_in(top)}`, arc spans `{bottom_arc, nullopt, top_arc, nullopt}`, UV `{(0,0),(sweep,0),(sweep,height),(0,height)}` (follow the exact orientation pattern of `detail/fillet_box.hpp:196-226`).
- Assemble with `assemble_curved_solid(faces, tolerance)`; on failure return `invalid_topology "fillet_edges: fillet could not be assembled"`.

In `fillet.hpp`, after the options checks, try `fillet_box_bounds`; if it succeeds run the existing box body unchanged. Otherwise try `fillet_prism_bounds`; on success:

- reject an empty selection with the existing `"fillet_edges: at least one edge is required"`;
- for each selected `EdgeId`, if dead return `index_out_of_range "fillet_edges: edge does not belong to the solid"`; if it is not in `prism.vertical` return `unsupported "fillet_edges: only vertical prism edges are supported"`;
- deduplicate indices, call `fillet_prism_edges`, and return `FilletResult{brep, unique_count}`.

- [ ] **Step 4: Run the filtered tests in both configurations**

Run: `build\tests\Debug\geom_tests.exe --test-case="*prism*"` and the Release binary.
Expected: 2 cases, 22 assertions, 0 failed; the existing `--test-case="*fillet*"` suite passes 9/898.

- [ ] **Step 5: Commit**

```bash
git add include/geom/brep/detail/fillet_prism.hpp include/geom/brep/fillet.hpp tests/unit/fillet_prism_test.cpp tests/CMakeLists.txt
git commit -m "feat(brep): fillet one vertical edge of a convex prism"
```

---

### Task 2: Round any subset of the vertical edges

**Files:**
- Modify: `include/geom/brep/detail/fillet_prism.hpp`
- Modify: `tests/unit/fillet_prism_test.cpp`

**Interfaces:**
- Consumes: Task 1's `FilletPrism`, `fillet_prism_edges`, cap/side/cylinder builders.
- Produces: the same `fillet_prism_edges` signature, now accepting any subset of `0..n-1`, with duplicates rejected or collapsed before the call.

- [ ] **Step 1: Add failing subset/all cases**

Append to `tests/unit/fillet_prism_test.cpp`:

```cpp
std::vector<geom::Point2d> regular_polygon(int n, double radius) {
  std::vector<geom::Point2d> polygon;
  for (int i = 0; i < n; ++i) {
    const double angle = 2.0 * kPi * static_cast<double>(i) / static_cast<double>(n);
    polygon.push_back(geom::Point2d{radius * std::cos(angle), radius * std::sin(angle)});
  }
  return polygon;
}

TEST_CASE("fillet_edges rounds every vertical edge of a triangular prism") {
  const geom::Result<geom::Brep> prism = make_prism(triangle(), 0.0, 3.0, kTolerance);
  REQUIRE(prism.has_value());
  const std::vector<geom::EdgeId> edges = vertical_edges(prism.value());
  REQUIRE(edges.size() == 3u);
  const geom::Result<geom::FilletResult> result =
      geom::fillet_edges(prism.value(), edges, geom::FilletOptions{0.2, kTolerance});
  REQUIRE(result.has_value());
  CHECK(result.value().filleted_edges == 3u);
  CHECK(volume_of(result.value().brep) ==
        doctest::Approx(3.0 * std::sqrt(3.0) - 9.0 * corner_area(0.2, kPi / 3.0))
            .epsilon(1e-9));
  CHECK(result.value().brep.vertex_count() == 12u);
  CHECK(result.value().brep.edge_count() == 18u);
  CHECK(result.value().brep.face_count() == 8u);
}

TEST_CASE("fillet_edges rounds an adjacent pair on a hexagonal prism") {
  const std::vector<geom::Point2d> polygon = regular_polygon(6, 1.0);
  const geom::Result<geom::Brep> prism = make_prism(polygon, 0.0, 2.0, kTolerance);
  REQUIRE(prism.has_value());
  std::vector<geom::EdgeId> edges = vertical_edges(prism.value());
  REQUIRE(edges.size() == 6u);
  const std::vector<geom::EdgeId> selected{edges[0], edges[1]};
  const geom::Result<geom::FilletResult> result =
      geom::fillet_edges(prism.value(), selected, geom::FilletOptions{0.2, kTolerance});
  REQUIRE(result.has_value());
  const double area = std::abs(geom::signed_area2(polygon));
  const double expected =
      area * 2.0 - 2.0 * 2.0 * corner_area(0.2, 2.0 * kPi / 3.0);
  CHECK(volume_of(result.value().brep) == doctest::Approx(expected).epsilon(1e-9));
  CHECK(result.value().brep.vertex_count() == 16u);
  CHECK(result.value().brep.edge_count() == 24u);
  CHECK(result.value().brep.face_count() == 10u);
}
```

Add the includes `<geom/cgeom/polygon2.hpp>` (already present) and use `geom::signed_area2`.

- [ ] **Step 2: Run and verify RED**

Run the two new cases.
Expected: the first returns only one rounded edge (or fails the count/volume), the second fails because only the first index builds a corner.

- [ ] **Step 3: Generalise the builders**

- Normalise the selected vertex list: sort, unique; reject an empty list before calling.
- Build `selected[n]` flags. Every cap ring uses the flag per vertex; every side face trims its two ends by `t_i` only when that end is selected.
- Emit one cylinder per selected vertex.
- Guard each side before assembly: `L - (selected_i ? t_i : 0) - (selected_{i+1} ? t_{i+1} : 0) > tolerance`; if not, return `unsupported "fillet_edges: radius leaves a prism face narrower than the tolerance"`.
- Keep the Task 1 exact volume/lexicographic assumptions; the code must not special-case `k == 1`.

- [ ] **Step 4: Run filtered tests**

Run `--test-case="*prism*"` in Debug and Release.
Expected: 3 cases, with the triangle one-edge case still 8/12/6 and the full/partial counts above.

- [ ] **Step 5: Commit**

```bash
git add include/geom/brep/detail/fillet_prism.hpp tests/unit/fillet_prism_test.cpp
git commit -m "feat(brep): fillet any vertical edge subset of a convex prism"
```

---

### Task 3: Diagnostics, tilted prism, immutability and watertightness

**Files:**
- Modify: `include/geom/brep/detail/fillet_prism.hpp`
- Modify: `tests/unit/fillet_prism_test.cpp`

**Interfaces:**
- Consumes: everything from Tasks 1鈥?.
- Produces: stable rejection messages; no API change.

- [ ] **Step 1: Add failing diagnostic cases**

Append cases that assert exact messages:

```cpp
TEST_CASE("fillet_edges rejects a non-prism with the combined message") {
  const geom::Result<geom::Brep> cylinder = geom::make_cylinder(
      geom::Point3d{0, 0, 0}, geom::Vec3d{0, 0, 1}, 1.0, 2.0, kTolerance);
  REQUIRE(cylinder.has_value());
  const std::vector<geom::EdgeId> edges{cylinder.value().edge_ids().front()};
  const geom::Result<geom::FilletResult> result =
      geom::fillet_edges(cylinder.value(), edges, geom::FilletOptions{0.2, kTolerance});
  REQUIRE(!result.has_value());
  CHECK(result.error().code == geom::Errc::unsupported);
  CHECK(result.error().message == "fillet_edges: only a box or a convex prism is supported");
}

TEST_CASE("fillet_edges rejects a horizontal prism edge") {
  const geom::Result<geom::Brep> prism = make_prism(triangle(), 0.0, 3.0, kTolerance);
  REQUIRE(prism.has_value());
  geom::EdgeId horizontal{};
  for (const geom::EdgeId edge_id : prism.value().edge_ids()) {
    const geom::Edge& edge = prism.value().edge(edge_id);
    if (std::abs(prism.value().vertex(edge.v0).point[2] - 3.0) <= kTolerance &&
        std::abs(prism.value().vertex(edge.v1).point[2] - 3.0) <= kTolerance) {
      horizontal = edge_id;
      break;
    }
  }
  REQUIRE(horizontal != geom::EdgeId{});
  const std::vector<geom::EdgeId> edges{horizontal};
  const geom::Result<geom::FilletResult> result =
      geom::fillet_edges(prism.value(), edges, geom::FilletOptions{0.2, kTolerance});
  REQUIRE(!result.has_value());
  CHECK(result.error().code == geom::Errc::unsupported);
  CHECK(result.error().message == "fillet_edges: only vertical prism edges are supported");
}

TEST_CASE("fillet_edges rejects a radius that leaves a prism face too narrow") {
  const geom::Result<geom::Brep> prism = make_prism(triangle(), 0.0, 3.0, kTolerance);
  REQUIRE(prism.has_value());
  const std::vector<geom::EdgeId> edges = vertical_edges(prism.value());
  const geom::Result<geom::FilletResult> result =
      geom::fillet_edges(prism.value(), edges, geom::FilletOptions{1.0, kTolerance});
  REQUIRE(!result.has_value());
  CHECK(result.error().code == geom::Errc::unsupported);
  CHECK(result.error().message ==
        "fillet_edges: radius leaves a prism face narrower than the tolerance");
}
```

Also add: dead `EdgeId` 鈫?`index_out_of_range`; empty span 鈫?`"fillet_edges: at least one edge is required"`; duplicate/order invariance; a short-side radius 鈫?`"radius leaves a prism face narrower than the tolerance"`.

- [ ] **Step 2: Run and verify RED**

Expect the non-prism case to leak the old box-only message, and the radius case to reach assembly or return a wrong code.

- [ ] **Step 3: Implement the guards and the tilted fixture**

- Do not impose a global inradius limit: a single selected corner may be filleted with `r` above the polygon inradius as long as its tangent length fits both adjacent sides. The per-vertex tangent check (`t <` both adjacent side lengths) returns `unsupported "fillet_edges: radius exceeds the prism cross-section"`; the every-side width guard returns `"fillet_edges: radius leaves a prism face narrower than the tolerance"`.
- Add a positive partial case: triangle prism `h=3`, only the edge at `(0,0)`, `r=0.6` (above the inradius 0.577), volume `4.45651090582474`, V/E/F 8/12/6.
- For every side, apply the Task 2 width guard before building any face.
- Map an assembly failure to `invalid_topology "fillet_edges: fillet could not be assembled"`; never leak `assemble_curved_solid` internals.
- Add `make_prism_in_frame(polygon2d, origin, axis, u, tolerance)`: build the same faces in the `(u, cross(axis,u))` plane and assemble. Add one case with `axis = normalize({1,1,1})`, all three triangle edges, `r = 0.2`, comparing volume to the same `3*sqrt(3) - 9*corner_area` anchor and V/E/F 12/18/8.
- Add the immutability case: snapshot `vertex_count`, `edge_count`, `face_count`, `volume_of(input)` before and after a successful call and a rejected radius; require exact equality. Call the same input twice and require exact volume/count equality and identical tessellation `positions`/`triangles` vectors.

- [ ] **Step 4: Run Debug and Release filtered suites**

Run `--test-case="*prism*"` and `--test-case="*fillet*"` in both configurations.
Expected: all prism cases pass; box fillet remains 8/887 (or higher only if this task adds a case to that file, which it must not).

- [ ] **Step 5: Commit**

```bash
git add include/geom/brep/detail/fillet_prism.hpp tests/unit/fillet_prism_test.cpp
git commit -m "test(brep): pin prism fillet diagnostics and invariants"
```

---

### Task 4: Benchmark, ADR, architecture and Results

**Files:**
- Create: `bench/bench_m4c3b.cpp`
- Modify: `CMakeLists.txt` (target `geom_bench_m4c3b`, guarded by `GEOM_BUILD_BENCH`)
- Create: `docs/adr/0043-prism-fillets.md`
- Create: `docs/benchmarks/m4c3b-fillet.md`
- Modify: `docs/architecture/geometry-library-architecture.md` 搂8.9 and 搂16.4
- Modify: `.github/workflows/ci.yml`
- Modify: this plan's Results section

**Interfaces:**
- Consumes: `fillet_edges`, `make_prism`-equivalent in the bench, `curved_solid_volume`, `validate_topology`.
- Produces: reproducible timing and volume evidence; no library API change.

- [ ] **Step 1: Add the benchmark**

`bench/bench_m4c3b.cpp` builds a 64-gon prism (`radius 2`, height 10) with the same public `assemble_solid` fixture, fillets all 64 vertical edges with `radius 0.05`, checks `validate_topology`, one solid/one shell and volume against `area*10 - 64*10*corner_area(0.05, 2*pi/64)`, then times 5 runs and prints `min/median/max` plus the volume. Exit 2 on any mismatch or validation failure, like the other benches.

Add to the root `CMakeLists.txt`:

```cmake
if(GEOM_BUILD_BENCH)
  add_executable(geom_bench_m4c3b bench/bench_m4c3b.cpp)
  target_link_libraries(geom_bench_m4c3b PRIVATE geom::core geom::warnings)
  target_compile_features(geom_bench_m4c3b PRIVATE cxx_std_20)
endif()
```

- [ ] **Step 2: Build and run the bench**

Configure with `-DGEOM_BUILD_BENCH=ON`, build Release, run `geom_bench_m4c3b.exe`.
Expected: exit 0, volume within `1e-9` relative, timings printed.

- [ ] **Step 3: Write the evidence docs**

- ADR-0043: decision (direct prism fragments, exact cylinders, cap arcs), scope (convex prism, vertical edges, equal radius, subset), rejected alternatives (generic sweep/vertex sphere patch), failure contract, limitations.
- `docs/benchmarks/m4c3b-fillet.md`: machine, compiler, workload, min/median/max, volume, reproduction command.
- 搂8.9: append the M4c-3b contract (prism only; adjacent vertex blends/variable radius remain M4c-3c).
- 搂16.4: evidence line with commit range and counts.
- Fill this plan's Results table with measured values; do not claim a PASS for anything not run.

- [ ] **Step 4: Full Release suite**

Run the complete Release suite (`geom_tests.exe`) and the new bench.
Expected: 0 failures; record exact case/assertion counts. Debug full suite is not a gate; run the Debug `*prism*` and `*fillet*` lanes.

- [ ] **Step 5: Commit**

```bash
git add bench/bench_m4c3b.cpp CMakeLists.txt docs/adr/0043-prism-fillets.md docs/benchmarks/m4c3b-fillet.md docs/architecture/geometry-library-architecture.md docs/superpowers/plans/2026-11-27-m4c3b-prism-fillet.md .github/workflows/ci.yml
git commit -m "docs(brep): record the M4c-3b prism fillet contract and benchmark"
```

---

## Results

| Item | Value |
|---|---|
| Triangle one edge, r=0.2, h=3 | expected `5.113970031941959`; actual `5.1139700319419079` (relative `1.0e-14`) |
| Triangle all edges | expected `4.949605250412612`; actual `4.9496052504125752` (relative `7.4e-15`) |
| Hexagon adjacent pair | expected `5.187552183732021`; actual `5.18755218373202` |
| Tilted triangle all edges | expected same as axis-aligned triangle; actual `4.9496052504125903` (same anchor, relative `4.3e-15`) |
| Debug/Release `*prism*` cases/assertions | Debug 11/98, Release 11/98, 0 failed (the filter matches one pre-existing chamfer case and misses three prism cases whose names lack "prism"; the prism file itself adds 13 cases/127 assertions) |
| Full Release cases/assertions | 1,172 cases / 683,738 assertions, 0 failed, 956 s (baseline before M4c-3b: 1,159/683,611; the prism file adds 13/127) |
| Bench volume / min-median-max | `125.461876478546259` vs analytic `125.461876478546088`; `11.683 / 14.477 / 15.385 ms` (five timed runs, timed region is `fillet_edges` only) |
| Commit range | code under test `2f5a258..b77cd88` on `sdd/m4c3b`; the Task 4 benchmark, docs and CI commit sits on top |

Deviations recorded during execution:

- The in-tree target is `geom_core` (which already links `geom_warnings`); `geom::core`/`geom::warnings` exist only as installed-package export names, so the bench links `geom_core` and requests `cxx_std_20` explicitly. Using the plan's literal names would fail configure.
- `corner_area` takes the polygon interior angle, as in `tests/unit/fillet_prism_test.cpp`; the expected volume therefore uses `theta = pi - 2*pi/64`, not the literal `2*pi/64` from the brief. The latter would give ~95.3 instead of ~125.46 and fail the check.
- The brief's "11 cases/98 assertions" is a filter artifact: `*prism*` also matches `chamfer one box edge removes a triangular prism` and misses three prism test cases whose names do not contain "prism". The file adds 13 cases/127 assertions, and the full Release count moves from 1,159/683,611 to 1,172/683,738.
- Debug full suite was not run (not a gate); the Debug `*prism*` and `*fillet*` lanes were run and passed.

## Self-Review Notes

- The plan uses only existing public/detail APIs; `assemble_solid` and `assemble_curved_solid` are already proven in M2/M4.
- The radius guard uses the exact tangent length `r*cot(theta/2)`, not `r`, so acute polygon corners are handled.
- Convexity and ring orientation use `sign_orient2d`; no new floating orientation test.
- Deliberately out of scope: adjacent non-vertical edges, spherical vertex blends, variable radius, non-convex polygons, multi-shell solids, and general curved prisms. These stay in M4c-3c.




