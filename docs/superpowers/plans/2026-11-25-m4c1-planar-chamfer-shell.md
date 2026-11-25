# Planar Chamfer and Shell (M4c-1) Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use subagent-driven-development (recommended) or executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Add exact constant-distance chamfer and uniform-thickness shell operations for convex, planar-faced, single-shell B-rep solids, with deterministic diagnostics for every unsupported input.

**Architecture:** Both operations reduce to repeated half-space cuts. Chamfering an edge cuts the solid with the plane through the two face setback points and keeps the side away from the edge; shelling cuts the solid with every face plane shifted inward and subtracts the resulting inner solid from the original. A requested open face is shifted **outward** so the inner solid protrudes through the opening and the difference cuts a clean pocket instead of leaving two coplanar caps. Cutting reuses `split_solid_by_plane`; the only new geometry code is frame construction, convexity validation, and distance-vs-face-width checks. No new tolerance regime and no global epsilon are introduced.

**Tech Stack:** C++20, header-only `geom_core`, `geom::v1`, doctest 2.4.11, CMake 4.x, MSVC 19.44 `/W4 /WX /utf-8`, PowerShell on Windows.

**Spec:** `docs/architecture/geometry-library-architecture.md` §8.9 (blend operations) and §16.4 (milestone evidence); `docs/quality/testing-strategy.md` (test classes, negative paths, determinism).

## Global Constraints

- C++20, header-only: every new function is `inline` in a header under `include/geom/brep/`; namespace `geom { inline namespace v1 { ... } }`.
- Public entry points return `Result<T>`; no exceptions cross the API; every failure carries a stable `Errc` and an exact message prefix.
- Inputs are `const Brep&` and never mutated. Failures leave the caller's document untouched and return no partial solid.
- Tolerances are explicit arguments; no global epsilon. A geometric guard uses `max(64 * tolerance, 1e-12 * scale)` and the plan states it where used.
- The only accepted operands in this milestone are a **single solid with a single closed shell, all faces planar, and the solid convex**. Everything else is `Errc::unsupported` with a diagnostic naming the reason.
- The only accepted edges are straight `LineCurve3` edges with exactly two coedges on two distinct faces.
- Tests use doctest, fixed literals, and independent expected values (analytic volumes), not values copied from the implementation. Every public error path has a negative test.
- Build with `-DGEOM_WERROR=ON`, run both Debug and Release, and keep `/W4 /WX /utf-8` clean on MSVC.
- One task = one commit. Record Debug and Release test counts in the task report.
- ASCII only in source and tests; documents may use the repository's existing style.

## Scope and Deferrals

This plan deliberately splits M4's blend work. It ships exact chamfer and shell for convex planar solids. It does **not** ship fillet, because an exact rolling-ball fillet needs an analytic cylindrical patch plus vertex-blend surgery, and a faceted approximation would violate the library's "no half-finished solid" rule. Fillet (and variable-radius blends, chamfer-by-angle, non-convex solids, curved faces, multi-shell inputs, and offset/sweep) is M4c-2 and needs its own plan.

The convex restriction is a correctness boundary, not an implementation shortcut: with all edges convex, every selected chamfer plane cuts off exactly the corner cap it names, and the inner half-space intersection is exactly the eroded solid. A non-convex solid would make a single plane cut remove material elsewhere, so it is rejected rather than guessed at.

## Review Focus

1. **Chamfer distance at or beyond the adjacent face width** — a caller can ask for `distance = 2.0` on a 2-wide face. Expected: `Errc::unsupported` naming the width limit, never a closed solid with a collapsed or flipped face.
2. **Three chamfers meeting at one vertex** — sequential half-space cuts must leave a closed solid with a planar vertex facet, `chi == 2`, and no non-manifold edge. Expected: valid solid and the analytic volume for a fully chamfered cube.
3. **Open-face shell** — the inner solid is shifted outward through the omitted face. Expected: one open pocket, the omitted face keeps an outer and an inner loop, and the tessellated mesh has zero boundary edges (no coincident cap, no zero-thickness wall).
4. **Wall thickness at or above the local thickness** — the inner solid vanishes. Expected: `Errc::degenerate` ("wall thickness leaves no interior"), never an empty Brep or a self-intersecting one.
5. **Curved or non-convex operand** — a cylinder or an L-shaped union is easy to pass by mistake. Expected: `Errc::unsupported` before any cutting, with the input counts unchanged.

## File Map

- Create `include/geom/brep/detail/blend_common.hpp` — shape/convexity validation, outward plane frames, straight-edge and inward-direction helpers, face-width ray query.
- Create `include/geom/brep/chamfer.hpp` — `ChamferOptions`, `ChamferResult`, `chamfer_edges`.
- Create `include/geom/brep/shell.hpp` — `ShellOptions`, `ShellResult`, `shell_solid`.
- Create `tests/unit/chamfer_test.cpp` — chamfer geometry, multi-edge and diagnostic tests.
- Create `tests/unit/shell_test.cpp` — closed, open and diagnostic shell tests.
- Create `tests/unit/blend_qa_test.cpp` — watertight tessellation, determinism and input-immutability cross-checks.
- Create `bench/bench_m4c1.cpp` — chamfer-all-edges and closed-shell timings plus volume checksums.
- Modify `tests/CMakeLists.txt` — register the three test files.
- Modify `CMakeLists.txt` — register `geom_bench_m4c1` under `GEOM_BUILD_BENCH`.
- Create `docs/benchmarks/m4c1-blends.md` — measured release numbers.
- Create `docs/adr/0041-planar-blends-via-half-space-cuts.md` — decision and limits.
- Modify `docs/architecture/geometry-library-architecture.md` — §8.9 scope line and §16.4 M4c-1 evidence.
- Modify this plan — tick steps and fill the Results table.

---

### Task 1: Shared planar-convex validation and single-edge chamfer

**Files:**
- Create: `include/geom/brep/detail/blend_common.hpp`
- Create: `include/geom/brep/chamfer.hpp`
- Test: `tests/unit/chamfer_test.cpp`
- Modify: `tests/CMakeLists.txt`

**Interfaces:**
- Consumes: `Brep`, `SolidId`, `FaceId`, `EdgeId`, `Plane3`, `LineCurve3`, `split_solid_by_plane(const Brep&, SolidId, const Plane3&, double)` returning `SplitSolidResult{below, above}`, `validate_topology`, `mass_properties`.
- Produces:
  - `struct detail::BlendPlanarSolid { SolidId solid; std::vector<FaceId> faces; std::vector<Vec3d> outward; std::vector<double> offset; double scale; };`
  - `Result<detail::BlendPlanarSolid> detail::blend_planar_convex_solid(const Brep&, double tolerance, std::string_view caller);`
  - `Result<Vec3d> detail::blend_edge_direction(const Brep&, EdgeId, std::string_view caller);`
  - `Result<Vec3d> detail::blend_inward_direction(const Brep&, FaceId, EdgeId, const Vec3d& edge_dir, const Vec3d& edge_mid, std::string_view caller);`
  - `Result<double> detail::blend_face_width(const Brep&, FaceId, const Vec3d& origin, const Vec3d& direction, std::string_view caller);`
  - `struct ChamferOptions { double distance = 0.0; double tolerance = 1e-9; };`
  - `struct ChamferResult { Brep brep; std::uint32_t chamfered_edges = 0; };`
  - `Result<ChamferResult> chamfer_edges(const Brep&, std::span<const EdgeId>, const ChamferOptions& = {});`

- [ ] **Step 1: Write the failing single-edge test**

```cpp
#include <doctest/doctest.h>

#include <array>
#include <cmath>
#include <cstdint>
#include <limits>
#include <span>
#include <string>
#include <vector>

#include <geom/brep/boolean_solid.hpp>
#include <geom/brep/chamfer.hpp>
#include <geom/brep/mass_properties.hpp>
#include <geom/brep/primitives.hpp>
#include <geom/brep/validate.hpp>
#include <geom/core/error.hpp>

using namespace geom;

namespace {

Brep make_box_solid(Point3d lo, Point3d hi) {
  const Result<Brep> box = make_box(lo, hi, 1e-9);
  REQUIRE(box.has_value());
  return box.value();
}

double volume_of(const Brep& brep) {
  const Result<MassProperties> mass = mass_properties(brep);
  REQUIRE(mass.has_value());
  return mass.value().volume;
}

// Finds the straight edge whose endpoints are the two given corners.
EdgeId edge_between(const Brep& brep, Point3d a, Point3d b) {
  for (const EdgeId id : brep.edge_ids()) {
    const Point3d p0 = brep.vertex(brep.edge(id).v0).point;
    const Point3d p1 = brep.vertex(brep.edge(id).v1).point;
    const auto same = [](Point3d p, Point3d q) {
      return std::abs(p[0] - q[0]) + std::abs(p[1] - q[1]) +
                 std::abs(p[2] - q[2]) <
             1e-12;
    };
    if ((same(p0, a) && same(p1, b)) || (same(p0, b) && same(p1, a))) return id;
  }
  REQUIRE_MESSAGE(false, "edge not found");
  return {};
}

}  // namespace

TEST_CASE("chamfer one box edge removes a triangular prism") {
  const Brep box = make_box_solid(Point3d{0.0, 0.0, 0.0}, Point3d{2.0, 3.0, 4.0});
  const EdgeId edge =
      edge_between(box, Point3d{0.0, 0.0, 0.0}, Point3d{0.0, 0.0, 4.0});
  const std::array<EdgeId, 1> edges{edge};

  ChamferOptions options;
  options.distance = 0.25;
  const Result<ChamferResult> result = chamfer_edges(box, edges, options);
  REQUIRE(result.has_value());
  CHECK(result.value().chamfered_edges == 1u);
  CHECK(result.value().brep.solid_count() == 1u);
  CHECK(result.value().brep.shell_count() == 1u);
  CHECK(validate_topology(result.value().brep).ok());
  // Removed cross-section is d^2/2 = 0.03125 over the 4-long edge.
  CHECK(std::abs(volume_of(result.value().brep) - 23.875) <= 1e-9);
  // The input document is untouched.
  CHECK(box.face_count() == 6u);
  CHECK(std::abs(volume_of(box) - 24.0) <= 1e-12);
}

TEST_CASE("chamfer rejects a non-positive distance") {
  const Brep box = make_box_solid(Point3d{0.0, 0.0, 0.0}, Point3d{1.0, 1.0, 1.0});
  const std::array<EdgeId, 1> edges{box.edge_ids().front()};
  ChamferOptions options;
  options.distance = 0.0;
  const Result<ChamferResult> result = chamfer_edges(box, edges, options);
  REQUIRE_FALSE(result.has_value());
  CHECK(result.error().code == Errc::degenerate);
}
```

- [ ] **Step 2: Run the test to verify it fails**

Run: `cmake --build build --config Debug --target geom_tests && .\build\Debug\geom_tests.exe --test-case="*chamfer*"`
Expected: compile error `Cannot open include file: 'geom/brep/chamfer.hpp'` (the header and source list entry do not exist yet).

- [ ] **Step 3: Implement the shared helpers and `chamfer_edges`**

`blend_planar_convex_solid` checks, in order: exactly one solid and one shell (`Errc::unsupported`), a closed shell (`Errc::invalid_topology`), every face surface holds `Plane3` (`Errc::unsupported`), and every vertex lies on the inner side of every outward face plane within `max(64 * tolerance, 1e-12 * scale)` (`Errc::unsupported`, "solid must be convex"). It stores unit outward normals (`face.reversed ? plane.normal * -1.0 : plane.normal`; `Vec` has no unary minus) and normalized offsets.

`blend_inward_direction` builds the in-face perpendicular `cross(outward, edge_dir)`, orients it toward the face centroid, and rejects a zero result. `blend_face_width` projects the face loop vertices into the edge frame `(direction, cross(outward, direction))`, intersects the ray with each boundary segment, skips the selected edge's own endpoints, and returns the smallest positive distance.

`chamfer_edges` validates the options (`distance` finite and positive, `tolerance` finite and positive), builds the planar-convex view, requires a non-empty distinct edge list, and for each edge computes the two setback points `A = mid + d * inward0`, `B = mid + d * inward1`. Reject coplanar faces when `norm2(inward0 + inward1) <= 1e-12`. The cut plane is `Plane3{normalize(inward0 + inward1), dot(unit, A)}`; the edge lies on its negative side, so keep `above`. Check `d < blend_face_width(...)` for each adjacent face first. Apply the cuts sequentially to a working copy (`current = std::move(split->above)`, then `current.solid_ids().front()`), require one solid after each cut, and run `validate_topology` before returning.

- [ ] **Step 4: Run Debug and Release tests**

Run:
`cmake --build build --config Debug --target geom_tests && .\build\Debug\geom_tests.exe --test-case="*chamfer*"`
then
`cmake --build build --config Release --target geom_tests && .\build\Release\geom_tests.exe --test-case="*chamfer*"`
Expected: both report 2 test cases and 11 assertions passing, 0 failures, 0 warnings.

- [ ] **Step 5: Commit**

```powershell
git add include/geom/brep/detail/blend_common.hpp include/geom/brep/chamfer.hpp tests/unit/chamfer_test.cpp tests/CMakeLists.txt
git commit -m "feat(brep): chamfer a straight edge of a convex planar solid"
```

---

### Task 2: Multi-edge chamfer and vertex facets

**Files:**
- Modify: `include/geom/brep/chamfer.hpp`
- Test: `tests/unit/chamfer_test.cpp`

**Interfaces:**
- Consumes: Task 1's `chamfer_edges`, `detail::BlendPlanarSolid`, `detail::blend_inward_direction`, `detail::blend_face_width`.
- Produces: the same `chamfer_edges` signature with all planes computed from the original solid before any cut, duplicate `EdgeId`s counted once, and `chamfered_edges` equal to the number of distinct edges applied.

- [ ] **Step 1: Add the failing multi-edge tests**

```cpp
TEST_CASE("chamfer two edges sharing a vertex subtracts their overlap once") {
  const Brep box = make_box_solid(Point3d{0.0, 0.0, 0.0}, Point3d{2.0, 2.0, 2.0});
  const std::array<EdgeId, 2> edges{
      edge_between(box, Point3d{0.0, 0.0, 0.0}, Point3d{0.0, 0.0, 2.0}),
      edge_between(box, Point3d{0.0, 0.0, 0.0}, Point3d{0.0, 2.0, 0.0})};
  ChamferOptions options;
  options.distance = 0.25;
  const Result<ChamferResult> result = chamfer_edges(box, edges, options);
  REQUIRE(result.has_value());
  CHECK(result.value().chamfered_edges == 2u);
  CHECK(validate_topology(result.value().brep).ok());
  // Two prisms of 0.0625 each, their d^3/3 corner overlap removed once.
  CHECK(std::abs(volume_of(result.value().brep) - 7.880208333333333) <= 1e-9);
}

TEST_CASE("chamfer all twelve cube edges leaves the analytic volume") {
  const Brep box = make_box_solid(Point3d{0.0, 0.0, 0.0}, Point3d{2.0, 2.0, 2.0});
  std::vector<EdgeId> edges = box.edge_ids();
  REQUIRE(edges.size() == 12u);
  ChamferOptions options;
  options.distance = 0.25;
  const Result<ChamferResult> result = chamfer_edges(box, edges, options);
  REQUIRE(result.has_value());
  CHECK(result.value().chamfered_edges == 12u);
  CHECK(validate_topology(result.value().brep).ok());
  // Inclusion-exclusion over the twelve prisms: a^3 - 6*a*d^2 + 6*d^3.
  CHECK(std::abs(volume_of(result.value().brep) - 7.34375) <= 1e-9);
}

TEST_CASE("chamfer counts a repeated edge once") {
  const Brep box = make_box_solid(Point3d{0.0, 0.0, 0.0}, Point3d{1.0, 1.0, 1.0});
  const EdgeId edge = edge_between(box, Point3d{0.0, 0.0, 0.0}, Point3d{0.0, 0.0, 1.0});
  const std::array<EdgeId, 2> edges{edge, edge};
  ChamferOptions options;
  options.distance = 0.1;
  const Result<ChamferResult> result = chamfer_edges(box, edges, options);
  REQUIRE(result.has_value());
  CHECK(result.value().chamfered_edges == 1u);
}
```

- [ ] **Step 2: Run to verify failure**

Run: `.\build\Debug\geom_tests.exe --test-case="*chamfer*"`
Expected: the all-edges case fails on `chamfered_edges`/volume; a repeated edge either double-cuts or reports 2.

- [ ] **Step 3: Implement the multi-edge pass**

Deduplicate handles by `EdgeId` (stable first-occurrence order), compute every cut plane and both width checks from the original `Brep` into a vector, then apply the planes in that order. Do not recompute faces or edges from the shrinking solid. After each cut, require `current.solid_count() == 1u`; otherwise return `Errc::invalid_topology` ("chamfer_edges: cut split the solid").

- [ ] **Step 4: Run Debug and Release tests**

Run both configurations with `--test-case="*chamfer*"`.
Expected: 5 cases, 22 assertions, 0 failures.

- [ ] **Step 5: Commit**

```powershell
git add include/geom/brep/chamfer.hpp tests/unit/chamfer_test.cpp
git commit -m "feat(brep): chamfer a chain of edges and their shared vertices"
```

---

### Task 3: Chamfer diagnostics and width guard

**Files:**
- Modify: `include/geom/brep/chamfer.hpp`
- Modify: `include/geom/brep/detail/blend_common.hpp`
- Test: `tests/unit/chamfer_test.cpp`

**Interfaces:**
- Consumes: Task 2's implementation.
- Produces: stable `Errc`/message ordering — options, shape, planarity, convexity, edge list, per-edge topology, width, cut; each negative test asserts both code and a message prefix.

- [ ] **Step 1: Add the failing diagnostic tests**

```cpp
TEST_CASE("chamfer rejects a curved operand") {
  const Result<Brep> cylinder = make_cylinder(Point3d{0.0, 0.0, 0.0},
                                              Vec3d{0.0, 0.0, 1.0}, 1.0, 2.0, 1e-9);
  REQUIRE(cylinder.has_value());
  const std::array<EdgeId, 1> edges{cylinder.value().edge_ids().front()};
  ChamferOptions options;
  options.distance = 0.1;
  const Result<ChamferResult> result = chamfer_edges(cylinder.value(), edges, options);
  REQUIRE_FALSE(result.has_value());
  CHECK(result.error().code == Errc::unsupported);
  CHECK(result.error().message.find("all faces must be planar") != std::string::npos);
}

TEST_CASE("chamfer rejects a non-convex solid") {
  const Brep a = make_box_solid(Point3d{0.0, 0.0, 0.0}, Point3d{2.0, 1.0, 1.0});
  const Brep b = make_box_solid(Point3d{0.0, 0.5, 0.0}, Point3d{1.0, 1.5, 1.0});
  const Result<SolidBooleanResult> formed =
      boolean_solid(a, b, SolidBooleanOp::union_);
  REQUIRE(formed.has_value());
  const std::array<EdgeId, 1> edges{formed.value().brep.edge_ids().front()};
  ChamferOptions options;
  options.distance = 0.1;
  const Result<ChamferResult> result = chamfer_edges(formed.value().brep, edges, options);
  REQUIRE_FALSE(result.has_value());
  CHECK(result.error().code == Errc::unsupported);
  CHECK(result.error().message.find("solid must be convex") != std::string::npos);
}

TEST_CASE("chamfer rejects a distance at the face width") {
  const Brep box = make_box_solid(Point3d{0.0, 0.0, 0.0}, Point3d{2.0, 2.0, 2.0});
  const std::array<EdgeId, 1> edges{
      edge_between(box, Point3d{0.0, 0.0, 0.0}, Point3d{0.0, 0.0, 2.0})};
  ChamferOptions options;
  options.distance = 2.0;
  const Result<ChamferResult> result = chamfer_edges(box, edges, options);
  REQUIRE_FALSE(result.has_value());
  CHECK(result.error().code == Errc::unsupported);
  CHECK(result.error().message.find("exceeds the adjacent face width") != std::string::npos);
}

TEST_CASE("chamfer rejects an empty edge list and a dead edge") {
  const Brep box = make_box_solid(Point3d{0.0, 0.0, 0.0}, Point3d{1.0, 1.0, 1.0});
  ChamferOptions options;
  options.distance = 0.1;
  const std::array<EdgeId, 0> none{};
  const Result<ChamferResult> empty = chamfer_edges(box, none, options);
  REQUIRE_FALSE(empty.has_value());
  CHECK(empty.error().code == Errc::unsupported);

  const std::array<EdgeId, 1> dead{EdgeId{}};
  const Result<ChamferResult> bad = chamfer_edges(box, dead, options);
  REQUIRE_FALSE(bad.has_value());
  CHECK(bad.error().code == Errc::index_out_of_range);
}

TEST_CASE("chamfer rejects an option error before touching the solid") {
  const Brep box = make_box_solid(Point3d{0.0, 0.0, 0.0}, Point3d{1.0, 1.0, 1.0});
  const std::array<EdgeId, 1> edges{box.edge_ids().front()};
  ChamferOptions options;
  options.distance = std::numeric_limits<double>::quiet_NaN();
  const Result<ChamferResult> non_finite = chamfer_edges(box, edges, options);
  REQUIRE_FALSE(non_finite.has_value());
  CHECK(non_finite.error().code == Errc::non_finite);

  options.distance = 0.1;
  options.tolerance = 0.0;
  const Result<ChamferResult> bad_tolerance = chamfer_edges(box, edges, options);
  REQUIRE_FALSE(bad_tolerance.has_value());
  CHECK(bad_tolerance.error().code == Errc::tolerance_violation);
}
```

- [ ] **Step 2: Run to verify failure**

Run: `.\build\Debug\geom_tests.exe --test-case="*chamfer*"`
Expected: at least the non-convex and width cases fail or report a different code/message.

- [ ] **Step 3: Implement the checks in the documented order**

Options first, then `blend_planar_convex_solid`, then the edge list (non-empty, alive, `LineCurve3`, two coedges on distinct faces), then the width guard. Use `std::string_view` caller prefixes so `blend_common.hpp` can serve shell. Do not fall back to a smaller distance and do not return the input on failure.

- [ ] **Step 4: Run Debug and Release tests**

Run both configurations with `--test-case="*chamfer*"`.
Expected: 10 cases, 39 assertions, 0 failures.

- [ ] **Step 5: Commit**

```powershell
git add include/geom/brep/chamfer.hpp include/geom/brep/detail/blend_common.hpp tests/unit/chamfer_test.cpp
git commit -m "fix(brep): diagnose unsupported chamfer operands and oversize distances"
```

---

### Task 4: Closed shell

**Files:**
- Create: `include/geom/brep/shell.hpp`
- Test: `tests/unit/shell_test.cpp`
- Modify: `tests/CMakeLists.txt`

**Interfaces:**
- Consumes: `detail::BlendPlanarSolid`, `split_solid_by_plane`, `boolean_solid(const Brep&, const Brep&, SolidBooleanOp, const SolidBooleanOptions&)`, `validate_topology`.
- Produces:
  - `struct ShellOptions { double thickness = 0.0; double tolerance = 1e-9; };`
  - `struct ShellResult { Brep brep; std::uint32_t cavity_shells = 0; };`
  - `Result<ShellResult> shell_solid(const Brep&, std::span<const FaceId> open_faces, const ShellOptions& = {});`

- [ ] **Step 1: Write the failing closed-shell test**

```cpp
#include <doctest/doctest.h>

#include <array>
#include <cmath>
#include <cstdint>
#include <limits>
#include <span>
#include <string>
#include <variant>

#include <geom/brep/mass_properties.hpp>
#include <geom/brep/primitives.hpp>
#include <geom/brep/shell.hpp>
#include <geom/brep/validate.hpp>
#include <geom/core/error.hpp>

using namespace geom;

namespace {

Brep box_solid(Point3d lo, Point3d hi) {
  const Result<Brep> box = make_box(lo, hi, 1e-9);
  REQUIRE(box.has_value());
  return box.value();
}

double volume_of(const Brep& brep) {
  const Result<MassProperties> mass = mass_properties(brep);
  REQUIRE(mass.has_value());
  return mass.value().volume;
}

}  // namespace

TEST_CASE("closed shell leaves a uniform cavity") {
  const Brep box = box_solid(Point3d{0.0, 0.0, 0.0}, Point3d{2.0, 3.0, 4.0});
  ShellOptions options;
  options.thickness = 0.25;
  const std::array<FaceId, 0> open{};
  const Result<ShellResult> result = shell_solid(box, open, options);
  REQUIRE(result.has_value());
  CHECK(result.value().cavity_shells == 1u);
  CHECK(result.value().brep.solid_count() == 1u);
  const auto shells = result.value().brep.solid_shells(result.value().brep.solid_ids().front());
  CHECK(shells.size() == 2u);
  CHECK(validate_topology(result.value().brep).ok());
  // 24 - (1.5 * 2.5 * 3.5).
  CHECK(std::abs(volume_of(result.value().brep) - 10.875) <= 1e-9);
}

TEST_CASE("shell rejects a thickness that leaves no interior") {
  const Brep box = box_solid(Point3d{0.0, 0.0, 0.0}, Point3d{2.0, 3.0, 4.0});
  ShellOptions options;
  options.thickness = 1.0;
  const std::array<FaceId, 0> open{};
  const Result<ShellResult> result = shell_solid(box, open, options);
  REQUIRE_FALSE(result.has_value());
  CHECK(result.error().code == Errc::degenerate);
  CHECK(result.error().message.find("leaves no interior") != std::string::npos);
}
```

- [ ] **Step 2: Run to verify failure**

Run: `cmake --build build --config Debug --target geom_tests && .\build\Debug\geom_tests.exe --test-case="*shell*"`
Expected: compile error for the missing `shell.hpp`.

- [ ] **Step 3: Implement `shell_solid`**

Build `detail::blend_planar_convex_solid`. Validate options. Validate and deduplicate `open_faces` (dead or out-of-range is `Errc::index_out_of_range`; an empty closed set is `Errc::unsupported`). Build the cavity as the intersection of one half-space per face, seeded by a box that contains the input with a numeric margin `max(2 * thickness, 1e-6 * max(1, extent), 1000 * tolerance)` on every side: a closed face cuts with `Plane3{outward[i], offset[i] - thickness}` and keeps `below`; an open face cuts with `Plane3{outward[i], offset[i] + thickness}` and keeps `below`, so the cavity protrudes through the opening. Shifting an open plane outward while starting from the input solid would be an identity cut and leave a coincident cap, so the seed box is required. Cut axis-aligned faces before oblique ones (stable) so an oblique plane never meets a face an earlier oblique cut already trimmed. Require exactly one solid after each cut. If a plane is at or below every vertex, or a cut leaves no solid, return `Errc::degenerate` ("shell_solid: wall thickness leaves no interior"). If a split fails, test whether the shifted half-spaces still contain a point more than `tolerance` inside every plane by averaging the feasible triple-intersection vertices; if none, return `Errc::degenerate` ("shell_solid: wall thickness leaves no measurable interior at this tolerance"), otherwise propagate the splitter's error code with the `shell_solid:` prefix ("shell_solid: cavity intersection could not be assembled"). Then `boolean_solid(brep, inner, SolidBooleanOp::difference, {options.tolerance})`; a non-finite result is `Errc::non_finite` ("shell_solid: cavity difference failed a numeric check") and any other failure is the "no measurable interior" degenerate case. Require one solid, `validate_topology`, and set `cavity_shells = result_shells - 1` (closed open-list) or `0` when any face is open.

- [ ] **Step 4: Run Debug and Release tests**

Run both configurations with `--test-case="*shell*"`.
Expected: eleven added cases (closed, no interior, open top, opposite open, wall below resolution, cavity thinner than the tolerance, all faces open, oblique prism, classifier witness, all-oblique cavity, sub-ulp seed); the `*shell*` filter also matches unrelated pre-existing tests, so check the new case names rather than the total.

- [ ] **Step 5: Commit**

```powershell
git add include/geom/brep/shell.hpp tests/unit/shell_test.cpp tests/CMakeLists.txt
git commit -m "feat(brep): shell a convex planar solid with a uniform wall"
```

---

### Task 5: Open-face shell and shell diagnostics

**Files:**
- Modify: `include/geom/brep/shell.hpp`
- Test: `tests/unit/shell_test.cpp`

**Interfaces:**
- Consumes: Task 4's `shell_solid`.
- Produces: open-face semantics — every requested face is omitted from the inner half-space intersection, the pocket opens through it, and `cavity_shells` is 0; duplicate open faces are applied once.

- [ ] **Step 1: Add diagnostic and duplicate-face tests**

```cpp
TEST_CASE("shell rejects bad options and unsupported operands") {
  const Brep box = box_solid(Point3d{0.0, 0.0, 0.0}, Point3d{2.0, 2.0, 2.0});
  const std::array<FaceId, 0> open{};
  ShellOptions options;
  options.thickness = 0.0;
  CHECK(shell_solid(box, open, options).error().code == Errc::degenerate);

  options.thickness = 0.25;
  options.tolerance = std::numeric_limits<double>::infinity();
  CHECK(shell_solid(box, open, options).error().code == Errc::non_finite);

  const std::array<FaceId, 1> dead{FaceId{}};
  options.tolerance = 1e-9;
  CHECK(shell_solid(box, dead, options).error().code == Errc::index_out_of_range);
}

TEST_CASE("a duplicate open face is applied once") {
  const Brep box = box_solid(Point3d{0.0, 0.0, 0.0}, Point3d{2.0, 3.0, 4.0});
  FaceId top{};
  for (const FaceId face : box.face_ids()) {
    const auto& plane = std::get<Plane3>(box.geometry().surfaces[box.face(face).surface]);
    const Vec3d normal = box.face(face).reversed ? plane.normal * -1.0 : plane.normal;
    if (normal[2] > 0.5) top = face;
  }
  REQUIRE_FALSE(top.is_null());
  ShellOptions options;
  options.thickness = 0.25;
  const std::array<FaceId, 2> open{top, top};
  const Result<ShellResult> result = shell_solid(box, open, options);
  REQUIRE(result.has_value());
  CHECK(result.value().cavity_shells == 0u);
  CHECK(std::abs(volume_of(result.value().brep) - 9.9375) <= 1e-9);
}
```

- [ ] **Step 2: Run the new cases**

Run: `.\build\Debug\geom_tests.exe --test-case="*shell*"`.
Expected: both cases pass; the open-face and all-open behavior already landed in Task 4, so this step is test-only. If the duplicate case fails, fix deduplication in `shell_solid` before continuing.

- [ ] **Step 3: Keep the guards covered**

The Task 4 implementation already deduplicates open faces in first-occurrence order, rejects a dead or out-of-range `FaceId` with `Errc::index_out_of_range`, rejects an empty closed set with `Errc::unsupported` ("shell_solid: at least one face must stay closed"), and maps unrepresentable thin walls to `Errc::degenerate`. A `FaceId` from a *different* `Brep` cannot be detected: handles are only `(index, generation)`. Do not change those contracts here; only adjust the code if Step 2 finds a real defect.

- [ ] **Step 4: Run Debug and Release tests**

Run both configurations with `--test-case="*shell*"`.
Expected: the two new cases pass in both configurations alongside the eight Task 4 cases.

- [ ] **Step 5: Commit**

```powershell
git add tests/unit/shell_test.cpp
git commit -m "test(brep): pin shell diagnostics and duplicate open faces"
```

---

### Task 6: Watertightness, determinism and input immutability

**Files:**
- Create: `tests/unit/blend_qa_test.cpp`
- Modify: `tests/CMakeLists.txt`

**Interfaces:**
- Consumes: `chamfer_edges`, `shell_solid`, `tessellate`, `TriMesh`, `validate_topology`, `mass_properties`.
- Produces: QA tests that every output is watertight (no triangle edge used once), deterministic (two runs produce identical positions and triangles), and that the inputs are unchanged.

- [ ] **Step 1: Write the failing QA tests**

```cpp
#include <doctest/doctest.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <map>
#include <span>
#include <utility>
#include <vector>

#include <geom/brep/chamfer.hpp>
#include <geom/brep/mass_properties.hpp>
#include <geom/brep/primitives.hpp>
#include <geom/brep/shell.hpp>
#include <geom/brep/tessellate.hpp>
#include <geom/brep/validate.hpp>
#include <geom/core/error.hpp>

using namespace geom;

namespace {

Brep box_solid(Point3d lo, Point3d hi) {
  const Result<Brep> box = make_box(lo, hi, 1e-9);
  REQUIRE(box.has_value());
  return box.value();
}

// Returns the number of triangle edges used exactly once.
std::size_t boundary_edges(const TriMesh& mesh) {
  std::map<std::pair<std::uint32_t, std::uint32_t>, int> uses;
  for (const auto& tri : mesh.triangles) {
    for (int i = 0; i < 3; ++i) {
      std::uint32_t a = tri[static_cast<std::size_t>(i)];
      std::uint32_t b = tri[static_cast<std::size_t>((i + 1) % 3)];
      if (a > b) std::swap(a, b);
      ++uses[{a, b}];
    }
  }
  std::size_t count = 0;
  for (const auto& [edge, n] : uses) {
    (void)edge;
    if (n == 1) ++count;
  }
  return count;
}

}  // namespace

TEST_CASE("chamfer and shell outputs tessellate watertight") {
  const Brep box = box_solid(Point3d{0.0, 0.0, 0.0}, Point3d{2.0, 3.0, 4.0});
  ChamferOptions chamfer;
  chamfer.distance = 0.25;
  std::vector<EdgeId> edges = box.edge_ids();
  const Result<ChamferResult> beveled = chamfer_edges(box, edges, chamfer);
  REQUIRE(beveled.has_value());
  const Result<TriMesh> bevel_mesh = tessellate(beveled.value().brep, TessellationOptions{});
  REQUIRE(bevel_mesh.has_value());
  CHECK_FALSE(bevel_mesh.value().triangles.empty());
  CHECK(boundary_edges(bevel_mesh.value()) == 0u);

  ShellOptions shell;
  shell.thickness = 0.25;
  const std::array<FaceId, 0> open{};
  const Result<ShellResult> hollow = shell_solid(box, open, shell);
  REQUIRE(hollow.has_value());
  const Result<TriMesh> shell_mesh = tessellate(hollow.value().brep, TessellationOptions{});
  REQUIRE(shell_mesh.has_value());
  CHECK(boundary_edges(shell_mesh.value()) == 0u);
}

TEST_CASE("blend outputs are deterministic and leave inputs unchanged") {
  const Brep box = box_solid(Point3d{0.0, 0.0, 0.0}, Point3d{2.0, 2.0, 2.0});
  ChamferOptions chamfer;
  chamfer.distance = 0.25;
  const std::vector<EdgeId> edges = box.edge_ids();
  const Result<ChamferResult> first = chamfer_edges(box, edges, chamfer);
  const Result<ChamferResult> second = chamfer_edges(box, edges, chamfer);
  REQUIRE(first.has_value());
  REQUIRE(second.has_value());
  CHECK(first.value().brep.vertex_count() == second.value().brep.vertex_count());
  CHECK(first.value().brep.edge_count() == second.value().brep.edge_count());
  CHECK(first.value().brep.face_count() == second.value().brep.face_count());
  const Result<MassProperties> a = mass_properties(first.value().brep);
  const Result<MassProperties> b = mass_properties(second.value().brep);
  REQUIRE(a.has_value());
  REQUIRE(b.has_value());
  CHECK(a.value().volume == b.value().volume);

  const Result<TriMesh> mesh_a = tessellate(first.value().brep, TessellationOptions{});
  const Result<TriMesh> mesh_b = tessellate(second.value().brep, TessellationOptions{});
  REQUIRE(mesh_a.has_value());
  REQUIRE(mesh_b.has_value());
  CHECK(mesh_a.value().positions.size() == mesh_b.value().positions.size());
  CHECK(mesh_a.value().triangles == mesh_b.value().triangles);

  CHECK(box.face_count() == 6u);
  const Result<MassProperties> original = mass_properties(box);
  REQUIRE(original.has_value());
  CHECK(std::abs(original.value().volume - 8.0) <= 1e-12);
}
```

- [ ] **Step 2: Run to verify failure**

Run: `cmake --build build --config Debug --target geom_tests && .\build\Debug\geom_tests.exe --test-case="*watertight*,*deterministic*"`
Expected: compile error until the file is registered, then a real boundary or determinism failure if the implementation has one.

- [ ] **Step 3: Fix implementation if a QA check fails**

A boundary edge means the cut result has a T-junction or a missing cap: do not weaken the test; return to Task 1/4 and fix the plane choice or validate the result before returning. Determinism failures mean unordered iteration; use stable first-occurrence order everywhere.

- [ ] **Step 4: Run Debug and Release tests and the full suite**

Run both configurations with `--test-case="*watertight*,*deterministic*"`, then the full suite.
Expected: 2 QA cases pass; the full Debug and Release suites report 0 failures.

- [ ] **Step 5: Commit**

```powershell
git add tests/unit/blend_qa_test.cpp tests/CMakeLists.txt
git commit -m "test(brep): pin blend watertightness, determinism and input immutability"
```

---

### Task 7: Benchmark, ADR and milestone evidence

**Files:**
- Create: `bench/bench_m4c1.cpp`
- Modify: `CMakeLists.txt`
- Create: `docs/benchmarks/m4c1-blends.md`
- Create: `docs/adr/0041-planar-blends-via-half-space-cuts.md`
- Modify: `docs/architecture/geometry-library-architecture.md`
- Modify: this plan

**Interfaces:**
- Consumes: `chamfer_edges`, `shell_solid`, `make_box`, `mass_properties`.
- Produces: reproducible timing evidence and the documented scope boundary.

- [ ] **Step 1: Write the benchmark**

`bench/bench_m4c1.cpp` builds a `100 x 100 x 100` box with tolerance `1e-9`, uses `brep.edge_ids()` (12 edges) and `ChamferOptions{1.0, 1e-9}`, times five runs with `std::chrono::steady_clock`, prints min/median/max milliseconds and the resulting volume. It then shells the same box with `ShellOptions{1.0, 1e-9}` and no open faces, printing the same statistics and `volume == 58808`. A failed operation prints the `Errc` and exits `2`. The volume checks are checksums: if `chamfered` volume is not `999406` within `1e-6` relative, or shell volume is not `58808` within `1e-6`, exit `3`.

- [ ] **Step 2: Register and run it**

Add under `if(GEOM_BUILD_BENCH)`:

```cmake
  add_executable(geom_bench_m4c1 bench/bench_m4c1.cpp)
  target_link_libraries(geom_bench_m4c1 PRIVATE geom_core)
```

Run: `cmake -S . -B build-bench -DGEOM_BUILD_BENCH=ON -DGEOM_WERROR=ON -DCMAKE_BUILD_TYPE=Release && cmake --build build-bench --config Release --target geom_bench_m4c1 && .\build-bench\Release\geom_bench_m4c1.exe`
Expected: exit 0, volumes `999406` and `58808`, and three timing lines.

- [ ] **Step 3: Record the evidence**

Write `docs/benchmarks/m4c1-blends.md` with: compiler and CPU, build type, the five-run min/median/max for both operations, the two volume checksums, the test counts from Task 6, and the explicit scope ("convex, planar, single-shell; fillet deferred to M4c-2"). Add ADR-0041 with the decision (half-space cuts, convex restriction, no faceted fillet) and the alternatives rejected (polyhedral fillet approximation, boolean tool solids). Add a §8.9 scope paragraph and a §16.4 M4c-1 entry with the same counts.

- [ ] **Step 4: Run the complete suite in both configurations**

Debug and Release, full `geom_tests`, plus `ctest`.
Expected: 0 failures; record the exact case/assertion counts in the plan's Results table.

- [ ] **Step 5: Commit**

```powershell
git add bench/bench_m4c1.cpp CMakeLists.txt docs/benchmarks/m4c1-blends.md docs/adr/0041-planar-blends-via-half-space-cuts.md docs/architecture/geometry-library-architecture.md docs/superpowers/plans/2026-11-25-m4c1-planar-chamfer-shell.md
git commit -m "docs(m4c1): record planar chamfer and shell evidence"
```

---

## M4c-1 Exit Criteria

- `chamfer_edges` accepts a convex single-shell planar solid and a non-empty list of straight two-coedge edges, returns a valid closed solid, and matches the analytic volume for one edge, two adjacent edges, and all twelve cube edges.
- `shell_solid` returns a closed cavity for an empty open list and one open pocket for one or more open faces, with the analytic volumes for the box fixtures.
- Every unsupported input (curved face, non-convex solid, dead face/edge, empty list, coplanar edge, oversized distance, vanishing wall) returns a stable `Errc` and message and leaves the input unchanged.
- Tessellation of every output has zero boundary edges; two runs are bit-identical; `validate_topology` is clean.
- `docs/benchmarks/m4c1-blends.md`, ADR-0041 and the architecture §8.9/§16.4 text carry the measured numbers and the explicit convex-planar limit.
- Full Debug and Release suites pass; CI is green before merge.

## Results

| Item | Value |
| --- | --- |
| Branch / base | `sdd/m4c1`, base `9e39182` |
| Chamfer one edge volume | TBD |
| Chamfer all edges volume | TBD |
| Shell closed volume | TBD |
| Shell open-top volume | TBD |
| Bench chamfer all edges (min/median) | TBD |
| Bench closed shell (min/median) | TBD |
| Debug cases / assertions | TBD |
| Release cases / assertions | TBD |
| Merged commit | TBD |

## Self-Review

- **Spec coverage:** The plan implements the two blend operations that reduce exactly to half-space cuts. Fillet, variable radius, chamfer angle, non-convex and curved operands, offset and sweep are explicitly deferred to M4c-2, so §8.9 is only partially covered by design; the deferral is written into ADR-0041 and §16.4 rather than implied.
- **Step scan:** Each step names the test or the exact function to write; no step says "handle edge cases" or "add validation" without the check and its error.
- **Type consistency:** `ChamferOptions`/`ChamferResult`/`chamfer_edges` and `ShellOptions`/`ShellResult`/`shell_solid` are used with the same names in tests, QA and bench; the shared detail frame is declared once in Task 1 and consumed by Tasks 2–5.
- **Review Focus:** face-width overshoot (Task 3), vertex facet (Task 2), open pocket (Task 5), vanishing wall (Task 4), curved/non-convex rejection (Tasks 3 and 5) each have a named test.
- **Proportion:** The plan lists signatures and tests; code bodies are left to the implementer except the test assertions and the exact volume formulas.
