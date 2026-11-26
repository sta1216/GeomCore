# M4c-2 Half-Space Shell Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use subagent-driven-development (recommended) or executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Let `shell_solid` open arbitrary convex planar solids through selected faces, including oblique faces, without relying on caller face order.

**Architecture:** The cavity is the intersection of half-spaces: each closed shell face moves inward by the wall thickness, each open face moves outward, and a seed box closes the unbounded directions. The cavity polyhedron is built directly from plane-triple vertices and assembled with `assemble_solid`, then subtracted from the input; no sequential plane splitting is involved, so the result cannot depend on face order. The existing inset-box fast path stays for closed axis-aligned boxes.

**Tech Stack:** C++20, CMake 4.4.4, MSVC 19.44 `/W4 /WX /utf-8`, doctest 2.4.11, existing `geom/brep` planar B-rep and `boolean_solid` difference.

**Spec:** `docs/architecture/geometry-library-architecture.md` §8.9 and §4.6; the accepted M4c-1 contract and its deviations are in `docs/superpowers/plans/2026-11-25-m4c1-planar-chamfer-shell.md` and ADR-0041.

## Scope

M4c-2 completes the shell operation:

- any single-shell, closed, convex, all-planar solid (including oblique faces);
- any subset of faces may be open, including two opposite faces;
- wall thickness and tolerance are explicit; no global epsilon.

Deferred to later milestones, and not claimed here: fillet (constant or variable radius), edge chains, offset/sweep, non-convex or multi-shell operands, curved faces, and the general `split_solid_by_plane` face-order hardening. The reason fillet is not folded in: an exact rolling-ball fillet needs trimmed cylindrical and spherical patches with pcurves and vertex blends, which is a separate subsystem, and a faceted approximation is forbidden by §8.9.

## Global Constraints

- C++20; public headers compile standalone with `/W4 /WX /utf-8`.
- No new third-party dependencies.
- No global epsilon. All tolerances are the caller's `ShellOptions::tolerance` scaled by the input extent, and each use is named in the source.
- `shell_solid` either returns a validated closed solid or an `Error`; it never returns a half-built solid and never mutates the input.
- Determinism: permuting `Brep` face/vertex order must not change the geometry of the result.
- Validation before modification: options, solid/shell count, planarity, convexity, open-face ids, and face width run before any builder is called.
- TDD for every task: failing test first, then minimal implementation, then Debug and Release `*shell*` plus the relevant lane.
- One commit per task; do not amend a reviewed task without a new review.

## Review Focus

1. A closed box whose faces are presented in any of the 720 orders still shells to the same volume and a valid solid.
2. An open box (one or two opposite open faces) produces a pocket that actually breaks the surface, with no coincident cap left behind.
3. A convex oblique solid (the right tetrahedron in the tests) shells with the analytically expected inner volume.
4. A thickness at or beyond the local half-thickness returns `degenerate`, not a leaked assembler or boolean error.
5. The cavity is wound outward (`assemble_solid` sees positive volume); a normal or ring-order sign error must fail a test, not silently invert the cavity.

---

### Task 1: Half-space vertices

**Files:**
- Create: `include/geom/brep/detail/half_space.hpp`
- Test: `tests/unit/half_space_test.cpp`
- Modify: `tests/CMakeLists.txt` (append `unit/half_space_test.cpp`)

**Interfaces:**
- Produces: `struct HalfSpacePlane { Vec3d normal; double offset; };` with inside meaning `dot(normal, x) <= offset`, and `Result<std::vector<Point3d>> half_space_vertices(std::span<const HalfSpacePlane> planes, double tolerance);`. `normal` may be non-unit on input; the helper normalizes it and divides `offset` by the same length. A plane with a zero or non-finite normal is rejected. Empty intersection is an empty vector, not an error.

- [ ] **Step 1: Write the failing tests**

Add `tests/unit/half_space_test.cpp`:

```cpp
#include <doctest/doctest.h>

#include <array>
#include <cmath>
#include <span>
#include <string_view>
#include <vector>

#include <geom/brep/detail/half_space.hpp>

using namespace geom;

namespace {

std::vector<detail::HalfSpacePlane> unit_box_planes() {
  return {
      {{1.0, 0.0, 0.0}, 2.0},
      {{-1.0, 0.0, 0.0}, 0.0},
      {{0.0, 1.0, 0.0}, 3.0},
      {{0.0, -1.0, 0.0}, 0.0},
      {{0.0, 0.0, 1.0}, 4.0},
      {{0.0, 0.0, -1.0}, 0.0},
  };
}

}  // namespace

TEST_CASE("half-space vertices of a box are its eight corners") {
  const auto planes = unit_box_planes();
  const Result<std::vector<Point3d>> vertices =
      detail::half_space_vertices(planes, 1e-9);
  REQUIRE(vertices.has_value());
  CHECK(vertices.value().size() == 8u);
  for (const Point3d& p : vertices.value()) {
    CHECK((p[0] == 0.0 || p[0] == 2.0));
    CHECK((p[1] == 0.0 || p[1] == 3.0));
    CHECK((p[2] == 0.0 || p[2] == 4.0));
  }
}

TEST_CASE("half-space vertex order does not depend on plane order") {
  auto planes = unit_box_planes();
  const Result<std::vector<Point3d>> reference =
      detail::half_space_vertices(planes, 1e-9);
  REQUIRE(reference.has_value());
  std::swap(planes[0], planes[3]);
  std::swap(planes[1], planes[4]);
  const Result<std::vector<Point3d>> shuffled =
      detail::half_space_vertices(planes, 1e-9);
  REQUIRE(shuffled.has_value());
  CHECK(shuffled.value().size() == reference.value().size());
  for (const Point3d& p : shuffled.value()) {
    bool found = false;
    for (const Point3d& q : reference.value()) {
      found = found || (p[0] == q[0] && p[1] == q[1] && p[2] == q[2]);
    }
    CHECK(found);
  }
}

TEST_CASE("half-space vertices ignore a redundant plane") {
  auto planes = unit_box_planes();
  planes.push_back({{0.0, 0.0, 1.0}, 10.0});
  const Result<std::vector<Point3d>> vertices =
      detail::half_space_vertices(planes, 1e-9);
  REQUIRE(vertices.has_value());
  CHECK(vertices.value().size() == 8u);
}

TEST_CASE("half-space vertices report an empty intersection") {
  const std::array<detail::HalfSpacePlane, 2> planes{{
      {{1.0, 0.0, 0.0}, 0.0},
      {{-1.0, 0.0, 0.0}, -1.0},
  }};
  const Result<std::vector<Point3d>> vertices =
      detail::half_space_vertices(planes, 1e-9);
  REQUIRE(vertices.has_value());
  CHECK(vertices.value().empty());
}

TEST_CASE("half-space vertices reject a zero normal") {
  const std::array<detail::HalfSpacePlane, 1> planes{{
      {{0.0, 0.0, 0.0}, 1.0},
  }};
  const Result<std::vector<Point3d>> vertices =
      detail::half_space_vertices(planes, 1e-9);
  REQUIRE_FALSE(vertices.has_value());
  CHECK(vertices.error().code == Errc::degenerate);
}
```

- [ ] **Step 2: Create the registration and run the test to verify it fails**

Append `unit/half_space_test.cpp` to `GEOM_TEST_SOURCES` in `tests/CMakeLists.txt`. Re-run CMake and build; expected failure is `C1083: Cannot open include file: 'geom/brep/detail/half_space.hpp'`.

- [ ] **Step 3: Implement the vertices function**

Create `include/geom/brep/detail/half_space.hpp`:

```cpp
#pragma once

#include <cmath>
#include <cstddef>
#include <span>
#include <vector>

#include <geom/core/error.hpp>
#include <geom/core/result.hpp>
#include <geom/primitives/point.hpp>
#include <geom/primitives/vec.hpp>

namespace geom {
inline namespace v1 {
namespace detail {

struct HalfSpacePlane {
  Vec3d normal{};
  double offset = 0.0;
};

Result<std::vector<Point3d>> half_space_vertices(
    std::span<const HalfSpacePlane> planes, double tolerance);

}  // namespace detail
}  // namespace v1
}  // namespace geom
```

Implementation rules (all functions in this header are `inline`):

- Reject non-finite `tolerance <= 0` with `Errc::tolerance_violation`, message `"half_space: tolerance must be finite and positive"`.
- Normalize each plane: `length = sqrt(dot(n,n))`; non-finite or `length <= 0` → `Errc::degenerate`, message `"half_space: plane normal is zero"`; `offset /= length`.
- Enumerate every triple `(i, j, k)`. Solve with Cramer's rule using `double`: `det = dot(n_i, cross(n_j, n_k))`, `p = (d_i*cross(n_j,n_k) + d_j*cross(n_k,n_i) + d_k*cross(n_i,n_j)) / det`; skip when `|det| <= 1e-14` or the determinant or any coordinate is not finite.
- Accept the point when every plane satisfies `dot(n, p) <= offset + slack`, where `slack = max(tolerance, 1e-12 * scale)` and `scale = max(1, max |offset|, max |p|)`.
- Deduplicate accepted points with a linear scan using `slack`; keep the first representative so iteration order is deterministic for a fixed plane order.
- Return the vector; empty is valid.

- [ ] **Step 4: Run the test to verify it passes**

Run: `cmake --build build --config Debug --target geom_tests -- /v:minimal` then `.\build\tests\Debug\geom_tests.exe --source-file="*half_space*"`. Expected: 5 cases pass.

- [ ] **Step 5: Commit**

```bash
git add include/geom/brep/detail/half_space.hpp tests/unit/half_space_test.cpp tests/CMakeLists.txt
git commit -m "feat(brep): enumerate convex half-space vertices"
```

---

### Task 2: Half-space faces

**Files:**
- Modify: `include/geom/brep/detail/half_space.hpp`
- Modify: `tests/unit/half_space_test.cpp`

**Interfaces:**
- Consumes: `HalfSpacePlane`, `half_space_vertices`.
- Produces: `Result<std::vector<AssembledFace>> half_space_faces(std::span<const HalfSpacePlane> planes, std::span<const Point3d> vertices, double tolerance);` Each returned face has an outward `normal` and a counter-clockwise `outer` ring. Planes with fewer than three distinct incident vertices are skipped (redundant planes).

- [ ] **Step 1: Write the failing tests**

Append:

```cpp
TEST_CASE("half-space faces of a box are six outward quads") {
  const auto planes = unit_box_planes();
  const Result<std::vector<Point3d>> vertices =
      detail::half_space_vertices(planes, 1e-9);
  REQUIRE(vertices.has_value());
  const Result<std::vector<AssembledFace>> faces =
      detail::half_space_faces(planes, vertices.value(), 1e-9);
  REQUIRE(faces.has_value());
  CHECK(faces.value().size() == 6u);
  for (const AssembledFace& face : faces.value()) {
    CHECK(face.outer.size() == 4u);
    CHECK(face.holes.empty());
  }
  const Result<Brep> solid =
      assemble_solid(faces.value(), 1e-9);
  REQUIRE(solid.has_value());
  const Result<MassProperties> mass = mass_properties(solid.value());
  REQUIRE(mass.has_value());
  CHECK(std::abs(mass.value().volume - 24.0) <= 1e-9);
}
```

Add the includes `<geom/brep/assemble.hpp>`, `<geom/brep/mass_properties.hpp>`, and `<geom/brep/validate.hpp>` to the test file.

- [ ] **Step 2: Run the test to verify it fails**

Expected: `half_space_faces` is not a member of `detail`.

- [ ] **Step 3: Implement face extraction**

In `half_space.hpp`, add:

- For each plane `i`, collect indices of vertices with `|dot(n_i, p) - offset_i| <= slack`.
- If fewer than three, skip.
- Build an orthonormal basis in the plane: `u` is the normalized component of any non-parallel axis (prefer `+x`, then `+y`, then `+z`), `v = cross(n_i, u)`.
- Sort the ids by `atan2(dot(p - c, v), dot(p - c, u))` around the centroid `c`; this order is counter-clockwise about `n_i`.
- Emit `AssembledFace{{points in sorted order}, {}, n_i}`.
- Return `Errc::degenerate` with `"half_space: fewer than three planes bound a cell"` if no face was emitted.

- [ ] **Step 4: Run the tests**

Run the same filtered command. Expected: 6 cases pass, box volume 24.

- [ ] **Step 5: Commit**

```bash
git add include/geom/brep/detail/half_space.hpp tests/unit/half_space_test.cpp
git commit -m "feat(brep): extract outward half-space face rings"
```

---

### Task 3: Half-space solid

**Files:**
- Modify: `include/geom/brep/detail/half_space.hpp`
- Modify: `tests/unit/half_space_test.cpp`

**Interfaces:**
- Consumes: `half_space_vertices`, `half_space_faces`.
- Produces: `Result<Brep> half_space_solid(std::span<const HalfSpacePlane> planes, double tolerance);` Empty intersection is `Errc::degenerate` with `"half_space: planes have no common point"`; an assembly failure is returned unchanged.

- [ ] **Step 1: Write the failing tests**

Append:

```cpp
TEST_CASE("half-space solid of a right tetrahedron has the analytic volume") {
  const double s = std::sqrt(3.0);
  const std::array<detail::HalfSpacePlane, 4> planes{{
      {{-1.0, 0.0, 0.0}, 0.0},
      {{0.0, -1.0, 0.0}, 0.0},
      {{0.0, 0.0, -1.0}, 0.0},
      {{1.0, 1.0, 1.0}, 2.0},
  }};
  const Result<Brep> solid = detail::half_space_solid(planes, 1e-9);
  REQUIRE(solid.has_value());
  CHECK(validate_topology(solid.value()).ok());
  const Result<MassProperties> mass = mass_properties(solid.value());
  REQUIRE(mass.has_value());
  CHECK(std::abs(mass.value().volume - 4.0 / 3.0) <= 1e-9);
}

TEST_CASE("half-space solid rejects an empty cell") {
  const std::array<detail::HalfSpacePlane, 2> planes{{
      {{1.0, 0.0, 0.0}, 0.0},
      {{-1.0, 0.0, 0.0}, -1.0},
  }};
  const Result<Brep> solid = detail::half_space_solid(planes, 1e-9);
  REQUIRE_FALSE(solid.has_value());
  CHECK(solid.error().code == Errc::degenerate);
  CHECK(std::string_view(solid.error().message) ==
        "half_space: planes have no common point");
}
```

The tetrahedron expected volume is `2^3 / 6 = 4/3`.

- [ ] **Step 2: Run the test to verify it fails**

Expected: `half_space_solid` is not a member of `detail`.

- [ ] **Step 3: Implement the solid function**

Call `half_space_vertices`; if empty return the degenerate error; call `half_space_faces`; return `assemble_solid(faces, tolerance)`.

- [ ] **Step 4: Run the tests**

Expected: 8 cases pass.

- [ ] **Step 5: Commit**

```bash
git add include/geom/brep/detail/half_space.hpp tests/unit/half_space_test.cpp
git commit -m "feat(brep): assemble a convex half-space cell"
```

---

### Task 4: Shell through a half-space cavity

**Files:**
- Modify: `include/geom/brep/shell.hpp`
- Modify: `include/geom/brep/detail/half_space.hpp` (only if a helper is missing)
- Modify: `tests/unit/shell_test.cpp`

**Interfaces:**
- Public API unchanged: `ShellOptions{thickness, tolerance}`, `ShellResult{Brep brep; uint32_t cavity_shells}`, `shell_solid(const Brep&, std::span<const FaceId>, const ShellOptions&)`.
- New internal step: `BlendPlanarSolid` from `detail::blend_planar_convex_solid` supplies `faces`, `outward`, `offset`, `scale`.

- [ ] **Step 1: Write the failing tests**

In `tests/unit/shell_test.cpp`, delete the five cases that expect the M4c-1 `unsupported` messages and replace them with the cases below:

- `"shell refuses an open face"` → `"shell opens a box through one face"`
- `"shell rejects opening every face"` → `"shell rejects opening every face"` with the new message
- `"shell refuses an oblique face"` → delete (the tetrahedron case covers oblique faces)
- `"shell refuses any oblique face in any face order"` → delete (Task 5 covers face order)
- `"a duplicate open face request is still refused"` → `"a duplicate open face request is accepted once"`

```cpp
TEST_CASE("shell opens a box through one face") {
  const Brep box = box_solid(Point3d{0.0, 0.0, 0.0}, Point3d{2.0, 3.0, 4.0});
  const std::array<FaceId, 1> open{top_face(box)};
  ShellOptions options;
  options.thickness = 0.25;
  options.tolerance = 1e-9;
  const Result<ShellResult> result = shell_solid(box, open, options);
  REQUIRE(result.has_value());
  CHECK(result.value().cavity_shells == 0u);
  CHECK(result.value().brep.shell_count() == 1u);
  CHECK(validate_topology(result.value().brep).ok());
  const Result<MassProperties> mass = mass_properties(result.value().brep);
  REQUIRE(mass.has_value());
  CHECK(std::abs(mass.value().volume - 9.9375) <= 1e-9);
}

TEST_CASE("shell opens a box through two opposite faces") {
  const Brep box = box_solid(Point3d{0.0, 0.0, 0.0}, Point3d{2.0, 3.0, 4.0});
  const std::array<FaceId, 2> open{top_face(box), bottom_face(box)};
  ShellOptions options;
  options.thickness = 0.25;
  options.tolerance = 1e-9;
  const Result<ShellResult> result = shell_solid(box, open, options);
  REQUIRE(result.has_value());
  CHECK(result.value().cavity_shells == 0u);
  CHECK(std::abs(volume_of(result.value().brep) - 9.0) <= 1e-9);
}

TEST_CASE("shell of an oblique tetrahedron matches the inset volume") {
  const Brep tetra = right_tetrahedron(2.0);
  const std::array<FaceId, 0> open{};
  ShellOptions options;
  options.thickness = 0.1;
  options.tolerance = 1e-9;
  const Result<ShellResult> result = shell_solid(tetra, open, options);
  REQUIRE(result.has_value());
  CHECK(result.value().cavity_shells == 1u);
  const double outer = 4.0 / 3.0;
  const double s = 2.0 - options.thickness * (3.0 + std::sqrt(3.0));
  const double inner = s * s * s / 6.0;
  CHECK(std::abs(volume_of(result.value().brep) - (outer - inner)) <= 1e-9);
}

TEST_CASE("shell rejects opening every face") {
  const Brep box = box_solid(Point3d{0.0, 0.0, 0.0}, Point3d{2.0, 3.0, 4.0});
  const auto faces = box.face_ids();
  ShellOptions options;
  options.thickness = 0.25;
  options.tolerance = 1e-9;
  const Result<ShellResult> result = shell_solid(box, faces, options);
  REQUIRE_FALSE(result.has_value());
  CHECK(result.error().code == Errc::unsupported);
  CHECK(std::string_view(result.error().message) ==
        "shell_solid: at least one face must stay closed");
}

TEST_CASE("a duplicate open face request is accepted once") {
  const Brep box = box_solid(Point3d{0.0, 0.0, 0.0}, Point3d{2.0, 3.0, 4.0});
  const FaceId top = top_face(box);
  const std::array<FaceId, 2> open{top, top};
  ShellOptions options;
  options.thickness = 0.25;
  options.tolerance = 1e-9;
  const Result<ShellResult> result = shell_solid(box, open, options);
  REQUIRE(result.has_value());
  CHECK(result.value().cavity_shells == 0u);
  CHECK(std::abs(volume_of(result.value().brep) - 9.9375) <= 1e-9);
}
```

Add these helpers in the existing anonymous namespace, before the first test:

```cpp
FaceId top_face(const Brep& brep) {
  for (const FaceId face : brep.face_ids()) {
    const auto& plane =
        std::get<Plane3>(brep.geometry().surfaces[brep.face(face).surface]);
    const Vec3d normal = brep.face(face).reversed ? plane.normal * -1.0 : plane.normal;
    if (normal[2] > 0.5) return face;
  }
  FAIL("box has no +Z face");
  return FaceId{};
}

FaceId bottom_face(const Brep& brep) {
  for (const FaceId face : brep.face_ids()) {
    const auto& plane =
        std::get<Plane3>(brep.geometry().surfaces[brep.face(face).surface]);
    const Vec3d normal = brep.face(face).reversed ? plane.normal * -1.0 : plane.normal;
    if (normal[2] < -0.5) return face;
  }
  FAIL("box has no -Z face");
  return FaceId{};
}

Brep right_tetrahedron(double a) {
  const auto p = [](double x, double y, double z) { return Point3d{x, y, z}; };
  const std::array<AssembledFace, 4> faces{{
      {{p(0, 0, 0), p(0, 0, a), p(0, a, 0)}, {}, {-1, 0, 0}},
      {{p(0, 0, 0), p(a, 0, 0), p(0, 0, a)}, {}, {0, -1, 0}},
      {{p(0, a, 0), p(a, 0, 0), p(0, 0, 0)}, {}, {0, 0, -1}},
      {{p(a, 0, 0), p(0, a, 0), p(0, 0, a)}, {}, {1 / std::sqrt(3.0), 1 / std::sqrt(3.0),
                                                  1 / std::sqrt(3.0)}},
  }};
  const Result<Brep> tetra = assemble_solid(faces, 1e-9);
  REQUIRE(tetra.has_value());
  return tetra.value();
}
```

The tetrahedron is the review case for an oblique face: its inset planes are `x=t`, `y=t`, `z=t`, `x+y+z=a-t*sqrt(3)`, so the inner volume is `(a - t*(3+sqrt(3)))^3 / 6`.

- [ ] **Step 2: Run the tests to verify they fail**

Expected: the open box and tetra cases fail with the current `unsupported` messages.

- [ ] **Step 3: Implement the general cavity**

In `shell.hpp`, keep the options and `blend_planar_convex_solid` validation, and:

1. Include `geom/brep/detail/half_space.hpp` and rewrite the file/`ShellResult` comments: M4c-2 supports any closed convex planar solid and open faces; `cavity_shells` is the number of enclosed pockets after the cut and is `0` when any face is open.
2. Collect `open` as a sorted unique set of face indices; a duplicate open face is accepted once.
3. If every face is open, return `Errc::unsupported`, message `"shell_solid: at least one face must stay closed"`.
4. If no face is open and every face normal is axis-aligned, keep the existing inset-box fast path.
5. Otherwise compute the input AABB from its live vertices. Expand it by `margin = max({thickness, 1e-6 * max(1.0, extent), 1000.0 * tolerance})` on every axis.
6. Build `std::vector<detail::HalfSpacePlane>`: six seed-box planes with outward normals `±x/±y/±z`, then one plane per shell face with `normal = outward[i]`, `offset = offset[i] + (open ? thickness : -thickness)`.
7. Call `detail::half_space_solid(planes, tolerance)`; on failure map to `Errc::degenerate`, message `"shell_solid: wall thickness leaves no measurable interior at this tolerance"` unless the failure is `non_finite`.
8. Subtract: `boolean_solid(brep, cavity, SolidBooleanOp::difference, boolean_options)` with `boolean_options.tolerance = tolerance`; on failure map the same way except `non_finite`.
9. Require one resulting solid; run `validate_topology`; set `cavity_shells = open.empty() ? shells - 1 : 0` and require `shells >= 2` when closed.

The seed-box margin is not cosmetic: with an open face the cavity half-space is relaxed outward, and the seed box is what keeps the cavity bounded. It must extend past the input by at least `thickness` so the pocket breaks the surface.

- [ ] **Step 4: Run the tests to verify they pass**

Run:

```
cmake --build build --config Debug --target geom_tests -- /v:minimal
.\build\tests\Debug\geom_tests.exe --test-case="*shell*"
cmake --build build --config Release --target geom_tests -- /v:minimal
.\build\tests\Release\geom_tests.exe --test-case="*shell*"
```

Expected: the new cases pass, closed-box volume remains 10.875, and no test still expects the old `unsupported` open/oblique messages.

- [ ] **Step 5: Commit**

```bash
git add include/geom/brep/shell.hpp tests/unit/shell_test.cpp
git commit -m "feat(brep): shell convex solids through a half-space cavity"
```

---

### Task 5: Order and diagnostic coverage

**Files:**
- Modify: `tests/unit/shell_test.cpp`

**Interfaces:**
- Consumes: the `shell_solid` behavior from Task 4.

- [ ] **Step 1: Extend the existing 720-order test to cover the open case**

The existing `"closed axis-aligned shell is independent of face order"` case already builds `solid` from the permuted `AssembledFace` list. Inside that loop, after the closed-shell checks, add:

```cpp
    const std::array<FaceId, 1> open{top_face(solid.value())};
    const Result<ShellResult> opened = shell_solid(solid.value(), open, options);
    if (!opened.has_value()) {
      first_failure = opened.error().message;
      break;
    }
    if (std::abs(volume_of(opened.value().brep) - 9.9375) > 1e-9) {
      first_failure = "wrong open volume";
      break;
    }
```

Keep the existing `target` of 720 in Release and 24 in Debug, and the existing `INFO` diagnostic. Do not add a `permute_faces` or `face_of` helper: the test already has the permuted solid in hand, and `top_face` finds its open face by outward normal.

- [ ] **Step 2: Add diagnostics that must not leak**

```cpp
TEST_CASE("shell rejects a wall at or beyond the half-thickness") {
  const Brep box = box_solid(Point3d{0.0, 0.0, 0.0}, Point3d{2.0, 3.0, 4.0});
  ShellOptions options;
  options.thickness = 1.0;
  options.tolerance = 1e-9;
  const Result<ShellResult> result = shell_solid(box, {}, options);
  REQUIRE_FALSE(result.has_value());
  CHECK(result.error().code == Errc::degenerate);
}

TEST_CASE("shell rejects an oblique wall that leaves no interior") {
  const Brep tetra = right_tetrahedron(2.0);
  ShellOptions options;
  options.thickness = 0.7;
  options.tolerance = 1e-9;
  const Result<ShellResult> result = shell_solid(tetra, {}, options);
  REQUIRE_FALSE(result.has_value());
  CHECK(result.error().code == Errc::degenerate);
}
```

The tetrahedron's inradius is `a*sqrt(6)/12 = 0.408`, so `t = 0.7` must fail.

- [ ] **Step 3: Run the tests and record the actual counts**

Run `*shell*` in Debug and Release. Update any count expectations in the plan only after observing the real numbers; do not guess.

- [ ] **Step 4: Commit**

```bash
git add tests/unit/shell_test.cpp
git commit -m "test(brep): pin shell order independence and diagnostics"
```

---

### Task 6: Watertight and deterministic open shell

**Files:**
- Modify: `tests/unit/blend_qa_test.cpp`

**Interfaces:**
- Consumes: `shell_solid`, `tessellate`, `mass_properties`.

- [ ] **Step 1: Extend the QA case**

In the second QA case, after the closed shell checks, add an open-top box:

```cpp
const std::array<FaceId, 1> open_top{top_face(box)};
const Result<ShellResult> pocket_a = shell_solid(box, open_top, shell);
const Result<ShellResult> pocket_b = shell_solid(box, open_top, shell);
REQUIRE(pocket_a.has_value());
REQUIRE(pocket_b.has_value());
CHECK(pocket_a.value().cavity_shells == 0u);
CHECK(std::abs(volume_of(pocket_a.value().brep) - 4.0625) <= 1e-9);
const Result<TriMesh> pocket_mesh_a =
    tessellate(pocket_a.value().brep, TessellationOptions{});
const Result<TriMesh> pocket_mesh_b =
    tessellate(pocket_b.value().brep, TessellationOptions{});
REQUIRE(pocket_mesh_a.has_value());
REQUIRE(pocket_mesh_b.has_value());
CHECK(boundary_edges(pocket_mesh_a.value()) == 0u);
CHECK(same_positions(pocket_mesh_a.value().positions,
                     pocket_mesh_b.value().positions));
CHECK(pocket_mesh_a.value().triangles == pocket_mesh_b.value().triangles);
```

The QA box is `2x2x2`; an open top with `t=0.25` removes `1.5*1.5*1.75 = 3.9375` and leaves `8 - 3.9375 = 4.0625`. The test must assert `4.0625`, not the removed volume. Add `#include <variant>`, `#include <geom/surfaces/plane.hpp>`, and a local `top_face` helper in the anonymous namespace (same outward-normal scan as `shell_test.cpp`) because `blend_qa_test.cpp` does not have one.

- [ ] **Step 2: Run the QA file in both configurations**

Expected: 2 cases pass; assertions increase from the M4c-1 count.

- [ ] **Step 3: Commit**

```bash
git add tests/unit/blend_qa_test.cpp
git commit -m "test(brep): prove open shells are watertight and deterministic"
```

---

### Task 7: Benchmark, ADR, architecture and Results

**Files:**
- Modify: `bench/bench_m4c1.cpp`
- Create: `docs/benchmarks/m4c2-shell.md`
- Modify: `docs/adr/0041-planar-blends-via-half-space-cuts.md`
- Modify: `docs/architecture/geometry-library-architecture.md` (§8.9 and §16.4)
- Modify: `docs/superpowers/plans/2026-11-26-m4c2-open-shell.md` (Results table)

- [ ] **Step 1: Add open and oblique shell workloads to the bench**

The 100³ box keeps the closed case. Add `#include <geom/brep/assemble.hpp>`, `#include <geom/surfaces/plane.hpp>`, a local `top_face(const geom::Brep&)` helper (outward-normal scan), a local `right_tetrahedron(double)` helper (the same four `AssembledFace`s as `shell_test.cpp`), and:

```cpp
// Open top: the cavity protrudes through +Z and the outer box clips it.
const std::array<geom::FaceId, 1> open{top_face(box.value())};
geom::ShellOptions open_options;
open_options.thickness = 1.0;
open_options.tolerance = 1e-9;
const auto opened = geom::shell_solid(box.value(), open, open_options);
if (!opened.has_value()) return 3;
const double open_volume =
    geom::mass_properties(opened.value().brep).value().volume;
// Removed: 98 * 98 * 99 = 950,796; shell = 1,000,000 - 950,796 = 49,204.
if (!close_to(open_volume, 49204.0)) return 4;

// Oblique tetrahedron, t = 0.1: analytic shell volume.
const geom::Brep tetra = right_tetrahedron(2.0);
geom::ShellOptions tetra_options;
tetra_options.thickness = 0.1;
tetra_options.tolerance = 1e-9;
const auto tet = geom::shell_solid(tetra, {}, tetra_options);
if (!tet.has_value()) return 5;
const double expected =
    4.0 / 3.0 -
    std::pow(2.0 - 0.1 * (3.0 + std::sqrt(3.0)), 3) / 6.0;
if (!close_to(geom::mass_properties(tet.value().brep).value().volume, expected)) return 6;
```

Report min/median/max for each workload; fail with a non-zero exit on a geometry mismatch.

- [ ] **Step 2: Write the benchmark record**

`docs/benchmarks/m4c2-shell.md` must state: machine and compiler version, the three workloads, observed volumes, min/median/max milliseconds, and the reproduction command. Do not copy the M4c-1 numbers.

- [ ] **Step 3: Update the contract documents**

- ADR-0041: add an "M4c-2 extension" section: open faces and oblique convex faces are supported by a half-space cavity; the M4c-1 axis-aligned limitation and the `unsupported` messages are superseded; fillet remains deferred.
- Architecture §8.9: replace the M4c-1-only sentence with the M4c-2 contract and keep the fillet deferral.
- Architecture §16.4: append the M4c-2 evidence paragraph with the actual commit hashes, counts, and volumes.
- Plan Results: fill the measured counts, Debug/Release status, and the known Debug PENDING caveat.

- [ ] **Step 4: Run the full Release suite**

Run `.\build\tests\Release\geom_tests.exe`; expected 0 failures. Record the exact case/assertion count. Full Debug remains PENDING for the eight known long `brep_boolean_curved_test.cpp` cases; run the Debug lane `*shell*,*half_space*,*blend*,*editor*` instead and record it.

- [ ] **Step 5: Commit**

```bash
git add bench/bench_m4c1.cpp docs/benchmarks/m4c2-shell.md docs/adr/0041-planar-blends-via-half-space-cuts.md docs/architecture/geometry-library-architecture.md docs/superpowers/plans/2026-11-26-m4c2-open-shell.md
git commit -m "docs(m4c2): record open-shell evidence"
```

---

## Results

| Item | Value |
| --- | --- |
| Debug `*shell*` | 51 cases / 524 assertions, 0 failed, 1,099 skipped |
| Debug lane `*shell*,*half_space*,*blend*,*editor*` | 91 cases / 946 assertions, 0 failed, 1,059 skipped |
| Release `*shell*` | 51 cases / 1,916 assertions, 0 failed, 1,099 skipped |
| Full Release | 1,150 cases / 682,717 assertions, 0 failed, 0 skipped (exit 0, ~12m40s) |
| Open top box `2x3x4`, t=0.25 | 9.9375 |
| Opposite open box, t=0.25 | 9.0 |
| Oblique tetra, a=2, t=0.1 | 0.740147367097487 (analytic) |
| Bench 100^3 closed / open | 58,808 / 49,204 |

The QA file now holds 54 assertions (46 in
`blend outputs are deterministic and leave inputs unchanged` + 8 in
`chamfer and shell outputs tessellate watertight`); this supersedes the
52-assertion figure in `task-6-report.md`.
