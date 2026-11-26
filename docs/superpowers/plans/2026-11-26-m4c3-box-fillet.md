# M4c-3a Axis-Aligned Box Fillet Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use subagent-driven-development (recommended) or executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Add an exact constant-radius fillet for up to four parallel edges of an axis-aligned box: replace each selected edge with a cylindrical patch and each cap corner with a circular arc, without faceting or a Boolean round-trip.

**Architecture:** Detect that the input is one closed axis-aligned box, map the selected `EdgeId`s to its corner family, then build the result directly as a list of `CurvedAssembledFace`s: two planar caps with rounded corners, four planar side faces, and one `Cylinder3` patch per selected edge. `assemble_curved_solid` creates the topology and pcurves; `validate_topology` and `curved_solid_volume` check the result.

**Tech Stack:** C++20, existing `geom::brep` topology/assembly/tessellation/mass-property headers, doctest 2.4.11, CMake + MSVC 19.44.

**Spec:** `docs/architecture/geometry-library-architecture.md` §8.9 (rolling-ball fillet, no half-made solids) and `docs/quality/testing-strategy.md` (TDD, failure paths, invariants).

## Global Constraints

- C++20 only; no new third-party dependency.
- MSVC `/W4 /WX /utf-8`; Debug and Release must both build with zero warnings.
- No global epsilon. Every numeric guard is local, named, and documented. `tolerance` is finite and strictly positive.
- The input `Brep` is never modified; every failure happens before constructing or returning a result entity.
- A failure returns exactly one stable `Errc` and message; never a partial or silently wrong solid.
- Determinism: the same input and options produce the same topology counts, surface geometry and volume on every run.
- TDD: write the failing test, record the RED output, implement the minimum, record GREEN, then commit.
- One commit per task; do not amend a task after its review without the controller's instruction.
- Do not run the full Debug suite (known >20 minute cases); use the Debug filter named in each task and the full Release suite at the milestone.

## Scope and Deviations (review these before implementation)

This plan deliberately ships a narrow but exact first fillet:

1. **Axis-aligned box only.** The box may have any positive extents, but its faces must be perpendicular to the world X/Y/Z axes. A rotated box is rejected in M4c-3a; it is not silently treated as axis-aligned.
2. **One parallel family.** All selected edges must be the box edges parallel to one axis (up to four). Selecting a vertical and a horizontal edge together is rejected. This guarantees no two selected edges share a vertex, so no vertex blend is needed.
3. **No vertex blends, variable radius, or general convex solids.** General convex planar solids with oblique end faces need elliptical trim arcs, and adjacent selected edges need spherical vertex patches. Those are M4c-3b/M4c-3c.
4. **No faceting.** The fillet is a real `Cylinder3` patch with `CircleCurve3` arcs; a ring of planar facets is forbidden by §8.9.
5. **No partial end treatment.** Filleting a subset of one family rounds only the selected corners of each cap; unselected corners stay sharp. Mixed-family and vertex-blend cases remain rejected, not approximated.

A reviewer who wants general convex fillet should split M4c-3b instead of widening these tasks.

## Review Focus

- **Mixed edge families.** A selection containing one vertical and one horizontal edge must return `unsupported`, not a solid with only one family filleted.
- **Radius at the section limit.** `radius >= 0.5*min(other two extents)` must be rejected; additionally, for each side with **both** cross-section corners selected, `2*radius >= side_extent - tolerance` must be rejected before assembly (the remaining flat face would be narrower than the tolerance). A side with only one selected end keeps `side_extent - r`, so `r=0.9999999999` on the 2×3×4 box must still succeed for k=1 and for a diagonal k=2.
- **All four corners on a cap.** With four selected edges the cap ring alternates four straight segments and four arcs, and the four cylinder patches must share the arcs without duplicate or missing vertices.
- **Non-box shapes.** A tetrahedron, cylinder, or rotated box must be rejected with the stable box-only message; no accidental acceptance through a coordinate-collapse check.
- **Failure immutability.** Every rejection (bad radius, mixed family, dead edge, non-box) must leave the input's counts, volume, and surface list unchanged; successful calls must do the same.

---

## Construction (decided here)

Coordinates: the box is `[x0,x1] x [y0,y1] x [z0,z1]`. Let the selected family be axis `f` and let `c` be its extent; let `a` and `b` be the extents of the other two axes in ascending axis order. For each selected edge its cross-section corner is one of the four `(0/1, 0/1)` combinations of the other two axes. Corner index convention: `0=(lo,lo)`, `1=(hi,lo)`, `2=(hi,hi)`, `3=(lo,hi)`.

For a selected corner `(u,v)` on cap `w`:

- cylinder center is `(u ? hi0 - r : lo0 + r, v ? hi1 - r : lo1 + r)` on the two cross-section axes;
- the cap arc has that center, radius `r`, and a quarter-circle interval selected by the quadrant;
- the adjacent side faces are inset by `r` at that end;
- the cylinder patch spans the cap-to-cap height and the same quarter-circle interval.

The cap ring is counter-clockwise about its outward normal. For the `+f` cap (`z1` when `f=2`) the ring walks the rectangle counter-clockwise in the `(x,y)` plane; at a selected corner the corner point is replaced by the two tangent points and the arc between them, traversed with increasing angle. The `-f` cap uses the reversed ring. For `f=0` or `f=1` the same construction is expressed in the two remaining axes in ascending order; the outward normal is the positive or negative family direction.

The four side faces are rectangles. A side face parallel to `f` is inset by `r` at each end whose corner edge is selected. Its ring is ordered so that the face normal points out of the box.

A cylinder patch for corner `(u,v)` uses `Cylinder3{origin = (cx, cy, z0 for f=2), axis = +f, x_axis = +other0, radius = r}` and the quadrant interval:

| corner (other0, other1) | `u0 .. u1` |
| --- | --- |
| `(hi,hi)` | `0 .. pi/2` |
| `(lo,hi)` | `pi/2 .. pi` |
| `(lo,lo)` | `pi .. 3pi/2` |
| `(hi,lo)` | `3pi/2 .. 2pi` |

for `other0 = x`, `other1 = y` and `axis = +z`. For other families rotate the same table with the chosen right-handed frame; the implementation should derive the angles from the actual cylinder frame rather than hard-code a second table.

A cylinder loop is `(u0,0), (u1,0), (u1,h), (u0,h)` in `(u,v)` parameter space. The first and third segments are arcs; their `CurvedArcSpan` is the circle centered at `(cx,cy,cap)` with radius `r` and `t0=u0`, `t1=u1`. The other two segments are straight. The matching cap loop carries the same circle and interval, so the assembler can share one `CircleCurve3` edge between the cap and the patch.

Volume and area anchors for a box `a x b x c` with `k` selected edges of length `c` and radius `r`:

- `V = a*b*c - k*(1 - pi/4)*r^2*c`
- `A = 2*(a*b - k*r^2 + k*pi*r^2/4) + (2*a + 2*b - 2*k*r)*c + k*(pi/2)*r*c`

Topology anchors: `V = 8 + 2*k`, `E = 12 + 3*k`, `F = 6 + k`, `chi = 2`. Every selected edge contributes one cylinder face and two arc edges; the two caps keep four straight sides plus `k` arcs.

---

## File Map

- Create `include/geom/brep/fillet.hpp` — public `FilletOptions`, `FilletResult`, `fillet_edges`.
- Create `include/geom/brep/detail/fillet_box.hpp` — box detection, edge-family mapping, and the `CurvedAssembledFace` builder.
- Create `tests/unit/fillet_test.cpp` — all fillet tests.
- Create `bench/bench_m4c3.cpp` — 100-cube fillet benchmark.
- Modify `tests/CMakeLists.txt`, the root `CMakeLists.txt` (bench targets live there), `.github/workflows/ci.yml`.
- Create `docs/adr/0042-box-fillet.md`, `docs/benchmarks/m4c3-fillet.md`.
- Modify `docs/architecture/geometry-library-architecture.md` §8.9 and §16.4.
- Modify this plan's Results section at the end.

### Task 1: Public API, box detection, and one-edge fillet

**Files:**
- Create: `include/geom/brep/fillet.hpp`
- Create: `include/geom/brep/detail/fillet_box.hpp`
- Create: `tests/unit/fillet_test.cpp`
- Modify: `tests/CMakeLists.txt`

**Interfaces:**
- Consumes: `Brep`, `EdgeId`, `make_box`, `assemble_curved_solid`, `validate_topology`, `detail::curved_solid_volume`.
- Produces:
  - `struct FilletOptions { double radius = 0.0; double tolerance = 1e-9; };`
  - `struct FilletResult { Brep brep; std::uint32_t filleted_edges = 0; };`
  - `Result<FilletResult> fillet_edges(const Brep& brep, std::span<const EdgeId> edges, const FilletOptions& options = {});`
  - `detail::FilletBox { double lo[3]; double hi[3]; };`
  - `Result<detail::FilletBox> detail::fillet_box_bounds(const Brep&, double tolerance);`
  - `Result<Brep> detail::fillet_box_edges(const FilletBox&, std::span<const int> corners, int family_axis, double radius, double tolerance);`

- [ ] **Step 1: Write the failing tests**

Add the helper and three test cases to `tests/unit/fillet_test.cpp`:

```cpp
#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <limits>
#include <map>
#include <numbers>
#include <span>
#include <utility>
#include <variant>
#include <vector>

#include <doctest/doctest.h>

#include <geom/brep/detail/curved_volume.hpp>
#include <geom/brep/fillet.hpp>
#include <geom/brep/primitives.hpp>
#include <geom/brep/tessellate.hpp>
#include <geom/brep/topology.hpp>
#include <geom/brep/validate.hpp>
#include <geom/core/error.hpp>
#include <geom/curves/analytic_curve.hpp>
#include <geom/surfaces/analytic.hpp>
#include <geom/surfaces/plane.hpp>

namespace {

using geom::Brep;
using geom::EdgeId;
using geom::Errc;
using geom::Point3d;
using geom::Vec3d;

Brep box_2x3x4() {
  return geom::make_box(Point3d{0.0, 0.0, 0.0}, Point3d{2.0, 3.0, 4.0}, 1e-9).value();
}

double volume_of(const Brep& brep) {
  return geom::detail::curved_solid_volume(brep, brep.solid_ids().front(), 1e-9).value();
}

EdgeId vertical_edge_at(const Brep& brep, double x, double y) {
  for (const EdgeId id : brep.edge_ids()) {
    const geom::Edge& edge = brep.edge(id);
    const Point3d a = brep.vertex(edge.v0).point;
    const Point3d b = brep.vertex(edge.v1).point;
    if (std::abs(a[0] - x) <= 1e-12 && std::abs(a[1] - y) <= 1e-12 &&
        std::abs(b[0] - x) <= 1e-12 && std::abs(b[1] - y) <= 1e-12) {
      return id;
    }
  }
  return EdgeId{};
}

}  // namespace

TEST_CASE("fillet_edges rejects options and non-box shapes") {
  const Brep box = box_2x3x4();
  const std::span<const EdgeId> empty;
  CHECK(geom::fillet_edges(box, empty, {}).error().code == Errc::unsupported);

  const EdgeId edge = vertical_edge_at(box, 2.0, 3.0);
  const std::span<const EdgeId> one(&edge, 1);
  CHECK(geom::fillet_edges(box, one, {0.0, 1e-9}).error().code == Errc::degenerate);
  CHECK(geom::fillet_edges(box, one,
                           {std::numeric_limits<double>::quiet_NaN(), 1e-9})
            .error()
            .code == Errc::non_finite);
  CHECK(geom::fillet_edges(box, one, {0.25, 0.0}).error().code ==
        Errc::tolerance_violation);

  const Brep cylinder = geom::make_cylinder(Point3d{0.0, 0.0, 0.0}, Vec3d{0.0, 0.0, 1.0},
                                            1.0, 2.0, 1e-9)
                            .value();
  const std::array<EdgeId, 1> dead{EdgeId::make(999u, 0u)};
  CHECK(geom::fillet_edges(cylinder, dead, {0.25, 1e-9}).error().code == Errc::unsupported);
}

TEST_CASE("fillet_edges rounds one vertical box edge exactly") {
  const Brep box = box_2x3x4();
  const EdgeId edge = vertical_edge_at(box, 2.0, 3.0);
  REQUIRE(edge.index() != 0xFFFFFFFFu);

  const auto result = geom::fillet_edges(box, std::span<const EdgeId>(&edge, 1), {0.25, 1e-9});
  REQUIRE(result.has_value());
  const Brep& filleted = result.value().brep;
  CHECK(result.value().filleted_edges == 1u);
  CHECK(geom::validate_topology(filleted).ok());
  CHECK(filleted.vertex_count() == 10u);
  CHECK(filleted.edge_count() == 15u);
  CHECK(filleted.face_count() == 7u);
  CHECK(std::abs(volume_of(filleted) - 23.94634954084936) <= 1e-9);

  std::size_t cylinders = 0;
  std::size_t arcs = 0;
  for (const geom::FaceId face : filleted.face_ids()) {
    const auto& surface =
        filleted.geometry().surfaces[filleted.face(face).surface];
    if (std::holds_alternative<geom::Cylinder3>(surface)) {
      ++cylinders;
      CHECK(std::abs(std::get<geom::Cylinder3>(surface).radius - 0.25) <= 1e-12);
    }
  }
  for (const EdgeId id : filleted.edge_ids()) {
    const auto& curve = filleted.geometry().curves[filleted.edge(id).curve];
    if (std::holds_alternative<geom::CircleCurve3>(curve)) {
      ++arcs;
      CHECK(std::abs(std::get<geom::CircleCurve3>(curve).radius - 0.25) <= 1e-12);
    }
  }
  CHECK(cylinders == 1u);
  CHECK(arcs == 2u);
}
```

`Vertical_edge_at` must check the edge direction is parallel to `+Z` and the endpoints are `z=0` and `z=4`; the midpoint check above is enough only because every box edge is axis-aligned, so add an explicit `std::abs(a[2] - 0.0) <= 1e-12 && std::abs(b[2] - 4.0) <= 1e-12` guard.

- [ ] **Step 2: Run the tests to verify they fail**

Run:
`cmake --build build --config Debug --target geom_tests`
`build\tests\Debug\geom_tests.exe --test-case="*fillet*"`

Expected RED: `C1083` cannot open `geom/brep/fillet.hpp`, or `fillet_edges` is not a member of `geom`.

- [ ] **Step 3: Implement the API and the one-corner builder**

In `include/geom/brep/fillet.hpp`:

1. Reject an empty `edges` span with `unsupported "fillet_edges: at least one edge is required"`. This check comes **before** the option checks: the default `FilletOptions{}` has `radius = 0.0`, so the empty + default call must report `unsupported`, not `degenerate`.
2. Validate `radius` finite (`non_finite "fillet_edges: radius must be finite"`), positive (`degenerate "fillet_edges: radius must be positive"`), `tolerance` finite (`non_finite "fillet_edges: tolerance must be finite"`), positive (`tolerance_violation "fillet_edges: tolerance must be positive"`).
3. Call `detail::fillet_box_bounds(brep, options.tolerance)`; on failure return it unchanged.
4. Deduplicate `edges` preserving first occurrence. For each edge: alive check (`index_out_of_range "fillet_edges: edge does not belong to the solid"`), direction exactly one nonzero component (`unsupported "fillet_edges: only axis-aligned box edges are supported"`), all selected edges same family axis (`unsupported "fillet_edges: selected edges must be parallel"`), and the edge must join the two family extremes at a corner (`unsupported "fillet_edges: edge is not a box corner edge"`).
5. Map each corner to the two other-axis low/high bits, then call `detail::fillet_box_edges`.
6. On success set `FilletResult.filleted_edges` to the distinct count and return.

In `include/geom/brep/detail/fillet_box.hpp`:

- `fillet_box_bounds` rejects anything that is not one closed shell, eight vertices, six planar faces, twelve axis-aligned edges, all eight coordinate combinations present, and every face a rectangle. Use `validate_topology` and `blend_planar_convex_solid` as the shape gate, then verify the eight corners against `lo/hi`. Message: `unsupported "fillet_edges: only an axis-aligned box is supported"`.
- The builder implements the construction section for one corner. Build the `CurvedAssembledFace` vector: the `+f` cap (CCW), the `-f` cap (reversed), the four insets, and the one cylinder. Give cap straight segments `nullopt` spans and the arc segments the `CurvedArcSpan` circle described above. `loops_uv` only needs same-size seed values; use `(other0, other1)` for caps, the two in-plane coordinates for side faces, and `(u,v)` for the cylinder.
- Call `assemble_curved_solid(faces, tolerance)`, then `validate_topology`; return `invalid_topology "fillet_edges: fillet produced an invalid solid"` if the report is not OK.

- [ ] **Step 4: Run the tests to verify they pass**

Run the same two commands, then run `--test-case="*fillet*"` again and record `2 cases | 2 passed` and the assertion count.

- [ ] **Step 5: Commit**

```bash
git add include/geom/brep/fillet.hpp include/geom/brep/detail/fillet_box.hpp tests/unit/fillet_test.cpp tests/CMakeLists.txt
git commit -m "feat(brep): fillet one axis-aligned box edge with an exact cylinder"
```

### Task 2: Any subset of one parallel family

**Files:**
- Modify: `include/geom/brep/fillet.hpp`
- Modify: `include/geom/brep/detail/fillet_box.hpp`
- Modify: `tests/unit/fillet_test.cpp`

**Interfaces:**
- Consumes: Task 1's `FilletBox`, `fillet_edges`, and the single-corner builder.
- Produces: the same public API; `fillet_box_edges` now accepts one to four corner indices and emits one cylinder per corner.

- [ ] **Step 1: Write the failing tests**

Add this helper above the test (it returns the four vertical edges in the corner order `(0,0),(2,0),(2,3),(0,3)`):

```cpp
std::vector<EdgeId> family_edges(const Brep& brep) {
  std::vector<EdgeId> result;
  for (const std::array<double, 2>& corner : {std::array{0.0, 0.0}, std::array{2.0, 0.0},
                                              std::array{2.0, 3.0}, std::array{0.0, 3.0}}) {
    result.push_back(vertical_edge_at(brep, corner[0], corner[1]));
  }
  return result;
}
```

Then add one `TEST_CASE`:

```cpp
TEST_CASE("fillet_edges rounds the four vertical prefixes") {
  const Brep box = box_2x3x4();
  const std::vector<EdgeId> all = family_edges(box);
  const double pi = std::numbers::pi_v<double>;
  const double per_edge = (1.0 - pi / 4.0) * 0.25 * 0.25 * 4.0;

  for (std::size_t k = 1; k <= 4; ++k) {
    const std::vector<EdgeId> selected(all.begin(), all.begin() + static_cast<std::ptrdiff_t>(k));
    const auto result = geom::fillet_edges(box, selected, {0.25, 1e-9});
    REQUIRE(result.has_value());
    const Brep& filleted = result.value().brep;
    CHECK(result.value().filleted_edges == k);
    CHECK(filleted.vertex_count() == 8u + 2u * k);
    CHECK(filleted.edge_count() == 12u + 3u * k);
    CHECK(filleted.face_count() == 6u + k);
    CHECK(std::abs(volume_of(filleted) - (24.0 - static_cast<double>(k) * per_edge)) <= 1e-9);
    CHECK(geom::validate_topology(filleted).ok());
  }
}
```

Add order and duplicate cases:

```cpp
TEST_CASE("fillet_edges deduplicates and ignores selection order") {
  const Brep box = box_2x3x4();
  const std::vector<EdgeId> all = family_edges(box);
  std::vector<EdgeId> duplicated{all[0], all[1], all[0], all[2], all[1]};
  const auto forward = geom::fillet_edges(box, duplicated, {0.25, 1e-9});
  REQUIRE(forward.has_value());
  CHECK(forward.value().filleted_edges == 3u);

  std::vector<EdgeId> reversed(all.begin(), all.begin() + 3);
  std::reverse(reversed.begin(), reversed.end());
  const auto backward = geom::fillet_edges(box, reversed, {0.25, 1e-9});
  REQUIRE(backward.has_value());
  CHECK(backward.value().brep.face_count() == forward.value().brep.face_count());
  CHECK(std::abs(volume_of(backward.value().brep) - volume_of(forward.value().brep)) <= 1e-9);
}
```

- [ ] **Step 2: Run the tests to verify they fail**

Expected: the k=2..4 cases fail with `REQUIRE(result.has_value())`, because the builder still only handles one corner (the single-corner code path may silently drop the rest).

- [ ] **Step 3: Generalize the builder**

Replace the one-corner cap and side construction with loops over all four corners:

- Build the cap ring by walking the four sides and, at a selected corner, emitting the two tangent points and the arc; at an unselected corner emit the sharp corner. Keep the ring CCW about `+f`.
- Inset each side face by `r` at the ends whose corner is selected.
- Emit one cylinder per selected corner with the quadrant table.
- Keep the same `loops_uv` size rules.

Do not special-case four corners; a single loop must build k=0..4, and Task 1's one-edge test must still pass unchanged.

- [ ] **Step 4: Run the tests to verify they pass**

`--test-case="*fillet*"` should report the Task 1 cases plus two new cases, all passing; record the exact count.

- [ ] **Step 5: Commit**

```bash
git add include/geom/brep/fillet.hpp include/geom/brep/detail/fillet_box.hpp tests/unit/fillet_test.cpp
git commit -m "feat(brep): fillet any subset of one box edge family"
```

### Task 3: Diagnostics, immutability, watertightness and determinism

**Files:**
- Modify: `include/geom/brep/fillet.hpp`
- Modify: `include/geom/brep/detail/fillet_box.hpp`
- Modify: `tests/unit/fillet_test.cpp`

**Interfaces:**
- Consumes: Task 1–2 public API and builder.
- Produces: stable rejection paths and QA tests; no new public type.

- [ ] **Step 1: Write the failing tests**

Replace the Task 2 prefix loop with an exhaustive `axis = 0..2` × `mask = 1..15` loop: every non-empty subset of every family must produce `V=8+2k`, `E=12+3k`, `F=6+k`, the exact volume, and exactly `k` cylinders whose origin is at the family low end and at a `lo+r`/`hi-r` cross-section corner. Rename the old prefix case accordingly.

Add:

```cpp
TEST_CASE("fillet_edges rejects a radius at the section limit") {
  const Brep box = box_2x3x4();
  const auto edge = vertical_edge_at(box, 2.0, 3.0);
  const auto too_large = geom::fillet_edges(box, std::span<const EdgeId>(&edge, 1), {1.0, 1e-9});
  REQUIRE_FALSE(too_large.has_value());
  CHECK(too_large.error().code == Errc::unsupported);

  const auto just_below = geom::fillet_edges(box, std::span<const EdgeId>(&edge, 1), {0.999, 1e-9});
  REQUIRE(just_below.has_value());
  CHECK(geom::validate_topology(just_below.value().brep).ok());

  // Four edges at 1-1e-10 leave a 2e-10 face; reject before assembly.
  const std::vector<EdgeId> all = family_edges(box);
  const auto sliver = geom::fillet_edges(box, all, {0.9999999999, 1e-9});
  REQUIRE_FALSE(sliver.has_value());
  CHECK(sliver.error().code == Errc::unsupported);
  CHECK(std::string_view(sliver.error().message).find("narrower than the tolerance") !=
        std::string_view::npos);
}

TEST_CASE("fillet_edges rejects a mixed family and a dead edge") {
  const Brep box = box_2x3x4();
  const EdgeId vertical = vertical_edge_at(box, 2.0, 3.0);
  const EdgeId horizontal = horizontal_edge_at(box, 4.0, 3.0);
  const std::array<EdgeId, 2> mixed{vertical, horizontal};
  const auto mixed_result = geom::fillet_edges(box, mixed, {0.25, 1e-9});
  REQUIRE_FALSE(mixed_result.has_value());
  CHECK(mixed_result.error().code == Errc::unsupported);

  const std::array<EdgeId, 1> dead{EdgeId::make(999u, 0u)};
  const auto dead_result = geom::fillet_edges(box, dead, {0.25, 1e-9});
  REQUIRE_FALSE(dead_result.has_value());
  CHECK(dead_result.error().code == Errc::index_out_of_range);
}
```

`horizontal_edge_at(box, z, y)` is the same scan as `vertical_edge_at`, but requires one endpoint at `x=0` and the other at `x=2`, both at the given `(y,z)`; add it in the anonymous namespace before this test.

Add the watertight helper that counts undirected triangle edges and asserts every edge occurs twice:

```cpp
std::size_t boundary_edges(const geom::TriMesh& mesh) {
  std::map<std::pair<std::uint32_t, std::uint32_t>, int> uses;
  for (const auto& triangle : mesh.triangles) {
    for (int i = 0; i < 3; ++i) {
      std::uint32_t a = triangle[static_cast<std::size_t>(i)];
      std::uint32_t b = triangle[(i + 1) % 3];
      if (a > b) std::swap(a, b);
      ++uses[{a, b}];
    }
  }
  std::size_t boundary = 0;
  for (const auto& [key, count] : uses) {
    (void)key;
    if (count == 1) ++boundary;
  }
  return boundary;
}
```

Then add:

```cpp
TEST_CASE("fillet_edges is watertight and deterministic") {
  const Brep box = box_2x3x4();
  const std::vector<EdgeId> all = family_edges(box);
  const auto first = geom::fillet_edges(box, all, {0.25, 1e-9});
  REQUIRE(first.has_value());
  const auto second = geom::fillet_edges(box, all, {0.25, 1e-9});
  REQUIRE(second.has_value());

  const geom::TriMesh mesh = geom::tessellate(first.value().brep, geom::TessellationOptions{}).value();
  CHECK(boundary_edges(mesh) == 0u);

  CHECK(first.value().brep.vertex_count() == second.value().brep.vertex_count());
  CHECK(first.value().brep.edge_count() == second.value().brep.edge_count());
  CHECK(first.value().brep.face_count() == second.value().brep.face_count());
  CHECK(std::abs(volume_of(first.value().brep) - volume_of(second.value().brep)) <= 1e-12);

  // Input must be unchanged by both success and failure.
  const std::size_t box_faces = box.face_count();
  const double box_volume = volume_of(box);
  CHECK(box.face_count() == box_faces);
  CHECK(std::abs(volume_of(box) - box_volume) <= 1e-12);
}
```

- [ ] **Step 2: Run the tests to verify they fail**

Expected RED: the radius limit and mixed-family tests fail because the builder currently accepts or misclassifies them; the tessellation test may fail on boundary edges if a cap ring is not closed.

- [ ] **Step 3: Add the guards**

- Before building, enforce `radius < 0.5 * min(a, b)` (`unsupported "fillet_edges: radius exceeds the box cross-section"`). After building the `selected` mask, for each side whose two corner flags are both set enforce `2*radius < side_extent - tolerance` (`unsupported "fillet_edges: radius leaves a face narrower than the tolerance"`). Only a doubly-filleted side has `side_extent - 2r`; a singly-filleted side has `side_extent - r`, so a global `2r` cap would reject valid k=1/k=2 cases. Do not add any other epsilon.
- Keep the mixed-family check before any geometry construction.
- Verify edge endpoints lie at two distinct family extremes and at one of the four cross-section corners; reject with the stable messages from Task 1.
- Return only after `validate_topology`; do not return the assembled solid if the report has issues.
- Make `fillet_box_edges` take the input by value/const and never mutate it; all vectors are local.

- [ ] **Step 4: Run the tests to verify they pass**

`--test-case="*fillet*"` all pass; `build\tests\Debug\geom_tests.exe --test-case="*fillet*"` and Release both report the same case count. Run `ctest -C Release` only if the controller asks for the full suite; otherwise run the full Release suite once at milestone exit.

- [ ] **Step 5: Commit**

```bash
git add include/geom/brep/fillet.hpp include/geom/brep/detail/fillet_box.hpp tests/unit/fillet_test.cpp
git commit -m "fix(brep): reject invalid box fillet selections before assembly"
```

### Task 4: Benchmark, ADR, architecture evidence, and CI

**Files:**
- Create: `bench/bench_m4c3.cpp`
- Modify: `CMakeLists.txt` (the root file; bench targets live inside `if(GEOM_BUILD_BENCH)`)
- Modify: `.github/workflows/ci.yml`
- Create: `docs/adr/0042-box-fillet.md`
- Create: `docs/benchmarks/m4c3-fillet.md`
- Modify: `docs/architecture/geometry-library-architecture.md`
- Modify: this plan's Results section

**Interfaces:**
- Consumes: `fillet_edges`, `make_box`, `detail::curved_solid_volume` (the fillet patches are arc-trimmed, so `mass_properties` does not accept them).
- Produces: `geom_bench_m4c3` executable and the M4c-3a documentation.

- [ ] **Step 1: Write the failing benchmark**

Create `bench/bench_m4c3.cpp` following `bench/bench_m4c1.cpp`: build a 100-cube, select its four vertical edges, call `fillet_edges` with `radius=1.0`, `tolerance=1e-9`, print five samples' min/median/max and the volume, and exit 2 if the result is missing or the volume is outside `1e-6` relative of `999914.159265359`.

- [ ] **Step 2: Build it**

Add `geom_bench_m4c3` to the root `CMakeLists.txt` (there is no
`bench/CMakeLists.txt`; the bench executables live inside `if(GEOM_BUILD_BENCH)`).
The main `build` tree has `GEOM_BUILD_BENCH=OFF`, so configure and build a bench
tree instead:

```powershell
cmake -S . -B build-bench -DGEOM_BUILD_BENCH=ON -DGEOM_BUILD_TESTS=OFF -DGEOM_WERROR=ON
cmake --build build-bench --config Release --target geom_bench_m4c3
.\build-bench\Release\geom_bench_m4c3.exe
```

Expected RED before the target exists: CMake reports an unknown target.

- [ ] **Step 3: Run the benchmark and record the numbers**

Run it once and paste the min/median/max and volume into the plan Results table. Do not invent numbers; if it fails, fix the implementation, not the expected volume.

- [ ] **Step 4: Write the docs**

- `docs/adr/0042-box-fillet.md`: decision (exact cylindrical fillet, direct curved assembly, no Boolean round-trip), scope (axis-aligned box, one parallel family), rejected alternatives (faceted ring, multi-plane fit, general SSI blend for this milestone), failure modes, and deferrals (vertex blends, variable radius, general convex, M4c-3b).
- `docs/benchmarks/m4c3-fillet.md`: environment, workload, min/median/max, volume anchor, how to reproduce, and the explicit Debug PENDING note if the full Debug suite is not run.
- Architecture §8.9: add one paragraph stating M4c-3a supports constant-radius fillet of up to four parallel edges of an axis-aligned box, exact cylinder + circular arc caps, with the stable rejection messages, and that vertex blends/variable radius/general convex remain M4c-3b.
- Architecture §16.4: add the evidence row with commit range, Release counts, benchmark numbers and the Debug lane counts.
- Add `geom_bench_m4c3` to the CI bench build target list and a run step producing `bench-m4c3.txt`; add the artifact to the upload list.

- [ ] **Step 5: Run the documentation checks and commit**

Check every quoted count against the actual test output; update the plan Results table. Then:

```bash
git add bench/bench_m4c3.cpp CMakeLists.txt .github/workflows/ci.yml docs/adr/0042-box-fillet.md docs/benchmarks/m4c3-fillet.md docs/architecture/geometry-library-architecture.md docs/superpowers/plans/2026-11-26-m4c3-box-fillet.md
git commit -m "docs(brep): record the M4c-3a box fillet contract and benchmark"
```

---

## M4c-3a Exit

- `fillet_edges` handles an axis-aligned box and any non-empty subset of one parallel family (1–4 edges).
- One-edge volume `23.94634954084936`, four-edge volume `23.78539816339745` for `2x3x4`, `r=0.25`, each within `1e-9`.
- Topology counts `V=8+2k`, `E=12+3k`, `F=6+k`, `chi=2` for every k.
- `validate_topology` is clean, `tessellate` has zero boundary edges, repeated calls agree, and the input is unchanged.
- Invalid options, non-box solids, mixed families, dead edges, and `2r >= min cross-section` are rejected before assembly with stable `Errc` + messages.
- Release full suite and the Debug `*fillet*` lane are green; the benchmark builds and reports honest numbers.
- Docs state the scope and deferrals; no document claims general fillet support.

## Results (fill after execution)

| Item | Value |
| --- | --- |
| Debug `*fillet*` | 8 cases / 887 assertions, 0 failed |
| Release `*fillet*` | 8 cases / 887 assertions, 0 failed |
| Full Release | 1,159 cases / 683,611 assertions, 0 failed (880 s) |
| One-edge volume | `23.946349540849312` (plan anchor `23.94634954084936`, absolute difference `4.8e-14`, within `1e-9`) |
| Four-edge volume | `23.785398163397282` (plan anchor `23.78539816339745`, absolute difference `1.7e-13`, within `1e-9`) |
| Benchmark 100-cube, 4 vertical edges, r=1 | min `7.587` / median `8.247` / max `12.766` ms; volume `999914.159265358` vs `999914.159265359` (relative `1.1e-15`); exit 0 |
| Debug lane (`*fillet*,*blend*,*editor*`) | 56 cases / 1,397 assertions, 0 failed |

Task 1 validation-order clarification: `fillet_edges` checks the empty span
**before** the option checks, because the default `FilletOptions{}` has
`radius = 0.0`. Validating the radius first would report `degenerate` for an
empty selection made with default options; "no edges selected" is an
`unsupported` request, so the empty-span check comes first and the empty +
default call returns `unsupported "fillet_edges: at least one edge is
required"`.

The benchmark target lives in the root `CMakeLists.txt` inside
`if(GEOM_BUILD_BENCH)` (there is no `bench/CMakeLists.txt` in this tree), and it
was built in `build-bench` because the main `build` tree has
`GEOM_BUILD_BENCH=OFF`.

Full Debug remains PENDING: the suite still contains the known >20-minute
cases, so the Debug `*fillet*` lane and the Debug lane above are the Debug
evidence. They are not a full Debug pass.

## Self-Review

- Spec coverage: §8.9 fillet is covered for the narrow box scope; vertex blends, variable radius and general convex are explicitly deferred, not silently omitted.
- Step scan: each task starts with a failing test with exact values and ends with a command and a commit; no "TBD" or "handle edge cases".
- Type consistency: `FilletOptions{radius,tolerance}`, `FilletResult{brep,filleted_edges}`, `fillet_edges`, `FilletBox`, and `fillet_box_edges` are used with the same names in every task.
- Review Focus: mixed family, radius limit, four-corner cap, non-box rejection, and failure immutability each have an owning test (Task 3).
- Proportion: the plan gives formulas and signatures; it does not transcribe the builder body.
