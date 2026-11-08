# M3a: Planar B-rep Boolean and Multi-Shell Solids

## Goal

Turn two valid, closed, all-planar B-rep solids into a third valid, closed solid
(or a diagnosed failure), for union, intersection, difference and symmetric
difference, including coplanar contact and results with cavities, through-holes
and multiple disconnected bodies. This is the M3 slice that proves the §8.8
pipeline — face pairing, face splitting, classification, stitching, validation —
end to end without needing surface-surface intersection. Curved faces are
rejected with `Errc::unsupported`; SSI, curved booleans and QA-9 are M3b.

## Architecture

The boolean is a **cut-and-classify** pipeline over planar faces:

1. **Pair faces** of A and B by bounding-box overlap. For each pair, compute the
   intersection of the two planar face regions: a set of straight 3D segments
   (`intersect_planar_faces`). Coplanar pairs contribute the boundary of their
   2D overlap polygon.
2. **Cut every face** by the segments that lie in its plane
   (`cut_face_by_segments`): project the face loops and the cuts into a stable
   plane frame, build a planar arrangement with the existing exact
   `detail::split_segments` + `Arrangement2`, and keep the arrangement faces
   that lie inside the original face. Each kept face is a fragment: one outer
   ring plus zero or more hole rings, in 2D.
3. **Classify every fragment** by sampling the point at `sample + n*delta` and
   `sample - n*delta`, testing membership in A and B with
   `point_in_solid` (ray parity against the planar face loops), and evaluating
   the boolean membership function on both sides. The fragment survives iff the
   two sides differ; its normal is flipped when the material is on the `+n`
   side. This one rule handles coplanar duplicate faces, interior faces that
   must disappear, and orientation, without special-casing the operation.
4. **Assemble** the surviving fragments (`assemble_solid`): weld vertices on a
   `tol` grid, build edges from ring adjacency, require every edge to be used
   exactly twice with opposite directions, group faces into shells by signed
   volume, assign negative shells as cavities of the containing positive shell,
   and emit one solid per positive component.
5. **Validate** the result with `validate_topology` and `validate_tolerances`
   before returning it.

All predicates that decide topology are exact or derived from exact results
(`sign_orient2d`, `Expansion`, the arrangement). Floating-point arithmetic is
only used to produce candidate geometry and sample points.

## Tech Stack

- C++20, MSVC 19.44 /W4 /WX /utf-8 (GCC 13 and Clang 18 in CI)
- CMake 4.4.4, doctest v2.4.11
- Existing layers: `geom/core`, `geom/primitives`, `geom/surfaces/plane.hpp`,
  `geom/cgeom` (`boolean2`, `arrangement2`), `geom/brep`
  (`topology`, `solid_editor`, `validate`, `mass_properties`, `tessellate`)

## Spec

- `docs/architecture/geometry-library-architecture.md` v2.1 §5.8 (solid =
  outer shell + inner shells), §8.2, §8.3, §8.4, §8.7, §8.8, §12.3,
  §16.2/§16.3 (M3 scope).
- `docs/quality/testing-strategy.md` (coverage gates, exact/float dual paths,
  negative tests for every `Errc`).
- QA-8: boolean of two ~10k-face solids, volume relative error ≤ 1e-9,
  ≤ 5 s single-threaded. This plan measures QA-8 and records PASS/FAIL
  honestly; it does not weaken the gate.

## Global Constraints

- Header-only, `inline namespace v1`, `#pragma once`? Follow the existing
  headers (they use include guards in `geom/core`, `#pragma once` elsewhere;
  match the neighbouring file in each directory).
- No new third-party dependency. No `throw` in library code; errors are
  `Result<T>` with `Errc`.
- No global epsilon. Every tolerance is an explicit parameter; all new public
  functions take `double tol` or a `BooleanOptions` with `tolerance`.
- Every new `Errc` value used here already exists (`invalid_topology`,
  `not_closed`, `unsupported`, `degenerate`, `numeric_failure`,
  `tolerance_violation`, `index_out_of_range`, `non_finite`). If a new value is
  added, update `to_string` and `tests/unit/coverage_gaps_test.cpp` in the same
  commit.
- Every test file must be added to `GEOM_TEST_SOURCES` in
  `tests/CMakeLists.txt`.
- TDD: write the test, run it, show the RED failure, implement, run Debug and
  Release, then commit once with a `feat(m3a): ...` message.
- No weak assertions: `CHECK(x.has_value())` alone is not a test; assert the
  concrete value, error code and message.
- Determinism: no iteration over `std::unordered_map` without sorting the
  result; welding and assembly must be stable across runs. Tests run twice and
  compare.
- Public headers must stay under 400 lines where practical; split
  `boolean_solid.hpp` into `inside.hpp`, `planar_cut.hpp`, `assemble.hpp` as
  planned below.

## Scope and Deviations (for review)

1. **Planar faces only.** Any face whose surface variant is not `Plane3`
   returns `Errc::unsupported` with message
   `"boolean_solid: curved faces are deferred to M3b"`. No SSI, no QA-9.
2. **`Pcurve2d` stays a straight 2D segment.** Planar boolean only produces
   straight edges, so a line pcurve is exact. Curved pcurves arrive with M3b.
3. **No general healing / gap bridging.** Inputs must already validate; a
   `1e-6` gap is not healed (QA-12 deferred). `assemble_solid` fails with
   `not_closed` when an edge is used once.
4. **Coplanar support is region-level.** Two coplanar faces are handled by
   their 2D polygon overlap; partially overlapping coplanar faces, opposite
   normals, and stacked boxes are supported. Exact coincident vertices are
   welded; a coplanar contact that leaves a non-manifold edge is a diagnosed
   failure, not a repair.
5. **Multi-shell results are supported; multi-shell inputs are not.** A solid
   with cavities is a valid result, but `boolean_solid` requires both operands
   to have exactly one shell each (message
   `"boolean_solid: multi-shell operands are deferred to M3c"`).
6. **No face merging.** Coplanar fragments stay split (a stacked-box union has
   10 faces, not 6). Volume and validity are the contract; face count is not.

## Review Focus

Review these five areas in that order:

1. **Classification rule.** Verify the `f(p+) != f(p-)` rule and the flip
   direction on four cases: transversal union, stacked coplanar union,
   difference that opens a cavity, and intersection of touching boxes
   (empty result).
2. **Exactness boundary.** The segment computation and the arrangement use
   exact predicates; sample points and welding use doubles. Check that no
   topology decision (inside/outside, edge pairing) depends on a float
   comparison with an ad-hoc epsilon.
3. **Stitching.** Every edge must be used exactly twice with opposite
   directions. Check the welding tolerance, the deterministic ordering, and
   that the shell signed-volume rule cannot misclassify a cavity.
4. **Failure modes.** Non-planar face, multi-shell operand, non-manifold
   contact, empty result and open input must return a specific `Errc` and
   message; no crash, no silent empty `Brep`.
5. **QA-8 honesty.** The bench must print the measured time and volume error;
   a FAIL is recorded as FAIL in `docs/benchmarks/m3a.md`, never rounded into
   a pass.

## File Map

New public headers under `include/geom/brep/`:

- `plane_frame.hpp` — stable per-plane (u, v) frame
- `face_query.hpp` — point-in-face and point-in-shell queries
- `face_intersect.hpp` — planar face-pair intersection segments
- `inside.hpp` — point-in-solid ray parity
- `planar_cut.hpp` — face cutting into fragments
- `assemble.hpp` — fragments to validated `Brep`
- `boolean_solid.hpp` — the driver, ops, options, report

New tests under `tests/unit/`:

- `plane_frame_test.cpp`
- `brep_face_query_test.cpp`
- `brep_face_intersect_test.cpp`
- `brep_inside_test.cpp`
- `brep_planar_cut_test.cpp`
- `brep_assemble_test.cpp`
- `brep_boolean_solid_test.cpp`
- `brep_boolean_shell_test.cpp`
- `brep_boolean_solid_errors_test.cpp`

Bench and docs:

- `bench/bench_m3.cpp`
- `docs/benchmarks/m3a.md`
- `docs/adr/0023-planar-boolean-first.md`
- architecture §16.3 and this plan's checkboxes

## Task 1: Stable plane frame

**Files:**
- Create: `include/geom/brep/plane_frame.hpp`
- Create: `tests/unit/plane_frame_test.cpp`
- Modify: `tests/CMakeLists.txt`

**Interfaces:**
- Consumes: `Plane3`, `Point3d`, `Vec3d`, `Point2d`, `dot`, `cross`, `Result`.
- Produces:
  - `struct PlaneFrame { Point3d origin; Vec3d u; Vec3d v; Vec3d normal; };`
  - `Result<PlaneFrame> make_plane_frame(const Plane3& plane);`
  - `Point2d plane_frame_to_2d(const PlaneFrame& frame, const Point3d& p) noexcept;`
  - `Point3d plane_frame_to_3d(const PlaneFrame& frame, const Point2d& q) noexcept;`
  - `Plane3 plane_frame_plane(const PlaneFrame& frame) noexcept;`
  - `void detail::plane_frame_basis(const Vec3d& normal, Vec3d& u, Vec3d& v);`

**Step 1: Write the failing test**

```cpp
#include <doctest/doctest.h>
#include <geom/brep/plane_frame.hpp>

using namespace geom;

TEST_CASE("plane_frame: orthonormal right-handed frame on a tilted plane") {
  const Vec3d n{1.0, 2.0, 2.0};  // |n| = 3
  const Plane3 plane{n, 7.0};
  const auto frame = make_plane_frame(plane);
  REQUIRE(frame.has_value());
  const PlaneFrame& f = frame.value();
  CHECK(dot(f.normal, f.normal) == doctest::Approx(1.0).epsilon(1e-15));
  CHECK(dot(f.u, f.u) == doctest::Approx(1.0).epsilon(1e-15));
  CHECK(dot(f.v, f.v) == doctest::Approx(1.0).epsilon(1e-15));
  CHECK(dot(f.u, f.v) == doctest::Approx(0.0).epsilon(1e-15));
  CHECK(dot(cross(f.u, f.v), f.normal) == doctest::Approx(1.0).epsilon(1e-15));
  CHECK(f.origin[0] + 2.0 * f.origin[1] + 2.0 * f.origin[2] ==
        doctest::Approx(7.0).epsilon(1e-15));
}

TEST_CASE("plane_frame: 3d to 2d round trip is exact on the plane") {
  const Plane3 plane{Vec3d{0.0, 1.0, 0.0}, -2.0};  // y = -2
  const auto frame = make_plane_frame(plane);
  REQUIRE(frame.has_value());
  for (const Point3d p : {Point3d{1.0, -2.0, 0.5}, Point3d{-3.25, -2.0, 4.0},
                          Point3d{0.0, -2.0, 0.0}}) {
    const Point3d back = plane_frame_to_3d(frame.value(), plane_frame_to_2d(frame.value(), p));
    CHECK(back[0] == doctest::Approx(p[0]).epsilon(1e-14));
    CHECK(back[1] == doctest::Approx(p[1]).epsilon(1e-14));
    CHECK(back[2] == doctest::Approx(p[2]).epsilon(1e-14));
  }
}

TEST_CASE("plane_frame: degenerate and non-finite planes are rejected") {
  const auto zero = make_plane_frame(Plane3{Vec3d{0.0, 0.0, 0.0}, 0.0});
  REQUIRE_FALSE(zero.has_value());
  CHECK(zero.error().code == Errc::degenerate);
  const auto nan = make_plane_frame(Plane3{Vec3d{0.0, 0.0, std::numeric_limits<double>::quiet_NaN()}, 0.0});
  REQUIRE_FALSE(nan.has_value());
  CHECK(nan.error().code == Errc::non_finite);
}

TEST_CASE("plane_frame: the basis only depends on the normal direction") {
  const Vec3d n{3.0, 4.0, 0.0};
  const auto a = make_plane_frame(Plane3{n, 1.0});
  const auto b = make_plane_frame(Plane3{n, 9.0});
  REQUIRE(a.has_value());
  REQUIRE(b.has_value());
  CHECK(dot(a.value().u, b.value().u) == doctest::Approx(1.0).epsilon(1e-15));
  CHECK(dot(a.value().v, b.value().v) == doctest::Approx(1.0).epsilon(1e-15));
  CHECK(plane_frame_plane(a.value()).offset == doctest::Approx(1.0).epsilon(1e-15));
}
```

**Step 2: Run it**

```powershell
cmake --build build --config Debug --target geom_tests
build\tests\Debug\geom_tests.exe --test-case="plane_frame*"
```

Expected RED: `plane_frame_test.cpp(3): error C1083: Cannot open include file:
'geom/brep/plane_frame.hpp'`.

**Step 3: Implement**

- Validate `normal` components with `std::isfinite` (→ `non_finite`) and
  `|normal| > 0` (→ `degenerate`).
- Unit normal `n = normal / |normal|`. Choose the reference axis with the
  largest `|n[i]|`; `u = normalize(cross(ref, n))`; `v = cross(n, u)`. Then
  `u × v = n`.
- `origin = n * (offset / |normal|)`.
- `to_2d(p) = {dot(p - origin, u), dot(p - origin, v)}`;
  `to_3d(q) = origin + q[0]*u + q[1]*v`.
- `plane_frame_plane` returns `{normal, dot(normal, origin)}`.

**Step 4: Run**

Debug and Release full suite must pass; standalone-compile the new header with
`/W4 /WX /utf-8`.

**Step 5: Commit**

`feat(m3a): add stable plane frames for planar B-rep work`

## Task 2: Point-in-face and point-in-shell queries

**Files:**
- Create: `include/geom/brep/face_query.hpp`
- Create: `tests/unit/brep_face_query_test.cpp`
- Modify: `tests/CMakeLists.txt`

**Interfaces:**
- Consumes: `PlaneFrame` (Task 1), `Brep`, `FaceId`, `ShellId`, `boolean_point_in_ring`,
  `make_box`.
- Produces:
  - `Result<std::vector<std::vector<Point2d>>> face_loops_2d(const Brep&, FaceId, const PlaneFrame&);`
  - `Result<bool> point_in_face(const Brep&, FaceId, const Point3d&, double tol);`
  - `Result<bool> point_on_face_boundary(const Brep&, FaceId, const Point3d&, double tol);`
  - `Result<bool> point_in_shell(const Brep&, ShellId, const Point3d&, double tol);`
  - `Result<Plane3> face_plane(const Brep&, FaceId);`

**Step 1: Write the failing test**

```cpp
#include <doctest/doctest.h>
#include <geom/brep/face_query.hpp>
#include <geom/brep/primitives.hpp>

using namespace geom;

namespace {
// Finds the face whose plane has the given outward unit normal and offset.
// Copy this helper into every test file that needs to name a face.
FaceId plane_face(const Brep& brep, Vec3d want_normal, double want_offset) {
  for (const FaceId id : brep.face_ids()) {
    const auto& variant = brep.geometry().surfaces[brep.face(id).surface];
    const auto* plane = std::get_if<Plane3>(&variant);
    if (plane == nullptr) continue;
    const double len = std::sqrt(dot(plane->normal, plane->normal));
    const Vec3d unit{plane->normal[0] / len, plane->normal[1] / len,
                     plane->normal[2] / len};
    if (dot(unit, want_normal) > 0.999 &&
        std::fabs(plane->offset / len - want_offset) <= 1e-9) {
      return id;
    }
  }
  return FaceId::null();
}
}

TEST_CASE("face_query: a box top face contains its centre, not its rim") {
  const Brep box = make_box(Point3d{0.0, 0.0, 0.0}, Point3d{2.0, 2.0, 2.0}, 1e-12).value();
  const FaceId top = plane_face(box, Vec3d{0.0, 0.0, 1.0}, 2.0);
  REQUIRE_FALSE(top.is_null());
  REQUIRE(point_in_face(box, top, Point3d{1.0, 1.0, 2.0}, 1e-9).has_value());
  CHECK(point_in_face(box, top, Point3d{1.0, 1.0, 2.0}, 1e-9).value());
  CHECK_FALSE(point_in_face(box, top, Point3d{0.0, 0.0, 2.0}, 1e-9).value());
  CHECK_FALSE(point_in_face(box, top, Point3d{3.0, 1.0, 2.0}, 1e-9).value());
  CHECK(point_on_face_boundary(box, top, Point3d{0.0, 1.0, 2.0}, 1e-9).value());
}

TEST_CASE("face_query: a shell contains an interior point and not an exterior one") {
  const Brep box = make_box(Point3d{0.0, 0.0, 0.0}, Point3d{2.0, 2.0, 2.0}, 1e-12).value();
  const ShellId shell = box.shell_ids().front();
  CHECK(point_in_shell(box, shell, Point3d{1.0, 1.0, 1.0}, 1e-9).value());
  CHECK_FALSE(point_in_shell(box, shell, Point3d{5.0, 1.0, 1.0}, 1e-9).value());
}
```

**Step 2: Run it** — expected RED: missing `face_query.hpp`.

**Step 3: Implement**

- `face_plane`: map `face.surface` through `geometry().surfaces`; return
  `unsupported` for a non-`Plane3` variant.
- `face_loops_2d`: walk `face.first_loop`, then `loop.next`; for each loop walk
  `coedge.next` from `loop.start`; append the start vertex of each coedge
  (`edge.v0` or `edge.v1` when `coedge.reversed`), projected with the frame.
  Drop a repeated closing point.
- `point_in_face`: project `p`; `point_in_face` is true when the point is
  inside the outer loop and outside every hole. Outer loop = the loop with
  `outer == true`; use `boolean_point_in_ring` for each ring. Boundary is
  decided separately by `point_on_face_boundary` (distance to every boundary
  segment in 2D ≤ `tol`; returns false when the point is not on the plane
  within `tol`).
- `point_in_shell`: ray parity over the shell's faces, using a direction chosen
  from a fixed table `{(0.5773, 0.5774, 0.5775), ...}`; for each face intersect
  the ray with its plane, reject `t <= tol`, test `point_in_face`. If any hit is
  within `tol` of a face boundary, retry the next direction; after the table is
  exhausted return `numeric_failure`. Odd parity = inside.

**Step 4: Run** Debug/Release plus standalone header compile.

**Step 5: Commit**

`feat(m3a): add exact point-in-face and point-in-shell queries`

## Task 3: Planar face-pair intersection

**Files:**
- Create: `include/geom/brep/face_intersect.hpp`
- Create: `tests/unit/brep_face_intersect_test.cpp`
- Modify: `tests/CMakeLists.txt`

**Interfaces:**
- Consumes: Tasks 1–2, `boolean2`, `Polygon2d`, `Segment3`.
- Produces:
  - `struct Segment3 { Point3d a; Point3d b; };`
  - `struct FacePairCut { std::vector<Segment3> segments; bool coplanar = false; };`
  - `Result<FacePairCut> intersect_planar_faces(const Brep& a, FaceId fa, const Brep& b, FaceId fb, double tol);`

**Step 1: Write the failing test**

```cpp
#include <doctest/doctest.h>
#include <geom/brep/face_intersect.hpp>
#include <geom/brep/primitives.hpp>

using namespace geom;

TEST_CASE("face_intersect: perpendicular boxes cut a top face along one segment") {
  const Brep a = make_box(Point3d{0, 0, 0}, Point3d{2, 2, 2}, 1e-12).value();
  const Brep b = make_box(Point3d{1, 1, 0}, Point3d{3, 3, 3}, 1e-12).value();
  const FaceId top = plane_face(a, Vec3d{0.0, 0.0, 1.0}, 2.0);
  const FaceId side = plane_face(b, Vec3d{-1.0, 0.0, 0.0}, 1.0);
  REQUIRE_FALSE(top.is_null());
  REQUIRE_FALSE(side.is_null());
  const auto cut = intersect_planar_faces(a, top, b, side, 1e-12);
  REQUIRE(cut.has_value());
  CHECK_FALSE(cut.value().coplanar);
  REQUIRE(cut.value().segments.size() == 1);
  CHECK(cut.value().segments[0].a[0] == doctest::Approx(1.0).epsilon(1e-12));
  CHECK(cut.value().segments[0].a[2] == doctest::Approx(2.0).epsilon(1e-12));
  CHECK(cut.value().segments[0].b[1] == doctest::Approx(2.0).epsilon(1e-12));
}

TEST_CASE("face_intersect: coplanar faces return the overlap boundary") {
  const Brep a = make_box(Point3d{0, 0, 0}, Point3d{2, 2, 2}, 1e-12).value();
  const Brep b = make_box(Point3d{0, 0, 2}, Point3d{1, 1, 4}, 1e-12).value();
  const FaceId top_a = plane_face(a, Vec3d{0.0, 0.0, 1.0}, 2.0);
  const FaceId bottom_b = plane_face(b, Vec3d{0.0, 0.0, -1.0}, 2.0);
  REQUIRE_FALSE(top_a.is_null());
  REQUIRE_FALSE(bottom_b.is_null());
  const auto cut = intersect_planar_faces(a, top_a, b, bottom_b, 1e-12);
  REQUIRE(cut.has_value());
  CHECK(cut.value().coplanar);
  CHECK(cut.value().segments.size() == 4);
}

TEST_CASE("face_intersect: disjoint faces produce no segments") {
  const Brep a = make_box(Point3d{0, 0, 0}, Point3d{1, 1, 1}, 1e-12).value();
  const Brep b = make_box(Point3d{5, 5, 5}, Point3d{6, 6, 6}, 1e-12).value();
  const auto cut = intersect_planar_faces(
      a, plane_face(a, Vec3d{0, 0, 1}, 1.0), b, plane_face(b, Vec3d{0, 0, -1}, 5.0), 1e-12);
  REQUIRE(cut.has_value());
  CHECK(cut.value().segments.empty());
  CHECK_FALSE(cut.value().coplanar);
}
```

**Step 2: Run it** — expected RED: missing header.

**Step 3: Implement**

- Get both planes and frames. Compute `d = cross(nA, nB)`; `|d| <= 1e-12`
  means parallel: if the offset difference along a unit normal is ≤ `tol`,
  handle the coplanar case, else return empty.
- Non-parallel: line point `p0` from solving `nA·x = bA`, `nB·x = bB` with the
  third equation `d·x = 0`; direction `d`.
- Clip the line to face A: project the line into A's frame; for every boundary
  edge compute the intersection parameter with the infinite line; sort the
  parameters; keep intervals whose midpoint is inside the face
  (`point_in_face`). Do the same for B, then intersect the two interval sets
  (map every endpoint to the global line parameter `t = dot(x - p0, d)`).
- Emit one `Segment3` per non-empty interval intersection.
- Coplanar: project B's face loops into A's frame as a `Polygon2d`, run
  `boolean2({faceA, faceB}, BooleanOp::intersection)`, and emit each result
  ring edge as a 3D segment (`plane_frame_to_3d`). Set `coplanar = true`.

**Step 4: Run** Debug/Release.

**Step 5: Commit**

`feat(m3a): intersect planar B-rep faces into cut segments`

## Task 4: Point-in-solid ray parity

**Files:**
- Create: `include/geom/brep/inside.hpp`
- Create: `tests/unit/brep_inside_test.cpp`
- Modify: `tests/CMakeLists.txt`

**Interfaces:**
- Consumes: Task 2 `point_in_face`, `face_plane`.
- Produces:
  - `Result<bool> point_in_solid(const Brep&, SolidId, const Point3d&, double tol);`

**Step 1: Write the failing test**

```cpp
TEST_CASE("inside: a box classifies interior and exterior points") {
  const Brep box = make_box(Point3d{0, 0, 0}, Point3d{2, 2, 2}, 1e-12).value();
  const SolidId solid = box.solid_ids().front();
  CHECK(point_in_solid(box, solid, Point3d{1, 1, 1}, 1e-9).value());
  CHECK(point_in_solid(box, solid, Point3d{0.1, 0.1, 0.1}, 1e-9).value());
  CHECK_FALSE(point_in_solid(box, solid, Point3d{2.5, 1, 1}, 1e-9).value());
  CHECK_FALSE(point_in_solid(box, solid, Point3d{-1, -1, -1}, 1e-9).value());
}
```

Add a translated box (offset 1e4) and a ray that grazes an edge to exercise
the direction retry.

**Step 2: Run it** — expected RED.

**Step 3: Implement**

- For each direction in the retry table, count intersections over every face of
  the solid's shells; parity decides. Use `face_plane` + ray-plane `t`; ignore
  `t <= tol`; test `point_in_face`; if the hit is within `tol` of the face
  boundary, mark the direction ambiguous and retry.
- If all directions are ambiguous, return
  `Error{Errc::numeric_failure, "point_in_solid: every probe ray is degenerate", ...}`.
- An empty solid list or a solid with no shells → `invalid_topology`.

**Step 4: Run** Debug/Release.

**Step 5: Commit**

`feat(m3a): classify points against planar solids by ray parity`

## Task 5: Cut a planar face into fragments

**Files:**
- Create: `include/geom/brep/planar_cut.hpp`
- Create: `tests/unit/brep_planar_cut_test.cpp`
- Modify: `tests/CMakeLists.txt`

**Interfaces:**
- Consumes: `detail::split_segments`, `Arrangement2`,
  `triangulate_polygon2`, `Polygon2d`, Tasks 1–3.
- Produces:
  - `struct FaceFragment2d { std::vector<Point2d> outer; std::vector<std::vector<Point2d>> holes; Point2d sample; };`
  - `Result<std::vector<FaceFragment2d>> cut_face_by_segments(const Brep&, FaceId, std::span<const Segment3>, double tol);`

**Step 1: Write the failing test**

```cpp
namespace {
double shoelace(std::span<const Point2d> ring) {
  double sum = 0.0;
  for (std::size_t i = 0; i < ring.size(); ++i) {
    const Point2d& a = ring[i];
    const Point2d& b = ring[(i + 1) % ring.size()];
    sum += a[0] * b[1] - b[0] * a[1];
  }
  return 0.5 * sum;
}
double fragment_area(const FaceFragment2d& f) {
  double area = std::fabs(shoelace(f.outer));
  for (const auto& hole : f.holes) area -= std::fabs(shoelace(hole));
  return area;
}
using Segments = std::vector<Segment3>;
}

TEST_CASE("planar_cut: one transverse cut splits a square into two fragments") {
  const Brep box = make_box(Point3d{0, 0, 0}, Point3d{2, 2, 2}, 1e-12).value();
  const FaceId top = plane_face(box, Vec3d{0, 0, 1}, 2.0);
  const Segments cuts{{Point3d{1, 0, 2}, Point3d{1, 2, 2}}};
  const auto fragments = cut_face_by_segments(box, top, cuts, 1e-12);
  REQUIRE(fragments.has_value());
  REQUIRE(fragments.value().size() == 2);
  double total = 0.0;
  for (const FaceFragment2d& f : fragments.value()) {
    CHECK(f.holes.empty());
    total += fragment_area(f);
  }
  CHECK(total == doctest::Approx(4.0).epsilon(1e-12));
}

TEST_CASE("planar_cut: no cuts returns the whole face") {
  const Brep box = make_box(Point3d{0, 0, 0}, Point3d{2, 2, 2}, 1e-12).value();
  const FaceId top = plane_face(box, Vec3d{0, 0, 1}, 2.0);
  const Segments none;
  const auto fragments = cut_face_by_segments(box, top, none, 1e-12);
  REQUIRE(fragments.has_value());
  REQUIRE(fragments.value().size() == 1);
  CHECK(fragments.value()[0].outer.size() == 4);
  CHECK(fragment_area(fragments.value()[0]) == doctest::Approx(4.0).epsilon(1e-12));
}

TEST_CASE("planar_cut: a cut that misses the face returns the whole face") {
  const Brep box = make_box(Point3d{0, 0, 0}, Point3d{2, 2, 2}, 1e-12).value();
  const FaceId top = plane_face(box, Vec3d{0, 0, 1}, 2.0);
  const Segments outside{{Point3d{5, -1, 2}, Point3d{5, 3, 2}}};
  const auto fragments = cut_face_by_segments(box, top, outside, 1e-12);
  REQUIRE(fragments.has_value());
  REQUIRE(fragments.value().size() == 1);
  CHECK(fragment_area(fragments.value()[0]) == doctest::Approx(4.0).epsilon(1e-12));
}

TEST_CASE("planar_cut: crossing cuts make four fragments whose areas sum to four") {
  const Brep box = make_box(Point3d{0, 0, 0}, Point3d{2, 2, 2}, 1e-12).value();
  const FaceId top = plane_face(box, Vec3d{0, 0, 1}, 2.0);
  const Segments cuts{{Point3d{1, 0, 2}, Point3d{1, 2, 2}},
                      {Point3d{0, 1, 2}, Point3d{2, 1, 2}}};
  const auto fragments = cut_face_by_segments(box, top, cuts, 1e-12);
  REQUIRE(fragments.has_value());
  REQUIRE(fragments.value().size() == 4);
  double total = 0.0;
  for (const FaceFragment2d& f : fragments.value()) total += fragment_area(f);
  CHECK(total == doctest::Approx(4.0).epsilon(1e-12));
}
```

**Step 2: Run it** — expected RED.

**Step 3: Implement**

- If `segments` is empty, return one fragment built from the face loops
  (projected and normalised: outer CCW, holes CW); this fast path keeps the
  boolean cheap for untouched faces.
- Project every loop edge and every cut segment to the frame.
- Build `SplitInput`s with `source = 0` for loop edges and `1..n` for cuts;
  call `detail::split_segments`; build `Arrangement2`.
- For each bounded arrangement face, collect **every** boundary cycle whose
  half-edges carry that face id (a face may have several cycles: one outer and
  several hole rings; do not start from a single half-edge and assume one
  cycle). Convert each cycle to 2D; the largest-|area| cycle is the outer ring,
  the rest are holes; normalise orientations.
- Individual cut segments may end inside the face when the partner face ends
  there; the segments from all partner faces form the closed cut curve
  together, so the arrangement is expected to be complete. If a face still has
  a dangling cut, it is a diagnosed `invalid_topology` error rather than a
  silent spike in a ring.
- Build a `Polygon2d` from the rings and call `triangulate_polygon2`; the
  centroid of the first triangle is a strictly interior sample (if
  triangulation fails, try the next triangle; if none, skip the face).
- Keep the arrangement face only when the sample is inside the original face
  (`point_in_face` semantics on the projected loops): inside outer and outside
  every hole.
- Store that interior sample in `FaceFragment2d::sample`.
- Return fragments sorted by (area desc, lexicographic first point) for
  determinism.

**Step 4: Run** Debug/Release.

**Step 5: Commit**

`feat(m3a): cut planar B-rep faces into fragments`

## Task 6: Assemble fragments into a validated Brep

**Files:**
- Create: `include/geom/brep/assemble.hpp`
- Create: `tests/unit/brep_assemble_test.cpp`
- Modify: `tests/CMakeLists.txt`

**Interfaces:**
- Consumes: `SolidEditor`, `Plane3`, `LineCurve3`, `Pcurve2d`, `validate_*`.
- Produces:
  - `struct AssembledFace { std::vector<Point3d> outer; std::vector<std::vector<Point3d>> holes; Vec3d normal; };`
  - `Result<Brep> assemble_solid(std::span<const AssembledFace>, double tol);`

**Step 1: Write the failing test**

```cpp
namespace {
AssembledFace quad(std::array<Point3d, 4> p, Vec3d n) {
  return AssembledFace{{p[0], p[1], p[2], p[3]}, {}, n};
}
std::vector<AssembledFace> box_faces(Point3d lo, Point3d hi) {
  const double x0 = lo[0], y0 = lo[1], z0 = lo[2];
  const double x1 = hi[0], y1 = hi[1], z1 = hi[2];
  return {
      quad({{x0, y0, z1}, {x1, y0, z1}, {x1, y1, z1}, {x0, y1, z1}}, {0, 0, 1}),
      quad({{x0, y0, z0}, {x0, y1, z0}, {x1, y1, z0}, {x1, y0, z0}}, {0, 0, -1}),
      quad({{x0, y1, z0}, {x0, y1, z1}, {x1, y1, z1}, {x1, y1, z0}}, {0, 1, 0}),
      quad({{x0, y0, z0}, {x1, y0, z0}, {x1, y0, z1}, {x0, y0, z1}}, {0, -1, 0}),
      quad({{x1, y0, z0}, {x1, y1, z0}, {x1, y1, z1}, {x1, y0, z1}}, {1, 0, 0}),
      quad({{x0, y0, z0}, {x0, y0, z1}, {x0, y1, z1}, {x0, y1, z0}}, {-1, 0, 0}),
  };
}
}

TEST_CASE("assemble: six quads rebuild a closed box") {
  const auto faces = box_faces({0, 0, 0}, {2, 3, 4});
  const auto brep = assemble_solid(faces, 1e-9);
  REQUIRE(brep.has_value());
  CHECK(brep.value().vertex_count() == 8);
  CHECK(brep.value().edge_count() == 12);
  CHECK(brep.value().face_count() == 6);
  CHECK(brep.value().shell_count() == 1);
  CHECK(brep.value().solid_count() == 1);
  CHECK(mass_properties(brep.value()).value().volume ==
        doctest::Approx(24.0).epsilon(1e-12));
}

TEST_CASE("assemble: a missing face is a diagnosed not_closed failure") {
  auto faces = box_faces({0, 0, 0}, {1, 1, 1});
  faces.pop_back();
  const auto brep = assemble_solid(faces, 1e-9);
  REQUIRE_FALSE(brep.has_value());
  CHECK(brep.error().code == Errc::not_closed);
  CHECK(std::string_view(brep.error().message).find("edge") != std::string_view::npos);
}

TEST_CASE("assemble: welding joins independently computed intersection points") {
  // Two quads share an edge whose endpoints differ by 1e-12: build four
  // vertices and two triangles forming a square; the result must weld to four
  // vertices and fail closed (a single quad pair is an open sheet).
  const Vec3d n{0, 0, 1};
  const Point3d a{0, 0, 0}, b{1, 0, 0};
  const Point3d c{1 + 1e-12, 1, 0}, d{0, 1, 0};
  const std::vector<AssembledFace> faces{
      AssembledFace{{a, b, c}, {}, n}, AssembledFace{{a, c, d}, {}, n}};
  const auto brep = assemble_solid(faces, 1e-9);
  REQUIRE_FALSE(brep.has_value());
  CHECK(brep.error().code == Errc::not_closed);
}
```

**Step 2: Run it** — expected RED.

**Step 3: Implement**

- `VertexWelder`: snap each coordinate to `llround(x / tol)` and search the 27
  neighbouring cells for an existing vertex within `tol`; insert otherwise.
  Reject non-finite coordinates and `tol <= 0`.
- For every face: weld the outer and hole rings; create directed edges
  (previous → next). Count each undirected edge's uses and directions.
- Any edge with one use → `not_closed`, message
  `"assemble_solid: edge is used once"`, `context_a = edge index`.
  Any edge with more than two uses, or two uses in the same direction →
  `invalid_topology` with `"assemble_solid: non-manifold edge"`.
- Add vertices, line curves (`origin = p0`, `direction = (p1-p0)/|p1-p0|`,
  `t1 = |p1-p0|`), edges, planes (dedupe equal `(normal, offset)` within tol),
  faces and loops to a `SolidEditor`; `plane.reversed = false` because the
  assembled normal is already outward.
- Group faces into connected components through shared edges. Per component,
  sum the oriented signed volume; a positive component is an outer shell, a
  negative one is a cavity. Assign each cavity to the outer shell containing a
  sample point of the cavity; if none or more than one contains it → error.
- Emit one solid per outer shell (`add_solid({outer, cavities...}, tol)`) and
  `commit()`; return its `Result<Brep>`.
- Pcurve for each coedge: the edge endpoints in that face's `(u, v)` frame,
  ordered along the coedge direction.

**Step 4: Run** Debug/Release.

**Step 5: Commit**

`feat(m3a): assemble fragments into a validated Brep solid`

## Task 7: Boolean driver (transversal and coplanar)

**Files:**
- Create: `include/geom/brep/boolean_solid.hpp`
- Create: `tests/unit/brep_boolean_solid_test.cpp`
- Modify: `tests/CMakeLists.txt`

**Interfaces:**
- Consumes: Tasks 1–6.
- Produces:
  - `enum class SolidBooleanOp { union_, intersection, difference, symmetric_difference };`
  - `struct SolidBooleanOptions { double tolerance = 1e-9; };`
  - `struct SolidBooleanReport { std::uint32_t face_pairs = 0, cut_segments = 0, fragments = 0, kept = 0, dropped = 0, flipped = 0; double max_cut_residual = 0.0; };`
  - `struct SolidBooleanResult { Brep brep; SolidBooleanReport report; };`
  - `Result<SolidBooleanResult> boolean_solid(const Brep& a, const Brep& b, SolidBooleanOp op, const SolidBooleanOptions& options = {});`
  - `bool detail::solid_boolean_member(bool in_a, bool in_b, SolidBooleanOp op) noexcept;`

**Step 1: Write the failing test**

```cpp
TEST_CASE("boolean_solid: overlapping boxes, four operations") {
  const Brep a = make_box(Point3d{0, 0, 0}, Point3d{2, 2, 2}, 1e-12).value();
  const Brep b = make_box(Point3d{1, 1, 1}, Point3d{3, 3, 3}, 1e-12).value();
  struct Case { SolidBooleanOp op; double volume; };
  for (const Case c : {Case{SolidBooleanOp::union_, 15.0},
                       Case{SolidBooleanOp::intersection, 1.0},
                       Case{SolidBooleanOp::difference, 7.0},
                       Case{SolidBooleanOp::symmetric_difference, 14.0}}) {
    const auto result = boolean_solid(a, b, c.op);
    REQUIRE(result.has_value());
    const Brep& out = result.value().brep;
    CHECK(validate_topology(out).ok());
    CHECK(validate_tolerances(out).ok());
    CHECK(mass_properties(out).value().volume ==
          doctest::Approx(c.volume).epsilon(1e-9));
    CHECK(out.solid_count() == 1);
  }
}

TEST_CASE("boolean_solid: stacked boxes merge their shared face away") {
  const Brep a = make_box(Point3d{0, 0, 0}, Point3d{1, 1, 1}, 1e-12).value();
  const Brep b = make_box(Point3d{0, 0, 1}, Point3d{1, 1, 2}, 1e-12).value();
  const auto result = boolean_solid(a, b, SolidBooleanOp::union_);
  REQUIRE(result.has_value());
  CHECK(mass_properties(result.value().brep).value().volume ==
        doctest::Approx(2.0).epsilon(1e-9));
  CHECK(result.value().brep.shell_count() == 1);
}
```

**Step 2: Run it** — expected RED.

**Step 3: Implement**

- Validate options (`tol > 0`, finite) and both inputs
  (`validate_topology`, `validate_tolerances`); require exactly one solid and
  one shell each; require every face to be `Plane3`.
- Bounding-box prefilter over faces (AABBs expanded by `tol`); for each
  candidate pair call `intersect_planar_faces` and append the segments to both
  faces' cut lists.
- Cut every face with its list (empty list → original face fragment). Collect
  `(fragment, frame, face_normal)` records. Map each fragment ring back to 3D;
  the outward normal is `face_normal` flipped when `face.reversed` is true, and
  the rings are reversed when the flip changes the normal, so every assembled
  outer ring is CCW about its outward normal and every hole ring is CW.
- Sample: `p3 = plane_frame_to_3d(frame, fragment.sample)`;
  `delta = max(64 * tol, 1e-7 * bbox_diagonal)`. Compute
  `in_a_plus`, `in_b_plus` at `p3 + n*delta`, and `in_a_minus`, `in_b_minus` at
  `p3 - n*delta` via `point_in_solid`.
- `f_plus = member(in_a_plus, in_b_plus, op)`, likewise minus. Drop when equal.
  Keep when `f_plus == false` (normal stays); when `f_plus == true` flip the
  ring orientation and the normal, and increment `flipped`.
- Deduplicate coincident fragments (same plane within tol and same ring
  multiset); keep the first in deterministic order.
- Feed the kept `AssembledFace`s to `assemble_solid`; validate; return.
- Empty kept set → `Error{Errc::degenerate, "boolean_solid: result is empty", ...}`.

**Step 4: Run** Debug/Release.

**Step 5: Commit**

`feat(m3a): add planar B-rep boolean for union, intersection and difference`

## Task 8: Multi-shell results, cavities and through-holes

**Files:**
- Modify: `include/geom/brep/validate.hpp` (allow multi-shell solids,
  orientation and containment checks)
- Modify: `include/geom/brep/assemble.hpp` (cavity assignment; already required
  by Task 6, wire it to the validator)
- Create: `tests/unit/brep_boolean_shell_test.cpp`
- Modify: `tests/CMakeLists.txt`

**Interfaces:**
- Consumes: Task 7, `point_in_shell`, `euler_characteristic`.
- Produces:
  - `validate_topology` accepts solids with one positive shell and any number
    of negative (cavity) shells; rejects two positive shells in one solid.
  - `euler_characteristic` unchanged in shape; documented as
    `V - E + F - (L - F) = 2(S - H)` with `S` the shell count.

**Step 1: Write the failing test**

```cpp
TEST_CASE("boolean_shell: a fully interior difference opens a cavity") {
  const Brep a = make_box(Point3d{0, 0, 0}, Point3d{10, 10, 10}, 1e-12).value();
  const Brep b = make_box(Point3d{3, 3, 3}, Point3d{7, 7, 7}, 1e-12).value();
  const auto result = boolean_solid(a, b, SolidBooleanOp::difference);
  REQUIRE(result.has_value());
  const Brep& out = result.value().brep;
  CHECK(out.shell_count() == 2);
  CHECK(out.solid_count() == 1);
  CHECK(validate_topology(out).ok());
  CHECK(mass_properties(out).value().volume == doctest::Approx(936.0).epsilon(1e-9));
  const SolidId solid = out.solid_ids().front();
  CHECK(euler_characteristic(out, solid) == 4);  // 2 * (2 shells - 0 handles)
}

TEST_CASE("boolean_shell: a difference that punches a through-hole has genus one") {
  const Brep plate = make_box(Point3d{0, 0, 0}, Point3d{10, 10, 2}, 1e-12).value();
  const Brep punch = make_box(Point3d{4, 4, -1}, Point3d{6, 6, 3}, 1e-12).value();
  const auto result = boolean_solid(plate, punch, SolidBooleanOp::difference);
  REQUIRE(result.has_value());
  const Brep& out = result.value().brep;
  CHECK(out.shell_count() == 1);
  CHECK(mass_properties(out).value().volume == doctest::Approx(192.0).epsilon(1e-9));
  CHECK(euler_characteristic(out, out.solid_ids().front()) == 0);  // 2 - 2*1
}
```

**Step 2: Run it** — expected RED: the difference returns `unsupported` from
`commit()` (M2 multi-shell guard), or the cavity is mis-grouped.

**Step 3: Implement**

- Remove the `shell_count != 1` rejection. For each solid: exactly one shell
  with positive signed volume (outer), every other shell negative, every cavity
  contained in the outer shell (`point_in_shell` of a sample point on the
  cavity side), cavities pairwise disjoint. Two positive shells →
  `invalid_topology` with `"solid has more than one outer shell"`.
- Keep `mass_properties` summing all oriented face contributions; add a test
  that the cavity shell contributes a negative volume.
- In `assemble_solid`, group components by signed volume, assign cavities by
  containment, and emit one solid per outer shell.

**Step 4: Run** Debug/Release.

**Step 5: Commit**

`feat(m3a): support solids with cavities and through-holes`

## Task 9: Diagnostics and degenerate inputs

**Files:**
- Modify: `include/geom/brep/boolean_solid.hpp`
- Create: `tests/unit/brep_boolean_solid_errors_test.cpp`
- Modify: `tests/CMakeLists.txt`

**Interfaces:**
- Consumes: Task 7.
- Produces: stable `Errc` + message contract for each failure class:
  - curved operand → `unsupported` / `"boolean_solid: curved faces are deferred to M3b"`
  - multi-shell operand → `unsupported` / `"boolean_solid: multi-shell operands are deferred to M3c"`
  - empty result → `degenerate` / `"boolean_solid: result is empty"`
  - open/non-manifold assembly → `not_closed` or `invalid_topology`
  - invalid tolerance → `tolerance_violation` / `"boolean_solid: tolerance must be positive and finite"`

**Step 1: Write the failing test**

```cpp
TEST_CASE("boolean_solid errors: curved operands are unsupported") {
  const Brep box = make_box(Point3d{0, 0, 0}, Point3d{1, 1, 1}, 1e-12).value();
  const Brep cylinder = make_cylinder(Point3d{0, 0, 0}, Vec3d{0, 0, 1}, 1.0, 2.0, 1e-12).value();
  const auto result = boolean_solid(box, cylinder, SolidBooleanOp::union_);
  REQUIRE_FALSE(result.has_value());
  CHECK(result.error().code == Errc::unsupported);
  CHECK(std::string(result.error().message) == "boolean_solid: curved faces are deferred to M3b");
}

TEST_CASE("boolean_solid errors: subtracting an identical box is empty") {
  const Brep a = make_box({0, 0, 0}, {1, 1, 1}, 1e-12).value();
  const Brep b = make_box({0, 0, 0}, {1, 1, 1}, 1e-12).value();
  const auto result = boolean_solid(a, b, SolidBooleanOp::difference);
  REQUIRE_FALSE(result.has_value());
  CHECK(result.error().code == Errc::degenerate);
}

TEST_CASE("boolean_solid errors: boxes that meet only at an edge are not silent") {
  const Brep a = make_box({0, 0, 0}, {1, 1, 1}, 1e-12).value();
  const Brep b = make_box({1, 1, 0}, {2, 2, 1}, 1e-12).value();
  const auto result = boolean_solid(a, b, SolidBooleanOp::union_);
  REQUIRE_FALSE(result.has_value());
  CHECK((result.error().code == Errc::not_closed ||
         result.error().code == Errc::invalid_topology));
}
```

**Step 2: Run it** — some cases will fail with a wrong code or a crash.

**Step 3: Implement** the guards in the listed order (options → operands →
faces → empty → assembly). Never let an `ASSERT`/`std::terminate` path be the
diagnosis.

**Step 4: Run** Debug/Release.

**Step 5: Commit**

`feat(m3a): diagnose unsupported and degenerate boolean inputs`

## Task 10: QA-8 benchmark

**Files:**
- Create: `bench/bench_m3.cpp`
- Modify: `CMakeLists.txt` (`GEOM_BUILD_BENCH` target)
- Create: `docs/benchmarks/m3a.md`

**Interfaces:**
- Consumes: `boolean_solid`, `mass_properties`, `SolidEditor`.
- Produces: a program printing
  `QA8 faces_a=<n> faces_b=<n> ms=<t> volume_rel_err=<e> target=5000ms/1e-9 PASS|FAIL`.

**Step 1: Write the benchmark**

- Build a `prism(center, radius, height, sides, tol)` solid with
  `sides = 5000`: two cap faces with `sides` coedges each and `sides` side
  quads (≈ 5002 faces per solid).
- Operands: two prisms with the same axis, offset by `0.25 * radius`.
- Run `union_`, `intersection`, `difference`; compute
  `|V(union) + V(intersection) - Va - Vb| / (Va + Vb)` and
  `|V(difference) + V(intersection) - Va| / Va`.
- Time only the boolean calls; print the numbers; exit 1 on FAIL, but do not
  throw.

**Step 2: Run it**

```powershell
cmake --build build-rel --config Release --target geom_bench_m3
build-rel\bench\Release\geom_bench_m3.exe
```

Record the output in `docs/benchmarks/m3a.md`, including the machine and the
fact that FAIL is FAIL.

**Step 3: Record** — nothing to implement beyond the bench; if the time is over
5 s or the error over 1e-9, state it explicitly.

**Step 4: Run** Debug build of the library must stay green; the bench is
Release-only.

**Step 5: Commit**

`feat(m3a): add QA-8 boolean benchmark and record measurements`

## Task 11: Documentation and ADR

**Files:**
- Create: `docs/adr/0023-planar-boolean-first.md`
- Modify: `docs/architecture/geometry-library-architecture.md` §16.3
- Modify: `docs/superpowers/plans/2026-11-09-m3a-planar-boolean.md` (checkboxes)
- Modify: `docs/benchmarks/m3a.md`

**Content requirements:**

- ADR-0023: decision to ship planar boolean first; alternatives (SSI-first,
  mesh boolean); consequences (curved faces rejected, QA-9 deferred, pcurves
  straight); revisit trigger (M3b).
- §16.3: M3a status line with test count, QA-8 PASS/FAIL, deformation scope.
- Plan: tick completed steps; leave QA-8 unchecked if it failed.

**Step 1–5:** normal TDD does not apply; verify the docs match the measured
numbers, then commit:

`docs(m3a): record planar boolean scope and QA-8 evidence`

## M3a Exit Acceptance

- Debug and Release full suite green, zero C++ warnings, `/W4 /WX /utf-8`.
- All new headers compile standalone.
- `validate_topology` and `validate_tolerances` pass on every boolean result
  the tests expect to succeed.
- Volume identity holds within 1e-9 for the test corpus (transversal, stacked,
  cavity, through-hole, disjoint bodies).
- Failure cases return the documented `Errc` and message, never a crash and
  never a silent empty `Brep`.
- QA-8 measured and recorded; if it FAILs, `docs/benchmarks/m3a.md` and §16.3
  say FAIL and the gap is listed as follow-up work.
- `docs/adr/0023-*.md` exists and §16.3 links it.

## Self-Check

- [x] Every task has Files, Interfaces, five steps and one commit.
- [x] Every test snippet compiles against the interfaces listed in the same
      task and the existing headers (`Handle`, `Result`, `Plane3`, `Brep`).
- [x] No task depends on a later task's header.
- [x] `face_intersect` consumes `face_query`; `inside` consumes `face_query`;
      `planar_cut` consumes `face_intersect` + `arrangement2`; `assemble`
      consumes `solid_editor`; `boolean_solid` consumes all of them.
- [x] Every new public function has a failure test (Task 9 collects them).
- [x] No global epsilon; every tolerance is a parameter.
- [ ] Deviations 1–6 are acceptable to the reviewer.
