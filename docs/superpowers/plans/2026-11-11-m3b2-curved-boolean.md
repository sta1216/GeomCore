# M3b-2 Curved Boolean Implementation Plan

> **For the implementing agent:** follow this plan task by task, in order. Every task starts with a failing test, then the minimal implementation that makes it pass, then a single commit. Do not skip the RED step and do not batch two tasks into one commit.

## Goal

Turn the M3b-1 SSI engine into trimmed curved faces and a boolean of two analytic B-rep solids:

1. `Pcurve2d` becomes a polyline so trimmed edges can carry multi-segment parameter curves.
2. A face's loops can be read in surface parameter space (`face_loops_uv`), including periodic seams and edges that have no stored pcurve.
3. A face can be cut by SSI curves into fragments (`split_face_curved`).
4. Fragments from both operands are welded into a validated `Brep` (`assemble_curved_solid`).
5. A point-in-mesh classifier and the boolean driver (`boolean_solid_curved`) keep or drop fragments by the regularized membership rule.

M3b-2 stops at analytic surfaces. NURBS faces return `unsupported`; cone operands are deferred because their lateral loops pass through the apex (see deviation 4); coplanar face contact and results with more than one closed component are out of scope and must be diagnosed, never silently accepted.

## Architecture

- **Layer:** L5 B-rep algorithms, same layer as `boolean_solid` (M3a). `geom/brep/detail/*` holds non-public helpers; public headers are `pcurve2d.hpp`, `face_uv.hpp`, `trim_face.hpp`, `assemble_curved.hpp`, `split_solid.hpp`, `mesh_inside.hpp`, `boolean_solid_curved.hpp`.
- **Trimming is parameter-domain 2D boolean.** Each 3D cut curve produced by `ssi()` is turned into a polyline in the face's `(u, v)` domain. The face's loop polylines plus the cut polylines are fed to `Arrangement2` as sources 0..N. The arrangement gives exact-noded rings; a containment tree over those rings reconstructs outer loops and holes. Every vertex is mapped back to 3D with `surface_eval`, so the 3D fragment loops are always on the surface to within the tolerance.
- **Assembly is weld-by-geometry.** M3b-2 does not mutate either operand's topology. Fragments are collected into `CurvedAssembledFace` values; `assemble_curved_solid` welds their 3D vertices within `tol`, creates one polyline NURBS edge per shared polyline, stores pcurves in edge `v0 -> v1` direction, groups faces that share edges into one shell, and calls `validate_topology`.
- **Classification is explicit and off-surface.** Fragment representatives are centroids of the first ear of the outer loop. They are classified against the other operand's tessellation. A representative within `tol` of the other solid's surface is a `numeric_failure`, never a guess.
- **Diagnostics over silent output.** Every unsupported configuration (NURBS, coplanar contact, multi-component result, non-manifold weld, ambiguous classification) returns an `Error` with a stable `Errc` and message. An empty but valid result is allowed only for a true disjoint operation and is reported through the result counters.

## Tech Stack

- C++20, header-only public API, CMake 3.25+ with the existing `geom_core` interface target.
- doctest v2.4.11; tests registered in `tests/CMakeLists.txt` under `GEOM_TEST_SOURCES`.
- MSVC 19.44 `/W4 /WX /utf-8`; GCC/Clang `-Wall -Wextra -Wpedantic -Werror` already configured.
- Existing primitives: `Point3d`, `Vec3d`, `Plane3`, `Cylinder3`, `Cone3`, `Sphere3`, `Torus3`, `NurbsCurve`, `NurbsSurface`, `Brep`, `SolidEditor`, `validate_topology`, `ssi()`, `Arrangement2`, `tessellate`, `polygon_area2`, `sign_orient2d`, `triangulate_polygon2`.

## Spec

Source of truth: `docs/architecture/geometry-library-architecture.md` v2.1

- §8.7 "修剪/分割": intersection curve -> parameter-domain 2D boolean -> new loops -> new coedges/pcurves.
- §8.8 "曲面布尔": face pairing -> SSI -> trim -> classify -> keep/drop -> stitch -> validate.
- §16.3 M3b-2 scope: curved boolean, multi-segment pcurves, mixed-scale duplicate-curve residuals.
- `docs/quality/testing-strategy.md`: every public API needs positive, boundary and failure tests; every `Errc` path needs at least one test; no silent success on invalid input; degeneracies must be asserted, not skipped.

## Global Constraints

- All new code is header-only and inside `namespace geom { inline namespace v1 {`.
- No new third-party dependency.
- Tolerance is always a parameter; no global epsilon, no hard-coded absolute thresholds except documented machine-epsilon floors.
- Topological decisions (inside/outside, crossing) use exact predicates (`sign_orient2d`, `sign_incircle`) or the existing exact arrangement; floating point is used for geometry construction and for rejecting clearly separated cases.
- Public functions return `Result<T>` and never throw.
- Existing M1b/M2/M3a behaviour and messages must not change; the full suite must stay green.
- Every new test file is added to `GEOM_TEST_SOURCES`; every task ends with Debug and Release `ctest` green and one commit.
- Do not leave commented-out code, TODO markers, or tests that assert `>= 0` on a counter that is always non-negative.

## Scope, splits and approved deviations

M3b is split into three plans:

- **M3b-2 (this plan):** polyline pcurves, parameter-domain trimming, analytic curved boolean.
- **M3b-3:** healing, shell repair, QA-12.
- **M3c:** NURBS operand booleans and exact rational SSI (blocked on the M3b-1 NURBS coincident predicate).

Deviations from a literal reading of the architecture that require review:

1. **Closed edges are sampled, not preserved analytically.** A circular cap edge is sampled into a 64-segment polyline in the parameter domain; the resulting trimmed face has many short line edges instead of one analytic circle. This is exact within `tol` for trimming purposes but changes edge count. It is documented in the header and measured by the chord-error test.
2. **No coplanar face merging.** A pair of coplanar face patches that overlap in area returns `unsupported "boolean_solid_curved: coplanar face contact is not supported yet"`; M3a's planar boolean still handles all-planar operands. This keeps M3b-2 from re-implementing M3a's coplanar splitter in parameter space.
3. **Single connected result only.** A union whose result is disjoint, or a difference that creates a cavity, returns `unsupported "boolean_solid_curved: result is not a single connected component"`. M3a keeps its multi-solid support; M3b-2 defers it until assembly can prove component orientation for curved shells without tessellating masses.
4. **Cone operands are deferred.** Every `make_cone` lateral face has a loop that runs through the apex, where the surface parameterization has no unique `u`; `face_loops_uv` returns `unsupported` for that loop and the boolean driver rejects `Cone3` before splitting. Apex-free cone patches are supported by `face_loops_uv` and tested; full cone booleans need a degenerate-vertex representation in the parameter domain and are left to M3c.

These four deviations are also listed in ADR-0025.

## Review Focus

The reviewer must attack these five areas first:

1. **Periodic seams.** A face whose loop crosses `u = 0 / 2π` (full cylinder lateral, sphere) must unroll into one continuous loop, not two; pcurves must be stored unwrapped. Probe a full cylinder and a sphere cap with a cut at `u = 0`.
2. **Cut chains ending on a face boundary.** A cut that exits through a loop edge or a vertex must be snapped to that boundary without creating zero-area fragments or duplicate rings. Probe plane cuts through a box corner and through a cylinder cap seam.
3. **Pcurve direction.** Every pcurve must follow the edge from `v0` to `v1`; a reversed coedge sets `reversed` but never reverses the stored pcurve. Probe `split_edge` on an edge created by `assemble_curved_solid`, then re-weld and validate.
4. **Classification near the other surface.** Representatives within `tol` of the other solid must produce a diagnostic; a tangency must not be classified by luck. Probe sphere tangent to a plane and box touching a cylinder.
5. **Sliver and zero-area fragments.** Fragments whose normalized area is below `tol²` must be dropped and counted, not emitted as faces that later fail `validate_topology`. Probe a cut grazing a face within `tol`.

## File Map

| Path | Purpose |
| --- | --- |
| `include/geom/brep/pcurve2d.hpp` | polyline pcurve operations |
| `include/geom/brep/detail/curve_polyline.hpp` | degree-1 NURBS polyline curve factory |
| `include/geom/brep/face_uv.hpp` | face loops in surface parameter space |
| `include/geom/brep/trim_face.hpp` | parameter-domain face splitting |
| `include/geom/brep/assemble_curved.hpp` | weld curved fragments into a `Brep` |
| `include/geom/brep/split_solid.hpp` | split a solid by a plane (integration path) |
| `include/geom/brep/mesh_inside.hpp` | point-in-tessellation classifier |
| `include/geom/brep/boolean_solid_curved.hpp` | boolean driver |
| `tests/unit/brep_pcurve2d_test.cpp` | Task 1 |
| `tests/unit/brep_polyline_curve_test.cpp` | Task 2 |
| `tests/unit/brep_face_uv_test.cpp` | Task 3 |
| `tests/unit/brep_trim_face_test.cpp` | Task 4 |
| `tests/unit/brep_assemble_curved_test.cpp` | Task 5 |
| `tests/unit/brep_split_solid_test.cpp` | Task 6 |
| `tests/unit/brep_mesh_inside_test.cpp` | Task 7 |
| `tests/unit/brep_boolean_curved_test.cpp` | Task 8 |
| `tests/unit/brep_curved_fixture.hpp` | shared quarter-cylinder fixture (test-only) |
| `bench/bench_boolean_curved.cpp` | Task 9 |
| `docs/benchmarks/m3b2-curved-boolean.md` | Task 9 |
| `docs/adr/0025-curved-boolean-analytic.md` | Task 9 |

---

## Task dependencies

| Task | Consumes (from earlier tasks) |
| --- | --- |
| 1 | M2 `Pcurve2d`, `SolidEditor` transaction |
| 2 | `NurbsCurve` |
| 3 | `Brep`, `SolidEditor`, `CircleCurve3`, `make_plane_frame` |
| 4 | Task 3 `face_loops_uv`; M1b `Arrangement2`, `split_segments` |
| 5 | Task 1 `Pcurve2d`, Task 2 `make_polyline_curve` |
| 6 | Tasks 3-5; M3b-1 `ssi` |
| 7 | M2 `tessellate`, `TriMesh` |
| 8 | Tasks 4-7; M3a `SolidBooleanOp`, `solid_boolean_member` |
| 9 | Tasks 1-8 |

---

## Task 1: Polyline pcurves

**Files**

- Modify: `include/geom/brep/topology.hpp`
- Create: `include/geom/brep/pcurve2d.hpp`
- Modify: `include/geom/brep/solid_editor.hpp` (`split_edge`, `join_vertex`)
- Create: `tests/unit/brep_pcurve2d_test.cpp`
- Modify: `tests/unit/brep_euler_ops_test.cpp`, `tests/unit/brep_assemble_test.cpp` (mechanical p0/p1 -> points.front()/points.back())
- Modify: `tests/CMakeLists.txt`

**Interfaces**

```cpp
// topology.hpp
struct Pcurve2d {
  std::vector<Point2d> points;
  Pcurve2d() = default;
  Pcurve2d(std::initializer_list<Point2d> pts) : points(pts) {}
  explicit Pcurve2d(std::vector<Point2d> pts) : points(std::move(pts)) {}
};
```

```cpp
// pcurve2d.hpp
double pcurve_length(const Pcurve2d& pcurve) noexcept;   // 0 for a malformed curve
Result<Pcurve2d> pcurve_split(const Pcurve2d& pcurve, double t);  // t clamped to [0, 1]
Pcurve2d pcurve_merge(const Pcurve2d& lower, const Pcurve2d& upper);
bool pcurve_well_formed(const Pcurve2d& pcurve) noexcept;  // >= 2 finite points
```

Semantics:

- `pcurve_split` walks chord length; `t` is clamped to `[0, 1]`. The split point is inserted once on both sides, so `pcurve_merge(lower, upper)` drops the duplicate junction when the two ends are exactly equal.
- A malformed pcurve (fewer than two points, non-finite coordinate) makes `pcurve_split` return `non_finite` / `degenerate`; `pcurve_length` returns 0 rather than NaN.
- `SolidEditor::split_edge` splits every attached pcurve with `pcurve_split(source, t)`; `join_vertex` rebuilds with `pcurve_merge(lower, upper)`. A failed split aborts `split_edge` with that error and leaves the editor unchanged (the existing transaction guarantee).

**Steps**

1. Add this test (RED first):

```cpp
TEST_CASE("pcurve_split follows chord length and clamps t") {
  geom::Pcurve2d pc{{0.0, 0.0}, {1.0, 0.0}, {1.0, 2.0}};
  CHECK(geom::pcurve_length(pc) == doctest::Approx(3.0));
  const auto half = geom::pcurve_split(pc, 0.5);
  REQUIRE(half.has_value());
  CHECK(half->points.size() == 3);
  CHECK(half->points.back()[0] == doctest::Approx(1.0));
  CHECK(half->points.back()[1] == doctest::Approx(0.5));
  const auto over = geom::pcurve_split(pc, 4.0);
  REQUIRE(over.has_value());
  CHECK(over->points.back() == pc.points.back());
}
```

Also test `pcurve_merge(half, upper)` reproduces the three original points plus the inserted junction (4 points), equal geometry, and non-finite input rejection.

2. Build: `cmake --build build --config Debug` must fail because `Pcurve2d` has no `points` yet. Capture that failure in the report.
3. Change `Pcurve2d`, add `pcurve2d.hpp`, update `solid_editor.hpp`.
4. Run the new test, then the whole suite in Debug and Release. Existing `Pcurve2d{{0,0},{1,1}}` call sites must still compile through `initializer_list`.
5. Commit: `feat(m3b2): make pcurves polylines`.

## Task 2: Polyline curve factory

**Files**

- Create: `include/geom/brep/detail/curve_polyline.hpp`
- Create: `tests/unit/brep_polyline_curve_test.cpp`
- Modify: `tests/CMakeLists.txt`

**Interfaces**

```cpp
Result<NurbsCurve> make_polyline_curve(std::span<const Point3d> points);
```

Semantics: degree 1, clamped uniform knots in `[0, 1]`, control points are the input; `eval(0)` and `eval(1)` are the first and last points; the curve passes through every control point at `i / (n - 1)`. Fewer than 2 points -> `degenerate`; non-finite point -> `non_finite`.

**Steps**

1. Test:

```cpp
TEST_CASE("a polyline curve interpolates every vertex") {
  const std::array<geom::Point3d, 3> pts{
      geom::Point3d{0, 0, 0}, geom::Point3d{1, 0, 0}, geom::Point3d{1, 1, 0}};
  const auto curve = geom::make_polyline_curve(pts);
  REQUIRE(curve.has_value());
  CHECK(curve->degree == 1);
  const auto p = curve->eval(0.5);
  CHECK(p[0] == doctest::Approx(1.0));
  CHECK(p[1] == doctest::Approx(0.0));
  const auto q = curve->eval(1.0);
  CHECK(q[1] == doctest::Approx(1.0));
}
```

2. RED: missing header.
3. Implement; verify `eval` at `i/(n-1)` equals `pts[i]` for all i (add the loop assertion).
4. Full suite Debug/Release.
5. Commit: `feat(m3b2): add polyline curve factory`.

## Task 3: Face loops in parameter space

**Files**

- Create: `include/geom/brep/face_uv.hpp`
- Create: `tests/unit/brep_face_uv_test.cpp`
- Modify: `tests/CMakeLists.txt`

**Interfaces**

```cpp
struct FaceUvLoop {
  std::vector<Point2d> points;
  bool outer = false;
};

// Returns one loop per face boundary. `outer` mirrors Loop::outer. Loops are
// unwrapped across periodic seams: consecutive u values differ by less than pi.
Result<std::vector<FaceUvLoop>> face_loops_uv(const Brep& brep, FaceId face,
                                              double tolerance);
```

Implementation:

- Walk each loop via `Loop::start` and `Coedge::next`.
- For a coedge with a stored pcurve, copy `Pcurve2d::points`; otherwise sample the edge curve: `LineCurve3` -> 2 endpoints, `CircleCurve3` -> 32 segments, `NurbsCurve` -> 32 uniform samples. Reverse the sequence when `Coedge::reversed`.
- Drop the duplicated junction point between consecutive coedges (within `tolerance` in 3D, exact in uv when stored pcurves exist).
- If the edge has no pcurve, project sampled 3D points with a per-surface inverse: plane frame projection; cylinder `u = atan2(dot(d,y), dot(d,x))`, `v = dot(d,axis)`; cone/sphere/torus use their own `du/dv` framing; a `NurbsSurface` face returns `unsupported`.
- Unwrap `u` by ±2π so consecutive values are within π; unwrap `v` the same way for a torus. Cylinder/cone/sphere `v` is not periodic.
- Close the loop by checking the last point against the first within `max(tolerance, 1e-12 * loop_extent)`; otherwise `numeric_failure "face_loops_uv: loop is not closed"`.
- Reject non-finite coordinates with `non_finite`; reject a face with no loop with `invalid_topology`.
- Apex-touching cone loops: a `make_cone` lateral face is bounded by the base circle plus a seam that runs base -> apex -> base, so its final point has `rho = 0` and no unique `u`. Such a loop returns `unsupported "face_loops_uv: cone projection is undefined at or beyond the apex"`; `Cone3` is consequently outside the M3b-2 boolean operand set (see Task 8). An apex-free cone patch (a closed circle at `rho > 0`) is supported and tested.

**Steps**

1. Test with `make_box(Point3d{-1,-1,-1}, Point3d{1,1,1}, 1e-9)`:

```cpp
TEST_CASE("box faces expose four uv corners") {
  const auto box = geom::make_box(geom::Point3d{-1, -1, -1}, geom::Point3d{1, 1, 1}, 1e-9);
  REQUIRE(box.has_value());
  const auto loops = geom::face_loops_uv(*box, geom::FaceId{0}, 1e-9);
  REQUIRE(loops.has_value());
  REQUIRE(loops->size() == 1);
  CHECK((*loops)[0].outer);
  CHECK((*loops)[0].points.size() == 4);
}
```

2. Add a cylinder lateral test: `make_cylinder` face 2 is the lateral face; its uv loop has four runs and a total u span of `2π`; assert `max_u - min_u == doctest::Approx(2π)` and that no consecutive jump exceeds π.
3. RED; implement; Debug/Release.
4. Commit: `feat(m3b2): expose face loops in surface parameters`.

## Task 4: Split a face by SSI cuts

**Files**

- Create: `include/geom/brep/trim_face.hpp`
- Create: `tests/unit/brep_trim_face_test.cpp`
- Modify: `tests/CMakeLists.txt`

**Interfaces**

```cpp
struct FaceCut {
  std::vector<Point2d> uv;   // polyline in this face's parameters
};

struct CurvedFragment {
  std::vector<std::vector<Point3d>> loops;     // outer first, then holes
  std::vector<std::vector<Point2d>> loops_uv;
  bool reversed = false;                        // orientation flipped by the caller
};

Result<std::vector<CurvedFragment>> split_face_curved(
    const Brep& brep, FaceId face, std::span<const FaceCut> cuts, double tolerance);
```

Algorithm:

1. `face_loops_uv(brep, face, tolerance)`; require exactly one outer loop.
2. Build `ArrangementSegment`s: loop points as source 0 (orientation preserved); each cut polyline as sources 1..N. Snap a cut endpoint that lies within `tolerance` of a loop segment by replacing its coordinate with the exact projection onto that segment (clamped to its ends), so `split_segments` sees an exact touch and the loop splits there; a cut endpoint farther away stays free (it must meet another cut or close a loop).
3. Drop cut segments shorter than `tolerance`; drop fragments whose absolute area is below `tolerance * tolerance`.
4. `Arrangement2::build`. For each bounded face ring with winding(0) != 0, collect its ring polygon and area. Also collect ring vertices.
5. Representative point per ring: triangulate the ring with `triangulate_polygon2` (outer only); the centroid of the first triangle is inside. Classify `winding(0)` at that point (0 = void ring, >0 = material ring).
6. Build the containment tree: parent = smallest-area ring that contains the representative point of the child. Fragments = material rings; holes = direct void children. A material ring contained in a void ring is its own fragment (island).
7. For each fragment, order holes deterministically (area descending, then lexicographic first point). Emit 3D loops with `surface_eval(surface, uv)`; reject non-finite evaluation with `numeric_failure`.
8. `reversed` is false: the arrangement outer ring is already counter-clockwise in uv, which maps to the face's natural normal. The boolean driver sets it when the fragment must be flipped.

**Steps**

1. Fixture: create `tests/unit/brep_curved_fixture.hpp` with `test::make_quarter_cylinder(double radius, double height, double tol)` and `test::quarter_cylinder_lateral(const Brep&)`. The solid is the quarter cylinder `u in [0, pi/2]`, `v in [0, height]` around the `+z` axis at the origin:
   - lateral face: `Cylinder3{origin, +z, +x, radius}` with an outer loop of four polylines (u = 0 line, top arc at v = height, u = pi/2 line, bottom arc at v = 0); arcs use 16 chords each;
   - top and bottom caps: `Plane3` faces, quarter-disk outer loops (two radii + the same 16-chord arc);
   - side faces at u = 0 and u = pi/2: `Plane3` rectangles;
   - every coedge gets a `Pcurve2d` in its own surface parameters (bottom arc `(t, 0)`, top arc `(t, height)`, side lines `(u, t)`), so the fixture does not depend on Task 3's projection path;
   - `validate_topology` must be clean. The 16-chord arcs make the exact volume `height * (n / 2) * radius^2 * sin(pi / (2n))` with `n = 16`, which is about 0.16% below the analytic `pi * radius^2 * height / 4`; the analytic value is *not* the volume of this B-rep. Trimmed curved faces cannot be tessellated by the M2 `tessellate` (it handles full analytic primitives only), so Task 4b adds the volume query.
   Put the helper in the header only (no CMake entry) so Tasks 4–8 include it.
2. Test:

```cpp
TEST_CASE("a lateral cut splits a cylindrical face into two fragments") {
  const auto solid = geom::test::make_quarter_cylinder(1.0, 2.0, 1e-9);
  REQUIRE(solid.has_value());
  const geom::FaceId lateral = geom::test::quarter_cylinder_lateral(*solid);
  const std::array<geom::FaceCut, 1> cuts{geom::FaceCut{{{0.0, 0.5}, {1.5707963267948966, 0.5}}}};
  const auto parts = geom::split_face_curved(*solid, lateral, cuts, 1e-9);
  REQUIRE(parts.has_value());
  REQUIRE(parts->size() == 2);
  CHECK((*parts)[0].loops.size() == 1);
  CHECK((*parts)[1].loops.size() == 1);
}
```

3. Add a hole test: cut a closed square loop strictly inside a box face; expect one fragment with one hole; assert `loops.size() == 2` and both ring counts.
4. Add a boundary-touching cut test: a cut from one loop edge to the opposite edge splits the face into two fragments, no holes.
5. RED; implement; Debug/Release.
6. Commit: `feat(m3b2): split curved faces in parameter space`.

## Task 4b: Volume by parameter-space quadrature

The M2 `tessellate` only meshes full analytic primitives (a trimmed cylinder patch with polyline boundaries fails with `tessellate: cylinder face is missing a cap circle`) and `mass_properties` only integrates untrimmed rectangles. A Boolean result has trimmed curved faces, so both are unavailable for verification. The divergence theorem gives the volume from the face loops alone: `V = (1/3) * sum_faces integral(S . (Su x Sv) du dv)` over the exact uv polygon. The polygon boundary is piecewise linear in uv (pcurves are polylines), so the domain is exact and a degree-5 rule per triangle is accurate for analytic surfaces.

**Files**

- Create: `include/geom/brep/detail/curved_volume.hpp`
- Create: `tests/unit/brep_curved_volume_test.cpp`
- Modify: `tests/CMakeLists.txt`

**Interfaces**

```cpp
namespace geom::v1::detail {
// Signed volume of one solid, integrated over each face's exact uv polygon.
// Analytic surfaces only (Plane3, Cylinder3, Cone3, Sphere3, Torus3).
Result<double> curved_solid_volume(const Brep& brep, SolidId solid, double tolerance);
}
```

Algorithm:

1. Reject a null/foreign `SolidId` (`index_out_of_range`), a missing/zero shell, an open shell (`invalid_topology`), `tolerance` not finite or `<= 0` (`tolerance_violation`), and any face whose surface is `NurbsSurface` (`unsupported "curved_solid_volume: NURBS surfaces are not supported"`).
2. For each face: `face_loops_uv(brep, face, tolerance)`; require one outer loop and build `Polygon2d{outer, holes}` with the returned orientation (outer counter-clockwise, holes clockwise).
3. `triangulate_polygon2(polygon)`; failure is `degenerate "curved_solid_volume: face triangulation failed"` carrying the face index.
4. For each triangle with uv vertices `a, b, c`: map a 7-point Dunavant rule (degree 5: barycentric `(1/3,1/3,1/3)` weight `0.225`; the three `(0.05971587,0.47014206,0.47014206)` permutations weight `0.13239415`; the three `(0.79742699,0.10128651,0.10128651)` permutations weight `0.12593918`) to uv, evaluate `surface_eval`, `surface_du`, `surface_dv`, and accumulate `2 * area_uv * w * dot(S, cross(Su, Sv)) / 3`.
5. Negate a face's contribution when `face.reversed` is set (the loop orientation then follows the reversed normal).
6. Sum in deterministic face order; require a finite positive total, otherwise `numeric_failure "curved_solid_volume: volume is not finite"` / `"curved_solid_volume: volume is not positive"`.

**Steps**

1. Test `make_box(2, 2, 2)` -> 8 within `1e-12` relative (planar faces take the same quadrature path).
2. Test `make_cylinder(2.0, 3.0, 1e-9)` -> the 32-gon prism volume `h * (n / 2) * r^2 * sin(2 * pi / n)` with `n` equal to the pcurve sample count used by `face_loops_uv`, within `1e-9` relative; also assert it is within `1e-3` relative of `pi * r^2 * h` (the sampling error).
3. Test the quarter-cylinder fixture -> `height * (n / 2) * radius^2 * sin(pi / (2n))` with `n = 16` within `1e-9`, and within `2e-3` of `pi * radius^2 * height / 4`.
4. Error paths: null `SolidId`, empty `Brep`, `tolerance` 0/NaN, and (if constructible without a NURBS face) NURBS rejection.
5. RED; implement; Debug/Release.
6. Commit: `feat(m3b2): integrate curved solid volume in parameter space`.

Note for Task 9 docs: trimmed curved faces remain non-tessellatable and `mass_properties` remains limited to untrimmed rectangles; record this as an M3c item. Volume is available through this helper (and is what Task 5/6/8 and the bench use for verification).

## Task 5: Assemble curved fragments

**Files**

- Create: `include/geom/brep/assemble_curved.hpp`
- Create: `tests/unit/brep_assemble_curved_test.cpp`
- Modify: `tests/CMakeLists.txt`

**Interfaces**

```cpp
struct CurvedAssembledFace {
  SurfaceGeometry surface;
  std::vector<std::vector<Point3d>> loops;      // outer first
  std::vector<std::vector<Point2d>> loops_uv;   // same sizes as loops
  std::vector<std::vector<std::optional<CurvedArcSpan>>> arc_spans;  // same shape; may be empty
};

Result<Brep> assemble_curved_solid(std::span<const CurvedAssembledFace> faces,
                                   double tolerance);
```

Algorithm:

1. Validate: at least one face; every loop has >= 3 points; `loops.size() == loops_uv.size()`; matching sizes; finite points; tolerance finite and > 0.
2. Split every loop into pieces: a maximal run of segments whose `arc_spans` entries are present and lie on the same circle (center, radius and plane normal within `tolerance`, allowing a different angular phase) becomes one arc piece; every other segment is a straight piece. The run's signed sweeps must be continuous at the shared junctions (the end angle of one segment equals the start angle of the next); a reversed or discontinuous run is `invalid_topology "assemble_curved_solid: arc run is inconsistent"`. Piece endpoints are the only vertices.
3. Weld piece endpoints within `tolerance` (grid with cell size `2 * tolerance`, plus neighbour-cell search). Reject a weld that merges endpoints farther than `tolerance`.
4. Build each undirected edge key from the welded vertex ids. Every edge must have exactly two uses, and both uses must be the same kind: arcs on the same circle, or straight segments. A closed arc (start vertex == end vertex, e.g. a full rim circle) is legal: its key is `(v, v)` and it must still be used twice. Otherwise `invalid_topology "assemble_curved_solid: edge is used once"` / `"shared edge curve mismatch"`.
5. Create one edge curve per key in the direction of the first use. An arc run becomes one `CircleCurve3` whose center/radius come from the span. Anchor `x_axis` at the run's start radial (unit, in the span's circle plane) and set `y_axis = unit(cross(plane_normal, x_axis))`; when the run's signed sweep is negative, negate `y_axis` instead of swapping the axes (a swap would move the start point). The edge stores `t0 = 0` and `t1 = |sweep| > 0`, with a closed full rim forced to exactly `2*pi`; increasing parameter then follows the run's signed sweep and `eval(t0)` is the welded start vertex. Storing the exact circle (rather than a NURBS approximation) preserves the `face_loops_uv` arc provenance and keeps `curved_solid_volume` exact. A straight piece uses `make_polyline_curve({a, b})`.
6. For each use, sample the edge at `n = clamp(ceil(span / (2*acos(1 - min(tolerance, 1e-6)/r))), 32, 4096)` angles from `t0` to `t1` (span = `t1 - t0`, r = circle radius; for a straight piece use its chord length and two endpoints), project each sample with `detail::face_uv_project` and unwrap a periodic u against the previous sample with `detail::face_uv_unwrap`; store the dense pcurve in edge direction `v0 -> v1` and record a backwards use with `coedge.reversed` (the planar assembler convention).
7. Deduplicate surfaces by exact `SurfaceGeometry` comparison within `tolerance` (plane normal/offset, cylinder frame/radius, cone, sphere, torus). Reuse the surface id.
8. `SolidEditor`: add vertices, edges (pcurve = uv polyline), faces, loops (`outer = index 0`), one shell with all faces and `closed = true`, one solid. Commit.
9. Run `validate_topology` on the result; any finding becomes an `invalid_topology` error carrying the first finding's entity index and message.

**Steps**

1. Exact-arc test: build `make_cylinder(2, 3)`, read its three faces with `face_loops_uv` (the lateral face is a full-period rectangle, the caps carry circular provenance), and assemble them back into a solid; expect `validate_topology(result).findings.empty()`, exactly 3 faces, `euler_characteristic == 2`, both rim edges stored as `CircleCurve3` with `|t1 - t0| == 2*pi`, every pcurve well-formed and closed, and `detail::curved_solid_volume(result, tol)` within `1e-9` of `pi * r^2 * h` (arcs stay exact).
2. Chordal test: assemble the five faces of `make_quarter_cylinder` (its 16-chord `LineCurve3` caps carry no provenance, but the lateral face is a smooth cylindrical patch, so the enclosed volume is the smooth-lateral value, not the prism value); expect a valid solid and `curved_solid_volume` within `1e-9` of `h * (pi*r^2/6 + n*r^2*sin(pi/(2*n))/6)`, and within `2e-3` relative of the prism value `h * (n/2) * r^2 * sin(pi/(2*n))`, `n = 16` (the difference is the lateral bulge).
3. Add a mismatch test: two faces sharing an edge with incompatible curves (arc vs chordal polyline, or two different circles) -> `invalid_topology`, message contains `"curve mismatch"`.
4. Add a single-use edge test -> `"edge is used once"`.
5. RED; implement; Debug/Release.
6. Commit: `feat(m3b2): assemble curved fragments into a Brep`.

## Task 6: Split a solid by a plane

**Files**

- Create: `include/geom/brep/split_solid.hpp`
- Create: `include/geom/brep/detail/arc_heal.hpp`
- Create: `tests/unit/brep_split_solid_test.cpp`
- Modify: `tests/CMakeLists.txt`

**Interfaces**

```cpp
struct SplitSolidResult {
  Brep below;   // negative side of the plane normal
  Brep above;   // positive side
};
Result<SplitSolidResult> split_solid_by_plane(const Brep& brep, SolidId solid,
                                              const Plane3& plane, double tolerance);
```

Algorithm:

1. For every face of the solid, call `face_loops_uv`; build the face patch (`SsiPatch{surface, uv_bounds}`) from the uv loop bounds (expanded by `tolerance`).
2. Build the plane patch: frame from `make_plane_frame(plane)`, uv box covering the solid's bounding box corners projected into that frame.
3. `ssi(face_patch, plane_patch, options)`. A `tangent_contact` or `coincident` issue returns `unsupported "split_solid_by_plane: plane is tangent or coincident with a face"`. Map each returned segment's `ua/va` to a `FaceCut` for that face.
4. `split_face_curved` for every face; classify each fragment by `surface_eval` at its outer-loop representative point (first ear centroid, as in Task 4). Side = `dot(normal, p) - offset` sign.
4b. Collapse known circular cut chains into exact arcs (`detail::collapse_circular_chains` in `include/geom/brep/detail/arc_heal.hpp`). A cut circle is known analytically, without fitting: for a `Cylinder3` face whose axis is parallel to the plane normal (`|dot(n, axis)| >= 1 - 1e-12`), the centre is the intersection of the axis line with the plane and the radius is the cylinder radius; for a `Sphere3` face with centre distance `d = dot(n, centre) - offset` and `|d| < radius`, the centre is the projection of the sphere centre and the radius is `sqrt(radius^2 - d^2)`. `Cone3`, `Torus3` and `NurbsSurface` faces produce no circle. For every fragment loop and cap loop, scan maximal runs of at least two consecutive segments whose vertices all lie on one of those circles within `max(10*tolerance, 1e-9*radius)`; replace an accepted run by its two endpoint vertices plus one `CurvedArcSpan` (`CircleCurve3` with the plane normal as axis and the start radial as `x_axis`, `t0/t1` from the endpoints, orientation following the loop). If the run closes on itself (a full turn), split it at the vertex nearest the half-turn into two spans shorter than a half turn: `assemble_curved_solid` rejects a single `2*pi` span as direction-ambiguous. Points that do not lie on a known circle stay polylines (box cuts, torus and cone cuts). SSI vertices are Newton-refined onto both surfaces, so a cylinder cut perpendicular to its axis is exact; this is what makes `pi*r^2*h/8` attainable. A sphere cut is exact in 3D but its uv image is curved, so the parameter-domain region (and therefore `curved_solid_volume` on that face) stays chord-limited; ADR-0025 records that exact volume is only claimed for cylinder cuts perpendicular to the axis plus their planar caps.
5. Build the cut cap: project every returned cut polyline onto the plane frame, weld and stitch the pieces into closed loops; reject an unclosed loop with `numeric_failure "split_solid_by_plane: cut loops do not close"`. Each side gets one planar cap face on the same plane with its own outward normal: `+n` for the negative side and `-n` for the positive side (so its plane offset flips sign too), with the loops projected into that normal's uv frame and kept with the outer loop CCW. Heal the cap loop with the same step-4b helper so the cap and the trimmed lateral face share identical exact arc edges.
6. Faces with no cut are assigned whole to one side; the two cap faces are assigned one per side.
7. Assemble each side with `assemble_curved_solid`.
8. If either side has no faces inside the plane (no fragments and no cap area), return `degenerate "split_solid_by_plane: plane does not split the solid"`.

**Steps**

1. Test the quarter cylinder split at `z = 1`: two solids, `validate_topology` clean, each `detail::curved_solid_volume` equal to `pi * r^2 * h / 8` within `1e-9` (exact arcs, so the analytic half is attainable). Across both results, the lateral/cap boundary on `z = 1` must consist of `CircleCurve3` edges (not polylines), each used exactly twice.
2. Test a box split at `z = 0.5` (all-planar integration path), volumes 4 and 4 within `1e-12`.
3. Test a full cylinder split at `z = h/2`: two solids, each volume `pi*r^2*h/2` within `1e-9`, and no `2*pi` arc span anywhere (the full cut circle must be split into two spans per loop).
4. Test tangent plane at `z = 2` (touching the top face) -> `unsupported` with the documented message.
5. RED; implement; Debug/Release.
6. Commit: `feat(m3b2): split solids by a plane through trimmed faces`.

## Task 7: Point-in-tessellation classifier

**Files**

- Create: `include/geom/brep/mesh_inside.hpp`
- Create: `tests/unit/brep_mesh_inside_test.cpp`
- Modify: `tests/CMakeLists.txt`

**Interfaces**

```cpp
// Ray-parity test over a closed triangulation. Returns true when `p` is inside.
// If `p` is within `tolerance` of the mesh surface, or every probe ray is
// degenerate, returns numeric_failure instead of guessing.
Result<bool> point_in_tessellation(const TriMesh& mesh, const Point3d& p,
                                   double tolerance);
```

Algorithm: four direction vectors (`+x`, `+y`, `+z`, normalized `(1,1,1)`). For each, count proper ray-triangle crossings with Moller-Trumbore; a hit with barycentric coordinate within `tolerance` of an edge/vertex marks the ray ambiguous and the next direction is tried. Odd/even parity decides. All four ambiguous -> `numeric_failure "point_in_tessellation: every probe ray is degenerate"`. Non-finite point, empty mesh, or non-positive tolerance -> `non_finite` / `degenerate` / `tolerance_violation`.

**Steps**

1. Build a unit cube with `tessellate(make_box(...))`; test the center (true), `(2, 0, 0)` (false), and `(0.5, 0.5, 0.5 + 1e-12)` at tol `1e-9` (numeric_failure).
2. Add an edge-grazing direction test: point whose `+x` ray hits the shared diagonal of two triangles; the classifier must fall through to another direction and still answer correctly.
3. RED; implement; Debug/Release.
4. Commit: `feat(m3b2): classify points against a tessellated solid`.

## Task 8: Curved boolean driver

**Files**

- Create: `include/geom/brep/boolean_solid_curved.hpp`
- Create: `tests/unit/brep_boolean_curved_test.cpp`
- Modify: `tests/CMakeLists.txt`

**Interfaces**

```cpp
struct CurvedBooleanOptions {
  double tolerance = 1e-9;
};

struct CurvedBooleanReport {
  std::uint32_t face_pairs = 0;
  std::uint32_t cuts = 0;
  std::uint32_t fragments = 0;
  std::uint32_t kept = 0;
  std::uint32_t dropped = 0;
  std::uint32_t flipped = 0;
  double max_cut_residual = 0.0;
};

struct CurvedBooleanResult {
  Brep brep;
  CurvedBooleanReport report;
};

Result<CurvedBooleanResult> boolean_solid_curved(
    const Brep& a, SolidId a_solid, const Brep& b, SolidId b_solid,
    SolidBooleanOp op, const CurvedBooleanOptions& options = {});
```

Algorithm:

1. Validate both operands: exactly one solid each; closed shells; `validate_topology` clean; every surface is `Plane3`, `Cylinder3`, `Sphere3` or `Torus3`; `NurbsSurface` -> `unsupported "boolean_solid_curved: NURBS operands are not supported"`; `Cone3` -> `unsupported "boolean_solid_curved: cone operands are not supported (apex-touching loops are deferred)"`. Tolerance must be finite and > 0.
2. Compute each face's 3D AABB by sampling its uv domain on a 5x5 grid. Pair faces whose AABBs overlap within `tolerance`.
3. For each pair, build patches from the uv loops and call `ssi`. `coincident` -> `unsupported "boolean_solid_curved: coplanar face contact is not supported yet"`. `tangent_contact` -> `unsupported "boolean_solid_curved: tangent face contact is not supported yet"`. `degenerate_patch` -> skip the pair. `no_convergence` -> `numeric_failure "boolean_solid_curved: cut curve did not converge"`.
4. Convert each returned segment into a cut for the A side (`ua/va`) and the B side (`ub/vb`); merge cuts per face and drop segments shorter than `tolerance`.
5. `split_face_curved` for every face of both operands. A face with no cut is one whole fragment.
6. For each A fragment: representative point at the first-ear centroid; if within `tolerance` of B's tessellation surface -> `numeric_failure "boolean_solid_curved: fragment representative lies on the other solid"`; else `point_in_tessellation` gives `in_b`. Keep when `solid_boolean_member(true, in_b, op)`.
7. For each B fragment: classify against A; keep when `solid_boolean_member(in_a, true, op)` and set `reversed = true` so its normals point out of the result.
8. If nothing is kept -> `numeric_failure "boolean_solid_curved: result is empty"`.
9. Assemble; if assembly reports more than one component -> `unsupported "boolean_solid_curved: result is not a single connected component"`.
10. Validate the result and fill the report counters deterministically (faces sorted by operand then index).

**Steps**

1. Test `split_solid_by_plane` as the reference: `difference` of the quarter cylinder with a box occupying `z >= 1` must equal the lower half within `1e-6` relative volume (both sides measured with `detail::curved_solid_volume`) and the same face count.
2. Test union of two overlapping boxes via the curved driver (all-planar path) and compare volume with `boolean_solid`; assert same volume within `1e-9` and no `validate_topology` findings.
3. Test NURBS operand -> `unsupported`.
4. Test cone operand (`make_cone`) -> `unsupported "boolean_solid_curved: cone operands are not supported (apex-touching loops are deferred)"`.
5. Test coplanar face contact: two boxes sharing a face -> `unsupported`.
6. Test disjoint union -> `unsupported "not a single connected component"`.
7. RED; implement; Debug/Release.
8. Commit: `feat(m3b2): boolean analytic B-rep solids with trimmed faces`.

## Task 9: QA, bench, docs

**Files**

- Create: `bench/bench_boolean_curved.cpp`
- Modify: `CMakeLists.txt` (`GEOM_BUILD_BENCH`)
- Create: `docs/benchmarks/m3b2-curved-boolean.md`
- Create: `docs/adr/0025-curved-boolean-analytic.md`
- Modify: `docs/architecture/geometry-library-architecture.md` §16.3
- Modify: this plan (tick completed steps)

**Steps**

1. Bench: a quarter-cylinder and a box with 5000 lateral samples; measure `boolean_solid_curved` union/intersection/difference with `detail::curved_solid_volume` verification. Report min/median/max over 5 runs. Target: median <= 2 s per operation; record FAIL honestly. Exit non-zero on validation failure, not on slowness.
2. QA corpus in `brep_boolean_curved_test.cpp`: twelve cases covering cylinder-box, sphere-box, torus-box, seam cuts, boundary cuts, tangency rejection, coplanar rejection, NURBS rejection, and cone rejection (documented `unsupported`); each result must be either a validated solid whose `detail::curved_solid_volume` matches the independently expected volume within `1e-6` relative, or a diagnosed error.
3. ADR-0025 records: analytic-only scope, sampled closed edges, no coplanar contact, cone operands deferred (apex-touching loops), single-component result, seam unwrapping, classification-tolerance policy, and the fact that trimmed curved faces are not yet tessellatable or accepted by `mass_properties` (volume is available through `detail::curved_solid_volume`; trimmed tessellation and mass properties are M3c work).
4. Update §16.3 with the M3b-2 status and the measured bench number; leave QA-8 and QA-9 unchanged and do not claim M3b-2 exit until the bench and QA corpus pass.
5. Full suite Debug/Release, `ctest`, warnings clean.

## M3b-2 exit criteria

- `Pcurve2d` is a polyline everywhere; `split_edge`/`join_vertex` preserve multi-segment pcurves.
- `face_loops_uv` unrolls periodic seams and samples edges without pcurves.
- `split_face_curved` produces outer+holes fragments with correct nesting; cuts touching a boundary never create a sliver.
- `assemble_curved_solid` produces a `validate_topology`-clean `Brep` from curved fragments with pcurves in edge direction.
- `split_solid_by_plane` splits a quarter cylinder and a box with parameter-space quadrature volumes.
- Trimmed curved faces remain outside `tessellate`/`mass_properties`; the limitation and the `curved_solid_volume` substitute are documented in ADR-0025 and §16.3.
- `boolean_solid_curved` passes the twelve-case corpus with diagnosed unsupported configurations.
- Bench and QA recorded in `docs/benchmarks/m3b2-curved-boolean.md`; FAILs are reported as FAILs.
- No regression: full Debug/Release suites green, 0 warnings.

## Self-review checklist (to run before hand-off)

- [ ] Every task has Files, Interfaces, Steps, and an explicit RED command.
- [ ] Every new public function has at least one failure-path test and one boundary test.
- [ ] The five Review Focus probes are each pinned by a named test.
- [ ] The three deviations are documented in ADR-0025 and the plan.
- [ ] No task depends on a function introduced by a later task.
- [ ] No placeholder, TODO, or "fix later" text remains.
