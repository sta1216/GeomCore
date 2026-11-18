# M3c-3c-4b1 — Cone chart, apex projection and open-chain domain cuts

## Goal

Make the seam-only domain splitter usable on a cone lateral face: project points
at the apex, describe the cone chart (`u` periodic, `v` from the base circle to
the collapsed apex), dispatch synthetic cone faces to
`trim_face_split_domain`, and accept open cut chains whose two ends land on
domain sides (real, seam or collapsed). The boolean gate for cone operands stays
closed: no boolean operation changes result in this slice.

## Architecture

`face_loops_uv` projects a face onto its surface chart; a cone lateral face is
the canonical case where one domain side (the apex) collapses to a single 3D
point, exactly like a sphere pole but with only one pole. The splitter's
geometric core already collapses a material cycle that owns a
`TrimSideKind::kCollapsed` side into one pole point
(`trim_face_split_domain`). What is missing is upstream:

- `face_uv_project_cone` refuses every point at or beyond the apex, so the
  stored loop cannot be read and every downstream consumer inherits the error.
- `trim_face_domain` knows sphere and torus charts only; there is no cone chart.
- `split_face_curved` has no cone dispatch, so a cone face falls through to the
  generic clip path with a zero-area stored loop.
- `trim_face_split_domain` still pairs open cut pieces as seam-to-seam only; a
  meridian from the base circle to the apex is refused with
  `"open cuts on a seam-only face"`.

Downstream representation is unchanged: fragments still carry
`loops`/`loops_uv` of equal length, and the apex is one uv point on the back
seam copy (`u = 2π`) produced by collapsing the top side, matching the sphere
pole convention.

## Tech Stack

C++20, CMake 4.4.4, MSVC 19.44 (`/W4 /WX /utf-8`), doctest 2.4.11. No new
dependency. Existing helpers: `trim_face_split_domain`,
`trim_face_window_chains`, `trim_face_domain_edge_side`,
`trim_face_walk_cycle`, `Arrangement2`, `split_segments`.

## Spec

- `docs/architecture/geometry-library-architecture.md` v2.1 §8.5–§8.6 (SSI and
  arrangement-based face splitting), §16.3 status table.
- `docs/adr/0032-general-domain-splitter.md` (general splitter contract).
- Probe report `D:\ai_geometry\.superpowers\sdd\m3c3c4b-probe-report.md`
  (cone anatomy, bypassed splitter results, apex-ring analysis).

## Global Constraints

- Exact predicates (`sign_orient2d`) decide topology; floating point only
  samples. No new epsilon constants beyond the existing relative parameter
  tolerance.
- No silent empty or self-touching result: every unsupported input returns an
  `Error` with a stable message and context ids.
- The sphere and torus paths must stay bitwise identical. The canonical cone
  mesh, mass properties and tessellation must stay unchanged.
- The boolean gate for cone operands stays as it is in this slice
  (`"boolean_solid_curved: cone operands need closed-cut and apex topology;
  booleans are deferred to M3c-3c-3"`).
- TDD: test first, observe RED, minimal implementation, Debug and Release
  suites, one commit per task.
- MSBuild header dependency tracking is unreliable in this repository: after a
  header edit, `touch` one `.cpp` or build with `/t:Rebuild`.
- Test files are already registered in `tests/CMakeLists.txt`; no CMake change
  is expected.

## Review Focus

1. Apex projection: only points within a rounding-level band of the apex are
   accepted; points on the other nappe (beyond the apex) still error. The
   returned `u` is stable (`previous_u`), so the ring closes without a jump, and
   the apex appears once after consecutive-duplicate removal.
2. Cone chart: bounds `[0, 2π] × [0, radius/tan(half_angle)]`, `u` periodic,
   `v = v_lo` real, `v = v_hi` collapsed. A genuinely trimmed cone face must
   not be reopened as the full chart.
3. Dispatch: only a synthetic cone face (stored loop area zero, one loop) goes
   to `trim_face_split_domain`; a trimmed cone face keeps the existing
   `unsupported` behaviour.
4. Open chains: both ends must lie within the parameter tolerance of a domain
   side; an end strictly inside the domain is refused; a chain lying on a side
   is refused as a degenerate cut. The arrangement must produce two material
   cycles, each keeping the real boundary edge and collapsing the apex side to
   one point.
5. Fragment geometry: a horizontal closed cut on the cone yields a cap with uv
   area `π` (one apex point) and a band with area `2π`; a meridian open cut
   yields two half-charts of area `2π` each, both with one apex point.
6. Regression: sphere cap and torus offset tests untouched; boolean gate
   message unchanged. Three temporary cone refusals keep the pre-fan contracts
   until 4b2/4b3: `split_face_curved: cone lateral faces are not supported yet`,
   `curved_solid_volume: cone faces are not supported yet`, and
   `trimmed_face_boundaries: cone faces are not supported yet`. No other
   assembly, mass or volume change.

## File Map

- `include/geom/brep/face_uv.hpp` — apex-tolerant cone projection; pass
  `previous_u` through `face_uv_project`.
- `include/geom/brep/trim_domain.hpp` — cone chart in `trim_face_domain`.
- `include/geom/brep/trim_face.hpp` — cone dispatch in `split_face_curved`;
  open ends on any domain side in `trim_face_split_domain`.
- `tests/unit/brep_face_uv_test.cpp` — cone loop projection (rewrites the
  pinned apex-error test).
- `tests/unit/brep_trim_face_test.cpp` — cone closed cut, cone meridian, generic
  rectangle open cut, rejection of an interior end.
- `docs/adr/0034-cone-chart-open-cuts.md` — new ADR.
- `docs/architecture/geometry-library-architecture.md` — §16.3 paragraph.

## Task 1 — Apex-tolerant cone projection

**Files:** `include/geom/brep/face_uv.hpp`,
`tests/unit/brep_face_uv_test.cpp`, plus the temporary cone refusals in
`include/geom/brep/trim_face.hpp`, `include/geom/brep/detail/curved_volume.hpp`,
`include/geom/brep/trimmed_domain.hpp` and their tests.

**Interfaces:**

```cpp
inline Result<Point2d> face_uv_project_cone(const Cone3& surface,
                                            const Point3d& point,
                                            double previous_u,
                                            bool has_previous_u);
```

**Consumes:** `Cone3` (`origin`, `axis`, `x_axis`, `radius`, `half_angle`),
`face_uv_project` dispatch, `kFaceUvPi`.

**Produces:** apex accepted as `{has_previous_u ? previous_u : 0.0, v_apex}`
where `v_apex = radius / tan(half_angle)`; the band is symmetric
(`|v - v_apex| <= 1e-9*max(1, v_apex)`) and only points strictly above it
(`v > v_apex + tol`, on the other nappe) return
`"face_loops_uv: cone projection is undefined at or beyond the apex"`. A cone
thinner than its own band (`v_apex <= tol`) keeps the old refusal.
`face_loops_uv` on a canonical cone lateral face returns one loop.

**Steps:**

1. Add the failing test (replace the current apex-error assertion):

   ```cpp
   TEST_CASE("face loops uv: the cone apex projects to one chart point") {
     const auto cone = geom::make_cone(geom::Point3d{0, 0, 0},
                                       geom::Vec3d{0, 0, 1}, 1.0, 2.0, 1e-9);
     REQUIRE(cone.has_value());
     const geom::FaceId lateral = face_with_surface<geom::Cone3>(cone.value());
     const auto loops = geom::face_loops_uv(cone.value(), lateral, 1e-9);
     REQUIRE(loops.has_value());
     REQUIRE(loops.value().size() == 1u);
     const auto& points = loops.value().front().points;
     std::size_t apex = 0, rim = 0;
     for (const geom::Point2d& p : points) {
       CHECK(p[0] >= -1e-12);
       CHECK(p[0] <= 2.0 * geom::kFaceUvPi + 1e-12);
       if (p[1] == 2.0) ++apex;
       if (p[1] == 0.0) ++rim;
     }
     CHECK(apex == 1u);
     CHECK(rim > 1u);
   }
   ```

   Keep an explicit error test for a point beyond the apex (construct with
   `face_uv_eval` at `v = v_apex + 1e-3` and project it).

2. Run the focused suite; expect the new case to fail on the old error.
3. Implement: compute `v_apex`; accept when `v <= v_apex + tol` and
   `v >= v_apex - tol`, returning the previous azimuth and `v_apex`; otherwise
   keep the existing error. Extend the call in `face_uv_project` to pass
   `previous_u`/`has_previous_u`.
4. Run focused and full Debug/Release suites; expect zero warnings.
5. Commit `fix(face_uv): project the cone apex onto its chart point`.

## Task 2 — Cone chart and synthetic dispatch

**Files:** `include/geom/brep/trim_domain.hpp`,
`include/geom/brep/trim_face.hpp`, `tests/unit/brep_trim_face_test.cpp`.

**Interfaces:**

```cpp
if (std::holds_alternative<Cone3>(surface)) {
  domain.bounds = {0.0, 2.0 * kFaceUvPi, 0.0,
                   radius / std::tan(half_angle)};
  domain.periodic_u = true;
  domain.periodic_v = false;
  domain.sides = {kSeamCopy, kSeamCopy, kRealBoundary, kCollapsed};
  return domain;
}
```

**Consumes:** `trim_face_domain`, `detail::trim_face_split_domain`,
`split_face_curved` dispatch, `signed_area2`.

**Produces:** `split_face_curved` on a canonical cone lateral face routes to
`trim_face_split_domain` and splits closed cuts. A cone face whose stored loop
is not the full chart keeps the existing refusal.

**Steps:**

1. Add the failing test: `split_face_curved(cone, lateral, {FaceCut{u,v}})` with
   a closed horizontal cut at `v = 1` (65 points, endpoints `(0,1)` and
   `(2π,1)`) returns 2 fragments; the cap contains exactly one point at
   `v = 2`; and the absolute uv areas are `π` and `2π` (1e-12).
2. Run; expect `"split_face_curved: open cuts on a seam-only face"` or
   `"no synthetic domain"`.
3. Implement: replace the temporary Task 1 guard in `split_face_curved` with
   the cone dispatch. A canonical cone's stored slit loop has area
   `π * v_apex` (half the chart), not zero, so the synthetic test is: one
   stored loop, `| |signed_area2| - 0.5 * chart | <= 1e-9 * chart`, and the
   loop's points lie only on `v_lo`/`v_hi` with exactly one point on the
   collapsed apex row. Synthetic + cuts goes to `trim_face_split_domain` with
   empty holes; synthetic + no cuts falls through unchanged; any other cone
   returns `unsupported "split_face_curved: trimmed cone faces are not
   supported yet"` with `(face.index(), 0)`. Do not route trimmed cones
   through the generic clip path.
4. Run focused and full suites (expect 988 cases / about 668,860 assertions).
5. Commit `feat(trim_face): split closed cuts on a cone chart`.

## Task 3 — Open cut ends on any domain side

**Files:** `include/geom/brep/trim_face.hpp`,
`tests/unit/brep_trim_face_test.cpp`.

**Interfaces:** no new helper: `detail::trim_face_point_on_domain_side`
(`trim_domain.hpp:154`) already reports whether a point lies within the
parameter tolerance of any side.

**Consumes:** `trim_face_window_chains`, `trim_face_seam_only_seam_pair`,
`trim_face_domain_edge_side`, `trim_face_point_on_domain_side`,
`split_segments`.

**Produces:** an open chain whose two ends lie on domain sides becomes a
non-island chain; both ends are snapped to the side value; the arrangement
splits the synthetic side segment at the endpoint. An end strictly inside the
domain is already refused upstream by `trim_face_window_chains` with
`unsupported "trim_face: cut chain ends inside the parameter domain"`; a piece
with one end on a periodic seam and the other elsewhere is refused with
`unsupported "split_face_curved: open cut crosses the seam"` (re-origining the
chart is deferred); a chain lying on a side returns
`"split_face_curved: cut lies on the parameter domain boundary"` as today.

**Steps:**

1. Add the failing tests:

   ```cpp
   // Cone lateral face, meridian from the base rim to the apex at u = 0.7
   // (deliberately not one of the pi/32 side samples).
   std::vector<geom::Point2d> meridian;
   meridian.push_back({0.7, 0.0});
   meridian.push_back({0.7, 2.0});
   const std::vector<geom::FaceCut> cuts{geom::FaceCut{meridian}};
   const auto split = geom::split_face_curved(cone.value(), lateral, cuts, 1e-9);
   REQUIRE(split.has_value());
   REQUIRE(split.value().size() == 2u);
   // Each half chart has area 2*pi and exactly one apex point.
   ```

   Also split a rectangle domain (`{0,2π,0,2}`, sides
   `{seam,seam,real,real}`) with the meridian from `(0.7,0)` to `(0.7,2)` and
   check two fragments of area `2π`; and a rejected chain ending at `(0.7,1)`
   inside the chart (assert the upstream message).

2. Run; expect the seam-pair refusal.
3. Implement: in the pairing loop, after seam-pair matching fails, accept an
   open piece when both ends pass `trim_face_point_on_domain_side` and neither
   end is on a periodic seam side; snap each end exactly to its side
   (`u_lo/u_hi/v_lo/v_hi`). Reject pieces that mix a seam end with a non-seam
   end using the new message. Keep the seam-pair path, closed-island path and
   partition guard unchanged. Ensure the arrangement has a vertex at the cut
   endpoint: if `split_segments` does not split the side segment there, split
   the side polyline at the endpoint before building the arrangement.
4. Run focused and full suites.
5. Commit `feat(trim_face): accept open cuts that end on the domain boundary`.

## Task 4 — QA-21, ADR and the deferred contract

**Files:** `tests/unit/brep_trim_face_test.cpp`,
`tests/unit/brep_face_uv_test.cpp`,
`docs/adr/0034-cone-chart-open-cuts.md`,
`docs/architecture/geometry-library-architecture.md`, this plan.

**Steps:**

1. Add QA-21 as one test case that repeats the Task 1–3 checks on the canonical
   cone: apex loop, closed horizontal cut (`π` cap with one apex, `2π` band),
   meridian split (two `2π` halves, one apex each), and asserts every fragment
   has `loops.size() == loops_uv.size()`.
2. Write ADR-0034: cone chart, apex projection rule, open-end acceptance,
   representations (apex at the back seam copy), refusal matrix (interior end,
   side-lying chain, trimmed cone face, and the temporary
   `split_face_curved`/`curved_solid_volume`/`trimmed_face_boundaries` cone
   refusals), and what stays deferred (boolean gate, apex-fan volume, trimmed
   cone tessellation, general twisted crossings).
3. Add the §16.3 paragraph and cross-reference ADR-0032/ADR-0034; leave
   line endings as found (CRLF).
4. Run Debug and Release suites; record the exact counts.
5. Commit `docs(m3c3c4b1): record the cone chart and open-cut contract`.

## Exit criteria

- `face_loops_uv` on a canonical cone lateral face succeeds and yields one
  loop with exactly one apex point at `v = 2`.
- `split_face_curved` splits a closed horizontal cone cut into a `π` cap and a
  `2π` band, and a meridian open cut into two `2π` halves.
- Open ends strictly inside the chart, chains lying on a side, and trimmed cone
  faces are refused with the stable messages above.
- The boolean gate for cone operands and all sphere/torus/cylinder outputs are
  unchanged (bitwise regression tests already cover them).
- Debug and Release: all tests pass, zero warnings, `git status` clean after
  each task.

## Self-review

- [ ] Every task has a RED test that fails on the parent commit.
- [ ] Apex projection limit is relative and documented; no new absolute epsilon.
- [ ] Open-end detection uses the same relative parameter tolerance as the
      splitter, and endpoints are snapped so the arrangement has exact vertices.
- [ ] Assembly code is untouched; mass and volume paths only add the three
      temporary cone refusals recorded above, with no numerical change.
- [ ] ADR records the refusal matrix and the deferred boolean/volume work.
- [ ] Test counts are recorded in the exit report and in ledger.
