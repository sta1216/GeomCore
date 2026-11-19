# M3c-3c-4b2: Cone operands in the curved boolean

## Goal

Enable cone operands in `boolean_solid_curved` for plane cuts: remove the cone
gate, integrate the apex fan so a cone fragment's volume and mass are correct,
weld the duplicate apex copies produced by the seam-only splitter so the two
halves assemble, and pin the end-to-end volumes with a QA test.

The boolean gate is opened only for cone operands cut by planes. Cone/cone,
cone/torus and curved-curved pairs stay deferred.

## Architecture

- `trim_face_domain` already describes the cone chart (`u` periodic, `v` from the
  base rim to the collapsed apex at `radius / tan(half_angle)`), and
  `split_face_curved` already dispatches synthetic cone faces and accepts open
  cuts ending on real, seam or collapsed sides.
- `curved_solid_volume` refuses cone faces ("cone faces are not supported yet")
  because the apex fan is not integrated. `trimmed_face_boundaries` refuses cone
  faces for the same reason. Both guards are temporary.
- `mass_properties` already has `mass_cone_domain` for the canonical cone; the
  trimmed path refuses cone boundaries through `trimmed_face_boundaries`.
- `assemble_curved_solid` pairs rims by undirected endpoint identity. A cone
  fragment carries two uv copies of the apex (cut end and back-seam copy) that
  evaluate to the same 3D point, so the ring contains a zero-length 3D edge and
  the assembly reports `edge is used once`.
- `boolean_solid_curved` rejects cone operands at its surface gate with
  "cone operands need closed-cut and apex topology".

## Tech Stack

C++20, CMake + VS2022 (MSVC 19.44) and GCC/Clang in CI, doctest 2.4.11,
`/W4 /WX /utf-8`, exact predicates for topology decisions.

## Spec

- Architecture v2.1 §8.5–8.10 (curved boolean, mesh, volume), §16.3.
- ADR-0034 (cone chart and open cuts), ADR-0035 (this slice).

## Global Constraints

- No global epsilon; exact predicates for orientation/inside decisions.
- Whole-cone operands and plane cutters only; no cone/cone, cone/torus or
  NURBS pairs.
- Sphere and torus results must not change bitwise on existing corpora.
- Every failure returns a `Result` with a message; no silent empty result.
- MSBuild header dependency tracking is broken: touch a `.cpp` or `/t:Rebuild`
  after header edits.
- Tests are registered in `tests/CMakeLists.txt`; no production code in a docs
  task.

## Scope and deviations (approval needed)

1. Only plane cutters; `boolean_solid_curved` accepts a cone operand only when
   every face of the other operand is planar. Curved/curved pairs stay refused.
2. Only whole-cone operands (single shell, single lateral face, canonical
   seam-only chart). A cone that is one face of a larger B-rep stays refused.
3. The apex weld collapses coincident 3D apex copies inside a ring before
   pcurve creation; it does not merge the two uv copies in `loops_uv`.
4. trimmed cone tessellation (4b3) stays out of scope.

## Review Focus

- Apex fan volume: whole cone 2*pi/3, upper small cone pi/6, lower frustum
  7*pi/12 (measured values below), and the fan must not double count the apex.
- Apex weld: no zero-length 3D edge, ring sizes preserved, exact 1e-12 match,
  deterministic, and no change to plane/cylinder/sphere/torus paths.
- Boolean completeness: union/intersection/difference produce one closed shell,
  outward orientation, validate ok, and the reported volume matches the analytic
  values within 1e-9 relative.
- Refusal paths: cone/cone, cone/torus, trimmed cone, non-plane cutter.
- Cross-platform: MSVC and GCC/Clang (the apex row is `radius/tan(half_angle)`,
  not the literal height).

## File Map

- `include/geom/brep/detail/curved_volume.hpp` — remove cone guard, apex fan.
- `include/geom/brep/trimmed_domain.hpp` — cone boundary support.
- `include/geom/brep/mass_properties.hpp` — trimmed cone path.
- `include/geom/brep/detail/assemble_curved.hpp` — apex weld.
- `include/geom/brep/boolean_solid_curved.hpp` — gate, plane-only guard.
- `tests/unit/brep_curved_volume_test.cpp`, `brep_boolean_curved_test.cpp`,
  `brep_trim_face_test.cpp`, `tessellate_trimmed_test.cpp`.
- `docs/adr/0035-cone-boolean.md`, `docs/benchmarks/m3c3c4b2.md`,
  `docs/architecture/geometry-library-architecture.md` §16.3.

## Task 1: cone volume fan

Files: `include/geom/brep/detail/curved_volume.hpp`,
`include/geom/brep/trimmed_domain.hpp`,
`tests/unit/brep_curved_volume_test.cpp`.

Interfaces:
- `detail::curved_solid_volume(const Brep&, SolidId, double tolerance) -> Result<double>` (existing).
- `detail::trimmed_face_boundaries(const Brep&, FaceId, const detail::FaceUvProjection&, double) -> Result<std::vector<TrimmedBoundary>>` (existing).

Steps:
1. Write failing tests: whole cone volume `2*pi/3` (1e-12 relative), upper small
   cone from `cone - box z>=1` is Task 3, here test canonical cone only and a
   cone fragment with a plane cut at z=1 (upper part) if the splitter is
   reachable through the public API; otherwise test the canonical cone.
2. Remove the `Cone3` guard in `curved_solid_volume` and remove the cone guard in
   `trimmed_face_boundaries`.
3. Ensure the uv triangulation drops zero-area apex triangles from the fan and
   integrates the remaining fan triangles exactly once; the apex row is a single
   uv row per fragment.
4. Run Debug/Release; the canonical cone must equal `2*pi/3` and the existing
   sphere/torus volumes must be unchanged bitwise.

Acceptance: whole cone `2.0943951023931953`, upper cone `0.2617993877991494`
(after Task 3 split) and no change to existing volume tests.

## Task 2: apex weld in assembly

Files: `include/geom/brep/detail/assemble_curved.hpp`,
`tests/unit/brep_trim_face_test.cpp`, `tests/unit/brep_boolean_curved_test.cpp`.

Steps:
1. Write a failing test: split the canonical cone by the plane z=1 through
   `split_solid_by_plane` (or the curved boolean internals if the public path is
   not available) and assemble the two fragments; expect one closed shell.
2. In the harvest/assemble path, collapse consecutive 3D-coincident points
   within 1e-12 tolerance so the two apex copies become one vertex; keep
   `loops_uv` untouched and keep the ring sizes consistent.
3. Skip pcurve creation for the collapsed zero-length segment (the apex has no
   tangent); every other edge keeps its pcurve.
4. Run the two-half assembly and the existing plane/cylinder/sphere/torus
   assembly tests.

Acceptance: two-half cone assembles to one closed shell, volume `2*pi/3`,
`validate_topology` ok, no zero-length 3D edge in any ring.

## Task 3: open the cone gate for plane cuts

Files: `include/geom/brep/boolean_solid_curved.hpp`,
`include/geom/brep/split_solid.hpp`,
`tests/unit/brep_boolean_curved_test.cpp`.

Steps:
1. Failing e2e tests: cone r=1,h=2 against box `{-2,-2,1}..{2,2,4}`:
   intersection volume `0.2617993877991494` (upper small cone, apex present),
   difference `1.8325957145940461` (lower frustum), union with box
   `{-2,-2,1}..{2,2,3}` volume `32 + 1.8325957145940461` (box plus lower
   frustum). Assert one shell, closed, outward, validate ok, bitwise determinism.
2. Replace the cone gate with a plane-only guard: allow Cone3 operands when the
   other operand's faces are all `Plane3`; otherwise keep an explicit
   `unsupported "boolean_solid_curved: cone operands only support plane cutters"`
   message.
4. Add the analytic cone/plane cut circle to `split_solid_face_circle` for planes
   perpendicular to the cone axis: `rho = radius - t*tan(half_angle)` at axial
   distance `t = (offset - dot(normal, origin)) / dot(normal, axis)`; center on
   the axis, frame from `cone.x_axis` projected onto the plane. This replaces the
   inscribed SSI polyline and is what lets the 1e-12 acceptance pass.
5. For cone operands, refuse cutters whose planes intersect the cone's bounding
   sphere but are not perpendicular to the axis:
   `unsupported "boolean_solid_curved: cone cutters must be perpendicular to the
   cone axis; oblique and axis-containing planes are deferred to M3c-3c-4b3"`.
   A plane through the axis currently fails assembly with half-turn rim arcs, so
   it must be refused deterministically instead.
6. Keep cone/cone and cone/torus refused before pairing.
7. Run the full curved boolean corpus (sphere/torus/cylinder/plane) and confirm
   no result changes.

Acceptance: the three volumes above, refusals exact, existing corpus unchanged.

## Task 4: QA-22, ADR-0035, bench, §16.3

Files: `tests/unit/brep_boolean_curved_test.cpp`,
`docs/adr/0035-cone-boolean.md`, `docs/adr/0034-cone-chart-open-cuts.md`,
`docs/benchmarks/m3c3c4b1.md`, `docs/benchmarks/m3c3c4b2.md`,
`docs/architecture/geometry-library-architecture.md`,
`docs/superpowers/plans/2026-11-19-m3c3c4b2-cone-boolean.md`.

Steps:
1. Add QA-22 as a single test that runs the three operations and asserts the
   analytic volumes, closure, validate and determinism.
2. ADR-0035 records the plane-only scope, apex weld, apex fan, refusals and
   deferred cone/cone, tessellation and curved cutters.
3. Benchmark note records the measured volumes and no timing gate; §16.3 gets an
   M3c-3c-4b2 paragraph. ADR-0034 and the M3c-3c-4b1 benchmark still say the cone
   volume and trimmed-boundary paths refuse cone faces; update them to the new
   support (the apex fan and the expanded trimmed boundary).
4. Tick the plan boxes; final counts must match the measured suite.

## Exit criteria

- Debug and Release: all tests pass, zero warnings.
- Cone operands work for plane cutters with the analytic volumes above.
- Cone/cone, cone/torus, trimmed cone and non-plane cutters are refused with
  exact messages.
- No sphere/torus/plane regression.

## Self-review

- [ ] Every task has a RED test and a concrete acceptance value.
- [ ] No task changes a public signature.
- [ ] The apex conventions match ADR-0034 / QA-21.
- [ ] The cone gate message and the plane-only guard are distinct.
- [ ] Counts and document references are updated together.
