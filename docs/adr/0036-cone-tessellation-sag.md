# ADR-0036: Cone chordal sag from the surface

Status: accepted for M3c-3c-4b3a (2026-11-20); supersedes the trimmed-cone
tessellation limitation in ADR-0035.

## Context

M3c-3c-4b2 (ADR-0035) opened whole-cone/plane booleans and `split_solid_by_plane`
started producing trimmed cone frusta, but those faces still could not tessellate.
The canonical `z = 1` lower frustum (cone radius `1`, height `2`) builds 121 strip
triangles and then fails with `Errc::numeric_failure` and
`"tessellate: trimmed face refinement budget is exhausted"` (context `(1,0)`).

The trimmed refiner measured chord sag in parameter space: it evaluated
`surface_middle = face_uv_eval(projection, middle_uv)` at the mean of the endpoint
parameters and compared that point with the 3D chord midpoint. The two rims of a
frustum are sampled independently (72 samples at the base and 51 at the `z = 1`
cut), so an interior chord can join vertices on different `v` rows. Its parameter
midpoint is the mean of those two rows, which is not the chord's surface point on a
cone. The metric therefore reports a large error for a chord that is already close
to the surface, and the refiner keeps splitting until the budget is exhausted.

The probe measured, at `TessellationOptions{1e-3, 0.35}` on the `z = 1` lower
frustum:

- maximum parameter sag at the final check: `9.1028e-3`, i.e. `9.1x` the `1e-3`
  deflection;
- distance from the same triangle centroids (the points the final check measures)
  to the cone surface at the same azimuth and
  axial height: `1.09421e-3`, i.e. `1.09x` the deflection and about `8.3x` smaller;
- the parameter predicate splits `1032` edges, exactly the `8*121 + 64` budget, and
  still fails;
- the surface predicate splits `13` edges and succeeds;
- the final parameter check still reports `6.54071e-3` after the surface-based
  refinement, so the acceptance check needed the same change, not only the split
  predicate.

The other surface kinds keep their existing parameter metrics; the suites pin their
tessellation counts and positions as unchanged, and changing them would be an
unrelated regression. Only the cone showed the mismatch between "mean parameter"
and "chord midpoint".

## Decision

- **Surface distance, cone only.** Add
  `detail::tessellate_cone_radial_sag_square(const FaceUvProjection&, const Point3d&,
  double hint_u)` (`include/geom/brep/tessellate.hpp:1488`). It returns `nullopt` when
  the projection has no surface, the surface is not a `Cone3`, or
  `face_uv_project(projection, point, hint_u, false)` fails. Otherwise it returns the
  squared 3D distance from `point` to `face_uv_eval(projection, uv)`. The cone
  projection keeps the point's azimuth and axial height, so the evaluated point is
  the cone point at those coordinates. This is not necessarily the exact Euclidean
  nearest point; the distance to it is an upper bound on the distance to the nearest
  surface point. The bound is conservative: a mesh accepted with this value is at
  least as close as the number says, and the predicate may split slightly more than
  a true nearest-point metric would. It never under-splits or accepts a mesh that the
  nearest-surface check would reject.
- **Split decision.** In `tessellate_refine_trimmed_edges`, the existing parameter
  sag is still computed first, then overridden when the helper returns a value
  (`tessellate.hpp:1602-1605`). The split location remains `surface_middle` at the
  parameter midpoint, not the projection of `chord_middle`; only the decision to
  split changes. The limit remains `sag_limit = 0.9 * deflection`
  (`tessellate.hpp:1514`), and the split budget remains `8*triangles + 64`
  (`tessellate.hpp:1543`). Boundary edges are still never split.
- **Final acceptance.** The final check in `tessellate_trimmed_face` keeps the
  parameter centroid sag and overrides it with
  `tessellate_cone_radial_sag_square(projection, centroid, centroid_uv[0])`
  (`tessellate.hpp:2039-2042`). The rejection threshold remains `deflection`
  (`tessellate.hpp:2043`), with the same `Errc::unsupported` and the same message
  `"tessellate: trimmed face does not follow the surface"`.
- **Scope is bitwise unchanged elsewhere.** Non-`Cone3` surfaces, a null surface and
  a failed projection all keep the parameter metric. No signature of `tessellate`,
  `TessellationOptions` or `FaceUvProjection` changes. QA-13 keeps the cut cylinder at
  894 triangles / 449 positions at `{2e-3, 0.35}`, and the canonical box, cylinder
  and torus keep 12/8, 88/46 and 1152/576.
- **The strip zipper is sufficient for frusta.** No stitching code is added. The
  strip builder already samples each rim independently and `add_zipper`
  (`tessellate.hpp:1281`) merges the bottom and top chains by the smaller next `u`,
  emitting `nb + nt - 2` triangles (`tessellate.hpp:1285`). The `72`/`51` rim sample
  counts therefore meet with no T-junction, and the rim vertices remain the exact
  vertices the neighbouring faces sampled. The only thing blocking the frustum was
  the sag metric.
- **Apex-bearing fragments keep the canonical fan.** `sample_curved_face` tries the
  analytic builder first; a cone face goes to `tessellate_cone_face`
  (`tessellate.hpp:717`), which fans the base ring to one apex vertex. That path
  still handles the whole cone and every boolean fragment that retains the apex, and
  the measured matrix shows exactly one apex vertex for those cases. If an apex ring
  ever reached the trimmed path, it is not rectilinear and the existing gate rejects
  it with `"tessellate: trimmed curved face boundary is not rectilinear"`
  (`tessellate.hpp:1771-1774`). This slice does not rewrite the fan, the seam or the
  uv ring; that work stays deferred to M3c-3c-4b3.
- **No new tolerance.** The metric is an exact double computation, not an epsilon.
  The deflection budget is still the only tolerance: `0.9 * deflection` for
  refinement and `deflection` for final acceptance.

## Evidence

Fixture: `make_cone({0,0,0}, {0,0,1}, 1, 2, 1e-9)`, boolean tolerance `1e-9`,
`TessellationOptions{1e-3, 0.35}`. The `z = 1` case cuts the cone into the upper
small cone (apex at `z = 2`) and the lower frustum; the box is
`{-2,-2,1}..{2,2,4}`, and the union uses `{-2,-2,1}..{2,2,3}`. The translated cone
sits at `(1,2,3)` and is split at `z = 4`. `V`/`T` are positions/triangles. Volume
error is mesh volume minus the analytic value, in absolute units.

| case | V | T | mesh volume | analytic | volume error | apex vertices |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| QA-22 intersection | 51 | 98 | `0.261110903259` | `pi/12` | `-6.885e-4` | 1 |
| QA-22 difference | 134 | 264 | `1.830287422577` | `7*pi/12` | `-2.308e-3` | 0 |
| QA-22 union (low box) | 142 | 280 | `33.830287422577` | `32 + 7*pi/12` | `-2.308e-3` | 0 |
| split `z = 1`, below | 134 | 264 | `1.830287422577` | `7*pi/12` | `-2.308e-3` | 0 |
| split `z = 1`, above | 51 | 98 | `0.261110903259` | `pi/12` | `-6.885e-4` | 1 |
| split `z = 0.5`, below | 140 | 276 | `1.209417839998` | `pi*0.5/3*(1 + 0.75 + 0.75^2)` | `-1.404e-3` | 0 |
| split `z = 1.5`, below | 124 | 244 | `2.058752793525` | `pi*1.5/3*(1 + 0.25 + 0.25^2)` | `-2.917e-3` | 0 |
| translated split `z = 4`, below | 134 | 264 | `1.830287422577` | `7*pi/12` | `-2.308e-3` | 0 |
| translated split `z = 4`, above | 51 | 98 | `0.261110903259` | `pi/12` | `-6.885e-4` | 1 |
| translated boolean difference | 134 | 264 | `1.830287422577` | `7*pi/12` | `-2.308e-3` | 0 |
| apex intersection `z = 0.5` | 62 | 120 | `0.882011366069` | `0.28125*pi` | `-1.562e-3` | 1 |
| apex intersection `z = 1.5` | 37 | 70 | `0.032559033313` | `0.0104166667*pi` | `-1.659e-4` | 1 |
| whole cone | 72 | 140 | `2.091662475470` | `2*pi/3` | `-2.733e-3` | 1 |

Every row has 0 boundary edges, 0 non-manifold edges, `triangle_face.size() ==
triangles.size()` with every id inside `face_count()`, and a second tessellation that
is bitwise identical in positions, triangles and `triangle_face`. The apex-bearing
rows have exactly one apex vertex: `(0,0,2)` for the origin cases and `(1,2,5)` for
the translated upper cone. The `z = 0.5` cap is an ear-clipped polygon with no vertex
at `(0,0,0.5)`; that is the measured geometry, and the unique apex is the cone's
analytic `(0,0,2)`.

The maximum absolute volume error is `2.917e-3` (split `z = 1.5`, below). The
largest relative error is `5.07e-3` at apex intersection `z = 1.5`
(`1.659e-4 / 0.0327249`), because that volume is small; the case is pinned by its
12-decimal literal at `1e-9` relative, not by the `5e-3` band. Task 1's five-row
QA-22 census uses `doctest::Approx(...).epsilon(5e-3)` and its largest relative
error is `2.63e-3` at the intersection. The original brief and plan text call `2.917e-3` the
"maximum relative error"; the measured probe column is absolute, so this ADR
records both quantities rather than repeating the label.

Focused matrix (four cases: the Task 1 trimmed-frusta case and the three Task 2
cases): 2,682 assertions pass in Debug and Release. Full suites: Debug and Release
each report 1013 cases / 673,296 assertions, 0 failed, zero warnings under
`/W4 /WX /utf-8` with `GEOM_WERROR=ON`. Task 1's baseline was 1010 / 671,694; Task 2
added exactly 3 cases and 1,602 assertions.

## Consequences

- Trimmed cone frusta produced by `boolean_solid_curved` and
  `split_solid_by_plane` now tessellate watertight at the requested deflection, with
  the counts and volumes in the table. The `72`/`51` rim pair is the ordinary case,
  not a special case.
- The cone is the only surface that takes the geometric sag path. The parameter
  metric remains for planes, cylinders, spheres, tori and every non-cone trimmed
  face, including its fallback when a cone projection fails, so those results stay
  bitwise unchanged.
- The metric is cheap and local: one projection plus one surface evaluation per
  candidate edge in the refiner and per triangle in the final check. It adds no
  state, no tolerance and no public API.
- The value is a conservative upper bound, so it can split slightly more than the
  exact nearest-surface distance would. That is intentional: the contract promises
  the deflection bound, and over-splitting is preferable to accepting a mesh whose
  true error is unknown.
- Apex-bearing fragments are still canonical-only. The apex fan, seam rewrite and
  uv-ring apex insertion remain in M3c-3c-4b3, together with oblique/axis-containing
  cone planes, general curved/curved cone booleans and re-centring of translated
  cones whose seam copies collide.
