# M3b-2b: exact circular arcs in curved boolean results

> Addendum to `2026-11-11-m3b2-curved-boolean.md`. User decision (m11256, refined):
> preserve exact circular edges. This addendum replaces the "all boundary edges
> become polylines" rule of Task 5 and adds the two tasks needed to carry arc
> provenance from `face_loops_uv` to assembly.

## Scope

Deliver exact circular arcs for every result edge that follows an original
`CircleCurve3`, and dense pcurves for those arcs. Parameter-domain trimming keeps
using the linear arrangement (Task 4): the exact circle/cut intersections are
inserted into the boundary loops first, so the linear arrangement splits at the
true crossing points. The chord between inserted points is only used to decide
containment/topology; the emitted edge geometry is the exact arc. Residual risk:
a cut that passes entirely inside the thin lens between a chord and its arc is
resolved by the chord, not the arc; this is documented as the topology
resolution limit (it does not affect edge geometry or volume of the emitted
solid).

Out of scope: a curved 2D arrangement with algebraic intersections (M3c), cone
operands (still unsupported), NURBS operands.

## Task 4c: arc provenance in face loops

Files: `include/geom/brep/face_uv.hpp`, `tests/unit/brep_face_uv_test.cpp`.

```cpp
struct FaceUvArcRef {
  CircleCurve3 circle;   // full circle the source edge lies on
  double t0 = 0.0;       // directed start angle of the segment
  double t1 = 0.0;       // directed end angle; sweep = t1 - t0 in (-pi, pi]
};

struct FaceUvLoop {
  std::vector<Point2d> points;
  bool outer = false;
  // Same size as points: arcs[i] describes points[i] -> points[(i+1) % n]
  // (the last entry is the closing segment). Engaged only for segments that
  // come from a CircleCurve3 coedge.
  std::vector<std::optional<FaceUvArcRef>> arcs;
};
```

`face_loops_uv` fills `arcs` per segment for coedges whose curve is a
`CircleCurve3` (t = atan2 in the circle frame; the segment's directed sweep is
taken from the coedge orientation) and leaves straight, NURBS and sampled
segments disengaged. A loop whose size disagrees is `invalid_topology`; a
segment whose sweep is ambiguous (a two-point semicircle pcurve, or a chord and
tangent whose directions disagree) is `numeric_failure "face_loops_uv:
circular arc direction is ambiguous"` rather than being guessed.

Junction semantics: when a duplicate point is removed, its two neighbouring
segments are spliced and the zero-length segment is dropped; the closing point
is popped but the closing segment is kept. Two legitimately different circles
meeting at one point (the torus seam) are a **corner between arcs**, not an
error; circle-mismatch detection lives in the assembler (Task 5), not here.

Tests: a cylinder cap loop has engaged refs with radius/center equal to the
edge and increasing t; the cylinder lateral loop has engaged refs only on the
bottom and top rim runs; a box face has none; the half-cylinder fixture reports
17-point semicircle rims whose spans are exact half turns; a half-disk
(`CircleCurve3` with negative `t0`) normalizes every t into `[0, 2*pi)` and its
straight diameter segment stays disengaged; the quarter-cylinder fixture stores
its arcs as 16 `LineCurve3` chords, so its loop is entirely disengaged
(regression).

## Task 4d: exact arc spans in fragments

Files: `include/geom/brep/trim_face.hpp`, `tests/unit/brep_trim_face_test.cpp`.

```cpp
struct CurvedArcSpan {
  CircleCurve3 circle;
  double t0 = 0.0;       // 0 <= t0 < t1 <= t0 + 2*pi
  double t1 = 0.0;
};

struct CurvedFragment {
  std::vector<std::vector<Point3d>> loops;
  std::vector<std::vector<Point2d>> loops_uv;
  bool reversed = false;
  // arc_spans[l][i] describes segment i -> (i+1) % points of loop l.
  std::vector<std::vector<std::optional<CurvedArcSpan>>> arc_spans;
};
```

Algorithm:

1. For every face loop and every boundary segment whose source is a circle
   (`arcs[i]` and `arcs[i+1]` on the same circle), intersect that circle exactly
   with every cut segment (3D points from `surface_eval` of the cut uv). Solve
   `|p0 + s*d - c|^2 = r^2` in a scaled frame (`m = max(|p0-c|, r)`), keep
   `s in [0,1]`, and require the point to be within `tolerance` of the circle.
2. Insert each exact intersection into the loop as a new point: 3D point on the
   circle, uv from `face_uv_project` of the face surface; record its arc ref
   `(circle, t)`. Split the cut polyline at the same 3D parameter so the cut
   ends exactly on the arc.
3. Run the existing linear arrangement and containment-tree split on the
   augmented loops and cuts. Map every output ring vertex back to its loop
   point by exact 3D coordinates (bitwise; arrangement copies input vertices)
   and read its arc ref.
4. For each ring segment whose two endpoints carry refs on the same circle
   (center/axes/radius equal within `tolerance`), and whose traversal follows
   the circle (the ring direction fixes t0 -> t1; add 2*pi when t decreases),
   set `arc_spans[l][i]`. Segments involving a cut endpoint keep no span unless
   that endpoint was an inserted intersection (then it has a ref).
5. Invariant tests: every engaged span's endpoints satisfy
   `|p - circle.eval(t)| <= tolerance`; spans never overlap and never exceed
   2*pi; a full quarter-cylinder lateral loop split by a plane at z = h/2 has
   the bottom and top arcs engaged and the two cut lines not engaged.

Tests: cylinder split by the plane `x = 0` (lateral halves keep two `pi` spans
each, cap halves one `pi` span each, complementary); cylinder split at
`z = h/2` (full-circle spans); cap disk cut by a chord (major/minor spans); a
cut that touches an arc at one point (single inserted vertex); a polygonal face
has no spans; the chordal quarter-cylinder fixture has no spans; mismatch of
circle/tolerance is `invalid_topology` in the assembler.

## Task 4e: arc-exact planar area for volume

Files: `include/geom/brep/detail/curved_volume.hpp`,
`tests/unit/brep_curved_volume_test.cpp`.

The quadrature integrates over the chord polygon, so a planar face bounded by a
circle loses the circular segments (a 32-gon cap is 0.64% small). For a face
whose surface is `Plane3`, replace the triangle-sum area by the exact signed
area of its loops using Green's theorem per segment:

- line segment: `0.5 * cross(p_i, p_{i+1})` in uv;
- circular arc (`arcs[i]`/`arcs[i+1]` engaged on the same circle **and** the
  circle's plane is parallel to the face plane within `tolerance`): with the
  circle center projected into uv and `delta = t1 - t0` (add `2*pi` once when
  needed),
  `0.5 * (r^2 * delta + r * (cx*(sin t1 - sin t0) - cy*(cos t1 - cos t0)))`.

Sum the loop areas with their orientation (outer CCW positive, holes CW
negative). `area_uv` becomes `|sum|` for the planar path; curved faces keep the
existing subdivided quadrature. This makes cylinder caps exact while lateral
rims (circle plane perpendicular to the face) stay line segments and are
unaffected.

Tests: `make_cylinder(2,3)` volume equals `pi*r^2*h` within 1e-9 relative; a
half cylinder `x >= 0` built by splitting equals `pi*r^2*h/2` within 1e-9; box
unchanged; a cap with a hole (two coaxial circle loops) equals
`pi*(R^2-r^2)*h` within 1e-9; the chordal quarter fixture keeps its polygonal
volume (no arc provenance).

## Task 5 revision: exact arcs in assembly

`assemble_curved_solid` uses `arc_spans` when present:

1. Both uses of an edge must agree: same circle within `tolerance`, spans equal
   up to reversal. Otherwise `invalid_topology "assemble_curved_solid: arc span
   mismatch"`.
2. Build the edge curve as an exact `CircleCurve3` on the span's circle: `x_axis`
   is the unit radial at the start angle and `y_axis = cross(plane_normal,
   x_axis)`, with the edge's `t0/t1` chosen so increasing parameter follows the
   run's signed sweep (`t1 > t0`; swap the axes and negate the angles when the
   sweep is negative; a closed full rim keeps `t1 - t0 == 2*pi`). Storing the
   exact circle rather than a NURBS approximation preserves `face_loops_uv` arc
   provenance, so `curved_solid_volume` stays exact instead of falling back to
   the sampled pcurve. Endpoints replace the welded vertices' positions exactly.
3. Generate both pcurves by sampling the arc at
   `n = clamp(ceil(span / (2*acos(1 - min(tolerance, 1e-6)/r))), 32, 4096)`
   angles from `t0` to `t1` and projecting each into the adjacent face's surface
   with `face_uv_project`; unwrap periodic u against the previous sample. The
   resulting pcurve is dense enough for parameter-domain trimming and area
   estimation.
4. All other edges keep `make_polyline_curve` and the existing uv polyline.

Tests: a full cylinder split at `z = h/2` then assembled; every arc edge is a
`CircleCurve3` with `|t1 - t0|` equal to its span and endpoints exactly on the
original circle; `curved_solid_volume` of each half equals `pi*r^2*h/2` within
1e-9 relative (Task 4e makes the caps exact); a mismatched span pair is
`invalid_topology`.
