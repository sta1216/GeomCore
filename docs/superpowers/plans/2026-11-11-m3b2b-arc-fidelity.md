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
  double t = 0.0;        // angle of this point, [0, 2*pi)
};

struct FaceUvLoop {
  std::vector<Point2d> points;
  bool outer = false;
  // Same size as points; engaged for points sampled from a CircleCurve3 edge.
  std::vector<std::optional<FaceUvArcRef>> arcs;
};
```

`face_loops_uv` fills `arcs` for every coedge whose curve is a `CircleCurve3`
(t = atan2 in the circle frame, normalized), and leaves it empty for lines,
NURBS, and sampled approximations. A loop whose size disagrees is
`invalid_topology`.

Tests: a cylinder cap loop has engaged refs with radius/center equal to the
edge and strictly increasing t; a polyline edge has none; a mixed loop (arc +
line) has engaged refs only on the arc run; t is normalized to [0, 2*pi) and
matches the 3D point within 1e-12.

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

Tests: quarter-cylinder lateral face cut by `z = h/2` (arcs engaged, exact t
ends at 0, pi/2); cap face cut by a chord (arc spans on both sides); a cut that
touches an arc at one point (single inserted vertex); a polygonal face has no
spans; mismatch of circle/tolerance is `invalid_topology` in the assembler.

## Task 5 revision: exact arcs in assembly

`assemble_curved_solid` uses `arc_spans` when present:

1. Both uses of an edge must agree: same circle within `tolerance`, spans equal
   up to reversal. Otherwise `invalid_topology "assemble_curved_solid: arc span
   mismatch"`.
2. Build the edge curve as a rational quadratic NURBS arc
   (`make_circle_arc_nurbs`): split the span into pieces of at most `pi/2`,
   each piece is a quadratic Bezier with control points
   `P0 = C + r*x(t0)`, `P2 = C + r*x(t1)`,
   `P1 = C + (r / cos(theta/2)) * x((t0+t1)/2)` and weights
   `{1, cos(theta/2), 1}`; concatenate pieces into one `NurbsCurve` (degree 2,
   shared endpoints, knot multiplicity 2 at joins). Endpoints replace the
   welded vertices' positions exactly.
3. Generate both pcurves by sampling the arc at
   `n = clamp(ceil(span / (2*acos(1 - min(tolerance, 1e-6)/r))), 32, 4096)`
   points and projecting each into the adjacent face's surface with
   `face_uv_project`; unwrap periodic u against the previous sample. The
   resulting pcurve is dense enough that its chord error is <= 1e-6 for volume
   quadrature.
4. All other edges keep `make_polyline_curve` and the existing uv polyline.

Tests: quarter-cylinder fixture split at `z = h/2` then assembled; every arc
edge is a `NurbsCurve` with degree 2 and endpoints exactly on the original
circle; `curved_solid_volume` equals `pi*r^2*h/4` within `1e-9` relative (the
32-chord fixture is healed); full-cylinder cap edges stay `CircleCurve3` when
the whole edge is untouched; a mismatched span pair is `invalid_topology`.
