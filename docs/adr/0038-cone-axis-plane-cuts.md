# ADR-0038: Cutting cones with planes through the axis

Status: accepted for M3c-3c-4b3c-2 (2026-11-22); narrows the axis-plane deferral recorded in ADR-0036 and ADR-0037. Oblique cone cutters remain deferred to M3c-3c-4b3.

## Context

`boolean_solid_curved` has accepted whole cones cut by planes since M3c-3c-4b2, but only when the plane is perpendicular to the cone axis (the section is a circle) or when its orientation is harmless because the plane misses the cone's conservative bounding sphere. Planes through the axis were refused by `curved_boolean_cone_plane_cut_supported` (`include/geom/brep/boolean_solid_curved.hpp:1445-1466`) and reported by the gate as deferred to M3c-3c-4b3.

That refusal was conservative, not a representational limit. A plane containing the cone axis meets the lateral surface in two rulings: two generator lines that run from the apex to the base circle. The resulting cut on the lateral face is two meridian chains, which the existing surface/surface intersector already returns exactly. No ellipse has to be traced, no synthetic chain or vertex fan is introduced, and the split/trim/assemble path used by perpendicular cuts consumes the chains unchanged.

Two details needed explicit treatment. Whole-cone operands only: the canonical chart is `[0, 2π] x [0, radius/tan(half_angle)]`, whose apex row collapses to a point, so trimmed cone fragments and non-canonical charts remain outside the contract. At the apex the chart carries two parameter copies of the same 3D point; for an axis-plane cut those copies are the legitimate closure of the two meridians, not the pathological coincident seam copies that the placement guard rejects.

## Decision

1. Whole cones may be cut by a plane that is perpendicular to the axis or contains it. `detail::curved_boolean_cone_plane_contains_axis` (`boolean_solid_curved.hpp:1428-1443`) normalizes the plane normal and cone axis, rejects non-finite or zero-length inputs, requires `|dot(unit_normal, unit_axis)| <= kCurvedBooleanConePerpendicularSlack` (`1e-9`), and requires the cone origin to lie on the plane within the caller's `tolerance`. The origin test is a placement coincidence owned by the caller, exactly like the other placement checks in this file. `curved_boolean_cone_plane_cut_supported` (`:1445-1466`) accepts the axis case at `:1458`, after the perpendicular branch and before the bounding-sphere test.
2. The cutter must still be a plane. A cone paired with a non-planar operand is refused at `:2101-2105`; a cone that is not a whole single-shell, single-side, single-loop canonical cone is refused at `:296-325`.
3. Face pairs are pruned by AABB before surface dispatch. The first statement of the inner pair loop is `curved_boolean_bounds_overlap(face_a.bounds, face_b.bounds, tolerance)` at `:2117`; non-overlapping pairs are skipped. This is a cost filter only: an overlapping pair still reaches the gate.
4. The seam-copy guard admits apex copies. It computes `apex = origin + axis * (radius / tan(half_angle))` at `:2682-2683` and permits a coincident consecutive pair when both points are within the Euclidean `tolerance` of that apex (`:2688-2689`). The exemption is per pair; it does not widen the whole-cone, cutter-type, or any other guard, and coincident seam copies away from the apex remain refused.

## Evidence

Fixture: canonical cone `make_cone({0,0,0}, +z, radius 1, height 2)` and box `{0,-1.5,-0.5}..{1.5,1.5,2.5}`, boolean tolerance `1e-9`. QA-22 already covers the perpendicular contract; the six cut cases and the predicate pin live in `tests/unit/brep_boolean_curved_test.cpp:936-1065`. Volumes are measured from the returned B-rep by `detail::curved_solid_volume`, not copied from the analytic formula.

| Operation | Faces | Edges | Vertices | Shells | Zero-length edges | Measured volume | Analytic | Relative error |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| difference | 3 | 64 | 63 | 1 | 0 | 1.0471975511965941 | π/3 = 1.0471975511965976 | 3.4e-15 |
| intersection | 4 | 66 | 64 | 1 | 0 | 1.0471975511965979 | π/3 | 2.1e-16 |
| union | 8 | 76 | 71 | 1 | 0 | 14.547197551196593 | 13.5 + π/3 = 14.547197551196598 | 3.7e-16 |

All three results are one shell, pass `validate_topology` and `shell_is_closed_by_uses`, have positive volume, and reproduce bitwise on a second difference run. A cone rotated to the `+x` axis (box `{-0.5,0,-1.5}..{2.5,1.5,1.5}`) and a cone translated to `(0,5,0)` (box `{0,3.5,-0.5}..{1.5,6.5,2.5}`) both give difference `1.0471975511965941` with `3/64/63`, one closed shell, no zero-length edges, and a bitwise-identical repeat.

Refusals stay deterministic and leave no entity: a box just below the apex (`{0,-1.5,1.5}..{1.5,1.5,2.5}`) is `unsupported` from `split_face_curved`; a cutter face narrower than the cone diameter (`{0,-1.5,0.2}..{1.5,1.5,1.0}`) is `unsupported` with `"trim_face: cut chain ends inside the parameter domain"`; a cutter whose plane contains a seam ray (`{-1.5,0,-0.5}..{1.5,1.5,2.5}`) is `unsupported` with `"boolean_solid_curved: result touches itself along an edge"`. The near-apex case pins no message because it is a downstream splitter refusal, not a gate contract.

Release focused `-tc="curved boolean:*"` runs 51/51 cases and 1492/1492 assertions; Debug focused matches. Both full configurations run 1038/1038 cases and 674,758/674,758 assertions with zero failures.

## Consequences

- Axis-containing planes are now a supported whole-cone/plane cut class, with exact two-ruling sections and no synthetic chains. The positive path's volume agrees with `π/3` and `13.5 + π/3` to a few ulps, and the topology is closed.
- The gate message only reports oblique planes: `"boolean_solid_curved: cone cutters must be perpendicular to or contain the cone axis; oblique planes are deferred to M3c-3c-4b3"` (`:2132-2134`). The predicate comment documents the same split.
- Limits remain: oblique planes (ellipse sections) and offset-parallel planes, reversed difference, cone/cone, cone/torus, cone/cylinder, cone/sphere pairs, trimmed cone fragments as boolean operands, a cutter narrower than the cone diameter, seam-ray cutters, a dedicated apex-only diagnostic, and same-face double pole-run merging are not supported. The apex exemption does not make translated cones with truly coincident seam copies legal; those are still refused.
- The AABB filter changes cost, not semantics. Worst-case gate work is still one predicate per overlapping face pair.
- Any future ellipse-capable splitter should replace the gate rather than bypass it, so the supported class stays explicit and testable.