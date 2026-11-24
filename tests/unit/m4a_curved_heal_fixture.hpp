#pragma once

// Fixtures for the M4a curved-shell healing tests.
//
// Every positive fixture is ONE solid whose faces live in ONE open shell; the
// shell is the union of two components separated by a narrow gap along a pair
// of coincident analytic carriers. Faces carry no pcurves, so face_loops_uv
// projects the edge samples on demand -- the state heal_brep leaves a model in
// before its pcurve-rebuild phase rewrites them.
//
// All circles are split into two half-turn arcs at u = 0 and u = pi, matching
// what split_solid produces for a full-circle rim, so heal_plan never sees a
// closed (zero chord length) boundary edge.
//
// The negative fixtures are used by later tasks: a radius mismatch, an
// off-surface band, a closed-rim representation and a plane lid. The NURBS
// fixtures are the ones this task's kind guards must reject.

#include <array>
#include <cmath>
#include <cstdint>
#include <numbers>
#include <span>
#include <utility>

#include <geom/brep/solid_editor.hpp>
#include <geom/brep/topology.hpp>
#include <geom/curves/analytic_curve.hpp>
#include <geom/curves/nurbs_curve.hpp>
#include <geom/surfaces/analytic.hpp>
#include <geom/surfaces/nurbs_surface.hpp>
#include <geom/surfaces/plane.hpp>

namespace geom::test::m4a {

constexpr double kPi = std::numbers::pi_v<double>;
constexpr double kTwoPi = 2.0 * kPi;

inline Point3d m4a_point(double x, double y, double z) { return {x, y, z}; }

// Adds a vertex and returns its handle; aborts on failure (fixtures are
// test-only, so a build failure is a test bug).
inline VertexId m4a_vertex(SolidEditor& editor, Point3d p) {
  return editor.add_vertex(p, 0.0).value();
}

inline EdgeId m4a_edge(SolidEditor& editor, VertexId a, VertexId b,
                       std::uint32_t curve, double t0, double t1) {
  return editor.add_edge(a, b, curve, t0, t1, 0.0).value();
}

inline FaceId m4a_face(SolidEditor& editor, std::uint32_t surface) {
  return editor.add_face(surface, false, 0.0).value();
}

// Same, with the face's surface normal flipped relative to the loop winding.
// Needed when the chart orientation (not the loop) puts the outer side on the
// other side, as on the inner branch of the torus.
inline FaceId m4a_face_reversed(SolidEditor& editor, std::uint32_t surface) {
  return editor.add_face(surface, true, 0.0).value();
}

inline SolidId m4a_add_solid(SolidEditor& editor, std::span<const FaceId> faces) {
  const ShellId shell = editor.add_shell(faces, false).value();
  return editor.add_solid(std::array<ShellId, 1>{shell}, 0.0).value();
}

struct CurvedSolid {
  Brep brep;
  SolidId solid{};
};

// ---------------------------------------------------------------------------
// Case 1: cylinder radius r, height 2, radial wall split at z = 1, upper band
// moved +gap. Both lateral faces use the SAME Cylinder3 surface (origin at the
// base) unless `upper_radius` differs; then the upper band gets its own
// cylinder, so the two carriers disagree while the rim vertices stay within
// max_gap.
// ---------------------------------------------------------------------------
// `lower_top_circle_z` can lift the kept lower rim's curve without moving its
// vertices, which models a corrupt kept edge for the kept-curve guard.
inline CurvedSolid make_cylinder_bands(double radius, double gap,
                                       double upper_radius, bool off_surface,
                                       bool shifted_chart,
                                       double lower_top_circle_z = 1.0) {
  SolidEditor editor;
  const std::uint32_t cylinder = editor.add_surface(Cylinder3{
      m4a_point(0.0, 0.0, 0.0), {0.0, 0.0, 1.0}, {1.0, 0.0, 0.0}, radius});
  // A different-radius upper band is a genuinely different carrier. When
  // `off_surface` is set the upper band keeps the base surface while its rims
  // sit at `upper_radius`, which models corrupt input rather than a second
  // carrier. `shifted_chart` gives the upper band an equal carrier whose
  // origin (and so its v chart) sits at z = 1, as a split fixture would
  // produce.
  const std::uint32_t upper_cylinder =
      (off_surface || (upper_radius == radius && !shifted_chart))
          ? cylinder
          : editor.add_surface(Cylinder3{
                m4a_point(0.0, 0.0, shifted_chart ? 1.0 : 0.0),
                {0.0, 0.0, 1.0}, {1.0, 0.0, 0.0}, upper_radius});
  const std::uint32_t bottom_plane =
      editor.add_surface(Plane3{{0.0, 0.0, -1.0}, 0.0});
  const std::uint32_t top_plane =
      editor.add_surface(Plane3{{0.0, 0.0, 1.0}, 2.0 + gap});

  const std::uint32_t circle_bottom = editor.add_curve(
      CircleCurve3{m4a_point(0.0, 0.0, 0.0), {1.0, 0.0, 0.0}, {0.0, 1.0, 0.0}, radius});
  const std::uint32_t circle_lower_top = editor.add_curve(CircleCurve3{
      m4a_point(0.0, 0.0, lower_top_circle_z), {1.0, 0.0, 0.0}, {0.0, 1.0, 0.0},
      radius});
  const std::uint32_t circle_upper_bottom = editor.add_curve(CircleCurve3{
      m4a_point(0.0, 0.0, 1.0 + gap), {1.0, 0.0, 0.0}, {0.0, 1.0, 0.0}, upper_radius});
  const std::uint32_t circle_upper_top = editor.add_curve(CircleCurve3{
      m4a_point(0.0, 0.0, 2.0 + gap), {1.0, 0.0, 0.0}, {0.0, 1.0, 0.0}, upper_radius});
  const std::uint32_t seam_lower = editor.add_curve(
      LineCurve3{m4a_point(radius, 0.0, 0.0), {0.0, 0.0, 1.0}});
  const std::uint32_t seam_upper = editor.add_curve(
      LineCurve3{m4a_point(upper_radius, 0.0, 1.0 + gap), {0.0, 0.0, 1.0}});

  const VertexId b0 = m4a_vertex(editor, m4a_point(radius, 0.0, 0.0));
  const VertexId b1 = m4a_vertex(editor, m4a_point(-radius, 0.0, 0.0));
  const VertexId t0 = m4a_vertex(editor, m4a_point(radius, 0.0, 1.0));
  const VertexId t1 = m4a_vertex(editor, m4a_point(-radius, 0.0, 1.0));
  const VertexId u0 = m4a_vertex(editor, m4a_point(upper_radius, 0.0, 1.0 + gap));
  const VertexId u1 = m4a_vertex(editor, m4a_point(-upper_radius, 0.0, 1.0 + gap));
  const VertexId v0 = m4a_vertex(editor, m4a_point(upper_radius, 0.0, 2.0 + gap));
  const VertexId v1 = m4a_vertex(editor, m4a_point(-upper_radius, 0.0, 2.0 + gap));

  const EdgeId bottom0 = m4a_edge(editor, b0, b1, circle_bottom, 0.0, kPi);
  const EdgeId bottom1 = m4a_edge(editor, b1, b0, circle_bottom, kPi, kTwoPi);
  const EdgeId lower_top0 = m4a_edge(editor, t0, t1, circle_lower_top, 0.0, kPi);
  const EdgeId lower_top1 = m4a_edge(editor, t1, t0, circle_lower_top, kPi, kTwoPi);
  const EdgeId upper_bottom0 =
      m4a_edge(editor, u0, u1, circle_upper_bottom, 0.0, kPi);
  const EdgeId upper_bottom1 =
      m4a_edge(editor, u1, u0, circle_upper_bottom, kPi, kTwoPi);
  const EdgeId upper_top0 = m4a_edge(editor, v0, v1, circle_upper_top, 0.0, kPi);
  const EdgeId upper_top1 = m4a_edge(editor, v1, v0, circle_upper_top, kPi, kTwoPi);
  const EdgeId lower_seam = m4a_edge(editor, b0, t0, seam_lower, 0.0, 1.0);
  const EdgeId upper_seam = m4a_edge(editor, u0, v0, seam_upper, 0.0, 1.0);

  const FaceId lower_lateral = m4a_face(editor, cylinder);
  const std::array<CoedgeSpec, 6> lower_loop{{
      {bottom0, false}, {bottom1, false}, {lower_seam, false},
      {lower_top1, true}, {lower_top0, true}, {lower_seam, true}}};
  editor.add_loop(lower_lateral, lower_loop, true).value();

  const FaceId bottom_cap = m4a_face(editor, bottom_plane);
  const std::array<CoedgeSpec, 2> bottom_loop{{
      {bottom0, true}, {bottom1, true}}};
  editor.add_loop(bottom_cap, bottom_loop, true).value();

  const FaceId upper_lateral = m4a_face(editor, upper_cylinder);
  const std::array<CoedgeSpec, 6> upper_loop{{
      {upper_bottom0, false}, {upper_bottom1, false}, {upper_seam, false},
      {upper_top1, true}, {upper_top0, true}, {upper_seam, true}}};
  editor.add_loop(upper_lateral, upper_loop, true).value();

  const FaceId top_cap = m4a_face(editor, top_plane);
  const std::array<CoedgeSpec, 2> top_loop{{
      {upper_top0, false}, {upper_top1, false}}};
  editor.add_loop(top_cap, top_loop, true).value();

  const std::array<FaceId, 4> faces{
      {lower_lateral, bottom_cap, upper_lateral, top_cap}};
  const SolidId solid = m4a_add_solid(editor, faces);
  return CurvedSolid{editor.commit().value(), solid};
}

inline CurvedSolid make_two_band_cylinder(double radius, double gap) {
  return make_cylinder_bands(radius, gap, radius, false, false);
}

// Case 1-negative-b: same construction but the upper band sits on a cylinder
// with a different radius; `upper_radius` is 1.001 for the distance-rejected
// case and 1 + 2e-7 for the carrier-rejected case.
inline CurvedSolid make_cylinder_radius_gap(double gap, double upper_radius) {
  return make_cylinder_bands(1.0, gap, upper_radius, false, false);
}

// Case 4e: both carriers agree and every rim curve lies on the shared cylinder,
// but the kept lower rim curve is lifted to z = 1.5 while its vertices stay at
// z = 1. Only the kept-curve-versus-welded-vertices guard can catch this.
inline CurvedSolid make_cylinder_kept_curve_gap(double gap) {
  return make_cylinder_bands(1.0, gap, 1.0, false, false, 1.5);
}

// Case 4d: both bands claim the same base cylinder but the upper band's rims
// sit at 1 + 2e-7, i.e. its boundary already floats off its own face. This is
// corrupt input: the surface predicate alone accepts it, so the
// edge-on-surface guard must reject it before welding.
inline CurvedSolid make_cylinder_off_surface_gap(
    double gap, double upper_radius = 1.0 + 2e-7) {
  return make_cylinder_bands(1.0, gap, upper_radius, true, false);
}

// Case 1e: the same split, but each band's rims are single CLOSED circle edges
// (v0 == v1) instead of two half-turn arcs. This is the representation
// make_cylinder produces; the gaps are real circles, so the healer must say
// why it cannot pair them.
inline CurvedSolid make_cylinder_closed_rims_gap(double gap) {
  SolidEditor editor;
  const std::uint32_t cylinder = editor.add_surface(Cylinder3{
      m4a_point(0.0, 0.0, 0.0), {0.0, 0.0, 1.0}, {1.0, 0.0, 0.0}, 1.0});
  const std::uint32_t bottom_plane =
      editor.add_surface(Plane3{{0.0, 0.0, -1.0}, 0.0});
  const std::uint32_t top_plane =
      editor.add_surface(Plane3{{0.0, 0.0, 1.0}, 2.0 + gap});
  const std::uint32_t circle_bottom = editor.add_curve(
      CircleCurve3{m4a_point(0.0, 0.0, 0.0), {1.0, 0.0, 0.0}, {0.0, 1.0, 0.0}, 1.0});
  const std::uint32_t circle_lower_top = editor.add_curve(
      CircleCurve3{m4a_point(0.0, 0.0, 1.0), {1.0, 0.0, 0.0}, {0.0, 1.0, 0.0}, 1.0});
  const std::uint32_t circle_upper_bottom = editor.add_curve(CircleCurve3{
      m4a_point(0.0, 0.0, 1.0 + gap), {1.0, 0.0, 0.0}, {0.0, 1.0, 0.0}, 1.0});
  const std::uint32_t circle_upper_top = editor.add_curve(CircleCurve3{
      m4a_point(0.0, 0.0, 2.0 + gap), {1.0, 0.0, 0.0}, {0.0, 1.0, 0.0}, 1.0});
  const std::uint32_t seam_lower = editor.add_curve(
      LineCurve3{m4a_point(1.0, 0.0, 0.0), {0.0, 0.0, 1.0}});
  const std::uint32_t seam_upper = editor.add_curve(
      LineCurve3{m4a_point(1.0, 0.0, 1.0 + gap), {0.0, 0.0, 1.0}});

  const VertexId b0 = m4a_vertex(editor, m4a_point(1.0, 0.0, 0.0));
  const VertexId t0 = m4a_vertex(editor, m4a_point(1.0, 0.0, 1.0));
  const VertexId u0 = m4a_vertex(editor, m4a_point(1.0, 0.0, 1.0 + gap));
  const VertexId v0 = m4a_vertex(editor, m4a_point(1.0, 0.0, 2.0 + gap));

  const EdgeId bottom = m4a_edge(editor, b0, b0, circle_bottom, 0.0, kTwoPi);
  const EdgeId lower_top = m4a_edge(editor, t0, t0, circle_lower_top, 0.0, kTwoPi);
  const EdgeId upper_bottom =
      m4a_edge(editor, u0, u0, circle_upper_bottom, 0.0, kTwoPi);
  const EdgeId upper_top =
      m4a_edge(editor, v0, v0, circle_upper_top, 0.0, kTwoPi);
  const EdgeId lower_seam = m4a_edge(editor, b0, t0, seam_lower, 0.0, 1.0);
  const EdgeId upper_seam = m4a_edge(editor, u0, v0, seam_upper, 0.0, 1.0);

  const FaceId lower_lateral = m4a_face(editor, cylinder);
  const std::array<CoedgeSpec, 4> lower_loop{{
      {bottom, false}, {lower_seam, false},
      {lower_top, true}, {lower_seam, true}}};
  editor.add_loop(lower_lateral, lower_loop, true).value();

  const FaceId bottom_cap = m4a_face(editor, bottom_plane);
  const std::array<CoedgeSpec, 1> bottom_loop{{{bottom, true}}};
  editor.add_loop(bottom_cap, bottom_loop, true).value();

  const FaceId upper_lateral = m4a_face(editor, cylinder);
  const std::array<CoedgeSpec, 4> upper_loop{{
      {upper_bottom, false}, {upper_seam, false},
      {upper_top, true}, {upper_seam, true}}};
  editor.add_loop(upper_lateral, upper_loop, true).value();

  const FaceId top_cap = m4a_face(editor, top_plane);
  const std::array<CoedgeSpec, 1> top_loop{{{upper_top, false}}};
  editor.add_loop(top_cap, top_loop, true).value();

  const std::array<FaceId, 4> faces{
      {lower_lateral, bottom_cap, upper_lateral, top_cap}};
  const SolidId solid = m4a_add_solid(editor, faces);
  return CurvedSolid{editor.commit().value(), solid};
}

// Case 4c: lower cylinder band closed at z = 0, upper lid is a PLANE disk at
// z = 1 + gap whose boundary circle has the same radius. The carriers are a
// cylinder and a plane; the healer must reject instead of welding.
inline CurvedSolid make_cylinder_plane_gap(double gap) {
  SolidEditor editor;
  const std::uint32_t cylinder = editor.add_surface(Cylinder3{
      m4a_point(0.0, 0.0, 0.0), {0.0, 0.0, 1.0}, {1.0, 0.0, 0.0}, 1.0});
  const std::uint32_t bottom_plane =
      editor.add_surface(Plane3{{0.0, 0.0, -1.0}, 0.0});
  const std::uint32_t lid_plane =
      editor.add_surface(Plane3{{0.0, 0.0, 1.0}, 1.0 + gap});

  const std::uint32_t circle_bottom = editor.add_curve(
      CircleCurve3{m4a_point(0.0, 0.0, 0.0), {1.0, 0.0, 0.0}, {0.0, 1.0, 0.0}, 1.0});
  const std::uint32_t circle_top = editor.add_curve(
      CircleCurve3{m4a_point(0.0, 0.0, 1.0), {1.0, 0.0, 0.0}, {0.0, 1.0, 0.0}, 1.0});
  const std::uint32_t circle_lid = editor.add_curve(CircleCurve3{
      m4a_point(0.0, 0.0, 1.0 + gap), {1.0, 0.0, 0.0}, {0.0, 1.0, 0.0}, 1.0});
  const std::uint32_t seam = editor.add_curve(
      LineCurve3{m4a_point(1.0, 0.0, 0.0), {0.0, 0.0, 1.0}});

  const VertexId b0 = m4a_vertex(editor, m4a_point(1.0, 0.0, 0.0));
  const VertexId b1 = m4a_vertex(editor, m4a_point(-1.0, 0.0, 0.0));
  const VertexId t0 = m4a_vertex(editor, m4a_point(1.0, 0.0, 1.0));
  const VertexId t1 = m4a_vertex(editor, m4a_point(-1.0, 0.0, 1.0));
  const VertexId l0 = m4a_vertex(editor, m4a_point(1.0, 0.0, 1.0 + gap));
  const VertexId l1 = m4a_vertex(editor, m4a_point(-1.0, 0.0, 1.0 + gap));

  const EdgeId bottom0 = m4a_edge(editor, b0, b1, circle_bottom, 0.0, kPi);
  const EdgeId bottom1 = m4a_edge(editor, b1, b0, circle_bottom, kPi, kTwoPi);
  const EdgeId top0 = m4a_edge(editor, t0, t1, circle_top, 0.0, kPi);
  const EdgeId top1 = m4a_edge(editor, t1, t0, circle_top, kPi, kTwoPi);
  const EdgeId lid0 = m4a_edge(editor, l0, l1, circle_lid, 0.0, kPi);
  const EdgeId lid1 = m4a_edge(editor, l1, l0, circle_lid, kPi, kTwoPi);
  const EdgeId seam_edge = m4a_edge(editor, b0, t0, seam, 0.0, 1.0);

  const FaceId lateral = m4a_face(editor, cylinder);
  const std::array<CoedgeSpec, 6> lateral_loop{{
      {bottom0, false}, {bottom1, false}, {seam_edge, false},
      {top1, true}, {top0, true}, {seam_edge, true}}};
  editor.add_loop(lateral, lateral_loop, true).value();

  const FaceId cap = m4a_face(editor, bottom_plane);
  const std::array<CoedgeSpec, 2> cap_loop{{
      {bottom0, true}, {bottom1, true}}};
  editor.add_loop(cap, cap_loop, true).value();

  const FaceId lid = m4a_face(editor, lid_plane);
  const std::array<CoedgeSpec, 2> lid_loop{{{lid0, false}, {lid1, false}}};
  editor.add_loop(lid, lid_loop, true).value();

  const std::array<FaceId, 3> faces{{lateral, cap, lid}};
  const SolidId solid = m4a_add_solid(editor, faces);
  return CurvedSolid{editor.commit().value(), solid};
}

// ---------------------------------------------------------------------------
// Case 2: sphere radius 1 split at the equator, boundary circles at
// v = +/-eps on the SAME sphere. Each hemisphere is two lune faces (u in
// [0, pi] and [pi, 2pi]) so both uv polygons are proper rectangles and the
// pole fans close.
// ---------------------------------------------------------------------------
inline CurvedSolid make_sphere_gap(double eps) {
  SolidEditor editor;
  const double radius = std::sqrt(1.0 - eps * eps);
  const std::uint32_t sphere =
      editor.add_surface(Sphere3{m4a_point(0.0, 0.0, 0.0), {0.0, 0.0, 1.0},
                                 {1.0, 0.0, 0.0}, 1.0});

  const std::uint32_t circle_north = editor.add_curve(CircleCurve3{
      m4a_point(0.0, 0.0, eps), {1.0, 0.0, 0.0}, {0.0, 1.0, 0.0}, radius});
  const std::uint32_t circle_south = editor.add_curve(CircleCurve3{
      m4a_point(0.0, 0.0, -eps), {1.0, 0.0, 0.0}, {0.0, 1.0, 0.0}, radius});
  const std::uint32_t meridian_n0 = editor.add_curve(CircleCurve3{
      m4a_point(0.0, 0.0, 0.0), {1.0, 0.0, 0.0}, {0.0, 0.0, 1.0}, 1.0});
  const std::uint32_t meridian_npi = editor.add_curve(CircleCurve3{
      m4a_point(0.0, 0.0, 0.0), {-1.0, 0.0, 0.0}, {0.0, 0.0, 1.0}, 1.0});
  const std::uint32_t meridian_s0 = editor.add_curve(CircleCurve3{
      m4a_point(0.0, 0.0, 0.0), {1.0, 0.0, 0.0}, {0.0, 0.0, -1.0}, 1.0});
  const std::uint32_t meridian_spi = editor.add_curve(CircleCurve3{
      m4a_point(0.0, 0.0, 0.0), {-1.0, 0.0, 0.0}, {0.0, 0.0, -1.0}, 1.0});

  const VertexId n0 = m4a_vertex(editor, m4a_point(radius, 0.0, eps));
  const VertexId n1 = m4a_vertex(editor, m4a_point(-radius, 0.0, eps));
  const VertexId np = m4a_vertex(editor, m4a_point(0.0, 0.0, 1.0));
  const VertexId s0 = m4a_vertex(editor, m4a_point(radius, 0.0, -eps));
  const VertexId s1 = m4a_vertex(editor, m4a_point(-radius, 0.0, -eps));
  const VertexId sp = m4a_vertex(editor, m4a_point(0.0, 0.0, -1.0));

  const double theta = std::asin(eps);  // latitude of the equator circle
  const EdgeId north_east = m4a_edge(editor, n0, n1, circle_north, 0.0, kPi);
  const EdgeId north_west = m4a_edge(editor, n1, n0, circle_north, kPi, kTwoPi);
  const EdgeId south_east = m4a_edge(editor, s0, s1, circle_south, 0.0, kPi);
  const EdgeId south_west = m4a_edge(editor, s1, s0, circle_south, kPi, kTwoPi);
  const EdgeId nmer0 = m4a_edge(editor, n0, np, meridian_n0, theta, kPi / 2.0);
  const EdgeId nmerpi = m4a_edge(editor, n1, np, meridian_npi, theta, kPi / 2.0);
  const EdgeId smer0 = m4a_edge(editor, s0, sp, meridian_s0, theta, kPi / 2.0);
  const EdgeId smerpi = m4a_edge(editor, s1, sp, meridian_spi, theta, kPi / 2.0);

  const FaceId north_east_face = m4a_face(editor, sphere);
  const std::array<CoedgeSpec, 3> north_east_loop{{
      {north_east, false}, {nmerpi, false}, {nmer0, true}}};
  editor.add_loop(north_east_face, north_east_loop, true).value();

  const FaceId north_west_face = m4a_face(editor, sphere);
  const std::array<CoedgeSpec, 3> north_west_loop{{
      {north_west, false}, {nmer0, false}, {nmerpi, true}}};
  editor.add_loop(north_west_face, north_west_loop, true).value();

  const FaceId south_east_face = m4a_face(editor, sphere);
  const std::array<CoedgeSpec, 3> south_east_loop{{
      {south_east, true}, {smer0, false}, {smerpi, true}}};
  editor.add_loop(south_east_face, south_east_loop, true).value();

  const FaceId south_west_face = m4a_face(editor, sphere);
  const std::array<CoedgeSpec, 3> south_west_loop{{
      {south_west, true}, {smerpi, false}, {smer0, true}}};
  editor.add_loop(south_west_face, south_west_loop, true).value();

  const std::array<FaceId, 4> faces{
      {north_east_face, north_west_face, south_east_face, south_west_face}};
  const SolidId solid = m4a_add_solid(editor, faces);
  return CurvedSolid{editor.commit().value(), solid};
}

// ---------------------------------------------------------------------------
// Case 3: cone radius 1, height 2 (half angle atan(0.5), apex at v = 2) with
// the band v in [1, 1+gap] removed. The upper band ends at the apex, so its
// face loop is [two base arcs, seam to apex, seam back].
// ---------------------------------------------------------------------------
inline CurvedSolid make_cone_gap(double gap) {
  SolidEditor editor;
  const double tan_half = 0.5;
  const double half_angle = std::atan(tan_half);
  const std::uint32_t cone = editor.add_surface(Cone3{
      m4a_point(0.0, 0.0, 0.0), {0.0, 0.0, 1.0}, {1.0, 0.0, 0.0}, 1.0, half_angle});
  const std::uint32_t bottom_plane =
      editor.add_surface(Plane3{{0.0, 0.0, -1.0}, 0.0});

  const double r_lower = 1.0 - tan_half * 1.0;          // 0.5 at v = 1
  const double r_upper = 1.0 - tan_half * (1.0 + gap);  // 0.5 - 0.5*gap
  const std::uint32_t circle_bottom = editor.add_curve(
      CircleCurve3{m4a_point(0.0, 0.0, 0.0), {1.0, 0.0, 0.0}, {0.0, 1.0, 0.0}, 1.0});
  const std::uint32_t circle_lower_top = editor.add_curve(CircleCurve3{
      m4a_point(0.0, 0.0, 1.0), {1.0, 0.0, 0.0}, {0.0, 1.0, 0.0}, r_lower});
  const std::uint32_t circle_upper_bottom = editor.add_curve(CircleCurve3{
      m4a_point(0.0, 0.0, 1.0 + gap), {1.0, 0.0, 0.0}, {0.0, 1.0, 0.0}, r_upper});
  // Lower seam follows the cone slant from (1, 0, 0) to (r_lower, 0, 1).
  const Vec3d lower_delta{r_lower - 1.0, 0.0, 1.0};
  const double lower_len =
      std::hypot(lower_delta[0], lower_delta[1], lower_delta[2]);
  const Vec3d lower_dir{lower_delta[0] / lower_len, lower_delta[1] / lower_len,
                        lower_delta[2] / lower_len};
  const std::uint32_t seam_lower =
      editor.add_curve(LineCurve3{m4a_point(1.0, 0.0, 0.0), lower_dir});
  // Upper seam: from (r_upper, 0, 1+gap) to the apex (0, 0, 2).
  const Vec3d upper_delta{-r_upper, 0.0, 1.0 - gap};
  const double upper_len =
      std::hypot(upper_delta[0], upper_delta[1], upper_delta[2]);
  const Vec3d upper_dir{upper_delta[0] / upper_len, upper_delta[1] / upper_len,
                        upper_delta[2] / upper_len};
  const std::uint32_t seam_upper = editor.add_curve(
      LineCurve3{m4a_point(r_upper, 0.0, 1.0 + gap), upper_dir});

  const VertexId b0 = m4a_vertex(editor, m4a_point(1.0, 0.0, 0.0));
  const VertexId b1 = m4a_vertex(editor, m4a_point(-1.0, 0.0, 0.0));
  const VertexId t0 = m4a_vertex(editor, m4a_point(r_lower, 0.0, 1.0));
  const VertexId t1 = m4a_vertex(editor, m4a_point(-r_lower, 0.0, 1.0));
  const VertexId u0 = m4a_vertex(editor, m4a_point(r_upper, 0.0, 1.0 + gap));
  const VertexId u1 = m4a_vertex(editor, m4a_point(-r_upper, 0.0, 1.0 + gap));
  const VertexId apex = m4a_vertex(editor, m4a_point(0.0, 0.0, 2.0));

  const EdgeId bottom0 = m4a_edge(editor, b0, b1, circle_bottom, 0.0, kPi);
  const EdgeId bottom1 = m4a_edge(editor, b1, b0, circle_bottom, kPi, kTwoPi);
  const EdgeId lower_top0 = m4a_edge(editor, t0, t1, circle_lower_top, 0.0, kPi);
  const EdgeId lower_top1 = m4a_edge(editor, t1, t0, circle_lower_top, kPi, kTwoPi);
  const EdgeId upper_bottom0 =
      m4a_edge(editor, u0, u1, circle_upper_bottom, 0.0, kPi);
  const EdgeId upper_bottom1 =
      m4a_edge(editor, u1, u0, circle_upper_bottom, kPi, kTwoPi);
  const EdgeId lower_seam = m4a_edge(editor, b0, t0, seam_lower, 0.0, lower_len);
  const EdgeId upper_seam =
      m4a_edge(editor, u0, apex, seam_upper, 0.0, upper_len);

  const FaceId lower_lateral = m4a_face(editor, cone);
  const std::array<CoedgeSpec, 6> lower_loop{{
      {bottom0, false}, {bottom1, false}, {lower_seam, false},
      {lower_top1, true}, {lower_top0, true}, {lower_seam, true}}};
  editor.add_loop(lower_lateral, lower_loop, true).value();

  const FaceId cap = m4a_face(editor, bottom_plane);
  const std::array<CoedgeSpec, 2> cap_loop{{
      {bottom0, true}, {bottom1, true}}};
  editor.add_loop(cap, cap_loop, true).value();

  const FaceId upper_lateral = m4a_face(editor, cone);
  const std::array<CoedgeSpec, 4> upper_loop{{
      {upper_bottom0, false}, {upper_bottom1, false},
      {upper_seam, false}, {upper_seam, true}}};
  editor.add_loop(upper_lateral, upper_loop, true).value();

  const std::array<FaceId, 3> faces{{lower_lateral, cap, upper_lateral}};
  const SolidId solid = m4a_add_solid(editor, faces);
  return CurvedSolid{editor.commit().value(), solid};
}

// ---------------------------------------------------------------------------
// M4b-1: torus band fixtures. Each is ONE solid in ONE open shell: two bands of
// the SAME torus (R = 1, r = 0.5) separated by a narrow gap along a pair of open
// rim circles, with one planar cap per band. No rim is a closed edge, so the
// T1/T2 pre-split pass has nothing to split here: the planner pairs the open
// halves directly.
//
// Every rim is split into two half-turn arcs at parameter 0 and pi. The healer
// is expected to weld each pair and keep the union of the two bands, i.e. the
// gap becomes material: the healed volumes are the volumes of the un-gapped
// shapes, computed analytically (Pappus for the sector, a meridian integral for
// the barrel bands) and written as the brief's quadrature constants.
//
// The open rim vertex chord (the reported gap) is chosen by construction:
// 2*(R + r)*sin(delta/2) = gap for a u gap and 2*r*sin(delta/2) = gap for a
// v gap. The other endpoint of each pair is closer (the inner equator chord is
// (R - r)/(R + r) times the outer one), so `gap` is the pair's maximum.
// ---------------------------------------------------------------------------
namespace torus_fixture {

inline constexpr double kMajorRadius = 1.0;
inline constexpr double kMinorRadius = 0.5;

// Rim at constant u, parameterized by the chart v. The canonical frame
// (center at the origin, axis +z, x_axis +x) makes radial(u) = (cos u, sin u,
// 0), so the rim is a circle of radius r around R*radial(u) in the plane that
// radial(u) and the axis span. Requires v_origin == 0.
inline CircleCurve3 u_rim_circle(const Torus3& torus, double u) {
  const Vec3d radial{std::cos(u), std::sin(u), 0.0};
  return CircleCurve3{torus.center + radial * kMajorRadius, radial,
                      Vec3d{0.0, 0.0, 1.0}, kMinorRadius};
}

// Rim at constant chart v, parameterized by u. Torus3::eval maps chart v to the
// physical vv = v + v_origin, so the rim is the circle in the plane
// z = r*sin(vv) with radius R + r*cos(vv).
inline CircleCurve3 v_rim_circle(const Torus3& torus, double v) {
  const double vv = torus.v_origin == 0.0 ? v : v + torus.v_origin;
  return CircleCurve3{torus.center + Vec3d{0.0, 0.0, kMinorRadius * std::sin(vv)},
                      Vec3d{1.0, 0.0, 0.0}, Vec3d{0.0, 1.0, 0.0},
                      kMajorRadius + kMinorRadius * std::cos(vv)};
}

// The u = 0 meridian, parameterized by the chart v: its frame is rotated by
// v_origin so that eval(t) == torus.eval(0, t) for every chart v_origin.
inline CircleCurve3 seam_meridian(const Torus3& torus) {
  const double v0 = torus.v_origin == 0.0 ? 0.0 : torus.v_origin;
  const Vec3d radial{1.0, 0.0, 0.0};
  const Vec3d axis{0.0, 0.0, 1.0};
  const Vec3d x = radial * std::cos(v0) + axis * std::sin(v0);
  const Vec3d y = radial * (-std::sin(v0)) + axis * std::cos(v0);
  return CircleCurve3{torus.center + radial * kMajorRadius, x, y, kMinorRadius};
}

struct Rim {
  VertexId p0{};
  VertexId p1{};
  EdgeId arc0{};
  EdgeId arc1{};
};

// Adds the rim carrier plus its two half-turn arcs (0..pi, pi..2pi). Vertices
// are taken from the surface so they lie exactly on it; the carrier differs
// from them by at most a few ulps, far below the healer's 1e-9 check
// tolerance.
inline Rim add_rim(SolidEditor& editor, const CircleCurve3& circle,
                   Point3d p0, Point3d p1) {
  const std::uint32_t curve = editor.add_curve(circle);
  const VertexId v0 = m4a_vertex(editor, p0);
  const VertexId v1 = m4a_vertex(editor, p1);
  return Rim{v0, v1, m4a_edge(editor, v0, v1, curve, 0.0, kPi),
             m4a_edge(editor, v1, v0, curve, kPi, kTwoPi)};
}

}  // namespace torus_fixture

// Case 5a: a 3.0 rad sector around the axis, split in the middle. Both bands
// span every v, so the open rims are u-circles. Healed volume by Pappus:
//   V = theta * R * pi * r^2 = 3.0 * 1 * pi * 0.25 = 3*pi/4
//     = 2.356194490192345...
// (the brief's constant is the same value from quadrature).
inline CurvedSolid make_torus_u_gap(double gap = 1e-6) {
  using namespace torus_fixture;
  const Torus3 torus{m4a_point(0.0, 0.0, 0.0), {0.0, 0.0, 1.0}, {1.0, 0.0, 0.0},
                     kMajorRadius, kMinorRadius, 0.0};
  const double u_cap_low = 0.0;
  const double u_cap_high = 3.0;
  // The open rims sit delta apart, where the outer-equator chord
  // 2*(R + r)*sin(delta/2) is `gap` minus a 1e-8 relative margin. The margin
  // absorbs the rounding of 1.5 +/- delta/2 (about 1e-10 relative here) so the
  // measured vertex chord stays below the default 1e-6 max_gap; the inner
  // equator chord is (R - r)/(R + r) times smaller, so it is the outer chord
  // that the planner reports.
  const double target_chord = gap * (1.0 - 1e-8);
  const double delta = 2.0 * std::asin(target_chord /
                                       (2.0 * (kMajorRadius + kMinorRadius)));
  const double u_open_low = 0.5 * (u_cap_low + u_cap_high) - 0.5 * delta;
  const double u_open_high = 0.5 * (u_cap_low + u_cap_high) + 0.5 * delta;

  SolidEditor editor;
  const std::uint32_t surface = editor.add_surface(torus);
  // The cap at u = 0 faces -y, the one at u = 3 faces +tangent(3).
  const std::uint32_t cap_low_plane =
      editor.add_surface(Plane3{{0.0, -1.0, 0.0}, 0.0});
  const std::uint32_t cap_high_plane = editor.add_surface(
      Plane3{{-std::sin(u_cap_high), std::cos(u_cap_high), 0.0}, 0.0});
  // The v = 0 seam is an arc of the outer equator, parameterized by u.
  const std::uint32_t seam = editor.add_curve(
      CircleCurve3{torus.center, {1.0, 0.0, 0.0}, {0.0, 1.0, 0.0},
                   kMajorRadius + kMinorRadius});

  const Rim rim_cap_low = add_rim(editor, u_rim_circle(torus, u_cap_low),
                                  torus.eval(u_cap_low, 0.0),
                                  torus.eval(u_cap_low, kPi));
  const Rim rim_open_low = add_rim(editor, u_rim_circle(torus, u_open_low),
                                   torus.eval(u_open_low, 0.0),
                                   torus.eval(u_open_low, kPi));
  const Rim rim_open_high = add_rim(editor, u_rim_circle(torus, u_open_high),
                                    torus.eval(u_open_high, 0.0),
                                    torus.eval(u_open_high, kPi));
  const Rim rim_cap_high = add_rim(editor, u_rim_circle(torus, u_cap_high),
                                   torus.eval(u_cap_high, 0.0),
                                   torus.eval(u_cap_high, kPi));
  const EdgeId seam_low =
      m4a_edge(editor, rim_cap_low.p0, rim_open_low.p0, seam, u_cap_low,
               u_open_low);
  const EdgeId seam_high =
      m4a_edge(editor, rim_open_high.p0, rim_cap_high.p0, seam, u_open_high,
               u_cap_high);

  const FaceId lower_lateral = m4a_face(editor, surface);
  const std::array<CoedgeSpec, 6> lower_loop{{
      {rim_open_low.arc0, false}, {rim_open_low.arc1, false},
      {seam_low, true}, {rim_cap_low.arc1, true}, {rim_cap_low.arc0, true},
      {seam_low, false}}};
  editor.add_loop(lower_lateral, lower_loop, true).value();

  const FaceId cap_low = m4a_face(editor, cap_low_plane);
  const std::array<CoedgeSpec, 2> cap_low_loop{{
      {rim_cap_low.arc0, false}, {rim_cap_low.arc1, false}}};
  editor.add_loop(cap_low, cap_low_loop, true).value();

  const FaceId upper_lateral = m4a_face(editor, surface);
  const std::array<CoedgeSpec, 6> upper_loop{{
      {rim_cap_high.arc0, false}, {rim_cap_high.arc1, false},
      {seam_high, true}, {rim_open_high.arc1, true},
      {rim_open_high.arc0, true}, {seam_high, false}}};
  editor.add_loop(upper_lateral, upper_loop, true).value();

  const FaceId cap_high = m4a_face(editor, cap_high_plane);
  const std::array<CoedgeSpec, 2> cap_high_loop{{
      {rim_cap_high.arc0, true}, {rim_cap_high.arc1, true}}};
  editor.add_loop(cap_high, cap_high_loop, true).value();

  const std::array<FaceId, 4> faces{
      {lower_lateral, cap_low, upper_lateral, cap_high}};
  const SolidId solid = m4a_add_solid(editor, faces);
  return CurvedSolid{editor.commit().value(), solid};
}

// Case 5b: the OUTER branch barrel, v in [-1, 1] split at v = 0, full u
// revolution, capped by two disks at z = +/- c with c = r*sin(1). For each z
// the barrel radius is rho(z) = R + sqrt(r^2 - z^2), so the healed volume is
//   V = pi * integral_{-c}^{c} rho(z)^2 dz
//     = pi * [2*R^2*c + 2*R*(c*s + r^2*asin(c/r)) + 2*r^2*c - (2/3)*c^3]
// with s = sqrt(r^2 - c^2) and asin(c/r) = 1, i.e. 5.433419527331613...
inline CurvedSolid make_torus_outer_v_gap(double gap = 1e-6) {
  using namespace torus_fixture;
  const Torus3 torus{m4a_point(0.0, 0.0, 0.0), {0.0, 0.0, 1.0}, {1.0, 0.0, 0.0},
                     kMajorRadius, kMinorRadius, 0.0};
  const double v_cap_low = -1.0;
  const double v_cap_high = 1.0;
  // Same 1e-8 relative margin as the u-gap fixture: it keeps the measured
  // 2*r*sin(delta/2) chord strictly below the default 1e-6 max_gap.
  const double target_chord = gap * (1.0 - 1e-8);
  const double delta = 2.0 * std::asin(target_chord / (2.0 * kMinorRadius));
  const double v_open_low = 0.5 * (v_cap_low + v_cap_high) - 0.5 * delta;
  const double v_open_high = 0.5 * (v_cap_low + v_cap_high) + 0.5 * delta;
  const double cap_half_height = kMinorRadius * std::sin(v_cap_high);

  SolidEditor editor;
  const std::uint32_t surface = editor.add_surface(torus);
  const std::uint32_t cap_low_plane =
      editor.add_surface(Plane3{{0.0, 0.0, -1.0}, cap_half_height});
  const std::uint32_t cap_high_plane =
      editor.add_surface(Plane3{{0.0, 0.0, 1.0}, cap_half_height});
  const std::uint32_t seam = editor.add_curve(seam_meridian(torus));

  const Rim rim_cap_low = add_rim(editor, v_rim_circle(torus, v_cap_low),
                                  torus.eval(0.0, v_cap_low),
                                  torus.eval(kPi, v_cap_low));
  const Rim rim_open_low = add_rim(editor, v_rim_circle(torus, v_open_low),
                                   torus.eval(0.0, v_open_low),
                                   torus.eval(kPi, v_open_low));
  const Rim rim_open_high = add_rim(editor, v_rim_circle(torus, v_open_high),
                                    torus.eval(0.0, v_open_high),
                                    torus.eval(kPi, v_open_high));
  const Rim rim_cap_high = add_rim(editor, v_rim_circle(torus, v_cap_high),
                                   torus.eval(0.0, v_cap_high),
                                   torus.eval(kPi, v_cap_high));
  const EdgeId seam_low = m4a_edge(editor, rim_cap_low.p0, rim_open_low.p0, seam,
                                   v_cap_low, v_open_low);
  const EdgeId seam_high =
      m4a_edge(editor, rim_open_high.p0, rim_cap_high.p0, seam, v_open_high,
               v_cap_high);

  const FaceId lower_lateral = m4a_face(editor, surface);
  const std::array<CoedgeSpec, 6> lower_loop{{
      {rim_cap_low.arc0, false}, {rim_cap_low.arc1, false},
      {seam_low, false}, {rim_open_low.arc1, true}, {rim_open_low.arc0, true},
      {seam_low, true}}};
  editor.add_loop(lower_lateral, lower_loop, true).value();

  const FaceId cap_low = m4a_face(editor, cap_low_plane);
  const std::array<CoedgeSpec, 2> cap_low_loop{{
      {rim_cap_low.arc0, true}, {rim_cap_low.arc1, true}}};
  editor.add_loop(cap_low, cap_low_loop, true).value();

  const FaceId upper_lateral = m4a_face(editor, surface);
  const std::array<CoedgeSpec, 6> upper_loop{{
      {rim_open_high.arc0, false}, {rim_open_high.arc1, false},
      {seam_high, false}, {rim_cap_high.arc1, true},
      {rim_cap_high.arc0, true}, {seam_high, true}}};
  editor.add_loop(upper_lateral, upper_loop, true).value();

  const FaceId cap_high = m4a_face(editor, cap_high_plane);
  const std::array<CoedgeSpec, 2> cap_high_loop{{
      {rim_cap_high.arc0, false}, {rim_cap_high.arc1, false}}};
  editor.add_loop(cap_high, cap_high_loop, true).value();

  const std::array<FaceId, 4> faces{
      {lower_lateral, cap_low, upper_lateral, cap_high}};
  const SolidId solid = m4a_add_solid(editor, faces);
  return CurvedSolid{editor.commit().value(), solid};
}

// Case 5c: the same barrel split, but around the INNER branch: the surface's
// chart seam (v_origin = pi) sits on the inner equator, so the band's chart v
// straddles 0 with the two halves on opposite sides. face_uv_project_torus
// returns v = atan2(axial, q - R) - v_origin without wrapping, so the upper
// half projects to v around -2*pi while the lower half projects to v around 0:
// a rebuilt pcurve must not jump between the two.
//
// Both caps face away from the band, which runs from z = +c down to z = -c, so
// the chart v direction is opposite to +z. The band's material is inside the
// tube, which on the inner branch puts the surface's outward (tube) normal on
// the material's inner side: the two lateral faces are reversed too.
// Volume, with rho(z) = R - sqrt(r^2 - z^2) on the inner branch:
//   V = pi * [2*R^2*c - 2*R*(c*s + r^2*asin(c/r)) + 2*r^2*c - (2/3)*c^3]
//     = 0.863505815719987... (the brief's constant is the same value from
// quadrature).
inline CurvedSolid make_torus_inner_v_seam_gap(double gap = 1e-6) {
  using namespace torus_fixture;
  const Torus3 torus{m4a_point(0.0, 0.0, 0.0), {0.0, 0.0, 1.0}, {1.0, 0.0, 0.0},
                     kMajorRadius, kMinorRadius, kPi};
  const double v_cap_low = -1.0;
  const double v_cap_high = 1.0;
  // Same 1e-8 relative margin as the outer-branch fixture, for the same reason.
  const double target_chord = gap * (1.0 - 1e-8);
  const double delta = 2.0 * std::asin(target_chord / (2.0 * kMinorRadius));
  const double v_open_low = 0.5 * (v_cap_low + v_cap_high) - 0.5 * delta;
  const double v_open_high = 0.5 * (v_cap_low + v_cap_high) + 0.5 * delta;
  // The chart -1 cap sits at z = r*sin(pi - 1) = +c; the chart +1 cap at -c.
  const double cap_half_height = kMinorRadius * std::sin(kPi - 1.0);

  SolidEditor editor;
  const std::uint32_t surface = editor.add_surface(torus);
  const std::uint32_t cap_low_plane =
      editor.add_surface(Plane3{{0.0, 0.0, 1.0}, cap_half_height});
  const std::uint32_t cap_high_plane =
      editor.add_surface(Plane3{{0.0, 0.0, -1.0}, cap_half_height});
  const std::uint32_t seam = editor.add_curve(seam_meridian(torus));

  const Rim rim_cap_low = add_rim(editor, v_rim_circle(torus, v_cap_low),
                                  torus.eval(0.0, v_cap_low),
                                  torus.eval(kPi, v_cap_low));
  const Rim rim_open_low = add_rim(editor, v_rim_circle(torus, v_open_low),
                                   torus.eval(0.0, v_open_low),
                                   torus.eval(kPi, v_open_low));
  const Rim rim_open_high = add_rim(editor, v_rim_circle(torus, v_open_high),
                                    torus.eval(0.0, v_open_high),
                                    torus.eval(kPi, v_open_high));
  const Rim rim_cap_high = add_rim(editor, v_rim_circle(torus, v_cap_high),
                                   torus.eval(0.0, v_cap_high),
                                   torus.eval(kPi, v_cap_high));
  const EdgeId seam_low = m4a_edge(editor, rim_cap_low.p0, rim_open_low.p0, seam,
                                   v_cap_low, v_open_low);
  const EdgeId seam_high =
      m4a_edge(editor, rim_open_high.p0, rim_cap_high.p0, seam, v_open_high,
               v_cap_high);

  const FaceId lower_lateral = m4a_face_reversed(editor, surface);
  const std::array<CoedgeSpec, 6> lower_loop{{
      {rim_cap_low.arc0, false}, {rim_cap_low.arc1, false},
      {seam_low, false}, {rim_open_low.arc1, true}, {rim_open_low.arc0, true},
      {seam_low, true}}};
  editor.add_loop(lower_lateral, lower_loop, true).value();

  const FaceId cap_low = m4a_face_reversed(editor, cap_low_plane);
  const std::array<CoedgeSpec, 2> cap_low_loop{{
      {rim_cap_low.arc0, true}, {rim_cap_low.arc1, true}}};
  editor.add_loop(cap_low, cap_low_loop, true).value();

  const FaceId upper_lateral = m4a_face_reversed(editor, surface);
  const std::array<CoedgeSpec, 6> upper_loop{{
      {rim_open_high.arc0, false}, {rim_open_high.arc1, false},
      {seam_high, false}, {rim_cap_high.arc1, true},
      {rim_cap_high.arc0, true}, {seam_high, true}}};
  editor.add_loop(upper_lateral, upper_loop, true).value();

  const FaceId cap_high = m4a_face_reversed(editor, cap_high_plane);
  const std::array<CoedgeSpec, 2> cap_high_loop{{
      {rim_cap_high.arc0, false}, {rim_cap_high.arc1, false}}};
  editor.add_loop(cap_high, cap_high_loop, true).value();

  const std::array<FaceId, 4> faces{
      {lower_lateral, cap_low, upper_lateral, cap_high}};
  const SolidId solid = m4a_add_solid(editor, faces);
  return CurvedSolid{editor.commit().value(), solid};
}

// Case 5d: two torus halves cut by the planes z = +/-gap/2. Unlike the band
// fixtures above, every boundary rim here is a full-turn circle with a single
// vertex (v0 == v1, parameter span 0..2*pi), so this case exercises the
// automatic closed-rim split before planning. The outer and inner rims are
// paired across the gap and welded. The split reuses the freed closed-rim edge
// handles, so the surviving (kept) outer rim is the lower one while the kept
// inner rim is the upper one; the healed upper face spans [-alpha, pi-alpha] and
// the lower face spans [pi-alpha, 2*pi-alpha], together covering the whole
// torus exactly once. Healed volume is the full torus
//   V = 2*pi^2*R*r^2 = pi^2/2 = 4.934802200544679...
//
// Both halves lie on the same Torus3 surface and their loops are the canonical
// (u, v) rectangles, so the shared rims are traversed in opposite directions
// across the gap.
inline CurvedSolid make_torus_closed_rims_gap(double gap = 1e-6) {
  using namespace torus_fixture;
  const Torus3 torus{m4a_point(0.0, 0.0, 0.0), {0.0, 0.0, 1.0}, {1.0, 0.0, 0.0},
                     kMajorRadius, kMinorRadius, 0.0};
  // Same 1e-8 relative margin as the open torus band fixtures: the cut height
  // is recovered through asin/sin around the tube, so the margin keeps the
  // cross-gap vertex chord strictly below the default 1e-6 max_gap.
  const double target_half_height = 0.5 * gap * (1.0 - 1e-8);
  const double alpha = std::asin(target_half_height / kMinorRadius);
  const double v_outer_high = alpha;
  const double v_inner_high = kPi - alpha;
  const double v_inner_low = kPi + alpha;
  const double v_outer_low = kTwoPi - alpha;

  SolidEditor editor;
  const std::uint32_t surface = editor.add_surface(torus);
  const std::uint32_t seam = editor.add_curve(seam_meridian(torus));

  // Closed rims: one vertex at u = 0 and one closed edge per rim. The split
  // pass later adds the opposite u = pi vertex and the two half-turn arcs.
  struct ClosedRim {
    VertexId vertex{};
    EdgeId edge{};
  };
  const auto add_closed_rim = [&](double chart_v) {
    const std::uint32_t curve = editor.add_curve(v_rim_circle(torus, chart_v));
    const VertexId vertex = m4a_vertex(editor, torus.eval(0.0, chart_v));
    return ClosedRim{vertex,
                     m4a_edge(editor, vertex, vertex, curve, 0.0, kTwoPi)};
  };
  const ClosedRim rim_outer_high = add_closed_rim(v_outer_high);
  const ClosedRim rim_inner_high = add_closed_rim(v_inner_high);
  const ClosedRim rim_inner_low = add_closed_rim(v_inner_low);
  const ClosedRim rim_outer_low = add_closed_rim(v_outer_low);

  const EdgeId seam_high =
      m4a_edge(editor, rim_outer_high.vertex, rim_inner_high.vertex, seam,
               v_outer_high, v_inner_high);
  const EdgeId seam_low =
      m4a_edge(editor, rim_inner_low.vertex, rim_outer_low.vertex, seam,
               v_inner_low, v_outer_low);

  const FaceId high = m4a_face(editor, surface);
  const std::array<CoedgeSpec, 4> high_loop{{
      {rim_outer_high.edge, false}, {seam_high, false},
      {rim_inner_high.edge, true}, {seam_high, true}}};
  editor.add_loop(high, high_loop, true).value();

  const FaceId low = m4a_face(editor, surface);
  const std::array<CoedgeSpec, 4> low_loop{{
      {rim_inner_low.edge, false}, {seam_low, false},
      {rim_outer_low.edge, true}, {seam_low, true}}};
  editor.add_loop(low, low_loop, true).value();

  const std::array<FaceId, 2> faces{{high, low}};
  const SolidId solid = m4a_add_solid(editor, faces);
  return CurvedSolid{editor.commit().value(), solid};
}

// Case 5e: the upper torus half of case 5d capped by a plane disk at the same
// height, with the disk rim coincident with the half's outer closed rim. The
// two rims pair at zero gap, so the surface guard must reject the pair because
// a Torus3 face and a Plane3 face cannot be the same carrier. The inner rim has
// no counterpart, which would report gap_too_large; the surface guard runs
// first on the coincident outer pair.
inline CurvedSolid make_torus_plane_lid_gap(double gap = 1e-6) {
  using namespace torus_fixture;
  const Torus3 torus{m4a_point(0.0, 0.0, 0.0), {0.0, 0.0, 1.0}, {1.0, 0.0, 0.0},
                     kMajorRadius, kMinorRadius, 0.0};
  const double target_half_height = 0.5 * gap * (1.0 - 1e-8);
  const double alpha = std::asin(target_half_height / kMinorRadius);
  const double outer_radius = kMajorRadius + kMinorRadius * std::cos(alpha);

  SolidEditor editor;
  const std::uint32_t torus_surface = editor.add_surface(torus);
  const std::uint32_t plane_surface =
      editor.add_surface(Plane3{{0.0, 0.0, 1.0}, target_half_height});
  const std::uint32_t seam = editor.add_curve(seam_meridian(torus));

  const auto add_closed_rim = [&](double chart_v) {
    const std::uint32_t curve = editor.add_curve(v_rim_circle(torus, chart_v));
    const VertexId vertex = m4a_vertex(editor, torus.eval(0.0, chart_v));
    return std::pair<VertexId, EdgeId>{
        vertex, m4a_edge(editor, vertex, vertex, curve, 0.0, kTwoPi)};
  };
  const auto [outer_vertex, outer_rim] = add_closed_rim(alpha);
  const auto [inner_vertex, inner_rim] = add_closed_rim(kPi - alpha);
  const EdgeId seam_edge = m4a_edge(editor, outer_vertex, inner_vertex, seam,
                                    alpha, kPi - alpha);

  const FaceId upper = m4a_face(editor, torus_surface);
  const std::array<CoedgeSpec, 4> upper_loop{{
      {outer_rim, false}, {seam_edge, false}, {inner_rim, true},
      {seam_edge, true}}};
  editor.add_loop(upper, upper_loop, true).value();

  // The plane lid is a disk whose boundary is the torus half's outer rim.
  const std::uint32_t lid_curve = editor.add_curve(
      CircleCurve3{m4a_point(0.0, 0.0, target_half_height),
                   {1.0, 0.0, 0.0}, {0.0, 1.0, 0.0}, outer_radius});
  const VertexId lid_vertex =
      m4a_vertex(editor, m4a_point(outer_radius, 0.0, target_half_height));
  const EdgeId lid_rim =
      m4a_edge(editor, lid_vertex, lid_vertex, lid_curve, 0.0, kTwoPi);
  const FaceId lid = m4a_face(editor, plane_surface);
  const std::array<CoedgeSpec, 1> lid_loop{{{lid_rim, false}}};
  editor.add_loop(lid, lid_loop, true).value();

  const std::array<FaceId, 2> faces{{upper, lid}};
  const SolidId solid = m4a_add_solid(editor, faces);
  return CurvedSolid{editor.commit().value(), solid};
}

// Case 5f: the two halves of case 5d, but the lower outer closed rim is a
// circle of radius R + r*cos(alpha) + radial_shift, so its carrier no longer
// lies on the shared Torus3 surface. The cut height keeps a 1e-2 relative
// margin instead of 1e-8 so the cross-gap chord stays below the default 1e-6
// max_gap even with the radial shift; planning then pairs the rims and reaches
// the off-surface guard instead of reporting gap_too_large.
inline CurvedSolid make_torus_off_surface_gap(double gap = 1e-6,
                                              double radial_shift = 1e-7) {
  using namespace torus_fixture;
  const Torus3 torus{m4a_point(0.0, 0.0, 0.0), {0.0, 0.0, 1.0}, {1.0, 0.0, 0.0},
                     kMajorRadius, kMinorRadius, 0.0};
  const double target_half_height = 0.5 * gap * (1.0 - 1e-2);
  const double alpha = std::asin(target_half_height / kMinorRadius);
  const double outer_radius = kMajorRadius + kMinorRadius * std::cos(alpha);

  SolidEditor editor;
  const std::uint32_t surface = editor.add_surface(torus);
  const std::uint32_t seam = editor.add_curve(seam_meridian(torus));

  const auto add_closed_rim = [&](const CircleCurve3& circle,
                                  const Point3d& at) {
    const std::uint32_t curve = editor.add_curve(circle);
    const VertexId vertex = m4a_vertex(editor, at);
    return std::pair<VertexId, EdgeId>{
        vertex, m4a_edge(editor, vertex, vertex, curve, 0.0, kTwoPi)};
  };
  const auto [outer_high_vertex, outer_high] =
      add_closed_rim(v_rim_circle(torus, alpha), torus.eval(0.0, alpha));
  const auto [inner_high_vertex, inner_high] = add_closed_rim(
      v_rim_circle(torus, kPi - alpha), torus.eval(0.0, kPi - alpha));
  const auto [inner_low_vertex, inner_low] = add_closed_rim(
      v_rim_circle(torus, kPi + alpha), torus.eval(0.0, kPi + alpha));
  // The lower outer rim is pushed radially outward, off the torus surface.
  const CircleCurve3 outer_low_circle{
      m4a_point(0.0, 0.0, -target_half_height), {1.0, 0.0, 0.0},
      {0.0, 1.0, 0.0}, outer_radius + radial_shift};
  const auto [outer_low_vertex, outer_low] =
      add_closed_rim(outer_low_circle,
                     m4a_point(outer_radius + radial_shift, 0.0,
                               -target_half_height));

  const EdgeId seam_high = m4a_edge(editor, outer_high_vertex, inner_high_vertex,
                                    seam, alpha, kPi - alpha);
  const EdgeId seam_low = m4a_edge(editor, inner_low_vertex, outer_low_vertex,
                                   seam, kPi + alpha, kTwoPi - alpha);

  const FaceId high = m4a_face(editor, surface);
  const std::array<CoedgeSpec, 4> high_loop{{
      {outer_high, false}, {seam_high, false}, {inner_high, true},
      {seam_high, true}}};
  editor.add_loop(high, high_loop, true).value();

  const FaceId low = m4a_face(editor, surface);
  const std::array<CoedgeSpec, 4> low_loop{{
      {inner_low, false}, {seam_low, false}, {outer_low, true},
      {seam_low, true}}};
  editor.add_loop(low, low_loop, true).value();

  const std::array<FaceId, 2> faces{{high, low}};
  const SolidId solid = m4a_add_solid(editor, faces);
  return CurvedSolid{editor.commit().value(), solid};
}

// ---------------------------------------------------------------------------
// Negative kind fixtures: a closed unit cube.
//
// The bottom face carries a flat NurbsSurface in make_nurbs_face_solid(); in
// make_nurbs_edge_solid() every face is planar and edge b0 carries a NURBS
// curve from (0,0,0) to (1,0,0). One edge is traversed in each direction, so
// the topology commits as a closed shell.
// ---------------------------------------------------------------------------
namespace detail {

struct CubeFixture {
  Brep brep;
  SolidId solid{};
  FaceId bottom{};
};

inline CubeFixture build_unit_cube(bool nurbs_face, bool nurbs_edge) {
  SolidEditor editor;
  const std::uint32_t plane =
      editor.add_surface(Plane3{{0.0, 0.0, 1.0}, 0.0});
  const std::uint32_t bottom_surface =
      nurbs_face
          ? editor.add_surface(NurbsSurface(
                1, 1, {0.0, 0.0, 1.0, 1.0}, {0.0, 0.0, 1.0, 1.0},
                {{0.0, 0.0, 0.0}, {1.0, 0.0, 0.0}, {0.0, 1.0, 0.0},
                 {1.0, 1.0, 0.0}},
                {1.0, 1.0, 1.0, 1.0}))
          : plane;
  const std::uint32_t b0_curve =
      nurbs_edge
          ? editor.add_curve(NurbsCurve(1, {0.0, 0.0, 1.0, 1.0},
                                        {{0.0, 0.0, 0.0}, {1.0, 0.0, 0.0}},
                                        {1.0, 1.0}))
          : kNoGeometry;

  const VertexId v000 = m4a_vertex(editor, m4a_point(0.0, 0.0, 0.0));
  const VertexId v100 = m4a_vertex(editor, m4a_point(1.0, 0.0, 0.0));
  const VertexId v110 = m4a_vertex(editor, m4a_point(1.0, 1.0, 0.0));
  const VertexId v010 = m4a_vertex(editor, m4a_point(0.0, 1.0, 0.0));
  const VertexId v001 = m4a_vertex(editor, m4a_point(0.0, 0.0, 1.0));
  const VertexId v101 = m4a_vertex(editor, m4a_point(1.0, 0.0, 1.0));
  const VertexId v111 = m4a_vertex(editor, m4a_point(1.0, 1.0, 1.0));
  const VertexId v011 = m4a_vertex(editor, m4a_point(0.0, 1.0, 1.0));

  const EdgeId b0 = m4a_edge(editor, v000, v100, b0_curve, 0.0, 1.0);
  const EdgeId b1 = m4a_edge(editor, v100, v110, kNoGeometry, 0.0, 1.0);
  const EdgeId b2 = m4a_edge(editor, v110, v010, kNoGeometry, 0.0, 1.0);
  const EdgeId b3 = m4a_edge(editor, v010, v000, kNoGeometry, 0.0, 1.0);
  const EdgeId t0 = m4a_edge(editor, v001, v101, kNoGeometry, 0.0, 1.0);
  const EdgeId t1 = m4a_edge(editor, v101, v111, kNoGeometry, 0.0, 1.0);
  const EdgeId t2 = m4a_edge(editor, v111, v011, kNoGeometry, 0.0, 1.0);
  const EdgeId t3 = m4a_edge(editor, v011, v001, kNoGeometry, 0.0, 1.0);
  const EdgeId s0 = m4a_edge(editor, v000, v001, kNoGeometry, 0.0, 1.0);
  const EdgeId s1 = m4a_edge(editor, v100, v101, kNoGeometry, 0.0, 1.0);
  const EdgeId s2 = m4a_edge(editor, v110, v111, kNoGeometry, 0.0, 1.0);
  const EdgeId s3 = m4a_edge(editor, v010, v011, kNoGeometry, 0.0, 1.0);

  const FaceId top = m4a_face(editor, plane);
  const std::array<CoedgeSpec, 4> top_loop{{
      {t0, false}, {t1, false}, {t2, false}, {t3, false}}};
  editor.add_loop(top, top_loop, true).value();

  const FaceId x1 = m4a_face(editor, plane);
  const std::array<CoedgeSpec, 4> x1_loop{{
      {b1, false}, {s2, false}, {t1, true}, {s1, true}}};
  editor.add_loop(x1, x1_loop, true).value();

  const FaceId y1 = m4a_face(editor, plane);
  const std::array<CoedgeSpec, 4> y1_loop{{
      {s3, false}, {t2, true}, {s2, true}, {b2, false}}};
  editor.add_loop(y1, y1_loop, true).value();

  const FaceId x0 = m4a_face(editor, plane);
  const std::array<CoedgeSpec, 4> x0_loop{{
      {s0, false}, {t3, true}, {s3, true}, {b3, false}}};
  editor.add_loop(x0, x0_loop, true).value();

  const FaceId y0 = m4a_face(editor, plane);
  const std::array<CoedgeSpec, 4> y0_loop{{
      {b0, false}, {s1, false}, {t0, true}, {s0, true}}};
  editor.add_loop(y0, y0_loop, true).value();

  // Added last so the face index a guard reports is non-zero; the shell order
  // below still puts it first, so the surface guard reaches it first.
  const FaceId bottom = m4a_face(editor, bottom_surface);
  const std::array<CoedgeSpec, 4> bottom_loop{{
      {b3, true}, {b2, true}, {b1, true}, {b0, true}}};
  editor.add_loop(bottom, bottom_loop, true).value();

  const std::array<FaceId, 6> faces{{bottom, top, x1, y1, x0, y0}};
  const SolidId solid = m4a_add_solid(editor, faces);
  return CubeFixture{editor.commit().value(), solid, bottom};
}

}  // namespace detail

struct UnsupportedFaceSolid {
  Brep brep;
  SolidId solid{};
  FaceId face{};
};

struct UnsupportedCurveSolid {
  Brep brep;
  SolidId solid{};
};

inline UnsupportedFaceSolid make_nurbs_face_solid() {
  detail::CubeFixture cube = detail::build_unit_cube(true, false);
  return UnsupportedFaceSolid{std::move(cube.brep), cube.solid, cube.bottom};
}

inline UnsupportedCurveSolid make_nurbs_edge_solid() {
  detail::CubeFixture cube = detail::build_unit_cube(false, true);
  return UnsupportedCurveSolid{std::move(cube.brep), cube.solid};
}

}  // namespace geom::test::m4a
