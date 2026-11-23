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
inline CurvedSolid make_cylinder_bands(double radius, double gap,
                                       double upper_radius, bool off_surface,
                                       bool shifted_chart) {
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
  const std::uint32_t circle_lower_top = editor.add_curve(
      CircleCurve3{m4a_point(0.0, 0.0, 1.0), {1.0, 0.0, 0.0}, {0.0, 1.0, 0.0}, radius});
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
      {s2, false}, {t1, true}, {s1, true}, {b1, true}}};
  editor.add_loop(x1, x1_loop, true).value();

  const FaceId y1 = m4a_face(editor, plane);
  const std::array<CoedgeSpec, 4> y1_loop{{
      {s3, false}, {t2, true}, {s2, true}, {b2, true}}};
  editor.add_loop(y1, y1_loop, true).value();

  const FaceId x0 = m4a_face(editor, plane);
  const std::array<CoedgeSpec, 4> x0_loop{{
      {s0, false}, {t3, true}, {s3, true}, {b3, true}}};
  editor.add_loop(x0, x0_loop, true).value();

  const FaceId y0 = m4a_face(editor, plane);
  const std::array<CoedgeSpec, 4> y0_loop{{
      {s1, false}, {t0, true}, {s0, true}, {b0, true}}};
  editor.add_loop(y0, y0_loop, true).value();

  // Added last so the face index a guard reports is non-zero; the shell order
  // below still puts it first, so the surface guard reaches it first.
  const FaceId bottom = m4a_face(editor, bottom_surface);
  const std::array<CoedgeSpec, 4> bottom_loop{{
      {b0, false}, {b1, false}, {b2, false}, {b3, false}}};
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
