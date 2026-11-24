#pragma once

// Fixtures for the M4b-2 sphere-pole healing tests.
//
// The fixture is ONE solid: the upper hemisphere of the unit sphere, closed by
// the z = 0 disk, with the narrow band v in [v_up, v_cap] removed at the north
// pole. The band's pole rim sits `delta` below the pole and the cap's rim sits
// `gap` above it, so the two rims are a small-latitude pair on the SAME sphere.
// The cap ends exactly at the pole (v = pi/2) and its meridian seam is used
// twice, once in each direction, which is the collapsed-pole chart the M4b-2
// projection must survive.
//
// All circles are split into two half-turn arcs at u = 0 and u = pi, so
// heal_plan never sees a closed (zero chord length) boundary edge and no
// closed-rim split is needed. The solid's volume is 2*pi/3.

#include <array>
#include <cmath>
#include <cstdint>

#include "m4a_curved_heal_fixture.hpp"

namespace geom::test::m4b {

// Sphere radius 1 with the band [v_up, v_up + gap] removed at the north pole,
// where v_up = pi/2 - delta. Vertices are created e0, e1, p0, p1, c0, c1, pole;
// the gap pair is p0/p1 on the band rim and c0/c1 on the cap rim.
inline m4a::CurvedSolid make_sphere_pole_gap(double delta, double gap) {
  using namespace geom::test::m4a;

  SolidEditor editor;
  const double radius = 1.0;
  const std::uint32_t sphere =
      editor.add_surface(Sphere3{m4a_point(0.0, 0.0, 0.0), {0.0, 0.0, 1.0},
                                 {1.0, 0.0, 0.0}, radius});
  const std::uint32_t eq_low = editor.add_curve(
      CircleCurve3{m4a_point(0.0, 0.0, 0.0), {1.0, 0.0, 0.0}, {0.0, 1.0, 0.0}, radius});
  const double v_up = kPi * 0.5 - delta;
  const std::uint32_t pole_rim = editor.add_curve(CircleCurve3{
      m4a_point(0.0, 0.0, radius * std::sin(v_up)), {1.0, 0.0, 0.0}, {0.0, 1.0, 0.0},
      radius * std::cos(v_up)});
  const double v_cap = v_up + gap;
  const std::uint32_t cap_rim = editor.add_curve(CircleCurve3{
      m4a_point(0.0, 0.0, radius * std::sin(v_cap)), {1.0, 0.0, 0.0}, {0.0, 1.0, 0.0},
      radius * std::cos(v_cap)});
  const std::uint32_t mer = editor.add_curve(CircleCurve3{
      m4a_point(0.0, 0.0, 0.0), {1.0, 0.0, 0.0}, {0.0, 0.0, 1.0}, radius});
  const std::uint32_t plane = editor.add_surface(Plane3{{0.0, 0.0, -1.0}, 0.0});

  const VertexId e0 = m4a_vertex(editor, m4a_point(radius, 0.0, 0.0));
  const VertexId e1 = m4a_vertex(editor, m4a_point(-radius, 0.0, 0.0));
  const VertexId p0 = m4a_vertex(editor,
                                 m4a_point(radius * std::cos(v_up), 0.0, radius * std::sin(v_up)));
  const VertexId p1 = m4a_vertex(editor,
                                 m4a_point(-radius * std::cos(v_up), 0.0, radius * std::sin(v_up)));
  const VertexId c0 = m4a_vertex(editor,
                                 m4a_point(radius * std::cos(v_cap), 0.0, radius * std::sin(v_cap)));
  const VertexId c1 = m4a_vertex(editor,
                                 m4a_point(-radius * std::cos(v_cap), 0.0, radius * std::sin(v_cap)));
  const VertexId pole = m4a_vertex(editor, m4a_point(0.0, 0.0, radius));

  const EdgeId eqa = m4a_edge(editor, e0, e1, eq_low, 0.0, kPi);
  const EdgeId eqb = m4a_edge(editor, e1, e0, eq_low, kPi, kTwoPi);
  const EdgeId pola = m4a_edge(editor, p0, p1, pole_rim, 0.0, kPi);
  const EdgeId polb = m4a_edge(editor, p1, p0, pole_rim, kPi, kTwoPi);
  const EdgeId capa = m4a_edge(editor, c0, c1, cap_rim, 0.0, kPi);
  const EdgeId capb = m4a_edge(editor, c1, c0, cap_rim, kPi, kTwoPi);
  const EdgeId seam_low = m4a_edge(editor, e0, p0, mer, 0.0, v_up);
  const EdgeId seam_high = m4a_edge(editor, c0, pole, mer, v_cap, kPi * 0.5);

  const FaceId band = m4a_face(editor, sphere);
  const std::array<CoedgeSpec, 6> band_loop{{
      {eqa, false}, {eqb, false}, {seam_low, false},
      {polb, true}, {pola, true}, {seam_low, true}}};
  editor.add_loop(band, band_loop, true).value();
  const FaceId cap = m4a_face(editor, plane);
  const std::array<CoedgeSpec, 2> cap_loop{{{eqa, true}, {eqb, true}}};
  editor.add_loop(cap, cap_loop, true).value();
  const FaceId pole_cap = m4a_face(editor, sphere);
  const std::array<CoedgeSpec, 4> pole_loop{{
      {capa, false}, {capb, false}, {seam_high, false}, {seam_high, true}}};
  editor.add_loop(pole_cap, pole_loop, true).value();

  const std::array<FaceId, 3> faces{{band, cap, pole_cap}};
  const SolidId solid = m4a_add_solid(editor, faces);
  return m4a::CurvedSolid{editor.commit().value(), solid};
}

// ---------------------------------------------------------------------------
// Cone fixtures for the dissolve tests.
//
// Both builders use the same frame: base radius 1, half angle atan(1/2), so the
// apex sits two units above z_base. Only the apex handling differs.
// ---------------------------------------------------------------------------

// Cone with a pinhole at the apex: the top rim has radius r0 (r0 > 0) and is
// either two half-turn arcs (split_rim) or one closed arc. The lateral loop
// uses `seam` twice, once in each direction; the base face closes the cone at
// z_base. The shell stays open, so the zero-length-by-limit top rim is a legal
// boundary of the sheet.
struct ConePinhole {
  EdgeId base0{}, base1{}, seam{}, top_a{}, top_b{};
  VertexId b0{}, b1{}, t0{}, t1{};
  FaceId lateral{}, base{};
  LoopId lateral_loop{}, base_loop{};
  ShellId shell{};
  SolidId solid{};
  bool split_rim = true;
};

inline ConePinhole add_cone_pinhole(SolidEditor& editor, double z_base, double r0,
                                    bool split_rim) {
  using namespace geom::test::m4a;

  ConePinhole cone;
  cone.split_rim = split_rim;
  const std::uint32_t cone_surface = editor.add_surface(Cone3{
      m4a_point(0.0, 0.0, z_base), {0.0, 0.0, 1.0}, {1.0, 0.0, 0.0}, 1.0, std::atan(0.5)});
  const std::uint32_t base_plane = editor.add_surface(Plane3{{0.0, 0.0, -1.0}, z_base});
  const std::uint32_t base_circle = editor.add_curve(CircleCurve3{
      m4a_point(0.0, 0.0, z_base), {1.0, 0.0, 0.0}, {0.0, 1.0, 0.0}, 1.0});
  // v = radius / tan(half_angle) = 2 is the apex; trimming at r0 keeps v_top.
  const double z_top = z_base + 2.0 - 2.0 * r0;
  const std::uint32_t top_circle = editor.add_curve(CircleCurve3{
      m4a_point(0.0, 0.0, z_top), {1.0, 0.0, 0.0}, {0.0, 1.0, 0.0}, r0});
  const double dx = r0 - 1.0;
  const double dz = z_top - z_base;
  const double seam_length = std::sqrt(dx * dx + dz * dz);
  const std::uint32_t seam_curve = editor.add_curve(LineCurve3{
      m4a_point(1.0, 0.0, z_base), {dx / seam_length, 0.0, dz / seam_length}});

  cone.b0 = m4a_vertex(editor, m4a_point(1.0, 0.0, z_base));
  cone.b1 = m4a_vertex(editor, m4a_point(-1.0, 0.0, z_base));
  cone.t0 = m4a_vertex(editor, m4a_point(r0, 0.0, z_top));
  cone.t1 = m4a_vertex(editor, m4a_point(-r0, 0.0, z_top));

  cone.base0 = m4a_edge(editor, cone.b0, cone.b1, base_circle, 0.0, kPi);
  cone.base1 = m4a_edge(editor, cone.b1, cone.b0, base_circle, kPi, kTwoPi);
  cone.seam = m4a_edge(editor, cone.b0, cone.t0, seam_curve, 0.0, seam_length);
  cone.top_a = split_rim
                   ? m4a_edge(editor, cone.t0, cone.t1, top_circle, 0.0, kPi)
                   : m4a_edge(editor, cone.t0, cone.t0, top_circle, 0.0, kTwoPi);
  if (split_rim) {
    cone.top_b = m4a_edge(editor, cone.t1, cone.t0, top_circle, kPi, kTwoPi);
  }

  cone.lateral = m4a_face(editor, cone_surface);
  cone.base = m4a_face(editor, base_plane);
  if (split_rim) {
    const std::array<CoedgeSpec, 6> loop{{
        {cone.base0, false}, {cone.base1, false}, {cone.seam, false},
        {cone.top_a, false}, {cone.top_b, false}, {cone.seam, true}}};
    cone.lateral_loop = editor.add_loop(cone.lateral, loop, true).value();
  } else {
    const std::array<CoedgeSpec, 5> loop{{
        {cone.base0, false}, {cone.base1, false}, {cone.seam, false},
        {cone.top_a, false}, {cone.seam, true}}};
    cone.lateral_loop = editor.add_loop(cone.lateral, loop, true).value();
  }
  const std::array<CoedgeSpec, 2> base_loop{{{cone.base0, true}, {cone.base1, true}}};
  cone.base_loop = editor.add_loop(cone.base, base_loop, true).value();
  cone.shell =
      editor.add_shell(std::array<FaceId, 2>{cone.lateral, cone.base}, false).value();
  cone.solid = editor.add_solid(std::array<ShellId, 1>{cone.shell}, 0.0).value();
  return cone;
}

// Cone pinched shut at the apex: the zero-length top rim is gone, so the apex is
// ONE vertex and the seam is used twice in the same lateral loop. The shell is
// closed and the result is a valid two-face solid of volume 2*pi/3 (the form
// dissolve + point weld must produce).
struct ConePinch {
  EdgeId base0{}, base1{}, seam{};
  VertexId b0{}, b1{}, apex{};
  FaceId lateral{}, base{};
  LoopId lateral_loop{}, base_loop{};
  ShellId shell{};
  SolidId solid{};
};

inline ConePinch add_cone_pinch(SolidEditor& editor, double z_base) {
  using namespace geom::test::m4a;

  ConePinch cone;
  const std::uint32_t cone_surface = editor.add_surface(Cone3{
      m4a_point(0.0, 0.0, z_base), {0.0, 0.0, 1.0}, {1.0, 0.0, 0.0}, 1.0, std::atan(0.5)});
  const std::uint32_t base_plane = editor.add_surface(Plane3{{0.0, 0.0, -1.0}, z_base});
  const std::uint32_t base_circle = editor.add_curve(CircleCurve3{
      m4a_point(0.0, 0.0, z_base), {1.0, 0.0, 0.0}, {0.0, 1.0, 0.0}, 1.0});
  const std::uint32_t seam_curve =
      editor.add_curve(LineCurve3{m4a_point(1.0, 0.0, z_base),
                                  {-1.0 / std::sqrt(5.0), 0.0, 2.0 / std::sqrt(5.0)}});

  cone.b0 = m4a_vertex(editor, m4a_point(1.0, 0.0, z_base));
  cone.b1 = m4a_vertex(editor, m4a_point(-1.0, 0.0, z_base));
  cone.apex = m4a_vertex(editor, m4a_point(0.0, 0.0, z_base + 2.0));

  cone.base0 = m4a_edge(editor, cone.b0, cone.b1, base_circle, 0.0, kPi);
  cone.base1 = m4a_edge(editor, cone.b1, cone.b0, base_circle, kPi, kTwoPi);
  cone.seam = m4a_edge(editor, cone.b0, cone.apex, seam_curve, 0.0, std::sqrt(5.0));

  cone.lateral = m4a_face(editor, cone_surface);
  cone.base = m4a_face(editor, base_plane);
  const std::array<CoedgeSpec, 4> loop{{
      {cone.base0, false}, {cone.base1, false}, {cone.seam, false}, {cone.seam, true}}};
  cone.lateral_loop = editor.add_loop(cone.lateral, loop, true).value();
  const std::array<CoedgeSpec, 2> base_loop{{{cone.base0, true}, {cone.base1, true}}};
  cone.base_loop = editor.add_loop(cone.base, base_loop, true).value();
  cone.shell =
      editor.add_shell(std::array<FaceId, 2>{cone.lateral, cone.base}, true).value();
  cone.solid = editor.add_solid(std::array<ShellId, 1>{cone.shell}, 0.0).value();
  return cone;
}

}  // namespace geom::test::m4b
