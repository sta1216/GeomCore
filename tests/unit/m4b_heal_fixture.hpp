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

}  // namespace geom::test::m4b
