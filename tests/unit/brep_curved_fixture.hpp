#pragma once

// Curved B-rep fixtures for the M3b-2 curved boolean tasks. The quarter
// cylinder is a closed solid with five faces: one lateral Cylinder3 face and
// four planar faces (bottom and top quarter disks, and the two flat sides at
// u = 0 and u = pi/2). Both circular arcs are approximated with sixteen chords,
// and every coedge carries a pcurve in its own surface parameters.

#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <numbers>
#include <vector>

#include <geom/brep/solid_editor.hpp>
#include <geom/brep/topology.hpp>
#include <geom/core/error.hpp>
#include <geom/curves/analytic_curve.hpp>
#include <geom/primitives/point.hpp>
#include <geom/surfaces/analytic.hpp>
#include <geom/surfaces/plane.hpp>

namespace geom::test {

inline Brep make_quarter_cylinder(double radius, double height, double tol) {
  constexpr std::size_t kSegments = 16;
  const double quarter = std::numbers::pi_v<double> / 2.0;
  const double delta = quarter / static_cast<double>(kSegments);

  SolidEditor editor;
  const std::uint32_t lateral_surface = editor.add_surface(
      Cylinder3{Point3d{0.0, 0.0, 0.0}, Vec3d{0.0, 0.0, 1.0},
                Vec3d{1.0, 0.0, 0.0}, radius});
  const std::uint32_t bottom_surface =
      editor.add_surface(Plane3{Vec3d{0.0, 0.0, -1.0}, 0.0});
  const std::uint32_t top_surface =
      editor.add_surface(Plane3{Vec3d{0.0, 0.0, 1.0}, height});
  const std::uint32_t side0_surface =
      editor.add_surface(Plane3{Vec3d{0.0, -1.0, 0.0}, 0.0});
  const std::uint32_t side1_surface =
      editor.add_surface(Plane3{Vec3d{-1.0, 0.0, 0.0}, 0.0});

  std::array<Point3d, kSegments + 1> bottom{};
  std::array<Point3d, kSegments + 1> top{};
  std::array<double, kSegments + 1> angle{};
  for (std::size_t i = 0; i <= kSegments; ++i) {
    angle[i] = delta * static_cast<double>(i);
    const double c = std::cos(angle[i]);
    const double s = std::sin(angle[i]);
    bottom[i] = Point3d{radius * c, radius * s, 0.0};
    top[i] = Point3d{radius * c, radius * s, height};
  }
  const Point3d origin0{0.0, 0.0, 0.0};
  const Point3d origin1{0.0, 0.0, height};

  std::array<VertexId, kSegments + 1> bottom_vertex{};
  std::array<VertexId, kSegments + 1> top_vertex{};
  for (std::size_t i = 0; i <= kSegments; ++i) {
    bottom_vertex[i] = editor.add_vertex(bottom[i], 0.0).value();
    top_vertex[i] = editor.add_vertex(top[i], 0.0).value();
  }
  const VertexId origin0_vertex = editor.add_vertex(origin0, 0.0).value();
  const VertexId origin1_vertex = editor.add_vertex(origin1, 0.0).value();

  // Parameter maps of the four planar faces, expressed in the same frames that
  // make_plane_frame() derives from their normals.
  const auto uv_bottom = [](const Point3d& p) {
    return Point2d{p[1], p[0]};
  };
  const auto uv_top = [](const Point3d& p) {
    return Point2d{-p[1], p[0]};
  };
  const auto uv_side0 = [](const Point3d& p) {
    return Point2d{-p[2], p[0]};
  };
  const auto uv_side1 = [](const Point3d& p) {
    return Point2d{p[2], p[1]};
  };

  std::array<EdgeId, kSegments> bottom_arc{};
  std::array<EdgeId, kSegments> top_arc{};
  std::array<std::uint32_t, kSegments> bottom_arc_bottom_pcurve{};
  std::array<std::uint32_t, kSegments> bottom_arc_lateral_pcurve{};
  std::array<std::uint32_t, kSegments> top_arc_top_pcurve{};
  std::array<std::uint32_t, kSegments> top_arc_lateral_pcurve{};
  for (std::size_t i = 0; i < kSegments; ++i) {
    bottom_arc[i] =
        editor
            .add_edge(bottom_vertex[i], bottom_vertex[i + 1],
                      editor.add_curve(LineCurve3{
                          bottom[i], Vec3d{bottom[i + 1][0] - bottom[i][0],
                                           bottom[i + 1][1] - bottom[i][1],
                                           bottom[i + 1][2] - bottom[i][2]}}),
                      0.0, 1.0, 0.0)
            .value();
    top_arc[i] =
        editor
            .add_edge(top_vertex[i], top_vertex[i + 1],
                      editor.add_curve(LineCurve3{
                          top[i], Vec3d{top[i + 1][0] - top[i][0],
                                        top[i + 1][1] - top[i][1],
                                        top[i + 1][2] - top[i][2]}}),
                      0.0, 1.0, 0.0)
            .value();
    bottom_arc_bottom_pcurve[i] = editor.add_pcurve(Pcurve2d{
        {uv_bottom(bottom[i]), uv_bottom(bottom[i + 1])}});
    bottom_arc_lateral_pcurve[i] = editor.add_pcurve(
        Pcurve2d{{Point2d{angle[i], 0.0}, Point2d{angle[i + 1], 0.0}}});
    top_arc_top_pcurve[i] = editor.add_pcurve(
        Pcurve2d{{uv_top(top[i]), uv_top(top[i + 1])}});
    top_arc_lateral_pcurve[i] = editor.add_pcurve(
        Pcurve2d{{Point2d{angle[i], height}, Point2d{angle[i + 1], height}}});
  }

  const auto radial_edge = [&editor](const Point3d& from, const Point3d& to,
                                     VertexId v0, VertexId v1) {
    const std::uint32_t curve = editor.add_curve(LineCurve3{
        from, Vec3d{to[0] - from[0], to[1] - from[1], to[2] - from[2]}});
    return editor.add_edge(v0, v1, curve, 0.0, 1.0, 0.0).value();
  };

  const EdgeId radial0 =
      radial_edge(origin0, bottom[0], origin0_vertex, bottom_vertex[0]);
  const EdgeId radial1 =
      radial_edge(origin0, bottom[kSegments], origin0_vertex,
                  bottom_vertex[kSegments]);
  const EdgeId radial0_top =
      radial_edge(origin1, top[0], origin1_vertex, top_vertex[0]);
  const EdgeId radial1_top =
      radial_edge(origin1, top[kSegments], origin1_vertex,
                  top_vertex[kSegments]);
  const EdgeId side0_line =
      radial_edge(bottom[0], top[0], bottom_vertex[0], top_vertex[0]);
  const EdgeId side1_line =
      radial_edge(bottom[kSegments], top[kSegments], bottom_vertex[kSegments],
                  top_vertex[kSegments]);
  const EdgeId axis_edge =
      radial_edge(origin0, origin1, origin0_vertex, origin1_vertex);

  const std::uint32_t radial0_bottom_pcurve = editor.add_pcurve(
      Pcurve2d{{uv_bottom(origin0), uv_bottom(bottom[0])}});
  const std::uint32_t radial1_bottom_pcurve = editor.add_pcurve(
      Pcurve2d{{uv_bottom(origin0), uv_bottom(bottom[kSegments])}});
  const std::uint32_t radial0_top_pcurve = editor.add_pcurve(
      Pcurve2d{{uv_top(origin1), uv_top(top[0])}});
  const std::uint32_t radial1_top_pcurve = editor.add_pcurve(
      Pcurve2d{{uv_top(origin1), uv_top(top[kSegments])}});
  const std::uint32_t radial0_side0_pcurve = editor.add_pcurve(
      Pcurve2d{{uv_side0(origin0), uv_side0(bottom[0])}});
  const std::uint32_t radial1_side1_pcurve = editor.add_pcurve(
      Pcurve2d{{uv_side1(origin0), uv_side1(bottom[kSegments])}});
  const std::uint32_t radial0_top_side0_pcurve = editor.add_pcurve(
      Pcurve2d{{uv_side0(origin1), uv_side0(top[0])}});
  const std::uint32_t radial1_top_side1_pcurve = editor.add_pcurve(
      Pcurve2d{{uv_side1(origin1), uv_side1(top[kSegments])}});
  const std::uint32_t side0_line_pcurve = editor.add_pcurve(Pcurve2d{
      {uv_side0(bottom[0]), uv_side0(top[0])}});
  const std::uint32_t side1_line_pcurve = editor.add_pcurve(Pcurve2d{
      {uv_side1(bottom[kSegments]), uv_side1(top[kSegments])}});
  const std::uint32_t side0_line_lateral_pcurve = editor.add_pcurve(
      Pcurve2d{{Point2d{0.0, 0.0}, Point2d{0.0, height}}});
  const std::uint32_t side1_line_lateral_pcurve = editor.add_pcurve(
      Pcurve2d{{Point2d{quarter, 0.0}, Point2d{quarter, height}}});
  const std::uint32_t axis_side0_pcurve = editor.add_pcurve(
      Pcurve2d{{uv_side0(origin0), uv_side0(origin1)}});
  const std::uint32_t axis_side1_pcurve = editor.add_pcurve(
      Pcurve2d{{uv_side1(origin0), uv_side1(origin1)}});

  const FaceId lateral = editor.add_face(lateral_surface, false, tol).value();
  const FaceId bottom_cap =
      editor.add_face(bottom_surface, false, tol).value();
  const FaceId top_cap = editor.add_face(top_surface, false, tol).value();
  const FaceId side0 = editor.add_face(side0_surface, false, tol).value();
  const FaceId side1 = editor.add_face(side1_surface, false, tol).value();

  std::vector<CoedgeSpec> specs;
  const auto push = [&specs](EdgeId edge, bool reversed, std::uint32_t pcurve) {
    CoedgeSpec spec;
    spec.edge = edge;
    spec.reversed = reversed;
    spec.pcurve = pcurve;
    spec.tol = 0.0;
    specs.push_back(spec);
  };

  // Lateral face: bottom arc forward, u = pi/2 line up, top arc backward,
  // u = 0 line down.
  specs.clear();
  for (std::size_t i = 0; i < kSegments; ++i) {
    push(bottom_arc[i], false, bottom_arc_lateral_pcurve[i]);
  }
  push(side1_line, false, side1_line_lateral_pcurve);
  for (std::size_t i = kSegments; i-- > 0;) {
    push(top_arc[i], true, top_arc_lateral_pcurve[i]);
  }
  push(side0_line, true, side0_line_lateral_pcurve);
  editor.add_loop(lateral, specs, true).value();

  // Bottom cap: from the origin to the far radius, back along the arc, home.
  specs.clear();
  push(radial1, false, radial1_bottom_pcurve);
  for (std::size_t i = kSegments; i-- > 0;) {
    push(bottom_arc[i], true, bottom_arc_bottom_pcurve[i]);
  }
  push(radial0, true, radial0_bottom_pcurve);
  editor.add_loop(bottom_cap, specs, true).value();

  // Top cap: mirror of the bottom cap.
  specs.clear();
  push(radial0_top, false, radial0_top_pcurve);
  for (std::size_t i = 0; i < kSegments; ++i) {
    push(top_arc[i], false, top_arc_top_pcurve[i]);
  }
  push(radial1_top, true, radial1_top_pcurve);
  editor.add_loop(top_cap, specs, true).value();

  // Flat side at u = 0.
  specs.clear();
  push(radial0, false, radial0_side0_pcurve);
  push(side0_line, false, side0_line_pcurve);
  push(radial0_top, true, radial0_top_side0_pcurve);
  push(axis_edge, true, axis_side0_pcurve);
  editor.add_loop(side0, specs, true).value();

  // Flat side at u = pi/2.
  specs.clear();
  push(axis_edge, false, axis_side1_pcurve);
  push(radial1_top, false, radial1_top_side1_pcurve);
  push(side1_line, true, side1_line_pcurve);
  push(radial1, true, radial1_side1_pcurve);
  editor.add_loop(side1, specs, true).value();

  const std::array<FaceId, 5> faces{lateral, bottom_cap, top_cap, side0, side1};
  editor.add_shell(faces, true).value();
  return editor.commit().value();
}

// The single lateral Cylinder3 face of the quarter cylinder fixture.
inline FaceId quarter_cylinder_lateral(const Brep& brep) {
  for (const FaceId id : brep.face_ids()) {
    const Face& face = brep.face(id);
    if (face.surface == kNoGeometry) continue;
    if (std::holds_alternative<Cylinder3>(
            brep.geometry().surfaces[face.surface])) {
      return id;
    }
  }
  return FaceId{};
}

}  // namespace geom::test
