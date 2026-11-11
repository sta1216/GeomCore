#include <doctest/doctest.h>

#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <iomanip>
#include <limits>
#include <numbers>
#include <string_view>
#include <variant>
#include <vector>

#include <geom/brep/detail/curved_volume.hpp>
#include <geom/brep/face_uv.hpp>
#include <geom/brep/primitives.hpp>
#include <geom/brep/solid_editor.hpp>
#include <geom/brep/topology.hpp>
#include <geom/core/error.hpp>
#include <geom/curves/analytic_curve.hpp>
#include <geom/surfaces/nurbs_surface.hpp>
#include <geom/surfaces/plane.hpp>

#include "brep_curved_fixture.hpp"

namespace {

constexpr double kPi = std::numbers::pi_v<double>;
constexpr double kTwoPi = 2.0 * kPi;
constexpr double kTolerance = 1e-9;

double relative_error(double actual, double expected) {
  return std::abs(actual - expected) / std::abs(expected);
}

void check_error(const geom::Result<double>& result, geom::Errc code, std::string_view message) {
  REQUIRE_FALSE(result.has_value());
  CHECK(result.error().code == code);
  CHECK(std::string_view(result.error().message) == message);
}

geom::FaceId top_plane_face(const geom::Brep& brep) {
  for (const geom::FaceId id : brep.face_ids()) {
    const geom::Face& face = brep.face(id);
    if (face.surface == geom::kNoGeometry) continue;
    const auto& surface = brep.geometry().surfaces[face.surface];
    if (!std::holds_alternative<geom::Plane3>(surface)) continue;
    const geom::Plane3& plane = std::get<geom::Plane3>(surface);
    if (plane.normal[2] > 0.5) return id;
  }
  return geom::FaceId{};
}

// A one-face, topologically closed shell whose only face is a bilinear NURBS
// patch. curved_solid_volume must reject the face before it attempts to
// integrate it.
geom::Brep make_nurbs_face_solid() {
  geom::SolidEditor editor;
  const std::uint32_t surface = editor.add_surface(geom::NurbsSurface(
      1, 1, {0.0, 0.0, 1.0, 1.0}, {0.0, 0.0, 1.0, 1.0},
      {{0.0, 0.0, 0.0}, {1.0, 0.0, 0.0}, {0.0, 1.0, 0.0}, {1.0, 1.0, 0.0}},
      {1.0, 1.0, 1.0, 1.0}));

  const std::array<geom::Point3d, 3> points{
      geom::Point3d{0.0, 0.0, 0.0},
      geom::Point3d{1.0, 0.0, 0.0},
      geom::Point3d{0.0, 1.0, 0.0},
  };
  std::array<geom::VertexId, 3> vertices{};
  for (std::size_t i = 0; i < vertices.size(); ++i) {
    vertices[i] = editor.add_vertex(points[i], 0.0).value();
  }
  std::array<geom::EdgeId, 3> edges{};
  for (std::size_t i = 0; i < edges.size(); ++i) {
    edges[i] = editor
                   .add_edge(vertices[i], vertices[(i + 1) % vertices.size()],
                             geom::kNoGeometry, 0.0, 1.0, 0.0)
                   .value();
  }

  const geom::FaceId face = editor.add_face(surface, false, 0.0).value();
  std::array<geom::CoedgeSpec, 3> coedges{};
  for (std::size_t i = 0; i < coedges.size(); ++i) {
    coedges[i] = geom::CoedgeSpec{edges[i], false, geom::kNoGeometry, 0.0};
  }
  editor.add_loop(face, coedges, true).value();

  const geom::ShellId shell = editor.add_shell(std::array{face}, true).value();
  editor.add_solid(std::array{shell}, 0.0).value();
  return editor.brep();
}

// make_box() always builds a positively oriented solid and SolidEditor has no
// setter for Face::reversed, so this repeats the make_box construction with
// the top face marked reversed. Topology and geometry are otherwise identical
// to make_box({0,0,0}, {2,2,2}, kTolerance).
geom::Brep make_reversed_top_box() {
  constexpr double side = 2.0;
  const geom::Point3d corners[8] = {
      {0.0, 0.0, 0.0}, {side, 0.0, 0.0}, {side, side, 0.0}, {0.0, side, 0.0},
      {0.0, 0.0, side}, {side, 0.0, side}, {side, side, side}, {0.0, side, side},
  };

  geom::SolidEditor editor;
  geom::VertexId vertices[8];
  for (int i = 0; i < 8; ++i) {
    vertices[i] = editor.add_vertex(corners[i], kTolerance).value();
  }

  struct EdgeSpec {
    int v0;
    int v1;
  };
  const EdgeSpec edge_specs[12] = {
      {0, 1}, {1, 2}, {2, 3}, {3, 0},  // bottom ring
      {4, 5}, {5, 6}, {6, 7}, {7, 4},  // top ring
      {1, 5}, {2, 6}, {3, 7}, {0, 4},  // verticals
  };
  geom::EdgeId edges[12];
  for (int i = 0; i < 12; ++i) {
    const geom::Point3d from = corners[edge_specs[i].v0];
    const geom::Point3d to = corners[edge_specs[i].v1];
    const geom::Vec3d delta = to - from;
    const double length = std::hypot(delta[0], delta[1], delta[2]);
    const geom::Vec3d direction{delta[0] / length, delta[1] / length, delta[2] / length};
    const std::uint32_t curve = editor.add_curve(geom::LineCurve3{from, direction});
    edges[i] = editor
                   .add_edge(vertices[edge_specs[i].v0], vertices[edge_specs[i].v1], curve,
                             0.0, length, kTolerance)
                   .value();
  }

  struct CoedgeUse {
    int edge;
    bool reversed;
  };
  struct FaceSpec {
    geom::Vec3d normal;
    int corner;
    CoedgeUse loop[4];
  };
  // Loop order is counter-clockwise seen from outside along `normal`.
  const FaceSpec face_specs[6] = {
      {{0.0, 0.0, -1.0}, 0, {{3, true}, {2, true}, {1, true}, {0, true}}},    // bottom
      {{0.0, 0.0, 1.0}, 4, {{4, false}, {5, false}, {6, false}, {7, false}}},  // top
      {{0.0, -1.0, 0.0}, 0, {{0, false}, {8, false}, {4, true}, {11, true}}},  // front
      {{0.0, 1.0, 0.0}, 2, {{2, false}, {10, false}, {6, true}, {9, true}}},   // back
      {{-1.0, 0.0, 0.0}, 0, {{11, false}, {7, true}, {10, true}, {3, false}}}, // left
      {{1.0, 0.0, 0.0}, 1, {{1, false}, {9, false}, {5, true}, {8, true}}},    // right
  };

  std::array<geom::FaceId, 6> faces{};
  for (int f = 0; f < 6; ++f) {
    const FaceSpec& spec = face_specs[f];
    const geom::Point3d corner = corners[spec.corner];
    const std::uint32_t surface = editor.add_surface(geom::Plane3{
        spec.normal,
        geom::dot(spec.normal, geom::Vec3d{corner[0], corner[1], corner[2]})});
    // Only the top face is marked reversed.
    faces[static_cast<std::size_t>(f)] =
        editor.add_face(surface, f == 1, kTolerance).value();

    std::array<geom::CoedgeSpec, 4> coedges{};
    for (int i = 0; i < 4; ++i) {
      coedges[static_cast<std::size_t>(i)] =
          geom::CoedgeSpec{edges[spec.loop[i].edge], spec.loop[i].reversed};
    }
    editor.add_loop(faces[static_cast<std::size_t>(f)], coedges, true).value();
  }

  const geom::ShellId shell = editor.add_shell(faces, true).value();
  editor.add_solid(std::array{shell}, kTolerance).value();
  return editor.commit().value();
}

// A closed half cylinder: a lateral Cylinder3 patch over u in [0, pi], two
// semicircular caps whose rims carry exact CircleCurve3 provenance (sixteen arc
// segments each) and the flat rectangle on the diameter plane. The cap rims are
// true arcs rather than chords, so curved_volume_planar_area must integrate
// pi*r^2/2 per cap and the solid volume must be pi*r^2*h/2.
geom::Brep make_half_cylinder(double radius, double height, double tol) {
  constexpr std::size_t kSegments = 16;
  const double delta = kPi / static_cast<double>(kSegments);

  geom::SolidEditor editor;
  const std::uint32_t lateral_surface = editor.add_surface(geom::Cylinder3{
      geom::Point3d{0.0, 0.0, 0.0}, geom::Vec3d{0.0, 0.0, 1.0},
      geom::Vec3d{1.0, 0.0, 0.0}, radius});
  const std::uint32_t bottom_surface =
      editor.add_surface(geom::Plane3{geom::Vec3d{0.0, 0.0, -1.0}, 0.0});
  const std::uint32_t top_surface =
      editor.add_surface(geom::Plane3{geom::Vec3d{0.0, 0.0, 1.0}, height});
  const std::uint32_t flat_surface =
      editor.add_surface(geom::Plane3{geom::Vec3d{0.0, -1.0, 0.0}, 0.0});

  std::array<geom::Point3d, kSegments + 1> bottom{};
  std::array<geom::Point3d, kSegments + 1> top{};
  std::array<double, kSegments + 1> angle{};
  for (std::size_t i = 0; i <= kSegments; ++i) {
    angle[i] = delta * static_cast<double>(i);
    const double c = std::cos(angle[i]);
    const double s = std::sin(angle[i]);
    bottom[i] = geom::Point3d{radius * c, radius * s, 0.0};
    top[i] = geom::Point3d{radius * c, radius * s, height};
  }

  std::array<geom::VertexId, kSegments + 1> bottom_vertex{};
  std::array<geom::VertexId, kSegments + 1> top_vertex{};
  for (std::size_t i = 0; i <= kSegments; ++i) {
    bottom_vertex[i] = editor.add_vertex(bottom[i], 0.0).value();
    top_vertex[i] = editor.add_vertex(top[i], 0.0).value();
  }

  // Parameter maps of the planar faces, in the frames make_plane_frame() derives
  // from their plane normals.
  const auto uv_bottom = [](const geom::Point3d& p) {
    return geom::Point2d{p[1], p[0]};
  };
  const auto uv_top = [](const geom::Point3d& p) {
    return geom::Point2d{-p[1], p[0]};
  };
  const auto uv_flat = [](const geom::Point3d& p) {
    return geom::Point2d{-p[2], p[0]};
  };

  const geom::CircleCurve3 bottom_circle{geom::Point3d{0.0, 0.0, 0.0},
                                         geom::Vec3d{1.0, 0.0, 0.0},
                                         geom::Vec3d{0.0, 1.0, 0.0}, radius};
  const geom::CircleCurve3 top_circle{geom::Point3d{0.0, 0.0, height},
                                      geom::Vec3d{1.0, 0.0, 0.0},
                                      geom::Vec3d{0.0, 1.0, 0.0}, radius};

  std::array<geom::EdgeId, kSegments> bottom_arc{};
  std::array<geom::EdgeId, kSegments> top_arc{};
  std::array<std::uint32_t, kSegments> bottom_arc_bottom_pcurve{};
  std::array<std::uint32_t, kSegments> bottom_arc_lateral_pcurve{};
  std::array<std::uint32_t, kSegments> top_arc_top_pcurve{};
  std::array<std::uint32_t, kSegments> top_arc_lateral_pcurve{};
  for (std::size_t i = 0; i < kSegments; ++i) {
    bottom_arc[i] =
        editor
            .add_edge(bottom_vertex[i], bottom_vertex[i + 1],
                      editor.add_curve(bottom_circle), angle[i], angle[i + 1], 0.0)
            .value();
    top_arc[i] =
        editor
            .add_edge(top_vertex[i], top_vertex[i + 1], editor.add_curve(top_circle),
                      angle[i], angle[i + 1], 0.0)
            .value();
    bottom_arc_bottom_pcurve[i] = editor.add_pcurve(
        geom::Pcurve2d{{uv_bottom(bottom[i]), uv_bottom(bottom[i + 1])}});
    bottom_arc_lateral_pcurve[i] = editor.add_pcurve(geom::Pcurve2d{
        {geom::Point2d{angle[i], 0.0}, geom::Point2d{angle[i + 1], 0.0}}});
    top_arc_top_pcurve[i] = editor.add_pcurve(
        geom::Pcurve2d{{uv_top(top[i]), uv_top(top[i + 1])}});
    top_arc_lateral_pcurve[i] = editor.add_pcurve(geom::Pcurve2d{
        {geom::Point2d{angle[i], height}, geom::Point2d{angle[i + 1], height}}});
  }

  const auto line_edge = [&editor](const geom::Point3d& from, const geom::Point3d& to,
                                   geom::VertexId v0, geom::VertexId v1) {
    const std::uint32_t curve = editor.add_curve(geom::LineCurve3{
        from, geom::Vec3d{to[0] - from[0], to[1] - from[1], to[2] - from[2]}});
    return editor.add_edge(v0, v1, curve, 0.0, 1.0, 0.0).value();
  };

  const geom::EdgeId bottom_diameter = line_edge(bottom[kSegments], bottom[0],
                                                 bottom_vertex[kSegments],
                                                 bottom_vertex[0]);
  const geom::EdgeId top_diameter =
      line_edge(top[kSegments], top[0], top_vertex[kSegments], top_vertex[0]);
  const geom::EdgeId side0_line = line_edge(bottom[0], top[0], bottom_vertex[0],
                                            top_vertex[0]);
  const geom::EdgeId side1_line =
      line_edge(bottom[kSegments], top[kSegments], bottom_vertex[kSegments],
                top_vertex[kSegments]);

  const std::uint32_t bottom_diameter_bottom_pcurve = editor.add_pcurve(
      geom::Pcurve2d{{uv_bottom(bottom[kSegments]), uv_bottom(bottom[0])}});
  const std::uint32_t bottom_diameter_flat_pcurve = editor.add_pcurve(
      geom::Pcurve2d{{uv_flat(bottom[kSegments]), uv_flat(bottom[0])}});
  const std::uint32_t top_diameter_top_pcurve = editor.add_pcurve(
      geom::Pcurve2d{{uv_top(top[kSegments]), uv_top(top[0])}});
  const std::uint32_t top_diameter_flat_pcurve = editor.add_pcurve(
      geom::Pcurve2d{{uv_flat(top[kSegments]), uv_flat(top[0])}});
  const std::uint32_t side0_line_lateral_pcurve = editor.add_pcurve(
      geom::Pcurve2d{{geom::Point2d{0.0, 0.0}, geom::Point2d{0.0, height}}});
  const std::uint32_t side1_line_lateral_pcurve = editor.add_pcurve(
      geom::Pcurve2d{{geom::Point2d{kPi, 0.0}, geom::Point2d{kPi, height}}});
  const std::uint32_t side0_line_flat_pcurve = editor.add_pcurve(
      geom::Pcurve2d{{uv_flat(bottom[0]), uv_flat(top[0])}});
  const std::uint32_t side1_line_flat_pcurve = editor.add_pcurve(
      geom::Pcurve2d{{uv_flat(bottom[kSegments]), uv_flat(top[kSegments])}});

  const geom::FaceId lateral = editor.add_face(lateral_surface, false, tol).value();
  const geom::FaceId bottom_cap = editor.add_face(bottom_surface, false, tol).value();
  const geom::FaceId top_cap = editor.add_face(top_surface, false, tol).value();
  const geom::FaceId flat = editor.add_face(flat_surface, false, tol).value();

  std::vector<geom::CoedgeSpec> specs;
  const auto push = [&specs](geom::EdgeId edge, bool reversed, std::uint32_t pcurve) {
    specs.push_back(geom::CoedgeSpec{edge, reversed, pcurve, 0.0});
  };

  // Lateral face: the bottom half rim forward, the u = pi line up, the top half
  // rim backward, the u = 0 line down (counter-clockwise in (u, v)).
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

  // Bottom cap (outward normal -z): diameter first, then the rim backward. This
  // is clockwise in world xy, hence counter-clockwise in the (y, x) frame.
  specs.clear();
  push(bottom_diameter, true, bottom_diameter_bottom_pcurve);
  for (std::size_t i = kSegments; i-- > 0;) {
    push(bottom_arc[i], true, bottom_arc_bottom_pcurve[i]);
  }
  editor.add_loop(bottom_cap, specs, true).value();

  // Top cap (outward normal +z): rim forward, then the diameter home.
  specs.clear();
  for (std::size_t i = 0; i < kSegments; ++i) {
    push(top_arc[i], false, top_arc_top_pcurve[i]);
  }
  push(top_diameter, false, top_diameter_top_pcurve);
  editor.add_loop(top_cap, specs, true).value();

  // Flat face on the diameter plane (outward normal -y).
  specs.clear();
  push(side0_line, false, side0_line_flat_pcurve);
  push(top_diameter, true, top_diameter_flat_pcurve);
  push(side1_line, true, side1_line_flat_pcurve);
  push(bottom_diameter, false, bottom_diameter_flat_pcurve);
  editor.add_loop(flat, specs, true).value();

  const std::array<geom::FaceId, 4> faces{lateral, bottom_cap, top_cap, flat};
  const geom::ShellId shell = editor.add_shell(faces, true).value();
  editor.add_solid(std::array{shell}, tol).value();
  return editor.commit().value();
}

}  // namespace

TEST_CASE("curved_solid_volume: box") {
  const geom::Brep brep = geom::make_box({0.0, 0.0, 0.0}, {2.0, 2.0, 2.0}, kTolerance).value();
  const geom::Result<double> volume =
      geom::detail::curved_solid_volume(brep, brep.solid_ids().front(), kTolerance);
  REQUIRE(volume.has_value());
  CHECK(relative_error(volume.value(), 8.0) <= 1e-12);
}

TEST_CASE("curved_solid_volume: cylinder") {
  constexpr double radius = 2.0;
  constexpr double height = 3.0;
  const geom::Brep brep =
      geom::make_cylinder({0.0, 0.0, 0.0}, {0.0, 0.0, 1.0}, radius, height, kTolerance).value();
  const geom::FaceId cap = top_plane_face(brep);
  REQUIRE_FALSE(cap.is_null());

  const geom::Result<std::vector<geom::FaceUvLoop>> loops =
      geom::face_loops_uv(brep, cap, kTolerance);
  REQUIRE(loops.has_value());
  REQUIRE(loops.value().size() == 1);
  REQUIRE(loops.value().front().outer);
  const std::size_t segments = loops.value().front().points.size();
  CHECK(segments == 32);

  const geom::Result<double> volume =
      geom::detail::curved_solid_volume(brep, brep.solid_ids().front(), kTolerance);
  REQUIRE(volume.has_value());

  // The caps carry exact circle provenance, so the arc-aware planar area uses
  // pi*r^2 rather than the 32-gon; the lateral face is integrated on the exact
  // cylinder. The result must therefore equal the analytic volume.
  CHECK(relative_error(volume.value(), kPi * radius * radius * height) <= 1e-9);
}

TEST_CASE("curved_solid_volume: half cylinder keeps exact semicircle rims") {
  constexpr double radius = 2.0;
  constexpr double height = 3.0;
  const geom::Brep brep = make_half_cylinder(radius, height, kTolerance);

  geom::FaceId bottom_cap;
  geom::FaceId top_cap;
  for (const geom::FaceId id : brep.face_ids()) {
    const geom::Face& face = brep.face(id);
    if (face.surface == geom::kNoGeometry) continue;
    const auto& surface = brep.geometry().surfaces[face.surface];
    if (!std::holds_alternative<geom::Plane3>(surface)) continue;
    const geom::Plane3& plane = std::get<geom::Plane3>(surface);
    if (plane.normal[2] < -0.5) bottom_cap = id;
    if (plane.normal[2] > 0.5) top_cap = id;
  }
  REQUIRE_FALSE(bottom_cap.is_null());
  REQUIRE_FALSE(top_cap.is_null());

  const geom::Result<std::vector<geom::FaceUvLoop>> loops =
      geom::face_loops_uv(brep, bottom_cap, kTolerance);
  REQUIRE(loops.has_value());
  REQUIRE(loops.value().size() == 1);
  const geom::FaceUvLoop& rim = loops.value().front();
  CHECK(rim.outer);
  // Sixteen arc segments plus the straight diameter, which carries no arc
  // reference even though both of its endpoints lie on the rim circle.
  REQUIRE(rim.points.size() == 17);
  REQUIRE(rim.arcs.size() == rim.points.size());
  // The loop starts with the straight diameter, then follows the sixteen rim
  // segments.
  CHECK_FALSE(rim.arcs.front().has_value());
  for (std::size_t i = 1; i < rim.points.size(); ++i) {
    REQUIRE(rim.arcs[i].has_value());
    CHECK(std::fabs(std::fabs(rim.arcs[i]->t1 - rim.arcs[i]->t0) - kPi / 16.0) <=
          1e-12);
  }

  const geom::Result<geom::detail::FaceUvProjection> projection =
      geom::detail::face_uv_projection(brep, brep.face(bottom_cap), kTolerance);
  REQUIRE(projection.has_value());
  const geom::Result<double> area = geom::detail::curved_volume_planar_area(
      loops.value(), projection.value().frame, kTolerance);
  REQUIRE(area.has_value());
  CHECK(relative_error(area.value(), 0.5 * kPi * radius * radius) <= 1e-12);

  const geom::Result<double> volume =
      geom::detail::curved_solid_volume(brep, brep.solid_ids().front(), kTolerance);
  REQUIRE(volume.has_value());
  CHECK(relative_error(volume.value(), 0.5 * kPi * radius * radius * height) <= 1e-9);
}

TEST_CASE("curved_solid_volume: a concave arc keeps its negative contribution") {
  // A 2 x 2 square with a half-disk notch cut out of its bottom edge. The notch
  // arc is travelled clockwise about its own circle (sweep -pi/16 per segment),
  // so the exact area is the square minus the half disk, not plus it.
  constexpr double radius = 1.0;
  const geom::Result<geom::PlaneFrame> frame =
      geom::make_plane_frame(geom::Plane3{geom::Vec3d{0.0, 0.0, 1.0}, 0.0});
  REQUIRE(frame.has_value());

  const geom::CircleCurve3 circle{geom::Point3d{1.0, 0.0, 0.0},
                                  geom::Vec3d{1.0, 0.0, 0.0},
                                  geom::Vec3d{0.0, 1.0, 0.0}, radius};
  constexpr std::size_t kSegments = 16;
  const double delta = kPi / static_cast<double>(kSegments);

  geom::FaceUvLoop loop;
  loop.outer = true;
  // The plane frame for normal +z maps world (x, y) to uv (-y, x), so the notch
  // arc is a circle of radius 1 centred at uv (0, 1) swept from -pi/2 down to
  // -3*pi/2 (clockwise about its circle).
  const geom::Point2d center{0.0, 1.0};
  for (std::size_t i = 0; i <= kSegments; ++i) {
    const double angle = -0.5 * kPi - delta * static_cast<double>(i);
    loop.points.push_back(geom::Point2d{center[0] + radius * std::cos(angle),
                                        center[1] + radius * std::sin(angle)});
    if (i < kSegments) {
      loop.arcs.push_back(geom::FaceUvArcRef{circle, angle, angle - delta});
    }
  }
  loop.points.push_back(geom::Point2d{-2.0, 2.0});
  loop.points.push_back(geom::Point2d{-2.0, 0.0});
  loop.arcs.push_back(std::nullopt);  // (0, 2) -> (-2, 2)
  loop.arcs.push_back(std::nullopt);  // (-2, 2) -> (-2, 0)
  loop.arcs.push_back(std::nullopt);  // (-2, 0) -> (0, 0), closing the notch
  REQUIRE(loop.arcs.size() == loop.points.size());

  const geom::Result<double> area =
      geom::detail::curved_volume_planar_area({loop}, frame.value(), kTolerance);
  REQUIRE(area.has_value());
  CHECK(relative_error(area.value(), 4.0 - 0.5 * kPi * radius * radius) <= 1e-12);
}

TEST_CASE("curved_solid_volume: an arc wider than pi is sampled per segment") {
  // A 270-degree pie of radius 1: the rim is split into 24 segments of 3*pi/48,
  // each well below pi, so the signed per-segment sweeps integrate to 3*pi/4.
  constexpr std::size_t kSegments = 24;
  constexpr double radius = 1.0;
  const double sweep = 1.5 * kPi / static_cast<double>(kSegments);
  const geom::Result<geom::PlaneFrame> frame =
      geom::make_plane_frame(geom::Plane3{geom::Vec3d{0.0, 0.0, 1.0}, 0.0});
  REQUIRE(frame.has_value());

  const geom::CircleCurve3 circle{geom::Point3d{0.0, 0.0, 0.0},
                                  geom::Vec3d{1.0, 0.0, 0.0},
                                  geom::Vec3d{0.0, 1.0, 0.0}, radius};
  geom::FaceUvLoop loop;
  loop.outer = true;
  for (std::size_t i = 0; i <= kSegments; ++i) {
    const double angle = sweep * static_cast<double>(i);
    loop.points.push_back(
        geom::Point2d{radius * std::cos(angle), radius * std::sin(angle)});
    if (i < kSegments) {
      loop.arcs.push_back(geom::FaceUvArcRef{circle, angle, angle + sweep});
    }
  }
  loop.points.push_back(geom::Point2d{0.0, 0.0});
  loop.arcs.push_back(std::nullopt);  // last rim point -> center
  loop.arcs.push_back(std::nullopt);  // center -> first rim point
  REQUIRE(loop.arcs.size() == loop.points.size());

  const geom::Result<double> area =
      geom::detail::curved_volume_planar_area({loop}, frame.value(), kTolerance);
  REQUIRE(area.has_value());
  CHECK(relative_error(area.value(), 0.75 * kPi * radius * radius) <= 1e-12);
}

TEST_CASE("curved_solid_volume: quarter cylinder") {
  constexpr double radius = 1.0;
  constexpr double height = 1.0;
  // The fixture commits a closed shell but no solid (and must not be edited),
  // so attach the shell to a solid here.
  geom::SolidEditor editor(geom::test::make_quarter_cylinder(radius, height, kTolerance));
  REQUIRE(editor.brep().shell_count() == 1);
  const geom::ShellId shell = editor.brep().shell_ids().front();
  editor.add_solid(std::array{shell}, kTolerance).value();
  const geom::Brep brep = editor.commit().value();
  const geom::FaceId lateral = geom::test::quarter_cylinder_lateral(brep);
  REQUIRE_FALSE(lateral.is_null());

  const geom::Result<double> volume =
      geom::detail::curved_solid_volume(brep, brep.solid_ids().front(), kTolerance);
  REQUIRE(volume.has_value());

  // The fixture has 16 chords per quarter arc. The lateral quarter stays a
  // smooth cylinder patch, so the result is the half-cylinder lateral term
  // plus the inscribed 16-gon top-cap term.
  constexpr double n = 16.0;
  const double cap_area = 0.5 * n * radius * radius * std::sin((kPi / 2.0) / n);
  const double expected = height * (kPi * radius * radius / 2.0 + cap_area) / 3.0;
  CHECK(relative_error(volume.value(), expected) <= 1e-9);
  CHECK(relative_error(volume.value(), kPi * radius * radius * height / 4.0) <= 2e-3);
}

TEST_CASE("curved_solid_volume: translated cylinder") {
  constexpr double radius = 2.0;
  constexpr double height = 3.0;
  const geom::Brep origin_brep =
      geom::make_cylinder({0.0, 0.0, 0.0}, {0.0, 0.0, 1.0}, radius, height, kTolerance).value();
  const geom::Brep translated_brep =
      geom::make_cylinder({1000.0, 0.0, 0.0}, {0.0, 0.0, 1.0}, radius, height, kTolerance).value();

  const geom::Result<double> origin_volume =
      geom::detail::curved_solid_volume(origin_brep, origin_brep.solid_ids().front(), kTolerance);
  const geom::Result<double> translated_volume = geom::detail::curved_solid_volume(
      translated_brep, translated_brep.solid_ids().front(), kTolerance);
  REQUIRE(origin_volume.has_value());
  REQUIRE(translated_volume.has_value());

  MESSAGE("origin cylinder volume: " << std::setprecision(17) << origin_volume.value());
  MESSAGE("translated cylinder volume: " << std::setprecision(17) << translated_volume.value());

  // Quadrature over the raw uv triangles used to depend on where the solid sat.
  CHECK(relative_error(translated_volume.value(), origin_volume.value()) <= 1e-9);
  CHECK(relative_error(translated_volume.value(), kPi * radius * radius * height) <= 1e-9);
}

TEST_CASE("curved_solid_volume: box far from the origin") {
  constexpr double base = 1e9;
  const geom::Brep brep =
      geom::make_box({base, base, base}, {base + 2.0, base + 2.0, base + 2.0}, kTolerance)
          .value();
  const geom::Result<double> volume =
      geom::detail::curved_solid_volume(brep, brep.solid_ids().front(), kTolerance);
  REQUIRE(volume.has_value());
  MESSAGE("box at 1e9 volume: " << std::setprecision(17) << volume.value());

  // The absolute plane offsets are 1e9; the uv polygon shift plus the reference
  // point of the divergence integrand keep the two large opposite face terms
  // from cancelling the 8 away.
  CHECK(relative_error(volume.value(), 8.0) <= 1e-9);
}

TEST_CASE("curved_solid_volume: torus") {
  constexpr double major_radius = 2.0;
  constexpr double minor_radius = 0.5;
  const geom::Brep brep =
      geom::make_torus({0.0, 0.0, 0.0}, {0.0, 0.0, 1.0}, major_radius, minor_radius, kTolerance)
          .value();
  const geom::Result<double> volume =
      geom::detail::curved_solid_volume(brep, brep.solid_ids().front(), kTolerance);
  REQUIRE(volume.has_value());
  MESSAGE("torus volume: " << std::setprecision(17) << volume.value());

  const double expected =
      2.0 * kPi * kPi * major_radius * minor_radius * minor_radius;
  CHECK(relative_error(volume.value(), expected) <= 1e-8);
}

TEST_CASE("curved_solid_volume: reversed top face") {
  const geom::Brep brep = make_reversed_top_box();
  const geom::Result<double> volume =
      geom::detail::curved_solid_volume(brep, brep.solid_ids().front(), kTolerance);
  REQUIRE(volume.has_value());
  MESSAGE("reversed top box volume: " << std::setprecision(17) << volume.value());

  // The six outward faces give 8; reversing the top face (which contributes
  // 2 * 4 / 3 = 8/3) subtracts it twice.
  CHECK(relative_error(volume.value(), 8.0 / 3.0) <= 1e-9);
}

TEST_CASE("curved_solid_volume: error paths") {
  const geom::Brep brep = geom::make_box({0.0, 0.0, 0.0}, {2.0, 2.0, 2.0}, kTolerance).value();
  const geom::SolidId solid = brep.solid_ids().front();

  check_error(geom::detail::curved_solid_volume(brep, geom::SolidId{}, kTolerance),
              geom::Errc::index_out_of_range,
              "curved_solid_volume: solid is not alive");

  const geom::Brep empty;
  check_error(geom::detail::curved_solid_volume(empty, geom::SolidId::make(0, 1), kTolerance),
              geom::Errc::index_out_of_range,
              "curved_solid_volume: solid is not alive");

  check_error(geom::detail::curved_solid_volume(brep, solid, 0.0),
              geom::Errc::tolerance_violation,
              "curved_solid_volume: tolerance must be finite and positive");
  check_error(geom::detail::curved_solid_volume(brep, solid,
                                                std::numeric_limits<double>::quiet_NaN()),
              geom::Errc::tolerance_violation,
              "curved_solid_volume: tolerance must be finite and positive");
}

TEST_CASE("curved_solid_volume: rejects NURBS faces") {
  const geom::Brep brep = make_nurbs_face_solid();
  REQUIRE(brep.solid_count() == 1);
  check_error(geom::detail::curved_solid_volume(brep, brep.solid_ids().front(), kTolerance),
              geom::Errc::unsupported,
              "curved_solid_volume: NURBS surfaces are not supported");
}
