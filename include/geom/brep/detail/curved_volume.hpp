#pragma once

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <span>
#include <utility>
#include <variant>
#include <vector>

#include <geom/brep/detail/ssi_surface.hpp>
#include <geom/brep/face_uv.hpp>
#include <geom/brep/topology.hpp>
#include <geom/cgeom/polygon2.hpp>
#include <geom/cgeom/triangulate2.hpp>
#include <geom/core/error.hpp>
#include <geom/core/result.hpp>
#include <geom/primitives/point.hpp>
#include <geom/primitives/vec.hpp>
#include <geom/surfaces/nurbs_surface.hpp>
#include <geom/surfaces/plane.hpp>

namespace geom {
inline namespace v1 {
namespace detail {

// Seven-point degree-five Dunavant rule on the reference triangle. The
// barycentric coordinates are (a, b, c) and the weights sum to one, so scaling
// by a triangle's area integrates over that triangle. The weights below are
// the full-precision Dunavant values; the eight-decimal truncations quoted in
// the task brief sum to 0.99999999 and would shrink every volume by 1e-8.
struct CurvedVolumeQuadraturePoint {
  double a;
  double b;
  double c;
  double weight;
};

inline constexpr std::array<CurvedVolumeQuadraturePoint, 7> kCurvedVolumeQuadraturePoints{{
    {1.0 / 3.0, 1.0 / 3.0, 1.0 / 3.0, 0.225},
    {0.059715871789770, 0.470142064105115, 0.470142064105115, 0.132394152788506},
    {0.470142064105115, 0.059715871789770, 0.470142064105115, 0.132394152788506},
    {0.470142064105115, 0.470142064105115, 0.059715871789770, 0.132394152788506},
    {0.797426985353087, 0.101286507323456, 0.101286507323456, 0.125939180544827},
    {0.101286507323456, 0.797426985353087, 0.101286507323456, 0.125939180544827},
    {0.101286507323456, 0.101286507323456, 0.797426985353087, 0.125939180544827},
}};

// A curved uv triangle is split into 2^level pieces per edge so that no piece
// spans more than roughly kCurvedVolumeTargetEdge of parameter space. The cap
// of eight levels (256 pieces per edge) keeps the cost bounded while still
// resolving a full 2*pi turn to about 0.025, which converges the cylinder and
// torus integrals to roughly 1e-9 relative.
inline constexpr double kCurvedVolumeTargetEdge = 0.025;

inline int curved_volume_subdivision_level(double max_edge) noexcept {
  constexpr double kMaxPanels = 256.0;
  if (!(max_edge > kCurvedVolumeTargetEdge)) return 0;
  const double panels = max_edge / kCurvedVolumeTargetEdge;
  if (!(panels < kMaxPanels)) return 8;
  return std::min(8, static_cast<int>(std::ceil(std::log2(panels))));
}

// Signed volume of a closed solid obtained by integrating
// (1/3) S . (S_u x S_v) over every face's uv polygon.
//
// Scope: planar faces are integrated exactly by the seven-point rule. Cylinder
// lateral faces and caps, and trimmed analytic patches in general, are
// integrated on the exact surface with the smooth quadrature refined until it
// converges. Sphere and cone seam and apex cases are not handled, and NURBS
// surfaces are rejected. The result is the volume of the polygonal (chord)
// B-rep that the sampled face boundaries describe, not the analytic ideal: the
// curved faces are approximated by the chords and parameter-domain polygons
// that face_loops_uv produces, so a cylinder with n cap chords comes out
// slightly smaller than pi * r^2 * h.
//
// Conditioning: every uv polygon is translated by its first vertex before
// triangulation and before the uv-area computation, and each sample is
// translated back before surface_eval(), so the parameterization is unchanged
// while the shoelace sums stay well conditioned even when the face parameters
// are of order 1e9. The divergence integrand is additionally measured from a
// fixed point of the solid: a closed shell has zero normal flux, so this is an
// exact change of origin, but it keeps a solid placed far from the world origin
// from losing its volume to cancellation between the large face contributions.
// Signed area of one face loop using Green's theorem, with exact circular-arc
// segments when the provenance is available. `frame` is the face's parameter
// frame; `shift` keeps the sums well conditioned for large parameters.
inline double curved_volume_arc_term(const Point2d& a, const Point2d& b,
                                     const Point2d& center, double radius, double turn) {
  const double ax = a[0] - center[0];
  const double ay = a[1] - center[1];
  const double bx = b[0] - center[0];
  const double by = b[1] - center[1];
  const double angle_a = std::atan2(ay, ax);
  const double angle_b = std::atan2(by, bx);
  double delta = angle_b - angle_a;
  if (turn > 0.0) {
    while (delta <= 0.0) delta += 2.0 * kFaceUvPi;
  } else {
    while (delta >= 0.0) delta -= 2.0 * kFaceUvPi;
  }
  return 0.5 * (radius * radius * delta +
                radius * (center[0] * (std::sin(angle_b) - std::sin(angle_a)) -
                          center[1] * (std::cos(angle_b) - std::cos(angle_a))));
}

inline bool curved_volume_same_circle(const FaceUvArcRef& a, const FaceUvArcRef& b,
                                      double tolerance) {
  const double slack = std::max(tolerance, 1e-12 * std::max(a.circle.radius, 1.0));
  const double center_distance = std::hypot(
      a.circle.center[0] - b.circle.center[0], a.circle.center[1] - b.circle.center[1],
      a.circle.center[2] - b.circle.center[2]);
  return center_distance <= slack && std::fabs(a.circle.radius - b.circle.radius) <= slack;
}

// Exact signed area of a planar face's loops in parameter space. Circle arcs
// whose plane is parallel to the face plane contribute their circular-segment
// area; everything else is a straight chord. This is what makes a cylinder cap
// integrate to pi*r^2 rather than to the 32-gon polygon.
inline Result<double> curved_volume_planar_area(const std::vector<FaceUvLoop>& loops,
                                                const PlaneFrame& frame, double tolerance) {
  const Point2d shift = loops.front().points.front();
  double area = 0.0;
  for (const FaceUvLoop& loop : loops) {
    const std::size_t count = loop.points.size();
    double twice_area = 0.0;
    for (std::size_t i = 0; i < count; ++i) {
      const Point2d& a = loop.points[i];
      const Point2d& b = loop.points[(i + 1) % count];
      twice_area += (a[0] - shift[0]) * (b[1] - shift[1]) -
                    (b[0] - shift[0]) * (a[1] - shift[1]);
    }
    const double turn = twice_area >= 0.0 ? 1.0 : -1.0;
    for (std::size_t i = 0; i < count; ++i) {
      const std::size_t next = (i + 1) % count;
      const Point2d& a = loop.points[i];
      const Point2d& b = loop.points[next];
      bool arc = false;
      double radius = 0.0;
      Point2d center{0.0, 0.0};
      if (loop.arcs.size() == count && loop.arcs[i].has_value() &&
          loop.arcs[next].has_value() &&
          curved_volume_same_circle(*loop.arcs[i], *loop.arcs[next], tolerance)) {
        const CircleCurve3& circle = loop.arcs[i]->circle;
        const Vec3d circle_normal = cross(circle.x_axis, circle.y_axis);
        const double circle_normal_length =
            std::hypot(circle_normal[0], circle_normal[1], circle_normal[2]);
        if (circle_normal_length > 0.0) {
          const double alignment =
              std::fabs(dot(circle_normal, frame.normal)) / circle_normal_length;
          const double plane_distance =
              std::fabs(dot(circle.center - frame.origin, frame.normal));
          const double slack =
              std::max(tolerance, 1e-12 * std::max(circle.radius, 1.0));
          if (alignment >= 1.0 - 1e-12 && plane_distance <= slack) {
            const Vec3d offset = circle.center - frame.origin;
            center[0] = dot(offset, frame.u) - shift[0];
            center[1] = dot(offset, frame.v) - shift[1];
            radius = circle.radius;
            arc = true;
          }
        }
      }
      if (arc) {
        area += curved_volume_arc_term(Point2d{a[0] - shift[0], a[1] - shift[1]},
                                       Point2d{b[0] - shift[0], b[1] - shift[1]}, center,
                                       radius, turn);
      } else {
        area += 0.5 * ((a[0] - shift[0]) * (b[1] - shift[1]) -
                       (b[0] - shift[0]) * (a[1] - shift[1]));
      }
    }
  }
  return area;
}

inline Result<double> curved_solid_volume(const Brep& brep, SolidId solid, double tolerance) {
  if (!std::isfinite(tolerance) || tolerance <= 0.0) {
    return Error{Errc::tolerance_violation,
                 "curved_solid_volume: tolerance must be finite and positive", 0u, 0u};
  }
  if (solid.is_null() || !brep.solid_alive(solid)) {
    return Error{Errc::index_out_of_range, "curved_solid_volume: solid is not alive",
                 solid.index(), 0u};
  }

  const std::span<const ShellId> shells = brep.solid_shells(solid);
  if (shells.empty()) {
    return Error{Errc::invalid_topology, "curved_solid_volume: solid has no shells",
                 solid.index(), 0u};
  }

  // Validate the whole solid before doing any integration so an unsupported
  // face or an open shell is reported regardless of face order.
  for (const ShellId shell_id : shells) {
    const Shell& shell = brep.shell(shell_id);
    if (!shell.closed) {
      return Error{Errc::invalid_topology, "curved_solid_volume: shell is not closed",
                   shell_id.index(), solid.index()};
    }
    const std::span<const FaceId> faces = brep.shell_faces(shell_id);
    if (faces.empty()) {
      return Error{Errc::invalid_topology, "curved_solid_volume: shell has no faces",
                   shell_id.index(), solid.index()};
    }
    for (const FaceId face_id : faces) {
      const Face& face = brep.face(face_id);
      if (face.surface == kNoGeometry ||
          static_cast<std::size_t>(face.surface) >= brep.geometry().surfaces.size()) {
        return Error{Errc::index_out_of_range,
                     "curved_solid_volume: face surface is not valid", face_id.index(),
                     solid.index()};
      }
      if (std::holds_alternative<NurbsSurface>(brep.geometry().surfaces[face.surface])) {
        return Error{Errc::unsupported,
                     "curved_solid_volume: NURBS surfaces are not supported",
                     face_id.index(), solid.index()};
      }
    }
  }

  double total = 0.0;
  bool has_reference = false;
  Point3d reference{0.0, 0.0, 0.0};
  for (const ShellId shell_id : shells) {
    for (const FaceId face_id : brep.shell_faces(shell_id)) {
      const Face& face = brep.face(face_id);

      const Result<std::vector<FaceUvLoop>> loops_result =
          face_loops_uv(brep, face_id, tolerance);
      if (!loops_result.has_value()) return loops_result.error();
      const std::vector<FaceUvLoop>& loops = loops_result.value();

      std::size_t outer_index = loops.size();
      for (std::size_t i = 0; i < loops.size(); ++i) {
        if (!loops[i].outer) continue;
        if (outer_index != loops.size()) {
          return Error{Errc::invalid_topology,
                       "curved_solid_volume: face must have exactly one outer loop",
                       face_id.index(), 0u};
        }
        outer_index = i;
      }
      if (outer_index == loops.size()) {
        return Error{Errc::invalid_topology,
                     "curved_solid_volume: face must have exactly one outer loop",
                     face_id.index(), 0u};
      }

      const SurfaceGeometry& surface = brep.geometry().surfaces[face.surface];
      const bool planar = std::holds_alternative<Plane3>(surface);
      if (planar) {
        // Planar faces have a constant divergence integrand, so only the exact
        // (arc-aware) area is needed; skip the chord triangulation entirely.
        const Result<FaceUvProjection> projection =
            face_uv_projection(brep, face, tolerance);
        if (!projection.has_value()) return projection.error();
        const Result<double> area =
            curved_volume_planar_area(loops, projection.value().frame, tolerance);
        if (!area.has_value()) return area.error();
        const Point2d& first = loops[outer_index].points.front();
        const Result<SurfacePoint> sample = surface_eval(surface, first[0], first[1]);
        if (!sample.has_value()) return sample.error();
        if (!has_reference) {
          reference = sample.value().p;
          has_reference = true;
        }
        const Vec3d normal = cross(sample.value().du, sample.value().dv);
        const Vec3d position{sample.value().p[0] - reference[0],
                             sample.value().p[1] - reference[1],
                             sample.value().p[2] - reference[2]};
        const double face_volume = area.value() * dot(position, normal) / 3.0;
        total += face.reversed ? -face_volume : face_volume;
        continue;
      }

      Polygon2d polygon;
      polygon.ring_starts.push_back(0u);
      const auto append_ring = [&polygon](const FaceUvLoop& loop) {
        polygon.vertices.insert(polygon.vertices.end(), loop.points.begin(),
                                loop.points.end());
        polygon.ring_starts.push_back(static_cast<std::uint32_t>(polygon.vertices.size()));
      };
      append_ring(loops[outer_index]);
      for (std::size_t i = 0; i < loops.size(); ++i) {
        if (i != outer_index) append_ring(loops[i]);
      }

      // Translate the polygon so that its first vertex is the origin. The
      // shoelace sums and the triangulation predicates then work on the small
      // offsets instead of the absolute parameters, which may be of order 1e9.
      const Point2d shift = polygon.vertices.front();
      for (Point2d& point : polygon.vertices) {
        point[0] -= shift[0];
        point[1] -= shift[1];
      }

      const Triangulation2 triangulation = triangulate_polygon2(polygon);
      if (!triangulation.ok) {
        return Error{Errc::degenerate, "curved_solid_volume: face triangulation failed",
                     face_id.index(), 0u};
      }

      double face_volume = 0.0;

      // Integrates the seven-point rule over one uv subtriangle. `area` is the
      // subtriangle's parameter-domain area. The shift is added back only for
      // the surface evaluation, so the parameterization is unchanged.
      const auto integrate = [&](const Point2d& p0, const Point2d& p1, const Point2d& p2,
                                 double area) -> Result<double> {
        double sum = 0.0;
        for (const CurvedVolumeQuadraturePoint& point : kCurvedVolumeQuadraturePoints) {
          const double u = point.a * p0[0] + point.b * p1[0] + point.c * p2[0] + shift[0];
          const double v = point.a * p0[1] + point.b * p1[1] + point.c * p2[1] + shift[1];
          const Result<SurfacePoint> sample = surface_eval(surface, u, v);
          if (!sample.has_value()) return sample.error();

          const SurfacePoint& evaluated = sample.value();
          if (!has_reference) {
            reference = evaluated.p;
            has_reference = true;
          }
          const Vec3d position{evaluated.p[0] - reference[0], evaluated.p[1] - reference[1],
                               evaluated.p[2] - reference[2]};
          const Vec3d normal = cross(evaluated.du, evaluated.dv);
          sum += area * point.weight * dot(position, normal) / 3.0;
        }
        return sum;
      };

      for (const std::array<std::uint32_t, 3>& triangle : triangulation.triangles) {
        const Point2d& a = polygon.vertices[triangle[0]];
        const Point2d& b = polygon.vertices[triangle[1]];
        const Point2d& c = polygon.vertices[triangle[2]];
        const double twice_area =
            (b[0] - a[0]) * (c[1] - a[1]) - (b[1] - a[1]) * (c[0] - a[0]);
        if (!(twice_area > 0.0)) {
          return Error{Errc::degenerate, "curved_solid_volume: face triangulation failed",
                       face_id.index(), 0u};
        }
        const double area_uv = 0.5 * twice_area;

        const auto edge_length = [](const Point2d& p, const Point2d& q) {
          return std::hypot(q[0] - p[0], q[1] - p[1]);
        };
        const int level =
            planar ? 0
                   : curved_volume_subdivision_level(std::max({edge_length(a, b),
                                                               edge_length(b, c),
                                                               edge_length(c, a)}));
        if (level == 0) {
          const Result<double> contribution = integrate(a, b, c, area_uv);
          if (!contribution.has_value()) return contribution.error();
          face_volume += contribution.value();
          continue;
        }

        // Barycentric grid with 2^level parts per edge. Each subtriangle has
        // area_uv / parts^2 and the same orientation as its parent.
        const int parts = 1 << level;
        const double step = 1.0 / static_cast<double>(parts);
        const double sub_area = area_uv * step * step;
        const auto grid_point = [&](int i, int j) {
          const double fi = static_cast<double>(i) * step;
          const double fj = static_cast<double>(j) * step;
          return Point2d{a[0] + fi * (b[0] - a[0]) + fj * (c[0] - a[0]),
                         a[1] + fi * (b[1] - a[1]) + fj * (c[1] - a[1])};
        };
        for (int i = 0; i < parts; ++i) {
          for (int j = 0; j < parts - i; ++j) {
            const Point2d p00 = grid_point(i, j);
            const Point2d p10 = grid_point(i + 1, j);
            const Point2d p01 = grid_point(i, j + 1);
            const Result<double> first = integrate(p00, p10, p01, sub_area);
            if (!first.has_value()) return first.error();
            face_volume += first.value();
            if (i + j + 1 < parts) {
              const Point2d p11 = grid_point(i + 1, j + 1);
              const Result<double> second = integrate(p10, p11, p01, sub_area);
              if (!second.has_value()) return second.error();
              face_volume += second.value();
            }
          }
        }
      }

      total += face.reversed ? -face_volume : face_volume;
    }
  }

  if (!std::isfinite(total)) {
    return Error{Errc::numeric_failure, "curved_solid_volume: volume is not finite",
                 solid.index(), 0u};
  }
  if (!(total > 0.0)) {
    return Error{Errc::numeric_failure, "curved_solid_volume: volume is not positive",
                 solid.index(), 0u};
  }
  return total;
}

}  // namespace detail
}  // namespace v1
}  // namespace geom
