#pragma once

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

// Signed volume of a closed solid obtained by integrating
// (1/3) S . (S_u x S_v) over every face's exact uv polygon. Planar faces are
// integrated exactly by the quadrature; curved faces use their analytic
// parameterization evaluated at the same seven points, so a trimmed face is
// handled without ever tessellating it in three dimensions.
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

      const Triangulation2 triangulation = triangulate_polygon2(polygon);
      if (!triangulation.ok) {
        return Error{Errc::degenerate, "curved_solid_volume: face triangulation failed",
                     face_id.index(), 0u};
      }

      const SurfaceGeometry& surface = brep.geometry().surfaces[face.surface];
      double face_volume = 0.0;
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

        for (const CurvedVolumeQuadraturePoint& point : kCurvedVolumeQuadraturePoints) {
          const double u = point.a * a[0] + point.b * b[0] + point.c * c[0];
          const double v = point.a * a[1] + point.b * b[1] + point.c * c[1];
          const Result<SurfacePoint> sample = surface_eval(surface, u, v);
          if (!sample.has_value()) return sample.error();

          const SurfacePoint& evaluated = sample.value();
          const Vec3d position{evaluated.p[0], evaluated.p[1], evaluated.p[2]};
          const Vec3d normal = cross(evaluated.du, evaluated.dv);
          face_volume += area_uv * point.weight * dot(position, normal) / 3.0;
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
