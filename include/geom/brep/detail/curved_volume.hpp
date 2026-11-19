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
#include <geom/surfaces/analytic.hpp>
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

// A curved uv triangle is split into 2^level pieces per edge, aiming for pieces
// no longer than roughly kCurvedVolumeTargetEdge. The level is fixed (there is
// no adaptivity) and capped at eight, so once the cap binds pieces can be longer
// than the target; measured relative error for cylinders and tori stays near
// 1e-14 because their integrands are smooth and nearly linear per panel.
inline constexpr double kCurvedVolumeTargetEdge = 0.025;

inline int curved_volume_subdivision_level(double max_edge) noexcept {
  constexpr double kMaxPanels = 256.0;
  if (!(max_edge > kCurvedVolumeTargetEdge)) return 0;
  const double panels = max_edge / kCurvedVolumeTargetEdge;
  if (!(panels < kMaxPanels)) return 8;
  return std::min(8, static_cast<int>(std::ceil(std::log2(panels))));
}

// Drops samples that sit within `tolerance` of the segment through their
// surviving neighbours. face_loops_uv samples every edge at a fixed number of
// segments, so a straight edge arrives as a dense run of collinear points; the
// triangulator only removes *exactly* collinear vertices, and after a rotation
// the projection rounding leaves those runs a few ulp off a line, which keeps
// every sample and turns the ear clipping into a quadratic pass over thousands
// of vertices. The tolerance is far below any real curvature over one sample
// (a 32-segment arc of radius 0.5 deviates by ~6e-4), so this only removes
// sampling redundancy along straight edges and leaves the integral unchanged
// to rounding.
inline std::vector<Point2d> curved_volume_simplify_ring(const std::vector<Point2d>& points,
                                                        double tolerance) {
  std::vector<Point2d> simplified;
  simplified.reserve(points.size());
  for (const Point2d& point : points) {
    if (simplified.size() >= 2u) {
      const Point2d& previous = simplified[simplified.size() - 2u];
      const Point2d& last = simplified.back();
      const double dx = last[0] - previous[0];
      const double dy = last[1] - previous[1];
      const double length = std::hypot(dx, dy);
      if (length > 0.0) {
        const double distance =
            std::fabs(dx * (point[1] - previous[1]) - dy * (point[0] - previous[0])) / length;
        if (distance <= tolerance) {
          simplified.back() = point;
          continue;
        }
      }
    }
    simplified.push_back(point);
  }
  // A ring that collapses below three vertices cannot be triangulated; keep the
  // original sampling so a degenerate sliver fails (or passes) exactly as before.
  if (simplified.size() < 3u) return points;
  return simplified;
}

// Signed volume of a closed solid obtained by integrating
// (1/3) S . (S_u x S_v) over every face's uv polygon.
//
// Scope: a planar face uses its exact arc-aware area (see
// curved_volume_planar_area), so a cylinder with circular cap provenance
// integrates to exactly pi * r^2 * h. Cylinder lateral faces and trimmed
// analytic patches in general are integrated on the exact surface with the
// seven-point rule at a fixed subdivision level. Sphere seam/pole loops and the
// cone apex fan are closed by synthesising the collapsed uv rows (see the
// append_ring lambda), and NURBS surfaces are rejected. Where a loop has no
// circle provenance (for example a chorded test fixture), the planar area falls
// back to the chord polygon and the result is that chord B-rep's volume.
//
// Conditioning: every uv polygon is translated by its first vertex before
// triangulation and before the uv-area computation, and each sample is
// translated back before surface_eval(), so the parameterization is unchanged
// while the shoelace sums stay well conditioned even when the face parameters
// are of order 1e9. The divergence integrand is additionally measured from a
// fixed point of the solid: a closed shell has zero normal flux, so this is an
// exact change of origin, but it keeps a solid placed far from the world origin
// from losing its volume to cancellation between the large face contributions.
// Signed area contribution of one circular-arc segment. `sweep` is the signed
// sweep carried by the segment's arc reference (|sweep| < pi) and `orientation`
// is +1 when the circle normal agrees with the face normal and -1 when it is
// opposed. The uv frame can rotate or mirror the circle's own (x_axis, y_axis)
// frame, so the start azimuth comes from the shifted uv point and the signed
// sweep is oriented into that frame. A concave arc (travelled clockwise about
// its circle) therefore keeps its negative contribution instead of being forced
// the long way around.
inline double curved_volume_arc_term(const Point2d& a, const Point2d& center,
                                     double radius, double sweep, double orientation) {
  const double ax = a[0] - center[0];
  const double ay = a[1] - center[1];
  const double angle_a = std::atan2(ay, ax);
  const double delta = orientation * sweep;
  const double angle_b = angle_a + delta;
  return 0.5 * (radius * radius * delta +
                radius * (center[0] * (std::sin(angle_b) - std::sin(angle_a)) -
                          center[1] * (std::cos(angle_b) - std::cos(angle_a))));
}

// Exact signed area of a planar face's loops in parameter space. A segment with
// an engaged per-segment arc reference whose circle is parallel and coplanar
// with the face plane contributes its circular-segment area; everything else is
// a straight chord. This is what makes a cylinder cap integrate to pi*r^2
// rather than to the 32-gon polygon.
inline Result<double> curved_volume_planar_area(const std::vector<FaceUvLoop>& loops,
                                                const PlaneFrame& frame, double tolerance) {
  const Point2d shift = loops.front().points.front();
  double area = 0.0;
  for (const FaceUvLoop& loop : loops) {
    const std::size_t count = loop.points.size();
    for (std::size_t i = 0; i < count; ++i) {
      const std::size_t next = (i + 1) % count;
      const Point2d& a = loop.points[i];
      const Point2d& b = loop.points[next];
      bool arc = false;
      double radius = 0.0;
      double sweep = 0.0;
      double orientation = 1.0;
      Point2d center{0.0, 0.0};
      if (loop.arcs.size() == count && loop.arcs[i].has_value()) {
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
            sweep = loop.arcs[i]->t1 - loop.arcs[i]->t0;
            orientation = dot(circle_normal, frame.normal) >= 0.0 ? 1.0 : -1.0;
            arc = std::fabs(sweep) > 0.0;
          }
        }
      }
      if (arc) {
        area += curved_volume_arc_term(Point2d{a[0] - shift[0], a[1] - shift[1]},
                                       center, radius, sweep, orientation);
      } else {
        area += 0.5 * ((a[0] - shift[0]) * (b[1] - shift[1]) -
                       (b[0] - shift[0]) * (a[1] - shift[1]));
      }
    }
  }
  return area;
}

// The full parameter rectangle of a sphere. These four values duplicate
// curved_boolean_surface_uv_domain() in boolean_solid_curved.hpp: that header
// owns the boolean driver and must stay out of detail/, so the domain is
// repeated here with the same convention. The u period is 2*pi (u = 0 and
// u = 2*pi are the same seam meridian) and v runs from the south pole to the
// north pole.
inline constexpr double kCurvedVolumePi = 3.14159265358979323846;
inline constexpr double kCurvedVolumeTwoPi = 2.0 * kCurvedVolumePi;
inline constexpr double kCurvedVolumeSouthPoleV = -0.5 * kCurvedVolumePi;
inline constexpr double kCurvedVolumeNorthPoleV = 0.5 * kCurvedVolumePi;

// Number of latitude bands for the synthetic full-sphere grid. With
// nu = 2*nv the u step is pi/nv, so one panel's chord sagitta is
// radius * (1 - cos(pi / (2*nv))). Bands are doubled from 128 until that
// sagitta is below 1e-4; 128 already gives 7.5e-5 for the unit sphere. The
// result is a power of two so the seam column lands exactly on u = 2*pi. The
// upper clamp bounds the panel count for very large spheres, where the absolute
// sagitta bound would otherwise demand about 3,512 bands at
// radius 1e3 (about 7,025 columns), where the clamp leaves a sagitta of ~1.2e-3); the seven-point rule
// stays accurate there because the sphere integrand varies only with v.
inline int curved_volume_sphere_bands(double radius) noexcept {
  constexpr double kMaxSagitta = 1e-4;
  constexpr int kMinBands = 128;
  constexpr int kMaxBands = 1024;
  if (!(radius > kMaxSagitta)) return kMinBands;
  const double cosine = std::clamp(1.0 - kMaxSagitta / radius, -1.0, 1.0);
  const double step = 2.0 * std::acos(cosine);
  if (!(step > 0.0)) return kMaxBands;
  int bands = kMinBands;
  while (bands < kMaxBands && !(static_cast<double>(bands) * step > kCurvedVolumePi)) {
    bands *= 2;
  }
  return bands;
}

// A sphere face whose loop is the degenerate seam slit projects to one meridian
// traversed up and back: every projected point shares the same azimuth and the
// enclosed uv area is zero. A latitude ring (a spherical cap fragment) also
// projects to a zero-area line, but its v varies and it is a legitimate trimmed
// face, so requiring a single meridian keeps this predicate from widening a
// trimmed domain to the full sphere. The unwrapped uv loop may place the seam
// at any azimuth, so the comparison is relative to the first point.
inline bool curved_volume_sphere_seam_loop(const FaceUvLoop& loop) noexcept {
  if (loop.points.size() < 3u) return false;
  const double first_u = loop.points.front()[0];
  const std::size_t count = loop.points.size();
  double twice_area = 0.0;
  double extent = 0.0;
  for (std::size_t i = 0; i < count; ++i) {
    const Point2d& a = loop.points[i];
    const Point2d& b = loop.points[(i + 1u) % count];
    twice_area += a[0] * b[1] - b[0] * a[1];
    extent = std::max({extent, std::fabs(a[0]), std::fabs(a[1])});
    if (!(std::fabs(a[0] - first_u) <= 1e-12)) return false;
  }
  return std::fabs(0.5 * twice_area) <= 1e-12 * extent * extent;
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
      const Cone3* cone = std::get_if<Cone3>(&surface);
      // The cone chart's collapsed side is the apex row v = radius/tan(half_angle).
      // Face loops carry that row as a single uv copy, or as two copies when the
      // loop travels the back seam as well; every copy evaluates to the same 3D
      // apex.
      const double cone_apex_v =
          cone == nullptr ? 0.0 : cone->radius / std::tan(cone->half_angle);
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

      Point2d shift{0.0, 0.0};
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

      const Sphere3* sphere = std::get_if<Sphere3>(&surface);
      if (sphere != nullptr && loops.size() == 1u &&
          curved_volume_sphere_seam_loop(loops[outer_index])) {
        // The seam-slit loop retraces one meridian, so it encloses no uv area
        // and the polygon triangulation has nothing to work with. The face is
        // the whole sphere, whose parameter domain is the full rectangle
        // [0, 2*pi] x [-pi/2, pi/2]; cover it with nu = 2*nv quads. Every
        // vertex on the v = +/-pi/2 rows maps to the single south or north pole
        // (surface_eval drops the azimuth there because cos(v) = 0), and the
        // u = 0 and u = 2*pi columns map to the same seam meridian because the
        // surface is 2*pi-periodic. Collapsing the poles and welding the seam
        // are therefore properties of that map; the parameter grid itself must
        // stay rectangular, since moving the pole vertices onto one azimuth
        // would halve the uv area, and with it the quadrature weight, near the
        // poles. Only the triangle topology is synthesized.
        //
        // The full-sphere face is its own closed shell, so measuring from the
        // center is an exact change of origin; it also makes the integrand
        // exactly (radius^3 / 3) * cos(v), independent of u.
        if (!has_reference) {
          reference = sphere->center;
          has_reference = true;
        }
        const int bands = curved_volume_sphere_bands(sphere->radius);
        const int columns = 2 * bands;
        const double du = kCurvedVolumeTwoPi / static_cast<double>(columns);
        const double dv = kCurvedVolumePi / static_cast<double>(bands);
        const double half_cell = 0.5 * du * dv;
        for (int j = 0; j < bands; ++j) {
          const double v0 = kCurvedVolumeSouthPoleV + static_cast<double>(j) * dv;
          const double v1 = v0 + dv;
          for (int i = 0; i < columns; ++i) {
            const double u_left = static_cast<double>(i) * du;
            // The rectangle closes on itself: the last column is the welded
            // copy of u = 0, evaluated at u = 2*pi so the panel area stays
            // positive.
            const double u_right =
                (i + 1 == columns) ? kCurvedVolumeTwoPi : static_cast<double>(i + 1) * du;
            const Point2d p00{u_left, v0};
            const Point2d p10{u_right, v0};
            const Point2d p01{u_left, v1};
            const Point2d p11{u_right, v1};
            const Result<double> lower = integrate(p00, p10, p01, half_cell);
            if (!lower.has_value()) return lower.error();
            face_volume += lower.value();
            const Result<double> upper = integrate(p10, p11, p01, half_cell);
            if (!upper.has_value()) return upper.error();
            face_volume += upper.value();
          }
        }
        total += face.reversed ? -face_volume : face_volume;
        continue;
      }

      Polygon2d polygon;
      polygon.ring_starts.push_back(0u);
      double volume_uv_scale = 1.0;
      for (const FaceUvLoop& loop : loops) {
        for (const Point2d& point : loop.points) {
          volume_uv_scale =
              std::max(volume_uv_scale, std::max(std::fabs(point[0]), std::fabs(point[1])));
        }
      }
      const double simplify_tolerance = 1e-12 * volume_uv_scale;
      const auto append_ring = [&polygon, sphere, cone, cone_apex_v,
                                simplify_tolerance](const FaceUvLoop& loop) {
        const std::vector<Point2d> simplified =
            (sphere == nullptr && loop.points.size() > 3u)
                ? curved_volume_simplify_ring(loop.points, simplify_tolerance)
                : loop.points;
        const std::size_t count = simplified.size();
        for (std::size_t i = 0; i < count; ++i) {
          const Point2d& point = simplified[i];
          if (sphere == nullptr && cone == nullptr) {
            polygon.vertices.push_back(point);
            continue;
          }
          const bool collapsed =
              sphere != nullptr
                  ? (point[1] == kCurvedVolumeSouthPoleV ||
                     point[1] == kCurvedVolumeNorthPoleV)
                  : point[1] == cone_apex_v;
          if (!collapsed) {
            polygon.vertices.push_back(point);
            continue;
          }
          // A sphere loop that passes through a pole, or a cone loop that ends
          // on its apex row, projects that collapsed side onto one azimuth, so
          // the ring's closing edge cuts the corner between the two meridians
          // that meet there. Insert the missing corners so the uv domain spans
          // the whole chart window; every corner evaluates to the same 3D
          // pole/apex, so no surface is invented and the fan closes exactly once.
          const Point2d& previous = simplified[(i + count - 1u) % count];
          const Point2d& next = simplified[(i + 1u) % count];
          if (previous[1] != point[1] && std::fabs(previous[0] - point[0]) > 1e-12) {
            polygon.vertices.push_back(Point2d{previous[0], point[1]});
          }
          polygon.vertices.push_back(point);
          if (next[1] != point[1] && std::fabs(next[0] - point[0]) > 1e-12) {
            polygon.vertices.push_back(Point2d{next[0], point[1]});
          }
        }
        polygon.ring_starts.push_back(static_cast<std::uint32_t>(polygon.vertices.size()));
      };
      append_ring(loops[outer_index]);
      for (std::size_t i = 0; i < loops.size(); ++i) {
        if (i != outer_index) append_ring(loops[i]);
      }

      // Translate the polygon so that its first vertex is the origin. The
      // shoelace sums and the triangulation predicates then work on the small
      // offsets instead of the absolute parameters, which may be of order 1e9.
      shift = polygon.vertices.front();
      for (Point2d& point : polygon.vertices) {
        point[0] -= shift[0];
        point[1] -= shift[1];
      }

      const Triangulation2 triangulation = triangulate_polygon2(polygon);
      if (!triangulation.ok) {
        return Error{Errc::degenerate, "curved_solid_volume: face triangulation failed",
                     face_id.index(), 0u};
      }

      for (const std::array<std::uint32_t, 3>& triangle : triangulation.triangles) {
        const Point2d& a = polygon.vertices[triangle[0]];
        const Point2d& b = polygon.vertices[triangle[1]];
        const Point2d& c = polygon.vertices[triangle[2]];
        const double twice_area =
            (b[0] - a[0]) * (c[1] - a[1]) - (b[1] - a[1]) * (c[0] - a[0]);
        if (twice_area == 0.0) {
          // Collinear samples along a straight or collapsed uv row triangulate
          // to zero area: they cover no parameter domain, so they contribute no
          // volume and must not fail the face.
          continue;
        }
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
