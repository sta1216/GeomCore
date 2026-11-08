#pragma once

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <numbers>
#include <vector>

#include <geom/brep/tessellate.hpp>
#include <geom/brep/topology.hpp>
#include <geom/core/error.hpp>
#include <geom/core/result.hpp>
#include <geom/curves/analytic_curve.hpp>
#include <geom/primitives/point.hpp>
#include <geom/primitives/vec.hpp>
#include <geom/surfaces/analytic.hpp>
#include <geom/surfaces/plane.hpp>

namespace geom {
inline namespace v1 {

// Unit-density mass properties of a closed B-rep, about the centroid.
//
// The volume, centroid and inertia tensor are surface integrals obtained from the
// divergence theorem:
//   V   = (1/3)  ∮ x·n dA
//   C   = (1/(4V)) ∮ x (x·n) dA
//   I_O = ∮ (|x|² 1 - x⊗x) (x·n) / 5 dA,  I_C = I_O - V (|C|² 1 - C⊗C)
// The (x·n) factor is folded into the element area, so the integrals are
// orientation-sensitive: reversed faces flip the sign of dA.
//
// Every shell must be flagged closed (otherwise Errc::not_closed), and the
// accumulated volume must be positive (otherwise Errc::numeric_failure), which
// rejects open or inward-oriented shells instead of returning a plausible number.
//
// Planar faces are sampled through their loop boundary at `sampling_deflection`
// and integrated with the exact polygon area moments (shoelace plus first and
// second moments); holes are the loops stored with the opposite winding, so
// their signed contribution is negative. Curved faces are integrated over their
// (u,v) parameter rectangle with `quadrature_order`-point Gauss-Legendre, using
// dA = cross(du, dv) du dv. A curved face is integrated only when its boundary
// is provably the untrimmed rectangle of the analytic surface: the loop must
// have the exact seam/circle coedge pattern the M2 primitives produce (full
// period in u, two constant-v circles for a cylinder band, a meridian for a
// sphere, two closed circles for a torus). Any other curved face returns
// Errc::unsupported rather than integrating the wrong region.
//
// All output values are finite; a non-finite integral is Errc::numeric_failure.
// Options: a non-positive or non-finite sampling_deflection and a
// quadrature_order below 1 are Errc::tolerance_violation.

struct MassProperties {
  double volume = 0.0;
  Point3d centroid{};
  double inertia[3][3]{};  // About the centroid, unit density.
};

struct MassPropertiesOptions {
  double sampling_deflection = 1e-9;  // Planar face boundary chord height.
  int quadrature_order = 24;          // Gauss-Legendre points per direction.
};

namespace detail {

struct MassIntegrals {
  double volume = 0.0;
  Vec3d first{};  // ∮ x (x·n) dA
  double second[3][3]{};  // I_O = ∮ (|x|² 1 - x⊗x) (x·n) dA
};

struct MassUvDomain {
  double u0 = 0.0;
  double u1 = 0.0;
  double v0 = 0.0;
  double v1 = 0.0;
};

inline double mass_pi() noexcept { return std::numbers::pi_v<double>; }

inline double mass_two_pi() noexcept { return 2.0 * std::numbers::pi_v<double>; }

inline Error mass_unsupported_rectangle(FaceId id) {
  return Error{Errc::unsupported,
               "mass_properties: curved face is not an untrimmed full rectangle", id.index(), 0u};
}

inline Error mass_not_closed(ShellId id) {
  return Error{Errc::not_closed, "mass_properties: shell is not closed", id.index(), 0u};
}

inline Error mass_numeric_failure(const char* message, std::uint32_t entity = 0u) {
  return Error{Errc::numeric_failure, message, entity, 0u};
}

inline double mass_geometry_tolerance(double deflection, double scale) noexcept {
  const double relative = 64.0 * std::numeric_limits<double>::epsilon() * (1.0 + std::fabs(scale));
  return std::max(deflection, relative);
}

inline const CircleCurve3* mass_circle_curve(const Brep& brep, EdgeId id) noexcept {
  const Edge& edge = brep.edge(id);
  if (edge.curve == kNoGeometry || edge.curve >= brep.geometry().curves.size()) return nullptr;
  return std::get_if<CircleCurve3>(&brep.geometry().curves[edge.curve]);
}

// One distinct edge of a face loop together with how often and in which
// directions the loop traverses it. Seams are used once per direction, cap
// circles once in total.
struct MassEdgeUse {
  EdgeId edge{};
  int count = 0;
  int forward = 0;
  int reversed = 0;
};

inline std::vector<MassEdgeUse> mass_edge_uses(const Brep& brep,
                                               const std::vector<CoedgeId>& coedges) {
  std::vector<MassEdgeUse> uses;
  for (const CoedgeId id : coedges) {
    const Coedge& coedge = brep.coedge(id);
    MassEdgeUse* use = nullptr;
    for (MassEdgeUse& candidate : uses) {
      if (candidate.edge == coedge.edge) {
        use = &candidate;
        break;
      }
    }
    if (use == nullptr) {
      uses.push_back(MassEdgeUse{coedge.edge, 0, 0, 0});
      use = &uses.back();
    }
    ++use->count;
    if (coedge.reversed) {
      ++use->reversed;
    } else {
      ++use->forward;
    }
  }
  return uses;
}

// Coedges of the face's single loop, in traversal order. A face with more than
// one loop is not the untrimmed rectangle of an analytic surface.
inline Result<std::vector<CoedgeId>> mass_single_loop_coedges(const Brep& brep, FaceId id) {
  const Face& face = brep.face(id);
  if (face.first_loop.is_null()) {
    return Error{Errc::invalid_topology, "mass_properties: face has no loop", id.index(), 0u};
  }
  if (!brep.loop(face.first_loop).next.is_null()) return mass_unsupported_rectangle(id);
  std::vector<CoedgeId> coedges;
  const std::size_t limit = brep.coedges().size();
  CoedgeId current = brep.loop(face.first_loop).start;
  if (current.is_null()) {
    return Error{Errc::invalid_topology, "mass_properties: loop has no coedges", id.index(), 0u};
  }
  do {
    coedges.push_back(current);
    if (coedges.size() > limit) {
      return Error{Errc::invalid_topology, "mass_properties: coedge cycle is cyclic", id.index(),
                   0u};
    }
    current = brep.coedge(current).next;
  } while (!current.is_null() && !(current == coedges.front()));
  if (current.is_null()) {
    return Error{Errc::invalid_topology, "mass_properties: loop is not closed", id.index(), 0u};
  }
  return coedges;
}

inline Vec3d mass_negate(Vec3d v) noexcept { return Vec3d{-v[0], -v[1], -v[2]}; }

// Adds one oriented surface element dA (already including du dv and the face
// reversal) with quadrature weight w.
inline void mass_add_surface_element(MassIntegrals& acc, Point3d point, Vec3d area,
                                     double weight) noexcept {
  const Vec3d x{point[0], point[1], point[2]};
  const double flux = dot(x, area);
  acc.volume += flux * weight / 3.0;
  acc.first = acc.first + x * (flux * weight);
  const double x2 = dot(x, x);
  for (int i = 0; i < 3; ++i) {
    for (int j = 0; j < 3; ++j) {
      const double delta = (i == j) ? x2 : 0.0;
      acc.second[i][j] += (delta - x[i] * x[j]) * flux * weight / 5.0;
    }
  }
}

inline void mass_gauss_legendre(int order, std::vector<double>& nodes,
                                std::vector<double>& weights) {
  nodes.assign(static_cast<std::size_t>(order), 0.0);
  weights.assign(static_cast<std::size_t>(order), 0.0);
  const int half = (order + 1) / 2;
  for (int i = 0; i < half; ++i) {
    double x = std::cos(std::numbers::pi_v<double> * (static_cast<double>(i) + 0.75) /
                        (static_cast<double>(order) + 0.5));
    for (int iteration = 0; iteration < 100; ++iteration) {
      double previous = 1.0;
      double current = x;
      for (int degree = 2; degree <= order; ++degree) {
        const double next = ((2.0 * degree - 1.0) * x * current - (degree - 1.0) * previous) /
                            static_cast<double>(degree);
        previous = current;
        current = next;
      }
      const double derivative = static_cast<double>(order) * (x * current - previous) /
                                (x * x - 1.0);
      const double step = current / derivative;
      x -= step;
      if (std::fabs(step) <= 1e-15) break;
    }
    double previous = 1.0;
    double current = x;
    for (int degree = 2; degree <= order; ++degree) {
      const double next = ((2.0 * degree - 1.0) * x * current - (degree - 1.0) * previous) /
                          static_cast<double>(degree);
      previous = current;
      current = next;
    }
    const double derivative = static_cast<double>(order) * (x * current - previous) /
                              (x * x - 1.0);
    const double weight = 2.0 / ((1.0 - x * x) * derivative * derivative);
    const std::size_t low = static_cast<std::size_t>(i);
    const std::size_t high = static_cast<std::size_t>(order - 1 - i);
    nodes[low] = -x;
    nodes[high] = x;
    weights[low] = weight;
    weights[high] = weight;
  }
}

template <class Surface>
inline void mass_integrate_parametric(const Surface& surface, const MassUvDomain& domain,
                                      bool reversed, const std::vector<double>& nodes,
                                      const std::vector<double>& weights, MassIntegrals& acc) {
  const double middle_u = 0.5 * (domain.u0 + domain.u1);
  const double middle_v = 0.5 * (domain.v0 + domain.v1);
  const double half_u = 0.5 * (domain.u1 - domain.u0);
  const double half_v = 0.5 * (domain.v1 - domain.v0);
  for (std::size_t i = 0; i < nodes.size(); ++i) {
    const double u = middle_u + half_u * nodes[i];
    for (std::size_t j = 0; j < nodes.size(); ++j) {
      const double v = middle_v + half_v * nodes[j];
      const double weight = half_u * half_v * weights[i] * weights[j];
      Vec3d area = cross(surface.du(u, v), surface.dv(u, v));
      if (reversed) area = mass_negate(area);
      mass_add_surface_element(acc, surface.eval(u, v), area, weight);
    }
  }
}

inline Result<MassUvDomain> mass_cylinder_domain(const Brep& brep, FaceId id,
                                                 const Cylinder3& cylinder,
                                                 double deflection) {
  const Result<std::vector<CoedgeId>> coedges = mass_single_loop_coedges(brep, id);
  if (!coedges.has_value()) return coedges.error();
  if (coedges.value().size() != 4) return mass_unsupported_rectangle(id);
  const std::vector<MassEdgeUse> uses = mass_edge_uses(brep, coedges.value());
  const MassEdgeUse* seam = nullptr;
  std::vector<const MassEdgeUse*> circles;
  for (const MassEdgeUse& use : uses) {
    if (use.count == 2 && use.forward == 1 && use.reversed == 1) {
      if (seam != nullptr) return mass_unsupported_rectangle(id);
      seam = &use;
    } else if (use.count == 1) {
      circles.push_back(&use);
    } else {
      return mass_unsupported_rectangle(id);
    }
  }
  if (seam == nullptr || circles.size() != 2) return mass_unsupported_rectangle(id);

  double levels[2] = {0.0, 0.0};
  for (int i = 0; i < 2; ++i) {
    const Edge& edge = brep.edge(circles[static_cast<std::size_t>(i)]->edge);
    if (!(edge.v0 == edge.v1)) return mass_unsupported_rectangle(id);
    const CircleCurve3* circle = mass_circle_curve(brep, circles[static_cast<std::size_t>(i)]->edge);
    if (circle == nullptr) return mass_unsupported_rectangle(id);
    const Vec3d relative = circle->center - cylinder.origin;
    const double axial = dot(relative, cylinder.axis);
    const double tolerance = mass_geometry_tolerance(
        deflection, std::fabs(axial) + std::fabs(circle->radius) + std::fabs(cylinder.radius));
    if (norm2(relative - cylinder.axis * axial) > tolerance) return mass_unsupported_rectangle(id);
    if (std::fabs(circle->radius - cylinder.radius) > tolerance) {
      return mass_unsupported_rectangle(id);
    }
    levels[i] = axial;
  }

  const double v0 = std::min(levels[0], levels[1]);
  const double v1 = std::max(levels[0], levels[1]);
  const double tolerance = mass_geometry_tolerance(
      deflection, std::fabs(v0) + std::fabs(v1) + std::fabs(cylinder.radius));
  if (!(v1 - v0 > tolerance)) return mass_unsupported_rectangle(id);

  const Edge& seam_edge = brep.edge(seam->edge);
  if (seam_edge.v0 == seam_edge.v1) return mass_unsupported_rectangle(id);
  const double start = dot(brep.vertex(seam_edge.v0).point - cylinder.origin, cylinder.axis);
  const double end = dot(brep.vertex(seam_edge.v1).point - cylinder.origin, cylinder.axis);
  const bool start_on_ring = std::fabs(start - v0) <= tolerance || std::fabs(start - v1) <= tolerance;
  const bool end_on_ring = std::fabs(end - v0) <= tolerance || std::fabs(end - v1) <= tolerance;
  if (!start_on_ring || !end_on_ring) return mass_unsupported_rectangle(id);
  return MassUvDomain{0.0, mass_two_pi(), v0, v1};
}

inline Result<MassUvDomain> mass_cone_domain(const Brep& brep, FaceId id, const Cone3& cone,
                                             double deflection) {
  const Result<std::vector<CoedgeId>> coedges = mass_single_loop_coedges(brep, id);
  if (!coedges.has_value()) return coedges.error();
  if (coedges.value().size() != 3) return mass_unsupported_rectangle(id);
  const std::vector<MassEdgeUse> uses = mass_edge_uses(brep, coedges.value());
  const MassEdgeUse* seam = nullptr;
  const MassEdgeUse* base = nullptr;
  for (const MassEdgeUse& use : uses) {
    if (use.count == 2 && use.forward == 1 && use.reversed == 1) {
      if (seam != nullptr) return mass_unsupported_rectangle(id);
      seam = &use;
    } else if (use.count == 1) {
      if (base != nullptr) return mass_unsupported_rectangle(id);
      base = &use;
    } else {
      return mass_unsupported_rectangle(id);
    }
  }
  if (seam == nullptr || base == nullptr) return mass_unsupported_rectangle(id);

  const Edge& base_edge = brep.edge(base->edge);
  if (!(base_edge.v0 == base_edge.v1)) return mass_unsupported_rectangle(id);
  const CircleCurve3* base_circle = mass_circle_curve(brep, base->edge);
  if (base_circle == nullptr) return mass_unsupported_rectangle(id);
  const Vec3d base_relative = base_circle->center - cone.origin;
  const double v0 = dot(base_relative, cone.axis);
  const double base_tolerance = mass_geometry_tolerance(
      deflection, std::fabs(v0) + std::fabs(base_circle->radius) + std::fabs(cone.radius));
  if (norm2(base_relative - cone.axis * v0) > base_tolerance) {
    return mass_unsupported_rectangle(id);
  }
  if (std::fabs(base_circle->radius - cone.radius) > base_tolerance) {
    return mass_unsupported_rectangle(id);
  }

  const Edge& seam_edge = brep.edge(seam->edge);
  if (seam_edge.v0 == seam_edge.v1) return mass_unsupported_rectangle(id);
  const double first = dot(brep.vertex(seam_edge.v0).point - cone.origin, cone.axis);
  const double second = dot(brep.vertex(seam_edge.v1).point - cone.origin, cone.axis);
  const double tolerance = mass_geometry_tolerance(deflection, std::fabs(first) + std::fabs(second));
  const bool first_is_base = std::fabs(first - v0) <= tolerance;
  const bool second_is_base = std::fabs(second - v0) <= tolerance;
  if (first_is_base == second_is_base) return mass_unsupported_rectangle(id);
  const double v1 = first_is_base ? second : first;
  if (!(v1 > v0 + tolerance)) return mass_unsupported_rectangle(id);

  // The face must stop at or before the apex; beyond it the parameterization
  // continues on the opposite nappe and a rectangle integral would be wrong.
  const double tangent = std::tan(cone.half_angle);
  if (!(tangent > 0.0) || !std::isfinite(tangent)) return mass_unsupported_rectangle(id);
  const double apex = cone.radius / tangent;
  if (!(v1 <= apex + tolerance)) return mass_unsupported_rectangle(id);
  return MassUvDomain{0.0, mass_two_pi(), v0, v1};
}

inline Result<MassUvDomain> mass_sphere_domain(const Brep& brep, FaceId id, const Sphere3& sphere,
                                               double deflection) {
  const Result<std::vector<CoedgeId>> coedges = mass_single_loop_coedges(brep, id);
  if (!coedges.has_value()) return coedges.error();
  if (coedges.value().size() != 2) return mass_unsupported_rectangle(id);
  const std::vector<MassEdgeUse> uses = mass_edge_uses(brep, coedges.value());
  if (uses.size() != 1) return mass_unsupported_rectangle(id);
  if (uses[0].count != 2 || uses[0].forward != 1 || uses[0].reversed != 1) {
    return mass_unsupported_rectangle(id);
  }

  const Edge& edge = brep.edge(uses[0].edge);
  if (edge.v0 == edge.v1) return mass_unsupported_rectangle(id);
  const CircleCurve3* circle = mass_circle_curve(brep, uses[0].edge);
  if (circle == nullptr) return mass_unsupported_rectangle(id);
  const double scale = std::fabs(sphere.radius);
  const double tolerance = mass_geometry_tolerance(deflection, scale);
  if (norm2(circle->center - sphere.center) > tolerance) return mass_unsupported_rectangle(id);
  if (std::fabs(circle->radius - sphere.radius) > tolerance) return mass_unsupported_rectangle(id);
  const double span = std::fabs(edge.t1 - edge.t0);
  if (std::fabs(span - mass_pi()) > tolerance) return mass_unsupported_rectangle(id);

  const Point3d poles[2] = {brep.vertex(edge.v0).point, brep.vertex(edge.v1).point};
  for (const Point3d& pole : poles) {
    const Vec3d relative = pole - sphere.center;
    const double axial = dot(relative, sphere.axis);
    if (std::fabs(axial) < sphere.radius - tolerance) return mass_unsupported_rectangle(id);
    if (norm2(relative - sphere.axis * axial) > tolerance) return mass_unsupported_rectangle(id);
  }
  return MassUvDomain{0.0, mass_two_pi(), -mass_pi() / 2.0, mass_pi() / 2.0};
}

inline Result<MassUvDomain> mass_torus_domain(const Brep& brep, FaceId id, const Torus3& torus,
                                              double deflection) {
  const Result<std::vector<CoedgeId>> coedges = mass_single_loop_coedges(brep, id);
  if (!coedges.has_value()) return coedges.error();
  if (coedges.value().size() != 4) return mass_unsupported_rectangle(id);
  const std::vector<MassEdgeUse> uses = mass_edge_uses(brep, coedges.value());
  if (uses.size() != 2) return mass_unsupported_rectangle(id);
  for (const MassEdgeUse& use : uses) {
    if (use.count != 2 || use.forward != 1 || use.reversed != 1) {
      return mass_unsupported_rectangle(id);
    }
    const Edge& edge = brep.edge(use.edge);
    if (!(edge.v0 == edge.v1)) return mass_unsupported_rectangle(id);
    if (mass_circle_curve(brep, use.edge) == nullptr) return mass_unsupported_rectangle(id);
  }

  const double outer = torus.major_radius + torus.minor_radius;
  const double scale = std::fabs(outer) + std::fabs(torus.minor_radius);
  const double tolerance = mass_geometry_tolerance(deflection, scale);
  const Point3d tube_center = torus.center + torus.x_axis * torus.major_radius;
  bool found_outer = false;
  bool found_inner = false;
  for (const MassEdgeUse& use : uses) {
    const CircleCurve3* circle = mass_circle_curve(brep, use.edge);
    const bool is_outer = norm2(circle->center - torus.center) <= tolerance &&
                          std::fabs(circle->radius - outer) <= tolerance;
    const bool is_inner = norm2(circle->center - tube_center) <= tolerance &&
                          std::fabs(circle->radius - torus.minor_radius) <= tolerance;
    if (is_outer == is_inner) return mass_unsupported_rectangle(id);
    if (is_outer) {
      if (found_outer) return mass_unsupported_rectangle(id);
      found_outer = true;
    } else {
      if (found_inner) return mass_unsupported_rectangle(id);
      found_inner = true;
    }
  }
  if (!found_outer || !found_inner) return mass_unsupported_rectangle(id);
  return MassUvDomain{0.0, mass_two_pi(), 0.0, mass_two_pi()};
}

// Boundary polygon of one loop, in traversal order. Consecutive duplicates (the
// shared vertex between two coedges) and the repeated closing point are dropped.
inline Result<std::vector<Point3d>> mass_loop_polygon(const Brep& brep, LoopId id,
                                                      double deflection) {
  std::vector<Point3d> points;
  const std::size_t limit = brep.coedges().size();
  CoedgeId current = brep.loop(id).start;
  if (current.is_null()) {
    return Error{Errc::invalid_topology, "mass_properties: loop has no coedges", id.index(), 0u};
  }
  const CoedgeId start = current;
  std::size_t steps = 0;
  do {
    if (steps++ > limit) {
      return Error{Errc::invalid_topology, "mass_properties: coedge cycle is cyclic", id.index(),
                   0u};
    }
    const Coedge& coedge = brep.coedge(current);
    const Result<std::vector<Point3d>> samples =
        tessellate_sample_edge(brep, coedge.edge, deflection);
    if (!samples.has_value()) return samples.error();
    const std::vector<Point3d>& ordered = samples.value();
    const auto append = [&points](const Point3d& point) {
      if (!points.empty()) {
        const Point3d& last = points.back();
        if (last[0] == point[0] && last[1] == point[1] && last[2] == point[2]) return;
      }
      points.push_back(point);
    };
    if (coedge.reversed) {
      for (std::size_t i = ordered.size(); i > 0; --i) append(ordered[i - 1]);
    } else {
      for (const Point3d& point : ordered) append(point);
    }
    current = coedge.next;
  } while (!current.is_null() && !(current == start));
  if (current.is_null()) {
    return Error{Errc::invalid_topology, "mass_properties: loop is not closed", id.index(), 0u};
  }
  if (points.size() >= 2) {
    const Point3d& first = points.front();
    const Point3d& last = points.back();
    if (first[0] == last[0] && first[1] == last[1] && first[2] == last[2]) points.pop_back();
  }
  return points;
}

inline Result<void> mass_add_planar_face(const Brep& brep, FaceId id, const Plane3& plane,
                                         const MassPropertiesOptions& options,
                                         MassIntegrals& acc) {
  const double normal_length = norm2(plane.normal);
  if (!std::isfinite(normal_length) || !(normal_length > 0.0)) {
    return mass_numeric_failure("mass_properties: face plane normal is not usable", id.index());
  }
  const Vec3d normal = plane.normal * (1.0 / normal_length);
  const double offset = plane.offset / normal_length;
  // x·n is constant on the plane: the signed plane distance, flipped when the
  // face normal is opposite the stored plane normal.
  const double distance = brep.face(id).reversed ? -offset : offset;

  const Vec3d reference = std::fabs(normal[0]) < 0.9 ? Vec3d{1.0, 0.0, 0.0}
                                                     : Vec3d{0.0, 1.0, 0.0};
  const Vec3d raw_u = cross(reference, normal);
  const double raw_length = norm2(raw_u);
  if (!std::isfinite(raw_length) || !(raw_length > 0.0)) {
    return mass_numeric_failure("mass_properties: face plane frame is degenerate", id.index());
  }
  const Vec3d u_axis = raw_u * (1.0 / raw_length);
  const Vec3d v_axis = cross(normal, u_axis);
  const Vec3d origin = normal * offset;

  const std::size_t loop_limit = brep.loops().size();
  LoopId loop_id = brep.face(id).first_loop;
  if (loop_id.is_null()) {
    return Error{Errc::invalid_topology, "mass_properties: face has no loop", id.index(), 0u};
  }
  std::size_t loop_steps = 0;
  while (!loop_id.is_null()) {
    if (loop_steps++ > loop_limit) {
      return Error{Errc::invalid_topology, "mass_properties: loop chain is cyclic", id.index(), 0u};
    }
    const Result<std::vector<Point3d>> polygon =
        mass_loop_polygon(brep, loop_id, options.sampling_deflection);
    if (!polygon.has_value()) return polygon.error();

    const std::vector<Point3d>& points = polygon.value();
    const std::size_t count = points.size();
    if (count >= 3) {
      double area = 0.0;
      double moment_a = 0.0;
      double moment_b = 0.0;
      double moment_aa = 0.0;
      double moment_bb = 0.0;
      double moment_ab = 0.0;
      std::vector<std::array<double, 2>> projected(count);
      for (std::size_t i = 0; i < count; ++i) {
        const Vec3d relative{points[i][0] - origin[0], points[i][1] - origin[1],
                             points[i][2] - origin[2]};
        projected[i] = {dot(relative, u_axis), dot(relative, v_axis)};
      }
      for (std::size_t i = 0; i < count; ++i) {
        const std::array<double, 2>& p = projected[i];
        const std::array<double, 2>& q = projected[(i + 1) % count];
        const double cross2 = p[0] * q[1] - q[0] * p[1];
        area += cross2;
        moment_a += (p[0] + q[0]) * cross2;
        moment_b += (p[1] + q[1]) * cross2;
        moment_aa += (p[0] * p[0] + p[0] * q[0] + q[0] * q[0]) * cross2;
        moment_bb += (p[1] * p[1] + p[1] * q[1] + q[1] * q[1]) * cross2;
        moment_ab += (2.0 * p[0] * p[1] + p[0] * q[1] + q[0] * p[1] + 2.0 * q[0] * q[1]) * cross2;
      }
      area *= 0.5;
      moment_a /= 6.0;
      moment_b /= 6.0;
      moment_aa /= 12.0;
      moment_bb /= 12.0;
      moment_ab /= 24.0;

      // ∫ x dA and ∫ x⊗x dA in world coordinates.
      Vec3d first_moment = origin * area + u_axis * moment_a + v_axis * moment_b;
      double second_moment[3][3]{};
      for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) {
          second_moment[i][j] =
              area * origin[i] * origin[j] +
              (origin[i] * u_axis[j] + u_axis[i] * origin[j]) * moment_a +
              (origin[i] * v_axis[j] + v_axis[i] * origin[j]) * moment_b +
              u_axis[i] * u_axis[j] * moment_aa +
              (u_axis[i] * v_axis[j] + v_axis[i] * u_axis[j]) * moment_ab +
              v_axis[i] * v_axis[j] * moment_bb;
        }
      }
      const double trace = second_moment[0][0] + second_moment[1][1] + second_moment[2][2];
      acc.volume += distance * area / 3.0;
      acc.first = acc.first + first_moment * distance;
      for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) {
          const double delta = (i == j) ? trace : 0.0;
          acc.second[i][j] += distance * (delta - second_moment[i][j]) / 5.0;
        }
      }
    }
    loop_id = brep.loop(loop_id).next;
  }
  return {};
}

inline Result<void> mass_add_face(const Brep& brep, FaceId id, const MassPropertiesOptions& options,
                                  const std::vector<double>& nodes,
                                  const std::vector<double>& weights, MassIntegrals& acc) {
  const Face& face = brep.face(id);
  if (face.surface == kNoGeometry || face.surface >= brep.geometry().surfaces.size()) {
    return Error{Errc::index_out_of_range,
                 "mass_properties: face surface index out of range", id.index(), 0u};
  }
  const SurfaceGeometry& surface = brep.geometry().surfaces[face.surface];
  if (const Plane3* plane = std::get_if<Plane3>(&surface)) {
    return mass_add_planar_face(brep, id, *plane, options, acc);
  }
  if (const Cylinder3* cylinder = std::get_if<Cylinder3>(&surface)) {
    const Result<MassUvDomain> domain =
        mass_cylinder_domain(brep, id, *cylinder, options.sampling_deflection);
    if (!domain.has_value()) return domain.error();
    mass_integrate_parametric(*cylinder, domain.value(), face.reversed, nodes, weights, acc);
    return {};
  }
  if (const Cone3* cone = std::get_if<Cone3>(&surface)) {
    const Result<MassUvDomain> domain =
        mass_cone_domain(brep, id, *cone, options.sampling_deflection);
    if (!domain.has_value()) return domain.error();
    mass_integrate_parametric(*cone, domain.value(), face.reversed, nodes, weights, acc);
    return {};
  }
  if (const Sphere3* sphere = std::get_if<Sphere3>(&surface)) {
    const Result<MassUvDomain> domain =
        mass_sphere_domain(brep, id, *sphere, options.sampling_deflection);
    if (!domain.has_value()) return domain.error();
    mass_integrate_parametric(*sphere, domain.value(), face.reversed, nodes, weights, acc);
    return {};
  }
  if (const Torus3* torus = std::get_if<Torus3>(&surface)) {
    const Result<MassUvDomain> domain =
        mass_torus_domain(brep, id, *torus, options.sampling_deflection);
    if (!domain.has_value()) return domain.error();
    mass_integrate_parametric(*torus, domain.value(), face.reversed, nodes, weights, acc);
    return {};
  }
  return Error{Errc::unsupported, "mass_properties: unsupported surface type", id.index(), 0u};
}

}  // namespace detail

inline Result<MassProperties> mass_properties(
    const Brep& brep, const MassPropertiesOptions& options = {}) {
  if (!std::isfinite(options.sampling_deflection) || !(options.sampling_deflection > 0.0) ||
      options.quadrature_order <= 0) {
    return Error{Errc::tolerance_violation, "mass_properties: options are out of range", 0u, 0u};
  }
  if (brep.shell_count() == 0) {
    return detail::mass_numeric_failure("mass_properties: brep has no shells");
  }

  detail::MassIntegrals integrals;
  std::vector<double> nodes;
  std::vector<double> weights;
  detail::mass_gauss_legendre(options.quadrature_order, nodes, weights);
  for (const ShellId shell_id : brep.shell_ids()) {
    if (!brep.shell(shell_id).closed) return detail::mass_not_closed(shell_id);
    for (const FaceId face_id : brep.shell_faces(shell_id)) {
      const Result<void> added =
          detail::mass_add_face(brep, face_id, options, nodes, weights, integrals);
      if (!added.has_value()) return added.error();
    }
  }

  if (!std::isfinite(integrals.volume) || !(integrals.volume > 0.0)) {
    return detail::mass_numeric_failure("mass_properties: volume is not positive");
  }
  const Vec3d centroid = integrals.first * (0.25 / integrals.volume);
  if (!std::isfinite(centroid[0]) || !std::isfinite(centroid[1]) ||
      !std::isfinite(centroid[2])) {
    return detail::mass_numeric_failure("mass_properties: centroid is not finite");
  }

  MassProperties result;
  result.volume = integrals.volume;
  result.centroid = Point3d{centroid[0], centroid[1], centroid[2]};
  const double centroid_squared = dot(centroid, centroid);
  for (int i = 0; i < 3; ++i) {
    for (int j = 0; j < 3; ++j) {
      const double delta = (i == j) ? centroid_squared : 0.0;
      const double shift = delta - centroid[i] * centroid[j];
      result.inertia[i][j] = integrals.second[i][j] - integrals.volume * shift;
      if (!std::isfinite(result.inertia[i][j])) {
        return detail::mass_numeric_failure("mass_properties: inertia is not finite");
      }
    }
  }
  return result;
}

}  // namespace v1
}  // namespace geom
