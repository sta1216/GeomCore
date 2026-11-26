#pragma once

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>
#include <span>
#include <utility>
#include <vector>

#include <geom/brep/assemble.hpp>
#include <geom/core/error.hpp>
#include <geom/core/result.hpp>
#include <geom/primitives/point.hpp>
#include <geom/primitives/vec.hpp>

namespace geom {
inline namespace v1 {
namespace detail {

// One half-space dot(normal, x) <= offset. `normal` need not be unit length:
// the enumeration normalizes it and divides `offset` by the same length, so
// callers may pass plane equations in any convenient scaling.
struct HalfSpacePlane {
  Vec3d normal{};
  double offset = 0.0;
};

// Enumerates the vertices of the convex cell that is the intersection of the
// given half-spaces. Every candidate is the solution of three plane equations,
// kept only when it satisfies all the other planes within `tolerance`, and
// deduplicated. The result is double-precision by design (the open-shell
// operation is allowed to model approximately), so no exact predicates are
// used here. An empty intersection is an empty vector, not an error.
inline Result<std::vector<Point3d>> half_space_vertices(
    std::span<const HalfSpacePlane> planes, double tolerance) {
  if (!std::isfinite(tolerance) || tolerance <= 0.0) {
    return Error{Errc::tolerance_violation,
                 "half_space: tolerance must be finite and positive"};
  }

  const std::size_t count = planes.size();
  std::vector<Vec3d> normals(count);
  std::vector<double> offsets(count);
  for (std::size_t i = 0; i < count; ++i) {
    const Vec3d n = planes[i].normal;
    const double length = std::sqrt(dot(n, n));
    if (!std::isfinite(length) || length <= 0.0) {
      return Error{Errc::degenerate,
                   "half_space: plane normal is zero or non-finite"};
    }
    if (!std::isfinite(planes[i].offset)) {
      return Error{Errc::non_finite, "half_space: plane offset is not finite"};
    }
    normals[i] = Vec3d{n[0] / length, n[1] / length, n[2] / length};
    offsets[i] = planes[i].offset / length;
  }

  std::vector<Point3d> accepted;
  for (std::size_t i = 0; i < count; ++i) {
    for (std::size_t j = i + 1; j < count; ++j) {
      for (std::size_t k = j + 1; k < count; ++k) {
        // Cramer's rule: p solves n_i . p = d_i for the three planes, with the
        // determinant written in mixed-product form.
        const Vec3d c0 = cross(normals[j], normals[k]);
        const Vec3d c1 = cross(normals[k], normals[i]);
        const Vec3d c2 = cross(normals[i], normals[j]);
        const double det = dot(normals[i], c0);
        // Exactly parallel planes have no single intersection point. Skip
        // rather than fail, since other triples may describe the same vertex.
        if (!std::isfinite(det) || det == 0.0) continue;
        const Vec3d p{
            (offsets[i] * c0[0] + offsets[j] * c1[0] + offsets[k] * c2[0]) / det,
            (offsets[i] * c0[1] + offsets[j] * c1[1] + offsets[k] * c2[1]) / det,
            (offsets[i] * c0[2] + offsets[j] * c1[2] + offsets[k] * c2[2]) / det};
        if (!std::isfinite(p[0]) || !std::isfinite(p[1]) ||
            !std::isfinite(p[2])) {
          continue;
        }

        // Relative slack for this candidate, with no absolute floor: it scales
        // with the triple's offsets and the solution's own magnitude so that
        // tiny and huge models keep the same relative behaviour.
        const double point_scale =
            std::max({std::fabs(offsets[i]), std::fabs(offsets[j]),
                      std::fabs(offsets[k]), std::fabs(p[0]),
                      std::fabs(p[1]), std::fabs(p[2])});
        const double slack = std::max(tolerance, 1e-12 * point_scale);

        // Conditioning check: if the intersection is so ill-conditioned that
        // the rounding error bound of the Cramer solve (about
        // eps * point_scale / |det|) exceeds the slack, the recovered point
        // says nothing reliable. Skip the triple instead of admitting a
        // spurious vertex.
        if (32.0 * std::numeric_limits<double>::epsilon() * point_scale /
                std::fabs(det) >
            slack) {
          continue;
        }

        bool inside = true;
        for (std::size_t l = 0; l < count && inside; ++l) {
          inside = dot(normals[l], p) <= offsets[l] + slack;
        }
        if (!inside) continue;

        // Duplicates arise when several plane triples meet at the same vertex.
        // Keep the first representative so the output order is deterministic
        // for a fixed plane order.
        bool duplicate = false;
        for (const Point3d& q : accepted) {
          duplicate = duplicate ||
                      (std::fabs(p[0] - q[0]) <= slack &&
                       std::fabs(p[1] - q[1]) <= slack &&
                       std::fabs(p[2] - q[2]) <= slack);
        }
        if (!duplicate) accepted.push_back(Point3d{p[0], p[1], p[2]});
      }
    }
  }
  return accepted;
}

// Extracts one outward-oriented boundary face per plane that supports at least
// three of `vertices`. `vertices` is expected to be the output of
// half_space_vertices for the same planes; a plane with fewer than three
// incident vertices does not bound the cell and is skipped, which is how
// redundant planes disappear. The returned ring is an open ring: the last
// point is not repeated.
//
// Coincident or duplicate planes are not detected: every plane that supports at
// least three vertices contributes its own face, so two equal planes produce two
// equal faces and the result cannot be assembled into a solid. Callers must
// deduplicate coincident planes before calling.
//
// A vertex belongs to plane i when `|dot(n_i, p) - offset_i| <= slack`, with the
// per-(plane, vertex) slack `max(tolerance, 1e-12 * max(|offset_i|, |p|
// components))`. The relative term scales with that plane's own equation and
// that point's coordinates, so no global model scale and no absolute floor are
// needed: a far redundant plane cannot capture the cell's vertices, and huge or
// tiny cells keep the same relative behaviour.
//
// Each ring is sorted by the polar angle of its points around their centroid in
// an orthonormal in-plane basis (u, v), with v = cross(n, u). Because (u, v, n)
// is right-handed, that order is counter-clockwise about the plane's unit
// normal, so the returned `normal` is outward and the ring winds
// counter-clockwise about it, as assemble_solid requires. Ties are broken by
// the incident order, so the rings are deterministic for a fixed plane order.
inline Result<std::vector<AssembledFace>> half_space_faces(
    std::span<const HalfSpacePlane> planes, std::span<const Point3d> vertices,
    double tolerance) {
  if (!std::isfinite(tolerance) || tolerance <= 0.0) {
    return Error{Errc::tolerance_violation,
                 "half_space: tolerance must be finite and positive"};
  }

  const std::size_t count = planes.size();
  std::vector<Vec3d> normals(count);
  std::vector<double> offsets(count);
  for (std::size_t i = 0; i < count; ++i) {
    const Vec3d n = planes[i].normal;
    const double length = std::sqrt(dot(n, n));
    if (!std::isfinite(length) || length <= 0.0) {
      return Error{Errc::degenerate,
                   "half_space: plane normal is zero or non-finite"};
    }
    if (!std::isfinite(planes[i].offset)) {
      return Error{Errc::non_finite, "half_space: plane offset is not finite"};
    }
    normals[i] = Vec3d{n[0] / length, n[1] / length, n[2] / length};
    offsets[i] = planes[i].offset / length;
  }

  struct incident_point {
    double angle = 0.0;
    std::size_t index = 0;
  };

  std::vector<AssembledFace> faces;
  faces.reserve(count);
  std::vector<std::size_t> incident;
  std::vector<incident_point> ordered;
  for (std::size_t i = 0; i < count; ++i) {
    incident.clear();
    for (std::size_t v = 0; v < vertices.size(); ++v) {
      const Point3d& p = vertices[v];
      const double plane_scale =
          std::max({std::fabs(offsets[i]), std::fabs(p[0]), std::fabs(p[1]),
                    std::fabs(p[2])});
      const double slack = std::max(tolerance, 1e-12 * plane_scale);
      const Vec3d point{p[0], p[1], p[2]};
      if (std::fabs(dot(normals[i], point) - offsets[i]) <= slack) {
        incident.push_back(v);
      }
    }
    // A plane touching fewer than three vertices cannot bound a polygon: it is
    // a redundant half-space or a plane that only grazes an edge or a corner.
    if (incident.size() < 3u) continue;

    // In-plane orthonormal basis. The axis with the largest component is the
    // most parallel to the normal and the least reliable, so the first axis
    // leaving a residual longer than 1/2 is taken; because a unit normal cannot
    // have two components of magnitude at least sqrt(3)/2, at least one axis
    // always qualifies. The order +x, +y, +z keeps the basis deterministic.
    Vec3d u{};
    bool have_u = false;
    for (int axis = 0; axis < 3 && !have_u; ++axis) {
      Vec3d e{};
      e[axis] = 1.0;
      const Vec3d residual = e - normals[i] * dot(normals[i], e);
      const double residual_length = std::sqrt(dot(residual, residual));
      if (residual_length > 0.5) {
        u = Vec3d{residual[0] / residual_length, residual[1] / residual_length,
                  residual[2] / residual_length};
        have_u = true;
      }
    }
    if (!have_u) continue;
    const Vec3d v_axis = cross(normals[i], u);

    Vec3d centroid{};
    for (const std::size_t index : incident) {
      centroid = centroid + Vec3d{vertices[index][0], vertices[index][1],
                                  vertices[index][2]};
    }
    const double inverse_count = 1.0 / static_cast<double>(incident.size());
    centroid = Vec3d{centroid[0] * inverse_count, centroid[1] * inverse_count,
                     centroid[2] * inverse_count};

    ordered.clear();
    ordered.reserve(incident.size());
    for (const std::size_t index : incident) {
      const Point3d& p = vertices[index];
      const Vec3d delta{p[0] - centroid[0], p[1] - centroid[1],
                        p[2] - centroid[2]};
      ordered.push_back(incident_point{
          std::atan2(dot(delta, v_axis), dot(delta, u)), index});
    }
    std::sort(ordered.begin(), ordered.end(),
              [](const incident_point& a, const incident_point& b) {
                if (a.angle != b.angle) return a.angle < b.angle;
                return a.index < b.index;
              });

    AssembledFace face;
    face.outer.reserve(ordered.size());
    for (const incident_point& entry : ordered) {
      face.outer.push_back(vertices[entry.index]);
    }
    face.normal = normals[i];
    faces.push_back(std::move(face));
  }

  if (faces.empty()) {
    return Error{Errc::degenerate,
                 "half_space: fewer than three planes bound a cell"};
  }
  return faces;
}

// Assembles the convex cell defined by `planes` into a closed B-rep solid.
// Enumerates the cell's vertices first, so an empty intersection is reported
// as degenerate before any assembly is attempted; then extracts one outward
// ring per bounding plane and hands them to assemble_solid. An assembly failure
// (for example coincident duplicate planes producing a non-manifold edge) is
// returned unchanged.
inline Result<Brep> half_space_solid(std::span<const HalfSpacePlane> planes,
                                     double tolerance) {
  const Result<std::vector<Point3d>> vertices =
      half_space_vertices(planes, tolerance);
  if (!vertices) return vertices.error();
  if (vertices.value().empty()) {
    return Error{Errc::degenerate,
                 "half_space: planes have no common point"};
  }
  const Result<std::vector<AssembledFace>> faces =
      half_space_faces(planes, vertices.value(), tolerance);
  if (!faces) return faces.error();
  return assemble_solid(faces.value(), tolerance);
}

}  // namespace detail
}  // namespace v1
}  // namespace geom
