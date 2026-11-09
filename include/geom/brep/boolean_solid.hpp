#pragma once

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <span>
#include <utility>
#include <vector>

#include <geom/brep/assemble.hpp>
#include <geom/brep/face_intersect.hpp>
#include <geom/brep/face_query.hpp>
#include <geom/brep/inside.hpp>
#include <geom/brep/plane_frame.hpp>
#include <geom/brep/planar_cut.hpp>
#include <geom/brep/topology.hpp>
#include <geom/brep/validate.hpp>
#include <geom/core/error.hpp>
#include <geom/core/result.hpp>
#include <geom/primitives/point.hpp>
#include <geom/primitives/vec.hpp>

namespace geom {
inline namespace v1 {

enum class SolidBooleanOp {
  union_,
  intersection,
  difference,
  symmetric_difference
};

struct SolidBooleanOptions {
  double tolerance = 1e-9;
};

struct SolidBooleanReport {
  std::uint32_t face_pairs = 0;
  std::uint32_t cut_segments = 0;
  std::uint32_t fragments = 0;
  std::uint32_t kept = 0;
  std::uint32_t dropped = 0;
  std::uint32_t flipped = 0;
  double max_cut_residual = 0.0;
};

struct SolidBooleanResult {
  Brep brep;
  SolidBooleanReport report;
};

namespace detail {

// Membership rule of the regularized Boolean operation.
inline bool solid_boolean_member(bool in_a, bool in_b,
                                 SolidBooleanOp op) noexcept {
  switch (op) {
    case SolidBooleanOp::union_:
      return in_a || in_b;
    case SolidBooleanOp::intersection:
      return in_a && in_b;
    case SolidBooleanOp::difference:
      return in_a && !in_b;
    case SolidBooleanOp::symmetric_difference:
      return in_a != in_b;
  }
  return false;
}

struct solid_boolean_bounds {
  Point3d lo{};
  Point3d hi{};
  bool valid = false;
};

inline void solid_boolean_extend(solid_boolean_bounds& bounds,
                                 const Point3d& point) noexcept {
  if (!bounds.valid) {
    bounds.lo = point;
    bounds.hi = point;
    bounds.valid = true;
    return;
  }
  for (int axis = 0; axis < 3; ++axis) {
    bounds.lo[axis] = std::min(bounds.lo[axis], point[axis]);
    bounds.hi[axis] = std::max(bounds.hi[axis], point[axis]);
  }
}

inline bool solid_boolean_overlap(const solid_boolean_bounds& first,
                                  const solid_boolean_bounds& second,
                                  double tol) noexcept {
  if (!first.valid || !second.valid) return false;
  for (int axis = 0; axis < 3; ++axis) {
    if (first.hi[axis] + tol < second.lo[axis]) return false;
    if (second.hi[axis] + tol < first.lo[axis]) return false;
  }
  return true;
}

inline bool solid_boolean_finite_point(const Point3d& point) noexcept {
  return std::isfinite(point[0]) && std::isfinite(point[1]) &&
         std::isfinite(point[2]);
}

inline double solid_boolean_dot_point(const Vec3d& normal,
                                     const Point3d& point) noexcept {
  return normal[0] * point[0] + normal[1] * point[1] + normal[2] * point[2];
}

inline Error solid_boolean_topology_error(const ValidationReport& report) {
  for (const ValidationIssue& issue : report.issues) {
    if (issue.code == Errc::unsupported) {
      return Error{Errc::unsupported, issue.message, issue.entity, 0u};
    }
  }
  const ValidationIssue& issue = report.issues.front();
  return Error{issue.code, issue.message, issue.entity, 0u};
}

inline Result<void> solid_boolean_check_operand(const Brep& brep) {
  const ValidationReport topology = validate_topology(brep);
  if (!topology.ok()) return solid_boolean_topology_error(topology);
  const ValidationReport tolerances = validate_tolerances(brep);
  if (!tolerances.ok()) {
    const ValidationIssue& issue = tolerances.issues.front();
    return Error{Errc::tolerance_violation, issue.message, issue.entity, 0u};
  }
  return {};
}

struct solid_boolean_face {
  FaceId id{};
  Plane3 plane{};
  PlaneFrame frame{};
  solid_boolean_bounds bounds{};
  std::vector<Segment3> cuts;
};

// Collects the shell faces of one operand together with the cuts that the
// other operand imposes on them.
inline Result<std::vector<solid_boolean_face>> solid_boolean_build_faces(
    const Brep& brep, std::span<const FaceId> face_ids) {
  std::vector<solid_boolean_face> faces;
  faces.reserve(face_ids.size());
  for (const FaceId id : face_ids) {
    const Result<Plane3> plane = face_plane(brep, id);
    if (!plane) {
      if (plane.error().code == Errc::unsupported) {
        return Error{Errc::unsupported,
                     "boolean_solid: curved faces are deferred to M3b",
                     id.index(), 0u};
      }
      return plane.error();
    }
    const Result<PlaneFrame> frame = make_plane_frame(plane.value());
    if (!frame) return frame.error();
    const Result<std::vector<std::vector<Point2d>>> loops =
        face_loops_2d(brep, id, frame.value());
    if (!loops) return loops.error();
    solid_boolean_face face;
    face.id = id;
    face.plane = plane.value();
    face.frame = frame.value();
    for (const std::vector<Point2d>& loop : loops.value()) {
      for (const Point2d& point : loop) {
        solid_boolean_extend(face.bounds,
                             plane_frame_to_3d(face.frame, point));
      }
    }
    faces.push_back(std::move(face));
  }
  return faces;
}

inline bool solid_boolean_same_ring(const std::vector<Point3d>& first,
                                    const std::vector<Point3d>& second,
                                    double tol) noexcept {
  if (first.size() != second.size()) return false;
  if (first.empty()) return true;
  const std::size_t count = first.size();
  for (std::size_t shift = 0; shift < count; ++shift) {
    bool matches = true;
    for (std::size_t i = 0; i < count && matches; ++i) {
      const Point3d& other = second[(shift + i) % count];
      for (int axis = 0; axis < 3; ++axis) {
        if (std::fabs(first[i][axis] - other[axis]) > tol) {
          matches = false;
          break;
        }
      }
    }
    if (matches) return true;
  }
  return false;
}

// Two fragments are coincident when their normals agree, their planes agree
// within tol, and their ring multisets agree within tol up to rotation.
inline bool solid_boolean_same_fragment(const AssembledFace& first,
                                        const AssembledFace& second,
                                        double tol) noexcept {
  for (int axis = 0; axis < 3; ++axis) {
    if (std::fabs(first.normal[axis] - second.normal[axis]) > tol) {
      return false;
    }
  }
  if (first.outer.empty() || second.outer.empty()) return false;
  const double first_offset = solid_boolean_dot_point(first.normal, first.outer.front());
  const double second_offset = solid_boolean_dot_point(second.normal, second.outer.front());
  if (std::fabs(first_offset - second_offset) > tol) return false;
  if (!solid_boolean_same_ring(first.outer, second.outer, tol)) return false;
  if (first.holes.size() != second.holes.size()) return false;
  std::vector<bool> used(second.holes.size(), false);
  for (const std::vector<Point3d>& hole : first.holes) {
    bool found = false;
    for (std::size_t i = 0; i < second.holes.size(); ++i) {
      if (!used[i] && solid_boolean_same_ring(hole, second.holes[i], tol)) {
        used[i] = true;
        found = true;
        break;
      }
    }
    if (!found) return false;
  }
  return true;
}

inline bool solid_boolean_fragment_less(const AssembledFace& first,
                                        const AssembledFace& second) noexcept {
  for (int axis = 0; axis < 3; ++axis) {
    if (first.normal[axis] != second.normal[axis]) {
      return first.normal[axis] < second.normal[axis];
    }
  }
  const double first_offset = solid_boolean_dot_point(first.normal, first.outer.front());
  const double second_offset = solid_boolean_dot_point(second.normal, second.outer.front());
  if (first_offset != second_offset) return first_offset < second_offset;
  const std::size_t count = std::min(first.outer.size(), second.outer.size());
  for (std::size_t i = 0; i < count; ++i) {
    for (int axis = 0; axis < 3; ++axis) {
      if (first.outer[i][axis] != second.outer[i][axis]) {
        return first.outer[i][axis] < second.outer[i][axis];
      }
    }
  }
  return first.outer.size() < second.outer.size();
}

// Keeps the first fragment of every coincident group in deterministic order.
inline std::vector<AssembledFace> solid_boolean_deduplicate(
    std::vector<AssembledFace> faces, double tol) {
  std::sort(faces.begin(), faces.end(), solid_boolean_fragment_less);
  std::vector<AssembledFace> unique;
  unique.reserve(faces.size());
  for (AssembledFace& face : faces) {
    bool duplicate = false;
    for (const AssembledFace& other : unique) {
      if (solid_boolean_same_fragment(other, face, tol)) {
        duplicate = true;
        break;
      }
    }
    if (!duplicate) unique.push_back(std::move(face));
  }
  return unique;
}

inline void solid_boolean_reverse_rings(std::vector<Point3d>& outer,
                                       std::vector<std::vector<Point3d>>& holes) {
  std::reverse(outer.begin(), outer.end());
  for (std::vector<Point3d>& hole : holes) {
    std::reverse(hole.begin(), hole.end());
  }
}

// Cuts one face, classifies every fragment by the operation, and appends the
// surviving outward-oriented fragments to `kept`.
inline Result<void> solid_boolean_classify_face(
    const Brep& owner, const solid_boolean_face& face, const Brep& a,
    SolidId solid_a, const Brep& b, SolidId solid_b, SolidBooleanOp op,
    double tol, double delta, SolidBooleanReport& report,
    std::vector<AssembledFace>& kept) {
  const Result<std::vector<FaceFragment2d>> fragments =
      cut_face_by_segments(owner, face.id, face.cuts, tol);
  if (!fragments) return fragments.error();

  Vec3d outward = face.frame.normal;
  const bool reversed = owner.face(face.id).reversed;
  if (reversed) {
    outward = Vec3d{-outward[0], -outward[1], -outward[2]};
  }

  for (const FaceFragment2d& fragment : fragments.value()) {
    ++report.fragments;
    const Point3d sample = plane_frame_to_3d(face.frame, fragment.sample);
    if (!solid_boolean_finite_point(sample)) {
      return Error{Errc::non_finite, "boolean_solid: fragment sample is not finite",
                   face.id.index(), 0u};
    }

    std::vector<Point3d> outer;
    outer.reserve(fragment.outer.size());
    for (const Point2d& point : fragment.outer) {
      outer.push_back(plane_frame_to_3d(face.frame, point));
    }
    std::vector<std::vector<Point3d>> holes;
    holes.reserve(fragment.holes.size());
    for (const std::vector<Point2d>& hole : fragment.holes) {
      std::vector<Point3d> ring;
      ring.reserve(hole.size());
      for (const Point2d& point : hole) {
        ring.push_back(plane_frame_to_3d(face.frame, point));
      }
      holes.push_back(std::move(ring));
    }

    Vec3d normal = outward;
    if (reversed) solid_boolean_reverse_rings(outer, holes);

    const Point3d plus{sample[0] + normal[0] * delta,
                       sample[1] + normal[1] * delta,
                       sample[2] + normal[2] * delta};
    const Point3d minus{sample[0] - normal[0] * delta,
                        sample[1] - normal[1] * delta,
                        sample[2] - normal[2] * delta};

    const Result<bool> in_a_plus = point_in_solid(a, solid_a, plus, tol);
    if (!in_a_plus) return in_a_plus.error();
    const Result<bool> in_b_plus = point_in_solid(b, solid_b, plus, tol);
    if (!in_b_plus) return in_b_plus.error();
    const Result<bool> in_a_minus = point_in_solid(a, solid_a, minus, tol);
    if (!in_a_minus) return in_a_minus.error();
    const Result<bool> in_b_minus = point_in_solid(b, solid_b, minus, tol);
    if (!in_b_minus) return in_b_minus.error();

    const bool keep_plus =
        solid_boolean_member(in_a_plus.value(), in_b_plus.value(), op);
    const bool keep_minus =
        solid_boolean_member(in_a_minus.value(), in_b_minus.value(), op);
    if (keep_plus == keep_minus) {
      ++report.dropped;
      continue;
    }
    if (keep_plus) {
      solid_boolean_reverse_rings(outer, holes);
      normal = Vec3d{-normal[0], -normal[1], -normal[2]};
      ++report.flipped;
    }
    kept.push_back(
        AssembledFace{std::move(outer), std::move(holes), normal});
    ++report.kept;
  }
  return {};
}

}  // namespace detail

// Computes the regularized planar Boolean of two single-shell solids. Both
// operands must be valid, single-solid, single-shell Breps whose faces are all
// planar; curved operands are rejected as unsupported. The result is assembled
// with `assemble_solid`, so a result that cannot form a closed manifold solid is
// reported as an error instead of being returned malformed.
inline Result<SolidBooleanResult> boolean_solid(
    const Brep& a, const Brep& b, SolidBooleanOp op,
    const SolidBooleanOptions& options = {}) {
  const double tol = options.tolerance;
  if (!std::isfinite(tol)) {
    return Error{Errc::non_finite,
                 "boolean_solid: tolerance must be positive and finite"};
  }
  if (!(tol > 0.0)) {
    return Error{Errc::degenerate,
                 "boolean_solid: tolerance must be positive and finite"};
  }

  const Result<void> check_a = detail::solid_boolean_check_operand(a);
  if (!check_a) return check_a.error();
  const Result<void> check_b = detail::solid_boolean_check_operand(b);
  if (!check_b) return check_b.error();

  const std::vector<SolidId> solids_a = a.solid_ids();
  const std::vector<SolidId> solids_b = b.solid_ids();
  if (solids_a.size() != 1u || solids_b.size() != 1u) {
    return Error{Errc::unsupported,
                 "boolean_solid: expected exactly one solid per operand"};
  }
  const std::span<const ShellId> shells_a = a.solid_shells(solids_a.front());
  const std::span<const ShellId> shells_b = b.solid_shells(solids_b.front());
  if (shells_a.size() != 1u || shells_b.size() != 1u) {
    return Error{Errc::unsupported,
                 "boolean_solid: multi-shell operands are deferred to M3c"};
  }
  const ShellId shell_a = shells_a.front();
  const ShellId shell_b = shells_b.front();
  if (!a.shell(shell_a).closed || !b.shell(shell_b).closed) {
    return Error{Errc::unsupported,
                 "boolean_solid: operand shell must be closed"};
  }

  Result<std::vector<detail::solid_boolean_face>> faces_a =
      detail::solid_boolean_build_faces(a, a.shell_faces(shell_a));
  if (!faces_a) return faces_a.error();
  Result<std::vector<detail::solid_boolean_face>> faces_b =
      detail::solid_boolean_build_faces(b, b.shell_faces(shell_b));
  if (!faces_b) return faces_b.error();

  SolidBooleanReport report;
  for (detail::solid_boolean_face& face : faces_a.value()) {
    for (detail::solid_boolean_face& other : faces_b.value()) {
      if (!detail::solid_boolean_overlap(face.bounds, other.bounds, tol)) {
        continue;
      }
      ++report.face_pairs;
      const Result<FacePairCut> cut = intersect_planar_faces(
          a, face.id, b, other.id, tol);
      if (!cut) return cut.error();
      for (const Segment3& segment : cut.value().segments) {
        if (!detail::solid_boolean_finite_point(segment.a) ||
            !detail::solid_boolean_finite_point(segment.b)) {
          return Error{Errc::numeric_failure,
                       "boolean_solid: cut segment is not finite"};
        }
        for (const Point3d& point : {segment.a, segment.b}) {
          const double first_residual =
              std::fabs(dot(face.frame.normal, point - face.frame.origin));
          const double second_residual =
              std::fabs(dot(other.frame.normal, point - other.frame.origin));
          if (!std::isfinite(first_residual) ||
              !std::isfinite(second_residual)) {
            return Error{Errc::numeric_failure,
                         "boolean_solid: cut residual is not finite"};
          }
          report.max_cut_residual =
              std::max(report.max_cut_residual,
                       std::max(first_residual, second_residual));
        }
        face.cuts.push_back(segment);
        other.cuts.push_back(segment);
        ++report.cut_segments;
      }
    }
  }

  detail::solid_boolean_bounds bounds;
  for (const detail::solid_boolean_face& face : faces_a.value()) {
    if (face.bounds.valid) {
      detail::solid_boolean_extend(bounds, face.bounds.lo);
      detail::solid_boolean_extend(bounds, face.bounds.hi);
    }
  }
  for (const detail::solid_boolean_face& face : faces_b.value()) {
    if (face.bounds.valid) {
      detail::solid_boolean_extend(bounds, face.bounds.lo);
      detail::solid_boolean_extend(bounds, face.bounds.hi);
    }
  }
  if (!bounds.valid) {
    return Error{Errc::degenerate, "boolean_solid: operands have no geometry"};
  }
  const Vec3d span = bounds.hi - bounds.lo;
  const double diagonal = std::sqrt(dot(span, span));
  const double delta = std::max(64.0 * tol, 1e-7 * diagonal);

  std::vector<AssembledFace> kept;
  for (const detail::solid_boolean_face& face : faces_a.value()) {
    const Result<void> status = detail::solid_boolean_classify_face(
        a, face, a, solids_a.front(), b, solids_b.front(), op, tol, delta,
        report, kept);
    if (!status) return status.error();
  }
  for (const detail::solid_boolean_face& face : faces_b.value()) {
    const Result<void> status = detail::solid_boolean_classify_face(
        b, face, a, solids_a.front(), b, solids_b.front(), op, tol, delta,
        report, kept);
    if (!status) return status.error();
  }

  std::vector<AssembledFace> unique =
      detail::solid_boolean_deduplicate(std::move(kept), tol);
  if (unique.empty()) {
    return Error{Errc::degenerate, "boolean_solid: result is empty"};
  }
  Result<Brep> assembled = assemble_solid(unique, tol);
  if (!assembled) return assembled.error();
  return SolidBooleanResult{std::move(assembled).value(), report};
}

}  // namespace v1
}  // namespace geom
