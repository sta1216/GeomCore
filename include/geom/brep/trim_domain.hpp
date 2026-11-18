#pragma once

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <span>
#include <utility>
#include <vector>

#include <geom/brep/face_uv.hpp>
#include <geom/brep/topology.hpp>
#include <geom/core/error.hpp>
#include <geom/core/result.hpp>
#include <geom/primitives/point.hpp>
#include <geom/surfaces/analytic.hpp>

namespace geom {
inline namespace v1 {

// One curve along which a curved face is split. The polyline lives in the
// surface parameter domain reported by face_loops_uv() for the same face. A
// closed cut repeats its first point as its last point.
struct FaceCut {
  std::vector<Point2d> uv;
};

namespace detail {

// How one side of a face's uv rectangle behaves when the face is split. A seam
// copy is the repeated parameter edge of a periodic axis: both copies name the
// same physical curve, so fragments along it are dropped. A collapsed side is a
// whole parameter edge that maps to one physical point (a sphere pole), so a
// run along it collapses to that point. A real boundary is the face's own
// stored edge and stays as geometry.
enum class TrimSideKind { kSeamCopy, kCollapsed, kRealBoundary };

// The full uv rectangle a synthetic-domain splitter works on. `bounds` is
// {u_lo, u_hi, v_lo, v_hi}; `sides` classifies the four sides in the same
// order.
struct TrimFaceDomain {
  std::array<double, 4> bounds{};
  bool periodic_u = false;
  bool periodic_v = false;
  std::array<TrimSideKind, 4> sides{};

  bool collapsed_v() const noexcept {
    return sides[2] == TrimSideKind::kCollapsed ||
           sides[3] == TrimSideKind::kCollapsed;
  }
};

// The analytic uv rectangle of a surface whose face has no real boundary loop.
// Only the synthetic-domain surfaces of this milestone are described: a sphere
// is periodic in u and collapses at both poles, a torus is periodic on both
// axes, and a cone is periodic in u with a real base edge at v = 0 and a
// collapsed apex at v = radius / tan(half_angle). Every other surface is not a
// synthetic domain yet.
//
// Precondition: the caller has already established that the face is a
// synthetic-domain face (its stored loop is the seam copy or the full
// rectangle). This function only inspects the surface: it does not look at the
// face loops, so a genuinely trimmed sphere would receive the full rectangle
// as well. The projection and tolerance parameters are reserved for that
// loop-level validation in the caller.
inline Result<TrimFaceDomain> trim_face_domain(const Brep& brep, FaceId face,
                                               const FaceUvProjection& projection,
                                               double tolerance) {
  (void)projection;
  (void)tolerance;
  if (face != FaceId{}) {
    const Face& stored = brep.face(face);
    if (stored.surface != kNoGeometry) {
      const SurfaceGeometry& surface = brep.geometry().surfaces[stored.surface];
      TrimFaceDomain domain;
      if (std::holds_alternative<Sphere3>(surface)) {
        domain.bounds = {0.0, 2.0 * kFaceUvPi, -0.5 * kFaceUvPi,
                         0.5 * kFaceUvPi};
        domain.periodic_u = true;
        domain.periodic_v = false;
        domain.sides = {TrimSideKind::kSeamCopy, TrimSideKind::kSeamCopy,
                        TrimSideKind::kCollapsed, TrimSideKind::kCollapsed};
        return domain;
      }
      if (std::holds_alternative<Torus3>(surface)) {
        domain.bounds = {0.0, 2.0 * kFaceUvPi, 0.0, 2.0 * kFaceUvPi};
        domain.periodic_u = true;
        domain.periodic_v = true;
        domain.sides = {TrimSideKind::kSeamCopy, TrimSideKind::kSeamCopy,
                        TrimSideKind::kSeamCopy, TrimSideKind::kSeamCopy};
        return domain;
      }
      if (const auto* cone = std::get_if<Cone3>(&surface)) {
        domain.bounds = {0.0, 2.0 * kFaceUvPi, 0.0,
                         cone->radius / std::tan(cone->half_angle)};
        domain.periodic_u = true;
        domain.periodic_v = false;
        domain.sides = {TrimSideKind::kSeamCopy, TrimSideKind::kSeamCopy,
                        TrimSideKind::kRealBoundary,
                        TrimSideKind::kCollapsed};
        return domain;
      }
    }
  }
  return Error{Errc::unsupported,
               "trim_face: no synthetic domain for this surface"};
}

// Builds the descriptor for a caller that already knows the parameter
// rectangle (the boolean driver stores it in its face info). Every periodic
// axis gets seam copies on both of its sides; every other side is a real
// boundary. A collapsed side is only implied by a known surface, so callers
// that have one use trim_face_domain() instead.
inline TrimFaceDomain trim_face_domain_from_ranges(
    const FaceUvProjection& projection, const std::array<double, 2>& u,
    const std::array<double, 2>& v) noexcept {
  TrimFaceDomain domain;
  domain.bounds = {u[0], u[1], v[0], v[1]};
  domain.periodic_u = projection.periodic_u();
  domain.periodic_v = projection.periodic_v();
  domain.sides = {TrimSideKind::kRealBoundary, TrimSideKind::kRealBoundary,
                  TrimSideKind::kRealBoundary, TrimSideKind::kRealBoundary};
  if (domain.periodic_u) {
    domain.sides[0] = TrimSideKind::kSeamCopy;
    domain.sides[1] = TrimSideKind::kSeamCopy;
  }
  if (domain.periodic_v) {
    domain.sides[2] = TrimSideKind::kSeamCopy;
    domain.sides[3] = TrimSideKind::kSeamCopy;
  }
  return domain;
}

// Wraps a periodic parameter into [base, base + period).
inline double trim_face_wrap_periodic(double value, double base,
                                      double period) noexcept {
  const double offset = value - base;
  const double wrapped = offset - period * std::floor(offset / period);
  return base + wrapped;
}

// True when the value lies within the standard relative parameter bound of one
// domain side value.
inline bool trim_face_near_side(double value, double side,
                                double tolerance = 1e-12) noexcept {
  return std::fabs(value - side) <=
         std::max(tolerance, 1e-12 * std::max({1.0, std::fabs(value),
                                               std::fabs(side)}));
}

// True when either endpoint of the two parameter copies names the same physical
// point within the relative bound.
inline bool trim_face_same_xy_rel(const Point2d& a, const Point2d& b,
                                  double tolerance = 1e-12) noexcept {
  const double bound_x =
      std::max(tolerance, 1e-12 * std::max({1.0, std::fabs(a[0]),
                                            std::fabs(b[0])}));
  const double bound_y =
      std::max(tolerance, 1e-12 * std::max({1.0, std::fabs(a[1]),
                                            std::fabs(b[1])}));
  return std::fabs(a[0] - b[0]) <= bound_x &&
         std::fabs(a[1] - b[1]) <= bound_y;
}

// True when the point sits on any of the four domain sides.
inline bool trim_face_point_on_domain_side(const Point2d& point,
                                           const TrimFaceDomain& domain,
                                           double tolerance = 1e-12) noexcept {
  return trim_face_near_side(point[0], domain.bounds[0], tolerance) ||
         trim_face_near_side(point[0], domain.bounds[1], tolerance) ||
         trim_face_near_side(point[1], domain.bounds[2], tolerance) ||
         trim_face_near_side(point[1], domain.bounds[3], tolerance);
}

// Drops the points that only duplicate their predecessor: a segment shorter
// than the tolerance carries no cutting information.
inline void trim_face_dedupe_chain(std::vector<Point2d>& points,
                                   double tolerance) {
  std::vector<Point2d> unique;
  unique.reserve(points.size());
  for (const Point2d& point : points) {
    if (!unique.empty()) {
      const double du = point[0] - unique.back()[0];
      const double dv = point[1] - unique.back()[1];
      if (du * du + dv * dv <= tolerance * tolerance) continue;
    }
    unique.push_back(point);
  }
  points = std::move(unique);
}

// Unwraps one parameter axis of a cut chain so consecutive samples differ by
// less than half a period; a non-periodic axis is only clamped to the face.
inline void trim_face_unwrap_axis(std::vector<Point2d>& points, int axis,
                                  bool periodic,
                                  const std::array<double, 2>& range) {
  if (points.empty()) return;
  if (!periodic) {
    for (Point2d& point : points) {
      point[axis] = std::min(std::max(point[axis], range[0]), range[1]);
    }
    return;
  }
  const double period = range[1] - range[0];
  if (!(period > 0.0)) return;
  const double half = 0.5 * period;
  for (std::size_t i = 1; i < points.size(); ++i) {
    const double delta = points[i][axis] - points[i - 1][axis];
    const double wrapped = delta - period * std::floor((delta + half) / period);
    points[i][axis] = points[i - 1][axis] + wrapped;
  }
}

// The traversal direction of a chain along one axis, judged from its short
// steps: +1 when the chain mostly increases, -1 when it mostly decreases. Steps
// longer than half a period are the seam crossings themselves and are ignored,
// so a chain dominated by one direction is not mistaken for a wrap. When every
// step is long (a two-point chain), assume a forward traversal; the caller can
// still disambiguate by ordering the points.
inline int trim_face_dominant_step_sign(const std::vector<Point2d>& points,
                                        int axis, double half) noexcept {
  std::size_t forward = 0;
  std::size_t backward = 0;
  for (std::size_t i = 1; i < points.size(); ++i) {
    const double delta = points[i][axis] - points[i - 1][axis];
    if (std::fabs(delta) > half) continue;
    if (delta > 0.0) {
      ++forward;
    } else if (delta < 0.0) {
      ++backward;
    }
  }
  return backward > forward ? -1 : 1;
}

// Splits a chain whose values on `axis` are wrapped into [base, base + period)
// where it steps across the periodic seam, inserting the seam point itself so
// both new ends land exactly on the seam edge they cross. A step from the high
// edge down to the low one (a forward traversal) is always split. A step from
// low to high is only split when `split_backward` is set: without it, the step
// is treated as part of the polyline inside the window, which is what the
// synthetic-domain splitter expects from a cut that is already windowed.
// When `dominant_sign` is +1 or -1 the split is direction-aware: only steps
// that wrap against the chain's dominant traversal direction are split, so a
// long step that stays inside the window is not mistaken for a seam crossing.
inline std::vector<std::vector<Point2d>> trim_face_split_at_seam(
    const std::vector<Point2d>& points, int axis, double base, double high,
    double tolerance, bool split_backward, int dominant_sign = 0) {
  const double period = high - base;
  const double half = 0.5 * period;
  const int other = 1 - axis;
  std::vector<std::vector<Point2d>> chains;
  if (points.size() < 2) return chains;
  std::vector<Point2d> current{points.front()};
  const auto finish = [&chains, &current, tolerance]() {
    if (current.size() < 2) {
      current.clear();
      return;
    }
    std::vector<Point2d> chain;
    chain.reserve(current.size());
    for (const Point2d& point : current) {
      if (!chain.empty()) {
        const double du = point[0] - chain.back()[0];
        const double dv = point[1] - chain.back()[1];
        if (du * du + dv * dv <= tolerance * tolerance) continue;
      }
      chain.push_back(point);
    }
    if (chain.size() >= 2) chains.push_back(std::move(chain));
    current.clear();
  };
  const auto seam_point = [axis, other](double along, double across) {
    Point2d point{};
    point[axis] = along;
    point[other] = across;
    return point;
  };
  for (std::size_t i = 1; i < points.size(); ++i) {
    const Point2d& previous = points[i - 1];
    const Point2d& point = points[i];
    const double delta = point[axis] - previous[axis];
    if (std::fabs(delta) <= half) {
      current.push_back(point);
      continue;
    }
    // The step wraps across the seam: `rising` marks a jump from the high seam
    // edge down to the low one.
    const bool rising = delta < 0.0;
    if (dominant_sign != 0) {
      // A forward chain wraps high->low; a backward chain wraps low->high.
      const bool crossing = dominant_sign > 0 ? rising : !rising;
      if (!crossing) {
        current.push_back(point);
        continue;
      }
    } else if (!rising && !split_backward) {
      current.push_back(point);
      continue;
    }
    const double low_part = previous[axis] - base;
    const double high_part = high - point[axis];
    const double span = rising ? (high - previous[axis]) + (point[axis] - base)
                               : low_part + high_part;
    double t = 0.5;
    if (span > 0.0) t = rising ? (high - previous[axis]) / span : low_part / span;
    const double crossing = rising ? high : base;
    const double across = previous[other] + t * (point[other] - previous[other]);
    current.push_back(seam_point(crossing, across));
    finish();
    current.push_back(seam_point(rising ? base : high, across));
    current.push_back(point);
  }
  finish();
  return chains;
}

// Wraps every periodic axis of a chain into the domain window and splits it at
// each seam it crosses, so every returned piece lies inside the window with its
// open ends exactly on a seam edge.
inline std::vector<std::vector<Point2d>> trim_face_window_polyline(
    const std::vector<Point2d>& chain, const TrimFaceDomain& domain,
    double tolerance, bool split_backward, bool directional = false) {
  std::vector<std::vector<Point2d>> pieces{chain};
  for (int axis = 0; axis < 2; ++axis) {
    const bool periodic = axis == 0 ? domain.periodic_u : domain.periodic_v;
    if (!periodic) continue;
    const double base = domain.bounds[2 * static_cast<std::size_t>(axis)];
    const double high = domain.bounds[2 * static_cast<std::size_t>(axis) + 1];
    if (!(high > base)) continue;
    std::vector<std::vector<Point2d>> split;
    split.reserve(pieces.size());
    for (const std::vector<Point2d>& piece : pieces) {
      if (piece.size() < 2) continue;
      std::vector<Point2d> wrapped = piece;
      for (Point2d& point : wrapped) {
        point[axis] = trim_face_wrap_periodic(point[axis], base, high - base);
      }
      const int dominant = directional
                               ? trim_face_dominant_step_sign(
                                     wrapped, axis, 0.5 * (high - base))
                               : 0;
      for (std::vector<Point2d>& run :
           trim_face_split_at_seam(wrapped, axis, base, high, tolerance,
                                   split_backward, dominant)) {
        split.push_back(std::move(run));
      }
    }
    pieces.swap(split);
  }
  return pieces;
}

// True when every sample of a chain sits on one of the periodic parameter
// boundaries. Such a chain traces the seam edge the face already has: splitting
// the face on it would only shave slivers off the loop copy instead of cutting
// the face, so the caller drops it.
inline bool trim_face_chain_on_domain_side(
    const std::vector<Point2d>& points, const TrimFaceDomain& domain,
    double tolerance) {
  if (points.size() < 2) return false;
  for (int axis = 0; axis < 2; ++axis) {
    const bool periodic = axis == 0 ? domain.periodic_u : domain.periodic_v;
    if (!periodic) continue;
    const double base = domain.bounds[2 * static_cast<std::size_t>(axis)];
    const double high = domain.bounds[2 * static_cast<std::size_t>(axis) + 1];
    const double period = high - base;
    if (!(period > 0.0)) continue;
    const double slack = std::max(64.0 * tolerance, 1e-6 * period);
    bool on_side = true;
    for (const Point2d& point : points) {
      double offset = point[axis] - base;
      offset -= period * std::floor(offset / period);
      const double gap = std::min(offset, period - offset);
      if (gap > slack) {
        on_side = false;
        break;
      }
    }
    if (on_side) return true;
  }
  return false;
}

// Windows the cut polylines of one face onto its parameter rectangle. Point
// order is preserved and every seam crossing inserts the seam point itself, so
// a returned piece either stays closed or ends exactly on a domain side. An
// open piece must start and end on the face boundary, except that a piece
// touching a periodic seam may keep one end inside: its continuation lives on
// the other seam copy. A piece that touches no boundary at all is rejected
// instead of being silently passed through.
inline Result<std::vector<std::vector<Point2d>>> trim_face_window_chains(
    std::span<const FaceCut> cuts, const TrimFaceDomain& domain,
    double tolerance) {
  std::vector<std::vector<Point2d>> chains;
  for (const FaceCut& cut : cuts) {
    if (cut.uv.size() < 2) {
      return Error{Errc::invalid_topology,
                   "trim_face: a cut has fewer than two points"};
    }
    std::vector<Point2d> points = cut.uv;
    trim_face_dedupe_chain(points, tolerance);
    if (points.size() < 2) {
      return Error{Errc::invalid_topology,
                   "trim_face: a cut has fewer than two points"};
    }
    for (std::vector<Point2d>& piece : trim_face_window_polyline(
             points, domain, tolerance, false, true)) {
      if (piece.size() < 2) continue;
      const bool closed =
          trim_face_same_xy_rel(piece.front(), piece.back(), tolerance);
      const bool front_on_side =
          trim_face_point_on_domain_side(piece.front(), domain, tolerance);
      const bool back_on_side =
          trim_face_point_on_domain_side(piece.back(), domain, tolerance);
      const auto on_identified = [&domain, tolerance](const Point2d& p) {
        for (int side = 0; side < 4; ++side) {
          if (domain.sides[side] == TrimSideKind::kRealBoundary) continue;
          const double value = side < 2 ? p[0] : p[1];
          if (trim_face_near_side(value, domain.bounds[side], tolerance)) {
            return true;
          }
        }
        return false;
      };
      const bool touches_identified =
          on_identified(piece.front()) || on_identified(piece.back());
      if (!closed && (!front_on_side || !back_on_side) &&
          !touches_identified) {
        return Error{Errc::unsupported,
                     "trim_face: cut chain ends inside the parameter domain"};
      }
      chains.push_back(std::move(piece));
    }
  }
  return chains;
}

}  // namespace detail
}  // namespace v1
}  // namespace geom
