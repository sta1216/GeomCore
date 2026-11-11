#include <doctest/doctest.h>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <numbers>
#include <optional>
#include <string_view>
#include <variant>
#include <vector>

#include <geom/brep/assemble_curved.hpp>
#include <geom/brep/detail/curved_volume.hpp>
#include <geom/brep/face_uv.hpp>
#include <geom/brep/plane_frame.hpp>
#include <geom/brep/primitives.hpp>
#include <geom/brep/topology.hpp>
#include <geom/brep/validate.hpp>
#include <geom/core/error.hpp>
#include <geom/curves/analytic_curve.hpp>
#include <geom/primitives/point.hpp>
#include <geom/primitives/vec.hpp>
#include <geom/surfaces/plane.hpp>

#include "brep_curved_fixture.hpp"

namespace {

constexpr double kPi = std::numbers::pi_v<double>;
constexpr double kTwoPi = 2.0 * kPi;
constexpr double kTolerance = 1e-9;

// Reads the faces of an existing solid back into the curved fragment form the
// assembler consumes: 3D rings, their parameter rings, and the per-segment arc
// provenance derived from the uv loops. The fixture pcurves were authored
// against exactly these plane frames, so face_uv_eval round-trips the rings.
std::vector<geom::CurvedAssembledFace> fragments_from(const geom::Brep& brep, double tolerance) {
  std::vector<geom::CurvedAssembledFace> fragments;
  fragments.reserve(brep.face_count());
  for (const geom::FaceId face_id : brep.face_ids()) {
    const auto projection =
        geom::detail::face_uv_projection(brep, brep.face(face_id), tolerance);
    REQUIRE(projection.has_value());
    const auto loops = geom::face_loops_uv(brep, face_id, tolerance);
    REQUIRE(loops.has_value());

    geom::CurvedAssembledFace fragment;
    fragment.surface = brep.geometry().surfaces[brep.face(face_id).surface];
    for (const geom::FaceUvLoop& loop : loops.value()) {
      std::vector<geom::Point3d> points;
      points.reserve(loop.points.size());
      for (const geom::Point2d& uv : loop.points) {
        points.push_back(geom::detail::face_uv_eval(projection.value(), uv));
      }
      fragment.loops.push_back(std::move(points));
      fragment.loops_uv.push_back(loop.points);
    }
    for (const geom::FaceUvLoop& loop : loops.value()) {
      std::vector<std::optional<geom::CurvedArcSpan>> spans(loop.points.size());
      for (std::size_t i = 0; i < loop.arcs.size(); ++i) {
        if (!loop.arcs[i].has_value()) continue;
        const geom::FaceUvArcRef& arc = loop.arcs[i].value();
        spans[i] = geom::CurvedArcSpan{arc.circle, std::min(arc.t0, arc.t1),
                                       std::max(arc.t0, arc.t1)};
      }
      fragment.arc_spans.push_back(std::move(spans));
    }
    fragments.push_back(std::move(fragment));
  }
  return fragments;
}

geom::Point2d plane_uv(const geom::PlaneFrame& frame, const geom::Point3d& point) {
  return geom::plane_frame_to_2d(frame, point);
}

// Evaluates a stored pcurve ring on the face surface and returns the 3D points.
std::vector<geom::Point3d> pcurve_positions(const geom::Brep& brep, const geom::Face& face,
                                            const std::vector<geom::Point2d>& points) {
  const auto projection = geom::detail::face_uv_projection(brep, face, kTolerance);
  REQUIRE(projection.has_value());
  std::vector<geom::Point3d> positions;
  positions.reserve(points.size());
  for (const geom::Point2d& uv : points) {
    positions.push_back(geom::detail::face_uv_eval(projection.value(), uv));
  }
  return positions;
}

}  // namespace

TEST_CASE("assemble_curved_solid rebuilds a cylinder from its three faces") {
  const geom::Brep cylinder =
      geom::make_cylinder({0.0, 0.0, 0.0}, {0.0, 0.0, 1.0}, 2.0, 3.0, kTolerance).value();
  const std::vector<geom::CurvedAssembledFace> faces = fragments_from(cylinder, kTolerance);
  REQUIRE(faces.size() == 3u);

  const auto result = geom::assemble_curved_solid(faces, kTolerance);
  REQUIRE(result.has_value());
  const geom::Brep& solid = result.value();

  const geom::ValidationReport report = geom::validate_topology(solid);
  CHECK(report.ok());
  CHECK(solid.face_count() == 3u);
  REQUIRE(solid.solid_ids().size() == 1u);
  const geom::SolidId solid_id = solid.solid_ids().front();
  CHECK(geom::euler_characteristic(solid, solid_id) == 2);

  std::size_t rim_edges = 0;
  for (const geom::EdgeId edge_id : solid.edge_ids()) {
    const geom::Edge& edge = solid.edge(edge_id);
    REQUIRE(edge.curve != geom::kNoGeometry);
    const geom::CurveGeometry& curve = solid.geometry().curves[edge.curve];
    if (!std::holds_alternative<geom::CircleCurve3>(curve)) continue;
    ++rim_edges;
    const geom::CircleCurve3& circle = std::get<geom::CircleCurve3>(curve);
    CHECK(circle.is_valid());
    CHECK(std::abs(circle.radius - 2.0) <= 1e-12);
    CHECK(edge.t0 == 0.0);
    CHECK(std::abs((edge.t1 - edge.t0) - kTwoPi) <= 1e-12);
    CHECK(edge.v0 == edge.v1);
  }
  CHECK(rim_edges == 2u);

  // Every coedge carries a dense pcurve along the edge direction: it starts at
  // the edge's v0, ends at its v1, is finite, and a full rim pcurve is closed.
  for (const geom::CoedgeId coedge_id : solid.coedge_ids()) {
    const geom::Coedge& coedge = solid.coedge(coedge_id);
    REQUIRE(coedge.pcurve != geom::kNoGeometry);
    const std::vector<geom::Point2d>& points = solid.geometry().pcurves[coedge.pcurve].points;
    CHECK(points.size() >= 2u);
    for (const geom::Point2d& uv : points) {
      CHECK(std::isfinite(uv[0]));
      CHECK(std::isfinite(uv[1]));
    }
    const geom::Edge& edge = solid.edge(coedge.edge);
    const std::vector<geom::Point3d> positions =
        pcurve_positions(solid, solid.face(coedge.face), points);
    const geom::Point3d v0 = solid.vertex(edge.v0).point;
    const geom::Point3d v1 = solid.vertex(edge.v1).point;
    CHECK(geom::norm2(positions.front() - v0) <= kTolerance);
    CHECK(geom::norm2(positions.back() - v1) <= kTolerance);
  }

  const auto volume = geom::detail::curved_solid_volume(solid, solid_id, kTolerance);
  REQUIRE(volume.has_value());
  CHECK(std::abs(volume.value() - kPi * 4.0 * 3.0) <= 1e-9);
}

TEST_CASE("assemble_curved_solid rebuilds a faceted quarter cylinder") {
  const geom::Brep quarter = geom::test::make_quarter_cylinder(2.0, 3.0, kTolerance);
  const std::vector<geom::CurvedAssembledFace> faces = fragments_from(quarter, kTolerance);
  REQUIRE(faces.size() == 5u);

  const auto result = geom::assemble_curved_solid(faces, kTolerance);
  REQUIRE(result.has_value());
  const geom::Brep& solid = result.value();
  CHECK(geom::validate_topology(solid).ok());
  CHECK(solid.face_count() == 5u);
  REQUIRE(solid.solid_ids().size() == 1u);
  const geom::SolidId solid_id = solid.solid_ids().front();
  CHECK(geom::euler_characteristic(solid, solid_id) == 2);

  for (const geom::EdgeId edge_id : solid.edge_ids()) {
    const geom::Edge& edge = solid.edge(edge_id);
    REQUIRE(edge.curve != geom::kNoGeometry);
    CHECK_FALSE(std::holds_alternative<geom::CircleCurve3>(
        solid.geometry().curves[edge.curve]));
  }

  // The lateral face keeps its exact cylinder patch, so the divergence-theorem
  // value is the smooth lateral integral plus the inscribed 16-gon cap term
  // (see tests/unit/brep_curved_volume_test.cpp). The faceted prism value
  // h * (n / 2) * r^2 * sin(pi / (2 n)) is the same solid's faceted
  // idealization, so it can only be matched to the faceting error below.
  constexpr double kSegments = 16.0;
  const double cap_area =
      0.5 * kSegments * 4.0 * std::sin((kPi / 2.0) / kSegments);
  const double expected = 3.0 * (kPi * 4.0 / 2.0 + cap_area) / 3.0;
  const auto volume = geom::detail::curved_solid_volume(solid, solid_id, kTolerance);
  REQUIRE(volume.has_value());
  CHECK(std::abs(volume.value() - expected) <= 1e-9);

  const double prism = 3.0 * (kSegments / 2.0) * 4.0 * std::sin(kPi / (2.0 * kSegments));
  CHECK(std::abs(volume.value() - prism) <= 2e-3 * prism);
}

TEST_CASE("assemble_curved_solid rejects an arc shared with a chord") {
  const geom::Plane3 plane{{0.0, 0.0, 1.0}, 0.0};
  const auto frame = geom::make_plane_frame(plane);
  REQUIRE(frame.has_value());

  const geom::Point3d p0{1.0, 0.0, 0.0};
  const geom::Point3d p1{-1.0, 0.0, 0.0};
  const geom::Point3d p2{0.0, -2.0, 0.0};

  geom::CurvedAssembledFace arc_face;
  arc_face.surface = plane;
  arc_face.loops = {{p0, p1, p2}};
  arc_face.loops_uv = {{plane_uv(frame.value(), p0), plane_uv(frame.value(), p1),
                        plane_uv(frame.value(), p2)}};
  geom::CurvedArcSpan span;
  span.circle = geom::CircleCurve3{{0.0, 0.0, 0.0}, {1.0, 0.0, 0.0}, {0.0, 1.0, 0.0}, 1.0};
  span.t0 = 0.0;
  span.t1 = kPi;
  arc_face.arc_spans = {{{span}, std::nullopt, std::nullopt}};

  geom::CurvedAssembledFace chord_face;
  chord_face.surface = plane;
  chord_face.loops = {{p0, p2, p1}};
  chord_face.loops_uv = {{plane_uv(frame.value(), p0), plane_uv(frame.value(), p2),
                          plane_uv(frame.value(), p1)}};
  chord_face.arc_spans = {};

  const std::vector<geom::CurvedAssembledFace> faces{arc_face, chord_face};
  const auto result = geom::assemble_curved_solid(faces, kTolerance);
  REQUIRE_FALSE(result.has_value());
  CHECK(result.error().code == geom::Errc::invalid_topology);
  CHECK(std::string_view(result.error().message) ==
        "assemble_curved_solid: shared edge curve mismatch");
}

TEST_CASE("assemble_curved_solid rejects an edge with a single use") {
  const geom::Plane3 plane{{0.0, 0.0, 1.0}, 0.0};
  const auto frame = geom::make_plane_frame(plane);
  REQUIRE(frame.has_value());
  const geom::Point3d p0{0.0, 0.0, 0.0};
  const geom::Point3d p1{1.0, 0.0, 0.0};
  const geom::Point3d p2{0.0, 1.0, 0.0};

  geom::CurvedAssembledFace face;
  face.surface = plane;
  face.loops = {{p0, p1, p2}};
  face.loops_uv = {{plane_uv(frame.value(), p0), plane_uv(frame.value(), p1),
                    plane_uv(frame.value(), p2)}};
  face.arc_spans = {};

  const std::vector<geom::CurvedAssembledFace> faces{face};
  const auto result = geom::assemble_curved_solid(faces, kTolerance);
  REQUIRE_FALSE(result.has_value());
  CHECK(result.error().code == geom::Errc::invalid_topology);
  CHECK(std::string_view(result.error().message) ==
        "assemble_curved_solid: edge is used once");
}

TEST_CASE("assemble_curved_solid rejects an arc run whose points do not match the span") {
  const geom::Plane3 plane{{0.0, 0.0, 1.0}, 0.0};
  const auto frame = geom::make_plane_frame(plane);
  REQUIRE(frame.has_value());
  const geom::Point3d p0{1.0, 0.0, 0.0};
  const geom::Point3d p1{0.0, 1.0, 0.0};
  const geom::Point3d p2{0.0, 0.0, 0.0};

  geom::CurvedAssembledFace face;
  face.surface = plane;
  face.loops = {{p0, p1, p2}};
  face.loops_uv = {{plane_uv(frame.value(), p0), plane_uv(frame.value(), p1),
                    plane_uv(frame.value(), p2)}};
  geom::CurvedArcSpan span;
  span.circle = geom::CircleCurve3{{0.0, 0.0, 0.0}, {1.0, 0.0, 0.0}, {0.0, 1.0, 0.0}, 1.0};
  span.t0 = 0.0;
  span.t1 = kPi;  // a half turn, but the loop points span a quarter turn
  face.arc_spans = {{{span}, std::nullopt, std::nullopt}};

  const std::vector<geom::CurvedAssembledFace> faces{face};
  const auto result = geom::assemble_curved_solid(faces, kTolerance);
  REQUIRE_FALSE(result.has_value());
  CHECK(result.error().code == geom::Errc::invalid_topology);
  CHECK(std::string_view(result.error().message) ==
        "assemble_curved_solid: arc run is inconsistent");
}
