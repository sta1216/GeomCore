#include <doctest/doctest.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <limits>
#include <numbers>
#include <string>
#include <string_view>
#include <vector>

#include <geom/brep/heal.hpp>
#include <geom/brep/detail/curved_volume.hpp>

#include "brep_heal_fixture.hpp"
#include "m4a_curved_heal_fixture.hpp"
#include "m4b_heal_fixture.hpp"

using namespace geom;

namespace {

constexpr std::string_view kUnimplementedMessage =
    "heal_brep: curved geometry is not supported yet";
constexpr std::string_view kSurfaceKindMessage =
    "heal_brep: surface type is not supported yet";
constexpr std::string_view kCurveKindMessage =
    "heal_brep: curve type is not supported yet";
constexpr std::string_view kCarrierMessage =
    "heal_brep: gap faces do not share the same surface";
constexpr std::string_view kOffSurfaceMessage =
    "heal_brep: gap edge is not on its face surface";
constexpr std::string_view kKeptCurveMessage =
    "heal_brep: kept edge curve does not match its welded vertices";
constexpr std::string_view kInvalidSurfaceMessage =
    "heal_brep: surface is invalid";

bool contains(std::string_view text, std::string_view needle) {
  return text.find(needle) != std::string_view::npos;
}

// Counts plus every vertex position, so a rejection that mutated the input
// would be visible even if the entity counts stayed the same.
std::vector<double> brep_fingerprint(const Brep& brep) {
  std::vector<double> fingerprint{
      static_cast<double>(brep.vertex_count()),
      static_cast<double>(brep.edge_count()),
      static_cast<double>(brep.face_count())};
  for (const VertexId id : brep.vertex_ids()) {
    const Point3d point = brep.vertex(id).point;
    fingerprint.push_back(point[0]);
    fingerprint.push_back(point[1]);
    fingerprint.push_back(point[2]);
  }
  return fingerprint;
}

// Full topological fingerprint for rejection immutability: entity/store
// counts, every vertex position and tolerance, every edge's carrier and
// parameter range with both endpoints, every coedge's wiring/pcurve/tolerance,
// every loop, and every face's surface/reversal/tolerance. Handles are recorded
// by index (null as the sentinel) so the copy must also agree on which slots are
// live, not only on the surviving geometry values.
std::vector<double> full_brep_fingerprint(const Brep& brep) {
  std::vector<double> fingerprint{
      static_cast<double>(brep.vertex_count()),
      static_cast<double>(brep.edge_count()),
      static_cast<double>(brep.coedge_count()),
      static_cast<double>(brep.loop_count()),
      static_cast<double>(brep.face_count()),
      static_cast<double>(brep.shell_count()),
      static_cast<double>(brep.solid_count()),
      static_cast<double>(brep.geometry().curves.size()),
      static_cast<double>(brep.geometry().surfaces.size()),
      static_cast<double>(brep.geometry().pcurves.size())};
  const auto index_of = [](const Handle& id) {
    return static_cast<double>(id.is_null() ? Handle::kNullIndex : id.index());
  };
  for (const VertexId id : brep.vertex_ids()) {
    const Vertex& vertex = brep.vertex(id);
    fingerprint.push_back(static_cast<double>(id.index()));
    fingerprint.push_back(vertex.point[0]);
    fingerprint.push_back(vertex.point[1]);
    fingerprint.push_back(vertex.point[2]);
    fingerprint.push_back(vertex.tol);
  }
  for (const EdgeId id : brep.edge_ids()) {
    const Edge& edge = brep.edge(id);
    fingerprint.push_back(static_cast<double>(id.index()));
    fingerprint.push_back(static_cast<double>(edge.curve));
    fingerprint.push_back(edge.t0);
    fingerprint.push_back(edge.t1);
    fingerprint.push_back(index_of(edge.v0));
    fingerprint.push_back(index_of(edge.v1));
    fingerprint.push_back(edge.tol);
  }
  for (const CoedgeId id : brep.coedge_ids()) {
    const Coedge& coedge = brep.coedge(id);
    fingerprint.push_back(static_cast<double>(id.index()));
    fingerprint.push_back(index_of(coedge.edge));
    fingerprint.push_back(index_of(coedge.face));
    fingerprint.push_back(index_of(coedge.loop));
    fingerprint.push_back(coedge.reversed ? 1.0 : 0.0);
    fingerprint.push_back(static_cast<double>(coedge.pcurve));
    fingerprint.push_back(index_of(coedge.next));
    fingerprint.push_back(index_of(coedge.prev));
    fingerprint.push_back(coedge.tol);
  }
  for (const LoopId id : brep.loop_ids()) {
    const Loop& loop = brep.loop(id);
    fingerprint.push_back(static_cast<double>(id.index()));
    fingerprint.push_back(index_of(loop.face));
    fingerprint.push_back(index_of(loop.start));
    fingerprint.push_back(index_of(loop.next));
    fingerprint.push_back(loop.outer ? 1.0 : 0.0);
  }
  for (const FaceId id : brep.face_ids()) {
    const Face& face = brep.face(id);
    fingerprint.push_back(static_cast<double>(id.index()));
    fingerprint.push_back(static_cast<double>(face.surface));
    fingerprint.push_back(index_of(face.first_loop));
    fingerprint.push_back(face.reversed ? 1.0 : 0.0);
    fingerprint.push_back(face.tol);
  }
  return fingerprint;
}

std::size_t count_issues(const detail::HealPlan& plan, HealIssueKind kind) {
  return static_cast<std::size_t>(std::count_if(
      plan.issues.begin(), plan.issues.end(),
      [kind](const HealIssue& issue) { return issue.kind == kind; }));
}

// Signed uv area of a closed polygon (shoelace formula).
double shoelace_area(const std::vector<Point2d>& points) {
  double twice = 0.0;
  for (std::size_t i = 0; i < points.size(); ++i) {
    const Point2d& a = points[i];
    const Point2d& b = points[(i + 1) % points.size()];
    twice += a[0] * b[1] - b[0] * a[1];
  }
  return 0.5 * twice;
}

struct UvBounds {
  double u_min = std::numeric_limits<double>::infinity();
  double u_max = -std::numeric_limits<double>::infinity();
  double v_min = std::numeric_limits<double>::infinity();
  double v_max = -std::numeric_limits<double>::infinity();
};

UvBounds uv_bounds(const std::vector<Point2d>& points) {
  UvBounds bounds;
  for (const Point2d& point : points) {
    bounds.u_min = std::min(bounds.u_min, point[0]);
    bounds.u_max = std::max(bounds.u_max, point[0]);
    bounds.v_min = std::min(bounds.v_min, point[1]);
    bounds.v_max = std::max(bounds.v_max, point[1]);
  }
  return bounds;
}

// A sphere lune is a parameter rectangle, so its area is the exact integral of
// the surface metric over the loop's uv bounds. Summing the uv polygon instead
// would drop the polar sliver where the pcurve chord meets the pole.
double sphere_area(const std::vector<FaceUvLoop>& loops, double radius) {
  double area = 0.0;
  for (const FaceUvLoop& loop : loops) {
    const UvBounds bounds = uv_bounds(loop.points);
    area += radius * radius * (bounds.u_max - bounds.u_min) *
            (std::sin(bounds.v_max) - std::sin(bounds.v_min));
  }
  return area;
}

// Concatenates the stored pcurves of the face's first loop in coedge order and
// returns the uv polygon those pcurves describe. Junction points are kept: the
// chart reader merges a pole/apex revisit by 3D position, which cuts the uv
// image of a closed lune, so the stored pcurves are the honest uv chart.
std::vector<Point2d> pcurve_loop_points(const Brep& brep, FaceId face_id) {
  const Face& face = brep.face(face_id);
  if (face.first_loop.is_null()) return {};
  const Loop& loop = brep.loop(face.first_loop);
  std::vector<Point2d> points;
  CoedgeId coedge_id = loop.start;
  for (std::size_t step = 0; step <= brep.coedge_count(); ++step) {
    if (!brep.coedge_alive(coedge_id)) return {};
    const Coedge& coedge = brep.coedge(coedge_id);
    if (coedge.pcurve == kNoGeometry) return {};
    const Pcurve2d& pcurve = brep.geometry().pcurves[coedge.pcurve];
    if (coedge.reversed) {
      for (auto it = pcurve.points.rbegin(); it != pcurve.points.rend(); ++it) {
        points.push_back(*it);
      }
    } else {
      points.insert(points.end(), pcurve.points.begin(), pcurve.points.end());
    }
    if (coedge.next == loop.start) break;
    coedge_id = coedge.next;
  }
  return points;
}

// Every coedge of `face_id` must carry a pcurve, circle edges must be sampled
// with interior points, and every stored uv point must evaluate back onto its
// carrier curve and re-project to the same parameter.
void check_face_pcurves(const Brep& brep, FaceId face_id) {
  const Result<detail::FaceUvProjection> projection =
      detail::face_uv_projection(brep, brep.face(face_id), 1e-9);
  REQUIRE(projection.has_value());
  const Result<std::vector<FaceUvLoop>> loops =
      face_loops_uv(brep, face_id, 1e-9);
  REQUIRE(loops.has_value());

  std::size_t coedges = 0;
  for (const CoedgeId coedge_id : brep.coedge_ids()) {
    const Coedge& coedge = brep.coedge(coedge_id);
    if (coedge.face != face_id) continue;
    ++coedges;
    REQUIRE(coedge.pcurve != kNoGeometry);
    REQUIRE(coedge.pcurve < brep.geometry().pcurves.size());
    const Pcurve2d& pcurve = brep.geometry().pcurves[coedge.pcurve];
    REQUIRE(pcurve.points.size() >= 2u);

    const Edge& edge = brep.edge(coedge.edge);
    REQUIRE(edge.curve != kNoGeometry);
    REQUIRE(edge.curve < brep.geometry().curves.size());
    const CurveGeometry& curve = brep.geometry().curves[edge.curve];
    if (std::get_if<CircleCurve3>(&curve) != nullptr) {
      CHECK(pcurve.points.size() > 2u);
      for (const Point2d& uv : pcurve.points) {
        const Point3d point = detail::face_uv_eval(projection.value(), uv);
        CHECK(detail::heal_curve_residual(curve, point) <= 1e-9);
        const Result<Point2d> back =
            detail::face_uv_project(projection.value(), point, uv[0], true);
        REQUIRE(back.has_value());
        double u = back.value()[0];
        if (projection.value().periodic_u()) {
          u = detail::face_uv_unwrap(u, uv[0]);
        }
        CHECK(std::fabs(u - uv[0]) <= 1e-9);
        CHECK(std::fabs(back.value()[1] - uv[1]) <= 1e-9);
      }
    } else {
      CHECK(std::get_if<LineCurve3>(&curve) != nullptr);
      CHECK(pcurve.points.size() == 2u);
    }
  }
  CHECK(coedges > 0u);
}

}  // namespace

TEST_CASE("heal: an unsupported NURBS face names the face index") {
  const test::m4a::UnsupportedFaceSolid fixture =
      test::m4a::make_nurbs_face_solid();

  const auto healed = heal_brep(fixture.brep, fixture.solid, {});
  REQUIRE_FALSE(healed.has_value());
  CHECK(healed.error().code == Errc::unsupported);
  CHECK(std::string_view(healed.error().message) == kSurfaceKindMessage);
  CHECK(healed.error().context_a == fixture.face.index());
}

TEST_CASE("heal: an invalid analytic surface names the surface index") {
  // probe8 built the sphere as Sphere3{center, radius}; brace elision put the
  // radius in `axis`, leaving axis == x_axis. That chart is degenerate, so the
  // surface is invalid even though the topology below is a valid B-rep. A
  // leading valid plane pushes the malformed sphere to surface index 1, while
  // the first face (the band) sits on it, so a face-index context would fail.
  const test::m4a::CurvedSolid fixture =
      test::m4b::make_malformed_sphere_pole_gap(1e-6, 5e-7, true);
  REQUIRE(fixture.brep.geometry().surfaces.size() == 3u);
  CHECK(std::get_if<Plane3>(&fixture.brep.geometry().surfaces[0]) != nullptr);
  const auto* sphere =
      std::get_if<Sphere3>(&fixture.brep.geometry().surfaces[1]);
  REQUIRE(sphere != nullptr);
  REQUIRE_FALSE(sphere->is_valid());
  const ShellId shell = fixture.brep.solid_shells(fixture.solid).front();
  const FaceId band = fixture.brep.shell_faces(shell).front();
  CHECK(fixture.brep.face(band).surface == 1u);
  const std::vector<double> input_fingerprint =
      full_brep_fingerprint(fixture.brep);

  const auto healed = heal_brep(fixture.brep, fixture.solid, {});
  if (healed.has_value()) {
    MESSAGE("heal_brep accepted an invalid analytic surface");
  }
  REQUIRE_FALSE(healed.has_value());
  CHECK(healed.error().code == Errc::unsupported);
  CHECK(std::string_view(healed.error().message) == kInvalidSurfaceMessage);
  // The context is the surface store index, not the face that uses it.
  CHECK(healed.error().context_a == 1u);
  CHECK(full_brep_fingerprint(fixture.brep) == input_fingerprint);

  // Repeating the heal is deterministic and still leaves the input untouched.
  const auto repeated = heal_brep(fixture.brep, fixture.solid, {});
  REQUIRE_FALSE(repeated.has_value());
  CHECK(repeated.error().code == Errc::unsupported);
  CHECK(std::string_view(repeated.error().message) == kInvalidSurfaceMessage);
  CHECK(repeated.error().context_a == 1u);
  CHECK(full_brep_fingerprint(fixture.brep) == input_fingerprint);
}

TEST_CASE("heal: an unsupported NURBS edge curve names the curve kind") {
  const test::m4a::UnsupportedCurveSolid fixture =
      test::m4a::make_nurbs_edge_solid();

  const auto healed = heal_brep(fixture.brep, fixture.solid, {});
  REQUIRE_FALSE(healed.has_value());
  CHECK(healed.error().code == Errc::unsupported);
  CHECK(std::string_view(healed.error().message) == kCurveKindMessage);
}

TEST_CASE("heal: a two-band cylinder gap is no longer unimplemented") {
  const test::m4a::CurvedSolid cylinder =
      test::m4a::make_two_band_cylinder(1.0, 1e-6);

  const auto healed = heal_brep(cylinder.brep, cylinder.solid, {});
  if (!healed.has_value()) {
    CHECK_FALSE(contains(healed.error().message, kUnimplementedMessage));
  } else {
    for (const HealIssue& issue : healed.value().report.issues) {
      CHECK_FALSE(contains(issue.message, kUnimplementedMessage));
    }
  }
}

TEST_CASE("heal: a cylinder gap tie-break picks the seam orientation") {
  const test::m4a::CurvedSolid cylinder =
      test::m4a::make_two_band_cylinder(1.0, 1e-6);
  const std::vector<EdgeId> boundary =
      detail::heal_boundary_edges(cylinder.brep, cylinder.solid);
  REQUIRE(boundary.size() == 4u);

  const auto plan_result = detail::heal_plan(cylinder.brep, boundary, 1e-6);
  REQUIRE(plan_result.has_value());
  const detail::HealPlan& plan = plan_result.value();
  CHECK(plan.pairs.size() == 2u);
  CHECK(count_issues(plan, HealIssueKind::ambiguous_gap) == 0u);
  for (const detail::HealEdgePair& pair : plan.pairs) {
    CHECK_FALSE(pair.same_direction);
  }
}

TEST_CASE("heal: a sphere lune gap tie-break picks the seam orientation") {
  const test::m4a::CurvedSolid sphere = test::m4a::make_sphere_gap(5e-7);
  const std::vector<EdgeId> boundary =
      detail::heal_boundary_edges(sphere.brep, sphere.solid);
  REQUIRE(boundary.size() == 4u);

  const auto plan_result = detail::heal_plan(sphere.brep, boundary, 1e-6);
  REQUIRE(plan_result.has_value());
  const detail::HealPlan& plan = plan_result.value();
  CHECK(plan.pairs.size() == 2u);
  CHECK(count_issues(plan, HealIssueKind::ambiguous_gap) == 0u);
  for (const detail::HealEdgePair& pair : plan.pairs) {
    CHECK_FALSE(pair.same_direction);
  }
}

TEST_CASE("heal: planar ambiguous rims keep ambiguous_gap through the tie-break") {
  const test::TornSheet sheet = test::make_ambiguous_rims(1e-6);
  const std::vector<EdgeId> boundary =
      detail::heal_boundary_edges(sheet.brep, sheet.solid);
  REQUIRE(boundary.size() == 12u);

  const auto plan_result = detail::heal_plan(sheet.brep, boundary, 1.5e-6);
  REQUIRE(plan_result.has_value());
  const detail::HealPlan& plan = plan_result.value();
  CHECK(plan.pairs.empty());
  REQUIRE(plan.issues.size() == 4u);
  std::vector<std::uint32_t> ambiguous_edges;
  for (const HealIssue& issue : plan.issues) {
    CHECK(issue.kind == HealIssueKind::ambiguous_gap);
    ambiguous_edges.push_back(issue.entity_a);
  }
  CHECK(ambiguous_edges == std::vector<std::uint32_t>{0, 1, 2, 3});
}

TEST_CASE("heal: a cylinder gap with mismatched carriers is rejected") {
  const test::m4a::CurvedSolid cylinder =
      test::m4a::make_cylinder_radius_gap(2e-7, 1.0 + 2e-7);
  const std::vector<double> before = brep_fingerprint(cylinder.brep);

  const auto healed = heal_brep(cylinder.brep, cylinder.solid, {});
  REQUIRE_FALSE(healed.has_value());
  CHECK(healed.error().code == Errc::unsupported);
  CHECK(std::string_view(healed.error().message) == kCarrierMessage);
  CHECK(healed.error().context_a == 2u);
  CHECK(healed.error().context_b == 4u);
  CHECK(brep_fingerprint(cylinder.brep) == before);
}

TEST_CASE("heal: a plane lid over a cylinder wall is rejected") {
  const test::m4a::CurvedSolid cylinder =
      test::m4a::make_cylinder_plane_gap(2e-7);
  const std::vector<double> before = brep_fingerprint(cylinder.brep);

  const auto healed = heal_brep(cylinder.brep, cylinder.solid, {});
  REQUIRE_FALSE(healed.has_value());
  CHECK(healed.error().code == Errc::unsupported);
  CHECK(std::string_view(healed.error().message) == kCarrierMessage);
  CHECK(healed.error().context_a == 2u);
  CHECK(healed.error().context_b == 4u);
  CHECK(brep_fingerprint(cylinder.brep) == before);
}

TEST_CASE("heal: gap rims floating off their cylinder are rejected") {
  const test::m4a::CurvedSolid cylinder =
      test::m4a::make_cylinder_off_surface_gap(2e-7);
  const std::vector<double> before = brep_fingerprint(cylinder.brep);

  const auto healed = heal_brep(cylinder.brep, cylinder.solid, {});
  REQUIRE_FALSE(healed.has_value());
  CHECK(healed.error().code == Errc::unsupported);
  CHECK(std::string_view(healed.error().message) == kOffSurfaceMessage);
  CHECK(healed.error().context_a == 2u);
  CHECK(healed.error().context_b == 4u);
  CHECK(brep_fingerprint(cylinder.brep) == before);
}

TEST_CASE("heal: a kept curve must pass through its welded vertices") {
  const test::m4a::CurvedSolid cylinder =
      test::m4a::make_cylinder_kept_curve_gap(2e-7);
  const std::vector<double> before = brep_fingerprint(cylinder.brep);

  const auto healed = heal_brep(cylinder.brep, cylinder.solid, {});
  REQUIRE_FALSE(healed.has_value());
  CHECK(healed.error().code == Errc::unsupported);
  CHECK(std::string_view(healed.error().message) == kKeptCurveMessage);
  CHECK(healed.error().context_a == 2u);
  CHECK(healed.error().context_b == 4u);
  CHECK(brep_fingerprint(cylinder.brep) == before);
}

TEST_CASE("heal: a radius mismatch beyond max_gap stays a gap_too_large report") {
  const test::m4a::CurvedSolid cylinder =
      test::m4a::make_cylinder_radius_gap(1e-3, 1.001);

  const auto healed = heal_brep(cylinder.brep, cylinder.solid, {});
  REQUIRE(healed.has_value());
  CHECK_FALSE(healed.value().report.changed);
  REQUIRE_FALSE(healed.value().report.issues.empty());
  const HealIssue& issue = healed.value().report.issues.front();
  CHECK(issue.kind == HealIssueKind::gap_too_large);
  CHECK(issue.entity_a == 2u);
  CHECK(issue.entity_b == 4u);
}

TEST_CASE("heal: a curved gap rebuilds arc pcurves") {
  constexpr double gap = 1e-6;

  SUBCASE("cylinder") {
    const test::m4a::CurvedSolid fixture =
        test::m4a::make_two_band_cylinder(1.0, gap);
    const auto healed = heal_brep(fixture.brep, fixture.solid, {});
    REQUIRE(healed.has_value());
    CHECK(healed.value().report.changed);
    CHECK(healed.value().report.issues.empty());
    const Brep& brep = healed.value().brep;

    std::size_t rebuilt = 0;
    for (const FaceId face_id : brep.face_ids()) {
      bool any = false;
      bool all = true;
      for (const CoedgeId coedge_id : brep.coedge_ids()) {
        if (brep.coedge(coedge_id).face != face_id) continue;
        const bool has = brep.coedge(coedge_id).pcurve != kNoGeometry;
        any = any || has;
        all = all && has;
      }
      CHECK(any == all);
      if (any) {
        ++rebuilt;
        check_face_pcurves(brep, face_id);
      }
    }
    CHECK(rebuilt == 2u);

    std::size_t circles = 0;
    std::size_t lines = 0;
    for (const EdgeId edge_id : brep.edge_ids()) {
      const Edge& edge = brep.edge(edge_id);
      REQUIRE(edge.curve != kNoGeometry);
      const CurveGeometry& curve = brep.geometry().curves[edge.curve];
      if (std::get_if<CircleCurve3>(&curve) != nullptr) {
        ++circles;
      } else {
        REQUIRE(std::get_if<LineCurve3>(&curve) != nullptr);
        ++lines;
      }
    }
    CHECK(circles == 6u);
    CHECK(lines == 2u);


    std::size_t lateral_points = 0;
    std::size_t cap_points = 0;
    double radius = 0.0;
    double lateral_uv_area = 0.0;
    for (const FaceId face_id : brep.face_ids()) {
      const auto loops = face_loops_uv(brep, face_id, 1e-9);
      REQUIRE(loops.has_value());
      REQUIRE(loops.value().size() == 1u);
      const std::vector<Point2d>& points = loops.value().front().points;
      const SurfaceGeometry& surface =
          brep.geometry().surfaces[brep.face(face_id).surface];
      if (const auto* cylinder = std::get_if<Cylinder3>(&surface)) {
        lateral_points = points.size();
        radius = cylinder->radius;
        lateral_uv_area += std::fabs(shoelace_area(points));
      } else {
        cap_points = points.size();
      }
    }
    // Collinear simplification leaves one interior sample per circle coedge, so
    // the rebuilt lateral chart is compact (3+3+2+3+3+2 minus five junction
    // drops and the closing duplicate). The cap stores no pcurves and keeps the
    // sampled resolution of its two half circles.
    CHECK(lateral_points == 10u);
    CHECK(cap_points == 64u);
    CHECK(std::fabs(radius - 1.0) <= 1e-12);

    const double area = radius * lateral_uv_area;
    const double expected = test::m4a::kTwoPi * 2.0;
    CHECK(std::fabs(area - expected) <= 1e-6 * expected);
    CHECK(std::fabs(area - test::m4a::kTwoPi * (2.0 + gap)) <= 1e-9);
  }

  SUBCASE("sphere") {
    constexpr double eps = 5e-7;
    const test::m4a::CurvedSolid fixture = test::m4a::make_sphere_gap(eps);
    const auto healed = heal_brep(fixture.brep, fixture.solid, {});
    REQUIRE(healed.has_value());
    CHECK(healed.value().report.changed);
    CHECK(healed.value().report.issues.empty());
    const Brep& brep = healed.value().brep;
    REQUIRE(brep.face_count() == 4u);

    double area = 0.0;
    double pcurve_area = 0.0;
    double reader_area_total = 0.0;
    for (const FaceId face_id : brep.face_ids()) {
      check_face_pcurves(brep, face_id);
      const auto loops = face_loops_uv(brep, face_id, 1e-9);
      REQUIRE(loops.has_value());
      REQUIRE(loops.value().size() == 1u);
      const std::vector<Point2d>& points = loops.value().front().points;
      // Both pole corners of a lune must survive the chart reader: with one of
      // them dropped the rectangle loses a side and the uv area falls to about
      // 4.7806 instead of pi^2/2.
      const double reader_area = shoelace_area(points);
      CHECK(std::fabs(reader_area - 0.5 * test::m4a::kPi * test::m4a::kPi) <= 1e-5);
      reader_area_total += reader_area;
      const auto* sphere = std::get_if<Sphere3>(
          &brep.geometry().surfaces[brep.face(face_id).surface]);
      REQUIRE(sphere != nullptr);
      area += sphere_area(loops.value(), sphere->radius);
      // A lune chart is a rectangle u in [0, pi], v in [eps, pi/2] (or the
      // mirrored south band): the stored pcurves must cover it, including both
      // corners of the pole, so the uv area is pi^2/2 per lune.
      const std::vector<Point2d> pcurve_points = pcurve_loop_points(brep, face_id);
      CHECK(pcurve_points.size() >= 3u);
      const double lune_area = shoelace_area(pcurve_points);
      CHECK(std::fabs(lune_area - 0.5 * test::m4a::kPi * test::m4a::kPi) <= 1e-5);
      pcurve_area += lune_area;
    }
    CHECK(std::fabs(area - 4.0 * test::m4a::kPi) <= 1e-6);
    CHECK(std::fabs(pcurve_area - 2.0 * test::m4a::kPi * test::m4a::kPi) <= 1e-5);
    CHECK(std::fabs(reader_area_total - 2.0 * test::m4a::kPi * test::m4a::kPi) <= 1e-5);

    // Closing the pole corners must also leave the closed solid integrable: the
    // old cut chart made this fail with Errc::degenerate.
    const auto volume =
        detail::curved_solid_volume(brep, fixture.solid, 1e-9);
    REQUIRE(volume.has_value());
    const double expected_volume = 4.0 / 3.0 * test::m4a::kPi;
    CHECK(std::fabs(volume.value() - expected_volume) <= 1e-12 * expected_volume);

    // The south-west lune sits on the negative-u side of the sphere: after the
    // gap is sewn its equator coedge belongs to the kept north circle, and the
    // chart must unwrap to [-pi, 0] instead of [0, 2*pi].
    const std::vector<FaceId> faces = brep.face_ids();
    REQUIRE(faces.size() == 4u);
    const auto southwest = face_loops_uv(brep, faces[3], 1e-9);
    REQUIRE(southwest.has_value());
    REQUIRE(southwest.value().size() == 1u);
    const UvBounds bounds = uv_bounds(southwest.value().front().points);
    CHECK(bounds.u_min >= -test::m4a::kPi - 1e-9);
    CHECK(bounds.u_max <= 1e-9);
  }

  SUBCASE("cone") {
    // The upper cone band ends exactly at the apex, so this subcase exercises
    // the apex fallback of face_uv_project: the apex projection is undefined
    // from its own coordinates and reuses the previous u.
    constexpr double cone_gap = 1e-7;
    const test::m4a::CurvedSolid fixture = test::m4a::make_cone_gap(cone_gap);
    const auto healed = heal_brep(fixture.brep, fixture.solid, {});
    REQUIRE(healed.has_value());
    CHECK(healed.value().report.changed);
    CHECK(healed.value().report.issues.empty());
    const Brep& brep = healed.value().brep;

    std::size_t rebuilt = 0;
    for (const FaceId face_id : brep.face_ids()) {
      bool any = false;
      bool all = true;
      for (const CoedgeId coedge_id : brep.coedge_ids()) {
        if (brep.coedge(coedge_id).face != face_id) continue;
        const bool has = brep.coedge(coedge_id).pcurve != kNoGeometry;
        any = any || has;
        all = all && has;
      }
      CHECK(any == all);
      if (any) {
        ++rebuilt;
        check_face_pcurves(brep, face_id);
      }
    }
    CHECK(rebuilt == 2u);

    for (const FaceId face_id : brep.face_ids()) {
      const auto loops = face_loops_uv(brep, face_id, 1e-9);
      REQUIRE(loops.has_value());
      REQUIRE(loops.value().size() == 1u);
      const auto* cone = std::get_if<Cone3>(
          &brep.geometry().surfaces[brep.face(face_id).surface]);
      if (cone == nullptr) continue;
      // The apex fan stays on the upper cone chart: every projected point is
      // inside the cone parameter box v in [0, height].
      for (const Point2d& uv : loops.value().front().points) {
        CHECK(uv[1] >= -1e-9);
        CHECK(uv[1] <= 2.0 + 1e-9);
      }
      // The seam runs up to the apex and back down the same column, so the
      // stored pcurve polygon closes as the full 2*pi x band rectangle. The
      // chart reader must keep both apex corners as well, or it reports only
      // the cut pi-wide chart.
      const std::vector<Point2d> pcurve_points = pcurve_loop_points(brep, face_id);
      CHECK(pcurve_points.size() >= 3u);
      CHECK(std::fabs(shoelace_area(pcurve_points) - test::m4a::kTwoPi) <= 1e-9);
      const double chart_area = shoelace_area(loops.value().front().points);
      CHECK(std::fabs(chart_area - test::m4a::kTwoPi) <= 1e-9);
    }
  }
}

namespace {

// Distance between two fixture vertices, in creation order. The healer reports
// these distances as HealReport::max_closed_gap, so the acceptance checks can
// compare the report against the geometry instead of the nominal gap.
double vertex_gap(const Brep& brep, std::size_t first, std::size_t second) {
  const std::vector<VertexId> ids = brep.vertex_ids();
  REQUIRE(first < ids.size());
  REQUIRE(second < ids.size());
  return norm2(brep.vertex(ids[first]).point - brep.vertex(ids[second]).point);
}

struct CurvedAcceptance {
  Brep brep;
  double volume = 0.0;
};

// Shared acceptance gate for a closed curved solid: the report must describe a
// healed, valid, chi == 2 solid whose volume matches the analytic value, and a
// second heal of the input must be deterministic while a heal of the result is
// a no-op.
CurvedAcceptance check_closed_curved_heal(
    const test::m4a::CurvedSolid& fixture, const HealOptions& options,
    double expected_volume, double relative_tolerance, double measured_gap) {
  const std::vector<double> input_fingerprint = brep_fingerprint(fixture.brep);

  const auto healed = heal_brep(fixture.brep, fixture.solid, options);
  if (!healed.has_value()) {
    MESSAGE("heal_brep rejected the closed solid: "
            << std::string(healed.error().message));
  }
  REQUIRE(healed.has_value());
  const HealResult& result = healed.value();
  CHECK(result.report.changed);
  CHECK(result.report.issues.empty());
  CHECK(result.report.boundary_edges == 4u);
  CHECK(result.report.welded_vertices == 2u);
  CHECK(result.report.sewn_edges == 2u);
  CHECK(result.report.max_closed_gap == measured_gap);
  CHECK(result.report.max_tolerance <= 1e-6);
  CHECK(brep_fingerprint(fixture.brep) == input_fingerprint);
  CHECK(euler_characteristic(result.brep, fixture.solid) == 2);
  CHECK(validate_topology(result.brep).ok());

  const auto volume =
      detail::curved_solid_volume(result.brep, fixture.solid, 1e-9);
  REQUIRE(volume.has_value());
  CHECK(std::fabs(volume.value() - expected_volume) <=
        relative_tolerance * expected_volume);

  // Repeating the heal on the same input must produce the same brep bitwise
  // (counts and every vertex position) and the identical volume.
  const auto repeated = heal_brep(fixture.brep, fixture.solid, options);
  REQUIRE(repeated.has_value());
  CHECK(brep_fingerprint(repeated.value().brep) ==
        brep_fingerprint(result.brep));
  CHECK(repeated.value().report.max_closed_gap == result.report.max_closed_gap);
  const auto repeated_volume =
      detail::curved_solid_volume(repeated.value().brep, fixture.solid, 1e-9);
  REQUIRE(repeated_volume.has_value());
  CHECK(repeated_volume.value() == volume.value());

  // A closed result has no boundary left to stitch: the second heal is a
  // no-op and leaves the geometry bitwise unchanged.
  const auto second = heal_brep(result.brep, fixture.solid, options);
  REQUIRE(second.has_value());
  CHECK_FALSE(second.value().report.changed);
  CHECK(second.value().report.issues.empty());
  CHECK(second.value().report.welded_vertices == 0u);
  CHECK(second.value().report.sewn_edges == 0u);
  CHECK(brep_fingerprint(second.value().brep) == brep_fingerprint(result.brep));
  const auto second_volume =
      detail::curved_solid_volume(second.value().brep, fixture.solid, 1e-9);
  REQUIRE(second_volume.has_value());
  CHECK(second_volume.value() == volume.value());
  CHECK(validate_topology(second.value().brep).ok());

  return CurvedAcceptance{result.brep, volume.value()};
}

// Acceptance gate for the two closed-rim torus halves. Every cross-gap rim is
// split into two arcs here, so the report has 8 boundary edges, 4 welded
// vertices and 4 sewn edges. A torus has Euler characteristic 0, not 2: after
// welding, V=4, E=6, F=2. The rest of the contract matches
// check_closed_curved_heal: deterministic repeat, unchanged input and an
// idempotent second heal.
CurvedAcceptance check_torus_halves_heal(
    const test::m4a::CurvedSolid& fixture, const HealOptions& options,
    double expected_volume, double relative_tolerance, double measured_gap) {
  const std::vector<double> input_fingerprint = brep_fingerprint(fixture.brep);

  const auto healed = heal_brep(fixture.brep, fixture.solid, options);
  if (!healed.has_value()) {
    MESSAGE("heal_brep rejected the torus halves: "
            << std::string(healed.error().message));
  }
  REQUIRE(healed.has_value());
  const HealResult& result = healed.value();
  CHECK(result.report.changed);
  CHECK(result.report.issues.empty());
  CHECK(result.report.boundary_edges == 8u);
  CHECK(result.report.welded_vertices == 4u);
  CHECK(result.report.sewn_edges == 4u);
  CHECK(result.report.max_closed_gap == measured_gap);
  CHECK(result.report.max_tolerance <= 1e-6);
  CHECK(brep_fingerprint(fixture.brep) == input_fingerprint);
  // A torus is not a sphere: V=4, E=6, F=2 gives chi = 0. The brief's chi == 2
  // was copied from the closed-rim cylinder; the healer's Euler guard compares
  // against the post-split working topology, which is also 0 here.
  CHECK(euler_characteristic(result.brep, fixture.solid) == 0);
  CHECK(validate_topology(result.brep).ok());

  const auto volume =
      detail::curved_solid_volume(result.brep, fixture.solid, 1e-9);
  REQUIRE(volume.has_value());
  CHECK(std::fabs(volume.value() - expected_volume) <=
        relative_tolerance * expected_volume);

  // Repeating the heal on the same input must produce the same brep bitwise
  // (counts and every vertex position) and the identical volume.
  const auto repeated = heal_brep(fixture.brep, fixture.solid, options);
  REQUIRE(repeated.has_value());
  CHECK(brep_fingerprint(repeated.value().brep) ==
        brep_fingerprint(result.brep));
  CHECK(repeated.value().report.max_closed_gap == result.report.max_closed_gap);
  const auto repeated_volume =
      detail::curved_solid_volume(repeated.value().brep, fixture.solid, 1e-9);
  REQUIRE(repeated_volume.has_value());
  CHECK(repeated_volume.value() == volume.value());

  // A closed result has no boundary left to stitch: the second heal is a
  // no-op and leaves the geometry bitwise unchanged.
  const auto second = heal_brep(result.brep, fixture.solid, options);
  REQUIRE(second.has_value());
  CHECK_FALSE(second.value().report.changed);
  CHECK(second.value().report.issues.empty());
  CHECK(second.value().report.welded_vertices == 0u);
  CHECK(second.value().report.sewn_edges == 0u);
  CHECK(brep_fingerprint(second.value().brep) == brep_fingerprint(result.brep));
  const auto second_volume =
      detail::curved_solid_volume(second.value().brep, fixture.solid, 1e-9);
  REQUIRE(second_volume.has_value());
  CHECK(second_volume.value() == volume.value());
  CHECK(validate_topology(second.value().brep).ok());

  return CurvedAcceptance{result.brep, volume.value()};
}

}  // namespace

TEST_CASE("heal: closed rim gaps are split and healed") {
  constexpr double gap = 1e-6;
  const test::m4a::CurvedSolid fixture =
      test::m4a::make_cylinder_closed_rims_gap(gap);
  // Vertices are created b0, t0, u0, v0; each closed rim is split at its
  // opposite-seam midpoint, then t0/t1 are welded to u0/u1 across the gap.
  const double measured_gap = vertex_gap(fixture.brep, 1, 2);

  const CurvedAcceptance accepted = check_closed_curved_heal(
      fixture, {}, test::m4a::kPi * (2.0 + gap), 1e-12, measured_gap);

  // Healed topology observed on the Release build: V=4, E=6, F=4 (chi=2).
  // Only the two cross-gap rims are boundary edges, so the split pass adds one
  // vertex+edge to each (working copy V=6, E=8); the two cross-gap welds then
  // drop those two vertices and merging the paired halves drops two edges.
  // The bottom/top rings are interior (lateral + cap) and are not split.
  CHECK(accepted.brep.vertex_count() == 4u);
  CHECK(accepted.brep.edge_count() == 6u);
  CHECK(accepted.brep.face_count() == 4u);

  // A successful closed-rim heal leaves one closed shell and no boundary edge.
  const auto shells = accepted.brep.solid_shells(fixture.solid);
  REQUIRE(shells.size() == 1u);
  CHECK(accepted.brep.shell(shells.front()).closed);
  CHECK(detail::heal_boundary_edges(accepted.brep, fixture.solid).empty());

  // The healer keeps the input height 2 + gap instead of shrinking the gap, as
  // in the open two-band case.
  CHECK(std::fabs(accepted.volume - test::m4a::kTwoPi) <=
        1e-6 * test::m4a::kTwoPi);
  CHECK(std::fabs(accepted.volume - 6.2831884487723997) <=
        1e-12 * 6.2831884487723997);
}

TEST_CASE("heal: a rejected closed-rim heal returns the caller's input") {
  // Both rejection classes split the two cross-gap rim edges on a working copy
  // first (post-split boundary 4 from working V=6/E=8), then reject before any
  // SolidEditor mutation. The returned brep must be the caller's input bit for
  // bit, while report.boundary_edges names the post-split working count (4),
  // not the input's two closed rim edges.
  constexpr double gap = 1e-6;
  const test::m4a::CurvedSolid fixture =
      test::m4a::make_cylinder_closed_rims_gap(gap);
  const std::vector<double> input = full_brep_fingerprint(fixture.brep);

  SUBCASE("max_gap below the rim seam") {
    HealOptions options;
    options.max_gap = 1e-12;
    const auto healed = heal_brep(fixture.brep, fixture.solid, options);
    REQUIRE(healed.has_value());
    const HealResult& result = healed.value();
    CHECK_FALSE(result.report.changed);
    CHECK(result.report.boundary_edges == 4u);
    CHECK(result.report.welded_vertices == 0u);
    CHECK(result.report.sewn_edges == 0u);
    REQUIRE(result.report.issues.size() == 1u);
    CHECK(result.report.issues.front().kind == HealIssueKind::gap_too_large);
    CHECK(full_brep_fingerprint(result.brep) == input);
    CHECK(full_brep_fingerprint(fixture.brep) == input);
  }

  SUBCASE("max_tolerance below the rim seam") {
    HealOptions options;
    options.max_tolerance = 1e-9;
    const auto healed = heal_brep(fixture.brep, fixture.solid, options);
    REQUIRE(healed.has_value());
    const HealResult& result = healed.value();
    CHECK_FALSE(result.report.changed);
    CHECK(result.report.boundary_edges == 4u);
    CHECK(result.report.welded_vertices == 0u);
    CHECK(result.report.sewn_edges == 0u);
    REQUIRE(result.report.issues.size() == 1u);
    CHECK(result.report.issues.front().kind == HealIssueKind::invalid_result);
    CHECK(std::string_view(result.report.issues.front().message) ==
          "heal_brep: tolerance cap exceeded");
    CHECK(full_brep_fingerprint(result.brep) == input);
    CHECK(full_brep_fingerprint(fixture.brep) == input);
  }
}

TEST_CASE("heal: analytic curved gaps are accepted") {
  SUBCASE("cylinder") {
    constexpr double gap = 1e-6;
    const test::m4a::CurvedSolid fixture =
        test::m4a::make_two_band_cylinder(1.0, gap);
    // Vertices are created b0, b1, t0, t1, u0, u1, v0, v1; the sewn pair joins
    // t0/ t1 to u0/u1.
    const double measured_gap = vertex_gap(fixture.brep, 2, 4);

    const CurvedAcceptance accepted = check_closed_curved_heal(
        fixture, {}, test::m4a::kPi * (2.0 + gap), 1e-12, measured_gap);

    // The healer repeats the input height 2+gap instead of shrinking the gap,
    // which is the 5e-7 relative difference from the ideal two-unit cylinder.
    CHECK(std::fabs(accepted.volume - test::m4a::kTwoPi) <=
          1e-6 * test::m4a::kTwoPi);
    // The brief's decimal is the measured quadrature value; it equals the
    // double product kPi*(2+gap) to 2.5e-14 relative.
    CHECK(std::fabs(accepted.volume - 6.2831884487723997) <=
          1e-12 * 6.2831884487723997);
  }

  SUBCASE("sphere") {
    // The fixture puts the two boundary circles at v = +/-eps, so the equator
    // vertex separation is 2*eps = 1e-6.
    constexpr double eps = 5e-7;
    const test::m4a::CurvedSolid fixture = test::m4a::make_sphere_gap(eps);
    // Vertices are created n0, n1, np, s0, s1, sp; the sewn pair joins n0/s0
    // and n1/s1.
    const double measured_gap = vertex_gap(fixture.brep, 0, 3);

    const CurvedAcceptance accepted = check_closed_curved_heal(
        fixture, {}, 4.0 / 3.0 * test::m4a::kPi, 1e-12, measured_gap);
    REQUIRE(accepted.brep.face_count() == 4u);

    // Each lune chart is the uv rectangle u in [0, pi], v in [eps, pi/2] (or
    // the mirrored south band). Its area differs from pi^2/2 by +/- pi*eps; the
    // north and south shifts cancel, so the uv sum is 2*pi^2 (the 4*pi*asin(eps)
    // sliver is an equatorial surface-area deficit, not a uv-area deficit).
    double uv_area_sum = 0.0;
    for (const FaceId face_id : accepted.brep.face_ids()) {
      const auto loops = face_loops_uv(accepted.brep, face_id, 1e-9);
      REQUIRE(loops.has_value());
      REQUIRE(loops.value().size() == 1u);
      const double lune_area = shoelace_area(loops.value().front().points);
      CHECK(lune_area > 0.0);
      CHECK(std::fabs(lune_area - 0.5 * test::m4a::kPi * test::m4a::kPi) <=
            1e-5);
      uv_area_sum += lune_area;
    }
    CHECK(std::fabs(uv_area_sum - 2.0 * test::m4a::kPi * test::m4a::kPi) <=
          1e-5);
  }

  SUBCASE("cone") {
    // The band-to-apex fixture's slant gap is gap*sqrt(1.25) = 8.944e-7, which
    // stays inside the default 1e-6 max_gap and max_tolerance caps.
    constexpr double gap = 8e-7;
    const test::m4a::CurvedSolid fixture = test::m4a::make_cone_gap(gap);
    // Vertices are created b0, b1, t0, t1, u0, u1, apex; the sewn pair joins
    // t0/t1 to u0/u1.
    const double measured_gap = vertex_gap(fixture.brep, 2, 4);

    const CurvedAcceptance accepted = check_closed_curved_heal(
        fixture, {}, 2.0943951023931864, 1e-12, measured_gap);
    CHECK(std::fabs(accepted.volume -
                    2.0 / 3.0 * test::m4a::kPi) <=
          1e-12 * (2.0 / 3.0 * test::m4a::kPi));
    REQUIRE(accepted.brep.face_count() == 3u);

    // Seven input vertices minus the two welded rim vertices; the apex is one
    // vertex, never split into a fan.
    CHECK(accepted.brep.vertex_count() == 5u);
    std::size_t apex_vertices = 0;
    for (const VertexId vertex_id : accepted.brep.vertex_ids()) {
      const Point3d point = accepted.brep.vertex(vertex_id).point;
      if (norm2(point - Point3d{0.0, 0.0, 2.0}) <= 1e-12) {
        ++apex_vertices;
      }
    }
    CHECK(apex_vertices == 1u);
  }
}

TEST_CASE("heal: a sphere pole gap closes at the pole") {
  // The cap ends exactly at the pole (v = pi/2) and its meridian seam is used
  // twice; the band's rim sits `delta` below the pole on the same unit sphere,
  // so the solid is the upper hemisphere and its volume is 2*pi/3. Both rims
  // are already split at u = 0 and u = pi, so no closed-rim split is involved.
  // The last pair needs an explicit max_gap: its 5e-5 material gap is larger
  // than the default 1e-6, and rejecting it would be correct behavior.
  const auto run_case = [](double delta, double gap, HealOptions options,
                           double tolerance_bound) {
    CAPTURE(delta);
    CAPTURE(gap);
    const test::m4a::CurvedSolid fixture =
        test::m4b::make_sphere_pole_gap(delta, gap);
    // Vertices are created e0, e1, p0, p1, c0, c1, pole; the sewn pair joins
    // p0/c0 and p1/c1 across the gap.
    const double measured_gap = vertex_gap(fixture.brep, 2, 4);
    const std::vector<double> input_fingerprint = brep_fingerprint(fixture.brep);

    const auto healed = heal_brep(fixture.brep, fixture.solid, options);
    if (!healed.has_value()) {
      MESSAGE("heal_brep rejected the near-pole sphere gap: "
              << std::string(healed.error().message));
    }
    REQUIRE(healed.has_value());
    const HealResult& result = healed.value();
    CHECK(result.report.changed);
    CHECK(result.report.issues.empty());
    CHECK(result.report.boundary_edges == 4u);
    CHECK(result.report.welded_vertices == 2u);
    CHECK(result.report.sewn_edges == 2u);
    CHECK(result.report.max_closed_gap == measured_gap);
    CHECK(result.report.max_tolerance <= tolerance_bound);
    CHECK(brep_fingerprint(fixture.brep) == input_fingerprint);
    CHECK(euler_characteristic(result.brep, fixture.solid) == 2);
    CHECK(validate_topology(result.brep).ok());

    const auto volume =
        detail::curved_solid_volume(result.brep, fixture.solid, 1e-9);
    if (!volume.has_value()) {
      MESSAGE("the healed sphere pole gap has no volume: "
              << std::string(volume.error().message));
    }
    REQUIRE(volume.has_value());
    CHECK(std::fabs(volume.value() - 2.0 / 3.0 * test::m4a::kPi) <=
          1e-12 * (2.0 / 3.0 * test::m4a::kPi));

    // A closed result is one shell with no boundary edge left.
    const auto shells = result.brep.solid_shells(fixture.solid);
    REQUIRE(shells.size() == 1u);
    CHECK(result.brep.shell(shells.front()).closed);
    CHECK(detail::heal_boundary_edges(result.brep, fixture.solid).empty());

    // Repeating the heal on the same input is deterministic.
    const auto repeated = heal_brep(fixture.brep, fixture.solid, options);
    REQUIRE(repeated.has_value());
    CHECK(brep_fingerprint(repeated.value().brep) ==
          brep_fingerprint(result.brep));
    const auto repeated_volume =
        detail::curved_solid_volume(repeated.value().brep, fixture.solid, 1e-9);
    REQUIRE(repeated_volume.has_value());
    CHECK(repeated_volume.value() == volume.value());
  };

  SUBCASE("the rim is one nanometre below the pole") {
    run_case(1e-8, 5e-9, HealOptions{}, 1e-6);
  }
  SUBCASE("the rim is fifty nanometres below the pole") {
    run_case(5e-8, 2.5e-8, HealOptions{}, 1e-6);
  }
  SUBCASE("a hundred-micrometre rim gap needs an explicit max_gap") {
    run_case(1e-4, 5e-5, HealOptions{1e-9, 1e-4, 1e-4}, 1e-4);
  }
}

TEST_CASE("heal: torus band gaps are accepted") {
  // M4b-1 acceptance for torus bands whose gap runs along u (around the axis)
  // or v (around the tube), including a band pair that straddles the torus v
  // seam. The gap becomes material, so each healed volume is the volume of the
  // un-gapped shape: Pappus for the sector, a meridian integral for the barrels.
  // Every reference below is the brief's quadrature decimal; each differs from
  // the closed form by at most 7e-15 relative, far inside the 1e-9 pin.
  SUBCASE("u gap on a 3 rad sector") {
    constexpr double gap = 1e-6;
    // Closed form: V = theta * R * pi * r^2 = 3*pi/4 = 2.356194490192345...
    constexpr double expected = 2.3561944901923608;
    const test::m4a::CurvedSolid fixture = test::m4a::make_torus_u_gap(gap);
    // Vertices are created a0,a1,b0,b1,c0,c1,d0,d1: b0/c0 are the v = 0 ends of
    // the two open rims; their chord is exactly `gap` at the outer equator while
    // the v = pi ends sit (R - r)/(R + r) times closer.
    const double measured_gap = vertex_gap(fixture.brep, 2, 4);

    const CurvedAcceptance accepted =
        check_closed_curved_heal(fixture, {}, expected, 1e-9, measured_gap);
    // The two rim welds pull the upper band's rim onto the lower one at
    // u = 1.5 - delta/2, so the sector spans exactly [0, 3].
    CHECK(accepted.brep.vertex_count() == 6u);
    CHECK(accepted.brep.edge_count() == 8u);
    CHECK(accepted.brep.face_count() == 4u);
    const auto shells = accepted.brep.solid_shells(fixture.solid);
    REQUIRE(shells.size() == 1u);
    CHECK(accepted.brep.shell(shells.front()).closed);
    CHECK(detail::heal_boundary_edges(accepted.brep, fixture.solid).empty());
  }

  SUBCASE("outer branch v gap") {
    constexpr double gap = 1e-6;
    // Closed form, rho(z) = R + sqrt(r^2 - z^2), c = r*sin(1):
    //   V = pi * [2*R^2*c + 2*R*(c*s + r^2*asin(c/r)) + 2*r^2*c - (2/3)*c^3]
    // with s = sqrt(r^2 - c^2) and asin(c/r) = 1.
    constexpr double expected = 5.4334195273316128;
    const test::m4a::CurvedSolid fixture = test::m4a::make_torus_outer_v_gap(gap);
    // Vertices are created a0,a1,b0,b1,c0,c1,d0,d1: b0/c0 are the u = 0 ends of
    // the two open rims, 2*r*sin(delta/2) = gap apart.
    const double measured_gap = vertex_gap(fixture.brep, 2, 4);

    const CurvedAcceptance accepted =
        check_closed_curved_heal(fixture, {}, expected, 1e-9, measured_gap);
    CHECK(accepted.brep.vertex_count() == 6u);
    CHECK(accepted.brep.edge_count() == 8u);
    CHECK(accepted.brep.face_count() == 4u);
    const auto shells = accepted.brep.solid_shells(fixture.solid);
    REQUIRE(shells.size() == 1u);
    CHECK(accepted.brep.shell(shells.front()).closed);
    CHECK(detail::heal_boundary_edges(accepted.brep, fixture.solid).empty());
  }

  SUBCASE("inner branch band across the v seam") {
    constexpr double gap = 1e-6;
    // Closed form, rho(z) = R - sqrt(r^2 - z^2), c = r*sin(1):
    //   V = pi * [2*R^2*c - 2*R*(c*s + r^2*asin(c/r)) + 2*r^2*c - (2/3)*c^3].
    constexpr double expected = 0.8635058157199863;
    const test::m4a::CurvedSolid fixture =
        test::m4a::make_torus_inner_v_seam_gap(gap);
    // Vertices are created a0,a1,b0,b1,c0,c1,d0,d1: b0/c0 are the u = 0 ends of
    // the two open rims. The chart v of the two halves is on opposite sides of
    // the surface's chart seam (v_origin = pi), so the projection puts one half
    // near v = 0 and the other near v = -2*pi.
    const double measured_gap = vertex_gap(fixture.brep, 2, 4);

    const CurvedAcceptance accepted =
        check_closed_curved_heal(fixture, {}, expected, 1e-9, measured_gap);
    CHECK(accepted.brep.vertex_count() == 6u);
    CHECK(accepted.brep.edge_count() == 8u);
    CHECK(accepted.brep.face_count() == 4u);
    const auto shells = accepted.brep.solid_shells(fixture.solid);
    REQUIRE(shells.size() == 1u);
    CHECK(accepted.brep.shell(shells.front()).closed);
    CHECK(detail::heal_boundary_edges(accepted.brep, fixture.solid).empty());

    // A pcurve rebuild that let the chart seam wrap mid-face would still give
    // the right volume (the torus is periodic in v), so pin the chart window
    // per face instead: every healed face spans less than a half turn in v.
    for (const FaceId face_id : accepted.brep.face_ids()) {
      const auto loops = face_loops_uv(accepted.brep, face_id, 1e-9);
      REQUIRE(loops.has_value());
      double v_min = 0.0;
      double v_max = 0.0;
      bool first = true;
      for (const FaceUvLoop& loop : loops.value()) {
        for (const Point2d& point : loop.points) {
          if (first) {
            v_min = point[1];
            v_max = point[1];
            first = false;
          } else {
            v_min = std::min(v_min, point[1]);
            v_max = std::max(v_max, point[1]);
          }
        }
      }
      CHECK(v_max - v_min < 2.0);
    }
  }
}

TEST_CASE("heal: torus halves with closed rims are split, welded and sewn") {
  constexpr double gap = 1e-6;
  const test::m4a::CurvedSolid fixture =
      test::m4a::make_torus_closed_rims_gap(gap);
  // Vertices are created outer-high, inner-high, inner-low, outer-low. Each
  // closed rim is split at u = pi before planning, so vertex_gap(0, 3) is the
  // outer cross-gap pair and is the report's maximum closed gap.
  const double measured_gap = vertex_gap(fixture.brep, 0, 3);

  // The two halves are glued across the cut planes, so the healed solid is the
  // whole torus: 2*pi^2*R*r^2 = pi^2/2.
  constexpr double expected_volume = 4.934802200544679;
  const CurvedAcceptance accepted = check_torus_halves_heal(
      fixture, {}, expected_volume, 1e-12, measured_gap);

  // Healed topology: V=4, E=6, F=2 (chi=0). The four closed rims split into
  // eight boundary arcs (working V=8, E=10); the four cross-gap pairs then weld
  // the lower arc endpoints onto the upper ones and merge each pair, dropping
  // four vertices and four edges.
  CHECK(accepted.brep.vertex_count() == 4u);
  CHECK(accepted.brep.edge_count() == 6u);
  CHECK(accepted.brep.face_count() == 2u);
  CHECK(euler_characteristic(accepted.brep, fixture.solid) == 0);

  const auto shells = accepted.brep.solid_shells(fixture.solid);
  REQUIRE(shells.size() == 1u);
  CHECK(accepted.brep.shell(shells.front()).closed);
  CHECK(detail::heal_boundary_edges(accepted.brep, fixture.solid).empty());

  // The split reuses freed edge handles, so the kept cross-gap rims are the
  // lower outer rim and the upper inner rim. The healed upper face therefore
  // spans v in [-alpha, pi-alpha] and the lower face spans [pi-alpha,
  // 2*pi-alpha]; together they cover exactly one turn, which is why the volume
  // is the full torus despite the chart window shifting by alpha.
  const double alpha = std::asin(0.5 * gap * (1.0 - 1e-8) /
                                 test::m4a::torus_fixture::kMinorRadius);
  constexpr double pi = std::numbers::pi_v<double>;
  const std::array<std::pair<double, double>, 2> expected_windows{{
      {-alpha, pi - alpha}, {pi - alpha, 2.0 * pi - alpha}}};
  std::size_t face_index = 0;
  for (const FaceId face_id : accepted.brep.face_ids()) {
    REQUIRE(face_index < expected_windows.size());
    const auto loops = face_loops_uv(accepted.brep, face_id, 1e-9);
    REQUIRE(loops.has_value());
    double v_min = 0.0;
    double v_max = 0.0;
    bool first = true;
    for (const FaceUvLoop& loop : loops.value()) {
      for (const Point2d& point : loop.points) {
        if (first) {
          v_min = point[1];
          v_max = point[1];
          first = false;
        } else {
          v_min = std::min(v_min, point[1]);
          v_max = std::max(v_max, point[1]);
        }
      }
    }
    CHECK(std::fabs(v_min - expected_windows[face_index].first) <= 1e-9);
    CHECK(std::fabs(v_max - expected_windows[face_index].second) <= 1e-9);
    ++face_index;
  }
  CHECK(face_index == 2u);
}

TEST_CASE("heal: torus closed rims beyond max_gap report gap_too_large") {
  // The same closed-rim halves cut 1e-3 apart: every rim is split, but no pair
  // fits inside the 1e-6 max_gap, so the healer reports the nearest miss.
  constexpr double gap = 1e-3;
  const test::m4a::CurvedSolid fixture =
      test::m4a::make_torus_closed_rims_gap(gap);

  const auto healed = heal_brep(fixture.brep, fixture.solid, {});
  REQUIRE(healed.has_value());
  const HealResult& result = healed.value();
  CHECK_FALSE(result.report.changed);
  CHECK(result.report.boundary_edges == 8u);
  CHECK(result.report.welded_vertices == 0u);
  CHECK(result.report.sewn_edges == 0u);
  REQUIRE(!result.report.issues.empty());
  const HealIssue& issue = result.report.issues.front();
  CHECK(issue.kind == HealIssueKind::gap_too_large);
  // The split reuses the freed closed-rim edge handles, so the nearest miss is
  // the inner rim pair: one half of the upper inner rim and one half of the
  // lower inner rim.
  CHECK(issue.entity_a == 0u);
  CHECK(issue.entity_b == 1u);
  CHECK(issue.gap > 1e-6);
}

TEST_CASE("heal: a torus rim against a plane lid reports a carrier mismatch") {
  // The upper torus half is capped by a plane disk whose rim is the half's
  // outer closed rim. The rims pair at zero gap, so the first pair the surface
  // guard sees joins a Torus3 face to a Plane3 face.
  const test::m4a::CurvedSolid fixture =
      test::m4a::make_torus_plane_lid_gap(1e-6);
  const std::vector<double> before = brep_fingerprint(fixture.brep);

  const auto healed = heal_brep(fixture.brep, fixture.solid, {});
  REQUIRE_FALSE(healed.has_value());
  CHECK(healed.error().code == Errc::unsupported);
  CHECK(std::string_view(healed.error().message) == kCarrierMessage);
  // The closed-rim split reuses freed edge handles: edge 1 is one plane-disk
  // arc (created after the torus outer rim was split) and edge 4 is the
  // coincident torus outer arc it was paired with.
  CHECK(healed.error().context_a == 1u);
  CHECK(healed.error().context_b == 4u);
  CHECK(brep_fingerprint(fixture.brep) == before);
}

TEST_CASE("heal: an off-surface torus rim reports the off-surface error") {
  // The lower outer rim is pushed 1e-7 radially outward, so its circle is not
  // on the shared Torus3 surface even though the cross-gap chord still fits
  // inside max_gap and planning pairs the halves.
  const test::m4a::CurvedSolid fixture =
      test::m4a::make_torus_off_surface_gap(1e-6, 1e-7);
  const std::vector<double> before = brep_fingerprint(fixture.brep);

  const auto healed = heal_brep(fixture.brep, fixture.solid, {});
  REQUIRE_FALSE(healed.has_value());
  CHECK(healed.error().code == Errc::unsupported);
  CHECK(std::string_view(healed.error().message) == kOffSurfaceMessage);
  // Handle reuse after the closed-rim split makes the kept edge (2) one shifted
  // lower outer arc and the dropped edge (6) its upper counterpart.
  CHECK(healed.error().context_a == 2u);
  CHECK(healed.error().context_b == 6u);
  CHECK(brep_fingerprint(fixture.brep) == before);
}

TEST_CASE("heal: an inherited carrier must still pass through its welded vertices") {
  const test::m4a::CurvedSolid fixture =
      test::m4a::make_two_band_cylinder(1.0, 1e-6);
  SolidEditor editor(fixture.brep);
  // Edges 2 and 3 are the lower rims (lower_top0/lower_top1); the fixture
  // allocates them in creation order.
  const std::vector<EdgeId> edge_ids = fixture.brep.edge_ids();
  REQUIRE(edge_ids.size() == 10u);
  const Result<void> cleared_lower =
      editor.set_edge_geometry(edge_ids[2], kNoGeometry, 0.0, 1.0);
  const Result<void> cleared_upper =
      editor.set_edge_geometry(edge_ids[3], kNoGeometry, 0.0, 1.0);
  REQUIRE(cleared_lower.has_value());
  REQUIRE(cleared_upper.has_value());
  const Result<Brep> committed = editor.commit();
  REQUIRE(committed.has_value());

  // merge_edges hands the dropped rim's circle to the curve-less keep rim, but
  // the inherited circle sits at the dropped rim's height: it no longer passes
  // through the keep vertices and the heal must reject instead of storing a
  // pcurve that leaves the chart.
  const auto healed = heal_brep(committed.value(), fixture.solid, {});
  REQUIRE_FALSE(healed.has_value());
  CHECK(healed.error().code == Errc::unsupported);
  CHECK(std::string_view(healed.error().message) ==
        "heal_brep: kept edge curve does not match its welded vertices");
}

TEST_CASE("heal: a reversed inherited carrier keeps the edge direction") {
  Brep brep;
  SolidEditor editor(brep);
  const Result<VertexId> start = editor.add_vertex(Point3d{1.0, 0.0, 0.0}, 0.0);
  const Result<VertexId> end = editor.add_vertex(Point3d{-1.0, 0.0, 0.0}, 0.0);
  REQUIRE(start.has_value());
  REQUIRE(end.has_value());
  CircleCurve3 circle;
  circle.center = Point3d{0.0, 0.0, 0.0};
  circle.radius = 1.0;
  const std::uint32_t curve = editor.add_curve(circle);
  const Result<EdgeId> edge = editor.add_edge(start.value(), end.value(), curve,
                                              0.0, test::m4a::kPi, 0.0);
  REQUIRE(edge.has_value());

  const Result<std::uint32_t> reversed =
      detail::heal_reverse_curve(editor, edge.value());
  REQUIRE(reversed.has_value());
  const auto& carrier = std::get<CircleCurve3>(
      editor.brep().geometry().curves[reversed.value()]);
  CHECK(norm2(carrier.eval(0.0) - circle.eval(test::m4a::kPi)) <= 1e-12);
  CHECK(norm2(carrier.eval(test::m4a::kPi) - circle.eval(0.0)) <= 1e-12);
  // Interior samples follow the same reversed parameterization, not a mirrored
  // arc: eval_reversed(t) == eval(t0 + t1 - t).
  CHECK(norm2(carrier.eval(0.5) - circle.eval(test::m4a::kPi - 0.5)) <= 1e-12);
}

TEST_CASE("heal: a flipped inherited carrier is reversed into the keep direction") {
  const test::FlippedRimBox box = test::make_flipped_rim_box();

  // The fixture stores the upper rim edges opposite to the lower ones, so the
  // plan must pick the crossed weld mapping for every rim pair.
  const std::vector<EdgeId> boundary =
      detail::heal_boundary_edges(box.brep, box.solid);
  REQUIRE(boundary.size() == 8u);
  const Result<detail::HealPlan> plan =
      detail::heal_plan(box.brep, boundary, 1e-6);
  REQUIRE(plan.has_value());
  REQUIRE(plan.value().pairs.size() == 4u);
  for (const detail::HealEdgePair& pair : plan.value().pairs) {
    CHECK_FALSE(pair.same_direction);
    const Edge& keep_edge = box.brep.edge(pair.keep);
    const Edge& drop_edge = box.brep.edge(pair.drop);
    const Point3d k0 = box.brep.vertex(keep_edge.v0).point;
    const Point3d k1 = box.brep.vertex(keep_edge.v1).point;
    const Point3d d0 = box.brep.vertex(drop_edge.v0).point;
    const Point3d d1 = box.brep.vertex(drop_edge.v1).point;
    CHECK(norm2(k0 - d1) == 0.0);
    CHECK(norm2(k1 - d0) == 0.0);
  }

  const auto healed = heal_brep(box.brep, box.solid, {});
  REQUIRE(healed.has_value());
  const HealResult& healed_result = healed.value();
  CHECK(healed_result.report.changed);
  CHECK(healed_result.report.issues.empty());

  // The inherited line starts at the dropped edge's v0. Without the flip repair
  // that is the keep edge's v1, so eval(t0) would land on the wrong endpoint.
  for (const EdgeId rim : box.lower_rim) {
    const Edge& edge = healed_result.brep.edge(rim);
    CHECK(edge.curve != kNoGeometry);
    if (edge.curve == kNoGeometry) continue;
    const auto* line = std::get_if<LineCurve3>(
        &healed_result.brep.geometry().curves[edge.curve]);
    REQUIRE(line != nullptr);
    CHECK(norm2(line->eval(edge.t0) -
                healed_result.brep.vertex(edge.v0).point) <= 1e-12);
    CHECK(norm2(line->eval(edge.t1) -
                healed_result.brep.vertex(edge.v1).point) <= 1e-12);
    CHECK(std::fabs(norm2(healed_result.brep.vertex(edge.v0).point -
                          healed_result.brep.vertex(edge.v1).point) -
                    1.0) <= 1e-12);
  }

  // The welded box is closed and has the unit volume of the two half boxes.
  const auto volume =
      detail::curved_solid_volume(healed_result.brep, box.solid, 1e-9);
  REQUIRE(volume.has_value());
  CHECK(std::fabs(volume.value() - 1.0) <= 1e-12);
}

TEST_CASE("heal: a cone pinhole closes at a welded point") {
  // The rim radius is zero, so the two split halves start and end on the same
  // point. Pairing cannot repair that rim: after a weld both ends of the gap
  // are one vertex. The point-gap weld pass has to dissolve both halves and
  // merge the two rim vertices before the planner sees a zero-length edge.
  SolidEditor editor;
  const test::m4b::ConePinhole cone =
      test::m4b::add_cone_pinhole(editor, 0.0, 0.0, true);
  const auto committed = editor.commit();
  REQUIRE(committed.has_value());
  const Brep brep = committed.value();
  const std::vector<double> input_fingerprint = brep_fingerprint(brep);
  const HealOptions options{1e-9, 1e-9, 1e-9};

  const auto healed = heal_brep(brep, cone.solid, options);
  if (!healed.has_value()) {
    MESSAGE("heal_brep rejected the zero-radius pinhole: "
            << std::string(healed.error().message));
  }
  REQUIRE(healed.has_value());
  const HealResult& result = healed.value();
  CHECK(result.report.changed);
  CHECK(result.report.issues.empty());
  CHECK(result.report.boundary_edges == 0u);
  CHECK(result.report.welded_vertices == 2u);
  CHECK(result.report.sewn_edges == 0u);
  CHECK(result.report.collapsed_edges == 2u);
  // The pass may not invent a tolerance the caller did not allow.
  CHECK(result.report.max_tolerance <= 1e-9);
  CHECK(brep_fingerprint(brep) == input_fingerprint);

  // V = b0, b1, t0; E = base0, base1, seam; F = lateral, base.
  CHECK(result.brep.vertex_count() == 3u);
  CHECK(result.brep.edge_count() == 3u);
  CHECK(result.brep.face_count() == 2u);
  CHECK(euler_characteristic(result.brep, cone.solid) == 2);
  CHECK(validate_topology(result.brep).ok());

  const auto volume =
      detail::curved_solid_volume(result.brep, cone.solid, 1e-9);
  if (!volume.has_value()) {
    MESSAGE("the point-welded pinhole has no volume: "
            << std::string(volume.error().message));
  }
  REQUIRE(volume.has_value());
  const double expected_volume = 2.0 / 3.0 * test::m4a::kPi;
  CHECK(std::fabs(volume.value() - expected_volume) <=
        1e-12 * expected_volume);

  const auto shells = result.brep.solid_shells(cone.solid);
  REQUIRE(shells.size() == 1u);
  CHECK(result.brep.shell(shells.front()).closed);
  CHECK(detail::heal_boundary_edges(result.brep, cone.solid).empty());

  // Repeating the heal on the same input is deterministic.
  const auto repeated = heal_brep(brep, cone.solid, options);
  REQUIRE(repeated.has_value());
  CHECK(brep_fingerprint(repeated.value().brep) == brep_fingerprint(result.brep));
  const auto repeated_volume =
      detail::curved_solid_volume(repeated.value().brep, cone.solid, 1e-9);
  REQUIRE(repeated_volume.has_value());
  CHECK(repeated_volume.value() == volume.value());
}

TEST_CASE("heal: a tiny pinhole rim welds only within max_gap") {
  // The rim radius is 1e-6, so the split halves are two arcs 2e-6 apart. With
  // max_gap = 1e-5 the point-gap pass may collapse them; with max_gap = 1e-9 it
  // must not, and the caller must be told no pair could be made instead of
  // receiving a silently open shell.
  SolidEditor editor;
  const test::m4b::ConePinhole cone =
      test::m4b::add_cone_pinhole(editor, 0.0, 1e-6, true);
  const auto committed = editor.commit();
  REQUIRE(committed.has_value());
  const Brep brep = committed.value();
  const std::vector<double> input_fingerprint = brep_fingerprint(brep);

  const auto healed = heal_brep(brep, cone.solid, HealOptions{1e-9, 1e-5, 1e-4});
  if (!healed.has_value()) {
    MESSAGE("heal_brep rejected the tiny pinhole: "
            << std::string(healed.error().message));
  }
  REQUIRE(healed.has_value());
  const HealResult& result = healed.value();
  CHECK(result.report.changed);
  CHECK(result.report.issues.empty());
  CHECK(result.report.boundary_edges == 0u);
  CHECK(result.report.welded_vertices == 2u);
  CHECK(result.report.sewn_edges == 0u);
  CHECK(result.report.collapsed_edges == 2u);
  // The merged rim vertices carry the 2e-6 rim diameter as tolerance.
  CHECK(result.report.max_tolerance <= 1e-4);
  CHECK(brep_fingerprint(brep) == input_fingerprint);

  CHECK(result.brep.vertex_count() == 3u);
  CHECK(result.brep.edge_count() == 3u);
  CHECK(result.brep.face_count() == 2u);
  CHECK(euler_characteristic(result.brep, cone.solid) == 2);
  CHECK(validate_topology(result.brep).ok());

  const auto volume =
      detail::curved_solid_volume(result.brep, cone.solid, 1e-9);
  if (!volume.has_value()) {
    MESSAGE("the point-welded tiny pinhole has no volume: "
            << std::string(volume.error().message));
  }
  REQUIRE(volume.has_value());
  const double expected_volume = 2.0 / 3.0 * test::m4a::kPi;
  CHECK(std::fabs(volume.value() - expected_volume) <=
        1e-9 * expected_volume);

  const auto shells = result.brep.solid_shells(cone.solid);
  REQUIRE(shells.size() == 1u);
  CHECK(result.brep.shell(shells.front()).closed);
  CHECK(detail::heal_boundary_edges(result.brep, cone.solid).empty());

  // Below the rim diameter nothing may weld; the solid stays open and the
  // caller gets the honest "no pair" answer instead of a false success.
  const auto rejected =
      heal_brep(brep, cone.solid, HealOptions{1e-9, 1e-9, 1e-4});
  REQUIRE(rejected.has_value());
  const HealResult& rejected_result = rejected.value();
  CHECK_FALSE(rejected_result.report.changed);
  CHECK(rejected_result.report.boundary_edges == 2u);
  CHECK(rejected_result.report.welded_vertices == 0u);
  CHECK(rejected_result.report.sewn_edges == 0u);
  CHECK(rejected_result.report.collapsed_edges == 0u);
  CHECK(brep_fingerprint(rejected_result.brep) == input_fingerprint);
  REQUIRE(rejected_result.report.issues.size() == 1u);
  CHECK(rejected_result.report.issues.front().kind ==
        HealIssueKind::invalid_result);
  CHECK(rejected_result.report.issues.front().message ==
        "heal_brep: no mutual pair");
  CHECK(detail::heal_boundary_edges(rejected_result.brep, cone.solid).size() ==
        2u);
}
