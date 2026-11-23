#include <doctest/doctest.h>

#include <algorithm>
#include <cstdint>
#include <string_view>
#include <vector>

#include <geom/brep/heal.hpp>

#include "brep_heal_fixture.hpp"
#include "m4a_curved_heal_fixture.hpp"

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

std::size_t count_issues(const detail::HealPlan& plan, HealIssueKind kind) {
  return static_cast<std::size_t>(std::count_if(
      plan.issues.begin(), plan.issues.end(),
      [kind](const HealIssue& issue) { return issue.kind == kind; }));
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
