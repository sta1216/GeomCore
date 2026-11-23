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

bool contains(std::string_view text, std::string_view needle) {
  return text.find(needle) != std::string_view::npos;
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
