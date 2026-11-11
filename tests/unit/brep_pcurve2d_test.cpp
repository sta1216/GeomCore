#include <doctest/doctest.h>

#include <array>
#include <cmath>
#include <cstdint>
#include <limits>
#include <string_view>
#include <utility>
#include <vector>

#include <geom/brep/pcurve2d.hpp>
#include <geom/brep/solid_editor.hpp>
#include <geom/brep/validate.hpp>
#include <geom/core/error.hpp>

using namespace geom;

namespace {

bool same_point(const Point2d& a, const Point2d& b) {
  return a[0] == b[0] && a[1] == b[1];
}

}  // namespace

TEST_CASE("pcurve2d stores a polyline") {
  const Pcurve2d empty;
  CHECK(empty.points.empty());

  const Pcurve2d two{{0.0, 0.0}, {1.0, 1.0}};
  REQUIRE(two.points.size() == 2u);
  CHECK(two.points.front()[0] == 0.0);
  CHECK(two.points.front()[1] == 0.0);
  CHECK(two.points.back()[0] == 1.0);
  CHECK(two.points.back()[1] == 1.0);

  const Point2d a{2.0, 3.0};
  const Point2d b{5.0, 7.0};
  const Pcurve2d from_points{a, b};
  REQUIRE(from_points.points.size() == 2u);
  CHECK(same_point(from_points.points.front(), a));
  CHECK(same_point(from_points.points.back(), b));

  const Pcurve2d three{{0.0, 0.0}, {1.0, 0.0}, {1.0, 2.0}};
  REQUIRE(three.points.size() == 3u);
  CHECK(same_point(three.points[1], Point2d{1.0, 0.0}));

  std::vector<Point2d> source{{0.0, 0.0}, {1.0, 1.0}, {2.0, 0.0}};
  const Pcurve2d from_vector(std::move(source));
  CHECK(from_vector.points.size() == 3u);
}

TEST_CASE("pcurve length sums chord lengths") {
  const Pcurve2d three_points{{0.0, 0.0}, {1.0, 0.0}, {1.0, 2.0}};
  CHECK(pcurve_length(three_points) == 3.0);

  const Pcurve2d triangle_leg{{0.0, 0.0}, {3.0, 4.0}};
  CHECK(pcurve_length(triangle_leg) == 5.0);

  const Pcurve2d repeated{{0.0, 0.0}, {0.0, 0.0}, {3.0, 4.0}};
  CHECK(pcurve_length(repeated) == 5.0);

  CHECK(pcurve_length(Pcurve2d{}) == 0.0);
  CHECK(pcurve_length(Pcurve2d{{1.0, 2.0}}) == 0.0);
  CHECK(pcurve_length(Pcurve2d{{0.0, 0.0},
                               {std::numeric_limits<double>::quiet_NaN(), 1.0}}) == 0.0);
  CHECK(pcurve_length(Pcurve2d{{0.0, 0.0},
                               {std::numeric_limits<double>::infinity(), 1.0}}) == 0.0);
}

TEST_CASE("pcurve well formed requires two finite points") {
  CHECK_FALSE(pcurve_well_formed(Pcurve2d{}));
  CHECK_FALSE(pcurve_well_formed(Pcurve2d{{0.0, 0.0}}));
  CHECK_FALSE(pcurve_well_formed(Pcurve2d{{0.0, 0.0},
                                          {std::numeric_limits<double>::quiet_NaN(), 0.0}}));
  CHECK_FALSE(pcurve_well_formed(Pcurve2d{{0.0, 0.0},
                                          {0.0, std::numeric_limits<double>::infinity()}}));
  CHECK(pcurve_well_formed(Pcurve2d{{0.0, 0.0}, {1.0, 1.0}}));
  CHECK(pcurve_well_formed(Pcurve2d{{0.0, 0.0}, {0.0, 0.0}}));
}

TEST_CASE("pcurve split walks chord length and clamps the parameter") {
  const Pcurve2d pc{{0.0, 0.0}, {1.0, 0.0}, {1.0, 2.0}};
  const auto half = pcurve_split(pc, 0.5);
  REQUIRE(half.has_value());
  CHECK(half.value().points.size() == 3u);
  CHECK(same_point(half.value().points.front(), Point2d{0.0, 0.0}));
  CHECK(same_point(half.value().points[1], Point2d{1.0, 0.0}));
  CHECK(half.value().points.back()[0] == 1.0);
  CHECK(half.value().points.back()[1] == doctest::Approx(0.5));

  const auto over = pcurve_split(pc, 4.0);
  REQUIRE(over.has_value());
  CHECK(over.value().points.size() == pc.points.size());
  CHECK(same_point(over.value().points.back(), pc.points.back()));

  const auto under = pcurve_split(pc, -3.0);
  REQUIRE(under.has_value());
  REQUIRE(under.value().points.size() == 2u);
  CHECK(same_point(under.value().points.front(), pc.points.front()));
  CHECK(same_point(under.value().points.back(), pc.points.front()));

  // A split that lands on an existing vertex keeps that vertex exactly once.
  const auto at_vertex = pcurve_split(pc, 1.0 / 3.0);
  REQUIRE(at_vertex.has_value());
  REQUIRE(at_vertex.value().points.size() == 2u);
  CHECK(same_point(at_vertex.value().points.back(), pc.points[1]));
}

TEST_CASE("pcurve split at one keeps the endpoint bit identical") {
  // The endpoint is not reachable by interpolating from the start: a + (b - a)
  // rounds to a different double than b, so a split at t = 1 must snap to the
  // existing vertex instead of interpolating.
  const Pcurve2d pc{{0.90474681086128, 0.0}, {-0.950247999257524, 0.0}};
  const auto split = pcurve_split(pc, 1.0);
  REQUIRE(split.has_value());
  REQUIRE(split.value().points.size() == 2u);
  CHECK(split.value().points.back()[0] == pc.points.back()[0]);
  CHECK(split.value().points.back()[1] == pc.points.back()[1]);
}

TEST_CASE("pcurve split rejects a non-finite total chord length") {
  // Every coordinate is finite, but the chord overflows to +inf, so the split
  // fraction would otherwise be NaN and produce wrong halves.
  const Pcurve2d overflow{{-1e308, -1e308}, {1e308, 1e308}};
  CHECK(pcurve_well_formed(overflow));
  CHECK(std::isinf(pcurve_length(overflow)));
  CHECK(pcurve_length(overflow) > 0.0);

  const auto split = pcurve_split(overflow, 0.5);
  REQUIRE_FALSE(split.has_value());
  CHECK(split.error().code == Errc::non_finite);
  CHECK(std::string_view(split.error().message) == "pcurve_split: chord length is not finite");
}

TEST_CASE("pcurve split rejects malformed input") {
  const auto no_points = pcurve_split(Pcurve2d{}, 0.5);
  REQUIRE_FALSE(no_points.has_value());
  CHECK(no_points.error().code == Errc::degenerate);

  const auto one_point = pcurve_split(Pcurve2d{{0.0, 0.0}}, 0.5);
  REQUIRE_FALSE(one_point.has_value());
  CHECK(one_point.error().code == Errc::degenerate);

  const auto zero_length = pcurve_split(Pcurve2d{{1.0, 1.0}, {1.0, 1.0}}, 0.5);
  REQUIRE_FALSE(zero_length.has_value());
  CHECK(zero_length.error().code == Errc::degenerate);

  const auto nan_point = pcurve_split(
      Pcurve2d{{0.0, 0.0}, {std::numeric_limits<double>::quiet_NaN(), 1.0}}, 0.5);
  REQUIRE_FALSE(nan_point.has_value());
  CHECK(nan_point.error().code == Errc::non_finite);

  const auto inf_point = pcurve_split(
      Pcurve2d{{0.0, 0.0}, {std::numeric_limits<double>::infinity(), 1.0}}, 0.5);
  REQUIRE_FALSE(inf_point.has_value());
  CHECK(inf_point.error().code == Errc::non_finite);

  const Pcurve2d pc{{0.0, 0.0}, {1.0, 1.0}};
  const auto nan_t = pcurve_split(pc, std::numeric_limits<double>::quiet_NaN());
  REQUIRE_FALSE(nan_t.has_value());
  CHECK(nan_t.error().code == Errc::non_finite);
}

TEST_CASE("pcurve merge drops exactly one duplicate junction") {
  const Pcurve2d pc{{0.0, 0.0}, {1.0, 0.0}, {1.0, 2.0}};
  const auto half = pcurve_split(pc, 0.5);
  REQUIRE(half.has_value());
  const Pcurve2d upper{{1.0, 0.5}, {1.0, 2.0}};

  const Pcurve2d merged = pcurve_merge(half.value(), upper);
  REQUIRE(merged.points.size() == 4u);
  CHECK(same_point(merged.points[0], pc.points[0]));
  CHECK(same_point(merged.points[1], pc.points[1]));
  CHECK(same_point(merged.points[2], Point2d{1.0, 0.5}));
  CHECK(same_point(merged.points[2], half.value().points.back()));
  CHECK(same_point(merged.points[3], pc.points[2]));

  // A junction that is not exactly shared is kept on both sides.
  const Pcurve2d shifted{{1.0, 0.6}, {1.0, 2.0}};
  const Pcurve2d joined = pcurve_merge(half.value(), shifted);
  REQUIRE(joined.points.size() == 5u);
  CHECK(same_point(joined.points[2], Point2d{1.0, 0.5}));
  CHECK(same_point(joined.points[3], Point2d{1.0, 0.6}));
}

TEST_CASE("pcurve merge keeps a one ulp junction") {
  const double x = 1.0;
  const double next = std::nextafter(x, std::numeric_limits<double>::infinity());
  REQUIRE(next != x);
  const Pcurve2d lower{{0.0, 0.0}, {x, 1.0}};
  const Pcurve2d upper{{next, 1.0}, {2.0, 0.0}};

  // One ulp apart is not exactly equal, so both junction points are kept.
  const Pcurve2d merged = pcurve_merge(lower, upper);
  REQUIRE(merged.points.size() == lower.points.size() + upper.points.size());
  CHECK(same_point(merged.points[1], Point2d{x, 1.0}));
  CHECK(same_point(merged.points[2], Point2d{next, 1.0}));
}

TEST_CASE("split edge splits each pcurve polyline and join rebuilds it") {
  SolidEditor ed;
  const VertexId v0 = ed.add_vertex({0.0, 0.0, 0.0}, 0.0).value();
  const VertexId v1 = ed.add_vertex({1.0, 0.0, 0.0}, 0.0).value();
  const std::uint32_t line = ed.add_curve(LineCurve3{{0.0, 0.0, 0.0}, {1.0, 0.0, 0.0}});
  const EdgeId edge = ed.add_edge(v0, v1, line, 0.0, 1.0, 0.0).value();
  const std::uint32_t plane = ed.add_surface(Plane3{{0.0, 0.0, 1.0}, 0.0});
  const FaceId face = ed.add_face(plane, false, 0.0).value();
  // A three-point polyline of length three so the split lands mid-segment.
  const Pcurve2d source{{0.0, 0.0}, {1.0, 0.0}, {1.0, 2.0}};
  const std::uint32_t pcurve = ed.add_pcurve(source);
  const std::array<CoedgeSpec, 1> specs{CoedgeSpec{edge, false, pcurve, 0.0}};
  const LoopId loop = ed.add_loop(face, specs, true).value();

  const Result<VertexId> split = ed.split_edge(edge, 0.5);
  REQUIRE(split.has_value());
  {
    const Brep& staged = ed.brep();
    REQUIRE(staged.geometry().pcurves.size() == 3u);
    const Coedge& head = staged.coedge(staged.loop(loop).start);
    const Coedge& tail = staged.coedge(head.next);
    const Pcurve2d& low = staged.geometry().pcurves[head.pcurve];
    const Pcurve2d& high = staged.geometry().pcurves[tail.pcurve];
    REQUIRE(low.points.size() == 3u);
    REQUIRE(high.points.size() == 2u);
    CHECK(same_point(low.points.back(), high.points.front()));
    CHECK(low.points.back()[0] == 1.0);
    CHECK(low.points.back()[1] == doctest::Approx(0.5));
    CHECK(pcurve_length(low) == doctest::Approx(1.5));
    CHECK(pcurve_length(high) == doctest::Approx(1.5));
  }

  REQUIRE(ed.join_vertex(split.value()).has_value());
  const auto committed = ed.commit();
  REQUIRE(committed.has_value());
  const Brep& b = committed.value();
  const Coedge& survivor = b.coedge(b.loop(loop).start);
  REQUIRE(survivor.pcurve != kNoGeometry);
  const Pcurve2d& rebuilt = b.geometry().pcurves[survivor.pcurve];
  // The original three points plus the junction, with the junction stored once.
  REQUIRE(rebuilt.points.size() == 4u);
  CHECK(same_point(rebuilt.points[0], source.points[0]));
  CHECK(same_point(rebuilt.points[1], source.points[1]));
  CHECK(same_point(rebuilt.points[2], Point2d{1.0, 0.5}));
  CHECK(same_point(rebuilt.points[3], source.points[2]));
  CHECK(validate_topology(b).ok());
}

TEST_CASE("split edge aborts unchanged when a pcurve cannot be split") {
  SolidEditor ed;
  const VertexId v0 = ed.add_vertex({0.0, 0.0, 0.0}, 0.0).value();
  const VertexId v1 = ed.add_vertex({1.0, 0.0, 0.0}, 0.0).value();
  const std::uint32_t line = ed.add_curve(LineCurve3{{0.0, 0.0, 0.0}, {1.0, 0.0, 0.0}});
  const EdgeId edge = ed.add_edge(v0, v1, line, 0.0, 1.0, 0.0).value();
  const std::uint32_t plane = ed.add_surface(Plane3{{0.0, 0.0, 1.0}, 0.0});
  const FaceId face = ed.add_face(plane, false, 0.0).value();
  const std::uint32_t pcurve = ed.add_pcurve(
      Pcurve2d{{0.0, 0.0}, {std::numeric_limits<double>::quiet_NaN(), 1.0}});
  const std::array<CoedgeSpec, 1> specs{CoedgeSpec{edge, false, pcurve, 0.0}};
  const LoopId loop = ed.add_loop(face, specs, true).value();
  const std::uint32_t pcurves_before =
      static_cast<std::uint32_t>(ed.brep().geometry().pcurves.size());

  const Result<VertexId> split = ed.split_edge(edge, 0.5);
  REQUIRE_FALSE(split.has_value());
  CHECK(split.error().code == Errc::non_finite);

  const Brep& staged = ed.brep();
  CHECK(staged.edge_alive(edge));
  CHECK(staged.edge_count() == 1u);
  CHECK(staged.vertex_count() == 2u);
  CHECK(staged.coedge_count() == 1u);
  CHECK(staged.geometry().pcurves.size() == pcurves_before);
  CHECK(staged.loop(loop).start != CoedgeId{});
}
