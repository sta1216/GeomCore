#include <doctest/doctest.h>

#include <geom/cgeom/arrangement2.hpp>
#include <geom/cgeom/polygon2.hpp>

#include <algorithm>
#include <cmath>
#include <limits>
#include <utility>
#include <vector>

using geom::v1::ArrangementSegment;
using geom::v1::Point2d;

namespace {

using geom::v1::detail::SplitResult;
using geom::v1::detail::split_segments;

int count_with_source(const SplitResult& r, std::uint32_t source) {
  int n = 0;
  for (const auto& s : r.segments) if (s.source == source) ++n;
  return n;
}

// Same normalized perpendicular-distance definition as the header, but written
// independently so the test does not simply re-run the implementation.
long double perpendicular_distance(const Point2d& p, const Point2d& a, const Point2d& b) {
  long double dx = static_cast<long double>(b[0]) - static_cast<long double>(a[0]);
  long double dy = static_cast<long double>(b[1]) - static_cast<long double>(a[1]);
  const long double dm = std::max(std::fabs(dx), std::fabs(dy));
  if (!(dm > 0.0L)) return std::numeric_limits<long double>::infinity();
  dx /= dm;
  dy /= dm;
  long double px = static_cast<long double>(p[0]) - static_cast<long double>(a[0]);
  long double py = static_cast<long double>(p[1]) - static_cast<long double>(a[1]);
  const long double pm = std::max(std::fabs(px), std::fabs(py));
  if (pm == 0.0L) return 0.0L;
  px /= pm;
  py /= pm;
  return std::fabs(dx * py - dy * px) / std::hypot(dx, dy) * pm;
}

void check_pairs_sorted(const SplitResult& r) {
  for (std::size_t i = 0; i < r.input_pairs.size(); ++i) {
    CHECK(r.input_pairs[i].first < r.input_pairs[i].second);
    if (i > 0) CHECK(r.input_pairs[i - 1] < r.input_pairs[i]);
  }
}

// Parameter of p on ab, using the dominant axis so huge and tiny mixed spans
// do not overflow or underflow the test helper itself.
long double parameter_on(const Point2d& p, const Point2d& a, const Point2d& b) {
  const long double dx = static_cast<long double>(b[0]) - static_cast<long double>(a[0]);
  const long double dy = static_cast<long double>(b[1]) - static_cast<long double>(a[1]);
  if (std::fabs(dx) >= std::fabs(dy)) {
    return (static_cast<long double>(p[0]) - static_cast<long double>(a[0])) / dx;
  }
  return (static_cast<long double>(p[1]) - static_cast<long double>(a[1])) / dy;
}

// Sine of the angle between p - a and b - a: zero when p is exactly on the
// supporting line, independent of the coordinate scale.
long double collinear_sine(const Point2d& p, const Point2d& a, const Point2d& b) {
  const long double dx = static_cast<long double>(b[0]) - static_cast<long double>(a[0]);
  const long double dy = static_cast<long double>(b[1]) - static_cast<long double>(a[1]);
  const long double dm = std::max(std::fabs(dx), std::fabs(dy));
  if (!(dm > 0.0L)) return std::numeric_limits<long double>::infinity();
  const long double ux = dx / dm;
  const long double uy = dy / dm;
  const long double px = static_cast<long double>(p[0]) - static_cast<long double>(a[0]);
  const long double py = static_cast<long double>(p[1]) - static_cast<long double>(a[1]);
  const long double pm = std::max(std::fabs(px), std::fabs(py));
  if (pm == 0.0L) return 0.0L;
  return std::fabs(ux * (py / pm) - uy * (px / pm));
}

}  // namespace

TEST_CASE("split: crossing X splits both segments") {
  const ArrangementSegment segs[2] = {
      {Point2d{0, 0}, Point2d{2, 2}, 0},
      {Point2d{0, 2}, Point2d{2, 0}, 1},
  };
  const SplitResult r = split_segments(segs);
  CHECK(r.segments.size() == 4);
  CHECK(count_with_source(r, 0) == 2);
  CHECK(count_with_source(r, 1) == 2);
  CHECK(r.input_pairs.size() == 1);
  check_pairs_sorted(r);
  CHECK(r.max_residual <= 1e-12);
}

TEST_CASE("split: collinear overlap and T junction") {
  const ArrangementSegment overlap[2] = {
      {Point2d{0, 0}, Point2d{2, 0}, 0},
      {Point2d{1, 0}, Point2d{3, 0}, 1},
  };
  const SplitResult r = split_segments(overlap);
  CHECK(r.segments.size() == 4);  // 0-1,1-2 (src0); 1-2,2-3 (src1)
  CHECK(r.input_pairs.size() == 1);

  const ArrangementSegment t[2] = {
      {Point2d{0, 0}, Point2d{2, 0}, 0},
      {Point2d{1, 0}, Point2d{1, 1}, 1},
  };
  const SplitResult rt = split_segments(t);
  CHECK(rt.segments.size() == 3);
  CHECK(rt.input_pairs.size() == 1);
  CHECK(rt.input_pairs[0] == std::make_pair(0u, 1u));
  check_pairs_sorted(rt);
  // The contact is the exact endpoint (1, 0), so the T has no residual.
  CHECK(rt.max_residual == 0.0);
}

// b - a overflows to +inf for this horizontal segment (MSVC long double is
// double), so the old parameter evaluation divided inf by inf: the endpoint
// touch at (0, 0) became NaN, was clamped to t = 0 and deduped, and the
// vertical segment's endpoint was left sitting on the horizontal's interior
// with max_residual = +inf.
TEST_CASE("split: overflowing endpoint difference splits a T junction at the origin") {
  const ArrangementSegment segs[2] = {
      {Point2d{-1e308, 0.0}, Point2d{1e308, 0.0}, 0},
      {Point2d{0.0, 0.0}, Point2d{0.0, 1.0}, 1},
  };
  const SplitResult r = split_segments(segs);
  REQUIRE(r.segments.size() == 3);
  REQUIRE(r.input_pairs.size() == 1);
  CHECK(r.input_pairs[0] == std::make_pair(0u, 1u));
  CHECK(count_with_source(r, 0) == 2);
  CHECK(count_with_source(r, 1) == 1);
  CHECK(r.max_residual == 0.0L);
  CHECK(static_cast<double>(geom::v1::detail::split_perp_distance(
            Point2d{0.0, 1.0}, Point2d{-1e308, 0.0}, Point2d{1e308, 0.0})) ==
        doctest::Approx(1.0));

  CHECK(r.segments[0].a[0] == -1e308);
  CHECK(r.segments[0].a[1] == 0.0);
  CHECK(r.segments[0].b[0] == 0.0);
  CHECK(r.segments[0].b[1] == 0.0);
  CHECK(r.segments[1].a[0] == 0.0);
  CHECK(r.segments[1].a[1] == 0.0);
  CHECK(r.segments[1].b[0] == 1e308);
  CHECK(r.segments[1].b[1] == 0.0);
}

TEST_CASE("split: triple point and duplicate segments") {
  const ArrangementSegment triple[3] = {
      {Point2d{-1, -1}, Point2d{1, 1}, 0},
      {Point2d{-1, 1}, Point2d{1, -1}, 0},
      {Point2d{-2, 0}, Point2d{2, 0}, 0},
  };
  const SplitResult r = split_segments(triple);
  CHECK(r.segments.size() == 6);
  CHECK(r.input_pairs.size() == 3);
  check_pairs_sorted(r);

  const ArrangementSegment dup[2] = {
      {Point2d{0, 0}, Point2d{1, 0}, 0},
      {Point2d{0, 0}, Point2d{1, 0}, 1},
  };
  const SplitResult rd = split_segments(dup);
  CHECK(rd.segments.size() == 2);
  CHECK(rd.input_pairs.size() == 1);
}

TEST_CASE("split: disjoint and vertical segments") {
  const ArrangementSegment disjoint[3] = {
      {Point2d{0, 0}, Point2d{1, 0}, 0},
      {Point2d{0, 1}, Point2d{1, 1}, 1},
      {Point2d{2, 0}, Point2d{2, 2}, 2},
  };
  const SplitResult r = split_segments(disjoint);
  CHECK(r.segments.size() == 3);
  CHECK(r.input_pairs.empty());
  CHECK(r.max_residual == 0.0);

  const ArrangementSegment cross[2] = {
      {Point2d{1, 0}, Point2d{1, 2}, 0},
      {Point2d{0, 1}, Point2d{2, 1}, 1},
  };
  CHECK(split_segments(cross).segments.size() == 4);
}

TEST_CASE("split: near-parallel crossing uses exact determinants") {
  // Both direction products are about 2^60 and cancel to a determinant of 2.
  // A floating denominator test sees 0 here; the exact orientation signs are
  // strict (-, +, +, -), so the pair must still be split.
  const double k = 1073741824.0;  // 2^30
  const ArrangementSegment segs[2] = {
      {Point2d{0, 0}, Point2d{k, k - 1}, 0},
      {Point2d{-1, -1}, Point2d{k + 1, k}, 1},
  };
  const SplitResult r = split_segments(segs);
  REQUIRE(r.segments.size() == 4);
  REQUIRE(r.input_pairs.size() == 1);
  CHECK(r.input_pairs[0] == std::make_pair(0u, 1u));
  check_pairs_sorted(r);

  const Point2d expected{k / 2.0, k / 2.0 - 0.5};
  const Point2d splits[4] = {
      r.segments[0].b, r.segments[1].a, r.segments[2].b, r.segments[3].a};
  for (const Point2d& p : splits) {
    CHECK(std::fabs(p[0] - expected[0]) <= 1e-9);
    CHECK(std::fabs(p[1] - expected[1]) <= 1e-9);
    CHECK(perpendicular_distance(p, segs[0].a, segs[0].b) <= 1e-9);
    CHECK(perpendicular_distance(p, segs[1].a, segs[1].b) <= 1e-9);
  }
  CHECK(r.max_residual <= 1e-9);
}

TEST_CASE("split: residual captures an off-line contact") {
  // A long segment crossed by a vertical one at a coordinate where neither
  // parameterization lands exactly on the other supporting line: the exact
  // intersection has a y offset that double cannot represent at this scale.
  const ArrangementSegment segs[2] = {
      {Point2d{0, 0}, Point2d{2.0e15, 2.0e15 + 0.25}, 0},
      {Point2d{3.0e14, 0}, Point2d{3.0e14, 2.0e15}, 1},
  };
  const SplitResult r = split_segments(segs);
  REQUIRE(r.segments.size() == 4);
  REQUIRE(r.input_pairs.size() == 1);
  CHECK(r.max_residual > 0.0);
  CHECK(r.max_residual <= 0.5);

  const Point2d splits[4] = {
      r.segments[0].b, r.segments[1].a, r.segments[2].b, r.segments[3].a};
  for (const Point2d& p : splits) {
    CHECK(perpendicular_distance(p, segs[0].a, segs[0].b) <= r.max_residual + 1e-12);
    CHECK(perpendicular_distance(p, segs[1].a, segs[1].b) <= r.max_residual + 1e-12);
  }
}

TEST_CASE("split: mixed-span 2e200 x 1e-160 uses the unscaled frame") {
  // A global 2^-e frame collapses the 1e-160 offset to zero and used to fall
  // back to t = u = 0.5, bending seg1 to (5e199, 0) and giving t = 0.5.
  const ArrangementSegment segs[2] = {
      {Point2d{-1e200, 0}, Point2d{2e200, 0}, 0},
      {Point2d{1e200, -1e-160}, Point2d{1e200, 1e-160}, 1},
  };
  const SplitResult r = split_segments(segs);
  REQUIRE(r.segments.size() == 4);
  REQUIRE(r.input_pairs.size() == 1);
  CHECK(r.input_pairs[0] == std::make_pair(0u, 1u));
  check_pairs_sorted(r);
  CHECK(r.max_residual == 0.0);

  const Point2d contact = r.segments[0].b;
  CHECK(std::fabs(parameter_on(contact, segs[0].a, segs[0].b) - 2.0L / 3.0L) <= 1e-9L);
  CHECK(std::fabs(parameter_on(contact, segs[1].a, segs[1].b) - 0.5L) <= 1e-9L);
  for (const auto& piece : r.segments) {
    const ArrangementSegment& in = piece.source == 0 ? segs[0] : segs[1];
    CHECK(collinear_sine(piece.a, in.a, in.b) <= 1e-12L);
    CHECK(collinear_sine(piece.b, in.a, in.b) <= 1e-12L);
  }
}

TEST_CASE("split: mixed-span 1e-160 x 1e200 uses the unscaled frame") {
  // The dual span: the global frame collapses to a zero-length seg1 and used
  // to return contact (0, 0) with t = 0.5 and max_residual = 0.
  const ArrangementSegment segs[2] = {
      {Point2d{-1e-160, 0}, Point2d{1e-160, 0}, 0},
      {Point2d{0.4e-160, -1e200}, Point2d{0.4e-160, 1e200}, 1},
  };
  const SplitResult r = split_segments(segs);
  REQUIRE(r.segments.size() == 4);
  REQUIRE(r.input_pairs.size() == 1);
  CHECK(r.input_pairs[0] == std::make_pair(0u, 1u));
  check_pairs_sorted(r);
  CHECK(r.max_residual == 0.0);

  const Point2d contact = r.segments[0].b;
  CHECK(std::fabs(parameter_on(contact, segs[0].a, segs[0].b) - 0.7L) <= 1e-9L);
  CHECK(std::fabs(parameter_on(contact, segs[1].a, segs[1].b) - 0.5L) <= 1e-9L);
  for (const auto& piece : r.segments) {
    const ArrangementSegment& in = piece.source == 0 ? segs[0] : segs[1];
    CHECK(collinear_sine(piece.a, in.a, in.b) <= 1e-12L);
    CHECK(collinear_sine(piece.b, in.a, in.b) <= 1e-12L);
  }
}

TEST_CASE("split: near-limit 1e307 crossing uses the scaled fallback") {
  // The unscaled products are ~1e614 and overflow, while the 2^k frame keeps
  // the ratios representable. This exercises the fallback path, not the error
  // path, and must still emit four collinear pieces.
  const ArrangementSegment segs[2] = {
      {Point2d{-1e307, 0}, Point2d{1e307, 0}, 0},
      {Point2d{0, -1e307}, Point2d{0, 1e307}, 1},
  };
  const SplitResult r = split_segments(segs);
  REQUIRE(r.segments.size() == 4);
  REQUIRE(r.input_pairs.size() == 1);
  CHECK(r.input_pairs[0] == std::make_pair(0u, 1u));
  check_pairs_sorted(r);
  CHECK(r.max_residual == 0.0);

  const Point2d contact = r.segments[0].b;
  CHECK(contact[0] == 0.0);
  CHECK(contact[1] == 0.0);
  for (const auto& piece : r.segments) {
    const ArrangementSegment& in = piece.source == 0 ? segs[0] : segs[1];
    CHECK(collinear_sine(piece.a, in.a, in.b) <= 1e-12L);
    CHECK(collinear_sine(piece.b, in.a, in.b) <= 1e-12L);
  }
}

TEST_CASE("split: uniform 1e200 and 1e-200 scaling") {
  const double scales[2] = {1e200, 1e-200};
  for (const double s : scales) {
    const ArrangementSegment segs[2] = {
        {Point2d{0, 0}, Point2d{2 * s, 2 * s}, 0},
        {Point2d{0, 2 * s}, Point2d{2 * s, 0}, 1},
    };
    const SplitResult r = split_segments(segs);
    REQUIRE(r.segments.size() == 4);
    REQUIRE(r.input_pairs.size() == 1);
    CHECK(std::isfinite(r.max_residual));
    check_pairs_sorted(r);

    const Point2d expected{s, s};
    const Point2d splits[4] = {
        r.segments[0].b, r.segments[1].a, r.segments[2].b, r.segments[3].a};
    for (const Point2d& p : splits) {
      CHECK(std::fabs(p[0] - expected[0]) <= 1e-9 * s);
      CHECK(std::fabs(p[1] - expected[1]) <= 1e-9 * s);
    }
  }
}

TEST_CASE("split: underflowing exact ratio is reported as unrepresentable") {
  // The exact crossing parameter is 1e-360, which rounds to zero. The contact
  // (1e-160, 0) lies on both lines, but adding it at t = 0 makes the t-only
  // dedupe drop it, leaving cd's vertex on ab's interior without a split there.
  const ArrangementSegment segs[2] = {
      {Point2d{0, 0}, Point2d{1e200, 0}, 0},
      {Point2d{1e-160, -1e-160}, Point2d{1e-160, 1e-160}, 1},
  };
  const SplitResult r = split_segments(segs);
  REQUIRE(r.segments.size() == 2);
  REQUIRE(r.input_pairs.size() == 1);
  CHECK(r.input_pairs[0] == std::make_pair(0u, 1u));
  CHECK(std::isinf(r.max_residual));
  for (const auto& piece : r.segments) {
    const ArrangementSegment& in = piece.source == 0 ? segs[0] : segs[1];
    CHECK(piece.a[0] == in.a[0]);
    CHECK(piece.a[1] == in.a[1]);
    CHECK(piece.b[0] == in.b[0]);
    CHECK(piece.b[1] == in.b[1]);
  }
}

TEST_CASE("split: beyond-span crossing is deterministic and unsplit") {
  // Both the unscaled and scaled determinants underflow here; the pair must
  // still be recorded, but neither input may be split or bent.
  const ArrangementSegment segs[2] = {
      {Point2d{0, 0}, Point2d{1e308, 0}, 0},
      {Point2d{1e-300, -1e-300}, Point2d{1e-300, 1e-300}, 1},
  };
  const SplitResult r = split_segments(segs);
  REQUIRE(r.segments.size() == 2);
  REQUIRE(r.input_pairs.size() == 1);
  CHECK(r.input_pairs[0] == std::make_pair(0u, 1u));
  CHECK(std::isinf(r.max_residual));
  for (const auto& piece : r.segments) {
    const ArrangementSegment& in = piece.source == 0 ? segs[0] : segs[1];
    CHECK(piece.a[0] == in.a[0]);
    CHECK(piece.a[1] == in.a[1]);
    CHECK(piece.b[0] == in.b[0]);
    CHECK(piece.b[1] == in.b[1]);
  }
}

TEST_CASE("split: zero-length input is filtered without renumbering outputs") {
  const ArrangementSegment segs[3] = {
      {Point2d{5, 5}, Point2d{5, 5}, 9},  // dropped
      {Point2d{0, 0}, Point2d{2, 2}, 0},
      {Point2d{0, 2}, Point2d{2, 0}, 1},
  };
  const SplitResult r = split_segments(segs);
  REQUIRE(r.segments.size() == 4);
  REQUIRE(r.input_pairs.size() == 1);
  CHECK(r.input_pairs[0] == std::make_pair(1u, 2u));
  check_pairs_sorted(r);
  for (const auto& piece : r.segments) {
    CHECK(piece.source != 9);
  }
}

using geom::v1::Arrangement2;
using geom::v1::Polygon2d;

namespace {

Polygon2d square(double x0, double y0, double x1, double y1) {
  Polygon2d p;
  p.vertices = {{x0, y0}, {x1, y0}, {x1, y1}, {x0, y1}};
  p.ring_starts = {0, 4};
  return p;
}

std::vector<ArrangementSegment> ring_segments(const Polygon2d& p, std::uint32_t source) {
  std::vector<ArrangementSegment> out;
  const auto r = p.ring(0);
  for (std::uint32_t i = 0; i < r.size(); ++i) {
    out.push_back({r[i], r[(i + 1) % r.size()], source});
  }
  return out;
}

std::size_t bounded_faces(const Arrangement2& arr) {
  std::size_t count = 0;
  for (const auto& face : arr.faces()) {
    if (face.bounded) ++count;
  }
  return count;
}

}  // namespace

TEST_CASE("arrangement: single square") {
  const Polygon2d p = square(0, 0, 1, 1);
  const auto segs = ring_segments(p, 0);
  const Arrangement2 arr = Arrangement2::build(segs);
  CHECK(arr.vertex_count() == 4);
  CHECK(arr.edge_count() == 4);
  CHECK(arr.face_count() == 2);
  CHECK(arr.source_count() == 1);
  CHECK(arr.winding(arr.locate(Point2d{0.5, 0.5}), 0) == 1);
  CHECK(arr.locate(Point2d{2, 2}) == 0xFFFFFFFFu);
  CHECK(arr.max_residual() == 0.0);
}

TEST_CASE("arrangement: two overlapping squares") {
  std::vector<ArrangementSegment> segs = ring_segments(square(0, 0, 1, 1), 0);
  for (auto s : ring_segments(square(0.5, 0.5, 1.5, 1.5), 1)) segs.push_back(s);
  const Arrangement2 arr = Arrangement2::build(segs);
  CHECK(arr.vertex_count() == 10);
  CHECK(arr.edge_count() == 12);
  CHECK(arr.face_count() == 4);
  CHECK(arr.winding(arr.locate(Point2d{0.25, 0.75}), 0) == 1);
  const std::uint32_t overlap = arr.locate(Point2d{0.75, 0.75});
  CHECK(arr.winding(overlap, 0) == 1);
  CHECK(arr.winding(overlap, 1) == 1);
  CHECK(arr.faces()[overlap].winding == 2);
  const std::uint32_t right = arr.locate(Point2d{1.25, 0.75});
  CHECK(arr.winding(right, 1) == 1);
  CHECK(arr.winding(right, 0) == 0);
  CHECK(arr.locate(Point2d{2, 2}) == 0xFFFFFFFFu);
}

TEST_CASE("arrangement: shared edge keeps equal winding on both sides") {
  std::vector<ArrangementSegment> segs = ring_segments(square(0, 0, 1, 1), 0);
  for (auto s : ring_segments(square(1, 0, 2, 1), 1)) segs.push_back(s);
  const Arrangement2 arr = Arrangement2::build(segs);
  CHECK(arr.vertex_count() == 6);
  CHECK(arr.edge_count() == 7);
  CHECK(arr.face_count() == 3);
  CHECK(arr.winding(arr.locate(Point2d{0.5, 0.5}), 0) == 1);
  CHECK(arr.winding(arr.locate(Point2d{1.5, 0.5}), 1) == 1);
}

TEST_CASE("arrangement: vertical chain and point locate on boundary") {
  const Polygon2d p = square(-1, -1, 1, 1);
  const Arrangement2 arr = Arrangement2::build(ring_segments(p, 0));
  CHECK(arr.winding(arr.locate(Point2d{0, 0}), 0) == 1);
  CHECK(arr.locate(Point2d{1, 0}) == 0xFFFFFFFFu);  // boundary has no unique face
  // Bottom-edge midpoint: the even-odd ray test alone reports this point as
  // inside, so only the explicit on-edge rejection returns the sentinel.
  CHECK(arr.locate(Point2d{0, -1}) == 0xFFFFFFFFu);
}

TEST_CASE("arrangement: three nested squares accumulate per-source winding") {
  std::vector<ArrangementSegment> segs = ring_segments(square(0, 0, 10, 10), 0);
  for (auto s : ring_segments(square(2, 2, 8, 8), 1)) segs.push_back(s);
  for (auto s : ring_segments(square(3, 3, 7, 7), 2)) segs.push_back(s);
  const Arrangement2 arr = Arrangement2::build(segs);
  CHECK(arr.source_count() == 3);
  CHECK(bounded_faces(arr) == 3);
  const std::uint32_t core = arr.locate(Point2d{5, 5});
  CHECK(core != Arrangement2::kNoFace);
  CHECK(arr.faces()[core].bounded);
  CHECK(arr.winding(core, 0) == 1);
  CHECK(arr.winding(core, 1) == 1);
  CHECK(arr.winding(core, 2) == 1);
  CHECK(arr.faces()[core].winding == 3);
  const std::uint32_t ring0 = arr.locate(Point2d{1, 1});
  CHECK(arr.winding(ring0, 0) == 1);
  CHECK(arr.winding(ring0, 1) == 0);
  CHECK(arr.faces()[ring0].winding == 1);
  const std::uint32_t ring1 = arr.locate(Point2d{9, 9});  // inside 0..10 only
  CHECK(arr.winding(ring1, 0) == 1);
  CHECK(arr.winding(ring1, 1) == 0);
  CHECK(arr.faces()[ring1].winding == 1);
  CHECK(arr.locate(Point2d{11, 11}) == Arrangement2::kNoFace);
}

TEST_CASE("arrangement: extreme coordinate nested squares") {
  const double scales[2] = {1e150, 1e200};
  for (const double s : scales) {
    std::vector<ArrangementSegment> segs = ring_segments(square(s, s, 4 * s, 4 * s), 0);
    for (auto seg : ring_segments(square(2 * s, 2 * s, 3 * s, 3 * s), 1)) {
      segs.push_back(seg);
    }
    const Arrangement2 arr = Arrangement2::build(segs);
    CHECK(arr.vertex_count() == 8);
    CHECK(arr.edge_count() == 8);
    CHECK(bounded_faces(arr) == 2);
    const std::uint32_t inner = arr.locate(Point2d{2.5 * s, 2.5 * s});
    CHECK(inner != Arrangement2::kNoFace);
    CHECK(arr.faces()[inner].bounded);
    CHECK(arr.winding(inner, 0) == 1);
    CHECK(arr.winding(inner, 1) == 1);
    CHECK(arr.faces()[inner].winding == 2);
    const std::uint32_t ring = arr.locate(Point2d{1.5 * s, 1.5 * s});
    CHECK(ring != Arrangement2::kNoFace);
    CHECK(arr.winding(ring, 0) == 1);
    CHECK(arr.winding(ring, 1) == 0);
    CHECK(arr.locate(Point2d{5 * s, 5 * s}) == Arrangement2::kNoFace);

    std::vector<ArrangementSegment> negative =
        ring_segments(square(-4 * s, -4 * s, -s, -s), 0);
    for (auto seg : ring_segments(square(-3 * s, -3 * s, -2 * s, -2 * s), 1)) {
      negative.push_back(seg);
    }
    const Arrangement2 neg = Arrangement2::build(negative);
    CHECK(bounded_faces(neg) == 2);
    const std::uint32_t neg_inner = neg.locate(Point2d{-2.5 * s, -2.5 * s});
    CHECK(neg_inner != Arrangement2::kNoFace);
    CHECK(neg.winding(neg_inner, 0) == 1);
    CHECK(neg.winding(neg_inner, 1) == 1);
    CHECK(neg.faces()[neg_inner].winding == 2);
  }
}

TEST_CASE("arrangement: hole in one ring winds down") {
  Polygon2d hole;
  hole.vertices = {{4, 4}, {4, 6}, {6, 6}, {6, 4}};  // clockwise
  hole.ring_starts = {0, 4};
  std::vector<ArrangementSegment> segs = ring_segments(square(0, 0, 10, 10), 0);
  for (auto s : ring_segments(hole, 0)) segs.push_back(s);
  const Arrangement2 arr = Arrangement2::build(segs);
  CHECK(arr.face_count() == 4);
  CHECK(bounded_faces(arr) == 2);
  const std::uint32_t interior = arr.locate(Point2d{5, 5});
  CHECK(interior != Arrangement2::kNoFace);
  CHECK(arr.faces()[interior].bounded);
  CHECK(arr.faces()[interior].winding == 0);
  CHECK(arr.winding(interior, 0) == 0);
  const std::uint32_t body = arr.locate(Point2d{1, 1});
  CHECK(body != Arrangement2::kNoFace);
  CHECK(arr.faces()[body].winding == 1);
  CHECK(arr.winding(body, 0) == 1);
}

namespace {

Polygon2d polygon_of(std::vector<Point2d> vertices) {
  Polygon2d p;
  p.vertices = std::move(vertices);
  p.ring_starts = {0, static_cast<std::uint32_t>(p.vertices.size())};
  return p;
}

std::uint32_t unbounded_face(const Arrangement2& arr) {
  for (std::uint32_t f = 0; f < arr.faces().size(); ++f) {
    if (!arr.faces()[f].bounded) return f;
  }
  return Arrangement2::kNoFace;
}

}  // namespace

// The lexicographically smallest vertex of this kite has one outgoing edge
// pointing into the lower half-plane. The old angular seed picked that edge and
// started the component from the interior face, inverting every winding.
TEST_CASE("arrangement: kite seeds winding from its unbounded face") {
  const Polygon2d kite = polygon_of({{0, 0}, {2, 1}, {0, 4}, {-2, 1}});
  const Arrangement2 arr = Arrangement2::build(ring_segments(kite, 0));
  REQUIRE(arr.faces().size() == 2);
  const std::uint32_t inside = arr.locate(Point2d{0, 1});
  REQUIRE(inside != Arrangement2::kNoFace);
  CHECK(arr.faces()[inside].bounded);
  CHECK(arr.winding(inside, 0) == 1);
  CHECK(arr.faces()[inside].winding == 1);
  CHECK(arr.locate(Point2d{0, 5}) == Arrangement2::kNoFace);
  const std::uint32_t outside = unbounded_face(arr);
  REQUIRE(outside != Arrangement2::kNoFace);
  CHECK(arr.winding(outside, 0) == 0);
  CHECK(arr.faces()[outside].winding == 0);
}

// Same accidental dependency, one step more general: the leftmost vertex
// (-1, 1) has an outgoing edge with a negative y component.
TEST_CASE("arrangement: pentagon with a lower-half edge at its leftmost vertex") {
  const Polygon2d pentagon = polygon_of({{0, 0}, {3, -1}, {4, 1}, {2, 3}, {-1, 1}});
  const Arrangement2 arr = Arrangement2::build(ring_segments(pentagon, 0));
  const std::uint32_t inside = arr.locate(Point2d{1.6, 0.8});
  REQUIRE(inside != Arrangement2::kNoFace);
  CHECK(arr.faces()[inside].bounded);
  CHECK(arr.winding(inside, 0) == 1);
  CHECK(arr.faces()[inside].winding == 1);
  CHECK(arr.locate(Point2d{1.6, 5.0}) == Arrangement2::kNoFace);
}

// A spur traversed in both directions is a zero-area component. Zero area must
// compare smaller than every nonzero face; otherwise the spur is processed
// before the small square and its seed winding leaks into the square's face.
TEST_CASE("arrangement: zero-area spur does not outrank a small square") {
  std::vector<ArrangementSegment> segs = ring_segments(square(0.0, 0.0, 0.5, 0.5), 0);
  segs.push_back(ArrangementSegment{Point2d{0.1, 0.1}, Point2d{0.2, 0.1}, 1});
  segs.push_back(ArrangementSegment{Point2d{0.2, 0.1}, Point2d{0.1, 0.1}, 1});
  const Arrangement2 arr = Arrangement2::build(segs);
  const std::uint32_t face = arr.locate(Point2d{0.25, 0.25});
  REQUIRE(face != Arrangement2::kNoFace);
  CHECK(arr.faces()[face].bounded);
  CHECK(arr.winding(face, 0) == 1);
  CHECK(arr.winding(face, 1) == 0);
  CHECK(arr.faces()[face].winding == 1);
  // The degenerate spur orbit sits inside the square and must inherit its
  // source-0 winding; at the pre-fix ordering its seed resolved to 0 instead.
  int inside_faces = 0;
  for (std::uint32_t f = 0; f < arr.faces().size(); ++f) {
    if (arr.winding(f, 0) == 1) ++inside_faces;
  }
  CHECK(inside_faces == 2);
}

// One global vertex exponent collapses the unit square's area to zero next to
// the 1e200 square, so the small component loses its bounded face and locate()
// returns the sentinel. Per-cycle normalization keeps both components usable.
TEST_CASE("arrangement: disjoint components keep independent area scales") {
  std::vector<ArrangementSegment> segs = ring_segments(square(0.0, 0.0, 1.0, 1.0), 0);
  for (const auto& s : ring_segments(square(0.5e200, 0.5e200, 1e200, 1e200), 1)) {
    segs.push_back(s);
  }
  const Arrangement2 arr = Arrangement2::build(segs);
  const std::uint32_t unit = arr.locate(Point2d{0.5, 0.5});
  REQUIRE(unit != Arrangement2::kNoFace);
  CHECK(arr.faces()[unit].bounded);
  CHECK(arr.winding(unit, 0) == 1);
  const std::uint32_t huge = arr.locate(Point2d{0.75e200, 0.75e200});
  REQUIRE(huge != Arrangement2::kNoFace);
  CHECK(arr.faces()[huge].bounded);
  CHECK(arr.winding(huge, 1) == 1);
  CHECK(arr.locate(Point2d{-1.0, -1.0}) == Arrangement2::kNoFace);
  CHECK(arr.locate(Point2d{2e200, 2e200}) == Arrangement2::kNoFace);
}
// Zero-area faces are degenerate records (mantissa 0). Their magnitude must be
// below every nonzero face, even when the nonzero face's area is < 1 and its
// exponent is negative. The old helpers compared 0 * 2^0 against 0.5 * 2^-1 as
// if 1 < 0.5, which put a zero-area component ahead of a small square.
TEST_CASE("arrangement: zero-area faces sort below every nonzero face") {
  CHECK(geom::v1::detail::area_abs_less(0.0, 0, 0.5, -1));
  CHECK_FALSE(geom::v1::detail::area_abs_less(0.5, -1, 0.0, 0));
  CHECK_FALSE(geom::v1::detail::area_abs_less(0.0, 0, 0.0, 0));
  CHECK(geom::v1::detail::area_compare(0.0, 0, 0.5, -1) < 0);
  CHECK(geom::v1::detail::area_compare(0.5, -1, 0.0, 0) > 0);
  CHECK(geom::v1::detail::area_compare(0.0, 0, 0.0, 0) == 0);
  CHECK(geom::v1::detail::area_compare(0.0, 0, -0.5, -1) > 0);
}