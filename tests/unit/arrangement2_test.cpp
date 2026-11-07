#include <doctest/doctest.h>

#include <geom/cgeom/arrangement2.hpp>

#include <algorithm>
#include <cmath>

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
