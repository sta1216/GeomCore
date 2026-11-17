// QA-18a: the public synthetic-domain dispatch splits a whole sphere on closed
// latitude cuts and is deterministic across runs.
#include <cstddef>
#include <vector>

#include <doctest/doctest.h>

#include <geom/brep/primitives.hpp>
#include <geom/brep/trim_face.hpp>
#include <geom/core/error.hpp>
#include <geom/primitives/point.hpp>

namespace {

constexpr double kPi = 3.14159265358979323846;

geom::FaceCut latitude_cut(double v) {
  constexpr int kPoints = 64;
  geom::FaceCut cut;
  cut.uv.reserve(kPoints + 1);
  for (int i = 0; i < kPoints; ++i) {
    const double u = 2.0 * kPi * static_cast<double>(i) / kPoints;
    cut.uv.push_back(geom::Point2d{u, v});
  }
  cut.uv.push_back(cut.uv.front());
  return cut;
}

}  // namespace

TEST_CASE("QA-18a: a whole sphere splits on two latitude cuts") {
  constexpr double tolerance = 1e-9;
  const geom::Brep sphere =
      geom::make_sphere(geom::Point3d{0.0, 0.0, 0.0}, 1.0, tolerance).value();
  const geom::FaceId face = sphere.face_ids().front();
  const std::vector<geom::FaceCut> cuts{latitude_cut(0.6), latitude_cut(1.0)};

  const geom::Result<std::vector<geom::CurvedFragment>> first =
      geom::split_face_curved(sphere, face, cuts, tolerance);
  REQUIRE(first.has_value());
  const std::vector<geom::CurvedFragment>& fragments = first.value();
  REQUIRE(fragments.size() == 3u);
  // A pole-free band is emitted as one ring joined through the seam columns,
  // so every fragment of this chart is single-looped.
  for (const geom::CurvedFragment& fragment : fragments) {
    CHECK(fragment.loops_uv.size() == 1u);
    CHECK(fragment.loops.size() == 1u);
    CHECK_FALSE(fragment.loops_uv.front().empty());
  }

  const geom::Result<std::vector<geom::CurvedFragment>> second =
      geom::split_face_curved(sphere, face, cuts, tolerance);
  REQUIRE(second.has_value());
  REQUIRE(second.value().size() == fragments.size());
  for (std::size_t i = 0; i < fragments.size(); ++i) {
    CHECK(second.value()[i].loops_uv.size() == fragments[i].loops_uv.size());
  }
}
