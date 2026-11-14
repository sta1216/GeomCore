// QA-14: the M3c-2 export corpus for torus operands of `boolean_solid_curved`.
//
// Three rows over one torus (R = 1, r = 0.5, axis +z):
//   - intersection with the half space z >= 0: an axis-perpendicular cut through
//     the seam-only torus band at v = pi,
//   - intersection with the half space x >= 0: a plane that contains the torus
//     axis, so the cut crosses the v seam,
//   - intersection with a box that contains the whole torus: no cut at all.
//
// Every row requires a success, a `validate_topology`-clean single solid whose
// shells are closed, the analytic volume (1e-6 relative) and a second run that
// reproduces the volume bit for bit. Cut rows additionally require at least one
// cut, at least one kept fragment and `max_cut_residual <= 1e-3`; the no-cut row
// requires zero cuts and skips the residual bound. Cone and sphere rows are
// intentionally absent: both are pinned as `unsupported` diagnostics by Tasks 5
// and 7, not measured here.
#include <doctest/doctest.h>

#include <array>
#include <cstddef>
#include <string_view>
#include <vector>

#include <geom/brep/boolean_solid_curved.hpp>
#include <geom/brep/detail/curved_volume.hpp>
#include <geom/brep/primitives.hpp>
#include <geom/brep/validate.hpp>
#include <geom/core/error.hpp>

namespace {

constexpr double kTolerance = 1e-9;
// (2 * pi^2 * R * r^2) / 2 with R = 1, r = 0.5: half of the torus volume.
constexpr double kHalfVolume = 2.4674011002723395;
// 2 * pi^2 * R * r^2 with R = 1, r = 0.5.
constexpr double kWholeVolume = 4.934802200544679;
constexpr double kResidualBound = 1e-3;
constexpr double kVolumeRelative = 1e-6;

struct HalvedBox {
  const char* name;
  geom::Point3d lower;
  geom::Point3d upper;
  double expected;
  bool cut;
};

std::string_view message_of(const geom::Error& error) { return std::string_view(error.message); }

geom::Brep make_box_or_fail(geom::Point3d lower, geom::Point3d upper) {
  const geom::Result<geom::Brep> box = geom::make_box(lower, upper, kTolerance);
  REQUIRE_MESSAGE(box.has_value(),
                  (box.has_value() ? std::string_view{} : message_of(box.error())));
  return box.value();
}

geom::SolidId only_solid(const geom::Brep& brep) {
  const std::vector<geom::SolidId> solids = brep.solid_ids();
  REQUIRE(solids.size() == 1);
  return solids.front();
}

// One measured run: the validated volume and the driver report. The Brep itself
// is not returned because the second run only needs its volume.
struct Measured {
  double volume = 0.0;
  geom::CurvedBooleanReport report;
};

Measured measure(const geom::Brep& torus, const geom::Brep& box, geom::SolidBooleanOp op) {
  const geom::Result<geom::CurvedBooleanResult> result = geom::boolean_solid_curved(
      torus, only_solid(torus), box, only_solid(box), op, geom::CurvedBooleanOptions{kTolerance});
  REQUIRE_MESSAGE(result.has_value(),
                  (result.has_value() ? std::string_view{} : message_of(result.error())));
  const geom::Brep& brep = result.value().brep;
  CHECK(geom::validate_topology(brep).ok());
  const std::vector<geom::SolidId> solids = brep.solid_ids();
  REQUIRE(solids.size() == 1);
  const std::span<const geom::ShellId> shells = brep.solid_shells(solids.front());
  REQUIRE(shells.size() == 1);
  for (const geom::ShellId shell : shells) {
    CHECK(brep.shell(shell).closed);
  }
  const geom::Result<double> volume =
      geom::detail::curved_solid_volume(brep, solids.front(), kTolerance);
  REQUIRE_MESSAGE(volume.has_value(),
                  (volume.has_value() ? std::string_view{} : message_of(volume.error())));
  Measured measured;
  measured.volume = volume.value();
  measured.report = result.value().report;
  return measured;
}

}  // namespace

TEST_CASE("QA-14: curved booleans of analytic primitives") {
  const geom::Result<geom::Brep> torus_result =
      geom::make_torus(geom::Point3d{0.0, 0.0, 0.0}, geom::Vec3d{0.0, 0.0, 1.0}, 1.0, 0.5,
                       kTolerance);
  REQUIRE_MESSAGE(torus_result.has_value(),
                  (torus_result.has_value() ? std::string_view{}
                                            : message_of(torus_result.error())));
  const geom::Brep& torus = torus_result.value();

  const std::array<HalvedBox, 3> rows{{
      {"torus intersection z >= 0", {-2.0, -2.0, 0.0}, {2.0, 2.0, 2.0}, kHalfVolume, true},
      {"torus intersection x >= 0", {0.0, -2.0, -2.0}, {2.0, 2.0, 2.0}, kHalfVolume, true},
      {"torus inside a containing box", {-2.0, -2.0, -2.0}, {2.0, 2.0, 2.0}, kWholeVolume, false},
  }};

  for (const HalvedBox& row : rows) {
    CAPTURE(std::string(row.name));
    const geom::Brep box = make_box_or_fail(row.lower, row.upper);

    const Measured first =
        measure(torus, box, geom::SolidBooleanOp::intersection);
    CHECK(first.volume == doctest::Approx(row.expected).epsilon(kVolumeRelative));
    if (row.cut) {
      CHECK(first.report.cuts >= 1);
      CHECK(first.report.kept >= 1);
      CHECK(first.report.max_cut_residual <= kResidualBound);
    } else {
      CHECK(first.report.cuts == 0);
    }

    const Measured second =
        measure(torus, box, geom::SolidBooleanOp::intersection);
    CHECK(second.volume == first.volume);
    CHECK(second.report.cuts == first.report.cuts);
    CHECK(second.report.kept == first.report.kept);
    CHECK(second.report.max_cut_residual == first.report.max_cut_residual);
  }
}
