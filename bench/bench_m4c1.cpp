#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <variant>
#include <vector>

#include <geom/brep/assemble.hpp>
#include <geom/brep/chamfer.hpp>
#include <geom/brep/mass_properties.hpp>
#include <geom/brep/primitives.hpp>
#include <geom/brep/shell.hpp>
#include <geom/core/error.hpp>
#include <geom/surfaces/plane.hpp>

namespace {

using Clock = std::chrono::steady_clock;

double seconds_since(Clock::time_point start) {
  return std::chrono::duration<double>(Clock::now() - start).count();
}

struct Stats {
  double min_ms = 0.0;
  double median_ms = 0.0;
  double max_ms = 0.0;
};

Stats summarize(std::array<double, 5> samples_ms) {
  std::sort(samples_ms.begin(), samples_ms.end());
  return Stats{samples_ms.front(), samples_ms[2], samples_ms.back()};
}

void print_stats(const char* label, const Stats& stats, std::size_t runs) {
  std::printf("%s: min %.3f ms, median %.3f ms, max %.3f ms (%zu runs)\n", label,
              stats.min_ms, stats.median_ms, stats.max_ms, runs);
}

bool close_to(double actual, double expected) {
  return std::abs(actual - expected) <= 1e-6 * std::abs(expected);
}

geom::FaceId top_face(const geom::Brep& brep) {
  for (const geom::FaceId face : brep.face_ids()) {
    const auto& plane =
        std::get<geom::Plane3>(brep.geometry().surfaces[brep.face(face).surface]);
    const geom::Vec3d normal =
        brep.face(face).reversed ? plane.normal * -1.0 : plane.normal;
    if (normal[2] > 0.5) return face;
  }
  return geom::FaceId{};
}

// The four faces of the right tetrahedron with legs `a` along +X, +Y and +Z,
// copied from the acceptance test: outward normals and CCW rings about them.
geom::Brep right_tetrahedron(double a) {
  const auto p = [](double x, double y, double z) { return geom::Point3d{x, y, z}; };
  const std::array<geom::AssembledFace, 4> faces{{
      {{p(0, 0, 0), p(0, 0, a), p(0, a, 0)}, {}, {-1, 0, 0}},
      {{p(0, 0, 0), p(a, 0, 0), p(0, 0, a)}, {}, {0, -1, 0}},
      {{p(0, a, 0), p(a, 0, 0), p(0, 0, 0)}, {}, {0, 0, -1}},
      {{p(a, 0, 0), p(0, a, 0), p(0, 0, a)}, {},
       {1 / std::sqrt(3.0), 1 / std::sqrt(3.0), 1 / std::sqrt(3.0)}},
  }};
  const geom::Result<geom::Brep> tetra = geom::assemble_solid(faces, 1e-9);
  if (!tetra.has_value()) {
    std::fprintf(stderr, "assemble_solid(tetra) failed: %s\n",
                 geom::to_string(tetra.error().code));
    return {};
  }
  return tetra.value();
}

}  // namespace

int main() {
  const geom::Result<geom::Brep> box =
      geom::make_box(geom::Point3d{0.0, 0.0, 0.0}, geom::Point3d{100.0, 100.0, 100.0}, 1e-9);
  if (!box.has_value()) {
    std::fprintf(stderr, "make_box failed: %s\n", geom::to_string(box.error().code));
    return 2;
  }

  const std::vector<geom::EdgeId> edges = box.value().edge_ids();
  const geom::ChamferOptions chamfer_options{1.0, 1e-9};

  std::array<double, 5> chamfer_samples{};
  double chamfer_volume = 0.0;
  for (std::size_t i = 0; i < chamfer_samples.size(); ++i) {
    const Clock::time_point start = Clock::now();
    const geom::Result<geom::ChamferResult> result =
        geom::chamfer_edges(box.value(), edges, chamfer_options);
    chamfer_samples[i] = seconds_since(start) * 1000.0;
    if (!result.has_value()) {
      std::fprintf(stderr, "chamfer_edges failed: %s: %s\n",
                   geom::to_string(result.error().code), result.error().message);
      return 2;
    }
    const geom::Result<geom::MassProperties> props = geom::mass_properties(result.value().brep);
    if (!props.has_value()) {
      std::fprintf(stderr, "mass_properties(chamfer) failed: %s\n",
                   geom::to_string(props.error().code));
      return 2;
    }
    chamfer_volume = props.value().volume;
    if (i == 0 && result.value().chamfered_edges != edges.size()) {
      std::fprintf(stderr, "chamfered %u of %zu edges\n", result.value().chamfered_edges,
                   edges.size());
      return 2;
    }
  }

  const geom::ShellOptions shell_options{1.0, 1e-9};
  const std::vector<geom::FaceId> open_faces;
  std::array<double, 5> shell_samples{};
  double shell_volume = 0.0;
  for (std::size_t i = 0; i < shell_samples.size(); ++i) {
    const Clock::time_point start = Clock::now();
    const geom::Result<geom::ShellResult> result =
        geom::shell_solid(box.value(), open_faces, shell_options);
    shell_samples[i] = seconds_since(start) * 1000.0;
    if (!result.has_value()) {
      std::fprintf(stderr, "shell_solid failed: %s: %s\n",
                   geom::to_string(result.error().code), result.error().message);
      return 2;
    }
    const geom::Result<geom::MassProperties> props = geom::mass_properties(result.value().brep);
    if (!props.has_value()) {
      std::fprintf(stderr, "mass_properties(shell) failed: %s\n",
                   geom::to_string(props.error().code));
      return 2;
    }
    shell_volume = props.value().volume;
    if (result.value().cavity_shells != 1) {
      std::fprintf(stderr, "expected 1 cavity shell, got %u\n", result.value().cavity_shells);
      return 2;
    }
  }

  // Open top: the cavity protrudes through +Z and the outer box clips it.
  const std::array<geom::FaceId, 1> open{top_face(box.value())};
  geom::ShellOptions open_options;
  open_options.thickness = 1.0;
  open_options.tolerance = 1e-9;
  std::array<double, 5> open_samples{};
  double open_volume = 0.0;
  for (std::size_t i = 0; i < open_samples.size(); ++i) {
    const Clock::time_point start = Clock::now();
    const geom::Result<geom::ShellResult> opened =
        geom::shell_solid(box.value(), open, open_options);
    open_samples[i] = seconds_since(start) * 1000.0;
    if (!opened.has_value()) {
      std::fprintf(stderr, "shell_solid(open) failed: %s: %s\n",
                   geom::to_string(opened.error().code), opened.error().message);
      return 3;
    }
    const geom::Result<geom::MassProperties> props =
        geom::mass_properties(opened.value().brep);
    if (!props.has_value()) {
      std::fprintf(stderr, "mass_properties(open shell) failed: %s\n",
                   geom::to_string(props.error().code));
      return 3;
    }
    open_volume = props.value().volume;
    if (opened.value().cavity_shells != 0) {
      std::fprintf(stderr, "open shell reported %u cavity shells\n",
                   opened.value().cavity_shells);
      return 3;
    }
  }
  // Removed: 98 * 98 * 99 = 950,796; shell = 1,000,000 - 950,796 = 49,204.
  if (!close_to(open_volume, 49204.0)) return 4;

  // Oblique tetrahedron, t = 0.1: analytic shell volume.
  const geom::Brep tetra = right_tetrahedron(2.0);
  if (tetra.face_ids().empty()) return 5;
  geom::ShellOptions tetra_options;
  tetra_options.thickness = 0.1;
  tetra_options.tolerance = 1e-9;
  const double expected =
      4.0 / 3.0 -
      std::pow(2.0 - 0.1 * (3.0 + std::sqrt(3.0)), 3) / 6.0;
  std::array<double, 5> tetra_samples{};
  double tetra_volume = 0.0;
  for (std::size_t i = 0; i < tetra_samples.size(); ++i) {
    const Clock::time_point start = Clock::now();
    const geom::Result<geom::ShellResult> tet =
        geom::shell_solid(tetra, {}, tetra_options);
    tetra_samples[i] = seconds_since(start) * 1000.0;
    if (!tet.has_value()) {
      std::fprintf(stderr, "shell_solid(tetra) failed: %s: %s\n",
                   geom::to_string(tet.error().code), tet.error().message);
      return 5;
    }
    const geom::Result<geom::MassProperties> props =
        geom::mass_properties(tet.value().brep);
    if (!props.has_value()) {
      std::fprintf(stderr, "mass_properties(tetra shell) failed: %s\n",
                   geom::to_string(props.error().code));
      return 5;
    }
    tetra_volume = props.value().volume;
  }
  if (!close_to(tetra_volume, expected)) return 6;

  std::printf("chamfer volume: %.12f (expected 999406)\n", chamfer_volume);
  std::printf("shell volume:   %.12f (expected 58808)\n", shell_volume);
  std::printf("open top volume: %.12f (expected 49204)\n", open_volume);
  std::printf("tetra volume:    %.12f (expected %.12f)\n", tetra_volume, expected);
  print_stats("chamfer all 12 edges", summarize(chamfer_samples), chamfer_samples.size());
  print_stats("closed shell", summarize(shell_samples), shell_samples.size());
  print_stats("open top shell", summarize(open_samples), open_samples.size());
  print_stats("oblique tetra shell", summarize(tetra_samples), tetra_samples.size());

  if (!close_to(chamfer_volume, 999406.0) || !close_to(shell_volume, 58808.0)) {
    std::fprintf(stderr, "volume checksum mismatch\n");
    return 3;
  }
  return 0;
}
