// QA-8: planar Boolean on two faceted prisms.
//
// Builds two prisms with `assemble_solid` (caps plus one quad per side), checks
// each operand, then times `boolean_solid` for union / intersection / difference
// and reports the volume identity. A run passes an op only when its wall-clock
// time is at most 5 s, the result passes `validate_topology`, and the identity
// holds to 1e-9. Construction and validation are outside the timed region; at
// the official 5000 sides assembly alone takes about 15 s per operand.
//
// Exit codes: 0 all executed ops pass, 1 any op fails, 2 setup or Boolean
// error. Progress goes to stderr so the stdout QA8 lines stay a single
// machine-readable block; the three-op 5000-side run takes tens of minutes.
#include <geom/brep/assemble.hpp>
#include <geom/brep/boolean_solid.hpp>
#include <geom/brep/mass_properties.hpp>
#include <geom/brep/validate.hpp>

#include <charconv>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <numbers>
#include <string_view>
#include <utility>
#include <vector>

using namespace geom::v1;

namespace {

constexpr double kRadius = 10.0;
constexpr double kHeight = 4.0;
constexpr double kTol = 1e-9;
constexpr double kOffsetFactor = 0.25;
constexpr double kTargetMs = 5000.0;
constexpr double kTargetErr = 1e-9;
constexpr std::uint32_t kDefaultSides = 5000;
constexpr std::uint32_t kMinSides = 8;
constexpr std::uint32_t kMaxSides = 200000;

struct Op {
  const char* name;
  SolidBooleanOp code;
};

constexpr Op kOps[3] = {{"union", SolidBooleanOp::union_},
                        {"intersection", SolidBooleanOp::intersection},
                        {"difference", SolidBooleanOp::difference}};

struct OpRun {
  const char* name = "";
  double ms = 0.0;
  bool topology_ok = false;
  bool volume_ok = false;
  double volume = 0.0;
};

double ms_since(std::chrono::steady_clock::time_point start) {
  return std::chrono::duration<double, std::milli>(
             std::chrono::steady_clock::now() - start)
      .count();
}

void print_error(const Error& error) {
  std::printf("error=%s: %s\n", to_string(error.code), error.message);
  std::fflush(stdout);
}

// Outward-oriented faceted prism: reversed bottom cap, forward top cap and one
// quad per side, each quad winding counter-clockwise about its radial normal.
std::vector<AssembledFace> prism_faces(double center_x, std::uint32_t sides) {
  const double step = 2.0 * std::numbers::pi / static_cast<double>(sides);
  std::vector<Point3d> bottom(sides);
  std::vector<Point3d> top(sides);
  for (std::uint32_t i = 0; i < sides; ++i) {
    const double angle = step * static_cast<double>(i);
    const double x = center_x + kRadius * std::cos(angle);
    const double y = kRadius * std::sin(angle);
    bottom[i] = Point3d{x, y, 0.0};
    top[i] = Point3d{x, y, kHeight};
  }

  std::vector<AssembledFace> faces;
  faces.reserve(static_cast<std::size_t>(sides) + 2u);

  AssembledFace bottom_cap;
  bottom_cap.outer.assign(bottom.rbegin(), bottom.rend());
  bottom_cap.normal = Vec3d{0.0, 0.0, -1.0};
  faces.push_back(std::move(bottom_cap));

  AssembledFace top_cap;
  top_cap.outer = top;
  top_cap.normal = Vec3d{0.0, 0.0, 1.0};
  faces.push_back(std::move(top_cap));

  for (std::uint32_t i = 0; i < sides; ++i) {
    const std::uint32_t j = (i + 1u) % sides;
    const double mid = step * (static_cast<double>(i) + 0.5);
    AssembledFace side;
    side.outer = {bottom[i], bottom[j], top[j], top[i]};
    side.normal = Vec3d{std::cos(mid), std::sin(mid), 0.0};
    faces.push_back(std::move(side));
  }
  return faces;
}

// Checks one assembled prism and returns its measured volume. A failure is a
// setup error: the operands must be exactly the solid the benchmark describes,
// otherwise the QA8 time is not comparable with the recorded baseline.
bool verify_operand(const Brep& brep, std::uint32_t sides, const char* label,
                    double* volume_out) {
  const ValidationReport topology = validate_topology(brep);
  if (!topology.ok()) {
    const ValidationIssue& issue = topology.issues.front();
    std::printf("error=%s: operand %s: %s\n", to_string(issue.code), label,
                issue.message);
    return false;
  }
  const ValidationReport tolerances = validate_tolerances(brep);
  if (!tolerances.ok()) {
    const ValidationIssue& issue = tolerances.issues.front();
    std::printf("error=%s: operand %s: %s\n", to_string(issue.code), label,
                issue.message);
    return false;
  }
  const std::uint32_t expected_faces = sides + 2u;
  if (brep.face_count() != expected_faces) {
    std::printf("error=invalid_topology: operand %s: face_count %u != %u\n",
                label, static_cast<unsigned>(brep.face_count()),
                static_cast<unsigned>(expected_faces));
    return false;
  }

  const Result<MassProperties> mass = mass_properties(brep);
  if (!mass.has_value()) {
    std::printf("error=%s: operand %s: %s\n", to_string(mass.error().code),
                label, mass.error().message);
    return false;
  }
  const double volume = mass.value().volume;
  const double step = 2.0 * std::numbers::pi / static_cast<double>(sides);
  const double analytic = 0.5 * static_cast<double>(sides) * kRadius * kRadius *
                          std::sin(step) * kHeight;
  const double relative = std::fabs(volume - analytic) / analytic;
  if (!std::isfinite(volume) || !(relative <= 1e-12)) {
    std::printf("error=numeric_failure: operand %s: volume %.17g analytic "
                "%.17g relative %.3e\n",
                label, volume, analytic, relative);
    return false;
  }
  *volume_out = volume;
  return true;
}

// Runs one Boolean op and records the timing plus the result's topology and
// volume. The result's volume is not needed for a single-op run, where the
// identity is unavailable.
bool run_op(const Brep& a, const Brep& b, const Op& spec, bool measure_volume,
            OpRun* out) {
  std::fprintf(stderr, "[bench] %s start\n", spec.name);
  std::fflush(stderr);
  const auto start = std::chrono::steady_clock::now();
  const Result<SolidBooleanResult> result = boolean_solid(a, b, spec.code, {});
  const double ms = ms_since(start);
  std::fprintf(stderr, "[bench] %s done ms=%.3f\n", spec.name, ms);
  std::fflush(stderr);
  if (!result.has_value()) {
    print_error(result.error());
    return false;
  }

  const Brep& brep = result.value().brep;
  out->name = spec.name;
  out->ms = ms;
  out->topology_ok = validate_topology(brep).ok();
  if (measure_volume) {
    const Result<MassProperties> mass = mass_properties(brep);
    if (mass.has_value()) {
      out->volume = mass.value().volume;
      out->volume_ok = std::isfinite(out->volume);
    }
  }
  return true;
}

// One QA8 line. `have_error` is false when a result volume could not be
// measured: the identity is then printed as nan and the op cannot pass.
void print_qa8(const OpRun& run, bool have_error, double error, bool single_op,
               std::uint32_t faces_a, std::uint32_t faces_b) {
  const bool pass = run.ms <= kTargetMs && run.topology_ok &&
                    (!have_error || error <= kTargetErr) &&
                    (single_op || have_error);
  std::printf(
      "QA8 op=%s faces_a=%u faces_b=%u ms=%.3f volume_rel_err=", run.name,
      static_cast<unsigned>(faces_a), static_cast<unsigned>(faces_b), run.ms);
  if (have_error) {
    std::printf("%.3e", error);
  } else {
    std::printf("nan");
  }
  if (single_op) std::printf(" identity=unavailable");
  std::printf(" target=5000ms/1e-9 %s\n", pass ? "PASS" : "FAIL");
  std::fflush(stdout);
}

bool parse_sides(std::string_view text, std::uint32_t* out) {
  if (text.empty()) return false;
  unsigned long long value = 0;
  const auto [end, error] =
      std::from_chars(text.data(), text.data() + text.size(), value, 10);
  if (error != std::errc{} || end != text.data() + text.size()) return false;
  if (value < static_cast<unsigned long long>(kMinSides) ||
      value > static_cast<unsigned long long>(kMaxSides)) {
    return false;
  }
  *out = static_cast<std::uint32_t>(value);
  return true;
}

}  // namespace

int main(int argc, char** argv) {
  if (argc > 3) {
    std::printf("error=usage: geom_bench_m3 [sides] "
                "[union|intersection|difference|all]\n");
    return 2;
  }

  std::uint32_t sides = kDefaultSides;
  if (argc > 1) {
    if (!parse_sides(argv[1], &sides)) {
      std::printf("error=usage: sides must be an integer in [%u, %u]\n",
                  static_cast<unsigned>(kMinSides),
                  static_cast<unsigned>(kMaxSides));
      return 2;
    }
  }

  std::string_view op_name = argc > 2 ? std::string_view{argv[2]} : "all";
  bool single_op = false;
  const Op* single = nullptr;
  if (op_name == "all") {
    single_op = false;
  } else {
    for (const Op& candidate : kOps) {
      if (op_name == candidate.name) single = &candidate;
    }
    if (single == nullptr) {
      std::printf("error=usage: op must be union, intersection, difference or "
                  "all\n");
      return 2;
    }
    single_op = true;
  }

  const auto build_start = std::chrono::steady_clock::now();
  const Result<Brep> made_a = assemble_solid(prism_faces(0.0, sides), kTol);
  if (!made_a.has_value()) {
    print_error(made_a.error());
    return 2;
  }
  const Result<Brep> made_b =
      assemble_solid(prism_faces(kOffsetFactor * kRadius, sides), kTol);
  if (!made_b.has_value()) {
    print_error(made_b.error());
    return 2;
  }
  const double build_ms = ms_since(build_start);
  std::fprintf(stderr, "[bench] operands built ms=%.3f\n", build_ms);
  std::fflush(stderr);

  const Brep& a = made_a.value();
  const Brep& b = made_b.value();
  double volume_a = 0.0;
  double volume_b = 0.0;
  if (!verify_operand(a, sides, "A", &volume_a)) return 2;
  if (!verify_operand(b, sides, "B", &volume_b)) return 2;
  std::fprintf(stderr, "[bench] volume_a=%.17g volume_b=%.17g\n", volume_a,
               volume_b);
  std::fflush(stderr);

  const std::uint32_t faces_a = a.face_count();
  const std::uint32_t faces_b = b.face_count();

  if (single_op) {
    OpRun run;
    if (!run_op(a, b, *single, false, &run)) return 2;
    print_qa8(run, false, 0.0, true, faces_a, faces_b);
    return run.ms <= kTargetMs && run.topology_ok ? 0 : 1;
  }

  OpRun runs[3];
  for (std::size_t i = 0; i < 3; ++i) {
    if (!run_op(a, b, kOps[i], true, &runs[i])) return 2;
  }

  const OpRun& union_run = runs[0];
  const OpRun& intersection_run = runs[1];
  const OpRun& difference_run = runs[2];
  const bool have_volumes = union_run.volume_ok && intersection_run.volume_ok &&
                            difference_run.volume_ok;
  const double sum_error =
      have_volumes
          ? std::fabs(union_run.volume + intersection_run.volume - volume_a -
                      volume_b) /
                (volume_a + volume_b)
          : 0.0;
  const double difference_error =
      have_volumes ? std::fabs(difference_run.volume + intersection_run.volume -
                               volume_a) /
                         volume_a
                   : 0.0;

  print_qa8(union_run, have_volumes, sum_error, false, faces_a, faces_b);
  print_qa8(intersection_run, have_volumes, sum_error, false, faces_a, faces_b);
  print_qa8(difference_run, have_volumes, difference_error, false, faces_a,
            faces_b);

  bool any_failed = false;
  for (const OpRun& run : runs) {
    if (run.ms > kTargetMs) any_failed = true;
    if (!run.topology_ok) any_failed = true;
  }
  if (!have_volumes || sum_error > kTargetErr ||
      difference_error > kTargetErr) {
    any_failed = true;
  }
  return any_failed ? 1 : 0;
}
