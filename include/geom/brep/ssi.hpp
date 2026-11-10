#pragma once

// Public surface-surface intersection diagnostics.
//
// ssi() is declared here and implemented by the intersection driver; this
// header only pins the shared vocabulary so the diagnostics predicates and the
// driver agree on result types. The predicates themselves are inline functions
// in detail/ssi_solver.hpp, next to the point/curve tracing helpers they are
// used with.

#include <array>
#include <cstdint>
#include <vector>

#include <geom/brep/detail/ssi_patch.hpp>
#include <geom/brep/detail/ssi_solver.hpp>
#include <geom/core/result.hpp>

namespace geom {
inline namespace v1 {

// Why an SSI request returned no usable curve, or why its input was rejected
// before tracing. none is the zero value, so an uninitialized issue can never
// masquerade as a classified diagnostic.
enum class SsiDiagnostic : std::uint8_t {
  none = 0,
  tangent_contact,
  coincident,
  degenerate_patch,
  no_convergence,
};

// One diagnostic attached to an SsiResult. params carries the four patch
// parameters most relevant to the issue, in (u_a, v_a, u_b, v_b) order, or
// zeros when the issue is purely structural. message is a static, non-owning
// explanation owned by the implementation.
struct SsiIssue {
  SsiDiagnostic code = SsiDiagnostic::none;
  std::array<double, 4> params{};
  const char* message = "";
};

// Result of one SSI request. curves are the traced intersection curves (empty
// when the request was classified as degenerate before tracing). The counters
// describe the work done and are filled by the driver; the defaults are the
// zero-work state.
struct SsiResult {
  std::vector<SsiCurve> curves;
  std::vector<SsiIssue> issues;
  std::uint64_t candidate_pairs = 0;
  std::uint64_t accepted_segments = 0;
  std::uint64_t rejected_segments = 0;
  double max_residual = 0.0;
};

// Budget for one SSI request. tolerance is the geometric acceptance band,
// max_depth bounds the adaptive patch subdivision, and max_segments bounds the
// traced output. An invalid tolerance or a non-positive budget is rejected by
// ssi() rather than silently clamped.
struct SsiOptions {
  double tolerance = 1e-9;
  int max_depth = 12;
  std::uint64_t max_segments = 100000;
};

// Traces the intersection of two patches. Defined by the intersection driver.
Result<SsiResult> ssi(const SsiPatch& a, const SsiPatch& b, const SsiOptions& options);

}  // namespace v1
}  // namespace geom
