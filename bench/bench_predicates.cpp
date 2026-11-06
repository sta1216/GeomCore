#include <algorithm>
#include <array>
#include <chrono>
#include <cstddef>
#include <cstdint>
#include <cstdio>

#include <geom/numeric/predicates2.hpp>

namespace {

constexpr std::size_t kN = 2000000;
constexpr int kRepetitions = 5;
constexpr double kReleaseBudgetNs = 200.0;

// Triangle (0,0), (1,0), (x + 0.5, y - 511.5) over a 1024 x 1024 grid. The
// vertical offset is centred so that both orientations occur; a predicate that
// returns a constant (or all zeros) then fails the checksum below. The sample
// is computed inline so that no per-call store traffic pollutes the timing.
void run_pass(long long& acc) {
  for (std::size_t i = 0; i < kN; ++i) {
    const double cx = static_cast<double>(i & 1023u) + 0.5;
    const double cy = static_cast<double>((i >> 10) & 1023u) - 511.5;
    acc += geom::v1::sign_orient2d(0.0, 0.0, 1.0, 0.0, cx, cy);
  }
}

// Independent reference sign: doubling every coordinate turns the samples into
// integers, so the orientation sign is the sign of an int64 determinant that
// is computed without the predicate under test.
int reference_sign(std::size_t i) {
  const std::int64_t ax2 = 0;
  const std::int64_t ay2 = 0;
  const std::int64_t bx2 = 2;
  const std::int64_t by2 = 0;
  const std::int64_t cx2 = static_cast<std::int64_t>(2 * (i & 1023u) + 1);
  const std::int64_t cy2 = static_cast<std::int64_t>(2 * ((i >> 10) & 1023u) - 1023);
  const std::int64_t det = (bx2 - ax2) * (cy2 - ay2) - (by2 - ay2) * (cx2 - ax2);
  return det > 0 ? 1 : (det < 0 ? -1 : 0);
}

// One full pass over the fixed input set; the checksum lands in `acc`.

}  // namespace

int main() {
  std::int64_t expected = 0;
  for (std::size_t i = 0; i < kN; ++i) expected += reference_sign(i);

  // Warm-up pass: page-in the code, prime the branch predictor and the caches.
  long long acc = 0;
  run_pass(acc);
  acc = 0;

  std::array<double, kRepetitions> ns{};
  for (int r = 0; r < kRepetitions; ++r) {
    const auto t0 = std::chrono::steady_clock::now();
    run_pass(acc);
    const auto t1 = std::chrono::steady_clock::now();
    ns[static_cast<std::size_t>(r)] =
        std::chrono::duration<double, std::nano>(t1 - t0).count() / static_cast<double>(kN);
  }

  // Read the accumulated checksum once; a volatile read per call is not needed
  // because the comparison below keeps the calls alive.
  volatile long long sink = acc;
  const long long checksum = sink;

  std::array<double, kRepetitions> sorted = ns;
  std::sort(sorted.begin(), sorted.end());
  const double min_ns = sorted.front();
  const double median_ns = sorted[kRepetitions / 2];
  const double max_ns = sorted.back();

  std::printf("orient2d: %.0f calls x %d, min %.2f ns/call, median %.2f, max %.2f\n",
              static_cast<double>(kN), kRepetitions, min_ns, median_ns, max_ns);
  std::printf("checksum: %lld (expected %lld)\n", checksum,
              static_cast<long long>(kRepetitions * expected));

  if (checksum != static_cast<long long>(kRepetitions) * expected) {
    std::printf("FAIL: predicate result does not match the exact reference\n");
    return 2;
  }

#ifdef NDEBUG
  if (min_ns >= kReleaseBudgetNs) {
    std::printf("FAIL: %.2f ns/call exceeds the %.1f ns Release budget\n", min_ns,
                kReleaseBudgetNs);
    return 1;
  }
#else
  std::printf("note: Debug build, Release performance budget not enforced\n");
#endif
  return 0;
}
