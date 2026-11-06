#include <chrono>
#include <cstdio>

#include <geom/numeric/predicates2.hpp>

int main() {
  constexpr int kN = 2000000;
  volatile int sink = 0;
  const auto t0 = std::chrono::steady_clock::now();
  for (int i = 0; i < kN; ++i) {
    const double x = static_cast<double>(i & 1023);
    const double y = static_cast<double>((i >> 10) & 1023);
    sink += geom::v1::sign_orient2d(0.0, 0.0, 1.0, 0.0, x + 0.5, y + 0.5);
  }
  const auto t1 = std::chrono::steady_clock::now();
  const double ns = std::chrono::duration<double, std::nano>(t1 - t0).count() / kN;
  std::printf("orient2d: %.2f ns/call (sink=%d)\n", ns, sink);
  return ns < 200.0 ? 0 : 1;
}
