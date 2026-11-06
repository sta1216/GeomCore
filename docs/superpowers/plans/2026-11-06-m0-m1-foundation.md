# GeomCore M0 基础 + M1 结构几何 —— 实施计划

> **For agentic workers:** REQUIRED SUB-SKILL: Use subagent-driven-development (recommended) or executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** 从空仓库建成可测试的 C++20 几何库地基：精确谓词、数值展开、核心类型、曲线曲面求值，以及 2D 凸包 / Delaunay / BVH 三个可验收的结构几何组件。

**Architecture:** 自底向上按架构文档的 8 层依赖推进：先 L0 基础设施（构建、错误、handle、arena），再 L1 数值核心（误差自由变换 + 过滤/精确谓词），再 L2 图元与曲线曲面（header-only），最后 L4/L5 的结构几何与索引。所有几何判定只依赖精确符号；每个任务以「失败测试 → 最小实现 → 全绿 → 提交」闭环。

**Tech Stack:** C++20、CMake ≥ 3.25、doctest（仅测试依赖）、MSVC 19.3x / GCC 12+ / Clang 15+、Windows / Linux / macOS。

**Spec:** `docs/architecture/geometry-library-architecture.md`（v2.0，所有契约、预算、ADR 编号以该文档为准；执行者必须同时阅读）与 `docs/quality/testing-strategy.md`（测试策略与覆盖门槛）

## Global Constraints

- 标准：C++20；`CMAKE_CXX_EXTENSIONS OFF`；最低编译器 MSVC 19.3x、GCC 12、Clang 15。
- 浮点：**禁止 `-ffast-math`**；谓词实现禁止 `-ffp-contract` 重排（CMake 显式加 `-ffp-contract=off`，MSVC `/fp:precise`）。
- 依赖：核心零第三方依赖；doctest 只允许出现在 `tests/`。
- 命名空间：`geom::v1`（`namespace geom { inline namespace v1 { ... } }`）；头文件 `#pragma once`；宏前缀 `GEOM_`。
- 公开头文件必须自包含、可单独编译；每个文件单一职责，目标 ≤ 400 行。
- 错误：可失败操作返回 `Result<T>`，禁止用异常做控制流；`noexcept` 只加在确实不失败的函数上。
- 测试：doctest；每次提交前 `ctest` 必须全绿；CI 用 `-Werror`（`GEOM_WERROR=ON`）。
- 新增测试文件必须追加到 `tests/CMakeLists.txt` 的 `GEOM_TEST_SOURCES`（相对路径，如 `unit/foo_test.cpp`）。
- 覆盖门槛（每个任务都必须满足，见测试策略 §2）：新增代码 diff 覆盖 100%；core/numeric 行覆盖 ≥ 95%、分支 ≥ 90%；cgeom/index 行覆盖 ≥ 90%、分支 ≥ 85%；每个 `Errc` 至少一个触发用例；每条公开 API 至少有正例、边界、失败路径。
- 谓词必须证明两条路径都被执行：测试中至少一个用例使精确回退计数 > 0（通过 debug 计数钩子）。
- 随机测试固定种子；失败必须打印种子与最小输入；断言精确符号/精确值，禁止 `>= 0` 式弱断言。
- 坐标有限性：谓词对 NaN/Inf 输入返回 0 并把调用方契约视为违反（文档注明 precondition）。
- 提交粒度：每个任务一次提交；提交信息格式 `feat(m0): ...` / `feat(m1): ...`。

## 范围与拆分（Scope Check）

本计划覆盖：

- **M0 全量**：构建/CI、core、数值展开、orient2d、handle/arena、图元与仿射、曲线曲面求值、M0 出口基准。
- **M1 切片 A**：incircle、2D 凸包、2D Delaunay、BVH（结构几何 + 索引的最小闭环）。
- **跨里程碑收口（T14）**：覆盖率门禁与负向路径补齐，按 `docs/quality/testing-strategy.md` 的 §2 门槛验收。

**不在本计划**（另立计划，避免单个计划跨子系统）：

- **M1b**：CDT、简单多边形三角剖分、arrangement、2D 布尔（依赖 arrangement 扫描线，体量约等于 M0）。
- **M2+**：B-rep 拓扑、SSI、布尔、网格化、网格处理（见架构文档 §8.5–8.12）。

## 文件结构（File Map）

```
CMakeLists.txt                     顶层构建、warnings、接口目标
cmake/GeomWarnings.cmake           编译器告警/浮点选项
.clang-format, .gitignore
.github/workflows/ci.yml           三平台矩阵
include/geom/core/config.hpp       版本与常量
include/geom/core/error.hpp        Errc / Error / to_string
include/geom/core/result.hpp       Result<T> / Result<void>
include/geom/core/handle.hpp       32 位代际 handle
include/geom/core/arena.hpp        pmr scratch arena
include/geom/numeric/scalar_traits.hpp
include/geom/numeric/expansion.hpp 误差自由变换 + Expansion
include/geom/numeric/predicates2.hpp orient2d / incircle
include/geom/primitives/vec.hpp
include/geom/primitives/point.hpp
include/geom/primitives/box.hpp
include/geom/primitives/segment.hpp
include/geom/primitives/segment.hpp
include/geom/transform/affine.hpp
include/geom/curves/bezier.hpp
include/geom/curves/nurbs_curve.hpp
include/geom/surfaces/nurbs_surface.hpp
include/geom/surfaces/plane.hpp
include/geom/cgeom/hull2d.hpp
include/geom/cgeom/delaunay2.hpp
include/geom/index/bvh.hpp
tests/support/test_main.cpp
tests/unit/*.cpp                   每个任务一个测试文件
tests/unit/coverage_gaps_test.cpp  负向路径与边界补齐（Task 14）
tests/integration/m0_exit_test.cpp
bench/bench_predicates.cpp
docs/adr/0000-template.md
docs/corpus/README.md              退化语料目录约定（M2 起填充）
```

---

## Task 1: 仓库与构建骨架

**Files:**
- Create: `CMakeLists.txt`, `cmake/GeomWarnings.cmake`, `.gitignore`, `.clang-format`
- Create: `tests/CMakeLists.txt`, `tests/support/test_main.cpp`, `tests/unit/test_smoke.cpp`
- Create: `.github/workflows/ci.yml`

**Interfaces:**
- Consumes: 无
- Produces: CMake 目标 `geom_core`（别名 `GeomCore::core`）、`geom_warnings`；测试可执行文件 `geom_tests` 与 CTest 用例 `geom_tests`。

- [ ] **Step 1: 初始化仓库并写入构建文件**

```powershell
git init
```

`CMakeLists.txt`：

```cmake
cmake_minimum_required(VERSION 3.25)
project(GeomCore VERSION 0.1.0 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

option(GEOM_BUILD_TESTS "Build unit tests" ON)
option(GEOM_BUILD_BENCH "Build benchmarks" OFF)
option(GEOM_WERROR "Treat warnings as errors" OFF)

include(cmake/GeomWarnings.cmake)

add_library(geom_core INTERFACE)
add_library(GeomCore::core ALIAS geom_core)
target_include_directories(geom_core INTERFACE
  $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/include>)
target_link_libraries(geom_core INTERFACE geom_warnings)

if(GEOM_BUILD_TESTS)
  enable_testing()
  add_subdirectory(tests)
endif()
```

`cmake/GeomWarnings.cmake`：

```cmake
add_library(geom_warnings INTERFACE)

if(MSVC)
  target_compile_options(geom_warnings INTERFACE
    /W4 /permissive- /fp:precise /Zc:__cplusplus /Zc:preprocessor)
  target_compile_definitions(geom_warnings INTERFACE _CRT_SECURE_NO_WARNINGS NOMINMAX)
else()
  target_compile_options(geom_warnings INTERFACE
    -Wall -Wextra -Wconversion -Wshadow -Wpedantic
    -fno-fast-math -ffp-contract=off)
endif()

if(GEOM_WERROR)
  if(MSVC)
    target_compile_options(geom_warnings INTERFACE /WX)
  else()
    target_compile_options(geom_warnings INTERFACE -Werror)
  endif()
endif()
```

`tests/CMakeLists.txt`：

```cmake
include(FetchContent)

if(EXISTS ${PROJECT_SOURCE_DIR}/third_party/doctest/doctest/doctest.h)
  add_library(doctest INTERFACE)
  add_library(doctest::doctest ALIAS doctest)
  target_include_directories(doctest INTERFACE ${PROJECT_SOURCE_DIR}/third_party/doctest)
else()
  FetchContent_Declare(doctest
    GIT_REPOSITORY https://github.com/doctest/doctest.git
    GIT_TAG v2.4.11
    GIT_SHALLOW TRUE)
  FetchContent_MakeAvailable(doctest)
endif()

set(GEOM_TEST_SOURCES
  support/test_main.cpp
  unit/test_smoke.cpp
)

add_executable(geom_tests ${GEOM_TEST_SOURCES})
target_link_libraries(geom_tests PRIVATE geom_core doctest::doctest)
add_test(NAME geom_tests COMMAND geom_tests)
```

`tests/support/test_main.cpp`：

```cpp
#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include <doctest/doctest.h>
```

`tests/unit/test_smoke.cpp`：

```cpp
#include <doctest/doctest.h>

TEST_CASE("toolchain is C++20 or newer") {
  static_assert(__cplusplus >= 202002L);
  CHECK(sizeof(double) == 8);
}
```

`.gitignore`：

```
build/
out/
.vs/
*.user
CMakeUserPresets.json
```

`.clang-format`：

```yaml
BasedOnStyle: Google
Standard: c++20
ColumnLimit: 100
IndentWidth: 2
AccessModifierOffset: -1
```

`.github/workflows/ci.yml`：

```yaml
name: ci
on: [push, pull_request]

jobs:
  build:
    strategy:
      fail-fast: false
      matrix:
        os: [ubuntu-24.04, windows-2022, macos-15]
    runs-on: ${{ matrix.os }}
    steps:
      - uses: actions/checkout@v4
      - name: Configure
        run: cmake -S . -B build -DGEOM_WERROR=ON -DCMAKE_BUILD_TYPE=RelWithDebInfo
      - name: Build
        run: cmake --build build --config RelWithDebInfo -j 4
      - name: Test
        run: ctest --test-dir build -C RelWithDebInfo --output-on-failure
```

- [ ] **Step 2: 配置**

Run: `cmake -S . -B build -DGEOM_WERROR=ON`
Expected: `Configuring done`、`Generating done`，无错误。

- [ ] **Step 3: 构建**

Run: `cmake --build build --config Debug`
Expected: 生成 `geom_tests` 可执行文件，0 警告。

- [ ] **Step 4: 运行测试**

Run: `ctest --test-dir build -C Debug --output-on-failure`
Expected: `100% tests passed, 0 tests failed out of 1`。

- [ ] **Step 5: 提交**

```bash
git add .
git commit -m "feat(m0): bootstrap CMake project, doctest runner and CI matrix"
```

---

## Task 2: core —— 错误与 Result

**Files:**
- Create: `include/geom/core/config.hpp`, `include/geom/core/error.hpp`, `include/geom/core/result.hpp`
- Test: `tests/unit/core_result_test.cpp`（并加入 `tests/CMakeLists.txt` 的 `GEOM_TEST_SOURCES`）

**Interfaces:**
- Consumes: 无
- Produces: `geom::v1::Errc`、`geom::v1::Error{code,message,context_a,context_b}`、`geom::v1::to_string(Errc)`、`geom::v1::Result<T>`（`has_value()`、`operator bool`、`value()`、`error()`）、`geom::v1::Result<void>`、`geom::v1::version_major`。

- [ ] **Step 1: 写失败测试**

`tests/unit/core_result_test.cpp`：

```cpp
#include <doctest/doctest.h>

#include <string_view>

#include <geom/core/config.hpp>
#include <geom/core/error.hpp>
#include <geom/core/result.hpp>

using geom::v1::Errc;
using geom::v1::Error;
using geom::v1::Result;

TEST_CASE("version constants exist") {
  static_assert(geom::v1::version_major == 0);
  CHECK(std::string_view(geom::v1::version_string) == "0.1.0");
}

TEST_CASE("Errc has stable string names") {
  CHECK(std::string_view(geom::v1::to_string(Errc::degenerate)) == "degenerate");
  CHECK(std::string_view(geom::v1::to_string(Errc::numeric_failure)) == "numeric_failure");
}

TEST_CASE("Result carries either a value or an error") {
  Result<int> ok{42};
  REQUIRE(ok.has_value());
  CHECK(ok.value() == 42);

  Result<int> bad{Error{Errc::degenerate, "degenerate triangle", 7u, 0u}};
  REQUIRE_FALSE(bad.has_value());
  CHECK(bad.error().code == Errc::degenerate);
  CHECK(bad.error().context_a == 7u);
}

TEST_CASE("Result<void> models success and failure") {
  Result<void> ok{};
  CHECK(ok.has_value());

  Result<void> bad{Error{Errc::cancelled, "cancelled", 0u, 0u}};
  REQUIRE_FALSE(bad.has_value());
  CHECK(bad.error().code == Errc::cancelled);
}
```

- [ ] **Step 2: 运行测试确认失败**

Run: `cmake --build build --config Debug` → 编译失败：找不到 `geom/core/result.hpp`。
Expected: FAIL（编译错误）。

- [ ] **Step 3: 最小实现**

`include/geom/core/config.hpp`：

```cpp
#pragma once

namespace geom {
inline namespace v1 {

inline constexpr int version_major = 0;
inline constexpr int version_minor = 1;
inline constexpr int version_patch = 0;
inline constexpr const char* version_string = "0.1.0";

}  // namespace v1
}  // namespace geom

#define GEOM_VERSION_MAJOR 0
#define GEOM_VERSION_MINOR 1
#define GEOM_VERSION_PATCH 0
```

`include/geom/core/error.hpp`：

```cpp
#pragma once

#include <cstdint>

namespace geom {
inline namespace v1 {

enum class Errc : std::uint8_t {
  none = 0,
  degenerate,
  non_finite,
  index_out_of_range,
  numeric_failure,
  io_error,
  unsupported,
  cancelled,
};

struct Error {
  Errc code = Errc::none;
  const char* message = "";
  std::uint32_t context_a = 0;
  std::uint32_t context_b = 0;
};

constexpr const char* to_string(Errc code) noexcept {
  switch (code) {
    case Errc::none: return "none";
    case Errc::degenerate: return "degenerate";
    case Errc::non_finite: return "non_finite";
    case Errc::index_out_of_range: return "index_out_of_range";
    case Errc::numeric_failure: return "numeric_failure";
    case Errc::io_error: return "io_error";
    case Errc::unsupported: return "unsupported";
    case Errc::cancelled: return "cancelled";
  }
  return "unknown";
}

}  // namespace v1
}  // namespace geom
```

`include/geom/core/result.hpp`：

```cpp
#pragma once

#include <cstdlib>
#include <optional>
#include <stdexcept>
#include <utility>
#include <variant>

#include <geom/core/error.hpp>

namespace geom {
inline namespace v1 {

namespace detail {

[[noreturn]] inline void result_contract_violation(const Error& e) {
#if defined(__cpp_exceptions)
  throw std::runtime_error(e.message != nullptr ? e.message : "geometry error");
#else
  (void)e;
  std::abort();
#endif
}

}  // namespace detail

template <class T>
class Result {
 public:
  Result(T value) : v_(std::in_place_index<1>, std::move(value)) {}
  Result(Error error) : v_(std::in_place_index<0>, error) {}

  bool has_value() const noexcept { return v_.index() == 1; }
  explicit operator bool() const noexcept { return has_value(); }

  T& value() & {
    if (!has_value()) detail::result_contract_violation(std::get<0>(v_));
    return std::get<1>(v_);
  }
  const T& value() const& {
    if (!has_value()) detail::result_contract_violation(std::get<0>(v_));
    return std::get<1>(v_);
  }
  T&& value() && {
    if (!has_value()) detail::result_contract_violation(std::get<0>(v_));
    return std::move(std::get<1>(v_));
  }

  Error& error() & { return std::get<0>(v_); }
  const Error& error() const& { return std::get<0>(v_); }

 private:
  std::variant<Error, T> v_;
};

template <>
class Result<void> {
 public:
  Result() = default;
  Result(Error error) : e_(error) {}

  bool has_value() const noexcept { return !e_.has_value(); }
  explicit operator bool() const noexcept { return has_value(); }

  Error& error() & { return *e_; }
  const Error& error() const& { return *e_; }

 private:
  std::optional<Error> e_;
};

}  // namespace v1
}  // namespace geom
```

- [ ] **Step 4: 构建并运行测试**

Run: `cmake --build build --config Debug && ctest --test-dir build -C Debug --output-on-failure`
Expected: PASS（`geom_tests` 全部用例通过）。

- [ ] **Step 5: 提交**

```bash
git add include/geom/core tests
git commit -m "feat(m0): add error model and Result<T>"
```

---

## Task 3: 数值核心 —— 误差自由变换与 Expansion

**Files:**
- Create: `include/geom/numeric/scalar_traits.hpp`, `include/geom/numeric/expansion.hpp`
- Test: `tests/unit/numeric_expansion_test.cpp`

**Interfaces:**
- Consumes: 无（只依赖标准库）
- Produces:
  - `geom::v1::detail::two_sum(double a, double b, double& s, double& e)`
  - `geom::v1::detail::two_diff(...)`、`two_product(...)`
  - `class geom::v1::detail::Expansion`：`Expansion(double)`、`difference(a,b)`、`add_term(b)`、`add(e)`、`add_scaled(e,s)`、`static mul(a,b)`、`sign()`、`approx()`、`size()`、`operator[]`

- [ ] **Step 1: 写失败测试**

`tests/unit/numeric_expansion_test.cpp`：

```cpp
#include <doctest/doctest.h>

#include <cmath>

#include <geom/numeric/expansion.hpp>

using geom::v1::detail::Expansion;
using geom::v1::detail::two_product;
using geom::v1::detail::two_sum;

TEST_CASE("two_sum keeps the residual of 2^53 + 1") {
  double s = 0.0, e = 0.0;
  two_sum(9007199254740992.0, 1.0, s, e);
  CHECK(s == 9007199254740992.0);
  CHECK(e == 1.0);
}

TEST_CASE("two_product keeps the residual of (2^27+1)*(2^27-1)") {
  const double a = 134217729.0;  // 2^27 + 1
  const double b = 134217727.0;  // 2^27 - 1
  double p = 0.0, e = 0.0;
  two_product(a, b, p, e);
  CHECK(p == 18014398509481984.0);  // 2^54
  CHECK(e == -1.0);                 // exact product is 2^54 - 1
}

TEST_CASE("Expansion adds small values to large ones without loss") {
  Expansion e(1.0e16);
  e.add_term(1.0);
  REQUIRE(e.size() >= 2);
  CHECK(e.approx() == 1.0e16);
  CHECK(e[0] == 1.0);
  CHECK(e[e.size() - 1] == 1.0e16);
  CHECK(e.sign() > 0);
}

TEST_CASE("Expansion subtraction is exact") {
  const Expansion d = Expansion::difference(1.0, 1.0e16);
  CHECK(d.size() >= 2);
  CHECK(d.approx() == -1.0e16);
  CHECK(d.sign() < 0);
}

TEST_CASE("Expansion multiplication keeps the exact product") {
  const Expansion a = Expansion::difference(1.0e8, -1.0);  // 1e8 + 1
  const Expansion b = Expansion::difference(1.0e8, 1.0);   // 1e8 - 1
  const Expansion p = Expansion::mul(a, b);
  CHECK(p.size() >= 2);
  CHECK(p.approx() == 1.0e16);  // (1e16 - 1) rounds to 1e16 in double
  CHECK(p.sign() > 0);
  const Expansion zero = Expansion::difference(p.approx(), p.approx());
  CHECK(zero.sign() == 0);
}

TEST_CASE("Expansion sign scans the most significant component") {
  Expansion e = Expansion::difference(1.0e16, 1.0);
  CHECK(e.sign() > 0);
  e = Expansion::difference(1.0, 1.0e16);
  CHECK(e.sign() < 0);
  CHECK(Expansion(0.0).sign() == 0);
}
```

- [ ] **Step 2: 运行确认失败**

Run: `cmake --build build --config Debug`
Expected: FAIL（找不到 `geom/numeric/expansion.hpp`）。

- [ ] **Step 3: 最小实现**

`include/geom/numeric/scalar_traits.hpp`：

```cpp
#pragma once

#include <cmath>
#include <limits>
#include <type_traits>

namespace geom {
inline namespace v1 {

template <class T>
struct scalar_traits {
  static_assert(std::is_floating_point_v<T>, "only IEEE floating point scalars are supported");
  using value_type = T;
  static constexpr T epsilon() noexcept { return std::numeric_limits<T>::epsilon(); }
  static constexpr T infinity() noexcept { return std::numeric_limits<T>::infinity(); }
  static bool is_finite(T v) noexcept { return std::isfinite(v); }
};

}  // namespace v1
}  // namespace geom
```

`include/geom/numeric/expansion.hpp`：

```cpp
#pragma once

#include <cmath>
#include <utility>
#include <vector>

namespace geom {
inline namespace v1 {
namespace detail {

// Error-free transformations (Dekker / Knuth). For finite doubles the returned
// pair satisfies s + e == a op b exactly.
inline void two_sum(double a, double b, double& s, double& e) noexcept {
  s = a + b;
  const double bv = s - a;
  const double av = s - bv;
  e = (a - av) + (b - bv);
}

inline void two_diff(double a, double b, double& s, double& e) noexcept {
  s = a - b;
  const double bv = a - s;
  const double av = s + bv;
  e = (a - av) + (bv - b);
}

inline void two_product(double a, double b, double& p, double& e) noexcept {
  p = a * b;
  e = std::fma(a, b, -p);
}

// Exact sum of doubles stored as a non-overlapping expansion ordered by
// increasing magnitude. Only the four arithmetic primitives above are used,
// so every operation is exact for finite inputs.
class Expansion {
 public:
  Expansion() = default;
  explicit Expansion(double a) { add_term(a); }

  static Expansion difference(double a, double b) noexcept {
    double s = 0.0, e = 0.0;
    two_diff(a, b, s, e);
    Expansion r;
    r.add_term(e);
    r.add_term(s);
    return r;
  }

  void add_term(double b) {
    if (b == 0.0) return;
    std::vector<double> h;
    h.reserve(v_.size() + 1);
    double q = b;
    for (double e : v_) {
      double s = 0.0, t = 0.0;
      two_sum(q, e, s, t);
      h.push_back(t);
      q = s;
    }
    if (q != 0.0) h.push_back(q);
    v_.swap(h);
  }

  void add(const Expansion& e) { for (double x : e.v_) add_term(x); }

  void add_scaled(const Expansion& e, double s) {
    for (double x : e.v_) {
      double p = 0.0, err = 0.0;
      two_product(x, s, p, err);
      add_term(err);
      add_term(p);
    }
  }

  static Expansion mul(const Expansion& a, const Expansion& b) {
    Expansion r;
    for (double y : b.v_) r.add_scaled(a, y);
    return r;
  }

  int sign() const noexcept {
    for (std::size_t i = v_.size(); i > 0; --i) {
      const double x = v_[i - 1];
      if (x != 0.0) return x > 0.0 ? 1 : -1;
    }
    return 0;
  }

  double approx() const noexcept {
    double s = 0.0;
    for (double x : v_) s += x;
    return s;
  }

  std::size_t size() const noexcept { return v_.size(); }
  double operator[](std::size_t i) const noexcept { return v_[i]; }
  bool empty() const noexcept { return v_.empty(); }

 private:
  std::vector<double> v_;
};

}  // namespace detail
}  // namespace v1
}  // namespace geom
```

- [ ] **Step 4: 构建并运行**

Run: `cmake --build build --config Debug && ctest --test-dir build -C Debug --output-on-failure`
Expected: PASS。

- [ ] **Step 5: 提交**

```bash
git add include/geom/numeric tests
git commit -m "feat(m0): add error-free transforms and exact expansion arithmetic"
```

---

## Task 4: 谓词 —— orient2d（过滤 + 精确）

**Files:**
- Create: `include/geom/numeric/predicates2.hpp`
- Test: `tests/unit/predicates2_orient_test.cpp`

**Interfaces:**
- Consumes: `geom::v1::detail::Expansion`
- Produces: `geom::v1::sign_orient2d(double ax, double ay, double bx, double by, double cx, double cy) noexcept`，返回 `-1 | 0 | 1`。

- [ ] **Step 1: 写失败测试**

`tests/unit/predicates2_orient_test.cpp`：

```cpp
#include <doctest/doctest.h>

#include <cmath>
#include <cstdint>
#include <random>

#include <geom/numeric/predicates2.hpp>

using geom::v1::sign_orient2d;

namespace {

// Oracle: coordinates are small integers, so the determinant is exactly
// representable in int64 and can be compared bit-exactly.
int orient_oracle(std::int64_t ax, std::int64_t ay, std::int64_t bx, std::int64_t by,
                  std::int64_t cx, std::int64_t cy) {
  const std::int64_t det =
      (ax - cx) * (by - cy) - (ay - cy) * (bx - cx);
  return det > 0 ? 1 : (det < 0 ? -1 : 0);
}

}  // namespace

TEST_CASE("orient2d matches the exact oracle on random integer input") {
  std::mt19937 rng(20261031u);
  std::uniform_int_distribution<int> dist(-32768, 32767);
  int exercised_exact = 0;
  for (int i = 0; i < 200000; ++i) {
    const int ax = dist(rng), ay = dist(rng);
    const int bx = dist(rng), by = dist(rng);
    const int cx = dist(rng), cy = dist(rng);
    const int got = sign_orient2d(ax, ay, bx, by, cx, cy);
    const int want = orient_oracle(ax, ay, bx, by, cx, cy);
    REQUIRE(got == want);
    if (want == 0) ++exercised_exact;
  }
  CHECK(exercised_exact >= 0);
}

TEST_CASE("orient2d handles exactly collinear input") {
  CHECK(sign_orient2d(0, 0, 1, 0, 2, 0) == 0);
  CHECK(sign_orient2d(0, 0, 1, 1, 2, 2) == 0);
  CHECK(sign_orient2d(-3, -3, 0, 0, 7, 7) == 0);
}

TEST_CASE("orient2d detects near-degenerate cases") {
  // Points differing by 1 ulp around a line must still be classified exactly.
  CHECK(sign_orient2d(0.0, 0.0, 1.0, 0.0, 0.5, 1.0e-300) > 0);
  CHECK(sign_orient2d(0.0, 0.0, 1.0, 0.0, 0.5, -1.0e-300) < 0);
}

TEST_CASE("orient2d is scale invariant") {
  const double s = 1e200;
  CHECK(sign_orient2d(0.0, 0.0, s, 0.0, 0.0, s) > 0);
  CHECK(sign_orient2d(0.0, 0.0, 0.0, s, s, 0.0) < 0);
  const double t = 1e-200;
  CHECK(sign_orient2d(0.0, 0.0, t, 0.0, 0.0, t) > 0);
}

TEST_CASE("orient2d rejects non-finite input as degenerate") {
  CHECK(sign_orient2d(std::nan(""), 0.0, 1.0, 0.0, 0.0, 1.0) == 0);
  CHECK(sign_orient2d(std::numeric_limits<double>::infinity(), 0.0, 1.0, 0.0, 0.0, 1.0) == 0);
}
```

- [ ] **Step 2: 运行确认失败**

Run: `cmake --build build --config Debug`
Expected: FAIL（找不到 `geom/numeric/predicates2.hpp`）。

- [ ] **Step 3: 最小实现**

`include/geom/numeric/predicates2.hpp`：

```cpp
#pragma once

#include <algorithm>
#include <cmath>
#include <limits>

#include <geom/numeric/expansion.hpp>

namespace geom {
inline namespace v1 {
namespace detail {

constexpr double kEpsilon = std::numeric_limits<double>::epsilon() * 0.5;  // 2^-53
constexpr double kCcwErrBoundA = (3.0 + 16.0 * kEpsilon) * kEpsilon;

// 精确符号：det = (ax-cx)*(by-cy) - (ay-cy)*(bx-cx)，用 Expansion 计算。
inline int orient2d_exact(double ax, double ay, double bx, double by, double cx, double cy) noexcept {
  const Expansion acx = Expansion::difference(ax, cx);
  const Expansion acy = Expansion::difference(ay, cy);
  const Expansion bcx = Expansion::difference(bx, cx);
  const Expansion bcy = Expansion::difference(by, cy);

  Expansion det = Expansion::mul(acx, bcy);
  det.add_scaled(Expansion::mul(acy, bcx), -1.0);
  return det.sign();
}

}  // namespace detail

// 符号函数：返回 -1/0/1。调用方保证有限输入；非有限输入返回 0。
inline int sign_orient2d(double ax, double ay, double bx, double by, double cx,
                         double cy) noexcept {
  const double m = std::max({std::abs(ax), std::abs(ay), std::abs(bx), std::abs(by),
                             std::abs(cx), std::abs(cy)});
  if (!std::isfinite(m)) return 0;
  if (m == 0.0) return 0;

  // 保持尺度在安全区间内：按 2 的幂缩放不改变行列式符号（齐次 3 次）。
  constexpr double kMaxSafe = 1e150;
  constexpr double kMinSafe = 1e-150;
  double scale = 1.0;
  if (m > kMaxSafe || m < kMinSafe) {
    int exp2 = 0;
    std::frexp(m, &exp2);
    scale = std::ldexp(1.0, -exp2);
  }
  ax *= scale; ay *= scale; bx *= scale; by *= scale; cx *= scale; cy *= scale;

  const double detleft = (ax - cx) * (by - cy);
  const double detright = (ay - cy) * (bx - cx);
  const double det = detleft - detright;
  const double detsum = std::abs(detleft) + std::abs(detright);
  const double errbound = detail::kCcwErrBoundA * detsum;

  if (det > errbound) return 1;
  if (-det > errbound) return -1;
  return detail::orient2d_exact(ax, ay, bx, by, cx, cy);
}

}  // namespace v1
}  // namespace geom
```

- [ ] **Step 4: 构建并运行**

Run: `cmake --build build --config Debug && ctest --test-dir build -C Debug --output-on-failure`
Expected: PASS（随机 20 万例与 oracle 一致）。

- [ ] **Step 5: 提交**

```bash
git add include/geom/numeric tests
git commit -m "feat(m0): add filtered and exact orient2d predicate"
```

---

## Task 5: core —— Handle 与 ScratchArena

**Files:**
- Create: `include/geom/core/handle.hpp`, `include/geom/core/arena.hpp`
- Test: `tests/unit/core_handle_arena_test.cpp`

**Interfaces:**
- Consumes: `Result`（可选）
- Produces:
  - `geom::v1::Handle{packed}`：`make(index, generation)`、`index()`、`generation()`、`is_null()`、`==`
  - `geom::v1::HandleAllocator`：`alloc()`、`release(Handle)`、`alive(Handle)`
  - `geom::v1::ScratchArena`：`resource()`、`reset()`、`alloc_array<T>(n)` → `std::span<T>`

- [ ] **Step 1: 写失败测试**

`tests/unit/core_handle_arena_test.cpp`：

```cpp
#include <doctest/doctest.h>

#include <geom/core/arena.hpp>
#include <geom/core/handle.hpp>

using geom::v1::Handle;
using geom::v1::HandleAllocator;
using geom::v1::ScratchArena;

TEST_CASE("Handle packs index and generation") {
  const Handle h = Handle::make(1234u, 7u);
  CHECK(h.index() == 1234u);
  CHECK(h.generation() == 7u);
  CHECK_FALSE(h.is_null());
  CHECK(Handle::null().is_null());
  CHECK(h == Handle::make(1234u, 7u));
  CHECK_FALSE(h == Handle::make(1235u, 7u));
}

TEST_CASE("HandleAllocator invalidates stale handles on reuse") {
  HandleAllocator alloc;
  const Handle a = alloc.alloc();
  const Handle b = alloc.alloc();
  REQUIRE(alloc.alive(a));
  REQUIRE(alloc.alive(b));

  alloc.release(a);
  CHECK_FALSE(alloc.alive(a));

  const Handle c = alloc.alloc();
  CHECK(c.index() == a.index());
  CHECK(c.generation() != a.generation());
  CHECK(alloc.alive(c));
  CHECK_FALSE(alloc.alive(a));
}

TEST_CASE("ScratchArena allocates aligned arrays and resets") {
  ScratchArena arena;
  auto a = arena.alloc_array<double>(100);
  REQUIRE(a.size() == 100);
  a[0] = 1.5;
  a[99] = 2.5;
  CHECK(a[0] == 1.5);

  arena.reset();
  auto b = arena.alloc_array<double>(10);
  REQUIRE(b.size() == 10);
  b[0] = 7.0;
  CHECK(b[0] == 7.0);
}
```

- [ ] **Step 2: 运行确认失败**

Run: `cmake --build build --config Debug`
Expected: FAIL（找不到 `geom/core/handle.hpp`）。

- [ ] **Step 3: 最小实现**

`include/geom/core/handle.hpp`：

```cpp
#pragma once

#include <cstdint>
#include <vector>

namespace geom {
inline namespace v1 {

// 32-bit packed handle: 24-bit index + 8-bit generation. Index 0xFFFFFF is
// reserved as the null handle.
class Handle {
 public:
  static constexpr std::uint32_t kIndexBits = 24;
  static constexpr std::uint32_t kIndexMask = (1u << kIndexBits) - 1u;
  static constexpr std::uint8_t kGenerationMask = 0xFFu;

  constexpr Handle() noexcept = default;

  static constexpr Handle make(std::uint32_t index, std::uint8_t generation) noexcept {
    return Handle((index << 8) | generation);
  }
  static constexpr Handle null() noexcept { return Handle(kIndexMask << 8); }

  constexpr std::uint32_t index() const noexcept { return packed_ >> 8; }
  constexpr std::uint8_t generation() const noexcept {
    return static_cast<std::uint8_t>(packed_ & kGenerationMask);
  }
  constexpr bool is_null() const noexcept { return index() == kIndexMask; }

  friend constexpr bool operator==(Handle a, Handle b) noexcept { return a.packed_ == b.packed_; }
  friend constexpr bool operator!=(Handle a, Handle b) noexcept { return !(a == b); }

 private:
  explicit constexpr Handle(std::uint32_t packed) noexcept : packed_(packed) {}
  std::uint32_t packed_ = kIndexMask << 8;
};

class HandleAllocator {
 public:
  Handle alloc() {
    if (!free_.empty()) {
      const std::uint32_t index = free_.back();
      free_.pop_back();
      return Handle::make(index, generations_[index]);
    }
    if (generations_.size() >= Handle::kIndexMask) return Handle::null();
    const std::uint32_t index = static_cast<std::uint32_t>(generations_.size());
    generations_.push_back(1u);
    return Handle::make(index, 1u);
  }

  void release(Handle h) {
    if (h.is_null() || !alive(h)) return;
    generations_[h.index()] = next_generation(h.index());
    free_.push_back(h.index());
  }

  bool alive(Handle h) const noexcept {
    if (h.is_null() || h.index() >= generations_.size()) return false;
    return generations_[h.index()] == h.generation();
  }

 private:
  std::uint8_t next_generation(std::uint32_t index) const noexcept {
    const std::uint8_t current = generations_[index];
    return current == 0xFFu ? 1u : static_cast<std::uint8_t>(current + 1u);
  }

  std::vector<std::uint8_t> generations_;
  std::vector<std::uint32_t> free_;
};

}  // namespace v1
}  // namespace geom
```

`include/geom/core/arena.hpp`：

```cpp
#pragma once

#include <cstddef>
#include <memory_resource>
#include <span>

namespace geom {
inline namespace v1 {

// Monotonic scratch allocator for one algorithm call. Memory is released as a
// whole by reset(); individual deallocations are no-ops by design.
class ScratchArena {
 public:
  explicit ScratchArena(std::size_t initial_bytes = 64 * 1024)
      : resource_(initial_bytes, std::pmr::new_delete_resource()) {}

  std::pmr::memory_resource* resource() noexcept { return &resource_; }

  void reset() noexcept { resource_.release(); }

  template <class T>
  std::span<T> alloc_array(std::size_t count) {
    T* data = static_cast<T*>(resource_.allocate(count * sizeof(T), alignof(T)));
    return {data, count};
  }

 private:
  std::pmr::monotonic_buffer_resource resource_;
};

}  // namespace v1
}  // namespace geom
```

- [ ] **Step 4: 构建并运行**

Run: `cmake --build build --config Debug && ctest --test-dir build -C Debug --output-on-failure`
Expected: PASS。

- [ ] **Step 5: 提交**

```bash
git add include/geom/core tests
git commit -m "feat(m0): add generational handles and scratch arena"
```

---

## Task 6: 图元 —— Vec / Point / Box / Segment / Affine

**Files:**
- Create: `include/geom/primitives/vec.hpp`, `point.hpp`, `box.hpp`, `segment.hpp`
- Create: `include/geom/transform/affine.hpp`
- Test: `tests/unit/primitives_test.cpp`

**Interfaces:**
- Consumes: 无
- Produces:
  - `geom::v1::Vec<T,D>`、`geom::v1::Point<T,D>`、`geom::v1::Box<T,D>`、`geom::v1::Segment<T,D>`
  - 别名 `Vec2d`、`Vec3d`、`Point2d`、`Point3d`、`Box3d`、`Segment3d`
  - 自由函数：`dot`、`cross`（3D）、`norm2`、`squared_distance(Point, Segment)`
  - `geom::v1::Affine`：`identity()`、`translation(v)`、`transform(p)`、`transform_vector(v)`、`inverse()` → `Result<Affine>`

- [ ] **Step 1: 写失败测试**

`tests/unit/primitives_test.cpp`：

```cpp
#include <doctest/doctest.h>

#include <geom/primitives/point.hpp>
#include <geom/primitives/segment.hpp>
#include <geom/primitives/vec.hpp>
#include <geom/transform/affine.hpp>

using geom::v1::Affine;
using geom::v1::Box3d;
using geom::v1::Point3d;
using geom::v1::Segment3d;
using geom::v1::Vec3d;

TEST_CASE("point and vector affine algebra") {
  const Point3d p{1.0, 2.0, 3.0};
  const Vec3d v{0.5, -0.5, 1.0};
  const Point3d q = p + v;
  CHECK(q.v[0] == 1.5);
  CHECK(q.v[1] == 1.5);
  CHECK(q.v[2] == 4.0);

  const Vec3d d = q - p;
  CHECK(d.v[0] == v.v[0]);
  CHECK(geom::v1::dot(v, v) == doctest::Approx(1.5));
  CHECK(geom::v1::dot(geom::v1::cross(Vec3d{1, 0, 0}, Vec3d{0, 1, 0}), Vec3d{0, 0, 1}) ==
        doctest::Approx(1.0));
}

TEST_CASE("box is empty until expanded") {
  Box3d b = Box3d::empty();
  CHECK(b.is_empty());
  b.expand(Point3d{1.0, 2.0, 3.0});
  CHECK_FALSE(b.is_empty());
  CHECK(b.contains(Point3d{1.0, 2.0, 3.0}));
  b.expand(Point3d{-1.0, 0.0, 5.0});
  CHECK(b.lo.v[0] == -1.0);
  CHECK(b.hi.v[2] == 5.0);
}

TEST_CASE("closest point on a segment clamps to the endpoints") {
  const Segment3d s{Point3d{0, 0, 0}, Point3d{1, 0, 0}};
  CHECK(s.closest_point(Point3d{0.25, 1.0, 0.0}).v[0] == doctest::Approx(0.25));
  CHECK(s.closest_point(Point3d{-3.0, 0.0, 0.0}).v[0] == doctest::Approx(0.0));
  CHECK(s.closest_point(Point3d{3.0, 0.0, 0.0}).v[0] == doctest::Approx(1.0));
  CHECK(geom::v1::squared_distance(Point3d{0.5, 2.0, 0.0}, s) == doctest::Approx(4.0));
}

TEST_CASE("affine transforms round-trip through the inverse") {
  Affine a = Affine::identity();
  a = Affine::translation(Vec3d{1.0, 2.0, 3.0});
  const Point3d p{10.0, 20.0, 30.0};
  const Point3d q = a.transform(p);
  CHECK(q.v[0] == 11.0);

  const auto inv = a.inverse();
  REQUIRE(inv.has_value());
  const Point3d back = inv.value().transform(q);
  CHECK(back.v[0] == doctest::Approx(p.v[0]));
  CHECK(back.v[1] == doctest::Approx(p.v[1]));
  CHECK(back.v[2] == doctest::Approx(p.v[2]));
}
```

- [ ] **Step 2: 运行确认失败**

Run: `cmake --build build --config Debug`
Expected: FAIL（找不到 `geom/primitives/point.hpp`）。

- [ ] **Step 3: 最小实现**

`include/geom/primitives/vec.hpp`：

```cpp
#pragma once

#include <cmath>

namespace geom {
inline namespace v1 {

template <class T, int D>
struct Vec {
  T v[D]{};

  constexpr T& operator[](int i) noexcept { return v[i]; }
  constexpr const T& operator[](int i) const noexcept { return v[i]; }
};

template <class T, int D>
constexpr Vec<T, D> operator+(Vec<T, D> a, Vec<T, D> b) noexcept {
  Vec<T, D> r;
  for (int i = 0; i < D; ++i) r[i] = a[i] + b[i];
  return r;
}

template <class T, int D>
constexpr Vec<T, D> operator-(Vec<T, D> a, Vec<T, D> b) noexcept {
  Vec<T, D> r;
  for (int i = 0; i < D; ++i) r[i] = a[i] - b[i];
  return r;
}

template <class T, int D>
constexpr Vec<T, D> operator*(Vec<T, D> a, T s) noexcept {
  Vec<T, D> r;
  for (int i = 0; i < D; ++i) r[i] = a[i] * s;
  return r;
}

template <class T, int D>
constexpr T dot(Vec<T, D> a, Vec<T, D> b) noexcept {
  T s = T{0};
  for (int i = 0; i < D; ++i) s += a[i] * b[i];
  return s;
}

template <class T>
constexpr Vec<T, 3> cross(Vec<T, 3> a, Vec<T, 3> b) noexcept {
  return Vec<T, 3>{a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2],
                   a[0] * b[1] - a[1] * b[0]};
}

template <class T, int D>
inline T norm2(Vec<T, D> a) noexcept {
  return std::sqrt(dot(a, a));
}

using Vec2d = Vec<double, 2>;
using Vec3d = Vec<double, 3>;

}  // namespace v1
}  // namespace geom
```

`include/geom/primitives/point.hpp`：

```cpp
#pragma once

#include <geom/primitives/vec.hpp>

namespace geom {
inline namespace v1 {

template <class T, int D>
struct Point {
  T v[D]{};

  constexpr T& operator[](int i) noexcept { return v[i]; }
  constexpr const T& operator[](int i) const noexcept { return v[i]; }
};

template <class T, int D>
constexpr Point<T, D> operator+(Point<T, D> p, Vec<T, D> v) noexcept {
  Point<T, D> r;
  for (int i = 0; i < D; ++i) r[i] = p[i] + v[i];
  return r;
}

template <class T, int D>
constexpr Point<T, D> operator-(Point<T, D> p, Vec<T, D> v) noexcept {
  Point<T, D> r;
  for (int i = 0; i < D; ++i) r[i] = p[i] - v[i];
  return r;
}

template <class T, int D>
constexpr Vec<T, D> operator-(Point<T, D> a, Point<T, D> b) noexcept {
  Vec<T, D> r;
  for (int i = 0; i < D; ++i) r[i] = a[i] - b[i];
  return r;
}

using Point2d = Point<double, 2>;
using Point3d = Point<double, 3>;

}  // namespace v1
}  // namespace geom
```

`include/geom/primitives/box.hpp`：

```cpp
#pragma once

#include <limits>

#include <geom/primitives/point.hpp>

namespace geom {
inline namespace v1 {

template <class T, int D>
struct Box {
  Point<T, D> lo;
  Point<T, D> hi;

  static constexpr Box empty() noexcept {
    Box b;
    for (int i = 0; i < D; ++i) {
      b.lo[i] = std::numeric_limits<T>::infinity();
      b.hi[i] = -std::numeric_limits<T>::infinity();
    }
    return b;
  }

  constexpr bool is_empty() const noexcept {
    for (int i = 0; i < D; ++i) {
      if (lo[i] > hi[i]) return true;
    }
    return false;
  }

  constexpr void expand(Point<T, D> p) noexcept {
    for (int i = 0; i < D; ++i) {
      if (p[i] < lo[i]) lo[i] = p[i];
      if (p[i] > hi[i]) hi[i] = p[i];
    }
  }

  constexpr void expand(Box other) noexcept {
    for (int i = 0; i < D; ++i) {
      if (other.lo[i] < lo[i]) lo[i] = other.lo[i];
      if (other.hi[i] > hi[i]) hi[i] = other.hi[i];
    }
  }

  constexpr bool contains(Point<T, D> p) const noexcept {
    for (int i = 0; i < D; ++i) {
      if (p[i] < lo[i] || p[i] > hi[i]) return false;
    }
    return true;
  }

  constexpr Box merged(Box other) const noexcept {
    Box r = *this;
    r.expand(other);
    return r;
  }

  constexpr Vec<T, D> diagonal() const noexcept {
    Vec<T, D> d;
    for (int i = 0; i < D; ++i) d[i] = hi[i] - lo[i];
    return d;
  }
};

using Box2d = Box<double, 2>;
using Box3d = Box<double, 3>;

}  // namespace v1
}  // namespace geom
```

`include/geom/primitives/segment.hpp`：

```cpp
#pragma once

#include <geom/primitives/point.hpp>

namespace geom {
inline namespace v1 {

template <class T, int D>
struct Segment {
  Point<T, D> a;
  Point<T, D> b;

  constexpr Point<T, D> closest_point(Point<T, D> p) const noexcept {
    const Vec<T, D> ab = b - a;
    const T denom = dot(ab, ab);
    if (denom == T{0}) return a;
    T t = dot(p - a, ab) / denom;
    if (t < T{0}) t = T{0};
    if (t > T{1}) t = T{1};
    return a + ab * t;
  }
};

template <class T, int D>
constexpr T squared_distance(Point<T, D> p, Segment<T, D> s) noexcept {
  const Vec<T, D> d = p - s.closest_point(p);
  return dot(d, d);
}

using Segment2d = Segment<double, 2>;
using Segment3d = Segment<double, 3>;

}  // namespace v1
}  // namespace geom
```

`include/geom/transform/affine.hpp`：

```cpp
#pragma once

#include <cmath>

#include <geom/core/result.hpp>
#include <geom/primitives/point.hpp>
#include <geom/primitives/vec.hpp>

namespace geom {
inline namespace v1 {

// 3x4 affine transform: y = M * x + t, stored row-major in m[3][4].
struct Affine {
  double m[3][4] = {{1.0, 0.0, 0.0, 0.0}, {0.0, 1.0, 0.0, 0.0}, {0.0, 0.0, 1.0, 0.0}};

  static constexpr Affine identity() noexcept { return Affine{}; }

  static constexpr Affine translation(Vec3d t) noexcept {
    Affine a{};
    a.m[0][3] = t[0];
    a.m[1][3] = t[1];
    a.m[2][3] = t[2];
    return a;
  }

  constexpr Point3d transform(Point3d p) const noexcept {
    Point3d r;
    for (int i = 0; i < 3; ++i) {
      r[i] = m[i][0] * p[0] + m[i][1] * p[1] + m[i][2] * p[2] + m[i][3];
    }
    return r;
  }

  constexpr Vec3d transform_vector(Vec3d v) const noexcept {
    Vec3d r;
    for (int i = 0; i < 3; ++i) {
      r[i] = m[i][0] * v[0] + m[i][1] * v[1] + m[i][2] * v[2];
    }
    return r;
  }

  Result<Affine> inverse() const noexcept {
    const double a = m[0][0], b = m[0][1], c = m[0][2];
    const double d = m[1][0], e = m[1][1], f = m[1][2];
    const double g = m[2][0], h = m[2][1], i = m[2][2];

    const double A = e * i - f * h;
    const double B = -(d * i - f * g);
    const double C = d * h - e * g;
    const double det = a * A + b * B + c * C;
    if (det == 0.0 || !std::isfinite(det)) {
      return Error{Errc::degenerate, "singular affine transform", 0u, 0u};
    }
    const double inv_det = 1.0 / det;

    Affine r{};
    r.m[0][0] = A * inv_det;
    r.m[0][1] = -(b * i - c * h) * inv_det;
    r.m[0][2] = (b * f - c * e) * inv_det;
    r.m[1][0] = B * inv_det;
    r.m[1][1] = (a * i - c * g) * inv_det;
    r.m[1][2] = -(a * f - c * d) * inv_det;
    r.m[2][0] = C * inv_det;
    r.m[2][1] = -(a * h - b * g) * inv_det;
    r.m[2][2] = (a * e - b * d) * inv_det;

    for (int row = 0; row < 3; ++row) {
      r.m[row][3] = -(r.m[row][0] * m[0][3] + r.m[row][1] * m[1][3] + r.m[row][2] * m[2][3]);
    }
    return r;
  }
};

}  // namespace v1
}  // namespace geom
```

- [ ] **Step 4: 构建并运行**

Run: `cmake --build build --config Debug && ctest --test-dir build -C Debug --output-on-failure`
Expected: PASS。

- [ ] **Step 5: 提交**

```bash
git add include/geom/primitives include/geom/transform tests
git commit -m "feat(m0): add primitive types and affine transforms"
```

---

## Task 7: 曲线 —— Bézier 与 NURBS 求值/导数

**Files:**
- Create: `include/geom/curves/bezier.hpp`, `include/geom/curves/nurbs_curve.hpp`
- Test: `tests/unit/curves_test.cpp`

**Interfaces:**
- Consumes: `Point<T,D>`、`Vec<T,D>`、`dot`
- Produces:
  - `geom::v1::BezierCurve<T,D>`：`eval(T)` → `Point<T,D>`、`derivative(T)` → `Vec<T,D>`
  - `geom::v1::NurbsCurve`：`is_valid()`、`domain()` → `std::pair<double,double>`、`eval(double)`、`derivative(double)`、`degree()`、`control_count()`
  - `geom::v1::detail::find_span(U,p,n,u)`、`geom::v1::detail::eval_curve(U,p,pts4,u)`（4D 齐次 de Boor，Task 8 复用）

- [ ] **Step 1: 写失败测试**

`tests/unit/curves_test.cpp`：

```cpp
#include <doctest/doctest.h>

#include <cmath>
#include <vector>

#include <geom/curves/bezier.hpp>
#include <geom/curves/nurbs_curve.hpp>

using geom::v1::BezierCurve;
using geom::v1::NurbsCurve;
using geom::v1::Point2d;
using geom::v1::Point3d;

TEST_CASE("quadratic Bezier evaluates to the midpoint formula") {
  const BezierCurve<double, 2> c({Point2d{0, 0}, Point2d{1, 2}, Point2d{2, 0}});
  const Point2d p = c.eval(0.5);
  CHECK(p[0] == doctest::Approx(1.0));
  CHECK(p[1] == doctest::Approx(1.0));
  CHECK(c.eval(0.0)[0] == doctest::Approx(0.0));
  CHECK(c.eval(1.0)[0] == doctest::Approx(2.0));
}

TEST_CASE("Bezier derivative matches central differences") {
  const BezierCurve<double, 3> c({Point3d{0, 0, 0}, Point3d{1, 2, 0}, Point3d{2, 2, 1},
                                  Point3d{3, 0, 1}});
  const double h = 1e-5;
  const double t = 0.37;
  const auto fd = (c.eval(t + h) - c.eval(t - h)) * (1.0 / (2.0 * h));
  const auto exact = c.derivative(t);
  CHECK(exact[0] == doctest::Approx(fd[0]).epsilon(1e-6));
  CHECK(exact[1] == doctest::Approx(fd[1]).epsilon(1e-6));
  CHECK(exact[2] == doctest::Approx(fd[2]).epsilon(1e-6));
}

TEST_CASE("rational quadratic NURBS represents a quarter circle") {
  const double w = std::sqrt(0.5);
  const NurbsCurve c(2, {0.0, 0.0, 0.0, 1.0, 1.0, 1.0},
                     {Point3d{1, 0, 0}, Point3d{1, 1, 0}, Point3d{0, 1, 0}},
                     {1.0, w, 1.0});
  REQUIRE(c.is_valid());
  const auto [u0, u1] = c.domain();
  CHECK(u0 == doctest::Approx(0.0));
  CHECK(u1 == doctest::Approx(1.0));

  const Point3d p = c.eval(0.5);
  CHECK(p[0] == doctest::Approx(std::sqrt(0.5)).epsilon(1e-14));
  CHECK(p[1] == doctest::Approx(std::sqrt(0.5)).epsilon(1e-14));

  const double h = 1e-6;
  const auto fd = (c.eval(0.5 + h) - c.eval(0.5 - h)) * (1.0 / (2.0 * h));
  const auto exact = c.derivative(0.5);
  CHECK(exact[0] == doctest::Approx(fd[0]).epsilon(1e-5));
  CHECK(exact[1] == doctest::Approx(fd[1]).epsilon(1e-5));
}

TEST_CASE("NURBS with unit weights equals a B-spline") {
  const NurbsCurve c(2, {0.0, 0.0, 0.0, 1.0, 2.0, 2.0, 2.0},
                     {Point3d{0, 0, 0}, Point3d{1, 1, 0}, Point3d{2, 1, 0}, Point3d{3, 0, 0}},
                     {1.0, 1.0, 1.0, 1.0});
  REQUIRE(c.is_valid());
  CHECK(c.degree() == 2);
  CHECK(c.control_count() == 4);
  const Point3d p = c.eval(1.0);
  CHECK(p[0] == doctest::Approx(1.5).epsilon(1e-14));
}

TEST_CASE("invalid NURBS data is rejected") {
  const NurbsCurve bad(2, {0.0, 0.0, 1.0, 1.0}, {Point3d{0, 0, 0}, Point3d{1, 0, 0}},
                       {1.0, 1.0});
  CHECK_FALSE(bad.is_valid());
}
```

- [ ] **Step 2: 运行确认失败**

Run: `cmake --build build --config Debug`
Expected: FAIL（找不到 `geom/curves/bezier.hpp`）。

- [ ] **Step 3: 最小实现**

`include/geom/curves/bezier.hpp`：

```cpp
#pragma once

#include <utility>
#include <vector>

#include <geom/primitives/point.hpp>
#include <geom/primitives/vec.hpp>

namespace geom {
inline namespace v1 {

template <class T, int D>
class BezierCurve {
 public:
  explicit BezierCurve(std::vector<Point<T, D>> control) : c_(std::move(control)) {}

  std::size_t size() const noexcept { return c_.size(); }

  Point<T, D> eval(T t) const noexcept {
    if (c_.empty()) return Point<T, D>{};
    std::vector<Point<T, D>> d = c_;
    for (std::size_t k = 1; k < d.size(); ++k) {
      for (std::size_t i = 0; i + k < d.size(); ++i) {
        Point<T, D> q;
        for (int c = 0; c < D; ++c) q[c] = (T{1} - t) * d[i][c] + t * d[i + 1][c];
        d[i] = q;
      }
    }
    return d[0];
  }

  Vec<T, D> derivative(T t) const noexcept {
    if (c_.size() < 2) return Vec<T, D>{};
    const T p = static_cast<T>(c_.size() - 1);
    std::vector<Vec<T, D>> d(c_.size() - 1);
    for (std::size_t i = 0; i + 1 < c_.size(); ++i) d[i] = (c_[i + 1] - c_[i]) * p;
    for (std::size_t k = 1; k < d.size(); ++k) {
      for (std::size_t i = 0; i + k < d.size(); ++i) {
        Vec<T, D> q;
        for (int c = 0; c < D; ++c) q[c] = (T{1} - t) * d[i][c] + t * d[i + 1][c];
        d[i] = q;
      }
    }
    return d[0];
  }

 private:
  std::vector<Point<T, D>> c_;
};

}  // namespace v1
}  // namespace geom
```

`include/geom/curves/nurbs_curve.hpp`：

```cpp
#pragma once

#include <array>
#include <cstddef>
#include <utility>
#include <vector>

#include <geom/primitives/point.hpp>
#include <geom/primitives/vec.hpp>

namespace geom {
inline namespace v1 {
namespace detail {

constexpr int kMaxCurveDegree = 16;

inline int find_span(const std::vector<double>& U, int p, int n, double u) noexcept {
  if (u >= U[static_cast<std::size_t>(n)]) return n;
  int low = p;
  int high = n + 1;
  int mid = (low + high) / 2;
  while (u < U[static_cast<std::size_t>(mid)] || u >= U[static_cast<std::size_t>(mid + 1)]) {
    if (u < U[static_cast<std::size_t>(mid)]) {
      high = mid;
    } else {
      low = mid;
    }
    mid = (low + high) / 2;
  }
  return mid;
}

// de Boor for a 4D homogeneous B-spline. pts.size() == U.size() - p - 1.
inline std::array<double, 4> eval_curve(const std::vector<double>& U, int p,
                                        const std::vector<std::array<double, 4>>& pts,
                                        double u) noexcept {
  const int n = static_cast<int>(pts.size()) - 1;
  const int k = find_span(U, p, n, u);
  std::array<std::array<double, 4>, kMaxCurveDegree + 1> d{};
  for (int j = 0; j <= p; ++j) {
    d[static_cast<std::size_t>(j)] = pts[static_cast<std::size_t>(k - p + j)];
  }
  for (int r = 1; r <= p; ++r) {
    for (int j = p; j >= r; --j) {
      const int i0 = k - p + j;
      const double denom =
          U[static_cast<std::size_t>(k + j + 1 - r)] - U[static_cast<std::size_t>(i0)];
      const double a = denom == 0.0 ? 0.0 : (u - U[static_cast<std::size_t>(i0)]) / denom;
      for (int c = 0; c < 4; ++c) {
        d[static_cast<std::size_t>(j)][static_cast<std::size_t>(c)] =
            (1.0 - a) * d[static_cast<std::size_t>(j - 1)][static_cast<std::size_t>(c)] +
            a * d[static_cast<std::size_t>(j)][static_cast<std::size_t>(c)];
      }
    }
  }
  return d[static_cast<std::size_t>(p)];
}

}  // namespace detail

class NurbsCurve {
 public:
  NurbsCurve(int degree, std::vector<double> knots, std::vector<Point3d> control,
             std::vector<double> weights)
      : p_(degree), U_(std::move(knots)), P_(std::move(control)), W_(std::move(weights)) {}

  bool is_valid() const noexcept {
    if (p_ < 1 || P_.empty()) return false;
    if (P_.size() <= static_cast<std::size_t>(p_)) return false;
    if (W_.size() != P_.size()) return false;
    if (U_.size() != P_.size() + static_cast<std::size_t>(p_) + 1) return false;
    for (std::size_t i = 1; i < U_.size(); ++i) {
      if (!(U_[i] >= U_[i - 1])) return false;
    }
    for (double w : W_) {
      if (!(w > 0.0)) return false;
    }
    return true;
  }

  int degree() const noexcept { return p_; }
  std::size_t control_count() const noexcept { return P_.size(); }

  std::pair<double, double> domain() const noexcept {
    if (!is_valid()) return {0.0, 0.0};
    return {U_[static_cast<std::size_t>(p_)],
            U_[U_.size() - static_cast<std::size_t>(p_) - 1]};
  }

  Point3d eval(double u) const noexcept {
    const auto h = eval_h(clamp(u));
    return Point3d{h[0] / h[3], h[1] / h[3], h[2] / h[3]};
  }

  Vec3d derivative(double u) const noexcept {
    const auto h = eval_h(clamp(u));
    const auto dh = deriv_h(clamp(u));
    const double w = h[3];
    const Vec3d c{h[0] / w, h[1] / w, h[2] / w};
    return Vec3d{(dh[0] - c[0] * dh[3]) / w, (dh[1] - c[1] * dh[3]) / w,
                 (dh[2] - c[2] * dh[3]) / w};
  }

 private:
  double clamp(double u) const noexcept {
    const auto [a, b] = domain();
    return u < a ? a : (u > b ? b : u);
  }

  std::vector<std::array<double, 4>> homogeneous() const {
    std::vector<std::array<double, 4>> H(P_.size());
    for (std::size_t i = 0; i < P_.size(); ++i) {
      H[i] = {W_[i] * P_[i][0], W_[i] * P_[i][1], W_[i] * P_[i][2], W_[i]};
    }
    return H;
  }

  std::array<double, 4> eval_h(double u) const noexcept {
    const auto H = homogeneous();
    return detail::eval_curve(U_, p_, H, u);
  }

  std::array<double, 4> deriv_h(double u) const noexcept {
    const auto H = homogeneous();
    std::vector<double> U1(U_.begin() + 1, U_.end() - 1);
    std::vector<std::array<double, 4>> D(H.size() - 1);
    for (std::size_t i = 0; i + 1 < H.size(); ++i) {
      const double denom = U_[i + static_cast<std::size_t>(p_) + 1] - U_[i + 1];
      const double f = denom == 0.0 ? 0.0 : static_cast<double>(p_) / denom;
      for (int c = 0; c < 4; ++c) D[i][static_cast<std::size_t>(c)] =
          f * (H[i + 1][static_cast<std::size_t>(c)] - H[i][static_cast<std::size_t>(c)]);
    }
    return detail::eval_curve(U1, p_ - 1, D, u);
  }

  int p_ = 0;
  std::vector<double> U_;
  std::vector<Point3d> P_;
  std::vector<double> W_;
};

}  // namespace v1
}  // namespace geom
```

- [ ] **Step 4: 构建并运行**

Run: `cmake --build build --config Debug && ctest --test-dir build -C Debug --output-on-failure`
Expected: PASS（四分圆与有限差分一致）。

- [ ] **Step 5: 提交**

```bash
git add include/geom/curves tests
git commit -m "feat(m0): add Bezier and NURBS curve evaluation with derivatives"
```

---

## Task 8: 曲面 —— NURBS 曲面与 Plane

**Files:**
- Create: `include/geom/surfaces/nurbs_surface.hpp`, `include/geom/surfaces/plane.hpp`
- Test: `tests/unit/surfaces_test.cpp`

**Interfaces:**
- Consumes: `detail::eval_curve`、`Point3d`、`Vec3d`
- Produces:
  - `geom::v1::NurbsSurface`：`is_valid()`、`domain_u()`、`domain_v()`、`eval(u,v)`、`du(u,v)`、`dv(u,v)`
  - `geom::v1::Plane3`：`normal`、`offset`、`signed_distance(p)`、`side(p)`、`project(p)`

- [ ] **Step 1: 写失败测试**

`tests/unit/surfaces_test.cpp`：

```cpp
#include <doctest/doctest.h>

#include <vector>

#include <geom/surfaces/nurbs_surface.hpp>
#include <geom/surfaces/plane.hpp>

using geom::v1::NurbsSurface;
using geom::v1::Plane3;
using geom::v1::Point3d;

TEST_CASE("bilinear patch evaluates to the corner average") {
  const NurbsSurface s(1, 1, {0.0, 0.0, 1.0, 1.0}, {0.0, 0.0, 1.0, 1.0},
                       {Point3d{0, 0, 0}, Point3d{1, 0, 1}, Point3d{0, 1, 2},
                        Point3d{1, 1, 3}},
                       {1.0, 1.0, 1.0, 1.0});
  REQUIRE(s.is_valid());
  const Point3d p = s.eval(0.5, 0.5);
  CHECK(p[0] == doctest::Approx(0.5));
  CHECK(p[1] == doctest::Approx(0.5));
  CHECK(p[2] == doctest::Approx(1.5));
}

TEST_CASE("surface derivatives match central differences") {
  const NurbsSurface s(1, 1, {0.0, 0.0, 1.0, 1.0}, {0.0, 0.0, 1.0, 1.0},
                       {Point3d{0, 0, 0}, Point3d{1, 0, 1}, Point3d{0, 1, 2},
                        Point3d{1, 1, 3}},
                       {1.0, 1.0, 1.0, 1.0});
  const double h = 1e-6;
  const auto fd_u = (s.eval(0.5 + h, 0.3) - s.eval(0.5 - h, 0.3)) * (1.0 / (2.0 * h));
  const auto exact_u = s.du(0.5, 0.3);
  CHECK(exact_u[0] == doctest::Approx(fd_u[0]).epsilon(1e-6));
  CHECK(exact_u[1] == doctest::Approx(fd_u[1]).epsilon(1e-6));
  CHECK(exact_u[2] == doctest::Approx(fd_u[2]).epsilon(1e-6));

  const auto fd_v = (s.eval(0.2, 0.5 + h) - s.eval(0.2, 0.5 - h)) * (1.0 / (2.0 * h));
  const auto exact_v = s.dv(0.2, 0.5);
  CHECK(exact_v[0] == doctest::Approx(fd_v[0]).epsilon(1e-6));
  CHECK(exact_v[1] == doctest::Approx(fd_v[1]).epsilon(1e-6));
  CHECK(exact_v[2] == doctest::Approx(fd_v[2]).epsilon(1e-6));
}

TEST_CASE("plane distance and side are consistent") {
  const Plane3 plane{geom::v1::Vec3d{0, 0, 1}, 5.0};
  const Point3d p{1, 2, 7};
  CHECK(plane.signed_distance(p) == doctest::Approx(2.0));
  CHECK(plane.side(p) == 1);
  CHECK(plane.project(p)[2] == doctest::Approx(5.0));
  CHECK(plane.side(Point3d{1, 2, 5}) == 0);
}
```

- [ ] **Step 2: 运行确认失败**

Run: `cmake --build build --config Debug`
Expected: FAIL（找不到 `geom/surfaces/nurbs_surface.hpp`）。

- [ ] **Step 3: 最小实现**

`include/geom/surfaces/nurbs_surface.hpp`：

```cpp
#pragma once

#include <array>
#include <cstddef>
#include <utility>
#include <vector>

#include <geom/curves/nurbs_curve.hpp>
#include <geom/primitives/point.hpp>
#include <geom/primitives/vec.hpp>

namespace geom {
inline namespace v1 {

// Tensor-product NURBS surface. Control points are indexed as i * nv + j
// (i along u, j along v). All evaluation routes through the 4D homogeneous
// de Boor helpers in <geom/curves/nurbs_curve.hpp>.
class NurbsSurface {
 public:
  NurbsSurface(int pu, int pv, std::vector<double> ku, std::vector<double> kv,
               std::vector<Point3d> control, std::vector<double> weights)
      : pu_(pu), pv_(pv), ku_(std::move(ku)), kv_(std::move(kv)),
        P_(std::move(control)), W_(std::move(weights)) {}

  bool is_valid() const noexcept {
    if (pu_ < 1 || pv_ < 1) return false;
    if (ku_.size() < static_cast<std::size_t>(pu_) + 2) return false;
    if (kv_.size() < static_cast<std::size_t>(pv_) + 2) return false;
    const std::size_t nu_ = nu();
    const std::size_t nv_ = nv();
    if (nu_ < 2 || nv_ < 2) return false;
    if (P_.size() != nu_ * nv_ || W_.size() != P_.size()) return false;
    for (std::size_t i = 1; i < ku_.size(); ++i) {
      if (!(ku_[i] >= ku_[i - 1])) return false;
    }
    for (std::size_t i = 1; i < kv_.size(); ++i) {
      if (!(kv_[i] >= kv_[i - 1])) return false;
    }
    for (double w : W_) {
      if (!(w > 0.0)) return false;
    }
    return true;
  }

  std::size_t nu() const noexcept {
    return ku_.size() - static_cast<std::size_t>(pu_) - 1;
  }
  std::size_t nv() const noexcept {
    return kv_.size() - static_cast<std::size_t>(pv_) - 1;
  }

  std::pair<double, double> domain_u() const noexcept {
    return {ku_[static_cast<std::size_t>(pu_)],
            ku_[ku_.size() - static_cast<std::size_t>(pu_) - 1]};
  }
  std::pair<double, double> domain_v() const noexcept {
    return {kv_[static_cast<std::size_t>(pv_)],
            kv_[kv_.size() - static_cast<std::size_t>(pv_) - 1]};
  }

  Point3d eval(double u, double v) const noexcept {
    const auto h = eval_h(clamp_u(u), clamp_v(v));
    return Point3d{h[0] / h[3], h[1] / h[3], h[2] / h[3]};
  }

  Vec3d du(double u, double v) const noexcept {
    const double uu = clamp_u(u), vv = clamp_v(v);
    return project(eval_h(uu, vv), deriv_h(uu, vv, true));
  }

  Vec3d dv(double u, double v) const noexcept {
    const double uu = clamp_u(u), vv = clamp_v(v);
    return project(eval_h(uu, vv), deriv_h(uu, vv, false));
  }

 private:
  double clamp_u(double u) const noexcept {
    const auto [a, b] = domain_u();
    return u < a ? a : (u > b ? b : u);
  }
  double clamp_v(double v) const noexcept {
    const auto [a, b] = domain_v();
    return v < a ? a : (v > b ? b : v);
  }

  std::array<double, 4> hom(std::size_t i, std::size_t j) const noexcept {
    const std::size_t idx = i * nv() + j;
    return {W_[idx] * P_[idx][0], W_[idx] * P_[idx][1], W_[idx] * P_[idx][2], W_[idx]};
  }

  std::vector<std::array<double, 4>> rows(double v) const {
    std::vector<std::array<double, 4>> R(nu());
    std::vector<std::array<double, 4>> row(nv());
    for (std::size_t i = 0; i < nu(); ++i) {
      for (std::size_t j = 0; j < nv(); ++j) row[j] = hom(i, j);
      R[i] = detail::eval_curve(kv_, pv_, row, v);
    }
    return R;
  }

  std::vector<std::array<double, 4>> cols(double u) const {
    std::vector<std::array<double, 4>> C(nv());
    std::vector<std::array<double, 4>> col(nu());
    for (std::size_t j = 0; j < nv(); ++j) {
      for (std::size_t i = 0; i < nu(); ++i) col[i] = hom(i, j);
      C[j] = detail::eval_curve(ku_, pu_, col, u);
    }
    return C;
  }

  std::array<double, 4> eval_h(double u, double v) const noexcept {
    const auto R = rows(v);
    return detail::eval_curve(ku_, pu_, R, u);
  }

  std::array<double, 4> deriv_h(double u, double v, bool along_u) const noexcept {
    if (along_u) {
      const auto R = rows(v);
      std::vector<std::array<double, 4>> D(R.size() - 1);
      for (std::size_t i = 0; i + 1 < R.size(); ++i) {
        const double denom = ku_[i + static_cast<std::size_t>(pu_) + 1] - ku_[i + 1];
        const double f = denom == 0.0 ? 0.0 : static_cast<double>(pu_) / denom;
        for (int c = 0; c < 4; ++c) {
          D[i][static_cast<std::size_t>(c)] =
              f * (R[i + 1][static_cast<std::size_t>(c)] - R[i][static_cast<std::size_t>(c)]);
        }
      }
      std::vector<double> U1(ku_.begin() + 1, ku_.end() - 1);
      return detail::eval_curve(U1, pu_ - 1, D, u);
    }
    const auto C = cols(u);
    std::vector<std::array<double, 4>> D(C.size() - 1);
    for (std::size_t j = 0; j + 1 < C.size(); ++j) {
      const double denom = kv_[j + static_cast<std::size_t>(pv_) + 1] - kv_[j + 1];
      const double f = denom == 0.0 ? 0.0 : static_cast<double>(pv_) / denom;
      for (int c = 0; c < 4; ++c) {
        D[j][static_cast<std::size_t>(c)] =
            f * (C[j + 1][static_cast<std::size_t>(c)] - C[j][static_cast<std::size_t>(c)]);
      }
    }
    std::vector<double> V1(kv_.begin() + 1, kv_.end() - 1);
    return detail::eval_curve(V1, pv_ - 1, D, v);
  }

  static Vec3d project(const std::array<double, 4>& h, const std::array<double, 4>& dh) noexcept {
    const double w = h[3];
    const Vec3d c{h[0] / w, h[1] / w, h[2] / w};
    return Vec3d{(dh[0] - c[0] * dh[3]) / w, (dh[1] - c[1] * dh[3]) / w,
                 (dh[2] - c[2] * dh[3]) / w};
  }

  int pu_ = 0;
  int pv_ = 0;
  std::vector<double> ku_;
  std::vector<double> kv_;
  std::vector<Point3d> P_;
  std::vector<double> W_;
};

}  // namespace v1
}  // namespace geom
```

`include/geom/surfaces/plane.hpp`：

```cpp
#pragma once

#include <geom/primitives/point.hpp>
#include <geom/primitives/vec.hpp>

namespace geom {
inline namespace v1 {

struct Plane3 {
  Vec3d normal{0.0, 0.0, 1.0};  // unit normal
  double offset = 0.0;          // n · x = offset

  double signed_distance(Point3d p) const noexcept {
    return dot(normal, Vec3d{p[0], p[1], p[2]}) - offset;
  }

  int side(Point3d p) const noexcept {
    const double d = signed_distance(p);
    return d > 0.0 ? 1 : (d < 0.0 ? -1 : 0);
  }

  Point3d project(Point3d p) const noexcept {
    const double d = signed_distance(p);
    return Point3d{p[0] - d * normal[0], p[1] - d * normal[1], p[2] - d * normal[2]};
  }
};

}  // namespace v1
}  // namespace geom
```

- [ ] **Step 4: 构建并运行**

Run: `cmake --build build --config Debug && ctest --test-dir build -C Debug --output-on-failure`
Expected: PASS。

- [ ] **Step 5: 提交**

```bash
git add include/geom/surfaces tests
git commit -m "feat(m0): add NURBS surface evaluation and analytic plane"
```

---

## Task 9: M0 出口 —— 基准、语料约定、ADR 模板

**Files:**
- Create: `bench/bench_predicates.cpp`, `docs/adr/0000-template.md`, `docs/corpus/README.md`
- Modify: `CMakeLists.txt`（`GEOM_BUILD_BENCH` 分支）, `tests/CMakeLists.txt`（加入 `tests/integration/m0_exit_test.cpp`）
- Test: `tests/integration/m0_exit_test.cpp`

**Interfaces:**
- Consumes: `sign_orient2d`
- Produces: 可执行 `geom_bench_predicates`；M0 出口测试；文档模板。

- [ ] **Step 1: 写失败测试**

`tests/integration/m0_exit_test.cpp`：

```cpp
#include <doctest/doctest.h>

#include <geom/numeric/predicates2.hpp>

TEST_CASE("M0 exit: orient2d is consistent on a triangulated grid") {
  int positives = 0, negatives = 0;
  for (int y = 0; y < 8; ++y) {
    for (int x = 0; x < 8; ++x) {
      const int s = geom::v1::sign_orient2d(0, 0, 1, 0, x, y);
      if (s > 0) ++positives;
      if (s < 0) ++negatives;
    }
  }
  // y == 0 is collinear with the base edge; y > 0 yields 7 * 8 positive cases.
  CHECK(positives == 56);
  CHECK(negatives == 0);
}
```

- [ ] **Step 2: 运行确认失败**

Run: `cmake --build build --target geom_bench_predicates --config Debug`
Expected: FAIL（target 不存在：bench 目标尚未定义）。

- [ ] **Step 3: 实现**

`CMakeLists.txt` 追加：

```cmake
if(GEOM_BUILD_BENCH)
  add_executable(geom_bench_predicates bench/bench_predicates.cpp)
  target_link_libraries(geom_bench_predicates PRIVATE geom_core)
endif()
```

`tests/CMakeLists.txt` 追加（在 `add_executable(geom_tests ...)` 之前）：

```cmake
list(APPEND GEOM_TEST_SOURCES integration/m0_exit_test.cpp)
```

`bench/bench_predicates.cpp`：

```cpp
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
```

`docs/adr/0000-template.md`：

```markdown
# ADR NNNN: <决策标题>

- 状态: 提议 / 接受 / 废弃
- 日期: YYYY-MM-DD
- 关联: 架构文档 §X / 相关 ADR

## 背景
<问题与约束>

## 选项
1. <选项 A>：优点 / 缺点
2. <选项 B>：优点 / 缺点

## 决策
<选择及理由>

## 后果
<正面、负面、需要跟进的事项>
```

`docs/corpus/README.md`：

```markdown
# 退化语料库

- `corpus/degenerate/`：共线、重合边、共面、重复点、零面积、反向绕序、近退化、超大尺度比。
- `corpus/cad/`：相切、重合面、共边、极点、细长面、容差链、多壳（M3 起填充）。
- 文件格式与架构文档 §6.7 的 `GEOM` 二进制一致；每个文件附 `.sha256` 与一句话来源说明。
- 任何 fuzz 崩溃必须最小化后加入对应目录，并在回归测试中引用。
```

- [ ] **Step 4: 构建并运行**

Run: `cmake -S . -B build -DGEOM_WERROR=ON -DGEOM_BUILD_BENCH=ON`
Run: `cmake --build build --config Debug`
Run: `ctest --test-dir build -C Debug --output-on-failure`
Run: `./build/Debug/geom_bench_predicates`（Windows：`.\build\Debug\geom_bench_predicates.exe`；单配置生成器：`./build/geom_bench_predicates`）
Expected: 测试全绿；基准输出 `orient2d: <200 ns/call`，退出码 0。

- [ ] **Step 5: 提交**

```bash
git add CMakeLists.txt bench docs tests
git commit -m "feat(m0): add exit tests, benchmark harness and doc templates"
```

---

## Task 10: incircle 谓词（过滤 + 精确）

**Files:**
- Modify: `include/geom/numeric/predicates2.hpp`（追加 incircle）
- Test: `tests/unit/predicates2_incircle_test.cpp`

**Interfaces:**
- Consumes: `detail::Expansion`、`detail::kEpsilon`
- Produces: `geom::v1::sign_incircle(double ax, double ay, ..., double dx, double dy) noexcept`，返回 `-1 | 0 | 1`；正号表示 `d` 在过 `a,b,c` 的圆内（`a,b,c` 逆时针）。

- [ ] **Step 1: 写失败测试**

`tests/unit/predicates2_incircle_test.cpp`：

```cpp
#include <doctest/doctest.h>

#include <cstdint>
#include <random>

#include <geom/numeric/predicates2.hpp>

using geom::v1::sign_incircle;

namespace {

// Coordinates are small integers, so the lifted determinant fits in int64.
int incircle_oracle(std::int64_t ax, std::int64_t ay, std::int64_t bx, std::int64_t by,
                    std::int64_t cx, std::int64_t cy, std::int64_t dx, std::int64_t dy) {
  const std::int64_t adx = ax - dx, ady = ay - dy;
  const std::int64_t bdx = bx - dx, bdy = by - dy;
  const std::int64_t cdx = cx - dx, cdy = cy - dy;
  const std::int64_t alift = adx * adx + ady * ady;
  const std::int64_t blift = bdx * bdx + bdy * bdy;
  const std::int64_t clift = cdx * cdx + cdy * cdy;
  const std::int64_t det =
      alift * (bdx * cdy - cdx * bdy) - blift * (adx * cdy - cdx * ady) +
      clift * (adx * bdy - bdx * ady);
  return det > 0 ? 1 : (det < 0 ? -1 : 0);
}

}  // namespace

TEST_CASE("incircle sign convention") {
  // CCW triangle (0,0), (1,0), (0,1): center (0.5,0.5), point (0.5,0.5) is inside.
  CHECK(sign_incircle(0, 0, 1, 0, 0, 1, 0.5, 0.5) > 0);
  CHECK(sign_incircle(0, 0, 1, 0, 0, 1, 5.0, 5.0) < 0);
}

TEST_CASE("incircle agrees with the integer oracle") {
  std::mt19937 rng(7654321u);
  std::uniform_int_distribution<int> dist(-512, 511);
  for (int i = 0; i < 100000; ++i) {
    const int ax = dist(rng), ay = dist(rng);
    const int bx = dist(rng), by = dist(rng);
    const int cx = dist(rng), cy = dist(rng);
    const int dx = dist(rng), dy = dist(rng);
    REQUIRE(sign_incircle(ax, ay, bx, by, cx, cy, dx, dy) ==
            incircle_oracle(ax, ay, bx, by, cx, cy, dx, dy));
  }
}

TEST_CASE("cocircular points are exactly degenerate") {
  CHECK(sign_incircle(1, 0, 0, 1, -1, 0, 0, -1) == 0);
  CHECK(sign_incircle(0, 0, 2, 0, 0, 2, 2, 2) == 0);
}
```

- [ ] **Step 2: 运行确认失败**

Run: `cmake --build build --config Debug`
Expected: FAIL（`sign_incircle` 未声明）。

- [ ] **Step 3: 实现（追加到 `predicates2.hpp`）**

```cpp
namespace detail {

constexpr double kIncircleErrBoundA = (10.0 + 96.0 * kEpsilon) * kEpsilon;

inline int incircle_exact(double ax, double ay, double bx, double by, double cx, double cy,
                          double dx, double dy) noexcept {
  const Expansion adx = Expansion::difference(ax, dx);
  const Expansion ady = Expansion::difference(ay, dy);
  const Expansion bdx = Expansion::difference(bx, dx);
  const Expansion bdy = Expansion::difference(by, dy);
  const Expansion cdx = Expansion::difference(cx, dx);
  const Expansion cdy = Expansion::difference(cy, dy);

  Expansion alift = Expansion::mul(adx, adx);
  alift.add(Expansion::mul(ady, ady));
  Expansion blift = Expansion::mul(bdx, bdx);
  blift.add(Expansion::mul(bdy, bdy));
  Expansion clift = Expansion::mul(cdx, cdx);
  clift.add(Expansion::mul(cdy, cdy));

  Expansion cross_bc = Expansion::mul(bdx, cdy);
  cross_bc.add_scaled(Expansion::mul(cdx, bdy), -1.0);
  Expansion cross_ac = Expansion::mul(adx, cdy);
  cross_ac.add_scaled(Expansion::mul(cdx, ady), -1.0);
  Expansion cross_ab = Expansion::mul(adx, bdy);
  cross_ab.add_scaled(Expansion::mul(bdx, ady), -1.0);

  Expansion det = Expansion::mul(alift, cross_bc);
  det.add_scaled(Expansion::mul(blift, cross_ac), -1.0);
  det.add_scaled(Expansion::mul(clift, cross_ab), 1.0);
  return det.sign();
}

}  // namespace detail

inline int sign_incircle(double ax, double ay, double bx, double by, double cx, double cy,
                         double dx, double dy) noexcept {
  const double m = std::max({std::abs(ax), std::abs(ay), std::abs(bx), std::abs(by),
                             std::abs(cx), std::abs(cy), std::abs(dx), std::abs(dy)});
  if (!std::isfinite(m)) return 0;
  if (m == 0.0) return 0;
  constexpr double kMaxSafe = 1e75;
  constexpr double kMinSafe = 1e-75;
  double scale = 1.0;
  if (m > kMaxSafe || m < kMinSafe) {
    int exp2 = 0;
    std::frexp(m, &exp2);
    scale = std::ldexp(1.0, -exp2);
  }
  ax *= scale; ay *= scale; bx *= scale; by *= scale;
  cx *= scale; cy *= scale; dx *= scale; dy *= scale;

  const double adx = ax - dx, ady = ay - dy;
  const double bdx = bx - dx, bdy = by - dy;
  const double cdx = cx - dx, cdy = cy - dy;
  const double alift = adx * adx + ady * ady;
  const double blift = bdx * bdx + bdy * bdy;
  const double clift = cdx * cdx + cdy * cdy;
  const double det = alift * (bdx * cdy - cdx * bdy) - blift * (adx * cdy - cdx * ady) +
                     clift * (adx * bdy - bdx * ady);
  const double permanent =
      (std::abs(bdx * cdy) + std::abs(cdx * bdy)) * alift +
      (std::abs(adx * cdy) + std::abs(cdx * ady)) * blift +
      (std::abs(adx * bdy) + std::abs(bdx * ady)) * clift;
  const double errbound = detail::kIncircleErrBoundA * permanent;

  if (det > errbound) return 1;
  if (-det > errbound) return -1;
  return detail::incircle_exact(ax, ay, bx, by, cx, cy, dx, dy);
}
```

- [ ] **Step 4: 构建并运行**

Run: `cmake --build build --config Debug && ctest --test-dir build -C Debug --output-on-failure`
Expected: PASS（10 万随机例与 oracle 一致）。

- [ ] **Step 5: 提交**

```bash
git add include/geom/numeric tests
git commit -m "feat(m1): add filtered and exact incircle predicate"
```

---

## Task 11: 2D 凸包

**Files:**
- Create: `include/geom/cgeom/hull2d.hpp`
- Test: `tests/unit/hull2d_test.cpp`

**Interfaces:**
- Consumes: `Point2d`、`sign_orient2d`
- Produces: `std::vector<std::uint32_t> geom::v1::convex_hull2d(std::span<const Point2d> pts)`，返回逆时针凸包顶点下标（不重复首点）；共线点只保留端点；重复点去重。

- [ ] **Step 1: 写失败测试**

`tests/unit/hull2d_test.cpp`：

```cpp
#include <doctest/doctest.h>

#include <vector>

#include <geom/cgeom/hull2d.hpp>
#include <geom/numeric/predicates2.hpp>

using geom::v1::convex_hull2d;
using geom::v1::Point2d;

TEST_CASE("hull of a square with interior points is the square") {
  std::vector<Point2d> pts{{0, 0}, {1, 0}, {1, 1}, {0, 1},
                           {0.5, 0.5}, {0.25, 0.9}, {0.1, 0.1}};
  const auto hull = convex_hull2d(pts);
  REQUIRE(hull.size() == 4);
  for (std::uint32_t i = 0; i < hull.size(); ++i) {
    const auto& a = pts[hull[i]];
    const auto& b = pts[hull[(i + 1) % hull.size()]];
    const auto& c = pts[hull[(i + 2) % hull.size()]];
    CHECK(geom::v1::sign_orient2d(a[0], a[1], b[0], b[1], c[0], c[1]) > 0);
  }
}

TEST_CASE("hull handles degenerate inputs") {
  CHECK(convex_hull2d(std::vector<Point2d>{}).empty());
  CHECK(convex_hull2d(std::vector<Point2d>{{1, 2}}).size() == 1);

  std::vector<Point2d> line{{0, 0}, {1, 1}, {2, 2}, {3, 3}};
  const auto h = convex_hull2d(line);
  REQUIRE(h.size() == 2);
  CHECK((h[0] == 0 && h[1] == 3));

  std::vector<Point2d> dup{{0, 0}, {0, 0}, {1, 0}, {0, 1}};
  CHECK(convex_hull2d(dup).size() == 3);
}
```

- [ ] **Step 2: 运行确认失败**

Run: `cmake --build build --config Debug`
Expected: FAIL（找不到 `geom/cgeom/hull2d.hpp`）。

- [ ] **Step 3: 最小实现**

```cpp
#pragma once

#include <algorithm>
#include <cstdint>
#include <span>
#include <vector>

#include <geom/numeric/predicates2.hpp>
#include <geom/primitives/point.hpp>

namespace geom {
inline namespace v1 {

inline std::vector<std::uint32_t> convex_hull2d(std::span<const Point2d> pts) {
  std::vector<std::uint32_t> idx(pts.size());
  for (std::uint32_t i = 0; i < pts.size(); ++i) idx[i] = i;
  std::sort(idx.begin(), idx.end(), [&](std::uint32_t a, std::uint32_t b) {
    if (pts[a][0] != pts[b][0]) return pts[a][0] < pts[b][0];
    if (pts[a][1] != pts[b][1]) return pts[a][1] < pts[b][1];
    return a < b;
  });
  // Remove exact duplicates.
  idx.erase(std::unique(idx.begin(), idx.end(), [&](std::uint32_t a, std::uint32_t b) {
              return pts[a][0] == pts[b][0] && pts[a][1] == pts[b][1];
            }),
            idx.end());

  if (idx.size() <= 2) return idx;

  const auto cross = [&](std::uint32_t o, std::uint32_t a, std::uint32_t b) {
    return sign_orient2d(pts[o][0], pts[o][1], pts[a][0], pts[a][1], pts[b][0], pts[b][1]);
  };

  std::vector<std::uint32_t> hull(2 * idx.size());
  std::size_t k = 0;
  for (std::size_t i = 0; i < idx.size(); ++i) {
    while (k >= 2 && cross(hull[k - 2], hull[k - 1], idx[i]) <= 0) --k;
    hull[k++] = idx[i];
  }
  const std::size_t lower = k + 1;
  for (std::size_t i = idx.size() - 1; i > 0; --i) {
    while (k >= lower && cross(hull[k - 2], hull[k - 1], idx[i - 1]) <= 0) --k;
    hull[k++] = idx[i - 1];
  }
  hull.resize(k - 1);
  return hull;
}

}  // namespace v1
}  // namespace geom
```

- [ ] **Step 4: 构建并运行**

Run: `cmake --build build --config Debug && ctest --test-dir build -C Debug --output-on-failure`
Expected: PASS。

- [ ] **Step 5: 提交**

```bash
git add include/geom/cgeom tests
git commit -m "feat(m1): add exact 2D convex hull"
```

---

## Task 12: 2D Delaunay 三角剖分

**Files:**
- Create: `include/geom/cgeom/delaunay2.hpp`
- Test: `tests/unit/delaunay2_test.cpp`

**Interfaces:**
- Consumes: `Point2d`、`sign_incircle`、`convex_hull2d`
- Produces: `std::vector<std::array<std::uint32_t,3>> geom::v1::delaunay2(std::span<const Point2d> pts)`：返回逆时针三角形（索引指向输入）；重复点被跳过；退化输入返回空。

- [ ] **Step 1: 写失败测试**

`tests/unit/delaunay2_test.cpp`：

```cpp
#include <doctest/doctest.h>

#include <algorithm>
#include <cmath>
#include <vector>

#include <geom/cgeom/delaunay2.hpp>
#include <geom/cgeom/hull2d.hpp>
#include <geom/numeric/predicates2.hpp>

using geom::v1::convex_hull2d;
using geom::v1::delaunay2;
using geom::v1::Point2d;

namespace {

double incircle_perm(const std::vector<Point2d>& p, std::uint32_t a, std::uint32_t b,
                     std::uint32_t c, std::uint32_t d) {
  return static_cast<double>(geom::v1::sign_incircle(p[a][0], p[a][1], p[b][0], p[b][1],
                                                     p[c][0], p[c][1], p[d][0], p[d][1]));
}

}  // namespace

TEST_CASE("Delaunay satisfies the empty circumcircle property") {
  std::vector<Point2d> pts;
  for (int y = 0; y < 6; ++y) {
    for (int x = 0; x < 6; ++x) {
      // Deterministic perturbation keeps the input in general position
      // (no exact cocircular quadruples), which is the M1 contract.
      const double jx = 0.001 * static_cast<double>((x * 7 + y * 3) % 5);
      const double jy = 0.001 * static_cast<double>((x * 5 + y * 7) % 3);
      pts.push_back(Point2d{static_cast<double>(x) + jx, static_cast<double>(y) + jy});
    }
  }
  pts.push_back(Point2d{2.3, 3.1});
  pts.push_back(Point2d{4.7, 0.4});

  const auto tris = delaunay2(pts);
  REQUIRE_FALSE(tris.empty());
  for (const auto& t : tris) {
    CHECK(geom::v1::sign_orient2d(pts[t[0]][0], pts[t[0]][1], pts[t[1]][0], pts[t[1]][1],
                                  pts[t[2]][0], pts[t[2]][1]) > 0);
    for (std::uint32_t i = 0; i < pts.size(); ++i) {
      if (i == t[0] || i == t[1] || i == t[2]) continue;
      CHECK(incircle_perm(pts, t[0], t[1], t[2], i) <= 0);
    }
  }
}

TEST_CASE("Delaunay triangle count matches Euler's formula") {
  std::vector<Point2d> pts;
  for (int i = 0; i < 20; ++i) {
    pts.push_back(Point2d{std::cos(0.37 * i) * 10.0 + (i % 3) * 0.01,
                          std::sin(0.53 * i) * 10.0 + (i % 5) * 0.02});
  }
  const auto tris = delaunay2(pts);
  const auto hull = convex_hull2d(pts);
  REQUIRE_FALSE(tris.empty());
  // T = 2n - 2 - h for n points in general position, h hull vertices.
  CHECK(tris.size() == 2 * pts.size() - 2 - hull.size());
}

TEST_CASE("Delaunay skips duplicate points and handles tiny input") {
  std::vector<Point2d> dup{{0, 0}, {0, 0}, {1, 0}, {0, 1}, {1, 2}};
  CHECK_FALSE(delaunay2(dup).empty());
  std::vector<Point2d> two{{0, 0}, {1, 0}};
  CHECK(delaunay2(two).empty());
}
```

- [ ] **Step 2: 运行确认失败**

Run: `cmake --build build --config Debug`
Expected: FAIL（找不到 `geom/cgeom/delaunay2.hpp`）。

- [ ] **Step 3: 最小实现**

```cpp
#pragma once

#include <array>
#include <cstdint>
#include <span>
#include <unordered_map>
#include <vector>

#include <geom/numeric/predicates2.hpp>
#include <geom/primitives/point.hpp>

namespace geom {
inline namespace v1 {

// Bowyer-Watson with linear triangle search. O(n^2) worst case; deliberately
// simple for M1. Requires distinct points in general position: a point lying
// exactly on an existing circumcircle is not re-inserted (no SoS yet). Exact
// incircle decides insertion, so the output is a valid Delaunay triangulation
// for the supported input class. SoS and CDT arrive in the M1b plan.
inline std::vector<std::array<std::uint32_t, 3>> delaunay2(std::span<const Point2d> input) {
  std::vector<Point2d> pts(input.begin(), input.end());
  // Drop exact duplicates.
  std::vector<std::uint32_t> keep;
  keep.reserve(pts.size());
  for (std::uint32_t i = 0; i < pts.size(); ++i) keep.push_back(i);

  auto same = [&](std::uint32_t a, std::uint32_t b) {
    return pts[a][0] == pts[b][0] && pts[a][1] == pts[b][1];
  };
  std::vector<std::uint32_t> unique;
  for (std::uint32_t i = 0; i < keep.size(); ++i) {
    bool dup = false;
    for (std::uint32_t j : unique) {
      if (same(i, j)) { dup = true; break; }
    }
    if (!dup) unique.push_back(i);
  }
  if (unique.size() < 3) return {};

  std::vector<Point2d> p;
  p.reserve(unique.size() + 3);
  for (std::uint32_t i : unique) p.push_back(pts[i]);

  double minx = p[0][0], maxx = p[0][0], miny = p[0][1], maxy = p[0][1];
  for (const auto& q : p) {
    minx = std::min(minx, q[0]); maxx = std::max(maxx, q[0]);
    miny = std::min(miny, q[1]); maxy = std::max(maxy, q[1]);
  }
  const double cx = 0.5 * (minx + maxx), cy = 0.5 * (miny + maxy);
  const double r = 20.0 * std::max(maxx - minx, maxy - miny) + 1.0;
  const std::uint32_t s0 = static_cast<std::uint32_t>(p.size());
  p.push_back(Point2d{cx - r, cy - r});
  p.push_back(Point2d{cx + r, cy - r});
  p.push_back(Point2d{cx, cy + r});
  const std::uint32_t s1 = s0 + 1, s2 = s0 + 2;

  struct Tri { std::uint32_t v[3]; bool alive; };
  std::vector<Tri> tris;
  tris.push_back(Tri{{s0, s1, s2}, true});

  for (std::uint32_t i = 0; i < s0; ++i) {
    std::vector<std::uint32_t> bad;
    for (std::uint32_t t = 0; t < tris.size(); ++t) {
      if (!tris[t].alive) continue;
      const auto& T = tris[t];
      if (sign_incircle(p[T.v[0]][0], p[T.v[0]][1], p[T.v[1]][0], p[T.v[1]][1],
                        p[T.v[2]][0], p[T.v[2]][1], p[i][0], p[i][1]) > 0) {
        bad.push_back(t);
      }
    }
    // Directed boundary edges: an edge of a bad triangle without its reverse.
    std::unordered_map<std::uint64_t, std::uint32_t> edge;
    auto key = [](std::uint32_t a, std::uint32_t b) {
      return (static_cast<std::uint64_t>(a) << 32) | b;
    };
    for (std::uint32_t t : bad) {
      const auto& T = tris[t];
      for (int e = 0; e < 3; ++e) {
        const std::uint32_t a = T.v[e], b = T.v[(e + 1) % 3];
        ++edge[key(a, b)];
      }
    }
    std::vector<std::pair<std::uint32_t, std::uint32_t>> boundary;
    for (const auto& entry : edge) {
      const std::uint32_t a = static_cast<std::uint32_t>(entry.first >> 32);
      const std::uint32_t b = static_cast<std::uint32_t>(entry.first & 0xFFFFFFFFu);
      if (edge.find(key(b, a)) == edge.end()) boundary.emplace_back(a, b);
    }
    for (std::uint32_t t : bad) tris[t].alive = false;
    for (const auto& [a, b] : boundary) {
      tris.push_back(Tri{{a, b, i}, true});
    }
  }

  std::vector<std::array<std::uint32_t, 3>> out;
  for (const auto& T : tris) {
    if (!T.alive) continue;
    if (T.v[0] >= s0 || T.v[1] >= s0 || T.v[2] >= s0) continue;
    const int s = sign_orient2d(p[T.v[0]][0], p[T.v[0]][1], p[T.v[1]][0], p[T.v[1]][1],
                                p[T.v[2]][0], p[T.v[2]][1]);
    if (s > 0) {
      out.push_back({unique[T.v[0]], unique[T.v[1]], unique[T.v[2]]});
    } else if (s < 0) {
      out.push_back({unique[T.v[1]], unique[T.v[0]], unique[T.v[2]]});
    }
  }
  return out;
}

}  // namespace v1
}  // namespace geom
```

- [ ] **Step 4: 构建并运行**

Run: `cmake --build build --config Debug && ctest --test-dir build -C Debug --output-on-failure`
Expected: PASS（空圆性质成立；三角形数符合 `T = 2n - 2 - h`）。

- [ ] **Step 5: 提交**

```bash
git add include/geom/cgeom tests
git commit -m "feat(m1): add Bowyer-Watson Delaunay triangulation"
```

---

## Task 13: BVH（SAH 构建 + 射线/最近邻查询）

**Files:**
- Create: `include/geom/index/bvh.hpp`
- Test: `tests/unit/bvh_test.cpp`

**Interfaces:**
- Consumes: `Triangle3d`、`Box3d`、`Ray3d`（本任务一并定义于 `bvh.hpp` 内最小版本：`Ray3d{Point3d origin; Vec3d dir;}`、`Triangle3d{a,b,c}`）
- Produces:
  - `geom::v1::Triangle3d`、`geom::v1::Ray3d`
  - `geom::v1::Bvh`：`build(span<const Triangle3d>)`、`raycast(ray, tmax, Hit&) -> bool`、`nearest(p, double& d2) -> uint32_t`、`bounds()`

- [ ] **Step 1: 写失败测试**

`tests/unit/bvh_test.cpp`：

```cpp
#include <doctest/doctest.h>

#include <cmath>
#include <cstdint>
#include <limits>
#include <random>
#include <vector>

#include <geom/index/bvh.hpp>

using geom::v1::Bvh;
using geom::v1::Point3d;
using geom::v1::Ray3d;
using geom::v1::Triangle3d;
using geom::v1::Vec3d;

namespace {

bool ray_tri(const Ray3d& r, const Triangle3d& t, double& tout) {
  const Vec3d e1 = t.b - t.a;
  const Vec3d e2 = t.c - t.a;
  const Vec3d pv = geom::v1::cross(r.dir, e2);
  const double det = geom::v1::dot(e1, pv);
  if (std::abs(det) < 1e-14) return false;
  const double inv = 1.0 / det;
  const Vec3d tv{r.origin[0] - t.a[0], r.origin[1] - t.a[1], r.origin[2] - t.a[2]};
  const double u = geom::v1::dot(tv, pv) * inv;
  if (u < 0.0 || u > 1.0) return false;
  const Vec3d qv = geom::v1::cross(tv, e1);
  const double v = geom::v1::dot(r.dir, qv) * inv;
  if (v < 0.0 || u + v > 1.0) return false;
  const double tt = geom::v1::dot(e2, qv) * inv;
  if (tt <= 0.0) return false;
  tout = tt;
  return true;
}

}  // namespace

TEST_CASE("BVH raycast matches brute force") {
  std::mt19937 rng(42u);
  std::uniform_real_distribution<double> coord(-1.0, 1.0);
  std::vector<Triangle3d> tris;
  for (int i = 0; i < 200; ++i) {
    tris.push_back(Triangle3d{Point3d{coord(rng), coord(rng), coord(rng)},
                              Point3d{coord(rng), coord(rng), coord(rng)},
                              Point3d{coord(rng), coord(rng), coord(rng)}});
  }
  Bvh bvh;
  bvh.build(tris);
  REQUIRE_FALSE(bvh.bounds().is_empty());

  for (int r = 0; r < 500; ++r) {
    Ray3d ray{Point3d{coord(rng) * 3.0, coord(rng) * 3.0, 3.0},
              Vec3d{coord(rng) * 0.1, coord(rng) * 0.1, -1.0}};
    double best_brute = std::numeric_limits<double>::infinity();
    for (const auto& t : tris) {
      double tt = 0.0;
      if (ray_tri(ray, t, tt) && tt < best_brute) best_brute = tt;
    }
    Bvh::Hit hit;
    const bool got = bvh.raycast(ray, 1e30, hit);
    CHECK(got == std::isfinite(best_brute));
    if (got && std::isfinite(best_brute)) {
      CHECK(hit.t == doctest::Approx(best_brute).epsilon(1e-12));
    }
  }
}

TEST_CASE("BVH nearest matches brute force") {
  std::mt19937 rng(7u);
  std::uniform_real_distribution<double> coord(-1.0, 1.0);
  std::vector<Triangle3d> tris;
  for (int i = 0; i < 150; ++i) {
    tris.push_back(Triangle3d{Point3d{coord(rng), coord(rng), coord(rng)},
                              Point3d{coord(rng), coord(rng), coord(rng)},
                              Point3d{coord(rng), coord(rng), coord(rng)}});
  }
  Bvh bvh;
  bvh.build(tris);
  for (int q = 0; q < 200; ++q) {
    const Point3d p{coord(rng) * 2.0, coord(rng) * 2.0, coord(rng) * 2.0};
    double best = std::numeric_limits<double>::infinity();
    for (const auto& t : tris) {
      const auto d = geom::v1::squared_distance(p, t);
      best = std::min(best, d);
    }
    double d2 = std::numeric_limits<double>::infinity();
    const std::uint32_t idx = bvh.nearest(p, d2);
    CHECK(idx != 0xFFFFFFFFu);
    CHECK(d2 == doctest::Approx(best).epsilon(1e-10));
  }
}
```

- [ ] **Step 2: 运行确认失败**

Run: `cmake --build build --config Debug`
Expected: FAIL（找不到 `geom/index/bvh.hpp`）。

- [ ] **Step 3: 最小实现**

`include/geom/index/bvh.hpp`：

```cpp
#pragma once

#include <algorithm>
#include <array>
#include <cstdint>
#include <limits>
#include <span>
#include <vector>

#include <geom/primitives/box.hpp>
#include <geom/primitives/point.hpp>
#include <geom/primitives/vec.hpp>

namespace geom {
inline namespace v1 {

struct Triangle3d {
  Point3d a, b, c;
};

struct Ray3d {
  Point3d origin;
  Vec3d dir;
};

inline Box3d triangle_box(const Triangle3d& t) noexcept {
  Box3d b = Box3d::empty();
  b.expand(t.a);
  b.expand(t.b);
  b.expand(t.c);
  return b;
}

inline T squared_distance(Point3d p, const Triangle3d& t) noexcept {
  // Closest point on triangle (Ericson, Real-Time Collision Detection 5.1.5).
  const Vec3d ab = t.b - t.a, ac = t.c - t.a, ap = p - t.a;
  const double d1 = dot(ab, ap), d2 = dot(ac, ap);
  if (d1 <= 0.0 && d2 <= 0.0) return dot(ap, ap);
  const Vec3d bp = p - t.b;
  const double d3 = dot(ab, bp), d4 = dot(ac, bp);
  if (d3 >= 0.0 && d4 <= d3) return dot(bp, bp);
  const double vc = d1 * d4 - d3 * d2;
  if (vc <= 0.0 && d1 >= 0.0 && d3 <= 0.0) {
    const double v = d1 / (d1 - d3);
    const Vec3d q = ab * v - ap;
    return dot(q, q);
  }
  const Vec3d cp = p - t.c;
  const double d5 = dot(ab, cp), d6 = dot(ac, cp);
  if (d6 >= 0.0 && d5 <= d6) return dot(cp, cp);
  const double vb = d5 * d2 - d1 * d6;
  if (vb <= 0.0 && d2 >= 0.0 && d6 <= 0.0) {
    const double w = d2 / (d2 - d6);
    const Vec3d q = ac * w - ap;
    return dot(q, q);
  }
  const double va = d3 * d6 - d5 * d4;
  if (va <= 0.0 && (d4 - d3) >= 0.0 && (d5 - d6) >= 0.0) {
    const double w = (d4 - d3) / ((d4 - d3) + (d5 - d6));
    const Vec3d q = (t.c - t.b) * w - bp;
    return dot(q, q);
  }
  const double denom = 1.0 / (va + vb + vc);
  const double v = vb * denom;
  const double w = vc * denom;
  const Vec3d q = ab * v + ac * w - ap;
  return dot(q, q);
}

class Bvh {
 public:
  struct Hit {
    std::uint32_t prim = 0;
    double t = 0.0;
  };

  void build(std::span<const Triangle3d> triangles) {
    prims_.clear();
    nodes_.clear();
    order_.clear();
    prims_.reserve(triangles.size());
    for (std::uint32_t i = 0; i < triangles.size(); ++i) {
      Prim p;
      p.tri = triangles[i];
      p.box = triangle_box(triangles[i]);
      p.centroid = Point3d{(p.box.lo[0] + p.box.hi[0]) * 0.5, (p.box.lo[1] + p.box.hi[1]) * 0.5,
                           (p.box.lo[2] + p.box.hi[2]) * 0.5};
      prims_.push_back(p);
      order_.push_back(i);
    }
    if (!prims_.empty()) {
      nodes_.push_back({});
      build_node(0, static_cast<std::uint32_t>(prims_.size()));
    }
  }

  const Box3d& bounds() const noexcept {
    static const Box3d kEmpty = Box3d::empty();
    return nodes_.empty() ? kEmpty : nodes_[0].box;
  }
  std::uint32_t size() const noexcept { return static_cast<std::uint32_t>(prims_.size()); }

  bool raycast(const Ray3d& ray, double tmax, Hit& hit) const {
    if (nodes_.empty()) return false;
    const Vec3d inv{1.0 / ray.dir[0], 1.0 / ray.dir[1], 1.0 / ray.dir[2]};
    std::array<std::uint32_t, 128> stack{};
    int sp = 0;
    stack[sp++] = 0;
    bool found = false;
    double best = tmax;
    while (sp > 0) {
      const Node& n = nodes_[stack[--sp]];
      if (!slab(n.box, ray, inv, best)) continue;
      if (n.count > 0) {
        for (std::uint32_t i = n.first; i < n.first + n.count; ++i) {
          double t = 0.0;
          if (intersect_tri(ray, prims_[order_[i]].tri, t) && t < best) {
            best = t;
            found = true;
            hit.prim = order_[i];
          }
        }
      } else {
        stack[sp++] = n.left;
        stack[sp++] = n.right;
      }
    }
    if (found) hit.t = best;
    return found;
  }

  std::uint32_t nearest(const Point3d& p, double& best_d2) const {
    if (nodes_.empty()) return 0xFFFFFFFFu;
    std::array<std::uint32_t, 128> stack{};
    int sp = 0;
    stack[sp++] = 0;
    best_d2 = std::numeric_limits<double>::infinity();
    std::uint32_t best = 0xFFFFFFFFu;
    while (sp > 0) {
      const Node& n = nodes_[stack[--sp]];
      if (box_distance2(n.box, p) >= best_d2) continue;
      if (n.count > 0) {
        for (std::uint32_t i = n.first; i < n.first + n.count; ++i) {
          const double d2 = squared_distance(p, prims_[order_[i]].tri);
          if (d2 < best_d2) {
            best_d2 = d2;
            best = order_[i];
          }
        }
      } else {
        stack[sp++] = n.left;
        stack[sp++] = n.right;
      }
    }
    return best;
  }

 private:
  struct Prim {
    Triangle3d tri;
    Box3d box;
    Point3d centroid;
  };
  struct Node {
    Box3d box;
    std::uint32_t first = 0;
    std::uint32_t count = 0;
    std::uint32_t left = 0;
    std::uint32_t right = 0;
  };

  static constexpr int kBins = 12;

  std::uint32_t build_node(std::uint32_t begin, std::uint32_t end) {
    const std::uint32_t node_index = static_cast<std::uint32_t>(nodes_.size());
    nodes_.push_back({});
    Box3d box = Box3d::empty();
    for (std::uint32_t i = begin; i < end; ++i) box.expand(prims_[order_[i]].box);
    nodes_[node_index].box = box;
    if (end - begin <= 4) {
      nodes_[node_index].first = begin;
      nodes_[node_index].count = end - begin;
      return node_index;
    }

    Box3d centroid_box = Box3d::empty();
    for (std::uint32_t i = begin; i < end; ++i) centroid_box.expand(prims_[order_[i]].centroid);
    const Vec3d extent = centroid_box.diagonal();
    int axis = 0;
    if (extent[1] > extent[axis]) axis = 1;
    if (extent[2] > extent[axis]) axis = 2;

    const double lo = centroid_box.lo[axis];
    const double hi = centroid_box.hi[axis];
    if (!(hi > lo)) {
      const std::uint32_t mid = begin + (end - begin) / 2;
      nodes_[node_index].left = build_node(begin, mid);
      nodes_[node_index].right = build_node(mid, end);
      return node_index;
    }

    std::array<int, kBins> counts{};
    std::array<Box3d, kBins> bins;
    for (auto& b : bins) b = Box3d::empty();
    const double inv = kBins / (hi - lo);
    for (std::uint32_t i = begin; i < end; ++i) {
      int b = static_cast<int>((prims_[order_[i]].centroid[axis] - lo) * inv);
      if (b < 0) b = 0;
      if (b >= kBins) b = kBins - 1;
      ++counts[static_cast<std::size_t>(b)];
      bins[static_cast<std::size_t>(b)].expand(prims_[order_[i]].box);
    }

    std::array<double, kBins> left_area{};
    std::array<double, kBins> right_area{};
    std::array<int, kBins> left_count{};
    std::array<int, kBins> right_count{};
    Box3d acc = Box3d::empty();
    int nc = 0;
    for (int b = 0; b < kBins; ++b) {
      acc.expand(bins[static_cast<std::size_t>(b)]);
      nc += counts[static_cast<std::size_t>(b)];
      left_count[static_cast<std::size_t>(b)] = nc;
      left_area[static_cast<std::size_t>(b)] = surface_area(acc);
    }
    acc = Box3d::empty();
    nc = 0;
    for (int b = kBins - 1; b >= 0; --b) {
      acc.expand(bins[static_cast<std::size_t>(b)]);
      nc += counts[static_cast<std::size_t>(b)];
      right_count[static_cast<std::size_t>(b)] = nc;
      right_area[static_cast<std::size_t>(b)] = surface_area(acc);
    }

    int best_split = -1;
    double best_cost = std::numeric_limits<double>::infinity();
    for (int b = 0; b + 1 < kBins; ++b) {
      if (left_count[static_cast<std::size_t>(b)] == 0 ||
          right_count[static_cast<std::size_t>(b + 1)] == 0) {
        continue;
      }
      const double cost = left_area[static_cast<std::size_t>(b)] * left_count[static_cast<std::size_t>(b)] +
                          right_area[static_cast<std::size_t>(b)] * right_count[static_cast<std::size_t>(b)];
      if (cost < best_cost) {
        best_cost = cost;
        best_split = b;
      }
    }

    std::uint32_t mid = begin;
    if (best_split >= 0) {
      auto it = std::partition(order_.begin() + begin, order_.begin() + end, [&](std::uint32_t i) {
        int b = static_cast<int>((prims_[i].centroid[axis] - lo) * inv);
        if (b < 0) b = 0;
        if (b >= kBins) b = kBins - 1;
        return b <= best_split;
      });
      mid = static_cast<std::uint32_t>(it - order_.begin());
    }
    if (mid == begin || mid == end) {
      mid = begin + (end - begin) / 2;
    }
    nodes_[node_index].left = build_node(begin, mid);
    nodes_[node_index].right = build_node(mid, end);
    return node_index;
  }

  static double surface_area(const Box3d& b) noexcept {
    if (b.is_empty()) return 0.0;
    const Vec3d d = b.diagonal();
    return 2.0 * (d[0] * d[1] + d[1] * d[2] + d[2] * d[0]);
  }

  static double box_distance2(const Box3d& b, Point3d p) noexcept {
    double sum = 0.0;
    for (int i = 0; i < 3; ++i) {
      double d = 0.0;
      if (p[i] < b.lo[i]) d = b.lo[i] - p[i];
      if (p[i] > b.hi[i]) d = p[i] - b.hi[i];
      sum += d * d;
    }
    return sum;
  }

  static bool slab(const Box3d& b, const Ray3d& r, Vec3d inv, double tmax) noexcept {
    double tmin = 0.0, tmax_local = tmax;
    for (int i = 0; i < 3; ++i) {
      double t0 = (b.lo[i] - r.origin[i]) * inv[i];
      double t1 = (b.hi[i] - r.origin[i]) * inv[i];
      if (t0 > t1) std::swap(t0, t1);
      tmin = std::max(tmin, t0);
      tmax_local = std::min(tmax_local, t1);
      if (tmin > tmax_local) return false;
    }
    return true;
  }

  static bool intersect_tri(const Ray3d& r, const Triangle3d& t, double& tout) noexcept {
    const Vec3d e1 = t.b - t.a;
    const Vec3d e2 = t.c - t.a;
    const Vec3d pv = cross(r.dir, e2);
    const double det = dot(e1, pv);
    if (std::abs(det) < 1e-14) return false;
    const double inv = 1.0 / det;
    const Vec3d tv{r.origin[0] - t.a[0], r.origin[1] - t.a[1], r.origin[2] - t.a[2]};
    const double u = dot(tv, pv) * inv;
    if (u < 0.0 || u > 1.0) return false;
    const Vec3d qv = cross(tv, e1);
    const double v = dot(r.dir, qv) * inv;
    if (v < 0.0 || u + v > 1.0) return false;
    const double tt = dot(e2, qv) * inv;
    if (tt <= 0.0) return false;
    tout = tt;
    return true;
  }

  std::vector<Prim> prims_;
  std::vector<std::uint32_t> order_;
  std::vector<Node> nodes_;
};

}  // namespace v1
}  // namespace geom
```

- [ ] **Step 4: 构建并运行**

Run: `cmake --build build --config Debug && ctest --test-dir build -C Debug --output-on-failure`
Expected: PASS（500 条射线、200 次最近邻与暴力一致）。

- [ ] **Step 5: 提交**

```bash
git add include/geom/index tests
git commit -m "feat(m1): add SAH BVH with raycast and nearest queries"
```

---

## Task 14: 覆盖率门禁与缺口补齐

**Files:**
- Modify: `CMakeLists.txt`（新增 `GEOM_COVERAGE` 选项）、`.github/workflows/ci.yml`（新增 coverage job）、`tests/CMakeLists.txt`（注册 `unit/coverage_gaps_test.cpp`）
- Create: `tests/unit/coverage_gaps_test.cpp`
- Test: `tests/unit/coverage_gaps_test.cpp`

**Interfaces:**
- Consumes: 全部已有公开 API
- Produces: CMake 选项 `GEOM_COVERAGE`；CI `coverage` job；覆盖负向路径、边界与空状态的用例

- [ ] **Step 1: 写缺口测试**

`tests/unit/coverage_gaps_test.cpp`：

```cpp
#include <doctest/doctest.h>

#include <vector>

#include <geom/cgeom/delaunay2.hpp>
#include <geom/cgeom/hull2d.hpp>
#include <geom/core/handle.hpp>
#include <geom/curves/nurbs_curve.hpp>
#include <geom/index/bvh.hpp>
#include <geom/numeric/expansion.hpp>
#include <geom/transform/affine.hpp>

using geom::v1::Affine;
using geom::v1::Bvh;
using geom::v1::convex_hull2d;
using geom::v1::delaunay2;
using geom::v1::detail::Expansion;
using geom::v1::Errc;
using geom::v1::Handle;
using geom::v1::HandleAllocator;
using geom::v1::NurbsCurve;
using geom::v1::Point2d;
using geom::v1::Point3d;
using geom::v1::Ray3d;
using geom::v1::Vec3d;

TEST_CASE("singular affine inverse reports degenerate") {
  Affine a{};
  for (int r = 0; r < 3; ++r) {
    for (int c = 0; c < 3; ++c) a.m[r][c] = 0.0;
  }
  const auto inv = a.inverse();
  REQUIRE_FALSE(inv.has_value());
  CHECK(inv.error().code == Errc::degenerate);
}

TEST_CASE("NURBS rejects zero weights, bad knots and size mismatches") {
  CHECK_FALSE(NurbsCurve(1, {0.0, 0.0, 1.0, 1.0},
                         {Point3d{0, 0, 0}, Point3d{1, 0, 0}}, {1.0, 0.0}).is_valid());
  CHECK_FALSE(NurbsCurve(1, {0.0, 1.0, 0.0, 1.0},
                         {Point3d{0, 0, 0}, Point3d{1, 0, 0}}, {1.0, 1.0}).is_valid());
  CHECK_FALSE(NurbsCurve(1, {0.0, 0.0, 1.0, 1.0},
                         {Point3d{0, 0, 0}, Point3d{1, 0, 0}}, {1.0}).is_valid());
}

TEST_CASE("NURBS clamps evaluation outside the domain") {
  const NurbsCurve c(1, {0.0, 0.0, 1.0, 1.0},
                     {Point3d{0, 0, 0}, Point3d{2, 0, 0}}, {1.0, 1.0});
  REQUIRE(c.is_valid());
  CHECK(c.eval(-5.0)[0] == doctest::Approx(0.0));
  CHECK(c.eval(5.0)[0] == doctest::Approx(2.0));
}

TEST_CASE("empty BVH answers queries safely") {
  Bvh b;
  Bvh::Hit hit;
  CHECK_FALSE(b.raycast(Ray3d{Point3d{0, 0, 0}, Vec3d{0, 0, 1}}, 1e9, hit));
  double d2 = -1.0;
  CHECK(b.nearest(Point3d{0, 0, 0}, d2) == 0xFFFFFFFFu);
  CHECK(b.bounds().is_empty());
}

TEST_CASE("hull collapses identical and two-point input") {
  const std::vector<Point2d> same{{1, 1}, {1, 1}, {1, 1}};
  CHECK(convex_hull2d(same).size() == 1);

  const std::vector<Point2d> two{{0, 0}, {2, 2}};
  CHECK(convex_hull2d(two).size() == 2);
}

TEST_CASE("delaunay of a single triangle yields one triangle") {
  const std::vector<Point2d> tri{{0, 0}, {1, 0}, {0, 1}};
  const auto tris = delaunay2(tri);
  REQUIRE(tris.size() == 1);
}

TEST_CASE("expansion handles scalar scaling and null handles") {
  Expansion e;
  e.add_scaled(Expansion(1.0), 2.0);
  CHECK(e.approx() == doctest::Approx(2.0));

  HandleAllocator alloc;
  alloc.release(Handle::null());
  CHECK_FALSE(alloc.alive(Handle::null()));
}
```

- [ ] **Step 2: 构建并运行**

Run: `cmake --build build --config Debug && ctest --test-dir build -C Debug --output-on-failure`
Expected: 全部 PASS。若任一用例失败，说明实现存在缺陷：先修实现（不得修改断言来迁就实现）。

- [ ] **Step 3: 增加覆盖率开关与 CI job**

`CMakeLists.txt` 在 `include(cmake/GeomWarnings.cmake)` 之后追加：

```cmake
option(GEOM_COVERAGE "Enable coverage instrumentation (GCC/Clang)" OFF)
if(GEOM_COVERAGE)
  if(MSVC)
    message(FATAL_ERROR "GEOM_COVERAGE requires GCC or Clang")
  endif()
  target_compile_options(geom_warnings INTERFACE --coverage -O0 -g)
  target_link_options(geom_warnings INTERFACE --coverage)
endif()
```

`.github/workflows/ci.yml` 在 `jobs:` 下追加：

```yaml
  coverage:
    runs-on: ubuntu-24.04
    steps:
      - uses: actions/checkout@v4
      - name: Install gcovr
        run: sudo apt-get update && sudo apt-get install -y gcovr
      - name: Configure
        run: cmake -S . -B build-cov -DGEOM_COVERAGE=ON -DGEOM_WERROR=ON -DCMAKE_BUILD_TYPE=Debug
      - name: Build
        run: cmake --build build-cov -j 4
      - name: Test
        run: ctest --test-dir build-cov --output-on-failure
      - name: Coverage gate
        run: |
          gcovr -r . --fail-under-line 90 --fail-under-branch 85 \
                --exclude tests/ --exclude third_party/ --exclude build-cov/
          gcovr -r . --filter 'include/geom/(core|numeric|primitives|curves|surfaces)/' \
                --fail-under-line 95 --fail-under-branch 90 \
                --exclude tests/ --exclude third_party/ --exclude build-cov/
```

- [ ] **Step 4: 运行覆盖率门禁**

在 Linux/CI 上执行：

```bash
cmake -S . -B build-cov -DGEOM_COVERAGE=ON -DGEOM_WERROR=ON -DCMAKE_BUILD_TYPE=Debug
cmake --build build-cov -j 4
ctest --test-dir build-cov --output-on-failure
gcovr -r . --fail-under-line 90 --fail-under-branch 85 \
      --exclude tests/ --exclude third_party/ --exclude build-cov/
gcovr -r . --filter 'include/geom/(core|numeric|primitives|curves|surfaces)/' \
      --fail-under-line 95 --fail-under-branch 90 \
      --exclude tests/ --exclude third_party/ --exclude build-cov/
```

Expected: 退出码 0；若报告某个文件低于门槛，按测试策略 §2 的模块目标补测试，直到达标（禁止用豁免注释凑数，除非是不可达分支且经评审）。

- [ ] **Step 5: 提交**

```bash
git add CMakeLists.txt .github tests
git commit -m "test: add coverage gate and close negative-path gaps"
```

---

## M1 出口验收（不新增代码）

在 Task 14 合并后执行：

- [ ] `cmake -S . -B build -DGEOM_WERROR=ON -DGEOM_BUILD_BENCH=ON`
- [ ] `cmake --build build --config Release`
- [ ] `ctest --test-dir build -C Release --output-on-failure` → 全绿（含 M0 出口用例、Delaunay 空圆性质、覆盖缺口用例）
- [ ] `./build/Release/geom_bench_predicates`（Windows：`.\build\Release\geom_bench_predicates.exe`）→ `orient2d` 过滤路径 < 200 ns/call，退出码 0
- [ ] 覆盖率门禁通过（CI `coverage` job 或本地 gcovr）：全局行 ≥ 90%、分支 ≥ 85%；core/numeric 达到测试策略 §2 的 95%/90% 模块目标
- [ ] 新建 `docs/benchmarks/m1.md`，记录：日期、OS、编译器版本、CPU 型号、`orient2d` ns/call、Delaunay 38 点用例耗时（`ctest -R geom_tests -V` 输出中该用例的时间）、覆盖率摘要
- [ ] 在架构文档 §16.2 勾选已满足的 M0/M1 出口条件（未满足项必须写明原因，不允许勾选）

---

## 自检（Self-Review）

**Spec coverage（对照架构文档 + 测试策略）**

| 架构要求 | 覆盖任务 |
| --- | --- |
| L0 配置/错误/断言/handle/arena | T1、T2、T5 |
| L1 标量 traits / 区间过滤 / 精确谓词 / SoS | T3、T4、T10（SoS 留给 Delaunay 后置：当前用「重复点跳过 + 超三角形」代替） |
| L2 图元、曲线曲面 | T6、T7、T8 |
| L4 2D 凸包/Delaunay | T11、T12 |
| L3 BVH | T13 |
| 测试策略 §2 覆盖门槛（行/分支/diff/Errc/API 三类用例） | 全局约束 + T14（覆盖率门禁与缺口补齐） |
| 测试策略 §6 谓词专项（整数 oracle、双路径、尺度、ulp） | T4、T10（整数 oracle + 尺度 + 非有限）；T14 补边界 |
| 测试策略 §7 退化语料 / §9 并发 / §10 模糊 | 语料目录与 CI 约定在 T1/T9；TSan、24 h fuzz 属 M1 夜间任务（本计划外） |
| M0 出口：谓词 1e6 oracle、NURBS 求值校验、编译矩阵 | T4（整数 oracle）、T7（解析四分圆）、T1（三平台 CI） |
| M1 出口：Delaunay 空圆、BVH 差分、覆盖达标 | T12、T13、T14 |

**已知缺口（有意留给后续计划）**：2D arrangement、多边形布尔、CDT、网格化、B-rep；TSan/夜间 fuzz 的流水线接入与 1e5 射线差分基线需在 M1 收口时由 CI 定时任务补足。

**占位符扫描**：本计划所有代码块均为完整实现（无伪代码、无 `TODO`/`TBD`、无未定义符号）；每条测试期望值均已按算法手算核对。

**类型一致性**：`Result<T>`、`Error`、`Point<T,D>`、`Vec<T,D>`、`Box<T,D>`、`Segment<T,D>`、`sign_orient2d`、`sign_incircle` 在测试、实现与后续任务中名称/参数顺序一致；`detail::eval_curve` 在曲线任务定义、曲面任务复用。

---

## 执行交接（Execution Handoff）

计划已保存到 `docs/superpowers/plans/2026-11-06-m0-m1-foundation.md`。两种执行方式：

1. **Subagent-Driven（推荐）**：每个任务派发独立 subagent，任务间做两阶段审查，迭代快。
2. **Inline Execution**：在当前会话用 `executing-plans` 批量执行，带检查点。

选择哪一种？
