# GeomCore M1b 计算几何补全 —— 实施计划

> **For agentic workers:** REQUIRED SUB-SKILL: Use subagent-driven-development (recommended) or executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** 补齐 M1 延期的四项能力——约束 Delaunay（CDT）、简单多边形三角剖分（含洞）、基于 arrangement 的 2D 多边形布尔、KD 树——并让 QA-1/QA-2/QA-4 成为可复现的验收门禁。

**Architecture:** 在既有 `geom::v1` header-only 体系内自底向上推进：先定义共享的 2D 多边形数据模型，再实现互不依赖的 KD 树、三角剖分与 arrangement（边求交 → 平面细分 → 逐面 winding），最后由 arrangement 驱动布尔运算；CDT 复用 M1a 的 `delaunay2` 与精确谓词。所有拓扑判定只用 `sign_orient2d` / `sign_incircle` 的符号；坐标输出允许 double 近似。

**Tech Stack:** C++20、CMake ≥ 3.25、doctest（仅测试）、MSVC 19.3x / GCC 12+ / Clang 15+。

**Spec:** `docs/architecture/geometry-library-architecture.md`（v2.1，§8.1–8.4、§12.3、§16.2–16.3 为本计划的上位契约）与 `docs/quality/testing-strategy.md`（覆盖门槛、性质测试、性能门禁）。

## Global Constraints

- 标准：C++20；`CMAKE_CXX_EXTENSIONS OFF`；最低编译器 MSVC 19.3x、GCC 12、Clang 15。
- 浮点：禁止 `-ffast-math`；谓词路径保持 `/fp:precise` 与 `-ffp-contract=off`。
- 依赖：核心零第三方依赖；doctest 只允许出现在 `tests/`。
- 命名空间：`geom::v1`；头文件 `#pragma once`；宏前缀 `GEOM_`。
- 公开头文件必须自包含、可单独编译；单文件目标 ≤ 400 行，超过时按职责拆到 `detail/` 头。
- 错误：可失败操作返回 `Result<T>`；禁止异常做控制流；`noexcept` 只用于确实不失败的函数。
- 坐标契约：所有新 API 要求有限坐标；谓词可解析的指数跨度上限沿用 `kMaxSpan = 1314`（`predicates2.hpp`），超出时返回 `Errc::unsupported` 或文档化的空结果，绝不返回部分拓扑。
- 退化契约：不实现 SoS。CDT 的非约束点仍要求一般位置（无三点共线、无共圆）；布尔允许自交输入，但输出必须规则化。
- 测试：每个新测试文件追加到 `tests/CMakeLists.txt` 的 `GEOM_TEST_SOURCES`；提交前 Debug 与 Release `ctest` 全绿；CI 用 `GEOM_WERROR=ON`。
- 覆盖门槛（测试策略 §2）：新增代码 diff 覆盖 100%；`cgeom`/`index` 行 ≥ 90%、分支 ≥ 85%；每个 `Errc` 至少一个触发用例；每条公开 API 正例、边界、失败路径齐全。
- 随机测试固定种子并打印种子；断言精确值或 `doctest::Approx` 的明确容差，禁止 `>= 0` 式弱断言。
- 提交粒度：每个任务一次提交，信息格式 `feat(m1b): ...` / `test(m1b): ...` / `docs(m1b): ...`。

## 范围与拆分（Scope Check）

本计划覆盖 M1b 全部四项交付物与出口验收：

- **T1** 共享 `Polygon2d` 模型、精确环方向、面积与结构校验；
- **T2** `KdTree`（build / nearest / k-nearest / range）；
- **T3** 简单多边形三角剖分（含洞，桥接 + 耳切）；
- **T4** 线段求交与基本边拆分（arrangement 的前半段）；
- **T5** 平面细分 DCEL、逐面 winding 与点定位；
- **T6** `boolean2` 四种运算、三种 fill rule、输出规则化；
- **T7** `cdt2` 约束插入与恢复；
- **T8** QA-1/QA-2/QA-4 基准与并发测试、文档收口。

**不在本计划**：B-rep、SSI、网格布尔（M2+）；Bentley–Ottmann 全事件扫描线；snap rounding / 整数模式；SoS；3D Delaunay；R-tree；夜间 fuzz 流水线（沿用 CI 定时任务，另立计划）。

### 与架构文档的偏差（需评审）

1. **三角剖分**：架构 §8.2 写「单调多边形分解 + 耳切」。本计划先实现**桥接 + 耳切（O(n²)）**。理由：M1b 无三角剖分性能门禁，先以可证明正确的简单算法落地；API 不变，单调分解可在 M2 面网格化出现性能压力时替换实现。风险：超大面（>5 万顶点）三角剖分偏慢。
2. **arrangement 求交**：架构 §8.3 写「扫描线建 arrangement」。本计划用**按 x 排序的活动集宽相位 + 精确谓词求交**（典型 O((n+c) log n)，最坏 O(n²)），再建 DCEL。理由：全事件 Bentley–Ottmann 的退化处理量远大于本里程碑可验证范围；精确符号仍是唯一判定依据。QA-1 的 80 ms 门禁就是该选择是否成立的判据；若不达标，升级为扫描线另立计划。
3. **自交检测**：`allow_self_intersections=false` 时用求交阶段记录的原始段对判定，不做二次扫描。

## Review Focus

以下五类输入是规格没有逐条展开、但最可能让使用者踩坑的地方；每个都在对应任务中有断言：

1. **共线/重合边**（两个多边形共享整条边）：布尔输出不得出现内部边、重复边或零面积环——T6 共享边用例。
2. **垂直边、T 型接点、三线共点**：arrangement 顶点合并与角度排序必须稳定——T4/T5 用例。
3. **洞的嵌套与归属**：差集产生洞、洞必须挂到最小包含外环——T6 嵌套用例。
4. **重复点与零长度边**：不得产生退化三角形或除零——T1/T3/T4 用例。
5. **极端指数**（1e200 与 1e-200 同现）：必须显式拒绝（`Errc::unsupported`）而不是静默产生错误拓扑——T6 用例。

## 文件结构（File Map）

```
include/geom/cgeom/polygon2.hpp       Polygon2d / PolygonSet / 面积 / 环方向 / 校验
include/geom/cgeom/triangulate2.hpp   简单多边形（含洞）三角剖分
include/geom/cgeom/arrangement2.hpp   ArrangementSegment / 线段拆分 / Arrangement2
include/geom/cgeom/boolean2.hpp       BooleanOp / FillRule / boolean2
include/geom/cgeom/cdt2.hpp           CdtConstraint / cdt2
include/geom/index/kdtree.hpp         KdTree
tests/unit/polygon2_test.cpp
tests/unit/kdtree_test.cpp
tests/unit/triangulate2_test.cpp
tests/unit/arrangement2_test.cpp
tests/unit/boolean2_test.cpp
tests/unit/cdt2_test.cpp
tests/unit/qa_concurrency_test.cpp    QA-2 一致性（TSan 下同一用例）
bench/bench_m1b.cpp                   QA-1 / QA-4
docs/benchmarks/m1b.md
docs/adr/0020-m1b-arrangement-and-triangulation.md
```

修改：`tests/CMakeLists.txt`（测试注册）、`CMakeLists.txt`（bench 目标、`GEOM_TSAN`）、`include/geom/core/error.hpp`（新增 `Errc::self_intersecting`）、`docs/architecture/geometry-library-architecture.md`（§16.3 状态）。

---

## Task 1: 2D 多边形模型与精确环方向

**Files:**
- Create: `include/geom/cgeom/polygon2.hpp`
- Test: `tests/unit/polygon2_test.cpp`
- Modify: `tests/CMakeLists.txt`

**Interfaces:**
- Consumes: `geom::Point2d`、`geom::Errc`、`sign_orient2d`。
- Produces:
  - `struct Polygon2d { std::vector<Point2d> vertices; std::vector<std::uint32_t> ring_starts; std::uint32_t ring_count() const noexcept; std::span<const Point2d> ring(std::uint32_t i) const noexcept; };`
  - `struct PolygonSet { std::vector<Polygon2d> polygons; };`
  - `Errc validate_polygon2(const Polygon2d&) noexcept;`
  - `int ring_orientation2(std::span<const Point2d> ring) noexcept;`（+1 CCW / -1 CW / 0 退化）
  - `double signed_area2(std::span<const Point2d> ring) noexcept;`
  - `double polygon_area2(const Polygon2d&) noexcept;`
  - `double polygon_set_area2(const PolygonSet&) noexcept;`

- [x] **Step 1: 写失败测试**

`tests/unit/polygon2_test.cpp`：

```cpp
#include <doctest/doctest.h>

#include <geom/cgeom/polygon2.hpp>

#include <limits>

using geom::v1::Errc;
using geom::v1::Polygon2d;
using geom::v1::Point2d;

namespace {

Polygon2d square_with_hole() {
  Polygon2d p;
  // 4x4 outer CCW, 2x2 hole CW.
  p.vertices = {{0, 0}, {4, 0}, {4, 4}, {0, 4},
                {1, 1}, {1, 3}, {3, 3}, {3, 1}};
  p.ring_starts = {0, 4, 8};
  return p;
}

}  // namespace

TEST_CASE("polygon2: ring iteration, orientation and area") {
  const Polygon2d poly = square_with_hole();
  CHECK(poly.ring_count() == 2);
  CHECK(poly.ring(0).size() == 4);
  CHECK(poly.ring(1).size() == 4);
  CHECK(geom::v1::ring_orientation2(poly.ring(0)) == 1);
  CHECK(geom::v1::ring_orientation2(poly.ring(1)) == -1);
  CHECK(geom::v1::signed_area2(poly.ring(0)) == doctest::Approx(16.0));
  CHECK(geom::v1::signed_area2(poly.ring(1)) == doctest::Approx(-4.0));
  CHECK(geom::v1::polygon_area2(poly) == doctest::Approx(12.0));
}

TEST_CASE("polygon2: area is independent of winding") {
  Polygon2d p;
  p.vertices = {{0, 0}, {0, 1}, {1, 1}, {1, 0}};  // CW outer
  p.ring_starts = {0, 4};
  CHECK(geom::v1::ring_orientation2(p.ring(0)) == -1);
  CHECK(geom::v1::polygon_area2(p) == doctest::Approx(1.0));
  CHECK(geom::v1::ring_orientation2(std::vector<Point2d>{{0, 0}, {1, 1}}) == 0);
}

TEST_CASE("polygon2: duplicate consecutive points do not change orientation") {
  Polygon2d p;
  p.vertices = {{0, 0}, {0, 0}, {1, 0}, {1, 1}, {0, 1}};
  p.ring_starts = {0, 5};
  CHECK(geom::v1::ring_orientation2(p.ring(0)) == 1);
}

TEST_CASE("polygon2: structural validation") {
  Polygon2d p;
  p.vertices = {{0, 0}, {1, 0}, {0, 1}};
  p.ring_starts = {0, 3};
  CHECK(geom::v1::validate_polygon2(p) == Errc::none);

  p.ring_starts = {0, 2, 3};
  CHECK(geom::v1::validate_polygon2(p) == Errc::degenerate);

  p.ring_starts = {0, 4};
  CHECK(geom::v1::validate_polygon2(p) == Errc::degenerate);

  p.ring_starts = {1, 3};
  CHECK(geom::v1::validate_polygon2(p) == Errc::degenerate);

  p.ring_starts = {0, 3};
  p.vertices[0][0] = std::numeric_limits<double>::quiet_NaN();
  CHECK(geom::v1::validate_polygon2(p) == Errc::non_finite);

  p.vertices[0][0] = std::numeric_limits<double>::infinity();
  CHECK(geom::v1::validate_polygon2(p) == Errc::non_finite);
}

TEST_CASE("polygon2: set area sums polygons") {
  geom::v1::PolygonSet set;
  set.polygons.push_back(square_with_hole());
  set.polygons.push_back(square_with_hole());
  CHECK(geom::v1::polygon_set_area2(set) == doctest::Approx(24.0));
}
```

- [x] **Step 2: 运行确认失败**

```powershell
cmake -S . -B build -DGEOM_WERROR=ON
cmake --build build --config Debug
```

Expected: 编译失败，`geom/cgeom/polygon2.hpp: No such file or directory`。

- [x] **Step 3: 实现 `polygon2.hpp`**

实现要点（签名见 Interfaces，全部 `inline`/`noexcept`）：

- `ring_count()`：`ring_starts.size()` 为 0 时返回 0，否则 `size()-1`；`ring(i)` 返回 `std::span{vertices.data()+ring_starts[i], ring_starts[i+1]-ring_starts[i]}`。
- `validate_polygon2`：先逐点 `std::isfinite`（任一非有限 → `Errc::non_finite`）；再查 `ring_starts.size() >= 2`、`front()==0`、`back()==vertices.size()`、严格递增、每个环 ≥3 点；任一不满足 → `Errc::degenerate`；否则 `none`。
- `ring_orientation2`：找字典序最小（x 相同比 y）的顶点 `k`；分别沿环向前/向后跳过与 `p[k]` 完全相等的重复点得到 `prev/next`；若三者不能构成两个不同方向 → 0；否则返回 `sign_orient2d(prev, k, next)`。该符号在极值顶点处与整环方向一致（简单环；自交环未定义）。
- `signed_area2`：`long double` 累加 `x_i*y_{i+1} - x_{i+1}*y_i` 后乘 0.5 转 `double`；这不是精确谓词，仅用于报告与测试容差。
- `polygon_area2`：`|outer| - Σ|hole|`（用 `signed_area2` 的绝对值），保证与环绕方向无关。
- `polygon_set_area2`：累加。

- [x] **Step 4: 运行确认通过**

```powershell
cmake --build build --config Debug
ctest --test-dir build -C Debug --output-on-failure
```

Expected: `geom_tests` 通过；`polygon2` 用例 5 个全绿。

- [x] **Step 5: 提交**

```powershell
git add include/geom/cgeom/polygon2.hpp tests/unit/polygon2_test.cpp tests/CMakeLists.txt
git commit -m "feat(m1b): add 2D polygon model with exact ring orientation"
```

---

## Task 2: KD 树

**Files:**
- Create: `include/geom/index/kdtree.hpp`
- Test: `tests/unit/kdtree_test.cpp`
- Modify: `tests/CMakeLists.txt`

**Interfaces:**
- Consumes: `Point3d`、`Box3d`、`Box3d::empty()` / `expand` / `lo` / `hi`。
- Produces:

```cpp
class KdTree {
 public:
  void build(std::span<const Point3d> points);
  std::uint32_t size() const noexcept;
  const Box3d& bounds() const noexcept;
  // 返回最近点索引；空树返回 0xFFFFFFFF 且 best_d2 = +inf。
  // 距离相等时取最小索引。
  std::uint32_t nearest(const Point3d& p, double& best_d2) const;
  // 写入至多 min(k, out.size()) 个索引，按 (d2, index) 升序；返回写入数。
  std::uint32_t nearest_k(const Point3d& p, std::uint32_t k,
                          std::span<std::uint32_t> out) const;
  // 盒内点（含边界）的索引，升序；空树返回空。
  std::vector<std::uint32_t> range(const Box3d& box) const;
};
```

- [x] **Step 1: 写失败测试**

`tests/unit/kdtree_test.cpp`：

```cpp
#include <doctest/doctest.h>

#include <geom/index/kdtree.hpp>

#include <algorithm>
#include <cmath>
#include <limits>
#include <random>
#include <vector>

using geom::v1::Box3d;
using geom::v1::KdTree;
using geom::v1::Point3d;

namespace {

double dist2(const Point3d& a, const Point3d& b) {
  const double dx = a[0] - b[0];
  const double dy = a[1] - b[1];
  const double dz = a[2] - b[2];
  return dx * dx + dy * dy + dz * dz;
}

}  // namespace

TEST_CASE("kdtree: empty tree") {
  KdTree tree;
  double d2 = 0.0;
  CHECK(tree.size() == 0);
  CHECK(tree.bounds().lo[0] > tree.bounds().hi[0]);  // empty box
  CHECK(tree.nearest(Point3d{0, 0, 0}, d2) == 0xFFFFFFFFu);
  CHECK(std::isinf(d2));
  std::uint32_t out[3]{};
  CHECK(tree.nearest_k(Point3d{0, 0, 0}, 3, out) == 0);
  CHECK(tree.range(Box3d::empty()).empty());
}

TEST_CASE("kdtree: single point and duplicate coordinates") {
  std::vector<Point3d> pts{{1, 1, 1}, {1, 1, 1}};
  KdTree tree;
  tree.build(pts);
  double d2 = 0.0;
  CHECK(tree.nearest(Point3d{0, 0, 0}, d2) == 0);
  CHECK(d2 == doctest::Approx(3.0));
  CHECK(tree.size() == 2);

  std::vector<std::uint32_t> out(2);
  CHECK(tree.nearest_k(Point3d{0, 0, 0}, 2, out) == 2);
  CHECK(out[0] == 0);
  CHECK(out[1] == 1);
}

TEST_CASE("kdtree: nearest, k-nearest and range match brute force") {
  std::mt19937 rng(12345);
  std::uniform_real_distribution<double> uni(-100.0, 100.0);
  std::vector<Point3d> pts(2000);
  for (Point3d& p : pts) p = Point3d{uni(rng), uni(rng), uni(rng)};

  KdTree tree;
  tree.build(pts);
  CHECK(tree.size() == pts.size());

  for (int q = 0; q < 200; ++q) {
    const Point3d query{uni(rng), uni(rng), uni(rng)};
    double best = std::numeric_limits<double>::infinity();
    std::uint32_t bi = 0xFFFFFFFFu;
    for (std::uint32_t i = 0; i < pts.size(); ++i) {
      const double d = dist2(query, pts[i]);
      if (d < best || (d == best && i < bi)) {
        best = d;
        bi = i;
      }
    }
    double got = 0.0;
    CHECK(tree.nearest(query, got) == bi);
    CHECK(got == best);

    std::vector<std::uint32_t> all(pts.size());
    for (std::uint32_t i = 0; i < pts.size(); ++i) all[i] = i;
    std::sort(all.begin(), all.end(), [&](std::uint32_t a, std::uint32_t b) {
      const double da = dist2(query, pts[a]);
      const double db = dist2(query, pts[b]);
      return da != db ? da < db : a < b;
    });
    std::vector<std::uint32_t> out(9);
    const std::uint32_t n = tree.nearest_k(query, 9, out);
    REQUIRE(n == 9);
    for (std::uint32_t i = 0; i < n; ++i) CHECK(out[i] == all[i]);
  }

  const Box3d box{Point3d{-50, -50, -50}, Point3d{50, 50, 50}};
  std::vector<std::uint32_t> want;
  for (std::uint32_t i = 0; i < pts.size(); ++i) {
    const bool in = pts[i][0] >= box.lo[0] && pts[i][0] <= box.hi[0] &&
                    pts[i][1] >= box.lo[1] && pts[i][1] <= box.hi[1] &&
                    pts[i][2] >= box.lo[2] && pts[i][2] <= box.hi[2];
    if (in) want.push_back(i);
  }
  CHECK(tree.range(box) == want);
}
```

- [x] **Step 2: 运行确认失败**

```powershell
cmake --build build --config Debug
```

Expected: 找不到 `geom/index/kdtree.hpp`。

- [x] **Step 3: 实现 `kdtree.hpp`**

- 存储：`points_`（构建时拷贝）、`order_`（索引排列）、`nodes_`；`Node { Box3d box; std::int32_t left, right; std::uint32_t first, count; }`，叶子 `count>0`，内部 `count==0`。
- `build`：空输入清空并返回；否则递归构建：节点盒 = 点盒并集；叶子阈值 8；按盒最长轴（相等取 x→y→z）取中位点分割，`nth_element` 比较键 `(coord[axis], index)` 保证确定性；深度上限 64，到达上限强制叶子。
- `nearest` / `nearest_k`：显式栈（深度 ≤64，数组 128）；`box_distance2` 剪枝；`nearest` 的平局按索引；`nearest_k` 用大小为 k 的最大堆（`std::push_heap`，比较 `(d2, index)` 的“更差”），剪枝条件 `box_d2 > worst_d2`（相等时继续下降，因为更小索引可能胜出）。
- `range`：遍历与盒相交的节点，收集后 `std::sort` 升序。
- 文档注明：坐标有限；空树语义；`range` 线性于输出大小 + 访问节点数。

- [x] **Step 4: 运行确认通过**

```powershell
cmake --build build --config Debug
ctest --test-dir build -C Debug --output-on-failure
```

Expected: 3 个 `kdtree` 用例全绿。

- [x] **Step 5: 提交**

```powershell
git add include/geom/index/kdtree.hpp tests/unit/kdtree_test.cpp tests/CMakeLists.txt
git commit -m "feat(m1b): add KD tree with nearest, k-nearest and range queries"
```

---

## Task 3: 简单多边形三角剖分（含洞）

**Files:**
- Create: `include/geom/cgeom/triangulate2.hpp`
- Test: `tests/unit/triangulate2_test.cpp`
- Modify: `tests/CMakeLists.txt`

**Interfaces:**
- Consumes: `Polygon2d`、`validate_polygon2`、`ring_orientation2`、`sign_orient2d`、`polygon_area2`。
- Produces:

```cpp
struct Triangulation2 {
  std::vector<std::array<std::uint32_t, 3>> triangles;  // CCW，索引到 poly.vertices
  bool ok = false;
};
Triangulation2 triangulate_polygon2(const Polygon2d& poly);
```

前置条件（文档化）：`validate_polygon2(poly) == Errc::none`；各环简单、互不相交；洞严格位于外环内部且互不嵌套；环可任意方向；允许重复的相邻点（内部去重）。违反前置条件或校验失败时返回 `{ {}, false }`。

- [x] **Step 1: 写失败测试**

`tests/unit/triangulate2_test.cpp`：

```cpp
#include <doctest/doctest.h>

#include <geom/cgeom/polygon2.hpp>
#include <geom/cgeom/triangulate2.hpp>
#include <geom/numeric/predicates2.hpp>

#include <cmath>
#include <limits>

using geom::v1::Polygon2d;
using geom::v1::Point2d;
using geom::v1::Triangulation2;
using geom::v1::polygon_area2;
using geom::v1::triangulate_polygon2;

namespace {

double triangle_area2(const Polygon2d& p, const std::array<std::uint32_t, 3>& t) {
  return std::fabs(geom::v1::signed_area2(
      std::vector<Point2d>{p.vertices[t[0]], p.vertices[t[1]], p.vertices[t[2]]}));
}

double covered_area(const Polygon2d& p, const Triangulation2& tri) {
  double sum = 0.0;
  for (const auto& t : tri.triangles) sum += triangle_area2(p, t);
  return sum;
}

Polygon2d single(std::vector<Point2d> ring) {
  Polygon2d p;
  p.vertices = std::move(ring);
  p.ring_starts = {0, static_cast<std::uint32_t>(p.vertices.size())};
  return p;
}

}  // namespace

TEST_CASE("triangulate2: square") {
  const Polygon2d p = single({{0, 0}, {1, 0}, {1, 1}, {0, 1}});
  const Triangulation2 tri = triangulate_polygon2(p);
  REQUIRE(tri.ok);
  CHECK(tri.triangles.size() == 2);
  for (const auto& t : tri.triangles) {
    CHECK(geom::v1::sign_orient2d(p.vertices[t[0]][0], p.vertices[t[0]][1],
                                  p.vertices[t[1]][0], p.vertices[t[1]][1],
                                  p.vertices[t[2]][0], p.vertices[t[2]][1]) > 0);
  }
  CHECK(covered_area(p, tri) == doctest::Approx(1.0));
}

TEST_CASE("triangulate2: concave L-shape") {
  const Polygon2d p = single({{0, 0}, {2, 0}, {2, 1}, {1, 1}, {1, 2}, {0, 2}});
  const Triangulation2 tri = triangulate_polygon2(p);
  REQUIRE(tri.ok);
  CHECK(tri.triangles.size() == 4);
  CHECK(covered_area(p, tri) == doctest::Approx(3.0));
}

TEST_CASE("triangulate2: polygon with a hole") {
  Polygon2d p;
  p.vertices = {{0, 0}, {4, 0}, {4, 4}, {0, 4},
                {1, 1}, {1, 3}, {3, 3}, {3, 1}};
  p.ring_starts = {0, 4, 8};
  const Triangulation2 tri = triangulate_polygon2(p);
  REQUIRE(tri.ok);
  CHECK(tri.triangles.size() == 8);  // n + 2h - 2
  CHECK(covered_area(p, tri) == doctest::Approx(12.0));
}

TEST_CASE("triangulate2: reversed rings and collinear vertices") {
  Polygon2d reversed;
  reversed.vertices = {{0, 0}, {0, 1}, {1, 1}, {1, 0}};  // CW outer
  reversed.ring_starts = {0, 4};
  const Triangulation2 tri_cw = triangulate_polygon2(reversed);
  REQUIRE(tri_cw.ok);
  CHECK(covered_area(reversed, tri_cw) == doctest::Approx(1.0));

  const Polygon2d collinear = single({{0, 0}, {1, 0}, {2, 0}, {2, 2}, {0, 2}});
  const Triangulation2 tri_col = triangulate_polygon2(collinear);
  REQUIRE(tri_col.ok);
  CHECK(tri_col.triangles.size() == 2);  // collinear vertex (1,0) is dropped
  CHECK(covered_area(collinear, tri_col) == doctest::Approx(4.0));
}

TEST_CASE("triangulate2: rejects malformed input") {
  Polygon2d p;
  p.ring_starts = {0, 2};
  p.vertices = {{0, 0}, {1, 0}};
  CHECK(triangulate_polygon2(p).ok == false);

  p.vertices[0][0] = std::numeric_limits<double>::quiet_NaN();
  p.ring_starts = {0, 3};
  CHECK(triangulate_polygon2(p).ok == false);
}
```

- [x] **Step 2: 运行确认失败**

```powershell
cmake --build build --config Debug
```

Expected: 找不到 `triangulate2.hpp`。

- [x] **Step 3: 实现 `triangulate2.hpp`**

实现要点：

1. 校验；把每个环转成工作顶点 `{Point2d p; std::uint32_t original;}`，去掉相邻完全重复点，并用 `sign_orient2d(prev, cur, next) == 0` 去掉共线点（保留原索引映射）。
2. 归一化：外环 CCW、洞 CW（`ring_orientation2`），否则反转工作序列。
3. 洞桥接：按洞的最右点（x 最大，平局取 y 最大）向 +x 射线，找当前外环上最近交点所在边；选该边 x 较大的端点作为桥接点 `b`，用「桥接段不与任何环边真相交（端点除外）」的精确测试确认可见；不可见时尝试另一端点，再失败则按环序尝试洞的其它顶点，全部失败返回 `ok=false`。把 `b` 与洞最右点 `m` 各复制一次接入外环（`b` 在环中出现两次）。
4. 耳切：在合并后的工作环上循环找耳——`sign_orient2d(prev, cur, next) > 0`，且没有其它顶点严格位于三角形内或在三角形边界上（重复桥接点按坐标相等跳过）；裁掉耳、输出三角形（CCW，原始索引），直到剩 3 点输出最后一个。无耳可裁时返回 `ok=false`。
5. `ok=true` 当且仅当输出的每个三角形非退化且总面积与 `polygon_area2` 的相对误差 ≤ 1e-9（内部自检）。
6. 复杂度 O(n²)，文档注明；`triangles` 不排序（按裁耳顺序），调用方自行规范化。

- [x] **Step 4: 运行确认通过**

```powershell
cmake --build build --config Debug
ctest --test-dir build -C Debug --output-on-failure
```

Expected: 5 个 `triangulate2` 用例全绿。

- [x] **Step 5: 提交**

```powershell
git add include/geom/cgeom/triangulate2.hpp tests/unit/triangulate2_test.cpp tests/CMakeLists.txt
git commit -m "feat(m1b): triangulate simple polygons with holes"
```

---

## Task 4: 线段求交与基本边拆分

**Files:**
- Create: `include/geom/cgeom/arrangement2.hpp`
- Test: `tests/unit/arrangement2_test.cpp`（本任务先写拆分部分）
- Modify: `tests/CMakeLists.txt`

**Interfaces:**
- Consumes: `Point2d`、`sign_orient2d`、`kMaxSpan` 语义。
- Produces:

```cpp
struct ArrangementSegment {
  Point2d a{};
  Point2d b{};
  std::uint32_t source = 0;  // 环/操作数编号，用于 winding 聚合
};

namespace detail {
struct SplitSegment { Point2d a; Point2d b; std::uint32_t source; };
struct SplitResult {
  std::vector<SplitSegment> segments;                     // 拆分后的基本段
  std::vector<std::pair<std::uint32_t, std::uint32_t>> input_pairs;  // 相交的原段对 i<j
  double max_residual = 0.0;                              // 交点距支撑线的最大偏差
};
SplitResult split_segments(std::span<const ArrangementSegment> segments);
}  // namespace detail
```

`split_segments` 契约：输入段有限、非零长；输出按原段方向保序拆分；`input_pairs` 含共享端点/共线重叠/真相交的全部原段对；`max_residual` 是交点相对两条支撑线的最大垂距（坐标单位）。

- [x] **Step 1: 写失败测试**

`tests/unit/arrangement2_test.cpp`（本任务只包含这些用例；后续任务追加）：

```cpp
#include <doctest/doctest.h>

#include <geom/cgeom/arrangement2.hpp>

#include <algorithm>
#include <cmath>

using geom::v1::ArrangementSegment;
using geom::v1::Point2d;

namespace {

using geom::v1::detail::split_segments;
using geom::v1::detail::SplitResult;

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
```

- [x] **Step 2: 运行确认失败**

```powershell
cmake --build build --config Debug
```

Expected: 找不到 `arrangement2.hpp`。

- [x] **Step 3: 实现拆分**

- 过滤零长段（`a == b`）。宽相位：按 `min(x)` 排序，活动集 + `max(x)` 与 y 区间剪枝，枚举候选对。
- 精确判定：用 `sign_orient2d` 的四个符号区分相离、端点接触、真相交、共线重叠。真相交由 `long double` 计算参数 `t,u`，夹到 `[0,1]`，点 = `a + t*(b-a)` 取 double；接触点直接用已有端点坐标（保证精确相等）。
- 每条原段收集参数 `t`（含 0、1 与所有交点），排序去重（`t` 用 `long double` 比较，端点强制 0/1），按序生成子段，保留 `source` 与原方向。
- 共线重叠：把重叠端点（落在对方内部的）加入双方参数表；不合并同向重复段（保留 source 语义）。
- `max_residual`：每个新交点 `p` 到两条支撑线的垂距取大者，累计最大；用 `long double` 计算 `|cross(dir, p-a)| / |dir|`。
- `input_pairs`：每对相交（含端点接触/共线重叠）记录一次，排序去重。

- [x] **Step 4: 运行确认通过**

```powershell
cmake --build build --config Debug
ctest --test-dir build -C Debug --output-on-failure
```

Expected: 4 个 `split:` 用例全绿（段数/对数逐项精确）。

- [x] **Step 5: 提交**

```powershell
git add include/geom/cgeom/arrangement2.hpp tests/unit/arrangement2_test.cpp tests/CMakeLists.txt
git commit -m "feat(m1b): split segments at exact-predicate intersections"
```

---

## Task 5: 平面细分与逐面 winding

**Files:**
- Modify: `include/geom/cgeom/arrangement2.hpp`
- Test: `tests/unit/arrangement2_test.cpp`（追加）

**Interfaces:**
- Consumes: `detail::split_segments`、`sign_orient2d`。
- Produces:

```cpp
class Arrangement2 {
 public:
  struct HalfEdge {
    std::uint32_t origin = 0;
    std::uint32_t twin = 0;
    std::uint32_t next = 0;
    std::uint32_t face = 0;
  };
  struct Face {
    std::uint32_t edge = 0;   // 该面上的任一半边
    std::int32_t winding = 0; // 所有 source 的 winding 之和
    bool bounded = false;     // 面环面积为正（内部面）
  };

  static Arrangement2 build(std::span<const ArrangementSegment> segments);

  std::uint32_t vertex_count() const noexcept;
  std::uint32_t edge_count() const noexcept;    // 无向边数 = half_edges().size()/2
  std::uint32_t face_count() const noexcept;
  std::uint32_t source_count() const noexcept;
  std::span<const Point2d> vertices() const noexcept;
  std::span<const HalfEdge> half_edges() const noexcept;
  std::span<const Face> faces() const noexcept;
  std::int32_t winding(std::uint32_t face, std::uint32_t source) const noexcept;
  // 点所在最内层有界面；不在任何有界面（或在边界上）返回 0xFFFFFFFF。
  std::uint32_t locate(const Point2d& p) const noexcept;
  std::span<const std::pair<std::uint32_t, std::uint32_t>> input_pairs() const noexcept;
  double max_residual() const noexcept;
};
```

- [x] **Step 1: 写失败测试（追加到 `arrangement2_test.cpp`）**

```cpp
#include <geom/cgeom/polygon2.hpp>  // for the Polygon2d test fixtures

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
  CHECK(arr.winding(arr.locate(Point2d{0.75, 0.75}), 0) == 2);
  CHECK(arr.winding(arr.locate(Point2d{1.25, 0.75}), 1) == 1);
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
  CHECK(arr.locate(Point2d{1, 0}) == 0xFFFFFFFFu);  // 边界上无唯一面
}
```

- [x] **Step 2: 运行确认失败**

```powershell
cmake --build build --config Debug
```

Expected: `Arrangement2` 未定义。

- [x] **Step 3: 实现 DCEL 与 winding**

1. 调 `detail::split_segments`；把每个子段端点按坐标精确去重成顶点 id（排序 + 相等合并）。
2. 聚合无向边 `(min(u,v), max(u,v))`，保存 `delta_low_to_high[source]`：方向 `u→v` 时 low→high 记 +1，否则 -1。
3. 每个顶点收集出边，按极角排序：比较器先用上半平面/下半平面分组，再 `sign_orient2d(v, a, b) > 0`；同向重复出边合并（同一无向边只保留一条）。
4. 半边：每条无向边生成两条。`next(h)`：在 `h` 的终点处，取 `twin(h)` 在出边环序中的前一条（顺时针前驱），即标准 DCEL 的 face-left 约定。
5. 面：沿 `next` 遍历分配 face id，同时累加有向面积（鞋带，`long double`）；面积 > 0 记 bounded。
6. winding：所有面积 < 0 的环（无界侧）初始化为 0；沿边 BFS：对半边 `h`（u→v），`w[left] = w[right] + s * delta_low_to_high`，其中 `s = (h.origin < twin(h).origin) ? +1 : -1`。每个 source 独立。所有 source 的 winding 之和写入 `Face::winding`。
7. `locate`：对每个 bounded 环做精确射线穿越测试（用 `sign_orient2d` 判交点侧），收集包含 `p` 的环，取 |面积| 最小者；边界返回 `0xFFFFFFFF`。
8. `source_count()` = 最大 source + 1（空 arrangement 为 0）。

文档注明：输入应是闭合链（每个 source 的段首尾相接），否则 winding 只相对于初始化的无界面；面环模型把多连通面的每个边界环记为独立 face 记录（winding 相同），布尔只使用 bounded 面。

- [x] **Step 4: 运行确认通过**

```powershell
cmake --build build --config Debug
ctest --test-dir build -C Debug --output-on-failure
```

Expected: 8 个 `arrangement`/`split` 用例全绿；顶点/边/面计数与 winding 精确相等。

- [x] **Step 5: 提交**

```powershell
git add include/geom/cgeom/arrangement2.hpp tests/unit/arrangement2_test.cpp
git commit -m "feat(m1b): build planar arrangement with per-source winding"
```

---

## Task 6: 2D 多边形布尔

**Files:**
- Create: `include/geom/cgeom/boolean2.hpp`
- Modify: `include/geom/core/error.hpp`（新增 `Errc::self_intersecting`）
- Test: `tests/unit/boolean2_test.cpp`
- Modify: `tests/CMakeLists.txt`

**Interfaces:**
- Consumes: `Polygon2d`、`PolygonSet`、`validate_polygon2`、`Arrangement2`、`sign_orient2d`、`Result<T>`。
- Produces:

```cpp
enum class BooleanOp { union_, intersection, difference, symmetric_difference };
enum class FillRule { non_zero, even_odd, positive };

struct BooleanOptions {
  FillRule fill_rule = FillRule::non_zero;
  bool allow_self_intersections = true;
};

struct BooleanReport {
  std::uint32_t input_edges = 0;
  std::uint32_t output_polygons = 0;
  std::uint32_t merged_vertices = 0;  // 输出时因共线合并删除的顶点数
  double max_residual = 0.0;
};

struct BooleanResult {
  PolygonSet polygons;
  BooleanReport report;
};

Result<BooleanResult> boolean2(std::span<const Polygon2d> operands, BooleanOp op,
                               const BooleanOptions& options = {});
```

语义：`union_`/`intersection`/`symmetric_difference` 对全部操作数；`difference` 为 `operands[0] - union(其余)`。fill rule 逐操作数作用于其 winding。错误：`non_finite`、`degenerate`（空操作数/<3 点）、`unsupported`（指数跨度 > 1314）、`self_intersecting`（`allow_self_intersections=false` 且同一环的非相邻边相交）。

- [x] **Step 1: 写失败测试**

`tests/unit/boolean2_test.cpp`：

```cpp
#include <doctest/doctest.h>

#include <geom/cgeom/boolean2.hpp>
#include <geom/cgeom/polygon2.hpp>

#include <cmath>
#include <limits>

using geom::v1::BooleanOp;
using geom::v1::BooleanOptions;
using geom::v1::boolean2;
using geom::v1::Errc;
using geom::v1::FillRule;
using geom::v1::Polygon2d;
using geom::v1::PolygonSet;
using geom::v1::Point2d;
using geom::v1::polygon_set_area2;

namespace {

Polygon2d rect(double x0, double y0, double x1, double y1) {
  Polygon2d p;
  p.vertices = {{x0, y0}, {x1, y0}, {x1, y1}, {x0, y1}};
  p.ring_starts = {0, 4};
  return p;
}

double area_of(const PolygonSet& set) { return polygon_set_area2(set); }

}  // namespace

TEST_CASE("boolean2: overlapping squares, four operations") {
  const Polygon2d a = rect(0, 0, 1, 1);
  const Polygon2d b = rect(0.5, 0.5, 1.5, 1.5);
  const Polygon2d ops[2] = {a, b};

  auto uni = boolean2(ops, BooleanOp::union_);
  REQUIRE(uni.has_value());
  CHECK(area_of(uni.value().polygons) == doctest::Approx(1.75));
  CHECK(uni.value().report.output_polygons == 1);

  auto inter = boolean2(ops, BooleanOp::intersection);
  REQUIRE(inter.has_value());
  CHECK(area_of(inter.value().polygons) == doctest::Approx(0.25));

  auto diff = boolean2(ops, BooleanOp::difference);
  REQUIRE(diff.has_value());
  CHECK(area_of(diff.value().polygons) == doctest::Approx(0.75));

  auto xr = boolean2(ops, BooleanOp::symmetric_difference);
  REQUIRE(xr.has_value());
  CHECK(area_of(xr.value().polygons) == doctest::Approx(1.5));
}

TEST_CASE("boolean2: disjoint and nested") {
  const Polygon2d a = rect(0, 0, 1, 1);
  const Polygon2d far = rect(3, 3, 4, 4);
  const Polygon2d ops[2] = {a, far};
  auto uni = boolean2(ops, BooleanOp::union_);
  REQUIRE(uni.has_value());
  CHECK(uni.value().polygons.size() == 2);
  CHECK(area_of(uni.value().polygons) == doctest::Approx(2.0));
  auto inter = boolean2(ops, BooleanOp::intersection);
  REQUIRE(inter.has_value());
  CHECK(inter.value().polygons.polygons.empty());

  const Polygon2d outer = rect(0, 0, 4, 4);
  const Polygon2d inner = rect(1, 1, 3, 3);
  const Polygon2d nested[2] = {outer, inner};
  auto diff = boolean2(nested, BooleanOp::difference);
  REQUIRE(diff.has_value());
  CHECK(area_of(diff.value().polygons) == doctest::Approx(12.0));
  REQUIRE(diff.value().polygons.polygons.size() == 1);
  CHECK(diff.value().polygons.polygons[0].ring_count() == 2);
}

TEST_CASE("boolean2: shared edge produces one clean rectangle") {
  const Polygon2d a = rect(0, 0, 1, 1);
  const Polygon2d b = rect(1, 0, 2, 1);
  const Polygon2d ops[2] = {a, b};
  auto uni = boolean2(ops, BooleanOp::union_);
  REQUIRE(uni.has_value());
  REQUIRE(uni.value().polygons.polygons.size() == 1);
  const Polygon2d& p = uni.value().polygons.polygons[0];
  CHECK(p.ring_count() == 1);
  CHECK(p.ring(0).size() == 4);  // collinear shared-edge endpoints merged away
  CHECK(area_of(uni.value().polygons) == doctest::Approx(2.0));
}

TEST_CASE("boolean2: fill rules and self-intersection") {
  Polygon2d bowtie;
  bowtie.vertices = {{0, 0}, {2, 2}, {2, 0}, {0, 2}};
  bowtie.ring_starts = {0, 4};
  const Polygon2d ops[1] = {bowtie};
  auto ok = boolean2(ops, BooleanOp::union_);
  REQUIRE(ok.has_value());
  CHECK(std::isfinite(area_of(ok.value().polygons)));

  BooleanOptions strict;
  strict.allow_self_intersections = false;
  auto bad = boolean2(ops, BooleanOp::union_, strict);
  REQUIRE(!bad.has_value());
  CHECK(bad.error().code == Errc::self_intersecting);

  const Polygon2d a = rect(0, 0, 1, 1);
  const Polygon2d b = rect(0, 0, 1, 1);
  const Polygon2d same[2] = {a, b};
  BooleanOptions even;
  even.fill_rule = FillRule::even_odd;
  auto x = boolean2(same, BooleanOp::union_, even);
  REQUIRE(x.has_value());
  CHECK(x.value().polygons.polygons.empty());
}

TEST_CASE("boolean2: rejects unsupported input") {
  Polygon2d nan = rect(0, 0, 1, 1);
  nan.vertices[0][0] = std::numeric_limits<double>::quiet_NaN();
  const Polygon2d bad[1] = {nan};
  auto r = boolean2(bad, BooleanOp::union_);
  REQUIRE(!r.has_value());
  CHECK(r.error().code == Errc::non_finite);

  Polygon2d huge = rect(0, 0, 1e200, 1);
  huge.vertices[1][1] = 1e-200;
  const Polygon2d wide[1] = {huge};
  auto w = boolean2(wide, BooleanOp::union_);
  REQUIRE(!w.has_value());
  CHECK(w.error().code == Errc::unsupported);

  std::vector<Polygon2d> none;
  auto e = boolean2(none, BooleanOp::union_);
  REQUIRE(!e.has_value());
  CHECK(e.error().code == Errc::degenerate);
}
```

- [x] **Step 2: 运行确认失败**

```powershell
cmake --build build --config Debug
```

Expected: 找不到 `boolean2.hpp`；`Errc::self_intersecting` 不存在。

- [x] **Step 3: 实现**

1. `error.hpp` 在 `Errc` 末尾加 `self_intersecting`，`to_string` 加 `"self_intersecting"`；`coverage_gaps_test.cpp` 补一条 `to_string` 断言。
2. 校验与跨度：空操作数 → `degenerate`；每个 `validate_polygon2`；跨度复用 `delaunay2` 的 `kMaxSpan=1314` 逻辑（`frexp` 扫描）→ 超限 `unsupported`。
3. 归一化：外环 CCW、洞 CW（反转顶点序列，`ring_starts` 重建）。每段生成 `ArrangementSegment{ p_i, p_{i+1}, source=操作数索引 }`。
4. `Arrangement2::build`；若 `!allow_self_intersections`，遍历 `input_pairs()`，把原段索引映射到 `(ring, position)`，同环且非相邻（含首尾相邻）→ `self_intersecting`。
5. 面分类：对每个 bounded 面，`w` = 各操作数 winding；`accept(w)` 按 fill rule（`non_zero: w!=0`、`even_odd: w%2!=0`、`positive: w>0`）；`union_`: 任一 accept；`intersection`: 全部 accept；`difference`: accept(0) 且其余全部 !accept；`symmetric_difference`: accept 的操作数个数为奇数。
6. 边界提取：对每条半边，`kept[face] != kept[twin.face]` 即边界；从边界半边出发，沿 `next` 前进并在遇到内部边时跨到 twin 继续（直到下一条边界边），收集闭环；若某环退化（<3 点）丢弃。合并共线顶点：`sign_orient2d(prev, cur, next) == 0` 时删除 `cur` 并计数。
7. 定向：环有向面积 >0 为外环，<0 为洞；洞挂到包含其一个顶点的最小 |面积| 外环（精确射线测试）。每个环旋转到字典序最小顶点开头；外环按首顶点排序，洞按首顶点排序。
8. 输出 `Polygon2d`（`vertices` 为各环拼接，`ring_starts` 累加）；`report.input_edges` = 原始段数，`output_polygons`，`merged_vertices`，`max_residual = arrangement.max_residual()`。
9. 空结果合法（例如不相交的交集）。

- [x] **Step 4: 运行确认通过**

```powershell
cmake --build build --config Debug
ctest --test-dir build -C Debug --output-on-failure
```

Expected: 6 个 `boolean2` 用例全绿；面积误差在 `doctest::Approx` 默认容差内。

- [x] **Step 5: 提交**

```powershell
git add include/geom/cgeom/boolean2.hpp include/geom/core/error.hpp tests/unit/boolean2_test.cpp tests/unit/coverage_gaps_test.cpp tests/CMakeLists.txt
git commit -m "feat(m1b): add arrangement-based 2D polygon boolean operations"
```

---

## Task 7: 约束 Delaunay（CDT）

**Files:**
- Create: `include/geom/cgeom/cdt2.hpp`
- Test: `tests/unit/cdt2_test.cpp`
- Modify: `tests/CMakeLists.txt`

**Interfaces:**
- Consumes: `delaunay2`、`sign_orient2d`、`sign_incircle`、`Result<T>`、`Errc`。
- Produces:

```cpp
struct CdtConstraint {
  std::uint32_t a = 0;
  std::uint32_t b = 0;
};

struct CdtResult {
  std::vector<std::array<std::uint32_t, 3>> triangles;  // CCW，索引到 points
};

Result<CdtResult> cdt2(std::span<const Point2d> points,
                       std::span<const CdtConstraint> constraints);
```

契约：非约束点需一般位置（同 `delaunay2`）；约束边可位于凸包内或边界上，端点必须是输入点；约束之间只允许在共享端点相接；重复/反向约束去重。错误：`non_finite`、`index_out_of_range`、`degenerate`（点数不足/全共线/`a==b`）、`self_intersecting`（约束相交或与非相邻约束共线重叠）。

- [x] **Step 1: 写失败测试**

`tests/unit/cdt2_test.cpp`：

```cpp
#include <doctest/doctest.h>

#include <geom/cgeom/cdt2.hpp>
#include <geom/numeric/predicates2.hpp>

#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <vector>

using geom::v1::CdtConstraint;
using geom::v1::cdt2;
using geom::v1::Errc;
using geom::v1::Point2d;

namespace {

bool has_edge(const std::vector<std::array<std::uint32_t, 3>>& tris,
              std::uint32_t a, std::uint32_t b) {
  for (const auto& t : tris) {
    for (int e = 0; e < 3; ++e) {
      const std::uint32_t u = t[e], v = t[(e + 1) % 3];
      if ((u == a && v == b) || (u == b && v == a)) return true;
    }
  }
  return false;
}

int edge_use(const std::vector<std::array<std::uint32_t, 3>>& tris,
             std::uint32_t a, std::uint32_t b) {
  int n = 0;
  for (const auto& t : tris) {
    for (int e = 0; e < 3; ++e) {
      const std::uint32_t u = t[e], v = t[(e + 1) % 3];
      if ((u == a && v == b) || (u == b && v == a)) ++n;
    }
  }
  return n;
}

}  // namespace

TEST_CASE("cdt2: recovers either diagonal of a convex quadrilateral") {
  const std::vector<Point2d> pts{{0, 0}, {4, 0}, {3, 2}, {0, 3}};
  const CdtConstraint c02[1] = {{0, 2}};
  auto r02 = cdt2(pts, c02);
  REQUIRE(r02.has_value());
  CHECK(has_edge(r02.value().triangles, 0, 2));
  CHECK(edge_use(r02.value().triangles, 0, 2) == 2);

  const CdtConstraint c13[1] = {{1, 3}};
  auto r13 = cdt2(pts, c13);
  REQUIRE(r13.has_value());
  CHECK(has_edge(r13.value().triangles, 1, 3));
  CHECK(edge_use(r13.value().triangles, 1, 3) == 2);
}

TEST_CASE("cdt2: constrained edge on the hull and duplicate constraints") {
  const std::vector<Point2d> pts{{0, 0}, {4, 0}, {3, 2}, {0, 3}};
  const CdtConstraint cs[3] = {{0, 1}, {1, 0}, {0, 1}};
  auto r = cdt2(pts, cs);
  REQUIRE(r.has_value());
  CHECK(edge_use(r.value().triangles, 0, 1) == 1);  // hull edge
}

TEST_CASE("cdt2: rejects crossing constraints and bad indices") {
  const std::vector<Point2d> pts{{0, 0}, {4, 0}, {3, 2}, {0, 3}};
  const CdtConstraint cross[2] = {{0, 2}, {1, 3}};
  auto r = cdt2(pts, cross);
  REQUIRE(!r.has_value());
  CHECK(r.error().code == Errc::self_intersecting);

  const CdtConstraint bad[1] = {{0, 9}};
  auto b = cdt2(pts, bad);
  REQUIRE(!b.has_value());
  CHECK(b.error().code == Errc::index_out_of_range);

  const CdtConstraint degenerate[1] = {{2, 2}};
  auto d = cdt2(pts, degenerate);
  REQUIRE(!d.has_value());
  CHECK(d.error().code == Errc::degenerate);
}

TEST_CASE("cdt2: non-constrained edges stay Delaunay") {
  std::vector<Point2d> pts;
  for (int y = 0; y < 6; ++y) {
    for (int x = 0; x < 6; ++x) {
      pts.push_back(Point2d{static_cast<double>(x) + 0.013 * y,
                            static_cast<double>(y) + 0.017 * x});
    }
  }
  const CdtConstraint cs[2] = {{0, 7}, {2, 9}};
  auto r = cdt2(pts, cs);
  REQUIRE(r.has_value());
  std::vector<bool> constrained(pts.size() * pts.size(), false);
  constrained[0 * pts.size() + 7] = constrained[7 * pts.size() + 0] = true;
  constrained[2 * pts.size() + 9] = constrained[9 * pts.size() + 2] = true;
  for (const auto& t : r.value().triangles) {
    for (int e = 0; e < 3; ++e) {
      const std::uint32_t a = t[e], b = t[(e + 1) % 3], c = t[(e + 2) % 3];
      if (constrained[a * pts.size() + b]) continue;
      // Opposite vertex across edge (a,b): find the other triangle.
      std::uint32_t d = 0xFFFFFFFFu;
      for (const auto& u : r.value().triangles) {
        for (int k = 0; k < 3; ++k) {
          const std::uint32_t p = u[k], q = u[(k + 1) % 3];
          if ((p == b && q == a) || (p == a && q == b)) d = u[(k + 2) % 3];
        }
      }
      if (d == 0xFFFFFFFFu) continue;  // hull edge
      CHECK(geom::v1::sign_incircle(pts[a][0], pts[a][1], pts[b][0], pts[b][1],
                                    pts[c][0], pts[c][1], pts[d][0], pts[d][1]) <= 0);
    }
  }
}
```

- [x] **Step 2: 运行确认失败**

```powershell
cmake --build build --config Debug
```

Expected: 找不到 `cdt2.hpp`。

- [x] **Step 3: 实现 `cdt2.hpp`**

1. 校验：坐标有限、约束索引 < n、`a != b`、跨度 ≤1314；每对约束用精确谓词检查是否只在共享端点相接，否则 `self_intersecting`；去重（无向边）。
2. 基三角化：`delaunay2(points)`；空且 n≥3 → `degenerate`。
3. 建立三角形数组与边→三角形邻接；标记约束边集合。
4. 约束恢复（Sloan）：若约束边已存在，标记；否则从 `u` 出发找出所有被 `uv` 穿过的三角形（精确 orient 判定），删除这些三角形得到腔体边界链；把腔体按 `uv` 分成两个简单多边形，用耳切（可复用 `triangulate2` 的耳切核心，或内联同样的精确判凸/包含判定）重新三角化；`uv` 标记为约束。
5. Delaunay 恢复：对非约束内部边做 Lawson 翻转，`sign_incircle > 0` 时翻，直到无翻转或达到迭代上限（`4 * triangles.size()`，超出返回 `numeric_failure`）。
6. 输出：每个三角形用 `sign_orient2d` 保证 CCW，按 `(a,b,c)` 排序；`ok` 由 `Result` 表达。

- [x] **Step 4: 运行确认通过**

```powershell
cmake --build build --config Debug
ctest --test-dir build -C Debug --output-on-failure
```

Expected: 4 个 `cdt2` 用例全绿；两条对角线的约束各自被恢复。

- [x] **Step 5: 提交**

```powershell
git add include/geom/cgeom/cdt2.hpp tests/unit/cdt2_test.cpp tests/CMakeLists.txt
git commit -m "feat(m1b): constrained Delaunay triangulation with edge recovery"
```

---

## Task 8: QA-1/QA-2/QA-4 与出口文档

**Files:**
- Create: `bench/bench_m1b.cpp`, `tests/unit/qa_concurrency_test.cpp`, `docs/benchmarks/m1b.md`, `docs/adr/0020-m1b-arrangement-and-triangulation.md`
- Modify: `CMakeLists.txt`, `tests/CMakeLists.txt`, `docs/architecture/geometry-library-architecture.md`（§16.3）

**Interfaces:**
- Consumes: `KdTree`、`Bvh`、`boolean2`、`polygon_area2`。
- Produces: `geom_bench_m1b` 可执行文件；QA-2 并发一致性测试；M1b 状态文档。

- [ ] **Step 1: 写 QA-2 测试与基准骨架**

`tests/unit/qa_concurrency_test.cpp`：

```cpp
#include <doctest/doctest.h>

#include <geom/index/bvh.hpp>

#include <array>
#include <cmath>
#include <cstdint>
#include <thread>
#include <vector>

using geom::v1::Bvh;
using geom::v1::Point3d;
using geom::v1::Ray3d;
using geom::v1::Triangle3d;

TEST_CASE("qa2: 4 threads x 1e6 rays agree with serial") {
  std::vector<Triangle3d> tris;
  for (int y = 0; y < 100; ++y) {
    for (int x = 0; x < 100; ++x) {
      const double fx = static_cast<double>(x), fy = static_cast<double>(y);
      tris.push_back(Triangle3d{Point3d{fx, fy, 0}, Point3d{fx + 0.9, fy, 0},
                                Point3d{fx, fy + 0.9, 0}});
    }
  }
  Bvh bvh;
  bvh.build(tris);

  constexpr std::uint32_t kRays = 1000000;
  auto hash_of = [&](std::uint32_t seed) {
    std::uint64_t h = 1469598103934665603ull;
    for (std::uint32_t i = 0; i < kRays; ++i) {
      seed = seed * 1664525u + 1013904223u;
      const double ox = (seed % 10000) * 0.01;
      seed = seed * 1664525u + 1013904223u;
      const double oy = (seed % 10000) * 0.01;
      const Ray3d ray{Point3d{ox, oy, -1}, geom::v1::Vec3d{0, 0, 1}};
      Bvh::Hit hit;
      const bool ok = bvh.raycast(ray, 100.0, hit);
      h ^= ok ? (static_cast<std::uint64_t>(hit.prim) + 1) : 0;
      h *= 1099511628211ull;
    }
    return h;
  };

  std::array<std::uint64_t, 4> parallel{};
  std::array<std::thread, 4> workers;
  for (std::uint32_t t = 0; t < 4; ++t) {
    workers[t] = std::thread([&, t] { parallel[t] = hash_of(1000 + t); });
  }
  for (auto& w : workers) w.join();
  for (std::uint32_t t = 0; t < 4; ++t) CHECK(parallel[t] == hash_of(1000 + t));
}
```

基准 `bench/bench_m1b.cpp`（QA-1 + QA-4，Release 下超预算返回 2）：

```cpp
#include <geom/cgeom/boolean2.hpp>
#include <geom/cgeom/polygon2.hpp>
#include <geom/index/bvh.hpp>

#include <chrono>
#include <cmath>
#include <cstdio>
#include <random>
#include <vector>

using namespace geom::v1;

namespace {

double ms_since(std::chrono::steady_clock::time_point t0) {
  return std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - t0)
      .count();
}

Polygon2d star(std::uint32_t n, double r0, double r1) {
  Polygon2d p;
  p.vertices.resize(n);
  for (std::uint32_t i = 0; i < n; ++i) {
    const double a = 2.0 * 3.14159265358979323846 * i / n;
    const double r = (i % 2 == 0) ? r0 : r1;
    p.vertices[i] = Point2d{std::cos(a) * r, std::sin(a) * r};
  }
  p.ring_starts = {0, n};
  return p;
}

}  // namespace

int main() {
  // QA-1: two 100k-vertex polygons, area identity |A| + |B| = |A u B| + |A n B|.
  Polygon2d a = star(50000, 100.0, 99.0);
  Polygon2d b = a;
  for (auto& q : b.vertices) q[0] += 0.35;
  std::vector<Polygon2d> ops{a, b};

  const auto t0 = std::chrono::steady_clock::now();
  auto uni = boolean2(ops, BooleanOp::union_);
  const double union_ms = ms_since(t0);
  auto inter = boolean2(ops, BooleanOp::intersection);
  if (!uni || !inter) {
    std::printf("QA-1 failed: boolean error\n");
    return 2;
  }
  const double area_a = polygon_area2(a);
  const double area_b = polygon_area2(b);
  const double lhs = area_a + area_b;
  const double rhs = polygon_set_area2(uni.value().polygons) +
                     polygon_set_area2(inter.value().polygons);
  const double rel = std::fabs(lhs - rhs) / std::fabs(lhs);
  std::printf("QA-1 union %.2f ms, area rel err %.3e (limit 1e-12, 80 ms)\n", union_ms,
              rel);

  // QA-4: 1M triangles, build and 8M rays.
  std::vector<Triangle3d> tris;
  tris.reserve(1000000);
  for (int y = 0; y < 500; ++y) {
    for (int x = 0; x < 1000; ++x) {
      const double fx = x, fy = y;
      tris.push_back(Triangle3d{Point3d{fx, fy, 0}, Point3d{fx + 1, fy, 0},
                                Point3d{fx, fy + 1, 0}});
      tris.push_back(Triangle3d{Point3d{fx + 1, fy, 0}, Point3d{fx + 1, fy + 1, 0},
                                Point3d{fx, fy + 1, 0}});
    }
  }
  Bvh bvh;
  const auto tb = std::chrono::steady_clock::now();
  bvh.build(tris);
  const double build_ms = ms_since(tb);

  std::mt19937 rng(20261107);
  std::uniform_real_distribution<double> uni_x(0.0, 1000.0);
  std::uniform_real_distribution<double> uni_y(0.0, 500.0);
  constexpr std::uint32_t kRays = 8000000;
  std::uint64_t hits = 0;
  const auto tr = std::chrono::steady_clock::now();
  for (std::uint32_t i = 0; i < kRays; ++i) {
    const Ray3d ray{Point3d{uni_x(rng), uni_y(rng), -1}, Vec3d{0, 0, 1}};
    Bvh::Hit hit;
    if (bvh.raycast(ray, 10.0, hit)) ++hits;
  }
  const double ray_s = ms_since(tr);
  const double mrays = kRays / (ray_s / 1000.0) / 1e6;
  std::printf("QA-4 build %.1f ms (limit 350), %.2f Mray/s (limit 8), hits %llu\n",
              build_ms, mrays, static_cast<unsigned long long>(hits));

  bool ok = rel <= 1e-12 && union_ms <= 80.0 && build_ms <= 350.0 && mrays >= 8.0;
#ifndef NDEBUG
  ok = true;  // Debug builds report only.
#endif
  return ok ? 0 : 2;
}
```

- [ ] **Step 2: 注册构建与测试**

`CMakeLists.txt`：

```cmake
option(GEOM_TSAN "Build with ThreadSanitizer (GCC/Clang)" OFF)
if(GEOM_TSAN)
  if(MSVC)
    message(FATAL_ERROR "GEOM_TSAN requires GCC or Clang")
  endif()
  target_compile_options(geom_warnings INTERFACE -fsanitize=thread -g)
  target_link_options(geom_warnings INTERFACE -fsanitize=thread)
endif()
```

`GEOM_BUILD_BENCH` 分支追加：

```cmake
add_executable(geom_bench_m1b bench/bench_m1b.cpp)
target_link_libraries(geom_bench_m1b PRIVATE geom_core)
```

`tests/CMakeLists.txt` 的 `GEOM_TEST_SOURCES` 追加 `unit/qa_concurrency_test.cpp`。

- [ ] **Step 3: 运行**

```powershell
cmake -S . -B build -DGEOM_WERROR=ON -DGEOM_BUILD_BENCH=ON
cmake --build build --config Release
ctest --test-dir build -C Release --output-on-failure
.\build\Release\geom_bench_m1b.exe
```

Expected: ctest 全绿（含 QA-2 一致性）；基准退出码 0，QA-1 ≤80 ms 且面积相对误差 ≤1e-12，QA-4 构建 ≤350 ms、≥8 Mray/s。若任一超预算，记录实测值并在本任务报告中标记 FAIL，不得修改门禁。

- [ ] **Step 4: 写文档**

- `docs/benchmarks/m1b.md`：日期、OS、CPU、编译器、`QA-1` 顶点数/耗时/面积误差、`QA-4` 三角形数/构建/吞吐、`ctest` 用例数与断言数、未达标项（如有）。
- `docs/adr/0020-m1b-arrangement-and-triangulation.md`：记录本计划「与架构文档的偏差」三条、选择理由、QA-1/QA-4 判据与升级触发条件。
- 架构 §16.3：M1b 改为「代码与测试完成 / 验收状态」两段；QA-1/2/4 未全达标前不得写「已满足」；M1a 覆盖率若首次 CI 已有数据则回填实际百分比并把 M1a 改为验收完成，否则保留待回填说明。

- [ ] **Step 5: 提交**

```powershell
git add bench/bench_m1b.cpp tests/unit/qa_concurrency_test.cpp CMakeLists.txt tests/CMakeLists.txt docs/benchmarks/m1b.md docs/adr/0020-m1b-arrangement-and-triangulation.md docs/architecture/geometry-library-architecture.md
git commit -m "test(m1b): add QA-1/QA-2/QA-4 gates and M1b exit evidence"
```

---

## M1b 出口验收（不新增代码）

在 Task 8 之后执行，全部满足才可声称 M1b 完成：

- [x] `cmake -S . -B build -DGEOM_WERROR=ON -DGEOM_BUILD_BENCH=ON`
- [x] `cmake --build build --config Debug` 与 `--config Release` 均 0 警告（`GEOM_WERROR=ON`）
- [x] `ctest --test-dir build -C Debug --output-on-failure` 与 `-C Release` 全绿
- [ ] `geom_bench_m1b`：QA-1 并集 ≤80 ms、面积相对误差 ≤1e-12；QA-4 构建 ≤350 ms、≥8 Mray/s/核（2026-11-07 实测 **FAIL**：winding 播种修复后 QA-1 面积误差 1.783e-13 已通过，但并集 4.8–4.9 s；QA-4 构建 2.6–2.7 s、0.42–0.44 Mray/s；见 `docs/benchmarks/m1b.md`）
- [ ] QA-2：`GEOM_TSAN=ON` 的 Linux/Clang 构建下 `ctest -R qa_concurrency` 无 TSan 报告（本地 MSVC 只跑一致性断言，通过；TSan 待 CI）
- [ ] 覆盖率：`GEOM_COVERAGE=ON` + gcovr，`cgeom`/`index` 行 ≥90%、分支 ≥85%；diff 覆盖 100%（本地 MSVC 无 gcovr，待首次 CI `coverage` job 回填）
- [x] `docs/benchmarks/m1b.md` 记录实测环境与数值
- [x] 架构 §16.3 的 M1b 状态与实际一致（未达标项明确写「未声称」）

---

## 自检（Self-Review）

**Spec coverage**

| 架构/策略要求 | 覆盖任务 |
| --- | --- |
| §8.2 简单多边形三角剖分（含洞） | T3 |
| §8.2 CDT（约束插入 + 腔体重新三角化） | T7 |
| §8.3 2D 布尔（arrangement + winding + fill rule + 输出规则化 + BooleanReport） | T4/T5/T6 |
| §8.4 arrangement（DCEL、面/边/顶点、点定位） | T5 |
| §16.3 KD 树 | T2 |
| §8.3 自交输入默认允许 / `allow_self_intersections=false` 报错 | T6 |
| QA-1（10 万顶点布尔 ≤80 ms、面积 ≤1e-12） | T8 |
| QA-2（4 线程、TSan、结果一致） | T8 |
| QA-4（100 万三角形构建 ≤350 ms、≥8 Mray/s） | T8 |
| 测试策略 §2 覆盖门槛、Errc 用例、API 三类用例 | 全局约束 + T6（新 Errc）+ T8 |
| 测试策略退化语料：共线/重合边、垂直边、T 接点、重复点 | T1/T3/T4/T6 |

**占位符扫描**：无 `TODO`/`TBD`；所有测试期望值为手算结果；T4/T5/T6 的计数已按欧拉公式核算（单连通 V−E+F=2）。

**已知风险**：QA-1 的 80 ms 与 QA-4 的吞吐是最可能不达标的两项；若失败，按 ADR-0020 的升级触发条件另立「扫描线/性能」计划，不修改门禁、不声称 M1b 完成。

**类型一致性**：`Polygon2d::ring_starts`（`uint32_t`）、`Triangulation2::triangles`、`ArrangementSegment{Point2d a,b; uint32_t source}`、`Arrangement2::HalfEdge/Face`、`BooleanOp`/`FillRule`/`BooleanOptions`、`CdtConstraint` 在测试与后续任务中名称与顺序一致；`detail::split_segments` 只在 T4/T5 使用。
