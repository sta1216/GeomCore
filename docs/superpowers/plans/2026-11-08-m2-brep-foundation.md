# M2 B-rep 基础实施计划

> **For agentic workers:** REQUIRED SUB-SKILL: Use subagent-driven-development (recommended) or executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** 建立可编辑、可校验、可网格化的 B-rep 基础——扁平拓扑存储、Euler 算子、解析基本体、显式容差链、水密网格化与质量属性——并用 QA-10/QA-11 验收。

**Architecture:** 拓扑是扁平数组 + 32 位代际 handle（§6.4/§6.8），所有变更经 `SolidEditor` 暂存并在提交前校验，失败回滚；几何保存在 `GeometryStore` 的 arena 中，实体只存句柄与参数域；解析曲面/曲线保留类型身份（ADR-016）。网格化是 B-rep 与网格世界的唯一契约（ADR-017）：边采样一次、相邻面复用同一批采样点，输出水密三角网。质量属性用散度定理把体积/质心/惯性矩化为边界上的面积分：平面面用精确多边形公式，解析曲面面用 Gauss–Legendre 求积。

**Tech Stack:** C++20、CMake 4.x、MSVC `/W4 /WX`、doctest v2.4.11、既有 `geom::v1` 的 primitives/numeric/cgeom/index/curves/surfaces。

**Spec:** `docs/architecture/geometry-library-architecture.md`（v2.1）——§4.6 容忍建模、§5.7 曲线曲面类型、§5.8 B-rep 拓扑实体、§6.4/§6.8 存储与事务、§8.7 Euler 算子与基本体、§8.10 网格化、§8.11 质量属性、§12.3 不变量校验器、§16.2 M2 出口；测试策略见 `docs/quality/testing-strategy.md`。

## Global Constraints

- 命名空间 `geom::v1`；公开头文件 `#pragma once`；不得在 `include/` 中出现 doctest。
- C++20；MSVC 19.44 `/W4 /WX /utf-8 /fp:precise`；GCC/Clang `-Wall -Wextra -Werror`；新增头文件必须能单独编译通过。
- 无全局 epsilon：所有容差来自实体字段或调用参数；禁止用放宽容差「修复」失败。
- 拓扑只能经 `SolidEditor` 修改；提交前必须通过拓扑与容差校验，失败返回 `Result` 并保持调用前状态。
- 新增 `Errc` 枚举值必须同时更新 `to_string` 与 `tests/unit/coverage_gaps_test.cpp` 的字符串断言。
- 每个测试文件都要加入 `tests/CMakeLists.txt` 的 `GEOM_TEST_SOURCES`。
- 断言禁止 `CHECK(x >= 0)` 式弱断言；数值比较写明确容差。
- 每个任务先写失败测试并运行确认失败（RED 证据），再实现，Debug 与 Release 都通过，最后单 commit（`feat(m2): ...`）。
- 严禁 `throw` 穿越公开边界；公开 API 可失败处返回 `Result<T>`；实体容量（2^24）耗尽可能时返回 `Errc::unsupported` 并附 message，而不是抛异常。
- 不改变既有公开 API 的语义；本计划只做新增。

## 范围与偏差（需评审）

1. **Euler 算子收敛为 MEV/KEV**：M2 只实现「分裂边 / 合并度为 2 的顶点」；MEF/KEF 需要面分割与 pcurve 重建，依赖 M3 的 SSI 与面参数域布尔，推迟到 M3。
2. **pcurve 只用线性段**：M2 的基本体所有 face 边界在参数域内都是直线段（平面多边形、圆柱/圆锥的 u=const 或 v=const 边、球/环面的网格边界），因此 `Pcurve2d` 定义为 `{p0, p1}`；一般曲线 pcurve 留给 M3。
3. **网格化只覆盖解析曲面**：平面、圆柱、圆锥、球、环面按参数域规则采样；NURBS 面的自适应细分推迟到 M3。
4. **QA-11 的参考值用解析闭式而非 OCCT**：CI 不引入 OCCT 依赖；对 box/cylinder/sphere/torus 用闭式体积、质心、惯性矩做 oracle，并在 `docs/benchmarks/m2.md` 记录该替代。
5. **不做序列化、属性系统、Solid 布尔**：这些分别属于 §6.7、§6.8 属性表与 §8.8，留给后续里程碑。

## Review Focus

1. **面朝向与 coedge 顺序**：外法向 = `reversed` 取反后的 `du × dv`；内腔 shell 必须整体反向。测试须覆盖外壳与内腔两种朝向，而不是只看体积符号。
2. **容差链单调性**：`edge.tol ≥ max(两端 vertex.tol)`、`face.tol ≥ 相邻 edge.tol`、`coedge.tol ≤ face.tol`、`solid.tol ≥ max(子实体 tol)`；违反时返回 `Errc::tolerance_violation` 而不是静默通过。
3. **接缝与闭合曲面**：圆柱侧面/圆锥/球/环面的接缝边在同一条 loop 里出现两次（正反各一），coedge 配对校验不能把它当成「边出现两次」而误报。
4. **退化极点**：球的两极由接缝边端点表示；网格化必须在极点合并为单个顶点且不产生零面积三角形。
5. **网格化边界复用**：相邻面共享边的采样点必须逐位相同（同一 `edge_points` 来源），反向后顺序相反；任何一面自行重采样都会破坏水密性。

---

## File Map

| 文件 | 责任 |
| --- | --- |
| `include/geom/curves/analytic_curve.hpp` | `LineCurve3`、`CircleCurve3`：求值、导数、包围盒、合法性 |
| `include/geom/surfaces/analytic.hpp` | `Cylinder3`、`Cone3`、`Sphere3`、`Torus3`：求值、偏导、单位法向、合法性 |
| `include/geom/brep/topology.hpp` | id 类型、实体、`GeometryStore`、`Pcurve2d`、`Brep` 只读访问 |
| `include/geom/brep/solid_editor.hpp` | `SolidEditor`：创建/删除、暂存、`commit`/`rollback`、MEV/KEV |
| `include/geom/brep/validate.hpp` | `ValidationIssue`、`ValidationReport`、拓扑与容差校验器 |
| `include/geom/brep/primitives.hpp` | `make_box`、`make_cylinder`、`make_cone`、`make_sphere`、`make_torus` |
| `include/geom/brep/tessellate.hpp` | `TriMesh`、`TessellationOptions`、`tessellate` |
| `include/geom/brep/mass_properties.hpp` | `MassProperties`、`mass_properties` |
| `bench/bench_m2.cpp` | QA-10 网格化吞吐、QA-11 质量属性误差 |
| `docs/benchmarks/m2.md` | 实测数据与 FAIL/PASS 记录 |
| `docs/adr/0021-brep-storage-and-euler-operators.md` | 扁平存储、handle、事务、MEV/KEV 决策 |
| `docs/adr/0022-tessellation-and-mass-properties.md` | 采样契约、偏差控制、积分方法 |

---

### Task 1: 解析曲线（直线与圆）

**Files:**
- Create: `include/geom/curves/analytic_curve.hpp`
- Test: `tests/unit/analytic_curve_test.cpp`
- Modify: `tests/CMakeLists.txt`（在 `GEOM_TEST_SOURCES` 中加一行 `unit/analytic_curve_test.cpp`）

**Interfaces:**
- Consumes: `Point3d`、`Vec3d`、`Box3d`（`geom/primitives/{point,vec,box}.hpp`）；`dot`、`cross`。
- Produces: `LineCurve3`、`CircleCurve3`，签名如下：

```cpp
struct LineCurve3 {
  Point3d origin{};
  Vec3d direction{1.0, 0.0, 0.0};
  bool is_valid() const noexcept;
  Point3d eval(double t) const noexcept;        // t 为弧长
  Vec3d derivative(double) const noexcept;      // 恒为 direction
  Box3d bounds(double t0, double t1) const noexcept;
};

struct CircleCurve3 {
  Point3d center{};
  Vec3d x_axis{1.0, 0.0, 0.0};
  Vec3d y_axis{0.0, 1.0, 0.0};
  double radius = 1.0;
  bool is_valid() const noexcept;
  Point3d eval(double angle) const noexcept;
  Vec3d derivative(double angle) const noexcept;
  Box3d bounds(double a0, double a1) const noexcept;
};
```

- [ ] **Step 1: 写失败测试**

```cpp
TEST_CASE("line curve evaluates and bounds") {
  const geom::LineCurve3 line{{1, 2, 3}, {0, 1, 0}};
  CHECK(line.is_valid());
  const geom::Point3d p = line.eval(4.0);
  CHECK(p[0] == doctest::Approx(1.0));
  CHECK(p[1] == doctest::Approx(6.0));
  CHECK(p[2] == doctest::Approx(3.0));
  CHECK(line.derivative(4.0)[1] == doctest::Approx(1.0));
  const geom::Box3d b = line.bounds(0.0, 4.0);
  CHECK(b.lo[0] == doctest::Approx(1.0));
  CHECK(b.lo[1] == doctest::Approx(2.0));
  CHECK(b.hi[0] == doctest::Approx(1.0));
  CHECK(b.hi[1] == doctest::Approx(6.0));
}

TEST_CASE("circle curve quarter arc bounds are exact") {
  const geom::CircleCurve3 c{{0, 0, 0}, {1, 0, 0}, {0, 1, 0}, 2.0};
  CHECK(c.is_valid());
  const geom::Point3d p0 = c.eval(0.0);
  CHECK(p0[0] == doctest::Approx(2.0));
  CHECK(p0[1] == doctest::Approx(0.0));
  const geom::Box3d b = c.bounds(0.0, 3.14159265358979323846 / 2);
  CHECK(b.lo[0] == doctest::Approx(0.0).epsilon(1e-12));
  CHECK(b.hi[0] == doctest::Approx(2.0).epsilon(1e-12));
  CHECK(b.lo[1] == doctest::Approx(0.0).epsilon(1e-12));
  CHECK(b.hi[1] == doctest::Approx(2.0).epsilon(1e-12));
}

TEST_CASE("analytic curves reject non-finite and unnormalized data") {
  geom::LineCurve3 line{{0, 0, 0}, {0, 2, 0}};
  CHECK_FALSE(line.is_valid());
  line.direction = {0, 0, 0};
  CHECK_FALSE(line.is_valid());
  line.direction = {std::numeric_limits<double>::quiet_NaN(), 0, 0};
  CHECK_FALSE(line.is_valid());
  geom::CircleCurve3 c{};
  c.radius = std::numeric_limits<double>::infinity();
  CHECK_FALSE(c.is_valid());
  c = geom::CircleCurve3{};
  c.y_axis = {1, 0, 0};  // 与 x_axis 重合
  CHECK_FALSE(c.is_valid());
}
```

- [ ] **Step 2: 运行确认失败**

Run: `cmake --build build --config Debug --target geom_tests` → 期望 `Cannot open include file: 'geom/curves/analytic_curve.hpp'`。

- [ ] **Step 3: 实现**

`LineCurve3::is_valid`：坐标/方向全部有限，且 `|1 - dot(dir,dir)| <= 1e-12`。
`CircleCurve3::is_valid`：全部有限，`radius > 0`，两轴单位且正交（`|dot| <= 1e-12`，`| |x|-1 | <= 1e-12`、`| |y|-1 | <= 1e-12`）。
`CircleCurve3::bounds(a0,a1)`：每个轴在区间内的极值由 `center + r(cos t x + sin t y)` 的解析驻点确定——候选角为 `a0`、`a1` 以及各轴分量导数零点（`atan2` 解出的 `t` 加上 `kπ`），落在区间内的才参与；不要用采样近似。前置条件 `std::isfinite(a0) && std::isfinite(a1) && a0 <= a1`，debug 下 `assert`。
`LineCurve3::bounds(t0,t1)`：两端点构成的盒（`t0 <= t1` 前置）。

- [ ] **Step 4: 运行确认通过**

Run: `cmake --build build --config Debug --target geom_tests; .\build\Debug\geom_tests.exe --test-case="*analytic curve*"` → 期望全过；再跑 Release 与全套 `ctest`。

- [ ] **Step 5: 提交**

```bash
git add include/geom/curves/analytic_curve.hpp tests/unit/analytic_curve_test.cpp tests/CMakeLists.txt
git commit -m "feat(m2): add analytic line and circle curves"
```

---

### Task 2: 解析曲面（圆柱、圆锥、球、环面）

**Files:**
- Create: `include/geom/surfaces/analytic.hpp`
- Test: `tests/unit/analytic_surface_test.cpp`
- Modify: `tests/CMakeLists.txt`

**Interfaces:**
- Consumes: `Point3d`、`Vec3d`、`dot`、`cross`、`std::sqrt`。
- Produces:

```cpp
struct Cylinder3 {
  Point3d origin{};
  Vec3d axis{0.0, 0.0, 1.0};     // 单位
  Vec3d x_axis{1.0, 0.0, 0.0};   // 单位且 ⊥ axis
  double radius = 1.0;
  bool is_valid() const noexcept;
  Point3d eval(double u, double v) const noexcept;   // u 角度，v 沿轴距离
  Vec3d du(double u, double v) const noexcept;
  Vec3d dv(double u, double v) const noexcept;       // = axis
  Vec3d normal(double u, double v) const noexcept;   // 单位外法向
};

struct Cone3 {
  Point3d origin{};              // v = 0 处的圆心
  Vec3d axis{0.0, 0.0, 1.0};
  Vec3d x_axis{1.0, 0.0, 0.0};
  double radius = 1.0;           // v = 0 处半径
  double half_angle = 0.5;       // ∈ (0, π/2)
  bool is_valid() const noexcept;
  Point3d eval(double u, double v) const noexcept;
  Vec3d du(double u, double v) const noexcept;
  Vec3d dv(double u, double v) const noexcept;
  Vec3d normal(double u, double v) const noexcept;
};

struct Sphere3 {
  Point3d center{};
  Vec3d axis{0.0, 0.0, 1.0};
  Vec3d x_axis{1.0, 0.0, 0.0};
  double radius = 1.0;
  bool is_valid() const noexcept;
  Point3d eval(double u, double v) const noexcept;   // v ∈ [-π/2, π/2]
  Vec3d du(double u, double v) const noexcept;
  Vec3d dv(double u, double v) const noexcept;
  Vec3d normal(double u, double v) const noexcept;
};

struct Torus3 {
  Point3d center{};
  Vec3d axis{0.0, 0.0, 1.0};
  Vec3d x_axis{1.0, 0.0, 0.0};
  double major_radius = 2.0;     // R
  double minor_radius = 0.5;     // r，0 < r < R
  bool is_valid() const noexcept;
  Point3d eval(double u, double v) const noexcept;
  Vec3d du(double u, double v) const noexcept;
  Vec3d dv(double u, double v) const noexcept;
  Vec3d normal(double u, double v) const noexcept;
};
```

- [ ] **Step 1: 写失败测试**

```cpp
TEST_CASE("cylinder evaluation, derivatives and outward normal") {
  const geom::Cylinder3 cyl{{0, 0, 0}, {0, 0, 1}, {1, 0, 0}, 2.0};
  CHECK(cyl.is_valid());
  const geom::Point3d p = cyl.eval(0.0, 3.0);
  CHECK(p[0] == doctest::Approx(2.0));
  CHECK(p[1] == doctest::Approx(0.0));
  CHECK(p[2] == doctest::Approx(3.0));
  CHECK(cyl.eval(3.14159265358979323846 / 2, 0.0)[1] == doctest::Approx(2.0).epsilon(1e-12));
  const geom::Vec3d d = cyl.du(0.0, 0.0);
  CHECK(d[0] == doctest::Approx(0.0));
  CHECK(d[1] == doctest::Approx(2.0));
  CHECK(cyl.dv(0.0, 0.0)[2] == doctest::Approx(1.0));
  const geom::Vec3d n0 = cyl.normal(0.0, 0.0);
  CHECK(n0[0] == doctest::Approx(1.0).epsilon(1e-12));
  CHECK(geom::dot(cyl.normal(1.2, 5.0), cyl.normal(1.2, 5.0)) == doctest::Approx(1.0).epsilon(1e-12));
}

TEST_CASE("sphere normal is unit and points outward at poles") {
  const geom::Sphere3 s{{0, 0, 0}, {0, 0, 1}, {1, 0, 0}, 3.0};
  CHECK(s.is_valid());
  const geom::Point3d q = s.eval(0.0, 0.0);
  CHECK(q[0] == doctest::Approx(3.0));
  CHECK(q[1] == doctest::Approx(0.0));
  const geom::Vec3d np = s.normal(0.0, 3.14159265358979323846 / 2);
  CHECK(np[2] == doctest::Approx(1.0).epsilon(1e-12));
  CHECK(geom::dot(np, np) == doctest::Approx(1.0).epsilon(1e-12));
}

TEST_CASE("torus geometry is rejected when r >= R") {
  geom::Torus3 t{};
  t.major_radius = 1.0;
  t.minor_radius = 1.0;
  CHECK_FALSE(t.is_valid());
  t.minor_radius = 0.5;
  CHECK(t.is_valid());
}

TEST_CASE("cone derivative follows the slant") {
  const geom::Cone3 cone{{0, 0, 0}, {0, 0, 1}, {1, 0, 0}, 1.0, 0.25};
  CHECK(cone.is_valid());
  CHECK(cone.eval(0.0, 2.0)[2] == doctest::Approx(2.0).epsilon(1e-12));
  CHECK(cone.eval(0.0, 2.0)[0] ==
        doctest::Approx(1.0 - 2.0 * std::tan(0.25)).epsilon(1e-12));
  CHECK(cone.dv(0.0, 2.0)[0] == doctest::Approx(-std::tan(0.25)).epsilon(1e-12));
}
```

- [ ] **Step 2: 运行确认失败**

Run: `cmake --build build --config Debug --target geom_tests` → 期望找不到 `analytic.hpp`。

- [ ] **Step 3: 实现**

公式：圆柱 `p = origin + v·axis + r(cos u·x + sin u·y)`，`y = cross(axis, x)`；圆锥 `p = origin + v·axis + (r − v·tan α)(cos u·x + sin u·y)`（`v` 增大半径减小，`v = r/tan α` 处收成锥顶）；球 `p = center + r(cos v cos u·x + cos v sin u·y + sin v·axis)`；环面 `p = center + (R + r cos v)(cos u·x + sin u·y) + r sin v·axis`。
法向 = `normalize(cross(du, dv))`（圆锥注意 `dv` 含径向项；球在极点 `cross(du,dv)` 退化，按 `v` 的符号返回 `±axis`，并保证单位长度）。`is_valid` 检查有限性、轴单位与正交、半径/半角范围（球半径 > 0；环面 `0 < r < R`）；退化输入返回 false，不抛异常。

- [ ] **Step 4: 运行确认通过**

Run: `cmake --build build --config Debug --target geom_tests; .\build\Debug\geom_tests.exe --test-case="*analytic surface*"`，再跑 Release 与全套 `ctest`。

- [ ] **Step 5: 提交**

```bash
git add include/geom/surfaces/analytic.hpp tests/unit/analytic_surface_test.cpp tests/CMakeLists.txt
git commit -m "feat(m2): add analytic cylinder, cone, sphere and torus surfaces"
```

### Task 3: 扁平拓扑存储与只读 Brep

**Files:**
- Create: `include/geom/brep/topology.hpp`
- Modify: `include/geom/core/error.hpp`（新增三个 Errc）
- Test: `tests/unit/brep_topology_test.cpp`
- Modify: `tests/CMakeLists.txt`、`tests/unit/coverage_gaps_test.cpp`

**Interfaces:**
- Consumes: `Handle`、`HandleAllocator`、`Point2d/Point3d`、`Plane3`、`LineCurve3`、`CircleCurve3`、`Cylinder3`、`Cone3`、`Sphere3`、`Torus3`、`NurbsCurve<double>`、`NurbsSurface<double>`。
- Produces:

```cpp
namespace geom::v1 {

inline constexpr std::uint32_t kNoGeometry = 0xFFFFFFFFu;

using VertexId = Handle;
using EdgeId = Handle;
using CoedgeId = Handle;
using LoopId = Handle;
using FaceId = Handle;
using ShellId = Handle;
using SolidId = Handle;

struct Vertex { Point3d point{}; double tol = 0.0; };
struct Edge {
  std::uint32_t curve = kNoGeometry;
  double t0 = 0.0;
  double t1 = 1.0;
  VertexId v0{};
  VertexId v1{};
  double tol = 0.0;
};
struct Coedge {
  EdgeId edge{};
  FaceId face{};
  LoopId loop{};
  bool reversed = false;
  std::uint32_t pcurve = kNoGeometry;
  CoedgeId next{};
  CoedgeId prev{};
  double tol = 0.0;
};
struct Loop { FaceId face{}; CoedgeId start{}; LoopId next{}; bool outer = false; };
struct Face { std::uint32_t surface = kNoGeometry; LoopId first_loop{}; bool reversed = false; double tol = 0.0; };
struct Shell { std::uint32_t first_face = 0; std::uint32_t face_count = 0; bool closed = false; };
struct SolidEntity { std::uint32_t first_shell = 0; std::uint32_t shell_count = 0; double tol = 0.0; };

using CurveGeometry = std::variant<LineCurve3, CircleCurve3, NurbsCurve<double>>;
using SurfaceGeometry =
    std::variant<Plane3, Cylinder3, Cone3, Sphere3, Torus3, NurbsSurface<double>>;

struct Pcurve2d { Point2d p0{}; Point2d p1{}; };
struct GeometryStore {
  std::vector<CurveGeometry> curves;
  std::vector<SurfaceGeometry> surfaces;
  std::vector<Pcurve2d> pcurves;
};

class Brep {
 public:
  Brep() = default;

  bool vertex_alive(VertexId id) const noexcept;
  const Vertex& vertex(VertexId id) const noexcept;   // precondition: alive
  std::uint32_t vertex_count() const noexcept;
  bool edge_alive(EdgeId id) const noexcept;
  const Edge& edge(EdgeId id) const noexcept;
  std::uint32_t edge_count() const noexcept;
  bool coedge_alive(CoedgeId id) const noexcept;
  const Coedge& coedge(CoedgeId id) const noexcept;
  std::uint32_t coedge_count() const noexcept;
  bool loop_alive(LoopId id) const noexcept;
  const Loop& loop(LoopId id) const noexcept;
  std::uint32_t loop_count() const noexcept;
  bool face_alive(FaceId id) const noexcept;
  const Face& face(FaceId id) const noexcept;
  std::uint32_t face_count() const noexcept;
  bool shell_alive(ShellId id) const noexcept;
  const Shell& shell(ShellId id) const noexcept;
  std::uint32_t shell_count() const noexcept;
  bool solid_alive(SolidId id) const noexcept;
  const SolidEntity& solid(SolidId id) const noexcept;
  std::uint32_t solid_count() const noexcept;

  const std::vector<Vertex>& vertices() const noexcept;
  const std::vector<Edge>& edges() const noexcept;
  const std::vector<Coedge>& coedges() const noexcept;
  const std::vector<Loop>& loops() const noexcept;
  const std::vector<Face>& faces() const noexcept;
  const std::vector<Shell>& shells() const noexcept;
  const std::vector<SolidEntity>& solids() const noexcept;
  std::vector<SolidId> solid_ids() const noexcept;  // alive solid 的 id（构造顺序）
  const GeometryStore& geometry() const noexcept;
  std::span<const FaceId> shell_faces(ShellId id) const noexcept;
  std::span<const ShellId> solid_shells(SolidId id) const noexcept;

 private:
  friend class SolidEditor;
  // allocators + 每个实体的扁平数组；shell_faces_ / solid_shells_ 为全局池，
  // Shell::first_face / face_count 与 SolidEntity::first_shell / shell_count 指向它。
};

}  // namespace geom::v1
```

新增错误码：`Errc::invalid_topology`、`Errc::tolerance_violation`、`Errc::not_closed`（追加在 `self_intersecting` 之后，不改变既有值）。

- [ ] **Step 1: 写失败测试**

```cpp
TEST_CASE("empty brep has no live entities") {
  geom::Brep b;
  CHECK(b.vertex_count() == 0);
  CHECK(b.edge_count() == 0);
  CHECK(b.face_count() == 0);
  CHECK(b.solid_count() == 0);
  CHECK_FALSE(b.vertex_alive(geom::VertexId::null()));
  CHECK(b.geometry().curves.empty());
}

TEST_CASE("brep ids are generation-checked handles") {
  static_assert(std::is_same_v<geom::VertexId, geom::Handle>);
  const geom::VertexId v = geom::VertexId::make(3, 2);
  CHECK(v.index() == 3);
  CHECK(v.generation() == 2);
  CHECK_FALSE(v.is_null());
  CHECK(geom::VertexId::null().is_null());
}

TEST_CASE("geometry store preserves analytic type identity") {
  geom::GeometryStore g;
  g.surfaces.emplace_back(geom::Plane3{{0, 0, 1}, 0.0});
  CHECK(std::holds_alternative<geom::Plane3>(g.surfaces.back()));
  g.surfaces.emplace_back(geom::Cylinder3{});
  CHECK(std::holds_alternative<geom::Cylinder3>(g.surfaces.back()));
  g.curves.emplace_back(geom::CircleCurve3{});
  CHECK(std::holds_alternative<geom::CircleCurve3>(g.curves.back()));
}

TEST_CASE("new error codes stringify") {
  CHECK(std::string_view(geom::to_string(geom::Errc::invalid_topology)) == "invalid_topology");
  CHECK(std::string_view(geom::to_string(geom::Errc::tolerance_violation)) == "tolerance_violation");
  CHECK(std::string_view(geom::to_string(geom::Errc::not_closed)) == "not_closed");
}
```

- [ ] **Step 2: 运行确认失败**

Run: `cmake --build build --config Debug --target geom_tests` → 期望找不到 `geom/brep/topology.hpp`，且 `Errc::invalid_topology` 未定义。

- [ ] **Step 3: 实现**

Brep 内部持有 7 个 `HandleAllocator` 与 7 个 `std::vector<Entity>`；`add` 时 `alloc()`，若返回 null（容量 2^24 耗尽）不抛异常——`SolidEditor` 把它转成 `Errc::unsupported`（message "entity capacity exhausted"）。数组按 `handle.index()` 增长，空洞槽位保持默认值；`alive` 判定必须同时检查 index < size 与该类型 allocator 的 generation。`shell_faces` / `solid_shells` 返回 `std::span`，越界或 dangling handle 返回空 span（debug 下 `assert`）。`Brep` 的拷贝是值语义（向量与 allocator 一起复制），提交时返回副本。

**注意**：`Errc::not_closed` 与 `tolerance_violation` 本任务只加枚举与字符串；`coverage_gaps_test.cpp` 的枚举遍历断言同步补上。

- [ ] **Step 4: 运行确认通过**

Run: Debug/Release 全套 + `ctest`；`geom_tests.exe --test-case="*brep*"` 全过。

- [ ] **Step 5: 提交**

```bash
git add include/geom/brep/topology.hpp include/geom/core/error.hpp tests/unit/brep_topology_test.cpp tests/unit/coverage_gaps_test.cpp tests/CMakeLists.txt
git commit -m "feat(m2): add flat B-rep topology storage and read-only Brep"
```

---

### Task 4: SolidEditor 暂存编辑

**Files:**
- Create: `include/geom/brep/solid_editor.hpp`
- Test: `tests/unit/solid_editor_test.cpp`
- Modify: `tests/CMakeLists.txt`

**Interfaces:**
- Consumes: `Brep`、`Vertex/Edge/Coedge/Loop/Face/Shell/SolidEntity`、`GeometryStore`、`Result`。
- Produces:

```cpp
struct CoedgeSpec {
  EdgeId edge{};
  bool reversed = false;
  std::uint32_t pcurve = kNoGeometry;
  double tol = 0.0;
};

class SolidEditor {
 public:
  SolidEditor() = default;
  explicit SolidEditor(const Brep& base);     // 从已有 Brep 继续编辑
  const Brep& brep() const noexcept;          // 暂存视图

  std::uint32_t add_curve(CurveGeometry curve);
  std::uint32_t add_surface(SurfaceGeometry surface);
  std::uint32_t add_pcurve(Pcurve2d pcurve);

  Result<VertexId> add_vertex(Point3d p, double tol);
  Result<EdgeId> add_edge(VertexId v0, VertexId v1, std::uint32_t curve,
                          double t0, double t1, double tol);
  Result<FaceId> add_face(std::uint32_t surface, bool reversed, double tol);
  Result<LoopId> add_loop(FaceId face, std::span<const CoedgeSpec> coedges,
                          bool outer);
  Result<ShellId> add_shell(std::span<const FaceId> faces, bool closed);
  Result<SolidId> add_solid(std::span<const ShellId> shells, double tol);

  Result<Brep> commit() const;   // 校验通过返回新 Brep；失败不改变编辑器
  void rollback();               // 丢弃暂存，回到构造时的 base

 private:
  Brep staged_;
};
```

- [ ] **Step 1: 写失败测试**

```cpp
namespace {
struct SquareFace {
  geom::SolidEditor editor;
  geom::FaceId face;
  geom::LoopId loop;
  std::array<geom::VertexId, 4> vertices;
};
SquareFace make_square_face() {
  geom::SolidEditor ed;
  SquareFace s{std::move(ed), {}, {}, {}};
  const geom::Point3d pts[4] = {{0, 0, 0}, {1, 0, 0}, {1, 1, 0}, {0, 1, 0}};
  for (int i = 0; i < 4; ++i) s.vertices[i] = s.editor.add_vertex(pts[i], 0.0).value();
  const std::uint32_t plane = s.editor.add_surface(geom::Plane3{{0, 0, 1}, 0.0});
  const std::uint32_t line = s.editor.add_curve(
      geom::LineCurve3{{0, 0, 0}, {1, 0, 0}});
  std::array<geom::EdgeId, 4> edges{};
  for (int i = 0; i < 4; ++i) {
    edges[i] = s.editor.add_edge(s.vertices[i], s.vertices[(i + 1) % 4], line, 0.0, 1.0, 0.0)
                   .value();
  }
  s.face = s.editor.add_face(plane, false, 0.0).value();
  std::array<geom::CoedgeSpec, 4> specs{};
  for (int i = 0; i < 4; ++i) specs[i] = geom::CoedgeSpec{edges[i], false};
  s.loop = s.editor.add_loop(s.face, specs, true).value();
  auto shell = s.editor.add_shell(std::array<geom::FaceId, 1>{s.face}, false).value();
  s.editor.add_solid(std::array<geom::ShellId, 1>{shell}, 0.0).value();
  return s;
}
}  // namespace

TEST_CASE("editor commits a single-face open shell") {
  SquareFace s = make_square_face();
  const auto r = s.editor.commit();
  REQUIRE(r.has_value());
  const geom::Brep& b = r.value();
  CHECK(b.vertex_count() == 4);
  CHECK(b.edge_count() == 4);
  CHECK(b.coedge_count() == 4);
  CHECK(b.face_count() == 1);
  CHECK(b.solid_count() == 1);
  CHECK(b.shell(b.solid_shells(b.solid_ids().front()).front()).closed == false);
}

TEST_CASE("loop coedges form a closed cycle") {
  SquareFace s = make_square_face();
  const auto r = s.editor.commit();
  REQUIRE(r.has_value());
  const geom::Brep& b = r.value();
  const geom::Loop& loop = b.loop(s.loop);
  geom::CoedgeId id = loop.start;
  int steps = 0;
  do {
    const geom::Coedge& c = b.coedge(id);
    CHECK(c.loop == s.loop);
    CHECK(c.face == s.face);
    id = c.next;
    ++steps;
  } while (id != loop.start && steps <= 4);
  CHECK(steps == 4);
}

TEST_CASE("commit rejects a loop with no coedges") {
  geom::SolidEditor ed;
  const std::uint32_t plane = ed.add_surface(geom::Plane3{{0, 0, 1}, 0.0});
  const geom::FaceId f = ed.add_face(plane, false, 0.0).value();
  const auto r = ed.add_loop(f, std::span<const geom::CoedgeSpec>{}, true);
  CHECK_FALSE(r.has_value());
  CHECK(r.error().code == geom::Errc::invalid_topology);
  CHECK(ed.brep().face_count() == 1);  // 失败的 add_loop 没有破坏已有实体
  CHECK(ed.brep().loop_count() == 0);
}

TEST_CASE("rollback discards staged entities") {
  geom::SolidEditor ed;
  CHECK(ed.add_vertex({1, 2, 3}, 0.0).has_value());
  CHECK(ed.brep().vertex_count() == 1);
  ed.rollback();
  CHECK(ed.brep().vertex_count() == 0);
}
```

- [ ] **Step 2: 运行确认失败**

Run: `cmake --build build --config Debug --target geom_tests` → 期望找不到 `geom/brep/solid_editor.hpp`。

- [ ] **Step 3: 实现**

`add_vertex`：拒绝非有限坐标 → `Errc::non_finite`；`tol < 0` → `Errc::tolerance_violation`。`add_edge`：两个顶点必须 alive（否则 `index_out_of_range`），`t0 < t1` 且有限（否则 `degenerate`），curve 索引 < `curves.size()` 或 `kNoGeometry`。`add_face`：surface 索引合法或 `kNoGeometry`。`add_loop`：coedge 数量 ≥ 1（否则 `invalid_topology`），每条 edge alive，同一 loop 内 edge 出现次数 ≤ 2，`next/prev` 按给定顺序成环。`add_shell`：faces 非空且全部 alive；`add_solid`：shells 非空且 alive。所有失败返回 `Result` 的错误分支，不修改 `staged_`。

`commit()`：在副本上检查（1）所有 alive 实体引用的 handle 都 alive；（2）每个 loop 的 coedge 链首尾相接且 `coedge.loop/face` 一致；（3）每个 face 至多一个 `outer == true` 的 loop。通过则返回副本，否则 `Errc::invalid_topology`；`T5/T6` 会把完整校验器接进来。

- [ ] **Step 4: 运行确认通过**

Run: Debug/Release 全套 + `ctest`；`--test-case="*editor*"` 全过。

- [ ] **Step 5: 提交**

```bash
git add include/geom/brep/solid_editor.hpp tests/unit/solid_editor_test.cpp tests/CMakeLists.txt
git commit -m "feat(m2): add staged SolidEditor with commit and rollback"
```

### Task 5: 拓扑校验器与 Euler 特征数

**Files:**
- Create: `include/geom/brep/validate.hpp`
- Modify: `include/geom/brep/solid_editor.hpp`（接入校验；新增低层 `erase`）
- Test: `tests/unit/brep_validate_test.cpp`
- Modify: `tests/CMakeLists.txt`

**Interfaces:**
- Consumes: `Brep`、`Errc`。
- Produces:

```cpp
struct ValidationIssue {
  Errc code = Errc::none;
  std::uint32_t entity = 0xFFFFFFFFu;  // handle.index()，UINT32_MAX 表示全局
  const char* message = "";
};

struct ValidationReport {
  std::vector<ValidationIssue> issues;
  bool ok() const noexcept { return issues.empty(); }
};

ValidationReport validate_topology(const Brep& brep);
std::int64_t euler_characteristic(const Brep& brep, SolidId solid);

// SolidEditor 新增（注意 EdgeId/VertexId 都是 Handle 别名，不能用重载区分）：
void erase_edge(EdgeId id);      // 只释放 handle，不级联修复引用
void erase_vertex(VertexId id);
```

- [ ] **Step 1: 写失败测试**

先在本测试文件顶部复制 T4 的 `SquareFace` / `make_square_face()` 辅助（各测试目标独立编译，不能跨 TU 共享匿名命名空间符号）。

```cpp
TEST_CASE("square face passes topology validation") {
  SquareFace s = make_square_face();
  const auto r = s.editor.commit();
  REQUIRE(r.has_value());
  CHECK(geom::validate_topology(r.value()).ok());
  CHECK(geom::euler_characteristic(r.value(), r.value().solid_ids().front()) == 1);
}

TEST_CASE("dangling edge reference is rejected at commit") {
  SquareFace s = make_square_face();
  const geom::EdgeId e0 = s.editor.brep().coedge(s.editor.brep().loop(s.loop).start).edge;
  s.editor.erase_edge(e0);  // 低层破坏：coedge 仍引用已释放的边
  const auto r = s.editor.commit();
  REQUIRE_FALSE(r.has_value());
  CHECK(r.error().code == geom::Errc::invalid_topology);
}

TEST_CASE("edge referenced by three coedges is rejected") {
  SquareFace s = make_square_face();
  const geom::EdgeId e0 = s.editor.brep().coedge(s.editor.brep().loop(s.loop).start).edge;
  const std::array<geom::CoedgeSpec, 3> specs{
      geom::CoedgeSpec{e0, false}, geom::CoedgeSpec{e0, true},
      geom::CoedgeSpec{e0, false}};
  const auto r = s.editor.add_loop(s.face, specs, false);
  CHECK_FALSE(r.has_value());
  CHECK(r.error().code == geom::Errc::invalid_topology);
}
```

- [ ] **Step 2: 运行确认失败**

Run: `cmake --build build --config Debug --target geom_tests` → 期望找不到 `geom/brep/validate.hpp`。

- [ ] **Step 3: 实现**

`validate_topology` 逐项检查并把第一处失败的 `Errc` 与实体 index 填入 issue（可继续收集多处，但最多 8 条）：

1. 每个 alive 实体引用的 handle 必须 alive：`Coedge.edge/face/loop/next/prev`、`Loop.face/start/next`、`Face.first_loop`、`Edge.v0/v1`。
2. 每个 loop：从 `start` 出发沿 `next` 必须回到 `start`，`prev` 反向一致，`coedge.loop == loop` 且 `coedge.face == loop.face`。
3. 每个 face：至多一个 `outer == true` 的 loop；loop 链无环。
4. 每条 edge：coedge 数 ∈ {1, 2}；若为 2，两个 coedge 的 `reversed` 必须相反（圆柱接缝正是这种「同边两向」）。
5. 每个 `closed == true` 的 shell：其所有面的边都必须是 2 个 coedge，否则报 `Errc::not_closed`；`closed == false` 的 shell 允许边界边（open sheet），但 `mass_properties` 会拒绝这样的实体。
6. `euler_characteristic`：遍历 solid 的 shell→face→loop→coedge 收集去重后的 V/E/F 集合，返回 `V - E + F - (L - F)`（L 为 loop 总数，L−F 为内环数）。**不**在 validator 中断言等于 2：M2 的环面是 0、单面片是 1，断言由每个基本体的测试负责。

`SolidEditor::commit` 依次调用 `validate_topology` 与 `validate_tolerances`（T6 加入），任一失败返回 `Result` 错误（`code` 取第一条 issue，`context_a` 取其 entity）。

- [ ] **Step 4: 运行确认通过**

Run: Debug/Release 全套 + `ctest`；`--test-case="*validate*"` 全过。

- [ ] **Step 5: 提交**

```bash
git add include/geom/brep/validate.hpp include/geom/brep/solid_editor.hpp tests/unit/brep_validate_test.cpp tests/CMakeLists.txt
git commit -m "feat(m2): validate B-rep topology and report Euler characteristic"
```

---

### Task 6: 容差链校验（§4.6）

**Files:**
- Modify: `include/geom/brep/validate.hpp`（新增 `validate_tolerances`）
- Modify: `include/geom/brep/solid_editor.hpp`（commit 接入）
- Test: `tests/unit/brep_tolerance_test.cpp`
- Modify: `tests/CMakeLists.txt`

**Interfaces:**
- Produces: `ValidationReport validate_tolerances(const Brep& brep);`
- 规则常量：`inline constexpr double kToleranceGrowthFactor = 10.0;`（§4.6 的「增长 >10× 视为可疑」）。

- [ ] **Step 1: 写失败测试**

```cpp
namespace {
geom::SolidEditor make_square_editor(double vertex_tol, double edge_tol,
                                     double face_tol, double solid_tol) {
  geom::SolidEditor ed;
  const geom::Point3d pts[4] = {{0, 0, 0}, {1, 0, 0}, {1, 1, 0}, {0, 1, 0}};
  std::array<geom::VertexId, 4> v{};
  for (int i = 0; i < 4; ++i) v[i] = ed.add_vertex(pts[i], vertex_tol).value();
  const std::uint32_t plane = ed.add_surface(geom::Plane3{{0, 0, 1}, 0.0});
  const std::uint32_t line = ed.add_curve(geom::LineCurve3{{0, 0, 0}, {1, 0, 0}});
  std::array<geom::EdgeId, 4> e{};
  for (int i = 0; i < 4; ++i) {
    e[i] = ed.add_edge(v[i], v[(i + 1) % 4], line, 0.0, 1.0, edge_tol).value();
  }
  const geom::FaceId f = ed.add_face(plane, false, face_tol).value();
  std::array<geom::CoedgeSpec, 4> specs{};
  for (int i = 0; i < 4; ++i) specs[i] = geom::CoedgeSpec{e[i], false};
  ed.add_loop(f, specs, true).value();
  const auto sh = ed.add_shell(std::array<geom::FaceId, 1>{f}, false).value();
  ed.add_solid(std::array<geom::ShellId, 1>{sh}, solid_tol).value();
  return ed;
}
}  // namespace

TEST_CASE("edge tolerance must dominate its vertices") {
  auto ed = make_square_editor(/*vertex=*/0.1, /*edge=*/0.05, /*face=*/0.2,
                               /*solid=*/0.2);
  const auto r = ed.commit();
  REQUIRE_FALSE(r.has_value());
  CHECK(r.error().code == geom::Errc::tolerance_violation);
}

TEST_CASE("solid tolerance must dominate its faces") {
  auto ed = make_square_editor(0.0, 0.0, /*face=*/1.0, /*solid=*/0.1);
  const auto r = ed.commit();
  REQUIRE_FALSE(r.has_value());
  CHECK(r.error().code == geom::Errc::tolerance_violation);
}

TEST_CASE("tolerance growth beyond ten times is a violation") {
  auto ed = make_square_editor(0.0, /*edge=*/0.1, /*face=*/2.0, /*solid=*/2.0);
  const auto r = ed.commit();
  REQUIRE_FALSE(r.has_value());
  CHECK(r.error().code == geom::Errc::tolerance_violation);
  CHECK(std::string_view(r.error().message).find("growth") != std::string_view::npos);
}

TEST_CASE("valid tolerance chain is accepted") {
  auto ed = make_square_editor(0.01, 0.02, 0.03, 0.03);
  const auto r = ed.commit();
  REQUIRE(r.has_value());
  CHECK(geom::validate_tolerances(r.value()).ok());
}

TEST_CASE("negative tolerance is rejected at construction") {
  auto ed = make_square_editor(0.0, 0.0, 0.0, 0.0);
  CHECK_FALSE(ed.add_vertex({0, 0, 0}, -1.0).has_value());
  CHECK_FALSE(ed.add_face(0, false, -0.5).has_value());
}
```

- [ ] **Step 2: 运行确认失败**

Run: `cmake --build build --config Debug --target geom_tests` → 期望 `validate_tolerances` 未声明。

- [ ] **Step 3: 实现**

规则（每条失败返回一条 issue，`entity` 为子实体 index）：

1. 所有 tol 必须有限且 ≥ 0（`tolerance_violation`）。
2. `edge.tol ≥ max(vertex(v0).tol, vertex(v1).tol)`。
3. 对每条用到该边的 coedge：`face.tol ≥ edge.tol`；`coedge.tol ≤ face.tol`（pcurve 容差属于面内对象，不得大于所属面），coedge 默认 0 因此恒满足。
4. `solid.tol ≥ max(shell 内所有 face.tol)`。
5. 增长检查：`parent.tol > 0 && child.tol > 0 && parent.tol > kToleranceGrowthFactor * child.tol` 时记一条 `tolerance_violation`，`message = "tolerance growth exceeds 10x"`（不引入新 Errc）。

- [ ] **Step 4: 运行确认通过**

Run: Debug/Release 全套 + `ctest`；`--test-case="*tolerance*"` 全过。

- [ ] **Step 5: 提交**

```bash
git add include/geom/brep/validate.hpp include/geom/brep/solid_editor.hpp tests/unit/brep_tolerance_test.cpp tests/CMakeLists.txt
git commit -m "feat(m2): enforce the B-rep tolerance chain"
```

---

### Task 7: Euler 算子 MEV / KEV

**Files:**
- Modify: `include/geom/brep/solid_editor.hpp`
- Test: `tests/unit/brep_euler_ops_test.cpp`
- Modify: `tests/CMakeLists.txt`

**Interfaces:**
- Produces:

```cpp
Result<VertexId> split_edge(EdgeId edge, double t);  // MEV
Result<bool> join_vertex(VertexId vertex);           // KEV
```

- [ ] **Step 1: 写失败测试**

先在本测试文件顶部复制 T4 的 `SquareFace` / `make_square_face()` 辅助（放在匿名命名空间，避免与其他测试 TU 的同名符号冲突）。

```cpp
TEST_CASE("split edge inserts a vertex and keeps the loop closed") {
  SquareFace s = make_square_face();
  const geom::Brep before = s.editor.commit().value();
  const geom::EdgeId e0 = before.coedge(before.loop(s.loop).start).edge;
  REQUIRE(s.editor.split_edge(e0, 0.5).has_value());
  const auto after = s.editor.commit();
  REQUIRE(after.has_value());
  const geom::Brep& b = after.value();
  CHECK(b.vertex_count() == before.vertex_count() + 1);
  CHECK(b.edge_count() == before.edge_count() + 1);
  CHECK(b.coedge_count() == before.coedge_count() + 1);
  CHECK(geom::validate_topology(b).ok());
  CHECK(geom::euler_characteristic(b, b.solid_ids().front()) == 1);
}

TEST_CASE("join vertex restores the original counts") {
  SquareFace s = make_square_face();
  const geom::EdgeId e0 = s.editor.brep().coedge(s.editor.brep().loop(s.loop).start).edge;
  const geom::VertexId mid = s.editor.split_edge(e0, 0.5).value();
  CHECK(s.editor.join_vertex(mid).has_value());
  const auto r = s.editor.commit();
  REQUIRE(r.has_value());
  CHECK(r.value().vertex_count() == 4);
  CHECK(r.value().edge_count() == 4);
  CHECK(geom::validate_topology(r.value()).ok());
}

TEST_CASE("split at an endpoint is rejected") {
  SquareFace s = make_square_face();
  const geom::EdgeId e0 = s.editor.brep().coedge(s.editor.brep().loop(s.loop).start).edge;
  CHECK_FALSE(s.editor.split_edge(e0, 0.0).has_value());
  CHECK_FALSE(s.editor.split_edge(e0, 1.0).has_value());
  CHECK(s.editor.split_edge(e0, 0.25).has_value());
}

TEST_CASE("join rejects a vertex with three incident edges") {
  SquareFace s = make_square_face();
  const geom::VertexId extra =
      s.editor.add_vertex({0.5, -1.0, 0.0}, 0.0).value();
  const std::uint32_t line =
      s.editor.add_curve(geom::LineCurve3{{1, 0, 0}, {0, 1, 0}});
  REQUIRE(s.editor.add_edge(s.vertices[1], extra, line, 0.0, 1.0, 0.0).has_value());
  const auto joined = s.editor.join_vertex(s.vertices[1]);
  REQUIRE_FALSE(joined.has_value());
  CHECK(joined.error().code == geom::Errc::invalid_topology);
}
```

- [ ] **Step 2: 运行确认失败**

Run: `cmake --build build --config Debug --target geom_tests` → `split_edge` 未声明。

- [ ] **Step 3: 实现**

`split_edge(e, t)`：边需 alive，`t0 < t < t1`，曲线可求值；新顶点 `p = curve.eval(t)`、tol = `edge.tol`；创建 e0(v0→vnew, t0..t)、e1(vnew→v1, t..t1)，两者 tol 与 pcurve 继承；遍历该边的每个 coedge，在原 loop 位置用两个新 coedge 替换（方向继承，pcurve 按线性段插值拆分或置 `kNoGeometry`）；释放旧边与旧 coedge。返回新顶点。

`join_vertex(v)`：v 必须 alive 且恰好 2 条 incident edge（按 coedge 统计；否则 `Errc::invalid_topology`），两条边共享同一曲线且参数区间首尾相接（`|e0.t1 - e1.t0| <= 0` 精确或同端点）、在 v 处切线方向一致（点积 > 0，否则 `degenerate`）；合并为一条边（t0 = min、t1 = max），重写两个 coedge 的 edge 引用并重新接线，释放 v 与其中一条边。

- [ ] **Step 4: 运行确认通过**

Run: Debug/Release 全套 + `ctest`；`--test-case="*split*"`、`--test-case="*join*"` 全过。

- [ ] **Step 5: 提交**

```bash
git add include/geom/brep/solid_editor.hpp tests/unit/brep_euler_ops_test.cpp tests/CMakeLists.txt
git commit -m "feat(m2): add MEV and KEV Euler operators"
```

### Task 8: 长方体基本体

**Files:**
- Create: `include/geom/brep/primitives.hpp`
- Test: `tests/unit/brep_box_test.cpp`
- Modify: `tests/CMakeLists.txt`

**Interfaces:**
- Consumes: `SolidEditor`、`LineCurve3`、`Plane3`、`Point3d`。
- Produces:

```cpp
Result<Brep> make_box(Point3d lo, Point3d hi, double tol);
```

- [ ] **Step 1: 写失败测试**

```cpp
TEST_CASE("box has the expected topology") {
  const auto r = geom::make_box({0, 0, 0}, {2, 3, 4}, 0.0);
  REQUIRE(r.has_value());
  const geom::Brep& b = r.value();
  CHECK(b.vertex_count() == 8);
  CHECK(b.edge_count() == 12);
  CHECK(b.face_count() == 6);
  CHECK(b.coedge_count() == 24);
  CHECK(b.solid_count() == 1);
  CHECK(geom::validate_topology(b).ok());
  CHECK(geom::validate_tolerances(b).ok());
  const geom::SolidId s = b.solid_ids().front();
  CHECK(geom::euler_characteristic(b, s) == 2);
  CHECK(b.shell(b.solid_shells(s).front()).closed);
}

TEST_CASE("box face normals point outward") {
  const auto r = geom::make_box({-1, -2, -3}, {1, 2, 3}, 0.0);
  REQUIRE(r.has_value());
  const geom::Brep& b = r.value();
  for (const geom::Face& f : b.faces()) {
    const auto& plane = std::get<geom::Plane3>(b.geometry().surfaces[f.surface]);
    const double h = std::hypot(plane.normal[0], plane.normal[1], plane.normal[2]);
    const geom::Vec3d unit{plane.normal[0] / h, plane.normal[1] / h,
                           plane.normal[2] / h};
    const geom::CoedgeId c0 = b.loop(f.first_loop).start;
    const geom::Point3d sample = b.vertex(b.edge(b.coedge(c0).edge).v0).point;
    CHECK(geom::dot(unit, geom::Vec3d{sample[0], sample[1], sample[2]}) > 0.0);
  }
}

TEST_CASE("degenerate and non-finite boxes are rejected") {
  CHECK_FALSE(geom::make_box({0, 0, 0}, {0, 1, 1}, 0.0).has_value());
  CHECK_FALSE(geom::make_box({0, 0, 0}, {1, 1, 1}, -1.0).has_value());
  const double nan = std::numeric_limits<double>::quiet_NaN();
  CHECK_FALSE(geom::make_box({0, 0, 0}, {1, nan, 1}, 0.0).has_value());
}
```

- [ ] **Step 2: 运行确认失败**

Run: `cmake --build build --config Debug --target geom_tests` → `make_box` 未声明。

- [ ] **Step 3: 实现**

8 个角点（`tol` 全部相同）、12 条边（每条用 `LineCurve3{v0.point, unit_dir}`，`t0 = 0`、`t1 = |v1 - v0|`，方向为单位向量；边 tol 与顶点相同）。6 个面按「外法向 = 面内环的右手方向」构造：底面 `-z`、顶面 `+z`、前面 `-y`、后面 `+y`、左面 `-x`、右面 `+x`；每个 loop 的 coedge 顺序必须让平面法向与环方向一致（在平面上按外法向看是逆时针）。面的 `reversed = false`，surface 用 `Plane3{normal, dot(normal, corner)}`。最后 `add_shell(6 faces, closed=true)`、`add_solid`、`commit()`。任何 `add_*` 失败立即返回其错误（`Result` 传播）。

- [ ] **Step 4: 运行确认通过**

Run: Debug/Release 全套 + `ctest`；`--test-case="*box*"` 全过。

- [ ] **Step 5: 提交**

```bash
git add include/geom/brep/primitives.hpp tests/unit/brep_box_test.cpp tests/CMakeLists.txt
git commit -m "feat(m2): build box solids through the editor"
```

---

### Task 9: 圆柱基本体

**Files:**
- Modify: `include/geom/brep/primitives.hpp`
- Test: `tests/unit/brep_cylinder_test.cpp`
- Modify: `tests/CMakeLists.txt`

**Interfaces:**
- Produces: `Result<Brep> make_cylinder(Point3d base, Vec3d axis, double radius, double height, double tol);`

- [ ] **Step 1: 写失败测试**

```cpp
TEST_CASE("cylinder topology uses a seam edge twice") {
  const auto r = geom::make_cylinder({0, 0, 0}, {0, 0, 1}, 2.0, 5.0, 0.0);
  REQUIRE(r.has_value());
  const geom::Brep& b = r.value();
  CHECK(b.vertex_count() == 2);   // 底/顶接缝点
  CHECK(b.edge_count() == 3);     // 底圆、顶圆、接缝
  CHECK(b.face_count() == 3);     // 底盖、顶盖、侧面
  CHECK(b.coedge_count() == 6);
  CHECK(geom::validate_topology(b).ok());
  CHECK(geom::euler_characteristic(b, b.solid_ids().front()) == 2);
  CHECK(b.shell(b.solid_shells(b.solid_ids().front()).front()).closed);
}

TEST_CASE("cylinder faces have consistent outward normals") {
  const auto r = geom::make_cylinder({0, 0, 0}, {0, 0, 1}, 2.0, 5.0, 0.0);
  REQUIRE(r.has_value());
  const geom::Brep& b = r.value();
  for (const geom::Face& f : b.faces()) {
    if (!std::holds_alternative<geom::Cylinder3>(b.geometry().surfaces[f.surface])) continue;
    const auto sample = b.vertex(b.edge(b.coedge(b.loop(f.first_loop).start).edge).v0).point;
    const auto& cyl = std::get<geom::Cylinder3>(b.geometry().surfaces[f.surface]);
    const double u = std::atan2(sample[1], sample[0]);
    CHECK(geom::dot(cyl.normal(u, sample[2]), geom::Vec3d{sample[0], sample[1], 0}) > 0.0);
  }
}

TEST_CASE("cylinder rejects bad radius, height and axis") {
  CHECK_FALSE(geom::make_cylinder({0, 0, 0}, {0, 0, 1}, 0.0, 1.0, 0.0).has_value());
  CHECK_FALSE(geom::make_cylinder({0, 0, 0}, {0, 0, 1}, 1.0, -1.0, 0.0).has_value());
  CHECK_FALSE(geom::make_cylinder({0, 0, 0}, {0, 0, 2}, 1.0, 1.0, 0.0).has_value());
}
```

- [ ] **Step 2: 运行确认失败**

Run: 期望 `make_cylinder` 未声明。

- [ ] **Step 3: 实现**

构造规范：`x = 任一与 axis 不平行的单位向量投影后归一化`（取 `|axis[0]| < 0.9 ? (1,0,0) : (0,1,0)`，再 `x = normalize(cross(cross(axis, cand), axis))`… 实现只需 `x` 单位且 ⊥ axis），`y = cross(axis, x)`，`top = base + height*axis`，接缝点 `p0 = base + radius*x`、`p1 = top + radius*x`。

三条边：
- 底圆：`CircleCurve3{base, x, y, radius}`，`v0 = v1 = p0`，`t0=0, t1=2π`；
- 顶圆：同上平移到 `top`，`v0 = v1 = p1`；
- 接缝：`LineCurve3{p0, axis}`，`t0=0, t1=height`，`v0=p0, v1=p1`。

三个面与三个 loop（这是本任务最容易错的地方，coedge 方向必须满足「边恰有两个反向 coedge」）：
- 底盖 `Plane3{-axis, dot(-axis, base)}`：loop 只含底圆，`reversed = true`（外法向朝下，圆参数方向朝上时需反向）。
- 顶盖 `Plane3{axis, dot(axis, top)}`：loop 只含顶圆，`reversed = false`。
- 侧面 `Cylinder3{base, axis, x, radius}`：loop 依次为 `底圆(false) → 接缝(false) → 顶圆(true) → 接缝(true)`，共 4 个 coedge；两处接缝必须是同一 `EdgeId` 的正反两个 coedge。

`add_shell(3 faces, closed=true)`，`add_solid`。

- [ ] **Step 4: 运行确认通过**

Run: Debug/Release 全套 + `ctest`；`--test-case="*cylinder*"` 全过。

- [ ] **Step 5: 提交**

```bash
git add include/geom/brep/primitives.hpp tests/unit/brep_cylinder_test.cpp tests/CMakeLists.txt
git commit -m "feat(m2): build cylinder solids with a seam edge"
```

---

### Task 10: 圆锥与球基本体

**Files:**
- Modify: `include/geom/brep/primitives.hpp`
- Test: `tests/unit/brep_cone_sphere_test.cpp`
- Modify: `tests/CMakeLists.txt`

**Interfaces:**
- Produces:
  - `Result<Brep> make_cone(Point3d base, Vec3d axis, double base_radius, double height, double tol);`（`half_angle = atan2(base_radius, height)` 由参数导出）
  - `Result<Brep> make_sphere(Point3d center, double radius, double tol);`

- [ ] **Step 1: 写失败测试**

```cpp
TEST_CASE("cone topology shares the base circle and the seam") {
  const auto r = geom::make_cone({0, 0, 0}, {0, 0, 1}, 2.0, 3.0, 0.0);
  REQUIRE(r.has_value());
  const geom::Brep& b = r.value();
  CHECK(b.vertex_count() == 2);   // 底接缝点、锥顶
  CHECK(b.edge_count() == 2);     // 底圆、接缝
  CHECK(b.face_count() == 2);     // 底盖、侧面
  CHECK(b.coedge_count() == 4);
  CHECK(geom::validate_topology(b).ok());
  CHECK(geom::euler_characteristic(b, b.solid_ids().front()) == 2);
}

TEST_CASE("sphere is a single face with a two-sided seam") {
  const auto r = geom::make_sphere({0, 0, 0}, 3.0, 0.0);
  REQUIRE(r.has_value());
  const geom::Brep& b = r.value();
  CHECK(b.vertex_count() == 2);   // 南北极
  CHECK(b.edge_count() == 1);     // 子午线接缝
  CHECK(b.face_count() == 1);
  CHECK(b.coedge_count() == 2);
  CHECK(geom::validate_topology(b).ok());
  CHECK(geom::euler_characteristic(b, b.solid_ids().front()) == 2);
}

TEST_CASE("cone and sphere reject degenerate parameters") {
  CHECK_FALSE(geom::make_cone({0, 0, 0}, {0, 0, 1}, 0.0, 1.0, 0.0).has_value());
  CHECK_FALSE(geom::make_cone({0, 0, 0}, {0, 0, 1}, 1.0, 0.0, 0.0).has_value());
  CHECK_FALSE(geom::make_cone({0, 0, 0}, {0, 0, 2}, 1.0, 1.0, 0.0).has_value());
  CHECK_FALSE(geom::make_sphere({0, 0, 0}, 0.0, 0.0).has_value());
  CHECK_FALSE(geom::make_sphere({0, 0, 0}, -1.0, 0.0).has_value());
}
```

- [ ] **Step 2: 运行确认失败**

Run: 期望 `make_cone` / `make_sphere` 未声明。

- [ ] **Step 3: 实现**

圆锥：底圆 `v0 = v1 = p_base`、接缝 `p_base → apex`（`apex = base + height*axis`，`base_radius > 0`、`height > 0`、`axis` 单位；否则 `Errc::degenerate` / `non_finite`）。`Cone3.half_angle = atan2(base_radius, height)`，这样顶点与底圆在几何上严格一致，不引入第二个可矛盾的参数。侧面 `Cone3{base, axis, x, base_radius, half_angle}`：loop = `底圆(false) → 接缝(false) → 接缝(true)`；底盖 `Plane3{-axis, ...}`：底圆 `reversed = true`。共 2 边、4 coedge。

球：接缝是半条大圆（从南极到北极，`CircleCurve3` 参数 `[-π/2, π/2]`，平面含 axis 与 x）；`v_south`、`v_north` 两个顶点；单面 `Sphere3`，loop 为 `接缝(false) → 接缝(true)`；`add_shell(1 face, closed=true)`。球面 `Sphere3` 的 `eval/du/dv/normal` 在极点按 T2 的规则取极限。

- [ ] **Step 4: 运行确认通过**

Run: Debug/Release 全套 + `ctest`；`--test-case="*cone*"`、`--test-case="*sphere*"` 全过。

- [ ] **Step 5: 提交**

```bash
git add include/geom/brep/primitives.hpp tests/unit/brep_cone_sphere_test.cpp tests/CMakeLists.txt
git commit -m "feat(m2): build cone and sphere solids"
```

---

### Task 11: 环面基本体

**Files:**
- Modify: `include/geom/brep/primitives.hpp`
- Test: `tests/unit/brep_torus_test.cpp`
- Modify: `tests/CMakeLists.txt`

**Interfaces:**
- Produces: `Result<Brep> make_torus(Point3d center, Vec3d axis, double major_radius, double minor_radius, double tol);`

- [ ] **Step 1: 写失败测试**

```cpp
TEST_CASE("torus has Euler characteristic zero") {
  const auto r = geom::make_torus({0, 0, 0}, {0, 0, 1}, 3.0, 1.0, 0.0);
  REQUIRE(r.has_value());
  const geom::Brep& b = r.value();
  CHECK(b.vertex_count() == 1);   // 两个接缝的交点
  CHECK(b.edge_count() == 2);     // u 接缝、v 接缝（都是闭合边）
  CHECK(b.face_count() == 1);
  CHECK(b.coedge_count() == 4);
  CHECK(geom::validate_topology(b).ok());
  CHECK(geom::euler_characteristic(b, b.solid_ids().front()) == 0);
}

TEST_CASE("torus rejects r >= R and non-unit axis") {
  CHECK_FALSE(geom::make_torus({0, 0, 0}, {0, 0, 1}, 1.0, 1.0, 0.0).has_value());
  CHECK_FALSE(geom::make_torus({0, 0, 0}, {0, 0, 2}, 3.0, 1.0, 0.0).has_value());
  CHECK_FALSE(geom::make_torus({0, 0, 0}, {0, 0, 1}, 0.0, 0.1, 0.0).has_value());
}
```

- [ ] **Step 2: 运行确认失败**

Run: 期望 `make_torus` 未声明。

- [ ] **Step 3: 实现**

接缝顶点 `p = center + (R + r)*x`（u=0, v=0）。u 接缝：`CircleCurve3` 在平面 z = 0 内绕轴，半径 `R + r`，`v0 = v1 = p`，`t0=0,t1=2π`；v 接缝：`CircleCurve3` 在含 `x` 与 `axis` 的平面内、圆心 `center + R*x`、半径 `r`，`v0 = v1 = p`，`t0=0,t1=2π`。单面 `Torus3`，loop = `u 接缝(false) → v 接缝(false) → u 接缝(true) → v 接缝(true)`。`shell(1 face, closed=true)`。

- [ ] **Step 4: 运行确认通过**

Run: Debug/Release 全套 + `ctest`；`--test-case="*torus*"` 全过。

- [ ] **Step 5: 提交**

```bash
git add include/geom/brep/primitives.hpp tests/unit/brep_torus_test.cpp tests/CMakeLists.txt
git commit -m "feat(m2): build torus solids with two seam edges"
```

### Task 12: 网格化——边采样与平面面

**Files:**
- Create: `include/geom/brep/tessellate.hpp`
- Test: `tests/unit/tessellate_planar_test.cpp`
- Modify: `tests/CMakeLists.txt`

**Interfaces:**
- Consumes: `Brep`、`triangulate_polygon2`、`Polygon2d`、`Plane3`、`LineCurve3`、`CircleCurve3`、`RingOrientation`。
- Produces:

```cpp
struct TriMesh {
  std::vector<Point3d> positions;
  std::vector<std::array<std::uint32_t, 3>> triangles;
  std::vector<std::uint32_t> triangle_face;  // 与 triangles 等长
};

struct TessellationOptions {
  double deflection = 1e-3;        // 弦高上限，> 0
  double angle_tolerance = 0.35;   // 相邻三角形法向夹角上限（弧度）
};

Result<TriMesh> tessellate(const Brep& brep, const TessellationOptions& options);
```

- [ ] **Step 1: 写失败测试**

```cpp
// 需要 T4 的 SquareFace / make_square_face()（复制到本文件匿名命名空间）。
TEST_CASE("single square face tessellates into two triangles") {
  SquareFace s = make_square_face();
  const auto b = s.editor.commit();
  REQUIRE(b.has_value());
  const auto result = geom::tessellate(b.value(), {});
  REQUIRE(result.has_value());
  const geom::TriMesh& mesh = result.value();
  CHECK(mesh.positions.size() == 4);
  CHECK(mesh.triangles.size() == 2);
  CHECK(mesh.triangle_face.size() == 2);
}

TEST_CASE("box mesh is closed and edge-manifold") {
  const auto b = geom::make_box({0, 0, 0}, {2, 3, 4}, 0.0);
  REQUIRE(b.has_value());
  const auto result = geom::tessellate(b.value(), {});
  REQUIRE(result.has_value());
  const geom::TriMesh& mesh = result.value();
  CHECK(mesh.positions.size() == 8);
  CHECK(mesh.triangles.size() == 12);
  std::map<std::pair<std::uint32_t, std::uint32_t>, int> edge_use;
  for (const auto& t : mesh.triangles) {
    for (int i = 0; i < 3; ++i) {
      auto a = t[i];
      auto c = t[(i + 1) % 3];
      if (a > c) std::swap(a, c);
      ++edge_use[{a, c}];
    }
  }
  CHECK(edge_use.size() == 18);
  for (const auto& [edge, count] : edge_use) CHECK(count == 2);
}

TEST_CASE("tessellation rejects non-positive deflection") {
  const auto b = geom::make_box({0, 0, 0}, {1, 1, 1}, 0.0);
  REQUIRE(b.has_value());
  CHECK_FALSE(geom::tessellate(b.value(), {0.0, 0.35}).has_value());
  CHECK_FALSE(geom::tessellate(b.value(), {-1.0, 0.35}).has_value());
}
```

- [ ] **Step 2: 运行确认失败**

Run: 期望找不到 `geom/brep/tessellate.hpp`。

- [ ] **Step 3: 实现**

1. **边采样缓存**：按 `(EdgeId, deflection)` 生成采样点；`LineCurve3` 取两端点；`CircleCurve3` 取 `n = ceil(|t1-t0| / dθ)`，`dθ = 2*acos(max(-1.0, 1.0 - deflection/radius))`（`deflection >= 2r` 时取 π）；其他曲线返回 `Errc::unsupported`。采样点必须按边的参数方向保存。
2. **焊接表**：端点按 `VertexId` 焊接（同一拓扑顶点永远映射到同一个全局下标，圆柱/圆锥的接缝端点因此与端盖边界自动重合）；边的内部采样按 `(edge.index(), 采样序号)` 焊接。相邻面因此复用同一批顶点（水密的基础）。
3. **平面面**：用面所在平面构造正交基 `(u_axis, v_axis = cross(n, u_axis))`，把 loop 的边界采样点投影为 `Point2d`；按 `outer` 决定环方向（外环 CCW、内环 CW），组装 `Polygon2d{vertices, ring_starts}`；调 `triangulate_polygon2`；三角形法向与 `face_normal = (reversed ? -1 : 1) * unit(plane.normal)` 同向，反向则交换两个顶点。
4. 面上每个 loop 的第一个 coedge 从采样序号 0 开始；相邻 coedge 共享的端点只保留一次（`reversed` 时从 `n` 递减到 0）。空 loop、零面积面返回 `Errc::degenerate`。

- [ ] **Step 4: 运行确认通过**

Run: Debug/Release 全套 + `ctest`；`--test-case="*tessellat*"` 全过。

- [ ] **Step 5: 提交**

```bash
git add include/geom/brep/tessellate.hpp tests/unit/tessellate_planar_test.cpp tests/CMakeLists.txt
git commit -m "feat(m2): tessellate planar B-rep faces with shared edge samples"
```

---

### Task 13: 曲面网格化与 QA-10

**Files:**
- Modify: `include/geom/brep/tessellate.hpp`
- Create: `bench/bench_m2.cpp`
- Modify: `CMakeLists.txt`（`GEOM_BUILD_BENCH` 下新增 `geom_bench_m2`；定义 `GEOM_TSAN` 已存在则复用）
- Test: `tests/unit/tessellate_curved_test.cpp`
- Modify: `tests/CMakeLists.txt`

**Interfaces:**
- 扩展 `tessellate` 支持 `Cylinder3`、`Cone3`、`Sphere3`、`Torus3` 面；新增内部函数 `detail::sample_curved_face`。
- QA-10 基准入口：`bench/bench_m2.cpp` 构造 16,666 个 box（6 面 = 100k 面）并计时 `tessellate`。

- [ ] **Step 1: 写失败测试**

```cpp
TEST_CASE("cylinder mesh is watertight across the seam and caps") {
  const auto b = geom::make_cylinder({0, 0, 0}, {0, 0, 1}, 2.0, 5.0, 0.0);
  REQUIRE(b.has_value());
  const auto result = geom::tessellate(b.value(), {0.02, 0.35});
  REQUIRE(result.has_value());
  const geom::TriMesh& mesh = result.value();
  // 每条无向边恰好两次（闭合 2-manifold），且没有零面积三角形。
  std::map<std::pair<std::uint32_t, std::uint32_t>, int> use;
  for (const auto& t : mesh.triangles) {
    const geom::Point3d a = mesh.positions[t[0]];
    const geom::Point3d c = mesh.positions[t[1]];
    const geom::Point3d d = mesh.positions[t[2]];
    const geom::Vec3d n = geom::cross(c - a, d - a);
    CHECK(geom::dot(n, n) > 0.0);
    for (int i = 0; i < 3; ++i) {
      auto p = t[i];
      auto q = t[(i + 1) % 3];
      if (p > q) std::swap(p, q);
      ++use[{p, q}];
    }
  }
  for (const auto& [e, count] : use) CHECK(count == 2);
}

TEST_CASE("sphere and torus meshes have Euler characteristic 2 and 0") {
  const auto chi = [](const geom::TriMesh& m) {
    std::set<std::pair<std::uint32_t, std::uint32_t>> edges;
    for (const auto& t : m.triangles) {
      for (int i = 0; i < 3; ++i) {
        auto a = t[i];
        auto b = t[(i + 1) % 3];
        if (a > b) std::swap(a, b);
        edges.insert({a, b});
      }
    }
    return static_cast<std::int64_t>(m.positions.size()) -
           static_cast<std::int64_t>(edges.size()) +
           static_cast<std::int64_t>(m.triangles.size());
  };
  const auto s = geom::make_sphere({0, 0, 0}, 2.0, 0.0);
  REQUIRE(s.has_value());
  const auto sm = geom::tessellate(s.value(), {0.02, 0.35});
  REQUIRE(sm.has_value());
  const geom::TriMesh& sphere_mesh = sm.value();
  CHECK(sphere_mesh.positions.size() > 100);
  CHECK(chi(sphere_mesh) == 2);
  const auto t = geom::make_torus({0, 0, 0}, {0, 0, 1}, 3.0, 1.0, 0.0);
  REQUIRE(t.has_value());
  const auto tm = geom::tessellate(t.value(), {0.02, 0.35});
  REQUIRE(tm.has_value());
  CHECK(chi(tm.value()) == 0);
}

TEST_CASE("chord height stays within deflection on a circle") {
  const auto b = geom::make_cylinder({0, 0, 0}, {0, 0, 1}, 2.0, 1.0, 0.0);
  REQUIRE(b.has_value());
  const auto result = geom::tessellate(b.value(), {0.01, 0.35});
  REQUIRE(result.has_value());
  const geom::TriMesh& mesh = result.value();
  std::vector<geom::Point3d> ring;
  for (const geom::Point3d& p : mesh.positions) {
    if (std::abs(p[2]) < 1e-12) ring.push_back(p);
  }
  REQUIRE(ring.size() > 8);
  std::sort(ring.begin(), ring.end(), [](const geom::Point3d& a, const geom::Point3d& c) {
    return std::atan2(a[1], a[0]) < std::atan2(c[1], c[0]);
  });
  for (std::size_t i = 0; i < ring.size(); ++i) {
    const geom::Point3d& a = ring[i];
    const geom::Point3d& c = ring[(i + 1) % ring.size()];
    const double mx = 0.5 * (a[0] + c[0]);
    const double my = 0.5 * (a[1] + c[1]);
    CHECK(2.0 - std::hypot(mx, my) <= 0.01 + 1e-12);  // 矢高
  }
}
```

- [ ] **Step 2: 运行确认失败**

Run: 期望 `tessellate` 对圆柱返回 `unsupported`（或测试失败于水密断言）。

- [ ] **Step 3: 实现**

- **圆柱侧面**：u 采样数取与该圆柱端圆相同的 `n`（保证与端盖边界共用采样点），v 取 2 行；每对相邻 u 生成两个三角形；方向由 `du × dv` 决定。接缝列 u=0 / u=2π 必须落在同一条接缝边的两个采样上（端点已按 `VertexId` 与端盖共享），不能各自新建顶点。
- **圆锥侧面**：底圆采样 `n` 个点 + 顶点，扇形三角化；顶点只出现一次。
- **球**：u 方向 `n_u = ceil(2π / dθ)`，v 方向 `n_v = ceil(π / dθ)`；两极各自合并为单点，极点行用扇形；接缝列焊接。
- **环面**：u、v 都用 `dθ` 采样并首尾焊接；每个四边形拆两个三角形，按 `du × dv` 定向。
- 球与环面的 u/v 采样数必须由接缝边自身的采样数导出（接缝边按 T12 的 `dθ` 规则采样），保证接缝列/行逐点一致，而不是各自独立取整。
- 每个三角形的 `triangle_face` 记录来源面索引；顶点焊接统一走 T12 的 `(edge, sample)` 表，曲面内部顶点用 `(face.index(), grid_u, grid_v)` 键。

**QA-10 基准**（`bench/bench_m2.cpp`）：16,666 个 1×1×1 box 平铺（每 1000 个一组，位置错开），共 100,000 面；网格化 3 次取最小值，打印面数、三角形数、毫秒、Mtri/s；断言水密（每条无向边两次）与零翻转；输出 `QA10 faces=100000 tris=... ms=... Mtri/s=...`。目标 ≥ 1 Mtri/s，未达标时打印 `FAIL` 但进程返回 0（基准不设失败退出码，验收在文档与 CI 步骤中判定）。

- [ ] **Step 4: 运行确认通过**

Run: Debug/Release 全套 + `ctest`；`--test-case="*curved*"`、`--test-case="*sphere*"`、`--test-case="*torus*"`；`geom_bench_m2.exe` 打印 QA-10 数据。

- [ ] **Step 5: 提交**

```bash
git add include/geom/brep/tessellate.hpp bench/bench_m2.cpp CMakeLists.txt tests/unit/tessellate_curved_test.cpp tests/CMakeLists.txt
git commit -m "feat(m2): tessellate analytic curved faces and add QA-10 benchmark"
```

---

### Task 14: 质量属性、QA-11 与文档收尾

**Files:**
- Create: `include/geom/brep/mass_properties.hpp`
- Test: `tests/unit/mass_properties_test.cpp`
- Modify: `bench/bench_m2.cpp`（QA-11 误差）
- Modify: `tests/CMakeLists.txt`
- Create: `docs/benchmarks/m2.md`
- Create: `docs/adr/0021-brep-storage-and-euler-operators.md`
- Create: `docs/adr/0022-tessellation-and-mass-properties.md`
- Modify: `docs/architecture/geometry-library-architecture.md`（§16.3 状态、M2 表）
- Modify: 本计划（勾选完成步骤）

**Interfaces:**
- Produces:

```cpp
struct MassProperties {
  double volume = 0.0;
  Point3d centroid{};
  double inertia[3][3]{};   // 关于质心，单位密度
};

struct MassPropertiesOptions {
  double sampling_deflection = 1e-9;  // 平面面边界采样弦高
  int quadrature_order = 24;          // 曲面 Gauss–Legendre 每方向点数
};

Result<MassProperties> mass_properties(const Brep& brep,
                                       const MassPropertiesOptions& options = {});
```

- [x] **Step 1: 写失败测试**

```cpp
TEST_CASE("box mass properties match closed forms") {
  const auto b = geom::make_box({0, 0, 0}, {2, 3, 4}, 0.0);
  REQUIRE(b.has_value());
  const auto m = geom::mass_properties(b.value());
  REQUIRE(m.has_value());
  const geom::MassProperties& mp = m.value();
  CHECK(mp.volume == doctest::Approx(24.0).epsilon(1e-9));
  CHECK(mp.centroid[0] == doctest::Approx(1.0).epsilon(1e-9));
  CHECK(mp.centroid[1] == doctest::Approx(1.5).epsilon(1e-9));
  CHECK(mp.centroid[2] == doctest::Approx(2.0).epsilon(1e-9));
  // I_xx = m/12 (b^2 + c^2) = 24/12 (9 + 16) = 50
  CHECK(mp.inertia[0][0] == doctest::Approx(50.0).epsilon(1e-6));
}

TEST_CASE("sphere and torus volumes match closed forms") {
  const auto s = geom::make_sphere({0, 0, 0}, 3.0, 0.0);
  REQUIRE(s.has_value());
  const auto sm = geom::mass_properties(s.value());
  REQUIRE(sm.has_value());
  CHECK(sm.value().volume ==
        doctest::Approx(36.0 * 3.14159265358979323846).epsilon(1e-6));
  CHECK(sm.value().centroid[0] == doctest::Approx(0.0).epsilon(1e-9));
  CHECK(sm.value().centroid[1] == doctest::Approx(0.0).epsilon(1e-9));
  CHECK(sm.value().centroid[2] == doctest::Approx(0.0).epsilon(1e-9));

  const auto t = geom::make_torus({0, 0, 0}, {0, 0, 1}, 3.0, 1.0, 0.0);
  REQUIRE(t.has_value());
  const auto tm = geom::mass_properties(t.value());
  REQUIRE(tm.has_value());
  const double expected = 2.0 * 3.14159265358979323846 * 3.14159265358979323846 * 3.0;
  CHECK(tm.value().volume == doctest::Approx(expected).epsilon(1e-6));
  CHECK(tm.value().centroid[0] == doctest::Approx(0.0).epsilon(1e-9));
  CHECK(tm.value().centroid[1] == doctest::Approx(0.0).epsilon(1e-9));
  CHECK(tm.value().centroid[2] == doctest::Approx(0.0).epsilon(1e-9));
}

// 需要 T4 的 SquareFace / make_square_face()（复制到本文件匿名命名空间）。
TEST_CASE("open shells have no mass properties") {
  SquareFace s = make_square_face();
  const auto b = s.editor.commit();
  REQUIRE(b.has_value());
  const auto m = geom::mass_properties(b.value());
  REQUIRE_FALSE(m.has_value());
  CHECK(m.error().code == geom::Errc::not_closed);
}
```

- [x] **Step 2: 运行确认失败**

Run: 期望找不到 `geom/brep/mass_properties.hpp`。

- [ ] **Step 3: 实现**

前置：所有 shell 必须 `closed`（否则 `not_closed`）、法向一致朝外（体积必须 > 0，否则 `numeric_failure`）。

散度定理（单位密度）：
- `V = (1/3) ∮ x·n dA`
- `C = (1/(4V)) ∮ x (x·n) dA`
- `I = ∮ [(|x|² δ − x⊗x)(x·n) / 5] dA`，最后平移到质心（平行轴定理）。

平面面：把 loop 边界按 `sampling_deflection` 采样后投影到 2D，用多边形精确公式（鞋带 + 二阶矩）计算 `∮ P(x) (x·n) dA`，`n` 为常数；内环取负贡献。
曲面面：对 `(u,v)` 参数矩形做 `order × order` Gauss–Legendre（节点/权重用 `std::array` 常量表或 `gauss_legendre` 现场计算），`dA_vec = cross(du,dv) dudv`，`|x|²`/`x⊗x` 在求积点求值。曲面面假定未裁剪矩形（M2 基本体满足；裁剪面返回 `unsupported`）。

`bench/bench_m2.cpp` 增加 QA-11：对 box/cylinder/sphere/torus 计算体积、质心、惯性矩与闭式比较，打印最大相对误差；阈值 1e-6，超限打印 `FAIL`（进程仍返回 0，判定在 CI 文档）。

**文档**：
- `docs/benchmarks/m2.md`：环境、QA-10 面数/三角形数/ms/Mtri/s、QA-11 各基本体误差、是否达标、已知限制（M2 网格化只覆盖解析面）。
- ADR-0021：扁平数组 + 代际 handle + `SolidEditor` 事务 + 仅 MEV/KEV 的理由。
- ADR-0022：边采样复用契约、偏差公式、质量属性积分方法与 QA-11 解析 oracle 的替代说明。
- 架构 §16.3：M2 状态改为「代码与测试完成，QA-10/QA-11 实测见 m2.md」；若未达标必须写 FAIL，不得声称达成。

- [x] **Step 4: 运行确认通过**

Run: Debug/Release 全套 + `ctest`；`--test-case="*mass*"`；`geom_bench_m2.exe` 打印 QA-10/QA-11。

- [x] **Step 5: 提交**

```bash
git add include/geom/brep/mass_properties.hpp tests/unit/mass_properties_test.cpp bench/bench_m2.cpp tests/CMakeLists.txt docs/benchmarks/m2.md docs/adr/0021-brep-storage-and-euler-operators.md docs/adr/0022-tessellation-and-mass-properties.md docs/architecture/geometry-library-architecture.md docs/superpowers/plans/2026-11-08-m2-brep-foundation.md
git commit -m "feat(m2): compute B-rep mass properties and record QA-10/QA-11 evidence"
```

---

## M2 出口验收

- [ ] QA-10：100,000 面 B-rep 网格化——水密（每条无向边恰好两次）、无零面积/翻转三角形、弦高 ≤ deflection、共享边采样逐点一致；`docs/benchmarks/m2.md` 记录实测面数、三角形数、耗时与 Mtri/s。**（2026-11-08 实测 FAIL：99,996 面 / 199,992 三角形，0.19–0.21 Mtri/s vs 1.0 Mtri/s；水密与朝外卷绕通过，详见 `docs/benchmarks/m2.md`）**
- [x] QA-11：box/cylinder/sphere/torus 的体积、质心、惯性矩相对误差 ≤ 1e-6，记录最大值。**（实测最大相对误差 5.33e-10）**
- [x] 所有基本体 `validate_topology` 与 `validate_tolerances` 通过；Euler 特征数断言：box/cylinder/cone/sphere = 2、torus = 0。
- [ ] Debug 与 Release 全套 `ctest` 通过；GCC/Clang CI 通过；无新增警告。**（本地 MSVC Debug/Release 均 366/366 用例、450,255 断言、无新增警告；GCC/Clang CI 未在本任务中复跑）**
- [x] 架构 §16.3、`docs/benchmarks/m2.md`、ADR-0021/0022 已更新，数字与 CI 一致。

## 自检

- [ ] Spec 覆盖：§4.6（T6）、§5.7（T1/T2）、§5.8（T3）、§6.4/6.8（T3/T4）、§8.7（T4/T5/T7–T11）、§8.10（T12/T13）、§8.11（T14）、§12.3（T5/T6）、QA-10/QA-11（T13/T14）。
- [ ] 每个任务的 Interfaces 与前一任务的实际签名一致；没有「后续补上」的占位符。
- [ ] 每个 Review Focus 至少有一个对应测试：朝向（T8/T9）、容差链（T6）、接缝（T9/T11）、极点（T10/T13）、边界复用（T12/T13）。
- [ ] 计划中的代码块可直接编译（除依赖前序任务的类型外）；测试通过下一条命令即可运行。
- [ ] 没有弱断言、没有 `TODO`、没有注释掉的测试代码。

