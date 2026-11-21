# M3c-3c-4b3c-2：锥 apex 环闭包与过轴平面切 Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use subagent-driven-development (recommended) or executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** 修复锥面修剪环在 apex 处的欠覆盖，并让整锥对「包含锥轴的平面」完成差/交/并布尔。

**Architecture:** 三处定点修复 + 一处文档：①`trim_face_split_domain` 的 `touches_pole` 分支在环沿塌缩行闭合时补上 front 的 apex 参数副本；②`split_face_curved` 的圆弧径向投影若落在既有 loop 顶点 `snap_radius` 内，复用该顶点；③布尔驱动的锥/平面 gate 接受过轴平面、按 AABB 跳过不相交的锥/平面对、并把 apex 参数副本排除在「seam copies coincide」守卫之外。SSI 已经返回全部 6 条切链，不新增解析切链。

**Tech Stack:** C++20；CMake 4.4.4 + VS2022（MSVC 19.44）；doctest 2.4.11；`/W4 /WX /utf-8`。

**Spec:** `.superpowers/sdd/m3c3c4b3c2-probe-report.md`（探针报告，含每项测量值与原型）；架构 §8.10/§8.12、ADR-0037。

## Global Constraints

- 不新增依赖；不改公共 API（`Torus3`/`Cone3`/`split_face_curved`/`boolean_solid_curved` 签名不变）。
- 拓扑判定只用精确谓词；浮点只用于采样、容差匹配与拒绝预处理。
- 不得静默返回错误实体：拒绝必须是精确 `Errc` + 消息；`failure_report` 在拒绝前记录。
- 锥/球/环/柱/盒的既有体积、计数、逐位确定性不得回归；除 Task 1 明确列出的 5 个 pin 外，不改既有测试期望。
- 每任务 TDD：先改测试→构建确认 RED→最小实现→Debug/Release 全量→单 commit。
- 改头文件后 MSBuild 头依赖失效：`touch` 任一 `.cpp` 或 `/t:Rebuild` 后再跑测试。
- 测试放在既有文件（`tests/unit/brep_trim_face_test.cpp`、`tests/unit/brep_boolean_curved_test.cpp`），已注册，无需改 CMake。

## Review Focus

1. **过轴但只擦到 apex 的平面**：必须响亮拒绝，绝不能返回半个锥体。Task 3 的 apex-only 用例钉住。
2. **切割面只覆盖部分直径**（box 面比锥截面窄）：切链在参数域内终止，应报 `"trim_face: cut chain ends inside the parameter domain"`，不得半装配。Task 3 用例钉住。
3. **过轴平面同时包含一条缝射线**（如 `y >= 0`）：环自触，应报 `"boolean_solid_curved: result touches itself along an edge"`。Task 3 用例钉住。
4. **旋转/平移的锥**：方位角与 apex 必须用 `Cone3::axis/x_axis/y_axis`，不能用世界 z；Task 3 用 +x 轴与 `(0,5,0)` 平移各跑一次差集。
5. **apex 参数副本 vs 真正的缝重合**：前者必须焊接后成功；后者（平移锥）仍拒绝。Task 3 保留既有拒绝用例并新增 apex 用例。

---

## File Map

| 文件 | 责任 |
|---|---|
| `include/geom/brep/trim_face.hpp` | Task 1：`trim_face_split_domain` 的 pole 环闭合；Task 2：圆弧径向投影后的顶点复用 |
| `include/geom/brep/boolean_solid_curved.hpp` | Task 3：锥/平面 gate、AABB 过滤、apex 副本豁免 |
| `tests/unit/brep_trim_face_test.cpp` | Task 1 的 3 个 pin 更新 + 1 个新用例；Task 2 的 1 个新用例 |
| `tests/unit/brep_boolean_curved_test.cpp` | Task 3 的过轴成功用例、谓词 pin、4 个拒绝用例 |
| `docs/adr/0038-cone-axis-plane-cuts.md` | Task 4：契约、证据、延期项 |
| `docs/benchmarks/m3c3c4b3c2.md` | Task 4：实测矩阵 |
| `docs/architecture/geometry-library-architecture.md` | Task 4：§16.3 段与计数 |

---

### Task 1: 锥 apex 环沿塌缩行闭合

**Files:**
- Modify: `include/geom/brep/trim_face.hpp:1774-1785`
- Test: `tests/unit/brep_trim_face_test.cpp:776-824`（cap）、`:826-935`（meridian）、`:1064`（QA-21）

**Interfaces:**
- Consumes: `detail::trim_face_seam_only_near_side(double value, double side) -> bool`（`trim_face.hpp:767-769`，1e-12 相对）。
- Produces: pole-touching run 的环在起点位于底边（`v == v_lo`）时，末尾依次出现 `(u_back, v_apex)`、`(u_front, v_apex)`、再重复首点；环的所有显式边轴对齐，uv 面积恰好是扇形面积。

- [ ] **Step 1: 更新 3 个失效 pin 并新增扇形用例**

在 `tests/unit/brep_trim_face_test.cpp`：

1. `"trim_face: a closed horizontal cut splits a cone cap from a band"`：cap 面积由 π 改为 `2*kPi`，`apex_points` 由 1 改为 2；band 仍 `2*kPi`；断言两片绝对面积和 `4*kPi`（`epsilon(1e-9)`）。
2. `"trim_face: a meridian open cut splits the cone chart in two"`：`areas[0] == 1.4`、`apex_points[0] == 2`、`areas[1] == 2.0*(2.0*kPi-0.7)`，和 `4*kPi`；注释删去「half its chart area because the pole closure replaces the collapsed top edge with a chord」。
3. `"QA-21: a cone meridian split is deterministic and closes at the apex"`：同样的 1.4 / 2 修正。
4. 新增：

```cpp
TEST_CASE("trim_face: the apex ring of a seam-adjacent wedge is the full sector") {
  const Brep cone = make_cone(Point3d{0,0,0}, Vec3d{0,0,1}, 1.0, 2.0, kTol).value();
  const FaceId lateral = face_with_surface<Cone3>(cone);
  const double v_apex = 1.0 / std::tan(std::atan2(1.0, 2.0));
  const auto split = [&](double u_cut) {
    std::vector<Point2d> meridian{{u_cut, 0.0}, {u_cut, v_apex}};
    return split_face_curved(cone, lateral, std::vector<FaceCut>{FaceCut{meridian}}, 1e-9);
  };
  // u = pi/2: seam-adjacent wedge, 19-point ring ending with both apex copies.
  const auto half = split(kPi / 2.0);
  REQUIRE(half.has_value());
  REQUIRE(half.value().size() == 2u);
  const std::vector<Point2d>& wedge = half.value()[0].loops_uv.front();
  REQUIRE(wedge.size() == 19u);
  CHECK(wedge[17][0] == doctest::Approx(kPi / 2.0));
  CHECK(wedge[17][1] == doctest::Approx(v_apex));
  CHECK(wedge[18][0] == doctest::Approx(0.0));
  CHECK(wedge[18][1] == doctest::Approx(v_apex));
  CHECK(std::fabs(signed_area2(wedge)) == doctest::Approx(kPi).epsilon(1e-9));
  // u = 0.7: the same closure at a non-side sample.
  const auto small = split(0.7);
  REQUIRE(small.has_value());
  REQUIRE(small.value().size() == 2u);
  const std::vector<Point2d>& sector = small.value()[0].loops_uv.front();
  REQUIRE(sector.size() == 11u);
  CHECK(sector[9][0] == doctest::Approx(0.7));
  CHECK(sector[9][1] == doctest::Approx(v_apex));
  CHECK(sector[10][0] == doctest::Approx(0.0));
  CHECK(sector[10][1] == doctest::Approx(v_apex));
  CHECK(std::fabs(signed_area2(sector)) == doctest::Approx(1.4).epsilon(1e-9));
  for (const auto& result : {half, small}) {
    double sum = 0.0;
    for (const CurvedFragment& fragment : result.value()) {
      sum += std::fabs(signed_area2(fragment.loops_uv.front()));
    }
    CHECK(sum == doctest::Approx(4.0 * kPi).epsilon(1e-9));
  }
}
```

- [ ] **Step 2: 运行确认 RED**

Run: `.\build\tests\Debug\geom_tests.exe -tc=*cone*`

Expected: 5 个用例失败——上面 3 个被改的 pin、新扇形用例，以及 cap 用例；失败值为 `0.7`/`π/2` 与 `1.4`/`π` 的差。

- [ ] **Step 3: 在 `touches_pole` 分支补 front 的 apex 副本**

`include/geom/brep/trim_face.hpp` 的

```cpp
        points.push_back(Point2d{points.back()[0], pole_v});
        points.push_back(points.front());
```

之间插入：

```cpp
        if (!trim_face_seam_only_near_side(points.front()[1], pole_v)) {
          // The run started on the opposite side row, so the ring must still
          // travel the whole collapsed row to the front's meridian: the row is
          // one 3D point, and its copy at the front's u carries that travel.
          // Without it the ring closes with a uv diagonal from the pole to the
          // front, which under-covers the chart by the front's u extent and
          // leaves the pcurve off the ring.
          points.push_back(Point2d{points.front()[0], pole_v});
        }
```

- [ ] **Step 4: 运行确认 GREEN（定向 + 两配置全量）**

Run: `.\build\tests\Debug\geom_tests.exe -tc=*cone*` → 46 用例全过。
Run: Debug 与 Release 全量 → 记录实际用例/断言数（基线 1031/674,604）。

- [ ] **Step 5: 提交**

```bash
git add include/geom/brep/trim_face.hpp tests/unit/brep_trim_face_test.cpp
git commit -m "fix(trim_face): close the apex ring along the collapsed row"
```

---

### Task 2: 圆弧径向投影后复用既有 rim 顶点

**Files:**
- Modify: `include/geom/brep/trim_face.hpp`（`split_face_curved` 的 arc 径向投影分支，基提交约 `:2668`）
- Test: `tests/unit/brep_trim_face_test.cpp`（文件末尾新增）

**Interfaces:**
- Produces: `namespace detail { inline bool trim_face_reuse_loop_vertex(std::span<const std::vector<Point2d>> loops, Point2d& point, double radius); }`——若 `loops` 中存在与 `point` 距离 ≤ `radius` 的顶点，把 `point` 改为该顶点并返回 true；否则返回 false 且不改 `point`。

- [ ] **Step 1: 写失败测试（底圆盘直径切）**

```cpp
TEST_CASE("trim_face: a diameter cut on the base disk reuses the rim vertex") {
  const Brep cone = make_cone(Point3d{0,0,0}, Vec3d{0,0,1}, 1.0, 2.0, kTol).value();
  const FaceId base = face_with_surface<Plane3>(cone);
  // Base chart uv = (y, x); the cut is the diameter y in [-1, 1] at x = 0.
  const std::vector<Point2d> diameter{{-1.0, 0.0}, {1.0, 0.0}};
  const auto split = split_face_curved(cone, base, std::vector<FaceCut>{FaceCut{diameter}}, 1e-9);
  REQUIRE(split.has_value());
  REQUIRE(split.value().size() == 2u);
  double area = 0.0;
  for (const CurvedFragment& fragment : split.value()) {
    REQUIRE(fragment.loops_uv.size() == 1u);
    const std::vector<Point2d>& ring = fragment.loops_uv.front();
    area += std::fabs(signed_area2(ring));
    for (std::size_t i = 0; i < ring.size(); ++i) {
      const Point2d& a = ring[i];
      const Point2d& b = ring[(i + 1) % ring.size()];
      CHECK(std::hypot(a[0] - b[0], a[1] - b[1]) > 1e-12);
    }
  }
  CHECK(area == doctest::Approx(kPi).epsilon(1e-9));
}
```

- [ ] **Step 2: 运行确认 RED**

Run: `.\build\tests\Debug\geom_tests.exe -tc="trim_face: a diameter cut*"`

Expected: FAIL——环里有一对相邻点相距约 `1.8e-16`（径向投影落在 `(cos(-pi/2), sin(-pi/2))` 的末位，而 loop 顶点是精确 `(-1,0)`）。

- [ ] **Step 3: 实现 helper 并在 arc 分支调用**

在 `trim_face.hpp` 的 detail 区（`trim_face_find_snap` 附近）加：

```cpp
inline bool trim_face_reuse_loop_vertex(std::span<const std::vector<Point2d>> loops,
                                        Point2d& point, double radius) {
  for (const std::vector<Point2d>& loop : loops) {
    for (const Point2d& vertex : loop) {
      if (std::hypot(vertex[0] - point[0], vertex[1] - point[1]) <= radius) {
        point = vertex;
        return true;
      }
    }
  }
  return false;
}
```

在 `split_face_curved` 的 arc 径向投影分支里，投影得到 `snap.point` 后调用：

```cpp
if (trim_face_reuse_loop_vertex(loops, snap.point, snap_radius)) {
  snap.vertex = true;
}
```

仅在该分支调用（不要在通用边界 snap 里用），避免把普通边界吸附改成顶点吸附。

- [ ] **Step 4: 运行确认 GREEN（定向 + 两配置全量）**

Run: 定向用例 1 个通过；Debug/Release 全量通过，记录计数。

- [ ] **Step 5: 提交**

```bash
git add include/geom/brep/trim_face.hpp tests/unit/brep_trim_face_test.cpp
git commit -m "fix(trim_face): reuse the rim vertex after radial projection"
```

---

### Task 3: 接受并装配过轴锥/平面切

**Files:**
- Modify: `include/geom/brep/boolean_solid_curved.hpp`（谓词基提交 `:1426-1445`；锥/平面对循环 `:2093-2126`；apex 副本守卫 `:2654-2668`）
- Test: `tests/unit/brep_boolean_curved_test.cpp`（替换 `:932` 的拒绝用例、更新谓词 pin `:955-993`、新增负例）

**Interfaces:**
- Produces: `namespace detail { inline bool curved_boolean_cone_plane_contains_axis(const Cone3& cone, const Plane3& plane, double tolerance); }`——`|dot(unit_normal, axis)| <= kCurvedBooleanConePerpendicularSlack` 且 `|dot(unit_normal, origin) - offset/|normal|| <= tolerance`。
- Consumes: Task 1 的 apex 双副本环；Task 2 的底圆盘顶点复用。

- [ ] **Step 1: 写失败测试**

1. 谓词 pin：`supported(cone, Plane3{{1,0,0}, 0.0}) == true`（过轴）；`Plane3{{1,0,0}, 0.5} == false`（平行偏移、双曲线族）；`Plane3{{1,0,0}, 2.0} == true`（远到错过保守包围球）；既有垂直/斜切断言不变。
2. 替换 `"curved boolean: a cone refuses an axis-containing plane cut"`：

```cpp
TEST_CASE("curved boolean: an axis-containing plane cuts the cone in two rulings") {
  const Brep cone = make_cone(Point3d{0,0,0}, Vec3d{0,0,1}, 1.0, 2.0, kTol).value();
  const Brep box = make_box(Point3d{0,-1.5,-0.5}, Point3d{1.5,1.5,2.5}, kTol).value();
  const auto run = [&](SolidBooleanOp op) {
    return boolean_solid_curved(cone, cone.solid_ids().front(), box, box.solid_ids().front(),
                                op, CurvedBooleanOptions{1e-9});
  };
  const auto difference = run(SolidBooleanOp::difference);
  const auto intersection = run(SolidBooleanOp::intersection);
  const auto united = run(SolidBooleanOp::union);
  REQUIRE(difference.has_value());
  REQUIRE(intersection.has_value());
  REQUIRE(united.has_value());
  const auto check = [](const CurvedBooleanResult& result, double volume,
                        std::size_t faces) {
    CHECK(validate_topology(result.brep).ok());
    CHECK(result.brep.face_ids().size() == faces);
    CHECK(std::fabs(volume_of(result.brep) - volume) <= 1e-12);
    CHECK(result.brep.solid_shells(only_solid(result.brep)).size() == 1u);
  };
  check(difference.value(), kPi / 3.0, 3u);
  check(intersection.value(), kPi / 3.0, 4u);
  check(united.value(), 13.5 + kPi / 3.0, 8u);
  // Determinism: the difference repeats bitwise.
  const auto again = run(SolidBooleanOp::difference);
  REQUIRE(again.has_value());
  CHECK(volume_of(again.value().brep) == volume_of(difference.value().brep));
}
```

（`volume_of`、`only_solid`、`kTol` 已在测试文件中；若 `CurvedBooleanResult` 字段名不同，以文件既有用法为准。）

3. 放置与拒绝用例：
   - 锥旋转到 +x（`make_cone(origin, Vec3d{1,0,0}, 1.0, 2.0, kTol)`）对 `box y>=0`：difference 体积 `π/3`。
   - 锥平移 `(0,5,0)`：difference 体积 `π/3`。
   - apex-only：`box {0,-1.5,1.5}..{1.5,1.5,2.5}` → `Errc::unsupported`，无实体。
   - 切割面窄于直径：`box {0,-1.5,0.2}..{1.5,1.5,1.0}` → `unsupported`，消息 `"trim_face: cut chain ends inside the parameter domain"`。
   - 切割面含缝射线：`box {0,0,-0.5}..{1.5,1.5,2.5}` → `unsupported`，消息 `"boolean_solid_curved: result touches itself along an edge"`。
   - 保留 `"curved boolean: a translated cone refuses coincident seam copies"` 不变。

- [ ] **Step 2: 运行确认 RED**

Run: `.\build\tests\Debug\geom_tests.exe -tc="curved boolean: a*n axis-containing*"`

Expected: FAIL——当前返回 `unsupported "cone cutters must be perpendicular..."`。

- [ ] **Step 3: 三处最小修改**

1. 谓词 `curved_boolean_cone_plane_cut_supported`：在垂直分支之后、包围球分支之前：

```cpp
  if (curved_boolean_cone_plane_contains_axis(cone, plane, tolerance)) return true;
```

2. 锥/平面对循环内，在调用谓词前跳过不相交的盒：

```cpp
      if (!curved_boolean_bounds_overlap(face_a.bounds, face_b.bounds, tolerance)) continue;
```

3. apex 副本守卫：若两个逐位相同的点都在 `apex = origin + axis * (radius / tan(half_angle))` 的 `tolerance` 内，`continue`（不记录失败）；真正的缝重合仍拒绝。

- [ ] **Step 4: 运行确认 GREEN（定向 + 两配置全量）**

Run: 定向 `-tc="curved boolean:*"` 与 `-tc=*cone*`；Debug/Release 全量 → 记录计数。

- [ ] **Step 5: 提交**

```bash
git add include/geom/brep/boolean_solid_curved.hpp tests/unit/brep_boolean_curved_test.cpp
git commit -m "feat(boolean): cut cones with planes through the axis"
```

---

### Task 4: 契约文档与基准

**Files:**
- Create: `docs/adr/0038-cone-axis-plane-cuts.md`
- Create: `docs/benchmarks/m3c3c4b3c2.md`
- Modify: `docs/architecture/geometry-library-architecture.md`（§16.3 末尾加 M3c-3c-4b3c-2 段；ADR-0037 的 cone 延期句改为指向 ADR-0038）
- Modify: `docs/superpowers/plans/2026-11-22-m3c3c4b3c2-cone-axis-planes.md`（勾选自检）

**Interfaces:**
- 无代码接口；文档必须与 Task 3 的实测值一致。

- [ ] **Step 1: 写 ADR-0038（≤120 行）**

记录：过轴平面切的两条经线由既有 SSI 提供（零合成链）；gate 的两条新条件；AABB 过滤；apex 副本豁免；apex 环闭包修复；验收矩阵（3 个运算、旋转/平移、拒绝类）；延期项（反序差集、偏移平行/斜切、切割面窄于直径的参数域裁剪、apex-only 专用诊断、含缝射线、非整锥、同面双 pole-run 合并）。

- [ ] **Step 2: 写 benchmark 与 §16.3 段**

`docs/benchmarks/m3c3c4b3c2.md`：测试矩阵（difference 1.0471975511965941、intersection 1.0471975511965979、union 14.547197551196593；faces 3/4/8；edges 64/66/76；vertices 63/64/71；0 zero-length；逐位确定性）、与 π/3 的相对误差、运行命令、环境。

- [ ] **Step 3: 全量两配置并统一计数**

Run: Debug/Release 全量；把实际 `N cases / M assertions` 写进 ADR、bench、§16.3。

- [ ] **Step 4: 提交**

```bash
git add docs/adr/0038-cone-axis-plane-cuts.md docs/benchmarks/m3c3c4b3c2.md \
        docs/architecture/geometry-library-architecture.md \
        docs/superpowers/plans/2026-11-22-m3c3c4b3c2-cone-axis-planes.md
git commit -m "docs(m3c3c4b3c2): record the axis-plane cone contract"
```

---

## Exit Criteria

- `-tc=*cone*` 46 用例全过；`-tc="curved boolean:*"` 全过；Debug/Release 全量 0 失败、0 警告。
- 过轴差/交/并体积与解析值相对误差 ≤1e-12；单壳、`validate_topology` 通过、0 零长 3D 边、两次运行逐位相同。
- 旋转 +x、平移 (0,5,0) 的差集同为 `π/3`。
- apex-only、窄切割面、含缝射线、平移锥缝重合全部确定性拒绝，无错误实体。
- 除 Task 1 列出的 3 个 pin 外无既有测试期望改动；`git diff --check` 干净。

## Self-Review

- [ ] 每个 Task 的测试名、数值、文件行号与探针报告一致。
- [ ] Task 3 依赖 Task 1/2；三者按序执行，不并行改同一头文件。
- [ ] Review Focus 5 条各自有任务测试（1→Task 3 apex-only，2→Task 3 窄切割面，3→Task 3 含缝射线，4→Task 3 旋转/平移，5→Task 3 apex 用例 + 既有平移锥拒绝）。
- [ ] 无占位符、无未定义 helper、无跨任务类型漂移。
- [ ] 允许的 pin 变更只有 3 个；文档计数在 Task 4 填入实测值。
