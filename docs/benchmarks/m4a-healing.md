# M4a 曲面缝隙愈合证据（2026-11-24）

**结论：三个解析曲面缝隙在 Release 下均闭合成通过校验的实体。** 圆柱（缝隙
`1e-6`）体积 `6.2831884487723997`，球（赤道分离 `1e-6`）`4.1887902047863808`，
锥（带终止于顶点，斜缝 `8.9442719102563569e-07`）`2.0943951023931864`；三者
`chi == 2`、`validate_topology` 干净、二次 heal 幂等、重复运行逐位一致。本文是
正确性证据，heal 是离线操作，**不设时间门禁**。

## 环境

| 项 | 值 |
| --- | --- |
| 日期 | 2026-11-24，Asia/Shanghai |
| 系统 | Windows 10 企业版 19045 x64 |
| CPU / 内存 | Intel Core i7-8550U @ 1.80 GHz / 15.9 GB |
| 编译器 | MSVC 19.44.35217.0（VS2022 17.14），`/W4 /WX /utf-8 /fp:precise` |
| 构建 | CMake 4.4.4，`GEOM_WERROR=ON`，Release `/O2`，单线程 |
| 代码基线 | `e0fcc6e docs(heal): correct the sphere uv-area rationale` + Task 6 闭合 rim 诊断 |

## 测量定义

- **gap**：`HealReport::max_closed_gap`，单位与模型一致。
- **welded / sewn**：`HealReport::welded_vertices` / `sewn_edges`；M4a 语料只有
  配对焊接，`welded_vertices` 按唯一存活焊点根计数；点隙焊接的簇成员计数与
  `collapsed_edges` 见 M4b-2，M4a 语料不涉及。
- **tolerance**：`HealReport::max_tolerance`；三个接受用例都等于 `max_closed_gap`，
  且不超过默认 `1e-6`。
- **V/E/F**：`brep.vertex_count()` / `edge_count()` / `face_count()`，箭头左侧为
  输入、右侧为 `HealResult::brep`。`chi = V - E + F`。
- **volume**：`detail::curved_solid_volume(brep, solid, 1e-9)` 的返回值；解析值用
  同一 `gap` 的闭式体积。相对误差按 `|实测-解析| / |解析|`。
- **fingerprint**：测试侧 `brep_fingerprint`（实体计数 + 每个顶点 xyz），用于
  断言拒绝路径输入不变；接受路径另断言重复调用与二次 heal 的结果一致。

## 验收结果

| 用例 | 输入（边界边） | V/E/F（前 → 后） | welded / sewn | max_closed_gap | 体积（实测） | 解析值 | 相对误差 | chi |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| 圆柱 | `make_two_band_cylinder(1, 1e-6)`，4 条 | 8/10/4 → 6/8/4 | 2 / 2 | `9.9999999991773336e-07` | `6.2831884487723997` | `π*(2+1e-6) = 6.2831884487722407` | 2.5e-14 | 2 |
| 球（四片 lune） | `make_sphere_gap(5e-7)`，赤道分离 `1e-6`，4 条 | 6/8/4 → 4/6/4 | 2 / 2 | `9.9999999999999995e-07` | `4.1887902047863808` | `4π/3 = 4.1887902047863905` | 2.3e-15 | 2 |
| 锥（带终止于顶点） | `make_cone_gap(8e-7)`，轴向 `8e-7` / 斜缝 `8.9442719102563569e-07`，4 条 | 7/8/3 → 5/6/3 | 2 / 2 | `8.9442719102563569e-07` | `2.0943951023931864` | `2π/3 = 2.0943951023931953` | 4.2e-15 | 2 |

三个用例共用同一验收门（`check_closed_curved_heal`）：`report.changed == true`、
`report.issues.empty()`、`report.boundary_edges == 4`、`welded_vertices == 2`、
`sewn_edges == 2`、`max_closed_gap` 等于夹具实测缝宽、`max_tolerance <= 1e-6`、
输入指纹不变、闭合壳通过 `validate_topology`、`chi == 2`、体积在相对门内、重复
运行逐位一致、二次 heal 为 `changed == false` 且 0 焊点/0 缝合边。

补充观察：

- **圆柱**：heal 保留输入高度 `2+gap`（不收缩缝隙），所以相对理想 `2π` 的差为
  5.0e-7（门禁 1e-6）；表中解析值取 `π*(2+gap)`，相对误差 2.5e-14。切向半径
  `2π` 的检查用 1e-6 相对门。
- **球**：四片 lune 的 uv 面积各自为正，总和为 `2π²`（1e-5 内）；赤道处
  `4π*asin(eps)` 的球面 sliver 是面积亏欠，不是 uv 图卡亏欠。焊接后 4 顶点、
  6 边、4 面。
- **锥**：7 个输入顶点减去 2 个焊接顶点得 5 个，恰好 1 个顶点在
  `(0, 0, 2)`；锥顶不展开成扇形。斜缝 `8.944e-7` 留在默认 `max_gap` /
  `max_tolerance = 1e-6` 内。

## 负例（同一语料）

| 用例 | 观察 |
| --- | --- |
| 闭合整圆 rim（`make_cylinder_closed_rims_gap(1e-6)`） | `Errc::unsupported`，消息 `"heal_brep: closed gap edge cannot be paired; split the rim first"`，context `(1, 0)`，输入指纹不变；旧树上是 `invalid_topology "heal_plan: boundary edge has zero length"` |
| 半径差 `2e-7` / 圆柱上加平面盖 / rim 浮在半径 `1+2e-7` 上 | `unsupported`，分别报 carrier / edge-on-surface 固定消息，context 为 (keep 边, drop 边)，输入指纹不变 |
| 继承曲线不穿过焊后端点 | `unsupported "heal_brep: kept edge curve does not match its welded vertices"` |
| 半径差 `1e-3`，默认 `max_gap=1e-6` | 不抛错，报告 `gap_too_large`（`changed == false`，边 2/4） |

## 复现

```powershell
cmake --build build-rel --config Release --target geom_tests
.\build-rel\tests\Release\geom_tests.exe --test-case="heal: analytic curved gaps are accepted"
.\build-rel\tests\Release\geom_tests.exe --test-case="heal:*,curved volume:*,tessellate:*,validate:*"
```

当前目标过滤（Release）为 68 个用例、4,274 条断言、0 失败、1,000 跳过；闭合 rim
用例单独运行 1 个用例、6 条断言全绿。全量 Debug/Release 只在里程碑评审运行。

## 范围与限制

- 支持类别见 ADR-0026 补遗：面限 `Plane3`/`Cylinder3`/`Cone3`/`Sphere3`/`Torus3`，
  边限 `kNoGeometry`/`LineCurve3`/`CircleCurve3`；NURBS 与其他类型在修改前拒绝。
- torus 缝隙愈合（v 周期、双半径）延后到 M4b，本里程碑不声称支持；`Torus3`
  仅通过面类型白名单。
- 顶点/极点缩成一点的 apex/pole gap 不在范围内；锥用例是“带终止于顶点”，边界
  仍是有限半径圆。
- 本文不记录逐次耗时，也不构成吞吐基准；heal 是离线操作，无时间门禁。
- 平面 QA-12 的数值与消息未变化（裂盒体积 1、误差 `1.11e-16`、4/4 焊点/缝合边、
  1e-3 拒绝于 `(0,0,0.5005)`）。契约与取舍见 ADR-0026 补遗，架构状态见 §16.4。
