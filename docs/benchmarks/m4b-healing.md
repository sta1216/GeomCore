# M4b-1 闭合 rim 与 torus 缝隙证据（2026-11-24）

**结论：闭合边界 rim 在私有副本上自动拆分后可与普通缝隙一样愈合，torus 的开缝带
与闭合 rim 半体也在 Release 下闭合成通过校验的实体。** 圆柱闭合 rim 体积
`6.2831884487723997`（对 `π*(2+gap)` 相对 `2.5e-14`），三个 torus 开缝带为
`2.3561944901923608`、`5.4334195273316128`、`0.8635058157199863`，闭合 rim 半体为
`4.934802200544679`（`π²/2`）；圆柱与开缝带 `chi == 2`，闭合 torus 半体
**`chi == 0`**（计划与简报曾误写为 2）。本文是正确性证据，heal 是离线操作，
**不设时间门禁**。

## 环境

| 项 | 值 |
| --- | --- |
| 日期 | 2026-11-24，Asia/Shanghai |
| 系统 | Windows 10 企业版 19045 x64 |
| CPU / 内存 | Intel Core i7-8550U @ 1.80 GHz / 15.9 GB |
| 编译器 | MSVC 19.44.35217.0（VS2022 17.14），`/W4 /WX /utf-8 /fp:precise` |
| 构建 | CMake 4.4.4，`GEOM_WERROR=ON`，Release `/O2`，单线程 |
| 代码基线 | `19773bc test(heal): accept closed-rim torus halves with negatives` + Task 6 文档与拒绝锁 |

## 测量定义

- **gap**：`HealReport::max_closed_gap`，单位与模型一致。torus 夹具把目标弦长取为
  `gap*(1-1e-8)`，保证实测缝宽严格落在默认 `max_gap=1e-6` 内。
- **welded / sewn**：`HealReport::welded_vertices` / `sewn_edges`；配对焊接按
  唯一存活焊点根计数（torus 半体 4 个根即 4），点隙焊接按簇内成员计数
  （M4b-2 的针孔把 2 个顶点合成 1 个存活焊点，记 2）；`sewn_edges` 是实际
  `merge_edges` 的配对数。
- **tolerance**：`HealReport::max_tolerance`；所有接受用例都等于 `max_closed_gap`，
  不超过默认 `1e-6`。
- **V/E/F**：`brep.vertex_count()` / `edge_count()` / `face_count()`，箭头左侧为
  调用者输入、右侧为 `HealResult::brep`；括号内的“拆分后”是 `heal_brep` 私有
  工作副本的规模，即计划器实际看到的拓扑。`chi = V - E + F`（torus 半体按
  `euler_characteristic` 的 `V-E+F-(L-F)` 计算亦为 0）。
- **boundary_edges**：`HealReport::boundary_edges` 是闭合 rim 拆分后、再经点隙
  焊接后的工作边界边数（圆柱闭合 rim 4、torus 半体 8），不是输入的闭合 rim
  条数（2/4）；M4b-2 的 r0=0 针孔焊接后为 0，被 `max_gap` 拒绝时为 2。
- **collapsed_edges**：`HealReport::collapsed_edges` 是点隙焊接消解掉的边界边数
  （针孔 2 条 rim 边）；配对焊接不产生该项，M4b-1 语料恒为 0。
- **volume**：`detail::curved_solid_volume(brep, solid, 1e-9)` 的返回值；锚相对闭式误差按
  `|锚-解析| / |解析|`，torus 半体为绝对误差。
- **fingerprint**：拒绝路径用完整指纹（实体/几何存储计数、顶点位置/容差、边端点/
  曲线/参数域/容差、面曲面/反向/容差、coedge 边/面/环/反向/pcurve/next/prev/
  容差、loop 连接），断言返回 Brep 与调用者输入逐项相同。

## 验收结果

| 用例 | 输入（边界边） | V/E/F（前 → 后） | welded / sewn | max_closed_gap | 体积（验收锚） | 解析值 | 锚相对闭式误差 | chi |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| 圆柱闭合 rim | `make_cylinder_closed_rims_gap(1e-6)`，2 条闭合缝边（拆分后 4 条） | 4/6/4 → 4/6/4（拆分后 6/8/4） | 2 / 2 | `9.9999999991773336e-07` | `6.2831884487723997` | `π*(2+1e-6) = 6.2831884487722407` | 2.5e-14 | 2 |
| torus u 向开缝带 | `make_torus_u_gap(1e-6)`，2 条闭合缝边（拆分后 4 条边界边） | 8/10/4 → 6/8/4 | 2 / 2 | `9.999999901566576e-07` | `2.3561944901923608` | `3π/4 = 2.356194490192345` | 6.8e-15 | 2 |
| torus 外支 v 缝带 | `make_torus_outer_v_gap(1e-6)`，2 条闭合缝边（拆分后 4 条） | 8/10/4 → 6/8/4 | 2 / 2 | `9.9999998999999992e-07` | `5.4334195273316128` | `5.433419527331613` | 3.7e-17 | 2 |
| torus 内支跨 v 缝带 | `make_torus_inner_v_seam_gap(1e-6)`，2 条闭合缝边（拆分后 4 条） | 8/10/4 → 6/8/4 | 2 / 2 | `9.9999998992555949e-07` | `0.8635058157199863` | `0.8635058157199873` | 1.2e-15 | 2 |
| torus 闭合 rim 半体 | `make_torus_closed_rims_gap(1e-6)`，4 条闭合缝边（拆分后 8 条） | 4/6/2 → 4/6/2（拆分后 8/10/2） | 4 / 4 | `9.9999999030728907e-07` | `4.934802200544679` | `π²/2 = 4.934802200544679` | ≤1e-12（绝对） | **0** |

Measured returns from the same Release binary (`detail::curved_solid_volume(brep, solid, 1e-9)`):
cylinder `6.2831884487724015` (2.56e-14 vs `pi*(2+1e-6)`), torus u band
`2.3561944901923297` (6.41e-15), outer v band `5.433419527331667` (9.97e-15),
inner v-seam band `0.86350581571999063` (3.9e-15), closed-rim halves
`4.9348022005446222` (5.7e-14 absolute vs `pi^2/2`). The table column above is the
**acceptance anchor** used by the tests (brief quadrature reference), not the
implementation return; both stay inside the 1e-9 gate.

所有接受用例共用验收门：`report.changed == true`、`report.issues.empty()`、
`report.boundary_edges` 等于拆分后的工作计数、焊点/缝合边数符合上表、
`max_closed_gap` 等于夹具实测缝宽、输入完整指纹不变、闭合壳通过
`validate_topology`、体积在误差门内、重复运行逐位一致、二次 heal 为
`changed == false` 且 0 焊点/0 缝合边。

补充观察：

- **圆柱闭合 rim**：两条整圆缝边各在参数中点拆成半圆弧，拆分后工作副本为
  `6/8/4` 和 4 条边界边；heal 保留输入高度 `2+gap`（不收缩缝隙），表中解析值取
  `π*(2+gap)`。相对理想 `2π` 的差为 5.0e-7（门禁 1e-6）。
- **torus 开缝带**：三个夹具分别在 u 向、外支 v 向、跨 v 缝的内支切开；内支的
  lateral 反向且 `Torus3(v_origin=π)`，内支图卡 v 跨度小于 2.0。体积对独立参考值
  的相对误差最大 `6.7e-15`（门限 1e-9）。
- **torus 闭合 rim 半体**：两半体在 `z = ±h`、`h = (gap/2)*(1-1e-8)` 处切开，
  `alpha = asin(h/r)`；四段图卡窗口 `[-alpha, π-alpha]` 与
  `[π-alpha, 2π-alpha]` 各在 `1e-9` 内。`V=4, E=6, F=2` 的闭合 torus 的欧拉示性数
  是 0；把它写成 2 会把“两个圆盘盖”的拓扑误套到环面上。

## 负例（同一语料）

| 用例 | 观察 |
| --- | --- |
| 圆柱闭合 rim，`max_gap=1e-12` | 返回 `gap_too_large`；Brep 是调用者输入且完整指纹逐项不变，`changed == false`，`welded == sewn == 0`，`report.boundary_edges == 4`（拆分后） |
| 圆柱闭合 rim，`max_tolerance=1e-9` | 返回 `invalid_result "heal_brep: tolerance cap exceeded"`；同样返回输入、完整指纹不变、`changed == false`、`boundary_edges == 4` |
| torus 半体半径差 `1e-3`，默认 `max_gap=1e-6` | 不抛错，报告 `gap_too_large`（拆分后边界 8，entity_a/b = 0/1，`changed == false`） |
| 圆柱上加平面盖 | `unsupported` carrier 消息，context 为 (keep 边, drop 边) = `(1, 4)`，输入指纹不变 |
| rim 浮在半径 `1+2e-7` 上 | `unsupported "heal_brep: gap edge is not on its face surface"`，context `(2, 6)`，输入指纹不变 |
| NURBS 面/边、歧义缝 | 复用 M4a 用例；消息、context 与输入不变性均未改 |

Task 6 新增的锁测试把前两行钉在一起：同一个闭合 rim 夹具分别用
`max_gap=1e-12` 与 `max_tolerance=1e-9` 触发两类拒绝，断言两者都返回调用者输入
的完整指纹、`changed == false`、`report.boundary_edges == 4`。

## 复现

```powershell
cmake --build build-rel --config Release --target geom_tests
.\build-rel\tests\Release\geom_tests.exe --test-case="*heal*,*torus*"
.\build-rel\tests\Release\geom_tests.exe --test-case="heal: a rejected closed-rim heal returns the caller's input"
```

当前快速通道（Release）为 149 个用例、12,576 条断言、0 失败、936 跳过；其中
`*heal*` 95 用例/2,002 断言、`*torus*` 59 用例/10,754 断言。M4b-1 收尾（T5）的
基线是 `*heal*` 91/1,866、`*torus*` 59/10,754、合并 145/12,440；M4b-2 新增的
4 个 heal 用例/136 条断言来自球极点、两个锥针孔与无效解析面用例（见
`docs/benchmarks/m4b2-healing.md`）。

## 范围与限制

- 支持类别见 ADR-0026 补遗：面限 `Plane3`/`Cylinder3`/`Cone3`/`Sphere3`/`Torus3`，
  边限 `kNoGeometry`/`LineCurve3`/`CircleCurve3`；NURBS 与其他类型在修改前拒绝。
- apex/pole gap 已在 M4b-2 关闭，见 `docs/benchmarks/m4b2-healing.md`；多壳点隙
  延后。fillet/chamfer/shell 在更晚的里程碑；NURBS 曲面/曲线与曲面重拟合仍不支持。
- partial-u 不做：本里程碑总是整条闭合边在参数中点拆分，不保留部分整圆边，也不
  处理 torus 的 u 周期多圈缠绕。
- 拆分只发生在私有工作副本；拒绝时返回调用者输入，绝不把半拆分的拓扑泄漏给
  调用方。
- 本文不记录逐次耗时，也不构成吞吐基准；heal 是离线操作，无时间门禁。
- 平面 QA-12 的数值与消息未变化。契约与取舍见 ADR-0026 补遗，架构状态见 §16.4；
  M4a 的实测见 `docs/benchmarks/m4a-healing.md`。
