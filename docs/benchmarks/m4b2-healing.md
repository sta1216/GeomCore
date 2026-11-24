# M4b-2 极点与点隙愈合证据（2026-11-24）

**结论：球极点缝隙与锥顶针孔在 Release 下闭合成通过校验的实体，超出 `max_gap`
的针孔保持开口并得到诚实的 `no mutual pair`，无效解析曲面在修改前被拒绝。**
球极点三例（`1e-8/5e-9`、`5e-8/2.5e-8`、`1e-4/5e-5`）与锥针孔两例
（r0=0、r0=1e-6）体积均回到 `2π/3`；针孔把两个边界顶点焊成一个存活焊点，
`welded_vertices` 记 2，`collapsed_edges` 记 2。本文是正确性证据，heal 是离线
操作，**不设时间门禁**。

## 环境

| 项 | 值 |
| --- | --- |
| 日期 | 2026-11-24，Asia/Shanghai |
| 系统 | Windows 10 企业版 19045 x64 |
| CPU / 内存 | Intel Core i7-8550U @ 1.80 GHz / 15.9 GB |
| 编译器 | MSVC 19.44.35217.0（VS2022 17.14），`/W4 /WX /utf-8 /fp:precise` |
| 构建 | CMake 4.4.4，`GEOM_WERROR=ON`，Debug 无优化、Release `/O2`，单线程 |
| 代码基线 | `c41c6e2 fix(brep): restrict dissolve_edge to a single loop` |

## 测量定义

- **gap**：`HealReport::max_closed_gap`。极点用例的输入参数 `gap` 是两条 rim
  的纬度差，报告值是 `p0`/`c0` 两顶点的**弦距**
  `2 * sin(gap / 2)`，因此 `5e-5` 的弧长对应
  `4.9999999994791671e-05` 的弦距；针孔不产生 `max_closed_gap`。
- **welded / sewn**：`HealReport::welded_vertices` / `sewn_edges`。配对焊接按
  唯一存活焊点根计数；点隙焊接还计入**簇内全部成员**，所以 r0=0 与 r0=1e-6
  的针孔都是 2（两个边界顶点合成一个存活焊点），而 `sewn_edges == 0`。
- **boundary_edges**：`HealReport::boundary_edges` 是闭合 rim 拆分后、**再经
  点隙焊接后**的工作边界边数；针孔愈合后为 0，针孔被 `max_gap` 拒绝时为 2。
- **collapsed_edges**：`HealReport::collapsed_edges` 是点隙焊接消解掉的边界
  边数；每个针孔的两条 split rim 边都被消解，故为 2。配对焊接不产生该项。
- **tolerance**：`HealReport::max_tolerance`；点隙焊点的容差取簇直径。r0=1e-6
  的 rim 直径为 `2e-6`，测试直接断言 `max_tolerance <= 2e-6`；r0=0 的簇直径
  为零，上限 `1e-9`。把上限压到 `1e-9` 时整个计划必须被拒，见“拒绝路径”。
- **V/E/F**：`brep.vertex_count()` / `edge_count()` / `face_count()`，箭头左侧
  为输入、右侧为 `HealResult::brep`；`chi = V - E + F`。
- **volume**：`detail::curved_solid_volume(brep, solid, 1e-9)` 的返回值。解析
  锚为半球/满锥的闭式 `2π/3`，`|实测 - 锚|` 由 doctest `-s` 的断言展开读出；
  用例分别断言相对误差不超过 `1e-12`（极点与 r0=0）或 `1e-9`（r0=1e-6）。

## 解析锚与实测返回

解析锚（解析几何闭式）：半球（极点夹具）与底半径 1、高 2 的满锥（针孔夹具）
体积均为

`V = 2π/3 = 2.0943951023931953...`

下表右侧是**实测返回与解析锚的绝对差**，不是第二个解析值，也不是用锚回填的
数字；数值来自同一 Release 二进制运行 `geom_tests.exe -s` 时断言展开
`std::fabs(volume.value() - expected) <= ...` 的左端。

| 用例 | 解析锚 | 绝对差 |
| --- | --- | --- |
| 球极点 `(1e-8, 5e-9)` | `2π/3` | `1.19904e-14` |
| 球极点 `(5e-8, 2.5e-8)` | `2π/3` | `4.88498e-15` |
| 球极点 `(1e-4, 5e-5)` | `2π/3` | `1.95399e-14` |
| 锥针孔 r0=0 | `2π/3` | `1.33227e-15` |
| 锥针孔 r0=1e-6 | `2π/3` | `4.44089e-15` |

## 验收结果

### 球极点缝隙

夹具是单位球上半球（`z = 0` 圆盘封底），在北极挖掉 `[v_up, v_up + gap]` 的窄
带；带 rim 在极点下方 `delta`，帽 rim 在极点上方 `gap`，帽终止于极点。所有圆
已拆成两条半圆，因此不触发闭合 rim 拆分。三例都断言输入指纹不变、闭合单壳、
愈合后无边界边、`validate_topology` 干净、重复 heal 逐位一致。

| `(delta, gap)` | 选项 | `boundary_edges` | welded / sewn | `max_closed_gap` | 容差上限 | chi |
| --- | --- | --- | --- | --- | --- | --- |
| `(1e-8, 5e-9)` | 默认 | 4 | 2 / 2 | `5e-9` | `1e-6` | 2 |
| `(5e-8, 2.5e-8)` | 默认 | 4 | 2 / 2 | `2.5e-8` | `1e-6` | 2 |
| `(1e-4, 5e-5)` | `{1e-9, 1e-4, 1e-4}` | 4 | 2 / 2 | `4.9999999994791671e-05` | `1e-4` | 2 |

第三例的 `5e-5` 材料间隙大于默认 `max_gap=1e-6`，必须显式放宽；否则正确行为
是拒绝。它的存在同时说明极点修复不依赖“极小缝”假设。

### 锥顶针孔

锥底半径 1、半角 `atan(1/2)`（顶点在 `z_base + 2`），顶部 rim 半径 `r0`，rim
已拆成两条半圆；shell 初始开口。点隙 pass 消解两条 rim 边并把两个 rim 顶点
焊成一个存活顶点。

| 用例 | 选项 | `boundary_edges` | welded / sewn | `collapsed_edges` | V/E/F | chi | 体积相对误差 |
| --- | --- | --- | --- | --- | --- | --- | --- |
| r0=0 | `{1e-9, 1e-9, 1e-9}` | 0 | 2 / 0 | 2 | 3/3/2 | 2 | ≤ `1e-12` |
| r0=1e-6 | `{1e-9, 1e-5, 1e-4}` | 0 | 2 / 0 | 2 | 3/3/2 | 2 | ≤ `1e-9` |

r0=1e-6 的 rim 顶点间距为 `2e-6`，`max_gap=1e-5` 允许塌缩，而
`max_tolerance=1e-4` 覆盖 `2e-6` 的簇直径。两例都断言输入指纹不变、单一闭合
壳、愈合后无边界边、`validate_topology` 干净、重复 heal 逐位一致。

## 拒绝路径

| 用例 | 选项 | 结果 |
| --- | --- | --- |
| r0=1e-6 针孔，`max_gap` 低于 rim 直径 | `{1e-9, 1e-9, 1e-4}` | `changed == false`、`boundary_edges == 2`、welded/sewn/collapsed 均为 0、恰一条 `invalid_result "heal_brep: no mutual pair"`、返回 Brep 指纹等于输入、实体仍留 2 条边界边 |
| r0=1e-6 针孔，`max_tolerance` 低于簇直径 | `{1e-9, 1e-5, 1e-9}` | `changed == false`、`boundary_edges == 2`、welded/collapsed 均为 0、`max_tolerance == 0`、恰一条 `invalid_result "heal_brep: tolerance cap exceeded"`、返回 Brep 指纹等于输入 |
| 无效球面 frame（axis 与 x_axis 重合） | 默认 | `Err(unsupported, "heal_brep: surface is invalid")`、`context_a == 1`（**曲面存储索引**，不是面索引）、完整指纹不变、重复调用结果一致 |

无效球面例用 `make_malformed_sphere_pole_gap(1e-6, 5e-7, true)`：前面加一个有效
`Plane3` 把坏球推到曲面索引 1，而首个面仍坐在该球上，因此用面索引作 context
的实现会失败。校验只覆盖提供 `is_valid()` 的
`Cylinder3`/`Cone3`/`Sphere3`/`Torus3`；`Plane3` 没有 `is_valid()`，平面不做
该校验。

## 计数与安全泳道

- Release 快速通道 `-tc="*heal*"`：95 用例 / 2,013 断言，0 失败；其中新用例
  `heal: an invalid analytic surface names the surface index` 为 15 条断言。
- Release UV 泳道 `-tc="*uv*"`：9 用例 / 388 断言，0 失败。
- Debug 全量：**未完成（PENDING，移交里程碑）**。本机 `/Od` 下球体布尔用例极
  慢，分片运行 5 小时完成 1,077/1,085 用例、553,598 断言、0 失败；未跑完的 8 个
  用例全在 `brep_boolean_curved_test.cpp`（球带、球冠、补集、平板与 torus 半空间）。
  Debug 安全泳道 `-tc="*heal*"` 为 95/2,013、`-tc="*uv*"` 为 9/388，0 失败。
- 全量 Release：1085 用例 / 680,739 断言，0 失败（1074→1085 用例、680,503→
  680,739 断言，+236 条）。

Debug/Release 均为强制全量重编译后的测试二进制（`geom_tests`），不接受增量
对象文件。
