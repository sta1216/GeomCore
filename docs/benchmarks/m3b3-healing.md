# M3b-3 Healing QA-12 证据（2026-11-13）

**结论：QA-12 正反例均通过。** 1e-6 缝隙的裂盒闭合成通过校验的实体，体积为 1
（误差 `1.11e-16`，门禁 1e-9），4 个唯一焊点、4 条缝合边，最大容差 1e-6；
1e-3 缝隙被拒绝，量值 0.001、位置 `(0, 0, 0.5005)`，输入逐字段指纹不变；
对已闭合结果再次 heal 为幂等。上述结果在 Debug 与 Release 一致。

## 环境

| 项 | 值 |
| --- | --- |
| 日期 | 2026-11-13，Asia/Shanghai |
| 系统 | Windows 10 企业版 19045 x64 |
| CPU / 内存 | Intel Core i7-8550U @ 1.80 GHz / 15.9 GB |
| 编译器 | MSVC 19.44.35217.0（VS2022 17.14），`/W4 /WX /utf-8 /fp:precise` |
| 构建 | CMake 4.4.4，`GEOM_WERROR=ON`，Release `/O2`，单线程 |
| 代码基线 | `6e370f8 test(m3b3): cover healing rejection, idempotence and ambiguity` |

## 测量定义

- **gap**：成功时 `HealReport::max_closed_gap`；拒绝时
  `HealIssue::gap`（最近候选的实际端点距离）。
- **welded / sewn**：`HealReport::welded_vertices` / `sewn_edges`；前者按唯一
  存活焊点根计数，后者是实际 `merge_edges` 的配对数。
- **tolerance**：`HealReport::max_tolerance`，即提升操作返回值与最大焊点容差的
  较大者；另检查 `validate_tolerances` 与“所有顶点 tol ≤ cap”。
- **volume / volume error**：`detail::curved_solid_volume(brep, solid, 1e-9)`
  的返回值，以及 `|volume - 1|`（裂盒解析体积为 1）。
- **rejection location**：`HealIssue::location`，即最近端点对的中点；
  **magnitude** 即同一 issue 的 `gap`。
- **fingerprint**：测试侧 `heal_fingerprint`，按索引输出每张实体表的存活标志、
  代数、几何/容差字段、loop/shell/solid 成员列表，以及 geometry store 的
  curve/surface/pcurve 清单与 pcurve 点；用于“输入不变”和幂等断言。

## QA-12 结果

| 用例 | 输入 | 结果 | gap | welded / sewn | max tol | volume error | 拒绝位置 / 量值 |
| --- | --- | --- | --- | --- | --- | --- | --- |
| QA-12 正例 | 裂盒，缝隙 1e-6，`max_gap=1e-6` | 闭合，`validate_topology` 通过 | 1e-6 | 4 / 4 | 1e-6 | 1.11e-16 | — |
| 幂等 | 对上一结果再次 heal | `changed=false`，`issues.empty()`，指纹不变 | 0 | 0 / 0 | 1e-6（不变） | — | — |
| QA-12 反例 | 裂盒，缝隙 1e-3，`max_gap=1e-6` | 拒绝，`changed=false`，指纹不变 | 0.001 | 0 / 0 | 0 | — | `(0, 0, 0.5005)` / 0.001，边 16 与 20 |

正例闭合后边界边从 8 变为 0，壳 `closed=true`，顶点/边/面计数为 12/20/10；
反例的 `gap_too_large` 是唯一 issue，`entity_b > entity_a` 且两者都在边界集合中。

## 辅助诊断（同一语料）

| 用例 | 观察 |
| --- | --- |
| 三条平行边、A 到 B/C 等距（`g=0.9e-6`，`max_gap=1e-6`） | 4 条 `ambiguous_gap`，不缝合，指纹不变 |
| 裂盒两半同向遍历（`same_direction=true`） | 4 条 `inconsistent_orientation`，不缝合 |
| 裂盒 1e-6，`max_tolerance=1e-7` | 一条 `invalid_result`，消息 `"heal_brep: tolerance cap exceeded"` |
| 裂盒 1e-6，输入链已预抬到 1e-6，`max_tolerance=1e-7` | 预检查仍拒绝；没有顶点 tol 高于 1e-7 |
| 圆柱实体 / 平面裂盒上挂 `CircleCurve3` | `Errc::unsupported`，消息 `"heal_brep: curved geometry is not supported yet"` |
| 单三角形（三条边界边，两两共享顶点） | 一条 `invalid_result`，消息 `"heal_brep: no mutual pair"` |

## 复现

```powershell
cmake --build build --config Release
build\tests\Release\geom_tests.exe --test-case="heal*" --no-skip
build\tests\Release\geom_tests.exe --test-case="heal: QA-12*" --no-skip
build\tests\Release\geom_tests.exe --count=true
```

全量回归：

```powershell
ctest --test-dir build -C Release --output-on-failure
```

Debug 配置把 `Release` 换成 `Debug` 即可；两种配置当前均为 859 个用例、
623,973 条断言全绿，`/W4 /WX /utf-8` 零警告。

## 范围与限制

- 只覆盖平面面、直线边、单实体（可多壳）的初始缝合；曲面缝隙、多实体缝合、
  小边/小面移除与曲面重拟合是 M4 范围。
- 本文是 QA-12 正确性证据，不是吞吐基准；heal 是离线操作，未记录逐次耗时。
- 容差 provenance 只到 `max_closed_gap` / `max_tolerance` 粒度，未记录每次
  放大的来源与量值（§4.6 invariant 5 的偏差，见 ADR-0026）。
