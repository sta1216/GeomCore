# M3c-3c-1 一般位置球/球透镜：QA-17 与实测

本文件记录 M3c-3c-1 的验收数据：QA-17（37° 放置）与 15° 已知难例的体积/拓扑实测、
33 次旋转语料与全套的墙钟耗时，以及该里程碑的适用范围。取舍与设计理由见
[ADR-0031](../adr/0031-general-sphere-lenses.md)。

本任务**不新增 bench 可执行文件**（计划未要求）。下面的耗时全部来自 doctest
聚焦运行的墙钟测量，不是统计基准；它们只用来记录 QA 用例与语料的实际量级。

## 环境

| 项 | 值 |
| --- | --- |
| CPU / 内存 | Intel Core i7-8550U（4C/8T）/ 15.9 GB |
| 操作系统 | Windows 10 19045 (x64，企业版) |
| 编译器 | MSVC 19.44.35217.0（`/W4 /WX /utf-8`，GEOM_WERROR=ON） |
| CMake | 4.4.4，Visual Studio 17 2022 generator |
| 构建配置 | Debug（`/Od`）与 Release（`/O2`），两套树分别为 `build`、`build-rel` |
| 容差 | `1e-9`（`CurvedBooleanOptions::tolerance`） |
| 计时工具 | PowerShell 5.1 `Measure-Command`（墙钟，含进程启动） |

## 语料

- **QA-17**（`tests/unit/curved_boolean_qa17_test.cpp`）：两半径 1 的整球，球心
  `(0,0,0)` 与 `1.5 * (cos37°, sin37°, 0)`，即 `d = 1.5`、中心线方向 37°
  （不在 15°/30°/10° 任何网格上）。交集，2 次布尔调用，10 条断言。
- **15° 已知难例**（`tests/unit/brep_boolean_curved_test.cpp` 的
  `curved boolean: an off-axis sphere lens is built`）：同一 `d = 1.5`、方向 15°，
  交集与并集各一次，21 条断言（含两个操作数的 `validate_topology` 与共享 rim）。
- **旋转语料**（`curved boolean: sphere lenses accept rotated placements`）：
  30° 全圈交集 12 次、0°/45°/90° 的并/差 6 次、0–90° 的 10° 细网格交集 10 次、
  `mt19937` 固定种子 5 个随机抽样，共 33 次布尔调用、396 条断言。
- 结构判据：QA-17 断言 `validate_topology().ok()`、恰好 1 个 shell、chart 体积在
  `1e-6` 相对内等于解析值、同参数两次运行体积逐位相同；15° 与旋转语料额外断言
  交/并/差的解析体积、恰好一条被两个 coedge 反向使用的共享 rim 边（每条 6 条
  断言）。

## 口径

- **体积度量**：lens rim 在 `1e-9` 下约 4098 个 uv 点，`detail::curved_solid_volume`
  的 ear clipping 在数千点轮廓上不收敛（M3c-3b Task 4 观察），因此 QA-17 与既有
  lens 用例一样用测试共享头
  `tests/unit/brep_boolean_curved_fixture.hpp` 的 `sphere_face_volume`：对每个
  `Sphere3` 面取 `face_loops_uv` 的环，用 shoelace 定符号、从 `±π/2` 极点求
  cap 高 `h`，累加 `sign*πh²(3r−h)/3`。
- **计时口径**：每个聚焦用例在无其他构建/测试负载时连跑 3 次，报中位数与
  min–max；Debug 的 33 次语料只跑 1 次（175 s 量级，取单次墙钟）。墙钟包含
  进程启动（Release 实测约 `0.05 s`，Debug `0.05–0.20 s`），因此每个用例的
  "单次布尔"是 `(墙钟 − 启动) / 布尔次数` 的近似。计时是 informational，
  没有阈值，也不进 CI 断言。
- **重复性**：QA-17 的第二个断言组要求两次布尔调用的 chart 体积**逐位相同**；
  本页的耗时数字是单机单次观察，不声称跨机器稳定性。

## QA-17 实测

`tests/unit/curved_boolean_qa17_test.cpp`，37° 放置、`d = 1.5`：

| 项 | 值 |
| --- | --- |
| 期望交集体积 | `11π/96 = 0.359974158223830` |
| 相对容差 | `1e-6`（chart 求和） |
| 拓扑 | `validate_topology().ok()`，恰好 1 个 shell |
| 重复性 | 同参数两次运行，chart 体积逐位相同 |
| 断言 | 10 条全部通过，Debug 与 Release 一致 |
| 墙钟（Debug，3 次） | 中位数 `10.35 s`（`10.08–10.50`） |
| 墙钟（Release，3 次） | 中位数 `1.03 s`（`0.94–1.31`） |

## 15° 已知难例实测

| 项 | 值 |
| --- | --- |
| 期望交集体积 | `11π/96 = 0.359974158223830`（`1e-6`） |
| 期望并集体积 | `245π/96 = 8.017606251348951`（`1e-6`） |
| 拓扑 | 交、并各自 `validate_topology().ok()`、1 个 shell、单条共享 rim |
| 断言 | 21 条全部通过 |
| 墙钟（Debug，3 次） | 中位数 `11.62 s`（`11.31–11.65`） |
| 墙钟（Release，3 次） | 中位数 `1.02 s`（`1.01–1.04`） |

## 耗时汇总

| 用例 | 布尔次数 | 配置 | 墙钟中位数 | 范围 | 单次布尔（约） | 断言 |
| --- | --- | --- | --- | --- | --- | --- |
| QA-17（37°） | 2 | Debug | `10.35 s` | `10.08–10.50` | `5.1 s` | 10 |
| QA-17（37°） | 2 | Release | `1.03 s` | `0.94–1.31` | `0.49 s` | 10 |
| 15°（交 + 并） | 2 | Debug | `11.62 s` | `11.31–11.65` | `5.8 s` | 21 |
| 15°（交 + 并） | 2 | Release | `1.02 s` | `1.01–1.04` | `0.49 s` | 21 |
| 旋转语料 | 33 | Debug | `175.04 s`（单次） | — | `5.3 s` | 396 |
| 旋转语料 | 33 | Release | `14.97 s` | `14.77–18.60` | `0.45 s` | 396 |

Debug 比 Release 慢约 `10×`（`/Od` 与 container 检查），但两者比例稳定；15° 与
37° 的 Release 墙钟相差不到 `0.02 s`，与"任意方向走同一条重画/切分路径、没有
按角度分支"一致。

## 全套

| 配置 | 用例 | 断言 | 结果 | 全套墙钟 |
| --- | --- | --- | --- | --- |
| Debug（`build`） | 952 | 660,012 | 全部通过，0 失败，0 跳过 | `761 s` |
| Release（`build-rel`） | 952 | 660,012 | 全部通过，0 失败，0 跳过 | `74.7 s` |

基线为 951 个用例 / 660,002 条断言；本任务新增 1 个用例（QA-17）与 10 条断言，
其余用例与断言数不变。两套构建在 `/W4 /WX /utf-8` 下零警告（GEOM_WERROR=ON）。
Debug 全套中旋转语料占约 `175 s`（约 23%），QA-17 约 `10 s`（约 1.3%）。

## 结论

- QA-17 在非网格角度 37° 通过：单壳、拓扑干净、chart 体积 `11π/96` 在 `1e-6`
  内、两次逐位相同 → 一般位置球/球已被接受，不再是"斜置拒绝"。
- 成本仍由布尔本身主导：Release 约 `0.5 s/次`，Debug 约 `5.3 s/次`，主要是重画
  两个整球 + 约 4098 个 rim 采样点的弧愈合 + 两次切分/装配；chart 体积度量不是
  瓶颈。
- 37° 的 QA-17 与 15° 的既有难例在 Release 下几乎同速，说明"一般位置"没有引入
  额外的慢路径。

## 范围与限制

- 只接受**整球**操作数：球帽、修剪球面、球 + 额外面都拒绝
  （`"...need whole-sphere operands"`）。
- 只接受**严格相交**的两球：包含、外切、内切、分离拒绝
  （`"...need two crossing spheres"`）；M3c-3c-1 没有放开这些。
- 中心线方向不再受限（任意方向都接受）；8 ulp 比较只用于重画后的不动点终止，
  不是放置门。
- 没有 rim snapping：Task 3 实测两片 cap 的 rim 本就相同，投影到交圆的补丁是
  惰性的，未进入生产代码。
- 锥、torus、NURBS、三元布尔、布尔结果修剪面的网格化（`tessellate`）仍未支持。
- lens rim 在 `1e-9` 下约 4098 点，`curved_solid_volume` 仍不收敛；QA-17 继续用
  chart cap 求和度量。斜置图卡上 `face_loops_uv` 有约 `2e-8` 的极点短差，仅测试
  helper 用 `1e-6` 容忍。
- 旋转语料中 5 个 `mt19937` 抽样的具体角度/距离依赖标准库
  `uniform_real_distribution` 的实现；换工具链位置会移动，但仍在 `[0,360)` 与
  `[1.05,1.95]` 内。

## 复现

```powershell
cmake --build build --config Debug --target geom_tests
.\build\tests\Debug\geom_tests.exe --test-case="*QA-17*"
.\build\tests\Debug\geom_tests.exe --test-case="*off-axis sphere lens*"
.\build\tests\Debug\geom_tests.exe --test-case="*rotated placements*"

cmake --build build-rel --config Release --target geom_tests
.\build-rel\tests\Release\geom_tests.exe --test-case="*QA-17*"
.\build-rel\tests\Release\geom_tests.exe --test-case="*off-axis sphere lens*"
.\build-rel\tests\Release\geom_tests.exe --test-case="*rotated placements*"
```

计时（同一命令包一层 `Measure-Command`，无其他负载时连跑 3 次取中位数）：

```powershell
Measure-Command { .\build-rel\tests\Release\geom_tests.exe --test-case="*QA-17*" }
```

全套：

```powershell
.\build\tests\Debug\geom_tests.exe
.\build-rel\tests\Release\geom_tests.exe
```
