# M3b-2 曲面布尔 QA 语料与基准（2026-11-12）

本文件记录 M3b-2 出口的两部分证据：十二例 QA 语料（六个正例 + 六个诊断反例）与
`boolean_solid_curved` 单位圆柱-盒体的逐操作耗时。数字为单机单次测量，门禁按计划
固定值判定，不是跨机器承诺。
**结论：QA 语料 12/12 通过（功能 PASS）；基准并/交/差中位数 156.9 / 156.3 / 152.2 ms，
目标 2000 ms/操作，PASS。**

## 环境

| 项 | 值 |
| --- | --- |
| 日期 | 2026-11-12（Asia/Shanghai） |
| OS | Windows 10 企业版 10.0.19045，x64 |
| CPU | Intel Core i7-8550U @ 1.80 GHz（4 核 8 线程，15 W 移动处理器） |
| 内存 | 15.9 GB |
| 编译器 | MSVC 19.44.35217.0（VS 2022 17.14，工具集 14.44.35207），`/W4 /WX /utf-8 /fp:precise` |
| 构建 | CMake 4.4.4，VS 2022 generator，`-DGEOM_WERROR=ON -DGEOM_BUILD_BENCH=ON` |
| 配置 | Release（`/O2`），单线程 |
| 提交 | 基线 `fba2de2`；Task 9 的基准、语料与文档在同一工作树上实测后提交 |

## 口径

- **操作数**：`make_cylinder({0,0,0}, {0,0,1}, 1.0, 2.0, 1e-9)`（r = 1、h = 2，
  下帽在 z = 0）与 `make_box({0.5,-0.5,0.5}, {2.0,0.5,1.5}, 1e-9)`（偏移帽盒，
  体积 1.5）。圆柱体积 2π，盒体高 1 且底面在圆盘内，交集是一个高 1 的圆盘弓形。
- **解析体积**：令 `I = π/6 + √3/4 − 1/2`（盒体足印与单位圆盘的公共面积），
  则 intersection = I ≈ 0.456611477490518、union = 1.5 + 2π − I ≈ 7.326573829689、
  difference = 2π − I ≈ 5.826573829689。
- **基准**：每个操作先跑一次预热（不计入样本），再用 `std::chrono::steady_clock`
  计时 5 次；计时区间是一次完整的 `boolean_solid_curved` 调用（两次分类网格、
  SSI、切面、分类、装配）。每次结果都用 `geom::validate_topology` 与
  `detail::curved_solid_volume` 对照解析值，相对容差 1e-9。
- **采样**：分类网格由默认 `TessellationOptions`（弦高 1e-3）构建，单位整圆
  71 段，圆柱分类网格 142 顶点 / 280 三角形。驱动没有弦高旋钮，计划里原先的
  "5000 侧面采样"不可表达，这里记录**实际**段数，不虚构更细的采样。
- 通过条件取**中位数** ≤ 2000 ms/操作，min/median/max 全部打印。退出码：
  0 全部校验通过（时间目标即使未达也打印 FAIL 但退出 0），1 校验或体积失败，
  2 构造/驱动错误。

## QA 语料（12 例）

`tests/unit/brep_boolean_curved_test.cpp` 的
`brep: boolean_solid_curved QA corpus` 用一张表跑十二例；每例跑两次，正例要求
体积完全相同、反例要求 (code, message) 完全相同。

| # | 用例 | 操作 | 解析体积 | 实测体积 | 面数 | 结果 |
| --- | --- | --- | --- | --- | --- | --- |
| 1 | 偏移帽盒 | union | 1.5 + 2π − I = 7.326573829689 | 7.326573829689 | 8 | PASS |
| 2 | 偏移帽盒 | intersection | I = 0.456611477491 | 0.456611477490 | 7 | PASS |
| 3 | 偏移帽盒 | difference | 2π − I = 5.826573829689 | 5.826573829689 | 8 | PASS |
| 4 | 缝上盒（z ∈ [1, 10]） | union | 144 + π = 147.141592653590 | 147.141592653590 | 8 | PASS |
| 5 | 平板盒（x ≥ 0.5） | intersection | 2(π/3 − √3/4) = 1.228369698608 | 1.228369698609 | 5 | PASS |
| 6 | 贯穿帽盒（z ∈ [1.5, 3]） | intersection | 0.587720000938 | 0.587720000939 | 6 | PASS |
| 7 | 球 `r = 2` + 盒 | union | — | `unsupported`：sphere operands are deferred to M3c | — | PASS（诊断） |
| 8 | 圆环 R = 2, r = 0.5 + 盒 | union | — | `unsupported`：torus operands are deferred to M3c | — | PASS（诊断） |
| 9 | 锥 + 盒 | union | — | `unsupported`：cone operands are not supported（apex-touching loops are deferred） | — | PASS（诊断） |
| 10 | NURBS 底面的盒 + 盒 | union | — | `unsupported`：NURBS operands are not supported | — | PASS（诊断） |
| 11 | 盒面与圆柱相切（x = ±1、y = ±1） | union | — | `unsupported`：tangent face contact is not supported yet | — | PASS（诊断） |
| 12 | 盒底落在上帽平面 z = 2 | union | — | `unsupported`：coplanar face contact is not supported yet | — | PASS（诊断） |

正例的六个体积都与独立推导的解析值相差 ≤ 1e-6 相对（实测残差在 1e-12 量级），
且 `validate_topology` 无发现；反例六个都是 `Errc::unsupported`，消息与上表逐字
一致且非空。语料用例本身 152 条断言全通过；`--test-case=*boolean_solid_curved*`
共 17 个用例、416 条断言通过。

## 基准

`build-m3b2-rel\Release\geom_bench_m3b2.exe`（5 次采样，退出码 0）：

```
m3b2-curved-boolean op=union runs=5 min_ms=154.132 median_ms=156.857 max_ms=180.794 target_ms=2000 faces=8 PASS
m3b2-curved-boolean op=intersection runs=5 min_ms=151.417 median_ms=156.328 max_ms=162.705 target_ms=2000 faces=7 PASS
m3b2-curved-boolean op=difference runs=5 min_ms=145.203 median_ms=152.209 max_ms=180.953 target_ms=2000 faces=8 PASS
m3b2-curved-boolean circle_segments=71 deflection=1e-03 tessellation_vertices=142 tessellation_triangles=280
m3b2-curved-boolean compiler=msvc os=windows pointer_bits=64 target_met=yes verified=yes
```

- 三个操作的中位数都在 152–157 ms，远低于 2000 ms 门禁；同机第一次运行为
  154.4 / 161.2 / 163.4 ms（min 144.5、max 172.4），批次间波动约 10%，量级不变。
- 每次调用的 5 个样本全部通过 `validate_topology` 与 1e-9 相对体积校验；
  面数与解析结构一致（8/7/8）。
- 成本结构：每次调用都重新构建两个分类网格（圆柱 280 三角形 + 盒体 12 三角形）、
  枚举 6 个面片对做 SSI、切面并装配；没有跨调用缓存。相对 M3a 平面布尔
  （5002 面棱柱 636 s）小三个数量级，是因为这里只有 3 + 6 = 9 个面。

## 结论

1. **QA 语料 PASS**：六正例（含缝、边界、双帽三类切割）全部是
   `validate_topology` 干净的实体，体积相对误差 ≤ 1e-6；六反例全部是带稳定
   消息的 `unsupported`，无静默成功、无空消息；12/12 确定性通过。
2. **基准 PASS**：并/交/差中位数 156.9 / 156.3 / 152.2 ms，目标 2000 ms/操作。
3. 正确性、性能两项同时达成，按计划可以勾选 M3b-2 出口。
4. 全套测试 Debug 与 Release 各 797 个用例、623,289 条断言全绿；`ctest`
   Debug 189.1 s、Release 25.8 s 通过；`/W4 /WX /utf-8` 零警告。

## 范围与限制

- **解析作用域**：只支持 `Plane3`/`Cylinder3` 操作数；球、圆环、锥、NURBS 是
  显式 `unsupported`（M3c）。
- **共面与相切不支持**：两者都诊断为 `unsupported`，不猜测空/满分类。
- **翻转曲面片段不支持**：差集/对称差中需要翻转朝外的曲面片段在切面前拒绝。
- **单连通结果**：不相交并集与腔体结果诊断为 `unsupported`。
- **修剪曲面面不能网格化/算质量属性**：`tessellate` 需要未修剪的闭合帽圆，
  `mass_properties` 需要未修剪整矩形；体积只能走 `detail::curved_solid_volume`
  （参数域积分 + 散度定理），修剪网格化是 M3c 工作。
- **分类网格弦高固定**：1e-3（整圆 71 段），驱动无旋钮；更细采样要等 M3c 的
  分类网格预算接口。
- **三链共点**：`distinct_cuts` 组内饱和为 2，三链共点按两链处理（未触发）。

## 复现

```powershell
cmake -S . -B build-m3b2-rel -DGEOM_WERROR=ON -DGEOM_BUILD_BENCH=ON
cmake --build build-m3b2-rel --config Release --target geom_tests geom_bench_m3b2
ctest --test-dir build-m3b2-rel -C Release --output-on-failure
build-m3b2-rel\Release\geom_bench_m3b2.exe       # 3 操作 x 6 次，约 3 s
build-m3b2-rel\tests\Release\geom_tests.exe --test-case="*boolean_solid_curved*"
build-m3b2-rel\tests\Release\geom_tests.exe --test-case="brep: boolean_solid_curved QA corpus"
```
