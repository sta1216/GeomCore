# M3b-1 QA-9 SSI 基准与退化语料（2026-11-11）

本文件记录 M3b-1 出口 QA-9 的两部分证据：200 例退化语料的诊断结果，以及单位球-平面
横截求交的逐对耗时。数字为单机单次测量，门禁按计划固定值判定，不是跨机器承诺。
**结论：QA-9 语料 200/200 全部被诊断（无效输出 0），正确性通过；基准中位数
15.91 ms/对，远超 100 µs/对的目标（约 159 倍），性能门禁 FAIL。**

## 环境

| 项 | 值 |
| --- | --- |
| 日期 | 2026-11-11（Asia/Shanghai） |
| OS | Windows 10 企业版 10.0.19045，x64 |
| CPU | Intel Core i7-8550U @ 1.80 GHz（4 核 8 线程，15 W 移动处理器） |
| 内存 | 15.9 GB |
| 编译器 | MSVC 19.44.35217.0（VS 2022 17.14，工具集 14.44.35207），`/W4 /WX /utf-8 /fp:precise` |
| 构建 | CMake 4.4.4，VS 2022 generator，`-DGEOM_WERROR=ON -DGEOM_BUILD_BENCH=ON` |
| 配置 | Release（`/O2`），单线程 |

## 口径

- **QA-9 语料**：`tests/unit/ssi_degenerate_test.cpp` 的
  `ssi: QA-9 degenerate corpus is either stitched or diagnosed`，四族各 50 例，共
  200 例：球-平面南极点相切、不同写法/尺度的重合平面、同轴等半径圆柱、圆环-平面
  相切。每例调用 `ssi(..., SsiOptions{1e-8})`；若返回曲线，要求每段
  `residual ≤ 1e-8`；若无曲线，要求至少一条 `SsiIssue` 且每个诊断码非 `none`。
  用例末尾断言 `cases == 200` 且 `diagnosed == 200`。
- **基准**：单位球补丁 `r = 1`，`u ∈ [0, 2π]`、`v ∈ [-1.5, 1.5]`；平面
  `z = 0.1`，参数域 `[-2, 2]²`；`SsiOptions{tolerance = 1e-9}`（`max_depth = 12`、
  `max_segments = 100000`）。重复 1000 次，每次只对 `ssi()` 计时；第一次调用不计入
  样本（预热）。任一次计时调用无曲线即记为结果不可用（FAIL，退出码 1）；构造或
  `ssi()` 返回错误才退出码 2。
- 计时用 `std::chrono::steady_clock`，微秒保留三位小数。通过条件取**中位数**
  ≤ 100 µs/对；min/median/max 全部打印，便于判断是否只是尾部长。
- 退出码：0 门禁通过，1 门禁未过或结果不可用，2 构造/SSI 错误。

## QA-9 语料

四族均为按构造退化，因此每例都必须落到诊断分支：

| 族 | 夹具 | 例数 | 实测诊断 |
| --- | --- | --- | --- |
| 1 | 球 `r = 0.5 + 0.01·i`，平面 `z = -r` 在南极点相切；球补丁 `v ∈ [-π/2, π/2]` | 50 | `tangent_contact` ×50 |
| 2 | 平面 `z = d` 与 `z = d`（法向写作 `{0,0,-3}`、偏移 `-3d`），`d = 0.1·i` | 50 | `coincident` ×50 |
| 3 | 同轴等半径圆柱，轴心沿 z 偏移 `0.01·i`，半径 1 | 50 | `coincident` ×50 |
| 4 | 圆环 `R = 2.0 + 0.01·i`、`r = 0.5`，与平面 `z = 0.5` 相切 | 50 | `tangent_contact` ×50 |

实测：`cases = 200`、`diagnosed = 200`（100 例 `tangent_contact` + 100 例
`coincident`，`none` 为 0），无曲线分支被触发；该用例 1002 条断言全部通过（200 例
各 5 条，加末尾 2 条计数断言）。说明相切/重合没有被误当横截线输出，也没有静默空集。

## 基准

`build-rel\Release\geom_bench_ssi.exe`（默认 1000 次），退出码 1：

```
QA9 calls=1000 min_us=9393.900 median_us=15910.550 max_us=24480.700 target_us=100 FAIL
QA9 curves=1000 segments=31000 candidate_pairs=200000 rejected_segments=65000 max_residual=6.441e-10
```

- 每次调用都返回 1 条闭合曲线，1000 次共 31,000 段，`max_residual = 6.441e-10`
  （门禁 1e-9），正确性通过。
- 每次约 200 个 AABB 候选三角形对、96 个提取线段（31 段输出 + 65 段饱和差额）。
- **中位数 15.91 ms/对，目标 100 µs/对，约 159 倍；min 9.39 ms、max 24.48 ms。
  PASS/FAIL 判定为 FAIL，不改门禁。**
- 同机此前一次运行的中位数为 12.97 ms（min 8.23、max 21.71），本次高约 23%。15 W
  移动处理器在持续满载下降频会带来这一量级的波动；两次都远未接近 100 µs。
- 计时区间只有 `ssi()`；补丁构造、结果校验与统计都在区间之外。每次调用都重新
  细分两张补丁、重新播种并追踪，没有跨调用缓存，这是当前实现的结构性成本。

## 结论

1. **QA-9 语料 PASS（正确性）**：200/200 被诊断（100 `tangent_contact` +
   100 `coincident`，含重合平面的反法向写法），无曲线输出、无静默空集；残差断言未
   被触发。
2. **QA-9 基准 FAIL（性能）**：1000 次单位球-平面横截求交的中位数 15.91 ms/对，
   目标 100 µs/对，约 159 倍；min/median/max 分别 9.39 / 15.91 / 24.48 ms。
3. 正确性侧全部通过：每次调用 1 条曲线，共 31,000 段，最大残差 6.441e-10（≤1e-9）。
4. 全套测试 Debug 与 Release 各 673 个用例、549,443 条断言全绿；`/W4 /WX /utf-8`
   零警告；QA-9 用例是回归通过（功能在 Task 1–Task 9 已完成，未观察到 RED）。

## 范围与限制

- **只覆盖未修剪参数补丁**：入口是 `SsiPatch`，不是 B-rep 面；面级调用、多段 pcurve、
  曲面面分割与曲面布尔属于 M3b-2。
- **数值交线，不是精确有理**：残差由 `tolerance` 约束；ADR-0015 已排除全局精确有理
  SSI。相切/重合是诊断，不是精确拓扑。
- **细分是启发式**：弦高容差 `max(1e-3·diagonal, 64·tolerance)`，每格 9 个采样点；
  不同深度的叶子不边共形，存在悬挂节点；最坏成本按 `4^max_depth` 增长。
- **相切检测下限**：几何带 `max(1e-12·(1+max|p|), 0.05·tolerance)`；更浅的切入按
  契约报 `tangent_contact`。
- **近相切分支可能合并**：种子聚类半径约 0.007×较小补丁对角线，间距更小的分支会被
  并成一条曲线；由 Task 6 谓词诊断，是 M3b-1 的已知限制。
- **混合尺度重复曲线残差**推迟到 M3b-2；M3b-1 只保证单次调用的残差与确定性。
- **性能缺口**：每次调用无缓存地重新细分/播种/追踪，中位数约 16 ms/对；100 µs/对
  目标在当前实现下未达成。
- **M3b-3 healing 未交付**；QA-8（M3a 平面布尔性能）仍为 FAIL，不由本任务改变。

## 复现

```powershell
cmake -S . -B build-rel -DGEOM_WERROR=ON -DGEOM_BUILD_BENCH=ON
cmake --build build-rel --config Release --target geom_tests geom_bench_ssi
ctest --test-dir build-rel -C Release --output-on-failure
build-rel\Release\geom_bench_ssi.exe            # 1000 次，约 16 s
build-rel\Release\geom_bench_ssi.exe 100        # 冒烟
```

QA-9 语料也可单独运行：

```powershell
build-rel\tests\Release\geom_tests.exe --test-case="ssi: QA-9*"
```
