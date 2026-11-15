# M3c-3a 球面帽曲面布尔 QA-15 与基准（2026-11-15）

本文件记录 M3c-3a 出口的两部分证据：QA-15 三行球帽矩阵，以及
`boolean_solid_curved` 在三个球操作上的逐操作耗时。数字为单机单次测量，基准
没有时间门禁，只记录实测值，不是跨机器承诺。
**结论：QA-15 3/3 通过（功能 PASS，两次运行体积逐位相同）；基准中位数
布尔调用 21.1 / 22.1 / 31.2 ms，体积积分 232.0 / 232.0 / 214.8 ms，
其中体积积分是主要成本，没有时间门禁。**

## 环境

| 项 | 值 |
| --- | --- |
| 日期 | 2026-11-15（Asia/Shanghai） |
| OS | Windows 10 企业版 10.0.19045，x64 |
| CPU | Intel Core i7-8550U @ 1.80 GHz（4 核 8 线程，15 W 移动处理器） |
| 内存 | 15.9 GB |
| 编译器 | MSVC 19.44.35217.0（VS 2022 17.14，工具集 14.44.35207），`/W4 /WX /utf-8 /fp:precise` |
| 构建 | CMake 4.4.4，VS 2022 generator，`-DGEOM_WERROR=ON -DGEOM_BUILD_BENCH=ON -DGEOM_BUILD_TESTS=OFF` |
| 配置 | Release（`/O2`），单线程 |
| 提交 | 基线 `fc0b076`（M3c-3a Task 5）；QA-15、基准与文档在同一工作树上实测后提交 |

## 口径

- **操作数**：`make_sphere({0,0,0}, 1.0, 1e-9)`（单位球）与三个轴对齐盒体：
  帽盒 `(-2,-2,0.5)..(2,2,2)`（`z >= 0.5` 的球帽，帽高 `h = 0.5`）、
  同一帽盒用于差集（球去掉帽）、平板盒 `(-2,-2,-0.8)..(2,2,0.8)`
  （`|z| <= 0.8` 的 band）。
- **解析体积**：球体积 `4π/3 = 4.1887902047863905`；球帽体积
  `πh²(3R−h)/3`，`h = 0.5` 时为 `0.6544984694978736`，补集为
  `3.534291735288517`；band 为球减去两个 `h = 0.2` 的帽，
  `3.9542179533184135`。
- **QA-15**（`tests/unit/brep_boolean_curved_test.cpp`）：三行每行跑两次，要求
  `boolean_solid_curved` 成功、`validate_topology().ok()`、恰好一个 shell、
  体积相对误差 `<= 1e-6`；第二次运行体积与第一次逐位相同（`==`）。
  cap 与 complement 还要求把每个面投影成 `CurvedBooleanFaceInfo`、用
  `detail::curved_boolean_solid_mesh` 建分类网格并通过
  `detail::curved_boolean_mesh_closure(...).closed()`、`used_once == 0`、
  `nonmanifold == 0`。slab 的 band 面轮廓是两条跨缝合线的圆（8194 个 uv 点），
  `curved_boolean_mesh_face` 的 `triangulate_polygon2` 拒绝该轮廓
  （`"boolean_solid_curved: a face outline could not be triangulated"`），因此
  该行改用 `geom::detail::heal_detail::shell_is_closed_by_uses`：shell 中每条边
  恰被两个 coedge 使用。这是本里程碑已知限制，记录在 ADR-0029 与 §范围与限制。
- **基准**（`bench/bench_m3c3.cpp`）：每个操作先跑一次不计时的预热（校验
  拓扑与闭合 shell），再用 `std::chrono::steady_clock` 计时 100 次；计时区间
  分成一次完整的 `boolean_solid_curved` 调用与一次
  `detail::curved_solid_volume` 调用两部分，分别给出 min/median/p95/max。
  每次迭代都按 1e-6 相对校验体积，并与第一次的体积位型逐位比较；100 次体积
  的 FNV-1a 64 位 checksum 打印在行内与结尾。拓扑只在第一次与最后一次校验。
  退出码：0 全部校验通过（含 checksum 一致），2 构造/驱动/校验错误或
  checksum 不一致。与 `bench_m3c2.cpp` 的 1/2 约定不同，本基准按任务要求把
  校验失败统一为 2。
- **单位**：毫秒；无时间门禁，`PASS` 只代表校验通过。

## QA-15（3 行）

| # | 用例 | 操作 | 解析体积 | 实测体积 | 面数 | cuts | kept | 残差 | 确定性 | 结果 |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| 1 | `sphere ∩ z >= 0.5`（球帽） | intersection | 0.6544984694978736 | 0.654498469497871 | 2 | 2 | 2 | 7.8e-13 | 逐位相同 | PASS |
| 2 | `sphere − z >= 0.5`（帽补集） | difference | 3.534291735288517 | 3.534291735288498 | 2 | 2 | 2 | 7.8e-13 | 逐位相同 | PASS |
| 3 | `sphere ∩ |z| <= 0.8`（band） | intersection | 3.9542179533184135 | 3.954217953318337 | 3 | 4 | 3 | 1.5e-12 | 逐位相同 | PASS |

三行体积残差都在 `1e-12` 量级，`validate_topology` 无发现，实体是单壳闭合
shell。行 1 的 cut 是纬度圆 `v = asin(0.5)`，球面 hat 触及北极，按“恰好一个
精确极点”闭合；行 2 是同一纬度圆下方的补集，触及南极；行 3 有两条纬度圆
`v = ±asin(0.8)`，在 `u = 0` 缝合列拼成一个 band 环，`cuts = 4` 对应两条球面
切链加两个平面 cap 圆。行 1/2 的分类网格闭合；行 3 的 band 轮廓无法三角化，
按 coedge 使用计数校验闭合（见上）。

QA-15 单用例 83 条断言；`geom_tests.exe --test-case="QA-15*"` 1 个用例、
83 条断言通过，同一次运行打印 "939 skipped"，与全量 940 个用例一致。

## 基准

`build-m3c3-bench\Release\geom_bench_m3c3.exe`（每操作 100 次采样，退出码 0）：

```
m3c3-sphere-cap case=cap runs=100 boolean_min_ms=17.794 boolean_median_ms=21.095 boolean_p95_ms=29.894 boolean_max_ms=35.656 volume_min_ms=196.109 volume_median_ms=232.044 volume_p95_ms=322.898 volume_max_ms=380.303 faces=2 expected_volume=0.654498469497874 measured_volume=0.654498469497871 cuts=2 kept=2 max_cut_residual=7.83373e-13 checksum=0xfe892d7046ac8d33 verified=yes
m3c3-sphere-cap case=complement runs=100 boolean_min_ms=16.836 boolean_median_ms=22.108 boolean_p95_ms=36.526 boolean_max_ms=44.890 volume_min_ms=185.272 volume_median_ms=231.961 volume_p95_ms=350.442 volume_max_ms=373.614 faces=2 expected_volume=3.534291735288517 measured_volume=3.534291735288498 cuts=2 kept=2 max_cut_residual=7.83373e-13 checksum=0x7b90ebb9ee71755b verified=yes
m3c3-sphere-cap case=slab runs=100 boolean_min_ms=26.680 boolean_median_ms=31.217 boolean_p95_ms=38.632 boolean_max_ms=49.554 volume_min_ms=196.084 volume_median_ms=214.838 volume_p95_ms=242.443 volume_max_ms=271.046 faces=3 expected_volume=3.954217953318413 measured_volume=3.954217953318337 cuts=4 kept=3 max_cut_residual=1.47691e-12 checksum=0xf1dc2ecedb553d0b verified=yes
m3c3-sphere-cap compiler=msvc os=windows pointer_bits=64 runs=100 checksum=0x723cfe20657534f3 verified=yes
```

- 三个操作的布尔调用中位数 21.1 / 22.1 / 31.2 ms；体积积分中位数
  232.0 / 232.0 / 214.8 ms。体积积分是布尔调用的约 7–11 倍，因为修剪面环在
  uv 上按每个 coedge 32 段采样，积分器要把 uv 多边形三角化并在每个子三角形上
  打 7 个点，子三角形数量由 uv 边长决定（上限 8 级）。
- band 一行的布尔调用最贵（31.2 ms，约 cap 的 1.5 倍）：三个面、四条 cut、
  两个平面 cap 圆都要解析愈合，球面 band 还要拼环；cap 与 complement 的面数、
  cut 数相同，耗时也接近（21.1 vs 22.1 ms）。
- **波动**：本次测量与 M3c-2 在同机（i7-8550U，15 W）完成；体积积分的 p95 与
  max 明显高于中位数（例如 cap 的 322.9 / 380.3 ms），主要来自同机其他进程
  争用。量级稳定，绝对数字只作本机参考。
- 100 次体积逐位相同，三个行内 checksum 与结尾语料 checksum
  `0x723cfe20657534f3` 一致；`verified=yes`、退出码 0。

## 结论

1. **QA-15 PASS**：三行（球帽、帽补集、两个纬度环之间的 band）都是
   `validate_topology` 干净的闭合单壳实体，体积相对误差 `<= 1e-6`，两次运行
   体积逐位相同；cap/complement 的分类网格闭合，band 行按每条边恰被两个
   coedge 使用校验闭合。
2. **基准无时间门禁**：记录的是本机实测值，布尔调用中位数 21.1 / 22.1 /
   31.2 ms，体积积分中位数 232.0 / 232.0 / 214.8 ms，三项全部校验通过
   （`verified=yes`，退出码 0）。
3. 球面闭合纬度环操作数与 `u = 0` 缝合列拼接都可以端到端使用，结果可重复；
   球/球 lens、锥、偏移/倾斜 torus 切平面、一般修剪网格化、NURBS 仍是显式
   `unsupported` 或未实现。
4. 全套测试 Debug 与 Release 各 940 个用例、659,153 条断言全绿；
   `/W4 /WX /utf-8` 零警告。

## 范围与限制

- **只支持闭合纬度环**：开放切链、完全位于域内部的闭环、跨缝合线两次的闭环
  继续返回 `unsupported`（消息见 ADR-0029）；球/球 lens（`d = 1.5`）不在本
  里程碑，需要跨链结合与精确球/球圆。
- **band 分类网格不可用**：slab 的 band 面轮廓是两条跨缝合线的圆，`mesh_face`
  的 `triangulate_polygon2` 拒绝该轮廓；QA-15 该行改用 coedge 使用计数，一般
  修剪网格化（含跨缝轮廓与孔洞）留给 M3c-3c。这不影响布尔结果、体积或
  `validate_topology`。
- **锥、偏移/倾斜 torus、精确 torus 求交**：锥需要顶点投影与闭合/开放切链
  拓扑；torus 的偏移平行/倾斜切平面与 torus/cylinder 配对未做，留给 M3c-3b。
- **NURBS 操作数**：仍是显式 `unsupported`，留给 M3c-3c。
- **修剪结果不能网格化**：布尔输出的修剪面边界是非直边/含孔洞，`tessellate`
  按 ADR-0027 确定性拒绝；体积只能走 `detail::curved_solid_volume`。
- **分类网格弦高固定**：默认 1e-3，驱动无旋钮；基准数字只对应这一采样。
- **单线程单机**：Release `/O2`、单线程、单机单次测量，不与 M3c-2 的圆环
  数字直接比较。

## 复现

```powershell
cmake -S . -B build-m3c3-bench -DGEOM_WERROR=ON -DGEOM_BUILD_BENCH=ON -DGEOM_BUILD_TESTS=OFF
cmake --build build-m3c3-bench --config Release --target geom_bench_m3c3
build-m3c3-bench\Release\geom_bench_m3c3.exe       # 3 操作 x 101 次，约数分钟
.\build\tests\Debug\geom_tests.exe --test-case="QA-15*"
.\build\tests\Debug\geom_tests.exe                 # 全量 940 个用例
```
