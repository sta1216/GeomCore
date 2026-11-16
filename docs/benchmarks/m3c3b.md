# M3c-3b 球/球透镜布尔：QA-16 与基准实测

本文件记录 M3c-3b（球/球透镜）的验收数据：QA-16 的体积/拓扑实测、`geom_bench_m3c3b`
的 Release 运行时间，以及该里程碑的适用范围。取舍与设计理由见
[ADR-0030](../adr/0030-sphere-lens-booleans.md)。

## 环境

| 项 | 值 |
| --- | --- |
| CPU / 内存 | Intel Core i7-8550U / 15.9 GB |
| 操作系统 | Windows 10 19045 (x64) |
| 编译器 | MSVC 19.44.35217.0（`/W4 /WX /utf-8`，GEOM_WERROR=ON） |
| CMake | 4.4.4，Visual Studio 17 2022 generator |
| 构建配置 | Release（`/O2`），Debug 仅用于全用例验收 |
| 容差 | `1e-9`（`CurvedBooleanOptions::tolerance`） |

语料：两半径 1 的整球，球心 `(0,0,0)` 与 `(1.5,0,0)`，即 `d = 1.5`、交圆半径
`r = sqrt(1 - 0.75²) = 0.6614`。

## 口径

- **QA-16 的体积度量**：lens rim 在 `1e-9` 下约 4098 个 uv 点，
  `detail::curved_solid_volume` 的 ear clipping 在数千点轮廓上不收敛（Task 4
  观察：分钟级不返回），因此 QA-16 不用它度量 lens。改用测试本地
  `sphere_face_volume`：对每个 `Sphere3` 面取 `face_loops_uv` 的环，用 shoelace
  定符号、从 `±π/2` 极点求 cap 高 `h`，累加 `sign*πh²(3r−h)/3`。
- **交叉校验（一次性）**：同一 helper 先在 M3c-3a 的 cap 语料
  （`sphere ∩ box z >= 0.5`，`h = 0.5`）上与库体积比对：库值
  `0.6544984694978736`（相对 `1e-6` 内），helper 与库值相差 `<= 1e-12`。
  校验通过后才用它度量 lens；该项也写进 ADR-0030，metric 不是空口断言。
- **每行结构判据**：2 个面、1 个 shell、`validate_topology().ok()`、
  `geom::detail::heal_detail::shell_is_closed_by_uses` 通过、恰有 1 条被两个
  coedge 反向使用的 rim 边；同参数连跑两次，chart 体积逐位相同，`cuts` / `kept`
  与拓扑计数相同。
- **基准口径**：每操作 1 次预热 + 20 次计时，两个区间——整个
  `boolean_solid_curved` 调用、结果的有符号 chart cap 体积。每轮都把体积与闭式解
  比对（交 `1e-9` 相对、并/差 `1e-6` 相对）并与首轮逐位比较；首末轮跑
  `validate_topology`、单 solid/单 shell/2 面/闭合壳/每条边恰两用。输出附带
  64 位 FNV-1a checksum 与退出码（0 全部校验通过，2 失败）。基准是 informational，
  没有时间阈值，CI 的 bench job 也是 `continue-on-error`。

## QA-16 实测

`tests/unit/brep_boolean_curved_test.cpp` 的
`QA-16: sphere lens booleans are closed and deterministic`，单位球对
（`d = 1.5`）：

| 操作 | 期望体积 | 实测 chart 体积 | 残差 | 相对容差 |
| --- | --- | --- | --- | --- |
| `A ∩ B`（lens） | `11π/96 = 0.359974158223830` | `0.359974158224599` | `7.69e-13` | `1e-9` |
| `A ∪ B` | `245π/96 = 8.017606251348951` | `8.017606251349720` | `7.69e-13` | `1e-6` |
| `A − B` | `39π/32 = 3.828816046562560` | `3.828816046562561` | `4.44e-16` | `1e-6` |

通用行（不等半径 + 非零原点 + 中心线沿 `+x`）：`ra = 1.2`、`rb = 0.8`、
`ca = (0.3, -0.2, 0.1)`、`cb = ca + (1.4, 0, 0)`；cap 高
`ha = 3/14 = 0.214285714285714`、`hb = 27/70 = 0.385714285714286`，参考 lens 体积
由 cap 公式（不读实现）给出 `0.476624485444623`，交集行在 `1e-9` 相对内一致。
四行都满足上述结构判据与两次逐位相同。单用例 166 条断言通过。

全套：Debug 与 Release 各 **949 个用例、659,595 条断言**全绿，
`/W4 /WX /utf-8` 零警告（GEOM_WERROR=ON）。

## 基准实测（Release，20 次/操作）

`build-m3c3-bench/Release/geom_bench_m3c3b.exe` 的 stdout：

```text
m3c3b-sphere-lens case=intersection runs=20 boolean_min_ms=403.897 boolean_median_ms=421.005 boolean_p95_ms=446.873 boolean_max_ms=457.559 volume_min_ms=4.354 volume_median_ms=5.339 volume_p95_ms=7.266 volume_max_ms=8.040 faces=2 expected_volume=0.359974158223830 measured_volume=0.359974158224599 residual=7.688e-13 cuts=2 kept=2 max_cut_residual=5.5933e-13 checksum=0x38adc65c64aac4a3 verified=yes
m3c3b-sphere-lens case=union runs=20 boolean_min_ms=390.567 boolean_median_ms=416.365 boolean_p95_ms=442.518 boolean_max_ms=448.676 volume_min_ms=2.808 volume_median_ms=3.561 volume_p95_ms=5.628 volume_max_ms=5.756 faces=2 expected_volume=8.017606251348951 measured_volume=8.017606251349720 residual=7.692e-13 cuts=2 kept=2 max_cut_residual=5.5933e-13 checksum=0xf86cbfbe4ffcbca3 verified=yes
m3c3b-sphere-lens case=difference runs=20 boolean_min_ms=401.798 boolean_median_ms=421.284 boolean_p95_ms=445.960 boolean_max_ms=460.860 volume_min_ms=2.832 volume_median_ms=3.699 volume_p95_ms=5.524 volume_max_ms=6.730 faces=2 expected_volume=3.828816046562560 measured_volume=3.828816046562561 residual=4.441e-16 cuts=2 kept=2 max_cut_residual=5.5933e-13 checksum=0xcb6a8a8897ade05b verified=yes
m3c3b-sphere-lens compiler=msvc os=windows pointer_bits=64 runs=20 checksum=0xa9d320a11993421b verified=yes
```

| 操作 | 布尔中位数 | 布尔 p95 | chart 体积中位数 | `cuts` | `kept` | `max_cut_residual` |
| --- | --- | --- | --- | --- | --- | --- |
| 交 | `421.0 ms` | `446.9 ms` | `5.34 ms` | 2 | 2 | `5.59e-13` |
| 并 | `416.4 ms` | `442.5 ms` | `3.56 ms` | 2 | 2 | `5.59e-13` |
| 差 | `421.3 ms` | `446.0 ms` | `3.70 ms` | 2 | 2 | `5.59e-13` |

退出码 0，20 次体积逐位相同（checksum 与上表逐行一致）。

## 结论

- lens 体积与闭式解一致（残差 `<= 7.7e-13`，远小于 `1e-9`/`1e-6` 门限），
  拓扑每行都是 2 面、1 壳、单 rim、闭合；两次运行逐位相同 → 结果确定。
- 成本主导在布尔本身：每次约 `0.42 s`，其中重画两个整球 + 约 4098 个 rim 采样点的精确
  愈合 + 两次切分/装配；chart 体积只占 `3.6–5.3 ms`（约 1%），因为 lens 已不再
  走 `curved_solid_volume`。与 M3c-3a 的"体积积分是主要成本"（布尔 21–31 ms、
  体积 232 ms）相反，这是度量口径变化的结果，不是性能回归对比。
- 三个操作的中位数差异在 `1.2%` 以内，交/并/差共用同一条切分路径。

## 范围与限制

- 只接受**整球**操作数：球帽、修剪球面、球 + 额外面都拒绝
  （`"...need whole-sphere operands"`）。
- 只接受**相交**的两球：包含、外切、内切、分离拒绝
  （`"...need two crossing spheres"`）。
- 中心线必须在 `8*eps` 内与某坐标轴对齐；斜置拒绝
  （`"...placement is not axis-aligned; general placements are deferred to M3c-3c"`）。
  一般位置球/球、包含/相切/分离、球/锥、球/torus、三元布尔、一般修剪网格化与
  NURBS 都留给 M3c-3c。
- lens rim 在 `1e-9` 下约 4098 点，`curved_solid_volume` 会卡住（ear clipping）；
  QA-16 因此度量 chart cap 体积，并在 3a cap 语料上交叉校验一次。

## 复现

```powershell
cmake -S . -B build-m3c3-bench -G "Visual Studio 17 2022" -A x64 `
  -DGEOM_BUILD_BENCH=ON -DGEOM_WERROR=ON
cmake --build build-m3c3-bench --config Release --target geom_bench_m3c3b
.\build-m3c3-bench\Release\geom_bench_m3c3b.exe
# 退出码 0 = 全部体积/拓扑/逐位校验通过；2 = 失败
.\build-m3c3-bench\Release\geom_bench_m3c3b.exe; $LASTEXITCODE
```
