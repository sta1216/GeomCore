# ADR 0031: M3c-3c-1 一般位置球/球透镜

- Status: Accepted
- Date: 2026-11-16
- Related: architecture §8.10、§16.3；M3c-3c-1 plan Task 1–Task 5；ADR-0029、
  ADR-0030（其"中心线必须轴对齐"的前置条件由本 ADR 取代）

## 背景

ADR-0030 把球/球 lens 关闭为"两个整球 + 相交 + 中心线在 8 ulp 内与坐标轴对齐"，
并把一般位置留给 M3c-3c。M3c-3c-1 把可接受的球心连线扩到任意方向后，通用放置
暴露出两个数值边界缺陷，二者都与"斜置图卡"有关：

1. **缝端残差。** SSI 采样点投影回图卡后可能落在参数域内几个 ulp。缝对闭合逻辑
   容忍这种漂移，但 arrangement 不容忍：切链悬空，产生一个零宽度面。Task 1 实测
   一个 `v = asin(0.5)` 的 25 点纬度环，端点 `u = 2.9101132668452026e-18`（缝合线
   内侧 2.9e-18）在修复前得到 52 = 2×25+2 点、`signed_area2 == 0` 的单个片段。
2. **seam-only 判定依赖精确零。** `split_face_curved` 旧代码用
   `fabs(signed_area2(loop)) == 0.0` 识别"整球缝环"。斜置图卡上 64 点缝环的
   shoelace 面积是舍入噪声 `6.1629758220391547e-33`，不是零；于是一个整球面被送进
   常规裁剪路径，其面积过滤把全部片段丢掉，装配只剩另一个球的两片 cap，报
   `assemble_curved_solid: edge is used once`（对外映射为
   `boolean_solid_curved: result touches itself along an edge`）。

原计划还假设第二个缺陷是"两片 cap 的 rim 各自经不同图卡抬起，顶点/弧表示不同"，
准备把 rim 点投影到球/球交圆上。Task 3 的探针推翻了该假设：两片 cap 的 rim 在本
路径上已经逐点相同，投影补丁完全惰性。该补丁已删除，**本里程碑没有 rim
snapping**；本 ADR 只记录实际进入生产代码的规则。

## 选项

1. **一般球面/球面裁剪域 + 多环网格化。** 对球面片、锥面片统一做 NURBS 轮廓、
   多环域三角化与精确求交。终态，但属于 M3c-3c 后续里程碑，无法在单一出口验证
   lens；ADR-0030 已拒绝把它与球/球混在一起。
2. **保留 M3c-3b 的重画路径，修两个数值边界并放开中心线方向。** 球/球交圆在
   共享图卡下仍是精确纬度，切分/分类/装配复用 M3c-3a 的闭合纬度环路径；只需让
   "整球缝环"的识别不依赖轴对齐。选用。
3. **继续要求中心线轴对齐。** 即 ADR-0030 的现状。任意一对 `make_sphere` 默认
   输入几乎都不满足，等于只支持特制输入；M3c-3c-1 的目的正是去掉它。拒绝。
4. **rim 投影到交圆。** 针对"两片 rim 不同"的诊断，实测两片 rim 已经相同，
   补丁不改变任何结果（惰性）。不采用，也不进入生产代码。

## 决策

采用选项 2。

- **接受输入。** 两个操作数必须各是单个整球（单壳、单面、面为 `Sphere3`、
  `synthetic_domain`、恰 1 条缝边与 2 个顶点），否则报
  `"boolean_solid_curved: sphere/sphere booleans need whole-sphere operands"`
  （`boolean_solid_curved.hpp:1790-1794`）。两球必须严格相交：
  `d > tolerance && d < ra + rb - tolerance && d > |ra - rb| + tolerance`，
  否则报 `"boolean_solid_curved: sphere/sphere booleans need two crossing
  spheres"`（`boolean_solid_curved.hpp:1801-1805`）。包含、外切、内切、分离仍
  没有交圆，统一在此拒绝。求交、并、差三个操作走同一条切分/装配路径。
- **任意中心线。** `curved_boolean_rechart_whole_spheres` 直接取
  `axis = delta / distance`（`boolean_solid_curved.hpp:1807`），不再计算主分量、
  不再有 `kAxisAlignment` 门，也不再返回 "placement is not axis-aligned"。
  第二球用 `-axis`；`candidate = |axis[0]| >= 0.9 ? (0,1,0) : (1,0,0)`，
  `x_axis = normalize(cross(axis, candidate))`（`boolean_solid_curved.hpp:1815-1820`）。
  8 ulp 比较（`kAxisAgreement`，`boolean_solid_curved.hpp:1809-1813`）只用于
  重画后的**不动点终止**：重建后的 `axis`/`x_axis` 精确等于同一段确定性代码算出
  的值，所以递归深度恒为 2（一次重画 + 一次不动点），不加递归计数器。它不是
  放置门，任何方向的斜置都接受。重画只改图卡朝向，圆心、半径、容差语义与报告
  计数器不变。
- **缝端点吸附。** `detail::trim_face_split_seam_only` 在整链 `u` 平移之后、
  域边界检查之前，把距 `u_lo`/`u_hi` 不超过
  `seam_eps = 1e-12 * max({1.0, |u_lo|, |u_hi|})` 的坐标精确吸附到该侧；当
  `projection.periodic_v()` 时对 `v_lo`/`v_hi` 用同一 `seam_eps`
  （`trim_face.hpp:988-1003`）。随后从吸附后的点重算 `min_u`/`max_u`，再做
  `"split_face_curved: cut leaves the parameter domain"` 检查
  （`trim_face.hpp:1004-1015`），因为旧值早于吸附。该规则只改变距域边界
  1e-12 相对量级的窄带，把"缝对已经当作在缝上"的点真正接到域侧，避免悬空切链
  与零宽度面；它不改变任何拓扑谓词，也不是全局 epsilon。
- **seam-only 相对面积判定。** 面边界环的 |uv 面积| 满足
  `|signed_area2(loops.front())| <= 1e-12 * extent^2` 时，才把该面当作"缝环"
  并打开整张参数矩形（`trim_face.hpp:1511-1526`）；`extent` 是面边界的 uv
  对角线 `hypot(max_u - min_u, max_v - min_v)`（`trim_face.hpp:1493`）。为什么
  用相对面积而不是类型门：被真正修剪过的球面片同样是 `Sphere3`，仅凭
  `trim_face_seam_only_domain(surface)` 的类型判断会把它误送进"整张矩形"路径；
  而真实修剪边界的面积是 `extent^2` 的显著份额，比 `1e-12 * extent^2` 高若干
  个数量级，因此面积量级能把"缝环的舍入噪声"和"真实区域"分开——斜置图卡缝环
  实测 `6.16e-33`，轴对齐缝环恰为 `0`，被修剪的面远高于界。该常数是判定用的
  经验界，不是模型容差：它必须远大于缝环的舍入面积（约 1e-32 量级），又必须
  远小于任何真实 uv 区域面积。
- **没有 rim snapping。** 本里程碑不把 rim 点投影到交圆，也不新增任何 rim 吸附
  代码。Task 3 的探针显示两片 cap 的 rim 在共享图卡重画下已经相同，投影补丁
  不改变结果；rim 一致性由"两球共轴、共享 `x_axis`、交圆是两张图卡的同一条
  纬度"这一构造保证，并由既有装配的顶点/弧焊接收口。这是对原计划的偏离，
  必须在文档中如实记录，不能声称存在 rim snapping。
- **测试 helper 的极点容差。** 测试本地 `sphere_face_cap_volume`（现位于
  `tests/unit/brep_boolean_curved_fixture.hpp`）用 `face_loops_uv` 的环做
  shoelace 定符号、从 `±π/2` 极点求 cap 高。斜置图卡上极点通过
  `asin(dot(offset, axis) / r)` 恢复，`|axis|^2 = 1 ± 1e-16` 的舍入被 `asin`
  放大成约 `2e-8` 的 `v` 偏差；因此 helper 接受距极点 `1e-6` 以内的环，诱导的
  cap 高误差是 `O(1e-12)`，远低于 `1e-6` 的体积断言。这是**测试代码**的容差；
  生产侧 `face_loops_uv` 仍有约 `2e-8` 的极点短差，是本里程碑记录的限制。

## 后果

- **QA-17（新增）。** `tests/unit/curved_boolean_qa17_test.cpp`，单位球对、
  `d = 1.5`、中心线方向 37°（不在 15°/30°/10° 网格上）：交集体积
  `11π/96 = 0.359974158223830...`（`1e-6` 相对，chart 求和度量）、恰好 1 个
  shell、恰好 2 个面（两片 cap）、单条被两个 coedge 反向使用的 rim、
  `validate_topology().ok()`、同参数连跑两次 chart 体积逐位相同；17 条断言
  （含 `face_ids().size() == 2` 与 `check_one_shared_rim` 的 6 条 rim 断言），
  Debug 与 Release 都通过。QA-17 使用测试共享头
  `tests/unit/brep_boolean_curved_fixture.hpp` 的 `sphere_face_volume` 与
  `check_one_shared_rim`（由 `brep_boolean_curved_test.cpp` 迁出，均未改名），
  不引入新的体积积分。
- **放置语料（实测）。** 33 次旋转语料由四部分构成：30° 全圈交集 12 次、
  0°/45°/90° 的并/差 6 次、0–90° 以 10° 为步长的细网格交集 10 次、固定种子
  随机抽样 5 次，共 396 条断言。15° 已知难例是**独立用例**（既有的
  `curved boolean: an off-axis sphere lens is built`，交、并各一次，2 次布尔、
  21 条断言，现复用 `check_one_shared_rim`），不计入这 33 次。交集相对 `1e-9`、
  并/差相对 `1e-6`，全部通过；每条结果都是单壳、`validate_topology` 干净、
  单条被两个 coedge 反向使用的 rim。随机抽样的具体角度/距离依赖标准库的
  `uniform_real_distribution` 实现，换工具链位置会移动，但仍落在 `[0,360)` 与
  `[1.05,1.95]` 内。
- **全套计数（实测）。** Debug 与 Release 各 **952 个用例、660,019 条断言**，
  0 失败、0 跳过，`/W4 /WX /utf-8` 零警告；基线 951/660,002，本任务新增 1 个
  用例、17 条断言（QA-17 含 6 条共享 rim 断言）。`sphere_face_volume` 与
  `check_one_shared_rim` 迁出只改测试组织，不改变原有用例的断言数。
- **基准与耗时（实测）。** 本任务不新增 bench 可执行文件；QA-17 与 15° 用例
  的聚焦墙钟时间、33 次语料与全套耗时见 `docs/benchmarks/m3c3c.md`（同一台
  i7-8550U / Windows 10 / MSVC 19.44，Debug 与 Release 分列，informational、
  无阈值）。
- **限制（M3c-3c-1 出口）。** 仍只支持整球操作数；仍只支持严格相交；包含、
  相切、分离、非整球/部分球一律 `unsupported`。锥、torus、NURBS、三元布尔、
  一般修剪网格化未做。布尔结果的修剪面仍不能 `tessellate`。lens rim 在
  `1e-9` 下约 4098 个 uv 点，`curved_solid_volume` 的 ear clipping 仍不收敛，
  因此 QA-17 继续用 chart cap 求和度量；`face_loops_uv` 在斜置图卡上有约
  `2e-8` 的极点短差，仅测试 helper 用 `1e-6` 容忍。`seam_eps` 与
  `1e-12 * extent^2` 都是判定用的经验界，前者只覆盖 u 域宽度（periodic-v
  图卡沿用同一值），以后新增 periodic-v seam-only 切分时需要重新论证。
  本 ADR 只关闭 M3c-3c-1，不得据此声称 M3c-3c 或 M3c-3 完成。

## 重新评估的触发条件

M3c-3c 支持包含/相切/分离、锥/球或球/torus 或三元布尔、或布尔结果修剪面网格化
时；`curved_solid_volume` 能处理数千点 rim（届时可移除 QA-17 的 chart 度量，
直接用库体积）；`make_sphere` 的图卡约定（缝合线、极点、`x_axis`）改变，或
`kAxisAgreement` 需要放宽到超过 8 ulp 时；`trim_face_split_seam_only` 的
seam-only 面积界或 `seam_eps` 在更大/更小尺度的图卡上被实测推翻时。
