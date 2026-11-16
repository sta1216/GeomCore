# ADR 0030: M3c-3b 球/球透镜（lens）布尔

- Status: Accepted
- Date: 2026-11-16
- Related: architecture §8.12、§16.3；M3c-3b plan Task 1–Task 5；ADR-0028、ADR-0029

## 背景

M3c-3a 关闭了球帽（闭合且过缝合线一点的纬度环），并把球/球 lens（两半径 1、球心
距 1.5 的两球）显式留给 M3c-3b（ADR-0029 的限制清单）。lens 的困难不在求交：两球
边界交于一个精确圆，圆心在两球心连线上、平面垂直于该连线，`curved_boolean_pair_circles`
在 M3c-2 已能给出这个圆。真正的问题有三件：

1. **缝与交圆不共面。** `make_sphere` 的缝合线是一条经线（子午线），而球/球交圆
   位于垂直于中心线的平面内。只有两球的 `axis` 恰好沿中心线（且 `x_axis` 相同）
   时，交圆才是两张图卡共有的纬度；否则交圆与缝斜交，`split_face_curved` 的
   seam-only 分支按"切链必须终止在缝上"裁剪，会报
   `"split_face_curved: open cuts on a seam-only face"`，lens 无法切分。
2. **两个球面帽不等于球帽。** lens 的结果壳由两个球面 cap 组成，共享一条 rim 边；
   每个 cap 是"球面片"而不是球，M3c-3a 的合成矩形域与纬度环逻辑都不直接适用。
3. **体积积分退化。** lens rim 在容差 `1e-9` 下约 4098 个 uv 点，`detail::curved_solid_volume`
   的 ear clipping 在数千点轮廓上不收敛（Task 4 观察：分钟级不返回），因此不能用
   它验证 QA-16 的 lens 体积。

## 选项

1. **通用球面/球面裁剪域。** 对球面片、lens、锥面片统一做 NURBS 轮廓 + 多环域
   三角化 + 精确求交。优点是终态；缺点是一般修剪网格化是 M3c-3c 的独立里程碑，
   M3c-3b 无法在单一出口验证 lens。
2. **把整球重画到共享图卡（re-chart），复用闭合纬度环路径。** 两球都是整球时，
   用 `make_sphere(center, axis = ±n, 共享 x_axis, radius)` 重建两个操作数（n 为
   中心线单位向量）。这样交圆变成两图卡的纬度、缝是同一条经线，lens 正好落进
   M3c-3a 已支持的"闭合且过缝一点的纬度环"切分/分类/装配路径，零新增切分代码。
   代价是需要显式的前置条件：操作数必须是整球，两球必须相交，中心线必须轴对齐。
   优点是把"球/球"从"新曲面类型 + 新裁剪域"降级为"图卡朝向变换"，可在
   M3c-3b 内端到端验证，且不触碰 M3c-3a 的 cap 分支。
3. **只在 axis 已沿中心线时支持，否则拒绝。** 改动最小，但 `make_sphere` 默认把
   极轴沿 z、`x_axis` 沿 x，任意一对球几乎都不能直接满足，等于只支持特制输入。
   评审已拒绝。
4. **继续拒绝球/球，等 M3c-3c 与一般修剪网格化一起做。** 缺点：lens 是最容易
   验证的球/球结果，推迟它会让 M3c-3c 同时承担一般位置球/球、锥、修剪网格化三块
   风险，没有中间验证点。

## 决策

采用选项 2。

- **整球判定。** 新增 `detail::CurvedBooleanSphereOperand` 与
  `detail::curved_boolean_sphere_operand`（`boolean_solid_curved.hpp:1703-1744`）：
  单壳、单面、面是 `Sphere3`、`curved_boolean_face_info` 报 `synthetic_domain`、
  文档恰有 1 条边与 2 个顶点（缝合线的一次往返）。任一条不满足就不是整球，配对
  返回 `"boolean_solid_curved: sphere/sphere booleans need whole-sphere operands"`
  （`boolean_solid_curved.hpp:1791-1796`）——部分球、修剪球面、球加额外面都不接受。
- **相交判定。** 球心距 `d = |cb - ca|` 必须满足
  `d > tolerance && d < ra + rb - tolerance && d > |ra - rb| + tolerance`
  （`boolean_solid_curved.hpp:1802-1807`），否则报
  `"boolean_solid_curved: sphere/sphere booleans need two crossing spheres"`。
  包含（`d <= |ra-rb|`）、相切（`d == ra+rb` 或 `d == |ra-rb|`）、分离
  （`d >= ra+rb`）都没有交圆，统一在此拒绝，不进入切分。
- **轴对齐。** 取中心线的主分量方向（`dominant`），要求
  `|unit[dominant]| >= 1 - 8*eps`（`kAxisAlignment`，
  `boolean_solid_curved.hpp:1817-1824`）；否则报
  `"boolean_solid_curved: sphere/sphere placement is not axis-aligned; general
  placements are deferred to M3c-3c"`。8 ulp 只吸收浮点舍入，任何有意的斜置都不
  接受——重画后的缝必须落在坐标轴上的极点经线，否则缝愈合跟不上。
- **共享 frame。** `axis` 取中心线的主分量符号单位向量，第二球用 `-axis`；
  `candidate = |axis[0]| >= 0.9 ? (0,1,0) : (1,0,0)`，
  `x_axis = normalize(cross(axis, candidate))`，两球共用，使两张图卡的缝是同一条
  经线。此时交圆（纬度）恰好在一个经线点接触缝，走 M3c-3a 的支持路径。
- **终止性（结构性，不计数）。** 重画后用
  `|dot(first.axis, axis) - 1| <= 8*eps`、`|dot(second.axis, -axis) - 1| <= 8*eps`、
  `|dot(first.x_axis, second.x_axis) - 1| <= 8*eps`（`kAxisAgreement`，
  `boolean_solid_curved.hpp:1827-1832`）判定"已是不动点"，返回 `nullopt` 让驱动
  走原路径。重建后的轴精确等于 `±n`（同一段确定性代码算出），所以递归深度恒为
  2：一次重画，一次不动点。不加递归计数器，以免掩盖破坏该论证的 bug。
- **驱动接入点。** `boolean_solid_curved` 在算完两侧 `face_info` 之后、torus 门
  与 flipped-fragment 预检之前调用 `detail::curved_boolean_rechart_whole_spheres`
  （`boolean_solid_curved.hpp:1923-1938`）；成功则用重建的两个单球文档递归调用
  自身，并把 `failure_report` 原样传下去。放在 torus 门之前是因为球对不含 torus；
  放在 flipped 预检之前是为了让"非整球/不相交"的诊断不被"切链无法修剪"掩盖。
- **图卡透明。** 重画只改 `axis` / `x_axis`（图卡朝向），圆心与半径原样保留，
  因而体积、容差语义、报告计数器与操作数几何都不变；结果 Brep 的每个面仍是
  `Sphere3` + 纬度环 + 极点行，`validate_topology` 与三次布尔（交/并/差）都复用
  现有代码。
- **体积度量与交叉校验。** lens rim 约 4098 点，`curved_solid_volume` 在其上退化，
  QA-16 改用测试本地 `sphere_face_volume`：对每个 `Sphere3` 面用 `face_loops_uv`
  的环做 shoelace 定符号、从 `±π/2` 极点求 cap 高 `h`，累加 `sign*πh²(3r−h)/3`。
  该 helper 先在 M3c-3a cap 语料（`sphere ∩ box z >= 0.5`，`h = 0.5`）上与库体积
  交叉校验一次：库值 `0.6544984694978736`（相对 `1e-6` 内），helper 与库值相差
  `<= 1e-12`；校验通过才用它度量 lens。这一项写进 ADR 是为了让"lens 体积由
  chart 求和得出"不是空口断言。

## 后果

- **QA-16（实测）。** `tests/unit/brep_boolean_curved_test.cpp` 的
  `QA-16: sphere lens booleans are closed and deterministic` 三行（两半径 1、
  球心距 1.5、容差 `1e-9`）：交 `11π/96 = 0.359974158223830`（相对 `1e-9`）、
  并 `245π/96 = 8.017606251348951`（`1e-6`）、差 `39π/32 = 3.828816046562561`
  （`1e-6`）；实测 chart 体积 `0.359974158224599` / `8.017606251349720` /
  `3.828816046562561`，残差 `7.7e-13` / `7.7e-13` / `4.4e-16`。每行：2 个面、
  1 个 shell、`validate_topology().ok()`、`shell_is_closed_by_uses` 通过、恰有
  1 条被两个 coedge 反向使用的 rim 边；同参数连跑两次，chart 体积逐位相同，
  拓扑计数、`cuts`、`kept` 也相同。
- **通用行（实测）。** 半径 1.2 / 0.8、球心距 1.4、`ca = {0.3, -0.2, 0.1}`、
  `cb = ca + {1.4, 0, 0}`（非零原点、中心线沿 +x）：cap 高
  `ha = 3/14 = 0.214285714285714`、`hb = 27/70 = 0.385714285714286`，cap 体积
  `0.162804109089092` 与 `0.313820376355531`，参考 lens 体积由 cap 公式
  （`πh²(3r−h)/3`，不读实现）给出 `0.476624485444623`；交集行在 `1e-9` 相对内
  一致，同一行也验证 2 面/1 壳/闭合/单 rim/两次逐位相同。该行固定了"不等半径 +
  非零原点 + 中心线沿坐标轴"这条一般性，而不只是单位球在原点。
- **方向性。** `A ∩ B`、`A ∪ B`、`A − B` 三个操作各自命中一条结果：交与并对称，
  差 `A − B` 取 A 的补集减 B 的小帽，helper 的 shoelace 符号（`reversed` 已烘进
  几何）保证带符号求和为正；若方向错，体积会变成负值或补集值，QA-16 会在
  `1e-9`/`1e-6` 内失败。
- **负向路径（已有用例）。** 非轴对齐放置（沿 `unit(1,2,2)`、同一 1.2/0.8/1.4
  配对）报 `"boolean_solid_curved: sphere/sphere placement is not axis-aligned;
  general placements are deferred to M3c-3c"`；非整球、包含/相切/分离各报上面
  两条消息。这些都在 QA-16 之外由既有用例覆盖，本任务只补 lens 语料。
- **基准（`geom_bench_m3c3b`，Release，20 次/操作，informational 无阈值）。**
  布尔调用中位数 `421.0` / `416.4` / `421.3` ms，chart 体积中位数
  `5.34` / `3.56` / `3.70` ms；20 次体积逐位相同，checksum 与退出码 0 见
  `docs/benchmarks/m3c3b.md`。布尔调用是主要成本（重画 + 224 段 rim 愈合 +
  两次切分），体积度量已不再是瓶颈——与 M3c-3a 的"体积积分是主要成本"相反，
  因为 lens 换成了 chart 求和。
- **限制清单（M3c-3b 出口）。** 只支持整球操作数（球帽、修剪球面、球+额外面
  拒绝）；只支持两球相交（包含、外切、内切、分离拒绝）；中心线必须在 `8*eps`
  内轴对齐（一般位置 deferred to M3c-3c）；lens rim 在 `1e-9` 下约 4098 点，
  `curved_solid_volume` 会卡住（ear clipping），QA-16 因此用 chart cap 体积 +
  3a 语料上的一次交叉校验；球/锥、球/torus、三元布尔、一般修剪网格化与 NURBS
  仍未做。本 ADR 只关闭球/球 lens，不得据此声称 M3c-3 完成。
- **复用。** 整球判定与共享 frame 的构造对 M3c-3c 的一般球/球可直接复用；一般
  位置需要额外的"缝绕交圆平面重画"（让缝与交圆共面），这是 M3c-3c 的新工作，
  不是本 ADR 的扩展；若 `make_sphere` 的图卡约定（缝、极点、`x_axis`）改变，
  `kAxisAlignment` / `kAxisAgreement` 与共享 frame 的推导需重新论证。

## 重新评估的触发条件

M3c-3c 支持一般球/球位置（非轴对齐中心线）、包含/相切/分离、或球/锥/球三元时；
`curved_solid_volume` 能处理数千点 rim（届时可移除 QA-16 的 chart 度量折衷，直接
用库体积）；`make_sphere` 的图卡约定（缝合线、极点、`x_axis`）改变时；
`kAxisAlignment` / `kAxisAgreement` 需要放宽到超过 8 ulp 时。
