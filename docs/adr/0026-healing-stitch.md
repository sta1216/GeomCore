# ADR 0026: M3b-3 的初始缝合与容差链

- Status: Accepted
- Date: 2026-11-13
- Related: architecture §4.6、§8.12、§12.6、§16.3；M3b-3 plan Task 3–Task 5；ADR-0021、ADR-0022、ADR-0025

## 背景

QA-12 要求：1e-6 量级的缝隙装配必须要么闭合成通过校验的实体，要么带缝隙位置与
量值拒绝，不能静默返回一个看起来闭合的 `Brep`。M3b-3 的输入是单实体、单壳或
多壳、平面面且边为直线（`kNoGeometry` 或 `LineCurve3`）的开口壳；典型场景是
盒体上盖与侧壁之间留下一条均匀缝隙。

现有 `SolidEditor` 已经提供 `merge_vertices`、`merge_edges`、`set_shell_closed`、
`set_coedge_pcurve`、`set_edge_geometry` 与 `raise_tolerances(cap)`，但还没有
“找配对 → 执行 → 校验 → 拒绝或发布”的事务层。直接把最近顶点焊掉会在
T 型/竞争配对、方向翻转、pcurve 失效和容差爆炸上产生不可诊断的坏结果。

## 选项

1. **全局 epsilon 或最近点全焊。** 实现最短，但会把非互近邻、同面边和只是
   共享一个顶点的边也焊在一起；无法解释拒绝了哪条缝，也无法保证容差链单调，
   与 §4.6 的反模式直接冲突。
2. **互近邻配对 + 显式容差链 + staged 事务。** 先测量全部候选，再按
   `max_gap` 接受；只缝双方都首选彼此的边；每次修改都走 `SolidEditor`，
   `commit()` 通过 `validate_topology`/`validate_tolerances` 后才发布。实现量
   可控，拒绝路径可定位、可量化。
3. **同时做小边折叠、小面删除与曲面重拟合。** 能修复更多模型，但会引入
   局部重拟合、误差上界与 `lossy` 标记，超出 M3b-3 的验证范围。
4. **多实体与曲面缝隙一起支持。** 需要跨实体容差带与曲面 pcurve 重建，
   当前没有足够的 QA 语料，容易把“未验证”误报成“已支持”。

## 决策

采用选项 2，并把契约固定为：

- **范围。** `heal_brep` 只处理目标实体的边界边；目标实体的每个面必须是
  `Plane3`，每条被使用的边曲线必须是 `kNoGeometry` 或 `LineCurve3`，否则在
  修改前返回 `Errc::unsupported`。其他实体可以存在于同一 `Brep` 中，但必须
  满足结构前置条件且不会被修改；目标实体的边界边不得被另一个实体共用。
- **配对规则。** `max_gap` 扮演 §8.12 中“两侧容差带相交”的判据：每条边界边
  先测量所有候选，按端点两种对应方式中字典序更优的 `(最大端点距离, 距离和)`
  取 `gap`；只有 `gap ≤ max_gap` 的候选才进入接受集。接受集按
  `(gap, sum, keep, drop)` 排序，双方都首选彼此的边才缝合。若一条边的两个
  接受候选 `gap` 相差不超过 `1e-12`，报 `ambiguous_gap` 并排除它的全部候选；
  只共享顶点或同属一个面的边不作为候选。
- **方向规则。** 不使用 coedge 的原始 `reversed` 标志直接比较。先把 coedge 的
  遍历起点按配对焊点映射到存活顶点，再与 `keep.v0` 比较，得到焊后方向；
  两侧方向必须相反，否则报 `inconsistent_orientation`，不缝。
- **容差规则。** 焊点容差取 `max(keep.tol, drop.tol, gap)`；边容差支配其端点，
  面容差支配其边与 coedge，实体容差支配其面。所有提升只增不减，并由
  `raise_tolerances(max_tolerance)` 统一执行；任何会超过 `max_tolerance` 的
  提升都拒绝。预检查还覆盖“输入链本身已高于 cap、提升函数因此不会再报错”的
  情况。`max_gap` 与 `max_tolerance` 都必须是有限正数。
- **事务性。** 全部修改在 `SolidEditor` 的 staged 副本上完成；pcurve 修复、
  壳闭合、容差提升之后调用 `commit()`。`commit()` 的拓扑/容差校验失败、壳
  仍未闭合、Euler 特征变化或任何 primitive 失败都返回**原输入** `Brep`
  （`changed=false`，计数为 0），绝不发布部分缝合结果。成功时报告
  `max_closed_gap`、`max_tolerance`、唯一焊点数与缝合边数。
- **pcurve 修复。** 焊接会移动端点，因此**每个脏面（拥有与焊后存活顶点相邻的边的面）**上的所有 coedge 都要重投影：按升序 coedge 索引，在面自己的 uv 中写入逐段折线 pcurve（圆弧保留内点、周期缝展开、极点/锥顶借用邻点方位），三维载体保留；投影失败即拒绝整个 heal。（本节的两点 pcurve + `kNoGeometry` 契约已被 2026-11-24 补遗取代。）
- **拒绝诊断。** 无接受候选且最近候选超过 `max_gap` 时，全局只报一条
  `gap_too_large`，带最近距离、最近端点对中点与两条边的索引；竞争候选只产生
  `ambiguous_gap`，不会同时产生 `gap_too_large`。只有共享顶点、没有互近邻时
  报 `invalid_result "heal_brep: no mutual pair"`，而不是静默成功。

## 后果

- **QA-12 正例（实测）。** 1e-6 缝隙的裂盒：8 条边界边 → 4 个唯一焊点、
  4 条缝合边，`max_closed_gap=1e-6`、`max_tolerance=1e-6`，闭合壳通过
  `validate_topology`；`detail::curved_solid_volume` 体积为 1，误差
  `1.11e-16`（门禁 1e-9）。同一结果再次 heal 为 `changed=false`、
  `issues.empty()`，且逐字段指纹不变。
- **QA-12 反例（实测）。** 1e-3 缝隙在 `max_gap=1e-6` 下被拒绝：
  `gap=0.001`，位置 `(0, 0, 0.5005)`（最近端点对中点），边索引 16/20，
  `changed=false`，输入指纹不变。
- **可诊断性。** 竞争配对、方向不一致、容差超限、曲线输入、未知实体分别有
  稳定的 `Errc` 或 `HealIssueKind` 与固定消息；拒绝路径不会返回半成品。
- **已知偏差（§4.6 invariant 5）。** 报告只记录 `max_closed_gap` 与
  `max_tolerance`，不记录每次容差放大的来源与量值；逐次 provenance 延后到
  M4。当前测试通过“输入指纹不变”和“所有顶点 tol ≤ cap”覆盖可观察部分。
- **范围限制。** 小边/小面移除、曲面重拟合、多实体缝合、曲面缝隙桥接仍延后
  到 M4；M3b-3 只交付初始焊接 + 缝合 + 容差提升。平面/直线限制不是静默降级：
  曲面面或非直线边在入口返回 `unsupported`。
- **测试证据。** Debug 与 Release 各 859 个用例、623,976 条断言全绿；
  `/W4 /WX /utf-8` 零警告。QA-12 的表格与复现命令见
  `docs/benchmarks/m3b3-healing.md`。

## 重新评估的触发条件

M4 加入小边/小面折叠、曲面重拟合、多实体缝合或曲面缝隙桥接时；或出现需要
逐次容差 provenance 的审计/交换场景时；或 QA-12 语料扩展到非均匀缝隙、
T 型接缝与多壳实体时。

---

## 补遗（2026-11-24）：M4a 曲面缝隙愈合

M4a 保持上文平面契约（配对、方向、容差与 staged 事务）不变，把同一
`SolidEditor` 事务扩展到解析曲面与圆弧边；以下是新增的契约。

- **支持类别（白名单）。** 目标实体的每个面必须是 `Plane3`、`Cylinder3`、
  `Cone3`、`Sphere3` 或 `Torus3`；每条被使用的边曲线必须是 `kNoGeometry`、
  `LineCurve3` 或 `CircleCurve3`。其他曲面/曲线（含 NURBS）在修改前分别以
  `unsupported "heal_brep: surface type is not supported yet"` /
  `"heal_brep: curve type is not supported yet"` 拒绝。`Torus3` 只通过类型门
  （允许同一实体中的 torus 邻面），torus 缝隙的 v 周期与双半径配对尚未验收，
  愈合本身延后到 M4b；不得据此声称支持 torus 缝隙。
- **载体判据。** `detail::heal_surfaces_coincide` 要求两个候选面是同一解析类型，
  并在 `check_tolerance = max(options.tolerance, 1e-12 * surface_scale)` 内比较：
  平面之间恒相容（平面缝合不依赖放置）；圆柱比较半径、轴方向与轴线偏移；球比较
  半径与球心；锥比较半径、半角、轴方向与顶点；torus 比较主/次半径、轴与中心。
  类型不同、或已支持解析量超出该带，报
  `unsupported "heal_brep: gap faces do not share the same surface"`，context 为
  (keep 边, drop 边)。每条边界边还要在 `t0`、参数中点与 `t1` 三点落在自己所属面
  的曲面上（残差 ≤ 同一带），否则报
  `"heal_brep: gap edge is not on its face surface"`；被保留边的三维曲线必须穿过
  焊后两个端点，否则报
  `"heal_brep: kept edge curve does not match its welded vertices"`。这些门都在
  `SolidEditor` 第一次修改之前运行。
- **方向 tie-break（仅非平面）。** 半径/纬线关于中心对称的载体上，两个候选的端点
  距离可以逐位相同，仅靠 `(gap, sum, keep, drop)` 无法分辨真缝与对称诱饵。对差距
  在 `1e-12 * max(1, best_gap)` 内的候选，只有当该边界边首个 coedge 所在面不是
  `Plane3` 时，才用 `SolidEditor::merge_edges` 将来采用的同一方向谓词（焊后遍历
  起点相同）挑出唯一方向相容的候选；方向相容候选不是恰好一个时仍报
  `ambiguous_gap`。缺失几何按平面处理，因此四条 `make_ambiguous_rims` 的
  `ambiguous_gap` 期望与消息完全不变。
- **逐面 pcurve 契约。** 焊接后对目标实体的每个面（不止缝隙面）重建 pcurve：
  按面/环/coedge 顺序取三维路径（直线两端点；圆按不超过 0.1 rad 的密集采样，
  闭合边取 ±2π 全周），用 `previous_u` 提示投影到该面图卡并展开周期 u/v；反向
  coedge 反向遍历采样；极点/锥顶样本保持精确 v、u 借用相邻非奇异样本，避免扇形
  积分成零；周期缝的整周回描在 |Δu| > π 时平移整数个 2π 以闭合图卡；去掉共线
  内点，但圆弧段至少保留一个内点使弧向可判。每条 coedge 写入至少两点，否则报
  `Errc::degenerate "heal_brep: pcurve has fewer than two points"`（context 为
  coedge 索引）。三维边载体不清除，`CircleCurve3`/`LineCurve3` 仍可供后续操作使用。
- **pcurve 读取规则。** `face_uv` 优先使用存储的 pcurve（反向 coedge 反向读取），
  否则采样三维边再投影；每段圆弧的扫角必须落在 ±π 内（|sweep| ≥ π-1e-6 报
  `numeric_failure "face_loops_uv: circular arc direction is ambiguous"`）；点偏离
  圆的程度超过 `max(tolerance, 1e-12*max(1,radius))` 报
  `"face_loops_uv: point is off its circle"`；环在三维中判闭合、去掉闭合重复点后
  不足三点报 `"face_loops_uv: loop has fewer than three points"`。
- **闭合 rim 诊断。** 若边界边的 `v0 == v1`（`make_cylinder` 式的单条整圆边），
  它无法配对：焊接会把两端塌到同一点，方向检查也没有遍历起点。`heal_brep` 在
  边界表与白名单已知之后、`heal_plan`/任何修改之前返回
  `Errc::unsupported "heal_brep: closed gap edge cannot be paired; split the rim
  first"`，context 为 (edge 索引, 0)。`detail::heal_plan` 对直接调用者仍保留
  `invalid_topology "heal_plan: boundary edge has zero length"`（两个不同顶点
  恰好重合的零长边仍走该门）。`heal_brep` 内部不拆分 rim；拆分是 `split_solid`
  的职责（`collapse_circular_chains` 是 M4b 的复用件）。
- **无限直线限制。** `heal_curve_residual` 对 `LineCurve3` 取点到无限载体直线的
  距离，不是到线段 `[t0,t1]` 的距离；`heal_edge_lies_on_surface` 在 `t0`、中点、
  `t1` 三点取样，因此它验证载体直线与这三个参数点，边端与顶点的一致性由
  kept-curve 门禁另行保证。不要把它当作线段范围校验。
- **延期项。** torus 缝隙愈合（v 周期、双半径、缝配对）延后到 M4b；顶点/极点
  缩成一点的 apex/pole gap 不在本里程碑范围（锥带终止于顶点、边界仍是有限半径
  圆的用例已接受）；NURBS 曲面/曲线、曲面重拟合与多实体缝合仍不支持。
- **平面 QA-12 不变。** 本补遗没有改变任何平面数值：裂盒体积 1（误差
  `1.11e-16`）、4 焊点/4 缝合边、`max_closed_gap=max_tolerance=1e-6`、1e-3 缝隙
  仍在 `(0,0,0.5005)` 拒绝（边 16/20）、幂等指纹不变；四条 `ambiguous_gap` 与
  既有消息均保持。新增的唯一行为是闭合 rim 的显式 `unsupported`。
- **测试证据（截至 Task 6）。** 新增
  `heal: a closed gap rim is reported instead of a zero-length edge`：在旧树上
  RED（`invalid_topology` 与 `"heal_plan: boundary edge has zero length"`），实现后
  GREEN，且输入指纹不变。Release 目标过滤
  `heal:*,curved volume:*,tessellate:*,validate:*` 实测 68 个用例、4,274 条断言、
  0 失败、1,000 跳过；全量 Debug/Release 只在里程碑评审运行。三个验收值
  （Release）：圆柱 `π*(2+1e-6)` rel 2.5e-14、球 `4π/3` rel 2.3e-15、锥 `2π/3`
  rel 4.2e-15；`chi == 2`；细节见 `docs/benchmarks/m4a-healing.md` 与
  architecture §16.4。
- **缺少曲面。** `surface == kNoGeometry` 的面在修改前返回 `invalid_topology "heal_brep: face has no surface"`（Task 1 起的诊断，未单独测试）。

---

## 补遗（2026-11-24）：M4b-1 闭合 rim 与 torus 缝隙

M4b-1 取代 M4a 补遗中的两条结论：闭合整圆 rim 不再以 `unsupported` 拒绝，而是在
私有副本上自动拆分后愈合；torus 缝隙（v 周期、双半径、闭合 rim 半体）进入支持
范围。平面契约与 M4a 的载体/方向/pcurve 契约不变；以下是新增与修正。

- **闭合 rim 自动拆分（私有副本）。** `heal_brep` 在边界表、白名单与载体门禁
  通过之后、`heal_plan` 与任何修改之前调用
  `detail::heal_split_closed_rims(brep, solid)`：目标实体没有 `v0 == v1` 的边界边
  时返回空 optional（零开销）；存在时用 `SolidEditor` 复制输入，对每条此类边在
  参数中点 `edge.t0 + 0.5*(edge.t1-edge.t0)` 调 `split_edge` 后提交。拆分只发生在
  工作副本，调用者的 `Brep` 不被修改；`split_edge` 失败时返回该错误、staging 副本
  保持未动，`heal_brep` 同样不改动输入。M4a 的
  `unsupported "closed gap edge cannot be paired; split the rim
  first"` 诊断随之删除，`detail::heal_plan` 对直接调用者仍保留零长边的
  `invalid_topology`。
- **中点拆分后的计数。** 每条闭合圆边拆成两段参数各半的圆弧，工作副本的边界边
  数暂时翻倍：圆柱闭合 rim 输入 `V/E/F = 4/6/4`，其中 2 条闭合缝边是边界边；
  拆分后 `6/8/4`、4 条边界边。torus 闭合 rim 半体输入 `4/6/2`、4 条闭合边界边，
  拆分后 `8/10/2`、8 条边界边。拆分顶点容差为 0，随后由普通互近邻配对焊回，
  因此愈合结果回到输入规模（圆柱 `4/6/4`、torus 半体 `4/6/2`）。报告中的
  `boundary_edges` 记录的是拆分后的工作计数（圆柱 4、torus 半体 8），不是输入
  的 2/4。
- **端点定向焊接。** `heal_plan` 不再按顶点索引排序后归一化焊点，而是以保留边
  （`merge_edges` 的幸存载体）为基准：端点同向时焊 `{keep.v0→drop.v0,
  keep.v1→drop.v1}`，反向时焊 `{keep.v0→drop.v1, keep.v1→drop.v0}`。闭合 rim
  拆开后，较小顶点索引可能属于被丢弃的 rim，旧的 `normalize_weld` 会把保留载体
  的端点换到对侧，导致 kept-curve 门禁在大缝隙下误报；端点定向消除了该缺陷，
  `normalize_weld` 已删除。
- **同缝方向归一化。** 首对候选固定一条缝的焊接方向；同一缝的后续配对（拆开的
  rim 的另一半）若其焊点已被反向占用，则交换该对的 keep/drop 并同时交换每个
  焊点的两端，使所有配对共享同一幸存端点。距离并列仍只在首个 coedge 邻面非
  `Plane3` 时用焊后遍历方向挑选唯一相容候选，否则报 `ambiguous_gap`；平面行为
  不变。
- **条件三分 pcurve。** 闭合边按全周采样；当共线简化后只剩两端点、且中间采样点
  与任一端点的 u 或 v 跨过图卡精确中点（`|Δ| > π`）时，
  `face_uv_unwrap` 会把略超 π 的前向扫掠折回上一列，把 rim pcurve 压成尖刺并让
  环的 uv 面积减半。此时 `heal_rebuild_face_pcurves` 改在 1/3 与 2/3 采样处存点
  （而非中点），使弧向与面积都稳定；不触发该条件的紧凑 pcurve 采样完全不变
  （既有 `lateral_points == 10` 用例保持）。
- **圆柱闭合 rim 验收。** `make_cylinder_closed_rims_gap(1e-6)`：拆分后工作副本
  `6/8/4`，愈合 `4/6/4`、单一闭合壳、无边界边、`chi == 2`、
  `max_closed_gap = max_tolerance = 9.9999999991773336e-07`、2 焊点/2 缝合边，
  体积 `6.2831884487723997` 对 `π*(2+1e-6)` 相对 `2.5e-14`，重复运行逐位一致，
  二次 heal 为 `changed == false`。
- **torus 开缝带验收。** 三个开缝带（u 缝、外支 v 缝、跨 v 缝的内支）输入均
  `8/10/4`，愈合 `6/8/4`、`chi == 2`、2 焊点/2 缝合边、愈合后无边界边；体积分别为
  `2.3561944901923608`（对 `3π/4`）、`5.4334195273316128`、
  `0.8635058157199863`，相对参考均不超过 `1.0e-14`；实测缝宽
  `9.999999901566576e-07`、`9.9999998999999992e-07`、
  `9.9999998992555949e-07`（夹具按 `gap*(1-1e-8)` 留余量，稳定落在默认
  `max_gap` 内）；内支图卡 v 跨度小于 2.0，反向 lateral 与 `Torus3(v_origin=π)`
  均被正确处理。
- **torus 闭合 rim 半体验收与 chi 修正。** 两个半体在 `z = ±h`（`h = (gap/2) *
  (1-1e-8)`）切开，`alpha = asin(h/r)`；输入 `4/6/2`、4 条闭合边界边，拆分后
  `8/10/2`、8 条边界边，愈合 `4/6/2`、单一闭合壳、愈合后无边界边、4 焊点/4 缝合边，
  `max_closed_gap = 9.9999999030728907e-07`，体积 `4.934802200544679` 对
  `π²/2` 绝对误差不超过 `1e-12`，四段图卡窗口 `[-alpha, π-alpha]` 与
  `[π-alpha, 2π-alpha]` 各在 `1e-9` 内。**`chi == 0`**：计划与简报曾把 torus
  半体写成 `chi == 2`，那是错的；`euler_characteristic`
  （`validate.hpp:905-959`）按 `V - E + F - (L - F)` 计算，闭合 torus 的
  `V=4, E=6, F=2` 给出 0。所有 M4b-1 文档、计划与 benchmark 以 0 为准。
- **拒绝后的返回语义（重要）。** 拆分发生在工作副本，因此愈合失败时返回的
  `HealResult::brep` 是**调用者的原始输入**（拓扑、几何存储、容差逐位不变），
  而 `report.boundary_edges` 是**拆分后的工作计数**，不是输入的闭合 rim 条数：
  圆柱为 4、torus 半体为 8。两类拒绝都满足 `changed == false`、
  `welded_vertices == sewn_edges == 0`：`max_gap` 小于实测缝宽时报
  `gap_too_large`；焊接容差超过 `max_tolerance` 时报
  `invalid_result "heal_brep: tolerance cap exceeded"`（容差取保留/丢弃顶点容差
  与配对 gap 的最大值）。
- **测试锁。** `heal: a rejected closed-rim heal returns the caller's input` 用
  `max_gap=1e-12` 触发 `gap_too_large`、用 `max_tolerance=1e-9` 触发
  `invalid_result "heal_brep: tolerance cap exceeded"`。两个分支都断言返回 Brep
  与输入的**完整指纹**逐项相同：实体/几何存储计数、顶点位置与容差、边端点/曲线
  索引/参数域/容差、面曲面/反向/容差、coedge 的边/面/环/反向/pcurve/next/prev/
  容差、loop 的 face/start/next/outer，以及 `report.boundary_edges == 4`。这样
  即使实现把工作副本的一部分泄漏进返回值，或把报告计数错写成输入值，测试都会
  失败。
- **负例。** `gap=1e-3` 仍报 `gap_too_large`（拆分后边界 8，entity_a/b = 0/1）；
  平面盖报 carrier 消息（context 1/4）；off-surface rim 报 (2/6)；NURBS 面/边与
  歧义缝复用 M4a 用例；所有拒绝的输入指纹不变。
- **延期项。** 顶点/极点缩成一点的 apex/pole gap 延后到 M4b-2；fillet/chamfer/
  shell 在更晚的里程碑；NURBS 曲面/曲线、曲面重拟合与多实体缝合仍不支持；
  partial-u（只拆闭合 rim 的一段、保留部分整圆边）不做，本里程碑总是整条边拆分；
  torus 的 u 周期多圈缠绕仍不在验收范围。
- **测试证据。** 新增/扩展用例见 `tests/unit/brep_heal_curved_test.cpp` 与
  `tests/unit/m4a_curved_heal_fixture.hpp`（简报中的
  `brep_heal_curved_fixture.hpp` 不存在）。Release 快速通道 `-tc="*heal*"` 为
  91 用例/1,866 断言、`-tc="*torus*"` 为 59/10,754、合并 145/12,440，0 失败；
  全量 Debug/Release 已运行：各 1074 用例 / 680,503 断言，0 失败。实测表见
  `docs/benchmarks/m4b-healing.md`，架构状态见 §16.4。
