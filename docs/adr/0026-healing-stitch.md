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
- **pcurve 修复。** 焊接会移动端点，因此不仅缝合边，**所有与焊后存活顶点
  相邻的边**都要重投影：按升序 coedge 索引，在每个面自己的 uv 中写入
  `{project(v0), project(v1)}` 两点 pcurve，并把边的曲线清为 `kNoGeometry`、
  `t0/t1` 复位为 0/1。投影失败即拒绝整个 heal。
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
