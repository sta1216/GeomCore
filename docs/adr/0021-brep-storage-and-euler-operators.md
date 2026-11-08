# ADR 0021: B-rep 扁平存储、代际 handle 与 SolidEditor 事务

- 状态: 接受
- 日期: 2026-11-08
- 关联: 架构文档 §4.6、§8.7、§12.3、§16.3；M2 计划 Task 1–Task 11

## 背景

M2 需要一套可测试、可回滚、且不依赖异常与全局 epsilon 的 B-rep 内核。架构
§8.7 要求所有拓扑变更经受控算子并维护欧拉-庞加莱不变量；§4.6 要求容差沿
vertex/edge/face/solid 显式传播。实现必须同时满足：

1. 引用稳定：拓扑实体被删除后，旧 id 不得静默指向新实体（否则测试与缓存会
   在 use-after-free 式错误上静默通过）。
2. 失败原子：任何编辑（含提交时的校验失败）不得留下半成品；失败返回
   `Errc`，不抛异常。
3. 可复制：`commit()` 需要产出独立的 `Brep` 值，编辑器可继续复用。

## 选项

1. **指针式半边图**（节点含裸指针/`unique_ptr`）：导航快；但所有权、删除后
   悬挂、深拷贝与序列化都复杂，失败回滚需要额外的 undo 日志。
2. **扁平数组 + 代际 handle + 事务副本**：实体存放在 `std::vector` 中，id 是
   `(index, generation)`；删除只递增代际；编辑在 `Brep` 副本上进行，提交时
   校验并返回该副本。
3. **直接原地修改 + 校验失败即中止**：实现最省；但失败后拓扑已部分改变，
   无法回滚，违反 §8.7 的事务要求。

## 决策

选择选项 2，并在 `include/geom/brep/topology.hpp` 中落地：

- **扁平存储**：`Vertex`/`Edge`/`Coedge`/`Loop`/`Face`/`Shell`/`SolidEntity`
  各占一个 `std::vector` 槽位；`GeometryStore` 另存 curves/surfaces/pcurves。
  遍历是连续内存，缓存友好；索引即槽位，序列化只需数组。
- **代际 handle**：`Handle` 为 `index + generation`，`kNullIndex = 0xFFFFFF`；
  `HandleAllocator::alloc()` 复用空闲槽时递增代际，`release()` 也递增代际。
  因此删除后残留的旧 id 在 `alive()` 上返回 false，不会别名到新实体；容量
  耗尽时 `alloc()` 返回 `Handle::null()`。
- **`SolidEditor` 事务**：构造函数复制输入 `Brep` 到私有 `base_` 与 `staged_`；
  所有编辑只改 `staged_`。`commit()` 先复制一份 `staged_`，再跑
  `validate_topology` 与 `validate_tolerances`，任一失败返回
  `Errc::invalid_topology`/`tolerance_violation` 且编辑器状态不变；成功则
  返回该副本。`rollback()` 用 `base_` 覆盖 `staged_`，丢弃全部未提交编辑。
- **仅实现 MEV/KEV**：M2 的基本体全部由 `add_vertex`/`add_edge`/… 直接构造，
  欧拉算子主要用于证明「拓扑变更必须经受控算子」并支撑后续编辑。`split_edge`
  （MEV）与 `join_vertex`（KEV）覆盖了增加/删除顶点与边的最小闭环；MEF/KEM
  等需要新建/合并环与面，M2 没有消费方，留待布尔与高级建模里程碑。
- **单壳实体（M2 范围限制）**：M2 不建模壳朝向，无法区分外表面与腔体内表面，
  因此每个 `solid` 恰好引用一个 shell；多壳实体（腔体、内壳）与布尔运算一起
  延后到 M3。`validate_topology` 以 `Errc::unsupported` 拒绝
  `shell_count != 1`，消息为 "multi-shell solids are deferred to M3"。

## 后果

- 正面：id 稳定且可检测失效；提交失败零副作用；无异常、无全局 epsilon；
  `commit()` 返回的 `Brep` 与编辑器完全解耦，测试可直接比较前后计数。
- 负面：每次提交复制整套数组（M2 规模可接受，大模型需增量校验）；代际
  handle 占两个字（index + generation），比裸索引大；`HandleAllocator` 的
  空闲链表使分配不是严格 O(1) 最坏情形，但 M2 未观测到热点。
- 已知限制：KEV 只接受恰好两条入射 coedge 的顶点，闭合壳上的四度顶点与
  「同一条边被一个环使用两次」的透镜环会被拒绝并返回 `invalid_topology`；
  这不是通用 KEV，后续里程碑需扩展为按环合并。相关测试
  `join rejects a vertex with four coedge uses` 固定了该边界。
- 已知限制（M2 单壳）：闭合体的欧拉特征退化为 `2 - 2g`（box 2、torus 0）；
  多壳 S/H 校验随腔体/内壳延后到 M3，届时应与布尔运算一起恢复。
- 跟进：若大模型提交成为瓶颈，引入脏区标记与增量校验；扩展欧拉算子族时
  保持同一事务与代际语义。
