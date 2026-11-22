# M3.5 性能专项（第一切片：BVH 与网格化）

## Goal

在保持现有正确性契约与可复现输出的前提下，把两个可独立验证的门禁拉回目标：

- QA-4：1M 三角形 BVH 构建 ≤ 350 ms，8M 射线 ≥ 8 Mray/s。
- QA-10：平面盒网格化 ≥ 1.0 Mtri/s（水密、朝外不变）。

QA-1（2D 布尔 ≤ 80 ms）与 QA-9（SSI ≤ 100 µs）作为第二切片，等本切片给出 profile 方法后再排。

## 基线（2026-11-22，VS2022 19.44 Release，4 核笔记本）

| 项 | 实测 | 目标 | 差距 |
| --- | --- | --- | --- |
| QA-1 union（2×50k 星形） | 1,888 ms | 80 ms | 23.6× |
| QA-4 build（1M 三角） | 1,372 ms（空载中位） | 350 ms | 3.9× |
| QA-4 raycast（8M 射线） | 0.45 Mray/s（空载中位） | 8 Mray/s | 17.8× |
| QA-10（16,666 盒） | 0.22 Mtri/s（55.6 µs/盒） | 1.0 Mtri/s | 4.5× |
| QA-9（sphere/plane） | 9.6–19.5 ms | 0.1 ms | 96–195× |

补充微基准（`qa10probe`，5000 盒）：`make_box` 82 µs/盒，`tessellate` 115 µs/盒，重复同一盒 126 µs/次——耗时在网格化路径本身，不在 fixture 构造；盒的 12 三角说明每条通用平面面路径（环/采样/三角化/展开）约 15–20 µs。

## Task 1 profile 结论（2026-11-22，Release，3 次中位）

报告 `.superpowers/perf/m3-5-profile.md`（20.9 KB）。要点：

- **QA-4 build 1372 ms**：prim 3.8% · SAH 分箱 24.1% · median split（nth_element）34.9% · 子包围盒重扫 23.2% · 其它 13.9%。单次分箱扫描同时累积 per-bin 包围盒、用前缀/后缀并集得到子盒、按 bin 划分，可去掉 nth_element 与两次重扫，预计 −45%（~720–760 ms）；SoA prim 再省 15–25%；顶层并行可到 350 ms 量级。
- **QA-4 raycast 0.45 Mray/s**：36.9 node pops/ray、4.0 prim tests/ray；traverse 占 78%。命中后 `best==1` 时剩余 908,672 次 pops **全部被 slab 拒绝**，tie 规则实际上从不触发（hits_at_tie=0）。有序子节点 + tie-aware 剪枝实测 **变慢 8–12%**，放弃。有效方向：2/4-wide 节点或内联子 AABB（每层一次 cache line，36.9→~19 touches）→ ~0.72–0.78 Mray/s；按叶子顺序重排 prims → +12–17%。报告估算 1M 点查询的 DRAM 下限约 0.7–1.0 Mray/s，8 Mray/s 需 packet/SIMD 与更好的数据布局，本轮目标改为“显著提升并记录下限”。
- **QA-10 55.6 µs/盒**：setup 1.0% · 环/coedge/采样 27.5%（edge_positions 17.2%）· uv/定向 17.7% · triangulate_polygon2 27.6% · expand/append 18.6% · cone 后处理 0.1%。每盒 **288 次 malloc/free、8688 B**；bump-arena 试验 54.9→28.8 µs（−47%）。方向：arena/out-param + 凸四边形 fan 快路径 + 盒专用发射；12 µs/盒才到 1.0 Mtri/s，需同时消掉三角化与采样开销。

## 诊断（代码阅读，Task 1 已用计时确认）

### QA-4

- `Bvh::build` 逐层对全部图元做 12-bin SAH 扫描，再对左右半区各重扫一次算包围盒；1M 图元 × 深度 ~20 意味着多次全量扫描，且 `Prim`（三角形+盒+质心）约 100 B，1M 条即 ~100 MB 的搬运与缓存压力。
- `raycast` 无子节点排序、无 tie 感知剪枝；本 benchmark 的射线全部落在同一平面（t 相同），一旦命中 `best == t`，所有 `tmin == best` 的盒仍会通过 slab 测试，遍历量接近整棵树。这是 0.25 Mray/s 的主因，而不是三角形求交本身。
- 修复方向：命中后按「盒的 `tmin == best` 且节点内最小图元索引 ≥ best_prim」剪枝（保持 tie→最小索引契约）；子节点按 `tmin` 排序；节点/图元 SoA 化以降低带宽。

### QA-10

- `tessellate_planar_face` 每个面都做：环/coedge 遍历、`tessellate_edge_positions`（每条边返回新 vector）、uv 投影、`triangulate_polygon2`、`tessellate_expand_trimmed_runs`、多次 vector 拷贝；对 4 边矩形盒这些全是固定开销。
- `tessellate()` 还做 `brep.face_ids()` 全量复制、`state.vertices/edges` 预分配与最后的 cone-rim 扫描。
- 修复方向：为「全平面 + 每面单外环 4 边 + 4 顶点」的盒加速（直接从拓扑发射 12 三角，保持现有三角形顺序与顶点索引），并把通用路径中的重复分配改为 reserve/out-param；先用计时确认哪一段占 115 µs。

## Global Constraints

- 不改公开 API、不改既有测试期望；输出必须逐位保持（canonical mesh hash、box 12/8、cylinder/torus 计数、sphere/torus 体积）。
- 所有剪枝必须保持 raycast 契约：区间 `(0, tmax]`、`t == tmax` 接受、`t` 相同取最小 prim 索引。
- 多线程只允许用于构建（子树并行），查询保持单线程、确定性；线程数不影响树形与结果。
- 新代码只在 Release 下评估性能；Debug 不要求达标，但必须编译与测试通过。
- 性能数字必须来自本机 Release 可执行文件，记录命令、构建类型与中位数（≥3 次）。

## Review Focus

1. tie-break 与 tmax 边界是否仍正确（新增共面/等距射线测试）。
2. BVH 构建并行是否确定性（同一输入多次构建的节点数组逐位相同）。
3. 盒加速路径是否逐位复现现有 12 三角形与顶点索引；非盒/带孔/曲面必须落回通用路径。
4. 是否引入未定义行为（越界、对齐、数据竞争）与 Debug 断言失败。
5. 性能收益是否可重复（报告给出 3 次中位数与基线对比）。

## Task 1：Profile harness 与热点报告（只读探针）

- Files：`.superpowers/perf/` 下的临时探针与报告（不提交产品代码）。
- 工作：用 chrono 拆解 QA-4 的 build（初始化/分箱/划分/包围盒）与 raycast（节点访问/三角形测试/命中次数）；拆解 QA-10 的单盒网格化（环构建/采样/投影/三角化/展开/追加）；产出 `.superpowers/perf/m3-5-profile.md`。
- 验收：每个阶段有毫秒数与占比，能回答「优化哪一段收益最大」。
- 门禁：Release 定向运行；不改产品代码。

## Task 2：QA-4 查询（宽节点 + 叶子序 prims）

- Files：`include/geom/index/bvh.hpp`、`tests/unit/bvh_test.cpp`。
- 工作：
  - 节点改为 2-wide：内部节点内联两个子 AABB（每层一次 64 B 行读取），叶子保持 `first/count`；或等价的 SoA 子盒数组。
  - `prims_` 按叶子顺序重排（构建完成后一次置换），使叶子测试的三角形连续。
  - 不做有序子节点/tie-aware 剪枝（profile 证明负收益）；保留 `tmin <= best` 语义与 tie 规则。
  - 可选：4 条同向射线打包的 SIMD 叶测试，只有实测有增益才保留。
- 测试：既有 bvh 全组 + 共面/等距/`t==tmax`/共享边用例；1M 共面网格基准。
- 验收：≥0.9 Mray/s（相对 0.45 至少 2×），报告给出实测与 DRAM 下限；8 Mray/s 若未达，记录差距与后续 packet 方案，不伪造。
- 结果（2026-11-23）：宽节点 + 叶子序 prims 实测 0.66–0.67 Mray/s（约 1.5×），未达 0.9；单射线延迟下限实测 0.87–0.91，packet 方案留后续；8 Mray/s 维持 gap。

## Task 3：QA-4 构建（单次分箱 + SoA + 顶层并行）

- Files：`include/geom/index/bvh.hpp`、`tests/unit/bvh_test.cpp`。
- 工作：
  - 一次扫描同时填 counts/bins，用前缀/后缀并集得到左右子盒；SAH 有效时按 bin 划分，彻底去掉 median `nth_element` 与两次子盒重扫（保留 median 回退仅当 SAH 无有效划分）。
  - prim 拆为 SoA（box/centroid/tri 三个数组），扫描只触碰所需字段。
  - 顶层若干层并行构建子树（确定性拼接 + 子索引修正）；线程数不影响 `nodes_` 字节。
- 测试：同一路径重复构建 `nodes_` 逐位相同（含跨对象重建）；串行与并行的公开查询结果必须一致（节点编号允许不同）；既有 build/raycast/nearest 全组不变。
- 验收：≤350 ms（3 次中位）；确定性测试通过。

## Task 3：QA-4 构建加速（确定性并行 + 减少重扫）

- Files：`include/geom/index/bvh.hpp`、`tests/unit/bvh_test.cpp`。
- 工作：
  - 记录每次划分的左右包围盒，避免再全扫两半；或合并到分箱扫描。
  - 顶层若干层按子树并行构建（`std::async` 或固定线程池），子树节点写入本地数组后按确定顺序拼接并修正子索引；线程数只影响速度，不影响结果。
  - 保留 12-bin SAH 与现有中位数回退规则。
- 测试：同一输入两次构建的 `nodes_` 逐位相同（新增内部可见性测试或哈希）；既有 build/raycast/nearest 全组不变；1M 三角构建基准。
- 验收：QA-4 build ≤ 350 ms（3 次中位数）；并行/串行结果逐位相同。

## Task 4：QA-10 盒网格化快路径 + arena

- Files：`include/geom/brep/tessellate.hpp`、`tests/unit/tessellate_*`、`tests/unit/brep_boolean_curved_test.cpp`（回归）。
- 工作：
  - 识别「全平面、每面单外环、4 coedge、4 顶点」的盒，直接按现有三角形顺序发射 12 三角，复用已焊接顶点索引，保持 `triangle_face` 与 `reversed`。
  - `TessellateState`/ring/polygon 的临时缓冲改为调用内复用；`tessellate_edge_positions` 返回引用或写入 out-param，消除每盒 288 次分配。
  - 不满足盒条件的输入走通用路径（同样享受缓冲复用）。
- 测试：盒输出与通用路径逐位相同（位置/三角形/triangle_face）；带孔/非矩形/曲面面仍走通用路径；canonical 哈希不变。
- 验收：≥1.0 Mtri/s（12 µs/盒）；水密、朝外、8 位置/12 三角。

## Task 5：通用平面路径的分配削减（若 Task 4 未达 1.0）

- Files：`include/geom/brep/tessellate.hpp`。
- 工作：`tessellate_edge_positions` 改 out-param 复用、`append_ring` 复用缓冲、`tessellate()` 避免 `face_ids()` 复制；保持输出逐位不变。
- 验收：QA-10 达标；所有 canonical 哈希不变。

## 退出准则

- QA-4 与 QA-10 达标（本机 Release，3 次中位数）；QA-1/QA-9 记录新基线并给出下一片计划。
- Debug/Release 全量 0 失败、0 警告；CI 六 job 全绿。
- ADR-0039（性能策略与 tie-aware 剪枝契约）、`docs/benchmarks/m3-5.md`、架构 §16.3 更新。

## 自检

- [ ] 性能数字有命令、构建类型、3 次中位数。
- [ ] tie-break/tmax 契约有新增测试。
- [ ] 盒路径与通用路径逐位比对。
- [ ] 并行构建确定性测试。
- [ ] 无 API/测试期望变更。
