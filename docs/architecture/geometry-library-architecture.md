# GeomCore：通用计算几何与 CAD 曲面内核的 C++ 库 —— 架构设计

| 项 | 内容 |
| --- | --- |
| 文档版本 | 2.1（架构评审稿） |
| 日期 | 2026-11-06 |
| 项目代号 | GeomCore（仓库 `ai_geometry`，命名空间 `geom`，宏前缀 `GEOM_`） |
| 状态 | 待评审 |
| 范围 | 架构、分层、关键决策（ADR）、契约与质量属性；**不含**实现代码与项目排期 |
| 读者 | 架构组、核心开发、构建/测试基础设施负责人、下游集成方（Python/机器人/CAD 前端/曲面建模） |

| 版本 | 日期 | 变更 |
| --- | --- | --- |
| 1.0 | 2026-10-31 | 初稿：通用计算几何 + 网格处理库 |
| 2.0 | 2026-11-06 | 主战场调整为「通用计算几何 + **CAD 曲面内核**」：新增曲线曲面、B-rep 拓扑与建模、SSI、网格化、容忍建模容差体系；能力矩阵、里程碑、风险与 ADR 重排 |
| 2.1 | 2026-11-06 | 质量保障独立成体系：新增测试策略文档，§12 增加分层、覆盖率/变异门槛、负向用例与里程碑测试出口 |

---

## 0. 摘要（TL;DR）

1. **双主战场**：通用计算几何（精确谓词、三角剖分、多边形布尔、空间索引、网格处理）与 **CAD 曲面内核**（B-rep 拓扑、NURBS/解析曲面、曲面求交 SSI、B-rep 布尔、圆角/抽壳、网格化、质量属性）。
2. **8 层单向依赖**（L0 基础设施 → L1 数值核心 → L2 图元与曲线曲面 → L3 空间索引 → L4 结构几何与网格 → L5 建模内核与组合算法 → L6 互操作 → L7 工具）；B-rep 与网格布尔同在 L5、互不依赖，仅通过「网格化契约」交互（ADR-013）。
3. **精度分层**：符号判定走「过滤 → 自适应精确 → 精确数域」漏斗；**谓词只承诺符号正确**。排列/三角剖分/多边形布尔等组合结构的拓扑一致性只依赖符号，坐标输出允许 double 近似。
4. **容忍建模是 B-rep 的一等公民**：顶点/边/面/实体带显式容差，容差参与缝合、网格化与布尔判定，但**不与精确谓词混用**；两者边界是架构契约（ADR-014）。
5. **图元/曲线曲面按标量 + 维度泛型；精度内核只在算法入口注入**，不在类型上病毒式扩散（ADR-002）。
6. **批量数据用 SoA 视图、拓扑编辑用 handle 索引、所有权显式**；查询路径零分配、零虚函数、零异常。
7. **核心零强制依赖**，C++20 基线；`geom::core` 纯头文件，B-rep/算法以编译库交付并显式实例化常见标量；**只有 C API 承诺 ABI 稳定**。
8. 质量策略是「性质测试 + 变形测试 + 差分参考 + 结构化模糊 + 不变量校验」五件套；**曲面求交与 B-rep 布尔是最高风险项**，以可诊断失败替代静默错误；网格化以「水密 + 偏差有界」为硬契约。

---

## 1. 驱动因素、目标与非目标

### 1.1 架构驱动因素（按优先级）

| 优先级 | 驱动因素 | 含义 | 落地手段 |
| --- | --- | --- | --- |
| 1 | 正确性 / 鲁棒性 | 退化输入不得产生拓扑自相矛盾的输出，也不得静默给出错误结果 | 精确谓词、显式退化策略、`Result` 诊断、不变量校验 |
| 2 | 性能 | 与专用手写代码同量级；尾延迟可预测；可扩展到大输入 | 过滤谓词、SoA、BVH/SAH、零分配热路径、确定性并行 |
| 3 | 可组合性 | 算法可被复用、组合、替换；用户可注入自定义标量/索引 | concepts 扩展点、算法为自由函数模板、显式 context |
| 4 | 可移植 / 可构建 | Windows / Linux / macOS / WASM；无缝接入 CMake/vcpkg/Conan | 零强制依赖、CMake 导出配置、编译选项矩阵 |
| 5 | 可演化 | 十年尺度的 API 与 ABI 演进 | 分层隔离、inline 版本命名空间、C API 稳定层、ADR 流程 |

### 1.2 质量属性场景（可验收）

| # | 场景 | 期望响应 | 度量 |
| --- | --- | --- | --- |
| QA-1 | 输入 10 万个顶点、含共线与重合边的多边形执行布尔并 | 输出合法多边形（无自交、洞结构正确），面积误差在阈值内 | 相对面积误差 ≤ 1e-12；单线程 ≤ 80 ms |
| QA-2 | 4 线程并发对同一只读 `BVH` 做 100 万次射线查询 | 无数据竞争，结果与串行一致 | TSan 干净；每线程结果哈希一致 |
| QA-3 | 输入含 NaN/Inf 或索引越界的三角网 | 在 API 边界被拒绝，返回 `Errc::non_finite` / `Errc::index_out_of_range` | 不崩溃、不产生部分写入 |
| QA-4 | 100 万三角形网格构建 BVH 并做 800 万条射线查询 | 构建与查询达到预算 | 构建 ≤ 350 ms（单线程）；≥ 8 M ray/s/核 |
| QA-5 | 输入 closed 2-manifold 的两网格做布尔差 | 输出 closed 2-manifold、定向一致 | 流形性/水密性校验通过；欧拉特征符合预期 |
| QA-6 | 升级 minor 版本 | 源码兼容，废弃项有迁移说明 | 编译全部示例与下游快照；`-Werror` 下无新告警 |
| QA-7 | 同一输入在 Windows/Linux/macOS + MSVC/GCC/Clang 运行 | 声明为 D0 的算法输出逐位一致 | 规范化输出哈希一致 |
| QA-8 | 两个各 1 万面的 B-rep 实体执行布尔差 | 输出有效 B-rep；网格化后水密；体积守恒 | 体积相对误差 ≤ 1e-9；单线程 ≤ 5 s |
| QA-9 | 相切/重合面/共边等退化 SSI 语料（≥ 200 例） | 要么给出拓扑可缝合的交线，要么返回可诊断失败 | 无效输出 0 例；每例附诊断码 |
| QA-10 | 10 万面 B-rep 网格化 | 水密、无翻转、偏差有界 | 弦高 ≤ 用户容差；共享边采样逐点一致 |
| QA-11 | 质量属性（体积/质心/惯性矩）对比 OCCT 参考 | 数值一致 | 相对误差 ≤ 1e-6 |
| QA-12 | 带 1e-6 缝隙的装配做缝合 | 缝隙被桥接或明确拒绝 | 缝合后容差链自洽；拒绝时给 gap 位置与量值 |

### 1.3 非目标（明确不做，避免范围蔓延）

- 不做渲染、GPU 管线、着色器；不做运行时反射与插件 ABI。
- 做 B-rep 曲面内核，但**不做参数化特征树/历史回放、装配约束求解、STEP/IGES 解析**（属于应用层或后续适配器）；不做「全局精确有理 SSI」的承诺，曲面求交以「区间认证 + 容差」交付。
- 不提供「全局 epsilon」；不承诺浮点几何在任意退化输入下的精确坐标（只承诺符号与拓扑）。
- 不强制依赖 Eigen / GMP / TBB / Boost；不做隐式精度降级。
- 不当「算法百科全书」：每个能力面必须通过性能与维护成本门槛（见 §16 里程碑出口条件）。

### 1.4 可验证假设（影响面标注）

| # | 假设 | 若不成立的影响 | 验证方式（里程碑） |
| --- | --- | --- | --- |
| A1 | 主战场是「通用计算几何 + CAD 曲面内核」双主线（2026-11-06 已确认） | 若退回单主线，可裁剪 §8.5–8.12 与 curves/surfaces/brep 模块，路线图收缩两个里程碑 | 已完成（用户确认） |
| A2 | 2D 与 3D 是主战场，N 维仅凸包/最近邻等有限能力 | L2 维度泛型需扩张，算法成本显著上升 | M1 API 评审 |
| A3 | 用户接受「拓扑精确 + 坐标近似 + 显式容差」的组合契约 | 若要求坐标也逐位精确，需引入精确坐标输出与区间类型，API 与性能大幅变化 | M1 文档评审 |
| A4 | 主流目标平台支持 SSE2/NEON，且不启用 `-ffast-math` / FTZ | 谓词误差界失效，需关闭 SIMD 路径或改精确后端 | M0 CI 矩阵 |
| A5 | 测试期允许引入 CGAL/Manifold 作差分参考（不随产品分发） | 差分测试需自建参考实现，测试成本上升 | M0 法务/合规确认 |
| A6 | 接受容忍建模：B-rep 实体带显式容差，几何精度由容差契约保证而非逐位精确 | 若要求严格精确 B-rep，需引入精确有理曲面求交，性能与工期数量级上升 | M2 设计评审 |
| A7 | 测试期可用 OCCT 作 B-rep 差分参考（LGPL，仅测试，不进产物） | 需自建参考或降低 B-rep 差分覆盖 | M0 合规确认 |

---

## 2. 架构原则

| 编号 | 原则 | 说明与强制方式 |
| --- | --- | --- |
| P1 | 单向分层 | 只能依赖更低层；`tools/check_layers` 解析 include 图，CI 阻塞违规 |
| P2 | 精度与容差显式化 | 谓词零容差（符号精确）；建模容差（B-rep 实体容差、焊接、网格化）是一等显式参数，按尺度归一化并写入结果报告 |
| P3 | 核心零依赖 | L0–L4 只依赖 C++ 标准库；可选后端经编译开关隔离 |
| P4 | 视图优先、所有权显式 | 输入用 `std::span`/view，输出用调用方缓冲区或返回值；禁止隐式拷贝大对象 |
| P5 | 热路径约束 | 查询/谓词路径：无虚函数、无堆分配、无异常抛出、无 `std::function` |
| P6 | 确定性 | 算法声明确定性等级 D0/D1/D2；默认 D1（拓扑一致），D0 需显式开启 |
| P7 | 失败可诊断 | 可失败操作返回 `Result<T, Error>`；错误携带类别与最小上下文，不抛裸异常 |
| P8 | 组合优于继承 | 算法是自由函数模板；状态放进显式 `Options`/`Context`；无类层次继承 |
| P9 | 编译期成本预算 | 单模块头文件编译 ≤ 0.5 s/TU；聚合头 ≤ 2 s/TU；算法库显式实例化 |
| P10 | 每个模块有 fitness function | 依赖规则、不变量、性能预算、ABI 检查均有可执行守护 |
| P11 | 单一几何真相源 | B-rep 是 CAD 形态的权威表示；网格由网格化派生。反向的网格→B-rep 仅作近似重建，必须标记 lossy 并给误差 |

---

## 3. 总体架构

### 3.1 逻辑视图：分层

```
┌────────────────────────────────────────────────────────────────────────────┐
│ L7  工具与验证      bench · fuzz · oracle(差分参考) · viz-debug · dashboard  │
├────────────────────────────────────────────────────────────────────────────┤
│ L6  互操作          io(STL/OBJ/PLY/glTF) · C API(geom_c) · pybind11 绑定     │
│                     Eigen 零拷贝适配 · WASM 构建产物                          │
├────────────────────────────────────────────────────────────────────────────┤
│ L5  建模内核与组合算法 B-rep：NURBS/解析曲面 · SSI · 布尔 · 圆角/抽壳/偏移 ·  │
│                     网格化 · 质量属性 · 缝合 │ mesh boolean · SDF/MC · 体素化  │
├────────────────────────────────────────────────────────────────────────────┤
│ L4  结构几何与网格  convex hull · Delaunay/CDT · polygon triangulation ·     │
│                     2D arrangement(DCEL) · polygon boolean · halfedge mesh ·  │
│                     simplification · remeshing · winding number              │
├────────────────────────────────────────────────────────────────────────────┤
│ L3  空间索引与查询  BVH · KD-tree · R-tree(STR) · uniform grid ·             │
│                     nearest/raycast/overlap · broad-phase pair reporting     │
├────────────────────────────────────────────────────────────────────────────┤
│ L2  图元与曲线曲面  Point/Vector/Box/Segment/Triangle · NURBS 曲线/曲面 ·     │
│                     解析曲面(平面/柱/锥/球/环面) · Transform · predicates      │
├────────────────────────────────────────────────────────────────────────────┤
│ L1  数值核心        scalar traits · interval/affine filter · expansions ·    │
│                     exact rational/bigfloat adapter · SoS · scale context    │
├────────────────────────────────────────────────────────────────────────────┤
│ L0  基础设施        config/compiler detect · Result/Error · precondition ·   │
│                     span helpers · handle · arena/allocator · hash · meta    │
└────────────────────────────────────────────────────────────────────────────┘
                           依赖方向：只允许向下（↓），禁止反向与跨层回路
```

### 3.2 模块目录与职责

| 模块 | 层 | 职责 | 关键类型/接口 | 交付形态 |
| --- | --- | --- | --- | --- |
| `geom/core` | L0 | 配置、错误、断言、handle、arena、元编程 | `Result<T>`, `Error`, `Handle`, `Arena` | header-only |
| `geom/numeric` | L1 | 标量 traits、区间过滤、精确谓词、精确数域、SoS | `sign_orient2d/3d`, `ExactReal` concept, `ScaleContext` | header-only（精确后端可选编译） |
| `geom/primitives` | L2 | 图元、仿射变换、基础几何函数 | `Point<T,D>`, `Box<T,D>`, `Triangle3<T>`, `Affine<T,3>` | header-only |
| `geom/curves` | L2 | 直线/圆弧/椭圆/NURBS 曲线：求值、导数、反求、弧长 | `NurbsCurve<T>`, `Arc3<T>` | header-only |
| `geom/surfaces` | L2 | 解析曲面（平面/柱/锥/球/环面）、NURBS 曲面、Bézier 片 | `NurbsSurface<T>`, `Plane`, `Cylinder`, `Cone`, `Sphere`, `Torus` | header-only |
| `geom/index` | L3 | 空间索引与查询 | `Bvh`, `KdTree`, `RTree`, `Grid` | 编译库（显式实例化） |
| `geom/cgeom` | L4 | 计算几何与网格结构 | `delaunay2`, `polygon_boolean`, `HalfedgeMesh`, `simplify` | 编译库 |
| `geom/brep` | L5 | B-rep 拓扑、建模算子、SSI、布尔、圆角/抽壳、网格化、质量属性、缝合 | `Solid`, `Face`, `Edge`, `tessellate`, `boolean_op` | 编译库 |
| `geom/algo/*` | L5 | 网格布尔、SDF/Marching Cubes、体素化、偏移 | `mesh_boolean`, `sdf`, `voxelize` | 编译库 |
| `geom/io` | L6 | 网格/多边形序列化 | `read_stl`, `write_obj`, `read_ply` | 编译库 |
| `geom/capi` | L6 | C ABI 门面 | `geom_mesh*`, `geom_solid*` 不透明句柄 | 编译库（唯一 ABI 承诺） |
| `geom/adapters` | L6 | Eigen/Python/其他后端适配 | `Eigen::Map` 互转、pybind 模块 | 编译库（可选依赖） |
| `geom/experimental` | — | 未承诺稳定的 API | 任意 | 编译库（不保证兼容） |

### 3.3 依赖矩阵（只允许 ✔）

| 使用方 \ 被依赖 | core | numeric | primitives | curves | surfaces | index | cgeom | brep | algo | io | capi | adapters |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| core | ✔ | | | | | | | | | | | |
| numeric | ✔ | ✔ | | | | | | | | | | |
| primitives | ✔ | ✔ | ✔ | | | | | | | | | |
| curves | ✔ | ✔ | ✔ | ✔ | | | | | | | | |
| surfaces | ✔ | ✔ | ✔ | ✔ | ✔ | | | | | | | |
| index | ✔ | ✔ | ✔ | | | ✔ | | | | | | |
| cgeom | ✔ | ✔ | ✔ | | | ✔ | ✔ | | | | | |
| brep | ✔ | ✔ | ✔ | ✔ | ✔ | ✔ | ✔ | ✔ | | | | |
| algo | ✔ | ✔ | ✔ | | | ✔ | ✔ | | ✔ | | | |
| io | ✔ | ✔ | ✔ | | | | ✔ | | | ✔ | | |
| capi | ✔ | ✔ | ✔ | ✔ | ✔ | ✔ | ✔ | ✔ | ✔ | ✔ | ✔ | |
| adapters | ✔ | ✔ | ✔ | | | | | ✔ | | | | ✔ |

`algo` 内部再分 `detail/`，允许同层互调，但禁止反向 include 上层。禁止在公开头文件中 include 任何第三方头（唯一例外：`adapters/` 与可选精确后端，且用编译开关隔离）。

### 3.4 运行视图：三种典型数据流

```
(1) 批处理流水线（Python/C++ 离线）
    文件 → io 解析 → 校验(有限性/索引/规模) → 索引构建 → 算法 → 不变量校验 → 序列化
(2) 交互查询（机器人/前端）
    加载一次 → 只读索引 → 每帧 N 次 raycast/overlap/nearest → 无分配、无锁
(3) 网格编辑（前端/工艺）
    halfedge 网格 → 局部编辑(handle 失效检测) → 重建索引(refit/rebuild) → 布尔 → 结果校验
(4) 曲面建模（CAD 内核）
    B-rep Solid → 建模算子(布尔/圆角/抽壳/偏移) → SSI + 面分割 → 修复/缝合
              → 校验(欧拉-庞加莱/水密/定向) → 网格化(水密) → 渲染/仿真/下游
```

### 3.5 部署视图

| 形态 | 内容 | 说明 |
| --- | --- | --- |
| 头文件库 | `geom::core`（L0–L2） | 轻量集成，编译期即可用；不含重量算法 |
| 静态/动态库 | `geom::brep` + `geom::algo` + `geom::io` + `geom::capi` | 显式实例化 float/double；控制编译时间与二进制体积 |
| CAD 宿主集成 | `geom::capi` 暴露 B-rep 不透明句柄；STEP/IGES 由外部适配器转换 | 内核不内置 CAD 交换格式解析（见 §1.3） |
| WASM | 同上（无异常、无线程变体） | 供 Web 端几何预览 |
| Python wheel | pybind11 绑定 + buffer protocol | 零拷贝接 NumPy |

### 3.6 4+1 视图映射

| 视图 | 关注 | 章节 |
| --- | --- | --- |
| 逻辑 | 分层、模块、依赖 | §3.1–3.3 |
| 运行 | 数据流、线程、内存所有权 | §3.4, §6, §9 |
| 开发 | 目录、构建、编码规范 | §11, 附录 A |
| 部署 | 库形态、包管理、跨平台 | §3.5, §11 |
| 场景 | 质量属性验收 | §1.2, §12 |

---

## 4. 数值与鲁棒性架构（核心）

几何库的成败首先取决于「退化输入下是否还讲道理」。本章是全部算法的公共地基。

### 4.1 三层精度漏斗

```
sign(orient2d(a,b,c))
 ├─ A  float 过滤：float 计算 + 静态误差界 │err│ > bound ? ── 符号可信 ──▶ 返回   （经验占比最大）
 ├─ B  double 过滤：1× 求值 + 误差界                        ── 符号可信 ──▶ 返回
 ├─ C  自适应精确：Shewchuk 展开式（误差自由变换）           ── 符号可信 ──▶ 返回
 └─ D  精确数域：有理数/大浮点（可选后端）                    ────────────▶ 返回
```

| Tier | 实现 | 触发条件 | 相对成本 | 说明 |
| --- | --- | --- | --- | --- |
| A | `float` 求值 + `static errbound` | 绝大多数输入 | 1× | 以 float 求值（需处理上溢/下溢）；标量本身为 `float` 时跳过本层 |
| B | `double` 求值 + 误差界 | A 失败 | ~1.5× | 主流路径，通常覆盖 >99% |
| C | 误差自由变换 + 展开式求和 | B 失败（近退化） | ~10–50× | 只计算符号，按需扩展 |
| D | `ExactReal` 适配（有理数/大浮点） | C 失败或用户强制 | 100–1000× | 编译期开关，默认仅在 arrangement/布尔分类启用 |

设计约束：
- 误差界必须由谓词实现自己声明，不依赖编译器优化假设；实现文件强制 `#pragma`/编译选项关闭浮点收缩（`-ffp-contract=off`，MSVC `/fp:precise`）。
- 32 位 x86 必须走 SSE2（`-mfpmath=sse`），禁止 x87 扩展精度进入误差界推导。
- 禁止在库内全局修改舍入模式/FTZ/DAZ；若用户外部启用，文档明确「精确性保证失效」。
- FMA 允许使用，但误差界按「有 FMA」与「无 FMA」两套常量分别推导，并在 CI 两种编译下验证符号一致。

### 4.2 谓词契约

| 谓词 | 语义 | 承诺 |
| --- | --- | --- |
| `sign_orient2d(a,b,c)` | 三角形有向面积符号 | 符号 100% 正确；量值不保证 |
| `sign_incircle(a,b,c,d)` | d 相对过 a,b,c 圆的内外 | 同上 |
| `sign_orient3d(a,b,c,d)` | 四面体有向体积符号 | 同上 |
| `sign_insphere(a,b,c,d,e)` | e 相对四面体外接球内外 | 同上 |
| 派生：线段相交符号、点在线/面侧、共面判定 | 由上述组合 | 同上 |
| 距离/最近点/量值类 | 浮点近似 | 相对误差 ≤ 8 ulp（单次运算链，文档给出界） |

规则：**符号函数绝不接收 epsilon**；量值函数绝不参与拓扑决策。

### 4.3 退化处理策略：两种手段，分开使用

| 场景 | 手段 | 理由 |
| --- | --- | --- |
| 组合结构（arrangement、Delaunay、三角剖分） | 符号扰动 SoS（Edelsbrunner–Mücke），按元素 id 的字典序打破平局 | 结构必须无歧义；坐标不可移动，但「谁在前」可由确定性规则决定 |
| 用户可见几何（布尔、偏移、网格） | 显式退化分类 + 规则化（regularized）输出 + 诊断报告 | 不得静默移动用户顶点；退化必须体现为「边/面被合并」或明确错误 |

SoS 的扰动深度上限固定（如 3 层），并记录到 `SolveReport.sos_triggered`，便于复现与线上定位。

### 4.4 容差契约（Tolerance Contract）

1. **谓词容差 = 0**：不存在隐式 epsilon。
2. **合并/焊接容差**（`weld(eps)`、`snap_rounding`）：必须显式传参；内部按输入包围盒对角线归一化，绝对 epsilon 只作为用户显式请求的语义。
3. **输出简化容差**（`simplify(tol)`、`offset(miter_limit)`）：显式，且写入结果报告（实际最大偏差）。
4. **点集包含的边界策略** `ContainsPolicy { closed, open, boundary_touching }`：用枚举表达语义，不用 epsilon 模拟。
5. 所有接受容差的 API 在文档中给出「容差如何相对模型尺度解释」的公式与示例。

### 4.5 一致性模型：拓扑精确 + 坐标近似

```
输入(近似坐标) ──▶ 精确符号 ──▶ 组合结构(精确、可信) ──▶ 坐标输出(double 近似)
                        ▲                                        │
                        └──── 任何拓扑决策都不得使用近似坐标 ────────┘
```

明确否决的路径：
- 先舍入到整数网格再做几何（Clipper 式）：会改变输入语义，且整数溢出/尺度依赖难以通用。
- 先算 double 交点再判符号：误差会级联成自交与裂缝。
- 用统一 epsilon 吸收近退化：在尺度差异大的模型上必然失效。

### 4.6 容忍建模（B-rep 的容差体系）

B-rep 不是「精确几何」，而是**带显式容差的近似模型**。容差不是全局 epsilon，而是绑定到实体上的精度声明：

| 字段 | 语义 | 典型初始化 |
| --- | --- | --- |
| `Vertex::tol` | 顶点到其所连边/曲线端点的最大允许偏差 | 1e-7 × 模型尺度 |
| `Edge::tol` | 边到其三维曲线/pcurve 的最大允许偏差 | ≥ 两端点 tol |
| `Face::tol` | 面到其支撑曲面的最大允许偏差（含环内裁剪曲线偏差） | ≥ 相邻边 tol |
| `Solid::tol` | 实体级精度上界，用于装配与交换 | = max(子实体 tol) |

一致性不变量（建模操作后必须校验，失败即回滚到操作前快照）：
1. 顶点到所连边的距离 ≤ `max(vertex.tol, edge.tol)`；
2. 边曲线与两侧面的 pcurve 在各自参数域上的偏差 ≤ `face.tol`；
3. 相邻面的公共边两侧容差带必须相交（可缝合）；
4. `Solid::tol ≥ max(所有子实体 tol)`；
5. 容差只能被显式操作放大，且每次放大记录来源与量值。

容差传播：布尔/求交/圆角产生的新实体，容差 = `max(输入容差, 数值误差界) + 受控余量`；增长超过输入容差 **10×** 时上报 `ToleranceGrowth` 诊断，禁止静默放大。

反模式（评审必拒）：用放大容差「修复」布尔失败、全局 epsilon 缝合、静默移动顶点。

### 4.7 曲面求交的数值基础

目标：求交结果不追求逐位精确，但必须**可认证、可缝合**。

- **区间/仿射算术**：对 patch 做包含界与符号判定，区间分离即剪枝；
- **Bézier 提取 + 凸包界**：控制点包围盒比参数域盒更紧；有理曲面在齐次坐标（4D）中取凸包再投影；
- **递归细分**：平坦性判据 = 控制点偏离平均平面 ≤ `tol_sub`；落到候选片对后再精化；
- **精化**：区间 Newton / 二分，收敛到目标容差；残差 `r = max(dist(p, S_a), dist(p, S_b))` 必须记录；
- **退化分类**：切向接触（区间无法分离但梯度秩 < 2）、重合面、边界穿越，分别走「保留接触边 / 2D 参数域布尔 / 诊断」；
- **认证输出**：每段交线输出 `{pcurve_a, pcurve_b, tol, residual}`，上层缝合只依赖这些字段，不重新猜测几何；
- **与精确后端的关系**：允许用有理区间做符号认证，但**不承诺全局精确有理交线**（ADR-015）。

### 4.8 数值失败的可观测性

```cpp
struct NumericStats {          // 由调用方选择收集（默认零开销）
  uint64_t filtered   = 0;     // 走过滤路径的谓词数
  uint64_t adaptive   = 0;     // 走自适应精确
  uint64_t exact      = 0;     // 走精确后端
  uint64_t sos        = 0;     // 触发符号扰动
  uint64_t give_up    = 0;     // 精度耗尽（应恒为 0，否则是缺陷）
};
```

- 每个算法可选返回 `SolveReport`（含 `NumericStats`、退化计数、最大残差）。
- `give_up > 0` 一律视为缺陷：转成 `Errc::numeric_failure` 并附带复现用例的最小要素。

---

## 5. 核心类型与 API 设计

### 5.1 代数空间：Point / Vector / Direction / Normal

| 类型 | 语义 | 关键约束 |
| --- | --- | --- |
| `Vector<T,D>` | 自由向量 | 可加、可缩放、可点积/叉积（D=3） |
| `Point<T,D>` | 仿射点 | `Point + Vector = Point`；`Point − Point = Vector`；**禁止 Point + Point** |
| `Direction<T,D>` | 单位方向 | 只能由 `normalize(v)` 与 `Direction::from(v)` 构造，构造后模长 = 1（误差界内） |
| `Normal<T,D>` | 单位法向 | 与 Direction 同构但语义不同（避免方向/法向混用） |
| `Frame<T,3>` | 正交标架（原点 + 三轴） | 用于局部坐标与 OBB |
| `Affine<T,3>` | 4×3 仿射（平移+线性） | 逆变换返回 `Result`（奇异矩阵显式失败） |

设计取舍：区分 `Point`/`Vector` 能静态拦截一类经典错误（点加点、法向当点用），代价是少量样板；`Direction`/`Normal` 分离仅为防语义误用，两者可 `reinterpret` 级别零成本互转（显式函数）。

### 5.2 图元与不变量

| 类型 | 表示 | sizeof（double, 64-bit） | 允许退化 | 不变量 |
| --- | --- | --- | --- | --- |
| `Point2d/3d` | 2/3 × T | 16 / 24 | — | 坐标有限（边界处校验） |
| `Box<T,D>` | lo, hi | 32 / 48 | ✔（空盒 `lo > hi` 表示空） | 构造后 `lo ≤ hi` 或明确为空 |
| `Segment<T,D>` | 2 点 | 32 / 48 | ✔（零长线段） | 无 |
| `Line<T,D>` | 点 + 方向 | 32 / 48 | ✘ | 方向非零 |
| `Ray<T,D>` | 原点 + 方向 | 32 / 48 | ✘ | 方向非零 |
| `Plane<T,3>` | 法向 + 偏移（单位法向） | 32 | ✘ | 法向单位化 |
| `Triangle<T,3>` | 3 点 | 72 | ✔（面积 0） | 无；退化由算法显式处理 |
| `Tetrahedron<T,3>` | 4 点 | 96 | ✔ | 无 |
| `AABB` 即 `Box<T,D>` | | | | |
| `Sphere<T,3>` | 球心 + 半径 | 32 | ✔(r=0) | r ≥ 0 |
| `Circle<T,2>` | 圆心 + 半径 | 24 | ✔ | r ≥ 0 |

退化图元是**合法值**，不是错误；错误只出现在「算法要求非退化」时，且必须是 `Result` 而非断言崩溃（未定义行为仅限文档明确声明的 precondition 违反）。

### 5.3 泛型参数轴

```cpp
namespace geom::v1 {

template <class T> struct scalar_traits;     // 内置 float/double/long double；自定义标量经特化接入

template <class T> concept Scalar =
    requires(T a, T b) {
        { a + b } -> std::convertible_to<T>;
        { a / b } -> std::convertible_to<T>;
        { std::isfinite(a) } -> std::convertible_to<bool>;
        typename scalar_traits<T>;
        typename std::numeric_limits<T>;
    };

template <class T, int D> struct Point;      // D = 2/3 为主；N 维算法复用 D > 3 的实例化
template <class T, int D> struct Vec;
template <class T, int D> struct Box;

// 精度策略只在算法入口注入
struct KernelFast;      // 纯 double，不启用精确回退（明确不保证退化正确性）
struct KernelRobust;    // 默认：过滤 + 自适应精确
struct KernelExact;     // 全精确（慢，用于对照与离线校验）

} // namespace geom::v1
```

三条泛型轴互不耦合：标量（`T`）→ 图元；维度（`D`）→ 图元与算法；精度策略（`K`）→ 仅算法。**图元不携带内核参数**，避免 CGAL 式类型爆炸（ADR-002）。

### 5.4 API 命名与返回语义

| 语义 | 命名 | 返回 | 说明 |
| --- | --- | --- | --- |
| 布尔判定 | `intersects(a,b)`, `contains(x,p)`, `is_degenerate(x)` | `bool` | 永不抛异常（输入已校验时） |
| 构造交 | `intersection(a,b)` | `IntersectionResult`（variant + kind 枚举） | 「无交」与「退化交」必须可区分，不能用 `optional` 糊在一起 |
| 距离 | `distance(a,b)`, `squared_distance(a,b)` | `T` | 优先 `squared_distance`，避免无谓开方 |
| 最近点 | `closest_point(a, p)` | `Point` | |
| 投影 | `project(seg, poly, mesh, p)` | `Point` / 求交结果 | |
| 包围盒 | `bounds(x)` | `Box` | |
| 变换 | `transform(x, M)`；`inverse(M)` | 值 / `Result` | 变换后自动更新退化状态 |
| 归一化 | `normalize(v)` | `optional<Direction>` | 零向量返回 nullopt，不产生 NaN |
| 精确化 | `to_exact(x)` / `to_double(x)` | 显式转换 | 禁止隐式窄化 |

API 纪律：
- 所有公开函数 `const`、`noexcept`（除非明确可能失败）；
- 返回 `Result` 的函数不写输出参数；返回大对象的函数支持调用方缓冲区重载；
- 重载族必须语义一致（同一名字在不同维度上含义相同）；
- 每个头文件自包含、可单独编译（include-what-you-use 检查）。

### 5.5 错误模型

```cpp
enum class Errc : uint8_t {
  none, degenerate, not_manifold, self_intersecting, non_finite,
  index_out_of_range, numeric_failure, io_error, unsupported, cancelled
};

struct Error {
  Errc code;
  const char* message;        // 静态字符串，零分配
  uint32_t context_a;         // 如元素 id / 行号 / 谓词编号
  uint32_t context_b;
};

template <class T> using Result = /* std::expected<T, Error>（C++23）或等价 shim */;
```

- 异常策略编译期可选：`GEOM_ERROR_MODE = expected | throw | abort`；库内部**只在入口**转换，核心路径不抛。
- 断言分两级：`GEOM_PRECONDITION`（内部逻辑不变量，debug 触发）与 `GEOM_EXPECT`（用户输入契约，debug/release 都检查热点外的边界）。
- 用户输入校验集中在 I/O、C API、pybind 边界与接受 `span` 的批处理入口（每 span 一次，向量化实现）。

### 5.6 版本化命名空间

```cpp
namespace geom { inline namespace v1 { ... } }
```

- 破坏性变更 → `v2`，旧版以 `geom::v1` 别名保留至少两个 minor；
- 废弃：`[[deprecated("use X")]]` + 迁移文档 + 在 CHANGELOG 标注移除版本；
- 实验 API 一律放 `geom::experimental`，不承诺兼容。

### 5.7 曲线与曲面类型

| 类型 | 表示 | 关键不变量 | 存储 |
| --- | --- | --- | --- |
| `Segment` / `Arc` / `Ellipse` / `Line` | 解析 | 参数域有效、半径 > 0 | 定长（≤ 64 B） |
| `NurbsCurve<T>` | 阶数 + 控制点 + 权重（有理时）+ 单调节点向量 + 定义域 | 节点非减；端点重数 ≤ p+1；权重 > 0；无零长度节点区间 | 类型 32 B 视图 + arena 中的控制点/节点 |
| `BezierCurve<T>` | 控制点 | — | 视图 |
| `Plane` | 单位法向 + 偏移 | 法向单位 | 32 B |
| `Cylinder` / `Cone` / `Sphere` / `Torus` | 轴/原点/半径/半角 | 轴单位、半径 > 0、半角 ∈ (0, π) | ≤ 64 B |
| `NurbsSurface<T>` | 双阶数 + 控制网格 + 权重 + 双节点向量 | 同曲线；允许退化边（极点）但必须显式标记 | 视图 + arena |
| `BezierPatch<T>` | 控制网格 | — | 视图 |

原则：**解析曲面永远保留类型身份**，不被泛化为 NURBS（求值与求交都更精确、更快）；NURBS 是统一的后备表示（ADR-016）。

### 5.8 B-rep 拓扑实体

| 实体 | 关键字段 | 不变量 |
| --- | --- | --- |
| `Vertex` | point, tol | tol ≥ 0；坐标有限 |
| `Edge` | curve, [t0,t1], v0, v1, tol | 端点落在 curve 上（≤ tol）；t0 < t1；无零长度参数区间 |
| `Coedge` | edge, sense, loop, pcurve, tol | 闭合实体中每条边恰有 2 个 coedge；pcurve 与 3D 曲线在 face 容差内一致 |
| `Loop` | coedge 环, face, 外/内标记 | 环闭合；每个 face 恰有一个外环；内环互不相交 |
| `Face` | surface, loops, sense, tol | 环在参数域内合法；法向与 sense 一致；面积有限 |
| `Shell` | face 集合, 开/闭 | 闭合 shell 每条边恰有 2 个 coedge |
| `Solid` | 外壳 + 内腔 shell 集合, tol | 外壳外向、内腔反向；壳之间不相交 |

存储：拓扑是**扁平数组 + 32 位代际 handle**（§6.4），不是指针图；曲线/曲面/pcurve 在 arena 中按实体索引存放（§6.8）。尺寸估算（含对齐，不含可变长几何）：Vertex 40 B、Coedge 40 B、Edge 72 B、Face 96 B、Loop 32 B、Shell 48 B、Solid 64 B。

## 6. 数据结构与内存架构

### 6.1 值语义与布局

- 图元是 **trivially copyable 的聚合体**，可 `memcpy`、可 `constexpr`（在标量支持时）；
- 热路径按值传递小对象（≤ 64 B），大对象（网格、多边形集）用视图；
- 对齐：`Point3d` 8 字节对齐；SIMD 批量视图要求 16/32 字节对齐（运行时检查，失败走标量路径）。

### 6.2 SoA 视图与批量接口

```cpp
template <class T> struct PointCloudView { std::span<const T> x, y, z; };
template <class T> struct TriangleSoupView {
  PointCloudView<T> pos;                    // 顶点位置
  std::span<const std::array<uint32_t,3>> tris;
};

// 批量查询：回调式，避免中间分配
template <class K, class F>
void bvh_query(const Bvh<K>& b, std::span<const Ray<typename K::scalar,3>> rays,
               F&& on_hit);
```

- AoS ↔ SoA 互转是显式操作（`to_soa`/`to_aos`），不在算法内部隐式转换；
- 批量 API 以「一次校验 + 向量化内核」为原则，单个 API 调用不做 N 次校验。

### 6.3 容器策略

| 需求 | 结构 | 理由 |
| --- | --- | --- |
| 多边形集合（带洞） | `PolygonSet`：CSR（`rings` 偏移数组 + 顶点数组 + `ring→polygon` 映射） | 一次分配、缓存友好、无 `vector<vector<>>` 碎片 |
| 多段线集合 | `PolylineSet`（同上 CSR） | 同上 |
| 小规模临时点集 | `SmallVector<Point,16>` | 避免布尔运算中的反复堆分配 |
| 网格 | `HalfedgeMesh` + 独立位置数组 | 拓扑与几何分离，支持只读几何的零拷贝替换 |
| 索引 | 扁平数组 + 显式空闲表（可选） | 可序列化、无指针 |

### 6.4 拓扑句柄（handle）设计

```cpp
// 32 位打包：24 位索引 + 8 位代数（防 ABA）
struct Handle {
  uint32_t packed;
  static constexpr uint32_t kIndexBits = 24;      // 单网格上限 16,777,215 个元素
  static constexpr uint8_t  kIndexMask  = 0xFFu;
  uint32_t index() const { return packed >> 8; }
  uint8_t  generation() const { return packed & 0xFFu; }
};
static_assert(sizeof(Handle) == 4);
```

- 元素删除时代数 +1；`is_valid(h)` 检出悬垂句柄；
- 溢出（超过 24 位）在插入时返回 `Errc::index_out_of_range`，绝不静默回绕；
- 句柄不暴露给 C API 之外的持久化层；序列化用「索引 + 代数」并重建。

### 6.5 Halfedge 网格布局

| 字段 | 偏移 | 类型 | 说明 |
| --- | --- | --- | --- |
| next | 0 | u32 | 同面环上的下一条半边 |
| twin | 4 | u32 | 对偶半边 |
| origin | 8 | u32 | 起点顶点 |
| face | 12 | u32 | 所属面 |
| flags | 16 | u16 | 边界/删除/属性位 |
| reserved | 18 | 6 B | 对齐填充，未来扩展（如属性句柄） |

- `sizeof(Halfedge) = 24`；100 万三角形 ≈ 300 万半边 ≈ 72 MB；
- 顶点：位置（SoA 主存储）+ `halfedge`（4 B）+ padding → 32 B/顶点；
- 面：`halfedge` + `flags` = 8 B；
- 总内存（1M 三角形）≈ 72 + 16 + 8 + 索引开销 ≈ **100 MB**，满足 §13 预算；
- 优化：叶级顶点位置可外置（只读网格用 `TriangleSoupView` 零拷贝），编辑时再物化。

### 6.6 分配器

- 所有批量临时内存走 `std::pmr::monotonic_buffer_resource`（`geom::Arena`），单个算法调用一个 arena，出作用域整体释放；
- 用户可注入 `std::pmr::memory_resource`；默认 `new_delete_resource`；
- 查询/谓词路径无分配（静态断言 + 测试挂钩统计分配次数，CI 门禁）。

### 6.7 序列化

- 二进制格式 `GEOM`：magic + 版本 + 端序标记 + 段表（`[type, offset, count, stride]`）+ 对齐填充；
- 只序列化值/索引，不序列化指针与句柄代数；
- 读取时做「长度、偏移、索引范围、有限性」全量校验（防不可信输入）；
- 可选 glTF/PLY 互转由 `io` 层完成，核心格式保持最小。

### 6.8 B-rep 存储、属性与事务

- **拓扑**：扁平数组 + 32 位代际 handle（§6.4）；删除走空闲表 + 代数递增，数组可整体迁移/序列化，无指针。
- **几何**：曲线/曲面的控制点、节点向量存放在 arena，按几何句柄索引；实体只保存「几何句柄 + 参数域」，避免拷贝。
- **pcurve**：每个 coedge 一条 2D 参数曲线；与 3D 边的参数化不要求一致，用段级映射表建立对应，缝合时按参数区间比对。
- **属性系统**：`Attributes` 为稀疏表（`(handle, type_id) → value`），支持 int/double/string/POD；**属性不参与几何运算**，仅随实体做增删改查与序列化。
- **事务**：建模操作在 `SolidEditor` 上进行（staged changes），提交时生成新 `Solid`；失败自动回滚，保证「要么有效、要么原样」。
- **序列化**：B-rep 段在 §6.7 格式上扩展（拓扑数组、几何 arena、pcurve 表、属性表）；读取时强制校验容差链（§4.6）与索引范围。

## 7. 空间索引与查询子系统

### 7.1 统一索引契约

以下以 `double`/3D 为例；实际 concept 对 `T`、`D` 泛型，且查询回调是模板参数（避免 `std::function` 间接调用）。

```cpp
template <class I>
concept SpatialIndex = requires(const I& i, Box3d b, Ray3d r, Point3d p) {
  { i.bounds() }        -> std::same_as<Box3d>;
  { i.size() }          -> std::convertible_to<size_t>;
  i.query(b, /* callback */);   // 范围
  i.query(r, /* callback */);   // 射线
  i.closest(p);                 // 最近图元
};
```

实现按「构建一次、只读并发查询」优化；动态更新用独立的 `refit`/`rebuild` 策略，不做隐式增量维护。

### 7.2 BVH（主力）

| 决策 | 选择 | 理由 |
| --- | --- | --- |
| 构建 | 分箱 SAH（16 箱，沿最长轴），O(n log n) | 质量/成本平衡；可并行分区 |
| 并行构建 | 按 Morton 序确定性分区，两阶段（分区 → 并行子树） | 结果与线程数无关（D0） |
| 节点 | 4 宽（BVH4）为默认，标量二叉为调试/对照 | 4 宽在遍历上常有 1.5–2× 收益 |
| 射线-三角形 | 水密（watertight）求交，与遍历的 slab 测试共用同一舍入语义 | 避免网格裂缝导致的漏检 |
| 查询 | 显式栈（固定 128 深度）+ 有序子节点近邻优先 | 无递归、可预测 |
| 最近邻 | 带优先队列的 best-first + 剪枝 | 精确最近，非近似 |
| 更新 | `refit()`（O(n) 自底向上）与 `rebuild()`；调用方选择 | 显式成本，避免隐式抖动 |
| 序列化 | 扁平节点数组 + 索引 | 可 mmap 只读加载 |

### 7.3 其他索引

| 索引 | 适用 | 关键参数 | 备注 |
| --- | --- | --- | --- |
| KD-tree | 点云最近邻、范围查询 | 叶大小 16–32、SAH 切分 | 3D 点云优于 BVH |
| R-tree（STR 批量装载） | 2D 多边形/矩形集合、GIS 式查询 | 节点填充率 70% | 动态插入走 R*-tree 变体（后置） |
| Uniform grid / 空间哈希 | 均匀尺度、稠密场景 | 单元格 = 平均尺寸 × 2 | 构建 O(n)，最坏退化（尺度差异）需回退 BVH |
| Broad-phase | 两集合相交对、自交检测 | 复用 BVH 双树遍历 | 输出对按 (id_a, id_b) 排序，保证确定性 |

### 7.4 复杂度与预算

| 操作 | 复杂度 | 预算（单线程，1M 三角形） |
| --- | --- | --- |
| SAH 构建 | O(n log n) | ≤ 350 ms |
| refit | O(n) | ≤ 25 ms |
| 射线查询 | O(log n) 期望 | ≥ 8 M ray/s/核 |
| 最近邻 | O(log n) 期望 | ≤ 2 µs/查询（k=1） |
| 自交 broad-phase（1M 三角） | 输出敏感 | ≤ 400 ms（无交时） |

### 7.5 面向曲面的索引与内外判定

| 需求 | 结构 | 说明 |
| --- | --- | --- |
| 面/片配对（SSI 候选） | 面级 BVH + 容差扩展盒 | 用控制凸包与曲面包围盒（含 `tol`）建树；双树遍历输出候选片对 |
| 最近点查询 | 面 BVH + 逐面 Newton | 剪枝到候选面后做参数域反求（§8.5），用残差认证 |
| 点内外判定 | 网格化 BVH 射线奇偶 + 精确符号 | 先网格化（§8.10）求射线奇偶；近边界/相切时用 winding number + 精确符号复核 |
| 干涉检查 | 双 Solid 的 BVH | 输出面/边对与最小距离，不做布尔 |

关键决策：索引只存包围盒与图元 id，不持有几何；几何变化必须显式 `rebuild`/`refit`；容差扩展是建树参数（显式传参）。

---

## 8. 算法子系统

每个算法必须以「输入契约 → 输出契约 → 鲁棒性策略 → 复杂度 → 诊断」五段式定义；以下是各子系统的架构决策。

### 8.1 凸包（2D/3D/N 维）

- 2D：Andrew 单调链，O(n log n)，共线点策略显式参数（保留/剔除/仅端点）；
- 3D：增量法 + 精确 `orient3d`，退化共面输入走「先投影到平面再做 2D」的快速路径；
- N 维：增量 + 精确 `insphere` 可选；输出为面索引（`vector<array<uint32_t,D>>`）；
- 无精确后端时，3D 凸包在近共面输入上**必须返回诊断**而非错误拓扑。

### 8.2 三角剖分

| 算法 | 输入 | 关键决策 | 鲁棒性 |
| --- | --- | --- | --- |
| 简单多边形三角剖分（含洞） | 简单多边形（可带洞） | 单调多边形分解 + 耳切；洞用桥接边合并 | 精确 `orient2d` 判凸；退化（共线）走「合并顶点」规则 |
| 2D Delaunay | 点集 | 增量插入 + 行走定位（last-face）+ SoS | 精确 `incircle`；重复点先按容差策略显式去重 |
| 约束 Delaunay（CDT） | 点集 + 约束边 | 约束插入 + 腔体重新三角化（Sloan） | 约束相交必须显式报错或先求交拆分 |
| 3D Delaunay | 点集 | 后置（M4+） | 退化（共面/共球）风险高，默认要求精确后端 |

### 8.3 2D 布尔与排列

架构选择：**基于 arrangement（DCEL）的扫描线 + 精确 winding 分类**，而非 Vatti/Greiner–Hormann 类直接裁剪。

```
输入多边形集(可自交) 
   → 边片段化(数值精确求交, 坐标近似) 
   → 扫描线建 arrangement(精确 orient 与 SoS 打破平局)
   → 每个面计算 winding number(精确符号累加)
   → 按 fill rule(non_zero / even_odd / positive) 分类保留面
   → 提取边界 → 合并共线边 → 输出 PolygonSet
```

- 支持 `union/intersection/difference/xor` 与自定义谓词（用户可表达「A 且非 B」的任意组合）；
- 自交输入默认允许（arrangement 天然处理），`allow_self_intersections=false` 时先做 O(n log n) 自交检测并报 `Errc::self_intersecting`；
- 整数模式：`snap_rounding` 到 `int64` 网格后走精确整数运算（对标 Clipper 语义），作为显式选项而非默认；
- 输出规则化：相邻共线边合并、零面积环剔除、洞归属校验；报告 `BooleanReport{input_edges, output_polygons, merged_vertices, max_residual}`。

否决：Vatti（边界情形多、难证明正确）、Greiner–Hormann（退化与洞支持差）、纯 double 扫描线（近退化必然错）。

### 8.4 2D Arrangement（DCEL）

- 半边结构，面、边、顶点三层；支持 overlay（两个 arrangement 的叠加）；
- 只读查询：面定位（walk + 随机化加速）、包含判定、遍历；
- 与布尔共用同一套扫描线与精确符号，避免两套退化逻辑。

### 8.5 曲线/曲面求值

- **求值**：NURBS 用 de Boor（O(p²)）；导数用解析公式；解析曲面走闭式；Bézier 片用 de Casteljau。求值路径零分配、标量模板化。
- **包围**：控制点为凸包的紧包围；有理情形在齐次 4D 中取包围后投影；细分后包围更紧（§4.7 用其剪枝）。
- **反求（point inversion）**：采样/网格给初值 + Newton 迭代，区间认证收敛；失败返回候选参数集而非错误点。
- **弧长/面积**：自适应 Gauss–Legendre + 误差界；按导数零点分割单调区间。
- **剪裁对应**：pcurve（2D 参数）与 surface 参数化必须显式一致；跨面比较统一在「面参数域 ↔ 3D」双射上做。

### 8.6 曲线与曲面求交（SSI）

```
SSI(S_a, S_b, tol):
  面片 BVH 配对(容差扩展) → 候选片对
  → 区间/仿射包含界剪枝（分离即丢弃）
  → 递归细分至平坦/区间可分离
  → 参数域 Newton/二分精化 → 交点 + 残差 r
  → 交线追踪（自适应步长，切向时降速）
  → 端点匹配 / 边界求交 / 拓扑拼接
  → 输出: 有序交线段 + {pcurve_a, pcurve_b, tol, r}
  失败 → 诊断（退化类型、位置、候选集）
```

| 情形 | 处理 | 输出 |
| --- | --- | --- |
| 横截（transversal） | 细分 + Newton，主路径 | 交线段 |
| 切向接触 | 区间不可分离；保留接触边/点 | 接触曲线 + 容差标注 |
| 面相重合 | 同曲面或参数重映射识别 → 2D 参数域布尔 | 面片边界 |
| 边界穿越 | 与 face 环求交，端点落在边上 | 边界交点 + 边拆分请求 |

保证：每段交线残差 ≤ `tol`、两端点一致、可被缝合消费；**不承诺精确有理交线**（ADR-015）。曲线-曲线用 3D 最近点迭代 + 区间；曲线-曲面把曲线代入曲面反求参数域。

### 8.7 B-rep 拓扑与建模算子

- **Euler 算子**（MEV/MEF/KEF…）维护欧拉-庞加莱不变量（§12.3）；**所有拓扑变更必须经受控算子**，禁止直接改数组。
- **基本体**：box/cylinder/cone/sphere/torus 的解析构造（带容差初始化）。
- **扫掠类**：extrude（方向/直线）、revolve（轴 + 角度）、loft（截面序列）、sweep（路径 + 标架传播 + 扭转）；扫掠前做自交检测（网格化粗检 + SSI 确认）。
- **修剪/分割**：交线 → 面参数域 2D 布尔（§8.3 arrangement）→ 新环 → 新 coedge/pcurve；环方向由面 sense 决定。
- **变换**：刚体、均匀缩放、镜像（镜像必须翻转 sense 并重算环方向）。
- **事务**：每个算子在 `SolidEditor` 上执行，提交前 `validate`；失败回滚并返回 `ModelingReport{new_faces, tol_growth, diagnostics}`。

### 8.8 B-rep 布尔

```
boolean(A, B, op):
  面配对(BVH + tol 扩展) → 面-面 SSI(§8.6) → 交线
  → 参数域分割面(§8.7) → 新环/coedge
  → 分类: 网格化射线奇偶 + winding/精确符号（边界模糊时）
  → 保留/丢弃面片 → 缝合（coedge 配对，容差带相交）
  → 修复(§8.12) → 校验(§12.3)
  → Solid 或 诊断
```

- **保证**：输入为有效 closed solid → 输出为有效 closed solid，且 `Solid::tol` 显式；失败必须可诊断，**绝不输出拓扑非法的 B-rep**。
- **共面重合**：走 2D 参数域布尔；**相切**：保留接触边并标注；**多壳/内腔**：逐壳分类后重组。
- **不承诺**：与商业内核相同的面片细分/拓扑；只承诺几何在容差内等价。
- 性能目标：1 万面输入 ≤ 5 s 单线程（§13.1）。

### 8.9 圆角、倒角、抽壳与偏移

| 操作 | 方法 | 失败模式（必须诊断） |
| --- | --- | --- |
| 圆角 | rolling-ball 近似；等/变半径；边链 | 半径超过局部曲率半径、相邻圆角干涉、尖角 |
| 倒角 | 距离/角度；与圆角共用边处理框架 | 面退化、相邻倒角冲突 |
| 抽壳 | 面偏移 + SSI + 修剪 + 封口 | 壁厚大于局部厚度、自交、封口失败 |
| 偏移 | 曲面偏移（NURBS 偏移不可精确表示 → 近似 + 误差报告）；实体偏移 | 偏移距离超过曲率半径 |

统一 blend 引擎：**边链 → 局部面片 → SSI → 缝合**；任何失败不得产生半成品，返回诊断与可回退的输入快照。

### 8.10 网格化（tessellation）

B-rep 与网格世界的**唯一契约**（ADR-017）：

| 契约项 | 要求 |
| --- | --- |
| 水密 | 相邻面共享边的采样点逐点相同；闭合实体无边界边 |
| 偏差 | 弦高 ≤ 用户容差；角度偏差 ≤ 用户角度容差 |
| 无翻转 | 所有三角形法向与面 sense 一致；无零面积 |
| 覆盖 | 面参数域被完整三角化；内环挖空正确 |
| 可复现 | 同输入同参数 → 同网格（D0） |

```
tessellate(Solid, deflection):
  边采样（所有边，容差驱动）→ 面参数域 CDT(§8.2，边界点复用)
  → 3D 抬高 + 中点/曲率偏差估计 → 超差则细分 → 焊接 → 定向
```

输出 `TriangleSoup`（共享顶点）或直接 `HalfedgeMesh`；开壳/多壳/退化边（极点）显式处理。目标吞吐 ≥ 1M 三角形/s 单线程（中等复杂度曲面）。

### 8.11 质量属性

- 体积/面积/质心/惯性矩：散度定理 + 裁剪面 Gauss–Legendre 积分，误差界显式给出；
- 与网格化结果交叉验证（容差内一致），不一致则报 `Errc::numeric_failure`；
- 零面积/零体积面贡献为 0，但计入诊断计数。

### 8.12 修复与缝合（heal）

- **缝隙桥接**：找容差带相交的 coedge 对，合并顶点或插入桥接边，记录 `gap` 量值；
- **容差重分配**：显式、单调、有上限；禁止「放大到能过」；
- **小边/小面移除**：折叠 + 局部重拟合，标记 `lossy` 并给误差上界；
- **重拟合曲面**：近似（最小二乘 NURBS），误差上界写入报告；
- **拒绝策略**：`gap > 阈值` 或修复会破坏 Euler 不变量时，返回诊断而非静默修改。

### 8.13 网格处理

| 能力 | 算法 | 关键约束 |
| --- | --- | --- |
| Halfedge 构建 | 索引三角网 → 半边，O(n) | 非流形边（>2 面）与边界必须显式区分；构建失败返回 `not_manifold` 并给位置 |
| 网格校验 | 流形性、定向一致、水密性、欧拉特征、自交 | 作为公共工具，所有布尔/简化输出前强制执行 |
| 简化 | QEM（二次误差度量）+ 边界/法向约束 | 保持流形与特征边；误差报告实际 Hausdorff 上界 |
| 重网格 | 各向同性（Botsch–Kobbelt 式分裂/坍缩/翻转/松弛） | 目标边长按尺度归一化；禁止翻转 |
| 平滑 | Laplacian / Taubin | 体积损失可量化 |
| 定向 | 基于有向体积符号的连通分量翻转 | 多分量分别定向 |

### 8.14 网格布尔（网格侧最高风险项）

```
A, B (closed 2-manifold 三角网)
  → 显式焊接/吸附(可选, 记录被合并顶点数)
  → 双 BVH 候选对(精确包围盒剪枝)
  → 精确符号的三角形-三角形相交(交点坐标近似, 拓扑关系精确)
  → 共面重叠: 投影到平面做 2D 精确布尔, 结果回投
  → 交点注册 → 边细分 → 重三角化(局部约束 Delaunay, 保形)
  → 区域分类: 射线奇偶 + 精确符号(或 winding number)
  → 保留/丢弃 → 缝合 → 一致性修复(去重、合并共线、闭合边界)
  → 不变量校验: manifold / watertight / orientation / 无自交
```

保证声明（写入文档与测试）：
- 输入为 **closed 2-manifold 且无自交** → 输出为 closed 2-manifold，拓扑一致；坐标误差 ≤ 输入尺度 × 1e-12；
- 输入不满足条件 → 返回 `Result` 诊断（含失败元素 id），**绝不输出「看起来像但拓扑非法」的网格**；
- 不承诺与商业内核的数值级一致（这是浮点几何的本质限制）。

### 8.15 SDF / Marching Cubes / 体素化

- 窄带 SDF：BVH 最近距离 + 符号（winding number 或射线奇偶，精确符号）；
- Marching Cubes 用标准表 + 顶点焊接；等值面在网格边界处需外扩 1 层；
- 体素化（表面/实体）：保守面片法 + 内部填充，输出可逆（误差报告）。

### 8.16 其他查询级算法

最近点对、Minkowski 和（2D 凸多边形）、线段/三角形相交、点集直径（旋转卡壳 2D）、多边形面积/周长/质心（精确求和：Kahan/补偿求和）。均要求「结果 + 误差报告」。

---

## 9. 并发与执行模型

### 9.1 线程安全契约

| 对象 | 契约 |
| --- | --- |
| 只读索引/网格/多边形集 | 多线程并发读安全（无内部可变缓存） |
| 只读 B-rep `Solid` | 并发读安全；建模算子不修改输入，返回新 `Solid`（事务语义，§6.8） |
| 算法调用 | 无全局可变状态；所有临时内存来自调用方 arena 或线程局部 |
| 用户回调 | 由调用方保证线程安全；库不串行化回调 |
| 数值环境 | 库不修改舍入模式/FP 环境标志 |

### 9.2 执行策略

```cpp
struct Serial {};
struct Parallel { unsigned threads; size_t grain; };
template <class Exec, class F>
void parallel_for(Exec, size_t n, size_t grain, F&& f);
```

- 默认串行；并行必须显式传入执行策略，不做「偷偷开线程」；
- 并行分区基于确定性划分（Morton/索引区间），确保 D0 算法与线程数无关；
- B-rep：面配对与 SSI 可按面并行，**拼接/修复/校验串行**（D1）；输出顺序按面 id 归一化；
- 取消：`CancelToken` 可选传入长算法（B-rep 布尔、网格布尔、重网格、网格化），在检查点轮询。

### 9.3 确定性等级

| 等级 | 含义 | 适用 |
| --- | --- | --- |
| D0 | 逐位一致（跨平台/跨线程数） | 谓词、索引构建、凸包、Delaunay 结构 |
| D1 | 拓扑一致，坐标允许 ulp 级差异 | 并行浮点归约、部分网格算法 |
| D2 | 仅几何等价（容差内） | 迭代式几何优化（重网格、平滑） |

每个算法头文件显式标注等级，测试按等级选取断言方式。

---

## 10. 扩展机制与互操作

### 10.1 concepts 扩展点

| 扩展点 | concept | 用户可替换 | 示例 |
| --- | --- | --- | --- |
| 标量 | `Scalar` | ✔（内置 float/double/long double 经 traits；自定义标量特化 `scalar_traits`） | 嵌入式用 `float` 省内存 |
| 精确后端 | `ExactReal` | ✔（有理数/大浮点） | 离线高精度校验 |
| 空间索引 | `SpatialIndex` | ✔ | 用户已有加速结构接入查询算法 |
| 分配器 | `std::pmr::memory_resource` | ✔ | 游戏引擎 arena |
| 属性系统 | `Attributes` | ✔（自定义 `type_id` + POD） | 材质、颜色、工艺标签 |
| I/O | `MeshReader`/`MeshWriter` | ✔ | 私有格式 |

### 10.2 互操作

- **Eigen**：`Eigen::Map` 零拷贝互转（`adapters/eigen.hpp`，仅显式 include）；
- **Python**：pybind11 + `buffer protocol` 零拷贝接 NumPy；批量接口一次跨界；
- **C API**：不透明句柄（网格/多边形/`Solid`）+ `geom_status` 返回码 + 显式释放函数；唯一承诺 ABI 稳定；
- **STEP/IGES**：内核**不内置**解析（§1.3）；由适配器或应用层转换为 B-rep API，预留 `BrepImporter` 概念；
- **WASM**：`-fno-exceptions` 变体 + 单线程执行策略；
- **序列化**：见 §6.7/§6.8；glTF 仅作为可选 I/O，不进入核心依赖。

---

## 11. 构建、依赖与打包

### 11.1 CMake 目标

| target | 类型 | 内容 | 依赖 |
| --- | --- | --- | --- |
| `geom::core` | INTERFACE | L0–L2 全部头文件 | 无 |
| `geom::brep` | STATIC | L5 B-rep 内核（SSI / 布尔 / 圆角 / 网格化 / 质量属性 / healing） | core, index, cgeom |
| `geom::algo` | STATIC/SHARED | L5 网格算法（网格布尔 / SDF / 体素化 / 偏移） | core, index, cgeom |
| `geom::io` | STATIC | L6 I/O（含输出前网格/多边形校验） | core, cgeom |
| `geom::capi` | SHARED | C ABI | brep, algo, io |
| `geom::tests` / `geom::bench` / `geom::fuzz` | 内部 | 测试基础设施 | 全部 |

### 11.2 编译选项矩阵

| 项 | 默认 | 说明 |
| --- | --- | --- |
| C++ 标准 | C++20 | concepts/ranges/`std::span`；C++23 `expected` 可用时自动启用 |
| 浮点 | `-O2`；**禁止 `-ffast-math`**；`-ffp-contract=off`（谓词 TU） | 误开 fast-math 会破坏误差界，CI 用编译检查拦截 |
| 异常 | 开（可关） | `GEOM_NO_EXCEPTIONS` 时错误仅走 `Result`/abort |
| RTTI | 关（默认） | 无运行时类型识别需求 |
| 并行 | 关（默认） | 可选 TBB/自建线程池后端 |
| 精确后端 | 关（默认） | 可选 GMP/MPFR/boost::multiprecision |
| 警告 | `-Wall -Wextra -Wconversion`（/W4） | 库自身零警告 |

### 11.3 依赖策略

| 依赖 | 角色 | 强制 | 隔离方式 |
| --- | --- | --- | --- |
| C++ 标准库 | 基础 | 是 | — |
| GMP / MPFR | 精确数域 | 否 | `GEOM_USE_GMP` + `numeric/exact_gmp.hpp` |
| TBB / 线程池 | 并行 | 否 | `GEOM_USE_TBB` + 执行策略适配 |
| Eigen | 互操作 | 否 | 仅 `adapters/eigen.hpp` |
| zlib/zstd | 压缩序列化 | 否 | 仅 `io` 层可选 |
| CGAL / Manifold / Clipper2 / OCCT | **测试期差分参考** | 否 | 仅 `tests/oracle/`，不进入安装产物 |

### 11.4 分发

- CMake 安装导出：`geomConfig.cmake` + `geomTargets.cmake` + `pkg-config`；
- vcpkg / Conan 配方；`FetchContent` 支持；
- 语义化版本；`CHANGELOG.md` 记录 API/行为变更；
- CI 矩阵：{MSVC 19.3x, GCC 12+, Clang 15+} × {Windows, Linux, macOS} × {Debug, Release, RelWithDebInfo} × {exceptions on/off} × {sanitizers}。

---

## 12. 质量保障体系

测试策略的完整定义（覆盖目标、测试设计方法、工具链、缺陷流程、反模式）见 **[testing-strategy.md](../quality/testing-strategy.md)**；本章只固化架构级约束。

### 12.1 测试分层与覆盖门槛

| 层级 | 内容 | 运行时机 | 时间预算 |
| --- | --- | --- | --- |
| L1 单元 | 单函数/单类型：等价类、边界、错误路径 | 每次提交 | < 30 s |
| L2 组件 | 模块内组合（如 hull + Delaunay 计数一致） | 每次提交 | < 60 s |
| L3 性质/变形/差分 | 随机性质、变换不变性、与 oracle 对比 | PR（快速档）/ 夜间（全量） | 5 min / 1 h |
| L4 系统/集成 | C API、序列化往返、语言绑定 | PR | < 10 min |
| L5 非功能 | 性能、TSan、fuzz、变异 | 夜间 / 每周 | 按任务 |

| 门槛 | core / numeric / primitives / curves / surfaces | cgeom / index | brep（M2+） |
| --- | --- | --- | --- |
| 行覆盖 | ≥ 95% | ≥ 90% | ≥ 85% |
| 分支覆盖 | ≥ 85% | ≥ 70% | ≥ 80% |
| diff 覆盖（本 PR 改动） | 100% | 100% | 100% |
| 变异分数 | ≥ 80% | ≥ 70% | ≥ 60% |

规则：
- diff 覆盖 100% 为硬门禁；未覆盖分支必须标注 `// GEOM_COVERAGE_EXEMPT: <理由>` 并评审；
- 每个 `Errc` 至少一个触发用例；每条公开 API 至少有正例、边界、失败路径三类用例；
- 谓词的过滤路径与精确路径都必须被测试执行（计数器证明，而非口头声明）；
- 覆盖率以 Linux + GCC + gcovr 结果为准，避免工具链差异；Windows 不做覆盖率统计。
- 分支覆盖按 gcov 语义计数（每个 `&&`/`||` 操作数与防御守卫各算一个分支）；M1b 实测
  全局 77.6%、cgeom 73.1%，门禁据此校准为全局 75%、core 85%、cgeom 70%（详见
  `docs/quality/testing-strategy.md` §2）。

### 12.2 测试五件套

| 类型 | 工具 | 覆盖目标 | 判据 |
| --- | --- | --- | --- |
| 单元测试 | doctest（低编译成本） | 每个公开函数 + 每个失败路径 | 覆盖门槛见 §12.1；断言精确符号/精确值，禁止 `>= 0` 式弱断言 |
| 性质测试 | rapidcheck（固定种子） | 代数律、不变量、守恒量 | 随机 10k 例无违反；失败打印种子与最小输入 |
| 变形测试 | 自研 generator | 平移/旋转/缩放/输入置换/细分不变性 | 几何结果在容差内一致；1e-200…1e200 尺度不改变符号 |
| 差分测试 | 整数 oracle + CGAL / Manifold / Clipper2 / OCCT（仅测试） | 谓词 1e6 例；小规模（n ≤ 200）几何对比；B-rep 对 OCCT | 符号 100% 一致；拓扑同构 + Hausdorff ≤ 阈值；体积/面积一致 |
| 模糊测试 | libFuzzer（结构感知 mutator） | 解析器与算法入口 | 夜间 1 h；里程碑出口要求 24 h 无崩溃/超时/UB |

### 12.3 不变量校验器（公共工具）

| 模块 | 不变量 |
| --- | --- |
| 多边形集 | 环闭合、无自交（可选）、洞在壳内、绕序一致、面积 ≥ 0 |
| 三角网 | 索引合法、无退化面、定向一致、流形性、水密性、欧拉特征 |
| arrangement | 每条边恰有两个面/边界、面绕序与 winding 一致 |
| BVH | 节点包围盒包含子节点、叶节点与图元一一对应、无空内部节点 |
| 布尔输出 | 输入并集/差集的面积（体积）守恒在容差内 |
| B-rep 拓扑 | 欧拉-庞加莱 `V − E + F − (L − F) = 2(S − H)`；每条边恰 2 个 coedge；外环唯一；壳定向一致；无悬垂 handle |
| B-rep 容差链 | 顶点/边/面容差满足 §4.6 的 5 条不变量；`Solid::tol` ≥ 子实体最大容差 |
| pcurve | 每条 coedge 的 pcurve 与 3D 边在 face 容差内逐点一致（采样校验） |
| 网格化 | 水密（无边界边）、无翻转、弦高/角度偏差 ≤ 参数、共享边采样逐点一致 |
| SSI 输出 | 每段残差 ≤ tol；段间端点一致；交线位于两曲面上（采样校验） |

每个校验器自身必须有负向测试：喂入故意损坏的结构，断言返回具体 `Errc` 而非崩溃或静默通过。

### 12.4 退化语料库

- `tests/corpus/degenerate/`：共线、重合边、共面、重复点、零面积、反向绕序、近退化（1 ulp 级扰动）、超大尺度比（1e9:1），每类 ≥ 20 例；
- `tests/corpus/cad/`：相切接触、重合面、共边、退化边（极点）、细长面、容差链装配、多壳/内腔、自交输入，每类 ≥ 20 例；
- 每个模糊测试崩溃都固化为语料 + 单测；语料带 `.sha256` 与来源说明，删除需评审；
- 语料序列化格式与 §6.7 一致，便于跨平台复现。

### 12.5 静态检查与质量门禁

- clang-tidy（bugprone/performance/modernize 子集，`-Werror`）；
- include 依赖检查（§3.3 矩阵）；
- ABI 检查：仅对 `geom::capi` 使用 abi-compliance-checker；
- 分配计数测试：查询路径 0 分配；
- 覆盖率门禁：CI 运行 gcovr，模块阈值见 §12.1，未达标阻塞合并；
- 变异测试：每周夜间，核心模块分数低于门槛时生成缺口清单并在下一迭代清零；
- 并发：TSan 夜间任务覆盖只读索引/网格的并发查询与线程数扫描；
- 确定性检查：D0 算法在多平台输出哈希一致（Windows/Linux/macOS × MSVC/GCC/Clang）；
- 模糊测试：夜间 1 h，里程碑出口要求 24 h 无崩溃；
- 性能门禁：金标准基准回归 > 5% 阻塞合并（见 §13）。

### 12.6 里程碑测试出口

| 里程碑 | 测试出口（与 testing-strategy §14 一致） |
| --- | --- |
| M0 | 谓词整数 oracle 1e6 通过；退化枚举全过；core/numeric 覆盖达标；跨平台符号哈希一致 |
| M1 | Delaunay 空圆性质 1e3 用例；BVH 差分 1e5 射线；覆盖达标；夜间 fuzz 24 h 无崩溃 |
| M2 | B-rep 不变量负向用例齐全；网格化水密率 100%；OCCT 差分通过 |
| M3 | 网格布尔体积守恒；SSI 语料通过率 100%（失败可诊断）；接入 OSS-Fuzz |
| M4 | 性能预算达标；变异分数达标；三平台矩阵全绿 |

---

## 13. 性能工程

### 13.1 初始预算（待 M1 基线校准，校准后作为 CI 门禁）

| 指标 | 目标值 | 条件 |
| --- | --- | --- |
| `sign_orient2d` 过滤路径 | ≤ 2.5 ns/call | x86-64, -O2, 无 AVX 特化 |
| `sign_orient3d` 过滤路径 | ≤ 5 ns/call | 同上 |
| 1M 点 2D Delaunay | ≤ 1.2 s 单线程；≤ 250 ms（16 线程） | 均匀分布 |
| 1M 三角形 BVH 构建 | ≤ 350 ms 单线程 | SAH 分箱 |
| 射线查询吞吐 | ≥ 8 M ray/s/核 | 网格 1M 三角 |
| 10 万顶点多边形并 | ≤ 80 ms 单线程 | 含 arrangement + 分类 |
| NURBS 曲面求值（p=q=3） | ≤ 80 ns/call | 单 patch，双精度 |
| B-rep 网格化吞吐 | ≥ 1M tri/s 单线程 | 中等复杂度曲面 |
| SSI 单对面（横截） | ≤ 100 µs | 细分深度 ≤ 12 |
| 1 万面 B-rep 布尔 | ≤ 5 s 单线程 | 单一操作，含修复与校验 |
| 质量属性（1 万面） | ≤ 50 ms | 体积/质心/惯性矩 |
| 1M 三角网格内存（halfedge） | ≤ 120 MB | 含位置与拓扑 |
| 单模块头文件编译 | ≤ 0.5 s/TU | MSVC/GCC/Clang 取最差 |
| 聚合头编译 | ≤ 2 s/TU | 同上 |
| 静态库体积 | ≤ 8 MB | Release，无调试信息，仅 double 实例化 |

### 13.2 热点路径设计检查表（PR 评审强制）

- [ ] 无堆分配（析构函数/RAII 也不得分配）
- [ ] 无虚函数、无 `std::function`、无 `shared_ptr` 拷贝
- [ ] 无异常抛出；错误以返回值表达
- [ ] 无 `pow`/`sqrt` 滥用（优先平方距离/比较）
- [ ] 无隐式 `double`↔`float` 转换（`-Wconversion` 把关）
- [ ] 分支可预测性：退化判定前置，热路径分支少而稳定
- [ ] 数据布局：热字段相邻、避免 false sharing（并行分区按 64 B 对齐）
- [ ] 不修改全局 FP 环境

### 13.3 方法

- Google Benchmark + 固定随机种子 + 统计显著性（中位数 + 置信区间）；
- 性能门禁：3 个金标准基准回归 > 5% 阻塞合并；
- SIMD：`geom::simd` 最小封装（SSE2/AVX2/NEON）+ 标量 oracle 对拍（逐元素语义一致，非逐位）；
- 谓词保持标量（分支多），批量查询/包含测试可 4/8 宽向量化；
- 定期做「无过滤 vs 过滤」对比，防止过滤层本身成为开销。

---

## 14. 版本、兼容性与治理

| 项 | 策略 |
| --- | --- |
| 版本 | SemVer；`MAJOR` = 破坏性 API；`MINOR` = 兼容新增；`PATCH` = 修复 |
| API 兼容 | 源码兼容优先；破坏性变更走 `v2` 命名空间 + 保留期 ≥ 2 个 minor |
| ABI | 仅 `geom::capi` 承诺稳定；C++ 侧不承诺（显式文档化） |
| 废弃 | `[[deprecated]]` + CHANGELOG + 迁移指南；至少一个 minor 后才移除 |
| 实验 API | `geom::experimental`，无兼容承诺 |
| 决策记录 | `docs/adr/NNNN-title.md`，一页纸：背景/选项/决策/后果 |
| 评审 | API 变更需 RFC + 两名维护者批准；架构变更需更新本文档并版本+1 |

---

## 15. 风险登记册

| # | 风险 | 影响 | 概率 | 缓解 | 触发信号 |
| --- | --- | --- | --- | --- | --- |
| R1 | 布尔运算（B-rep 与网格）在退化输入上不可靠 | 高（核心卖点受损） | 高 | 严格不变量校验 + 可诊断失败 + 语料库持续扩充；必要时限制输入契约 | 差分测试失败率 > 0.1% |
| R2 | 精确算术性能悬崖（Tier D 比例升高） | 中 | 中 | 监控 `NumericStats.exact` 比例；优化过滤层；分级降级策略 | 精确率 > 0.1% |
| R3 | 模板/头文件导致的编译期爆炸 | 中 | 中 | 显式实例化、类型擦除门面、分层头文件、`extern template` | 聚合头 > 2 s |
| R4 | API 膨胀与命名不一致 | 中 | 中 | API 评审清单 + 面向任务的能力矩阵 + 删除非必要重载 | 公开符号数季度增长 > 20% |
| R5 | 用户误用容差导致「看起来对」的结果 | 中 | 高 | 显式容差 + 尺度归一化 + 报告实际偏差 + 文档示例 | 用户反馈「结果偏移」 |
| R6 | 平台浮点差异（FMA 收缩、x87、FTZ） | 高 | 中 | 编译选项矩阵 + 双误差界 + 跨平台确定性测试 | 哈希不一致 |
| R7 | 并行引入非确定性 | 中 | 中 | 确定性分区 + D0/D1 分级 + 线程数扫描测试 | 同输入不同结果 |
| R8 | 参考实现许可污染 | 高 | 低 | 仅测试期依赖、独立 `tests/oracle` 目标、产物扫描 | 产物含 GPL 符号 |
| R9 | 性能回归悄悄累积 | 中 | 高 | 基准 CI 门禁 + 仪表盘 + 季度性能审计 | 5% 回归告警 |
| R10 | 范围蔓延（特征树 / STEP / 装配约束） | 高 | 中 | §1.3 非目标 + 里程碑出口条件 + 架构评审 | 需求池出现参数化历史/装配求解 |
| R11 | SSI 在切向/重合输入上失败率超预期 | 高 | 高 | 退化分类 + CAD 语料库 + 诊断回退；必要时限制契约（要求横截相交） | SSI 诊断率 > 1% |
| R12 | 容差链不一致导致缝合失败或模型漂移 | 高 | 中 | 容差传播规则（§4.6）+ 每步校验 + 增长上限告警 | `tol_growth > 10×` |
| R13 | 圆角/抽壳失败率高，用户信任下降 | 中 | 高 | 明确失败模式 + 诊断 + 不做半成品 | 失败率 > 5% |
| R14 | 高精度算术在 SSI 上性能悬崖 | 中 | 中 | 区间过滤 + 仅认证阶段升精度 + 按面并行 | 单对面 > 1 ms |

---

## 16. 演进路线（架构演进，非项目排期）

### 16.1 能力矩阵（"全面"的可验收边界）

| 能力域 | 具体能力 | 落地里程碑 |
| --- | --- | --- |
| 谓词与图元 | 2D/3D 精确谓词、全套图元、仿射变换、包围盒 | M0 |
| 曲线曲面 | 解析曲面、NURBS/Bézier 求值/导数/反求/包围、弧长 | M0 |
| 2D 计算几何 | 凸包、Delaunay、CDT、简单多边形三角剖分、布尔、arrangement | M1 |
| 空间索引 | BVH/KD/R-tree、射线、最近邻、broad-phase 相交对、自交检测 | M1 |
| B-rep 基础 | 拓扑、Euler 算子、基本体、网格化、质量属性、不变量校验器 | M2 |
| 曲面求交 | SSI（横截/切向/重合/边界）、曲线-曲线、曲线-曲面 | M3 |
| B-rep 布尔 | 面配对、面分割、分类、缝合、修复、容差模型 | M3 |
| 高级建模 | 圆角/倒角/抽壳/偏移、扫掠类、healing | M4 |
| 网格算法 | 网格处理、网格布尔、SDF/Marching Cubes、体素化 | M4 |
| 互操作 | STL/OBJ/PLY/glTF、C API、pybind11、Eigen、WASM | M5 |
| 首版明确不做 | 参数化特征树、STEP/IGES 解析、装配约束求解、GPU 加速 | — |

### 16.2 里程碑

| 里程碑 | 内容 | 出口条件（必须全部满足） |
| --- | --- | --- |
| M0 地基 | L0/L1/L2、谓词、曲线曲面求值、Build/CI、测试骨架、语料库雏形 | 谓词与 CGAL 符号 100% 一致（1e6 随机 + 退化集）；NURBS 求值/导数对解析解校验通过；编译矩阵全绿；依赖检查生效 |
| M1 计算几何与索引 | hull、Delaunay、CDT、多边形三角剖分、多边形布尔、KD/BVH 完整（拆分为 M1a/M1b，见 §16.3） | QA-1/QA-2/QA-4 达标；性质与差分测试通过；性能预算校准并入库 |
| M2 B-rep 基础 | 拓扑、Euler 算子、基本体、网格化、质量属性、B-rep 校验器 | QA-10/QA-11 达标；网格化水密率 100%；欧拉-庞加莱校验全绿；容差链校验生效 |
| M3 求交与布尔 | SSI、面分割、B-rep 布尔、容差模型、healing 初版 | QA-8/QA-9/QA-12 达标；R11/R12 缓解措施验证；失败全部可诊断 |
| M4 高级建模与网格 | 圆角/倒角/抽壳/偏移、扫掠类、healing 完整、网格处理/网格布尔/SDF/体素化 | QA-5 达标；R13 失败率 ≤ 5% 且全部带诊断；网格与 B-rep 双侧回归通过 |
| M5 产品化 | I/O、C API、pybind、Eigen、WASM、性能打磨、1.0 | QA-3/QA-6/QA-7 达标；C ABI 检查通过；文档与示例完备 |

依赖关系：M1 依赖 M0 的谓词、曲线曲面与测试设施；M2 依赖 M1 的 arrangement/CDT/索引；M3 强依赖 M2 的拓扑、网格化与容差模型；M4 在 M3 之上；M5 依赖前五者的 API 稳定。

### 16.3 M0/M1 出口状态（2026-11-07）

**M0 地基：已满足。** 谓词以独立整数/Fraction oracle 对拍（1e6 随机 + 共线、1 ulp、
抵消、1e±200、次正规语料），NURBS 求值与导数对解析解校验，三平台 CI 矩阵已配置；
`orient2d` 过滤路径实测约 48–76 ns/call（不同批次；预算 200 ns），见
`docs/benchmarks/m1.md`。

**M1 按两段验收，避免把未实现项算作达标：**

- **M1a（代码与测试完成；覆盖率已实测）**：`hull2d`、`delaunay2`（一般位置
  前置）、`Bvh`（raycast/nearest）及其性质与暴力差分测试；安装包与 consumer 工程
  验证。对应实施计划 T11–T14。CI `coverage` job（Release 覆盖率构建）实测：全局
  行 96.6%、分支 77.6%；core 模块分支 86.7%、行覆盖 ≥95%。按 §12.1 校准后的门禁
  通过。
- **M1b（代码与测试完成；验收未达成）**：`cdt2`、`triangulate2`、
  `arrangement2`/`boolean2`、KD 树与 BVH 的代码及性质/差分测试已完成；QA-2
  一致性测试与 `geom_bench_m1b` 已入库（实现取舍见 ADR-020）。**验收状态
  （2026-11-08 实测，见 `docs/benchmarks/m1b.md`）**：QA-1 FAIL，仅为性能——
  winding 播种修复后面积相对误差 1.78e-13（门禁 1e-12）已通过；拆分阶段优化后
  并集 1.08 s（较 4.8–4.9 s 降低约 4.5×），仍超 80 ms 门禁约 13×；QA-4 FAIL
  （构建 3.2 s vs 350 ms，0.34 Mray/s vs 8）；QA-2 的 `tsan` job 已在
  ubuntu-24.04 CI 通过；覆盖率实测 cgeom 分支
  73.1%（门禁 70），行 ≥90。因此 **M1b 代码质量门禁通过，但 QA-1/QA-4 性能验收
  未达成**，不得据 §16.2 的 M1 行声称已达标。

依赖检查：M0/M1a 的公开头文件不依赖 doctest；`Result` 只用于可失败构造与算法
错误返回（`Affine::inverse`、`boolean2`、`split_segments`、`cdt2` 等）；谓词之外
无全局 epsilon。

---

## 17. 附录

### A. 目录结构

```
ai_geometry/
├── CMakeLists.txt
├── cmake/                     # 安装导出、presets、编译器检查
├── include/geom/
│   ├── core/                  # config, result, error, handle, arena, span, meta
│   ├── numeric/               # traits, interval/affine, predicates2/3, exact, sos
│   ├── primitives/            # point, vector, box, segment, ray, plane, triangle, ...
│   ├── curves/                # line, arc, ellipse, nurbs_curve, bezier_curve
│   ├── surfaces/              # plane, cylinder, cone, sphere, torus, nurbs_surface, bezier_patch
│   ├── transform/             # affine, frame, quaternion
│   ├── index/                 # bvh, kdtree, rtree, grid
│   ├── cgeom/                 # hull, triangulation, polygon boolean, arrangement, mesh
│   ├── brep/                  # topology, euler, model, ssi, boolean, blend, tessellate, heal, massprops
│   ├── algo/                  # mesh_boolean, sdf, voxelize, offset
│   ├── query/                 # nearest, raycast, overlap
│   ├── io/                    # stl, obj, ply, gltf
│   ├── capi/                  # C ABI
│   ├── adapters/              # eigen, python
│   └── experimental/
├── src/                       # 编译单元（显式实例化、算法实现）
├── tests/
│   ├── unit/  property/  metamorphic/  oracle/  corpus/
│   ├── corpus/degenerate/     # 计算几何退化
│   └── corpus/cad/            # 相切/重合/容差链/多壳
├── bench/
├── fuzz/
├── tools/                     # check_layers, alloc_counter, corpus_minimizer
├── docs/
│   ├── architecture/          # 本文档
│   ├── adr/
│   └── user-guide/
└── examples/
```

### B. 术语表

| 术语 | 含义 |
| --- | --- |
| 过滤谓词 | 带静态误差界的快速浮点判定，符号可信则直接返回 |
| 自适应精确 | 按需扩展精度的展开式算术，只求符号 |
| SoS | Simulation of Simplicity，符号扰动打破退化平局 |
| 规则化布尔 | 输出为闭集，零面积/零体积部分被吸收 |
| arrangement | 由边、顶点、面构成的平面细分（DCEL 表示） |
| winding number | 点相对多边形的环绕数，决定内外 |
| halfedge | 半边，每个无向边拆为两条有向半边，O(1) 邻接遍历 |
| D0/D1/D2 | 确定性等级：逐位一致 / 拓扑一致 / 几何等价 |
| 水密射线求交 | 保证相邻三角形共享边处不漏检/不重复的浮点一致求交 |
| B-rep | 边界表示：用「面（曲面 + 环）/ 边 / 顶点」描述实体，CAD 内核的标准表示 |
| coedge | 半边（有向边）在 B-rep 中的称谓：携带所在环、朝向与 pcurve |
| pcurve | 曲线在某个面参数域中的 2D 表示，用于面上求交、裁剪与网格化 |
| 容忍建模 | 实体带显式容差的建模范式（Parasolid/ACIS 系）；容差是一等数据，不是全局 epsilon |
| SSI | Surface–Surface Intersection，曲面-曲面求交 |
| Bézier clipping | 用 Bézier 控制多边形/区间裁剪参数域的细分求交技术 |
| Euler 算子 | 保持欧拉-庞加莱不变量的原子拓扑修改（MEV/MEF/KEF…） |
| 欧拉-庞加莱公式 | `V − E + F − (L − F) = 2(S − H)`，B-rep 有效性校验的基准 |
| 网格化 / 弦高 | 用三角网逼近 B-rep；偏差由弦高（deflection）与角度容差控制 |
| 双树遍历 | 两个 BVH 同步下降求候选对（宽相位）的遍历技术 |

### C. 与主流库的定位对比

| 维度 | GeomCore | CGAL | Boost.Geometry | Clipper2 | Manifold | OpenCASCADE |
| --- | --- | --- | --- | --- | --- | --- |
| 定位 | 通用计算几何 + CAD 曲面内核 | 算法大全、精确几何 | 2D 泛型几何 | 2D 整数裁剪 | 网格布尔 | B-rep 内核 + 建模 |
| 强制依赖 | 无 | Boost（+GMP 可选） | Boost | 无 | 无（TBB 可选） | 自带；可选 TBB/FreeType |
| 精度模型 | 分层过滤 + 容忍 B-rep | 精确内核（EPICK/EPECK） | 重缩放 | 整数坐标 | 浮点 + 鲁棒构造 | 浮点 + 容忍建模 |
| B-rep 曲面内核 | 是（M2–M4） | 无完整 B-rep（有 Nef/多面体） | 否 | 否 | 否 | 是（标杆） |
| 3D 网格布尔 | 是（M4） | 有（接口重） | 否 | 否 | 是（标杆） | 有（网格化后） |
| 编译成本 | 低（分层 + 显式实例化） | 高 | 低–中 | 低 | 低 | 高 |
| 许可 | Apache-2.0（建议） | GPL/LGPL 混合 | BSL-1.0 | BSL-1.0 | Apache-2.0 | LGPL-2.1 + 例外 |

定位目标：**CGAL 常用计算几何的 80% + OCCT 级 B-rep 常用建模能力 + Manifold 级网格布尔 + Clipper2 级 2D 吞吐**，同时保持零强制依赖与低编译成本。

### D. ADR 索引

| ADR | 决策 | 状态 |
| --- | --- | --- |
| 001 | 单向 8 层分层，依赖由工具强制 | 提议 |
| 002 | 图元标量/维度泛型，精度内核仅在算法入口注入（否决 CGAL 式全参数化） | 提议 |
| 003 | 谓词分层：过滤 → 自适应 → 精确；符号精确、量值近似 | 提议 |
| 004 | 无全局 epsilon；容差显式并按尺度归一化 | 提议 |
| 005 | 拓扑用 32 位代际 handle；批量用 SoA 视图 | 提议 |
| 006 | 查询路径零分配、零虚函数、零异常 | 提议 |
| 007 | 核心零强制依赖；可选后端用编译开关隔离 | 提议 |
| 008 | 错误模型 `Result<T, Error>` + 可配置异常模式 | 提议 |
| 009 | 确定性分级 D0/D1/D2；并行必须显式选择且结果与线程数无关 | 提议 |
| 010 | 仅 C API 承诺 ABI 稳定 | 提议 |
| 011 | ~~首版不做 B-rep/NURBS~~ → 已被 ADR-013–019 取代（主战场调整） | 废弃 |
| 012 | 2D 布尔基于 arrangement + winding，否决 Vatti/Greiner–Hormann | 提议 |
| 013 | L5 更名为「建模内核与组合算法」：B-rep 与网格算法同层、互不依赖，经网格化契约交互 | 提议 |
| 014 | 容忍建模：每实体容差 + 精确谓词分层；容差显式传播并校验 | 提议 |
| 015 | 曲面求交：细分 + 区间认证 + 容差交付；不做全局精确有理 SSI | 提议 |
| 016 | NURBS 为统一后备表示，解析曲面保留类型身份与闭式快速路径 | 提议 |
| 017 | 网格化是 B-rep 与网格世界的唯一契约（水密 / 偏差 / 共享边采样） | 提议 |
| 018 | 内核只提供建模算子；特征树/历史、STEP/IGES、装配约束不做 | 提议 |
| 019 | B-rep 存储：扁平数组 + 代际 handle + arena 几何 + 属性稀疏表 | 提议 |
| 020 | M1b arrangement/三角剖分取舍：桥接+耳切、x 排序宽相位、`input_pairs` 自交判定 | 接受 |

### E. 需产品侧确认的决策（不阻塞架构）

1. **许可**：文档按 Apache-2.0 设计（便于商业集成）；若需 GPL 兼容或双许可，需在 M0 前确定。
2. **特征树边界**：内核只提供无历史的建模算子（ADR-018）。参数化历史若要做，放应用层还是内核新增历史层？直接决定 `Solid` 是否携带操作记录。
3. **STEP/IGES 适配器**：内核不解析（§1.3）。若首年需要交换能力，应单独立项（解析器 + 单位/容差映射），不并入内核工期。
4. **精确后端形态**：默认仅编译期可选；若线上需要运行时动态切换精度，需要额外设计 ABI 与性能隔离。
5. **容忍建模作为对外契约**：是否把「B-rep 容差语义」写进产品承诺（A6）。若坚持逐位精确，需重估 SSI 与布尔的技术路线（成本数量级上升）。

---

*本文档只描述架构与决策，不含实现代码；任何实现必须遵循本文件的分层与契约，偏离需通过 ADR 记录并经架构评审。*
