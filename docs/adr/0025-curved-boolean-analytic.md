# ADR 0025: M3b-2 的解析曲面布尔

- Status: Accepted
- Date: 2026-11-12
- Related: architecture §8.8, §16.2, §16.3; M3b-2 plan Task 1–Task 9; ADR-0004、ADR-0021、ADR-0022、ADR-0023、ADR-0024

## 背景

M3a 交付了平面面实体的布尔（ADR-0023），M3b-1 交付了未修剪参数补丁上的
SSI（ADR-0024）。M3b-2 要把两者接起来：把 SSI 交线变成面上的修剪曲线，
把两个解析曲面实体的面切成片段，按正则化成员关系保留/丢弃片段，再焊接成
一个通过 `validate_topology` 的 `Brep`。这条管线要处理周期缝、采样边界、
分类容差与装配焊接四类数值问题，而 NURBS 操作数、共面接触、翻转曲面片段等
子问题在 M3b 阶段没有可靠的判据。

## 选项

1. **全解析曲面 + NURBS + 共面合并一次做完。** 一次到位，但把 NURBS 有理
   SSI（M3b-1 已阻塞）、共面参数域合并与曲面片段翻转同时拉进来，风险与
   验证面都过大。
2. **只做 Plane3/Cylinder3 操作数，其余显式诊断。** 管线的拓扑、容差、
   焊接与分类机器与最终形态一致，但把未验证的配置挡在入口之外。
3. **网格布尔。** 丢弃 B-rep 身份、容差语义与精确拓扑，与 M2/M3a 的方向冲突。
4. **参数域精确有理运算。** 超出 M3b 范围。

## 决策

M3b-2 采用选项 2，并明确记录以下契约：

- **解析作用域。** 操作数面只允许 `Plane3` 与 `Cylinder3`。`Sphere3`、`Torus3`
  与锥面在入口处返回 `Errc::unsupported`（`"boolean_solid_curved: sphere operands
  are deferred to M3c"`、`"... torus operands are deferred to M3c"`、
  `"boolean_solid_curved: cone operands are not supported (apex-touching loops are
  deferred)"`）；NURBS 面返回 `"boolean_solid_curved: NURBS operands are not
  supported"`。球与圆环的 uv 域一个周期折叠、一个双周期，当前分类网格无法
  自洽三角化，因此是"诚实延后"而不是"部分支持"。
- **闭合边采样 + 精确圆愈合。** 闭合圆边在参数域按弦高采样（默认
  `TessellationOptions::deflection = 1e-3` 时单位圆 71 段），面上的环因此是多段
  折线；但交线落到的圆弧在 `split_face_curved`/`split_solid_heal_loop` 中被识别为
  `CurvedArcSpan`，片段环保留精确圆弧而不是弦，体积误差由圆面积给出而不是
  由采样密度给出。
- **无共面接触。** 两个面片的平面部分在面积上重合时返回
  `"boolean_solid_curved: coplanar face contact is not supported yet"`；M3a 的
  平面布尔仍可处理纯平面操作数。相切同样诊断为
  `"boolean_solid_curved: tangent face contact is not supported yet"`，不猜测
  空/满的分类。
- **翻转曲面片段不支持。** 差集/对称差中若某个曲面片段必须以翻转朝向外，
  在切面之前就返回 `"boolean_solid_curved: flipped curved fragments are not
  supported yet"`：曲面片的环反向需要重排圆弧跨度，当前没有验证过。
- **结果只允许单连通分量。** 不相交并集、产生腔体的差集返回
  `"boolean_solid_curved: result is not a single connected component"`；M3a 的
  多实体/腔体支持不变。
- **周期缝解缠。** 交线的参数坐标经 `face_uv_project` 落在 (-π, π]，而面的 uv
  窗口是 [u₀, u₀+2π]。驱动在 `curved_boolean_clip_pair` 里把每个裁剪片段重新
  wrap 进面窗口，并在 `u = 0`/`u = 2π` 处切开；闭合一整圈的交线旋转到缝边，
  两端精确落在两条缝边上，装配后成为一条被使用两次的闭合边。
- **分类容差与接合策略。** 端点接合的容差是
  `junction_tol = max(64·tolerance, 1e-5·uv 对角线)`：SSI 追踪在采样步长的一小
  部分内停止（单位圆柱上实测 2.5e-5 rad ≈ 4e-6 的 uv 对角线），而"真正停在
  面内"的负例（约 1e-2）比它大 150 倍以上，因此不会被吸收。修复时：
  同一角点上的两条开链用各自端段所在直线的交点重建共享角点（归一化后
  接受 `|along| ≤ 2·junction_tol` 的 uv 位移），而不是焊到任一原始端点；
  圆弧边界只在能落到精确圆点时接受（径向投影回圆），弦投影会被拒绝，
  否则 `map_ring` 会报"片段顶点离开了它的环"。圆弧台阶阈值 `0.5 rad`
  区分"圆弧段"与"边界弦"，使愈合不会把弦误当成圆弧。
- **修剪曲面面不参与网格化与质量属性。** `tessellate` 要求圆柱面恰好有两条
  闭合圆边（cap circle），修剪后不满足；`mass_properties` 对非"未修剪整矩形"
  的曲面面返回 `"mass_properties: curved face is not an untrimmed full rectangle"`。
  体积只能通过 `detail::curved_solid_volume`（参数域面积分 + 散度定理）得到，
  这是 M3b-2 的正式体积入口；修剪后网格化与质量属性留给 M3c。
- **分类网格水密。** 分类网格优先用 `tessellate()` 的网格（完整基本体）；
  修剪操作数回退到"uv 三角化 + 按容差焊接 + 边界共形扇 + 闭合校验"。回退网格
  不是水密时返回 `numeric_failure "boolean_solid_curved: classification mesh is
  not watertight"`，绝不拿开口网格做内外分类。

## 后果

- **功能验收（12/12 QA 语料）。** 六个正例（偏移帽盒并/交/差、缝上盒并、
  x ≥ 0.5 平板交、贯穿双帽盒交）全部返回通过 `validate_topology` 的实体，体积
  相对解析值 ≤ 1e-6；六个反例（球、圆环、锥、NURBS、相切、共面）全部返回
  `Errc::unsupported` 且消息与上表逐字一致；每个输入跑两次，体积或
  (code, message) 完全相同。
- **性能（`geom_bench_m3b2`，5 次，Release，i7-8550U）。** 并集中位数
  156.9 ms（8 面）、交 156.3 ms（7 面）、差 152.2 ms（8 面），目标是中位数
  ≤ 2000 ms/操作，**PASS**；分类网格为默认 1e-3 弦高（整圆 71 段，142 顶点、
  280 三角形），驱动没有弦高旋钮，因此不虚构更细的采样。
- 曲面布尔的能力边界是显式的：球/圆环/锥/NURBS、共面接触、翻转曲面片段、
  多连通结果全部是带稳定消息的 `unsupported`，不会静默给出错误实体。
- 修剪曲面的体积可用但网格化/质量属性不可用，调用方必须接受这个非对称性，
  直到 M3c 补上修剪网格化。
- 缝感知的裁剪与接合只覆盖"平面-圆柱"这一对；纯直线交线对走普通裁剪路径
  （没有精确圆弧需要落点），行为与 M3a 的平面裁剪一致。
- `distinct_cuts` 在端点组内饱和为 2，三链共点的组按两链处理；当前语料没有
  触发，作为已知限制记录。

## 重新评估的触发条件

M3c 的 NURBS 操作数布尔、球/圆环的 uv 折叠处理、修剪曲面网格化与质量属性、
曲面片段翻转，或出现三链以上共点的真实用例时。
