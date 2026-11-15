# ADR 0029: M3c-3a 球面帽的闭合纬度环布尔

- Status: Accepted
- Date: 2026-11-15
- Related: architecture §8.12、§16.3；M3c-3a plan Task 1–Task 6；ADR-0025、ADR-0027、ADR-0028

## 背景

M3c-2 结束时 `boolean_solid_curved` 已支持平面、圆柱面与圆环面操作数，球仍是显式
`unsupported`（ADR-0028）。球的困难与圆环不同，不是切平面限制，而是三件事叠在
一起：

1. `make_sphere` 与 `make_torus` 一样只保存一条面积为零的缝合线环（u 方向整圆在
   两极塌缩成点），切分器如果把环当作参数域，`curved_boolean_clip_pair` 会把每条
   切链都裁掉，因为切链不落在零面积环内。
2. 球的 v 图卡是 `[-π/2, π/2]`，极点行（合成矩形的上下边）在 3D 中塌缩成一个点。
   一个触及极点的 cap 片段若仍按“缝合线两端闭合”处理，会在极点附近留下一整行
   重复点或让环自交。
3. 两个纬度环之间的 band 片段在 uv 上是两条闭合环，但它们不是两个独立边界：
   同一条 band 的上下 rim 要在缝合列处拼成一个环，否则装配会看到两个零面积边界。

M3c-3a 只切“闭合且经过缝合线一点的纬度环”：球 ∩ 平面（球帽）、球 − 平面
（帽的补集）、球 ∩ 平板（两个纬度环之间的 band）。球/球 lens、锥、偏移 torus 与
一般修剪网格化都不在本里程碑内，先把这三类结果端到端固定下来，给球操作数一个
可验证的出口。

## 选项

1. **通用球面/锥面裁剪域。** 对球面片、球/球透镜、锥面片统一做 NURBS 轮廓 +
   多环域三角化 + 精确求交。优点是最终形态；缺点是球/球精确圆、跨链结合、孔洞
   域网格化各自都是独立里程碑，M3c-3a 无法在单一出口验证。
2. **只支持“闭合且过缝一点的纬度环”，复用合成矩形域。** 把 seam-only 面的
   解析参数矩形暴露给驱动，裁剪在矩形上做；极点 cap 用恰好一个精确极点 uv 点
   表示；band 的两个纬度环在缝合列拼接。优点是复用 M3b-2/M3c-2 已验证的切分、
   分类与装配；缺点是内部闭环、跨缝两次的闭环、开放切链仍拒绝，球/球 lens 与锥
   仍不能支持。
3. **用缝合线对闭合极点行。** 把触及极点的合成边当作普通缝合线两端（`u=0` 与
   `u=2π` 两个点）闭合。优点是改动最小；缺点是极点处的环会含一整行退化点，体积
   积分与装配都会看到零面积三角形，产出错误几何。评审已拒绝。
4. **继续把球标为 unsupported，等 M3c-3b 一次做完球/锥。** 缺点：球帽是最容易
   验证的球操作数，推迟它会让 M3c-3b 同时承担 cap 拓扑、球/球透镜与锥顶点三块
   风险，没有中间验证点。

## 决策

采用选项 2。

- **暴露合成参数域。** `CurvedBooleanFaceInfo` 增加 `synthetic_domain` 与
  `synthetic_boundary{u_lo, v_lo, u_hi, v_hi}`（`boolean_solid_curved.hpp:215-219`）；
  `curved_boolean_face_info` 在面的存储环是零面积缝合线时置位
  （`boolean_solid_curved.hpp:278-279`）。球的域是
  `[0, 2π] × [-π/2, π/2]`，圆环是 `[0, 2π]²`，锥是 `[0, 2π] × [0, r/tan θ]`。
  驱动的其余路径仍只读 `info.loops`，新分支全部由该标志门控，因此圆环/圆柱/平面
  行为不变。
- **在矩形上裁剪，不在缝合环上裁剪。** `curved_boolean_region_loops` 对
  `synthetic_domain` 面返回合成矩形的四角多边形，普通面仍返回真实 uv 环
  （`boolean_solid_curved.hpp:1149-1167`）。`curved_boolean_clip_pair` 用一个链
  采样在“本面区域”和“对面投影区域”上的成员判定裁剪切链：全在内部时保留整链，
  否则按段保留，并在边界穿越处调用 `curved_boolean_boundary_crossing` 求精确
  3D 交点再投影回本面 uv，使两个面在同一个 3D 顶点相接，而不是各自用中点近似
  （`boolean_solid_curved.hpp:1169-1260`）。
- **恰好一个精确极点。** 球图卡是唯一 v 域端点恰为 `-π/2`/`+π/2` 的 seam-only
  图卡，`trim_face_seam_only_pole_v` 据此识别极点行（`trim_face.hpp:891-907`）。
  一个 material cycle 拥有顶边或底边的一段时，环在 rim 末端与重复首点之间插入
  一个 `v` 取自域端点的极点 uv 点，再按 `rim → pole → rim` 闭合
  （`trim_face.hpp:1145-1200`）。`v` 直接来自域值，所以极点精确等于 `±π/2`，不会
  出现投影近似；两个不同的极点不会被闭合成一条经线。
- **band 在缝合列拼接。** 无极点且恰有两个环的片段就是两个纬度环之间的 band：
  两个环的端点分别是缝合线的两份参数副本（相差一整周），因此不能靠 seam 标志
  判断。拼接把两个环按顺序合并、丢掉各自的重复首点、在 `u=0` 缝合列处连成一个
  闭合环并重新归一化绕向；只有 `pole_chart && !touches_pole && rings.size() == 2`
  才执行，cap（一个环或触及极点）不会被误拼。两个零面积环会把同一区域写成两段
  边界，拼成一个环才是 band 的真正拓扑。
- **内部闭环仍然拒绝。** 切分器只连接终止在缝合线上的切链：一条完全位于域内部
  的闭环需要把两段弧拼成一个环（否则丢失一侧），一条跨缝合线两次的闭环需要连接
  多段 run；两者都在 arrangement 之前拒绝，调用方不会拿到部分切分。相关判据在
  `trim_face.hpp:1041-1076`，消息分别为
  `"split_face_curved: closed cut does not reach the seam"` 与
  `"split_face_curved: closed cut crosses the seam twice"`；落在参数域边界上的切链
  更早被拒绝为 `"split_face_curved: cut lies on the parameter domain boundary"`
  （`trim_face.hpp:1017-1027`）。这些早退是正确性保证，不得删除；球帽的纬度环
  恰好只有一个缝合线接触点（由两个参数副本表示），所以走的是支持路径。
- **球门提升。** `curved_boolean_surface_gate` 不再拒绝 `Sphere3`；入口消息
  `"boolean_solid_curved: sphere operands need cap topology; booleans are deferred to M3c-3"`
  随 M3c-2 的球缺口说明一起失效。门的范围仍是“闭合纬度环”：开放切链、内部闭环、
  跨缝两次的闭环继续返回 `unsupported`，球/球配对不在本里程碑。
- **体积与容差。** 体积仍走 `detail::curved_solid_volume`；球面在 M3c-2 已实现按
  极点分带的解析积分，本里程碑直接复用。全程容差参数 `1e-9`，QA-15 只要求体积
  相对误差 `<= 1e-6`；极点、缝合线与域边界都用精确值或已有恒等判据，没有新增
  全局 epsilon。

## 后果

- **QA-15（实测）。** `tests/unit/brep_boolean_curved_test.cpp` 的
  `QA-15: sphere cap booleans are closed and deterministic` 三行：
  `sphere ∩ box z ≥ 0.5` 体积 `0.6544984694978736`、
  `sphere − box` 体积 `3.534291735288517`、
  `sphere ∩ |z| <= 0.8` 体积 `3.9542179533184135`。每行跑两次，均满足
  `validate_topology().ok()`、恰好一个 shell、体积相对误差 `<= 1e-6`，第二次
  运行体积与第一次逐位相同。cap/complement 的球面 hat 与平面 cap 面可三角化，
  分类网格 `curved_boolean_mesh_closure(...).closed()` 通过；slab 的 band 面
  轮廓由两条跨缝合线的圆组成（8194 个 uv 点），`curved_boolean_mesh_face` 的
  `triangulate_polygon2` 拒绝该轮廓（`"boolean_solid_curved: a face outline could
  not be triangulated"`），因此该行改用结构判据
  `geom::detail::heal_detail::shell_is_closed_by_uses`：shell 中每条边恰被两个
  coedge 使用。单用例 83 条断言通过。
- **全套测试（实测）。** Debug 与 Release 各 940 个用例、659,153 条断言全绿；
  MSVC `/W4 /WX /utf-8` 零警告（GEOM_WERROR=ON）。
- **基准（`geom_bench_m3c3`，Release，100 次/操作）。** 三个操作的布尔调用中位数
  21.1 / 22.1 / 31.2 ms，体积积分中位数 232.0 / 232.0 / 214.8 ms；100 次体积
  逐位相同，checksum 与退出码 0 见 `docs/benchmarks/m3c3a.md`。体积积分仍是主要
  成本（约 7–11 倍），与 M3c-2 的结论一致。
- **限制清单（M3c-3b/3c）。** 球/球 lens（`d = 1.5`）需要跨链结合两条缝接触与
  精确球/球圆；锥操作数需要顶点投影与闭合/开放切链拓扑；偏移/倾斜 torus 切平面
  与精确 torus/cylinder 求交未做；布尔结果的一般修剪网格化（孔洞、非直边边界、
  跨缝轮廓）未做；NURBS 操作数未做。内部闭环、跨缝两次闭环与开放切链继续显式
  拒绝。任何没有真正切到球面的操作数都会被拒绝并返回错误，包括“盒完全在球内”这类
  合法交集——不会把未经修剪的整球或整盒当作结果返回。厚度低于约 `1e-2` 的球冠在
  当前弧段拼接与牛顿迭代下会失败（`1e-3` 起出现 arc run 不一致或收敛失败），可用
  下限以 QA-15 的 0.2 厚度语料为准。本 ADR 只关闭球面闭合纬度环，不得据此声称
  M3c-3 完成。
- **复用。** 合成矩形暴露、矩形裁剪与合成边剥离对锥可直接复用；极点处理只适用
  于 v 域端点为 `±π/2` 的球图卡；band 的“两环拼一环”判据以 `rings.size() == 2`
  为界，后续若出现三环及以上的无极点片段必须重新设计而不是扩展该分支。

## 重新评估的触发条件

M3c-3b 加入球/球 lens、锥操作数或偏移/倾斜 torus 切平面时；M3c-3c 加入布尔
结果的一般修剪网格化或 NURBS 操作数时；内部闭环/跨缝两次闭环需要支持时；极点
行与缝合线接触的判据需要放宽容差时；band 拼接遇到两个以上环或非纬度环时。
