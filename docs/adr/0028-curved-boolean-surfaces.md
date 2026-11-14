# ADR 0028: M3c-2 圆环操作数与缝合线锚定的曲面布尔

- Status: Accepted
- Date: 2026-11-14
- Related: architecture §16.3；M3c-2 plan Task 1–Task 8；ADR-0025、ADR-0027

## 背景

M3b-2 的 `boolean_solid_curved` 只接受 `Plane3`/`Cylinder3` 面：圆柱面有解析 cap
圆，平面有直线边界，切链与边界都能用直边/圆弧表示。M3c-2 要把 `Torus3` 纳入
操作数，但圆环面有两个 M3b-2 没有的性质：

1. `make_torus` 只生成一条缝合线环（外赤道整圆 + 管截面整圆），按环推导的 uv
   域面积为零；若直接沿用“域来自环”的切分路径，切链会落在退化域上。
2. torus ∩ plane 的交线一般不是圆弧（含轴平面之外），平面帽上的切圆虽然可以
   解析愈合，环面上的切链只能是 SSI 折线的弦逼近。

同一里程碑里球与锥各缺不同能力：球需要零面积环 + 极点 cap 环与极点扇形体积；
锥的切链只有内部闭环或终止在底圆的开放经线链，还需要顶点焊接。把三者一起做会
把 M3c-2 拆成三个独立的大任务，且没有中间的验证点。因此本 ADR 先固定“支持哪
类裁剪”，再固定操作数矩阵。

## 选项

1. **通用精确方案。** 精确 torus/plane、torus/cylinder 求交 + NURBS 轮廓 + 多环
   域三角化。优点是最完整；缺点是四次曲线求交、孔洞域网格化与精确愈合各自
   都是独立里程碑，M3c-2 无法在单一出口验证。
2. **缝合线锚定 + 合成参数域。** 只接受“闭合且经过缝合线一点的切链”，对
   seam-only 面（球/环/锥侧面）合成解析整矩形域，在合成域上跑既有 arrangement，
   再把合成边链剥掉。优点：复用 M3b-2 已验证的切面、分类、装配；缺点：开放
   切链与内部闭环不支持，因此锥不能在本里程碑支持。
3. **只加球/锥、不加环。** 缺点：球与锥都比环难（零面积环、极点 cap、顶点
   焊接、开放链），把最容易验证的操作数留到最后没有收益。
4. **近似/吸附对齐。** 把近轴、近垂直平面吸附到精确平面。缺点：产出与真实
   平面不符的几何，评审已拒绝。

## 决策

采用选项 2，操作数矩阵固定为：

- **圆环支持，附平面限制。** 切割平面必须垂直于环轴（归一化后 `|n·axis|` 在
  8 ulp 内等于 1）或包含环轴（`|n·axis|` 在 8 ulp 内为 0，且
  `|n·center − offset| <= tolerance`）。`make_torus` 只保证轴在
  `kAnalyticSurfaceUnitTolerance = 1e-12` 内为单位长且原样保存，单位轴也可能让
  归一化点积舍入到 `1 ± 1 ulp`，所以 `curved_boolean_torus_plane_cut_supported`
  与 `curved_boolean_pair_circles` 都先归一化轴，再用 8 ulp 带比较；带外一律
  在配对前拒绝，消息
  `"boolean_solid_curved: torus cut plane must be perpendicular to or contain the torus axis"`，
  不做几何吸附。由于点积在 1 处是驻点，该带也放进约 `6e-8 rad` 的真实倾斜；
  这类输入若进入精确切分路径会在 `trim_face` 明确失败
  `"split_face_curved: segment intersection is not representable"`，不会产出错误
  几何。偏移平行平面与倾斜平面留给 M3c-3。
- **合成参数边界。** seam-only 面在切分前合成解析参数矩形（环为 `[0, 2π]²`，
  球为 `[0, 2π] × [-π/2, π/2]`），在合成域上走既有 arrangement；随后剥掉
  落在合成矩形边上的链（`curved_boolean_chain_on_domain_side`，容差
  `max(64·tol, 1e-6·2π)`），闭合切链按 3D 端点恒等闭合。合成边不进入任何
  `CurvedFragment::loops`，装配只看到真实切线。
- **v 图卡包装。** 环的 `face_uv_project_torus` 返回 `(−π, π]`，而域是
  `[0, 2π]`；切链先在每个周期轴上展开，再包装进 `[base, base + 2π)` 并在该轴
  的缝合线处切分（`curved_boolean_window_chains`）。跨 `v = 0/2π` 的闭链因此
  不会留下 `2π` 跳变或开放端点。
- **只支持缝合线锚定的闭合切链。** 开放切链返回
  `unsupported "split_face_curved: open cuts on a seam-only face"`；不接触缝合线
  的内部闭环返回 `"closed cut crosses the seam twice"`。这条限制是锥不能在本
  里程碑支持的直接原因。
- **圆环切链是弦逼近。** 环面交线按 SSI 采样折线切分，没有精确圆弧愈合；
  残差由 `CurvedBooleanReport::max_cut_residual` 报告，QA-14 要求 `<= 1e-3`。
  平面帽上的切圆仍可解析愈合：垂直轴平面给出 `R ± sqrt(r² − d²)` 两圆，含轴
  平面给出两个管截面圆；这是 QA-14 体积能到 `1e-6` 相对的原因。
- **体积积分。** 修剪曲面面沿 uv 多边形三角化后按 7 点 Dunavant 在精确曲面上
  积分；环的采样冗余（每个 coedge 32 段）先用 `1e-12 · uv 尺度` 的容差简化，
  否则旋转后仅差几 ulp 的“共线”点会让耳切从线性退化到分钟级；若简化会把环
  压到少于 3 点则保留原采样，退化行为与简化前一致。球面走极点分带解析积分，
  供 M3c-3 使用，本里程碑不测球布尔。
- **球明确不支持。** 缺口是零面积缝合环的裁剪、极点绕行 cap 环、
  `split_solid_heal_fragment` 把绕行经线折进切割弧跨度、`球 ∩ 球` 的跨链结合；
  消息 `"boolean_solid_curved: sphere operands need cap topology; booleans are deferred to M3c-3"`，
  `face_pairs == 0`。
- **锥明确不支持。** 锥/平面切只产生内部闭环（被判 closed cut crosses the seam
  twice）或终止在底圆的开放经线链（被判 open cuts），且需要顶点焊接与顶点
  扇形体积，`curved_solid_volume(cone)` 目前也在顶点投影处失败；消息
  `"boolean_solid_curved: cone operands need closed-cut and apex topology; booleans are deferred to M3c-3"`，
  `face_pairs == fragments == kept == 0`。

## 后果

- **QA-14（实测）。** `tests/unit/curved_boolean_qa14_test.cpp` 三行：
  `torus ∩ z ≥ 0` 与 `torus ∩ x ≥ 0`（跨 v 缝）体积各
  `2.4674011002723395`，残差 `<= 1e-3`，两次运行体积逐位相同；
  `torus ⊂ 盒` 体积 `4.934802200544679`、`cuts == 0`。单用例 67 条断言通过；
  全量 Debug 与 Release 各 932 个用例、659,156 条断言全绿，`/W4 /WX /utf-8`
  零警告。
- **基准（`geom_bench_m3c2`，Release，100 次/操作）。** 三个操作分别计时
  `boolean_solid_curved` 与 `detail::curved_solid_volume`，无时间门禁；实测
  中位数与完整环境见 `docs/benchmarks/m3c2.md`。
- **修剪网格化仍不支持布尔结果。** 布尔输出的修剪面边界是非直边/含孔洞，
  `tessellate` 会按 ADR-0027 的域形状确定性拒绝；体积只能走
  `detail::curved_solid_volume`。通用修剪网格化是 M3c-3 的工作。
- **限制清单。** 偏移平行/倾斜环切平面、内部闭环、开放切链、锥顶点、球 cap、
  NURBS 操作数、精确 torus/cone 愈合都留给 M3c-3；环面的分类网格仍按默认
  弦高构建，驱动没有采样预算旋钮。
- **复用。** 合成域与剥边逻辑在球/锥具备裁剪能力后可直接复用；8 ulp 归一化带
  对所有环切入口一致，包括 `curved_boolean_pair_circles`。

## 重新评估的触发条件

M3c-3 加入球/锥操作数、一般布尔结果网格化、偏移平行或倾斜环切平面、内部闭环
或开放切链、NURBS 操作数时；环面切链需要精确圆弧愈合（当前是弦逼近，残差
`<= 1e-3`）时；体积积分的环简化容差需要随几何尺度调整时。
