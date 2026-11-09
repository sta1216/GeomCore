# M3a QA-8 平面布尔基准（2026-11-10）

本文件记录 M3a 出口 QA-8 的机器、构建与实测数据：两个 5000 边棱柱的并、交、差
布尔运算耗时与体积恒等式。数字为单机单次测量，门禁按计划固定值判定，不是跨机器
承诺。**结论：体积恒等式在 5.9e-14 内成立，但三次运算的耗时都远超 5 s 门禁
（约 110–127 倍），QA-8 判定为 FAIL。**

## 环境

| 项 | 值 |
| --- | --- |
| 日期 | 2026-11-10（Asia/Shanghai） |
| OS | Windows 10 企业版 10.0.19045，x64 |
| CPU | Intel Core i7-8550U @ 1.80 GHz（4 核 8 线程，15 W 移动处理器） |
| 内存 | 15.9 GB |
| 编译器 | MSVC 19.44.35217.0（VS 2022 17.14，工具集 14.44.35207），`/W4 /WX /utf-8 /fp:precise` |
| 构建 | CMake 4.4.4，VS 2022 generator，`-DGEOM_WERROR=ON -DGEOM_BUILD_BENCH=ON` |
| 配置 | Release（`/O2`），单线程 |

## 口径

- 操作数：两个正 5000 边形棱柱，半径 10、高 4、焊接容差 1e-9；A 轴心在
  `x = 0`，B 轴心在 `x = 0.25 * radius`，同轴、z 对齐。每个棱柱由
  `assemble_solid` 装配：底盖（逆序环、法向 `(0,0,-1)`）、顶盖（正序环、法向
  `(0,0,1)`）与 5000 个侧面四边形，共 5002 个面。
- 装配后逐个校验 `validate_topology`、`validate_tolerances`、面数
  `sides + 2`，以及 `mass_properties` 体积与解析值
  `0.5 * n * r^2 * sin(2π/n) * h` 的相对误差 ≤ 1e-12；任一不满足即退出码 2。
- 计时只覆盖 `boolean_solid(a, b, op, {})`，用 `std::chrono::steady_clock`，
  毫秒保留三位小数；装配、校验与结果体积测量都在计时区间之外。
- `all` 模式在相同操作数上依次运行并、交、差。并集与交集行使用
  `|Vu + Vi - Va - Vb| / (Va + Vb)`，差集行使用 `|Vd + Vi - Va| / Va`。
- 通过条件分两种：单个 op 只看 `ms <= 5000` 与 `validate_topology`，恒等式在
  单 op 下不可用（该行打印 `identity=unavailable`、`volume_rel_err=nan`），
  1e-9 不是单 op 门禁；`all` 模式的每一行还要求对应恒等式
  `volume_rel_err <= 1e-9`。任一行 FAIL，进程退出码为 1。

## 冒烟：200 边

`geom_bench_m3.exe 200 all`，退出码 0，总墙钟 3.2 s（其中两个操作数装配
65.8 ms，不计入计时）。输出：

```
QA8 op=union faces_a=202 faces_b=202 ms=1064.064 volume_rel_err=1.810e-16 target=5000ms/1e-9 PASS
QA8 op=intersection faces_a=202 faces_b=202 ms=991.887 volume_rel_err=1.810e-16 target=5000ms/1e-9 PASS
QA8 op=difference faces_a=202 faces_b=202 ms=1053.113 volume_rel_err=1.086e-15 target=5000ms/1e-9 PASS
```

200 边时三次运算都在 1.3 s 内，恒等式误差 1e-16 量级，说明装配、布尔与体积
测量路径本身正确。冒烟行里的 PASS 只证明这条路径能跑通，不构成 QA-8 结论；
QA-8 只由 5000 边官方运行判定。

## 官方：5000 边

`geom_bench_m3.exe 5000 all`，三次运算依次耗时 636.0 s、550.4 s、559.3 s，
布尔总计 1745.8 s（29.1 min）；加上两个操作数装配（35.4 s）与校验，进程总墙钟
29.7 min。输出：

```
QA8 op=union faces_a=5002 faces_b=5002 ms=636040.265 volume_rel_err=5.926e-14 target=5000ms/1e-9 FAIL
QA8 op=intersection faces_a=5002 faces_b=5002 ms=550446.162 volume_rel_err=5.926e-14 target=5000ms/1e-9 FAIL
QA8 op=difference faces_a=5002 faces_b=5002 ms=559296.597 volume_rel_err=5.880e-14 target=5000ms/1e-9 FAIL
```

三次布尔调用都返回了结果：`run_op` 在 `boolean_solid` 失败时会打印
`error=...` 并以退出码 2 结束，不会走到 QA8 行（bench/bench_m3.cpp:164-171、
:288、:295）。结果拓扑也在同一进程内、打印之前检查：`run_op` 用
`validate_topology(brep).ok()` 设置 `topology_ok`（bench/bench_m3.cpp:176），
`print_qa8` 把它计入 PASS/FAIL（bench/bench_m3.cpp:191-193），因此行内判定已经
包含拓扑校验，不需要读者另行推断。两个操作数的实测体积为
`Va = 1256.6367307022517`、`Vb = 1256.6367307022515`，恒等式相对误差 5.9e-14，
满足 1e-9。**耗时全部超 5 s：并集约 127.2 倍、交集约 110.1 倍、差集约 111.9 倍，
因此 QA-8 为 FAIL。**

## 结论

1. **QA-8 FAIL（性能）**：5000 边并/交/差分别 636.0 s、550.4 s、559.3 s，
   门禁 5 s，超出约 110–127 倍；恒等式 5.9e-14 满足 1e-9。FAIL 的直接证据是
   时间越界（拓扑校验也计入该行判定，但行本身不区分失败来源）。与 M1b 的
   QA-1 结论一致。
2. 与先前的二次增长测量相符：1000 边并集 29.8 s，5000 边 636 s，边数 5 倍而
   耗时 21.3 倍（纯二次预测 25 倍），说明成本仍由逐对候选/切割主导。
3. 同机此前的探针运行测得并集 563,564 ms，本次为 636,040 ms（+13%）。15 W
   移动处理器在长时间满载下降频会带来这一量级的波动；门禁按计划固定值判定，
   仍为 FAIL。
4. 本次只新增基准程序，未改动库与测试；Debug 构建在本机通过，
   `ctest -C Debug` 1/1 passed，35.96 s。

## 范围与限制

M3a 的布尔只覆盖平面面；以下限制是已记录的契约，不是待修 bug。曲面面布尔、
SSI 与 QA-9 随 M3b 交付。

- **仅平面面**：任一面的 surface 不是 `Plane3` 时返回 `Errc::unsupported`，
  消息为 `"boolean_solid: curved faces are deferred to M3b"`；圆柱/圆锥/球/
  圆环实体的布尔与 QA-9 一并延后。
- **pcurve 为直线**：M3a 只产生直边，`Pcurve2d` 是精确的 2D 线段；曲面 pcurve
  随 SSI 一起引入。
- **共面补丁不得触及宿主面边界**：与宿主面环相接的共面补丁会留下 T 型接点，
  诊断为 `Errc::not_closed`（`"assemble_solid: edge is used once"`），不做
  T 型传播或静默缝合。
- **腔体归属用 AABB 预检 + 点/接触探测**：多壳结果（外壳 + 腔体壳）受支持，
  但腔体包含性与不相交判定先用保守的 AABB 预检，再做点/接触探测，不做精确
  壳-壳求交。
- **非流形接触是诊断失败**：边相接触的 L 形对称差诊断为
  `Errc::invalid_topology`，记录为失败而不是缝合。
- **只接受单实体、单壳操作数**：每个输入必须恰好一个实体且该实体只有一个壳，
  结果可为多壳。多实体输入返回 `Errc::unsupported`，消息为
  `"boolean_solid: expected exactly one solid per operand"`；多壳输入返回
  `Errc::unsupported`，消息为
  `"boolean_solid: multi-shell operands are deferred to M3c"`。
- **门禁状态**：QA-8 是 M3a 唯一失败的出口门禁，失败原因仅是性能——并、交、
  差分别 636,040 ms、550,446 ms、559,297 ms，比 5,000 ms 门禁慢 110–127 倍；
  体积恒等式 5.93e-14 / 5.93e-14 / 5.88e-14 远优于 1e-9，正确性通过。Debug 与
  Release 各 516 个用例、455,084 条断言全绿，`/W4 /WX /utf-8` 下零警告。

## 复现

```powershell
cmake -S . -B build-rel -DGEOM_WERROR=ON -DGEOM_BUILD_BENCH=ON
cmake --build build-rel --config Release --target geom_bench_m3
build-rel\Release\geom_bench_m3.exe 200 all   # 冒烟，约 3 s
build-rel\Release\geom_bench_m3.exe 5000      # 官方，约 30 min
```

基准在 stderr 打印每个 op 的 `start` / `done` 进度与装配耗时，stdout 只保留上述
QA8 行，便于长时间运行时观察进度而不污染记录。
