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
- 单个 op 通过的条件：`ms <= 5000` 且 `volume_rel_err <= 1e-9` 且结果通过
  `validate_topology`；任一行 FAIL，进程退出码为 1。

## 冒烟：200 边

`geom_bench_m3.exe 200 all`，退出码 0，总墙钟 4.4 s（其中两个操作数装配
66 ms，不计入计时）。输出：

```
QA8 op=union faces_a=202 faces_b=202 ms=1196.089 volume_rel_err=1.810e-16 target=5000ms/1e-9 PASS
QA8 op=intersection faces_a=202 faces_b=202 ms=1075.615 volume_rel_err=1.810e-16 target=5000ms/1e-9 PASS
QA8 op=difference faces_a=202 faces_b=202 ms=1296.292 volume_rel_err=1.086e-15 target=5000ms/1e-9 PASS
```

200 边时三次运算都在 1.3 s 内，恒等式误差 1e-16 量级，说明装配、布尔与体积
测量路径本身正确。

## 官方：5000 边

`geom_bench_m3.exe 5000 all`，三次运算依次耗时 636.0 s、550.4 s、559.3 s，
布尔总计 1745.8 s（29.1 min）；加上两个操作数装配（35.4 s）与校验，进程总墙钟
29.7 min。输出：

```
QA8 op=union faces_a=5002 faces_b=5002 ms=636040.265 volume_rel_err=5.926e-14 target=5000ms/1e-9 FAIL
QA8 op=intersection faces_a=5002 faces_b=5002 ms=550446.162 volume_rel_err=5.926e-14 target=5000ms/1e-9 FAIL
QA8 op=difference faces_a=5002 faces_b=5002 ms=559296.597 volume_rel_err=5.880e-14 target=5000ms/1e-9 FAIL
```

三次布尔调用都正常返回（否则按约定退出码为 2）；两个操作数的实测体积为
`Va = 1256.6367307022517`、`Vb = 1256.6367307022515`，恒等式相对误差 5.9e-14，
满足 1e-9。**耗时全部超 5 s：并集约 127.2 倍、交集约 110.1 倍、差集约 111.9 倍，
因此 QA-8 为 FAIL。**

## 结论

1. **QA-8 FAIL（仅性能）**：5000 边并/交/差分别 636.0 s、550.4 s、559.3 s，
   门禁 5 s；恒等式 5.9e-14 已满足 1e-9，拓扑有效。失败是纯性能问题，
   与 M1b 的 QA-1 结论一致。
2. 与先前的二次增长测量相符：1000 边并集 29.8 s，5000 边 636 s，边数 5 倍而
   耗时 21.3 倍（纯二次预测 25 倍），说明成本仍由逐对候选/切割主导。
3. 同机此前的探针运行测得并集 563,564 ms，本次为 636,040 ms（+13%）。15 W
   移动处理器在长时间满载下降频会带来这一量级的波动；门禁按计划固定值判定，
   仍为 FAIL。
4. 本次只新增基准程序，未改动库与测试；Debug 全量用例保持通过（见提交说明）。

## 复现

```powershell
cmake -S . -B build-rel -DGEOM_WERROR=ON -DGEOM_BUILD_BENCH=ON
cmake --build build-rel --config Release --target geom_bench_m3
build-rel\Release\geom_bench_m3.exe 200 all   # 冒烟，约 4 s
build-rel\Release\geom_bench_m3.exe 5000      # 官方，约 30 min
```

基准在 stderr 打印每个 op 的 `start` / `done` 进度与装配耗时，stdout 只保留上述
QA8 行，便于长时间运行时观察进度而不污染记录。
