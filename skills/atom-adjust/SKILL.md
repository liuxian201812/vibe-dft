---
name: atom-adjust
description: "编辑 POSCAR（结构文件）的坐标、选择性动力学标志或显式 F 全晶格仿射形变；键操作考虑周期边界，缺陷基态搜索与批量扫描转交 ShakeNBreak（缺陷结构搜索工具）。"
---
# atom-adjust

用于两类结构操作：

- `adjust.py`：对现有 VASP（第一性原理计算软件）`POSCAR` 做一次性的、可审计的坐标或
  Selective dynamics（选择性动力学）编辑。
- `snb_workflow.py`：为 ShakeNBreak（SnB，缺陷结构搜索工具）写配置并打印
  `snb-*` 命令；它不调用 SnB API（应用程序接口），也不直接提交作业。

## 路由边界

| 需求 | 入口 |
|---|---|
| 单原子/小群组的平移、键长比例、旋转、镜像、绝对定位、冻结/释放 | `adjust.py` |
| 显式变形梯度 F、固定分数坐标的全晶格主动形变 | `adjust.py affine`，见下方独立契约 |
| 缺陷基态搜索、批量键畸变、rattle（随机扰动）、二聚体畸变和跨电荷态分析 | SnB，经 `snb_workflow.py` 路由 |

不在本 skill 中重新实现 SnB 搜索或分析逻辑。按场景选取命令时读取
[references/examples.md](references/examples.md)。

## `adjust.py` 契约

### `affine` 独立契约

```bash
python3 "$SKILL/scripts/adjust.py" affine \
  --F 1 0 0 0 1 0 0 0 1 -o NEW_POSCAR SOURCE_POSCAR
```

`$SKILL` 为本 skill 根目录，示例仅为单位矩阵。F 为固定笛卡尔坐标系中的
无量纲 3×3 变形梯度，九个参数按行输入；A 的三行为晶格矢量，单位为 Å（埃）。
以列向量表示笛卡尔位置 r，主动形变为：

$$
A' = A F^{T}, \qquad r' = F r, \qquad f' = f,\qquad V'/V=\det F.
$$

f 为无量纲分数坐标，V 为晶胞体积。这不是保持实际位置不变的被动换基，
不自动从应变张量推断 F；调用者负责提供完整 F。接受一般正行列式 F，
包括其中的几何旋转，但不旋转电子/自旋态，不做原子松弛、应变能或计算提交。

- 复用 `read_poscar(..., strict_affine=True)`；不另设 POSCAR 解析器。
  仅支持单个有限正 scale（缩放因子）、显式物种及匹配的非负整数计数、
  Direct（分数）或 Cartesian（笛卡尔，含 K 标记）坐标和可选 T/F 标志。
  **负 scale 体积模式、三个 scale、无物种行和未支持尾块在此入口明确拒绝**，
  不走旧命令的 negative-scale（负缩放）警告后继续路径。
- 静态几何来源可携带精确全零的离子速度：坐标之后恰为一行速度模式加 N 行
  三分量数值（N 为原子数）。空白模式行表示 Cartesian；也接受大小写不敏感的
  `Cartesian/C/K` 或 `Direct/D`。不跳过空白行、不猜任意未知模式。
  支持普通十进制及 `E/e` 科学计数法，以 `Decimal`（十进制精确数）判断精确零，
  不用浮点容差；非零（包括浮点下溢非零）、非有限、坏数值、格矢速度、预测器、
  截断或多余行均拒绝。纯空白结尾仍视为没有速度块。
  合法全零尾块的原始字节保留，不由调用者截尾；零在任意线性 F 下仍为零。
  这只扩展静态几何兼容，不是分子动力学状态续算或速度变换接口。
  格式依据：VASP Wiki（官方文档）POSCAR 的 Ion velocities、Lattice velocities、
  MD extra 段，核读版本 `oldid=37632`（北京时间2026-09-15）：
  `https://vasp.at/wiki/index.php?title=POSCAR&oldid=37632#poscar9`。
- 输出统一 scale=1、Direct，保持原分数坐标，不包裹、不排序。Cartesian 源按
  原始晶格正确转译，原正 scale 同时作用于晶格和笛卡尔位置；保留原注释、
  物种/计数行、原子顺序、坐标行尾标志及注释。
- F、源晶格与目标晶格必须有限、右手且可逆；条件数大于 `1e12` 时按数值退化拒绝。
  不支持的格式、奇异/反向晶胞、非有限坐标和不完整标志均在输出前失败。
- 源须为普通非符号链接文件；源目标不得同路径，目标及悬空链接已存在均拒绝。
  在目标目录写临时文件并以硬链接排他发布，不提供原地写回/覆盖/降级选项；
  文件系统不支持硬链接时失败，不覆盖兜底。
- 成功时标准输出给出 JSON（结构化数据）审计：F、前后晶格/体积、物种/计数、
  源目标 SHA256（内容摘要）、坐标模式与保持项。`velocity_tail` 另记精确零判断、
  速度模式、行数、源文件行区间及原文保留政策；整文件摘要绑定尾块来源和输出。
  它只证明结构操作，不认证科学用途。

以下坐标操作契约不适用于 `affine` 的输出和失败边界。

### 输入与坐标

- `--atoms` 使用 1-based（从 1 开始）索引，支持逗号分隔的单点和包含端点
  的范围，例如 `5`、`5,7,9`、`5-10`；空、非整数和越界输入必须失败，
  重复索引会被排序并去重。
- 读取 Direct（分数坐标）和 Cartesian（笛卡尔坐标）POSCAR，接受正、负 scale（缩放因子）；负 scale 按当前晶格
  原样处理并发出警告，不把它静默解释成另一套晶格。
- 内部统一用分数坐标计算，输出保留输入的坐标模式和 scale；写出的分数
  坐标包裹到 `[0, 1)`。
- 注释、元素行、坐标行和已有 Selective-dynamics 标志按行保留，不经过
  `pymatgen`（材料结构库）往返重排。

键方向和最近邻距离使用周期性最短镜像（minimum image）Cartesian（笛卡尔）
距离；倾斜晶胞会搜索邻近晶格像，不能用逐分量取整替代。`plane`、旋转轴和
旋转 pivot（枢轴）使用选定的晶胞内 Cartesian 点，不自动把跨边界片段展开。

### 操作与失败边界

| 操作 | 关键语义 |
|---|---|
| `translate` | `--dist` 与 `--ratio` 互斥；`--ratio` 只能配 `bond:j`，按每个移动原子相对参考原子 `j` 的原始最短镜像向量缩放。 |
| `--dir a\|b\|c` | 沿晶格基矢正向移动。 |
| `--dir bond:j` | 从移动原子组的质心指向 `j` 的周期性最短镜像方向；零向量失败。 |
| `--dir vec:x,y,z` | 使用非零 Cartesian 向量并归一化。 |
| `--dir plane:i,j,k` | 用三点右手叉积求平面法向；三点共线失败。 |
| `--dir perp:j,ref:k` | 将质心到 `k` 的向量对 `bond:j` 做 Gram-Schmidt（格拉姆-施密特）正交化；零垂直分量失败。 |
| `rotate` | Rodrigues（罗德里格斯）公式；正角遵循 `axis_i -> axis_j` 的右手定则；默认绕组质心，也可指定 `--pivot k`。 |
| `mirror` | 关于三点定义的平面镜像；平面退化失败。 |
| `set` | `--unit frac` 或 `--unit ang`；组内第一个原子到目标位置，其余原子整体平移以保持组内几何。 |
| `fix` / `free` | 目标原子设为 `F/F/F` 或 `T/T/T`；没有 Selective dynamics 行时插入该行，并令未指定原子默认为 `T/T/T`。 |

所有操作都必须在输入解析、索引检查和几何向量检查通过后才写文件。解析
错误、零长度方向、退化平面/轴、非法单位或互斥参数冲突返回非零状态。

### 输出与副作用

```text
python3 "$SKILL/scripts/adjust.py" translate --atoms 5 --dir bond:6 --ratio 0.8 --dry-run POSCAR
```

- 无 `-o`/`-c` 时默认原地写回；原地写回默认只在相邻 `.bak` 不存在时创建
  备份。`--no-backup` 可关闭，`-o` 写新文件，`-c` 只向标准输出。
- `--dry-run` 只报告位移和可选的邻居键长，不修改输入或备份。
- `--print-bonds N` 保留操作前邻居身份，再报告操作后的同一对距离；不能把
  变化后的邻居排序位置误当作同一原子。

## `snb_workflow.py` 契约

```text
python3 "$SKILL/scripts/snb_workflow.py" prep --defect DEFECT_POSCAR --bulk BULK_POSCAR
```

- `prep` 写 `snb_config.yaml` 并打印 `snb-generate`；已有配置须显式使用
  `--force` 才能覆盖。
- 在执行 `snb-generate` 前必须按实际母相填写 `oxidation_states`；该字段
  决定 SnB 的配位邻居畸变数量，不能用空字典当作默认化学信息。
- `submit`、`parse`、`regenerate` 和 `groundstate` 只打印对应的 SnB 命令；
  `--submit-command` 只是命令文本的一部分，不授予本 skill 直接提交权限。
- SnB 的命令行、默认畸变网格、随机扰动参数和跨电荷态规则以 SnB 自身版本为准；
  本辅助脚本只负责配置骨架和命令路由。

## 依赖与验证

- `adjust.py` 只需要 Python（编程语言）和 `numpy`（数值计算库）；`snb_workflow.py` 只用
  标准库；SnB 仅在第三方搜索流程中需要安装。
- 运行 `tests/` 验证 POSCAR 格式保持、周期边界、键长比例、旋转/镜像、
  Selective dynamics 和输出模式；测试夹具不是普通示例，不要删除。
