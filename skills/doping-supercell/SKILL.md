---
name: doping-supercell
description: "两阶段生成替位掺杂超胞：先验收未掺杂超胞，再按精确位点或中心近邻替位，保留父映射、物种顺序与哈希溯源。用于替位掺杂及不等价位点结构。"
---
# 替位掺杂超胞生成

该 skill 只处理替位掺杂的两阶段结构生成：

1. `build-supercell`：母相 POSCAR（结构文件）-> pristine（未掺杂）超胞 + manifest
   （清单）。
2. `dope`：已接受的 pristine 超胞 + manifest -> 掺杂 POSCAR + provenance
   （来源记录）。

第二阶段永远不重新搜索矩阵、不扩展结构；必须先检查并接受第一阶段的
POSCAR 和 manifest。运行时脚本与 Schema（结构模式）是当前公开契约：
[Schema/兼容说明](references/schema-and-compatibility.md) 和 `schemas/` 不得删改。
当前公开快照未随附专用 `tests/` 夹具，因此不能声称已有该目录的测试覆盖。

## 命令行入口

```text
SKILL=<DOPING_SUPERCELL_SKILL_DIR>
python3 "$SKILL/scripts/dope_supercell.py" build-supercell <BUILD_CONFIG_JSON>
python3 "$SKILL/scripts/dope_supercell.py" dope <DOPE_CONFIG_JSON>
```

`build` 是 `build-supercell` 的别名，`substitute` 是 `dope` 的别名；子命令
必填，不存在组合式“扩展并掺杂”命令。

## 第一阶段：构建未掺杂超胞

### 两种互斥模式

| 模式 | 必须有 | 禁止有 |
|---|---|---|
| `auto` | `min_atoms`、`min_image_distance` | `supercell_matrix` |
| `explicit` | 非奇异 3x3 整数 `supercell_matrix` | `min_atoms`、`min_image_distance` |

母相路径使用 `parent_poscar`。自动模式调用
`doped.generation.get_ideal_supercell_matrix`，且始终以用户提交的精确母相
为输入，不静默切换原胞/标准化结构，因为矩阵必须可逆向映射到该母相。

显式矩阵直接采用：

```text
L_supercell = M @ L_parent
```

其中 `M` 是 3x3 整数超胞矩阵，`L_parent` 和 `L_supercell` 是按当前
`pymatgen`（材料结构库）约定存储的晶格矩阵。矩阵形状、整数性和行列式
均需通过检查。

### 第一阶段接受门

构建器始终验证：

- `N_supercell = abs(det(M)) * N_parent`；
- 化学计量等于母相组成乘以 `abs(det(M))`；
- species-block order（连续物种块顺序）与母相一致；
- 写出的 POSCAR 重新读取后仍通过同一组检查；
- 自动模式满足两个用户阈值。

`checks.formula`、`checks.atoms` 和 `checks.species_order` 可声明精确期望；
未提供时由母相和矩阵派生。成功清单的结构模式为
`doping-supercell.supercell-manifest/v1`，至少包含母相/超胞路径、SHA256（内容校验值）、
组成、原子数、物种顺序、模式、请求参数、矩阵、行列式、最小镜像距离、
`parent_preserved: true` 和验证结果。

只有 `status == "passed"` 且 `validation.passed == true` 才能进入第二阶段。

## 第二阶段：父结构约束掺杂

配置要求 `supercell_manifest` 和恰好一次 `substitution.from -> to` 替换；
`substitution.parent_site` 指定母相编号 `index_1based` 与 `wyckoff`，
按需给出 `mapping_tolerance`、`symprec`。`output` 指定掺杂结构，
`provenance` 指定来源记录；完整字段见 `schemas/dope.schema.json`，
旧字段迁移时再读[兼容说明](references/schema-and-compatibility.md)。

`supercell_poscar` 可指向搬迁后的 pristine 副本，但 SHA256 必须与 manifest
一致；`parent_poscar` 也只允许使用哈希一致的搬迁副本。manifest 路径是
默认来源。

第二阶段在替换前验证清单的结构模式/类型/状态、清单哈希、母相和
未掺杂结构哈希、原子数、化学式、晶格、矩阵及物种顺序。`checks.formula`、
`checks.atoms`、`checks.species_order` 三项均为必填，并在写入前后各检查一次。

## 站点选择

### 自动模式

未提供任何精确超胞索引时，默认使用
`nearest_supercell_center`：

1. 枚举满足 `from`、母相站点编号、Wyckoff（晶体学等价位点标记）字符串和可选母相分数坐标约束
   的所有平移等价候选；
2. 计算候选在写出 POSCAR 平行六面体内到分数坐标
   `(0.5, 0.5, 0.5)` 的直接笛卡尔欧氏距离；
3. 周期性最短镜像距离只作为诊断，不参与排序；
4. 距离相差不超过 `1e-10 Å` 时视为并列，选择最低 pristine
   `site_1based`。

该距离定义适用于正交和倾斜晶胞，不能换成逐分量分数坐标范数。
provenance 必须记录候选全集、映射结果、两种距离、并列容差和最终选择。

### 精确模式

可使用 `site_1based` 或 `site_0based`，只能提供一种。精确索引优先级最高：
它只验证指定原子，不枚举其他候选，也不会被更靠近中心的候选替换。两种
索引同时存在、与 `site_selection` 合用、布尔/负值/非整数/越界输入都必须失败。

`parent_site.wyckoff` 必须使用实际 `SpacegroupAnalyzer`/`spglib`（空间群分析库）返回的完整
字符串，不要自行把诸如 `4c` 归一化成 `c`。父站点编号和 Wyckoff 条件缺失
时，严格第二阶段拒绝执行。

## 父站点映射与输出顺序

选择后使用行向量分数坐标映射：

```text
f_parent = f_supercell @ M mod 1
```

映射必须唯一，并通过目标物种、母相编号、母相物种、Wyckoff、周期性映射
距离和可选母相分数坐标检查。不同文件哈希、超胞形状或局部原子号本身
不能证明是不同晶体学位点；父结构映射是权威依据。

替换前，目标宿主原子移到其连续宿主物种块的最前端，再替换为掺杂原子，
因此不会产生分裂的 `Host / Dopant / Host` 块。若宿主物种不是一个连续块，
拒绝猜测并失败。最终物种顺序必须由 `checks.species_order` 固定。

## 失败关闭与兼容边界

- 第二阶段配置出现 `mode`、`supercell_matrix`、`min_atoms` 或
  `min_image_distance` 时失败；不能重新定义或生成超胞。
- 严格上下文的未知字段、拼写错误和未支持的选择字段失败，不能静默回退
  到自动模式。
- 新 provenance 使用
  `doping-supercell.doping-provenance/v2`，并保留源 manifest、母相、
  未掺杂结构、矩阵、站点选择、父映射、最终检查和掺杂 POSCAR 哈希。

读取旧清单或调用旧接口时，查阅[结构模式与兼容接口](references/schema-and-compatibility.md)；
旧接口不属于严格第二阶段，不能用其绕过接受门。发布前至少运行顶层隐私扫描，
并对 Stage 1/Stage 2 各使用一份本地占位配置做显式验证；当前公开快照不宣称
存在未随仓库发布的专用测试夹具。
