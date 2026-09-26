---
name: defect-complex-builder
description: "从明确替位点准备掺杂与复合缺陷审核清单，可显式生成默认关闭的孤立填隙掺杂待审候选；或在已验收替位超胞上构建 1-3 个附加组分并保留原子映射和组成差分。"
---
# 复合缺陷结构构建

复合候选构建器只消费 `doping-supercell` 已通过验收的替位结果，在同一超胞中增加
1-3 个 vacancy（空位）、interstitial（填隙）或 antisite（反位）组分。它不重新生成超胞、不重新
选择替位位点、不做结构弛豫，也不把少量初始构型宣称为全局基态搜索；
下述结构阶段编排器则依次调用已有的独立结构生产者。

## 从宿主结构准备可计算的外掺杂结构清单

`scripts/prepare_extrinsic_inventory.py` 是结构阶段编排器，依次**调用**
`doping-supercell` 的两阶段构建和本 skill 的复合候选构建；它不直接改写原子、
不执行 VASP（第一性原理软件）作业。用户指定母相 `host_poscar`、正对角
`supercell_matrix`、替位原子及未掺杂超胞中的**精确**
`substitution.site_1based`（从 1 起），程序通过几何映射读取母相站点及完整
Wyckoff（晶体学位点标记）来核对，不默认从多个站点猜一个。可选
`complex.operations`；默认不生成任何附加组分，也不生成**孤立外源填隙掺杂**。
显式设置 `isolated_interstitial.enabled=true` 后，编排器从
`defect-screening` 仅采集通过筛选的外源填隙首批候选（默认至多 2 个），
逐个验证其与替位使用的**同一宿主、同一超胞**，记录逐位点原生清单、
POSCAR 哈希、几何坐标及相对纯净超胞的组成差。生成器自带的 pristine 副本
只在几何与位点身份核验通过后，把报告 pristine 身份重绑到替位使用的冻结
POSCAR，并保留生成器原始路径与哈希供审计；氧化态仅是生成候选所需的
假设，建议电荷不能代替实际 `charge_audit`／`spin_audit`。这一步只产出待审
候选；审核后才能进入 VASP 可接收清单，不能把“替位＋填隙”复合物当作孤立
填隙掺杂。

```json
{
  "schema": "extrinsic-structure-recipe-1",
  "host_poscar": "/absolute/host/POSCAR",
  "supercell_matrix": [[2,0,0],[0,2,0],[0,0,2]],
  "substitution": {"from": "A", "to": "D", "site_1based": 5},
  "complex": {
    "operations": [{"id": "vacancy_B", "kind": "vacancy", "species": "B"}],
    "selection": {"max_structures": 2}
  }
}
```

仅当需要孤立填隙掺杂时，在上述配方顶层另加：

```json
"isolated_interstitial": {
  "enabled": true,
  "host_oxidation_states": {"A": 1, "B": -1},
  "dopant_oxidation_state": 3,
  "max_candidates": 2
}
```

`host_oxidation_states` 可省略，但必须能从母相唯一确定中性氧化态，否则
上游会要求审查且无合格候选时失败。`dopant_oxidation_state` 必须显式填写。
草稿中的 `isolated_interstitial_source` 和 `isolated_interstitials` 是待审
结构数据。未明确启用时既不生成也不进入清单；不得从有填隙候选的草稿删除
字段后冒称通过填隙接入。

```bash
python3 "$SKILL/scripts/prepare_extrinsic_inventory.py" \
  prepare /absolute/recipe.json /absolute/new-structure-directory
```

此步仅输出原生结构记录与 `review-draft.json`，状态是
`PENDING_SCIENTIFIC_REVIEW`。研究者从草稿选结构，对每个待算
`state_id` 提供实际批准的 `charge_audit`、`spin_audit`；选复合结构还须
提供精确结构哈希绑定的 `complex_candidate_audits`。草稿含孤立填隙时，还须
在 `interstitial_candidate_audits` 中按候选结构 ID **逐个**给出明确
`status=approved` 或 `status=rejected`：批准项绑定来源报告哈希、候选 POSCAR
哈希、位点 ID、物种、筛选记录、组成差与分数坐标；拒绝项须记录原因。
缺审核、来源或候选哈希变化、点位/母超胞不一致、获批候选未被 `states[]`
选择时一律拒绝，且不会写出清单。评论不能只给 `charge_hint`、不能自填默认
批准。审核文件是
`schema=extrinsic-inventory-review-1`，带 `draft_sha256`、`states[]`
（各含 `structure_id/defect_id/geometry_id/state_id` 与完整两种审核）
、`complex_candidate_audits`（按所选结构 ID 索引）以及启用填隙时的
`interstitial_candidate_audits`。两类审核包含
`status=approved`、`reviewer`、`evidence_ref`、`structure_sha256`，
电荷审核还含经审核整数 `q`，自旋审核含
`spin_state: {"nupdown": 整数}` 或
`{"constraint": "unconstrained"}`。填隙获批后会写入清单顶层
`interstitial_structures`（含 `structure_id/provenance_path/
provenance_sha256/candidate_path/candidate_sha256/site_id/selection/screen/
composition_delta_from_pristine/site_audit`），与替位、复合条目共享同一
pristine 来源和冻结哈希清单。审核后运行：

```bash
python3 "$SKILL/scripts/prepare_extrinsic_inventory.py" \
  freeze /absolute/new-structure-directory/review-draft.json \
  /absolute/approved-review.json /absolute/new-inventory.json
```

冻结器验证宿主/纯净/替位/复合/填隙结构及来源没有自草稿漂移，然后产生
`vasp-extrinsic-structure-inventory-1`。将其作为
`vasp-calculation` 的 `workflow.source_inventory`；原生计算端会重新验收
每个来源、审核及冻结文件。旧本征作业的跨请求复用不是此结构清单的功能，
须按交接报告的身份/方法一致性门另做，不能仅按目录或组成拼接结果。

## 高频流程

1. 先用 `doping-supercell` 完成替位，并取得
   `doping-supercell.doping-provenance/v2` 来源记录。
2. 为每个附加组分提供候选：
   - 空位可给出已掺杂结构中的 1-based（从 1 开始）原子编号；省略时枚举指定物种，
     优先选择离掺杂中心较近者；不会删除替位中心。
   - 反位给出 `from_species`（被替换的宿主物种）及 `species`（替入物种，须已在超胞中），
     可显式给候选编号，或从近邻自动枚举；禁止覆盖掺杂中心。
   - 填隙可显式给分数/笛卡尔坐标及 `geometry_source.method/details`，也可改用
     `generate` 在已掺杂结构上测量中心—配体键长，用邻近配体平均距离略缩短的
     球壳做周期性空隙搜索。生成器只给 1–2 个非碰撞起始位置，记录实测键长、
     目标距离、最小间距和搜索方向；缺少可测参考配体时失败，不凭元素半径猜位置。
3. 运行候选排序；先审计分项得分和未知信息，再构建少量代表结构。
4. 对生成结构做后续平移、键长比例、旋转、镜像或冻结操作时，调用
   `atom-adjust/scripts/adjust.py`，不要在本 skill 中另写坐标编辑逻辑。

```text
SKILL=/path/to/defect-complex-builder
python3 "$SKILL/scripts/build_defect_complex.py" rank config.json
python3 "$SKILL/scripts/build_defect_complex.py" build config.json
```

`rank` 不写结构；`build` 默认只写排名前 3 个，`selection.max_structures` 最大为
12。候选组合数默认上限 500，超过时要求先缩小候选集，不静默展开批量初筛。

## 配置契约

顶层必须提供：

```json
{
  "doping_provenance": "doping_provenance.json",
  "center_charge_hint": 1,
  "operations": [
    {
      "id": "cl_i",
      "kind": "interstitial",
      "species": "Cl",
      "charge_hint": -1,
      "electrostatic_basis": "Sb_Zn donor compensation hypothesis",
      "candidates": [
        {
          "id": "near_sb",
          "fractional_coords": [0.45, 0.50, 0.47],
          "geometry_source": {
            "method": "coordination-guided",
            "details": "candidate placed in an open Sb coordination direction"
          },
          "bonding": {
            "target_sites_1based": [1],
            "ideal_distance_A": 2.5,
            "tolerance_A": 0.4,
            "basis": "representative Sb-Cl distance hypothesis"
          }
        }
      ]
    }
  ],
  "selection": {"max_structures": 3},
  "output_dir": "complex_candidates",
  "provenance": "complex_candidates/provenance.json"
}
```

`base_poscar` 可指向搬迁后的相同 POSCAR（结构文件），但其 SHA256（内容校验值）
必须与替位来源记录一致。`center_site_1based` 省略时使用替位输出站点。
`charge_hint` 仅允许 `-1`、`0`、`1` 或 `null`；它是排序提示，不是缺陷真实电荷证明。
未报告的电性或成键信息保留为 `unknown` / `not_reported`，不得补猜。

空位操作使用 `candidate_sites_1based`；填隙候选的
`fractional_coords` 与 `cartesian_coords_A` 二选一。每个操作和候选 ID 必须唯一。
与显式 `candidates` 互斥的填隙自动模式示例：

```json
{
  "id": "cl_i", "kind": "interstitial", "species": "Cl",
  "generate": {
    "reference_species": "Cl",
    "reference_neighbor_count": 4,
    "shortening_fraction": 0.05,
    "max_candidates": 2,
    "minimum_clearance_A": 1.2
  }
}
```

`reference_neighbor_count` 取实际测到的邻近配体数，默认 4；若该中心没有足够多
该物种配体，改由研究者提供有来源的显式位置。`shortening_fraction` 是起始几何
建议而非普适键长定律，不能凭它证明构型稳定。反位配置举例：
`{"id":"cl_cs","kind":"antisite","from_species":"Cs","species":"Cl"}`。

## 候选选择

排序总分越低越优，默认权重为：

```text
score = 0.50 * center_distance
      + 0.25 * electrostatic
      + 0.25 * possible_bonding
```

- `center_distance`：候选到替位中心的周期最短距离，经
  `ranking.distance_scale_A` 归一化。
- `electrostatic`：显式电荷提示异号为有利、同号为不利；零或未知为中性/未知。
  这是候选启发式，不替代正式 charge state（电荷态）定义。
- `possible_bonding`：按候选提供的目标原子、理想键长和容差计算；缺失时记为未知。

每项都写入排序记录，平分时按候选 ID 稳定排序。组合中重复删除同一原子、
在同一原子上删除和反位、填隙彼此或与保留原子距离小于
`selection.minimum_separation_A` 的构型会被拒绝。
自动枚举空位/反位按 `selection.environment_tolerance_A`（默认 0.12 Å）比较
候选到掺杂中心的距离和逐物种近邻距离，几何非常相近者暂缓并记录
`representative_id`；研究者**显式指定**的位点不被此阈值自动合并。
此外自动枚举只将离中心最近距离加
`selection.nearest_shell_window_A`（默认 0.35 Å，可调整）以内的候选纳入
组合；更远且环境不同的点记录为 `deferred_distinct_distant_unrepresented`，
不是物理等价或永久排除。显式指定候选不受近邻窗口筛除。此窗口只控制
首轮探索成本，不代表结合能或激发发射能的普适阈值。
暂缓并不代表弛豫后形成能或光谱等价，发现差异时应恢复计算。
默认只取总分最高的少量构型；需要更广的畸变或基态搜索时转交
ShakeNBreak（缺陷结构搜索工具），不在这里扩展批量扫描。

## 验收门

构建前必须满足：

- 替位 provenance 的 schema、状态、阶段和 `validation.passed` 正确；
- 输入 POSCAR 哈希、原子数和组成与 provenance 一致；
- 附加组分恰为 1-3 个，索引、物种、坐标和几何来源完整；
- 候选组合未超过上限，且至少有一个通过距离与重复位点检查。

写出后必须重新读取并验证：

- 原子数和组成差分由实际结构自动计算；
- 若上游替位 provenance 含可验收的纯净超胞，须验其 POSCAR 哈希与晶格；
  除“相对已掺杂基底”差分，还逐候选记录
  `composition_delta_from_pristine`（形成能的有符号 `n_i` 来源）。
  老来源缺此记录时写 `null`，后续形成能输入组装必须另取已审计的
  纯净超胞结构差分，不能误将相对掺杂基底的差分当作 `n_i`；
- 原有原子保留 `source_site_1based` 映射，新增填隙明确标为无源站点；
- 物种块顺序稳定，新增同种填隙并入原物种块；
- 每个输出记录候选 ID、组分、几何来源、分项得分、路径和 SHA256。

成功来源记录使用 `defect-complex-builder.provenance/v1`。该记录证明结构构建
过程可追溯，不证明所选构型是最低能，也不授权 VASP（第一性原理软件）提交。
