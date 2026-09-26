# 计算索引适配契约

## 边界

`scripts/build_formation_energy_bundle.py` 是离线结构化数据适配器，不是 VASP 解析器或任务索引客户端。它只打开 `--calculation-index` 与 `--phase-analysis` 指定的 JSON；`result_ref`、`source_ref.file` 和 `locator` 只作为可追溯来源保留，绝不据此遍历目录或读取 `OUTCAR`、`OSZICAR`、`vasprun.xml`、POSCAR 等原始文件。

VASP 工程师或上游适配器负责提供并验收下列数值。结构来源可引用 `doping-supercell`、`defect-complex-builder` 或 `defect-screening` 产物；本适配器不创建或修改结构，不从 POSCAR 推导 `n_i`。任务身份必须来自已接受的计算索引和完成证据，不根据目录名、路径、用户名或作业输出推断归属。

## JSON 契约

根对象必须包含：

| 字段 | 要求 |
|---|---|
| `schema_version` | 整数 `1` |
| `status` | 严格为 `accepted` |
| `acceptance` | `status=accepted`、`accepted_at` 和索引本身的 `source_ref` |
| `material_id` | 与纯净宿主、缺陷态一致的稳定材料标识 |
| `calculation_method` | 含 `HSE06` 的能量方法 |
| `pristine` | 共享宿主能、超胞晶格、带边和各自来源 |
| `states` | 一个或多个已验收缺陷态记录 |
| `elemental_references` | 所有参与元素的 HSE06 eV/atom 参考及来源 |

每个 `source_ref` 至少为：

```json
{
  "task_id": "accepted-task-id",
  "file": "result-or-manifest-reference",
  "locator": "field/table/record locator"
}
```

`task_id` 必须和引用所属 task 完全一致。顶层 `acceptance.source_ref` 是索引级来源，可不含 `task_id`。

每个 task（任务）对象必须包含非空的 `project_id`、`run_id`、`cluster`、`task_id`、`step_id`、`attempt_id`、`result_ref`、`accepted_at`，并且 `status=success`、`acceptance_status=accepted`，以及含相同 `task_id` 的 `acceptance_ref`。`job_id` 可同时提供；若索引代表 SLURM（作业调度系统）提交，应使用真实提交回执中的 ID。纯净宿主、每个缺陷态及每个单质参考都各自提供 task 身份，任一缺失或未验收即拒绝整批组装。

`pristine` 至少包含：

```json
{
  "task": {},
  "E_host_ev": -100.0,
  "E_host_source_ref": {},
  "supercell_lattice_angstrom": [[8, 0, 0], [0, 9, 0], [0, 0, 10]],
  "lattice_source_ref": {},
  "band_edges": {
    "calculation_method": "HSE06",
    "calculation_type": "self_consistent",
    "kpoints": {"scheme": "Gamma", "mesh": [1, 1, 1]},
    "converged": true,
    "vbm_reference_ev": -5.0,
    "cbm_reference_ev": -3.0,
    "gap_ev": 2.0,
    "source_refs": {
      "vbm_reference_ev": {},
      "cbm_reference_ev": {},
      "gap_ev": {},
      "calculation_type": {},
      "converged": {},
      "kpoints": {}
    }
  }
}
```

每个带边来源引用指向同一已验收的纯净宿主 task。`gap_ev` 必须为正且与 `cbm_reference_ev - vbm_reference_ev` 一致。宿主 `E_host_ev` 是组装后所有缺陷态共用的参考能；其来源引用会随每态保存在输出的 `origin.calculation_index.state_provenance`。

每个 `states[]` 条目必须提供：

```json
{
  "state_id": "state-id",
  "defect_id": "defect-id",
  "geometry_id": "audited-geometry-id",
  "charge": 0,
  "n_i": {"A": -1, "B": 1},
  "E_def_ev": -96.0,
  "E_corr_ev": 0.0,
  "E_corr_status": "justified_zero",
  "E_corr_description": "Explicit reason for the zero correction.",
  "material_id": "host-id",
  "calculation_method": "HSE06",
  "supercell_lattice_angstrom": [[8, 0, 0], [0, 9, 0], [0, 0, 10]],
  "task": {},
  "source_refs": {
    "state_id": {},
    "defect_id": {},
    "geometry_id": {},
    "charge": {},
    "n_i": {},
    "E_def_ev": {},
    "E_corr_ev": {},
    "supercell_lattice_angstrom": {}
  }
}
```

`n_i` 是结构审核者提供的带符号整数原子差：加入为正、移除为负；适配器只验证和转交，不计算结构差。`E_corr_ev` 必须显式给出并附来源。`E_corr_status` 只能是 `applied` 或 `justified_zero`；后者要求数值为零且有非空理由。缺失修正不得编码成零。

`elemental_references` 按元素符号索引，每项包含 `energy_ev`、`method`、`unit=eV/atom`、稳定的 `reference_id`、已验收的 `task` 和 `source_ref`。必须覆盖宿主全部元素、所有 `n_i` 物种及相图中出现的外来固定化学势物种；任何缺项或非 HSE06 参考都会阻断组装。

## 相图输入与输出

相图结果必须为 `schema_version=1`、`status=accepted`，并带有完整宿主元素列表、约束集合和 `selected_points`。批处理默认只消费 `selection_source=automatic_extreme_face_convex_interpolation` 的点。每个点都必须给出精确覆盖宿主所有元素的 `delta_mu_eV`、控制元素/分数（自动点）、逐条 `constraint_margins_eV` 和 `constraints_source`。适配器独立重算宿主等式、单质上界及每条竞争相不等式，并逐项比对记录余量；不接受不完整向量、投影坐标或伪造的余量。

O 优先、唯一卤素默认及歧义后的 `requires_explicit_control_element` 状态由相图程序记录。歧义状态没有自动点时，批处理不生成空白或猜测点；用户需指定控制元素重跑相图，或以 `--point-name` 明确选择有效人工审计点。

成功后，每个选中点单独输出 `formation-input.json`、formation `validate_input` 接受的分析审计 JSON、SVG 和 PNG。输入的 `origin` 保留索引路径与验收、task 身份、每个数值的来源定位、元素参考 ID，以及相图文件路径、点向量、算法和约束余量。混合化学势明确为 PBE 相对点加 HSE06 单质参考；它不代表 HSE06 稳定域。

## 失败关闭

以下情况必须在写任何 bundle 文件前失败：索引或顶层验收不通过；任一 task/attempt/step 未验收或缺接受回执；宿主能量、晶格、VBM、CBM 或带隙缺失/矛盾；缺陷 `E_def_ev`、`E_corr_ev`、修正状态/理由、`charge`、`n_i`、几何或来源引用缺失；元素参考缺失或方法不匹配；相图点缺少完整向量、约束源或一致余量；点违反任一稳定域约束；或者输出目录已含文件。

这里的验收针对结构化输入契约和已声明的来源关系，不替代科研人员对相图完整性、结构/电子态物理可比性、有限尺寸模型、修正方法或实验可达性的判断。
