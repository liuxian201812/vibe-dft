# 计算端 `calculation-results-1` 索引接入

## 双层索引

`vasp-calculation/scripts/calculation_results.py` 定义的 `calculation-results-1`
包含 PBE 宿主/竞争相/元素参考、HSE06 单质/纯净超胞/缺陷态记录以及原生任务、
attempt（尝试）、调度与完成来源。完整输入仍优先接受**顶层**
`status=READY`、`ready=true`、无拒绝理由的索引，且每行必须为 `READY`、
producer gate（计算端验收门）`VERIFIED`、numerical audit（数值审计）`TRUSTED`、
电子/离子收敛布尔值均为真及收敛来源精确绑定同一个 task/attempt。
另允许极窄的**修正待审路径**：顶层仍为 `status=INCOMPLETE`、`ready=false`，
但必须有完整且无重复的 `expected_records`、无 `adapter_diagnostics`，每行
producer gate 为 `VERIFIED`、numerical audit 为 `TRUSTED`，且顶层及逐行
`rejection_reasons` 只能精确涉及 `CORRECTION_DECLARATION_MISSING`（PBE 修正声明缺项）
和/或 `E_CORR_MISSING`（缺陷修正缺项）。拒绝理由必须对应到确实缺失的修正字段；
任一方法、带边零点/对齐、来源、收敛或其他缺口都会拒绝。适配只生成新的后处理
索引，不改写原生文件，也不把它重标为 `READY`；输出保留
`source_index_status=INCOMPLETE` 和源文件 SHA-256（安全哈希算法）以供追溯。
本分支不从服务器或 `OUTCAR` 收集索引；计算端若尚未从真实任务导出
这一文件，不能把测试构造的 READY 样例当成生产结果。

`scripts/prepare_correction_drafts.py` 可从上述完整可信的原生索引只读生成
精确身份绑定的 PBE 和缺陷修正**待审草稿**。它不会计算修正、声称有真实
审计来源或自动批准任何决议；缺值为 `null`、批准为 `PENDING`。带电原生
`E_corr` 已有数值时，只预置 `alignment_only` 草稿并保留原值，不构造
第二份修正。草稿不是本节下文的已批准输入，使用前须由独立科学审查者
补齐模型、数值、符号及双侧电势/对齐证据。
若任一 PBE 行的原生 `correction_policy_id` 缺失，草稿保留
`policy_id: null`、`policy_review.status=PENDING`，并列出确实存在的
`known_native_policy_ids`；由独立审查者确认统一策略并批准绑定原生索引的
逐行决议，不能把已知行的 ID 自动填给缺失行。多个不同的已知原生策略 ID
仍是口径冲突，草稿生成失败，待上游解释。原生修正单位也缺失时，
草稿的 `unit` 为 `null`，须由审查者显式确定，不能猜成每计算胞修正。

`scripts/adapt_calculation_results.py` 只读转换上述原生索引，产生：

- `phase-input.json`：已按化学式单元/原子归一并明确加过一次修正的 PBE 总能输入；
- `formation-calculation-index.json`：映射至
  [组装器输入](calculation-index-adapter.md)的已验收任务、HSE06 单质参考、超胞带边和缺陷态。

两个文件都保存原生索引路径和逐字段 locator（定位点）；相图分析会转存
`phase-input.json` 中的修正决议来源。本转换不重算 `n_i`、电荷修正或结构，
不改变计算端的任务归属。原生 `pristine_hse` 角色按其单 Γ 已收敛带边契约映射
为形成能输入的自洽纯净超胞；它不独立读取 `vasprun.xml` 核验原始计算类型。
外部决议的 `source_refs.E_corr_ev` 指向独立缺陷修正决议文件；同时保留原生
`E_corr` 声明及原生 `E_def_ev` 能量来源作为双重证据。

## PBE 修正决议

原生索引允许 PBE 行明确声明 `correction.status=not_applied` 且数值为空。
本 skill **不能**把该状态默认为零。使用独立的修正决议 JSON，逐条给出
在该材料/相参考中采用的**数值、单位、来源及理由**：

```json
{
  "schema_version": 1,
  "policy_id": "与原生 PBE correction_policy_id 完全相同",
  "rows": {
    "host_pbe/host-record-id": {
      "native_status": "not_applied",
      "value_ev": 0.0,
      "unit": "eV/cell",
      "source": "相应修正方案及审计记录",
      "rationale": "说明该数值的适用性"
    }
  }
}
```

键为 `kind/record_id`；每个 `not_applied` PBE 行都须有决议，不能只补 Cl₂。
也可用 `eV/formula_unit`；元素分子/固体的每原子修正按其真实
`formula_units/atom_count` 换算，不能把每个 Cl₂ 分子的数值直接当成每个 Cl
原子的数值。已在原生索引明确 `provided`/`zero_justified` 的行**禁止**
再给决议，以免重复修正。修正决议属于后处理，不是给 VASP 计算 skill
增加替换 PBE 总能或自动套经验参数的权限。

修正待审的 `INCOMPLETE` 索引还要求每条实际使用的 PBE 决议增加
`"decision_status": "APPROVED"`；决议文件顶层还须有精确源文件的
`native_index_sha256`。每个逐行决议须含完整 `identity_binding` 和
`energy_source` 对象，身份字段与该行完全一致，能量来源的 `ref`、
`locator`、`binding_basis` 和完整身份绑定须与原生 `numerical_audit.energy_source`
相同。这样同一 `record_id` 的旧 attempt 决议不能重放到新 attempt 或另一份索引。
`native_status` 必须精确等于原生字段状态。
原生修正对象缺失时可用 `native_status: "missing"`，显式未应用时用
`not_applied`。缺少批准状态、行缺项、多余行或重复 JSON 键都会失败。
原有 READY 路径继续兼容不含 `decision_status` 的既有决议格式。
原生 PBE 策略 ID 全部或部分缺失时，外部决议须明确给出统一
`policy_id`，不得与任何已知原生 ID 冲突；全部缺失还须逐条给出与源索引
字节哈希和 task/attempt 绑定的已批准修正决议。若从待审草稿沿用
`policy_review`，必须将其标为 `APPROVED`，填写可核的 `source` 与
`rationale`，并保留与原生索引一致的 `known_native_policy_ids`；
仅填写 `policy_id`、保留 `PENDING` 的草稿不能用于正式作图。
直接编制、未携带 `policy_review` 的已批准且逐行绑定的决议继续兼容。

纯净及缺陷 HSE06 方法在 READY 和修正待审入口均须有逐元素的
POTCAR（赝势数据文件）标题与 SHA-256 身份。同种元素的赝势须一致，
而因缺陷改变组成导致的整份 POTCAR 哈希不同不构成错误；同一参考集合
的家族策略、泛函、自旋轨道耦合、共同超胞和 k 点设置仍须兼容。

## 缺陷态修正与带边对齐

缺陷态的外部修正必须使用独立 JSON 文件和 CLI（命令行接口）参数
`--defect-correction-resolutions`，不能混入 PBE 相修正。文件用源原生索引的
**原始字节 SHA-256** 绑定，并以完整 identity（身份）、`state_id` 和整数 `q`
唯一定位每态。`identity_binding` 必须包含 `project_id`、`run_id`、`cluster`、
`task_id`、`package_task_id`、`step_id`、`attempt_id`、`physical_attempt_id`、
`selector_attempt` 和 `job_id`。决议的核心字段如下：

```json
{
  "schema_version": 1,
  "native_index_sha256": "源 calculation-results.json 的 64 位小写 SHA-256",
  "decisions": [{
    "identity_binding": {
      "project_id": "project-id",
      "run_id": "run-id",
      "cluster": "cluster-id",
      "task_id": "task-id",
      "package_task_id": "package-task-id",
      "step_id": "step-id",
      "attempt_id": "attempt-id",
      "physical_attempt_id": "physical-attempt-id",
      "selector_attempt": "root",
      "job_id": "job-id"
    },
    "state_id": "defect-state-q-minus-one",
    "q": -1,
    "native_status": "not_applied",
    "decision_status": "APPROVED",
    "correction_status": "applied",
    "value_ev": 0.12,
    "unit": "eV",
    "model": "经审核的有限尺寸修正模型",
    "source": "模型实现、文献或本地审计产物",
    "rationale": "说明数值、符号和适用性",
    "alignment": {
      "status": "APPROVED",
      "treatment": "included_in_model",
      "reference_id": "与带边同参考的稳定标识",
      "model": "经审核的有限尺寸修正模型",
      "source": "实际势能/对齐来源及审计产物",
      "rationale": "明确该模型已覆盖 potential alignment（势能对齐）",
      "evidence": {
        "pristine_energy_source": {
          "ref": "原生纯净超胞 EIGENVAL 的 SHA-256",
          "locator": "原生带边来源定位"
        },
        "defect_energy_source": {
          "ref": "原生缺陷能量来源的 SHA-256",
          "locator": "原生缺陷能量来源定位"
        },
        "pristine_potential": {
          "source_type": "LOCPOT",
          "source_ref": {
            "ref": "纯净超胞原始电势来源的 SHA-256",
            "locator": "LOCPOT 或逐站点电势定位",
            "binding_basis": "target_probe_selected_result_and_exact_identity",
            "identity_binding": {"task_id": "完整纯净超胞任务身份"}
          },
          "structure": {
            "supercell_id": "与原生纯净超胞一致",
            "lattice_angstrom": [[1, 0, 0], [0, 1, 0], [0, 0, 1]]
          },
          "parameters": {"potential_tag": "方法实际使用的参数"}
        },
        "defect_potential": {
          "source_type": "LOCPOT",
          "source_ref": {
            "ref": "缺陷超胞原始电势来源的 SHA-256",
            "locator": "LOCPOT 或逐站点电势定位",
            "binding_basis": "target_probe_selected_result_and_exact_identity",
            "identity_binding": {"task_id": "完整缺陷态任务身份"}
          },
          "structure": {
            "supercell_id": "与原生缺陷超胞一致",
            "lattice_angstrom": [[1, 0, 0], [0, 1, 0], [0, 0, 1]]
          },
          "parameters": {"potential_tag": "方法实际使用的参数"}
        },
        "audit_source_ref": {
          "sha256": "独立对齐审计产物的 64 位小写 SHA-256",
          "locator": "审计产物中的记录定位",
          "binding_basis": "independent_alignment_audit_bound_to_native_sources",
          "potential_source_sha256": {
            "pristine": "纯净超胞原始电势 SHA-256",
            "defect": "缺陷超胞原始电势 SHA-256"
          },
          "identity_bindings": {
            "pristine": {"project_id": "完整纯净超胞身份字段"},
            "defect": {"project_id": "完整缺陷态身份字段"}
          }
        }
      }
    }
  }]
}
```

示例中的所有 identity 对象都须展开为上文列出的全部身份字段。纯净与缺陷
能量来源的 `ref`/`locator` 必须匹配原生索引；两份原始电势来源还须分别以
SHA-256、定位点和完整 task/attempt 身份绑定到其结构，**同时与原生索引
`data.locpot` 中已验收的选中尝试 LOCPOT 的路径、SHA-256、大小和身份逐项一致**，
并记录所用势能类型和非空方法参数。选用逐站点电势时，同样需要原生索引
`data.potential_sources` 提供对应类型的已验收来源；当前计算端默认仅导出
LOCPOT 元数据，因此不能用决议文件自填的逐站点电势来源绕过原生证据。
两侧超胞 ID 与原生记录一致，且晶格须互相可比；独立审计来源
必须列出并绑定两侧原始电势 SHA-256 及 task/attempt。`treatment` 可为
`included_in_model`，此时
`value_ev` 已包含势能对齐，不得再提供独立对齐数值；也可为
`separate_same_reference`，此时必须在 `alignment` 中提供经审核的
`value_ev`、`unit: "eV"`、`model`、`source` 和 `rationale`，并与缺陷修正相加。
两种形式均须明确 `reference_id`。错哈希、错 attempt/state/q、重复/多余/缺失
决议、非 `APPROVED` 状态或来源不匹配都会失败。
`correction_status` 可为 `applied` 或零状态 `zero_justified` /
`justified_zero`；零状态必须对应数值零，输出索引统一使用组装器接受的
`justified_zero`。

原生索引若已经 `READY` 且包含 `q != 0` 的缺陷态，**原生 `E_corr`
数值本身仍不足以证明已做同参考电势对齐**。另交一份同样绑定原生字节
SHA-256、精确任务/attempt/state/q 的缺陷决议文件；该态的决议须声明
`purpose: "alignment_only"`、`correction_status: "native_value_unchanged"`，
并以 `value_ev`、`unit: "eV"`、`model` 和 `native_status` 逐项核对原生
`E_corr`，`decision_status: "APPROVED"`。其 `alignment` 使用上面的双侧原始
电势/能量及独立对齐审计契约：`included_in_model` 不重复加数值，
`separate_same_reference` 才将另行审定的对齐值加入。不能通过另写一个
E_corr 数值来替换原生值；缺对齐决议的 READY 带电态失败关闭。没有带电态的
原有 READY 索引不需要增加此文件。带电 READY 若附此文件，PBE 决议也须
增加与修正待审分支相同的原生索引 SHA-256、逐行 `APPROVED`、完整
task/attempt 身份及原生能量来源绑定，不能以旧参考尝试的相修正重放。
原生索引若为 `INCOMPLETE`、缺项只在 PBE，但某个带电态已有原生
`E_corr`，也可对此态用相同的 `alignment_only` 决议；不能因为顶层
未就绪就把已有数值替换为另一份缺陷修正。

纯净超胞 `band_edges.energy_zero_provenance.status=DERIVED` 即使绑定了真实
EIGENVAL，也只证明纯净超胞 VBM（价带顶）来源；若
`defect_alignment_status=NOT_COMPUTED`，它不代表带电缺陷对齐已完成。修正待审
输入仍须提供身份绑定的真实纯净带边来源。对 `q != 0`，必须再有明确覆盖势能
对齐的批准模型，或另行审计的同参考对齐证据；一个 `E_corr` 数值本身不能补齐
能量零点/对齐缺口。对 `q = 0`，不要求虚构 `q·VBM` 对齐，决议中应省略
`alignment`，但纯净带边来源仍须有效。
适配器只校验决议字段、哈希与来源绑定，不计算修正、不读取 LOCPOT
（局域静电势文件）或逐站点电势文件，也不替研究者审核模型是否科学适用。
真实输入仍须由研究者提供批准决议及原始能量/电势和审计来源；缺失时拒绝输出。

## 一次运行三个条件

```bash
python3 scripts/run_from_calculation_results.py \
  --native-index "/absolute/calculation-results.json" \
  --correction-resolutions "/absolute/correction-resolutions.json" \
  --defect-correction-resolutions "/absolute/defect-correction-resolutions.json" \
  --output-dir "/absolute/new-output"
```

修正待审索引中只要有缺陷态的原生 `E_corr` 缺失，或原生索引虽为 READY
但有带电态需要独立对齐审计，就须指定
`--defect-correction-resolutions`；仅 PBE 修正待审、没有待审缺陷态时可省略。
已有中性 READY 输入路径及其 PBE 决议格式保持兼容。

可用 `--control-element Cl` 指定化学势控制元素，
`--separate-defect-id <已审计ID>` 给复合缺陷单独成图。默认由相图程序按
卤素/氧规则取三个内部点；该命令分别写出相图分析、每个点的形成能输入和
独立图，三个点不触发三组 VASP 计算。输出目录必须不存在。

若原生索引同时含本征、替位掺杂及复合缺陷态，需要在**每个化学势点**将
所有缺陷画在一张图时，添加 `--all-in-one`。每个 `defect_id` 仍只与自身
的价态/几何取最低包络；高能缺陷过多时可不启用此选项，沿用默认每图至多
六种缺陷的拆图。`--all-in-one` 与 `--separate-defect-id` 不可同时使用，
输入中缺少合格掺杂/复合态时不会自动补造曲线。

掺杂元素通常选一个常见**二元相上界**，例如计算端索引含
`competitive_phase_pbe` 的 `record_id=sbcl3`、式量 `SbCl3`，
且含 PBE/HSE06 Sb 单质参考时，使用
`--dopant-binary-phase-id sbcl3`。该相不挤进宿主稳定域的竞争相约束，
而是在每个富/中/贫点按已修正的 PBE 形成能和该点的 Cl 化学势算
`Δμ_Sb = min[0,(ΔHf(SbCl3)−3Δμ_Cl)]`；一般式按二元相的真实化学计量除以
掺杂元素系数。输出分别记录二元相来源、未截断上限、是否改由 Sb 单质上界
主导及约束余量。这个定义是一种已声明的储库边界选择，不等于已知实验掺杂
浓度；二元相 PBE 修正和 HSE06 Sb 单质参考缺一不可。

若不使用逐点二元相上界，也可使用 `--external-fixed-mu` 传入
带 PBE 参考集合与单位声明的**显式固定**相对化学势 JSON，例如
`{"Sb":{"value":-1.0,"unit":"eV","meaning":"delta_mu","elemental_reference":"PBE参考集合ID"}}`。
这一模式只是研究者选定的外元素化学势条件；同一掺杂元素不得同时用固定值
和二元相逐点上界。即使没有外元素竞争相，固定外元素 Δμ 仍须有同口径的
PBE 外元素单质参考以及形成能所需的 HSE06 单质参考；固定值不能替代
单质计算。没有外来相时两者均省略。

## 失败范围

任一任务未验收、索引行不完整、原生未就绪原因超出允许的修正缺项、修正状态
不明、PBE 元素参考缺失、带边来源/对齐未审计、带边未收敛、缺少 HSE06 元素
参考、缺陷电荷修正缺数值或错误混用参考时拒绝输出正式图。
单命令入口在创建输出目录前检查所有选中点的形成能数值、元素参考和纵轴
显示界限；输入本身无效时不会留下相图半成品。输出目录新建后若磁盘或
绘图库发生异常，仍可能留下部分产物；只有相图分析和
`formation/formation-energy-bundle.json` 均存在且审计通过才视作本次运行完成。
