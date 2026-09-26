---
name: chemical-potential-and-formation-energy
description: "Post-process audited JSON into chemical-potential stability domains, feasible full-vector points, and per-point defect-formation-energy envelopes from accepted structured calculation indices. Never extract from raw VASP outputs or imply experimental accessibility."
---
# 化学势稳定域与缺陷形成能

本 skill 处理已审计并显式整理的数值 JSON（结构化数据文件）。它不扫描计算目录、不从 VASP（第一性原理软件）原始输出提取数值、不计算带电缺陷修正，也不验收计算是否收敛。字段校验通过不代表来源、修正模型或物理可比性已获科学验收。

## 流程

- **化学势稳定域**：使用 `scripts/chemical_potential_phase_diagram.py`。在形成同一能量参考和修正口径的稳定域后，默认控制元素优先为宿主中的 O；没有 O 时，仅在存在一个经 `pymatgen.core.Element(symbol).is_halogen`（材料元素分类接口）识别的卤素时默认选该卤素。用户显式给出的 `point_selection.control_element` 始终优先。默认从该元素贫到富按 1/4、1/2、3/4 插值；分数可改为 `[1/3, 2/3]` 或包含 0、1 的任意唯一 `[0,1]` 值。每个自动点包含全部宿主元素的 $\Delta\mu$，并核验宿主等式、单质上界和全部竞争相约束。
- **无默认控制元素及高维宿主**：多卤素宿主或既无 O 又无唯一卤素的宿主仍正常计算并输出相稳定域；`point_selection.selection_status` 标记 `requires_explicit_control_element`，此时不生成自动代表点。用户必须显式指定控制元素后重跑相图，才能让批处理生成相应形成能输入。显式 `mode=manual` 仍可提供完整 `selected_points`；程序逐点核验全部坐标和约束余量。四元固定一个、五元固定两个元素的二维切片须引用已计算的真实截面；五元完整稳定域为四维，图中仅画明确标注的二维**投影**，不是完整域或切片。每个选点始终保留全部宿主元素坐标。
- **形成能**：单点分析使用 `scripts/defect_formation_energy.py`。批量组装使用 `scripts/build_formation_energy_bundle.py`，输入为已验收的结构化计算索引、相图分析 JSON 和全元素 HSE06（Heyd–Scuseria–Ernzerhof 杂化泛函）单质参考。脚本先校验任务验收、来源引用、能量、带边、化学计量差和逐态 `E_corr`，再按每个已选化学势点写一份符合 formation `validate_input` 的 JSON，并分别生成审计记录和图。每个 `defect_id` 独立取最低包络；费米能横轴为 VBM=0 到纯净超胞 CBM。
- **混合缺陷合图**：同一个已验收索引含本征、替位掺杂、孤立填隙掺杂或复合缺陷时，`run_from_calculation_results.py --all-in-one` 将每个化学势点的全部 `defect_id` 画在同一张图；默认仍按每图至多六个缺陷拆图。孤立填隙须先由计算端以审核后的同一纯净超胞、原生身份和 `n_i` 导出到**这一份**索引；本 skill 不从结构草稿或两份旧索引拼接其能量。`--all-in-one` 与 `--separate-defect-id` 互斥；合图只改变展示，不改变各缺陷独立包络或三点化学势。过密图可继续使用默认拆图。
- **独立修正审查准备**：对可信的原生结果索引，`scripts/prepare_correction_drafts.py` 可只读写出 PBE 和逐态电荷/对齐决议的身份绑定**草稿**；数值、模型、来源与批准状态均需研究者独立审核，草稿保持 `PENDING` 且不能用于正式作图。它只减少手工抄写 task/attempt/来源的错误，不替代气相、介电或电荷修正计算。
- **计算端索引接入**：`calculation-results-1` 整体为 `READY`，或虽为 `INCOMPLETE` 但计算/数值已验收且**仅缺独立修正决议**时，使用 `scripts/run_from_calculation_results.py` 做显式 PBE 修正决议、格式转换、相图分析与逐点形成能图；仅需转换而暂不作图时用 `scripts/adapt_calculation_results.py`。后一分支须以原生索引字节哈希及精确任务/attempt 绑定外部 PBE/逐态修正决议；**已有原生 E_corr 的 READY 带电态**仍须另附只审计对齐、不替换原数值的独立决议。所有带电态均须有与原生索引匹配的双侧原始电势及势能对齐审计，中性态不虚构对齐。它不更改原生 `INCOMPLETE` 状态，也不自动填零。指定常见掺杂二元相的 `record_id` 可对各稳定点求外来元素的二元相上界，且不把稀掺杂相当成宿主相图的新坐标维度。两脚本只消费显式 JSON，不替计算端从原始输出搜集结果。运行前读[原生索引与修正决议](references/native-calculation-results.md)；字段检查不能代替修正方法的科学审定。
- **结构请求**：结构生成不属于本 skill。按结构类型分别转交 `doping-supercell`、`defect-complex-builder` 或 `defect-screening`；本 skill 只消费它们及计算端交付的已审核结构/数值引用。
- **运行与复核**：可从任意工作目录调用脚本，输入和新输出位置由参数显式指定。索引字段与来源关联见[计算索引适配契约](references/calculation-index-adapter.md)，各脚本参数用 `--help` 查询。

富/贫端点和内部插值点只是给定模型稳定域中的几何条件，不默认等同于实验可达的生长条件。缺少经审计的 `E_corr`、共同 `E_host`、VBM/CBM、`n_i`、完整 HSE06 元素参考或任一源任务验收时，批处理失败关闭。脚本不会推断或补齐这些量。
