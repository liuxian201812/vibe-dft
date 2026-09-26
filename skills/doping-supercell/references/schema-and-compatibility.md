# Schema 与兼容接口

本文件只记录低频但仍有效的字段契约。当前公开运行时行为以
`scripts/core/dope_supercell.py` 和 `schemas/` 为准。仓库未随附专用
`tests/` 夹具，因此文档不以未发布测试作为契约依据。

## 活动 Schema

| 文件 | 用途 | 状态 |
|---|---|---|
| `schemas/build_supercell.schema.json` | Stage 1 配置 | 活动 |
| `schemas/supercell_manifest.schema.json` | Stage 1 manifest | 活动，v1 |
| `schemas/dope.schema.json` | Stage 2 配置 | 活动 |
| `schemas/doping_provenance_v2.schema.json` | Stage 2 provenance | 活动，v2 |
| `schemas/doping_provenance.schema.json` | 旧 provenance | 仅兼容，不生成新结果 |

Schema 文件是运行资产，不要移入普通参考目录或删除。

## 字段别名

运行时对可迁移配置接受下列别名；同一逻辑字段出现多个不同值时直接失败。

| 逻辑字段 | 接受的字段 |
|---|---|
| Stage 1 parent POSCAR | `parent_poscar`、`poscar` |
| Stage 1 manifest 输出 | `manifest`、`manifest_json` |
| Stage 2 manifest 输入 | `supercell_manifest`、`manifest`、`manifest_json` |
| Stage 2 pristine POSCAR | `supercell_poscar`、`supercell`、`poscar` |
| Stage 2 substitution | `substitution` 或单元素的 `substitutions` |
| Stage 2 doped 输出 | `provenance`、`report`、`report_json` |
| parent site index | `parent_site.index_1based`、`index_0based`，或对应的 `parent_index_*` 别名 |
| target site index | `substitution.site_*`、`target_site_*`，或 `target_site.supercell_index_*` |
| target selection | `substitution.site_selection`、顶层 `site_selection`、`target_site.selection` |
| parent fractional coordinate | `fractional_coordinate`、`fractional_coords` |
| output checks | `checks.*`；兼容 `validation.*` 和顶层 `expected_*`/`required_*` |

规范新配置优先使用主入口文档中的字段，不要为了兼容而同时填写别名。

## 严格 Stage 2

- 顶层只能有运行时允许的字段；`target_site`、`parent_site` 和
  `substitution` 子对象也有独立的未知字段检查。
- `mode`、`supercell_matrix`、`min_atoms` 和 `min_image_distance` 即使位于
  兼容白名单中，出现在 Stage 2 配置仍会被拒绝；Stage 2 不重新生成超胞。
- `substitution` 与 `substitutions` 不能同时出现，后者必须恰好包含一个对象。
- 精确站点只能提供一种索引基准；布尔值、负值、非整数和越界值失败。
- `site_selection` 只能是 `nearest_supercell_center`，且不能与精确站点合用。
- 选择自动模式时，provenance 必须含全部匹配平移等价候选、固定
  `1e-10 Å` 并列容差和最低 pristine 站点编号的确定性决策。
- 选择精确模式时，provenance 必须标记
  `candidate_enumeration_performed: false`，不能因为其他候选异常而改变旧行为。

## 旧 Python 接口

以下名称继续保留给旧调用者：

- `generate_supercell(...)`
- `apply_substitutions(...)`
- `map_supercell_site_to_parent(...)`
- `dope_supercell(config_path)`

`generate_supercell(..., allow_primitive=True)` 仍可用于旧的 primitive
比较路径；可审计的 Stage 1 builder 固定使用 `allow_primitive=False`。
`apply_substitutions()` 仍接受零基 `sites`、`count` 和 `all` 形式，但这些
形式不属于严格 Stage 2。`dope_supercell()` 现在只执行 Stage 2，不恢复旧的
“扩展并掺杂”组合动作。

## Provenance 最低字段

新输出必须声称
`doping-supercell.doping-provenance/v2`，并保留：

- source manifest、parent POSCAR、pristine POSCAR 的路径和 SHA256；
- `supercell_matrix`、行向量映射公式和 determinant；
- `from`/`to`、输入或解析出的 supercell 站点、输出站点；
- parent 站点编号、物种、分数坐标、Wyckoff 字符串和 mapping distance；
- 自动选择的几何中心、直接 Cartesian 距离、周期性诊断距离及候选；
- 最终 formula、atom count、species-block order 和 doped POSCAR SHA256。

`doping_provenance.schema.json` 仅用于识别旧 v1 数据；不要让新结果伪装成
v1，也不要用 v1 校验器替代 v2 校验。
