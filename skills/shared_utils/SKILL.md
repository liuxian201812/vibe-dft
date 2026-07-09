---
name: shared-utils
description: Shared utility library for VASP workflow scripts: POTCAR generation/ENCUT calculation, ENCUT three-tier rule, DFT+U auto-detection, restart-loop bash templates, and generic runtime command retrieval. Imported by other skills in this repository.
---

# Shared Utils Library

## Modules

| Module | Exports | Used By |
|--------|---------|---------|
| `server.py` | `get_vasp_cmd()`, `default_ncore()` | vasp-input, vasp-submit |
| `potcar.py` | `generate_potcar()`, `calc_encut()`, `get_potcar_lib_dir()` | vasp-input |
| `restart.py` | `restart_loop_bash()`, `contcar_restart_bash()` | vasp-input common.py |
| `lda_u.py` | `auto_lda_u()`, `merge_u_config()` | vasp-input relax phase |

## Environment Variables

- `VIBEDFT_VASP_CMD` — VASP execution command template with `{vasp_bin}` placeholder
- `VIBEDFT_POTCAR_DIR` — POTCAR library directory for pure-Python fallback

## ENCUT Rule

```
raw = max(ENMAX) * 1.5
tier in [400, 520, 680]; capped at 680
```

## DFT+U Auto-Detection

- Anion detected by electronegativity
- U tables: MP (primary) + MIT (supplement)
- Cl/Br/I/N fallback to O table
- LMAXMIX: 4 for d-electrons, 6 for f-electrons
