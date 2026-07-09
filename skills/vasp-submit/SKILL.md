---
name: vasp-submit
description: Generate independent SLURM submission scripts for VASP jobs with explicit user-provided scheduler settings. Use after vasp-input has generated inputs.
---

# VASP Submit Basic

## Scope

This skill generates one `run_<step>.sh` per step and one `submit.sh` that submits them independently.

It intentionally does not generate dependency chains, `afterok`, daemon monitoring, queue fallback, or task status export.

## Usage

```bash
python3 scripts/vasp_submit_suite.py --config task_config.json
cd <workdir> && bash submit.sh
```

The config can be a simple `steps` config or the same `tasks` config used by `vasp-input`.

## Basic Steps Config

```json
{
  "workdir": "work/manual_submit",
  "runtime": "local-profile",
  "partition": "<scheduler-partition>",
  "ncpus": 32,
  "steps": [
    {
      "name": "relax",
      "inputs": {
        "poscar": "calc/POSCAR",
        "incar": "calc/INCAR",
        "potcar": "calc/POTCAR",
        "kpoints": "calc/KPOINTS"
      },
      "vasp_bin": "vasp_std"
    }
  ]
}
```
