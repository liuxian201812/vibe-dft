---
name: vasp-calculation
description: Generate local VASP inputs and independent SLURM submit scripts for relax, static/scf, band, and CCD calculations. Use when preparing VASP inputs or submission files without a private scheduling service.
---

# VASP Calculation

This skill is the public, local-only VASP preparation entry point. It combines
input generation and submit-script assembly in one skill and does not include
remote task registration, queue discovery, private cluster profiles, or
automatic job control.

## Scope

- Input generation for `relax`, `static`, `scf`, `band`, and `ccd`.
- Independent SLURM submit-script generation from the generated task manifest.
- Explicit runtime information supplied by the user or a private local profile.

Not included:

- Automatic server selection or remote status ownership.
- DeltaSCF, chemical-potential, phonon, NEB, MD, GW/BSE, optics, or elastic
  workflow generation. Use the corresponding public post-processing or
  structure skill when available.

## Generate Inputs

```bash
python3 scripts/vasp_calculation.py prepare -- --phase relax --config <RELAX_CONFIG_JSON>
python3 scripts/vasp_calculation.py prepare -- --phase band --config <BAND_CONFIG_JSON>
python3 scripts/vasp_calculation.py prepare -- --phase ccd --config <CCD_CONFIG_JSON>
```

The input generator writes `INCAR`, `KPOINTS`, `POTCAR`, `POSCAR`, and
`submit.sh` under the configured work directory and emits the task manifest
used by the submit-script assembler.

## Generate Submit Scripts

```bash
python3 scripts/vasp_calculation.py submit-scripts -- --config <TASK_CONFIG_JSON>
```

This creates one `run_<step>.sh` per step and one `submit.sh` that submits the
independent steps. It intentionally does not create `afterok` chains, daemon
monitoring, queue fallback, or remote task status export.

## Runtime Configuration

The repository contains no site-specific defaults. Set only the values needed
for the current machine and keep real values outside tracked files:

```bash
export VIBEDFT_VASP_CMD='<VASP_LAUNCH_COMMAND_WITH_{vasp_bin}>'
export VIBEDFT_POTCAR_DIR='<POTCAR_LIBRARY_PATH>'
```

The command template must contain the `{vasp_bin}` placeholder. Scheduler
partition, account, QoS, walltime, and other site-specific values belong in
the task config or in an untracked local profile.

## Input Examples

Relax:

```json
{
  "workdir": "<LOCAL_WORK_DIRECTORY>",
  "runtime": "<LOCAL_RUNTIME_PROFILE>",
  "partition": "<SCHEDULER_PARTITION>",
  "ncpus": 32,
  "host_structures": {"host": "POSCAR"},
  "tasks": [
    {
      "id": "01",
      "name": "relax",
      "type": "relax",
      "host": "host",
      "template": "LR",
      "functional": "PBE",
      "isif": 3
    }
  ]
}
```

CCD:

```json
{
  "workdir": "<LOCAL_CCD_WORK_DIRECTORY>",
  "runtime": "<LOCAL_RUNTIME_PROFILE>",
  "partition": "<SCHEDULER_PARTITION>",
  "ncpus": 32,
  "poscar_gs": "<GROUND_STATE_CONTCAR>",
  "poscar_es": "<EXCITED_STATE_CONTCAR>",
  "functional": "PBE",
  "n_images": 41,
  "nupdown": 0
}
```

CCD endpoints must already exist and use the same site order. Keep any config
containing real scheduler, account, path, or runtime values untracked.
