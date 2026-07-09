---
name: ccd
description: Generate and post-process VASP configuration-coordinate diagram (CCD) calculations. Use when interpolating two relaxed endpoint structures into CCD single-point image jobs and plotting alpha-energy curves.
---

# CCD Basic

## Scope

This skill handles CCD calculations:

1. Interpolate between two existing endpoint structures.
2. Generate standalone VASP single-point directories for each image.
3. Plot image energies after VASP finishes.

It does not submit dependencies, wait for relax jobs, infer endpoints from ΔSCF workflows, or manage multi-state physics automatically.

## Requirements

- `poscar_gs` and `poscar_es` must already exist.
- Endpoint structures must have the same site count and site species order.
- POTCAR generation requires `VIBEDFT_POTCAR_DIR` or `PMG_VASP_PSP_DIR`.

## Generate CCD Inputs

Example `ccd_config.json`:

```json
{
  "workdir": "work/ccd_example",
  "runtime": "local-profile",
  "partition": "<scheduler-partition>",
  "ncpus": 32,
  "poscar_gs": "gs/CONTCAR",
  "poscar_es": "es/CONTCAR",
  "functional": "PBE",
  "n_images": 41,
  "alpha_start": -0.5,
  "alpha_stop": 1.5,
  "nupdown": 0
}
```

Run:

```bash
python3 scripts/ccd_generate.py ccd_config.json
```

Output:

```text
work/ccd_example/ccd/
├── ccd_manifest.json
├── image_000/
│   ├── POSCAR INCAR KPOINTS POTCAR submit.sh
├── image_001/
└── ...
```

## Submit Images

Each image has its own `submit.sh`. Submit manually or with a simple shell loop:

```bash
for d in work/ccd_example/ccd/image_*; do (cd "$d" && sbatch submit.sh); done
```

## Plot CCD

After VASP finishes:

```bash
python3 scripts/ccd_plot.py work/ccd_example/ccd --output ccd.png
```
