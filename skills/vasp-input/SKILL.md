---
name: vasp-input
description: Generate VASP inputs for relax, static/scf, band, and CCD calculations with site-specific runtime details supplied out of band. Use when writing INCAR/KPOINTS/POTCAR/POSCAR/submit.sh for standard VASP jobs.
---

# VASP Input Basic

## Supported Phases

- `relax`: structure relaxation or static/scf via templates `SR`, `LR`, `ST`
- `band`: NSCF band-path input using seekpath
- `ccd`: interpolate two existing endpoint structures and write standalone single-point image directories
- `auto`: dispatch task types `relax`, `static`, `scf`, `band`, `ccd`

## Not Included

No ΔSCF, chemical potentials, phonons, NEB, MD, GW/BSE, SOC, D3, optics, or elastic workflows.

## Usage

```bash
python3 scripts/vasp_input_suite.py --phase relax --config relax_config.json
python3 scripts/vasp_input_suite.py --phase band --config band_config.json
python3 scripts/vasp_input_suite.py --phase ccd --config ccd_config.json
```

## Relax Config

```json
{
  "workdir": "work/example",
  "runtime": "local-profile",
  "partition": "<scheduler-partition>",
  "ncpus": 32,
  "host_structures": {"host": "POSCAR"},
  "tasks": [
    {"id": "01", "name": "relax", "type": "relax", "host": "host", "template": "LR", "functional": "PBE", "isif": 3}
  ]
}
```

## CCD Config

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
  "nupdown": 0
}
```

CCD endpoints must already exist and have the same site order.
