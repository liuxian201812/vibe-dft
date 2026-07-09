# QUICKSTART

This document is for users who want to try `vibedft` without embedding any
site-specific runtime facts in the repository.

## 1. Clone

```bash
git clone https://github.com/liuxian201812/vibe-dft.git
cd vibe-dft
```

## 2. Install Python Dependencies

```bash
python3 -m venv .venv
.venv/bin/pip install -r requirements.txt
```

This installs the local dependencies needed to run the public scripts in this
repository. It does not describe or manage any remote runtime environment.

## 3. Configure Environment Variables

Copy relevant lines from `.env.example` into your shell profile, then replace
placeholders with your own runtime values.

```bash
cp .env.example /tmp/vibedft.env.example
```

Typical workflow:

1. Open `.env.example`
2. Copy the generic section into `~/.bashrc`
3. Replace every `<...>` placeholder with your own real value
4. `source ~/.bashrc`

Check:

```bash
python3 skills/runtime-config/scripts/check_env.py
```

Expected result: all required variables show `set`.

If you do not need remote submit-script generation, you can still use the repo
for POSCAR validation, input-file generation, and local post-processing.

## 4. Run Smoke Test First

```bash
python3 smoke_tests/make_test_inputs.py
```

Then follow `SMOKE_TEST.md`.

## 5. First Real Workflow

### Search and download structure

```bash
python3 skills/poscar-generation/scripts/search_cod.py search "Cs2NaInCl6"
python3 skills/poscar-generation/scripts/search_cod.py download <cod_id> -o host.cif
```

### Convert CIF to POSCAR

```python
from pymatgen.core import Structure
Structure.from_file("host.cif").to(filename="POSCAR", fmt="poscar")
```

### Validate POSCAR

```bash
python3 skills/poscar-generation/scripts/validate_poscar.py POSCAR
```

### Generate relax input

Prepare `relax_config.json`:

```json
{
  "workdir": "work/example",
  "runtime": "local-profile",
  "partition": "<scheduler-partition>",
  "qos": "<optional-qos>",
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
      "isif": 3,
      "kpoints_density": 0.04
    }
  ]
}
```

Generate input files:

```bash
python3 skills/vasp-input/scripts/vasp_input_suite.py --phase relax --config relax_config.json
```

### Generate submit scripts

```bash
python3 skills/vasp-submit/scripts/vasp_submit_suite.py --config relax_config.json
```

### Submit

```bash
cd work/example
bash submit.sh
```

## 6. CCD Workflow

If you already have two relaxed endpoint structures:

```bash
python3 skills/ccd/scripts/ccd_generate.py ccd_config.json
for d in work/ccd_example/ccd/image_*; do (cd "$d" && sbatch submit.sh); done
python3 skills/ccd/scripts/ccd_plot.py work/ccd_example/ccd --output ccd.png
```

## Notes

- Submission is independent-only. No `afterok` dependency chain is generated.
- Submit scripts depend on explicit scheduler settings in your config JSON.
- POTCAR generation depends on the configured POTCAR library path.
- It is for stable standard workflows, not every advanced VASP mode.
