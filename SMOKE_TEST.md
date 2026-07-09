# SMOKE TEST

This document verifies that `vibedft` is usable before you trust it for a real calculation.

## Goal

Check these paths end-to-end:

1. environment variables
2. POSCAR validation
3. relax input generation
4. band input generation
5. submit script generation
6. CCD image generation
7. CCD plotting

The smoke test is local. It uses a fake POTCAR library and fake CCD `OSZICAR` files, so it does not consume cluster time.

## Step 1. Create Test Inputs

```bash
python3 smoke_tests/make_test_inputs.py
```

This writes files under `smoke_tests/tmp/`:

- `POSCAR_host`
- `POSCAR_gs`
- `POSCAR_es`
- `fake_potcars/`
- `relax_config.json`
- `band_config.json`
- `ccd_config.json`

## Step 2. Export Temporary Environment Variables

```bash
export VIBEDFT_VASP_CMD='printf "Simulated %s\n" "{vasp_bin}"'
export VIBEDFT_POTCAR_DIR="$PWD/smoke_tests/tmp/fake_potcars"
```

## Step 3. Check Environment Variables

```bash
python3 skills/runtime-config/scripts/check_env.py
```

Expected: all variables show `set`.

## Step 4. Validate POSCAR

```bash
python3 skills/poscar-generation/scripts/validate_poscar.py smoke_tests/tmp/POSCAR_host
```

Expected: `Verdict: PASS`.

## Step 5. Generate Relax Input

```bash
python3 skills/vasp-input/scripts/vasp_input_suite.py --phase relax --config smoke_tests/tmp/relax_config.json
```

Expected output directory:

```text
smoke_tests/tmp/relax_case/
├── 01_relax/
├── tasks/01_relax/
├── run_01_relax.sh
└── submit.sh
```

Check:

```bash
ls smoke_tests/tmp/relax_case/tasks/01_relax
```

Expected files:

- `POSCAR`
- `INCAR`
- `KPOINTS`
- `POTCAR`
- `submit.sh`

## Step 6. Generate Band Input

```bash
python3 skills/vasp-input/scripts/vasp_input_suite.py --phase band --config smoke_tests/tmp/band_config.json
```

Expected output directory:

```text
smoke_tests/tmp/band_case/band/
```

Expected files:

- `POSCAR`
- `INCAR`
- `KPOINTS`
- `POTCAR`
- `submit.sh`

## Step 7. Generate Submit Scripts

```bash
python3 skills/vasp-submit/scripts/vasp_submit_suite.py --config smoke_tests/tmp/relax_config.json
```

Expected:

- `run_01_relax.sh`
- `submit.sh`

## Step 8. Generate CCD Images

```bash
python3 skills/ccd/scripts/ccd_generate.py smoke_tests/tmp/ccd_config.json
```

Expected output directory:

```text
smoke_tests/tmp/ccd_case/ccd/
├── ccd_manifest.json
├── image_000/
├── image_001/
├── image_002/
├── image_003/
└── image_004/
```

## Step 9. Add Fake CCD Energies

The plotter expects real `OSZICAR` files. For a local smoke test, write fake ones:

```bash
for d in smoke_tests/tmp/ccd_case/ccd/image_*; do
  n=$(basename "$d" | sed 's/image_//')
  e=$(awk -v n="$n" 'BEGIN { printf "%.6f", -10.0 + (n + 0) * 0.05 }')
  printf ' 1 F=  %s E0=  %s  d E =0.000000\n' "$e" "$e" > "$d/OSZICAR"
done
```

## Step 10. Plot CCD

```bash
python3 skills/ccd/scripts/ccd_plot.py smoke_tests/tmp/ccd_case/ccd --output smoke_tests/tmp/ccd_case/ccd.png
```

Expected:

- terminal prints `Wrote smoke_tests/tmp/ccd_case/ccd.png`
- output PNG exists

## Pass Criteria

The smoke test passes if all of the following are true:

1. `check_env.py` shows `set`
2. POSCAR validation passes
3. relax input files are generated
4. band input files are generated
5. submit scripts are generated
6. CCD image directories are generated
7. CCD plot is written successfully

## After the Smoke Test

Delete the generated local test files if you do not want to keep them:

```bash
rm -rf smoke_tests/tmp
```
