#!/bin/bash
set -euo pipefail

ROOT=$(cd "$(dirname "$0")/.." && pwd)
cd "$ROOT"

PYTHON=${PYTHON:-python3}

"$PYTHON" smoke_tests/make_test_inputs.py >/dev/null

export VIBEDFT_VASP_CMD='printf "Simulated %s\n" "{vasp_bin}"'
export VIBEDFT_POTCAR_DIR="$ROOT/smoke_tests/tmp/fake_potcars"

"$PYTHON" skills/runtime-config/scripts/check_env.py
"$PYTHON" skills/poscar-generation/scripts/validate_poscar.py smoke_tests/tmp/POSCAR_host
"$PYTHON" skills/vasp-input/scripts/vasp_input_suite.py --phase relax --config smoke_tests/tmp/relax_config.json
"$PYTHON" skills/vasp-input/scripts/vasp_input_suite.py --phase band --config smoke_tests/tmp/band_config.json
"$PYTHON" skills/vasp-submit/scripts/vasp_submit_suite.py --config smoke_tests/tmp/relax_config.json
"$PYTHON" skills/ccd/scripts/ccd_generate.py smoke_tests/tmp/ccd_config.json

for d in smoke_tests/tmp/ccd_case/ccd/image_*; do
  n=$(basename "$d" | sed 's/image_//')
  e=$(awk -v n="$n" 'BEGIN { printf "%.6f", -10.0 + (n + 0) * 0.05 }')
  printf ' 1 F=  %s E0=  %s  d E =0.000000\n' "$e" "$e" > "$d/OSZICAR"
done

"$PYTHON" skills/ccd/scripts/ccd_plot.py smoke_tests/tmp/ccd_case/ccd --output smoke_tests/tmp/ccd_case/ccd.png

echo "Smoke test passed."
