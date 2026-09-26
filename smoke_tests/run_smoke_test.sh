#!/bin/bash
set -euo pipefail

ROOT=$(cd "$(dirname "$0")/.." && pwd)
cd "$ROOT"

PYTHON=${PYTHON:-python3}

"$PYTHON" smoke_tests/make_test_inputs.py >/dev/null

export VIBEDFT_VASP_CMD='printf "Simulated %s\n" "{vasp_bin}"'
export VIBEDFT_POTCAR_DIR="$ROOT/smoke_tests/tmp/fake_potcars"

"$PYTHON" skills/poscar-generation/scripts/validate_poscar.py smoke_tests/tmp/POSCAR_host \
  --template smoke_tests/tmp/POSCAR_host --task-id smoke-test --source synthetic
"$PYTHON" skills/vasp-calculation/scripts/vasp_calculation.py prepare -- \
  --phase relax --config smoke_tests/tmp/relax_config.json
"$PYTHON" skills/vasp-calculation/scripts/vasp_calculation.py prepare -- \
  --phase band --config smoke_tests/tmp/band_config.json
"$PYTHON" skills/vasp-calculation/scripts/vasp_calculation.py submit-scripts -- \
  --config smoke_tests/tmp/relax_config.json

echo "Smoke test passed."
