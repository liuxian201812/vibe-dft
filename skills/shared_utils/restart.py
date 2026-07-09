"""Auto-restart loop templates for VASP submit scripts.

Max 2 loops per project decision (CONTCAR -> POSCAR restart).
"""


def restart_loop_bash(vasp_cmd: str, max_loops: int = 2) -> str:
    """Return bash snippet for auto-restart loop with max_loops iterations.

    Args:
        vasp_cmd: Full VASP command line (e.g. "mpirun vasp_std").
        max_loops: Maximum restart iterations (default 2).

    Returns:
        Bash code string to be inserted into submit.sh.
    """
    return f"""\
# Auto-restart loop (max {max_loops})
MAX_LOOPS={max_loops}
CONVERGENCE="reached required accuracy"
IS_CONVERGED=false
for i in $(seq 1 $MAX_LOOPS); do
  {vasp_cmd} > vasp_loop.out
  if grep -q "$CONVERGENCE" vasp_loop.out; then
    IS_CONVERGED=true
    break
  fi
  if [ -s CONTCAR ]; then
    mv POSCAR "POSCAR_$(date +%Y%m%d_%H%M%S)_loop$i"
    cp CONTCAR POSCAR
  else
    echo "ERROR: CONTCAR missing, cannot restart" >&2
    exit 1
  fi
done
if [ "$IS_CONVERGED" = false ]; then
  echo "ERROR: not converged after $MAX_LOOPS attempts" >&2
  exit 1
fi
# Backup original POSCAR, then keep CONTCAR as POSCAR (use cp to preserve CONTCAR)
mv POSCAR "POSCAR_$(date +%Y%m%d_%H%M%S)_final"
cp CONTCAR POSCAR
"""


def contcar_restart_bash(vasp_cmd: str) -> str:
    """Return bash snippet for one-shot CONTCAR->POSCAR restart, then run.

    Useful for HSE06 single-point or NSCF runs that just need to copy
    relaxed structure before running (no convergence loop).
    """
    return f"""\
# Copy relaxed structure if available
if [ -s CONTCAR ] && [ CONTCAR -nt POSCAR ]; then
  mv POSCAR POSCAR"$SLURM_JOB_ID"
  cp CONTCAR POSCAR
fi
{vasp_cmd}
"""
