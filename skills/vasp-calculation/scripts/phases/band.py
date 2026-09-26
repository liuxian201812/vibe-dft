"""Phase: band (能带计算 NSCF)

Migrated from band-structure/scripts/band_pipeline.py.
Uses common.make_incar() with ST template instead of hardcoded Incar(dict).
"""

from __future__ import annotations

import json
import subprocess
import sys
from pathlib import Path

from pymatgen.core import Structure
from pymatgen.io.vasp import Kpoints, Poscar

_SCRIPT_DIR = Path(__file__).resolve().parent
_PARENT_SCRIPT = _SCRIPT_DIR.parent
if str(_PARENT_SCRIPT) not in sys.path:
    sys.path.insert(0, str(_PARENT_SCRIPT))

from common import (
    make_incar, get_encut, write_potcar,
    write_submit_sh, replace_tag,
)


def _band_ncore(server_label: str, ncpus: int, natoms: int | None = None) -> int:
    try:
        from shared_utils.server import default_ncore
        return default_ncore(ncpus, natoms=natoms)
    except ImportError:
        return max(ncpus // 4, 1)


def _generate_band_path(poscar_path: str) -> dict:
    """Generate explicit k-path via seekpath.

    seekpath expects the spglib tuple `(cell, positions, numbers)` rather than a
    pymatgen `Structure` instance directly.
    """
    import seekpath

    struct = Structure.from_file(poscar_path)
    seekpath_input = (
        struct.lattice.matrix,
        [site.frac_coords for site in struct],
        [site.specie.number for site in struct],
    )
    return seekpath.get_explicit_k_path(
        seekpath_input, reference_distance=0.025, symprec=1e-5
    )


def _write_band_kpoints(task_dir: Path, poscar_path: str):
    """Write band-path KPOINTS."""
    result = _generate_band_path(poscar_path)
    kpts = result["explicit_kpoints_rel"]
    labels = result["explicit_kpoints_labels"]
    lines = [str(len(kpts)), "Reciprocal"]
    for i in range(0, len(kpts), 2):
        if i + 1 < len(kpts):
            lines.append(
                f"  {kpts[i][0]:.10f}  {kpts[i][1]:.10f}  {kpts[i][2]:.10f}  !{labels[i]}"
            )
            lines.append(
                f"  {kpts[i+1][0]:.10f}  {kpts[i+1][1]:.10f}  {kpts[i+1][2]:.10f}  !{labels[i+1]}"
            )
            lines.append("")
    (task_dir / "KPOINTS").write_text("\n".join(lines))


def generate(config: dict, workdir: Path) -> dict:
    """Generate band structure VASP inputs.

    Config:
        poscar: path to relaxed CONTCAR/POSCAR
        functional: PBE / R2SCAN / HSE06
        encut: optional (auto from POTCAR if omitted)
        ncpus, partition, qos: SLURM params
    """
    poscar_path = config["poscar"]
    functional = config.get("functional", "R2SCAN")
    _server_label = config.get("runtime", config.get("server", ""))
    ncpus = int(config.get("ncpus", 1))
    partition = config.get("partition", "")
    qos = config.get("qos", "")
    vasp_cmd = config.get("vasp_cmd", "")
    if not partition:
        raise ValueError("partition is required for band submit script generation")

    task_dir = workdir / "band"
    task_dir.mkdir(parents=True, exist_ok=True)

    # POSCAR
    struct = Structure.from_file(poscar_path)
    Poscar(struct).write_file(str(task_dir / "POSCAR"))

    # KPOINTS (band path)
    _write_band_kpoints(task_dir, str(task_dir / "POSCAR"))

    # POTCAR
    try:
        write_potcar(task_dir)
    except Exception as e:
        print(f"  POTCAR skip ({e})")

    encut = config.get("encut")
    if not encut and (task_dir / "POTCAR").exists():
        encut = get_encut(str(task_dir / "POTCAR"))
    encut = encut or 400

    # INCAR: ST template + functional + ICHARG=11 + LORBIT=10
    incar = make_incar(
        template_key="ST",
        functional=functional,
        encut=encut,
        ncore=_band_ncore(_server_label, ncpus, natoms=len(struct)),
    )
    incar = replace_tag(incar, "ICHARG", "11")
    if "LORBIT" not in incar:
        incar += "\nLORBIT = 10"
    (task_dir / "INCAR").write_text(incar)

    # submit.sh
    write_submit_sh(
        workdir, "band",
        vasp_bin="vasp_std", ncores=ncpus, partition=partition,
        qos=qos, vasp_cmd=vasp_cmd,
    )

    manifest = {
        "phase": "band",
        "functional": functional,
        "encut": encut,
        "workdir": str(task_dir),
        "steps": [{"name": "band", "dir": str(task_dir), "depends_on": []}],
    }
    (task_dir / "band_manifest.json").write_text(json.dumps(manifest, indent=2))
    print(f"  band: OK ({functional}, encut={encut})")
    return manifest
