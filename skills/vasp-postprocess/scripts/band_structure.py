#!/usr/bin/env python3
"""Band structure post-processing: extract band gap and band data from VASP output.

Usage:
  python3 post_band.py /path/to/task_dir
"""

from __future__ import annotations

import json, sys
from pathlib import Path


def extract_band_gap(vasprun_path: str | Path) -> dict:
    """Extract band gap, VBM, CBM from a converged vasprun.xml.

    Args:
        vasprun_path: Path to vasprun.xml (NSCF with band structure)

    Returns:
        Dict with keys: band_gap, vbm, cbm, is_direct, fermi_level
    """
    from pymatgen.io.vasp.outputs import Vasprun
    v = Vasprun(str(vasprun_path))
    if not v.converged:
        raise ValueError(f"Calculation not converged: {vasprun_path}")
    props = v.eigenvalue_band_properties
    return {
        "band_gap": props[0],
        "vbm": float(props[1]),
        "cbm": float(props[2]),
        "is_direct": bool(props[3]),
        "fermi_level": v.efermi if hasattr(v, 'efermi') else 0.0,
    }


def extract_band_structure(vasprun_path: str | Path) -> dict:
    """Extract full band structure (k-points, eigenvalues) from vasprun.xml.

    Args:
        vasprun_path: Path to vasprun.xml

    Returns:
        Dict with keys: kpoints (list of 3-tuples), labels (list),
        bands (list of bands: list of eigenvalues per kpoint), nbands, nkpoints
    """
    from pymatgen.io.vasp.outputs import Vasprun
    from pymatgen.electronic_structure.bandstructure import get_reconstructed_band_structure
    v = Vasprun(str(vasprun_path))
    bs = v.get_band_structure(kpoints_filename=None)
    return {
        "nkpoints": len(bs.kpoints),
        "nbands": bs.nb_bands,
        "efermi": bs.efermi,
        "is_metal": bs.is_metal(),
        "direct_gap": bs.get_direct_band_gap(),
    }


def call_vaspkit_211(task_dir: str | Path) -> str:
    """Run vaspkit -task 211 to extract BAND.dat and KLINES.dat.

    Returns:
        Path to BAND.dat
    """
    import subprocess
    td = Path(task_dir)
    subprocess.run(["vaspkit", "-task", "211"], cwd=str(td), check=True,
                   timeout=60, capture_output=True)
    band_dat = td / "BAND.dat"
    if not band_dat.exists():
        raise FileNotFoundError(f"vaspkit 211 did not produce BAND.dat in {td}")
    return str(band_dat)


def main():
    if len(sys.argv) < 2:
        print(__doc__, file=sys.stderr)
        sys.exit(1)
    td = Path(sys.argv[1])
    if not td.is_dir():
        print(f"ERROR: {td} not found", file=sys.stderr)
        sys.exit(1)
    vasprun = td / "vasprun.xml"
    if not vasprun.exists():
        print(f"ERROR: {vasprun} not found", file=sys.stderr)
        sys.exit(1)
    try:
        gap = extract_band_gap(vasprun)
        print(json.dumps(gap, indent=2))
    except Exception as e:
        print(f"Band gap extraction failed: {e}", file=sys.stderr)
        sys.exit(1)
    try:
        bs = extract_band_structure(vasprun)
        print(f"Bands: {bs['nkpoints']} k-points × {bs['nbands']} bands")
        print(f"Is metal: {bs['is_metal']}")
    except Exception as e:
        print(f"Band structure extraction failed: {e}", file=sys.stderr)


if __name__ == "__main__":
    main()
