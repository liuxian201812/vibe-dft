"""POTCAR and ENCUT utilities.

POTCAR generation:
  1. Try pymatgen ``Potcar``.
  2. Fallback: pure-Python concatenation from a user-provided POTCAR library.

ENCUT rule: max_enmax * 1.5, round up to nearest tier in [400, 520, 680].
"""
from __future__ import annotations

import os
import re
from pathlib import Path


# Recommended POTCAR symbol for each element (PBE 5.4), from MPRelaxSet.CONFIG.
# Matches vaspkit RECOMMENDED_POTCAR=TRUE behaviour.
_RECOMMENDED_SYMBOL: dict[str, str] = {
    "Ac": "Ac", "Ag": "Ag", "Al": "Al", "Ar": "Ar", "As": "As",
    "Au": "Au", "B": "B", "Ba": "Ba_sv", "Be": "Be_sv", "Bi": "Bi",
    "Br": "Br", "C": "C", "Ca": "Ca_sv", "Cd": "Cd", "Ce": "Ce",
    "Cl": "Cl", "Co": "Co", "Cr": "Cr_pv", "Cs": "Cs_sv", "Cu": "Cu_pv",
    "Dy": "Dy_3", "Er": "Er_3", "Eu": "Eu", "F": "F", "Fe": "Fe_pv",
    "Ga": "Ga_d", "Gd": "Gd", "Ge": "Ge_d", "H": "H", "He": "He",
    "Hf": "Hf_pv", "Hg": "Hg", "Ho": "Ho_3", "I": "I", "In": "In_d",
    "Ir": "Ir", "K": "K_sv", "Kr": "Kr", "La": "La", "Li": "Li_sv",
    "Lu": "Lu_3", "Mg": "Mg_pv", "Mn": "Mn_pv", "Mo": "Mo_pv",
    "N": "N", "Na": "Na_pv", "Nb": "Nb_pv", "Nd": "Nd_3", "Ne": "Ne",
    "Ni": "Ni_pv", "Np": "Np", "O": "O", "Os": "Os_pv", "P": "P",
    "Pa": "Pa", "Pb": "Pb_d", "Pd": "Pd", "Pm": "Pm_3", "Pr": "Pr_3",
    "Pt": "Pt", "Pu": "Pu", "Rb": "Rb_sv", "Re": "Re_pv", "Rh": "Rh_pv",
    "Ru": "Ru_pv", "S": "S", "Sb": "Sb", "Sc": "Sc_sv", "Se": "Se",
    "Si": "Si", "Sm": "Sm_3", "Sn": "Sn_d", "Sr": "Sr_sv", "Ta": "Ta_pv",
    "Tb": "Tb_3", "Tc": "Tc_pv", "Te": "Te", "Th": "Th", "Ti": "Ti_pv",
    "Tl": "Tl_d", "Tm": "Tm_3", "U": "U", "V": "V_pv", "W": "W_pv",
    "Xe": "Xe", "Y": "Y_sv", "Yb": "Yb_2", "Zn": "Zn", "Zr": "Zr_sv",
}


def get_potcar_lib_dir(server: str | None = None) -> str:
    """Get POTCAR library path from generic environment variables.

    `server` is accepted for compatibility but ignored by the public repo.
    """
    value = (
        os.environ.get("VIBEDFT_POTCAR_DIR")
        or os.environ.get("PMG_VASP_PSP_DIR")
        or os.environ.get("VASP_POTCAR_DIR")
    )
    if not value:
        raise ValueError(
            "No POTCAR library configured. Export VIBEDFT_POTCAR_DIR or PMG_VASP_PSP_DIR."
        )
    return value


def _read_elements(poscar_path: Path) -> list[str]:
    """Extract element symbols from a POSCAR file (header only, no pymatgen).

    The element-symbols line sits at index 5 (after comment, scale, 3 lattice
    vectors). Coordinate-type markers ("Direct"/"Cartesian") come *after* the
    atom-counts line, so they never collide with index 5. We must NOT skip a
    line merely because its first character is 'C'/'S'/'D' -- that wrongly
    discards element lines starting with Cs, Cl, Ca, Co, Cr, Cu, Sb, Sc, ...
    Instead, return the first non-numeric, non-coordinate-marker line at/after
    index 5.
    """
    coord_markers = {"direct", "cartesian", "d", "c", "selectivedynamics"}
    with open(poscar_path) as f:
        lines = f.readlines()
    for i, line in enumerate(lines):
        if i < 5:
            continue
        toks = line.split()
        if not toks:
            continue
        if line.strip().lower() in coord_markers:
            continue
        if all(t.replace(".", "", 1).lstrip("-").isdigit() for t in toks):
            continue  # atom-counts line (numeric)
        return toks
    raise ValueError(f"Cannot parse element line from {poscar_path}")


def generate_potcar(out_dir: str | Path, server: str | None = None) -> Path:
    """Generate POTCAR.

    Priority:
      1. pymatgen ``Potcar``.
      2. Pure-Python concatenation.
    """
    od = Path(out_dir)
    poscar_path = od / "POSCAR"

    if not poscar_path.exists():
        raise FileNotFoundError(f"POTCAR generation requires POSCAR: {poscar_path}")

    elements = _read_elements(poscar_path)

    # 1. pymatgen path
    try:
        from pymatgen.io.vasp import Potcar
        from pymatgen.io.vasp.sets import MPRelaxSet  # noqa: F401

        potcar_map = MPRelaxSet.CONFIG["POTCAR"]
        symbols = [potcar_map[el] for el in elements]
        Potcar(symbols=symbols, functional="PBE").write_file(str(od / "POTCAR"))
        return od / "POTCAR"
    except (ImportError, Exception):
        pass

    # 2. Pure-Python concatenation
    lib_dir = Path(get_potcar_lib_dir(server))
    symbols = [_RECOMMENDED_SYMBOL.get(el, el) for el in elements]

    with open(od / "POTCAR", "w") as out:
        for sym in symbols:
            src = lib_dir / sym / "POTCAR"
            if not src.exists():
                raise FileNotFoundError(
                    f"POTCAR not found: {src}. "
                    "Set VIBEDFT_POTCAR_DIR or check POTCAR library."
                )
            out.write(src.read_text())
    return od / "POTCAR"


def calc_encut(potcar_path: str | Path) -> int:
    """Calculate ENCUT from POTCAR ENMAX values.

    Rule: max_enmax * 1.5, round up to nearest tier in [400, 520, 680].
    Capped at 680 even if raw exceeds it.
    """
    max_enmax = 0.0
    with open(potcar_path) as f:
        for line in f:
            m = re.search(r"ENMAX\s*=\s*([\d.]+)", line)
            if m:
                max_enmax = max(max_enmax, float(m.group(1)))
    raw = max_enmax * 1.5
    for tier in (400, 520, 680):
        if raw <= tier:
            return tier
    return 680
