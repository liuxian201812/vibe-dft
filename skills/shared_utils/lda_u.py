"""DFT+U automatic parameter lookup.

Rules:
  - U values: MP (primary) + MIT (supplement for S environment and extra TM)
  - Anion detection: by electronegativity (highest = anion)
  - Anion lookup priority: O, F, S have own tables; Cl/Br/I/N fallback to O
  - LMAXMIX: 4 for d-electron, 6 for f-electron
  - Auto + manual override: auto-detected values can be overridden by manual config

Usage:
  from shared_utils.lda_u import auto_lda_u, merge_u_config
  u_config = auto_lda_u(structure)  # auto-detect from POSCAR
  u_config = merge_u_config(u_config, manual_config)  # apply manual overrides
"""

from pymatgen.core import Structure, Element
from typing import Optional

# DFT+U table: {anion: {element: {"L": l, "U": u, "J": j}}}
# Sources:
#   MP (MPRelaxSet): Co, Cr, Fe, Mn, Mo, Ni, V, W in O/F environment
#   MIT (MITRelaxSet): Ag, Cu, Nb, Re, Ta + S environment for Fe/Mn
#   Cl/Br/I/N: fallback to O table
_U_TABLE = {
    "O": {
        "Co": {"L": 2, "U": 3.32, "J": 0.0},
        "Cr": {"L": 2, "U": 3.7,  "J": 0.0},
        "Fe": {"L": 2, "U": 5.3,  "J": 0.0},
        "Mn": {"L": 2, "U": 3.9,  "J": 0.0},
        "Mo": {"L": 2, "U": 4.38, "J": 0.0},
        "Ni": {"L": 2, "U": 6.2,  "J": 0.0},
        "V":  {"L": 2, "U": 3.25, "J": 0.0},
        "W":  {"L": 2, "U": 6.2,  "J": 0.0},
        # MIT supplements (not in MP):
        "Ag": {"L": 2, "U": 1.5,  "J": 0.0},
        "Cu": {"L": 2, "U": 4.0,  "J": 0.0},
        "Nb": {"L": 2, "U": 1.5,  "J": 0.0},
        "Re": {"L": 2, "U": 2.0,  "J": 0.0},
        "Ta": {"L": 2, "U": 2.0,  "J": 0.0},
    },
    "F": {
        # MP: F uses same U values as O
        "Co": {"L": 2, "U": 3.32, "J": 0.0},
        "Cr": {"L": 2, "U": 3.7,  "J": 0.0},
        "Fe": {"L": 2, "U": 5.3,  "J": 0.0},
        "Mn": {"L": 2, "U": 3.9,  "J": 0.0},
        "Mo": {"L": 2, "U": 4.38, "J": 0.0},
        "Ni": {"L": 2, "U": 6.2,  "J": 0.0},
        "V":  {"L": 2, "U": 3.25, "J": 0.0},
        "W":  {"L": 2, "U": 6.2,  "J": 0.0},
        "Ag": {"L": 2, "U": 1.5,  "J": 0.0},
        "Cu": {"L": 2, "U": 4.0,  "J": 0.0},
        "Nb": {"L": 2, "U": 1.5,  "J": 0.0},
        "Re": {"L": 2, "U": 2.0,  "J": 0.0},
        "Ta": {"L": 2, "U": 2.0,  "J": 0.0},
    },
    "S": {
        # MIT: S environment, smaller U (more covalent)
        "Fe": {"L": 2, "U": 1.9, "J": 0.0},
        "Mn": {"L": 2, "U": 2.5, "J": 0.0},
    },
}

# Elements with known d/f electrons that may need +U
# 3d: Sc~Zn, 4d: Y~Cd, 5d: La~Hg (excluding main-group like In, Sn, etc.)
# 4f: La~Lu (lanthanides)
_D_ELEMENTS = {
    # 3d
    "Sc", "Ti", "V", "Cr", "Mn", "Fe", "Co", "Ni", "Cu", "Zn",
    # 4d
    "Y", "Zr", "Nb", "Mo", "Tc", "Ru", "Rh", "Pd", "Ag", "Cd",
    # 5d (TM, not main group)
    "Hf", "Ta", "W", "Re", "Os", "Ir", "Pt", "Au", "Hg",
}

_F_ELEMENTS = {
    # 4f lanthanides
    "La", "Ce", "Pr", "Nd", "Pm", "Sm", "Eu", "Gd",
    "Tb", "Dy", "Ho", "Er", "Tm", "Yb", "Lu",
    # 5f actinides
    "Ac", "Th", "Pa", "U", "Np", "Pu", "Am", "Cm",
}

# Anion fallback: Cl/Br/I/N/P/As -> O table (most similar with table)
_FALLBACK_ANION = "O"


def _detect_anion(structure: Structure) -> str:
    """Detect the primary anion by electronegativity.

    Returns the symbol of the most electronegative element.
    """
    max_en = -1.0
    anion = None
    for elem in structure.composition.elements:
        try:
            en = elem.X
        except (AttributeError, TypeError):
            continue
        if en is not None and en > max_en:
            max_en = en
            anion = str(elem)
    return anion or "O"


def _resolve_anion_table(anion: str) -> dict:
    """Get the U-value table for a given anion environment.

    Cl/Br/I/N/P/As fallback to O table (no dedicated fitting available).
    """
    if anion in _U_TABLE:
        return _U_TABLE[anion]
    return _U_TABLE[_FALLBACK_ANION]


def auto_lda_u(structure: Structure) -> Optional[dict]:
    """Auto-detect DFT+U parameters from structure.

    Args:
        structure: pymatgen Structure (from POSCAR)

    Returns:
        u_config dict {element: {"L": l, "U": u, "J": j}} or None if no TM found.
        LMAXMIX is included under key "_LMAXMIX".
    """
    anion = _detect_anion(structure)
    table = _resolve_anion_table(anion)

    # Find TM/f-electron elements in structure that have U values
    elements_in_struct = [str(e) for e in structure.composition.elements]
    u_config = {}
    has_f = False

    for elem in elements_in_struct:
        if elem in _F_ELEMENTS:
            has_f = True
        if elem in table:
            u_config[elem] = dict(table[elem])
        # f-elements not in standard table: use default U=4.0 if no entry
        elif elem in _F_ELEMENTS:
            u_config[elem] = {"L": 3, "U": 4.0, "J": 0.0}

    if not u_config:
        return None

    # LMAXMIX: 6 for f-electrons, 4 for d-electrons only
    u_config["_LMAXMIX"] = 6 if has_f else 4
    return u_config


def merge_u_config(auto_config: Optional[dict], manual_config: Optional[dict]) -> Optional[dict]:
    """Merge auto-detected U config with manual overrides.

    Manual values override auto-detected ones.
    Elements only in manual are added.
    Elements only in auto are kept.

    Args:
        auto_config: from auto_lda_u() or None
        manual_config: user-provided {element: {"L": l, "U": u, "J": j}}

    Returns:
        Merged u_config or None if both are None/empty.
    """
    if not auto_config and not manual_config:
        return None

    merged = {}
    if auto_config:
        # Preserve _LMAXMIX from auto detection
        lmaxmix = auto_config.pop("_LMAXMIX", 4)
        merged.update(auto_config)
        merged["_LMAXMIX"] = lmaxmix
    else:
        # No auto config: start with default LMAXMIX
        merged["_LMAXMIX"] = 4

    if manual_config:
        for elem, vals in manual_config.items():
            if elem in merged:
                # Override individual fields
                merged[elem]["L"] = vals.get("L", merged[elem]["L"])
                merged[elem]["U"] = vals.get("U", merged[elem]["U"])
                merged[elem]["J"] = vals.get("J", merged[elem]["J"])
            else:
                merged[elem] = {
                    "L": vals.get("L", -1),
                    "U": vals.get("U", 0.0),
                    "J": vals.get("J", 0.0),
                }
            # Update LMAXMIX if manual element is f-electron
            if elem in _F_ELEMENTS:
                merged["_LMAXMIX"] = max(merged.get("_LMAXMIX", 4), 6)

    return merged if merged else None
