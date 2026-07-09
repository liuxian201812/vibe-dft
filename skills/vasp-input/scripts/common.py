"""Common utilities for vasp-input-suite.

All shared INCAR template logic, POTCAR generation, ENCUT calculation, and
submit.sh generation are centralized here.
Each phase module imports from this module.
"""
from __future__ import annotations

import re
import sys
import json
import shutil
import subprocess
from pathlib import Path

from pymatgen.core import Structure
from pymatgen.io.vasp import Kpoints, Poscar

# ── sys.path setup for shared_utils ──────────────────────────────
_SKILLS_ROOT = Path(__file__).resolve().parents[3] / "skills"  # repo root → skills/
if str(_SKILLS_ROOT) not in sys.path:
    sys.path.insert(0, str(_SKILLS_ROOT))
from shared_utils import (
    calc_encut, generate_potcar, get_vasp_cmd,
    restart_loop_bash, auto_lda_u, merge_u_config,
)

# ── Templates ────────────────────────────────────────────────────
_SCRIPT_DIR = Path(__file__).resolve().parent
TEMPLATES_DIR = _SCRIPT_DIR / "templates"


def read_template(key: str, templates_dir: Path = TEMPLATES_DIR) -> str:
    """Read a template file (base or modifier) by its key (e.g. 'SR' → INCAR_SR)."""
    fp = templates_dir / f"INCAR_{key}"
    if not fp.exists():
        raise FileNotFoundError(f"Template not found: {fp}")
    return fp.read_text()


def replace_tag(text: str, tag: str, value: str) -> str:
    """Replace parameter value in INCAR text, preserving inline comment."""
    return re.sub(
        rf'^({tag}\s*=\s*)(\S+)',
        rf'\g<1>{value}',
        text, flags=re.MULTILINE,
    )


def fill_placeholders(text: str, params: dict) -> str:
    """Replace {key} placeholders with values from params dict."""
    for k, v in params.items():
        text = text.replace(f"{{{k}}}", str(v))
    unfilled = re.findall(r"\{(\w+)\}", text)
    if unfilled:
        missing = ", ".join(sorted(set(unfilled)))
        raise ValueError(f"Missing required parameter(s): {missing}.")
    return text


def deduplicate_params(incar: str) -> str:
    """Remove duplicate parameters, keeping only the LAST occurrence."""
    lines = incar.split('\n')
    keep = [True] * len(lines)
    seen = {}
    for i, line in enumerate(lines):
        s = line.strip()
        if '=' in s and not s.startswith('#') and not s.startswith('!') and not s.startswith('*'):
            tag = s.split('=')[0].strip().split()[0]
            if tag:
                if tag in seen:
                    keep[seen[tag]] = False
                seen[tag] = i
    return '\n'.join(line for i, line in enumerate(lines) if keep[i])


def apply_functional(incar: str, functional: str) -> str:
    """Apply or switch XC functional in INCAR text.

    Supported: PBE, R2SCAN, HSE06.
    """
    functional = functional.upper()
    lines = incar.splitlines()
    filtered = []
    for line in lines:
        s = line.strip()
        if s.startswith(("GGA", "METAGGA", "LHFCALC", "AEXX", "HFSCREEN")):
            continue
        filtered.append(line)
    incar = "\n".join(filtered)

    if functional == "PBE":
        incar += "\nGGA = PE\nGGA_COMPAT = .FALSE.\n"
    elif functional == "SCAN":
        incar += "\nMETAGGA = SCAN\nLASPH = .TRUE.\nGGA_COMPAT = .FALSE.\nLMIXTAU = .TRUE.\n"
    elif functional == "R2SCAN":
        incar += "\nMETAGGA = R2SCAN\nLASPH = .TRUE.\nGGA_COMPAT = .FALSE.\nLMIXTAU = .TRUE.\n"
    elif functional in ("HSE06", "PBE0", "HYBRID"):
        incar += (
            "\nLHFCALC = .TRUE.\n"
            "AEXX = 0.25\n"
            "HFSCREEN = 0.2\n"
            "ALGO = Normal\n"
        )
    else:
        raise ValueError(f"Unsupported functional: {functional}")
    return incar


def make_pu_modifier(u_config: dict, poscar_path: str = "") -> str | None:
    """Build DFT+U modifier fragment from u_config dict.

    u_config may contain a special key "_LMAXMIX" (int: 4 or 6).
    """
    if not u_config:
        return None
    lmaxmix = u_config.pop("_LMAXMIX", 4)
    elements = list(u_config.keys())
    if poscar_path and Path(poscar_path).exists():
        poscar = Poscar.from_file(poscar_path)
        elements = list(poscar.site_symbols)
    ldau_l = [str(int(u_config.get(e, {}).get("L", -1))) for e in elements]
    ldau_u = [str(u_config.get(e, {}).get("U", 0.0)) for e in elements]
    ldau_j = [str(u_config.get(e, {}).get("J", 0.0)) for e in elements]
    pu = read_template("PU")
    pu = fill_placeholders(pu, {
        "ldaul": " ".join(ldau_l),
        "ldauu": " ".join(ldau_u),
        "ldauj": " ".join(ldau_j),
        "lmaxmix": str(lmaxmix),
    })
    return pu


def validate_combination(template_key: str, modifiers: list[str]):
    """Raise ValueError if template + modifier combination is physically meaningless."""
    allowed_templates = {"SR", "LR", "ST"}
    allowed_modifiers = {"H6", "PU"}
    if template_key not in allowed_templates:
        raise ValueError(f"Unsupported template: {template_key}. Use SR, LR, or ST.")
    unknown = set(modifiers) - allowed_modifiers
    if unknown:
        raise ValueError(f"Unsupported modifier(s): {sorted(unknown)}. Use H6 and/or PU.")


def make_incar(
    template_key: str,
    modifiers: list[str] | None = None,
    encut: int = 400,
    nupdown: int | None = None,
    ncore: int = 12,
    params: dict | None = None,
    isif: int | None = None,
    functional: str = "PBE",
    nelect: int | None = None,
    overrides: dict | None = None,
    ferwe: str | None = None,
    ferdo: str | None = None,
) -> str:
    """Generate INCAR from base template + modifiers with overrides.

    Unified entry point for all phases.
    """
    modifiers = modifiers or []
    validate_combination(template_key, modifiers)

    # 1. Read base template
    incar = read_template(template_key)

    # 2. Append modifiers
    for mod in modifiers:
        incar += "\n" + read_template(mod)

    # 3. Apply functional
    incar = apply_functional(incar, functional)

    # 4. Deduplicate
    incar = deduplicate_params(incar)

    # 5. Fill placeholders
    fill = {"encut": str(encut), "ncore": str(ncore)}
    if params:
        fill.update(params)
    incar = fill_placeholders(incar, fill)

    # 5b. Remove commented ENCUT residue (templates have # ENCUT = {encut})
    incar = re.sub(r'^\s*#\s+ENCUT\s*=\s*\d+.*\n?', '', incar, flags=re.MULTILINE)

    # 6. Override runtime parameters
    incar = replace_tag(incar, "EDIFF", "1E-05")

    if nupdown is not None:
        if re.search(r'^NUPDOWN\s*=', incar, flags=re.MULTILINE):
            incar = replace_tag(incar, "NUPDOWN", str(nupdown))
        else:
            incar += f"\nNUPDOWN = {nupdown}"

    if nelect is not None:
        # Template may have "# NELECT = ..." (commented out); uncomment + set value
        import re as _re
        incar = _re.sub(
            r'^#\s*NELECT\s*=\s*\S*[^\n]*',
            f'NELECT = {nelect}',
            incar, flags=_re.MULTILINE,
        )
        if f"NELECT = {nelect}" not in incar:
            incar += f"\nNELECT = {nelect}"

    # Constrained occupation (FERWE/FERDO) — for ΔSCF singlet/triplet excitations
    if ferwe is not None:
        incar += f"\nFERWE = {ferwe}"
    if ferdo is not None:
        incar += f"\nFERDO = {ferdo}"

    if not re.search(r'^ENCUT\s*=', incar, flags=re.MULTILINE):
        incar += f"\nENCUT = {encut}"

    if "NBLOCK_FOCK" not in incar and any(m in incar for m in ["LHFCALC", "AEXX"]):
        incar += "\nNBLOCK_FOCK = 1"

    if isif is not None:
        if re.search(r'^ISIF\s*=', incar, flags=re.MULTILINE):
            incar = replace_tag(incar, "ISIF", str(isif))
        else:
            incar += f"\nISIF = {isif}"

    if overrides:
        for tag, val in overrides.items():
            if re.search(rf'^{tag}\s*=', incar, flags=re.MULTILINE):
                incar = replace_tag(incar, tag, str(val))
            else:
                incar += f"\n{tag} = {val}"

    return incar.strip() + "\n"


# ── POTCAR / ENCUT ───────────────────────────────────────────────

def write_potcar(workdir: Path, elements: list[str] | None = None, server: str | None = None):
    """Generate POTCAR (pymatgen first, pure-Python fallback)."""
    generate_potcar(str(workdir), server=server)


def get_encut(potcar_path: str | Path) -> int:
    """Calculate ENCUT via three-tier rule (400/520/680)."""
    return calc_encut(str(potcar_path))


def get_auto_u_config(poscar_path: str | Path, manual_u: dict | None = None) -> dict | None:
    """Auto-detect DFT+U params + manual override."""
    struct = Structure.from_file(str(poscar_path)) if Path(poscar_path).exists() else None
    auto_u = auto_lda_u(struct) if struct else None
    return merge_u_config(auto_u, manual_u) if (auto_u or manual_u) else None


def slurm_alloc_lines(ncpus: int, server: str = "") -> list[str]:
    """Return SBATCH allocation lines.

    Request only the total MPI task count. Do not constrain node count or
    tasks-per-node; SLURM should decide how many nodes are needed for ``-n``.
    """
    return [f"#SBATCH -n {ncpus}"]


def write_submit_sh(
    workdir: Path,
    task_dir: str,
    vasp_bin: str = "vasp_std",
    ncores: int = 1,
    partition: str = "",
    node: int = 1,
    qos: str = "",
    vasp_cmd: str = "",
    parent_dirs: list[str] | None = None,
    restart_loop: bool = False,
    max_loops: int = 5,
):
    """Write submit.sh for a single VASP task.

    If restart_loop=True, wraps VASP run in an auto-restart loop that copies
    CONTCAR→POSCAR until ``reached required accuracy`` appears in the output.
    Use for structure relaxation only (not single-point).
    """
    parent_dirs = parent_dirs or []
    if not partition:
        raise ValueError("partition is required for submit.sh generation")
    if not ncores:
        raise ValueError("ncores must be a positive integer for submit.sh generation")
    task_path = workdir / task_dir
    lines = [
        "#!/bin/bash",
        f"#SBATCH -J {task_dir}",
        f"#SBATCH -p {partition}",
        f"#SBATCH -n {ncores}",
    ]
    lines += [
        "#SBATCH --time=36:00:00",
        "#SBATCH -o %j.log",
    ]
    if qos:
        lines.append(f"#SBATCH --qos={qos}")
    lines.append("ulimit -s unlimited")

    lines.append(f"cd {task_path} || exit 1")

    for pd in parent_dirs:
        src = workdir / pd
        lines.append(f'if [[ -s "{src}/CONTCAR" && ! -s POSCAR ]]; then')
        lines.append(f'  cp "{src}/CONTCAR" POSCAR')
        lines.append(f'  cp "{src}/POTCAR" . 2>/dev/null || true')
        lines.append(f'  cp "{src}/KPOINTS" . 2>/dev/null || true')
        lines.append('fi')

    lines.append('if [[ -s CONTCAR ]] && [[ CONTCAR -nt POSCAR ]]; then')
    lines.append('  mv POSCAR POSCAR"$SLURM_JOB_ID"; cp CONTCAR POSCAR')
    lines.append('fi')

    # Resolve full VASP command
    if vasp_cmd:
        full_cmd = vasp_cmd.replace("{vasp_bin}", vasp_bin)
    else:
        full_cmd = get_vasp_cmd().replace("{vasp_bin}", vasp_bin)

    if restart_loop:
        from shared_utils.restart import restart_loop_bash
        lines.append(restart_loop_bash(full_cmd, max_loops=max_loops))
    else:
        for p in full_cmd.split("&&"):
            lines.append(p.strip())

    content = "\n".join(lines) + "\n"
    submit_path = task_path / "submit.sh"
    submit_path.write_text(content)
    submit_path.chmod(0o755)


def scf_convergence_wrapper(vasp_cmd: str) -> str:
    """Generate bash script wrapping VASP with ALGO fallback for SCF convergence.

    Convergence check (for single-point calculations, NSW=0):
      - Extract EDIFF from INCAR
      - Read last DAV/CGA/SDA line from OSZICAR
      - Check if abs(dE) <= EDIFF
      - If not converged, fall back to ALGO=All (max 1 retry)

    VASP stdout/stderr goes directly to SLURM log (no temp file).
    """
    return (
        '# ALGO=Normal first attempt\n'
        f'sed -i \'/^ALGO/c\\ALGO = Normal\' INCAR\n'
        f'{vasp_cmd}\n'
        '# Check electronic convergence via OSZICAR\n'
        'TARGET_EDIFF=$(grep -iw "EDIFF" INCAR | awk -F"=" \'{print $2}\' | xargs)\n'
        'LAST_DE=$(grep -E "SDA:|CGA:|DAV:" OSZICAR 2>/dev/null | tail -1 | awk \'{print $4}\')\n'
        'IS_CONVERGED=$(awk -v de="$LAST_DE" -v ediff="$TARGET_EDIFF" \'BEGIN { gsub("E", "e", de); if ( (de>0?de+0:-de+0) <= ediff+0 ) print "true"; else print "false" }\')\n'
        'if [ "$IS_CONVERGED" = "true" ]; then\n'
        '  echo "SCF converged with ALGO=Normal (dE=$LAST_DE <= EDIFF=$TARGET_EDIFF)"\n'
        'else\n'
        '  echo "SCF not converged with ALGO=Normal (dE=$LAST_DE > EDIFF=$TARGET_EDIFF), retrying with ALGO=All"\n'
        f'  sed -i \'/^ALGO/c\\ALGO = All\' INCAR\n'
        f'  {vasp_cmd}\n'
        '  LAST_DE=$(grep -E "SDA:|CGA:|DAV:" OSZICAR 2>/dev/null | tail -1 | awk \'{print $4}\')\n'
        '  IS_CONVERGED=$(awk -v de="$LAST_DE" -v ediff="$TARGET_EDIFF" \'BEGIN { gsub("E", "e", de); if ( (de>0?de+0:-de+0) <= ediff+0 ) print "true"; else print "false" }\')\n'
        '  if [ "$IS_CONVERGED" = "true" ]; then\n'
        '    echo "SCF converged with ALGO=All (dE=$LAST_DE <= EDIFF=$TARGET_EDIFF)"\n'
        '  else\n'
        '    echo "ERROR: SCF did not converge with Normal or All (dE=$LAST_DE > EDIFF=$TARGET_EDIFF)"\n'
        '    exit 1\n'
        '  fi\n'
        'fi\n'
    )
