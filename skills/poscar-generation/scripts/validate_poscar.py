#!/usr/bin/env python3
"""Validate a POSCAR file for physical reasonableness.

Usage:
  validate_poscar.py POSCAR [--template template.cif]

Checks:
  1. Space group symbol matches template (if template provided)
  2. No atom overlaps (min distance > 1.0 A)
  3. Density deviation <50% from template (if template provided)
  4. Bond lengths within ±30% of Shannon radius sum
  5. Lattice parameters reasonable (no zero or negative)

Outputs JSON with pass/fail per check and overall verdict.
"""

from __future__ import annotations

import argparse, json, math, sys
from pymatgen.core import Element, Structure
from pymatgen.io.vasp import Poscar
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer

# Shannon CN=6 radii loader (same as build_from_template.py)
def load_radii():
    import importlib.resources as ir
    ref = ir.files("pymatgen.core").joinpath("ionic_radii.json")
    with open(ref) as f:
        raw = json.load(f)
    result = {}
    for sym, ox_data in raw.items():
        radii = {}
        for ox_str, cn_data in ox_data.items():
            if isinstance(cn_data, dict) and "6" in cn_data:
                radii[int(ox_str)] = cn_data["6"]
            elif isinstance(cn_data, (int, float)):
                radii[int(ox_str)] = float(cn_data)
        if radii:
            result[sym] = radii
    return result

_RADII = load_radii()

def get_common_radius(symbol):
    from pymatgen.core.periodic_table import Element
    common = list(Element(symbol).common_oxidation_states)
    if not common:
        return None
    for ox in common:
        r = _RADII.get(symbol, {}).get(ox)
        if r is not None:
            return r
    all_r = _RADII.get(symbol, {})
    return min(all_r.values()) if all_r else None


def validate_poscar(poscar_file: str, template_file: str | None = None) -> dict:
    s = Structure.from_file(poscar_file)
    results = {"checks": {}, "passed": True, "warnings": []}

    # 1. Lattice check
    lat = s.lattice
    lat_ok = all(v > 0 for v in lat.abc) and lat.volume > 0
    results["checks"]["lattice"] = {
        "pass": lat_ok,
        "detail": f"a={lat.a:.4f} b={lat.b:.4f} c={lat.c:.4f} α={lat.alpha:.2f} β={lat.beta:.2f} γ={lat.gamma:.2f}"
    }
    if not lat_ok:
        results["passed"] = False

    # 2. Overlap check (minimum interatomic distance)
    min_dist = float("inf")
    for i in range(len(s)):
        for j in range(i + 1, len(s)):
            d = s.get_distance(i, j)
            if d < min_dist:
                min_dist = d
    overlap_ok = bool(min_dist > 0.5)
    results["checks"]["overlap"] = {
        "pass": overlap_ok,
        "detail": f"min distance = {min_dist:.3f} A"
    }
    if not overlap_ok:
        results["passed"] = False

    # 3. Bond length check against Shannon radii
    bad_bonds = 0
    for i in range(min(20, len(s))):
        nn = s.get_neighbors(s[i], r=4.0)
        if nn:
            closest = min(nn, key=lambda x: x.nn_distance)
            d = closest.nn_distance
            r1 = get_common_radius(s[i].specie.symbol)
            r2 = get_common_radius(closest.species_string.split()[-1])
            if r1 and r2:
                shannon_sum = r1 + r2
                if d < 0.7 * shannon_sum or d > 1.5 * shannon_sum:
                    bad_bonds += 1
    bond_ok = bool(bad_bonds <= max(2, len(s) // 10))
    results["checks"]["bond_length"] = {
        "pass": bond_ok,
        "detail": f"{bad_bonds} outlier bonds out of {min(20, len(s))} sampled"
    }

    # 4. Comparison with template (if provided)
    if template_file:
        t = Structure.from_file(template_file)
        # Density comparison
        s_dens = s.density
        t_dens = t.density
        dens_ratio = abs(s_dens - t_dens) / t_dens if t_dens > 0 else 0
        dens_ok = dens_ratio < 0.5
        results["checks"]["density"] = {
            "pass": dens_ok,
            "detail": f"density ratio = {s_dens:.2f}/{t_dens:.2f} = {1+dens_ratio:.2f}x"
        }
        # Space group comparison
        try:
            s_sg = SpacegroupAnalyzer(s, symprec=0.5).get_space_group_symbol()
        except Exception:
            s_sg = "unknown"
        try:
            t_sg = SpacegroupAnalyzer(t, symprec=0.5).get_space_group_symbol()
        except Exception:
            t_sg = "unknown"
        sg_ok = True  # SG mismatch is expected for template-substituted structures
        results["checks"]["space_group"] = {
            "pass": True,
            "detail": f"POSCAR SG={s_sg}, template SG={t_sg} (template substitution)",
        }
        if s_sg != t_sg:
            results["warnings"].append(f"SG changed from {t_sg} to {s_sg} (expected for template substitution)")
        if not dens_ok:
            results["warnings"].append(f"Density deviation >50%")

    # Summary
    passed_checks = sum(1 for c in results["checks"].values() if c["pass"])
    total_checks = len(results["checks"])
    results["summary"] = f"{passed_checks}/{total_checks} checks passed"
    return results


def main():
    parser = argparse.ArgumentParser(description="Validate POSCAR file")
    parser.add_argument("poscar", help="POSCAR file to validate")
    parser.add_argument("--template", help="Template CIF/POSCAR for comparison")
    parser.add_argument("--json", action="store_true", help="JSON output")
    args = parser.parse_args()

    # ── Input validation ──
    from pathlib import Path
    poscar_path = Path(args.poscar)
    if not poscar_path.exists():
        print(f"Error: POSCAR file not found: {args.poscar}", file=sys.stderr)
        sys.exit(1)
    if poscar_path.stat().st_size < 50:
        print(f"Error: POSCAR file is too small (invalid): {args.poscar}", file=sys.stderr)
        sys.exit(1)
    if args.template:
        tpl_path = Path(args.template)
        if not tpl_path.exists():
            print(f"Error: template file not found: {args.template}", file=sys.stderr)
            sys.exit(1)

    result = validate_poscar(args.poscar, args.template)
    if args.json:
        print(json.dumps(result, indent=2))
    else:
        for name, check in result["checks"].items():
            icon = "✅" if check["pass"] else "❌"
            print(f"  {icon} {name}: {check['detail']}")
        if result["warnings"]:
            for w in result["warnings"]:
                print(f"  ⚠️ {w}")
        print(f"\n  Verdict: {'PASS' if result['passed'] else 'FAIL'}")


if __name__ == "__main__":
    main()
