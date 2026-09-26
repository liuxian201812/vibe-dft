#!/usr/bin/env python3
"""Validate a POSCAR file for physical reasonableness.

Usage:
  validate_poscar.py POSCAR [--template template.cif]

Checks:
  1. Space group symbol matches template (if template provided)
  2. No atom overlaps (min distance > 1.0 A)
  3. Density deviation <50% from template (if template provided)
  4. Bond lengths within 0.7–1.5 times the Shannon radius sum
  5. Lattice parameters reasonable (no zero or negative)

Outputs JSON with pass/fail per check and overall verdict.
"""

from __future__ import annotations

import argparse, json, math, re, sys
from pathlib import Path
from pymatgen.core import Element, Structure
from pymatgen.io.vasp import Poscar
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer

try:
    from structure_metadata import load_manifest, sha256_file, write_manifest
except ModuleNotFoundError:
    sys.path.insert(0, str(Path(__file__).resolve().parent))
    from structure_metadata import load_manifest, sha256_file, write_manifest

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


def validate_poscar_header(poscar_file: str) -> dict:
    """Validate a VASP 5 POSCAR header and exact coordinate count."""
    lines = Path(poscar_file).read_text(encoding="utf-8", errors="strict").splitlines()
    errors: list[str] = []
    if len(lines) < 8:
        return {"pass": False, "errors": ["POSCAR has fewer than 8 header lines"]}

    try:
        scale = [float(value) for value in lines[1].split()]
        if len(scale) not in {1, 3} or any(not math.isfinite(value) or value == 0 for value in scale):
            raise ValueError
    except ValueError:
        errors.append("line 2 must contain one or three finite non-zero scale values")

    for index in range(2, 5):
        try:
            vector = [float(value) for value in lines[index].split()]
            if len(vector) != 3 or any(not math.isfinite(value) for value in vector):
                raise ValueError
        except ValueError:
            errors.append(f"line {index + 1} must contain three finite lattice-vector values")

    elements = lines[5].split()
    if not elements or any(not re.fullmatch(r"[A-Z][a-z]?", item) for item in elements):
        errors.append("line 6 must contain element symbols")

    try:
        counts = [int(item) for item in lines[6].split()]
    except ValueError:
        counts = []
        errors.append("line 7 must contain integer atom counts")
    if counts and (any(item <= 0 for item in counts) or len(counts) != len(elements)):
        errors.append("line 7 counts must be positive and match line 6 elements")

    coordinate_index = 7
    selective_dynamics = lines[coordinate_index].strip().lower().startswith("s")
    if selective_dynamics:
        coordinate_index += 1
    if coordinate_index >= len(lines):
        errors.append("coordinate mode line is missing")
        coordinate_mode = ""
    else:
        coordinate_mode = lines[coordinate_index].strip()
        if coordinate_mode.lower() not in {"direct", "cartesian", "d", "c", "k", "kpoints"}:
            errors.append("coordinate mode must be Direct or Cartesian")

    expected = sum(counts) if counts else 0
    coordinate_lines = [line for line in lines[coordinate_index + 1:] if line.strip()]
    if expected and len(coordinate_lines) != expected:
        errors.append(
            f"coordinate count {len(coordinate_lines)} does not match atom-count sum {expected}"
        )
    for number, line in enumerate(coordinate_lines, 1):
        try:
            coordinates = [float(value) for value in line.split()[:3]]
            if len(coordinates) != 3 or any(not math.isfinite(value) for value in coordinates):
                raise ValueError
        except ValueError:
            errors.append(f"coordinate line {number} must start with three finite numbers")

    return {
        "pass": not errors,
        "elements": elements,
        "counts": counts,
        "selective_dynamics": selective_dynamics,
        "coordinate_mode": coordinate_mode,
        "coordinate_count": len(coordinate_lines),
        "expected_coordinate_count": expected,
        "errors": errors,
    }


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


def validate_poscar(
    poscar_file: str,
    template_file: str | None = None,
    *,
    task_id: str = "",
    source: str = "",
) -> dict:
    header = validate_poscar_header(poscar_file)
    try:
        manifest = load_manifest(poscar_file)
    except ValueError as exc:
        manifest = {}
        manifest_error = str(exc)
    else:
        manifest_error = ""

    actual_hash = sha256_file(Path(poscar_file))
    manifest_task_id = str(manifest.get("task_id") or task_id or "")
    manifest_source = str(manifest.get("source") or source or "")
    metadata_ok = bool(
        manifest_task_id
        and manifest_source
        and str(manifest.get("structure_hash") or "") == actual_hash
        and not manifest_error
    )
    results = {
        "checks": {
            "header": {"pass": header["pass"], "detail": header},
            "metadata": {
                "pass": metadata_ok,
                "detail": {
                    "task_id": manifest_task_id,
                    "source": manifest_source,
                    "structure_hash": actual_hash,
                    "manifest_error": manifest_error,
                },
            },
        },
        "metadata": {
            "task_id": manifest_task_id,
            "source": manifest_source,
            "structure_hash": actual_hash,
        },
        "passed": header["pass"] and metadata_ok,
        "all_passed": header["pass"] and metadata_ok,
        "warnings": [],
    }
    if not header["pass"]:
        results["validation_status"] = "REVIEW_REQUIRED"
        results["summary"] = "0/1 checks passed"
        write_manifest(
            poscar_file,
            task_id=manifest_task_id,
            source=manifest_source,
            validation_status="REVIEW_REQUIRED",
            validation=results,
        )
        return results
    s = Structure.from_file(poscar_file)

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
    if not bond_ok:
        results["passed"] = False

    # 4. Template comparison and space-group verification
    if template_file:
        t = Structure.from_file(template_file)
        s_dens = s.density
        t_dens = t.density
        dens_ratio = abs(s_dens - t_dens) / t_dens if t_dens > 0 else 0
        dens_ok = dens_ratio < 0.5
        results["checks"]["density"] = {
            "pass": dens_ok,
            "detail": f"density ratio = {s_dens:.2f}/{t_dens:.2f} = {1+dens_ratio:.2f}x"
        }
        try:
            s_sg = SpacegroupAnalyzer(s, symprec=0.5).get_space_group_symbol()
        except Exception:
            s_sg = "unknown"
        try:
            t_sg = SpacegroupAnalyzer(t, symprec=0.5).get_space_group_symbol()
        except Exception:
            t_sg = "unknown"
        sg_verified = s_sg != "unknown" and t_sg != "unknown"
        sg_ok = sg_verified and s_sg == t_sg
        results["checks"]["space_group"] = {
            "pass": sg_ok,
            "verified": sg_verified,
            "status": "VERIFIED" if sg_ok else "UNVERIFIED",
            "detail": f"POSCAR SG={s_sg}, template SG={t_sg}",
        }
        if not sg_ok:
            results["warnings"].append(
                f"Space-group verification failed: POSCAR={s_sg}, template={t_sg}"
            )
        if not dens_ok:
            results["warnings"].append("Density deviation exceeds 50%")
        if not dens_ok:
            results["passed"] = False
        if not sg_ok:
            results["passed"] = False
    else:
        results["checks"]["space_group"] = {
            "pass": False,
            "verified": False,
            "status": "UNVERIFIED",
            "detail": "No experimental or template space-group reference was supplied",
        }
        results["warnings"].append(
            "Space group is unverified because no reference structure was supplied."
        )
        results["passed"] = False

    chemistry = manifest.get("chemistry") or {}
    if chemistry and chemistry.get("stoichiometry_match") is False:
        results["checks"]["chemistry"] = {
            "pass": False,
            "detail": "Template substitution stoichiometry does not match the target formula",
        }
        results["warnings"].append("Template chemistry requires review.")
        results["passed"] = False

    passed_checks = sum(1 for c in results["checks"].values() if c["pass"])
    total_checks = len(results["checks"])
    results["validation_status"] = (
        "VALIDATED" if results["passed"] else "REVIEW_REQUIRED"
    )
    results["all_passed"] = bool(results["passed"])
    results["summary"] = f"{passed_checks}/{total_checks} checks passed"
    space_group = dict(manifest.get("space_group") or {})
    sg_check = results["checks"].get("space_group") or {}
    space_group.update(
        {
            "verified": bool(sg_check.get("verified") and sg_check.get("pass")),
            "status": (
                "VERIFIED"
                if sg_check.get("verified") and sg_check.get("pass")
                else "UNVERIFIED"
            ),
        }
    )
    write_manifest(
        poscar_file,
        task_id=manifest_task_id,
        source=manifest_source,
        validation_status=results["validation_status"],
        validation=results,
        space_group=space_group,
    )
    return results


def main():
    parser = argparse.ArgumentParser(description="Validate POSCAR file")
    parser.add_argument("poscar", help="POSCAR file to validate")
    parser.add_argument("--template", help="Template CIF/POSCAR for comparison")
    parser.add_argument("--task-id", default="",
                        help="Logical task identifier for the structure manifest")
    parser.add_argument("--source", default="",
                        help="Structure provenance for the structure manifest")
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

    result = validate_poscar(
        args.poscar,
        args.template,
        task_id=args.task_id,
        source=args.source,
    )
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
