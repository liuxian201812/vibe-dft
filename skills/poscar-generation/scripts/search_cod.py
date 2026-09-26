#!/usr/bin/env python3
"""Search for crystal structure templates from MP + COD.

Three modes:
  search <formula>         Level 1 — search MP (experimental only) + COD
  template <SG> <a> <b> <c> <natoms>  Level 2 — search by space group + lattice + atom count
  download <cod_id>        Download a CIF from COD by its ID

Output: JSON array of candidates sorted by match score.
  Each candidate:
    {source, id, formula, sg, a, b, c, volume, atoms, score, cif_url}
  source is "mp" (Materials Project) or "cod" (Crystallography Open Database).
"""

from __future__ import annotations

import argparse, json, math, os, sys, urllib.request, urllib.parse

COD_CIF = "https://www.crystallography.net/cod/{}.cif"
COD_RESULT = "https://www.crystallography.net/cod/result.php"
MP_API_KEY = os.environ.get("MP_API_KEY")


def _parse_elements(formula: str) -> list[str]:
    """Extract element symbols from a chemical formula."""
    from pymatgen.core import Composition
    try:
        comp = Composition(formula)
        return [str(el) for el in comp.elements]
    except Exception:
        return []


def _post_cod(params: dict) -> list[dict]:
    """POST to COD result.php and return parsed JSON list."""
    params["format"] = "json"
    params.setdefault("include_theoretical", "1")
    data = urllib.parse.urlencode(params).encode()
    try:
        req = urllib.request.Request(COD_RESULT, data=data)
        with urllib.request.urlopen(req, timeout=30) as resp:
            return json.loads(resp.read().decode())
    except Exception as e:
        return [{"error": f"COD search failed: {e}"}]


def search_cod_by_elements(elements: list[str]) -> list[dict]:
    """Search COD by element list (up to 8 elements)."""
    params = {}
    for i, el in enumerate(elements[:8]):
        params[f"el{i+1}"] = el
    raw = _post_cod(params)
    results = []
    for entry in raw:
        if not isinstance(entry, dict):
            continue
        cod_id = entry.get("file", "")
        if not cod_id:
            continue
        results.append({
            "source": "cod",
            "id": cod_id,
            "formula": entry.get("formula", ""),
            "sg": entry.get("sg", ""),
            "a": float(entry.get("a", 0) or 0),
            "b": float(entry.get("b", 0) or 0),
            "c": float(entry.get("c", 0) or 0),
            "volume": float(entry.get("volume", 0) or 0),
            "z": int(entry.get("Z", 0) or 0),
            "score": 1.0,
            "cif_url": COD_CIF.format(cod_id),
        })
    return results


def search_cod_by_template(sg_name: str, a: float, b: float, c: float,
                            natoms: int, formula: str | None = None) -> list[dict]:
    """Search COD by space group + lattice parameters + atom count.

    natoms = total number of atoms in conventional cell.
    COD's Z = formula units per cell. Z = natoms / (atoms per formula).
    If formula is provided, Z is computed exactly and used for filtering.
    """
    vol_target = a * b * c
    params = {"space_group_number": sg_name}

    # Lattice parameter ranges (±10%)
    tol = 0.10
    params["amin"] = f"{a * (1 - tol):.2f}"
    params["amax"] = f"{a * (1 + tol):.2f}"
    params["bmin"] = f"{b * (1 - tol):.2f}"
    params["bmax"] = f"{b * (1 + tol):.2f}"
    params["cmin"] = f"{c * (1 - tol):.2f}"
    params["cmax"] = f"{c * (1 + tol):.2f}"

    # Compute Z from formula + natoms for precise filtering
    if formula:
        from pymatgen.core import Composition
        try:
            comp = Composition(formula)
            atoms_per_fu = sum(comp.element_composition.values())
            z = natoms // atoms_per_fu
            params["minZ"] = str(z)
            params["maxZ"] = str(z)
        except Exception:
            pass  # fall back to no Z filter

    raw = _post_cod(params)
    candidates = []
    for entry in raw:
        if not isinstance(entry, dict):
            continue
        cod_id = entry.get("file", "")
        if not cod_id:
            continue
        ea = float(entry.get("a", 0) or 0)
        eb = float(entry.get("b", 0) or 0)
        ec = float(entry.get("c", 0) or 0)
        ev = float(entry.get("volume", 0) or 0) or (ea * eb * ec)
        e_sg = entry.get("sg", "")
        if not ea or not eb or not ec:
            continue

        # Score by lattice match only (Z may not correlate with atom count here)
        la = abs(ea - a) / a
        lb = abs(eb - b) / b
        lc = abs(ec - c) / c
        lv = abs(ev - vol_target) / vol_target
        if lv > 0.15 or la > 0.15 or lb > 0.15 or lc > 0.15:
            continue

        lat_score = 1.0 - (la + lb + lc) / 3
        vol_score = 1.0 - lv
        score = lat_score * 0.6 + vol_score * 0.4

        candidates.append({
            "source": "cod",
            "id": cod_id,
            "formula": entry.get("formula", ""),
            "sg": e_sg,
            "a": round(ea, 4), "b": round(eb, 4), "c": round(ec, 4),
            "volume": round(ev, 1),
            "z": int(entry.get("Z", 0) or 0),
            "score": round(score, 4),
            "cif_url": COD_CIF.format(cod_id),
        })
    candidates.sort(key=lambda x: -x["score"])
    return candidates[:15]


def search_mp(formula: str, max_results: int = 10) -> list[dict]:
    """Search Materials Project for experimental-only structures (theoretical=False)."""
    if not MP_API_KEY:
        raise RuntimeError("MP_API_KEY environment variable is required for Materials Project fallback search")
    from pymatgen.ext.matproj import MPRester
    results = []
    try:
        with MPRester(MP_API_KEY) as mpr:
            docs = list(mpr.summary.search(
                formula=formula,
                theoretical=False,
            ))
            for d in docs:
                sg = d.get("symmetry", {}).get("symbol", "")
                s = d.get("structure")
                if not s:
                    continue
                results.append({
                    "source": "mp",
                    "id": d.get("material_id", ""),
                    "formula": d.get("formula_pretty", formula),
                    "sg": sg,
                    "a": round(s.lattice.a, 4),
                    "b": round(s.lattice.b, 4),
                    "c": round(s.lattice.c, 4),
                    "volume": round(s.volume, 1),
                    "atoms": d.get("nsites", 0),
                    "score": 1.0,
                    "cif_url": f"mp:{d.get('material_id', '')}",
                    "theoretical": d.get("theoretical", True),
                })
    except Exception as e:
        raise RuntimeError(f"Materials Project search failed for formula={formula}: {e}") from e
    return results


def search_by_formula(formula: str) -> list[dict]:
    """Level 1: search COD first (faster), then MP (experimental only) as fallback."""
    elements = _parse_elements(formula)
    if elements:
        cod_results = search_cod_by_elements(elements)
        if cod_results:
            return cod_results
    # Fallback to MP (experimental only)
    mp_results = search_mp(formula)
    if mp_results:
        return mp_results
    return []


def search_by_template(sg: str, a: float, b: float, c: float,
                        natoms: int, formula: str | None = None) -> list[dict]:
    """Level 2: search COD by space group + lattice + atoms."""
    return search_cod_by_template(sg, a, b, c, natoms, formula=formula)


def download_cif(cod_id: str, output: str = "") -> str:
    """Download a CIF from COD by its ID."""
    url = COD_CIF.format(cod_id)
    try:
        with urllib.request.urlopen(url, timeout=30) as resp:
            cif_text = resp.read().decode()
    except Exception as e:
        return json.dumps([{"error": f"Download failed: {e}"}])

    if not output:
        output = f"{cod_id}.cif"
    with open(output, "w") as f:
        f.write(cif_text)
    return json.dumps([{"cod_id": cod_id, "cif_file": output}])


def main():
    parser = argparse.ArgumentParser(description="Search crystal structure templates")
    sub = parser.add_subparsers(dest="mode", required=True)

    p1 = sub.add_parser("search", help="Level 1: search by formula (MP + COD)")
    p1.add_argument("formula", help="Chemical formula (e.g. AB2C3)")

    p2 = sub.add_parser("template", help="Level 2: search by SG + lattice + atoms")
    p2.add_argument("sg", help="Space group")
    p2.add_argument("a", type=float, help="a (Angstrom)")
    p2.add_argument("b", type=float, help="b (Angstrom)")
    p2.add_argument("c", type=float, help="c (Angstrom)")
    p2.add_argument("natoms", type=int, help="Atoms in conventional cell")
    p2.add_argument("--formula", help="Target formula (e.g. AB2C3), used to compute Z for precise COD filtering")

    p3 = sub.add_parser("download", help="Download CIF from COD")
    p3.add_argument("cod_id", help="COD ID (e.g. 4003575)")
    p3.add_argument("-o", "--output", help="Output file")

    args = parser.parse_args()

    # ── Input validation ──
    import re
    if args.mode == "search":
        if not re.match(r'^[A-Z][a-z]?\d*([A-Z][a-z]?\d*)*$', args.formula.replace('₀₁₂₃₄₅₆₇₈₉', '')):
            print(json.dumps([{"error": f"Invalid formula format: {args.formula}"}]))
            sys.exit(1)

    elif args.mode == "template":
        if args.a <= 0 or args.b <= 0 or args.c <= 0:
            print(json.dumps([{"error": f"Lattice parameters must be positive: a={args.a} b={args.b} c={args.c}"}]))
            sys.exit(1)
        if args.natoms <= 0:
            print(json.dumps([{"error": f"Atom count must be positive: {args.natoms}"}]))
            sys.exit(1)
        if not args.sg or len(args.sg) < 2:
            print(json.dumps([{"error": f"Invalid space group: {args.sg}"}]))
            sys.exit(1)

    elif args.mode == "download":
        if not re.match(r'^\d+$', args.cod_id):
            print(json.dumps([{"error": f"COD ID must be numeric: {args.cod_id}"}]))
            sys.exit(1)

    if args.mode == "search":
        results = search_by_formula(args.formula)
    elif args.mode == "template":
        results = search_by_template(args.sg, args.a, args.b, args.c, args.natoms, formula=args.formula)
    elif args.mode == "download":
        results = download_cif(args.cod_id, args.output or "")
        print(results)
        return

    print(json.dumps(results, indent=2))


if __name__ == "__main__":
    main()
