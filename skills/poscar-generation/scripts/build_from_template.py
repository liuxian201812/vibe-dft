#!/usr/bin/env python3
"""Build a POSCAR from a template CIF by Wyckoff-based element substitution.

Usage:
  build_from_template.py <template.cif> <target_formula> [-o POSCAR]

Logic:
  1. Read template CIF → get space group, lattice, atom sites with Wyckoff labels
  2. Parse target formula → determine which elements replace which template elements
  3. Group sites by Wyckoff label → group by distance to nearest anion
  4. Assign target elements: larger ionic radius → longer bonds (sites farther from anion)
  5. Write POSCAR with substituted elements

Requires: pymatgen
"""

from __future__ import annotations

import argparse, json, sys
from collections import defaultdict
from pathlib import Path

from pymatgen.core import Element, Structure, Lattice
from pymatgen.io.cif import CifParser
from pymatgen.io.vasp import Poscar

try:
    from structure_metadata import write_manifest
except ModuleNotFoundError:
    sys.path.insert(0, str(Path(__file__).resolve().parent))
    from structure_metadata import write_manifest


def load_ionic_radii() -> dict[str, dict[int, float]]:
    """Load Shannon CN=6 radii from pymatgen's bundled data."""
    import importlib.resources as ir
    import json
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


_RADII = load_ionic_radii()


def get_preferred_radius(symbol: str) -> float:
    """Get the most physically reasonable Shannon radius for an element.

    For metals, prefers positive oxidation states.
    For halogens/chalcogens, prefers negative states.
    Returns the first matching radius or 0.0 if none found.
    """
    radii = _RADII.get(symbol, {})
    if not radii:
        return 0.0
    if symbol in {"F", "Cl", "Br", "I", "O", "S", "Se", "Te"}:
        neg = [r for o, r in radii.items() if o < 0]
        if neg:
            return neg[0]
    pos = [r for o, r in radii.items() if o > 0]
    if pos:
        return pos[0]
    return list(radii.values())[0]


def parse_formula(formula: str) -> dict[str, int]:
    """Parse a chemical formula string into {element: count} dict.

    Handles formats like: Cs2NaInCl6, Cs₂NaInCl₆, Cs2 Na1 In1 Cl6
    """
    from pymatgen.core import Composition
    c = Composition(formula)
    return {str(el): int(amt) for el, amt in c.element_composition.items()}


def build_poscar(template_cif: str, target_formula: str,
                 output: str = "POSCAR",
                 cation_elements: list[str] | None = None,
                 lattice_params: tuple[float, float, float, float, float, float] | None = None,
                 task_id: str = "",
                 source: str = "") -> str:
    """Build POSCAR from template CIF with element substitution.

    Args:
        template_cif: Path to template CIF file.
        target_formula: Target chemical formula (e.g. A2BC6).
        output: Output POSCAR file path.
        cation_elements: List of cation elements (e.g. ['Cs', 'K', 'In']).
        lattice_params: Optional (a, b, c, alpha, beta, gamma) to override template lattice.
        task_id: Logical task identifier recorded in the structure manifest.
        source: Human-readable provenance for the structure.

    Returns path to generated POSCAR file.
    """
    # Parse template
    parser = CifParser(template_cif)
    template = parser.parse_structures(primitive=False)[0]
    from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
    sga = SpacegroupAnalyzer(template)
    sg = sga.get_space_group_symbol()

    # Get target composition
    target_comp = parse_formula(target_formula)
    target_elements = sorted(target_comp.keys())

    # Get template composition
    template_comp = template.composition.element_composition
    template_elements = sorted(template_comp.keys())

    # Compute site sizes for each Wyckoff group
    from collections import defaultdict
    from pymatgen.symmetry.analyzer import SpacegroupAnalyzer

    sga = SpacegroupAnalyzer(template, symprec=0.3)
    try:
        ds = sga.get_symmetry_dataset()
        wyckoff_labels = list(ds.wyckoffs) if hasattr(ds, 'wyckoffs') else []
        equiv_atoms = list(ds.equivalent_atoms) if hasattr(ds, 'equivalent_atoms') else list(range(len(template)))
        use_wyckoff = len(wyckoff_labels) == len(template) and len(set(wyckoff_labels)) > 1
    except Exception:
        use_wyckoff = False

    wyckoff_sites = defaultdict(list)  # group_key → list of (index, avg_distance_to_anion)

    for i, site in enumerate(template):
        # Compute nearest-neighbor distance (site size proxy)
        min_dist = 5.0
        for j, other in enumerate(template):
            if i == j:
                continue
            d = site.distance(other)
            if d < min_dist and d > 0.01:
                min_dist = d
        if min_dist >= 5.0:
            min_dist = 0.0

        # Group key: Wyckoff letter + equivalent atom index, or element + distance if unavailable
        if use_wyckoff:
            key = f"{wyckoff_labels[i]}_{equiv_atoms[i]}"
        else:
            key = f"{site.species.elements[0].symbol}_{min_dist:.2f}"
        wyckoff_sites[key].append((i, min_dist, site.species.elements[0].symbol))

    # Build Wyckoff groups with average site size
    wyckoff_groups = []
    for key, sites in wyckoff_sites.items():
        avg_dist = sum(s[1] for s in sites) / len(sites)
        sym = sites[0][2]
        wyckoff_groups.append({
            "key": key,
            "symbol": sym,
            "count": len(sites),
            "avg_dist": avg_dist,
            "indices": [s[0] for s in sites],
        })
    wyckoff_groups.sort(key=lambda x: -x["avg_dist"])

    # Scale only by an exact formula-unit ratio. Never pad a mismatch with an
    # arbitrary element because that would create a false stoichiometry claim.
    template_n = len(template)
    target_n = sum(target_comp.values())
    scale = (
        template_n // target_n
        if target_n and template_n % target_n == 0
        else None
    )
    scaled_comp = (
        {el: count * scale for el, count in target_comp.items()}
        if scale is not None
        else dict(target_comp)
    )
    stoichiometry_match = (
        scale is not None and sum(scaled_comp.values()) == template_n
    )

    # Determine which elements are cations vs anions
    if not cation_elements:
        cation_elements = []
    cation_set = set(cation_elements)
    # Everything else in scaled_comp that's NOT a cation is an anion
    anion_elements = [el for el in scaled_comp if el not in cation_set]

    # Preserve shared species before assigning replacements. Target cation
    # names do not describe a different cation that is still in the template.
    preserved = {
        el for el, count in scaled_comp.items()
        if float(template_comp.get(el, 0)) == count
    }
    source_cations = set(cation_set)
    for site in template:
        # Recent pymatgen removed PeriodicSite.specie; use the species element.
        symbol = site.species.elements[0].symbol
        oxidation = getattr(site.species.elements[0], "oxi_state", None)
        common_states = Element(symbol).common_oxidation_states
        if (oxidation is not None and oxidation > 0) or (
            oxidation is None and common_states and all(q > 0 for q in common_states)
        ):
            source_cations.add(symbol)
    all_groups = [g for g in wyckoff_groups if g["symbol"] not in preserved]
    cation_groups = [g for g in all_groups if g["symbol"] in source_cations]
    anion_groups = [g for g in all_groups if g["symbol"] not in source_cations]

    # Sort cations by radius descending, anion groups by avg_dist descending
    cation_rank = sorted(
        [(el, get_preferred_radius(el)) for el in cation_elements],
        key=lambda x: -x[1],
    )
    anion_rank = sorted(
        [(el, get_preferred_radius(el)) for el in anion_elements],
        key=lambda x: -x[1],
    )
    cation_groups.sort(key=lambda x: -x["avg_dist"])
    anion_groups.sort(key=lambda x: -x["avg_dist"])

    # Assign cations → cation groups, anions → anion groups
    new_species = [
        site.species.elements[0].symbol
        if site.species.elements[0].symbol in preserved
        else ""
        for site in template
    ]
    remaining = {
        el: 0 if el in preserved else cnt for el, cnt in scaled_comp.items()
    }

    def assign_elements(ranked, groups):
        for target_el, r_el in ranked:
            need = remaining[target_el]
            while need > 0 and groups:
                g = groups.pop(0)
                take = min(need, g["count"])
                for idx in g["indices"][:take]:
                    new_species[idx] = target_el
                need -= take
                if take < g["count"]:
                    g["count"] -= take
                    g["indices"] = g["indices"][take:]
                    groups.insert(0, g)
                remaining[target_el] = need

    assign_elements(cation_rank, cation_groups)
    assign_elements(anion_rank, anion_groups)

    # Keep the POSCAR syntactically usable for review, but retain an explicit
    # mismatch flag if target counts could not be assigned exactly.
    unassigned = [i for i, s in enumerate(new_species) if not s]
    if unassigned:
        all_els = [el for el in scaled_comp if remaining.get(el, 0) > 0]
        if all_els:
            filler_el = all_els[0]
            for i in unassigned:
                new_species[i] = filler_el
        else:
            for i in unassigned:
                new_species[i] = template[i].species.elements[0].symbol

    # Report changes
    from collections import Counter as _Counter
    old_counts = _Counter(site.species.elements[0].symbol for site in template)
    new_counts = _Counter(new_species)
    stoichiometry_match = stoichiometry_match and dict(new_counts) == scaled_comp
    changes = []
    for el in sorted(set(list(old_counts.keys()) + list(new_counts.keys()))):
        oc = old_counts.get(el, 0)
        nc = new_counts.get(el, 0)
        if oc != nc:
            changes.append(f"{el}: {oc}→{nc}")
    print(f"Element changes: {', '.join(changes) if changes else 'none'}")

    # Build structure and sort by element for clean POSCAR output
    if lattice_params:
        a, b, c, alpha, beta, gamma = lattice_params
        new_lattice = Lattice.from_parameters(a, b, c, alpha, beta, gamma)
    else:
        new_lattice = template.lattice
    new_struct = Structure(
        lattice=new_lattice,
        species=new_species,
        coords=[site.frac_coords for site in template],
        coords_are_cartesian=False,
    )
    # Sort by atomic number to group identical elements together
    new_struct.sort()  # pymatgen default: by element Z, deterministic for any system

    # Write POSCAR
    poscar = Poscar(new_struct)
    poscar.write_file(output)

    try:
        generated_sg = SpacegroupAnalyzer(new_struct, symprec=0.3).get_space_group_symbol()
    except Exception:
        generated_sg = "unknown"
    validation_status = "UNVALIDATED" if stoichiometry_match else "REVIEW_REQUIRED"
    warnings = []
    if not stoichiometry_match:
        warnings.append(
            "Generated element counts do not match the target formula scaled to the template."
        )
    warnings.append("Space-group preservation has not been experimentally validated.")
    write_manifest(
        output,
        task_id=task_id,
        source=source or f"COD template: {template_cif}",
        method="cod_template",
        validation_status=validation_status,
        validation={"all_passed": False, "status": validation_status, "warnings": warnings},
        chemistry={
            "template_atom_count": template_n,
            "target_formula": target_formula,
            "target_scaled_counts": dict(scaled_comp),
            "generated_counts": dict(new_counts),
            "stoichiometry_match": stoichiometry_match,
        },
        space_group={
            "template": sg,
            "generated": generated_sg,
            "verified": False,
            "status": "UNVERIFIED",
        },
        template={"source": "COD", "path": str(Path(template_cif).resolve())},
    )
    if not stoichiometry_match:
        print("WARNING: stoichiometry mismatch; structure requires review.")

    # Print summary
    comp_before = template.composition.reduced_formula
    comp_after = new_struct.composition.reduced_formula
    print(f"Template: {comp_before} → Target: {comp_after}")
    print(f"Atoms: {len(new_struct)}")
    print(f"Lattice: {new_struct.lattice}")

    return output


def main():
    parser = argparse.ArgumentParser(
        description="Build POSCAR from template CIF via element substitution")
    parser.add_argument("template_cif", help="Template CIF file")
    parser.add_argument("target_formula", help="Target formula (e.g. Cs2KInCl6)")
    parser.add_argument("-o", "--output", default="POSCAR",
                        help="Output POSCAR file (default: POSCAR)")
    parser.add_argument("--cations", nargs="+", default=[],
                        help="Cation elements (e.g. Cs K In). Elements not in this list are treated as anions.")
    parser.add_argument("--lattice",
                        help="Override lattice: a,b,c,alpha,beta,gamma (e.g. 16.97,16.97,10.99,90,90,90)")
    parser.add_argument("--task-id", default="",
                        help="Logical task identifier for the structure manifest")
    parser.add_argument("--source", default="",
                        help="Structure provenance recorded in the manifest")
    args = parser.parse_args()

    # ── Input validation ──
    from pathlib import Path
    cif_path = Path(args.template_cif)
    if not cif_path.exists():
        print(f"Error: template CIF not found: {args.template_cif}", file=sys.stderr)
        sys.exit(1)
    if cif_path.stat().st_size < 100:
        print(f"Error: template CIF is too small (possibly empty/download failed): {args.template_cif}", file=sys.stderr)
        sys.exit(1)

    from pymatgen.core import Element
    # Validate target formula
    try:
        comp = parse_formula(args.target_formula)
        if not comp:
            print(f"Error: could not parse target formula: {args.target_formula}", file=sys.stderr)
            sys.exit(1)
    except Exception as e:
        print(f"Error: invalid target formula '{args.target_formula}': {e}", file=sys.stderr)
        sys.exit(1)

    # Validate cation elements
    from pymatgen.core.periodic_table import Element
    for el in args.cations:
        try:
            Element(el)
        except Exception as e:
            print(f"Error: invalid cation element '{el}': {e}", file=sys.stderr)
            sys.exit(1)

    # Validate that all elements in formula exist
    for el in comp:
        try:
            Element(el)
        except Exception as e:
            print(f"Error: unknown element '{el}' in formula: {e}", file=sys.stderr)
            sys.exit(1)

    # Parse optional lattice override
    lattice_params = None
    if args.lattice:
        parts = [float(x) for x in args.lattice.split(",")]
        if len(parts) != 6:
            print(f"Error: --lattice requires 6 comma-separated values (a,b,c,alpha,beta,gamma), got {len(parts)}",
                  file=sys.stderr)
            sys.exit(1)
        lattice_params = tuple(parts)

    build_poscar(
        args.template_cif,
        args.target_formula,
        args.output,
        cation_elements=args.cations,
        lattice_params=lattice_params,
        task_id=args.task_id,
        source=args.source,
    )


if __name__ == "__main__":
    main()
