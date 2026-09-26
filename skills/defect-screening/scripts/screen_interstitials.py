#!/usr/bin/env python3
"""Screen interstitial defects using bulk shortest-bond-length criterion.

Usage:
  python3 screen_interstitials.py defect_POSCAR bulk_POSCAR [--coeff 0.75]

Logic:
  1. Read bulk POSCAR, find the shortest interatomic bond length d_min.
  2. Compare bulk and defect compositions:
     - Same total atom count → substitution/vacancy → PASS
     - More atoms in defect → find extra atoms (interstitials)
  3. For each interstitial atom, find its nearest neighbor distance l.
  4. Check: l > coeff × d_min
  5. If not → report UNREASONABLE

Rationale:
  Shannon effective ionic radii are averages over many structures, not hard
  minima. The actual shortest bond in the bulk (d_min) is a better reference
  for whether an interstitial is placed too close.
"""

import argparse
import json

from pymatgen.core import Structure


def get_composition(structure: Structure) -> dict[str, int]:
    counts = {}
    for site in structure:
        el = site.specie.symbol
        counts[el] = counts.get(el, 0) + 1
    return counts


def find_min_bond_length(bulk: Structure) -> float:
    """Compute the minimum center-to-center distance between ANY pair of atoms."""
    min_d = float("inf")
    for i in range(len(bulk)):
        for j in range(i + 1, len(bulk)):
            d = bulk.get_distance(i, j)
            if d < min_d:
                min_d = d
    return min_d


def find_interstitial_indices(defect: Structure, bulk: Structure) -> list[int]:
    if len(defect) <= len(bulk):
        return []

    defect_comp = get_composition(defect)
    bulk_comp = get_composition(bulk)

    # Find elements that are only in defect (e.g. dopant in doped system)
    new_elements = [el for el in defect_comp if el not in bulk_comp]
    if not new_elements:
        return []

    # The extra atoms are at the end of the structure (doped convention)
    extra_count = sum(defect_comp[el] for el in new_elements)
    total = len(defect)
    return list(range(total - extra_count, total))


def screen_interstitial(defect_path: str, bulk_path: str,
                        coeff: float = 0.75) -> list[dict]:
    defect = Structure.from_file(defect_path)
    bulk = Structure.from_file(bulk_path)

    d_min = find_min_bond_length(bulk)
    threshold = coeff * d_min

    inter_indices = find_interstitial_indices(defect, bulk)
    if not inter_indices:
        return []

    results = []
    for idx in inter_indices:
        site = defect[idx]
        species = site.species.elements[0].symbol

        nn_list = defect.get_neighbors(site, r=5.0)
        if not nn_list:
            results.append({
                "index": idx, "species": species,
                "verdict": "PASS", "reason": "no neighbors found",
            })
            continue

        nn = min(nn_list, key=lambda x: x.nn_distance)
        nn_species = nn.species_string.split()[-1]
        l_dist = nn.nn_distance

        if l_dist > threshold:
            verdict = "PASS"
        else:
            verdict = "UNREASONABLE"

        results.append({
            "index": idx, "species": species, "neighbor": nn_species,
            "distance": round(l_dist, 3), "d_min_bulk": round(d_min, 3),
            "threshold": round(threshold, 3), "coeff": coeff,
            "verdict": verdict,
        })

    return results


def main():
    parser = argparse.ArgumentParser(
        description="Screen interstitial defects using bulk shortest-bond-length criterion"
    )
    parser.add_argument("defect_poscar", help="Defect POSCAR")
    parser.add_argument("bulk_poscar", help="Bulk (pristine) supercell POSCAR")
    parser.add_argument("--coeff", type=float, default=0.75,
                        help="Fraction of d_min_bulk to use as threshold (default: 0.75)")
    parser.add_argument("--json", action="store_true")
    args = parser.parse_args()

    results = screen_interstitial(args.defect_poscar, args.bulk_poscar, args.coeff)

    if args.json:
        print(json.dumps(results, indent=2))
        return

    if not results:
        print("No interstitial atoms detected (substitution/vacancy) → PASS")
        return

    all_pass = all(r["verdict"] == "PASS" for r in results)
    for r in results:
        icon = "✅" if r["verdict"] == "PASS" else "❌"
        print(f"  {icon} {r['species']:4s} → {r['neighbor']:4s}  "
              f"l={r['distance']:.3f}A  d_min={r['d_min_bulk']:.3f}A  "
              f"thresh={r['threshold']:.3f}A  {r['verdict']}")

    if all_pass:
        print(f"\nAll {len(results)} interstitial(s) reasonable.")
    else:
        n_fail = sum(1 for r in results if r["verdict"] != "PASS")
        print(f"\n{n_fail}/{len(results)} interstitial(s) UNREASONABLE.")


if __name__ == "__main__":
    main()
