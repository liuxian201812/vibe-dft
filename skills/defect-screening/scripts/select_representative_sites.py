#!/usr/bin/env python3
"""Generate doped defect sites and plan a structure-only first calculation batch.

Site similarity is an *initial-geometry scheduling heuristic*, not a statement
of equal relaxed formation or optical energies. No VASP inputs are produced.
"""

from __future__ import annotations

import argparse
import json
import math
from collections import defaultdict
from pathlib import Path

from doped.core import Interstitial, Substitution, Vacancy
from doped.generation import DefectsGenerator
from pymatgen.core import Structure
from pymatgen.io.vasp import Poscar


def _identity(entry):
    defect = entry.defect
    if isinstance(defect, Vacancy):
        kind = "vacancy"
    elif isinstance(defect, Substitution):
        kind = "substitution"
    elif isinstance(defect, Interstitial):
        kind = "interstitial"
    else:
        raise ValueError(f"Unsupported defect type: {type(defect).__name__}")
    # doped's uncharged base name (not the point-group-decorated entry name)
    # identifies the replacement species AND the original sublattice.
    return kind, defect.name


def _environment(entry):
    defect = entry.defect
    structure = defect.structure  # pristine host, not the already-defective cell
    center = defect.site.coords
    # All symmetry-equivalent sites have the same environment. Include images
    # but remove the unmodified central atom for vacancies/substitutions.
    near = sorted(structure.get_sites_in_sphere(center, 10.0), key=lambda n: n.nn_distance)
    if not isinstance(defect, Interstitial):
        near = [n for n in near if n.nn_distance > 0.05]
    if not near:
        raise ValueError(f"No pristine neighbors for {entry.name}")
    r0 = near[0].nn_distance
    if r0 < 0.05:
        raise ValueError(f"Candidate overlaps an atom: {entry.name}")
    cutoff = min(9.0, max(5.0, 2.4 * r0))
    shells = defaultdict(list)
    for neighbor in near:
        if neighbor.nn_distance <= cutoff:
            shells[neighbor.specie.symbol].append(round(float(neighbor.nn_distance), 5))
    return {"nearest_angstrom": round(r0, 5),
            "cutoff_angstrom": round(cutoff, 5),
            "neighbors_angstrom": {k: sorted(v) for k, v in sorted(shells.items())}}


def _distance(a, b):
    """Species-resolved soft radial comparison including unmatched neighbors.

    The farther coordination shells contribute less, but no element is
    ignored (in particular, cation neighbors around anion vacancies).
    """
    r0 = min(a["nearest_angstrom"], b["nearest_angstrom"])
    scale = max(0.3, r0 * 0.25)
    weighted_difference = total_weight = 0.0
    for species in sorted(set(a["neighbors_angstrom"]) | set(b["neighbors_angstrom"])):
        aa = a["neighbors_angstrom"].get(species, [])
        bb = b["neighbors_angstrom"].get(species, [])
        for i in range(max(len(aa), len(bb))):
            da = aa[i] if i < len(aa) else None
            db = bb[i] if i < len(bb) else None
            r = min(x for x in (da, db) if x is not None)
            weight = math.exp(-(r - r0) / max(0.3, 0.8 * r0))
            total_weight += weight
            weighted_difference += weight * (
                min(1.0, abs(da - db) / scale) if da is not None and db is not None else 1.0
            )
    return round(weighted_difference / total_weight, 6) if total_weight else 1.0


def plan(entries, threshold=0.15, interstitial_budget=2):
    """Keep every candidate and charge state; mark medoids for first batch.

    Complete-link grouping avoids transitive chaining of dissimilar sites.
    Results depend on pristine geometry and are not a bound on energy errors.
    """
    if not 0 <= threshold <= 1:
        raise ValueError("threshold must be in [0,1]")
    if interstitial_budget < 1:
        raise ValueError("interstitial_budget must be at least 1")
    sites = {}
    for key, entry in sorted(entries.items()):
        family = _identity(entry)
        name = key.rsplit("_", 1)[0]  # only strip a charge suffix
        if name not in sites:
            sites[name] = {
                "family": list(family),
                "host_fractional_coordinates": [round(float(x % 1), 7) for x in entry.defect.site.frac_coords],
                "wyckoff": entry.wyckoff,
                "environment": _environment(entry),
                "charges": [],
                "entry_keys": [],
            }
        elif sites[name]["family"] != list(family):
            raise ValueError(f"Conflicting defect identities: {name}")
        sites[name]["charges"].append(int(entry.charge_state))
        sites[name]["entry_keys"].append(key)

    by_family = defaultdict(list)
    for name, item in sites.items():
        item["charges"].sort()
        by_family[tuple(item["family"])].append(name)
    groups = []
    for family in sorted(by_family):
        names = sorted(by_family[family])
        distances = {
            (a, b): _distance(sites[a]["environment"], sites[b]["environment"])
            for i, a in enumerate(names) for b in names[i + 1:]
        }
        def distance(a, b):
            return distances[tuple(sorted((a, b)))] if a != b else 0.0

        clusters = [[n] for n in names]
        while True:
            candidates = [
                (max(distance(a, b) for a in c for b in d), tuple(c), tuple(d), i, j)
                for i, c in enumerate(clusters) for j, d in enumerate(clusters[i + 1:], i + 1)
            ]
            valid = [x for x in candidates if x[0] <= threshold]
            if not valid:
                break
            _, _, _, i, j = min(valid)
            clusters[i] = sorted(clusters[i] + clusters[j])
            del clusters[j]
        for members in clusters:
            medoid = min(members, key=lambda n: (
                max(distance(n, other) for other in members),
                sum(distance(n, other) for other in members), n
            ))
            groups.append({
                "family": list(family), "representative": medoid,
                "members": members,
                "deferred": [n for n in members if n != medoid],
                "max_pairwise_distance": max(
                    (distance(a, b) for a in members for b in members), default=0
                ),
            })
    # Void-site generation can yield many genuinely different environments.
    # Prioritize a small *exploratory* batch, not pretend the other sites are
    # represented. Both must have near-maximal host clearance; within that
    # window the second maximizes environment diversity. The window is a
    # geometry heuristic, not an ionic-radius or formation-energy threshold.
    selected = set()
    for family in sorted(by_family):
        family_groups = [g for g in groups if tuple(g["family"]) == family]
        if family[0] != "interstitial" or len(family_groups) <= interstitial_budget:
            chosen = family_groups
        else:
            first = max(family_groups, key=lambda g: (
                sites[g["representative"]]["environment"]["nearest_angstrom"],
                g["representative"]
            ))
            chosen = [first]
            largest_void = sites[first["representative"]]["environment"]["nearest_angstrom"]
            clearance_window = min(0.4, 0.15 * largest_void)
            eligible = [g for g in family_groups
                        if sites[g["representative"]]["environment"]["nearest_angstrom"]
                        >= largest_void - clearance_window]
            while len(chosen) < min(interstitial_budget, len(eligible)):
                remaining = [g for g in eligible if g not in chosen]
                chosen.append(max(remaining, key=lambda g: (
                    min(_distance(sites[g["representative"]]["environment"],
                                  sites[c["representative"]]["environment"]) for c in chosen),
                    sites[g["representative"]]["environment"]["nearest_angstrom"],
                    g["representative"]
                )))
        for group in family_groups:
            group["first_batch"] = group in chosen
            if group["first_batch"]:
                selected.add(group["representative"])
            for member in group["members"]:
                sites[member]["status"] = ("first_batch" if group["first_batch"] and
                                           member == group["representative"] else
                                           "deferred_similar" if group["first_batch"] else
                                           "deferred_distinct_unrepresented")
    return {"schema_version": 1, "metric": "species_resolved_weighted_radial_neighborhood",
            "threshold": threshold, "sites": sites, "groups": groups,
            "summary": {"all_sites": len(sites), "first_batch_sites": len(selected),
                        "deferred_sites": len(sites) - len(selected),
                        "deferred_distinct_unrepresented": sum(
                            item["status"] == "deferred_distinct_unrepresented"
                            for item in sites.values()
                        ),
                        "all_charge_entries": len(entries),
                        "first_batch_charge_entries": sum(
                            len(sites[n]["entry_keys"]) for n in selected
                        )}}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("poscar", type=Path)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--extrinsic", nargs="*", help="Optional dopant elements, e.g. Sb")
    parser.add_argument("--threshold", type=float, default=0.15)
    parser.add_argument("--interstitial-budget", type=int, default=2,
                        help="Exploratory void environments per interstitial species (default: 2)")
    parser.add_argument("--min-dist", type=float, default=None, help="Override doped void exclusion distance (Å)")
    parser.add_argument("--min-atoms", type=int, default=80)
    parser.add_argument("--min-image-distance", type=float, default=12.0)
    args = parser.parse_args()
    if args.output.exists() and any(args.output.iterdir()):
        parser.error("output directory must be new or empty; existing structures will not be overwritten")
    structure = Structure.from_file(args.poscar)
    interstitial_kwargs = {"include_unique_wyckoffs": True}
    if args.min_dist is not None:
        interstitial_kwargs["min_dist"] = args.min_dist
    generator = DefectsGenerator(
        structure, extrinsic=args.extrinsic or None,
        interstitial_gen_kwargs=interstitial_kwargs,
        supercell_gen_kwargs={"min_atoms": args.min_atoms,
                              "min_image_distance": args.min_image_distance},
        processes=1,
    )
    report = plan(generator.defect_entries, args.threshold, args.interstitial_budget)
    report["source_poscar"] = str(args.poscar.resolve())
    report["policy"] = (
        "Pristine geometric similarity schedules only the first calculation batch. "
        "No energy or optical-error guarantee; all charge states of each selected site "
        "are kept. Distinct, unrepresented interstitial environments are explicitly "
        "flagged; a two-site exploratory budget is not coverage. Compare deferred "
        "sites if relaxation, formation energies, electronic localization or "
        "spectroscopy contradict the approximation."
    )
    selected = {g["representative"] for g in report["groups"] if g["first_batch"]}
    args.output.mkdir(parents=True, exist_ok=True)
    for name in selected:
        for key in report["sites"][name]["entry_keys"]:
            entry = generator.defect_entries[key]
            if entry.defect_supercell is None:
                raise ValueError(f"Missing defect supercell for {key}")
            target = args.output / "first_batch" / key
            target.mkdir(parents=True, exist_ok=True)
            Poscar(entry.defect_supercell).write_file(str(target / "POSCAR"))
            report["sites"][name].setdefault("first_batch_poscars", []).append(str(target / "POSCAR"))
    (args.output / "site_selection.json").write_text(
        json.dumps(report, indent=2, ensure_ascii=False) + "\n", encoding="utf-8"
    )
    print(json.dumps(report["summary"], ensure_ascii=False))


if __name__ == "__main__":
    main()
