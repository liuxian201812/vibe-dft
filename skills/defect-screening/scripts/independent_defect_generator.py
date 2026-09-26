#!/usr/bin/env python3
"""Auditable, structure-only point-defect generation without doped at runtime.

Geometric grouping schedules an exploratory first batch; it does not establish
equal relaxed energies. Oxidation states and charge intervals are hypotheses,
not a prediction of thermodynamic stability.
"""

from __future__ import annotations

import argparse
import hashlib
import itertools
import json
import math
from collections import defaultdict
from pathlib import Path

import numpy as np
from pymatgen.analysis.defects.generators import (
    SubstitutionGenerator,
    VacancyGenerator,
    VoronoiInterstitialGenerator,
)
from pymatgen.core import Element, Lattice, PeriodicSite, Structure
from pymatgen.io.vasp import Poscar
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer


def _integer_states(values, label):
    if not isinstance(values, dict):
        raise ValueError(f"{label}: expected a JSON object mapping elements to integers")
    states = {}
    for element, value in values.items():
        symbol = Element(element).symbol
        if isinstance(value, bool) or not isinstance(value, (int, float)) or not math.isfinite(value):
            raise ValueError(f"{label}: {symbol} must have a finite integer oxidation state")
        if int(value) != value:
            raise ValueError(f"{label}: {symbol} requires an integer oxidation state")
        states[symbol] = int(value)
    return states


def resolve_oxidation_states(structure, explicit=None):
    """Use a complete explicit map, site decorations, or one unique neutral guess."""
    elements = set(structure.symbol_set)
    if explicit is not None:
        states = _integer_states(explicit, "explicit host oxidation states")
        if set(states) != elements:
            raise ValueError(f"host oxidation states must specify exactly {sorted(elements)}")
        source = "explicit"
    elif all(hasattr(site.specie, "oxi_state") for site in structure):
        by_element = defaultdict(set)
        for site in structure:
            by_element[site.specie.symbol].add(site.specie.oxi_state)
        if any(len(values) != 1 for values in by_element.values()):
            return None, "REVIEW_REQUIRED: mixed site oxidation states need an explicit per-site policy"
        states = _integer_states({key: next(iter(value)) for key, value in by_element.items()}, "site")
        source = "site_decorations"
    else:
        guesses = structure.composition.oxi_state_guesses()
        try:
            unique = {tuple(sorted(_integer_states(guess, "guess").items())) for guess in guesses}
        except ValueError:
            return None, "REVIEW_REQUIRED: inferred mean oxidation states are non-integer"
        if len(unique) != 1:
            return None, f"REVIEW_REQUIRED: {len(unique)} distinct neutral oxidation-state guesses"
        states = dict(unique.pop())
        source = "unique_neutral_composition_guess"
    if abs(sum(structure.composition.get_atomic_fraction(e) * states[e] for e in elements)) > 1e-6:
        return None, "REVIEW_REQUIRED: host oxidation states do not balance charge"
    return states, source


def _env(structure, frac, remove_center):
    center = structure.lattice.get_cartesian_coords(frac)
    neighbors = sorted(structure.get_sites_in_sphere(center, 10.0), key=lambda n: n.nn_distance)
    if remove_center:
        neighbors = [n for n in neighbors if n.nn_distance > 0.05]
    if not neighbors:
        # Avoid a radius proportional to the lattice length for sparse, large cells.
        distances = [
            (float(structure.lattice.get_distance_and_image(frac, site.frac_coords)[0]), site.specie.symbol)
            for site in structure
        ]
        if remove_center:
            distances = [(dist, el) for dist, el in distances if dist > 0.05]
        if not distances:
            distances = [(min(structure.lattice.get_niggli_reduced_lattice().abc),
                          structure[0].specie.symbol)]
        distance, species = min(distances)
        return {"nearest_angstrom": round(distance, 5),
                "neighbors_angstrom": {species: [round(distance, 5)]}}
    if not neighbors or neighbors[0].nn_distance < 0.05:
        raise ValueError("Cannot establish an unoccupied local environment")
    clearance = float(neighbors[0].nn_distance)
    cutoff = min(9.0, max(5.0, 2.4 * clearance))
    shells = defaultdict(list)
    for neighbor in neighbors:
        if neighbor.nn_distance <= cutoff:
            shells[neighbor.specie.symbol].append(round(float(neighbor.nn_distance), 5))
    return {
        "nearest_angstrom": round(clearance, 5),
        "neighbors_angstrom": {key: sorted(val) for key, val in sorted(shells.items())},
    }


def _env_distance(a, b):
    nearest = min(a["nearest_angstrom"], b["nearest_angstrom"])
    scale = max(0.3, 0.25 * nearest)
    difference = weight_sum = 0.0
    for element in sorted(set(a["neighbors_angstrom"]) | set(b["neighbors_angstrom"])):
        left, right = a["neighbors_angstrom"].get(element, []), b["neighbors_angstrom"].get(element, [])
        for index in range(max(len(left), len(right))):
            x = left[index] if index < len(left) else None
            y = right[index] if index < len(right) else None
            distance = min(v for v in (x, y) if v is not None)
            weight = math.exp(-(distance - nearest) / max(0.3, 0.8 * nearest))
            weight_sum += weight
            difference += weight * (min(1.0, abs(x - y) / scale) if x is not None and y is not None else 1.0)
    return round(difference / weight_sum, 6) if weight_sum else 1.0


def _site_symmetry_count(operations, lattice, frac, tolerance):
    return sum(bool(lattice.get_distance_and_image(frac, op.operate(frac))[0] <= tolerance)
               for op in operations)


def _charge_plan(kind, old_state, new_state, max_abs_charge, alternative_states=()):
    if (kind != "interstitial" and old_state is None) or (kind != "vacancy" and new_state is None):
        return [], {"status": "REVIEW_REQUIRED", "reason": "missing oxidation-state evidence"}
    if kind == "vacancy":
        nominal = -old_state
        endpoints = [0, nominal]
        if abs(old_state) == 1:
            endpoints.append(-nominal)
        rule = "vacancy: neutral through nominal ionization; one opposite-side state only for |host valence|=1"
    else:
        nominal = new_state - (old_state if kind == "substitution" else 0)
        offset = old_state if kind == "substitution" else 0
        endpoints = [0, nominal, *(state - offset for state in alternative_states)]
        if nominal == 0:
            endpoints.extend([-1, 1])
        rule = ("neutral through nominal and tabulated common oxidation-state differences; "
                "+/-1 for nominally isovalent defects; alternatives are not stability predictions")
    low, high = min(endpoints), max(endpoints)
    audit = {"status": "OK", "rule": rule, "nominal_charge": nominal,
             "alternative_oxidation_states": list(alternative_states),
             "bounds": [low, high], "max_abs_charge": max_abs_charge}
    if max(abs(low), abs(high)) > max_abs_charge:
        audit.update(status="REVIEW_REQUIRED", reason="proposed charge exceeds configured safety bound; none discarded")
    return list(range(low, high + 1)), audit


def _supercell(structure, min_atoms, min_image_distance, max_atoms, explicit_matrix=None):
    if len(structure) > max_atoms:
        raise ValueError("input already exceeds max_atoms")
    if explicit_matrix is not None:
        if (
            not isinstance(explicit_matrix, (list, tuple))
            or len(explicit_matrix) != 3
            or any(not isinstance(row, (list, tuple)) or len(row) != 3 for row in explicit_matrix)
            or any(
                isinstance(value, bool) or not isinstance(value, int)
                for row in explicit_matrix for value in row
            )
            or any(
                explicit_matrix[i][j] != (explicit_matrix[i][i] if i == j else 0)
                for i in range(3) for j in range(3)
            )
            or any(explicit_matrix[i][i] < 1 for i in range(3))
        ):
            raise ValueError("supercell_matrix must be a positive diagonal 3x3 integer matrix")
        factors = [explicit_matrix[i][i] for i in range(3)]
        atom_count = len(structure) * math.prod(factors)
        if not min_atoms <= atom_count <= max_atoms:
            raise ValueError("explicit supercell_matrix violates min_atoms/max_atoms")
        trial = structure.lattice.matrix * np.array(factors)[:, None]
        image = min(Lattice(trial).get_niggli_reduced_lattice().abc)
        if image + 1e-7 < min_image_distance:
            raise ValueError("explicit supercell_matrix violates min_image_distance")
        return factors, round(float(image), 6)
    limit = max(1, max_atoms // len(structure))
    side = max(math.ceil(limit ** (1 / 3)) + 2,
               *(math.ceil(min_image_distance / length) + 2 for length in structure.lattice.abc))
    choices = (
        factors for factors in itertools.product(range(1, min(side, limit) + 1), repeat=3)
        if max(min_atoms, 2) <= len(structure) * math.prod(factors) <= max_atoms
    )
    for factors in sorted(choices, key=lambda f: (math.prod(f), max(f) / min(f), f)):
        trial = structure.lattice.matrix * np.array(factors)[:, None]
        image = min(Lattice(trial).get_niggli_reduced_lattice().abc)
        if image + 1e-7 >= min_image_distance:
            return list(factors), round(float(image), 6)
    raise ValueError("no diagonal supercell satisfies min_atoms, min_image_distance and max_atoms")


def _mapped_index(host_supercell, original_site, target_frac=None):
    lattice = host_supercell.lattice
    # Convert every supercell atom back into the source lattice to identify all
    # translated replicas, not just the atom at the original Cartesian corner.
    matches = [
        index for index, site in enumerate(host_supercell)
        if site.specie.symbol == original_site.specie.symbol
        and original_site.lattice.get_distance_and_image(
            original_site.frac_coords,
            original_site.lattice.get_fractional_coords(site.coords),
        )[0] < 0.01
    ]
    if not matches:
        raise ValueError("cannot map defect position onto the pristine supercell")
    center = np.asarray(target_frac if target_frac is not None else [0.5, 0.5, 0.5])
    return min(matches, key=lambda index: (
        np.linalg.norm(lattice.get_cartesian_coords(host_supercell[index].frac_coords - center)), index
    ))


def _build_defect(host_supercell, matrix, item):
    result = host_supercell.copy()
    frac = np.asarray(item["fractional_coordinates"])
    if item["kind"] == "interstitial":
        target = np.mod(frac / matrix, 1)
        # Any supercell image is equivalent; select the one nearest the cell center.
        images = [np.mod(target + np.array(shift) / matrix, 1)
                  for shift in itertools.product(*(range(n) for n in matrix))]
        placed = min(images, key=lambda pos: (
            np.linalg.norm(result.lattice.get_cartesian_coords(pos - 0.5)), tuple(pos)
        ))
        result.append(item["new_element"], placed, validate_proximity=False)
        mapped = len(result) - 1
    else:
        source = PeriodicSite(item["old_element"], frac, item["_source_lattice"])
        mapped = _mapped_index(host_supercell, source)
        placed = np.mod(result[mapped].frac_coords, 1)
        if item["kind"] == "vacancy":
            result.remove_sites([mapped])
        else:
            result.replace(mapped, item["new_element"])
    return result, mapped, [round(float(x), 8) for x in placed]


def _select(sites, threshold, budget):
    by_family = defaultdict(list)
    for item in sites:
        if item["screen"]["status"] == "PASS":
            by_family[tuple(item["family"])].append(item)
    groups = []
    for family, members in sorted(by_family.items()):
        distance = lambda a, b: _env_distance(a["environment"], b["environment"]) if a is not b else 0.0
        clusters = [[m] for m in members]
        while True:
            pairs = [
                (max(distance(a, b) for a in left for b in right), i, j)
                for i, left in enumerate(clusters) for j, right in enumerate(clusters)
                if i < j
            ]
            valid = [pair for pair in pairs if pair[0] <= threshold]
            if not valid:
                break
            _, i, j = min(valid)
            clusters[i].extend(clusters.pop(j))
        family_groups = []
        for cluster in clusters:
            representative = min(cluster, key=lambda x: (
                max(distance(x, other) for other in cluster),
                sum(distance(x, other) for other in cluster), x["id"]
            ))
            family_groups.append({"family": list(family), "representative": representative["id"],
                                  "members": [x["id"] for x in cluster], "_rep": representative})
        if family[0] == "interstitial" and len(family_groups) > budget:
            first = max(family_groups, key=lambda g: (
                g["_rep"]["environment"]["nearest_angstrom"], g["representative"]
            ))
            chosen = [first]
            largest = first["_rep"]["environment"]["nearest_angstrom"]
            window = min(0.4, 0.15 * largest)
            eligible = [g for g in family_groups
                        if g["_rep"]["environment"]["nearest_angstrom"] >= largest - window]
            while len(chosen) < min(budget, len(eligible)):
                remaining = [g for g in eligible if g not in chosen]
                chosen.append(max(remaining, key=lambda g: (
                    min(distance(g["_rep"], selected["_rep"]) for selected in chosen),
                    g["_rep"]["environment"]["nearest_angstrom"], g["representative"]
                )))
        else:
            chosen = family_groups
        for group in family_groups:
            selected = group in chosen
            group["first_batch"] = selected
            for member in members:
                if member["id"] in group["members"]:
                    member["selection"] = (
                        "first_batch" if selected and member["id"] == group["representative"] else
                        "deferred_similar" if selected else "deferred_distinct_unrepresented"
                    )
                    member["represented_by"] = group["representative"] if selected else None
            del group["_rep"]
        groups.extend(family_groups)
    for item in sites:
        if item["screen"]["status"] != "PASS":
            item["selection"] = item["screen"]["status"]
            item["represented_by"] = None
    return groups


def generate(poscar, output, *, extrinsic=(), oxidation_states=None, dopant_oxidation_states=None,
             min_atoms=80, max_atoms=512, min_image_distance=10.0, min_void_dist=0.9,
             min_sym_ops=2, threshold=0.15, interstitial_budget=2, max_abs_charge=8,
             symprec=0.01, coverage="first_batch", supercell_matrix=None):
    """Create an audited manifest and only selected, screened defect POSCARs."""
    source, output = Path(poscar).resolve(), Path(output).absolute()
    if output.is_symlink() or (output.exists() and (not output.is_dir() or any(output.iterdir()))):
        raise ValueError("output directory must be new or empty (no overwrites)")
    if coverage not in {"first_batch", "all_pass"}:
        raise ValueError("coverage must be 'first_batch' or 'all_pass'")
    if not (0 <= threshold <= 1 and interstitial_budget >= 1 and min_sym_ops >= 1
            and min_atoms >= 1 and max_atoms >= min_atoms and min_image_distance >= 0
            and min_void_dist > 0 and max_abs_charge >= 1 and symprec > 0):
        raise ValueError("invalid screening, charge, symmetry, void or supercell settings")
    structure = Structure.from_file(source)
    if not structure.is_ordered or not structure.is_valid(tol=0.1):
        raise ValueError("requires an ordered, non-overlapping 3D crystalline host")
    host_elements = set(structure.symbol_set)
    extrinsic = sorted(set(Element(x).symbol for x in extrinsic) - host_elements)
    oxi, oxi_source = resolve_oxidation_states(structure, oxidation_states)
    dopant = _integer_states(dopant_oxidation_states or {}, "dopant oxidation states")
    if set(dopant) - set(extrinsic):
        raise ValueError("dopant oxidation states must refer to listed extrinsic elements")
    matrix, image_distance = _supercell(
        structure, min_atoms, min_image_distance, max_atoms, supercell_matrix,
    )
    clean = structure.copy()
    clean.remove_oxidation_states()
    pristine_clean = clean.copy()
    operations = SpacegroupAnalyzer(clean, symprec=symprec).get_symmetry_operations()
    vacancies = [("vacancy", defect, None) for defect in VacancyGenerator(symprec=symprec).generate(clean)]
    sub_map = {element: sorted((host_elements | set(extrinsic)) - {element}) for element in host_elements}
    substitutions = [
        ("substitution", defect, None)
        for defect in SubstitutionGenerator(symprec=symprec).generate(clean, sub_map)
    ]
    inserted = sorted(host_elements | set(extrinsic))
    void_generator = VoronoiInterstitialGenerator(min_dist=min_void_dist)
    try:
        interstitials = [("interstitial", defect, None)
                         for defect in void_generator.generate(clean, inserted)]
    except ValueError as exc:
        if "not enough values to unpack" not in str(exc):
            raise
        interstitials = []
    sites = []
    raw = vacancies + substitutions + interstitials
    raw.sort(key=lambda entry: (
        entry[0], entry[1].site.specie.symbol,
        entry[1].structure[entry[1].defect_site_index].specie.symbol
        if entry[0] == "substitution" else "",
        tuple(round(float(x % 1), 7) for x in entry[1].site.frac_coords)
    ))
    for index, (kind, defect, _) in enumerate(raw):
        frac = np.mod(defect.site.frac_coords, 1)
        if kind == "substitution":
            old = defect.structure[defect.defect_site_index].specie.symbol
            new = defect.site.specie.symbol
        elif kind == "vacancy":
            old, new = defect.site.specie.symbol, None
        else:
            old, new = None, defect.site.specie.symbol
        family = [kind, old, new]
        reasons = []
        review = []
        if oxi is None:
            review.append(oxi_source)
        if new in extrinsic and new not in dopant:
            review.append(f"extrinsic {new} needs an explicit oxidation state")
        old_state = oxi.get(old) if oxi and old else None
        new_state = (oxi.get(new) if oxi and new in host_elements else dopant.get(new))
        if kind == "substitution":
            en_old, en_new = Element(old).X, Element(new).X
            if en_old is None or en_new is None or math.isnan(en_old) or math.isnan(en_new):
                review.append("Pauling electronegativity unavailable")
            elif abs(en_old - en_new) > 1.0 + 1e-9:
                reasons.append("electronegativity_difference_gt_1.0")
            if old_state is not None and new_state is not None:
                if old_state * new_state > 0 and abs(old_state - new_state) > 3:
                    reasons.append("same_sign_oxidation_difference_gt_3")
        count = _site_symmetry_count(operations, clean.lattice, frac, symprec)
        if kind == "interstitial" and count < min_sym_ops:
            reasons.append(f"site_symmetry_operations_lt_{min_sym_ops}")
        alternatives = (
            sorted(set(Element(new).common_oxidation_states) - {new_state})
            if kind != "vacancy" and new in host_elements and new_state is not None else ()
        )
        charges, charge_audit = _charge_plan(
            kind, old_state, new_state, max_abs_charge, alternatives
        )
        if charge_audit["status"] == "REVIEW_REQUIRED":
            review.append(charge_audit["reason"])
        screen = {"status": "EXCLUDED" if reasons else "REVIEW_REQUIRED" if review else "PASS",
                  "reasons": reasons, "review_reasons": review,
                  "electronegativity_difference": (
                      round(abs(Element(old).X - Element(new).X), 5)
                      if kind == "substitution" and Element(old).X is not None
                      and Element(new).X is not None else None
                  ), "site_symmetry_operations": count}
        original_index = None
        if old is not None:
            original_index = min(
                (i for i, site in enumerate(clean) if site.specie.symbol == old),
                key=lambda i: clean.lattice.get_distance_and_image(frac, clean[i].frac_coords)[0]
            )
        sites.append({
            "id": f"{kind}_{old or 'void'}_{new or 'removed'}_{index:04d}",
            "kind": kind, "family": family, "old_element": old, "new_element": new,
            "fractional_coordinates": [round(float(x), 8) for x in frac],
            "equivalent_fractional_coordinates": [
                [round(float(x % 1), 8) for x in site.frac_coords]
                for site in (defect.equivalent_sites or [defect.site])
            ],
            "parent_site_index": original_index,
            "multiplicity": len(defect.equivalent_sites or [defect.site]),
            "environment": _env(clean, frac, kind != "interstitial"),
            "oxidation_states": {"removed": old_state, "inserted": new_state},
            "charges": charges, "charge_audit": charge_audit, "screen": screen,
        })
    groups = _select(sites, threshold, interstitial_budget)
    first_batch_charge_entries = sum(
        len(site["charges"]) for site in sites if site["selection"] == "first_batch"
    )
    host_supercell = pristine_clean.copy()
    host_supercell.make_supercell(np.diag(matrix))
    pristine_target = output / "pristine_supercell" / "POSCAR"
    pristine_bytes = Poscar(host_supercell).get_str().encode("utf-8")
    pristine_metadata = {
        "poscar": str(pristine_target.resolve()),
        "sha256": hashlib.sha256(pristine_bytes).hexdigest(),
        "atom_count": len(host_supercell),
        "composition": {
            element: int(amount)
            for element, amount in sorted(host_supercell.composition.get_el_amt_dict().items())
        },
        "lattice_vectors": host_supercell.lattice.matrix.tolist(),
    }
    pass_sites = [site for site in sites if site["screen"]["status"] == "PASS"]
    covered_sites = (
        pass_sites if coverage == "all_pass"
        else [site for site in sites if site["selection"] == "first_batch"]
    )
    covered_site_ids = {site["id"] for site in covered_sites}
    coverage_charge_entries = sum(len(site["charges"]) for site in covered_sites)
    pending = []
    for item in sites:
        item["coverage_poscars"] = []
        is_first_batch = item["selection"] == "first_batch"
        is_covered = item["id"] in covered_site_ids
        if not is_first_batch and not is_covered:
            continue
        item["_source_lattice"] = clean.lattice
        defect_cell, mapped, placed = _build_defect(host_supercell, matrix, item)
        del item["_source_lattice"]
        item["supercell_parent_index"] = mapped if item["kind"] != "interstitial" else None
        item["supercell_defect_index"] = mapped if item["kind"] != "vacancy" else None
        item["supercell_fractional_coordinates"] = placed
        poscar_content = Poscar(defect_cell)
        poscar_sha256 = hashlib.sha256(
            poscar_content.get_str().encode("utf-8")
        ).hexdigest()
        if is_first_batch:
            item["first_batch_poscars"] = []
            item["first_batch_poscar_hashes"] = {}
            for charge in item["charges"]:
                target = output / "first_batch" / f"{item['id']}_q{charge:+d}" / "POSCAR"
                item["first_batch_poscars"].append(str(target))
                item["first_batch_poscar_hashes"][str(charge)] = poscar_sha256
                pending.append((target, poscar_content))
        if is_covered:
            coverage_root = "all_pass" if coverage == "all_pass" else "first_batch"
            item["coverage_poscar_hashes"] = {}
            for charge in item["charges"]:
                target = output / coverage_root / f"{item['id']}_q{charge:+d}" / "POSCAR"
                item["coverage_poscars"].append(str(target))
                item["coverage_poscar_hashes"][str(charge)] = poscar_sha256
                if not is_first_batch or coverage == "all_pass":
                    pending.append((target, poscar_content))
    report = {
        "schema_version": 1, "source_poscar": str(source), "host_formula": clean.composition.reduced_formula,
        "host_oxidation_states": oxi, "oxidation_source": oxi_source,
        "extrinsic_oxidation_states": dopant,
        "pristine_supercell": pristine_metadata,
        "supercell": {"matrix": np.diag(matrix).tolist(), "atom_count": len(clean) * math.prod(matrix),
                      "minimum_image_distance_angstrom": image_distance},
        "coverage": {
            "mode": coverage,
            "full_pass": coverage == "all_pass",
            "scope": "screen.status == PASS" if coverage == "all_pass" else "selection == first_batch",
            "site_count": len(covered_sites),
            "charge_entries": coverage_charge_entries,
        },
        "settings": {"min_atoms": min_atoms, "max_atoms": max_atoms,
                     "min_image_distance": min_image_distance, "min_void_dist": min_void_dist,
                     "min_sym_ops": min_sym_ops, "threshold": threshold,
                     "interstitial_budget": interstitial_budget, "max_abs_charge": max_abs_charge,
                     "symprec": symprec, "supercell_matrix": (
                         np.diag(matrix).tolist() if supercell_matrix is not None else None
                     )},
        "policy": "Geometric similarity is only first-batch scheduling, not an energy bound. "
                  "Deferred distinct voids are not represented or physically excluded. "
                  "Unrelaxed charge proposals require electronic validation. "
                  "Voronoi void enumeration and min_void_dist are geometric search limits. "
                  "all_pass coverage writes every PASS site and its existing charge states "
                  "without changing first-batch selection.",
        "sites": sites, "groups": groups,
        "summary": {"all_sites": len(sites), "first_batch_sites": sum(
            s["selection"] == "first_batch" for s in sites
        ), "deferred_distinct_unrepresented": sum(
            s["selection"] == "deferred_distinct_unrepresented" for s in sites
        ), "review_required": sum(s["screen"]["status"] == "REVIEW_REQUIRED" for s in sites),
                    "excluded": sum(s["screen"]["status"] == "EXCLUDED" for s in sites),
                    "all_charge_entries": sum(len(s["charges"]) for s in sites),
                    "first_batch_charge_entries": first_batch_charge_entries,
                    "coverage_sites": len(covered_sites),
                    "coverage_charge_entries": coverage_charge_entries},
    }
    manifest = json.dumps(report, indent=2, ensure_ascii=False) + "\n"
    pristine_target.parent.mkdir(parents=True, exist_ok=True)
    pristine_target.write_bytes(pristine_bytes)
    if hashlib.sha256(pristine_target.read_bytes()).hexdigest() != pristine_metadata["sha256"]:
        raise OSError("pristine supercell POSCAR bytes changed while writing")
    for target, poscar_content in pending:
        target.parent.mkdir(parents=True, exist_ok=True)
        poscar_content.write_file(str(target))
        if hashlib.sha256(target.read_bytes()).hexdigest() != hashlib.sha256(
            poscar_content.get_str().encode("utf-8")
        ).hexdigest():
            raise OSError("defect POSCAR bytes changed while writing")
    output.mkdir(parents=True, exist_ok=True)
    (output / "site_selection.json").write_text(manifest, encoding="utf-8")
    return report


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("poscar", type=Path)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--extrinsic", nargs="*", default=[])
    parser.add_argument("--oxidation-states", type=json.loads, help='Complete JSON map, e.g. {"Na":1,"Cl":-1}')
    parser.add_argument("--dopant-oxidation-states", type=json.loads, help="JSON map for foreign species")
    parser.add_argument("--min-atoms", type=int, default=80)
    parser.add_argument("--max-atoms", type=int, default=512)
    parser.add_argument("--min-image-distance", type=float, default=10.0)
    parser.add_argument(
        "--supercell-matrix", type=json.loads,
        help="Explicit positive diagonal 3x3 supercell matrix (JSON), checked against size/image constraints",
    )
    parser.add_argument("--min-void-dist", type=float, default=0.9)
    parser.add_argument("--min-sym-ops", type=int, default=2)
    parser.add_argument("--threshold", type=float, default=0.15)
    parser.add_argument("--interstitial-budget", type=int, default=2)
    parser.add_argument("--max-abs-charge", type=int, default=8)
    parser.add_argument("--symprec", type=float, default=0.01)
    parser.add_argument("--coverage", choices=("first_batch", "all_pass"), default="first_batch")
    args = parser.parse_args(argv)
    report = generate(
        args.poscar, args.output, extrinsic=args.extrinsic,
        oxidation_states=args.oxidation_states, dopant_oxidation_states=args.dopant_oxidation_states,
        min_atoms=args.min_atoms, max_atoms=args.max_atoms, min_image_distance=args.min_image_distance,
        min_void_dist=args.min_void_dist, min_sym_ops=args.min_sym_ops, threshold=args.threshold,
        interstitial_budget=args.interstitial_budget, max_abs_charge=args.max_abs_charge,
        symprec=args.symprec, coverage=args.coverage,
        supercell_matrix=args.supercell_matrix,
    )
    print(json.dumps(report["summary"], ensure_ascii=False))


if __name__ == "__main__":
    main()
