#!/usr/bin/env python3
"""Coordinate independent structure producers, then freeze reviewed extrinsic states.

The prepare stage makes structures and a review draft, never an approved inventory.
The freeze stage requires separately recorded, structure-hash-bound scientific reviews.
Neither stage creates a VASP request or controls a job.
"""

from __future__ import annotations

import argparse
from copy import deepcopy
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import sys
from typing import Any

import numpy as np
from pymatgen.io.vasp import Poscar

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "doping-supercell" / "scripts"))
sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "defect-screening" / "scripts"))
sys.path.insert(0, str(Path(__file__).resolve().parent))
from core.dope_supercell import (  # noqa: E402
    build_supercell,
    dope_existing_supercell,
    map_supercell_site_to_parent,
)
from build_defect_complex import build_config  # noqa: E402
from independent_defect_generator import generate as generate_defects  # noqa: E402


class InventoryError(ValueError):
    pass


def _object(value: Any, name: str) -> dict[str, Any]:
    if not isinstance(value, dict):
        raise InventoryError(f"{name} must be a JSON object")
    return value


def _unique_pairs(pairs: list[tuple[str, Any]]) -> dict[str, Any]:
    result = {}
    for key, value in pairs:
        if key in result:
            raise InventoryError(f"duplicate JSON key: {key}")
        result[key] = value
    return result


def _read(path: str | Path) -> dict[str, Any]:
    return _object(
        json.loads(Path(path).read_text(encoding="utf-8"), object_pairs_hook=_unique_pairs),
        str(path),
    )


def _write(path: Path, value: dict[str, Any]) -> None:
    path.write_text(
        json.dumps(value, indent=2, ensure_ascii=False, allow_nan=False) + "\n",
        encoding="utf-8",
    )


def _digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _absolute_file(value: Any, name: str) -> Path:
    if not isinstance(value, str):
        raise InventoryError(f"{name} must be an absolute file path")
    path = Path(value)
    if not path.is_absolute() or ".." in path.parts or not path.is_file():
        raise InventoryError(f"{name} must be an existing absolute file without '..'")
    return path.resolve()


def _isolated_interstitial_config(value: Any, dopant: str) -> dict[str, Any] | None:
    if value is None:
        return None
    config = _object(value, "recipe.isolated_interstitial")
    if set(config) - {"enabled", "host_oxidation_states", "dopant_oxidation_state", "max_candidates"}:
        raise InventoryError("isolated_interstitial contains unsupported fields")
    if type(config.get("enabled")) is not bool:
        raise InventoryError("isolated_interstitial.enabled must be explicitly true or false")
    if not config["enabled"]:
        if set(config) != {"enabled"}:
            raise InventoryError("disabled isolated_interstitial must not specify generation parameters")
        return None
    state = config.get("dopant_oxidation_state")
    if isinstance(state, bool) or not isinstance(state, int):
        raise InventoryError(f"isolated_interstitial requires reviewed oxidation-state hypothesis for {dopant}")
    count = config.get("max_candidates", 2)
    if isinstance(count, bool) or not isinstance(count, int) or not 1 <= count <= 2:
        raise InventoryError("isolated_interstitial.max_candidates must be 1 or 2")
    host_states = config.get("host_oxidation_states")
    if host_states is not None and (
        not isinstance(host_states, dict)
        or any(not isinstance(k, str) or isinstance(v, bool) or not isinstance(v, int)
               for k, v in host_states.items())
    ):
        raise InventoryError("isolated_interstitial.host_oxidation_states must map elements to integer states")
    return {"host_oxidation_states": host_states, "dopant_oxidation_state": state, "max_candidates": count}


def _same_pristine(left_path: Path, right_path: Path) -> bool:
    """Require an identical physical parent, including origin and atom order."""
    left = Poscar.from_file(str(left_path), check_for_potcar=False).structure
    right = Poscar.from_file(str(right_path), check_for_potcar=False).structure
    return (
        len(left) == len(right)
        and np.allclose(left.lattice.matrix, right.lattice.matrix, rtol=0, atol=1e-6)
        and all(
            a.specie == b.specie
            and left.lattice.get_distance_and_image(a.frac_coords, b.frac_coords)[0] < 1e-5
            for a, b in zip(left, right)
        )
    )


def _fractional_coordinates(value: Any, name: str) -> list[float]:
    if (
        not isinstance(value, list)
        or len(value) != 3
        or any(
            isinstance(component, bool)
            or not isinstance(component, (int, float))
            or not np.isfinite(component)
            for component in value
        )
    ):
        raise InventoryError(f"{name} must contain three finite fractional coordinates")
    # Periodic images are equivalent; normalize into [0, 1) as the consumer does.
    return [float(component) % 1.0 for component in value]


def _same_fractional_coordinates(left: Any, right: Any, name: str) -> bool:
    left_values = np.asarray(_fractional_coordinates(left, name), dtype=float)
    right_values = np.asarray(_fractional_coordinates(right, name), dtype=float)
    delta = left_values - right_values
    delta -= np.rint(delta)
    return bool(np.all(np.abs(delta) < 1e-5))


def _prepare_isolated_interstitial(
    host: Path, pristine_path: Path, matrix: list[list[int]], dopant: str,
    config: dict[str, Any], output: Path,
) -> tuple[dict[str, Any], list[dict[str, Any]]]:
    defect_output = output / "isolated_interstitial_candidates"
    report = generate_defects(
        host, defect_output, extrinsic=(dopant,),
        oxidation_states=config["host_oxidation_states"],
        dopant_oxidation_states={dopant: config["dopant_oxidation_state"]},
        supercell_matrix=matrix, min_atoms=1, min_image_distance=0,
        coverage="first_batch",
    )
    report_path = defect_output / "site_selection.json"
    report_pristine = _absolute_file(report["pristine_supercell"]["poscar"], "generator pristine")
    if (
        _read(report_path) != report
        or Path(report["source_poscar"]).resolve() != host
        or report["supercell"]["matrix"] != matrix
        or _digest(report_pristine) != report["pristine_supercell"]["sha256"]
        or not _same_pristine(pristine_path, report_pristine)
    ):
        raise InventoryError("interstitial generator parent or supercell differs from substitution parent")
    # The calculation-side interstitial consumer requires the producer report to
    # name the exact frozen pristine POSCAR used by the substitution.  Rebind
    # that identity only after the physical parent check above, and retain the
    # generator's original identity for audit.
    producer_pristine_poscar = str(report_pristine)
    producer_pristine_sha256 = _digest(report_pristine)
    bound_pristine_sha256 = _digest(pristine_path)
    report["pristine_supercell"]["poscar"] = str(pristine_path)
    report["pristine_supercell"]["sha256"] = bound_pristine_sha256
    report["producer_pristine_supercell"] = {
        "path": producer_pristine_poscar,
        "sha256": producer_pristine_sha256,
    }
    _write(report_path, report)
    rebound_report_sha256 = _digest(report_path)
    selected = [
        site for site in report["sites"]
        if site["kind"] == "interstitial" and site["new_element"] == dopant
        and site["screen"]["status"] == "PASS" and site["selection"] == "first_batch"
    ]
    selected = selected[:config["max_candidates"]]
    if not selected:
        raise InventoryError("no screened isolated dopant interstitial candidate; review generator output")
    records = []
    pristine_structure = Poscar.from_file(str(pristine_path), check_for_potcar=False).structure
    for site in selected:
        paths = site["coverage_poscars"]
        if len(paths) != len(site["charges"]) or not paths:
            raise InventoryError("isolated interstitial source lacks generated POSCARs")
        # One geometry per site: generator charge proposals do not authorize electronic states.
        poscar = _absolute_file(paths[0], "generator interstitial POSCAR")
        sha = site["coverage_poscar_hashes"][str(site["charges"][0])]
        if _digest(poscar) != sha:
            raise InventoryError("isolated interstitial POSCAR differs from generator record")
        candidate = Poscar.from_file(str(poscar), check_for_potcar=False).structure
        dopant_sites = [atom for atom in candidate if atom.specie.symbol == dopant]
        host_sites = [atom for atom in candidate if atom.specie.symbol != dopant]
        if (
            len(candidate) != len(pristine_structure) + 1
            or len(dopant_sites) != 1
            or len(host_sites) != len(pristine_structure)
            or not np.allclose(candidate.lattice.matrix, pristine_structure.lattice.matrix, rtol=0, atol=1e-6)
            or any(
                left.specie != right.specie
                or pristine_structure.lattice.get_distance_and_image(
                    left.frac_coords, right.frac_coords,
                )[0] >= 1e-5
                for left, right in zip(pristine_structure, host_sites)
            )
            or pristine_structure.lattice.get_distance_and_image(
                dopant_sites[0].frac_coords, site["supercell_fractional_coordinates"],
            )[0] >= 1e-5
        ):
            raise InventoryError("isolated interstitial is not one dopant inserted into the same pristine supercell")
        records.append({
            "structure_id": f"isolated_{site['id']}",
            "source_site_id": site["id"],
            "poscar": str(poscar),
            "sha256": sha,
            "fractional_coordinates": site["fractional_coordinates"],
            "supercell_fractional_coordinates": site["supercell_fractional_coordinates"],
            "source_selection": site["selection"],
            "source_screen": site["screen"],
            "composition_delta_from_pristine": {dopant: 1},
            "charge_proposals_unreviewed": site["charges"],
            "status": "PENDING_SCIENTIFIC_REVIEW",
        })
    return {
        "site_selection_path": str(report_path),
        "site_selection_sha256": rebound_report_sha256,
        "pristine_poscar": str(pristine_path),
        "pristine_sha256": bound_pristine_sha256,
        "producer_pristine_poscar": producer_pristine_poscar,
        "producer_pristine_sha256": producer_pristine_sha256,
        "parent_validation": (
            "physically identical lattice, ordered atom species and periodic positions; "
            "report pristine identity rebound to the frozen substitution parent"
        ),
    }, records


def prepare(recipe_path: str | Path, output_dir: str | Path) -> Path:
    """Make a fresh two-stage pristine/doped structure and optional complex candidates."""
    recipe_path = _absolute_file(str(recipe_path), "recipe")
    recipe = _read(recipe_path)
    if recipe.get("schema") != "extrinsic-structure-recipe-1":
        raise InventoryError("recipe.schema must be extrinsic-structure-recipe-1")
    if set(recipe) - {"schema", "host_poscar", "supercell_matrix", "substitution", "complex", "isolated_interstitial"}:
        raise InventoryError("recipe contains unsupported fields")
    host = _absolute_file(recipe.get("host_poscar"), "recipe.host_poscar")
    matrix = recipe.get("supercell_matrix")
    if (
        not isinstance(matrix, list) or len(matrix) != 3
        or any(not isinstance(row, list) or len(row) != 3 for row in matrix)
        or any(
            isinstance(matrix[i][j], bool) or not isinstance(matrix[i][j], int)
            or (i != j and matrix[i][j] != 0)
            or (i == j and matrix[i][j] <= 0)
            for i in range(3) for j in range(3)
        )
    ):
        raise InventoryError("supercell_matrix must be positive diagonal integers")
    substitution = _object(recipe.get("substitution"), "recipe.substitution")
    if set(substitution) - {"from", "to", "site_1based", "symprec"}:
        raise InventoryError("substitution contains unsupported fields")
    from_species, dopant = substitution.get("from"), substitution.get("to")
    parent = Poscar.from_file(str(host), check_for_potcar=False)
    if (
        not isinstance(from_species, str) or not isinstance(dopant, str)
        or from_species not in parent.site_symbols or dopant in parent.site_symbols
    ):
        raise InventoryError("substitution must name a host species and a new dopant species")
    isolated_config = _isolated_interstitial_config(recipe.get("isolated_interstitial"), dopant)
    requested_site = substitution.get("site_1based")
    if (
        isinstance(requested_site, bool) or not isinstance(requested_site, int)
        or requested_site < 1
    ):
        raise InventoryError("substitution.site_1based must explicitly select a pristine-supercell atom")
    symprec = substitution.get("symprec", 0.1)
    if isinstance(symprec, bool) or not isinstance(symprec, (int, float)) or not 0 < symprec < 1:
        raise InventoryError("substitution.symprec must be a positive symmetry tolerance below 1 Angstrom")
    complex_recipe = recipe.get("complex")
    if complex_recipe is not None:
        complex_recipe = _object(complex_recipe, "recipe.complex")
        if set(complex_recipe) - {"operations", "selection", "center_charge_hint"}:
            raise InventoryError("complex contains unsupported fields")
        if not isinstance(complex_recipe.get("operations"), list):
            raise InventoryError("complex.operations must be a list")
    output = Path(output_dir).expanduser().resolve()
    if output.exists():
        raise InventoryError(f"output directory must not exist: {output}")
    output.mkdir(parents=True)
    build_file = output / "build.json"
    manifest_path = output / "pristine-manifest.json"
    pristine_path = output / "POSCAR_pristine"
    _write(build_file, {
        "mode": "explicit", "parent_poscar": str(host),
        "supercell_matrix": matrix,
        "output": str(pristine_path), "manifest": str(manifest_path),
    })
    manifest = build_supercell(str(build_file))
    mapped = map_supercell_site_to_parent(
        parent.structure,
        Poscar.from_file(str(pristine_path), check_for_potcar=False).structure,
        matrix,
        requested_site,
        expected_species=from_species,
        symprec=float(symprec),
    )
    if not mapped["parent_wyckoff"]:
        raise InventoryError("target parent Wyckoff could not be determined; review the host structure")
    doping_substitution = {
        "from": from_species,
        "to": dopant,
        "site_1based": requested_site,
        "parent_site": {
            "index_1based": mapped["mapped_parent_site"],
            "wyckoff": mapped["parent_wyckoff"],
            "symprec": float(symprec),
        },
    }
    counts = {key: int(value) for key, value in manifest["pristine_supercell"]["formula"].items()}
    counts[from_species] -= 1
    counts[dopant] = 1
    order = list(manifest["pristine_supercell"]["species_order"])
    index = order.index(from_species)
    if counts[from_species] == 0:
        order[index] = dopant
        del counts[from_species]
    else:
        order.insert(index, dopant)
    dope_file = output / "dope.json"
    provenance_path = output / "doping-provenance.json"
    _write(dope_file, {
        "supercell_manifest": str(manifest_path),
        "substitution": doping_substitution,
        "checks": {
            "formula": counts,
            "atoms": sum(counts.values()),
            "species_order": order,
        },
        "output": str(output / "POSCAR_doped"),
        "provenance": str(provenance_path),
    })
    doping = dope_existing_supercell(str(dope_file))
    isolated_source = None
    isolated_records: list[dict[str, Any]] = []
    if isolated_config is not None:
        isolated_source, isolated_records = _prepare_isolated_interstitial(
            host, pristine_path, matrix, dopant, isolated_config, output,
        )
    complexes: list[dict[str, Any]] = []
    if complex_recipe is not None:
        complex_file = output / "complex.json"
        complex_path = output / "complex_candidates" / "provenance.json"
        _write(complex_file, {
            "doping_provenance": str(provenance_path),
            **complex_recipe,
            "output_dir": str(output / "complex_candidates"),
            "provenance": str(complex_path),
        })
        built = build_config(complex_file)
        complexes = [
            {
                "structure_id": f"complex_{row['rank']:03d}",
                "poscar": row["path"],
                "sha256": row["sha256"],
                "provenance_path": str(complex_path),
                "provenance_sha256": _digest(complex_path),
                "candidate_ids": row["candidate_ids"],
            }
            for row in built["outputs"]
        ]
    draft = output / "review-draft.json"
    _write(draft, {
        "schema": "extrinsic-inventory-review-draft-1",
        "status": "PENDING_SCIENTIFIC_REVIEW",
        "host_poscar": str(host),
        "host_sha256": _digest(host),
        "pristine_manifest": str(manifest_path),
        "pristine_manifest_sha256": _digest(manifest_path),
        "doping_provenance": str(provenance_path),
        "doping_provenance_sha256": _digest(provenance_path),
        "doped": {
            "structure_id": "dopant_substitution",
            "poscar": doping["doped_supercell"]["path"],
            "sha256": doping["doped_supercell"]["sha256"],
        },
        "complexes": complexes,
        "isolated_interstitial_source": isolated_source,
        "isolated_interstitials": isolated_records,
    })
    return draft


_INTERSTITIAL_SOURCE_FIELDS = {
    "site_selection_path", "site_selection_sha256",
    "pristine_poscar", "pristine_sha256",
    "producer_pristine_poscar", "producer_pristine_sha256",
    "parent_validation",
}
_INTERSTITIAL_RECORD_FIELDS = {
    "structure_id", "source_site_id", "poscar", "sha256",
    "fractional_coordinates", "supercell_fractional_coordinates",
    "source_selection", "source_screen", "composition_delta_from_pristine",
    "charge_proposals_unreviewed", "status",
}
_INTERSTITIAL_SITE_AUDIT_FIELDS = {
    "status", "reviewer", "evidence_ref", "source_report_sha256",
    "candidate_sha256", "site_id", "species", "selection", "screen",
    "composition_delta_from_pristine", "fractional_coordinates",
    "supercell_fractional_coordinates",
}


def _interstitial_enabled(draft: dict[str, Any]) -> bool:
    return bool(draft.get("isolated_interstitials")) or (
        draft.get("isolated_interstitial_source") is not None
    )


def _interstitial_delta(pristine_path: Path, candidate_path: Path) -> dict[str, int]:
    pristine = Poscar.from_file(str(pristine_path), check_for_potcar=False).structure
    candidate = Poscar.from_file(str(candidate_path), check_for_potcar=False).structure
    delta: dict[str, int] = {}
    for element, amount in candidate.composition.get_el_amt_dict().items():
        delta[element] = int(amount) - int(
            pristine.composition.get_el_amt_dict().get(element, 0)
        )
    for element, amount in pristine.composition.get_el_amt_dict().items():
        delta.setdefault(element, -int(amount))
    return {element: value for element, value in sorted(delta.items()) if value}


def _candidate_preserves_pristine(
    pristine_path: Path,
    candidate_path: Path,
    species: str,
    supercell_point: list[float],
) -> bool:
    pristine = Poscar.from_file(str(pristine_path), check_for_potcar=False).structure
    candidate = Poscar.from_file(str(candidate_path), check_for_potcar=False).structure
    if (
        len(candidate) != len(pristine) + 1
        or not np.allclose(
            candidate.lattice.matrix, pristine.lattice.matrix, rtol=0, atol=1e-6
        )
    ):
        return False
    unmatched = list(range(len(candidate)))
    for parent in pristine:
        matches = [
            index for index in unmatched
            if candidate[index].specie.symbol == parent.specie.symbol
            and pristine.lattice.get_distance_and_image(
                parent.frac_coords, candidate[index].frac_coords,
            )[0] < 1e-5
        ]
        if len(matches) != 1:
            return False
        unmatched.remove(matches[0])
    if len(unmatched) != 1:
        return False
    inserted = candidate[unmatched[0]]
    return (
        inserted.specie.symbol == species
        and pristine.lattice.get_distance_and_image(
            inserted.frac_coords, supercell_point,
        )[0] < 1e-5
    )


def _validated_interstitial_candidates(
    draft: dict[str, Any], manifest: dict[str, Any], dopant: str,
) -> list[dict[str, Any]]:
    """Rebind every draft candidate to the frozen producer report and pristine source."""
    source = _object(draft.get("isolated_interstitial_source"), "draft.isolated_interstitial_source")
    if set(source) != _INTERSTITIAL_SOURCE_FIELDS:
        raise InventoryError("isolated interstitial source has unsupported or missing fields")
    report_path = _absolute_file(source["site_selection_path"], "interstitial site selection")
    report = _read(report_path)
    if _digest(report_path) != source["site_selection_sha256"]:
        raise InventoryError("interstitial site selection changed after preparation")
    if report.get("schema_version") != 1:
        raise InventoryError("interstitial site selection schema_version must be 1")
    pristine_record = manifest.get("pristine_supercell")
    if not isinstance(pristine_record, dict):
        raise InventoryError("pristine manifest lacks pristine_supercell")
    pristine_path = _absolute_file(pristine_record.get("path"), "pristine manifest POSCAR")
    pristine_sha = _digest(pristine_path)
    if (
        pristine_sha != pristine_record.get("sha256")
        or Path(source["pristine_poscar"]).resolve() != pristine_path
        or source["pristine_sha256"] != pristine_sha
    ):
        raise InventoryError("isolated interstitial source does not bind the frozen pristine supercell")
    parent = _object(manifest.get("parent"), "pristine manifest.parent")
    report_parent = _absolute_file(report.get("source_poscar"), "interstitial report parent")
    if (
        report_parent != _absolute_file(parent.get("path"), "pristine manifest parent")
        or _digest(report_parent) != parent.get("sha256")
    ):
        raise InventoryError("interstitial report does not name the substitution parent POSCAR")
    report_pristine = _object(
        report.get("pristine_supercell"), "interstitial report pristine_supercell",
    )
    if (
        Path(report_pristine.get("poscar", "")).resolve() != pristine_path
        or report_pristine.get("sha256") != pristine_sha
    ):
        raise InventoryError("interstitial report pristine identity conflicts with the frozen source")
    pristine_structure = Poscar.from_file(
        str(pristine_path), check_for_potcar=False
    ).structure
    pristine_counts = {
        str(element): int(amount)
        for element, amount in pristine_structure.composition.get_el_amt_dict().items()
    }
    if (
        report_pristine.get("atom_count") != len(pristine_structure)
        or report_pristine.get("composition") != pristine_counts
        or not np.allclose(
            np.asarray(report_pristine.get("lattice_vectors"), dtype=float),
            pristine_structure.lattice.matrix,
            rtol=0,
            atol=1e-6,
        )
    ):
        raise InventoryError("interstitial report pristine metadata conflicts with its POSCAR")
    if (
        not isinstance(dopant, str)
        or not dopant
        or dopant in pristine_counts
    ):
        raise InventoryError("interstitial dopant must be an external element")
    extrinsic_states = report.get("extrinsic_oxidation_states")
    if not isinstance(extrinsic_states, dict) or dopant not in extrinsic_states:
        raise InventoryError("interstitial report does not declare the external dopant")
    producer = _object(
        report.get("producer_pristine_supercell"),
        "interstitial report producer_pristine_supercell",
    )
    producer_path = _absolute_file(producer.get("path"), "interstitial producer pristine")
    if producer.get("sha256") != _digest(producer_path) or not _same_pristine(
        producer_path, pristine_path,
    ):
        raise InventoryError("interstitial report generator pristine is not the substitution parent")
    generation = _object(manifest.get("generation"), "pristine manifest.generation")
    matrix = generation.get("supercell_matrix")
    if _object(report.get("supercell"), "interstitial report supercell").get("matrix") != matrix:
        raise InventoryError("interstitial report supercell differs from the substitution matrix")
    coverage = _object(report.get("coverage"), "interstitial report coverage")
    if coverage.get("mode") != "first_batch" or coverage.get("full_pass") is not False:
        raise InventoryError("interstitial report coverage must be first_batch")
    sites = report.get("sites")
    if not isinstance(sites, list):
        raise InventoryError("interstitial report sites must be a list")

    records = draft.get("isolated_interstitials")
    if not isinstance(records, list) or not records:
        raise InventoryError("enabled isolated interstitial draft must contain candidates")
    validated: list[dict[str, Any]] = []
    seen_ids: set[str] = set()
    for index, raw_record in enumerate(records):
        record = _object(raw_record, f"draft.isolated_interstitials[{index}]")
        if set(record) != _INTERSTITIAL_RECORD_FIELDS:
            raise InventoryError("isolated interstitial record has unsupported or missing fields")
        structure_id = record["structure_id"]
        if not isinstance(structure_id, str) or not structure_id or structure_id in seen_ids:
            raise InventoryError("isolated interstitial structure IDs must be unique strings")
        seen_ids.add(structure_id)
        site_id = record["source_site_id"]
        matches = [
            site for site in sites
            if isinstance(site, dict) and site.get("id") == site_id
        ]
        if len(matches) != 1:
            raise InventoryError("isolated interstitial draft site is not in the producer report")
        site = matches[0]
        if (
            site.get("kind") != "interstitial"
            or site.get("old_element") is not None
            or site.get("new_element") != dopant
            or site.get("selection") != record["source_selection"]
            or site.get("screen") != record["source_screen"]
            or not isinstance(site.get("screen"), dict)
            or site["screen"].get("status") != "PASS"
        ):
            raise InventoryError("isolated interstitial draft site no longer matches the producer review")
        fractional = _fractional_coordinates(
            record["fractional_coordinates"], "isolated interstitial fractional coordinates",
        )
        supercell_point = _fractional_coordinates(
            record["supercell_fractional_coordinates"],
            "isolated interstitial supercell fractional coordinates",
        )
        if (
            not _same_fractional_coordinates(
                fractional, site.get("fractional_coordinates"), "producer fractional coordinates",
            )
            or not _same_fractional_coordinates(
                supercell_point, site.get("supercell_fractional_coordinates"),
                "producer supercell fractional coordinates",
            )
        ):
            raise InventoryError("isolated interstitial draft coordinates changed after preparation")
        candidate_path = _absolute_file(record["poscar"], "isolated interstitial POSCAR")
        candidate_sha = _digest(candidate_path)
        coverage_poscars = site.get("coverage_poscars")
        charges = site.get("charges")
        coverage_hashes = site.get("coverage_poscar_hashes")
        if (
            not isinstance(coverage_poscars, list)
            or not isinstance(charges, list)
            or not isinstance(coverage_hashes, dict)
        ):
            raise InventoryError("producer interstitial site lacks coverage records")
        if (
            not charges
            or any(isinstance(charge, bool) or not isinstance(charge, int) for charge in charges)
            or len(set(charges)) != len(charges)
            or record["charge_proposals_unreviewed"] != charges
        ):
            raise InventoryError(
                "draft interstitial charge proposals must preserve the producer list"
            )
        indexes = [
            position for position, value in enumerate(coverage_poscars)
            if isinstance(value, str) and Path(value).resolve() == candidate_path
        ]
        if len(indexes) != 1:
            raise InventoryError("draft interstitial POSCAR is not a producer coverage POSCAR")
        producer_sha = coverage_hashes.get(str(charges[indexes[0]]))
        if producer_sha != candidate_sha or record["sha256"] != candidate_sha:
            raise InventoryError("isolated interstitial POSCAR hash conflicts with the producer")
        delta = _interstitial_delta(pristine_path, candidate_path)
        if (
            delta != {dopant: 1}
            or record["composition_delta_from_pristine"] != delta
            or not _candidate_preserves_pristine(
                pristine_path, candidate_path, dopant, supercell_point,
            )
        ):
            raise InventoryError("isolated interstitial is not one dopant in the frozen pristine supercell")
        validated.append({
            "structure_id": structure_id,
            "sha256": candidate_sha,
            "poscar": str(candidate_path),
            "provenance_path": str(report_path),
            "provenance_sha256": source["site_selection_sha256"],
            "site_id": site_id,
            "species": dopant,
            "selection": record["source_selection"],
            "screen": deepcopy(record["source_screen"]),
            "composition_delta_from_pristine": delta,
            "fractional_coordinates": fractional,
            "supercell_fractional_coordinates": supercell_point,
        })
    return validated


def _validated_interstitial_audits(
    review: dict[str, Any], candidates: list[dict[str, Any]],
) -> tuple[list[dict[str, Any]], set[str]]:
    audits = review.get("interstitial_candidate_audits")
    by_id = {candidate["structure_id"]: candidate for candidate in candidates}
    if not isinstance(audits, dict) or set(audits) != set(by_id):
        raise InventoryError(
            "interstitial_candidate_audits must explicitly review every drafted candidate"
        )
    approved_ids: set[str] = set()
    selected: list[dict[str, Any]] = []
    for structure_id in sorted(by_id):
        candidate = by_id[structure_id]
        audit = _object(audits[structure_id], f"interstitial audit {structure_id}")
        status = audit.get("status")
        if not audit.get("reviewer") or not audit.get("evidence_ref"):
            raise InventoryError("interstitial audit requires reviewer and evidence_ref")
        if audit.get("candidate_sha256") != candidate["sha256"]:
            raise InventoryError("interstitial audit does not bind the candidate POSCAR SHA256")
        if status == "rejected":
            if set(audit) != {
                "status", "reviewer", "evidence_ref", "candidate_sha256", "reason",
            } or not audit.get("reason"):
                raise InventoryError("rejected interstitial audit must record a reason")
            continue
        if status != "approved" or set(audit) != _INTERSTITIAL_SITE_AUDIT_FIELDS:
            raise InventoryError("interstitial audit must be an explicit approved or rejected review")
        if (
            audit["source_report_sha256"] != candidate["provenance_sha256"]
            or audit["site_id"] != candidate["site_id"]
            or audit["species"] != candidate["species"]
            or audit["selection"] != candidate["selection"]
            or audit["screen"] != candidate["screen"]
            or audit["composition_delta_from_pristine"]
            != candidate["composition_delta_from_pristine"]
            or not _same_fractional_coordinates(
                audit["fractional_coordinates"], candidate["fractional_coordinates"],
                "audited fractional coordinates",
            )
            or not _same_fractional_coordinates(
                audit["supercell_fractional_coordinates"],
                candidate["supercell_fractional_coordinates"],
                "audited supercell fractional coordinates",
            )
        ):
            raise InventoryError("interstitial audit does not bind the reviewed site and source")
        approved_ids.add(structure_id)
        selected.append({
            "structure_id": structure_id,
            "provenance_path": candidate["provenance_path"],
            "provenance_sha256": candidate["provenance_sha256"],
            "candidate_path": candidate["poscar"],
            "candidate_sha256": candidate["sha256"],
            "site_id": candidate["site_id"],
            "selection": candidate["selection"],
            "screen": deepcopy(candidate["screen"]),
            "composition_delta_from_pristine": deepcopy(
                candidate["composition_delta_from_pristine"]
            ),
            "site_audit": deepcopy(audit),
        })
    return selected, approved_ids


def freeze(draft_path: str | Path, review_path: str | Path, inventory_path: str | Path) -> Path:
    """Require independently approved state/geometry decisions before VASP inventory."""
    draft = _read(_absolute_file(str(draft_path), "draft"))
    review = _read(_absolute_file(str(review_path), "review"))
    if draft.get("schema") != "extrinsic-inventory-review-draft-1":
        raise InventoryError("unsupported review draft")
    interstitial_enabled = _interstitial_enabled(draft)
    if review.get("schema") != "extrinsic-inventory-review-1":
        raise InventoryError("review.schema must be extrinsic-inventory-review-1")
    review_fields = {"schema", "draft_sha256", "states", "complex_candidate_audits"}
    if interstitial_enabled:
        review_fields.add("interstitial_candidate_audits")
    if set(review) != review_fields:
        raise InventoryError(
            "review fields must match the drafted structure families: "
            + ", ".join(sorted(review_fields))
        )
    if review.get("draft_sha256") != _digest(Path(draft_path)):
        raise InventoryError("review must bind the exact draft SHA256")
    host = _absolute_file(draft.get("host_poscar"), "draft.host_poscar")
    if _digest(host) != draft.get("host_sha256"):
        raise InventoryError("host POSCAR changed after preparation")
    manifest = _read(_absolute_file(draft.get("pristine_manifest"), "draft.pristine_manifest"))
    doping = _read(_absolute_file(draft.get("doping_provenance"), "draft.doping_provenance"))
    if (
        _digest(Path(draft["pristine_manifest"])) != draft.get("pristine_manifest_sha256")
        or _digest(Path(draft["doping_provenance"])) != draft.get("doping_provenance_sha256")
    ):
        raise InventoryError("producer manifest or doping provenance changed after preparation")
    records = [draft["doped"], *draft["complexes"]]
    by_id = {item["structure_id"]: item for item in records}
    if len(by_id) != len(records):
        raise InventoryError("duplicate structure_id in draft")
    for record in records:
        if _digest(_absolute_file(record["poscar"], "draft candidate")) != record["sha256"]:
            raise InventoryError("candidate POSCAR changed after preparation")
    for record in draft["complexes"]:
        if (
            _digest(_absolute_file(record["provenance_path"], "draft complex provenance"))
            != record["provenance_sha256"]
        ):
            raise InventoryError("complex provenance changed after preparation")
    interstitial_candidates = []
    selected_interstitials = []
    approved_interstitial_ids: set[str] = set()
    if interstitial_enabled:
        substitution = doping.get("substitution")
        if not isinstance(substitution, dict) or not isinstance(
            substitution.get("to"), str
        ):
            raise InventoryError("doping provenance lacks its substitution dopant")
        interstitial_candidates = _validated_interstitial_candidates(
            draft, manifest, substitution["to"],
        )
        selected_interstitials, approved_interstitial_ids = _validated_interstitial_audits(
            review, interstitial_candidates,
        )
    interstitial_by_id = {
        candidate["structure_id"]: candidate for candidate in interstitial_candidates
    }
    all_by_id = {**by_id, **{
        structure_id: {
            "poscar": candidate["poscar"],
            "sha256": candidate["sha256"],
        }
        for structure_id, candidate in interstitial_by_id.items()
    }}
    states = review.get("states")
    if not isinstance(states, list) or not states:
        raise InventoryError("review.states must include an approved charge/spin state")
    selected_ids = {state.get("structure_id") for state in states if isinstance(state, dict)}
    if None in selected_ids or selected_ids - set(all_by_id):
        raise InventoryError("review names a missing candidate structure")
    if draft["doped"]["structure_id"] not in selected_ids:
        raise InventoryError("review must include the isolated dopant substitution")
    missing_interstitial_states = approved_interstitial_ids - selected_ids
    if missing_interstitial_states:
        raise InventoryError(
            "every approved interstitial must be explicitly selected by states: "
            + ", ".join(sorted(missing_interstitial_states))
        )
    unapproved_interstitial_states = (
        selected_ids & set(interstitial_by_id) - approved_interstitial_ids
    )
    if unapproved_interstitial_states:
        raise InventoryError(
            "review states select an interstitial without an approved site audit: "
            + ", ".join(sorted(unapproved_interstitial_states))
        )
    reviewed_states = []
    seen_state_ids: set[str] = set()
    for item in states:
        state = _object(item, "review.states[]")
        if (
            set(state) != {
                "structure_id", "defect_id", "geometry_id", "state_id",
                "charge_audit", "spin_audit",
            }
            or any(
                not isinstance(state[field], str) or not state[field]
                for field in ("structure_id", "defect_id", "geometry_id", "state_id")
            )
            or not isinstance(state["state_id"], str)
            or not state["state_id"]
            or state["state_id"] in seen_state_ids
        ):
            raise InventoryError("review.states[] must contain unique complete state IDs")
        seen_state_ids.add(state["state_id"])
        binding = all_by_id[state["structure_id"]]["sha256"]
        for field in ("charge_audit", "spin_audit"):
            audit = _object(state.get(field), f"review.{field}")
            if (
                audit.get("status") != "approved"
                or audit.get("structure_sha256") != binding
                or not audit.get("reviewer")
                or not audit.get("evidence_ref")
            ):
                raise InventoryError(f"review.{field} must be approved and bind the exact POSCAR")
        charge = state["charge_audit"].get("q")
        spin = state["spin_audit"].get("spin_state")
        if isinstance(charge, bool) or not isinstance(charge, int):
            raise InventoryError("charge_audit.q must be a reviewed integer")
        if not isinstance(spin, dict) or not (
            set(spin) == {"nupdown"}
            and isinstance(spin["nupdown"], int)
            and not isinstance(spin["nupdown"], bool)
            or spin == {"constraint": "unconstrained"}
        ):
            raise InventoryError("spin_audit.spin_state must set integer nupdown or unconstrained")
        reviewed_states.append(state)
    complex_audits = review.get("complex_candidate_audits", {})
    selected_complex_ids = selected_ids & set(by_id) - {draft["doped"]["structure_id"]}
    if not isinstance(complex_audits, dict) or set(complex_audits) != selected_complex_ids:
        raise InventoryError("complex candidate audits must cover exactly the selected complexes")
    selected_complexes = []
    for structure_id, audit in sorted(complex_audits.items()):
        record = by_id[structure_id]
        audit = _object(audit, f"candidate audit {structure_id}")
        if (
            audit.get("status") != "approved"
            or audit.get("structure_sha256") != record["sha256"]
            or not audit.get("reviewer")
            or not audit.get("evidence_ref")
        ):
            raise InventoryError("complex candidate audit must bind the selected POSCAR")
        selected_complexes.append({
            "structure_id": structure_id,
            "provenance_path": record["provenance_path"],
            "candidate_sha256": record["sha256"],
            "doped_structure_id": draft["doped"]["structure_id"],
            "candidate_audit": audit,
        })
    files = {
        str(host), draft["pristine_manifest"],
        manifest["pristine_supercell"]["path"],
        draft["doping_provenance"], doping["doped_supercell"]["path"],
        *(record["provenance_path"] for record in selected_complexes),
        *(all_by_id[structure_id]["poscar"] for structure_id in selected_ids),
        *(candidate["provenance_path"] for candidate in selected_interstitials),
        *(candidate["candidate_path"] for candidate in selected_interstitials),
    }
    inventory = {
        "schema": "vasp-extrinsic-structure-inventory-1",
        "pristine": {"manifest_path": draft["pristine_manifest"]},
        "doped_structures": [{
            "structure_id": draft["doped"]["structure_id"],
            "provenance_path": draft["doping_provenance"],
        }],
        "complex_structures": selected_complexes,
        **(
            {"interstitial_structures": selected_interstitials}
            if interstitial_enabled
            else {}
        ),
        "states": reviewed_states,
        "freeze": {
            "status": "frozen",
            "frozen_at_utc": datetime.now(timezone.utc).isoformat(),
            "files": [
                {"path": path, "sha256": _digest(Path(path))}
                for path in sorted(files)
            ],
        },
    }
    destination = Path(inventory_path).expanduser().resolve()
    if destination.exists():
        raise InventoryError(f"inventory already exists: {destination}")
    # The calculation consumer validates all provenance and frozen files again.
    destination.parent.mkdir(parents=True, exist_ok=True)
    _write(destination, inventory)
    return destination


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest="command", required=True)
    pre = sub.add_parser("prepare")
    pre.add_argument("recipe")
    pre.add_argument("output_dir")
    done = sub.add_parser("freeze")
    done.add_argument("draft")
    done.add_argument("review")
    done.add_argument("inventory")
    args = parser.parse_args(argv)
    try:
        result = (
            prepare(args.recipe, args.output_dir) if args.command == "prepare"
            else freeze(args.draft, args.review, args.inventory)
        )
    except (InventoryError, OSError, ValueError, KeyError, TypeError) as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 2
    print(result)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
