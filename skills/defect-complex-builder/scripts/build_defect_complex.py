#!/usr/bin/env python3
"""Rank and build small, auditable defect-complex candidate sets."""

from __future__ import annotations

import argparse
import hashlib
import itertools
import json
import math
import os
import sys
import tempfile
from copy import deepcopy
from pathlib import Path
from typing import Any

import numpy as np
from pymatgen.core import PeriodicSite, Structure
from pymatgen.io.vasp import Poscar


PROVENANCE_SCHEMA = "doping-supercell.doping-provenance/v2"
OUTPUT_SCHEMA = "defect-complex-builder.provenance/v1"
DEFAULT_WEIGHTS = {"center_distance": 0.50, "electrostatic": 0.25, "possible_bonding": 0.25}


class ConfigError(ValueError):
    """Raised when an input cannot satisfy the construction contract."""


def _load_json(path: Path) -> dict[str, Any]:
    try:
        data = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as exc:
        raise ConfigError(f"cannot read JSON {path}: {exc}") from exc
    if not isinstance(data, dict):
        raise ConfigError(f"JSON root must be an object: {path}")
    return data


def _resolve(raw: str, relative_to: Path) -> Path:
    path = Path(raw).expanduser()
    if not path.is_absolute():
        path = relative_to / path
    return path.resolve()


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _composition(structure: Structure) -> dict[str, int | float]:
    result: dict[str, int | float] = {}
    for element, amount in structure.composition.as_dict().items():
        rounded = round(float(amount))
        result[element] = int(rounded) if math.isclose(float(amount), rounded) else float(amount)
    return result


def _composition_delta(base: Structure, output: Structure) -> dict[str, int | float]:
    before = _composition(base)
    after = _composition(output)
    delta: dict[str, int | float] = {}
    for element in sorted(set(before) | set(after)):
        value = float(after.get(element, 0)) - float(before.get(element, 0))
        if not math.isclose(value, 0.0):
            rounded = round(value)
            delta[element] = int(rounded) if math.isclose(value, rounded) else value
    return delta


def _atomic_write_json(path: Path, payload: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.NamedTemporaryFile(
        mode="w", encoding="utf-8", dir=path.parent, prefix=f".{path.name}.", delete=False
    ) as handle:
        json.dump(payload, handle, indent=2, sort_keys=True)
        handle.write("\n")
        temporary = Path(handle.name)
    os.replace(temporary, path)


def _atomic_write_poscar(path: Path, structure: Structure) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.NamedTemporaryFile(dir=path.parent, prefix=f".{path.name}.", delete=False) as handle:
        temporary = Path(handle.name)
    try:
        Poscar(structure, sort_structure=False).write_file(temporary)
        os.replace(temporary, path)
    finally:
        temporary.unlink(missing_ok=True)


def _as_object(value: Any, label: str) -> dict[str, Any]:
    if not isinstance(value, dict):
        raise ConfigError(f"{label} must be an object")
    return value


def _as_nonempty_string(value: Any, label: str) -> str:
    if not isinstance(value, str) or not value.strip():
        raise ConfigError(f"{label} must be a non-empty string")
    return value.strip()


def _charge_hint(value: Any, label: str) -> int | None:
    if value is None:
        return None
    if isinstance(value, bool) or value not in {-1, 0, 1}:
        raise ConfigError(f"{label} must be -1, 0, 1, or null")
    return int(value)


def _positive_float(value: Any, label: str) -> float:
    if isinstance(value, bool):
        raise ConfigError(f"{label} must be a positive number")
    try:
        number = float(value)
    except (TypeError, ValueError) as exc:
        raise ConfigError(f"{label} must be a positive number") from exc
    if not math.isfinite(number) or number <= 0:
        raise ConfigError(f"{label} must be a positive finite number")
    return number


def _site_index(value: Any, size: int, label: str) -> int:
    if isinstance(value, bool) or not isinstance(value, int):
        raise ConfigError(f"{label} must be an integer")
    if value < 1 or value > size:
        raise ConfigError(f"{label}={value} is outside 1..{size}")
    return value - 1


def _distance(structure: Structure, frac_a: np.ndarray, frac_b: np.ndarray) -> float:
    return float(structure.lattice.get_distance_and_image(frac_a, frac_b)[0])


def _load_source(config: dict[str, Any], config_dir: Path) -> dict[str, Any]:
    provenance_raw = _as_nonempty_string(config.get("doping_provenance"), "doping_provenance")
    provenance_path = _resolve(provenance_raw, config_dir)
    provenance = _load_json(provenance_path)
    if provenance.get("schema") != PROVENANCE_SCHEMA:
        raise ConfigError(f"doping provenance schema must be {PROVENANCE_SCHEMA}")
    if provenance.get("artifact_type") != "doped_supercell":
        raise ConfigError("doping provenance artifact_type must be doped_supercell")
    if provenance.get("status") != "passed" or provenance.get("stage") != "dope":
        raise ConfigError("doping provenance must be a passed dope-stage artifact")
    validation = provenance.get("validation")
    if not isinstance(validation, dict) or validation.get("passed") is not True:
        raise ConfigError("doping provenance validation.passed must be true")

    doped = _as_object(provenance.get("doped_supercell"), "doped_supercell")
    expected_hash = _as_nonempty_string(doped.get("sha256"), "doped_supercell.sha256")
    default_base = _as_nonempty_string(doped.get("path"), "doped_supercell.path")
    if config.get("base_poscar") is None:
        base_path = _resolve(default_base, provenance_path.parent)
    else:
        base_path = _resolve(
            _as_nonempty_string(config.get("base_poscar"), "base_poscar"), config_dir
        )
    if not base_path.is_file():
        raise ConfigError(f"base POSCAR does not exist: {base_path}")
    actual_hash = _sha256(base_path)
    if actual_hash != expected_hash:
        raise ConfigError("base POSCAR SHA256 does not match doping provenance")

    poscar = Poscar.from_file(base_path)
    structure = poscar.structure
    if doped.get("atoms") != len(structure):
        raise ConfigError("base atom count does not match doping provenance")
    if doped.get("formula") != _composition(structure):
        raise ConfigError("base composition does not match doping provenance")

    pristine_structure = None
    pristine_source = None
    pristine_record = provenance.get("pristine_supercell")
    if pristine_record is not None:
        pristine = _as_object(pristine_record, "pristine_supercell")
        pristine_path = _resolve(
            _as_nonempty_string(pristine.get("path"), "pristine_supercell.path"),
            provenance_path.parent,
        )
        pristine_hash = _as_nonempty_string(
            pristine.get("sha256"), "pristine_supercell.sha256"
        )
        if not pristine_path.is_file() or _sha256(pristine_path) != pristine_hash:
            raise ConfigError("pristine supercell POSCAR SHA256 does not match doping provenance")
        pristine_structure = Poscar.from_file(pristine_path).structure
        if pristine.get("atoms") != len(pristine_structure) or pristine.get("formula") != _composition(pristine_structure):
            raise ConfigError("pristine supercell count or composition disagrees with doping provenance")
        if not np.allclose(
            pristine_structure.lattice.matrix, structure.lattice.matrix, atol=1e-5, rtol=1e-7
        ):
            raise ConfigError("pristine and doped supercell lattice disagree")
        pristine_source = {
            "path": str(pristine_path),
            "sha256": pristine_hash,
            "composition": _composition(pristine_structure),
        }

    substitution = _as_object(provenance.get("substitution"), "substitution")
    default_center = substitution.get("output_site_1based")
    center_index = _site_index(
        config.get("center_site_1based", default_center),
        len(structure),
        "center_site_1based",
    )
    return {
        "provenance_path": provenance_path,
        "provenance": provenance,
        "base_path": base_path,
        "base_hash": actual_hash,
        "poscar": poscar,
        "structure": structure,
        "pristine_structure": pristine_structure,
        "pristine_source": pristine_source,
        "center_index": center_index,
    }


def _geometry_source(value: Any, label: str) -> dict[str, str]:
    source = _as_object(value, label)
    return {
        "method": _as_nonempty_string(source.get("method"), f"{label}.method"),
        "details": _as_nonempty_string(source.get("details"), f"{label}.details"),
    }


def _bonding_score(
    candidate: dict[str, Any],
    operation: dict[str, Any],
    structure: Structure,
    frac_coords: np.ndarray,
) -> dict[str, Any]:
    raw = candidate.get("bonding", operation.get("bonding"))
    if raw is None:
        return {
            "status": "unknown",
            "basis": "not_reported",
            "score": 0.5,
            "targets": [],
        }
    bonding = _as_object(raw, f"{operation['id']}.bonding")
    target_values = bonding.get("target_sites_1based")
    if not isinstance(target_values, list) or not target_values:
        raise ConfigError(f"{operation['id']}.bonding.target_sites_1based must be non-empty")
    ideal = _positive_float(bonding.get("ideal_distance_A"), "bonding.ideal_distance_A")
    tolerance = _positive_float(bonding.get("tolerance_A"), "bonding.tolerance_A")
    basis = _as_nonempty_string(bonding.get("basis"), "bonding.basis")
    distances = []
    for raw_index in target_values:
        index = _site_index(raw_index, len(structure), "bonding target site")
        distance = _distance(structure, frac_coords, structure[index].frac_coords)
        distances.append(
            {
                "site_1based": index + 1,
                "species": structure[index].specie.symbol,
                "distance_A": distance,
                "absolute_deviation_A": abs(distance - ideal),
            }
        )
    best = min(distances, key=lambda item: (item["absolute_deviation_A"], item["site_1based"]))
    score = min(best["absolute_deviation_A"] / tolerance, 2.0) / 2.0
    return {
        "status": "reported",
        "basis": basis,
        "ideal_distance_A": ideal,
        "tolerance_A": tolerance,
        "best_target_site_1based": best["site_1based"],
        "score": score,
        "targets": distances,
    }


def _electrostatic_score(
    center_hint: int | None, candidate_hint: int | None, basis: str
) -> dict[str, Any]:
    if center_hint is None or candidate_hint is None:
        assessment, score = "unknown", 0.5
    elif center_hint == 0 or candidate_hint == 0:
        assessment, score = "neutral", 0.5
    elif center_hint == -candidate_hint:
        assessment, score = "favorable_opposite_sign", 0.0
    else:
        assessment, score = "unfavorable_same_sign", 1.0
    return {
        "assessment": assessment,
        "center_charge_hint": center_hint,
        "component_charge_hint": candidate_hint,
        "basis": basis,
        "score": score,
    }


def _weights(config: dict[str, Any]) -> tuple[dict[str, float], float]:
    ranking = config.get("ranking", {})
    ranking = _as_object(ranking, "ranking")
    raw_weights = ranking.get("weights", DEFAULT_WEIGHTS)
    raw_weights = _as_object(raw_weights, "ranking.weights")
    if set(raw_weights) != set(DEFAULT_WEIGHTS):
        raise ConfigError(
            "ranking.weights must contain center_distance, electrostatic, and possible_bonding"
        )
    weights = {key: _positive_float(raw_weights[key], f"ranking.weights.{key}") for key in raw_weights}
    total = sum(weights.values())
    weights = {key: value / total for key, value in weights.items()}
    scale = _positive_float(ranking.get("distance_scale_A", 5.0), "ranking.distance_scale_A")
    return weights, scale


def _rank_candidate(
    candidate: dict[str, Any],
    operation: dict[str, Any],
    structure: Structure,
    center_frac: np.ndarray,
    center_hint: int | None,
    weights: dict[str, float],
    distance_scale: float,
) -> dict[str, Any]:
    result = deepcopy(candidate)
    frac = np.array(candidate["fractional_coords"], dtype=float) % 1.0
    center_distance = _distance(structure, center_frac, frac)
    distance_score = min(center_distance / distance_scale, 2.0) / 2.0
    candidate_hint = _charge_hint(
        candidate.get("charge_hint", operation.get("charge_hint")),
        f"{operation['id']} candidate charge_hint",
    )
    basis = candidate.get("electrostatic_basis", operation.get("electrostatic_basis"))
    if basis is None:
        basis = "not_reported"
    else:
        basis = _as_nonempty_string(basis, f"{operation['id']}.electrostatic_basis")
    electrostatic = _electrostatic_score(center_hint, candidate_hint, basis)
    bonding = _bonding_score(candidate, operation, structure, frac)
    components = {
        "center_distance": distance_score,
        "electrostatic": electrostatic["score"],
        "possible_bonding": bonding["score"],
    }
    total = sum(weights[key] * components[key] for key in weights)
    result.update(
        {
            "fractional_coords": [float(value) for value in frac],
            "center_distance_A": center_distance,
            "electrostatic": electrostatic,
            "possible_bonding": bonding,
            "score_components": components,
            "total_score": total,
        }
    )
    return result


def _auto_interstitial_candidates(
    operation: dict[str, Any],
    structure: Structure,
    center_index: int,
    minimum_separation: float,
) -> list[dict[str, Any]]:
    """Place a small number of auditable points on a ligand-derived shell."""
    spec = _as_object(operation.get("generate"), f"{operation['id']}.generate")
    reference_species = _as_nonempty_string(
        spec.get("reference_species", operation["species"]),
        f"{operation['id']}.generate.reference_species",
    )
    neighbor_count = spec.get("reference_neighbor_count", 4)
    count = spec.get("max_candidates", 2)
    if isinstance(neighbor_count, bool) or not isinstance(neighbor_count, int) or neighbor_count < 1:
        raise ConfigError("generate.reference_neighbor_count must be a positive integer")
    if isinstance(count, bool) or not isinstance(count, int) or not 1 <= count <= 2:
        raise ConfigError("generate.max_candidates must be 1 or 2")
    shortening = spec.get("shortening_fraction", 0.05)
    if isinstance(shortening, bool) or not isinstance(shortening, int | float) or not 0 <= shortening < 0.25:
        raise ConfigError("generate.shortening_fraction must be in [0, 0.25)")
    clearance = _positive_float(
        spec.get("minimum_clearance_A", max(minimum_separation, 1.2)),
        "generate.minimum_clearance_A",
    )
    center = structure[center_index]
    neighbors = sorted(
        (
            _distance(structure, center.frac_coords, site.frac_coords)
            for index, site in enumerate(structure)
            if index != center_index and site.specie.symbol == reference_species
        )
    )
    if len(neighbors) < neighbor_count:
        raise ConfigError(
            f"{operation['id']} needs {neighbor_count} measured {reference_species} "
            f"neighbors around the dopant; found {len(neighbors)}"
        )
    reference_distances = neighbors[:neighbor_count]
    mean_distance = sum(reference_distances) / neighbor_count
    target_distance = mean_distance * (1.0 - float(shortening))
    if target_distance <= clearance:
        raise ConfigError(
            f"{operation['id']} target distance {target_distance:.4f} A "
            "does not exceed the minimum clearance"
        )

    # Deterministic Fibonacci sphere; optimize empty space at a fixed, measured
    # center distance instead of pushing the inserted atom onto an existing ligand.
    directions = []
    golden_angle = math.pi * (3.0 - math.sqrt(5.0))
    for index in range(128):
        z = 1.0 - 2.0 * (index + 0.5) / 128.0
        radius = math.sqrt(1.0 - z * z)
        theta = index * golden_angle
        directions.append(np.array([radius * math.cos(theta), radius * math.sin(theta), z]))
    proposals = []
    for index, direction in enumerate(directions):
        cart = center.coords + direction * target_distance
        frac = structure.lattice.get_fractional_coords(cart) % 1.0
        actual_center_distance = _distance(structure, center.frac_coords, frac)
        # In a very small supercell, a periodic image can be closer than the
        # intended radius: reject it rather than silently changing the model.
        if abs(actual_center_distance - target_distance) > 1e-5:
            continue
        nearest = min(
            _distance(structure, frac, site.frac_coords)
            for j, site in enumerate(structure) if j != center_index
        )
        if nearest + 1e-8 >= clearance:
            proposals.append((nearest, index, frac))
    proposals.sort(key=lambda item: (-item[0], item[1]))
    chosen: list[tuple[float, int, np.ndarray]] = []
    for proposal in proposals:
        if all(
            _distance(structure, proposal[2], other[2]) >= max(clearance, target_distance * 0.5)
            for other in chosen
        ):
            chosen.append(proposal)
        if len(chosen) == count:
            break
    if not chosen:
        raise ConfigError(
            f"{operation['id']} found no collision-free interstitial near "
            f"the dopant at {target_distance:.4f} A; review the measured shell or relax constraints"
        )
    return [
        {
            "id": f"{operation['id']}:auto-{index + 1}",
            "fractional_coords": [float(value) for value in frac],
            "geometry_source": {
                "method": "measured-dopant-ligand-shell-void-search",
                "details": (
                    f"reference {reference_species} distances (A): "
                    f"{[round(value, 6) for value in reference_distances]}; "
                    f"mean={mean_distance:.6f} A, shortening={shortening}, "
                    f"target={target_distance:.6f} A; Fibonacci-128 direction {direction_id}; "
                    f"nearest noncenter atom={nearest:.6f} A"
                ),
            },
            "bonding": {
                "target_sites_1based": [center_index + 1],
                "ideal_distance_A": target_distance,
                "tolerance_A": max(0.1, target_distance * 0.1),
                "basis": f"measured {reference_species} shell around substituted center",
            },
            "auto_geometry_audit": {
                "reference_species": reference_species,
                "neighbor_distances_A": reference_distances,
                "mean_distance_A": mean_distance,
                "shortening_fraction": float(shortening),
                "target_distance_A": target_distance,
                "nearest_noncenter_A": nearest,
                "minimum_clearance_A": clearance,
                "direction_index": direction_id,
            },
        }
        for index, (nearest, direction_id, frac) in enumerate(chosen)
    ]


def _local_environment(structure: Structure, index: int) -> dict[str, tuple[float, ...]]:
    result: dict[str, list[float]] = {}
    for other_index, site in enumerate(structure):
        if other_index == index:
            continue
        result.setdefault(site.specie.symbol, []).append(
            _distance(structure, structure[index].frac_coords, site.frac_coords)
        )
    return {species: tuple(sorted(distances)[:4]) for species, distances in result.items()}


def _near_environment(
    left: dict[str, tuple[float, ...]],
    right: dict[str, tuple[float, ...]],
    tolerance: float,
) -> bool:
    return left.keys() == right.keys() and all(
        len(left[species]) == len(right[species])
        and all(abs(a - b) <= tolerance for a, b in zip(left[species], right[species], strict=True))
        for species in left
    )


def _operation_candidates(
    raw_operation: Any,
    operation_number: int,
    structure: Structure,
    center_index: int,
    center_hint: int | None,
    weights: dict[str, float],
    distance_scale: float,
    minimum_separation: float,
    environment_tolerance: float,
    nearest_shell_window: float,
) -> dict[str, Any]:
    operation = _as_object(raw_operation, f"operations[{operation_number}]")
    operation_id = _as_nonempty_string(
        operation.get("id", f"component_{operation_number + 1}"),
        f"operations[{operation_number}].id",
    )
    kind = operation.get("kind")
    if kind not in {"vacancy", "interstitial", "antisite"}:
        raise ConfigError(f"{operation_id}.kind must be vacancy, interstitial or antisite")
    species = _as_nonempty_string(operation.get("species"), f"{operation_id}.species")
    normalized = dict(operation)
    normalized.update({"id": operation_id, "kind": kind, "species": species})
    normalized["charge_hint"] = _charge_hint(operation.get("charge_hint"), f"{operation_id}.charge_hint")

    candidates: list[dict[str, Any]] = []
    if kind in {"vacancy", "antisite"}:
        if kind == "antisite":
            from_species = _as_nonempty_string(
                operation.get("from_species"), f"{operation_id}.from_species"
            )
            if from_species == species or species not in structure.symbol_set:
                raise ConfigError(
                    f"{operation_id} antisite must replace a different host species "
                    "with a species already present in the doped supercell"
                )
            normalized["from_species"] = from_species
        else:
            from_species = species
        requested = operation.get("candidate_sites_1based")
        if requested is None:
            indices = [
                index
                for index, site in enumerate(structure)
                if site.specie.symbol == from_species and index != center_index
            ]
        else:
            if not isinstance(requested, list) or not requested:
                raise ConfigError(f"{operation_id}.candidate_sites_1based must be non-empty")
            indices = [
                _site_index(value, len(structure), f"{operation_id}.candidate_sites_1based")
                for value in requested
            ]
        if len(indices) != len(set(indices)):
            raise ConfigError(f"{operation_id} contains duplicate source candidate sites")
        for index in indices:
            if index == center_index:
                raise ConfigError(f"{operation_id} cannot remove the substitution center or replace it")
            if structure[index].specie.symbol != from_species:
                raise ConfigError(
                    f"{operation_id} site {index + 1} is {structure[index].specie.symbol}, not {from_species}"
                )
            candidates.append(
                {
                    "id": f"{operation_id}:site-{index + 1}",
                    "kind": kind,
                    "species": species,
                    **({"from_species": from_species} if kind == "antisite" else {}),
                    "source_site_1based": index + 1,
                    "fractional_coords": [float(value) for value in structure[index].frac_coords],
                    "geometry_source": {
                        "method": "existing-supercell-site",
                        "details": (
                            f"{'remove' if kind == 'vacancy' else 'replace'} source site "
                            f"{index + 1} from accepted doped supercell"
                        ),
                    },
                }
            )
    else:
        if "generate" in operation and "candidates" in operation:
            raise ConfigError(f"{operation_id} accepts either generate or candidates, not both")
        raw_candidates = (
            _auto_interstitial_candidates(operation, structure, center_index, minimum_separation)
            if "generate" in operation else operation.get("candidates")
        )
        if not isinstance(raw_candidates, list) or not raw_candidates:
            raise ConfigError(f"{operation_id}.candidates must be a non-empty list")
        for candidate_number, raw_candidate in enumerate(raw_candidates):
            candidate = _as_object(raw_candidate, f"{operation_id}.candidates[{candidate_number}]")
            candidate_id = _as_nonempty_string(
                candidate.get("id", f"{operation_id}:candidate-{candidate_number + 1}"),
                f"{operation_id}.candidates[{candidate_number}].id",
            )
            has_fractional = "fractional_coords" in candidate
            has_cartesian = "cartesian_coords_A" in candidate
            if has_fractional == has_cartesian:
                raise ConfigError(
                    f"{candidate_id} must provide exactly one of fractional_coords or cartesian_coords_A"
                )
            key = "fractional_coords" if has_fractional else "cartesian_coords_A"
            coords = candidate[key]
            if (
                not isinstance(coords, list)
                or len(coords) != 3
                or any(isinstance(value, bool) for value in coords)
            ):
                raise ConfigError(f"{candidate_id}.{key} must be a three-number list")
            try:
                vector = np.array(coords, dtype=float)
            except (TypeError, ValueError) as exc:
                raise ConfigError(f"{candidate_id}.{key} must be numeric") from exc
            if not np.all(np.isfinite(vector)):
                raise ConfigError(f"{candidate_id}.{key} must contain finite numbers")
            frac = vector if has_fractional else structure.lattice.get_fractional_coords(vector)
            normalized_candidate = dict(candidate)
            normalized_candidate.update(
                {
                    "id": candidate_id,
                    "kind": kind,
                    "species": species,
                    "fractional_coords": [float(value) for value in frac % 1.0],
                    "geometry_source": _geometry_source(
                        candidate.get("geometry_source"), f"{candidate_id}.geometry_source"
                    ),
                }
            )
            candidates.append(normalized_candidate)

    if not candidates:
        raise ConfigError(f"{operation_id} has no eligible candidates")
    ids = [candidate["id"] for candidate in candidates]
    if len(ids) != len(set(ids)):
        raise ConfigError(f"{operation_id} candidate IDs must be unique")
    ranked = [
        _rank_candidate(
            candidate,
            normalized,
            structure,
            structure[center_index].frac_coords,
            center_hint,
            weights,
            distance_scale,
        )
        for candidate in candidates
    ]
    ranked.sort(key=lambda item: (item["total_score"], item["id"]))
    normalized["candidates"] = ranked
    selected: list[dict[str, Any]] = []
    deferred: list[dict[str, Any]] = []
    # Only auto-enumerated host sites are thinned; explicitly supplied sites
    # remain explicit user hypotheses even if their unrelaxed fingerprints match.
    if kind in {"vacancy", "antisite"} and operation.get("candidate_sites_1based") is None:
        for candidate in ranked:
            fingerprint = _local_environment(
                structure, candidate["source_site_1based"] - 1
            )
            candidate["environment_fingerprint_A"] = fingerprint
            representative = next(
                (
                    other for other in selected
                    if abs(other["center_distance_A"] - candidate["center_distance_A"])
                    <= environment_tolerance
                    and _near_environment(
                        fingerprint,
                        other["environment_fingerprint_A"],
                        environment_tolerance,
                    )
                ),
                None,
            )
            if representative is None:
                selected.append(candidate)
            else:
                deferred.append({
                    "candidate_id": candidate["id"],
                    "representative_id": representative["id"],
                    "reason": "deferred_similar_unrelaxed_environment",
                    "tolerance_A": environment_tolerance,
                    "center_distance_difference_A": abs(
                        candidate["center_distance_A"] - representative["center_distance_A"]
                    ),
                    "max_neighbor_distance_difference_A": max(
                        abs(a - b)
                        for species in fingerprint
                        for a, b in zip(
                            fingerprint[species],
                            representative["environment_fingerprint_A"][species],
                            strict=True,
                        )
                    ),
                })
    else:
        selected = ranked
    if kind in {"vacancy", "antisite"} and operation.get("candidate_sites_1based") is None:
        nearest = min(item["center_distance_A"] for item in ranked)
        near_shell = []
        for candidate in selected:
            if candidate["center_distance_A"] <= nearest + nearest_shell_window + 1e-8:
                near_shell.append(candidate)
            else:
                deferred.append({
                    "candidate_id": candidate["id"],
                    "representative_id": None,
                    "reason": "deferred_distinct_distant_unrepresented",
                    "center_distance_A": candidate["center_distance_A"],
                    "nearest_center_distance_A": nearest,
                    "nearest_shell_window_A": nearest_shell_window,
                })
        selected = near_shell
    normalized["selected_candidate_ids"] = [item["id"] for item in selected]
    normalized["deferred_candidates"] = deferred
    normalized["_selected_candidates"] = selected
    return normalized


def _selection(config: dict[str, Any]) -> dict[str, Any]:
    raw = _as_object(config.get("selection", {}), "selection")
    max_structures = raw.get("max_structures", 3)
    max_combinations = raw.get("max_combinations", 500)
    if isinstance(max_structures, bool) or not isinstance(max_structures, int):
        raise ConfigError("selection.max_structures must be an integer")
    if max_structures < 1 or max_structures > 12:
        raise ConfigError("selection.max_structures must be within 1..12")
    if isinstance(max_combinations, bool) or not isinstance(max_combinations, int):
        raise ConfigError("selection.max_combinations must be an integer")
    if max_combinations < 1:
        raise ConfigError("selection.max_combinations must be positive")
    minimum = _positive_float(
        raw.get("minimum_separation_A", 0.7), "selection.minimum_separation_A"
    )
    tolerance = _positive_float(
        raw.get("environment_tolerance_A", 0.12), "selection.environment_tolerance_A"
    )
    window_raw = raw.get("nearest_shell_window_A", 0.35)
    if isinstance(window_raw, bool) or not isinstance(window_raw, int | float):
        raise ConfigError("selection.nearest_shell_window_A must be a nonnegative finite number")
    window = float(window_raw)
    if not math.isfinite(window) or window < 0:
        raise ConfigError("selection.nearest_shell_window_A must be a nonnegative finite number")
    return {
        "max_structures": max_structures,
        "max_combinations": max_combinations,
        "minimum_separation_A": minimum,
        "environment_tolerance_A": tolerance,
        "nearest_shell_window_A": window,
    }


def _combination_rejection(
    combination: tuple[dict[str, Any], ...],
    structure: Structure,
    minimum_separation: float,
) -> str | None:
    source_indices = [
        int(candidate["source_site_1based"]) - 1
        for candidate in combination
        if candidate["kind"] in {"vacancy", "antisite"}
    ]
    if len(source_indices) != len(set(source_indices)):
        return "duplicate vacancy/antisite source site"
    vacancy_indices = {
        int(candidate["source_site_1based"]) - 1
        for candidate in combination if candidate["kind"] == "vacancy"
    }
    retained = [index for index in range(len(structure)) if index not in vacancy_indices]
    interstitials = [candidate for candidate in combination if candidate["kind"] == "interstitial"]
    for candidate in interstitials:
        frac = np.array(candidate["fractional_coords"], dtype=float)
        nearest = min(
            (_distance(structure, frac, structure[index].frac_coords) for index in retained),
            default=math.inf,
        )
        if nearest < minimum_separation:
            return (
                f"interstitial {candidate['id']} is {nearest:.6f} A from a retained atom, "
                f"below {minimum_separation:.6f} A"
            )
    for first, second in itertools.combinations(interstitials, 2):
        separation = _distance(
            structure,
            np.array(first["fractional_coords"], dtype=float),
            np.array(second["fractional_coords"], dtype=float),
        )
        if separation < minimum_separation:
            return (
                f"interstitials {first['id']} and {second['id']} are {separation:.6f} A apart, "
                f"below {minimum_separation:.6f} A"
            )
    return None


def rank_config(config_path: str | Path) -> dict[str, Any]:
    path = Path(config_path).expanduser().resolve()
    config = _load_json(path)
    source = _load_source(config, path.parent)
    operations_raw = config.get("operations")
    if not isinstance(operations_raw, list) or not 1 <= len(operations_raw) <= 3:
        raise ConfigError("operations must contain 1..3 vacancy/interstitial/antisite components")
    center_hint = _charge_hint(config.get("center_charge_hint"), "center_charge_hint")
    weights, distance_scale = _weights(config)
    selection = _selection(config)
    operations = [
        _operation_candidates(
            operation,
            number,
            source["structure"],
            source["center_index"],
            center_hint,
            weights,
            distance_scale,
            selection["minimum_separation_A"],
            selection["environment_tolerance_A"],
            selection["nearest_shell_window_A"],
        )
        for number, operation in enumerate(operations_raw)
    ]
    operation_ids = [operation["id"] for operation in operations]
    if len(operation_ids) != len(set(operation_ids)):
        raise ConfigError("operation IDs must be unique")

    combination_count = math.prod(len(operation["_selected_candidates"]) for operation in operations)
    if combination_count > selection["max_combinations"]:
        raise ConfigError(
            f"candidate product {combination_count} exceeds selection.max_combinations="
            f"{selection['max_combinations']}; narrow the candidate sets"
        )
    ranked: list[dict[str, Any]] = []
    rejected: list[dict[str, Any]] = []
    for combination in itertools.product(*(operation["_selected_candidates"] for operation in operations)):
        ids = [candidate["id"] for candidate in combination]
        reason = _combination_rejection(
            combination, source["structure"], selection["minimum_separation_A"]
        )
        if reason:
            rejected.append({"candidate_ids": ids, "reason": reason})
            continue
        score = sum(float(candidate["total_score"]) for candidate in combination) / len(combination)
        component_means = {
            key: sum(float(candidate["score_components"][key]) for candidate in combination)
            / len(combination)
            for key in DEFAULT_WEIGHTS
        }
        ranked.append(
            {
                "candidate_ids": ids,
                "total_score": score,
                "mean_score_components": component_means,
                "components": [deepcopy(candidate) for candidate in combination],
            }
        )
    ranked.sort(key=lambda item: (item["total_score"], tuple(item["candidate_ids"])))
    if not ranked:
        raise ConfigError("no candidate combination passed geometry checks")
    for rank, item in enumerate(ranked, start=1):
        item["rank"] = rank

    center = source["structure"][source["center_index"]]
    return {
        "schema": "defect-complex-builder.ranking/v1",
        "status": "passed",
        "config": str(path),
        "source": {
            "doping_provenance": str(source["provenance_path"]),
            "doping_provenance_sha256": _sha256(source["provenance_path"]),
            "base_poscar": str(source["base_path"]),
            "base_poscar_sha256": source["base_hash"],
        },
        "center": {
            "site_1based": source["center_index"] + 1,
            "species": center.specie.symbol,
            "fractional_coords": [float(value) for value in center.frac_coords],
            "charge_hint": center_hint,
        },
        "ranking": {
            "lower_is_better": True,
            "weights": weights,
            "distance_scale_A": distance_scale,
            "unknown_dimension_score": 0.5,
            "formula": "weighted mean of center_distance, electrostatic, possible_bonding",
        },
        "selection": selection,
        "operations": [
            {key: value for key, value in operation.items() if key != "_selected_candidates"}
            for operation in operations
        ],
        "candidate_combination_count": combination_count,
        "valid_combination_count": len(ranked),
        "rejected_combinations": rejected,
        "ranked_combinations": ranked,
    }


def _species_order(poscar: Poscar, components: list[dict[str, Any]]) -> list[str]:
    order = list(poscar.site_symbols)
    for component in components:
        species = component["species"]
        if component["kind"] in {"interstitial", "antisite"} and species not in order:
            order.append(species)
    return order


def _build_structure(
    base: Structure,
    poscar: Poscar,
    components: list[dict[str, Any]],
) -> tuple[Structure, list[dict[str, Any]], list[dict[str, Any]], list[str]]:
    vacancies = {
        int(component["source_site_1based"]) - 1
        for component in components
        if component["kind"] == "vacancy"
    }
    antisites = {
        int(component["source_site_1based"]) - 1: component
        for component in components if component["kind"] == "antisite"
    }
    entries: list[dict[str, Any]] = []
    for index, site in enumerate(base):
        if index not in vacancies:
            antisite = antisites.get(index)
            species = antisite["species"] if antisite else site.specie.symbol
            entries.append(
                {
                    "site": (
                        PeriodicSite(species, site.frac_coords, base.lattice, to_unit_cell=True)
                        if antisite else site
                    ),
                    "species": species,
                    "source_site_1based": index + 1,
                    "component_id": antisite["id"] if antisite else None,
                }
            )
    for component in components:
        if component["kind"] == "interstitial":
            entries.append(
                {
                    "site": PeriodicSite(
                        component["species"],
                        component["fractional_coords"],
                        base.lattice,
                        to_unit_cell=True,
                    ),
                    "species": component["species"],
                    "source_site_1based": None,
                    "component_id": component["id"],
                }
            )
    order = _species_order(poscar, components)
    order_index = {species: index for index, species in enumerate(order)}
    entries.sort(
        key=lambda item: (
            order_index[item["species"]],
            item["source_site_1based"] is None,
            item["source_site_1based"] or 0,
            item["component_id"] or "",
        )
    )
    output = Structure.from_sites([entry["site"] for entry in entries], charge=base.charge)
    mapping = []
    component_output_sites: dict[str, int] = {}
    for output_index, entry in enumerate(entries, start=1):
        mapping.append(
            {
                "output_site_1based": output_index,
                "species": entry["species"],
                "source_site_1based": entry["source_site_1based"],
                "component_id": entry["component_id"],
            }
        )
        if entry["component_id"] is not None:
            component_output_sites[entry["component_id"]] = output_index
    built_components = []
    for component in components:
        record = deepcopy(component)
        if component["kind"] in {"interstitial", "antisite"}:
            record["output_site_1based"] = component_output_sites[component["id"]]
        else:
            record["output_site_1based"] = None
        built_components.append(record)
    return output, mapping, built_components, order


def build_config(config_path: str | Path) -> dict[str, Any]:
    path = Path(config_path).expanduser().resolve()
    config = _load_json(path)
    ranking = rank_config(path)
    source = _load_source(config, path.parent)
    selection = ranking["selection"]
    selected = ranking["ranked_combinations"][: selection["max_structures"]]
    output_dir = _resolve(str(config.get("output_dir", "defect_complexes")), path.parent)
    provenance_path = _resolve(
        str(config.get("provenance", output_dir / "provenance.json")), path.parent
    )
    overwrite = config.get("overwrite", False)
    if not isinstance(overwrite, bool):
        raise ConfigError("overwrite must be true or false")
    output_paths = [output_dir / f"POSCAR_complex_{rank:03d}" for rank in range(1, len(selected) + 1)]
    protected = [provenance_path, *output_paths]
    if not overwrite:
        existing = [str(target) for target in protected if target.exists()]
        if existing:
            raise ConfigError(f"refusing to overwrite existing output(s): {', '.join(existing)}")
    if source["base_path"] in protected or source["provenance_path"] in protected or path in protected:
        raise ConfigError("output paths must be distinct from config and source artifacts")

    outputs = []
    for output_path, combination in zip(output_paths, selected, strict=True):
        structure, mapping, components, species_order = _build_structure(
            source["structure"], source["poscar"], combination["components"]
        )
        _atomic_write_poscar(output_path, structure)
        written = Poscar.from_file(output_path)
        if len(written.structure) != len(mapping):
            raise ConfigError(f"written atom count mismatch: {output_path}")
        if written.site_symbols != species_order:
            raise ConfigError(f"written species order mismatch: {output_path}")
        expected_delta = _composition_delta(source["structure"], structure)
        written_delta = _composition_delta(source["structure"], written.structure)
        if written_delta != expected_delta:
            raise ConfigError(f"written composition delta mismatch: {output_path}")
        outputs.append(
            {
                "rank": combination["rank"],
                "total_score": combination["total_score"],
                "mean_score_components": combination["mean_score_components"],
                "candidate_ids": combination["candidate_ids"],
                "path": str(output_path),
                "sha256": _sha256(output_path),
                "atoms": len(written.structure),
                "composition": _composition(written.structure),
                "composition_delta_from_doped_base": written_delta,
                "composition_delta_from_pristine": (
                    _composition_delta(source["pristine_structure"], written.structure)
                    if source["pristine_structure"] is not None else None
                ),
                "species_order": written.site_symbols,
                "components": components,
                "site_mapping": mapping,
                "validation": {
                    "passed": True,
                    "atom_count": True,
                    "composition_delta": True,
                    "species_order": True,
                    "minimum_separation_A": selection["minimum_separation_A"],
                },
            }
        )

    payload = {
        "schema": OUTPUT_SCHEMA,
        "artifact_type": "defect_complex_candidate_set",
        "status": "passed",
        "config": str(path),
        "source": {
            "doping_provenance": str(source["provenance_path"]),
            "doping_provenance_sha256": _sha256(source["provenance_path"]),
            "doping_provenance_schema": PROVENANCE_SCHEMA,
            "base_poscar": str(source["base_path"]),
            "base_poscar_sha256": source["base_hash"],
            "pristine_supercell": source["pristine_source"],
            "atoms": len(source["structure"]),
            "composition": _composition(source["structure"]),
            "species_order": source["poscar"].site_symbols,
        },
        "center": ranking["center"],
        "ranking_contract": ranking["ranking"],
        "selection": selection,
        "candidate_inventory": ranking["operations"],
        "candidate_combination_count": ranking["candidate_combination_count"],
        "valid_combination_count": ranking["valid_combination_count"],
        "rejected_combinations": ranking["rejected_combinations"],
        "selected_count": len(outputs),
        "unsearched_scope": {
            "exhaustive_search_performed": False,
            "structure_relaxation_performed": False,
            "note": "Only the requested top-ranked representative starting geometries were written.",
        },
        "outputs": outputs,
        "provenance": str(provenance_path),
    }
    _atomic_write_json(provenance_path, payload)
    return payload


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        description="Rank and build 1-3-component defect complexes from an accepted doped supercell"
    )
    subparsers = parser.add_subparsers(dest="command", required=True)
    for command in ("rank", "build"):
        child = subparsers.add_parser(command)
        child.add_argument("config", help="JSON configuration path")
    args = parser.parse_args(argv)
    try:
        payload = rank_config(args.config) if args.command == "rank" else build_config(args.config)
    except ConfigError as exc:
        print(f"error: {exc}", file=sys.stderr)
        return 2
    print(json.dumps(payload, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
