#!/usr/bin/env python3
"""Assemble and plot audited formation-energy inputs at selected phase points."""

from __future__ import annotations

import argparse
import importlib.util
import json
import math
import re
import sys
from copy import deepcopy
from pathlib import Path
from typing import Any, Sequence


SCHEMA_VERSION = 1
_TOL = 1e-8


def _load_formation_module():
    script = Path(__file__).with_name("defect_formation_energy.py")
    spec = importlib.util.spec_from_file_location("defect_formation_energy", script)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"cannot load formation-energy processor: {script}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


formation = _load_formation_module()


def _fail(message: str) -> None:
    raise ValueError(message)


def _text(value: Any, label: str) -> str:
    if not isinstance(value, str) or not value.strip():
        _fail(f"{label} must be a non-empty string")
    return value.strip()


def _number(value: Any, label: str) -> float:
    if isinstance(value, bool) or not isinstance(value, int | float):
        _fail(f"{label} must be a finite number")
    result = float(value)
    if not math.isfinite(result):
        _fail(f"{label} must be a finite number")
    return result


def _source_ref(
    raw: Any,
    label: str,
    *,
    expected_task_id: str | None = None,
) -> dict[str, str]:
    if not isinstance(raw, dict):
        _fail(f"{label} must identify its source file and locator")
    result = {
        "file": _text(raw.get("file"), f"{label}.file"),
        "locator": _text(raw.get("locator"), f"{label}.locator"),
    }
    task_id = raw.get("task_id")
    if expected_task_id is not None:
        task_id = _text(task_id, f"{label}.task_id")
        if task_id != expected_task_id:
            _fail(f"{label}.task_id does not match its accepted task")
        result["task_id"] = task_id
    elif task_id is not None:
        result["task_id"] = _text(task_id, f"{label}.task_id")
    return result


def _accepted_task(raw: Any, label: str) -> dict[str, Any]:
    if not isinstance(raw, dict):
        _fail(f"{label} must contain an accepted task identity")
    required = (
        "project_id",
        "run_id",
        "cluster",
        "task_id",
        "step_id",
        "attempt_id",
        "result_ref",
        "status",
        "acceptance_status",
        "accepted_at",
        "acceptance_ref",
    )
    task = {key: _text(raw.get(key), f"{label}.{key}") for key in required if key != "acceptance_ref"}
    if task["status"] != "success":
        _fail(f"{label}.status must be success")
    if task["acceptance_status"] != "accepted":
        _fail(f"{label}.acceptance_status must be accepted")
    task["acceptance_ref"] = _source_ref(
        raw.get("acceptance_ref"),
        f"{label}.acceptance_ref",
        expected_task_id=task["task_id"],
    )
    if "job_id" in raw:
        task["job_id"] = _text(raw["job_id"], f"{label}.job_id")
    return task


def _accepted_index(raw: Any) -> dict[str, Any]:
    if not isinstance(raw, dict):
        _fail("calculation index root must be an object")
    if raw.get("schema_version") != SCHEMA_VERSION:
        _fail(f"calculation index schema_version must be {SCHEMA_VERSION}")
    if raw.get("status") != "accepted":
        _fail("calculation index status must be accepted")
    acceptance = raw.get("acceptance")
    if not isinstance(acceptance, dict) or acceptance.get("status") != "accepted":
        _fail("calculation index requires an accepted top-level acceptance record")
    acceptance = {
        "status": "accepted",
        "accepted_at": _text(acceptance.get("accepted_at"), "acceptance.accepted_at"),
        "source_ref": _source_ref(acceptance.get("source_ref"), "acceptance.source_ref"),
    }
    material_id = _text(raw.get("material_id"), "material_id")
    method = _text(raw.get("calculation_method"), "calculation_method")
    if "HSE06" not in method.upper():
        _fail("calculation index calculation_method must identify HSE06")

    pristine_raw = raw.get("pristine")
    if not isinstance(pristine_raw, dict):
        _fail("pristine must provide the accepted host calculation")
    pristine_task = _accepted_task(pristine_raw.get("task"), "pristine.task")
    host_energy = _number(pristine_raw.get("E_host_ev"), "pristine.E_host_ev")
    host_energy_ref = _source_ref(
        pristine_raw.get("E_host_source_ref"),
        "pristine.E_host_source_ref",
        expected_task_id=pristine_task["task_id"],
    )
    lattice = pristine_raw.get("supercell_lattice_angstrom")
    lattice_ref = _source_ref(
        pristine_raw.get("lattice_source_ref"),
        "pristine.lattice_source_ref",
        expected_task_id=pristine_task["task_id"],
    )
    band_raw = pristine_raw.get("band_edges")
    if not isinstance(band_raw, dict):
        _fail("pristine.band_edges is required")
    band_method = _text(
        band_raw.get("calculation_method"), "pristine.band_edges.calculation_method"
    )
    if band_method != method:
        _fail("pristine band-edge method must match the calculation index method")
    vbm = _number(band_raw.get("vbm_reference_ev"), "pristine.band_edges.vbm_reference_ev")
    cbm = _number(band_raw.get("cbm_reference_ev"), "pristine.band_edges.cbm_reference_ev")
    gap = _number(band_raw.get("gap_ev"), "pristine.band_edges.gap_ev")
    if gap <= 0 or not math.isclose(cbm - vbm, gap, rel_tol=1e-7, abs_tol=1e-5):
        _fail("pristine band edges require a positive gap equal to CBM minus VBM")
    if band_raw.get("calculation_type") != "self_consistent":
        _fail("pristine band edges must be from a self-consistent calculation")
    if band_raw.get("converged") is not True:
        _fail("pristine band-edge calculation must be accepted as converged")
    kpoints = band_raw.get("kpoints")
    if not isinstance(kpoints, dict) or kpoints != {
        "scheme": "Gamma",
        "mesh": [1, 1, 1],
    }:
        _fail("pristine band edges must use Gamma-only 1x1x1")
    band_refs_raw = band_raw.get("source_refs")
    required_band_sources = (
        "vbm_reference_ev",
        "cbm_reference_ev",
        "gap_ev",
        "calculation_type",
        "converged",
        "kpoints",
    )
    if not isinstance(band_refs_raw, dict):
        _fail("pristine.band_edges.source_refs are required")
    band_refs = {
        field: _source_ref(
            band_refs_raw.get(field),
            f"pristine.band_edges.source_refs.{field}",
            expected_task_id=pristine_task["task_id"],
        )
        for field in required_band_sources
    }

    states_raw = raw.get("states")
    if not isinstance(states_raw, list) or not states_raw:
        _fail("states must be a non-empty list")
    states = []
    state_ids: set[str] = set()
    required_state_sources = (
        "state_id",
        "defect_id",
        "geometry_id",
        "charge",
        "n_i",
        "E_def_ev",
        "E_corr_ev",
        "supercell_lattice_angstrom",
    )
    for index, item in enumerate(states_raw):
        label = f"states[{index}]"
        if not isinstance(item, dict):
            _fail(f"{label} must be an object")
        state_id = _text(item.get("state_id"), f"{label}.state_id")
        if state_id in state_ids:
            _fail(f"duplicate state_id: {state_id}")
        state_ids.add(state_id)
        task = _accepted_task(item.get("task"), f"{label}.task")
        if _text(item.get("material_id"), f"{label}.material_id") != material_id:
            _fail(f"{label}.material_id does not match the pristine host")
        if _text(item.get("calculation_method"), f"{label}.calculation_method") != method:
            _fail(f"{label}.calculation_method does not match the index method")

        n_i = item.get("n_i")
        if not isinstance(n_i, dict):
            _fail(f"{label}.n_i is required")
        parsed_n_i = {}
        for element, count in n_i.items():
            element = _text(element, f"{label}.n_i element")
            amount = _number(count, f"{label}.n_i.{element}")
            if not math.isclose(amount, round(amount), rel_tol=0.0, abs_tol=_TOL):
                _fail(f"{label}.n_i.{element} must be an integer atom count")
            parsed_n_i[element] = int(round(amount))

        correction = _number(item.get("E_corr_ev"), f"{label}.E_corr_ev")
        correction_status = item.get("E_corr_status")
        correction_description = _text(
            item.get("E_corr_description"), f"{label}.E_corr_description"
        )
        if correction_status not in {"applied", "justified_zero"}:
            _fail(f"{label}.E_corr_status must be applied or justified_zero")
        if correction_status == "justified_zero" and not math.isclose(
            correction, 0.0, rel_tol=0.0, abs_tol=_TOL
        ):
            _fail(f"{label}.E_corr_status=justified_zero requires E_corr_ev=0")
        source_refs_raw = item.get("source_refs")
        if not isinstance(source_refs_raw, dict):
            _fail(f"{label}.source_refs are required")
        source_refs = {
            field: _source_ref(
                source_refs_raw.get(field),
                f"{label}.source_refs.{field}",
                expected_task_id=task["task_id"],
            )
            for field in required_state_sources
        }
        state = {
            "state_id": state_id,
            "defect_id": _text(item.get("defect_id"), f"{label}.defect_id"),
            "geometry_id": _text(item.get("geometry_id"), f"{label}.geometry_id"),
            "charge": item.get("charge"),
            "n_i": parsed_n_i,
            "E_def_ev": _number(item.get("E_def_ev"), f"{label}.E_def_ev"),
            "E_corr_ev": correction,
            "E_corr_status": correction_status,
            "E_corr_description": correction_description,
            "material_id": material_id,
            "calculation_method": method,
            "supercell_lattice_angstrom": item.get("supercell_lattice_angstrom"),
            "task": task,
            "source_refs": source_refs,
        }
        states.append(state)

    references_raw = raw.get("elemental_references")
    if not isinstance(references_raw, dict):
        _fail("elemental_references must provide HSE06 references for every required element")
    references = {}
    for element, item in references_raw.items():
        element = _text(element, "elemental_references element")
        if not isinstance(item, dict):
            _fail(f"elemental_references.{element} must be an object")
        if item.get("method") != method or item.get("unit") != "eV/atom":
            _fail(f"elemental_references.{element} must be {method} eV/atom")
        task = _accepted_task(item.get("task"), f"elemental_references.{element}.task")
        references[element] = {
            "energy_ev": _number(
                item.get("energy_ev"), f"elemental_references.{element}.energy_ev"
            ),
            "method": method,
            "unit": "eV/atom",
            "reference_id": _text(
                item.get("reference_id"),
                f"elemental_references.{element}.reference_id",
            ),
            "task": task,
            "source_ref": _source_ref(
                item.get("source_ref"),
                f"elemental_references.{element}.source_ref",
                expected_task_id=task["task_id"],
            ),
        }

    return {
        "schema_version": SCHEMA_VERSION,
        "status": "accepted",
        "acceptance": acceptance,
        "material_id": material_id,
        "calculation_method": method,
        "pristine": {
            "task": pristine_task,
            "E_host_ev": host_energy,
            "E_host_source_ref": host_energy_ref,
            "supercell_lattice_angstrom": lattice,
            "lattice_source_ref": lattice_ref,
            "band_edges": {
                "calculation_method": band_method,
                "calculation_type": "self_consistent",
                "kpoints": {"scheme": "Gamma", "mesh": [1, 1, 1]},
                "converged": True,
                "vbm_reference_ev": vbm,
                "cbm_reference_ev": cbm,
                "gap_ev": gap,
                "source_refs": band_refs,
            },
        },
        "states": states,
        "elemental_references": references,
    }


def _validate_phase_point(
    phase_analysis: dict[str, Any], point: dict[str, Any]
) -> tuple[list[str], dict[str, float], dict[str, float]]:
    if phase_analysis.get("status") != "accepted" or phase_analysis.get("schema_version") != 1:
        _fail("phase-analysis input must be an accepted schema-version-1 result")
    host = phase_analysis.get("host")
    if not isinstance(host, dict):
        _fail("phase-analysis result has no host record")
    formula = host.get("formula")
    elements = host.get("elements")
    if (
        not isinstance(formula, dict)
        or not isinstance(elements, list)
        or not elements
        or set(formula) != set(elements)
    ):
        _fail("phase-analysis host formula and element list are incomplete")
    values = point.get("delta_mu_eV")
    if not isinstance(values, dict) or set(values) != set(elements):
        _fail(f"selected point {point.get('name')} must contain the full host delta_mu vector")
    delta_mu = {
        element: _number(values[element], f"selected_points.{point.get('name')}.{element}")
        for element in elements
    }
    if any(value > _TOL for value in delta_mu.values()):
        _fail(f"selected point {point.get('name')} exceeds an elemental upper bound")
    host_hf = _number(
        host.get("formation_energy_eV_per_formula"),
        "phase-analysis host formation energy",
    )
    host_residual = sum(
        _number(formula[element], f"host.formula.{element}") * delta_mu[element]
        for element in elements
    ) - host_hf
    if abs(host_residual) > _TOL * max(1.0, abs(host_hf)):
        _fail(f"selected point {point.get('name')} violates the host equality")

    constraints = phase_analysis.get("constraints")
    if not isinstance(constraints, list) or not constraints:
        _fail("phase-analysis result must include its complete constraint list")
    expected_margins = {}
    for index, constraint in enumerate(constraints):
        label = f"constraints[{index}]"
        if not isinstance(constraint, dict):
            _fail(f"{label} must be an object")
        coefficients = constraint.get("coefficients")
        if not isinstance(coefficients, dict) or set(coefficients) != set(elements):
            _fail(f"{label}.coefficients must cover every host element")
        lhs = sum(
            _number(coefficients[element], f"{label}.coefficients.{element}")
            * delta_mu[element]
            for element in elements
        )
        rhs = _number(constraint.get("rhs_eV"), f"{label}.rhs_eV")
        constraint_id = _text(constraint.get("id"), f"{label}.id")
        if constraint.get("equality"):
            if abs(lhs - rhs) > _TOL * max(1.0, abs(rhs)):
                _fail(f"selected point {point.get('name')} violates {constraint_id}")
        else:
            margin = rhs - lhs
            if margin < -_TOL * max(1.0, abs(rhs)):
                _fail(f"selected point {point.get('name')} violates {constraint_id}")
            expected_margins[constraint_id] = margin

    recorded_margins = point.get("constraint_margins_eV")
    if not isinstance(recorded_margins, dict) or set(recorded_margins) != set(
        expected_margins
    ):
        _fail(f"selected point {point.get('name')} has incomplete constraint-margin provenance")
    for constraint_id, margin in expected_margins.items():
        recorded = _number(
            recorded_margins[constraint_id],
            f"selected_points.{point.get('name')}.constraint_margins_eV.{constraint_id}",
        )
        if not math.isclose(recorded, margin, rel_tol=1e-7, abs_tol=_TOL):
            _fail(f"selected point {point.get('name')} has inconsistent {constraint_id} margin")
    if set(point.get("constraints_source", [])) != {
        constraint["id"] for constraint in constraints
    }:
        _fail(f"selected point {point.get('name')} has incomplete constraint sources")

    external = phase_analysis.get("external_fixed_mu_eV", {})
    if not isinstance(external, dict):
        _fail("phase-analysis external_fixed_mu_eV must be an object")
    external_mu = {
        _text(element, "external element"): _number(value, f"external_fixed_mu_eV.{element}")
        for element, value in external.items()
    }
    binary = phase_analysis.get("dopant_reservoir")
    per_point = point.get("external_fixed_mu_eV", external_mu)
    if not isinstance(per_point, dict):
        _fail("selected point external_fixed_mu_eV must be an object")
    if binary is not None:
        if not isinstance(binary, dict):
            _fail("dopant_reservoir must be an object")
        dopant = _text(binary.get("dopant"), "dopant_reservoir.dopant")
        ligand = _text(binary.get("ligand"), "dopant_reservoir.ligand")
        formula = binary.get("formula")
        if (
            not isinstance(formula, dict) or set(formula) != {dopant, ligand}
            or ligand not in elements or dopant in elements or dopant in external_mu
        ):
            _fail("dopant_reservoir has an inconsistent binary host/dopant formula")
        formation_energy = _number(
            binary.get("formation_energy_eV_per_formula"),
            "dopant_reservoir.formation_energy_eV_per_formula",
        )
        dopant_count = _number(formula[dopant], "dopant_reservoir.dopant_count")
        ligand_count = _number(formula[ligand], "dopant_reservoir.ligand_count")
        if dopant_count <= 0 or ligand_count <= 0:
            _fail("dopant_reservoir formula counts must be positive")
        uncapped = (formation_energy - ligand_count * delta_mu[ligand]) / dopant_count
        expected = min(0.0, uncapped)
        if set(per_point) != set(external_mu) | {dopant} or not math.isclose(
            _number(per_point.get(dopant), f"selected_points.{dopant}"),
            expected, rel_tol=1e-7, abs_tol=_TOL,
        ):
            _fail(f"selected point has an inconsistent {dopant} binary-reservoir limit")
        audit = point.get("dopant_reservoir_audit")
        margin = formation_energy - dopant_count * expected - ligand_count * delta_mu[ligand]
        if (
            not isinstance(audit, dict)
            or audit.get("phase_name") != binary.get("phase_name")
            or audit.get("source") != binary.get("source")
            or audit.get("dopant") != dopant
            or audit.get("ligand") != ligand
            or audit.get("active_bound") != (
                "elemental" if uncapped >= 0 else "binary_phase"
            )
            or not math.isclose(
                _number(
                    audit.get("uncapped_binary_limit_eV"),
                    "dopant_reservoir_audit.uncapped_binary_limit_eV",
                ),
                uncapped, rel_tol=1e-7, abs_tol=_TOL,
            )
            or not math.isclose(
                _number(audit.get("selected_delta_mu_eV"), "dopant_reservoir_audit.selected_delta_mu_eV"),
                expected, rel_tol=1e-7, abs_tol=_TOL,
            )
            or not math.isclose(
                _number(
                    audit.get("binary_constraint_margin_eV"),
                    "dopant_reservoir_audit.binary_constraint_margin_eV",
                ),
                margin, rel_tol=1e-7, abs_tol=_TOL,
            )
        ):
            _fail("selected point lacks consistent binary-reservoir audit")
        external_mu[dopant] = expected
    elif set(per_point) != set(external_mu):
        _fail("selected point adds an undeclared external chemical potential")
    for element, expected in external_mu.items():
        if not math.isclose(
            _number(per_point.get(element), f"selected_points.external_fixed_mu_eV.{element}"),
            expected, rel_tol=1e-7, abs_tol=_TOL,
        ):
            _fail(f"selected point external {element} chemical potential is inconsistent")
    if any(value > _TOL for value in external_mu.values()):
        _fail("external fixed chemical potentials must not exceed the elemental upper bound")
    return elements, delta_mu, external_mu


def _select_phase_points(
    phase_analysis: dict[str, Any], point_names: Sequence[str] | None
) -> list[dict[str, Any]]:
    points = phase_analysis.get("selected_points")
    if not isinstance(points, list):
        _fail("phase-analysis selected_points must be a list")
    if not points:
        selection_audit = phase_analysis.get("point_selection", {})
        if (
            isinstance(selection_audit, dict)
            and selection_audit.get("selection_status")
            == "requires_explicit_control_element"
        ):
            _fail(
                "phase domain is available but has no selected points; set an explicit "
                "point_selection.control_element and rerun phase analysis, or add "
                "complete manual selected_points and select them with --point-name"
            )
        if (
            isinstance(selection_audit, dict)
            and selection_audit.get("selection_status") == "manual_points_only"
        ):
            _fail(
                "manual selection mode has no selected points; add complete manual "
                "selected_points and select them with --point-name"
            )
        _fail("phase-analysis result contains no selected points")
    if point_names:
        requested = list(point_names)
        if len(requested) != len(set(requested)):
            _fail("point names must be unique")
        by_name = {point.get("name"): point for point in points if isinstance(point, dict)}
        missing = set(requested) - set(by_name)
        if missing:
            _fail(f"unknown selected point names: {', '.join(sorted(missing))}")
        selected = [by_name[name] for name in requested]
    else:
        selected = [
            point
            for point in points
            if isinstance(point, dict)
            and point.get("selection_source")
            == "automatic_extreme_face_convex_interpolation"
        ]
        if not selected:
            selection_audit = phase_analysis.get("point_selection", {})
            if (
                isinstance(selection_audit, dict)
                and selection_audit.get("selection_status") == "manual_points_only"
            ):
                _fail(
                    "phase analysis has only manual points; select an audited point "
                    "with --point-name"
                )
            reason = (
                selection_audit.get("control_element_reason")
                if isinstance(selection_audit, dict)
                else None
            )
            if reason:
                _fail(
                    "phase domain is available but automatic formation points were not "
                    f"selected ({reason}); rerun phase analysis with an explicit "
                    "point_selection.control_element or specify audited --point-name values"
                )
            _fail(
                "no automatic phase points; specify point_selection.control_element "
                "and rerun phase analysis, or specify audited --point-name values"
            )
    for point in selected:
        if not isinstance(point, dict):
            _fail("selected phase points must be objects")
    return selected


def _assemble_input(
    index: dict[str, Any],
    phase_analysis: dict[str, Any],
    point: dict[str, Any],
    elements: list[str],
    delta_mu: dict[str, float],
    external_mu: dict[str, float],
    *,
    calculation_index_path: str,
    phase_analysis_path: str,
) -> dict[str, Any]:
    pristine = index["pristine"]
    required_elements = set(elements) | set(external_mu)
    for state in index["states"]:
        required_elements.update(state["n_i"])
    missing_refs = required_elements - set(index["elemental_references"])
    if missing_refs:
        _fail(
            "missing accepted HSE06 elemental references for: "
            + ", ".join(sorted(missing_refs))
        )

    cp_values = {}
    for element in sorted(required_elements):
        if element in delta_mu:
            relative = delta_mu[element]
        elif element in external_mu:
            relative = external_mu[element]
        else:
            _fail(f"no selected chemical potential is available for {element}")
        cp_values[element] = {
            "delta_mu_ev": relative,
            "elemental_reference_ev": index["elemental_references"][element][
                "energy_ev"
            ],
        }

    states = []
    state_provenance = {}
    for item in index["states"]:
        state = {
            key: item[key]
            for key in (
                "state_id",
                "defect_id",
                "geometry_id",
                "charge",
                "n_i",
                "E_def_ev",
                "E_corr_ev",
                "material_id",
                "calculation_method",
                "supercell_lattice_angstrom",
            )
        }
        state["E_host_ev"] = pristine["E_host_ev"]
        state["source_refs"] = deepcopy(item["source_refs"])
        state["task"] = deepcopy(item["task"])
        state["E_corr_status"] = item["E_corr_status"]
        state["E_corr_description"] = item["E_corr_description"]
        states.append(state)
        state_provenance[item["state_id"]] = {
            "task": deepcopy(item["task"]),
            "source_refs": deepcopy(item["source_refs"]),
            "E_host_source_ref": deepcopy(pristine["E_host_source_ref"]),
            "correction_status": item["E_corr_status"],
            "correction_description": item["E_corr_description"],
        }

    band = pristine["band_edges"]
    band_edges = {
        "material_id": index["material_id"],
        "supercell_lattice_angstrom": pristine["supercell_lattice_angstrom"],
        **{
            key: band[key]
            for key in (
                "calculation_method",
                "calculation_type",
                "kpoints",
                "converged",
                "vbm_reference_ev",
                "cbm_reference_ev",
                "gap_ev",
            )
        },
        "provenance": json.dumps(
            {
                "task": pristine["task"],
                "source_refs": band["source_refs"],
            },
            ensure_ascii=False,
            sort_keys=True,
        ),
        "source_refs": deepcopy(band["source_refs"]),
        "task": deepcopy(pristine["task"]),
    }
    point_record = deepcopy(point)
    point_record["delta_mu_eV"] = delta_mu
    origin = {
        "calculation_index": {
            "path": calculation_index_path,
            "schema_version": SCHEMA_VERSION,
            "status": "accepted",
            "acceptance": deepcopy(index["acceptance"]),
            "pristine_task": deepcopy(pristine["task"]),
            "E_host_source_ref": deepcopy(pristine["E_host_source_ref"]),
            "lattice_source_ref": deepcopy(pristine["lattice_source_ref"]),
            "state_provenance": state_provenance,
            "elemental_reference_provenance": {
                element: {
                    "reference_id": index["elemental_references"][element][
                        "reference_id"
                    ],
                    "task": deepcopy(
                        index["elemental_references"][element]["task"]
                    ),
                    "source_ref": deepcopy(
                        index["elemental_references"][element]["source_ref"]
                    ),
                }
                for element in sorted(required_elements)
            },
        },
        "phase_analysis": {
            "path": phase_analysis_path,
            "status": phase_analysis["status"],
            "point": point_record,
            "constraint_sources": point["constraints_source"],
            "constraint_margins_eV": point["constraint_margins_eV"],
            "external_fixed_mu_eV": external_mu,
        },
    }
    return {
        "schema_version": SCHEMA_VERSION,
        "origin": origin,
        "material_id": index["material_id"],
        "calculation_method": index["calculation_method"],
        "pristine_band_edges": band_edges,
        "chemical_potentials": {
            "mode": "relative_plus_elemental_reference",
            "relative_method": phase_analysis["chemical_potential_functional"],
            "elemental_reference_method": index["calculation_method"],
            "values": cp_values,
        },
        "states": states,
    }


def _safe_name(value: str) -> str:
    slug = re.sub(r"[^A-Za-z0-9_-]+", "_", value).strip("_")
    return slug[:64] or "point"


def preflight_bundle(
    calculation_index: dict[str, Any],
    phase_analysis: dict[str, Any],
    *,
    calculation_index_path: str | Path,
    phase_analysis_path: str | Path,
    point_names: Sequence[str] | None = None,
    separate_defect_ids: Sequence[str] = (),
    max_defects_per_figure: int = 6,
    y_max_ev: float | None = None,
) -> tuple[list[tuple[str, dict[str, Any], dict[str, Any]]], list[str]]:
    """Validate every selected point and plot limit without writing artifacts."""
    if (
        isinstance(max_defects_per_figure, bool)
        or not isinstance(max_defects_per_figure, int)
        or max_defects_per_figure < 1
    ):
        _fail("max_defects_per_figure must be a positive integer")
    index = _accepted_index(calculation_index)
    if not isinstance(phase_analysis, dict):
        _fail("phase-analysis input must be an object")
    if phase_analysis.get("chemical_potential_functional") != "PBE":
        _fail("phase-analysis chemical potentials must be relative to PBE references")
    points = _select_phase_points(phase_analysis, point_names)
    prepared = []
    point_names_seen = set()
    for point in points:
        name = _text(point.get("name"), "selected point name")
        if name in point_names_seen:
            _fail(f"duplicate selected point name: {name}")
        point_names_seen.add(name)
        elements, delta_mu, external_mu = _validate_phase_point(phase_analysis, point)
        data = _assemble_input(
            index,
            phase_analysis,
            point,
            elements,
            delta_mu,
            external_mu,
            calculation_index_path=str(Path(calculation_index_path).expanduser().resolve()),
            phase_analysis_path=str(Path(phase_analysis_path).expanduser().resolve()),
        )
        audit = formation.analyze(data)
        formation._y_bounds(audit["defects"], y_max_ev)
        prepared.append((name, data, audit))

    separate_ids = list(separate_defect_ids)
    if len(separate_ids) != len(set(separate_ids)):
        _fail("separate_defect_ids must not contain duplicates")
    known_defects = {item["defect_id"] for item in prepared[0][2]["defects"]}
    unknown_defects = set(separate_ids) - known_defects
    if unknown_defects:
        _fail(f"unknown separate defect IDs: {', '.join(sorted(unknown_defects))}")
    return prepared, separate_ids


def build_bundle(
    calculation_index: dict[str, Any],
    phase_analysis: dict[str, Any],
    output_dir: str | Path,
    *,
    calculation_index_path: str | Path,
    phase_analysis_path: str | Path,
    point_names: Sequence[str] | None = None,
    separate_defect_ids: Sequence[str] = (),
    max_defects_per_figure: int = 6,
    y_max_ev: float | None = None,
    title: str | None = None,
) -> dict[str, Any]:
    prepared, separate_ids = preflight_bundle(
        calculation_index,
        phase_analysis,
        calculation_index_path=calculation_index_path,
        phase_analysis_path=phase_analysis_path,
        point_names=point_names,
        separate_defect_ids=separate_defect_ids,
        max_defects_per_figure=max_defects_per_figure,
        y_max_ev=y_max_ev,
    )

    root = Path(output_dir).expanduser().resolve()
    if root.exists() and (not root.is_dir() or any(root.iterdir())):
        _fail(f"output directory must be new or empty; refusing to overwrite: {root}")
    if not root.exists():
        root.mkdir(parents=True)

    point_outputs = []
    for index_number, (name, data, audit) in enumerate(prepared, start=1):
        point_dir = root / f"point_{index_number:02d}_{_safe_name(name)}"
        point_dir.mkdir()
        input_path = point_dir / "formation-input.json"
        input_path.write_text(
            json.dumps(data, indent=2, ensure_ascii=False, allow_nan=False) + "\n",
            encoding="utf-8",
        )
        result = formation.write_outputs(
            audit,
            point_dir / "formation",
            max_defects_per_figure=max_defects_per_figure,
            separate_defect_ids=separate_ids,
            y_max_ev=y_max_ev,
            title=title,
            input_path=input_path,
        )
        point_outputs.append(
            {
                "point_name": name,
                "formation_input": str(input_path),
                "audit_json": result["output"]["audit_json"],
                "figures": result["output"]["figures"],
                "selection_source": data["origin"]["phase_analysis"]["point"][
                    "selection_source"
                ],
                "control_element": data["origin"]["phase_analysis"]["point"].get(
                    "control_element"
                ),
                "fraction_from_poor_to_rich": data["origin"]["phase_analysis"][
                    "point"
                ].get("fraction_from_poor_to_rich"),
            }
        )

    manifest = {
        "schema_version": SCHEMA_VERSION,
        "status": "accepted",
        "calculation_index": str(
            Path(calculation_index_path).expanduser().resolve()
        ),
        "phase_analysis": str(Path(phase_analysis_path).expanduser().resolve()),
        "chemical_potential_convention": "PBE relative delta_mu plus HSE06 elemental reference",
        "points": point_outputs,
        "separate_defect_ids": separate_ids,
    }
    manifest_path = root / "formation-energy-bundle.json"
    manifest_path.write_text(
        json.dumps(manifest, indent=2, ensure_ascii=False, allow_nan=False) + "\n",
        encoding="utf-8",
    )
    manifest["manifest"] = str(manifest_path)
    return manifest


def main(argv: Sequence[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        description=(
            "Build and plot formation-energy inputs from an accepted structured "
            "calculation index and a phase-analysis result."
        )
    )
    parser.add_argument("--calculation-index", required=True)
    parser.add_argument("--phase-analysis", required=True)
    parser.add_argument("--output-dir", required=True)
    parser.add_argument(
        "--point-name",
        action="append",
        help="explicitly select an audited point name; repeatable; default is automatic points",
    )
    parser.add_argument(
        "--separate-defect-id",
        action="append",
        default=[],
        help="plot this explicitly named defect separately; repeatable",
    )
    parser.add_argument("--max-defects-per-figure", type=int, default=6)
    parser.add_argument("--y-max-ev", type=float)
    parser.add_argument("--title")
    args = parser.parse_args(argv)
    try:
        index_path = Path(args.calculation_index).expanduser().resolve()
        phase_path = Path(args.phase_analysis).expanduser().resolve()
        calculation_index = formation.load_json(index_path)
        phase_analysis = formation.load_json(phase_path)
        manifest = build_bundle(
            calculation_index,
            phase_analysis,
            args.output_dir,
            calculation_index_path=index_path,
            phase_analysis_path=phase_path,
            point_names=args.point_name,
            separate_defect_ids=args.separate_defect_id,
            max_defects_per_figure=args.max_defects_per_figure,
            y_max_ev=args.y_max_ev,
            title=args.title,
        )
    except (OSError, ValueError, TypeError) as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 2
    for point in manifest["points"]:
        print(f"Built {point['point_name']}: {point['formation_input']}")
        print(f"Audit: {point['audit_json']}")
        for figure in point["figures"]:
            print(f"Figure: {figure['files']['svg']} / {figure['files']['png']}")
    print(f"Bundle manifest: {manifest['manifest']}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
