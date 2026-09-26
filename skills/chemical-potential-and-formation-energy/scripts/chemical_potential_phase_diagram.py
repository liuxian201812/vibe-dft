#!/usr/bin/env python3
"""Build chemical-potential stability domains for two- through five-element hosts.

Input formulas are explicit element-to-count mappings. Energies must be either
consistent formation energies per formula unit or total energies referenced to
one consistent set of elemental energies. All plotted and reported chemical
potentials are relative values, delta_mu_i = mu_i - mu_i(element), in eV.

Example structure (formation-energy mode):
{
  "energy_basis": "formation_energy",
  "correction_policy": {
    "scheme": "none",
    "mode": "already_applied"
  },
  "host": {
    "formula": {"Cs": 2, "Zn": 1, "Cl": 4},
    "energy": {
      "value": -8.0, "unit": "eV/formula", "functional": "PBE",
      "elemental_reference": "PBE-2026-set-A",
      "correction": {
        "value": 0.0, "unit": "eV/formula",
        "source": "none", "applied": true
      }
    }
  },
  "phases": [
    {
      "name": "ZnCl2",
      "formula": {"Zn": 1, "Cl": 2},
      "energy": {
        "value": -2.0, "unit": "eV/formula", "functional": "PBE",
        "elemental_reference": "PBE-2026-set-A",
        "correction": {
          "value": 0.0, "unit": "eV/formula",
          "source": "none", "applied": true
        }
      }
    }
  ],
  "plot": {"axes": ["Cl", "Zn"]},
  "slices": [],
  "selected_points": []
}

For total-energy mode set energy_basis to "total_energy", use unit "eV" for
host/phases, and supply elemental_energies as per-atom values in eV/atom with
matching functional and elemental_reference fields. Every energy row also
declares a correction object with value, unit, source, and applied. The common
correction_policy.mode is "already_applied" (audit only; never add again) or
"apply" (input values are raw and the declared correction is added once).
External elements in a competitor phase must have an explicit external_fixed_mu
value in eV relative to the same elemental reference; they are never silently
dropped. An explicitly fixed dopant need not occur in a host competitor phase:
its typed delta_mu <= 0 carries through every selected point without changing
the host stability domain.

The host equality is eliminated analytically. The remaining one- through
four-dimensional bounded half-space is solved by vertex enumeration, including
lower-dimensional degenerate domains. For five-element hosts, the plot is an
explicitly labeled two-dimensional projection or a fixed two-dimensional slice;
the analysis retains the complete four-dimensional vertices and five-element
chemical-potential vectors.
"""

from __future__ import annotations

import argparse
import itertools
import json
import math
import sys
from pathlib import Path
from typing import Any

import numpy as np


TOL = 1.0e-8
DEFAULT_POINT_FRACTIONS = (0.25, 0.5, 0.75)


class DomainError(ValueError):
    """Invalid input or a chemical-potential domain with no feasible points."""


def _as_finite_number(value: Any, label: str) -> float:
    if isinstance(value, bool) or not isinstance(value, (int, float)):
        raise DomainError(f"{label} must be a finite number")
    number = float(value)
    if not math.isfinite(number):
        raise DomainError(f"{label} must be a finite number")
    return number


def _formula(raw: Any, label: str) -> dict[str, float]:
    if not isinstance(raw, dict) or not raw:
        raise DomainError(f"{label} must be a non-empty element-to-count object")
    result: dict[str, float] = {}
    for element, count in raw.items():
        if not isinstance(element, str) or not element.strip():
            raise DomainError(f"{label} contains an invalid element name")
        amount = _as_finite_number(count, f"{label}.{element}")
        if amount <= 0:
            raise DomainError(f"{label}.{element} must be positive")
        result[element] = amount
    return result


def _correction_policy(data: dict[str, Any]) -> dict[str, str]:
    policy = data.get("correction_policy")
    if not isinstance(policy, dict):
        raise DomainError("correction_policy must declare a common scheme and mode")
    scheme = policy.get("scheme")
    mode = policy.get("mode")
    if not isinstance(scheme, str) or not scheme.strip():
        raise DomainError("correction_policy.scheme must be a non-empty string")
    if mode not in {"already_applied", "apply"}:
        raise DomainError("correction_policy.mode must be 'already_applied' or 'apply'")
    return {"scheme": scheme, "mode": mode}


def _energy_metadata(
    record: Any,
    expected_unit: str,
    label: str,
    correction_policy: dict[str, str],
) -> tuple[float, dict[str, str], dict[str, Any]]:
    if not isinstance(record, dict):
        raise DomainError(f"{label} must be an energy object")
    raw_value = _as_finite_number(record.get("value"), f"{label}.value")
    unit = record.get("unit")
    functional = record.get("functional")
    reference = record.get("elemental_reference")
    correction = record.get("correction")
    if unit != expected_unit:
        raise DomainError(f"{label}.unit must be {expected_unit!r}, got {unit!r}")
    if not all(isinstance(item, str) and item.strip() for item in (functional, reference)):
        raise DomainError(
            f"{label} must explicitly declare functional and elemental_reference"
        )
    if not isinstance(correction, dict):
        raise DomainError(
            f"{label}.correction must declare value, unit, source, and applied"
        )
    correction_value = _as_finite_number(
        correction.get("value"), f"{label}.correction.value"
    )
    correction_unit = correction.get("unit")
    correction_source = correction.get("source")
    applied = correction.get("applied")
    if correction_unit != expected_unit:
        raise DomainError(
            f"{label}.correction.unit must be {expected_unit!r}, got {correction_unit!r}"
        )
    if not isinstance(correction_source, str) or not correction_source.strip():
        raise DomainError(f"{label}.correction.source must be a non-empty string")
    if not isinstance(applied, bool):
        raise DomainError(f"{label}.correction.applied must be true or false")
    if correction_policy["mode"] == "already_applied" and not applied:
        raise DomainError(
            f"{label} declares correction not applied, but correction_policy.mode "
            "is 'already_applied'"
        )
    if correction_policy["mode"] == "apply" and applied:
        raise DomainError(
            f"{label} correction is already applied; refusing to apply it twice"
        )
    if correction_source.lower() == "none" and abs(correction_value) > TOL:
        raise DomainError(
            f"{label} has a non-zero correction but source is 'none'"
        )
    if correction_policy["scheme"].lower() == "none" and abs(correction_value) > TOL:
        raise DomainError(
            f"{label} has a non-zero correction while correction_policy.scheme is 'none'"
        )
    effective_value = (
        raw_value
        if correction_policy["mode"] == "already_applied"
        else raw_value + correction_value
    )
    audit = {
        "input_value": raw_value,
        "correction_value": correction_value,
        "correction_source": correction_source,
        "applied_in_input": applied,
        "application_action": (
            "not_reapplied" if applied else "applied_once"
        ),
        "effective_value": effective_value,
        "unit": expected_unit,
    }
    return effective_value, {
        "unit": unit,
        "functional": functional,
        "elemental_reference": reference,
    }, audit


def _read_energy_rows(data: dict[str, Any], host_formula: dict[str, float]) -> tuple[
    float, list[dict[str, Any]], dict[str, Any]
]:
    basis = data.get("energy_basis")
    if basis not in {"formation_energy", "total_energy"}:
        raise DomainError("energy_basis must be 'formation_energy' or 'total_energy'")
    phases_raw = data.get("phases")
    if not isinstance(phases_raw, list):
        raise DomainError("phases must be a list")
    host_raw = data.get("host")
    if not isinstance(host_raw, dict):
        raise DomainError("host must be an object")
    expected_unit = "eV/formula" if basis == "formation_energy" else "eV"
    energy_key = "formation_energy" if basis == "formation_energy" else "total_energy"
    correction_policy = _correction_policy(data)

    entities: list[dict[str, Any]] = []
    host_energy, host_metadata, host_correction = _energy_metadata(
        host_raw.get("energy"), expected_unit, f"host.{energy_key}", correction_policy
    )
    entities.append(
        {
            "name": str(host_raw.get("name", "host")),
            "formula": host_formula,
            "raw_energy": host_energy,
            "metadata": host_metadata,
            "correction_audit": host_correction,
            "is_host": True,
        }
    )
    for index, phase_raw in enumerate(phases_raw):
        if not isinstance(phase_raw, dict):
            raise DomainError(f"phases[{index}] must be an object")
        name = phase_raw.get("name")
        if not isinstance(name, str) or not name.strip():
            raise DomainError(f"phases[{index}].name must be a non-empty string")
        formula = _formula(phase_raw.get("formula"), f"phases[{index}].formula")
        value, metadata, correction_audit = _energy_metadata(
            phase_raw.get("energy"),
            expected_unit,
            f"phases[{index}].{energy_key}",
            correction_policy,
        )
        entities.append(
            {
                "name": name,
                "formula": formula,
                "raw_energy": value,
                "metadata": metadata,
                "correction_audit": correction_audit,
                "is_host": False,
            }
        )

    common_metadata = entities[0]["metadata"]
    if common_metadata["functional"] != "PBE":
        raise DomainError("chemical-potential stability energies must use functional='PBE'")
    for entity in entities[1:]:
        for key in ("unit", "functional", "elemental_reference"):
            if entity["metadata"][key] != common_metadata[key]:
                raise DomainError(
                    f"energy provenance mismatch for {entity['name']}: {key} differs "
                    "from the host"
                )

    reference_values: dict[str, float] = {}
    reference_corrections: dict[str, Any] = {}
    all_elements = set(host_formula)
    for entity in entities:
        all_elements.update(entity["formula"])
    if basis == "formation_energy":
        for entity in entities:
            entity["formation_energy"] = entity["raw_energy"]
    else:
        raw_refs = data.get("elemental_energies")
        if not isinstance(raw_refs, dict):
            raise DomainError("total_energy mode requires elemental_energies")
        for element in sorted(all_elements):
            record = raw_refs.get(element)
            value, metadata, correction_audit = _energy_metadata(
                record,
                "eV/atom",
                f"elemental_energies.{element}",
                correction_policy,
            )
            if metadata["functional"] != common_metadata["functional"]:
                raise DomainError(f"elemental reference functional mismatch for {element}")
            if metadata["elemental_reference"] != common_metadata["elemental_reference"]:
                raise DomainError(f"elemental reference set mismatch for {element}")
            reference_values[element] = value
            reference_corrections[element] = correction_audit
        for entity in entities:
            entity["formation_energy"] = entity["raw_energy"] - sum(
                count * reference_values[element]
                for element, count in entity["formula"].items()
            )

    provenance = {
        "basis": basis,
        "energy_key": energy_key,
        "common_metadata": common_metadata,
        "correction_policy": correction_policy,
        "corrections": {
            entity["name"]: entity["correction_audit"] for entity in entities
        },
        "elemental_reference_corrections": reference_corrections,
        "elemental_energies_eV_per_atom": reference_values if basis == "total_energy" else None,
        "formation_energies_eV_per_formula": {
            entity["name"]: entity["formation_energy"] for entity in entities
        },
    }
    return entities[0]["formation_energy"], entities[1:], provenance


def _external_mu(
    data: dict[str, Any],
    phases: list[dict[str, Any]],
    host_elements: list[str],
    expected_reference: str,
) -> dict[str, float]:
    needed = {
        element
        for phase in phases
        for element in phase["formula"]
        if element not in host_elements
    }
    raw = data.get("external_fixed_mu", {})
    if not isinstance(raw, dict):
        raise DomainError("external_fixed_mu must be an object")
    missing = sorted(needed - set(raw))
    if missing:
        raise DomainError(
            "competitor phases contain elements outside the host; explicitly provide "
            f"external_fixed_mu for: {', '.join(missing)}"
        )
    overlap = sorted(set(raw) & set(host_elements))
    if overlap:
        raise DomainError(
            "external_fixed_mu must refer to elements outside the host: "
            + ", ".join(overlap)
        )
    result: dict[str, float] = {}
    for element in raw:
        if not isinstance(element, str) or not element.strip():
            raise DomainError("external_fixed_mu contains an invalid element name")
        item = raw[element]
        if (
            not isinstance(item, dict)
            or item.get("unit") != "eV"
            or item.get("elemental_reference") != expected_reference
            or item.get("meaning") != "delta_mu"
        ):
            raise DomainError(
                f"external_fixed_mu.{element} must declare unit='eV', "
                "meaning='delta_mu', and the matching elemental_reference"
            )
        value = _as_finite_number(item.get("value"), f"external_fixed_mu.{element}")
        if value > 0:
            raise DomainError(
                f"external_fixed_mu.{element} exceeds the elemental upper bound (0 eV)"
            )
        result[element] = value
    return result


def _make_full_constraints(
    host_formula: dict[str, float],
    host_hf: float,
    phases: list[dict[str, Any]],
    external_mu: dict[str, float],
) -> list[dict[str, Any]]:
    elements = list(host_formula)
    constraints: list[dict[str, Any]] = []
    for element in elements:
        row = {item: 0.0 for item in elements}
        row[element] = 1.0
        constraints.append(
            {
                "id": f"elemental:{element}",
                "label": f"{element} elemental limit",
                "source": {"type": "elemental_limit", "element": element, "relation": "delta_mu <= 0"},
                "coefficients": row,
                "rhs_eV": 0.0,
            }
        )
    host_row = {element: host_formula[element] for element in elements}
    constraints.append(
        {
            "id": "host:formation_energy",
            "label": "host formation-energy equality",
            "source": {
                "type": "host_stability",
                "formula": host_formula,
                "formation_energy_eV_per_formula": host_hf,
                "relation": "sum(n_i * delta_mu_i) = delta_Hf(host)",
            },
            "coefficients": host_row,
            "rhs_eV": host_hf,
            "equality": True,
        }
    )
    for index, phase in enumerate(phases):
        formula = phase["formula"]
        external_term = sum(
            count * external_mu[element]
            for element, count in formula.items()
            if element not in host_formula
        )
        row = {element: formula.get(element, 0.0) for element in elements}
        constraints.append(
            {
                "id": f"phase:{index}:{phase['name']}",
                "label": phase["name"],
                "source": {
                    "type": "competing_phase",
                    "name": phase["name"],
                    "formula": formula,
                    "formation_energy_eV_per_formula": phase["formation_energy"],
                    "external_fixed_mu_eV": {
                        element: external_mu[element]
                        for element in formula
                        if element in external_mu
                    },
                    "relation": "sum(m_i * delta_mu_i) <= delta_Hf(phase)",
                },
                "coefficients": row,
                "rhs_eV": phase["formation_energy"] - external_term,
            }
        )
    return constraints


def _parameterized_inequalities(
    elements: list[str],
    host_formula: dict[str, float],
    host_hf: float,
    constraints: list[dict[str, Any]],
    axes: list[str],
    fixed: dict[str, float] | None = None,
) -> tuple[np.ndarray, np.ndarray, list[dict[str, Any]], dict[str, Any]]:
    fixed = fixed or {}
    if len(set(axes)) != len(axes) or any(element not in elements for element in axes):
        raise DomainError("coordinate axes must be distinct host elements")
    if any(element not in elements for element in fixed):
        raise DomainError("fixed slice elements must belong to the host")
    if set(axes) & set(fixed):
        raise DomainError("a fixed element cannot also be a coordinate axis")
    remaining = [element for element in elements if element not in axes and element not in fixed]
    if len(remaining) != 1:
        raise DomainError("axes and fixed elements must leave exactly one host element")
    pivot = remaining[0]
    if host_formula[pivot] <= 0:
        raise DomainError("the eliminated host element must have positive stoichiometry")
    for element, value in fixed.items():
        if element not in elements:
            raise DomainError(f"fixed slice element {element!r} is not in the host")
        _as_finite_number(value, f"fixed[{element}]")

    dimension = len(axes)
    offset = np.zeros(len(elements), dtype=float)
    transform = np.zeros((len(elements), dimension), dtype=float)
    index = {element: position for position, element in enumerate(elements)}
    for element, value in fixed.items():
        offset[index[element]] = value
    for axis_index, element in enumerate(axes):
        transform[index[element], axis_index] = 1.0
    fixed_host_sum = sum(host_formula[element] * value for element, value in fixed.items())
    for axis_index, element in enumerate(axes):
        transform[index[pivot], axis_index] = -host_formula[element] / host_formula[pivot]
    offset[index[pivot]] = (host_hf - fixed_host_sum) / host_formula[pivot]

    matrix: list[np.ndarray] = []
    bounds: list[float] = []
    reduced_constraints: list[dict[str, Any]] = []
    for constraint in constraints:
        if constraint.get("equality"):
            continue
        full_coeff = np.array(
            [constraint["coefficients"][element] for element in elements], dtype=float
        )
        coefficients = full_coeff @ transform
        rhs = constraint["rhs_eV"] - float(full_coeff @ offset)
        norm = float(np.linalg.norm(coefficients))
        if norm > TOL:
            coefficients = coefficients / norm
            rhs /= norm
        reduced = {
            "id": constraint["id"],
            "label": constraint["label"],
            "source": constraint["source"],
            "coefficients": coefficients.tolist(),
            "rhs_eV": float(rhs),
        }
        reduced_constraints.append(reduced)
        matrix.append(coefficients)
        bounds.append(float(rhs))
    parameterization = {
        "axes": axes,
        "fixed": fixed,
        "eliminated_element": pivot,
        "host_equality": {
            "formula": host_formula,
            "formation_energy_eV_per_formula": host_hf,
        },
        "dimension_before_constraints": dimension,
    }
    return np.asarray(matrix, dtype=float), np.asarray(bounds, dtype=float), reduced_constraints, parameterization


def _enumerate_vertices(matrix: np.ndarray, bounds: np.ndarray, tolerance: float = TOL) -> list[np.ndarray]:
    dimension = matrix.shape[1]
    if dimension == 0:
        if np.all(-bounds <= tolerance):
            return [np.empty(0, dtype=float)]
        raise DomainError("chemical-potential domain is empty")
    vertices: list[np.ndarray] = []
    for selected in itertools.combinations(range(len(bounds)), dimension):
        submatrix = matrix[list(selected)]
        if np.linalg.matrix_rank(submatrix, tol=tolerance) < dimension:
            continue
        try:
            point = np.linalg.solve(submatrix, bounds[list(selected)])
        except np.linalg.LinAlgError:
            continue
        residual = matrix @ point - bounds
        scale = np.maximum(1.0, np.abs(bounds))
        if np.all(residual <= tolerance * scale):
            if not any(np.linalg.norm(point - existing, ord=np.inf) <= tolerance * 10 for existing in vertices):
                vertices.append(point)
    if not vertices:
        raise DomainError("chemical-potential stability domain is empty")
    vertices.sort(key=lambda item: tuple(float(value) for value in item))
    return vertices


def _four_dimensional_halfspace_prefilter(
    matrix: np.ndarray,
    bounds: np.ndarray,
    constraints: list[dict[str, Any]],
) -> tuple[np.ndarray, np.ndarray, dict[str, Any]]:
    """Discard only inequalities already implied by the retained half-spaces."""
    elemental = [
        index for index, item in enumerate(constraints)
        if item["source"].get("type") == "elemental_limit"
    ]
    if len(elemental) != 5:
        raise DomainError("five-element host requires all five elemental upper bounds")
    try:
        from scipy.optimize import linprog
    except ImportError:
        return matrix, bounds, {
            "algorithm": "all_constraints_vertex_enumeration",
            "retained_count": len(bounds),
            "redundant_count": 0,
            "redundant_constraint_ids": [],
        }
    retained = list(elemental)
    redundant = []
    for index, constraint in enumerate(constraints):
        if index in elemental:
            continue
        result = linprog(
            -matrix[index],
            A_ub=matrix[retained],
            b_ub=bounds[retained],
            bounds=[(None, None)] * 4,
            method="highs",
        )
        # A failed, unbounded, or infeasible optimization cannot prove
        # redundancy: keep the inequality and let vertex enumeration decide.
        if not result.success or result.fun is None or not math.isfinite(float(result.fun)):
            retained.append(index)
            continue
        maximum = -float(result.fun)
        rhs = float(bounds[index])
        if maximum <= rhs + TOL * max(1.0, abs(maximum), abs(rhs)):
            redundant.append(constraint["id"])
        else:
            retained.append(index)
    return matrix[retained], bounds[retained], {
        "algorithm": "sequential_linear_programming_redundancy_check",
        "retained_count": len(retained),
        "redundant_count": len(redundant),
        "redundant_constraint_ids": redundant,
    }


def _affine_dimension(vertices: list[np.ndarray], tolerance: float = TOL) -> int:
    if len(vertices) <= 1:
        return 0
    centered = np.asarray(vertices[1:]) - np.asarray(vertices[0])
    return int(np.linalg.matrix_rank(centered, tol=tolerance * 10))


def _two_dimensional_projection(
    vertices: list[dict[str, Any]], axes: list[str]
) -> dict[str, Any]:
    """Project the complete vertex set; the filled hull is not a 4D slice."""
    points = sorted({
        tuple(float(vertex["delta_mu_eV"][axis]) for axis in axes)
        for vertex in vertices
    })

    def cross(first: tuple[float, float], second: tuple[float, float],
              third: tuple[float, float]) -> float:
        return (
            (second[0] - first[0]) * (third[1] - first[1])
            - (second[1] - first[1]) * (third[0] - first[0])
        )

    if len(points) <= 2:
        hull = points
    else:
        lower: list[tuple[float, float]] = []
        upper: list[tuple[float, float]] = []
        for sequence, destination in ((points, lower), (list(reversed(points)), upper)):
            for point in sequence:
                while len(destination) >= 2 and cross(
                    destination[-2], destination[-1], point
                ) <= 0:
                    destination.pop()
                destination.append(point)
        hull = lower[:-1] + upper[:-1]
    return {
        "type": "2d_linear_projection_of_4d_ambient_domain",
        "axes": axes,
        "is_full_stability_domain": False,
        "hull_coordinates_eV": [list(point) for point in hull],
        "projected_vertex_count": len(points),
        "interpretation": (
            "The displayed hull contains feasible projections, not complete chemical-potential "
            "vectors; use domain.vertices or selected_points for all five elements."
        ),
    }


def _full_mu(
    point: np.ndarray,
    parameterization: dict[str, Any],
    elements: list[str],
    host_formula: dict[str, float],
    host_hf: float,
) -> dict[str, float]:
    axes = parameterization["axes"]
    fixed = parameterization["fixed"]
    pivot = parameterization["eliminated_element"]
    result = dict(fixed)
    for axis, value in zip(axes, point):
        result[axis] = float(value)
    result[pivot] = (
        host_hf
        - sum(host_formula[element] * value for element, value in result.items())
    ) / host_formula[pivot]
    return {element: float(result[element]) for element in elements}


def _active_constraints(point: np.ndarray, constraints: list[dict[str, Any]], tolerance: float = TOL) -> list[str]:
    active = []
    for constraint in constraints:
        lhs = float(np.dot(constraint["coefficients"], point))
        if abs(lhs - constraint["rhs_eV"]) <= tolerance * max(1.0, abs(constraint["rhs_eV"])):
            active.append(constraint["id"])
    return active


def _ordered_indices(points: list[np.ndarray], normal: np.ndarray | None = None) -> list[int]:
    array = np.asarray(points)
    if array.shape[1] == 2:
        center = np.mean(array, axis=0)
        angles = np.arctan2(array[:, 1] - center[1], array[:, 0] - center[0])
        return [int(index) for index in np.argsort(angles)]
    if array.shape[1] == 3 and normal is not None and len(array) >= 3:
        center = np.mean(array, axis=0)
        direction = array[0] - center
        norm = np.linalg.norm(direction)
        if norm <= TOL:
            return list(range(len(points)))
        u = direction / norm
        v = np.cross(normal, u)
        vnorm = np.linalg.norm(v)
        if vnorm <= TOL:
            return list(range(len(points)))
        v /= vnorm
        projected = np.column_stack(((array - center) @ u, (array - center) @ v))
        return _ordered_indices([row for row in projected])
    return list(range(len(points)))


def _geometry(
    vertices: list[np.ndarray],
    constraints: list[dict[str, Any]],
    parameterization: dict[str, Any],
) -> dict[str, Any]:
    dimension = len(parameterization["axes"])
    affine_dimension = _affine_dimension(vertices)
    vertex_rows = []
    for point in vertices:
        vertex_rows.append(
            {
                "coordinates": {
                    element: float(value)
                    for element, value in zip(parameterization["axes"], point)
                },
                "delta_mu_eV": _full_mu(
                    point,
                    parameterization,
                    list(parameterization["host_equality"]["formula"]),
                    parameterization["host_equality"]["formula"],
                    parameterization["host_equality"]["formation_energy_eV_per_formula"],
                ),
                "active_constraints": _active_constraints(point, constraints),
            }
        )
    result: dict[str, Any] = {
        "status": "full_dimensional" if affine_dimension == dimension else "degenerate",
        "ambient_dimension": dimension,
        "affine_dimension": affine_dimension,
        "vertex_count": len(vertex_rows),
        "vertices": vertex_rows,
        "parameterization": parameterization,
    }
    if dimension == 2 and affine_dimension == 2:
        coordinates = [np.asarray([row["coordinates"][axis] for axis in parameterization["axes"]]) for row in vertex_rows]
        result["ordered_vertex_indices"] = _ordered_indices(coordinates)
        edges = []
        order = result["ordered_vertex_indices"]
        for position, first in enumerate(order):
            second = order[(position + 1) % len(order)]
            common = sorted(set(vertex_rows[first]["active_constraints"]) & set(vertex_rows[second]["active_constraints"]))
            edges.append(
                {
                    "vertices": [first, second],
                    "active_constraints": common,
                    "competing_phases": [
                        item["label"]
                        for item in constraints
                        if item["id"] in common and item["id"].startswith("phase:")
                    ],
                }
            )
        result["edges"] = edges
    if dimension == 3 and affine_dimension == 3:
        faces = []
        for constraint in constraints:
            active_indices = [
                index
                for index, point in enumerate(vertices)
                if abs(float(np.dot(constraint["coefficients"], point)) - constraint["rhs_eV"])
                <= TOL * max(1.0, abs(constraint["rhs_eV"]))
            ]
            if len(active_indices) < 3:
                continue
            face_points = np.asarray([vertices[index] for index in active_indices])
            if np.linalg.matrix_rank(face_points[1:] - face_points[0], tol=TOL * 10) < 2:
                continue
            normal = np.asarray(constraint["coefficients"], dtype=float)
            ordered_local = _ordered_indices([vertices[index] for index in active_indices], normal)
            faces.append(
                {
                    "constraint": constraint["id"],
                    "label": constraint["label"],
                    "vertex_indices": [active_indices[index] for index in ordered_local],
                }
            )
        result["faces"] = faces
    return result


def _validate_selected_points(
    raw_points: Any,
    elements: list[str],
    host_formula: dict[str, float],
    host_hf: float,
    constraints: list[dict[str, Any]],
) -> list[dict[str, Any]]:
    if raw_points is None:
        return []
    if not isinstance(raw_points, list):
        raise DomainError("selected_points must be a list of user-specified points")
    result = []
    for index, raw in enumerate(raw_points):
        if not isinstance(raw, dict) or not isinstance(raw.get("delta_mu_eV"), dict):
            raise DomainError(f"selected_points[{index}] requires a delta_mu_eV object")
        values_raw = raw["delta_mu_eV"]
        if set(values_raw) != set(elements):
            raise DomainError(
                f"selected_points[{index}].delta_mu_eV must contain exactly: "
                + ", ".join(elements)
            )
        values = {element: _as_finite_number(values_raw[element], f"selected_points[{index}].{element}") for element in elements}
        fixed_raw = raw.get("fixed_slice_values_eV", {})
        if not isinstance(fixed_raw, dict) or len(fixed_raw) > max(0, len(elements) - 3):
            raise DomainError(
                f"selected_points[{index}].fixed_slice_values_eV contains too many "
                "fixed host elements for a two-dimensional slice"
            )
        fixed_values = {
            element: _as_finite_number(value, f"selected_points[{index}].fixed_slice_values_eV.{element}")
            for element, value in fixed_raw.items()
        }
        if any(element not in elements for element in fixed_values):
            raise DomainError(f"selected_points[{index}] fixes an element outside the host")
        for element, value in fixed_values.items():
            if abs(values[element] - value) > TOL * max(1.0, abs(value)):
                raise DomainError(
                    f"selected point {index} does not match its fixed slice value for {element}"
                )
        host_residual = sum(host_formula[element] * values[element] for element in elements) - host_hf
        violated = []
        margins = {}
        for constraint in constraints:
            lhs = sum(constraint["coefficients"][element] * values[element] for element in elements)
            if constraint.get("equality"):
                continue
            margins[constraint["id"]] = constraint["rhs_eV"] - lhs
            if lhs > constraint["rhs_eV"] + TOL * max(1.0, abs(constraint["rhs_eV"])):
                violated.append(constraint["id"])
        if abs(host_residual) > TOL * max(1.0, abs(host_hf)):
            raise DomainError(f"selected point {index} does not satisfy the host formation-energy equality")
        if violated:
            raise DomainError(f"selected point {index} violates constraints: {', '.join(violated)}")
        name = raw.get("name", f"user_point_{index + 1}")
        if not isinstance(name, str) or not name.strip():
            raise DomainError(f"selected_points[{index}].name must be a non-empty string")
        result.append(
            {
                "name": name.strip(),
                "selection_source": "explicit_user_input",
                "delta_mu_eV": values,
                "constraints_source": [item["id"] for item in constraints],
                "domain_dimension": len(elements) - 1 - len(fixed_values),
                "domain_ambient_dimension": len(elements) - 1 - len(fixed_values),
                "fixed_slice_values_eV": fixed_values,
                "host_equality_residual_eV": host_residual,
                "constraint_margins_eV": margins,
                "minimum_constraint_margin_eV": min(margins.values()) if margins else None,
            }
        )
    return result


def _default_control_element(
    elements: list[str],
) -> tuple[str | None, str, str | None]:
    if "O" in elements:
        return "O", "oxygen_priority", None
    try:
        from pymatgen.core import Element
    except ImportError:
        return (
            None,
            "unresolved_requires_explicit_control_element",
            "automatic halogen detection requires pymatgen",
        )

    halogens = []
    for element in elements:
        try:
            if Element(element).is_halogen:
                halogens.append(element)
        except (ValueError, KeyError):
            continue
    if len(halogens) == 1:
        return halogens[0], "sole_pymatgen_halogen", None
    if len(halogens) > 1:
        return (
            None,
            "unresolved_requires_explicit_control_element",
            "host contains multiple recognized halogens",
        )
    return (
        None,
        "unresolved_requires_explicit_control_element",
        "host has no oxygen or uniquely recognized halogen",
    )


def _point_selection_config(
    raw: Any,
    manual_points: Any,
    elements: list[str],
    domain: dict[str, Any],
    slices: list[dict[str, Any]],
) -> tuple[list[dict[str, Any]], dict[str, Any]]:
    if raw is None:
        raw = {}
    if not isinstance(raw, dict):
        raise DomainError("point_selection must be an object")
    if manual_points is None:
        manual_points = []
    if not isinstance(manual_points, list):
        raise DomainError("selected_points must be a list of complete manual audit points")
    allowed = {"mode", "control_element", "fractions", "fixed_slice_values_eV"}
    extra = set(raw) - allowed
    if extra:
        raise DomainError(f"point_selection has unsupported fields: {', '.join(sorted(extra))}")

    mode = raw.get("mode", "auto_and_manual")
    if not isinstance(mode, str) or mode not in {"auto", "manual", "auto_and_manual"}:
        raise DomainError("point_selection.mode must be auto, manual, or auto_and_manual")
    manual_count = len(manual_points)
    if mode == "manual" and manual_count == 0:
        raise DomainError("manual point-selection mode requires selected_points")
    if mode == "auto" and manual_count:
        raise DomainError(
            "selected_points were supplied with mode=auto; use auto_and_manual or manual"
        )
    auto_enabled = mode in {"auto", "auto_and_manual"}
    if not auto_enabled:
        return [], {
            "mode": mode,
            "selection_status": "manual_points_only",
            "automatic_points_created": 0,
            "manual_points_included": manual_count,
            "default_control_rule": "oxygen_then_unique_pymatgen_halogen",
            "selection_domain": None,
        }

    if "control_element" in raw:
        control_element = raw["control_element"]
        control_element_source = "explicit_input"
        control_element_reason = None
    else:
        (
            control_element,
            control_element_source,
            control_element_reason,
        ) = _default_control_element(elements)
    if control_element is not None and (
        not isinstance(control_element, str) or control_element not in elements
    ):
        raise DomainError(
            "point_selection.control_element must be one of the host elements"
        )
    fractions_raw = raw.get("fractions", list(DEFAULT_POINT_FRACTIONS))
    if not isinstance(fractions_raw, list) or not fractions_raw:
        raise DomainError("point_selection.fractions must be a non-empty list")
    fractions = [
        _as_finite_number(value, f"point_selection.fractions[{index}]")
        for index, value in enumerate(fractions_raw)
    ]
    if any(value < 0 or value > 1 for value in fractions):
        raise DomainError("point_selection fractions must lie in [0, 1]")
    if any(
        abs(left - right) <= TOL
        for index, left in enumerate(fractions)
        for right in fractions[index + 1 :]
    ):
        raise DomainError("point_selection fractions must be unique")

    fixed_raw = raw.get("fixed_slice_values_eV", {})
    if not isinstance(fixed_raw, dict) or len(fixed_raw) > max(0, len(elements) - 3):
        raise DomainError(
            "point_selection.fixed_slice_values_eV must identify an actual 2D slice"
        )
    fixed = {
        element: _as_finite_number(
            value, f"point_selection.fixed_slice_values_eV.{element}"
        )
        for element, value in fixed_raw.items()
    }
    if any(element not in elements for element in fixed):
        raise DomainError("point-selection slice fixes an element outside the host")

    geometry = domain
    selection_domain = {"type": "full_stability_domain"}
    if fixed:
        matching = [
            item
            for item in slices
            if set(item.get("fixed_slice_values_eV", {})) == set(fixed)
            and all(
                math.isclose(
                    item["fixed_slice_values_eV"][element],
                    value,
                    rel_tol=0.0,
                    abs_tol=TOL * max(1.0, abs(value)),
                )
                for element, value in fixed.items()
            )
        ]
        if len(matching) != 1:
            raise DomainError(
                "point-selection slice must match exactly one computed fixed-chemical-potential cross-section"
            )
        geometry = matching[0]
        if geometry.get("status") == "empty" or not geometry.get("vertices"):
            raise DomainError("cannot select points from an empty fixed slice")
        selection_domain = {
            "type": "fixed_chemical_potential_cross_section",
            "fixed_slice_values_eV": fixed,
            "geometry_type": geometry.get("geometry_type"),
            "slice_label": geometry.get("slice_label"),
        }

    if control_element is None:
        return [], {
            "mode": mode,
            "selection_status": "requires_explicit_control_element",
            "automatic_points_created": 0,
            "manual_points_included": manual_count if mode == "auto_and_manual" else 0,
            "control_element": None,
            "control_element_source": control_element_source,
            "control_element_reason": control_element_reason,
            "default_control_rule": "oxygen_then_unique_pymatgen_halogen",
            "fractions_from_poor_to_rich": fractions,
            "control_interval_eV": None,
            "selection_domain": selection_domain,
            "physical_interpretation": (
                "stability domain is calculated; representative points require an "
                "explicit control element"
            ),
        }

    vertices = geometry.get("vertices")
    if not isinstance(vertices, list) or not vertices:
        raise DomainError("cannot select points from a domain without vertices")
    vectors = []
    for index, vertex in enumerate(vertices):
        values = vertex.get("delta_mu_eV")
        if not isinstance(values, dict) or set(values) != set(elements):
            raise DomainError(
                f"domain vertex {index} does not contain the complete host chemical-potential vector"
            )
        vectors.append(
            {element: _as_finite_number(values[element], f"domain.vertices[{index}].{element}") for element in elements}
        )

    control_values = [vector[control_element] for vector in vectors]
    minimum = min(control_values)
    maximum = max(control_values)
    scale = max(1.0, abs(minimum), abs(maximum))
    face_tolerance = TOL * scale
    poor_indices = [
        index for index, value in enumerate(control_values) if abs(value - minimum) <= face_tolerance
    ]
    rich_indices = [
        index for index, value in enumerate(control_values) if abs(value - maximum) <= face_tolerance
    ]
    if math.isclose(minimum, maximum, rel_tol=0.0, abs_tol=face_tolerance):
        return [], {
            "mode": mode,
            "selection_status": "degenerate_control_range",
            "automatic_points_created": 0,
            "manual_points_included": manual_count if mode == "auto_and_manual" else 0,
            "control_element": control_element,
            "control_element_source": control_element_source,
            "default_control_rule": "oxygen_then_unique_pymatgen_halogen",
            "fractions_from_poor_to_rich": fractions,
            "control_interval_eV": {"minimum": minimum, "maximum": maximum},
            "poor_extreme_face_vertex_indices": poor_indices,
            "rich_extreme_face_vertex_indices": rich_indices,
            "selection_domain": selection_domain,
            "control_range_degenerate": True,
            "physical_interpretation": (
                "the control element has zero feasible range; no distinct poor, "
                "middle, and rich points can be selected automatically"
            ),
        }

    def face_representative(indices: list[int]) -> dict[str, float]:
        ordered = sorted(
            (vectors[index] for index in indices),
            key=lambda vector: tuple(vector[element] for element in elements),
        )
        return {
            element: sum(vector[element] for vector in ordered) / len(ordered)
            for element in elements
        }

    poor = face_representative(poor_indices)
    rich = face_representative(rich_indices)
    generated = []
    point_metadata = []
    for fraction in fractions:
        values = {
            element: (1.0 - fraction) * poor[element] + fraction * rich[element]
            for element in elements
        }
        if fraction == 0:
            label = "poor_limit"
        elif fraction == 1:
            label = "rich_limit"
        elif math.isclose(fraction, 0.5, rel_tol=0.0, abs_tol=TOL):
            label = "middle"
        elif fraction < 0.5:
            label = "poor_side"
        else:
            label = "rich_side"
        generated.append(
            {
                "name": f"{control_element}-{label}-f{fraction:g}",
                "delta_mu_eV": values,
                "fixed_slice_values_eV": fixed,
            }
        )
        point_metadata.append(
            {
                "selection_source": "automatic_extreme_face_convex_interpolation",
                "selection_algorithm": (
                    "deterministic_centroids_of_control_element_extreme_face_vertices_"
                    "then_convex_interpolation_from_poor_to_rich"
                ),
                "control_element": control_element,
                "control_element_source": control_element_source,
                "control_min_eV": minimum,
                "control_max_eV": maximum,
                "fraction_from_poor_to_rich": fraction,
                "poor_extreme_face_vertex_indices": poor_indices,
                "rich_extreme_face_vertex_indices": rich_indices,
                "poor_extreme_face_representative_delta_mu_eV": poor,
                "rich_extreme_face_representative_delta_mu_eV": rich,
            }
        )

    return generated, {
        "mode": mode,
        "selection_status": "resolved",
        "automatic_points_created": len(generated),
        "manual_points_included": manual_count if mode == "auto_and_manual" else 0,
        "control_element": control_element,
        "control_element_source": control_element_source,
        "default_control_rule": "oxygen_then_unique_pymatgen_halogen",
        "fractions_from_poor_to_rich": fractions,
        "control_interval_eV": {"minimum": minimum, "maximum": maximum},
        "poor_extreme_face_vertex_indices": poor_indices,
        "rich_extreme_face_vertex_indices": rich_indices,
        "poor_extreme_face_representative_delta_mu_eV": poor,
        "rich_extreme_face_representative_delta_mu_eV": rich,
        "selection_domain": selection_domain,
        "control_range_degenerate": math.isclose(
            minimum, maximum, rel_tol=0.0, abs_tol=face_tolerance
        ),
        "physical_interpretation": (
            "geometric points inside the computed stability domain; not evidence that "
            "the extreme chemical-potential conditions are experimentally accessible"
        ),
        "point_metadata": point_metadata,
    }


def _dopant_binary_reservoir(
    raw: Any,
    host_elements: list[str],
    pbe_reference: str,
    points: list[dict[str, Any]],
    globally_fixed: dict[str, float],
) -> dict[str, Any] | None:
    if raw is None:
        return None
    if not isinstance(raw, dict) or set(raw) != {
        "dopant", "phase_name", "formula", "formation_energy_eV_per_formula",
        "elemental_reference", "source",
    }:
        raise DomainError(
            "dopant_reservoir requires dopant, phase_name, binary formula, "
            "formation_energy_eV_per_formula, elemental_reference, and source"
        )
    dopant = raw["dopant"]
    formula = _formula(raw["formula"], "dopant_reservoir.formula")
    ligand = next((element for element in formula if element != dopant), None)
    if (
        not isinstance(dopant, str)
        or not dopant.strip()
        or dopant in host_elements
        or len(formula) != 2
        or dopant not in formula
        or ligand not in host_elements
    ):
        raise DomainError("dopant_reservoir must be one dopant and one host element")
    if dopant in globally_fixed:
        raise DomainError(
            f"{dopant} cannot be both an explicitly fixed external chemical "
            "potential and an automatically bounded binary reservoir"
        )
    if raw["elemental_reference"] != pbe_reference:
        raise DomainError("dopant_reservoir uses a different PBE elemental reference")
    if not isinstance(raw["phase_name"], str) or not raw["phase_name"].strip():
        raise DomainError("dopant_reservoir.phase_name must be nonempty")
    if not isinstance(raw["source"], str) or not raw["source"].strip():
        raise DomainError("dopant_reservoir.source must identify the corrected binary phase")
    formation_energy = _as_finite_number(
        raw["formation_energy_eV_per_formula"],
        "dopant_reservoir.formation_energy_eV_per_formula",
    )
    for point in points:
        limit = (
            formation_energy - formula[ligand] * point["delta_mu_eV"][ligand]
        ) / formula[dopant]
        selected = min(0.0, limit)
        point["external_fixed_mu_eV"] = {**globally_fixed, dopant: selected}
        point["dopant_reservoir_audit"] = {
            "phase_name": raw["phase_name"],
            "source": raw["source"],
            "dopant": dopant,
            "ligand": ligand,
            "uncapped_binary_limit_eV": limit,
            "selected_delta_mu_eV": selected,
            "active_bound": "elemental" if limit >= 0 else "binary_phase",
            "binary_constraint_margin_eV": (
                formation_energy
                - formula[dopant] * selected
                - formula[ligand] * point["delta_mu_eV"][ligand]
            ),
        }
    return {
        **raw,
        "formula": formula,
        "ligand": ligand,
        "selection_rule": "min(0, (PBE_binary_formation_energy - ligand_count * host_delta_mu) / dopant_count)",
        "host_stability_domain_unchanged": True,
    }


def analyze(data: dict[str, Any]) -> dict[str, Any]:
    if not isinstance(data, dict):
        raise DomainError("input JSON root must be an object")
    host_raw = data.get("host")
    if not isinstance(host_raw, dict):
        raise DomainError("host must be an object")
    host_formula = _formula(host_raw.get("formula"), "host.formula")
    elements = list(host_formula)
    if len(elements) not in (2, 3, 4, 5):
        raise DomainError("host must contain exactly 2, 3, 4, or 5 elements")
    host_hf, phases, energy_provenance = _read_energy_rows(data, host_formula)
    external_mu = _external_mu(
        data, phases, elements, energy_provenance["common_metadata"]["elemental_reference"]
    )
    constraints = _make_full_constraints(host_formula, host_hf, phases, external_mu)

    ambient = len(elements) - 1
    if ambient == 1:
        axes = [elements[0]]
    elif ambient == 2:
        plot = data.get("plot", {})
        axes = plot.get("axes", []) if isinstance(plot, dict) else []
        if not axes:
            axes = elements[-2:]
        if not isinstance(axes, list) or len(axes) != 2:
            raise DomainError("plot.axes must name two host elements for a ternary host")
    else:
        axes = elements[:-1]

    matrix, bounds, reduced_constraints, parameterization = _parameterized_inequalities(
        elements, host_formula, host_hf, constraints, axes
    )
    if ambient == 4:
        solve_matrix, solve_bounds, prefilter = _four_dimensional_halfspace_prefilter(
            matrix, bounds, reduced_constraints
        )
    else:
        solve_matrix, solve_bounds, prefilter = matrix, bounds, None
    vertices = _enumerate_vertices(solve_matrix, solve_bounds)
    base_geometry = _geometry(vertices, reduced_constraints, parameterization)
    if prefilter is not None:
        base_geometry["halfspace_prefilter"] = prefilter
    if base_geometry["ambient_dimension"] != ambient:
        base_geometry["status"] = "degenerate"
    if ambient == 3:
        base_geometry["geometry_type"] = "3d_polytope"
    elif ambient == 4:
        plot = data.get("plot", {})
        if not isinstance(plot, dict):
            raise DomainError("plot must be an object for a five-element host")
        projection_axes = plot.get("projection_axes", elements[:2])
        if (
            not isinstance(projection_axes, list)
            or len(projection_axes) != 2
            or any(not isinstance(axis, str) for axis in projection_axes)
            or len(set(projection_axes)) != 2
            or any(axis not in elements for axis in projection_axes)
        ):
            raise DomainError("plot.projection_axes must name two distinct host elements")
        base_geometry["geometry_type"] = "4d_stability_domain"
        base_geometry["projection"] = _two_dimensional_projection(
            base_geometry["vertices"], projection_axes
        )
    elif ambient == 2:
        base_geometry["geometry_type"] = "2d_stability_domain"
    else:
        base_geometry["geometry_type"] = "1d_stability_interval"

    slices_result = []
    raw_slices = data.get("slices", [])
    if not isinstance(raw_slices, list):
        raise DomainError("slices must be a list")
    for slice_index, spec in enumerate(raw_slices):
        if ambient not in (3, 4):
            raise DomainError(
                "fixed-element 2D slices are supported for four- and five-element hosts"
            )
        fixed_count = ambient - 2
        if (
            not isinstance(spec, dict)
            or not isinstance(spec.get("fixed"), dict)
            or len(spec["fixed"]) != fixed_count
        ):
            raise DomainError(
                f"slices[{slice_index}].fixed must contain exactly "
                f"{'one' if fixed_count == 1 else 'two'} element/value"
            )
        fixed = {
            element: _as_finite_number(value, f"slices[{slice_index}].fixed.{element}")
            for element, value in spec["fixed"].items()
        }
        slice_axes = spec.get("axes")
        if not isinstance(slice_axes, list) or len(slice_axes) != 2:
            raise DomainError(f"slices[{slice_index}].axes must contain two host elements")
        smat, sbounds, sconstraints, sparam = _parameterized_inequalities(
            elements, host_formula, host_hf, constraints, slice_axes, fixed
        )
        try:
            slice_vertices = _enumerate_vertices(smat, sbounds)
            slice_geometry = _geometry(slice_vertices, sconstraints, sparam)
            slice_geometry["geometry_type"] = "2d_fixed_chemical_potential_slice"
            slice_geometry["slice_label"] = "fixed-chemical-potential cross-section, not projection"
        except DomainError as exc:
            slice_geometry = {
                "status": "empty",
                "geometry_type": "2d_fixed_chemical_potential_slice",
                "slice_label": "fixed-chemical-potential cross-section, not projection",
                "fixed_slice_values_eV": fixed,
                "axes": slice_axes,
                "error": str(exc),
                "ambient_dimension": 2,
                "affine_dimension": None,
                "vertex_count": 0,
                "vertices": [],
            }
        slice_geometry["fixed_slice_values_eV"] = fixed
        slices_result.append(slice_geometry)

    automatic_raw, point_selection_audit = _point_selection_config(
        data.get("point_selection", {}),
        data.get("selected_points", []),
        elements,
        base_geometry,
        slices_result,
    )
    automatic_points = _validate_selected_points(
        automatic_raw, elements, host_formula, host_hf, constraints
    )
    for point, metadata in zip(
        automatic_points, point_selection_audit.get("point_metadata", []), strict=True
    ):
        point.update(metadata)
    manual_raw = data.get("selected_points", [])
    if point_selection_audit["mode"] == "auto":
        manual_raw = []
    selected_points = automatic_points + _validate_selected_points(
        manual_raw, elements, host_formula, host_hf, constraints
    )
    names = [item["name"] for item in selected_points]
    if len(names) != len(set(names)):
        raise DomainError("selected point names must be unique")
    for item in selected_points:
        item["external_fixed_mu_eV"] = external_mu
        if item["fixed_slice_values_eV"]:
            matching_slice = next(
                (
                    geometry
                    for geometry in slices_result
                    if geometry.get("fixed_slice_values_eV") == item["fixed_slice_values_eV"]
                ),
                None,
            )
            if matching_slice is None:
                raise DomainError(
                    f"selected point {item['name']!r} declares a fixed slice that is "
                    "not present in slices"
                )
            item["domain_dimension"] = matching_slice["affine_dimension"]
            item["domain_ambient_dimension"] = matching_slice["ambient_dimension"]
        else:
            item["domain_dimension"] = base_geometry["affine_dimension"]
            item["domain_ambient_dimension"] = base_geometry["ambient_dimension"]

    binary_reservoir = _dopant_binary_reservoir(
        data.get("dopant_reservoir"),
        elements,
        energy_provenance["common_metadata"]["elemental_reference"],
        selected_points,
        external_mu,
    )
    return {
        "schema_version": 1,
        "status": "accepted",
        "origin": data.get("origin", {}),
        "chemical_potential_convention": "delta_mu_i = mu_i - mu_i(elemental reference), in eV",
        "chemical_potential_functional": "PBE",
        "absolute_chemical_potentials_reported": False,
        "downstream_absolute_mu_contract": {
            "constructed_by_this_script": False,
            "reference_functional": "HSE06",
            "rule": "mu_i_absolute_HSE06 = mu_i_reference_HSE06 + delta_mu_i_PBE",
        },
        "host": {
            "formula": host_formula,
            "elements": elements,
            "formation_energy_eV_per_formula": host_hf,
        },
        "energy_provenance": energy_provenance,
        "external_fixed_mu_eV": external_mu,
        "dopant_reservoir": binary_reservoir,
        "dimension": {
            "chemical_potential_dimension": ambient,
            "affine_dimension": base_geometry["affine_dimension"],
            "status": base_geometry["status"],
        },
        "constraints": constraints,
        "domain": base_geometry,
        "slices": slices_result,
        "selected_points": selected_points,
        "point_selection": {
            key: value
            for key, value in point_selection_audit.items()
            if key != "point_metadata"
        },
        "representative_point_selected_automatically": (
            point_selection_audit["automatic_points_created"] > 0
        ),
        "tolerance": {"absolute": TOL, "relative": TOL},
    }


def _plot_geometry(analysis: dict[str, Any], geometry: dict[str, Any], output: Path, title: str) -> list[Path]:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from mpl_toolkits.mplot3d.art3d import Poly3DCollection

    output.mkdir(parents=True, exist_ok=True)
    written: list[Path] = []
    dimension = geometry["ambient_dimension"]
    vertices = geometry["vertices"]
    fig = plt.figure(figsize=(6.4, 5.2), constrained_layout=True)
    if dimension == 4:
        ax = fig.add_subplot(111)
        projection = geometry["projection"]
        axes = projection["axes"]
        hull = np.asarray(projection["hull_coordinates_eV"], dtype=float)
        points = np.asarray([
            [vertex["delta_mu_eV"][axis] for axis in axes] for vertex in vertices
        ])
        if len(hull) >= 3:
            ax.fill(hull[:, 0], hull[:, 1], color="#b9dce8", alpha=0.68)
        if len(hull) >= 2:
            boundary = np.vstack((hull, hull[0])) if len(hull) >= 3 else hull
            ax.plot(boundary[:, 0], boundary[:, 1], color="#167d9a", linewidth=1.6)
        ax.scatter(points[:, 0], points[:, 1], color="#202b35", s=17)
        ax.set_xlabel(f"$\\Delta\\mu_{{{axes[0]}}}$ (eV)")
        ax.set_ylabel(f"$\\Delta\\mu_{{{axes[1]}}}$ (eV)")
        ax.set_aspect("equal", adjustable="box")
        ax.grid(color="#d8dde0", linestyle=":", linewidth=0.7)
        ax.text(
            0.02, 0.02,
            f"2D projection of 4D ambient domain (affine dim "
            f"{geometry['affine_dimension']})",
            transform=ax.transAxes, fontsize=8, va="bottom",
        )
    elif dimension == 3 and geometry.get("affine_dimension") == 3:
        ax = fig.add_subplot(111, projection="3d")
        axes = geometry["parameterization"]["axes"]
        points = np.asarray([[row["coordinates"][element] for element in axes] for row in vertices])
        for face in geometry.get("faces", []):
            polygon = points[face["vertex_indices"]]
            ax.add_collection3d(
                Poly3DCollection([polygon], alpha=0.24, edgecolor="black", linewidth=0.7)
            )
        ax.scatter(points[:, 0], points[:, 1], points[:, 2], color="#c43c39", s=18)
        ax.set_xlabel(f"$\\Delta\\mu_{{{axes[0]}}}$ (eV)")
        ax.set_ylabel(f"$\\Delta\\mu_{{{axes[1]}}}$ (eV)")
        ax.set_zlabel(f"$\\Delta\\mu_{{{axes[2]}}}$ (eV)")
    elif dimension == 2 and geometry.get("affine_dimension") == 2:
        ax = fig.add_subplot(111)
        axes = geometry["parameterization"]["axes"]
        order = geometry["ordered_vertex_indices"]
        points = np.asarray(
            [[vertices[index]["coordinates"][axis] for axis in axes] for index in order]
        )
        ax.fill(points[:, 0], points[:, 1], color="#b9dce8", alpha=0.68, zorder=1)
        ax.plot(
            np.r_[points[:, 0], points[0, 0]],
            np.r_[points[:, 1], points[0, 1]],
            color="#167d9a",
            linewidth=1.6,
            zorder=2,
        )
        ax.scatter(points[:, 0], points[:, 1], color="#202b35", s=17, zorder=3)
        for edge in geometry.get("edges", []):
            names = edge["competing_phases"]
            if not names:
                continue
            first, second = edge["vertices"]
            p1 = [vertices[first]["coordinates"][axis] for axis in axes]
            p2 = [vertices[second]["coordinates"][axis] for axis in axes]
            midpoint = (np.asarray(p1) + np.asarray(p2)) / 2
            ax.annotate(", ".join(names), midpoint, fontsize=8, xytext=(3, 3), textcoords="offset points")
        ax.set_xlabel(f"$\\Delta\\mu_{{{axes[0]}}}$ (eV)")
        ax.set_ylabel(f"$\\Delta\\mu_{{{axes[1]}}}$ (eV)")
        ax.set_aspect("equal", adjustable="box")
        ax.grid(color="#d8dde0", linestyle=":", linewidth=0.7)
    elif dimension == 2:
        ax = fig.add_subplot(111)
        axes = geometry["parameterization"]["axes"]
        points = np.asarray(
            [[row["coordinates"][axis] for axis in axes] for row in vertices]
        )
        if geometry["affine_dimension"] == 1:
            direction = points[-1] - points[0]
            ordered = sorted(
                points, key=lambda point: float(np.dot(point - points[0], direction))
            )
            ax.plot(
                [point[0] for point in ordered],
                [point[1] for point in ordered],
                color="#167d9a", linewidth=1.6,
            )
        ax.scatter(points[:, 0], points[:, 1], color="#202b35", s=20, zorder=3)
        ax.set_xlabel(f"$\\Delta\\mu_{{{axes[0]}}}$ (eV)")
        ax.set_ylabel(f"$\\Delta\\mu_{{{axes[1]}}}$ (eV)")
        ax.set_aspect("equal", adjustable="box")
        ax.grid(color="#d8dde0", linestyle=":", linewidth=0.7)
    elif dimension == 1:
        ax = fig.add_subplot(111)
        element = geometry["parameterization"]["axes"][0]
        values = sorted(row["coordinates"][element] for row in vertices)
        if len(values) == 1:
            ax.scatter([values[0]], [0], color="#167d9a")
        else:
            ax.plot(values, [0] * len(values), color="#167d9a", linewidth=2)
            ax.scatter(values, [0] * len(values), color="#202b35", zorder=3)
        ax.set_xlabel(f"$\\Delta\\mu_{{{element}}}$ (eV)")
        ax.set_yticks([])
        ax.grid(axis="x", color="#d8dde0", linestyle=":", linewidth=0.7)
    elif dimension == 3:
        ax = fig.add_subplot(111, projection="3d")
        axes = geometry["parameterization"]["axes"]
        points = np.asarray(
            [[row["coordinates"][element] for element in axes] for row in vertices]
        )
        if geometry["affine_dimension"] == 2 and len(points) >= 3:
            _, _, basis = np.linalg.svd(points - points.mean(axis=0))
            projected = (points - points.mean(axis=0)) @ basis[:2].T
            order = _ordered_indices([point for point in projected])
            ax.add_collection3d(
                Poly3DCollection([points[order]], alpha=0.24, edgecolor="black")
            )
        elif geometry["affine_dimension"] == 1 and len(points) >= 2:
            direction = points[-1] - points[0]
            ordered = sorted(
                points, key=lambda point: float(np.dot(point - points[0], direction))
            )
            ax.plot(*np.asarray(ordered).T, color="#167d9a", linewidth=1.6)
        ax.scatter(points[:, 0], points[:, 1], points[:, 2], color="#c43c39", s=20)
        ax.set_xlabel(f"$\\Delta\\mu_{{{axes[0]}}}$ (eV)")
        ax.set_ylabel(f"$\\Delta\\mu_{{{axes[1]}}}$ (eV)")
        ax.set_zlabel(f"$\\Delta\\mu_{{{axes[2]}}}$ (eV)")
    else:
        ax = fig.add_subplot(111)
        ax.text(0.5, 0.5, f"Degenerate domain (affine dimension {geometry.get('affine_dimension')})", ha="center", va="center")
        ax.set_axis_off()
    # Relative chemical potentials cannot exceed their elemental reference.
    # Fix the displayed upper bounds as well: Matplotlib's default margins
    # otherwise extend the picture into physically forbidden positive values.
    ax.set_xlim(right=0.0)
    if dimension >= 2:
        ax.set_ylim(top=0.0)
    if dimension == 3:
        ax.set_zlim(top=0.0)
    ax.set_title(title, fontsize=11)
    for suffix in ("svg", "png"):
        path = output / f"phase_diagram.{suffix}"
        fig.savefig(path, dpi=220)
        written.append(path)
    plt.close(fig)
    return written


def write_outputs(analysis: dict[str, Any], output_dir: Path) -> list[Path]:
    if output_dir.exists():
        if not output_dir.is_dir() or any(output_dir.iterdir()):
            raise DomainError(f"output directory must be new or empty; refusing to overwrite: {output_dir}")
    else:
        output_dir.mkdir(parents=True)
    paths = [output_dir / "analysis.json"]
    paths[0].write_text(json.dumps(analysis, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    title = "".join(
        f"{element}{count:g}" if count != 1 else element
        for element, count in analysis["host"]["formula"].items()
    )
    paths.extend(_plot_geometry(analysis, analysis["domain"], output_dir / "domain", f"Host stability domain: {title}"))
    for index, geometry in enumerate(analysis["slices"]):
        if geometry.get("status") == "empty":
            continue
        paths.extend(
            _plot_geometry(
                analysis,
                geometry,
                output_dir / f"slice_{index + 1}",
                f"Fixed chemical-potential slice: {geometry['fixed_slice_values_eV']}",
            )
        )
    return paths


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        description="Build a 2–5 element chemical-potential stability domain."
    )
    parser.add_argument("--input", required=True, help="Explicit JSON energy/stoichiometry input")
    parser.add_argument("--output-dir", required=True, help="New output directory; existing files are never overwritten")
    args = parser.parse_args(argv)
    try:
        input_path = Path(args.input)
        if not input_path.is_file():
            raise DomainError(f"input file does not exist: {input_path}")
        data = json.loads(input_path.read_text(encoding="utf-8"))
        analysis = analyze(data)
        written = write_outputs(analysis, Path(args.output_dir))
    except (DomainError, OSError, json.JSONDecodeError) as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 2
    for path in written:
        print(f"Written: {path}")
    print(
        f"Accepted {analysis['host']['formula']}: "
        f"{analysis['dimension']['affine_dimension']}D stability domain"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
