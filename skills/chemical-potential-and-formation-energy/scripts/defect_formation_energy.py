#!/usr/bin/env python3
"""Plot audited defect formation-energy lower envelopes without doped."""

from __future__ import annotations

import argparse
import json
import math
import sys
from collections import defaultdict
from datetime import datetime, timezone, timedelta
from pathlib import Path
from typing import Any, Sequence

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np


SCHEMA_VERSION = 1
_NUMBER_TOL = 1e-8
_ENERGY_TOL = 1e-9
_X_TOL = 1e-10
_PALETTE = ("#e31a99", "#f58220", "#a64b00", "#32b94a", "#00b9c7")


def _fail(message: str) -> None:
    raise ValueError(message)


def _number(value: Any, path: str) -> float:
    if isinstance(value, bool) or not isinstance(value, int | float):
        _fail(f"{path} must be a finite number")
    result = float(value)
    if not math.isfinite(result):
        _fail(f"{path} must be a finite number")
    return result


def _integer(value: Any, path: str) -> int:
    number = _number(value, path)
    rounded = round(number)
    if abs(number - rounded) > _NUMBER_TOL:
        _fail(f"{path} must be an integer")
    return int(rounded)


def _text(value: Any, path: str) -> str:
    if not isinstance(value, str) or not value.strip():
        _fail(f"{path} must be a non-empty string")
    return value.strip()


def _lattice(value: Any, path: str) -> list[list[float]]:
    if not isinstance(value, list) or len(value) != 3:
        _fail(f"{path} must be a 3x3 lattice matrix in angstrom")
    matrix = []
    for row_index, row in enumerate(value):
        if not isinstance(row, list) or len(row) != 3:
            _fail(f"{path} must be a 3x3 lattice matrix in angstrom")
        matrix.append([_number(cell, f"{path}[{row_index}]") for cell in row])
    if abs(float(np.linalg.det(np.asarray(matrix)))) < 1e-8:
        _fail(f"{path} must describe a non-singular supercell")
    return matrix


def _close(a: float, b: float, *, atol: float = 1e-6) -> bool:
    return math.isclose(a, b, rel_tol=1e-7, abs_tol=atol)


def _same_lattice(a: list[list[float]], b: list[list[float]]) -> bool:
    return bool(np.allclose(a, b, rtol=1e-7, atol=1e-5))


def _duplicate_rejecting_object(pairs: list[tuple[str, Any]]) -> dict[str, Any]:
    result: dict[str, Any] = {}
    for key, value in pairs:
        if key in result:
            _fail(f"duplicate JSON key: {key}")
        result[key] = value
    return result


def _reject_nonstandard_constant(value: str) -> None:
    _fail(f"non-finite JSON constant is not allowed: {value}")


def load_json(path: str | Path) -> dict[str, Any]:
    try:
        data = json.loads(
            Path(path).read_text(encoding="utf-8"),
            object_pairs_hook=_duplicate_rejecting_object,
            parse_constant=_reject_nonstandard_constant,
        )
    except json.JSONDecodeError as exc:
        raise ValueError(f"invalid JSON: {exc}") from exc
    if not isinstance(data, dict):
        _fail("input root must be a JSON object")
    return data


def _parse_chemical_potentials(
    raw: Any, calculation_method: str
) -> tuple[dict[str, float], dict[str, Any]]:
    if not isinstance(raw, dict):
        _fail("chemical_potentials must be an object")
    mode = raw.get("mode")
    if mode == "absolute":
        required = {"mode", "absolute_method", "values_ev"}
        if set(raw) != required:
            _fail(
                "absolute chemical_potentials must contain exactly mode, "
                "absolute_method, and values_ev"
            )
        method = _text(raw["absolute_method"], "chemical_potentials.absolute_method")
        if method != calculation_method:
            _fail("absolute chemical potentials must use the total-energy method")
        values = raw["values_ev"]
        if not isinstance(values, dict):
            _fail("chemical_potentials.values_ev must be an object")
        parsed = {
            _text(element, "chemical-potential element"): _number(
                value, f"chemical_potentials.values_ev.{element}"
            )
            for element, value in values.items()
        }
        return parsed, {
            "mode": mode,
            "absolute_method": method,
            "values_ev": parsed,
            "resolved_absolute_mu_ev": parsed,
        }

    if mode == "relative_plus_elemental_reference":
        required = {
            "mode",
            "relative_method",
            "elemental_reference_method",
            "values",
        }
        if set(raw) != required:
            _fail(
                "relative chemical potentials must explicitly pair relative_method, "
                "elemental_reference_method, and values"
            )
        relative_method = _text(
            raw["relative_method"], "chemical_potentials.relative_method"
        )
        reference_method = _text(
            raw["elemental_reference_method"],
            "chemical_potentials.elemental_reference_method",
        )
        if reference_method != calculation_method:
            _fail(
                "elemental reference method must match the defect and host "
                "total-energy method"
            )
        values = raw["values"]
        if not isinstance(values, dict):
            _fail("chemical_potentials.values must be an object")
        parsed: dict[str, float] = {}
        components: dict[str, dict[str, float]] = {}
        for element, entry in values.items():
            element = _text(element, "chemical-potential element")
            if not isinstance(entry, dict) or set(entry) != {
                "delta_mu_ev",
                "elemental_reference_ev",
            }:
                _fail(
                    f"chemical_potentials.values.{element} must contain exactly "
                    "delta_mu_ev and elemental_reference_ev"
                )
            delta = _number(
                entry["delta_mu_ev"],
                f"chemical_potentials.values.{element}.delta_mu_ev",
            )
            reference = _number(
                entry["elemental_reference_ev"],
                f"chemical_potentials.values.{element}.elemental_reference_ev",
            )
            parsed[element] = reference + delta
            components[element] = {
                "delta_mu_ev": delta,
                "elemental_reference_ev": reference,
                "resolved_absolute_mu_ev": parsed[element],
            }
        return parsed, {
            "mode": mode,
            "relative_method": relative_method,
            "elemental_reference_method": reference_method,
            "values": components,
            "resolved_absolute_mu_ev": parsed,
        }

    _fail(
        "chemical_potentials.mode must be 'absolute' or "
        "'relative_plus_elemental_reference'"
    )


def validate_input(raw: dict[str, Any]) -> dict[str, Any]:
    if (
        isinstance(raw.get("schema_version"), bool)
        or raw.get("schema_version") != SCHEMA_VERSION
    ):
        _fail(f"schema_version must be {SCHEMA_VERSION}")

    calculation_method = _text(raw.get("calculation_method"), "calculation_method")
    if "HSE06" not in calculation_method.upper():
        _fail("calculation_method must identify HSE06")
    band = raw.get("pristine_band_edges")
    if not isinstance(band, dict):
        _fail("pristine_band_edges must be an object")

    material_id = _text(band.get("material_id"), "pristine_band_edges.material_id")
    band_method = _text(
        band.get("calculation_method"),
        "pristine_band_edges.calculation_method",
    )
    if band_method != calculation_method:
        _fail("pristine band-edge method must match calculation_method")

    lattice = _lattice(
        band.get("supercell_lattice_angstrom"),
        "pristine_band_edges.supercell_lattice_angstrom",
    )
    kpoints = band.get("kpoints")
    if not isinstance(kpoints, dict):
        _fail("pristine_band_edges.kpoints must be an object")
    if kpoints.get("scheme") != "Gamma" or kpoints.get("mesh") != [1, 1, 1]:
        _fail("pristine band edges must come from a Gamma-only 1x1x1 calculation")
    if band.get("calculation_type") != "self_consistent":
        _fail("pristine band edges must come from a self-consistent calculation")
    if band.get("converged") is not True:
        _fail("pristine band-edge calculation must be explicitly converged")
    _text(band.get("provenance"), "pristine_band_edges.provenance")
    vbm = _number(band.get("vbm_reference_ev"), "pristine_band_edges.vbm_reference_ev")
    cbm = _number(band.get("cbm_reference_ev"), "pristine_band_edges.cbm_reference_ev")
    gap = _number(band.get("gap_ev"), "pristine_band_edges.gap_ev")
    if gap <= 0:
        _fail("gap_ev must be positive")
    if not _close(cbm - vbm, gap, atol=1e-5):
        _fail("gap_ev must equal cbm_reference_ev - vbm_reference_ev")

    mu, chemical_potential_audit = _parse_chemical_potentials(
        raw.get("chemical_potentials"), calculation_method
    )
    raw_states = raw.get("states")
    if not isinstance(raw_states, list) or not raw_states:
        _fail("states must be a non-empty list")

    states: list[dict[str, Any]] = []
    state_ids: set[str] = set()
    host_energies: list[float] = []
    required_state_fields = {
        "state_id",
        "defect_id",
        "geometry_id",
        "charge",
        "E_def_ev",
        "E_host_ev",
        "n_i",
        "E_corr_ev",
        "material_id",
        "calculation_method",
        "supercell_lattice_angstrom",
    }
    for index, item in enumerate(raw_states):
        path = f"states[{index}]"
        if not isinstance(item, dict):
            _fail(f"{path} must be an object")
        missing = required_state_fields - item.keys()
        if missing:
            _fail(f"{path} missing required fields: {', '.join(sorted(missing))}")
        state_id = _text(item["state_id"], f"{path}.state_id")
        if state_id in state_ids:
            _fail(f"state_id must be unique: {state_id}")
        state_ids.add(state_id)
        defect_id = _text(item["defect_id"], f"{path}.defect_id")
        geometry_id = _text(item["geometry_id"], f"{path}.geometry_id")
        state_material = _text(item["material_id"], f"{path}.material_id")
        state_method = _text(item["calculation_method"], f"{path}.calculation_method")
        state_lattice = _lattice(
            item["supercell_lattice_angstrom"],
            f"{path}.supercell_lattice_angstrom",
        )
        if state_material != material_id:
            _fail(f"{path} material_id is inconsistent with pristine band edges")
        if state_method != calculation_method:
            _fail(f"{path} calculation_method is inconsistent")
        if not _same_lattice(state_lattice, lattice):
            _fail(f"{path} supercell lattice is inconsistent with pristine band edges")

        n_i = item["n_i"]
        if not isinstance(n_i, dict):
            _fail(f"{path}.n_i must map species to signed atom counts")
        parsed_n_i: dict[str, int] = {}
        for element, count in n_i.items():
            element = _text(element, f"{path}.n_i species")
            parsed_n_i[element] = _integer(count, f"{path}.n_i.{element}")
            if element not in mu:
                _fail(f"no chemical potential supplied for {element}")

        charge = _integer(item["charge"], f"{path}.charge")
        e_def = _number(item["E_def_ev"], f"{path}.E_def_ev")
        e_host = _number(item["E_host_ev"], f"{path}.E_host_ev")
        e_corr = _number(item["E_corr_ev"], f"{path}.E_corr_ev")
        host_energies.append(e_host)
        chemical_term = sum(parsed_n_i[element] * mu[element] for element in parsed_n_i)
        intercept = e_def - e_host + e_corr - chemical_term + charge * vbm
        states.append(
            {
                "state_id": state_id,
                "defect_id": defect_id,
                "geometry_id": geometry_id,
                "charge": charge,
                "E_def_ev": e_def,
                "E_host_ev": e_host,
                "n_i": parsed_n_i,
                "E_corr_ev": e_corr,
                "chemical_potential_term_ev": chemical_term,
                "vbm_reference_ev": vbm,
                "formation_energy_intercept_ev": intercept,
                "slope_charge_e": charge,
                "material_id": state_material,
                "calculation_method": state_method,
                "supercell_lattice_angstrom": state_lattice,
            }
        )

    if any(not _close(value, host_energies[0], atol=1e-5) for value in host_energies[1:]):
        _fail("all states must use the same pristine E_host_ev")

    return {
        "schema_version": SCHEMA_VERSION,
        "material_id": material_id,
        "calculation_method": calculation_method,
        "supercell_lattice_angstrom": lattice,
        "band_edges": {
            "vbm_reference_ev": vbm,
            "cbm_reference_ev": cbm,
            "gap_ev": gap,
            "provenance": band["provenance"],
            "calculation_method": band_method,
            "calculation_type": band["calculation_type"],
            "kpoints": {"scheme": "Gamma", "mesh": [1, 1, 1]},
            "converged": True,
        },
        "chemical_potentials": chemical_potential_audit,
        "states": states,
        "origin": raw.get("origin", {}),
    }


def _energy(state: dict[str, Any], fermi_level_ev: float) -> float:
    return state["formation_energy_intercept_ev"] + state["charge"] * fermi_level_ev


def _unique_sorted(values: list[float], tolerance: float = _X_TOL) -> list[float]:
    result: list[float] = []
    for value in sorted(values):
        if not result or abs(value - result[-1]) > tolerance:
            result.append(value)
    return result


def _active_states(
    states: list[dict[str, Any]], fermi_level_ev: float
) -> list[dict[str, Any]]:
    energies = [_energy(state, fermi_level_ev) for state in states]
    minimum = min(energies)
    return [
        state
        for state, energy in zip(states, energies, strict=True)
        if energy <= minimum + _ENERGY_TOL
    ]


def _crossing_points(
    states: list[dict[str, Any]], gap: float
) -> list[dict[str, Any]]:
    result = []
    for left_index, left in enumerate(states):
        for right in states[left_index + 1 :]:
            q_left = left["charge"]
            q_right = right["charge"]
            if q_left == q_right:
                continue
            x = (
                right["formation_energy_intercept_ev"]
                - left["formation_energy_intercept_ev"]
            ) / (q_left - q_right)
            if -_X_TOL <= x <= gap + _X_TOL:
                x = min(gap, max(0.0, x))
                active = _active_states(states, x)
                minimum = min(_energy(state, x) for state in states)
                result.append(
                    {
                        "state_ids": [left["state_id"], right["state_id"]],
                        "charge_pair": [q_left, q_right],
                        "fermi_level_ev": x,
                        "formation_energy_ev": minimum,
                        "in_gap": 0.0 < x < gap,
                        "on_lower_envelope": left in active and right in active,
                    }
                )
    return sorted(result, key=lambda point: (point["fermi_level_ev"], point["state_ids"]))


def _defect_envelope(states: list[dict[str, Any]], gap: float) -> dict[str, Any]:
    defect_id = states[0]["defect_id"]
    candidates = [0.0, gap]
    for left_index, left in enumerate(states):
        for right in states[left_index + 1 :]:
            q_left = left["charge"]
            q_right = right["charge"]
            if q_left == q_right:
                continue
            x = (
                right["formation_energy_intercept_ev"]
                - left["formation_energy_intercept_ev"]
            ) / (q_left - q_right)
            if _X_TOL < x < gap - _X_TOL:
                candidates.append(x)
    breakpoints = _unique_sorted(candidates)

    raw_intervals: list[dict[str, Any]] = []
    for left_x, right_x in zip(breakpoints[:-1], breakpoints[1:], strict=True):
        if right_x - left_x <= _X_TOL:
            continue
        midpoint = (left_x + right_x) / 2.0
        winners = _active_states(states, midpoint)
        state_ids = sorted(state["state_id"] for state in winners)
        charges = sorted({state["charge"] for state in winners})
        geometries = sorted({state["geometry_id"] for state in winners})
        interval = {
            "fermi_level_start_ev": left_x,
            "fermi_level_end_ev": right_x,
            "formation_energy_start_ev": min(
                _energy(state, left_x) for state in states
            ),
            "formation_energy_end_ev": min(
                _energy(state, right_x) for state in states
            ),
            "state_ids": state_ids,
            "charge_states": charges,
            "geometry_ids": geometries,
            "slope_charge_e": charges[0],
        }
        if raw_intervals and raw_intervals[-1]["state_ids"] == state_ids:
            raw_intervals[-1]["fermi_level_end_ev"] = right_x
            raw_intervals[-1]["formation_energy_end_ev"] = interval[
                "formation_energy_end_ev"
            ]
        else:
            raw_intervals.append(interval)

    transitions = []
    for left, right in zip(raw_intervals[:-1], raw_intervals[1:], strict=True):
        x = left["fermi_level_end_ev"]
        if abs(x - right["fermi_level_start_ev"]) > _X_TOL:
            continue
        left_charges = left["charge_states"]
        right_charges = right["charge_states"]
        if set(left_charges) == set(right_charges):
            continue
        active = _active_states(states, x)
        active_charges = sorted({state["charge"] for state in active})
        charge_pairs = sorted(
            {
                (q_left, q_right)
                for q_left in left_charges
                for q_right in right_charges
                if q_left != q_right
            }
        )
        transitions.append(
            {
                "fermi_level_ev": x,
                "formation_energy_ev": min(_energy(state, x) for state in states),
                "left_charge_states": left_charges,
                "right_charge_states": right_charges,
                "charge_pairs": [list(pair) for pair in charge_pairs],
                "state_ids_at_crossing": sorted(
                    state["state_id"] for state in active
                ),
                "in_gap": 0.0 < x < gap,
                "negative_u_jump": any(
                    abs(q_right - q_left) > 1 for q_left, q_right in charge_pairs
                ),
            }
        )

    return {
        "defect_id": defect_id,
        "state_ids": [state["state_id"] for state in states],
        "charge_states_present": sorted({state["charge"] for state in states}),
        "geometry_ids_present": sorted({state["geometry_id"] for state in states}),
        "segments": raw_intervals,
        "crossing_points": _crossing_points(states, gap),
        "transition_levels": [item for item in transitions if item["in_gap"]],
    }


def analyze(raw: dict[str, Any]) -> dict[str, Any]:
    data = validate_input(raw)
    grouped: dict[str, list[dict[str, Any]]] = defaultdict(list)
    for state in data["states"]:
        grouped[state["defect_id"]].append(state)
    envelopes = [
        _defect_envelope(grouped[defect_id], data["band_edges"]["gap_ev"])
        for defect_id in sorted(grouped)
    ]
    return {
        "schema_version": SCHEMA_VERSION,
        "material_id": data["material_id"],
        "calculation_method": data["calculation_method"],
        "band_edges": data["band_edges"],
        "chemical_potentials": data["chemical_potentials"],
        "origin": data["origin"],
        "states": data["states"],
        "defects": envelopes,
    }


def _y_bounds(
    defects: list[dict[str, Any]], requested_max: float | None
) -> tuple[float, float]:
    values = [
        energy
        for defect in defects
        for segment in defect["segments"]
        for energy in (
            segment["formation_energy_start_ev"],
            segment["formation_energy_end_ev"],
        )
    ]
    if not values:
        _fail("no formation-energy envelope could be constructed")
    low = min(values)
    high = max(values)
    energy_span = high - low
    # A flat neutral envelope must not yield a visually exaggerated 0.03 eV
    # y-axis from padding alone; retain at least 0.5 eV unless overridden.
    padding = max(
        0.06 * max(energy_span, 0.25),
        0.5 * max(0.5 - energy_span, 0.0),
    )
    bottom = low - padding
    top = high + padding
    if requested_max is not None:
        requested_max = _number(requested_max, "y_max_ev")
        if requested_max <= bottom:
            _fail("y_max_ev would hide every low-energy envelope")
        top = min(top, requested_max)
    if top <= bottom:
        _fail("invalid y-axis range")
    return bottom, top


def _plot_defects(
    audit: dict[str, Any],
    defects: list[dict[str, Any]],
    output_path: Path,
    *,
    y_max_ev: float | None = None,
    title: str | None = None,
) -> tuple[float, float]:
    gap = audit["band_edges"]["gap_ev"]
    # Automatic limits are local to the plotted group. An explicit cutoff
    # instead uses one common lower bound, so clipping a high-only group
    # cannot fail midway through the output batch.
    display_defects = audit["defects"] if y_max_ev is not None else defects
    y_bounds = _y_bounds(display_defects, y_max_ev)
    fig, ax = plt.subplots(figsize=(4.4, 3.4), constrained_layout=True)
    for index, defect in enumerate(defects):
        color = _PALETTE[index % len(_PALETTE)]
        segments = defect["segments"]
        if not segments:
            continue
        x_values = [segments[0]["fermi_level_start_ev"]]
        y_values = [segments[0]["formation_energy_start_ev"]]
        for segment in segments:
            x_values.append(segment["fermi_level_end_ev"])
            y_values.append(segment["formation_energy_end_ev"])
        ax.plot(
            x_values,
            y_values,
            color=color,
            linewidth=1.55,
            label=defect["defect_id"],
            solid_capstyle="round",
        )
        for transition in defect["transition_levels"]:
            ax.plot(
                transition["fermi_level_ev"],
                transition["formation_energy_ev"],
                marker="o",
                markersize=3.3,
                markerfacecolor="white",
                markeredgecolor=color,
                markeredgewidth=0.9,
                linestyle="none",
                zorder=3,
            )

    ax.axhline(0.0, color="#555555", linewidth=0.65, linestyle=(0, (3, 3)), zorder=0)
    ax.set_xlim(0.0, gap)
    ax.set_ylim(*y_bounds)
    ax.set_xlabel(r"Fermi level, $E_F$ (eV)")
    ax.set_ylabel("Formation energy (eV)")
    ax.tick_params(direction="out", length=3, width=0.7)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    for spine in ("left", "bottom"):
        ax.spines[spine].set_linewidth(0.75)
    if title:
        ax.set_title(title, fontsize=10)
    if len(defects) <= 6:
        ax.legend(
            frameon=False,
            loc="upper left",
            bbox_to_anchor=(1.01, 1.0),
            borderaxespad=0.0,
            fontsize=8,
        )
    else:
        ax.legend(
            frameon=False,
            loc="upper left",
            bbox_to_anchor=(1.01, 1.0),
            borderaxespad=0.0,
            fontsize=7,
            ncol=2,
        )
    output_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(
        output_path,
        format=output_path.suffix.lstrip("."),
        bbox_inches="tight",
    )
    plt.close(fig)
    return y_bounds


def write_outputs(
    audit: dict[str, Any],
    output_prefix: str | Path,
    *,
    max_defects_per_figure: int = 6,
    separate_defect_ids: Sequence[str] = (),
    y_max_ev: float | None = None,
    title: str | None = None,
    input_path: str | Path | None = None,
) -> dict[str, Any]:
    if (
        isinstance(max_defects_per_figure, bool)
        or not isinstance(max_defects_per_figure, int)
        or max_defects_per_figure < 1
    ):
        _fail("max_defects_per_figure must be a positive integer")
    prefix = Path(output_prefix).expanduser().resolve()
    defects = audit["defects"]
    if not defects:
        _fail("no defects available to plot")
    if isinstance(separate_defect_ids, (str, bytes)) or not isinstance(
        separate_defect_ids, Sequence
    ):
        _fail("separate_defect_ids must be a sequence of defect IDs")
    requested_separate = list(separate_defect_ids)
    if len(requested_separate) != len(set(requested_separate)):
        _fail("separate_defect_ids must not contain duplicates")
    known_ids = {defect["defect_id"] for defect in defects}
    unknown_ids = set(requested_separate) - known_ids
    if unknown_ids:
        _fail(f"unknown separate defect IDs: {', '.join(sorted(unknown_ids))}")

    figure_specs: list[tuple[list[dict[str, Any]], str, str | None]] = []
    if not requested_separate:
        chunks = [
            defects[index : index + max_defects_per_figure]
            for index in range(0, len(defects), max_defects_per_figure)
        ]
        figure_specs = [
            (chunk, "all", None if len(chunks) == 1 else f"{index:02d}")
            for index, chunk in enumerate(chunks, start=1)
        ]
    else:
        remaining = [
            defect for defect in defects if defect["defect_id"] not in requested_separate
        ]
        remaining_chunks = [
            remaining[index : index + max_defects_per_figure]
            for index in range(0, len(remaining), max_defects_per_figure)
        ]
        figure_specs.extend(
            (chunk, "remaining", f"remaining_{index:02d}")
            for index, chunk in enumerate(remaining_chunks, start=1)
        )
        for index, defect_id in enumerate(requested_separate, start=1):
            defect = next(item for item in defects if item["defect_id"] == defect_id)
            figure_specs.append(([defect], "separate", f"separate_{index:02d}"))
    stems = [
        prefix
        if suffix is None
        else prefix.with_name(f"{prefix.name}_{suffix}")
        for _, _, suffix in figure_specs
    ]
    audit_path = Path(f"{prefix}.audit.json")
    planned_paths = [
        audit_path,
        *(Path(f"{stem}.{extension}") for stem in stems for extension in ("svg", "png")),
    ]
    existing = [str(path) for path in planned_paths if path.exists()]
    if existing:
        _fail(f"refusing to overwrite existing outputs: {', '.join(existing)}")
    prefix.parent.mkdir(parents=True, exist_ok=True)

    figures = []
    for stem, (chunk, group, _) in zip(stems, figure_specs, strict=True):
        paths = {
            "svg": f"{stem}.svg",
            "png": f"{stem}.png",
        }
        y_bounds = None
        for extension in ("svg", "png"):
            y_bounds = _plot_defects(
                audit,
                chunk,
                Path(paths[extension]),
                y_max_ev=y_max_ev,
                title=title,
            )
        figures.append(
            {
                "defect_ids": [defect["defect_id"] for defect in chunk],
                "group": group,
                "files": paths,
                "y_limits_ev": list(y_bounds or ()),
                "x_limits_ev": [0.0, audit["band_edges"]["gap_ev"]],
                "high_energy_clipped": bool(
                    y_max_ev is not None
                    and any(
                        segment[key] > y_bounds[1]
                        for defect in chunk
                        for segment in defect["segments"]
                        for key in (
                            "formation_energy_start_ev",
                            "formation_energy_end_ev",
                        )
                    )
                ),
            }
        )

    audit["output"] = {
        "figures": figures,
        "max_defects_per_figure": max_defects_per_figure,
        "separate_defect_ids": requested_separate,
        "y_max_ev": y_max_ev,
        "input_path": str(Path(input_path).expanduser().resolve())
        if input_path is not None
        else None,
        "generated_at": datetime.now(timezone(timedelta(hours=8))).isoformat(),
    }
    audit["output"]["audit_json"] = str(audit_path)
    audit_path.write_text(
        json.dumps(audit, indent=2, ensure_ascii=False, allow_nan=False) + "\n",
        encoding="utf-8",
    )
    return audit


def main(argv: Sequence[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        description="Plot defect formation-energy lower envelopes from audited JSON."
    )
    parser.add_argument("--input", required=True, help="audited defect-state JSON")
    parser.add_argument(
        "--output-prefix",
        required=True,
        help="output path prefix for SVG, PNG, and audit JSON",
    )
    parser.add_argument(
        "--max-defects-per-figure",
        type=int,
        default=6,
        help="split crowded diagrams after this many defects (default: 6)",
    )
    parser.add_argument(
        "--separate-defect-id",
        action="append",
        default=[],
        help="plot this explicitly named defect in its own figure (repeatable)",
    )
    parser.add_argument(
        "--y-max-ev",
        type=float,
        help="optional upper formation-energy limit; low-energy envelopes remain visible",
    )
    parser.add_argument("--title", help="optional figure title")
    args = parser.parse_args(argv)
    try:
        raw = load_json(args.input)
        audit = analyze(raw)
        result = write_outputs(
            audit,
            args.output_prefix,
            max_defects_per_figure=args.max_defects_per_figure,
            separate_defect_ids=args.separate_defect_id,
            y_max_ev=args.y_max_ev,
            title=args.title,
            input_path=args.input,
        )
    except (OSError, ValueError) as exc:
        print(f"错误: {exc}", file=sys.stderr)
        return 2

    print(f"已写入审计 JSON: {result['output']['audit_json']}")
    for figure in result["output"]["figures"]:
        print(f"已写入图: {figure['files']['svg']} / {figure['files']['png']}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
