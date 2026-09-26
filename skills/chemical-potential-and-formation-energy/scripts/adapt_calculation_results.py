#!/usr/bin/env python3
"""Convert a READY calculation-results-1 export into audited plot inputs.

The calculation producer owns task completion and native provenance. This
adapter only consumes explicit JSON; it never scans VASP directories or
interprets a directory name as a successful calculation.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import re
import sys
from collections import Counter
from pathlib import Path
from typing import Any

from pymatgen.core import Composition


_IDENTITY_FIELDS = (
    "project_id", "run_id", "cluster", "task_id", "package_task_id",
    "step_id", "attempt_id", "physical_attempt_id",
    "selector_attempt", "job_id",
)
_SOURCE_BINDING_BASIS = "target_probe_selected_result_and_exact_identity"
_ALIGNMENT_BINDING_BASIS = "independent_alignment_audit_bound_to_native_sources"
_PENDING_REASON_CODES = {
    "CORRECTION_DECLARATION_MISSING",
    "E_CORR_MISSING",
}
_SHA256_RE = re.compile(r"^[0-9a-f]{64}$")


def _fail(message: str) -> None:
    raise ValueError(message)


def _mapping(value: Any, label: str) -> dict[str, Any]:
    if not isinstance(value, dict):
        _fail(f"{label} must be an object")
    return value


def _text(value: Any, label: str) -> str:
    if not isinstance(value, str) or not value.strip():
        _fail(f"{label} must be a non-empty string")
    return value.strip()


def _number(value: Any, label: str) -> float:
    if isinstance(value, bool) or not isinstance(value, (int, float)):
        _fail(f"{label} must be a finite number")
    result = float(value)
    if not math.isfinite(result):
        _fail(f"{label} must be a finite number")
    return result


def _positive_int(value: Any, label: str) -> int:
    if isinstance(value, bool) or not isinstance(value, int) or value < 1:
        _fail(f"{label} must be a positive integer")
    return value


def _integer(value: Any, label: str) -> int:
    if isinstance(value, bool) or not isinstance(value, int):
        _fail(f"{label} must be an integer")
    return value


def _unique_object(pairs: list[tuple[str, Any]]) -> dict[str, Any]:
    result: dict[str, Any] = {}
    for key, value in pairs:
        if key in result:
            _fail(f"duplicate JSON object key: {key}")
        result[key] = value
    return result


def load_json_file(path: str | Path) -> Any:
    return json.loads(
        Path(path).expanduser().read_text(encoding="utf-8"),
        object_pairs_hook=_unique_object,
    )


def _sha256(value: Any, label: str) -> str:
    text = _text(value, label)
    if not _SHA256_RE.fullmatch(text):
        _fail(f"{label} must be a lowercase SHA-256 digest")
    return text


def _pseudopotential_element_map(value: Any, label: str) -> dict[str, tuple[str, str]]:
    raw = _mapping(value, label)
    if not raw:
        _fail(f"{label} must contain trusted element-level POTCAR identities")
    result = {}
    for element, entry_value in raw.items():
        symbol = _text(element, f"{label} element")
        if symbol in result:
            _fail(f"{label} contains duplicate normalized element {symbol}")
        entry = _mapping(entry_value, f"{label}.{symbol}")
        result[symbol] = (
            _text(entry.get("title"), f"{label}.{symbol}.title"),
            _sha256(entry.get("sha256"), f"{label}.{symbol}.sha256"),
        )
    return result


def _validate_hse_methods(grouped: dict[str, list[tuple[int, dict[str, Any]]]]) -> None:
    signatures: set[tuple[Any, ...]] = set()
    element_maps: dict[int, dict[str, tuple[str, str]]] = {}
    for kind in ("pristine_hse", "defect_state_hse"):
        for row_index, row in grouped.get(kind, []):
            data = _mapping(row.get("data"), f"records[{row_index}].data")
            method = _mapping(
                data.get("method"), f"records[{row_index}].data.method"
            )
            if not isinstance(method.get("soc"), bool):
                _fail(f"records[{row_index}] HSE06 method has no explicit SOC setting")
            _text(
                method.get("pseudopotential_set_id"),
                f"records[{row_index}].method.pseudopotential_set_id",
            )
            raw_policy = method.get("pseudopotential_policy_id")
            family_policy = (
                None
                if raw_policy is None
                else _text(
                    raw_policy,
                    f"records[{row_index}].method.pseudopotential_policy_id",
                )
            )
            element_maps[row_index] = _pseudopotential_element_map(
                method.get("pseudopotentials"),
                f"records[{row_index}].data.method.pseudopotentials",
            )
            signatures.add((
                _text(method.get("method_id"), f"records[{row_index}].method_id"),
                _text(
                    method.get("reference_set_id"),
                    f"records[{row_index}].reference_set_id",
                ),
                family_policy,
                method["soc"],
                _text(
                    method.get("comparison_id"),
                    f"records[{row_index}].comparison_id",
                ),
                _text(
                    method.get("supercell_id"),
                    f"records[{row_index}].supercell_id",
                ),
                json.dumps(method.get("kpoints"), sort_keys=True),
            ))
    if len(signatures) != 1:
        _fail("native index has inconsistent HSE06 methods")
    pristine_index = grouped["pristine_hse"][0][0]
    pristine_map = element_maps[pristine_index]
    for row_index, _ in grouped["defect_state_hse"]:
        defect_map = element_maps[row_index]
        shared = set(pristine_map) & set(defect_map)
        if not shared:
            _fail(
                f"records[{row_index}] HSE06 defect and pristine methods "
                "have no shared element-level POTCAR evidence"
            )
        for element in shared:
            if pristine_map[element] != defect_map[element]:
                _fail(
                    f"records[{row_index}] HSE06 shared-element POTCAR "
                    f"identity differs for {element}"
                )


def _identity_values(identity: Any, label: str) -> tuple[str, ...]:
    source = _mapping(identity, label)
    return tuple(_text(source.get(field), f"{label}.{field}") for field in _IDENTITY_FIELDS)


def _bound_source(
    value: Any,
    identity: Any,
    label: str,
    *,
    require_sha256: bool = False,
) -> dict[str, Any]:
    source = _mapping(value, label)
    ref = _text(source.get("ref"), f"{label}.ref")
    _text(source.get("locator"), f"{label}.locator")
    if require_sha256:
        _sha256(ref, f"{label}.ref")
    if source.get("binding_basis") != _SOURCE_BINDING_BASIS:
        _fail(f"{label} does not declare exact selected-result identity binding")
    expected = _identity_values(identity, f"{label}.expected_identity")
    binding = _mapping(source.get("identity_binding"), f"{label}.identity_binding")
    actual = tuple(
        _text(binding.get(field), f"{label}.identity_binding.{field}")
        for field in _IDENTITY_FIELDS
    )
    if actual != expected:
        _fail(f"{label} is not bound to the exact task and attempt")
    return source


def _record_keys(values: Any, label: str) -> set[tuple[str, str]]:
    if not isinstance(values, list):
        _fail(f"{label} must be a list")
    keys: list[tuple[str, str]] = []
    for index, value in enumerate(values):
        item = _mapping(value, f"{label}[{index}]")
        keys.append((
            _text(item.get("kind"), f"{label}[{index}].kind"),
            _text(item.get("record_id"), f"{label}[{index}].record_id"),
        ))
    if len(keys) != len(set(keys)):
        _fail(f"{label} contains duplicate records")
    return set(keys)


def _pending_reason_code(value: Any, label: str) -> str:
    issue = _mapping(value, label)
    code = _text(issue.get("code"), f"{label}.code")
    if code not in _PENDING_REASON_CODES:
        _fail(f"{label} is not a declared correction-missing reason")
    return code


def _formula(value: Any, label: str) -> dict[str, float]:
    text = _text(value, label)
    try:
        formula = Composition(text).as_dict()
    except (ValueError, TypeError) as exc:
        raise ValueError(f"{label} is not a valid chemical formula: {text}") from exc
    if not formula or any(amount <= 0 for amount in formula.values()):
        _fail(f"{label} must have positive element counts")
    return {element: float(count) for element, count in formula.items()}


def _source(native_path: str, row_index: int, field: str, task_id: str) -> dict[str, str]:
    return {
        "task_id": task_id,
        "file": native_path,
        "locator": f"records[{row_index}].{field}",
    }


def _task(native_path: str, index: int, row: dict[str, Any]) -> dict[str, str]:
    identity = _mapping(row.get("identity"), f"records[{index}].identity")
    completion = _mapping(
        _mapping(row.get("provenance"), f"records[{index}].provenance").get("completion"),
        f"records[{index}].provenance.completion",
    )
    task = {
        key: _text(identity.get(key), f"records[{index}].identity.{key}")
        for key in ("project_id", "run_id", "cluster", "task_id", "step_id", "attempt_id")
    }
    task.update(
        {
            "job_id": _text(identity.get("job_id"), f"records[{index}].identity.job_id"),
            "result_ref": _text(
                identity.get("result_directory"),
                f"records[{index}].identity.result_directory",
            ),
            "status": "success",
            "acceptance_status": "accepted",
            "accepted_at": _text(
                completion.get("verified_at"),
                f"records[{index}].provenance.completion.verified_at",
            ),
            "acceptance_ref": _source(
                native_path, index, "producer_gate.status", task["task_id"]
            ),
        }
    )
    return task


def _normal_correction(
    index: int,
    row: dict[str, Any],
    resolutions: dict[str, Any],
    *,
    divisor: int,
    unit: str,
    phase_unit: str,
    pending_index: bool,
    strict_resolution_binding: bool,
) -> tuple[float, dict[str, Any]]:
    data = row["data"]
    key = f"{row['kind']}/{row['record_id']}"
    raw_native = data.get("correction")
    if raw_native is None:
        native = {}
        native_status = "missing"
    else:
        native = _mapping(raw_native, f"records[{index}].data.correction")
        native_status = native.get("status") or "missing"
    if native_status in {"not_applied", "missing"}:
        if native_status == "missing" and not pending_index:
            _fail(f"{key} has no native correction declaration")
        if native_status == "not_applied" and (
            native.get("applied") is not False or native.get("value_ev") is not None
        ):
            _fail(f"{key} not_applied correction contains an applied value")
        if key not in resolutions:
            _fail(f"{key} has no numerical PBE correction; provide an audited resolution")
        resolved = _mapping(resolutions.get(key), f"correction_resolutions.rows.{key}")
        if resolved.get("native_status") != native_status:
            _fail(f"{key} resolution must identify its native {native_status} status")
        decision_status = resolved.get("decision_status")
        if strict_resolution_binding and decision_status != "APPROVED":
            _fail(f"{key} PBE correction resolution is not APPROVED")
        if decision_status is not None and decision_status != "APPROVED":
            _fail(f"{key} PBE correction resolution is not APPROVED")
        if strict_resolution_binding:
            identity = _mapping(row.get("identity"), f"records[{index}].identity")
            if (
                _identity_values(
                    resolved.get("identity_binding"),
                    f"correction_resolutions.rows.{key}.identity_binding",
                )
                != _identity_values(identity, f"records[{index}].identity")
            ):
                _fail(f"{key} PBE correction resolution has a mismatched task or attempt")
            native_source = _bound_source(
                _mapping(
                    row.get("numerical_audit"),
                    f"records[{index}].numerical_audit",
                ).get("energy_source"),
                identity,
                f"records[{index}].numerical_audit.energy_source",
                require_sha256=True,
            )
            resolved_source = _bound_source(
                resolved.get("energy_source"),
                identity,
                f"correction_resolutions.rows.{key}.energy_source",
                require_sha256=True,
            )
            if (
                resolved_source.get("ref") != native_source.get("ref")
                or resolved_source.get("locator") != native_source.get("locator")
            ):
                _fail(f"{key} PBE correction resolution has a mismatched energy source")
        value = _number(resolved.get("value_ev"), f"{key}.value_ev")
        correction_unit = resolved.get("unit")
        source = _text(resolved.get("source"), f"{key}.source")
        _text(resolved.get("rationale"), f"{key}.rationale")
    elif native_status in {"provided", "zero_justified"}:
        if key in resolutions:
            _fail(f"{key} already has a numerical correction; refusing a second resolution")
        value = _number(native.get("value_ev"), f"{key}.correction.value_ev")
        correction_unit = native.get("unit")
        source = _text(native.get("source"), f"{key}.correction.source")
        if native_status == "zero_justified" and value != 0:
            _fail(f"{key} justified-zero correction is not zero")
    else:
        _fail(f"{key} PBE correction has no accepted status")
    if correction_unit not in {"eV/cell", "eV/formula_unit"}:
        _fail(f"{key} PBE correction has an unsupported unit")
    amount = value / divisor if correction_unit == "eV/cell" else value * unit
    return amount, {
        "value": amount,
        "unit": phase_unit,
        "source": source,
        "applied": True,
    }


def _file_document(
    path: str | Path,
    expected: Any,
    label: str,
) -> tuple[Path, bytes, str]:
    resolved = Path(path).expanduser().resolve()
    raw = resolved.read_bytes()
    parsed = json.loads(raw.decode("utf-8"), object_pairs_hook=_unique_object)
    if parsed != expected:
        _fail(f"{label} object does not match its source JSON file")
    return resolved, raw, hashlib.sha256(raw).hexdigest()


def _decision_identity_key(
    identity: Any,
    state_id: Any,
    charge: Any,
    label: str,
) -> tuple[tuple[str, ...], str, int]:
    return (
        _identity_values(identity, f"{label}.identity_binding"),
        _text(state_id, f"{label}.state_id"),
        _integer(charge, f"{label}.q"),
    )


def _defect_decisions(
    value: Any,
    *,
    native_sha256: str,
) -> dict[tuple[tuple[str, ...], str, int], tuple[int, dict[str, Any]]]:
    document = _mapping(value, "defect correction resolutions")
    if document.get("schema_version") != 1:
        _fail("defect correction resolutions require schema_version=1")
    if document.get("native_index_sha256") != native_sha256:
        _fail("defect correction resolutions are not bound to this native index SHA-256")
    decisions = document.get("decisions")
    if not isinstance(decisions, list) or not decisions:
        _fail("defect correction resolutions.decisions must be a non-empty list")
    result = {}
    for index, value in enumerate(decisions):
        decision = _mapping(value, f"defect correction resolutions.decisions[{index}]")
        key = _decision_identity_key(
            decision.get("identity_binding"),
            decision.get("state_id"),
            decision.get("q"),
            f"decisions[{index}]",
        )
        if key in result:
            _fail(f"duplicate defect correction decision for {key[1]} q={key[2]}")
        result[key] = (index, decision)
    return result


def _same_source(value: Any, expected: dict[str, Any], label: str) -> None:
    source = _mapping(value, label)
    for field in ("ref", "locator"):
        if source.get(field) != expected.get(field):
            _fail(f"{label}.{field} does not match the native source")


def _lattice_matrix(value: Any, label: str) -> list[list[float]]:
    if not isinstance(value, (list, tuple)) or len(value) != 3:
        _fail(f"{label} must be a 3x3 lattice")
    matrix = []
    for row_index, row in enumerate(value):
        if not isinstance(row, (list, tuple)) or len(row) != 3:
            _fail(f"{label}[{row_index}] must have three values")
        matrix.append([
            _number(item, f"{label}[{row_index}][{column_index}]")
            for column_index, item in enumerate(row)
        ])
    return matrix


def _potential_evidence(
    value: Any,
    *,
    identity: dict[str, Any],
    native_data: dict[str, Any],
    label: str,
) -> tuple[dict[str, Any], list[list[float]], str]:
    evidence = _mapping(value, label)
    source_type = _text(evidence.get("source_type"), f"{label}.source_type")
    if source_type not in {"LOCPOT", "site_potential", "site_potentials"}:
        _fail(f"{label}.source_type must identify raw potential evidence")
    source = _bound_source(
        evidence.get("source_ref"),
        identity,
        f"{label}.source_ref",
        require_sha256=True,
    )
    if source_type == "LOCPOT":
        native_potential = _mapping(native_data.get("locpot"), f"{label}.native_locpot")
    else:
        native_sources = native_data.get("potential_sources")
        if not isinstance(native_sources, dict) or source_type not in native_sources:
            _fail(f"{label} has no native source for {source_type}")
        native_potential = _mapping(
            native_sources[source_type], f"{label}.native_potential_sources.{source_type}"
        )
    if (
        native_potential.get("status") != "available"
        or native_potential.get("binding_basis")
        != "selected_attempt_artifact_and_exact_identity"
        or _identity_values(
            native_potential.get("identity_binding"),
            f"{label}.native_identity_binding",
        ) != _identity_values(identity, f"{label}.expected_identity")
        or _sha256(native_potential.get("sha256"), f"{label}.native_sha256") != source["ref"]
        or _text(native_potential.get("path"), f"{label}.native_path") != source["locator"]
        or isinstance(native_potential.get("size_bytes"), bool)
        or not isinstance(native_potential.get("size_bytes"), int)
        or native_potential["size_bytes"] <= 0
    ):
        _fail(f"{label} does not match the native selected potential artifact")
    structure = _mapping(evidence.get("structure"), f"{label}.structure")
    native_method = _mapping(native_data.get("method"), f"{label}.native_method")
    expected_supercell = _text(
        native_data.get("supercell_id") or native_method.get("supercell_id"),
        f"{label}.expected_supercell_id",
    )
    if _text(structure.get("supercell_id"), f"{label}.structure.supercell_id") != expected_supercell:
        _fail(f"{label} potential evidence belongs to a different supercell")
    matrix = _lattice_matrix(
        structure.get("lattice_angstrom"), f"{label}.structure.lattice_angstrom"
    )
    native_matrix = _lattice_matrix(
        native_data.get("lattice_angstrom"), f"{label}.native_lattice_angstrom"
    )
    if any(
        not math.isclose(actual, expected, rel_tol=1e-7, abs_tol=1e-5)
        for actual_row, expected_row in zip(matrix, native_matrix)
        for actual, expected in zip(actual_row, expected_row)
    ):
        _fail(f"{label} potential evidence lattice differs from the native structure")
    parameters = _mapping(evidence.get("parameters"), f"{label}.parameters")
    if not parameters or any(
        not isinstance(key, str) or not key.strip() or parameter is None
        or (isinstance(parameter, str) and not parameter.strip())
        for key, parameter in parameters.items()
    ):
        _fail(f"{label}.parameters must explicitly record non-empty source parameters")
    return source, matrix, expected_supercell


def _validate_alignment(
    value: Any,
    *,
    charge: int,
    correction_model: str,
    pristine_identity: dict[str, Any],
    defect_identity: dict[str, Any],
    pristine_data: dict[str, Any],
    defect_data: dict[str, Any],
    pristine_source: dict[str, Any],
    defect_source: dict[str, Any],
    native_energy_zero_id: Any,
    label: str,
) -> tuple[float, dict[str, Any] | None]:
    if charge == 0:
        if value is not None:
            _fail(f"{label} must omit charge alignment for q=0")
        return 0.0, None
    alignment = _mapping(value, label)
    if alignment.get("status") != "APPROVED":
        _fail(f"{label} is not APPROVED")
    treatment = alignment.get("treatment")
    if treatment not in {"included_in_model", "separate_same_reference"}:
        _fail(f"{label}.treatment must explicitly cover potential alignment")
    model = _text(alignment.get("model"), f"{label}.model")
    _text(alignment.get("source"), f"{label}.source")
    _text(alignment.get("rationale"), f"{label}.rationale")
    reference_id = _text(alignment.get("reference_id"), f"{label}.reference_id")
    if native_energy_zero_id and reference_id != native_energy_zero_id:
        _fail(f"{label}.reference_id differs from the native pristine energy zero")
    if treatment == "included_in_model" and model != correction_model:
        _fail(f"{label}.model must match the E_corr model when alignment is included")

    evidence = _mapping(alignment.get("evidence"), f"{label}.evidence")
    _same_source(
        evidence.get("pristine_energy_source"),
        pristine_source,
        f"{label}.evidence.pristine_energy_source",
    )
    _same_source(
        evidence.get("defect_energy_source"),
        defect_source,
        f"{label}.evidence.defect_energy_source",
    )
    pristine_potential_source, pristine_lattice, pristine_supercell = _potential_evidence(
        evidence.get("pristine_potential"),
        identity=pristine_identity,
        native_data=pristine_data,
        label=f"{label}.evidence.pristine_potential",
    )
    defect_potential_source, defect_lattice, defect_supercell = _potential_evidence(
        evidence.get("defect_potential"),
        identity=defect_identity,
        native_data=defect_data,
        label=f"{label}.evidence.defect_potential",
    )
    if pristine_supercell != defect_supercell or any(
        not math.isclose(actual, expected, rel_tol=1e-7, abs_tol=1e-5)
        for actual_row, expected_row in zip(pristine_lattice, defect_lattice)
        for actual, expected in zip(actual_row, expected_row)
    ):
        _fail(f"{label} raw potential sources are not comparable in one supercell")
    audit_source = _mapping(
        evidence.get("audit_source_ref"),
        f"{label}.evidence.audit_source_ref",
    )
    _sha256(audit_source.get("sha256"), f"{label}.evidence.audit_source_ref.sha256")
    _text(audit_source.get("locator"), f"{label}.evidence.audit_source_ref.locator")
    if audit_source.get("binding_basis") != _ALIGNMENT_BINDING_BASIS:
        _fail(f"{label} audit source is not bound to the exact native sources")
    potential_refs = _mapping(
        audit_source.get("potential_source_sha256"),
        f"{label}.evidence.audit_source_ref.potential_source_sha256",
    )
    if (
        potential_refs.get("pristine") != pristine_potential_source["ref"]
        or potential_refs.get("defect") != defect_potential_source["ref"]
    ):
        _fail(f"{label} audit source is not bound to both raw potential sources")
    bindings = _mapping(
        audit_source.get("identity_bindings"),
        f"{label}.evidence.audit_source_ref.identity_bindings",
    )
    if (
        _identity_values(bindings.get("pristine"), f"{label}.pristine_identity_binding")
        != _identity_values(pristine_identity, f"{label}.expected_pristine_identity")
        or _identity_values(bindings.get("defect"), f"{label}.defect_identity_binding")
        != _identity_values(defect_identity, f"{label}.expected_defect_identity")
    ):
        _fail(f"{label} audit source has a mismatched task or attempt")

    if treatment == "included_in_model":
        if "value_ev" in alignment or "unit" in alignment:
            _fail(f"{label} must not add a second alignment value already included in E_corr")
        return 0.0, alignment
    if alignment.get("unit") != "eV":
        _fail(f"{label}.unit must be eV for separately audited alignment")
    value_ev = _number(alignment.get("value_ev"), f"{label}.value_ev")
    return value_ev, alignment


def adapt(
    native: dict[str, Any],
    correction_resolutions: dict[str, Any],
    *,
    native_path: str | Path,
    defect_correction_resolutions: dict[str, Any] | None = None,
    defect_correction_path: str | Path | None = None,
    point_selection: dict[str, Any] | None = None,
    external_fixed_mu: dict[str, Any] | None = None,
    dopant_binary_phase_id: str | None = None,
) -> tuple[dict[str, Any], dict[str, Any]]:
    """Return phase input and the calculation index for the bundle builder."""
    native = _mapping(native, "native calculation results")
    if native.get("schema_version") != "calculation-results-1":
        _fail("native schema_version must be calculation-results-1")
    pending_index = native.get("status") == "INCOMPLETE"
    if native.get("status") == "READY":
        if native.get("ready") is not True:
            _fail("READY native calculation results must have ready=true")
    elif pending_index:
        if native.get("ready") is not False:
            _fail("INCOMPLETE native calculation results must have ready=false")
    else:
        _fail("native calculation results must be READY or correction-only INCOMPLETE")

    audit = _mapping(native.get("audit"), "native.audit")
    top_reasons = audit.get("rejection_reasons")
    if not isinstance(top_reasons, list):
        _fail("native.audit.rejection_reasons must be a list")
    if pending_index:
        if native.get("adapter_diagnostics", []) != []:
            _fail("INCOMPLETE native index has adapter_diagnostics")
        top_codes = [_pending_reason_code(item, f"native.audit.rejection_reasons[{i}]")
                     for i, item in enumerate(top_reasons)]
        if not top_codes:
            _fail(
                "native calculation results must be READY unless an INCOMPLETE "
                "index has correction-only rejection reasons"
            )
    elif top_reasons != []:
        _fail("native calculation results contain unresolved rejection reasons")

    rows = native.get("records")
    if not isinstance(rows, list) or not rows:
        _fail("native.records must be a non-empty list")
    expected_keys = _record_keys(native.get("expected_records"), "native.expected_records")
    actual_keys = _record_keys(rows, "native.records")
    if expected_keys != actual_keys:
        _fail("native expected_records do not match the actual record set")

    metadata = _mapping(native.get("metadata"), "native.metadata")
    material_id = _text(metadata.get("material_id"), "native.metadata.material_id")
    pbe_reference = _text(
        metadata.get("pbe_reference_set_id"), "native.metadata.pbe_reference_set_id"
    )
    _text(metadata.get("hse_reference_set_id"), "native.metadata.hse_reference_set_id")
    host_formula = _formula(metadata.get("formula"), "native.metadata.formula")
    resolution_doc = _mapping(correction_resolutions, "correction resolutions")
    if resolution_doc.get("schema_version") != 1:
        _fail("correction resolutions require schema_version=1")
    resolutions = _mapping(resolution_doc.get("rows"), "correction_resolutions.rows")
    native_file = Path(native_path).expanduser().resolve()
    path = str(native_file)
    native_sha256 = None
    if pending_index or defect_correction_resolutions is not None:
        _, _, native_sha256 = _file_document(
            native_file, native, "native calculation results"
        )
    strict_resolution_binding = pending_index or defect_correction_resolutions is not None
    if strict_resolution_binding:
        if resolution_doc.get("native_index_sha256") != native_sha256:
            _fail("PBE correction resolutions are not bound to this native index SHA-256")

    defect_doc = defect_correction_resolutions
    defect_file: Path | None = None
    defect_sha256: str | None = None
    if defect_correction_path is not None:
        if defect_doc is None:
            defect_doc = load_json_file(defect_correction_path)
        defect_file, _, defect_sha256 = _file_document(
            defect_correction_path, defect_doc, "defect correction resolutions"
        )
    elif defect_doc is not None:
        _fail("defect correction resolutions require their source file path")
    defect_decisions = (
        _defect_decisions(defect_doc, native_sha256=native_sha256)
        if defect_doc is not None and native_sha256 is not None
        else {}
    )

    grouped: dict[str, list[tuple[int, dict[str, Any]]]] = {}
    row_reason_codes: Counter[str] = Counter()
    for index, value in enumerate(rows):
        row = _mapping(value, f"records[{index}]")
        producer_gate = _mapping(row.get("producer_gate"), f"records[{index}].producer_gate")
        numerical_audit = _mapping(
            row.get("numerical_audit"), f"records[{index}].numerical_audit"
        )
        if producer_gate.get("status") != "VERIFIED":
            _fail(f"records[{index}] producer_gate is not VERIFIED")
        if numerical_audit.get("status") != "TRUSTED":
            _fail(f"records[{index}] numerical_audit is not TRUSTED")
        row_status = row.get("status")
        row_reasons = row.get("rejection_reasons")
        if pending_index:
            if row_status not in {"READY", "INCOMPLETE"}:
                _fail(f"records[{index}] status is not READY or correction-only INCOMPLETE")
            if producer_gate.get("rejection_reasons", []) != []:
                _fail(f"records[{index}] VERIFIED producer_gate has rejection reasons")
            if not isinstance(row_reasons, list):
                _fail(f"records[{index}].rejection_reasons must be a list")
            if row_status == "READY" and row_reasons:
                _fail(f"records[{index}] READY row has rejection reasons")
            if row_status == "INCOMPLETE" and not row_reasons:
                _fail(f"records[{index}] INCOMPLETE row has no declared correction gap")
            kind = _text(row.get("kind"), f"records[{index}].kind")
            data_for_reasons = _mapping(row.get("data"), f"records[{index}].data")
            for reason_index, reason in enumerate(row_reasons):
                code = _pending_reason_code(
                    reason, f"records[{index}].rejection_reasons[{reason_index}]"
                )
                row_reason_codes[code] += 1
                if code == "CORRECTION_DECLARATION_MISSING":
                    correction = data_for_reasons.get("correction")
                    correction_status = (
                        correction.get("status")
                        if isinstance(correction, dict)
                        else "missing"
                    )
                    if (
                        kind not in {"host_pbe", "competitive_phase_pbe", "elemental_pbe"}
                        or correction_status not in {"missing", "not_applied", None}
                        or (
                            correction_status == "not_applied"
                            and (
                                correction.get("applied") is not False
                                or correction.get("value_ev") is not None
                            )
                        )
                    ):
                        _fail(
                            f"records[{index}] CORRECTION_DECLARATION_MISSING "
                            "does not identify a missing PBE correction"
                        )
                elif code == "E_CORR_MISSING":
                    correction = data_for_reasons.get("E_corr")
                    correction_status = (
                        correction.get("status")
                        if isinstance(correction, dict)
                        else "missing"
                    )
                    if (
                        kind != "defect_state_hse"
                        or correction_status not in {"missing", "not_applied", None}
                        or (
                            correction_status == "not_applied"
                            and (
                                correction.get("applied") is not False
                                or correction.get("value_ev") is not None
                            )
                        )
                    ):
                        _fail(
                            f"records[{index}] E_CORR_MISSING does not identify "
                            "a missing defect correction"
                        )
        else:
            if row_status != "READY":
                _fail(f"records[{index}] status is not READY")
            if row_reasons != []:
                _fail(f"records[{index}] contains unresolved rejection reasons")
        convergence = _mapping(
            _mapping(row.get("data"), f"records[{index}].data").get("convergence"),
            f"records[{index}].data.convergence",
        )
        if convergence.get("electronic") is not True or convergence.get("ionic") is not True:
            _fail(f"records[{index}] lacks verified electronic and ionic convergence")
        convergence_source = _mapping(
            convergence.get("source"),
            f"records[{index}].data.convergence.source",
        )
        identity = _mapping(row.get("identity"), f"records[{index}].identity")
        source_binding = _mapping(
            convergence_source.get("identity_binding"),
            f"records[{index}].data.convergence.source.identity_binding",
        )
        if (
            not convergence_source.get("ref")
            or not convergence_source.get("locator")
            or convergence_source.get("binding_basis") != _SOURCE_BINDING_BASIS
            or any(
                not identity.get(field) or source_binding.get(field) != identity.get(field)
                for field in _IDENTITY_FIELDS
            )
        ):
            _fail(f"records[{index}] convergence source is not bound to its task attempt")
        grouped.setdefault(_text(row.get("kind"), f"records[{index}].kind"), []).append((index, row))

    if pending_index and set(row_reason_codes) != set(top_codes):
        _fail("top-level and row rejection reasons do not match correction-only gaps")

    for kind in ("host_pbe", "pristine_hse"):
        if len(grouped.get(kind, [])) != 1:
            _fail(f"native index must contain exactly one {kind} record")
    if not grouped.get("defect_state_hse"):
        _fail("native index contains no defect_state_hse records")
    if not grouped.get("elemental_pbe") or not grouped.get("elemental_hse"):
        _fail("native index must contain PBE and HSE06 elemental references")
    _validate_hse_methods(grouped)

    pristine_identity: dict[str, Any] | None = None
    pristine_band_source: dict[str, Any] | None = None
    pristine_energy_zero_id: Any = None
    if pending_index or defect_doc is not None:
        pristine_index, pristine_row = grouped["pristine_hse"][0]
        pristine_identity = _mapping(
            pristine_row.get("identity"), f"records[{pristine_index}].identity"
        )
        pristine_data = _mapping(
            pristine_row.get("data"), f"records[{pristine_index}].data"
        )
        band_evidence = _mapping(
            pristine_data.get("band_edges"),
            f"records[{pristine_index}].data.band_edges",
        )
        pristine_band_source = _bound_source(
            band_evidence.get("source"),
            pristine_identity,
            f"records[{pristine_index}].data.band_edges.source",
            require_sha256=True,
        )
        zero_provenance = _mapping(
            band_evidence.get("energy_zero_provenance"),
            f"records[{pristine_index}].data.band_edges.energy_zero_provenance",
        )
        if zero_provenance.get("status") not in {"VERIFIED", "DERIVED"}:
            _fail("INCOMPLETE native index lacks verified pristine band-edge zero provenance")
        if zero_provenance.get("status") == "DERIVED":
            if (
                _sha256(
                    zero_provenance.get("eigenval_sha256"),
                    "pristine band-edge eigenval_sha256",
                )
                != pristine_band_source["ref"]
            ):
                _fail("derived pristine band-edge zero is not bound to its EIGENVAL source")
        pristine_energy_zero_id = band_evidence.get("energy_zero_id")
        method_zero_id = _mapping(
            pristine_data.get("method"), f"records[{pristine_index}].data.method"
        ).get("energy_zero_id")
        if (
            pristine_energy_zero_id
            and method_zero_id
            and pristine_energy_zero_id != method_zero_id
        ):
            _fail("pristine band edges and method have different energy-zero IDs")
        for row_index, row in grouped["defect_state_hse"]:
            _bound_source(
                _mapping(
                    row.get("numerical_audit"),
                    f"records[{row_index}].numerical_audit",
                ).get("energy_source"),
                _mapping(row.get("identity"), f"records[{row_index}].identity"),
                f"records[{row_index}].numerical_audit.energy_source",
                require_sha256=True,
            )

    native_policy_ids = []
    for kind in ("host_pbe", "competitive_phase_pbe", "elemental_pbe"):
        for row_index, row in grouped.get(kind, []):
            method_data = _mapping(
                _mapping(row.get("data"), f"records[{row_index}].data").get("method"),
                f"records[{row_index}].data.method",
            )
            native_policy_id = method_data.get("correction_policy_id")
            if native_policy_id is not None:
                native_policy_ids.append(
                    _text(
                        native_policy_id,
                        f"records[{row_index}] PBE correction policy",
                    )
                )
    if len(set(native_policy_ids)) > 1:
        _fail("PBE rows have different correction policies")
    resolution_policy = _text(
        resolution_doc.get("policy_id"),
        "correction_resolutions.policy_id",
    )
    if native_policy_ids and resolution_policy != native_policy_ids[0]:
        _fail("correction resolution policy does not match the native PBE policy")
    if not native_policy_ids and not strict_resolution_binding:
        _fail(
            "missing native PBE policies require SHA-256- and row-bound "
            "external correction resolutions"
        )
    if "policy_review" in resolution_doc:
        review = _mapping(
            resolution_doc["policy_review"],
            "correction_resolutions.policy_review",
        )
        if review.get("status") != "APPROVED":
            _fail(
                "correction_resolutions.policy_review is "
                f"{review.get('status')!r}, not APPROVED"
            )
        _text(review.get("source"), "correction_resolutions.policy_review.source")
        _text(
            review.get("rationale"),
            "correction_resolutions.policy_review.rationale",
        )
        known = review.get("known_native_policy_ids")
        if (
            not isinstance(known, list)
            or any(not isinstance(item, str) or not item.strip() for item in known)
            or len(set(known)) != len(known)
            or set(known) != set(native_policy_ids)
        ):
            _fail(
                "correction_resolutions.policy_review known native policy IDs "
                "do not match the source index"
            )
    policy = resolution_policy

    host_index, host = grouped["host_pbe"][0]
    if _formula(host["data"].get("formula"), "host_pbe.formula") != host_formula:
        _fail("host PBE formula does not match the package formula")
    pbe_rows = [
        (index, row)
        for kind in ("host_pbe", "competitive_phase_pbe", "elemental_pbe")
        for index, row in grouped.get(kind, [])
    ]
    used_resolutions = set()

    def phase_energy(index: int, row: dict[str, Any], *, elemental: bool) -> dict[str, Any]:
        data = _mapping(row.get("data"), f"records[{index}].data")
        raw = _number(data.get("total_energy_ev"), f"records[{index}].data.total_energy_ev")
        formula_units = _positive_int(data.get("formula_units"), "formula_units")
        atoms = _positive_int(data.get("atom_count"), "atom_count") if elemental else None
        divisor = atoms if atoms is not None else formula_units
        unit = formula_units / atoms if atoms is not None else 1.0
        reported = _number(
            data.get("energy_ev_per_atom" if elemental else "energy_ev_per_formula_unit"),
            f"records[{index}].data.energy_per_unit",
        )
        if not math.isclose(raw / divisor, reported, abs_tol=1e-6, rel_tol=1e-7):
            _fail(f"records[{index}] normalized PBE energy disagrees with raw total")
        correction, provenance = _normal_correction(
            index,
            row,
            resolutions,
            divisor=divisor,
            unit=unit,
            phase_unit="eV/atom" if elemental else "eV",
            pending_index=pending_index,
            strict_resolution_binding=strict_resolution_binding,
        )
        key = f"{row['kind']}/{row['record_id']}"
        if key in resolutions:
            used_resolutions.add(key)
        return {
            "value": reported + correction,
            "unit": provenance["unit"],
            "functional": "PBE",
            "elemental_reference": pbe_reference,
            "correction": provenance,
        }

    host_energy = phase_energy(host_index, host, elemental=False)
    external_fixed_mu_data = (
        {}
        if external_fixed_mu is None
        else _mapping(external_fixed_mu, "external_fixed_mu")
    )
    phases = []
    phase_names = set()
    selected_binary = None
    external_phase_rows = []
    external_phase_elements = set()
    for index, row in grouped.get("competitive_phase_pbe", []):
        data = row["data"]
        # calculation-results-1 currently verifies the phase label in the
        # source snapshot but exports only its stable record_id.
        name = _text(
            data.get("phase") or row.get("record_id"),
            f"records[{index}].data.phase or record_id",
        )
        if name in phase_names:
            _fail(f"duplicate competing phase name: {name}")
        phase_names.add(name)
        formula = _formula(data.get("formula"), f"records[{index}].data.formula")
        phase_record = {
            "name": name,
            "formula": formula,
            "energy": phase_energy(index, row, elemental=False),
        }
        external_elements = set(formula) - set(host_formula)
        if external_elements:
            host_elements = set(formula) & set(host_formula)
            if len(formula) != 2 or len(external_elements) != 1 or len(host_elements) != 1:
                _fail(
                    f"external phase {row['record_id']} must be a binary with "
                    "one dopant and one host element"
                )
            dopant = next(iter(external_elements))
            external_phase_rows.append((index, row, phase_record, dopant))
            external_phase_elements.update(formula)
            if row["record_id"] == dopant_binary_phase_id:
                selected_binary = (index, row, phase_record)
        else:
            phases.append(phase_record)
    if dopant_binary_phase_id is not None and selected_binary is None:
        _fail(f"no competing phase has record_id={dopant_binary_phase_id!r}")
    selected_dopant = None
    if selected_binary is not None:
        selected_formula = selected_binary[2]["formula"]
        selected_dopant = next(iter(set(selected_formula) - set(host_formula)))
        if selected_dopant in external_fixed_mu_data:
            _fail(
                f"{selected_dopant} cannot use both an explicit fixed chemical "
                "potential and a binary-reservoir upper bound"
            )
    for _, row, _, dopant in external_phase_rows:
        if (
            row["record_id"] != dopant_binary_phase_id
            and dopant != selected_dopant
            and dopant not in external_fixed_mu_data
        ):
            _fail(
                f"external phase {row['record_id']} requires either "
                "dopant_binary_phase_id selection or an explicit external_fixed_mu"
            )
    elemental = {}
    for index, row in grouped["elemental_pbe"]:
        element = _text(row["data"].get("element"), f"records[{index}].data.element")
        if element in elemental:
            _fail(f"duplicate PBE elemental reference for {element}")
        elemental[element] = phase_energy(index, row, elemental=True)
    if set(resolutions) != used_resolutions:
        _fail("unused PBE correction resolutions: " + ", ".join(sorted(set(resolutions) - used_resolutions)))
    if not native_policy_ids:
        pbe_row_keys = {
            f"{row['kind']}/{row['record_id']}"
            for kind in ("host_pbe", "competitive_phase_pbe", "elemental_pbe")
            for _, row in grouped.get(kind, [])
        }
        if used_resolutions != pbe_row_keys:
            _fail(
                "missing native PBE policies require an approved bound "
                "external correction resolution for every PBE row"
            )
    needed_pbe = set(host_formula)
    for phase in phases:
        needed_pbe.update(phase["formula"])
    needed_pbe.update(external_phase_elements)
    needed_pbe.update(external_fixed_mu_data)
    reservoir = None
    if selected_binary is not None:
        binary_index, binary_row, binary = selected_binary
        formula = binary["formula"]
        extrinsic = set(formula) - set(host_formula)
        ligand = set(formula) & set(host_formula)
        if len(formula) != 2 or len(extrinsic) != 1 or len(ligand) != 1:
            _fail("selected dopant reservoir must be one dopant plus one host element")
        dopant = next(iter(extrinsic))
        needed_pbe.update(formula)
        formation_energy = binary["energy"]["value"] - sum(
            formula[element] * elemental[element]["value"]
            for element in formula if element in elemental
        )
        if any(element not in elemental for element in formula):
            _fail(f"selected binary reservoir lacks PBE elemental reference for {dopant}")
        reservoir = {
            "dopant": dopant,
            "phase_name": binary["name"],
            "formula": formula,
            "formation_energy_eV_per_formula": formation_energy,
            "elemental_reference": pbe_reference,
            "source": f"{path}#records[{binary_index}].data",
        }
    if needed_pbe != set(elemental):
        _fail("PBE elemental references must cover exactly the host and competitor species")

    hse_index, pristine = grouped["pristine_hse"][0]
    hse_data = _mapping(pristine.get("data"), "pristine_hse.data")
    hse_method = _mapping(hse_data.get("method"), "pristine_hse.data.method")
    method = "HSE06|" + "|".join(
        _text(hse_method.get(field), f"pristine_hse.method.{field}")
        for field in ("method_id", "reference_set_id", "pseudopotential_set_id")
    ) + f"|soc={hse_method.get('soc')}"
    band = _mapping(hse_data.get("band_edges"), "pristine_hse.data.band_edges")
    if band.get("converged") is not True:
        _fail("pristine HSE band edges must be converged")
    hse_task = _task(path, hse_index, pristine)
    band_refs = {
        key: _source(path, hse_index, f"data.band_edges.{value}", hse_task["task_id"])
        for key, value in (
            ("vbm_reference_ev", "vbm_ev"),
            ("cbm_reference_ev", "cbm_ev"),
            ("gap_ev", "gap_ev"),
            ("calculation_type", "converged"),
            ("converged", "converged"),
            ("kpoints", "method.kpoints"),
        )
    }
    band_input = {
        "calculation_method": method,
        "calculation_type": "self_consistent",
        "kpoints": {"scheme": "Gamma", "mesh": [1, 1, 1]},
        "converged": True,
        "vbm_reference_ev": _number(band.get("vbm_ev"), "pristine.band_edges.vbm_ev"),
        "cbm_reference_ev": _number(band.get("cbm_ev"), "pristine.band_edges.cbm_ev"),
        "gap_ev": _number(band.get("gap_ev"), "pristine.band_edges.gap_ev"),
        "source_refs": band_refs,
    }
    if pending_index:
        band_input["energy_zero_id"] = pristine_energy_zero_id
        band_input["energy_zero_provenance"] = dict(
            _mapping(
                hse_data["band_edges"].get("energy_zero_provenance"),
                "pristine.band_edges.energy_zero_provenance",
            )
        )
        band_input["energy_zero_source_ref"] = _source(
            path, hse_index, "data.band_edges.source", hse_task["task_id"]
        )
    index = {
        "schema_version": 1,
        "status": "accepted",
        "source_index_status": native["status"],
        "source_index_ready": native["ready"],
        "acceptance": {
            "status": "accepted",
            "accepted_at": hse_task["accepted_at"],
            "source_ref": {"file": path, "locator": "status"},
        },
        "material_id": material_id,
        "calculation_method": method,
        "pristine": {
            "task": hse_task,
            "E_host_ev": _number(hse_data.get("total_energy_ev"), "pristine_hse.total_energy_ev"),
            "E_host_source_ref": _source(path, hse_index, "data.total_energy_ev", hse_task["task_id"]),
            "supercell_lattice_angstrom": hse_data.get("lattice_angstrom"),
            "lattice_source_ref": _source(path, hse_index, "data.lattice_angstrom", hse_task["task_id"]),
            "band_edges": band_input,
        },
        "states": [],
        "elemental_references": {},
    }
    if native_sha256 is not None:
        index["source_index_sha256"] = native_sha256

    alignment_reference_ids: set[str] = set()
    alignment_decision_refs: list[dict[str, Any]] = []
    used_defect_decisions = set()
    for row_index, row in grouped["defect_state_hse"]:
        data = _mapping(row.get("data"), f"records[{row_index}].data")
        task = _task(path, row_index, row)
        identity = _mapping(row.get("identity"), f"records[{row_index}].identity")
        state_id = _text(data.get("state_id"), f"records[{row_index}].data.state_id")
        charge = _integer(data.get("q"), f"records[{row_index}].data.q")
        raw_correction = data.get("E_corr")
        correction = (
            _mapping(raw_correction, f"records[{row_index}].data.E_corr")
            if raw_correction is not None
            else {}
        )
        native_correction_status = correction.get("status") or "missing"
        decision_key = _decision_identity_key(
            identity, state_id, charge, f"records[{row_index}]"
        )
        decision_entry = defect_decisions.get(decision_key)
        if native_correction_status in {"not_applied", "missing"}:
            if not pending_index:
                _fail(f"defect {state_id} lacks an accepted numerical E_corr")
            if decision_entry is None:
                _fail(f"defect {state_id} q={charge} has no approved E_corr decision")
            decision_index, decision = decision_entry
            used_defect_decisions.add(decision_key)
            if decision.get("decision_status") != "APPROVED":
                _fail(f"defect {state_id} q={charge} E_corr decision is not APPROVED")
            if decision.get("native_status") != native_correction_status:
                _fail(f"defect {state_id} E_corr decision has the wrong native status")
            decision_correction_status = decision.get("correction_status")
            if decision_correction_status not in {
                "applied", "zero_justified", "justified_zero"
            }:
                _fail(f"defect {state_id} E_corr decision status is not approved for use")
            correction_value = _number(
                decision.get("value_ev"), f"decisions[{decision_index}].value_ev"
            )
            if (
                decision_correction_status in {"zero_justified", "justified_zero"}
                and correction_value != 0
            ):
                _fail(f"defect {state_id} justified-zero E_corr decision is not zero")
            correction_status = (
                "justified_zero"
                if decision_correction_status in {"zero_justified", "justified_zero"}
                else "applied"
            )
            if decision.get("unit") != "eV":
                _fail(f"defect {state_id} E_corr decision unit must be eV")
            correction_model = _text(
                decision.get("model"), f"decisions[{decision_index}].model"
            )
            correction_source = _text(
                decision.get("source"), f"decisions[{decision_index}].source"
            )
            correction_rationale = _text(
                decision.get("rationale"), f"decisions[{decision_index}].rationale"
            )
            correction_file_ref = {
                "task_id": task["task_id"],
                "file": str(defect_file),
                "locator": f"decisions[{decision_index}]",
                "sha256": defect_sha256,
            }
            alignment_value, alignment_evidence = _validate_alignment(
                decision.get("alignment"),
                charge=charge,
                correction_model=correction_model,
                pristine_identity=pristine_identity or {},
                defect_identity=identity,
                pristine_data=hse_data,
                defect_data=data,
                pristine_source=pristine_band_source or {},
                defect_source=_mapping(
                    _mapping(
                        row.get("numerical_audit"),
                        f"records[{row_index}].numerical_audit",
                    ).get("energy_source"),
                    f"records[{row_index}].numerical_audit.energy_source",
                ),
                native_energy_zero_id=pristine_energy_zero_id,
                label=f"decisions[{decision_index}].alignment",
            )
            correction_value += alignment_value
            correction_native_source_ref = _source(
                path, row_index, "data.E_corr", task["task_id"]
            )
        elif native_correction_status in {"provided", "zero_justified"}:
            if decision_entry is not None and (
                charge == 0 or decision_entry[1].get("purpose") != "alignment_only"
            ):
                _fail(f"defect {state_id} already has a native E_corr; refusing a second correction")
            correction_status = (
                "justified_zero"
                if native_correction_status == "zero_justified"
                else "applied"
            )
            correction_value = _number(
                correction.get("value_ev"), f"records[{row_index}].data.E_corr.value_ev"
            )
            correction_model = _text(
                correction.get("model"), f"records[{row_index}].data.E_corr.model"
            )
            correction_source = _text(
                correction.get("source"), f"records[{row_index}].data.E_corr.source"
            )
            native_rationale = correction.get("rationale")
            correction_rationale = (
                native_rationale.strip()
                if isinstance(native_rationale, str)
                else ""
            )
            if correction.get("unit") != "eV":
                _fail(f"defect {state_id} native E_corr unit must be eV")
            if correction_status == "justified_zero" and correction_value != 0:
                _fail(f"defect {state_id} justified-zero E_corr is not zero")
            correction_native_source_ref = _source(
                path, row_index, "data.E_corr", task["task_id"]
            )
            correction_file_ref = correction_native_source_ref
            alignment_evidence = None
            if decision_entry is None and pending_index:
                alignment_value, alignment_evidence = _validate_alignment(
                    correction.get("alignment"),
                    charge=charge,
                    correction_model=correction_model,
                    pristine_identity=pristine_identity or {},
                    defect_identity=identity,
                    pristine_data=hse_data,
                    defect_data=data,
                    pristine_source=pristine_band_source or {},
                    defect_source=_mapping(
                        _mapping(
                            row.get("numerical_audit"),
                            f"records[{row_index}].numerical_audit",
                        ).get("energy_source"),
                        f"records[{row_index}].numerical_audit.energy_source",
                    ),
                    native_energy_zero_id=pristine_energy_zero_id,
                    label=f"records[{row_index}].data.E_corr.alignment",
                )
                correction_value += alignment_value
                if correction_status == "justified_zero" and correction_value != 0:
                    correction_status = "applied"
            elif decision_entry is not None:
                decision_index, decision = decision_entry
                used_defect_decisions.add(decision_key)
                if (
                    decision.get("decision_status") != "APPROVED"
                    or decision.get("native_status") != native_correction_status
                    or decision.get("correction_status") != "native_value_unchanged"
                    or decision.get("unit") != "eV"
                    or _number(decision.get("value_ev"), f"decisions[{decision_index}].value_ev")
                    != correction_value
                    or decision.get("model") != correction_model
                ):
                    _fail(f"defect {state_id} alignment-only decision changes native E_corr")
                _text(decision.get("source"), f"decisions[{decision_index}].source")
                _text(decision.get("rationale"), f"decisions[{decision_index}].rationale")
                alignment_value, alignment_evidence = _validate_alignment(
                    decision.get("alignment"),
                    charge=charge,
                    correction_model=correction_model,
                    pristine_identity=pristine_identity or {},
                    defect_identity=identity,
                    pristine_data=hse_data,
                    defect_data=data,
                    pristine_source=pristine_band_source or {},
                    defect_source=_mapping(
                        _mapping(
                            row.get("numerical_audit"),
                            f"records[{row_index}].numerical_audit",
                        ).get("energy_source"),
                        f"records[{row_index}].numerical_audit.energy_source",
                    ),
                    native_energy_zero_id=pristine_energy_zero_id,
                    label=f"decisions[{decision_index}].alignment",
                )
                correction_value += alignment_value
                if correction_status == "justified_zero" and correction_value != 0:
                    correction_status = "applied"
                correction_file_ref = {
                    "task_id": task["task_id"],
                    "file": str(defect_file),
                    "locator": f"decisions[{decision_index}]",
                    "sha256": defect_sha256,
                }
        else:
            _fail(f"defect {state_id} has no accepted numerical E_corr status")

        if charge != 0:
            if alignment_evidence is None:
                _fail(f"charged defect {state_id} lacks audited alignment evidence")
            alignment_reference_ids.add(
                _text(alignment_evidence.get("reference_id"), "alignment.reference_id")
            )
            alignment_decision_refs.append({
                "state_id": state_id,
                "source_ref": correction_file_ref,
            })
        fields = {
            "state_id": "state_id",
            "defect_id": "defect_id",
            "geometry_id": "geometry_id",
            "charge": "q",
            "n_i": "n_i",
            "E_def_ev": "E_def_ev",
            "E_corr_ev": "E_corr.value_ev",
            "supercell_lattice_angstrom": "lattice_angstrom",
        }
        state_item = {
            "state_id": state_id,
            "defect_id": data.get("defect_id"),
            "geometry_id": data.get("geometry_id"),
            "charge": charge,
            "n_i": data.get("n_i"),
            "E_def_ev": data.get("E_def_ev"),
            "E_corr_ev": correction_value,
            "E_corr_status": correction_status,
            "E_corr_description": "; ".join(
                value for value in (
                    correction_model,
                    correction_source,
                    correction_rationale,
                    (
                        "potential alignment " + alignment_evidence["treatment"]
                        if alignment_evidence else None
                    ),
                    alignment_evidence.get("model") if alignment_evidence else None,
                    alignment_evidence.get("source") if alignment_evidence else None,
                    alignment_evidence.get("rationale") if alignment_evidence else None,
                ) if value
            ),
            "material_id": material_id,
            "calculation_method": method,
            "supercell_lattice_angstrom": data.get("lattice_angstrom"),
            "task": task,
            "source_refs": {
                key: _source(path, row_index, f"data.{value}", task["task_id"])
                for key, value in fields.items()
            },
            "E_corr_native_source_ref": correction_native_source_ref,
            "E_corr_decision_source_ref": correction_file_ref,
        }
        if alignment_evidence is not None:
            state_item["E_corr_alignment_evidence"] = alignment_evidence
        state_item["source_refs"]["E_corr_ev"] = correction_file_ref
        index["states"].append(state_item)

    if pending_index and len(alignment_reference_ids) > 1:
        _fail("charged defects do not share one audited band-edge alignment reference")
    if alignment_reference_ids:
        band_input["energy_zero_id"] = next(iter(alignment_reference_ids))
        band_input["alignment_decision_source_refs"] = alignment_decision_refs
        index["pristine"]["band_edges"] = band_input

    if defect_decisions and set(defect_decisions) != used_defect_decisions:
        unused = sorted(
            (key[1], key[2]) for key in set(defect_decisions) - used_defect_decisions
        )
        _fail(f"unused defect correction decisions: {unused}")
    for row_index, row in grouped["elemental_hse"]:
        data = row["data"]
        element = _text(data.get("element"), "elemental_hse.element")
        if element in index["elemental_references"]:
            _fail(f"duplicate HSE06 elemental reference for {element}")
        task = _task(path, row_index, row)
        index["elemental_references"][element] = {
            "energy_ev": _number(data.get("energy_ev_per_atom"), f"elemental_hse.{element}.energy_ev_per_atom"),
            "method": method,
            "unit": "eV/atom",
            "reference_id": _text(data["method"].get("reference_set_id"), "HSE elemental reference"),
            "task": task,
            "source_ref": _source(path, row_index, "data.energy_ev_per_atom", task["task_id"]),
        }
    phase_input = {
        "energy_basis": "total_energy",
        "correction_policy": {"scheme": policy, "mode": "already_applied"},
        "origin": {
            "native_calculation_results": path,
            "native_status": native["status"],
            **({"native_index_sha256": native_sha256} if native_sha256 else {}),
            "correction_resolutions": resolution_doc,
            **({
                "defect_correction_resolutions": {
                    "file": str(defect_file),
                    "sha256": defect_sha256,
                }
            } if defect_file is not None else {}),
            "normalization": (
                "raw VASP total energy per formula unit or atom plus one "
                "unit-normalized, explicitly audited correction"
            ),
        },
        "host": {"name": "host", "formula": host_formula, "energy": host_energy},
        "phases": phases,
        "elemental_energies": elemental,
        "point_selection": point_selection or {},
        "external_fixed_mu": external_fixed_mu_data,
        **({"dopant_reservoir": reservoir} if reservoir is not None else {}),
    }
    return phase_input, index


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        description="Adapt a READY or correction-only INCOMPLETE calculation-results-1 index."
    )
    parser.add_argument("--native-index", required=True)
    parser.add_argument("--correction-resolutions", required=True)
    parser.add_argument("--defect-correction-resolutions")
    parser.add_argument("--output-dir", required=True)
    parser.add_argument("--control-element")
    parser.add_argument("--external-fixed-mu", help="optional audited JSON mapping of external elements")
    parser.add_argument(
        "--dopant-binary-phase-id",
        help="native competitive_phase_pbe record_id used as one per-point dopant binary limit",
    )
    args = parser.parse_args(argv)
    try:
        native_path = Path(args.native_index).expanduser().resolve()
        corrections = load_json_file(args.correction_resolutions)
        defect_corrections = (
            load_json_file(args.defect_correction_resolutions)
            if args.defect_correction_resolutions else None
        )
        external = (
            load_json_file(args.external_fixed_mu)
            if args.external_fixed_mu else None
        )
        phase_input, calculation_index = adapt(
            load_json_file(native_path),
            corrections,
            native_path=native_path,
            defect_correction_resolutions=defect_corrections,
            defect_correction_path=args.defect_correction_resolutions,
            point_selection={"control_element": args.control_element}
            if args.control_element else None,
            external_fixed_mu=external,
            dopant_binary_phase_id=args.dopant_binary_phase_id,
        )
        out = Path(args.output_dir).expanduser().resolve()
        if out.exists():
            _fail(f"output directory must not already exist: {out}")
        out.mkdir(parents=True)
        (out / "phase-input.json").write_text(
            json.dumps(phase_input, indent=2, ensure_ascii=False, allow_nan=False) + "\n",
            encoding="utf-8",
        )
        (out / "formation-calculation-index.json").write_text(
            json.dumps(calculation_index, indent=2, ensure_ascii=False, allow_nan=False) + "\n",
            encoding="utf-8",
        )
    except (ValueError, OSError, TypeError) as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 2
    print(f"Written: {out / 'phase-input.json'}")
    print(f"Written: {out / 'formation-calculation-index.json'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
