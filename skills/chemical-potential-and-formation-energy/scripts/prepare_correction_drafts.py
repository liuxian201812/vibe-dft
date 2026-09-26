#!/usr/bin/env python3
"""Prepare identity-bound, deliberately unapproved correction review drafts.

This reads one native JSON index. It neither reads VASP files nor calculates,
approves, or assigns a zero value to any missing scientific correction.
"""

from __future__ import annotations

import argparse
import copy
import hashlib
import importlib.util
import json
import sys
from pathlib import Path
from typing import Any

_ADAPTER_PATH = Path(__file__).with_name("adapt_calculation_results.py")
_ADAPTER_SPEC = importlib.util.spec_from_file_location("_correction_draft_adapter", _ADAPTER_PATH)
if _ADAPTER_SPEC is None or _ADAPTER_SPEC.loader is None:
    raise RuntimeError(f"cannot load {_ADAPTER_PATH}")
adapter = importlib.util.module_from_spec(_ADAPTER_SPEC)
_ADAPTER_SPEC.loader.exec_module(adapter)


PBE_KINDS = frozenset({"host_pbe", "competitive_phase_pbe", "elemental_pbe"})


def _fail(message: str) -> None:
    raise ValueError(message)


def _potential_hint(row: dict[str, Any]) -> dict[str, Any]:
    potential = row["data"].get("locpot")
    if not isinstance(potential, dict) or potential.get("status") != "available":
        return {"status": "missing_native_potential"}
    return {
        "status": "available",
        "path": potential.get("path"),
        "sha256": potential.get("sha256"),
        "size_bytes": potential.get("size_bytes"),
        "identity_binding": copy.deepcopy(potential.get("identity_binding")),
    }


def _review_drafts(native: dict[str, Any], native_sha: str) -> tuple[dict[str, Any], dict[str, Any]]:
    if native.get("schema_version") != "calculation-results-1":
        _fail("native index schema is not calculation-results-1")
    status = native.get("status")
    if status not in {"READY", "INCOMPLETE"} or native.get("ready") is not (status == "READY"):
        _fail("native index is neither READY nor correction-only INCOMPLETE")
    audit = adapter._mapping(native.get("audit"), "native.audit")
    top_reasons = audit.get("rejection_reasons")
    if not isinstance(top_reasons, list):
        _fail("native rejection reasons must be a list")
    if status == "READY" and top_reasons:
        _fail("READY native index has rejection reasons")
    if status == "INCOMPLETE":
        if native.get("adapter_diagnostics", []) != []:
            _fail("native index has unresolved adapter diagnostics")
        top_codes = {
            adapter._pending_reason_code(reason, "native.audit.rejection_reasons")
            for reason in top_reasons
        }
        if not top_codes:
            _fail("INCOMPLETE index is not limited to correction gaps")
    else:
        top_codes = set()

    rows = native.get("records")
    if not isinstance(rows, list) or not rows:
        _fail("native index has no records")
    if adapter._record_keys(native.get("expected_records"), "expected_records") != adapter._record_keys(
        rows, "records"
    ):
        _fail("native expected records do not match actual records")

    row_codes: set[str] = set()
    policies: set[str] = set()
    missing_policies = 0
    pbe_rows: dict[str, Any] = {}
    defect_decisions: list[dict[str, Any]] = []
    pristine = [
        row for row in rows if row.get("kind") == "pristine_hse"
    ]
    if len(pristine) != 1:
        _fail("native index needs one pristine HSE supercell")
    for index, row in enumerate(rows):
        row = adapter._mapping(row, f"records[{index}]")
        identity = adapter._mapping(row.get("identity"), f"records[{index}].identity")
        adapter._identity_values(identity, f"records[{index}].identity")
        gate = adapter._mapping(row.get("producer_gate"), f"records[{index}].producer_gate")
        numerical = adapter._mapping(row.get("numerical_audit"), f"records[{index}].numerical_audit")
        if gate.get("status") != "VERIFIED" or numerical.get("status") != "TRUSTED":
            _fail(f"records[{index}] lacks trusted production and numerical gates")
        data = adapter._mapping(row.get("data"), f"records[{index}].data")
        kind = adapter._text(row.get("kind"), f"records[{index}].kind")
        convergence = adapter._mapping(data.get("convergence"), f"records[{index}].convergence")
        if convergence.get("electronic") is not True or convergence.get("ionic") is not True:
            _fail(f"records[{index}] has unverified convergence")
        energy_source = adapter._bound_source(
            numerical.get("energy_source"), identity,
            f"records[{index}].energy_source",
            require_sha256=kind in PBE_KINDS or kind == "defect_state_hse",
        )
        adapter._bound_source(
            convergence.get("source"), identity,
            f"records[{index}].convergence.source",
        )
        reasons = row.get("rejection_reasons")
        if not isinstance(reasons, list):
            _fail(f"records[{index}] has no explicit rejection-reason list")
        codes = {
            adapter._pending_reason_code(reason, f"records[{index}].rejection_reasons")
            for reason in reasons
        }
        row_codes.update(codes)
        if row.get("status") == "READY" and reasons:
            _fail(f"records[{index}] is READY but has rejection reasons")
        if row.get("status") == "INCOMPLETE" and not reasons:
            _fail(f"records[{index}] is INCOMPLETE without correction gaps")
        if row.get("status") not in {"READY", "INCOMPLETE"}:
            _fail(f"records[{index}] has an unsupported state")
        if status == "READY" and row.get("status") != "READY":
            _fail(f"records[{index}] is not READY in a READY index")
        if "CORRECTION_DECLARATION_MISSING" in codes and kind not in PBE_KINDS:
            _fail(f"records[{index}] PBE correction reason does not identify a PBE row")
        if "E_CORR_MISSING" in codes and kind != "defect_state_hse":
            _fail(f"records[{index}] defect correction reason does not identify a defect")

        if kind in PBE_KINDS:
            method = adapter._mapping(data.get("method"), f"records[{index}].method")
            policy_id = method.get("correction_policy_id")
            if isinstance(policy_id, str) and policy_id.strip():
                policies.add(policy_id.strip())
            else:
                missing_policies += 1
            raw_corr = data.get("correction")
            corr = adapter._mapping(raw_corr, f"records[{index}].correction") if raw_corr is not None else {}
            corr_status = corr.get("status") or "missing"
            if corr_status in {"missing", "not_applied"}:
                if corr_status == "missing" and status == "READY":
                    _fail("READY PBE row lacks a native correction declaration")
                if corr_status == "missing" and "CORRECTION_DECLARATION_MISSING" not in codes:
                    _fail(f"records[{index}] lacks a declared PBE correction gap")
                key = f"{kind}/{adapter._text(row.get('record_id'), 'record_id')}"
                pbe_rows[key] = {
                    "native_status": corr_status,
                    "decision_status": "PENDING",
                    "identity_binding": {
                        field: identity[field] for field in adapter._IDENTITY_FIELDS
                    },
                    "energy_source": copy.deepcopy(energy_source),
                    "value_ev": None,
                    "unit": corr.get("unit") if corr.get("unit") in {"eV/cell", "eV/formula_unit"} else None,
                    "source": None,
                    "rationale": None,
                }
            elif corr_status not in {"provided", "zero_justified"}:
                _fail(f"records[{index}] has an unsupported PBE correction state")
            if "CORRECTION_DECLARATION_MISSING" in codes and corr_status not in {"missing", "not_applied"}:
                _fail(f"records[{index}] correction gap contradicts its native value")
        if kind == "defect_state_hse":
            correction = data.get("E_corr")
            corr = adapter._mapping(correction, f"records[{index}].E_corr") if correction is not None else {}
            corr_status = corr.get("status") or "missing"
            q = adapter._integer(data.get("q"), f"records[{index}].q")
            missing = corr_status in {"missing", "not_applied"}
            if missing and status == "READY":
                _fail("READY defect row has no native E_corr")
            if missing and "E_CORR_MISSING" not in codes:
                _fail(f"records[{index}] lacks a declared defect correction gap")
            if "E_CORR_MISSING" in codes and not missing:
                _fail(f"records[{index}] E_corr gap contradicts its native value")
            if missing or q != 0:
                decision = {
                    "identity_binding": {
                        field: identity[field] for field in adapter._IDENTITY_FIELDS
                    },
                    "state_id": adapter._text(data.get("state_id"), f"records[{index}].state_id"),
                    "q": q,
                    "native_status": corr_status,
                    "decision_status": "PENDING",
                    "correction_status": "unresolved" if missing else "native_value_unchanged",
                    "value_ev": None if missing else corr.get("value_ev"),
                    "unit": "eV",
                    "model": None if missing else corr.get("model"),
                    "source": None,
                    "rationale": None,
                    "alignment": None,
                }
                if not missing:
                    decision["purpose"] = "alignment_only"
                if q != 0:
                    decision["native_potential_hints"] = {
                        "pristine": _potential_hint(pristine[0]),
                        "defect": _potential_hint(row),
                    }
                defect_decisions.append(decision)

    if status == "INCOMPLETE" and row_codes != top_codes:
        _fail("top-level and row correction gaps differ")
    if len(policies) > 1:
        _fail("PBE rows do not share one native correction policy")
    if not policies and not missing_policies:
        _fail("native index contains no PBE rows")
    # A known ID on only some rows is evidence to review, not authority to
    # impute that policy to the rows whose native declaration is missing.
    policy = next(iter(policies)) if policies and not missing_policies else None
    return (
        {
            "schema_version": 1,
            "native_index_sha256": native_sha,
            "policy_id": policy,
            **({
                "policy_review": {
                    "status": "PENDING",
                    "source": None,
                    "rationale": None,
                    "known_native_policy_ids": sorted(policies),
                }
            } if policy is None else {}),
            "rows": pbe_rows,
        },
        {
            "schema_version": 1,
            "native_index_sha256": native_sha,
            "decisions": defect_decisions,
        },
    )


def prepare(native_path: str | Path, output_dir: str | Path) -> dict[str, Any]:
    source = Path(native_path).expanduser().resolve()
    output = Path(output_dir).expanduser().resolve()
    if output.exists():
        _fail(f"output directory already exists: {output}")
    raw = source.read_bytes()
    native = json.loads(raw.decode("utf-8"), object_pairs_hook=adapter._unique_object)
    pbe, defects = _review_drafts(native, hashlib.sha256(raw).hexdigest())
    summary = {
        "schema_version": 1,
        "status": "PENDING_SCIENTIFIC_REVIEW",
        "native_index": str(source),
        "native_status": native["status"],
        "pbe_decisions": len(pbe["rows"]),
        "defect_decisions": len(defects["decisions"]),
        "no_correction_values_or_approvals_inferred": True,
    }
    output.mkdir(parents=True)
    files = {
        "pbe-corrections-draft.json": pbe,
        "draft-manifest.json": summary,
    }
    if defects["decisions"]:
        files["defect-corrections-draft.json"] = defects
    for name, payload in files.items():
        (output / name).write_text(
            json.dumps(payload, indent=2, ensure_ascii=False, allow_nan=False) + "\n",
            encoding="utf-8",
        )
    return summary


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description="Write unapproved native-bound correction review drafts.")
    parser.add_argument("--native-index", required=True)
    parser.add_argument("--output-dir", required=True)
    args = parser.parse_args(argv)
    try:
        result = prepare(args.native_index, args.output_dir)
    except (OSError, ValueError, TypeError) as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 2
    print(json.dumps(result, ensure_ascii=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
