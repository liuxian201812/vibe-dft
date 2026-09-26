#!/usr/bin/env python3
"""Run all three audited postprocessing stages from one native result index."""

from __future__ import annotations

import argparse
import importlib.util
import json
import sys
from pathlib import Path
from typing import Any


def _load(name: str):
    file = Path(__file__).with_name(f"{name}.py")
    spec = importlib.util.spec_from_file_location(name, file)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"cannot load {file}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


adapter = _load("adapt_calculation_results")
phase = _load("chemical_potential_phase_diagram")
bundle = _load("build_formation_energy_bundle")


def _read_json(path: str | Path) -> dict[str, Any]:
    return adapter.load_json_file(path)


def run(
    native_path: str | Path,
    correction_path: str | Path,
    output_dir: str | Path,
    *,
    defect_correction_path: str | Path | None = None,
    control_element: str | None = None,
    external_fixed_mu: dict[str, Any] | None = None,
    dopant_binary_phase_id: str | None = None,
    separate_defect_ids: tuple[str, ...] = (),
    all_in_one: bool = False,
    y_max_ev: float | None = None,
) -> dict[str, Any]:
    if all_in_one and separate_defect_ids:
        raise ValueError("--all-in-one conflicts with --separate-defect-id")
    native_path = Path(native_path).expanduser().resolve()
    correction_path = Path(correction_path).expanduser().resolve()
    output = Path(output_dir).expanduser().resolve()
    if output.exists():
        raise ValueError(f"output directory must not already exist: {output}")
    phase_input, formation_index = adapter.adapt(
        _read_json(native_path),
        _read_json(correction_path),
        native_path=native_path,
        defect_correction_resolutions=(
            _read_json(defect_correction_path)
            if defect_correction_path is not None else None
        ),
        defect_correction_path=defect_correction_path,
        point_selection={"control_element": control_element} if control_element else None,
        external_fixed_mu=external_fixed_mu,
        dopant_binary_phase_id=dopant_binary_phase_id,
    )
    analysis = phase.analyze(phase_input)
    if not analysis["selected_points"]:
        raise ValueError(
            "no selected points; set --control-element for this host, or provide "
            "manually audited full-vector chemical potentials"
        )
    phase_path = output / "phase" / "analysis.json"
    index_path = output / "formation-calculation-index.json"
    max_defects_per_figure = (
        max(1, len(formation_index["states"])) if all_in_one else 6
    )

    # Fail before creating output even if the phase alone is valid but formation
    # needs a missing elemental reference, correction, or usable plotting bound.
    bundle.preflight_bundle(
        formation_index,
        analysis,
        calculation_index_path=index_path,
        phase_analysis_path=phase_path,
        separate_defect_ids=separate_defect_ids,
        max_defects_per_figure=max_defects_per_figure,
        y_max_ev=y_max_ev,
    )
    # No VASP job, task directory, or remote state is touched.
    output.mkdir(parents=True)
    (output / "phase-input.json").write_text(
        json.dumps(phase_input, indent=2, ensure_ascii=False, allow_nan=False) + "\n",
        encoding="utf-8",
    )
    index_path.write_text(
        json.dumps(formation_index, indent=2, ensure_ascii=False, allow_nan=False) + "\n",
        encoding="utf-8",
    )
    phase.write_outputs(analysis, output / "phase")
    figures = bundle.build_bundle(
        formation_index,
        analysis,
        output / "formation",
        calculation_index_path=index_path,
        phase_analysis_path=phase_path,
        separate_defect_ids=separate_defect_ids,
        max_defects_per_figure=max_defects_per_figure,
        y_max_ev=y_max_ev,
    )
    return {
        "phase_analysis": str(phase_path),
        "formation_bundle": figures["manifest"],
        "points": [item["point_name"] for item in figures["points"]],
    }


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        description="Run phase and formation-energy figures from one audited VASP results index."
    )
    parser.add_argument("--native-index", required=True)
    parser.add_argument("--correction-resolutions", required=True)
    parser.add_argument("--defect-correction-resolutions")
    parser.add_argument("--output-dir", required=True)
    parser.add_argument("--control-element")
    parser.add_argument("--external-fixed-mu", help="JSON file with explicit external relative chemical potentials")
    parser.add_argument("--dopant-binary-phase-id", help="binary competing-phase record_id for point-dependent dopant limit")
    parser.add_argument("--separate-defect-id", action="append", default=[])
    parser.add_argument(
        "--all-in-one",
        action="store_true",
        help="draw all accepted defect IDs in one formation-energy figure per chemical-potential point",
    )
    parser.add_argument("--y-max-ev", type=float)
    args = parser.parse_args(argv)
    try:
        result = run(
            args.native_index,
            args.correction_resolutions,
            args.output_dir,
            defect_correction_path=args.defect_correction_resolutions,
            control_element=args.control_element,
            external_fixed_mu=_read_json(args.external_fixed_mu) if args.external_fixed_mu else None,
            dopant_binary_phase_id=args.dopant_binary_phase_id,
            separate_defect_ids=tuple(args.separate_defect_id),
            all_in_one=args.all_in_one,
            y_max_ev=args.y_max_ev,
        )
    except (OSError, ValueError, TypeError, phase.DomainError) as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 2
    print(f"Phase analysis: {result['phase_analysis']}")
    print(f"Formation bundle: {result['formation_bundle']}")
    print("Points: " + ", ".join(result["points"]))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
