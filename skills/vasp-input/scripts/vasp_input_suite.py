#!/usr/bin/env python3
"""vibedft VASP input generation.

Usage:
    python3 vasp_input_suite.py --phase <phase> --config <config.json>
    python3 vasp_input_suite.py --phase auto --config vasp_tasks_config.json

Phases:
    relax     - Structure relaxation (ISIF=2/3), static, SCF
    band      - Band structure (NSCF, seekpath k-path)
    ccd       - Configuration coordinate diagram (structure interpolation)
    auto      - Auto-dispatch by task.type for relax/static/scf/band/ccd
"""

import argparse
import json
import sys
from datetime import datetime
from pathlib import Path

_SCRIPT_DIR = Path(__file__).resolve().parent
if str(_SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(_SCRIPT_DIR))


def _load_phase(phase_name: str):
    import importlib
    try:
        return importlib.import_module(f"phases.{phase_name}")
    except ImportError as exc:
        print(f"[error] cannot import phases.{phase_name}: {exc}", file=sys.stderr)
        sys.exit(1)


PHASE_MAP = {
    "relax": "relax",
    "band": "band",
    "ccd": "ccd",
}
# Types that map to phase name (task.type → phase)
TYPE_TO_PHASE = {
    "relax": "relax",
    "static": "relax",
    "ccd": "ccd",
    "band": "band",
    "scf": "relax",  # diagnostic SCF / DOS
}


def _auto(config: dict, workdir: Path) -> dict:
    """Auto-dispatch: read config['tasks'], route each by type."""
    tasks = config.get("tasks", [])
    if not tasks:
        print("[auto] No tasks defined")
        return {"steps": []}

    server_info = {
        "name": config.get("runtime", config.get("server", "")),
        "ncpus": int(config.get("ncpus", 1)),
        "partition": config.get("partition", ""),
        "qos": config.get("qos", ""),
    }

    # Build task lookup: id → task (for depends_on resolution)
    task_lookup = {}
    for t in tasks:
        t["__step_name"] = f"{t['id']}_{t['name']}"
        task_lookup[t["id"]] = t
    # Inject lookup into each task for child phases
    for t in tasks:
        t["__all_tasks"] = task_lookup
        t["host_structures"] = config.get("host_structures", {})

    all_steps = []
    for task in tasks:
        ttype = task.get("type", "")
        phase_name = TYPE_TO_PHASE.get(ttype)
        if phase_name is None:
            print(f"  [auto] SKIP {task.get('name')}: unknown type='{ttype}'")
            continue
        tid = task["id"]
        task_dir = workdir / "tasks" / f"{tid}_{task['name']}"
        task_dir.mkdir(parents=True, exist_ok=True)

        phase = _load_phase(phase_name)
        print(f"  [auto] {task['name']}: type={ttype} → {phase_name}")
        result = phase.generate_task(task, task_dir, server_info)
        steps = result.get("steps", [])
        all_steps.extend(steps)
        print(f"         → {len(steps)} steps")

    # Write merged tasks_manifest.json for vasp-submit-suite
    manifest = {
        "workdir": str(workdir),
        "runtime": server_info["name"],
        "ncpus": server_info["ncpus"],
        "steps": all_steps,
    }
    (workdir / "tasks_manifest.json").write_text(json.dumps(manifest, indent=2, ensure_ascii=False))
    print(f"\n  → tasks_manifest.json: {len(all_steps)} total steps")
    return manifest


def _write_task_manifest(workdir: Path, phase_name: str, config_name: str,
                         runtime: str, steps: list):
    """Write a per-phase task manifest to .task_manifests/.

    These manifests are consumed by vasp-submit-suite's --auto-submit mode
    to build a dependency DAG and manage job submission.
    """
    manifest_dir = workdir / ".task_manifests"
    manifest_dir.mkdir(exist_ok=True)
    manifest_path = manifest_dir / f"{phase_name}_{config_name}.json"

    phase_manifest = {
        "phase": phase_name,
        "config": config_name,
        "generated_at": datetime.now().isoformat(),
        "runtime": runtime,
        "tasks": [
            {
                "name": s["name"],
                "dir": s["dir"],
                "dependencies": s.get("depends_on", []),
                "type": phase_name,
            }
            for s in steps
        ],
    }
    manifest_path.write_text(json.dumps(phase_manifest, indent=2, ensure_ascii=False))
    print(f"  → manifest: {manifest_path.name} ({len(steps)} tasks)")


def main():
    parser = argparse.ArgumentParser(description="Unified VASP input generation suite")
    parser.add_argument(
        "--phase", required=True,
        choices=list(PHASE_MAP.keys()) + ["auto"],
        help="Calculation phase (auto = dispatch by task.type)"
    )
    parser.add_argument(
        "--config", required=True,
        help="Path to config JSON file"
    )
    args = parser.parse_args()

    config = json.loads(Path(args.config).read_text())
    workdir = Path(config.get("workdir", "."))
    workdir.mkdir(parents=True, exist_ok=True)

    if args.phase == "auto":
        result = _auto(config, workdir)
        _write_task_manifest(
            workdir, "auto", Path(args.config).stem,
            config.get("runtime", config.get("server", "")), result.get("steps", []),
        )
        return result

    phase_module = _load_phase(PHASE_MAP[args.phase])
    print(f"=== vasp-input-suite: phase={args.phase} ===")
    result = phase_module.generate(config, workdir)
    _write_task_manifest(
        workdir, args.phase, Path(args.config).stem,
        config.get("runtime", config.get("server", "")), result.get("steps", []),
    )
    print(f"=== Done ===")
    return result


if __name__ == "__main__":
    main()
