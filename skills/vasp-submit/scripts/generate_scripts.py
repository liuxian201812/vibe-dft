#!/usr/bin/env python3
"""Generate independent SLURM scripts for VASP tasks.

Accepted config formats:
1. Basic chain-style config:
   {"workdir": "...", "partition": "...", "ncpus": 32, "steps": [{"name": "relax", "inputs": {...}}]}

2. vasp-input task config:
   {"workdir": "...", "partition": "...", "ncpus": 32, "tasks": [...], "host_structures": {...}}

Dependency fields are intentionally ignored in vibedft. Each step is submitted
independently; no afterok chain and no CONTCAR handoff are generated.
"""
from __future__ import annotations

import json
import shutil
import sys
from pathlib import Path

_SKILLS_ROOT = Path(__file__).resolve().parents[3] / "skills"
if str(_SKILLS_ROOT) not in sys.path:
    sys.path.insert(0, str(_SKILLS_ROOT))

from shared_utils import get_vasp_cmd


def _server_info(config: dict) -> dict:
    partition = config.get("partition", "")
    ncpus = int(config.get("ncpus", 1))
    if not partition:
        raise ValueError("partition is required in submit config")
    return {
        "name": config.get("runtime", config.get("server", "")),
        "partition": partition,
        "qos": config.get("qos", ""),
        "ncpus": ncpus,
        "vasp_cmd": config.get("vasp_cmd") or get_vasp_cmd(),
    }


def _slurm_header(name: str, server: dict, time: str) -> str:
    qos_line = f"#SBATCH --qos={server['qos']}\n" if server.get("qos") else ""
    return (
        "#!/bin/bash\n"
        f"#SBATCH -J {name}\n"
        f"#SBATCH -p {server['partition']}\n"
        f"#SBATCH -n {server['ncpus']}\n"
        f"#SBATCH --time={time}\n"
        "#SBATCH -o %j.log\n"
        f"{qos_line}"
        "ulimit -s unlimited\n"
    )


def _vasp_tasks_to_steps(config: dict) -> list[dict]:
    tasks = config.get("tasks", [])
    host_structs = config.get("host_structures", {})
    workdir = Path(config["workdir"])
    steps = []
    for task in tasks:
        name = f"{task['id']}_{task['name']}"
        task_dir = workdir / "tasks" / name
        inputs = {
            "poscar": host_structs.get(task.get("host", ""), str(task_dir / "POSCAR")),
            "incar": str(task_dir / "INCAR"),
            "potcar": str(task_dir / "POTCAR"),
            "kpoints": str(task_dir / "KPOINTS"),
        }
        steps.append({
            "name": name,
            "inputs": inputs,
            "vasp_bin": task.get("vasp_bin", "vasp_gam" if task.get("kpoints_gamma") else "vasp_std"),
            "time": task.get("time", "36:00:00"),
        })
    return steps


def _copy_inputs(step: dict, step_dir: Path):
    step_dir.mkdir(parents=True, exist_ok=True)
    for key, dest_name in (("poscar", "POSCAR"), ("incar", "INCAR"), ("potcar", "POTCAR"), ("kpoints", "KPOINTS")):
        src = step.get("inputs", {}).get(key, "")
        if src and Path(src).exists():
            dst = step_dir / dest_name
            if Path(src).resolve() != dst.resolve():
                shutil.copy2(src, dst)


def _write_run_script(workdir: Path, step: dict, server: dict):
    name = step["name"]
    step_dir = workdir / name
    _copy_inputs(step, step_dir)
    vasp_bin = step.get("vasp_bin", "vasp_std")
    vasp_cmd = server["vasp_cmd"].replace("{vasp_bin}", vasp_bin)
    script = _slurm_header(name, server, step.get("time", "36:00:00"))
    script += f'cd "{step_dir}" || exit 1\n'
    script += 'if [ ! -s POSCAR ] || [ ! -s INCAR ] || [ ! -s KPOINTS ] || [ ! -s POTCAR ]; then\n'
    script += '  echo "ERROR: missing POSCAR/INCAR/KPOINTS/POTCAR" >&2\n  exit 1\nfi\n'
    for part in vasp_cmd.split("&&"):
        script += part.strip() + "\n"
    run_path = workdir / f"run_{name}.sh"
    run_path.write_text(script)
    run_path.chmod(0o755)
    print(f"  {run_path.name}")


def _write_submit(workdir: Path, steps: list[dict]):
    lines = ["#!/bin/bash", f'cd "{workdir}" || exit 1']
    for step in steps:
        name = step["name"]
        lines.append(f'sbatch "run_{name}.sh"')
    submit_path = workdir / "submit.sh"
    submit_path.write_text("\n".join(lines) + "\n")
    submit_path.chmod(0o755)
    print("  submit.sh")


def main():
    if len(sys.argv) < 2:
        print(__doc__)
        sys.exit(1)
    config = json.loads(Path(sys.argv[1]).read_text())
    workdir = Path(config["workdir"])
    workdir.mkdir(parents=True, exist_ok=True)
    server = _server_info(config)

    steps = config.get("steps") or _vasp_tasks_to_steps(config)
    if not steps:
        print("No steps or tasks defined", file=sys.stderr)
        sys.exit(1)

    for step in steps:
        _write_run_script(workdir, step, server)
    _write_submit(workdir, steps)
    print(f"\n  workdir: {workdir}")
    if server["name"]:
        print(f"  runtime: {server['name']}")
    print(f"  scheduler: partition={server['partition']} ncpus={server['ncpus']}")
    print(f"  steps:   {len(steps)} independent jobs")
    print(f"  submit:  cd {workdir} && bash submit.sh")


if __name__ == "__main__":
    main()
