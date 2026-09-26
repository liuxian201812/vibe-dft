"""Phase: relax/static/scf.

Basic vibedft input generation for normal VASP calculations. Supported task
types are structure relaxation, static single-point, and simple SCF/DOS-style
single-point inputs. Phonon, NEB, MD, and dependency chains are intentionally
not included.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

from pymatgen.core import Structure
from pymatgen.io.vasp import Kpoints, Poscar

_SCRIPT_DIR = Path(__file__).resolve().parent
_PARENT_SCRIPT = _SCRIPT_DIR.parent
if str(_PARENT_SCRIPT) not in sys.path:
    sys.path.insert(0, str(_PARENT_SCRIPT))

from common import get_auto_u_config, get_encut, make_incar, make_pu_modifier, write_potcar, write_submit_sh
from shared_utils.server import default_ncore


def _write_kpoints(task_dir: Path, gamma: bool, density: float = 0.04):
    if gamma:
        Kpoints.gamma_automatic([1, 1, 1]).write_file(str(task_dir / "KPOINTS"))
        return
    struct = Structure.from_file(str(task_dir / "POSCAR"))
    mesh = [max(1, round(1.0 / (a * density))) for a in struct.lattice.abc]
    Kpoints.monkhorst_automatic(mesh).write_file(str(task_dir / "KPOINTS"))


def _apply_substitutions(struct: Structure, substitutions: list[dict]) -> Structure:
    out = struct.copy()
    for sub in substitutions:
        from_el = sub.get("from", "")
        to_el = sub.get("to", "")
        if not from_el or not to_el:
            continue
        if sub.get("sites") is not None:
            for idx in sub["sites"]:
                out[int(idx)] = to_el
        elif sub.get("all", False):
            out.replace_species({from_el: to_el})
        else:
            count = int(sub.get("count", 0))
            replaced = 0
            for idx, site in enumerate(out):
                if site.specie.symbol == from_el and replaced < count:
                    out[idx] = to_el
                    replaced += 1
    return out


def _build_task(task: dict, task_dir: Path, server_info: dict) -> dict:
    ncpus = int(server_info.get("ncpus") or 1)
    partition = server_info.get("partition") or ""
    qos = server_info.get("qos", "")
    vasp_cmd = task.get("vasp_cmd", "")
    if not partition:
        raise ValueError("partition is required in config or task runtime info")

    host_structs = task.get("host_structures", {})
    host_key = task.get("host", "")
    host_path = task.get("poscar") or host_structs.get(host_key)
    if not host_path:
        raise ValueError(f"Task {task.get('name', '')}: missing poscar or host_structures[{host_key!r}]")

    task_dir.mkdir(parents=True, exist_ok=True)
    struct = Structure.from_file(host_path)
    struct = _apply_substitutions(struct, task.get("substitutions", []))
    Poscar(struct).write_file(str(task_dir / "POSCAR"))

    gamma = bool(task.get("kpoints_gamma", False))
    _write_kpoints(task_dir, gamma, float(task.get("kpoints_density", 0.04)))

    try:
        write_potcar(task_dir, server=server_info.get("name", ""))
    except Exception as exc:
        print(f"  POTCAR skip ({exc})")

    potcar = task_dir / "POTCAR"
    encut = int(task.get("encut") or (get_encut(potcar) if potcar.exists() else 400))
    template = task.get("template", "LR")
    if template not in {"SR", "LR", "ST"}:
        raise ValueError(f"Unsupported template '{template}'. Use SR, LR, or ST.")

    incar = make_incar(
        template_key=template,
        modifiers=task.get("modifiers", []),
        encut=encut,
        ncore=default_ncore(ncpus, natoms=len(struct)),
        nupdown=task.get("nupdown"),
        params=task.get("params", {}),
        isif=task.get("isif"),
        functional=task.get("functional", "PBE"),
        overrides=task.get("overrides", {}),
    )

    if task.get("lda_u", True):
        u_config = get_auto_u_config(task_dir / "POSCAR", task.get("u", {}))
        if u_config:
            pu = make_pu_modifier(u_config, str(task_dir / "POSCAR"))
            if pu:
                incar += "\n" + pu

    (task_dir / "INCAR").write_text(incar)

    is_relax = template in {"SR", "LR"}
    vasp_bin = task.get("vasp_bin") or ("vasp_gam" if gamma else "vasp_std")
    write_submit_sh(
        task_dir.parent,
        task_dir.name,
        vasp_bin=vasp_bin,
        ncores=ncpus,
        partition=partition,
        qos=qos,
        vasp_cmd=vasp_cmd,
        restart_loop=is_relax,
        max_loops=int(task.get("max_loops", 3)),
    )
    return {"name": task_dir.name, "dir": str(task_dir), "depends_on": []}


def generate_task(task: dict, task_dir: Path, server_info: dict | None = None) -> dict:
    step = _build_task(task, task_dir, server_info or {})
    return {"steps": [step], "bin": task.get("vasp_bin", "vasp_std")}


def generate(config: dict, workdir: Path) -> dict:
    server_info = {
        "name": config.get("runtime", config.get("server", "")),
        "ncpus": int(config.get("ncpus", 1)),
        "partition": config.get("partition", ""),
        "qos": config.get("qos", ""),
    }
    tasks_dir = workdir / "tasks"
    tasks_dir.mkdir(parents=True, exist_ok=True)
    steps = []
    for task in config.get("tasks", []):
        task = dict(task)
        task["host_structures"] = config.get("host_structures", {})
        name = f"{task['id']}_{task['name']}"
        steps.append(_build_task(task, tasks_dir / name, server_info))

    manifest = {"workdir": str(workdir), "steps": steps}
    (workdir / "tasks_manifest.json").write_text(json.dumps(manifest, indent=2))
    print(f"\nDone. {len(steps)} tasks in {workdir}")
    return manifest
