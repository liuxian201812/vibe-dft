"""Phase: ccd (configuration-coordinate single-point images).

Basic vibedft CCD support is intentionally direct: both endpoint structures
must already exist. This module interpolates between them and writes standalone
VASP single-point directories for each image. It does not manage relaxation
dependencies or deferred POSCAR generation.
"""
from __future__ import annotations

import json
import shutil
import sys
from pathlib import Path

import numpy as np
from pymatgen.core import Structure
from pymatgen.io.vasp import Kpoints, Poscar

_SCRIPT_DIR = Path(__file__).resolve().parent
_PARENT_SCRIPT = _SCRIPT_DIR.parent
if str(_PARENT_SCRIPT) not in sys.path:
    sys.path.insert(0, str(_PARENT_SCRIPT))

from common import get_encut, make_incar, replace_tag, write_potcar, write_submit_sh
from shared_utils.server import default_ncore


def _pbc_delta(frac1, frac2):
    delta = frac2 - frac1
    return delta - np.round(delta)


def _interpolate(s1: Structure, s2: Structure, alpha: float) -> Structure:
    if len(s1) != len(s2):
        raise ValueError("CCD endpoints must have the same number of sites")
    if [str(site.specie) for site in s1] != [str(site.specie) for site in s2]:
        raise ValueError("CCD endpoints must have the same site species order")
    out = s1.copy()
    for i in range(len(s1)):
        new_frac = (s1[i].frac_coords + alpha * _pbc_delta(s1[i].frac_coords, s2[i].frac_coords)) % 1.0
        out.replace(i, str(s1[i].specie), coords=new_frac)
    return out


def _parse_alphas(config: dict) -> list[float]:
    if "alphas" in config:
        return [float(x) for x in config["alphas"]]
    start = float(config.get("alpha_start", -0.5))
    stop = float(config.get("alpha_stop", 1.5))
    n_images = int(config.get("n_images", 41))
    if n_images < 2:
        raise ValueError("n_images must be >= 2")
    return [float(x) for x in np.linspace(start, stop, n_images)]


def _copy_or_make_potcar(ccd_dir: Path, config: dict, server: str):
    if config.get("potcar"):
        shutil.copy2(config["potcar"], ccd_dir / "POTCAR")
        return
    write_potcar(ccd_dir, server=server)


def generate(config: dict, workdir: Path) -> dict:
    """Generate CCD image directories.

    Required config keys:
      poscar_gs: ground-state endpoint POSCAR/CONTCAR
      poscar_es: excited-state endpoint POSCAR/CONTCAR

    Optional keys:
      workdir, functional, runtime, ncpus, partition, qos, nupdown, n_images,
      alpha_start, alpha_stop, alphas, encut, potcar, vasp_cmd
    """
    poscar_gs = Path(config["poscar_gs"])
    poscar_es = Path(config["poscar_es"])
    if not poscar_gs.exists():
        raise FileNotFoundError(f"poscar_gs not found: {poscar_gs}")
    if not poscar_es.exists():
        raise FileNotFoundError(f"poscar_es not found: {poscar_es}")

    runtime = config.get("runtime", config.get("server", ""))
    ncpus = int(config.get("ncpus", 1))
    partition = config.get("partition", "")
    qos = config.get("qos", "")
    functional = config.get("functional", "PBE")
    nupdown = config.get("nupdown")
    vasp_cmd = config.get("vasp_cmd", "")
    if not partition:
        raise ValueError("partition is required for CCD submit script generation")

    s_gs = Structure.from_file(str(poscar_gs))
    s_es = Structure.from_file(str(poscar_es))
    alphas = _parse_alphas(config)

    ccd_dir = workdir / config.get("output_subdir", "ccd")
    ccd_dir.mkdir(parents=True, exist_ok=True)
    Poscar(s_gs).write_file(str(ccd_dir / "POSCAR"))
    _copy_or_make_potcar(ccd_dir, config, runtime)

    encut = int(config.get("encut") or get_encut(ccd_dir / "POTCAR"))
    incar = make_incar(
        template_key="ST",
        functional=functional,
        encut=encut,
        ncore=default_ncore(ncpus, natoms=len(s_gs)),
        nupdown=nupdown,
        overrides={"NSW": 0, "IBRION": -1, "ISIF": 0, **config.get("overrides", {})},
    )
    if nupdown is not None:
        incar = replace_tag(incar, "NUPDOWN", str(int(nupdown)))
    (ccd_dir / "INCAR").write_text(incar)
    Kpoints.gamma_automatic([1, 1, 1]).write_file(str(ccd_dir / "KPOINTS"))

    steps = []
    images = []
    for idx, alpha in enumerate(alphas):
        img_name = f"image_{idx:03d}"
        img_dir = ccd_dir / img_name
        img_dir.mkdir(parents=True, exist_ok=True)
        Poscar(_interpolate(s_gs, s_es, alpha)).write_file(str(img_dir / "POSCAR"))
        for fname in ("INCAR", "KPOINTS", "POTCAR"):
            shutil.copy2(ccd_dir / fname, img_dir / fname)
        write_submit_sh(
            ccd_dir,
            img_name,
            vasp_bin="vasp_gam",
            ncores=ncpus,
            partition=partition,
            qos=qos,
            vasp_cmd=vasp_cmd,
            restart_loop=False,
        )
        entry = {"name": img_name, "dir": str(img_dir), "alpha": alpha}
        steps.append({"name": img_name, "dir": str(img_dir), "depends_on": []})
        images.append(entry)

    manifest = {
        "phase": "ccd",
        "poscar_gs": str(poscar_gs),
        "poscar_es": str(poscar_es),
        "functional": functional,
        "encut": encut,
        "n_images": len(images),
        "images": images,
        "steps": steps,
    }
    (ccd_dir / "ccd_manifest.json").write_text(json.dumps(manifest, indent=2))
    print(f"  ccd: OK ({len(images)} images, {functional}, encut={encut})")
    return manifest


def generate_task(task: dict, task_dir: Path, server_info: dict | None = None) -> dict:
    merged = dict(task)
    merged.setdefault("workdir", str(task_dir.parent))
    merged.setdefault("output_subdir", task_dir.name)
    if server_info:
        merged.setdefault("runtime", server_info.get("name") or "")
        merged.setdefault("ncpus", server_info.get("ncpus"))
        merged.setdefault("partition", server_info.get("partition"))
        merged.setdefault("qos", server_info.get("qos"))
    result = generate(merged, task_dir.parent)
    return {"steps": result.get("steps", [])}
