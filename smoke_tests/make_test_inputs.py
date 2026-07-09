#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path

from pymatgen.core import Lattice, Structure


def write_poscar(path: Path, a: float = 5.64):
    structure = Structure(Lattice.cubic(a), ["Na", "Cl"], [[0, 0, 0], [0.5, 0.5, 0.5]])
    structure.to(filename=str(path), fmt="poscar")


def main():
    root = Path(__file__).resolve().parent / "tmp"
    root.mkdir(parents=True, exist_ok=True)

    host = root / "POSCAR_host"
    gs = root / "POSCAR_gs"
    es = root / "POSCAR_es"
    write_poscar(host)
    write_poscar(gs, 5.64)
    write_poscar(es, 5.80)

    (root / "fake_potcars" / "Na_pv").mkdir(parents=True, exist_ok=True)
    (root / "fake_potcars" / "Cl").mkdir(parents=True, exist_ok=True)
    fake_na = "PAW_PBE Na_pv\n   ENMAX = 200.000; ENMIN = 150.000\n"
    fake_cl = "PAW_PBE Cl\n   ENMAX = 250.000; ENMIN = 180.000\n"
    (root / "fake_potcars" / "Na_pv" / "POTCAR").write_text(fake_na)
    (root / "fake_potcars" / "Cl" / "POTCAR").write_text(fake_cl)

    relax_config = {
        "workdir": str(root / "relax_case"),
        "runtime": "smoke-test",
        "partition": "debug",
        "ncpus": 4,
        "host_structures": {"host": str(host)},
        "tasks": [
            {
                "id": "01",
                "name": "relax",
                "type": "relax",
                "host": "host",
                "template": "LR",
                "functional": "PBE",
                "isif": 3,
                "kpoints_density": 0.20,
            }
        ],
    }
    (root / "relax_config.json").write_text(json.dumps(relax_config, indent=2))

    band_config = {
        "workdir": str(root / "band_case"),
        "runtime": "smoke-test",
        "partition": "debug",
        "ncpus": 4,
        "poscar": str(host),
        "functional": "PBE",
    }
    (root / "band_config.json").write_text(json.dumps(band_config, indent=2))

    ccd_config = {
        "workdir": str(root / "ccd_case"),
        "runtime": "smoke-test",
        "partition": "debug",
        "ncpus": 4,
        "poscar_gs": str(gs),
        "poscar_es": str(es),
        "functional": "PBE",
        "n_images": 5,
        "nupdown": 0,
    }
    (root / "ccd_config.json").write_text(json.dumps(ccd_config, indent=2))

    print(root)


if __name__ == "__main__":
    main()
