"""Generic runtime helpers for vibedft.

This repository does not embed any site-specific runtime defaults. Remote
execution details must come from user config JSON or private local environment
variables.
"""
from __future__ import annotations

import os


def get_vasp_cmd() -> str:
    """Return the VASP execution command template.

    The command must come from `VIBEDFT_VASP_CMD` and contain a `{vasp_bin}`
    placeholder so callers can substitute `vasp_std`, `vasp_gam`, etc.
    """
    cmd = os.environ.get("VIBEDFT_VASP_CMD", "").strip()
    if not cmd:
        raise ValueError(
            "VIBEDFT_VASP_CMD is not set. Export a site-specific command template "
            "such as a module-wrapper or scheduler launch command with a {vasp_bin} placeholder."
        )
    if "{vasp_bin}" not in cmd:
        raise ValueError("VIBEDFT_VASP_CMD must contain the placeholder {vasp_bin}.")
    return cmd


def default_ncore(ncpus: int | None, natoms: int | None = None) -> int:
    """Return a conservative generic NCORE guess.

    Public repo code avoids embedding site-specific node topology. Use one
    quarter of the requested MPI ranks when available, with a minimum of 1.
    """
    if not ncpus:
        return 1
    return max(int(ncpus) // 4, 1)


def get_python_env(server: str | None) -> str:
    return "python3"


def get_qos_count(qos: str) -> int:
    return 0


def count_my_jobs(qos: str, states: str = "RUNNING") -> int:
    return 0


def tres_exceeded(max_tres: dict) -> bool:
    return False


def partition_has_capacity(partition: str) -> bool:
    return True
