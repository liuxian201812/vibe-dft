#!/usr/bin/env python3
"""Local-only public entry point for VASP input and submit-script generation."""
from __future__ import annotations

import argparse
from pathlib import Path
import subprocess
import sys


SCRIPT_DIR = Path(__file__).resolve().parent


def _run(script: str, args: list[str]) -> int:
    if args[:1] == ["--"]:
        args = args[1:]
    return subprocess.call([sys.executable, str(SCRIPT_DIR / script), *args])


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Generate local VASP inputs or independent SLURM submit scripts."
    )
    subparsers = parser.add_subparsers(dest="action", required=True)
    for action, help_text in (
        ("prepare", "generate VASP inputs"),
        ("submit-scripts", "generate independent submit scripts"),
    ):
        child = subparsers.add_parser(action, help=help_text)
        child.add_argument("args", nargs=argparse.REMAINDER)
    namespace = parser.parse_args()

    if namespace.action == "prepare":
        return _run("vasp_input_suite.py", namespace.args)
    return _run("vasp_submit_suite.py", namespace.args)


if __name__ == "__main__":
    raise SystemExit(main())
