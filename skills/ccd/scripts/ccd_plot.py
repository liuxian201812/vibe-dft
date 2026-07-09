#!/usr/bin/env python3
"""Wrapper for vibedft CCD plotting."""
from __future__ import annotations

import subprocess
import sys
from pathlib import Path


def main():
    root = Path(__file__).resolve().parents[3]  # skills/ccd/scripts/ → repo root
    runner = root / "skills" / "vasp-postprocess" / "scripts" / "ccd_diagram.py"
    subprocess.run([sys.executable, str(runner), *sys.argv[1:]], check=True)


if __name__ == "__main__":
    main()
