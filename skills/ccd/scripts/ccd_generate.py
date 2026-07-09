#!/usr/bin/env python3
"""Wrapper for vibedft CCD input generation."""
from __future__ import annotations

import argparse
import subprocess
import sys
from pathlib import Path


def main():
    parser = argparse.ArgumentParser(description="Generate CCD image inputs from two endpoint structures")
    parser.add_argument("config", help="CCD config JSON")
    args = parser.parse_args()
    root = Path(__file__).resolve().parents[3]  # skills/ccd/scripts/ → repo root
    runner = root / "skills" / "vasp-input" / "scripts" / "vasp_input_suite.py"
    subprocess.run([sys.executable, str(runner), "--phase", "ccd", "--config", args.config], check=True)


if __name__ == "__main__":
    main()
