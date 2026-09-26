#!/usr/bin/env python3
"""vibedft submit helper.

Usage:
  python3 vasp_submit_suite.py --config task_config.json

This version only generates independent SLURM run scripts. It does not
build dependency chains, afterok links, daemons, queue fallbacks, or status
exporters.
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

_SCRIPT_DIR = Path(__file__).resolve().parent
if str(_SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(_SCRIPT_DIR))


def main():
    parser = argparse.ArgumentParser(description="Generate independent VASP submit scripts")
    parser.add_argument("--config", required=True, help="Path to JSON config")
    args = parser.parse_args()

    config_path = Path(args.config)
    if not config_path.exists():
        print(f"Config not found: {config_path}", file=sys.stderr)
        sys.exit(1)

    from generate_scripts import main as generate_main

    sys.argv = ["generate_scripts.py", str(config_path)]
    generate_main()


if __name__ == "__main__":
    main()
