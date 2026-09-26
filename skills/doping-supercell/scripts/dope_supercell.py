#!/usr/bin/env python3
"""Two-stage CLI for pristine supercell generation and exact-site doping.

Usage:
  dope_supercell.py build-supercell build_config.json
  dope_supercell.py dope doping_config.json
"""

import argparse
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

from core.dope_supercell import build_supercell, dope_existing_supercell


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        description="Build an auditable pristine supercell, then dope it in a separate step"
    )
    subparsers = parser.add_subparsers(dest="command", required=True)

    build_parser = subparsers.add_parser(
        "build-supercell",
        aliases=["build"],
        help="Stage 1: generate pristine supercell + manifest",
    )
    build_parser.add_argument("config", help="Stage-1 build config JSON")

    dope_parser = subparsers.add_parser(
        "dope",
        aliases=["substitute"],
        help="Stage 2: dope an accepted pristine supercell",
    )
    dope_parser.add_argument("config", help="Stage-2 doping config JSON")

    args = parser.parse_args(argv)
    if args.command in {"build-supercell", "build"}:
        build_supercell(args.config)
    else:
        dope_existing_supercell(args.config)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
