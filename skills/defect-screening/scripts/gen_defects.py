#!/usr/bin/env python3
"""Compatibility CLI for the unified defect-generation pipeline.

The historical positional arguments are preserved. Generation, screening,
input writing, manifests, and direct-submit preparation all live in
defect_generation_pipeline.py.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

SCRIPT_DIR = Path(__file__).resolve().parent
if str(SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPT_DIR))
from defect_generation_pipeline import generate_from_config


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Generate screened defect VASP inputs"
    )
    parser.add_argument("poscar", help="Primitive-cell POSCAR")
    parser.add_argument("host_name", help="Host identifier")
    parser.add_argument(
        "--extrinsic",
        nargs="+",
        default=[],
        help="Extrinsic dopants or host-to-dopant mapping values",
    )
    parser.add_argument("--output", default="defects", help="Output directory")
    parser.add_argument("--verbose", "-v", action="store_true")
    parser.add_argument(
        "--functional",
        choices=("SCAN", "R2SCAN"),
        default="R2SCAN",
        help="Explicit meta-GGA used for relaxation",
    )
    parser.add_argument("--task-id", default="")
    parser.add_argument("--approval-id", default="")
    parser.add_argument("--server", default="")
    parser.add_argument("--remote-workdir", default="")
    args = parser.parse_args()

    poscar_path = Path(args.poscar).resolve()
    output_dir = Path(args.output).resolve()
    submission = {}
    if args.task_id:
        submission["task_id"] = args.task_id
    if args.approval_id:
        submission["approval_id"] = args.approval_id
    if args.remote_workdir:
        submission["remote_workdir"] = args.remote_workdir

    config = {
        "poscar": str(poscar_path),
        "host_id": args.host_name,
        "extrinsic": args.extrinsic,
        "output_dir": str(output_dir / args.host_name),
        "functional": args.functional,
        "verbose": args.verbose,
        "submission": submission,
    }
    if args.server:
        config["server"] = args.server
    generate_from_config(config, output_dir)


if __name__ == "__main__":
    try:
        main()
    except Exception as exc:
        print(f"error: {exc}", file=sys.stderr)
        raise
