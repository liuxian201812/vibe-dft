#!/bin/bash
set -euo pipefail

ROOT=$(cd "$(dirname "$0")/.." && pwd)
cd "$ROOT"

python3 - <<'PY'
from pathlib import Path
import re
import sys

patterns = [
    r"\b6138\b",
    r"\bcompute\b",
    r"VIBEDFT_LAB",
    r"VIBEDFT_MICRO",
    r"HPC_LAB",
    r"HPC_MICRO",
    r"/mnt/",
    r"/share/",
    r"/home/",
    r"/opt/",
    r"ckduan",
    r"wenjiaxu",
    r"c606",
    r"oneAPI",
    r"Singularity",
    r"singularity",
]

allow_dirs = {".git", ".venv"}
skip_prefixes = {Path("smoke_tests/tmp")}
bad = []
for path in Path(".").rglob("*"):
    if not path.is_file():
        continue
    if any(part in allow_dirs for part in path.parts):
        continue
    if any(prefix in path.parents for prefix in skip_prefixes):
        continue
    if path == Path("tools/check_no_site_details.sh"):
        continue
    if path.suffix in {".png", ".pdf", ".pack", ".idx", ".rev", ".pyc"}:
        continue
    try:
        text = path.read_text(encoding="utf-8")
    except Exception:
        continue
    for pattern in patterns:
        if re.search(pattern, text):
            bad.append((str(path), pattern))

if bad:
    for path, pattern in bad:
        print(f"forbidden pattern {pattern!r} in {path}")
    sys.exit(1)

print("No site-specific text found.")
PY
