#!/usr/bin/env python3
"""Check generic vibedft environment variables without printing secrets."""
from __future__ import annotations

import os


REQUIRED = [
    "VIBEDFT_VASP_CMD",
    "VIBEDFT_POTCAR_DIR",
]


def main() -> None:
    ok = True
    for key in REQUIRED:
        present = bool(os.environ.get(key))
        ok = ok and present
        print(f"{key}: {'set' if present else 'missing'}")
    raise SystemExit(0 if ok else 1)


if __name__ == "__main__":
    main()
