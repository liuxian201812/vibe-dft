#!/usr/bin/env python3
"""Basic COD structure search and CIF download.

Usage:
  python3 search_cod.py search "SiO2"
  python3 search_cod.py download 9008567 -o SiO2.cif

Output is JSON. This version intentionally avoids Materials Project and
template substitution so it works without API keys.
"""
from __future__ import annotations

import argparse
import json
import urllib.parse
import urllib.request

COD_CIF = "https://www.crystallography.net/cod/{}.cif"
COD_RESULT = "https://www.crystallography.net/cod/result.php"


def _parse_elements(formula: str) -> list[str]:
    from pymatgen.core import Composition

    comp = Composition(formula)
    return [str(el) for el in comp.elements]


def _post_cod(params: dict) -> list[dict]:
    params = dict(params)
    params["format"] = "json"
    params.setdefault("include_theoretical", "0")
    data = urllib.parse.urlencode(params).encode()
    req = urllib.request.Request(COD_RESULT, data=data)
    with urllib.request.urlopen(req, timeout=30) as resp:
        return json.loads(resp.read().decode())


def search(formula: str) -> list[dict]:
    elements = _parse_elements(formula)
    params = {f"el{i + 1}": el for i, el in enumerate(elements[:8])}
    raw = _post_cod(params)
    results = []
    for entry in raw:
        if not isinstance(entry, dict) or not entry.get("file"):
            continue
        cod_id = entry["file"]
        results.append({
            "source": "cod",
            "id": cod_id,
            "formula": entry.get("formula", ""),
            "sg": entry.get("sg", ""),
            "a": float(entry.get("a", 0) or 0),
            "b": float(entry.get("b", 0) or 0),
            "c": float(entry.get("c", 0) or 0),
            "volume": float(entry.get("volume", 0) or 0),
            "z": int(entry.get("Z", 0) or 0),
            "cif_url": COD_CIF.format(cod_id),
        })
    return results


def download(cod_id: str, output: str | None = None) -> dict:
    output = output or f"{cod_id}.cif"
    with urllib.request.urlopen(COD_CIF.format(cod_id), timeout=30) as resp:
        cif_text = resp.read().decode()
    with open(output, "w") as f:
        f.write(cif_text)
    return {"cod_id": cod_id, "cif_file": output}


def main():
    parser = argparse.ArgumentParser(description="Basic COD search/download")
    sub = parser.add_subparsers(dest="mode", required=True)
    p_search = sub.add_parser("search", help="Search COD by formula elements")
    p_search.add_argument("formula")
    p_download = sub.add_parser("download", help="Download COD CIF by ID")
    p_download.add_argument("cod_id")
    p_download.add_argument("-o", "--output")
    args = parser.parse_args()

    if args.mode == "search":
        print(json.dumps(search(args.formula), indent=2))
    else:
        print(json.dumps(download(args.cod_id, args.output), indent=2))


if __name__ == "__main__":
    main()
