#!/usr/bin/env python3
"""Sidecar metadata for generated and validated structures."""

from __future__ import annotations

import hashlib
import json
import os
import tempfile
from datetime import datetime
from pathlib import Path
from typing import Any


MANIFEST_NAME = "structure_manifest.json"


def now_iso() -> str:
    return datetime.now().astimezone().isoformat(timespec="seconds")


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def manifest_path(poscar_path: str | Path) -> Path:
    return Path(poscar_path).resolve().with_name(MANIFEST_NAME)


def load_manifest(poscar_path: str | Path) -> dict[str, Any]:
    path = manifest_path(poscar_path)
    if not path.is_file():
        return {}
    payload = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(payload, dict):
        raise ValueError(f"structure manifest must be an object: {path}")
    return payload


def atomic_write_json(path: Path, payload: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.NamedTemporaryFile(
        "w",
        encoding="utf-8",
        dir=path.parent,
        prefix=f".{path.name}.",
        suffix=".tmp",
        delete=False,
    ) as handle:
        tmp_path = Path(handle.name)
        json.dump(payload, handle, indent=2, ensure_ascii=True)
        handle.write("\n")
        handle.flush()
        os.fsync(handle.fileno())
    os.replace(tmp_path, path)


def write_manifest(
    poscar_path: str | Path,
    *,
    task_id: str = "",
    source: str = "",
    method: str = "",
    validation_status: str = "UNVALIDATED",
    validation: dict[str, Any] | None = None,
    chemistry: dict[str, Any] | None = None,
    space_group: dict[str, Any] | None = None,
    template: dict[str, Any] | None = None,
    extra: dict[str, Any] | None = None,
) -> dict[str, Any]:
    poscar = Path(poscar_path).resolve()
    current = load_manifest(poscar) if manifest_path(poscar).exists() else {}
    payload: dict[str, Any] = {
        **current,
        "schema_version": 1,
        "task_id": str(task_id or current.get("task_id") or ""),
        "structure_hash": sha256_file(poscar) if poscar.is_file() else "",
        "source": str(source or current.get("source") or ""),
        "method": str(method or current.get("method") or ""),
        "validation_status": str(
            validation_status or current.get("validation_status") or "UNVALIDATED"
        ).upper(),
        "validation": validation if validation is not None else current.get("validation", {}),
        "chemistry": chemistry if chemistry is not None else current.get("chemistry", {}),
        "space_group": (
            space_group if space_group is not None else current.get("space_group", {})
        ),
        "template": template if template is not None else current.get("template", {}),
        "updated_at": now_iso(),
    }
    if extra:
        payload.update(extra)
    atomic_write_json(manifest_path(poscar), payload)
    return payload
