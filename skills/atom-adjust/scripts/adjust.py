#!/usr/bin/env python3
"""adjust.py — POSCAR atom position editor (SnB complement).

Pure numpy, zero pymatgen dependency. Designed for manual structural tweaks
that ShakeNBreak does not cover: single-atom directed translation, rotation,
mirroring, absolute placement, and selective-dynamics toggles.

Complements SnB which handles automated batch bond-distortion scanning,
rattling, ground-state search, and analysis.
"""
from __future__ import annotations

import argparse
from decimal import Decimal, InvalidOperation
import hashlib
import json
import os
import re
import shutil
import stat
import sys
import tempfile
from pathlib import Path

import numpy as np


# ---------------------------------------------------------------------------
# Errors
# ---------------------------------------------------------------------------

class PoscarError(Exception):
    """Raised on malformed POSCAR files."""


def _affine_velocity_tail(tail: list[str], n_atoms: int) -> dict:
    """Accept only a mode line plus N exactly zero ion velocities."""
    if not any(line.strip() for line in tail):
        return {"kind": "none", "policy": "no_velocity_block"}
    mode_line = tail[0].strip().lower()
    if mode_line in ("", "c", "cartesian", "k"):
        mode = "Cartesian"
    elif mode_line in ("d", "direct"):
        mode = "Direct"
    else:
        raise PoscarError("affine unsupported velocity header/lattice velocities/predictor block")
    if len(tail) != n_atoms + 1:
        raise PoscarError("affine velocity block must have exactly one mode line and N rows")
    number = r"[+-]?(?:[0-9]+(?:\.[0-9]*)?|\.[0-9]+)(?:[Ee][+-]?[0-9]+)?"
    for row in tail[1:]:
        tokens = row.split()
        if len(tokens) != 3:
            raise PoscarError("affine requires exactly three velocity components per ion")
        for token in tokens:
            if not re.fullmatch(number, token):
                raise PoscarError("affine unsupported velocity numeric token")
            try:
                # Decimal construction is exact and does not underflow like float.
                is_zero = Decimal(token).is_zero()
            except InvalidOperation as exc:
                raise PoscarError("affine invalid velocity numeric token") from exc
            if not is_zero:
                raise PoscarError("affine only supports exactly zero ion velocities")
    return {
        "kind": "zero_ion_velocities",
        "policy": "preserve_original_tail_bytes",
        "coordinate_mode": mode,
        "mode_line": "blank" if mode_line == "" else "explicit",
        "rows": n_atoms,
        "exact_zero": True,
    }


# ---------------------------------------------------------------------------
# POSCAR I/O (line-preserving; keeps comments, SD flags, element rows)
# ---------------------------------------------------------------------------

def read_poscar(path: Path, *, strict_affine: bool = False) -> dict:
    """Parse a VASP POSCAR into a structured dict.

    Supports Direct and Cartesian coordinate modes, optional Selective
    dynamics rows, positive and negative scale factors (negative scale →
    volume mode, treated as raw lattice with a warning).
    strict_affine rejects unsupported scale/format variants instead of
    applying legacy warning-and-continue behavior.
    """
    source_bytes = path.read_bytes()
    lines = source_bytes.decode("utf-8").splitlines()
    if len(lines) < 8:
        raise PoscarError(f"POSCAR too short ({len(lines)} lines, need >=8)")

    comment = lines[0]

    try:
        scale = float(lines[1].split()[0])
    except (ValueError, IndexError) as exc:
        raise PoscarError(f"Invalid scale line: {lines[1]!r}") from exc
    if strict_affine and (
        len(lines[1].split()) != 1 or not np.isfinite(scale) or scale <= 0
    ):
        raise PoscarError("affine requires one finite positive scale; volume mode unsupported")

    try:
        lattice_raw = np.array(
            [list(map(float, lines[i].split())) for i in range(2, 5)]
        )
    except (ValueError, IndexError) as exc:
        raise PoscarError(f"Invalid lattice vectors (lines 3-5): {exc}") from exc
    if lattice_raw.shape != (3, 3):
        raise PoscarError(f"Lattice block wrong shape: {lattice_raw.shape}")
    if strict_affine:
        _positive_matrix(lattice_raw, "source lattice")

    # Line 5 may be element symbols or directly the counts.
    line5 = lines[5].split()
    try:
        counts = list(map(int, line5))
        elements: list[str] | None = None
        counts_idx = 5
    except ValueError:
        elements = line5
        try:
            counts = list(map(int, lines[6].split()))
        except ValueError as exc:
            raise PoscarError("Cannot locate atom counts line") from exc
        counts_idx = 6

    if strict_affine and (
        not elements or len(elements) != len(counts)
        or any(not re.fullmatch(r"[A-Z][a-z]?", item) for item in elements)
        or any(count < 0 for count in counts) or sum(counts) <= 0
    ):
        raise PoscarError("affine requires explicit species and matching nonnegative counts")
    n_total = sum(counts)

    # Detect "Selective dynamics" line, then Direct/Cartesian.
    idx = counts_idx + 1
    has_sd = False
    if idx < len(lines) and lines[idx].strip().lower().startswith("s"):
        has_sd = True
        idx += 1
    coord_type_idx = idx
    if idx >= len(lines):
        raise PoscarError("Missing Direct/Cartesian line")
    type_line = lines[idx].strip().lower()
    if type_line.startswith("d"):
        is_direct = True
    elif type_line.startswith(("c", "k")):
        is_direct = False
    else:
        raise PoscarError(f"Unknown coord-type line: {lines[idx]!r}")
    coord_start = idx + 1

    if len(lines) < coord_start + n_total:
        raise PoscarError(
            f"Not enough coordinate lines: need {n_total}, "
            f"have {len(lines) - coord_start}"
        )

    coords_listed = np.zeros((n_total, 3))
    sd_flags: list[list[str]] = []
    for i in range(n_total):
        parts = lines[coord_start + i].split()
        if len(parts) < 3:
            raise PoscarError(
                f"Bad coord line {coord_start + i + 1}: "
                f"{lines[coord_start + i]!r}"
            )
        try:
            coords_listed[i] = [float(parts[0]), float(parts[1]), float(parts[2])]
        except ValueError as exc:
            raise PoscarError(
                f"Bad coord line {coord_start + i + 1}: {exc}"
            ) from exc
        if has_sd:
            if strict_affine and (
                len(parts) < 6 or any(flag.upper() not in ("T", "F") for flag in parts[3:6])
            ):
                raise PoscarError("affine requires three explicit T/F flags for every atom")
            flags = parts[3:6] if len(parts) >= 6 else ["T", "T", "T"]
            sd_flags.append([str(f).upper() for f in flags])

    if strict_affine:
        if not np.isfinite(coords_listed).all():
            raise PoscarError("affine requires finite coordinates")
        velocity_tail = _affine_velocity_tail(lines[coord_start + n_total:], n_total)

    # Work internally in fractional coordinates with the actual lattice.
    lattice = scale * lattice_raw if scale > 0 else lattice_raw
    if scale < 0:
        print(
            "warning: negative scale (volume mode) — treating lattice as-is",
            file=sys.stderr,
        )

    if is_direct:
        coords_frac = coords_listed.copy()
    else:
        # Cartesian listed values are scaled by `scale` to get Å; conversion to
        # fractional cancels scale, so use lattice_raw.
        coords_frac = coords_listed @ np.linalg.inv(lattice_raw)

    result = {
        "path": path,
        "lines": lines,
        "comment": comment,
        "scale": scale,
        "lattice_raw": lattice_raw,
        "lattice": lattice,
        "elements": elements,
        "counts": counts,
        "n_total": n_total,
        "has_sd": has_sd,
        "is_direct": is_direct,
        "coord_type_idx": coord_type_idx,
        "coord_start": coord_start,
        "coords_frac": coords_frac,
        "sd_flags": sd_flags,
    }
    if strict_affine:
        _positive_matrix(lattice, "scaled source lattice")
        if not np.isfinite(coords_frac).all():
            raise PoscarError("nonfinite converted coordinates")
        result["source_sha256"] = hashlib.sha256(source_bytes).hexdigest()
        result["velocity_tail"] = velocity_tail
        if velocity_tail["kind"] == "zero_ion_velocities":
            velocity_tail["source_lines_1based"] = [coord_start + n_total + 1, len(lines)]
            result["velocity_tail_bytes"] = b"".join(
                source_bytes.splitlines(keepends=True)[coord_start + n_total:]
            )
    return result


def _positive_matrix(value, label: str) -> np.ndarray:
    """Reject invalid orientation and numerically unresolved inversion."""
    array = np.asarray(value)
    if array.shape != (3, 3) or array.dtype.kind not in "fiu":
        raise PoscarError(f"{label} must be a real numeric 3x3 matrix")
    matrix = array.astype(float)
    if not np.isfinite(matrix).all():
        raise PoscarError(f"{label} must be finite")
    with np.errstate(over="raise", invalid="raise", under="ignore"):
        determinant = float(np.linalg.det(matrix))
        condition = float(np.linalg.cond(matrix))
    if not np.isfinite(determinant) or determinant <= 0:
        raise PoscarError(f"{label} must be nonsingular and right-handed")
    if not np.isfinite(condition) or condition > 1e12:
        raise PoscarError(f"{label} is numerically singular (condition > 1e12)")
    return matrix


def affine_lattice(source: Path, target: Path, deformation) -> dict:
    """Active r'=F r at fixed fractional coordinates; never replace a file."""
    source, target = Path(source), Path(target)
    if source.resolve() == target.resolve():
        raise PoscarError("affine source and target must differ")
    if not stat.S_ISREG(source.lstat().st_mode):
        raise PoscarError("affine source must be a regular, non-symlink file")
    if os.path.lexists(target):
        raise FileExistsError("affine target already exists")
    matrix = _positive_matrix(deformation, "F")
    data = read_poscar(source, strict_affine=True)
    with np.errstate(over="raise", invalid="raise"):
        lattice = _positive_matrix(data["lattice"] @ matrix.T, "target lattice")

    lines = data["lines"][:]
    lines[1] = "1.0"
    lines[2:5] = [" ".join(format(x, ".17g") for x in row) for row in lattice]
    lines[data["coord_type_idx"]] = "Direct"
    # Keep Direct lines verbatim. For Cartesian input, replace only the three
    # numeric fields, preserving the original flags and any atom annotation.
    if not data["is_direct"]:
        for i, frac in enumerate(data["coords_frac"]):
            old_line = lines[data["coord_start"] + i]
            suffix = re.match(r"\s*\S+\s+\S+\s+\S+(.*)$", old_line).group(1)
            lines[data["coord_start"] + i] = " ".join(format(x, ".17g") for x in frac) + suffix
    payload = ("\n".join(lines) + "\n").encode("utf-8")
    if "velocity_tail_bytes" in data:
        prefix = "\n".join(lines[:data["coord_start"] + data["n_total"]]) + "\n"
        payload = prefix.encode("utf-8") + data["velocity_tail_bytes"]
    report = {
        "schema": "affine_lattice_v1",
        "operation": "active_deformation_fixed_fractional",
        "convention": "A_out = A_in @ F.T; row_lattice; Cartesian r_out = F @ r_in",
        "F": matrix.tolist(),
        "lattice_before_angstrom": data["lattice"].tolist(),
        "lattice_after_angstrom": lattice.tolist(),
        "volume_before_angstrom3": float(np.linalg.det(data["lattice"])),
        "volume_after_angstrom3": float(np.linalg.det(lattice)),
        "det_F": float(np.linalg.det(matrix)),
        "source_sha256": data["source_sha256"],
        "target_sha256": hashlib.sha256(payload).hexdigest(),
        "species": data["elements"],
        "counts": data["counts"],
        "atom_order_preserved": True,
        "selective_dynamics_preserved": True,
        "input_coordinate_mode": "Direct" if data["is_direct"] else "Cartesian",
        "output_coordinate_mode": "Direct",
        "fractional_wrapping": False,
        "velocity_tail": data["velocity_tail"],
    }
    # Hard-link publication is atomic and fails if another writer wins.
    # The temporary file is on the target filesystem; no overwrite fallback.
    temporary = None
    try:
        with tempfile.NamedTemporaryFile(dir=target.parent, prefix=".affine-", delete=False) as stream:
            temporary = Path(stream.name)
            stream.write(payload)
            stream.flush()
            os.fsync(stream.fileno())
        os.link(temporary, target)
    finally:
        if temporary is not None:
            temporary.unlink()
    return report


def cmd_affine(args: argparse.Namespace) -> int:
    try:
        report = affine_lattice(
            Path(args.poscar), Path(args.output), np.array(args.F).reshape(3, 3)
        )
    except (PoscarError, OSError, ValueError, IndexError,
            np.linalg.LinAlgError, FloatingPointError) as exc:
        print(f"affine error: {exc}", file=sys.stderr)
        return 2
    print(json.dumps(report, allow_nan=False, sort_keys=True))
    return 0


def _format_coord(listed: np.ndarray, sd_flags: list[str] | None,
                  precision: int = 10) -> str:
    coord_str = " ".join(f"{x:.{precision}f}" for x in listed)
    if sd_flags is not None:
        return f"{coord_str} {' '.join(sd_flags)}"
    return coord_str


def build_lines(data: dict, new_frac: np.ndarray,
                new_sd_flags: list[list[str]] | None = None) -> list[str]:
    """Return a fresh list of lines with updated fractional coords.

    If ``new_sd_flags`` is None the original SD flags (if any) are preserved.
    """
    lines = data["lines"][:]
    for i, frac in enumerate(new_frac):
        if data["is_direct"]:
            listed = frac
        else:
            listed = frac @ data["lattice_raw"]
        if new_sd_flags is not None:
            sd = new_sd_flags[i]
        elif data["has_sd"]:
            sd = data["sd_flags"][i] if i < len(data["sd_flags"]) else None
        else:
            sd = None
        lines[data["coord_start"] + i] = _format_coord(listed, sd)
    return lines


# ---------------------------------------------------------------------------
# Geometry helpers
# ---------------------------------------------------------------------------

def frac_to_cart(frac: np.ndarray, lattice: np.ndarray) -> np.ndarray:
    return frac @ lattice


def cart_to_frac(cart: np.ndarray, lattice: np.ndarray) -> np.ndarray:
    return cart @ np.linalg.inv(lattice)


def wrap_frac(frac: np.ndarray) -> np.ndarray:
    return frac % 1.0


def min_image_delta(frac_diff: np.ndarray, lattice: np.ndarray) -> np.ndarray:
    """Minimum-image Cartesian displacement for a fractional difference.

    Component-wise fractional wrapping is only sufficient for orthogonal
    cells. Search neighbouring lattice images so skewed supercells also use
    the physically shortest displacement.
    """
    centered = frac_diff - np.round(frac_diff)
    offsets = [(0, 0, 0)] + [
        (i, j, k)
        for i in (-1, 0, 1)
        for j in (-1, 0, 1)
        for k in (-1, 0, 1)
        if (i, j, k) != (0, 0, 0)
    ]
    candidates = np.array([centered + offset for offset in offsets])
    cart = candidates @ lattice
    return cart[np.argmin(np.einsum("ij,ij->i", cart, cart))]


def distance_pbc(frac_a: np.ndarray, frac_b: np.ndarray,
                 lattice: np.ndarray) -> float:
    return float(np.linalg.norm(min_image_delta(frac_a - frac_b, lattice)))


def nearest_neighbors(frac_query: np.ndarray, coords_frac: np.ndarray,
                      lattice: np.ndarray, exclude: set[int] | None = None,
                      n: int = 4) -> list[tuple[int, float]]:
    """Return the ``n`` nearest atoms (index, distance) by PBC distance."""
    exclude = exclude or set()
    dists = []
    for i in range(len(coords_frac)):
        if i in exclude:
            continue
        d = distance_pbc(frac_query, coords_frac[i], lattice)
        dists.append((i, d))
    dists.sort(key=lambda t: t[1])
    return dists[:n]


# ---------------------------------------------------------------------------
# Spec parsing
# ---------------------------------------------------------------------------

def parse_atoms(spec: str, n_total: int) -> list[int]:
    """Parse '5', '5,7,9', '5-10' → sorted unique 0-based indices."""
    indices: list[int] = []
    for part in spec.split(","):
        part = part.strip()
        if not part:
            continue
        if "-" in part:
            a, b = part.split("-", 1)
            indices.extend(range(int(a) - 1, int(b)))
        else:
            indices.append(int(part) - 1)
    if not indices:
        raise ValueError(f"No atoms parsed from {spec!r}")
    for idx in indices:
        if idx < 0 or idx >= n_total:
            raise ValueError(
                f"Atom index {idx + 1} out of range (1-{n_total})"
            )
    return sorted(set(indices))


def parse_direction(spec: str, coords_frac: np.ndarray,
                    lattice: np.ndarray,
                    centroid_frac: np.ndarray) -> np.ndarray:
    """Return a unit cartesian vector for the given direction spec.

    Specs:
      a/b/c            — lattice basis vector
      bond:j           — from moving centroid toward atom j
      vec:x,y,z        — explicit cartesian vector
      plane:i,j,k      — normal of plane defined by atoms i,j,k (right-handed)
      perp:j,ref:k     — component of (k - centroid) perpendicular to bond:j,
                          pointing toward ref:k
    """
    s = spec.strip().lower()

    if s in ("a", "b", "c"):
        axis = ord(s) - ord("a")
        v = lattice[axis]
        norm = np.linalg.norm(v)
        if norm < 1e-10:
            raise ValueError(f"Lattice vector {s!r} has zero norm")
        return v / norm

    if s.startswith("bond:"):
        j = int(s[5:].split(",")[0]) - 1
        v = min_image_delta(coords_frac[j] - centroid_frac, lattice)
        norm = np.linalg.norm(v)
        if norm < 1e-10:
            raise ValueError(
                f"bond:{j + 1} direction is zero (atom coincides with movers)"
            )
        return v / norm

    if s.startswith("vec:"):
        v = np.array([float(x) for x in s[4:].split(",")])
        if v.shape != (3,):
            raise ValueError(f"vec: needs 3 components, got {v}")
        norm = np.linalg.norm(v)
        if norm < 1e-10:
            raise ValueError("vec: direction has zero length")
        return v / norm

    if s.startswith("plane:"):
        idxs = [int(x) - 1 for x in s[6:].split(",")]
        if len(idxs) != 3:
            raise ValueError("plane: needs 3 atom indices (i,j,k)")
        p1 = frac_to_cart(coords_frac[idxs[0]], lattice)
        p2 = frac_to_cart(coords_frac[idxs[1]], lattice)
        p3 = frac_to_cart(coords_frac[idxs[2]], lattice)
        normal = np.cross(p2 - p1, p3 - p1)
        norm = np.linalg.norm(normal)
        if norm < 1e-10:
            raise ValueError("plane: points are collinear; no normal")
        return normal / norm

    if s.startswith("perp:"):
        rest = s[5:]
        parts = rest.split(",")
        if len(parts) != 2 or not parts[1].strip().lower().startswith("ref:"):
            raise ValueError("perp: format is 'perp:j,ref:k'")
        j = int(parts[0]) - 1
        k = int(parts[1].strip()[4:]) - 1
        v_bond = (frac_to_cart(coords_frac[j], lattice)
                  - frac_to_cart(centroid_frac, lattice))
        v_ref = (frac_to_cart(coords_frac[k], lattice)
                 - frac_to_cart(centroid_frac, lattice))
        nb = np.linalg.norm(v_bond)
        if nb < 1e-10:
            raise ValueError(
                f"perp:{j + 1}: bond direction zero (atom coincides with movers)"
            )
        b_hat = v_bond / nb
        perp = v_ref - np.dot(v_ref, b_hat) * b_hat  # Gram-Schmidt
        np_ = np.linalg.norm(perp)
        if np_ < 1e-10:
            raise ValueError(
                f"perp:{j + 1},ref:{k + 1}: ref atom lies on bond axis; "
                "perpendicular is undefined"
            )
        return perp / np_

    raise ValueError(f"Unknown direction spec: {spec!r}")


# ---------------------------------------------------------------------------
# Operations
# ---------------------------------------------------------------------------

def op_translate(data: dict, atoms: list[int], direction_spec: str,
                 dist: float | None, ratio: float | None) -> np.ndarray:
    """Translate atoms. Either --dist (Å along direction) or --ratio (bond: only).

    --ratio R: new_pos_i = r_j + R * (r_i - r_j).
        R=0.8 → bond shrinks to 80%, R=1.2 → bond grows to 120%.
        Equivalent to SnB's Bond_Distortion_(R-1).
    """
    coords = data["coords_frac"].copy()
    centroid = np.mean(coords[atoms], axis=0)

    if ratio is not None:
        if not direction_spec.strip().lower().startswith("bond:"):
            raise ValueError("--ratio requires a 'bond:j' direction")
        j = int(direction_spec.split(":")[1].split(",")[0]) - 1
        for i in atoms:
            bond_cart = min_image_delta(
                coords[i] - coords[j], data["lattice"]
            )
            scaled_frac = cart_to_frac(ratio * bond_cart, data["lattice"])
            coords[i] = wrap_frac(coords[j] + scaled_frac)
        return coords

    if dist is None:
        raise ValueError("either --dist or --ratio is required")

    unit_vec = parse_direction(direction_spec, coords, data["lattice"], centroid)
    disp_cart = unit_vec * dist
    disp_frac = cart_to_frac(disp_cart, data["lattice"])
    for i in atoms:
        coords[i] = wrap_frac(coords[i] + disp_frac)
    return coords


def op_rotate(data: dict, atoms: list[int], axis_i: int, axis_j: int,
              angle_deg: float, pivot: int | None = None) -> np.ndarray:
    """Rotate a group of atoms around the axis (i→j) by angle_deg (degrees),
    about pivot atom k (defaults to group centroid).

    Uses Rodrigues' rotation formula."""
    coords = data["coords_frac"].copy()
    lattice = data["lattice"]

    p_i = frac_to_cart(coords[axis_i], lattice)
    p_j = frac_to_cart(coords[axis_j], lattice)
    axis = p_j - p_i
    norm = np.linalg.norm(axis)
    if norm < 1e-10:
        raise ValueError(
            f"rotation axis atoms {axis_i + 1}-{axis_j + 1} coincide"
        )
    k_hat = axis / norm

    if pivot is not None:
        pivot_cart = frac_to_cart(coords[pivot], lattice)
    else:
        pivot_cart = np.mean([frac_to_cart(coords[i], lattice) for i in atoms],
                             axis=0)

    theta = np.radians(angle_deg)
    c, s = np.cos(theta), np.sin(theta)

    for i in atoms:
        v = frac_to_cart(coords[i], lattice) - pivot_cart
        rotated = (v * c
                   + np.cross(k_hat, v) * s
                   + np.outer(v @ k_hat, k_hat) * (1 - c))
        coords[i] = wrap_frac(cart_to_frac(rotated + pivot_cart, lattice))
    return coords


def op_mirror(data: dict, atoms: list[int], p_i: int, p_j: int,
              p_k: int) -> np.ndarray:
    """Mirror atoms across the plane defined by points i, j, k."""
    coords = data["coords_frac"].copy()
    lattice = data["lattice"]

    pi = frac_to_cart(coords[p_i], lattice)
    pj = frac_to_cart(coords[p_j], lattice)
    pk = frac_to_cart(coords[p_k], lattice)
    normal = np.cross(pj - pi, pk - pi)
    norm = np.linalg.norm(normal)
    if norm < 1e-10:
        raise ValueError(
            f"mirror plane atoms {p_i + 1},{p_j + 1},{p_k + 1} are collinear"
        )
    n_hat = normal / norm

    for i in atoms:
        v = frac_to_cart(coords[i], lattice) - pi
        reflected = v - 2.0 * np.dot(v, n_hat) * n_hat
        coords[i] = wrap_frac(cart_to_frac(reflected + pi, lattice))
    return coords


def op_set(data: dict, atoms: list[int], target: tuple[float, float, float],
           unit: str) -> np.ndarray:
    """Set absolute position. unit=frac → target is fractional;
    unit=ang → target is cartesian Å (converted via current lattice)."""
    coords = data["coords_frac"].copy()
    lattice = data["lattice"]
    if unit == "frac":
        new_frac = np.array(target, dtype=float)
    elif unit == "ang":
        new_frac = cart_to_frac(np.array(target, dtype=float), lattice)
    else:
        raise ValueError(f"Unknown unit: {unit!r}")
    # Set the first moving atom exactly; offset the rest rigidly so the
    # group's internal geometry is preserved.
    base_old = coords[atoms[0]]
    delta = new_frac - base_old
    for i in atoms:
        coords[i] = wrap_frac(coords[i] + delta)
    return coords


def op_set_sd(data: dict, atoms: list[int], free: bool) -> tuple[np.ndarray, list[list[str]]]:
    """Toggle selective-dynamics flags. If the file has no SD row, one is
    inserted (all atoms default to free=T/T/T, then `atoms` set to fix=F/F/F
    when free=False)."""
    flag = "T" if free else "F"
    n = data["n_total"]
    sd = [list(data["sd_flags"][i]) if data["has_sd"] and i < len(data["sd_flags"])
          else ["T", "T", "T"]
          for i in range(n)]
    for i in atoms:
        sd[i] = [flag, flag, flag]
    return data["coords_frac"].copy(), sd


# ---------------------------------------------------------------------------
# Reporting
# ---------------------------------------------------------------------------

def describe_displacement(data: dict, atoms: list[int],
                          new_coords: np.ndarray) -> str:
    """Human-readable summary of what changed (for --dry-run or verbose)."""
    lattice = data["lattice"]
    lines_out = []
    for i in atoms:
        old = data["coords_frac"][i]
        new = new_coords[i]
        cart_old = frac_to_cart(old, lattice)
        cart_new = frac_to_cart(new, lattice)
        delta_cart = min_image_delta(new - old, lattice)
        d = np.linalg.norm(delta_cart)
        lines_out.append(
            f"  atom {i + 1:3d}  frac {old.round(6)} → {new.round(6)}  "
            f"(Δ={d:.4f} Å, dir={delta_cart.round(4)})"
        )
    return "\n".join(lines_out)


def describe_bonds(data: dict, atoms: list[int], new_coords: np.ndarray,
                   n: int = 4) -> str:
    """Show nearest-neighbour distances before and after the move."""
    lattice = data["lattice"]
    lines_out = []
    for i in atoms:
        before = nearest_neighbors(
            data["coords_frac"][i], data["coords_frac"], lattice,
            exclude={i}, n=n,
        )
        lines_out.append(f"  atom {i + 1} nearest neighbours (Å):")
        # Preserve the before-neighbour identity when reporting the after
        # distance.  Independently sorting both lists and zipping them can
        # silently relabel a different atom when neighbour ranks change.
        for j_b, d_b in before:
            d_a = distance_pbc(new_coords[i], new_coords[j_b], lattice)
            arrow = "→"
            lines_out.append(f"    atom {j_b + 1:3d}: {d_b:.4f} {arrow} {d_a:.4f}")
    return "\n".join(lines_out)


# ---------------------------------------------------------------------------
# Output dispatch
# ---------------------------------------------------------------------------

def emit(data: dict, new_lines: list[str], args: argparse.Namespace) -> int:
    text = "\n".join(new_lines) + "\n"
    if args.stdout:
        sys.stdout.write(text)
        return 0
    target = Path(args.output) if args.output else data["path"]
    if args.backup and not args.output and not args.no_backup is False:
        pass  # placeholder; handled below
    if args.backup and target == data["path"] and not args.output:
        bak = data["path"].with_suffix(data["path"].suffix + ".bak")
        if not bak.exists():
            shutil.copy2(data["path"], bak)
            print(f"backup → {bak}", file=sys.stderr)
    target.write_text(text)
    if not args.output:
        print(f"wrote {data['path']}", file=sys.stderr)
    else:
        print(f"wrote {target}", file=sys.stderr)
    return 0


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def _add_common_args(p: argparse.ArgumentParser) -> None:
    p.add_argument("-i", "--inplace", action="store_true",
                   help="edit POSCAR in place (default if no -o/-c)")
    p.add_argument("-o", "--output", default=None,
                   help="write to a new file instead of in-place")
    p.add_argument("-c", "--stdout", action="store_true",
                   help="write result to stdout, do not modify any file")
    p.add_argument("--backup", dest="backup", action="store_true",
                   default=True, help="write .bak next to source (default)")
    p.add_argument("--no-backup", dest="backup", action="store_false",
                   help="skip writing .bak")
    p.add_argument("--dry-run", action="store_true",
                   help="print planned changes, do not modify files")
    p.add_argument("--print-bonds", type=int, default=4, metavar="N",
                   help="show N nearest-neighbour distances before/after (default 4)")
    p.add_argument("-v", "--verbose", action="store_true")
    p.add_argument("poscar", nargs="?", default="POSCAR",
                   help="POSCAR path (default: ./POSCAR)")


def _resolve_output(args: argparse.Namespace) -> None:
    """Normalise so that exactly one of stdout/output/inplace is active."""
    if args.stdout:
        return
    if args.output:
        return
    # default: inplace
    args.inplace = True


def cmd_translate(args: argparse.Namespace) -> int:
    data = read_poscar(Path(args.poscar))
    atoms = parse_atoms(args.atoms, data["n_total"])
    new_coords = op_translate(data, atoms, args.dir, args.dist, args.ratio)

    if args.verbose or args.dry_run:
        print(f"translate atoms {args.atoms} dir={args.dir} "
              f"{'ratio=' + str(args.ratio) if args.ratio is not None else 'dist=' + str(args.dist)}",
              file=sys.stderr)
        print(describe_displacement(data, atoms, new_coords), file=sys.stderr)
        if args.print_bonds > 0:
            print(describe_bonds(data, atoms, new_coords, args.print_bonds),
                  file=sys.stderr)
    if args.dry_run:
        return 0

    new_lines = build_lines(data, new_coords)
    return _write_output(data, new_lines, args)


def cmd_rotate(args: argparse.Namespace) -> int:
    data = read_poscar(Path(args.poscar))
    atoms = parse_atoms(args.atoms, data["n_total"])
    axis_i, axis_j = [int(x) - 1 for x in args.axis.split(",")]
    pivot = int(args.pivot) - 1 if args.pivot else None
    new_coords = op_rotate(data, atoms, axis_i, axis_j, args.angle, pivot)

    if args.verbose or args.dry_run:
        print(f"rotate atoms {args.atoms} axis={args.axis} angle={args.angle}° "
              f"pivot={args.pivot or 'centroid'}", file=sys.stderr)
        print(describe_displacement(data, atoms, new_coords), file=sys.stderr)
    if args.dry_run:
        return 0

    new_lines = build_lines(data, new_coords)
    return _write_output(data, new_lines, args)


def cmd_mirror(args: argparse.Namespace) -> int:
    data = read_poscar(Path(args.poscar))
    atoms = parse_atoms(args.atoms, data["n_total"])
    p_i, p_j, p_k = [int(x) - 1 for x in args.plane.split(",")]
    new_coords = op_mirror(data, atoms, p_i, p_j, p_k)

    if args.verbose or args.dry_run:
        print(f"mirror atoms {args.atoms} plane={args.plane}", file=sys.stderr)
        print(describe_displacement(data, atoms, new_coords), file=sys.stderr)
    if args.dry_run:
        return 0

    new_lines = build_lines(data, new_coords)
    return _write_output(data, new_lines, args)


def cmd_set(args: argparse.Namespace) -> int:
    data = read_poscar(Path(args.poscar))
    atoms = parse_atoms(args.atoms, data["n_total"])
    target = tuple(float(x) for x in args.to.split(","))
    if len(target) != 3:
        raise ValueError("--to needs x,y,z")
    new_coords = op_set(data, atoms, target, args.unit)

    if args.verbose or args.dry_run:
        print(f"set atoms {args.atoms} to {target} ({args.unit})", file=sys.stderr)
        print(describe_displacement(data, atoms, new_coords), file=sys.stderr)
    if args.dry_run:
        return 0

    new_lines = build_lines(data, new_coords)
    return _write_output(data, new_lines, args)


def cmd_fix(args: argparse.Namespace) -> int:
    return _cmd_sd(args, free=False)


def cmd_free(args: argparse.Namespace) -> int:
    return _cmd_sd(args, free=True)


def _cmd_sd(args: argparse.Namespace, free: bool) -> int:
    data = read_poscar(Path(args.poscar))
    atoms = parse_atoms(args.atoms, data["n_total"])
    new_coords, new_sd = op_set_sd(data, atoms, free)

    # If we just inserted an SD row, we must rewrite the whole file with the
    # "Selective dynamics" line added.
    if not data["has_sd"]:
        lines = data["lines"][:]
        lines.insert(data["coord_type_idx"], "Selective dynamics")
        # Rebuild data with updated indices
        data = dict(data)
        data["lines"] = lines
        data["has_sd"] = True
        data["coord_type_idx"] += 1
        data["coord_start"] += 1
        data["sd_flags"] = [["T", "T", "T"]] * data["n_total"]
        # Now apply the fix/free to the freshly-seeded SD flags
        new_sd = [list(f) for f in new_sd]
        for i in atoms:
            flag = "T" if free else "F"
            new_sd[i] = [flag, flag, flag]

    if args.verbose or args.dry_run:
        action = "free" if free else "fix"
        print(f"{action} atoms {args.atoms}", file=sys.stderr)
    if args.dry_run:
        return 0

    new_lines = build_lines(data, new_coords, new_sd_flags=new_sd)
    return _write_output(data, new_lines, args)


def _write_output(data: dict, new_lines: list[str],
                  args: argparse.Namespace) -> int:
    text = "\n".join(new_lines) + "\n"
    if args.stdout:
        sys.stdout.write(text)
        return 0
    target = Path(args.output) if args.output else data["path"]
    if args.backup and target == data["path"] and not args.dry_run:
        bak = data["path"].parent / (data["path"].name + ".bak")
        if not bak.exists():
            shutil.copy2(data["path"], bak)
            print(f"backup → {bak}", file=sys.stderr)
    target.write_text(text)
    print(f"wrote {target}", file=sys.stderr)
    return 0


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="adjust.py",
        description="POSCAR atom-position editor (SnB complement).",
    )
    sub = parser.add_subparsers(dest="command", required=True)

    p_a = sub.add_parser("affine", help="active lattice deformation at fixed fractional coordinates")
    p_a.add_argument("--F", type=float, nargs=9, required=True,
                     help="nine row-major entries of dimensionless Cartesian F; A_out=A_in@F.T")
    p_a.add_argument("-o", "--output", required=True, help="new output file; never overwrite")
    p_a.add_argument("poscar", help="explicit source POSCAR")
    p_a.set_defaults(func=cmd_affine)

    p_t = sub.add_parser("translate", help="translate atoms along a direction")
    p_t.add_argument("--atoms", required=True,
                     help="atom spec, e.g. 5 / 5,7,9 / 5-10 (1-based)")
    p_t.add_argument("--dir", required=True,
                     help="a/b/c | bond:j | vec:x,y,z | plane:i,j,k | perp:j,ref:k")
    g = p_t.add_mutually_exclusive_group(required=True)
    g.add_argument("--dist", type=float, help="displacement magnitude (Å)")
    g.add_argument("--ratio", type=float,
                   help="target bond length / original (requires bond:j). "
                        "0.8 → 80%%, 1.2 → 120%%")
    _add_common_args(p_t)
    p_t.set_defaults(func=cmd_translate)

    p_r = sub.add_parser("rotate", help="rotate atoms about an axis")
    p_r.add_argument("--atoms", required=True)
    p_r.add_argument("--axis", required=True,
                     help="two atom indices defining the axis: i,j")
    p_r.add_argument("--angle", type=float, required=True,
                     help="rotation angle in degrees (right-hand rule)")
    p_r.add_argument("--pivot", default=None,
                     help="pivot atom index (default: group centroid)")
    _add_common_args(p_r)
    p_r.set_defaults(func=cmd_rotate)

    p_m = sub.add_parser("mirror", help="mirror atoms across a plane")
    p_m.add_argument("--atoms", required=True)
    p_m.add_argument("--plane", required=True,
                     help="three atom indices defining the plane: i,j,k")
    _add_common_args(p_m)
    p_m.set_defaults(func=cmd_mirror)

    p_s = sub.add_parser("set", help="set absolute atom position")
    p_s.add_argument("--atoms", required=True)
    p_s.add_argument("--to", required=True, help="target x,y,z")
    p_s.add_argument("--unit", choices=["ang", "frac"], default="ang",
                     help="units of --to (default: ang)")
    _add_common_args(p_s)
    p_s.set_defaults(func=cmd_set)

    p_fix = sub.add_parser("fix", help="freeze atoms (Selective dynamics F/F/F)")
    p_fix.add_argument("--atoms", required=True)
    _add_common_args(p_fix)
    p_fix.set_defaults(func=cmd_fix)

    p_free = sub.add_parser("free", help="free atoms (Selective dynamics T/T/T)")
    p_free.add_argument("--atoms", required=True)
    _add_common_args(p_free)
    p_free.set_defaults(func=cmd_free)

    return parser


def main(argv: list[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    if args.command != "affine":
        _resolve_output(args)
    try:
        return args.func(args)
    except (PoscarError, ValueError) as exc:
        print(f"error: {exc}", file=sys.stderr)
        return 2


if __name__ == "__main__":
    sys.exit(main())
