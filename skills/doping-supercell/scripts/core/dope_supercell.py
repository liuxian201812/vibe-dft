"""Generate doped supercells with optional strict site provenance.

The legacy path still supports doped's automatic supercell search and the
historical ``count``/``all``/zero-based ``sites`` substitution forms.  The
strict path maps the target back to the supplied parent structure and validates
the written POSCAR.  An exact supercell index remains authoritative when
provided; otherwise the nearest translation-equivalent site to the periodic
supercell center is selected deterministically.
"""

from __future__ import annotations

import hashlib
import json
from pathlib import Path
from typing import Any

import numpy as np
from pymatgen.core import Composition, Structure
from pymatgen.io.vasp import Poscar
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer


_SUPERCELL_CENTER_FRAC = np.array([0.5, 0.5, 0.5], dtype=float)
_CENTER_DISTANCE_TIE_TOLERANCE_A = 1e-10


def _validate_supercell_matrix(supercell_matrix: Any) -> np.ndarray:
    """Return a validated 3x3 nonsingular integer supercell matrix."""
    raw = np.asarray(supercell_matrix)
    if raw.shape != (3, 3):
        raise ValueError("supercell_matrix must be a 3x3 integer matrix")
    if not np.issubdtype(raw.dtype, np.number) or not np.all(np.isfinite(raw)):
        raise ValueError("supercell_matrix must contain finite numbers")
    rounded = np.rint(raw).astype(int)
    if not np.allclose(raw.astype(float), rounded, atol=0, rtol=0):
        raise ValueError("supercell_matrix entries must be exact integers")
    determinant = float(np.linalg.det(rounded))
    if abs(determinant) < 0.5:
        raise ValueError("supercell_matrix must be nonsingular")
    return rounded


def generate_supercell(
    structure: Structure,
    min_image_distance: float = 12.0,
    min_atoms: int = 80,
    supercell_matrix: Any | None = None,
    allow_primitive: bool = True,
) -> tuple[Structure, np.ndarray]:
    """Generate a supercell.

    When ``supercell_matrix`` is supplied, it is applied directly to the input
    parent structure.  No primitive-cell conversion or automatic search is
    performed, which makes parent-site provenance deterministic.

    Otherwise the legacy doped-based automatic search is used.
    """
    if supercell_matrix is not None:
        matrix = _validate_supercell_matrix(supercell_matrix)
        supercell = structure.copy()
        supercell.make_supercell(matrix)
        expected_atoms = len(structure) * abs(int(round(np.linalg.det(matrix))))
        if len(supercell) != expected_atoms:
            raise ValueError(
                "Explicit supercell atom-count mismatch: "
                f"expected {expected_atoms}, generated {len(supercell)}"
            )
        expected_lattice = matrix @ np.asarray(structure.lattice.matrix)
        if not np.allclose(supercell.lattice.matrix, expected_lattice, atol=1e-8, rtol=0):
            raise ValueError("Explicit supercell lattice does not match M @ parent_lattice")
        return supercell, matrix

    from doped.generation import get_ideal_supercell_matrix

    def try_expand(s: Structure) -> tuple[Structure, np.ndarray] | None:
        current_min_img = _compute_min_image(s.lattice.matrix)
        if len(s) >= min_atoms and current_min_img >= min_image_distance:
            return (s, np.eye(3, dtype=int))
        mat = get_ideal_supercell_matrix(
            s, min_image_distance=min_image_distance, min_atoms=min_atoms,
        )
        if mat is not None:
            return (s * mat, np.asarray(mat, dtype=int))
        return None

    best = try_expand(structure)

    # Legacy callers can still opt into primitive-cell comparison.  The
    # auditable two-stage builder disables it so M always maps the generated
    # supercell directly back to the exact supplied parent.
    if allow_primitive:
        try:
            sga = SpacegroupAnalyzer(structure, symprec=0.2)
            prim = sga.get_primitive_standard_structure()
            if len(prim) != len(structure):
                prim_result = try_expand(prim)
                if prim_result and (best is None or len(prim_result[0]) < len(best[0])):
                    best = prim_result
        except Exception:
            pass

    if best is None:
        raise ValueError(
            f"Could not find supercell meeting min_image_distance={min_image_distance} "
            f"and min_atoms={min_atoms}"
        )

    return best


def _compute_min_image(lattice_matrix: np.ndarray) -> float:
    """Compute minimum distance between periodic images."""
    min_d = float("inf")
    for i in range(-2, 3):
        for j in range(-2, 3):
            for k in range(-2, 3):
                if i == 0 and j == 0 and k == 0:
                    continue
                vector = i * lattice_matrix[0] + j * lattice_matrix[1] + k * lattice_matrix[2]
                min_d = min(min_d, float(np.linalg.norm(vector)))
    return min_d


def _require_1based_index(value: Any, *, field: str, size: int) -> int:
    """Validate an exact integer 1-based index and return its 0-based value."""
    if isinstance(value, bool) or not isinstance(value, (int, np.integer)):
        raise ValueError(f"{field} must be an integer using 1-based indexing")
    if value < 1 or value > size:
        raise ValueError(f"{field}={value} is outside the valid 1..{size} range")
    return int(value) - 1


def _parent_wyckoff_labels(parent: Structure, symprec: float) -> list[str] | None:
    """Return full per-site Wyckoff symbols (for example 4c), not bare letters."""
    try:
        symmetrized = SpacegroupAnalyzer(
            parent, symprec=symprec
        ).get_symmetrized_structure()
        labels: list[str | None] = [None] * len(parent)
        for symbol, group in zip(
            symmetrized.wyckoff_symbols,
            symmetrized.equivalent_indices,
            strict=True,
        ):
            for index in group:
                labels[index] = str(symbol)
        return [label for label in labels if label is not None] if all(labels) else None
    except Exception:
        return None


def _periodic_distance(lattice, frac_a: Any, frac_b: Any) -> float:
    return float(lattice.get_all_distances([frac_a], [frac_b])[0, 0])


def map_supercell_site_to_parent(
    parent: Structure,
    supercell: Structure,
    supercell_matrix: Any,
    site_1based: int,
    *,
    expected_species: str | None = None,
    expected_parent_index_1based: int | None = None,
    expected_parent_fractional_coordinate: Any | None = None,
    expected_parent_wyckoff: str | None = None,
    mapping_tolerance: float = 1e-5,
    symprec: float = 0.2,
    precomputed_wyckoff_labels: list[str] | None = None,
) -> dict:
    """Map one exact supercell site back to its parent site and validate it.

    Fractional coordinates are row vectors, so the inverse provenance mapping
    is exactly ``f_parent = f_supercell @ M mod 1``.
    """
    matrix = _validate_supercell_matrix(supercell_matrix)
    site_index = _require_1based_index(
        site_1based, field="site_1based", size=len(supercell)
    )
    if mapping_tolerance <= 0:
        raise ValueError("mapping_tolerance must be positive")

    selected = supercell[site_index]
    if expected_species is not None and selected.species_string != expected_species:
        raise ValueError(
            f"site_1based={site_1based} is {selected.species_string}, "
            f"expected {expected_species}"
        )

    f_supercell = np.asarray(selected.frac_coords, dtype=float)
    f_parent = np.mod(f_supercell @ matrix, 1.0)
    distances = np.asarray(
        parent.lattice.get_all_distances([f_parent], parent.frac_coords)[0],
        dtype=float,
    )
    parent_index = int(np.argmin(distances))
    mapping_distance = float(distances[parent_index])
    if mapping_distance > mapping_tolerance:
        raise ValueError(
            "Supercell target does not map to a parent atom within tolerance: "
            f"nearest distance={mapping_distance:.8g} Å, tolerance={mapping_tolerance:.8g} Å"
        )

    tied = np.flatnonzero(np.isclose(distances, mapping_distance, atol=1e-10, rtol=0))
    if len(tied) != 1:
        raise ValueError(
            "Supercell target maps ambiguously to multiple parent atoms: "
            + ", ".join(str(int(i) + 1) for i in tied)
        )

    parent_site = parent[parent_index]
    if parent_site.species_string != selected.species_string:
        raise ValueError(
            "Mapped parent species mismatch: "
            f"supercell site is {selected.species_string}, parent site "
            f"{parent_index + 1} is {parent_site.species_string}"
        )

    if expected_parent_index_1based is not None:
        expected_parent_index = _require_1based_index(
            expected_parent_index_1based,
            field="expected_parent_index_1based",
            size=len(parent),
        )
        if parent_index != expected_parent_index:
            raise ValueError(
                f"Mapped parent site is {parent_index + 1}, "
                f"expected {expected_parent_index_1based}"
            )

    if expected_parent_fractional_coordinate is not None:
        expected_frac = np.asarray(expected_parent_fractional_coordinate, dtype=float)
        if expected_frac.shape != (3,) or not np.all(np.isfinite(expected_frac)):
            raise ValueError("expected_parent_fractional_coordinate must contain 3 finite numbers")
        expected_distance = _periodic_distance(parent.lattice, f_parent, expected_frac)
        if expected_distance > mapping_tolerance:
            raise ValueError(
                "Mapped parent fractional coordinate differs from the expected coordinate: "
                f"distance={expected_distance:.8g} Å"
            )

    wyckoff_labels = (
        precomputed_wyckoff_labels
        if precomputed_wyckoff_labels is not None
        else _parent_wyckoff_labels(parent, symprec)
    )
    if wyckoff_labels is not None and len(wyckoff_labels) != len(parent):
        raise ValueError("Precomputed parent Wyckoff labels have the wrong length")
    parent_wyckoff = wyckoff_labels[parent_index] if wyckoff_labels else None
    if expected_parent_wyckoff is not None:
        if parent_wyckoff is None:
            raise ValueError("Could not determine parent Wyckoff labels for validation")
        if str(parent_wyckoff) != str(expected_parent_wyckoff):
            raise ValueError(
                f"Mapped parent Wyckoff is {parent_wyckoff}, "
                f"expected {expected_parent_wyckoff}"
            )

    return {
        "supercell_site_1based": site_1based,
        "supercell_species": selected.species_string,
        "supercell_fractional_coordinate": f_supercell.tolist(),
        "mapping_formula": "f_parent = f_supercell @ M mod 1",
        "mapped_parent_fractional_coordinate": f_parent.tolist(),
        "mapped_parent_site": parent_index + 1,
        "parent_atom_index": parent_index + 1,
        "parent_atom_index_1based": parent_index + 1,
        "parent_species": parent_site.species_string,
        "parent_fractional_coordinate": np.asarray(parent_site.frac_coords).tolist(),
        "parent_wyckoff": parent_wyckoff,
        "mapping_distance": mapping_distance,
        "mapping_distance_A": mapping_distance,
        "mapping_tolerance_A": mapping_tolerance,
        "validation": {
            "expected_parent_index_1based": expected_parent_index_1based,
            "expected_parent_fractional_coordinate": (
                list(expected_parent_fractional_coordinate)
                if expected_parent_fractional_coordinate is not None
                else None
            ),
            "expected_parent_wyckoff": expected_parent_wyckoff,
            "passed": True,
        },
    }


def _center_distance_record(supercell: Structure, site_index: int) -> dict:
    """Describe direct and periodic distances from one site to the cell center.

    Selection uses the direct Cartesian displacement inside the written POSCAR
    parallelepiped from the in-cell fractional coordinate to ``(0.5, 0.5,
    0.5)``. No periodic wrapping is applied to that ranking distance. A
    periodic shortest-image distance is recorded separately as a diagnostic.
    """
    site = supercell[site_index]
    frac_in_cell = np.mod(np.asarray(site.frac_coords, dtype=float), 1.0)

    direct_displacement_frac = frac_in_cell - _SUPERCELL_CENTER_FRAC
    direct_displacement_cart = supercell.lattice.get_cartesian_coords(
        direct_displacement_frac
    )
    direct_distance = float(np.linalg.norm(direct_displacement_cart))

    periodic_distance, periodic_image = supercell.lattice.get_distance_and_image(
        _SUPERCELL_CENTER_FRAC, frac_in_cell
    )
    periodic_image = np.asarray(periodic_image, dtype=int)
    periodic_displacement_frac = (
        frac_in_cell + periodic_image - _SUPERCELL_CENTER_FRAC
    )
    periodic_displacement_cart = supercell.lattice.get_cartesian_coords(
        periodic_displacement_frac
    )
    return {
        "site_1based": site_index + 1,
        "site_0based": site_index,
        "fractional_coordinate": frac_in_cell.tolist(),
        "direct_center_displacement_fractional": direct_displacement_frac.tolist(),
        "direct_center_displacement_cartesian_A": np.asarray(
            direct_displacement_cart
        ).tolist(),
        "direct_distance_to_center_A": direct_distance,
        "periodic_nearest_image_translation": periodic_image.tolist(),
        "periodic_center_displacement_fractional": periodic_displacement_frac.tolist(),
        "periodic_center_displacement_cartesian_A": np.asarray(
            periodic_displacement_cart
        ).tolist(),
        "periodic_shortest_image_distance_to_center_A": float(periodic_distance),
    }


def _matching_translation_candidates(
    parent: Structure,
    supercell: Structure,
    supercell_matrix: Any,
    *,
    expected_species: str,
    expected_parent_index_1based: int,
    expected_parent_fractional_coordinate: Any | None,
    expected_parent_wyckoff: str,
    mapping_tolerance: float,
    symprec: float,
) -> list[dict]:
    """Return every supercell translation satisfying the parent-site constraints."""
    candidates: list[dict] = []
    expected_frac = (
        np.asarray(expected_parent_fractional_coordinate, dtype=float)
        if expected_parent_fractional_coordinate is not None
        else None
    )
    if expected_frac is not None and (
        expected_frac.shape != (3,) or not np.all(np.isfinite(expected_frac))
    ):
        raise ValueError(
            "expected_parent_fractional_coordinate must contain 3 finite numbers"
        )
    wyckoff_labels = _parent_wyckoff_labels(parent, symprec)
    if wyckoff_labels is None:
        raise ValueError("Could not determine parent Wyckoff labels for site selection")
    for site_index, site in enumerate(supercell):
        if site.species_string != expected_species:
            continue
        mapping = map_supercell_site_to_parent(
            parent,
            supercell,
            supercell_matrix,
            site_index + 1,
            expected_species=expected_species,
            mapping_tolerance=mapping_tolerance,
            symprec=symprec,
            precomputed_wyckoff_labels=wyckoff_labels,
        )
        if mapping["mapped_parent_site"] != expected_parent_index_1based:
            continue
        if str(mapping["parent_wyckoff"]) != str(expected_parent_wyckoff):
            continue
        if expected_frac is not None:
            expected_distance = _periodic_distance(
                parent.lattice,
                mapping["mapped_parent_fractional_coordinate"],
                expected_frac,
            )
            if expected_distance > mapping_tolerance:
                continue

        record = _center_distance_record(supercell, site_index)
        record.update(
            {
                "species": site.species_string,
                "mapped_parent_site_1based": mapping["mapped_parent_site"],
                "mapped_parent_wyckoff": mapping["parent_wyckoff"],
                "mapping_distance_A": mapping["mapping_distance_A"],
            }
        )
        candidates.append(record)

    if not candidates:
        raise ValueError(
            "No pristine-supercell site satisfies the requested element, parent-site, "
            "Wyckoff, and optional parent-coordinate constraints"
        )
    return sorted(candidates, key=lambda item: item["site_0based"])


def _select_nearest_center_candidate(candidates: list[dict]) -> tuple[dict, list[int]]:
    """Select by direct in-parallelepiped Cartesian distance, then lowest index."""
    if not candidates:
        raise ValueError("At least one center-selection candidate is required")
    minimum_distance = min(
        candidate["direct_distance_to_center_A"] for candidate in candidates
    )
    tied_sites = sorted(
        candidate["site_1based"]
        for candidate in candidates
        if np.isclose(
            candidate["direct_distance_to_center_A"],
            minimum_distance,
            atol=_CENTER_DISTANCE_TIE_TOLERANCE_A,
            rtol=0,
        )
    )
    selected_site = tied_sites[0]
    selected = next(
        candidate for candidate in candidates if candidate["site_1based"] == selected_site
    )
    return selected, tied_sites


def _automatic_site_selection(
    parent: Structure,
    supercell: Structure,
    supercell_matrix: Any,
    *,
    selection_mode_defaulted: bool,
    expected_species: str,
    parent_validation: dict,
) -> tuple[int, dict]:
    """Enumerate matching translations and select the direct-center nearest one."""
    candidates = _matching_translation_candidates(
        parent,
        supercell,
        supercell_matrix,
        expected_species=expected_species,
        expected_parent_index_1based=parent_validation[
            "expected_parent_index_1based"
        ],
        expected_parent_fractional_coordinate=parent_validation[
            "expected_parent_fractional_coordinate"
        ],
        expected_parent_wyckoff=parent_validation["expected_parent_wyckoff"],
        mapping_tolerance=parent_validation["mapping_tolerance"],
        symprec=parent_validation["symprec"],
    )
    selected, tied_sites = _select_nearest_center_candidate(candidates)
    selected_site_1based = selected["site_1based"]
    selection = {
        "mode": "nearest_supercell_center",
        "mode_defaulted": bool(selection_mode_defaulted),
        "candidate_enumeration_performed": True,
        "explicit_site_preserved": False,
        "rationale": (
            "Selected the matching translation-equivalent site with the smallest "
            "direct Cartesian distance inside the written POSCAR parallelepiped "
            "to fractional center (0.5, 0.5, 0.5); no PBC wrapping is used for "
            "ranking, and ties use the lowest pristine-supercell 1-based index."
        ),
        "center_fractional_coordinate": _SUPERCELL_CENTER_FRAC.tolist(),
        "ranking_distance_definition": (
            "Direct Euclidean distance in supercell Cartesian space from the "
            "candidate's in-cell fractional coordinate to frac=(0.5,0.5,0.5), "
            "without periodic wrapping"
        ),
        "periodic_diagnostic_definition": (
            "Periodic shortest-image Euclidean distance in supercell Cartesian "
            "space; recorded for diagnostics only and not used for ranking"
        ),
        "tie_break": {
            "rule": (
                "lowest pristine-supercell site_1based among candidates whose "
                "direct center distances differ from the minimum by at most the tolerance"
            ),
            "distance_tolerance_A": _CENTER_DISTANCE_TIE_TOLERANCE_A,
            "tied_candidate_sites_1based": tied_sites,
        },
        "constraints": {
            "species": expected_species,
            "parent_site_1based": parent_validation[
                "expected_parent_index_1based"
            ],
            "parent_wyckoff": parent_validation["expected_parent_wyckoff"],
            "parent_fractional_coordinate": parent_validation[
                "expected_parent_fractional_coordinate"
            ],
        },
        "candidate_count": len(candidates),
        "candidates": candidates,
        "resolved_supercell_site_1based": selected_site_1based,
        "resolved_supercell_site_0based": selected_site_1based - 1,
        "resolved_direct_distance_to_center_A": selected[
            "direct_distance_to_center_A"
        ],
        "resolved_periodic_shortest_image_distance_to_center_A": selected[
            "periodic_shortest_image_distance_to_center_A"
        ],
    }
    return selected_site_1based, selection


def _explicit_site_selection(supercell: Structure, site_1based: int) -> dict:
    """Record center diagnostics for one explicit site without enumerating peers."""
    selected = _center_distance_record(supercell, site_1based - 1)
    return {
        "mode": "explicit_site",
        "mode_defaulted": False,
        "candidate_enumeration_performed": False,
        "explicit_site_preserved": True,
        "rationale": (
            "The user supplied an exact pristine-supercell site index; only that "
            "atom was strictly validated and it was preserved without enumerating "
            "or comparing other translation-equivalent sites."
        ),
        "center_fractional_coordinate": _SUPERCELL_CENTER_FRAC.tolist(),
        "ranking_distance_definition": (
            "No ranking was performed for explicit mode; the direct in-cell "
            "Cartesian center distance is recorded as a diagnostic"
        ),
        "periodic_diagnostic_definition": (
            "Periodic shortest-image Euclidean distance in supercell Cartesian "
            "space; recorded for diagnostics only"
        ),
        "resolved_supercell_site_1based": site_1based,
        "resolved_supercell_site_0based": site_1based - 1,
        "resolved_direct_distance_to_center_A": selected[
            "direct_distance_to_center_A"
        ],
        "resolved_periodic_shortest_image_distance_to_center_A": selected[
            "periodic_shortest_image_distance_to_center_A"
        ],
        "selected_site_diagnostics": selected,
    }


def _move_target_to_host_block_front(
    structure: Structure, target_index: int, host_species: str
) -> tuple[Structure, int]:
    """Move target to the first position of one contiguous host-species block."""
    host_indices = [
        index for index, site in enumerate(structure)
        if site.species_string == host_species
    ]
    if target_index not in host_indices:
        raise ValueError(
            f"Selected target site is not in the {host_species} host-species block"
        )
    if host_indices != list(range(host_indices[0], host_indices[-1] + 1)):
        raise ValueError(
            f"Host species {host_species} is not one contiguous block; refusing to reorder"
        )

    host_start = host_indices[0]
    if target_index == host_start:
        return structure.copy(), host_start

    order = list(range(len(structure)))
    order.pop(target_index)
    order.insert(host_start, target_index)
    reordered = Structure.from_sites([structure[index] for index in order])
    return reordered, host_start


def _species_block_order(
    structure: Structure, *, reject_split_blocks: bool = False
) -> list[str]:
    """Return species blocks in file order.

    Legacy substitutions may create split blocks, so those remain reportable.
    Strict output validation rejects them when an expected order is supplied.
    """
    order: list[str] = []
    seen: set[str] = set()
    previous: str | None = None
    for site in structure:
        species = site.species_string
        if species != previous:
            if reject_split_blocks and species in seen:
                raise ValueError(f"Species {species} occurs in multiple non-contiguous blocks")
            order.append(species)
            seen.add(species)
            previous = species
    return order


def apply_substitutions(structure: Structure, substitutions: list[dict]) -> Structure:
    """Apply legacy element substitutions to a structure.

    ``sites`` remains zero-based for backward compatibility.  New strict jobs
    should use ``site_1based`` through :func:`dope_supercell` instead.
    """
    s = structure.copy()

    for sub in substitutions:
        from_el = sub["from"]
        to_el = sub["to"]
        all_sites = sub.get("all", False)
        count = sub.get("count", None)
        target_sites = sub.get("sites", None)

        from_indices = [i for i, site in enumerate(s) if site.species_string == from_el]

        if not from_indices:
            raise ValueError(f"Element '{from_el}' not found in structure")

        if target_sites:
            for idx in target_sites:
                if idx >= len(s) or s[idx].species_string != from_el:
                    print(f"Warning: site {idx} is not {from_el}, skipping")
                    continue
                s[idx] = to_el

        elif all_sites:
            for idx in from_indices:
                s[idx] = to_el

        elif count is not None:
            if count > len(from_indices):
                raise ValueError(
                    f"Cannot replace {count} {from_el}: only {len(from_indices)} available"
                )
            try:
                sga = SpacegroupAnalyzer(s, symprec=0.2)
                wyckoff_labels = sga.get_symmetry_dataset().wyckoffs
            except Exception:
                wyckoff_labels = [f"site_{i}" for i in range(len(s))]

            wyckoff_groups: dict[str, list[int]] = {}
            for idx in from_indices:
                label = wyckoff_labels[idx] if idx < len(wyckoff_labels) else f"site_{idx}"
                wyckoff_groups.setdefault(label, []).append(idx)

            picked: list[int] = []
            group_idx = 0
            group_keys = list(wyckoff_groups.keys())
            while len(picked) < count:
                label = group_keys[group_idx % len(group_keys)]
                indices = wyckoff_groups[label]
                next_i = len(picked) // len(group_keys)
                if next_i < len(indices):
                    picked.append(indices[next_i])
                group_idx += 1

            for idx in picked:
                s[idx] = to_el

        else:
            s[from_indices[0]] = to_el

    return s


def _reject_unknown_fields(obj: dict, allowed: set[str], *, context: str) -> None:
    """Fail closed on misspelled or unsupported fields in strict stage 2."""
    unknown = sorted(set(obj) - allowed)
    if unknown:
        raise ValueError(
            f"Unknown field(s) in {context}: {', '.join(unknown)}"
        )


def _validate_stage2_config_fields(config: dict) -> None:
    _reject_unknown_fields(
        config,
        {
            "supercell_manifest", "manifest", "manifest_json",
            "supercell_poscar", "supercell", "poscar", "parent_poscar",
            "substitution", "substitutions", "checks", "output",
            "provenance", "report", "report_json", "target_site",
            "parent_site", "target_site_1based", "target_site_0based",
            "expected_parent_site_1based", "expected_parent_site_0based",
            "expected_parent_wyckoff", "expected_parent_fractional_coordinate",
            "mapping_tolerance", "symprec", "site_selection",
            "mode", "supercell_matrix", "min_atoms", "min_image_distance",
        },
        context="doping config",
    )


def _validate_substitution_fields(substitution: dict) -> None:
    _reject_unknown_fields(
        substitution,
        {
            "from", "to", "site_1based", "site_0based", "parent_site",
            "target_site_1based", "target_site_0based", "site_selection",
            "sites", "count", "all",
        },
        context="substitution",
    )


def _coalesce(values: list[tuple[str, Any]], *, name: str) -> Any:
    """Return one non-None alias value and reject conflicting aliases."""
    present = [(field, value) for field, value in values if value is not None]
    if not present:
        return None
    first_field, first_value = present[0]
    for field, value in present[1:]:
        if value != first_value:
            raise ValueError(
                f"Conflicting {name} values in {first_field} and {field}"
            )
    return first_value


def _require_0based_index(value: Any, *, field: str, size: int) -> int:
    """Validate an exact integer 0-based index and return it unchanged."""
    if isinstance(value, bool) or not isinstance(value, (int, np.integer)):
        raise ValueError(f"{field} must be an integer using 0-based indexing")
    if value < 0 or value >= size:
        raise ValueError(f"{field}={value} is outside the valid 0..{size - 1} range")
    return int(value)


def _strict_target_config(
    config: dict,
    substitution: dict,
    *,
    supercell_size: int,
    parent_size: int,
) -> tuple[int | None, dict, dict]:
    """Collect an exact or center-selected target and parent constraints.

    Exact ``site_1based``/``site_0based`` values remain authoritative. If no
    exact index is supplied, ``nearest_supercell_center`` is the default.
    """
    target_raw = config.get("target_site")
    if target_raw is not None and not isinstance(target_raw, dict):
        raise ValueError("target_site must be an object")
    target_obj = target_raw or {}
    _reject_unknown_fields(
        target_obj,
        {
            "supercell_index_1based", "supercell_index_0based",
            "site_1based", "site_0based", "parent_index_1based",
            "parent_index_0based", "parent_wyckoff",
            "parent_fractional_coordinate", "mapping_tolerance", "symprec",
            "selection",
        },
        context="target_site",
    )

    parent_obj = substitution.get("parent_site") or config.get("parent_site") or {}
    if not isinstance(parent_obj, dict):
        raise ValueError("parent_site must be an object")
    _reject_unknown_fields(
        parent_obj,
        {
            "index_1based", "index_0based", "parent_index_1based",
            "parent_index_0based", "wyckoff", "fractional_coordinate",
            "fractional_coords", "mapping_tolerance", "symprec",
        },
        context="parent_site",
    )

    site_1based_raw = _coalesce(
        [
            ("substitution.site_1based", substitution.get("site_1based")),
            ("substitution.target_site_1based", substitution.get("target_site_1based")),
            ("target_site_1based", config.get("target_site_1based")),
            ("target_site.supercell_index_1based", target_obj.get("supercell_index_1based")),
            ("target_site.site_1based", target_obj.get("site_1based")),
        ],
        name="1-based supercell target site",
    )
    site_0based_raw = _coalesce(
        [
            ("substitution.site_0based", substitution.get("site_0based")),
            ("substitution.target_site_0based", substitution.get("target_site_0based")),
            ("target_site_0based", config.get("target_site_0based")),
            ("target_site.supercell_index_0based", target_obj.get("supercell_index_0based")),
            ("target_site.site_0based", target_obj.get("site_0based")),
        ],
        name="0-based supercell target site",
    )
    selection_raw = _coalesce(
        [
            ("substitution.site_selection", substitution.get("site_selection")),
            ("target_site.selection", target_obj.get("selection")),
            ("site_selection", config.get("site_selection")),
        ],
        name="supercell site selection mode",
    )
    if selection_raw is not None and selection_raw != "nearest_supercell_center":
        raise ValueError(
            "site_selection must be 'nearest_supercell_center' when provided"
        )
    if site_1based_raw is not None and site_0based_raw is not None:
        raise ValueError("Specify exactly one of site_1based or site_0based, not both")
    if selection_raw is not None and (
        site_1based_raw is not None or site_0based_raw is not None
    ):
        raise ValueError(
            "Do not combine site_selection with explicit site_1based/site_0based"
        )
    if site_1based_raw is None and site_0based_raw is None:
        site_1based = None
        index_metadata = {
            "selection_mode": "nearest_supercell_center",
            "selection_mode_defaulted": selection_raw is None,
        }
    elif site_1based_raw is not None:
        site_index = _require_1based_index(
            site_1based_raw, field="site_1based", size=supercell_size
        )
        site_1based = site_index + 1
        index_metadata = {
            "input_index_base": 1,
            "input_supercell_site_1based": site_1based,
            "input_supercell_site_0based": site_index,
            "selection_mode": "explicit_site",
            "selection_mode_defaulted": False,
        }
    else:
        site_index = _require_0based_index(
            site_0based_raw, field="site_0based", size=supercell_size
        )
        site_1based = site_index + 1
        index_metadata = {
            "input_index_base": 0,
            "input_supercell_site_1based": site_1based,
            "input_supercell_site_0based": site_index,
            "selection_mode": "explicit_site",
            "selection_mode_defaulted": False,
        }

    parent_index_1based_raw = _coalesce(
        [
            ("parent_site.index_1based", parent_obj.get("index_1based")),
            ("parent_site.parent_index_1based", parent_obj.get("parent_index_1based")),
            ("target_site.parent_index_1based", target_obj.get("parent_index_1based")),
            ("expected_parent_site_1based", config.get("expected_parent_site_1based")),
        ],
        name="1-based expected parent site",
    )
    parent_index_0based_raw = _coalesce(
        [
            ("parent_site.index_0based", parent_obj.get("index_0based")),
            ("parent_site.parent_index_0based", parent_obj.get("parent_index_0based")),
            ("target_site.parent_index_0based", target_obj.get("parent_index_0based")),
            ("expected_parent_site_0based", config.get("expected_parent_site_0based")),
        ],
        name="0-based expected parent site",
    )
    if parent_index_1based_raw is not None and parent_index_0based_raw is not None:
        raise ValueError(
            "Specify exactly one of parent index_1based or index_0based, not both"
        )
    if parent_index_1based_raw is not None:
        parent_index_1based = _require_1based_index(
            parent_index_1based_raw,
            field="parent_site.index_1based",
            size=parent_size,
        ) + 1
    elif parent_index_0based_raw is not None:
        parent_index_1based = _require_0based_index(
            parent_index_0based_raw,
            field="parent_site.index_0based",
            size=parent_size,
        ) + 1
    else:
        parent_index_1based = None

    parent_wyckoff = _coalesce(
        [
            ("parent_site.wyckoff", parent_obj.get("wyckoff")),
            ("target_site.parent_wyckoff", target_obj.get("parent_wyckoff")),
            ("expected_parent_wyckoff", config.get("expected_parent_wyckoff")),
        ],
        name="expected parent Wyckoff",
    )
    parent_fractional = _coalesce(
        [
            ("parent_site.fractional_coordinate", parent_obj.get("fractional_coordinate")),
            ("parent_site.fractional_coords", parent_obj.get("fractional_coords")),
            ("target_site.parent_fractional_coordinate", target_obj.get("parent_fractional_coordinate")),
            ("expected_parent_fractional_coordinate", config.get("expected_parent_fractional_coordinate")),
        ],
        name="expected parent fractional coordinate",
    )
    tolerance = _coalesce(
        [
            ("parent_site.mapping_tolerance", parent_obj.get("mapping_tolerance")),
            ("target_site.mapping_tolerance", target_obj.get("mapping_tolerance")),
            ("mapping_tolerance", config.get("mapping_tolerance")),
        ],
        name="mapping tolerance",
    )
    symprec = _coalesce(
        [
            ("parent_site.symprec", parent_obj.get("symprec")),
            ("target_site.symprec", target_obj.get("symprec")),
            ("symprec", config.get("symprec")),
        ],
        name="symmetry tolerance",
    )

    return site_1based, index_metadata, {
        "expected_parent_index_1based": parent_index_1based,
        "expected_parent_fractional_coordinate": parent_fractional,
        "expected_parent_wyckoff": parent_wyckoff,
        "mapping_tolerance": 1e-5 if tolerance is None else float(tolerance),
        "symprec": 0.2 if symprec is None else float(symprec),
    }


def _validation_config(config: dict) -> dict:
    checks = config.get("checks") or config.get("validation") or {}
    if not isinstance(checks, dict):
        raise ValueError("checks/validation must be an object")
    return {
        "formula": _coalesce(
            [
                ("checks.formula", checks.get("formula")),
                ("checks.expected_formula", checks.get("expected_formula")),
                ("expected_formula", config.get("expected_formula")),
                ("required_formula", config.get("required_formula")),
            ],
            name="formula check",
        ),
        "atoms": _coalesce(
            [
                ("checks.atoms", checks.get("atoms")),
                ("checks.expected_atoms", checks.get("expected_atoms")),
                ("expected_atoms", config.get("expected_atoms")),
                ("required_atoms", config.get("required_atoms")),
                ("required_atom_count", config.get("required_atom_count")),
            ],
            name="atom-count check",
        ),
        "species_order": _coalesce(
            [
                ("checks.species_order", checks.get("species_order")),
                ("checks.expected_species_order", checks.get("expected_species_order")),
                ("expected_species_order", config.get("expected_species_order")),
                ("required_species_order", config.get("required_species_order")),
            ],
            name="species-order check",
        ),
    }


def _validate_output_structure(structure: Structure, checks: dict) -> dict:
    """Run exact formula, atom-count, and species-block order checks."""
    expected_formula = checks.get("formula")
    expected_atoms = checks.get("atoms")
    expected_order = checks.get("species_order")

    actual_order = _species_block_order(
        structure, reject_split_blocks=expected_order is not None
    )
    if expected_formula is not None:
        try:
            expected_composition = Composition(expected_formula)
        except Exception as exc:
            raise ValueError(f"Invalid expected formula {expected_formula!r}") from exc
        if structure.composition != expected_composition:
            raise ValueError(
                "Output formula check failed: "
                f"actual={structure.composition.as_dict()}, "
                f"expected={expected_composition.as_dict()}"
            )
    if expected_atoms is not None:
        if isinstance(expected_atoms, bool) or not isinstance(
            expected_atoms, (int, np.integer)
        ):
            raise ValueError("Expected atom count must be an integer")
        if len(structure) != int(expected_atoms):
            raise ValueError(
                f"Output atom-count check failed: actual={len(structure)}, "
                f"expected={int(expected_atoms)}"
            )
    if expected_order is not None:
        if not isinstance(expected_order, list) or not all(
            isinstance(item, str) and item for item in expected_order
        ):
            raise ValueError("Expected species order must be a non-empty string list")
        if actual_order != expected_order:
            raise ValueError(
                f"Output species-order check failed: actual={actual_order}, "
                f"expected={expected_order}"
            )

    return {
        "formula": structure.composition.as_dict(),
        "atoms": len(structure),
        "species_order": actual_order,
        "expected_formula": expected_formula,
        "expected_atoms": expected_atoms,
        "expected_species_order": expected_order,
        "passed": True,
    }


def _sha256(path: str | Path) -> str:
    digest = hashlib.sha256()
    with open(path, "rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _write_json(path: str | Path, payload: dict) -> None:
    json_path = Path(path)
    json_path.parent.mkdir(parents=True, exist_ok=True)
    json_path.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")


def _read_json(path: str | Path) -> dict:
    with open(path, encoding="utf-8") as handle:
        payload = json.load(handle)
    if not isinstance(payload, dict):
        raise ValueError(f"JSON document must be an object: {path}")
    return payload


def _path_alias(config: dict, aliases: list[str], *, required: bool, name: str) -> str | None:
    value = _coalesce([(alias, config.get(alias)) for alias in aliases], name=name)
    if value is None and required:
        raise ValueError(f"Missing required {name}: use one of {', '.join(aliases)}")
    if value is not None and not isinstance(value, str):
        raise ValueError(f"{name} must be a path string")
    return value


def _ensure_distinct_paths(named_paths: list[tuple[str, str | Path]]) -> None:
    resolved: dict[Path, str] = {}
    for name, path in named_paths:
        candidate = Path(path).expanduser().resolve()
        if candidate in resolved:
            raise ValueError(
                f"{name} path must differ from {resolved[candidate]} path: {candidate}"
            )
        resolved[candidate] = name


def _written_poscar_validation(path: str | Path, checks: dict) -> tuple[Poscar, dict]:
    written = Poscar.from_file(str(path))
    validation = _validate_output_structure(written.structure, checks)
    if written.site_symbols != validation["species_order"]:
        raise ValueError(
            "Written POSCAR symbol line disagrees with structure block order: "
            f"symbols={written.site_symbols}, blocks={validation['species_order']}"
        )
    return written, validation


def build_supercell(config_path: str) -> dict:
    """Stage 1: generate and validate one pristine supercell plus manifest.

    Modes are mutually exclusive:
      - ``auto`` requires min_atoms + min_image_distance and forbids a matrix.
      - ``explicit`` requires supercell_matrix and forbids auto-search fields.
    """
    config = _read_json(config_path)
    mode = config.get("mode")
    if mode not in {"auto", "explicit"}:
        raise ValueError("build-supercell config mode must be 'auto' or 'explicit'")

    parent_path = _path_alias(
        config, ["parent_poscar", "poscar"], required=True, name="parent POSCAR"
    )
    assert parent_path is not None
    output = str(config.get("output", "POSCAR_pristine_supercell"))
    manifest_path = str(
        config.get("manifest")
        or config.get("manifest_json")
        or (output + ".manifest.json")
    )
    _ensure_distinct_paths(
        [("parent POSCAR", parent_path), ("pristine output", output), ("manifest", manifest_path)]
    )

    parent_poscar = Poscar.from_file(parent_path)
    parent = parent_poscar.structure
    explicit_matrix = config.get("supercell_matrix")
    if mode == "explicit":
        if explicit_matrix is None:
            raise ValueError("explicit mode requires supercell_matrix")
        forbidden = [field for field in ("min_atoms", "min_image_distance") if field in config]
        if forbidden:
            raise ValueError(
                "explicit mode forbids auto-search fields: " + ", ".join(forbidden)
            )
        supercell, matrix = generate_supercell(
            parent, supercell_matrix=explicit_matrix, allow_primitive=False
        )
        requested = {"supercell_matrix": matrix.tolist()}
    else:
        if explicit_matrix is not None:
            raise ValueError("auto mode forbids supercell_matrix")
        missing = [field for field in ("min_atoms", "min_image_distance") if field not in config]
        if missing:
            raise ValueError("auto mode requires " + " and ".join(missing))
        min_atoms = config["min_atoms"]
        min_image_distance = config["min_image_distance"]
        if isinstance(min_atoms, bool) or not isinstance(min_atoms, (int, np.integer)) or min_atoms <= 0:
            raise ValueError("min_atoms must be a positive integer")
        if not isinstance(min_image_distance, (int, float)) or isinstance(min_image_distance, bool) or min_image_distance <= 0:
            raise ValueError("min_image_distance must be a positive number")
        supercell, matrix = generate_supercell(
            parent,
            min_image_distance=float(min_image_distance),
            min_atoms=int(min_atoms),
            allow_primitive=False,
        )
        requested = {
            "min_atoms": int(min_atoms),
            "min_image_distance": float(min_image_distance),
        }

    determinant = abs(int(round(np.linalg.det(matrix))))
    expected_atoms = len(parent) * determinant
    expected_formula = (parent.composition * determinant).as_dict()
    derived_checks = {
        "formula": expected_formula,
        "atoms": expected_atoms,
        "species_order": parent_poscar.site_symbols,
    }
    declared_checks = _validation_config(config)
    checks = {
        key: declared_checks[key] if declared_checks[key] is not None else derived_checks[key]
        for key in derived_checks
    }
    _validate_output_structure(supercell, checks)

    min_image = _compute_min_image(supercell.lattice.matrix)
    if mode == "auto":
        if len(supercell) < requested["min_atoms"]:
            raise ValueError(
                f"Automatic supercell has {len(supercell)} atoms, below min_atoms={requested['min_atoms']}"
            )
        if min_image + 1e-8 < requested["min_image_distance"]:
            raise ValueError(
                "Automatic supercell minimum image distance is below request: "
                f"actual={min_image:.8g} Å, requested={requested['min_image_distance']:.8g} Å"
            )

    Path(output).parent.mkdir(parents=True, exist_ok=True)
    Poscar(supercell, sort_structure=False).write_file(output)
    written, validation = _written_poscar_validation(output, checks)
    output_abs = str(Path(output).resolve())
    parent_abs = str(Path(parent_path).resolve())
    manifest_abs = str(Path(manifest_path).resolve())
    payload = {
        "schema": "doping-supercell.supercell-manifest/v1",
        "artifact_type": "pristine_supercell",
        "status": "passed",
        "stage": "build-supercell",
        "mode": mode,
        "config": str(Path(config_path).resolve()),
        "parent": {
            "path": parent_abs,
            "sha256": _sha256(parent_path),
            "atoms": len(parent),
            "formula": parent.composition.as_dict(),
            "species_order": parent_poscar.site_symbols,
        },
        "generation": {
            "mode": mode,
            "requested": requested,
            "supercell_matrix": matrix.tolist(),
            "supercell_det": determinant,
            "minimum_image_distance_A": min_image,
            "parent_preserved": True,
            "fractional_mapping": "f_parent = f_supercell @ M mod 1",
        },
        "pristine_supercell": {
            "path": output_abs,
            "sha256": _sha256(output),
            "atoms": len(written.structure),
            "formula": written.structure.composition.as_dict(),
            "species_order": written.site_symbols,
        },
        "validation": validation,
        "manifest": manifest_abs,
        # Flat compatibility/provenance keys.
        "parent_poscar": parent_abs,
        "parent_structure_hash": _sha256(parent_path),
        "supercell_matrix": matrix.tolist(),
        "supercell_det": determinant,
        "supercell_atoms": len(written.structure),
        "supercell_min_image_distance": round(min_image, 8),
        "output": output_abs,
        "output_sha256": _sha256(output),
    }
    _write_json(manifest_path, payload)
    print(json.dumps(payload, indent=2))
    return payload


def _load_supercell_manifest(path: str | Path) -> dict:
    manifest = _read_json(path)
    if manifest.get("schema") != "doping-supercell.supercell-manifest/v1":
        raise ValueError("Unsupported or missing pristine supercell manifest schema")
    if manifest.get("artifact_type") != "pristine_supercell":
        raise ValueError("Manifest is not a pristine_supercell artifact")
    if manifest.get("status") != "passed" or not manifest.get("validation", {}).get("passed"):
        raise ValueError("Pristine supercell manifest is not accepted/passed")
    return manifest


def _single_substitution(config: dict) -> dict:
    direct = config.get("substitution")
    listed = config.get("substitutions")
    if direct is not None and listed is not None:
        raise ValueError("Specify substitution or substitutions, not both")
    if direct is not None:
        if not isinstance(direct, dict):
            raise ValueError("substitution must be an object")
        return direct
    if not isinstance(listed, list) or len(listed) != 1 or not isinstance(listed[0], dict):
        raise ValueError("Doping stage requires exactly one substitution")
    return listed[0]


def dope_existing_supercell(config_path: str) -> dict:
    """Stage 2: dope an already generated pristine supercell.

    This function never calls :func:`generate_supercell`.  It verifies the
    pristine artifact and parent hashes from the stage-1 manifest before any
    site replacement.
    """
    config = _read_json(config_path)
    forbidden = [
        field
        for field in ("mode", "supercell_matrix", "min_atoms", "min_image_distance")
        if field in config
    ]
    if forbidden:
        raise ValueError(
            "Doping stage cannot regenerate/redefine the supercell; remove: "
            + ", ".join(forbidden)
        )
    _validate_stage2_config_fields(config)

    manifest_path = _path_alias(
        config,
        ["supercell_manifest", "manifest", "manifest_json"],
        required=True,
        name="supercell manifest",
    )
    assert manifest_path is not None
    manifest = _load_supercell_manifest(manifest_path)
    manifest_hash = _sha256(manifest_path)

    pristine_record = manifest.get("pristine_supercell") or {}
    parent_record = manifest.get("parent") or {}
    matrix = _validate_supercell_matrix(
        (manifest.get("generation") or {}).get("supercell_matrix")
        or manifest.get("supercell_matrix")
    )

    configured_supercell = _path_alias(
        config,
        ["supercell_poscar", "supercell", "poscar"],
        required=False,
        name="pristine supercell POSCAR",
    )
    supercell_path = configured_supercell or pristine_record.get("path")
    if not supercell_path:
        raise ValueError("Manifest does not provide a pristine supercell path")
    configured_parent = _path_alias(
        config, ["parent_poscar"], required=False, name="parent POSCAR"
    )
    parent_path = configured_parent or parent_record.get("path")
    if not parent_path:
        raise ValueError("Manifest does not provide a parent POSCAR path")

    expected_supercell_hash = pristine_record.get("sha256") or manifest.get("output_sha256")
    expected_parent_hash = parent_record.get("sha256") or manifest.get("parent_structure_hash")
    if not expected_supercell_hash or _sha256(supercell_path) != expected_supercell_hash:
        raise ValueError("Pristine supercell SHA256 does not match the accepted manifest")
    if not expected_parent_hash or _sha256(parent_path) != expected_parent_hash:
        raise ValueError("Parent structure SHA256 does not match the accepted manifest")

    parent_poscar = Poscar.from_file(parent_path)
    pristine_poscar = Poscar.from_file(supercell_path)
    parent = parent_poscar.structure
    pristine = pristine_poscar.structure
    determinant = abs(int(round(np.linalg.det(matrix))))
    if len(pristine) != len(parent) * determinant:
        raise ValueError("Pristine supercell atom count is inconsistent with parent and matrix")
    if pristine.composition != parent.composition * determinant:
        raise ValueError("Pristine supercell formula is inconsistent with parent and matrix")
    expected_lattice = matrix @ np.asarray(parent.lattice.matrix)
    if not np.allclose(pristine.lattice.matrix, expected_lattice, atol=1e-8, rtol=0):
        raise ValueError("Pristine supercell lattice is inconsistent with parent and matrix")
    if pristine_poscar.site_symbols != pristine_record.get("species_order"):
        raise ValueError("Pristine supercell species order differs from its manifest")

    substitution = _single_substitution(config)
    _validate_substitution_fields(substitution)
    if any(key in substitution for key in ("sites", "count", "all")):
        raise ValueError(
            "Strict doping stage requires site_1based or site_0based; "
            "legacy sites/count/all are not accepted"
        )
    site_1based, index_metadata, parent_validation = _strict_target_config(
        config,
        substitution,
        supercell_size=len(pristine),
        parent_size=len(parent),
    )
    if parent_validation["expected_parent_index_1based"] is None:
        raise ValueError("Doping stage requires expected parent atom index")
    if parent_validation["expected_parent_wyckoff"] is None:
        raise ValueError("Doping stage requires expected parent Wyckoff")

    checks = _validation_config(config)
    missing_checks = [key for key in ("formula", "atoms", "species_order") if checks[key] is None]
    if missing_checks:
        raise ValueError(
            "Doping stage requires formula/atoms/species_order checks; missing "
            + ", ".join(missing_checks)
        )

    from_element = substitution["from"]
    to_element = substitution["to"]
    if site_1based is not None:
        # Exact mode intentionally validates only the requested atom. It must
        # not enumerate other candidates or acquire new failure modes from them.
        mapping = map_supercell_site_to_parent(
            parent,
            pristine,
            matrix,
            site_1based,
            expected_species=from_element,
            **parent_validation,
        )
        selection = _explicit_site_selection(pristine, site_1based)
    else:
        site_1based, selection = _automatic_site_selection(
            parent,
            pristine,
            matrix,
            selection_mode_defaulted=index_metadata["selection_mode_defaulted"],
            expected_species=from_element,
            parent_validation=parent_validation,
        )
        mapping = map_supercell_site_to_parent(
            parent,
            pristine,
            matrix,
            site_1based,
            expected_species=from_element,
            **parent_validation,
        )
    original_target_index = site_1based - 1
    doped, output_target_index = _move_target_to_host_block_front(
        pristine, original_target_index, from_element
    )
    doped[output_target_index] = to_element
    mapping.update(index_metadata)
    mapping.update(
        {
            "from": from_element,
            "to": to_element,
            "resolved_supercell_site_1based": site_1based,
            "resolved_supercell_site_0based": site_1based - 1,
            "output_site_1based": output_target_index + 1,
            "output_site_0based": output_target_index,
            "moved_to_host_block_front": output_target_index != original_target_index,
            "site_selection": selection,
        }
    )
    _validate_output_structure(doped, checks)

    output = str(config.get("output", "POSCAR_doped"))
    provenance_path = str(
        config.get("provenance")
        or config.get("report")
        or config.get("report_json")
        or (output + ".provenance.json")
    )
    _ensure_distinct_paths(
        [
            ("parent POSCAR", parent_path),
            ("pristine supercell", supercell_path),
            ("supercell manifest", manifest_path),
            ("doped output", output),
            ("doping provenance", provenance_path),
        ]
    )

    Path(output).parent.mkdir(parents=True, exist_ok=True)
    Poscar(doped, sort_structure=False).write_file(output)
    written, validation = _written_poscar_validation(output, checks)
    payload = {
        "schema": "doping-supercell.doping-provenance/v2",
        "artifact_type": "doped_supercell",
        "status": "passed",
        "stage": "dope",
        "config": str(Path(config_path).resolve()),
        "source_manifest": {
            "path": str(Path(manifest_path).resolve()),
            "sha256": manifest_hash,
            "schema": manifest["schema"],
        },
        "parent": {
            "path": str(Path(parent_path).resolve()),
            "sha256": expected_parent_hash,
            "atoms": len(parent),
            "formula": parent.composition.as_dict(),
        },
        "pristine_supercell": {
            "path": str(Path(supercell_path).resolve()),
            "sha256": expected_supercell_hash,
            "atoms": len(pristine),
            "formula": pristine.composition.as_dict(),
            "species_order": pristine_poscar.site_symbols,
        },
        "supercell_matrix": matrix.tolist(),
        "supercell_det": determinant,
        "substitution": mapping,
        "validation": validation,
        "doped_supercell": {
            "path": str(Path(output).resolve()),
            "sha256": _sha256(output),
            "atoms": len(written.structure),
            "formula": written.structure.composition.as_dict(),
            "species_order": written.site_symbols,
        },
        "output": str(Path(output).resolve()),
        "output_sha256": _sha256(output),
        "provenance": str(Path(provenance_path).resolve()),
        "report": str(Path(provenance_path).resolve()),
        "report_json": str(Path(provenance_path).resolve()),
    }
    _write_json(provenance_path, payload)
    print(json.dumps(payload, indent=2))
    return payload


def dope_supercell(config_path: str) -> dict:
    """Backward-compatible Python name for stage-2 doping only.

    Combined expansion+doping was intentionally removed.  The config must
    reference a passed stage-1 supercell manifest.
    """
    return dope_existing_supercell(config_path)
