---
name: poscar-generation
description: Generate host POSCARs from COD CIFs and validate POSCAR format/geometry. Use when preparing pristine structures for VASP workflows without embedding site-specific runtime facts.
---

# POSCAR Generation Basic

## Scope

This skill only supports direct COD search/download and POSCAR validation.

It does not do template substitution, Materials Project fallback, literature structure recovery, or remote relaxation submission.

## Commands

```bash
python3 scripts/search_cod.py search "<formula>"
python3 scripts/search_cod.py download <cod_id> -o structure.cif
python3 scripts/validate_poscar.py POSCAR
```

Convert CIF to POSCAR with pymatgen:

```python
from pymatgen.core import Structure
Structure.from_file("structure.cif").to(filename="POSCAR", fmt="poscar")
```

## Validation

`validate_poscar.py` checks lattice sanity, atom overlap, sampled bond lengths, optional density/space-group comparison to a template.
