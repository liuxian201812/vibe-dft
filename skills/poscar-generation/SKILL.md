---
name: poscar-generation
description: "Generate and validate POSCAR structures from experimental or template sources. Use before vasp-calculation creates VASP inputs."
---

# POSCAR Generation

This skill creates and validates structure files. It does not submit jobs,
connect to private servers, or manage task ownership.

## Workflow

1. Confirm the chemical formula, space group, lattice, atom count, and source.
2. Search the Crystallography Open Database (COD) or use an experimental
   template explicitly supplied by the user.
3. Build `POSCAR` and write `structure_manifest.json`.
4. Validate formula, site order, lattice, and source consistency.
5. Pass the validated structure to `vasp-calculation`.

Set `SKILL` to this skill directory and replace angle-bracket placeholders only
in your local command or untracked configuration:

```bash
SKILL=<POSCAR_GENERATION_SKILL_DIR>
python3 "$SKILL/scripts/search_cod.py" template "<SPACE_GROUP>" \
  <LATTICE_A> <LATTICE_B> <LATTICE_C> <ATOM_COUNT> --formula "<FORMULA>"
python3 "$SKILL/scripts/build_from_template.py" <TEMPLATE_CIF> "<FORMULA>" \
  --cations <CATION_ELEMENTS> --task-id "<TASK_ID>" \
  --source "COD:<COD_ID>" -o <OUTPUT_POSCAR>
python3 "$SKILL/scripts/validate_poscar.py" <OUTPUT_POSCAR> \
  --template <TEMPLATE_CIF> --task-id "<TASK_ID>" \
  --source "COD:<COD_ID>" --json
```

The validation status is `VALIDATED`, `REVIEW_REQUIRED`, or `UNVERIFIED`.
Do not convert a review or unverified result into a successful structure.
