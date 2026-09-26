---
name: defect-screening
description: "Generate and screen defect structures locally. Use before vasp-calculation creates VASP inputs."
---

# Defect Screening

This public version handles local structure generation and screening only. It
does not submit jobs or interpret an unrelaxed structure as a final
formation-energy result.

## Independent Generator

`scripts/independent_defect_generator.py` generates crystallographically
inequivalent vacancies, antisites, and interstitials from a validated host
`POSCAR`. It writes selected structures, a pristine supercell, and
`site_selection.json` with hashes, composition, lattice, candidates, and
deferral reasons.

Set `SKILL` to this skill directory and replace the placeholders locally:

```bash
SKILL=<DEFECT_SCREENING_SKILL_DIR>
python3 "$SKILL/scripts/independent_defect_generator.py" <HOST_POSCAR> \
  --output <NEW_OUTPUT_DIRECTORY> --min-image-distance 12
```

Keep machine-specific paths and runtime values in local, untracked
configuration. Explicit dopants and oxidation-state overrides must be supplied
by the user. Unknown oxidation states remain `REVIEW_REQUIRED`; the script must
not invent them.

## Additional Structure Tools

- `scripts/select_representative_sites.py` uses `doped` to build and group
  candidate sites when that package is installed.
- `scripts/screen_interstitials.py` performs geometrical checks.
- `scripts/gen_defects.py` is kept as a compatibility wrapper and does not
  create VASP inputs or submit jobs.

All generated structures must pass the current structure-validation workflow
before they are passed to `vasp-calculation`.
