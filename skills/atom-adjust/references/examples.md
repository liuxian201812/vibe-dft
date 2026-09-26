# Generic command examples

All angle-bracket values are placeholders. Replace them only in your local
shell or an ignored configuration file; do not commit machine-specific paths,
account names, hostnames, or submission commands.

Set the skill directory once:

```bash
SKILL=<ATOM_ADJUST_SKILL_DIR>
```

## Inspect a translation without writing

```bash
python3 "$SKILL/scripts/adjust.py" translate \
  --atoms <ATOM_INDEXES> --dir bond:<REFERENCE_ATOM_INDEX> \
  --ratio <BOND_LENGTH_RATIO> --dry-run <SOURCE_POSCAR>
```

## Write an affine transformation

```bash
python3 "$SKILL/scripts/adjust.py" affine \
  --F <F11> <F12> <F13> <F21> <F22> <F23> <F31> <F32> <F33> \
  -o <OUTPUT_POSCAR> <SOURCE_POSCAR>
```

## Prepare a ShakeNBreak configuration

```bash
python3 "$SKILL/scripts/snb_workflow.py" prep \
  --defect <DEFECT_POSCAR> --bulk <BULK_POSCAR>
```

The helper prints third-party commands for review. It does not submit jobs or
supply private runtime settings.
