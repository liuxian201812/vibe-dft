# Quick start

This walkthrough keeps runtime details local and uses placeholders for every
site-specific value.

## 1. Clone and install

```bash
git clone https://github.com/liuxian201812/vibe-dft.git
cd vibe-dft
python3 -m venv .venv
.venv/bin/pip install -r requirements.txt
export PYTHONPATH="$PWD/skills:$PYTHONPATH"
```

## 2. Configure local runtime values

Copy the variable names from `.env.example` into your current shell or a private,
untracked shell profile. Replace the placeholders locally:

```bash
export VIBEDFT_VASP_CMD='<VASP_LAUNCH_COMMAND_WITH_{vasp_bin}>'
export VIBEDFT_POTCAR_DIR='<POTCAR_LIBRARY_PATH>'
```

`VIBEDFT_VASP_CMD` must contain the literal `{vasp_bin}` placeholder so the
scripts can substitute `vasp_std`, `vasp_gam`, or another executable name.

Never commit real usernames, account names, contact addresses, hostnames,
cluster names, network addresses, scheduler values, absolute site paths,
tokens, or keys. Scheduler partition, account, QoS, walltime, and work
directories belong in a local task configuration or an ignored private profile.

## 3. Validate a structure

```bash
python3 skills/poscar-generation/scripts/validate_poscar.py POSCAR \
  --task-id example --source local
```

For strict space-group comparison, also pass `--template <TEMPLATE_CIF>`.

## 4. Generate VASP inputs

Create a local `relax_config.json`. Values in angle brackets are placeholders
and must be replaced before execution:

```json
{
  "workdir": "<LOCAL_WORK_DIRECTORY>",
  "runtime": "<LOCAL_RUNTIME_PROFILE_NAME>",
  "partition": "<SCHEDULER_PARTITION>",
  "ncpus": 32,
  "host_structures": {"host": "POSCAR"},
  "tasks": [
    {
      "id": "01",
      "name": "relax",
      "type": "relax",
      "host": "host",
      "template": "LR",
      "functional": "PBE",
      "isif": 3
    }
  ]
}
```

Keep this file untracked when it contains real site-specific values. Generate
the input files with:

```bash
python3 skills/vasp-calculation/scripts/vasp_calculation.py prepare -- \
  --phase relax --config relax_config.json
```

## 5. Generate and inspect submit scripts

```bash
python3 skills/vasp-calculation/scripts/vasp_calculation.py submit-scripts -- \
  --config relax_config.json
```

The generator writes one independent job script per step and a top-level
`submit.sh`. It does not create dependency chains, select a remote server, or
track remote jobs. Review every generated script before submitting it.

## 6. Defect and formation-energy workflows

- Use `atom-adjust` for explicit atomic translation, rotation, mirroring,
  selective-dynamics changes, or lattice strain.
- Use `defect-screening` to generate and screen candidate defect structures.
- Use `doping-supercell` or `defect-complex-builder` for reviewed doped or
  complex-defect structures.
- Pass validated structures to `vasp-calculation`.
- Use `chemical-potential-and-formation-energy` only with audited JSON; it does
  not extract numbers from raw VASP output.

## 7. Verify before publishing

Check that commit attribution is public-safe:

```bash
git config user.name
git config user.email
```

Use a public or noreply address rather than a private contact address. Then run
both repository checks from the root:

```bash
bash tools/check_no_site_details.sh
PYTHON=.venv/bin/python bash smoke_tests/run_smoke_test.sh
```

The privacy scanner reports only the file path and rule name; it does not echo
the matched secret into CI logs. It checks the current working tree, not old Git
history.

Delete generated local test files when they are no longer needed:

```bash
rm -rf smoke_tests/tmp
```
