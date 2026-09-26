# vibe-dft — VASP workflow skill library

Public, local-first agent skills for VASP structure preparation, input
generation, submit-script assembly, defect workflows, and chemical-potential
or formation-energy post-processing.

## Quick install

```bash
git clone https://github.com/liuxian201812/vibe-dft.git
cd vibe-dft
bash install.sh              # OpenCode
bash install.sh claude-code  # Claude Code
```

Add the repository's `skills` directory to `PYTHONPATH`:

```bash
export VIBEDFT_HOME="$PWD"
export PYTHONPATH="$PYTHONPATH:$VIBEDFT_HOME/skills"
```

Set local runtime values using the placeholders documented in `.env.example`
before generating submit scripts. Do not put real cluster, account, host, path,
or credential values into tracked files.

## Included skills

| Skill | Description |
|---|---|
| `shared-utils` | POTCAR/ENCUT, DFT+U, restart helpers, and generic runtime commands |
| `poscar-generation` | Search COD or a user template, build POSCAR, and validate structure identity |
| `atom-adjust` | Translate, rotate, mirror, fix, or strain atomic structures |
| `vasp-calculation` | Generate VASP inputs and independent SLURM submit scripts |
| `defect-screening` | Generate and screen vacancy, antisite, and interstitial candidates |
| `doping-supercell` | Build doped supercells with explicit provenance and parent mapping |
| `defect-complex-builder` | Build reviewed dopant-defect complexes |
| `chemical-potential-and-formation-energy` | Convert audited JSON into stability domains and formation-energy envelopes |

## Repository layout

```text
vibe-dft/
├── README.md
├── QUICKSTART.md
├── SMOKE_TEST.md
├── install.sh
├── requirements.txt
├── .env.example
├── .github/
│   └── workflows/
│       └── smoke-test.yml
├── skills/
│   ├── shared_utils/
│   ├── poscar-generation/
│   ├── atom-adjust/
│   ├── vasp-calculation/
│   ├── defect-screening/
│   ├── doping-supercell/
│   ├── defect-complex-builder/
│   └── chemical-potential-and-formation-energy/
├── smoke_tests/
│   ├── make_test_inputs.py
│   └── run_smoke_test.sh
└── tools/
    └── check_no_site_details.sh
```

## Scope

This public repository contains local preparation and analysis code only. It
does not include private cluster profiles, remote orchestration, job ownership,
automatic server selection, or site-specific runtime facts.

The repository intentionally excludes real:

- usernames, personal identifiers, account names, and contact addresses;
- hostnames, cluster names, partitions, and network addresses;
- absolute site-specific paths and private runtime commands;
- API tokens, passwords, SSH keys, and other credentials.

Examples use descriptive placeholders. Real values must stay in the current
shell, a private configuration store, or an ignored local file.

Git commit metadata is public too. Before committing, verify that `user.name`
and `user.email` are suitable for public attribution; use a public or noreply
address rather than a private contact address.

## Quick start

```bash
python3 -m venv .venv
.venv/bin/pip install -r requirements.txt
```

Follow `QUICKSTART.md` for the local workflow. Run the local smoke test and the
privacy scan with:

```bash
PYTHON=.venv/bin/python bash smoke_tests/run_smoke_test.sh
bash tools/check_no_site_details.sh
```

The GitHub Actions workflow runs both checks on pushes and pull requests.

## Prerequisites

- Python 3.10+
- `pymatgen`, `numpy`, `seekpath`, and `matplotlib` from `requirements.txt`
- optional `doped` for the `doped`-based defect-selection workflow
- a local VASP command template in `VIBEDFT_VASP_CMD`
- a local POTCAR library path in `VIBEDFT_POTCAR_DIR`

## Before publishing changes

1. Keep all examples generic and replace environment-specific values with
   descriptive placeholders.
2. Do not add generated calculation outputs, local profiles, `.env` files, or
   credentials.
3. Check `git config user.name` and `git config user.email` for public-safe
   commit attribution.
4. Run the smoke test and privacy scan shown above.
5. Review the complete Git diff, including renamed or deleted documentation.
6. Treat broken links or references to directories that are not in the current
   public tree as documentation bugs.
