# vibedft — VASP workflow skill library

Agent skills for structure preparation, VASP input generation, SLURM submit-script assembly, and lightweight post-processing.

## Quick Install

```bash
git clone https://github.com/liuxian201812/vibe-dft.git
cd vibe-dft
bash install.sh              # OpenCode (default)
bash install.sh claude-code  # Claude Code
```

Add to your shell profile:

```bash
export VIBEDFT_HOME="$PWD"
export PYTHONPATH="$PYTHONPATH:$VIBEDFT_HOME"
```

Set runtime variables (see `.env.example`) and you are ready.

## Skills

| Skill | Description |
|-------|-------------|
| `shared-utils` | Common library: POTCAR/ENCUT, DFT+U auto-detection, restart loops, runtime command retrieval |
| `poscar-generation` | Search COD CIF by formula, download, validate POSCAR format/geometry |
| `vasp-input` | Generate INCAR/KPOINTS/POTCAR/POSCAR/submit.sh for relax, static, band, CCD |
| `vasp-submit` | Assemble independent SLURM submit scripts from task manifests |
| `ccd` | CCD wrappers: interpolate image inputs from two endpoint structures, plot alpha-energy curve |
| `vasp-postprocess` | Extract band gap, band edges, PDOS, CCD diagram after VASP finishes |
| `runtime-config` | Guidance for supplying site-specific runtime details without leaking secrets |

## Quick Start

```bash
python3 -m venv .venv
.venv/bin/pip install -r requirements.txt
```

Follow `QUICKSTART.md` for a walkthrough or run `bash smoke_tests/run_smoke_test.sh` to verify the toolchain.

## Repository Layout

```
vibedft/
├── install.sh               # one-command skill installer
├── shared_utils/            # Python library (add to PYTHONPATH)
├── skills/
│   ├── poscar-generation/
│   ├── vasp-input/
│   ├── vasp-submit/
│   ├── ccd/
│   ├── vasp-postprocess/
│   └── runtime-config/
├── smoke_tests/
├── .env.example
└── requirements.txt
```

## Prerequisites

- Python 3.10+
- pymatgen, numpy, seekpath, matplotlib (`requirements.txt`)
- A VASP command template in `VIBEDFT_VASP_CMD`
- A POTCAR library path in `VIBEDFT_POTCAR_DIR`
- (optional) A SLURM cluster for remote execution

## What's Included

- COD CIF search, download, and POSCAR validation
- Relax/static/scf, band, and CCD input generation
- Independent SLURM submit scripts (no `afterok` chains)
- Band gap, PDOS, and CCD diagram post-processing
- Generic runtime configuration (no site-specific secrets tracked)

## Smoke Test

```bash
PYTHON=.venv/bin/python bash smoke_tests/run_smoke_test.sh
```

Verifies POSCAR validation, relax/band/CCD input generation, submit scripts, and CCD plotting — all locally without a cluster.
