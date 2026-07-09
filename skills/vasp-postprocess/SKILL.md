---
name: vasp-postprocess
description: Extract VASP results: band gap, band edges, PDOS summaries, and CCD alpha-energy curves. Use after VASP calculations finish.
---

# VASP Postprocess Basic

## Tools

- `band_structure.py`: print band gap, VBM, CBM, direct/indirect status from `vasprun.xml`
- `band_edges.py`: Python helper for VBM/CBM extraction
- `pdos.py`: Python helper for element PDOS integrals
- `ccd_diagram.py`: plot CCD image energies from a `ccd_manifest.json` directory

## Usage

```bash
python3 scripts/band_structure.py /path/to/calc_dir
python3 scripts/ccd_diagram.py /path/to/ccd --output ccd.png
```

`pdos.py` and `band_edges.py` expose importable helper functions for scripts/notebooks.
