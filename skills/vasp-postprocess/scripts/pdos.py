from pymatgen.io.vasp.outputs import Vasprun


def extract_pdos(vr_path: str, elements: list = None) -> dict:
    import numpy as np

    vr = Vasprun(vr_path); cdos = vr.complete_dos
    energies = np.array(cdos.energies, dtype=float)
    totals = {}
    pdos = cdos.get_element_dos()
    for e in cdos.structure.composition.elements:
        name = str(e)
        if elements and name not in elements: continue
        if name in pdos:
            densities = np.array(pdos[name].get_densities(), dtype=float)
            totals[name] = {
                "integrated_states": round(float(np.trapz(densities, energies)), 4),
                "integration": "trapz DOS(E) dE",
            }
    return {"elements": list(totals.keys()), "element_pdos_integrals": totals, "efermi": round(cdos.efermi, 4), "band_gap": round(cdos.get_gap(), 4) if cdos.get_gap() else None}
