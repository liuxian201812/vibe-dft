from pymatgen.io.vasp.outputs import Vasprun
def extract_band_edges(vr_path: str) -> dict:
    vr = Vasprun(vr_path); bs = vr.get_band_structure()
    vbm = bs.get_vbm(); cbm = bs.get_cbm()
    return {"vbm_ev": round(vbm["energy"], 4), "cbm_ev": round(cbm["energy"], 4) if cbm else None, "band_gap_ev": round(bs.get_band_gap()["energy"], 4), "is_direct": bs.get_band_gap()["direct"], "efermi_ev": round(vr.efermi, 4)}
