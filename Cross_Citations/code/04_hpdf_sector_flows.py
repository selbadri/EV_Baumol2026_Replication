"""
==========================================================================================
 An Anatomy of the Baumol Cost Disease -- Elbadri & Vizcaino (2026)
 ------------------------------------------------------------------------------------------
 File:    04_hpdf_sector_flows.py
 Purpose: Issued patents by sector on the USPC route, 1790-2014, from the USPTO Historical
          Patent Data Files. Used only for the change in sector composition before 1976
          (Appendix C, "Years before 1976")
 ------------------------------------------------------------------------------------------
 Inputs:  ../Empirical/input/historical_masterfile.dta   USPTO historical masterfile (not
                                                         tracked; see ../Empirical/README.md)
          ../Empirical/input/ALP_USPC/NAICS_2007/uspc_Class_to_naics07_2.txt
                                                         USPC -> NAICS 2007 crosswalk
                                                         (Goldschlag, Lybbert & Zolas 2016)
 Outputs: output/patents/hpdf_sector_flows.csv           year, Goods, RI_Services, NRI_Services
 ------------------------------------------------------------------------------------------
 Notes:   The same construction as ../Empirical/code/10_patents.do: issued patents by
          numeric USPC class, spread over two-digit NAICS sectors with the probability
          crosswalk, public administration (NAICS 92) excluded. Its 1947-2014 values
          reproduce PAT_USPC_* in data_moments.xlsx.
==========================================================================================
"""

import pandas as pd

import paths as P

SAVE_TABLES = True
CHUNK = 2_000_000
# NAICS 2007 two-digit sectors, as in Empirical/code/01_concordance_tables.do (USPTO_2digits)
RI_NAICS = {51, 52, 54}


def sector_of(naics2):
    if naics2 <= 33:
        return "Goods"
    return "RI_Services" if naics2 in RI_NAICS else "NRI_Services"


def main():
    master = P.need(P.HPDF_MASTER)
    crosswalk = P.need(P.ALP_USPC_NAICS)
    out_dir = P.OUT_PATENTS
    out_dir.mkdir(parents=True, exist_ok=True)

    # Issued patents by numeric USPC class and year. Non-numeric classes (design, plant) are
    # dropped, as `real(uspc)` does in 10_patents.do.
    parts = []
    for ch in pd.read_stata(master, columns=["uspc", "disp_dt", "disp_ty"], chunksize=CHUNK):
        ch = ch[ch["disp_ty"] == "ISS"]
        d = pd.DataFrame({"uspc": pd.to_numeric(ch["uspc"], errors="coerce"),
                          "year": ch["disp_dt"].dt.year}).dropna()
        parts.append(d.groupby(["uspc", "year"]).size())
    counts = (pd.concat(parts).groupby(level=[0, 1]).sum()
                .rename("n").reset_index())

    cw = pd.read_csv(crosswalk)
    # Private economy, as in 10_patents.do: public administration (NAICS 92) is dropped. It is
    # 5 to 16 percent of NRI-S on this route before 1976 and rising, but 0.1 percent of patents
    # on the IPC route, so keeping it made the pre-1976 change in the NRI-S composition measure
    # a component the 1976 anchor barely contains.
    cw = cw[cw["naics07_2"].astype(int) != 92]
    cw["sector"] = cw["naics07_2"].astype(int).map(sector_of)
    m = counts.merge(cw, on="uspc", how="inner")
    m["w"] = m["n"] * m["probability_weight"]
    flows = (m.groupby(["year", "sector"])["w"].sum().unstack(fill_value=0.0)
               [["Goods", "RI_Services", "NRI_Services"]])
    flows.index = flows.index.astype(int)
    flows.index.name = "year"

    print(f"{flows.index.min()} to {flows.index.max()}; "
          f"1950 shares {(flows.loc[1950] / flows.loc[1950].sum()).round(4).to_dict()}")
    if SAVE_TABLES:
        flows.to_csv(out_dir / "hpdf_sector_flows.csv")
        print(f"wrote {out_dir / 'hpdf_sector_flows.csv'}")


if __name__ == "__main__":
    main()
