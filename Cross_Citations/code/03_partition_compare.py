"""
==========================================================================================
 An Anatomy of the Baumol Cost Disease -- Elbadri & Vizcaino (2026)
 ------------------------------------------------------------------------------------------
 File:    03_partition_compare.py
 Purpose: Check that the patent-side sector map (ISIC Rev 3.1 classes) and the NAICS
          partition used for R&D workers (Appendix A) agree, class by class (Appendix D,
          notes to Table D.1)
 ------------------------------------------------------------------------------------------
 Inputs:  input/concordances/isic31_4digit_to_sector.csv   ISIC Rev 3.1 -> sector (Table D.1)
          input/concordances/ipc4_to_isic_rev31_4.txt      ALP concordance, for the weights
 Outputs: output/patents/partition_comparison.csv          every service ISIC class, its NAICS
                                                           code, and both sector labels
          output/patents/partition_mismatches.tex          the disagreeing classes (7494, 9231)
          output/patents/isic31_4digit_to_sector_harmonised.csv
                                                           the map with them moved to RI-Services,
                                                           used by the 'harmonised' run of step 05
 ------------------------------------------------------------------------------------------
 Notes:   RI-Services is NAICS 51, 52 and 54 (Appendix A). The NAICS counterpart of each
          service ISIC class is recorded in NAICS_OF below, from the ISIC Rev 3.1 to NAICS
          2002 correspondence.
==========================================================================================
"""

import sys

import pandas as pd

import paths as P

SAVE_TABLES = True

SECTORS = ["Goods", "RI_Services", "NRI_Services"]
RI_NAICS2 = {"51", "52", "54"}          # Appendix D, Research Intensive Services

# ISIC Rev 3.1 four digit class -> NAICS code of the BEA industry it belongs to.
# Service classes only. Goods classes are ISIC below 5000 on both sides by construction.
NAICS_OF = {
    "2211": "5111", "2212": "5111", "2213": "5111", "2219": "5112",
    "5010": "441", "5020": "8111", "5030": "441", "5040": "441", "5050": "447",
    "5110": "425", "5122": "424", "5131": "424", "5139": "423", "5141": "423",
    "5142": "423", "5143": "423", "5151": "423", "5152": "423", "5159": "423",
    "5211": "445", "5220": "445", "5231": "446", "5233": "446", "5239": "453",
    "5259": "453", "5260": "8114",
    "5510": "721", "5520": "722",
    "6010": "482", "6021": "485", "6023": "485", "6030": "486", "6110": "483",
    "6210": "481", "6220": "481", "6304": "5615", "6309": "488",
    "6411": "491", "6412": "492", "6420": "517",
    "6511": "522", "6519": "522", "6592": "522", "6599": "523",
    "6601": "524", "6602": "525", "6603": "524",
    "6711": "523", "6712": "523", "6719": "523", "6720": "524",
    "7010": "531", "7020": "532",
    "7221": "5112", "7229": "5415", "7230": "5182", "7240": "5191", "7250": "8112",
    "7290": "5415",
    "7310": "5417",
    "7411": "5411", "7412": "5412", "7413": "5413", "7414": "5416",
    "7421": "5413", "7422": "5413", "7430": "5418",
    "7491": "5613", "7492": "5616", "7493": "5617", "7494": "5419",
    "7495": "5619", "7499": "5619",
    "7511": "921", "7512": "923", "7513": "926", "7514": "921",
    "7521": "928", "7523": "922", "7530": "923",
    "8010": "611", "8022": "611", "8030": "611", "8090": "611",
    "8511": "622", "8512": "621", "8519": "621", "8531": "623", "8532": "624",
    "9000": "562",
    "9211": "5121", "9212": "5121", "9213": "5151", "9214": "711", "9219": "711",
    "9220": "5191", "9231": "5191", "9232": "712", "9233": "712",
    "9241": "713", "9249": "713",
    "9302": "8121", "9303": "8122", "9309": "8129", "9500": "814", "9900": "928",
}

NAMES = {
    "7494": "Photographic activities",
    "9231": "Library and archive activities",
    "7414": "Business and management consultancy",
    "9220": "News agency activities",
    "6420": "Telecommunications",
    "7250": "Maintenance and repair of office and computing machinery",
}


def main():
    f_map = P.need(P.SECTOR_MAP)
    f_cw = P.need(P.ALP_IPC)
    out_dir = P.OUT_PATENTS
    out_dir.mkdir(parents=True, exist_ok=True)
    print(f"Sector map:  {f_map}\nConcordance: {f_cw}\nOutput:      {out_dir}\n")

    m = pd.read_csv(f_map, dtype=str)
    m.columns = [c.strip() for c in m.columns]
    m["isic"] = m["isic_rev31_4"].str.strip().str.zfill(4)
    m["patent_sector"] = (m["Sector"].str.strip()
                          .replace({"Manufacturing": "Goods", "Other_Goods": "Goods"}))
    svc = m[m["patent_sector"] != "Goods"].copy()

    svc["naics"] = svc["isic"].map(NAICS_OF)
    unmapped = svc[svc["naics"].isna()]
    if len(unmapped):
        print("Service ISIC classes with no NAICS entry in this script:")
        print(unmapped[["isic", "patent_sector"]].to_string(index=False))
        sys.exit("Add them to NAICS_OF before drawing conclusions.")

    svc["naics2"] = svc["naics"].str[:2]
    svc["allocation_sector"] = svc["naics2"].apply(
        lambda n: "RI_Services" if n in RI_NAICS2 else "NRI_Services")
    svc["agree"] = svc["patent_sector"] == svc["allocation_sector"]

    # patent-side weight reaching each ISIC class, from the concordance alone
    cw = pd.read_csv(f_cw, dtype=str)
    cw.columns = [c.strip() for c in cw.columns]
    cw["isic"] = cw["isic_rev31_4"].str.strip().str.zfill(4)
    cw["w"] = pd.to_numeric(cw["probability_weight"], errors="coerce").fillna(0.0)
    w_isic = cw.groupby("isic")["w"].sum()
    n_ipc = cw.groupby("isic")["ipc4"].nunique()
    svc["conc_weight"] = svc["isic"].map(w_isic).fillna(0.0)
    svc["n_ipc4"] = svc["isic"].map(n_ipc).fillna(0).astype(int)
    svc["conc_share"] = svc["conc_weight"] / cw["w"].sum()
    svc["name"] = svc["isic"].map(NAMES)

    cols = ["isic", "name", "naics", "patent_sector", "allocation_sector", "agree",
            "n_ipc4", "conc_weight", "conc_share"]
    out = svc[cols].sort_values("isic")
    if SAVE_TABLES:
        out.to_csv(out_dir / "partition_comparison.csv", index=False)

    bad = out[~out["agree"]]
    print(f"{len(out)} service ISIC classes compared, {len(bad)} disagree.\n")
    print("Patent-side RI_Services classes, and where NAICS puts them:")
    print(out[out.patent_sector == "RI_Services"]
          [["isic", "naics", "allocation_sector", "conc_share"]].to_string(index=False))
    print("\nDisagreements:")
    print(bad[["isic", "name", "naics", "patent_sector", "allocation_sector",
               "n_ipc4", "conc_share"]].to_string(index=False))
    print(f"\nConcordance weight in disagreeing classes: {bad.conc_share.sum():.5f} "
          f"of the total, in {bad.n_ipc4.sum()} IPC subclasses.")
    if len(bad):
        top = cw[cw["isic"].isin(bad["isic"])].sort_values("w", ascending=False)
        print("\nIPC subclasses feeding them:")
        print(top[["ipc4", "isic", "w"]].head(15).to_string(index=False))

    # ------------------------------------------------------------ harmonised map
    h = m[["isic_rev31_4", "Sector", "Note"]].copy()
    fix = dict(zip(bad["isic"], bad["allocation_sector"]))
    idx = m["isic"].isin(fix)
    h.loc[idx, "Sector"] = m.loc[idx, "isic"].map(fix)
    h.loc[idx, "Note"] = ("reassigned to match NAICS "
                          + m.loc[idx, "isic"].map(dict(zip(bad["isic"], bad["naics"]))))
    if SAVE_TABLES:
        h.to_csv(out_dir / "isic31_4digit_to_sector_harmonised.csv", index=False)
        with open(out_dir / "partition_mismatches.tex", "w") as fh:
            fh.write("\\begin{tabular}{llllr}\n\\toprule\n"
                     "ISIC & Activity & NAICS & Patent side & Allocation side \\\\\n"
                     "\\midrule\n")
            for _, r in bad.iterrows():
                fh.write(f"{r.isic} & {r['name'] or ''} & {r.naics} & "
                         f"{r.patent_sector.replace('_','-')} & "
                         f"{r.allocation_sector.replace('_','-')} \\\\\n")
            fh.write("\\bottomrule\n\\end{tabular}\n")
    print(f"\nWrote partition_comparison.csv, partition_mismatches.tex and "
          f"isic31_4digit_to_sector_harmonised.csv to {out_dir}")


if __name__ == "__main__":
    main()
