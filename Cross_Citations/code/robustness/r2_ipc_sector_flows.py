"""
==========================================================================================
 An Anatomy of the Baumol Cost Disease -- Elbadri & Vizcaino (2026)
 ------------------------------------------------------------------------------------------
 File:    robustness/r2_ipc_sector_flows.py
 Purpose: Patent flows by sector on the IPC route, 1976-2023, all patents and patents with
          US-resident inventors (input to r6)
 ------------------------------------------------------------------------------------------
 Inputs:  input/PatentsView/g_patent.tsv.zip, g_ipc_at_issue.tsv.zip
          input/PatentsView/g_inventor_disambiguated.tsv.zip   inventor variants
          input/PatentsView/g_location_disambiguated.tsv.zip  only if the inventor file
                                                               has no country column
          input/concordances/ipc4_to_isic_rev31_4.txt, isic31_4digit_to_sector.csv
 Outputs: output/robustness/ipc_sector_flows.csv   basis, year, sector, patents
 ------------------------------------------------------------------------------------------
 Notes:   Not in the paper. The share of US patents with a foreign inventor rose from about
          a tenth to over a half, so the US-inventor bases bracket the propensity correction.
          Each patent has total weight one on the `all` basis; mass is checked. The main
          output is written before the inventor variants are attempted.
==========================================================================================
"""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))   # code/, for paths.py
import paths as PATHS  # noqa: E402

import traceback
from pathlib import Path

import pandas as pd

CHUNK = 2_000_000
BLOCK = 4_000_000
YEAR_MIN = 1976
SECTORS = ["Goods", "RI_Services", "NRI_Services"]
RI_DIVISIONS = {"22", "64", "65", "66", "67", "72", "73", "74", "92"}


# ----------------------------------------------------------------- small helpers

def locate(name):
    """Input folders, as set in code/paths.py."""
    return {"PatentsView": PATHS.PATENTSVIEW, "Concordances": PATHS.CONCORDANCES,
            "tpf": PATHS.OECD_TPF, "quality": PATHS.OECD_QUALITY, "kpss": PATHS.KPSS}[name]


def pick(folder, stem):
    for ext in (".tsv.zip", ".tsv"):
        p = folder / (stem + ext)
        if p.exists():
            return p
    return None


def read_header(path):
    return list(pd.read_csv(path, sep="\t", nrows=0, dtype=str).columns)


def resolve(available, wanted, label):
    lookup = {str(c).lower(): c for c in available}
    for name in ([wanted] if isinstance(wanted, str) else wanted):
        if name.lower() in lookup:
            return lookup[name.lower()]
    sys.exit(f"No column for {label}. Tried {wanted}. File contains: {sorted(available)}")


def to_pid(series):
    """Patent identifiers as a nullable integer. Non-utility grants carry letter prefixes."""
    return pd.to_numeric(series, errors="coerce").astype("Int64")


def fallback_sector(code):
    c = str(code).zfill(4)
    if not c.isdigit():
        return None
    if int(c) < 5000:
        return "Goods"
    return "RI_Services" if c[:2] in RI_DIVISIONS else "NRI_Services"


# ----------------------------------------------------------------- concordance

def ipc_sector_table(conc):
    m = pd.read_csv(conc / "isic31_4digit_to_sector.csv", dtype=str)
    m.columns = [c.strip() for c in m.columns]
    m["isic_rev31_4"] = m["isic_rev31_4"].str.strip().str.zfill(4)
    to_sector = (m.set_index("isic_rev31_4")["Sector"].str.strip()
                  .replace({"Manufacturing": "Goods", "Other_Goods": "Goods"}))

    cw = pd.read_csv(conc / "ipc4_to_isic_rev31_4.txt", dtype=str)
    cw.columns = [c.strip() for c in cw.columns]
    c_cls = resolve(cw.columns, ["ipc4", "ipc"], "IPC4 code")
    c_isic = resolve(cw.columns, ["isic_rev31_4", "isic"], "ISIC code")
    c_w = resolve(cw.columns, ["probability_weight", "weight"], "probability weight")
    cw["cls"] = cw[c_cls].str.strip().str.strip('"').str.upper()
    cw["isic"] = cw[c_isic].str.strip().str.strip('"').str.zfill(4)
    cw["w"] = pd.to_numeric(cw[c_w], errors="coerce").fillna(0.0)
    cw["sector"] = cw["isic"].map(to_sector)
    gap = cw["sector"].isna()
    if gap.any():
        share = cw.loc[gap, "w"].sum() / cw["w"].sum()
        cw.loc[gap, "sector"] = cw.loc[gap, "isic"].map(fallback_sector)
        print(f"  {cw.loc[gap,'isic'].nunique()} ISIC codes absent from the sector map "
              f"({share:.5f} of weight), assigned by the division rule")
    cw = cw.dropna(subset=["sector"])
    S = cw.groupby(["cls", "sector"])["w"].sum().unstack(fill_value=0.0)
    for s in SECTORS:
        if s not in S.columns:
            S[s] = 0.0
    S = S[SECTORS]
    S = S.div(S.sum(axis=1), axis=0)      # row stochastic, so mass is preserved
    print(f"  {len(S):,} IPC subclasses mapped")
    return S


# ----------------------------------------------------------------- patent side

def grant_years(path):
    print(f"\nReading {path.name} ...")
    head = read_header(path)
    c_id = resolve(head, ["patent_id"], "patent id")
    c_dt = resolve(head, ["patent_date", "date", "grant_date"], "grant date")
    c_ty = resolve(head, ["patent_type", "type"], "patent type")
    pat = pd.read_csv(path, sep="\t", usecols=[c_id, c_dt, c_ty], dtype=str,
                      on_bad_lines="skip")
    pat = pat.rename(columns={c_id: "patent_id", c_dt: "date", c_ty: "type"})
    pat = pat[pat["type"].astype(str).str.lower() == "utility"]
    pat["patent_id"] = to_pid(pat["patent_id"])
    pat["year"] = pd.to_numeric(pat["date"].str[:4], errors="coerce")
    pat = pat.dropna(subset=["year", "patent_id"]).drop_duplicates("patent_id")
    pat["patent_id"] = pat["patent_id"].astype("int64")
    s = pat.set_index("patent_id")["year"].astype("int16")
    print(f"  {len(s):,} utility patents, grant years {s.min()} to {s.max()}")
    print("  grants per year, recent:")
    print(s.value_counts().sort_index().loc[2010:].to_string())
    return s


def ipc_sector_vectors(path, S):
    """One sector probability vector per patent, from its IPC subclasses."""
    print(f"\nStreaming {path.name} ...")
    head = read_header(path)
    c_id = resolve(head, ["patent_id"], "patent id")
    c_se = resolve(head, ["section"], "IPC section")
    c_cl = resolve(head, ["ipc_class", "class"], "IPC class")
    c_sb = resolve(head, ["subclass"], "IPC subclass")
    print(f"  using {c_id}, {c_se}, {c_cl}, {c_sb}")

    parts, n_rows, n_cov = [], 0, 0
    for i, ch in enumerate(pd.read_csv(path, sep="\t", usecols=[c_id, c_se, c_cl, c_sb],
                                       dtype=str, chunksize=CHUNK,
                                       on_bad_lines="skip"), 1):
        ch = ch.rename(columns={c_id: "patent_id", c_se: "section",
                                c_cl: "ipc_class", c_sb: "subclass"})
        n_rows += len(ch)
        ch = ch.dropna(subset=["patent_id", "section", "ipc_class", "subclass"])
        cls = (ch["section"].str.strip().str.upper()
               + ch["ipc_class"].str.strip().str.zfill(2)
               + ch["subclass"].str.strip().str.upper())
        pid = to_pid(ch["patent_id"])
        code = pd.Categorical(cls, categories=S.index).codes
        n_cov += int((code >= 0).sum())
        keep = (code >= 0) & pid.notna().to_numpy()
        parts.append(pd.DataFrame({"patent_id": pid.to_numpy()[keep].astype("int64"),
                                   "code": code[keep].astype("int16")}))
        print(f"  chunk {i}: {n_rows:,} rows read")

    # Deduplicate across the whole file, not within chunks, so a patent whose rows straddle
    # a chunk boundary is not counted twice.
    pairs = pd.concat(parts, ignore_index=True).drop_duplicates()
    del parts
    print(f"  {n_rows:,} rows, concordance covers {n_cov/max(n_rows,1):.4f} of them")
    print(f"  {len(pairs):,} distinct covered pairs, {pairs.patent_id.nunique():,} patents")

    pairs["w"] = 1.0 / pairs.groupby("patent_id")["code"].transform("size")
    Sm = S.to_numpy(dtype="float64")
    codes = pairs["code"].to_numpy()
    wv = pairs["w"].to_numpy()
    pid = pairs["patent_id"].to_numpy()
    frames = []
    for a in range(0, len(pairs), BLOCK):
        b = min(a + BLOCK, len(pairs))
        W = Sm[codes[a:b]] * wv[a:b, None]
        d = pd.DataFrame(W, columns=SECTORS)
        d["patent_id"] = pid[a:b]
        frames.append(d.groupby("patent_id", as_index=False)[SECTORS].sum())
    out = pd.concat(frames, ignore_index=True).groupby("patent_id",
                                                       as_index=False)[SECTORS].sum()
    total, n = out[SECTORS].to_numpy().sum(), len(out)
    print(f"  {n:,} patents retained, total weight {total:,.1f}")
    if abs(total - n) > max(1.0, 1e-9 * n):
        sys.exit(f"Mass check FAILED: {total:,.1f} against {n:,} patents.")
    print("  mass check passed")
    return out


def us_inventor_weights(pv):
    """Per patent: 1 if any inventor is in the US, and the US share of inventors."""
    f_inv = pick(pv, "g_inventor_disambiguated")
    if f_inv is None:
        print("\ng_inventor_disambiguated not found, skipping the inventor variants")
        return None
    print(f"\nStreaming {f_inv.name} ...")
    head = read_header(f_inv)
    print("  columns:", head)
    c_id = resolve(head, ["patent_id"], "patent id")
    country_col = next((c for c in head
                        if c.lower() in ("disambig_country", "country",
                                         "inventor_country", "country_code")), None)
    loc_col = next((c for c in head if c.lower() in ("location_id", "disambig_location_id")),
                   None)

    if country_col:
        print(f"  using the country column {country_col}")
        use = [c_id, country_col]
    elif loc_col:
        f_loc = pick(pv, "g_location_disambiguated")
        if f_loc is None:
            print("  no country column and no g_location_disambiguated.tsv.zip in the "
                  "PatentsView folder.\n  Download it from the PatentsView bulk page and "
                  "re-run; skipping the inventor variants.")
            return None
        print(f"  joining {f_loc.name} on {loc_col}")
        use = [c_id, loc_col]
    else:
        print("  no country or location column, skipping the inventor variants")
        return None

    parts, n_rows = [], 0
    for i, ch in enumerate(pd.read_csv(f_inv, sep="\t", usecols=use, dtype=str,
                                       chunksize=CHUNK, on_bad_lines="skip"), 1):
        n_rows += len(ch)
        ch = ch.rename(columns={c_id: "patent_id"})
        ch["patent_id"] = to_pid(ch["patent_id"])
        ch = ch.dropna(subset=["patent_id"])
        ch["patent_id"] = ch["patent_id"].astype("int64")
        parts.append(ch)
        print(f"  chunk {i}: {n_rows:,} rows read")
    inv = pd.concat(parts, ignore_index=True)
    del parts

    if country_col:
        inv["is_us"] = inv[country_col].astype(str).str.strip().str.upper().eq("US")
    else:
        loc = pd.read_csv(f_loc, sep="\t", dtype=str, on_bad_lines="skip")
        lc = resolve(loc.columns, ["location_id", "id"], "location id")
        cc = resolve(loc.columns, ["disambig_country", "country", "country_code"],
                     "country in g_location_disambiguated")
        mp = loc.set_index(lc)[cc].astype(str).str.strip().str.upper()
        inv["is_us"] = inv[loc_col].map(mp).eq("US")

    g = inv.groupby("patent_id")["is_us"].agg(["sum", "size"])
    out = pd.DataFrame({"any_us": (g["sum"] > 0).astype(float),
                        "share_us": g["sum"] / g["size"]})
    print(f"  {len(out):,} patents with inventor records, "
          f"{out.any_us.mean():.4f} have at least one US inventor")
    return out


# ----------------------------------------------------------------- main

def main():
    pv = locate("PatentsView")
    conc = locate("Concordances")
    out_dir = PATHS.OUT_ROBUSTNESS
    out_dir.mkdir(exist_ok=True, parents=True)
    print(f"PatentsView: {pv}\nConcordances: {conc}\nOutput: {out_dir}\n")

    f_pat = pick(pv, "g_patent")
    f_ipc = pick(pv, "g_ipc_at_issue")
    for f, n in ((f_pat, "g_patent"), (f_ipc, "g_ipc_at_issue")):
        if f is None:
            sys.exit(f"Not found in {pv}: {n}.tsv.zip")

    print("Building the IPC to sector table ...")
    S = ipc_sector_table(conc)

    year_of = grant_years(f_pat)
    V = ipc_sector_vectors(f_ipc, S)

    V["year"] = V["patent_id"].map(year_of)
    V = V.dropna(subset=["year"])
    V["year"] = V["year"].astype(int)
    V = V[V["year"] >= YEAR_MIN]
    print(f"\n{len(V):,} patents dated and in sample, {V.year.min()} to {V.year.max()}")

    frames = []
    g = V.groupby("year")[SECTORS].sum().reset_index().melt(
        id_vars="year", var_name="sector", value_name="patents")
    g["basis"] = "all"
    frames.append(g)

    # main output first, so a problem below does not cost the run
    pd.concat(frames)[["basis", "year", "sector", "patents"]].to_csv(
        out_dir / "ipc_sector_flows.csv", index=False)
    print(f"\nWrote {out_dir/'ipc_sector_flows.csv'} (all basis)")
    sh = g.pivot(index="year", columns="sector", values="patents")
    sh = sh.div(sh.sum(axis=1), axis=0)
    print("\nsector shares, every fifth year")
    print(sh.loc[sh.index % 5 == 0].round(4).to_string())

    try:
        W = us_inventor_weights(pv)
    except Exception:
        print("\nInventor variants failed:\n" + traceback.format_exc())
        W = None

    if W is not None:
        M = V.merge(W, left_on="patent_id", right_index=True, how="inner")
        print(f"\n{len(M):,} of {len(V):,} patents matched to inventor records")
        for name, col in (("us_any", "any_us"), ("us_share", "share_us")):
            d = M[SECTORS].mul(M[col], axis=0)
            d["year"] = M["year"].to_numpy()
            gg = (d.groupby("year")[SECTORS].sum().reset_index()
                    .melt(id_vars="year", var_name="sector", value_name="patents"))
            gg["basis"] = name
            frames.append(gg)
        pd.concat(frames)[["basis", "year", "sector", "patents"]].to_csv(
            out_dir / "ipc_sector_flows.csv", index=False)
        print(f"Rewrote {out_dir/'ipc_sector_flows.csv'} with all three bases")

        tot = (pd.concat(frames).groupby(["basis", "year"])["patents"].sum()
                 .unstack("basis"))
        print("\ntotal patents per year by basis, every fifth year")
        print(tot.loc[tot.index % 5 == 0].round(0).to_string())

    print("\nDone.")


if __name__ == "__main__":
    main()
