"""
==========================================================================================
 An Anatomy of the Baumol Cost Disease -- Elbadri & Vizcaino (2026)
 ------------------------------------------------------------------------------------------
 File:    robustness/r5_kpss_flows.py
 Purpose: Patent flows by sector weighted by the stock-market value of each patent (Kogan,
          Papanikolaou, Seru & Stoffman 2017, extended to 2023), input to r6
 ------------------------------------------------------------------------------------------
 Inputs:  input/KPSS/KPSS_2024.csv
          input/PatentsView/g_patent.tsv.zip, g_ipc_at_issue.tsv.zip
          input/concordances/ipc4_to_isic_rev31_4.txt, isic31_4digit_to_sector.csv
          output/patents/bm_ids_cpc.csv   with --drop (step 01)
 Outputs: output/robustness/kpss_flows{tag}.csv            value-weighted counts and total value
          output/robustness/kpss_value_by_sector{tag}.csv  mean and total value
          output/robustness/kpss_coverage{tag}.csv         match rate by sector and year
 ------------------------------------------------------------------------------------------
 Notes:   Not in the paper. Usage:
            python code/robustness/r5_kpss_flows.py
            python code/robustness/r5_kpss_flows.py --drop output/patents/bm_ids_cpc.csv --tag _exbm
          KPSS values exist only for patents of listed firms, so the match rate by sector
          is reported first. Each weight is value over the grant-year mean.
==========================================================================================
"""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))   # code/, for paths.py
import paths as PATHS  # noqa: E402

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd

CHUNK, BLOCK = 2_000_000, 4_000_000
YEAR_MIN = 1976
SECTORS = ["Goods", "RI_Services", "NRI_Services"]
RI_DIVISIONS = {"22", "64", "65", "66", "67", "72", "73", "74", "92"}

_ap = argparse.ArgumentParser()
_ap.add_argument("--out")
_ap.add_argument("--drop", help="CSV with a patent_id column; those patents are dropped "
                                "before any weighting. Use output/patents/bm_ids_cpc.csv to remove "
                                "business method patents.")
_ap.add_argument("--tag", default="", help="suffix for the output file names, so an "
                                           "exclusion run does not overwrite the baseline")
ARGS, _ = _ap.parse_known_args()


def locate(name):
    """Input folders, as set in code/paths.py."""
    return {"PatentsView": PATHS.PATENTSVIEW, "Concordances": PATHS.CONCORDANCES,
            "tpf": PATHS.OECD_TPF, "quality": PATHS.OECD_QUALITY, "kpss": PATHS.KPSS}[name]


def pick(folder, stem, exts=(".tsv.zip", ".tsv", ".csv", ".csv.zip", ".txt")):
    for ext in exts:
        p = folder / (stem + ext)
        if p.exists():
            return p
    hits = [h for h in sorted(folder.rglob(f"*{stem}*"))
            if h.is_file() and h.suffix.lower() not in (".7z", ".zip") or
            h.name.lower().endswith(".csv.zip")]
    return hits[0] if hits else None


def read_header(p, sep="\t"):
    return list(pd.read_csv(p, sep=sep, nrows=0, dtype=str).columns)


def resolve(cols, wanted, label):
    lk = {str(c).strip().lower(): c for c in cols}
    for w in ([wanted] if isinstance(wanted, str) else wanted):
        if w.lower() in lk:
            return lk[w.lower()]
    sys.exit(f"No column for {label}. Tried {wanted}. File has: {sorted(cols)}")


def to_pid(s):
    return pd.to_numeric(s, errors="coerce").astype("Int64")


def fallback_sector(c):
    c = str(c).zfill(4)
    if not c.isdigit():
        return None
    if int(c) < 5000:
        return "Goods"
    return "RI_Services" if c[:2] in RI_DIVISIONS else "NRI_Services"


def ipc_sector_table(conc):
    m = pd.read_csv(conc / "isic31_4digit_to_sector.csv", dtype=str)
    m.columns = [c.strip() for c in m.columns]
    m["isic_rev31_4"] = m["isic_rev31_4"].str.strip().str.zfill(4)
    to_sec = (m.set_index("isic_rev31_4")["Sector"].str.strip()
               .replace({"Manufacturing": "Goods", "Other_Goods": "Goods"}))
    cw = pd.read_csv(conc / "ipc4_to_isic_rev31_4.txt", dtype=str)
    cw.columns = [c.strip() for c in cw.columns]
    cw["cls"] = (cw[resolve(cw.columns, ["ipc4", "ipc"], "IPC4")]
                 .str.strip().str.strip('"').str.upper())
    cw["isic"] = (cw[resolve(cw.columns, ["isic_rev31_4", "isic"], "ISIC")]
                  .str.strip().str.strip('"').str.zfill(4))
    cw["w"] = pd.to_numeric(
        cw[resolve(cw.columns, ["probability_weight", "weight"], "weight")],
        errors="coerce").fillna(0.0)
    cw["sector"] = cw["isic"].map(to_sec)
    gap = cw["sector"].isna()
    if gap.any():
        cw.loc[gap, "sector"] = cw.loc[gap, "isic"].map(fallback_sector)
    cw = cw.dropna(subset=["sector"])
    S = cw.groupby(["cls", "sector"])["w"].sum().unstack(fill_value=0.0)
    for s in SECTORS:
        if s not in S.columns:
            S[s] = 0.0
    S = S[SECTORS]
    print(f"  {len(S):,} IPC subclasses mapped")
    return S.div(S.sum(axis=1), axis=0)


def grant_years(path):
    print(f"\nReading {path.name} ...")
    h = read_header(path)
    c_id = resolve(h, ["patent_id"], "patent id")
    c_dt = resolve(h, ["patent_date", "date", "grant_date"], "grant date")
    c_ty = resolve(h, ["patent_type", "type"], "patent type")
    p = pd.read_csv(path, sep="\t", usecols=[c_id, c_dt, c_ty], dtype=str,
                    on_bad_lines="skip")
    p = p.rename(columns={c_id: "patent_id", c_dt: "date", c_ty: "type"})
    p = p[p["type"].astype(str).str.lower() == "utility"]
    p["patent_id"] = to_pid(p["patent_id"])
    p["year"] = pd.to_numeric(p["date"].str[:4], errors="coerce")
    p = p.dropna(subset=["year", "patent_id"]).drop_duplicates("patent_id")
    s = pd.Series(p["year"].astype("int16").to_numpy(),
                  index=p["patent_id"].astype("int64").to_numpy())
    print(f"  {len(s):,} utility patents, {s.min()} to {s.max()}")
    return s


def sector_vectors(path, S):
    print(f"\nStreaming {path.name} ...")
    h = read_header(path)
    want = {"patent_id": ["patent_id"], "section": ["section"],
            "ipc_class": ["ipc_class", "class"], "subclass": ["subclass"]}
    cols = {k: resolve(h, v, k) for k, v in want.items()}
    parts, n, cov = [], 0, 0
    for ch in pd.read_csv(path, sep="\t", usecols=list(cols.values()), dtype=str,
                          chunksize=CHUNK, on_bad_lines="skip"):
        ch = ch.rename(columns={v: k for k, v in cols.items()})
        n += len(ch)
        ch = ch.dropna(subset=list(want))
        cls = (ch["section"].str.strip().str.upper()
               + ch["ipc_class"].str.strip().str.zfill(2)
               + ch["subclass"].str.strip().str.upper())
        pid = to_pid(ch["patent_id"])
        code = pd.Categorical(cls, categories=S.index).codes
        cov += int((code >= 0).sum())
        keep = (code >= 0) & pid.notna().to_numpy()
        parts.append(pd.DataFrame({"patent_id": pid.to_numpy()[keep].astype("int64"),
                                   "code": code[keep].astype("int16")}))
        print(f"  {n:,} rows read")
    pairs = pd.concat(parts, ignore_index=True).drop_duplicates()
    del parts
    print(f"  concordance covers {cov / max(n, 1):.4f} of rows")
    pairs["w"] = 1.0 / pairs.groupby("patent_id")["code"].transform("size")
    Sm = S.to_numpy(dtype="float64")
    codes, wv, pid = (pairs["code"].to_numpy(), pairs["w"].to_numpy(),
                      pairs["patent_id"].to_numpy())
    frames = []
    for a in range(0, len(pairs), BLOCK):
        b = min(a + BLOCK, len(pairs))
        d = pd.DataFrame(Sm[codes[a:b]] * wv[a:b, None], columns=SECTORS)
        d["patent_id"] = pid[a:b]
        frames.append(d.groupby("patent_id", as_index=False)[SECTORS].sum())
    out = pd.concat(frames, ignore_index=True).groupby("patent_id",
                                                       as_index=False)[SECTORS].sum()
    tot, npat = out[SECTORS].to_numpy().sum(), len(out)
    print(f"  {npat:,} patents, total weight {tot:,.1f}")
    if abs(tot - npat) > max(1.0, 1e-9 * npat):
        sys.exit(f"Mass check FAILED: {tot:,.1f} against {npat:,}")
    print("  mass check passed")
    return out


def kpss_values(path):
    print(f"\nReading {path.name} ...")
    h = read_header(path, sep=",")
    print(f"  columns: {h}")
    c_id = resolve(h, ["patent_num", "patnum", "patent_id"], "patent number")
    c_v = resolve(h, ["xi_real"], "real patent value")
    c_n = next((c for c in h if str(c).strip().lower() == "xi_nominal"), None)
    use = [c_id, c_v] + ([c_n] if c_n else [])
    k = pd.read_csv(path, usecols=use, on_bad_lines="skip")
    k = k.rename(columns={c_id: "patent_id", c_v: "xi_real",
                          **({c_n: "xi_nominal"} if c_n else {})})
    k["patent_id"] = to_pid(k["patent_id"])
    k["xi_real"] = pd.to_numeric(k["xi_real"], errors="coerce")
    k = k.dropna(subset=["patent_id", "xi_real"])
    k["patent_id"] = k["patent_id"].astype("int64")
    before = len(k)
    k = k.drop_duplicates(subset=["patent_id"])
    print(f"  {before:,} valued patents, {len(k):,} distinct")
    print(f"  xi_real: mean {k.xi_real.mean():.2f}, median {k.xi_real.median():.2f}, "
          f"max {k.xi_real.max():,.0f} (millions of 1982 dollars)")
    return k


def main():
    pv, conc, kp = locate("PatentsView"), locate("Concordances"), locate("kpss")
    out_dir = Path(ARGS.out) if ARGS.out else PATHS.OUT_ROBUSTNESS
    out_dir.mkdir(exist_ok=True, parents=True)
    print(f"PatentsView:  {pv}\nConcordances: {conc}\nKPSS:         {kp}\n"
          f"Output:       {out_dir}\n")

    f_pat, f_ipc = pick(pv, "g_patent"), pick(pv, "g_ipc_at_issue")
    f_k = pick(kp, "KPSS", exts=(".csv",))
    for f, nm in ((f_pat, "g_patent"), (f_ipc, "g_ipc_at_issue"), (f_k, "KPSS_2024.csv")):
        if f is None:
            sys.exit(f"Not found: {nm}")

    print("Building the IPC to sector table ...")
    S = ipc_sector_table(conc)
    year_of = grant_years(f_pat)
    V = sector_vectors(f_ipc, S)
    V["year"] = V["patent_id"].map(year_of)
    V = V.dropna(subset=["year"])
    V["year"] = V["year"].astype(int)
    V = V[V["year"] >= YEAR_MIN].reset_index(drop=True)

    if ARGS.drop:
        d = pd.read_csv(ARGS.drop)
        col = "patent_id" if "patent_id" in d.columns else d.columns[0]
        ids = set(pd.to_numeric(d[col], errors="coerce").dropna().astype("int64"))
        before = len(V)
        V = V[~V["patent_id"].isin(ids)].reset_index(drop=True)
        print(f"\nDropped {before - len(V):,} of {before:,} patents listed in "
              f"{Path(ARGS.drop).name}")

    K = kpss_values(f_k)
    M = V.merge(K[["patent_id", "xi_real"]], on="patent_id", how="left")
    has = M["xi_real"].notna()
    print(f"\n{has.mean():.4f} of patents in sample carry a KPSS value")

    # ---------------------------------------------------------------- coverage first
    num = M[SECTORS].mul(has.astype(float), axis=0)
    den = M[SECTORS]
    num["year"] = M["year"].to_numpy(); den["year"] = M["year"].to_numpy()
    cov = (num.groupby("year")[SECTORS].sum() / den.groupby("year")[SECTORS].sum())
    cov.to_csv(out_dir / f"kpss_coverage{ARGS.tag}.csv")
    print("\nShare of each sector's patents carrying a KPSS value, every fifth year.")
    print("Read this before anything else: KPSS covers public firms only, and if the")
    print("sectors are matched at very different rates the comparison is selected.\n")
    print(cov.loc[cov.index % 5 == 0].round(4).to_string())

    # ---------------------------------------------------------------- value by sector
    D = M[has].copy()
    vnum = D[SECTORS].mul(D["xi_real"].to_numpy(), axis=0)
    vden = D[SECTORS]
    vnum["year"] = D["year"].to_numpy(); vden["year"] = D["year"].to_numpy()
    total_v = vnum.groupby("year")[SECTORS].sum()
    mean_v = total_v / vden.groupby("year")[SECTORS].sum()
    pd.concat({"mean_value": mean_v, "total_value": total_v}, axis=1).to_csv(
        out_dir / f"kpss_value_by_sector{ARGS.tag}.csv")

    print("\n\nMean value per patent by sector, millions of 1982 dollars")
    p = mean_v.copy(); p["RI_over_Goods"] = p.RI_Services / p.Goods
    print(p.loc[p.index % 5 == 0].round(3).to_string())
    print("\nTotal value by sector, and the RI-Services share of it")
    q = total_v.copy(); q["RI_share"] = q.RI_Services / q.sum(axis=1)
    print(q.loc[q.index % 5 == 0].round(0).to_string())

    # ---------------------------------------------------------------- flows
    mu = D.groupby("year")["xi_real"].transform("mean")
    wt = np.where(mu > 0, D["xi_real"] / mu, 0.0)
    frames = []
    g = M.groupby("year")[SECTORS].sum().reset_index().melt(
        id_vars="year", var_name="sector", value_name="patents")
    g["measure"] = "unweighted_all"
    frames.append(g)
    g = D.groupby("year")[SECTORS].sum().reset_index().melt(
        id_vars="year", var_name="sector", value_name="patents")
    g["measure"] = "unweighted_matched"
    frames.append(g)
    e = D[SECTORS].mul(wt, axis=0); e["year"] = D["year"].to_numpy()
    g = (e.groupby("year")[SECTORS].sum().reset_index()
           .melt(id_vars="year", var_name="sector", value_name="patents"))
    g["measure"] = "value_weighted"
    frames.append(g)
    g = (total_v.reset_index().melt(id_vars="year", var_name="sector",
                                    value_name="patents"))
    g["measure"] = "total_value"
    frames.append(g)
    out = pd.concat(frames, ignore_index=True)[["measure", "year", "sector", "patents"]]
    out.to_csv(out_dir / f"kpss_flows{ARGS.tag}.csv", index=False)

    print(f"\nWrote {out_dir / f'kpss_coverage{ARGS.tag}.csv'}, "
          f"{out_dir / f'kpss_value_by_sector{ARGS.tag}.csv'}, "
          f"{out_dir / f'kpss_flows{ARGS.tag}.csv'}")

    print("\nRI-Services share of the flow, by measure, every fifth year")
    piv = {}
    for m in out.measure.unique():
        pp = out[out.measure == m].pivot(index="year", columns="sector", values="patents")
        piv[m] = pp.RI_Services / pp.sum(axis=1)
    P = pd.DataFrame(piv)
    print(P.loc[P.index % 5 == 0].round(4).to_string())


if __name__ == "__main__":
    main()
