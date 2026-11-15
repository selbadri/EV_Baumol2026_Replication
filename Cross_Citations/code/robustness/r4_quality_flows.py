"""
==========================================================================================
 An Anatomy of the Baumol Cost Disease -- Elbadri & Vizcaino (2026)
 ------------------------------------------------------------------------------------------
 File:    robustness/r4_quality_flows.py
 Purpose: Patent flows by sector weighted by the OECD patent quality indicators, input to r6
 ------------------------------------------------------------------------------------------
 Inputs:  input/OECD_Quality_Feb2025/202502_OECD_PATENT_QUALITY_USPTO_INDIC.txt
          input/PatentsView/g_patent.tsv.zip, g_ipc_at_issue.tsv.zip
          input/concordances/ipc4_to_isic_rev31_4.txt, isic31_4digit_to_sector.csv
 Outputs: output/robustness/quality_sector_flows.csv   flows by indicator, two normalisations
          output/robustness/quality_by_sector.csv      mean indicators by sector and year
          output/robustness/quality_coverage.csv       match rate by year
 ------------------------------------------------------------------------------------------
 Notes:   Not in the paper. Uses the raw indicators (fwd_cits5, fwd_cits7, family_size,
          renewal), each weight (1 + Q) over its grant-year mean so weights average one.
          The OECD composites are normalised within technology field, which removes the
          cross-sector variation of interest, so they are reported only as a benchmark.
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
RAW = ["fwd_cits5", "fwd_cits7", "family_size", "renewal"]
COMPOSITE = ["quality_index_4", "quality_index_6"]

_ap = argparse.ArgumentParser()
_ap.add_argument("--out")
ARGS, _ = _ap.parse_known_args()


def locate(name):
    """Input folders, as set in code/paths.py."""
    return {"PatentsView": PATHS.PATENTSVIEW, "Concordances": PATHS.CONCORDANCES,
            "tpf": PATHS.OECD_TPF, "quality": PATHS.OECD_QUALITY, "kpss": PATHS.KPSS}[name]


def pick(folder, stem, exts=(".tsv.zip", ".tsv", ".txt")):
    for ext in exts:
        p = folder / (stem + ext)
        if p.exists():
            return p
    hits = [h for h in sorted(folder.rglob(f"*{stem}*"))
            if h.is_file() and h.suffix.lower() != ".7z"]
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


def quality_table(path):
    """Per patent quality indicators, keyed on the US grant number in pub_nbr."""
    print(f"\nStreaming {path.name} ...")
    h = read_header(path, sep="|")
    print(f"  columns: {h}")
    c_pub = resolve(h, ["pub_nbr", "publn_nbr", "us_pub_nbr"], "US publication number")
    c_fil = resolve(h, ["filing"], "filing year")
    c_tf = resolve(h, ["tech_field"], "technology field")
    have = [c for c in RAW + COMPOSITE if c in h]
    missing = [c for c in RAW + COMPOSITE if c not in h]
    if missing:
        print(f"  not present, skipped: {missing}")
    use = [c_pub, c_fil, c_tf] + have

    parts, n = [], 0
    for ch in pd.read_csv(path, sep="|", usecols=use, dtype=str, chunksize=CHUNK,
                          on_bad_lines="skip"):
        n += len(ch)
        ch = ch.rename(columns={c_pub: "pub", c_fil: "filing", c_tf: "tech_field"})
        pid = to_pid(ch["pub"].str.strip().str.upper().str.removeprefix("US"))
        d = pd.DataFrame({"patent_id": pid})
        d["filing"] = pd.to_numeric(ch["filing"], errors="coerce")
        d["tech_field"] = pd.to_numeric(ch["tech_field"], errors="coerce")
        for c in have:
            d[c] = pd.to_numeric(ch[c], errors="coerce")
        parts.append(d.dropna(subset=["patent_id"]))
        print(f"  {n:,} rows read")
    q = pd.concat(parts, ignore_index=True)
    del parts
    q["patent_id"] = q["patent_id"].astype("int64")
    before = len(q)
    q = q.drop_duplicates(subset=["patent_id"])
    print(f"  {n:,} rows, {before:,} with a numeric grant number, "
          f"{len(q):,} distinct patents")
    return q, have


def main():
    pv, conc, qual = locate("PatentsView"), locate("Concordances"), locate("quality")
    out_dir = Path(ARGS.out) if ARGS.out else PATHS.OUT_ROBUSTNESS
    out_dir.mkdir(exist_ok=True, parents=True)
    print(f"PatentsView:  {pv}\nConcordances: {conc}\nQuality:      {qual}\n"
          f"Output:       {out_dir}\n")

    f_pat, f_ipc = pick(pv, "g_patent"), pick(pv, "g_ipc_at_issue")
    f_q = pick(qual, "202502_OECD_PATENT_QUALITY_USPTO_INDIC", exts=(".txt",))
    if f_q is not None and "COHORT" in f_q.name.upper():
        cands = [h for h in qual.rglob("*INDIC*")
                 if h.is_file() and "COHORT" not in h.name.upper()
                 and h.suffix.lower() != ".7z"]
        f_q = cands[0] if cands else None
    for f, nm in ((f_pat, "g_patent"), (f_ipc, "g_ipc_at_issue"),
                  (f_q, "the patent-level quality file")):
        if f is None:
            sys.exit(f"Not found: {nm}")

    print("Building the IPC to sector table ...")
    S = ipc_sector_table(conc)
    year_of = grant_years(f_pat)
    V = sector_vectors(f_ipc, S)
    V["year"] = V["patent_id"].map(year_of)
    V = V.dropna(subset=["year"])
    V["year"] = V["year"].astype(int)
    V = V[V["year"] >= YEAR_MIN]

    Q, have = quality_table(f_q)
    M = V.merge(Q.drop(columns=["filing", "tech_field"]), on="patent_id", how="left")
    matched = M[have[0]].notna() if have else pd.Series(False, index=M.index)
    cov = (pd.DataFrame({"year": M["year"], "m": matched.astype(float)})
             .groupby("year")["m"].agg(["mean", "size"])
             .rename(columns={"mean": "matched", "size": "patents"}))
    cov.to_csv(out_dir / "quality_coverage.csv")
    print("\nMatch rate to the OECD quality file, every fifth year")
    print(cov.loc[cov.index % 5 == 0].round(4).to_string())

    print("\nMean of each indicator by grant year, every fifth year")
    print("(watch fwd_cits5 and fwd_cits7 collapse at the end: that is window truncation)")
    print(M.groupby("year")[have].mean().pipe(
        lambda d: d.loc[d.index % 5 == 0]).round(3).to_string())

    # ------------------------------------------------------------- key diagnostic
    rows = []
    for c in have:
        d = M.dropna(subset=[c])
        w = d[SECTORS].to_numpy()                      # sector probability per patent
        num = pd.DataFrame(w * d[c].to_numpy()[:, None], columns=SECTORS)
        den = pd.DataFrame(w, columns=SECTORS)
        num["year"] = d["year"].to_numpy(); den["year"] = d["year"].to_numpy()
        mean_by = (num.groupby("year")[SECTORS].sum()
                   / den.groupby("year")[SECTORS].sum())
        mean_by = mean_by.reset_index().melt(id_vars="year", var_name="sector",
                                             value_name="mean_quality")
        mean_by["measure"] = c
        rows.append(mean_by)
    D = pd.concat(rows, ignore_index=True)
    D.to_csv(out_dir / "quality_by_sector.csv", index=False)
    print("\n\nMean quality per patent by sector. This is the diagnostic that decides")
    print("everything: if RI-Services patents are not worse, weighting cannot help.\n")
    for c in have:
        p = D[D.measure == c].pivot(index="year", columns="sector",
                                    values="mean_quality")[SECTORS]
        p["RI_over_Goods"] = p.RI_Services / p.Goods
        yrs = [y for y in (1980, 1990, 2000, 2005, 2010, 2015) if y in p.index]
        print(f"{c}")
        print(p.loc[yrs].round(3).to_string())
        print()

    # ------------------------------------------------------------- weighted flows
    frames = []
    g = M.groupby("year")[SECTORS].sum().reset_index().melt(
        id_vars="year", var_name="sector", value_name="patents")
    g["measure"] = "unweighted"
    frames.append(g)
    for c in have:
        d = M.dropna(subset=[c]).copy()
        for lab, x in (("plus1", d[c] + 1.0), ("raw", d[c].astype(float))):
            mu = x.groupby(d["year"]).transform("mean")
            wt = np.where(mu > 0, x / mu, 0.0)
            e = d[SECTORS].mul(wt, axis=0)
            e["year"] = d["year"].to_numpy()
            gg = (e.groupby("year")[SECTORS].sum().reset_index()
                    .melt(id_vars="year", var_name="sector", value_name="patents"))
            gg["measure"] = f"{c}_{lab}"
            frames.append(gg)
    out = pd.concat(frames, ignore_index=True)[["measure", "year", "sector", "patents"]]
    out.to_csv(out_dir / "quality_sector_flows.csv", index=False)
    print(f"Wrote {out_dir / 'quality_sector_flows.csv'}, "
          f"{out_dir / 'quality_by_sector.csv'}, {out_dir / 'quality_coverage.csv'}")

    print("\nRI-Services share of the flow, by measure, every fifth year")
    piv = {}
    for m in out.measure.unique():
        p = out[out.measure == m].pivot(index="year", columns="sector", values="patents")
        piv[m] = p.RI_Services / p.sum(axis=1)
    P = pd.DataFrame(piv)
    print(P.loc[P.index % 5 == 0].round(4).to_string())


if __name__ == "__main__":
    main()
