"""
==========================================================================================
 An Anatomy of the Baumol Cost Disease -- Elbadri & Vizcaino (2026)
 ------------------------------------------------------------------------------------------
 File:    robustness/r3_triadic_flows.py
 Purpose: Patent flows by sector restricted to triadic families (EPO, JPO and USPTO), input to r6
 ------------------------------------------------------------------------------------------
 Inputs:  input/OECD_TPF_Feb2025/202501_TPF_USPTO.txt   triadic membership by US grant number
          input/OECD_TPF_Feb2025/202501_TPF_Core.txt    family dates, coverage check only
          input/PatentsView/g_patent.tsv.zip, g_ipc_at_issue.tsv.zip
          input/PatentsView/g_inventor_disambiguated.tsv.zip   optional, US-inventor bases
          input/concordances/ipc4_to_isic_rev31_4.txt, isic31_4digit_to_sector.csv
 Outputs: output/robustness/tpf_sector_flows.csv   basis (all, triadic, us_any, triadic_us), year,
                                                   sector, patents
 ------------------------------------------------------------------------------------------
 Notes:   Not in the paper. Asks whether the RI-Services expansion survives the triadic
          filter, which removes about ninety percent of grants. Recent cohorts are
          truncated; coverage by year is printed. Mass is checked.
==========================================================================================
"""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))   # code/, for paths.py
import paths as PATHS  # noqa: E402

import argparse
import sys
import traceback
from pathlib import Path

import numpy as np
import pandas as pd

CHUNK, BLOCK = 2_000_000, 4_000_000
YEAR_MIN = 1976
SECTORS = ["Goods", "RI_Services", "NRI_Services"]
RI_DIVISIONS = {"22", "64", "65", "66", "67", "72", "73", "74", "92"}


_ap = argparse.ArgumentParser(add_help=True)
_ap.add_argument("--out", help="output folder, defaults to output/robustness")
ARGS, _ = _ap.parse_known_args()


# ----------------------------------------------------------------- locating things

def locate(name):
    """Input folders, as set in code/paths.py."""
    return {"PatentsView": PATHS.PATENTSVIEW, "Concordances": PATHS.CONCORDANCES,
            "tpf": PATHS.OECD_TPF, "quality": PATHS.OECD_QUALITY, "kpss": PATHS.KPSS}[name]


def pick(folder, stem):
    """Find one data file by stem, whatever extension it carries."""
    for ext in (".tsv.zip", ".tsv", ".txt", ".txt.zip"):
        p = folder / (stem + ext)
        if p.exists():
            return p
    hits = sorted(folder.glob(f"*{stem}*"))
    hits = [h for h in hits if h.is_file()]
    return hits[0] if hits else None


def read_header(p, sep=None):
    if sep is None:
        sep = "|" if "".join(p.suffixes).startswith(".txt") else "\t"
    return list(pd.read_csv(p, sep=sep, nrows=0, dtype=str).columns)


def resolve(cols, wanted, label):
    lk = {str(c).strip().strip('"').lower(): c for c in cols}
    for w in ([wanted] if isinstance(wanted, str) else wanted):
        if w.lower() in lk:
            return lk[w.lower()]
    sys.exit(f"No column for {label}. Tried {wanted}. File has: {sorted(cols)}")


def to_pid(s):
    """Patent identifiers as a nullable integer. Non-utility grants carry letter prefixes."""
    return pd.to_numeric(s, errors="coerce").astype("Int64")


# ----------------------------------------------------------------- concordance

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
    cw["cls"] = (cw[resolve(cw.columns, ["ipc4", "ipc"], "IPC4 code")]
                 .str.strip().str.strip('"').str.upper())
    cw["isic"] = (cw[resolve(cw.columns, ["isic_rev31_4", "isic"], "ISIC code")]
                  .str.strip().str.strip('"').str.zfill(4))
    cw["w"] = pd.to_numeric(
        cw[resolve(cw.columns, ["probability_weight", "weight"], "probability weight")],
        errors="coerce").fillna(0.0)

    cw["sector"] = cw["isic"].map(to_sec)
    gap = cw["sector"].isna()
    if gap.any():
        share = cw.loc[gap, "w"].sum() / cw["w"].sum()
        cw.loc[gap, "sector"] = cw.loc[gap, "isic"].map(fallback_sector)
        print(f"  {cw.loc[gap, 'isic'].nunique()} ISIC codes absent from the sector map "
              f"({share:.5f} of weight), assigned by the division rule")
    cw = cw.dropna(subset=["sector"])

    S = cw.groupby(["cls", "sector"])["w"].sum().unstack(fill_value=0.0)
    for s in SECTORS:
        if s not in S.columns:
            S[s] = 0.0
    S = S[SECTORS]
    print(f"  {len(S):,} IPC subclasses mapped")
    return S.div(S.sum(axis=1), axis=0)      # row stochastic, so mass is preserved


# ----------------------------------------------------------------- triadic membership

def triadic_ids(tpf_dir):
    f = pick(tpf_dir, "202501_TPF_USPTO")
    if f is None:
        sys.exit(f"Not found: 202501_TPF_USPTO.txt in {tpf_dir}")
    print(f"\nStreaming {f.name} ...")
    head = read_header(f, sep="|")
    c_nbr = resolve(head, ["USPTO_Nbr", "uspto_nbr"], "US grant number in TPF_USPTO")
    print(f"  columns: {head}")

    parts, n = [], 0
    for ch in pd.read_csv(f, sep="|", usecols=[c_nbr], dtype=str,
                          chunksize=CHUNK, on_bad_lines="skip"):
        n += len(ch)
        parts.append(ch[c_nbr].str.strip().str.strip('"').str.upper().str.removeprefix("US"))
        print(f"  {n:,} rows read")
    raw = pd.concat(parts, ignore_index=True).dropna()
    del parts

    num = raw.str.fullmatch(r"\d+")
    print(f"  {len(raw):,} US members of triadic families, "
          f"{num.mean():.4f} are plain numeric grant numbers")
    if (~num).any():
        print(f"    non-numeric examples: {raw[~num].drop_duplicates().head(6).tolist()}")
    ids = pd.to_numeric(raw[num], errors="coerce").dropna().astype("int64").unique()
    print(f"  {len(ids):,} distinct triadic USPTO grant numbers, "
          f"range {ids.min():,} to {ids.max():,}")
    return set(ids.tolist())


def coverage_by_year(tpf_dir):
    f = pick(tpf_dir, "202501_TPF_Core")
    if f is None:
        print("\nTPF_Core not found, skipping the coverage check")
        return
    head = read_header(f, sep="|")
    c_g = resolve(head, ["USPTO_grant", "uspto_grant"], "USPTO grant date in TPF_Core")
    core = pd.read_csv(f, sep="|", usecols=[c_g], dtype=str, on_bad_lines="skip")
    yr = pd.to_numeric(core[c_g].str[:4], errors="coerce").dropna().astype(int)
    t = yr[(yr >= 1976) & (yr <= 2025)].value_counts().sort_index()
    print("\nTriadic families by USPTO grant year, every second year")
    print(t[t.index % 2 == 0].to_string())
    print("\n  The fall-off at the end is truncation, since a family must clear three offices")
    print("  before it appears. It sets the usable window, not a decline in triadic patenting.")


# ----------------------------------------------------------------- patent side

def grant_years(path):
    print(f"\nReading {path.name} ...")
    h = read_header(path, sep="\t")
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
    print(f"  {len(s):,} utility patents, grant years {s.min()} to {s.max()}")
    return s


def sector_vectors(path, S):
    """One sector probability vector per patent, from its IPC subclasses."""
    print(f"\nStreaming {path.name} ...")
    h = read_header(path, sep="\t")
    want = {"patent_id": ["patent_id"], "section": ["section"],
            "ipc_class": ["ipc_class", "class"], "subclass": ["subclass"]}
    cols = {k: resolve(h, v, k) for k, v in want.items()}
    print(f"  using {', '.join(cols.values())}")

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

    # Deduplicate across the whole file, not within chunks, so a patent whose rows straddle
    # a chunk boundary is not counted twice.
    pairs = pd.concat(parts, ignore_index=True).drop_duplicates()
    del parts
    print(f"  concordance covers {cov / max(n, 1):.4f} of rows, "
          f"{len(pairs):,} covered pairs, {pairs.patent_id.nunique():,} patents")

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

    total, n_pat = out[SECTORS].to_numpy().sum(), len(out)
    print(f"  {n_pat:,} patents retained, total weight {total:,.1f}")
    if abs(total - n_pat) > max(1.0, 1e-9 * n_pat):
        sys.exit(f"Mass check FAILED: {total:,.1f} against {n_pat:,} patents.")
    print("  mass check passed")
    return out


def us_weights(pv):
    """One if the patent has at least one US-resident inventor, else zero."""
    f = pick(pv, "g_inventor_disambiguated")
    if f is None:
        print("\ng_inventor_disambiguated not found, skipping the US bases")
        return None
    h = read_header(f, sep="\t")
    print(f"\nStreaming {f.name} ...\n  columns: {h}")
    c_id = resolve(h, ["patent_id"], "patent id")
    cc = next((c for c in h if str(c).lower() in
               ("disambig_country", "country", "inventor_country", "country_code")), None)
    lc = next((c for c in h if str(c).lower() in
               ("location_id", "disambig_location_id")), None)
    f_loc = None
    if cc:
        use = [c_id, cc]
        print(f"  using the country column {cc}")
    elif lc:
        f_loc = pick(pv, "g_location_disambiguated")
        if f_loc is None:
            print("  no country column and no g_location_disambiguated in the PatentsView "
                  "folder, skipping the US bases")
            return None
        use = [c_id, lc]
        print(f"  joining {f_loc.name} on {lc}")
    else:
        print("  no country or location column, skipping the US bases")
        return None

    parts, n = [], 0
    for ch in pd.read_csv(f, sep="\t", usecols=use, dtype=str, chunksize=CHUNK,
                          on_bad_lines="skip"):
        n += len(ch)
        ch = ch.rename(columns={c_id: "patent_id"})
        ch["patent_id"] = to_pid(ch["patent_id"])
        parts.append(ch.dropna(subset=["patent_id"]))
        print(f"  {n:,} rows read")
    inv = pd.concat(parts, ignore_index=True)
    del parts
    inv["patent_id"] = inv["patent_id"].astype("int64")

    if cc:
        inv["is_us"] = inv[cc].astype(str).str.strip().str.upper().eq("US")
    else:
        loc = pd.read_csv(f_loc, sep="\t", dtype=str, on_bad_lines="skip")
        k = resolve(loc.columns, ["location_id", "id"], "location id")
        v = resolve(loc.columns, ["disambig_country", "country", "country_code"],
                    "country in g_location_disambiguated")
        mp = loc.set_index(k)[v].astype(str).str.strip().str.upper()
        inv["is_us"] = inv[lc].map(mp).eq("US")

    g = inv.groupby("patent_id")["is_us"].max().astype(float)
    print(f"  {len(g):,} patents with inventor records, "
          f"{g.mean():.4f} with at least one US inventor")
    return g


# ----------------------------------------------------------------- main

def main():
    pv = locate("PatentsView")
    conc = locate("Concordances")
    tpf = locate("tpf")
    out_dir = Path(ARGS.out) if ARGS.out else PATHS.OUT_ROBUSTNESS
    out_dir.mkdir(exist_ok=True, parents=True)
    print(f"PatentsView:  {pv}\nConcordances: {conc}\nTPF:          {tpf}\n"
          f"Output:       {out_dir}\n")

    f_pat, f_ipc = pick(pv, "g_patent"), pick(pv, "g_ipc_at_issue")
    for f, nm in ((f_pat, "g_patent"), (f_ipc, "g_ipc_at_issue")):
        if f is None:
            sys.exit(f"Not found in {pv}: {nm}")

    print("Building the IPC to sector table ...")
    S = ipc_sector_table(conc)

    tri = triadic_ids(tpf)
    coverage_by_year(tpf)

    year_of = grant_years(f_pat)
    V = sector_vectors(f_ipc, S)

    V["year"] = V["patent_id"].map(year_of)
    V = V.dropna(subset=["year"])
    V["year"] = V["year"].astype(int)
    V = V[V["year"] >= YEAR_MIN].reset_index(drop=True)
    V["triadic"] = V["patent_id"].isin(tri).astype(float)
    print(f"\n{len(V):,} patents dated and in sample, {V.year.min()} to {V.year.max()}")
    print(f"  {V.triadic.mean():.4f} belong to a triadic family")

    bases = {"all": np.ones(len(V)), "triadic": V["triadic"].to_numpy()}
    try:
        W = us_weights(pv)
        if W is not None:
            us = V["patent_id"].map(W).fillna(0.0).to_numpy()
            bases["us_any"] = us
            bases["triadic_us"] = us * V["triadic"].to_numpy()
    except Exception:
        print("\nUS bases failed, continuing without them:\n" + traceback.format_exc())

    frames = []
    for name, w in bases.items():
        d = V[SECTORS].mul(w, axis=0)
        d["year"] = V["year"].to_numpy()
        g = (d.groupby("year")[SECTORS].sum().reset_index()
               .melt(id_vars="year", var_name="sector", value_name="patents"))
        g["basis"] = name
        frames.append(g)
    out = pd.concat(frames, ignore_index=True)[["basis", "year", "sector", "patents"]]
    out.to_csv(out_dir / "tpf_sector_flows.csv", index=False)
    print(f"\nWrote {out_dir / 'tpf_sector_flows.csv'}")

    tot = out.groupby(["basis", "year"])["patents"].sum().unstack("basis")
    print("\ntotal patents per year by basis, every fifth year")
    print(tot.loc[tot.index % 5 == 0].round(0).to_string())

    for b in ("all", "triadic"):
        sh = out[out.basis == b].pivot(index="year", columns="sector", values="patents")
        sh = sh.div(sh.sum(axis=1), axis=0)
        print(f"\n{b} sector shares, every fifth year")
        print(sh.loc[sh.index % 5 == 0].round(4).to_string())


if __name__ == "__main__":
    main()
