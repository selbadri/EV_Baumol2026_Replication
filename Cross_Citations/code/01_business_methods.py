"""
==========================================================================================
 An Anatomy of the Baumol Cost Disease -- Elbadri & Vizcaino (2026)
 ------------------------------------------------------------------------------------------
 File:    01_business_methods.py
 Purpose: Flag business-method patents (current CPC subclass G06Q) and count utility
          patents by IPC subclass and grant year without them (Appendix C)
 ------------------------------------------------------------------------------------------
 Inputs:  input/PatentsView/g_patent.tsv.zip         grant year, utility filter
          input/PatentsView/g_ipc_at_issue.tsv.zip   IPC subclasses at issue
          input/PatentsView/g_cpc_current.tsv.zip    current CPC classification
 Outputs: output/patents/bm_ids_cpc.csv              flagged patent ids (also used by step 06)
          output/patents/ipc4_year_counts_excpc.csv  ipc4, year, patents, excluding G06Q
          output/patents/bm_cpc_coverage.csv         flagged patents and share of grants, by year
 ------------------------------------------------------------------------------------------
 Notes:   The current CPC applies one definition to every grant year, including patents
          reclassified after issue. The at-issue IPC flag does not: G06Q was created in the
          2006 IPC reform, and before then business methods sat in G06F 17/60. A patent
          contributes total weight one, split equally across its distinct IPC subclasses,
          deduplicated across the whole file; the script exits if the mass check fails.
==========================================================================================
"""

import re
import sys

import pandas as pd

import paths as P

SAVE_TABLES = True

CHUNK = 2_000_000
YEAR_MIN = 1976
TARGET = "G06Q"

FULL_CODE = re.compile(r"^[A-HY]\d{2}[A-Z]$")


# ----------------------------------------------------------------- small helpers

def read_header(path):
    return list(pd.read_csv(path, sep="\t", nrows=0, dtype=str).columns)


def resolve(available, wanted, label, required=True):
    lookup = {str(c).lower(): c for c in available}
    for name in ([wanted] if isinstance(wanted, str) else wanted):
        if name.lower() in lookup:
            return lookup[name.lower()]
    if not required:
        return None
    sys.exit(f"No column for {label}. Tried {wanted}. File contains: {sorted(available)}")


def to_pid(series):
    """Patent identifiers as a nullable integer. Non-utility grants carry letter prefixes."""
    return pd.to_numeric(series, errors="coerce").astype("Int64")


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
    return s


# ----------------------------------------------------------------- the CPC flag

def cpc_subclass_builder(head, sample):
    """CPC files store the subclass either as one code, G06Q, or as three columns."""
    c_id = resolve(head, ["patent_id"], "patent id")
    c_sb = resolve(head, ["cpc_subclass", "subclass", "cpc_subclass_id"], "CPC subclass")
    looks_full = sample[c_sb].dropna().astype(str).str.strip().str.upper() \
                            .str.match(FULL_CODE).mean() > 0.9
    if looks_full:
        print(f"  {c_sb} already holds the full subclass code")
        return c_id, (lambda ch: ch[c_sb].astype(str).str.strip().str.upper())
    c_se = resolve(head, ["cpc_section", "section", "cpc_section_id"], "CPC section")
    c_cl = resolve(head, ["cpc_class", "class", "cpc_class_id"], "CPC class")
    print(f"  building the subclass from {c_se}, {c_cl}, {c_sb}")
    return c_id, (lambda ch: (ch[c_se].astype(str).str.strip().str.upper()
                              + ch[c_cl].astype(str).str.strip().str.zfill(2)
                              + ch[c_sb].astype(str).str.strip().str.upper()))


def flagged_ids(path):
    """Patent identifiers whose current CPC classification contains the target subclass."""
    print(f"\nStreaming {path.name} for current CPC {TARGET} ...")
    head = read_header(path)
    sample = pd.read_csv(path, sep="\t", nrows=2000, dtype=str, on_bad_lines="skip")
    c_id, build = cpc_subclass_builder(head, sample)

    ids, n_rows = set(), 0
    for i, ch in enumerate(pd.read_csv(path, sep="\t", dtype=str, chunksize=CHUNK,
                                       on_bad_lines="skip"), 1):
        n_rows += len(ch)
        sub = build(ch)
        hit = ch.loc[sub == TARGET, c_id]
        ids.update(to_pid(hit).dropna().astype("int64").tolist())
        print(f"  chunk {i}: {n_rows:,} rows read, {len(ids):,} flagged so far")
    print(f"  {n_rows:,} rows, {len(ids):,} patents carry current CPC {TARGET}")
    return ids


# ----------------------------------------------------------------- IPC counts

def ipc_pairs(path):
    print(f"\nStreaming {path.name} ...")
    head = read_header(path)
    c_id = resolve(head, ["patent_id"], "patent id")
    c_se = resolve(head, ["section"], "IPC section")
    c_cl = resolve(head, ["ipc_class", "class"], "IPC class")
    c_sb = resolve(head, ["subclass"], "IPC subclass")
    parts, n_rows = [], 0
    for i, ch in enumerate(pd.read_csv(path, sep="\t", usecols=[c_id, c_se, c_cl, c_sb],
                                       dtype=str, chunksize=CHUNK,
                                       on_bad_lines="skip"), 1):
        n_rows += len(ch)
        ch = ch.dropna(subset=[c_id, c_se, c_cl, c_sb])
        cls = (ch[c_se].str.strip().str.upper()
               + ch[c_cl].str.strip().str.zfill(2)
               + ch[c_sb].str.strip().str.upper())
        pid = to_pid(ch[c_id])
        keep = pid.notna().to_numpy()
        parts.append(pd.DataFrame({"patent_id": pid.to_numpy()[keep].astype("int64"),
                                   "ipc4": cls.to_numpy()[keep]}))
        print(f"  chunk {i}: {n_rows:,} rows read")
    pairs = pd.concat(parts, ignore_index=True).drop_duplicates()
    print(f"  {n_rows:,} rows, {len(pairs):,} distinct pairs, "
          f"{pairs.patent_id.nunique():,} patents")
    return pairs


def counts_by_year(pairs, year_of):
    pairs = pairs.copy()
    pairs["w"] = 1.0 / pairs.groupby("patent_id")["ipc4"].transform("size")
    pairs["year"] = pairs["patent_id"].map(year_of)
    pairs = pairs.dropna(subset=["year"])
    pairs["year"] = pairs["year"].astype(int)
    pairs = pairs[pairs["year"] >= YEAR_MIN]
    out = (pairs.groupby(["ipc4", "year"])["w"].sum()
                .reset_index().rename(columns={"w": "patents"}))
    gap = (out.groupby("year")["patents"].sum()
           - pairs.groupby("year")["patent_id"].nunique()).abs().max()
    print(f"  annual totals reproduce the patent count to {gap:.6f}")
    if gap > 1e-6:
        sys.exit("Mass check FAILED.")
    return out


# ----------------------------------------------------------------- main

def main():
    out_dir = P.OUT_PATENTS
    out_dir.mkdir(parents=True, exist_ok=True)
    f_pat = P.pv_file("g_patent")
    f_ipc = P.pv_file("g_ipc_at_issue")
    f_cpc = P.pv_file("g_cpc_current")

    year_of = grant_years(f_pat)
    bm = flagged_ids(f_cpc)
    pairs = ipc_pairs(f_ipc)

    kept = pairs[~pairs["patent_id"].isin(bm)]
    print("\n=== Counts excluding current CPC G06Q ===")
    out = counts_by_year(kept, year_of)
    if SAVE_TABLES:
        out.to_csv(out_dir / "ipc4_year_counts_excpc.csv", index=False)
        print(f"  wrote {out_dir / 'ipc4_year_counts_excpc.csv'}  ({len(out):,} rows)")
        pd.DataFrame({"patent_id": sorted(bm)}).to_csv(out_dir / "bm_ids_cpc.csv",
                                                       index=False)

    # ------------------------------------------------------------- coverage report
    yr = pd.Series(sorted(bm), name="patent_id").map(year_of).dropna().astype(int)
    flagged = yr[yr >= YEAR_MIN].value_counts().sort_index()
    total = (pd.Series(year_of.to_numpy(), name="year")
             .value_counts().sort_index().loc[YEAR_MIN:])
    cov = pd.DataFrame({"flagged": flagged, "grants": total}).fillna(0)
    cov["share"] = cov["flagged"] / cov["grants"]
    if SAVE_TABLES:
        cov.to_csv(out_dir / "bm_cpc_coverage.csv")
    print("\nPatents flagged as business methods by current CPC, and their share of grants:\n")
    for y in [1976, 1980, 1985, 1990, 1995, 1998, 2000, 2005, 2010, 2014, 2019, 2023]:
        if y in cov.index:
            print(f"  {y}  {int(cov.loc[y, 'flagged']):>7,}   {cov.loc[y, 'share']:.4f}")
    pre = cov.loc[:2005, "flagged"].sum()
    print(f"\n  {int(pre):,} flagged patents were granted before 2006. These are the ones the "
          f"at-issue\n  IPC route could not see, and they are why that route was not "
          f"comparable over time.")


if __name__ == "__main__":
    main()
