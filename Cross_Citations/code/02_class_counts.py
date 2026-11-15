"""
==========================================================================================
 An Anatomy of the Baumol Cost Disease -- Elbadri & Vizcaino (2026)
 ------------------------------------------------------------------------------------------
 File:    02_class_counts.py
 Purpose: Utility patents by technology class and grant year, 1976 onward. Everything
          downstream (sector aggregation, exclusions, partition check) runs on these small
          tables in step 05, without the bulk PatentsView files
 ------------------------------------------------------------------------------------------
 Inputs:  input/PatentsView/g_patent.tsv.zip          grant year, utility filter
          input/PatentsView/g_ipc_at_issue.tsv.zip    IPC subclasses at issue
          input/PatentsView/g_uspc_at_issue.tsv.zip   USPC main class, optional, 1976-2015
 Outputs: output/patents/ipc4_year_counts.csv         ipc4, year, patents
          output/patents/ipc4_year_counts_exbm.csv    the same, excluding patents with an
                                                      at-issue G06Q subclass or USPC class 705
          output/patents/uspc_year_counts.csv         uspc_mainclass, year, patents (if the
                                                      USPC file is present)
 ------------------------------------------------------------------------------------------
 Notes:   A patent contributes total weight one, split equally across its distinct
          subclasses, so summing over classes within a year reproduces the annual grant
          count. Deduplication is across the whole file, not within chunks.
==========================================================================================
"""

import sys

import pandas as pd

import paths as P

SAVE_TABLES = True

CHUNK = 2_000_000
YEAR_MIN = 1976


# ----------------------------------------------------------------- small helpers

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


# ----------------------------------------------------------------- grant years

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


# ----------------------------------------------------------------- class counts

def class_pairs(path, build_code, label):
    """Distinct (patent_id, class) pairs over the whole file."""
    print(f"\nStreaming {path.name} ...")
    parts, n_rows = [], 0
    for i, ch in enumerate(pd.read_csv(path, sep="\t", dtype=str, chunksize=CHUNK,
                                       on_bad_lines="skip"), 1):
        n_rows += len(ch)
        pid, cls = build_code(ch)
        keep = pid.notna().to_numpy() & cls.notna().to_numpy()
        parts.append(pd.DataFrame({"patent_id": pid.to_numpy()[keep].astype("int64"),
                                   label: cls.to_numpy()[keep]}))
        print(f"  chunk {i}: {n_rows:,} rows read")
    pairs = pd.concat(parts, ignore_index=True).drop_duplicates()
    del parts
    print(f"  {n_rows:,} rows, {len(pairs):,} distinct pairs, "
          f"{pairs.patent_id.nunique():,} patents, {pairs[label].nunique():,} classes")
    return pairs


def counts_by_year(pairs, label, year_of):
    """Each patent contributes total weight one, split equally across its distinct classes."""
    pairs = pairs.copy()
    pairs["w"] = 1.0 / pairs.groupby("patent_id")[label].transform("size")
    pairs["year"] = pairs["patent_id"].map(year_of)
    pairs = pairs.dropna(subset=["year"])
    pairs["year"] = pairs["year"].astype(int)
    pairs = pairs[pairs["year"] >= YEAR_MIN]

    out = (pairs.groupby([label, "year"])["w"].sum()
                .reset_index().rename(columns={"w": "patents"}))
    tot = out.groupby("year")["patents"].sum()
    n = pairs.groupby("year")["patent_id"].nunique()
    gap = (tot - n).abs().max()
    print(f"  annual totals reproduce the patent count to {gap:.6f}")
    if gap > 1e-6:
        sys.exit("Mass check FAILED.")
    print(f"  {out.year.min()} to {out.year.max()}, "
          f"{int(round(tot.loc[out.year.max()])):,} patents in the last year")
    return out


def ipc_code(ch):
    head = list(ch.columns)
    c_id = resolve(head, ["patent_id"], "patent id")
    c_se = resolve(head, ["section"], "IPC section")
    c_cl = resolve(head, ["ipc_class", "class"], "IPC class")
    c_sb = resolve(head, ["subclass"], "IPC subclass")
    ch = ch.dropna(subset=[c_id, c_se, c_cl, c_sb])
    cls = (ch[c_se].str.strip().str.upper()
           + ch[c_cl].str.strip().str.zfill(2)
           + ch[c_sb].str.strip().str.upper())
    return to_pid(ch[c_id]), cls


def uspc_code(ch):
    head = list(ch.columns)
    c_id = resolve(head, ["patent_id"], "patent id")
    c_mc = resolve(head, ["uspc_mainclass_id", "mainclass_id", "mainclass"],
                   "USPC main class")
    ch = ch.dropna(subset=[c_id, c_mc])
    return to_pid(ch[c_id]), ch[c_mc].str.strip().str.upper()


# ----------------------------------------------------------------- main

def main():
    out_dir = P.OUT_PATENTS
    out_dir.mkdir(parents=True, exist_ok=True)
    f_pat = P.pv_file("g_patent")
    f_ipc = P.pv_file("g_ipc_at_issue")

    year_of = grant_years(f_pat)

    print("\n=== IPC subclasses ===")
    ipc_pairs = class_pairs(f_ipc, ipc_code, "ipc4")
    ipc = counts_by_year(ipc_pairs, "ipc4", year_of)
    if SAVE_TABLES:
        ipc.to_csv(out_dir / "ipc4_year_counts.csv", index=False)
        print(f"  wrote {out_dir / 'ipc4_year_counts.csv'}  ({len(ipc):,} rows)")

    # ------------------------------------------------------------- USPC main classes
    uspc_pairs = None
    f_uspc = next((p for p in (P.PATENTSVIEW / "g_uspc_at_issue.tsv.zip",
                               P.PATENTSVIEW / "g_uspc_at_issue.tsv") if p.exists()), None)
    if f_uspc is None:
        print("\ng_uspc_at_issue not found. The G06Q exclusion below runs on the IPC route "
              "alone; USPC 705 is a cross check only.")
    else:
        print("\n=== USPC main classes ===")
        try:
            uspc_pairs = class_pairs(f_uspc, uspc_code, "uspc_mainclass")
            uspc = counts_by_year(uspc_pairs, "uspc_mainclass", year_of)
            if SAVE_TABLES:
                uspc.to_csv(out_dir / "uspc_year_counts.csv", index=False)
                print(f"  wrote {out_dir / 'uspc_year_counts.csv'}  ({len(uspc):,} rows)")
        except Exception as e:  # the USPC file is superseded, so a failure here is not fatal
            print(f"  USPC pass failed: {e}\n  The IPC output above is unaffected.")
            uspc_pairs = None

    # ------------------------------------------------------- business method exclusion
    # A patent is dropped whole if it carries any G06Q subclass or sits in USPC main class
    # 705. Dropping the class rather than the patent would leave a business method patent
    # in the count through its secondary classes, which is not the counterfactual wanted.
    print("\n=== Business method exclusion ===")
    drop = set(ipc_pairs.loc[ipc_pairs["ipc4"] == "G06Q", "patent_id"].unique())
    print(f"  {len(drop):,} patents carry a G06Q subclass")
    if uspc_pairs is not None:
        u705 = set(uspc_pairs.loc[uspc_pairs["uspc_mainclass"].astype(str).str.lstrip("0")
                                  == "705", "patent_id"].unique())
        print(f"  {len(u705):,} patents sit in USPC main class 705, "
              f"{len(drop & u705):,} of them also in G06Q")
        drop |= u705
    print(f"  {len(drop):,} patents dropped in total")

    kept = ipc_pairs[~ipc_pairs["patent_id"].isin(drop)]
    exbm = counts_by_year(kept, "ipc4", year_of)
    if SAVE_TABLES:
        exbm.to_csv(out_dir / "ipc4_year_counts_exbm.csv", index=False)
        print(f"  wrote {out_dir / 'ipc4_year_counts_exbm.csv'}  ({len(exbm):,} rows)")

    a = ipc.groupby("year")["patents"].sum()
    b = exbm.groupby("year")["patents"].sum()
    sh = (1 - b / a).dropna()
    print("\n  share of grants excluded, selected years:")
    for y in [1976, 1990, 1998, 2005, 2010, 2014, 2019, 2023]:
        if y in sh.index:
            print(f"    {y}  {sh.loc[y]:.4f}")


if __name__ == "__main__":
    main()
