"""
==========================================================================================
 An Anatomy of the Baumol Cost Disease -- Elbadri & Vizcaino (2026)
 ------------------------------------------------------------------------------------------
 File:    robustness/r1_concordance_overlap.py
 Purpose: IPC route against USPC route: how much of the sector composition comes from the
          choice of source classification, on patents carrying both, 1976-2014
 ------------------------------------------------------------------------------------------
 Inputs:  input/PatentsView/g_patent.tsv.zip, g_ipc_at_issue.tsv.zip, g_uspc_at_issue.tsv.zip
          input/concordances/ipc4_to_isic_rev31_4.txt
          input/concordances/uspc_Class_to_isic_rev31_4.txt   ALP USPC -> ISIC Rev 3.1
          input/concordances/isic31_4digit_to_sector.csv
 Outputs: output/robustness/overlap_annual_shares.csv      sector shares by year and route
          output/robustness/overlap_agreement_matrix.csv   pooled 3x3 agreement matrix
          output/robustness/overlap_annual_agreement.csv   agreement rate by year
 ------------------------------------------------------------------------------------------
 Notes:   Not in the paper. The two routes agree on sector margins to within about 0.01,
          which supports using the IPC route after the USPC was retired in 2015. Each
          patent has total weight one on each route; mass is checked at every step.
==========================================================================================
"""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))   # code/, for paths.py
import paths as PATHS  # noqa: E402

import sys
from pathlib import Path

import numpy as np
import pandas as pd

CHUNK = 2_000_000
BLOCK = 4_000_000
YEAR_MIN, YEAR_MAX = 1976, 2014
SECTORS = ["Goods", "RI_Services", "NRI_Services"]


# ----------------------------------------------------------------- small helpers



def locate(name):
    """Input folders, as set in code/paths.py."""
    return {"PatentsView": PATHS.PATENTSVIEW, "Concordances": PATHS.CONCORDANCES,
            "tpf": PATHS.OECD_TPF, "quality": PATHS.OECD_QUALITY, "kpss": PATHS.KPSS}[name]


def read_header(path):
    return list(pd.read_csv(path, sep="\t", nrows=0, dtype=str).columns)


def resolve(available, wanted, label):
    lookup = {str(c).lower(): c for c in available}
    for name in ([wanted] if isinstance(wanted, str) else wanted):
        if name.lower() in lookup:
            return lookup[name.lower()]
    sys.exit(f"No column for {label}. Tried {wanted}. File contains: {sorted(available)}")


def to_pid(series):
    """Patent identifiers as int32.

    Utility grant numbers are numeric. Reissues, plant patents and design patents carry a
    letter prefix, are not in the sample, and become missing here so the pairs frames stay
    small enough to hold in memory.
    """
    return pd.to_numeric(series, errors="coerce").astype("Int64")


def normalise_uspc(series):
    """Put USPC class codes into one format: 2 and 02 both become 002."""
    s = series.astype(str).str.strip().str.upper().str.strip('"')
    return s.where(~s.str.fullmatch(r"\d+"), s.str.zfill(3))


# ----------------------------------------------------------------- concordances

# The research intensive service divisions, at the two digit level of ISIC Rev 3.1, that
# the sector map assigns to RI-Services. Used only as a fallback for codes the map file
# does not list, so that no citation mass is silently discarded.
RI_DIVISIONS = {"22", "64", "65", "66", "67", "72", "73", "74", "92"}


def fallback_sector(code):
    """Assign a four digit ISIC Rev 3.1 code the sector map does not list."""
    c = str(code).zfill(4)
    if not c.isdigit():
        return None
    if int(c) < 5000:
        return "Goods"
    return "RI_Services" if c[:2] in RI_DIVISIONS else "NRI_Services"


def sector_map(path):
    """Return a Series from four digit ISIC Rev 3.1 code to one of three sectors."""
    m = pd.read_csv(path, dtype=str)
    m.columns = [c.strip() for c in m.columns]
    m["isic_rev31_4"] = m["isic_rev31_4"].str.strip().str.zfill(4)
    s = m.set_index("isic_rev31_4")["Sector"].str.strip()
    # The paper works with three sectors. Manufacturing and Other_Goods are both Goods.
    return s.replace({"Manufacturing": "Goods", "Other_Goods": "Goods"})


def class_to_sector(path, class_names, isic_names, isic_to_sector, uspc=False):
    """Collapse a class to ISIC concordance into a row stochastic class by sector table."""
    cw = pd.read_csv(path, dtype=str)
    cw.columns = [c.strip() for c in cw.columns]
    c_cls = resolve(cw.columns, class_names, f"class code in {Path(path).name}")
    c_isic = resolve(cw.columns, isic_names, f"ISIC code in {Path(path).name}")
    c_w = resolve(cw.columns, ["probability_weight", "weight"], "probability weight")

    cw["cls"] = cw[c_cls].str.strip().str.strip('"').str.upper()
    if uspc:
        cw["cls"] = normalise_uspc(cw["cls"])
    cw["isic"] = cw[c_isic].str.strip().str.strip('"').str.zfill(4)
    cw["w"] = pd.to_numeric(cw[c_w], errors="coerce").fillna(0.0)

    cw["sector"] = cw["isic"].map(isic_to_sector)
    gap = cw["sector"].isna()
    if gap.any():
        codes = sorted(cw.loc[gap, "isic"].unique())
        share = cw.loc[gap, "w"].sum() / cw["w"].sum()
        cw.loc[gap, "sector"] = cw.loc[gap, "isic"].map(fallback_sector)
        print(f"    {len(codes)} ISIC codes are absent from the sector map "
              f"({share:.5f} of concordance weight), assigned by the division rule: "
              f"{', '.join(codes[:10])}")
    cw = cw.dropna(subset=["sector"])

    S = cw.groupby(["cls", "sector"])["w"].sum().unstack(fill_value=0.0)
    for s in SECTORS:
        if s not in S.columns:
            S[s] = 0.0
    S = S[SECTORS]
    S = S.div(S.sum(axis=1), axis=0)          # row stochastic, so mass is preserved
    print(f"    {len(S):,} classes mapped")
    return S


# ----------------------------------------------------------------- patent side

def ipc4_builder(chunk):
    """Assemble the four character IPC subclass from its three PatentsView columns."""
    return (chunk["section"].str.strip().str.upper()
            + chunk["ipc_class"].str.strip().str.zfill(2)
            + chunk["subclass"].str.strip().str.upper())


def read_pairs(path, id_names, col_spec, build, S, label, uspc=False):
    """Stream a PatentsView class table and return distinct patent by class pairs.

    Classes the concordance does not cover are dropped inside the loop, which keeps the
    held frame small. The share of pairs dropped is reported.
    """
    print(f"\nStreaming {Path(path).name} ...")
    head = read_header(path)
    c_id = resolve(head, id_names, f"patent id in {label}")
    cols = {name: resolve(head, cands, f"{name} in {label}")
            for name, cands in col_spec.items()}
    print(f"  using {c_id}, " + ", ".join(cols.values()))

    usecols = [c_id] + list(cols.values())
    parts, n_rows, n_kept = [], 0, 0
    for i, chunk in enumerate(pd.read_csv(path, sep="\t", usecols=usecols, dtype=str,
                                          chunksize=CHUNK, on_bad_lines="skip"), 1):
        chunk = chunk.rename(columns={v: k for k, v in cols.items()})
        chunk = chunk.rename(columns={c_id: "patent_id"})
        n_rows += len(chunk)
        chunk = chunk.dropna(subset=["patent_id"] + list(cols))
        cls = build(chunk)
        if uspc:
            cls = normalise_uspc(cls)
        pid = to_pid(chunk["patent_id"])
        code = pd.Categorical(cls, categories=S.index).codes
        keep = (code >= 0) & pid.notna().to_numpy()
        n_kept += int((code >= 0).sum())
        parts.append(pd.DataFrame({"patent_id": pid.to_numpy()[keep].astype("int32"),
                                   "code": code[keep].astype("int16")}))
        print(f"  chunk {i}: {n_rows:,} rows read")

    # Deduplicate across the whole file, not within chunks, so a patent whose rows
    # straddle a chunk boundary is not counted twice.
    pairs = pd.concat(parts, ignore_index=True).drop_duplicates()
    del parts
    print(f"  {n_rows:,} rows, concordance covers {n_kept / max(n_rows, 1):.4f} of them")
    print(f"  {len(pairs):,} distinct covered patent-class pairs, "
          f"{pairs.patent_id.nunique():,} patents")
    return pairs


def sector_vectors(pairs, S, label):
    """Turn patent by class pairs into one sector probability vector per patent."""
    pairs = pairs.copy()
    # Weight after restricting to covered classes, so each retained patent still carries
    # total weight one and dropped classes do not shrink a patent's contribution.
    pairs["w"] = 1.0 / pairs.groupby("patent_id")["code"].transform("size")

    codes = pairs["code"].to_numpy()
    Sm = S.to_numpy(dtype="float64")
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

    total, n_pat = out[SECTORS].to_numpy().sum(), len(out)
    print(f"  {n_pat:,} patents retained, total weight {total:,.1f}")
    if abs(total - n_pat) > max(1.0, 1e-9 * n_pat):
        sys.exit(f"Mass check FAILED for {label}: {total:,.1f} against {n_pat:,} patents.")
    print("  mass check passed")
    return out


def grant_years(path):
    print(f"\nReading {Path(path).name} ...")
    head = read_header(path)
    c_id = resolve(head, ["patent_id"], "patent id in g_patent")
    c_dt = resolve(head, ["patent_date", "date", "grant_date"], "grant date")
    c_ty = resolve(head, ["patent_type", "type"], "patent type")
    pat = pd.read_csv(path, sep="\t", usecols=[c_id, c_dt, c_ty], dtype=str,
                      on_bad_lines="skip")
    pat = pat.rename(columns={c_id: "patent_id", c_dt: "date", c_ty: "type"})
    pat = pat[pat["type"].astype(str).str.lower() == "utility"]
    pat["patent_id"] = to_pid(pat["patent_id"])
    pat["year"] = pd.to_numeric(pat["date"].str[:4], errors="coerce")
    pat = pat.dropna(subset=["year", "patent_id"]).drop_duplicates("patent_id")
    pat["patent_id"] = pat["patent_id"].astype("int32")
    s = pat.set_index("patent_id")["year"].astype("int16")
    print(f"  {len(s):,} utility patents, grant years {s.min()} to {s.max()}")
    return s


# ----------------------------------------------------------------- main

def main():
    pv = locate("PatentsView")
    conc = locate("Concordances")
    out_dir = PATHS.OUT_ROBUSTNESS
    out_dir.mkdir(exist_ok=True, parents=True)

    f_ipc = pv / "g_ipc_at_issue.tsv.zip"
    if not f_ipc.exists():
        f_ipc = pv / "g_ipc_at_issue.tsv"
    f_uspc = pv / "g_uspc_at_issue.tsv.zip"
    if not f_uspc.exists():
        f_uspc = pv / "g_uspc_at_issue.tsv"
    f_pat = pv / "g_patent.tsv.zip"
    if not f_pat.exists():
        f_pat = pv / "g_patent.tsv"
    f_c_ipc = PATHS.ALP_IPC
    f_c_uspc = PATHS.ALP_USPC_ISIC
    f_sec = PATHS.SECTOR_MAP
    for f in (f_ipc, f_uspc, f_pat, f_c_ipc, f_c_uspc, f_sec):
        if not f.exists():
            sys.exit(f"Not found: {f}")
    print(f"PatentsView: {pv}\nConcordances: {conc}\n")

    print("Building class to sector tables ...")
    isic_to_sector = sector_map(f_sec)
    print("  IPC route:")
    S_ipc = class_to_sector(f_c_ipc, ["ipc4", "ipc"], ["isic_rev31_4", "isic"],
                            isic_to_sector)
    print("  USPC route:")
    S_uspc = class_to_sector(f_c_uspc, ["uspc", "uspc_class", "uspc_Class"],
                             ["isic_rev31_4", "isic"], isic_to_sector, uspc=True)

    year_of = grant_years(f_pat)

    p_ipc = read_pairs(f_ipc, ["patent_id"],
                       {"section": ["section"], "ipc_class": ["ipc_class", "class"],
                        "subclass": ["subclass"]},
                       ipc4_builder, S_ipc, "g_ipc_at_issue")
    v_ipc = sector_vectors(p_ipc, S_ipc, "IPC route")
    del p_ipc

    p_uspc = read_pairs(f_uspc, ["patent_id"],
                        {"cls_raw": ["uspc_mainclass_id", "mainclass_id", "uspc_class",
                                     "class_id", "uspc_mainclass"]},
                        lambda c: c["cls_raw"], S_uspc, "g_uspc_at_issue", uspc=True)
    v_uspc = sector_vectors(p_uspc, S_uspc, "USPC route")
    del p_uspc

    # ------------------------------------------------------------- common sample
    m = v_ipc.merge(v_uspc, on="patent_id", suffixes=("_ipc", "_uspc"))
    m["year"] = m["patent_id"].map(year_of)
    m = m.dropna(subset=["year"])
    m["year"] = m["year"].astype(int)
    m = m[(m["year"] >= YEAR_MIN) & (m["year"] <= YEAR_MAX)]
    print(f"\nCommon sample: {len(m):,} patents, {YEAR_MIN} to {YEAR_MAX}, "
          f"from {len(v_ipc):,} on the IPC route and {len(v_uspc):,} on the USPC route")

    A = m[[f"{s}_ipc" for s in SECTORS]].to_numpy()
    B = m[[f"{s}_uspc" for s in SECTORS]].to_numpy()

    # ------------------------------------------------------------- annual shares
    rows = []
    for route, X in (("ipc", A), ("uspc", B)):
        d = pd.DataFrame(X, columns=SECTORS)
        d["year"] = m["year"].to_numpy()
        g = d.groupby("year")[SECTORS].sum()
        g = g.div(g.sum(axis=1), axis=0)
        g["route"] = route
        rows.append(g.reset_index())
    shares = pd.concat(rows, ignore_index=True)[["route", "year"] + SECTORS]
    shares.to_csv(out_dir / "overlap_annual_shares.csv", index=False)

    for s in SECTORS:
        piv = shares.pivot(index="year", columns="route", values=s)
        piv["gap"] = piv["uspc"] - piv["ipc"]
        print(f"\n{s} share of patenting by route")
        print(piv.loc[piv.index % 5 == 0].round(4).to_string())

    # ------------------------------------------------------------- agreement
    # Row is the USPC assignment, column is the IPC assignment. Each patent contributes
    # the outer product of its two probability vectors, so total mass is the patent count.
    M = pd.DataFrame(B.T @ A, index=SECTORS, columns=SECTORS)
    M.index.name, M.columns.name = "uspc_route", "ipc_route"
    Mn = M / M.to_numpy().sum()
    Mn.to_csv(out_dir / "overlap_agreement_matrix.csv")
    print("\nPooled agreement matrix, shares of total patent mass")
    print(Mn.round(4).to_string())
    print(f"\n  probabilistic agreement rate (trace): {np.trace(Mn.to_numpy()):.4f}")

    # Chance agreement under the independence benchmark, so the rate can be read against
    # something. Cohen's kappa rescales the trace by the agreement two independent
    # assignments with these margins would produce.
    chance = float(Mn.sum(axis=1) @ Mn.sum(axis=0))
    trace = float(np.trace(Mn.to_numpy()))
    print(f"  chance agreement under independent margins: {chance:.4f}")
    print(f"  Cohen's kappa: {(trace - chance) / (1 - chance):.4f}")

    # Modal assignment: give each patent whole to its highest probability sector.
    hi = np.asarray(SECTORS)[A.argmax(axis=1)]
    hu = np.asarray(SECTORS)[B.argmax(axis=1)]
    H = pd.crosstab(pd.Series(hu, name="uspc_route"), pd.Series(hi, name="ipc_route"),
                    normalize=True).reindex(index=SECTORS, columns=SECTORS,
                                            fill_value=0.0)
    print("\nModal assignment agreement matrix")
    print(H.round(4).to_string())
    print(f"\n  modal agreement rate: {np.trace(H.to_numpy()):.4f}")

    # ------------------------------------------------------------- by year
    d = pd.DataFrame({"year": m["year"].to_numpy(),
                      "agree": (A * B).sum(axis=1),
                      "modal": (hi == hu).astype(float)})
    ann = d.groupby("year").agg(prob_agreement=("agree", "mean"),
                                modal_agreement=("modal", "mean"),
                                patents=("agree", "size"))
    ann.to_csv(out_dir / "overlap_annual_agreement.csv")
    print("\nAgreement by year, every fifth year")
    print(ann.loc[ann.index % 5 == 0].round(4).to_string())

    print(f"\nWrote three files to {out_dir}")


if __name__ == "__main__":
    main()
