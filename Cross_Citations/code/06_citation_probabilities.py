"""
==========================================================================================
 An Anatomy of the Baumol Cost Disease -- Elbadri & Vizcaino (2026)
 ------------------------------------------------------------------------------------------
 File:    06_citation_probabilities.py
 Purpose: Cross-sector citation probabilities delta_jl(t) (Section 2, Figure 3; Section 4,
          knowledge-spillover elasticities; Appendix D, Table D.2)
 ------------------------------------------------------------------------------------------
 Usage:   python code/06_citation_probabilities.py               full run (reads ~25 GB)
          python code/06_citation_probabilities.py --from-cache  start from the cached
                 technology-level citation counts, skipping steps 1-4 (a few seconds)
 ------------------------------------------------------------------------------------------
 Inputs:  input/OECD_Citations_Feb2025/202501_USPTO_CITATIONS_1..5.txt (or .7z)
          input/OECD_Citations_Feb2025/202501_USPTO_Equivalent.txt      (or .7z)
          input/PatentsView/g_ipc_at_issue.tsv.zip
          input/concordances/ipc4_to_isic_rev31_4.txt
          input/concordances/isic31_4digit_to_sector.csv
          output/patents/bm_ids_cpc.csv                       (step 01)
 Outputs: output/citations/cell_counts_ipc4_exbm.parquet      citation counts by IPC4 pair,
                                                               year, origin and category
          output/citations/
            spillover_backward1.csv        annual delta_jl(t) by citing grant year, with the
                                           citation count (Figure 3; read by Calibration)
            spillover_robustness.csv       pooled delta_jl, benchmark and alternative samples
            Table_D2.tex                   Table D.2
            spillover_vs_null.csv          against independent citing (Section 2, fn. 9)
            spillover_matrix.csv           pooled 3x3 matrix
            spillover_parameters.csv       the same in delta_jl notation
            spillover_fixed5.csv           five-year blocks
            spillover_delta_path.csv       annual path: constant, backward and fixed
            spillover_time_summary.csv     range of each delta_jl over the sample
            spillover_citation_shares.csv  all nine cells by year, two normalisations
            spillover_citation_shares_raw.csv  citations made by each citing sector, by year
            spillover_by_destination.pdf, spillover_by_source.pdf   diagnostic figures
 ------------------------------------------------------------------------------------------
 Notes:   Citations are those made by USPTO-granted patents in the OECD Citations database,
          February 2025 edition. Each citation carries total weight one: the citing patent's
          weight is split equally across its distinct IPC4 subclasses, and a cited
          application takes the average technology distribution of its granted US
          equivalents. The IPC4 x IPC4 counts C_t are mapped to sectors with the
          row-stochastic IPC4 x sector matrix S from the Lybbert-Zolas (ALP) concordance,
          M_t = S' C_t S, and each citing sector's row is normalised to sum to one
          (Appendix D, equations D.1-D.4). Mass is asserted after every step.
          Business-method patents (current CPC G06Q, step 01) are dropped on both the citing
          and the cited side. The model uses the mean of the annual probabilities in
          spillover_backward1.csv over citing years 1978-2020 (Calibration/code).
==========================================================================================
"""

import gc
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

import paths as P

FROM_CACHE = "--from-cache" in sys.argv

EXCLUDE_BUSINESS_METHODS = True
CHUNK_SIZE = 2_000_000

YEAR_MIN = 1976            # first grant year with IPC classes at issue in PatentsView
YEAR_MAX = 2025            # every citing year in the extract; the model uses 1978-2020
WINDOW_YEARS = 1           # backward window: citations made in that grant year only
WINDOW_FIXED = 5           # length of the non-overlapping blocks

CITATION_FILES = [f"202501_USPTO_CITATIONS_{i}" for i in range(1, 6)]
EQUIVALENT_FILE = "202501_USPTO_Equivalent"

SECTOR_CLASSES = ["Manufacturing", "Other_Goods", "RI_Services", "NRI_Services"]
SECTORS = ["Goods", "RI_Services", "NRI_Services"]
SECTOR_COLOUR = {"Goods": "black", "RI_Services": "forestgreen",
                 "NRI_Services": "firebrick", "Manufacturing": "black"}
SECTOR_LABEL = {"Goods": "Goods", "RI_Services": "RI Services",
                "NRI_Services": "NRI Services", "Manufacturing": "Manufacturing"}
CITING_NAME = "Citing sector (destination)"
CITED_NAME = "Cited sector (source)"
LABEL = {"Goods": "m", "Manufacturing": "m", "RI_Services": "r", "NRI_Services": "n"}
DELTA_COLS = [f"delta_{LABEL[a]}{LABEL[b]}" for a in SECTORS for b in SECTORS]

OUT = P.OUT_CITATIONS


# =========================================================================================
# Steps 1-4: technology-level citation counts (skipped with --from-cache)
# =========================================================================================

def oecd_text(stem):
    """The extracted text file of one OECD extract, extracting the .7z only if needed."""
    for folder in (P.OECD_CITATIONS, P.INTERMEDIATE):
        for ext in (".txt", ".csv"):
            if (folder / (stem + ext)).exists():
                return folder / (stem + ext)
    archive = P.need(P.OECD_CITATIONS / (stem + ".7z"))
    import py7zr                                  # only needed when extracting
    print(f"  extracting {archive.name} into {P.INTERMEDIATE} ...")
    with py7zr.SevenZipFile(archive, mode="r") as z:
        names = [n for n in z.getnames() if n.lower().endswith((".txt", ".csv"))]
        if len(names) != 1:
            raise ValueError(f"Expected one text file in {archive.name}, found {names}")
        z.extractall(path=P.INTERMEDIATE)
    return P.INTERMEDIATE / names[0]


def read_header(path, sep="|"):
    """Column names on the first line of a delimited text file."""
    with open(path, "r", encoding="utf-8", errors="replace") as handle:
        return [c.strip() for c in handle.readline().rstrip("\n").split(sep)]


def resolve_columns(available, wanted, required=True):
    """Map wanted column names onto the names present, ignoring case (the OECD extracts
    capitalise the same field differently)."""
    lookup = {c.lower(): c for c in available}
    resolved, missing = {}, []
    for name in wanted:
        actual = lookup.get(name.lower())
        if actual is None:
            missing.append(name)
        else:
            resolved[name] = actual
    if required and missing:
        raise KeyError(f"Columns not found: {missing}. The file contains: {sorted(available)}")
    return resolved


def citing_distribution():
    """Step 2. IPC4 distribution of each granted US patent, weights summing to one."""
    usecols = ["patent_id", "section", "ipc_class", "subclass"]
    parts = []
    for chunk in pd.read_csv(P.pv_file("g_ipc_at_issue"), sep="\t", usecols=usecols,
                             dtype=str, on_bad_lines="skip", chunksize=CHUNK_SIZE):
        chunk = chunk.dropna(subset=usecols)
        chunk["ipc4"] = (chunk["section"].str.strip()
                         + chunk["ipc_class"].str.strip().str.zfill(2)
                         + chunk["subclass"].str.strip())
        parts.append(chunk[["patent_id", "ipc4"]].drop_duplicates())
    df_ipc = pd.concat(parts, ignore_index=True).drop_duplicates()
    del parts
    gc.collect()

    # Dropping flagged patents here removes them as citing patents and, through their
    # granted US equivalents, as cited patents.
    if EXCLUDE_BUSINESS_METHODS:
        bm_ids = set(pd.read_csv(P.need(P.OUT_PATENTS / "bm_ids_cpc.csv"),
                                 dtype=str)["patent_id"])
        n_before = df_ipc.patent_id.nunique()
        df_ipc = df_ipc[~df_ipc["patent_id"].isin(bm_ids)]
        print(f"Business-method patents dropped: {n_before - df_ipc.patent_id.nunique():,} "
              f"of {len(bm_ids):,} flagged")

    df_ipc["w"] = 1.0 / df_ipc.groupby("patent_id")["ipc4"].transform("size")
    assert np.allclose(df_ipc.groupby("patent_id")["w"].sum(), 1.0)
    print(f"Patents with at least one IPC code: {df_ipc.patent_id.nunique():,}")
    print(f"Distinct IPC4 codes: {df_ipc.ipc4.nunique():,}")
    return df_ipc


def cited_distribution(path_equiv, df_ipc):
    """Step 3. IPC4 distribution of each cited application, via its granted US equivalents."""
    cols = resolve_columns(read_header(path_equiv), ["Cited_appln_id", "Eqv_Pub_nbr"])
    eq = (pd.read_csv(path_equiv, sep="|", usecols=list(cols.values()), dtype=str)
            .rename(columns={v: k for k, v in cols.items()}))
    eq["Cited_appln_id"] = pd.to_numeric(eq["Cited_appln_id"], errors="coerce")
    eq = eq.dropna(subset=["Cited_appln_id", "Eqv_Pub_nbr"])
    eq["Cited_appln_id"] = eq["Cited_appln_id"].astype("int64")
    eq["patent_id"] = eq["Eqv_Pub_nbr"].str.removeprefix("US")

    # Grants only: pre-grant publications are US + four-digit year + six digits.
    is_grant = ~eq["patent_id"].str.match(r"^(19|20)\d{2}\d{6}$", na=False)
    eq = eq.loc[is_grant, ["Cited_appln_id", "patent_id"]].drop_duplicates()
    eq["w_eqv"] = 1.0 / eq.groupby("Cited_appln_id")["patent_id"].transform("size")

    dist = (eq.merge(df_ipc, on="patent_id", how="inner")
              .assign(w=lambda d: d["w_eqv"] * d["w"])
              .groupby(["Cited_appln_id", "ipc4"], as_index=False)["w"].sum())
    mass = dist.groupby("Cited_appln_id")["w"].sum()
    dist = dist.merge(mass.rename("total"), left_on="Cited_appln_id", right_index=True)
    dist["w"] = dist["w"] / dist["total"]
    dist = dist[["Cited_appln_id", "ipc4", "w"]]
    assert np.allclose(dist.groupby("Cited_appln_id")["w"].sum(), 1.0)
    print(f"Cited applications with a technology distribution: "
          f"{dist.Cited_appln_id.nunique():,}")
    return dist


def citation_cells(paths_citations, df_ipc, cited_dist):
    """Step 4. Weighted citations by citing IPC4, cited IPC4, citing grant year, citation
    origin (SEA examiner, APP applicant) and category (X, Y, ...)."""
    req = ["Citing_pub_nbr", "Citing_pub_date", "Cited_appln_id"]
    opt = ["Citn_origin", "Citn_category"]
    acc, n_rows, n_kept = [], 0, 0.0
    for path in paths_citations:
        header = read_header(path)
        c_req = resolve_columns(header, req, required=True)
        c_opt = resolve_columns(header, opt, required=False)
        rename = {v: k for k, v in {**c_req, **c_opt}.items()}
        for chunk in pd.read_csv(path, sep="|", usecols=list(c_req.values()) + list(c_opt.values()),
                                 dtype=str, chunksize=CHUNK_SIZE, on_bad_lines="skip"):
            n_rows += len(chunk)
            chunk = chunk.rename(columns=rename)
            for name in opt:
                if name not in chunk.columns:
                    chunk[name] = pd.NA
            chunk["Cited_appln_id"] = pd.to_numeric(chunk["Cited_appln_id"], errors="coerce")
            chunk = chunk.dropna(subset=["Cited_appln_id", "Citing_pub_nbr"])
            if chunk.empty:
                continue
            chunk["Cited_appln_id"] = chunk["Cited_appln_id"].astype("int64")
            chunk["patent_id"] = chunk["Citing_pub_nbr"].str.removeprefix("US")
            chunk["year"] = pd.to_numeric(chunk["Citing_pub_date"].str[:4], errors="coerce")
            chunk = chunk[chunk["year"].between(YEAR_MIN, YEAR_MAX)]
            if chunk.empty:
                continue
            chunk["year"] = chunk["year"].astype("int16")
            chunk = chunk.drop_duplicates(subset=["Citing_pub_nbr", "Cited_appln_id"])

            merged = (chunk.merge(df_ipc.rename(columns={"ipc4": "ipc4_citing", "w": "w_citing"}),
                                  on="patent_id", how="inner")
                           .merge(cited_dist.rename(columns={"ipc4": "ipc4_cited", "w": "w_cited"}),
                                  on="Cited_appln_id", how="inner"))
            if merged.empty:
                continue
            merged["w"] = merged["w_citing"] * merged["w_cited"]
            n_kept += float(merged["w"].sum())
            acc.append(merged.groupby(["ipc4_citing", "ipc4_cited", "year", "Citn_origin",
                                       "Citn_category"], dropna=False, as_index=False)["w"].sum())
            del merged
            gc.collect()
        print(f"Finished {path.name}: running weighted citation total {n_kept:,.0f}")

    cells = (pd.concat(acc, ignore_index=True)
               .groupby(["ipc4_citing", "ipc4_cited", "year", "Citn_origin", "Citn_category"],
                        dropna=False, as_index=False)["w"].sum())
    print(f"\nCitation rows read: {n_rows:,}")
    print(f"Citations retained: {cells['w'].sum():,.0f} "
          f"(share of rows {n_kept / n_rows:.3f})")
    return cells


def build_cells():
    """Steps 1-4, writing the cache."""
    print("Step 1. OECD extracts")
    path_equiv = oecd_text(EQUIVALENT_FILE)
    paths_citations = [oecd_text(s) for s in CITATION_FILES]
    for p in [path_equiv] + paths_citations:
        print(f"  {p}")
    print("\nStep 2. Technology distribution of each US patent")
    df_ipc = citing_distribution()
    print("\nStep 3. Technology distribution of each cited application")
    cited = cited_distribution(path_equiv, df_ipc)
    print("\nStep 4. Technology-by-technology citation counts")
    cells = citation_cells(paths_citations, df_ipc, cited)
    P.CITATION_CELLS.parent.mkdir(parents=True, exist_ok=True)
    cells.to_parquet(P.CITATION_CELLS, index=False)
    print(f"Wrote {P.CITATION_CELLS}")
    return cells


# =========================================================================================
# Steps 5-11: from technologies to sectors
# =========================================================================================

def sector_matrix():
    """Step 5. Row-stochastic IPC4 x industry-class matrix S from the ALP concordance."""
    alp = (pd.read_csv(P.need(P.ALP_IPC))
             .rename(columns={"isic_rev31_4": "isic4", "probability_weight": "p"}))
    assert np.allclose(alp.groupby("ipc4")["p"].sum(), 1.0)
    sector_map = pd.read_csv(P.need(P.SECTOR_MAP))
    alp = alp.merge(sector_map[["isic_rev31_4", "Sector"]],
                    left_on="isic4", right_on="isic_rev31_4", how="left")
    assert alp["Sector"].notna().all(), "Some ISIC classes are not assigned to a sector"
    S = (alp.groupby(["ipc4", "Sector"])["p"].sum()
            .unstack(fill_value=0.0)
            .reindex(columns=SECTOR_CLASSES, fill_value=0.0))
    assert np.allclose(S.sum(axis=1), 1.0)
    print(f"IPC4 codes in the concordance: {len(S):,}")
    return S


def class_matrix(cells, S):
    """M = S' C S across the four industry classes. Codes outside the concordance drop."""
    codes = S.index
    cells = cells[cells.ipc4_citing.isin(codes) & cells.ipc4_cited.isin(codes)]
    C = (cells.pivot_table(index="ipc4_citing", columns="ipc4_cited", values="w",
                           aggfunc="sum", fill_value=0.0)
              .reindex(index=codes, columns=codes, fill_value=0.0))
    # numpy 2 with Apple Accelerate can warn spuriously on matmul; the result is checked
    with np.errstate(all="ignore"):
        M = S.to_numpy().T @ C.to_numpy() @ S.to_numpy()
    assert np.isfinite(M).all()
    return pd.DataFrame(M, index=SECTOR_CLASSES, columns=SECTOR_CLASSES)


def combine_goods(M4):
    """Manufacturing and the other Goods industries form one Goods sector."""
    M = M4.copy()
    M.loc["Goods"] = M.loc["Manufacturing"] + M.loc["Other_Goods"]
    M["Goods"] = M["Manufacturing"] + M["Other_Goods"]
    return M.loc[SECTORS, SECTORS]


def restrict_to_manufacturing(M4):
    """Drop the Goods industries outside manufacturing and condition on what remains."""
    keep = ["Manufacturing", "RI_Services", "NRI_Services"]
    return M4.loc[keep, keep]


def intensities(M):
    """Normalise each citing sector to one. Columns are citing, rows cited."""
    d = M.T.div(M.T.sum(axis=0), axis=1)
    d.index.name, d.columns.name = CITED_NAME, CITING_NAME
    return d


def parameter_table(d):
    """The intensities in the model's delta_jl notation."""
    cols = list(d.columns)
    return pd.DataFrame([{"parameter": f"delta_{LABEL[citing]}{LABEL[cited]}",
                          "source (cited)": SECTOR_LABEL[cited],
                          "destination (citing)": SECTOR_LABEL[citing],
                          "value": round(d.loc[cited, citing], 4)}
                         for citing in cols for cited in cols])


def make_windows(kind):
    if kind == "backward":
        return [(t - WINDOW_YEARS + 1, t)
                for t in range(YEAR_MIN + WINDOW_YEARS - 1, YEAR_MAX + 1)]
    if kind == "fixed":
        return [(y, min(y + WINDOW_FIXED - 1, YEAR_MAX))
                for y in range(YEAR_MIN, YEAR_MAX + 1, WINDOW_FIXED)]
    raise ValueError(kind)


def window_intensities(cells, S, windows, reduce_fn=combine_goods):
    rows = []
    for lo, hi in windows:
        sub = cells[cells.year.between(lo, hi)]
        if sub.empty or sub["w"].sum() == 0:
            print(f"Skipping window {lo}-{hi}: no citations")
            continue
        m4 = class_matrix(sub, S)
        d = intensities(reduce_fn(m4))
        if not np.isfinite(d.to_numpy()).all():
            print(f"Skipping window {lo}-{hi}: at least one sector has no citations")
            continue
        rec = {"start": lo, "end": hi, "citations": m4.to_numpy().sum()}
        for citing in d.columns:
            for cited in d.index:
                rec[f"delta_{LABEL[citing]}{LABEL[cited]}"] = d.loc[cited, citing]
        rows.append(rec)
    return pd.DataFrame(rows)


TABLE_D2_COLUMNS = [("Benchmark", "Benchmark"), ("1976-1990", "1976--1990"),
                    ("1991-2005", "1991--2005"), ("2006-2020", "2006--2020"),
                    ("Examiner citations", "Examiner"), ("Applicant citations", "Applicant"),
                    ("X or Y category", "X or Y")]
TABLE_D2_ROWS = [("delta_mm", r"$\delta_{m,m}$"), ("delta_mr", r"$\delta_{m,rs}$"),
                 ("delta_mn", r"$\delta_{m,ns}$"), ("delta_rm", r"$\delta_{rs,m}$"),
                 ("delta_rr", r"$\delta_{rs,rs}$"), ("delta_rn", r"$\delta_{rs,ns}$"),
                 ("delta_nm", r"$\delta_{ns,m}$"), ("delta_nr", r"$\delta_{ns,rs}$"),
                 ("delta_nn", r"$\delta_{ns,ns}$")]


def table_d2(robustness):
    """Table D.2 as a booktabs tabular."""
    lines = [r"\begin{tabular}{l" + "c" * len(TABLE_D2_COLUMNS) + "}", r"\toprule",
             "Parameter & " + " & ".join(h for _, h in TABLE_D2_COLUMNS) + r" \\",
             r"\midrule"]
    for key, tex in TABLE_D2_ROWS:
        vals = " & ".join(f"{robustness.loc[key, c]:.3f}" for c, _ in TABLE_D2_COLUMNS)
        lines.append(f"{tex} & {vals} " + r"\\")
    lines += [r"\bottomrule", r"\end{tabular}", ""]
    return "\n".join(lines)


def plot_intensities(backward, delta, by_source=False):
    """Annual intensities against the pooled values, by citing or by cited sector."""
    fig, axes = plt.subplots(1, 3, figsize=(13, 4.2), sharey=True)
    for ax, panel in zip(axes, SECTORS):
        for other in SECTORS:
            citing, cited = (other, panel) if by_source else (panel, other)
            col = f"delta_{LABEL[citing]}{LABEL[cited]}"
            colour = SECTOR_COLOUR[other]
            ax.plot(backward["end"], backward[col], color=colour, linewidth=1.6,
                    label=SECTOR_LABEL[other])
            ax.axhline(delta.loc[cited, citing], color=colour, linestyle="--", linewidth=0.8)
        ax.set_title(SECTOR_LABEL[panel])
        ax.grid(alpha=0.3)
    axes[1].set_xlabel("Grant year of the citing patent")
    axes[0].set_ylabel("Share of citations made")
    axes[0].legend(frameon=False, title=(CITING_NAME if by_source else CITED_NAME))
    fig.tight_layout()
    return fig


def main():
    OUT.mkdir(parents=True, exist_ok=True)

    if FROM_CACHE:
        cells = pd.read_parquet(P.need(P.CITATION_CELLS))
        print(f"Read cached citation counts from {P.CITATION_CELLS}")
    else:
        cells = build_cells()
    print(f"Citations: {cells['w'].sum():,.0f}, technology pairs: {len(cells):,}")

    print("\nStep 5. Technology-to-sector concordance")
    S = sector_matrix()

    print("\nStep 6. Sector-by-sector citation matrix")
    in_conc = cells.ipc4_citing.isin(S.index) & cells.ipc4_cited.isin(S.index)
    print(f"Share of citation mass with both codes in the concordance: "
          f"{cells.loc[in_conc, 'w'].sum() / cells['w'].sum():.4f}")
    M4 = class_matrix(cells, S)
    M = combine_goods(M4)

    print("\nStep 7. Pooled intensities")
    delta = intensities(M)
    assert np.allclose(delta.sum(axis=0), 1.0)
    print(delta.round(4).to_string())
    parameter_table(delta).to_csv(OUT / "spillover_parameters.csv", index=False)
    delta.to_csv(OUT / "spillover_matrix.csv")

    print("\nStep 8. Against independent citing (no sorting)")
    codes = S.index
    C = (cells[in_conc].pivot_table(index="ipc4_citing", columns="ipc4_cited", values="w",
                                    aggfunc="sum", fill_value=0.0)
                       .reindex(index=codes, columns=codes, fill_value=0.0).to_numpy())
    with np.errstate(all="ignore"):           # see class_matrix
        C_null = C.sum(axis=1, keepdims=True) @ C.sum(axis=0, keepdims=True) / C.sum()
        M_null = S.to_numpy().T @ C_null @ S.to_numpy()
    assert np.isfinite(M_null).all()
    M4_null = pd.DataFrame(M_null, index=SECTOR_CLASSES, columns=SECTOR_CLASSES)
    delta_null = intensities(combine_goods(M4_null))
    comparison = pd.concat({"Estimated": delta, "No sorting": delta_null,
                            "Difference": delta - delta_null}, axis=1).round(3)
    comparison.to_csv(OUT / "spillover_vs_null.csv")
    for citing in SECTORS:
        gap = float(np.abs(delta[citing] - delta_null[citing]).sum() / 2)
        print(f"  total variation distance, {citing} citing: {gap:.3f}")

    print("\nStep 9. Alternative samples (Table D.2)")
    variants = {"Benchmark": (cells, combine_goods)}
    for lo, hi, name in [(1976, 1990, "1976-1990"), (1991, 2005, "1991-2005"),
                         (2006, 2020, "2006-2020")]:
        variants[name] = (cells[cells.year.between(lo, hi)], combine_goods)
    variants["Examiner citations"] = (cells[cells.Citn_origin.eq("SEA")], combine_goods)
    variants["Applicant citations"] = (cells[cells.Citn_origin.eq("APP")], combine_goods)
    variants["X or Y category"] = (
        cells[cells.Citn_category.fillna("").str.contains("X|Y", regex=True)], combine_goods)
    variants["Manufacturing only"] = (cells, restrict_to_manufacturing)
    rows = []
    for name, (sub, reduce_fn) in variants.items():
        if sub.empty or sub["w"].sum() == 0:
            print(f"Skipping {name}: no citations in this variant")
            continue
        d = intensities(reduce_fn(class_matrix(sub, S)))
        for r in parameter_table(d).itertuples():
            rows.append({"variant": name, "parameter": r.parameter, "value": r.value})
    robustness = (pd.DataFrame(rows)
                    .pivot(index="parameter", columns="variant", values="value")
                    .reindex(columns=list(variants)).round(3))
    robustness.to_csv(OUT / "spillover_robustness.csv")
    (OUT / "Table_D2.tex").write_text(table_d2(robustness), encoding="utf-8")
    print(robustness.to_string())

    print("\nStep 10. Intensities over time")
    backward = window_intensities(cells, S, make_windows("backward"))
    fixed = window_intensities(cells, S, make_windows("fixed"))
    for name, frame in [("backward", backward), ("fixed", fixed)]:
        for citing in SECTORS:
            cols = [f"delta_{LABEL[citing]}{LABEL[c]}" for c in SECTORS]
            assert np.allclose(frame[cols].sum(axis=1), 1.0), f"{name}: {citing}"
    backward.to_csv(OUT / f"spillover_backward{WINDOW_YEARS}.csv", index=False)
    fixed.to_csv(OUT / f"spillover_fixed{WINDOW_FIXED}.csv", index=False)
    print(f"  {len(backward)} annual windows, {backward.end.min()}-{backward.end.max()}")

    years = np.arange(YEAR_MIN, YEAR_MAX + 1)
    constant = pd.DataFrame({"year": years, "method": "constant"})
    for citing in SECTORS:
        for cited in SECTORS:
            constant[f"delta_{LABEL[citing]}{LABEL[cited]}"] = delta.loc[cited, citing]
    back = pd.DataFrame({"year": years, "method": "backward"})
    by_end = backward.set_index("end")
    for col in DELTA_COLS:
        back[col] = [by_end[col].get(y, by_end[col].iloc[0]) for y in years]
    step = pd.DataFrame({"year": years, "method": "fixed"})
    block_of = {y: i for i, (lo, hi) in enumerate(make_windows("fixed"))
                for y in range(lo, hi + 1)}
    for col in DELTA_COLS:
        step[col] = [fixed[col].iloc[block_of[y]] for y in years]
    delta_path = pd.concat([constant, back, step], ignore_index=True)
    for citing in SECTORS:
        cols = [f"delta_{LABEL[citing]}{LABEL[c]}" for c in SECTORS]
        assert np.allclose(delta_path[cols].sum(axis=1), 1.0)
    delta_path.to_csv(OUT / "spillover_delta_path.csv", index=False)

    for by_source, name in [(False, "spillover_by_destination.pdf"),
                            (True, "spillover_by_source.pdf")]:
        fig = plot_intensities(backward, delta, by_source)
        fig.savefig(OUT / name, bbox_inches="tight")
        plt.close(fig)

    # All nine cells by year; share_of_citing_sector reproduces spillover_backward1.csv.
    share_rows = []
    for lo, hi in make_windows("backward"):
        window = cells[cells.year.between(lo, hi)]
        if window.empty or window["w"].sum() == 0:
            continue
        matrix = combine_goods(class_matrix(window, S))      # rows citing, columns cited
        total = float(matrix.to_numpy().sum())
        citing_totals = matrix.sum(axis=1)
        if total <= 0 or (citing_totals <= 0).any():
            continue
        for citing in SECTORS:
            for cited in SECTORS:
                cell = float(matrix.loc[citing, cited])
                share_rows.append({
                    "year": hi, "window_start": lo, "window_end": hi,
                    "citing_sector": citing, "citing_sector_label": SECTOR_LABEL[citing],
                    "cited_sector": cited, "cited_sector_label": SECTOR_LABEL[cited],
                    "citations": cell,
                    "share_of_citing_sector": cell / float(citing_totals.loc[citing]),
                    "share_of_all_citations": cell / total,
                    "citing_sector_citations": float(citing_totals.loc[citing]),
                    "citing_sector_share_of_total": float(citing_totals.loc[citing]) / total,
                })
    shares = pd.DataFrame(share_rows)
    assert np.allclose(shares.groupby(["year", "citing_sector"])["share_of_citing_sector"].sum(), 1.0)
    assert np.allclose(shares.groupby("year")["share_of_all_citations"].sum(), 1.0)
    check = shares.pivot_table(index="year", columns=["citing_sector", "cited_sector"],
                               values="share_of_citing_sector")
    for citing in SECTORS:
        for cited in SECTORS:
            col = f"delta_{LABEL[citing]}{LABEL[cited]}"
            assert np.allclose(check[(citing, cited)].to_numpy(), by_end[col].to_numpy()), col
    shares.to_csv(OUT / "spillover_citation_shares.csv", index=False)
    raw = (shares.groupby(["year", "window_start", "window_end", "citing_sector",
                           "citing_sector_label"], as_index=False)["citations"].sum()
                 .rename(columns={"citing_sector": "sector", "citing_sector_label": "sector_label"}))
    raw["share"] = raw["citations"] / raw.groupby("year")["citations"].transform("sum")
    raw.to_csv(OUT / "spillover_citation_shares_raw.csv", index=False)

    summary = pd.DataFrame(index=DELTA_COLS)
    summary["constant"] = [delta.loc[cited, citing] for citing in SECTORS for cited in SECTORS]
    summary["backward minimum"] = backward[DELTA_COLS].min()
    summary["backward maximum"] = backward[DELTA_COLS].max()
    summary["first window"] = backward[DELTA_COLS].iloc[0]
    summary["last window"] = backward[DELTA_COLS].iloc[-1]
    summary["change"] = summary["last window"] - summary["first window"]
    summary.to_csv(OUT / "spillover_time_summary.csv")

    print("\nStep 11. Values used in the model: mean over citing years 1978-2020, "
          "each citing sector renormalised to one")
    b = backward[backward.start.between(1978, 2020)][DELTA_COLS].mean()
    for s in "mrn":
        cols = [f"delta_{s}{x}" for x in "mrn"]
        print(f"  {s}: " + "  ".join(f"{v:.3f}" for v in b[cols] / b[cols].sum()))

    print(f"\nWrote {OUT}")


if __name__ == "__main__":
    main()
