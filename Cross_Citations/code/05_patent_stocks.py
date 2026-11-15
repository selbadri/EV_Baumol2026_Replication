"""
==========================================================================================
 An Anatomy of the Baumol Cost Disease -- Elbadri & Vizcaino (2026)
 ------------------------------------------------------------------------------------------
 File:    05_patent_stocks.py
 Purpose: Sector shares of issued patents and knowledge stocks by sector, 1947-2023
          (Section 2, Figures 1-2 right panels; Table 1, initial knowledge stocks;
          Appendix C), plus alternative aggregations for robustness
 ------------------------------------------------------------------------------------------
 Inputs:  output/patents/ipc4_year_counts.csv                       (step 02)
          output/patents/ipc4_year_counts_exbm.csv                  (step 02)
          output/patents/ipc4_year_counts_excpc.csv                 (step 01)
          output/patents/isic31_4digit_to_sector_harmonised.csv     (step 03)
          output/patents/hpdf_sector_flows.csv                      (step 04)
          input/concordances/ipc4_to_isic_rev31_4.txt
          input/concordances/isic31_4digit_to_sector.csv
          input/USPTO_HPDF/annual.csv          USPTO annual patent counts, 1840-2014
 Outputs: output/patents/
            ex_cpc_hist_flow_shares.csv  annual sector shares of grants, 1947-2023 (Figures
                                         1-2; read by ../Empirical/code/10_patents.do)
            ex_cpc_hist_stocks.csv       knowledge stocks N_j(t), 1950-2023 (Table 1 N-bar,
                                         Goods stock target; read by 10b_patent_stocks.do)
            class_stocks.csv             stocks under every aggregation (read by Calibration
                                         when PATENT_SERIES is not the default)
            class_sector_flows.csv, class_growth.csv, class_contributions.csv,
            tab_class_contributions.tex, fig_class_stocks.pdf, fig_class_decomposition.pdf
                                         diagnostics, not in the paper
 ------------------------------------------------------------------------------------------
 Notes:   Aggregations built from the same class counts:
            baseline      the sector map of Appendix D, all utility patents
            ex_bm         excluding patents with at-issue G06Q or USPC 705 (not time-consistent)
            ex_cpc        excluding patents whose current CPC contains G06Q (step 01)
            harmonised    baseline with the two classes of step 03 moved to RI-Services
            ex_cpc_hist   THE PAPER'S SERIES. ex_cpc, accumulated without obsolescence
                          (N_j(t) = N_j(t-1) + P_j(t)), and before 1976 the 1976 IPC shares
                          scaled by the change in the USPC-route shares (step 04), then
                          renormalised (Appendix C). Levels and the sector partition stay
                          those of the IPC route.
          The others accumulate with geometric obsolescence at 0.15 and hold the 1976 shares
          before 1976. Before 1976 the level is the USPTO annual grant count scaled to each
          aggregation's own 1976 total. The script also decomposes the change in each
          service sector's flow into the contribution of individual IPC subclasses.
==========================================================================================
"""

import sys

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

import paths as P

SAVE_FIGURES = SAVE_TABLES = True

SECTORS = ["Goods", "RI_Services", "NRI_Services"]
LAB = {"Goods": "Goods", "RI_Services": "RI Services", "NRI_Services": "NRI Services"}
COL = {"Goods": "black", "RI_Services": "forestgreen", "NRI_Services": "firebrick"}
RI_DIVISIONS = {"22", "64", "65", "66", "67", "72", "73", "74", "92"}

DELTA = 0.15
Y0, Y1 = 1950, 2023
SPLICE = 1976
SMOOTH = 3                      # centred window on interval endpoints
INTERVALS = [(1976, 1990), (1990, 2000), (2000, 2010), (2010, 2019), (2019, 2023)]
TOP = 10


def optional(path):
    """path if it exists, else None (the run that needs it is skipped)."""
    return path if path.exists() else None


def sector_matrix(f_cw, f_map):
    """Row stochastic IPC subclass by sector matrix, so patent mass is preserved."""
    m = pd.read_csv(f_map, dtype=str)
    m.columns = [c.strip() for c in m.columns]
    to_sector = (m.assign(isic=m["isic_rev31_4"].str.strip().str.zfill(4))
                  .set_index("isic")["Sector"].str.strip()
                  .replace({"Manufacturing": "Goods", "Other_Goods": "Goods"}))

    def fallback(code):
        c = str(code).zfill(4)
        if not c.isdigit():
            return None
        if int(c) < 5000:
            return "Goods"
        return "RI_Services" if c[:2] in RI_DIVISIONS else "NRI_Services"

    cw = pd.read_csv(f_cw, dtype=str)
    cw.columns = [c.strip() for c in cw.columns]
    cw["cls"] = cw["ipc4"].str.strip().str.upper()
    cw["isic"] = cw["isic_rev31_4"].str.strip().str.zfill(4)
    cw["w"] = pd.to_numeric(cw["probability_weight"], errors="coerce").fillna(0.0)
    cw["sector"] = cw["isic"].map(to_sector)
    gap = cw["sector"].isna()
    if gap.any():
        cw.loc[gap, "sector"] = cw.loc[gap, "isic"].map(fallback)
    cw = cw.dropna(subset=["sector"])
    S = cw.groupby(["cls", "sector"])["w"].sum().unstack(fill_value=0.0)
    for s in SECTORS:
        if s not in S.columns:
            S[s] = 0.0
    S = S[SECTORS]
    return S.div(S.sum(axis=1), axis=0)


def class_sector_flows(counts, S):
    """Patents by class and year, spread over sectors. Index year, columns (class, sector)."""
    C = counts.pivot(index="year", columns="ipc4", values="patents").fillna(0.0)
    keep = [c for c in C.columns if c in S.index]
    lost = C.drop(columns=keep).to_numpy().sum() / C.to_numpy().sum()
    if lost > 0.02:
        print(f"  warning: {lost:.4f} of patent mass sits in classes absent from the "
              f"concordance")
    C = C[keep]
    Sm = S.loc[keep]
    out = {s: C.mul(Sm[s], axis=1) for s in SECTORS}
    return C, out


def stock(flow, delta=DELTA):
    """Perpetual inventory with geometric obsolescence."""
    out, s = [], 0.0
    for v in flow.to_numpy():
        s = (1.0 - delta) * s + v
        out.append(s)
    return pd.Series(out, index=flow.index)


def spliced_shares(F, ref, hist=None):
    """Sector shares of the annual flow: measured from 1976, extended backward before that.

    Before the first measured year the sector shares are held at their first-year values,
    unless `hist` (sector flows from another route covering those years) is given. Then the
    first-year shares are scaled by hist's share in year t relative to its first-year share,
    and renormalised, so only the change in composition is borrowed from `hist`."""
    share = F.div(F.sum(axis=1), axis=0)
    y1 = share.index.min()
    pre = [y for y in ref.index if Y0 - 110 <= y < y1]
    if hist is None:
        pre_share = pd.DataFrame({c: pd.Series(share.iloc[0][c], index=pre) for c in SECTORS})
    else:
        h = hist[SECTORS].div(hist[SECTORS].sum(axis=1), axis=0)
        missing = [y for y in pre + [y1] if y not in h.index]
        if missing:
            sys.exit(f"historical sector flows lack years {missing[:5]}...")
        pre_share = h.loc[pre].div(h.loc[y1]).mul(share.loc[y1])
        pre_share = pre_share.div(pre_share.sum(axis=1), axis=0)
    return pd.concat([pre_share, share]).sort_index()


def spliced_stock(F, ref, delta=DELTA, hist=None):
    """Sector flows put on the reference aggregate with spliced_shares, then accumulated."""
    share = spliced_shares(F, ref, hist)
    full = share.mul(ref.reindex(share.index), axis=0).dropna()
    N = pd.DataFrame({c: stock(full[c], delta) for c in SECTORS})
    N.index.name = "year"
    return N.loc[Y0:F.index.max()]


def endpoint(x, y):
    """Centred mean of a series around year y, to take grant-processing noise out."""
    h = SMOOTH // 2
    w = x.loc[max(x.index.min(), y - h):min(x.index.max(), y + h)]
    return float(w.mean())


def main():
    out_dir = P.OUT_PATENTS
    f_cw = P.need(P.ALP_IPC)
    f_map = P.need(P.SECTOR_MAP)
    f_cnt = P.need(out_dir / "ipc4_year_counts.csv")
    f_ann = P.need(P.HPDF_ANNUAL)
    print(f"Counts:      {f_cnt}\nConcordance: {f_cw}\nOutput:      {out_dir}\n")

    S_base = sector_matrix(f_cw, f_map)
    maps = {"baseline": (f_cnt, S_base)}

    f_exbm = optional(out_dir / "ipc4_year_counts_exbm.csv")
    if f_exbm is not None:
        maps["ex_bm"] = (f_exbm, S_base)
    else:
        print("ipc4_year_counts_exbm.csv not found, skipping the exclusion run")

    f_excpc = optional(out_dir / "ipc4_year_counts_excpc.csv")
    if f_excpc is not None:
        maps["ex_cpc"] = (f_excpc, S_base)
    else:
        print("ipc4_year_counts_excpc.csv not found, skipping the current-CPC exclusion")

    # Stock construction per aggregation; anything not listed uses DELTA and constant
    # pre-1976 shares.
    spec = {}
    f_hist = optional(out_dir / "hpdf_sector_flows.csv")
    if f_excpc is not None and f_hist is not None:
        maps["ex_cpc_hist"] = (f_excpc, S_base)
        spec["ex_cpc_hist"] = {"delta": 0.0,
                               "hist": pd.read_csv(f_hist).set_index("year")}
    else:
        print("hpdf_sector_flows.csv or the ex_cpc counts not found, skipping ex_cpc_hist; "
              "run steps 01 and 04 first")

    f_harm = optional(out_dir / "isic31_4digit_to_sector_harmonised.csv")
    if f_harm is not None:
        maps["harmonised"] = (f_cnt, sector_matrix(f_cw, f_harm))
    else:
        print("isic31_4digit_to_sector_harmonised.csv not found, "
              "run step 03 first for the partition run")

    # Reference aggregate flow, as in robustness/r6_stock_variants.py, built from each aggregation's own
    # counts: an exclusion must remove the patents from the level as well as from the shares.
    F_agg = pd.read_csv(f_ann).set_index("year")["total_iss"].astype(float)
    counts0 = pd.read_csv(f_cnt)

    def reference_flow(counts):
        ref_post = counts.groupby("year")["patents"].sum().sort_index()
        k = ref_post.loc[SPLICE] / F_agg.loc[SPLICE]
        return pd.concat([F_agg.loc[:SPLICE - 1] * k, ref_post]).sort_index(), k

    flows, stocks, growth, contrib = [], [], [], []
    for key, (path, S) in maps.items():
        counts = pd.read_csv(path)
        REF, k = reference_flow(counts)
        print(f"{key:11s} reference flow {REF.index.min()} to {REF.index.max()}, "
              f"pre-{SPLICE} scaled by {k:.4f}")
        C, byclass = class_sector_flows(counts, S)
        F = pd.DataFrame({s: byclass[s].sum(axis=1) for s in SECTORS})
        print(f"{key:11s} {F.index.min()}-{F.index.max()}, "
              f"{int(F.loc[F.index.max()].sum()):,} patents in the last year")

        f = F.reset_index().melt(id_vars="year", var_name="sector", value_name="patents")
        f["aggregation"] = key
        flows.append(f)

        N = spliced_stock(F, REF, **spec.get(key, {}))
        if key == "ex_cpc_hist" and SAVE_TABLES:
            # Annual sector shares of grants for the Section 2 figures, on the route the
            # model's stocks use; the stock file starts in 1950, so its differences cannot.
            sh = spliced_shares(F, REF, spec[key]["hist"]).loc[1947:Y1]
            (sh.reset_index().rename(columns={"index": "year"})
               .melt(id_vars="year", var_name="sector", value_name="share")
               .to_csv(out_dir / "ex_cpc_hist_flow_shares.csv", index=False))
        n = N.reset_index().melt(id_vars="year", var_name="sector", value_name="N")
        if key == "ex_cpc_hist" and SAVE_TABLES:
            # The knowledge stocks read by ../Empirical/code/10b_patent_stocks.do
            n.to_csv(out_dir / "ex_cpc_hist_stocks.csv", index=False)
        n["lnN"] = np.log(n["N"])
        n["aggregation"] = key
        stocks.append(n)

        for t0, t1 in INTERVALS:
            if t1 > F.index.max():
                continue
            for s in SECTORS:
                a, b = endpoint(F[s], t0), endpoint(F[s], t1)
                growth.append({"aggregation": key, "sector": s, "t0": t0, "t1": t1,
                               "flow_t0": a, "flow_t1": b,
                               "g_flow": np.log(b / a) / (t1 - t0),
                               "g_stock": (np.log(N.loc[t1, s] / N.loc[t0, s]) / (t1 - t0)
                                           if t1 in N.index else np.nan)})
                if key != "baseline":
                    continue
                base = a
                for c in byclass[s].columns:
                    ca, cb = endpoint(byclass[s][c], t0), endpoint(byclass[s][c], t1)
                    if cb - ca == 0:
                        continue
                    contrib.append({"sector": s, "t0": t0, "t1": t1, "ipc4": c,
                                    "flow_t0": ca, "flow_t1": cb,
                                    "contribution": (cb - ca) / base})

    FL = pd.concat(flows, ignore_index=True)
    ST = pd.concat(stocks, ignore_index=True)
    GR = pd.DataFrame(growth)
    CO = pd.DataFrame(contrib)
    if SAVE_TABLES:
        FL.to_csv(out_dir / "class_sector_flows.csv", index=False)
        ST.to_csv(out_dir / "class_stocks.csv", index=False)
        GR.to_csv(out_dir / "class_growth.csv", index=False)
        CO.to_csv(out_dir / "class_contributions.csv", index=False)

    # ------------------------------------------------------------------ growth table
    print("\n\nAnnualised log growth of the sector flow, by interval\n")
    g = GR.pivot_table(index=["sector", "t0"], columns="aggregation", values="g_flow")
    print(g.round(4).to_string())
    print("\nAnnualised log growth of the stock, by interval\n")
    print(GR.pivot_table(index=["sector", "t0"], columns="aggregation",
                         values="g_stock").round(4).to_string())

    # ------------------------------------------------------------ contributions
    print("\n\nTop contributing IPC subclasses, share of the interval's proportional "
          "change in the sector flow\n")
    lines = []
    for s in ["RI_Services", "NRI_Services"]:
        for t0, t1 in INTERVALS:
            d = CO[(CO.sector == s) & (CO.t0 == t0)]
            if d.empty:
                continue
            tot = d.contribution.sum()
            d = d.reindex(d.contribution.abs().sort_values(ascending=False).index).head(TOP)
            print(f"{LAB[s]}, {t0} to {t1}: total {tot:+.3f}")
            for _, r in d.iterrows():
                print(f"    {r.ipc4}  {r.contribution:+.4f}   "
                      f"{100 * r.contribution / tot:5.1f}% of the change")
                lines.append((LAB[s], f"{t0}-{t1}", r.ipc4, r.contribution,
                              100 * r.contribution / tot))
            gq = CO[(CO.sector == s) & (CO.t0 == t0) & (CO.ipc4 == "G06Q")]
            if len(gq):
                print(f"    G06Q share of the change: "
                      f"{100 * gq.contribution.iloc[0] / tot:.1f}%")
            hh = ((d.contribution / tot) ** 2).sum()
            print(f"    concentration of the top {TOP}: {hh:.3f}\n")

    if SAVE_TABLES:
        with open(out_dir / "tab_class_contributions.tex", "w") as fh:
            fh.write("\\begin{tabular}{lllrr}\n\\toprule\nSector & Interval & IPC & "
                     "Contribution & Share of change \\\\\n\\midrule\n")
            for a, b, c, d_, e in lines:
                fh.write(f"{a} & {b} & {c} & {d_:+.4f} & {e:.1f}\\% \\\\\n")
            fh.write("\\bottomrule\n\\end{tabular}\n")

    # ------------------------------------------------------ the partition mass check
    C0 = counts0.pivot(index="year", columns="ipc4", values="patents").fillna(0.0)
    for c in ["G03D", "G09B", "G06Q"]:
        if c in C0.columns:
            tot = C0.sum(axis=1)
            print(f"{c}: patents and share of all grants")
            for y in [1976, 1985, 1990, 2000, 2010, 2019, 2023]:
                if y in C0.index:
                    print(f"    {y}  {C0.loc[y, c]:8.1f}   {C0.loc[y, c] / tot.loc[y]:.5f}")
            print()

    # ------------------------------------------------------------------ figures
    if SAVE_FIGURES:
        keys = list(maps)
        fig, axes = plt.subplots(1, len(keys), figsize=(3.5 * len(keys), 3.4),
                                 sharex=True, sharey=True, squeeze=False)
        for a, key in zip(axes[0], keys):
            d = ST[ST.aggregation == key]
            for s in SECTORS:
                x = d[d.sector == s]
                a.plot(x.year, x.lnN, color=COL[s], lw=1.5, label=LAB[s])
            a.axvline(SPLICE, color="0.6", lw=0.8, ls=":")
            a.set_title(key, fontsize=9)
            a.set_xlabel("grant year", fontsize=8)
            a.spines[["top", "right"]].set_visible(False)
        axes[0][0].set_ylabel(r"$\ln N_j$", fontsize=9)
        axes[0][0].legend(fontsize=7, frameon=False, loc="lower right")
        fig.tight_layout()
        fig.savefig(out_dir / "fig_class_stocks.pdf")
        fig.savefig(out_dir / "fig_class_stocks.png", dpi=170)

        fig, axes = plt.subplots(1, 2, figsize=(10, 3.8), sharey=True)
        for a, s in zip(axes, ["RI_Services", "NRI_Services"]):
            labs, gq, oth = [], [], []
            for t0, t1 in INTERVALS:
                d = CO[(CO.sector == s) & (CO.t0 == t0)]
                if d.empty:
                    continue
                labs.append(f"{t0}-{t1}")
                q = d.loc[d.ipc4 == "G06Q", "contribution"].sum()
                gq.append(q)
                oth.append(d.contribution.sum() - q)
            x = np.arange(len(labs))
            a.bar(x, oth, color="0.75", label="all other classes")
            a.bar(x, gq, bottom=oth, color="firebrick", label="G06Q")
            a.set_xticks(x); a.set_xticklabels(labs, fontsize=8)
            a.axhline(0, color="black", lw=0.7)
            a.set_title(LAB[s], fontsize=10)
            a.spines[["top", "right"]].set_visible(False)
        axes[0].set_ylabel("proportional change in the sector flow", fontsize=8)
        axes[0].legend(fontsize=8, frameon=False)
        fig.tight_layout()
        fig.savefig(out_dir / "fig_class_decomposition.pdf")
        fig.savefig(out_dir / "fig_class_decomposition.png", dpi=170)

    print(f"Wrote class_sector_flows.csv, class_stocks.csv, class_growth.csv, "
          f"class_contributions.csv,\ntab_class_contributions.tex, "
          f"fig_class_stocks.pdf and fig_class_decomposition.pdf to {out_dir}")


if __name__ == "__main__":
    main()
