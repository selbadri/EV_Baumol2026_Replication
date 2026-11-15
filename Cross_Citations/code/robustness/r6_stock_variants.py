"""
==========================================================================================
 An Anatomy of the Baumol Cost Disease -- Elbadri & Vizcaino (2026)
 ------------------------------------------------------------------------------------------
 File:    robustness/r6_stock_variants.py
 Purpose: Patent stocks by sector under eight weightings (count, KPSS value, OECD quality,
          triadic, US inventor), 1950-2023
 ------------------------------------------------------------------------------------------
 Inputs:  output/robustness/kpss_flows.csv, kpss_flows_exbm.csv (r5), quality_sector_flows.csv (r4),
          tpf_sector_flows.csv (r3), ipc_sector_flows.csv (r2)
          input/USPTO_HPDF/annual.csv   USPTO annual grant counts, before 1976
 Outputs: output/robustness/stocks_all_variants.csv   N and ln N by variant, sector and year
          output/robustness/stocks_summary.csv        ln N at selected dates (read by
                                                      Calibration, KPSS rows)
          output/robustness/fig_stocks_headline.pdf, fig_stocks_robustness.pdf
 ------------------------------------------------------------------------------------------
 Notes:   Not in the paper. N_j(t) = (1 - 0.15) N_j(t-1) + F_j(t). Every weight averages one
          within the grant year and every variant is put on the true annual grant total,
          so a variant changes the sector composition and not the number of varieties.
          Before 1976 the aggregate flow is split at each variant's first measured shares.
          Usable windows: KPSS and counts 1976-2023, fwd_cits5 to 2019, fwd_cits7 to 2017,
          renewal 1985-2005, triadic from 1985.
==========================================================================================
"""

import sys
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))   # code/, for paths.py
import paths as PATHS  # noqa: E402

SAVE_FIGURES = SAVE_TABLES = True

SECTORS = ["Goods", "RI_Services", "NRI_Services"]
LAB = {"Goods": "Goods", "RI_Services": "RI Services", "NRI_Services": "NRI Services"}
COL = {"Goods": "black", "RI_Services": "forestgreen", "NRI_Services": "firebrick"}
DELTA = 0.15
Y0, Y1 = 1950, 2023
SPLICE = 1976

# variant: (source file, measure column value, first measured year, last usable year, label)
# The last usable year is set by truncation in the underlying indicator, not by the data end.
VARIANTS = [
    ("value_weighted",  "kpss",    "value_weighted",     1976, 2023,
     "KPSS market value"),
    ("unweighted",      "ipc",     "all",                1976, 2023,
     "Unweighted count"),
    ("fwd_cits5",       "quality", "fwd_cits5_plus1",    1976, 2019,
     "Forward citations, 5 year window"),
    ("fwd_cits7",       "quality", "fwd_cits7_plus1",    1976, 2017,
     "Forward citations, 7 year window"),
    ("family_size",     "quality", "family_size_plus1",  1976, 2023,
     "Patent family size"),
    ("renewal",         "quality", "renewal_plus1",      1985, 2005,
     "Renewal to full term"),
    ("triadic",         "tpf",     "triadic",            1985, 2023,
     "Triadic families only"),
    ("us_inventor",     "ipc",     "us_any",             1976, 2023,
     "US resident inventor only"),
    ("value_weighted_exbm", "kpss_exbm", "value_weighted",  1976, 2023,
     "KPSS value, excl. business methods"),
]

FILES = {"kpss": "kpss_flows.csv", "quality": "quality_sector_flows.csv",
         "tpf": "tpf_sector_flows.csv", "ipc": "ipc_sector_flows.csv"}

# Optional. Written by r5_kpss_flows.py with --drop bm_ids_cpc.csv --tag _exbm, which is
# the headline value weighting applied to the sample with business method patents removed.
OPTIONAL = {"kpss_exbm": "kpss_flows_exbm.csv"}


def find_inputs():
    d = PATHS.OUT_ROBUSTNESS
    found = {k: PATHS.need(d / f) for k, f in FILES.items()}
    for k, f in OPTIONAL.items():
        if (d / f).exists():
            found[k] = d / f
        else:
            print(f"  optional input {f} not present, its variant is skipped")
    return found, PATHS.need(PATHS.HPDF_ANNUAL), d


def stock(flow, delta=DELTA):
    """Perpetual inventory with geometric obsolescence, initialised at zero."""
    out, s = [], 0.0
    for v in flow.to_numpy():
        s = (1.0 - delta) * s + v
        out.append(s)
    return pd.Series(out, index=flow.index)


def main():
    found, f_agg, out_dir = find_inputs()
    print("Inputs")
    for k, p in found.items():
        print(f"  {k:8s} {p}")
    print(f"  aggregate {f_agg}\n")

    src = {k: pd.read_csv(p) for k, p in found.items()}
    F_agg = (pd.read_csv(f_agg).set_index("year")["total_iss"].astype(float))

    # Reference annual flow: true utility grants. From 1976 this is the full IPC sample; before
    # that the USPTO historical file, scaled to meet it at the splice year.
    ipc = src["ipc"]
    col_i = "measure" if "measure" in ipc.columns else "basis"
    ref_post = (ipc[ipc[col_i] == "all"].groupby("year")["patents"].sum().sort_index())
    k = ref_post.loc[SPLICE] / F_agg.loc[SPLICE]
    ref_pre = F_agg.loc[:SPLICE - 1] * k
    REF = pd.concat([ref_pre, ref_post]).sort_index()
    print(f"Reference flow: {REF.index.min()} to {REF.index.max()}, "
          f"pre-{SPLICE} scaled by {k:.4f}\n")

    rows, series = [], {}
    for key, s, measure, first, last, label in VARIANTS:
        if s not in src:
            continue
        d = src[s]
        # the collection scripts name this column `measure` or `basis`
        col = "measure" if "measure" in d.columns else "basis"
        if measure not in set(d[col]):
            print(f"  skipping {key}: no {col} named {measure} in {FILES[s]}, "
                  f"available: {sorted(set(d[col]))}")
            continue
        P = (d[d[col] == measure]
             .pivot(index="year", columns="sector", values="patents")[SECTORS]
             .sort_index())
        P = P.loc[first:last]

        # Every variant is put on the same annual total, the true count of utility grants,
        # so that a variant changes only the composition of varieties across sectors and not
        # how many there are. Without this the quality variants would inherit the coverage
        # ramp of the OECD file, which matches only 1.5 percent of 1976 grants, and the
        # subsample variants would be read as a smaller economy rather than a reweighted one.
        share = P.div(P.sum(axis=1), axis=0)
        pre_years = [y for y in REF.index if Y0 - 30 <= y < share.index.min()]
        pre_share = pd.DataFrame({c: pd.Series(share.iloc[0][c], index=pre_years)
                                  for c in SECTORS})
        share = pd.concat([pre_share, share]).sort_index()
        full = share.mul(REF.reindex(share.index), axis=0).dropna()

        N = pd.DataFrame({c: stock(full[c]) for c in SECTORS})
        N = N.loc[Y0:last]
        series[key] = (N, label, first, last)

        for y in N.index:
            for c in SECTORS:
                rows.append({"variant": key, "label": label, "sector": c, "year": y,
                             "N": N.loc[y, c], "lnN": np.log(N.loc[y, c]),
                             "measured_from": first, "usable_to": last})
        print(f"  {key:14s} {label:34s} {first}-{last}, "
              f"{'measured' if first <= SPLICE else 'measured from ' + str(first)}")

    A = pd.DataFrame(rows)
    if SAVE_TABLES:
        A.to_csv(out_dir / "stocks_all_variants.csv", index=False)

    # ---------------------------------------------------------------- summary table
    marks = [1950, 1960, 1976, 1990, 2000, 2010, 2019, 2023]
    sm = []
    for key, (N, label, first, last) in series.items():
        for c in SECTORS:
            r = {"variant": key, "label": label, "sector": c}
            for y in marks:
                r[str(y)] = round(np.log(N.loc[y, c]), 3) if y in N.index else np.nan
            sm.append(r)
    SM = pd.DataFrame(sm)
    if SAVE_TABLES:
        SM.to_csv(out_dir / "stocks_summary.csv", index=False)
    print("\n\nln stock of varieties by sector, at selected dates\n")
    print(SM.to_string(index=False))

    # ---------------------------------------------------------------- figures
    N, label, first, last = series["value_weighted"]
    fig, ax = plt.subplots(figsize=(6.2, 4.0))
    for c in SECTORS:
        ax.plot(N.index, np.log(N[c]), color=COL[c], lw=1.8, label=LAB[c])
    ax.axvline(SPLICE, color="0.6", lw=0.8, ls=":")
    ax.annotate("measured sector split begins", xy=(SPLICE, ax.get_ylim()[0]),
                xytext=(SPLICE + 1.5, ax.get_ylim()[0] + 0.4), fontsize=7, color="0.35")
    ax.set_xlabel("grant year"); ax.set_ylabel(r"$\ln N_j$")
    ax.set_title(label, fontsize=10)
    ax.legend(fontsize=8, frameon=False, loc="lower right")
    ax.spines[["top", "right"]].set_visible(False)
    fig.tight_layout()
    if SAVE_FIGURES:
        fig.savefig(out_dir / "fig_stocks_headline.pdf")
        fig.savefig(out_dir / "fig_stocks_headline.png", dpi=170)

    keys = [k for k, *_ in VARIANTS if k in series]
    n = len(keys)
    ncol = 4
    nrow = int(np.ceil(n / ncol))
    fig, axes = plt.subplots(nrow, ncol, figsize=(3.1 * ncol, 2.7 * nrow),
                             sharex=True, sharey=True)
    for a, key in zip(axes.ravel(), keys):
        N, label, first, last = series[key]
        for c in SECTORS:
            a.plot(N.index, np.log(N[c]), color=COL[c], lw=1.4)
        if first > SPLICE:
            a.axvspan(Y0, first, color="0.92", zorder=0)
        a.set_title(label, fontsize=8.5)
        a.spines[["top", "right"]].set_visible(False)
    for a in axes.ravel()[n:]:
        a.set_visible(False)
    for a in axes[-1]:
        a.set_xlabel("grant year", fontsize=8)
    for a in axes[:, 0]:
        a.set_ylabel(r"$\ln N_j$", fontsize=9)
    axes.ravel()[0].legend([plt.Line2D([], [], color=COL[c], lw=1.4) for c in SECTORS],
                           [LAB[c] for c in SECTORS], fontsize=7, frameon=False,
                           loc="lower right")
    fig.tight_layout()
    if SAVE_FIGURES:
        fig.savefig(out_dir / "fig_stocks_robustness.pdf")
        fig.savefig(out_dir / "fig_stocks_robustness.png", dpi=170)

    print(f"\nWrote to {out_dir}:")
    print("  stocks_all_variants.csv, stocks_summary.csv")
    print("  fig_stocks_headline.pdf, fig_stocks_robustness.pdf")
    print("\nShaded panels mark variants whose measure does not exist before the year")
    print("shown, so their pre-sample carries the aggregate path at a fixed composition.")


if __name__ == "__main__":
    main()
