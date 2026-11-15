# An Anatomy of the Baumol Cost Disease: Patents and Cross-Sector Citations

**Selim Elbadri & Juan I. Vizcaino (2026)**

## Description

This folder builds the patent data used in *An Anatomy of the Baumol Cost Disease*. A set of Python scripts assigns US patents and patent citations to the paper's three sectors, Goods, research-intensive services (RI-Services) and non-research-intensive services (NRI-Services), and produces:

1. **Patent shares and knowledge stocks**, 1947–2023. These are the sector shares of issued patents in Section 2 and the knowledge stocks of the model (Appendix C).
2. **Cross-sector citation probabilities** δ<sub>jℓ</sub>(t), 1977–2023. These are the citation shares in Figure 3, the knowledge-spillover elasticities of Section 4 and Table D.2 (Appendix D).
3. **Alternative patent stocks** (optional). These are weighted by stock-market value, OECD quality indicators, triadic filing or inventor residence. They are robustness checks and do not appear in the paper.

Patents are assigned to sectors through their IPC subclasses with the Algorithmic Links with Probabilities concordance of Lybbert & Zolas (2014)[^1] to ISIC Rev. 3.1. ISIC classes are then mapped to sectors as in Table D.1. Business-method patents (current CPC subclass G06Q) are excluded throughout.

The two sibling folders read the results directly from `output/`:

| File | Read by |
|------|---------|
| `output/patents/ex_cpc_hist_flow_shares.csv` | `../Empirical/code/10_patents.do` (`PATSH_*`) |
| `output/patents/ex_cpc_hist_stocks.csv` | `../Empirical/code/10b_patent_stocks.do` (`PATSTOCK_*`) |
| `output/citations/spillover_backward1.csv` | `../Calibration/code/paths.jl` (Figure 3, spillover elasticities) |
| `output/robustness/stocks_summary.csv` | `../Calibration/code/paths.jl` (KPSS patent-stock summary) |
| `output/patents/class_stocks.csv` | `../Calibration/code/paths.jl` (only when `PATENT_SERIES` is not the default) |

In turn, step 04 reads the USPTO historical masterfile and the USPC–NAICS crosswalk from `../Empirical/input/`.

## Attribution & License

If you use this code, please cite the paper:

> Elbadri, S. & Vizcaino, J. I. (2026). *An Anatomy of the Baumol Cost Disease*. Working paper.

The code in this folder is released under the MIT License (see `LICENSE`). Third-party data remain under their providers' terms and are not redistributed here, apart from the USPTO annual counts in `input/USPTO_HPDF/`. The ISIC-to-sector concordance in `input/concordances/` is the authors' own.

## Repository Structure

```
Cross_Citations/
├── run_all.py                      Runs every step below, in order
├── requirements.txt                Python packages
├── README.md
├── LICENSE
├── code/
│   ├── paths.py                    All input and output paths
│   ├── 01_business_methods.py      Business-method flag (current CPC G06Q)
│   ├── 02_class_counts.py          Utility patents by IPC subclass and grant year
│   ├── 03_partition_compare.py     ISIC vs NAICS sector partition (Table D.1 notes)
│   ├── 04_hpdf_sector_flows.py     USPC-route sector flows, 1790–2014 (before 1976)
│   ├── 05_patent_stocks.py         Patent shares and knowledge stocks
│   ├── 06_citation_probabilities.py  Cross-sector citation probabilities
│   └── robustness/                 Alternative stocks, not in the paper
│       ├── r1_concordance_overlap.py   IPC route vs USPC route, 1976–2014
│       ├── r2_ipc_sector_flows.py      Flows, all patents and US-resident inventors
│       ├── r3_triadic_flows.py         Flows, triadic families
│       ├── r4_quality_flows.py         Flows, OECD quality indicators
│       ├── r5_kpss_flows.py            Flows, KPSS stock-market value
│       └── r6_stock_variants.py        Stocks under every weighting
├── input/
│   ├── concordances/
│   │   ├── ipc4_to_isic_rev31_4.txt        ALP concordance, IPC4 → ISIC Rev. 3.1 [1], not included
│   │   ├── isic31_4digit_to_sector.csv     ISIC Rev. 3.1 class → sector (Table D.1)
│   │   └── uspc_Class_to_isic_rev31_4.txt  ALP concordance, USPC → ISIC Rev. 3.1 [2] (r1 only), not included
│   ├── USPTO_HPDF/annual.csv       USPTO annual patent counts, 1840–2014 [3]
│   ├── PatentsView/                PatentsView tables [4], not included
│   ├── OECD_Citations_Feb2025/     OECD Citations database [5], not included
│   ├── OECD_TPF_Feb2025/           OECD Triadic Patent Families [6], not included (r3)
│   ├── OECD_Quality_Feb2025/       OECD Patent Quality Indicators [7], not included (r4)
│   └── KPSS/                       KPSS patent values [8], not included (r5)
├── intermediate/                   OECD archives unpacked by step 06 (not tracked)
└── output/
    ├── patents/                    Steps 01–05
    ├── citations/                  Step 06
    ├── robustness/                 code/robustness/
    └── logs/                       One log per step (not tracked)
```

### Steps

| Step | File | Builds |
|------|------|--------|
| 01 | `01_business_methods.py` | Patents whose current CPC contains G06Q (`bm_ids_cpc.csv`). IPC-subclass counts without them, and their share of grants |
| 02 | `02_class_counts.py` | Utility patents by IPC subclass and grant year, 1976–2023. Each patent has weight one, split across its subclasses |
| 03 | `03_partition_compare.py` | Compares the patent-side ISIC sector map with the NAICS partition of Appendix A, class by class |
| 04 | `04_hpdf_sector_flows.py` | Patents by sector on the USPC route, 1790–2014, as in `../Empirical/code/10_patents.do` |
| 05 | `05_patent_stocks.py` | Sector shares and stocks. The paper's series is the `ex_cpc_hist` aggregation: IPC route without business methods, no obsolescence, and before 1976 the change in composition from step 04 (Appendix C). Also writes four alternative aggregations |
| 06 | `06_citation_probabilities.py` | Citations made by USPTO-granted patents, each with weight one, mapped to sectors as M = S′CS and normalised by citing sector (Appendix D, eqs. D.1–D.4). Annual, pooled and alternative-sample probabilities |

Steps 01, 02, 04 and the first part of 06 read the large downloads. They write small tables in `output/`: the IPC-subclass counts (`ipc4_year_counts*.csv`, `uspc_year_counts.csv`), the business-method patent list (`bm_ids_cpc.csv`), the USPC-route flows, and `output/citations/cell_counts_ipc4_exbm.parquet` (citations by citing and cited subclass, year, origin and category). With those tables, everything after them reruns in about a minute (`--quick`, below). They are built from licensed data, so the public package does not include them; it includes the sector-level results read by `../Empirical` and `../Calibration`.

## Outputs and the paper

| Paper item | File | Step |
|------------|------|------|
| Figures 1–2, right panels (sector shares of issued patents) | `output/patents/ex_cpc_hist_flow_shares.csv`, via `../Empirical` (`PATSH_*`). Drawn by `../Calibration/code/motivation_plots.jl` | 05 |
| Figure 3 (citation probabilities by citing year) | `output/citations/spillover_backward1.csv`. Drawn by `../Calibration/code/motivation_plots.jl` | 06 |
| Table 1, initial knowledge stocks N̄ = (0.96, 0.01, 0.03); Goods patent-stock target (1.50 in 2023) | `output/patents/ex_cpc_hist_stocks.csv`, via `../Empirical` (`PATSTOCK_*`) | 05 |
| Section 4.1, spillover elasticities δ<sub>jℓ</sub> (0.955, 0.792, 0.811, …) | Mean of `spillover_backward1.csv` over citing years 1978–2020, each citing sector renormalised (computed in `../Calibration`; also printed at the end of step 06) | 06 |
| Section 2, footnote 9 (citations against independence) | `output/citations/spillover_vs_null.csv` | 06 |
| Appendix C, sector shares 1950/2014/2023 and N̄ = (0.961, 0.008, 0.032) | `output/patents/ex_cpc_hist_flow_shares.csv`, `ex_cpc_hist_stocks.csv` | 05 |
| Appendix C and D, business-method share of grants (0.9–4.8%) | `output/patents/bm_cpc_coverage.csv` | 01 |
| Appendix D, citations per year (4,178 in 1980 … 632,035 in 2023) | `citations` column of `output/citations/spillover_backward1.csv` | 06 |
| Table D.1 and the complete ISIC-to-sector concordance | `input/concordances/isic31_4digit_to_sector.csv`. The two classes discussed in the notes are in `output/patents/partition_mismatches.tex` | 03 |
| Table D.2 | `output/citations/Table_D2.tex` (all columns in `spillover_robustness.csv`) | 06 |

The Appendix C comparison with the USPC route (RI-Services 1.4% of grants in 1950, 15.6% in 2014) comes from `../Empirical/code/10_patents.do`. The other files in `output/` are diagnostics that do not appear in the paper. These include the four alternative aggregations in `class_stocks.csv`, the class decomposition of the services surge, the five-year citation blocks and the two `spillover_by_*.pdf` figures.

## Data availability

The files below are not included. Save them in the folder shown; the top-level README gives step-by-step download instructions. PatentsView (the three main tables), the OECD Citations database and the IPC concordance are needed to rerun this folder; the rest only for the robustness checks.

| Folder | Files | Source |
|--------|-------|--------|
| `input/PatentsView/` | `g_patent.tsv.zip`, `g_ipc_at_issue.tsv.zip`, `g_cpc_current.tsv.zip`. Optional: `g_uspc_at_issue.tsv.zip` (step 02, r1). Robustness only: `g_inventor_disambiguated.tsv.zip`, `g_location_disambiguated.tsv.zip` (r2, r3) | PatentsView final release (12/31/2024), USPTO record on Zenodo[^4], about 2 GB. Keep the files zipped |
| `input/OECD_Citations_Feb2025/` | `202501_USPTO_CITATIONS_1.7z` … `_5.7z`, `202501_USPTO_Equivalent.7z` (or the unpacked `.txt` files) | OECD Citations database, February 2025[^5]. Distributed by the OECD on request. About 3 GB packed and 24 GB unpacked; step 06 unpacks into `intermediate/` |
| `input/OECD_TPF_Feb2025/` | `202501_TPF_USPTO.txt`, `202501_TPF_Core.txt` | OECD Triadic Patent Families database, February 2025[^6] (r3 only) |
| `input/OECD_Quality_Feb2025/` | `202502_OECD_PATENT_QUALITY_USPTO_INDIC.txt` | OECD Patent Quality Indicators database, February 2025[^7] (r4 only) |
| `input/concordances/` | `ipc4_to_isic_rev31_4.txt`; `uspc_Class_to_isic_rev31_4.txt` (r1 only) | ALP concordances of Lybbert & Zolas (2014)[^1] and Goldschlag, Lybbert & Zolas (2016)[^2], version 2209, from https://sites.google.com/site/nikolaszolas/PatentCrosswalk (*IPC Crosswalks* and *USPC Crosswalks*, folders `ISIC_Rev31/`). Posted publicly without a stated license, so not redistributed |
| `input/KPSS/` | `KPSS_2024.csv` | Kogan, Papanikolaou, Seru & Stoffman extended data[^8]. Unzip `KPSS_2024.zip` from the repository (r5 only) |

## Instructions

**Requirements:** Python 3.9 or later (run with 3.9.6) and the packages pinned in `requirements.txt`. The full run reads about 30 GB of text and takes about 25 minutes; the quick run takes about a minute.

1. Open a terminal in this folder and install the packages once:
   ```
   python -m pip install -r requirements.txt
   ```
2. Save the downloads as described above and run:
   ```
   python run_all.py
   python run_all.py --robustness     # also the alternative stocks in code/robustness/
   ```
   Once steps 01, 02, 04 and 06 have written their tables, `python run_all.py --quick` reruns the rest in about a minute.

Its results are included in `output/`, so `../Empirical` and `../Calibration` run without rerunning this folder; if you rerun it, do so before them. Each script can also be run on its own, e.g. `python code/06_citation_probabilities.py --from-cache`. Inputs and outputs are listed in each file's header. Logs are written to `output/logs/`.

## Contact

Selim Elbadri: selim.elbadri@nottingham.ac.uk  
Juan I. Vizcaino: juan.vizcaino@nottingham.ac.uk

## References

[^1]: Lybbert, T. J. & Zolas, N. J. (2014). Getting Patents and Economic Data to Speak to Each Other: An "Algorithmic Links with Probabilities" Approach for Joint Analyses of Patenting and Economic Activity. *Research Policy* 43(3), 530–542. Concordances: https://sites.google.com/site/nikolaszolas/PatentCrosswalk
[^2]: Goldschlag, N., Lybbert, T. J. & Zolas, N. J. (2016). An "Algorithmic Links with Probabilities" Crosswalk for USPC and CPC Patent Classifications. US Census Bureau CES Working Paper 16-15. https://sites.google.com/site/nikolaszolas/PatentCrosswalk
[^3]: Marco, A. C., Carley, M., Jackson, S. & Myers, A. (2015). The USPTO Historical Patent Data Files: Two Centuries of Innovation. https://www.uspto.gov/ip-policy/economic-research/research-datasets/historical-patent-data-files
[^4]: United States Patent and Trademark Office (2025). *Final Release of PatentsView Metadata, Pre-Grant and Granted (12/31/2024)*. Zenodo. https://doi.org/10.5281/zenodo.15783125
[^5]: OECD Citations database, February 2025. Distributed free of charge on request by the OECD Directorate for Science, Technology and Innovation; derived from EPO PATSTAT, Autumn 2024. Method: Webb, C., Dernis, H., Harhoff, D. & Hoisl, K. (2005). Analysing European and International Patent Citations: A Set of EPO Patent Database Building Blocks. OECD STI Working Paper 2005/9.
[^6]: OECD Triadic Patent Families database, February 2025. Distributed free of charge on request by the OECD Directorate for Science, Technology and Innovation; derived from EPO PATSTAT, Autumn 2024.
[^7]: OECD Patent Quality Indicators database, February 2025. Distributed free of charge on request by the OECD Directorate for Science, Technology and Innovation; derived from EPO PATSTAT, Autumn 2024. Indicators: Squicciarini, M., Dernis, H. & Criscuolo, C. (2013). Measuring Patent Quality: Indicators of Technological and Economic Value. OECD STI Working Paper 2013/03. https://doi.org/10.1787/5k4522wkw1r8-en
[^8]: Kogan, L., Papanikolaou, D., Seru, A. & Stoffman, N. (2017). Technological Innovation, Resource Allocation, and Growth. *Quarterly Journal of Economics* 132(2), 665–712. Data: https://github.com/KPSS2017/Technological-Innovation-Resource-Allocation-and-Growth-Extended-Data
