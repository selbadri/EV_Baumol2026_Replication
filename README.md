# An Anatomy of the Baumol Cost Disease: Replication Package

**Selim Elbadri & Juan I. Vizcaino (2026)**

## Overview

This package reproduces every data series, figure, table and quantitative result of *An Anatomy of the Baumol Cost Disease*. The paper splits the US private economy into three sectors: Goods, research-intensive services (RI-Services) and non-research-intensive services (NRI-Services).

The package has three parts. Each part is a self-contained folder with its own README, LICENSE and one script that runs everything in it.

| Folder | Language | What it does | Main output |
|--------|----------|--------------|-------------|
| [`Cross_Citations/`](Cross_Citations/README.md) | Python | Assigns US patents and patent citations to the three sectors. Builds the patent shares and knowledge stocks, the cross-sector citation probabilities and Table D.2 | `output/patents/`, `output/citations/` |
| [`Empirical/`](Empirical/README.md) | Stata | Builds the annual sector data for 1947–2023: labor and value-added shares, prices, productivity, capital, markups, patents and R&D workers | `output/dataset/data_moments.xlsx` |
| [`Calibration/`](Calibration/README.md) | Julia | Draws the data figures, calibrates the model, and runs the quantitative exercises and robustness checks | `output/figures/`, `output/calibration/`, `output/counterfactuals/`, `output/robustness/` |

```
Cross_Citations  ──►  Empirical  ──►  Calibration
   patent shares        data_moments.xlsx     figures, Table 2,
   and stocks                                 exercises, Tables J.1, H.2
        └──── citation probabilities, stock summary ────►
```

- `Empirical` reads the patent shares and stocks from `../Cross_Citations/output/patents/`.
- `Calibration` reads `../Empirical/output/dataset/data_moments.xlsx`, the NSF R&D benchmarks in `../Empirical/input/manual/NSF/`, and the citation and stock files in `../Cross_Citations/output/`.
- `Cross_Citations` step 04 reads the USPTO historical masterfile and the USPC–NAICS crosswalk from `../Empirical/input/`.

All paths are relative, so the three folders must stay side by side. Nothing needs to be edited to run the package on another machine. The outputs of every step are included, so each part can also be run on its own from the stored results of the part before it.

## Data Availability and Provenance Statements

### Statement about rights

- [x] We certify that the authors of the manuscript have legitimate access to and permission to use the data used in this manuscript.
- [x] We certify that the authors of the manuscript have documented permission to redistribute/publish the data contained within this replication package. Data that cannot be redistributed are listed below as not included, with instructions to obtain them.

### Summary of availability

- [ ] All data **are** publicly available.
- [x] Some data **cannot be made** publicly available.
- [ ] **No data can be made** publicly available.

Six inputs used for the paper's results are not included: the IPUMS USA extract and the IPUMS occupation crosswalk (IPUMS terms of use do not allow redistribution), the Caunedo, Jaume & Keller (2023) quality-adjusted equipment prices (distributed in their replication package), the U.S. KLEMS file (its terms do not state that it may be redistributed), the USPTO historical masterfile (public domain, but at 758 MB too large for GitHub), and the ALP USPC–NAICS patent crosswalk (posted by its authors without a stated license). All six are free to download. The patent data (PatentsView, OECD and KPSS) are not included either, nor are the patent-class tables built from them. They are needed only to rerun `Cross_Citations`. Its sector-level results, which `Empirical` and `Calibration` read, are included in `Cross_Citations/output/`, so the other two parts run without them.

### Data sources

**Included in the package.**

| Data | File(s) | Folder | Provider and access | License / terms |
|------|---------|--------|---------------------|-----------------|
| US industry production accounts, 1947–2023 | `EV_production_accounts_1947to2023.xlsx` | `Empirical/input/` | Elbadri & Vizcaino (2025)[^1], Zenodo | See the Zenodo record |
| BEA fixed assets, detailed estimates by industry and asset | `detailnonres_stk1.xlsx`, `_stk2.xlsx`, `_inv2.xlsx` | `Empirical/input/` | BEA[^2] | US government work, public domain |
| BLS capital detail (rental prices, stocks, deflators) | `rental-prices-major-industries.xlsx` | `Empirical/input/` | BLS[^3] | US government work, public domain |
| Penn World Table 11.0 | `pwt110.dta` | `Empirical/input/` | Feenstra, Inklaar & Timmer (2015)[^5] | CC BY 4.0 |
| Jordà–Schularick–Taylor Macrohistory Database, R6 | `JSTdatasetR6.dta` | `Empirical/input/manual/JST/` | Jordà, Schularick & Taylor (2017)[^6] | CC BY-NC-SA 4.0 |
| NSCG 2023 public-use files and Table 3-1 | `NSCG_2023/college_grads_2023.zip`, `manual/NSCG/` | `Empirical/input/` | NCSES[^11] | US government work, public domain |
| O*NET STEM occupations | `All_STEM_Occupations.xlsx` | `Empirical/input/` | O*NET Resource Center[^12] | CC BY 4.0 |
| NSF R&D employment benchmarks and FRED USPRIV | `manual/NSF/` | `Empirical/input/` | NSF / NCSES; FRED (BLS CES) | Public domain |
| USPTO annual patent counts, 1840–2014 | `input/USPTO_HPDF/annual.csv` | `Cross_Citations/` | Marco et al. (2015)[^8] | US government work, public domain |
| Concordances and assignments built by the authors | `input/manual/*.xlsx`, `input/concordances/isic31_4digit_to_sector.csv` | `Empirical/`, `Cross_Citations/` | This paper | MIT (with the code) |

**Not included, needed for the paper's results.**

| Data | Save as | Needed by | How to obtain |
|------|---------|-----------|---------------|
| IPUMS USA extract | `Empirical/input/usa_00026.dta` (4.5 GB) | `Empirical` step 11 | Free registration at IPUMS USA[^10]; build the extract described below |
| IPUMS OCCSOC crosswalk, 2000 onward | `Empirical/input/occ_occsoc_crosswalk_2000_onward.xlsx` | `Empirical` step 11 | Free download from IPUMS USA[^10] (`occsoc_crosswalk_2000_onward.xlsx`, renamed; see below) |
| U.S. KLEMS, March 2017 release | `Empirical/input/usa_wk_mar_2017.xlsx` | `Empirical` steps 06, 07b, 08 | Jorgenson, Ho & Samuels (2017)[^4]. Free download from World KLEMS, analytical data, https://www.worldklems.net/data/analytical/ (file `usa_wk_mar_2017.xlsx`) |
| Quality-adjusted equipment and software prices | `Empirical/input/manual/CJK/Prices_by_asset_type_with_software_091720.dta` | `Empirical` step 08 | Caunedo, Jaume & Keller (2023)[^7], replication package, Zenodo record 7779807, file `empirics/rawdata/PriceKquality/Prices_by_asset_type_with_software_091720.dta` |
| USPTO Historical Patent Data Files, patent masterfile | `Empirical/input/historical_masterfile.dta` (758 MB) | `Empirical` step 10, `Cross_Citations` step 04 | Marco et al. (2015)[^8]. US government work, public domain, but too large for GitHub. Free download from the USPTO (step 5 below) |
| ALP USPC–NAICS 2007 crosswalk | `Empirical/input/ALP_USPC/NAICS_2007/uspc_Class_to_naics07_2.txt` | `Empirical` step 10, `Cross_Citations` step 04 | Goldschlag, Lybbert & Zolas (2016)[^9]. Posted publicly by the authors without a stated license, so not redistributed. Free download (step 6 below) |

*IPUMS USA extract.* At https://usa.ipums.org, create an extract with:

- **Samples (10):** 1940 1%, 1950 1%, 1960 1%, 1970 1% state fm1, 1980 5% state, 1990 5% state, 2000 5%, 2010 ACS, 2020 ACS (with the experimental weights IPUMS supplies for the COVID-19 period), 2023 ACS.
- **Variables:** `SEX`, `AGE`, `EDUC`, `EMPSTAT`, `OCC1990`, `OCCSOC`, `IND1990`, `INDNAICS` (the detailed versions `EDUCD` and `EMPSTATD` are added automatically), plus the household and person variables IPUMS preselects (`YEAR`, `SAMPLE`, `SERIAL`, `CBSERIAL`, `HHWT`, `CLUSTER`, `STRATA`, `GQ`, `PERNUM`, `PERWT`).
- **Cases:** all persons, no case selection (54,138,391 records).
- **Format:** `.dat` (fixed-width text), rectangular (person), about 5.9 GB, with the Stata command file. Run the IPUMS command file in the folder of the `.dat` file and save the result as `usa_00026.dta`.

Step 11 keeps employed persons aged 21–65 in private industries from 1950 on, so the 1940 sample and the unused variables do not affect the results.

The USPTO historical masterfile is public domain but too large for GitHub, so it is not included; step 5 below explains how to download it from the USPTO Historical Patent Data Files page[^8]. The code uses three variables: `disp_ty` (disposal type), `disp_dt` (disposal date) and `uspc` (USPC class).

**Not included, needed only to rerun `Cross_Citations`.** Its sector-level results are included in `Cross_Citations/output/`.

| Data | Folder | Provider and access |
|------|--------|---------------------|
| PatentsView granted-patent tables (about 2 GB) | `Cross_Citations/input/PatentsView/` | PatentsView final release (12/31/2024), USPTO record on Zenodo[^14], free, CC BY 4.0 |
| OECD Citations database, February 2025 (about 3 GB packed) | `Cross_Citations/input/OECD_Citations_Feb2025/` | OECD[^15], free on request from the OECD Directorate for Science, Technology and Innovation |
| OECD Triadic Patent Families and Patent Quality Indicators, February 2025 | `Cross_Citations/input/OECD_TPF_Feb2025/`, `OECD_Quality_Feb2025/` | OECD[^15], free on request from the OECD Directorate for Science, Technology and Innovation (robustness only) |
| ALP IPC–ISIC and USPC–ISIC concordances | `Cross_Citations/input/concordances/` (`ipc4_to_isic_rev31_4.txt`; `uspc_Class_to_isic_rev31_4.txt`, robustness only) | Lybbert & Zolas (2014)[^13]; Goldschlag, Lybbert & Zolas (2016)[^9]. Free download from the authors' page (step 6 below) |
| KPSS patent values, 2024 extension | `Cross_Citations/input/KPSS/` | Kogan, Papanikolaou, Seru & Stoffman (2017)[^16], public GitHub repository (robustness only) |

The README of each folder lists every file, its source and the step that reads it.

### Downloading the data that are not included

All of these files are free. Save each one under the exact name and folder shown; the code looks for these names.

**For the paper's results (`Empirical`)**

1. **IPUMS USA extract** → `Empirical/input/usa_00026.dta`. Register at https://usa.ipums.org, choose *Select Data*, and build the extract described above (samples, variables, all persons, `.dat` format). When the extract is ready, download the `.dat.gz` file and the Stata command file (`.do`), unzip the `.dat` file into the same folder, run the command file in Stata, and save the result as `usa_00026.dta` in `Empirical/input/`.
2. **IPUMS OCCSOC crosswalk, 2000 onward** → `Empirical/input/occ_occsoc_crosswalk_2000_onward.xlsx`. On the IPUMS page *Codes for Occupation SOC (OCCSOC) in the 2000 Census and the ACS/PRCS Samples from 2000 Onward* (https://usa.ipums.org/usa/volii/occsoc18.shtml), click *Excel format* (file `occsoc_crosswalk_2000_onward.xlsx`, https://usa.ipums.org/usa/resources/volii/occsoc_crosswalk_2000_onward.xlsx) and rename it. The code reads the sheet *Codes without Descriptions*.
3. **U.S. KLEMS, March 2017 release** → `Empirical/input/usa_wk_mar_2017.xlsx`. On the World KLEMS analytical data page (https://www.worldklems.net/data/analytical/), download the U.S. file `usa_wk_mar_2017.xlsx`.
4. **Caunedo, Jaume & Keller (2023) equipment and software prices** → `Empirical/input/manual/CJK/Prices_by_asset_type_with_software_091720.dta`. Download their replication package from Zenodo (https://doi.org/10.5281/zenodo.7779807, version 1.3) and copy the file `empirics/rawdata/PriceKquality/Prices_by_asset_type_with_software_091720.dta`. The same file is in the authors' GitHub repository, https://github.com/julicaunedo/Occupational_Exposure_to_CETC.
5. **USPTO historical patent masterfile, 2014 release** → `Empirical/input/historical_masterfile.dta`. On the USPTO Historical Patent Data Files page (https://www.uspto.gov/ip-policy/economic-research/research-datasets/historical-patent-data-files), click the DTA link *2014 historical_masterfile.dta* (758 MB). The page returns a temporary download link to paste into the browser. The same files are in the USPTO Bulk Data Directory, https://data.uspto.gov/bulkdata/datasets/histexc.
6. **ALP patent crosswalks, version 2209** → `Empirical/input/ALP_USPC/NAICS_2007/` and `Cross_Citations/input/concordances/`. On Nikolas Zolas's *Patent Crosswalks* page (https://sites.google.com/site/nikolaszolas/PatentCrosswalk), download *USPC Crosswalks, Version 2209* and *IPC Crosswalks, Version 2209* (the full sets at the top of the page, not the versions excluding service industries) and unzip them. Copy:
   - `ALP_USPC/NAICS_2007/uspc_Class_to_naics07_2.txt` to `Empirical/input/ALP_USPC/NAICS_2007/` (used by `Empirical` step 10 and `Cross_Citations` step 04);
   - `ALP_IPC/ISIC_Rev31/ipc4_to_isic_rev31_4.txt` to `Cross_Citations/input/concordances/` (only to rerun `Cross_Citations`, steps 03, 05 and 06);
   - `ALP_USPC/ISIC_Rev31/uspc_Class_to_isic_rev31_4.txt` to `Cross_Citations/input/concordances/` (robustness only, r1).
   The authors post these files publicly but state no license, so they are not redistributed here.

**Only to rerun `Cross_Citations`** (its results are included, so this is optional; it also needs the IPC crosswalk from step 6)

7. **PatentsView tables** → `Cross_Citations/input/PatentsView/`. From the USPTO's final PatentsView release on Zenodo (https://doi.org/10.5281/zenodo.15783125), download `g_patent.tsv.zip` (223 MB), `g_ipc_at_issue.tsv.zip` (372 MB) and `g_cpc_current.tsv.zip` (485 MB). Optional: `g_uspc_at_issue.tsv.zip` (612 MB, step 02 and r1). Robustness only: `g_inventor_disambiguated.tsv.zip` (673 MB) and `g_location_disambiguated.tsv.zip` (3 MB, r2 and r3). Keep the files zipped.
8. **OECD Citations database, February 2025** → `Cross_Citations/input/OECD_Citations_Feb2025/`. The OECD Directorate for Science, Technology and Innovation distributes its patent databases free of charge on request; registered users receive download links for each new edition by email. Save `202501_USPTO_CITATIONS_1.7z` to `_5.7z` and `202501_USPTO_Equivalent.7z` (about 3 GB). Step 06 unpacks them into `Cross_Citations/intermediate/` (about 24 GB); the unpacked `.txt` files can also be saved directly.
9. **OECD Triadic Patent Families database and Patent Quality Indicators database, February 2025** (robustness only, r3 and r4) → `Cross_Citations/input/OECD_TPF_Feb2025/` (`202501_TPF_USPTO.txt`, `202501_TPF_Core.txt`) and `Cross_Citations/input/OECD_Quality_Feb2025/` (`202502_OECD_PATENT_QUALITY_USPTO_INDIC.txt`). Same request as step 8.
10. **KPSS patent values, extended to 2024** (robustness only, r5) → `Cross_Citations/input/KPSS/KPSS_2024.csv`. Download `KPSS_2024.zip` from https://github.com/KPSS2017/Technological-Innovation-Resource-Allocation-and-Growth-Extended-Data and unzip it.

## Computational Requirements

### Software

- **Python** 3.9 or later (run with Python 3.9.6), with the packages pinned in `Cross_Citations/requirements.txt`.
- **Stata** 18 or later (run with StataNow 19). No user-written packages.
- **Julia** 1.12, with every package version pinned in `Calibration/Manifest.toml`.

### Hardware and run time

The results in the paper were produced on an Apple MacBook Pro (Apple M5 Pro, 15 cores, 24 GB of memory) running macOS 26.6.

| Part | Command | Run time | Notes |
|------|---------|----------|-------|
| Cross_Citations (optional) | `python run_all.py` | about 25 minutes | Needs the PatentsView and OECD downloads and reads about 30 GB of text. Its results are included, so this step can be skipped |
| Empirical | `do master.do` | about 12 minutes | Step 11 loads the 4.5 GB IPUMS file; at least 16 GB of memory is recommended |
| Calibration | `julia --project=. run_all.jl --skip-robustness` | about 5 minutes | Figures, calibration, Table 2 and the exercises |
| Calibration robustness | `julia --project=. code/robustness.jl --redo` | about 3.5 hours | Reruns the ten recalibrations for Tables J.1 and H.2. Without `--redo`, `run_all.jl` reuses the stored runs in `output/robustness/` and only rebuilds the tables |

About 40 GB of free disk space is needed for a full rebuild of the patent tables; otherwise 10 GB is enough. No step uses random numbers, so every run gives the same results.

## Instructions to Replicators

Run the three parts in this order. Each README has the details.

1. **Cross_Citations** (optional; its results are included). After saving the patent downloads listed above:
   ```
   cd Cross_Citations
   python -m pip install -r requirements.txt
   python run_all.py
   ```
2. **Empirical.** Save the files of download steps 1–6 where shown above, then in Stata:
   ```stata
   cd "path/to/Empirical"
   do master.do
   ```
3. **Calibration:**
   ```
   cd Calibration
   julia --project=. -e 'using Pkg; Pkg.instantiate()'
   julia --project=. run_all.jl
   ```
   Use `--skip-calibration` to keep the stored calibration and `--skip-robustness` to leave out Tables J.1 and H.2. Step 5 reuses the alternative calibrations stored in `output/robustness/`; to rerun them, run `julia --project=. code/robustness.jl --redo` (about 3.5 hours).

## List of Tables and Figures

| Paper item | Folder | File | Script |
|------------|--------|------|--------|
| Figures 1–5, Figure B.1 | Calibration | `output/figures/Figure_1.pdf` … `Figure_5.pdf`, `Figure_B1.pdf` | `code/motivation_plots.jl` |
| Figures 6–7 | Calibration | `output/figures/Figure_6.pdf`, `Figure_7.pdf` | `code/calibration_plots.jl` |
| Table 1 | — | Externally set parameters, constants in `Calibration/code/smm_ge_unified.jl`; initial knowledge stocks from `Cross_Citations/output/patents/ex_cpc_hist_stocks.csv` | — |
| Table 2 | Calibration | `output/calibration/Table_2.tex` | `code/calibration_table.jl` |
| Table 3, Table H.1 | Calibration | `output/counterfactuals/Table_3.tex`, `Table_H1.tex` | `code/counterfactuals/Exercise1.jl` |
| Table 4 | Calibration | `output/counterfactuals/Table_4.tex` | `code/counterfactuals/Exercise2.jl` |
| Tables 5, 6, H.3, H.4 | Calibration | `output/counterfactuals/Table_5.tex`, `Table_6.tex`, `Table_H3.tex`, `Table_H4.tex` | `code/counterfactuals/Exercise3.jl` |
| Table 7 | Calibration | `output/counterfactuals/Table_7.tex` | `code/counterfactuals/Exercise4.jl` |
| Table 8 | Calibration | `output/counterfactuals/Table_8.tex` | `code/counterfactuals/Exercise5.jl` |
| Tables J.1, H.2 | Calibration | `output/robustness/Table_J1.tex`, `Table_H2.tex` | `code/robustness.jl` |
| Table A.1 | Empirical | `input/manual/concordance_tables.xlsx` | — (input) |
| Table D.1 | Cross_Citations | `input/concordances/isic31_4digit_to_sector.csv` | — (input) |
| Table D.2 | Cross_Citations | `output/citations/Table_D2.tex` | `code/06_citation_probabilities.py` |
| Table E.1 | Empirical | `output/dataset/RD_intensity.xlsx` | `code/11_scientists.do` |
| Data series behind the figures and targets | Empirical | `output/dataset/data_moments.xlsx` | `master.do` |

The README of each folder maps every number quoted in the text to the script that produces it.

## License and Citation

The code is released under the MIT License (see the `LICENSE` file at the top level and in each folder). Third-party data remain under their providers' terms, listed above.

If you use this code or data, please cite the paper and the replication package:

> Elbadri, S. & Vizcaino, J. I. (2026). *An Anatomy of the Baumol Cost Disease*. Working paper.

> Elbadri, S. & Vizcaino, J. I. (2026). *Replication package for "An Anatomy of the Baumol Cost Disease"*. Zenodo. https://doi.org/10.5281/zenodo.24265690. Code: https://github.com/selbadri/EV_Baumol2026_Replication

## Contact

Selim Elbadri: selim.elbadri@nottingham.ac.uk  
Juan I. Vizcaino: juan.vizcaino@nottingham.ac.uk

## References

[^1]: Elbadri, S. & Vizcaino, J. I. (2025). *Measuring US Industry-Level Productivity between 1947 and 2023*. Dataset, Zenodo. https://doi.org/10.5281/zenodo.17600449
[^2]: U.S. Bureau of Economic Analysis. *Fixed Assets Accounts Tables*, detailed estimates by industry and asset (accessed 9 May 2026). https://apps.bea.gov/national/FA2004/Details/Index.htm
[^3]: U.S. Bureau of Labor Statistics (2025). *Rental Prices*, Productivity program tables. https://www.bls.gov/productivity/tables/
[^4]: Jorgenson, D. W., Ho, M. S. & Samuels, J. D. (2017). Educational Attainment and the Revival of US Economic Growth. In *Education, Skills, and Technical Change*, University of Chicago Press, 23–60. Data: https://www.worldklems.net/data
[^5]: Feenstra, R. C., Inklaar, R. & Timmer, M. P. (2015). The Next Generation of the Penn World Table. *American Economic Review* 105(10), 3150–3182. https://www.rug.nl/ggdc/productivity/pwt/
[^6]: Jordà, Ò., Schularick, M. & Taylor, A. M. (2017). Macrofinancial History and the New Business Cycle Facts. *NBER Macroeconomics Annual* 31(1), 213–263; Jordà, Ò., Knoll, K., Kuvshinov, D., Schularick, M. & Taylor, A. M. (2019). The Rate of Return on Everything, 1870–2015. *Quarterly Journal of Economics* 134(3), 1225–1298. Release 6: https://www.macrohistory.net/database/
[^7]: Caunedo, J., Jaume, D. & Keller, E. (2023). Occupational Exposure to Capital-Embodied Technical Change. *American Economic Review* 113(6), 1642–1685. Replication package: https://doi.org/10.5281/zenodo.7779807
[^8]: Marco, A. C., Carley, M., Jackson, S. & Myers, A. (2015). The USPTO Historical Patent Data Files: Two Centuries of Innovation. https://www.uspto.gov/ip-policy/economic-research/research-datasets/historical-patent-data-files
[^9]: Goldschlag, N., Lybbert, T. J. & Zolas, N. J. (2016). An "Algorithmic Links with Probabilities" Crosswalk for USPC and CPC Patent Classifications. US Census Bureau CES Working Paper 16-15. https://sites.google.com/site/nikolaszolas/PatentCrosswalk
[^10]: Ruggles, S. et al. (2025). IPUMS USA: Version 16.0. Minneapolis, MN: IPUMS. https://usa.ipums.org/usa/
[^11]: National Center for Science and Engineering Statistics. *National Survey of College Graduates 2023*, public-use files and Table 3-1 (NSF 25-322). https://ncses.nsf.gov/surveys/national-survey-college-graduates/2023
[^12]: National Center for O*NET Development. *Browse STEM Occupations*. O*NET OnLine. https://www.onetonline.org/find/stem?t=0
[^13]: Lybbert, T. J. & Zolas, N. J. (2014). Getting Patents and Economic Data to Speak to Each Other: An "Algorithmic Links with Probabilities" Approach for Joint Analyses of Patenting and Economic Activity. *Research Policy* 43(3), 530–542. Concordances: https://sites.google.com/site/nikolaszolas/PatentCrosswalk
[^14]: United States Patent and Trademark Office (2025). *Final Release of PatentsView Metadata, Pre-Grant and Granted (12/31/2024)*. Zenodo. https://doi.org/10.5281/zenodo.15783125
[^15]: OECD Citations database, February 2025; OECD Triadic Patent Families database, February 2025; OECD Patent Quality Indicators database, February 2025. Distributed free of charge on request by the OECD Directorate for Science, Technology and Innovation; derived from EPO PATSTAT, Autumn 2024. Method: Webb, C., Dernis, H., Harhoff, D. & Hoisl, K. (2005), Analysing European and International Patent Citations, OECD STI Working Paper 2005/9; Squicciarini, M., Dernis, H. & Criscuolo, C. (2013), Measuring Patent Quality, OECD STI Working Paper 2013/03, https://doi.org/10.1787/5k4522wkw1r8-en.
[^16]: Kogan, L., Papanikolaou, D., Seru, A. & Stoffman, N. (2017). Technological Innovation, Resource Allocation, and Growth. *Quarterly Journal of Economics* 132(2), 665–712. Data: https://github.com/KPSS2017/Technological-Innovation-Resource-Allocation-and-Growth-Extended-Data
