# An Anatomy of the Baumol Cost Disease: Empirical Data

**Selim Elbadri & Juan I. Vizcaino (2026)**

## Description

This folder builds the data used in *An Anatomy of the Baumol Cost Disease*. A Stata pipeline constructs annual series for 1947–2023 for the US private economy, split into three sectors:

| Code | Sector |
|------|--------|
| `GOODS` | Goods-producing industries |
| `RI` | Research-intensive services |
| `NRI` | Non-research-intensive services |

The series cover labor shares, value added shares, relative prices, labor productivity, labor quality, machine-producer markups, capital prices, population, patent shares and stocks, and R&D worker shares. They are compiled into `output/dataset/data_moments.xlsx`, which the calibration code reads to draw the paper's figures and set the model targets. The industry-to-sector mapping is in `input/manual/concordance_tables.xlsx` (Appendix A, Table A.1).

The main sources are the industry production accounts of Elbadri & Vizcaino (2025)[^1], BEA fixed assets[^2], BLS capital detail[^3], U.S. KLEMS[^4], the Penn World Table[^5], the Jordà–Schularick–Taylor Macrohistory Database[^6], the quality-adjusted equipment prices of Caunedo, Jaume & Keller (2023)[^7], USPTO patent data[^8] with the Goldschlag, Lybbert & Zolas (2016) concordance[^9], IPUMS USA[^10], the NSCG 2023[^11] and O*NET[^12].

## Attribution & License

If you use this data or code, please cite the paper:

> Elbadri, S. & Vizcaino, J. I. (2026). *An Anatomy of the Baumol Cost Disease*. Working paper.

The industry production accounts in `input/EV_production_accounts_1947to2023.xlsx` come from a separate dataset, which should be cited on its own[^1].

The code in this folder is released under the MIT License (see `LICENSE`). Third-party data in `input/` remain under their providers' terms; in particular, the JST database is CC BY-NC-SA 4.0 (see `input/manual/JST/README`), and the IPUMS files, the U.S. KLEMS file and the Caunedo–Jaume–Keller price file are not redistributed here.

## Repository Structure

```
Empirical/
├── master.do                  Runs the full pipeline
├── README.md
├── LICENSE
├── code/                      Steps 00–12 (see table below)
├── input/                     Raw data, never modified by the pipeline
│   ├── EV_production_accounts_1947to2023.xlsx   Industry production accounts [1]
│   ├── detailnonres_stk1.xlsx, _stk2.xlsx, _inv2.xlsx   BEA fixed assets by industry and asset [2]
│   ├── rental-prices-major-industries.xlsx      BLS rental prices, stocks, deflators [3]
│   ├── usa_wk_mar_2017.xlsx                     U.S. KLEMS, March 2017 release [4], not included
│   ├── pwt110.dta                               Penn World Table 11.0 [5]
│   ├── historical_masterfile.dta                USPTO historical patent masterfile [8], not included (size)
│   ├── ALP_USPC/NAICS_2007/                     USPC–NAICS 2007 crosswalk [9], not included
│   ├── usa_00026.dta                            IPUMS USA extract [10], not included
│   ├── occ_occsoc_crosswalk_2000_onward.xlsx    IPUMS OCCSOC crosswalk [10], not included
│   ├── NSCG_2023/college_grads_2023.zip         NSCG 2023 public-use files [11]
│   ├── All_STEM_Occupations.xlsx                O*NET STEM occupations [12]
│   └── manual/                                  Files curated by the authors, each folder with a README
│       ├── concordance_tables.xlsx              Industry-to-sector concordances
│       ├── occupations_manual_assignment.xlsx   R&D occupations (Appendix B, Table B.3)
│       ├── scientists_industry_crosswalk.xlsx   INDNAICS to industry groups
│       ├── CJK/                                 Quality-adjusted equipment prices [7] (price file not included)
│       ├── JST/                                 Macrohistory Database, R6 [6]
│       ├── NSCG/                                NSCG Table 3-1 and occupation mapping [11]
│       └── NSF/                                 NSF R&D benchmarks (Appendix B; read by ../Calibration, not by Stata)
├── intermediate/              Created by the pipeline (not tracked)
└── output/
    ├── master_log.log
    └── dataset/
        ├── data_moments.xlsx  Compiled series
        └── RD_intensity.xlsx  R&D workers per 100 employees by industry (Appendix E, Table E.1)
```

### Steps

| Step | File | Builds | Series in `data_moments.xlsx` |
|------|------|--------|------|
| 00 | `00_programs.do` | Shared programs (`ev_load`, `ev_tornqvist_weights`) | |
| 01 | `01_concordance_tables.do` | Industry-to-sector concordances | |
| 02 | `02_labor_shares.do` | Labor share of value added | `LABVASHARE_*` |
| 03 | `03_value_added_shares.do` | Share of nominal value added | `VASHARE_*` |
| 04 | `04_value_added_prices.do` | Value added prices relative to the private economy | `VAPRICE_*`, `VAP_AGGREGATE` |
| 05 | `05_labor_productivity_growth.do` | Labor productivity | `LP_*`, `LP_AGGREGATE` |
| 06 | `06_labor_quality_growth.do` | Labor quality | `LABQLTY_*` |
| 07a | `07a_capital_assets.do` | BEA capital asset panel | |
| 07b | `07b_rates_of_return.do` | Internal and external rates of return | |
| 07c | `07c_markups.do` | Machine-producer profit shares (Appendix F) | `profit_share_*` |
| 08 | `08_capital_prices.do` | Capital prices, quality-adjusted equipment (Appendix G) | `CAPPRICE_*` |
| 09 | `09_demographics.do` | US population | `pop`, `pop_g` |
| 10 | `10_patents.do` | Patent shares, IPC route and USPC route (Appendix C) | `PATSH_*`, `PAT_USPC_*`, `PATSH_USPC_*` |
| 10b | `10b_patent_stocks.do` | Patent stocks | `PATSTOCK_*`, `PATSTOCKSH_*` |
| 11 | `11_scientists.do` | R&D worker shares and the RI/NRI classification (Appendices B, E) | `RDSH_*`, `SCI_SH` |
| 12 | `12_compilation.do` | Merges all series and writes `data_moments.xlsx` | |

Indices equal 1 in 1947 (the paper renormalizes to 1950). Shares are fractions, not percentages. Sheet `README` in `data_moments.xlsx` describes each variable.

### Not produced here

The paper's figures and model results come from the calibration code, which reads `data_moments.xlsx`. Steps 10 and 10b read the patent shares and stocks of the IPC route from `../Cross_Citations/output/patents/` (`ex_cpc_hist_flow_shares.csv`, `ex_cpc_hist_stocks.csv`), which are included in that folder. The citation probabilities also come from there.

## Instructions

**Requirements:** Stata 18 or later; no user-written packages. A full run takes about 12 minutes (StataNow 19, 2026). Step 11 loads a 4.5 GB microdata file, so allow enough memory.

1. Obtain the six files that are not included and save them where shown:
   - `usa_00026.dta`: an IPUMS USA extract[^10]. Samples: 1940 1%, 1950 1%, 1960 1%, 1970 1% state fm1, 1980 5% state, 1990 5% state, 2000 5%, 2010 ACS, 2020 ACS, 2023 ACS. Variables: `SEX`, `AGE`, `EDUC`, `EMPSTAT`, `OCC1990`, `OCCSOC`, `IND1990`, `INDNAICS` (with `EDUCD`, `EMPSTATD` and the variables IPUMS preselects). All persons, no case selection; `.dat` format, rectangular (person). Run the Stata command file IPUMS provides and save the result as `usa_00026.dta`.
   - `occ_occsoc_crosswalk_2000_onward.xlsx`: the IPUMS USA OCCSOC crosswalk for 2000 onward[^10]: on https://usa.ipums.org/usa/volii/occsoc18.shtml click *Excel format* (`occsoc_crosswalk_2000_onward.xlsx`) and save it under this name in `input/`.
   - `usa_wk_mar_2017.xlsx`: U.S. KLEMS, March 2017 release[^4], from https://www.worldklems.net/data/analytical/, saved in `input/`.
   - `Prices_by_asset_type_with_software_091720.dta`: from the replication package of Caunedo, Jaume & Keller (2023)[^7] (Zenodo record 7779807, folder `empirics/rawdata/PriceKquality/`), saved in `input/manual/CJK/`.
   - `ALP_USPC/NAICS_2007/uspc_Class_to_naics07_2.txt`: the USPC–NAICS 2007 crosswalk of Goldschlag, Lybbert & Zolas (2016)[^9]. On https://sites.google.com/site/nikolaszolas/PatentCrosswalk download *USPC Crosswalks, Version 2209* (the full set, not the version excluding service industries), unzip it and copy `ALP_USPC/NAICS_2007/uspc_Class_to_naics07_2.txt` to `input/ALP_USPC/NAICS_2007/`. The authors post it publicly but state no license, so it is not redistributed here.
   - `historical_masterfile.dta`: the USPTO historical patent masterfile, 2014 release[^8], public domain but too large for GitHub. On https://www.uspto.gov/ip-policy/economic-research/research-datasets/historical-patent-data-files click the DTA link *2014 historical_masterfile.dta* (758 MB) and save it in `input/`.
2. Open Stata, change to this folder and run the master file:
   ```stata
   cd "path/to/Empirical"
   do master.do
   ```
   No edits are needed. `master.do` takes the root from the working directory and creates `intermediate/` and `output/`. To run it from elsewhere, set `global root "path/to/Empirical"` first.

The log is written to `output/master_log.log`.

## Contact

Selim Elbadri: selim.elbadri@nottingham.ac.uk  
Juan I. Vizcaino: juan.vizcaino@nottingham.ac.uk

## References

[^1]: Elbadri, S. & Vizcaino, J. I. (2025). *Measuring US Industry-Level Productivity between 1947 and 2023*. Dataset, Zenodo. https://doi.org/10.5281/zenodo.17600449. Code: https://github.com/selbadri/Measuring-US-Industry-Level-Productivity-1947-to-2023
[^2]: U.S. Bureau of Economic Analysis. *Fixed Assets Accounts Tables*, detailed estimates by industry and asset (accessed 9 May 2026). https://apps.bea.gov/national/FA2004/Details/Index.htm
[^3]: U.S. Bureau of Labor Statistics (2025). *Rental Prices*, Productivity program tables. https://www.bls.gov/productivity/tables/
[^4]: Jorgenson, D. W., Ho, M. S. & Samuels, J. D. (2017). Educational Attainment and the Revival of US Economic Growth. In *Education, Skills, and Technical Change*, University of Chicago Press, 23–60. Data: https://www.worldklems.net/data
[^5]: Feenstra, R. C., Inklaar, R. & Timmer, M. P. (2015). The Next Generation of the Penn World Table. *American Economic Review* 105(10), 3150–3182. https://www.rug.nl/ggdc/productivity/pwt/
[^6]: Jordà, Ò., Schularick, M. & Taylor, A. M. (2017). Macrofinancial History and the New Business Cycle Facts. *NBER Macroeconomics Annual* 31(1), 213–263; Jordà, Ò., Knoll, K., Kuvshinov, D., Schularick, M. & Taylor, A. M. (2019). The Rate of Return on Everything, 1870–2015. *Quarterly Journal of Economics* 134(3), 1225–1298. Release 6: https://www.macrohistory.net/database/
[^7]: Caunedo, J., Jaume, D. & Keller, E. (2023). Occupational Exposure to Capital-Embodied Technical Change. *American Economic Review* 113(6), 1642–1685.
[^8]: Marco, A. C., Carley, M., Jackson, S. & Myers, A. (2015). The USPTO Historical Patent Data Files: Two Centuries of Innovation. https://www.uspto.gov/ip-policy/economic-research/research-datasets/historical-patent-data-files
[^9]: Goldschlag, N., Lybbert, T. J. & Zolas, N. J. (2016). An "Algorithmic Links with Probabilities" Crosswalk for USPC and CPC Patent Classifications. US Census Bureau CES Working Paper 16-15. https://sites.google.com/site/nikolaszolas/PatentCrosswalk
[^10]: Ruggles, S. et al. (2025). IPUMS USA: Version 16.0. Minneapolis, MN: IPUMS. https://usa.ipums.org/usa/
[^11]: National Center for Science and Engineering Statistics. *National Survey of College Graduates 2023*, public-use files. https://ncses.nsf.gov/surveys/national-survey-college-graduates/2023
[^12]: National Center for O*NET Development. *Browse STEM Occupations*. O*NET OnLine. https://www.onetonline.org/find/stem?t=0
