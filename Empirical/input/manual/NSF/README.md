# NSF business R&D benchmarks

NSF benchmarks for Appendix B (Table B.2 and Figure B.1). No step of this package reads them;
the figure is drawn by the calibration code from `nsf_rd_benchmarks.csv`. No series in
`data_moments.xlsx` depends on these files.

## Source files

| File | Source | Content |
|---|---|---|
| `h-19.xls` | NSF, Research and Development in Industry, historical Table H-19 | FTE R&D scientists and engineers in R&D-performing companies, January of each year, by SIC industry, 1954, 1957-58, 1963-64, 1967-68, 1971-99 (thousands) |
| `sird2007_tab67.xls` | NSF, Survey of Industrial Research and Development 2007, Table 67 | FTE R&D scientists and engineers, January 2008, by NAICS industry (thousands) |
| `berd2023_tab046.xlsx` | NSF, Business Enterprise Research and Development Survey 2023, Table 46 | Worldwide, domestic and foreign R&D employees by NAICS industry (thousands) |
| `USPRIV_monthly.csv` | FRED series USPRIV (BLS CES), downloaded 30 September 2026 | All employees, total private, monthly, 1939-01 to 2026-08 (thousands) |

## Processed file `nsf_rd_benchmarks.csv`

One row per NSF year. `rd_thous` is R&D employment in thousands; `usp_thous` is the annual
average of USPRIV (from `USPRIV_monthly.csv`); `rd_pct_priv = 100 * rd_thous / usp_thous`.

- `mfg_share` (H-19 years): 100 x (total - nonmanufacturing) / total. Nonmanufacturing is
  reported from 1963; it is suppressed in 1990 and 1991. Before the 1990s the survey sample
  covered few service firms, so the level of this share overstates manufacturing.
- Three-sector shares, our sector definitions:
  - 1999 (H-19, SIC). Total 997.7. Goods = manufacturing 659.1 + utilities + other
    nonmanufacturing (agriculture, mining, construction) 43.4. Transportation and utilities
    is 16.1, of which communications 14.1, so utilities and transportation are 2.0; we assign
    1.0 to Goods and 1.0 to NRI (the table does not split them). RI = communications 14.1 +
    finance, insurance and real estate 17.5 + computer and data processing 122.8 +
    engineering and management services 57.6. NRI = the remainder (trade, health and other
    services, transportation). Shares 70.5 / 21.2 / 8.2.
  - 2007 (SIRD, NAICS, January 2008 count). Total 1,130.5. Goods = manufacturing 717.1 +
    utilities 0.7 + the residual suppressed in the table 12.6. RI = information 139.1 +
    finance and insurance 19.3 + professional, scientific and technical services 196.1.
    NRI = wholesale 15.1 + retail 14.7 + health care 7.7 + other nonmanufacturing 8.1.
    Shares 64.6 / 31.4 / 4.0.
  - 2023 (BERD, domestic R&D employees). Total 2,126. Goods = manufacturing 1,047 + mining 4
    + utilities 6. RI = information 497 + finance and insurance 90 + professional services
    348. NRI = wholesale 22 + electronic shopping 31 + transportation 27 + real estate 5 +
    health care 10 + other nonmanufacturing (construction, retail, management, admin) 40.
    Shares 49.7 / 44.0 / 6.3.
- Concepts: 1954-2007 count FTE R&D scientists and engineers; 2023 counts all domestic R&D
  employees (including technicians and support staff), so the figure plots it as a separate
  series. All NSF surveys cover R&D-performing companies only and exclude universities,
  hospitals and nonprofit institutions.
