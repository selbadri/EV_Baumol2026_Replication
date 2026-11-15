R&D activity rates by occupation, used by code/11_scientists.do (step 4).

nsf25322-tab003-001.xlsx
  NCSES, National Survey of College Graduates 2023, Table 3-1: employed college graduates, by
  level of highest degree, minor occupation, and primary work activity (NSF 25-322).
  https://ncses.nsf.gov/pubs/nsf25322/assets/data-tables/tables/nsf25322-tab003-001.xlsx
  Downloaded 26 September 2026.

nscg_rd_rates_occ1990.csv
  Maps each occ1990 code with R&D occupations to an NSCG minor occupation. nscg_group is the
  NSCG minor occupation group code (N3OCPRNG; ENG pools all engineering groups) used by
  step 4b of code/11_scientists.do, which computes from the NSCG 2023 public-use microdata the
  share of college graduates whose primary work activity is R&D, by group and sector. The
  columns nscg_rd_primary and nscg_total give, for all degree levels, the number (thousands)
  whose primary work activity is research and development and the total, from Table 3-1; they
  are kept as a check (the microdata reproduce them within 0.009) and were the weights until
  30 September 2026.

input/NSCG_2023/college_grads_2023.zip
  NCSES, National Survey of College Graduates 2023, public-use data files and documentation
  (codebook Ppcg23, data dictionary Dpcg23, annotated questionnaire).
  https://ncses.nsf.gov/822/assets/0/files/college_grads_2023.zip
  Downloaded 30 September 2026 (64 MB). Step 4b reads N3OCPRNG, WAPRI, WTSURVY and INDCODE. Architects, foresters and conservation scientists, and urban planners carry zero:
  their work is design, land management and planning. Industrial engineers are suppressed in the
  NSCG and take the rate for all engineers.
