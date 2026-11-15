Quality-adjusted equipment and software prices, used by code/08_capital_prices.do.

Prices_by_asset_type_with_software_091720.dta
  Source: Caunedo, J., D. Jaume and E. Keller, "Occupational Exposure to Capital-Embodied
  Technical Change", replication package, Zenodo record 7779807 (v1.3, identical in v3),
  file empirics/rawdata/PriceKquality/Prices_by_asset_type_with_software_091720.dta.
  Downloaded 25 September 2026 from https://github.com/julicaunedo/Occupational_Exposure_to_CETC
  Contents: price levels (2009 = 100) for 22 NIPA equipment types (nipa_code) and software
  (nipa_code 99, from 1960), 1947-2016. The authors note that the prices build on Cummins and
  Violante (2002).

cjk_bea_assets.csv   BEA fixed-asset codes (BEA_Panel assetcode) -> CJK nipa_code. Taken from
                     the "By detailed equipment" table of the Prices sheet of
                     Capital_embodiment_measure.xlsx in the same package; EO71 (household
                     appliances, absent there) is assigned to 41, and the three software codes
                     to 99.
cjk_bls_assets.csv   BLS capital-account asset categories -> CJK nipa_code, by asset name.

Assets with no CJK counterpart (structures, R&D, entertainment originals, software before
1960) keep their official deflators.
