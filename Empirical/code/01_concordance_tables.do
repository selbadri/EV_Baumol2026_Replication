
*=======================================================================================
* An Anatomy of the Baumol Cost Disease -- Elbadri & Vizcaino (2026)
*---------------------------------------------------------------------------------------
* File:    01_concordance_tables.do
* Purpose: Convert industry-to-sector concordance tables from Excel to .dta format
*---------------------------------------------------------------------------------------
* Inputs:  input/manual/concordance_tables.xlsx      -- manually curated concordances
*
* Outputs: intermediate/concordances/EV_44Ind_concordance.dta        -- EV (2025) 44-industry        -> broad sectors
*          intermediate/concordances/EV_63Ind_concordance.dta        -- EV (2025) 63-industry        -> broad sectors
*          intermediate/concordances/USPTO_2digits.dta               -- NAICS 2-digit                -> broad sectors
*          intermediate/concordances/BEA_sector_map.dta              -- BEA sector   -> industry_id
*          intermediate/concordances/BLS_naics_map.dta               -- BLS naics    -> industry_id
*          intermediate/concordances/BroadSector_map.dta             -- industry_id  -> BroadSector
*=======================================================================================

set varabbrev off

*=======================================================================================
* STEP 1: Turn each concordance table into .dta
*=======================================================================================

*Import bridge and turn it into .dta
import excel "$raw_manual/concordance_tables.xlsx", sheet("EV_44Ind") firstrow clear
assert _N == 44
save "$conc/EV_44Ind_concordance.dta", replace

*Import bridge and turn it into .dta
import excel "$raw_manual/concordance_tables.xlsx", sheet("EV_63Ind") firstrow clear
assert _N == 63
save "$conc/EV_63Ind_concordance.dta", replace

*Import bridge and turn it into .dta
import excel "$raw_manual/concordance_tables.xlsx", sheet("USPTO_2digits") firstrow clear
*Rows with empty Sector are subtotal/header rows, not industries
assert _N == 24
save "$conc/USPTO_2digits.dta", replace

*Import BEA sector → industry_id crosswalk (used by 07b, 07c, 08)
import excel "$raw_manual/concordance_tables.xlsx", sheet("BEA_sector_map") firstrow clear
rename IndustryDescription indname
save "$conc/BEA_sector_map.dta", replace

*Import BLS naics → industry_id crosswalk (used by 07b, 07c, 08)
import excel "$raw_manual/concordance_tables.xlsx", sheet("BLS_naics_map") firstrow clear
rename (IndustryDescription code sector) (indname naics_code naics_sector)
save "$conc/BLS_naics_map.dta", replace

*Import industry_id → BroadSector crosswalk (used by 07c, 08)
import excel "$raw_manual/concordance_tables.xlsx", sheet("BroadSector_map") cellrange(A2) clear
rename A industry_id
rename C broadsector
keep industry_id broadsector
drop if missing(industry_id)
destring industry_id, replace force
save "$conc/BroadSector_map.dta", replace

*=======================================================================================
*=======================================================================================
