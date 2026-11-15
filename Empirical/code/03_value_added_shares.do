*=======================================================================================
* An Anatomy of the Baumol Cost Disease -- Elbadri & Vizcaino (2026)
*---------------------------------------------------------------------------------------
* File:    03_value_added_shares.do
* Purpose: Compute value added shares by broad sector
*---------------------------------------------------------------------------------------
* Inputs:  input/EV_production_accounts_1947to2023.xlsx  -- EV (2025) industry accounts
*          intermediate/concordances/EV_44Ind_concordance.dta  -- industry-to-sector concordance
*
* Outputs: intermediate/03_value_added_shares.dta
*            VASHARE_{GOODS,RI,NRI}  -- share of total value added (annual, 1947-2023)
*=======================================================================================

set varabbrev off

*=======================================================================================
* STEP 1: Compute value added shares by broad sector
*=======================================================================================

ev_load
collapse (sum) VA, by(SectorCode year)
bys year: egen total_VA = total(VA)
gen VASHARE_ = VA/total_VA
drop VA total_VA
reshape wide VASHARE_, i(year) j(SectorCode) string
rename year Year
save "$interm/03_value_added_shares.dta", replace

*=======================================================================================
*=======================================================================================

