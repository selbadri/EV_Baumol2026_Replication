*=======================================================================================
* An Anatomy of the Baumol Cost Disease -- Elbadri & Vizcaino (2026)
*---------------------------------------------------------------------------------------
* File:    02_labor_shares.do
* Purpose: Compute Tornqvist-weighted labor share of value added by broad sector
*---------------------------------------------------------------------------------------
* Inputs:  input/EV_production_accounts_1947to2023.xlsx  -- EV (2025) industry accounts
*          intermediate/concordances/EV_44Ind_concordance.dta  -- industry-to-sector concordance
*
* Outputs: intermediate/02_labor_shares.dta
*            LABVASHARE_{GOODS,RI,NRI}  -- labor share of value added (annual, 1947-2023)
*=======================================================================================

set varabbrev off

*=======================================================================================
* STEP 1: Compute labor compensation share of value added by broad sector
*=======================================================================================

ev_load
collapse (sum) LAB VA, by(SectorCode year)
gen LABVASHARE_ = LAB/VA
drop LAB VA
reshape wide LABVASHARE_, i(year) j(SectorCode) string
rename year Year
save "$interm/02_labor_shares.dta", replace

*=======================================================================================
*=======================================================================================



