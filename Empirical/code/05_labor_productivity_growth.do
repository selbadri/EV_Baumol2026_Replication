*=======================================================================================
* An Anatomy of the Baumol Cost Disease -- Elbadri & Vizcaino (2026)
*---------------------------------------------------------------------------------------
* File:    05_labor_productivity_growth.do
* Purpose: Compute Tornqvist-weighted labor productivity index by broad sector and aggregate
*---------------------------------------------------------------------------------------
* Inputs:  input/EV_production_accounts_1947to2023.xlsx  -- EV (2025) industry accounts
*          intermediate/concordances/EV_44Ind_concordance.dta  -- industry-to-sector concordance
*
* Outputs: intermediate/05_labor_productivity_growth.dta
*            LP_{GOODS,RI,NRI,AGGREGATE}  -- labor productivity index (1947 = 1.00)
*=======================================================================================

set varabbrev off

*=======================================================================================
* STEP 1: Compute value added (Tornqvist) weights within each sector
*=======================================================================================

ev_load
ev_tornqvist_weights

*=======================================================================================
* STEP 2: Compute sector-wide labor productivity growth rate and nominal value added 
*=======================================================================================

gen dln_LP = ln(LPVA)-ln(L.LPVA)
gen LPg_contr_ = dln_LP*tornqv_weight
rename VA VA_
collapse (sum) LPg_contr_ VA_, by(year SectorCode)
reshape wide LPg_contr_ VA_, i(year) j(SectorCode) string

*=======================================================================================
* STEP 3: Compute aggregate labor productivity growth rate
*=======================================================================================

tsset year
egen TOTAL_VA = rowtotal(VA_GOODS VA_NRI VA_RI)
local sect "GOODS RI NRI"
foreach x of local sect{
	gen VASHARE_`x' = VA_`x'/TOTAL_VA
	gen lag_VASHARE_`x' = L.VASHARE_`x'
	gen tornqv_weight_`x' = (VASHARE_`x'+lag_VASHARE_`x')/2
	gen contr_`x' = tornqv_weight_`x'*LPg_contr_`x'
	drop VASHARE_`x' lag_VASHARE_`x' tornqv_weight_`x'
}
egen LPg_contr_AGGREGATE = rowtotal(contr_GOODS contr_RI contr_NRI)

*=======================================================================================
* STEP 4: Turn sector and aggregate labor productivity growth into indices
*=======================================================================================

local sect_aggregate "GOODS RI NRI AGGREGATE"
foreach x of local sect_aggregate{
	gen LPg_cumsum_`x' = sum(LPg_contr_`x')
	gen LP_`x' = exp(LPg_cumsum_`x')
}
rename year Year
keep Year LP_*

save "$interm/05_labor_productivity_growth.dta", replace

*=======================================================================================
*=======================================================================================




