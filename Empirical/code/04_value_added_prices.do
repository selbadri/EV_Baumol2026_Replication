*=======================================================================================
* An Anatomy of the Baumol Cost Disease -- Elbadri & Vizcaino (2026)
*---------------------------------------------------------------------------------------
* File:    04_value_added_prices.do
* Purpose: Compute Tornqvist sectoral value added price indices and aggregate price growth
*---------------------------------------------------------------------------------------
* Inputs:  input/EV_production_accounts_1947to2023.xlsx  -- EV (2025) industry accounts
*          intermediate/concordances/EV_44Ind_concordance.dta  -- industry-to-sector concordance
*
* Outputs: intermediate/04_value_added_prices.dta
*            VAPRICE_{GOODS,RI,NRI}  -- relative sectoral VA price (1947 = 1.00)
*            VAP_AGGREGATE           -- private-economy VA price index, Tornqvist (1947 = 1.00; used by 07c, 08)
*            dln_VAP                 -- aggregate VA price growth rate
*=======================================================================================

set varabbrev off

*=======================================================================================
* STEP 1: Compute value added (Tornqvist) weights within each sector
*=======================================================================================

ev_load
ev_tornqvist_weights

*=======================================================================================
* STEP 2: Compute change in value added prices by industry-year
*=======================================================================================

egen VA_base = max(cond(year==1947, VA, .)), by(industry_id)
gen VA_I = VA/VA_base
gen VAP = VA_I/VA_QI
gen dln_VAP = (ln(VAP)-ln(L.VAP))

*=======================================================================================
* STEP 3: Compute value added and value added price growth by sector-year
*=======================================================================================

gen VAPg_contr_ = dln_VAP*tornqv_weight
rename VA VA_
collapse (sum) VAPg_contr_ VA_, by(SectorCode year)
reshape wide VAPg_contr_ VA_, i(year) j(SectorCode) string

*=======================================================================================
* STEP 4: Compute aggregate value added price change
*=======================================================================================

tsset year
egen TOTAL_VA = rowtotal(VA_GOODS VA_NRI VA_RI)
local sect "GOODS RI NRI"
foreach x of local sect{
	gen VASHARE_`x' = VA_`x'/TOTAL_VA
	gen lag_VASHARE_`x' = L.VASHARE_`x'
	gen tornqv_weight_`x' = (VASHARE_`x'+lag_VASHARE_`x')/2
	gen contr_`x' = tornqv_weight_`x'*VAPg_contr_`x'
	drop VASHARE_`x' lag_VASHARE_`x' tornqv_weight_`x'
}
egen VAPg_contr_AGGREGATE = rowtotal(contr_GOODS contr_RI contr_NRI)

*=======================================================================================
* STEP 5: Turn sector and aggregate value added price growth into indices
*=======================================================================================

local sect_aggregate "GOODS RI NRI AGGREGATE"
foreach x of local sect_aggregate {
	gen VAPg_cumsum_`x' = sum(VAPg_contr_`x')
	gen VAP_`x' = exp(VAPg_cumsum_`x')
}
rename year Year
keep Year VAP_* VAPg_contr_AGGREGATE

*=======================================================================================
* STEP 6: Divide sectoral value added prices by aggregate value added price 
*=======================================================================================

local sect "GOODS RI NRI"
foreach x of local sect{
	gen relVAP_`x' = VAP_`x'/VAP_AGGREGATE
	drop VAP_`x'
	rename relVAP_`x' VAPRICE_`x'
}
rename VAPg_contr_AGGREGATE dln_VAP

save "$interm/04_value_added_prices.dta", replace

*=======================================================================================
*=======================================================================================




