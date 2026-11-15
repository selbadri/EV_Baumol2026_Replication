*=======================================================================================
* An Anatomy of the Baumol Cost Disease -- Elbadri & Vizcaino (2026)
*---------------------------------------------------------------------------------------
* File:    06_labor_quality_growth.do
* Purpose: Compute Tornqvist-weighted labor quality index by broad sector
*---------------------------------------------------------------------------------------
* Inputs:  input/EV_production_accounts_1947to2023.xlsx  -- EV (2025) industry accounts
*          intermediate/concordances/EV_44Ind_concordance.dta  -- industry-to-sector concordance
*          input/usa_wk_mar_2017.xlsx  -- U.S. KLEMS, labor input of Other services to 1987
*
* Outputs: intermediate/06_labor_quality_growth.dta
*            LABQLTY_{GOODS,RI,NRI}  -- labor quality index (1947 = 1.00)
*=======================================================================================

set varabbrev off

*=======================================================================================
* STEP 1: Compute labor compensation (Tornqvist) weights within each sector
*=======================================================================================

* Labor quality is an index of labor input, so industries are weighted by their share of
* sector labor compensation. Value-added weights gave real estate, whose value added is
* mostly imputed rent on owner-occupied housing, about a quarter of NRI-S against 2-3
* percent of its labor, and its small workforce's labor-quality swings (+35 percent to
* 1968, -30 percent to 1976) dominated the sector index.
ev_load

* Other services, except government (61), through 1987: the EV accounts take labor input
* from the BEA-BLS experimental account for 1947-2016 (Eldridge et al. 2020), whose
* non-college labor quantity for this industry falls 33 percent over 1963-1980 while its
* hours rise 89 percent, so measured labor quality falls by two thirds. U.S. KLEMS (March
* 2017), with nearly identical nominal labor compensation, has real labor input growing
* 3.9 percent a year over 1963-1987. Its real labor input growth replaces the EV growth
* through 1987; EV hours are kept, and the two sources agree after 1987.
preserve
import excel "$raw/usa_wk_mar_2017.xlsx", sheet("KLEMdata") cellrange(A3) clear
rename (A B C D E F G H I J) ///
    (year industry go_nom cap_nom lab_nom int_nom go_con cap_con lab_con int_con)
destring year industry lab_con, replace force
keep if industry == 61 & year <= 1987
keep year lab_con
rename lab_con wk_lab61
tempfile wk61
save `wk61'
restore
merge m:1 year using `wk61', nogenerate keep(master match)

xtset industry_id year
bysort SectorCode year: egen Total_LAB = total(LAB)
gen LABshare = LAB / Total_LAB
sort industry_id year
gen tornqv_weight = (LABshare + L.LABshare) / 2

*=======================================================================================
* STEP 2: Compute sector-wide labor quality growth rate
*=======================================================================================

gen dln_LAB = ln(LAB_QI)-ln(L.LAB_QI)
replace dln_LAB = ln(wk_lab61)-ln(L.wk_lab61) if industry_id == 61 & year <= 1987
assert !missing(dln_LAB) if industry_id == 61 & inrange(year, 1948, 1987)
gen dln_LABQLTY = dln_LAB-(ln(HRS_QI)-ln(L.HRS_QI))
gen LABQLTYg_contr_ = dln_LABQLTY*tornqv_weight
collapse (sum) LABQLTYg_contr_, by(year SectorCode)
reshape wide LABQLTYg_contr_, i(year) j(SectorCode) string

*=======================================================================================
* STEP 3: Turn into labor quality index by sector 
*=======================================================================================

local sect "GOODS RI NRI"
foreach x of local sect{
	gen LABQLTYg_cumsum_`x' = sum(LABQLTYg_contr_`x')
	gen LABQLTY_`x' = exp(LABQLTYg_cumsum_`x')
}
rename year Year
keep Year LABQLTY_*

save "$interm/06_labor_quality_growth.dta", replace

*=======================================================================================
*=======================================================================================




