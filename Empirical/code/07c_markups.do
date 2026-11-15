*=======================================================================================
* An Anatomy of the Baumol Cost Disease -- Elbadri & Vizcaino (2026)
*---------------------------------------------------------------------------------------
* File:    07c_markups.do
* Purpose: Machine-producer profit shares by broad sector (Appendix F of the paper),
*          the markup input of the model
*---------------------------------------------------------------------------------------
* Construction
*   Machines are the produced nonresidential assets of the capital prices (step 08):
*   equipment, intellectual property products and structures. Land, inventories and the
*   tenant-occupied housing of real estate are excluded, so that markups and capital
*   prices cover the same assets. For machine assets m of industry i,
*       psi = 1 - sum_m RC_m / sum_m CC_m,
*       CC_m = r_i SV_m + B_m      (capital income at the industry's internal return)
*       RC_m = r*  SV_m + B_m      (required cost at the external return r*)
*       SV_m = T P(t-1) K,   B_m = T (delta P - EpiRel P(t-1)) K,
*   with EpiRel each asset's mean change in price relative to the private VA deflator. The
*   internal return r_i is re-solved over all assets of the industry so that its capital
*   income is unchanged: r_i = (CC_i - sum B) / sum SV. Revaluations then fall on every
*   asset per dollar of stock, and expected relative price changes enter both sides.
*   - Negative r_i is set to r* (BLS rule); r_i is capped at 30 percent.
*   - r* is a yield-based expected real return on wealth (JST R6, United States): dividend
*     yield + real GDP growth, rental yield + real house-price growth, real long rate and
*     real bill rate, weighted by regressing JST's realized return on wealth on the four
*     asset returns over 1950-2020 (weights rescaled to sum to one). Held at 2020 after 2020.
*   - BLS capital detail from 1988 (published capital costs, all assets, tax parameters);
*     before 1988, BEA fixed assets with the spliced internal rate of 07b, with capital
*     income and required cost each ratio-spliced to BLS at 1988 by industry.
*   - All private industries, real estate included.
*   The calibration uses a quadratic trend in time of each sector's gross markup 1/(1-psi)
*   over 1950-2023.
*
* Inputs:  intermediate/07_markups/BEA_Panel.dta           (07a)
*          intermediate/07_markups/IRR_ERR_1947to2024.dta  (07b)
*          intermediate/04_value_added_prices.dta          (04)
*          input/rental-prices-major-industries.xlsx       (BLS capital detail)
*          input/manual/JST/JSTdatasetR6.dta               (Jorda-Schularick-Taylor R6)
*          intermediate/concordances/{BEA_sector_map, BLS_naics_map, BroadSector_map}.dta
* Outputs: intermediate/07_markups.dta  (Year, profit_share_GOODS/RI/NRI, 1950-2023)
*          intermediate/07_markups/required_return.dta
*=======================================================================================

set varabbrev off

*---------------------------------------------------------------------------------------
* 1. Required real return on wealth (JST)
*---------------------------------------------------------------------------------------
use "$raw_manual/JST/JSTdatasetR6.dta", clear
keep if iso == "USA"
keep year eq_tr housing_tr bond_tr bill_rate capital_tr eq_dp housing_rent_yd ltrate cpi rgdpbarro pop hpnom
tsset year
gen double pi = ln(cpi) - ln(L.cpi)
* Expected inflation: mean of the five previous years (at least three)
gen double Epi = .
forvalues t = 2/`=_N' {
    quietly summarize pi in `=max(1, `t' - 5)'/`=`t' - 1'
    if r(N) >= 3 quietly replace Epi = r(mean) in `t'
}
regress capital_tr eq_tr housing_tr bond_tr bill_rate if inrange(year, 1950, 2020), noconstant
matrix b = e(b)
scalar wsum = b[1,1] + b[1,2] + b[1,3] + b[1,4]
scalar w_eq = b[1,1] / wsum
scalar w_h  = b[1,2] / wsum
scalar w_b  = b[1,3] / wsum
scalar w_s  = b[1,4] / wsum
gen double dlgdp = ln(rgdpbarro * pop) - ln(L.rgdpbarro * L.pop)
gen double dlhp  = (ln(hpnom) - ln(cpi)) - (ln(L.hpnom) - ln(L.cpi))
quietly summarize dlgdp if inrange(year, 1950, 2020)
scalar g_eq = r(mean)
quietly summarize dlhp if inrange(year, 1950, 2020)
scalar g_h = r(mean)
gen double rstar = w_eq * (eq_dp + g_eq) + w_h * (housing_rent_yd + g_h) ///
                 + w_b * (ltrate / 100 - Epi) + w_s * (bill_rate - Epi)
display as result "JST weights (eq, housing, bonds, bills): " %5.3f w_eq " " %5.3f w_h " " %5.3f w_b " " %5.3f w_s ///
    "   g_eq " %6.4f g_eq "   g_h " %6.4f g_h
keep year rstar
keep if inrange(year, 1945, 2023)
* Hold the last observed value (JST ends in 2020) through 2023
tsset year
tsfill
quietly summarize year
local add = 2023 - r(max)
if `add' > 0 {
    set obs `=_N + `add''
    replace year = year[_n-1] + 1 if missing(year)
}
sort year
replace rstar = rstar[_n-1] if missing(rstar)
label variable rstar "Required real return: yield-based expected return on wealth (JST R6)"
save "$interm/07_markups/required_return.dta", replace

*---------------------------------------------------------------------------------------
* Private VA price inflation, for expected relative price changes
*---------------------------------------------------------------------------------------
use "$interm/04_value_added_prices.dta", clear
rename Year year
sort year
gen double pi_va = ln(VAP_AGGREGATE) - ln(VAP_AGGREGATE[_n-1])
keep year pi_va
tempfile piva
save `piva'

*---------------------------------------------------------------------------------------
* 2. BEA fixed assets: capital income at the spliced internal rate, re-solved returns
*---------------------------------------------------------------------------------------
use "$interm/07_markups/BEA_Panel.dta", clear
keep if stk2 > 0 & P > 0 & !missing(stk2) & !missing(P)
drop LP
sort sector assetcode year
by sector assetcode: gen double LP = P[_n-1]
merge m:1 year using `piva', keep(master match) nogenerate
gen double dlrel = ln(P / LP) - pi_va
bysort sector assetcode: egen double EpiRel = mean(cond(inrange(year, 1948, 2023), dlrel, .))
drop if missing(LP)
merge m:1 sector using "$conc/BEA_sector_map.dta", keepusing(industry_id) keep(match) nogenerate
merge m:1 industry_id year using "$interm/07_markups/IRR_ERR_1947to2024.dta", keepusing(IRR_it) keep(match) nogenerate
* BEA nonresidential fixed assets: equipment, IPP and structures, all produced
gen byte machine = 1
gen double CC = (IRR_it * LP + delta * P - cap_gain_dollar) * stk2
gen double SV = LP * stk2
gen double B  = (delta * P - cond(missing(EpiRel), 0, EpiRel) * LP) * stk2
gen double SVm = SV * machine
gen double Bm  = B * machine
collapse (sum) CCi = CC SV B SVm Bm (max) anym = machine, by(industry_id year)
keep if anym == 1
gen double rt = (CCi - B) / SV
keep industry_id year rt SVm Bm
tempfile bea
save `bea'

*---------------------------------------------------------------------------------------
* 3. BLS capital detail, 1988-2023: published capital costs, re-solved returns
*---------------------------------------------------------------------------------------
import excel "$raw/rental-prices-major-industries.xlsx", sheet("Annual") cellrange(A3) firstrow clear
rename NAICS naics_code
rename Measure measure
rename AssetCategory asset
local year 1987
foreach var of varlist F-AQ {
    rename `var' v`year'
    local year = `year' + 1
}
foreach var of varlist v1987-v2024 {
    capture destring `var', replace force ignore(",")
}
keep if inlist(measure, "Capital asset price deflator", "Wealth stock depreciation rate", "Tax parameter", ///
    "Productive capital stock", "Capital costs based on internal rate of return")
gen str4 m = ""
replace m = "P"   if measure == "Capital asset price deflator"
replace m = "dep" if measure == "Wealth stock depreciation rate"
replace m = "T"   if measure == "Tax parameter"
replace m = "K"   if measure == "Productive capital stock"
replace m = "CCI" if measure == "Capital costs based on internal rate of return"
keep naics_code asset m v*
collapse (mean) v*, by(naics_code asset m)
reshape long v, i(naics_code asset m) j(year)
reshape wide v, i(naics_code asset year) j(m) string
rename v* *
preserve
    keep if asset == "No associated assets"
    keep naics_code year CCI
    merge m:1 naics_code using "$conc/BLS_naics_map.dta", keepusing(industry_id) keep(match) nogenerate
    collapse (sum) CCi = CCI, by(industry_id year)
    tempfile blscc
    save `blscc'
restore
drop if asset == "No associated assets"
keep if K > 0 & !missing(K)
replace T = 1 if missing(T)
replace dep = 0 if missing(dep)
sort naics_code asset year
by naics_code asset: gen double LP = P[_n-1]
merge m:1 year using `piva', keep(master match) nogenerate
gen double dlrel = ln(P / LP) - pi_va
bysort naics_code asset: egen double EpiRel = mean(cond(inrange(year, 1988, 2023), dlrel, .))
drop if missing(LP)
* Machines: produced nonresidential assets in the BLS asset list, as in step 08
gen byte machine = !inlist(asset, "Land stock", "Finished goods inventory", "Materials inventory", "Work-in-Process inventory") ///
    & strpos(lower(asset), "tenant-occupied") == 0
gen double SV = T * LP * K
gen double B  = T * (dep * P - cond(missing(EpiRel), 0, EpiRel) * LP) * K
gen double SVm = SV * machine
gen double Bm  = B * machine
merge m:1 naics_code using "$conc/BLS_naics_map.dta", keepusing(industry_id) keep(match) nogenerate
collapse (sum) SV B SVm Bm (max) anym = machine, by(industry_id year)
keep if anym == 1
merge 1:1 industry_id year using `blscc', keep(match) nogenerate
gen double rt = (CCi - B) / SV
keep industry_id year rt SVm Bm
tempfile bls
save `bls'

*---------------------------------------------------------------------------------------
* 4. Capital income and required cost on machines; BLS rule and cap
*---------------------------------------------------------------------------------------
foreach src in bea bls {
    use ``src'', clear
    merge m:1 industry_id using "$conc/BroadSector_map.dta", keepusing(broadsector) keep(match) nogenerate
    merge m:1 year using "$interm/07_markups/required_return.dta", keep(master match) nogenerate
    * Returns above 30 percent are capped; negative or missing returns take the required return
    gen double r_used = cond(rt >= 0 & !missing(rt), min(rt, 0.30), rstar)
    gen double CCm = r_used * SVm + Bm
    gen double RCm = rstar * SVm + Bm
    drop if missing(RCm)
    keep industry_id broadsector year CCm RCm
    tempfile `src'_m
    save ``src'_m'
}

* Splice BEA to BLS at 1988, industry by industry, separately for CC and RC
use `bls_m', clear
keep if year == 1988
rename (CCm RCm) (CC88_bls RC88_bls)
keep industry_id CC88_bls RC88_bls
tempfile b88
save `b88'
use `bea_m', clear
keep if year == 1988
rename (CCm RCm) (CC88_bea RC88_bea)
keep industry_id CC88_bea RC88_bea
merge 1:1 industry_id using `b88', keep(match) nogenerate
gen double k_cc = CC88_bls / CC88_bea
gen double k_rc = RC88_bls / RC88_bea
keep industry_id k_cc k_rc
tempfile kf
save `kf'
use `bea_m', clear
keep if year < 1988
merge m:1 industry_id using `kf', keep(match) nogenerate
replace CCm = CCm * k_cc
replace RCm = RCm * k_rc
drop if missing(CCm) | missing(RCm)
drop k_cc k_rc
append using `bls_m'
drop if year > 2023 | (year >= 1988 & missing(CCm))
drop if year < 1988 & missing(CCm)

collapse (sum) CCm RCm, by(broadsector year)
gen double psi = 1 - RCm / CCm
keep broadsector year psi
reshape wide psi, i(year) j(broadsector) string
rename psiGOODS profit_share_GOODS
rename psiRI    profit_share_RI
rename psiNRI   profit_share_NRI
rename year Year
keep if inrange(Year, 1950, 2023)
label variable profit_share_GOODS "Machine-producer profit share, Goods (Appendix G)"
label variable profit_share_RI    "Machine-producer profit share, RI-S (Appendix G)"
label variable profit_share_NRI   "Machine-producer profit share, NRI-S (Appendix G)"
sort Year
save "$interm/07_markups.dta", replace
foreach s in GOODS RI NRI {
    quietly gen double mk_`s' = 1 / (1 - profit_share_`s')
    quietly summarize mk_`s'
    display as result "`s' gross machine markup, 1950-2023 mean: " %6.4f r(mean)
}
