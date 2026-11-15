*=======================================================================================
* An Anatomy of the Baumol Cost Disease -- Elbadri & Vizcaino (2026)
*---------------------------------------------------------------------------------------
* File:    08_capital_prices.do
* Purpose: Sectoral capital price indices relative to the private-economy VA price
*          (Appendix G of the paper), 1947-2023
*---------------------------------------------------------------------------------------
* Inputs:  input/rental-prices-major-industries.xlsx  -- BLS asset deflators, 1987-2024
*          intermediate/07_markups/BEA_Panel.dta      -- BEA asset prices, 1946-2024 (07a)
*          intermediate/07_markups/IRR_ERR_1947to2024.dta -- rates of return for the asset weights (07b)
*          input/manual/CJK/                          -- Caunedo-Jaume-Keller quality-adjusted
*                                                        equipment prices and asset concordances
*          input/usa_wk_mar_2017.xlsx                 -- capital compensation weights, 1947-2014
*          input/EV_production_accounts_1947to2023.xlsx -- capital compensation weights, 2015-2023
*          intermediate/04_value_added_prices.dta     -- private-economy VA price
*          intermediate/concordances/{BEA_sector_map, BLS_naics_map, BroadSector_map}.dta
* Outputs: intermediate/08_capital_prices.dta
*            CAPPRICE_{GOODS,RI,NRI} -- relative capital price index (1947 = 1.00)
*---------------------------------------------------------------------------------------
* Coverage: produced nonresidential assets (equipment, IPP, structures); land, inventories
*   and tenant-occupied housing are dropped.
* Prices:   CJK quality-adjusted for equipment and software, official for the rest (Part 0).
* Method:   Tornqvist aggregation, asset -> industry with weights R x K (R the rental price
*   at the internal rate, or at the external rate where that is not positive), and
*   industry -> sector with capital-compensation weights. BEA growth rates cover 1948-1987
*   and BLS growth rates 1988-2023.
*=======================================================================================

*====================================================================
* PART 0: Quality-adjusted equipment and software prices
*   Caunedo, Jaume and Keller (CJK), built on Cummins and Violante (2002):
*   22 equipment types and software, 1947-2016 (input/manual/CJK, see README).
*   Through 2016 an asset's official price change is replaced by the CJK change for
*   its type. From 2017, where CJK ends, the official change is kept and each type's
*   average quality gap (CJK minus official) over 2001-2016 is added, so the
*   adjustment continues at its recent rate without a break in 2016. Structures,
*   R&D, entertainment originals and software before 1960 keep official prices.
*====================================================================

use "$raw_manual/CJK/Prices_by_asset_type_with_software_091720.dta", clear
keep nipa_code price_*
reshape long price_, i(nipa_code) j(year)
sort nipa_code year
by nipa_code: gen double d_ln_cjk = ln(price_) - ln(price_[_n-1])
keep nipa_code year d_ln_cjk
tempfile cjk_long
save `cjk_long'

import delimited "$raw_manual/CJK/cjk_bea_assets.csv", clear varnames(1) stringcols(1)
tempfile cjk_bea
save `cjk_bea'

import delimited "$raw_manual/CJK/cjk_bls_assets.csv", clear varnames(1) stringcols(1) bindquote(strict)
tempfile cjk_bls
save `cjk_bls'

* Official change by CJK type from the BEA asset deflators (common across industries;
* the median guards against strays), then the 2001-2016 average gap.
use "$interm/07_markups/BEA_Panel.dta", clear
keep if P > 0 & !missing(P)
collapse (median) P, by(assetcode year)
sort assetcode year
by assetcode: gen double d_off = ln(P) - ln(P[_n-1])
merge m:1 assetcode using `cjk_bea', keep(match) nogenerate
collapse (mean) d_off, by(nipa_code year)
merge 1:1 nipa_code year using `cjk_long', keep(match) nogenerate
keep if inrange(year, 2001, 2016)
gen double g = d_ln_cjk - d_off
collapse (mean) gap = g, by(nipa_code)
tempfile cjk_gap
save `cjk_gap'

*====================================================================
* PART A: BLS period (1988–2024) — asset → industry Törnqvist
*====================================================================

import excel "$raw/rental-prices-major-industries.xlsx", ///
    sheet("Annual") cellrange(A3) firstrow clear

rename NAICS         naics
rename Sector        sector
rename Measure       measure
rename AssetCategory assetcategory
rename Units         units

local year 1987
foreach var of varlist F-AQ {
    rename `var' _`year'
    local year = `year' + 1
}
foreach var of varlist _1987-_2024 {
    capture destring `var', replace force ignore(",")
}

*--- A1: Keep required measures -------------------------------------
keep if inlist(measure, ///
    "Capital asset price deflator",              ///
    "Productive capital stock",                  ///
    "Rental price with internal rate of return", ///
    "Rental price with external rate of return")
drop if assetcategory == "No associated assets"
drop units

* Produced assets only: land and inventories are dropped. The model's capital is the output
* of machine producers; land is not produced and inventories are working capital. BLS
* includes both from 1987 while the BEA assets used before 1988 include neither, so the
* restriction also gives both periods the same coverage.
drop if assetcategory == "Land stock" | strpos(assetcategory, "inventory") > 0
* Tenant-occupied housing (real estate; with its acquisition and disposal costs) is dropped
* too: the BEA nonresidential assets used before 1988 exclude it, so keeping it would change
* the NRI-S coverage in 1988 (about 75 percent of the real-estate weight; decided 6 Oct 2026).
drop if strpos(lower(assetcategory), "tenant-occupied") > 0

*--- A2: Reshape wide → long ----------------------------------------
reshape long _, i(naics sector measure assetcategory) j(year)
rename _ value

*--- A3: Pivot measures to columns ----------------------------------
gen byte mnum = .
replace mnum = 1 if measure == "Capital asset price deflator"
replace mnum = 2 if measure == "Productive capital stock"
replace mnum = 3 if measure == "Rental price with internal rate of return"
replace mnum = 4 if measure == "Rental price with external rate of return"
drop measure
reshape wide value, i(naics sector assetcategory year) j(mnum)
rename value1 P_deflator
rename value2 k_prod
rename value3 R_irr
rename value4 R_err
label variable P_deflator "Capital asset price deflator (BLS)"
label variable k_prod     "Productive capital stock (BLS)"
label variable R_irr      "Rental price with IRR (BLS)"
label variable R_err      "Rental price with ERR (BLS)"

*--- A4: Effective weight (same rule as files 03/04) ----------------
gen double w_eff = cond(R_irr > 0, R_irr * k_prod, R_err * k_prod)
* BLS R_ERR can be negative. A negative capital-compensation weight has no meaning as a
* Tornqvist share (securities, 2008, gave shares of -4.9 and +11.5), so floor at zero; the
* shares below renormalise over the remaining assets.
replace w_eff = 0 if w_eff < 0
label variable w_eff "Effective weight (R_IRR×k if R_IRR>0, else R_ERR×k), floored at 0"

*--- A5: Asset-level log price change -------------------------------
* Sort so _n-1 is the prior year for the same asset×sector
* Year 1987 is retained here: its d_ln_P will be missing (no 1986
* in BLS), but its w_eff is needed as the lagged weight for 1988.
sort naics sector assetcategory year
by naics sector assetcategory: ///
    gen double d_ln_P = ln(P_deflator) - ln(P_deflator[_n-1])
label variable d_ln_P "Log change in capital asset price deflator (asset)"

*--- A5b: Quality-adjusted prices for equipment and software --------
merge m:1 assetcategory using `cjk_bls', keep(master match) nogenerate
merge m:1 nipa_code year using `cjk_long', keep(master match) nogenerate
merge m:1 nipa_code using `cjk_gap', keep(master match) nogenerate
replace d_ln_P = d_ln_cjk if year <= 2016 & !missing(d_ln_cjk)
replace d_ln_P = d_ln_P + gap if year > 2016 & !missing(gap)
drop nipa_code d_ln_cjk gap
* The merges reorder the data; A6 takes lags by position.
sort naics sector assetcategory year

*--- A6: First-stage Törnqvist shares (asset level) -----------------
* Average share = 0.5 × (current share + lagged share)
by naics sector assetcategory: gen double w_eff_lag = w_eff[_n-1]

bysort naics sector year: egen double W_tot     = total(w_eff)
bysort naics sector year: egen double W_tot_lag = total(w_eff_lag)

gen double s_t   = w_eff     / W_tot
gen double s_lag = w_eff_lag / W_tot_lag
gen double s_bar = 0.5 * (s_t + s_lag)

gen double torn_a = s_bar * d_ln_P
label variable torn_a "Törnqvist contribution (asset, BLS)"

*--- A7: First-stage collapse: asset → naics×sector -----------------
* Keep year 1987 (d_ln_P_sec = 0 since torn_a missing → sum = 0);
* its w_sec value is needed as lagged weight in second stage.
* 1987 is dropped after second stage.
collapse (sum) d_ln_P_sec = torn_a  (sum) w_sec = w_eff, by(naics sector year)
label variable d_ln_P_sec "Törnqvist log ΔP (asset→sector, BLS)"
label variable w_sec      "Effective weight sum (sector-year, BLS)"

*--- A8: Map naics sectors → 61-industry classification -------------
rename naics naics_code
merge m:1 naics_code using "$conc/BLS_naics_map.dta", keep(master match) nogenerate

*--- A9: Second-stage Törnqvist: naics×sector → industry_id ---------
sort industry_id naics_code year
by industry_id naics_code: gen double w_sec_lag = w_sec[_n-1]

bysort industry_id year: egen double W_ind     = total(w_sec)
bysort industry_id year: egen double W_ind_lag = total(w_sec_lag)

gen double s_ind_t   = w_sec     / W_ind
gen double s_ind_lag = w_sec_lag / W_ind_lag
gen double s_bar_ind = 0.5 * (s_ind_t + s_ind_lag)

gen double torn_ind = s_bar_ind * d_ln_P_sec

collapse (sum) d_ln_P_new = torn_ind  (sum) W_new = w_sec  (firstnm) indname, ///
         by(industry_id year)

* Drop 1987: retained only to supply lagged weights for 1988
drop if year == 1987

label variable d_ln_P_new "Törnqvist log ΔP (BLS, industry-year, 1988–2024)"
label variable W_new      "Effective weight sum (industry-year, BLS)"

tempfile bls_prices
save `bls_prices'
display as result "Part A complete: BLS price changes, `=_N' obs."

*====================================================================
* PART B: BEA period (1948–1987) — asset → industry Törnqvist
*====================================================================

use "$interm/07_markups/BEA_Panel.dta", clear


*--- B1: Merge industry mapping (BEA sector → industry_id) ----------
merge m:1 sector using "$conc/BEA_sector_map.dta", keep(match) nogenerate

*--- B2: Merge IRR_it, ERR_t for weight construction ----------------
preserve
use "$interm/07_markups/IRR_ERR_1947to2024.dta", clear
keep industry_id year IRR_it ERR_t
tempfile irr
save `irr'
restore

merge m:1 industry_id year using `irr', keep(match) nogenerate

*--- B3: Effective weight at asset level (same rule as file 03 B3) --
gen double R_irr_kit = LP * IRR_it + P * delta - cap_gain_dollar
gen double R_err_kit = LP * ERR_t  + P * delta
gen double w_eff = cond(R_irr_kit > 0, R_irr_kit * stk2, R_err_kit * stk2)
* Floor at zero, as in A4.
replace w_eff = 0 if w_eff < 0
keep if year >= 1947

*--- B4: Asset-level log price change --------------------------------
* Year 1947: d_ln_P = missing (1946 absent after keep if year>=1947)
* Year 1947 row is retained to supply w_eff_lag for year 1948.
sort sector assetcode year
by sector assetcode: gen double d_ln_P = ln(P) - ln(P[_n-1])
by sector assetcode: gen double w_eff_lag = w_eff[_n-1]

*--- B4b: Quality-adjusted prices for equipment and software, as in A5b
merge m:1 assetcode using `cjk_bea', keep(master match) nogenerate
merge m:1 nipa_code year using `cjk_long', keep(master match) nogenerate
merge m:1 nipa_code using `cjk_gap', keep(master match) nogenerate
replace d_ln_P = d_ln_cjk if year <= 2016 & !missing(d_ln_cjk)
replace d_ln_P = d_ln_P + gap if year > 2016 & !missing(gap)
drop nipa_code d_ln_cjk gap

*--- B5: Törnqvist shares (asset → industry_id, one stage) ----------
bysort industry_id year: egen double W_tot     = total(w_eff)
bysort industry_id year: egen double W_tot_lag = total(w_eff_lag)

gen double s_t   = w_eff     / W_tot
gen double s_lag = w_eff_lag / W_tot_lag
gen double s_bar = 0.5 * (s_t + s_lag)
gen double torn_a = s_bar * d_ln_P

*--- B6: Collapse to industry × year; drop 1947 placeholder --------
* collapse (sum) returns 0 for all-missing (year 1947); mark as missing
* then drop — year 1947 served only to provide lagged weights for 1948.
collapse (sum) d_ln_P_it = torn_a  (firstnm) indname, ///
         by(industry_id year)
replace d_ln_P_it = . if year == 1947
drop if year == 1947
keep if year <= 1987

label variable d_ln_P_it "Törnqvist log ΔP (BEA, industry-year, 1948–1987)"

tempfile bea_prices
save `bea_prices'
display as result "Part B complete: BEA price changes, `=_N' obs."

*====================================================================
* PART C: Append BEA+BLS, industry → BroadSector Törnqvist, cumulate
*====================================================================

*--- C1: Combine BEA (1948–1987) and BLS (1988–2024) ---------------
* Data in memory is bea_prices (industry_id indname year d_ln_P_it)
preserve
use `bls_prices', clear
rename d_ln_P_new d_ln_P_it
keep industry_id indname year d_ln_P_it
tempfile bls_renamed
save `bls_renamed'
restore

keep industry_id indname year d_ln_P_it
append using `bls_renamed'
* Combined: 1948–2024, no splice needed — growth rates chain directly
label variable d_ln_P_it "Törnqvist log ΔP (industry-year, 1948–2024)"

*--- C2: Merge BroadSector from crosswalk Sheet3 -------------------
merge m:1 industry_id using "$conc/BroadSector_map.dta", keep(master match) nogenerate

*--- C3: Build capital compensation weights (1947–2023) -------------
* Identical construction to file 05 STEP 4

preserve
import excel "$raw/usa_wk_mar_2017.xlsx", ///
    sheet("KLEMdata") cellrange(A3) clear
rename (A B C D E F G H I J) ///
    (year industry go_nom cap_old lab_nom int_nom go_con cap_con lab_con int_con)
destring year industry cap_old, replace force
keep year industry cap_old
rename industry industry_id
keep if year <= 2014
tempfile capcomp_wk
save `capcomp_wk'
restore

preserve
import excel "$raw/EV_production_accounts_1947to2023.xlsx", ///
    sheet("Data_63Ind_1963to2023") firstrow clear
keep industry_id year CAP
destring industry_id year CAP, replace force
keep if industry_id <= 61
rename CAP cap_ev
tempfile capcomp_ev_all
save `capcomp_ev_all'
restore

preserve
use `capcomp_wk', clear
keep if year == 2014
rename cap_old cap_wk2014
keep industry_id cap_wk2014
tempfile anchor_wk
save `anchor_wk'
restore

preserve
use `capcomp_ev_all', clear
keep if year == 2014
rename cap_ev cap_ev2014
keep industry_id cap_ev2014
tempfile anchor_ev
save `anchor_ev'
restore

preserve
use `capcomp_ev_all', clear
keep if year >= 2015
merge m:1 industry_id using `anchor_wk', nogenerate keep(match)
merge m:1 industry_id using `anchor_ev', nogenerate keep(match)
gen double CC_it = cap_ev * (cap_wk2014 / cap_ev2014)
keep industry_id year CC_it
tempfile capcomp_post
save `capcomp_post'
restore

preserve
use `capcomp_wk', clear
rename cap_old CC_it
append using `capcomp_post'
* CC_it now covers 1947–2023
tempfile CC_weights
save `CC_weights'
restore

* Lagged CC: shift year by +1 so CC_1947 becomes CC_lag for year 1948, etc.
preserve
use `CC_weights', clear
rename CC_it CC_lag
replace year = year + 1
tempfile CC_lag_weights
save `CC_lag_weights'
restore

*--- C4: Merge CC current and lagged --------------------------------
* keep(match) naturally drops year 2024 (CC ends 2023) — consistent with step 05
merge 1:1 industry_id year using `CC_weights',     keep(match)        nogenerate
merge 1:1 industry_id year using `CC_lag_weights', keep(master match) nogenerate
label variable CC_it  "Capital compensation (current year)"
label variable CC_lag "Capital compensation (lagged year)"

*--- C5: Stage-2 Törnqvist: industry → BroadSector -----------------
bysort broadsector year: egen double CC_tot     = total(CC_it)
bysort broadsector year: egen double CC_tot_lag = total(CC_lag)

gen double s_CC_t   = CC_it  / CC_tot
gen double s_CC_lag = CC_lag / CC_tot_lag
gen double s_bar_CC = 0.5 * (s_CC_t + s_CC_lag)

gen double torn_bs = s_bar_CC * d_ln_P_it

collapse (sum) d_ln_P_bs = torn_bs, by(broadsector year)
label variable d_ln_P_bs "Törnqvist log ΔP (broad sector-year, 1948–2023)"

*--- C6: Add base year 1947 (P_bs = 1), cumulate -------------------
preserve
keep broadsector
duplicates drop
gen year = 1947
gen double d_ln_P_bs = 0
tempfile base47
save `base47'
restore

append using `base47'
sort broadsector year

* Running sum of log changes; base 1947 contributes 0 → P_bs_1947 = exp(0) = 1
by broadsector: gen double ln_P_bs = sum(d_ln_P_bs)
gen double P_bs = exp(ln_P_bs)

label variable broadsector "Broad sector (GOODS / RI / NRI)"
label variable year        "Year"
label variable d_ln_P_bs   "Törnqvist log change in capital price (broad sector)"
label variable ln_P_bs     "Cumulated log capital price (base 1947 = 0)"
label variable P_bs        "Capital price index (base 1947 = 1)"

*====================================================================
* Deflate capital price by aggregate VA price (from step 04)
*====================================================================

preserve
    use "$interm/04_value_added_prices.dta", clear
    rename Year year
    keep year VAP_AGGREGATE
    tempfile vap_agg
    save `vap_agg'
restore

merge m:1 year using `vap_agg', keep(master match) nogenerate
replace P_bs = P_bs / VAP_AGGREGATE
label variable P_bs "Capital price index deflated by aggregate VA price (base 1947 = 1)"
drop VAP_AGGREGATE

*====================================================================
* SAVE  (wide format, Year capitalised, to match 12_compilation.do)
*====================================================================

drop d_ln_P_bs ln_P_bs
reshape wide P_bs, i(year) j(broadsector) string
rename P_bsGOODS CAPPRICE_GOODS
rename P_bsRI    CAPPRICE_RI
rename P_bsNRI   CAPPRICE_NRI
rename year Year
sort Year
save "$interm/08_capital_prices.dta", replace
display as result "Saved 08_capital_prices.dta (`=_N' obs)."

*=======================================================================================
*=======================================================================================
