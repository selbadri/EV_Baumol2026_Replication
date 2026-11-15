*=======================================================================================
* An Anatomy of the Baumol Cost Disease -- Elbadri & Vizcaino (2026)
*---------------------------------------------------------------------------------------
* File:    07b_rates_of_return.do
* Purpose: Industry internal rates of return (IRR) and the external rate (ERR), 1947-2024,
*          61 private industries, used by steps 07c and 08
*---------------------------------------------------------------------------------------
* Inputs:  intermediate/07_markups/BEA_Panel.dta              -- from 07a
*          input/usa_wk_mar_2017.xlsx                         -- US KLEMS capital compensation, 1947-2016
*          input/rental-prices-major-industries.xlsx          -- BLS capital detail, 1987-2024
*          input/pwt110.dta                                   -- PWT 11.0 irr, USA (external rate)
*          intermediate/concordances/BEA_sector_map.dta       -- BEA sector -> industry_id
*          intermediate/concordances/BLS_naics_map.dta        -- BLS naics  -> industry_id
* Outputs: intermediate/07_markups/IRR_ERR_1947to2024.dta
*            industry_id indname year
*            IRR_it  -- internal rate; negative values replaced by ERR_t
*            ERR_t   -- external rate (PWT 11.0 irr, USA)
*            IRR_old -- BEA-based internal rate, longest series
*---------------------------------------------------------------------------------------
* Method (BLS Handbook of Methods)
*   SV_kit = k x tau^P x P                              (stock value)
*   OC_kit = k x [tau^P x (P x delta - Pdot) + P x tau^L] (ownership cost; Pdot = MA3_P - L.MA3_P)
*   IRR_it = (CC_it - OC_it) / SV_it                    (solved from capital compensation)
*   Before 1987 tau^P = 1, tau^L = 0. From 1987 the BLS internal rate is used directly;
*   before 1987 the BEA-based rate is ratio-spliced to the BLS level in 1987.
*=======================================================================================

*====================================================================
* PART A: BEA-based OC, SV, CC → 61-industry level (1947–2016)
*====================================================================

use "$interm/07_markups/BEA_Panel.dta", clear

*--- A1: Per-asset ownership cost and stock value (τ^P=1, τ^L=0) ----
* OC_kit = stk2 × (P×δ − Ṗ)       Ṗ = cap_gain_dollar from BEA_Panel
* SV_kit = stk1  (= stk2 × P)
gen double OC_kit = stk2 * (P * delta - cap_gain_dollar)
gen double SV_kit = stk1
label variable OC_kit "Per-asset OC: stk2×(P×δ−Ṗ)  [τ^P=1,τ^L=0]"
label variable SV_kit "Per-asset SV: stk1 (= stk2×P) [τ^P=1,τ^L=0]"

*--- A2: Map BEA sectors to 61-industry classification ---------------
merge m:1 sector using "$conc/BEA_sector_map.dta", keep(match) nogenerate

*--- A3: Collapse to industry × year ---------------------------------
keep if year >= 1947
collapse (sum) OC_old = OC_kit  SV_old = SV_kit, by(industry_id indname year)
label variable OC_old "BEA OC (τ^P=1,τ^L=0), Σ_k stk2×(P×δ−Ṗ)"
label variable SV_old "BEA SV (τ^P=1,τ^L=0), Σ_k stk2×P"

*--- A4: Load BEA KLEM nominal capital costs (usa_wk_mar_2017) -------
* cap_old = nominal capital costs, 1947–2016
preserve
import excel "$raw/usa_wk_mar_2017.xlsx", ///
    sheet("KLEMdata") cellrange(A3) clear
rename (A B C D E F G H I J) ///
    (year industry go_nom cap_old lab_nom int_nom go_con cap_con lab_con int_con)
destring year industry cap_old, replace force
keep year industry cap_old
rename industry industry_id
label variable cap_old "Nominal capital costs, BEA KLEM (usa_wk_mar_2017)"
tempfile klem
save `klem'
restore

merge 1:1 industry_id year using `klem', keep(master match) nogenerate

*--- A5: BEA-based IRR (longest available series) --------------------
* Available wherever cap_old, OC_old, and SV_old are all non-missing.
* cap_old ends at 2016; cap_gain_dollar unavailable at boundary years.
gen double IRR_old = (cap_old - OC_old) / SV_old
label variable IRR_old "BEA-based IRR = (cap_old − OC_old) / SV_old  [longest series]"

tempfile bea_series
save `bea_series'
display as result "Part A complete: BEA series, `=_N' obs."

*====================================================================
* PART B: BLS IRR → 61-industry level (1987–2024)
* Uses published IRR, aggregated with weight w = R_IRR × productive stock.
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

*--- B1: Keep IRR and components for weighting -----------------------
keep if inlist(measure, ///
    "Internal rate of return", ///
    "Rental price with internal rate of return", ///
    "Productive capital stock")
drop units

reshape long _, i(naics sector measure assetcategory) j(year)
rename _ value

* Sector-year IRR scalar
gen double irr_raw = value if measure == "Internal rate of return"
bysort naics sector year: egen double irr_bls = max(irr_raw)
drop irr_raw

* Note: `irr_bls` is already the BLS IRR for a given sector-year.
* Nothing below is used to recompute sector-year IRR.
* This block only builds a sector-year weight, `w_sec`, so that when
* multiple BLS sectors map into the same `industry_id`, the final
* industry-year IRR can be formed as a weighted average across those
* sector-year IRR values.
* Requested sector-year weight: w_sec = Σ_assets(R_IRR × k_prod)
preserve
keep if inlist(measure, "Rental price with internal rate of return", "Productive capital stock")

gen byte mnum = .
replace mnum = 1 if measure == "Rental price with internal rate of return"
replace mnum = 2 if measure == "Productive capital stock"
drop measure

reshape wide value, i(naics sector assetcategory year) j(mnum)
rename value1 R_irr
rename value2 k_prod

gen double w_irrk = R_irr * k_prod
collapse (sum) w_sec = w_irrk, by(naics sector year)

tempfile bls_w
save `bls_w'
restore

keep naics sector year irr_bls
duplicates drop
merge 1:1 naics sector year using `bls_w', keep(master match) nogenerate

*--- B2: Map to 61-industry classification --------------------------
rename naics naics_code
merge m:1 naics_code using "$conc/BLS_naics_map.dta", keep(master match) nogenerate

gen double irr_w = irr_bls * w_sec
collapse (sum) irr_w w_sec (firstnm) indname, by(industry_id year)
gen double irr_new = irr_w / w_sec
drop irr_w
label variable irr_new "BLS-published IRR, wtd by R_IRR×productive stock (industry-year, 1987+)"

tempfile bls_series
save `bls_series'
display as result "Part B complete: BLS IRR, `=_N' obs."

*====================================================================
* PART C: Splice IRR at 1987
*====================================================================

*--- C1: Merge BEA and BLS series ------------------------------------
use `bea_series', clear
merge 1:1 industry_id year using `bls_series', nogenerate

*--- C2: Ratio-splice IRR_old to BLS level at 1987 ------------------
* 1987+  : BLS irr_new directly (authoritative).
* pre-1987: IRR_old carries BEA time variation; rescaled to BLS level
*   at 1987 so the series is continuous at the splice point.
preserve
keep if year == 1987
rename IRR_old  IRR_old87
rename irr_new  irr_new87
keep industry_id IRR_old87 irr_new87
tempfile anch87
save `anch87'
restore

merge m:1 industry_id using `anch87', nogenerate keep(master match)

gen double IRR_it = irr_new                              if year >= 1987 & !missing(irr_new)
replace    IRR_it = irr_new87 * (IRR_old / IRR_old87)   if year <  1987 & !missing(IRR_old87) & IRR_old87 != 0
drop IRR_old87 irr_new87
label variable IRR_it "IRR: BLS 1987+; BEA ratio-spliced at 1987 pre-1987"

display as result "Part C complete: spliced series computed."

*====================================================================
* PART D: ERR from PWT 11.0 and replacement
*====================================================================

*--- D1: External rate of return: PWT 11.0 irr, United States ---------
preserve 
use "$raw/pwt110.dta", clear
keep countrycode year irr
keep if countrycode == "USA"
rename irr ERR_t
tempfile err_series
save `err_series'
restore

merge m:1 year using `err_series', nogenerate keep(master match)

*--- D2: Replace negative IRR with ERR ------------------------------
replace IRR_it = ERR_t if IRR_it < 0 & !missing(ERR_t)

display as result "Part D complete: `=sum(IRR_it == ERR_t & !missing(ERR_t))' obs replaced with ERR."

*====================================================================
* SAVE OUTPUT
*====================================================================

keep industry_id indname year IRR_it ERR_t IRR_old
label variable industry_id "Industry identifier (61-industry classification)"
label variable indname      "Industry description"
label variable year         "Year"
label variable IRR_it       "Internal rate of return (cleaned: negative IRR replaced by ERR_t)"
label variable ERR_t        "External rate of return (PWT 11.0 irr, USA)"
label variable IRR_old      "BEA-based IRR, longest series (cap_old−OC_old)/SV_old"

sort industry_id year
save "$interm/07_markups/IRR_ERR_1947to2024.dta", replace
display as result "Saved IRR_ERR_1947to2024.dta (`=_N' obs, `=r(k)' variables)."
