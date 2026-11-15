*=======================================================================================
* An Anatomy of the Baumol Cost Disease -- Elbadri & Vizcaino (2026)
*---------------------------------------------------------------------------------------
* File:    07a_capital_assets.do
* Purpose: BEA capital asset panel by industry and asset, 1946-2024: asset price, depreciation
*          rate and capital gains, used by steps 07b, 07c and 08
*---------------------------------------------------------------------------------------
* Inputs:  input/detailnonres_stk1.xlsx  -- BEA current-cost net stock (millions of $)
*          input/detailnonres_stk2.xlsx  -- BEA chained 2017-$ net stock
*          input/detailnonres_inv2.xlsx  -- BEA chained 2017-$ gross investment
* Outputs: intermediate/07_markups/BEA_Panel.dta
*            one row per BEA sector x asset x year: stk1 stk2 inv2,
*            P     = stk1 / stk2 (asset price deflator, 2017 = 1), LP = L.P,
*            delta = depreciation rate (from inv2, trimmed to [0,1]),
*            MA3_P = (L.P + P + F.P) / 3, cap_gain_dollar = MA3_P - L.MA3_P
* Note:    tax parameters tau^P = 1, tau^L = 0 before 1987 (no BLS tax data); ownership
*          costs and stock values are computed in 07b.
*=======================================================================================

capture mkdir "$interm/07_markups"

local stk1_file "$raw/detailnonres_stk1.xlsx"
local stk2_file "$raw/detailnonres_stk2.xlsx"
local inv2_file "$raw/detailnonres_inv2.xlsx"

*====================================================================
* PHASE 1: Import all sheets from each file and stack into one dataset
*====================================================================
* Sheet layout (each sheet = one BEA industry, e.g. "110C - Farms"):
*   Row 1  : sector label
*   Row 2  : description
*   Row 3  : units
*   Row 4  : update date
*   Row 5  : blank
*   Row 6  : column headers — "Asset Codes" | "NIPA Asset Types" | 1925 | … | 2024
*   Row 7  : blank spacer
*   Rows 8+: data rows
*
* Import strategy: cellrange(A6) so row 1 of imported data = the header
*   row. Read actual year values from that row; rename only 4-digit-year
*   columns to _YYYY; drop all other columns.

local measure_stk1 "Current-cost net stock"
local measure_stk2 "Chained 2017-dollar net stock"
local measure_inv2 "Chained 2017-dollar gross investment"
local units_stk1   "Millions of dollars"
local units_stk2   "Millions of chained 2017 dollars"
local units_inv2   "Millions of chained 2017 dollars"

tempfile master_all
local is_first_all 1

foreach stub in stk1 stk2 inv2 {

    local file "``stub'_file'"

    quietly import excel "`file'", describe
    local nsheets = r(N_worksheet)
    forvalues s = 1/`nsheets' {
        local shname_`s' = r(worksheet_`s')
    }

    forvalues s = 1/`nsheets' {
        local sh "`shname_`s''"

        if regexm(lower("`sh'"), "readme|read me|documentation|notes") {
            display as text "  Skipping sheet: `sh'"
            continue
        }

        import excel "`file'", sheet("`sh'") cellrange(A6) clear allstring

        ds A B, not
        foreach v of varlist `r(varlist)' {
            local hdr = strtrim(`v'[1])
            if regexm("`hdr'", "^[0-9][0-9][0-9][0-9]$") {
                rename `v' _`hdr'
            }
            else {
                drop `v'
            }
        }

        drop in 1/2

        drop if missing(A) | strtrim(A) == ""

        rename A assetcode
        rename B assetname

        capture ds _*
        if _rc == 0 {
            foreach v of varlist `r(varlist)' {
                capture destring `v', replace force ignore(",")
            }
        }
        else {
            display as text "  Warning: no year columns in sheet `sh' (`stub') — skipping."
            continue
        }

        gen str20 sector  = "`sh'"
        gen str60 measure = "`measure_`stub''"
        gen str50 units   = "`units_`stub''"

        if `is_first_all' {
            save `master_all', replace
            local is_first_all 0
        }
        else {
            append using `master_all'
            save `master_all', replace
        }
    }

    display as result "Phase 1 complete: `stub' — `nsheets' sheets appended."
}

use `master_all', clear
order sector assetcode assetname measure units
display as result "Phase 1 total: `=_N' observations."

*====================================================================
* PHASE 2: Clean, reshape, compute P, LP, δ, MA3_P, cap_gain_dollar
*====================================================================

*--- Step 2a: Drop aggregate rows ------------------------------------
drop if inlist(upper(strtrim(assetcode)), ///
    "EQUIPMENT", "STRUCTURES", "IPP", "ES00")

*--- Step 2b: Reshape wide-year → long --------------------------------
reshape long _, i(sector assetcode assetname measure units) j(year)
rename _ value

*--- Step 2c: Pivot measures from rows to columns --------------------
gen byte mnum = .
replace mnum = 1 if measure == "Current-cost net stock"
replace mnum = 2 if measure == "Chained 2017-dollar net stock"
replace mnum = 3 if measure == "Chained 2017-dollar gross investment"
drop measure units
reshape wide value, i(sector assetcode assetname year) j(mnum)
rename value1 stk1
rename value2 stk2
rename value3 inv2
label variable stk1 "Current-cost net stock (millions $)"
label variable stk2 "Chained 2017-dollar net stock (millions $)"
label variable inv2 "Chained 2017-dollar gross investment (millions $)"

*--- Step 2d: Panel setup --------------------------------------------
encode sector,    gen(sectorcode)
encode assetcode, gen(assetcd)
egen panelid = group(sectorcode assetcd)
xtset panelid year

*--- Step 2e: Capital asset price deflator ---------------------------
gen double P = stk1 / stk2
label variable P "Capital asset price deflator (= 1.0 in 2017)"

*--- Step 2f: Lagged price (for the rental price formula) -
gen double LP = L.P
label variable LP "Lagged capital asset price deflator (P_{t-1})"

*--- Step 2g: Wealth stock depreciation rate -------------------------
gen double delta = (inv2 - (stk2 - L.stk2)) / ((stk2 + L.stk2) / 2)
replace delta = 0 if delta <  0 & !missing(delta)
replace delta = 1 if delta >  1 & !missing(delta)
label variable delta "Wealth stock depreciation rate (inv2-based, trimmed [0,1])"

*--- Step 2h: 3-year centred MA of P and dollar capital gain ---------
* BLS HOM: capital gain = "change in the 3-year moving average of
*   the asset price" (dollar change, not a rate).
gen double MA3_P          = (L.P + P + F.P) / 3
gen double cap_gain_dollar = MA3_P - L.MA3_P
label variable MA3_P          "3-yr centred MA of price deflator"
label variable cap_gain_dollar "Capital gain: dollar change in 3-yr centred MA of P (Ṗ)"

*--- Step 2i: Trim to 1945+ -----------------------------------------------
* 1946 needed so L.P is valid at 1947.
* 1945 needed so L.MA3_P (= MA3_P_1946 = (P_1945+P_1946+P_1947)/3) is
* valid at 1947, giving a non-missing cap_gain_dollar in 1947.
keep if year >= 1945

*--- Step 2j: Drop panel bookkeeping variables ----------------------
drop sectorcode assetcd panelid

order sector assetcode assetname year stk1 stk2 inv2 P LP delta MA3_P cap_gain_dollar
sort sector assetcode year

save "$interm/07_markups/BEA_Panel.dta", replace
display as result "Phase 2 complete: BEA_Panel.dta saved (`=_N' obs, `=r(k)' vars)."
