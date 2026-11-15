
*=======================================================================================
* An Anatomy of the Baumol Cost Disease -- Elbadri & Vizcaino (2026)
*---------------------------------------------------------------------------------------
* File:    09_demographics.do
* Purpose: Extract US population levels and growth rate from Penn World Tables 11.0
*---------------------------------------------------------------------------------------
* Inputs:  input/pwt110.dta  -- Penn World Tables v11.0
*
* Outputs: intermediate/09_demographics.dta
*            pop    -- US population (millions)
*            pop_g  -- log population growth rate
*=======================================================================================

set varabbrev off

*=======================================================================================
* STEP 1: Extract population level and growth rate
*=======================================================================================

use "$raw/pwt110.dta", clear
keep if countrycode == "USA"
keep year pop
tsset year
gen pop_g = log(pop) - log(L.pop)
rename year Year

save "$interm/09_demographics.dta", replace

*=======================================================================================
*=======================================================================================
