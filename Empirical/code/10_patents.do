*=======================================================================================
* An Anatomy of the Baumol Cost Disease -- Elbadri & Vizcaino (2026)
*---------------------------------------------------------------------------------------
* File:    10_patents.do
* Purpose: Sector shares of issued patents by year (Section 2): the IPC route used by the
*          model's knowledge stocks, and the USPC route kept for comparison
*---------------------------------------------------------------------------------------
* Inputs:  input/ALP_USPC/NAICS_2007/uspc_Class_to_naics07_2.txt   -- USPC-NAICS 2007 (2-digit) crosswalk (Goldschlag, Lybbert & Zolas, 2016)
*          input/historical_masterfile.dta                         -- USPTO historical patent masterfile (https://www.uspto.gov/ip-policy/economic-research/research-datasets/historical-patent-data-files)
*          intermediate/concordances/USPTO_2digits.dta             -- NAICS 2-digit to sector concordance
*          ../Cross_Citations/output/patents/ex_cpc_hist_flow_shares.csv -- IPC-route sector shares of grants,
*            1947-2023: ex_cpc_hist in Cross_Citations/code/05_patent_stocks.py (IPC
*            subclasses, the partition of the citation matrix; business-method patents, current
*            CPC G06Q, excluded; before 1976 the composition moves with the USPC route)
*
* Outputs: intermediate/10_patents.dta
*            PATSH_{sector}       -- share of issued patents, IPC route (the series in Section 2)
*            PAT_USPC_{sector}    -- probability-weighted patent counts, USPC route, by issue year
*            PATSH_USPC_{sector}  -- shares on the USPC route (business methods included)
*=======================================================================================

set varabbrev off

*=======================================================================================
* STEP 1: Pre-save USPC-NAICS 2007 (2-digit) crosswalk
*=======================================================================================

import delimited "$raw/ALP_USPC/NAICS_2007/uspc_Class_to_naics07_2.txt", clear
save "$interm/10_patents/uspc_Class_to_naics07_2.dta", replace

*=======================================================================================
* STEP 2: Aggregate issued patents by USPC class and year
*=======================================================================================

* Load historical patent masterfile (USPTO: https://www.uspto.gov/ip-policy/economic-research/research-datasets/historical-patent-data-files)
use "$raw/historical_masterfile.dta", clear

* Keep issued patents only and extract issue year
keep if disp_ty == "ISS"
keep uspc disp_dt
gen year = year(disp_dt)
drop disp_dt

* Destring USPC class code
gen uspc_1 = real(uspc)
drop uspc
rename uspc_1 uspc
drop if year == .
drop if uspc == .

* Collapse to USPC-year counts, then reshape wide for merging
gen total = 1
collapse (sum) total, by(year uspc)
rename total _
reshape wide _, i(uspc) j(year)

*=======================================================================================
* STEP 3: Apply probability weights and aggregate to NAICS 2-digit level
*=======================================================================================

* Merge with USPC-NAICS crosswalk
merge 1:m uspc using "$interm/10_patents/uspc_Class_to_naics07_2.dta"
keep if _merge == 3
drop _merge
sort uspc

* Reshape back to time series format and compute probability-weighted counts
reshape long _, i(uspc naics07_2 probability_weight) j(year)
sort year uspc
drop if _ == .
rename _ total
gen patents = probability_weight * total
keep patents naics07_2 year
collapse (sum) patents, by(year naics07_2)

*=======================================================================================
* STEP 4: Map to broad sectors using concordance and save
*=======================================================================================

* Reshape wide by year for concordance merge
xtset naics07_2 year
reshape wide patents, i(naics07_2) j(year)
rename patents* _*

* Merge with sector concordance
tostring naics07_2, gen(naics)
drop naics07_2
rename naics naics07_2
merge 1:1 naics07_2 using "$conc/USPTO_2digits.dta"
* Private economy: public administration (NAICS 92) is excluded, as in every other series.
drop if naics07_2 == "92"
collapse (sum) _1790-_2014, by(SectorCode)

* Reshape to panel format and save
reshape long _, i(SectorCode) j(Year)
drop if Year < 1947
reshape wide _, i(Year) j(SectorCode) string
rename _* PAT_*
local sect "GOODS RI NRI"
foreach x of local sect {
	gen PATSH_`x' = PAT_`x'/(PAT_GOODS+PAT_RI+PAT_NRI)
}
keep Year PATSH* PAT_*
foreach x of local sect {
	rename PAT_`x' PAT_USPC_`x'
	rename PATSH_`x' PATSH_USPC_`x'
}
tempfile uspc
save `uspc'

* IPC route: the sector assignment of the model's knowledge stocks, so the Section 2 shares and
* the stocks share one construction and both exclude business methods.
import delimited "$crosscit/ex_cpc_hist_flow_shares.csv", clear varnames(1) asdouble
replace sector = "GOODS" if sector == "Goods"
replace sector = "RI"    if sector == "RI_Services"
replace sector = "NRI"   if sector == "NRI_Services"
rename year Year
reshape wide share, i(Year) j(sector) string
rename share* PATSH_*
merge 1:1 Year using `uspc', nogenerate
sort Year
save "$interm/10_patents.dta", replace

*=======================================================================================
*=======================================================================================

















