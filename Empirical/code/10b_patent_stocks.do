*=======================================================================================
* An Anatomy of the Baumol Cost Disease -- Elbadri & Vizcaino (2026)
*---------------------------------------------------------------------------------------
* File:    10b_patent_stocks.do
* Purpose: Sectoral patent stocks used by the model (knowledge stocks and 1950 N-bar)
*---------------------------------------------------------------------------------------
* Inputs:  ../Cross_Citations/output/patents/ex_cpc_hist_stocks.csv
*            Aggregation ex_cpc_hist from Cross_Citations/code/05_patent_stocks.py:
*            utility grants assigned to sectors through IPC subclasses (the partition of the
*            citation matrix), business-method patents (current CPC G06Q) excluded, accumulated
*            without obsolescence; before 1976 the sector composition moves with the USPC-based
*            shares of the USPTO Historical Patent Data Files.
*
* Outputs: intermediate/10b_patent_stocks.dta
*            PATSTOCK_{sector}    -- patent stock (number of patents)
*            PATSTOCKSH_{sector}  -- sector share of the total stock; the 1950 row is N-bar
*=======================================================================================

set varabbrev off

import delimited "$crosscit/ex_cpc_hist_stocks.csv", clear varnames(1) asdouble
replace sector = "GOODS" if sector == "Goods"
replace sector = "RI"    if sector == "RI_Services"
replace sector = "NRI"   if sector == "NRI_Services"
rename year Year
keep if Year <= 2023
reshape wide n, i(Year) j(sector) string
rename n* PATSTOCK_*

gen double PATSTOCK_TOTAL = PATSTOCK_GOODS + PATSTOCK_RI + PATSTOCK_NRI
foreach x in GOODS RI NRI {
	gen double PATSTOCKSH_`x' = PATSTOCK_`x' / PATSTOCK_TOTAL
}
drop PATSTOCK_TOTAL
keep Year PATSTOCK_* PATSTOCKSH_*
save "$interm/10b_patent_stocks.dta", replace
