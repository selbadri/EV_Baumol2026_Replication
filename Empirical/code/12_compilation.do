
*=======================================================================================
* An Anatomy of the Baumol Cost Disease -- Elbadri & Vizcaino (2026)
*---------------------------------------------------------------------------------------
* File:    12_compilation.do
* Purpose: Merge all constructed datasets; export the data moments
*---------------------------------------------------------------------------------------
* Inputs:  intermediate/02_labor_shares.dta
*          intermediate/03_value_added_shares.dta
*          intermediate/04_value_added_prices.dta
*          intermediate/05_labor_productivity_growth.dta
*          intermediate/06_labor_quality_growth.dta
*          intermediate/07_markups.dta
*          intermediate/08_capital_prices.dta
*          intermediate/09_demographics.dta
*          intermediate/10_patents.dta
*          intermediate/10b_patent_stocks.dta
*          intermediate/11_scientists.dta
*
* Outputs: output/dataset/data_moments.xlsx   -- compiled data moments table
*                                             (DataMoments sheet: data; README sheet: variable descriptions)
*          output/dataset/data_moments_markups_noRE.xlsx, data_moments_markups_BEA.xlsx
*                                          -- the same table with the markup series of the two
*                                             checks of Appendix F (intermediate/07_markups_*.dta)
*=======================================================================================

set varabbrev off

*=======================================================================================
* STEP 1: Merge all datasets together
*=======================================================================================

use "$interm/02_labor_shares.dta", clear

local mergesteps "03_value_added_shares 04_value_added_prices 05_labor_productivity_growth 06_labor_quality_growth 07_markups 08_capital_prices 09_demographics 10_patents 10b_patent_stocks 11_scientists"
foreach f of local mergesteps {
	merge 1:1 Year using "$interm/`f'.dta", nogenerate
}

*=======================================================================================
* STEP 2: Export all data moments to data_moments
*=======================================================================================

drop dln_VAP // no longer needed
order Year LABVASHARE_* VASHARE* VAP* CAPPRICE* LP* *LABQLTY* profit_share* pop* PAT_* PATSH* PATSTOCK_* PATSTOCKSH_* RDSH* SCI_SH

foreach var of varlist LABVASHARE_* VASHARE* VAP* CAPPRICE* LP* *LABQLTY* profit_share* pop* PAT_* PATSH* PATSTOCK_* PATSTOCKSH_* RDSH* SCI_SH {
    format `var' %9.2f
}
export excel using "$output/dataset/data_moments.xlsx", replace firstrow(variables) sheet("DataMoments") 

*=======================================================================================
* STEP 3: Write README sheet to data_moments.xlsx
*=======================================================================================

putexcel set "$output/dataset/data_moments.xlsx", sheet("README") modify

* --- Sector abbreviations ---
putexcel A1 = "SECTOR ABBREVIATIONS", bold
putexcel A2 = "Abbreviation"          B2 = "Full Name"
putexcel A3 = "GOODS"                 B3 = "Goods"
putexcel A4 = "RI"                    B4 = "Research Intensive Services"
putexcel A5 = "NRI"                   B5 = "Non-Research Intensive Services"

* --- Variables ---
putexcel A7 = "VARIABLES", bold
putexcel A8  = "Variable"              B8  = "Description"
putexcel A9  = "Year"                  B9  = "Calendar year"
putexcel A10 = "LABVASHARE_GOODS"      B10 = "Labor share of value added -- Goods"
putexcel A11 = "LABVASHARE_RI"         B11 = "Labor share of value added -- Research Intensive Services"
putexcel A12 = "LABVASHARE_NRI"        B12 = "Labor share of value added -- Non-Research Intensive Services"
putexcel A13 = "VASHARE_GOODS"         B13 = "Share of total nominal value added -- Goods"
putexcel A14 = "VASHARE_RI"            B14 = "Share of total nominal value added -- Research Intensive Services"
putexcel A15 = "VASHARE_NRI"           B15 = "Share of total nominal value added -- Non-Research Intensive Services"
putexcel A16 = "VAPRICE_GOODS"         B16 = "Relative sectoral value added price (1947 = 1.00) -- Goods"
putexcel A17 = "VAPRICE_RI"            B17 = "Relative sectoral value added price (1947 = 1.00) -- Research Intensive Services"
putexcel A18 = "VAPRICE_NRI"           B18 = "Relative sectoral value added price (1947 = 1.00) -- Non-Research Intensive Services"
putexcel A19 = "CAPPRICE_GOODS"        B19 = "Relative capital price index (1947 = 1.00) -- Goods"
putexcel A20 = "CAPPRICE_RI"           B20 = "Relative capital price index (1947 = 1.00) -- Research Intensive Services"
putexcel A21 = "CAPPRICE_NRI"          B21 = "Relative capital price index (1947 = 1.00) -- Non-Research Intensive Services"
putexcel A22 = "LP_GOODS"              B22 = "Labor productivity index (1947 = 1.00) -- Goods"
putexcel A23 = "LP_RI"                 B23 = "Labor productivity index (1947 = 1.00) -- Research Intensive Services"
putexcel A24 = "LP_NRI"                B24 = "Labor productivity index (1947 = 1.00) -- Non-Research Intensive Services"
putexcel A25 = "LP_AGGREGATE"          B25 = "Labor productivity index (1947 = 1.00) -- Aggregate"
putexcel A26 = "LABQLTY_GOODS"         B26 = "Labor quality index (1947 = 1.00) -- Goods"
putexcel A27 = "LABQLTY_RI"            B27 = "Labor quality index (1947 = 1.00) -- Research Intensive Services"
putexcel A28 = "LABQLTY_NRI"           B28 = "Labor quality index (1947 = 1.00) -- Non-Research Intensive Services"
putexcel A29 = "profit_share_GOODS"    B29 = "Machine-producer profit share (produced nonresidential assets; the model uses a quadratic trend of the gross markup 1/(1-share) over 1950-2023) -- Goods"
putexcel A30 = "profit_share_RI"       B30 = "Machine-producer profit share (produced nonresidential assets; the model uses a quadratic trend of the gross markup 1/(1-share) over 1950-2023) -- Research Intensive Services"
putexcel A31 = "profit_share_NRI"      B31 = "Machine-producer profit share (produced nonresidential assets; the model uses a quadratic trend of the gross markup 1/(1-share) over 1950-2023) -- Non-Research Intensive Services"
putexcel A32 = "PAT_USPC_GOODS"        B32 = "Number of issued patents, USPC route -- Goods"
putexcel A33 = "PAT_USPC_RI"           B33 = "Number of issued patents, USPC route -- Research Intensive Services"
putexcel A34 = "PAT_USPC_NRI"          B34 = "Number of issued patents, USPC route -- Non-Research Intensive Services"
putexcel A35 = "pop"                   B35 = "US population (millions)"
putexcel A36 = "pop_g"                 B36 = "US log population growth rate"
putexcel A37 = "PATSH_GOODS"           B37 = "Share of issued patents, IPC route, business methods excluded -- Goods"
putexcel A38 = "PATSH_RI"              B38 = "Share of issued patents, IPC route, business methods excluded -- Research Intensive Services"
putexcel A39 = "PATSH_NRI"             B39 = "Share of issued patents, IPC route, business methods excluded -- Non-Research Intensive Services"
putexcel A40 = "RDSH_GOODS"            B40 = "R&D scientist share of total R&D scientists -- Goods"
putexcel A41 = "RDSH_RI"               B41 = "R&D scientist share of total R&D scientists -- Research Intensive Services"
putexcel A42 = "RDSH_NRI"              B42 = "R&D scientist share of total R&D scientists -- Non-Research Intensive Services"
putexcel A43 = "SCI_SH"                B43 = "R&D scientists as share of total employment (private industries)"
putexcel A44 = "PATSTOCK_GOODS"        B44 = "Patent stock, no obsolescence, IPC route, business methods excluded -- Goods"
putexcel A45 = "PATSTOCK_RI"           B45 = "Patent stock, no obsolescence, IPC route, business methods excluded -- Research Intensive Services"
putexcel A46 = "PATSTOCK_NRI"          B46 = "Patent stock, no obsolescence, IPC route, business methods excluded -- Non-Research Intensive Services"
putexcel A47 = "PATSTOCKSH_GOODS"      B47 = "Share of the total patent stock -- Goods (1950 value = initial knowledge share N-bar)"
putexcel A48 = "PATSTOCKSH_RI"         B48 = "Share of the total patent stock -- Research Intensive Services"
putexcel A49 = "PATSTOCKSH_NRI"        B49 = "Share of the total patent stock -- Non-Research Intensive Services"
putexcel A50 = "PATSH_USPC_GOODS"      B50 = "Share of issued patents, USPC route (comparison) -- Goods"
putexcel A51 = "PATSH_USPC_RI"         B51 = "Share of issued patents, USPC route (comparison) -- Research Intensive Services"
putexcel A52 = "PATSH_USPC_NRI"        B52 = "Share of issued patents, USPC route (comparison) -- Non-Research Intensive Services"
putexcel A53 = "VAP_AGGREGATE"         B53 = "Private-economy value added price index, Tornqvist (1947 = 1.00)"

*=======================================================================================
* STEP 4: Data moments for the two markup checks of Appendix F
*         Same table, with the markup series built without real estate (footnote 41)
*         and from BEA data in every year (footnote 42); read by Calibration/code/robustness.jl
*=======================================================================================

foreach v in noRE BEA {
    preserve
    drop profit_share_*
    merge 1:1 Year using "$interm/07_markups_`v'.dta", nogenerate keepusing(profit_share_*)
    order Year LABVASHARE_* VASHARE* VAP* CAPPRICE* LP* *LABQLTY* profit_share* pop* PAT_* PATSH* PATSTOCK_* PATSTOCKSH_* RDSH* SCI_SH
    export excel using "$output/dataset/data_moments_markups_`v'.xlsx", replace firstrow(variables) sheet("DataMoments")
    restore
}

*=======================================================================================
*=======================================================================================
