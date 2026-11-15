*=======================================================================================
* An Anatomy of the Baumol Cost Disease -- Elbadri & Vizcaino (2026)
*---------------------------------------------------------------------------------------
* File:    00_programs.do
* Purpose: Define shared utility programs used across the EV data pipeline
*---------------------------------------------------------------------------------------
* Programs: ev_load              -- import EV 44-industry data and merge concordance
*           ev_tornqvist_weights -- compute within-sector Tornqvist VA weights
*=======================================================================================

capture program drop ev_load
program define ev_load
    * No arguments.
    * Requires globals: $raw (path to input/), $conc (path to intermediate/concordances/)
    * Loads: $raw/EV_production_accounts_1947to2023.xlsx, sheet "Data_44Ind_1947to2023"
    * Merges: $conc/EV_44Ind_concordance.dta on industry_id (m:1, must be exact)
    * Result: 44-industry annual panel in memory with SectorCode attached
    import excel "$raw/EV_production_accounts_1947to2023.xlsx", ///
        sheet("Data_44Ind_1947to2023") firstrow clear
    merge m:1 industry_id using "$conc/EV_44Ind_concordance.dta"
    assert _merge == 3
    drop _merge

    * The model economy is the private economy: government output is not sold, and its
    * researchers are not allocated by the value of patents. Drop Federal (62) and State
    * and local (63), so every series, including the aggregate VA price, is private.
    drop if inlist(industry_id, 62, 63)
end

capture program drop ev_tornqvist_weights
program define ev_tornqvist_weights
    * No arguments.
    * Requires in memory: industry_id (numeric panel id), year, SectorCode (string), VA (nominal value added)
    * Result: adds tornqv_weight to the dataset in memory
    *         tornqv_weight = average of current and lagged VA share within sector-year
    xtset industry_id year
    bysort SectorCode year: egen Total_VA = total(VA)
    gen VAshare = VA / Total_VA
    sort industry_id year
    gen lag_VAshare = L.VAshare
    gen tornqv_weight = (VAshare + lag_VAshare) / 2
end

*=======================================================================================
*=======================================================================================
