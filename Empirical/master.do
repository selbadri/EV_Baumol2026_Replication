*=======================================================================================
* An Anatomy of the Baumol Cost Disease -- Elbadri & Vizcaino (2026)
*---------------------------------------------------------------------------------------
* File:    master.do
* Purpose: Run the whole Empirical step in order (code/00-12) and write the data moments
*---------------------------------------------------------------------------------------
* Usage:   Run from the package root, with no edits:
*              cd "path/to/Empirical"
*              do master.do
*          To run from elsewhere, set the root first:
*              global root "path/to/Empirical"
*
* Inputs:  input/ (see README), ../Cross_Citations/output/patents/
* Outputs: output/dataset/data_moments.xlsx, output/master_log.log
*=======================================================================================

clear all
set more off
version 18

*---------------------------------------------------------------------------------------
* 1. Directories
*---------------------------------------------------------------------------------------

if `"$root"' == "" global root "`c(pwd)'"

capture confirm file "$root/master.do"
if _rc {
    display as error "Package root not found: $root"
    display as error "cd to the Empirical folder, or set global root, then run master.do."
    exit 601
}

global code        "$root/code"
global raw         "$root/input"
global raw_manual  "$root/input/manual"
global interm      "$root/intermediate"
global conc        "$root/intermediate/concordances"
global output      "$root/output"
global crosscit    "$root/../Cross_Citations/output/patents"

* Create the generated folders, so the pipeline runs from an empty intermediate/ and output/
foreach d in "$interm" "$conc" "$interm/07_markups" "$interm/10_patents" ///
             "$interm/11_scientists" "$output" "$output/dataset" {
    capture mkdir "`d'"
}

*---------------------------------------------------------------------------------------
* 2. Log
*---------------------------------------------------------------------------------------

capture log close
log using "$output/master_log.log", replace

*---------------------------------------------------------------------------------------
* 3. Pipeline
*---------------------------------------------------------------------------------------

display _n "[0] Shared programs"
do "$code/00_programs.do"

display _n "[1] Concordance tables"
do "$code/01_concordance_tables.do"

display _n "[2] Labor shares of value added"
do "$code/02_labor_shares.do"

display _n "[3] Value-added shares"
do "$code/03_value_added_shares.do"

display _n "[4] Value-added prices"
do "$code/04_value_added_prices.do"

display _n "[5] Labor productivity growth"
do "$code/05_labor_productivity_growth.do"

display _n "[6] Labor quality growth"
do "$code/06_labor_quality_growth.do"

display _n "[7a] BEA capital asset panel"
do "$code/07a_capital_assets.do"

display _n "[7b] Internal and external rates of return"
do "$code/07b_rates_of_return.do"

display _n "[7c] Machine-producer markups"
do "$code/07c_markups.do"

display _n "[8] Capital prices"
do "$code/08_capital_prices.do"

display _n "[9] Population growth"
do "$code/09_demographics.do"

display _n "[10] Patent shares"
do "$code/10_patents.do"

display _n "[10b] Patent stocks"
do "$code/10b_patent_stocks.do"

display _n "[11] Scientist shares"
do "$code/11_scientists.do"

display _n "[12] Data moments"
do "$code/12_compilation.do"

*---------------------------------------------------------------------------------------
* 4. Close log
*---------------------------------------------------------------------------------------

log close
display _n "Empirical finished."
