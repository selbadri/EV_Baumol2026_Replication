
*=======================================================================================
* An Anatomy of the Baumol Cost Disease -- Elbadri & Vizcaino (2026)
*---------------------------------------------------------------------------------------
* File:    11_scientists.do
* Purpose: R&D worker shares by sector and year, and the RI/NRI classification of service
*          industries (Appendices B and E)
*---------------------------------------------------------------------------------------
* Inputs:  input/usa_00026.dta                                    -- IPUMS USA extract, 1950-2023 (not included; see README)
*          input/All_STEM_Occupations.xlsx                        -- O*NET STEM occupation list
*          input/occ_occsoc_crosswalk_2000_onward.xlsx            -- IPUMS OCCSOC crosswalk, 2000 onward
*          input/manual/occupations_manual_assignment.xlsx        -- Manual assignments (sheets: STEM -- R&D occupations, occupations, ind1990_crosswalk, indnaics_crosswalk)
*          input/manual/scientists_industry_crosswalk.xlsx        -- Manual INDNAICS → industry group crosswalk
*          input/manual/NSCG/nscg_rd_rates_occ1990.csv            -- OCC1990 → NSCG occupation group
*          input/NSCG_2023/college_grads_2023.zip                 -- NSCG 2023 public-use file (NCSES)
*
* Outputs: intermediate/11_scientists.dta                         -- R&D scientist shares by sector and year
*          output/dataset/RD_intensity.xlsx                       -- R&D intensity and classification by industry
*=======================================================================================

set varabbrev off

*=======================================================================================
* STEP 1: Identify R&D Occupations from O*NET STEM list 
*=======================================================================================

* Pre-save inclusion list (manually coded Yes/No, sheet "STEM -- R&D occupations" in occupations_manual_assignment.xlsx)
import excel "$raw_manual/occupations_manual_assignment.xlsx", sheet("STEM -- R&D occupations") firstrow clear
keep if Included == "Yes"
keep Code
save "$interm/11_scientists/occ_codes_included.dta", replace

import excel "$raw/All_STEM_Occupations.xlsx", sheet("Browse STEM Occupations") cellrange(A4:C291) firstrow clear

* Keep R&D and design practitioners only
keep if OccupationTypes == "Research, Development, Design, and Practitioners" | ///
		OccupationTypes == "Research, Development, Design, and Practitioners; Technologists and Technicians"

* Further restrict to manually-approved occupations (see occupations_manual_assignment.xlsx, sheet "STEM -- R&D occupations")
merge m:1 Code using "$interm/11_scientists/occ_codes_included.dta", keep(match) nogen

* Clean SOC codes: strip decimal places and dashes
replace Code = substr(Code, 1, strpos(Code, ".") - 1) if strpos(Code, ".")
gen soc_clean = subinstr(Code, "-", "", .)
drop Code OccupationTypes
rename Occupation ONET_Title

save "$interm/11_scientists/rd_occupations.dta", replace

*=======================================================================================
* STEP 2: Clean IPUMS (O*NET → OCCSOC) crosswalk
*=======================================================================================

import excel "$raw/occ_occsoc_crosswalk_2000_onward.xlsx", ///
	sheet("Codes without Descriptions") firstrow clear
keep OnwardACSPRCS Occupationtitle
drop if OnwardACSPRCS == "" | missing(Occupationtitle)
rename OnwardACSPRCS soc_clean
rename Occupationtitle IPUMS_Title

save "$interm/11_scientists/occsoc_occ_crosswalk.dta", replace

*=======================================================================================
* STEP 3: Map R&D occupations to OCCSOC via cascade 6→5→4-digit merge + manual
*
* Strategy: attempt an exact 6-digit SOC match first. For remaining unmatched entries,
* progressively coarsen the match key: 5-digit (crosswalk entries with one trailing X),
* then 4-digit (entries with two trailing X's). Remaining unmatched are assigned manually.
*=======================================================================================

* Crosswalk subset: 5-digit wildcard entries, either one trailing X (e.g. "17207X") or a
* SOC broad-occupation code ending in 0 (e.g. "191010"), which IPUMS uses for aggregated
* occupations. Without the second case, detailed codes such as 19-1011 fall through to the
* 4-digit "1910XX" (other life scientists) and are misassigned.
use "$interm/11_scientists/occsoc_occ_crosswalk.dta", clear
keep if (substr(soc_clean, -1, 1) == "X" & ///
		 strpos(substr(soc_clean, 1, length(soc_clean) - 1), "X") == 0) | ///
		(substr(soc_clean, -1, 1) == "0" & strpos(soc_clean, "X") == 0)
gen soc_prefix = substr(soc_clean, 1, 5)
rename soc_clean _5d_soc
save "$interm/11_scientists/crosswalk_5digit.dta", replace

* Crosswalk subset: 4-digit wildcard entries (exactly two trailing X's, e.g. "1720XX")
use "$interm/11_scientists/occsoc_occ_crosswalk.dta", clear
keep if substr(soc_clean, -2, 2) == "XX" & ///
		strpos(substr(soc_clean, 1, length(soc_clean) - 2), "X") == 0
gen soc_prefix = substr(soc_clean, 1, 4)
rename soc_clean _4d_soc
save "$interm/11_scientists/crosswalk_4digit.dta", replace

* ---- 6-digit exact match ----
use "$interm/11_scientists/rd_occupations.dta", clear
merge m:1 soc_clean using "$interm/11_scientists/occsoc_occ_crosswalk.dta"
gen level_merge = 6 if _merge == 3
preserve
	keep if _merge == 3
	save "$interm/11_scientists/matched.dta", replace
restore
keep if _merge == 1
drop _merge IPUMS_Title level_merge

* ---- 5-digit match on 6-digit unmatched ----
gen soc_prefix = substr(soc_clean, 1, 5)
merge m:1 soc_prefix using "$interm/11_scientists/crosswalk_5digit.dta"
drop soc_prefix
preserve
	keep if _merge == 3
	gen level_merge = 5
	append using "$interm/11_scientists/matched.dta"
	save "$interm/11_scientists/matched.dta", replace
restore
keep if _merge == 1
drop _merge IPUMS_Title _5d_soc

* ---- 4-digit match on remaining unmatched ----
gen soc_prefix = substr(soc_clean, 1, 4)
merge m:1 soc_prefix using "$interm/11_scientists/crosswalk_4digit.dta"
drop soc_prefix
drop if _merge == 2   // discard crosswalk-only rows
preserve
	keep if _merge == 3
	gen level_merge = 4
	append using "$interm/11_scientists/matched.dta"
	save "$interm/11_scientists/matched.dta", replace
restore
keep if _merge == 1
drop _merge IPUMS_Title _4d_soc
save "$interm/11_scientists/unmatched.dta", replace

* ---- Manual match for remaining unmatched ----
import excel "$raw_manual/occupations_manual_assignment.xlsx", sheet("occupations") cellrange(A1:D56) firstrow clear
drop if missing(soc_clean)
tostring soc_clean, replace
save "$interm/11_scientists/manual_assignment.dta", replace

use "$interm/11_scientists/unmatched.dta", clear
merge 1:1 ONET_Title soc_clean using "$interm/11_scientists/manual_assignment.dta"
drop if _merge == 2
rename Crosswalk_Title IPUMS_Title
drop _merge
rename Code Manual

* ---- Combine all cohorts and build final OCCSOC code ----
append using "$interm/11_scientists/matched.dta"
gen occsoc = ""
replace occsoc = soc_clean if level_merge == 6
replace occsoc = _5d_soc   if level_merge == 5
replace occsoc = _4d_soc   if level_merge == 4
replace occsoc = Manual    if missing(level_merge)
keep ONET_Title soc_clean IPUMS_Title occsoc

* Manual corrections for recently split/renamed SOC codes
replace occsoc      = "172070" if inlist(soc_clean, "172071", "172072")
replace IPUMS_Title = "Electrical and Electronics Engineers" if occsoc == "172070"
replace occsoc      = "172110" if soc_clean == "172112"
replace IPUMS_Title = "Industrial Engineers, including Health and Safety" if occsoc == "172110"
replace occsoc      = "1721YY" if inlist(soc_clean, "172111", "172161", "172199")
replace IPUMS_Title = "Other Engineers" if occsoc == "1721YY"
replace occsoc      = "191020" if soc_clean == "191029"
replace IPUMS_Title = "Biological Scientists" if occsoc == "191020"
replace occsoc      = "299000" if soc_clean == "299091"
replace IPUMS_Title = "Other Healthcare Practitioners and Technical Occupations" if occsoc == "299000"

* Collapse to unique OCCSOC codes (multiple O*NET codes can map to the same OCCSOC)
sort occsoc
gen id = 1
collapse (sum) id, by(occsoc)
drop id
save "$interm/11_scientists/rd_occupations_occsoc.dta", replace

*=======================================================================================
* STEP 4: Map OCCSOC R&D occupations into OCC1990 codes and apply to all IPUMS years
*=======================================================================================

* Sample: employed workers, women and men, aged 21-65 in the private economy (public
* administration, ind1990 900 and above, excluded), census years from 1950. The model's
* researcher share is a share of all employment; the inherited men-only sample was dropped
* on 8 October 2026 (it overstated the share and understated NRI-S, where many women
* researchers work in health care and education).
local ipums_vars "year perwt sex age empstat occ1990 occsoc ind1990 indnaics"
local ipums_if   "year >= 1950 & empstat == 1 & age >= 21 & age <= 65 & ind1990 < 900"

* OCCSOC -> OCC1990 using 2023 ACS observations, where both codings exist. An occ1990 code
* can pool R&D and non-R&D occupations (e.g. economists with market research analysts), so
* each code carries its R&D share of 2023 employment, rd_share in [0,1], rather than a 0/1
* flag set when any member occupation is R&D.
use `ipums_vars' using "$raw/usa_00026.dta" if `ipums_if' & year == 2023, clear
drop if occ1990 == 999 | occ1990 == . | occ1990 >= 905
keep occ1990 occsoc perwt
merge m:1 occsoc using "$interm/11_scientists/rd_occupations_occsoc.dta", keep(master match)
gen byte rd_occ = (_merge == 3)
drop _merge
collapse (mean) rd_share = rd_occ [pw = perwt], by(occ1990)
save "$interm/11_scientists/stem_occ1990.dta", replace

* Apply R&D shares to all IPUMS years
use `ipums_vars' using "$raw/usa_00026.dta" if `ipums_if', clear
drop if occ1990 == 999 | occ1990 == . | occ1990 >= 905
drop if ind1990 == 0 | ind1990 == . | ind1990 >= 940

merge m:1 occ1990 using "$interm/11_scientists/stem_occ1990.dta", keep(master match)

* occ1990 codes absent from the 2023 ACS but R&D occupations in earlier years
replace rd_share = 1 if _merge == 1 & inlist(occ1990, 67, 168) // statisticians, sociologists
replace rd_share = 0 if missing(rd_share)
drop _merge

* Occupation titles overstate research: most engineers and computer workers design, build or
* maintain. Each code's R&D share is multiplied (end of step 5) by the share of college
* graduates in its NSCG occupation group and sector whose primary work activity is R&D (step
* 4b). Here each code only receives its NSCG group; urban planners, architects and foresters
* have none and carry zero.
preserve
	import delimited "$raw_manual/NSCG/nscg_rd_rates_occ1990.csv", clear varnames(1) bindquote(strict) stringcols(4)
	keep occ1990 nscg_group
	tempfile nscg
	save `nscg'
restore
merge m:1 occ1990 using `nscg', keep(master match) gen(_mnscg)
assert _mnscg == 3 if rd_share > 0          // every R&D code must be listed in the mapping
replace rd_share = 0 if rd_share > 0 & missing(nscg_group)
drop _mnscg

* Only these are used below; the saved microdata files are large.
keep year perwt occ1990 ind1990 indnaics rd_share nscg_group
compress

save "$interm/11_scientists/ipums_usa_stemocc.dta", replace

*=======================================================================================
* STEP 4b: R&D activity rates by NSCG occupation group and sector, NSCG 2023 microdata
*=======================================================================================

* R&D is a primary work activity (most hours in the principal job, WAPRI) of basic research,
* applied research or development; design is excluded. With this definition the microdata
* reproduce the published NSCG Table 3-1 rates within 0.009. Rates are computed on employers
* outside public administration and the military, by sector of the employer's Census industry
* code, because R&D intensity within an occupation differs across sectors (programmers and
* engineers in NRI do more R&D). Thin cells are shrunk toward the occupation group's rate
* across these employers with a weight of NSCG_POOL respondents. Engineers share one pooled
* group (ENG) for industrial engineers.
local NSCG_POOL = 30
preserve
	* The public-use file is read from its zip, extracted once into the intermediate folder.
	local here "`c(pwd)'"
	capture mkdir "$interm/11_scientists/nscg"
	quietly cd "$interm/11_scientists/nscg"
	capture confirm file "pcg23Public/epcg23.sas7bdat"
	if _rc unzipfile "$raw/NSCG_2023/college_grads_2023.zip", replace
	quietly cd "`here'"
	import sas N3OCPRNG WAPRI WTSURVY INDCODE using ///
		"$interm/11_scientists/nscg/pcg23Public/epcg23.sas7bdat", clear case(upper)
	foreach v in N3OCPRNG WAPRI INDCODE {
		capture confirm string variable `v'
		if _rc tostring `v', replace format(%02.0f)
		replace `v' = strtrim(`v')
	}
	drop if WAPRI == "L" | WAPRI == "" | WAPRI == "."
	gen byte rd = inlist(WAPRI, "02", "03", "04")
	gen ind = real(INDCODE)
	gen str5 S3 = ""
	replace S3 = "GOODS" if ind <= 3990
	replace S3 = "RI" if inrange(ind, 6470, 6992) | inrange(ind, 7270, 7490) | ind == 7570
	replace S3 = "NRI" if S3 == "" & ind < 9370
	drop if S3 == ""                              // public administration and military
	rename N3OCPRNG nscg_group
	expand 2 if inlist(nscg_group, "51", "52", "53", "54", "55", "56", "57"), gen(dup)
	replace nscg_group = "ENG" if dup
	gen double w_rd = WTSURVY * rd
	gen n = 1
	tempfile resp grp
	save `resp'
	collapse (sum) w_rd WTSURVY, by(nscg_group)
	gen double rate_grp = w_rd / WTSURVY
	keep nscg_group rate_grp
	save `grp'
	use `resp', clear
	collapse (sum) w_rd WTSURVY n, by(nscg_group S3)
	fillin nscg_group S3
	replace n = 0 if _fillin
	merge m:1 nscg_group using `grp', nogen assert(match)
	gen double nscg_rate = cond(n > 0, w_rd / WTSURVY, rate_grp)
	replace nscg_rate = (n * nscg_rate + `NSCG_POOL' * rate_grp) / (n + `NSCG_POOL')
	keep nscg_group S3 nscg_rate n rate_grp
	rename n nscg_n
	save "$interm/11_scientists/nscg_rates_sector.dta", replace
	list if inlist(nscg_group, "11", "ENG", "22"), noobs sep(3)
restore

*=======================================================================================
* STEP 5: Map IPUMS observations to BEA industries via IND1990 → INDNAICS crosswalk
*=======================================================================================

* Pre-save manual INDNAICS → industry crosswalk
import excel "$raw_manual/scientists_industry_crosswalk.xlsx", sheet("crosswalk") firstrow clear
* Utilities (NAICS 22) and retail bakeries (NAICS 311811) are Goods in the paper and in the
* value-added concordance; the manual crosswalk lists them as services.
replace Sector = "GOODS" if substr(indnaics, 1, 2) == "22" | indnaics == "311811"
save "$interm/11_scientists/scientists_industry_crosswalk.dta", replace

* Build a data-driven IND1990 → INDNAICS crosswalk (one INDNAICS per IND1990): among the
* INDNAICS codes listed in the manual industry crosswalk, the one holding the largest weighted
* employment among observations coded under both. Taking the alphabetically first code sent,
* for example, business services n.e.c. to Information. Every IND1990 code has at least one
* listed INDNAICS candidate; unlisted codes would leave workers without an industry.
use "$interm/11_scientists/ipums_usa_stemocc.dta", clear
keep if !missing(ind1990) & !missing(indnaics)
collapse (sum) perwt, by(ind1990 indnaics)
merge m:1 indnaics using "$interm/11_scientists/scientists_industry_crosswalk.dta", ///
	keep(match) keepusing(indnaics) nogen
gsort ind1990 -perwt indnaics
by ind1990: keep if _n == 1
keep ind1990 indnaics
save "$interm/11_scientists/ind1990_indnaics_crosswalk.dta", replace

* Merge crosswalks onto IPUMS observations
use "$interm/11_scientists/ipums_usa_stemocc.dta", clear
drop indnaics
merge m:1 ind1990 using "$interm/11_scientists/ind1990_indnaics_crosswalk.dta", ///
	keep(master match) nogen

* Merge to assign industry_group and industry_name to each observation
merge m:1 indnaics using "$interm/11_scientists/scientists_industry_crosswalk.dta", ///
	keep(master match) nogen

* IND1990 codes that exist only before 2000, when IPUMS has no INDNAICS, have no crosswalk
* pair. Previously these workers were dropped, which removed, for example, electrical
* machinery, photographic equipment and office machines from Goods R&D in 1950-1990.
replace industry_name = "Manufacturing" if inlist(ind1990, 150, 220, 241, 290, 321, 350) | ///
	inlist(ind1990, 380, 381)
replace Sector = "GOODS" if industry_name == "Manufacturing" & missing(indnaics)
replace industry_name = "Forestry, fishing, and related activities" if ind1990 == 30
replace Sector = "GOODS" if ind1990 == 30
replace industry_name = "Retail trade" if inlist(ind1990, 590, 592, 602)
replace industry_name = "Information" if ind1990 == 442
replace industry_name = "Transportation and warehousing" if ind1990 == 422
replace industry_name = "Health care and social assistance" if inlist(ind1990, 830, 863)
replace industry_name = "Educational services" if ind1990 == 851
replace industry_name = "Other services, except government" if ind1990 == 790
* IND1990 893 ("miscellaneous professional and related services") maps to performing arts
* (NAICS 711) from 2000 on, where that is its modal industry. In 1980 and 1990 it holds
* professional-service workers (geologists, psychologists, actuaries, engineers), so it
* belongs to Professional, scientific, and technical services in those years.
replace industry_name = "Professional, scientific, and technical services" if ind1990 == 893 & year <= 1990
replace Sector = "SERVICES" if inlist(ind1990, 422, 442, 590, 592, 602, 790, 830, 851) | ///
	ind1990 == 863

* Every worker must reach an industry and a sector, or the shares silently lose them.
assert !missing(industry_name) & industry_name != "" & Sector != ""

* Researchers in private universities, hospitals and nonprofit research institutes are
* counted in Educational services and in Health care and social assistance, both NRI.

* NSCG activity rate by occupation group and the worker's sector (step 4b). The RI list is the
* one fixed in step 6.
gen str5 S3 = cond(Sector == "GOODS", "GOODS", cond(inlist(industry_name, "Information", ///
	"Finance and insurance", "Professional, scientific, and technical services"), "RI", "NRI"))
merge m:1 nscg_group S3 using "$interm/11_scientists/nscg_rates_sector.dta", ///
	keep(master match) keepusing(nscg_rate) nogen
assert !missing(nscg_rate) if rd_share > 0
replace rd_share = rd_share * nscg_rate if rd_share > 0
drop S3 nscg_rate nscg_group

compress
save "$interm/11_scientists/ipums_usa_industries.dta", replace

*=======================================================================================
* STEP 6: Export R&D Intensity Table
*=======================================================================================

* Compute weighted STEM and non-STEM employment, collapse to industry × sector × year
gen stem_occw     = rd_share * perwt
gen non_stem_occw = (1 - rd_share) * perwt
collapse (sum) stem_occw non_stem_occw, by(industry_name Sector year)
drop if industry_name == "" | year == .
* R&D workers per 100 employees.
gen rd_int = 100 * stem_occw / (stem_occw + non_stem_occw)

* Benchmarks (Appendix E): Goods, Manufacturing, the private economy and all services, each
* computed per census year for the group as a whole and averaged over 1970-2020 and 2000-2020.
foreach b in goods mfg private services {
	preserve
		if "`b'" == "goods"    keep if Sector == "GOODS"
		if "`b'" == "mfg"      keep if industry_name == "Manufacturing"
		if "`b'" == "services" keep if Sector != "GOODS"
		collapse (sum) stem_occw non_stem_occw, by(year)
		gen b_int = 100 * stem_occw / (stem_occw + non_stem_occw)
		quietly summarize b_int if inrange(year, 1970, 2020)
		local `b'_7020 = r(mean)
		quietly summarize b_int if inlist(year, 2000, 2010, 2020)
		local `b'_0020 = r(mean)
	restore
	display as result "`b' R&D workers per 100 employees: 1970-2020 " %6.3f ``b'_7020' ", 2000-2020 " %6.3f ``b'_0020'
}

* Reshape to wide (one row per industry), compute average R&D intensity
reshape wide rd_int stem_occw non_stem_occw, i(industry_name Sector) j(year)
order industry_name Sector stem_occw* non_stem_occw* rd_int*
* Averages over 1970-2020 and 2000-2020 (Appendix E); 2023 is reported but not averaged.
egen rd_av = rowmean(rd_int1970 rd_int1980 rd_int1990 rd_int2000 rd_int2010 rd_int2020)
egen rd_av0020 = rowmean(rd_int2000 rd_int2010 rd_int2020)
sort rd_av

* Services are split into RI and NRI: Information, Finance and insurance, and Professional,
* scientific and technical services are RI: the service industries at or above the average
* of all services over 2000-2020 (the criterion of Appendix E since 8 October 2026; with all
* workers Finance, 0.86 per 100, is below the Goods average, 1.25, and approaches Goods by
* 2023, 1.50 against 1.53). Information and Professional services clear every benchmark;
* Finance clears the services average in both windows and in every census year from 2010. Educational services (0.66) and Health care (0.55) are
* below the services average (0.75) over 2000-2020 and over 1970-2020. The assignment is
* fixed to that list; the flags report it.
foreach b in goods mfg private services {
	gen byte Above_`b'_7020 = (rd_av     >= ``b'_7020') if !missing(rd_av)     & Sector != "GOODS"
	gen byte Above_`b'_0020 = (rd_av0020 >= ``b'_0020') if !missing(rd_av0020) & Sector != "GOODS"
}
replace Sector = cond(inlist(industry_name, "Information", "Finance and insurance", ///
	"Professional, scientific, and technical services"), "RI", "NRI") if Sector != "GOODS"

* Export intensity table (services + Manufacturing only)
preserve
	keep if inlist(Sector, "RI", "NRI") | industry_name == "Manufacturing"
	foreach var of varlist rd_int* rd_av rd_av0020 {
		format `var' %9.2f
	}
	rename rd_int*   RD_per100_*
	rename rd_av     RD_per100_Avg_1970_2020
	rename rd_av0020 RD_per100_Avg_2000_2020
	keep  industry_name RD_per100_* Above_* Sector
	rename industry_name Industry
	export excel using "$output/dataset/RD_intensity.xlsx", replace firstrow(variables) sheet("R&D Intensity Classification")
restore

*=======================================================================================
* STEP 7: Compute R&D Scientist Shares by Sector
*=======================================================================================

collapse (sum) stem_occw* non_stem_occw*, by(Sector)
reshape long  stem_occw non_stem_occw, i(Sector) j(Year)
reshape wide  stem_occw non_stem_occw, i(Year)  j(Sector) string

local sect "GOODS RI NRI"
foreach x of local sect {
	gen RDSH_`x' = stem_occw`x' / (stem_occwGOODS + stem_occwRI + stem_occwNRI)
}
gen SCI_SH = (stem_occwGOODS + stem_occwRI + stem_occwNRI) / ///
             (stem_occwGOODS + stem_occwRI + stem_occwNRI + non_stem_occwGOODS + non_stem_occwRI + non_stem_occwNRI)
keep Year RDSH_* SCI_SH

save "$interm/11_scientists.dta", replace

