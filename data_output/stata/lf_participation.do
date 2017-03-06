clear all
set more off
cd "/media/alex/Storage/Documents/Projects/Ageing and Education/fortran_code/seniors_new/data_output/stata"
import delimited "../life_labor_cd.txt", clear 
*change all (-1) values (values for dead) to missing
mvdecode v*, mv(-1)
forvalues i = 2/42 {
replace v`i' = 1 if v`i'>0 & v`i' ~= .
}

*transpose the data
xpose, clear
drop if inlist(_n,1)
gen age = _n+49
egen lfp = rowmean(v*)
egen sd_lfp = rowsd(v*)
twoway (line lfp age, sort)

drop v*

merge 1:1 age using labor_data.dta
keep if _merge == 3
drop _merge


twoway (line lfp age)(line lfp_data age), name(lfp)

*twoway (line sd_lfp age)(line sd_lfp_data age), name(sdlfp)
