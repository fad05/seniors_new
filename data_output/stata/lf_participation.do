clear all
set more off
cd "/media/alex/Storage/Documents/Projects/Ageing and Education/fortran_code/seniors_new/data_output/stata"
import delimited "../life_labor.txt", clear 
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
twoway (line lfp age, sort)
