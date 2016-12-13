clear all
set more off
cd "/media/alex/Storage/Documents/Projects/Ageing and Education/fortran_code/seniors_new/data_output/stata"
local name "life_labor"
import delimited "../`name'.txt", clear 
*change all (-1) values (values for dead) to missing
mvdecode v*, mv(-1)
*transpose the data
xpose, clear
drop if inlist(_n,1)
gen age = _n+49
egen mean_ls = rowmean(v*)
*twoway (line mean_ls age, sort), name(`name') title(`name')
keep age mean_ls
rename mean_ls `name'_eta2_8
merge 1:1 age using `name'.dta
drop _merge
save `name'.dta, replace
