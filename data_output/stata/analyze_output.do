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
egen hours = rowmean(v*)
egen sd_hours = rowsd(v*)
twoway (line hours age, sort), name(`name') title(`name')
drop v*
merge 1:1 age using labor_data.dta
keep if _merge == 3
drop _merge

twoway (line hours age)(line hours_data age), name(hours)

twoway (line sd_hours age)(line sd_hours_data age), name(sdhours)


/*
rename mean_ls `name'_ret65
merge 1:1 age using `name'.dta
drop _merge
save `name'.dta, replace
