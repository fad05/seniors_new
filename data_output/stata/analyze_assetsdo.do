clear all
set more off
cd "/media/alex/Storage/Documents/Projects/Ageing and Education/fortran_code/seniors_new/data_output/stata"
local name "life_assets_cd"
import delimited "../`name'.txt", clear 
*change all (-1) values (values for dead) to missing
mvdecode v*, mv(-1)
*transpose the data
xpose, clear
drop if inlist(_n,1)
gen age = _n+49
egen assets = rowmean(v*)
egen sd_assets = rowsd(v*)
twoway (line assets age, sort), name(`name') title(`name')
drop v*
merge 1:1 age using assets_data.dta
keep if _merge == 3
drop _merge

twoway (line assets age)(line assets_data age), name(assets)

twoway (line sd_assets age)(line sd_assets_data age), name(sdassets)
