clear all
set more off
cd "/media/alex/Storage/Documents/Projects/Ageing and Education/fortran_code/seniors_new/data_output/stata"
local name ""
import delimited "../life_labor_cd_6`name'.txt", clear 
*change all (-1) values (values for dead) to missing
mvdecode v*, mv(-1)
*transpose the data
xpose, clear
drop if inlist(_n,1)
gen age = _n+49
*Generate age categories
/*
gen agecat = .
replace agecat=52 if age>=50 & age<55
replace agecat=57 if age>=55 & age<60
replace agecat=62 if age>=60 & age<65
replace agecat=67 if age>=65 & age<70
replace agecat=72 if age>=70 & age<75
replace agecat=77 if age>=75 & age<80
replace agecat=82 if age>=80 & age<85
replace agecat=87 if age>=85 & age<=90
*replace agecat=92 if age>=90
*/
gen agecat = .
replace agecat=51 if age>=50 & age<53
replace agecat=54 if age>=53 & age<56
replace agecat=57 if age>=56 & age<59
replace agecat=60 if age>=59 & age<62
replace agecat=63 if age>=62 & age<65
replace agecat=66 if age>=65 & age<68
replace agecat=69 if age>=68 & age<71
replace agecat=72 if age>=71 & age<74
replace agecat=75 if age>=74 & age<77
replace agecat=78 if age>=77 & age<80
replace agecat=81 if age>=80 & age<83
replace agecat=84 if age>=83 & age<86
replace agecat=87 if age>=86 & age<=90

egen hours = rowmean(v*)
egen medhours = rowmedian(v*)
egen sd_hours = rowsd(v*)

/*
*preserve
ren hours hours_`name'
keep age hours_`name'
merge 1:1 age using hours_counterfact.dta
drop _merge
save hours_counterfact.dta, replace
*restore
*/
drop v*
merge 1:1 age using lfp_hours_data_6.dta
*keep if _merge == 3
*drop _merge

twoway (line hours age)(line mean_hours age),  yscale(range(0 2500)) name(hours) ytitle(Annual hours worked) xtitle(Age) legend(label(1 "Model") label(2 "Data"))
quietly graph export "model_fit/hours_model_data.png", width(3000) height(2500) replace

collapse hours hours_data, by(agecat)

twoway (line hours agecat)(line hours_data agecat), yscale(range(0 2500)) name(acecat_hours)


/*
rename mean_ls `name'_ret65
merge 1:1 age using `name'.dta
drop _merge
save `name'.dta, replace
