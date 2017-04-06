clear all
set more off
cd "/media/alex/Storage/Documents/Projects/Ageing and Education/fortran_code/seniors_new/data_output/stata"
*Load data for type 2
use lfp_hours_data_2.dta 
*merge with counterfactuals of model type 6 with inputs of type 2 (noncollege males of two different cohorts)
merge 1:1 age using lfp_counterfact.dta
drop _merge
merge 1:1 age using hours_counterfact.dta
drop _merge
drop if age>67
sort age
*twoway (line lfp_base age)(connected lfp_wage2 age )(line mean_lfp age), ylabel(0(0.2)1,grid) name(lfp_wage2)  yscale(range(0 1))  ytitle(Mean labor force participation) xtitle(Age) legend(label(1 "Baseline") label(2 "Counterfactual") label(3 "Data, cohort 1915-1934"))
*twoway (line hours_base age)(connected hours_wage2 age )(line mean_hours age), ylabel(0(0.2)1,grid) name(hours) yscale(range(1000 2000))  ytitle(Mean annual hours per worker) xtitle(Age) legend(label(1 "Baseline") label(2 "Counterfactual") label(3 "Data, cohort 1915-1934"))
twoway (line lfp_base age)(connected lfp_surv2 age )(line mean_lfp age), ylabel(0(0.2)1,grid) name(lfp__policy)  yscale(range(0 1))  ytitle(Mean labor force participation) xtitle(Age) legend(label(1 "Baseline") label(2 "Counterfactual") label(3 "Data, cohort 1915-1934"))

