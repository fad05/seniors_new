clear all
set more off
cd "/media/alex/Storage/Documents/Projects/Ageing and Education/fortran_code/seniors_new/data_output/stata"
*Load data for type 2
use lfp_hours_data_2.dta 
*merge with counterfactuals of model type 6 with inputs of type 2 (noncollege males of two different cohorts)
merge 1:1 age using lfp_counterfact.dta
drop _merge
merge 1:1 age using hours_counterfact.dta
