cd "/media/alex/Storage/Documents/Projects/Ageing and Education/fortran_code/seniors_new/data_output/stata"
use lfp_counterfact.dta, replace
drop if age>75
gsort -age
gen ones = 1
gen cumsum = sum(ones)
gen cumlfp_base = sum(lfp_base)/cumsum
gen cumlfp_all2 = sum(lfp_all2)/cumsum
drop if age>67
twoway (line cumlfp_base age)(line cumlfp_all2 age), ytitle(Cumulative LFP up from age) xtitle(Age) legend(label(1 "Baseline") label(2 "Counterfactual"))

