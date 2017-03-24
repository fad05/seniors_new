cd "/media/alex/Storage/Documents/Projects/Ageing and Education/fortran_code/seniors_new/data_output/stata"
import delimited "../app_age_cd.txt", clear 
*change all (-1) values (values for dead) to missing
mvdecode v*, mv(-1)
rename v2 ret_age
sort ret_age
gen a = 1
by ret_age: egen numret = sum(a)
collapse (mean) numret, by (ret_age)
line numret ret_age,sort
rename numret numret_eta2_8
merge 1:1 ret_age using retirement.dta
drop _merge
save retirement.dta, replace
exit, clear STATA
