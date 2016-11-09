cd "/media/alex/Storage/Documents/Projects/Ageing and Education/fortran_code/seniors_new/data_output/stata"
local name "app_age"
import delimited "../`name'.txt", clear 
*change all (-1) values (values for dead) to missing
mvdecode v*, mv(-1)
rename v2 ret_age
sort ret_age
gen a = 1
by ret_age: egen numret = sum(a)
collapse (mean) numret, by (ret_age)
line numret ret_age,sort
