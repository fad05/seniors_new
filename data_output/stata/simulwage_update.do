*This script mimics he procedure I use to obtain wage profiles and deviations from average profiles in the data.
*The difference is, I use model generated data here.
clear all
set more off
cd "/media/alex/Storage/Documents/Projects/Ageing and Education/fortran_code/seniors_new/data_output/stata"
import delimited "../life_wage.txt", clear 
mvdecode v*, mv(-1)
*transpose the data
xpose, clear
drop if inlist(_n,1)
gen age = _n+49
*reshape the data
reshape long v, i(age) j(pid)
rename v hourwage
sort age pid
*log hourly wage
gen loghwage = log(hourwage)
reg loghwage age
*predicted log wage
predict logwagepred
*"residuals" by age
gen zi = loghwage - logwagepred
sort age zi
by age: egen zi_mean = mean(zi)
by age: egen zi_sd = sd(zi)
*construct new fitted log wage mean and sd profiles by age
gen age2 = age^2
gen age3 = age^3
gen age4 = age^4
reg zi_sd age age2 age3 age4
predict zisd_pred
reg zi_mean age
predict zimean_pred
preserve
collapse logwagepred zimean_pred zisd_pred, by(age)
*Store counterfactual smoothed profile
export delimited logwagepred zimean_pred zisd_pred using "/media/alex/Storage/Documents/Projects/Ageing and Education/fortran_code/seniors_new/data_inputs/lwage_smooth_simul.csv", nolabel replace
restore

*Now, update the residuals.
*For each age, split residuals into 5 quintiles
egen zi_bin = xtile(zi), by (age) nq(5) 
*For each age and each zi quintile, find min and max value of zi
sort age zi_bin zi
by age zi_bin: egen bin_min = min(zi)
by age zi_bin: egen bin_max = max(zi)

*Save bin borders as a function of age into a file
sort zi_bin age
preserve
collapse bin_min bin_max, by (zi_bin age)
drop if zi_bin == .
reshape wide bin_min bin_max, i(age) j(zi_bin)
rename (bin_min1 bin_min2 bin_min3 bin_min4 bin_min5 bin_max5) (b1 b2 b3 b4 b5 b6)
*Smooth the borders
keep b? age
forvalues i = 1/6 { 
	reg b`i' age
	predict bpred`i'	
}
export delimited bpred? using "/media/alex/Storage/Documents/Projects/Ageing and Education/fortran_code/seniors_new/data_inputs/wagebin_borders_smooth_simul.csv", delimiter(";") nolabel replace 
restore

*wage shock transitions
sort pid age
gen zibin_next = zi_bin[_n+1]
gen delta_age = age[_n+1] - age[_n]
	*dropping last observation for each person; the information is still contained in healthy_next
by pid: gen numobs = _n
by pid: egen maxnumobs = max(numobs)
drop if numobs == maxnumobs
drop numobs maxnumobs

mlogit zibin_next zi_bin age, base(3)
predict zn1 zn2 zn3 zn4 zn5
sort zi_bin age zibin_next
*twoway (line zn1 age if zi_bin == 2 & zibin_next == 1)(line zn2 age if zi_bin == 2 & zibin_next == 2)(line zn3 age if zi_bin == 2 & zibin_next == 3)(line zn4 age if zi_bin == 2 & zibin_next == 4)(line zn5 age if zi_bin == 2 & zibin_next == 5)
*sort age zi_bin zibin_next
drop if zi_bin == . | zibin_next == .
collapse zn?, by(age zi_bin)

export delimited zn1 zn2 zn3 zn4 zn5 using "/media/alex/Storage/Documents/Projects/Ageing and Education/fortran_code/seniors_new/data_inputs/transmat_wage_simul.csv", delimiter(";") nolabel novarnames replace
exit, clear STATA
