
import excel "./monte_carlo_2sls_power_weak_het.xlsx", clear firstrow
//import excel "./monte_carlo_2sls_power_strong_het.xlsx", clear firstrow

sum simulation

sum F1
local F1_mean: di %3.2f `=r(mean)'
sum F2
local F2_mean: di %3.2f `=r(mean)' 
collapse (mean) r*, by(theta1) // theta1 is beta in the main text

label var r_lc_2sls "conventional lc_2sls"
label var r_k_2sls "conventional k_2sls"
label var r_lc_2sls_refined  "refined lc_2sls"
label var r_k_2sls_refined "refined k_2sls"

local norm =  (0.25^2 *5)^(1/2) // weak 
//local norm =  (0.5^2 *5)^(1/2) // strong
replace theta1 = theta1 / `norm' // rescale by the norm


line r* theta1 , ylabel(0(0.1)1) xlabel(-3(0.5)3) ytitle("pr(reject {&beta}=0)") xtitle("{&beta}{sub:0}/||{&pi}{sub:2}|| ") ///
	lpattern("l" "_" "-" "_-") note("Mean of First Stage F statistic: {&eta}: `F2_mean'; {&beta}: `F1_mean'") scheme(sj)


//graph export "./power_2sls_strong_het.eps", replace
graph export "./power_2sls_weak_het.eps", replace


