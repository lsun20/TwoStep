
import excel "./monte_carlo_2sls_power_weak_het.xlsx", clear firstrow
//import excel "./monte_carlo_2sls_power_strong_het.xlsx", clear firstrow

sum simulation

sum F1
local F2_mean: di %3.2f `=r(mean)'
sum F2
local F1_mean: di %3.2f `=r(mean)' // I flipped the parameters in write-up so X2 is X1
collapse (mean) r*, by(theta1)

label var r_lc_2sls "conventional lc_2sls"
label var r_k_2sls "conventional k_2sls"
label var r_lc_2sls_refined  "refined lc_2sls"
label var r_k_2sls_refined "refined k_2sls"
replace theta1 = theta1 / (0.25^2 *5)^(1/2) // rescale by the norm


line r* theta1 , ylabel(0(0.1)1) xlabel(-3(0.5)3) ytitle("pr(reject {&beta}=0)") xtitle("{&beta}{sub:0}") ///
	lpattern("l" "_" "-" "_-") note("Mean of First Stage F statistic: {&eta}: `F1_mean'; {&beta}: `F2_mean'") scheme(sj)


//graph export "./power_2sls_strong_het.eps", replace
graph export "./power_2sls_weak_het.eps", replace


