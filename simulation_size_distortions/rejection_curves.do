// Plot coverage of 2SLS confidence sets
use "./monte_carlo_p1_het_wald.dta", clear 
sum simulation
tab norm

rename first_rk  first // use the Kleibergen-Paap Wald F statistic
gen r_ar = arp < 0.05
gen r_wald = waldp < 0.05
gen weakF = first < 38.54 // SY critical value

gen r_2sF = .
replace r_2sF = r_ar if weakF == 1
replace r_2sF = r_wald if weakF == 0

rename lc_2sls r_lc2sls
gen weak2s = gamma_hat > 10


gen r_2sR = .
replace r_2sR = r_lc2sls if weak2s == 1
replace r_2sR = r_wald   if weak2s == 0



collapse (mean) first r_ar r_wald weakF r_2sF r_lc2sls weak2s r_2sR  gamma_hat, by(norm)

label var r_ar "AR test size"
label var r_wald "t-test size"
label var r_lc2sls  "LC_2sls test size"

gen a_ar = 1 - r_ar
gen a_wald = 1 - r_wald
gen a_2sF = 1 - r_2sF
gen a_2sR = 1- r_2sR


gen a_lc2sls = 1- r_lc2sls

label var a_ar "AR coverage"
label var a_wald "Wald coverage"
label var a_lc2sls  "LC_2sls coverage"

label var first "Mean Kleibergen-Paap rk Wald F Statistic"

line a_ar a_wald a_2sF first if norm<=0.11, ylabel(0(0.2)1 .85 .95, angle(0.2))   ytitle("pr(CS includes {&beta}=0)")  ///
	lpattern("l" "-." "-" "_-") scheme(sj) 	yline(.85 .95 , lstyle(grid)) ///
	legend(pos(5) ring(0) col(1) lab(3 "coverage after" "pretesting based on Stock and Yogo (2005) cutoff")) 
	
graph export "./coverage_p1_het.eps", replace

line a_lc2sls a_wald a_2sR first if norm<=0.11, ylabel(0(0.2)1 .85 .95, angle(0.2))    ytitle("pr(CS includes {&beta}=0)")   ///
	lpattern("l" "-." "-" "_-") scheme(sj) ///
	legend(pos(5) ring(0) col(1) lab(3 "coverage after" "pretesting based on whether distortion cutoff > 10")) ///
	yline(.85 .95, lstyle(grid))

graph export "./coverage_2sR_p1_het.eps", replace



