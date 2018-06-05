// Monte Carlo exercise
clear all
set rmsg on, permanently
set more off, permanently
cap log close
global monte_carlo_data "./monte_carlo_p1_data/"


matrix r = J(1,8,.)

mata: R = J(1,8,.) // initiate a mata matrix to store rejection results

forval j =  1/2500{

	use  "${monte_carlo_data}data`j'.dta", clear
	dis "file `j'"
	generate Y = epsilon
	generate X = u
	
	local norm_prev = 0
	foreach norm of numlist 0(0.01)0.06 0.061(0.001)0.079 0.8(0.01)0.11 {
	quietly {

	
		if `norm' > 0 {
			replace X = X + Z1*0.57*(`norm'-`norm_prev') + ///
			Z2*-0.223*(`norm'-`norm_prev') + Z3*-0.478*(`norm'-`norm_prev') + Z4*0.146*(`norm'-`norm_prev') + Z5*0.213*(`norm'-`norm_prev') + Z6*0.08*(`norm'-`norm_prev') + ///
			Z7*-0.167*(`norm'-`norm_prev') + Z8*-0.107*(`norm'-`norm_prev') + Z9*-0.356*(`norm'-`norm_prev') + Z10*-0.396*(`norm'-`norm_prev')
		}
			
			
			matrix r = `j', `norm'
			
			ivreg2 Y (X = Z*), robust first nocons
			
			matrix r = r, e(archi2p), e(cdf), e(rkf)
			
			test X = 0
			matrix r = r, r(p)
			
			
			twostepweakiv 2sls Y ( X = Z*),  citestlist(lc_2sls) nocons gridmin(-1) gridmax(1)
			matrix r = r, e(gamma_hat)
			twostepweakiv 2sls Y ( X = Z*),  citestlist(lc_2sls) gridmin(0) gridmax(0) nocons
			mata:citable = st_matrix("e(citable)")
			mata:table	= select(citable,citable[,1]:==0) 
			mata:r_lc_2sls	= table[.,3] :== 1 // reject indicator
			mata:st_matrix("lc_2sls",r_lc_2sls)
			matrix r = r, lc_2sls
			
			mata: rtable = st_matrix("r")
			mata: R = R \ rtable

			
			local norm_prev = `norm'
			
		}
	}

}


drop _all
mata: newvars = st_addvar("float", ("simulation", "norm", "arp", "first_cd", "first_rk", "waldp","gamma_hat","lc_2sls"))
mata: st_addobs(rows(R) - st_nobs() - 1) // set data matrix to right size
mata: st_store(., newvars, R[|2,1 \ .,.|] )
save .monte_carlo_p1_het_wald.dta, replace
