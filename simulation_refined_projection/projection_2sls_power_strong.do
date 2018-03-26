// Monte Carlo exercise
clear all
set rmsg on, permanently
set more off, permanently
cap log close
//global monte_carlo_data "/Users/lsun20/Dropbox (MIT)/2017Summer/TwoStep-local/monte_carlo_data/"
//global monte_carlo_data "E:\Dropbox (MIT)\TwoStep-local\monte_carlo_data\"
global monte_carlo_data "./monte_carlo_data/"

//putexcel set "monte_carlo_2sls_power_strong_het.xlsx", modify
putexcel set "monte_carlo_2sls_power_strong_het.xlsx", replace
matrix r = J(1,8,.)
mat colnames r = theta1 simulation r_lc_2sls r_k_2sls r_lc_2sls_refined r_k_2sls_refined F1 F2
putexcel A1 = matrix(r), colnames

 local norm =  (0.5^2 *5)^(1/2) 
// local norm = (0.03^2 * 5)^(1/2) // weak first stage
local i = 2
forval j = 1/500{
	use  "${monte_carlo_data}data`j'.dta", clear
	dis "file `j'"
	
	gen het = (exp(0.5*Z1))^(1/2) // scale by the heteroskedasticity term
	//gen het = 1
	generate X2_strong_het = Z1*0.1 + Z2*0.2 + Z3*0.3 + Z4*0.4 + Z5*0.5 + u1*het 
	generate X1_strong_het = Z1*0.5 + Z2*0.5 + Z3*0.5 + Z4*0.5 + Z5*0.5 + u2*het
	
	forval t = -25/25 { // use integer increment to avoid Stata's float point problem
		dis "t is `t'"

		quietly {

			local theta = `t'/25*`norm'*3 // over a range of [-3,3]
			local tl = -`norm'*6
			local tu = `norm'*6
			cap drop Y2_strong_het
			generate Y2_strong_het = X1_strong_het*`theta' + X2_strong_het*0 + epsilon*het
			
			twostepweakiv 2sls Y2_strong_het (X1_strong_het X2_strong_het = Z1 Z2 Z3 Z4 Z5) if _n <= 500,noconstant gridmin(0 `tl') gridmax(0 `tu') gridpoints(1 101) citestlist(lc_2sls k_2sls)
		
			mata:table = st_matrix("e(citable)")

		
			mata:rsum	= colsum(table)
			mata:r	= rsum[.,(4,7)]:== 101 // 4th column is lc_r and 7th is k_r - =101 points if every point rejected, then excluded in CI
			mata:st_matrix("r(reject)",r)
			matrix r_project = r(reject)
			
			twostepweakiv 2sls Y2_strong_het (X1_strong_het X2_strong_het = Z1 Z2 Z3 Z4 Z5) if _n <= 500 ,noconstant gridmin(0 `tl') gridmax(0 `tu') gridpoints(1 101) project(X1_strong_het) ptestlist (lc_2sls k_2sls)
			
			mata:citable = st_matrix("e(p1citable)")
			
			mata:table	= select(citable,citable[,1]:==0) 
			//mata: table	= select(citable,(0\1\0))
			mata:r	= table[.,(2,3)] :== 1 // reject indicator
			mata:st_matrix("r(reject)",r)
			matrix r_refined = r(reject)
		

			matrix r = `theta',`j',r_project, r_refined, e(F)'
			
			putexcel A`i' = matrix(r) 
			local i = `i' + 1
		}

	}
}
