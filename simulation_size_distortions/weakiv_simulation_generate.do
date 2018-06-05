// Generate simulation datasets (Z's, epsilon and u's)
// p=1 k=10: one endogenous variable and ten instruments
// simulation design follows the moderate endogeneity under heteroskedasticity
// in the supplementary appendix to 
// Andrews, Isaiah, "Valid Two-Step Identification-Robust Confidence Sets for GMM", Review of Economics and Statistics (2018 May).
clear
global monte_carlo_data "./monte_carlo_p1_data/"

forval j = 1/2500{
	dis "simulation is `j'"
	
	quietly {
	clear
	local s = `j' + 500
	set seed `s'
	set obs 10000
	
	mata: n = 10000
	mata: p = (1/10,1/10,1/10,1/10,1/10,1/10,1/10,1/10,1/10,1/10)
	mata: Z_tilde = rdiscrete(n, 1, p)
	mata: Z = J(n,10,0)
	
	mata: Z[.,1] = Z_tilde :== 1
	mata: Z[.,2] = Z_tilde :== 2
	mata: Z[.,3] = Z_tilde :== 3
	mata: Z[.,4] = Z_tilde :== 4
	mata: Z[.,5] = Z_tilde :== 5
	mata: Z[.,6] = Z_tilde :== 6
	mata: Z[.,7] = Z_tilde :== 7
	mata: Z[.,8] = Z_tilde :== 8
	mata: Z[.,9] = Z_tilde :== 9
	mata: Z[.,10] = Z_tilde :== 10
	//mata: sum(Z) // should be 500

	mata: Zvars = st_addvar("float", ("Z1", "Z2", "Z3", "Z4", "Z5", "Z6","Z7","Z8","Z9","Z10"))
	mata: st_store(., Zvars ,Z)

	matrix M = 0,0
	 
	matrix V1 = (1.102,0.004 \ 0.004,0.004)
	matrix V2 = (1.891,-2.045 \ -2.045,2.789)
	matrix V3 = (3.844,-3.969 \ -3.969,4.264)
	matrix V4 = (3.457,-1.369 \ -1.369,0.779)
	matrix V5 = (0.043,0.129 \ 0.129,0.395)
	matrix V6 = (0.832,-2.383 \ -2.383,7.026)
	matrix V7 = (1.111,-1.135 \ -1.135,1.226)
	matrix V8 = (2.254,-0.39 \ -0.39,0.308)
	matrix V9 = (10.786,-2.859 \ -2.859,1.709)
	matrix V10 = (0.419,-0.934 \ -0.934,6.099)


	gen epsilon = .
	gen u = .
	forval i = 1(1)10 {
		drawnorm epsilon`i' u`i' ,  cov(V`i') means(M) 
		replace epsilon = epsilon`i' if Z`i' == 1
		replace u = u`i' if Z`i' == 1
	}
	drop epsilon1 - u10
	gen s = `j'
	save "${monte_carlo_data}data`j'.dta", replace
	}
	
}

