// Generate simulation datasets (Z's, epsilon and u's)

clear
global monte_carlo_data "./monte_carlo_data/"
forval j = 1/500{
	dis "simulation is `j'"
	
	quietly {
	clear
	set seed `j'
	set obs 500
	
	mata: n = 500
	mata: Z = rnormal(n,5,0,1)
	mata: Zvars = st_addvar("float", ("Z1", "Z2", "Z3", "Z4", "Z5"))
	mata: st_store(., Zvars ,Z)

	matrix M = 0,0,0
	matrix V = (1,0.5,0.5 \ 0.5,1,0 \ 0.5,0,1)
	drawnorm epsilon u1 u2, n(500) cov(V) means(M)
	
	gen s = `j'
	save "${monte_carlo_data}data`j'.dta", replace
	}
	
}

