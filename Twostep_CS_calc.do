clear all
set rmsg on, permanently
set more off, permanently
discard
/*
We simulate a linear IV model of sample 10,000 observations with
Y_t = X_t\beta' + \epsilon_t
X_t = Z_t\pi' + u_t
for Y_t a Tx1 vector, X_t a Tx2 vector, and Z_t a Tx5 matrix.
We set the true parameter value \beta_0=(0,0), \epsilon_t and u_t are drawn from a 
bivariate normal distribution with mean 0, variance 1, and correlation \rho = 0.5.
We take Z_t~N(0,I) and \pi = (0,0,0,0,0)
*/
set seed 912
set obs 10000

mata: 
Z = rnormal(10000,5,0,1)
Zvars = st_addvar("float", ("Z1", "Z2", "Z3", "Z4", "Z5"))
st_store(., Zvars ,Z)
end

matrix M = 0,0,0
matrix V = (1,0.5,0.5 \ 0.5,1,0 \ 0.5,0,1)
drawnorm epsilon u1 u2, n(10000) cov(V) means(M)

generate X1_weak_hom = Z1*0 + Z2*0 + Z3*0 + Z4*0 + Z5*0 + u1
generate X2_weak_hom = Z1*0 + Z2*0 + Z3*0 + Z4*0 + Z5*0 + u2
generate X1_weak_het = Z1*0 + Z2*0 + Z3*0 + Z4*0 + Z5*0 + u1*exp(0.5*Z1)
generate X2_weak_het = Z1*0 + Z2*0 + Z3*0 + Z4*0 + Z5*0 + u2*exp(0.5*Z1)

generate X1_strong_hom = Z1*1 + Z2*4+ Z3*9 + Z4*16 + Z5*25 + u1
generate X2_strong_hom = Z1*1 + Z2*4+ Z3*9 + Z4*16 + Z5*25 + u2
generate X1_strong_het = Z1*1 + Z2*4+ Z3*9 + Z4*16 + Z5*25 + u1*exp(0.5*Z1)
generate X2_strong_het = Z1*1 + Z2*4+ Z3*9 + Z4*16 + Z5*25 + u2*exp(0.5*Z1)

generate Y1_weak_hom = X1_weak_hom*0 + epsilon
generate Y1_weak_het = X1_weak_het*0 + epsilon*exp(0.5*Z1)
generate Y1_strong_hom = X1_strong_hom*0 + epsilon
generate Y1_strong_het = X1_strong_het*0 + epsilon*exp(0.5*Z1)

generate Y2_weak_hom = X1_weak_hom*0 + X2_weak_hom*0 + epsilon
generate Y2_weak_het = X1_weak_het*0 + X2_weak_het*0 + epsilon*exp(0.5*Z1)
generate Y2_strong_hom = X1_strong_hom*0 + X2_strong_hom*0 + epsilon
generate Y2_strong_het = X1_strong_het*0 + X2_strong_het*0 + epsilon*exp(0.5*Z1)

*save Hom_IV_5_instruments.dta, replace
*outsheet using Hom_IV_5_instruments.dat, noname replace 
*save IV_5_instruments.dta, replace
outsheet using IV_5_instruments.dat, noname replace 

log using stata.log, text replace
*/

*use Hom_IV_5_instruments.dta, clear
// Test cases for k=p=1: only AR should be performed
weakestiv ivreg2 Y1_weak_hom (X1_weak_hom  = Z1),noconstant robust usegrid

// Test illegal cases:
weakestiv ivreg2 Y2_weak_hom (X1_weak_hom X2_weak_hom = Z1 Z2 Z3 Z4 Z5),noconstant robust
weakestiv ivreg2 Y2_weak_hom (X1_weak_hom X2_weak_hom = Z1 Z2 Z3 Z4 Z5),noconstant robust usegrid 
weakestiv ivreg2 Y2_weak_hom (X1_weak_hom X2_weak_hom = Z1 Z2 Z3 Z4 Z5),noconstant robust gridmin(-5 -5) gridmax(5 5) gridpoints(11 11)

// Test p=1, k=5
weakestiv ivreg2 Y1_weak_hom (X1_weak_hom  = Z1 Z2 Z3 Z4 Z5),noconstant robust gridmin(-5) gridmax(5) gridpoints(111)
weakestiv ivreg2 Y1_weak_het (X1_weak_het  = Z1 Z2 Z3 Z4 Z5),noconstant robust gridmin(-5) gridmax(5) gridpoints(111)
weakestiv ivreg2 Y1_strong_hom (X1_strong_hom  = Z1 Z2 Z3 Z4 Z5),noconstant robust gridmin(-5) gridmax(5) gridpoints(111)
weakestiv ivreg2 Y1_strong_het (X1_strong_het  = Z1 Z2 Z3 Z4 Z5),noconstant robust gridmin(-5) gridmax(5) gridpoints(111)
// default citestlist is k lc_2sls ar

// Test p=1, k=5
// only projection if no citestlsit is listed
weakestiv ivreg2 Y2_weak_hom (X1_weak_hom X2_weak_hom = Z1 Z2 Z3 Z4 Z5),noconstant robust project(X1_weak_hom) 
weakestiv ivreg2 Y2_strong_hom (X1_strong_hom X2_strong_hom = Z1 Z2 Z3 Z4 Z5),noconstant robust project(X1_strong_hom) 
weakestiv ivreg2 Y2_weak_het (X1_weak_het X2_weak_het = Z1 Z2 Z3 Z4 Z5),noconstant robust project(X1_weak_het) 
weakestiv ivreg2 Y2_strong_het (X1_strong_het X2_strong_het = Z1 Z2 Z3 Z4 Z5),noconstant robust project(X1_strong_het) 

// centers around Wald point estimate and width twice the Wald confidence interval 25x25

weakestiv ivreg2 Y (X1 X2 = Z1 Z2 Z3 Z4 Z5),noconstant robust gridmin(-5 -5) gridmax(5 5) gridpoints(11 11) project(X1)
weakestiv ivreg2 Y (X1 X2 = Z1 Z2 Z3 Z4 Z5),noconstant robust gridmin(-5 -5) gridmax(5 5) gridpoints(11 11) clrsims(0) project(X1) citestlist(AR LC_2sls)
// both full tests and projection

log close

/*
weakiv ivreg2 Y (X1 X2 = Z1 Z2 Z3 Z4 Z5),project(_all) gridmin(-50 -50) gridmax(50 50) gridpoints(10 10)
s

// System settings
// Somewhere I need to set a global seed
local state = c(seed)
set seed `state'
/* PACKAGES NEEDED
ssc install moremata (for the quantile function)
*/

*! TwoStep CS

*Calculate value of "a" (the linear combination weight) corresponding to
*gamma_min and different hypotheses

/* Pass variables into Mata environment - later!
local alpha	= 0.05
local crit_sims = 10
local nendog	= 2 // number of endogenous variables
local k 	= 5 // number of moments = number of excluded exogenous variables
local gamma_min = 0.05 // will introduce the vector later
*/

/*We starts with F = I(nendog) and project the confidence interval for now. 
We will calculate the marginal confidence interval later.
So p is just nendo. For simpilicity, start with one endogenous variab
*/

mata:
mata clear

alpha	= 0.05
gamma_min = 0.05
crit_sims = 10^6
nendog	= 2 // number of endogenous variables
k 	= 4 // number of moments = number of excluded exogenous variables
p = nendog

K_size=rchi2(crit_sims,1,p)
J_size=rchi2(crit_sims,1,k-p)
target_quant = invchi2(p,1-alpha)

// Definet he crit_obj function
real scalar crit_obj(real scalar a, 
			real matrix K_size, 
			real matrix J_size,
			real scalar target_quant,
			real scalar alpha,
			real scalar gamma_min
			)
{

	K_J = (1+a)*K_size+a*J_size

	crit_quant=mm_quantile(K_J, 1, 1-alpha-gamma_min)
	
	crit_obj=abs(target_quant-crit_quant)

	return(crit_obj)
}

// Minimize the crit_obj function
void eval_crit_obj(todo, a, K_size, J_size, target_quant, alpha, gamma_min,
			d, g, H)
{
	d = crit_obj(a, K_size, J_size, target_quant, alpha, gamma_min)
}

S = optimize_init() 
optimize_init_evaluator(S, &eval_crit_obj()) 
optimize_init_which(S, "min")
optimize_init_params(S, 0)
// We can only pass 5 parameters into the objective function
optimize_init_argument(S,1,K_size)
optimize_init_argument(S,2,J_size)
optimize_init_argument(S,3,target_quant)
optimize_init_argument(S,4,alpha)
optimize_init_argument(S,5,gamma_min)
// Specify the objective function has no derivative
optimize_init_evaluatortype(S, "d0")
// Use the Nelder-Mead method as Matlab's fminsearch
optimize_init_technique(S, "nm")

//default ptol (tolerance) is 1e-6 but we set it to matlab's default
optimize_init_nmsimplexdeltas(S, 1e-3)
optimize_init_conv_ptol(S, 1e-4)

a_min = optimize(S)
a_min

end
*/
use Hom_IV_5_instruments.dta

capture program drop compute_a_min
program compute_a_min, rclass
/*	syntax [,
		k_w_level(real 0)
		nexexog(integer 0)		///
		nendog(integer 0)		///
		nsendog(integer 0)		/// also boolean for strongly-IDed vars
		]
*/	
local nendog = 2
local nsendog = 0
local nexexog = 5
	tempname a_min k_df j_df
	local k_df		= `nendog' - `nsendog'
	local a = 1
	local g = 1
	

// The first column is p and the second column is k, the third column returns the a_min

	!awk -F"," '($1=="`a'")&&($2=="`g'")&&($3 =="`k_df'")&& ($4 == "`nexexog'") { print $5,$6 > "a.txt"}' "a_min.csv"
  	!awk -F"," '($1=="`a'")&&($2=="`g'")&&($3 =="1")&& ($4 == "`nexexog'") { print $5,$6 >> "a.txt"}' "a_min.csv"

 file open f using "a_min.txt", read
 file read f line
 while `=word("`line'",1)'!=1|`=word("`line'",2)'!=1|`=word("`line'",3)'!=2|`=word("`line'",4)'!=2 {
	file read f line	
 }
tokenize `line', parse(" ")
display "a_min is `5'"
display "lc_2sls_crit is `6'"
 
 display "`=word("`line'",1)'"
  display "`=word("`line'",2)'"
   display "`=word("`line'",3)'"
    display "`=word("`line'",4)'"
     display "`=word("`line'",5)'"
      display "`=word("`line'",6)'"
 tempname fh
 file open `fh' using "a.txt", read
 file read `fh' line
 local linenum = 0
 while r(eof)==0 {
    local linenum = `linenum' + 1
    scalar count = `linenum'
    local o`linenum' = `"`line'"'
    local no = `linenum'
    file read `fh' line
   }
  file close `fh'
  !rm a.tmp

tokenize `o1', parse(" ")
  return local a_min = `1'
  return loca lc_2sls_crit = `2'
  tokenize `o2', parse(" ")	
	return local a_min_p = `1'
	return local lc_2sls_crit_p = `1'
end
compute_a_min
return list
