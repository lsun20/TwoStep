
Package name:  twostepweakiv

Title:  Implementing valid two-step identification-robust confidence sets for linear instrumental-variables models

Author 1 name:  Liyang Sun
Author 1 from:  Liyang Sun, Department of Economics, MIT, Cambridge, US
Author 1 email: lsun20@mit.edu

Author 2 name:  
Author 2 from:  
Author 2 email: 

Author 3 name:  
Author 3 from:  
Author 3 email: 

Author 4 name:  
Author 4 from:  
Author 4 email: 

Author 5 name:  
Author 5 from:  
Author 5 email: 

Help keywords:  twostepweakiv

File list: twostepweakiv.ado twostepweakiv.sthlp 


Notes: 
To replicate the simulation results in Section 2, run the following in the same directory
1. weakiv_simulation_generate.do (outputs simulated datasets to ./monte_carlo_p1_data/)
2. weakiv_p1_rejection.do (outputs monte_carlo_p1_het_wald.dta)
3. rejection_curves.do
weakiv_simulation_generate.do generates random draws,weakiv_p1_rejection.do calculates test statistics, and rejection_curves.do plots the coverage curves using outputs from weakiv_p1_rejection.do.

To replicate the simulation results in Section 5, run the following in the same directory
1. Twostep_simulation_generate.do (outputs simulated datasets to ./monte_carlo_data/)
2. projection_2sls_power_strong.do (outputs monte_carlo_2sls_power_strong_het.xlsx)
3. projection_2sls_power_weak.do (outputs monte_carlo_2sls_power_weak_het.xlsx)
4. power_curves.do
Twostep_simulation_generate.do generates random draws, projection_2sls_power_xx.do calculates power, and power_curves.do plots the power curves using outputs from projection_2sls_power_xx.do (need to (un)comment the appropriate filename and norm local)

