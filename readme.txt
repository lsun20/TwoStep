
Package name:  twostepweakiv

Title:  Implementing valid two-step identification-robust confidence sets for linear instrumental-variables models

Author 1 name:  Liyang Sun
Author 1 from:  Liyang Sun, Department of Economics, MIT, Cambridge, US
Author 1 email: lsun20@mit.edu


Help keywords:  twostepweakiv

File list: twostepweakiv.ado twostepweakiv.sthlp 


Notes: 
To replicate the simulation results in Section 2, run the following files in ./simulation_size_distortions/
weakiv_simulation_generate.do
weakiv_p1_rejection.do (outputs monte_carlo_p1_het_wald.dta)
rejection_curves.do
weakiv_simulation_generate.do generates random draws,weakiv_p1_rejection.do calculates test statistics, and rejection_curves.do plots the coverage curves using outputs from weakiv_p1_rejection.do)

To replicate the simulation results in Section 5, run the following files in
./simulation_refined_projection/
Twostep_simulation_generate.do
projection_2sls_power_strong.do (outputs monte_carlo_2sls_power_strong_het.xlsx)
projection_2sls_power_weak.do (outputs monte_carlo_2sls_power_weak_het.xlsx)
power_curves.do
Twostep_simulation_generate.do generates random draws, projection_xx_power_xx.do calculates power, and power_curves.do plots the power curves using outputs from projection_xx_power_xx.do)
