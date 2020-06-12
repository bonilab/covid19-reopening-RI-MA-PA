## Starting Values for MCMC Chains

This file contains manually fit simulations that can be used as starting points for MCMC inference.  Please note the day of each manual fit as some older ones may be out of date if the code has been updated.


### June 9 - Rhode Island Fit

The following ODE call yields a good fit to the RI data:

```
        ./odesim none   -tf 160                 -loc RI                 -introday 55 
                        -prob-icu-vent 0.97     -mean-time-vent 6.0     -death-prob-home-80 0.22 
                        -time-symp-to-hosp 5.0  -dev-len-hospstay 0.65  -dev-icu-frac 0.6 
                        -hosp-frac-10 0.021     -hosp-frac-20 0.017     -hosp-frac-30 0.031 
                        -hosp-frac-40 0.05      -hosp-frac-50 0.10      -hosp-frac-60 0.15 
                        -hosp-frac-70 0.15      -hosp-frac-80 0.15 
                        -beta 2.3 1.70 0.41 0.195 0.195 0.19 0.17 0.14 0.12 0.07 0.06 0.06 0.06 0.06 0.06
```
