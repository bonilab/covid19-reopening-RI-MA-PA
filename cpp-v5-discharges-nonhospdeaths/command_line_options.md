### File describing the command-line usage of the executable 'odesim'

##### First argument is always the output filename

You can call `./odesim myoutput.tdl`

or `./odesim none`

and the 'none' keyword will skip outputting to any file, and everything gets sent to stdout instead.


##### -beta

You can call the function with an arbitrary number of beta parameters (transmission parameters) that allow for time-varying transmission in the model.  For example, you can call

   `./odesim none -beta 1.0 0.9 0.8 0.7`
   
and this will run the model with four different transmission rates, evenly distributed, in the order above, for the period between March 1 and the end of the simulation which is denoted by `tf = time final`.  


> WARNING
>
> For now, you have to set time final on the command line before you assign the betas
>
> e.g. `./odesim none -tf 120 -beta 3 4 5 6`




##### -death-prob-home-60, -death-prob-home-70, -death-prob-home-80

This is the probability of dying at home for individuals in their 60s, 70s, and 80s.

```diff
! Recommended prior distributions 

!  -death-prob-home-60 = [0.0 - 0.10]
!  -death-prob-home-70 = [0.0 - 0.20]
!  -death-prob-home-80 = [0.0 - 0.30]

```


##### -dev-icu-frac

This is the deviation in the parameter that determines the fraction of hospitalized individuals that are admitted to the ICU. This deviation will apply across all age groups.  Default is 1.0.  If you set this higher than 5.0, it will be reset back down to 5.0 

```diff
! Recommended prior distribution [0.5 - 1.5]
```

##### -dev-len-hospstay

This is the deviation in the parameter that determines the length of the typical hospital stay on the medical floor.  Default is 1.0.  If you set this higher than 2.5, it will be reset back down to 2.5 

```diff
! Recommended prior distribution [0.5 - 2.0]
```


##### -dev-ventdeath-mid

This is the deviation in the parameter that determines the fraction of 40-70 individuals that die after being on a ventilator. Default is 0.7.  Do not set this higher than 1.7. 

```diff
! Recommended prior distribution [0.5 - 1.3]
```


##### -diag

This is a flag that prints out some basic diagnostic info as to who was infected, who was hospitalized, and who died.  To be used like this:

    `./odesim none -diag`

##### -hosp-frac-10, -hosp-frac-20, -hosp-frac-30, -hosp-frac-40, -hosp-frac-50, -hosp-frac-60, -hosp-frac-70, -hosp-frac-80
    
This is the fraction of symptomatic individuals that are eventually hospitalized, by age group. `-hosp-frac-10` is the hospitalization probability for individuals in the 10-19 age group. Once this probability is set, it is also applied to the 0-9 age group. `-hosp-frac-20` is for the 20-29 age group, `-hosp-frac-30` is for the 30-39 age group, etc.

```diff
! Recommended prior distributions 

!  -hosp-frac-10 = [0.0 - 0.2]
!  -hosp-frac-20 = [0.0 - 0.2]
!  -hosp-frac-30 = [0.0 - 0.2]
!  -hosp-frac-40 = [0.0 - 0.5]
!  -hosp-frac-50 = [0.0 - 0.5]
!  -hosp-frac-60 = [0.0 - 1.0]
!  -hosp-frac-70 = [0.0 - 1.0]
!  -hosp-frac-80 = [0.0 - 1.0]

```


    
##### -introday

The introduction day of the first case. Feb 1 2020 is day 32, March 1 2020 is day 61, etc.


##### -loc

The location currently being modeled. E.g.

    `./odesim none -loc PA`
    
for Pennsylvania. Acceptable arguments currently are PA, MA, and RI.  Default is RI.


##### -mean-time-vent

This is the mean time that a surviving patient spends on a ventilator.  Default is set to 10.8 days.


##### -prob-icu-vent

This is the probability that a patient in the ICU will progress to mechanical ventilation.  Reasonable values are 0.75 to 0.99.


##### -rel-beta-hosp

This is the relative beta, i.e. the relative population mixing parameter, for a hospitalized individual (including ICU and ventilated). Default is 0.2. Must be between 0.0 and 1.0. Note that this beta does not change when social distancing is enforced.

```diff
! Recommended prior distribution [0.0 - 0.3]
```


##### -remove99

This option (or flag) removes the first 99 columns of the output.  These are the columns that keep track of the susceptible individuals, exposed individuals, and asymptomatic individuals.  These individuals are never enterred into the likelihood function, so it is not necessary to output them.



##### -self-isolation-factor

A factor describing the extent to which symptomatic individuals self-isolate by not leaving home and contacting others.  1.0 means that they behave identically to non-symptomatic individuals.  0.0 means that they self-isolate immediately after the first symptom appears and do not contact anyone until either (1) they need to be hospitalized or (2) all symptoms have disappeared.


##### -steps-per-day

This is the number of steps the Runge-Kutta integrator takes per day.  Default is set to 100.  Can set to 5 or 10 to improve speed in trial inference runs.



##### -susc-0-20, -susc-60-100

This is the relative susceptibility of certain age groups. Normally, the 30-39 age group is our reference group.  Here, the entire 20-59 age group is the reference group, since there is unlikely to be enough information to differentiate susceptibility among 10-year age bands.  The susceptibility of the 20-59 age group is 1.0, and the two other age-group susceptibilities may be able to be inferred.  Default settings are that the susceptibility of the 0-20 age group is 0.6 (slight adjustment from an odds-value of 0.51) and the susceptibility of the 60+ age group is 1.0.

```diff
! Recommended prior distributions 

!  -susc-0-20   = [0.2 - 1.5]
!  -susc-60-100 = [0.5 - 1.5]

```



##### -symp-frac

This is the probability that a 30-39 year-old infected individual will develop symptoms. Default value is 0.25. If you set this higher than 0.325, it will be reset back down to 0.325.

```diff
! Recommended prior distribution [0.1 - 0.3]
```

##### -tf
Time at which the ODEs are stopped. Day 1 is Jan 1 2020.  So, if you want to run the simulation through to April 30 2020, you would call

   `./odesim none -tf 121`
   
The 'time initial' right now is set to zero, which means that it is set to Jan 1 2020.  Cases are introduced via the `-introday` command-line argument.  


##### -time-symp-to-hosp

This is the average time (in days) between symptoms appearance and admission to hospital. It is the same for all ages right now.

```diff
! Recommended prior distribution [3.0 - 10.0]
```
