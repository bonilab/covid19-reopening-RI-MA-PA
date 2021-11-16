#!/usr/bin/env Rscript

### ms-figures.R
### last edited: 15 Nov 2021
### authors: Nathan Wikle
###
### Generates Figures 2 and 5 for manuscript and supplementary 
###   material, using CSV output from state-level runs. CAUTION: some 
###   state-level params are hard-coded into 'fig2panel'. These were current
###   for Sept 6 runs, as of 15 Nov 2021. They may need to be updated for 
###   future runs.




#################################################################################
### 1. Functions to help clean and process results before plotting
#################################################################################

### Process ODESIM trajectory for plotting.
traj.fine.process <- function(
  traj,   # output from ODESIM
  loc     # state used in the analysis
){

  ### relevant indices for plotting:
  sympt.cum.idx <- 290:298 # cumulative symptomatic cases
  hosp.cum.idx <- 299:307 # cumulative hospitalizations
  home.deaths.idx <- 254:262 # home deaths (cumul?)
  hosp.deaths.idx <- 263:271 # hosp. deaths
  
  hosp.curr.idx <- c(137:172,245:253) # hospitalized, but not in ICU
  icu.curr.idx <- c(173:181,236:244) # in ICU but not on ventilator
  vent.curr.idx <- 182:235 # in ICU and on ventilator
  
  # number of days simulated
  nr <- nrow(traj)

  ############################
  ### 1. Symptomatic cases ###
  ############################

  # total cumulative symptomatic cases
  traj.sympt <- traj[, sympt.cum.idx]
  
  # combine 0-9 and 10-19 age classes for MA
  if (loc == "MA") {
    traj.sympt.ma <- traj.sympt[, 2:ncol(traj.sympt)]
    traj.sympt.ma[, 1] <- traj.sympt[, 1] + traj.sympt[, 2]
    traj.sympt <- traj.sympt.ma
  }
  
  # calculating new cases for each day
  traj.sympt.new <-- traj.sympt
  traj.sympt.new[-1, ] <- traj.sympt[-1, ] - traj.sympt[-nr, ]
  tot.sympt.new <- rowSums(traj.sympt.new)
  
  ###########################
  ### 2. Hospitalizations ###
  ###########################

  # cumulative hospitalized by age:
  traj.hosp.cum <- traj[, hosp.cum.idx]
  
  # calculating new hosp from traj for each day
  traj.hosp.new <- traj.hosp.cum
  traj.hosp.new[-1, ] <- traj.hosp.cum[-1, ] - traj.hosp.cum[-nr, ]
  
  if (loc == "MA") {
    traj.hosp.new.ma <- traj.hosp.new[, 2:ncol(traj.hosp.new)]
    traj.hosp.new.ma[, 1] <- traj.hosp.new[, 1] + traj.hosp.new[, 2]
    traj.hosp.new <- traj.hosp.new.ma

    traj.hosp.cum.ma <- traj.hosp.cum[, 2:ncol(traj.hosp.cum)]
    traj.hosp.cum.ma[, 1] <- traj.hosp.cum[, 1] + traj.hosp.cum[, 2]
    traj.hosp.cum <- traj.hosp.cum.ma
  }

  tot.hosp.new <- rowSums(traj.hosp.new)
  
  #######################################
  ### 3. Current Hosps, Vent, and ICU ###
  #######################################

  # for proposing new s2.hosp parameter
  hosp.curr <- traj[, hosp.curr.idx]
  tot.hosp.curr <- rowSums(hosp.curr)
  
  # for proposing new s2.icu parameter
  icu.curr <- traj[, icu.curr.idx]
  tot.icu.curr <- rowSums(icu.curr)
  
  # for proposing new s2.vent parameter
  vent.curr <- traj[, vent.curr.idx]
  tot.vent.curr <- rowSums(vent.curr)
  
  #########################################
  ### 4. Deaths (home, hosp, and total) ###
  #########################################

  # deaths
  home.deaths <- traj[, home.deaths.idx]
  hosp.deaths <- traj[, hosp.deaths.idx]
  tot.home.deaths.cum <- rowSums(home.deaths, na.rm = T)
  tot.hosp.deaths.cum <- rowSums(hosp.deaths, na.rm = T)
  deaths.cum <- home.deaths + hosp.deaths
  
  home.deaths.new <- home.deaths
  home.deaths.new[-1, ] <- home.deaths[-1, ] - home.deaths[-nr, ]
  
  hosp.deaths.new <- hosp.deaths
  hosp.deaths.new[-1, ] <- hosp.deaths[-1, ] - hosp.deaths[-nr, ]
  
  # total (new) deaths
  deaths.new <- home.deaths.new + hosp.deaths.new
  
  tot.home.deaths.new <- rowSums(home.deaths.new)
  tot.hosp.deaths.new <- rowSums(hosp.deaths.new)
  tot.deaths.new <- tot.home.deaths.new + tot.hosp.deaths.new
  
  # combine 0-9 and 10-19 age classes for MA
  if (loc == "MA") {
    deaths.cum.ma <- deaths.cum[, 2:ncol(deaths.cum)]
    deaths.cum.ma[, 1] <- deaths.cum[, 1] + deaths.cum[, 2]
    deaths.cum <- deaths.cum.ma

    deaths.new.ma <- deaths.new[, 2:ncol(deaths.new)]
    deaths.new.ma[, 1] <- deaths.new[, 1] + deaths.new[, 2]
    deaths.new <- deaths.new.ma
  }
  
  ##############################
  ### 5. Hospital Discharges ###
  ##############################
  
  discharge.idx <- 281:289 
  hosp.discharges <- traj[,discharge.idx]
  hosp.discharges.new <- hosp.discharges
  hosp.discharges.new[-1, ] <- hosp.discharges[-1, ] - hosp.discharges[-nr, ]
  
  tot.hosp.discharges <- apply(traj[, discharge.idx], 1, sum, na.rm = T)
  tot.hosp.discharges.new <- tot.hosp.discharges
  tot.hosp.discharges.new[-1] <- tot.hosp.discharges[-1] - tot.hosp.discharges[-nr]

  # return list with processed results
  list(
    days.odesim = traj[, 1],
    # cases
    tot.sympt.new.odesim = as.numeric(tot.sympt.new),
    sympt.new.odesim = as.matrix(traj.sympt.new),
    tot.sympt.cum.odesim = as.numeric(rowSums(traj.sympt, na.rm = T)),
    sympt.cum.odesim = as.matrix(traj.sympt),
    # hosp.
    tot.hosp.new.odesim = as.numeric(tot.hosp.new),
    hosp.new.odesim = as.matrix(traj.hosp.new),
    tot.hosp.cum.odesim = as.numeric(rowSums(traj.hosp.cum, na.rm = T)),
    hosp.cum.odesim = as.matrix(traj.hosp.cum),
    # deaths
    tot.deaths.new.odesim = as.numeric(tot.deaths.new),
    deaths.new.odesim = as.matrix(deaths.new),
    tot.deaths.cum.odesim = as.numeric(apply(deaths.cum, 1, sum, na.rm = TRUE)),
    deaths.cum.odesim = as.matrix(deaths.cum),
    hosp.deaths.new.odesim = as.matrix(hosp.deaths.new),
    home.deaths.new.odesim = as.matrix(home.deaths.new),
    tot.hosp.deaths.new.odesim = as.numeric(tot.hosp.deaths.new),
    tot.home.deaths.new.odesim = as.numeric(tot.home.deaths.new),
    tot.hosp.deaths.cum.odesim = as.numeric(tot.hosp.deaths.cum),
    tot.home.deaths.cum.odesim = as.numeric(tot.home.deaths.cum),
    # current hosp, icu, vent
    tot.hosp.curr.noicu.odesim = tot.hosp.curr,
    tot.icu.curr.novent.odesim = tot.icu.curr,
    tot.vent.curr.odesim = tot.vent.curr,
    tot.hosp.curr.odesim = tot.hosp.curr + tot.icu.curr + tot.vent.curr,
    tot.icu.curr.odesim = tot.icu.curr + tot.vent.curr,
    # hosp. discharges
    tot.hosp.discharges.cum.odesim = as.numeric(tot.hosp.discharges),
    hosp.discharges.cum.odesim = as.matrix(hosp.discharges),
    tot.hosp.discharges.new.odesim = as.numeric(tot.hosp.discharges.new),
    hosp.discharges.new.odesim = as.matrix(hosp.discharges.new)
  )
}

### Process observed data for plotting.
data.fine.process <- function(
  df,   # observed data
  loc   # state used in the analysis
){

  # days after 2020 Jan 01  
  days <- df$daynum 
  n.days <- nrow(df)
  
  ### calculations for age-structured infill and pre-processing of data

  # cases
  cases.cum <- df$cumulative_confirmed_cases
  cases.new <- c(cases.cum[1], cases.cum[-1] - cases.cum[-n.days])
  
  if (loc == "MA") {
    cases.age.cum <- df[, 16:23]
    hosp.age.cum <- df[, 24:31]
    deaths.age.cum <- df[, 32:39]
  } else {
    cases.age.cum <- df[, 16:24]
    hosp.age.cum <- df[, 25:33]
    deaths.age.cum <- df[, 34:42]
  }
  
  cases.age.new <- as.matrix(rbind(
    cases.age.cum[1, ],
    cases.age.cum[-1, ] - cases.age.cum[-n.days, ]
  ))
  cases.age.new[cases.age.new < 0] <- 0
  
  # hospitalizations
  hosp.age.new <- as.matrix(rbind(
    hosp.age.cum[1, ],
    hosp.age.cum[-1, ] - hosp.age.cum[-n.days, ]
  ))
  hosp.age.new[hosp.age.new < 0] <- 0
  hosp.cum <- df$hospitalized_cumulative
  hosp.new <- c(hosp.cum[1], hosp.cum[-1] - hosp.cum[-n.days])

  if (loc == "CT") {
    hosp.curr <- df$current_hospitalized
  } else {
    hosp.curr <- df$hospitalized_currently
  }
  
  # icu
  icu.curr <- df$InICU_currently

  # # ventilator
  if (loc == "CT") {
    vent.curr <- df$OnVentilator_currently
  } else {
    vent.curr <- df$OnVentilator_Currently
  }
  
  # deaths
  tot.deaths.cum <- df$cumulative_deaths
  tot.deaths.new <- c(
    tot.deaths.cum[1],
    tot.deaths.cum[-1] - tot.deaths.cum[-n.days]
  )
  
  # deaths in hospital
  if (loc == "CT") {
    deaths.hosp.cum <- df$cumulative_hospital_deaths
  } else {
    deaths.hosp.cum <- df$Cumulative_hospital_deaths
  }
  deaths.home.cum <- df$cumulative_deaths - deaths.hosp.cum
  deaths.hosp.new <- c(
    deaths.hosp.cum[1],
    deaths.hosp.cum[-1] - deaths.hosp.cum[-n.days]
  )
  deaths.out.of.hosp.new <- tot.deaths.new - deaths.hosp.new
  if (is.null(deaths.hosp.cum[1])) {
    deaths.hosp.new <- NULL
    deaths.out.of.hosp.new <- NULL
  }
  deaths.age.new <- as.matrix(rbind(
    deaths.age.cum[1, ],
    deaths.age.cum[-1, ] - deaths.age.cum[-n.days, ]
  ))
  
  # hospital discharges
  tot.hosp.discharges <- df$Cumulative_hospital_discharges
  tot.new.hosp.discharges <- c(0, tot.hosp.discharges[-1] - tot.hosp.discharges[-n.days])
  
  # active surveillance
  if (loc == "MA") {
    tests.as <- df$nursinghome_tests_total
    pos.as <- df$nursinghome_tests_positive
  } else {
    tests.as <- rep(NA, nrow(df))
    pos.as <- rep(NA, nrow(df))
  }
  
  # return processed data
  list(
    days = days,
    # symptomatic cases
    tot.sympt.new = cases.new,
    tot.sympt.cum = cases.cum,
    sympt.new = cases.age.new,
    sympt.cum = cases.age.cum,
    # hospitalizations
    tot.hosp.new = hosp.new,
    hosp.new = hosp.age.new,
    tot.hosp.cum = hosp.cum,
    hosp.cum = hosp.age.cum,
    # current hosp, icu, and vent data
    tot.hosp.curr = hosp.curr,
    tot.icu.curr = icu.curr,
    tot.vent.curr = vent.curr,
    # death data
    tot.deaths.new = tot.deaths.new,
    tot.hosp.deaths.new = deaths.hosp.new,
    tot.hosp.deaths.cum = deaths.hosp.cum,
    tot.home.deaths.new = deaths.out.of.hosp.new,
    tot.home.deaths.cum = deaths.home.cum,
    deaths.cum = deaths.age.cum,
    deaths.new = deaths.age.new,
    tot.deaths.cum = tot.deaths.cum,
    # hospital discharges
    tot.hosp.discharges.new = tot.new.hosp.discharges,
    tot.hosp.discharges.cum = tot.hosp.discharges,
    tests.as = tests.as,
    pos.as = pos.as
  )
}

### Generate ODESIM trajectories for a given set of posterior samples.
traj.sim <- function(
  samples,        # list with posterior samples, used to generate ODESIM trajectories
  odepath,        # path to ODESIM directory
  csv = TRUE      # whether samples were generated using a CSV or Rdata file
){
  
  # number of samples available for traj sim
  n.samples <- nrow(samples$betas)
  
  for (k in 1:n.samples) {
    
    # last day of simulated output
    end.day <- max(samples$days) 
    
    if (csv) {
      n.b <- ncol(samples$betas)
      beta.daily <- as.vector(samples$betas[k, ])
    } else {
      # beta daily vector
      beta.daily <- samples$Z.beta %*% samples$betas[k, ]
    }
    
    if (k == 1) {
      beta.full <- matrix(NA, nrow = n.samples, ncol = length(beta.daily))
      beta.full[k, ] <- beta.daily
    } else {
      beta.full[k, ] <- beta.daily
    }
    
    # check for fitted hosp.report.rate parameter
    if (is.element("hosp.report.rate", colnames(samples$ode))) {
      idx.h <- which(colnames(samples$ode) == "hosp.report.rate")
      cumul.hosp.rr <- samples$ode[k, idx.h]
      ode.params <- samples$ode[k, -idx.h]
      names(ode.params) <- colnames(samples$ode)[-idx.h]
      constants <- samples$const
    } else {
      cumul.hosp.rr <- 1
      ode.params <- unlist(samples$ode[k, ])
      names(ode.params) <- colnames(samples$ode)
      constants <- samples$const
    }
    
    if (is.element("p.asympt", colnames(samples$ode))) {
      non.odesim <- c("p.asympt")
    } else {
      non.odesim  <- NULL
    }

    # generate trajectories
    if (samples$sf.choice) {
      sf.choice.names <- c(
        " ", "-symp-frac-davies", "-symp-frac-equal 0.3",
        "-symp-frac-equal 0.4", "-symp-frac-equal 0.5",
        "-symp-frac-equal 0.6", "-symp-frac-equal 0.7"
      )
      # generate trajectory
      traj.k <- traj.from.params(beta.daily,
        params = ode.params,
        tf = end.day,
        introday = samples$introday,
        const.params = constants,
        non.odesim.params = non.odesim,
        odepath = odepath,
        loc = samples$loc,
        symp = sf.choice.names[samples$sf.vals[k]]
      )
    } else {
      # generate trajectory
      traj.k <- traj.from.params(
        beta = beta.daily,
        params = ode.params,
        const.params = constants,
        non.odesim.params = non.odesim,
        introday = samples$introday,
        tf = end.day,
        odepath = odepath,
        loc = samples$loc,
        symp = NULL
      )
    }
    
    # process trajectory
    tp.k <- traj.fine.process(traj.k, loc = samples$loc)
    
    # create reporting rate vector of correct length
    if (!is.matrix(samples$Z.rr)) {
      rr.full <- as.vector(samples$Z.rr * samples$rr[k])
    } else if (ncol(samples$Z.rr) == 1) {
      rr.full <- as.vector(samples$Z.rr * samples$rr[k])
    } else {
      rr.full <- as.vector(samples$Z.rr %*% samples$rr[k, ])
    }
    
    if (k == 1) {
      rr.full.m <- matrix(NA, nrow = n.samples, ncol = length(rr.full))
      rr.full.m[k, ] <- rr.full
    } else {
      rr.full.m[k, ] <- rr.full
    }
    
    lower.rr <- rr.full[1]
    upper.rr <- rr.full[length(rr.full)]
    
    lower.rep <- which(tp.k$days.odesim < min(samples$df$daynum, na.rm = T))
    upper.rep <- which(tp.k$days.odesim > max(samples$df$daynum, na.rm = T))
    
    rr.full <- c(
      rep(lower.rr, length(lower.rep)),
      rr.full,
      rep(upper.rr, length(upper.rep))
    )
    
    if (csv) {
      rr.full <- rr.full[-length(rr.full)]
    }
    
    ###################
    ### FIXED DELAY ###
    ###################
    if (!is.null(samples$pres.delay)) {
      if (samples$pres.delay > 0) {
        P <- matrix(0, nrow = end.day + samples$pres.delay + 1, ncol = end.day + 1)
        for (j in 1:ncol(P)) {
          P[j:(j + samples$pres.delay), j] <- c(rep(0, samples$pres.delay), 1)
        }
      } else {
        P <- diag(end.day + 1)
      }
    } else {
      P <- diag(end.day + 1)
    }

    # ll.k <- loglik.odesim(
    #   traj = traj.k,
    #   df = samples$df, 
    #   P = P, 
    #   loc = samples$loc,
    #   report.rate = rr.full,
    #   nb.disp.params = c(0.9, 0.75, 0.6, 0.01, 0.8),
    #   s2.hosp = 30,
    #   s2.icu = 5,
    #   s2.vent = 4,
    #   lik.tot = FALSE, 
    #   lik.age = TRUE,
    #   lik.hosp.new = TRUE,
    #   lik.hosp.curr = TRUE,
    #   lik.icu.curr = TRUE,
    #   lik.vent.curr = TRUE,
    #   lik.tot.deaths = FALSE,
    #   lik.home.deaths = TRUE,
    #   lik.hosp.deaths = TRUE,
    #   lik.age.deaths = TRUE,
    #   lik.hosp.discharges = TRUE,
    #   lik.curr = FALSE,
    #   lik.old = TRUE,
    #   active.surv = FALSE,
    #   p.asympt = .4,
    #   total.size.constraint = FALSE
    # )
    
    # ll.full.k <- ll.k[1]
    
    if (k == 1){
      
      ### create structures
      
      # symptomatic cases (no age)
      tot.sympt.new <- matrix(NA, nrow = n.samples, ncol = length(tp.k$tot.sympt.new.odesim))
      tot.sympt.cum <- matrix(NA, nrow = n.samples, ncol = length(tp.k$tot.sympt.cum.odesim))
      
      # hospitalizations (no age)
      tot.hosp.new <- matrix(NA, nrow = n.samples, ncol = length(tp.k$tot.hosp.new.odesim))
      tot.hosp.cum <- matrix(NA, nrow = n.samples, ncol = length(tp.k$tot.hosp.cum.odesim))
      
      # deaths (no age)
      tot.deaths.new <- matrix(NA, nrow = n.samples, ncol = length(tp.k$tot.deaths.new.odesim))
      tot.deaths.cum <- matrix(NA, nrow = n.samples, ncol = length(tp.k$tot.deaths.cum.odesim))
      tot.hosp.deaths.new <- matrix(NA, nrow = n.samples, ncol = length(tp.k$tot.hosp.deaths.new.odesim))
      tot.home.deaths.new <- matrix(NA, nrow = n.samples, ncol = length(tp.k$tot.home.deaths.new.odesim))
      tot.hosp.deaths.cum <- matrix(NA, nrow = n.samples, ncol = length(tp.k$tot.hosp.deaths.cum.odesim))
      tot.home.deaths.cum <- matrix(NA, nrow = n.samples, ncol = length(tp.k$tot.home.deaths.cum.odesim))
      
      # current (no age)
      tot.hosp.curr.noicu <- matrix(NA, nrow = n.samples, ncol = length(tp.k$tot.hosp.curr.noicu.odesim))
      tot.icu.curr.novent <- matrix(NA, nrow = n.samples, ncol = length(tp.k$tot.icu.curr.novent.odesim))
      tot.vent.curr <- matrix(NA, nrow = n.samples, ncol = length(tp.k$tot.vent.curr.odesim))
      tot.hosp.curr <- matrix(NA, nrow = n.samples, ncol = length(tp.k$tot.hosp.curr.odesim))
      tot.icu.curr <- matrix(NA, nrow = n.samples, ncol = length(tp.k$tot.icu.curr.odesim))
      
      # discharges
      tot.hosp.discharges.new <- matrix(NA, nrow = n.samples, ncol = length(tp.k$tot.hosp.discharges.new.odesim))
      tot.hosp.discharges.cum <- matrix(NA, nrow = n.samples, ncol = length(tp.k$tot.hosp.discharges.cum.odesim))
      
      # ll.values <- rep(NA, length = n.samples)
      
      if (samples$loc == "MA"){
        
        ### symptomatic
        
        # new symptomatic cases (by age)
        sympt.new.01 <- matrix(NA, nrow = n.samples, ncol = length(tp.k$sympt.new.odesim[,1]))
        sympt.new.2 <- matrix(NA, nrow = n.samples, ncol = length(tp.k$sympt.new.odesim[,2]))
        sympt.new.3 <- matrix(NA, nrow = n.samples, ncol = length(tp.k$sympt.new.odesim[,3]))
        sympt.new.4 <- matrix(NA, nrow = n.samples, ncol = length(tp.k$sympt.new.odesim[,4]))
        sympt.new.5 <- matrix(NA, nrow = n.samples, ncol = length(tp.k$sympt.new.odesim[,5]))
        sympt.new.6 <- matrix(NA, nrow = n.samples, ncol = length(tp.k$sympt.new.odesim[,6]))
        sympt.new.7 <- matrix(NA, nrow = n.samples, ncol = length(tp.k$sympt.new.odesim[,7]))
        sympt.new.8 <- matrix(NA, nrow = n.samples, ncol = length(tp.k$sympt.new.odesim[,8]))
        
        # cum. symptomatic cases (by age)
        sympt.cum.01 <- matrix(NA, nrow = n.samples, ncol = length(tp.k$sympt.cum.odesim[,1]))
        sympt.cum.2 <- matrix(NA, nrow = n.samples, ncol = length(tp.k$sympt.cum.odesim[,2]))
        sympt.cum.3 <- matrix(NA, nrow = n.samples, ncol = length(tp.k$sympt.cum.odesim[,3]))
        sympt.cum.4 <- matrix(NA, nrow = n.samples, ncol = length(tp.k$sympt.cum.odesim[,4]))
        sympt.cum.5 <- matrix(NA, nrow = n.samples, ncol = length(tp.k$sympt.cum.odesim[,5]))
        sympt.cum.6 <- matrix(NA, nrow = n.samples, ncol = length(tp.k$sympt.cum.odesim[,6]))
        sympt.cum.7 <- matrix(NA, nrow = n.samples, ncol = length(tp.k$sympt.cum.odesim[,7]))
        sympt.cum.8 <- matrix(NA, nrow = n.samples, ncol = length(tp.k$sympt.cum.odesim[,8]))
        
        ### hospitalizations 
        
        # new hospitalization (by age)
        hosp.new.01 <- matrix(NA, nrow = n.samples, ncol = length(tp.k$hosp.new.odesim[,1]))
        hosp.new.2 <- matrix(NA, nrow = n.samples, ncol = length(tp.k$hosp.new.odesim[,2]))
        hosp.new.3 <- matrix(NA, nrow = n.samples, ncol = length(tp.k$hosp.new.odesim[,3]))
        hosp.new.4 <- matrix(NA, nrow = n.samples, ncol = length(tp.k$hosp.new.odesim[,4]))
        hosp.new.5 <- matrix(NA, nrow = n.samples, ncol = length(tp.k$hosp.new.odesim[,5]))
        hosp.new.6 <- matrix(NA, nrow = n.samples, ncol = length(tp.k$hosp.new.odesim[,6]))
        hosp.new.7 <- matrix(NA, nrow = n.samples, ncol = length(tp.k$hosp.new.odesim[,7]))
        hosp.new.8 <- matrix(NA, nrow = n.samples, ncol = length(tp.k$hosp.new.odesim[,8]))
        
        # cum. hospitalizations (by age)
        hosp.cum.01 <- matrix(NA, nrow = n.samples, ncol = length(tp.k$hosp.cum.odesim[,1]))
        hosp.cum.2 <- matrix(NA, nrow = n.samples, ncol = length(tp.k$hosp.cum.odesim[,2]))
        hosp.cum.3 <- matrix(NA, nrow = n.samples, ncol = length(tp.k$hosp.cum.odesim[,3]))
        hosp.cum.4 <- matrix(NA, nrow = n.samples, ncol = length(tp.k$hosp.cum.odesim[,4]))
        hosp.cum.5 <- matrix(NA, nrow = n.samples, ncol = length(tp.k$hosp.cum.odesim[,5]))
        hosp.cum.6 <- matrix(NA, nrow = n.samples, ncol = length(tp.k$hosp.cum.odesim[,6]))
        hosp.cum.7 <- matrix(NA, nrow = n.samples, ncol = length(tp.k$hosp.cum.odesim[,7]))
        hosp.cum.8 <- matrix(NA, nrow = n.samples, ncol = length(tp.k$hosp.cum.odesim[,8]))
        
        ### deaths
        
        ### new deaths
        deaths.new.01 <- matrix(NA, nrow = n.samples, ncol = length(tp.k$deaths.new.odesim[,1]))
        deaths.new.2 <- matrix(NA, nrow = n.samples, ncol = length(tp.k$deaths.new.odesim[,2]))
        deaths.new.3 <- matrix(NA, nrow = n.samples, ncol = length(tp.k$deaths.new.odesim[,3]))
        deaths.new.4 <- matrix(NA, nrow = n.samples, ncol = length(tp.k$deaths.new.odesim[,4]))
        deaths.new.5 <- matrix(NA, nrow = n.samples, ncol = length(tp.k$deaths.new.odesim[,5]))
        deaths.new.6 <- matrix(NA, nrow = n.samples, ncol = length(tp.k$deaths.new.odesim[,6]))
        deaths.new.7 <- matrix(NA, nrow = n.samples, ncol = length(tp.k$deaths.new.odesim[,7]))
        deaths.new.8 <- matrix(NA, nrow = n.samples, ncol = length(tp.k$deaths.new.odesim[,8]))
        
        ### cum. deaths
        deaths.cum.01 <- matrix(NA, nrow = n.samples, ncol = length(tp.k$deaths.cum.odesim[,1]))
        deaths.cum.2 <- matrix(NA, nrow = n.samples, ncol = length(tp.k$deaths.cum.odesim[,2]))
        deaths.cum.3 <- matrix(NA, nrow = n.samples, ncol = length(tp.k$deaths.cum.odesim[,3]))
        deaths.cum.4 <- matrix(NA, nrow = n.samples, ncol = length(tp.k$deaths.cum.odesim[,4]))
        deaths.cum.5 <- matrix(NA, nrow = n.samples, ncol = length(tp.k$deaths.cum.odesim[,5]))
        deaths.cum.6 <- matrix(NA, nrow = n.samples, ncol = length(tp.k$deaths.cum.odesim[,6]))
        deaths.cum.7 <- matrix(NA, nrow = n.samples, ncol = length(tp.k$deaths.cum.odesim[,7]))
        deaths.cum.8 <- matrix(NA, nrow = n.samples, ncol = length(tp.k$deaths.cum.odesim[,8]))
        
        ### hosp. deaths
        hosp.deaths.new.01 <- matrix(NA, nrow = n.samples, ncol = length(tp.k$hosp.deaths.new.odesim[,1]))
        hosp.deaths.new.2 <- matrix(NA, nrow = n.samples, ncol = length(tp.k$hosp.deaths.new.odesim[,2]))
        hosp.deaths.new.3 <- matrix(NA, nrow = n.samples, ncol = length(tp.k$hosp.deaths.new.odesim[,3]))
        hosp.deaths.new.4 <- matrix(NA, nrow = n.samples, ncol = length(tp.k$hosp.deaths.new.odesim[,4]))
        hosp.deaths.new.5 <- matrix(NA, nrow = n.samples, ncol = length(tp.k$hosp.deaths.new.odesim[,5]))
        hosp.deaths.new.6 <- matrix(NA, nrow = n.samples, ncol = length(tp.k$hosp.deaths.new.odesim[,6]))
        hosp.deaths.new.7 <- matrix(NA, nrow = n.samples, ncol = length(tp.k$hosp.deaths.new.odesim[,7]))
        hosp.deaths.new.8 <- matrix(NA, nrow = n.samples, ncol = length(tp.k$hosp.deaths.new.odesim[,8]))
        
        ### home deaths
        home.deaths.new.01 <- matrix(NA, nrow = n.samples, ncol = length(tp.k$home.deaths.new.odesim[,1]))
        home.deaths.new.2 <- matrix(NA, nrow = n.samples, ncol = length(tp.k$home.deaths.new.odesim[,2]))
        home.deaths.new.3 <- matrix(NA, nrow = n.samples, ncol = length(tp.k$home.deaths.new.odesim[,3]))
        home.deaths.new.4 <- matrix(NA, nrow = n.samples, ncol = length(tp.k$home.deaths.new.odesim[,4]))
        home.deaths.new.5 <- matrix(NA, nrow = n.samples, ncol = length(tp.k$home.deaths.new.odesim[,5]))
        home.deaths.new.6 <- matrix(NA, nrow = n.samples, ncol = length(tp.k$home.deaths.new.odesim[,6]))
        home.deaths.new.7 <- matrix(NA, nrow = n.samples, ncol = length(tp.k$home.deaths.new.odesim[,7]))
        home.deaths.new.8 <- matrix(NA, nrow = n.samples, ncol = length(tp.k$home.deaths.new.odesim[,8]))
        
        ### discharges
        
        ### new hosp discharges
        hosp.discharges.new.01 <- matrix(NA, nrow = n.samples, ncol = length(tp.k$hosp.discharges.new.odesim[,1]))
        hosp.discharges.new.2 <- matrix(NA, nrow = n.samples, ncol = length(tp.k$hosp.discharges.new.odesim[,2]))
        hosp.discharges.new.3 <- matrix(NA, nrow = n.samples, ncol = length(tp.k$hosp.discharges.new.odesim[,3]))
        hosp.discharges.new.4 <- matrix(NA, nrow = n.samples, ncol = length(tp.k$hosp.discharges.new.odesim[,4]))
        hosp.discharges.new.5 <- matrix(NA, nrow = n.samples, ncol = length(tp.k$hosp.discharges.new.odesim[,5]))
        hosp.discharges.new.6 <- matrix(NA, nrow = n.samples, ncol = length(tp.k$hosp.discharges.new.odesim[,6]))
        hosp.discharges.new.7 <- matrix(NA, nrow = n.samples, ncol = length(tp.k$hosp.discharges.new.odesim[,7]))
        hosp.discharges.new.8 <- matrix(NA, nrow = n.samples, ncol = length(tp.k$hosp.discharges.new.odesim[,8]))
        
        ### cum hosp discharges 
        hosp.discharges.cum.01 <- matrix(NA, nrow = n.samples, ncol = length(tp.k$hosp.discharges.cum.odesim[,1]))
        hosp.discharges.cum.2 <- matrix(NA, nrow = n.samples, ncol = length(tp.k$hosp.discharges.cum.odesim[,2]))
        hosp.discharges.cum.3 <- matrix(NA, nrow = n.samples, ncol = length(tp.k$hosp.discharges.cum.odesim[,3]))
        hosp.discharges.cum.4 <- matrix(NA, nrow = n.samples, ncol = length(tp.k$hosp.discharges.cum.odesim[,4]))
        hosp.discharges.cum.5 <- matrix(NA, nrow = n.samples, ncol = length(tp.k$hosp.discharges.cum.odesim[,5]))
        hosp.discharges.cum.6 <- matrix(NA, nrow = n.samples, ncol = length(tp.k$hosp.discharges.cum.odesim[,6]))
        hosp.discharges.cum.7 <- matrix(NA, nrow = n.samples, ncol = length(tp.k$hosp.discharges.cum.odesim[,7]))
        hosp.discharges.cum.8 <- matrix(NA, nrow = n.samples, ncol = length(tp.k$hosp.discharges.cum.odesim[,8]))
        
      } else {
        
        ### symptomatic
        
        # new symptomatic cases (by age)
        sympt.new.0 <- matrix(NA, nrow = n.samples, ncol = length(tp.k$sympt.new.odesim[,1]))
        sympt.new.1 <- matrix(NA, nrow = n.samples, ncol = length(tp.k$sympt.new.odesim[,2]))
        sympt.new.2 <- matrix(NA, nrow = n.samples, ncol = length(tp.k$sympt.new.odesim[,3]))
        sympt.new.3 <- matrix(NA, nrow = n.samples, ncol = length(tp.k$sympt.new.odesim[,4]))
        sympt.new.4 <- matrix(NA, nrow = n.samples, ncol = length(tp.k$sympt.new.odesim[,5]))
        sympt.new.5 <- matrix(NA, nrow = n.samples, ncol = length(tp.k$sympt.new.odesim[,6]))
        sympt.new.6 <- matrix(NA, nrow = n.samples, ncol = length(tp.k$sympt.new.odesim[,7]))
        sympt.new.7 <- matrix(NA, nrow = n.samples, ncol = length(tp.k$sympt.new.odesim[,8]))
        sympt.new.8 <- matrix(NA, nrow = n.samples, ncol = length(tp.k$sympt.new.odesim[,9]))
        
        # cum. symptomatic cases (by age)
        sympt.cum.0 <- matrix(NA, nrow = n.samples, ncol = length(tp.k$sympt.cum.odesim[,1]))
        sympt.cum.1 <- matrix(NA, nrow = n.samples, ncol = length(tp.k$sympt.cum.odesim[,2]))
        sympt.cum.2 <- matrix(NA, nrow = n.samples, ncol = length(tp.k$sympt.cum.odesim[,3]))
        sympt.cum.3 <- matrix(NA, nrow = n.samples, ncol = length(tp.k$sympt.cum.odesim[,4]))
        sympt.cum.4 <- matrix(NA, nrow = n.samples, ncol = length(tp.k$sympt.cum.odesim[,5]))
        sympt.cum.5 <- matrix(NA, nrow = n.samples, ncol = length(tp.k$sympt.cum.odesim[,6]))
        sympt.cum.6 <- matrix(NA, nrow = n.samples, ncol = length(tp.k$sympt.cum.odesim[,7]))
        sympt.cum.7 <- matrix(NA, nrow = n.samples, ncol = length(tp.k$sympt.cum.odesim[,8]))
        sympt.cum.8 <- matrix(NA, nrow = n.samples, ncol = length(tp.k$sympt.cum.odesim[,9]))
        
        ### hospitalizations 
        
        # new hospitalization (by age)
        hosp.new.0 <- matrix(NA, nrow = n.samples, ncol = length(tp.k$hosp.new.odesim[,1]))
        hosp.new.1 <- matrix(NA, nrow = n.samples, ncol = length(tp.k$hosp.new.odesim[,2]))
        hosp.new.2 <- matrix(NA, nrow = n.samples, ncol = length(tp.k$hosp.new.odesim[,3]))
        hosp.new.3 <- matrix(NA, nrow = n.samples, ncol = length(tp.k$hosp.new.odesim[,4]))
        hosp.new.4 <- matrix(NA, nrow = n.samples, ncol = length(tp.k$hosp.new.odesim[,5]))
        hosp.new.5 <- matrix(NA, nrow = n.samples, ncol = length(tp.k$hosp.new.odesim[,6]))
        hosp.new.6 <- matrix(NA, nrow = n.samples, ncol = length(tp.k$hosp.new.odesim[,7]))
        hosp.new.7 <- matrix(NA, nrow = n.samples, ncol = length(tp.k$hosp.new.odesim[,8]))
        hosp.new.8 <- matrix(NA, nrow = n.samples, ncol = length(tp.k$hosp.new.odesim[,9]))
        
        # cum. hospitalizations (by age)
        hosp.cum.0 <- matrix(NA, nrow = n.samples, ncol = length(tp.k$hosp.cum.odesim[,1]))
        hosp.cum.1 <- matrix(NA, nrow = n.samples, ncol = length(tp.k$hosp.cum.odesim[,2]))
        hosp.cum.2 <- matrix(NA, nrow = n.samples, ncol = length(tp.k$hosp.cum.odesim[,3]))
        hosp.cum.3 <- matrix(NA, nrow = n.samples, ncol = length(tp.k$hosp.cum.odesim[,4]))
        hosp.cum.4 <- matrix(NA, nrow = n.samples, ncol = length(tp.k$hosp.cum.odesim[,5]))
        hosp.cum.5 <- matrix(NA, nrow = n.samples, ncol = length(tp.k$hosp.cum.odesim[,6]))
        hosp.cum.6 <- matrix(NA, nrow = n.samples, ncol = length(tp.k$hosp.cum.odesim[,7]))
        hosp.cum.7 <- matrix(NA, nrow = n.samples, ncol = length(tp.k$hosp.cum.odesim[,8]))
        hosp.cum.8 <- matrix(NA, nrow = n.samples, ncol = length(tp.k$hosp.cum.odesim[,9]))
        
        ### deaths
        
        ### new deaths
        deaths.new.0 <- matrix(NA, nrow = n.samples, ncol = length(tp.k$deaths.new.odesim[,1]))
        deaths.new.1 <- matrix(NA, nrow = n.samples, ncol = length(tp.k$deaths.new.odesim[,2]))
        deaths.new.2 <- matrix(NA, nrow = n.samples, ncol = length(tp.k$deaths.new.odesim[,3]))
        deaths.new.3 <- matrix(NA, nrow = n.samples, ncol = length(tp.k$deaths.new.odesim[,4]))
        deaths.new.4 <- matrix(NA, nrow = n.samples, ncol = length(tp.k$deaths.new.odesim[,5]))
        deaths.new.5 <- matrix(NA, nrow = n.samples, ncol = length(tp.k$deaths.new.odesim[,6]))
        deaths.new.6 <- matrix(NA, nrow = n.samples, ncol = length(tp.k$deaths.new.odesim[,7]))
        deaths.new.7 <- matrix(NA, nrow = n.samples, ncol = length(tp.k$deaths.new.odesim[,8]))
        deaths.new.8 <- matrix(NA, nrow = n.samples, ncol = length(tp.k$deaths.new.odesim[,9]))
        
        ### cum. deaths
        deaths.cum.0 <- matrix(NA, nrow = n.samples, ncol = length(tp.k$deaths.cum.odesim[,1]))
        deaths.cum.1 <- matrix(NA, nrow = n.samples, ncol = length(tp.k$deaths.cum.odesim[,2]))
        deaths.cum.2 <- matrix(NA, nrow = n.samples, ncol = length(tp.k$deaths.cum.odesim[,3]))
        deaths.cum.3 <- matrix(NA, nrow = n.samples, ncol = length(tp.k$deaths.cum.odesim[,4]))
        deaths.cum.4 <- matrix(NA, nrow = n.samples, ncol = length(tp.k$deaths.cum.odesim[,5]))
        deaths.cum.5 <- matrix(NA, nrow = n.samples, ncol = length(tp.k$deaths.cum.odesim[,6]))
        deaths.cum.6 <- matrix(NA, nrow = n.samples, ncol = length(tp.k$deaths.cum.odesim[,7]))
        deaths.cum.7 <- matrix(NA, nrow = n.samples, ncol = length(tp.k$deaths.cum.odesim[,8]))
        deaths.cum.8 <- matrix(NA, nrow = n.samples, ncol = length(tp.k$deaths.cum.odesim[,9]))
        
        ### hosp. deaths
        hosp.deaths.new.0 <- matrix(NA, nrow = n.samples, ncol = length(tp.k$hosp.deaths.new.odesim[,1]))
        hosp.deaths.new.1 <- matrix(NA, nrow = n.samples, ncol = length(tp.k$hosp.deaths.new.odesim[,2]))
        hosp.deaths.new.2 <- matrix(NA, nrow = n.samples, ncol = length(tp.k$hosp.deaths.new.odesim[,3]))
        hosp.deaths.new.3 <- matrix(NA, nrow = n.samples, ncol = length(tp.k$hosp.deaths.new.odesim[,4]))
        hosp.deaths.new.4 <- matrix(NA, nrow = n.samples, ncol = length(tp.k$hosp.deaths.new.odesim[,5]))
        hosp.deaths.new.5 <- matrix(NA, nrow = n.samples, ncol = length(tp.k$hosp.deaths.new.odesim[,6]))
        hosp.deaths.new.6 <- matrix(NA, nrow = n.samples, ncol = length(tp.k$hosp.deaths.new.odesim[,7]))
        hosp.deaths.new.7 <- matrix(NA, nrow = n.samples, ncol = length(tp.k$hosp.deaths.new.odesim[,8]))
        hosp.deaths.new.8 <- matrix(NA, nrow = n.samples, ncol = length(tp.k$hosp.deaths.new.odesim[,9]))
        
        ### home deaths
        home.deaths.new.0 <- matrix(NA, nrow = n.samples, ncol = length(tp.k$home.deaths.new.odesim[,1]))
        home.deaths.new.1 <- matrix(NA, nrow = n.samples, ncol = length(tp.k$home.deaths.new.odesim[,2]))
        home.deaths.new.2 <- matrix(NA, nrow = n.samples, ncol = length(tp.k$home.deaths.new.odesim[,3]))
        home.deaths.new.3 <- matrix(NA, nrow = n.samples, ncol = length(tp.k$home.deaths.new.odesim[,4]))
        home.deaths.new.4 <- matrix(NA, nrow = n.samples, ncol = length(tp.k$home.deaths.new.odesim[,5]))
        home.deaths.new.5 <- matrix(NA, nrow = n.samples, ncol = length(tp.k$home.deaths.new.odesim[,6]))
        home.deaths.new.6 <- matrix(NA, nrow = n.samples, ncol = length(tp.k$home.deaths.new.odesim[,7]))
        home.deaths.new.7 <- matrix(NA, nrow = n.samples, ncol = length(tp.k$home.deaths.new.odesim[,8]))
        home.deaths.new.8 <- matrix(NA, nrow = n.samples, ncol = length(tp.k$home.deaths.new.odesim[,9]))
        
        ### discharges
        
        ### new hosp discharges
        hosp.discharges.new.0 <- matrix(NA, nrow = n.samples, ncol = length(tp.k$hosp.discharges.new.odesim[,1]))
        hosp.discharges.new.1 <- matrix(NA, nrow = n.samples, ncol = length(tp.k$hosp.discharges.new.odesim[,2]))
        hosp.discharges.new.2 <- matrix(NA, nrow = n.samples, ncol = length(tp.k$hosp.discharges.new.odesim[,3]))
        hosp.discharges.new.3 <- matrix(NA, nrow = n.samples, ncol = length(tp.k$hosp.discharges.new.odesim[,4]))
        hosp.discharges.new.4 <- matrix(NA, nrow = n.samples, ncol = length(tp.k$hosp.discharges.new.odesim[,5]))
        hosp.discharges.new.5 <- matrix(NA, nrow = n.samples, ncol = length(tp.k$hosp.discharges.new.odesim[,6]))
        hosp.discharges.new.6 <- matrix(NA, nrow = n.samples, ncol = length(tp.k$hosp.discharges.new.odesim[,7]))
        hosp.discharges.new.7 <- matrix(NA, nrow = n.samples, ncol = length(tp.k$hosp.discharges.new.odesim[,8]))
        hosp.discharges.new.8 <- matrix(NA, nrow = n.samples, ncol = length(tp.k$hosp.discharges.new.odesim[,9]))
        
        ### cum hosp discharges 
        hosp.discharges.cum.0 <- matrix(NA, nrow = n.samples, ncol = length(tp.k$hosp.discharges.cum.odesim[,1]))
        hosp.discharges.cum.1 <- matrix(NA, nrow = n.samples, ncol = length(tp.k$hosp.discharges.cum.odesim[,2]))
        hosp.discharges.cum.2 <- matrix(NA, nrow = n.samples, ncol = length(tp.k$hosp.discharges.cum.odesim[,3]))
        hosp.discharges.cum.3 <- matrix(NA, nrow = n.samples, ncol = length(tp.k$hosp.discharges.cum.odesim[,4]))
        hosp.discharges.cum.4 <- matrix(NA, nrow = n.samples, ncol = length(tp.k$hosp.discharges.cum.odesim[,5]))
        hosp.discharges.cum.5 <- matrix(NA, nrow = n.samples, ncol = length(tp.k$hosp.discharges.cum.odesim[,6]))
        hosp.discharges.cum.6 <- matrix(NA, nrow = n.samples, ncol = length(tp.k$hosp.discharges.cum.odesim[,7]))
        hosp.discharges.cum.7 <- matrix(NA, nrow = n.samples, ncol = length(tp.k$hosp.discharges.cum.odesim[,8]))
        hosp.discharges.cum.8 <- matrix(NA, nrow = n.samples, ncol = length(tp.k$hosp.discharges.cum.odesim[,9]))   
      }
    }  

    ### fill data structures

    # likelihood    
    # ll.values[k] <- ll.full.k

    # symptomatic cases (no age)
    traj.times <- round(traj.k[, 1])
    delay.rates.k <- P[, traj.times] %*% tp.k$tot.sympt.new.odesim
    tot.sympt.new[k, ] <- rr.full * delay.rates.k[traj.times]
    tot.sympt.cum[k, ] <- cumsum(tot.sympt.new[k, ])
      
    # hospitalizations (no age)
    tot.hosp.new[k, ] <- cumul.hosp.rr * tp.k$tot.hosp.new.odesim
    tot.hosp.cum[k, ] <- cumsum(tot.hosp.new[k, ])
      
    # deaths (no age)
    tot.deaths.new[k, ] <- tp.k$tot.deaths.new.odesim
    tot.deaths.cum[k, ] <- tp.k$tot.deaths.cum.odesim
    tot.hosp.deaths.new[k, ] <- tp.k$tot.hosp.deaths.new.odesim
    tot.home.deaths.new[k, ] <- tp.k$tot.home.deaths.new.odesim
    tot.hosp.deaths.cum[k, ] <- tp.k$tot.hosp.deaths.cum.odesim
    tot.home.deaths.cum[k, ] <- tp.k$tot.home.deaths.cum.odesim
      
    # current (no age)
    tot.hosp.curr.noicu[k, ] <- tp.k$tot.hosp.curr.noicu.odesim
    tot.icu.curr.novent[k, ] <- tp.k$tot.icu.curr.novent.odesim
    tot.vent.curr[k, ] <- tp.k$tot.vent.curr.odesim
    tot.hosp.curr[k, ] <- tp.k$tot.hosp.curr.odesim
    tot.icu.curr[k, ] <- tp.k$tot.icu.curr.odesim
    
    # discharges
    tot.hosp.discharges.new[k, ] <- tp.k$tot.hosp.discharges.new.odesim
    tot.hosp.discharges.cum[k, ] <- tp.k$tot.hosp.discharges.cum.odesim
     
    if (samples$loc == "MA"){
      
      ### symptomatic
      
      # new symptomatic cases (by age)
      sympt.new.01[k, ] <- rr.full * tp.k$sympt.new.odesim[, 1] 
      sympt.new.2[k, ] <- rr.full * tp.k$sympt.new.odesim[, 2] 
      sympt.new.3[k, ] <- rr.full * tp.k$sympt.new.odesim[, 3] 
      sympt.new.4[k, ] <- rr.full * tp.k$sympt.new.odesim[, 4] 
      sympt.new.5[k, ] <- rr.full * tp.k$sympt.new.odesim[, 5] 
      sympt.new.6[k, ] <- rr.full * tp.k$sympt.new.odesim[, 6] 
      sympt.new.7[k, ] <- rr.full * tp.k$sympt.new.odesim[, 7] 
      sympt.new.8[k, ] <- rr.full * tp.k$sympt.new.odesim[, 8] 
      
      # cum. symptomatic cases (by age)      
      sympt.cum.01[k, ] <- cumsum(sympt.new.01[k, ])
      sympt.cum.2[k, ] <- cumsum(sympt.new.2[k, ])
      sympt.cum.3[k, ] <- cumsum(sympt.new.3[k, ])
      sympt.cum.4[k, ] <- cumsum(sympt.new.4[k, ]) 
      sympt.cum.5[k, ] <- cumsum(sympt.new.5[k, ])
      sympt.cum.6[k, ] <- cumsum(sympt.new.6[k, ]) 
      sympt.cum.7[k, ] <- cumsum(sympt.new.7[k, ])
      sympt.cum.8[k, ] <- cumsum(sympt.new.8[k, ])
      
      ### hospitalizations 
      
      # new hospitalization (by age)
      hosp.new.01[k, ] <- cumul.hosp.rr * tp.k$hosp.new.odesim[, 1]
      hosp.new.2[k, ] <- cumul.hosp.rr * tp.k$hosp.new.odesim[, 2]
      hosp.new.3[k, ] <- cumul.hosp.rr * tp.k$hosp.new.odesim[, 3]
      hosp.new.4[k, ] <- cumul.hosp.rr * tp.k$hosp.new.odesim[, 4]
      hosp.new.5[k, ] <- cumul.hosp.rr * tp.k$hosp.new.odesim[, 5]
      hosp.new.6[k, ] <- cumul.hosp.rr * tp.k$hosp.new.odesim[, 6]
      hosp.new.7[k, ] <- cumul.hosp.rr * tp.k$hosp.new.odesim[, 7]
      hosp.new.8[k, ] <- cumul.hosp.rr * tp.k$hosp.new.odesim[, 8]
      
      # cum. hospitalizations (by age)
      hosp.cum.01[k, ] <- cumsum(hosp.new.01[k, ])
      hosp.cum.2[k, ] <- cumsum(hosp.new.2[k, ]) 
      hosp.cum.3[k, ] <- cumsum(hosp.new.3[k, ]) 
      hosp.cum.4[k, ] <- cumsum(hosp.new.4[k, ]) 
      hosp.cum.5[k, ] <- cumsum(hosp.new.5[k, ])  
      hosp.cum.6[k, ] <- cumsum(hosp.new.6[k, ]) 
      hosp.cum.7[k, ] <- cumsum(hosp.new.7[k, ])  
      hosp.cum.8[k, ] <- cumsum(hosp.new.8[k, ])  
      
      ### deaths
      
      ### new deaths
      deaths.new.01[k, ] <- tp.k$deaths.new.odesim[, 1]
      deaths.new.2[k, ] <- tp.k$deaths.new.odesim[, 2]
      deaths.new.3[k, ] <- tp.k$deaths.new.odesim[, 3]
      deaths.new.4[k, ] <- tp.k$deaths.new.odesim[, 4]
      deaths.new.5[k, ] <- tp.k$deaths.new.odesim[, 5]
      deaths.new.6[k, ] <- tp.k$deaths.new.odesim[, 6]
      deaths.new.7[k, ] <- tp.k$deaths.new.odesim[, 7]
      deaths.new.8[k, ] <- tp.k$deaths.new.odesim[, 8]
      
      ### cum. deaths
      deaths.cum.01[k, ] <- tp.k$deaths.cum.odesim[, 1]
      deaths.cum.2[k, ] <- tp.k$deaths.cum.odesim[, 2]
      deaths.cum.3[k, ] <- tp.k$deaths.cum.odesim[, 3]
      deaths.cum.4[k, ] <- tp.k$deaths.cum.odesim[, 4]
      deaths.cum.5[k, ] <- tp.k$deaths.cum.odesim[, 5]
      deaths.cum.6[k, ] <- tp.k$deaths.cum.odesim[, 6]
      deaths.cum.7[k, ] <- tp.k$deaths.cum.odesim[, 7]
      deaths.cum.8[k, ] <- tp.k$deaths.cum.odesim[, 8]
      
      ### hosp. deaths
      hosp.deaths.new.01[k, ] <- tp.k$hosp.deaths.new.odesim[, 1]
      hosp.deaths.new.2[k, ] <- tp.k$hosp.deaths.new.odesim[, 2]
      hosp.deaths.new.3[k, ] <- tp.k$hosp.deaths.new.odesim[, 3]
      hosp.deaths.new.4[k, ] <- tp.k$hosp.deaths.new.odesim[, 4]
      hosp.deaths.new.5[k, ] <- tp.k$hosp.deaths.new.odesim[, 5]
      hosp.deaths.new.6[k, ] <- tp.k$hosp.deaths.new.odesim[, 6]
      hosp.deaths.new.7[k, ] <- tp.k$hosp.deaths.new.odesim[, 7]
      hosp.deaths.new.8[k, ] <- tp.k$hosp.deaths.new.odesim[, 8]
      
      ### home deaths
      home.deaths.new.01[k, ] <- tp.k$home.deaths.new.odesim[, 1]
      home.deaths.new.2[k, ] <- tp.k$home.deaths.new.odesim[, 2]
      home.deaths.new.3[k, ] <- tp.k$home.deaths.new.odesim[, 3]
      home.deaths.new.4[k, ] <- tp.k$home.deaths.new.odesim[, 4]
      home.deaths.new.5[k, ] <- tp.k$home.deaths.new.odesim[, 5]
      home.deaths.new.6[k, ] <- tp.k$home.deaths.new.odesim[, 6]
      home.deaths.new.7[k, ] <- tp.k$home.deaths.new.odesim[, 7]
      home.deaths.new.8[k, ] <- tp.k$home.deaths.new.odesim[, 8]
      
      ### discharges
      
      ### new hosp discharges
      hosp.discharges.new.01[k, ] <- tp.k$hosp.discharges.new.odesim[, 1]
      hosp.discharges.new.2[k, ] <- tp.k$hosp.discharges.new.odesim[, 2]
      hosp.discharges.new.3[k, ] <- tp.k$hosp.discharges.new.odesim[, 3]
      hosp.discharges.new.4[k, ] <- tp.k$hosp.discharges.new.odesim[, 4]
      hosp.discharges.new.5[k, ] <- tp.k$hosp.discharges.new.odesim[, 5]
      hosp.discharges.new.6[k, ] <- tp.k$hosp.discharges.new.odesim[, 6]
      hosp.discharges.new.7[k, ] <- tp.k$hosp.discharges.new.odesim[, 7]
      hosp.discharges.new.8[k, ] <- tp.k$hosp.discharges.new.odesim[, 8]
      
      ### cum hosp discharges 
      hosp.discharges.cum.01[k, ] <- tp.k$hosp.discharges.cum.odesim[, 1]
      hosp.discharges.cum.2[k, ] <- tp.k$hosp.discharges.cum.odesim[, 2]
      hosp.discharges.cum.3[k, ] <- tp.k$hosp.discharges.cum.odesim[, 3]
      hosp.discharges.cum.4[k, ] <- tp.k$hosp.discharges.cum.odesim[, 4]
      hosp.discharges.cum.5[k, ] <- tp.k$hosp.discharges.cum.odesim[, 5]
      hosp.discharges.cum.6[k, ] <- tp.k$hosp.discharges.cum.odesim[, 6]
      hosp.discharges.cum.7[k, ] <- tp.k$hosp.discharges.cum.odesim[, 7]
      hosp.discharges.cum.8[k, ] <- tp.k$hosp.discharges.cum.odesim[, 8]
      
    } else {
      
      ### symptomatic
      
      # new symptomatic cases (by age)
      sympt.new.0[k, ] <- rr.full * tp.k$sympt.new.odesim[, 1] 
      sympt.new.1[k, ] <- rr.full * tp.k$sympt.new.odesim[, 2] 
      sympt.new.2[k, ] <- rr.full * tp.k$sympt.new.odesim[, 3] 
      sympt.new.3[k, ] <- rr.full * tp.k$sympt.new.odesim[, 4] 
      sympt.new.4[k, ] <- rr.full * tp.k$sympt.new.odesim[, 5] 
      sympt.new.5[k, ] <- rr.full * tp.k$sympt.new.odesim[, 6] 
      sympt.new.6[k, ] <- rr.full * tp.k$sympt.new.odesim[, 7] 
      sympt.new.7[k, ] <- rr.full * tp.k$sympt.new.odesim[, 8] 
      sympt.new.8[k, ] <- rr.full * tp.k$sympt.new.odesim[, 9] 
      
      # cum. symptomatic cases (by age)      
      sympt.cum.0[k, ] <- cumsum(sympt.new.0[k, ])
      sympt.cum.1[k, ] <- cumsum(sympt.new.1[k, ])
      sympt.cum.2[k, ] <- cumsum(sympt.new.2[k, ])
      sympt.cum.3[k, ] <- cumsum(sympt.new.3[k, ])
      sympt.cum.4[k, ] <- cumsum(sympt.new.4[k, ]) 
      sympt.cum.5[k, ] <- cumsum(sympt.new.5[k, ])
      sympt.cum.6[k, ] <- cumsum(sympt.new.6[k, ]) 
      sympt.cum.7[k, ] <- cumsum(sympt.new.7[k, ])
      sympt.cum.8[k, ] <- cumsum(sympt.new.8[k, ])
      
      ### hospitalizations 
      
      # new hospitalization (by age)
      hosp.new.0[k, ] <- cumul.hosp.rr * tp.k$hosp.new.odesim[, 1]
      hosp.new.1[k, ] <- cumul.hosp.rr * tp.k$hosp.new.odesim[, 2]
      hosp.new.2[k, ] <- cumul.hosp.rr * tp.k$hosp.new.odesim[, 3]
      hosp.new.3[k, ] <- cumul.hosp.rr * tp.k$hosp.new.odesim[, 4]
      hosp.new.4[k, ] <- cumul.hosp.rr * tp.k$hosp.new.odesim[, 5]
      hosp.new.5[k, ] <- cumul.hosp.rr * tp.k$hosp.new.odesim[, 6]
      hosp.new.6[k, ] <- cumul.hosp.rr * tp.k$hosp.new.odesim[, 7]
      hosp.new.7[k, ] <- cumul.hosp.rr * tp.k$hosp.new.odesim[, 8]
      hosp.new.8[k, ] <- cumul.hosp.rr * tp.k$hosp.new.odesim[, 9]
      
      # cum. hospitalizations (by age)      
      hosp.cum.0[k, ] <- cumsum(hosp.new.0[k, ])
      hosp.cum.1[k, ] <- cumsum(hosp.new.1[k, ])
      hosp.cum.2[k, ] <- cumsum(hosp.new.2[k, ])
      hosp.cum.3[k, ] <- cumsum(hosp.new.3[k, ])
      hosp.cum.4[k, ] <- cumsum(hosp.new.4[k, ])
      hosp.cum.5[k, ] <- cumsum(hosp.new.5[k, ])
      hosp.cum.6[k, ] <- cumsum(hosp.new.6[k, ])
      hosp.cum.7[k, ] <- cumsum(hosp.new.7[k, ])
      hosp.cum.8[k, ] <- cumsum(hosp.new.8[k, ])
      
      ### deaths
      
      ### new deaths
      deaths.new.0[k, ] <- tp.k$deaths.new.odesim[, 1]
      deaths.new.1[k, ] <- tp.k$deaths.new.odesim[, 2]
      deaths.new.2[k, ] <- tp.k$deaths.new.odesim[, 3]
      deaths.new.3[k, ] <- tp.k$deaths.new.odesim[, 4]
      deaths.new.4[k, ] <- tp.k$deaths.new.odesim[, 5]
      deaths.new.5[k, ] <- tp.k$deaths.new.odesim[, 6]
      deaths.new.6[k, ] <- tp.k$deaths.new.odesim[, 7]
      deaths.new.7[k, ] <- tp.k$deaths.new.odesim[, 8]
      deaths.new.8[k, ] <- tp.k$deaths.new.odesim[, 9]
      
      ### cum. deaths
      deaths.cum.0[k, ] <- tp.k$deaths.cum.odesim[, 1]
      deaths.cum.1[k, ] <- tp.k$deaths.cum.odesim[, 2]
      deaths.cum.2[k, ] <- tp.k$deaths.cum.odesim[, 3]
      deaths.cum.3[k, ] <- tp.k$deaths.cum.odesim[, 4]
      deaths.cum.4[k, ] <- tp.k$deaths.cum.odesim[, 5]
      deaths.cum.5[k, ] <- tp.k$deaths.cum.odesim[, 6]
      deaths.cum.6[k, ] <- tp.k$deaths.cum.odesim[, 7]
      deaths.cum.7[k, ] <- tp.k$deaths.cum.odesim[, 8]
      deaths.cum.8[k, ] <- tp.k$deaths.cum.odesim[, 9]
      
      ### hosp. deaths
      hosp.deaths.new.0[k, ] <- tp.k$hosp.deaths.new.odesim[, 1]
      hosp.deaths.new.1[k, ] <- tp.k$hosp.deaths.new.odesim[, 2]
      hosp.deaths.new.2[k, ] <- tp.k$hosp.deaths.new.odesim[, 3]
      hosp.deaths.new.3[k, ] <- tp.k$hosp.deaths.new.odesim[, 4]
      hosp.deaths.new.4[k, ] <- tp.k$hosp.deaths.new.odesim[, 5]
      hosp.deaths.new.5[k, ] <- tp.k$hosp.deaths.new.odesim[, 6]
      hosp.deaths.new.6[k, ] <- tp.k$hosp.deaths.new.odesim[, 7]
      hosp.deaths.new.7[k, ] <- tp.k$hosp.deaths.new.odesim[, 8]
      hosp.deaths.new.8[k, ] <- tp.k$hosp.deaths.new.odesim[, 9]
      
      ### home deaths
      home.deaths.new.0[k, ] <- tp.k$home.deaths.new.odesim[, 1]
      home.deaths.new.1[k, ] <- tp.k$home.deaths.new.odesim[, 2]
      home.deaths.new.2[k, ] <- tp.k$home.deaths.new.odesim[, 3]
      home.deaths.new.3[k, ] <- tp.k$home.deaths.new.odesim[, 4]
      home.deaths.new.4[k, ] <- tp.k$home.deaths.new.odesim[, 5]
      home.deaths.new.5[k, ] <- tp.k$home.deaths.new.odesim[, 6]
      home.deaths.new.6[k, ] <- tp.k$home.deaths.new.odesim[, 7]
      home.deaths.new.7[k, ] <- tp.k$home.deaths.new.odesim[, 8]
      home.deaths.new.8[k, ] <- tp.k$home.deaths.new.odesim[, 9]
      
      ### discharges
      
      ### new hosp discharges
      hosp.discharges.new.0[k, ] <- tp.k$hosp.discharges.new.odesim[, 1]
      hosp.discharges.new.1[k, ] <- tp.k$hosp.discharges.new.odesim[, 2]
      hosp.discharges.new.2[k, ] <- tp.k$hosp.discharges.new.odesim[, 3]
      hosp.discharges.new.3[k, ] <- tp.k$hosp.discharges.new.odesim[, 4]
      hosp.discharges.new.4[k, ] <- tp.k$hosp.discharges.new.odesim[, 5]
      hosp.discharges.new.5[k, ] <- tp.k$hosp.discharges.new.odesim[, 6]
      hosp.discharges.new.6[k, ] <- tp.k$hosp.discharges.new.odesim[, 7]
      hosp.discharges.new.7[k, ] <- tp.k$hosp.discharges.new.odesim[, 8]
      hosp.discharges.new.8[k, ] <- tp.k$hosp.discharges.new.odesim[, 9]
      
      ### cum hosp discharges 
      hosp.discharges.cum.0[k, ] <- tp.k$hosp.discharges.cum.odesim[, 1]
      hosp.discharges.cum.1[k, ] <- tp.k$hosp.discharges.cum.odesim[, 2]
      hosp.discharges.cum.2[k, ] <- tp.k$hosp.discharges.cum.odesim[, 3]
      hosp.discharges.cum.3[k, ] <- tp.k$hosp.discharges.cum.odesim[, 4]
      hosp.discharges.cum.4[k, ] <- tp.k$hosp.discharges.cum.odesim[, 5]
      hosp.discharges.cum.5[k, ] <- tp.k$hosp.discharges.cum.odesim[, 6]
      hosp.discharges.cum.6[k, ] <- tp.k$hosp.discharges.cum.odesim[, 7]
      hosp.discharges.cum.7[k, ] <- tp.k$hosp.discharges.cum.odesim[, 8]
      hosp.discharges.cum.8[k, ] <- tp.k$hosp.discharges.cum.odesim[, 9] 
    }
  }
    
  ### create lists for age-stratified results
  if (samples$loc == "MA"){
    
    sympt.new.age <- list(
      sympt.new.01,
      sympt.new.2,
      sympt.new.3,
      sympt.new.4,
      sympt.new.5,
      sympt.new.6,
      sympt.new.7,
      sympt.new.8
    )
    
    sympt.cum.age <- list(
      sympt.cum.01,
      sympt.cum.2,
      sympt.cum.3,
      sympt.cum.4,
      sympt.cum.5,
      sympt.cum.6,
      sympt.cum.7,
      sympt.cum.8
    )
    
    hosp.new.age <- list(
      hosp.new.01,
      hosp.new.2,
      hosp.new.3,
      hosp.new.4,
      hosp.new.5,
      hosp.new.6,
      hosp.new.7,
      hosp.new.8
    )
    
    hosp.cum.age <- list(
      hosp.cum.01,
      hosp.cum.2,
      hosp.cum.3,
      hosp.cum.4,
      hosp.cum.5,
      hosp.cum.6,
      hosp.cum.7,
      hosp.cum.8
    )
    
    deaths.new.age <- list(
      deaths.new.01,
      deaths.new.2,
      deaths.new.3,
      deaths.new.4,
      deaths.new.5,
      deaths.new.6,
      deaths.new.7,
      deaths.new.8
    )
        
    deaths.cum.age <- list(
      deaths.cum.01,
      deaths.cum.2,
      deaths.cum.3,
      deaths.cum.4,
      deaths.cum.5,
      deaths.cum.6,
      deaths.cum.7,
      deaths.cum.8
    )
    
    hosp.deaths.new.age <- list(
      hosp.deaths.new.01,
      hosp.deaths.new.2,
      hosp.deaths.new.3,
      hosp.deaths.new.4,
      hosp.deaths.new.5,
      hosp.deaths.new.6,
      hosp.deaths.new.7,
      hosp.deaths.new.8
    )
    
    home.deaths.new.age <- list(
      home.deaths.new.01,
      home.deaths.new.2,
      home.deaths.new.3,
      home.deaths.new.4,
      home.deaths.new.5,
      home.deaths.new.6,
      home.deaths.new.7,
      home.deaths.new.8
    )
    
    hosp.discharges.new.age <- list(
      hosp.discharges.new.01,
      hosp.discharges.new.2,
      hosp.discharges.new.3,
      hosp.discharges.new.4,
      hosp.discharges.new.5,
      hosp.discharges.new.6,
      hosp.discharges.new.7,
      hosp.discharges.new.8
    )
    
    hosp.discharges.cum.age <- list(
      hosp.discharges.cum.01,
      hosp.discharges.cum.2,
      hosp.discharges.cum.3,
      hosp.discharges.cum.4,
      hosp.discharges.cum.5,
      hosp.discharges.cum.6,
      hosp.discharges.cum.7,
      hosp.discharges.cum.8
    )
    
  } else {
    
    sympt.new.age <- list(
      sympt.new.0,
      sympt.new.1,
      sympt.new.2,
      sympt.new.3,
      sympt.new.4,
      sympt.new.5,
      sympt.new.6,
      sympt.new.7,
      sympt.new.8
    )
    
    sympt.cum.age <- list(
      sympt.cum.0,
      sympt.cum.1,
      sympt.cum.2,
      sympt.cum.3,
      sympt.cum.4,
      sympt.cum.5,
      sympt.cum.6,
      sympt.cum.7,
      sympt.cum.8
    )
    
    hosp.new.age <- list(
      hosp.new.0,
      hosp.new.1,
      hosp.new.2,
      hosp.new.3,
      hosp.new.4,
      hosp.new.5,
      hosp.new.6,
      hosp.new.7,
      hosp.new.8
    )
    
    hosp.cum.age <- list(
      hosp.cum.0,
      hosp.cum.1,
      hosp.cum.2,
      hosp.cum.3,
      hosp.cum.4,
      hosp.cum.5,
      hosp.cum.6,
      hosp.cum.7,
      hosp.cum.8
    )
    
    deaths.new.age <- list(
      deaths.new.0,
      deaths.new.1,
      deaths.new.2,
      deaths.new.3,
      deaths.new.4,
      deaths.new.5,
      deaths.new.6,
      deaths.new.7,
      deaths.new.8
    )
    
    deaths.cum.age <- list(
      deaths.cum.0,
      deaths.cum.1,
      deaths.cum.2,
      deaths.cum.3,
      deaths.cum.4,
      deaths.cum.5,
      deaths.cum.6,
      deaths.cum.7,
      deaths.cum.8
    )
    
    hosp.deaths.new.age <- list(
      hosp.deaths.new.0,
      hosp.deaths.new.1,
      hosp.deaths.new.2,
      hosp.deaths.new.3,
      hosp.deaths.new.4,
      hosp.deaths.new.5,
      hosp.deaths.new.6,
      hosp.deaths.new.7,
      hosp.deaths.new.8
    )
    
    home.deaths.new.age <- list(
      home.deaths.new.0,
      home.deaths.new.1,
      home.deaths.new.2,
      home.deaths.new.3,
      home.deaths.new.4,
      home.deaths.new.5,
      home.deaths.new.6,
      home.deaths.new.7,
      home.deaths.new.8
    )
    
    hosp.discharges.new.age <- list(
      hosp.discharges.new.0,
      hosp.discharges.new.1,
      hosp.discharges.new.2,
      hosp.discharges.new.3,
      hosp.discharges.new.4,
      hosp.discharges.new.5,
      hosp.discharges.new.6,
      hosp.discharges.new.7,
      hosp.discharges.new.8
    )
    
    hosp.discharges.cum.age <- list(
      hosp.discharges.cum.0,
      hosp.discharges.cum.1,
      hosp.discharges.cum.2,
      hosp.discharges.cum.3,
      hosp.discharges.cum.4,
      hosp.discharges.cum.5,
      hosp.discharges.cum.6,
      hosp.discharges.cum.7,
      hosp.discharges.cum.8
    ) 
  }
  
  # list to save results
  traj.vals <- list()
  
  # simulated days
  traj.vals$days.odesim <- tp.k$days.odesim
  
  # betas
  traj.vals$beta.full <- beta.full
  
  # rr
  traj.vals$rr.full <- rr.full.m
  
  # total results
  traj.vals$tot.hosp.curr <- tot.hosp.curr
  traj.vals$tot.vent.curr <- tot.vent.curr
  traj.vals$tot.icu.curr <- tot.icu.curr
  traj.vals$tot.hosp.new <- tot.hosp.new
  traj.vals$tot.sympt.new <- tot.sympt.new
  traj.vals$tot.deaths.new <- tot.deaths.new
  traj.vals$tot.sympt.new <- tot.sympt.new
  traj.vals$tot.sympt.cum <- tot.sympt.cum
  traj.vals$tot.hosp.new <- tot.hosp.new
  traj.vals$tot.hosp.cum <- tot.hosp.cum
  traj.vals$tot.deaths.new <- tot.deaths.new
  traj.vals$tot.deaths.cum <- tot.deaths.cum
  traj.vals$tot.hosp.deaths.new <- tot.hosp.deaths.new
  traj.vals$tot.home.deaths.new <- tot.home.deaths.new
  traj.vals$tot.hosp.deaths.cum <- tot.hosp.deaths.cum
  traj.vals$tot.home.deaths.cum <- tot.home.deaths.cum
  traj.vals$tot.hosp.curr.noicu <- tot.hosp.curr.noicu
  traj.vals$tot.icu.curr.novent <- tot.icu.curr.novent
  traj.vals$tot.vent.curr <- tot.vent.curr
  traj.vals$tot.hosp.curr <- tot.hosp.curr
  traj.vals$tot.icu.curr <- tot.icu.curr
  traj.vals$tot.hosp.discharges.new <- tot.hosp.discharges.new
  traj.vals$tot.hosp.discharges.cum <- tot.hosp.discharges.cum
  
  # age stratified results
  traj.vals$sympt.new.age <- sympt.new.age
  traj.vals$sympt.cum.age <- sympt.cum.age
  traj.vals$hosp.new.age <- hosp.new.age
  traj.vals$hosp.cum.age <- hosp.cum.age
  traj.vals$deaths.new.age <- deaths.new.age
  traj.vals$deaths.cum.age <- deaths.cum.age
  traj.vals$hosp.deaths.new.age <- hosp.deaths.new.age 
  traj.vals$home.deaths.new.age <- home.deaths.new.age
  traj.vals$hosp.discharges.new.age <- hosp.discharges.new.age
  traj.vals$hosp.discharges.cum.age <- hosp.discharges.cum.age
  
  # loglikelihood values
  # traj.vals$loglik <- ll.values
  
  # return list with all trajectory values
  return(traj.vals)
}


#################################################################################
### 2. Plotting functions
#################################################################################

### This function adds transparancy to a color. Define transparancy with an 
###   integer between 0 and 255, 0 being fully transparant and 255 being fully 
###   visable. Works with either color and trans a vector of equal length,
###   or one of the two of length 1.
addTrans <- function(
  color,  # vector of colors
  trans   # vector of transparency values
){
  
  if (length(color) != length(trans) & !any(c(length(color), length(trans)) == 1)) {
    stop("Vector lengths not correct")
  }
  
  if (length(color) == 1 & length(trans) > 1) {
    color <- rep(color, length(trans))
  }
  
  if (length(trans) == 1 & length(color) > 1) {
    trans <- rep(trans, length(color))
  }

  num2hex <- function(x) {
    hex <- unlist(strsplit("0123456789ABCDEF", split = ""))
    return(paste(hex[(x - x %% 16) / 16 + 1], hex[x %% 16 + 1], sep = ""))
  }
  
  rgb <- rbind(col2rgb(color), trans)
  res <- paste("#", apply(apply(rgb, 2, num2hex), 2, paste, collapse = ""), sep = "")
  return(res)
}


### Find the proportion of infecteds that enter the ICU, averaged across age classes.
icuTransformation <- function(
  samp,   # icu posterior parameter
  loc     # state used in the analysis
){
 
  if (loc == "RI") {
    pop.frac <- c(0.105, 0.123, 0.140, 0.127, 0.124, 0.135, 0.120, 0.074, 0.052)
  } else if (loc == "MA") {
    pop.frac <- c(
      0.10586466, 0.12243686, 0.14498786, 0.13384234,
      0.12230812, 0.14064248, 0.11801015, 0.06958116, 0.04232637
    )
  } else if (loc == "PA") {
    pop.frac <- c(
      0.11160395, 0.12229803, 0.13156525, 0.12581869,
      0.11809624, 0.13878546, 0.1270166, 0.07657303, 0.04824275
    )
  }

  lewis.prop <- c(
    0.304, 0.293, 0.2825, 0.301, 0.463, 0.4245, 0.46, 0.4835, 0.416
  )
  
  icu.prob <- sapply(samp, function(x) {
    sum(x * pop.frac * lewis.prop)
  })
  
  return(icu.prob)
}


fig2panel <- function(
  beta.directory, 
  ode.directory, 
  data.directory, 
  odepath, 
  loc, 
  const, 
  plot.name, 
  alpha, 
  subsample = NA, 
  axis.size = 1.75, lab.size = 1.75,
  ncol = 2, nrow = 4, 
  height = 12, width = 12, 
  png.true = F,
  font.family = "Helvetica",
  delay = 0
){
  
  
  ### 1. set-up samples of posterior ODESIM trajectories
  
  # grab samples from csv files
  betas <- read.csv(beta.directory)
  ode <- read.csv(ode.directory)
  
  # data
  df <- read.csv(data.directory)
  df <- data.frame(df)
  
  # clean up data and initialize spline basis functions
  if (loc == "PA"){
    
    # ode parameter output
    ode.inds <- 1:34
    
    # removing days before 61
    idx.remove <- which(df$daynum<61)
    if(length(idx.remove) > 0){
      df <- df[-idx.remove,]
    }
    
    # removing days after 250
    idx.remove <- which(df$daynum>250)
    if(length(idx.remove) > 0){
      df <- df[-idx.remove,]
    }
    
    # remove NAs
    df <- df[!is.na(df$daynum), ]
    
    # default values
    n.days <- nrow(df)
    days <- df$daynum
    end.day <- max(days, na.rm = TRUE) + 1
    num.days <- end.day - 60
    
    # cubic spline with one basis function every 7 days
    bspl <- create.bspline.basis(c(61, end.day), nbasis = round(num.days / 7))
    Z <- eval.basis(bspl, 61:end.day)
    
    # iSpline for reporting rate
    knots <- c(84, 92, 100, 108, 130, 160, 190)
    Z.rr <- iSpline(x = 61:end.day, knots = knots, degree = 2)
    Z.rr <- cbind(rep(1, nrow(Z.rr)), Z.rr)
    Z.rr <- Z.rr[, c(1, 4:7)]
    
    # introday
    id = 60
    
  } else if (loc == "RI"){
    
    # ode parameter output
    ode.inds <- 1:35
    
    # removing days before 61
    idx.remove <- which(df$daynum<61)
    if (length(idx.remove) > 0) {
      df = df[-idx.remove, ]
    }
    
    # default values
    n.days <- nrow(df)
    days <- df$daynum
    end.day <- max(days, na.rm = TRUE) + 1
    num.days <- end.day - 60 
    
    # cubic spline with one basis function every 7 days
    bspl <- create.bspline.basis(c(61, end.day), nbasis = round(num.days / 7))
    Z <- eval.basis(bspl, 61:end.day)
    
    # iSpline for reporting rate
    knots <- c(84, 92, 100, 108, 130, 160, 190)
    Z.rr <- iSpline(x = 61:end.day, knots = knots, degree = 2)
    Z.rr <- cbind(rep(1, nrow(Z.rr)), Z.rr)
    Z.rr <- Z.rr[, c(1, 4:7)]
    
    # introday
    id <- 55
    
  } else if (loc == "MA") {
    
    # ode parameter output
    ode.inds <- 1:35
    
    # removing days before 61
    idx.remove <- which(df$daynum<61)
    if(length(idx.remove)>0){
      df <- df[-idx.remove,]
    }
    
    # default values
    n.days <- nrow(df)
    days <- df$daynum
    end.day <- max(days, na.rm = TRUE) + 1
    num.days <- end.day - 60 
    
    ### Create spline expansions

    # Cubic spline with one basis function every 7 days
    bspl <- create.bspline.basis(c(61, end.day), nbasis = round(num.days / 7))
    Z <- eval.basis(bspl, 61:end.day)
    
    # zero-order spline with one rate
    bspl.rr <- create.bspline.basis(c(61, end.day), nbasis = 1, norder = 1)
    Z.rr <- eval.basis(bspl.rr, 61:end.day)
    Z.rr <- Z.rr[-nrow(Z.rr), ]
    
    # introday
    id <- 55
  }
  
  # list with required structures
  samples <- list(
    betas = as.matrix(betas),
    ode = as.matrix(ode[, ode.inds]),
    rr = as.matrix(ode[, -ode.inds]),
    days = days,
    Z.beta = Z,
    Z.rr = Z.rr,
    df = df,
    loc = loc,
    const = const,
    introday = id,
    pres.delay = delay,
    sf.choice = FALSE
  )
  
  # if subsampling, remove some samples for plotting
  if (!is.na(subsample)) {
    # remove some of the samples
    n.orig <- nrow(samples$betas)
    new.samples <- sample.int(n.orig, subsample, replace = FALSE)
    samples$betas <- samples$betas[new.samples, ]
    samples$ode <- samples$ode[new.samples, ]
    samples$rr <- samples$rr[new.samples, ]
  }
  
  # reformat column names
  colnames(samples$ode) <- gsub("\\.", "-", colnames(samples$ode))
  colnames(samples$ode) <- gsub(
    "hosp-report-rate",
    "hosp.report.rate",
    colnames(samples$ode)
  )
  
  # process trajectories from the samples
  tp <- traj.sim(
    samples = samples,
    odepath = odepath,
    csv = T
  )
  
  # process actual data
  dp <- data.fine.process(samples$df, loc = samples$loc)
  
  # number of simulated output
  S <- nrow(tp$tot.sympt.new)

  # vector of dates
  dates <- seq(as.Date("2020/1/1"), by = "day", length.out = 365)

  # last day of simulated output
  end.day <- max(samples$days)

  # save plots as pdf
  if (png.true){
    png(filename = plot.name,
        family = font.family, 
        width = width, 
        height = height, 
        units = "in",
        res = 100)
  } else {
    pdf(file = plot.name, family = font.family, h = height, w = width)
  }
  # grid of plots
  par(mfrow = c(nrow,ncol))
  par(mar = c(3, 2, 1, 1), oma = c(2, 2, 0.5, 0.5))

  ##########################################
  ### A. new symptomatic cases
  ##########################################

  y.max <- max(
    max(dp$tot.sympt.new, na.rm = T),
    max(tp$tot.sympt.new)
  )
  
  # plot total new symptomatic cases
  plot(dates[tp$days.odesim],
    tp$tot.sympt.new[1, ],
    type = "l", col = addTrans("grey", alpha),
    ylim = c(0, y.max), ylab = "", xlab = "",
    cex.axis = axis.size, cex.lab = lab.size, bty = "l"
  )
  title("New Symptomatic Cases (Reported)",
    adj = 0.95,
    line = -1.5, cex.main = 1.75
  )
  title("A", adj = 0.05, line = -1.5, cex.main = 2)

  for (j in 2:S) {
    lines(dates[tp$days.odesim],
      tp$tot.sympt.new[j, ],
      type = "l",
      col = addTrans("grey", alpha)
    )
  }
  
  # add median and 95% quantiles
  mn <- apply(tp$tot.sympt.new, 2, median)
  uq <- apply(tp$tot.sympt.new, 2, quantile, 0.975)
  lq <- apply(tp$tot.sympt.new, 2, quantile, 0.025)

  lines(dates[tp$days.odesim], mn, type = "l", col = "blue")
  lines(dates[tp$days.odesim], uq, type = "l", lty = 3, col = "blue")
  lines(dates[tp$days.odesim], lq, type = "l", lty = 3, col = "blue")

  # add observed data
  points(dates[dp$days], dp$tot.sympt.new, pch = 19, col = "black", cex = 0.5)


  #################################################
  ### B. New Hospitalizations
  #################################################

  y.max <- max(
    max(dp$tot.hosp.new, na.rm = T),
    max(tp$tot.hosp.new)
  )

  plot(dates[tp$days.odesim],
    tp$tot.hosp.new[1, ],
    type = "l", col = addTrans("grey", alpha),
    ylim = c(0, y.max), ylab = "", xlab = "",
    cex.axis = axis.size, cex.lab = lab.size, bty = "l"
  )

  title("New Hospitalizations",
    adj = 0.95,
    line = -1.5, cex.main = 1.75
  )
  
  title("B", adj = 0.05, line = -1.5, cex.main = 2)

  for (j in 2:S) {
    lines(dates[tp$days.odesim],
      tp$tot.hosp.new[j, ],
      type = "l",
      col = addTrans("grey", alpha)
    )
  }

  # add median and 95% quantiles
  mn <- apply(tp$tot.hosp.new, 2, median)
  uq <- apply(tp$tot.hosp.new, 2, quantile, 0.975)
  lq <- apply(tp$tot.hosp.new, 2, quantile, 0.025)

  lines(dates[tp$days.odesim], mn, type = "l", col = "blue")
  lines(dates[tp$days.odesim], uq, type = "l", lty = 3, col = "blue")
  lines(dates[tp$days.odesim], lq, type = "l", lty = 3, col = "blue")

  # add observed data (RI only)
  if (samples$loc == "RI") {
    points(dates[dp$days], dp$tot.hosp.new, pch = 19, col = "black", cex = 0.5)
  }
  

  ###################################
  ### C. Current Hospitalizations ###
  ###################################

  y.max <- max(
    max(dp$tot.hosp.curr, na.rm = T),
    max(tp$tot.hosp.curr)
  )
  
  # plot current hospitalizations
  plot(dates[tp$days.odesim],
    tp$tot.hosp.curr[1, ],
    type = "l", col = addTrans("grey", alpha),
    ylim = c(0, y.max), ylab = "", xlab = "",
    cex.axis = axis.size, cex.lab = lab.size, bty = "l"
  )
  title("Current Hospitalizations",
    adj = 0.95,
    line = -1.5, cex.main = 1.75
  )
  title("C", adj = 0.05, line = -1.5, cex.main = 2)

  for (j in 2:S) {
    lines(dates[tp$days.odesim],
      tp$tot.hosp.curr[j, ],
      type = "l",
      col = addTrans("grey", alpha)
    )
  }

  # add median and 95% quantiles
  mn <- apply(tp$tot.hosp.curr, 2, median)
  uq <- apply(tp$tot.hosp.curr, 2, quantile, 0.975)
  lq <- apply(tp$tot.hosp.curr, 2, quantile, 0.025)

  lines(dates[tp$days.odesim], mn, type = "l", col = "blue")
  lines(dates[tp$days.odesim], uq, type = "l", lty = 3, col = "blue")
  lines(dates[tp$days.odesim], lq, type = "l", lty = 3, col = "blue")

  # add observed data
  points(dates[dp$days], dp$tot.hosp.curr, pch = 19, col = "black", cex = 0.5)


  ###################################
  ### D. Current ICU Occupancy    ###
  ###################################

  y.max <- max(
    max(dp$tot.icu.curr, na.rm = T),
    max(tp$tot.icu.curr)
  )
  
  # plot current ICU
  plot(dates[tp$days.odesim],
    tp$tot.icu.curr[1, ],
    type = "l", col = addTrans("grey", alpha),
    ylim = c(0, y.max), ylab = "", xlab = "", yaxt = "none",
    cex.axis = axis.size, cex.lab = lab.size, bty = "l"
  )
  axis(2, seq(0, 120, 20),
    cex.axis = axis.size, cex.lab = lab.size,
    labels = c("0", "", "40", "", "80", "", "120")
  )
  title("Current ICU Occupancy",
    adj = 0.95,
    line = -1.5, cex.main = 1.75
  )
  title("D", adj = 0.05, line = -1.5, cex.main = 2)
  
  for (j in 2:S) {
    lines(dates[tp$days.odesim],
      tp$tot.icu.curr[j, ],
      type = "l",
      col = addTrans("grey", alpha)
    )
  }
  
  # add median and 95% quantiles
  mn <- apply(tp$tot.icu.curr, 2, median)
  uq <- apply(tp$tot.icu.curr, 2, quantile, 0.975)
  lq <- apply(tp$tot.icu.curr, 2, quantile, 0.025)

  lines(dates[tp$days.odesim], mn, type = "l", col = "blue")
  lines(dates[tp$days.odesim], uq, type = "l", lty = 3, col = "blue")
  lines(dates[tp$days.odesim], lq, type = "l", lty = 3, col = "blue")

  # observed data
  points(dates[dp$days], dp$tot.icu.curr, pch = 19, col = "black", cex = 0.5)


  #######################################
  ### E. Current Ventilator Occupancy ###
  #######################################

  y.max <- max(
    max(dp$tot.vent.curr, na.rm = T),
    max(tp$tot.vent.curr)
  )
  
  # plot current vent
  plot(dates[tp$days.odesim],
    tp$tot.vent.curr[1, ],
    type = "l", col = addTrans("grey", alpha),
    ylim = c(0, y.max), ylab = "", xlab = "",
    cex.axis = axis.size, cex.lab = lab.size, bty = "l"
  )
  title("Current Ventilator Occupancy",
    adj = 0.95, line = -1.5, cex.main = 1.75
  )
  title("E", adj = 0.05, line = -1.5, cex.main = 2)

  for (j in 2:S) {
    lines(dates[tp$days.odesim],
      tp$tot.vent.curr[j, ],
      type = "l",
      col = addTrans("grey", alpha)
    )
  }
  
  # add median and 95% quantiles
  mn <- apply(tp$tot.vent.curr, 2, median)
  uq <- apply(tp$tot.vent.curr, 2, quantile, 0.975)
  lq <- apply(tp$tot.vent.curr, 2, quantile, 0.025)

  lines(dates[tp$days.odesim], mn, type = "l", col = "blue")
  lines(dates[tp$days.odesim], uq, type = "l", lty = 3, col = "blue")
  lines(dates[tp$days.odesim], lq, type = "l", lty = 3, col = "blue")

  # observed data
  points(dates[dp$days], dp$tot.vent.curr, pch = 19, col = "black", cex = 0.5)


  ###############################################
  ### F. Deaths
  ###############################################

  y.max <- max(
    max(dp$tot.deaths.new, na.rm = T),
    max(tp$tot.deaths.new)
  )
  
  # plot new deaths
  plot(dates[tp$days.odesim],
    tp$tot.deaths.new[1, ],
    type = "l", col = addTrans("grey", alpha),
    ylim = c(0, y.max), ylab = "", xlab = "",
    cex.axis = axis.size, cex.lab = lab.size, bty = "l"
  )
  title("New Deaths", adj = 0.95, line = -1.5, cex.main = 1.75)
  title("F", adj = 0.05, line = -1.5, cex.main = 2)

  for (j in 2:S) {
    lines(dates[tp$days.odesim],
      tp$tot.deaths.new[j, ],
      type = "l",
      col = addTrans("grey", alpha)
    )
  }
  
  # add median and 95% quantiles
  mn <- apply(tp$tot.deaths.new, 2, median)
  uq <- apply(tp$tot.deaths.new, 2, quantile, 0.975)
  lq <- apply(tp$tot.deaths.new, 2, quantile, 0.025)

  lines(dates[tp$days.odesim], mn, type = "l", col = "blue")
  lines(dates[tp$days.odesim], uq, type = "l", lty = 3, col = "blue")
  lines(dates[tp$days.odesim], lq, type = "l", lty = 3, col = "blue")

  # observed data
  points(dates[dp$days], dp$tot.deaths.new, pch = 19, col = "black", cex = 0.5)

  ###############################################
  ### G. Hospital Deaths
  ###############################################

  if (samples$loc == "MA") {
    y.max <- max(tp$tot.hosp.deaths.new)
  } else {
    y.max <- max(
      max(dp$tot.hosp.deaths.new, na.rm = T),
      max(tp$tot.hosp.deaths.new)
    )
  }

  # plot new deaths
  plot(dates[tp$days.odesim],
    tp$tot.hosp.deaths.new[1, ],
    type = "l", col = addTrans("grey", alpha),
    ylim = c(0, y.max), ylab = "", xlab = "",
    cex.axis = axis.size, cex.lab = lab.size, bty = "l"
  )
  title("New Hospital Deaths",
    adj = 0.95,
    line = -1.5, cex.main = 1.75
  )
  title("G", adj = 0.05, line = -1.5, cex.main = 2)

  for (j in 2:S) {
    lines(dates[tp$days.odesim],
      tp$tot.hosp.deaths.new[j, ],
      type = "l",
      col = addTrans("grey", alpha)
    )
  }
  
  # add median and 95% quantiles
  mn <- apply(tp$tot.hosp.deaths.new, 2, median)
  uq <- apply(tp$tot.hosp.deaths.new, 2, quantile, 0.975)
  lq <- apply(tp$tot.hosp.deaths.new, 2, quantile, 0.025)

  lines(dates[tp$days.odesim], mn, type = "l", col = "blue")
  lines(dates[tp$days.odesim], uq, type = "l", lty = 3, col = "blue")
  lines(dates[tp$days.odesim], lq, type = "l", lty = 3, col = "blue")

  # observed data
  points(dates[dp$days], dp$tot.hosp.deaths.new, pch = 19, col = "black", cex = 0.5)


  ###############################################
  ### H. Discharges
  ###############################################

  # new home deaths
  if (samples$loc == "MA") {
    y.max <- max(tp$tot.hosp.discharges.new)
  } else {
    y.max <- max(
      max(dp$tot.hosp.discharges.new, na.rm = T),
      max(tp$tot.hosp.discharges.new)
    )
  }
  
  # plot new deaths
  plot(dates[tp$days.odesim],
    tp$tot.hosp.discharges.new[1, ],
    type = "l", col = addTrans("grey", alpha),
    ylim = c(0, y.max), ylab = "", xlab = "",
    cex.axis = axis.size, cex.lab = lab.size, bty = "l"
  )
  title("New Hospital Discharges",
    adj = 0.95,
    line = -1.5, cex.main = 1.75
  )
  title("H", adj = 0.05, line = -1.5, cex.main = 2)

  for (j in 2:S) {
    lines(dates[tp$days.odesim],
      tp$tot.hosp.discharges.new[j, ],
      type = "l",
      col = addTrans("grey", alpha)
    )
  }
  
  # add median and 95% quantiles
  mn <- apply(tp$tot.hosp.discharges.new, 2, median)
  uq <- apply(tp$tot.hosp.discharges.new, 2, quantile, 0.975)
  lq <- apply(tp$tot.hosp.discharges.new, 2, quantile, 0.025)

  lines(dates[tp$days.odesim], mn, type = "l", col = "blue")
  lines(dates[tp$days.odesim], uq, type = "l", lty = 3, col = "blue")
  lines(dates[tp$days.odesim], lq, type = "l", lty = 3, col = "blue")

  # observed data (not MA)
  if (samples$loc != "MA") {
    points(dates[dp$days], dp$tot.hosp.discharges.new, pch = 19, col = "black", cex = 0.5)
  }
  
  dev.off()
  embed_fonts(plot.name, outfile = plot.name)
}




fig5panel <- function(
  ## Rhode Island directories ##
  ri.beta, ri.ode, ri.data, ri.const,
  ## Pennsylvania directories ##
  pa.beta, pa.ode, pa.data, pa.const,
  ## Massachusetts directories ##
  ma.beta, ma.ode, ma.data, ma.const,
  ## ODESIM directory ##
  odepath, 
  ## plotting inputs ##
  plot.name, alpha, subsample = NA, 
  axis.size = 1, lab.size = 1,
  height = 16, width = 11, logplot = F,
  font.family = "Helvetica"
){
  
  ### make data structure for traj.sim
  
  ## i. Rhode Island ##
  
  # grab samples from csv files
  ri.betas <- read.csv(ri.beta)
  ri.ode <- read.csv(ri.ode)
  
  # data
  ri.df <- read.csv(ri.data)
  ri.df <- data.frame(ri.df)
  
  # ode index
  ri.ode.inds <- 1:35
  
  ## removing days before 61
  ri.idx.remove <- which(ri.df$daynum < 61)
  if (length(ri.idx.remove) > 0) {
    ri.df <- ri.df[-ri.idx.remove, ]
  }
  
  # default settings
  ri.n.days <- nrow(ri.df)
  ri.days <- ri.df$daynum
  ri.end.day <- max(ri.days, na.rm = TRUE) + 1
  ri.num.days <- ri.end.day - 60 ## number of days with data
  
  # cubic spline with one basis function every 7 days
  ri.bspl <- create.bspline.basis(c(61, ri.end.day), nbasis = round(ri.num.days / 7))
  ri.Z <- eval.basis(ri.bspl, 61:ri.end.day)
  
  # iSpline for reporting rate
  ri.knots <- c(84, 92, 100, 108, 130, 160, 190)
  ri.Z.rr <- iSpline(x = 61:ri.end.day, knots = ri.knots, degree = 2)
  ri.Z.rr <- cbind(rep(1, nrow(ri.Z.rr)), ri.Z.rr)
  ri.Z.rr <- ri.Z.rr[-nrow(ri.Z.rr), c(1, 4:7)]
  
  # introday
  ri.id <- 55
  
  # list with required structures
  ri.samples <- list(
    betas = as.matrix(ri.betas),
    ode = as.matrix(ri.ode[, ri.ode.inds]),
    rr = as.matrix(ri.ode[, -ri.ode.inds]),
    days = ri.days,
    Z.beta = ri.Z,
    Z.rr = ri.Z.rr,
    df = ri.df,
    loc = "RI",
    const = ri.const,
    introday = ri.id,
    pres.delay = 2,
    sf.choice = FALSE
  )
  
  # reformat column names
  colnames(ri.samples$ode) <- gsub('\\.', '-', colnames(ri.samples$ode))
  colnames(ri.samples$ode) <- gsub(
    "hosp-report-rate",
    "hosp.report.rate",
    colnames(ri.samples$ode)
  )
  
  # process trajectories from the samples
  ri.tp <- traj.sim(
    samples = ri.samples,
    odepath = odepath,
    csv = TRUE
  )
  
  # process actual data
  ri.dp <- data.fine.process(ri.samples$df, loc = ri.samples$loc)

  ## ii. Massachusetts ##
  
  # grab samples from csv files
  ma.betas <- read.csv(ma.beta)
  ma.ode <- read.csv(ma.ode)
  
  # data
  ma.df <- read.csv(ma.data)
  ma.df <- data.frame(ma.df)
  
  # ode index
  ma.ode.inds <- 1:35
  
  # removing days before 61
  ma.idx.remove <- which(ma.df$daynum < 61)
  if (length(ma.idx.remove) > 0) {
    ma.df <- ma.df[-ma.idx.remove, ]
  }
  
 # default settings
  ma.n.days <- nrow(ma.df)
  ma.days <- ma.df$daynum
  ma.end.day <- max(ma.days, na.rm = TRUE) + 1
  ma.num.days <- ma.end.day - 60 ## number of days with data
  
  # cubic spline with one basis function every 7 days
  ma.bspl <- create.bspline.basis(c(61, ma.end.day),
    nbasis = round(ma.num.days / 7)
  )
  ma.Z <- eval.basis(ma.bspl, 61:ma.end.day)
  
  # zero-order spline with one rate
  ma.bspl.rr <- create.bspline.basis(c(61, ma.end.day), nbasis = 1, norder = 1)
  ma.Z.rr <- eval.basis(ma.bspl.rr, 61:ma.end.day)
  ma.Z.rr <- ma.Z.rr[-nrow(ma.Z.rr),]
  ma.Z.rr <- as.vector(ma.Z.rr)
  
  # introday
  ma.id = 55
  
  # list with required structures
  ma.samples <- list(
    betas = as.matrix(ma.betas),
    ode = as.matrix(ma.ode[, ma.ode.inds]),
    rr = as.matrix(ma.ode[, -ma.ode.inds]),
    days = ma.days,
    Z.beta = ma.Z,
    Z.rr = ma.Z.rr,
    df = ma.df,
    loc = "MA",
    const = ma.const,
    pres.delay = 2,
    introday = ma.id,
    sf.choice = FALSE
  )
  
  # reformat column names
  colnames(ma.samples$ode) <- gsub('\\.', '-', colnames(ma.samples$ode))
  colnames(ma.samples$ode) <- gsub(
    "hosp-report-rate",
    "hosp.report.rate",
    colnames(ma.samples$ode)
  )
  
  # process trajectories from the samples
  ma.tp <- traj.sim(
    samples = ma.samples,
    odepath = odepath,
    csv = TRUE
  )
  
  # process actual data
  ma.dp <- data.fine.process(ma.samples$df, loc = ma.samples$loc)
  
  ## iii. Pennsylvania ##
  
  # ode index
  pa.ode.inds <- 1:34
  
  # grab samples from csv files
  pa.betas <- read.csv(pa.beta)
  pa.ode <- read.csv(pa.ode)
  
  # data
  pa.df <- read.csv(pa.data)
  pa.df <- data.frame(pa.df)
  
  ## removing days before 61
  pa.idx.remove <- which(pa.df$daynum < 61)
  if (length(pa.idx.remove) > 0) {
    pa.df <- pa.df[-pa.idx.remove, ]
  }
  
  # default settings
  pa.n.days <- nrow(pa.df)
  pa.days <- pa.df$daynum
  pa.end.day <- max(pa.days, na.rm = TRUE) + 1
  pa.num.days <- pa.end.day - 60 
  
  # cubic spline with one basis function every 7 days
  pa.bspl <- create.bspline.basis(c(61, pa.end.day), nbasis=round(pa.num.days / 7))
  pa.Z <- eval.basis(pa.bspl, 61:pa.end.day)
  
  # iSpline for reporting rate
  pa.knots <- c(84, 92, 100, 108, 130, 160, 190)
  pa.Z.rr <- iSpline(x = 61:pa.end.day, knots = pa.knots, degree = 2)
  pa.Z.rr <- cbind(rep(1, nrow(pa.Z.rr)), pa.Z.rr)
  pa.Z.rr <- pa.Z.rr[-nrow(pa.Z.rr), c(1, 4:7)]
  
  # introday
  pa.id <- 60

  # list with required structures
  pa.samples <- list(
    betas = as.matrix(pa.betas),
    ode = as.matrix(pa.ode[, pa.ode.inds]),
    rr = as.matrix(pa.ode[, -pa.ode.inds]),
    days = pa.days,
    Z.beta = pa.Z,
    Z.rr = pa.Z.rr,
    df = pa.df,
    loc = "PA",
    const = pa.const,
    pres.delay = 2,
    introday = pa.id,
    sf.choice = FALSE
  )
  
  # reformat column names
  colnames(pa.samples$ode) <- gsub('\\.', '-', colnames(pa.samples$ode))
  colnames(pa.samples$ode) <- gsub(
    "hosp-report-rate",
    "hosp.report.rate",
    colnames(pa.samples$ode)
  )
  
  # process trajectories from the samples
  pa.tp <- traj.sim(
    samples = pa.samples,
    odepath = odepath,
    csv = TRUE
  )
  
  # process actual data
  pa.dp <- data.fine.process(pa.samples$df, loc = pa.samples$loc)
  
  ri.color = "#733381"
  ma.color = "#DD571C"
  pa.color = "#589C48"
  
  ############################################################################
  ### PLOT FIGURE FIVE
  ############################################################################

  pdf(file = plot.name, family = font.family, h = height, w = width)
  
  ### left, right, bottom, top
  
  m0 <- matrix(c(0, 1, 0, 0.97), nrow = 1, ncol = 4, byrow = TRUE)
  
  m1 <- rbind(
    c(0.01, 0.99, 0.84, 0.99),
    c(0.01, 0.99, 0.67, 0.82),
    c(0.01, 0.99, 0.50, 0.65),
    c(0.01, 0.99, 0.33, 0.48),
    c(0.01, 0.99, 0.01, 0.31)
  )
  
  m2 <- rbind(
    c(0.01, 0.33, 0.02, 0.98),
    c(0.33, 0.66, 0.02, 0.98),
    c(0.66, 1, 0.02, 0.98)
  )
  
  m3 <- matrix(c(0.01, 0.99, 0.01, 0.99), nrow = 1, ncol = 4, byrow = TRUE)
  
  
  split.screen(m0)
  split.screen(m1, screen = 1)
  split.screen(m2, screen = 2)
  split.screen(m2, screen = 3)
  split.screen(m2, screen = 4)
  split.screen(m2, screen = 5)
  split.screen(m3, screen = 6)
  
  
  screen(1)
  par(
    mar = c(0, 0, 0, 0),
    oma = c(0, 0, 2, 0)
  )
  
  title("Rhode Island",
    adj = 0.15, line = 0, cex.main = 1.55, outer = TRUE
  )
  
  title("Massachusetts",
    adj = 0.53, line = 0, cex.main = 1.55, outer = TRUE
  )
  
  title("Pennsylvania",
    adj = 0.9, line = 0, cex.main = 1.55, outer = TRUE
  )
  
  ############################################################################
  ### A. REPORTING RATE 
  ############################################################################
  
  screen(2)
  # bottom, left, top, right
  par(
    mar = c(0, 0, 2, 0),
    oma = c(0, 0, 0, 0)
  )
  
  title("A. Symptomatic Reporting Rate Over Time",
    adj = 0, line = 1, cex.main = 1.25, outer = FALSE
  )
  
  #############################
  ###  A1: RI Reporting Rate
  #############################
  screen(7)
  
  # bottom, left, top, right
  par(
    mar = c(1, 2, 2, 0),
    oma = c(0, 0, 0, 0)
  )
  
  # some constants
  n.orig <- nrow(ri.samples$betas)
  new.samples <- sample.int(n.orig, 150, replace = FALSE)
  ri.rr <- ri.tp$rr.full[new.samples, ]
  
  # vector of dates
  dates <- seq(as.Date("2020/1/1"), by = "day", length.out = 365)
  
  # last day of simulated output
  ri.end.day <- max(ri.samples$days)
  
  # plot reporting rate over time
  plot(dates[c(ri.end.day - ncol(ri.rr) + 1):ri.end.day],
    ri.rr[1, ],
    type = "l", col = addTrans("grey", alpha),
    ylim = c(0, 1), ylab = "", xlab = "",
    cex.axis = axis.size, cex.lab = lab.size, bty = "l"
  )
  
  for (j in 2:nrow(ri.rr)) {
    lines(dates[c(ri.end.day - ncol(ri.rr) + 1):ri.end.day],
      ri.rr[j, ],
      type = "l",
      col = addTrans("grey", alpha)
    )
  }
  
  # add median and 95% quantiles
  mn <- apply(ri.rr, 2, median)
  uq <- apply(ri.rr, 2, quantile, 0.975)
  lq <- apply(ri.rr, 2, quantile, 0.025)
  
  lines(dates[c(ri.end.day - ncol(ri.rr) + 1):ri.end.day],
    mn, type = "l", col = ri.color, lwd = 2
  )
  lines(dates[c(ri.end.day - ncol(ri.rr) + 1):ri.end.day],
    uq, type = "l", lty = 3, col = ri.color, lwd = 2
  )
  lines(dates[c(ri.end.day - ncol(ri.rr) + 1):ri.end.day],
    lq, type = "l", lty = 3, col = ri.color, lwd = 2
  )
  
  #############################
  ###  A2: MA Reporting Rate
  #############################
  
  screen(8)
  # bottom, left, top, right
  par(
    mar = c(1, 2, 2, 0),
    oma = c(0, 0, 0, 0)
  )
  
  # vector of dates
  dates <- seq(as.Date("2020/1/1"), by = "day", length.out = 365)

  ma.rr <- ma.tp$rr.full[new.samples, ]
  
  # last day of simulated output
  ma.end.day <- max(ma.samples$days)
  
  
  # plot reporting rate over time
  plot(dates[c(ma.end.day - ncol(ma.rr) + 1):ma.end.day],
    ma.rr[1, ],
    type = "l", col = addTrans("grey", alpha),
    ylim = c(0, 1), ylab = "", xlab = "",
    cex.axis = axis.size, cex.lab = lab.size, bty = "l"
  )
  
  for (j in 2:nrow(ma.rr)) {
    lines(dates[c(ma.end.day - ncol(ma.rr) + 1):ma.end.day],
      ma.rr[j, ],
      type = "l",
      col = addTrans("grey", alpha)
    )
  }
  
  # add median and 95% quantiles
  mn <- apply(ma.rr, 2, median)
  uq <- apply(ma.rr, 2, quantile, 0.975)
  lq <- apply(ma.rr, 2, quantile, 0.025)
  
  lines(dates[c(ma.end.day - ncol(ma.rr) + 1):ma.end.day],
    mn,
    type = "l", col = ma.color, lwd = 2
  )
  lines(dates[c(ma.end.day - ncol(ma.rr) + 1):ma.end.day],
    uq,
    type = "l", lty = 3, col = ma.color, lwd = 2
  )
  lines(dates[c(ma.end.day - ncol(ma.rr) + 1):ma.end.day],
    lq,
    type = "l", lty = 3, col = ma.color, lwd = 2
  )
  
  #############################
  ###  A3: PA Reporting Rate
  #############################
  screen(9)
  # bottom, left, top, right
  par(
    mar = c(1, 2, 2, 0),
    oma = c(0, 0, 0, 0)
  )
  
  # vector of dates
  dates <- seq(as.Date("2020/1/1"), by = "day", length.out = 365)
  
  pa.rr <- pa.tp$rr.full[new.samples, ]
  
  # last day of simulated output
  pa.end.day <- max(pa.samples$days)
  
  # plot reporting rate over time
  plot(dates[c(pa.end.day - ncol(pa.rr) + 1):pa.end.day],
    pa.rr[1, ],
    type = "l", col = addTrans("grey", alpha),
    ylim = c(0, 1), ylab = "", xlab = "",
    cex.axis = axis.size, cex.lab = lab.size, bty = "l"
  )
  
  for (j in 2:nrow(pa.rr)) {
    lines(dates[c(pa.end.day - ncol(pa.rr) + 1):pa.end.day],
      pa.rr[j, ],
      type = "l",
      col = addTrans("grey", alpha)
    )
  }
  
  # add median and 95% quantiles
  mn <- apply(pa.rr, 2, median)
  uq <- apply(pa.rr, 2, quantile, 0.975)
  lq <- apply(pa.rr, 2, quantile, 0.025)
  
  lines(dates[c(pa.end.day - ncol(pa.rr) + 1):pa.end.day],
    mn,
    type = "l", col = pa.color, lwd = 2
  )
  lines(dates[c(pa.end.day - ncol(pa.rr) + 1):pa.end.day],
    uq,
    type = "l", lty = 3, col = pa.color, lwd = 2
  )
  lines(dates[c(pa.end.day - ncol(pa.rr) + 1):pa.end.day],
    lq,
    type = "l", lty = 3, col = pa.color, lwd = 2
  )
  
  ############################################################################
  ### B. HOSPITAL STAY
  ############################################################################
  
  screen(3)
  # bottom, left, top, right
  par(
    mar = c(0, 0, 2, 0),
    oma = c(0, 0, 0, 0)
  )
  
  title("B. Length of Hospital Stay (Days)",
    adj = 0, line = 1, cex.main = 1.25, outer = FALSE
  )
  
  # densities for each state
  ri.hospstay.ind <- which(colnames(ri.samples$ode) == "dev-len-hospstay")
  ri.hospstay.dens <- density(10 * ri.samples$ode[, ri.hospstay.ind])
  
  pa.hospstay.ind <- which(colnames(pa.samples$ode) == "dev-len-hospstay")
  pa.hospstay.dens <- density(10 * pa.samples$ode[, pa.hospstay.ind])
  
  ma.hospstay.ind <- which(colnames(ma.samples$ode) == "dev-len-hospstay")
  ma.hospstay.dens <- density(10 * ma.samples$ode[, ma.hospstay.ind])

  
  # plot limits
  max.x <- max(
    max(ri.hospstay.dens$x),
    max(pa.hospstay.dens$x),
    max(ma.hospstay.dens$x)
  )
  
  min.x <- min(
    min(ri.hospstay.dens$x),
    min(pa.hospstay.dens$x),
    min(ma.hospstay.dens$x)
  )
  
  #############################
  ###  B1: RI Hospital Stay
  #############################
  screen(10)
  # bottom, left, top, right
  par(
    mar = c(1, 2, 2, 0),
    oma = c(0, 0, 0, 0)
  )
  
  plot(ri.hospstay.dens,
    xlab = "", ylab = "", main = "", bty = "l",
    xlim = c(min.x, max.x), cex.axis = axis.size, cex.lab = lab.size
  )
  
  polygon(ri.hospstay.dens,
    col = addTrans(ri.color, 200),
    border = ri.color,
    lwd = 2
  )
  
  #############################
  ###  B2: MA Hospital Stay
  #############################
  screen(11)
  # bottom, left, top, right
  par(
    mar = c(1, 2, 2, 0),
    oma = c(0, 0, 0, 0)
  )
  
  plot(ma.hospstay.dens,
    xlab = "", ylab = "", main = "", bty = "l",
    xlim = c(min.x, max.x),
    cex.axis = axis.size, cex.lab = lab.size,
    col = addTrans(ma.color, 200)
  )
  
  polygon(ma.hospstay.dens,
    col = addTrans(ma.color, 200),
    border = addTrans(ma.color, 200),
    lwd = 2
  )
  
  #############################
  ###  B3: PA Hospital Stay
  #############################
  screen(12)
  # bottom, left, top, right
  par(
    mar = c(1, 2, 2, 0),
    oma = c(0, 0, 0, 0)
  )
  
  plot(pa.hospstay.dens,
    xlab = "", ylab = "", main = "", bty = "l",
    xlim = c(min.x, max.x),
    cex.axis = axis.size, cex.lab = lab.size,
    col = addTrans(pa.color, 200)
  )
  
  polygon(pa.hospstay.dens,
    col = addTrans(pa.color, 200),
    border = addTrans(pa.color, 200),
    lwd = 2
  )
  
  ############################################################################
  ### C. PROBABILITY OF DYING AT HOME
  ############################################################################
  
  screen(4)
  # bottom, left, top, right
  par(
    mar = c(0, 0, 2, 0),
    oma = c(0, 0, 0, 0)
  )
  
  title("C. Probability of Death Outside Hospital Setting",
    adj = 0, line = 1, cex.main = 1.25, outer = FALSE
  )
  
  # home 60 (RI only)
  ri.home60.ind <- which(colnames(ri.samples$ode) == "death-prob-home-60")
  ri.home60.dens <- density(ri.samples$ode[, ri.home60.ind])
  
  # home 70
  ri.home70.ind <- which(colnames(ri.samples$ode) == "death-prob-home-70")
  ri.home70.dens <- density(ri.samples$ode[, ri.home70.ind])
  
  ma.home70.ind <- which(colnames(ma.samples$ode) == "death-prob-home-70")
  ma.home70.dens <- density(ma.samples$ode[, ma.home70.ind])
  
  pa.home70.ind <- which(colnames(pa.samples$ode) == "death-prob-home-70")
  pa.home70.dens <- density(pa.samples$ode[,pa.home70.ind])
  
  # home 80
  ri.home80.ind <- which(colnames(ri.samples$ode) == "death-prob-home-80")
  ri.home80.dens <- density(ri.samples$ode[, ri.home80.ind])
  
  ma.home80.ind <- which(colnames(ma.samples$ode) == "death-prob-home-80")
  ma.home80.dens <- density(ma.samples$ode[, ma.home80.ind])
  
  pa.home80.ind <- which(colnames(pa.samples$ode) == "death-prob-home-80")
  pa.home80.dens <- density(pa.samples$ode[, pa.home80.ind])
  
  # plot limits
  max.x <- max(
    max(ri.home60.dens$x),
    max(ri.home70.dens$x),
    max(ma.home70.dens$x),
    max(pa.home70.dens$x),
    max(ri.home80.dens$x),
    max(ma.home80.dens$x),
    max(pa.home80.dens$x)
  )
  
  min.x <- min(
    min(ri.home60.dens$x),
    min(ri.home70.dens$x),
    min(ma.home70.dens$x),
    min(pa.home70.dens$x),
    min(ri.home80.dens$x),
    min(ma.home80.dens$x),
    min(pa.home80.dens$x)
  )
  
  ####################################
  ###  C1: RI Prob of Dying at Home 
  ####################################
  
  screen(13)
  # bottom, left, top, right
  par(
    mar = c(1, 2, 2, 0),
    oma = c(0, 0, 0, 0)
  )
  
  plot(ri.home60.dens,
    xlab = "", ylab = "", main = "", bty = "l",
    xlim = c(min.x, max.x),
    cex.axis = axis.size, cex.lab = lab.size,
    col = addTrans(ri.color, 25)
  )
  
  polygon(ri.home60.dens,
    col = addTrans(ri.color, 25),
    border = addTrans(ri.color, 25),
    lwd = 2
  )
  
  lines(ri.home70.dens,
    xlab = "", ylab = "", main = "", bty = "l",
    xlim = c(min.x, max.x),
    cex.axis = axis.size, cex.lab = lab.size,
    col = addTrans(ri.color, 100)
  )
  
  polygon(ri.home70.dens,
    col = addTrans(ri.color, 100),
    border = addTrans(ri.color, 100),
    lwd = 2
  )
  
  lines(ri.home80.dens,
    xlab = "", ylab = "", main = "", bty = "l",
    xlim = c(min.x, max.x),
    cex.axis = axis.size, cex.lab = lab.size,
    col = addTrans(ri.color, 250)
  )
  
  polygon(ri.home80.dens,
    col = addTrans(ri.color, 250),
    border = addTrans(ri.color, 250),
    lwd = 2
  )
  
  legend("topright",
    inset = 0.1,
    legend = c("60-69 y/o", "70-79 y/o", "80+ y/o"),
    col = c(
      addTrans(ri.color, 25),
      addTrans(ri.color, 100),
      addTrans(ri.color, 250)
    ),
    lty = rep(1, 3), cex = 1,
    box.lty = 0, lwd = rep(3, 3)
  )
  
  ####################################
  ###  C2: MA Prob of Dying at Home 
  ####################################
  screen(14)
  # bottom, left, top, right
  par(
    mar = c(1, 2, 2, 0),
    oma = c(0, 0, 0, 0)
  )
  
  plot(ma.home70.dens,
    xlab = "", ylab = "", main = "", bty = "l",
    xlim = c(min.x, max.x),
    cex.axis = axis.size, cex.lab = lab.size,
    col = addTrans(ma.color, 100)
  )
  
  polygon(ma.home70.dens,
    col = addTrans(ma.color, 100),
    border = addTrans(ma.color, 100),
    lwd = 2
  )
  
  lines(ma.home80.dens,
    xlab = "", ylab = "", main = "", bty = "l",
    xlim = c(min.x, max.x),
    cex.axis = axis.size, cex.lab = lab.size,
    col = addTrans(ma.color, 250)
  )
  
  polygon(ma.home80.dens,
    col = addTrans(ma.color, 250),
    border = addTrans(ma.color, 250),
    lwd = 2
  )
  
  legend("topright",
    inset = 0.1,
    legend = c("70-79 y/o", "80+ y/o"),
    col = c(
      addTrans(ma.color, 100),
      addTrans(ma.color, 250)
    ),
    lty = rep(1, 2), cex = 1,
    box.lty = 0, lwd = rep(3, 2)
  )
  
  ####################################
  ###  C3: PA Prob of Dying at Home 
  ####################################
  screen(15)
  # bottom, left, top, right
  par(
    mar = c(1, 2, 2, 0),
    oma = c(0, 0, 0, 0)
  )
  
  plot(pa.home70.dens,
    xlab = "", ylab = "", main = "", bty = "l",
    xlim = c(min.x, max.x),
    cex.axis = axis.size, cex.lab = lab.size,
    col = addTrans(pa.color, 100)
  )
  
  polygon(pa.home70.dens,
    col = addTrans(pa.color, 100),
    border = addTrans(pa.color, 100),
    lwd = 2
  )
  
  lines(pa.home80.dens,
    xlab = "", ylab = "", main = "", bty = "l",
    xlim = c(min.x, max.x),
    cex.axis = axis.size, cex.lab = lab.size,
    col = addTrans(pa.color, 250)
  )
  
  polygon(pa.home80.dens,
    col = addTrans(pa.color, 250),
    border = addTrans(pa.color, 250),
    lwd = 2
  )
  
  legend(0.1, 52,
    inset = 0.1,
    legend = c("70-79 y/o", "80+ y/o"),
    col = c(
      addTrans(pa.color, 100),
      addTrans(pa.color, 250)
    ),
    lty = rep(1, 2), cex = 1,
    box.lty = 0, lwd = rep(3, 2)
  )
  
  ############################################################################
  ### D. ICU ADMISSION PROBABILITY
  ############################################################################
  
  screen(5)
  # bottom, left, top, right
  par(
    mar = c(0, 0, 2, 0),
    oma = c(0, 0, 0, 0)
  )
  
  title("D. ICU Admission Probability",
    adj = 0, line = 1, cex.main = 1.25, outer = FALSE
  )
  
  # lockdown
  ri.icu1.ind <- which(colnames(ri.samples$ode) == "dev-icu-frac")
  ri.icu1.dens <- density(icuTransformation(ri.samples$ode[, ri.icu1.ind], "RI"))

  ma.icu1.ind <- which(colnames(ma.samples$ode) == "dev-icu-frac")
  ma.icu1.dens <- density(icuTransformation(ma.samples$ode[, ma.icu1.ind], "MA"))
  
  pa.icu1.ind <- which(colnames(pa.samples$ode) == "dev-icu-frac")
  pa.icu1.dens <- density(icuTransformation(pa.samples$ode[, pa.icu1.ind], "PA"))
  
  # post-lockdown
  ri.icu2.ind <- which(colnames(ri.samples$ode) == "dev-icu-frac-phase2")
  ri.icu2.dens <- density(icuTransformation(ri.samples$ode[, ri.icu2.ind], "RI"))
  
  ma.icu2.ind <- which(colnames(ma.samples$ode) == "dev-icu-frac-phase2")
  ma.icu2.dens <- density(icuTransformation(ma.samples$ode[, ma.icu2.ind], "MA"))
  
  pa.icu2.ind <- which(colnames(pa.samples$ode) == "dev-icu-frac-phase2")
  pa.icu2.dens <- density(icuTransformation(pa.samples$ode[, pa.icu2.ind], "PA"))
  
  # plot limits
  max.x <- max(
    max(ri.icu1.dens$x),
    max(ma.icu1.dens$x),
    max(pa.icu1.dens$x),
    max(ri.icu2.dens$x),
    max(ma.icu2.dens$x),
    max(pa.icu2.dens$x)
  )
  
  min.x <- min(
    min(ri.icu1.dens$x),
    min(ma.icu1.dens$x),
    min(pa.icu1.dens$x),
    min(ri.icu2.dens$x),
    min(ma.icu2.dens$x),
    min(pa.icu2.dens$x)
  )
  
  #################################
  ###  D1: RI ICU admission prob 
  #################################
  screen(16)
  # bottom, left, top, right
  par(
    mar = c(1, 2, 2, 0),
    oma = c(0, 0, 0, 0)
  )
  
  plot(ri.icu1.dens,
    xlab = "", ylab = "", main = "", bty = "l",
    xlim = c(min.x, max.x),
    ylim = c(0, max(max(ri.icu1.dens$y), max(ri.icu2.dens$y))),
    cex.axis = axis.size, cex.lab = lab.size,
    col = addTrans(ri.color, 100),
    yaxt = "n"
  )
  
  axis(2, at = c(0, 4, 8, 12), labels = c(0, 4, 8, 12), las = 2)
  
  polygon(ri.icu1.dens,
    col = addTrans(ri.color, 100),
    border = addTrans(ri.color, 100),
    lwd = 2
  )
  
  lines(ri.icu2.dens,
    xlim = c(min.x, max.x),
    ylim = c(0, max(max(ri.icu1.dens$y), max(ri.icu2.dens$y))),
    cex.axis = axis.size, cex.lab = lab.size,
    col = addTrans(ri.color, 250)
  )
  
  polygon(ri.icu2.dens,
    col = addTrans(ri.color, 250),
    border = addTrans(ri.color, 250),
    lwd = 2
  )
  
  legend(0.225, 15.5,
    legend = c("during lockdown", "post-lockdown"),
    col = c(
      addTrans(ri.color, 100),
      addTrans(ri.color, 250)
    ),
    lty = rep(1, 2), cex = 1,
    box.lty = 0, lwd = rep(3, 2)
  )
  
  #################################
  ###  D2: MA ICU admission prob 
  #################################
  screen(17)
  # bottom, left, top, right
  par(
    mar = c(1, 2, 2, 0),
    oma = c(0, 0, 0, 0)
  )
  
  plot(ma.icu1.dens,
    xlab = "", ylab = "", main = "", bty = "l",
    xlim = c(min.x, max.x),
    ylim = c(0, max(max(ma.icu1.dens$y), max(ma.icu2.dens$y))),
    cex.axis = axis.size, cex.lab = lab.size,
    col = addTrans(ma.color, 100)
  )
  
  polygon(ma.icu1.dens,
    col = addTrans(ma.color, 100),
    border = addTrans(ma.color, 100),
    lwd = 2
  )
  
  lines(ma.icu2.dens,
    xlim = c(min.x, max.x),
    ylim = c(0, max(max(ma.icu1.dens$y), max(ma.icu2.dens$y))),
    cex.axis = axis.size, cex.lab = lab.size,
    col = addTrans(ma.color, 250)
  )
  
  polygon(ma.icu2.dens,
    col = addTrans(ma.color, 250),
    border = addTrans(ma.color, 250),
    lwd = 2
  )
  
  legend(0.22, 19,
    legend = c("during lockdown", "post-lockdown"),
    col = c(
      addTrans(ma.color, 100),
      addTrans(ma.color, 250)
    ),
    lty = rep(1, 2), cex = 1,
    box.lty = 0, lwd = rep(3, 2)
  )
  
  #################################
  ###  D3: PA ICU admission prob 
  #################################
  screen(18)
  # bottom, left, top, right
  par(
    mar = c(1, 2, 2, 0),
    oma = c(0, 0, 0, 0)
  )
  
  plot(pa.icu1.dens,
    xlab = "", ylab = "", main = "", bty = "l",
    xlim = c(min.x, max.x),
    ylim = c(0, max(max(pa.icu1.dens$y), max(pa.icu2.dens$y))),
    cex.axis = axis.size, cex.lab = lab.size,
    col = addTrans(ma.color, 100)
  )
  
  polygon(pa.icu1.dens,
    col = addTrans(pa.color, 100),
    border = addTrans(pa.color, 100),
    lwd = 2
  )
  
  lines(pa.icu2.dens,
    xlim = c(min.x, max.x),
    ylim = c(0, max(max(pa.icu1.dens$y), max(pa.icu2.dens$y))),
    cex.axis = axis.size, cex.lab = lab.size,
    col = addTrans(pa.color, 250)
  )
  
  polygon(pa.icu2.dens,
    col = addTrans(pa.color, 250),
    border = addTrans(pa.color, 250),
    lwd = 2
  )
  
  legend(0.07, 27,
    legend = c("during lockdown", "post-", "lockdown"),
    col = c(
      addTrans(pa.color, 100),
      addTrans(pa.color, 250), NA
    ),
    lty = rep(1, 3), cex = 1,
    box.lty = 0, lwd = rep(3, 3)
  )
  
  ############################################################################
  ### E. Fraction of Symptomatic Individuals Hopsitalized, By Age
  ############################################################################
  
  screen(6)  ### bottom screen
  # bottom, left, top, right
  par(
    mar = c(0, 0, 2, 0),
    oma = c(0, 0, 0, 0)
  )
  
  title("E. Fraction of Symptomatic Individuals Hospitalized, By Age Group",
    adj = 0, line = 1, cex.main = 1.25, outer = FALSE
  )
  
  ### Plot with Credible Intervals ###
  screen(19)
  
  if (logplot) {
    # bottom, left, top, right
    par(
      mar = c(4, 3, 2, 0),
      oma = c(0, 0, 0, 0)
    )
  } else {
    # bottom, left, top, right
    par(
      mar = c(2, 3, 2, 0),
      oma = c(0, 0, 0, 0)
    )
  }
  
  ri.h10.ind <- which(colnames(ri.samples$ode) == "hosp-frac-10")               
  ri.hfrac.inds <- ri.h10.ind + 0:7   
  ri.hfrac <- ri.samples$ode[, ri.hfrac.inds]
  ri.hfrac.sum <- t(apply(ri.hfrac, 2, quantile, prob = c(0.5, 0.025, 0.975)))
  
  if (logplot) {
    ri.hfrac.sum <- log2(ri.hfrac.sum)
    
    plot(
      y = 1:8 + 0.1, x = ri.hfrac.sum[, 1], ylim = c(0.5, 8.5), xlim = c(-6.5, -1.5),
      xlab = expression(paste("fraction of symptomatic individualized who are hospitalized (", log[2], ")", sep = "")),
      ylab = "", main = "", col = addTrans(ri.color, 250),
      pch = 16, bty = "l", yaxt = "n", cex = 1.25,
      cex.axis = axis.size, cex.lab = lab.size
    )    
    axis(2,
      at = 1:8, labels = c(
        "10-19", "20-29", "30-39", "40-49", "50-59",
        "60-69", "70-79", "80+"
      ),
      cex.axis = axis.size, cex.lab = lab.size, las = 1
    )
    
    abline(h = 1, col = addTrans("grey", 50), lwd = 2)
    abline(h = 2, col = addTrans("grey", 50), lwd = 2)
    abline(h = 3, col = addTrans("grey", 50), lwd = 2)
    abline(h = 4, col = addTrans("grey", 50), lwd = 2)
    abline(h = 5, col = addTrans("grey", 50), lwd = 2)
    abline(h = 6, col = addTrans("grey", 50), lwd = 2)
    abline(h = 7, col = addTrans("grey", 50), lwd = 2)
    abline(h = 8, col = addTrans("grey", 50), lwd = 2)
    
    segments(
      y0 = 1:8 + 0.1, x0 = ri.hfrac.sum[, 2],
      y1 = 1:8 + 0.1, x1 = ri.hfrac.sum[, 3],
      col = addTrans(ri.color, 100),
      lty = 1, lwd = 3
    )
    
    points(
      y = 1:8 + 0.1, x = ri.hfrac.sum[, 1],
      col = addTrans(ri.color, 250),
      pch = 16, cex = 1.25
    )
    
    legend("bottomright",
      inset = 0.15,
      legend = c("Rhode Island"), # , "Massachusetts", "Pennsylvania"),
      col = c(addTrans(ri.color, 250)),
      # addTrans(ma.color, 250),
      # addTrans(pa.color, 250)),
      pch = rep(16, 1), cex = 1.25,
      box.lty = 0
    )
    
  } else {
    
    plot(
      y = 1:8 + 0.1, x = ri.hfrac.sum[, 1], ylim = c(0.5, 8.5), xlim = c(0, 0.35),
      xlab = "",
      ylab = "", main = "", col = addTrans(ri.color, 250),
      pch = 16, bty = "l", yaxt = "n", cex = 1.25,
      cex.axis = axis.size, cex.lab = lab.size
    )
    
    axis(2,
      at = 1:8, labels = c(
        "10-19", "20-29", "30-39", "40-49", "50-59",
        "60-69", "70-79", "80+"
      ),
      cex.axis = axis.size, cex.lab = lab.size, las = 1
    )
    
    abline(h = 1, col = addTrans("grey", 50), lwd = 2)
    abline(h = 2, col = addTrans("grey", 50), lwd = 2)
    abline(h = 3, col = addTrans("grey", 50), lwd = 2)
    abline(h = 4, col = addTrans("grey", 50), lwd = 2)
    abline(h = 5, col = addTrans("grey", 50), lwd = 2)
    abline(h = 6, col = addTrans("grey", 50), lwd = 2)
    abline(h = 7, col = addTrans("grey", 50), lwd = 2)
    abline(h = 8, col = addTrans("grey", 50), lwd = 2)
    
    segments(
      y0 = 1:8 + 0.1, x0 = ri.hfrac.sum[, 2],
      y1 = 1:8 + 0.1, x1 = ri.hfrac.sum[, 3],
      col = addTrans(ri.color, 100),
      lty = 1, lwd = 3
    )
    
    points(
      y = 1:8 + 0.1, x = ri.hfrac.sum[, 1],
      col = addTrans(ri.color, 250),
      pch = 16, cex = 1.25
    )
    
    legend("topleft",
      legend = c("Rhode Island"), # "Massachusetts", "Pennsylvania"),
      col = c(addTrans(ri.color, 250)),
      # addTrans(ma.color, 250),
      # addTrans(pa.color, 250)),
      pch = rep(16, 1), cex = 1.25,
      box.lty = 0
    )
  }

  close.screen(all = TRUE)
  
  dev.off()
  embed_fonts(plot.name, outfile = plot.name)
}


