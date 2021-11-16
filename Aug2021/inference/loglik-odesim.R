#!/usr/bin/env Rscript

### loglik-odesim.R
### last edited: 11 Nov 2021
### authors: Ephraim Hanks, Nathan Wikle


### 1. Determine which multinomial PMF evaluation to use:
###   If OPT_USE_PY_LL = TRUE, use a Python evaluation (much faster!)
###   Else, use the R default (dmultinom)

# Get global options (i.e., determine if Python should be used)
OPT_USE_PY_LL <- get0('OPT_USE_PY_LL', ifnotfound = TRUE)

# Helper function 1 - Get sums of time spans bounded by times
time_slicer <- function(vals, times) {
  out <- matrix(nrow = length(times) - 1, ncol = ncol(vals))
  for (kk in 2:length(times)) {
    this_slice <- vals[(times[kk - 1] + 1):times[[kk]], ]
    if (length(dim(this_slice)) == 2) {
      this_slice <- apply(this_slice, 2, sum)
    }
    out[kk - 1, ] <- this_slice
  }
  return(out)
}

# Helper function 2 - Sum up log-likelihood for several rows
r_ll <- function(means, counts){
  row_lls <- sapply(
    1:nrow(means),
    function(i) dmultinom(round(counts[i, ]), prob = means[i, ], log = TRUE)
  )
  return(sum(row_lls))
}

# Helper function 3 - Python-based multinomial logpmf
if (OPT_USE_PY_LL) {
  # Use Python multinomial log-pmf (MUCH FASTER)  
  tryCatch({
      library(reticulate)
      source_python("./py_faster_stats.py") # path is relative to this file, remember to source with chdir = TRUE
    }, error = function(e) {
      stop(c("Unable to import Python multinomial function. Set OPT_USE_PY_LL=FALSE to use R version.", "\n\n", e))
    }
  )
  py_ll <- py$py_ll
  multinom_log_pmf <- function(x, y) py_ll(as.matrix(x), as.matrix(y))
} else {
  # Use dmultinom() instead (slower)
  multinom_log_pmf <- r_ll
}

### 2. Define the loglikelihood for our model. See the Supplementary Materials
###       for a detailed description of our choice of likelihood.
loglik.odesim <- function(
  traj,               # ODESIM trajectory output
  df,                 # observed data
  dp = NULL,          # processed data
  odesim.ver = "v5",  # odesim version
  P = P,              # matrix with probablity of delay from symptoms to presentation
  loc = loc,          # state of interest ("RI", "MA", or "PA")
  report.rate,        # vector of reporting rate
  nb.disp.params = NULL,      # NB dispersion parameters
  extra.params = NULL,        # extra (non-ODESIM) parameters (estimated)
  extra.const.params = NULL,  # extra (non-ODESIM) parameters (constant) 
  s2.hosp = NULL,     # variance parameter, current hosps
  s2.icu = NULL,      # variance parameter, current ICU
  s2.vent = NULL,     # variance parameter, current vent
  lik.tot = TRUE,               # calculate likelihood for total cases only (Boolean)
  lik.age = FALSE,              # calculate likelihood for total AND age-stratified data (Boolean)
  lik.hosp.new = FALSE,         # calculate likelihood for new hospitalizaitons (Boolean)
  lik.tot.deaths = FALSE,       # calculate likelihood for total deaths ONLY (Boolean)
  lik.age.deaths = FALSE,       # calculate likelihood for total AND age-stratified deaths (Boolean)
  lik.home.deaths = FALSE,      # calculate likelihood for home deaths (Boolean)
  lik.hosp.deaths = FALSE,      # calculate likelihood for hosp. deaths (Boolean)
  lik.hosp.discharges = FALSE,  # calculate likelihood for new hosp. discharges (Boolean)
  lik.hosp.curr = TRUE,         # calculate likelihood for current hospitalizations (Boolean)
  lik.icu.curr = TRUE,          # calculate likelihood for current ICU data (Boolean)
  lik.vent.curr = TRUE,         # calculate likelihood for current vent data (Boolean)
  lik.curr = FALSE,             # use weekly, independent likelihood for current data (Boolean)
  lik.old = FALSE,              # if !lik.curr & lik.old: cond. likelihood based on increments
                                # if !lik.curr & !lik.old: cond. joint normal approximation
  active.surv = TRUE, # use active surveillance data (Boolean)
  p.asympt = .4,      # asymptomatic proportion (only used if active.surv = TRUE)
  total.size.constraint = FALSE,  # require estimated cumulative cases to be within 10% of observed data (Boolean)
  ...
){
  
  #######################################################################
  ### 1. handling "extra parameters" 
  #######################################################################

  # hospital reporting rate is assumed to be 1, unless otherwise specified in "extra.params"
  hosp.rr <- 1
  
  # read in values from "extra.params" (to be estimated!)
  if (length(extra.params) > 0) {
    extra.names <- names(extra.params)

    # hosp.report.rate
    if (is.element("hosp.report.rate", extra.names)) {
      hosp.rr <- extra.params["hosp.report.rate"]
    }

    # proportion of asymptomatics
    if (is.element("p.asympt", extra.names)) {
      p.asympt <- extra.params["p.asympt"]
    }
  }
  
  # read in values from "extra.params" (held constant!)
  if (length(extra.const.params) > 0) {
    extra.const.names <- names(extra.const.params)
    # hosp.report.rate
    if (is.element("hosp.report.rate", extra.const.names)) {
      hosp.rr <- as.numeric(extra.const.params["hosp.report.rate"])
    }
  }

  #######################################################################  
  ### 2. parse out NB dispersions params 
  #######################################################################
  
  nb.disp.new.cases <- NULL     # leave "NULL" to do Poisson lik
  nb.disp.new.hosp <- NULL      # leave "NULL" to do Poisson lik
  nb.disp.deaths <- NULL        # leave "NULL" to do Poisson lik
  nb.disp.home.deaths <- NULL   # leave "NULL" to do Poisson lik
  nb.disp.hosp.deaths <- NULL   # leave "NULL" to do Poisson lik

  # fill in NB dispersion parameters
  if (length(nb.disp.params) >= 1) {
    # new cases
    nb.disp.new.cases <- nb.disp.params[1]
  }
  if (length(nb.disp.params) >= 2) {
    # new hospitalizations
    nb.disp.new.hosp <- nb.disp.params[2]
  }
  if (length(nb.disp.params) >= 3 & (lik.tot.deaths | lik.age.deaths)) {
    # new deaths
    nb.disp.deaths <- nb.disp.params[3]
  }
  if (length(nb.disp.params) == 4 & lik.home.deaths & lik.hosp.deaths) {
    # new home and hosp deaths (no discharges)
    nb.disp.home.deaths <- nb.disp.params[3]
    nb.disp.hosp.deaths <- nb.disp.params[4]
  }
  if (length(nb.disp.params) == 4 & lik.hosp.discharges) {
    # new hospital discharges (no home/hosp deaths)
    nb.disp.hosp.discharges <- nb.disp.params[4]
  }
  if (length(nb.disp.params) == 5 & lik.home.deaths & lik.hosp.deaths & lik.hosp.discharges) {
    # new home/hosp deaths + new hospital discharges
    nb.disp.home.deaths <- nb.disp.params[3]
    nb.disp.hosp.deaths <- nb.disp.params[4]
    nb.disp.hosp.discharges <- nb.disp.params[5]
  }
    
  #######################################################################  
  ### 3. pre-processing of data 
  #######################################################################  
  
  # symptomatic reporting rate (for new cases)
  report.rate <- as.numeric(report.rate)
  
  # pre-process epidemic data
  if (is.null(dp)) {
    data.processed <- data.process(df, loc = loc, ...)
  } else {
    data.processed <- dp
  }
  
  # add elements of data.processed to global environment
  list2env(data.processed, globalenv())
  n.days <- length(days)
  
  #######################################################################
  ### 4. pre-processing of ODESIM output
  #######################################################################

  # pre-process ODESIM trajectory output (traj)
  traj.processed <- traj.process(traj, loc = loc, odesim.version = odesim.ver)
  # add elements of traj.processed to global environment
  list2env(traj.processed, globalenv())
  
  # link time indices of data and odesim output
  traj.times <- round(traj[,1])
  tmin <- min(days); tmax <- max(days)
  idx.min <- which(traj.times == tmin)
  idx.max <- which(traj.times == tmax)
  t.idx <- idx.min:idx.max
  
  # link odesim output with report.rate (which is daily starting at day 61)
  t.idx.rr <- t.idx - 10
  
  # linking to P and delay rate
  t.idx.delay <- tmin:tmax
  
  #######################################################################
  ### 5. Check that parameters are in-bounds, and that the total epidemic
  ###     size is within 10% of the observed total
  #######################################################################

  # cumulative number of cases (observed)
  tot.size.ep <- sum(tot.sympt.new, na.rm = T)
  # estimated (ODESIM) cumulative number of cases
  tot.size.ep.odesim <- cumsum(tot.sympt.new.odesim[t.idx] * report.rate[t.idx.rr])
  tot.size.ep.odesim <- max(tot.size.ep.odesim, na.rm = TRUE)
  # percent difference between observed and ODESIM cumulative cases
  tot.size.percent.off <- abs(tot.size.ep.odesim - tot.size.ep) / tot.size.ep
  
  # check if conditions are met
  if(
    min(report.rate, 
        nb.disp.new.cases, 
        nb.disp.new.hosp,
        s2.hosp, s2.icu, s2.vent) < 0 | (total.size.constraint & tot.size.percent.off > 0.10)
  ){ 
    # invalid likelihood (parameter(s) out of bounds or size constraint not met)
    ll <- -Inf
    ll.hosp.new <- NA
    ll.hosp <- NA
    ll.vent <- NA
    ll.icu <- NA
    ll.new <- -Inf
    ll.deaths <- NA
    ll.hosp.dis <- NA
    sympt.imputed <- NULL
    
  } else {
    # calculate likelihood!

    # initialize outputs
    ll.new <- integer()
    ll.hosp <- integer()
    ll.deaths <- integer()
    sympt.imputed <- integer()
    ll <- 0

    #######################################################################
    ### 6. likelihood for daily new cases
    #######################################################################

    if (lik.age & !active.surv) {
      ### likelihood for age-structured new cases (no active surveillance data)

      # account for testing delay:
      delay.rates <- as.vector(P[, traj.times] %*% tot.sympt.new.odesim) + .Machine$double.eps[1]
      # finding any na.values
      no.na <- which(!is.na(tot.sympt.new))
      # log-likelihood of total new cases each day:
      ll.tot <- sum(dnbinom(
        tot.sympt.new[no.na],
        mu = report.rate[t.idx.rr[no.na]] * delay.rates[t.idx.delay[no.na]],
        size = nb.disp.new.cases, log = TRUE
      ))
      # log-likelihood for age-structured data
      age.obs.times <- c(0, age.obs.times)

      # shift age means
      shifted.age.means <- P[, traj.times] %*% sympt.new.odesim + .Machine$double.eps[1]
      # account for reporting rate
      age.means <- apply(
        shifted.age.means, 2,
        function(X) {
          report.rate[t.idx.rr] * X[t.idx.delay]
        }
      )

      odesim.means.mat <- time_slicer(age.means, age.obs.times)

      zero.cases <- rowSums(age.new.cases) == 0
      ll.age <- multinom_log_pmf(
        odesim.means.mat[!zero.cases, ],
        age.new.cases[!zero.cases, ]
      )

      # add total and age-structured loglikelihoods
      ll.new <- ll.tot + ll.age
      ll <- ll + ll.new
    }
    
    if (lik.age & active.surv) {
      ### likelihood for age-structured new cases (WITH active surveillance data)

      # active surveillance constants
      sigma.e <- 0.8 
      sigma.a <- 0.9 
      p.ncas <- (1 - sigma.e * c.t[t.idx.rr - 1, ])^6
      
      # making SEIAR for all age classes
      i9 <- diag(9)
      S <- traj[t.idx.rr - 1, 2:10]
      E <- traj[t.idx.rr - 1, c(11:19, 20:28, 29:37, 38:46, 47:55, 56:64)] %*% rbind(i9, i9, i9, i9, i9, i9)
      A <- traj[t.idx.rr- 1, c(65:73, 74:82, 83:91, 92:100)] %*% rbind(i9, i9, i9, i9)
      I <- traj[t.idx.rr - 1, c(101:109, 110:118, 119:127, 128:136)] %*% rbind(i9, i9, i9, i9)
      R <- traj[t.idx.rr - 1, 272:280]
      prevalence <- (sigma.e * E + (1 - report.rate[t.idx.rr]) * I + sigma.a * A) /
        (S + E + (1 - report.rate[t.idx.rr]) * I + A + (p.asympt + (1 - p.asympt) * (1 - report.rate[t.idx.rr])) * R)
            
      if (loc == "RI"){
        # used to calculate expected number of positive tests from active surveillance
        pop.by.age <- c(111233, 130301, 148311, 134539, 131361, 143014, 127123, 78393, 55087)
      }
      
      # calculate mean new cases (passive and active)
      delay.rates <- P[,traj.times] %*% sympt.new.odesim + .Machine$double.eps[1]
      mean.passive.active <- report.rate[t.idx.rr] * delay.rates[t.idx.delay, ] * p.ncas +
        (c.t[t.idx.rr - 1, ] * prevalence) %*% diag(pop.by.age)
      tot.mean.pass.act <- apply(mean.passive.active, 1, sum)

      # finding any na.values
      no.na <- which(!is.na(tot.sympt.new))
      # log-likelihood of total new cases each day:
      ll.tot <- sum(dnbinom(tot.sympt.new[no.na], mu = tot.mean.pass.act[no.na], size = nb.disp.new.cases, log = TRUE))
      # log-likelihood for age-structured data
      age.obs.times <- c(0, age.obs.times)
      
      odesim.means.mat <- time_slicer(mean.passive.active, age.obs.times)
      ll.age <- multinom_log_pmf(odesim.means.mat, age.new.cases)
      
      # add total and age-structured loglikelihoods
      ll.new <- ll.tot + ll.age
      ll <- ll + ll.new
    }
    
    if (lik.tot) {
      ### likelihood for total new cases (not age-structured)

      # account for testing delay:
      delay.rates <- as.vector(P[, traj.times] %*% tot.sympt.new.odesim) + .Machine$double.eps[1]
      # finding any na.values
      no.na <- which(!is.na(tot.sympt.new))
      # log-likelihood given trajectory sim:
      ll.new <- sum(dnbinom(tot.sympt.new[no.na],
        mu = report.rate[t.idx.rr[no.na]] * delay.rates[t.idx.delay[no.na]],
        size = nb.disp.new.cases, log = TRUE
      ))
      ll <- ll + ll.new
    }
    
    #######################################################################
    ### 6. likelihood for new hospitalizations
    #######################################################################
    if (lik.hosp.new & lik.age){
      
      ### new hospitalizations (including age-stratified data)

      # find any na.values
      no.na <- which(!is.na(tot.hosp.new))
      # log-likelihood of total new hospitalizations each day:
      ll.hosp.tot <- sum(dnbinom(tot.hosp.new[no.na],
        mu = tot.hosp.new.odesim[no.na] * hosp.rr + .Machine$double.eps[1],
        size = nb.disp.new.hosp, log = TRUE
      ))
      
      # log-likelihood for age-structured data
      age.hosp.obs.times <- c(0, age.hosp.obs.times)
      # set negative values to zero
      hosp.new.odesim[hosp.new.odesim <= 0] <- .Machine$double.eps[1]
      ll.hosp.age <- 0
      
      for(kk in 2:length(age.hosp.obs.times)){
        odesim.means <- hosp.new.odesim[(age.hosp.obs.times[kk - 1] + 1):age.hosp.obs.times[[kk]], ]
        if (length(dim(odesim.means)) == 2) {
          odesim.means <- apply(odesim.means, 2, sum)
        }
        
        if (loc == "PA") {
          na.inds <- which(is.na(round(as.numeric(age.hosp.new[kk - 1, ]))))
          age.hosp.new.nona <- age.hosp.new[, -na.inds]
          odesim.means.nona <- odesim.means[-na.inds]
          ll.hosp.age <- ll.hosp.age + dmultinom(round(as.numeric(age.hosp.new.nona[kk - 1, ])),
            prob = as.numeric(odesim.means.nona), log = TRUE
          )
        } else {
          ll.hosp.age <- ll.hosp.age + dmultinom(round(as.numeric(age.hosp.new[kk - 1, ])),
            prob = as.numeric(odesim.means), log = TRUE
          )
        }
      }
      
      # combine age and total loglikelihood
      ll.hosp.new <- ll.hosp.tot + ll.hosp.age
      ll <- ll + ll.hosp.new

    } else if (lik.hosp.new & !lik.age) {
      
      ### new hospitalizations (WITHOUT age breakdown)
      
      # find any na.values
      no.na <- which(!is.na(tot.hosp.new))
      # log-likelihood of total new hospitalizations each day:
      ll.hosp.tot <- sum(dnbinom(tot.hosp.new[no.na],
        mu = tot.hosp.new.odesim[no.na] * hosp.rr + .Machine$double.eps[1], 
        size = nb.disp.new.hosp, log = TRUE
      ))
      
      # log-likelihood for age-structured data
      age.hosp.obs.times <- c(0, age.hosp.obs.times)
      # set negative values to zero
      hosp.new.odesim[hosp.new.odesim <= 0] <- .Machine$double.eps[1]
      
      odesim.means.mat <- time_slicer(hosp.new.odesim, age.hosp.obs.times)
      ll.hosp.age <- multinom_log_pmf(odesim.means.mat, age.hosp.new)
      
      # combine age and total likelihood
      ll.hosp.new <- ll.hosp.tot + ll.hosp.age
      ll <- ll + ll.hosp.new

    } else {

      # no new hospitalization data
      ll.hosp.new <- 0
      ll <- ll + ll.hosp.new
    }

    #######################################################################
    ### 7. likelihood for new deaths
    #######################################################################
    if (lik.tot.deaths) {
      ### new deaths (total, NO age data)

      # finding any na.values
      no.na <- which(!is.na(tot.deaths.new))
      # log-likelihood given trajectory sim:
      tot.deaths.new.odesim[tot.deaths.new.odesim <= 0] <- 0 + .Machine$double.eps[1]
      ll.deaths <- sum(dnbinom(tot.deaths.new[no.na],
        mu = tot.deaths.new.odesim[t.idx[no.na]], size = nb.disp.deaths,
        log = TRUE, na.rm = TRUE
      ))
      ll <- ll + ll.deaths
    }

    if (lik.age.deaths) {
      ### new deaths (total AND age data)

      # finding any na.values
      no.na <- which(!is.na(tot.deaths.new))
      # log-likelihood for total deaths:
      tot.deaths.new.odesim[tot.deaths.new.odesim <= 0] <- 0 +.Machine$double.eps[1]
      ll.tot.deaths <- sum(dnbinom(tot.deaths.new[no.na],
        mu = tot.deaths.new.odesim[t.idx[no.na]], size = nb.disp.deaths, log = TRUE
      ))

      # home and hosp deaths
      ll.home.hosp.deaths <- 0
      if (lik.hosp.deaths & lik.home.deaths & !(loc == "PA")) {
        tot.hosp.deaths.new[tot.hosp.deaths.new < 0] = 0
        tot.home.deaths.new[tot.home.deaths.new < 0] = 0
        no.na <- which(!is.na(tot.hosp.deaths.new + tot.home.deaths.new))
        ll.home.hosp.deaths <- sum(dmultinomial(as.matrix(cbind(
          tot.hosp.deaths.new,
          tot.home.deaths.new
        )[no.na, ]),
        prob = as.matrix(cbind(
          tot.hosp.deaths.new.odesim[t.idx],
          tot.home.deaths.new.odesim[t.idx]
        ))[no.na, ],
        log = TRUE
        ))
      }
    
      if (lik.hosp.deaths & lik.home.deaths & (loc == "PA")) {
        ll.home.hosp.deaths <- 0
        home.deaths.cum <- tot.deaths.cum[hosp.deaths.obs.times] - cumsum(hosp.deaths.new)
        home.deaths.new <- c(
          home.deaths.cum[1],
          home.deaths.cum[-1] - home.deaths.cum[-length(home.deaths.cum)]
        )

        home.deaths.new[home.deaths.new < 0] <- 0

        for (kk in 1:length(hosp.deaths.obs.times)) {
          if (kk == 1) {
            odesim.means <- c(
              sum(tot.hosp.deaths.new.odesim[1:hosp.deaths.obs.times[kk]]),
              sum(tot.home.deaths.new.odesim[1:hosp.deaths.obs.times[kk]])
            ) + .Machine$double.eps[1]
          } else {
            odesim.means <- c(
              sum(tot.hosp.deaths.new.odesim[(hosp.deaths.obs.times[kk - 1] + 1):hosp.deaths.obs.times[[kk]]]),
              sum(tot.home.deaths.new.odesim[(hosp.deaths.obs.times[kk - 1] + 1):hosp.deaths.obs.times[[kk]]])
            ) + .Machine$double.eps[1]
          }

          if (!is.na(hosp.deaths.new[[kk]] + home.deaths.new[[kk]])) {
            ll.home.hosp.deaths <- ll.home.hosp.deaths +
              dmultinom(round(c(hosp.deaths.new[kk], home.deaths.new[kk])),
                prob = as.numeric(odesim.means), log = TRUE
              )
          }
        }
      }

      # combine likelihoods
      ll.deaths.age <- 0
      # log-likelihood for age-structured data
      age.deaths.obs.times <- c(0, age.deaths.obs.times)
      # set negative values to zero
      deaths.new.odesim[deaths.new.odesim <= 0] <- .Machine$double.eps[1]
      odesim.means.mat <- time_slicer(deaths.new.odesim, age.deaths.obs.times)
      no.zero.deaths <- which(rowSums(age.deaths.new, na.rm = T) != 0)
      
      ll.deaths.age <- multinom_log_pmf(
        odesim.means.mat[no.zero.deaths, ],
        age.deaths.new[no.zero.deaths, ]
      )
      
      ll.deaths <- ll.tot.deaths + ll.deaths.age + ll.home.hosp.deaths
      ll <- ll + ll.deaths
      ll.vent <- NA
      ll.icu <- NA
    }

    
    
        
    
    #######################################################################
    ### 8. likelihood for current hospitalizations, ICU, and vent numbers
    #######################################################################

    if (lik.curr) { 
    
      ### A. likelihood for current (not increments/conditional) hosp/icu/vent
      ###      Note: this likelihood was found to be "best" in our paper!

      # current vent data
      if (lik.vent.curr) {
        # determine weekly times with data
        idx.curr.data <- which(tot.vent.curr > -1)
        idx.weekly.times <- min(idx.curr.data)
        keep.going <- TRUE

        while (keep.going) {
          idx.new.time <- which(idx.curr.data > max(idx.weekly.times) + 6)
          if (length(idx.new.time) > 0) {
            idx.weekly.times <- c(idx.weekly.times, min(idx.curr.data[idx.new.time]))
          } else {
            keep.going <- FALSE
          }
        }

        # calculate likelihood (normal, with constant variance)
        ll.vent <- sum(dnorm(tot.vent.curr[idx.weekly.times],
          tot.vent.curr.odesim[t.idx[idx.weekly.times]],
          sd = sqrt(s2.vent),
          log = TRUE
        ))
        # add to current ll
        ll <- ll + ll.vent
      } else {
        ll.vent <- integer(0)
      }

      # current ICU data
      if (lik.icu.curr) {
        # determine weekly times with data
        idx.curr.data <- which(tot.icu.curr > -1)
        idx.weekly.times <- min(idx.curr.data)
        keep.going <- TRUE

        while (keep.going) {
          idx.new.time <- which(idx.curr.data > max(idx.weekly.times) + 6)
          if (length(idx.new.time) > 0) {
            idx.weekly.times <- c(idx.weekly.times, min(idx.curr.data[idx.new.time]))
          } else {
            keep.going <- FALSE
          }
        }

        # calculate likelihood (normal, with constant variance)
        ll.icu <- sum(dnorm(tot.icu.curr[idx.weekly.times],
          tot.icu.curr.odesim[t.idx[idx.weekly.times]],
          sd = sqrt(s2.icu),
          log = TRUE
        ))
        # add to current ll
        ll <- ll + ll.icu
      } else {
        ll.icu <- integer(0)
      }
      

      # current hospitalization data
      if (lik.hosp.curr) {
        # determine weekly times with data
        idx.curr.data <- which(tot.hosp.curr > -1)
        idx.weekly.times <- min(idx.curr.data)
        keep.going <- TRUE

        while (keep.going) {
          idx.new.time <- which(idx.curr.data > max(idx.weekly.times) + 6)
          if (length(idx.new.time) > 0) {
            idx.weekly.times <- c(idx.weekly.times, min(idx.curr.data[idx.new.time]))
          } else {
            keep.going <- FALSE
          }
        }

        # calculate likelihood (normal, with constant variance)
        ll.hosp <- sum(dnorm(tot.hosp.curr[idx.weekly.times],
          tot.hosp.curr.odesim[t.idx[idx.weekly.times]],
          sd = sqrt(s2.hosp), log = TRUE
        ))
        # add to current ll
        ll <- ll + ll.hosp
      } else {
        ll.hosp <- integer(0)
      }
      
    } else { 
    
      ### B. likelihood for increments / conditional likelihood of hosp/icu/vent 
      
      ### B(i). likelihood used in the original submission 
      if (lik.old) {
        
        # daily change in vent numbers
        delta.vent <- c(tot.vent.curr[1], 
                        tot.vent.curr[-1] - tot.vent.curr[-n.days])
        delta.vent.odesim <- c(tot.vent.curr.odesim[1],
                               tot.vent.curr.odesim[-1] - tot.vent.curr.odesim[-n.days])
        # daily change in icu+vent numbers
        delta.icu <- c(tot.icu.curr[1],
                       tot.icu.curr[-1] - tot.icu.curr[-n.days])
        delta.icu.odesim <- c(tot.icu.curr.odesim[1],
                              tot.icu.curr.odesim[-1] - tot.icu.curr.odesim[-n.days])
        # daily change in icu (ignoring vent) numbers
        delta.icu.novent.odesim <- c(tot.icu.curr.novent.odesim[1],
                                     tot.icu.curr.novent.odesim[-1] - tot.icu.curr.novent.odesim[-n.days])
        # daily change in icu+vent+hosp numbers
        delta.hosp <- c(tot.hosp.curr[1], 
                        tot.hosp.curr[-1] - tot.hosp.curr[-n.days])
        delta.hosp.odesim <- c(tot.hosp.curr.odesim[1], 
                               tot.hosp.curr.odesim[-1] - tot.hosp.curr.odesim[-n.days])
        # daily change in hosp (ignorming icu/vent) numbers
        delta.hosp.noicu.odesim <- c(tot.hosp.curr.noicu.odesim[1], 
                                     tot.hosp.curr.noicu.odesim[-1] - tot.hosp.curr.noicu.odesim[-n.days])
        
        icu.no.na <- tot.icu.curr
        icu.no.na[tot.icu.curr > -1] <- 1
        hosp.no.na <- tot.hosp.curr
        hosp.no.na[tot.hosp.curr > -1] <- 1
        
        if (lik.vent.curr) {
          no.na.idx <- which(!is.na(delta.vent))
          # Gaussian likelihood
          ll.vent <- sum(dnorm(delta.vent[no.na.idx],
            mean = delta.vent.odesim[t.idx[no.na.idx]],
            sd = sqrt(s2.vent),
            log = TRUE
          ))
          # add to current ll
          ll <- ll + ll.vent
        } else {
          ll.vent <- integer(0)
        }
        
        if (lik.icu.curr) {
          # index for times with both icu and vent data
          idx.icu.vent <- which(!is.na(delta.icu) & !is.na(delta.vent))
          idx.icu.only <- which(!is.na(delta.icu) & is.na(delta.vent))
          # Gaussian likelihood
          ll.icu.vent <- 0
          ll.icu.only <- 0

          if (length(idx.icu.vent) > 0) {
            ll.icu.vent = sum(dnorm(delta.icu[idx.icu.vent] - delta.vent[idx.icu.vent],
              mean = delta.icu.novent.odesim[t.idx[idx.icu.vent]],
              sd = sqrt(s2.icu),
              log = TRUE
            ))
          }

          if (length(idx.icu.only) > 0) {
            ll.icu.only <- sum(dnorm(delta.icu[idx.icu.only],
              mean = delta.icu.odesim[t.idx[idx.icu.only]],
              sd = sqrt(s2.icu + s2.vent),
              log = TRUE
            ))
          }

          # add to current ll
          ll.icu <- ll.icu.vent + ll.icu.only
          ll <- ll + ll.icu
        } else {
          ll.icu <- integer(0)
        }
        
        if (lik.hosp.curr) {
          no.na.idx <- which(!is.na(tot.hosp.curr))
          # index for times with both hosp and icu or vent data
          idx.hosp.icu <- which(!is.na(delta.hosp + delta.icu))
          idx.hosp.vent <- which(!is.na(delta.hosp + delta.vent) & is.na(delta.icu))
          idx.hosp.only <- which(!is.na(delta.hosp) & is.na(delta.vent) & is.na(delta.icu))
          # Gaussian likelihood
          ll.hosp.vent <- 0
          ll.hosp.icu <- 0
          ll.hosp.only <- 0
          # case in which hosp and icu are observed
          if (length(idx.hosp.icu) > 0) {
            ll.hosp.icu = sum(dnorm(delta.hosp[idx.hosp.icu] - delta.icu[idx.hosp.icu],
              mean = delta.hosp.noicu.odesim[t.idx[idx.hosp.icu]],
              sd = sqrt(s2.hosp),
              log = TRUE
            ))
          }
          # case in which hosp and vent, but not icu, are observed
          if (length(idx.hosp.vent) > 0) {
            ll.hosp.vent <- sum(dnorm(delta.hosp[idx.hosp.vent] - delta.vent[idx.hosp.vent],
              mean = delta.hosp.odesim[t.idx[idx.hosp.vent]] - delta.vent.odesim[t.idx[idx.hosp.vent]],
              sd = sqrt(s2.hosp + s2.vent),
              log = TRUE
            ))
          }
          # case in which hosp is only observation
          if (length(idx.hosp.only) > 0) {
            ll.hosp.only <- sum(dnorm(delta.hosp[idx.hosp.only],
              mean = delta.hosp.odesim[t.idx[idx.hosp.only]],
              sd = sqrt(s2.hosp + s2.icu + s2.vent),
              log = TRUE
            ))
          }
          # add to current ll
          ll.hosp <- ll.hosp.icu + ll.hosp.vent + ll.hosp.only
          ll <- ll + ll.hosp
        } else {
          ll.hosp <- integer(0)
        }
          
      } else {

        ### B(ii). "new" likelihood, with dependence across current data streams  
        
        # daily change in vent numbers
        delta.vent <- c(tot.vent.curr[1], 
                        tot.vent.curr[-1] - tot.vent.curr[-n.days])
        delta.vent.odesim <- c(tot.vent.curr.odesim[1], 
                               tot.vent.curr.odesim[-1] - tot.vent.curr.odesim[-n.days])
        
        # daily change in icu+vent numbers
        delta.icu <- c(tot.icu.curr[1], 
                       tot.icu.curr[-1] - tot.icu.curr[-n.days])
        delta.icu.odesim <- c(tot.icu.curr.odesim[1], 
                              tot.icu.curr.odesim[-1] - tot.icu.curr.odesim[-n.days])
        
        # daily change in icu (ignoring vent) numbers
        delta.icu.novent.odesim <- c(tot.icu.curr.novent.odesim[1],
                                     tot.icu.curr.novent.odesim[-1] - tot.icu.curr.novent.odesim[-n.days])
        
        # daily change in icu+vent+hosp numbers
        delta.hosp <- c(tot.hosp.curr[1],
                        tot.hosp.curr[-1] - tot.hosp.curr[-n.days])
        delta.hosp.odesim <- c(tot.hosp.curr.odesim[1],
                               tot.hosp.curr.odesim[-1] - tot.hosp.curr.odesim[-n.days])
        
        # daily change in hosp (ignoring icu/vent) numbers
        delta.hosp.noicu.odesim <- c(tot.hosp.curr.noicu.odesim[1],
                                     tot.hosp.curr.noicu.odesim[-1] - tot.hosp.curr.noicu.odesim[-n.days])
          
        icu.no.na <- tot.icu.curr
        icu.no.na[tot.icu.curr > -1] <- 1
        hosp.no.na <- tot.hosp.curr
        hosp.no.na[tot.hosp.curr > -1] <-1
          
        if (lik.hosp.curr) {
              
          ### getting mean net movement from odesim output
          
          # net movement between icu and non-icu
          mu.ihhi <- delta.hosp.noicu.odesim - tot.hosp.new.odesim*hosp.rr + tot.hosp.discharges.new.odesim
          # net movement between icu and ventilated group
          mu.viiv <- -delta.vent.odesim + tot.deaths.new.odesim
              
          # making mean vector for all time points
          mu.hiv <- matrix(NA,nrow=length(t.idx),ncol=3)
              
          if (loc == "RI") {
            # set up data for RI:
            mu.hiv[, 1] <- tot.hosp.new - tot.hosp.discharges.new + mu.ihhi[t.idx]
            mu.hiv[, 2] <- -mu.ihhi[t.idx] + mu.viiv[t.idx]
            mu.hiv[, 3] <- -mu.viiv[t.idx] - tot.hosp.deaths.new
          } else {
            # for PA / MA, use odesim output instead of data
            mu.hiv[, 1] <- tot.hosp.new.odesim[t.idx] * hosp.rr -
              tot.hosp.discharges.new.odesim[t.idx] + mu.ihhi[t.idx]
            mu.hiv[, 2] <- -mu.ihhi[t.idx] + mu.viiv[t.idx]
            mu.hiv[, 3] <- -mu.viiv[t.idx] - tot.hosp.deaths.new.odesim[t.idx]
          }
            
          # make covariance matrix for all time points
          cov.hiv <- array(0,dim=c(3,3,length(t.idx)))
          cov.hiv[1,2,] <- -s2.icu
          cov.hiv[2,1,] <- -s2.icu
          cov.hiv[2,3,] <- -s2.vent
          cov.hiv[3,2,] <- -s2.vent
          cov.hiv[1,1,] <- s2.hosp * cumsum(delta.hosp.noicu.odesim)[t.idx] + 
                              s2.icu + .Machine$double.eps[1]
          cov.hiv[2,2,] <- s2.hosp * cumsum(delta.icu.novent.odesim)[t.idx] + 
                              s2.icu + s2.vent + .Machine$double.eps[1]
          cov.hiv[3,3,] <- s2.hosp * cumsum(delta.vent.odesim)[t.idx] + 
                              s2.vent + .Machine$double.eps[1]
          
          # change from model for hosp (without icu) and icu (without vent) and vent to
          #  a model for hosp (with icu and vent) and icu (with vent) and icu
          AA <- matrix(c(
            1, 1, 1,
            0, 1, 1,
            0, 0, 1
          ), nrow = 3, byrow = TRUE)

          mu.hiv <- t(AA %*% t(mu.hiv))

          for (ii in 1:nrow(mu.hiv)) {
            cov.hiv[, , ii] <- AA %*% cov.hiv[, , ii] %*% t(AA)
            cov.hiv[, , ii] <- 0.5 * (cov.hiv[, , ii] + t(cov.hiv[, , ii]))
          }
              
          if (loc == "RI") { 
            # likelihood for RI (all data present) 
            llvals.hosp.icu.vent <- rep(NA, nrow(mu.hiv))
                  
            for (ii in 1:nrow(mu.hiv)) {
              llvals.hosp.icu.vent[ii] <- dmvnorm(c(
                delta.hosp[ii],
                delta.icu[ii],
                delta.vent[ii]
              ),
              mean = mu.hiv[ii, ],
              sigma = cov.hiv[, , ii],
              log = TRUE
              )

              ll.hosp <- sum(llvals.hosp.icu.vent)
              ll.icu <- integer()
              ll.vent <- integer()

              # add to current likelihood
              ll <- ll + ll.hosp
            }
            
          } else { 
          
            ### likelihood for other states with lots of missing data ###
                  
            # find observed times
            t.obs.vent <- which(!is.na(tot.vent.curr))
            t.obs.icu <- which(!is.na(tot.icu.curr))
            t.obs.hosp <- which(!is.na(tot.hosp.curr))
                  
            ### vent likelihood
            if (lik.vent.curr) {
              delta.vent <- c(0, tot.vent.curr[t.obs.vent])
              delta.vent <- delta.vent[-1] - delta.vent[-length(delta.vent)]
              ll.vent <- dnorm(delta.vent[1],
                sum(mu.hiv[1:t.obs.vent[1], 3]),
                sqrt(sum(cov.hiv[3, 3, 1:t.obs.vent[1]])),
                log = TRUE
              )

              for (kk in 2:length(t.obs.vent)) {
                ll.vent <- ll.vent +
                  dnorm(delta.vent[kk],
                    sum(mu.hiv[(t.obs.vent[kk - 1] + 1):t.obs.vent[kk], 3]),
                    sqrt(sum(cov.hiv[3, 3, (t.obs.vent[kk - 1] + 1):t.obs.vent[1]])),
                    log = TRUE
                  )
              }

              if (length(ll.hosp) > 0) {
                ll.hosp <- ll.hosp + ll.vent
              } else {
                ll.hosp <- ll.vent
              }
            } else {
              ll.vent <- integer()
            }
                  
            ### icu likelihood
            if (lik.icu.curr) {
              delta.icu <- c(0, tot.icu.curr[t.obs.icu])
              delta.icu <- delta.icu[-1] - delta.icu[-length(delta.icu)]
              ll.icu <- dnorm(delta.icu[1],
                sum(mu.hiv[1:t.obs.icu[1], 3]),
                sqrt(sum(cov.hiv[3, 3, 1:t.obs.icu[1]])),
                log = TRUE
              )

              for (kk in 2:length(t.obs.icu)) {
                ll.icu <- ll.icu +
                  dnorm(delta.icu[kk],
                    sum(mu.hiv[(t.obs.icu[kk - 1] + 1):t.obs.icu[kk], 3]),
                    sqrt(sum(cov.hiv[3, 3, (t.obs.icu[kk - 1] + 1):t.obs.icu[1]])),
                    log = TRUE
                  )
              }

              if (length(ll.hosp) > 0) {
                ll.hosp <- ll.hosp + ll.icu
              } else {
                ll.hosp <- ll.icu
              }
            } else {
              ll.icu <- integer()
            }
                  
            ### hosp likelihood
            if (lik.hosp.curr) {
              delta.hosp <- c(0, tot.hosp.curr[t.obs.hosp])
              delta.hosp <- delta.hosp[-1] - delta.hosp[-length(delta.hosp)]
              ll.hosp.curr <- dnorm(delta.hosp[1],
                sum(mu.hiv[1:t.obs.hosp[1], 3]),
                sqrt(sum(cov.hiv[3, 3, 1:t.obs.hosp[1]])),
                log = TRUE
              )

              for (kk in 2:length(t.obs.hosp)) {
                ll.hosp.curr <- ll.hosp.curr +
                  dnorm(delta.hosp[kk],
                    sum(mu.hiv[(t.obs.hosp[kk - 1] + 1):t.obs.hosp[kk], 3]),
                    sqrt(sum(cov.hiv[3, 3, (t.obs.hosp[kk - 1] + 1):t.obs.hosp[1]])),
                    log = TRUE
                  )
              }

              if (length(ll.hosp) > 0) {
                ll.hosp <- ll.hosp + ll.hosp.curr
              } else {
                ll.hosp <- ll.hosp.curr
              }
            } else {
              ll.hosp.curr <- integer()
            }
                  
            if (length(ll.hosp) > 0) {
              ll <- ll + ll.hosp
            }
          }
        }
      }
    }
    
    #######################################################################
    ### 9. likelihood for new hospital discharges
    #######################################################################
    
    ll.hosp.dis <- NULL
    
    if (lik.hosp.discharges) {
      discharge.rates <- tot.hosp.discharges.new.odesim
      # finding any na.values
      no.na <- which(!is.na(tot.hosp.discharges.new))
      # log-likelihood given trajectory sim:
      ll.hosp.dis <- sum(dnbinom(tot.hosp.discharges.new[no.na],
        mu = discharge.rates[t.idx[no.na]],
        size = nb.disp.hosp.discharges,
        log = TRUE
      ))
      ll <- ll + ll.hosp.dis
    }
  } 
  
  ### combine likelihood results into a list
  list(
    ll = ll,
    ll.new = ll.new,
    ll.hosp.new = ll.hosp.new,
    ll.hosp = ll.hosp,
    ll.vent = ll.vent,
    ll.icu = ll.icu,
    sympt.imputed = sympt.imputed,
    ll.deaths = ll.deaths,
    ll.hosp.dis = ll.hosp.dis
  )
}
