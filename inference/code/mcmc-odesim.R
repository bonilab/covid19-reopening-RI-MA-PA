#!/usr/bin/env Rscript

# mcmc-odesim.R
# authors: Ephraim Hanks, Nathan Wikle, Emily Strong
# last edited: 20 Nov 2020 
#
# This file defines the mcmc.odesim function, which generates
#   samples from the posterior of all model parameters using
#   MCMC. A large number of function arguments are available
#   (see function definition); output is returned as a large list.

mcmc.odesim <- function(
  n.mcmc,                       # number of MCMC iterations to run
  df,                           # state-level covid data (must be a data frame) 
  odepath,                      # path to "odesim"
  odesim.ver = "v5",            # version of odesim (defaults to "v5")
  lik.tot = TRUE,               # evaluate likelihood for total new cases ONLY (ie, set F and use lik.age for age-strat. data; default = TRUE)
  lik.age = FALSE,              # evaluate likelihood for age-struc. new cases and hosps AND total new cases and hosps (default = FALSE)
  lik.hosp.new = TRUE,          # evaluate likelihood for new hosp. cases (default = FALSE)
  lik.hosp.curr = FALSE,        # evaluate likelihood for current hosp. (default = FALSE)
  lik.icu.curr = FALSE,         # evaluate likelihood for current icu admits (default = FALSE)
  lik.vent.curr = FALSE,        # evaluate likelihood for current vent admits (default = FALSE)
  lik.tot.deaths = FALSE,       # evaluate likelihood for tot. deaths ONLY (ie, set F and use lik.age.deaths for age-strat. data; default = FALSE)
  lik.age.deaths = FALSE,       # evaluate likelihood for age-struc. new deaths and total new deaths (default = FALSE)
  lik.home.deaths = FALSE,      # evaluate likelihood for new home deaths (default = FALSE)
  lik.hosp.deaths = FALSE,      # evaluate likelihood for new hosp. deaths (default = FALSE)
  lik.hosp.discharges = FALSE,  # evaluate likelihood for hospital discharges (default = FALSE)
  case.constraint = FALSE,      # constrain fit to cumulative cases to be within 10% of data (default = FALSE)
  active.surv = FALSE,          # include active surveillance data (default = FALSE)
  p.asympt = 0.4,               # proportion of asymptomatic individuals (default = 0.4; CAUTION: DO NOT CHANGE UNLESS ODESIM REFLECTS A DIFFERENT VALUE!)
  beta.start,                   # starting values for contact rate spline loadings
  spline.beta,                  # fda "spline" object spanning the time window
  report.rate.params.start,     # rate at which infecteds report (vector or scalar)
  spline.rr,                    # spline object (iSpline, bSpline, or constant vec) spanning time window
  ode.params.start = NULL,      # odesim param starting values (named vector matching odesim CLOs) 
  const.params = NULL,          # odesim CLOs, to be kept constant
  ode.params.prior.min = -Inf,  # vector of lower bounds for odesim params (uniform priors)
  ode.params.prior.max = Inf,   # vector of upper bounds for odesim params (uniform priors)
  non.odesim.params = NULL,     # non-odesim param starting values (eg, hosp.report.rate)
  lik.params.start = NULL,      # starting values for dispersion parameters in NB likelihood
  fixed.nb.disp = FALSE,        # Boolean indicating if NB dispersion params should be fixed (default = FALSE)
  start.day = 61,               # start day of odesim (default = 61, DON'T CHANGE)
  end.day,                      # end day of odesim
  introday = NULL,              # day of first infected
  loc = "RI",                   # US state used for analysis (one of "RI" (default), "MA", or "PA")
  s2.hosp.start = .01,          # initial value for current hosp variance hyperparam
  s2.icu.start = .01,           # initial value for current icu variance hyperparam
  s2.vent.start = .01,          # initial value for current vent variance hyperparam
  s2.beta.start = .01,          # initial value for beta prior (multivar. normal) marg. variance hyperparam
  s2.rr.start = .01,            # initial value for rr prior (multivar. normal) marg. variance hyperparam
  adapt.iter = 100,             # adaptive tuning update interval (log-adaptive tuning on MH based on Shaby and Wells, 2011)
  indep.mh = FALSE,             # if TRUE, propose beta separate from other params
  t.adapt.start = 0,            # number of times adaptive tuning var. has been updated (default = 0)
  prop.type = "tnorm",          # MH proposal type (default = "tnorm")
  adapt.type = "ShabyWells",    # adaptive tuning type (default = "ShabyWells")
  c0 = 1,                       # Shaby Wells adaptive tuning constant c0 
  c1 = 0.8,                     # Shaby Wells adaptive tuning constant c1
  var.tune = NULL,              # list of tuning param variances, order = (beta, ode.params, rr, s2, lik)
  Sigma.tune = NULL,            # list of tuning param covars, order = (beta, ode.params, rr, s2, lik)
  p.vecs,                       # weekly vector of delay probabilities (should be c(1, rep(0,6)) unless good reason otherwise)
  thin = 1,                     # rate to thin saved posterior samples (default = 1)
  plot.save = TRUE,             # plot trace plots while mcmc is running (default = TRUE)
  plot.rate = 10,               # refresh rate on trace plots (default = 10)
  plot.name = "traceplots.pdf", # name of trace plots (default = "traceplots.pdf")
  print.iter = FALSE,           # print iteration number everytime sample is saved (default = FALSE)
  sf.choice = FALSE,            # estimate proportion of symptomatics from 7 possible age-strat. combinations (default = FALSE)
) {

  ########################
  ### 1. Preliminaries ###
  ########################
  
  # adaptive structure
  t.adapt <- t.adapt.start

  # process data frame for likelihood evaluation
  data.processed <- data.process(df, loc = loc)
  
  # add data.processed results to global environment
  list2env(data.processed, globalenv())
    
  # useful constants
  n.days <- length(days)
  num.beta <- length(beta.start)
  num.rr <- length(report.rate.params.start)
  num.ode.params <- length(ode.params.start)
  num.lik.params <- length(lik.params.start)
    
  ### structures to save parameter samples:
  beta.save <- matrix(NA_real_, nrow = n.mcmc/thin, ncol = num.beta)
  rr.params.save <- matrix(NA_real_, nrow = n.mcmc/thin, ncol = num.rr)
  ode.params.save <- matrix(NA_real_, nrow = n.mcmc/thin, ncol = num.ode.params)
  ode.params.names <- names(ode.params.start)
  colnames(ode.params.save) <- ode.params.names
  
  if (num.lik.params > 0){
    lik.params.save <- matrix(NA_real_, nrow = n.mcmc/thin, ncol = num.lik.params)
  } else {
    lik.params.save <- NULL
  }
  
  s2.beta.save <- rep(NA_real_, n.mcmc/thin)
  s2.rr.save <- rep(NA_real_, n.mcmc/thin)
  s2.params.save <- matrix(NA_real_, n.mcmc/thin, 3)
  loglik.save <- rep(NA_real_, n.mcmc/thin)

  ### initialize parameters:
  
  # betas
  beta <- beta.start # vector of beta parameters
  Z <- eval.basis(spline.beta, start.day:end.day) # beta spline basis functions 
  beta.daily <- Z %*% beta
    
  # reporting rate
  if (is.matrix(spline.rr)){
    Z.rr <- spline.rr
  } else {
    Z.rr <- eval.basis(spline.rr, start.day:end.day)
  }
    
  rr.params <- report.rate.params.start
  rr.daily <- Z.rr %*% rr.params
  
  # hyperparameters
  s2.hosp <- s2.hosp.start
  s2.icu <- s2.icu.start
  s2.vent <- s2.vent.start
  s2.params <- c(s2.hosp, s2.icu, s2.vent)
  s2.beta <- s2.beta.start # marginal variance of beta prior
  s2.rr <- s2.rr.start
  ode.params <- ode.params.start
  lik.params <- lik.params.start
  lik.params.star <- lik.params ## needs to be NULL if lik.params = NULL so that loglik.odesim gets right values

  # extra parameters (like hospitalization reporting rate )
  extra.params <- NULL
  extra.const.params <- NULL
  extra.params.fitted.idx <- integer()
  extra.params.const.idx <- integer()
    
  if(length(non.odesim.params)>0){
    
    for(k in 1:length(non.odesim.params)){
      extra.params.fitted.idx <- c(extra.params.fitted.idx,
        which(names(ode.params) == non.odesim.params[k]))
    }
      
    if (length(extra.params.fitted.idx) > 0){
      extra.params <- ode.params[extra.params.fitted.idx]
    }
    
    for(k in 1:length(non.odesim.params)){
      extra.params.const.idx <- c(extra.params.const.idx, 
        which(names(const.params) == non.odesim.params[k]))
    }
    
    if (length(extra.params.const.idx) > 0){
      extra.const.params <- const.params[extra.params.const.idx]
    }
  }  
  
  ### parameter for symptomatic fractions
  
  if (sf.choice){
    # if true, create structure for sampling symp-frac settings (6 options)
    
    # list of model possibilities
    symp.vals <- c("", "-symp-frac-davies", "-symp-frac-equal 0.3",
                   "-symp-frac-equal 0.4", "-symp-frac-equal 0.5",
                   "-symp-frac-equal 0.6", "-symp-frac-equal 0.7")
    
    # number of symp-frac options
    n.sf <- length(symp.vals)
    
    # initialize to random starting symp-frac option
    K <- sample.int(n.sf, size = 1)
    
    # structure to store models
    K.vals <- rep(NA_integer_, n.mcmc/thin)
    symp.cur <- symp.vals[K]
    
  } else {
    symp.cur <- NULL
    K.vals <- NULL
  }
  
    
  ### create a matrix P which can be used to calculate Poisson rates after delay
  ###   note: this is no longer used for any state
  P <- matrix(0, nrow = end.day + max(sapply(p.vecs, length)), ncol = end.day+1)
  colnames(P) <- paste("Day", 1:(end.day+1), sep = " ")
  
  for(j in 1:ncol(P)){
    if (j < 74){
      ## period 1: beginning - March 13
      P[j:(j + length(p.vecs[[1]]) - 1), j] <- p.vecs[[1]]
    } else if ((j >= 74) & (j < 84)){
      ## period 2: March 14 - March 23
      P[j:(j + length(p.vecs[[2]]) - 1), j] <- p.vecs[[2]]
    } else if ((j >= 84) & (j < 88)){
      ## period 3: March 24 - March 27
      P[j:(j + length(p.vecs[[3]]) - 1), j] <- p.vecs[[3]]
    } else {
      ## period 4: March 28 - present
      P[j:(j + length(p.vecs[[4]]) - 1), j] <- p.vecs[[4]]
    }
  }


  #################
  ### 2. Priors ###
  #################

  ### beta prior: random walk
  ###   (penalized regression spline w/ 1st-order diffs)
  
  D <- diff(diag(num.beta), differences = 1)
  S <- crossprod(D)

  ### beta ~ N(0, s2.beta * S^-1)
  beta.prior.loglik <- function(beta, s2.beta, precision = S){
    if(min(beta) < 0){
      ll <- -Inf
    } else {
      ll <- -1 / 2 / s2.beta * (t(beta) %*% (precision %*% beta))
    }
    return(ll)
  }
  
  ### rr prior: random walk
  ###   (same prior for report.rate params as for beta...)
  
  D.rr <- diff(diag(num.rr), differences = 1)
  S.rr <- crossprod(D.rr)

  ### rr ~ N(0, s2.rr * S.rr^-1)
  rr.prior.loglik <- function(rr, s2.rr, precision = S.rr){
    if(min(rr) < 0 | max(rr) > 1){
      ll <- -Inf
    } else {
      ll <- -1 / 2 / s2.rr * (t(rr) %*% (precision %*% rr))
    }
    return(ll)
  }


  ### dispersion parameter: exponential prior
  ###   disp ~ Exp(lambda = 100)
  lik.params.prior.loglik <- function(lik.params, 
                                      lambda = rep(100, length(lik.params))){
    if(length(lik.params) < 1){
      ll <- 0
    } else {
      # exponential prior
      ll <- sum(dexp(lik.params, rate = lambda, log = TRUE))
      # improper uniform prior
      #   ll=0
      #   if(min(lik.params)<.1){
      #     ll=-Inf
      #   }
    }
    return(ll)
  }

  ### s2.beta prior: inverse gamma
  ###   s2.beta ~ IG(s = shape, r = rate)
  
  # initialize: s = 1, r = 1
  s <- 1; r <- 1

  # Uniform priors for odesim parameters
  ode.params.prior.loglik <- function(ode.params, ode.params.min, ode.params.max){
    ll <- 0
    if(sum(c(ode.params < ode.params.min , ode.params > ode.params.max)) > 0){
      ll <- -Inf
    }
    return(ll)
  }

  ### Uniform priors for s2 parameters
  s2.params.prior.loglik <- function(s2.params){
    ll <- 0
    if(min(s2.params) < 0){
      ll <- -Inf
    }
    return(ll)
  }

  ################################
  ### 3. Likelihood evaluation ###
  ################################

  ### simulate trajectory using current beta values:
  traj <- traj.from.params(beta = beta.daily, 
                           params = ode.params, 
                           tf = end.day, 
                           introday = introday,
                           const.params = const.params,
                           non.odesim.params = non.odesim.params,
                           odepath = odepath,
                           loc = loc,
                           symp = symp.cur)
    
  ### evaluate logliklihood under initial conditions
  llvals <- loglik.odesim(traj = traj,
                          df = df,
                          dp = data.processed,
                          odesim.ver = odesim.ver,
                          P = P,
                          loc = loc,
                          report.rate = rr.daily,
                          nb.disp.params = lik.params,
                          lik.tot = lik.tot,
                          lik.age = lik.age,
                          lik.hosp.new = lik.hosp.new,
                          lik.hosp.curr = lik.hosp.curr,
                          lik.icu.curr = lik.icu.curr,
                          lik.vent.curr = lik.vent.curr,
                          lik.tot.deaths = lik.tot.deaths,
                          lik.home.deaths = lik.home.deaths,
                          lik.hosp.deaths = lik.hosp.deaths,
                          lik.age.deaths = lik.age.deaths,
                          lik.hosp.discharges = lik.hosp.discharges,
                          active.surv = active.surv,
                          p.asympt = p.asympt,
                          case.constraint = case.constraint,
                          s2.hosp = s2.hosp,
                          s2.icu = s2.icu,
                          s2.vent = s2.vent,
                          extra.params = extra.params, ### if cumul. hosp. reporting rate is fitted 
                          extra.const.params = extra.const.params) ### if cumul. hosp. reporting rate is constant
  ll.current <- llvals$ll
  ll.new <- llvals$ll.new
  ll.hosp.new <- llvals$ll.new
  
  ###################################
  ### 4. Adaptive proposal set-up ###
  ###################################
  accept <- rep(0, length(var.tune))
  accept.tot <- accept
  never.adapt <- TRUE

  ##########################
  ### 5. MCMC iterations ###
  ##########################
  for(iter in 1:n.mcmc){

    # print-out iteration number every 100 iterations
    if (iter %% 100 == 0 & print.iter){
      cat(iter, " ")
    }

    #########################################################################
    ### Propose beta
    #########################################################################

    ### setting up error catching - reject beta if odesim fails
    beta.good <- 0
    
    while(beta.good < 1){
      # if(prop.type=="norm"){
      #     beta.star <- t(rmvnorm(1, c(beta), var.tune[1] * Sigma.tune[[1]]))
      # }
      # if(prop.type=="tnorm"){
      #     beta.star <- rtnorm(length(beta),beta,sqrt(var.tune[1]*diag(Sigma.tune[[1]])),0,Inf)
      # }
      
      # propose beta.star
      beta.star <- exp(rnorm(length(beta), log(beta), 
                             sqrt(var.tune[1] * diag(Sigma.tune[[1]]))))
      beta.daily.star <- Z %*% beta.star
      
      ### trajectory given beta.star
      traj.star <- try(traj.from.params(beta = beta.daily.star, 
                                        params = ode.params,
                                        tf = end.day,
                                        introday = introday,
                                        const.params = const.params,
                                        non.odesim.params = non.odesim.params,
                                        odepath = odepath,
                                        loc = loc,
                                        symp = symp.cur), silent = TRUE)
      
      # check if error was thrown. if not, leave while loop
      if(class(traj.star) != "try-error"){
        beta.good <- 1
      }
    }

    ### evaluate loglikelihood for beta.star
    llvals.star <- loglik.odesim(traj = traj.star,
                                 dp = data.processed,
                                 report.rate = rr.daily,
                                 nb.disp.params = lik.params,
                                 s2.hosp = s2.params[1],
                                 s2.icu = s2.params[2],
                                 s2.vent = s2.params[3],
                                 lik.tot = lik.tot,
                                 lik.age = lik.age,
                                 lik.hosp.new = lik.hosp.new,
                                 lik.hosp.curr = lik.hosp.curr,
                                 lik.icu.curr = lik.icu.curr, 
                                 lik.vent.curr = lik.vent.curr,
                                 lik.tot.deaths = lik.tot.deaths,
                                 lik.home.deaths = lik.home.deaths,
                                 lik.hosp.deaths = lik.hosp.deaths,
                                 lik.age.deaths = lik.age.deaths,
                                 lik.hosp.discharges = lik.hosp.discharges,
                                 active.surv = active.surv,
                                 p.asympt = p.asympt,
                                 case.constraint = case.constraint,
                                 df = df,
                                 odesim.ver = odesim.ver,
                                 P = P,
                                 loc = loc,
                                 extra.params = extra.params,  
                                 extra.const.params = extra.const.params) 
    ll.star <- llvals.star$ll
        
    ### accept/reject beta.star
        
    mh1 <- ll.star + beta.prior.loglik(beta.star, s2.beta) + sum(log(beta.star))
    mh2 <- ll.current + beta.prior.loglik(beta, s2.beta) + sum(log(beta))
    
    # if(prop.type=="tnorm"){
    #     mh1=mh1+sum(dtnorm(beta,beta.star,sqrt(var.tune[1]*diag(Sigma.tune[[1]])),0,Inf,log=TRUE))
    #     mh2=mh2+sum(dtnorm(beta.star,beta,sqrt(var.tune[1]*diag(Sigma.tune[[1]])),0,Inf,log=TRUE))
    # }
    
    if(is.na(mh1)){
      mh1 <- -Inf
    }

    ####if Unif(0,1) < mh1/mh2, accept new beta and disp parameters
    if (exp(mh1 - mh2) > runif(1)){
      
      ### accept beta.star
      beta <- beta.star
      beta.daily <- beta.daily.star
      traj <- traj.star
      llvals=llvals.star
      ll.current <- ll.star
      ## if(lik.age){
      ##     sympt.new.imputed=sympt.new.star
      ## }
      accept[1] <- accept[1] + 1
    }


    #########################################################################
    ### Propose ode.params
    #########################################################################

    ### setting up error catching - reject ode.params if odesim fails
    beta.good <- 0
    while(beta.good < 1){
      
      if(prop.type == "norm"){
        ode.params.star <- t(rmvnorm(1, ode.params, var.tune[2] * Sigma.tune[[2]]))
      }
      if(prop.type == "tnorm"){
        ode.params.star <- (rtnorm(length(ode.params), ode.params, 
          sqrt(var.tune[2] * diag(Sigma.tune[[2]])), ode.params.prior.min, ode.params.prior.max))
      }
      
      names(ode.params.star) <- ode.params.names
      extra.params.star <- ode.params.star[extra.params.fitted.idx]
      
      ### trajectory given ode.params.star
      traj.star <- try(traj.from.params(beta = beta.daily, 
                                        params = ode.params.star,
                                        tf = end.day, 
                                        introday = introday,
                                        const.params = const.params,
                                        non.odesim.params = non.odesim.params,
                                        odepath = odepath,
                                        loc = loc,
                                        symp = symp.cur), silent = TRUE)
      
      if(class(traj.star) != "try-error" & min(ode.params.star > 0)){
        beta.good <- 1
      }
    }

    ### evalulate loglikelihood for ode.params.star
    llvals.star <- loglik.odesim(traj = traj.star,
                                 dp = data.processed,
                                 report.rate = rr.daily,
                                 nb.disp.params = lik.params,
                                 s2.hosp = s2.params[1],
                                 s2.icu = s2.params[2],
                                 s2.vent = s2.params[3],
                                 lik.tot = lik.tot,
                                 lik.age = lik.age,
                                 lik.hosp.new = lik.hosp.new,
                                 lik.hosp.curr = lik.hosp.curr,
                                 lik.icu.curr = lik.icu.curr, 
                                 lik.vent.curr = lik.vent.curr,
                                 lik.tot.deaths = lik.tot.deaths,
                                 lik.home.deaths = lik.home.deaths,
                                 lik.hosp.deaths = lik.hosp.deaths,
                                 lik.age.deaths = lik.age.deaths,
                                 lik.hosp.discharges = lik.hosp.discharges,
                                 active.surv = active.surv,
                                 p.asympt = p.asympt,
                                 case.constraint = case.constraint,
                                 df = df,
                                 odesim.ver = odesim.ver,
                                 P = P,
                                 loc = loc,
                                 extra.params = extra.params.star, 
                                 extra.const.params = extra.const.params) 
    ll.star <- llvals.star$ll
    
    ### accept/reject ode.params.star

    mh1 <- ll.star + ode.params.prior.loglik(ode.params.star,
                                             ode.params.prior.min,
                                             ode.params.prior.max)
    mh2 <- ll.current +  ode.params.prior.loglik(ode.params,
                                                 ode.params.prior.min,
                                                 ode.params.prior.max) 
    if(prop.type == "tnorm"){
      mh1 <- mh1 + sum(dtnorm(ode.params,
                              ode.params.star,
                              sqrt(var.tune[2] * diag(Sigma.tune[[2]])),
                              ode.params.prior.min,
                              ode.params.prior.max, 
                              log = TRUE))
      
      mh2 <- mh2 + sum(dtnorm(ode.params.star,
                              ode.params,
                              sqrt(var.tune[2] * diag(Sigma.tune[[2]])),
                              ode.params.prior.min,
                              ode.params.prior.max,
                              log = TRUE))
    }
    
    if(is.na(mh1)){
      mh1 <- -Inf
    }

    ### if Unif(0,1) < mh1/mh2, accept new odesim params
    if (exp(mh1 - mh2) > runif(1)){
      ode.params <- ode.params.star
      extra.params <- extra.params.star
      traj <- traj.star
      llvals <- llvals.star
      ll.current <- ll.star
      # if(lik.age){
      #     sympt.new.imputed=sympt.new.star
      # }
      accept[2] <- accept[2] + 1
    }
    
    
    #########################################################################
    ### Propose symp-frac option (if sf.choice == TRUE)
    #########################################################################
    
    if (sf.choice){
      
      # propose model from uniform prior
      K.star <- sample(n.sf, size = 1)
      symp.star <- symp.vals[K.star]
      
      
      # calculate traj and likelihoods for each possible model:
      traj.vals <- list()
      
      good.k <- TRUE
      
        
      ## trajectory given beta.star
      traj.star <- try(traj.from.params(beta = beta.daily, 
                                        params = ode.params,
                                        tf = end.day,
                                        introday = introday,
                                        const.params = const.params,
                                        non.odesim.params = non.odesim.params,
                                        odepath = odepath,
                                        loc = loc,
                                        symp = symp.star),
                       silent = TRUE)
      if(class(traj.star) == "try-error"){
        good.k <- FALSE
      } 
      
      if (good.k){
        
        llvals.star <- loglik.odesim(traj = traj.star,
                                 dp = data.processed,
                                 report.rate = rr.daily,
                                 nb.disp.params = lik.params,
                                 s2.hosp = s2.params[1],
                                 s2.icu = s2.params[2],
                                 s2.vent = s2.params[3],
                                 lik.tot = lik.tot,
                                 lik.age = lik.age,
                                 lik.hosp.new = lik.hosp.new,
                                 lik.hosp.curr = lik.hosp.curr,
                                 lik.icu.curr = lik.icu.curr, 
                                 lik.vent.curr = lik.vent.curr,
                                 lik.tot.deaths = lik.tot.deaths,
                                 lik.home.deaths = lik.home.deaths,
                                 lik.hosp.deaths = lik.hosp.deaths,
                                 lik.age.deaths = lik.age.deaths,
                                 lik.hosp.discharges = lik.hosp.discharges,
                                 active.surv = active.surv,
                                 p.asympt = p.asympt,
                                 case.constraint = case.constraint,
                                 df = df,
                                 odesim.ver = odesim.ver,
                                 P = P,
                                 loc = loc,
                                 extra.params = extra.params, 
                                 extra.const.params = extra.const.params) 
        ll.star <- llvals.star$ll
        
        ### accept/reject new symp frac option
        mh1 <- ll.star # uniform prior and proposal...
        mh2 <- ll.current # uniform prior and proposal...
        
        if(is.na(mh1)){
          mh1 <- -Inf
        }
        
        ### if Unif(0,1) < mh1/mh2, accept symp.frac.star
        if (exp(mh1 - mh2) > runif(1)){
          K <- K.star 
          symp.cur <- symp.star
          llvals <- llvals.star
          ll.current <- ll.star
          traj <- traj.star
        }
      }
    }

    #########################################################################
    ### Propose rr.params
    #########################################################################

    # sample new rr params
    rr.params.star <- rtnorm(length(rr.params), rr.params, 
                            sqrt(diag(var.tune[3] * Sigma.tune[[3]])), 0, 1)
    rr.daily.star <- Z.rr %*% rr.params.star
        
        
    if (max(rr.daily.star, na.rm = T) < 1){
      # make sure rr < 1  
          
      ### evaluate logliklihood with rr.params.star
      llvals.star <- loglik.odesim(traj = traj,
                                   dp = data.processed,
                                   report.rate = rr.daily.star,
                                   nb.disp.params = lik.params,
                                   s2.hosp = s2.params[1],
                                   s2.icu = s2.params[2],
                                   s2.vent = s2.params[3],
                                   lik.tot = lik.tot,
                                   lik.age = lik.age,
                                   lik.hosp.new = lik.hosp.new,
                                   lik.hosp.curr = lik.hosp.curr,
                                   lik.icu.curr = lik.icu.curr, 
                                   lik.vent.curr = lik.vent.curr,
                                   lik.tot.deaths = lik.tot.deaths,
                                   lik.home.deaths = lik.home.deaths, 
                                   lik.hosp.deaths = lik.hosp.deaths,
                                   lik.age.deaths = lik.age.deaths,
                                   lik.hosp.discharges = lik.hosp.discharges,
                                   active.surv = active.surv,
                                   p.asympt = p.asympt,
                                   case.constraint = case.constraint,
                                   df = df,
                                   odesim.ver = odesim.ver,
                                   P = P,
                                   loc = loc,
                                   extra.params = extra.params, 
                                   extra.const.params = extra.const.params) 
      ll.star <- llvals.star$ll
          
      ### accept/reject rr.params.star
      
      mh1 <- ll.star + rr.prior.loglik(rr.params.star, 
                                       s2.rr, 
                                       precision = S.rr) 
          
      mh2 <- ll.current + rr.prior.loglik(rr.params, 
                                          s2.rr, 
                                          precision = S.rr)
        
          
      if(prop.type == "tnorm"){
        mh1 <- mh1 + sum(dtnorm(rr.params,
                                rr.params.star,
                                sqrt(var.tune[3] * diag(Sigma.tune[[3]])),
                                0, 1, log = TRUE))
        
        mh2 <- mh2 + sum(dtnorm(rr.params.star,
                                rr.params,
                                sqrt(var.tune[3] * diag(Sigma.tune[[3]])),
                                0, 1, log = TRUE))
      }
          
      if(is.na(mh1)){
        mh1 <- -Inf
      }
          
      ### if Unif(0,1) < mh1/mh2, accept rr.params
      if (exp(mh1 - mh2) > runif(1)){
        rr.params <- rr.params.star
        rr.daily <- rr.daily.star
        llvals <- llvals.star
        ll.current <- ll.star
        # if(lik.age){
        #     sympt.new.imputed=sympt.new.star
        # }
        accept[3] <- accept[3] + 1
      }
    }

        

    #########################################################################
    ### Propose s2.params
    #########################################################################      

    if(!lik.hosp.curr & !lik.vent.curr & !lik.icu.curr){
      # no need to propose new s2 parameters
    } else {
        
      # s2.params.star <- t(rtnorm(length(s2.params), s2.params, sqrt(diag(var.tune[4] * Sigma.tune[[4]])),0,Inf))
      s2.params.star <- exp(t(rnorm(length(s2.params), 
                                    log(s2.params), 
                                    sqrt(diag(var.tune[4] * Sigma.tune[[4]])))))
      
      ### evaluate logliklihood for s2.params.star
      llvals.star <- loglik.odesim(traj = traj,
                                   dp = data.processed,
                                   report.rate = rr.daily,
                                   nb.disp.params = lik.params,
                                   s2.hosp = s2.params.star[1],
                                   s2.icu = s2.params.star[2],
                                   s2.vent = s2.params.star[3],
                                   lik.tot = lik.tot,
                                   lik.age = lik.age,
                                   lik.hosp.new = lik.hosp.new,
                                   lik.hosp.curr = lik.hosp.curr,
                                   lik.icu.curr = lik.icu.curr, 
                                   lik.vent.curr = lik.vent.curr,
                                   lik.tot.deaths = lik.tot.deaths,
                                   lik.home.deaths = lik.home.deaths,
                                   lik.hosp.deaths = lik.hosp.deaths,
                                   lik.age.deaths = lik.age.deaths,
                                   lik.hosp.discharges = lik.hosp.discharges,
                                   active.surv = active.surv,
                                   p.asympt = p.asympt,
                                   case.constraint = case.constraint,
                                   df = df,
                                   odesim.ver = odesim.ver,
                                   P = P,
                                   loc = loc,
                                   extra.params = extra.params, 
                                   extra.const.params = extra.const.params) 
      ll.star <- llvals.star$ll

      ### accept/reject s2.params.star
        
      # determine which s2 values to include or not
      s2.true <- rep(FALSE, 3)
          
      mh1 <- 0; mh2 <- 0
      
      if(length(llvals.star$ll.hosp) > 0){
        mh1 <- mh1 + llvals.star$ll.hosp
        mh2 <- mh2 + llvals$ll.hosp
        s2.true[1] <- TRUE
      } 
      
      if(length(llvals.star$ll.icu) > 0){
        mh1 <- mh1 + llvals.star$ll.icu
        mh2 <- mh2 + llvals$ll.icu
        s2.true[2] <- TRUE
      } 
      
      if(length(llvals.star$ll.vent) > 0){
        mh1 <- mh1 + llvals.star$ll.vent
        mh2 <- mh2 + llvals$ll.vent
        s2.true[3] <- TRUE
      } 
        
      mh1 <- mh1 + s2.params.prior.loglik(s2.params.star[s2.true])
      mh2 <- mh2 + s2.params.prior.loglik(s2.params[s2.true])
          
      if(prop.type == "tnorm"){
        mh1 <- mh1 + sum(log(s2.params.star[s2.true]))
              #+sum(dtnorm(s2.params,s2.params.star,sqrt(var.tune[4]*diag(Sigma.tune[[4]])),0,Inf,log=TRUE))
        mh2 = mh2 + sum(log(s2.params[s2.true]))
            #+sum(dtnorm(s2.params.star,s2.params,sqrt(var.tune[4]*diag(Sigma.tune[[4]])),0,Inf,log=TRUE))
      }
          
      if(is.na(mh1)){
        mh1 <- -Inf
      }

      ### if Unif(0,1) < mh1/mh2, accept new s2.params
      if (exp(mh1 - mh2) > runif(1)){
        s2.params[s2.true] <- s2.params.star[s2.true]
        llvals <- llvals.star
        ll.current <- ll.star
        # if(lik.age){
        #     sympt.new.imputed=sympt.new.star
        # }
        accept[4] <- accept[4] + 1
      }
    }

    
    #########################################################################
    ### Propose lik.params (NB dispersion params)
    #########################################################################      

      if(!fixed.nb.disp){
          lik.params.star <- exp(t(rnorm(length(lik.params), 
                                         log(lik.params), 
                                         sqrt(diag(var.tune[5] * Sigma.tune[[5]])))))
          
          ### evaluate logliklihood with lik.params.star
          llvals.star <- loglik.odesim(traj = traj,
                                       dp = data.processed,
                                       report.rate = rr.daily,
                                       nb.disp.params = lik.params.star,
                                       s2.hosp = s2.params[1],
                                       s2.icu = s2.params[2],
                                       s2.vent = s2.params[3],
                                       lik.tot = lik.tot,
                                       lik.age = lik.age,
                                       lik.hosp.new = lik.hosp.new,
                                       lik.hosp.curr = lik.hosp.curr,
                                       lik.icu.curr = lik.icu.curr, 
                                       lik.vent.curr = lik.vent.curr,
                                       lik.tot.deaths = lik.tot.deaths,
                                       lik.home.deaths = lik.home.deaths,
                                       lik.hosp.deaths = lik.hosp.deaths,
                                       lik.age.deaths = lik.age.deaths,
                                       lik.hosp.discharges = lik.hosp.discharges,
                                       active.surv = active.surv,
                                       p.asympt = p.asympt,
                                       case.constraint = case.constraint,
                                       df = df,
                                       odesim.ver = odesim.ver,
                                       P = P,
                                       loc = loc,
                                       extra.params = extra.params, 
                                       extra.const.params=extra.const.params) 
          ll.star <- llvals.star$ll
          
          ## accept/reject lik.params.star
          
          mh1 <- ll.star +  
              lik.params.prior.loglik(lik.params.star) + 
              sum(log(lik.params.star))
          
          mh2 <- ll.current + 
              lik.params.prior.loglik(lik.params) + 
              sum(log(lik.params))
          
          if(is.na(mh1)){
              mh1 <- -Inf
          }
          
          ### if Unif(0,1) < mh1/mh2, accept new lik.params
          if (exp(mh1 - mh2) > runif(1)){
              lik.params <- lik.params.star
              llvals <- llvals.star
              ll.current <- ll.star
                                        # if(lik.age){
                                        #     sympt.new.imputed=sympt.new.star
                                        # }
              accept[5] <- accept[5] + 1
          }
      }


    #########################################################################
    ### Propose new s2.beta and s2.rr parameters
    #########################################################################          
      
    ### conjugate updates from inverse gamma marg. post.

    s2.beta <- 1/rgamma(1, shape = s + num.beta / 2,
                        rate = r + .5 * (t(beta) %*% (S %*% beta)))

    s2.rr <- 1/rgamma(1, shape = s + num.rr / 2,
                      rate = r + .5 * (t(rr.params) %*% (S.rr %*% rr.params)))

      
    #########################################################################
    ### Save parameters
    #########################################################################      

    if (iter %% thin == 0){
      
      beta.save[iter/thin,] <- beta
      ode.params.save[iter/thin,] <- ode.params
      rr.params.save[iter/thin,] <- rr.params
      s2.params.save[iter/thin,] <- s2.params
      lik.params.save[iter/thin,] <- lik.params
      s2.beta.save[iter/thin] <- s2.beta
      s2.rr.save[iter/thin] <- s2.rr
      loglik.save[iter/thin] <- ll.current
      # traj.save[iter/thin,,]=as.matrix(traj)
      
      if (sf.choice){
        K.vals[iter/thin] <- K
      }
    }

    
    #########################################################################
    ### Update adaptive tuning variances
    #########################################################################      

    if (iter %% adapt.iter == 0){
      
      never.adapt <- FALSE
      
      ### Using Shaby and Wells adaptive scheme for RWM...
      if(adapt.type == "ShabyWells"){
        # default constants
        #     c0 <- 1; c1 <- 0.8;
        r.opt <- 0.234
        r.hat <- accept / adapt.iter
        t.adapt <- iter/adapt.iter + t.adapt.start
        gamma.1 <- 1 / (t.adapt)^c1
        gamma.2 <- c0 * gamma.1
        var.tune <- exp(log(var.tune) + gamma.2 * (r.hat - r.opt))
                
        accept.tot <- accept.tot + accept
        accept <- rep(0, length(accept))
        cat("\n","current accept rate =", r.hat)
        cat("\n","new proposal var =", var.tune, "\n")
      }
      
    }

    #########################################################################
    ### plotting output to assess convergence
    #########################################################################    

    if(iter %% plot.rate == 0 & plot.save == TRUE){
      pdf(plot.name, width = 10, height = 7)
      matplot(beta.save[1:iter/thin,],
              type = "l",
              main = "beta")
      matplot(ode.params.save[1:iter/thin,],
              type = "l",
              main = "ode.params")
      matplot(rr.params.save[1:iter/thin,],
              type = "l",
              main = "report.rate")
      matplot(s2.params.save[1:iter/thin,],
              type = "l", 
              main = "s2.params")
      matplot(log(lik.params.save[1:iter/thin,]),
              type = "l",
              main = "log(NB dispersion params)")
      dp <- data.process(df, loc = loc)
      tp <- traj.process(traj, loc = loc, odesim.version = odesim.ver)
      plots.odesim(dp, tp, rr.daily[-length(rr.daily)])
      dev.off()
    }

  } ## end MCMC

  if(print.iter){
    cat("\n") # new line
  }

  
  #########################################
  ### 6. Adaptive structures for output ###
  #########################################

  accept.rate.final <- accept / n.mcmc
  Sigma.hat <- list(cov(beta.save),
                    cov(ode.params.save),
                    cov(rr.params.save),
                    cov(s2.params.save),
                    cov(lik.params.save))
    
  
  ############################
  ### 7. Save MCMC results ###
  ############################

  ### output MCMC samples in a list
  list(beta = beta.save, 
       rr.params = rr.params.save,
       ode.params = ode.params.save,
       lik.params = lik.params.save,
       s2.beta = s2.beta.save,
       s2.rr = s2.rr.save,
       s2.params = s2.params.save,
       traj.save = traj,
       sf.choice = sf.choice,
       sf.vals = K.vals,
       spline.beta = spline.beta,
       spline.rr = spline.rr,
       df = df,
       lik.tot = lik.tot,
       lik.age = lik.age,
       lik.hosp.new = lik.hosp.new,
       lik.hosp.curr = lik.hosp.curr,
       lik.icu.curr = lik.icu.curr,
       lik.vent.curr = lik.vent.curr,
       lik.tot.deaths = lik.tot.deaths,
       lik.home.deaths = lik.home.deaths,
       lik.hosp.deaths = lik.hosp.deaths,
       lik.age.deaths = lik.age.deaths,
       lik.hosp.discharges = lik.hosp.discharges,
       odesim.ver = odesim.ver,
       odepath = odepath,
       introday = introday,
       Sigma.hat = Sigma.hat,
       Sigma.tune = Sigma.tune,
       var.tune = var.tune,
       loglik.final.iter = ll.current,
       loglik.vals.final.iter = llvals,
       loglik = loglik.save,
       accept.rate = accept.rate.final,
       t.adapt.end = t.adapt,
       ode.params.prior.min = ode.params.prior.min, 
       ode.params.prior.max = ode.params.prior.max,  
       thin = thin,
       loc = loc,
       adapt.type = adapt.type,
       const.params = const.params,
       non.odesim.params = non.odesim.params,
       P=P,
       extra.params=extra.params,
       extra.const.params=extra.const.params,
       active.surv=active.surv,
       today = Sys.Date() )
}
