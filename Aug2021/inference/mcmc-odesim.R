#!/usr/bin/env Rscript

### mcmc-odesim.R
### last edited: 12 Nov 2021
### authors: Ephraim Hanks, Nathan Wikle


### Set up an MCMC sampler for our mechanistic model of SARS-CoV-2 in Rhode Island,
###   Massachusetts, or Pennsylvania. The output is a large list.
mcmc.odesim <- function(
  ### MCMC parameters ###
  n.mcmc,                         # number of MCMC iterations
  df,                             # observed data
  loc = "RI",                     # state of interest ("RI", "MA", or "PA")
  odepath,                        # path to "ODESIM" source code
  odesim.ver = "v5",              # odesim version
  start.day = 61,                 # start day for ODESIM
  end.day = 116,                  # end day for ODESIM
  introday = NULL,                # day of first infected (ODESIM)
  thin = 1,                       # rate to thin saved posterior samples
  print.iter = FALSE,             # print MCMC iterations (Boolean)
  plot.save = FALSE,              # iteratively save trace plots (Boolean)
  plot.rate = 10,                 # trace plot generation rate
  plot.name = "trace-plots.pdf",  # name of trace plot
  ### Likelihood specification ###
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
  p.vecs,               # weekly vector of delay probabilities
  pres.delay = 0,       # number of days from sypmtoms to clinical presenation         
  active.surv = FALSE,  # use active surveillance data (Boolean)
  p.asympt = .4,        # asymptomatic proportion (only used if active.surv = TRUE)
  total.size.constraint = FALSE,  # require estimated cumulative cases to be within 10% of observed data (Boolean)
  sf.choice = FALSE,              # estimate proportion of symptomatic infecteds (Boolean)
  ### Initialize parameters ###
  beta.start,               # initial values for spline loadings
  spline.beta,              # fda "spline" object spanning the time window in "days"
  report.rate.params.start, # initial reporting rate spline loadings
  spline.rr,                # reporting rate basis functions
  s2.hosp.start = .01,      # initial current hosp. variance 
  s2.icu.start = .01,       # initial current ICU variance 
  s2.vent.start = .01,      # initial current vent. variance 
  s2.beta.start = .01,      # initial beta prior variance
  s2.rr.start = .01,        # initial reporting rate prior variance
  const.params = NULL,      # vector of parameters to keep constant (optional)
  ode.params.start = NULL,      # ODESIM parameter initial values
  ode.params.prior.min = -Inf,  # lower bounds for ODESIM parameters (uniform prior)
  ode.params.prior.max = Inf,   # upper bounds for ODESIM parameters (uniform prior)
  non.odesim.params = NULL,     # vector of non-ODESIM parameters (such as "hosp.report.rate")
  lik.params.start = NULL,      # initial values for NB dispersion parameters
  fixed.nb.disp = FALSE,        # fix the NB dispersion parameters
  ### Adaptive proposal set-up ###
  adapt.iter = 100,           # interval to do log-adaptive tuning on M-H proposals
  indep.mh = FALSE,           # propose beta separate from other params. (Boolean)
  prop.type = "tnorm",        # M-H proposal type (should be "tnorm")
  adapt.type = "ShabyWells",  # type of adaptive tuning
  t.adapt.start = 0,          # adaptive iteration (should be 0)
  c0 = 1,                     # adaptive constant (see Shaby and Wells for details)
  c1 = 0.8,                   # adaptive constant (see Shaby and Wells for details)
  var.tune = NULL,            # list, order= beta,ode.params,rr,s2,lik
  Sigma.tune = NULL           # list, order= beta,ode.params,rr,s2,lik
){

  ########################################################################
  ### 1. Preliminaries 
  ########################################################################
  
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
    
  # structures to save parameter samples:
  beta.save <- matrix(NA_real_, nrow = n.mcmc/thin, ncol = num.beta)
  rr.params.save <- matrix(NA_real_, nrow = n.mcmc/thin, ncol = num.rr)
  ode.params.save <- matrix(NA_real_, nrow = n.mcmc/thin, ncol = num.ode.params)
  ode.params.names <- names(ode.params.start)
  colnames(ode.params.save) <- ode.params.names
  
  if (num.lik.params > 0) {
    lik.params.save <- matrix(NA_real_, nrow = n.mcmc / thin, ncol = num.lik.params)
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
  if (is.matrix(spline.rr)) {
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
  s2.beta <- s2.beta.start
  s2.rr <- s2.rr.start
  ode.params <- ode.params.start
  lik.params <- lik.params.start
  lik.params.star <- lik.params 

  # extra parameters (like hospitalization reporting rate )
  extra.params <- NULL
  extra.const.params <- NULL
  extra.params.fitted.idx <- integer()
  extra.params.const.idx <- integer()
    
  if (length(non.odesim.params) > 0) {
    for (k in 1:length(non.odesim.params)) {
      extra.params.fitted.idx <- c(
        extra.params.fitted.idx,
        which(names(ode.params) == non.odesim.params[k])
      )
    }

    if (length(extra.params.fitted.idx) > 0) {
      extra.params <- ode.params[extra.params.fitted.idx]
    }

    for (k in 1:length(non.odesim.params)) {
      extra.params.const.idx <- c(
        extra.params.const.idx,
        which(names(const.params) == non.odesim.params[k])
      )
    }

    if (length(extra.params.const.idx) > 0) {
      extra.const.params <- const.params[extra.params.const.idx]
    }
  }
  
  # parameter for estimation of prop. of symptomatic fractions
  if (sf.choice) {
    # if true, create structure for sampling symp-frac settings (6 options)

    # list of model possibilities
    symp.vals <- c(
      "", "-symp-frac-davies", "-symp-frac-equal 0.3",
      "-symp-frac-equal 0.4", "-symp-frac-equal 0.5",
      "-symp-frac-equal 0.6", "-symp-frac-equal 0.7"
    )

    # number of symp-frac options
    n.sf <- length(symp.vals)

    # initialize to random starting symp-frac option
    K <- sample.int(n.sf, size = 1)

    # structure to store models
    K.vals <- rep(NA_integer_, n.mcmc / thin)
    symp.cur <- symp.vals[K]
  } else {
    symp.cur <- NULL
    K.vals <- NULL
  }
    
  # create a matrix P to allow for a delay from symptoms to presentation
  P <- matrix(0, nrow = end.day + pres.delay + 1, ncol = end.day + 1)
  for (j in 1:ncol(P)) {
    P[j:(j + pres.delay), j] <- c(rep(0, pres.delay), 1)
  }
  
  ########################################################################
  ### 2. Define priors 
  ########################################################################

  ### beta prior: random walk
  ###   (penalized regression spline w/ 1st-order diffs)
  
  D <- diff(diag(num.beta), differences = 1)
  S <- crossprod(D)

  # beta ~ N(0, s2.beta * S^-1)
  beta.prior.loglik <- function(beta, s2.beta, precision = S) {
    if (min(beta) < 0) {
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
  rr.prior.loglik <- function(rr, s2.rr, precision = S.rr) {
    if (min(rr) < 0 | max(rr) > 1) {
      ll <- -Inf
    } else {
      ll <- -1 / 2 / s2.rr * (t(rr) %*% (precision %*% rr))
    }
    return(ll)
  }
  
  ### dispersion parameter: exponential prior
  ###   disp ~ Exp(lambda = 100)
  lik.params.prior.loglik <- function(lik.params,
                                      lambda = rep(100, length(lik.params))) {
    if (length(lik.params) < 1) {
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
  ode.params.prior.loglik <- function(ode.params, ode.params.min, ode.params.max) {
    ll <- 0
    if (sum(c(ode.params < ode.params.min, ode.params > ode.params.max)) > 0) {
      ll <- -Inf
    }
    return(ll)
  }
  
  ### Uniform priors for s2 parameters
  s2.params.prior.loglik <- function(s2.params) {
    ll <- 0
    if (min(s2.params) < 0) {
      ll <- -Inf
    }
    return(ll)
  }
  
  ########################################################################
  ### 3. Initial likelihood evaluation 
  ########################################################################

  # simulate trajectory using current beta values:
  traj <- traj.from.params(
    beta = beta.daily,
    params = ode.params,
    tf = end.day,
    introday = introday,
    const.params = const.params,
    non.odesim.params = non.odesim.params,
    odepath = odepath,
    loc = loc,
    symp = symp.cur
  )
  
  # evaluate logliklihood under initial conditions
  llvals <- loglik.odesim(
    traj = traj,
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
    lik.curr = lik.curr,
    lik.old = lik.old,
    active.surv = active.surv,
    p.asympt = p.asympt,
    total.size.constraint = total.size.constraint,
    s2.hosp = s2.hosp,
    s2.icu = s2.icu,
    s2.vent = s2.vent,
    extra.params = extra.params,
    extra.const.params = extra.const.params
  )
  
  ll.current <- llvals$ll
  ll.new <- llvals$ll.new
  ll.hosp.new <- llvals$ll.hosp.new
  
  ########################################################################
  ### 4. Adaptive proposal set-up 
  ########################################################################
  accept <- rep(0, length(var.tune))
  accept.tot <- accept
  never.adapt <- TRUE

  ########################################################################
  ### 5. MCMC iterations 
  ########################################################################
  for (iter in 1:n.mcmc) {

    # print-out iteration number every 100 iterations
    if (iter %% 100 == 0 & print.iter) {
      cat(iter, " ")
    }
    
    #########################################################################
    ### A. Propose beta
    #########################################################################

    # setting up error catching - reject beta if odesim fails
    beta.good <- 0
    
    while (beta.good < 1) {

      # propose beta.star
      beta.star <- exp(rnorm(
        length(beta), log(beta),
        sqrt(var.tune[1] * diag(Sigma.tune[[1]]))
      ))
      beta.daily.star <- Z %*% beta.star

      # trajectory given beta.star
      traj.star <- try(traj.from.params(
        beta = beta.daily.star,
        params = ode.params,
        tf = end.day,
        introday = introday,
        const.params = const.params,
        non.odesim.params = non.odesim.params,
        odepath = odepath,
        loc = loc,
        symp = symp.cur
      ), silent = TRUE)

      # check if error was thrown. if not, leave while loop
      if (class(traj.star) != "try-error") {
        beta.good <- 1
      }
    }

    # evaluate loglikelihood for beta.star
    llvals.star <- loglik.odesim(
      traj = traj.star,
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
      lik.curr = lik.curr,
      lik.old = lik.old,
      active.surv = active.surv,
      p.asympt = p.asympt,
      total.size.constraint = total.size.constraint,
      df = df,
      odesim.ver = odesim.ver,
      P = P,
      loc = loc,
      extra.params = extra.params,
      extra.const.params = extra.const.params
    ) 
    ll.star <- llvals.star$ll
        
    # accept/reject beta.star    
    mh1 <- ll.star + beta.prior.loglik(beta.star, s2.beta) + sum(log(beta.star))
    mh2 <- ll.current + beta.prior.loglik(beta, s2.beta) + sum(log(beta))
    
    if (is.na(mh1)) {
      mh1 <- -Inf
    }
    
    # if Unif(0,1) < mh1/mh2, accept new beta and disp parameters
    if (exp(mh1 - mh2) > runif(1)) {

      # accept beta.star
      beta <- beta.star
      beta.daily <- beta.daily.star
      traj <- traj.star
      llvals = llvals.star
      ll.current <- ll.star
      accept[1] <- accept[1] + 1
    }

    #########################################################################
    ### B. Propose ode.params
    #########################################################################

    # set up error catching - reject ode.params if odesim fails
    beta.good <- 0
    while (beta.good < 1) {
      if (prop.type == "norm") {
        ode.params.star <- t(rmvnorm(1, ode.params, var.tune[2] * Sigma.tune[[2]]))
      }
      if (prop.type == "tnorm") {
        ode.params.star <- (rtnorm(
          length(ode.params), ode.params,
          sqrt(var.tune[2] * diag(Sigma.tune[[2]])), ode.params.prior.min, ode.params.prior.max
        ))
      }

      names(ode.params.star) <- ode.params.names
      extra.params.star <- ode.params.star[extra.params.fitted.idx]

      # trajectory given ode.params.star
      traj.star <- try(traj.from.params(
        beta = beta.daily,
        params = ode.params.star,
        tf = end.day,
        introday = introday,
        const.params = const.params,
        non.odesim.params = non.odesim.params,
        odepath = odepath,
        loc = loc,
        symp = symp.cur
      ), silent = TRUE)

      if (class(traj.star) != "try-error" & min(ode.params.star > 0)) {
        beta.good <- 1
      }
    }
    
    # evalulate loglikelihood for ode.params.star
    llvals.star <- loglik.odesim(
      traj = traj.star,
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
      lik.curr = lik.curr,
      lik.old = lik.old,
      active.surv = active.surv,
      p.asympt = p.asympt,
      total.size.constraint = total.size.constraint,
      df = df,
      odesim.ver = odesim.ver,
      P = P,
      loc = loc,
      extra.params = extra.params.star,
      extra.const.params = extra.const.params
    )
    ll.star <- llvals.star$ll
        
    # accept/reject ode.params.star
    mh1 <- ll.star + ode.params.prior.loglik(
      ode.params.star,
      ode.params.prior.min,
      ode.params.prior.max
    )
    
    mh2 <- ll.current + ode.params.prior.loglik(
      ode.params,
      ode.params.prior.min,
      ode.params.prior.max
    )
    
    if (prop.type == "tnorm") {
      mh1 <- mh1 + sum(dtnorm(ode.params,
        ode.params.star,
        sqrt(var.tune[2] * diag(Sigma.tune[[2]])),
        ode.params.prior.min,
        ode.params.prior.max,
        log = TRUE
      ))

      mh2 <- mh2 + sum(dtnorm(ode.params.star,
        ode.params,
        sqrt(var.tune[2] * diag(Sigma.tune[[2]])),
        ode.params.prior.min,
        ode.params.prior.max,
        log = TRUE
      ))
    }
    
    if (is.na(mh1)) {
      mh1 <- -Inf
    }
    
    # if Unif(0,1) < mh1/mh2, accept new odesim params
    if (exp(mh1 - mh2) > runif(1)) {
      ode.params <- ode.params.star
      extra.params <- extra.params.star
      traj <- traj.star
      llvals <- llvals.star
      ll.current <- ll.star
      accept[2] <- accept[2] + 1
    }
    
    #########################################################################
    ### C. Propose symp-frac option (if sf.choice == TRUE)
    #########################################################################

    if (sf.choice) {

      # propose model from uniform prior
      K.star <- sample(n.sf, size = 1)
      symp.star <- symp.vals[K.star]

      # calculate traj and likelihoods for each possible model:
      traj.vals <- list()
      good.k <- TRUE

      # trajectory given beta.star
      traj.star <- try(traj.from.params(
        beta = beta.daily,
        params = ode.params,
        tf = end.day,
        introday = introday,
        const.params = const.params,
        non.odesim.params = non.odesim.params,
        odepath = odepath,
        loc = loc,
        symp = symp.star
      ),
      silent = TRUE
      )

      if (class(traj.star) == "try-error") {
        good.k <- FALSE
      }

      if (good.k) {
        llvals.star <- loglik.odesim(
          traj = traj.star,
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
          lik.curr = lik.curr,
          lik.old = lik.old,
          active.surv = active.surv,
          p.asympt = p.asympt,
          total.size.constraint = total.size.constraint,
          df = df,
          odesim.ver = odesim.ver,
          P = P,
          loc = loc,
          extra.params = extra.params,
          extra.const.params = extra.const.params
        )
        ll.star <- llvals.star$ll

        # accept/reject new symp frac option
        mh1 <- ll.star
        mh2 <- ll.current

        if (is.na(mh1)) {
          mh1 <- -Inf
        }

        # if Unif(0,1) < mh1/mh2, accept symp.frac.star
        if (exp(mh1 - mh2) > runif(1)) {
          K <- K.star
          symp.cur <- symp.star
          llvals <- llvals.star
          ll.current <- ll.star
          traj <- traj.star
        }
      }
    }
    
    #########################################################################
    ### D. Propose rr.params
    #########################################################################

    # sample new rr params
    rr.params.star <- rtnorm(
      length(rr.params), rr.params,
      sqrt(diag(var.tune[3] * Sigma.tune[[3]])), 0, 1
    )
    rr.daily.star <- Z.rr %*% rr.params.star
        
    # make sure rr < 1
    if (max(rr.daily.star, na.rm = T) < 1) {

      # evaluate logliklihood with rr.params.star
      llvals.star <- loglik.odesim(
        traj = traj,
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
        lik.curr = lik.curr,
        lik.old = lik.old,
        active.surv = active.surv,
        p.asympt = p.asympt,
        total.size.constraint = total.size.constraint,
        df = df,
        odesim.ver = odesim.ver,
        P = P,
        loc = loc,
        extra.params = extra.params,
        extra.const.params = extra.const.params
      )
      ll.star <- llvals.star$ll

      # accept/reject rr.params.star
      mh1 <- ll.star + rr.prior.loglik(rr.params.star,
        s2.rr,
        precision = S.rr
      )

      mh2 <- ll.current + rr.prior.loglik(rr.params,
        s2.rr,
        precision = S.rr
      )

      if (prop.type == "tnorm") {
        mh1 <- mh1 + sum(dtnorm(rr.params,
          rr.params.star,
          sqrt(var.tune[3] * diag(Sigma.tune[[3]])),
          0, 1,
          log = TRUE
        ))

        mh2 <- mh2 + sum(dtnorm(rr.params.star,
          rr.params,
          sqrt(var.tune[3] * diag(Sigma.tune[[3]])),
          0, 1,
          log = TRUE
        ))
      }

      if (is.na(mh1)) {
        mh1 <- -Inf
      }

      # if Unif(0,1) < mh1/mh2, accept rr.params
      if (exp(mh1 - mh2) > runif(1)) {
        rr.params <- rr.params.star
        rr.daily <- rr.daily.star
        llvals <- llvals.star
        ll.current <- ll.star
        accept[3] <- accept[3] + 1
      }
    }
    
    #########################################################################
    ### E. Propose s2.params
    #########################################################################      

    if(!lik.hosp.curr & !lik.vent.curr & !lik.icu.curr) {
      # no need to propose new s2 parameters
    } else {
        
      # s2.params.star <- t(rtnorm(length(s2.params), s2.params, sqrt(diag(var.tune[4] * Sigma.tune[[4]])),0,Inf))
      s2.params.star <- exp(t(rnorm(
        length(s2.params),
        log(s2.params),
        sqrt(diag(var.tune[4] * Sigma.tune[[4]]))
      )))
      
      # evaluate logliklihood for s2.params.star
      llvals.star <- loglik.odesim(
        traj = traj,
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
        lik.curr = lik.curr,
        lik.old = lik.old,
        active.surv = active.surv,
        p.asympt = p.asympt,
        total.size.constraint = total.size.constraint,
        df = df,
        odesim.ver = odesim.ver,
        P = P,
        loc = loc,
        extra.params = extra.params,
        extra.const.params = extra.const.params
      ) 
      ll.star <- llvals.star$ll

      # accept/reject s2.params.star
  
      # determine which s2 values to include or not
      s2.true <- rep(FALSE, 3)
      mh1 <- 0; mh2 <- 0
      
      if (length(llvals.star$ll.hosp) > 0) {
        mh1 <- mh1 + llvals.star$ll.hosp
        mh2 <- mh2 + llvals$ll.hosp
        s2.true[1] <- TRUE
      }
      
      if (length(llvals.star$ll.icu) > 0) {
        mh1 <- mh1 + llvals.star$ll.icu
        mh2 <- mh2 + llvals$ll.icu
        s2.true[2] <- TRUE
      }
      
      if (length(llvals.star$ll.vent) > 0) {
        mh1 <- mh1 + llvals.star$ll.vent
        mh2 <- mh2 + llvals$ll.vent
        s2.true[3] <- TRUE
      }
        
      mh1 <- mh1 + s2.params.prior.loglik(s2.params.star[s2.true])
      mh2 <- mh2 + s2.params.prior.loglik(s2.params[s2.true])
          
      if (prop.type == "tnorm") {
        mh1 <- mh1 + sum(log(s2.params.star[s2.true]))
        mh2 = mh2 + sum(log(s2.params[s2.true]))
      }
          
      if (is.na(mh1)) {
        mh1 <- -Inf
      }
      
      # if Unif(0,1) < mh1/mh2, accept new s2.params
      if (exp(mh1 - mh2) > runif(1)) {
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
    ### E. Propose lik.params (NB dispersion params)
    #########################################################################
    if (!fixed.nb.disp) {
      lik.params.star <- exp(t(rnorm(
        length(lik.params),
        log(lik.params),
        sqrt(diag(var.tune[5] * Sigma.tune[[5]]))
      )))

      # evaluate logliklihood with lik.params.star
      llvals.star <- loglik.odesim(
        traj = traj,
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
        lik.curr = lik.curr,
        lik.old = lik.old,
        active.surv = active.surv,
        p.asympt = p.asympt,
        total.size.constraint = total.size.constraint,
        df = df,
        odesim.ver = odesim.ver,
        P = P,
        loc = loc,
        extra.params = extra.params,
        extra.const.params = extra.const.params
      )
      ll.star <- llvals.star$ll

      # accept/reject lik.params.star
      mh1 <- ll.star +
        lik.params.prior.loglik(lik.params.star) +
        sum(log(lik.params.star))

      mh2 <- ll.current +
        lik.params.prior.loglik(lik.params) +
        sum(log(lik.params))

      if (is.na(mh1)) {
        mh1 <- -Inf
      }

      # if Unif(0,1) < mh1/mh2, accept new lik.params
      if (exp(mh1 - mh2) > runif(1)) {
        lik.params <- lik.params.star
        llvals <- llvals.star
        ll.current <- ll.star
        accept[5] <- accept[5] + 1
      }
    }
    
    #########################################################################
    ### F. Propose new s2.beta and s2.rr parameters
    #########################################################################

    # conjugate updates from inverse gamma marg. post.
    s2.beta <- 1 / rgamma(1,
      shape = s + num.beta / 2,
      rate = r + .5 * (t(beta) %*% (S %*% beta))
    )
    
    s2.rr <- 1 / rgamma(1,
      shape = s + num.rr / 2,
      rate = r + .5 * (t(rr.params) %*% (S.rr %*% rr.params))
    )
      
    #########################################################################
    ### G. Save parameters
    #########################################################################

    if (iter %% thin == 0) {
      beta.save[iter / thin, ] <- beta
      ode.params.save[iter / thin, ] <- ode.params
      rr.params.save[iter / thin, ] <- rr.params
      s2.params.save[iter / thin, ] <- s2.params
      lik.params.save[iter / thin, ] <- lik.params
      s2.beta.save[iter / thin] <- s2.beta
      s2.rr.save[iter / thin] <- s2.rr
      loglik.save[iter / thin] <- ll.current

      if (sf.choice) {
        K.vals[iter / thin] <- K
      }
    }
    

    
    #########################################################################
    ### H. Update adaptive tuning variances
    #########################################################################
    if (iter %% adapt.iter == 0) {
      never.adapt <- FALSE

      # Using Shaby and Wells adaptive scheme for RWM...
      if (adapt.type == "ShabyWells") {
        # default constants
        r.opt <- 0.234
        r.hat <- accept / adapt.iter
        t.adapt <- iter / adapt.iter + t.adapt.start
        gamma.1 <- 1 / (t.adapt)^c1
        gamma.2 <- c0 * gamma.1
        var.tune <- exp(log(var.tune) + gamma.2 * (r.hat - r.opt))

        accept.tot <- accept.tot + accept
        accept <- rep(0, length(accept))
        cat("\n", "current accept rate =", r.hat)
        cat("\n", "new proposal var =", var.tune, "\n")
      }
    }
    
    #########################################################################
    ### I. Plot output to assess convergence
    #########################################################################

    if (iter %% plot.rate == 0 & plot.save == TRUE) {
      pdf(plot.name, width = 10, height = 7)
      matplot(beta.save[1:iter / thin, ],
        type = "l",
        main = "beta"
      )
      matplot(ode.params.save[1:iter / thin, ],
        type = "l",
        main = "ode.params"
      )
      matplot(rr.params.save[1:iter / thin, ],
        type = "l",
        main = "report.rate"
      )
      matplot(s2.params.save[1:iter / thin, ],
        type = "l",
        main = "s2.params"
      )
      matplot(log(lik.params.save[1:iter / thin, ]),
        type = "l",
        main = "log(NB dispersion params)"
      )
      dp <- data.process(df, loc = loc)
      tp <- traj.process(traj, loc = loc, odesim.version = odesim.ver)
      plots.odesim(dp, tp, rr.daily[-length(rr.daily)])
      dev.off()
    }
  } 

  # new line
  if (print.iter) {
    cat("\n") 
  }

  ########################################################################
  ### 6. Adaptive structures for output 
  ########################################################################
  accept.rate.final <- accept / n.mcmc
  Sigma.hat <- list(
    cov(beta.save),
    cov(ode.params.save),
    cov(rr.params.save),
    cov(s2.params.save),
    cov(lik.params.save)
  )

  ########################################################################
  ### 7. Save MCMC results
  ########################################################################

  # output MCMC samples in a list
  list(
    beta = beta.save,
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
    lik.curr = lik.curr,
    lik.old = lik.old,
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
    P = P,
    extra.params = extra.params,
    extra.const.params = extra.const.params,
    active.surv = active.surv,
    today = Sys.Date()
  )
}
