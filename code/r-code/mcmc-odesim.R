##############################################
### Most Recent MCMC Function              ###
### Last edited: 18 August 2020            ###
##############################################

mcmc.odesim <- function(df, ## data.frame object
                    lik.tot=TRUE, ## If TRUE, evaluate the likelihood for the total new cases
                    lik.age=FALSE, ## if TRUE, evaluate the likelihood for the age stratified new cases
                    lik.hosp.new=TRUE, ## if TRUE, evaluate the likelihood for the new hosp. cases
                    lik.hosp.curr=FALSE,
                    lik.icu.curr=FALSE,
                    lik.vent.curr=FALSE,
                    lik.tot.deaths=FALSE,
                    lik.home.deaths=FALSE,
                    lik.hosp.deaths=FALSE,
                    lik.age.deaths=FALSE,
                    lik.hosp.discharges=FALSE,
                    active.surv=FALSE,
                    p.asympt=0.4,
                    total.size.constraint=FALSE,
                    odesim.ver="v4",
                    beta.start, ## starting values for spline loadings
                    spline.beta, ## fda "spline" object spanning the time window in "days"
                    s2.hosp.start=.01,
                    s2.icu.start=.01,
                    s2.vent.start=.01,
                    s2.beta.start=.01,
                    s2.rr.start=.01, 
                    const.params = NULL, ## vector of parameters to keep constant, if desired
                    report.rate.params.start, ## rate at which infecteds report
                    spline.rr,
                    ode.params.start=NULL, ## additional inputs to odesim.  See below for details
                    ode.params.prior.min=-Inf, ## bounds for uniform priors
                    ode.params.prior.max=Inf,  ## bounds for uniform priors
                    non.odesim.params=NULL, ## string vector of names of params in ode.params that are NOT used in odesim.  First case is "hosp.report.rate"
                    lik.params.start=NULL, ## starting values for dispersion parameters in NB likelihood
                    start.day=61,end.day=116, ## start and end day to the odesim
                    introday=NULL, ## day of first infected
                    loc="RI",
                    odepath, ## path to "odesim"
                    n.mcmc, ## number of MCMC iterations to run
                    adapt.iter = 100, ## interval to do log-adaptive tuning on MH
                    indep.mh=FALSE, ## if true, propose beta separate from other params.
                    t.adapt.start=0,
                    prop.type="tnorm",
                    adapt.type="ShabyWells",
                    c0=1, ## Shaby Wells adaptive params
                    c1=0.8, ## Shaby Wells adaptive params
                    var.tune=NULL, ## list, order= beta,ode.params,rr,s2,lik
                    Sigma.tune=NULL, ## list, order= beta,ode.params,rr,s2,lik
                    p.vecs, ## weekly vector of delay probabilities
                    thin=1, ## rate to thin saved posterior samples
                    plot.save=TRUE,
                    plot.rate=sqrt(2),
                    plot.name="trace.plots.pdf",
                    print.iter=FALSE,
                    sf.choice = FALSE, ## if true, includes parameter for symp-frac,
                    fixed.nb.disp = FALSE,
                    mult = 1
                    ){


  ########################
  ### 1. Preliminaries ###
  ########################
  
  # adaptive structure
  t.adapt <- t.adapt.start

  # process data frame for likelihood evaluation
  data.processed <- data.process(df, loc = loc, mult = mult)
  
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

  ## ### function to fill in missing age-structured data
  ## impute <- function(total.vals, age.vals,rates){
  ##     diff.num=total.vals-apply(age.vals,1,sum)
  ##     idx=which(diff.num>0)
  ##     ## for(i in idx){
  ##     ##     age.vals[i,]=age.vals[i,]+rmultinom(1,diff.num[i],rates[i,])
  ##     ## }
  ##     age.vals[idx,]=age.vals[idx,]+rmultinomial(length(idx),diff.num[idx],rates[idx,])
  ##     as.matrix(age.vals)
  ## }

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

  # traj.save <- array(NA,dim=c(n.mcmc/thin,dim(as.matrix(traj))))
  
 # browser()
  
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
                          total.size.constraint = total.size.constraint,
                          s2.hosp = s2.hosp,
                          s2.icu = s2.icu,
                          s2.vent = s2.vent,
                          extra.params = extra.params, ### if cumul. hosp. reporting rate is fitted 
                          extra.const.params = extra.const.params,
                          mult = mult) ### if cumul. hosp. reporting rate is constant
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
                                 total.size.constraint = total.size.constraint,
                                 df = df,
                                 odesim.ver = odesim.ver,
                                 P = P,
                                 loc = loc,
                                 extra.params = extra.params,  
                                 extra.const.params = extra.const.params,
                                 mult = mult) 
    ll.star <- llvals.star$ll

    ## ### resample imputations 
    ## llvals <- loglik.odesim(traj = traj,
    ##                         dp = data.processed,
    ##                         report.rate = rr.daily,
    ##                         nb.disp.params = lik.params,
    ##                         s2.hosp = s2.params[1],
    ##                         s2.icu = s2.params[2],
    ##                         s2.vent = s2.params[3],
    ##                         lik.tot = lik.tot,
    ##                         lik.age = lik.age,
    ##                         lik.hosp.new = lik.hosp.new,
    ##                         lik.hosp.curr = lik.hosp.curr,
    ##                         lik.icu.curr = lik.icu.curr, 
    ##                         lik.vent.curr = lik.vent.curr,
    ##                         lik.tot.deaths = lik.tot.deaths,
    ##                         lik.home.deaths = lik.home.deaths,
    ##                         lik.hosp.deaths = lik.hosp.deaths,
    ##                         lik.age.deaths = lik.age.deaths,
    ##                         lik.hosp.discharges = lik.hosp.discharges,
    ##                         active.surv = active.surv,
    ##                         p.asympt = p.asympt,
    ##                         total.size.constraint = total.size.constraint,
    ##                         df = df,
    ##                         odesim.ver = odesim.ver,
    ##                         P = P,
    ##                         loc = loc,
    ##                         extra.params = extra.params) 
    ## ll.current <- llvals$ll
        
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
                                 total.size.constraint = total.size.constraint,
                                 df = df,
                                 odesim.ver = odesim.ver,
                                 P = P,
                                 loc = loc,
                                 extra.params = extra.params.star, 
                                 extra.const.params = extra.const.params,
                                 mult = mult) 
    ll.star <- llvals.star$ll

    ## ### resample imputations and calc. likelihood for current odesim params
    ## llvals <- loglik.odesim(traj = traj,
    ##                         dp = data.processed,
    ##                         report.rate = rr.daily,
    ##                         nb.disp.params = lik.params,
    ##                         s2.hosp = s2.params[1],
    ##                         s2.icu = s2.params[2],
    ##                         s2.vent = s2.params[3],
    ##                         lik.tot = lik.tot,
    ##                         lik.age = lik.age,
    ##                         lik.hosp.new = lik.hosp.new,
    ##                         lik.hosp.curr = lik.hosp.curr,
    ##                         lik.icu.curr = lik.icu.curr, 
    ##                         lik.vent.curr = lik.vent.curr,
    ##                         lik.tot.deaths = lik.tot.deaths,
    ##                         lik.home.deaths = lik.home.deaths,
    ##                         lik.hosp.deaths = lik.hosp.deaths,
    ##                         lik.age.deaths = lik.age.deaths,
    ##                         lik.hosp.discharges = lik.hosp.discharges,
    ##                         active.surv = active.surv,
    ##                         p.asympt = p.asympt,
    ##                         total.size.constraint = total.size.constraint,
    ##                         df = df,
    ##                         odesim.ver = odesim.ver,
    ##                         P = P,
    ##                         loc = loc,
    ##                         extra.params = extra.params, 
    ##                         extra.const.params = extra.const.params) 
    ## ll.current <- llvals$ll

        
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
                                 total.size.constraint = total.size.constraint,
                                 df = df,
                                 odesim.ver = odesim.ver,
                                 P = P,
                                 loc = loc,
                                 extra.params = extra.params, 
                                 extra.const.params = extra.const.params,
                                 mult = mult) 
        ll.star <- llvals.star$ll
        
        
        ## ### resample imputations, reevaluate likelihood for current rr.params
        ## llvals <- loglik.odesim(traj = traj,
        ##                         dp = data.processed,
        ##                         report.rate = rr.daily,
        ##                         nb.disp.params = lik.params,
        ##                         s2.hosp = s2.params[1],
        ##                         s2.icu = s2.params[2],
        ##                         s2.vent = s2.params[3],
        ##                         lik.tot = lik.tot,
        ##                         lik.age = lik.age,
        ##                         lik.hosp.new = lik.hosp.new,
        ##                         lik.hosp.curr = lik.hosp.curr,
        ##                         lik.icu.curr = lik.icu.curr, 
        ##                         lik.vent.curr = lik.vent.curr,
        ##                         lik.tot.deaths = lik.tot.deaths,
        ##                         lik.home.deaths = lik.home.deaths,
        ##                         lik.hosp.deaths = lik.hosp.deaths,
        ##                         lik.age.deaths = lik.age.deaths,
        ##                         lik.hosp.discharges = lik.hosp.discharges,
        ##                         active.surv = active.surv,
        ##                         p.asympt = p.asympt,
        ##                         total.size.constraint = total.size.constraint,
        ##                         df = df,
        ##                         odesim.ver = odesim.ver,
        ##                         P = P,
        ##                         loc = loc,
        ##                         extra.params = extra.params, 
        ##                         extra.const.params = extra.const.params) 
        ## ll.current <- llvals$ll
      
      
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
                                   total.size.constraint = total.size.constraint,
                                   df = df,
                                   odesim.ver = odesim.ver,
                                   P = P,
                                   loc = loc,
                                   extra.params = extra.params, 
                                   extra.const.params = extra.const.params,
                                   mult = mult) 
      ll.star <- llvals.star$ll

      ## ### resample imputations, reevaluate likelihood for current rr.params
      ## llvals <- loglik.odesim(traj = traj,
      ##                         dp = data.processed,
      ##                         report.rate = rr.daily,
      ##                         nb.disp.params = lik.params,
      ##                         s2.hosp = s2.params[1],
      ##                         s2.icu = s2.params[2],
      ##                         s2.vent = s2.params[3],
      ##                         lik.tot = lik.tot,
      ##                         lik.age = lik.age,
      ##                         lik.hosp.new = lik.hosp.new,
      ##                         lik.hosp.curr = lik.hosp.curr,
      ##                         lik.icu.curr = lik.icu.curr, 
      ##                         lik.vent.curr = lik.vent.curr,
      ##                         lik.tot.deaths = lik.tot.deaths,
      ##                         lik.home.deaths = lik.home.deaths,
      ##                         lik.hosp.deaths = lik.hosp.deaths,
      ##                         lik.age.deaths = lik.age.deaths,
      ##                         lik.hosp.discharges = lik.hosp.discharges,
      ##                         active.surv = active.surv,
      ##                         p.asympt = p.asympt,
      ##                         total.size.constraint = total.size.constraint,
      ##                         df = df,
      ##                         odesim.ver = odesim.ver,
      ##                         P = P,
      ##                         loc = loc,
      ##                         extra.params = extra.params, 
      ##                         extra.const.params = extra.const.params) 
      ## ll.current <- llvals$ll

          
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
                                   total.size.constraint = total.size.constraint,
                                   df = df,
                                   odesim.ver = odesim.ver,
                                   P = P,
                                   loc = loc,
                                   extra.params = extra.params, 
                                   extra.const.params = extra.const.params,
                                   mult = mult) 
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

      # mh1 <- llvals.star$ll.hosp+llvals.star$ll.vent+llvals.star$ll.icu +  s2.params.prior.loglik(s2.params.star)
      # mh2 <- llvals$ll.hosp+llvals$ll.vent+llvals$ll.icu + s2.params.prior.loglik(s2.params)
          
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
                                       total.size.constraint = total.size.constraint,
                                       df = df,
                                       odesim.ver = odesim.ver,
                                       P = P,
                                       loc = loc,
                                       extra.params = extra.params, 
                                       extra.const.params=extra.const.params,
                                       mult = mult) 
          ll.star <- llvals.star$ll
          
          ## accept/reject lik.params.star
          
          mh1 <- ll.star +  
              lik.params.prior.loglik(lik.params.star) + 
              sum(log(lik.params.star))
          
          mh2 <- ll.current + 
              lik.params.prior.loglik(lik.params) + 
              sum(log(lik.params))
          
                                        # if(prop.type=="tnorm"){
                                        #   mh1=mh1+sum(dtnorm(lik.params,lik.params.star,sqrt(var.tune[5]*diag(Sigma.tune[[5]])),0,Inf,log=TRUE))
                                        #   mh2=mh2+sum(dtnorm(lik.params.star,lik.params,sqrt(var.tune[5]*diag(Sigma.tune[[5]])),0,Inf,log=TRUE))
                                        # }
          
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
                
        # Sigma.hat=Sigma.tune
        # Sigma.hat[[1]]=cov(beta.save[(iter/thin - adapt.iter/thin + 1):iter/thin,])
            
        # diag(Sigma.hat[[2]])=apply(ode.params.save[(iter/thin - adapt.iter/thin + 1):iter/thin,],2,var)
        # diag(Sigma.hat[[3]])=apply(rr.params.save[(iter/thin - adapt.iter/thin + 1):iter/thin,],2,var)
        # diag(Sigma.hat[[4]])=apply(s2.params.save[(iter/thin - adapt.iter/thin + 1):iter/thin,],2,var)
        # diag(Sigma.hat[[5]])=apply(lik.params.save[(iter/thin - adapt.iter/thin + 1):iter/thin,],2,var)diag(Sigma.hat[[1]])=apply(beta.save[(iter/thin - adapt.iter/thin + 1):iter/thin,],2,var)
        # diag(Sigma.hat[[2]])=apply(ode.params.save[(iter/thin - adapt.iter/thin + 1):iter/thin,],2,var)
        # diag(Sigma.hat[[3]])=apply(rr.params.save[(iter/thin - adapt.iter/thin + 1):iter/thin,],2,var)
        # diag(Sigma.hat[[4]])=apply(s2.params.save[(iter/thin - adapt.iter/thin + 1):iter/thin,],2,var)
        # diag(Sigma.hat[[5]])=apply(lik.params.save[(iter/thin - adapt.iter/thin + 1):iter/thin,],2,var)
        # # diag(Sigma.hat[[2]])=var(ode.params.save[(iter/thin - adapt.iter/thin + 1):iter/thin,])
        # # diag(Sigma.hat[[3]])=var(rr.params.save[(iter - adapt.iter + 1):iter,])
        # # diag(Sigma.hat[[4]])=var(s2.params.save[(iter - adapt.iter + 1):iter,])
        # # diag(Sigma.hat[[5]])=var(lik.params.save[(iter - adapt.iter + 1):iter,])
        # for(iii in 1:length(Sigma.hat)){
        #  Sigma.tune[[iii]] <- Sigma.tune[[iii]] + gamma.1 * (Sigma.hat[[iii]] - Sigma.tune[[iii]])
        # }
                
        accept.tot <- accept.tot + accept
        accept <- rep(0, length(accept))
        cat("\n","current accept rate =", r.hat)
        cat("\n","new proposal var =", var.tune, "\n")
      }

      # ### update based on empirical covar of most recent block
      #          
      # ### this block gives beta an independent update from ode.params
      # if(indep.mh){
      #   Sigma.hat.tmp=Sigma.hat
      #   Sigma.hat[1:nrow(Sigma.hat),1:nrow(Sigma.hat)]=0
      #   Sigma.hat[1:num.beta,1:num.beta]=Sigma.hat.tmp[1:num.beta,1:num.beta]
      #   Sigma.hat[num.beta+(1:num.ode.params),num.beta+(1:num.ode.params)]=Sigma.hat.tmp[num.beta+(1:num.ode.params),num.beta+(1:num.ode.params)]
      #   if(num.lik.params>0){
      #       Sigma.hat[num.beta+num.ode.params+(1:num.lik.params),num.beta+num.ode.params+(1:num.lik.params)]=Sigma.hat.tmp[num.beta+num.ode.params+(1:num.lik.params),num.beta+num.ode.params+(1:num.lik.params)]
      #   }
      # }
      #     
      # ### this block makes proposal independent
      # if(indep.mh){
      #   Sigma.hat.tmp=Sigma.hat
      #   Sigma.hat[1:nrow(Sigma.hat),1:nrow(Sigma.hat)]=0
      #   diag(Sigma.hat)=diag(Sigma.hat.tmp)
      # }
      #
      # r.hat <- accept / adapt.iter
      # gamma.1 <- 1 / (iter / adapt.iter)^c1
      # gamma.2 <- c0 * gamma.1
      # var.tune <- exp(log(var.tune) + gamma.2 * (r.hat - r.opt))
      # Sigma.tune <- Sigma.tune + gamma.1 * (Sigma.hat - Sigma.tune)
      # reset acceptance
      # accept.tot=accept.tot+accept
      # accept <- rep(0,length(accept))
      # 
      # if(adapt.type=="Covar"){
      #   Sigma.hat <- cov( cbind(beta.save,ode.params.save,s2.hosp.save,s2.icu.save,s2.vent.save,lik.params.save)[1:iter,])
      #   Sigma.tune=2.4^2/nrow(Sigma.hat)*Sigma.hat
      # }
      
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
      dp <- data.process(df, loc = loc, mult = mult)
      tp <- traj.process(traj, loc = loc, odesim.version = odesim.ver)
      plots.odesim(dp, tp, rr.daily[-length(rr.daily)])
      dev.off()
    }

    # # DIC calculations
    # if(iter/n.mcmc>.5){
    #   dbar <- dbar + (2 / n.mcmc) * (-2 * ll.current)
    #   if(inf.type=="age"){
    #     sympt.new.hat <- sympt.new.hat + (2 / n.mcmc) * sympt.new.imputed
    #   }
    # }

  } ## end MCMC

  if(print.iter){
    cat("\n") # new line
  }

  # ### DIC calculations (using 2nd half of chain)
  # ###   emh: DIC needs to be changed for age-structured data. Not comparable yet.
  #
  # beta.hat=apply(beta.save[-c(1:(n.mcmc/2)),],2,mean)
  # beta.daily.hat=Z%*%beta.hat
  # traj.hat=traj.from.params.v4(beta.daily.hat,tf=end.day,introday=introday,odepath=odepath)
  # delay.rates.hat <- (P[, traj.hat$times[delta.idx] ] %*% traj.hat$sympt.new[delta.idx,])+.Machine$double.eps[1]
  # ### log-likelihood given trajectory sim:
  # ll.hat <- sum(dpois(as.numeric(round(as.matrix(sympt.new.hat))), lambda = report.rate * as.numeric(delay.rates.hat[delta.idx,]), log = TRUE))
  #
  # dhat=-2*ll.hat
  # dic=2*dbar-dhat
  
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
       # dic = dic,
       # dbar = dbar,
       # dhat = dhat,
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
