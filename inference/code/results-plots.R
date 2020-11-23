#!/usr/bin/env Rscript

# results-plots.R
# authors: Nathan Wikle
# last edited: 23 Nov 2020 
#
# Function to generate output plots from fitted model, with 
#   multiple comparisons of fit to the data.

library(ggplot2)
library(stringr)
library(ggstance)

### A. Grab S samples from posterior
sample.grab <- function(S, directory, burnin = 0){
  # Input: 
  #   S: number of samples to grab from posterior output.
  #   directory: location of posterior output.
  #   burnin: number of burn-in samples to leave out.
  # Output:
  #   A sample of S log-likelihood evaluations saved with posterior output.
  
  # grab mcmc output
  mcmc.output <- mcmc.params(out.folder = directory, burnin = burnin)
  inds <- dim(mcmc.output$beta)
  
  # subsample the mcmc output
  subsamples <- cbind(sample(1:inds[1], size = S),
                      sample(1:inds[3], size = S, replace = T))
  
  # store results as a list
  res <- list()
  
  res$betas <- t(sapply(1:S, function(x){mcmc.output$beta[subsamples[x,1],,subsamples[x,2]]}))
  res$ode <- t(sapply(1:S, function(x){mcmc.output$ode[subsamples[x,1],,subsamples[x,2]]}))
  if (length(dim(mcmc.output$rr)) > 2){
    res$rr <- t(sapply(1:S, function(x){mcmc.output$rr[subsamples[x,1],,subsamples[x,2]]}))
  } else {
    res$rr <- t(sapply(1:S, function(x){mcmc.output$rr[subsamples[x,1],subsamples[x,2]]}))
  }
  res$lik.params <- t(sapply(1:S, function(x){mcmc.output$lik.params[subsamples[x,1],,subsamples[x,2]]}))
  res$loglik <- unlist(sapply(1:S, function(x){mcmc.output$loglik[subsamples[x,1],1,subsamples[x,2]]}))
  
  if (mcmc.output$sf.choice){
    res$sf.vals <- t(sapply(1:S, function(x){mcmc.output$sf.vals[subsamples[x,1],subsamples[x,2]]}))
  }
  
  res$days <- mcmc.output$days
  res$Z.beta <- mcmc.output$Z.beta
  res$Z.rr <- mcmc.output$Z.rr
  colnames(res$ode) <- mcmc.output$col.names
  res$df <- mcmc.output$df
  res$loc <- mcmc.output$loc
  res$const <- mcmc.output$const
  res$introday <- mcmc.output$introday
  res$sf.choice <- mcmc.output$sf.choice
  
  return(res)
}

mcmc.params <- function(out.folder, burnin = 0){
  # Input:
  #   out.folder: location of posterior output.
  #   burnin: number of burn-in samples to leave out.
  # Output:
  #   All MCMC parameters saved after burnin.
  
  # save output files in a list
  files.list <- list.files(out.folder,"*.Rdata")
  files.list <- files.list[order(nchar(files.list), files.list)]
  max.iter=length(files.list)
  
  # load mcmc output
  load(paste(out.folder, files.list[1], sep=""))
  data <- out
  
  # initialize parameters
  n.chains <- length(data)
  n.mcmc <- nrow(data[[1]]$beta)
  n.beta <- ncol(data[[1]]$beta)
  n.ode.params <- ncol(data[[1]]$ode.params)
  n.rr.params <- ncol(data[[1]]$rr.params)
  n.lik.params <- ncol(data[[1]]$lik.params)
  spline <- data[[1]]$spline
  df <- data[[1]]$df
  
  if (length(data[[1]]$sf.choice) > 0){
    sf.choice <- data[[1]]$sf.choice
  } else {
    sf.choice <- FALSE
  }
  
  
  # data structures
  beta.chains=array(NA,dim=c(n.mcmc*max.iter,n.beta,n.chains))
  ode.chains=array(NA,dim=c(n.mcmc*max.iter,n.ode.params,n.chains))
  rr.chains=array(NA,dim=c(n.mcmc*max.iter,n.rr.params,n.chains))
  lik.chains=array(NA,dim=c(n.mcmc*max.iter,n.lik.params,n.chains))
  lik.vals=array(NA, c(n.mcmc*max.iter,n.lik.params,n.chains))
  
  if (sf.choice){
    sf.vals = matrix(NA, nrow = n.mcmc * max.iter, ncol = n.chains)
  }
  
  for(iter in 1:max.iter){
    ## load in data and call it "data"
    load(paste(out.folder,files.list[iter],sep=""))
    data=out
    ## get betas and ode.params
    for(i in 1:n.chains){
      beta.chains[(iter-1)*n.mcmc+(1:n.mcmc),,i] <- data[[i]]$beta
      ode.chains[(iter-1)*n.mcmc+(1:n.mcmc),,i] <- data[[i]]$ode.params
      rr.chains[(iter-1)*n.mcmc+(1:n.mcmc),,i] <- data[[i]]$rr.params
      lik.chains[(iter-1)*n.mcmc+(1:n.mcmc),,i] <- data[[i]]$lik.params
      lik.vals[(iter-1)*n.mcmc+(1:n.mcmc),,i] <- data[[i]]$loglik
      
      if (sf.choice){
        sf.vals[(iter-1) * n.mcmc + (1:n.mcmc), i] = data[[i]]$sf.vals
      }
    }
  }
  
  # return list of values
  results <- list()
  if (burnin > 0){
    results$beta <- beta.chains[-c(1:burnin),,]
    results$ode <- ode.chains[-c(1:burnin),,]
    results$rr <- rr.chains[-c(1:burnin),,]
    results$lik.params <- lik.chains[-c(1:burnin),,]
    results$loglik <- lik.vals[-c(1:burnin),,]
    if (sf.choice){
      results$sf.vals <- sf.vals[-c(1:burnin),]
    }
  } else {
    results$beta <- beta.chains
    results$ode <- ode.chains
    results$rr <- rr.chains
    results$lik.params <- lik.chains
    results$loglik <- lik.vals
    if (sf.choice){
      results$sf.vals <- sf.vals
    }
  }
  
  library(fda)
  df=data[[1]]$df
  results$days <- df$daynum
  results$Z.beta <- eval.basis(61:max(results$days), data[[1]]$spline.beta)
  
  if (is.matrix(data[[1]]$spline.rr)){
    results$Z.rr <- data[[1]]$spline.rr[-nrow(data[[1]]$spline.rr),]
  } else {
    results$Z.rr <- eval.basis(61:max(results$days),data[[1]]$spline.rr) 
  }
  
  results$col.names <- colnames(data[[1]]$ode.params)
  results$df <- df
  results$loc <- data[[1]]$loc
  results$const <- data[[1]]$const.params
  results$introday <- data[[1]]$introday
  results$sf.choice <- sf.choice
  
  return(results)
}

### B. Generating trajectories

traj.sim <- function(samples, odepath, csv = FALSE){
  
  # number of samples available for traj sim
  n.samples <- nrow(samples$betas)
  
  for (k in 1:n.samples){
    
    # last day of simulated output
    end.day <- max(samples$days) 
    
    if (csv){
      n.b <- ncol(samples$betas)
      beta.daily <- as.vector(samples$betas[k,])
    } else {
      # beta daily vector
      beta.daily <- samples$Z.beta %*% samples$betas[k,]
    }
    
    if (k == 1){
      beta.full <- matrix(NA, nrow = n.samples, ncol = length(beta.daily))
      beta.full[k,] <- beta.daily
    } else {
      beta.full[k,] <- beta.daily
    }
    
    # check for fitted hosp.report.rate parameter
    if(is.element("hosp.report.rate", colnames(samples$ode))){
      idx.h <- which(colnames(samples$ode) == "hosp.report.rate")
      cumul.hosp.rr <- samples$ode[k, idx.h]
      ode.params <- samples$ode[k, -idx.h]
      names(ode.params) <- colnames(samples$ode)[-idx.h]
      constants <- samples$const
    } else {
      cumul.hosp.rr = 1
      ode.params <- unlist(samples$ode[k,])
      names(ode.params) <- colnames(samples$ode)
      constants <- samples$const
    }
    
    if(is.element("p.asympt", colnames(samples$ode))){
      non.odesim <- c("p.asympt")
    } else {
      non.odesim  <- NULL
    }
    
    if (samples$sf.choice){
      
      sf.choice.names <- c(" ", "-symp-frac-davies", "-symp-frac-equal 0.3",
                           "-symp-frac-equal 0.4", "-symp-frac-equal 0.5",
                           "-symp-frac-equal 0.6", "-symp-frac-equal 0.7")
      # generate trajectory
      traj.k <- traj.from.params(beta.daily, 
                               params=ode.params, 
                               tf = end.day, 
                               introday = samples$introday,
                               const.params = constants,
                               non.odesim.params = non.odesim,
                               odepath=odepath,
                               loc= samples$loc,
                               symp = sf.choice.names[samples$sf.vals[k]])
    } else {
      # generate trajectory
      traj.k <- traj.from.params(beta = beta.daily, 
                                 params = ode.params, 
                                 const.params = constants,
                                 non.odesim.params = non.odesim,
                                 introday = samples$introday,
                                 tf = end.day,
                                 odepath = odepath,
                                 loc = samples$loc,
                                 symp = NULL)
    }
    
    tp.k <- traj.fine.process(traj.k, loc = samples$loc)
    
    # create reporting rate vector of correct length
    
    if (!is.matrix(samples$Z.rr)){
      rr.full <- as.vector(samples$Z.rr * samples$rr[k])
    } else if (ncol(samples$Z.rr) == 1){
      rr.full <- as.vector(samples$Z.rr * samples$rr[k])
    } else {
      rr.full <- as.vector(samples$Z.rr %*% samples$rr[k,])
    }
    
    if (k == 1){
      rr.full.m <- matrix(NA, nrow = n.samples, ncol = length(rr.full))
      rr.full.m[k,] <- rr.full
    } else {
      rr.full.m[k,] <- rr.full
    }
    
    lower.rr <- rr.full[1]
    upper.rr <- rr.full[length(rr.full)]
    
    lower.rep <- which(tp.k$days.odesim < min(samples$df$daynum, na.rm = T))
    upper.rep <- which(tp.k$days.odesim > max(samples$df$daynum, na.rm = T))
    
    rr.full <- c(rep(lower.rr, length(lower.rep)),
                 rr.full,
                 rep(upper.rr, length(upper.rep)))
    
    
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
      
    # symptomatic cases (no age)
    tot.sympt.new[k,] <- rr.full * tp.k$tot.sympt.new.odesim
    tot.sympt.cum[k,] <- cumsum(tot.sympt.new[k,])
      
    # hospitalizations (no age)
    tot.hosp.new[k,] <- cumul.hosp.rr * tp.k$tot.hosp.new.odesim 
    tot.hosp.cum[k,] <- cumsum(tot.hosp.new[k,])
      
    # deaths (no age)
    tot.deaths.new[k,] <- tp.k$tot.deaths.new.odesim
    tot.deaths.cum[k,] <- tp.k$tot.deaths.cum.odesim
    tot.hosp.deaths.new[k,] <- tp.k$tot.hosp.deaths.new.odesim
    tot.home.deaths.new[k,] <- tp.k$tot.home.deaths.new.odesim
    tot.hosp.deaths.cum[k,] <- tp.k$tot.hosp.deaths.cum.odesim
    tot.home.deaths.cum[k,] <- tp.k$tot.home.deaths.cum.odesim
      
    # current (no age)
    tot.hosp.curr.noicu[k,] <- tp.k$tot.hosp.curr.noicu.odesim
    tot.icu.curr.novent[k,] <- tp.k$tot.icu.curr.novent.odesim
    tot.vent.curr[k,] <- tp.k$tot.vent.curr.odesim
    tot.hosp.curr[k,] <- tp.k$tot.hosp.curr.odesim
    tot.icu.curr[k,] <- tp.k$tot.icu.curr.odesim
    
    # discharges
    tot.hosp.discharges.new[k,] <- tp.k$tot.hosp.discharges.new.odesim
    tot.hosp.discharges.cum[k,] <- tp.k$tot.hosp.discharges.cum.odesim
     
    if (samples$loc == "MA"){
      
      ### symptomatic
      
      # new symptomatic cases (by age)
      sympt.new.01[k,] <- rr.full * tp.k$sympt.new.odesim[,1] 
      sympt.new.2[k,] <- rr.full * tp.k$sympt.new.odesim[,2] 
      sympt.new.3[k,] <- rr.full * tp.k$sympt.new.odesim[,3] 
      sympt.new.4[k,] <- rr.full * tp.k$sympt.new.odesim[,4] 
      sympt.new.5[k,] <- rr.full * tp.k$sympt.new.odesim[,5] 
      sympt.new.6[k,] <- rr.full * tp.k$sympt.new.odesim[,6] 
      sympt.new.7[k,] <- rr.full * tp.k$sympt.new.odesim[,7] 
      sympt.new.8[k,] <- rr.full * tp.k$sympt.new.odesim[,8] 
      
      # cum. symptomatic cases (by age)
      #sympt.cum.01[k,] <- rr.full * tp.k$sympt.cum.odesim[,1]
      # sympt.cum.2[k,] <- rr.full * tp.k$sympt.cum.odesim[,2] 
      # sympt.cum.3[k,] <- rr.full * tp.k$sympt.cum.odesim[,3] 
      # sympt.cum.4[k,] <- rr.full * tp.k$sympt.cum.odesim[,4] 
      # sympt.cum.5[k,] <- rr.full * tp.k$sympt.cum.odesim[,5] 
      # sympt.cum.6[k,] <- rr.full * tp.k$sympt.cum.odesim[,6] 
      # sympt.cum.7[k,] <- rr.full * tp.k$sympt.cum.odesim[,7] 
      # sympt.cum.8[k,] <- rr.full * tp.k$sympt.cum.odesim[,8]
      
      sympt.cum.01[k,] <- cumsum(sympt.new.01[k,])
      sympt.cum.2[k,] <- cumsum(sympt.new.2[k,])
      sympt.cum.3[k,] <- cumsum(sympt.new.3[k,])
      sympt.cum.4[k,] <- cumsum(sympt.new.4[k,]) 
      sympt.cum.5[k,] <- cumsum(sympt.new.5[k,])
      sympt.cum.6[k,] <- cumsum(sympt.new.6[k,]) 
      sympt.cum.7[k,] <- cumsum(sympt.new.7[k,])
      sympt.cum.8[k,] <- cumsum(sympt.new.8[k,])
      
      ### hospitalizations 
      
      # new hospitalization (by age)
      hosp.new.01[k,] <- cumul.hosp.rr * tp.k$hosp.new.odesim[,1]
      hosp.new.2[k,] <- cumul.hosp.rr * tp.k$hosp.new.odesim[,2]
      hosp.new.3[k,] <- cumul.hosp.rr * tp.k$hosp.new.odesim[,3]
      hosp.new.4[k,] <- cumul.hosp.rr * tp.k$hosp.new.odesim[,4]
      hosp.new.5[k,] <- cumul.hosp.rr * tp.k$hosp.new.odesim[,5]
      hosp.new.6[k,] <- cumul.hosp.rr * tp.k$hosp.new.odesim[,6]
      hosp.new.7[k,] <- cumul.hosp.rr * tp.k$hosp.new.odesim[,7]
      hosp.new.8[k,] <- cumul.hosp.rr * tp.k$hosp.new.odesim[,8]
      
      # cum. hospitalizations (by age)
      hosp.cum.01[k,] <- cumsum(hosp.new.01[k,])
      hosp.cum.2[k,] <- cumsum(hosp.new.2[k,]) 
      hosp.cum.3[k,] <- cumsum(hosp.new.3[k,]) 
      hosp.cum.4[k,] <- cumsum(hosp.new.4[k,]) 
      hosp.cum.5[k,] <- cumsum(hosp.new.5[k,])  
      hosp.cum.6[k,] <- cumsum(hosp.new.6[k,]) 
      hosp.cum.7[k,] <- cumsum(hosp.new.7[k,])  
      hosp.cum.8[k,] <- cumsum(hosp.new.8[k,])  
      
      ### deaths
      
      ### new deaths
      deaths.new.01[k,] <- tp.k$deaths.new.odesim[,1]
      deaths.new.2[k,] <- tp.k$deaths.new.odesim[,2]
      deaths.new.3[k,] <- tp.k$deaths.new.odesim[,3]
      deaths.new.4[k,] <- tp.k$deaths.new.odesim[,4]
      deaths.new.5[k,] <- tp.k$deaths.new.odesim[,5]
      deaths.new.6[k,] <- tp.k$deaths.new.odesim[,6]
      deaths.new.7[k,] <- tp.k$deaths.new.odesim[,7]
      deaths.new.8[k,] <- tp.k$deaths.new.odesim[,8]
      
      ### cum. deaths
      deaths.cum.01[k,] <- tp.k$deaths.cum.odesim[,1]
      deaths.cum.2[k,] <- tp.k$deaths.cum.odesim[,2]
      deaths.cum.3[k,] <- tp.k$deaths.cum.odesim[,3]
      deaths.cum.4[k,] <- tp.k$deaths.cum.odesim[,4]
      deaths.cum.5[k,] <- tp.k$deaths.cum.odesim[,5]
      deaths.cum.6[k,] <- tp.k$deaths.cum.odesim[,6]
      deaths.cum.7[k,] <- tp.k$deaths.cum.odesim[,7]
      deaths.cum.8[k,] <- tp.k$deaths.cum.odesim[,8]
      
      ### hosp. deaths
      hosp.deaths.new.01[k,] <- tp.k$hosp.deaths.new.odesim[,1]
      hosp.deaths.new.2[k,] <- tp.k$hosp.deaths.new.odesim[,2]
      hosp.deaths.new.3[k,] <- tp.k$hosp.deaths.new.odesim[,3]
      hosp.deaths.new.4[k,] <- tp.k$hosp.deaths.new.odesim[,4]
      hosp.deaths.new.5[k,] <- tp.k$hosp.deaths.new.odesim[,5]
      hosp.deaths.new.6[k,] <- tp.k$hosp.deaths.new.odesim[,6]
      hosp.deaths.new.7[k,] <- tp.k$hosp.deaths.new.odesim[,7]
      hosp.deaths.new.8[k,] <- tp.k$hosp.deaths.new.odesim[,8]
      
      ### home deaths
      home.deaths.new.01[k,] <- tp.k$home.deaths.new.odesim[,1]
      home.deaths.new.2[k,] <- tp.k$home.deaths.new.odesim[,2]
      home.deaths.new.3[k,] <- tp.k$home.deaths.new.odesim[,3]
      home.deaths.new.4[k,] <- tp.k$home.deaths.new.odesim[,4]
      home.deaths.new.5[k,] <- tp.k$home.deaths.new.odesim[,5]
      home.deaths.new.6[k,] <- tp.k$home.deaths.new.odesim[,6]
      home.deaths.new.7[k,] <- tp.k$home.deaths.new.odesim[,7]
      home.deaths.new.8[k,] <- tp.k$home.deaths.new.odesim[,8]
      
      ### discharges
      
      ### new hosp discharges
      hosp.discharges.new.01[k,] <- tp.k$hosp.discharges.new.odesim[,1]
      hosp.discharges.new.2[k,] <- tp.k$hosp.discharges.new.odesim[,2]
      hosp.discharges.new.3[k,] <- tp.k$hosp.discharges.new.odesim[,3]
      hosp.discharges.new.4[k,] <- tp.k$hosp.discharges.new.odesim[,4]
      hosp.discharges.new.5[k,] <- tp.k$hosp.discharges.new.odesim[,5]
      hosp.discharges.new.6[k,] <- tp.k$hosp.discharges.new.odesim[,6]
      hosp.discharges.new.7[k,] <- tp.k$hosp.discharges.new.odesim[,7]
      hosp.discharges.new.8[k,] <- tp.k$hosp.discharges.new.odesim[,8]
      
      ### cum hosp discharges 
      hosp.discharges.cum.01[k,] <- tp.k$hosp.discharges.cum.odesim[,1]
      hosp.discharges.cum.2[k,] <- tp.k$hosp.discharges.cum.odesim[,2]
      hosp.discharges.cum.3[k,] <- tp.k$hosp.discharges.cum.odesim[,3]
      hosp.discharges.cum.4[k,] <- tp.k$hosp.discharges.cum.odesim[,4]
      hosp.discharges.cum.5[k,] <- tp.k$hosp.discharges.cum.odesim[,5]
      hosp.discharges.cum.6[k,] <- tp.k$hosp.discharges.cum.odesim[,6]
      hosp.discharges.cum.7[k,] <- tp.k$hosp.discharges.cum.odesim[,7]
      hosp.discharges.cum.8[k,] <- tp.k$hosp.discharges.cum.odesim[,8]
      
      
    } else {
      
      ### symptomatic
      
      # new symptomatic cases (by age)
      sympt.new.0[k,] <- rr.full * tp.k$sympt.new.odesim[,1] 
      sympt.new.1[k,] <- rr.full * tp.k$sympt.new.odesim[,2] 
      sympt.new.2[k,] <- rr.full * tp.k$sympt.new.odesim[,3] 
      sympt.new.3[k,] <- rr.full * tp.k$sympt.new.odesim[,4] 
      sympt.new.4[k,] <- rr.full * tp.k$sympt.new.odesim[,5] 
      sympt.new.5[k,] <- rr.full * tp.k$sympt.new.odesim[,6] 
      sympt.new.6[k,] <- rr.full * tp.k$sympt.new.odesim[,7] 
      sympt.new.7[k,] <- rr.full * tp.k$sympt.new.odesim[,8] 
      sympt.new.8[k,] <- rr.full * tp.k$sympt.new.odesim[,9] 
      
      # cum. symptomatic cases (by age)
      # sympt.cum.0[k,] <- rr.full * tp.k$sympt.cum.odesim[,1] 
      # sympt.cum.1[k,] <- rr.full * tp.k$sympt.cum.odesim[,2]
      # sympt.cum.2[k,] <- rr.full * tp.k$sympt.cum.odesim[,3] 
      # sympt.cum.3[k,] <- rr.full * tp.k$sympt.cum.odesim[,4] 
      # sympt.cum.4[k,] <- rr.full * tp.k$sympt.cum.odesim[,5] 
      # sympt.cum.5[k,] <- rr.full * tp.k$sympt.cum.odesim[,6] 
      # sympt.cum.6[k,] <- rr.full * tp.k$sympt.cum.odesim[,7] 
      # sympt.cum.7[k,] <- rr.full * tp.k$sympt.cum.odesim[,8] 
      # sympt.cum.8[k,] <- rr.full * tp.k$sympt.cum.odesim[,9] 
      
      sympt.cum.0[k,] <- cumsum(sympt.new.0[k,])
      sympt.cum.1[k,] <- cumsum(sympt.new.1[k,])
      sympt.cum.2[k,] <- cumsum(sympt.new.2[k,])
      sympt.cum.3[k,] <- cumsum(sympt.new.3[k,])
      sympt.cum.4[k,] <- cumsum(sympt.new.4[k,]) 
      sympt.cum.5[k,] <- cumsum(sympt.new.5[k,])
      sympt.cum.6[k,] <- cumsum(sympt.new.6[k,]) 
      sympt.cum.7[k,] <- cumsum(sympt.new.7[k,])
      sympt.cum.8[k,] <- cumsum(sympt.new.8[k,])
      
      ### hospitalizations 
      
      # new hospitalization (by age)
      hosp.new.0[k,] <- cumul.hosp.rr * tp.k$hosp.new.odesim[,1]
      hosp.new.1[k,] <- cumul.hosp.rr * tp.k$hosp.new.odesim[,2]
      hosp.new.2[k,] <- cumul.hosp.rr * tp.k$hosp.new.odesim[,3]
      hosp.new.3[k,] <- cumul.hosp.rr * tp.k$hosp.new.odesim[,4]
      hosp.new.4[k,] <- cumul.hosp.rr * tp.k$hosp.new.odesim[,5]
      hosp.new.5[k,] <- cumul.hosp.rr * tp.k$hosp.new.odesim[,6]
      hosp.new.6[k,] <- cumul.hosp.rr * tp.k$hosp.new.odesim[,7]
      hosp.new.7[k,] <- cumul.hosp.rr * tp.k$hosp.new.odesim[,8]
      hosp.new.8[k,] <- cumul.hosp.rr * tp.k$hosp.new.odesim[,9]
      
      # cum. hospitalizations (by age)
      # hosp.cum.0[k,] <- cumul.hosp.rr * tp.k$hosp.cum.odesim[,1]
      # hosp.cum.1[k,] <- cumul.hosp.rr * tp.k$hosp.cum.odesim[,2]
      # hosp.cum.2[k,] <- cumul.hosp.rr * tp.k$hosp.cum.odesim[,3] 
      # hosp.cum.3[k,] <- cumul.hosp.rr * tp.k$hosp.cum.odesim[,4] 
      # hosp.cum.4[k,] <- cumul.hosp.rr * tp.k$hosp.cum.odesim[,5] 
      # hosp.cum.5[k,] <- cumul.hosp.rr * tp.k$hosp.cum.odesim[,6] 
      # hosp.cum.6[k,] <- cumul.hosp.rr * tp.k$hosp.cum.odesim[,7] 
      # hosp.cum.7[k,] <- cumul.hosp.rr * tp.k$hosp.cum.odesim[,8] 
      # hosp.cum.8[k,] <- cumul.hosp.rr * tp.k$hosp.cum.odesim[,9] 
      
      hosp.cum.0[k,] <- cumsum(hosp.new.0[k,])
      hosp.cum.1[k,] <- cumsum(hosp.new.1[k,])
      hosp.cum.2[k,] <- cumsum(hosp.new.2[k,])
      hosp.cum.3[k,] <- cumsum(hosp.new.3[k,])
      hosp.cum.4[k,] <- cumsum(hosp.new.4[k,])
      hosp.cum.5[k,] <- cumsum(hosp.new.5[k,])
      hosp.cum.6[k,] <- cumsum(hosp.new.6[k,])
      hosp.cum.7[k,] <- cumsum(hosp.new.7[k,])
      hosp.cum.8[k,] <- cumsum(hosp.new.8[k,])
      
      ### deaths
      
      ### new deaths
      deaths.new.0[k,] <- tp.k$deaths.new.odesim[,1]
      deaths.new.1[k,] <- tp.k$deaths.new.odesim[,2]
      deaths.new.2[k,] <- tp.k$deaths.new.odesim[,3]
      deaths.new.3[k,] <- tp.k$deaths.new.odesim[,4]
      deaths.new.4[k,] <- tp.k$deaths.new.odesim[,5]
      deaths.new.5[k,] <- tp.k$deaths.new.odesim[,6]
      deaths.new.6[k,] <- tp.k$deaths.new.odesim[,7]
      deaths.new.7[k,] <- tp.k$deaths.new.odesim[,8]
      deaths.new.8[k,] <- tp.k$deaths.new.odesim[,9]
      
      ### cum. deaths
      deaths.cum.0[k,] <- tp.k$deaths.cum.odesim[,1]
      deaths.cum.1[k,] <- tp.k$deaths.cum.odesim[,2]
      deaths.cum.2[k,] <- tp.k$deaths.cum.odesim[,3]
      deaths.cum.3[k,] <- tp.k$deaths.cum.odesim[,4]
      deaths.cum.4[k,] <- tp.k$deaths.cum.odesim[,5]
      deaths.cum.5[k,] <- tp.k$deaths.cum.odesim[,6]
      deaths.cum.6[k,] <- tp.k$deaths.cum.odesim[,7]
      deaths.cum.7[k,] <- tp.k$deaths.cum.odesim[,8]
      deaths.cum.8[k,] <- tp.k$deaths.cum.odesim[,9]
      
      ### hosp. deaths
      hosp.deaths.new.0[k,] <- tp.k$hosp.deaths.new.odesim[,1]
      hosp.deaths.new.1[k,] <- tp.k$hosp.deaths.new.odesim[,2]
      hosp.deaths.new.2[k,] <- tp.k$hosp.deaths.new.odesim[,3]
      hosp.deaths.new.3[k,] <- tp.k$hosp.deaths.new.odesim[,4]
      hosp.deaths.new.4[k,] <- tp.k$hosp.deaths.new.odesim[,5]
      hosp.deaths.new.5[k,] <- tp.k$hosp.deaths.new.odesim[,6]
      hosp.deaths.new.6[k,] <- tp.k$hosp.deaths.new.odesim[,7]
      hosp.deaths.new.7[k,] <- tp.k$hosp.deaths.new.odesim[,8]
      hosp.deaths.new.8[k,] <- tp.k$hosp.deaths.new.odesim[,9]
      
      ### home deaths
      home.deaths.new.0[k,] <- tp.k$home.deaths.new.odesim[,1]
      home.deaths.new.1[k,] <- tp.k$home.deaths.new.odesim[,2]
      home.deaths.new.2[k,] <- tp.k$home.deaths.new.odesim[,3]
      home.deaths.new.3[k,] <- tp.k$home.deaths.new.odesim[,4]
      home.deaths.new.4[k,] <- tp.k$home.deaths.new.odesim[,5]
      home.deaths.new.5[k,] <- tp.k$home.deaths.new.odesim[,6]
      home.deaths.new.6[k,] <- tp.k$home.deaths.new.odesim[,7]
      home.deaths.new.7[k,] <- tp.k$home.deaths.new.odesim[,8]
      home.deaths.new.8[k,] <- tp.k$home.deaths.new.odesim[,9]
      
      ### discharges
      
      ### new hosp discharges
      hosp.discharges.new.0[k,] <- tp.k$hosp.discharges.new.odesim[,1]
      hosp.discharges.new.1[k,] <- tp.k$hosp.discharges.new.odesim[,2]
      hosp.discharges.new.2[k,] <- tp.k$hosp.discharges.new.odesim[,3]
      hosp.discharges.new.3[k,] <- tp.k$hosp.discharges.new.odesim[,4]
      hosp.discharges.new.4[k,] <- tp.k$hosp.discharges.new.odesim[,5]
      hosp.discharges.new.5[k,] <- tp.k$hosp.discharges.new.odesim[,6]
      hosp.discharges.new.6[k,] <- tp.k$hosp.discharges.new.odesim[,7]
      hosp.discharges.new.7[k,] <- tp.k$hosp.discharges.new.odesim[,8]
      hosp.discharges.new.8[k,] <- tp.k$hosp.discharges.new.odesim[,9]
      
      ### cum hosp discharges 
      hosp.discharges.cum.0[k,] <- tp.k$hosp.discharges.cum.odesim[,1]
      hosp.discharges.cum.1[k,] <- tp.k$hosp.discharges.cum.odesim[,2]
      hosp.discharges.cum.2[k,] <- tp.k$hosp.discharges.cum.odesim[,3]
      hosp.discharges.cum.3[k,] <- tp.k$hosp.discharges.cum.odesim[,4]
      hosp.discharges.cum.4[k,] <- tp.k$hosp.discharges.cum.odesim[,5]
      hosp.discharges.cum.5[k,] <- tp.k$hosp.discharges.cum.odesim[,6]
      hosp.discharges.cum.6[k,] <- tp.k$hosp.discharges.cum.odesim[,7]
      hosp.discharges.cum.7[k,] <- tp.k$hosp.discharges.cum.odesim[,8]
      hosp.discharges.cum.8[k,] <- tp.k$hosp.discharges.cum.odesim[,9]
      
    }  
  }
    
    
  ### create lists for age-stratified results
  if (samples$loc == "MA"){
    
    sympt.new.age <- list(sympt.new.01,
                          sympt.new.2,
                          sympt.new.3,
                          sympt.new.4,
                          sympt.new.5,
                          sympt.new.6,
                          sympt.new.7,
                          sympt.new.8)
    
    sympt.cum.age <- list(sympt.cum.01,
                          sympt.cum.2,
                          sympt.cum.3,
                          sympt.cum.4,
                          sympt.cum.5,
                          sympt.cum.6,
                          sympt.cum.7,
                          sympt.cum.8)
    
    hosp.new.age <- list(hosp.new.01,
                         hosp.new.2,
                         hosp.new.3,
                         hosp.new.4,
                         hosp.new.5,
                         hosp.new.6,
                         hosp.new.7,
                         hosp.new.8)
    
    hosp.cum.age <- list(hosp.cum.01,
                         hosp.cum.2,
                         hosp.cum.3,
                         hosp.cum.4,
                         hosp.cum.5,
                         hosp.cum.6,
                         hosp.cum.7,
                         hosp.cum.8)
    
    deaths.new.age <- list(deaths.new.01,
                           deaths.new.2,
                           deaths.new.3,
                           deaths.new.4,
                           deaths.new.5,
                           deaths.new.6,
                           deaths.new.7,
                           deaths.new.8)
    
    
    deaths.cum.age <- list(deaths.cum.01,
                           deaths.cum.2,
                           deaths.cum.3,
                           deaths.cum.4,
                           deaths.cum.5,
                           deaths.cum.6,
                           deaths.cum.7,
                           deaths.cum.8)
    
    hosp.deaths.new.age <- list(hosp.deaths.new.01,
                                hosp.deaths.new.2,
                                hosp.deaths.new.3,
                                hosp.deaths.new.4,
                                hosp.deaths.new.5,
                                hosp.deaths.new.6,
                                hosp.deaths.new.7,
                                hosp.deaths.new.8)
    
    
    home.deaths.new.age <- list(home.deaths.new.01,
                                home.deaths.new.2,
                                home.deaths.new.3,
                                home.deaths.new.4,
                                home.deaths.new.5,
                                home.deaths.new.6,
                                home.deaths.new.7,
                                home.deaths.new.8)
    
    hosp.discharges.new.age <- list(hosp.discharges.new.01,
                                    hosp.discharges.new.2,
                                    hosp.discharges.new.3,
                                    hosp.discharges.new.4,
                                    hosp.discharges.new.5,
                                    hosp.discharges.new.6,
                                    hosp.discharges.new.7,
                                    hosp.discharges.new.8)
    
    
    hosp.discharges.cum.age <- list(hosp.discharges.cum.01,
                                    hosp.discharges.cum.2,
                                    hosp.discharges.cum.3,
                                    hosp.discharges.cum.4,
                                    hosp.discharges.cum.5,
                                    hosp.discharges.cum.6,
                                    hosp.discharges.cum.7,
                                    hosp.discharges.cum.8)
  } else {
    
    sympt.new.age <- list(sympt.new.0,
                          sympt.new.1,
                          sympt.new.2,
                          sympt.new.3,
                          sympt.new.4,
                          sympt.new.5,
                          sympt.new.6,
                          sympt.new.7,
                          sympt.new.8)
    
    sympt.cum.age <- list(sympt.cum.0,
                          sympt.cum.1,
                          sympt.cum.2,
                          sympt.cum.3,
                          sympt.cum.4,
                          sympt.cum.5,
                          sympt.cum.6,
                          sympt.cum.7,
                          sympt.cum.8)
    
    hosp.new.age <- list(hosp.new.0,
                         hosp.new.1,
                         hosp.new.2,
                         hosp.new.3,
                         hosp.new.4,
                         hosp.new.5,
                         hosp.new.6,
                         hosp.new.7,
                         hosp.new.8)
    
    hosp.cum.age <- list(hosp.cum.0,
                         hosp.cum.1,
                         hosp.cum.2,
                         hosp.cum.3,
                         hosp.cum.4,
                         hosp.cum.5,
                         hosp.cum.6,
                         hosp.cum.7,
                         hosp.cum.8)
    
    deaths.new.age <- list(deaths.new.0,
                           deaths.new.1,
                           deaths.new.2,
                           deaths.new.3,
                           deaths.new.4,
                           deaths.new.5,
                           deaths.new.6,
                           deaths.new.7,
                           deaths.new.8)
    
    
    deaths.cum.age <- list(deaths.cum.0,
                           deaths.cum.1,
                           deaths.cum.2,
                           deaths.cum.3,
                           deaths.cum.4,
                           deaths.cum.5,
                           deaths.cum.6,
                           deaths.cum.7,
                           deaths.cum.8)
    
    hosp.deaths.new.age <- list(hosp.deaths.new.0,
                                hosp.deaths.new.1,
                                hosp.deaths.new.2,
                                hosp.deaths.new.3,
                                hosp.deaths.new.4,
                                hosp.deaths.new.5,
                                hosp.deaths.new.6,
                                hosp.deaths.new.7,
                                hosp.deaths.new.8)
    
    
    home.deaths.new.age <- list(home.deaths.new.0,
                                home.deaths.new.1,
                                home.deaths.new.2,
                                home.deaths.new.3,
                                home.deaths.new.4,
                                home.deaths.new.5,
                                home.deaths.new.6,
                                home.deaths.new.7,
                                home.deaths.new.8)
  
    hosp.discharges.new.age <- list(hosp.discharges.new.0,
                                    hosp.discharges.new.1,
                                    hosp.discharges.new.2,
                                    hosp.discharges.new.3,
                                    hosp.discharges.new.4,
                                    hosp.discharges.new.5,
                                    hosp.discharges.new.6,
                                    hosp.discharges.new.7,
                                    hosp.discharges.new.8)
  
  
    hosp.discharges.cum.age <- list(hosp.discharges.cum.0,
                                    hosp.discharges.cum.1,
                                    hosp.discharges.cum.2,
                                    hosp.discharges.cum.3,
                                    hosp.discharges.cum.4,
                                    hosp.discharges.cum.5,
                                    hosp.discharges.cum.6,
                                    hosp.discharges.cum.7,
                                    hosp.discharges.cum.8)
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
  
  return(traj.vals)
}

traj.fine.process <- function(traj, loc){
  # Input:
  #   traj: simulated trajectory output.
  #   loc: location (either PA, MA, or RI).
  # Output:
  #   Parsed trajectory output, with case statistics of interest.
  
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
  traj.sympt <- traj[ ,sympt.cum.idx]
  
  # combine 0-9 and 10-19 age classes for MA
  if(loc=="MA"){
    traj.sympt.ma <- traj.sympt[,2:ncol(traj.sympt)]
    traj.sympt.ma[,1] <- traj.sympt[,1] + traj.sympt[,2]
    traj.sympt <- traj.sympt.ma
  }
  
  ## calculating new cases for each day
  traj.sympt.new = traj.sympt
  traj.sympt.new[-1,] = traj.sympt[-1,] - traj.sympt[-nr,]
  tot.sympt.new=rowSums(traj.sympt.new)
  
  
  ###########################
  ### 2. Hospitalizations ###
  ###########################
  
  ## cumulative hospitalized by age:
  traj.hosp.cum=traj[,hosp.cum.idx]
  
  ## calculating new hosp from traj for each day
  traj.hosp.new = traj.hosp.cum
  traj.hosp.new[-1,] = traj.hosp.cum[-1,] - traj.hosp.cum[-nr,]
  
  if(loc=="MA"){
    traj.hosp.new.ma = traj.hosp.new[,2:ncol(traj.hosp.new)]
    traj.hosp.new.ma[,1] = traj.hosp.new[,1] + traj.hosp.new[,2]
    traj.hosp.new = traj.hosp.new.ma
    
    traj.hosp.cum.ma = traj.hosp.cum[,2:ncol(traj.hosp.cum)]
    traj.hosp.cum.ma[,1] = traj.hosp.cum[,1] + traj.hosp.cum[,2]
    traj.hosp.cum = traj.hosp.cum.ma
  }
  
  tot.hosp.new=rowSums(traj.hosp.new)
  
  #######################################
  ### 3. Current Hosps, Vent, and ICU ###
  #######################################
  
  ## for proposing new s2.hosp parameter
  hosp.curr=traj[,hosp.curr.idx]
  tot.hosp.curr=rowSums(hosp.curr)
  
  ## for proposing new s2.icu parameter
  icu.curr=traj[,icu.curr.idx]
  tot.icu.curr=rowSums(icu.curr)
  
  ## for proposing new s2.vent parameter
  vent.curr=traj[,vent.curr.idx]
  tot.vent.curr=rowSums(vent.curr)
  
  #########################################
  ### 4. Deaths (home, hosp, and total) ###
  #########################################
  
  ## deaths
  home.deaths=traj[,home.deaths.idx]
  hosp.deaths=traj[,hosp.deaths.idx]
  
  tot.home.deaths.cum=rowSums(home.deaths, na.rm = T)
  tot.hosp.deaths.cum=rowSums(hosp.deaths, na.rm = T)
  
  deaths.cum=home.deaths+hosp.deaths
  
  #home.deaths.new=rbind(home.deaths[1,],as.matrix(home.deaths[-1,]-home.deaths[-nrow(traj),]))
  home.deaths.new = home.deaths
  home.deaths.new[-1,] = home.deaths[-1,] - home.deaths[-nr,]
  
  #hosp.deaths.new=rbind(hosp.deaths[1,],as.matrix(hosp.deaths[-1,]-hosp.deaths[-nrow(traj),]))
  hosp.deaths.new = hosp.deaths
  hosp.deaths.new[-1,] = hosp.deaths[-1,] - hosp.deaths[-nr,]
  
  # total (new) deaths
  deaths.new = home.deaths.new + hosp.deaths.new
  
  tot.home.deaths.new=rowSums(home.deaths.new)
  tot.hosp.deaths.new=rowSums(hosp.deaths.new)
  tot.deaths.new= tot.home.deaths.new + tot.hosp.deaths.new
  
  # combine 0-9 and 10-19 age classes for MA
  if(loc=="MA"){
    deaths.cum.ma <- deaths.cum[,2:ncol(deaths.cum)]
    deaths.cum.ma[,1] <- deaths.cum[,1] + deaths.cum[,2]
    deaths.cum <- deaths.cum.ma
    
    deaths.new.ma <- deaths.new[,2:ncol(deaths.new)]
    deaths.new.ma[,1] <- deaths.new[,1] + deaths.new[,2]
    deaths.new <- deaths.new.ma
  }
  
  ##############################
  ### 5. Hospital Discharges ###
  ##############################
  
  discharge.idx = 281:289 ## for v5
  hosp.discharges <- traj[,discharge.idx]
  hosp.discharges.new = hosp.discharges
  hosp.discharges.new[-1,] = hosp.discharges[-1,] - hosp.discharges[-nr,]
  
  tot.hosp.discharges = apply(traj[,discharge.idx],1,sum,na.rm = T)
  tot.hosp.discharges.new = tot.hosp.discharges
  tot.hosp.discharges.new[-1] = tot.hosp.discharges[-1] - tot.hosp.discharges[-nr]
  
  ## return
  list(days.odesim = traj[,1],
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
       tot.deaths.cum.odesim = as.numeric(apply(deaths.cum,1,sum,na.rm=TRUE)),
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
       tot.hosp.curr.odesim = tot.hosp.curr+tot.icu.curr+tot.vent.curr,
       tot.icu.curr.odesim = tot.icu.curr+tot.vent.curr,
       # hosp. discharges
       tot.hosp.discharges.cum.odesim = as.numeric(tot.hosp.discharges),
       hosp.discharges.cum.odesim = as.matrix(hosp.discharges),
       tot.hosp.discharges.new.odesim = as.numeric(tot.hosp.discharges.new),
       hosp.discharges.new.odesim = as.matrix(hosp.discharges.new)
  )
}

data.fine.process <- function(df, loc="RI"){
  # Input
  #   df: data frame object.
  #   loc: "RI", "PA", or "MA".
  # Output
  #   Named list with the following objects:
  #     days, cases.cum, cases.new, cases.age.new,
  #     hosp.cum, hosp.new, hosp.age.new, hosp.curr,
  #     icu.curr, vent.curr, deaths.cum, deaths.new
  
  days=df$daynum # days after Jan 1, 2020
  
  ## calculations for age-structured infill and pre-processing of data
  n.days=nrow(df)
  # cases
  cases.cum=df$cumulative_confirmed_cases
  cases.new=c(cases.cum[1],cases.cum[-1]-cases.cum[-n.days])
  if(loc=="MA"){
    cases.age.cum=df[,16:23]
    hosp.age.cum=df[,24:31]
    deaths.age.cum=df[,32:39]
  }else{
    cases.age.cum=df[,16:24]
    hosp.age.cum=df[,25:33]
    deaths.age.cum=df[,34:42]
  }
  cases.age.new=as.matrix(rbind(cases.age.cum[1,],cases.age.cum[-1,]-cases.age.cum[-n.days,]))
  cases.age.new[cases.age.new<0]=0
  # hospitalizations
  hosp.age.new=as.matrix(rbind(hosp.age.cum[1,],hosp.age.cum[-1,]-hosp.age.cum[-n.days,]))
  hosp.age.new[hosp.age.new<0]=0
  hosp.cum=df$hospitalized_cumulative
  hosp.new=c(hosp.cum[1],hosp.cum[-1]-hosp.cum[-n.days])
  hosp.curr=df$hospitalized_currently
  
  # icu
  icu.curr=df$InICU_currently
  # if(loc=="MA"){
  #   icu.curr=df$InICU_Currently
  # }
  # # ventilator
  vent.curr=df$OnVentilator_Currently
  # deaths
  tot.deaths.cum=df$cumulative_deaths
  tot.deaths.new=c(tot.deaths.cum[1],tot.deaths.cum[-1]-tot.deaths.cum[-n.days])
  # ## deaths in hospital
  # if(loc=="PA"){
  #   deaths.hosp.cum=NULL
  #   deaths.hosp.new=NULL
  #   deaths.out.of.hosp.new=NULL
  # }else{
  deaths.hosp.cum=df$Cumulative_hospital_deaths
  deaths.home.cum=df$cumulative_deaths - deaths.hosp.cum
  deaths.hosp.new=c(deaths.hosp.cum[1],deaths.hosp.cum[-1]-deaths.hosp.cum[-n.days])
  deaths.out.of.hosp.new=tot.deaths.new-deaths.hosp.new
  if(is.null(deaths.hosp.cum[1])){
    deaths.hosp.new=NULL
    deaths.out.of.hosp.new=NULL
  }
  
  deaths.age.new=as.matrix(rbind(deaths.age.cum[1,],deaths.age.cum[-1,]-deaths.age.cum[-n.days,]))
  
  # }
  
  ## hospital discharges
  tot.hosp.discharges=df$Cumulative_hospital_discharges
  tot.new.hosp.discharges=c(0,tot.hosp.discharges[-1]-tot.hosp.discharges[-n.days])
  
  ## active surveillance
  if(loc=="MA"){
    tests.as=df$nursinghome_tests_total
    pos.as=df$nursinghome_tests_positive
  }else{
    tests.as=rep(NA,nrow(df))
    pos.as=rep(NA,nrow(df))
  }
  
  list(days=days,
       # symptomatic cases
       tot.sympt.new=cases.new, 
       tot.sympt.cum=cases.cum,
       sympt.new=cases.age.new, 
       sympt.cum=cases.age.cum,
       # hospitalizations
       tot.hosp.new=hosp.new, 
       hosp.new=hosp.age.new, 
       tot.hosp.cum=hosp.cum,
       hosp.cum=hosp.age.cum,
       # current hosp, icu, and vent data
       tot.hosp.curr=hosp.curr, ##note: this is all in hosp (=hosp+icu+vent)
       tot.icu.curr=icu.curr, ## note: this is all in icu (=icu+vent)
       tot.vent.curr=vent.curr, ## current number on ventilators
       # death data
       tot.deaths.new=tot.deaths.new,
       tot.hosp.deaths.new=deaths.hosp.new,
       tot.hosp.deaths.cum=deaths.hosp.cum,
       tot.home.deaths.new=deaths.out.of.hosp.new,
       tot.home.deaths.cum=deaths.home.cum,
       deaths.cum=deaths.age.cum,
       deaths.new=deaths.age.new,
       tot.deaths.cum=tot.deaths.cum,
       # hospital discharges
       tot.hosp.discharges.new=tot.new.hosp.discharges,
       tot.hosp.discharges.cum=tot.hosp.discharges,
       tests.as=tests.as,
       pos.as=pos.as)
}

### C. Plotting trajectories
addTrans <- function(color,trans){
  # This function adds transparancy to a color.
  # Define transparancy with an integer between 0 and 255
  # 0 being fully transparant and 255 being fully visable
  # Works with either color and trans a vector of equal length,
  # or one of the two of length 1.
  
  if (length(color)!=length(trans)&!any(c(length(color),length(trans))==1)) stop("Vector lengths not correct")
  if (length(color)==1 & length(trans)>1) color <- rep(color,length(trans))
  if (length(trans)==1 & length(color)>1) trans <- rep(trans,length(color))
  
  num2hex <- function(x)
  {
    hex <- unlist(strsplit("0123456789ABCDEF",split=""))
    return(paste(hex[(x-x%%16)/16+1],hex[x%%16+1],sep=""))
  }
  rgb <- rbind(col2rgb(color),trans)
  res <- paste("#",apply(apply(rgb,2,num2hex),2,paste,collapse=""),sep="")
  return(res)
}

summary.plots <- function(S, directory, odepath, plot.name, alpha, burnin=0,
                          sf.choice.save = FALSE){
  
  # sample from posterior
  samples <- sample.grab(S = S, 
                         directory = directory,
                         burnin = burnin)
  
  # process trajectories from the samples
  tp <- traj.sim(samples = samples, 
                 odepath = odepath,
                 csv = F)
  # process actual data
  dp <- data.fine.process(samples$df, loc = samples$loc)
  
  # vector of dates
  dates <- seq(as.Date("2020/1/1"), by = "day", length.out = 365)
  
  # save plots as pdf
  pdf(file = plot.name, h = 8, w = 12)
  
  ##########################################
  ### 1. population mixing parameters 
  ##########################################
  par(mfrow = c(1,1))
  
  # max value for plots
  y.max <- min(max(tp$beta.full), 10)
  
  # last day of simulated output
  end.day <- max(samples$days)
  
  # plot betas over time
  plot(dates[c(end.day - ncol(tp$beta.full) + 1):end.day],
       tp$beta.full[1,], type = "l", col = addTrans("grey", alpha),
       ylim = c(0, y.max), ylab = "Beta", xlab = "", main = "Population Mixing Parameter")
  
  for(j in 2:S){
    lines(dates[c(end.day - ncol(tp$beta.full) + 1):end.day],
          tp$beta.full[j,], 
          type = "l", 
          col = addTrans("grey", alpha))
  }
  
  # add mean, and 95% quantiles
  beta.mn <- apply(tp$beta.full, 2, mean)
  beta.uq <- apply(tp$beta.full, 2, quantile, 0.975)
  beta.lq <- apply(tp$beta.full, 2, quantile, 0.025)
  
  lines(dates[c(end.day - ncol(tp$beta.full) + 1):end.day],
        beta.mn, type = "l", col = "blue")
  
  lines(dates[c(end.day - ncol(tp$beta.full) + 1):end.day],
        beta.uq, type = "l", lty = 3, col = "blue")
  
  lines(dates[c(end.day - ncol(tp$beta.full) + 1):end.day],
        beta.lq, type = "l", lty = 3, col = "blue")
  
  
  ### age-class mixing rates
  
  if (length(which(colnames(samples$ode) == "contact-rate-80")) > 0){
    
    # contact rate parameters
    cr.inds <- c(which(colnames(samples$ode) == "contact-rate-10"),
                 which(colnames(samples$ode) == "contact-rate-20"),
                 which(colnames(samples$ode) == "contact-rate-30"),
                 which(colnames(samples$ode) == "contact-rate-40"),
                 which(colnames(samples$ode) == "contact-rate-50"),
                 which(colnames(samples$ode) == "contact-rate-60"),
                 which(colnames(samples$ode) == "contact-rate-70"),
                 which(colnames(samples$ode) == "contact-rate-80"))
      
    cr.vals <- samples$ode[,cr.inds]
    
    
    median.cr.vals <- apply(cr.vals, 2, median)
    lower.bnd.cr.vals <- apply(cr.vals, 2, quantile, 0.025)
    upper.bnd.cr.vals <- apply(cr.vals, 2, quantile, 0.975)
    
    
    if (length(which(colnames(samples$ode) == "contact-rate-postld-80"))){
      
      # lockdown days
      lock.day <- floor(samples$ode[,which(colnames(samples$ode) == "firstlockdown-endday")])
      
      # post-lockdown contact rate parameters
      post.cr.inds <- c(which(colnames(samples$ode) == "contact-rate-postld-10"),
                   which(colnames(samples$ode) == "contact-rate-postld-20"),
                   which(colnames(samples$ode) == "contact-rate-postld-30"),
                   which(colnames(samples$ode) == "contact-rate-postld-40"),
                   which(colnames(samples$ode) == "contact-rate-postld-50"),
                   which(colnames(samples$ode) == "contact-rate-postld-60"),
                   which(colnames(samples$ode) == "contact-rate-postld-70"),
                   which(colnames(samples$ode) == "contact-rate-postld-80"))
      
      post.cr.vals <- samples$ode[,post.cr.inds]
      
      median.post.cr.vals <- apply(post.cr.vals, 2, median)
      lower.bnd.post.cr.vals <- apply(post.cr.vals, 2, quantile, 0.025)
      upper.bnd.post.cr.vals <- apply(post.cr.vals, 2, quantile, 0.975)
      
      
      ages <- c("10-19", "20-29", "30-39", "40-49", 
                "50-59", "60-69", "70-79", "80+")
      
      
      ### create data frame for plotting (NEW) ###
      
      df.all <- data.frame(ages = as.factor(c(ages, ages)),
                            median.pre = c(median.cr.vals, rep(NaN, 8)),
                            lower.pre = c(lower.bnd.cr.vals, rep(NaN, 8)),
                            upper.pre = c(upper.bnd.cr.vals, rep(NaN, 8)),
                            median.post = c(rep(NaN, 8), median.post.cr.vals),
                            lower.post = c(rep(NaN, 8), lower.bnd.post.cr.vals),
                            upper.post = c(rep(NaN, 8), upper.bnd.post.cr.vals),
                            cols = as.factor(c(rep("Lockdown", 8), 
                                               rep("Post-Lockdown", 8))))
      
      df.all$cols <- relevel(df.all$cols, "Lockdown")
      
      
      ### 1. Colors ###
      
      cr.plot <- ggplot(data = df.all,
                        mapping = aes(y = ages,
                                      x = median.post,
                                      color = cols)) + 
        geom_vline(xintercept = 1, linetype = "dashed",
                   color = "black", size = 0.5) +
        ggstance::geom_linerangeh(mapping = aes(xmin = lower.pre,
                                                xmax = upper.pre),
                                  size=1,
                                  na.rm = T,
                                  position = position_nudge(y = -0.05)) +
        geom_point(mapping = aes(x = median.pre), 
                   size = 2, na.rm = T,
                   position = position_nudge(y = -0.05))+
        ggstance::geom_linerangeh(mapping = aes(xmin = lower.post,
                                                xmax = upper.post),
                                  size=1,
                                  na.rm = T,
                                  position = position_nudge(y = 0.05)) +
        geom_point(mapping = aes(x = median.post),
                   size = 2, na.rm = T,
                   position = position_nudge(y = 0.05)) +
        labs(x = "Relative Contact Rates", y = "Age Class") +
        theme_bw() + 
        scale_color_manual(name = "Contact Rates",
                           values=c("Lockdown" = "gray55", "Post-Lockdown" = "black"))
       
      
      # ### 2. Transparency ###
      # cr.plot <- ggplot(data = df.test,
      #                   mapping = aes(y = ages,
      #                                 x = median.post,
      #                                 alpha = cols)) +
      #   geom_vline(xintercept = 1, linetype = "dashed",
      #              color = "black", size = 0.5) +
      #   ggstance::geom_linerangeh(mapping = aes(xmin = lower.pre,
      #                                           xmax = upper.pre),
      #                             size=1,
      #                             na.rm = T,
      #                             position = position_nudge(y = -0.05)) +
      #   geom_point(mapping = aes(x = median.pre),
      #              size = 2, na.rm = T,
      #              position = position_nudge(y = -0.05))+
      #   ggstance::geom_linerangeh(mapping = aes(xmin = lower.post,
      #                                           xmax = upper.post),
      #                             size=1,
      #                             na.rm = T,
      #                             position = position_nudge(y = 0.05)) +
      #   geom_point(mapping = aes(x = median.post),
      #              size = 2, na.rm = T,
      #              position = position_nudge(y = 0.05)) +
      #   labs(x = "Relative Contact Rates", y = "Age Class") +
      #   theme_bw() +
      #   scale_alpha_manual(name = "Contact Rates",
      #                      values=c("Lockdown" = 0.35, "Post-Lockdown" = 1))
      # 
      # 
      # ### 3. Transparency and Colors ###
      # cr.plot <- ggplot(data = df.test,
      #                   mapping = aes(y = ages,
      #                                 x = median.post,
      #                                 alpha = cols,
      #                                 color = ages)) +
      #   geom_vline(xintercept = 1, linetype = "dashed",
      #              color = "black", size = 0.5) +
      #   ggstance::geom_linerangeh(mapping = aes(xmin = lower.pre,
      #                                           xmax = upper.pre),
      #                             size=1,
      #                             na.rm = T,
      #                             position = position_nudge(y = -0.05)) +
      #   geom_point(mapping = aes(x = median.pre),
      #              size = 2, na.rm = T,
      #              position = position_nudge(y = -0.05))+
      #   ggstance::geom_linerangeh(mapping = aes(xmin = lower.post,
      #                                           xmax = upper.post),
      #                             size=1,
      #                             na.rm = T,
      #                             position = position_nudge(y = 0.05)) +
      #   geom_point(mapping = aes(x = median.post),
      #              size = 2, na.rm = T,
      #              position = position_nudge(y = 0.05)) +
      #   labs(x = "Relative Contact Rates", y = "Age Class") +
      #   theme_bw() +
      #   scale_color_discrete(guide = "none") +
      #   scale_alpha_manual(name = "Contact Rates",
      #                      values=c("Lockdown" = 0.35, "Post-Lockdown" = 1))

      
      
      ### 4. Pre/Post Facet (OLD) ###
      
      # df.long <- data.frame(ages = as.factor(c(ages, ages)),
      #                       median.post = c(median.cr.vals,
      #                                       median.post.cr.vals),
      #                       lower.post = c(lower.bnd.cr.vals,
      #                                      lower.bnd.post.cr.vals),
      #                       upper.post = c(upper.bnd.cr.vals,
      #                                      upper.bnd.post.cr.vals),
      #                       variable = c(rep("Lockdown", 8),
      #                                    rep("Post-Lockdown", 8)))
      # 
      
      # cr.plot <- ggplot(data = df.long, 
      #                   mapping = aes(x = ages, 
      #                                 y = median.post,
      #                                 color = ages)) + 
      #   geom_hline(yintercept = 1, linetype = "dashed", 
      #              color = "black", size = 0.5) + 
      #   geom_linerange(mapping = aes(ymin = lower.post,
      #                                ymax = upper.post),
      #                  size=1) +
      #   geom_point(mapping = aes(y = median.post), size = 2) + 
      #   facet_wrap(~variable) + #, scales = "free_x") + 
      #   labs(x = "Age Class", y = "Relative Contact Rates") + 
      #   coord_flip() + theme(legend.position = "none")
      
      # print to fancy plot
      print(cr.plot)
      
      # save individual png
      ggsave(filename = paste(str_remove(plot.name, ".pdf"), 
                              "_contactrates.png", sep = ""),
             plot = cr.plot,
             width = 7, height = 4, units = "in", dpi = 300)
             
  
      # beta.medians <- matrix(nrow = ncol(tp$beta.full), ncol = 9)
      # 
      # ### age-specific beta plots
      # 
      # par(mfrow = c(3,3))
      # 
      # # titles for plot
      # titles <- c("Population Mixing Parameter, Ages 0-9",
      #             "Population Mixing Parameter, Ages 10-19",
      #             "Population Mixing Parameter, Ages 20-29",
      #             "Population Mixing Parameter, Ages 30-39",
      #             "Population Mixing Parameter, Ages 40-49",
      #             "Population Mixing Parameter, Ages 50-59",
      #             "Population Mixing Parameter, Ages 60-69",
      #             "Population Mixing Parameter, Ages 70-79",
      #             "Population Mixing Parameter, Ages 80+")
      # 
      # # titles for plot
      # legend.title <- c("Median, Ages 0-9",
      #             "Median, Ages 10-19",
      #             "Median, Ages 20-29",
      #             "Median, Ages 30-39",
      #             "Median, Ages 40-49",
      #             "Median, Ages 50-59",
      #             "Median, Ages 60-69",
      #             "Median, Ages 70-79",
      #             "Median, Ages 80+")
      # 
      # for(i in 1:9){
      #   
      #   if (i == 1){
      #     ###  baseline rate
      #     age.betas <- tp$beta.full
      #   
      #   } else {
      #     
      #     # when does post-lockdown come into effect?
      #     beta.days <- c(end.day - ncol(tp$beta.full) + 1):end.day
      #     beta.lock <- sapply(lock.day, function(x) {which(beta.days == x)})
      #     
      #     # re-calculate betas based on contact rate parameters
      #     age.betas <- matrix(data = NA_real_, 
      #                         nrow = S,
      #                         ncol = ncol(tp$beta.full))
      #     
      #     for (j in 1:S){
      #       age.betas[j,] <- c(tp$beta.full[j, 1:beta.lock[j]] * cr.vals[j, i - 1],
      #                          tp$beta.full[j, (beta.lock[j] + 1):length(beta.days)] * post.cr.vals[j, i - 1])
      #     }
      #     
      #   }
      #   
      #   beta.medians[,i] <- apply(age.betas, 2, mean)
      #   
      #   # max value for plots
      #   y.max <- max(age.betas, beta.medians[,1])
      #   
      #   # last day of simulated output
      #   end.day <- max(samples$days)
      #   
      #   # plot betas over time
      #   plot(dates[c(end.day - ncol(tp$beta.full) + 1):end.day],
      #        beta.medians[,1], type = "l", lty = 1, lwd = 4,
      #        col = "black",
      #        ylim = c(0, y.max), ylab = "Beta", xlab = "", 
      #        main = titles[i])
      # 
      #   
      #   for(j in 1:S){
      #     lines(dates[c(end.day - ncol(tp$beta.full) + 1):end.day],
      #           age.betas[j,], 
      #           type = "l", 
      #           col = addTrans("grey", alpha))
      #   }
      #   
      #   # add mean, and 95% quantiles
      #   beta.mn <- apply(age.betas, 2, mean)
      #   beta.uq <- apply(age.betas, 2, quantile, 0.975)
      #   beta.lq <- apply(age.betas, 2, quantile, 0.025)
      #   
      #   lines(dates[c(end.day - ncol(tp$beta.full) + 1):end.day],
      #         beta.mn, type = "l", col = "blue", lwd = 2)
      #   
      #   lines(dates[c(end.day - ncol(tp$beta.full) + 1):end.day],
      #         beta.uq, type = "l", lty = 3, col = "blue")
      #   
      #   lines(dates[c(end.day - ncol(tp$beta.full) + 1):end.day],
      #         beta.lq, type = "l", lty = 3, col = "blue")
      #   
      #   legend("topleft", legend=c(legend.title[i], "Baseline (Ages 0-9)"),
      #          col=c("blue", "black"), lty=c(1,1), lwd = c(2,2), cex = 0.8,
      #          box.lty=0)
      #   
      # }
      
    } else {
      
      ages <- c("10-19", "20-29", "30-39", "40-49", 
                "50-59", "60-69", "70-79", "80+")
      
      
      df.short <- data.frame(ages = as.factor(ages),
                            median.post = median.cr.vals,
                            lower.post = lower.bnd.cr.vals,
                            upper.post = upper.bnd.cr.vals)
      
      cr.plot <- ggplot(data = df.short, 
                        mapping = aes(x = ages, 
                                      y = median.post,
                                      color = ages)) + 
        geom_hline(yintercept = 1, linetype = "dashed", 
                   color = "black", size = 0.5) + 
        geom_linerange(mapping = aes(ymin = lower.post,
                                     ymax = upper.post),
                       size=1) +
        geom_point(mapping = aes(y = median.post), size = 2) + 
        labs(x = "Age Class", y = "Relative Contact Rates") + 
        coord_flip() + theme(legend.position = "none")
      
      # save to fancy-plots
      print(cr.plot)    
      
      # save individual png
      ggsave(filename = paste(str_remove(plot.name, ".pdf"), 
                              "_contactrates.png", sep = ""),
             plot = cr.plot,
             width = 7, height = 4, units = "in", dpi = 300)
  
    }
    
  }

  
  ### symptomatic cases reporting rate
  par(mfrow = c(1,1))
  
  # max value for plots
  y.max <- max(tp$rr.full)
  
  # last day of simulated output
  end.day <- max(samples$days)
  
  # plot betas over time
  plot(dates[c(end.day - ncol(tp$rr.full) + 1):end.day],
       tp$rr.full[1,], type = "l", col = addTrans("grey", alpha),
       ylim = c(0, y.max), ylab = "Reporting Rate", xlab = "", main = "Symptomatic Cases Reporting Rate")
  
  for(j in 2:S){
    lines(dates[c(end.day - ncol(tp$rr.full) + 1):end.day],
          tp$rr.full[j,], 
          type = "l", 
          col = addTrans("grey", alpha))
  }
  
  # add mean, and 95% quantiles
  rr.mn <- apply(tp$rr.full, 2, mean)
  rr.uq <- apply(tp$rr.full, 2, quantile, 0.975)
  rr.lq <- apply(tp$rr.full, 2, quantile, 0.025)
  
  lines(dates[c(end.day - ncol(tp$rr.full) + 1):end.day],
        rr.mn, type = "l", col = "blue")
  
  lines(dates[c(end.day - ncol(tp$rr.full) + 1):end.day],
        rr.uq, type = "l", lty = 3, col = "blue")
  
  lines(dates[c(end.day - ncol(tp$rr.full) + 1):end.day],
        rr.lq, type = "l", lty = 3, col = "blue")
  
  #################################################
  ### 1.5. symp-frac posterior
  #################################################
  
  if (samples$sf.choice){
  
    sf.choice.names <- c(" ", "-symp-frac-davies", "-symp-frac-equal 0.3",
                         "-symp-frac-equal 0.4", "-symp-frac-equal 0.5",
                         "-symp-frac-equal 0.6", "-symp-frac-equal 0.7")
    
    counts <- rep(0, length(sf.choice.names))
    
    for (m in 1:length(sf.choice.names)){
      counts[m] <- sum(samples$sf.vals == m)
    }
    
    probs <- counts / sum(counts)
    
    print(probs)
    
    if (sf.choice.save){
      write.csv(as.matrix(samples$sf.vals), 
              file = paste(samples$loc, "_sfchoice.csv", sep = ""),
              col.names = sf.choice.names)
    }
    
    # bar plot with symp-frac posterior probabilities
    par(mfrow = c(1,1))
    barplot(probs, main="symp-frac posterior probabilities", 
            names.arg = sf.choice.names)
    
  }
  
  #################################################
  ### 2. symptomatic cases
  #################################################
  
  # total new symptomatic cases
  par(mfrow = c(1,1))
  
  y.max <- max(max(dp$tot.sympt.new, na.rm = T),
               max(tp$tot.sympt.new))
  
  # plot total new symptomatic cases
  plot(dates[tp$days.odesim],
       tp$tot.sympt.new[1,], type = "l", col = addTrans("grey", alpha),
       ylim = c(0, y.max), ylab = "New Cases", xlab = "",
       main = "New Symptomatic Cases (Observed)")
  
  for(j in 2:S){
    lines(dates[tp$days.odesim],
          tp$tot.sympt.new[j,], 
          type = "l", 
          col = addTrans("grey", alpha))
  }
  
  
  # add mean and 95% quantiles
  mn <- apply(tp$tot.sympt.new, 2, mean)
  uq <- apply(tp$tot.sympt.new, 2, quantile, 0.975)
  lq <- apply(tp$tot.sympt.new, 2, quantile, 0.025)
  
  lines(dates[tp$days.odesim], mn, type = "l", col = "blue")
  lines(dates[tp$days.odesim], uq, type = "l", lty = 3, col = "blue")
  lines(dates[tp$days.odesim], lq, type = "l", lty = 3, col = "blue")
  
  points(dates[dp$days], dp$tot.sympt.new, pch = 19, col = "black", cex = 0.5)
  
  
  # age-strat new sympt. cases
  par(mfrow = c(3,3))
  
  if (samples$loc == "MA"){
    titles <- c("New Cases, Ages 0-19",
                "New Cases, Ages 20-29",
                "New Cases, Ages 30-39",
                "New Cases, Ages 40-49",
                "New Cases, Ages 50-59",
                "New Cases, Ages 60-69",
                "New Cases, Ages 70-79",
                "New Cases, Ages 80+")
    
    for(i in 1:8){
      y.max <- max(max(dp$sympt.new[,i], na.rm = T),
                   max(tp$sympt.new.age[[i]]))
      
      # plot total new symptomatic cases
      plot(dates[tp$days.odesim],
           tp$sympt.new.age[[i]][1,], type = "l", col = addTrans("grey", alpha),
           ylim = c(0, y.max), ylab = "New Cases", xlab = "",
           main = titles[i])
      
      for(j in 2:S){
        lines(dates[tp$days.odesim],
              tp$sympt.new.age[[i]][j,], 
              type = "l", 
              col = addTrans("grey", alpha))
      }
      
      
      # add mean and 95% quantiles
      mn <- apply(tp$sympt.new.age[[i]], 2, mean)
      uq <- apply(tp$sympt.new.age[[i]], 2, quantile, 0.975)
      lq <- apply(tp$sympt.new.age[[i]], 2, quantile, 0.025)
      
      lines(dates[tp$days.odesim], mn, type = "l", col = "blue")
      lines(dates[tp$days.odesim], uq, type = "l", lty = 3, col = "blue")
      lines(dates[tp$days.odesim], lq, type = "l", lty = 3, col = "blue")
      
      points(dates[dp$days], dp$sympt.new[,i], pch = 19, col = "black", cex = 0.5)
      
    }
  } else {
    
    titles <- c("New Cases, Ages 0-9",
                "New Cases, Ages 10-19",
                "New Cases, Ages 20-29",
                "New Cases, Ages 30-39",
                "New Cases, Ages 40-49",
                "New Cases, Ages 50-59",
                "New Cases, Ages 60-69",
                "New Cases, Ages 70-79",
                "New Cases, Ages 80+")
    
    for(i in 1:9){
      y.max <- max(max(dp$sympt.new[,i], na.rm = T),
                   max(tp$sympt.new.age[[i]]))
      
      # plot total new symptomatic cases
      plot(dates[tp$days.odesim],
           tp$sympt.new.age[[i]][1,], type = "l", col = addTrans("grey", alpha),
           ylim = c(0, y.max), ylab = "New Cases", xlab = "",
           main = titles[i])
      
      for(j in 2:S){
        lines(dates[tp$days.odesim],
              tp$sympt.new.age[[i]][j,], 
              type = "l", 
              col = addTrans("grey", alpha))
      }
      
      
      # add mean and 95% quantiles
      mn <- apply(tp$sympt.new.age[[i]], 2, mean)
      uq <- apply(tp$sympt.new.age[[i]], 2, quantile, 0.975)
      lq <- apply(tp$sympt.new.age[[i]], 2, quantile, 0.025)
      
      lines(dates[tp$days.odesim], mn, type = "l", col = "blue")
      lines(dates[tp$days.odesim], uq, type = "l", lty = 3, col = "blue")
      lines(dates[tp$days.odesim], lq, type = "l", lty = 3, col = "blue")
      
      points(dates[dp$days], dp$sympt.new[,i], pch = 19, col = "black", cex = 0.5)
    }
    
  }
  
  
  ### total cumulative symptomatic cases
  par(mfrow = c(1,1))
  
  y.max <- max(max(dp$tot.sympt.cum, na.rm = T),
               max(tp$tot.sympt.cum))
  
  # plot total new symptomatic cases
  plot(dates[tp$days.odesim],
       tp$tot.sympt.cum[1,], type = "l", col = addTrans("grey", alpha),
       ylim = c(0, y.max), ylab = "Cases", xlab = "",
       main = "Cumulative Symptomatic Cases (Observed)")
  
  for(j in 2:S){
    lines(dates[tp$days.odesim],
          tp$tot.sympt.cum[j,], 
          type = "l", 
          col = addTrans("grey", alpha))
  }
  
  
  # add mean and 95% quantiles
  mn <- apply(tp$tot.sympt.cum, 2, mean)
  uq <- apply(tp$tot.sympt.cum, 2, quantile, 0.975)
  lq <- apply(tp$tot.sympt.cum, 2, quantile, 0.025)
  
  lines(dates[tp$days.odesim], mn, type = "l", col = "blue")
  lines(dates[tp$days.odesim], uq, type = "l", lty = 3, col = "blue")
  lines(dates[tp$days.odesim], lq, type = "l", lty = 3, col = "blue")
  
  points(dates[dp$days], dp$tot.sympt.cum, pch = 19, col = "black", cex = 0.5)
  
  ### age-strat cumulative sympt. cases
  par(mfrow = c(3,3))
  
  if (samples$loc == "MA"){
    titles <- c("Cumulative Cases, Ages 0-19",
                "Cumulative Cases, Ages 20-29",
                "Cumulative Cases, Ages 30-39",
                "Cumulative Cases, Ages 40-49",
                "Cumulative Cases, Ages 50-59",
                "Cumulative Cases, Ages 60-69",
                "Cumulative Cases, Ages 70-79",
                "Cumulative Cases, Ages 80+")
    
    for(i in 1:8){
      y.max <- max(max(dp$sympt.cum[,i], na.rm = T),
                   max(tp$sympt.cum.age[[i]]))
      
      # plot total new symptomatic cases
      plot(dates[tp$days.odesim],
           tp$sympt.cum.age[[i]][1,], type = "l", col = addTrans("grey", alpha),
           ylim = c(0, y.max), ylab = "Cases", xlab = "",
           main = titles[i])
      
      for(j in 2:S){
        lines(dates[tp$days.odesim],
              tp$sympt.cum.age[[i]][j,], 
              type = "l", 
              col = addTrans("grey", alpha))
      }
      
      
      # add mean and 95% quantiles
      mn <- apply(tp$sympt.cum.age[[i]], 2, mean)
      uq <- apply(tp$sympt.cum.age[[i]], 2, quantile, 0.975)
      lq <- apply(tp$sympt.cum.age[[i]], 2, quantile, 0.025)
      
      lines(dates[tp$days.odesim], mn, type = "l", col = "blue")
      lines(dates[tp$days.odesim], uq, type = "l", lty = 3, col = "blue")
      lines(dates[tp$days.odesim], lq, type = "l", lty = 3, col = "blue")
      
      points(dates[dp$days], dp$sympt.cum[,i], pch = 19, col = "black", cex = 0.5)
    }
  } else {
    
    titles <- c("Cumulative Cases, Ages 0-9",
                "Cumulative Cases, Ages 10-19",
                "Cumulative Cases, Ages 20-29",
                "Cumulative Cases, Ages 30-39",
                "Cumulative Cases, Ages 40-49",
                "Cumulative Cases, Ages 50-59",
                "Cumulative Cases, Ages 60-69",
                "Cumulative Cases, Ages 70-79",
                "Cumulative Cases, Ages 80+")
    
    for(i in 1:9){
      y.max <- max(max(dp$sympt.cum[,i], na.rm = T),
                   max(tp$sympt.cum.age[[i]]))
      
      # plot total new symptomatic cases
      plot(dates[tp$days.odesim],
           tp$sympt.cum.age[[i]][1,], type = "l", col = addTrans("grey", alpha),
           ylim = c(0, y.max), ylab = "Cases", xlab = "",
           main = titles[i])
      
      for(j in 2:S){
        lines(dates[tp$days.odesim],
              tp$sympt.cum.age[[i]][j,], 
              type = "l", 
              col = addTrans("grey", alpha))
      }
      
      
      # add mean and 95% quantiles
      mn <- apply(tp$sympt.cum.age[[i]], 2, mean)
      uq <- apply(tp$sympt.cum.age[[i]], 2, quantile, 0.975)
      lq <- apply(tp$sympt.cum.age[[i]], 2, quantile, 0.025)
      
      lines(dates[tp$days.odesim], mn, type = "l", col = "blue")
      lines(dates[tp$days.odesim], uq, type = "l", lty = 3, col = "blue")
      lines(dates[tp$days.odesim], lq, type = "l", lty = 3, col = "blue")
      
      points(dates[dp$days], dp$sympt.cum[,i], pch = 19, col = "black", cex = 0.5)
    }
    
  }
  
  #################################################
  ### 2. hospitalization cases
  #################################################
  
  # total new hospitalizations
  par(mfrow = c(1,1))
  
  y.max <- max(max(dp$tot.hosp.new, na.rm = T),
               max(tp$tot.hosp.new))
  
  # plot total new symptomatic cases
  plot(dates[tp$days.odesim],
       tp$tot.hosp.new[1,], type = "l", col = addTrans("grey", alpha),
       ylim = c(0, y.max), ylab = "New Hospitalizations", xlab = "",
       main = "New Hospitalizations (Observed)")
  
  for(j in 2:S){
    lines(dates[tp$days.odesim],
          tp$tot.hosp.new[j,], 
          type = "l", 
          col = addTrans("grey", alpha))
  }
  
  
  # add mean and 95% quantiles
  mn <- apply(tp$tot.hosp.new, 2, mean)
  uq <- apply(tp$tot.hosp.new, 2, quantile, 0.975)
  lq <- apply(tp$tot.hosp.new, 2, quantile, 0.025)
  
  lines(dates[tp$days.odesim], mn, type = "l", col = "blue")
  lines(dates[tp$days.odesim], uq, type = "l", lty = 3, col = "blue")
  lines(dates[tp$days.odesim], lq, type = "l", lty = 3, col = "blue")
  
  points(dates[dp$days], dp$tot.hosp.new, pch = 19, col = "black", cex = 0.5)
  
  # age-strat new sympt. cases
  par(mfrow = c(3,3))
  
  if (samples$loc == "MA"){
    titles <- c("New Hospitalizations, Ages 0-19",
                "New Hospitalizations, Ages 20-29",
                "New Hospitalizations, Ages 30-39",
                "New Hospitalizations, Ages 40-49",
                "New Hospitalizations, Ages 50-59",
                "New Hospitalizations, Ages 60-69",
                "New Hospitalizations, Ages 70-79",
                "New Hospitalizations, Ages 80+")
    
    for(i in 1:8){
      y.max <- max(max(dp$hosp.new[,i], na.rm = T),
                   max(tp$hosp.new.age[[i]]))
      
      # plot total new symptomatic cases
      plot(dates[tp$days.odesim],
           tp$hosp.new.age[[i]][1,], type = "l", col = addTrans("grey", alpha),
           ylim = c(0, y.max), ylab = "New Hospitalizations", xlab = "",
           main = titles[i])
      
      for(j in 2:S){
        lines(dates[tp$days.odesim],
              tp$hosp.new.age[[i]][j,], 
              type = "l", 
              col = addTrans("grey", alpha))
      }
      
      
      # add mean and 95% quantiles
      mn <- apply(tp$hosp.new.age[[i]], 2, mean)
      uq <- apply(tp$hosp.new.age[[i]], 2, quantile, 0.975)
      lq <- apply(tp$hosp.new.age[[i]], 2, quantile, 0.025)
      
      lines(dates[tp$days.odesim], mn, type = "l", col = "blue")
      lines(dates[tp$days.odesim], uq, type = "l", lty = 3, col = "blue")
      lines(dates[tp$days.odesim], lq, type = "l", lty = 3, col = "blue")
      
      points(dates[dp$days], dp$hosp.new[,i], pch = 19, col = "black", cex = 0.5)
      
    }
  } else {
    
    titles <- c("New Hospitalizations, Ages 0-9",
                "New Hospitalizations, Ages 10-19",
                "New Hospitalizations, Ages 20-29",
                "New Hospitalizations, Ages 30-39",
                "New Hospitalizations, Ages 40-49",
                "New Hospitalizations, Ages 50-59",
                "New Hospitalizations, Ages 60-69",
                "New Hospitalizations, Ages 70-79",
                "New Hospitalizations, Ages 80+")
    
    for(i in 1:9){
      y.max <- max(max(dp$hosp.new[,i], na.rm = T),
                   max(tp$hosp.new.age[[i]]))
      
      # plot total new symptomatic cases
      plot(dates[tp$days.odesim],
           tp$hosp.new.age[[i]][1,], type = "l", col = addTrans("grey", alpha),
           ylim = c(0, y.max), ylab = "New Hospitalizations", xlab = "",
           main = titles[i])
      
      for(j in 2:S){
        lines(dates[tp$days.odesim],
              tp$hosp.new.age[[i]][j,], 
              type = "l", 
              col = addTrans("grey", alpha))
      }
      
      
      # add mean and 95% quantiles
      mn <- apply(tp$hosp.new.age[[i]], 2, mean)
      uq <- apply(tp$hosp.new.age[[i]], 2, quantile, 0.975)
      lq <- apply(tp$hosp.new.age[[i]], 2, quantile, 0.025)
      
      lines(dates[tp$days.odesim], mn, type = "l", col = "blue")
      lines(dates[tp$days.odesim], uq, type = "l", lty = 3, col = "blue")
      lines(dates[tp$days.odesim], lq, type = "l", lty = 3, col = "blue")
      
      points(dates[dp$days], dp$hosp.new[,i], pch = 19, col = "black", cex = 0.5)
    }
    
  }
  
  ### total cumulative hospitalizations
  par(mfrow = c(1,1))
  
  y.max <- max(max(dp$tot.hosp.cum, na.rm = T),
               max(tp$tot.hosp.cum))
  
  # plot total new symptomatic cases
  plot(dates[tp$days.odesim],
       tp$tot.hosp.cum[1,], type = "l", col = addTrans("grey", alpha),
       ylim = c(0, y.max), ylab = "Hospitalizations", xlab = "",
       main = "Cumulative Hospitalizations (Observed)")
  
  for(j in 2:S){
    lines(dates[tp$days.odesim],
          tp$tot.hosp.cum[j,], 
          type = "l", 
          col = addTrans("grey", alpha))
  }
  
  # add mean and 95% quantiles
  mn <- apply(tp$tot.hosp.cum, 2, mean)
  uq <- apply(tp$tot.hosp.cum, 2, quantile, 0.975)
  lq <- apply(tp$tot.hosp.cum, 2, quantile, 0.025)
  
  lines(dates[tp$days.odesim], mn, type = "l", col = "blue")
  lines(dates[tp$days.odesim], uq, type = "l", lty = 3, col = "blue")
  lines(dates[tp$days.odesim], lq, type = "l", lty = 3, col = "blue")
  
  points(dates[dp$days], dp$tot.hosp.cum, pch = 19, col = "black", cex = 0.5)
  
  
  ### age-strat cumulative hospitalizations
  par(mfrow = c(3,3))
  
  if (samples$loc == "MA"){
    titles <- c("Cumul. Hospitalizations, Ages 0-19",
                "Cumul. Hospitalizations, Ages 20-29",
                "Cumul. Hospitalizations, Ages 30-39",
                "Cumul. Hospitalizations, Ages 40-49",
                "Cumul. Hospitalizations, Ages 50-59",
                "Cumul. Hospitalizations, Ages 60-69",
                "Cumul. Hospitalizations, Ages 70-79",
                "Cumul. Hospitalizations, Ages 80+")
    
    for(i in 1:8){
      y.max <- max(max(dp$hosp.cum[,i], na.rm = T),
                   max(tp$hosp.cum.age[[i]]))
      
      # plot total new symptomatic cases
      plot(dates[tp$days.odesim],
           tp$hosp.cum.age[[i]][1,], type = "l", col = addTrans("grey", alpha),
           ylim = c(0, y.max), ylab = "Hospitalizations", xlab = "",
           main = titles[i])
      
      for(j in 2:S){
        lines(dates[tp$days.odesim],
              tp$hosp.cum.age[[i]][j,], 
              type = "l", 
              col = addTrans("grey", alpha))
      }
      
      
      # add mean and 95% quantiles
      mn <- apply(tp$hosp.cum.age[[i]], 2, mean)
      uq <- apply(tp$hosp.cum.age[[i]], 2, quantile, 0.975)
      lq <- apply(tp$hosp.cum.age[[i]], 2, quantile, 0.025)
      
      lines(dates[tp$days.odesim], mn, type = "l", col = "blue")
      lines(dates[tp$days.odesim], uq, type = "l", lty = 3, col = "blue")
      lines(dates[tp$days.odesim], lq, type = "l", lty = 3, col = "blue")
      
      points(dates[dp$days], dp$hosp.cum[,i], pch = 19, col = "black", cex = 0.5)
    }
  } else {
    
    titles <- c("Cumul. Hospitalizations, Ages 0-9",
                "Cumul. Hospitalizations, Ages 10-19",
                "Cumul. Hospitalizations, Ages 20-29",
                "Cumul. Hospitalizations, Ages 30-39",
                "Cumul. Hospitalizations, Ages 40-49",
                "Cumul. Hospitalizations, Ages 50-59",
                "Cumul. Hospitalizations, Ages 60-69",
                "Cumul. Hospitalizations, Ages 70-79",
                "Cumul. Hospitalizations, Ages 80+")
    
    for(i in 1:9){
      y.max <- max(max(dp$hosp.cum[,i], na.rm = T),
                   max(tp$hosp.cum.age[[i]]))
      
      # plot total new symptomatic cases
      plot(dates[tp$days.odesim],
           tp$hosp.cum.age[[i]][1,], type = "l", col = addTrans("grey", alpha),
           ylim = c(0, y.max), ylab = "Hospitalizations", xlab = "",
           main = titles[i])
      
      for(j in 2:S){
        lines(dates[tp$days.odesim],
              tp$hosp.cum.age[[i]][j,], 
              type = "l", 
              col = addTrans("grey", alpha))
      }
      
      
      # add mean and 95% quantiles
      mn <- apply(tp$hosp.cum.age[[i]], 2, mean)
      uq <- apply(tp$hosp.cum.age[[i]], 2, quantile, 0.975)
      lq <- apply(tp$hosp.cum.age[[i]], 2, quantile, 0.025)
      
      lines(dates[tp$days.odesim], mn, type = "l", col = "blue")
      lines(dates[tp$days.odesim], uq, type = "l", lty = 3, col = "blue")
      lines(dates[tp$days.odesim], lq, type = "l", lty = 3, col = "blue")
      
      points(dates[dp$days], dp$hosp.cum[,i], pch = 19, col = "black", cex = 0.5)
    }
  }
    
  ########################################
  ### 4. Current Hosps, Vents, and ICU ###
  ########################################
  
  
  ### current hospitalizations
  par(mfrow = c(1,1))
  
  y.max <- max(max(dp$tot.hosp.curr, na.rm = T),
               max(tp$tot.hosp.curr))
  
  # plot current hospitalizations
  plot(dates[tp$days.odesim],
       tp$tot.hosp.curr[1,], type = "l", col = addTrans("grey", alpha),
       ylim = c(0, y.max), ylab = "Hospitalizations", xlab = "",
       main = "Current Hospitalizations")
  
  for(j in 2:S){
    lines(dates[tp$days.odesim],
          tp$tot.hosp.curr[j,], 
          type = "l", 
          col = addTrans("grey", alpha))
  }
  
  
  # add mean and 95% quantiles
  mn <- apply(tp$tot.hosp.curr, 2, mean)
  uq <- apply(tp$tot.hosp.curr, 2, quantile, 0.975)
  lq <- apply(tp$tot.hosp.curr, 2, quantile, 0.025)
  
  lines(dates[tp$days.odesim], mn, type = "l", col = "blue")
  lines(dates[tp$days.odesim], uq, type = "l", lty = 3, col = "blue")
  lines(dates[tp$days.odesim], lq, type = "l", lty = 3, col = "blue")
  
  points(dates[dp$days], dp$tot.hosp.curr, pch = 19, col = "black", cex = 0.5)
  
  ### current ICU
  par(mfrow = c(1,1))
  
  y.max <- max(max(dp$tot.icu.curr, na.rm = T),
               max(tp$tot.icu.curr))
  
  # plot current ICU
  plot(dates[tp$days.odesim],
       tp$tot.icu.curr[1,], type = "l", col = addTrans("grey", alpha),
       ylim = c(0, y.max), ylab = "ICU beds", xlab = "",
       main = "Current ICU occupancy")
  
  for(j in 2:S){
    lines(dates[tp$days.odesim],
          tp$tot.icu.curr[j,], 
          type = "l", 
          col = addTrans("grey", alpha))
  }
  
  
  # add mean and 95% quantiles
  mn <- apply(tp$tot.icu.curr, 2, mean)
  uq <- apply(tp$tot.icu.curr, 2, quantile, 0.975)
  lq <- apply(tp$tot.icu.curr, 2, quantile, 0.025)
  
  lines(dates[tp$days.odesim], mn, type = "l", col = "blue")
  lines(dates[tp$days.odesim], uq, type = "l", lty = 3, col = "blue")
  lines(dates[tp$days.odesim], lq, type = "l", lty = 3, col = "blue")
  
  points(dates[dp$days], dp$tot.icu.curr, pch = 19, col = "black", cex = 0.5)
  
  ### current vent
  par(mfrow = c(1,1))
  
  y.max <- max(max(dp$tot.vent.curr, na.rm = T),
               max(tp$tot.vent.curr))
  
  # plot current vent
  plot(dates[tp$days.odesim],
       tp$tot.vent.curr[1,], type = "l", col = addTrans("grey", alpha),
       ylim = c(0, y.max), ylab = "ventilators in use", xlab = "",
       main = "Current ventilator occupancy")
  
  for(j in 2:S){
    lines(dates[tp$days.odesim],
          tp$tot.vent.curr[j,], 
          type = "l", 
          col = addTrans("grey", alpha))
  }
  
  
  # add mean and 95% quantiles
  mn <- apply(tp$tot.vent.curr, 2, mean)
  uq <- apply(tp$tot.vent.curr, 2, quantile, 0.975)
  lq <- apply(tp$tot.vent.curr, 2, quantile, 0.025)
  
  lines(dates[tp$days.odesim], mn, type = "l", col = "blue")
  lines(dates[tp$days.odesim], uq, type = "l", lty = 3, col = "blue")
  lines(dates[tp$days.odesim], lq, type = "l", lty = 3, col = "blue")
  
  points(dates[dp$days], dp$tot.vent.curr, pch = 19, col = "black", cex = 0.5)
  
  ###############################################
  ### 5. Deaths
  ###############################################
  
  ### total new deaths
  par(mfrow = c(1,1))
  
  y.max <- max(max(dp$tot.deaths.new, na.rm = T),
               max(tp$tot.deaths.new))
  
  # plot new deaths
  plot(dates[tp$days.odesim],
       tp$tot.deaths.new[1,], type = "l", col = addTrans("grey", alpha),
       ylim = c(0, y.max), ylab = "New Deaths", xlab = "",
       main = "New Deaths")
  
  for(j in 2:S){
    lines(dates[tp$days.odesim],
          tp$tot.deaths.new[j,], 
          type = "l", 
          col = addTrans("grey", alpha))
  }
  
  
  # add mean and 95% quantiles
  mn <- apply(tp$tot.deaths.new, 2, mean)
  uq <- apply(tp$tot.deaths.new, 2, quantile, 0.975)
  lq <- apply(tp$tot.deaths.new, 2, quantile, 0.025)
  
  lines(dates[tp$days.odesim], mn, type = "l", col = "blue")
  lines(dates[tp$days.odesim], uq, type = "l", lty = 3, col = "blue")
  lines(dates[tp$days.odesim], lq, type = "l", lty = 3, col = "blue")
  
  points(dates[dp$days], dp$tot.deaths.new, pch = 19, col = "black", cex = 0.5)
  
  ### age-strat new deaths
  par(mfrow = c(3,3))
  
  if (samples$loc == "MA"){
    titles <- c("New Deaths, Ages 0-19",
                "New Deaths, Ages 20-29",
                "New Deaths, Ages 30-39",
                "New Deaths, Ages 40-49",
                "New Deaths, Ages 50-59",
                "New Deaths, Ages 60-69",
                "New Deaths, Ages 70-79",
                "New Deaths, Ages 80+")
    
    for(i in 1:8){
      y.max <- max(max(dp$deaths.new[,i], na.rm = T),
                    max(tp$deaths.new.age[[i]]))
      
      # plot total new symptomatic cases
      plot(dates[tp$days.odesim],
           tp$deaths.new.age[[i]][1,], type = "l", col = addTrans("grey", alpha),
           ylim = c(0, y.max), ylab = "New Deaths", xlab = "",
           main = titles[i])
      
      for(j in 2:S){
        lines(dates[tp$days.odesim],
              tp$deaths.new.age[[i]][j,], 
              type = "l", 
              col = addTrans("grey", alpha))
      }
      
      
      # add mean and 95% quantiles
      mn <- apply(tp$deaths.new.age[[i]], 2, mean)
      uq <- apply(tp$deaths.new.age[[i]], 2, quantile, 0.975)
      lq <- apply(tp$deaths.new.age[[i]], 2, quantile, 0.025)
      
      lines(dates[tp$days.odesim], mn, type = "l", col = "blue")
      lines(dates[tp$days.odesim], uq, type = "l", lty = 3, col = "blue")
      lines(dates[tp$days.odesim], lq, type = "l", lty = 3, col = "blue")
      
      points(dates[dp$days], dp$deaths.new[,i], pch = 19, col = "black", cex = 0.5)
      
    }
  } else {
    
    titles <- c("New Deaths, Ages 0-9",
                "New Deaths, Ages 10-19",
                "New Deaths, Ages 20-29",
                "New Deaths, Ages 30-39",
                "New Deaths, Ages 40-49",
                "New Deaths, Ages 50-59",
                "New Deaths, Ages 60-69",
                "New Deaths, Ages 70-79",
                "New Deaths, Ages 80+")
    
    for(i in 1:9){
      y.max <- max(max(dp$deaths.new[,i], na.rm = T),
                    max(tp$deaths.new.age[[i]]))
      
      # plot total new symptomatic cases
      plot(dates[tp$days.odesim],
           tp$deaths.new.age[[i]][1,], type = "l", col = addTrans("grey", alpha),
           ylim = c(0, y.max), ylab = "New Deaths", xlab = "",
           main = titles[i])
      
      for(j in 2:S){
        lines(dates[tp$days.odesim],
              tp$deaths.new.age[[i]][j,], 
              type = "l", 
              col = addTrans("grey", alpha))
      }
      
      
      # add mean and 95% quantiles
      mn <- apply(tp$deaths.new.age[[i]], 2, mean)
      uq <- apply(tp$deaths.new.age[[i]], 2, quantile, 0.975)
      lq <- apply(tp$deaths.new.age[[i]], 2, quantile, 0.025)
      
      lines(dates[tp$days.odesim], mn, type = "l", col = "blue")
      lines(dates[tp$days.odesim], uq, type = "l", lty = 3, col = "blue")
      lines(dates[tp$days.odesim], lq, type = "l", lty = 3, col = "blue")
      
      points(dates[dp$days], dp$deaths.new[,i], pch = 19, col = "black", cex = 0.5)
      
    }
  }  
  
  ### total cumulative deaths
  par(mfrow = c(1,1))
  
  y.max <- max(max(dp$tot.deaths.cum, na.rm = T),
               max(tp$tot.deaths.cum))
  
  # plot new deaths
  plot(dates[tp$days.odesim],
       tp$tot.deaths.cum[1,], type = "l", col = addTrans("grey", alpha),
       ylim = c(0, y.max), ylab = "Deaths", xlab = "",
       main = "Cumulative Deaths")
  
  for(j in 2:S){
    lines(dates[tp$days.odesim],
          tp$tot.deaths.cum[j,], 
          type = "l", 
          col = addTrans("grey", alpha))
  }
  
  
  # add mean and 95% quantiles
  mn <- apply(tp$tot.deaths.cum, 2, mean)
  uq <- apply(tp$tot.deaths.cum, 2, quantile, 0.975)
  lq <- apply(tp$tot.deaths.cum, 2, quantile, 0.025)
  
  lines(dates[tp$days.odesim], mn, type = "l", col = "blue")
  lines(dates[tp$days.odesim], uq, type = "l", lty = 3, col = "blue")
  lines(dates[tp$days.odesim], lq, type = "l", lty = 3, col = "blue")
  
  points(dates[dp$days], dp$tot.deaths.cum, pch = 19, col = "black", cex = 0.5)
  
  ### age-strat cum. deaths
  par(mfrow = c(3,3))
  
  if (samples$loc == "MA"){
    titles <- c("Cumulative Deaths, Ages 0-19",
                "Cumulative Deaths, Ages 20-29",
                "Cumulative Deaths, Ages 30-39",
                "Cumulative Deaths, Ages 40-49",
                "Cumulative Deaths, Ages 50-59",
                "Cumulative Deaths, Ages 60-69",
                "Cumulative Deaths, Ages 70-79",
                "Cumulative Deaths, Ages 80+")
    
    for(i in 1:8){
      y.max <- max(max(dp$deaths.cum[,i], na.rm = T),
                  max(tp$deaths.cum.age[[i]]))
      
      # plot total new symptomatic cases
      plot(dates[tp$days.odesim],
           tp$deaths.cum.age[[i]][1,], type = "l", col = addTrans("grey", alpha),
           ylim = c(0, y.max), ylab = "Deaths", xlab = "",
           main = titles[i])
      
      for(j in 2:S){
        lines(dates[tp$days.odesim],
              tp$deaths.cum.age[[i]][j,], 
              type = "l", 
              col = addTrans("grey", alpha))
      }
      
      
      # add mean and 95% quantiles
      mn <- apply(tp$deaths.cum.age[[i]], 2, mean)
      uq <- apply(tp$deaths.cum.age[[i]], 2, quantile, 0.975)
      lq <- apply(tp$deaths.cum.age[[i]], 2, quantile, 0.025)
      
      lines(dates[tp$days.odesim], mn, type = "l", col = "blue")
      lines(dates[tp$days.odesim], uq, type = "l", lty = 3, col = "blue")
      lines(dates[tp$days.odesim], lq, type = "l", lty = 3, col = "blue")
      
      points(dates[dp$days], dp$deaths.cum[,i], pch = 19, col = "black", cex = 0.5)
    }
  } else {
    
    titles <- c("Cumulative Deaths, Ages 0-9",
                "Cumulative Deaths, Ages 10-19",
                "Cumulative Deaths, Ages 20-29",
                "Cumulative Deaths, Ages 30-39",
                "Cumulative Deaths, Ages 40-49",
                "Cumulative Deaths, Ages 50-59",
                "Cumulative Deaths, Ages 60-69",
                "Cumulative Deaths, Ages 70-79",
                "Cumulative Deaths, Ages 80+")
    
    for(i in 1:9){
      y.max <- max(max(dp$deaths.cum[,i], na.rm = T),
                   max(tp$deaths.cum.age[[i]]))
      
      # plot total new symptomatic cases
      plot(dates[tp$days.odesim],
           tp$deaths.cum.age[[i]][1,], type = "l", col = addTrans("grey", alpha),
           ylim = c(0, y.max), ylab = "Deaths", xlab = "",
           main = titles[i])
      
      for(j in 2:S){
        lines(dates[tp$days.odesim],
              tp$deaths.cum.age[[i]][j,], 
              type = "l", 
              col = addTrans("grey", alpha))
      }
      
      
      # add mean and 95% quantiles
      mn <- apply(tp$deaths.cum.age[[i]], 2, mean)
      uq <- apply(tp$deaths.cum.age[[i]], 2, quantile, 0.975)
      lq <- apply(tp$deaths.cum.age[[i]], 2, quantile, 0.025)
      
      lines(dates[tp$days.odesim], mn, type = "l", col = "blue")
      lines(dates[tp$days.odesim], uq, type = "l", lty = 3, col = "blue")
      lines(dates[tp$days.odesim], lq, type = "l", lty = 3, col = "blue")
      
      points(dates[dp$days], dp$deaths.cum[,i], pch = 19, col = "black", cex = 0.5)
    }
  }  
  
  
  ### home and hosp deaths
  par(mfrow = c(1,1))
  
  if (samples$loc == "MA"){
    y.max <- max(tp$tot.home.deaths.new)
  } else {
    # new home deaths
    y.max <- max(max(dp$tot.home.deaths.new, na.rm = T),
                 max(tp$tot.home.deaths.new))
  }
  
  # plot new deaths
  plot(dates[tp$days.odesim],
       tp$tot.home.deaths.new[1,], type = "l", col = addTrans("grey", alpha),
       ylim = c(0, y.max), ylab = "Deaths", xlab = "",
       main = "New Home Deaths")
  
  for(j in 2:S){
    lines(dates[tp$days.odesim],
          tp$tot.home.deaths.new[j,], 
          type = "l", 
          col = addTrans("grey", alpha))
  }
  
  
  # add mean and 95% quantiles
  mn <- apply(tp$tot.home.deaths.new, 2, mean)
  uq <- apply(tp$tot.home.deaths.new, 2, quantile, 0.975)
  lq <- apply(tp$tot.home.deaths.new, 2, quantile, 0.025)
  
  lines(dates[tp$days.odesim], mn, type = "l", col = "blue")
  lines(dates[tp$days.odesim], uq, type = "l", lty = 3, col = "blue")
  lines(dates[tp$days.odesim], lq, type = "l", lty = 3, col = "blue")
  
  points(dates[dp$days], dp$tot.home.deaths.new, pch = 19, col = "black", cex = 0.5)
  
  
  # new hosp deaths
  par(mfrow = c(1,1))
  
  if (samples$loc == "MA"){
    y.max <- max(tp$tot.hosp.deaths.new)
  } else {
    y.max <- max(max(dp$tot.hosp.deaths.new, na.rm = T),
               max(tp$tot.hosp.deaths.new))
  }
  
  # plot new deaths
  plot(dates[tp$days.odesim],
       tp$tot.hosp.deaths.new[1,], type = "l", col = addTrans("grey", alpha),
       ylim = c(0, y.max), ylab = "Deaths", xlab = "",
       main = "New Hospital Deaths")
  
  for(j in 2:S){
    lines(dates[tp$days.odesim],
          tp$tot.hosp.deaths.new[j,], 
          type = "l", 
          col = addTrans("grey", alpha))
  }
  
  
  # add mean and 95% quantiles
  mn <- apply(tp$tot.hosp.deaths.new, 2, mean)
  uq <- apply(tp$tot.hosp.deaths.new, 2, quantile, 0.975)
  lq <- apply(tp$tot.hosp.deaths.new, 2, quantile, 0.025)
  
  lines(dates[tp$days.odesim], mn, type = "l", col = "blue")
  lines(dates[tp$days.odesim], uq, type = "l", lty = 3, col = "blue")
  lines(dates[tp$days.odesim], lq, type = "l", lty = 3, col = "blue")
  
  points(dates[dp$days], dp$tot.hosp.deaths.new, pch = 19, col = "black", cex = 0.5)
  
  
  ### cumulative home and hosp deaths
  par(mfrow = c(1,1))
  
  # cumul home deaths
  
  if (samples$loc == "MA"){
    y.max <- max(tp$tot.home.deaths.cum)
  } else {
    y.max <- max(max(dp$tot.home.deaths.cum, na.rm = T),
               max(tp$tot.home.deaths.cum))
  }
  
  # plot new deaths
  plot(dates[tp$days.odesim],
       tp$tot.home.deaths.cum[1,], type = "l", col = addTrans("grey", alpha),
       ylim = c(0, y.max), ylab = "Deaths", xlab = "",
       main = "Cumulative Home Deaths")
  
  for(j in 2:S){
    lines(dates[tp$days.odesim],
          tp$tot.home.deaths.cum[j,], 
          type = "l", 
          col = addTrans("grey", alpha))
  }
  
  
  # add mean and 95% quantiles
  mn <- apply(tp$tot.home.deaths.cum, 2, mean)
  uq <- apply(tp$tot.home.deaths.cum, 2, quantile, 0.975)
  lq <- apply(tp$tot.home.deaths.cum, 2, quantile, 0.025)
  
  lines(dates[tp$days.odesim], mn, type = "l", col = "blue")
  lines(dates[tp$days.odesim], uq, type = "l", lty = 3, col = "blue")
  lines(dates[tp$days.odesim], lq, type = "l", lty = 3, col = "blue")
  
  if (samples$loc != "MA"){
    points(dates[dp$days], dp$tot.home.deaths.cum, pch = 19, col = "black", cex = 0.5)
  }
  
  # cumul hosp deaths
  par(mfrow = c(1,1))
  
  if (samples$loc == "MA"){
    y.max <- max(tp$tot.hosp.deaths.cum)
  } else {
    y.max <- max(max(dp$tot.hosp.deaths.cum, na.rm = T),
               max(tp$tot.hosp.deaths.cum))
  }
  
  # plot new deaths
  plot(dates[tp$days.odesim],
       tp$tot.hosp.deaths.cum[1,], type = "l", col = addTrans("grey", alpha),
       ylim = c(0, y.max), ylab = "Deaths", xlab = "",
       main = "Cumulative Hospital Deaths")
  
  for(j in 2:S){
    lines(dates[tp$days.odesim],
          tp$tot.hosp.deaths.cum[j,], 
          type = "l", 
          col = addTrans("grey", alpha))
  }
  
  
  # add mean and 95% quantiles
  mn <- apply(tp$tot.hosp.deaths.cum, 2, mean)
  uq <- apply(tp$tot.hosp.deaths.cum, 2, quantile, 0.975)
  lq <- apply(tp$tot.hosp.deaths.cum, 2, quantile, 0.025)
  
  lines(dates[tp$days.odesim], mn, type = "l", col = "blue")
  lines(dates[tp$days.odesim], uq, type = "l", lty = 3, col = "blue")
  lines(dates[tp$days.odesim], lq, type = "l", lty = 3, col = "blue")
  
  points(dates[dp$days], dp$tot.hosp.deaths.cum, pch = 19, col = "black", cex = 0.5)
  
  
  ###############################################
  ### 6. Discharges
  ###############################################
  
  
  ### new discharges
  par(mfrow = c(1,1))
  
  # new home deaths
  if (samples$loc == "MA"){
    y.max <- max(tp$tot.hosp.discharges.new)
  } else {
    y.max <- max(max(dp$tot.hosp.discharges.new, na.rm = T),
               max(tp$tot.hosp.discharges.new))
  }
  
  # plot new deaths
  plot(dates[tp$days.odesim],
       tp$tot.hosp.discharges.new[1,], type = "l", col = addTrans("grey", alpha),
       ylim = c(0, y.max), ylab = "New Discharges", xlab = "",
       main = "New Hospital Discharges")
  
  for(j in 2:S){
    lines(dates[tp$days.odesim],
          tp$tot.hosp.discharges.new[j,], 
          type = "l", 
          col = addTrans("grey", alpha))
  }
  
  
  # add mean and 95% quantiles
  mn <- apply(tp$tot.hosp.discharges.new, 2, mean)
  uq <- apply(tp$tot.hosp.discharges.new, 2, quantile, 0.975)
  lq <- apply(tp$tot.hosp.discharges.new, 2, quantile, 0.025)
  
  lines(dates[tp$days.odesim], mn, type = "l", col = "blue")
  lines(dates[tp$days.odesim], uq, type = "l", lty = 3, col = "blue")
  lines(dates[tp$days.odesim], lq, type = "l", lty = 3, col = "blue")
  
  if (samples$loc != "MA"){
    points(dates[dp$days], dp$tot.hosp.discharges.new, pch = 19, col = "black", cex = 0.5)
  }
  
  
  ### cum discharges
  par(mfrow = c(1,1))
  
  if (samples$loc == "MA"){
    y.max <- max(tp$tot.hosp.discharges.cum)
  } else {
    y.max <- max(max(dp$tot.hosp.discharges.cum, na.rm = T),
               max(tp$tot.hosp.discharges.cum))
  }
  # plot cum discharges
  plot(dates[tp$days.odesim],
       tp$tot.hosp.discharges.cum[1,], type = "l", col = addTrans("grey", alpha),
       ylim = c(0, y.max), ylab = "Discharges", xlab = "",
       main = "Cumulative Hospital Discharges")
  
  for(j in 2:S){
    lines(dates[tp$days.odesim],
          tp$tot.hosp.discharges.cum[j,], 
          type = "l", 
          col = addTrans("grey", alpha))
  }
  
  
  # add mean and 95% quantiles
  mn <- apply(tp$tot.hosp.discharges.cum, 2, mean)
  uq <- apply(tp$tot.hosp.discharges.cum, 2, quantile, 0.975)
  lq <- apply(tp$tot.hosp.discharges.cum, 2, quantile, 0.025)
  
  lines(dates[tp$days.odesim], mn, type = "l", col = "blue")
  lines(dates[tp$days.odesim], uq, type = "l", lty = 3, col = "blue")
  lines(dates[tp$days.odesim], lq, type = "l", lty = 3, col = "blue")
  
  if (samples$loc != "MA"){
    points(dates[dp$days], dp$tot.hosp.discharges.cum, pch = 19, col = "black", cex = 0.5)
  }
  
  dev.off()
}




