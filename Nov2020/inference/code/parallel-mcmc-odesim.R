#!/usr/bin/env Rscript

# parallel-mcmc-odesim.R
# authors: Ephraim Hanks, Nathan Wikle, Emily Strong
# last edited: 20 Nov 2020 
#
# 


multichain.mcmc.odesim <- function(
  n.chains,                 # number of mcmc chains to run in parallel
  n.cores,                  # number of cores to request (should == n.chains)
  n.iter,                   # number of Rdata objects (with n.mcmc.per.iter samples each) to run
  n.mcmc.per.iter,          # number of mcmc samples per iteration
  resample = FALSE,         # whether resampling should be used (default = FALSE)
  parallel.type = "doMC",   # type of parallelization ("doMC" (default) or "psock")
  inf.file.dir = "./",      # location of inference files (eg, mcmc-odesim.R, loglik.odesim.R, etc)
  save.file.name = "out",   # prefix used when saving Rdata output 
  betas.lst,                # list of initial beta values, length = n.chains
  ode.params.lst,           # list of initial ode.param starting values, length = n.chains
  lik.params.lst,           # list of nb dispersion param starting values, length = n.chains
  rr.params.lst,            # list of reporrting rate loading staring values, length = n.chains
  s2.hosp.lst,              # list of curr hosp var starting values, length = n.chains
  s2.vent.lst,              # list of curr vent var starting values, length = n.chains
  s2.icu.lst,               # list of curr icu var starting values, length = n.chains
  Sigma.tune.start,         # covariance matrix for proposal dist.
  var.tune.start,           # variance for proposal dist.
  adapt = "ShabyWells",     # type of adaptive proposal scheme used (default = ShabyWells)
  thin.rt = 1,              # rate to thin saved posterior samples (default = 1)
  ...                       # additional arguments to mcmc.odesim (see documentation in mcmc-odesim.R for details)
){
    
  # set up proposal density structures
  Sigma.t <- Sigma.tune.start
  var.t <- list()
  for(i in 1:n.chains){
    var.t[[i]] <- var.tune.start
  }
  t.adpt <- 0
    
  if(parallel.type=="doMC"){
    ## start parallel
    require(doMC)
    registerDoMC(cores=n.cores)
   }
  
  if(parallel.type=="psock"){
    # Parallelization using Foreach
        
    # set up Foreach (with number of processors = n.chs!)
    mp_type = "PSOCK"
    cl <- parallel::makeCluster(n.chains, type = mp_type)
    doParallel::registerDoParallel(cl)
        
    clusterCall(cl, function() {
      ## library for vectorized multinomial
      library(mc2d)
      ## library for multivariate normal distribution
      library(mvtnorm)
      ## library for spline-based expansions
      library(fda)
      ## library for truncated normal
      library(msm)
      ## functions for local use
      source(paste(inf.file.dir, "loglik-odesim.R", sep = ""), chdir = TRUE)
      source(paste(inf.file.dir, "mcmc-odesim.R", sep = ""))
      source(paste(inf.file.dir, "traj-from-params.R", sep = ""))
      source(paste(inf.file.dir, "traj-process.R", sep = ""))
      source(paste(inf.file.dir, "data-process.R", sep = ""))
      source(paste(inf.file.dir, "plots-odesim.R", sep = ""))
    })    
  }
    
  for(iter in 1:n.iter){
    cat(iter, "\n")
    out <- foreach(i=1:n.chains) %dopar% mcmc.odesim(n.mcmc = n.mcmc.per.iter,
                                                         beta.start=betas.lst[[i]],
                                                         ode.params.start=ode.params.lst[[i]],
                                                         report.rate.params.start=rr.params.lst[[i]],
                                                         lik.params.start=lik.params.lst[[i]],
                                                         s2.hosp.start=s2.hosp.lst[[i]],
                                                         s2.icu.start=s2.icu.lst[[i]],
                                                         s2.vent.start=s2.vent.lst[[i]],
                                                         Sigma.tune=Sigma.t,
                                                         var.tune=var.t[[i]],
                                                         adapt.type= adapt,
                                                         thin=thin.rt,
                                                         t.adapt.start=t.adpt,
                                                         ...)
        
        
    ## resampling
    if (resample) {

      ll.vec=out[[1]]$loglik.final.iter
      for(k in 2:n.chains){
        ll.vec <- c(ll.vec, out[[2]]$loglik.final.iter)
      }
    
      resample.probs <- exp(ll.vec - max(ll.vec))
      resamp.idx <- sample(1:n.chains,n.chains, replace = TRUE, prob = resample.probs)
    
    } else {
      resamp.idx <- 1:n.chains
    }

    # save iteration outpuf from each chain
    last.idx <- n.mcmc.per.iter / thin.rt
    for(k in 1:n.chains){
      betas.lst[[k]] <- out[[resamp.idx[k]]]$beta[last.idx, ]
      ode.params.lst[[k]] <- out[[resamp.idx[k]]]$ode.params[last.idx, ]
      lik.params.lst[[k]] <- out[[resamp.idx[k]]]$lik.params[last.idx, ]
      rr.params.lst[[k]] <- out[[resamp.idx[k]]]$rr.params[last.idx, ]
      s2.hosp.lst[[k]] <- out[[resamp.idx[k]]]$s2.params[last.idx, 1]
      s2.icu.lst[[k]] <- out[[resamp.idx[k]]]$s2.params[last.idx, 2]
      s2.vent.lst[[k]] <- out[[resamp.idx[k]]]$s2.params[last.idx, 3]
      var.t[[k]] <- out[[resamp.idx[k]]]$var.tune
    }
    
    t.adpt <- out[[1]]$t.adapt.end
    
    # save output from all chains in .Rdata object
    iter.char <- as.character(iter)
    num.leading.zeros <- nchar(as.character(n.iter)) - nchar(iter.char)
    
    # determine number of leading zeros to include when naming file
    if(num.leading.zeros == 0){
      leading.zeros <- ""
    }
    
    if(num.leading.zeros == 1){
      leading.zeros <- "0"
    }
    
    if(num.leading.zeros == 2){
      leading.zeros <- "00"
    }
    
    if(num.leading.zeros == 3){
      leading.zeros <- "000"
    }
    
    # save output from most recent iteration
    save(out, file = paste(save.file.name, "-", leading.zeros, as.character(iter), ".Rdata", sep=""))
  }

  # return output
  return(out)
}

