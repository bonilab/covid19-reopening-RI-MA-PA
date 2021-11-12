#!/usr/bin/env Rscript

### multichain-mcmc.R
### last edited: 12 Nov 2021
### authors: Ephraim Hanks, Nathan Wikle

# This function sets up n.chains independent MCMC chains in parallel. It 
#   requires access to multiple processors (either through the R 'doMC' 
#   parallel backend or a Linux PSOCK cluster). MCMC output is periodically saved
#   in Rdata files.
multichain.mcmc.odesim <- function(
  n.chains,               # number of MCMC chains to run in parallel
  n.cores,                # number of processors (should == n.chains)
  n.iter,                 # number of iterations (N = n.iter * n.mcmc.per.iter)
  n.mcmc.per.iter,        # number of mcmc samples per iteration
  parallel.type = "doMC", # type of parallelization (either "doMC" or "psock")
  inf.file.dir = "./",    # path to source file directory
  save.file.name = "out", # name of output files
  betas.lst,              # initial beta parameter values
  ode.params.lst,         # initial ode parameter values
  lik.params.lst,         # initial neg. bin. dispersion values
  rr.params.lst,          # initial reporting rate parameters
  s2.hosp.lst,            # initial hospitalization variance parameter
  s2.vent.lst,            # initial ventilator variance parameter
  s2.icu.lst,             # initial ICU variance parameter
  Sigma.tune.start,       # MCMC proposal distribution covariance matrix
  var.tune.start,         # MCMC proposal distribution variance scaling param
  adapt = "ShabyWells",   # type of MCMC adaptive proposal scheume
  thin.rt = 1,            # MCMC thinning rate
  ...
){

  # set up proposal distribution tuning parameters    
  Sigma.t <- Sigma.tune.start
  var.t <- list()
  for (i in 1:n.chains) {
    var.t[[i]] <- var.tune.start
  }
  
  # adaptive proposal parameter
  t.adpt <- 0
    
  # parallelize, either using doMC or Foreach
  if (parallel.type == "doMC") {
    # parallelize using doMC
    require(doMC)
    registerDoMC(cores = n.cores)
  } else if (parallel.type == "psock") {

    # parallelization using PSOCK

    # set up Foreach (with number of processors = n.chains)
    mp_type <- "PSOCK"
    cl <- parallel::makeCluster(n.chains, type = mp_type)
    doParallel::registerDoParallel(cl)

    # set up source code and requisite libraries on each processor
    clusterCall(cl, function() {
      # library for vectorized multinomial
      library(mc2d)
      # library for multivariate normal distribution
      library(mvtnorm)
      # library for spline-based expansions
      library(fda)
      # library for truncated normal
      library(msm)
      # functions for local use
      source(paste(inf.file.dir, "loglik-odesim.R", sep = ""), chdir = TRUE)
      source(paste(inf.file.dir, "mcmc-odesim.R", sep = ""))
      source(paste(inf.file.dir, "traj-from-params.R", sep = ""))
      source(paste(inf.file.dir, "traj-process.R", sep = ""))
      source(paste(inf.file.dir, "data-process.R", sep = ""))
      source(paste(inf.file.dir, "plots-odesim.R", sep = ""))
    })
  }
  
  # run mcmc sampler on each processor, N = n.iter * n.mcmc.per.iter iterations
  for(iter in 1:n.iter){
    cat(iter, "\n")
    # run mcmc sampler on each processor, n = n.mcmc.per.iter
    out <- foreach(i = 1:n.chains) %dopar% mcmc.odesim(
      n.mcmc = n.mcmc.per.iter,
      beta.start = betas.lst[[i]],
      ode.params.start = ode.params.lst[[i]],
      report.rate.params.start = rr.params.lst[[i]],
      lik.params.start = lik.params.lst[[i]],
      s2.hosp.start = s2.hosp.lst[[i]],
      s2.icu.start = s2.icu.lst[[i]],
      s2.vent.start = s2.vent.lst[[i]],
      Sigma.tune = Sigma.t,
      var.tune = var.t[[i]],
      adapt.type = adapt,
      thin = thin.rt,
      t.adapt.start = t.adpt,
      ...
    )
        
    # determine which samples to keep for analysis 
    resamp.idx <- 1:n.chains
    
    # add output from each cluster into a list
    last.idx <- n.mcmc.per.iter / thin.rt
    for (k in 1:n.chains) {
      betas.lst[[k]] <- out[[resamp.idx[k]]]$beta[last.idx, ]
      ode.params.lst[[k]] <- out[[resamp.idx[k]]]$ode.params[last.idx, ]
      lik.params.lst[[k]] <- out[[resamp.idx[k]]]$lik.params[last.idx, ]
      rr.params.lst[[k]] <- out[[resamp.idx[k]]]$rr.params[last.idx, ]
      s2.hosp.lst[[k]] <- out[[resamp.idx[k]]]$s2.params[last.idx, 1]
      s2.icu.lst[[k]] <- out[[resamp.idx[k]]]$s2.params[last.idx, 2]
      s2.vent.lst[[k]] <- out[[resamp.idx[k]]]$s2.params[last.idx, 3]
      var.t[[k]] <- out[[resamp.idx[k]]]$var.tune
    }

    # update adaptive proposal parameter
    t.adpt <- out[[1]]$t.adapt.end
    
    # save output as Rdata file
    iter.char <- as.character(iter)
    num.leading.zeros <- nchar(as.character(n.iter)) - nchar(iter.char)
    if (num.leading.zeros == 0) {
      leading.zeros <- ""
    } else if (num.leading.zeros == 1) {
      leading.zeros <- "0"
    } else if (num.leading.zeros == 2) {
      leading.zeros <- "00"
    } else if (num.leading.zeros == 3) {
      leading.zeros <- "000"
    }
    save(out, file = paste(save.file.name, "-", leading.zeros,
      as.character(iter), ".Rdata",
      sep = ""
    ))
  }

  # return output
  return(out)
}



