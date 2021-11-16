#!/usr/bin/env Rscript

### process-results.R
### last edited: 14 Nov 2021
### authors: Ephraim Hanks, Nathan Wikle


### Process large MCMC output (.Rdata files) into small CSVs containing
###   1000 samples from the posterior. In addition, a README file with relevant
###   metadata for the analysis is also produced.
process.results <- function(
  out.folder,                   # location of .Rdata files
  burnin = 1,                   # burnin size
  loc = "RI",                   # location
  odepath = "../cpp-v5-discharges-nonhospdeaths",  # ODESIM location
  odesim.version = "v5",        # ODESIM version
  name = "output",              # output file name
  readme = TRUE,                # create README file
  cleanFiles = TRUE,            # remove .Rdata files after saving CSV (Boolean) 
  non.odesim.params = NULL,     # vector of non-odesim parameters
  const.params = NULL,          # vector of constant parameters
  lik.hosp.discharges = FALSE,  # hospital discharge data (Boolean)
  pres.delay = 2,               # delay from symptoms to presentation 
  active.surv = FALSE           # use active surveillance data (Boolean)
){
  
  ###############################################################################
  ### 1. Combine output from all .Rdata files
  ###############################################################################

  # grab file names
  files.list <- list.files(out.folder, "*.Rdata")
  files.list <- files.list[order(nchar(files.list), files.list)]
  max.iter <- length(files.list)

  # load all .Rdata files into the environment
  load(paste(out.folder, files.list[1], sep = ""))
  
  # data structures  
  data <- out
  n.chains <- length(data)
  n.mcmc <- nrow(data[[1]]$beta)
  n.beta <- ncol(data[[1]]$beta)
  n.ode.params <- ncol(data[[1]]$ode.params)
  n.rr.params <- ncol(data[[1]]$rr.params)
  n.lik.params <- ncol(data[[1]]$lik.params)
  n.s2.params <- ncol(data[[1]]$s2.params)
  spline <- data[[1]]$spline
  df <- data[[1]]$df
  beta.chains <- array(NA, dim = c(n.mcmc * max.iter, n.beta, n.chains))
  ode.chains <- array(NA, dim = c(n.mcmc * max.iter, n.ode.params, n.chains))
  rr.chains <- array(NA, dim = c(n.mcmc * max.iter, n.rr.params, n.chains))
  lik.chains <- array(NA, dim = c(n.mcmc * max.iter, n.lik.params, n.chains))
  s2.chains <- array(NA, dim = c(n.mcmc * max.iter, n.s2.params, n.chains))
  loglik.chains <- matrix(NA, nrow = n.mcmc * max.iter, ncol = n.chains)

  # save results in data objects
  for (iter in 1:max.iter) {
    ## load in data and call it "data"
    load(paste(out.folder, files.list[iter], sep = ""))
    data <- out
    ## get betas and ode.params
    for (i in 1:n.chains) {
      beta.chains[(iter - 1) * n.mcmc + (1:n.mcmc), , i] = data[[i]]$beta
      ode.chains[(iter - 1) * n.mcmc + (1:n.mcmc), , i] = data[[i]]$ode.params
      rr.chains[(iter - 1) * n.mcmc + (1:n.mcmc), , i] = data[[i]]$rr.params
      lik.chains[(iter - 1) * n.mcmc + (1:n.mcmc), , i] = data[[i]]$lik.params
      s2.chains[(iter - 1) * n.mcmc + (1:n.mcmc), , i] = data[[i]]$s2.params
      loglik.chains[(iter - 1) * n.mcmc + (1:n.mcmc), i] = data[[i]]$loglik
    }
  }


  df <- data[[1]]$df
  days <- df$daynum
  Z.beta <- eval.basis(61:max(days), data[[1]]$spline.beta)
  
  if (is.matrix(data[[1]]$spline.rr)) {
    Z.rr <- data[[1]]$spline.rr
    Z.rr <- Z.rr[-nrow(Z.rr), ]
  } else {
    Z.rr <- eval.basis(61:max(days), data[[1]]$spline.rr)
  }


  # read in data and save trajectories after burnin
  dp <- data.process(df, loc = loc)
  traj.list <- list()
  idx.remove <- ceiling(burnin / n.mcmc)
  loglik.save <- matrix(NA, nrow = max.iter - idx.remove + 1, ncol = n.chains)
  for (i in (idx.remove:max.iter)) {
    load(paste(out.folder, files.list[i], sep = ""))
    for (k in 1:n.chains) {
      traj.list[[length(traj.list) + 1]] <- traj.process(
        out[[k]]$traj,
        loc = loc,
        odesim.ver = odesim.version
      )
      loglik.save[i - idx.remove + 1, k] <- out[[k]]$loglik.final.iter
    }
  }
  n.traj <- length(traj.list)
  

  ###############################################################################
  ### 2. Make CSV output files
  ###############################################################################
  
  # number of samples to save
  M <- 1000
  idx <- round(seq(burnin + 1, dim(beta.chains)[1], length.out = M / n.chains))

  ### save betas 
  beta.idx <- beta.chains[idx, , ]
  beta.plot <- beta.idx[, , 1]
  
  for (i in 2:n.chains) {
    beta.plot <- rbind(beta.plot, beta.idx[, , i])
  }
  
  # each row of daily.beta.idx is a daily beta value from day = 61 to 
  #   the last day of data
  daily.beta <- t(Z.beta %*% t(beta.plot))
  
  write.table(rbind(61:(60 + ncol(daily.beta)), daily.beta),
    file = paste(out.folder, name, "_daily-betas_day", max(days), ".csv", sep = ""), 
    sep = ",", row.names = FALSE, col.names = FALSE
  )
  
  
  ### save ode parameters
  params.idx <- ode.chains[idx, , 1]
  for (i in 2:n.chains) {
    params.idx <- rbind(params.idx, ode.chains[idx, , i])
  }
  
  rr.idx <- rr.chains[idx, , 1]
  if (dim(rr.chains)[2] == 1) {
    rr.idx <- matrix(rr.idx, ncol = 1, nrow = length(rr.idx))
    for (i in 2:n.chains) {
      rr.idx <- rbind(rr.idx, matrix(rr.chains[idx, , i], ncol = 1))
    }
  } else {
    for (i in 2:n.chains) {
      rr.idx <- rbind(rr.idx, rr.chains[idx, , i])
    }
  }
  
  colnames(rr.idx) <- paste("rr.", as.character(1:ncol(rr.idx)), sep = "")
  colnames(params.idx) <- colnames(data[[1]]$ode.params)
  colnames(daily.beta) <- 61:(60 + ncol(daily.beta))
    
  write.table(cbind(params.idx, rr.idx),
    file = paste(out.folder, name, "_ode-params_day", max(days), ".csv", sep = ""), sep = ",",
    row.names = FALSE, col.names = TRUE
  )
  
  
  ###############################################################################
  ### 3. DIC calculations
  ###############################################################################

  # get posterior mean values
  sample.index <- cbind(sample.int(
    n = nrow(beta.chains[-c(1:burnin), , ]),
    size = 2500
  ))
  
  # mean betas
  small.beta <- beta.chains[burnin + sample.index, , ]
  full.beta <- apply(small.beta, c(1, 3), function(X) {
    Z.beta %*% X
  })
  daily.beta.hat <- rowMeans(apply(full.beta, 2, rowMeans))

  # mean reporting rate
  small.rr <- rr.chains[burnin + sample.index, , ]
  
  if (dim(rr.chains)[2] > 1) {
    full.rr <- apply(small.rr, c(1, 3), function(X) {
      Z.rr %*% X
    })
    daily.rr.hat <- rowMeans(apply(full.rr, 2, rowMeans))
  } else {
    daily.rr.hat <- Z.rr * mean(small.rr)
  }
  
  # odesim parameters
  ode.params.hat <- apply(ode.chains[burnin + sample.index, , ], 2, mean)
  names(ode.params.hat) <- colnames(data[[1]]$ode.params)

  # likelihood variance
  lik.params.hat <- apply(lik.chains[burnin + sample.index, , ], 2, mean)
  
  # other params
  s2.params.hat <- apply(s2.chains[burnin + sample.index, , ], 2, mean)
  
  # get loglikelihood of posterior mean
  if (length(which(names(data[[1]]) == "introday")) == 1) {
    introday <- data[[1]]$introday
  } else {
    introday <- 60
  }
  
  # get non.odesim.params (backwards compatibility)
  if (length(which(names(data[[1]]) == "non.odesim.params")) > 0) {
    n.o.p <- data[[1]]$non.odesim.params
  } else {
    n.o.p <- non.odesim.params
  }
  
  # get const.params (backwards compatibility)
  if (length(which(names(data[[1]]) == "const.params")) > 0) {
    c.p <- data[[1]]$const.params
  } else {
    c.p <- const.params
  }
  
  # end.day
  end.day <- max(data[[1]]$df$daynum, na.rm = T) + 1

  # trajectory using posterior mean
  traj.hat <- traj.from.params(
    beta = as.numeric(daily.beta.hat),
    params = ode.params.hat,
    tf = end.day,
    introday = introday,
    const.params = c.p,
    non.odesim.params = n.o.p,
    odepath = odepath,
    loc = loc,
    symp = NULL
  )
  
  # extra parameters (like hospitalization reporting rate )
  extra.params <- NULL
  extra.const.params <- NULL
  extra.params.fitted.idx <- integer()
  extra.params.const.idx <- integer()
  
  if (length(non.odesim.params) > 0) {
    for (k in 1:length(non.odesim.params)) {
      extra.params.fitted.idx <- c(
        extra.params.fitted.idx,
        which(names(ode.params.hat) == non.odesim.params[k])
      )
    }

    if (length(extra.params.fitted.idx) > 0) {
      extra.params <- ode.params.hat[extra.params.fitted.idx]
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
  

  # create P matrix
  P <- matrix(0, nrow = end.day + pres.delay + 1, ncol = end.day + 1)
  for (j in 1:ncol(P)) {
    P[j:(j + pres.delay), j] <- c(rep(0, pres.delay), 1)
  }

  # check that hospital discharge information is included
  if (is.null(data[[1]]$lik.hosp.discharges)) {
    data[[1]]$lik.hosp.discharges = lik.hosp.discharges
  }
  
  
  if (is.null(data[[1]]$active.surv)) {
    data[[1]]$active.surv = active.surv
  }
  
  # calculate loglikelihood for posterior mean estimate
  loglik.hat <- loglik.odesim(
    traj = traj.hat,
    df = df,
    dp = dp,
    odesim.ver = odesim.version,
    P = P,
    loc = loc,
    report.rate = c(daily.rr.hat, daily.rr.hat[length(daily.rr.hat)]),
    nb.disp.params = lik.params.hat,
    lik.tot = data[[1]]$lik.tot,
    lik.age = data[[1]]$lik.age,
    lik.hosp.new = data[[1]]$lik.hosp.new,
    lik.hosp.curr = data[[1]]$lik.hosp.curr,
    lik.icu.curr = data[[1]]$lik.icu.curr,
    lik.vent.curr = data[[1]]$lik.vent.curr,
    lik.tot.deaths = data[[1]]$lik.tot.deaths,
    lik.home.deaths = data[[1]]$lik.home.deaths,
    lik.hosp.deaths = data[[1]]$lik.hosp.deaths,
    lik.age.deaths = data[[1]]$lik.age.deaths,
    lik.hosp.discharges = data[[1]]$lik.hosp.discharges,
    lik.curr = data[[1]]$lik.curr,
    lik.old = data[[1]]$lik.old,
    active.surv = data[[1]]$active.surv,
    p.asympt = data[[1]]$p.asympt,
    total.size.constraint = FALSE,
    s2.hosp = s2.params.hat[1],
    s2.icu = s2.params.hat[2],
    s2.vent = s2.params.hat[3],
    extra.params = extra.params, 
    extra.const.params = extra.const.params
  )

  # calculate DIC
  Dhat <- -2 * loglik.hat$ll
  Dbar <- -2 * mean(loglik.chains[burnin + sample.index, ], na.rm = TRUE)
  pD <- Dbar - Dhat
  DIC <- pD + Dbar
  
  # compare with DIC using Gelman's penalization term:
  pG = mean(apply(loglik.chains[burnin + sample.index, ], 2, var))
  DIC.gelman <- Dhat + 2 * pG
  
  ###############################################################################
  ### 4. Make README file
  ###############################################################################
    
  if (readme) {

    # add useful metadata for future reference

    mcmc1 <- out[[1]]
    
    tab <- data.frame(input = integer(), value = integer())
    f <- 0
    tab[f + 1, 1] <- "Location"
    tab[f + 1, 2] <- mcmc1$loc
    tab[f + 2, 1] <- "Last Day of Data Used"
    tab[f + 2, 2] <- max(mcmc1$df$daynum)
    tab[f + 3, 1] <- "Day of MCMC Run"
    tab[f + 3, 2] <- as.character(mcmc1$today)
    tab[f + 4, ] <- " "
    
    j <- nrow(tab)
    tab[j + 1, ] <- " "
    tab[j + 2, 1] <- "Number of MCMC Chains"
    tab[j + 2, 2] <- length(out)
    tab[j + 3, 1] <- "Number of Iterations Per Chain"
    tab[j + 3, 2] <- n.mcmc * max.iter * mcmc1$thin
    tab[j + 4, 1] <- "Burnin (values before this discarded)"
    tab[j + 4, 2] <- burnin
    tab[j + 5, 1] <- "Adaptive Tuning Type"
    tab[j + 5, 2] <- mcmc1$adapt.type
    
    g <- nrow(tab)
    tab[g + 1, 1] <- "lik.tot"
    tab[g + 1, 2] <- mcmc1$lik.tot
    tab[g + 2, 1] <- "lik.age"
    tab[g + 2, 2] <- mcmc1$lik.age
    tab[g + 3, 1] <- "lik.hosp.new"
    tab[g + 3, 2] <- mcmc1$lik.hosp.new
    tab[g + 4, 1] <- "lik.hosp.curr"
    tab[g + 4, 2] <- mcmc1$lik.hosp.curr
    tab[g + 5, 1] <- "lik.icu.curr"
    tab[g + 5, 2] <- mcmc1$lik.icu.curr
    tab[g + 6, 1] <- "lik.vent.curr"
    tab[g + 6, 2] <- mcmc1$lik.vent.curr
    tab[g + 7, 1] <- "lik.tot.deaths"
    tab[g + 7, 2] <- mcmc1$lik.tot.deaths
    tab[g + 8, 1] <- "lik.home.deaths"
    tab[g + 8, 2] <- mcmc1$lik.home.deaths
    tab[g + 9, 1] <- "lik.age.deaths"
    tab[g + 9, 2] <- mcmc1$lik.age.deaths
    tab[g + 10, 1] <- "lik.hosp.deaths"
    tab[g + 10, 2] <- mcmc1$lik.hosp.deaths
    tab[g + 11, 1] <- "lik.hosp.discharges"
    tab[g + 11, 2] <- mcmc1$lik.hosp.discharges
    tab[g + 12, 1] <- "lik.curr"
    tab[g + 12, 2] <- mcmc1$lik.curr
    tab[g + 13, 1] <- "lik.old"
    tab[g + 13, 2] <- mcmc1$lik.old
    
    tab[g + 14, ] <- " "
    tab[g + 15, 1] <- "odesim version"
    tab[g + 15, 2] <- mcmc1$odesim.ver
    tab[g + 16, ] <- " "
    tab[g + 17, 1] <- "Mean Log-Likelihood"
    tab[g + 17, 2] <- mean(as.numeric(loglik.save))
    tab[g + 18, 1] <- "SD of Log-Likelihood"
    tab[g + 18, 2] <- sqrt(var(as.numeric(loglik.save)))
    tab[g + 19, ] <- " "
    tab[g + 20, 1] <- "Dbar"
    tab[g + 20, 2] <- Dbar
    tab[g + 21, 1] <- "Dhat"
    tab[g + 21, 2] <- Dhat
    tab[g + 22, 1] <- "DIC"
    tab[g + 22, 2] <- DIC
    tab[g + 23, 1] <- "DIC (Gelman)"
    tab[g + 23, 2] <- DIC.gelman
    tab[g + 24, ] <- " "
    
    h <- nrow(tab) - 1
    tab[h + 2, ] <- " "
    tab[h + 3, 1] <- "Contact Rate Betas (spline length)"
    tab[h + 3, 2] <- mcmc1$spline.beta$nbasis
    tab[h + 4, 1] <- "Contact Rate Betas (min day of spline)"
    tab[h + 4, 2] <- mcmc1$spline.beta$rangeval[1]
    tab[h + 5, 1] <- "Contact Rate Betas (max day of spline)"
    tab[h + 5, 2] <- mcmc1$spline.beta$rangeval[2]
    tab[h + 6, ] <- " "
    tab[h + 7, 1] <- "Reporting Rate (spline length)"
    
    if (is.matrix(mcmc1$spline.rr)) {
      tab[h + 7, 2] <- ncol(mcmc1$spline.rr)
      tab[h + 8, 1] <- "Reporting Rate (min day of spline)"
      tab[h + 8, 2] <- 61
      tab[h + 9, 1] <- "Reporting Rate (max day of spline)"
      tab[h + 9, 2] <- max(mcmc1$df$daynum, na.rm = T) + 1
      tab[h + 10, ] <- " "
    } else {
      tab[h + 7, 2] <- mcmc1$spline.rr$nbasis
      tab[h + 8, 1] <- "Reporting Rate (min day of spline)"
      tab[h + 8, 2] <- mcmc1$spline.rr$rangeval[1]
      tab[h + 9, 1] <- "Reporting Rate (max day of spline)"
      tab[h + 9, 2] <- mcmc1$spline.rr$rangeval[2]
      tab[h + 10, ] <- " "
    }
    
    tab$max <- " "
    i <- nrow(tab)
    tab[i + 1, ] <- " "
    tab[i + 2, 1] <- "Included ODESIM Params"
    tab[i + 2, 2] <- "Prior Min"
    tab[i + 2, 3] <- "Prior Max"
    n.params <- ncol(mcmc1$ode.params)
    tab[i + 2 + (1:n.params), ] <- cbind(
      colnames(mcmc1$ode.params),
      mcmc1$ode.params.prior.min, 
      mcmc1$ode.params.prior.max
    )

    if (length(mcmc1$const.params) > 0) {
      tab[i + 2 + n.params + 1, ] <- " "
      tab[i + 2 + n.params + 2, ] <- c("Constant ODESIM Params", "value", " ")
      n.const.params <- length(mcmc1$const.params)
      tab[i + 2 + n.params + 2 + (1:n.const.params), ] <- cbind(
        names(mcmc1$const.params), 
        mcmc1$const.params, 
        rep(" ", n.const.params)
      )
      if (length(mcmc1$introday) > 0) {
        tab[i + 2 + n.params + 2 + n.const.params + 1, ] <- c("introday", mcmc1$introday, " ")
      }
    } else {
      if (length(mcmc1$introday) > 0) {
        tab[i + 2 + n.params + 1, ] <- " "
        tab[i + 2 + n.params + 2, ] <- c("Constant ODESIM Params", "value", " ")
        tab[i + 2 + n.params + 3, ] <- c("introday", mcmc1$introday, " ")
      }
    }
    
    # write out README file
    write.table(tab,
      file = paste(out.folder, name, "-README", ".txt", sep = ""),
      sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE
    ) 
  }


  ###############################################################################
  ### 5. Remove .Rdata files
  ###############################################################################

  if (cleanFiles) {
    for (j in 1:length(files.list)) {
      # remove .Rdata files
      system(paste("rm ", out.folder, files.list[j], sep = ""))
    }
  }
}



