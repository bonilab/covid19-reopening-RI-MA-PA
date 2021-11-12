#!/usr/bin/env Rscript

### traj-from-params.R
### last edited: 12 Nov 2021
### authors: Ephraim Hanks, Nathan Wikle

### This wrapper function takes R objects and converts them into valid inputs and 
###   function calls to the ODESIM program.
traj.from.params <- function(
  beta,                     # vector of beta values for ODESIM
  tf,                       # integer day to stop sim (1 = 20200101)
  introday = NULL,          # day of first infected (typically = 55) 
  odepath,                  # path to ODESIM location
  loc = "RI",               # state used for the analysis
  params = NULL,            # vector of inferred ODESIM parameters
  const.params = NULL,      # vector of constant ODESIM parameters
  non.odesim.params = NULL, # parameters not used in ODESIM  
  symp = NULL               # symp-frac specification
){

  # remove non.odesim.params from params and const.params
  if (length(non.odesim.params) > 0) {
    idx.remove <- integer()

    for (k in 1:length(non.odesim.params)) {
      idx.remove <- c(idx.remove, which(names(params) == non.odesim.params[k]))
    }

    if (length(idx.remove) > 0) {
      params <- params[-idx.remove]
    } else {
      params <- NULL
    }
  }
    
  # make introday backwards compatible
  if (length(which(c(names(params), names(const.params)) == "introday")) + length(introday) == 0) {
    introday <- 55
  } else if (length(which(c(names(params), names(const.params)) == "introday")) + length(introday) > 1) {
    cat("WARNING: introday specified in two or more places.", "\n")
  }
    
  # create string to run odesim with parameters
  cmd1 <- paste(c(
    odepath, "odesim none -binary-output -tf ",
    as.character(tf)
  ), collapse = "")

  if (length(introday) == 1) {
    cmd1.5 <- paste(c(" -introday"), as.character(introday), collapse = " ")
    cmd2 <- paste(c(" -beta", as.character(beta)), collapse = " ")
    cmd <- paste(cmd1, cmd1.5, cmd2, sep = "")
  } else {
    cmd2 <- paste(c(" -beta", as.character(beta)), collapse = " ")
    cmd <- paste(cmd1, cmd2, sep = "")
  }
  
  # add constants, if any
  if (length(const.params) > 0) {
    const.names <- names(const.params)
    for (k in 1:length(const.params)) {
      cmd <- paste(c(
        cmd, " -", const.names[k], " ",
        as.character(const.params[k])
      ), collapse = "")
    }
  }
  
  # add parameters
  if (length(params) > 0) {
    names = names(params)
    for (k in 1:length(params)) {
      cmd = paste(c(
        cmd, " -", names(params)[k], " ",
        as.character(params[k])
      ), collapse = "")
    }
  }
  
  # specify symp-frac (used when sf.choice == TRUE)
  if (length(symp) > 0) {
    cmd <- paste(cmd, symp, collapse = "")
  }
  
  # add location parameter
  cmd <- paste(cmd, "-loc ", as.character(loc), "2>&1", collapse = "")
  
  # number of expected columns in output
  ncols <- 307
  # upper limit on number of output lines
  nlines <- as.integer(550) 

  # call the odesim program, with beta as input
  p <- pipe(cmd, 'rb')
  bin_result <- readBin(p, "double", nlines * ncols)
  close.connection(p)
  
  # reshape 1D result to (*-by-ncols)
  dim(bin_result) <- c(ncols, length(bin_result) / ncols)
  bin_result <- aperm(bin_result, c(2,1))

  # return trajectory results
  return(bin_result)
}
