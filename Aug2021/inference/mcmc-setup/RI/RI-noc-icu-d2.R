#!/usr/bin/env Rscript

### RI-noc-icu-d2.R
### last edited: 14 Nov 2021
### Ephraim Hanks, Nathan Wikle

#################################################################################
### 1. Preliminaries
#################################################################################

### create directory for output files

# file name
file.name <- "RI-noc-icu-d2"

# create directory
dir.create("../output/MCMC-redo")
dir.create("../output/MCMC-redo/RI")
dir.create(paste("../output/MCMC-redo/RI/", file.name, sep = ""))

### functions
source("./data-process.R")
source("./traj-process.R")
source("./traj-from-params.R")
source("./loglik-odesim.R", chdir = TRUE)
source("./mcmc-odesim.R")
source("./multichain-mcmc.R")
source("./plots-odesim.R")
source("./process-results.R")

### ODESIM location
odepath <- "../cpp-v5-discharges-nonhospdeaths/"

### load and clean data
ri.data <- read.csv("../data/RI-data_2021Sep06.csv")
ri.data <- data.frame(ri.data)

# removing days before 61
idx.remove <- which(ri.data$daynum < 61)
if (length(idx.remove) > 0) {
  ri.data <- ri.data[-idx.remove, ]
}

# number of days
n.days <- nrow(ri.data)
# vector of days
days <- ri.data$daynum

# end day (don't change!)
end.day <- max(days, na.rm = TRUE) + 1
# number of days with data (starts at 61 [01 Mar 2021])
num.days <- end.day - 60 

# process data
dp <- data.process(ri.data)

### contact rate (beta) basis functions 

# cubic spline with one basis function every 7 days
bspl <- create.bspline.basis(c(61, end.day), nbasis = round(num.days / 7))

# create basis evaluation function matrices
Z <- eval.basis(bspl, 61:end.day)

### reporting rate basis functions

# create I-spline basis functions, evaluated every day from 61:end.day
knots <- c(84, 92, 100, 108, 130, 160, 190)
Z.rr <- iSpline(x = 61:end.day, knots = knots, degree = 2)
Z.rr <- cbind(rep(1, nrow(Z.rr)), Z.rr)
Z.rr <- Z.rr[, c(1, 4:7)]

### create testing delay probabilities (obsolete)

# list of 4 vectors, each corresponding to a reporting period:
#   p1: Beginning - March 13,
#   p2: March 14 - March 23,
#   p3: March 24 - March 27,
#   p4: March 28 - present
# the 1st term = prob of 0 day delay, 2nd = prob. of 1 day delay, etc.
delay.probs <- list(
  p1 = c(1, rep(0, 7)),
  p2 = c(1, rep(0, 7)),
  p3 = c(1, rep(0, 7)),
  p4 = c(1, rep(0, 7))
)


#################################################################################
### 2. Initial values for parameters
#################################################################################

# starting values for reporting rate
rr.strt <- c(0.2, 0.1, 0.5, 0.1, 0.0005) 
rr.daily <- Z.rr %*% rr.strt

# starting values for beta (length should == ncol(Z))
beta.strt <- c(
  2, 2, 1, 0.5, 0.5, 0.35, 0.25, 0.25, 0.15, 0.15, 0.1, 0.1, 0.1, 0.1, 0.1,
  0.1, 0.1, 0.1, 0.1, 0.1, 0.2, .1, .2, .3, .3, .2, .1
)

# starting values for ODESIM and other parameters
prms.start <- c(
  9.12, 0.03, 0.15, 0.3,
  0.7, 0.7, 0.7, 150, 0.7, 0.6,
  0.027, 0.03, 0.044, 0.085, 0.12, 0.2, 0.29, 0.18,
  1, 1, 1, 1, 1, 1, 1, 1
)

names(prms.start) <- c(
  "mean-time-vent", "death-prob-home-60", "death-prob-home-70", "death-prob-home-80",
  "dev-len-hospstay", "dev-icu-frac", "dev-icu-frac-phase2", "dev-icu-frac-phase2beginday",
  "prob-icu-vent", "dev-ventdeath-mid",
  "hosp-frac-10", "hosp-frac-20", "hosp-frac-30", "hosp-frac-40", "hosp-frac-50",
  "hosp-frac-60", "hosp-frac-70", "hosp-frac-80",
  "contact-rate-10", "contact-rate-20", "contact-rate-30", "contact-rate-40",
  "contact-rate-50", "contact-rate-60", "contact-rate-70", "contact-rate-80"
)

# prior lower bounds
params.prior.min <- c(
  6, 0, 0, 0.1,
  0.25, 0.25, 0.25, 100, 0.5, 0.4,
  .001, .001, .02, .05, .08, .14, .15, .15,
  .1, .1, .1, .1, .1, .1, .1, .1
)

# prior upper bounds
params.prior.max <- c(
  12, 0.15, 0.25, 0.4,
  1.5, 1.0, 1.0, 200, 1.0, 1.3,
  .04, .04, .07, .12, .25, .3, .35, .3,
  10, 10, 10, 10, 10, 10, 10, 10
)

# names of parameters that are NOT odesim parameters, but are still to be estimated
non.odesim.params <- NULL

# inputs to odesim that are to be FIXED, and not estimated
const.prms.start <- c("", 2, 4)
names(const.prms.start) <- c("symp-frac-davies", "steps-per-day", "time-symp-to-hosp")

# rough calculation for tuning variance - make proposal variance smaller for 
#   parameters with smaller prior range
range <- params.prior.max - params.prior.min
ode.params.var.start <- range^2 / (mean(range^2))
names(ode.params.var.start) <- names(prms.start)

# starting values for NB dispersion parameters
nb.disp.start <- c(1, 1, 0.8, 0.01, 1.1) 

# create a "list" of starting values for all chains
betas.0 <- list()
ode.params.0 <- list()
const.params.0 <- list()
const.params.0 <- const.prms.start
rr.params.0 <- list()
lik.params.0 <- list()
s2.hosp.0 <- list()
s2.vent.0 <- list()
s2.icu.0 <- list()
for (k in 1:n.chains) {
  betas.0[[k]] <- beta.strt
  ode.params.0[[k]] <- prms.start
  rr.params.0[[k]] <- rr.strt
  lik.params.0[[k]] <- nb.disp.start
  s2.hosp.0[[k]] <- 1
  s2.vent.0[[k]] <- 1
  s2.icu.0[[k]] <- 1
}

# initialize proposal variance
var.tune.0 = c(.0001, .001, .007, 25, .01)

Sigma.tune.0 = list(
  Sigma.tune.beta = diag(length(betas.0[[1]])),
  Sigma.tune.ode.params = diag(ode.params.var.start),
  Sigma.tune.rr.params = diag(length(rr.params.0[[1]])),
  Sigma.tune.s2 = diag(3),
  Sigma.tune.ll.params = diag(length(lik.params.0[[1]]))
)


#################################################################################
### 3. MCMC chains (parallel)
#################################################################################

# generate posterior samples with MCMC
fit <- multichain.mcmc.odesim(
  # parallel set-up #
  parallel.type = parallelization, 
  n.chains = n.chains, 
  n.cores = n.proc, 
  n.iter = n.iter, 
  n.mcmc.per.iter = n.mcmc, 
  save.file.name = paste("../output/MCMC-redo/RI/", file.name, "/out", sep = ""),
  inf.file.dir = "./",
  # MCMC parameters #
  df = ri.data,
  loc = "RI",
  odepath = odepath, 
  odesim.ver = "v5",
  start.day = 61,
  end.day = max(ri.data$daynum, na.rm = TRUE) + 1,  
  introday = 55, 
  thin.rt = 10, 
  plot.save = FALSE, 
  # Likelihood specification #
  lik.tot = FALSE, 
  lik.age = TRUE, 
  lik.hosp.new = TRUE, 
  lik.tot.deaths = FALSE, 
  lik.age.deaths = TRUE, 
  lik.home.deaths = TRUE,
  lik.hosp.deaths = TRUE,
  lik.hosp.discharges = TRUE, 
  lik.hosp.curr = TRUE, 
  lik.icu.curr = TRUE, 
  lik.vent.curr = TRUE, 
  lik.curr = TRUE,
  lik.old = FALSE,
  p.vecs = delay.probs, 
  pres.delay = 2,
  active.surv = FALSE,
  p.asympt = 0.4, 
  total.size.constraint = TRUE, 
  sf.choice = FALSE,
  # Initialize parameters #
  betas.lst = betas.0, 
  spline.beta = bspl,
  rr.params.lst = rr.params.0, 
  spline.rr = Z.rr, 
  s2.hosp.lst = s2.hosp.0, 
  s2.vent.lst = s2.vent.0, 
  s2.icu.lst = s2.icu.0, 
  const.params = const.params.0,
  ode.params.lst = ode.params.0,
  ode.params.prior.min = params.prior.min, 
  ode.params.prior.max = params.prior.max,
  non.odesim.params = non.odesim.params,
  lik.params.lst = lik.params.0, 
  fixed.nb.disp = FALSE,  
  # Adaptive proposal set-up #
  adapt.iter = 100, 
  adapt = "ShabyWells", 
  var.tune.start = var.tune.0,
  Sigma.tune.start = Sigma.tune.0  
)


#################################################################################
### 4. Convert output to CSV
#################################################################################

# Process output Rdata files, saving CSV files and a summary README.
process.results(
  out.folder = paste("../output/MCMC-redo/RI/", file.name, "/", sep = ""),
  burnin = n.burnin,                  
  loc = "RI",                   
  odepath = odepath,  
  odesim.version = "v5",        
  name = file.name,             
  readme = TRUE,                
  cleanFiles = remove.rdata,           
  non.odesim.params = non.odesim.params,  
  const.params = const.prms.start,
  lik.hosp.discharges = TRUE,  
  pres.delay = 2,               
  active.surv = FALSE          
)

