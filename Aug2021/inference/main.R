#!/usr/bin/env Rscript

### main.R
### last edited: 16 Nov 2021
### authors: Nathan Wikle

### This file contains instructions and source code to reproduce the results
###   from Wikle, Tran, et al's (2021) manuscript, "SARS-CoV-2 epidemic after 
###   social and economic reopening in three US states reveals shifts in age 
###   structure and clinical characteristics."
###
### Software requirements are included in the README associated with this
###   repository. IMPORTANT: The working directory must be set to the location of
###   this file (main.R).


#################################################################################
### 1. Load packages  (version number)
#################################################################################
library(mc2d)           # 0.1.18
library(mvtnorm)        # 1.1.1
library(readxl)         # 1.3.1
library(fda)            # 5.1.9
library(snow)           # 0.4.3
library(doParallel)     # 1.0.16
library(foreach)        # 1.5.1
library(msm)            # 1.6.8
library(splines2)       # 1.18
library(doMC)           # 1.3.7
library(extrafont)      # 0.17
library(scales)         # 1.1.1
library(dplyr)          # 1.0.7
library(forcats)        # 0.5.1
library(viridis)        # 0.6.1
library(hrbrthemes)     # 0.8.0
library(ggplot2)        # 3.3.5
library(ggridges)       # 0.5.2


#################################################################################
### 2. Make ODESIM (deterministic solve of mechanistic model of SARS-CoV-2 epidemic)
#################################################################################

# path to ODESIM
odepath <- "../cpp-v5-discharges-nonhospdeaths/"
# compile ODESIM
system(paste("make -C", odepath))


#################################################################################
### 3. Inference
#################################################################################

### The following R scripts recreate the inference procedures used to obtain all
###   results in the manuscript and supplementary materials. This includes the
###   "best fitting models" shown in Figures 2, 5, and S13 and S14, as well as
###   the model comparison discussed in Table 3, the sensitivity analyses in
###   Supp. Materials (SM) Section 5.1 (time from symptoms to hospitalization),
###   SM Section 5.2 (time from symptoms to presentation), and SM Section 6
###   (excess deaths). The output from these scripts are stored in a directory,
###   "./output/MCMC-redo", that is automatically created when the file is souced.
###
### Note that output from these files is already included in 
###   "./output/manuscript-results"; it is not necessary that these files be 
###   sourced in order if interests is only in the existing results. 
###   
### If you choose to recreate all inference results, it is HIGHLY RECOMMENDED
###   that these files be run in PARALLEL. A single R script, distributed across 
###   five 2.2 GHz Intel Xeon processors from the PSU high-performance computing 
###   infrastructure (PSU Roar Supercomputer), takes ~2.5 days to run,.
###   If run sequentially, this section of the analysis would take ~65 days.

### computing set-up

# use python
OPT_USE_PY_LL = TRUE # set to FALSE if not using Python multinomial function
# number of processors
n.proc <- 5
# number of chains (defaults to n.proc)
n.chains <- n.proc
# type of parallelization
parallelization <- "doMC" # "psock" # must be either "psock" or "doMC"

# number of mcmc samples before saving
n.mcmc <- 1000
# number of times files are saved
n.iter <- 500
# note: total number of iterations: n.mcmc * n.iter

# burnin
#   note: by default, mcmc is thinned to every 10 samples. so,
#         burnin <- 100 indicates that the first 1000 samples are discarded.
n.burnin <- 30000
# remove .Rdata files
remove.rdata <- TRUE

### A. Rhode Island

# ## i. Model comparisons (Table 3)
# source("./mcmc-setup/RI/RI-c-icu-d2.R") # BEST FIT! Output used to make Fig 2
# source("./mcmc-setup/RI/RI-c-noicu-d2.R")
# source("./mcmc-setup/RI/RI-noc-icu-d2.R")
# source("./mcmc-setup/RI/RI-noc-noicu-d2.R")

# ## ii. Time-from-symptoms-to-hospitalization (Supp. Materials, Section 5.1)
# source("./mcmc-setup/RI/RI-sth-2.R")
# source("./mcmc-setup/RI/RI-sth-3.R")
# source("./mcmc-setup/RI/RI-sth-4.R")
# source("./mcmc-setup/RI/RI-sth-5.R")
# source("./mcmc-setup/RI/RI-sth-6.R")
# source("./mcmc-setup/RI/RI-sth-7.R")

# ## iii. Time-from-symptoms-to-presentation (Supp. Materials, Section 5.2)
# source("./mcmc-setup/RI/RI-stp-d0.R")
# source("./mcmc-setup/RI/RI-stp-d1.R")
# source("./mcmc-setup/RI/RI-stp-d2.R")
# source("./mcmc-setup/RI/RI-stp-d3.R")
# source("./mcmc-setup/RI/RI-stp-d4.R")

### B. Massachusetts

# ## i. Model comparisons (Table 3)
# source("./mcmc-setup/MA/MA-c-icu-d2.R") # BEST FIT!
# source("./mcmc-setup/MA/MA-c-noicu-d2.R")
# source("./mcmc-setup/MA/MA-noc-icu-d2.R")
# source("./mcmc-setup/MA/MA-noc-noicu-d2.R")

# ## ii. Excess deaths (Supp. Materials, Section 6)
# source("./mcmc-setup/MA/MA-ex-deaths.R")

### C. Pennsylvania

# ## i. Model comparisons (Table 3)
# source("./mcmc-setup/PA/PA-c-icu-d2.R") # BEST FIT!
# source("./mcmc-setup/PA/PA-c-noicu-d2.R")
# source("./mcmc-setup/PA/PA-noc-icu-d2.R")
# source("./mcmc-setup/PA/PA-noc-noicu-d2.R")

# ## ii. Excess deaths (Supp. Materials, Section 6)
# source("./mcmc-setup/PA/PA-ex-deaths.R")


#################################################################################
### 4. Figures
#################################################################################

### Create Figures 2 and 5 from the manuscript, as well as figures S3-S17 in the
###   supplementary materials. 
###
### Note: This script takes ~30 minutes to complete.
###
source("./create-figures.R")

