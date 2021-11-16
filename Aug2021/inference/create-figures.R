#!/usr/bin/env Rscript

### create-figures.R
### last edited: 15 Nov 2021
### authors: Nathan Wikle
###
### This script creates Figures 2 and 5 from the manuscript, as well as 
###   S3-S17 in the Supplementary Materials.


#################################################################################
### 1. Source functions
#################################################################################
source("./traj-from-params.R")
source("./loglik-odesim.R", chdir = TRUE)
source("./plots-odesim.R")
source("./process-results.R")
source("./ms-figures.R")
source("./supp-figures.R")

#################################################################################
### 2. load Ubuntu Mono fonts
#################################################################################
loadfonts()

#################################################################################
### 3. Manuscript (Figures 2 and 5)
#################################################################################

### Figure 2 (RI posterior trajectories)

set.seed(94)

const.vec.ri <- c("", 2, 4)
names(const.vec.ri) <- c("symp-frac-davies", "steps-per-day", "time-symp-to-hosp")

fig2panel(
  beta.directory = "../output/manuscript-results/RI/RI-bestfit_daily-betas_day250.csv",
  ode.directory = "../output/manuscript-results/RI/RI-bestfit_ode-params_day250.csv",
  data.directory = "../data/RI-data_2021Sep06.csv",
  odepath = "../cpp-v5-discharges-nonhospdeaths/",
  loc = "RI", const = const.vec.ri,
  plot.name = "../output/figures/manuscript/ms-figure2.pdf", 
  alpha = 50,
  subsample = 250,
  axis.size = 1.75, lab.size = 1.75,
  ncol = 2, nrow = 4, height = 12, width = 12,
  font.family = "Ubuntu Mono", delay = 2
)

### Figure 5 (posterior comparisons of reporting rate and ODE params)

const.vec.ri <- c("", 2, 4)
names(const.vec.ri) <- c("symp-frac-davies", "steps-per-day", "time-symp-to-hosp")

const.vec.ma <- c(2, "", 4)
names(const.vec.ma) = c("steps-per-day", "symp-frac-davies", "time-symp-to-hosp")

const.vec.pa <- c("", 0.8, 2, 4)
names(const.vec.pa) <- c("symp-frac-davies", "prob-icu-vent", "steps-per-day", "time-symp-to-hosp")

fig5panel(
  ri.beta = "../output/manuscript-results/RI/RI-bestfit_daily-betas_day250.csv",
  ri.ode = "../output/manuscript-results/RI/RI-bestfit_ode-params_day250.csv",
  ri.data = "../data/RI-data_2021Sep06.csv",
  ri.const = const.vec.ri,
  pa.beta = "../output/manuscript-results/PA/PA-bestfit_daily-betas_day250.csv",
  pa.ode = "../output/manuscript-results/PA/PA-bestfit_ode-params_day250.csv",
  pa.data = "../data/PA-data_2021Sep06.csv",
  pa.const = const.vec.pa,
  ma.beta = "../output/manuscript-results/MA/MA-bestfit_daily-betas_day250.csv",
  ma.ode = "../output/manuscript-results/MA/MA-bestfit_ode-params_day250.csv",
  ma.data = "../data/MA-data_2021Sep06.csv",
  ma.const = const.vec.ma,
  odepath = "../cpp-v5-discharges-nonhospdeaths/",
  plot.name = "../output/figures/manuscript/ms-figure5.pdf",
  alpha = 25, subsample = 250,
  axis.size = 1, lab.size = 1,
  height = 16, width = 11, logplot = F,
  font.family = "Ubuntu Mono"
)



#################################################################################
### 4. Supplementary Figures (S3-S17)
#################################################################################

#######################################
### Time-Symp-To-Hosp (Section 5.1) ###
#######################################

### Figure S3 (time-symp-to-hosp: loglike comps)

## i. prepare loglikelihood matrix for plotting

# read data
sth.ll.matrix <- read.csv("../output/manuscript-results/RI/symp-to-hosp/RI-STH_loglik.csv")

# create vector for plotting
sth.plot.data <- data.frame(
  text = as.character(
    c(rep(
      "sympt-to-hosp: 2 days",
      length(as.vector(sth.ll.matrix[ , 1]))
    ), rep(
      "sympt-to-hosp: 3 days",
      length(as.vector(sth.ll.matrix[ , 2]))
    ), rep(
      "sympt-to-hosp: 4 days",
      length(as.vector(sth.ll.matrix[ , 3]))
    ), rep(
      "sympt-to-hosp: 5 days",
      length(as.vector(sth.ll.matrix[ , 4]))
    ), rep(
      "sympt-to-hosp: 6 days",
      length(as.vector(sth.ll.matrix[ , 5]))
    ), rep(
      "sympt-to-hosp: 7 days",
      length(as.vector(sth.ll.matrix[ , 6]))
    ))
  ),
  value = c(
    as.vector(sth.ll.matrix[ , 1]), 
    as.vector(sth.ll.matrix[ , 2]),
    as.vector(sth.ll.matrix[ , 3]),
    as.vector(sth.ll.matrix[ , 4]),
    as.vector(sth.ll.matrix[ , 5]),
    as.vector(sth.ll.matrix[ , 6])
  )
)

# density plots for each model
sth.comp.plot <- sth.plot.data %>%
  mutate(text = fct_reorder(text, value)) %>%
  ggplot(aes(y=text, x=value,  fill=text)) +
  geom_density_ridges() +
  scale_fill_viridis(discrete=TRUE) +
  scale_color_viridis(discrete=TRUE) +
  theme_ipsum() +
  theme(
    legend.position="none",
    panel.spacing = unit(1, "lines"),
    strip.text.x = element_text(size = 8)
  ) +
  xlab("Posterior log-likelihood") +
  ylab("")

ggsave(
  filename = "../output/figures/supp-materials/supp-fig03.png",
  plot = sth.comp.plot,
  width = 7, 
  height = 5, 
  units = "in", 
  dpi = 300
)

### Figure S4 (time-symp-to-hosp: trajectory comps)

## i. read data

# ode parameters
ode.sth.2 <- read.csv("../output/manuscript-results/RI/symp-to-hosp/RI-sth-2_ode-params_day250.csv")
ode.sth.3 <- read.csv("../output/manuscript-results/RI/symp-to-hosp/RI-sth-3_ode-params_day250.csv")
ode.sth.4 <- read.csv("../output/manuscript-results/RI/symp-to-hosp/RI-sth-4_ode-params_day250.csv")
ode.sth.5 <- read.csv("../output/manuscript-results/RI/symp-to-hosp/RI-sth-5_ode-params_day250.csv")
ode.sth.6 <- read.csv("../output/manuscript-results/RI/symp-to-hosp/RI-sth-6_ode-params_day250.csv")
ode.sth.7 <- read.csv("../output/manuscript-results/RI/symp-to-hosp/RI-sth-7_ode-params_day250.csv")

# betas
beta.sth.2 <- read.csv("../output/manuscript-results/RI/symp-to-hosp/RI-sth-2_daily-betas_day250.csv")
beta.sth.3 <- read.csv("../output/manuscript-results/RI/symp-to-hosp/RI-sth-3_daily-betas_day250.csv")
beta.sth.4 <- read.csv("../output/manuscript-results/RI/symp-to-hosp/RI-sth-4_daily-betas_day250.csv")
beta.sth.5 <- read.csv("../output/manuscript-results/RI/symp-to-hosp/RI-sth-5_daily-betas_day250.csv")
beta.sth.6 <- read.csv("../output/manuscript-results/RI/symp-to-hosp/RI-sth-6_daily-betas_day250.csv")
beta.sth.7 <- read.csv("../output/manuscript-results/RI/symp-to-hosp/RI-sth-7_daily-betas_day250.csv")

# RI data
ri.data <- read.csv("../data/RI-data_2021Sep06.csv")
ri.data <- data.frame(ri.data)
# restrict data from March 1 to Sept 6 (daynum \in [61,250])
ri.data <- ri.data[(ri.data$daynum >= 61) & (ri.data$daynum <= 250), ]

## ii. store ode and beta parameters in a 3-dim array

# ode parameters
sth.ode.array <- array(c(
    as.matrix(ode.sth.2)[, -c(36:40)], as.matrix(ode.sth.3)[, -c(36:40)], as.matrix(ode.sth.4)[, -c(36:40)],
    as.matrix(ode.sth.5)[, -c(36:40)], as.matrix(ode.sth.6)[, -c(36:40)], as.matrix(ode.sth.7)[, -c(36:40)]
  ), dim = c(nrow(ode.sth.2), ncol(ode.sth.2) - 5, 6))

dimnames(sth.ode.array)[[2]] <- gsub("\\.", "-", colnames(ode.sth.2)[-c(36:40)])
dimnames(sth.ode.array)[[3]] <- c("sth-2", "sth-3", "sth-4", "sth-5", "sth-6", "sth-7")

# beta values
sth.beta.array <- array(c(
    as.matrix(beta.sth.2), as.matrix(beta.sth.3), as.matrix(beta.sth.4),
    as.matrix(beta.sth.5), as.matrix(beta.sth.6), as.matrix(beta.sth.7)
  ), dim = c(nrow(beta.sth.2), ncol(beta.sth.2), 6))


# reporting rate values
sth.rr.array <- array(c(
    as.matrix(ode.sth.2)[, c(36:40)], as.matrix(ode.sth.3)[, c(36:40)], as.matrix(ode.sth.4)[, c(36:40)],
    as.matrix(ode.sth.5)[, c(36:40)], as.matrix(ode.sth.6)[, c(36:40)], as.matrix(ode.sth.7)[, c(36:40)]
  ), dim = c(nrow(ode.sth.2), 5, 6))

## iii. plot trajectories

trajComps(
  beta.array = sth.beta.array, ode.array = sth.ode.array, rr.array = sth.rr.array,
  df = ri.data, loc = "RI", plot.type = "sth",
  plot.name = "../output/figures/supp-materials/supp-fig04.pdf",
  odepath = "../cpp-v5-discharges-nonhospdeaths/",
  col.pal = c("#C83B00", "#BF7100", "#B4982E", "#ACB864", "#ADD39B", "#BEEBCF"),
  height = 14, width = 8
)


### Figure S5 (time-symp-to-hosp: post. histogram comps)

# set-up legend info for plot
sth.rn.1 <- c(
  "Symp-To-Hosp = 2", "Symp-To-Hosp = 3", "Symp-To-Hosp = 4",
  "Symp-To-Hosp = 5", "Symp-To-Hosp = 6", "Symp-To-Hosp = 7"
)

sth.rn.2 <- c("2", "3", "4", "5", "6", "7")

# ode parameters

ode.sth.2 <- read.csv("../output/manuscript-results/RI/symp-to-hosp/RI-sth-2_ode-params_day250.csv")
ode.sth.3 <- read.csv("../output/manuscript-results/RI/symp-to-hosp/RI-sth-3_ode-params_day250.csv")
ode.sth.4 <- read.csv("../output/manuscript-results/RI/symp-to-hosp/RI-sth-4_ode-params_day250.csv")
ode.sth.5 <- read.csv("../output/manuscript-results/RI/symp-to-hosp/RI-sth-5_ode-params_day250.csv")
ode.sth.6 <- read.csv("../output/manuscript-results/RI/symp-to-hosp/RI-sth-6_ode-params_day250.csv")
ode.sth.7 <- read.csv("../output/manuscript-results/RI/symp-to-hosp/RI-sth-7_ode-params_day250.csv")

sth.ode.array <- array(c(
    as.matrix(ode.sth.2)[, -c(36:40)], as.matrix(ode.sth.3)[, -c(36:40)], as.matrix(ode.sth.4)[, -c(36:40)],
    as.matrix(ode.sth.5)[, -c(36:40)], as.matrix(ode.sth.6)[, -c(36:40)], as.matrix(ode.sth.7)[, -c(36:40)]
  ), dim = c(nrow(ode.sth.2), ncol(ode.sth.2) - 5, 6))

dimnames(sth.ode.array)[[2]] <- gsub("\\.", "-", colnames(ode.sth.2)[-c(36:40)])
dimnames(sth.ode.array)[[3]] <- c("sth-2", "sth-3", "sth-4", "sth-5", "sth-6", "sth-7")

# plot posterior histogram comparisons
odeComps(
  ode.array = sth.ode.array, 
  plot.name = "../output/figures/supp-materials/supp-fig05.pdf", 
  main.title = "Time-Symp-To-Hosp Comparison", 
  run.names1 = sth.rn.1, run.names2 = sth.rn.2,
  state = "(Rhode Island)", dens = FALSE,
  col.pal = c("#C83B00", "#BF7100", "#B4982E", "#ACB864", "#ADD39B", "#BEEBCF"), 
  h = 21, w = 15
)


#######################################
### Time-Symp-To-Pres (Section 5.2) ###
#######################################

### Figure S6 (time-symp-to-presentation: loglike comps)

## i. prepare loglikelihood matrix for plotting

# read data
stp.ll.matrix <- read.csv("../output/manuscript-results/RI/symp-to-presentation/RI-STP_loglik.csv")

# create vector for plotting
stp.plot.data <- data.frame(
  text = as.character(
    c(rep(
      "sympt-to-presentation: 0 days",
      length(as.vector(stp.ll.matrix[ , 1]))
    ), rep(
      "sympt-to-presentation: 1 days",
      length(as.vector(stp.ll.matrix[ , 2]))
    ), rep(
      "sympt-to-presentation: 2 days",
      length(as.vector(stp.ll.matrix[ , 3]))
    ), rep(
      "sympt-to-presentation: 3 days",
      length(as.vector(stp.ll.matrix[ , 4]))
    ), rep(
      "sympt-to-presentation: 4 days",
      length(as.vector(stp.ll.matrix[ , 5]))
    ))
  ),
  value = c(
    as.vector(stp.ll.matrix[ , 1]), 
    as.vector(stp.ll.matrix[ , 2]),
    as.vector(stp.ll.matrix[ , 3]),
    as.vector(stp.ll.matrix[ , 4]),
    as.vector(stp.ll.matrix[ , 5])
  )
)

# density plots for each model
stp.comp.plot <- stp.plot.data %>%
  mutate(text = fct_reorder(text, value)) %>%
  ggplot(aes(y=text, x=value,  fill=text)) +
  geom_density_ridges() +
  scale_fill_viridis(discrete=TRUE) +
  scale_color_viridis(discrete=TRUE) +
  theme_ipsum() +
  theme(
    legend.position="none",
    panel.spacing = unit(1, "lines"),
    strip.text.x = element_text(size = 8)
  ) +
  xlab("Posterior log-likelihood") +
  ylab("")
  

ggsave(
  filename = "../output/figures/supp-materials/supp-fig06.png",
  plot = stp.comp.plot,
  width = 7, 
  height = 5, 
  units = "in", 
  dpi = 300
)

### Figure S7 (time-symp-to-presentation: trajectory comps)

## i. read data

# ode parameters
ode.d0 <- read.csv("../output/manuscript-results/RI/symp-to-presentation/RI-stp-d0_ode-params_day250.csv")
ode.d1 <- read.csv("../output/manuscript-results/RI/symp-to-presentation/RI-stp-d1_ode-params_day250.csv")
ode.d2 <- read.csv("../output/manuscript-results/RI/symp-to-presentation/RI-stp-d2_ode-params_day250.csv")
ode.d3 <- read.csv("../output/manuscript-results/RI/symp-to-presentation/RI-stp-d3_ode-params_day250.csv")
ode.d4 <- read.csv("../output/manuscript-results/RI/symp-to-presentation/RI-stp-d4_ode-params_day250.csv")

# betas
beta.d0 <- read.csv("../output/manuscript-results/RI/symp-to-presentation/RI-stp-d0_daily-betas_day250.csv")
beta.d1 <- read.csv("../output/manuscript-results/RI/symp-to-presentation/RI-stp-d1_daily-betas_day250.csv")
beta.d2 <- read.csv("../output/manuscript-results/RI/symp-to-presentation/RI-stp-d2_daily-betas_day250.csv")
beta.d3 <- read.csv("../output/manuscript-results/RI/symp-to-presentation/RI-stp-d3_daily-betas_day250.csv")
beta.d4 <- read.csv("../output/manuscript-results/RI/symp-to-presentation/RI-stp-d4_daily-betas_day250.csv")

# actual RI covid data
ri.data <- read.csv("../data/RI-data_2021Sep06.csv")
ri.data <- data.frame(ri.data)
# restrict data from March 1 to Sept 6 (daynum \in [61,250]) 
ri.data <- ri.data[(ri.data$daynum >= 61) & (ri.data$daynum <= 250),]

## ii. store ode and beta parameters in a 3-dim array

# ode parameters
stp.ode.array <- array(c(
    as.matrix(ode.d0)[, -c(36:40)], as.matrix(ode.d1)[, -c(36:40)], as.matrix(ode.d2)[, -c(36:40)],
    as.matrix(ode.d3)[, -c(36:40)], as.matrix(ode.d4)[, -c(36:40)]
  ), dim = c(nrow(ode.d0), ncol(ode.d0) - 5, 5)
)

dimnames(stp.ode.array)[[2]] <- gsub("\\.", "-", colnames(ode.d0)[-c(36:40)])
dimnames(stp.ode.array)[[3]] <- c("delay-0", "delay-1", "delay-2", "delay-3", "delay-4")

# beta values
stp.beta.array <- array(c(
    as.matrix(beta.d0), as.matrix(beta.d1), as.matrix(beta.d2),
    as.matrix(beta.d3), as.matrix(beta.d4)
  ), dim = c(nrow(beta.d0), ncol(beta.d0), 5)
)

# reporting rate values
stp.rr.array <- array(c(
    as.matrix(ode.d0)[, c(36:40)], as.matrix(ode.d1)[, c(36:40)], as.matrix(ode.d2)[, c(36:40)],
    as.matrix(ode.d3)[, c(36:40)], as.matrix(ode.d4)[, c(36:40)]
  ), dim = c(nrow(ode.d0), 5, 5)
)

## iii. plot trajectories

trajComps(
  beta.array = stp.beta.array, ode.array = stp.ode.array, rr.array = stp.rr.array,
  df = ri.data, loc = "RI", plot.type = "stp",
  plot.name = "../output/figures/supp-materials/supp-fig07.pdf",
  odepath = "../cpp-v5-discharges-nonhospdeaths/",
  col.pal = c("#C83B00", "#BF7100", "#B4982E", "#ACB864", "#ADD39B", "#BEEBCF"),
  height = 14, width = 8
)

### Figure S8 (time-symp-to-presentation: post. histogram comps)

# set-up legend info for plot
stp.rn.1 <- c(
  "Symp-To-Pres = 0", "Symp-To-Pres = 1", "Symp-To-Pres = 2",
  "Symp-To-Pres = 3", "Symp-To-Pres = 4"
)

stp.rn.2 <- c("0", "1", "2", "3", "4")

# ode parameters
ode.d0 <- read.csv("../output/manuscript-results/RI/symp-to-presentation/RI-stp-d0_ode-params_day250.csv")
ode.d1 <- read.csv("../output/manuscript-results/RI/symp-to-presentation/RI-stp-d1_ode-params_day250.csv")
ode.d2 <- read.csv("../output/manuscript-results/RI/symp-to-presentation/RI-stp-d2_ode-params_day250.csv")
ode.d3 <- read.csv("../output/manuscript-results/RI/symp-to-presentation/RI-stp-d3_ode-params_day250.csv")
ode.d4 <- read.csv("../output/manuscript-results/RI/symp-to-presentation/RI-stp-d4_ode-params_day250.csv")

stp.ode.array <- array(c(
    as.matrix(ode.d0)[, -c(36:40)], as.matrix(ode.d1)[, -c(36:40)], as.matrix(ode.d2)[, -c(36:40)],
    as.matrix(ode.d3)[, -c(36:40)], as.matrix(ode.d4)[, -c(36:40)]
  ), dim = c(nrow(ode.d0), ncol(ode.d0) - 5, 5)
)

dimnames(stp.ode.array)[[2]] <- gsub("\\.", "-", colnames(ode.d0)[-c(36:40)])
dimnames(stp.ode.array)[[3]] <- c("delay-0", "delay-1", "delay-2", "delay-3", "delay-4")

# plot posterior histogram comparisons
odeComps(
  ode.array = stp.ode.array, 
  plot.name = "../output/figures/supp-materials/supp-fig08.pdf", 
  main.title = "Time-Symp-To-Pres Comparison", 
  run.names1 = stp.rn.1, run.names2 = stp.rn.2,
  state = "(Rhode Island)", dens = FALSE,
  col.pal = c("#C83B00", "#BF7100", "#B4982E", "#ACB864", "#ADD39B", "#BEEBCF"), 
  h = 21, w = 15
)

#################################
### Excess Deaths (Section 6) ###
#################################

### Figure S9 (excess death posterior histograms - MA)
summaryHists(
  ode.directory = "../output/manuscript-results/MA/excess-deaths/MA-ex-deaths_ode-params_day250.csv",
  ode.inds = 1:35,
  plot.name = "../output/figures/supp-materials/supp-fig09.pdf",
  loc = "MA", height = 16, width = 16,
  font.family = "Ubuntu Mono"
)

### Figure S10 (excess death posterior histograms - PA)
summaryHists(
  ode.directory = "../output/manuscript-results/PA/excess-deaths/PA-ex-deaths_ode-params_day250.csv",
  ode.inds = 1:34,
  plot.name = "../output/figures/supp-materials/supp-fig10.pdf",
  loc = "PA", height = 16, width = 16,
  font.family = "Ubuntu Mono"
)

### Figure S11 (excess deaths trajectories - MA)

## i. read in data

# ode parameters
ma.ode.best <- read.csv("../output/manuscript-results/MA/MA-bestfit_ode-params_day250.csv")
ma.ode.excess <- read.csv("../output/manuscript-results/MA/excess-deaths/MA-ex-deaths_ode-params_day250.csv")

# betas
ma.beta.best <- read.csv("../output/manuscript-results/MA/MA-bestfit_daily-betas_day250.csv")
ma.beta.excess <- read.csv("../output/manuscript-results/MA/excess-deaths/MA-ex-deaths_daily-betas_day250.csv")

# MA data
ma.data <- read.csv("../data/MA-data_2021Sep06.csv")
ma.data <- data.frame(ma.data)
# restrict data from March 1 to Sept 6 (daynum \in [61,250]) 
ma.data <- ma.data[(ma.data$daynum >= 61) & (ma.data$daynum <= 250),]

## ii. store ode and beta parameters in a 3-dim array

# ode parameters
ma.ode.array <- array(c(as.matrix(ma.ode.best)[, -36], as.matrix(ma.ode.excess)[, -36]),
  dim = c(nrow(ma.ode.best), ncol(ma.ode.best) - 1, 2)
)

dimnames(ma.ode.array)[[2]] <- gsub("\\.", "-", colnames(ma.ode.best)[-36])
dimnames(ma.ode.array)[[3]] <- c("Observed Deaths", "Excess Deaths")

# beta values
ma.beta.array <- array(c(as.matrix(ma.beta.best), as.matrix(ma.beta.excess)),
  dim = c(nrow(ma.beta.best), ncol(ma.beta.best), 2)
)

# reporting rate values
ma.rr.array <- array(c(as.matrix(ma.ode.best)[, 36], as.matrix(ma.ode.excess)[, 36]),
  dim = c(nrow(ma.ode.best), 1, 2)
)

## iii. generate trajectory plots

trajComps(
  beta.array = ma.beta.array, ode.array = ma.ode.array, rr.array = ma.rr.array,
  df = ma.data, loc = "MA", plot.type = "ex-deaths",
  plot.name = "../output/figures/supp-materials/supp-fig11.pdf",
  odepath = "../cpp-v5-discharges-nonhospdeaths/",
  col.pal = c("#C83B00", "#BEEBCF"),
  height = 14, width = 8
)

### Figure S12 (excess deaths trajectories - PA)

## i. read in data

# ode parameters
pa.ode.best <- read.csv("../output/manuscript-results/PA/PA-bestfit_ode-params_day250.csv")
pa.ode.excess <- read.csv("../output/manuscript-results/PA/excess-deaths/PA-ex-deaths_ode-params_day250.csv")

# betas
pa.beta.best <- read.csv("../output/manuscript-results/PA/PA-bestfit_daily-betas_day250.csv")
pa.beta.excess <- read.csv("../output/manuscript-results/PA/excess-deaths/PA-ex-deaths_daily-betas_day250.csv")

# MA data
pa.data <- read.csv("../data/PA-data_2021Sep06.csv")
pa.data <- data.frame(pa.data)
# restrict data from March 1 to Sept 6 (daynum \in [61,250]) 
pa.data <- pa.data[(pa.data$daynum >= 61) & (pa.data$daynum <= 250),]

## ii. store ode and beta parameters in a 3-dim array

# ode parameters
pa.ode.array <- array(c(as.matrix(pa.ode.best)[, -c(35:39)], as.matrix(pa.ode.excess)[, -c(35:39)]),
  dim = c(nrow(pa.ode.best), ncol(pa.ode.best) - 5, 2)
)

dimnames(pa.ode.array)[[2]] <- gsub("\\.", "-", colnames(pa.ode.best)[-c(35:39)])
dimnames(pa.ode.array)[[3]] <- c("Reported Deaths", "Excess Deaths")

# beta values
pa.beta.array <- array(c(as.matrix(pa.beta.best), as.matrix(pa.beta.excess)),
  dim = c(nrow(pa.beta.best), ncol(pa.beta.best), 2)
)

# reporting rate values
pa.rr.array <- array(c(as.matrix(pa.ode.best)[, c(35:39)], as.matrix(pa.ode.excess)[, c(35:39)]),
  dim = c(nrow(pa.ode.best), 5, 2)
)

## iii. generate trajectory plots
trajComps(
  beta.array = pa.beta.array, ode.array = pa.ode.array, rr.array = pa.rr.array,
  df = pa.data, loc = "PA", plot.type = "ex-deaths",
  plot.name = "../output/figures/supp-materials/supp-fig12.pdf",
  odepath = "../cpp-v5-discharges-nonhospdeaths/",
  col.pal = c("#C83B00", "#BEEBCF"),
  height = 14, width = 8
)


####################################
### Final Model Fits (Section 7) ###
####################################

### Figure S13 (MA posterior trajectories)

set.seed(31)

const.vec.ma <- c(2, "", 4)
names(const.vec.ma) = c("steps-per-day", "symp-frac-davies", "time-symp-to-hosp")

fig2panel(
  beta.directory = "../output/manuscript-results/MA/MA-bestfit_daily-betas_day250.csv",
  ode.directory = "../output/manuscript-results/MA/MA-bestfit_ode-params_day250.csv",
  data.directory = "../data/MA-data_2021Sep06.csv",
  odepath = "../cpp-v5-discharges-nonhospdeaths/",
  loc = "MA", const = const.vec.ma,
  plot.name = "../output/figures/supp-materials/supp-fig13.pdf", 
  alpha = 50,
  subsample = 250,
  axis.size = 1.75, lab.size = 1.75,
  ncol = 2, nrow = 4, height = 12, width = 12,
  font.family = "Ubuntu Mono", delay = 2
)


### Figure S14 (PA posterior trajectories)

set.seed(302)

const.vec.pa <- c("", 0.8, 2, 4)
names(const.vec.pa) <- c("symp-frac-davies", "prob-icu-vent", "steps-per-day", "time-symp-to-hosp")

fig2panel(
  beta.directory = "../output/manuscript-results/PA/PA-bestfit_daily-betas_day250.csv",
  ode.directory = "../output/manuscript-results/PA/PA-bestfit_ode-params_day250.csv",
  data.directory = "../data/PA-data_2021Sep06.csv",
  odepath = "../cpp-v5-discharges-nonhospdeaths/",
  loc = "PA", const = const.vec.pa,
  plot.name = "../output/figures/supp-materials/supp-fig14.pdf", 
  alpha = 50,
  subsample = 250,
  axis.size = 1.75, lab.size = 1.75,
  ncol = 2, nrow = 4, height = 12, width = 12,
  font.family = "Ubuntu Mono", delay = 2
)

### Figure S15 (RI posterior histograms)
summaryHists(
  ode.directory = "../output/manuscript-results/RI/RI-bestfit_ode-params_day250.csv",
  ode.inds = 1:35,
  plot.name = "../output/figures/supp-materials/supp-fig15.pdf",
  loc = "RI", height = 16, width = 16,
  font.family = "Ubuntu Mono"
)

### Figure S16 (MA posterior histograms)
summaryHists(
  ode.directory = "../output/manuscript-results/MA/MA-bestfit_ode-params_day250.csv",
  ode.inds = 1:35,
  plot.name = "../output/figures/supp-materials/supp-fig16.pdf",
  loc = "MA", height = 16, width = 16,
  font.family = "Ubuntu Mono"
)

### Figure S17 (PA posterior histograms)
summaryHists(
  ode.directory = "../output/manuscript-results/PA/PA-bestfit_ode-params_day250.csv",
  ode.inds = 1:34,
  plot.name = "../output/figures/supp-materials/supp-fig17.pdf",
  loc = "PA", height = 16, width = 16,
  font.family = "Ubuntu Mono"
)
