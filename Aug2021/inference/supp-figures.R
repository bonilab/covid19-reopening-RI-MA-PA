#!/usr/bin/env Rscript

### supp-figures.R
### last edited: 15 Nov 2021
### authors: Nathan Wikle
###
### Functions used to create figures found in the supplementary materials.

############################################################
### Compare ODESIM trajectories from various model fits. ###
############################################################
trajComps <- function(
  beta.array, ode.array, rr.array, df,
  plot.type, loc,
  plot.name,                    
  odepath = "./cpp-v5-discharges-nonhospdeaths/",
  col.pal = c("#C83B00", "#BEEBCF"),
  nrow = 5, ncol = 2, ...
){
  
  ### process data
  
  # number or days
  n.days <- nrow(df)
  days <- df$daynum
  end.day <- max(days, na.rm = TRUE) + 1
  num.days <- end.day - 60 
  
  # cubic spline with one basis function every 7 days
  bspl <- create.bspline.basis(c(61, end.day), nbasis = round(num.days / 7))
  Z <- eval.basis(bspl, 61:end.day)
  
  if (loc == "MA"){
    # zero-order spline with one rate
    bspl.rr <- create.bspline.basis(c(61, end.day), nbasis = 1, norder = 1)
    Z.rr <- eval.basis(bspl.rr, 61:end.day)
    Z.rr <- Z.rr[-nrow(Z.rr), ]
  } else {
    # iSpline for reporting rate
    knots <- c(84, 92, 100, 108, 130, 160, 190)
    Z.rr <- iSpline(x = 61:end.day, knots = knots, degree = 2)
    Z.rr <- cbind(rep(1, nrow(Z.rr)), Z.rr)
    Z.rr <- Z.rr[-nrow(Z.rr), c(1, 4:7)]
  }
  
  # introday
  if (loc == "PA") {
    id = 60
  } else {
    id = 55
  }
  
  # list structure with all relevant data
  structures <- list()
  for (k in 1:dim(beta.array)[3]){
    
    if (plot.type == "sth") {
      const <- c("", 2, k + 1)
      names(const) <- c("symp-frac-davies", "steps-per-day", "time-symp-to-hosp")
      delay <- 2
    } else if (plot.type == "stp") {
      const <- c("", 2, 4)
      names(const) <- c("symp-frac-davies", "steps-per-day", "time-symp-to-hosp")
      delay <- k - 1
    } else if (plot.type == "ex-deaths") {
      if (loc == "PA") {
        const <- c("", 0.8, 2, 4)
        names(const) <- c("symp-frac-davies", "prob-icu-vent", "steps-per-day", "time-symp-to-hosp")
      } else {
        const <- c("", 2, 4)
        names(const) <- c("symp-frac-davies", "steps-per-day", "time-symp-to-hosp")
      }
      delay <- 2
    }
    
    structures[[k]] <- list(
      betas = beta.array[, , k],
      ode = ode.array[, , k],
      rr = rr.array[, , k],
      days = days,
      Z.beta = Z,
      Z.rr = Z.rr,
      df = df,
      loc = loc,
      pres.delay = delay,
      const = const,
      introday = id,
      sf.choice = FALSE
    )
  }
  
  # generate trajectories using these structures
  traj <- list()
  for(k in 1:dim(beta.array)[3]){
    traj[[k]] <- traj.sim(
      samples = structures[[k]],
      odepath = odepath,
      csv = TRUE
    )
  }
  
  # process actual data
  dp <- data.fine.process(structures[[1]]$df, loc = structures[[1]]$loc)
  
  ### generate plots
  pdf(file = plot.name, ...)
  
  # grid of plots
  par(mfrow = c(nrow,ncol))
  par(mar = c(3, 2, 1, 1), oma = c(2, 2, 3, 0.5))

  # vector of dates
  dates <- seq(as.Date("2020/1/1"), by = "day", length.out = 365)
  
  # last day of simulated output
  end.day <- max(structures[[1]]$days)
  
  n.days <- ncol(traj[[1]]$beta.full)
  n.runs <- dim(rr.array)[3]
  
  ## A) Mixing Parameter (betas)
  
  beta.medians <- matrix(NA_real_, nrow = n.days, ncol = n.runs)
  for (k in 1:n.runs) {
    beta.medians[, k] <- apply(traj[[k]]$beta.full, 2, median)
  }
  
  reparam <- function(betas) {
    avg <- mean(betas[5:15])
    betas / avg
  }
  
  new.betas <- apply(beta.medians, 2, reparam)

  # max and min
  y.max <- min(c(max(new.betas), 4))

  # plot contact rate (beta) over time
  plot(dates[c(end.day - ncol(traj[[1]]$rr.full) + 1):end.day],
    new.betas[, 1],
    type = "l", col = col.pal[1], lwd = 1,
    ylim = c(0, y.max), ylab = "", xlab = "", bty = "l"
  )

  for (k in 2:dim(rr.array)[3]) {
    lines(dates[c(end.day - ncol(traj[[1]]$rr.full) + 1):end.day],
      new.betas[, k],
      col = col.pal[k], lwd = 1.5
    )
  }
  
  title(main = "A. ", adj = 0.05, line = -0.5, cex.main = 1.25, outer = FALSE)
  title(
    main = "Mixing Parameter (beta)", adj = 1,
    line = -0.5, cex.main = 1.25, outer = FALSE
  )
  

  if (plot.type == "sth") {
    legend(dates[80], 3.75,
      legend = c("2", "3", "4", "5", "6", "7"),
      bty = "n", lty = 1, lwd = 3, col = col.pal, title = "Time-Symp-To-Hosp:",
      cex = 1.25
    )
  } else if (plot.type == "stp") {
    legend(dates[80], 3.75,
      legend = c("0", "1", "2", "3", "4"),
      bty = "n", lty = 1, lwd = 3, col = col.pal, title = "Time-Symp-To-Pres:",
      cex = 1.25
    )
  } else if (plot.type == "ex-deaths") {
    legend(dates[80], 3.75,
      legend = c("Reported Deaths", "Reported + Excess Deaths"),
      bty = "n", lty = 1, lwd = 3, col = col.pal, title = "", cex = 1.25
    )
  }
  
  ## B) Reporting Rate
  
  # calculate median reporting rates
  rr.medians <- matrix(NA_real_, nrow = n.days, ncol = n.runs)
  for (k in 1:n.runs) {
    rr.medians[, k] <- apply(traj[[k]]$rr.full, 2, median)
  }
  
  # max and min
  y.max <- max(rr.medians); y.min <- min(rr.medians)
  
  # plot reporting rate over time
  plot(dates[c(end.day - ncol(traj[[1]]$rr.full) + 1):end.day],
    rr.medians[, 1],
    type = "l", col = col.pal[1], lwd = 2,
    ylim = c(0, 1), ylab = "", xlab = "", bty = "l"
  )
  
  for (k in 2:n.runs) {
    lines(dates[c(end.day - ncol(traj[[1]]$rr.full) + 1):end.day],
      rr.medians[, k], col = col.pal[k], lwd = 2
    )
  }
  
  title(main = "B. ", adj = 0.05, line = -0.5, cex.main = 1.25, outer = FALSE)
  title(
    main = "Symptomatic Reporting Rate Over Time", adj = 1,
    line = -0.5, cex.main = 1.25, outer = FALSE
  )
  
  ## C) New Cases
  
  # calculate median new cases
  rr.newcases <- matrix(NA_real_, nrow = ncol(traj[[k]]$tot.sympt.new), ncol = n.runs)
  for (k in 1:n.runs) {
    rr.newcases[, k] <- apply(traj[[k]]$tot.sympt.new, 2, median)
  }
  
  y.max <- max(rr.newcases)
  
  plot(dates[traj[[1]]$days.odesim],
       rr.newcases[,1], type = "l", col = col.pal[1], lwd = 2,
       ylim = c(0, y.max), ylab = "", xlab = "", bty = "l")
  
  for (k in 2:n.runs) {
    lines(dates[traj[[1]]$days.odesim],
      rr.newcases[, k], col = col.pal[k], lwd = 2
    )
  }
  
  points(dates[dp$days], dp$tot.sympt.new, pch = 19, col = "black", cex = 0.5)
  
  title(main = "C. ", adj = 0.05, line = -0.5, cex.main = 1.25, outer = FALSE)
  title(
    main = "New Symptomatic Cases (Reported)", adj = 1,
    line = -0.5, cex.main = 1.25, outer = FALSE
  )
  
  ## D) New Hospitalizations
  
  # calculate median new cases
  new.hosps.median <- matrix(NA_real_, nrow = ncol(traj[[k]]$tot.hosp.new), ncol = n.runs)
  for (k in 1:n.runs) {
    new.hosps.median[, k] <- apply(traj[[k]]$tot.hosp.new, 2, median)
  }
  
  y.max <- max(new.hosps.median)
  
  plot(dates[traj[[1]]$days.odesim],
    new.hosps.median[, 1],
    type = "l", col = col.pal[1], lwd = 2,
    ylim = c(0, y.max), ylab = "", xlab = "", bty = "l"
  )
  
  for (k in 2:n.runs) {
    lines(dates[traj[[1]]$days.odesim],
      new.hosps.median[, k],
      col = col.pal[k], lwd = 2
    )
  }
  
  if (loc == "RI") {
    points(dates[dp$days], dp$tot.hosp.new, pch = 19, col = "black", cex = 0.5)
  }
  
  title(main = "D. ", adj = 0.05, line = -1, cex.main = 1.25, outer = FALSE)
  title(
    main = "New Hospitalizations", adj = 1,
    line = -1, cex.main = 1.25, outer = FALSE
  )
  
  ## E) Current Hospitalizations
  
  # calculate median current hosps
  curr.hosps.median <- matrix(NA_real_, nrow = ncol(traj[[k]]$tot.hosp.curr), ncol = n.runs)
  for (k in 1:n.runs) {
    curr.hosps.median[, k] <- apply(traj[[k]]$tot.hosp.curr, 2, median)
  }
  
  y.max <- max(curr.hosps.median)
  
  plot(dates[traj[[1]]$days.odesim],
    curr.hosps.median[, 1],
    type = "l", col = col.pal[1], lwd = 2,
    ylim = c(0, y.max), ylab = "", xlab = "", bty = "l"
  )
  
  for (k in 2:n.runs) {
    lines(dates[traj[[1]]$days.odesim],
      curr.hosps.median[, k],
      col = col.pal[k], lwd = 2
    )
  }
  
  points(dates[dp$days], dp$tot.hosp.curr, pch = 19, col = "black", cex = 0.5)
  
  title(main = "E. ", adj = 0.05, line = -0.5, cex.main = 1.25, outer = FALSE)
  title(
    main = "Current Hospitalizations", adj = 1,
    line = -0.5, cex.main = 1.25, outer = FALSE
  )
  
  ## F) Current ICU Occupancy
  
  # calculate median current icu
  curr.icu.median <- matrix(NA_real_, nrow = ncol(traj[[k]]$tot.icu.curr), ncol = n.runs)
  for (k in 1:n.runs) {
    curr.icu.median[, k] <- apply(traj[[k]]$tot.icu.curr, 2, median)
  }
  
  y.max <- max(curr.icu.median)
  
  plot(dates[traj[[1]]$days.odesim],
    curr.icu.median[, 1],
    type = "l", col = col.pal[1], lwd = 2,
    ylim = c(0, y.max), ylab = "", xlab = "", bty = "l"
  )
  
  for (k in 2:n.runs) {
    lines(dates[traj[[1]]$days.odesim],
      curr.icu.median[, k],
      col = col.pal[k], lwd = 2
    )
  }
  
  if (loc != "PA") {
    points(dates[dp$days], dp$tot.icu.curr, pch = 19, col = "black", cex = 0.5)
  }
  
  title(main = "F. ", adj = 0.05, line = -0.5, cex.main = 1.25, outer = FALSE)
  title(
    main = "Current ICU Occupancy", adj = 1,
    line = -0.5, cex.main = 1.25, outer = FALSE
  )
  
  ## G) Current Vent. Occupancy
  
  # calculate median current hosps
  curr.vent.median <- matrix(NA_real_, nrow = ncol(traj[[k]]$tot.vent.curr), ncol = n.runs)
  for (k in 1:n.runs) {
    curr.vent.median[, k] <- apply(traj[[k]]$tot.vent.curr, 2, median)
  }
  
  y.max <- max(c(curr.vent.median, dp$tot.vent.curr), na.rm = T)

  plot(dates[traj[[1]]$days.odesim],
    curr.vent.median[, 1],
    type = "l", col = col.pal[1], lwd = 2,
    ylim = c(0, y.max), ylab = "", xlab = "", bty = "l"
  )
  
  for (k in 2:n.runs) {
    lines(dates[traj[[1]]$days.odesim],
      curr.vent.median[, k],
      col = col.pal[k], lwd = 2
    )
  }
  
  points(dates[dp$days], dp$tot.vent.curr, pch = 19, col = "black", cex = 0.5)
  
  title(main = "G. ", adj = 0.05, line = -0.5, cex.main = 1.25, outer = FALSE)
  title(
    main = "Current Ventilator Occupancy", adj = 1,
    line = -0.5, cex.main = 1.25, outer = FALSE
  )
  
  ## H) New Deaths
  
  # calculate median new deaths
  curr.deaths.median <- matrix(NA_real_, nrow = ncol(traj[[k]]$tot.deaths.new), ncol = n.runs)
  for (k in 1:n.runs) {
    curr.deaths.median[, k] <- apply(traj[[k]]$tot.deaths.new, 2, median)
  }
  
  y.max <- max(c(curr.deaths.median, dp$tot.deaths.new), na.rm = T)
  
  plot(dates[traj[[1]]$days.odesim],
    curr.deaths.median[, 1],
    type = "l", col = col.pal[1], lwd = 2,
    ylim = c(0, y.max), ylab = "", xlab = "", bty = "l"
  )
  
  for (k in 2:n.runs) {
    lines(dates[traj[[1]]$days.odesim],
      curr.deaths.median[, k],
      col = col.pal[k], lwd = 2
    )
  }
  
  points(dates[dp$days], dp$tot.deaths.new, pch = 19, col = "black", cex = 0.5)
  
  title(main = "H. ", adj = 0.05, line = -0.5, cex.main = 1.25, outer = FALSE)
  title(
    main = "New Deaths", adj = 1,
    line = -0.5, cex.main = 1.25, outer = FALSE
  )
  
  ## I) New Hospital Deaths
  
  # calculate median new hosp. deaths
  hosp.deaths.median <- matrix(NA_real_, nrow = ncol(traj[[k]]$tot.hosp.deaths.new), ncol = n.runs)
  for (k in 1:n.runs) {
    hosp.deaths.median[, k] <- apply(traj[[k]]$tot.hosp.deaths.new, 2, median)
  }
  
  y.max <- max(c(hosp.deaths.median, dp$tot.hosp.deaths.new), na.rm = T)
  
  plot(dates[traj[[1]]$days.odesim],
    hosp.deaths.median[, 1],
    type = "l", col = col.pal[1], lwd = 2,
    ylim = c(0, y.max), ylab = "", xlab = "", bty = "l"
  )
  
  for (k in 2:n.runs) {
    lines(dates[traj[[1]]$days.odesim],
      hosp.deaths.median[, k],
      col = col.pal[k], lwd = 2
    )
  }
  
  if (loc == "RI"){
    points(dates[dp$days], dp$tot.hosp.deaths.new, pch = 19, col = "black", cex = 0.5)
  }
  
  title(main = "I. ", adj = 0.05, line = -0.5, cex.main = 1.25, outer = FALSE)
  title(
    main = "New Hospital Deaths", adj = 1,
    line = -0.5, cex.main = 1.25, outer = FALSE
  )
  
  ## J. New Hospital Discharges
  
  # calculate median new hosp. deaths
  hosp.discharges.median <- matrix(NA_real_, nrow = ncol(traj[[k]]$tot.hosp.discharges.new), ncol = n.runs)
  for (k in 1:n.runs) {
    hosp.discharges.median[, k] <- apply(traj[[k]]$tot.hosp.discharges.new, 2, median)
  }
  
  y.max <- max(c(hosp.discharges.median, dp$tot.hosp.discharges.new), na.rm = T)
  
  plot(dates[traj[[1]]$days.odesim],
    hosp.discharges.median[, 1],
    type = "l", col = col.pal[1], lwd = 2,
    ylim = c(0, y.max), ylab = "", xlab = "", bty = "l"
  )
  
  for (k in 2:n.runs) {
    lines(dates[traj[[1]]$days.odesim],
      hosp.discharges.median[, k],
      col = col.pal[k], lwd = 2
    )
  }
  
  if (loc == "RI") {
    points(dates[dp$days], dp$tot.hosp.discharges.new, pch = 19, col = "black", cex = 0.5)
  }
  
  title(main = "J. ", adj = 0.05, line = -0.5, cex.main = 1.25, outer = FALSE)
  title(
    main = "New Hospital Discharges", adj = 1,
    line = -0.5, cex.main = 1.25, outer = FALSE
  )
  

  if (plot.type == "sth") {
    title(
      main = "Time-Symp-To-Hosp Model Fit Comparison", cex.main = 2,
      adj = 0, line = 1.5, outer = T
    )
    title(
      main = "(Rhode Island)", cex.main = 2,
      adj = 1, line = 1.5, outer = T
    )
  } else if (plot.type == "stp") {
    title(
      main = "Time-Symp-To-Presentation Model Fit Comparison", cex.main = 2,
      adj = 0, line = 1.5, outer = T
    )
    title(
      main = "(Rhode Island)", cex.main = 2,
      adj = 1, line = 1.5, outer = T
    )
  } else if (plot.type == "ex-deaths") {
    if (loc == "MA") {
      title(
        main = "Excess Deaths: Massachusetts", cex.main = 2,
        adj = 0, line = 1.5, outer = TRUE
      )
    } else if (loc == "PA") {
      title(
        main = "Excess Deaths: Pennsylvania", cex.main = 2,
        adj = 0, line = 1.5, outer = TRUE
      )
    }
  }
  dev.off()
}


#################################################
### Compare ODE histograms across model runs. ###
#################################################
odeComps <- function(
  ode.array, plot.name, main.title, run.names1, run.names2,
  state = "(Rhode Island)", dens = FALSE,
  col.pal = c("#C83B00", "#BF7100", "#B4982E", "#ACB864", "#ADD39B", "#BEEBCF"), 
  n.plot.per.page = 3, 
  ...
){

  # grab ode parameter names
  ode.param.names <- dimnames(ode.array)[[2]] 

  {pdf(file = plot.name, ...)

    ### left, right, bottom, top
    m0 <- matrix(c(0, 1, 0, 1), nrow = 1, ncol = 4, byrow = T)

    if (n.plot.per.page == 3){
      m1 <- rbind(c(0.01, 0.99, 0.65, 0.95),
                  c(0.01, 0.99, 0.335, 0.635),
                  c(0.01, 0.99, 0.02, 0.32))

      m2 <- rbind(c(0.07, 0.37, 0.5, 0.97),
                  c(0.37, 0.67, 0.5, 0.97),
                  c(0.67, 0.97, 0.5, 0.97),
                  c(0.07, 0.37, 0.03, 0.5),
                  c(0.37, 0.67, 0.03, 0.5),
                  c(0.67, 0.97, 0.03, 0.5))
    }

    for (j in 1:dim(ode.array)[2]){

      plot.num <- j %% n.plot.per.page

      if (plot.num == 0){
        plot.num <- n.plot.per.page
      }

      if (plot.num == 1){

        if (j > 1){
          close.screen(all = TRUE)
        }

        # first set of plots
        split.screen(m0)
        split.screen(m1, screen = 1)

        for (n.p in 1:n.plot.per.page){
          split.screen(m2, screen = 1 + n.p)
        } 
        # main screen (title)
        screen(1)
        par(mar = c(0,0,3,0),
            oma = c(0,0,0,0))

        if (j == 1){
          title(main = main.title, line = 0, adj = 0.05, outer = FALSE, cex.main = 3.5)
          title(main = state, line = 0, adj = 0.95, outer = FALSE, cex.main = 3.5)
        }
      } 

      # screen for set of plots
      screen(plot.num + 1) 

      # bottom, left, top, right
      par(mar = c(3,3,5,3), 
          oma = c(0,0,0,0))

      title(main = gsub("\\.", "-", ode.param.names[j]), adj = 0,
            line = 3.5, cex.main = 2.5, outer = FALSE)

      # kernel density estimates for each run
      den.str <- apply(ode.array[,j,], 2, density)

      # minimum and maximum x and y values
      min.x <- min(unlist(lapply(den.str, function(d){min(d[[1]])})))
      max.x <- max(unlist(lapply(den.str, function(d){max(d[[1]])})))
      min.y <- min(unlist(lapply(den.str, function(d){min(d[[2]])})))
      max.y <- max(unlist(lapply(den.str, function(d){max(d[[2]])})))

      for(k in 1:length(den.str)){
        screen(((plot.num - 1) * 6 + n.plot.per.page + 1) + k)
        par(mar = c(2,2,1,1),
            oma = c(0,0,1,0))

        if (dens){
          # density plots
          plot(den.str[[k]], xlim = c(min.x, max.x), ylim = c(0, max.y),
              bty = "L", xlab = "", ylab = "", col = alpha(col.pal[k], alpha = 0.95), lwd = 3, main = "")
          polygon(den.str[[k]], col= alpha(col.pal[k], alpha = 0.5),  border=col.pal[k]) 
        } else {
          # histograms
          hist(ode.array[,j,k], xlim = c(min.x, max.x), ylim = c(0, max.y),
               bty = "L", xlab = "", ylab = "", col = alpha(col.pal[k], alpha = 0.95), 
               main = "", probability = "TRUE")
        }
        
        if (j == 1){
          title(main = run.names1[k], line = -2, adj = 0.05, cex.main = 1.5)
        } else {
          title(main = run.names2[k], line = -2, adj = 0.05, cex.main = 1.5)
        }
        
      }
    }
    
    close.screen(all = TRUE)
    dev.off()
  }
}


###################################################################
### Create the histogram plots in Section 7 (Final Model Fits). ###
###################################################################
summaryHists <- function(
  ode.directory,
  ode.inds,
  plot.name,
  loc,
  axis.size = 1.75, lab.size = 1.75,
  ncol = 6, nrow = 6, height = 12, width = 12,
  font.family = "Helvetica"
) {

  # read ODE csv file
  ode <- read.csv(ode.directory)
  ode.final <- ode[, ode.inds]

  pdf(file = plot.name, family = font.family, h = height, w = width)

  # determine set-up
  par(mfrow = c(nrow, ncol))
  par(
    mar = c(2, 2, 3, 1),
    oma = c(2, 2, 3, 2)
  )

  for (j in 1:ncol(ode.final)) {
    hist(ode.final[, j], main = "", col = (j - 1) %% 8 + 1)
    title(colnames(ode.final)[j],
      adj = 0, line = 0, cex.main = 1.5
    )
  }

  if (loc == "RI") {
    title("Rhode Island", adj = 0, line = -0.5, cex.main = 3, outer = TRUE)
  } else if (loc == "MA") {
    title("Massachusetts", adj = 0, line = -0.5, cex.main = 3, outer = TRUE)
  } else if (loc == "PA") {
    title("Pennsylvania", adj = 0, line = -0.5, cex.main = 3, outer = TRUE)
  }

  dev.off()
  embed_fonts(plot.name, outfile = plot.name)
}






