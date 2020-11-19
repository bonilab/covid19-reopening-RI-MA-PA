### ms-figures.R
### Nathan Wikle
###
### Generates Figures 2, (part of) 4, and 5 for manuscript and supplementary 
###   material, using CSV output from state-level runs. CAUTION: some 
###   state-level params are hard-coded into 'fig2panel'. These were current
###   for Sept 6 runs, as of 12 Oct 2020. They may need to be updated for 
###   future runs.


fig2panel <- function(beta.directory, ode.directory, data.directory, odepath, 
                      loc, const, plot.name, alpha, subsample = NA, 
                      axis.size = 1.75, lab.size = 1.75,
                      ncol = 2, nrow = 4, height = 12, width = 12, png.true = F,
                      font.family = "Helvetica"){
  
  
  ### make samples data structure for traj.sim
  
  # grab samples from csv files
  betas <- read.csv(beta.directory)
  ode <- read.csv(ode.directory)
  
  # data
  df <- read.csv(data.directory)
  df <- data.frame(df)
  
  if (loc == "PA"){
    
    ode.inds <- 1:34
    
    ### clean up data (CHECK THIS IN FUTURE!)
    
    # cumulative value is wrong for age 80+
    #df[135:137,24] # ..., 11526, 13600, 11645, ...
    # change 13600 to 11600:
    df[136,24] <- 11600
    
    na.70 <- is.na(df[,41])
    df[na.70, 34:42] <- NA
    
    ## removing days before 61
    idx.remove <- which(df$daynum<61)
    if(length(idx.remove) > 0){
      df <- df[-idx.remove,]
    }
    
    ## removing days after 250
    idx.remove <- which(df$daynum>250)
    if(length(idx.remove) > 0){
      df <- df[-idx.remove,]
    }
    
    ## remove NAs
    df <- df[!is.na(df$daynum),]
    
    n.days <- nrow(df)
    days <- df$daynum
    
    end.day <- max(days, na.rm = TRUE) + 1
    num.days <- end.day - 60
    
    # cubic spline with one basis function every 7 days
    bspl <- create.bspline.basis(c(61, end.day), nbasis = round(num.days/7))
    Z <- eval.basis(bspl, 61:end.day)
    
    # iSpline for reporting rate
    knots <- c(84, 92, 100, 108, 130, 160, 190)
    Z.rr <- iSpline(x = 61:end.day, knots = knots, degree = 2)
    Z.rr <- cbind(rep(1, nrow(Z.rr)), Z.rr)
    Z.rr <- Z.rr[-nrow(Z.rr),c(1:3)]
    
    # introday
    id = 60
    
  } else if (loc == "RI"){
    
    ode.inds <- 1:36
    
    ## removing days before 61
    idx.remove=which(df$daynum<61)
    if(length(idx.remove)>0){
      df = df[-idx.remove,]
    }
    
    n.days <- nrow(df)
    days <- df$daynum
    
    ### mcmc initialization (don't change)
    end.day <- max(days, na.rm = TRUE) + 1
    num.days <- end.day - 60 ## number of days with data
    
    # cubic spline with one basis function every 7 days
    bspl <- create.bspline.basis(c(61, end.day), nbasis=round(num.days/7))
    Z <- eval.basis(bspl,61:end.day)
    
    # iSpline for reporting rate
    knots <- c(84, 92, 100, 108, 130, 160, 190)
    Z.rr <- iSpline(x = 61:end.day, knots = knots, degree = 2)
    Z.rr <- cbind(rep(1, nrow(Z.rr)), Z.rr)
    Z.rr <- Z.rr[-nrow(Z.rr),c(1,4:7)]
    
    # introday
    id <- 55
    
  } else if (loc == "MA") {
    
    ode.inds <- 1:35
    
    ## removing days before 61
    idx.remove <- which(df$daynum<61)
    if(length(idx.remove)>0){
      df <- df[-idx.remove,]
    }
    
    n.days <- nrow(df)
    days <- df$daynum
    
    ### mcmc initialization (don't change)
    end.day <- max(days, na.rm = TRUE) + 1
    num.days <- end.day - 60 ## number of days with data
    
    ### Create spline expansions
    
    # Cubic spline with one basis function every 7 days
    bspl <- create.bspline.basis(c(61, end.day), nbasis=round(num.days/7))
    Z <- eval.basis(bspl,61:end.day)
    
    # zero-order spline with one rate
    bspl.rr <- create.bspline.basis(c(61, end.day),nbasis=1,norder=1)
    Z.rr <- eval.basis(bspl.rr, 61:end.day)
    Z.rr <- Z.rr[-nrow(Z.rr),]
    
    # introday
    id <- 55
  }
  
  # list with required structures
  samples <- list(betas = as.matrix(betas),
                  ode = as.matrix(ode[,ode.inds]),
                  rr = as.matrix(ode[,-ode.inds]),
                  days = days,
                  Z.beta = Z,
                  Z.rr = Z.rr,
                  df = df,
                  loc = loc,
                  const = const,
                  introday = id,
                  sf.choice = FALSE)
  
  
  if (!is.na(subsample)){
    # remove some of the samples
    n.orig <- nrow(samples$betas)
    new.samples <- sample.int(n.orig, subsample, replace = FALSE)
    
    samples$betas <- samples$betas[new.samples,]
    samples$ode <- samples$ode[new.samples,]
    samples$rr <- samples$rr[new.samples,]
  }

  # reformat column names
  colnames(samples$ode) <- gsub('\\.', '-', colnames(samples$ode))
  colnames(samples$ode) <- gsub("hosp-report-rate", 
                                "hosp.report.rate", 
                                colnames(samples$ode))
  
  
  # process trajectories from the samples
  tp <- traj.sim(samples = samples,
                 odepath = odepath,
                 csv = T)
  # process actual data
  dp <- data.fine.process(samples$df, loc = samples$loc)
  
  # number of simulated output
  S <- nrow(tp$tot.sympt.new)

  # vector of dates
  dates <- seq(as.Date("2020/1/1"), by = "day", length.out = 365)

  # last day of simulated output
  end.day <- max(samples$days)

  # save plots as pdf
  if (png.true){
    png(filename = plot.name,
        family = font.family, 
        width = width, 
        height = height, 
        units = "in",
        res = 100)
  } else {
    pdf(file = plot.name, family = font.family, h = height, w = width)
  }
  # grid of plots
  par(mfrow = c(nrow,ncol))
  par(mar = c(3, 2, 1, 1), oma = c(2, 2, 0.5, 0.5))

  ##########################################
  ### A. new symptomatic cases
  ##########################################

  y.max <- max(max(dp$tot.sympt.new, na.rm = T),
               max(tp$tot.sympt.new))

  # plot total new symptomatic cases
  plot(dates[tp$days.odesim],
       tp$tot.sympt.new[1,], type = "l", col = addTrans("grey", alpha),
       ylim = c(0, y.max), ylab = "", xlab = "",
       cex.axis = axis.size, cex.lab = lab.size, bty = "l")
  #main = "New Symptomatic Cases (Reported)")
  title("New Symptomatic Cases (Reported)", adj = 0.95,
        line = -1.5, cex.main = 1.75)
  title("A", adj = 0.05, line = -1.5, cex.main = 2)

  for(j in 2:S){
    lines(dates[tp$days.odesim],
          tp$tot.sympt.new[j,],
          type = "l",
          col = addTrans("grey", alpha))
  }

  # add median and 95% quantiles
  mn <- apply(tp$tot.sympt.new, 2, median)
  uq <- apply(tp$tot.sympt.new, 2, quantile, 0.975)
  lq <- apply(tp$tot.sympt.new, 2, quantile, 0.025)

  lines(dates[tp$days.odesim], mn, type = "l", col = "blue")
  lines(dates[tp$days.odesim], uq, type = "l", lty = 3, col = "blue")
  lines(dates[tp$days.odesim], lq, type = "l", lty = 3, col = "blue")

  points(dates[dp$days], dp$tot.sympt.new, pch = 19, col = "black", cex = 0.5)

  #fig_label(LETTERS[1], cex=2)


  #################################################
  ### B. New Hospitalizations
  #################################################

  y.max <- max(max(dp$tot.hosp.new, na.rm = T),
               max(tp$tot.hosp.new))

  plot(dates[tp$days.odesim],
       tp$tot.hosp.new[1,], type = "l", col = addTrans("grey", alpha),
       ylim = c(0, y.max), ylab = "", xlab = "",
       cex.axis = axis.size, cex.lab = lab.size, bty = "l")
  title("New Hospitalizations", adj = 0.95,
        line = -1.5, cex.main = 1.75)
  title("B", adj = 0.05, line = -1.5, cex.main = 2)

  for(j in 2:S){
    lines(dates[tp$days.odesim],
          tp$tot.hosp.new[j,],
          type = "l",
          col = addTrans("grey", alpha))
  }

  # add median and 95% quantiles
  mn <- apply(tp$tot.hosp.new, 2, median)
  uq <- apply(tp$tot.hosp.new, 2, quantile, 0.975)
  lq <- apply(tp$tot.hosp.new, 2, quantile, 0.025)

  lines(dates[tp$days.odesim], mn, type = "l", col = "blue")
  lines(dates[tp$days.odesim], uq, type = "l", lty = 3, col = "blue")
  lines(dates[tp$days.odesim], lq, type = "l", lty = 3, col = "blue")

  if (samples$loc != "MA"){
    points(dates[dp$days], dp$tot.hosp.new, pch = 19, col = "black", cex = 0.5)
  }
  #fig_label(LETTERS[2], cex=2)

  ###################################
  ### C. Current Hospitalizations ###
  ###################################

  y.max <- max(max(dp$tot.hosp.curr, na.rm = T),
               max(tp$tot.hosp.curr))

  # plot current hospitalizations
  plot(dates[tp$days.odesim],
       tp$tot.hosp.curr[1,], type = "l", col = addTrans("grey", alpha),
       ylim = c(0, y.max), ylab = "", xlab = "",
       cex.axis = axis.size, cex.lab = lab.size, bty = "l")
  title("Current Hospitalizations", adj = 0.95,
        line = -1.5, cex.main = 1.75)
  title("C", adj = 0.05, line = -1.5, cex.main = 2)

  for(j in 2:S){
    lines(dates[tp$days.odesim],
          tp$tot.hosp.curr[j,],
          type = "l",
          col = addTrans("grey", alpha))
  }

  # add median and 95% quantiles
  mn <- apply(tp$tot.hosp.curr, 2, median)
  uq <- apply(tp$tot.hosp.curr, 2, quantile, 0.975)
  lq <- apply(tp$tot.hosp.curr, 2, quantile, 0.025)

  lines(dates[tp$days.odesim], mn, type = "l", col = "blue")
  lines(dates[tp$days.odesim], uq, type = "l", lty = 3, col = "blue")
  lines(dates[tp$days.odesim], lq, type = "l", lty = 3, col = "blue")

  points(dates[dp$days], dp$tot.hosp.curr, pch = 19, col = "black", cex = 0.5)

  #fig_label(LETTERS[3], cex=2)

  ###################################
  ### D. Current ICU Occupancy    ###
  ###################################

  y.max <- max(max(dp$tot.icu.curr, na.rm = T),
               max(tp$tot.icu.curr))

  # plot current ICU
  plot(dates[tp$days.odesim],
       tp$tot.icu.curr[1,], type = "l", col = addTrans("grey", alpha),
       ylim = c(0, y.max), ylab = "", xlab = "", yaxt = "none",
       cex.axis = axis.size, cex.lab = lab.size, bty = "l")
  axis(2, seq(0,120,20), cex.axis = axis.size, cex.lab = lab.size,
       labels = c("0", "", "40", "", "80", "", "120"))


  title("Current ICU Occupancy", adj = 0.95,
        line = -1.5, cex.main = 1.75)
  title("D", adj = 0.05, line = -1.5, cex.main = 2)

  for(j in 2:S){
    lines(dates[tp$days.odesim],
          tp$tot.icu.curr[j,],
          type = "l",
          col = addTrans("grey", alpha))
  }


  # add median and 95% quantiles
  mn <- apply(tp$tot.icu.curr, 2, median)
  uq <- apply(tp$tot.icu.curr, 2, quantile, 0.975)
  lq <- apply(tp$tot.icu.curr, 2, quantile, 0.025)

  lines(dates[tp$days.odesim], mn, type = "l", col = "blue")
  lines(dates[tp$days.odesim], uq, type = "l", lty = 3, col = "blue")
  lines(dates[tp$days.odesim], lq, type = "l", lty = 3, col = "blue")

  points(dates[dp$days], dp$tot.icu.curr, pch = 19, col = "black", cex = 0.5)

  #fig_label(LETTERS[4], cex=2)

  #######################################
  ### E. Current Ventilator Occupancy ###
  #######################################

  y.max <- max(max(dp$tot.vent.curr, na.rm = T),
               max(tp$tot.vent.curr))

  # plot current vent
  plot(dates[tp$days.odesim],
       tp$tot.vent.curr[1,], type = "l", col = addTrans("grey", alpha),
       ylim = c(0, y.max), ylab = "", xlab = "",
       cex.axis = axis.size, cex.lab = lab.size, bty = "l")
  title("Current Ventilator Occupancy",
        adj = 0.95, line = -1.5, cex.main = 1.75)
  title("E", adj = 0.05, line = -1.5, cex.main = 2)

  for(j in 2:S){
    lines(dates[tp$days.odesim],
          tp$tot.vent.curr[j,],
          type = "l",
          col = addTrans("grey", alpha))
  }


  # add median and 95% quantiles
  mn <- apply(tp$tot.vent.curr, 2, median)
  uq <- apply(tp$tot.vent.curr, 2, quantile, 0.975)
  lq <- apply(tp$tot.vent.curr, 2, quantile, 0.025)

  lines(dates[tp$days.odesim], mn, type = "l", col = "blue")
  lines(dates[tp$days.odesim], uq, type = "l", lty = 3, col = "blue")
  lines(dates[tp$days.odesim], lq, type = "l", lty = 3, col = "blue")

  points(dates[dp$days], dp$tot.vent.curr, pch = 19, col = "black", cex = 0.5)

  #  fig_label(LETTERS[5], cex=2)

  ###############################################
  ### F. Deaths
  ###############################################

  y.max <- max(max(dp$tot.deaths.new, na.rm = T),
               max(tp$tot.deaths.new))

  # plot new deaths
  plot(dates[tp$days.odesim],
       tp$tot.deaths.new[1,], type = "l", col = addTrans("grey", alpha),
       ylim = c(0, y.max), ylab = "", xlab = "",
       cex.axis = axis.size, cex.lab = lab.size, bty = "l")
  title("New Deaths", adj = 0.95, line = -1.5, cex.main = 1.75)
  title("F", adj = 0.05, line = -1.5, cex.main = 2)

  for(j in 2:S){
    lines(dates[tp$days.odesim],
          tp$tot.deaths.new[j,],
          type = "l",
          col = addTrans("grey", alpha))
  }


  # add median and 95% quantiles
  mn <- apply(tp$tot.deaths.new, 2, median)
  uq <- apply(tp$tot.deaths.new, 2, quantile, 0.975)
  lq <- apply(tp$tot.deaths.new, 2, quantile, 0.025)

  lines(dates[tp$days.odesim], mn, type = "l", col = "blue")
  lines(dates[tp$days.odesim], uq, type = "l", lty = 3, col = "blue")
  lines(dates[tp$days.odesim], lq, type = "l", lty = 3, col = "blue")

  points(dates[dp$days], dp$tot.deaths.new, pch = 19, col = "black", cex = 0.5)

  #  fig_label(LETTERS[6], cex=2)


  ###############################################
  ### G. Hospital Deaths
  ###############################################

  if (samples$loc == "MA"){
    y.max <- max(tp$tot.hosp.deaths.new)
  } else {
    y.max <- max(max(dp$tot.hosp.deaths.new, na.rm = T),
                 max(tp$tot.hosp.deaths.new))
  }

  # plot new deaths
  plot(dates[tp$days.odesim],
       tp$tot.hosp.deaths.new[1,], type = "l", col = addTrans("grey", alpha),
       ylim = c(0, y.max), ylab = "", xlab = "",
       cex.axis = axis.size, cex.lab = lab.size, bty = "l")
  title("New Hospital Deaths", adj = 0.95,
        line = -1.5, cex.main = 1.75)
  title("G", adj = 0.05, line = -1.5, cex.main = 2)

  for(j in 2:S){
    lines(dates[tp$days.odesim],
          tp$tot.hosp.deaths.new[j,],
          type = "l",
          col = addTrans("grey", alpha))
  }


  # add median and 95% quantiles
  mn <- apply(tp$tot.hosp.deaths.new, 2, median)
  uq <- apply(tp$tot.hosp.deaths.new, 2, quantile, 0.975)
  lq <- apply(tp$tot.hosp.deaths.new, 2, quantile, 0.025)

  lines(dates[tp$days.odesim], mn, type = "l", col = "blue")
  lines(dates[tp$days.odesim], uq, type = "l", lty = 3, col = "blue")
  lines(dates[tp$days.odesim], lq, type = "l", lty = 3, col = "blue")

  points(dates[dp$days], dp$tot.hosp.deaths.new, pch = 19, col = "black", cex = 0.5)

  #  fig_label(LETTERS[7], cex=2)


  ###############################################
  ### E. Discharges
  ###############################################

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
       ylim = c(0, y.max), ylab = "", xlab = "",
       cex.axis = axis.size, cex.lab = lab.size, bty = "l")
  title("New Hospital Discharges", adj = 0.95,
        line = -1.5, cex.main = 1.75)
  title("H", adj = 0.05, line = -1.5, cex.main = 2)

  for(j in 2:S){
    lines(dates[tp$days.odesim],
          tp$tot.hosp.discharges.new[j,],
          type = "l",
          col = addTrans("grey", alpha))
  }

  # add median and 95% quantiles
  mn <- apply(tp$tot.hosp.discharges.new, 2, median)
  uq <- apply(tp$tot.hosp.discharges.new, 2, quantile, 0.975)
  lq <- apply(tp$tot.hosp.discharges.new, 2, quantile, 0.025)

  lines(dates[tp$days.odesim], mn, type = "l", col = "blue")
  lines(dates[tp$days.odesim], uq, type = "l", lty = 3, col = "blue")
  lines(dates[tp$days.odesim], lq, type = "l", lty = 3, col = "blue")

  if (samples$loc != "MA"){
    points(dates[dp$days], dp$tot.hosp.discharges.new, pch = 19, col = "black", cex = 0.5)
  }

  #  fig_label(LETTERS[8], cex=2)

  dev.off()

  embed_fonts(plot.name, outfile = plot.name)
}

fig4 <- function(beta.directory, ode.directory, data.directory, odepath, 
                 loc, const, plot.name, subsample = NA, 
                 axis.size = 12, title.size = 14, 
                 legend.size = 12, legend.title.size = 14,
                 height = 4, width = 5,
                 font.family = "Helvetica"){
  
  ### make samples data structure for traj.sim
  
  # grab samples from csv files
  betas <- read.csv(beta.directory)
  ode <- read.csv(ode.directory)
  
  # data
  df <- read.csv(data.directory)
  df <- data.frame(df)
  
  if (loc == "PA"){
    
    ode.inds <- 1:34
    
    ### clean up data (CHECK THIS IN FUTURE!)
    
    # cumulative value is wrong for age 80+
    #df[135:137,24] # ..., 11526, 13600, 11645, ...
    # change 13600 to 11600:
    df[136,24] <- 11600
    
    na.70 <- is.na(df[,41])
    df[na.70, 34:42] <- NA
    
    ## removing days before 61
    idx.remove <- which(df$daynum<61)
    if(length(idx.remove) > 0){
      df <- df[-idx.remove,]
    }
    
    ## removing days after 250
    idx.remove <- which(df$daynum>250)
    if(length(idx.remove) > 0){
      df <- df[-idx.remove,]
    }
    
    ## remove NAs
    df <- df[!is.na(df$daynum),]
    
    n.days <- nrow(df)
    days <- df$daynum
    
    end.day <- max(days, na.rm = TRUE) + 1
    num.days <- end.day - 60
    
    # cubic spline with one basis function every 7 days
    bspl <- create.bspline.basis(c(61, end.day), nbasis = round(num.days/7))
    Z <- eval.basis(bspl, 61:end.day)
    
    # iSpline for reporting rate
    knots <- c(84, 92, 100, 108, 130, 160, 190)
    Z.rr <- iSpline(x = 61:end.day, knots = knots, degree = 2)
    Z.rr <- cbind(rep(1, nrow(Z.rr)), Z.rr)
    Z.rr <- Z.rr[-nrow(Z.rr),c(1:3)]
    
    # introday
    id = 60
    
  } else if (loc == "RI"){
    
    ode.inds <- 1:36
    
    ## removing days before 61
    idx.remove=which(df$daynum<61)
    if(length(idx.remove)>0){
      df = df[-idx.remove,]
    }
    
    n.days <- nrow(df)
    days <- df$daynum
    
    ### mcmc initialization (don't change)
    end.day <- max(days, na.rm = TRUE) + 1
    num.days <- end.day - 60 ## number of days with data
    
    # cubic spline with one basis function every 7 days
    bspl <- create.bspline.basis(c(61, end.day), nbasis=round(num.days/7))
    Z <- eval.basis(bspl,61:end.day)
    
    # iSpline for reporting rate
    knots <- c(84, 92, 100, 108, 130, 160, 190)
    Z.rr <- iSpline(x = 61:end.day, knots = knots, degree = 2)
    Z.rr <- cbind(rep(1, nrow(Z.rr)), Z.rr)
    Z.rr <- Z.rr[-nrow(Z.rr),c(1,4:7)]
    
    # introday
    id <- 55
    
  } else if (loc == "MA") {
    
    ode.inds <- 1:35
    
    ## removing days before 61
    idx.remove <- which(df$daynum<61)
    if(length(idx.remove)>0){
      df <- df[-idx.remove,]
    }
    
    n.days <- nrow(df)
    days <- df$daynum
    
    ### mcmc initialization (don't change)
    end.day <- max(days, na.rm = TRUE) + 1
    num.days <- end.day - 60 ## number of days with data
    
    ### Create spline expansions
    
    # Cubic spline with one basis function every 7 days
    bspl <- create.bspline.basis(c(61, end.day), nbasis=round(num.days/7))
    Z <- eval.basis(bspl,61:end.day)
    
    # zero-order spline with one rate
    bspl.rr <- create.bspline.basis(c(61, end.day),nbasis=1,norder=1)
    Z.rr <- eval.basis(bspl.rr, 61:end.day)
    Z.rr <- Z.rr[-nrow(Z.rr),]
    
    # introday
    id <- 55
  }
  
  # list with required structures
  samples <- list(betas = as.matrix(betas),
                  ode = as.matrix(ode[,ode.inds]),
                  rr = as.matrix(ode[,-ode.inds]),
                  days = days,
                  Z.beta = Z,
                  Z.rr = Z.rr,
                  df = df,
                  loc = loc,
                  const = const,
                  introday = id,
                  sf.choice = FALSE)
  
  # reformat column names
  colnames(samples$ode) <- gsub('\\.', '-', colnames(samples$ode))
  colnames(samples$ode) <- gsub("hosp-report-rate", 
                                "hosp.report.rate", 
                                colnames(samples$ode))
  
  ### plot contact rate parameters ###

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
      
      
  ### Plot Contact Rates ###
      
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
                       values=c("Lockdown" = "gray55", 
                                "Post-Lockdown" = "black")) + 
    theme(legend.justification=c(0.95,0.75), legend.position=c(0.95,0.75)) +
    theme(axis.text = element_text(size = axis.size),
          axis.title = element_text(size = title.size, face="bold")) + 
    theme(legend.text=element_text(size = legend.size),
          legend.title=element_text(size = legend.title.size, face="bold"),
          legend.background=element_blank()) + 
    theme(text = element_text(family="Ubuntu Mono")) 
    
    # save plot
    save_plot(plot.name,
              cr.plot, 
              base_height = height,
              base_width = width)
    
    # embed font in pdf
    embed_fonts(file = plot.name, outfile = plot.name)
}

fig5panel <- function(ri.beta, ri.ode, ri.data, ri.const,
                      pa.beta, pa.ode, pa.data, pa.const,
                      ma.beta, ma.ode, ma.data, ma.const,
                      odepath, 
                      plot.name, alpha, subsample = NA, 
                      axis.size = 1, lab.size = 1,
                      height = 16, width = 11, logplot = F,
                      font.family = "Helvetica"){
  
  ### make samples data structure for traj.sim
  
  ### Rhode Island ###
  
  # grab samples from csv files
  ri.betas <- read.csv(ri.beta)
  ri.ode <- read.csv(ri.ode)
  
  # data
  ri.df <- read.csv(ri.data)
  ri.df <- data.frame(ri.df)
  
  # ode index
  ri.ode.inds <- 1:36
  
  ## removing days before 61
  ri.idx.remove=which(ri.df$daynum<61)
  if(length(ri.idx.remove)>0){
    ri.df = ri.df[-ri.idx.remove,]
  }
  
  ri.n.days <- nrow(ri.df)
  ri.days <- ri.df$daynum
  
  ### mcmc initialization (don't change)
  ri.end.day <- max(ri.days, na.rm = TRUE) + 1
  ri.num.days <- ri.end.day - 60 ## number of days with data
  
  # cubic spline with one basis function every 7 days
  ri.bspl <- create.bspline.basis(c(61, ri.end.day), nbasis=round(ri.num.days / 7))
  ri.Z <- eval.basis(ri.bspl,61:ri.end.day)
  
  # iSpline for reporting rate
  ri.knots <- c(84, 92, 100, 108, 130, 160, 190)
  ri.Z.rr <- iSpline(x = 61:ri.end.day, knots = ri.knots, degree = 2)
  ri.Z.rr <- cbind(rep(1, nrow(ri.Z.rr)), ri.Z.rr)
  ri.Z.rr <- ri.Z.rr[-nrow(ri.Z.rr),c(1,4:7)]
  
  # introday
  ri.id = 55
  
  # list with required structures
  ri.samples <- list(betas = as.matrix(ri.betas),
                  ode = as.matrix(ri.ode[,ri.ode.inds]),
                  rr = as.matrix(ri.ode[,-ri.ode.inds]),
                  days = ri.days,
                  Z.beta = ri.Z,
                  Z.rr = ri.Z.rr,
                  df = ri.df,
                  loc = "RI",
                  const = ri.const,
                  introday = ri.id,
                  sf.choice = FALSE)
  
  # reformat column names
  colnames(ri.samples$ode) <- gsub('\\.', '-', colnames(ri.samples$ode))
  colnames(ri.samples$ode) <- gsub("hosp-report-rate", 
                                  "hosp.report.rate", 
                                  colnames(ri.samples$ode))
  
  
  # process trajectories from the samples
  ri.tp <- traj.sim(samples = ri.samples,
                 odepath = odepath,
                 csv = T)
  
  # process actual data
  ri.dp <- data.fine.process(ri.samples$df, loc = ri.samples$loc)
  
  ### Massachusetts ###
  
  # grab samples from csv files
  ma.betas <- read.csv(ma.beta)
  ma.ode <- read.csv(ma.ode)
  
  # data
  ma.df <- read.csv(ma.data)
  ma.df <- data.frame(ma.df)
  
  # ode index
  ma.ode.inds <- 1:35
  
  ## removing days before 61
  ma.idx.remove=which(ma.df$daynum<61)
  if(length(ma.idx.remove)>0){
    ma.df = ma.df[-ma.idx.remove,]
  }
  
  ma.n.days <- nrow(ma.df)
  ma.days <- ma.df$daynum
  
  ### mcmc initialization (don't change)
  ma.end.day <- max(ma.days, na.rm = TRUE) + 1
  ma.num.days <- ma.end.day - 60 ## number of days with data
  
  # cubic spline with one basis function every 7 days
  ma.bspl <- create.bspline.basis(c(61, ma.end.day), 
                                  nbasis = round(ma.num.days / 7))
  ma.Z <- eval.basis(ma.bspl, 61:ma.end.day)
  
  # zero-order spline with one rate
  ma.bspl.rr <- create.bspline.basis(c(61, ma.end.day), nbasis=1, norder=1)
  ma.Z.rr <- eval.basis(ma.bspl.rr, 61:ma.end.day)
  ma.Z.rr <- ma.Z.rr[-nrow(ma.Z.rr),]
  ma.Z.rr <- as.vector(ma.Z.rr)
  
  # introday
  ma.id = 55
  
  # list with required structures
  ma.samples <- list(betas = as.matrix(ma.betas),
                     ode = as.matrix(ma.ode[,ma.ode.inds]),
                     rr = as.matrix(ma.ode[,-ma.ode.inds]),
                     days = ma.days,
                     Z.beta = ma.Z,
                     Z.rr = ma.Z.rr,
                     df = ma.df,
                     loc = "MA",
                     const = ma.const,
                     introday = ma.id,
                     sf.choice = FALSE)
  
  # reformat column names
  colnames(ma.samples$ode) <- gsub('\\.', '-', colnames(ma.samples$ode))
  colnames(ma.samples$ode) <- gsub("hosp-report-rate", 
                                   "hosp.report.rate", 
                                   colnames(ma.samples$ode))
  
  
  # process trajectories from the samples
  ma.tp <- traj.sim(samples = ma.samples,
                    odepath = odepath,
                    csv = T)
  
  # process actual data
  ma.dp <- data.fine.process(ma.samples$df, loc = ma.samples$loc)
  
  ### Pennsylvania ###
  
  # ode index
  pa.ode.inds <- 1:34
  
  # grab samples from csv files
  pa.betas <- read.csv(pa.beta)
  pa.ode <- read.csv(pa.ode)
  
  # data
  pa.df <- read.csv(pa.data)
  pa.df <- data.frame(pa.df)
  
  # cumulative value is wrong for age 80+
  #   (change 13600 to 11600)
  pa.df[136,24] <- 11600
  
  pa.na.70 <- is.na(pa.df[,41])
  pa.df[pa.na.70, 34:42] <- NA
  
  ## removing days before 61
  pa.idx.remove <- which(pa.df$daynum<61)
  if(length(pa.idx.remove) > 0){
    pa.df <- pa.df[-pa.idx.remove,]
  }
  
  ## removing days after 250
  pa.idx.remove <- which(pa.df$daynum>250)
  if(length(pa.idx.remove) > 0){
    pa.df <- pa.df[-pa.idx.remove,]
  }
  
  ## remove NAs
  pa.df <- pa.df[!is.na(pa.df$daynum),]

  pa.n.days <- nrow(pa.df)
  pa.days <- pa.df$daynum
  
  ### mcmc initialization (don't change)
  pa.end.day <- max(pa.days, na.rm = TRUE) + 1
  pa.num.days <- pa.end.day - 60 ## number of days with data
  
  # cubic spline with one basis function every 7 days
  pa.bspl <- create.bspline.basis(c(61, pa.end.day), 
                                  nbasis = round(pa.num.days / 7))
  pa.Z <- eval.basis(pa.bspl, 61:ma.end.day)
  

  # iSpline for reporting rate
  pa.knots <- c(84, 92, 100, 108, 130, 160, 190)
  pa.Z.rr <- iSpline(x = 61:pa.end.day, knots = pa.knots, degree = 2)
  pa.Z.rr <- cbind(rep(1, nrow(pa.Z.rr)), pa.Z.rr)
  pa.Z.rr <- pa.Z.rr[-nrow(pa.Z.rr),c(1:3)]
  
  # introday
  pa.id = 60
  
  # list with required structures
  pa.samples <- list(betas = as.matrix(pa.betas),
                     ode = as.matrix(pa.ode[,pa.ode.inds]),
                     rr = as.matrix(pa.ode[,-pa.ode.inds]),
                     days = pa.days,
                     Z.beta = pa.Z,
                     Z.rr = pa.Z.rr,
                     df = pa.df,
                     loc = "PA",
                     const = pa.const,
                     introday = pa.id,
                     sf.choice = FALSE)
  
  # reformat column names
  colnames(pa.samples$ode) <- gsub('\\.', '-', colnames(pa.samples$ode))
  colnames(pa.samples$ode) <- gsub("hosp-report-rate", 
                                   "hosp.report.rate", 
                                   colnames(pa.samples$ode))
  
  
  # process trajectories from the samples
  pa.tp <- traj.sim(samples = pa.samples,
                    odepath = odepath,
                    csv = T)
  
  # process actual data
  pa.dp <- data.fine.process(pa.samples$df, loc = pa.samples$loc)
  
  ri.color = "#733381"
  ma.color = "#F58024"
  pa.color = "#589C48"
  
  ############################################################################
  ### PLOT FIGURE FIVE
  ############################################################################
  
  pdf(file = plot.name, family = font.family, h = height, w = width)
  
  ### left, right, bottom, top
  
  m0 <- matrix(c(0, 1, 0, 0.97), nrow = 1, ncol = 4, byrow = T)
  
  m1 <- rbind(c(0.01, 0.99, 0.84, 0.99),
              c(0.01, 0.99, 0.67, 0.82),
              c(0.01, 0.99, 0.50, 0.65),
              c(0.01, 0.99, 0.33, 0.48),
              c(0.01, 0.99, 0.01, 0.31))
  
  m2 <- rbind(c(0.01, 0.33, 0.02, 0.98),
              c(0.33, 0.66, 0.02, 0.98),
              c(0.66, 1, 0.02, 0.98))
  
  m3 <- matrix(c(0.01, 0.99, 0.01, 0.99), nrow = 1, ncol = 4, byrow = T)
  
  
  split.screen(m0)
  split.screen(m1, screen = 1)
  split.screen(m2, screen = 2)
  split.screen(m2, screen = 3)
  split.screen(m2, screen = 4)
  split.screen(m2, screen = 5)
  split.screen(m3, screen = 6)
  
  
  screen(1)
  par(mar = c(0, 0, 0, 0),
      oma = c(0, 0, 2, 0))
  
  title("Rhode Island", adj = 0.15,
        line = 0, cex.main = 1.55, outer = T)
  
  title("Massachusetts", adj = 0.53,
        line = 0, cex.main = 1.55, outer = T)
  
  title("Pennsylvania", adj = 0.9,
        line = 0, cex.main = 1.55, outer = T)
  
  ############################################################################
  ### A. REPORTING RATE 
  ############################################################################
  
  screen(2)
  # bottom, left, top, right
  par(mar = c(0, 0, 2, 0), 
      oma = c(0, 0, 0, 0))
  
  title("A. Symptomatic Reporting Rate Over Time", adj = 0,
        line = 1, cex.main = 1.25, outer = F)
  
  #############################
  ###  A1: RI Reporting Rate
  #############################
  screen(7)
  
  # bottom, left, top, right
  par(mar = c(1, 2, 2, 0), 
      oma = c(0, 0, 0, 0))
  
  # some constants
  n.orig <- nrow(ri.samples$betas)
  new.samples <- sample.int(n.orig, 150, replace = FALSE)
  ri.rr <- ri.tp$rr.full[new.samples,]
  
  # vector of dates
  dates <- seq(as.Date("2020/1/1"), by = "day", length.out = 365)
  
  # last day of simulated output
  ri.end.day <- max(ri.samples$days)
  
  
  # plot reporting rate over time
  plot(dates[c(ri.end.day - ncol(ri.rr) + 1):ri.end.day],
       ri.rr[1,], type = "l", col = addTrans("grey", alpha),
       ylim = c(0, 1), ylab = "", xlab = "",
       cex.axis = axis.size, cex.lab = lab.size, bty = "l")
  
  
  for(j in 2:nrow(ri.rr)){
    lines(dates[c(ri.end.day - ncol(ri.rr) + 1):ri.end.day],
          ri.rr[j,],
          type = "l",
          col = addTrans("grey", alpha))
  }
  
  # add median and 95% quantiles
  mn <- apply(ri.rr, 2, median)
  uq <- apply(ri.rr, 2, quantile, 0.975)
  lq <- apply(ri.rr, 2, quantile, 0.025)
  
  lines(dates[c(ri.end.day - ncol(ri.rr) + 1):ri.end.day], 
        mn, type = "l", col = ri.color, lwd = 2)
  lines(dates[c(ri.end.day - ncol(ri.rr) + 1):ri.end.day], 
        uq, type = "l", lty = 3, col = ri.color, lwd = 2)
  lines(dates[c(ri.end.day - ncol(ri.rr) + 1):ri.end.day], 
        lq, type = "l", lty = 3, col = ri.color, lwd = 2)
  
  #############################
  ###  A2: MA Reporting Rate
  #############################
  
  screen(8)
  # bottom, left, top, right
  par(mar = c(1, 2, 2, 0), 
      oma = c(0, 0, 0, 0))
  
  # vector of dates
  dates <- seq(as.Date("2020/1/1"), by = "day", length.out = 365)
  
  ma.rr <- ma.tp$rr.full[new.samples,]
  
  # last day of simulated output
  ma.end.day <- max(ma.samples$days)
  
  
  # plot reporting rate over time
  plot(dates[c(ma.end.day - ncol(ma.rr) + 1):ma.end.day],
       ma.rr[1,], type = "l", col = addTrans("grey", alpha),
       ylim = c(0, 1), ylab = "", xlab = "",
       cex.axis = axis.size, cex.lab = lab.size, bty = "l")
  
  for(j in 2:nrow(ma.rr)){
    lines(dates[c(ma.end.day - ncol(ma.rr) + 1):ma.end.day],
          ma.rr[j,],
          type = "l",
          col = addTrans("grey", alpha))
  }
  
  # add median and 95% quantiles
  mn <- apply(ma.rr, 2, median)
  uq <- apply(ma.rr, 2, quantile, 0.975)
  lq <- apply(ma.rr, 2, quantile, 0.025)
  
  lines(dates[c(ma.end.day - ncol(ma.rr) + 1):ma.end.day], 
        mn, type = "l", col = ma.color, lwd = 2)
  lines(dates[c(ma.end.day - ncol(ma.rr) + 1):ma.end.day], 
        uq, type = "l", lty = 3, col = ma.color, lwd = 2)
  lines(dates[c(ma.end.day - ncol(ma.rr) + 1):ma.end.day], 
        lq, type = "l", lty = 3, col = ma.color, lwd = 2)
  
  #############################
  ###  A3: PA Reporting Rate
  #############################
  screen(9)
  # bottom, left, top, right
  par(mar = c(1, 2, 2, 0), 
      oma = c(0, 0, 0, 0))
  
  # vector of dates
  dates <- seq(as.Date("2020/1/1"), by = "day", length.out = 365)
  
  pa.rr <- pa.tp$rr.full[new.samples,]
  
  # last day of simulated output
  pa.end.day <- max(pa.samples$days)
  
  # plot reporting rate over time
  plot(dates[c(pa.end.day - ncol(pa.rr) + 1):pa.end.day],
       pa.rr[1,], type = "l", col = addTrans("grey", alpha),
       ylim = c(0, 1), ylab = "", xlab = "",
       cex.axis = axis.size, cex.lab = lab.size, bty = "l")
  
  for(j in 2:nrow(pa.rr)){
    lines(dates[c(pa.end.day - ncol(pa.rr) + 1):pa.end.day],
          pa.rr[j,],
          type = "l",
          col = addTrans("grey", alpha))
  }
  
  # add median and 95% quantiles
  mn <- apply(pa.rr, 2, median)
  uq <- apply(pa.rr, 2, quantile, 0.975)
  lq <- apply(pa.rr, 2, quantile, 0.025)
  
  lines(dates[c(pa.end.day - ncol(pa.rr) + 1):pa.end.day], 
        mn, type = "l", col = pa.color, lwd = 2)
  lines(dates[c(pa.end.day - ncol(pa.rr) + 1):pa.end.day], 
        uq, type = "l", lty = 3, col = pa.color, lwd = 2)
  lines(dates[c(pa.end.day - ncol(pa.rr) + 1):pa.end.day], 
        lq, type = "l", lty = 3, col = pa.color, lwd = 2)
  
  
  ############################################################################
  ### B. HOSPITAL STAY
  ############################################################################
  
  screen(3)
  # bottom, left, top, right
  par(mar = c(0, 0, 2, 0), 
      oma = c(0, 0, 0, 0))
  
  title("B. Length of Hospital Stay (Days)", adj = 0,
        line = 1, cex.main = 1.25, outer = F)
  
  # densities for each state
  ri.hospstay.ind <- which(colnames(ri.samples$ode) == "dev-len-hospstay")
  ri.hospstay.dens <- density(10.7 * ri.samples$ode[,ri.hospstay.ind])
  
  pa.hospstay.ind <- which(colnames(pa.samples$ode) == "dev-len-hospstay")
  pa.hospstay.dens <- density(10.7 * pa.samples$ode[,pa.hospstay.ind])
  
  ma.hospstay.ind <- which(colnames(ma.samples$ode) == "dev-len-hospstay")
  ma.hospstay.dens <- density(10.7 * ma.samples$ode[, ma.hospstay.ind])
  
  # plot limits
  max.x <- max(max(ri.hospstay.dens$x),
               max(pa.hospstay.dens$x), 
               max(ma.hospstay.dens$x))
  
  min.x <- min(min(ri.hospstay.dens$x),
               min(pa.hospstay.dens$x), 
               min(ma.hospstay.dens$x))
  
  #############################
  ###  B1: RI Hospital Stay
  #############################
  screen(10)
  # bottom, left, top, right
  par(mar = c(1, 2, 2, 0), 
      oma = c(0, 0, 0, 0))
  
  plot(ri.hospstay.dens, xlab = "", ylab = "", main = "", bty = "l",
       xlim = c(min.x, max.x),
       cex.axis = axis.size, cex.lab = lab.size)
  
  polygon(ri.hospstay.dens, 
          col = addTrans(ri.color, 200), 
          border=ri.color,
          lwd = 2) 
  
  
  #############################
  ###  B2: MA Hospital Stay
  #############################
  screen(11)
  # bottom, left, top, right
  par(mar = c(1, 2, 2, 0), 
      oma = c(0, 0, 0, 0))
  
  plot(ma.hospstay.dens, xlab = "", ylab = "", main = "", bty = "l",
       xlim = c(min.x, max.x),
       cex.axis = axis.size, cex.lab = lab.size,
       col = addTrans(ma.color, 200))
  
  polygon(ma.hospstay.dens, 
          col = addTrans(ma.color, 200), 
          border = addTrans(ma.color, 200),
          lwd = 2) 
  
  #############################
  ###  B3: PA Hospital Stay
  #############################
  screen(12)
  # bottom, left, top, right
  par(mar = c(1, 2, 2, 0), 
      oma = c(0, 0, 0, 0))
  
  plot(pa.hospstay.dens, xlab = "", ylab = "", main = "", bty = "l",
       xlim = c(min.x, max.x),
       cex.axis = axis.size, cex.lab = lab.size,
       col = addTrans(pa.color, 200))
  
  polygon(pa.hospstay.dens, 
          col = addTrans(pa.color, 200), 
          border = addTrans(pa.color, 200),
          lwd = 2) 
  
  
  
  ############################################################################
  ### C. PROBABILITY OF DYING AT HOME
  ############################################################################
  
  screen(4)
  # bottom, left, top, right
  par(mar = c(0, 0, 2, 0), 
      oma = c(0, 0, 0, 0))
  
  title("C. Probability of Dying at Home", adj = 0,
        line = 1, cex.main = 1.25, outer = F)
  
  # home 60 (RI only)
  ri.home60.ind <- which(colnames(ri.samples$ode) == "death-prob-home-60")
  ri.home60.dens <- density(ri.samples$ode[,ri.home60.ind])
  
  # home 70
  ri.home70.ind <- which(colnames(ri.samples$ode) == "death-prob-home-70")
  ri.home70.dens <- density(ri.samples$ode[,ri.home70.ind])
  
  ma.home70.ind <- which(colnames(ma.samples$ode) == "death-prob-home-70")
  ma.home70.dens <- density(ma.samples$ode[,ma.home70.ind])
  
  pa.home70.ind <- which(colnames(pa.samples$ode) == "death-prob-home-70")
  pa.home70.dens <- density(pa.samples$ode[,pa.home70.ind])
  
  # home 80
  ri.home80.ind <- which(colnames(ri.samples$ode) == "death-prob-home-80")
  ri.home80.dens <- density(ri.samples$ode[,ri.home80.ind])
  
  ma.home80.ind <- which(colnames(ma.samples$ode) == "death-prob-home-80")
  ma.home80.dens <- density(ma.samples$ode[,ma.home80.ind])
  
  pa.home80.ind <- which(colnames(pa.samples$ode) == "death-prob-home-80")
  pa.home80.dens <- density(pa.samples$ode[,pa.home80.ind])
  
  
  # plot limits
  max.x <- max(max(ri.home60.dens$x),
               max(ri.home70.dens$x),
               max(ma.home70.dens$x),
               max(pa.home70.dens$x),
               max(ri.home80.dens$x),
               max(ma.home80.dens$x),
               max(pa.home80.dens$x))
  
  min.x <- min(min(ri.home60.dens$x),
               min(ri.home70.dens$x),
               min(ma.home70.dens$x),
               min(pa.home70.dens$x),
               min(ri.home80.dens$x),
               min(ma.home80.dens$x),
               min(pa.home80.dens$x))
  
  ####################################
  ###  C1: RI Prob of Dying at Home 
  ####################################
  
  screen(13)
  # bottom, left, top, right
  par(mar = c(1, 2, 2, 0), 
      oma = c(0, 0, 0, 0))
  
  plot(ri.home60.dens, xlab = "", ylab = "", main = "", bty = "l",
       xlim = c(min.x, max.x),
       cex.axis = axis.size, cex.lab = lab.size,
       col = addTrans(ri.color, 25))
  
  polygon(ri.home60.dens,
          col = addTrans(ri.color, 25),
          border = addTrans(ri.color, 25),
          lwd = 2)
  
  lines(ri.home70.dens, xlab = "", ylab = "", main = "", bty = "l",
        xlim = c(min.x, max.x),
        cex.axis = axis.size, cex.lab = lab.size,
        col = addTrans(ri.color, 100))
  
  polygon(ri.home70.dens,
          col = addTrans(ri.color, 100),
          border = addTrans(ri.color, 100),
          lwd = 2)
  
  lines(ri.home80.dens, xlab = "", ylab = "", main = "", bty = "l",
        xlim = c(min.x, max.x),
        cex.axis = axis.size, cex.lab = lab.size,
        col = addTrans(ri.color, 250))
  
  polygon(ri.home80.dens,
          col = addTrans(ri.color, 250),
          border = addTrans(ri.color, 250),
          lwd = 2)
  
  
  legend("topright", inset = 0.1,
         legend=c("60-69 y/o", "70-79 y/o", "80+ y/o"),
         col = c(addTrans(ri.color, 25),
                 addTrans(ri.color, 100),
                 addTrans(ri.color, 250)),
         lty=rep(1,3), cex=1,
         box.lty = 0, lwd = rep(3,3))
  
  
  ####################################
  ###  C2: MA Prob of Dying at Home 
  ####################################
  screen(14)
  # bottom, left, top, right
  par(mar = c(1, 2, 2, 0), 
      oma = c(0, 0, 0, 0))
  
  plot(ma.home70.dens, xlab = "", ylab = "", main = "", bty = "l",
       xlim = c(min.x, max.x),
       cex.axis = axis.size, cex.lab = lab.size,
       col = addTrans(ma.color, 100))
  
  polygon(ma.home70.dens,
          col = addTrans(ma.color, 100),
          border = addTrans(ma.color, 100),
          lwd = 2)
  
  lines(ma.home80.dens, xlab = "", ylab = "", main = "", bty = "l",
        xlim = c(min.x, max.x),
        cex.axis = axis.size, cex.lab = lab.size,
        col = addTrans(ma.color, 250))
  
  polygon(ma.home80.dens,
          col = addTrans(ma.color, 250),
          border = addTrans(ma.color, 250),
          lwd = 2)
  
  legend("topright", inset = 0.1,
         legend=c("70-79 y/o", "80+ y/o"),
         col = c(addTrans(ma.color, 100),
                 addTrans(ma.color, 250)),
         lty=rep(1,2), cex=1,
         box.lty = 0, lwd = rep(3,2))
  
  
  ####################################
  ###  C3: PA Prob of Dying at Home 
  ####################################
  screen(15)
  # bottom, left, top, right
  par(mar = c(1, 2, 2, 0), 
      oma = c(0, 0, 0, 0))
  
  plot(pa.home70.dens, xlab = "", ylab = "", main = "", bty = "l",
       xlim = c(min.x, max.x),
       cex.axis = axis.size, cex.lab = lab.size,
       col = addTrans(pa.color, 100))
  
  polygon(pa.home70.dens,
          col = addTrans(pa.color, 100),
          border = addTrans(pa.color, 100),
          lwd = 2)
  
  lines(pa.home80.dens, xlab = "", ylab = "", main = "", bty = "l",
        xlim = c(min.x, max.x),
        cex.axis = axis.size, cex.lab = lab.size,
        col = addTrans(pa.color, 250))
  
  polygon(pa.home80.dens,
          col = addTrans(pa.color, 250),
          border = addTrans(pa.color, 250),
          lwd = 2)
  
  legend("topright", inset = 0.1,
         legend=c("70-79 y/o", "80+ y/o"),
         col = c(addTrans(pa.color, 100),
                 addTrans(pa.color, 250)),
         lty=rep(1,2), cex=1,
         box.lty = 0, lwd = rep(3,2))
  
  
  ############################################################################
  ### D. ICU ADMISSION PROBABILITY
  ############################################################################
  
  screen(5)
  # bottom, left, top, right
  par(mar = c(0, 0, 2, 0), 
      oma = c(0, 0, 0, 0))
  
  title("D. ICU Admission Probability", adj = 0,
        line = 1, cex.main = 1.25, outer = F)
  
  # lockdown
  ri.icu1.ind <- which(colnames(ri.samples$ode) == "dev-icu-frac")
  ri.icu1.dens <- density(icuTransformation(ri.samples$ode[,ri.icu1.ind], "RI"))
  
  ma.icu1.ind <- which(colnames(ma.samples$ode) == "dev-icu-frac")
  ma.icu1.dens <- density(icuTransformation(ma.samples$ode[,ma.icu1.ind], "MA"))
  
  pa.icu1.ind <- which(colnames(pa.samples$ode) == "dev-icu-frac")
  pa.icu1.dens <- density(icuTransformation(pa.samples$ode[,pa.icu1.ind], "PA"))
  
  # post-lockdown
  ri.icu2.ind <- which(colnames(ri.samples$ode) == "dev-icu-frac-phase2")
  ri.icu2.dens <- density(icuTransformation(ri.samples$ode[,ri.icu2.ind], "RI"))
  
  ma.icu2.ind <- which(colnames(ma.samples$ode) == "dev-icu-frac-phase2")
  ma.icu2.dens <- density(icuTransformation(ma.samples$ode[,ma.icu2.ind], "MA"))
  
  pa.icu2.ind <- which(colnames(pa.samples$ode) == "dev-icu-frac-phase2")
  pa.icu2.dens <- density(icuTransformation(pa.samples$ode[,pa.icu2.ind], "PA"))
  
  # plot limits
  max.x <- max(max(ri.icu1.dens$x),
               max(ma.icu1.dens$x),
               max(pa.icu1.dens$x),
               max(ri.icu2.dens$x),
               max(ma.icu2.dens$x),
               max(pa.icu2.dens$x))
  
  min.x <- min(min(ri.icu1.dens$x),
               min(ma.icu1.dens$x),
               min(pa.icu1.dens$x),
               min(ri.icu2.dens$x),
               min(ma.icu2.dens$x),
               min(pa.icu2.dens$x))
  
  
  #################################
  ###  D1: RI ICU admission prob 
  #################################
  screen(16)
  # bottom, left, top, right
  par(mar = c(1, 2, 2, 0), 
      oma = c(0, 0, 0, 0))
  
  plot(ri.icu1.dens, xlab = "", ylab = "", main = "", bty = "l",
       xlim = c(min.x, max.x),
       ylim = c(0, max(max(ri.icu1.dens$y), max(ri.icu2.dens$y))),
       cex.axis = axis.size, cex.lab = lab.size,
       col = addTrans(ri.color, 100))
  
  polygon(ri.icu1.dens, 
          col = addTrans(ri.color, 100), 
          border = addTrans(ri.color, 100),
          lwd = 2)
  
  lines(ri.icu2.dens, 
        xlim = c(min.x, max.x),
        ylim = c(0, max(max(ri.icu1.dens$y), max(ri.icu2.dens$y))),
        cex.axis = axis.size, cex.lab = lab.size,
        col = addTrans(ri.color, 250))
  
  polygon(ri.icu2.dens, 
          col = addTrans(ri.color, 250), 
          border = addTrans(ri.color, 250),
          lwd = 2)
  
  legend("topright", inset = 0.05,
         legend=c("during lockdown", "post-lockdown"),
         col = c(addTrans(ri.color, 100),
                 addTrans(ri.color, 250)),
         lty=rep(1,2), cex=1,
         box.lty = 0, lwd = rep(3,2))  
  
  #################################
  ###  D2: MA ICU admission prob 
  #################################
  screen(17)
  # bottom, left, top, right
  par(mar = c(1, 2, 2, 0), 
      oma = c(0, 0, 0, 0))
  
  plot(ma.icu1.dens, xlab = "", ylab = "", main = "", bty = "l",
       xlim = c(min.x, max.x),
       ylim = c(0, max(max(ma.icu1.dens$y), max(ma.icu2.dens$y))),
       cex.axis = axis.size, cex.lab = lab.size,
       col = addTrans(ma.color, 100))
  
  polygon(ma.icu1.dens, 
          col = addTrans(ma.color, 100), 
          border = addTrans(ma.color, 100),
          lwd = 2)
  
  lines(ma.icu2.dens, 
        xlim = c(min.x, max.x),
        ylim = c(0, max(max(ma.icu1.dens$y), max(ma.icu2.dens$y))),
        cex.axis = axis.size, cex.lab = lab.size,
        col = addTrans(ma.color, 250))
  
  polygon(ma.icu2.dens, 
          col = addTrans(ma.color, 250), 
          border = addTrans(ma.color, 250),
          lwd = 2)
  
  legend("topleft", inset = 0.05,
         legend=c("during lockdown", "post-lockdown"),
         col = c(addTrans(ma.color, 100),
                 addTrans(ma.color, 250)),
         lty=rep(1,2), cex=1,
         box.lty = 0, lwd = rep(3,2))  
  
  #################################
  ###  D3: PA ICU admission prob 
  #################################
  screen(18)
  # bottom, left, top, right
  par(mar = c(1, 2, 2, 0), 
      oma = c(0, 0, 0, 0))
  
  plot(pa.icu1.dens, xlab = "", ylab = "", main = "", bty = "l",
       xlim = c(min.x, max.x),
       ylim = c(0, max(max(pa.icu1.dens$y), max(pa.icu2.dens$y))),
       cex.axis = axis.size, cex.lab = lab.size,
       col = addTrans(ma.color, 100))
  
  polygon(pa.icu1.dens, 
          col = addTrans(pa.color, 100), 
          border = addTrans(pa.color, 100),
          lwd = 2)
  
  lines(pa.icu2.dens, 
        xlim = c(min.x, max.x),
        ylim = c(0, max(max(pa.icu1.dens$y), max(pa.icu2.dens$y))),
        cex.axis = axis.size, cex.lab = lab.size,
        col = addTrans(pa.color, 250))
  
  polygon(pa.icu2.dens, 
          col = addTrans(pa.color, 250), 
          border = addTrans(pa.color, 250),
          lwd = 2)
  
  legend("topright", inset = 0.05,
         legend=c("during lockdown", "post-lockdown"),
         col = c(addTrans(pa.color, 100),
                 addTrans(pa.color, 250)),
         lty=rep(1,2), cex=1,
         box.lty = 0, lwd = rep(3,2))  
  
  
  
  ############################################################################
  ### E. Fraction of Symptomatic Individuals Hopsitalized, By Age
  ############################################################################
  
  screen(6)  ### bottom screen
  # bottom, left, top, right
  par(mar = c(0, 0, 2, 0), 
      oma = c(0, 0, 0, 0))
  
  title("E. Fraction of Symptomatic Individuals Hospitalized, By Age Group", adj = 0,
        line = 1, cex.main = 1.25, outer = F)
  
  ### Plot with Credible Intervals ###
  screen(19)
  
  if (logplot){
    # bottom, left, top, right
    par(mar = c(4, 3, 2, 0), 
        oma = c(0, 0, 0, 0))
  } else {
    # bottom, left, top, right
    par(mar = c(2, 3, 2, 0), 
        oma = c(0, 0, 0, 0))
  }
  
  
  ri.h10.ind <- which(colnames(ri.samples$ode) == "hosp-frac-10")               
  ri.hfrac.inds <- ri.h10.ind + 0:7   
  ri.hfrac <- ri.samples$ode[,ri.hfrac.inds]
  ri.hfrac.sum <- t(apply(ri.hfrac, 2, quantile, prob = c(0.5, 0.025, 0.975)))
  
  ma.h10.ind <- which(colnames(ma.samples$ode) == "hosp-frac-10")               
  ma.hfrac.inds <- ma.h10.ind + 0:7   
  ma.hfrac <- ma.samples$ode[,ma.hfrac.inds]
  ma.hfrac.sum <- t(apply(ma.hfrac, 2, quantile, prob = c(0.5, 0.025, 0.975)))
  
  pa.h10.ind <- which(colnames(pa.samples$ode) == "hosp-frac-10")               
  pa.hfrac.inds <- pa.h10.ind + 0:7   
  pa.hfrac <- pa.samples$ode[,pa.hfrac.inds]
  pa.hfrac.sum <- t(apply(pa.hfrac, 2, quantile, prob = c(0.5, 0.025, 0.975)))
  
  
  if (logplot){
    ri.hfrac.sum <- log2(ri.hfrac.sum)
    ma.hfrac.sum <- log2(ma.hfrac.sum)
    pa.hfrac.sum <- log2(pa.hfrac.sum)
    
    plot(y = 1:8 + 0.1, x = ri.hfrac.sum[,1], ylim = c(0.5, 8.5), xlim = c(-6.5,-1.5),
         xlab = expression(paste("fraction of symptomatic individualized who are hospitalized (", log[2], ")", sep = "")), 
         ylab = "", main = "", col = addTrans(ri.color, 250),
         pch = 16, bty = "l", yaxt = "n", cex = 1.25,
         cex.axis = axis.size, cex.lab = lab.size)
    
    axis(2, at=1:8, labels=c("10-19", "20-29", "30-39", "40-49", "50-59",
                             "60-69", "70-79", "80+"),
         cex.axis = axis.size, cex.lab = lab.size,
         las = 1)
    
    abline(h = 1, col = addTrans("grey", 50), lwd = 2)
    abline(h = 2, col = addTrans("grey", 50), lwd = 2)
    abline(h = 3, col = addTrans("grey", 50), lwd = 2)
    abline(h = 4, col = addTrans("grey", 50), lwd = 2)
    abline(h = 5, col = addTrans("grey", 50), lwd = 2)
    abline(h = 6, col = addTrans("grey", 50), lwd = 2)
    abline(h = 7, col = addTrans("grey", 50), lwd = 2)
    abline(h = 8, col = addTrans("grey", 50), lwd = 2)
    
    segments(y0 = 1:8 + 0.1, x0 = ri.hfrac.sum[,2], 
             y1 = 1:8 + 0.1, x1 = ri.hfrac.sum[,3],
             col = addTrans(ri.color, 100), 
             lty = 1, lwd = 3)
    
    points(y = 1:8 + 0.1, x = ri.hfrac.sum[,1],
           col = addTrans(ri.color, 250),
           pch = 16, cex = 1.25)
    
    
    segments(y0 = 1:8, x0 = ma.hfrac.sum[,2], 
             y1 = 1:8, x1 = ma.hfrac.sum[,3],
             col = addTrans(ma.color, 100), 
             lty = 1, lwd = 3)
    
    points(y = 1:8 , x = ma.hfrac.sum[,1], 
           col = addTrans(ma.color, 250),
           pch = 16, cex = 1.25)
    
    
    segments(y0 = 1:8 - 0.1, x0 = pa.hfrac.sum[,2], 
             y1 = 1:8 - 0.1, x1 = pa.hfrac.sum[,3],
             col = addTrans(pa.color, 100), 
             lty = 1, lwd = 3)
    
    points(y = 1:8 - 0.1, x = pa.hfrac.sum[,1], 
           col = addTrans(pa.color, 250),
           pch = 16, cex = 1.25)
    
    legend("bottomright", inset = 0.15,
           legend=c("Rhode Island", "Massachusetts", "Pennsylvania"),
           col = c(addTrans(ri.color, 250),
                   #addTrans(ma.color, 250),
                   addTrans(pa.color, 250)),
           pch = rep(16,3), cex=1.25,
           box.lty = 0)  
  } else {
    
    plot(y = 1:8 + 0.05, x = ri.hfrac.sum[,1], 
         ylim = c(0.5, 8.5), xlim = c(0,0.35),
         xlab = "",
         ylab = "", main = "", col = addTrans(ri.color, 250),
         pch = 16, bty = "l", yaxt = "n", cex = 1.25,
         cex.axis = axis.size, cex.lab = lab.size)
    
    axis(2, at=1:8, labels=c("10-19", "20-29", "30-39", "40-49", "50-59",
                             "60-69", "70-79", "80+"),
         cex.axis = axis.size, cex.lab = lab.size,
         las = 1)
    
    abline(h = 1, col = addTrans("grey", 50), lwd = 2)
    abline(h = 2, col = addTrans("grey", 50), lwd = 2)
    abline(h = 3, col = addTrans("grey", 50), lwd = 2)
    abline(h = 4, col = addTrans("grey", 50), lwd = 2)
    abline(h = 5, col = addTrans("grey", 50), lwd = 2)
    abline(h = 6, col = addTrans("grey", 50), lwd = 2)
    abline(h = 7, col = addTrans("grey", 50), lwd = 2)
    abline(h = 8, col = addTrans("grey", 50), lwd = 2)
    
    segments(y0 = 1:8 + 0.05, x0 = ri.hfrac.sum[,2], 
             y1 = 1:8 + 0.05, x1 = ri.hfrac.sum[,3],
             col = addTrans(ri.color, 100), 
             lty = 1, lwd = 3)
    
    points(y = 1:8 + 0.05, x = ri.hfrac.sum[,1],
           col = addTrans(ri.color, 250),
           pch = 16, cex = 1.25)
    
    
    # segments(y0 = 1:8, x0 = ma.hfrac.sum[,2], 
    #          y1 = 1:8, x1 = ma.hfrac.sum[,3],
    #          col = addTrans(ma.color, 100), 
    #          lty = 1, lwd = 3)
    # 
    # points(y = 1:8 , x = ma.hfrac.sum[,1], 
    #        col = addTrans(ma.color, 250),
    #        pch = 16, cex = 1.25)
    
    
    segments(y0 = 1:8 - 0.05, x0 = pa.hfrac.sum[,2], 
             y1 = 1:8 - 0.05, x1 = pa.hfrac.sum[,3],
             col = addTrans(pa.color, 100), 
             lty = 1, lwd = 3)
    
    points(y = 1:8 - 0.05, x = pa.hfrac.sum[,1], 
           col = addTrans(pa.color, 250),
           pch = 16, cex = 1.25)
    
    legend("bottomright", inset = 0.15,
           legend=c("Rhode Island", #"Massachusetts", "Pennsylvania"),
                    "Pennsylvania"),
           col = c(addTrans(ri.color, 250),
                   #addTrans(ma.color, 250),
                   addTrans(pa.color, 250)),
           pch = rep(16,2), cex=1.25,
           box.lty = 0)   
    
  }
  
  close.screen(all = TRUE)
  
  dev.off()
  embed_fonts(plot.name, outfile = plot.name)

}

### Functions used to determine (100 * alpha)% HPD interval, using the 
###   kernel density estimate for a given set of posterior samples.
###   Portions of this code (highest_alpha and as.pdf) represent slightly
###   modified functions presented by 'whuber' on crossvalidated:
###   (https://stats.stackexchange.com/questions/381520/how-can-i-estimate-
###     the-highest-posterior-density-interval-from-a-set-of-x-y-valu).

library(rootSolve)

highest_alpha <- function(alpha, df, x.min, x.max, y.min, y.max, ...) {
  p <- function(h) {
    g <- function(x) {y <- df(x); ifelse(y > h, y, 0)}
    integrate(g, x.min, x.max, ...)$value - alpha
  }
  uniroot(p, c(y.min, y.max), tol=1e-12)$root
}

as.pdf <- function(x, y, ...) {
  f <- approxfun(x, y, method="linear", yleft=0, yright=0, rule=2)
  const <- integrate(f, min(x), max(x), ...)$value
  approxfun(x, y/const, method="linear", yleft=0, yright=0, rule=2)
}

x.vals <- function(h, df, x.min, x.max){
  f <- function(x) {df(x) - h}
  uniroot.all(f, c(x.min, x.max), tol=1e-12)
}

HPDplot <- function(samp, alpha = 0.95, easy = T, ...){
  
  dens <- density(samp)
  
  if (easy){
    
    hpd.bounds <- hdi(samp, credMass = alpha)
    hpd.inds <- which((dens$x >= hpd.bounds[1]) & 
                        (dens$x <= hpd.bounds[2]))
    
    lower.ind <- which(dens$x == min(dens$x[hpd.inds]))
    upper.ind <- which(dens$x == max(dens$x[hpd.inds]))
    
    hpd.poly <- list(x = c(min(dens$x[hpd.inds]), 
                           dens$x[hpd.inds], 
                           max(dens$x[hpd.inds])),
                     y = c(0, dens$y[hpd.inds], 0))
    
    plot(dens, main = "", bty = "l", xlab = "", ylab = "", ...)
    # polygon(icufrac.dens, col=addTrans("grey", 100), border="black") 
    polygon(hpd.poly, col = addTrans("deepskyblue3", 50), 
            border = addTrans("deepskyblue3", 50))
    
    median.ind <- which.min(abs(dens$x - median(samp)))
    
    segments(x0 = dens$x[median.ind], y0 = 0, 
             x1 = dens$x[median.ind], 
             y1 = dens$y[median.ind], 
             col = "deepskyblue3", lwd = 2)
    
    segments(x0 = dens$x[lower.ind], y0 = 0, 
             x1 = dens$x[lower.ind], 
             y1 = dens$y[lower.ind], 
             col = "deepskyblue3", lty = 3, lwd = 2)
    
    segments(x0 = dens$x[upper.ind], y0 = 0, 
             x1 = dens$x[upper.ind], 
             y1 = dens$y[upper.ind], 
             col = "deepskyblue3", lty = 3, lwd = 2)
    
    lines(dens, lwd = 2)
  } else {
    
    h <- highest_alpha(alpha, as.pdf(dens$x, dens$y),
                       min(dens$x), max(dens$x),
                       min(dens$y), max(dens$y))
    
    x.bounds <- x.vals(h, as.pdf(dens$x, dens$y), 
                       min(dens$x), max(dens$x))
    
    if (length(x.bounds) == 2){
      
      hpd.inds <- dens$y >= h
      hpd.poly <- list(x = c(rep(x.bounds[1], 2), 
                             dens$x[hpd.inds], 
                             rep(x.bounds[2], 2)),
                       y = c(0, h, dens$y[hpd.inds], h, 0))
    } else {
      
      #stop("interval is multimodal; plotting code not yet implemented")
    }
    
    plot(dens, main = "", bty = "l", xlab = "", ylab = "", ...)
    # polygon(icufrac.dens, col=addTrans("grey", 100), border="black") 
    # polygon(hpd.poly, col = addTrans("deepskyblue3", 50), 
    #         border = addTrans("deepskyblue3", 50))
    
    abline(h = h, col = "red", lwd = 2)
    
    median.ind <- which.min(abs(dens$x - median(samp)))
    
    segments(x0 = dens$x[median.ind], y0 = 0, 
             x1 = dens$x[median.ind], 
             y1 = dens$y[median.ind], 
             col = "deepskyblue3", lwd = 2)
    
    segments(x0 = x.bounds[1], y0 = 0, 
             x1 = x.bounds[1], y1 = h, 
             col = "deepskyblue3", lty = 3, lwd = 2)
    
    segments(x0 = x.bounds[2], y0 = 0, 
             x1 = x.bounds[2], 
             y1 = h, 
             col = "deepskyblue3", lty = 3, lwd = 2)
    
    if (length(x.bounds) > 2){
      
      segments(x0 = x.bounds[3], y0 = 0, 
               x1 = x.bounds[3], y1 = h, 
               col = "deepskyblue3", lty = 3, lwd = 2)
      
      segments(x0 = x.bounds[4], y0 = 0, 
               x1 = x.bounds[4], 
               y1 = h, 
               col = "deepskyblue3", lty = 3, lwd = 2)
    }
    
    lines(dens, lwd = 2)
  }
}

hpdBounds <- function(samp, alpha = 0.95){
  
  # density plot estimate
  dens <- density(samp)
  
  
  h <- highest_alpha(alpha, as.pdf(dens$x, dens$y),
                     min(dens$x), max(dens$x),
                     min(dens$y), max(dens$y))
    
  x.bounds <- x.vals(h, as.pdf(dens$x, dens$y), 
                     min(dens$x), max(dens$x))
    
  return(x.bounds)
}

icuTransformation <- function(samp, loc){
 
  if (loc == "RI"){
    pop.frac <- c(0.105, 0.123, 0.140, 0.127, 0.124, 0.135, 0.120, 0.074, 0.052)
  } else if (loc == "MA"){
    pop.frac <- c(0.10586466, 0.12243686, 0.14498786, 0.13384234, 
                  0.12230812, 0.14064248, 0.11801015, 0.06958116, 0.04232637)
  } else if (loc == "PA"){
    pop.frac <- c(0.11160395, 0.12229803, 0.13156525, 0.12581869,
                  0.11809624, 0.13878546, 0.1270166, 0.07657303, 0.04824275)
  }
  
  lewis.prop <- c(0.304, 0.293, 0.2825, 0.301, 0.463, 0.4245,
                  0.46, 0.4835, 0.416)
  
  icu.prob <- sapply(samp, function(x) {sum(x * pop.frac * lewis.prop)}) 
  
  return(icu.prob)
}

paramSummary <- function(beta.directory, ode.directory, data.directory, 
  odepath, loc, const, subsample = NA){
  
  ### make samples data structure for traj.sim
  
  # grab samples from csv files
  betas <- read.csv(beta.directory)
  ode <- read.csv(ode.directory)
  
  # data
  df <- read.csv(data.directory)
  df <- data.frame(df)
  
  if (loc == "PA"){
    
    ode.inds <- 1:34
    
    ### clean up data (CHECK THIS IN FUTURE!)
    
    # cumulative value is wrong for age 80+
    #df[135:137,24] # ..., 11526, 13600, 11645, ...
    # change 13600 to 11600:
    df[136,24] <- 11600
    
    na.70 <- is.na(df[,41])
    df[na.70, 34:42] <- NA
    
    ## removing days before 61
    idx.remove <- which(df$daynum<61)
    if(length(idx.remove) > 0){
      df <- df[-idx.remove,]
    }
    
    ## removing days after 250
    idx.remove <- which(df$daynum>250)
    if(length(idx.remove) > 0){
      df <- df[-idx.remove,]
    }
    
    ## remove NAs
    df <- df[!is.na(df$daynum),]
    
    n.days <- nrow(df)
    days <- df$daynum
    
    end.day <- max(days, na.rm = TRUE) + 1
    num.days <- end.day - 60
    
    # cubic spline with one basis function every 7 days
    bspl <- create.bspline.basis(c(61, end.day), nbasis = round(num.days/7))
    Z <- eval.basis(bspl, 61:end.day)
    
    # iSpline for reporting rate
    knots <- c(84, 92, 100, 108, 130, 160, 190)
    Z.rr <- iSpline(x = 61:end.day, knots = knots, degree = 2)
    Z.rr <- cbind(rep(1, nrow(Z.rr)), Z.rr)
    Z.rr <- Z.rr[-nrow(Z.rr),c(1:3)]
    
    # introday
    id = 60
    
  } else if (loc == "RI"){
    
    ode.inds <- 1:36
    
    ## removing days before 61
    idx.remove=which(df$daynum<61)
    if(length(idx.remove)>0){
      df = df[-idx.remove,]
    }
    
    n.days <- nrow(df)
    days <- df$daynum
    
    ### mcmc initialization (don't change)
    end.day <- max(days, na.rm = TRUE) + 1
    num.days <- end.day - 60 ## number of days with data
    
    # cubic spline with one basis function every 7 days
    bspl <- create.bspline.basis(c(61, end.day), nbasis=round(num.days/7))
    Z <- eval.basis(bspl,61:end.day)
    
    # iSpline for reporting rate
    knots <- c(84, 92, 100, 108, 130, 160, 190)
    Z.rr <- iSpline(x = 61:end.day, knots = knots, degree = 2)
    Z.rr <- cbind(rep(1, nrow(Z.rr)), Z.rr)
    Z.rr <- Z.rr[-nrow(Z.rr),c(1,4:7)]
    
    # introday
    id <- 55
    
  } else if (loc == "MA") {
    
    ode.inds <- 1:35
    
    ## removing days before 61
    idx.remove <- which(df$daynum<61)
    if(length(idx.remove)>0){
      df <- df[-idx.remove,]
    }
    
    n.days <- nrow(df)
    days <- df$daynum
    
    ### mcmc initialization (don't change)
    end.day <- max(days, na.rm = TRUE) + 1
    num.days <- end.day - 60 ## number of days with data
    
    ### Create spline expansions
    
    # Cubic spline with one basis function every 7 days
    bspl <- create.bspline.basis(c(61, end.day), nbasis=round(num.days/7))
    Z <- eval.basis(bspl,61:end.day)
    
    # zero-order spline with one rate
    bspl.rr <- create.bspline.basis(c(61, end.day),nbasis=1,norder=1)
    Z.rr <- eval.basis(bspl.rr, 61:end.day)
    Z.rr <- Z.rr[-nrow(Z.rr),]
    
    # introday
    id <- 55
  }
  
  # list with required structures
  samples <- list(betas = as.matrix(betas),
                  ode = as.matrix(ode[,ode.inds]),
                  rr = as.matrix(ode[,-ode.inds]),
                  days = days,
                  Z.beta = Z,
                  Z.rr = Z.rr,
                  df = df,
                  loc = loc,
                  const = const,
                  introday = id,
                  sf.choice = FALSE)
  
  # reformat column names
  colnames(samples$ode) <- gsub('\\.', '-', colnames(samples$ode))
  colnames(samples$ode) <- gsub("hosp-report-rate", 
                                "hosp.report.rate", 
                                colnames(samples$ode))
  
  # process trajectories from the samples
  tp <- traj.sim(samples = samples,
                 odepath = odepath,
                 csv = T)
  
  
  ### determine ODE param values
  
  n.ode <- ncol(samples$ode)
  ode.summary <- data.frame(names = colnames(samples$ode),
                            median = rep(NA_real_, n.ode),
                            lower.ci = rep(NA_real_, n.ode),
                            upper.ci = rep(NA_real_, n.ode),
                            lower.hpd = rep(NA_real_, n.ode),
                            upper.hpd = rep(NA_real_, n.ode),
                            lower.hpd2 = rep(NA_real_, n.ode),
                            upper.hpd2 = rep(NA_real_, n.ode),
                            lower.hpd3 = rep(NA_real_, n.ode),
                            upper.hpd3 = rep(NA_real_, n.ode))
             
  for (k in 1:n.ode){
    
    if ((colnames(samples$ode)[k] == "dev-icu-frac") | 
        (colnames(samples$ode)[k] == "dev-icu-frac-phase2")){
      
      ode.k <- icuTransformation(samp = samples$ode[,k], loc = samples$loc)
      
    } else if (colnames(samples$ode)[k] == "dev-len-hospstay"){
      ode.k <- 10 * samples$ode[,k]
    } else {
      ode.k <- samples$ode[,k]
    }
    
    # median and credible interval
    ci <- quantile(ode.k, probs = c(0.025, 0.5, 0.975))
    ode.summary[k, 2:4] <- ci[c(2,1,3)]
    
    # hpd 
    hpd.k <- hpdBounds(ode.k, alpha = 0.95)
    ode.summary[k, 4 + 1:length(hpd.k)] <- hpd.k
  }
  
  ### determine reporting rate bounds
  
  rr.summary <- data.frame(date = c("March 1", "March 15", "June 6"),
                           median = rep(NA_real_, 3),
                           lower.ci = rep(NA_real_, 3),
                           uppder.ci = rep(NA_real_, 3))
  
  rr.summary[1, 2:4] <- quantile(tp$rr.full[,1], prob = c(0.5, 0.025, 0.975))
  rr.summary[2, 2:4] <- quantile(tp$rr.full[,15], prob = c(0.5, 0.025, 0.975))
  rr.summary[3, 2:4] <- quantile(tp$rr.full[,190], prob = c(0.5, 0.025, 0.975))
  
  
  return(list(ode = ode.summary,
              rr = rr.summary))
}



###
### 1. Read in functions.
###

OPT_USE_PY_LL=FALSE


source("./plot.chains.R")
source("./data.process.R")
source("./traj.process.R")
source("./traj.from.params.R")
source("./loglik.odesim.4.0.R", chdir = TRUE)
source("./results.plots.and.params.R")
source("./fancy-plots.R")

# register fonts with R's pdf output device
library('extrafont')
loadfonts()

library(HDInterval)
library(fda)
library(splines2)
library(cowplot)

### constant vectors:

const.vec.ri <- c("", 2)
names(const.vec.ri) <- c("symp-frac-davies", "steps-per-day")

const.vec.ma <- c(2, "")
names(const.vec.ma) = c("steps-per-day", "symp-frac-davies")

const.vec.pa <- c("", 0.8, 2)
names(const.vec.pa) <- c("symp-frac-davies", "prob-icu-vent", "steps-per-day")


### ri plots ###

set.seed(94)

fig2panel(beta.directory = "../../final_outputs_for_ms/pre-print/RI/RI_bestfit-results.daily.betas-day-250.csv", 
          ode.directory = "../../final_outputs_for_ms/pre-print/RI/RI_bestfit-results.ode.params-day-250.csv",
          data.directory = "../../data/Rhode_Island/RI_formatted_20200906.csv", 
          odepath = "../../cpp-v5-discharges-nonhospdeaths/", 
          loc = "RI", const = const.vec.ri, 
          plot.name = "../../final_outputs_for_ms/pre-print/RI/RI-ms-figure2.pdf", alpha = 50, 
          subsample = 250, 
          axis.size = 1.75, lab.size = 1.75,
          ncol = 2, nrow = 4, height = 12, width = 12,
          font.family = "Ubuntu Mono")

set.seed(1967)

fig4(beta.directory = "../../final_outputs_for_ms/pre-print/RI/RI_bestfit-results.daily.betas-day-250.csv", 
     ode.directory = "../../final_outputs_for_ms/pre-print/RI/RI_bestfit-results.ode.params-day-250.csv",
     data.directory = "../../data/Rhode_Island/RI_formatted_20200906.csv", 
     odepath = "../../cpp-v5-discharges-nonhospdeaths/",
     loc = "RI", const = const.vec.ri, 
     axis.size = 12, title.size = 14, 
     legend.size = 10, legend.title.size = 12,
     plot.name = "../../final_outputs_for_ms/pre-print/RI/RI-ms-fig4b.pdf",  
     height = 4, width = 6, font.family = "Ubuntu Mono")

set.seed(1255)

ri.summary <- paramSummary(beta.directory = "../../final_outputs_for_ms/pre-print/RI/RI_bestfit-results.daily.betas-day-250.csv", 
                           ode.directory = "../../final_outputs_for_ms/pre-print/RI/RI_bestfit-results.ode.params-day-250.csv",
                           data.directory = "../../data/Rhode_Island/RI_formatted_20200906.csv", 
                           odepath = "../../cpp-v5-discharges-nonhospdeaths/",
                           loc = "RI",
                           const = const.vec.ri,
                           subsample = NA)

### pa plots ###

set.seed(302)

fig2panel(beta.directory = "../../final_outputs_for_ms/pre-print/PA/PA_bestfit-results.daily.betas-day-250.csv", 
          ode.directory = "../../final_outputs_for_ms/pre-print/PA/PA_bestfit-results.ode.params-day-250.csv",
          data.directory = "../../data/Pennsylvania/PA_most_recent_data.csv", 
          odepath = "../../cpp-v5-discharges-nonhospdeaths/", 
          loc = "PA", const = const.vec.pa, 
          plot.name = "../../final_outputs_for_ms/pre-print/PA/PA-ms-figure2.pdf", 
          alpha = 50, 
          subsample = 250, 
          axis.size = 1.75, lab.size = 1.75,
          ncol = 2, nrow = 4, height = 12, width = 12,
          font.family = "Ubuntu Mono")

set.seed(92)

fig4(beta.directory = "../../final_outputs_for_ms/pre-print/PA/PA_bestfit-results.daily.betas-day-250.csv", 
     ode.directory = "../../final_outputs_for_ms/pre-print/PA/PA_bestfit-results.ode.params-day-250.csv",
     data.directory = "../../data/Pennsylvania/PA_most_recent_data.csv", 
     odepath = "../../cpp-v5-discharges-nonhospdeaths/", 
     loc = "PA", const = const.vec.pa, 
     axis.size = 12, title.size = 14, 
     legend.size = 10, legend.title.size = 12,
     plot.name = "../../final_outputs_for_ms/pre-print/PA/PA-ms-fig4b.pdf",  
     height = 4, width = 6, font.family = "Ubuntu Mono")

set.seed(93)

pa.summary <- paramSummary(beta.directory = "../../final_outputs_for_ms/pre-print/PA/PA_bestfit-results.daily.betas-day-250.csv", 
                           ode.directory = "../../final_outputs_for_ms/pre-print/PA/PA_bestfit-results.ode.params-day-250.csv",
                           data.directory = "../../data/Pennsylvania/PA_most_recent_data.csv", 
                           odepath = "../../cpp-v5-discharges-nonhospdeaths/", 
                           loc = "PA",
                           const = const.vec.pa,
                           subsample = NA)

### ma plots ###

set.seed(31)

fig2panel(beta.directory = "../../final_outputs_for_ms/pre-print/MA/MA-1104-cpi.daily.betas-day-250.csv",
          ode.directory = "../../final_outputs_for_ms/pre-print/MA/MA-1104-cpi.ode.params-day-250.csv",
          data.directory = "../../data/Massachusetts/Past_data/MA_20200906_data_7daysmoothed.csv",
          odepath = "../../cpp-v5-discharges-nonhospdeaths/",
          loc = "MA", const = const.vec.ma, 
          plot.name = "../../final_outputs_for_ms/pre-print/MA/MA-ms-figure2.pdf",
          alpha = 50, subsample = 250,
          axis.size = 1.75, lab.size = 1.75, ncol = 2, nrow = 4, height = 12, width = 12,
          font.family = "Ubuntu Mono")

set.seed(78)

fig4(beta.directory = "../../final_outputs_for_ms/pre-print/MA/MA-1104-cpi.daily.betas-day-250.csv",
     ode.directory = "../../final_outputs_for_ms/pre-print/MA/MA-1104-cpi.ode.params-day-250.csv",
     data.directory = "../../data/Massachusetts/Past_data/MA_20200906_data_7daysmoothed.csv",
     odepath = "../../cpp-v5-discharges-nonhospdeaths/",
     loc = "MA", const = const.vec.ma, 
     axis.size = 12, title.size = 14, 
     legend.size = 10, legend.title.size = 12,
     plot.name = "../../final_outputs_for_ms/pre-print/MA/MA-ms-fig4b.pdf",  
     height = 4, width = 6, font.family = "Ubuntu Mono")



set.seed(18)

ma.summary <- paramSummary(beta.directory = "../../final_outputs_for_ms/pre-print/MA/MA-1104-cpi.daily.betas-day-250.csv",
                           ode.directory = "../../final_outputs_for_ms/pre-print/MA/MA-1104-cpi.ode.params-day-250.csv",
                           data.directory = "../../data/Massachusetts/Past_data/MA_20200906_data_7daysmoothed.csv",
                           odepath = "../../cpp-v5-discharges-nonhospdeaths/",
                           loc = "MA",
                           const = const.vec.ma,
                           subsample = NA)


### FIGURE FIVE ###
set.seed(98366)

fig5panel(ri.beta = "../../final_outputs_for_ms/pre-print/RI/RI_bestfit-results.daily.betas-day-250.csv",
          ri.ode = "../../final_outputs_for_ms/pre-print/RI/RI_bestfit-results.ode.params-day-250.csv",
          ri.data = "../../data/Rhode_Island/RI_formatted_20200906.csv",
          ri.const = const.vec.ri,
          ma.beta = "../../final_outputs_for_ms/pre-print/MA/MA-1104-cpi.daily.betas-day-250.csv",
          ma.ode = "../../final_outputs_for_ms/pre-print/MA/MA-1104-cpi.ode.params-day-250.csv",
          ma.data = "../../data/Massachusetts/Past_data/MA_20200906_data_7daysmoothed.csv",
          ma.const = const.vec.ma,
          pa.beta = "../../final_outputs_for_ms/pre-print/PA/PA_bestfit-results.daily.betas-day-250.csv",
          pa.ode = "../../final_outputs_for_ms/pre-print/PA/PA_bestfit-results.ode.params-day-250.csv",
          pa.data = "../../data/Pennsylvania/PA_most_recent_data.csv",
          pa.const = const.vec.pa,
          plot.name = "../../final_outputs_for_ms/pre-print/fig5.pdf",
          alpha = 75,
          subsample = 250,
          odepath = "../../cpp-v5-discharges-nonhospdeaths/",
          axis.size = 1,
          lab.size = 1,
          font.family = "Ubuntu Mono")


ri.ode = "../../final_outputs_for_ms/pre-print/RI/RI_bestfit-results.ode.params-day-250.csv"
ode <- read.csv(ri.ode)

histFig(ode = ode[,1:36],
  nr = 6,
  nc = 6,
  plot.name = "test.pdf",
  font.family = "Ubuntu Mono",
  height = 16,
  width = 16)


histFig <- function(ode, nr = 6, nc = 6, plot.name, 
                    font.family = "Ubuntu Mono", height = 16, width = 16){
  
  pdf(file = plot.name, family = font.family, h = height, w = width)
  
  # grid of plots
  par(mfrow = c(nr,nc))
  par(mar = c(3, 2, 3, 2), oma = c(2, 2, 4, 2))
  
  for (i in 1:ncol(ode)){
    hist(ode[,i], main = "", xlab = "", ylab = "", col = i)
    title(colnames(ode)[i], adj = 0,
          line = 1, cex.main = 1.75)
  }
  
  title("Rhode Island", outer = T, adj = 0, line = 1, cex.main = 4)
  
  dev.off()
  embed_fonts(plot.name, outfile = plot.name)
}



