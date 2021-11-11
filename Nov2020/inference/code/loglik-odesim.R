#!/usr/bin/env Rscript

# loglik-odesim.R
# authors: Ephraim Hanks, Nathan Wikle, Emily Strong
# last edited: 20 Nov 2020 
#
# This file defines the functions used to calculate the ikelihoood for 
#   all odesim parameters. In an effort to speed up the base-r
#   implementation of 'dmultinom', a Python-based multinomial logpmf
#   function is defined -- this requires the embedding of a Python 
#   session within R via the 'reticulate' packaged. If Python is 
#   unavailable, a global variable called 'OPT_USE_PY_LL' must be 
#   set to FALSE. If not, an error will be thrown. If Python is 
#   enabled, this need not be defined.


### 1. Check if Python is enabled, write Python-based multinomial logpmf

# get global options, OPT_USE_PY_LL must == FALSE if Python is not enabled
OPT_USE_PY_LL = get0('OPT_USE_PY_LL', ifnotfound = TRUE)

# helper 1 - get sums of times spans bounded by times
time_slicer = function(vals, times){
  out = matrix(nrow=length(times)-1, ncol=ncol(vals))
  for(kk in 2:length(times)){
    this_slice=vals[(times[kk-1]+1):times[[kk]],]
    if(length(dim(this_slice))==2){
      this_slice = apply(this_slice,2,sum)
    }
    out[kk-1,] = this_slice
  }
  return(out)
}

# helper 2 - sum up log likelihood for several rows
r_ll = function(means, counts){
 row_lls = sapply(1:nrow(means),
                  function (i) dmultinom(round(counts[i,]), 
                                  prob=means[i,], log=TRUE))
 return(sum(row_lls))
}

# helper 3 - Python-based multinomial logpmf
if (OPT_USE_PY_LL){
  tryCatch({
    library(reticulate)
    source_python('tools/py_faster_stats.py') # path is relative to this file, remember to source with chdir=TRUE
   }, error=function(e){
     stop(c("Unable to import Python multinomial function.  Set OPT_USE_PY_LL=FALSE to use R version.", "\n\n", e))
     }
   )
  py_ll = py$py_ll
  multinom_log_pmf = function(x,y) py_ll(as.matrix(x), as.matrix(y))
} else {
  multinom_log_pmf = r_ll
}

### 2. Define loglik.odesim function to calculate loglikelihood. 
###       Outputs list of log-likelihood values.

loglik.odesim <- function(
  traj,                         # trajectories generated via odesim (using traj.from.params)
  df,                           # state-level covid data
  dp = NULL,                    # formatted state-level covid data (using data.process)
  odesim.ver = "v5",            # version of odesim (defaults to "v5")
  loc = loc,                    # US state used for analysis (one of "RI", "MA", or "PA")
  report.rate,                  # vector of daily symptomatic reporting rates
  nb.disp.params = NULL,        # vector neg. binomial dispersion parameters
  extra.params = NULL,          # fitted params not found in odesim (most likely hosp.report.rate)
  extra.const.params = NULL,    # constant params not found in odesim (most likely hosp.report.rate)
  s2.hosp = NULL,               # variance param of current hosp. likelihood
  s2.icu = NULL,                # variance param of current icu likelihood
  s2.vent = NULL,               # variance param of current vents likelihood
  lik.tot = TRUE,               # evaluate likelihood for total new cases ONLY (ie, set F and use lik.age for age-strat. data; default = TRUE)
  lik.age = FALSE,              # evaluate likelihood for age-struc. new cases and hosps AND total new cases and hosps (default = FALSE)
  lik.hosp.new = TRUE,          # evaluate likelihood for new hosp. cases (default = FALSE)
  lik.hosp.curr = FALSE,        # evaluate likelihood for current hosp. (default = FALSE)
  lik.icu.curr = FALSE,         # evaluate likelihood for current icu admits (default = FALSE)
  lik.vent.curr = FALSE,        # evaluate likelihood for current vent admits (default = FALSE)
  lik.tot.deaths = FALSE,       # evaluate likelihood for tot. deaths ONLY (ie, set F and use lik.age.deaths for age-strat. data; default = FALSE)
  lik.age.deaths = FALSE,       # evaluate likelihood for age-struc. new deaths and total new deaths (default = FALSE)
  lik.home.deaths = FALSE,      # evaluate likelihood for new home deaths (default = FALSE)
  lik.hosp.deaths = FALSE,      # evaluate likelihood for new hosp. deaths (default = FALSE)
  lik.hosp.discharges = FALSE,  # evaluate likelihood for hospital discharges (default = FALSE)
  active.surv = FALSE,          # include active surveillance data (default = FALSE)
  p.asympt = 0.4,               # proportion of asymptomatic individuals (default = 0.4; CAUTION: DO NOT CHANGE UNLESS ODESIM REFLECTS A 
  case.constraint = FALSE,      # constrain fit to cumulative cases to be within 10% of data (default = FALSE)
  P = P,                        # testing delay matrix, computed from 'p.vecs' input in mcmc-odesim.R
  ...
){     

  ##################################
  ## handling "extra parameters" ###
  ##################################

  ## default values go here
  hosp.rr=1
  
  ## reading in values from "extra.params"
  if(length(extra.params)>0){
    extra.names=names(extra.params)
    ## hosp.report.rate
    if(is.element("hosp.report.rate",extra.names)){
      hosp.rr=extra.params["hosp.report.rate"]
    }
    if(is.element("p.asympt",extra.names)){
      p.asympt=extra.params["p.asympt"]
    }
    ## add others here
  }

  ## reading in values from "extra.params"
  if(length(extra.const.params)>0){
    extra.const.names=names(extra.const.params)
    ## hosp.report.rate
    if(is.element("hosp.report.rate",extra.const.names)){
      hosp.rr=as.numeric(extra.const.params["hosp.report.rate"])
    }
    ## add others here
  }
    
  #########################################
  ### parsing out NB dispersions params ###
  #########################################
    
  nb.disp.new.cases=NULL ## leave "NULL" to do Poisson lik
  nb.disp.new.hosp=NULL ## leave "NULL" to do Poisson lik
  nb.disp.deaths=NULL ## leave "NULL" to do Poisson lik
  nb.disp.home.deaths=NULL ## leave "NULL" to do Poisson lik
  nb.disp.hosp.deaths=NULL ## leave "NULL" to do Poisson lik
  if(length(nb.disp.params)>=1){
    nb.disp.new.cases=nb.disp.params[1]
  }
  if(length(nb.disp.params)>=2){
    nb.disp.new.hosp=nb.disp.params[2]
  }
  if(length(nb.disp.params)>=3 & (lik.tot.deaths | lik.age.deaths)){
    nb.disp.deaths=nb.disp.params[3]
  }
  if(length(nb.disp.params)==4 & lik.home.deaths & lik.hosp.deaths){
    nb.disp.home.deaths=nb.disp.params[3]
    nb.disp.hosp.deaths=nb.disp.params[4]
  }
  if(length(nb.disp.params)==4 & lik.hosp.discharges){
    nb.disp.hosp.discharges=nb.disp.params[4]
  }
  if(length(nb.disp.params)==5 & lik.home.deaths & lik.hosp.deaths & lik.hosp.discharges){
    nb.disp.home.deaths=nb.disp.params[3]
    nb.disp.hosp.deaths=nb.disp.params[4]
    nb.disp.hosp.discharges=nb.disp.params[5]
  }
    
  ################################################################################
  ###                                                                          ###
  ### NOTE that this function REQUIRES a function called "data.process", which ### 
  ###   may take other parameters (via ...)                                    ###
  ###                                                                          ###  
  ################################################################################
    
  # symptomatic reporting rate for 
  report.rate <- as.numeric(report.rate)
   
  ##############################
  ### pre-processing of data ###
  ##############################
  if(is.null(dp)){  #Allow precomputed data.process
    data.processed=data.process(df,loc=loc,...)
  } else {
    data.processed=dp
  }
  list2env(data.processed,globalenv())
  n.days=length(days)
        
  #######################################
  ### pre-processing of odesim output ###
  #######################################
  
  traj.processed=traj.process(traj,loc=loc,odesim.version=odesim.ver)
  list2env(traj.processed,globalenv())

  ### linking times of data and odesim output
  traj.times=round(traj[,1])
  tmin <- min(days); tmax <- max(days)
  idx.min <- which(traj.times == tmin)
  idx.max <- which(traj.times == tmax)
  t.idx <- idx.min:idx.max
        
  ### linking to report.rate (which is daily starting at day 61)
  t.idx.rr=t.idx-10
    
  ### linking to P and delay rate
  t.idx.delay=tmin:tmax
        
  ############################################################################
  ###                                                                      ###
  ### Check if parameters are in bounds, and that case.constraint is ###
  ###   within 10% of observed total (return -Inf if violated)             ###
  ###                                                                      ###  
  ############################################################################

  tot.size.ep <- sum(tot.sympt.new, na.rm = T)  #cumsum(tot.sympt.new, na.rm = T)
    tot.size.ep.odesim <- cumsum(tot.sympt.new.odesim[t.idx] * report.rate[t.idx.rr]) # CHANGED TO MULTIPLICATION!
    tot.size.ep.odesim=max(tot.size.ep.odesim,na.rm=TRUE) ## get last, biggest, value
  tot.size.percent.off <- abs(tot.size.ep.odesim - tot.size.ep) / tot.size.ep
  # tot.size.percent.off <- tot.size.percent.off[length(tot.size.percent.off)]

   ## browser()
    
  ### check if conditions are met
  if(min(report.rate,nb.disp.new.cases,nb.disp.new.hosp,s2.hosp,s2.icu,s2.vent) < 0 | 
    (case.constraint & tot.size.percent.off > 0.10)){ # changed to being greater than 0.1
    
    # parameter out of bounds or total size constraint not met, return -Inf
    ll=-Inf
    ll.hosp.new=NA
    ll.hosp=NA
    ll.vent=NA
    ll.icu=NA
    ll.new=-Inf
    ll.deaths=NA
    ll.hosp.dis=NA
    sympt.imputed=NULL
    
  } else {
    
    ### active surveillance calculations
    
        
    # a_p,t  =  number of positive tests in active surveilance
    # a_n,t  =  number of negative tests in active surveilance
    # c_t  = coverage of the active surv program = (a_p,t  + a_n,t) / N
    # N = 6,897,200
    # p_NCAS     = probability you are not caught by active surveillance
    # = (1 - c_t)^6
    #    
    # pi_t   =  measured prevalence   = (E + I + A)  / ( S + E + (1-RHO)*I + A + (P_ASYMP + P_SYMP*(1-RHO)) *R )
    # 
    # S is columns 2-10
    # E is columns 11-64, inclusive
    # A is columns 65-100
    # I is columns 101-136
    # R is columns 272-280
    # 
    # P_SYMP = 0.6
    # P_ASYMP = 0.4
    # 
    # Like(Data|Model) = Prob[ Poisson( rho * DeltaJ_t * p_NCAS ) == NUM_CASES_TODAY - NUM_CASES_ACTIVE_SURV ) *  Prob[ Binomial( a_n,t+a_p,t , pi_t  ) == a_p,t ]
     
    # if(active.surv){
    #   N.total.pop.nursinghome=34360 
    #   c.t=.01 ## value for RI - see Supplement p7
    #   if(loc=="MA"){
    #     c.t=1497/34360 # = mean tests per day (from MA data 8/3/2020) / total nursinghome pop in MA
    #     # source for nuursinghome pop: https://www.kff.org/other/state-indicator/number-of-nursing-facility-residents/?currentTimeframe=0&sortModel=%7B%22colId%22:%22Location%22,%22sort%22:%22asc%22%7D
    #   }
    #   as.wts=c(.35,.65) ## relative probability of being tested if (70-80,80+)
    #   p.ncas=(1-c.t)^6 ## probability of not being tested for 6 days
    #   p.ncas[is.na(p.ncas)]=mean(p.ncas,na.rm=TRUE)
    #   c.t.mean=mean(c.t,na.rm=TRUE)
    # } else {
    #   p.ncas=rep(1, length(t.idx.rr))
    # }
        
    ##########################
    ### initialize outputs ###
    ##########################
    
    ll.new <- integer()
    ll.hosp <- integer()
    ll.deaths <- integer()
    sympt.imputed <- integer()
    ll <- 0
        
    ######################################
    ### likelihood for daily new cases ###
    ######################################
  
    
    ### likelihood for age-structured new cases
    if(lik.age & !active.surv){
          
      # account for testing delay:
      delay.rates <- as.vector(P[,traj.times] %*% tot.sympt.new.odesim) +.Machine$double.eps[1]
      # finding any na.values
      no.na=which(!is.na(tot.sympt.new))
      # log-likelihood of total new cases each day:
      ll.tot = sum(dnbinom(tot.sympt.new[no.na],mu=report.rate[t.idx.rr[no.na]] * delay.rates[t.idx.delay[no.na]], size=nb.disp.new.cases,log = TRUE))
      # log-likelihood for age-structured data
      age.obs.times=c(0,age.obs.times)

      odesim.means.mat = time_slicer(sympt.new.odesim, age.obs.times)
      
      zero.cases <- rowSums(age.new.cases) == 0
      ll.age = multinom_log_pmf(odesim.means.mat[!zero.cases,], 
                                age.new.cases[!zero.cases,])
      
      ll.new <- ll.tot + ll.age
      ll <- ll + ll.new
    }
    
    if(lik.age & active.surv){
      sigma.e=0.8 ## active surveillance test efficacy
      sigma.a=0.9 ## 
      p.ncas=(1-sigma.e*c.t[t.idx.rr-1,])^6
      # making SEIAR for all age classes 
      i9=diag(9)
      S=traj[t.idx.rr-1,2:10]
      E=traj[t.idx.rr-1,c(11:19,20:28,29:37,38:46,47:55,56:64)]%*%rbind(i9,i9,i9,i9,i9,i9)
      A=traj[t.idx.rr-1,c(65:73,74:82,83:91,92:100)]%*%rbind(i9,i9,i9,i9)
      I=traj[t.idx.rr-1,c(101:109,110:118,119:127,128:136)]%*%rbind(i9,i9,i9,i9)
      R=traj[t.idx.rr-1,272:280]
      prevalence=(sigma.e*E+(1-report.rate[t.idx.rr])*I+sigma.a*A)/(S+E+(1-report.rate[t.idx.rr])*I+A+(p.asympt + (1-p.asympt)*(1-report.rate[t.idx.rr]))*R)
      # # this line sums over 70-79 and 80+
      # measured.prevalence.summed=apply((E+I+A),1,sum)/apply(S+E+(1-report.rate[t.idx.rr])*I+A+(p.asympt + (1-p.asympt)*(1-report.rate[t.idx.rr]))*R,1,sum)
      # idx.as=which(!is.na(pos.as))
      # ll.active.surv=sum(dbinom(pos.as[idx.as],tests.as[idx.as],measured.prevalence.summed[idx.as],log=TRUE))
                        
      if (loc == "RI"){
        ### used to calculate expected number of positive tests from active surveillance
        pop.by.age <- c(111233, 130301, 148311, 134539, 131361, 
                         143014, 127123, 78393, 55087)
      }
      
      # calculate mean new cases (passive and active)
      delay.rates <- P[,traj.times] %*% sympt.new.odesim +.Machine$double.eps[1]
      mean.passive.active= report.rate[t.idx.rr] * delay.rates[t.idx.delay,]*p.ncas+(c.t[t.idx.rr-1,]*prevalence) %*% diag(pop.by.age)
      tot.mean.pass.act=apply(mean.passive.active,1,sum)
      # finding any na.values
      no.na=which(!is.na(tot.sympt.new))
      # log-likelihood of total new cases each day:
      ll.tot = sum(dnbinom(tot.sympt.new[no.na],mu=tot.mean.pass.act[no.na], size=nb.disp.new.cases,log = TRUE))
      # log-likelihood for age-structured data
      age.obs.times=c(0,age.obs.times)
      
      
      odesim.means.mat = time_slicer(mean.passive.active, age.obs.times)
      ll.age = multinom_log_pmf(odesim.means.mat, age.new.cases)
          
      ll.new <- ll.tot + ll.age 
      ll <- ll + ll.new
    }
    
        
    ###########################################################
    ### likelihood for total new cases (not age-structured) ###
    ###########################################################
    if(lik.tot){
      # account for testing delay:
      delay.rates <- as.vector(P[,traj.times] %*% tot.sympt.new.odesim) +.Machine$double.eps[1]
      # finding any na.values
      no.na <- which(!is.na(tot.sympt.new))
      # log-likelihood given trajectory sim:
      ll.new <- sum(dnbinom(tot.sympt.new[no.na],mu=report.rate[t.idx.rr[no.na]] * delay.rates[t.idx.delay[no.na]], size=nb.disp.new.cases,log = TRUE))
      ll <- ll + ll.new
    }

    #######################################################################
    ### additional likelihood for age structured new hospitalized cases ###
    #######################################################################

    if(lik.hosp.new & lik.age){
          
      # finding any na.values
      no.na=which(!is.na(tot.hosp.new))
      # log-likelihood of total new hospitalizations each day:
      ll.hosp.tot = sum(dnbinom(tot.hosp.new[no.na],mu=tot.hosp.new.odesim[no.na]*hosp.rr +.Machine$double.eps[1], size=nb.disp.new.hosp,log = TRUE))
      # log-likelihood for age-structured data
      age.hosp.obs.times=c(0,age.hosp.obs.times)
      # set negative values to zero
      hosp.new.odesim[hosp.new.odesim<=0]=.Machine$double.eps[1]
      ll.hosp.age = 0
      
      for(kk in 2:length(age.hosp.obs.times)){
        odesim.means=hosp.new.odesim[(age.hosp.obs.times[kk-1]+1):age.hosp.obs.times[[kk]],]
        if(length(dim(odesim.means))==2){
          odesim.means=apply(odesim.means,2,sum)
        }
        if (loc == "PA"){
          na.inds <- which(is.na(round(as.numeric(age.hosp.new[kk-1,]))))
          age.hosp.new.nona <- age.hosp.new[,-na.inds]
          odesim.means.nona <- odesim.means[-na.inds]
          ll.hosp.age=ll.hosp.age+dmultinom(round(as.numeric(age.hosp.new.nona[kk-1,])),prob=as.numeric(odesim.means.nona),log=TRUE)
        } else {
          ll.hosp.age=ll.hosp.age+dmultinom(round(as.numeric(age.hosp.new[kk-1,])),prob=as.numeric(odesim.means),log=TRUE)
        }
      }
          
      ll.hosp.new <- ll.hosp.tot + ll.hosp.age
      ll <- ll + ll.hosp.new
    } else if (lik.hosp.new & !lik.age) {
      ### AGE-STRATIFIED HOSPITALIZATIONS, NO AGE CASES...
        
      # finding any na.values
      no.na=which(!is.na(tot.hosp.new))
      # log-likelihood of total new hospitalizations each day:
      ll.hosp.tot = sum(dnbinom(tot.hosp.new[no.na],mu=tot.hosp.new.odesim[no.na]*hosp.rr +.Machine$double.eps[1], size=nb.disp.new.hosp,log = TRUE))
      # log-likelihood for age-structured data
      age.hosp.obs.times=c(0,age.hosp.obs.times)
      # set negative values to zero
      hosp.new.odesim[hosp.new.odesim<=0]=.Machine$double.eps[1]
    
      odesim.means.mat = time_slicer(hosp.new.odesim, age.hosp.obs.times)
      ll.hosp.age = multinom_log_pmf(odesim.means.mat, age.hosp.new)
          
      ll.hosp.new <- ll.hosp.tot + ll.hosp.age
      ll <- ll + ll.hosp.new
    } else {
      ll.hosp.new <- 0
      ll <- ll + ll.hosp.new
    }
      
    ############################################################
    ### likelihood for total new deaths (not age-structured) ###
    ############################################################
    if(lik.tot.deaths){
      # finding any na.values
      no.na=which(!is.na(tot.deaths.new))
      # log-likelihood given trajectory sim:
      tot.deaths.new.odesim[tot.deaths.new.odesim<=0]=0 +.Machine$double.eps[1]
      ll.deaths = sum(dnbinom(tot.deaths.new[no.na],mu=tot.deaths.new.odesim[t.idx[no.na]], size=nb.disp.deaths,log = TRUE,na.rm=TRUE))
      ll <- ll + ll.deaths
    }
        
    ################################################    
    ### likelihood for age-structured new deaths ###
    ################################################
    if(lik.age.deaths){
      # finding any na.values
      no.na <- which(!is.na(tot.deaths.new))
      # log-likelihood for total deaths:
      tot.deaths.new.odesim[tot.deaths.new.odesim <= 0] <- 0 +.Machine$double.eps[1]
      ll.tot.deaths <- sum(dnbinom(tot.deaths.new[no.na],mu=tot.deaths.new.odesim[t.idx[no.na]], size=nb.disp.deaths,log = TRUE))## tot.deaths.new.odesim[tot.deaths.new.odesim<=0]=0 +.Machine$double.eps[1]

      # home and hosp deaths
      ll.home.hosp.deaths <- 0
      if(lik.hosp.deaths & lik.home.deaths & !(loc=="PA") ){
        tot.hosp.deaths.new[tot.hosp.deaths.new < 0] = 0
        tot.home.deaths.new[tot.home.deaths.new < 0] = 0
        no.na <- which(!is.na(tot.hosp.deaths.new+tot.home.deaths.new))
        ll.home.hosp.deaths <- sum(dmultinomial(as.matrix(cbind(tot.hosp.deaths.new,
                                                                tot.home.deaths.new)[no.na,]),
                                                prob=as.matrix(cbind(tot.hosp.deaths.new.odesim[t.idx],
                                                                     tot.home.deaths.new.odesim[t.idx]))[no.na,],
                                                log=TRUE))
      }
      
      if(lik.hosp.deaths & lik.home.deaths & (loc=="PA") ){
        ll.home.hosp.deaths <- 0
        home.deaths.cum <- tot.deaths.cum[hosp.deaths.obs.times] - cumsum(hosp.deaths.new)
        home.deaths.new <- c(home.deaths.cum[1], 
                             home.deaths.cum[-1] - home.deaths.cum[-length(home.deaths.cum)])
        
        for(kk in 1:length(hosp.deaths.obs.times)){
          if (kk == 1){
            odesim.means <- c(sum(tot.hosp.deaths.new.odesim[1:hosp.deaths.obs.times[kk]]),
              sum(tot.home.deaths.new.odesim[1:hosp.deaths.obs.times[kk]])) + .Machine$double.eps[1]
          } else {
            odesim.means <- c(sum(tot.hosp.deaths.new.odesim[(hosp.deaths.obs.times[kk-1]+1):hosp.deaths.obs.times[[kk]]]), 
                              sum(tot.home.deaths.new.odesim[(hosp.deaths.obs.times[kk-1]+1):hosp.deaths.obs.times[[kk]]])) + .Machine$double.eps[1]
          }

          if(!is.na(hosp.deaths.new[[kk]]+home.deaths.new[[kk]])){
              ll.home.hosp.deaths <- ll.home.hosp.deaths + 
                  dmultinom(round(c(hosp.deaths.new[kk], home.deaths.new[kk])),
                            prob = as.numeric(odesim.means), log=TRUE)
          }
        }
      }
            
      ll.deaths.age = 0
      # log-likelihood for age-structured data
      age.deaths.obs.times = c(0,age.deaths.obs.times)
      # set negative values to zero
      deaths.new.odesim[deaths.new.odesim<=0]=.Machine$double.eps[1]
      odesim.means.mat = time_slicer(deaths.new.odesim, age.deaths.obs.times)
      
      no.zero.deaths <- which(rowSums(age.deaths.new, na.rm = T) != 0) 
      
      ll.deaths.age = multinom_log_pmf(odesim.means.mat[no.zero.deaths,], 
                                       age.deaths.new[no.zero.deaths,])
            
      ll.deaths <- ll.tot.deaths + ll.deaths.age + ll.home.hosp.deaths
      ll <- ll + ll.deaths
    }

        


# 
#         
#     ##
#     ## likelihood for hosp/home new deaths (not age-structured)
#     ##
#     if(lik.hosp.deaths & lik.home.deaths & !lik.age.deaths){
# 
#         ## finding any na.values
#         no.na=which(!is.na(tot.hosp.deaths.new))
#         ## log-likelihood given trajectory sim:
#         tot.hosp.deaths.new.odesim[tot.hosp.deaths.new.odesim<=0]=0 +.Machine$double.eps[1]
#         ll.hosp.deaths = sum(dnbinom(tot.hosp.deaths.new[no.na],mu=tot.hosp.deaths.new.odesim[t.idx[no.na]], size=nb.disp.hosp.deaths,log = TRUE))
#         ## finding any na.values
#         no.na=which(!is.na(tot.home.deaths.new))
#         ## log-likelihood given trajectory sim:
#         tot.home.deaths.new.odesim[tot.home.deaths.new.odesim<=0]=0 +.Machine$double.eps[1]
#         ll.home.deaths = sum(dnbinom(tot.home.deaths.new[no.na],mu=tot.home.deaths.new.odesim[t.idx[no.na]], size=nb.disp.home.deaths,log = TRUE))
# 
#         ll.deaths=ll.home.deaths+ll.hosp.deaths
#         ll=ll+ll.deaths
#     }
# 

    


    #######################################################################
    ### additional likelihood for total current hosp / vent / icu cases ###
    #######################################################################

    # daily change in vent numbers
    delta.vent=c(tot.vent.curr[1],tot.vent.curr[-1]-tot.vent.curr[-n.days])
    delta.vent.odesim=c(tot.vent.curr.odesim[1],tot.vent.curr.odesim[-1]-tot.vent.curr.odesim[-n.days])
    # daily change in icu+vent numbers
    delta.icu=c(tot.icu.curr[1],tot.icu.curr[-1]-tot.icu.curr[-n.days])
    delta.icu.odesim=c(tot.icu.curr.odesim[1],tot.icu.curr.odesim[-1]-tot.icu.curr.odesim[-n.days])
    # daily change in icu (ignoring vent) numbers
    delta.icu.novent.odesim=c(tot.icu.curr.novent.odesim[1],tot.icu.curr.novent.odesim[-1]-tot.icu.curr.novent.odesim[-n.days])
    # daily change in icu+vent+hosp numbers
    delta.hosp=c(tot.hosp.curr[1],tot.hosp.curr[-1]-tot.hosp.curr[-n.days])
    delta.hosp.odesim=c(tot.hosp.curr.odesim[1],tot.hosp.curr.odesim[-1]-tot.hosp.curr.odesim[-n.days])
    # daily change in hosp (ignorming icu/vent) numbers
    delta.hosp.noicu.odesim=c(tot.hosp.curr.noicu.odesim[1],tot.hosp.curr.noicu.odesim[-1]-tot.hosp.curr.noicu.odesim[-n.days])
  
    icu.no.na=tot.icu.curr
    icu.no.na[tot.icu.curr>-1]=1
    hosp.no.na=tot.hosp.curr
    hosp.no.na[tot.hosp.curr>-1]=1
    
    if(lik.vent.curr){
      no.na.idx=which(!is.na(delta.vent))
      # Gaussian likelihood
      ll.vent=sum(dnorm(delta.vent[no.na.idx],mean=delta.vent.odesim[t.idx[no.na.idx]],sd=sqrt(s2.vent),log=TRUE))
      # add to current ll
      ll <- ll + ll.vent
    }else{
      ll.vent <- integer(0)
    }

    if(lik.icu.curr){
      # index for times with both icu and vent data
      idx.icu.vent=which(!is.na(delta.icu & delta.vent))
      idx.icu.only=which(!is.na(delta.icu) & is.na(delta.vent)) 
      # Gaussian likelihood 
      ll.icu.vent=0
      ll.icu.only=0
      
      if(length(idx.icu.vent)>0){
        ll.icu.vent=sum(dnorm(delta.icu[idx.icu.vent]-delta.vent[idx.icu.vent],mean=delta.icu.novent.odesim[t.idx[idx.icu.vent]],sd=sqrt(s2.icu),log=TRUE))
      }
      
      if(length(idx.icu.only)>0){
        ll.icu.only=sum(dnorm(delta.icu[idx.icu.only],mean=delta.icu.odesim[t.idx[idx.icu.only]],sd=sqrt(s2.icu+s2.vent),log=TRUE))
      }
      
      # add to current ll
      ll.icu <- ll.icu.vent + ll.icu.only
      ll <- ll + ll.icu
      
    } else {
      ll.icu <- integer(0)
    }

    if(lik.hosp.curr){
      no.na.idx=which(!is.na(tot.hosp.curr))
      # index for times with both hosp and icu or vent data
      idx.hosp.icu=which(!is.na(delta.hosp+delta.icu))
      idx.hosp.vent=which(!is.na(delta.hosp+delta.vent) & is.na(delta.icu))
      idx.hosp.only=which(!is.na(delta.hosp) & is.na(delta.vent) & is.na(delta.icu)) 
      # Gaussian likelihood 
      ll.hosp.vent=0
      ll.hosp.icu=0
      ll.hosp.only=0
      # case in which hosp and icu are observed
      if(length(idx.hosp.icu)>0){
        ll.hosp.icu=sum(dnorm(delta.hosp[idx.hosp.icu]-delta.icu[idx.hosp.icu],mean=delta.hosp.noicu.odesim[t.idx[idx.hosp.icu]],sd=sqrt(s2.hosp),log=TRUE))
      }
      # case in which hosp and vent, but not icu, are observed
      if(length(idx.hosp.vent)>0){
        ll.hosp.vent=sum(dnorm(delta.hosp[idx.hosp.vent]-delta.vent[idx.hosp.vent],mean=delta.hosp.odesim[t.idx[idx.hosp.vent]]-delta.vent.odesim[t.idx[idx.hosp.vent]],sd=sqrt(s2.hosp+s2.vent),log=TRUE))
      }
      # case in which hosp is only observation
      if(length(idx.hosp.only)>0){
        ll.hosp.only=sum(dnorm(delta.hosp[idx.hosp.only],mean=delta.hosp.odesim[t.idx[idx.hosp.only]],sd=sqrt(s2.hosp+s2.icu+s2.vent),log=TRUE))
      }
      # add to current ll
      ll.hosp <- ll.hosp.icu + ll.hosp.vent + ll.hosp.only
      ll <- ll + ll.hosp
    } else {
      ll.hosp=integer(0)
    }

    ##########################################
    ### likelihood for hospital discharges ###
    ##########################################
    ll.hosp.dis=NULL
    if(lik.hosp.discharges){
      discharge.rates <- tot.hosp.discharges.new.odesim
      # finding any na.values
      no.na <- which(!is.na(tot.hosp.discharges.new))
      # log-likelihood given trajectory sim:
      ll.hosp.dis <- sum(dnbinom(tot.hosp.discharges.new[no.na],mu=discharge.rates[t.idx[no.na]] , size=nb.disp.hosp.discharges,log = TRUE))
      ll <- ll + ll.hosp.dis
    }

  } # end of loop checking if params are above zero

  
    ### results      
  list(ll = ll,
       ll.new = ll.new,
       ll.hosp.new = ll.hosp.new,
       ll.hosp = ll.hosp,
       ll.vent = ll.vent,
       ll.icu = ll.icu,
       sympt.imputed = sympt.imputed,
       ll.deaths = ll.deaths,
       ll.hosp.dis = ll.hosp.dis)
}
