#!/usr/bin/env Rscript

##
## example use
##
## plot.chains("./",nrow.plots=3,ncol.plots=5,pdf.plot=TRUE,pdf.name="chains.pdf")
## results.plots.and.params("./",name="MAtoday",loc="MA",odesim.version="v5",burnin=20000)

results.plots.and.params = function(out.folder,
                                    loc="RI",
                                    odesim.version="v5",
                                    burnin=1,
                                    name="output",
                                    readme=TRUE,
                                    odepath="../../",
                                    non.odesim.params=NULL,
                                    const.params=NULL,
                                    lik.hosp.discharges=FALSE,
                                    active.surv=FALSE){
  ##
  ## Function to produce uniform plots of output
  ## 
  
  files.list=list.files(out.folder,"*.Rdata")
  files.list <- files.list[order(nchar(files.list), files.list)]
  max.iter=length(files.list)
  ##
  load(paste(out.folder,files.list[1],sep=""))
  data=out
  n.chains=length(data)
  n.mcmc=nrow(data[[1]]$beta)
  n.beta=ncol(data[[1]]$beta)
  n.ode.params=ncol(data[[1]]$ode.params)
  n.rr.params=ncol(data[[1]]$rr.params)
  n.lik.params=ncol(data[[1]]$lik.params)
  n.s2.params=ncol(data[[1]]$s2.params)
  spline=data[[1]]$spline
  df=data[[1]]$df
  ##
  beta.chains=array(NA,dim=c(n.mcmc*max.iter,n.beta,n.chains))
  ode.chains=array(NA,dim=c(n.mcmc*max.iter,n.ode.params,n.chains))
  rr.chains=array(NA,dim=c(n.mcmc*max.iter,n.rr.params,n.chains))
  lik.chains=array(NA,dim=c(n.mcmc*max.iter,n.lik.params,n.chains))
  s2.chains=array(NA,dim=c(n.mcmc*max.iter,n.s2.params,n.chains))
  loglik.chains=matrix(NA,nrow=n.mcmc*max.iter,ncol=n.chains)
  for(iter in 1:max.iter){
    ## load in data and call it "data"
    load(paste(out.folder,files.list[iter],sep=""))
    data=out
    ## get betas and ode.params
    for(i in 1:n.chains){
      beta.chains[(iter-1)*n.mcmc+(1:n.mcmc),,i]=data[[i]]$beta
      ode.chains[(iter-1)*n.mcmc+(1:n.mcmc),,i]=data[[i]]$ode.params
      rr.chains[(iter-1)*n.mcmc+(1:n.mcmc),,i]=data[[i]]$rr.params
      lik.chains[(iter-1)*n.mcmc+(1:n.mcmc),,i]=data[[i]]$lik.params
      s2.chains[(iter-1)*n.mcmc+(1:n.mcmc),,i]=data[[i]]$s2.params
      loglik.chains[(iter-1)*n.mcmc+(1:n.mcmc),i]=data[[i]]$loglik
    }
  }
  
  
  
  library(fda)
  ## attach(data[[1]])
  df=data[[1]]$df
  days=df$daynum
  Z.beta=eval.basis(61:max(days),data[[1]]$spline.beta)
  
  if (is.matrix(data[[1]]$spline.rr)){
    Z.rr <- data[[1]]$spline.rr
    Z.rr <- Z.rr[-nrow(Z.rr),]
  } else {
    Z.rr <- eval.basis(61:max(days),data[[1]]$spline.rr)
  }
  #cumul.hosp.rr=data[[1]]$hosp.report.rate
  
  ## start pdf
  pdf(paste(name,".pdf",sep=""),width=10,height=7)
  
  ##
  ##
  ## 1. Estimated Contact rate parameter over time.
  ##
  ##
  
 ## browser()
  
  ## We first plot the estimated contact rate parameter (beta) over time.  This shows the changing population-level contact rate in the state.
  
  ## mean beta and 95% CIs
  par(mfrow=c(1,1))
  mn=Z.beta%*%apply(beta.chains[-(1:burnin),,], 2, mean)
  uq=Z.beta%*%apply(beta.chains[-(1:burnin),,], 2, quantile, 0.975)
  lq=Z.beta%*%apply(beta.chains[-(1:burnin),,], 2, quantile, 0.025)
  q75=Z.beta%*%apply(beta.chains[-(1:burnin),,], 2, quantile, 0.75)
  q25=Z.beta%*%apply(beta.chains[-(1:burnin),,], 2, quantile, 0.25)
  ## plots
  plot(days,mn, type = "l", lwd = 5, col="blue", ylim=c(min(lq),max(uq)),
       main = "Contact Rate Parameter Over Time", ylab = "beta",xlab="Julian Day (61=Mar 1, 92=Apr 1)")
  points(days,uq,type = "l", col = grey(.5),lty=3,lwd=3)
  points(days,lq,type = "l", col = grey(.5),lty=3,lwd=3)
  points(days,q75,type = "l", col = grey(.5),lty=1,lwd=3)
  points(days,q25,type = "l", col = grey(.5),lty=1,lwd=3)
  legend("topright",legend=c("Posterior Mean","Posterior 95% CI","Posterior 50% CI"),col=c("blue",grey(.5),grey(.5)),lty=c(2,3,1),lwd=3)
  
 
  ##  browser()
  
  ## mean report.rate and 95% CIs
  par(mfrow=c(1,1))
  if(dim(rr.chains)[2]>1){
      mn.rr=Z.rr%*%apply(rr.chains[-(1:burnin),,], 2, mean)
      uq=Z.rr%*%apply(rr.chains[-(1:burnin),,], 2, quantile, 0.975)
      lq=Z.rr%*%apply(rr.chains[-(1:burnin),,], 2, quantile, 0.025)
      q75=Z.rr%*%apply(rr.chains[-(1:burnin),,], 2, quantile, 0.75)
      q25=Z.rr%*%apply(rr.chains[-(1:burnin),,], 2, quantile, 0.25)
  }else{
      mn.rr=Z.rr*mean(rr.chains[-(1:burnin),,])
      uq=Z.rr%*%quantile(rr.chains[-(1:burnin),,], 0.975)
      lq=Z.rr%*%quantile(rr.chains[-(1:burnin),,], 0.025)
      q75=Z.rr%*%quantile(rr.chains[-(1:burnin),,], 0.75)
      q25=Z.rr%*%quantile(rr.chains[-(1:burnin),,], 0.25)
  }
    
  ## plots
  plot(days, mn.rr, type = "l", lwd = 5, col="blue", ylim=c(min(lq), max(uq)),
       main = "Reporting Rate Over Time", ylab = "beta",xlab="Julian Day (61=Mar 1, 92=Apr 1)")
  points(days,uq,type = "l", col = grey(.5),lty=3,lwd=3)
  points(days,lq,type = "l", col = grey(.5),lty=3,lwd=3)
  points(days,q75,type = "l", col = grey(.5),lty=1,lwd=3)
  points(days,q25,type = "l", col = grey(.5),lty=1,lwd=3)
  legend("topright",legend=c("Posterior Mean","Posterior 95% CI","Posterior 50% CI"),col=c("blue",grey(.5),grey(.5)),lty=c(2,3,1),lwd=3)
  
  
  
  
  ##
  ##
  ## Hospitalization params
  ##
  ##
  
  n.ode.params=ncol(data[[1]]$ode.params)
  par(mfrow=c(3,5))
  for(i in 1:n.ode.params){
    hist(ode.chains[-(1:burnin),i,],main=colnames(data[[1]]$ode.params)[i],col=i)
  }
  par(mfrow=c(1,1))
  
  ##
  ##
  ## Hospitalization variance params
  ##
  ##
  
  par(mfrow=c(2,2))
  for(i in 1:n.s2.params){
    hist(s2.chains[-(1:burnin),i,],main=paste("variance param ",i,sep=""),col=i)
  }
  par(mfrow=c(1,1))
  
  ##
  ##
  ## Hospitalization variance params
  ##
  ##
  
  par(mfrow=c(2,2))
  for(i in 1:n.lik.params){
    hist(lik.chains[-(1:burnin),i,],main=paste("NB disp param ",i,sep=""),col=i)
  }
  par(mfrow=c(1,1))
  
  
  ##
  ##
  ## 3 Estimated Epidemic Trajectories
  ##
  ##
  
  
  ##browser()
  
  ## read in data and saved trajectories after burnin
  dp=data.process(df,loc=loc)
  traj.list=list()
  idx.remove=ceiling(burnin/n.mcmc)
  loglik.save=matrix(NA,nrow=max.iter-idx.remove+1,ncol=n.chains)
  for(i in (idx.remove:max.iter)){
    load(paste(out.folder,files.list[i],sep=""))
    for(k in 1:n.chains){
      traj.list[[length(traj.list)+1]]=traj.process(out[[k]]$traj,loc=loc,odesim.ver=odesim.version)
      loglik.save[i-idx.remove+1,k]=out[[k]]$loglik.final.iter
    }
  }
  n.traj=length(traj.list)
  
  ## prep for plotting
  list2env(dp,globalenv())
  ## linking times of data and odesim output
  tmin <- min(days)
  tmax <- max(days)
  idx.min <- which(round(traj.list[[1]]$days.odesim) == tmin)
  idx.max <- which(round(traj.list[[1]]$days.odesim) == tmax)
  t.idx <- idx.min:idx.max
  days.odesim=traj.list[[1]]$days.odesim
  ##
  ## plot of new symptomatics
  ##
  
  plot(days,tot.sympt.new/mn.rr,col="blue",pch=20,cex=2,main="Daily Total New Reported Cases")
  for(i in 1:n.traj){
    points(days.odesim,traj.list[[i]]$tot.sympt.new.odesim,col="red",type="l")
  }
  points(days,tot.sympt.new/mn.rr,col="blue",pch=20,cex=2)
  
  ##
  ## plot of cumulative symptomatics
  ##
  tsn=tot.sympt.new
  tsn[is.na(tsn)]=0
  tsn=tsn/mn.rr
  plot(days,cumsum(tsn),col="blue",pch=20,cex=2,main="Cumulative Reported Cases")
  for(i in 1:n.traj){
    points(days.odesim,cumsum(traj.list[[i]]$tot.sympt.new.odesim),col="red",type="l")
  }
  points(days,cumsum(tsn),col="blue",pch=20,cex=2)
  
  
  
  ##
  ## plot of new hospitalizations
  ##
  
  if(is.element("hosp.report.rate",colnames(data[[1]]$ode.params))){
    idx.h=which(colnames(data[[1]]$ode.params)=="hosp.report.rate")
    cumul.hosp.rr=mean(ode.chains[-(1:burnin),idx.h,])
  } else {
    cumul.hosp.rr <- 1
  }
  
  
  if(sum(tot.hosp.new,na.rm=TRUE)>0){
    plot(days,tot.hosp.new/cumul.hosp.rr,col="blue",pch=20,cex=2,main="Daily Total New Hospitalized Cases")
    for(i in 1:n.traj){
      points(days.odesim,traj.list[[i]]$tot.hosp.new.odesim,col="red",type="l")
    }
    points(days,tot.hosp.new/cumul.hosp.rr,col="blue",pch=20,cex=2,main="Daily Total New Hospitalized Cases")
  }
  ##
  ## plot of current total hospitalized individuals
  ##
  if(sum(tot.hosp.curr,na.rm=TRUE)>0){
    plot(days,tot.hosp.curr,col="blue",pch=20,cex=2,main="Daily Current Hospitalized Cases")
    for(i in 1:n.traj){
      points(days.odesim,traj.list[[i]]$tot.hosp.curr.odesim,col="red",type="l")
    }
    points(days,tot.hosp.curr,col="blue",pch=20,cex=2)
  }
  ##
  ## plot of current total icu individuals
  ##
  if(sum(tot.icu.curr,na.rm=TRUE)>0){
    plot(days,tot.icu.curr,col="blue",pch=20,cex=2,main="Daily Current ICU Cases")
    for(i in 1:n.traj){
      points(days.odesim,traj.list[[i]]$tot.icu.curr.odesim,col="red",type="l")
    }
    points(days,tot.icu.curr,col="blue",pch=20,cex=2)
  }
  ##
  ## plot of current total individuals on ventilators
  ##
  if(sum(tot.vent.curr,na.rm=TRUE)>0){
    plot(days,tot.vent.curr,col="blue",pch=20,cex=2,main="Daily Current Cases on Ventilator")
    for(i in 1:n.traj){
      points(days.odesim,traj.list[[i]]$tot.vent.curr.odesim,col="red",type="l")
    }
    points(days,tot.vent.curr,col="blue",pch=20,cex=2)
  }
  ##
  ## plot of daily new deaths
  ##
  if(sum(tot.deaths.new,na.rm=TRUE)>0){
    plot(days,tot.deaths.new,col="blue",pch=20,cex=2,main="Daily Total New Deaths")
    for(i in 1:n.traj){
      points(days.odesim,traj.list[[i]]$tot.deaths.new.odesim,col="red",type="l")
    }
    points(days,tot.deaths.new,col="blue",pch=20,cex=2)
  }
  
  ##
  ## plot of cumulative new deaths
  ##
  if(sum(tot.deaths.new,na.rm=TRUE)>0){
    tdn=tot.deaths.new
    tdn[is.na(tdn)]=0
    plot(days,cumsum(tdn),col="blue",pch=20,cex=2,main=" Total Cumulative Deaths")
    for(i in 1:n.traj){
      points(days.odesim,cumsum(traj.list[[i]]$tot.deaths.new.odesim),col="red",type="l")
    }
    points(days,cumsum(tdn),col="blue",pch=20,cex=2)
  }
  
##  attach(data)
  
  plot(0,0,type="n",xlim=c(0,100),ylim=c(0,100))
  legend("topright",legend=c(paste("lik.tot = ",data[[1]]$lik.tot,sep=""),
                            paste("lik.age = ",data[[1]]$lik.age,sep=""),
                            paste("lik.hosp.new = ",data[[1]]$lik.hosp.new,sep=""),
                            paste("lik.hosp.curr = ",data[[1]]$lik.hosp.curr,sep=""),
                            paste("lik.icu.curr = ",data[[1]]$lik.icu.curr,sep=""),
                            paste("lik.vent.curr = ",data[[1]]$lik.vent.curr,sep=""),
                            paste("lik.tot.deaths = ",data[[1]]$lik.tot.deaths,sep=""),
                            paste("lik.hosp.deaths = ",data[[1]]$lik.hosp.deaths,sep=""),
                            paste("lik.home.deaths = ",data[[1]]$lik.home.deaths,sep="")))
         
 legend("topleft",legend=c(paste("State = ",loc,sep=""),
                            paste("Max Day of Data = ",max(days),sep=""),
                            paste("Day Plotted = ",Sys.Date(),sep="")))
                            

 legend("bottomleft"
        ,legend=c(paste("beta.hat=",paste(as.character(round(apply(beta.chains[-(1:burnin),,], 2, mean),2)),sep=" ",collapse=" ")),
                  paste("rr.hat=",paste(as.character(round(apply(rr.chains[-(1:burnin),,], 2, mean),2)),sep=" ",collapse=" ")),
                  paste("ode.params.hat=",paste(as.character(round(apply(ode.chains[-(1:burnin),,], 2, mean),2)),sep=" ",collapse=" ")),
                  paste("ode.params=",paste(colnames(data[[1]]$ode.params),collapse=", "))
                  ),cex=c(1,1,1,.75)
        )
 
 
 ## plotting loglikelihood
 ## matplot(loglik.save,type="l",main="Loglikelihood of MCMC iterations",xlab="MCMC iteration",ylab="logliklihood")
 hist(as.numeric(loglik.chains),col="yellow",main="Loglikelihood of MCMC iterations",xlab="loglik")
  
  dev.off()
  
  ##
  ## Make output files
  ##
  
  M=1000
  idx=round(seq(burnin+1,dim(beta.chains)[1],length.out=M/n.chains))
  beta.idx=beta.chains[idx,,]
  beta.plot=beta.idx[,,1]
  for(i in 2:n.chains){
    beta.plot=rbind(beta.plot,beta.idx[,,i])
  }
  daily.beta=t(Z.beta%*%t(beta.plot))
  ## each row of daily.beta.idx is a daily beta value from day=61 to the last day of data
  
  ##
    params.idx=ode.chains[idx,,1]
    for(i in 2:n.chains){
        params.idx=rbind(params.idx,ode.chains[idx,,i])
    }
    rr.idx=rr.chains[idx,,1]
    if(dim(rr.chains)[2]==1){
        rr.idx=matrix(rr.idx,ncol=1,nrow=length(rr.idx))
        for(i in 2:n.chains){
            
            rr.idx=rbind(rr.idx,matrix(rr.chains[idx,,i],ncol=1))
        }
    }else{
        for(i in 2:n.chains){
            rr.idx=rbind(rr.idx,rr.chains[idx,,i])
        }
    }
  
  colnames(rr.idx)=paste("rr.",as.character(1:ncol(rr.idx)),sep="")
  colnames(params.idx)=colnames(data[[1]]$ode.params)
  colnames(daily.beta)=61:(60+ncol(daily.beta))
  
  write.table(rbind(61:(60+ncol(daily.beta)),daily.beta),file=paste(name,".daily.betas-day-",max(days),".csv",sep=""),sep=",",row.names=FALSE,col.names=FALSE)
  
  write.table(beta.plot,file=paste(name,".spline.coeff.betas-day-",max(days),".csv",sep=""),sep=",",row.names=FALSE,col.names=FALSE)
  
  write.table(cbind(params.idx,rr.idx),file=paste(name,".ode.params-day-",max(days),".csv",sep=""),sep=",",row.names=FALSE,col.names=TRUE)
  
  ##
  ## DIC calculations
  ##
  
  ## browser()
  ## get posterior mean values
  daily.beta.hat=mn
  daily.rr.hat=mn.rr
  ode.params.hat=apply(ode.chains[-(1:burnin),,],2,mean)
  names(ode.params.hat)=colnames(data[[1]]$ode.params)
  lik.params.hat=apply(lik.chains[-(1:burnin),,],2,mean)
  s2.params.hat=apply(s2.chains[-(1:burnin),,],2,mean)
  ## get loglikelihood of posterior mean
  if(length(which(names(data[[1]])=="introday"))==1){
    introday=data[[1]]$introday
  }else{
    introday=60
  }
    ## get non.odesim.params (backwards compatibility)
    if(length(which(names(data[[1]])=="non.odesim.params"))>0){
        n.o.p=data[[1]]$non.odesim.params
    }else{
        n.o.p=non.odesim.params
    }
    ## get const.params (backwards compatibility)
    if(length(which(names(data[[1]])=="const.params"))>0){
        c.p=data[[1]]$const.params
    }else{
        c.p=const.params
    }
  
  # end.day
  end.day <- max(data[[1]]$df$daynum, na.rm = T) + 1
  
  traj.hat <- traj.from.params(beta=as.numeric(daily.beta.hat), 
                               params=ode.params.hat, 
                               const.params = c.p,
                               non.odesim.params=n.o.p,
                               introday = introday,
                               tf = end.day, 
                               odepath=odepath,
                               loc=loc,
                               symp = NULL)


    # extra parameters (like hospitalization reporting rate )
  extra.params <- NULL
  extra.const.params <- NULL
  extra.params.fitted.idx <- integer()
  extra.params.const.idx <- integer()
    
  if(length(non.odesim.params)>0){
    
    for(k in 1:length(non.odesim.params)){
      extra.params.fitted.idx <- c(extra.params.fitted.idx,
        which(names(ode.params.hat) == non.odesim.params[k]))
    }
      
    if (length(extra.params.fitted.idx) > 0){
      extra.params <- ode.params.hat[extra.params.fitted.idx]
    }
    
    for(k in 1:length(non.odesim.params)){
      extra.params.const.idx <- c(extra.params.const.idx, 
        which(names(const.params) == non.odesim.params[k]))
    }
    
    if (length(extra.params.const.idx) > 0){
      extra.const.params <- const.params[extra.params.const.idx]
    }
  }

    
  ### create a matrix P which can be used to calculate Poisson rates after delay
  ###   note: this is no longer used for any state

    p.vecs <- list(p1 = c(1, rep(0,7)),
                    p2 = c(1, rep(0,7)),
                    p3 = c(1, rep(0,7)),
                    p4 = c(1, rep(0,7)))

    P <- matrix(0, nrow = end.day + max(sapply(p.vecs, length)), ncol = end.day+1)
  colnames(P) <- paste("Day", 1:(end.day+1), sep = " ")
  
  for(j in 1:ncol(P)){
    if (j < 74){
      ## period 1: beginning - March 13
      P[j:(j + length(p.vecs[[1]]) - 1), j] <- p.vecs[[1]]
    } else if ((j >= 74) & (j < 84)){
      ## period 2: March 14 - March 23
      P[j:(j + length(p.vecs[[2]]) - 1), j] <- p.vecs[[2]]
    } else if ((j >= 84) & (j < 88)){
      ## period 3: March 24 - March 27
      P[j:(j + length(p.vecs[[3]]) - 1), j] <- p.vecs[[3]]
    } else {
      ## period 4: March 28 - present
      P[j:(j + length(p.vecs[[4]]) - 1), j] <- p.vecs[[4]]
    }
  }


     ##browser()

    if(is.null(data[[1]]$lik.hosp.discharges)){
        data[[1]]$lik.hosp.discharges=lik.hosp.discharges
    }
    if(is.null(data[[1]]$active.surv)){
        data[[1]]$active.surv=active.surv
    }
    
    loglik.hat=loglik.odesim(traj.hat,
                       df,
                       dp=dp,
                       odesim.ver=odesim.version,
                       P=P,
                       loc=loc,
                       report.rate=c(daily.rr.hat,daily.rr.hat[length(daily.rr.hat)]),
                       nb.disp.params=lik.params.hat,
                       lik.tot=data[[1]]$lik.tot,
                       lik.age=data[[1]]$lik.age,
                       lik.hosp.new=data[[1]]$lik.hosp.new,
                       lik.hosp.curr=data[[1]]$lik.hosp.curr,
                       lik.icu.curr=data[[1]]$lik.icu.curr,
                       lik.vent.curr=data[[1]]$lik.vent.curr,
                       lik.tot.deaths=data[[1]]$lik.tot.deaths,
                       lik.home.deaths=data[[1]]$lik.home.deaths,
                       lik.hosp.deaths=data[[1]]$lik.hosp.deaths,
                       lik.age.deaths=data[[1]]$lik.age.deaths,
                       lik.hosp.discharges=data[[1]]$lik.hosp.discharges,
                       active.surv=data[[1]]$active.surv,
                       p.asympt=data[[1]]$p.asympt,
                       total.size.constraint=FALSE,
                       s2.hosp=s2.params.hat[1],
                       s2.icu=s2.params.hat[2],
                       s2.vent=s2.params.hat[3],
                       extra.params=extra.params, ### if cumul. hosp. reporting rate is fitted 
                       extra.const.params=extra.const.params) ### if cumul. hosp. reporting rate is constant

    ## calculate DIC
    Dhat=-2*loglik.hat$ll
    Dbar=-2*mean(loglik.chains,na.rm=TRUE)
    pD=Dbar-Dhat
    DIC=pD+Dbar
    
  ##
  ## Make readme file
  ##
  
  ## browser()
  
  if(readme){
    mcmc1=(out[[1]])
    
    tab=data.frame(input=integer(),value=integer())
    f=0
    tab[f+1,1]="Location"
    tab[f+1,2]=mcmc1$loc
    tab[f+2,1]="Last Day of Data Used"
    tab[f+2,2]=max(mcmc1$df$daynum)
    tab[f+3,1]="Day of MCMC Run"
    tab[f+3,2]=as.character(mcmc1$today)
    tab[f+4,]=" "
    
    j=nrow(tab)
    tab[j+1,]=" "
    tab[j+2,1]="Number of MCMC Chains"
    tab[j+2,2]=length(out)
    tab[j+3,1]="Number of Iterations Per Chain"
    tab[j+3,2]=n.mcmc*max.iter*mcmc1$thin
    tab[j+4,1]="Burnin (values before this discarded)"
    tab[j+4,2]=burnin
    tab[j+5,1]="Adaptive Tuning Type"
    tab[j+5,2]=mcmc1$adapt.type
    
    g=nrow(tab)
    tab[g+1,1]="lik.tot"
    tab[g+1,2]=mcmc1$lik.tot
    tab[g+2,1]="lik.age"
    tab[g+2,2]=mcmc1$lik.age
    tab[g+3,1]="lik.hosp.new"
    tab[g+3,2]=mcmc1$lik.hosp.new
    tab[g+4,1]="lik.hosp.curr"
    tab[g+4,2]=mcmc1$lik.hosp.curr
    tab[g+5,1]="lik.icu.curr"
    tab[g+5,2]=mcmc1$lik.icu.curr
    tab[g+6,1]="lik.vent.curr"
    tab[g+6,2]=mcmc1$lik.vent.curr
    tab[g+7,1]="lik.tot.deaths"
    tab[g+7,2]=mcmc1$lik.tot.deaths
    tab[g+8,1]="lik.home.deaths"
    tab[g+8,2]=mcmc1$lik.home.deaths
    tab[g+9,1]="lik.age.deaths"
    tab[g+9,2]=mcmc1$lik.age.deaths
    tab[g+10,1]="lik.hosp.deaths"
    tab[g+10,2]=mcmc1$lik.hosp.deaths
    #tab[g+11,1]="lik.hosp.discharges"
    #tab[g+11,2]=mcmc1$lik.hosp.discharges
    tab[g+11,]=" "
    tab[g+12,1]="odesim version"
    tab[g+12,2]=mcmc1$odesim.ver
    tab[g+13,]=" "
    tab[g+14,1]="Mean Log-Likelihood"
    tab[g+14,2]=mean(as.numeric(loglik.save))
    tab[g+15,1]="SD of Log-Likelihood"
    tab[g+15,2]=sqrt(var(as.numeric(loglik.save)))
    tab[g+16,]=" "
    tab[g+17,1]="Dbar"
    tab[g+17,2]=Dbar
    tab[g+18,1]="Dhat"
    tab[g+18,2]=Dhat
    tab[g+19,1]="DIC"
    tab[g+19,2]=DIC
    tab[g+20,]=" "
                    
                     
    
    h=nrow(tab)-1
    # tab[h,]=" "
    # tab[h+1,1]="Params In Model"
    # tab[h+1,2]=" "
    tab[h+2,]=" "
    tab[h+3,1]="Contact Rate Betas (spline length)"
    tab[h+3,2]=mcmc1$spline.beta$nbasis
    tab[h+4,1]="Contact Rate Betas (min day of spline)"
    tab[h+4,2]=mcmc1$spline.beta$rangeval[1]
    tab[h+5,1]="Contact Rate Betas (max day of spline)"
    tab[h+5,2]=mcmc1$spline.beta$rangeval[2]
    tab[h+6,]=" "
    tab[h+7,1]="Reporting Rate (spline length)"
    
    if (is.matrix(mcmc1$spline.rr)){
      tab[h+7,2]=ncol(mcmc1$spline.rr)
      tab[h+8,1]="Reporting Rate (min day of spline)"
      tab[h+8,2]=61
      tab[h+9,1]="Reporting Rate (max day of spline)"
      tab[h+9,2]= max(mcmc1$df$daynum, na.rm = T) + 1
      tab[h+10,] = " "
    } else {
      tab[h+7,2]=mcmc1$spline.rr$nbasis
      tab[h+8,1]="Reporting Rate (min day of spline)"
      tab[h+8,2]=mcmc1$spline.rr$rangeval[1]
      tab[h+9,1]="Reporting Rate (max day of spline)"
      tab[h+9,2]=mcmc1$spline.rr$rangeval[2]
      tab[h+10,] = " "
    }
    ##tab[h+11,1]="Cumulative Hosp. Reporting Rate"
    ##tab[h+11,2]=mcmc1$hosp.report.rate
      
    tab$max=" "
    i=nrow(tab)
    tab[i+1,]=" "
    tab[i+2,1]="Included ODESIM Params"
    tab[i+2,2]="Prior Min"
    tab[i+2,3]="Prior Max"
    n.params=ncol(mcmc1$ode.params)
    tab[i+2+(1:n.params),]=cbind(colnames(mcmc1$ode.params),mcmc1$ode.params.prior.min,mcmc1$ode.params.prior.max)
    
    if(length(mcmc1$const.params) > 0){
      tab[i+2+n.params+1,]=" "
      tab[i+2+n.params+2,]= c("Constant ODESIM Params", "value", " ")
      n.const.params=length(mcmc1$const.params)
      tab[i+2+n.params+2+(1:n.const.params),] = cbind(names(mcmc1$const.params), mcmc1$const.params, rep(" ", n.const.params))
      if (length(mcmc1$introday) > 0){
        tab[i + 2 + n.params + 2 + n.const.params + 1, ] = c("introday", mcmc1$introday, " ")
      }
    } else {
      if (length(mcmc1$introday) > 0){
        tab[i+2+n.params+1,]=" "
        tab[i+2+n.params+2,]= c("Constant ODESIM Params", "value", " ")
        tab[i + 2 + n.params + 3, ] = c("introday", mcmc1$introday, " ")
      }
    }
    
    
    write.table(tab,file=paste(name,"-README",".txt",sep=""),sep="\t",row.names=FALSE,col.names=FALSE,quote=FALSE)
  }
  
}


