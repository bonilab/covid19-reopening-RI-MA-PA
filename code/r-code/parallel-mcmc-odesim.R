multichain.mcmc.odesim <- function(n.chains,
                                   n.cores,
                                   n.iter,
                                   resample=FALSE,
                                   parallel.type="doMC", ## "doMC" or "Foreach" parallelization
                                   inf.file.dir="./", ## location of mcmc.odesim.R / loglik.odesim.R / traj.from.params.R### remaining arguments go into mcmc.odesim
                                   n.mcmc.per.iter, ## number of mcmc samples per iteration
                                   save.file.name="out",
                                   betas.lst, ## starting values
                                   ode.params.lst, ## starting values
                                   lik.params.lst, ## starting values
                                   rr.params.lst, ## starting values
                                   s2.hosp.lst, #var# starting values
                                   s2.vent.lst, ## starting values
                                   s2.icu.lst, ## starting values
                                   Sigma.tune.start,
                                   var.tune.start,
                                   adapt = "ShabyWells",
                                   thin.rt=1,
                                   c0=1,
                                   c1=0.8,
                                   ...){

  ##browser()
    
    Sigma.t=Sigma.tune.start
    var.t=list()
    for(i in 1:n.chains){
        var.t[[i]]=var.tune.start
        
    }
    t.adpt=0
    
    
    if(parallel.type=="doMC"){
        ## start parallel
        require(doMC)
        registerDoMC(cores=n.cores)
    }
    if(parallel.type=="psock"){
### Parallelization using Foreach
        
                                        # set up Foreach (with number of processors = n.chs!)
        mp_type = "PSOCK"
        cl <- parallel::makeCluster(n.chains, type = mp_type)
        doParallel::registerDoParallel(cl)
        
        clusterCall(cl, function() {
            ## library for vectorized multinomial
            library(mc2d)
            ## library for multivariate normal distribution
            library(mvtnorm)
            ## library for spline-based expansions
            library(fda)
            ## library for truncated normal
            library(msm)
            ## functions for local use
            source(paste(inf.file.dir,"loglik.odesim.4.0.R",sep=""), chdir=TRUE)
            source(paste(inf.file.dir,"mcmc.odesim.2.0.R",sep=""))
            source(paste(inf.file.dir,"traj.from.params.R",sep=""))
            source(paste(inf.file.dir,"traj.process.R",sep=""))
            source(paste(inf.file.dir,"data.process.R",sep=""))
            source(paste(inf.file.dir,"plots.odesim.R",sep=""))
        })
        
    }
    
    i=1
    
    for(iter in 1:n.iter){
        cat(iter,"\n")
        out <- foreach(i=1:n.chains) %dopar% mcmc.odesim(n.mcmc = n.mcmc.per.iter,
                                                         beta.start=betas.lst[[i]],
                                                         ode.params.start=ode.params.lst[[i]],
                                                         report.rate.params.start=rr.params.lst[[i]],
                                                         lik.params.start=lik.params.lst[[i]],
                                                         s2.hosp.start=s2.hosp.lst[[i]],
                                                         s2.icu.start=s2.icu.lst[[i]],
                                                         s2.vent.start=s2.vent.lst[[i]],
                                                         Sigma.tune=Sigma.t,
                                                         var.tune=var.t[[i]],
                                                         adapt.type= adapt,
                                                         thin=thin.rt,
                                                         t.adapt.start=t.adpt,
                                        #c0=c0,
                                        #c1=c1,
                                                         ...)
        
        
        ## resampling
        if(resample){
            ll.vec=out[[1]]$loglik.final.iter
            for(k in 2:n.chains){
                ll.vec=c(ll.vec,out[[2]]$loglik.final.iter)
            }
            resample.probs=exp(ll.vec-max(ll.vec))
            resamp.idx=sample(1:n.chains,n.chains,replace=TRUE,prob=resample.probs)
        }else{
            resamp.idx=1:n.chains
        }

        last.idx=n.mcmc.per.iter/thin.rt
        for(k in 1:n.chains){
            betas.lst[[k]]=out[[resamp.idx[k]]]$beta[last.idx,]
            ode.params.lst[[k]]=out[[resamp.idx[k]]]$ode.params[last.idx,]
            lik.params.lst[[k]]=out[[resamp.idx[k]]]$lik.params[last.idx,]
            rr.params.lst[[k]]=out[[resamp.idx[k]]]$rr.params[last.idx,]
            s2.hosp.lst[[k]]=out[[resamp.idx[k]]]$s2.params[last.idx,1]
            s2.icu.lst[[k]]=out[[resamp.idx[k]]]$s2.params[last.idx,2]
            s2.vent.lst[[k]]=out[[resamp.idx[k]]]$s2.params[last.idx,3]
            var.t[[k]]=out[[resamp.idx[k]]]$var.tune
            ## Sigma.t[[k]]=out[[resamp.idx[k]]]$Sigma.tune
            
        }
        t.adpt=out[[1]]$t.adapt.end
        ## saving output
        iter.char=as.character(iter)
        num.leading.zeros=nchar(as.character(n.iter))-nchar(iter.char)
        if(num.leading.zeros==0){
          leading.zeros=""
        }
        if(num.leading.zeros==1){
          leading.zeros="0"
        }
        if(num.leading.zeros==2){
          leading.zeros="00"
        }
        if(num.leading.zeros==3){
          leading.zeros="000"
        }
        save(out,file=paste(save.file.name,"-",leading.zeros,as.character(iter),".Rdata",sep=""))
    }
    out
}


##### old code, saved just in case


          ## beta.start=betas.start[[i]],
        ##                                                  ode.params.start=ode.params.start[[i]],
        ##                                                  report.rate.params.start=report.rate.params.start[[i]],
        ##                                                  lik.params.start=lik.params.start[[i]],
        ##                                                  s2.hosp.start=s2.hosp.start[[i]],
        ##                                                  s2.icu.start=s2.icu.start[[i]],
        ##                                                  s2.vent.start=s2.vent.start[[i]],
        ##                                                  n.mcmc=n.mcmc,
        ##                                                  Sigma.tune=Sigma.tune,
        ##                                                  var.tune=var.tune,
        ##                                                  adapt.iter=sqrt(2),
        ##                                                  adapt.type="ShabyWells",
        ##                                                  c0=c0,
        ##                                                  c1=c1,
        ##                                                  ...)


    ##         Sigma.hat=out[[1]]$Sigma.hat
##         for(ii in 1:5){
##             Sigma.hat[[ii]]=Sigma.hat[[ii]]/n.chains
##         }
##         accept.rate=out[[1]]$accept.rate/n.chains
##         if(n.chains>1){
##             for(k in 2:n.chains){
##                 for(ii in 1:5)
##                     Sigma.hat[[ii]]=Sigma.hat[[ii]]+out[[k]]$Sigma.hat[[ii]]/n.chains
##             }
##             accept.rate=accept.rate+out[[k]]$accept.rate/n.chains
##         }
        
##         ##
##         ## tuning
##         ##
##         ## Using Shaby and Wells adaptive scheme for RWM...
##         if(adapt.type=="ShabyWells"){
##             for(ii in 1:5){
##                 ## default constants
##                 ## c0 <- 1; c1 <- 0.8;
##                 r.opt <- 0.234
##                 ## update based on empirical covar of most recent block
##                 r.hat <- accept.rate[[ii]]
##                 gamma.1 <- 1 / (iter)^c1
##                 gamma.2 <- c0 * gamma.1
##                 var.tune[[ii]] <- exp(log(var.tune[[ii]]) + gamma.2 * (r.hat - r.opt))
##                 Sigma.tune[[ii]] <- Sigma.tune[[ii]] + gamma.1 * (Sigma.hat[[ii]] - Sigma.tune[[ii]])
##             }
##         }
##         if(adapt.type=="Covar"){
##             ## Sigma.hat <- cov( cbind(beta.save,ode.params.save,lik.params.save)[1:iter,])
##             for(ii in 1:5){
##                 Sigma.tune[[ii]]=2.4^2/nrow(Sigma.hat[[ii]])*Sigma.hat[[ii]]
##             }
##         }
        
##         ## resampling
##         if(resample){
##             ll.vec=out[[1]]$loglik.final.iter
##             for(k in 2:n.chains){
##                 ll.vec=c(ll.vec,out[[2]]$loglik.final.iter)
##             }
##             resample.probs=exp(ll.vec-max(ll.vec))
##             resamp.idx=sample(1:n.chains,n.chains,replace=TRUE,prob=resample.probs)
##         }else{
##             resamp.idx=1:n.chains
##         }
##         for(k in 1:n.chains){
##             betas.start[[k]]=out[[resamp.idx[k]]]$beta[n.mcmc,]
##             ode.params.start[[k]]=out[[resamp.idx[k]]]$ode.params[n.mcmc,]
##             lik.params.start[[k]]=out[[resamp.idx[k]]]$lik.params[n.mcmc,]
##             rr.params.start[[k]]=out[[resamp.idx[k]]]$rr.params[n.mcmc,]
##             s2.hosp.start[[k]]=out[[resamp.idx[k]]]$s2.params[n.mcmc,1]
##             s2.icu.start[[k]]=out[[resamp.idx[k]]]$s2.params[n.mcmc,2]
##             s2.vent.start[[k]]=out[[resamp.idx[k]]]$s2.params[n.mcmc,3]
##         }
        
##         ## saving output
##         save(out,file=paste(save.file.name,"-",as.character(iter),".Rdata",sep=""))
##     }

##     out
## }

