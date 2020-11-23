#!/usr/bin/env Rscript

# results-chains.R
# authors: Nathan Wikle
# last edited: 23 Nov 2020 
#
# Function (plot.chains) to generate trace plots from MCMC results.

plot.chains = function(out.folder, 
                       nrow.plots = 3, 
                       ncol.plots = 5, 
                       pdf.plot = TRUE,
                       pdf.name = "chains.pdf"){
  # Input: 
  #     out.folder: location of MCMC output (.Rdata)
  #     nrow.plots: number of rows of plots (default = 3)
  #     ncol.plots: number of columns of plots (default = 5)
  #     pdf.plot: whether to make a pdf of the results (default = TRUE)
  #     pdf.name: name of pdf file (default = "chains.pdf")
  # Output: 
  #   Trace plots for all parameters estimated via MCMC.

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
  spline=data[[1]]$spline
  df=data[[1]]$df
  ##
  beta.chains=array(NA,dim=c(n.mcmc*max.iter,n.beta,n.chains))
  ode.chains=array(NA,dim=c(n.mcmc*max.iter,n.ode.params,n.chains))
  rr.chains=array(NA,dim=c(n.mcmc*max.iter,n.rr.params,n.chains))
  lik.chains=array(NA,dim=c(n.mcmc*max.iter,n.lik.params,n.chains))
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
    }
  }
  
  ## str(beta.chains)
  ##
  if(pdf.plot){
    pdf(pdf.name,width=10,height=7)
  }
  par(mfrow=c(nrow.plots,ncol.plots))
  for(k in 1:n.beta){
    matplot(beta.chains[,k,],type="l",col=k,main=paste("beta",k))
  }
  for(k in 1:n.ode.params){
    matplot(ode.chains[,k,],type="l",col=k,main=colnames(data[[i]]$ode.params)[k])
  }
  for(k in 1:n.rr.params){
    matplot(rr.chains[,k,],type="l",col=k,main=paste("rr",k))
  }
  for(k in 1:n.lik.params){
    matplot(lik.chains[,k,],type="l",col=k,main=paste("disp",k))
  }
  if(pdf.plot){
    dev.off()
  }
}
