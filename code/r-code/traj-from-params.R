
traj.from.params <- function(beta,
                             params = NULL,
                             const.params = NULL,
                             non.odesim.params=NULL,
                             introday = NULL,
                             tf,
                             odepath,
                             loc="RI",
                             symp = NULL
                             ){
  # A function that takes beta values as input and calls the odesim program.
  #   Input:
  #       1) beta: vector of beta values used in the sim
  #               should be a vector of length tf-60, with one value per day
  #       2) tf: integer day to stop sim (1 = 20200101)
  #       3) t0=jan 1, 2020
  #   Output:
  #       Data frame with times column and DeltaJ (simulated cases) column


    ##
    ## remove non.odesim.params from params and const.params
    ##

    if(length(non.odesim.params)>0){
      
      idx.remove <- integer()
      
      for(k in 1:length(non.odesim.params)){
        idx.remove <- c(idx.remove,which(names(params)==non.odesim.params[k]))
      }
        
      if (length(idx.remove) > 0){
        params <- params[-idx.remove]
      } else {
        params <- NULL
      }
    
      # if(length(const.params)>0){
      #   
      #   idx.remove <- integer()
      #   for(k in 1:length(non.odesim.params)){
      #     idx.remove <- c(idx.remove, which(names(const.params) == non.odesim.params[k]))
      #   }
      #   
      #   if (length(idx.remove) > 0){
      #     const.params <- const.params[-idx.remove]
      #   } else {
      #     const.params <- NULL
      #   }
      # }
    }
    
    
    ## make introday backwards compatible
    if(length(which(c(names(params),names(const.params))=="introday"))+length(introday)==0){
        introday=55
    }
    if(length(which(c(names(params),names(const.params))=="introday"))+length(introday)>1){
        cat("WARNING: introday specified in two or more places.","\n")
    }
    
  ## create string to run odesim with parameters
    cmd1 <- paste(c(odepath,"odesim none -binary-output -tf ",as.character(tf)),collapse="")
    if(length(introday)==1){
        cmd1.5 <- paste(c(" -introday"), as.character(introday), collapse = " ")
        cmd2 <- paste(c(" -beta", as.character(beta)), collapse = " ")
        cmd <- paste(cmd1, cmd1.5, cmd2, sep="")
    } else {
      cmd2 <- paste(c(" -beta", as.character(beta)), collapse = " ")
      cmd <- paste(cmd1, cmd2, sep="")
    }
    
  ## add constants, if any
  if(length(const.params) > 0){
    const.names = names(const.params)
    for(k in 1:length(const.params)){
      cmd=paste(c(cmd, " -", const.names[k], " ", as.character(const.params[k])), collapse = "")
    }
  }
    ## add parameters
  if(length(params)>0){
    names=names(params)
    for(k in 1:length(params)){
      cmd=paste(c(cmd," -",names(params)[k]," ",as.character(params[k])),collapse = "")
    }
  }
  
  if (length(symp) > 0){
    # specify symp-frac (used when sf.choice == TRUE)
    cmd <- paste(cmd, symp, collapse = "")
  }
  
  # add location parameter
  cmd=paste(cmd,"-loc ",as.character(loc), "2>&1",collapse="")
  
  
  
  
  ncols  = 307 #Number of expected columns in output
  nlines = as.integer(365) #Upper limit on number of output lines
  #nlines = as.integer(tf - introday + 1) # Output starts on introday
  #nlines = as.integer(tf) # Output starts on 0
  
  
  # call the odesim program, with beta as input
  p = pipe(cmd, 'rb')
  bin_result = readBin(p, 'double', nlines*ncols)
  close.connection(p)

  ##  browser()
  
  #Reshaping 1D result to (*-by-ncols)
  dim(bin_result) = c(ncols, length(bin_result)/ncols) # Because dim fills row-first :(
  bin_result = aperm(bin_result, c(2,1))

  bin_result
}
