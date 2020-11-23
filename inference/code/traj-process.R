#!/usr/bin/env Rscript

# traj-process.R
# authors: Ephraim Hanks, Nathan Wikle, Emily Strong
# last edited: 20 Nov 2020 
#
# This function (traj.process) takes output from odesim and returns a more 
#   interpretable list of outputs, which better match the available data.

traj.process <- function(traj, loc = "RI", odesim.version = "v5"){
  # Input:
  #   traj: odesim trajectories generated from 'traj.from.params'
  #   loc: US state ("RI" (default), "MA", "PA")
  #   odesim.version: version of odesim used for trajectories (default = "v5")
  # Output:
  #   Named list with the following objects:
  #     days, tot.sympt.new, tot.hosp.new, hosp.new, deaths.new, hosp.deaths.new, 
  #     home.deaths.new, tot.hosp.curr.noicu, tot.icu.curr.novent, tot.vent.curr, 
  #     tot.hosp.curr, tot.icu.curr, tot.deaths.new, tot.hosp.deaths.new,
  #     tot.home.deaths.new, tot.deaths.cum, deaths.cum, tot.hosp.discharges.new    


    if(odesim.version=="v4"){
        sympt.cum.idx=272:280
        hosp.cum.idx=281:289
    }
    if(odesim.version=="v5"){
        sympt.cum.idx=290:298
        hosp.cum.idx=299:307
        home.deaths.idx=254:262
        hosp.deaths.idx=263:271
    }
    hosp.curr.idx=c(137:172,245:253)
    icu.curr.idx=c(173:181,236:244)
    vent.curr.idx=182:235
    nr = nrow(traj)
    traj.sympt=traj[,sympt.cum.idx]
    if(loc=="MA"){
      traj.sympt.ma = traj.sympt[,2:ncol(traj.sympt)]
      traj.sympt.ma[,1] = traj.sympt[,1] + traj.sympt[,2]
      traj.sympt = traj.sympt.ma
    }
    ## calculating new cases for each day
    traj.sympt.new = traj.sympt
    traj.sympt.new[-1,] = traj.sympt[-1,] - traj.sympt[-nr,]
    tot.sympt.new=rowSums(traj.sympt.new)
    
    ## cumulative hospitalized by age:
    traj.hosp.cum=traj[,hosp.cum.idx]
    ## calculating new hosp from traj for each day
    traj.hosp.new = traj.hosp.cum
    traj.hosp.new[-1,] = traj.hosp.cum[-1,] - traj.hosp.cum[-nr,]

    if(loc=="MA"){
      traj.hosp.new.ma = traj.hosp.new[,2:ncol(traj.hosp.new)]
      traj.hosp.new.ma[,1] = traj.hosp.new[,1] + traj.hosp.new[,2]
      traj.hosp.new = traj.hosp.new.ma
    }
    
    tot.hosp.new=rowSums(traj.hosp.new)
    ## for proposing new s2.hosp parameter
    hosp.curr=traj[,hosp.curr.idx]
    tot.hosp.curr=rowSums(hosp.curr)
    ## for proposing new s2.icu parameter
    icu.curr=traj[,icu.curr.idx]
    tot.icu.curr=rowSums(icu.curr)
    ## for proposing new s2.vent parameter
    vent.curr=traj[,vent.curr.idx]
    tot.vent.curr=rowSums(vent.curr)
    
    ## deaths
    home.deaths=traj[,home.deaths.idx]
    hosp.deaths=traj[,hosp.deaths.idx]
    deaths.cum=home.deaths+hosp.deaths

    #home.deaths.new=rbind(home.deaths[1,],as.matrix(home.deaths[-1,]-home.deaths[-nrow(traj),]))
    home.deaths.new=home.deaths
    home.deaths.new[-1,] = home.deaths[-1,] - home.deaths[-nr,]

    #hosp.deaths.new=rbind(hosp.deaths[1,],as.matrix(hosp.deaths[-1,]-hosp.deaths[-nrow(traj),]))
    hosp.deaths.new = hosp.deaths
    hosp.deaths.new[-1,] = hosp.deaths[-1,] - hosp.deaths[-nr,]
    
    deaths.new=home.deaths.new+hosp.deaths.new
    
    
    tot.home.deaths.new=rowSums(home.deaths.new)
    tot.hosp.deaths.new=rowSums(hosp.deaths.new)
    tot.deaths.new= tot.home.deaths.new + tot.hosp.deaths.new

    
    if (loc == "MA"){
      # combine first two columns of deaths data
      deaths.new.ma <- deaths.new[,2:ncol(deaths.new)]
      deaths.new.ma[,1] <- deaths.new[,1] + deaths.new[,2]
      deaths.new <- deaths.new.ma
      
      deaths.cum.ma <- deaths.cum[,2:ncol(deaths.cum)]
      deaths.cum.ma[,1] <- deaths.cum[,1] + deaths.cum[,2]
      deaths.cum <- deaths.cum.ma
    }
    
    
    ## hospital discharges
    discharge.idx=281:289 ## for v5
    tot.hosp.discharges=apply(traj[,discharge.idx],1,sum)
    tot.hosp.discharges.new=tot.hosp.discharges
    tot.hosp.discharges.new[-1]=tot.hosp.discharges[-1]-tot.hosp.discharges[-nr]

    ## return
    list(days.odesim=traj[,1],
         tot.sympt.new.odesim=as.numeric(tot.sympt.new),
         sympt.new.odesim=as.matrix(traj.sympt.new),
         tot.hosp.new.odesim=as.numeric(tot.hosp.new),
         hosp.new.odesim=as.matrix(traj.hosp.new),
         deaths.new.odesim=as.matrix(deaths.new),
         hosp.deaths.new.odesim=as.matrix(hosp.deaths.new),
         home.deaths.new.odesim=as.matrix(home.deaths.new),
         tot.hosp.curr.noicu.odesim=tot.hosp.curr,
         tot.icu.curr.novent.odesim=tot.icu.curr,
         tot.vent.curr.odesim=tot.vent.curr,
         tot.hosp.curr.odesim=tot.hosp.curr+tot.icu.curr+tot.vent.curr,
         tot.icu.curr.odesim=tot.icu.curr+tot.vent.curr,
         tot.deaths.new.odesim=as.numeric(tot.deaths.new),
         tot.hosp.deaths.new.odesim=as.numeric(tot.hosp.deaths.new),
         tot.home.deaths.new.odesim=as.numeric(tot.home.deaths.new),
         tot.deaths.cum.odesim=as.numeric(apply(deaths.cum,1,sum,na.rm=TRUE)),
         deaths.cum.odesim=as.matrix(deaths.cum),
         tot.hosp.discharges.new.odesim=as.numeric(tot.hosp.discharges.new)
         )
}
