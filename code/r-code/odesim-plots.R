plots.odesim <- function(dp,tp,rr){
    ## dp =  output from data.process.R
    ## tp = output from traj.process.R
    ## rr = daily vector of reporting rates

    list2env(dp,globalenv())
    list2env(tp,globalenv())
    
    ##
    ## linking times of data and odesim output
    ##

    tmin <- min(days)
    tmax <- max(days)
    idx.min <- which(round(days.odesim) == tmin)
    idx.max <- which(round(days.odesim) == tmax)
    t.idx <- idx.min:idx.max

    ##
    ## linking times of data and odesim output
    ##
    
    
    
    ## linking to report.rate (which is daily starting at day 61)
    t.idx.rr=t.idx-10
    
    ## linking to P and delay rate
    t.idx.delay=tmin:tmax
    
    ##
    ## plot of new symptomatics
    ##

    plot(days,tot.sympt.new,col="blue",pch=20,cex=2,main="Daily Total New Reported Cases")
    points(days.odesim[t.idx],tot.sympt.new.odesim[t.idx]*rr[t.idx.rr],col="red",type="l")


    ##
    ## plot of new hospitalizations
    ##

    if(sum(tot.hosp.new,na.rm=TRUE)>0){
        plot(days,tot.hosp.new,col="blue",pch=20,cex=2,main="Daily Total New Hospitalized Cases")
        points(days.odesim,tot.hosp.new.odesim,col="red",type="l")
    }
    ##
    ## plot of current total hospitalized individuals
    ##

    if(sum(tot.hosp.curr,na.rm=TRUE)>0){
        plot(days,tot.hosp.curr,col="blue",pch=20,cex=2,main="Daily Current Hospitalized Cases")
        points(days.odesim,tot.hosp.curr.odesim,col="red",type="l")
    }
    ##
    ## plot of current total icu individuals
    ##

    if(sum(tot.icu.curr,na.rm=TRUE)>0){
        plot(days,tot.icu.curr,col="blue",pch=20,cex=2,main="Daily Current ICU Cases")
        points(days.odesim,tot.icu.curr.odesim,col="red",type="l")
    }

    ##
    ## plot of current total individuals on ventilators
    ##

    if(sum(tot.vent.curr,na.rm=TRUE)>0){
        plot(days,tot.vent.curr,col="blue",pch=20,cex=2,main="Daily Current Cases on Ventilator")
        points(days.odesim,tot.vent.curr.odesim,col="red",type="l")
    }
    ##
    ## plot of daily new deaths
    ##

    if(sum(tot.deaths.new,na.rm=TRUE)>0){
        plot(days,tot.deaths.new,col="blue",pch=20,cex=2,main="Daily Total New Deaths")
        points(days.odesim,tot.deaths.new.odesim,col="red",type="l")
    }

    ##
    ## plots of age-structured death data
    ##

    if(sum(deaths.cum,na.rm=TRUE)>0){
        matplot(days,deaths.cum,col=1:ncol(deaths.cum),type="p",pch=20,cex=2,main="Daily Cumulative Deaths By Age")
        matplot(days.odesim,deaths.cum.odesim,col=1:ncol(deaths.cum),type="l",add=TRUE,lty=1)
        legend("topleft",legend=(1:9)*10,col=1:9,pch=20)
    }
}
