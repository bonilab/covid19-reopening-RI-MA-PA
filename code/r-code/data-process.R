data.process <- function(df,loc="RI",data.params=NULL, mult = 1){
  # Input
  #   df: data frame object
  #   loc: "State"
  #   data.params: placeholder for parameters shifting the data
  #   mult: c.t multiplier, used to check if random screening code is correct
  # Output
  #   Named list with the following objects:
  #     days, cases.cum, cases.new, cases.age.new,
  #     hosp.cum, hosp.new, hosp.age.new, hosp.curr,
  #     icu.curr, vent.curr, deaths.cum, deaths.new
 
  days <- df$daynum # days after Jan 1, 2020
 
  # calculations for age-structured infill and pre-processing of data
  n.days <- nrow(df)
  # cases
  cases.cum <- df$cumulative_confirmed_cases
  cases.new <- c(cases.cum[1], cases.cum[-1] - cases.cum[-n.days])
  if(loc=="MA"){
    cases.age.cum=df[,16:23]
    hosp.age.cum=df[,24:31]
    deaths.age.cum=df[,32:39]
  }else{
    cases.age.cum=df[,16:24]
    hosp.age.cum=df[,25:33]
    deaths.age.cum=df[,34:42]
  }
  ## observation times of age structured new cases/hosp/deaths
  age.obs.times <- which(!is.na(apply(cases.age.cum,1,sum)))
  age.new.cases <- cases.age.cum[age.obs.times,]
  age.new.cases <- rbind(age.new.cases[1,], 
                         age.new.cases[-1,] - age.new.cases[-length(age.obs.times),])
  
  ### MODIFIED:
  # age.hosp.obs.times=which(!is.na(apply(hosp.age.cum,1,sum)))
  age.hosp.obs.times=which(rowSums(hosp.age.cum, na.rm = T) > 0)
  age.hosp.new.cases=hosp.age.cum[age.hosp.obs.times,]
  age.hosp.new.cases=rbind(age.hosp.new.cases[1,],age.hosp.new.cases[-1,]-age.hosp.new.cases[-length(age.hosp.obs.times),])
  # obs times for age structured deaths
  age.deaths.obs.times=which(!is.na(apply(deaths.age.cum,1,sum)))
  age.deaths.new.cases=deaths.age.cum[age.deaths.obs.times,]
  age.deaths.new.cases=rbind(age.deaths.new.cases[1,],age.deaths.new.cases[-1,]-age.deaths.new.cases[-length(age.deaths.obs.times),])
     
     
  cases.age.new=as.matrix(rbind(cases.age.cum[1,],cases.age.cum[-1,]-cases.age.cum[-n.days,]))
  cases.age.new[cases.age.new<0]=0
  tot.cases.age.new=apply(cases.age.new,1,sum)
  # hospitalizations
  ## browser()
  hosp.age.new=as.matrix(rbind(hosp.age.cum[1,],hosp.age.cum[-1,]-hosp.age.cum[-n.days,]))
  hosp.age.new[hosp.age.new<0]=0
  hosp.cum=df$hospitalized_cumulative
  hosp.new=c(hosp.cum[1],hosp.cum[-1]-hosp.cum[-n.days])
  hosp.curr=df$hospitalized_currently
  
  # icu
  icu.curr=df$InICU_currently
  # if(loc=="MA"){
  #   icu.curr=df$InICU_Currently
  # }
  # # ventilator
  vent.curr=df$OnVentilator_Currently
  # deaths
  tot.deaths.cum=df$cumulative_deaths
  tot.deaths.new=c(tot.deaths.cum[1],tot.deaths.cum[-1]-tot.deaths.cum[-n.days])
  # ## deaths in hospital
  # if(loc=="PA"){
  #   deaths.hosp.cum=NULL
  #   deaths.hosp.new=NULL
  #   deaths.out.of.hosp.new=NULL
  # }else{
  deaths.hosp.cum <- df$Cumulative_hospital_deaths
  deaths.hosp.new=c(deaths.hosp.cum[1],deaths.hosp.cum[-1]-deaths.hosp.cum[-n.days])
  deaths.out.of.hosp.new=tot.deaths.new-deaths.hosp.new
  if(is.null(deaths.hosp.cum[1])){
    deaths.hosp.new <- NULL
    deaths.out.of.hosp.new <- NULL
  }
  # obs times for home/hosp deaths
  hosp.deaths.obs.times <- which(!is.na(deaths.hosp.cum))
  hosp.deaths.new.cases <- deaths.hosp.cum[hosp.deaths.obs.times]
  hosp.deaths.new.cases <- c(hosp.deaths.new.cases[1], 
  hosp.deaths.new.cases[-1] - hosp.deaths.new.cases[-length(hosp.deaths.obs.times)])

   

  ## hospital discharges
  tot.hosp.discharges <- df$Cumulative_hospital_discharges
  tot.new.hosp.discharges <- c(0,tot.hosp.discharges[-1]-tot.hosp.discharges[-n.days])
     
  ## active surveillance (for MA nursinghome tests)
  if(loc=="MA"){
    tests.as=df$nursinghome_tests_total
    pos.as=df$nursinghome_tests_positive
    tests.as=c(tests.as[1], tests.as[-1]-tests.as[-length(tests.as)])
    pos.as=c(pos.as[1], pos.as[-1]-pos.as[-length(pos.as)])
  } else {
    tests.as=rep(NA,nrow(df))
    pos.as=rep(NA,nrow(df))
  }
  
  ## active surveillance for RI
  months.days.2020=list(jan=1:31,
                        feb=32:59,
                        mar=60:90,
                        apr=91:120,
                        may=121:151,
                        jun=152:181,
                        jul=182:212,
                        aug=213:243,
                        sep=244:273,
                        oct=274:304,
                        nov=305:334,
                        dec=335:365)
  idx.days=list()
  for(mmm in 1:12){
    idx.days[[mmm]]=which(is.element(days,months.days.2020[[mmm]]))
  }
  c.t = matrix(0,nrow=length(days),ncol=9)
  ## update below as needed
  ##browser()
  if(loc=="RI"){
    ri.pop.20to70=684347
    ri.pop.70plus=134539
    ## 20 - 70
    c.t[idx.days[[3]],3:7]=0/ri.pop.20to70 ## march
    c.t[idx.days[[4]],3:7]=46.2/ri.pop.20to70 ## april
    c.t[idx.days[[5]],3:7]=484.5/ri.pop.20to70 ## may
    c.t[idx.days[[6]],3:7]=639/ri.pop.20to70 ## june
    c.t[idx.days[[7]],3:7]=888.6/ri.pop.20to70 ## july
    c.t[idx.days[[8]],3:7]=888.6/ri.pop.20to70 ## august (NEED TO CHANGE EVENTUALLY)
    c.t[idx.days[[9]],3:7]=888.6/ri.pop.20to70 ## september (NEED TO CHANGE EVENTUALLY)
    ## 70+ 
    c.t[idx.days[[3]],8:9]=13.8/ri.pop.70plus ## march
    c.t[idx.days[[4]],8:9]=398.8/ri.pop.70plus ## april
    c.t[idx.days[[5]],8:9]=782.7/ri.pop.70plus ## may
    c.t[idx.days[[6]],8:9]=793.4/ri.pop.70plus ## june
    c.t[idx.days[[7]],8:9]=791.4/ri.pop.70plus ## july
    c.t[idx.days[[8]],8:9]=791.4/ri.pop.70plus ## august (NEED TO CHANGE EVENTUALLY)
    c.t[idx.days[[9]],8:9]=791.4/ri.pop.70plus ## september (NEED TO CHANGE EVENTUALLY)
    
    c.t <- c.t * mult
  }
 
  list(days=days,
       tot.sympt.new=cases.new, ## total daily new cases
       sympt.new=cases.age.new, ## age-structured daily new cases
       tot.hosp.new=hosp.new, ## total daily new hospitalizations
       hosp.new=hosp.age.new, ## age-structured daily new hospitalizations
       tot.hosp.curr=hosp.curr, ##note: this is all in hosp (=hosp+icu+vent)
       tot.icu.curr=icu.curr, ## note: this is all in icu (=icu+vent)
       tot.vent.curr=vent.curr, ## current number on ventilators
       tot.deaths.new=tot.deaths.new,
       tot.hosp.deaths.new=deaths.hosp.new,
       tot.home.deaths.new=deaths.out.of.hosp.new,
       deaths.cum=deaths.age.cum,
       tot.deaths.cum=tot.deaths.cum,
       tot.hosp.discharges.new=tot.new.hosp.discharges,
       tests.as=tests.as,
       pos.as=pos.as,
       age.obs.times=age.obs.times,
       age.new.cases=age.new.cases,
       age.hosp.obs.times=age.hosp.obs.times,
       age.hosp.new=age.hosp.new.cases,
       age.deaths.obs.times=age.deaths.obs.times,
       age.deaths.new=age.deaths.new.cases,
       hosp.deaths.obs.times=hosp.deaths.obs.times,
       hosp.deaths.new=hosp.deaths.new.cases,
       c.t=c.t)
 }

