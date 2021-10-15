rm(list=ls())

.libPaths("/global/home/users/cbrook/R/x86_64-pc-linux-gnu-library/4.0")
#setwd("/Users/caraebrook/Documents/R/R_repositories/Berkeley-Reopening/Dec-2020/all-runs/Re-Run-12-24/FigS1/")

#no group, no test, no trace

library(data.table)
library(plyr)
library(dplyr)
library(deSolve)
library(matrixStats)
library(fitdistrplus)

#load parameters including pre-run titer trajectories for each individual
load("titer.dat.20K.Rdata")
#load("titer.dat.2K.Rdata")
load("virus.par.12.15.Rdata")
load("pop.par.base.Rdata")



get.real.cases <- function(pop.dat, event.dat, titer.dat1, within.host.theta, group.limit){
  
  #if no cases caused, then ignore
  if((pop.dat$original_potential_cases_caused_UCB>0) & (pop.dat$num_infection_events>0)){
    
    
    #then allocate all the cases to the events
    #distribute cases at random amongst the events
    event.names <- 1:as.numeric(pop.dat$num_infection_events)
    actual.events <- sample(x=event.names, size=as.numeric(pop.dat$original_potential_cases_caused_UCB), replace = T)
    
    event.data <- cbind.data.frame(actual.events, event.dat[actual.events])
    names(event.data) <- c("event", "gentime")                              
    
    #and add the titer at the time of the event
    gen.tmp = as.list(event.data$gentime)
    event.data$titer <-  c(unlist(lapply(gen.tmp, grab.titer, dat.vir =titer.dat1)))
    
    #now that you have titer, here calculate the probability of transmission, given a certain viral load,
    #based off of the probabiliy model from the URT in Ke et al. 2020
    
    # in Ke et al. 2020, theta is fixed at 0.05 (could be modulated and/or fit to data)
    #draw Km from a normal disribution centered at the midpoint between the two values explored in Ke et al. 2020 (10^3 and 10^4)
    event.data$Km <- rnorm(nrow(event.data),mean=5500, sd=1000)
    
    event.data$prob_exposure = within.host.theta*(event.data$titer/(event.data$titer + event.data$Km))
    event.data$prob_exposure[event.data$prob_exposure<0] <- 0
    #probability is small: ~5% for a typical contact if theta = 0.05 as in Ke.
    #for theta = .7 here, up to 50% depending on theta
    
    #does the infection happen? make it a probabilistic outcome of the titer
    #then, you role a dice to see if this exposure causes an infection
    tmp.prob <- as.list(event.data$prob_exposure)
    
    event.data$InfectionYN = c(unlist(lapply(tmp.prob, test.titer)))
    
    
    
    #then total the events that actually happen to incorporate into the original data
    pop.dat$post_titer_potential_cases_caused_UCB <- sum(event.data$InfectionYN)
    
    #and then, if there is a group size limit, impose it here
    if((group.limit>0) & (pop.dat$obs_dist_limits==TRUE)){
      #gives you the number of successful transmissions per event
      event.sum <- ddply(event.data, .(event),summarize, N=sum(InfectionYN)) 
      event.sum$over_lim = event.sum$N-group.limit
      event.sum$over_lim[event.sum$over_lim<0] <- 0
      
      
      #truncate # of events for the IDs listed above to the group limit.
      event.data.list = dlply(subset(event.data, InfectionYN==1), .(event))
      
      
      
      new.event.list <- lapply(event.data.list, impose.group, group.limit=group.limit)
      
      #new.event.data <- do.call("rbind", new.event.list)
      new.event.data <-data.table::rbindlist(new.event.list)
      
      pop.dat$potential_cases_caused = sum(new.event.data$InfectionYN)
      
      #in this case, return the generation time table after the group intervention
      if(pop.dat$potential_cases_caused >0){
        dat.gen.tab <- cbind.data.frame(rep(unique(pop.dat$employ_ids), nrow(new.event.data)), new.event.data$gentime)
        names(dat.gen.tab) <- c("employ_ids", "generation_time")
      }else{
        dat.gen.tab <- cbind.data.frame(unique(pop.dat$employ_ids), NA)
        names(dat.gen.tab) <- c("employ_ids", "generation_time")
      }
      
      
    }else{
      
      pop.dat$potential_cases_caused <- pop.dat$post_titer_potential_cases_caused_UCB
      if(pop.dat$potential_cases_caused >0){
        event.data.out = subset(event.data, InfectionYN==1)
        dat.gen.tab <- cbind.data.frame(rep(unique(pop.dat$employ_ids), nrow(event.data.out)),   event.data.out$gentime)
        names(dat.gen.tab) <- c("employ_ids", "generation_time")
      }else{
        dat.gen.tab <- cbind.data.frame(unique(pop.dat$employ_ids), NA)
        names(dat.gen.tab) <- c("employ_ids", "generation_time")
      }
    }
    
  }else{ 
    #none take place
    
    #return the original data with 0s
    pop.dat$post_titer_potential_cases_caused_UCB <- 0
    pop.dat$potential_cases_caused <- 0
    
    #and return a table of generation times with nothing
    dat.gen.tab <- cbind.data.frame(unique(pop.dat$employ_ids), NA)
    names(dat.gen.tab) <- c("employ_ids", "generation_time")
    
  }
  
  
  return(list(pop.dat, dat.gen.tab))
  
}
test.titer <- function(prob1){
  Y_N =sample(c(0,1), size=1, prob  = c(1-prob1, prob1))
  return(Y_N)
}
impose.group <- function(event.dat1, group.limit){
  tot.transmissions = nrow(event.dat1)
  if(tot.transmissions>group.limit){
    choose.events <- sample(x=1:tot.transmissions, size=group.limit, replace = F)
    event.dat2 = event.dat1[choose.events,]
    return(event.dat2)
  }else{
    return(event.dat1)
  }
}
get.event.time <- function(dat, genTime){
  event.times = genTime(as.numeric(dat$num_infection_events))
  return(event.times)
}
grab.titer <- function(dat1, dat.vir){
  titer.out <- dat.vir$V[dat.vir$time>dat1][1]
  return(titer.out)
}
normal_fn <- function(meanpar=NULL, sdpar=NULL){
  out <- purrr::partial(rnorm,
                        mean = meanpar,
                        sd = sdpar)
  return(out)
} 
poisson_fn <- function(lambda=NULL){
  out <- purrr::partial(rpois,
                        lambda = lambda)
  return(out)
} 
lognormal_fn <- function(meanlogpar=NULL, sdlogpar=NULL){
  out <- purrr::partial(rlnorm,
                        meanlog = meanlogpar,
                        sdlog = sdlogpar)
  return(out)
} 
add.risk.cat <- function(dat, pop_dat){
  dat = data.table(dat)
  
  daily_new <- dat[, day := ceiling(time_isolation) #time_isolated
                   ][, .(daily_isolations = .N), by = day
                     ]
  
  pop_dat$add <- 0
  for (i in 1:length(daily_new$day)){
    pop_dat$add[pop_dat$day==daily_new$day[i]] <- daily_new$daily_isolations[i]
  }
  
  
  out.vect  <- as.data.frame(pop_dat$add)
  
  return(  out.vect)
}
add.risk.cat.exp <- function(dat, pop_dat, input_par){
  dat = data.table(dat)
  
  daily_new <- dat[, day := ceiling(exposure_time) 
                   ][, .(daily_exposures = .N), by = day
                     ]
  
  pop_dat$add <- 0
  for (i in 1:length(daily_new$day)){
    pop_dat$add[pop_dat$day==daily_new$day[i]] <- daily_new$daily_exposures[i]
  }
  
  
  out.vect  <- as.data.frame(pop_dat$add)
  
  
  #then add deaths based on each pop cat
  pop.cat = unique(dat$employ_cat)
  
  out.vect2 <- as.data.frame(as.numeric(input_par$par1[input_par$parameter=="CFR" & input_par$population==pop.cat])*out.vect)
  
  # dat.out = cbind.data.frame(out.vect, out.vect2)
  return(list(out.vect, out.vect2))
  #return(dat.out)
}
cross.infect <- function(dat, all.sus, input.par){
  pop.par = subset(input.par, population == unique(dat$infector_cat))
  
  #first, elim any populations for which there are no longer remaining susceptibles
  
  rem.cat = unique(all.sus$employ_cat)
  all.cat = unique(pop.par$par2[pop.par$parameter=="meta-pop"])
  
  missed.cat = setdiff(all.cat, rem.cat)
  
  pop.par$sub = 0
  for (i in 1: length(missed.cat)) {
    pop.par$sub[pop.par$parameter=="meta-pop" & pop.par$par2==missed.cat[i]] <- 1
  }
  
  pop.par = subset(pop.par, sub==0)
  
  #then allocate the population of the new cases based on the proportion within and without 
  tot.cases = nrow(dat)
  
  #then need to reallocate probabilities comparatively without the remaining
  possible.cat = unique(pop.par$par2[pop.par$parameter=="meta-pop"])
  
  old.cat = as.numeric(unique(input.par$par2[input.par$parameter=="meta-pop"]))
  old.prob = as.numeric(input.par$par1[input.par$parameter=="meta-pop"])[1:length(old.cat)]
  
  if(length(possible.cat)<length(old.cat)){
    if(length(possible.cat)==1){
      dat$new_cat = possible.cat
    }else{
      #if you've run out of probabilities, just, rellocate proportionally
      new.prob = rep((1/length(possible.cat)), length(possible.cat))
      dat$new_cat = sample(x=possible.cat, size = tot.cases, replace = TRUE, prob = new.prob)
      
    }
  }else{
    dat$new_cat = sample(x=old.cat, size = tot.cases, replace = TRUE, prob = old.prob)
    
  }
  
  
  
  
  return(dat)
  
  
}
assign.ID = function(sus.dat.sub, dat.new.sub){
  #at the very end of the time series, you may run out of susceptibles in the right category, in which case, these just become lost infections
  if(nrow(dat.new.sub)<=length(sus.dat.sub$employ_ids)){
    dat.new.sub$new_infected =    sample(sus.dat.sub$employ_ids, size=nrow(dat.new.sub), replace=FALSE)  
  }else{
    new.count = length(sus.dat.sub$employ_ids)
    new.missed = nrow(dat.new.sub) - new.count
    row.tmp = seq(1, nrow(dat.new.sub),1)
    row.take = sample(row.tmp, size = new.count, replace = FALSE)
    dat.new.sub <- dat.new.sub[row.take,] 
    dat.new.sub$new_infected =    sample(sus.dat.sub$employ_ids, size=nrow(dat.new.sub), replace=FALSE)  
  }
  
  return(dat.new.sub)
}
assign.infections <- function(pop.mat, gen_list, timestep, input.par){
  # assign new exposures (and times) based on 'actual cases caused' above
  # and move those that have transmitted to isolated/recovered state
  #(asymptomatics will be missed in iso time  unless tested)
  # timestep.prev = unique(pop.mat$timestep)
  
  
  #first, pair each case with its generation times
  new.mat <- dplyr::select(pop.mat, employ_ids, employ_cat, state, exposure_time, actual_cases_caused, time_isolation)
  new.mat <- new.mat[!is.na(new.mat$actual_cases_caused) & new.mat$state==1,]
  
  #only matters if it actually causes  cases.
  new.mat.zero = subset(new.mat, actual_cases_caused<1)
  new.mat <- subset(new.mat, actual_cases_caused>0)
  
  if(nrow(new.mat)>0){
    
    new.mat.list <- dlply(new.mat, .(employ_ids))
    #print("1")
    new.mat.list <-  lapply(new.mat.list, make.rows)
    
    #the new new mat - no longer includes those which cased 0 actual cases
    #should always have at least one row because of the if-statement above
    #new.mat <- do.call("rbind", new.mat.list)
    new.mat <- data.table::rbindlist(new.mat.list)
    
    #now attach a generation time with each of these cases and a random sample from the susceptibles
    new.mat$generation_time <- NA
    index.ids = unique(new.mat$employ_ids)
    
    
    for(i in 1:length(index.ids )){
      tmp = nrow(new.mat[new.mat$employ_ids == index.ids[i],])
      #print(index.ids[[i]])
      new.mat$generation_time[new.mat$employ_ids == index.ids[i]] <- gen_list[[index.ids[i]]]$generation_time[1:tmp]
    }
    
    
    
    #now, attach a place to infect (susceptible)
    #bias the sampling based on the proportion of infections within and without of your direct cohort 
    
    #first, pair the remaining susceptibles with their category
    
    
    all.sus <- cbind.data.frame(pop.mat$employ_ids[pop.mat$state==0],pop.mat$employ_cat[pop.mat$state==0])
    names(all.sus) = c("employ_ids", "employ_cat")
    new.list = dlply(new.mat, .(employ_ids))
    
    #cross infect by cat
    #print("2")
    new.list.out <- lapply(new.list, cross.infect, all.sus=all.sus, input.par=input.par)
    new.mat = data.table::rbindlist(new.list.out)
    #new.mat = do.call("rbind", new.list.out)
    rownames(new.mat) <- c()
    
    #then, assign names by category of new infections
    id.cat = data.frame(sort(unique(new.mat$new_cat)))
    all.sus = arrange(all.sus, employ_cat)
    
    names(id.cat) <- "employ_cat"
    tmp.sus = merge(x=all.sus, y=id.cat)
    tmp.sus.split = dlply(tmp.sus, .(employ_cat))
    new.mat.split <- dlply(new.mat, .(new_cat))
    
    #print("3")
    dat.new.split.out = mapply(FUN=assign.ID, sus.dat.sub= tmp.sus.split, dat.new.sub= new.mat.split, SIMPLIFY = FALSE)
    new.mat = data.table::rbindlist(dat.new.split.out)
    #new.mat = do.call("rbind", dat.new.split.out)
    
    new.mat$new_exposure_time = new.mat$exposure_time + new.mat$generation_time
    
    
    #and merge into pop.mat
    #new.merge <- dplyr::select(new.mat, new_infected, employ_ids, infector_iso_time, new_exposure_time)
    #names(new.merge) <- c("employ_ids", "infector", "infector_iso_time", "exposure_time")
    
    
    
    #now put them into pop.mat
    for(i in 1:nrow(new.mat)){
      
      #identify infector and iso time
      pop.mat$infector[pop.mat$employ_ids==new.mat$new_infected[i] ] <- new.mat$employ_ids[i]
      pop.mat$infector_iso_time[pop.mat$employ_ids==new.mat$new_infected[i]] <- new.mat$infector_iso_time[i] 
      
      pop.mat$time_of_tracing_iso[pop.mat$employ_ids==new.mat$new_infected[i] & pop.mat$traced==TRUE] <-  pop.mat$time_of_tracing_iso[pop.mat$employ_ids==new.mat$new_infected[i] & pop.mat$traced==TRUE] + new.mat$infector_iso_time[i]
      pop.mat$time_of_tracing_iso[pop.mat$employ_ids==new.mat$new_infected[i] & pop.mat$traced==FALSE] <- Inf #if traced==FALSE, this is never tracked
      
      #and exposure time
      pop.mat$exposure_time[pop.mat$employ_ids==new.mat$new_infected[i]] <- new.mat$new_exposure_time[i]
      #pop.mat$time_test_sensitive_start[pop.mat$employ_ids==new.mat$new_infected[i]] <-new.mat$new_exposure_time[i] + pop.mat$time_test_sensitive_start[pop.mat$employ_ids==new.mat$new_infected[i]] 
      #pop.mat$time_test_sensitive_end[pop.mat$employ_ids==new.mat$new_infected[i]] <-new.mat$new_exposure_time[i] + pop.mat$time_test_sensitive_end[pop.mat$employ_ids==new.mat$new_infected[i]] 
      pop.mat$case_source[pop.mat$employ_ids==new.mat$new_infected[i]] <- "UCB" #transmission within berkeley
      
      #change state - only if exposure time is already achieved
      pop.mat$state[pop.mat$employ_ids==new.mat$new_infected[i] & pop.mat$exposure_time<=timestep] <- 1
      
      
      #otherwise, they still stay suceptible - but you mark them 
      pop.mat$state[pop.mat$employ_ids==new.mat$new_infected[i] & pop.mat$exposure_time>timestep] <- 3
      
      
    } #else, just return pop mat
    
    
  }
  #now, make those that already transmitted recovered/isolated
  pop.mat$state[(pop.mat$state==1 & !is.na(pop.mat$actual_cases_caused))] <- 5
  #and, if any of the old "pre-exposed" have reached their exposure time, you can 
  
  #and return pop.mat
  return(pop.mat)
}
assign.last.infections <- function(pop.mat, gen_list, remaining.susceptibles, timestep){
  # assign new exposures (and times) based on 'actual cases caused' above
  # and move those that have transmitted to isolated/recovered state
  #(asymptomatics will be missed in iso time  unless tested)
  timestep.prev = unique(pop.mat$timestep)
  
  
  if(remaining.susceptibles>0){
    
    
    #first, pair each case with its generation times
    new.mat <- dplyr::select(pop.mat, employ_ids, state, exposure_time, actual_cases_caused, time_isolation)#, time_of_tracing_iso)
    new.mat <- new.mat[ new.mat$state==1 & !is.na(new.mat$actual_cases_caused) ,]
    
    #get rid of those that cause no cases
    new.mat <- new.mat[new.mat$actual_cases_caused>0,]
    #sum(new.mat$actual_cases_caused)>remaining susceptibles
    
    #so need to pick these at random to generate new infections instead
    all.possible = c(rep(new.mat$employ_ids, times=new.mat$actual_cases_caused))
    last.infector.ids = sample(all.possible, size=remaining.susceptibles,  replace=FALSE)
    
    
    last.infector.ids = data.frame(last.infector.ids)
    names( last.infector.ids) ="employ_ids"
    
    new.dat  = ddply(last.infector.ids,.(employ_ids), summarise, actual_cases_caused=length(employ_ids))
    
    #and new.mat becomes just these
    new.dat$state <- new.dat$time_isolation <- new.dat$exposure_time <- NA
    
    for (i in 1:nrow(new.mat)){
      new.dat$state[new.dat$employ_ids==new.mat$employ_ids[i]] <- new.mat$state[i]
      new.dat$time_isolation[new.dat$employ_ids==new.mat$employ_ids[i]] <- new.mat$time_isolation[i]
      new.dat$exposure_time[new.dat$employ_ids==new.mat$employ_ids[i]] <- new.mat$exposure_time[i]
    }
    
    #then, new dat takes the place of new mat
    
    new.dat.list <- dlply(new.dat, .(employ_ids))
    
    new.dat.list <-  lapply(new.dat.list, make.rows)
    new.dat <- data.table::rbindlist(new.dat.list)
    #new.dat <- do.call("rbind", new.dat.list)
    
    #now attach a generation time with each of these cases and a random sample from the susceptibles
    new.dat$generation_time <- NA
    index.ids = unique(new.dat$employ_ids)
    
    for(i in 1:length(index.ids )){
      #print(index.ids[[i]])
      new.mat$generation_time[new.mat$employ_ids == index.ids[i]] <- gen_list[[index.ids[i]]]$generation_time[1:length(new.mat$generation_time[new.mat$employ_ids == index.ids[i]])]
    }
    
    
    
    
    #now, attach a place to infect (susceptible) -- should be enough
    all.sus <- pop.mat$employ_ids[pop.mat$state==0]
    new.dat$new_infected <- sample(all.sus, size=nrow(new.dat), replace=FALSE)
    
    new.dat$new_exposure_time = new.dat$exposure_time + new.dat$generation_time 
    
    #now put them into pop.mat
    for(i in 1:nrow(new.dat)){
      
      #identify infector and iso time
      pop.mat$infector[pop.mat$employ_ids==new.dat$new_infected[i] ] <- new.dat$employ_ids[i]
      pop.mat$infector_iso_time[pop.mat$employ_ids==new.dat$new_infected[i]] <- new.dat$infector_iso_time[i] 
      
      pop.mat$time_of_tracing_iso[pop.mat$employ_ids==new.dat$new_infected[i] & pop.mat$traced==TRUE] <-  pop.mat$time_of_tracing_iso[pop.mat$employ_ids==new.dat$new_infected[i] & pop.mat$traced==TRUE] + new.dat$infector_iso_time[i]
      pop.mat$time_of_tracing_iso[pop.mat$employ_ids==new.dat$new_infected[i] & pop.mat$traced==FALSE] <- Inf #if traced==FALSE, this is never tracked
      
      #and exposure time
      pop.mat$exposure_time[pop.mat$employ_ids==new.dat$new_infected[i]] <- new.dat$new_exposure_time[i]
      #pop.mat$time_test_sensitive_start[pop.mat$employ_ids==new.mat$new_infected[i]] <-new.mat$new_exposure_time[i] + pop.mat$time_test_sensitive_start[pop.mat$employ_ids==new.mat$new_infected[i]] 
      #pop.mat$time_test_sensitive_end[pop.mat$employ_ids==new.mat$new_infected[i]] <-new.mat$new_exposure_time[i] + pop.mat$time_test_sensitive_end[pop.mat$employ_ids==new.mat$new_infected[i]] 
      
      pop.mat$case_source[pop.mat$employ_ids==new.dat$new_infected[i]] <- "UCB" #transmission within berkeley
      
      #change state - only if exposure time is already achieved
      pop.mat$state[pop.mat$employ_ids==new.dat$new_infected[i] & pop.mat$exposure_time<=timestep] <- 1
      
      
      #otherwise, they still stay suceptible - but you mark them 
      pop.mat$state[pop.mat$employ_ids==new.dat$new_infected[i] & pop.mat$exposure_time>timestep] <- 3
      
      
    } #else, just return pop mat
  }
  #other
  #now, make those that already transmitted recovered/isolated
  pop.mat$state[(pop.mat$state==1 & !is.na(pop.mat$actual_cases_caused))] <- 5
  #and, if any of the old "pre-exposed" have reached their exposure time, you can 
  
  #then, those active infections will cause no more new cases
  pop.mat$actual_cases_caused[pop.mat$state==1] <- 0
  
  #and return pop.mat
  return(pop.mat)
}
get.actual.cases = function(pop.dat, dat.gen, timestep){
  sub.gen =subset(dat.gen, employ_ids==unique(pop.dat$employ_ids))
  #count the number of exposure time + generation time combos that take place before the iso time
  sub.gen$new_exposures = sub.gen$generation_time + pop.dat$exposure_time
  n.out = length(sub.gen$new_exposures[sub.gen$new_exposures<=pop.dat$time_isolation])
  return(n.out)
}
get.symptom.onset <- function(dat, dat.vir, LOD){
  #get titer limit
  symptom.lim <- as.numeric(unique(dat$titer_lim_for_symptoms))
  #get the timing in the trajectory that first crosses this limit
  dat$time_of_symptom_onset  <- min(dat.vir$time[dat.vir$V>symptom.lim])
  
  dat$time_test_sensitive_start <- min(dat.vir$time[dat.vir$V>LOD])
  dat$time_test_sensitive_end <- max(dat.vir$time[dat.vir$V>LOD])
  
  #will return infinity if wrong
  
  return(dat)  
}
make.rows <- function(dat){
  n = dat$actual_cases_caused
  new.dat <- data.frame(matrix(NA, nrow=n, ncol=5) )
  names(new.dat) <- c("employ_ids", "exposure_time", "actual_cases_caused", "infector_iso_time", "infector_cat")#, "time_of_test_sensitivity")#, "time_of_tracing_iso")
  new.dat$employ_ids <- rep(dat$employ_ids, nrow(new.dat))
  new.dat$infector_iso_time <- rep(dat$time_isolation, nrow(new.dat))
  new.dat$infector_cat <- rep(dat$employ_cat, nrow(new.dat))
  new.dat$exposure_time <- rep(dat$exposure_time, nrow(new.dat))
  #new.dat$time_of_tracing_iso <- rep(dat$time_of_tracing_iso, nrow(new.dat))
  if(nrow(new.dat)>0){
    new.dat$actual_cases_caused <- 1
    return(new.dat)  
  } #else, return nothing
}
delayfn_surv <- function(delay_mean, delay_sd){
  
  out <- purrr::partial(rnorm,
                        mean = delay_mean,
                        sd = delay_sd)
  
  return(out)
}#symptomatic surveillance/TAT delay
generationTime_fn <- function(serial_dist=NULL, serial_shape = NULL, serial_scale = NULL) {
  if(serial_dist=="weibull"){
    out <- purrr::partial(rweibull,
                          shape = serial_shape,
                          scale = serial_scale)  
  }else if(serial_dist=="gamma"){
    out <- purrr::partial(rgamma,
                          shape = serial_shape,
                          scale = serial_scale) 
  }
  
  return(out)
} #weibull or gamma serial interval as the case may be
inc_fn <- function(n_inc_samp = NULL, meanInc=NULL, sdInc=NULL) {
  
  out= purrr::partial(rlnorm,
                      meanlog = log(meanInc), 
                      sdlog = log(sdInc))
  
  #out[out < 1] <-  1
  
  return(out)
} #lognormal incubation time draw
R0_fn <- function(meanR0=NULL, sdR0=NULL){
  out <- purrr::partial(rlnorm,
                        meanlog = log(meanR0),
                        sdlog = log(sdR0))
  return(out)
} #lognormal R0
R0_fn_nb <- function(muR0=NULL, sizeR0=NULL){
  out <- purrr::partial(rnbinom,
                        mu = muR0,
                        size = sizeR0)
  return(out)
} #nb R0
initiate.pop <- function(start.ID.employ, pop.UCB, n.init.exposed, pop.ID, within.host.theta, input.par,  R0fn, eventFn, titer.dat, LOD, virus.par){
  
  
  #sample serial interval
  genTime = generationTime_fn(serial_dist = virus.par$distribution[virus.par$parameter=="generation_time"],
                              serial_shape= virus.par$par1[virus.par$parameter=="generation_time"], 
                              serial_scale= virus.par$par2[virus.par$parameter=="generation_time"])
  
  
  pop.par = subset(input.par, population==pop.ID)
  
  #make table one
  pop.mat = cbind.data.frame(matrix(NA, nrow=pop.UCB, ncol =27))
  names(pop.mat) <- c( "employ_ids","employ_cat", "state",  "traced",  "testing", "obs_dist_limits",  "exposure_time",   "total_potential_cases_caused",  "original_potential_cases_caused_UCB", "num_infection_events", "post_titer_potential_cases_caused_UCB", "potential_cases_caused", "actual_cases_caused", "case_source", "infector", "time_test_sensitive_start", "time_test_sensitive_end", "infector_iso_time", "time_of_tracing_iso", "time_of_next_test", "time_of_testing_iso", "titer_lim_for_symptoms", "time_of_symptom_onset", "time_of_symptom_iso", "time_isolation", "reason_isolated",  "timestep")
  
  #and fill in all you can
  pop.mat$testing = pop.par$par1[pop.par$parameter=="test-on"]
  pop.mat$timestep = 0
  pop.mat$employ_cat = pop.ID
  
  #assign them all an employer ID 
  pop.mat$employ_ids = start.ID.employ:(pop.UCB+start.ID.employ-1)
  
  #assign a "first test" depending on how many days of testing per week...
  test_rotation = as.character(pop.par$par1[pop.par$parameter=="test-rotation"] )
  n.test.day.per.wk = as.numeric(pop.par$par1[pop.par$parameter=="n-test-days-per-week"])
  
  #then, if this is bigger than a weekly regime, half of them must happen one week and half the other
  if ((test_rotation=="biweekly" & n.test.day.per.wk==2) | (test_rotation=="weekly" & n.test.day.per.wk==2)){
    pop.mat$time_of_next_test = rep(c(3,7), length=pop.UCB) 
    
  }else if((test_rotation=="biweekly" & n.test.day.per.wk==5) | (test_rotation=="weekly" & n.test.day.per.wk==5)){
    
    pop.mat$time_of_next_test = rep(c(3,4,5,6,7), length=pop.UCB)
    
  }else if((test_rotation=="biweekly" & n.test.day.per.wk==7) | (test_rotation=="weekly" & n.test.day.per.wk==7)){
    pop.mat$time_of_next_test = rep(c(1,2,3,4,5,6,7), length=pop.UCB)
    
  }else if(test_rotation=="two-week" & n.test.day.per.wk==7){
    
    pop.mat$time_of_next_test = rep(c(1,2,3,4,5,6,7,8,9,10,11,12,13,14), length=pop.UCB)
    
  }else if(test_rotation=="two-week" & n.test.day.per.wk==5){
    
    pop.mat$time_of_next_test = rep(c(3,4,5,6,7,10,11,12,13,14), length=pop.UCB)
    
  }else if(test_rotation=="two-week" & n.test.day.per.wk==2){
    
    pop.mat$time_of_next_test = rep(c(3,7, 10,14), length=pop.UCB)
    
  }else if (test_rotation=="two-week-ThFri"){
    
    pop.mat$time_of_next_test = rep(c(seq((7-n.test.day.per.wk+1),7,1),seq((14-n.test.day.per.wk+1),14,1)), length=pop.UCB) #end of week. if 1 and 2 are weekend (sat/sun), then this is thursday/friday
    
  }else if (test_rotation=="two-week-daily"){
    
    pop.mat$time_of_next_test = rep(c(seq(1,14,1)), length=pop.UCB) #end of week. if 1 and 2 are weekend (sat/sun), then this is thursday/friday
    
  }else if (test_rotation=="two-week-MonTues"){
    
    pop.mat$time_of_next_test = rep(c(seq(3,4,1),seq(10,11,1)), length=pop.UCB) #end of week. if 1 and 2 are weekend (sat/sun), then this is thursday/friday
    
  }else if (test_rotation=="two-week-TuesWed"){
    
    pop.mat$time_of_next_test = rep(c(seq(4,5,1),seq(11,12,1)), length=pop.UCB) #end of week. if 1 and 2 are weekend (sat/sun), then this is thursday/friday
    
  }else if (test_rotation=="two-week-MonFri"){
    
    pop.mat$time_of_next_test = rep(c(3,7,10,14), length=pop.UCB) #end of week. if 1 and 2 are weekend (sat/sun), then this is thursday/friday
    
  }else if (test_rotation=="two-week-MonWed"){
    
    pop.mat$time_of_next_test = rep(c(3,5,10,12), length=pop.UCB) #end of week. if 1 and 2 are weekend (sat/sun), then this is thursday/friday
    
  }else if (test_rotation=="one-week-ThFri"){
    pop.mat$time_of_next_test = rep(seq((7-n.test.day.per.wk+1),7,1), length=pop.UCB)
  }else if (test_rotation=="one-week-MonTues"){
    pop.mat$time_of_next_test = rep(seq(3,4,1), length=pop.UCB)
  }else if (test_rotation=="one-week-MonFri"){
    pop.mat$time_of_next_test = rep(c(3,7), length=pop.UCB)
  }else if(test_rotation=="one-week-daily"){
    pop.mat$time_of_next_test = rep(c(seq(1,7,1)), length=pop.UCB)
  }else if(test_rotation=="none"){
    pop.mat$time_of_next_test=Inf
  }else if(test_rotation=="thrice-weekly-MonTues"){
    pop.mat$time_of_next_test = rep(c(3,4,10,11,17,18), length=pop.UCB)
  }else if (test_rotation=="two_day"){
    pop.mat$time_of_next_test = rep(c(3,7), length=pop.UCB)
  }else if (test_rotation=="four_week"){
    pop.mat$time_of_next_test = rep(c(seq((7-n.test.day.per.wk+1),7,1), seq((14-n.test.day.per.wk+1),14,1), seq((21-n.test.day.per.wk+1),21,1), seq((28-n.test.day.per.wk+1),28,1)), length=pop.UCB)
  }
  
  pop.mat$time_of_next_test = sample(pop.mat$time_of_next_test, size=length(pop.mat$time_of_next_test), replace = FALSE) #scramble
  
  
  
  prop.traced = as.numeric(pop.par$par1[pop.par$parameter=="prop.trace"])
  
  #for all, based on proportions, give whether traced or not
  pop.mat$traced = sample(x=c("TRUE", "FALSE"), size=length(pop.mat$traced), replace=TRUE, prob=c(prop.traced, 1-prop.traced))
  
  #and the same for proportion observign distancing limits
  prop.obs = as.numeric(pop.par$par1[pop.par$parameter=="percent-obs-dist-lim"])
  
  pop.mat$obs_dist_limits = sample(x=c("TRUE", "FALSE"), size=length(pop.mat$obs_dist_limits), replace=TRUE, prob=c(prop.obs, 1-prop.obs))
  
  # and whether asymp or not
  #pop.mat$stat_asymp = sample(x=c("TRUE", "FALSE"), size=length(pop.mat$stat_asymp), replace=TRUE, prob=c(prop.asym, 1-prop.asym))
  
  
  #make initial state variable
  pop.mat$state <- rep(as.integer(0),pop.UCB)
  
  #based on the proportion vaccinated, some get moved to recovered (state 5) right away 
  #for all of our model runs, this is 0, so this gets skipped
  if(as.numeric(pop.par$par1[pop.par$parameter=="prop-vaccinated"])>0){
    tot.vacc <- round(as.numeric(pop.par$par1[pop.par$parameter=="prop-vaccinated"])*pop.UCB,0)
    index.vacc = sample(1:pop.UCB, size=tot.vacc, replace=FALSE)
    pop.mat$state[index.vacc] <- 5
  }
  
  #then, regardless of vaccination, overwrite susceptibles with those initially exposed
  
  #initially exposed get distributed at random
  index.init = sample(as.numeric(rownames(pop.mat[pop.mat$state==0,])), size=n.init.exposed, replace=FALSE)
  pop.mat$state[index.init] <- 1
  
  #here, build distributions
  #symptomatic isolation delay
  delayfn_symp = delayfn_surv(delay_mean=as.numeric(pop.par$par1[pop.par$parameter=="iso-lag"]), 
                              delay_sd= as.numeric(pop.par$par2[pop.par$parameter=="iso-lag"]))
  
  #turnaround testing delay
  delayfn_TAT = delayfn_surv(delay_mean=as.numeric(pop.par$par1[pop.par$parameter=="TAT-lag"]), 
                             delay_sd= as.numeric(pop.par$par2[pop.par$parameter=="TAT-lag"]))
  
  
  #contact tracing lag
  delayfn_trace = delayfn_surv(delay_mean=as.numeric(pop.par$par1[pop.par$parameter=="trace-lag"]), 
                               delay_sd= as.numeric(pop.par$par2[pop.par$parameter=="trace-lag"]))
  
  #titer  limit for symptoms
  titer_lim = lognormal_fn(meanlogpar=log(as.numeric(pop.par$par1[pop.par$parameter=="symptom-lim"])), 
                           sdlogpar = log(as.numeric(pop.par$par2[pop.par$parameter=="symptom-lim"])))
  
  prop.cases.UCB = as.numeric(pop.par$par1[pop.par$parameter=="prop.cases.UCB"])
  
  #now generate potential new infections based on your status
  #this gives the weekend average number of infections
  #pop.mat$tot_potential_cases_caused[pop.mat$stat_asymp==TRUE] = floor(R0fn.asym.wk(length(pop.mat$tot_potential_cases_caused[pop.mat$stat_asymp==TRUE]))*prop.cases.UCB) 
  #pop.mat$tot_potential_cases_caused[pop.mat$stat_asymp==FALSE] = floor(R0fn.wk(length(pop.mat$tot_potential_cases_caused[pop.mat$stat_asymp==FALSE]))*prop.cases.UCB)
  
  
  #and during the week, fewer cases
  #pop.mat$wk_tot_potential_cases_caused[pop.mat$stat_asymp==TRUE] =  floor(R0fn.asym(length(pop.mat$tot_potential_cases_caused[pop.mat$stat_asymp==TRUE]))*prop.cases.UCB) 
  
  #here are all possible cases
  pop.mat$total_potential_cases_caused =  R0fn(length(pop.mat$employ_ids))
  
  #here are all possible at UC Berkeley - before the titer cull
  pop.mat$original_potential_cases_caused_UCB = floor(pop.mat$total_potential_cases_caused*prop.cases.UCB)
  
  
  #you should have already brought in a titer trajectory for everyone in your population 
  #choose a threshold titer for symptom onset 
  pop.mat$titer_lim_for_symptoms = titer_lim(pop.UCB)
  
  pop.split <- dlply(pop.mat, .(employ_ids))
  titer.split <- dlply(titer.dat, .(employ_ids))
  
  #now, based on this, go into each person's virus trajectory and calculate the timing of symptom onset
  #while you are at it, you can also look at their titer and the LOD and calculate the start/end times for which they are test sensitive
  
  pop.split.new <- mapply(get.symptom.onset, dat = pop.split, dat.vir=titer.split, MoreArgs = list(LOD=LOD), SIMPLIFY = FALSE)
  #when there is nothing that is under the limit, these infections become "asymptomatic" -- 
  #we can later play with the proportion that classify as this by modulating the mean value for the symptom onset limit
  
  pop.mat <- data.table::rbindlist(pop.split.new)
  #pop.mat <- do.call("rbind", pop.split.new)
  
  #and the delay to isolation
  pop.mat$time_of_symptom_onset[pop.mat$time_of_symptom_onset<0]<-0
  pop.mat$time_of_symptom_iso = delayfn_symp(pop.UCB)
  pop.mat$time_of_symptom_iso[pop.mat$time_of_symptom_iso<0]<- 0
  pop.mat$time_of_symptom_iso <- pop.mat$time_of_symptom_iso + pop.mat$time_of_symptom_onset
  pop.mat$time_of_testing_iso =  delayfn_TAT(pop.UCB)
  
  pop.mat$time_of_testing_iso[pop.mat$time_of_testing_iso<0] <- 0
  pop.mat$time_of_testing_iso <- pop.mat$time_of_testing_iso + pop.mat$time_of_next_test
  
  pop.mat$time_of_tracing_iso = delayfn_trace(pop.UCB)
  pop.mat$time_of_tracing_iso[pop.mat$time_of_tracing_iso<0] <- 0
  
  
  #now, if not traced, never:
  pop.mat$time_of_tracing_iso[pop.mat$traced==FALSE] <- Inf
  pop.mat$time_of_tracing_iso[pop.mat$state>0] <- Inf # new introductions cannot be traced
  pop.mat$infector[pop.mat$state>0] <- 0 # new introductions cannot be traced
  pop.mat$infector_iso_time[pop.mat$state>0] <- Inf # new introductions cannot be traced
  pop.mat$case_source[pop.mat$state>0] <- "alameda"
  
  #NOW, we generate new cases:
  #we break down each infectious individual based on that individual's: 
  #(a) within-host titer trajectory, (b) the selected value for within-host theta (how viral load translates toinfection probability), 
  #(c) the number of discrete transmission events that we draw for each person, and 
  #(d) the generation time of those contact events 
  #(for d, we currently use the Ferretti weibull, but we are hoping that a constant hazard of events
  # + the titer trajectory of the pathogen should roughly produce the expected generation time)
  
  
  #(1) First, for each person, we draw the number of possible cases from R0 - this equates to individual heterogeneity in infectiousness 
  # (one type of superspreading) and is already captured in the "total_potential_cases_caused" column, which then gets reduced down to the 
  # proportion in the UCB community in the "original_potential_cases_caused_UCB" column
  #(2) Then, we draw a number of contact events, among which the above cases get distributed. (this equates to event-based superspreading
  # - fewer event draws and a high number of transmissions from #1 generate the biggest superspreading events). Current, we draw this 
  # from a Poisson with lambda=3
  #(3) Then, for each "event", we draw a time that this event took place (here, represented from the generation time Weibull, though this could change) 
  #(4) Then, for each event + time, we go into individual's titer trajectory to determine if each transmission actually 
  # takes place, based on the person's titer load at the point of infection. Since our initial R0 is 2.5, we fix theta at .7, such that the max
  # probability of infection taking place is ~50% at peak viral load. If one 'event' generates multiple cases, each case is treated independently 
  # with this titer-based transmission probability.
  #(5) If there is a group size limit, it gets imposed here. Say that group limit is 6 and one event is supposed to generate 10 cases.
  # If this person abides by group limits (there is a parameter for this), we truncate the 10 person event to a 6 person event, and assume
  # as a worst-case scenario that all 6 of those people get infected 
  
  
  
  #first, draw number of transmission events per person
  pop.mat$num_infection_events <- eventFn(pop.UCB)
  
  #then get a list of event times per person for each of these events
  pop.list <- dlply(pop.mat, .(employ_ids))
  event.times.list <- lapply(pop.list, get.event.time, genTime=genTime)
  
  # now, each person has a number of cases, a number of events, a time for each event, 
  # and a virus titer trajectory.
  # take this information and determine which events actually take place and when they occur
  # also, if applicable, here impose the group size limit and record cases both before 
  # and after that limit occurs
  # return the data as well as the edited event times list that replaces each
  # failed case generation with NA
  
  
  
  double.list <- mapply(FUN=get.real.cases, pop.dat=pop.list, event.dat=event.times.list, titer.dat1 = titer.split, MoreArgs = list(within.host.theta=within.host.theta, group.limit=as.numeric(pop.par$par1[pop.par$parameter=="group-size-limit"])), SIMPLIFY = FALSE)
  
  pop.mat.list <- sapply(double.list, "[",1)
  pop.mat <- data.table::rbindlist(pop.mat.list)
  #pop.mat <- do.call("rbind", pop.mat.list)
  
  gen_time_list <- sapply(double.list, "[",2)
  dat.gen <- data.table::rbindlist(gen_time_list)
  #dat.gen = do.call("rbind",   gen_time_list)
  
  
  #now cases from potential get distributed among events
  #then we determine how many take place based on titers
  #then we remove those that don't take place based on group size limitation
  
  
  
  #then, we set an exposure time for those cases that actually occur
  pop.mat$exposure_time[pop.mat$state>0] <- 0
  
  
  #first, assume that isolation time is symptomatic
  pop.mat$time_isolation[pop.mat$state==1 ] <- as.numeric(pop.mat$time_of_symptom_iso[pop.mat$state==1 ])
  pop.mat$time_isolation = as.numeric(pop.mat$time_isolation)
  pop.mat$reason_isolated[pop.mat$state==1 ] <- "symptom_iso"
  
  #now, if testing (and, for other cases, tracing) comes first, we replace it
  #test needs to be AFTER start time of test sensitive and before end time of test sensitive
  pop.mat$reason_isolated[pop.mat$state==1 & pop.mat$time_isolation> pop.mat$time_of_testing_iso & pop.mat$testing==TRUE & pop.mat$time_of_next_test>=pop.mat$time_test_sensitive_start & pop.mat$time_of_next_test<pop.mat$time_test_sensitive_end] <- "testing_iso"
  pop.mat$time_isolation[pop.mat$state==1 & pop.mat$time_isolation> pop.mat$time_of_testing_iso & pop.mat$testing==TRUE & pop.mat$time_of_next_test>=pop.mat$time_test_sensitive_start & pop.mat$time_of_next_test<pop.mat$time_test_sensitive_end & complete.cases(pop.mat)] <- pop.mat$time_of_testing_iso[pop.mat$state==1 & pop.mat$time_isolation > pop.mat$time_of_testing_iso & pop.mat$testing==TRUE & pop.mat$time_of_next_test>=pop.mat$time_test_sensitive_start & pop.mat$time_of_next_test<pop.mat$time_test_sensitive_end]
  
  
  #and then, if any of these are Inf, change the reason to NA
  pop.mat$reason_isolated[pop.mat$time_isolation==Inf] <- NA
  
  #now, based on isolation time and the generation times in the vector, determine the actual number of cases caused,
  #then export, regroup with other half of population and assign those new infections in the next time step
  
  new.cases = dlply(pop.mat[pop.mat$state==1& pop.mat$potential_cases_caused>0,], .(employ_ids))
  
  
  #if potential cases were 0, then actual cases are too:
  pop.mat$actual_cases_caused[pop.mat$state==1 & pop.mat$potential_cases_caused ==0] <- 0
  pop.mat$actual_cases_caused[pop.mat$state==1 & pop.mat$potential_cases_caused > 0] <-  c(unlist(lapply(new.cases, get.actual.cases, dat.gen=dat.gen, timestep)))
  
  
  #now pop it back out, join with other sub-mats and assign those infections in time and space using your generation time vector.
  
  return(list(pop.mat,  dat.gen))
}
epidemic.step = function(pop.mat, timestep,  length_timestep, prob.out,  gen_list, input.par){
  #pop.mat <- as.data.frame(pop.mat)
  pop.par = subset(input.par, population ==unique(pop.mat$employ_cat))
  #advance timestep
  pop.mat$timestep = timestep
  
  #introduce outside infections into susceptible spaces (cannot misplace those "exposed" by UCB above since we are guaranteeing those transmissions to take place)
  #could easily modulate this for risk cohorts in future
  #check if weekend
  # we say days 1 and 2 are testing
  # days 
  ###MULTIPLE
  # if(timestep ==1 | timestep ==2 | (timestep%%7==1)| (timestep%%7==2)){
  #  n.outside.exposures = sum(sample(x=c(0,1), size=length(pop.mat$state[pop.mat$state==0]), replace=TRUE, prob = c(1-prob.out.wk, prob.out.wk)))
  #}else{
  n.outside.exposures = sum(sample(x=c(0,1), size=length(pop.mat$state[pop.mat$state==0]), replace=TRUE, prob = c(1-prob.out, prob.out)))
  #}
  
  
  if(n.outside.exposures>0){ 
    #if you find some, fill them in with an exposure time of now, distributed at random
    #could add in higher introduction rate for certain sub-groups in this case
    
    new.case.ids = sample(pop.mat$employ_ids[pop.mat$state==0], size = n.outside.exposures, replace=FALSE)
    # print(new.case.ids)
    #and assign
    for (i in 1:length(new.case.ids)){
      #print(i)
      #print(new.case.ids[i])
      #expose the new cases immediately - but only those that have reached the current timestep already
      #those "predestined" for exposure get passed over for now.
      pop.mat$state[pop.mat$employ_ids==new.case.ids[i]] <- 1
      #pop.mat$state[pop.mat$employ_ids==new.case.ids[i] & pop.mat$stat_asymp==TRUE] <- 2
      #exposure time is this timestep
      pop.mat$exposure_time[pop.mat$employ_ids==new.case.ids[i]] <- timestep #infection kicks off so you can now calculat symptom onset time
      #pop.mat$time_of_symptom_onset[pop.mat$employ_ids==new.case.ids[i]] <- pop.mat$time_of_symptom_onset[pop.mat$employ_ids==new.case.ids[i]] + timestep
      #tmp <- pop.mat$time_of_test_positivity[pop.mat$employ_ids==new.case.ids[i]]
      #pop.mat$time_test_sensitive_end[pop.mat$employ_ids==new.case.ids[i]] <- pop.mat$exposure_time[pop.mat$employ_ids==new.case.ids[i]] + pop.mat$time_test_sensitive_end[pop.mat$employ_ids==new.case.ids[i]] 
      #pop.mat$time_test_sensitive_start[pop.mat$employ_ids==new.case.ids[i]] <- pop.mat$exposure_time[pop.mat$employ_ids==new.case.ids[i]] + pop.mat$time_test_sensitive_start[pop.mat$employ_ids==new.case.ids[i]] 
      #infector is outside source that cannot be tracked
      pop.mat$infector[pop.mat$employ_ids==new.case.ids[i]] <- 0
      pop.mat$infector_iso_time[pop.mat$employ_ids==new.case.ids[i]] <- Inf
      pop.mat$case_source[pop.mat$employ_ids==new.case.ids[i]] <- "alameda"
      #introduced cases cannot be traced
      pop.mat$time_of_tracing_iso[pop.mat$employ_ids==new.case.ids[i]] <- Inf 
    }
  }
  
  #pop.mat <- subset(pop.mat, !is.na(employ_ids))
  #now, for those that are currently exposed (both from outside and UCB), 
  # compute distributions of iso time 
  # and new actual cases caused.
  #then, we can assign those times and move them to recovered status
  
  
  pop.mat.old <- pop.mat
  pop.mat <- dplyr::select(pop.mat, -(actual_cases_caused))
  #first, go ahead and move test postivity to the appropriate degree
  #pop.mat$time_of_test_positivity[pop.mat$state==1 | pop.mat$state==2] <- pop.mat$time_of_test_positivity[pop.mat$state==1 | pop.mat$state==2]  + pop.mat$exposure_time[pop.mat$state==1 | pop.mat$state==2] 
  #print("7")
  
  #first, assume that isolation time is symptomatic
  pop.mat$time_isolation[pop.mat$state==1 ] <- as.numeric(pop.mat$time_of_symptom_iso[pop.mat$state==1])
  pop.mat$time_isolation = as.numeric(pop.mat$time_isolation)
  pop.mat$reason_isolated[pop.mat$state==1 ] <- "symptom_iso"
  
  #now, if tracing comes first, we replace it
  #tracing only applicable within our community
  #print("8")
  pop.mat$reason_isolated[pop.mat$state==1 & pop.mat$time_isolation> pop.mat$time_of_tracing_iso & pop.mat$traced==TRUE  & complete.cases(pop.mat)  ] <- "tracing_iso"
  pop.mat$time_isolation[pop.mat$state==1 & pop.mat$time_isolation> pop.mat$time_of_tracing_iso & pop.mat$traced==TRUE  & complete.cases(pop.mat) ] <-  pop.mat$time_of_tracing_iso[pop.mat$state==1 & pop.mat$time_isolation> pop.mat$time_of_tracing_iso & pop.mat$traced==TRUE  & complete.cases(pop.mat)] 
  
  #or, finally, if testing comes first, we replace it here - IF the infection is test sensitive at the time of testing
  
  
  #print("9")
  pop.mat$reason_isolated[pop.mat$state==1 & pop.mat$time_isolation> pop.mat$time_of_testing_iso & pop.mat$testing==TRUE  & complete.cases(pop.mat) & pop.mat$time_of_next_test>=pop.mat$time_test_sensitive_start & pop.mat$time_of_next_test<pop.mat$time_test_sensitive_end] <- "testing_iso"
  pop.mat$time_isolation[pop.mat$state==1 & pop.mat$time_isolation> pop.mat$time_of_testing_iso & pop.mat$testing==TRUE & complete.cases(pop.mat) & pop.mat$time_of_next_test>=pop.mat$time_test_sensitive_start & pop.mat$time_of_next_test<pop.mat$time_test_sensitive_end] <- pop.mat$time_of_testing_iso[pop.mat$state==1 & pop.mat$time_isolation> pop.mat$time_of_testing_iso & pop.mat$testing==TRUE & complete.cases(pop.mat) & pop.mat$time_of_next_test>=pop.mat$time_test_sensitive_start & pop.mat$time_of_next_test<pop.mat$time_test_sensitive_end] 
  
  
  
  #and then, if any of these are Inf, change the reason to NA
  pop.mat$reason_isolated[pop.mat$time_isolation==Inf] <- NA
  
  #first, double-check that nothing was exposed after the isolation time  (would be based on tracing only) 
  #if that ever happens, that person becomes susceptible again because that infection was never generated
  #first flag
  #then, go in and find that person's infector and reduce their actual cases by one
  #based on this case that did not occur
  pop.mat$state[pop.mat$exposure_time>pop.mat$time_isolation & pop.mat$state==1  & pop.mat$reason_isolated=="tracing_iso" ] <- 7
  pop.mat$reason_isolated[pop.mat$state==7] <- NA
  pop.mat$time_isolation[pop.mat$state==7] <- NA
  pop.mat$case_source[pop.mat$state==7 ] <- NA
  
  #now remove a case from the infectors that "caused" these events
  infector.sub1 = pop.mat[pop.mat$state==7,]
  infector.sum1 = ddply(infector.sub1, .(infector), summarize, cases_removed = length(employ_ids)) #save this for the end
  
  
  pop.mat$infector[pop.mat$state==7] <- NA
  pop.mat$infector_iso_time[pop.mat$state==7] <- NA
  pop.mat$exposure_time[pop.mat$state==7]<- NA
  pop.mat$state[pop.mat$state==7] <- 0 
  
  
  #now, based on isolation time and the generation times in the vector, determine the actual number of cases caused,
  #then export, and assign those new infections in the next time step
  
  
  #now, advance forward all of the "time of etc." for susceptibles
  #and time of next testing for all
  
  
  new.cases = dlply(pop.mat[pop.mat$state==1& pop.mat$potential_cases_caused>0,], .(employ_ids))
  
  #if potential cases were 0, then actual cases are too:
  #pop.mat$actual_cases_caused[pop.mat$state==1 & pop.mat$potential_cases_caused ==0| pop.mat$state==2 & pop.mat$potential_cases_caused ==0] <- 0
  
  #but, if potential cases were greater than 0, then actual might be as well, depending on the isolation times
  #dat.gen.new = do.call("rbind", gen_list)
  dat.gen.new = data.table::rbindlist(gen_list)
  #pop.mat$actual_cases_caused[pop.mat$state==1 & pop.mat$potential_cases_caused > 0| pop.mat$state==2 & pop.mat$potential_cases_caused > 0] <-  c(unlist(lapply(new.cases, get.actual.cases, dat.gen=dat.gen.new, timestep=timestep, weekend.amp=weekend.amp)))
  new.actual.cases <-   c(unlist(lapply(new.cases, get.actual.cases, dat.gen=dat.gen.new, timestep=timestep)))
  
  
  #these have not kicked off, so let them kick forward  
  pop.mat$time_of_symptom_onset[pop.mat$state==0 | pop.mat$state==3] <- pop.mat$time_of_symptom_onset[pop.mat$state==0 | pop.mat$state==3] + length_timestep
  pop.mat$time_of_symptom_iso[pop.mat$state==0 | pop.mat$state==3 ] <- pop.mat$time_of_symptom_iso[pop.mat$state==0 | pop.mat$state==3] + length_timestep
  
  pop.mat$time_test_sensitive_start[pop.mat$state==0 | pop.mat$state==3 ] <- pop.mat$time_test_sensitive_start[pop.mat$state==0 | pop.mat$state==3] + length_timestep
  pop.mat$time_test_sensitive_end[pop.mat$state==0 | pop.mat$state==3 ] <- pop.mat$time_test_sensitive_end[pop.mat$state==0 | pop.mat$state==3] + length_timestep
  #tracing only gets started when infector iso time is assigned, so we don't touch it here
  
  #if you are at your current testing date, then next test is bumped into the future. 
  #Otherwise, you just advance in time until you reach it
  
  #but the lag time is maintained after the new test date, so deal with that first
  pop.mat$time_of_testing_iso = pop.mat$time_of_testing_iso - pop.mat$time_of_next_test #now this is just the lag time
  
  #now, compute actual next test day if today is the test day of the runs in question - add different frequencies depending on the type
  
  pop.mat$time_of_next_test[pop.mat$time_of_next_test==timestep] <- timestep + as.numeric(pop.par$par1[pop.par$parameter=="test-freq"])  
  
  
  
  #now put the lag back on to the new test day for isolation
  pop.mat$time_of_testing_iso <- pop.mat$time_of_testing_iso + pop.mat$time_of_next_test
  pop.mat$time_of_testing_iso[pop.mat$time_of_next_test==Inf] <- Inf
  
  
  #and, finally, check in on those that were "pre-exposed" up above.
  #move them up to their appropriate status if they should be exposed now
  #if they reach it, go ahead and assign their actual cases
  
  #first, eliminate if they should not occur
  #first flag
  #then, go in and find that person's infector and reduce their actual cases by one
  #based on this case that did not occur
  pop.mat$state[pop.mat$state==3 & pop.mat$exposure_time>pop.mat$infector_iso_time ] <- 8
  pop.mat$time_of_tracing_iso[pop.mat$state==8 & pop.mat$exposure_time>pop.mat$infector_iso_time ] <- pop.mat$time_of_tracing_iso[pop.mat$state==8 & pop.mat$exposure_time>pop.mat$infector_iso_time  ] - pop.mat$infector_iso_time[pop.mat$state==8 & pop.mat$exposure_time>pop.mat$infector_iso_time]
  pop.mat$case_source[pop.mat$state==8 & pop.mat$exposure_time>pop.mat$infector_iso_time] <- NA
  
  #now remove a case from the infectors that "caused" these events
  infector.sub2 = pop.mat[pop.mat$state==8,]
  infector.sum2 = ddply(infector.sub2, .(infector), summarize, cases_removed = length(employ_ids)) #save this for the end
  
  pop.mat$infector[pop.mat$state==8 & pop.mat$exposure_time>pop.mat$infector_iso_time ] <- NA
  pop.mat$state[pop.mat$state==8 & pop.mat$exposure_time>pop.mat$infector_iso_time ] <- 0
  pop.mat$infector_iso_time[pop.mat$state==0] <- NA
  pop.mat$exposure_time[pop.mat$state==0] <- NA
  
  if (exists('infector.sum1') &  exists('infector.sum2')){
    infector.sum <- rbind(infector.sum1, infector.sum2)  
  }else if(exists('infector.sum1')){
    infector.sum <- infector.sum1
  }else if(exists('infector.sum2')){
    infector.sum <- infector.sum2
  }
  
  
  #then, if they pass that test and still remain 'pre-exposed', check to see if they should be elevated in status to 1 or 2 
  #(meaning they have reached the exposure time)
  #if so, assign them isolation time and actual cases which get allocated in the next round. 
  #otherwise, they just keep current status as "pre-exposed"
  
  #first make them complete cases, so that R is not angry with new columns being filled in
  pop.mat$time_isolation[pop.mat$state==3  & pop.mat$exposure_time<=timestep] <- Inf
  pop.mat$reason_isolated[pop.mat$state==3  & pop.mat$exposure_time<=timestep] <- "in progress"
  
  
  #first, assume that isolation time is symptomatic
  #print("1")
  pop.mat$reason_isolated[pop.mat$state==3  & pop.mat$exposure_time<=timestep  & complete.cases(pop.mat)] <- "symptom_iso"
  # print("2")
  pop.mat$time_isolation[pop.mat$state==3  & pop.mat$exposure_time<=timestep & complete.cases(pop.mat)] <- pop.mat$time_of_symptom_iso[pop.mat$state==3  & pop.mat$exposure_time<=timestep & complete.cases(pop.mat)]
  pop.mat$time_isolation = as.numeric(pop.mat$time_isolation)
  
  
  #now, if tracing comes first, we replace it
  #print("3")
  pop.mat$reason_isolated[pop.mat$state==3  & pop.mat$exposure_time<=timestep & pop.mat$time_of_tracing_iso<pop.mat$time_isolation & pop.mat$traced==TRUE  & complete.cases(pop.mat)] <- "tracing_iso"
  # print("4")
  pop.mat$time_isolation[pop.mat$state==3  & pop.mat$exposure_time<=timestep & pop.mat$time_of_tracing_iso<pop.mat$time_isolation & pop.mat$traced==TRUE  & complete.cases(pop.mat)] <- pop.mat$time_of_tracing_iso[pop.mat$state==3  & pop.mat$exposure_time<=timestep & pop.mat$time_of_tracing_iso<pop.mat$time_isolation & pop.mat$traced==TRUE  & complete.cases(pop.mat)]
  
  
  #finally, if testing comes first, we replace it
  #print("5")
  pop.mat$reason_isolated[pop.mat$state==3  & pop.mat$exposure_time<=timestep & pop.mat$time_of_testing_iso<pop.mat$time_isolation & pop.mat$testing==TRUE & pop.mat$time_of_next_test>=pop.mat$time_test_sensitive_start & pop.mat$time_of_next_test<pop.mat$time_test_sensitive_end  & complete.cases(pop.mat)] <- "testing_iso"
  #print("6")
  
  #print(pop.mat[pop.mat$state==3  & pop.mat$exposure_time<=timestep & pop.mat$time_of_testing_iso<pop.mat$time_isolation & pop.mat$testing==TRUE  & complete.cases(pop.mat) ,])
  #print(pop.mat)
  pop.mat$time_isolation[pop.mat$state==3  & pop.mat$exposure_time<=timestep & pop.mat$time_of_testing_iso<pop.mat$time_isolation & pop.mat$testing==TRUE & pop.mat$time_of_next_test>=pop.mat$time_test_sensitive_start & pop.mat$time_of_next_test<pop.mat$time_test_sensitive_end  & complete.cases(pop.mat)] <- pop.mat$time_of_testing_iso[pop.mat$state==3  & pop.mat$exposure_time<=timestep & pop.mat$time_of_testing_iso<pop.mat$time_isolation & pop.mat$testing==TRUE &  pop.mat$time_of_next_test>=pop.mat$time_test_sensitive_start & pop.mat$time_of_next_test<pop.mat$time_test_sensitive_end  & complete.cases(pop.mat)]
  
  
  #and then, if any of these are Inf, change the reason to NA
  pop.mat$reason_isolated[pop.mat$time_isolation==Inf] <- NA
  
  #now, based on isolation time and the generation times in the vector, determine the actual number of cases caused,
  #then export, regroup with other half of population and assign those new infections in the next time step
  new.cases = dlply(pop.mat[pop.mat$state== 3& pop.mat$potential_cases_caused>0 & pop.mat$exposure_time<=timestep & complete.cases(pop.mat) ,], .(employ_ids))
  
  #if potential cases were 0, then actual cases are too:
  #pop.mat$actual_cases_caused[pop.mat$state==3 & pop.mat$potential_cases_caused ==0 & pop.mat$exposure_time<=timestep  | pop.mat$state==4 & pop.mat$potential_cases_caused ==0 & pop.mat$exposure_time<=timestep] <- 0
  
  tmp.dat =  pop.mat[pop.mat$state==3 & pop.mat$potential_cases_caused >0 & pop.mat$exposure_time<=timestep  ,]
  new.cases.2 <- dlply(tmp.dat, .(employ_ids))
  
  new.actual.cases.3  <-  c(unlist(lapply(new.cases.2, get.actual.cases, dat.gen=dat.gen.new, timestep=timestep)))
  
  #now add the actual cases back in 
  pop.mat <- cbind.data.frame(pop.mat, pop.mat.old$actual_cases_caused)
  names(pop.mat)[length(names(pop.mat))] <- "actual_cases_caused"
  #reorder
  pop.mat <- dplyr::select(pop.mat, names(pop.mat.old))
  
  #and add in the new actual cases
  pop.mat$actual_cases_caused[pop.mat$state==3 & pop.mat$potential_cases_caused ==0 & pop.mat$exposure_time<=timestep  | pop.mat$state==1 & pop.mat$potential_cases_caused ==0 & pop.mat$exposure_time<=timestep ] <- 0
  pop.mat$actual_cases_caused[pop.mat$state==1 & pop.mat$potential_cases_caused >0 & pop.mat$exposure_time<=timestep] <- new.actual.cases
  pop.mat$actual_cases_caused[pop.mat$state==3 & pop.mat$potential_cases_caused >0 & pop.mat$exposure_time<=timestep ] <- new.actual.cases.3
  
  
  #and, finally, change state so these cases can get allocated in the next round.
  pop.mat$state[pop.mat$state==3 & pop.mat$exposure_time<=timestep  ] <- 1
  
  
  #and remove any avoided cases if there were some
  if (exists('infector.sum')){
    if(nrow(infector.sum)>0){
      for(i in 1:length(infector.sum$infector)){
        pop.mat$actual_cases_caused[pop.mat$employ_ids==infector.sum$infector[i]] <-   pop.mat$actual_cases_caused[pop.mat$employ_ids==infector.sum$infector[i]] - infector.sum$cases_removed[i]  
      }
    }
  }
  
  
  #and return
  return(pop.mat)
  
}
get.mean.sd <- function(vector, name){
  
  #first, trim to same length
  min.length <- min(unlist(lapply(vector, length)))
  
  for (i in 1:length(vector)){
    vector[[i]] <- vector[[i]][1:min.length]
  }
  vec <- unlist(vector, use.names = FALSE)
  DIM <- c(length(vector[[1]]),1)
  n <- length(vector)
  
  list.mean <- tapply(vec, rep(1:prod(DIM),times = n), mean)
  attr(list.mean, "dim") <- DIM
  list.mean <- as.data.frame(list.mean)
  
  list.sd <- tapply(vec, rep(1:prod(DIM),times = n), sd)
  attr(list.sd, "dim") <- DIM
  list.sd <- as.data.frame(list.sd)
  
  list.uci = list.mean + 1.96*list.sd
  list.lci = list.mean - 1.96*list.sd
  
  list.lci[list.lci<0] <- 0
  list.uci[list.uci<0] <- 0
  
  dat= cbind.data.frame(list.mean, list.lci, list.uci)
  names(dat) = paste(c("mean", "lci", "uci"), name, sep="_")
  return(dat)
} 
get.mean.matrix <- function(mat){
  
  
  #first, trim to same length
  min.length <- min(unlist(lapply(mat, nrow)))
  
  n.cat = ncol(mat[[1]])/3
  
  
  for (i in 1:length(mat)){
    mat[[i]] <- mat[[i]][1:min.length,]
  }
  
  list.mean <- Reduce("+",mat) / length(mat)
  mat.2 <-  do.call("cbind", mat)
  #mat.2 <- data.table::rbindlist(mat)
  list.sd <- apply(mat.2, 1, sd)
  
  list.uci = list.mean + 1.96*list.sd
  list.lci = list.mean - 1.96*list.sd
  
  list.lci[list.lci<0] <- 0
  list.uci[list.uci<0] <- 0
  
  dat= cbind.data.frame(list.mean, list.lci, list.uci)
  
  names(dat) = c(paste0("mean_iso_cat_",seq(1,n.cat,1)), paste0("lci_iso_cat_",seq(1,n.cat,1)), paste0("uci_iso_cat_",seq(1,n.cat,1)), 
                 paste0("mean_exp_cat_",seq(1,n.cat,1)), paste0("lci_exp_cat_",seq(1,n.cat,1)), paste0("uci_exp_cat_",seq(1,n.cat,1)),
                 paste0("mean_deaths_cat_",seq(1,n.cat,1)), paste0("lci_deaths_cat_",seq(1,n.cat,1)), paste0("uci_deaths_cat_",seq(1,n.cat,1)))
  return(dat)
} 
convert.cat = function(dat){
  n.cat = ncol(dat)/9
  max.times = nrow(dat)
  iso.dat = dat[,1:(n.cat*3)]
  exp.dat = dat[,(n.cat*3+1):(n.cat*3*2)]
  death.dat = dat[,(n.cat*3*2+1):ncol(dat)]
  
  
  #then, sep by cat
  list.iso  <- list.exp <- list.deaths  <-   list() 
  
  
  for(i in 1:n.cat){
    list.iso[[i]] <- cbind.data.frame(iso.dat[,i], iso.dat[,i+n.cat],iso.dat[,i+(n.cat*2)])  
    list.exp[[i]] <- cbind.data.frame(exp.dat[,i],exp.dat[,i+n.cat],exp.dat[,i+(n.cat*2)])  
    list.deaths[[i]] <- cbind.data.frame(death.dat[,i], death.dat[,i+n.cat],death.dat[,i+(n.cat*2)])  
  }
  
  #iso.db <- do.call("rbind", list.iso)
  iso.db <- data.table::rbindlist(list.iso)
  iso.db$type = rep(1:n.cat, each = max.times)
  iso.db$type <- paste0("iso-pop-", iso.db$type)
  
  names(iso.db) <- c("mean", "lci", "uci", "type")
  
  #exp.db <- do.call("rbind", list.exp)
  exp.db <- data.table::rbindlist(list.exp)
  exp.db$type = rep(1:n.cat, each = max.times)
  exp.db$type <- paste0("exp-pop-", exp.db$type)
  
  names(exp.db) <- c("mean", "lci", "uci", "type")
  
  #death.db <- do.call("rbind", list.deaths)
  death.db <- data.table::rbindlist(list.deaths)
  death.db$type = rep(1:n.cat, each = max.times)
  death.db$type <- paste0("death-pop-", death.db$type)
  
  names(death.db) <- c("mean", "lci", "uci", "type")
  
  return(list(iso.db, exp.db, death.db))
}
R.fit.sum <- function(mat.df){
  
  #apply across all columns
  mean.all <- apply(mat.df, 2,mean)
  sd.all <- apply(mat.df, 2,sd)
  lci.all <- mean.all-1.96*sd.all
  lci.all[  lci.all < 0] <- 0
  uci.all <- mean.all+1.96*sd.all
  
  #and nbinom fit
  all.fit <-  apply(mat.df, 2, fitdist, distr="nbinom")
  
  #and return
  out.dat <- cbind.data.frame(mean.all, lci.all, uci.all)
  
  out.dat$class <-   names(mat.df)
  #names(out.dat) <- names(mat.df)
  #out.dat$estimate <- c("mean", "lci", "uci")
  #out.dat[out.dat<0] <- 0
  
  #and add fit
  size.out <- list()
  mu.out <- list()
  for(i in 1:length(all.fit)){
    size.out[[i]] <- all.fit[[i]]$estimate[1]
    mu.out[[i]] <- all.fit[[i]]$estimate[2]
    
  }
  
  
  size.out <- c(unlist(size.out))
  mu.out <- c(unlist(mu.out))
  out.dat$nb_mu <- mu.out
  out.dat$nb_size <- size.out
  
  
  # names(size.out) <- names(mu.out) <- names(out.dat)
  # out.dat <- rbind(out.dat, size.out, mu.out)
  # 
  # out.dat$total_potential_cases <- as.numeric(out.dat$total_potential_cases)
  # out.dat$UCB_potential_cases  <-  as.numeric(out.dat$UCB_potential_cases)
  # out.dat$UCB_post_group_potential_cases <-  as.numeric(out.dat$UCB_post_group_potential_cases)
  # out.dat$UCB_post_titer_potential_cases <-  as.numeric(out.dat$UCB_post_titer_potential_cases)
  # out.dat$UCB_post_isolations_actual_cases <-  as.numeric(out.dat$UCB_post_isolations_actual_cases)
  # 
  return(out.dat)  
}
R.fit.sum.lognorm <- function(mat.df){
  
  #apply across all columns
  mean.all <- apply(mat.df, 2,mean)
  sd.all <- apply(mat.df, 2,sd)
  lci.all <- mean.all-1.96*sd.all
  lci.all[  lci.all < 0] <- 0
  uci.all <- mean.all+1.96*sd.all
  
  #and return
  out.dat <- cbind.data.frame(mean.all, lci.all, uci.all)
  
  out.dat$class <-   names(mat.df)
  #names(out.dat) <- names(mat.df)
  #out.dat$estimate <- c("mean", "lci", "uci")
  #out.dat[out.dat<0] <- 0
  
  
  # 
  return(out.dat)  
}
simulate.epidemic <- function(input.pop, n.init.exposed.vector, employ.id.vector, times, virus.par, input.par, burnin,
                              test.freq, length_timestep, bay.area.prev, initial.R, within.host.theta,  titer.dat, LOD, test_rotation_name){
  
  
  
  
  if (virus.par$distribution[virus.par$parameter=="R0"]=="log-normal"){
    
    #sample R0 normal
    R0fn = R0_fn(meanR0=virus.par$par1[virus.par$parameter=="R0"],
                 sdR0=virus.par$par2[virus.par$parameter=="R0"])
    
    
  }else if(virus.par$distribution[virus.par$parameter=="R0"]=="negbinom"){
    
    #sample R0 normal
    R0fn = R0_fn_nb(muR0=virus.par$par1[virus.par$parameter=="R0"],
                    sizeR0=virus.par$par2[virus.par$parameter=="R0"])
    
  }
  
  #and the number of transmission events, from a negbinom
  #remember that fewer events = higher likelihood of a big superspreading event.
  #but the vast majority of people have both few events and few cases
  
  eventFn = poisson_fn(lambda =as.numeric(input.par$par1[input.par$parameter=="transmission-events"]))
  
  
  
  
  #and normal distribution of the detection limit
  
  #then, form your new populations
  #now split the population based on risk
  tot.pop = length(input.pop)
  pop.num = 1:tot.pop
  
  titer.dat$cat <- NA 
  for (i in 1:(length(pop.num)-1)){
    titer.dat$cat[titer.dat$employ_ids < employ.id.vector [i+1]  & titer.dat$employ_ids >= employ.id.vector [i]] <-   pop.num[i]
  }
  
  titer.dat$cat[is.na(titer.dat$cat)] <- pop.num[length(pop.num)]
  
  #and split 
  titer.dat.split <- dlply(titer.dat, .(cat))
  
  
  
  #make the proper number of pop.mat depending on the total number of subpopulations
  #populate each using the appropriate parameters
  out.list = mapply(FUN=initiate.pop, start.ID.employ = as.list(employ.id.vector), pop.UCB=as.list(input.pop), n.init.exposed= as.list(n.init.exposed.vector),  pop.ID = as.list(pop.num), titer.dat=titer.dat.split, 
                    MoreArgs= list(input.par=input.par, virus.par=virus.par,  R0fn=R0fn,   eventFn=eventFn, within.host.theta=within.host.theta, LOD=LOD))
  
  
  pop.list = out.list[1,]
  gen_list_long <- out.list[2,]
  #original.r0 <- out.list[3,][[1]]
  #gen_list_long_wkend <- out.list[3,]
  
  #pop.mat <- do.call("rbind", pop.list)
  pop.mat <- data.table::rbindlist(pop.list)
  #gen.dat.all <- do.call("rbind", gen_list_long)
  gen.dat.all <- data.table::rbindlist(gen_list_long)
  
  
  #now, double-check that the generation time dataframe is the same length as the number of unique employ ids
  if(sum(setdiff(pop.mat$employ_ids, gen.dat.all$employ_ids))>0){
    
    missing.ids <- setdiff(pop.mat$employ_ids, gen.dat.all$employ_ids)
    
    missing.cases <- list()
    for(i in 1:length(missing.ids)){
      missing.cases[[i]] <- pop.mat$potential_cases_caused[pop.mat$employ_ids==missing.ids[i]] 
    }
    missing.cases <- c(unlist(missing.cases))
    
    if(sum(missing.cases)>0){
      missing.gen <- genTime(missing.cases)  
      add.dat <- cbind.data.frame(rep(missing.ids, missing.cases), missing.gen)
    }else{
      missing.gen <- rep(NA, length(missing.cases))  
      add.dat <- cbind.data.frame(missing.ids, missing.gen)
    }
    
    
    names(add.dat) <- names(gen.dat.all)
    
    gen.dat.all <- rbind(gen.dat.all, add.dat)
    
    gen.dat.all <- arrange(gen.dat.all, employ_ids)
    
  }
  
  
  gen_list =  dlply(gen.dat.all, .(employ_ids))
  #gen_list_wk =  dlply(gen.dat.all.wk, .(employ_ids))
  
  foi.bay.area = initial.R*bay.area.prev*length_timestep #rate per day at which susceptibles become infected
  #foi.wkend = bay.area.R*bay.area.prev*length_timestep*weekend.amp
  prob.outside.exposure =1-(exp(-1*foi.bay.area)) #for each person in berkeley, this is the probability of getting exposed each day
  prob.outside.exposure[prob.outside.exposure<0] <- 0
  #prob.outside.exposure.wk =1-(exp(-1*foi.wkend))
  #could also be a vector
  
  times_vect = seq(length_timestep,times, by = length_timestep)
  
  for(i in 1: length(times_vect)){
    
    #print(i)
    timestep = times_vect[i]
    #could make other functions here if people mostly infect their own subgroups
    #here, we distribute the infections amongst new people and retire the old
    pop.mat = assign.infections(pop.mat = pop.mat, gen_list=gen_list, timestep = timestep, input.par = input.par)
    
    #now split it by population to introduce outside exposures
    pop.split = dlply(pop.mat, .(employ_cat))
    
    pop.mat.list = lapply(pop.split, FUN=epidemic.step, timestep= timestep, prob.out = prob.outside.exposure, gen_list=gen_list,  input.par=input.par, length_timestep = length_timestep)
    
    #then, rejoin
    #pop.mat = do.call("rbind", pop.mat.list)#print(i)
    pop.mat = data.table::rbindlist(pop.mat.list)
    
    #then, just keep tabs that there are enough susceptibles to fill the new cases in the next step
    remaining.susceptibles = length(pop.mat$state[pop.mat$state==0])
    future.cases = sum(pop.mat$actual_cases_caused[pop.mat$state==1])
    
    if(future.cases>remaining.susceptibles){ #if there are not enough susceptibles left for all of the assigned cases before you reach the end of the time series, then you go into the next step
      #print(i)
      pop.mat = assign.last.infections(pop.mat = pop.mat, gen_list = gen_list, remaining.susceptibles = remaining.susceptibles, timestep = timestep)
      #print(i)
    }
    
    
  }
  
  #collect all the "R" reduction info:
  R.mat <- dplyr::select(pop.mat, total_potential_cases_caused, original_potential_cases_caused_UCB, post_titer_potential_cases_caused_UCB, potential_cases_caused, actual_cases_caused)
  names(R.mat) <- c( "total_potential_cases", "UCB_potential_cases", "UCB_post_titer_potential_cases", "UCB_post_group_potential_cases", "UCB_post_isolations_actual_cases")
  
  R.mat <- arrange(R.mat, desc(total_potential_cases))
  R.mat$UCB_post_isolations_actual_cases[is.na(R.mat$UCB_post_isolations_actual_cases)] <- 0
  #R.mat <- as.matrix(R.mat)
  
  # #new R0
  # new.R0 = subset(pop.mat, !is.na(infector))
  # new.R0 = ddply(new.R0, .(infector), summarize, cases_caused=length(employ_ids))
  # tot.introductions = new.R0$cases_caused[new.R0$infector=="0"]
  # new.R0 = subset(new.R0, infector!="0")
  # 
  # maxID = max(pop.mat$employ_ids)
  # missing_ids <- (1:maxID)[!(1:maxID %in% new.R0$infector)]
  # 
  # # add in missing days if any are missing
  # if (length(missing_ids > 0)) {
  #   R0comp <- data.table::rbindlist(list(new.R0,
  #                                        data.table(infector = missing_ids,
  #                                                   cases_caused = 0)))
  # }
  # 
  # R0comp <- arrange(R0comp, infector)
  # 
  # #now add back in those cases not at UCB...
  # #original.r0$actual_cases_caused_UCB <- R0comp$cases_caused
  
  #get prop.asymptomatic at this cutoff
  prop.asym <- length(pop.mat$time_of_symptom_onset[pop.mat$time_of_symptom_onset==Inf])/length(pop.mat$time_of_symptom_iso)
  
  #from here, compute Reffective
  R.dat = dplyr::select(pop.mat, employ_ids, infector, time_isolation, case_source)
  
  R.dat = arrange(R.dat, time_isolation) #icidence will just be cases by time isolated
  #if not isolated, you don't count for incidence...
  R.dat = R.dat[!is.na(R.dat$time_isolation),]
  R.dat$time_isolation = ceiling(R.dat$time_isolation)
  
  #could add source. don't for now
  R.sum = ddply(R.dat, .(time_isolation), summarise, length(employ_ids))
  #R.sum = ddply(R.dat, .(time_isolated, source), summarise, length(employ_ids))
  names(R.sum) = c( "day", "incidence")
  
  #plot as incidence
  #plot(as.incidence(R.sum$incidence, dates = R.sum$day))
  
  #this will go in as your incidence data
  
  #now add in pairs to estimate the serial interval
  #T <- nrow(R.sum)
  #t_start <- seq(2, T-13) # starting at 2 as conditional on the past observations
  #t_end <- t_start + 13 
  # 
  # R.est = estimate_R(R.sum$incidence,
  #                    method="parametric_si",
  #                    config = make_config(list(#t_start = t_start,
  #                      #t_end =   t_end, 
  #                      mean_si = serial_mean, std_si = serial_sd)))
  # 
  # #plot(R.est, "R")
  # #get midpoint and R values and extract
  # R.out = cbind.data.frame(get.midpoint(par.low = R.est$R$t_start, par.hi = R.est$R$t_end), R.est$R$`Mean(R)`)
  # names(R.out) = c("day", "Reffective")
  # #and try it based on pairs
  
  pop.mat = data.table(pop.mat)
  #now, get broad incidence data to report
  UCB.mat = subset(pop.mat, case_source=="UCB")
  alameda.mat = subset(pop.mat, case_source=="alameda")
  
  
  
  symp.mat = subset(pop.mat, reason_isolated=="symptom_iso")
  trace.mat = subset(pop.mat, reason_isolated=="tracing_iso")
  test.mat = subset(pop.mat, reason_isolated=="testing_iso")
  
  
  
  daily_exposures <- pop.mat[, day := ceiling(exposure_time) #time_isolated
                             ][, .(daily_exposures = .N), by = day
                               ]
  
  # #daily isolations 
  daily_isolations <- pop.mat[, day := ceiling(time_isolation) #
                              ][, .(daily_isolations = .N), by = day
                                ]
  
  daily_cal <- UCB.mat[, day := ceiling(time_isolation) #time_isolated
                       ][, .(daily_isolations = .N), by = day
                         ]
  daily_alameda <- alameda.mat[, day := ceiling(time_isolation) #time_isolated
                               ][, .(daily_isolations = .N), by = day
                                 ]
  
  daily_symp <- symp.mat[, day := ceiling(time_isolation) #time_isolated
                         ][, .(daily_isolations = .N), by = day
                           ]
  
  daily_trace <- trace.mat[, day := ceiling(time_isolation) #time_isolated
                           ][, .(daily_isolations = .N), by = day
                             ]
  
  daily_test <- test.mat[, day := ceiling(time_isolation) #time_isolated
                         ][, .(daily_isolations = .N), by = day
                           ]
  
  
  # maximum outbreak day
  max_day <- ceiling(times)
  # days with 0 cases in 0:max_week
  #missing_days <- (0:max_day)[!(0:max_day %in% daily_isolations$day)]
  missing_days <- (0:max_day)[!(0:max_day %in% daily_exposures$day)]
  
  # add in missing days if any are missing
  if (length(missing_days > 0)) {
    daily_cases <- data.table::rbindlist(list(daily_exposures,
                                              data.table(day = missing_days,
                                                         daily_exposures = 0)))
  }
  
  #reorder as appropriate
  #daily_cases <- arrange(daily_cases, day)
  
  
  # order and sum up
  daily_cases <- daily_exposures[order(day)
                                 ][, cumulative := cumsum(daily_exposures)]
  
  
  # cut at max_week
  daily_cases <- daily_cases[day<=max_day]
  
  
  # and isoaltions
  daily_cases$daily_isolations <- 0
  for (i in 1:length(daily_isolations$day)){
    daily_cases$daily_isolations[daily_cases$day==daily_isolations$day[i]] <- daily_isolations$daily_isolations[i]
  }
  
  
  #and cumulative isolations
  daily_cases$cumulative_iso =  cumsum(daily_cases$daily_isolations)
  
  
  # #and cases in UCB vs out
  daily_cases$daily_UCB_isolations <- 0
  for (i in 1:length(daily_cal$day)){
    daily_cases$daily_UCB_isolations[daily_cases$day==daily_cal$day[i]] <- daily_cal$daily_isolations[i]
  }
  # 
  # #and cases in UCB vs out
  daily_cases$daily_alameda_isolations <- 0
  for (i in 1:length(daily_alameda$day)){
    daily_cases$daily_alameda_isolations[daily_cases$day==daily_alameda$day[i]] <- daily_alameda$daily_isolations[i]
  }
  
  daily_cases$daily_symptomatic_isolations <- 0
  for (i in 1:length(daily_symp$day)){
    daily_cases$daily_symptomatic_isolations[daily_cases$day==daily_symp$day[i]] <- daily_symp$daily_isolations[i]
  }
  
  daily_cases$daily_tracing_isolations <- 0
  for (i in 1:length(daily_trace$day)){
    daily_cases$daily_tracing_isolations[daily_cases$day==daily_trace$day[i]] <- daily_trace$daily_isolations[i]
  }
  
  
  daily_cases$daily_testing_isolations <- 0
  for (i in 1:length(daily_test$day)){
    daily_cases$daily_testing_isolations[daily_cases$day==daily_test$day[i]] <- daily_test$daily_isolations[i]
  }
  
  # 
  # #now attach R-effective 
  # daily_cases$Reffective = NA
  # 
  # for(i in 1:nrow(R.out)){
  #   daily_cases$Reffective[daily_cases$day==R.out$day[i]] <- R.out$Reffective[i]
  # }
  # 
  
  #add category
  pop.mat.cat= dlply(pop.mat, .(employ_cat))
  new_col <- lapply(pop.mat.cat, FUN=add.risk.cat, pop_dat=daily_cases)
  
  #and also the daily exposures
  new_col2 <- lapply(pop.mat.cat, FUN=add.risk.cat.exp, pop_dat=daily_cases, input_par=input.par)
  
  new_col_exp <- sapply(new_col2, "[", 1)
  new_col_deaths <- sapply(new_col2, "[", 2)
  
  #tmp = data.table::rbindlist(new_col)
  tmp = as.data.frame(do.call("cbind", new_col))
  names(tmp) <- paste0("isolations-employ-cat-", unique(input.par$population))
  
  tmp2 = as.data.frame(do.call("cbind", new_col_exp))
  #tmp2 = data.table::rbindlist(new_col_exp)
  names(tmp2) <- paste0("exposures-employ-cat-", unique(input.par$population))
  
  tmp3 = as.data.frame(do.call("cbind", new_col_deaths))
  #tmp3 = data.table::rbindlist(new_col_deaths)
  names(tmp3) <- paste0("deaths-employ-cat-", unique(input.par$population))
  
  #and attach to daily cases
  daily_cases <- cbind.data.frame(daily_cases, tmp, tmp2, tmp3)
  
  
  # #finally, calculate some summary statistics from the epidemic
  # tot.exposures = sum(daily_cases$daily_exposures, na.rm=T)
  # tot.isolations = sum(daily_cases$daily_isolations, na.rm=T)
  # #time.to.control = max(daily_cases$day[!is.na(daily_cases$Reffective)])
  # max.exposures.per.day = max(daily_cases$daily_exposures, na.rm=T)
  # mean.exposures.per.day = mean(daily_cases$daily_exposures, na.rm=T)
  # max.iso.per.day = max(daily_cases$daily_isolations, na.rm=T)
  # mean.iso.per.day = mean(daily_cases$daily_isolations, na.rm=T)
  # time.of.peak.iso = min(daily_cases$day[daily_cases$daily_isolations==max(daily_cases$daily_isolations, na.rm=T)])
  # time.of.peak.exposure = min(daily_cases$day[daily_cases$daily_exposures==max(daily_cases$daily_exposures, na.rm=T)])
  # 
  #and report out the max day before your cases are too few to calculate Reffective
  #out.stat <- c(tot.exposures, tot.isolations, max.exposures.per.day, mean.exposures.per.day, max.iso.per.day, mean.iso.per.day, time.of.peak.exposure, time.of.peak.iso) 
  #names(out.stat) <- c("total_exposures", "total_isolations",  "max_exp_per_day", "mean_exp_per_day", "max_iso_per_day", "mean_iso_per_day", "time_peak_exposure", "time_peak_isolation")
  
  pop.mat$LOD <- LOD
  #add TAT if this is a single population model, but if it is mixed in a multipop, note that
  if(length(unique(input.par$par1[input.par$parameter=="TAT-lag"]))==1){
    pop.mat$TAT <- unique(input.par$par1[input.par$parameter=="TAT-lag"])  
  }else{
    pop.mat$TAT <- "multiple"
  }
  
  pop.mat$test_rotation <- test_rotation_name
  
  
  return(list(daily_cases,pop.mat, prop.asym, R.mat))
}
replicate.epidemic = function(n.reps, input.pop, n.init.exposed.vector, employ.id.vector, times, virus.par, input.par, burnin, test.freq, length_timestep,
                              bay.area.prev, initial.R, within.host.theta, test_rotation_name, LOD, titer.dat){
  
  out = replicate(n.reps, simulate.epidemic(virus.par = virus.par,
                                            input.par = input.par,
                                            input.pop=input.pop, 
                                            n.init.exposed.vector=n.init.exposed.vector, 
                                            times=times, 
                                            bay.area.prev = bay.area.prev, 
                                            initial.R = initial.R, 
                                            within.host.theta = within.host.theta,
                                            burnin = burnin,  
                                            length_timestep=length_timestep,
                                            employ.id.vector =employ.id.vector,
                                            LOD = LOD,
                                            titer.dat = titer.dat,
                                            test_rotation_name = test_rotation_name), simplify = "array")
  #make list
  out.time<- out.daily <- out.cal <- out.iso <- out.cumulative <-  out.ala <- out.symp <- out.trace <- out.test  <- out.iso <-out.cum.iso <- pop.mat.chain <- out.prop.asym <- R.mat.out <-  list()
  
  #and make list of all the categories of sub-pop
  out.cat <- list()
  
  for (i in 1:ncol(out)){
    #tmp <- data.table::cbindlist(out[,i][[1]])
    tmp <- do.call("cbind", out[,i][[1]])
    out.time[[i]] <- c(unlist(tmp[,1]))
    out.daily[[i]] <- c(unlist(tmp[,2]))
    out.cumulative[[i]] <- c(unlist(tmp[,3]))
    out.iso[[i]] <- c(unlist(tmp[,4]))
    out.cum.iso[[i]] <- c(unlist(tmp[,5]))
    out.cal[[i]] <- c(unlist(tmp[,6]))
    out.ala[[i]] <- c(unlist(tmp[,7]))
    out.symp[[i]] <- c(unlist(tmp[,8]))
    out.trace[[i]] <- c(unlist(tmp[,9]))
    out.test[[i]] <- c(unlist(tmp[,10]))
    
    
    #out.R[[i]] <- c(unlist(tmp[,11]))
    
    out.cat[[i]] <- cbind(unlist(tmp[,11:(10+(length(unique(input.par$population)))*3)]))
    
    #and save a chain of pop.mat
    tmp2 <- out[,i][[2]]
    pop.mat.chain[[i]] <- tmp2
    
    
    #and the prop.asym
    tmp3 <- out[,i][[3]]
    out.prop.asym[[i]] <- tmp3
    
    tmp4 <- out[,i][[4]]
    rownames(tmp4) <- c()
    R.mat.out[[i]] <- tmp4
    #unique(input.par$population)
    
  }
  
  
  #now shorten them all to the same length and get mean + sd
  
  #print(out.time)
  mean.time = get.mean.sd(vector= out.time, name = "day")[,1]
  #print(out.daily)
  mean.daily = get.mean.sd(vector=out.daily, name = "exposures")
  #print(out.cumulative)
  mean.cumulative=  get.mean.sd(vector=out.cumulative, name = "cumulative")
  #print(out.cal)
  mean.cal = get.mean.sd(vector=out.cal, name="UCB")
  #print(out.ala)
  mean.ala =  get.mean.sd(vector=out.ala, name = "AlamedaCo")
  #print(out.low)
  mean.symp = get.mean.sd(vector=out.symp, name="symptomatic_iso")
  mean.trace = get.mean.sd(vector=out.trace, name="tracing_iso")
  mean.test = get.mean.sd(vector=out.test, name="testing_iso")
  #print(out.iso)
  mean.iso =  get.mean.sd(vector=out.iso, name = "isolations")
  #print(out.cum.iso)
  mean.cum.iso =  get.mean.sd(vector=out.cum.iso, name = "cumulative_isolations")
  #print(out.sum)
  #mean.sum = get.mean.sd.summary(out.sum)
  
  
  #and the employ-cat
  mean.cat = get.mean.matrix(mat=out.cat)
  #print(out.hi)
  
  
  mean.dat = cbind.data.frame(mean.time, mean.daily, mean.cumulative, mean.cal, mean.ala, mean.symp, mean.trace,mean.test, mean.iso, mean.cum.iso, mean.cat)#, mean.R)
  names(mean.dat)[1] = "day"   
  
  #all of the descriptors can now change within the pop
  mean.dat$LOD <- LOD
  
  if(length(unique(input.par$par1[input.par$parameter=="TAT-lag"]))==1){
    mean.dat$TAT <- unique(input.par$par1[input.par$parameter=="TAT-lag"])
  }else{
    mean.dat$TAT <- "multiple"
  }
  
  
  
  mean.dat$test_rotation <- test_rotation_name
  #mean.dat$prop_asym = prop.asym
  mean.dat$virus_par = unique(virus.par$version)
  mean.dat$distance_limit = unique(input.par$par1[input.par$parameter=="group-size-limit"])
  
  
  avg.prop.asym <- mean(c(unlist(out.prop.asym)))
  
  mean.dat$prop_asym= avg.prop.asym
  
  #and the long version
  mean.daily$type =  "all_exposures"
  mean.cumulative$type =   "cumulative"
  mean.cal$type = "UCB"
  mean.ala$type = "AlamedaCo"
  mean.symp$type = "symptomatic_iso"
  mean.trace$type = "tracing_iso"
  mean.test$type = "testing_iso"
  #mean.R$type = "Reffective"
  mean.iso$type= "isolations"
  #don't bother with employ-cat 00 can add later if needed
  mean.cat.long.list = convert.cat(mean.cat)
  mean.cat.long = data.table::rbindlist(mean.cat.long.list)
  #mean.cat.long = do.call("rbind", mean.cat.long.list)
  
  
  names(mean.daily) <- names(mean.cumulative) <- names(mean.cal) <- names(mean.ala) <- names(mean.symp) <- names(mean.trace) <- names(mean.test) <-  names(mean.iso)  <- c("mean", "lci", "uci", "type") #<- names(mean.R)
  
  mean.long <- rbind(mean.daily, mean.cumulative, mean.cal, mean.ala, mean.symp, mean.trace, mean.test, mean.iso, mean.cat.long)#, mean.R) 
  
  n.cat = length(input.pop)
  mean.long$day = c(rep(mean.time, (8+(3*n.cat))))#, mean.time[-1])
  
  
  mean.long$LOD <- LOD
  
  
  if(length(unique(input.par$par1[input.par$parameter=="TAT-lag"]))==1){
    mean.long$TAT <- unique(input.par$par1[input.par$parameter=="TAT-lag"])
  }else{
    mean.long$TAT <- "multiple"
  }
  
  mean.long$test_rotation <- test_rotation_name
  #mean.long$prop_asym = prop.asym
  mean.long$virus_par = unique(virus.par$version)
  
  mean.long$distance_limit = unique(input.par$par1[input.par$parameter=="group-size-limit"])
  
  
  mean.long$prop_asym = avg.prop.asym
  
  
  # mean.sum$sim_cat = sim_cat
  # #mean.sum$prop_asym = prop.asym
  # mean.sum$virus_par = unique(virus.par$version)
  # 
  # mean.sum$superspread = superspread
  # mean.sum$distance_limit = unique(input.par$par1[input.par$parameter=="group-size-limit"])
  
  # 
  # 
  # #and summarize R
  # mean.R = summarise.R(out.list.R=out.R, day.vec = mean.dat$day, n.reps=n.reps)
  # mean.R$LOD <- LOD
  # mean.R$TAT <- unique(input.par$par1[input.par$parameter=="TAT-lag"])
  # mean.R$test_rotation <- test_rotation_name
  # #mean.R$sim_cat = sim_cat
  # #mean.R$prop_asym = prop.asym
  # mean.R$virus_par = unique(virus.par$version)
  # 
  # mean.R$superspread = superspread
  # mean.R$distance_limit = unique(input.par$par1[input.par$parameter=="group-size-limit"])
  # 
  # mean.R$prop_asym <- avg.prop.asym
  # 
  # 
  # mean.R.mat = manage.R.matrix(mat.list=R.mat.out)
  
  
  #and do the best you can with the R-output
  #put it all together
  #R.mat.use <- do.call("rbind", R.mat.out)
  
  
  R.mat.use <- data.table::rbindlist(R.mat.out)
  R.mat.use <- arrange(R.mat.use, total_potential_cases)
  
  if(virus.par$distribution[virus.par$parameter=="R0"]=="negbinom"){
    mean.R.mat = R.fit.sum(R.mat.use)  
  }else{
    mean.R.mat = R.fit.sum.lognorm(R.mat.use) 
  }
  
  rownames(mean.R.mat) <- c()
  mean.R.mat$LOD <- LOD
  
  
  
  if(length(unique(input.par$par1[input.par$parameter=="TAT-lag"]))==1){
    mean.R.mat$TAT <- unique(input.par$par1[input.par$parameter=="TAT-lag"])
  }else{
    mean.R.mat$TAT <- "multiple"
  }
  
  mean.R.mat$test_rotation <- test_rotation_name
  mean.R.mat$virus_par = unique(virus.par$version)
  mean.R.mat$distance_limit = unique(input.par$par1[input.par$parameter=="group-size-limit"])
  mean.R.mat$prop_asym <- avg.prop.asym
  
  
  #return these summaries and the list of pop.mats
  return(list(mean.dat, mean.long, pop.mat.chain, mean.R.mat))  
}


pop.par.base$par1[pop.par.base$parameter=="TAT-lag"] <- 1
pop.par.base$par2[pop.par.base$parameter=="TAT-lag"] <- .5
pop.par.base$par1[pop.par.base$parameter=="test-rotation"] <- "two-week"
pop.par.base$par1[pop.par.base$parameter=="n-test-days-per-week"] <- 7
pop.par.base$par1[pop.par.base$parameter=="test-on"] <- TRUE
pop.par.base$par1[pop.par.base$parameter=="test-freq"] <- 14
pop.par.base$par1[pop.par.base$parameter=="prop.cases.UCB"] <- .5


out = replicate.epidemic(n.reps = 100,
                         virus.par = virus.par,
                         input.par = pop.par.base,
                         input.pop=c(20000),#2000
                         n.init.exposed.vector=c(100),#10
                         times=365*2,
                         bay.area.prev = .1/100,
                         initial.R = 2.5,
                         within.host.theta = .72,
                         burnin = 0,
                         length_timestep=1,
                         employ.id.vector = c(1),
                         LOD=(10^1),
                         titer.dat = titer.dat,
                         test_rotation_name = "two-week-7-test-days")

save(out, file = "two-week-7-test-days.Rdata")
