rm(list=ls())

#setwd("/Users/caraebrook/Documents/R/R_repositories/Berkeley-Reopening/final-model-runs/Fig3-Symptom-Iso-Lags")

#no group, no test, no trace

library(data.table)
library(plyr)
library(dplyr)
library(EpiEstim)
library(deSolve)
library(matrixStats)
library(fitdistrplus)

#load parameters including pre-run titer trajectories for each individual
load("titer.dat.20K.Rdata")
load("virus.par.9.4.Rdata")
load("pop.par.base.Rdata")


grab.titer <- function( gen.time, dat.vir){
  
  titer.inf <- dat.vir$V[dat.vir$time>gen.time$genTime][1]
  
  #we'll say that we need ~a dose of 1 million virions to cause an infection, so the FOI is
  FOI = titer.inf/1000000
  #and the probability of infection happening is
  ProbExposure = 1-(exp(-1*FOI))
  ProbExposure[ProbExposure<0] <- 0
  #does the infection happen? make it a probabilistic outcome of the titer
  #then, you role a dice to see if this exposure causes an infection
  InfectionYN = sample(c(0,1), size=1, prob = c(1-ProbExposure, ProbExposure))
  
  if(InfectionYN==1){
    return(gen.time)
  }else{
    gen.time$genTime <- NA
    return(gen.time)
  }}
simple.target <- function(t,x,parms){
  
  Tc = x[1]
  E = x[2]
  I = x[3]
  V = x[4]
  
  
  # The with() function gives access to the named values of parms within the
  # local environment created by the function
  with(as.list(parms),{
    
    dTcdt <- -parms$beta*Tc*V
    dEdt <- parms$beta*Tc*V - parms$k*E
    dIdt <- parms$k*E - parms$sigma*I
    dVdt <- parms$p*I - parms$c*V
    
    list(c(dTcdt,dEdt, dIdt,dVdt))
  })
}
assess.titer <- function(gen.dat2, dat.vir){
  #first, if the generation time is already NA, no need to assess titer!
  if(is.na(gen.dat2$generation_time)){
    return(gen.dat2)
  }else{
    #get the titer at the infection event
    titer.inf <- dat.vir$V[dat.vir$time>gen.dat2$generation_time][1]
    #we'll say that we need ~a dose of 1 million virions to cause an infection, so the FOI is
    FOI = titer.inf/1000000
    #and the probability of infection happening is
    ProbExposure = 1-(exp(-1*FOI))
    ProbExposure[ProbExposure<0] <- 0
    #does the infection happen? make it a probabilistic outcome of the titer
    #then, you role a dice to see if this exposure causes an infection
    InfectionYN = sample(c(0,1), size=1, prob = c(1-ProbExposure, ProbExposure))
    if (InfectionYN==1){
      return(gen.dat2)
    }else if (InfectionYN==0){
      gen.dat2$generation_time = NA
      return(gen.dat2)
    }
  }
}
normal_fn <- function(meanpar=NULL, sdpar=NULL){
  out <- purrr::partial(rnorm,
                        mean = meanpar,
                        sd = sdpar)
  return(out)
} 
lognormal_fn <- function(meanlogpar=NULL, sdlogpar=NULL){
  out <- purrr::partial(rlnorm,
                        meanlog = meanlogpar,
                        sdlog = sdlogpar)
  return(out)
} 
one.model.run <- function(beta, sigma, k, c, p, xstart, times){
  
  #build parameter matrix
  params = list(beta=beta,
                sigma=sigma,
                c=c,
                k=k,
                p=p)
  
  out <- as.data.frame(lsoda(y = xstart, times = times, func = simple.target, parms = params))
  out.rep <- dplyr::select(out, time, V)
  return(out.rep)
}
wrap.viral.loads <- function(beta.mean, sigma.mean, k.mean, c.mean, p.mean, p.sd, n.pop, times, xstart){
  
  #first set up draw distributions
  draw_beta = normal_fn(meanpar = exp(beta.mean), sdpar = beta.mean/3)
  draw_sigma = normal_fn(meanpar = exp(sigma.mean), sdpar = sigma.mean/3)
  draw_k = normal_fn(meanpar = exp(k.mean), sdpar = k.mean/3)
  draw_c = normal_fn(meanpar = exp(c.mean), sdpar = c.mean/3)
  draw_p = lognormal_fn(meanlogpar = log(p.mean), sdlogpar = log(p.sd))
  
  #p is drawn from your range
  
  #draw parameters from the lognormal distribution
  beta.list = as.list(log(draw_beta(n.pop)))
  sigma.list = as.list(log(draw_sigma(n.pop)))
  c.list = as.list(log(draw_c(n.pop)))
  k.list = as.list(log(draw_k(n.pop)))
  p.list = as.list(draw_p(n.pop))
  
  #and run over all the individuals in the pop
  out.list = mapply(one.model.run, beta=beta.list, k = k.list,
                    sigma=sigma.list, c=c.list, p=p.list, 
                    MoreArgs = list(xstart=xstart, times=times), SIMPLIFY = FALSE) 
  
  #each entry in list givens another individual
  employ.ID.list <- as.list(1:n.pop)
  
  out.list <- mapply(add.seq.ID, dat = out.list, employID = employ.ID.list, SIMPLIFY = FALSE)
  
  #and bind
  out.df <- do.call("rbind", out.list)
  head(out.df)
  
  #and attach par 
  out.df$beta <- rep(c(unlist(beta.list)), each = length(times))
  out.df$sigma <- rep(c(unlist(sigma.list)), each = length(times))
  out.df$c <- rep(c(unlist(c.list)), each = length(times))
  out.df$k <- rep(c(unlist(k.list)), each = length(times))
  out.df$p <- rep(c(unlist(p.list)), each = length(times))
  
  
  #and return
  return(out.df)
}
add.seq.ID <- function(dat, employID){
  dat$employIDs <- employID
  return(dat)
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
  
  
  out.vect  <- pop_dat$add 
  
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
  
  
  out.vect  <- pop_dat$add 
  
  
  #then add deaths based on each pop cat
  pop.cat = unique(dat$employ_cat)
  
  out.vect2 <- as.numeric(input_par$par1[input_par$parameter=="CFR" & input_par$population==pop.cat])*out.vect
  
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
add.noise<- function(par){
  par.new <-  rnorm(n=1, mean=par, sd=.25*par)    
  #not allowed to be negative, so if, so, give a small number:
  if(par.new<0){
    par.new = .0000001
  }
  return(par.new)
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
    new.mat <- do.call("rbind", new.mat.list)
    
    
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
    new.mat = do.call("rbind", new.list.out)
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
    
    new.mat = do.call("rbind", dat.new.split.out)
    
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
    new.dat <- do.call("rbind", new.dat.list)
    
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
cull.titer <- function(dat.vir, gen.dat){
  
  #in case there are multiple infections, split by row
  gen.dat.split <- dlply(gen.dat, .(rownames(gen.dat)))
  
  gen.dat.split2 <-  lapply(X = gen.dat.split,FUN= assess.titer, dat.vir=dat.vir)
  #and apply titer assessment
  
  #and recombine and return
  gen.dat.new = do.call("rbind", gen.dat.split2)
  
  return(gen.dat.new)
  
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
grab.next = function(dat){
  out.par = dat$generation_time
  return(out.par)
}
get.gen = function(dat, genTime, dat.vir){
  #rep.name = rep(dat$employ_ids, dat$potential_cases_caused)
  if(dat$potential_cases_caused>0){
    n.cases <- dat$potential_cases_caused
    #rep.gen = sort(genTime(dat$potential_cases_caused))  
    
    #first, how many cases would there have been??
    if(n.cases>0){
      rep.gen.dat = cbind.data.frame(genTime(n.cases))
    }else{
      rep.gen.dat = cbind.data.frame(NA)
    }
    
    names(rep.gen.dat) <-  "genTime"
    
    #and make list:
    rep.gen.list = dlply(rep.gen.dat, .(rownames(rep.gen.dat)))
    
    #now, lose some to the titer dynamics:
    rep.gen <- do.call("rbind",lapply(rep.gen.list, grab.titer, dat.vir=dat.vir))
    rep.gen$employ_ids <- dat$employ_ids
    rep.gen <- rep.gen[complete.cases(rep.gen),]
    
    #rep.gen is your cases post-titer dynamics
    
    if(nrow(rep.gen)>0){
      rep.gen <- arrange(rep.gen, genTime)
      rownames(rep.gen) <- c()
      rep.gen.int <- rep.gen
      
    }else{
      rep.gen <- cbind.data.frame(dat$employ_ids, NA)  
      names(rep.gen) <- c("employ_ids", "genTime")
      rep.gen.int <- rep.gen
      
    }
    
  }else{
    rep.name = dat$employ_ids
    rep.gen = NA
    dat.out = cbind.data.frame(rep.name, rep.gen)
    names(dat.out) = c("employ_ids", "genTime")
    rep.gen = rep.gen.int = dat.out
  }
  
  #for all, reorder
  rep.gen <- dplyr::select(rep.gen, employ_ids, genTime)
  rep.gen.int <- dplyr::select(rep.gen.int, employ_ids, genTime)
  return(list(rep.gen, rep.gen.int))
}
get.gen.super = function(dat, genTime, thresh, dist.limit, perc.inf.super, dat.vir){
  #we assume that all of the cases above some threshold 
  #case number are generated in a single superspreading event
  dist.limit.obs = dat$obs_distancing_limit
  #print(dat$employ_ids)
  #rep.name = rep(dat$employ_ids, dat$potential_cases_caused)
  
  if(dat$potential_cases_caused>0){
    n.case.super = dat$potential_cases_caused-thresh # all of these happen together
    n.case.super[n.case.super<0] <- 0
    n.case.disp =  dat$potential_cases_caused-n.case.super
    
    
    #first, how many cases would there have been??
    if(n.case.disp>0){
      rep.gen1 = cbind.data.frame("disp", genTime(n.case.disp))  
    }else(
      rep.gen1 = cbind.data.frame("disp", NA)  
    )
    if(n.case.super>0){
      rep.gen2 =  cbind.data.frame("super", rep(genTime(1), n.case.super))#here is where you would add in multiple superspreading events if you wanted  
    }else{
      rep.gen2 = cbind.data.frame("super", NA)  
    }
    
    names(rep.gen1) <- names(rep.gen2) <- c("type", "genTime")
    rep.gen.dat <- arrange(rbind(rep.gen1, rep.gen2), genTime)
    rep.gen.dat <- rep.gen.dat[complete.cases(rep.gen.dat),]
    
    rep.gen.list = dlply(rep.gen.dat, .(rownames(rep.gen.dat)))
    
    
    
    #now, lose some to the titer dynamics:
    rep.gen <- do.call("rbind",lapply(rep.gen.list, grab.titer, dat.vir=dat.vir))
    rep.gen$employ_ids <- dat$employ_ids
    rep.gen <- rep.gen[complete.cases(rep.gen),]
    
    if(nrow(rep.gen)>0){
      rep.gen <- arrange(rep.gen, genTime)
      rep.gen <- dplyr::select(rep.gen, employ_ids, genTime, type)
      rownames(rep.gen) <- c()
      rep.gen.int <- rep.gen
      
    }else{
      rep.gen <- cbind.data.frame(dat$employ_ids, NA,NA)  
      names(rep.gen) <- c("employ_ids", "genTime", "type")
      rep.gen.int <- rep.gen
      
    }
    
    
    
    #then, if there is an association limit that is abided
    #delete those superspreading events that happen beyond the threshold
    
    if(dist.limit.obs==TRUE & n.case.super>(perc.inf.super*dist.limit) & n.case.disp>0){ 
      #superspreading event avoided but other cases generated at discrete intervals
      if(!is.na(rep.gen.int$genTime[1])){
        rep.gen.int = subset(rep.gen.int, type!="super")  
      }else{
        rep.gen.int = rep.gen.int
      }
      
      
      
      
    }else if(dist.limit.obs==TRUE &  n.case.super>(perc.inf.super*dist.limit) & n.case.disp==0){ 
      #superspreading event avoided and no cases generated at all
      rep.gen.int <- rep.gen[1,]
      rep.gen.int$type <- NA
      rep.gen.int$genTime <- NA
      
      
    }
    #otherwise, there could be an intervention but there are no superspread cases or
    #the superspreading event is below the intervention threshold so it still occurs, so you return above
    
  }else{ #no cases anticipated anyway
    
    rep.gen <- cbind.data.frame(dat$employ_ids, NA,NA)  
    names(rep.gen) <- c("employ_ids", "genTime", "type")
    rep.gen.int <- rep.gen
  }
  
  
  return(list(rep.gen, rep.gen.int))
}
scale.back = function(dat){
  if(nrow(dat)>1){
    dat.new = dat[2:nrow(dat),]
  }
  else{
    dat.new = rbind.data.frame(c(NA,NA))
  }
  names(dat.new) <- c("employ_ids", "generation_time")
  return(dat.new)
  
}
grab.next.scale.back = function(dat){
  out.par = dat$generation_time[1]
  if(nrow(dat)>1){
    dat.new = dat[2:nrow(dat),]
  }
  else{
    dat.new = rbind.data.frame(c(NA,NA))
  }
  names(dat.new) <- c("employ_ids", "generation_time")
  return(list(out.par, dat.new))
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
get.midpoint <- function(par.low, par.hi){
  diff.add = (par.hi-par.low)/2
  mid=par.low+diff.add
  return(mid)
}
initiate.pop <- function(start.ID.employ, pop.UCB, n.init.exposed, pop.ID, input.par,  R0fn, genTime, superspread, titer.dat, LOD){
  
  pop.par = subset(input.par, population==pop.ID)
  
  #make table one
  pop.mat = cbind.data.frame(matrix(NA, nrow=pop.UCB, ncol =27))
  names(pop.mat) <- c( "employ_ids","employ_cat", "state",  "traced",  "testing", "superspread", "obs_dist_limits",  "exposure_time",   "total_potential_cases_caused",  "original_potential_cases_caused_UCB", "post_titer_potential_cases_caused_UCB", "potential_cases_caused", "actual_cases_caused", "case_source", "infector", "time_test_sensitive_start", "time_test_sensitive_end", "infector_iso_time", "time_of_tracing_iso", "time_of_next_test", "time_of_testing_iso", "titer_lim_for_symptoms", "time_of_symptom_onset", "time_of_symptom_iso", "time_isolation", "reason_isolated",  "timestep")
  
  #and fill in all you can
  pop.mat$testing = pop.par$par1[pop.par$parameter=="test-on"]
  pop.mat$superspread = superspread
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
  
  #initially exposed get distributed at random
  employ.id.init = sample(pop.mat$employ_ids, size=n.init.exposed, replace=FALSE)
  
  
  for(i in 1:n.init.exposed){
    pop.mat$state[pop.mat$employ_ids==employ.id.init[i]] <- 1
  }
  
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
  
  
  #you should have already brought in a titer for everyone in your population 
  #choose a threshold titer for symptom onset 
  pop.mat$titer_lim_for_symptoms = titer_lim(pop.UCB)
  
  pop.split <- dlply(pop.mat, .(employ_ids))
  titer.split <- dlply(titer.dat, .(employ_ids))
  
  #now, based on this, go into each person's virus trajectory and calculate the timing of symptom onset
  #while you are at it, you can also look at their titer and the LOD and calculate the start/end times for which they are test sensitive
  
  pop.split.new <-mapply(get.symptom.onset, dat = pop.split, dat.vir=titer.split, MoreArgs = list(LOD=LOD), SIMPLIFY = FALSE)
  #when there is nothing that is under the limit, these infections become "asymptomatic" -- 
  #we can later play with the proportion that classify as this by modulating the mean value for the symptom onset limit
  
  pop.mat <- do.call("rbind", pop.split.new)
  
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
  
  #now, make list of generation times for all future infections
  #once you have the generation time, go into each person's virus trajectory and extract the titer load and calculate,
  #based on FOI, whether this infection actually happens. If not, remove it here.
  
  
  #then, keep the trajectory so that you can later evaluate how it measures against the test sensitivity
  
  #for time to next case, you need a list of all of the cases
  tmp.dat = cbind.data.frame(pop.mat$employ_ids, pop.mat$original_potential_cases_caused_UCB, pop.mat$obs_dist_limits)
  
  names(tmp.dat) = c("employ_ids", "potential_cases_caused", "obs_distancing_limit")
  
  
  #any of these that cause > superspread threshold will get the same generation time
  tmp.list = dlply(tmp.dat, .(employ_ids))
  
  
  
  #depending on whether superspreading dynamics are allowed, 
  #the generations times are either clustered together or distributed
  #the superspread function also eliminates some cases if they are above the group association limit
  #otherwise the case number does not change at all
  if(superspread==FALSE){
    gen_time_list <- mapply(get.gen, dat=tmp.list, dat.vir = titer.split, MoreArgs = list(genTime=genTime), SIMPLIFY = FALSE)  
    pre_int_gen_time_list <- sapply(gen_time_list, "[",1)
    post_int_gen_time_list <- sapply(gen_time_list, "[",2)
  }else{
    gen_time_list <- mapply(get.gen.super,dat=tmp.list, dat.vir = titer.split, MoreArgs = list(genTime=genTime, thresh=as.numeric(pop.par$par1[pop.par$parameter=="superspread-threshold"]), dist.limit =  as.numeric(pop.par$par1[pop.par$parameter=="group-size-limit"]), perc.inf.super = as.numeric(pop.par$par1[pop.par$parameter=="percent-infected-super"])), SIMPLIFY = FALSE)
    pre_int_gen_time_list <- sapply(gen_time_list, "[",1)
    post_int_gen_time_list <- sapply(gen_time_list, "[",2)
  }
  
  
  #this gives potential cases both before and afer group size limits - record
  #this whole section is just to record the data point
  dat.gen.titer = do.call("rbind", pre_int_gen_time_list)
  dat.gen.group = do.call("rbind", post_int_gen_time_list)
  
  dat.gen.titer <- dat.gen.titer[complete.cases(dat.gen.titer),]
  dat.gen = dat.gen.group
  names(dat.gen)[2] <- "generation_time"
  dat.gen <- dat.gen[,1:2]
  dat.gen.group <- dat.gen.group[complete.cases(dat.gen.group),]
  
  titer.case = ddply(dat.gen.titer, .(employ_ids), summarize, tot_cases = length(genTime))
  group.case = ddply(dat.gen.group, .(employ_ids), summarize, tot_cases = length(genTime))
  
  
  
  max_id = max(pop.mat$employ_ids)
  missing_ids_titer <- (start.ID.employ:max_id)[!(start.ID.employ:max_id %in% titer.case$employ_ids)]
  missing_ids_group <- (start.ID.employ:max_id)[!(start.ID.employ:max_id %in% group.case$employ_ids)]
  
  titer.case = data.table(titer.case)
  group.case = data.table(group.case)
  
  # add in missing days if any are missing
  if (length(missing_ids_titer > 0)) {
    titer.case <- data.table::rbindlist(list(titer.case,
                                             data.table(employ_ids = missing_ids_titer,
                                                        tot_cases = 0)))
  }
  
  if (length(missing_ids_group > 0)) {
    group.case <- data.table::rbindlist(list(group.case,
                                             data.table(employ_ids = missing_ids_group,
                                                        tot_cases = 0)))
  }
  
  
  titer.case <- arrange(titer.case, employ_ids)
  group.case <- arrange(group.case, employ_ids)
  
  
  
  #and cases post titer cull
  pop.mat$post_titer_potential_cases_caused_UCB <- titer.case$tot_cases
  
  #and here are those infections that actually occur
  pop.mat$potential_cases_caused <- group.case$tot_cases
  
  pop.mat$exposure_time[pop.mat$state>0] <- 0
  
  
  #first, assume that isolation time is symptomatic
  pop.mat$time_isolation[pop.mat$state==1 ] <- as.numeric(pop.mat$time_of_symptom_iso[pop.mat$state==1 ])
  pop.mat$time_isolation = as.numeric(pop.mat$time_isolation)
  pop.mat$reason_isolated[pop.mat$state==1 ] <- "symptom_iso"
  
  #now, if testing (and, for other cases, tracing) comes first, we replace it
  #test needs to be AFTER start time of test sensitive and before end time of test sensitive
  pop.mat$reason_isolated[pop.mat$state==1 & pop.mat$time_isolation> pop.mat$time_of_testing_iso & pop.mat$testing==TRUE & pop.mat$time_of_next_test>=pop.mat$time_test_sensitive_start & pop.mat$time_of_next_test<pop.mat$time_test_sensitive_end] <- "testing_iso"
  pop.mat$time_isolation[pop.mat$state==1 & pop.mat$time_isolation> pop.mat$time_of_testing_iso & pop.mat$testing==TRUE & pop.mat$time_of_next_test>=pop.mat$time_test_sensitive_start & pop.mat$time_of_next_test<pop.mat$time_test_sensitive_end] <- pop.mat$time_of_testing_iso[pop.mat$state==1 & pop.mat$time_isolation > pop.mat$time_of_testing_iso & pop.mat$testing==TRUE & pop.mat$time_of_next_test>=pop.mat$time_test_sensitive_start & pop.mat$time_of_next_test<pop.mat$time_test_sensitive_end]
  
  
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
update.test <- function(dat, testing.seq){
  if(!is.na(dat$exposure_time)){
    #take the lowest testing seq that follows the latest onset
    if(dat$next_test < dat$exposure_time){
      dat$next_test <- min(testing.seq[testing.seq>dat$exposure_time])
    }#otherwise nothing changes
  }
  
  
  #return data
  return(dat)
  
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
  dat.gen.new = do.call("rbind", gen_list)
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
  mat.2 <- do.call("cbind", mat)
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
manage.R.matrix <- function(mat.list){
  
  
  #first, trim to same length
  min.length <- min(unlist(lapply(mat.list, nrow)))
  
  
  for (i in 1:length(mat.list)){
    mat.list[[i]] <- mat.list[[i]][1:min.length,]
  }
  
  list.mean <- Reduce("+",mat.list) / length(mat.list)
  mat.2 <- do.call("cbind", mat.list)
  list.sd <- apply(mat.2, 1, sd, na.rm=T)
  
  list.uci = list.mean + 1.96*list.sd
  list.lci = list.mean - 1.96*list.sd
  
  list.lci[list.lci<0] <- 0
  list.uci[list.uci<0] <- 0
  
  dat= cbind.data.frame(list.mean, list.lci, list.uci)
  
  names(dat) = c("mean_total_potential_cases", "mean_UCB_potential_cases", "mean_UCB_post_group_potential_cases", "mean_UCB_post_titer_potential_cases", "mean_UCB_post_isolations_actual_cases",
                 "lci_total_potential_cases", "lci_UCB_potential_cases", "lci_UCB_post_group_potential_cases", "lci_UCB_post_titer_potential_cases", "lci_UCB_post_isolations_actual_cases",
                 "uci_total_potential_cases", "uci_UCB_potential_cases", "uci_UCB_post_group_potential_cases", "uci_UCB_post_titer_potential_cases", "uci_UCB_post_isolations_actual_cases")
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
  
  iso.db <- do.call("rbind", list.iso)
  iso.db$type = rep(1:n.cat, each = max.times)
  iso.db$type <- paste0("iso-pop-", iso.db$type)
  
  names(iso.db) <- c("mean", "lci", "uci", "type")
  
  exp.db <- do.call("rbind", list.exp)
  exp.db$type = rep(1:n.cat, each = max.times)
  exp.db$type <- paste0("exp-pop-", exp.db$type)
  
  names(exp.db) <- c("mean", "lci", "uci", "type")
  
  death.db <- do.call("rbind", list.deaths)
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
get.mean.sd.summary <- function(vector){
  
  name.vec = names(vector[[1]])
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
  names(dat) = c("mean", "lci", "uci")
  dat$parameter = name.vec
  return(dat)
} 
summarise.R = function(out.list.R, day.vec, n.reps){
  
  #first, trim to same length
  min.length <- min(unlist(lapply(out.list.R, length)))
  
  for (i in 1:length(out.list.R)){
    out.list.R[[i]] <- out.list.R[[i]][1:min.length]
  }
  #in case it is longer than R-effective
  day.vec = day.vec[1:min.length]
  
  dat.R = c(unlist(out.list.R))
  time = rep(day.vec, n.reps)
  trial = rep(1:n.reps, each=length(day.vec))
  
  dat.R = cbind.data.frame(time, dat.R, trial)
  names(dat.R) = c("day", "Reffective", "trial")
  dat.R$lci <- dat.R$uci <- NA
  
  mean.dat = data.frame(do.call("cbind", out.list.R))
  mean.dat$means = rowSums(mean.dat)
  mean.dat$means = mean.dat$means/n.reps
  tmp.sd = mean.dat[,1:n.reps]
  tmp.sd = as.matrix(tmp.sd)
  out.sd =  rowSds(tmp.sd)
  mean.dat$uci = mean.dat$means + 1.96*out.sd 
  mean.dat$lci = mean.dat$means - 1.96*out.sd 
  mean.dat$time = day.vec
  mean.dat = dplyr::select(mean.dat, time, means, lci, uci)
  
  add.dat = dplyr::select(mean.dat, time, means, lci, uci)
  add.dat$trial = "mean"
  add.dat <- dplyr::select(add.dat, time, means, trial, lci, uci)
  names(add.dat) <- names(dat.R)
  
  dat.R = rbind(dat.R, add.dat)
  return(dat.R)
  
}
simulate.epidemic <- function(input.pop, n.init.exposed.vector, employ.id.vector, times, virus.par, input.par, burnin, serial_mean, serial_sd, test.freq, length_timestep, bay.area.prev, bay.area.R, superspread, titer.dat, LOD, test_rotation_name){
  
  
  #first, draw the virus par that do not change depending on cohort of people:
  
  #sample incubation time to symptoms - you will be quarantined in the next timestep (could explore delay introductions here)
  #incfn = inc_fn(n_inc_samp=NULL,
  #              meanInc= virus.par$par1[virus.par$parameter=="incubation_time"], 
  #             sdInc=virus.par$par2[virus.par$parameter=="incubation_time"]) 
  
  
  #sample serial interval
  genTime = generationTime_fn(serial_dist = virus.par$distribution[virus.par$parameter=="generation_time"],
                              serial_shape= virus.par$par1[virus.par$parameter=="generation_time"], 
                              serial_scale= virus.par$par2[virus.par$parameter=="generation_time"])
  
  #sample delay to test positivity
  #testPos = inc_fn(n_inc_samp=NULL,
  #                meanInc= virus.par$par1[virus.par$parameter=="test_positive_delay"], 
  #               sdInc=virus.par$par2[virus.par$parameter=="test_positive_delay"]) 
  
  
  if (virus.par$distribution[virus.par$parameter=="R0"]=="log-normal"){
    
    #sample R0 normal
    R0fn = R0_fn(meanR0=virus.par$par1[virus.par$parameter=="R0"],
                 sdR0=virus.par$par2[virus.par$parameter=="R0"])
    
    
  }else if(virus.par$distribution[virus.par$parameter=="R0"]=="negbinom"){
    
    #sample R0 normal
    R0fn = R0_fn_nb(muR0=virus.par$par1[virus.par$parameter=="R0"],
                    sizeR0=virus.par$par2[virus.par$parameter=="R0"])
    
  }
  
  
  #and normal distribution of the detection limit
  
  #then, form your new populations
  #now split the population based on risk
  tot.pop = length(input.pop)
  pop.num = 1:tot.pop
  
  #make the proper number of pop.mat depending on the total number of subpopulations
  #populate each using the appropriate parameters
  out.list = mapply(FUN=initiate.pop, start.ID.employ = as.list(employ.id.vector), pop.UCB=as.list(input.pop), n.init.exposed= as.list(n.init.exposed.vector),  pop.ID = as.list(pop.num), 
                    MoreArgs= list(input.par=input.par, titer.dat=titer.dat,  R0fn=R0fn,  genTime=genTime, superspread=superspread, LOD=LOD))
  
  
  pop.list = out.list[1,]
  gen_list_long <- out.list[2,]
  #original.r0 <- out.list[3,][[1]]
  #gen_list_long_wkend <- out.list[3,]
  
  pop.mat <- do.call("rbind", pop.list)
  gen.dat.all <- do.call("rbind", gen_list_long)
  #gen.dat.all.wk <- do.call("rbind", gen_list_long_wkend)
  
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
  
  foi.bay.area = bay.area.R*bay.area.prev*length_timestep #rate per day at which susceptibles become infected
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
    pop.mat = do.call("rbind", pop.mat.list)#print(i)
    
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
  
  tmp = as.data.frame(do.call("cbind", new_col))
  names(tmp) <- paste0("isolations-employ-cat-", unique(input.par$population))
  
  tmp2 = as.data.frame(do.call("cbind", new_col_exp))
  names(tmp2) <- paste0("exposures-employ-cat-", unique(input.par$population))
  
  tmp3 = as.data.frame(do.call("cbind", new_col_deaths))
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
  pop.mat$TAT <- unique(input.par$par1[input.par$parameter=="TAT-lag"])
  pop.mat$test_rotation <- test_rotation_name
  
  
  return(list(daily_cases,pop.mat, prop.asym, R.mat))
}
replicate.epidemic = function(n.reps, input.pop, n.init.exposed.vector, employ.id.vector, times, virus.par, input.par, burnin, serial_mean, serial_sd, test.freq, length_timestep,
                              bay.area.prev, bay.area.R, test_rotation_name, superspread, LOD, titer.dat){
  
  out = replicate(n.reps, simulate.epidemic(virus.par = virus.par,
                                            input.par = input.par,
                                            input.pop=input.pop, 
                                            n.init.exposed.vector=n.init.exposed.vector, 
                                            times=times, 
                                            bay.area.prev = bay.area.prev, 
                                            bay.area.R = bay.area.R, 
                                            burnin = burnin,  
                                            serial_mean = serial_mean, 
                                            serial_sd = serial_sd, 
                                            length_timestep=length_timestep,
                                            employ.id.vector =employ.id.vector,
                                            superspread = superspread,
                                            LOD = LOD,
                                            titer.dat = titer.dat,
                                            test_rotation_name = test_rotation_name), simplify = "array")
  #make list
  out.time<- out.daily <- out.cal <- out.iso <- out.cumulative <-  out.ala <- out.symp <- out.trace <- out.test  <- out.iso <-out.cum.iso <- pop.mat.chain <- out.prop.asym <- R.mat.out <-  list()
  
  #and make list of all the categories of sub-pop
  out.cat <- list()
  
  for (i in 1:ncol(out)){
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
  mean.dat$TAT <- unique(input.par$par1[input.par$parameter=="TAT-lag"])
  mean.dat$test_rotation <- test_rotation_name
  #mean.dat$prop_asym = prop.asym
  mean.dat$virus_par = unique(virus.par$version)
  mean.dat$superspread = superspread
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
  mean.cat.long = do.call("rbind", mean.cat.long.list)
  
  
  names(mean.daily) <- names(mean.cumulative) <- names(mean.cal) <- names(mean.ala) <- names(mean.symp) <- names(mean.trace) <- names(mean.test) <-  names(mean.iso)  <- c("mean", "lci", "uci", "type") #<- names(mean.R)
  
  mean.long <- rbind(mean.daily, mean.cumulative, mean.cal, mean.ala, mean.symp, mean.trace, mean.test, mean.iso, mean.cat.long)#, mean.R) 
  
  n.cat = length(input.pop)
  mean.long$day = c(rep(mean.time, (8+(3*n.cat))))#, mean.time[-1])
  
  
  mean.long$LOD <- LOD
  mean.long$TAT <- unique(input.par$par1[input.par$parameter=="TAT-lag"])
  mean.long$test_rotation <- test_rotation_name
  #mean.long$prop_asym = prop.asym
  mean.long$virus_par = unique(virus.par$version)
  
  mean.long$superspread = superspread
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
  R.mat.use <- do.call("rbind", R.mat.out)
  R.mat.use <- arrange(R.mat.use, total_potential_cases)
  
  mean.R.mat = R.fit.sum(R.mat.use)
  rownames(mean.R.mat) <- c()
  mean.R.mat$LOD <- LOD
  mean.R.mat$TAT <- unique(input.par$par1[input.par$parameter=="TAT-lag"])
  mean.R.mat$test_rotation <- test_rotation_name
  mean.R.mat$virus_par = unique(virus.par$version)
  mean.R.mat$superspread = superspread
  mean.R.mat$distance_limit = unique(input.par$par1[input.par$parameter=="group-size-limit"])
  mean.R.mat$prop_asym <- avg.prop.asym
  
  
  #return these summaries and the list of pop.mats
  return(list(mean.dat, mean.long, pop.mat.chain, mean.R.mat))  
}



pop.par.base$par1[pop.par.base$parameter=="TAT-lag"] <- 1
pop.par.base$par2[pop.par.base$parameter=="TAT-lag"] <- .5
pop.par.base$par1[pop.par.base$parameter=="test-rotation"] <- "biweekly"
pop.par.base$par1[pop.par.base$parameter=="n-test-days-per-week"] <- 2
pop.par.base$par1[pop.par.base$parameter=="test-on"] <- TRUE
pop.par.base$par1[pop.par.base$parameter=="test-freq"] <- 3



out = replicate.epidemic(n.reps = 100,
                         virus.par = virus.par,
                         input.par = pop.par.base,
                         input.pop=c(20000),#2000
                         n.init.exposed.vector=c(100),#10
                         times=365*2,
                         bay.area.prev = .1/100,
                         bay.area.R = 1.5,
                         burnin = 0,
                         serial_mean = 4.8,
                         serial_sd = 2.3,
                         length_timestep=1,
                         employ.id.vector = c(1),
                         superspread=TRUE,
                         LOD=(10^1),
                         titer.dat = titer.dat,
                         test_rotation_name = "biweekly-TAT1-LOD10")

save(out, file = "biweekly-test-TAT1-LOD10.Rdata")
