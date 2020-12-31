rm(list=ls())

library(plyr)
library(dplyr)


setwd("/Users/caraebrook/Documents/R/R_repositories/Berkeley-Reopening/Dec-2020/all-runs/FigS3-Contact-Trace-Enhancement/out/")

load("biweekly-test-TAT1-LOD10-trace.Rdata")
dat1 <- out[[1]]

load("biweekly-test-TAT2-LOD10-trace.Rdata")
dat2 <- out[[1]]

load("biweekly-test-TAT3-LOD10-trace.Rdata")
dat3 <- out[[1]]

load("biweekly-test-TAT4-LOD10-trace.Rdata")
dat4 <- out[[1]]

load("biweekly-test-TAT5-LOD10-trace.Rdata")
dat5 <- out[[1]]

load("biweekly-test-TAT10-LOD10-trace.Rdata")
dat6 <- out[[1]]

dat <- rbind(dat1,dat2, dat3, dat4, dat5, dat6)
rm(dat1,dat2, dat3, dat4, dat5, dat6)
unique(dat$LOD)
unique(dat$TAT)
unique(dat$test_rotation)
unique(dat$distance_limit)



load("weekly-test-TAT1-LOD10-trace.Rdata")
dat1 <- out[[1]]

load("weekly-test-TAT2-LOD10-trace.Rdata")
dat2 <- out[[1]]

load("weekly-test-TAT3-LOD10-trace.Rdata")
dat3 <- out[[1]]

load("weekly-test-TAT4-LOD10-trace.Rdata")
dat4 <- out[[1]]

load("weekly-test-TAT5-LOD10-trace.Rdata")
dat5 <- out[[1]]

load("weekly-test-TAT10-LOD10-trace.Rdata")
dat6 <- out[[1]]

dat <- rbind(dat, dat1,dat2, dat3, dat4, dat5, dat6)
rm(dat1,dat2, dat3, dat4, dat5, dat6)
unique(dat$LOD)
unique(dat$TAT)
unique(dat$test_rotation)
unique(dat$distance_limit)


load("two-week-test-TAT1-LOD10-trace.Rdata")
dat1 <- out[[1]]

load("two-week-test-TAT2-LOD10-trace.Rdata")
dat2 <- out[[1]]

load("two-week-test-TAT3-LOD10-trace.Rdata")
dat3 <- out[[1]]

load("two-week-test-TAT4-LOD10-trace.Rdata")
dat4 <- out[[1]]

load("two-week-test-TAT5-LOD10-trace.Rdata")
dat5 <- out[[1]]

load("two-week-test-TAT10-LOD10-trace.Rdata")
dat6 <- out[[1]]

dat <- rbind(dat, dat1,dat2, dat3, dat4, dat5, dat6)
rm(dat1,dat2, dat3, dat4, dat5, dat6)
unique(dat$LOD)
unique(dat$TAT)
unique(dat$test_rotation)
unique(dat$distance_limit)


dat$iso_lag = Inf


#and load in the symptom isolations
load("symptom-iso-1-trace.Rdata")
dat1 <- out[[1]]
dat1$iso_lag = 1
load("symptom-iso-2-trace.Rdata")
dat2 <- out[[1]]
dat2$iso_lag = 2
load("symptom-iso-3-trace.Rdata")
dat3 <- out[[1]]
dat3$iso_lag=3
load("symptom-iso-4-trace.Rdata")
dat4 <- out[[1]]
dat4$iso_lag=4
load("symptom-iso-5-trace.Rdata")
dat5 <- out[[1]]
dat5$iso_lag = 5


dat <- rbind(dat, dat1,dat2, dat3, dat4, dat5)
rm(dat1,dat2, dat3, dat4, dat5)
unique(dat$LOD)
unique(dat$TAT)
unique(dat$test_rotation)
unique(dat$distance_limit)


dat$intervention_class <- "testing"
dat$intervention_class[dat$iso_lag<Inf] <- "symptom-isolation"

#and a 0 for reference
load("/Users/caraebrook/Documents/R/R_repositories/Berkeley-Reopening/Dec-2020/all-runs/Fig2-Group-Limits/out/no-group-lim.Rdata")
dat1 <- out[[1]]
dat1$iso_lag = Inf
dat1$intervention_class <- "none"


dat <- rbind(dat, dat1)

dat$trace_lag = 1

dat.trace <- dat

#now load the previous, add trace lag of Inf and compare
load("/Users/caraebrook/Documents/R/R_repositories/Berkeley-Reopening/Dec-2020/all-runs/Fig3-Testing-TAT-LOD-Freq/out/dat.all.int.12.28.Rdata")
dat = subset(dat, LOD!="1e+03" & LOD!="1e+05")

head(dat)
head(dat.trace)
dat$trace_lag = Inf

dat.trace = rbind(dat.trace, dat)

unique(dat.trace$LOD)
dat.trace = subset(dat.trace, LOD ==10)
dat.trace$LOD[dat.trace$LOD=="10"] <- "1e+01"
#and save 
save(dat.trace, file = "dat.trace.12.28.Rdata")


### and for R0
rm(list=ls())

setwd("/Users/caraebrook/Documents/R/R_repositories/Berkeley-Reopening/Dec-2020/all-runs/FigS3-Contact-Trace-Enhancement/out/")

load("biweekly-test-TAT1-LOD10-trace.Rdata")
dat1 <- out[[4]]

load("biweekly-test-TAT2-LOD10-trace.Rdata")
dat2 <- out[[4]]

load("biweekly-test-TAT3-LOD10-trace.Rdata")
dat3 <- out[[4]]

load("biweekly-test-TAT4-LOD10-trace.Rdata")
dat4 <- out[[4]]

load("biweekly-test-TAT5-LOD10-trace.Rdata")
dat5 <- out[[4]]

load("biweekly-test-TAT10-LOD10-trace.Rdata")
dat6 <- out[[4]]

dat <- rbind(dat1,dat2, dat3, dat4, dat5, dat6)
rm(dat1,dat2, dat3, dat4, dat5, dat6)
unique(dat$LOD)
unique(dat$TAT)
unique(dat$test_rotation)
unique(dat$distance_limit)

load("weekly-test-TAT1-LOD10-trace.Rdata")
dat1 <- out[[4]]

load("weekly-test-TAT2-LOD10-trace.Rdata")
dat2 <- out[[4]]

load("weekly-test-TAT3-LOD10-trace.Rdata")
dat3 <- out[[4]]

load("weekly-test-TAT4-LOD10-trace.Rdata")
dat4 <- out[[4]]

load("weekly-test-TAT5-LOD10-trace.Rdata")
dat5 <- out[[4]]

load("weekly-test-TAT10-LOD10-trace.Rdata")
dat6 <- out[[4]]

dat <- rbind(dat, dat1,dat2, dat3, dat4, dat5, dat6)
rm(dat1,dat2, dat3, dat4, dat5, dat6)
unique(dat$LOD)
unique(dat$TAT)
unique(dat$test_rotation)
unique(dat$distance_limit)

load("two-week-test-TAT1-LOD10-trace.Rdata")
dat1 <- out[[4]]

load("two-week-test-TAT2-LOD10-trace.Rdata")
dat2 <- out[[4]]

load("two-week-test-TAT3-LOD10-trace.Rdata")
dat3 <- out[[4]]

load("two-week-test-TAT4-LOD10-trace.Rdata")
dat4 <- out[[4]]

load("two-week-test-TAT5-LOD10-trace.Rdata")
dat5 <- out[[4]]

load("two-week-test-TAT10-LOD10-trace.Rdata")
dat6 <- out[[4]]

dat <- rbind(dat, dat1,dat2, dat3, dat4, dat5, dat6)
rm(dat1,dat2, dat3, dat4, dat5, dat6)
unique(dat$LOD)
unique(dat$TAT)
unique(dat$test_rotation)
unique(dat$distance_limit)


dat$iso_lag = Inf
dat$intervention_class <- "testing"


#and load in the symptom isolations
load("symptom-iso-1-trace.Rdata")
dat1 <- out[[4]]
dat1$iso_lag = 1
dat1$intervention_class = "symptom-isolation"
load("symptom-iso-2-trace.Rdata")
dat2 <- out[[4]]
dat2$iso_lag = 2
dat2$intervention_class = "symptom-isolation"
load("symptom-iso-3-trace.Rdata")
dat3 <- out[[4]]
dat3$iso_lag=3
dat3$intervention_class = "symptom-isolation"
load("symptom-iso-4-trace.Rdata")
dat4 <- out[[4]]
dat4$iso_lag=4
dat4$intervention_class = "symptom-isolation"
load("symptom-iso-5-trace.Rdata")
dat5 <- out[[4]]
dat5$iso_lag = 5
dat5$intervention_class = "symptom-isolation"


dat <- rbind(dat, dat1,dat2, dat3, dat4, dat5)
rm(dat1,dat2, dat3, dat4, dat5)
unique(dat$LOD)
unique(dat$TAT)
unique(dat$test_rotation)
unique(dat$distance_limit)


#and load the comparative at 0

#and a 0 for reference
load("/Users/caraebrook/Documents/R/R_repositories/Berkeley-Reopening/Dec-2020/all-runs/Fig2-Group-Limits/out/no-group-lim.Rdata")
dat1 <- out[[4]]
dat1$iso_lag = Inf
dat1$intervention_class <- "none"


dat <- rbind(dat, dat1)

dat$trace_lag = 1

dat.trace.R0 <- dat


dat.sum.list <- dlply(dat.trace.R0, .(intervention_class, test_rotation, iso_lag, trace_lag))

get.rel.R0 <- function(dat1){
  start.R = dat1$mean.all[dat1$class=="UCB_post_titer_potential_cases"]
  start.R.lci = dat1$lci.all[dat1$class=="UCB_post_titer_potential_cases"]
  start.R.uci = dat1$uci.all[dat1$class=="UCB_post_titer_potential_cases"]
  
  finish.R = dat1$mean.all[dat1$class=="UCB_post_isolations_actual_cases"]
  finish.R.lci = dat1$lci.all[dat1$class=="UCB_post_isolations_actual_cases"]
  finish.R.uci = dat1$uci.all[dat1$class=="UCB_post_isolations_actual_cases"]
  
  out.R1 <- cbind.data.frame(start.R,  start.R.lci,  start.R.uci)
  out.R2 <- cbind.data.frame(finish.R,  finish.R.lci,  finish.R.uci)
  names(out.R1) <- names(out.R2) <- c("mean", "lci", "uci")  
  out.R <- rbind(out.R1, out.R2)
  out.R$intervention <- c("pre", "post")
  
  reduct.R = out.R1-out.R2
  reduct.R$LOD  <- unique(dat1$LOD)
  reduct.R$TAT  <- unique(dat1$TAT)
  reduct.R$test_rotation  <- unique(dat1$test_rotation)
  reduct.R$distance_limit  <- unique(dat1$distance_limit)
  #reduct.R$n_test_days  <- unique(dat1$n_test_days)
  reduct.R$iso_lag  <- unique(dat1$iso_lag)
  reduct.R$intervention_class  <- unique(dat1$intervention_class)
  #reduct.R$trace_lag  <- unique(dat1$trace_lag)
  reduct.R$mean[reduct.R$mean<0] <-0
  reduct.R$lci[reduct.R$lci<0] <-0
  reduct.R$uci[reduct.R$uci<0] <-0
  
  return(reduct.R)
  
}

dat.out.list <- lapply(dat.sum.list, get.rel.R0)

dat.trace.R0  <- do.call("rbind", dat.out.list)
rownames(dat.trace.R0) <- c()

dat.trace.R0$LOD[dat.trace.R0$LOD==10] <- "1e+01"
dat.trace.R0$LOD[dat.trace.R0$LOD==1000] <- "1e+03"

dat.trace.R0$trace_lag = 1

load("/Users/caraebrook/Documents/R/R_repositories/Berkeley-Reopening/Dec-2020/all-runs/Fig3-Testing-TAT-LOD-Freq/out/R0.all.int.12.28.Rdata")
dat.test.R0 = subset(dat.test.R0, LOD!="1e+03" & LOD!="1e+05")
head(dat.test.R0)
head(dat.trace.R0)

dat.test.R0$trace_lag <- Inf

head(dat.test.R0)
head(dat.trace.R0)
dat.trace.R0 <- rbind(dat.test.R0, dat.trace.R0)

unique(dat.trace.R0$LOD)
unique(dat.trace.R0$TAT)
unique(dat.trace.R0$test_rotation)
unique(dat.trace.R0$distance_limit)
unique(dat.trace.R0$iso_lag)
unique(dat.trace.R0$trace_lag)
dat.trace.R0$iso_lag = factor(dat.trace.R0$iso_lag, levels = c(1,2,3,4,5,Inf))
dat.trace.R0$TAT = factor(dat.trace.R0$TAT, levels = c(1,2,3,4,5,10, Inf))

#and save and plot.
setwd("/Users/caraebrook/Documents/R/R_repositories/Berkeley-Reopening/Dec-2020/all-runs/FigS3-Contact-Trace-Enhancement/out/")
save(dat.trace.R0, file = "R0.trace.12.28.Rdata")
