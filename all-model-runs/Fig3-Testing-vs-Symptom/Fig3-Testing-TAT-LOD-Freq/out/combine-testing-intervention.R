rm(list=ls())

library(plyr)
library(dplyr)

#load group testing data with LOD
setwd("/Users/caraebrook/Documents/R/R_repositories/Berkeley-Reopening/final-model-runs/Fig3-Testing-TAT-LOD-Freq/out/")

load("biweekly-test-TAT1-LOD10.Rdata")
dat1 <- out[[1]]

load("biweekly-test-TAT1-LOD103.Rdata")
dat2 <- out[[1]]

load("biweekly-test-TAT1-LOD105.Rdata")
dat3 <- out[[1]]

load("biweekly-test-TAT1-LOD107.Rdata")
dat4 <- out[[1]]

dat <- rbind(dat1,dat2, dat3, dat4)
rm(dat1,dat2, dat3, dat4)
unique(dat$LOD)
unique(dat$TAT)
unique(dat$test_rotation)
unique(dat$distance_limit)


load("biweekly-test-TAT2-LOD10.Rdata")
dat1 <- out[[1]]

load("biweekly-test-TAT2-LOD103.Rdata")
dat2 <- out[[1]]

load("biweekly-test-TAT2-LOD105.Rdata")
dat3 <- out[[1]]

load("biweekly-test-TAT2-LOD107.Rdata")
dat4 <- out[[1]]

dat <- rbind(dat, dat1,dat2, dat3, dat4)
rm(dat1,dat2, dat3, dat4)
unique(dat$LOD)
unique(dat$TAT)
unique(dat$test_rotation)
unique(dat$distance_limit)


load("biweekly-test-TAT3-LOD10.Rdata")
dat1 <- out[[1]]

load("biweekly-test-TAT3-LOD103.Rdata")
dat2 <- out[[1]]

load("biweekly-test-TAT3-LOD105.Rdata")
dat3 <- out[[1]]

load("biweekly-test-TAT3-LOD107.Rdata")
dat4 <- out[[1]]

dat <- rbind(dat, dat1,dat2, dat3, dat4)
rm(dat1,dat2, dat3, dat4)
unique(dat$LOD)
unique(dat$TAT)
unique(dat$test_rotation)
unique(dat$distance_limit)



load("biweekly-test-TAT4-LOD10.Rdata")
dat1 <- out[[1]]

load("biweekly-test-TAT4-LOD103.Rdata")
dat2 <- out[[1]]

load("biweekly-test-TAT4-LOD105.Rdata")
dat3 <- out[[1]]

load("biweekly-test-TAT4-LOD107.Rdata")
dat4 <- out[[1]]

dat <- rbind(dat, dat1,dat2, dat3, dat4)
rm(dat1,dat2, dat3, dat4)
unique(dat$LOD)
unique(dat$TAT)
unique(dat$test_rotation)
unique(dat$distance_limit)



load("biweekly-test-TAT5-LOD10.Rdata")
dat1 <- out[[1]]

load("biweekly-test-TAT5-LOD103.Rdata")
dat2 <- out[[1]]

load("biweekly-test-TAT5-LOD105.Rdata")
dat3 <- out[[1]]

load("biweekly-test-TAT5-LOD107.Rdata")
dat4 <- out[[1]]

dat <- rbind(dat, dat1,dat2, dat3, dat4)
rm(dat1,dat2, dat3, dat4)
unique(dat$LOD)
unique(dat$TAT)
unique(dat$test_rotation)
unique(dat$distance_limit)


load("biweekly-test-TAT10-LOD10.Rdata")
dat1 <- out[[1]]

load("biweekly-test-TAT10-LOD103.Rdata")
dat2 <- out[[1]]

load("biweekly-test-TAT10-LOD105.Rdata")
dat3 <- out[[1]]

load("biweekly-test-TAT10-LOD107.Rdata")
dat4 <- out[[1]]

dat <- rbind(dat, dat1,dat2, dat3, dat4)
rm(dat1,dat2, dat3, dat4)
unique(dat$LOD)
unique(dat$TAT)
unique(dat$test_rotation)
unique(dat$distance_limit)


#and now weekly
load("weekly-test-TAT1-LOD10.Rdata")
dat1 <- out[[1]]

load("weekly-test-TAT1-LOD103.Rdata")
dat2 <- out[[1]]

load("weekly-test-TAT1-LOD105.Rdata")
dat3 <- out[[1]]

load("weekly-test-TAT1-LOD107.Rdata")
dat4 <- out[[1]]

dat <- rbind(dat, dat1,dat2, dat3, dat4)
rm(dat1,dat2, dat3, dat4)
unique(dat$LOD)
unique(dat$TAT)
unique(dat$test_rotation)
unique(dat$distance_limit)


load("weekly-test-TAT2-LOD10.Rdata")
dat1 <- out[[1]]

load("weekly-test-TAT2-LOD103.Rdata")
dat2 <- out[[1]]

load("weekly-test-TAT2-LOD105.Rdata")
dat3 <- out[[1]]

load("weekly-test-TAT2-LOD107.Rdata")
dat4 <- out[[1]]

dat <- rbind(dat, dat1,dat2, dat3, dat4)
rm(dat1,dat2, dat3, dat4)
unique(dat$LOD)
unique(dat$TAT)
unique(dat$test_rotation)
unique(dat$distance_limit)


load("weekly-test-TAT3-LOD10.Rdata")
dat1 <- out[[1]]

load("weekly-test-TAT3-LOD103.Rdata")
dat2 <- out[[1]]

load("weekly-test-TAT3-LOD105.Rdata")
dat3 <- out[[1]]

load("weekly-test-TAT3-LOD107.Rdata")
dat4 <- out[[1]]

dat <- rbind(dat, dat2, dat1, dat3, dat4)#
rm(dat1, dat2, dat3, dat4)
unique(dat$LOD)
unique(dat$TAT)
unique(dat$test_rotation)
unique(dat$distance_limit)

load("weekly-test-TAT4-LOD10.Rdata")
dat1 <- out[[1]]

load("weekly-test-TAT4-LOD103.Rdata")
dat2 <- out[[1]]

load("weekly-test-TAT4-LOD105.Rdata")
dat3 <- out[[1]]

load("weekly-test-TAT4-LOD107.Rdata")
dat4 <- out[[1]]

dat <- rbind(dat, dat1,dat2, dat3, dat4)
rm(dat1,dat2, dat3, dat4)
unique(dat$LOD)
unique(dat$TAT)
unique(dat$test_rotation)
unique(dat$distance_limit)



load("weekly-test-TAT5-LOD10.Rdata")
dat1 <- out[[1]]

load("weekly-test-TAT5-LOD103.Rdata")
dat2 <- out[[1]]

load("weekly-test-TAT5-LOD105.Rdata")
dat3 <- out[[1]]

load("weekly-test-TAT5-LOD107.Rdata")
dat4 <- out[[1]]

dat <- rbind(dat, dat1,dat2, dat3, dat4)
rm(dat1,dat2, dat3, dat4)
unique(dat$LOD)
unique(dat$TAT)
unique(dat$test_rotation)
unique(dat$distance_limit)


load("weekly-test-TAT10-LOD10.Rdata")
dat1 <- out[[1]]

load("weekly-test-TAT10-LOD103.Rdata")
dat2 <- out[[1]]

load("weekly-test-TAT10-LOD105.Rdata")
dat3 <- out[[1]]

load("weekly-test-TAT10-LOD107.Rdata")
dat4 <- out[[1]]

dat <- rbind(dat, dat1,dat2, dat3, dat4)
rm(dat1,dat2, dat3, dat4)
unique(dat$LOD)
unique(dat$TAT)
unique(dat$test_rotation)
unique(dat$distance_limit)



#and two-week
load("two-week-test-TAT1-LOD10.Rdata")
dat1 <- out[[1]]

load("two-week-test-TAT1-LOD103.Rdata")
dat2 <- out[[1]]

load("two-week-test-TAT1-LOD105.Rdata")
dat3 <- out[[1]]

load("two-week-test-TAT1-LOD107.Rdata")
dat4 <- out[[1]]

dat <- rbind(dat, dat1,dat2, dat3, dat4)
rm(dat1,dat2, dat3, dat4)
unique(dat$LOD)
unique(dat$TAT)
unique(dat$test_rotation)
unique(dat$distance_limit)


load("two-week-test-TAT2-LOD10.Rdata")
dat1 <- out[[1]]

load("two-week-test-TAT2-LOD103.Rdata")
dat2 <- out[[1]]

load("two-week-test-TAT2-LOD105.Rdata")
dat3 <- out[[1]]

load("two-week-test-TAT2-LOD107.Rdata")
dat4 <- out[[1]]

dat <- rbind(dat, dat1,dat2, dat3, dat4)
rm(dat1,dat2, dat3, dat4)
unique(dat$LOD)
unique(dat$TAT)
unique(dat$test_rotation)
unique(dat$distance_limit)


load("two-week-test-TAT3-LOD10.Rdata")
dat1 <- out[[1]]

load("two-week-test-TAT3-LOD103.Rdata")
dat2 <- out[[1]]

load("two-week-test-TAT3-LOD105.Rdata")
dat3 <- out[[1]]

load("two-week-test-TAT3-LOD107.Rdata")
dat4 <- out[[1]]

dat <- rbind(dat, dat1,dat2, dat3, dat4)
rm(dat1,dat2, dat3, dat4)
unique(dat$LOD)
unique(dat$TAT)
unique(dat$test_rotation)
unique(dat$distance_limit)



load("two-week-test-TAT4-LOD10.Rdata")
dat1 <- out[[1]]

load("two-week-test-TAT4-LOD103.Rdata")
dat2 <- out[[1]]

load("two-week-test-TAT4-LOD105.Rdata")
dat3 <- out[[1]]

load("two-week-test-TAT4-LOD107.Rdata")
dat4 <- out[[1]]

dat <- rbind(dat, dat1,dat2, dat3, dat4)
rm(dat1,dat2, dat3, dat4)
unique(dat$LOD)
unique(dat$TAT)
unique(dat$test_rotation)
unique(dat$distance_limit)



load("two-week-test-TAT5-LOD10.Rdata")
dat1 <- out[[1]]

load("two-week-test-TAT5-LOD103.Rdata")
dat2 <- out[[1]]

load("two-week-test-TAT5-LOD105.Rdata")
dat3 <- out[[1]]

load("two-week-test-TAT5-LOD107.Rdata")
dat4 <- out[[1]]

dat <- rbind(dat, dat1,dat2, dat3, dat4)
rm(dat1,dat2, dat3, dat4)
unique(dat$LOD)
unique(dat$TAT)
unique(dat$test_rotation)
unique(dat$distance_limit)


load("two-week-test-TAT10-LOD10.Rdata")
dat1 <- out[[1]]

load("two-week-test-TAT10-LOD103.Rdata")
dat2 <- out[[1]]

load("two-week-test-TAT10-LOD105.Rdata")
dat3 <- out[[1]]

load("two-week-test-TAT10-LOD107.Rdata")
dat4 <- out[[1]]

dat <- rbind(dat, dat1,dat2, dat3, dat4)
rm(dat1,dat2, dat3, dat4)
unique(dat$LOD)
unique(dat$TAT)
unique(dat$test_rotation)
unique(dat$distance_limit)

#and add in the symptom lag
dat$iso_lag = Inf

#and symptom stuff
setwd("/Users/caraebrook/Documents/R/R_repositories/Berkeley-Reopening/final-model-runs/Fig3-Symptom-Iso-Lags/out/")
load("symptom-iso-1.Rdata")
dat1 <- out[[1]]
dat1$iso_lag = 1

load("symptom-iso-2.Rdata")
dat2 <- out[[1]]
dat2$iso_lag = 2

load("symptom-iso-3.Rdata")
dat3 <- out[[1]]
dat3$iso_lag = 3

load("symptom-iso-4.Rdata")
dat4 <- out[[1]]
dat4$iso_lag = 4

load("symptom-iso-5.Rdata")
dat5 <- out[[1]]
dat5$iso_lag = 5


dat <- rbind(dat, dat1,dat2, dat3, dat4, dat5)
rm(dat1,dat2, dat3, dat4, dat5)
unique(dat$LOD)
unique(dat$TAT)
unique(dat$test_rotation)
unique(dat$distance_limit)
unique(dat$iso_lag)

dat$intervention_class <- "testing"
dat$intervention_class[dat$iso_lag!=Inf] <- "symptom-isolation"

dat.hold = dat

#and add in the group size limits
#and the nothing intervention - not at high asym but should not matter
load("/Users/caraebrook/Documents/R/R_repositories/Berkeley-Reopening/final-model-runs/Fig2-Group-Limits/out/dat.group.lim.10.21.Rdata")
dat$iso_lag = Inf
dat$intervention_class <- "group-size-limit"

dat <- rbind(dat.hold, dat)
dat$intervention_class[dat$iso_lag==Inf & dat$distance_limit==Inf & dat$test_rotation=="none"] <- "none" 

setwd("/Users/caraebrook/Documents/R/R_repositories/Berkeley-Reopening/final-model-runs/Fig3-Testing-TAT-LOD-Freq/out/")
save(dat, file = "dat.all.int.10.21.Rdata")




### and for R0
rm(list=ls())


#load goup testing data with LOD
setwd("/Users/caraebrook/Documents/R/R_repositories/Berkeley-Reopening/final-model-runs/Fig3-Testing-TAT-LOD-Freq/out/")

load("biweekly-test-TAT1-LOD10.Rdata")
dat1 <- out[[4]]

load("biweekly-test-TAT1-LOD103.Rdata")
dat2 <- out[[4]]

load("biweekly-test-TAT1-LOD105.Rdata")
dat3 <- out[[4]]

load("biweekly-test-TAT1-LOD107.Rdata")
dat4 <- out[[4]]

dat <- rbind(dat1,dat2, dat3, dat4)
rm(dat1,dat2, dat3, dat4)
unique(dat$LOD)
unique(dat$TAT)
unique(dat$test_rotation)
unique(dat$distance_limit)


load("biweekly-test-TAT2-LOD10.Rdata")
dat1 <- out[[4]]

load("biweekly-test-TAT2-LOD103.Rdata")
dat2 <- out[[4]]

load("biweekly-test-TAT2-LOD105.Rdata")
dat3 <- out[[4]]

load("biweekly-test-TAT2-LOD107.Rdata")
dat4 <- out[[4]]

dat <- rbind(dat, dat1,dat2, dat3, dat4)
rm(dat1,dat2, dat3, dat4)
unique(dat$LOD)
unique(dat$TAT)
unique(dat$test_rotation)
unique(dat$distance_limit)


load("biweekly-test-TAT3-LOD10.Rdata")
dat1 <- out[[4]]

load("biweekly-test-TAT3-LOD103.Rdata")
dat2 <- out[[4]]

load("biweekly-test-TAT3-LOD105.Rdata")
dat3 <- out[[4]]

load("biweekly-test-TAT3-LOD107.Rdata")
dat4 <- out[[4]]

dat <- rbind(dat, dat1,dat2, dat3, dat4)
rm(dat1,dat2, dat3, dat4)
unique(dat$LOD)
unique(dat$TAT)
unique(dat$test_rotation)
unique(dat$distance_limit)



load("biweekly-test-TAT4-LOD10.Rdata")
dat1 <- out[[4]]

load("biweekly-test-TAT4-LOD103.Rdata")
dat2 <- out[[4]]

load("biweekly-test-TAT4-LOD105.Rdata")
dat3 <- out[[4]]

load("biweekly-test-TAT4-LOD107.Rdata")
dat4 <- out[[4]]

dat <- rbind(dat, dat1,dat2, dat3, dat4)
rm(dat1,dat2, dat3, dat4)
unique(dat$LOD)
unique(dat$TAT)
unique(dat$test_rotation)
unique(dat$distance_limit)



load("biweekly-test-TAT5-LOD10.Rdata")
dat1 <- out[[4]]

load("biweekly-test-TAT5-LOD103.Rdata")
dat2 <- out[[4]]

load("biweekly-test-TAT5-LOD105.Rdata")
dat3 <- out[[4]]

load("biweekly-test-TAT5-LOD107.Rdata")
dat4 <- out[[4]]

dat <- rbind(dat, dat1,dat2, dat3, dat4)
rm(dat1,dat2, dat3, dat4)
unique(dat$LOD)
unique(dat$TAT)
unique(dat$test_rotation)
unique(dat$distance_limit)


load("biweekly-test-TAT10-LOD10.Rdata")
dat1 <- out[[4]]

load("biweekly-test-TAT10-LOD103.Rdata")
dat2 <- out[[4]]

load("biweekly-test-TAT10-LOD105.Rdata")
dat3 <- out[[4]]

load("biweekly-test-TAT10-LOD107.Rdata")
dat4 <- out[[4]]

dat <- rbind(dat, dat1,dat2, dat3, dat4)
rm(dat1,dat2, dat3, dat4)
unique(dat$LOD)
unique(dat$TAT)
unique(dat$test_rotation)
unique(dat$distance_limit)


#and now weekly
load("weekly-test-TAT1-LOD10.Rdata")
dat1 <- out[[4]]

load("weekly-test-TAT1-LOD103.Rdata")
dat2 <- out[[4]]

load("weekly-test-TAT1-LOD105.Rdata")
dat3 <- out[[4]]

load("weekly-test-TAT1-LOD107.Rdata")
dat4 <- out[[4]]

dat <- rbind(dat, dat1,dat2, dat3, dat4)
rm(dat1,dat2, dat3, dat4)
unique(dat$LOD)
unique(dat$TAT)
unique(dat$test_rotation)
unique(dat$distance_limit)


load("weekly-test-TAT2-LOD10.Rdata")
dat1 <- out[[4]]

load("weekly-test-TAT2-LOD103.Rdata")
dat2 <- out[[4]]

load("weekly-test-TAT2-LOD105.Rdata")
dat3 <- out[[4]]

load("weekly-test-TAT2-LOD107.Rdata")
dat4 <- out[[4]]

dat <- rbind(dat, dat1,dat2, dat3, dat4)
rm(dat1,dat2, dat3, dat4)
unique(dat$LOD)
unique(dat$TAT)
unique(dat$test_rotation)
unique(dat$distance_limit)


load("weekly-test-TAT3-LOD10.Rdata")
dat1 <- out[[4]]

load("weekly-test-TAT3-LOD103.Rdata")
dat2 <- out[[4]]

load("weekly-test-TAT3-LOD105.Rdata")
dat3 <- out[[4]]

load("weekly-test-TAT3-LOD107.Rdata")
dat4 <- out[[4]]

dat <- rbind(dat, dat1, dat2, dat3, dat4)
rm(dat1, dat2, dat3, dat4)
unique(dat$LOD)
unique(dat$TAT)
unique(dat$test_rotation)
unique(dat$distance_limit)

load("weekly-test-TAT4-LOD10.Rdata")
dat1 <- out[[4]]

load("weekly-test-TAT4-LOD103.Rdata")
dat2 <- out[[4]]

load("weekly-test-TAT4-LOD105.Rdata")
dat3 <- out[[4]]

load("weekly-test-TAT4-LOD107.Rdata")
dat4 <- out[[4]]

dat <- rbind(dat, dat1,dat2, dat3, dat4)
rm(dat1,dat2, dat3, dat4)
unique(dat$LOD)
unique(dat$TAT)
unique(dat$test_rotation)
unique(dat$distance_limit)



load("weekly-test-TAT5-LOD10.Rdata")
dat1 <- out[[4]]

load("weekly-test-TAT5-LOD103.Rdata")
dat2 <- out[[4]]

load("weekly-test-TAT5-LOD105.Rdata")
dat3 <- out[[4]]

load("weekly-test-TAT5-LOD107.Rdata")
dat4 <- out[[4]]

dat <- rbind(dat, dat1,dat2, dat3, dat4)
rm(dat1,dat2, dat3, dat4)
unique(dat$LOD)
unique(dat$TAT)
unique(dat$test_rotation)
unique(dat$distance_limit)


load("weekly-test-TAT10-LOD10.Rdata")
dat1 <- out[[4]]

load("weekly-test-TAT10-LOD103.Rdata")
dat2 <- out[[4]]

load("weekly-test-TAT10-LOD105.Rdata")
dat3 <- out[[4]]

load("weekly-test-TAT10-LOD107.Rdata")
dat4 <- out[[4]]

dat <- rbind(dat, dat1,dat2, dat3, dat4)
rm(dat1,dat2, dat3, dat4)
unique(dat$LOD)
unique(dat$TAT)
unique(dat$test_rotation)
unique(dat$distance_limit)



#and two-week
load("two-week-test-TAT1-LOD10.Rdata")
dat1 <- out[[4]]

load("two-week-test-TAT1-LOD103.Rdata")
dat2 <- out[[4]]

load("two-week-test-TAT1-LOD105.Rdata")
dat3 <- out[[4]]

load("two-week-test-TAT1-LOD107.Rdata")
dat4 <- out[[4]]

dat <- rbind(dat, dat1,dat2, dat3, dat4)
rm(dat1,dat2, dat3, dat4)
unique(dat$LOD)
unique(dat$TAT)
unique(dat$test_rotation)
unique(dat$distance_limit)


load("two-week-test-TAT2-LOD10.Rdata")
dat1 <- out[[4]]

load("two-week-test-TAT2-LOD103.Rdata")
dat2 <- out[[4]]

load("two-week-test-TAT2-LOD105.Rdata")
dat3 <- out[[4]]

load("two-week-test-TAT2-LOD107.Rdata")
dat4 <- out[[4]]

dat <- rbind(dat, dat1,dat2, dat3, dat4)
rm(dat1,dat2, dat3, dat4)
unique(dat$LOD)
unique(dat$TAT)
unique(dat$test_rotation)
unique(dat$distance_limit)


load("two-week-test-TAT3-LOD10.Rdata")
dat1 <- out[[4]]

load("two-week-test-TAT3-LOD103.Rdata")
dat2 <- out[[4]]

load("two-week-test-TAT3-LOD105.Rdata")
dat3 <- out[[4]]

load("two-week-test-TAT3-LOD107.Rdata")
dat4 <- out[[4]]

dat <- rbind(dat, dat1,dat2, dat3, dat4)
rm(dat1,dat2, dat3, dat4)
unique(dat$LOD)
unique(dat$TAT)
unique(dat$test_rotation)
unique(dat$distance_limit)



load("two-week-test-TAT4-LOD10.Rdata")
dat1 <- out[[4]]

load("two-week-test-TAT4-LOD103.Rdata")
dat2 <- out[[4]]

load("two-week-test-TAT4-LOD105.Rdata")
dat3 <- out[[4]]

load("two-week-test-TAT4-LOD107.Rdata")
dat4 <- out[[4]]

dat <- rbind(dat, dat1,dat2, dat3, dat4)
rm(dat1,dat2, dat3, dat4)
unique(dat$LOD)
unique(dat$TAT)
unique(dat$test_rotation)
unique(dat$distance_limit)



load("two-week-test-TAT5-LOD10.Rdata")
dat1 <- out[[4]]

load("two-week-test-TAT5-LOD103.Rdata")
dat2 <- out[[4]]

load("two-week-test-TAT5-LOD105.Rdata")
dat3 <- out[[4]]

load("two-week-test-TAT5-LOD107.Rdata")
dat4 <- out[[4]]

dat <- rbind(dat, dat1,dat2, dat3, dat4)
rm(dat1,dat2, dat3, dat4)
unique(dat$LOD)
unique(dat$TAT)
unique(dat$test_rotation)
unique(dat$distance_limit)


load("two-week-test-TAT10-LOD10.Rdata")
dat1 <- out[[4]]

load("two-week-test-TAT10-LOD103.Rdata")
dat2 <- out[[4]]

load("two-week-test-TAT10-LOD105.Rdata")
dat3 <- out[[4]]

load("two-week-test-TAT10-LOD107.Rdata")
dat4 <- out[[4]]

dat <- rbind(dat, dat1,dat2, dat3, dat4)
rm(dat1,dat2, dat3, dat4)
unique(dat$LOD)
unique(dat$TAT)
unique(dat$test_rotation)
unique(dat$distance_limit)

#and add in the symptom lag
dat$iso_lag = Inf


#and symptom stuff
setwd("/Users/caraebrook/Documents/R/R_repositories/Berkeley-Reopening/final-model-runs/Fig3-Symptom-Iso-Lags/out/")

load("symptom-iso-1.Rdata")
dat1 <- out[[4]]
dat1$iso_lag = 1

load("symptom-iso-2.Rdata")
dat2 <- out[[4]]
dat2$iso_lag = 2

load("symptom-iso-3.Rdata")
dat3 <- out[[4]]
dat3$iso_lag = 3

load("symptom-iso-4.Rdata")
dat4 <- out[[4]]
dat4$iso_lag = 4

load("symptom-iso-5.Rdata")
dat5 <- out[[4]]
dat5$iso_lag = 5


dat <- rbind(dat, dat1,dat2, dat3, dat4, dat5)
rm(dat1,dat2, dat3, dat4, dat5)
unique(dat$LOD)
unique(dat$TAT)
unique(dat$test_rotation)
unique(dat$distance_limit)
unique(dat$iso_lag)

#add group limits in later

dat.test.R0 <- dat
head(dat.test.R0)


dat.sum.list <- dlply(dat.test.R0 , .(test_rotation, iso_lag))

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
  #reduct.R$trace_lag  <- unique(dat1$trace_lag)
  reduct.R$mean[reduct.R$mean<0] <-0
  reduct.R$lci[reduct.R$lci<0] <-0
  reduct.R$uci[reduct.R$uci<0] <-0
  
  return(reduct.R)
  
}

dat.out.list <- lapply(dat.sum.list, get.rel.R0)

dat.test.R0  <- do.call("rbind", dat.out.list)
rownames(dat.test.R0 ) <- c()

dat.test.R0$LOD[dat.test.R0$LOD==10] <- "1e+01"
dat.test.R0$LOD[dat.test.R0$LOD==1000] <- "1e+03"

dat.test.R0.hold <- dat.test.R0

#and add in the group size limits
load("/Users/caraebrook/Documents/R/R_repositories/Berkeley-Reopening/final-model-runs/Fig2-Group-Limits/out/R0.group.lim.10.21.Rdata")
head(dat.test.R0)
dat.test.R0$iso_lag = Inf
head(dat.test.R0.hold)


dat.test.R0 <- rbind(dat.test.R0.hold, dat.test.R0)

unique(dat.test.R0$LOD)
unique(dat.test.R0$TAT)
unique(dat.test.R0$test_rotation)
unique(dat.test.R0$distance_limit)
unique(dat.test.R0$iso_lag)

dat.test.R0$intervention_class <- "testing"
dat.test.R0$intervention_class[dat.test.R0$iso_lag!=Inf] <- "symptom-isolation"
dat.test.R0$intervention_class[dat.test.R0$distance_limit!=Inf] <- "group-size-limit"
dat.test.R0$intervention_class[dat.test.R0$distance_limit==Inf & dat.test.R0$iso_lag==Inf & dat.test.R0$test_rotation=="none"] <- "none"

head(dat.test.R0)
#and save and plot.
setwd("/Users/caraebrook/Documents/R/R_repositories/Berkeley-Reopening/final-model-runs/Fig3-Testing-TAT-LOD-Freq/out/")
save(dat.test.R0, file = "R0.all.int.10.21.Rdata")
