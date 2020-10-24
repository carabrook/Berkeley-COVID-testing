rm(list=ls())

library(plyr)
library(dplyr)

#load goup testing data with LOD
setwd("/Users/caraebrook/Documents/R/R_repositories/Berkeley-Reopening/final-model-runs/Fig4-Layered-Interventions/out/")

load("biweekly-2-test-days-trace-1-symp-group.Rdata")
dat1 <- out[[1]]
dat1$iso_lag = 1
dat1$trace_lag = 1
load("biweekly-2-test-days-trace-1-symp.Rdata")
dat2 <- out[[1]]
dat2$iso_lag = 1
dat2$trace_lag = 1
load("biweekly-2-test-days-trace-1.Rdata")
dat3 <- out[[1]]
dat3$iso_lag = Inf
dat3$trace_lag = 1
load("biweekly-2-test-days.Rdata")
dat4 <- out[[1]]
dat4$iso_lag = Inf
dat4$trace_lag = Inf
load("biweekly-5-test-days-trace-1-symp-group.Rdata")
dat5 <- out[[1]]
dat5$iso_lag = 1
dat5$trace_lag = 1
load("biweekly-5-test-days-trace-1-symp.Rdata")
dat6 <- out[[1]]
dat6$iso_lag = 1
dat6$trace_lag = 1
load("biweekly-5-test-days-trace-1.Rdata")
dat7 <- out[[1]]
dat7$iso_lag = Inf
dat7$trace_lag = 1
load("biweekly-5-test-days.Rdata")
dat8 <- out[[1]]
dat8$iso_lag = Inf
dat8$trace_lag = Inf
load("biweekly-7-test-days-trace-1-symp-group.Rdata")
dat9 <- out[[1]]
dat9$iso_lag = 1
dat9$trace_lag = 1
load("biweekly-7-test-days-trace-1-symp.Rdata")
dat10 <- out[[1]]
dat10$iso_lag = 1
dat10$trace_lag = 1
load("biweekly-7-test-days-trace-1.Rdata")
dat11 <- out[[1]]
dat11$iso_lag = Inf
dat11$trace_lag = 1
load("biweekly-7-test-days.Rdata")
dat12 <- out[[1]]
dat12$iso_lag = Inf
dat12$trace_lag = Inf


dat <- rbind(dat1,dat2, dat3, dat4, dat5, dat6,dat7, dat8, dat9, dat10, dat11, dat12)
rm(dat1,dat2, dat3, dat4, dat5, dat6,dat7,dat8, dat9, dat10, dat11, dat12)
unique(dat$LOD)
unique(dat$test_rotation)
unique(dat$distance_limit)
unique(dat$iso_lag)
unique(dat$trace_lag)

load("weekly-2-test-days-trace-1-symp-group.Rdata")
dat1 <- out[[1]]
dat1$iso_lag = 1
dat1$trace_lag = 1
load("weekly-2-test-days-trace-1-symp.Rdata")
dat2 <- out[[1]]
dat2$iso_lag = 1
dat2$trace_lag = 1
load("weekly-2-test-days-trace-1.Rdata")
dat3 <- out[[1]]
dat3$iso_lag = Inf
dat3$trace_lag = 1
load("weekly-2-test-days.Rdata")
dat4 <- out[[1]]
dat4$iso_lag = Inf
dat4$trace_lag = Inf
load("weekly-5-test-days-trace-1-symp-group.Rdata")
dat5 <- out[[1]]
dat5$iso_lag = 1
dat5$trace_lag = 1
load("weekly-5-test-days-trace-1-symp.Rdata")
dat6 <- out[[1]]
dat6$iso_lag = 1
dat6$trace_lag = 1
load("weekly-5-test-days-trace-1.Rdata")
dat7 <- out[[1]]
dat7$iso_lag = Inf
dat7$trace_lag = 1
load("weekly-5-test-days.Rdata")
dat8 <- out[[1]]
dat8$iso_lag = Inf
dat8$trace_lag = Inf
load("weekly-7-test-day-trace-1-symp-group.Rdata")
dat9 <- out[[1]]
dat9$iso_lag = 1
dat9$trace_lag = 1
load("weekly-7-test-day-trace-1-symp.Rdata")
dat10 <- out[[1]]
dat10$iso_lag = 1
dat10$trace_lag = 1
load("weekly-7-test-day-trace-1.Rdata")
dat11 <- out[[1]]
dat11$iso_lag = Inf
dat11$trace_lag = 1
load("weekly-7-test-days.Rdata")
dat12 <- out[[1]]
dat12$iso_lag = Inf
dat12$trace_lag = Inf


dat <- rbind(dat, dat1, dat2, dat3, dat4, dat5, dat6,dat7, dat8,dat9, dat10, dat11,  dat12)
rm(dat1,dat2, dat3, dat4, dat5, dat6,dat7,dat8, dat9, dat10, dat11, dat12)
unique(dat$LOD)
unique(dat$test_rotation)
unique(dat$distance_limit)


load("two-week-2-test-days-trace-1-symp-group.Rdata")
dat1 <- out[[1]]
dat1$iso_lag = 1
dat1$trace_lag = 1
load("two-week-2-test-days-trace-1-symp.Rdata")
dat2 <- out[[1]]
dat2$iso_lag = 1
dat2$trace_lag = 1
load("two-week-2-test-days-trace-1.Rdata")
dat3 <- out[[1]]
dat3$iso_lag = Inf
dat3$trace_lag = 1
load("two-week-2-test-days.Rdata")
dat4 <- out[[1]]
dat4$iso_lag = Inf
dat4$trace_lag = Inf
load("two-week-5-test-days-trace-symp-group.Rdata")
dat5 <- out[[1]]
dat5$iso_lag = 1
dat5$trace_lag = 1
load("two-week-5-test-days-trace-symp.Rdata")
dat6 <- out[[1]]
dat6$iso_lag = 1
dat6$trace_lag = 1
load("two-week-5-test-days-trace.Rdata")
dat7 <- out[[1]]
dat7$iso_lag = Inf
dat7$trace_lag = 1
load("two-week-5-test-days.Rdata")
dat8 <- out[[1]]
dat8$iso_lag = Inf
dat8$trace_lag = Inf
load("two-week-7-test-days-trace-1-symp-group.Rdata")
dat9 <- out[[1]]
dat9$iso_lag = 1
dat9$trace_lag = 1
load("two-week-7-test-days-trace-1-symp.Rdata")
dat10 <- out[[1]]
dat10$iso_lag = 1
dat10$trace_lag = 1
load("two-week-7-test-days-trace-1.Rdata")
dat11 <- out[[1]]
dat11$iso_lag = Inf
dat11$trace_lag = 1
load("two-week-7-test-days.Rdata")
dat12 <- out[[1]]
dat12$iso_lag = Inf
dat12$trace_lag = Inf



dat <- rbind(dat, dat1, dat2, dat3, dat4, dat5, dat6,dat7, dat8,dat9, dat10, dat11,  dat12)
rm(dat1,dat2, dat3, dat4, dat5, dat6,dat7,dat8, dat9, dat10, dat11, dat12)
unique(dat$LOD)
unique(dat$test_rotation)
unique(dat$distance_limit)


dat$n_test_days <- NA
dat$n_test_days[dat$test_rotation=="biweekly-2-test-days"] <- 2
dat$n_test_days[dat$test_rotation=="weekly-2-test-days"] <- 2
dat$n_test_days[dat$test_rotation=="two-week-2-test-days"]  <- 2
dat$n_test_days[dat$test_rotation=="biweekly-5-test-days"] <- 5
dat$n_test_days[dat$test_rotation=="weekly-5-test-days"] <- 5
dat$n_test_days[dat$test_rotation=="two-week-5-test-days"]  <- 5
dat$n_test_days[dat$test_rotation=="biweekly-7-test-days"] <- 7
dat$n_test_days[dat$test_rotation=="weekly-7-test-days"] <- 7
dat$n_test_days[dat$test_rotation=="two-week-7-test-days"]  <- 7


dat$test_rotation[dat$test_rotation=="biweekly-2-test-days"] <- "biweekly"
dat$test_rotation[dat$test_rotation=="weekly-2-test-days"] <- "weekly"
dat$test_rotation[dat$test_rotation=="two-week-2-test-days"]  <- "two-week"
dat$test_rotation[dat$test_rotation=="biweekly-5-test-days"] <- "biweekly"
dat$test_rotation[dat$test_rotation=="weekly-5-test-days"] <- "weekly"
dat$test_rotation[dat$test_rotation=="two-week-5-test-days"]  <- "two-week"
dat$test_rotation[dat$test_rotation=="biweekly-7-test-days"] <- "biweekly"
dat$test_rotation[dat$test_rotation=="weekly-7-test-days"] <- "weekly"
dat$test_rotation[dat$test_rotation=="two-week-7-test-days"]  <- "two-week"

dat$test_rotation = factor(dat$test_rotation, levels = c("biweekly", "weekly", "two-week"))

dat$n_test_days = factor(dat$n_test_days, levels = c(2,5,7))

dat.test.all <- dat


#and save and plot.
save(dat.test.all, file = "dat.test.all.10.12.Rdata")



rm(list=ls())

#and the R0 data
setwd("/Users/caraebrook/Documents/R/R_repositories/Berkeley-Reopening/final-model-runs/Fig4-Layered-Interventions/out/")

load("biweekly-2-test-days-trace-1-symp-group.Rdata")
dat1 <- out[[4]]
dat1$iso_lag = 1
dat1$trace_lag = 1
load("biweekly-2-test-days-trace-1-symp.Rdata")
dat2 <- out[[4]]
dat2$iso_lag = 1
dat2$trace_lag = 1
load("biweekly-2-test-days-trace-1.Rdata")
dat3 <- out[[4]]
dat3$iso_lag = Inf
dat3$trace_lag = 1
load("biweekly-2-test-days.Rdata")
dat4 <- out[[4]]
dat4$iso_lag = Inf
dat4$trace_lag = Inf
load("biweekly-5-test-days-trace-1-symp-group.Rdata")
dat5 <- out[[4]]
dat5$iso_lag = 1
dat5$trace_lag = 1
load("biweekly-5-test-days-trace-1-symp.Rdata")
dat6 <- out[[4]]
dat6$iso_lag = 1
dat6$trace_lag = 1
load("biweekly-5-test-days-trace-1.Rdata")
dat7 <- out[[4]]
dat7$iso_lag = Inf
dat7$trace_lag = 1
load("biweekly-5-test-days.Rdata")
dat8 <- out[[4]]
dat8$iso_lag = Inf
dat8$trace_lag = Inf
load("biweekly-7-test-days-trace-1-symp-group.Rdata")
dat9 <- out[[4]]
dat9$iso_lag = 1
dat9$trace_lag = 1
load("biweekly-7-test-days-trace-1-symp.Rdata")
dat10 <- out[[4]]
dat10$iso_lag = 1
dat10$trace_lag = 1
load("biweekly-7-test-days-trace-1.Rdata")
dat11 <- out[[4]]
dat11$iso_lag = Inf
dat11$trace_lag = 1
load("biweekly-7-test-days.Rdata")
dat12 <- out[[4]]
dat12$iso_lag = Inf
dat12$trace_lag = Inf


dat <- rbind(dat1,dat2, dat3, dat4, dat5, dat6,dat7, dat8, dat9, dat10, dat11, dat12)
rm(dat1,dat2, dat3, dat4, dat5, dat6,dat7,dat8, dat9, dat10, dat11, dat12)
unique(dat$LOD)
unique(dat$test_rotation)
unique(dat$distance_limit)
unique(dat$iso_lag)
unique(dat$trace_lag)

load("weekly-2-test-days-trace-1-symp-group.Rdata")
dat1 <- out[[4]]
dat1$iso_lag = 1
dat1$trace_lag = 1
load("weekly-2-test-days-trace-1-symp.Rdata")
dat2 <- out[[4]]
dat2$iso_lag = 1
dat2$trace_lag = 1
load("weekly-2-test-days-trace-1.Rdata")
dat3 <- out[[4]]
dat3$iso_lag = Inf
dat3$trace_lag = 1
load("weekly-2-test-days.Rdata")
dat4 <- out[[4]]
dat4$iso_lag = Inf
dat4$trace_lag = Inf
load("weekly-5-test-days-trace-1-symp-group.Rdata")
dat5 <- out[[4]]
dat5$iso_lag = 1
dat5$trace_lag = 1
load("weekly-5-test-days-trace-1-symp.Rdata")
dat6 <- out[[4]]
dat6$iso_lag = 1
dat6$trace_lag = 1
load("weekly-5-test-days-trace-1.Rdata")
dat7 <- out[[4]]
dat7$iso_lag = Inf
dat7$trace_lag = 1
load("weekly-5-test-days.Rdata")
dat8 <- out[[4]]
dat8$iso_lag = Inf
dat8$trace_lag = Inf
load("weekly-7-test-day-trace-1-symp-group.Rdata")
dat9 <- out[[4]]
dat9$iso_lag = 1
dat9$trace_lag = 1
load("weekly-7-test-day-trace-1-symp.Rdata")
dat10 <- out[[4]]
dat10$iso_lag = 1
dat10$trace_lag = 1
load("weekly-7-test-day-trace-1.Rdata")
dat11 <- out[[4]]
dat11$iso_lag = Inf
dat11$trace_lag = 1
load("weekly-7-test-days.Rdata")
dat12 <- out[[4]]
dat12$iso_lag = Inf
dat12$trace_lag = Inf


dat <- rbind(dat, dat1, dat2, dat3, dat4, dat5, dat6,dat7, dat8,dat9, dat10, dat11,  dat12)
rm(dat1,dat2, dat3, dat4, dat5, dat6,dat7,dat8, dat9, dat10, dat11, dat12)
unique(dat$LOD)
unique(dat$test_rotation)
unique(dat$distance_limit)


load("two-week-2-test-days-trace-1-symp-group.Rdata")
dat1 <- out[[4]]
dat1$iso_lag = 1
dat1$trace_lag = 1
load("two-week-2-test-days-trace-1-symp.Rdata")
dat2 <- out[[4]]
dat2$iso_lag = 1
dat2$trace_lag = 1
load("two-week-2-test-days-trace-1.Rdata")
dat3 <- out[[4]]
dat3$iso_lag = Inf
dat3$trace_lag = 1
load("two-week-2-test-days.Rdata")
dat4 <- out[[4]]
dat4$iso_lag = Inf
dat4$trace_lag = Inf
load("two-week-5-test-days-trace-symp-group.Rdata")
dat5 <- out[[4]]
dat5$iso_lag = 1
dat5$trace_lag = 1
load("two-week-5-test-days-trace-symp.Rdata")
dat6 <- out[[4]]
dat6$iso_lag = 1
dat6$trace_lag = 1
load("two-week-5-test-days-trace.Rdata")
dat7 <- out[[4]]
dat7$iso_lag = Inf
dat7$trace_lag = 1
load("two-week-5-test-days.Rdata")
dat8 <- out[[4]]
dat8$iso_lag = Inf
dat8$trace_lag = Inf
load("two-week-7-test-days-trace-1-symp-group.Rdata")
dat9 <- out[[4]]
dat9$iso_lag = 1
dat9$trace_lag = 1
load("two-week-7-test-days-trace-1-symp.Rdata")
dat10 <- out[[4]]
dat10$iso_lag = 1
dat10$trace_lag = 1
load("two-week-7-test-days-trace-1.Rdata")
dat11 <- out[[4]]
dat11$iso_lag = Inf
dat11$trace_lag = 1
load("two-week-7-test-days.Rdata")
dat12 <- out[[4]]
dat12$iso_lag = Inf
dat12$trace_lag = Inf



dat <- rbind(dat, dat1, dat2, dat3, dat4, dat5, dat6,dat7, dat8,dat9, dat10, dat11,  dat12)
rm(dat1,dat2, dat3, dat4, dat5, dat6,dat7,dat8, dat9, dat10, dat11, dat12)
unique(dat$LOD)
unique(dat$test_rotation)
unique(dat$distance_limit)



dat$n_test_days <- NA
dat$n_test_days[dat$test_rotation=="biweekly-2-test-days"] <- 2
dat$n_test_days[dat$test_rotation=="weekly-2-test-days"] <- 2
dat$n_test_days[dat$test_rotation=="two-week-2-test-days"]  <- 2
dat$n_test_days[dat$test_rotation=="biweekly-5-test-days"] <- 5
dat$n_test_days[dat$test_rotation=="weekly-5-test-days"] <- 5
dat$n_test_days[dat$test_rotation=="two-week-5-test-days"]  <- 5
dat$n_test_days[dat$test_rotation=="biweekly-7-test-days"] <- 7
dat$n_test_days[dat$test_rotation=="weekly-7-test-days"] <- 7
dat$n_test_days[dat$test_rotation=="two-week-7-test-days"]  <- 7


dat$test_rotation[dat$test_rotation=="biweekly-2-test-days"] <- "biweekly"
dat$test_rotation[dat$test_rotation=="weekly-2-test-days"] <- "weekly"
dat$test_rotation[dat$test_rotation=="two-week-2-test-days"]  <- "two-week"
dat$test_rotation[dat$test_rotation=="biweekly-5-test-days"] <- "biweekly"
dat$test_rotation[dat$test_rotation=="weekly-5-test-days"] <- "weekly"
dat$test_rotation[dat$test_rotation=="two-week-5-test-days"]  <- "two-week"
dat$test_rotation[dat$test_rotation=="biweekly-7-test-days"] <- "biweekly"
dat$test_rotation[dat$test_rotation=="weekly-7-test-days"] <- "weekly"
dat$test_rotation[dat$test_rotation=="two-week-7-test-days"]  <- "two-week"

dat$test_rotation = factor(dat$test_rotation, levels = c("biweekly", "weekly", "two-week"))

dat$n_test_days = factor(dat$n_test_days, levels = c(2,5,7))


dat.test.R0 <- dat
head(dat.test.R0)

dat.sum.list <- dlply(dat.test.R0 , .(test_rotation, n_test_days, distance_limit, iso_lag, trace_lag))

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
  reduct.R$n_test_days  <- unique(dat1$n_test_days)
  reduct.R$iso_lag  <- unique(dat1$iso_lag)
  reduct.R$trace_lag  <- unique(dat1$trace_lag)
  
  return(reduct.R)
  
}

dat.out.list <- lapply(dat.sum.list, get.rel.R0)

dat.test.R0  <- do.call("rbind", dat.out.list)
rownames(dat.test.R0 ) <- c()

dat.test.R0 $LOD[dat.test.R0 $LOD==10] <- "1e+01"

#dat.test$test_rotation[dat.test$test_rotation=="two_day"] <- "every-other-day"
dat.test.R0$test_rotation <- factor(dat.test.R0$test_rotation, levels = c("biweekly", "weekly", "two-week"))

unique(dat.test.R0$test_rotation)


#and save and plot.
save(dat.test.R0, file = "R0_test_all_group_10_12.Rdata")
