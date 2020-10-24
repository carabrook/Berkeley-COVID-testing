rm(list=ls())

library(plyr)
library(dplyr)

#load goup testing data with LOD
setwd("/Users/caraebrook/Documents/R/R_repositories/Berkeley-Reopening/final-model-runs/Fig5-multi-group/out/")

load("all-no-test.Rdata")
dat1 <- out[[1]]
load("biweekly-all.Rdata")
dat2 <- out[[1]]
load("biweekly-high-three-week-low.Rdata")
dat3 <- out[[1]]
load("biweekly-high-two-week-low.Rdata")
dat4 <- out[[1]]
load("biweekly-high-weekly-low.Rdata")
dat5 <- out[[1]]
load("two-week-all.Rdata")
dat6 <- out[[1]]
load("weekly-all.Rdata")
dat7 <- out[[1]]


dat <- rbind(dat1,dat2, dat3, dat4, dat5, dat6,dat7)
rm(dat1,dat2, dat3, dat4, dat5, dat6,dat7)
unique(dat$LOD)
unique(dat$test_rotation)
dat$test_rotation <- factor(dat$test_rotation, levels = c("biweekly-all", "weekly-all", "two-week-all", "biweekly-high-weekly-low", "biweekly-high-two-week-low", "biweekly-high-three-week-low", "all-no-test"))
unique(dat$distance_limit)

dat.multi <- dat

#and save and plot.
save(dat.multi, file = "dat_multi_10_12.Rdata")



rm(list=ls())

#and the R0 data
setwd("/Users/caraebrook/Documents/R/R_repositories/Berkeley-Reopening/final-model-runs/Fig5-multi-group/out/")

load("all-no-test.Rdata")
dat1 <- out[[4]]
load("biweekly-all.Rdata")
dat2 <- out[[4]]
load("biweekly-high-three-week-low.Rdata")
dat3 <- out[[4]]
load("biweekly-high-two-week-low.Rdata")
dat4 <- out[[4]]
load("biweekly-high-weekly-low.Rdata")
dat5 <- out[[4]]
load("two-week-all.Rdata")
dat6 <- out[[4]]
load("weekly-all.Rdata")
dat7 <- out[[4]]


dat <- rbind(dat1,dat2, dat3, dat4, dat5, dat6,dat7)
rm(dat1,dat2, dat3, dat4, dat5, dat6,dat7)
unique(dat$LOD)
unique(dat$test_rotation)
dat$test_rotation <- factor(dat$test_rotation, levels = c("biweekly-all", "weekly-all", "two-week-all", "biweekly-high-weekly-low", "biweekly-high-two-week-low", "biweekly-high-three-week-low", "all-no-test"))
unique(dat$distance_limit)

dat.multi.R0 <- dat

ref.R0 = dat.multi.R0$mean.all[dat.multi.R0$class=="UCB_post_isolations_actual_cases" & dat$test_rotation=="all-no-test"]
ref.R0.lci = dat.multi.R0$lci.all[dat.multi.R0$class=="UCB_post_isolations_actual_cases" & dat$test_rotation=="all-no-test"]
ref.R0.uci = dat.multi.R0$uci.all[dat.multi.R0$class=="UCB_post_isolations_actual_cases" & dat$test_rotation=="all-no-test"]

dat.sum.list <- dlply(dat.multi.R0, .(test_rotation))

get.rel.R0.multi <- function(dat1, refmean, reflci, refuci){
  start.R = dat1$mean.all[dat1$class=="UCB_post_titer_potential_cases"]
  start.R.lci = dat1$lci.all[dat1$class=="UCB_post_titer_potential_cases"]
  start.R.uci = dat1$uci.all[dat1$class=="UCB_post_titer_potential_cases"]
  
  finish.R = dat1$mean.all[dat1$class=="UCB_post_isolations_actual_cases"]
  finish.R.lci = dat1$lci.all[dat1$class=="UCB_post_isolations_actual_cases"]
  finish.R.uci = dat1$uci.all[dat1$class=="UCB_post_isolations_actual_cases"]
  
  out.R1 <- cbind.data.frame(refmean, reflci, refuci)
  out.R2 <- cbind.data.frame(finish.R,  finish.R.lci,  finish.R.uci)
  names(out.R1) <- names(out.R2) <- c("mean", "lci", "uci")  
  out.R <- rbind(out.R1, out.R2)
  out.R$intervention <- c("pre", "post")
  
  reduct.R = out.R1-out.R2
  reduct.R$LOD  <- unique(dat1$LOD)
  #reduct.R$TAT  <- unique(dat1$TAT)
  reduct.R$test_rotation  <- unique(dat1$test_rotation)
  reduct.R$distance_limit  <- unique(dat1$distance_limit)
  #reduct.R$n_test_days  <- unique(dat1$n_test_days)
  
  return(reduct.R)
  
}

dat.out.list <- lapply(dat.sum.list, get.rel.R0.multi, refmean=ref.R0 , reflci=ref.R0.lci, refuci = ref.R0.uci)

dat.multi.R0 <- do.call("rbind", dat.out.list)
rownames(dat.multi.R0) <- c()

dat.multi.R0$LOD[dat.multi.R0$LOD==10] <- "1e+01"

#dat.test$test_rotation[dat.test$test_rotation=="two_day"] <- "every-other-day"
dat.multi.R0$test_rotation <- factor(dat.multi.R0$test_rotation, levels = c("biweekly-all", "weekly-all", "two-week-all", "biweekly-high-weekly-low", "biweekly-high-two-week-low", "biweekly-high-three-week-low", "all-no-test"))

unique(dat.multi.R0$test_rotation)


#and save and plot.
save(dat.multi.R0, file = "R0_multi_group_10_12.Rdata")
