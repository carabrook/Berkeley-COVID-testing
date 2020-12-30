rm(list=ls())

library(ggplot2)
library(plyr)
library(dplyr)

setwd("~/Documents/R/R_repositories/Berkeley-Reopening/Dec-2020/all-runs/Fig2-Group-Limits/out/")

load("group-lim-6.Rdata")
dat1 <- out[[1]]
load("group-lim-12.Rdata")
dat2 <- out[[1]]
load("group-lim-16.Rdata")
dat3 <- out[[1]]
load("group-lim-20.Rdata")
dat4 <- out[[1]]
load("group-lim-50.Rdata")
dat5 <- out[[1]]
load("no-group-lim.Rdata")
dat6 <- out[[1]]

dat <- rbind(dat1, dat2, dat3, dat4, dat5, dat6)
head(dat)

#and save 
save(dat, file = "dat.group.lim.12.21.Rdata")


#and the R0 data
rm(list=ls())
setwd("~/Documents/R/R_repositories/Berkeley-Reopening/Dec-2020/all-runs/Fig2-Group-Limits/out/")

load("group-lim-6.Rdata")
dat1 <- out[[4]]

load("group-lim-12.Rdata")
dat2 <- out[[4]]

load("group-lim-16.Rdata")
dat3 <- out[[4]]

load("group-lim-20.Rdata")
dat4 <- out[[4]]

load("group-lim-50.Rdata")
dat5 <- out[[4]]

load("group-lim-none.Rdata")
dat6 <- out[[4]]



dat <- rbind(dat1,dat2, dat3, dat4, dat5, dat6)
unique(dat$LOD)
unique(dat$test_rotation)
unique(dat$distance_limit)
dat$distance_limit = factor(dat$distance_limit, levels = c(6,12,16,20,50,Inf))





dat.test.R0 <- dat
head(dat.test.R0)


dat.test.R0$class <- rep(c("total_potential_cases", "UCB_potential_cases", "UCB_post_titer_potential_cases", 
                           "UCB_post_group_potential_cases", "UCB_post_isolations_actual_cases"), 6)

dat.sum.list <- dlply(dat.test.R0 , .(distance_limit))


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
  #reduct.R$iso_lag  <- unique(dat1$iso_lag)
  #reduct.R$trace_lag  <- unique(dat1$trace_lag)
  reduct.R$mean[reduct.R$mean<0] <-0
  reduct.R$lci[reduct.R$lci<0] <-0
  reduct.R$uci[reduct.R$uci<0] <-0
  
  
  return(reduct.R)
  
}

dat.out.list <- lapply(dat.sum.list, get.rel.R0)

dat.test.R0  <- do.call("rbind", dat.out.list)
rownames(dat.test.R0 ) <- c()

dat.test.R0 $LOD[dat.test.R0 $LOD==10] <- "1e+01"

dat.test.R0$distance_limit = factor(dat.test.R0$distance_limit, levels = c(6,12,16,20,50,Inf))



#and save and plot.
save(dat.test.R0, file = "R0.group.lim.12.21.Rdata")

