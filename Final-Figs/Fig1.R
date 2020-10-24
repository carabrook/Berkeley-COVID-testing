#make final figures
library(plyr)
library(ggplot2)
library(dplyr)


rm(list=ls())

load("Final-Figs/titer.dat.20K.Rdata")
#load("titer.dat.20K.Rdata")
load("Final-Figs/single.pop.par.Rdata")


#Figure 1 within-host trajectories
lognormal_fn <- function(meanlogpar=NULL, sdlogpar=NULL){
  out <- purrr::partial(rlnorm,
                        meanlog = meanlogpar,
                        sdlog = sdlogpar)
  return(out)
} 
get.titer.lim <- function(input.par, titer.dat){
  
  titer_lim = lognormal_fn(meanlogpar=log(as.numeric(input.par$par1[input.par$parameter=="symptom-lim"])), 
                           sdlogpar = log(as.numeric(input.par$par2[input.par$parameter=="symptom-lim"])))
  titer.vect <- titer_lim(length(unique(titer.dat$employ_ids)))
  
  #and save them and plot with the actual lines
  titer.dat$symptom_lim <- rep(titer.vect, each = max(titer.dat$time))
  
  return(titer.dat)
}
get.titer.limits <- function(meanpar, sdpar, N){
  
  titer_lim = lognormal_fn( meanlogpar =  log(meanpar),
                            sdlogpar  = log(sdpar))
  
  out=titer_lim(N)

  yplot = median(out)
  yQ1 = quantile(out)[2]
  yQ2 = quantile(out)[3]
  yQ3 = quantile(out)[4]
  ymax=max(out)
  #and combine
  out.dat <- cbind.data.frame(yplot, yQ1, yQ2, yQ3)
  names(out.dat) <- c("median_titer_lim", "Q1_titer_lim", "Q2_titer_lim", "Q3_titer_lim")
  return(out.dat)
}


make_Figure_1 <- function(filename){
titer.lim <- get.titer.limits(meanpar = 1e+7, sdpar = 10000, N= 20000)

titer.lim$x = 0
titer.lim <- rbind(titer.lim, titer.lim)
titer.lim$x[2] = 100

#summarize titer dat
titer.sum <- ddply(titer.dat, .(time), summarize, meanV= mean(V), sdV=sd(V))
titer.sum$lciV = titer.sum$meanV-titer.sum$sdV
titer.sum$lciV[titer.sum$lciV<0] <-0
titer.sum$uciV = titer.sum$meanV+titer.sum$sdV

p1 <- ggplot(data=titer.sum)+ 
      #geom_ribbon(data=titer.lim, aes(x=x, ymin=Q1_titer_lim, ymax=Q3_titer_lim), alpha=.2, fill="darkgreen") +
      geom_line(aes(x=time, y=meanV), color="cornflowerblue") + 
      geom_ribbon(aes(x=time, ymin=lciV, ymax=uciV), fill = "cornflowerblue", alpha=.3) +
      theme_bw() + theme(panel.grid = element_blank(), axis.title = element_text(size=18), axis.text = element_text(size=16)) + 
      coord_cartesian(xlim=c(0,10), ylim = c(0, 5e9)) + ylab("viral titer") + xlab( "days post-exposure") +
      geom_line(data=titer.lim, aes(x=x, y = median_titer_lim), color="magenta", linetype=2) +
      geom_line(data=titer.lim, aes(x=x, y = Q3_titer_lim), color="magenta", linetype=3) +
      annotate("text", label="75th percentile symptom onset", x =2, y=4.8e9, color="magenta") +
      annotate("text", label="50th percentile\nsymptom onset", x = 1.1, y=4e8, color="magenta") 
      
print(p1)

ggsave(file = filename,
       units="mm",  
       width=60, 
       height=40, 
       scale=3, 
       dpi=200)

}


make_Figure_1("Final-Figs/Fig1-inset.png") 

