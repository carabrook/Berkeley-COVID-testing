#make final figures
library(plyr)
library(ggplot2)
library(dplyr)


rm(list=ls())

#setwd("/Users/caraebrook/Documents/R/R_repositories/Berkeley-Reopening/Dec-2020/all-runs/")

load("Final-Figs/titer.dat.20K.Rdata")
load("Final-Figs/pop.par.base.Rdata")


#Figure 1 within-host trajectories
lognormal_fn <- function(meanlogpar=NULL, sdlogpar=NULL){
  out <- purrr::partial(rlnorm,
                        meanlog = meanlogpar,
                        sdlog = sdlogpar)
  return(out)
} 

make_Figure_1 <- function(titer.dat, pop.par, filename){
  
  
#titer  limit for symptoms
titer_lim = lognormal_fn(meanlogpar=log(as.numeric(pop.par$par1[pop.par$parameter=="symptom-lim"])), 
                           sdlogpar = log(as.numeric(pop.par$par2[pop.par$parameter=="symptom-lim"])))

#and get the IQR for 20000 cases
titer.lim.dat <- titer_lim(nrow(titer.dat))
titer.lim.dat <- quantile(titer.lim.dat)
  
y.50 <-log10(titer.lim.dat[3])
y.50.text = y.50 + .7
y.75 <-log10(titer.lim.dat[4])
y.75.text = y.75 + .3

p2 <- ggplot(titer.dat) + geom_line(aes(x=time, y=log10(V), color=piT, group=employ_ids), show.legend = F) + 
  coord_cartesian(xlim=c(0,20), ylim =c(0,9)) +
  theme_bw() + geom_hline(yintercept = 3, color="navy", linetype=2, size=1) + 
  theme(panel.grid = element_blank(), axis.title = element_text(size=26), axis.text = element_text(size=22))+
  annotate("text", label="LOD=10^3", y = 3.3, x=18, color="navy", size=8) + 
  ylab("viral load (cp/ul RNA, log10-scale)") + xlab("days since exposure") + scale_color_gradient(low="aliceblue", high="dodgerblue4") +
  geom_hline(yintercept = 5, color="navy", linetype=2, size=1) + 
  scale_y_continuous(breaks=c(0,1,3,5,7), labels=c("0","10^1","10^3", "10^5", "10^7"))+
  geom_hline(yintercept = y.50, color="deeppink3", linetype=2, size=1) + 
  geom_hline(yintercept = y.75, color="deeppink3", linetype=2, size=1) + 
  annotate("text", label="LOD=10^5", y = 5.3, x=18, color="navy", size=8) +
  annotate("text", label="mean symptom onset/", y = y.50.text, x=18, color="deeppink3", size=8) +
  annotate("text", label="75th percentile symptom onset", y = y.75.text, x=16, color="deeppink3", size=8)
  
print(p2)

ggsave(file = filename,
       units="mm",  
       width=60, 
       height=40, 
       scale=6, 
       dpi=300)

}


make_Figure_1(titer.dat=titer.dat, pop.par=pop.par.base, filename= "Final-Figs/Fig1-inset.png") 

