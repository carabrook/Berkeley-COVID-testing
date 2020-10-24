rm(list = ls())
library(ggplot2)
library(plyr)
library(dplyr)


#Figure 3 shows gains in testing across a set distance limit - compare frequency/TAT/LOD

make_Fig_3_S2 <- function(filename){
  
load("Final-Figs/dat.trace.10.21.Rdata")

head(dat.trace)
dat.trace$test_rotation <- sapply(strsplit(dat.trace$test_rotation, "-"), '[', 1)
dat.trace$test_rotation[dat.trace$test_rotation=="biweekly"] <- "biweekly-testing"
dat.trace$test_rotation[dat.trace$test_rotation=="weekly"] <- "weekly-testing"
dat.trace$test_rotation[dat.trace$test_rotation=="two"] <- "every-two-weeks-testing"
dat.trace$test_rotation <- factor(dat.trace$test_rotation, levels = c("biweekly-testing", "weekly-testing", "every-two-weeks-testing"))


#dat.trace$trace_lag <- factor(dat.trace$trace_lag, levels=c(1, Inf))

dat.all = subset(dat.trace, day <=50)

#first, summarize by intervention
dat.sum <- ddply(subset(dat.all, day==50), .(trace_lag, intervention_class, test_rotation, TAT, iso_lag, trace_lag), summarize, total_cases=unique(mean_cumulative), lci_cases=unique(lci_cumulative), uci_cases=(uci_cumulative))

#and get cases saved
dat.sum$cases_saved = dat.sum$total_cases[dat.sum$intervention_class=="none" & dat.sum$trace_lag==1] - dat.sum$total_cases 
dat.sum$cases_saved_lci = dat.sum$lci_cases[dat.sum$intervention_class=="none" & dat.sum$trace_lag==1] - dat.sum$lci_cases 
dat.sum$cases_saved_uci = dat.sum$uci_cases[dat.sum$intervention_class=="none" & dat.sum$trace_lag==1] - dat.sum$uci_cases 



#dat.sum[dat.sum<0]=0

dat.sum$TAT = factor(dat.sum$TAT, levels = c(1,2,3,4,5,10, Inf))
dat.sum$iso_lag = factor(dat.sum$iso_lag, levels=c(1,2,3,4,5, Inf))
dat.sum$trace_lag = factor(dat.sum$trace_lag, levels = c(1, Inf))

dat.sum = subset(dat.sum, intervention_class!="none")


fillz <- c("1" = "tomato", "2" = "goldenrod", "3" = "springgreen3", "4"="turquoise3", "5"="royalblue", "10" = "magenta")



dat.sum$label = NA
dat.sum$label[dat.sum$trace_lag==Inf] <- "no contact tracing"
dat.sum$label[dat.sum$trace_lag==1] <- "with contact tracing"
dat.sum$label <- factor(dat.sum$label, levels = c("no contact tracing", "with contact tracing"))


dat.sum.sym = subset(dat.sum, intervention_class=="symptom-isolation")
dat.sum.test = subset(dat.sum, intervention_class=="testing")


p1a <- ggplot(dat.sum.sym)  + ylab("cases saved")   + facet_grid(label~intervention_class) +  scale_fill_manual(values=fillz) + scale_color_manual(values=fillz) + guides(fill="none") +
  geom_bar(aes(x=iso_lag, y=cases_saved, fill = iso_lag), stat = "identity", position=position_dodge(),  color="black", show.legend = F) + coord_cartesian(ylim=c(-1000,12000)) +
  geom_errorbar(aes(x=iso_lag, ymin=cases_saved_lci, ymax=cases_saved_uci, color=iso_lag), width=.3, position=position_dodge(width = .9), show.legend = F)  +
  theme_bw()+ theme(legend.position = c(.9,.85), plot.margin = unit(c(1.5,0,.5,.5), "lines"), panel.grid = element_blank(), strip.background = element_rect(fill="white"), #panel.spacing = unit(3, "lines"), 
                    axis.title = element_text(size=14),   axis.text = element_text(size=12), strip.text.x = element_text(size=11), strip.text.y = element_blank()) + xlab("lag to isolation (days)")

print(p1a)

dat.sum.test$test_rotation = factor(dat.sum.test$test_rotation, levels = c("biweekly-testing", "weekly-testing", "every-two-weeks-testing"))

p1b <- ggplot(dat.sum.test)  + ylab("cases saved")   + facet_grid(label~test_rotation) +  scale_fill_manual(values=fillz) + scale_color_manual(values=fillz) + guides(fill="none") +
  geom_bar(aes(x=TAT, y=cases_saved, fill = TAT), stat = "identity", position=position_dodge(),  color="black", show.legend = F) + coord_cartesian(ylim=c(-1000,12000)) +
  geom_errorbar(aes(x=TAT, ymin=cases_saved_lci, ymax=cases_saved_uci, color=TAT), width=.3, position=position_dodge(width = .9), show.legend = F)  +
  theme_bw()+ theme(legend.position = c(.9,.85), plot.margin = unit(c(1.5,.5,.5,.2), "lines"), panel.grid = element_blank(), strip.background = element_rect(fill="white"), #panel.spacing = unit(3, "lines"), 
                    axis.title.x = element_text(size=14),   axis.text.x = element_text(size=12), axis.text.y = element_blank(), axis.title.y = element_blank(), axis.ticks.y = element_blank(),
                    strip.text = element_text(size=11)) + xlab("testing mean turnaround time(days)")

print(p1b)



FigS1top <- cowplot::plot_grid(p1a, p1b, nrow=1, ncol=2, rel_widths = c(1.4,3))


#and the R0 reduction


print(FigS1top)
#and then behavioral interventions

#next plot R-effective across each intervention and compare

load("Final-Figs/R0.trace.10.21.Rdata")

dat.trace.R0$test_rotation <- sapply(strsplit(dat.trace.R0$test_rotation, "-"), '[', 1)
dat.trace.R0$test_rotation[dat.trace.R0$test_rotation=="biweekly"] <- "biweekly-testing"
dat.trace.R0$test_rotation[dat.trace.R0$test_rotation=="weekly"] <- "weekly-testing"
dat.trace.R0$test_rotation[dat.trace.R0$test_rotation=="two"] <- "every-two-weeks-testing"
dat.trace.R0$test_rotation <- factor(dat.trace.R0$test_rotation, levels = c("biweekly-testing", "weekly-testing", "every-two-weeks-testing"))

dat.trace.R0$label = "with contact tracing"
dat.trace.R0$label[dat.trace.R0$trace_lag==Inf] <- "no contact tracing"

dat.trace.R0$label = factor(dat.trace.R0$label, levels = c("no contact tracing", "with contact tracing"))

#now plot testing side-by-side with behavior
dat.test = subset(dat.trace.R0, intervention_class=="testing")
dat.symp = subset(dat.trace.R0, intervention_class=="symptom-isolation")




fillz <- c("1" = "tomato", "2" = "goldenrod", "3" = "springgreen3", "4"="turquoise3", "5"="royalblue", "10" = "magenta")

p.2a <- ggplot(data=dat.test) + geom_point(aes(x=TAT, y=mean, color=TAT), shape = 16, size=5,show.legend = F,  position=position_dodge(width=.5)) + 
  scale_color_manual(values=fillz) + 
  ylab("reduction in Reff from intervention") + xlab("testing mean turnaround time (days)") + facet_grid(label~test_rotation) +
  geom_errorbar(aes(x=TAT, ymin=lci, max=uci,  color=TAT), size=.1, show.legend = F,position=position_dodge(width=.5)) +
  theme_bw() + theme(panel.grid = element_blank(), axis.text.x= element_text(size = 12), axis.title.x = element_text(size=14), axis.title.y = element_blank(),
                     strip.text = element_text(size=11), axis.text.y = element_blank(), axis.ticks.y = element_blank(),
                     strip.background = element_rect(fill="white"), plot.margin = unit(c(1.5,.5,.5,0), "lines")) #+ geom_label(aes(x=behavioral_int, y=mean, label=behavioral_int))
print(p.2a)



p.2b <- ggplot(data=dat.symp) + geom_point(aes(x=iso_lag, y=mean, color=iso_lag), shape = 16, size=5,show.legend = F,  position=position_dodge(width=.5)) + 
  scale_color_manual(values=fillz) + 
  ylab("reduction in Reff from intervention") + xlab("testing mean turnaround time (days)") + facet_grid(label~intervention_class) +
  geom_errorbar(aes(x=iso_lag, ymin=lci, max=uci,  color=iso_lag), size=.1, show.legend = F,position=position_dodge(width=.5)) +
  theme_bw() + theme(panel.grid = element_blank(), axis.text = element_text(size = 12), axis.title = element_text(size=14), 
                     strip.text.x = element_text(size=11), strip.text.y = element_blank(),
                     strip.background = element_rect(fill="white"), plot.margin = unit(c(1.5,.2,.5,.5), "lines")) #+ geom_label(aes(x=behavioral_int, y=mean, label=behavioral_int))
print(p.2b)




FigS1bottom <- cowplot::plot_grid(p.2b, p.2a, nrow=1, ncol=2, rel_widths = c(1.2,3))

print(FigS1bottom)




#and together 
FigS1 <- cowplot::plot_grid(FigS1bottom, FigS1top, nrow=1, ncol=2, rel_widths = c(1,1), labels=c("A.", "B."), hjust = -.1)

print(FigS1)

ggsave(file = filename,
       units="mm",  
       width=140, 
       height=60, 
       scale=3, 
       dpi=300)
}

make_Fig_3_S2("Final-Figs/FigS3.png")


#and the second data set
make.Table.S2 <- function(filename){
  
  load("Final-Figs/dat.trace.10.21.Rdata")
  head(dat.trace)
  
  
  dat = subset(dat.trace, day <=50)
  
  #first, summarize by intervention
  dat.sum <- ddply(subset(dat, day==50), .(intervention_class, distance_limit, test_rotation, LOD, TAT, iso_lag, trace_lag), summarize, total_cases=unique(mean_cumulative), lci_cases=unique(lci_cumulative), uci_cases=(uci_cumulative))
  
  #and get cases saved
  dat.sum$cases_saved = dat.sum$total_cases[dat.sum$intervention_class=="none" & dat.sum$trace_lag==Inf] - dat.sum$total_cases 
  dat.sum$cases_saved_lci = dat.sum$lci_cases[dat.sum$intervention_class=="none"& dat.sum$trace_lag==Inf] - dat.sum$lci_cases 
  dat.sum$cases_saved_uci = dat.sum$uci_cases[dat.sum$intervention_class=="none"& dat.sum$trace_lag==Inf] - dat.sum$uci_cases 
  
  
  dat.sum$TAT = factor(dat.sum$TAT, levels = c(1,2,3,4,5,10, Inf))
  dat.sum$LOD = as.character(dat.sum$LOD)
  dat.sum$LOD[dat.sum$LOD=="10"] <- "1e+01"
  dat.sum$LOD[dat.sum$LOD=="1000"] <- "1e+03"
  dat.sum$LOD = factor(dat.sum$LOD, levels = c("1e+01", "1e+03", "1e+05", "1e+07"))
  
  dat.sum$iso_lag = factor(dat.sum$iso_lag, levels=c(1,2,3,4,5, Inf))
  
  dat.sum$distance_limit = factor(dat.sum$distance_limit, levels=c(6,12,16, 20,50,Inf))
  dat.sum$contact_tracing = "no contact tracing"
  dat.sum$contact_tracing[dat.sum$trace_lag==1] <- "with contact tracing"
  
  dat.sum$contact_tracing <- factor(dat.sum$contact_tracing, levels = c( "no contact tracing", "with contact tracing"))
  
  load("Final-Figs/R0.trace.10.21.Rdata")
  
  dat.trace.R0$TAT = factor(dat.trace.R0$TAT, levels = c(1,2,3,4,5,10, Inf))
  dat.trace.R0$LOD = as.character(dat.trace.R0$LOD)
  dat.trace.R0$LOD = factor(dat.trace.R0$LOD, levels = c("1e+01", "1e+03", "1e+05", "1e+07"))
  
  dat.trace.R0$iso_lag = factor(dat.trace.R0$iso_lag, levels=c(1,2,3,4,5, Inf))
  
  dat.trace.R0$distance_limit = factor(dat.trace.R0$distance_limit, levels=c(6,12,16, 20,50,Inf))
  
  
  #and save the raw data:
  head(dat.sum)
  head(dat.trace.R0)
  
  
  
  #and add in Reff
  dat.sum$uci_Reff_reduction  <- dat.sum$lci_Reff_reduction <- dat.sum$mean_Reff_reduction <- NA
  
  dat.sum$intervention = dat.sum$test_rotation
  dat.sum$intervention[dat.sum$intervention_class=="group-size-limit"] <- paste0("group-lim=",as.character(dat.sum$distance_limit[dat.sum$intervention_class=="group-size-limit"]))
  dat.sum$intervention[dat.sum$intervention_class=="symptom-isolation"] <- paste0("symptom-iso=",as.character(dat.sum$iso_lag[dat.sum$intervention_class=="symptom-isolation"]))
  
  
  dat.trace.R0$intervention = dat.trace.R0$test_rotation
  dat.trace.R0$intervention[dat.trace.R0$intervention_class=="group-size-limit"] <- paste0("group-lim=",as.character(dat.trace.R0$distance_limit[dat.trace.R0$intervention_class=="group-size-limit"]))
  dat.trace.R0$intervention[dat.trace.R0$intervention_class=="symptom-isolation"] <- paste0("symptom-iso=",as.character(dat.trace.R0$iso_lag[dat.trace.R0$intervention_class=="symptom-isolation"]))
  
  
  
  for (i in 1:length(dat.trace.R0$mean)){
    dat.sum$mean_Reff_reduction[dat.sum$intervention==dat.trace.R0$intervention[i] & dat.sum$trace_lag==dat.trace.R0$trace_lag[i]] <- dat.trace.R0$mean[i]
    dat.sum$lci_Reff_reduction[dat.sum$intervention==dat.trace.R0$intervention[i]& dat.sum$trace_lag==dat.trace.R0$trace_lag[i]] <- dat.trace.R0$lci[i]
    dat.sum$uci_Reff_reduction[dat.sum$intervention==dat.trace.R0$intervention[i]& dat.sum$trace_lag==dat.trace.R0$trace_lag[i]] <- dat.trace.R0$uci[i]
  }
  
  
  dat.sum$test_rotation <- sapply(strsplit(dat.sum$test_rotation, "-"), '[', 1)
  dat.sum$test_rotation[dat.sum$test_rotation=="biweekly"] <- "biweekly-testing"
  dat.sum$test_rotation[dat.sum$test_rotation=="weekly"] <- "weekly-testing"
  dat.sum$test_rotation[dat.sum$test_rotation=="two"] <- "every-two-weeks-testing"
  dat.sum$test_rotation <- factor(dat.sum$test_rotation, levels = c("biweekly-testing", "weekly-testing", "every-two-weeks-testing", "none"))
  
  dat.sum$intervention_class = factor(dat.sum$intervention_class, levels=c("testing", "group-size-limit", "symptom-isolation"))
  dat.sum <- arrange(dat.sum, contact_tracing, intervention_class, test_rotation, LOD, TAT, distance_limit)
  
  
  
  write.csv(dat.sum, file = filename)
  
}
make.Table.S2(filename = "Final-Figs/SI-Appendix-Table-S2.csv")
