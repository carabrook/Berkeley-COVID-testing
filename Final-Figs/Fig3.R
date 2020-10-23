rm(list = ls())
library(ggplot2)
library(plyr)
library(dplyr)

#Figure 3 shows gains in testing across a set distance limit - compare frequency/TAT/LOD


make_Fig_3 <- function(filename){
load("dat.all.int.10.21.Rdata")
head(dat)


dat = subset(dat, day <=50)

#first, summarize by intervention
dat.sum <- ddply(subset(dat, day==50), .(intervention_class, distance_limit, test_rotation, LOD, TAT, iso_lag), summarize, total_cases=unique(mean_cumulative), lci_cases=unique(lci_cumulative), uci_cases=(uci_cumulative))

#and get cases saved
dat.sum$cases_saved = dat.sum$total_cases[dat.sum$intervention_class=="none"] - dat.sum$total_cases 
dat.sum$cases_saved_lci = dat.sum$lci_cases[dat.sum$intervention_class=="none"] - dat.sum$lci_cases 
dat.sum$cases_saved_uci = dat.sum$uci_cases[dat.sum$intervention_class=="none"] - dat.sum$uci_cases 


dat.sum$TAT = factor(dat.sum$TAT, levels = c(1,2,3,4,5,10, Inf))
dat.sum$LOD = as.character(dat.sum$LOD)
dat.sum$LOD[dat.sum$LOD=="10"] <- "1e+01"
dat.sum$LOD[dat.sum$LOD=="1000"] <- "1e+03"
dat.sum$LOD = factor(dat.sum$LOD, levels = c("1e+01", "1e+03", "1e+05", "1e+07"))

dat.sum$iso_lag = factor(dat.sum$iso_lag, levels=c(1,2,3,4,5, Inf))

dat.sum$distance_limit = factor(dat.sum$distance_limit, levels=c(6,12,16, 20,50,Inf))

dat.sum.testing = subset(dat.sum, intervention_class=="testing")
dat.sum.symptom = subset(dat.sum, intervention_class=="symptom-isolation")
dat.sum.group = subset(dat.sum, intervention_class=="group-size-limit")

dat.sum.testing$test_rotation <- sapply(strsplit(dat.sum.testing$test_rotation, "-"), '[', 1)
dat.sum.testing$test_rotation[dat.sum.testing$test_rotation=="biweekly"] <- "biweekly-testing"
dat.sum.testing$test_rotation[dat.sum.testing$test_rotation=="weekly"] <- "weekly-testing"
dat.sum.testing$test_rotation[dat.sum.testing$test_rotation=="two"] <- "every-two-weeks-testing"
dat.sum.testing$test_rotation <- factor(dat.sum.testing$test_rotation, levels = c("biweekly-testing", "weekly-testing", "every-two-weeks-testing"))

#plus cases saved:
alphaz <- c("1e+01" = 1, "1e+03" = .7, "1e+05"=.5, "1e+07" =.3)
fillz <- c("1" = "tomato", "2" = "goldenrod", "3" = "springgreen3", "4"="turquoise3", "5"="royalblue", "10" = "magenta")
p1a2 <- ggplot(dat.sum.testing)  + ylab("cases saved") +  facet_grid(~test_rotation)  + scale_alpha_manual(values=alphaz) + #ggtitle("testing interventions")
  scale_fill_manual(values=fillz) + guides(fill="none") +  geom_bar(aes(x=TAT, y=cases_saved, alpha=LOD, fill = TAT), stat = "identity", position=position_dodge(),  color="black", show.legend = F) + coord_cartesian(ylim=c(-1000,12000)) +
  geom_errorbar(aes(x=TAT, ymin=cases_saved_lci, ymax=cases_saved_uci, group=LOD, color=TAT, alpha=LOD), width=.3, position=position_dodge(width = .9), show.legend = F)  +
  theme_bw()+ theme(legend.position = c(.9,.85), plot.margin = unit(c(.2,.2,.5,0), "lines"), panel.grid = element_blank(), strip.background = element_rect(fill="white"), axis.title.y = element_blank(), axis.ticks.y = element_blank(),
                    axis.title.x = element_text(size=12),  axis.text.y = element_blank(), axis.text.x = element_text(size=12), strip.text = element_blank()) + xlab("testing mean turnaround time (days)")

print(p1a2)

#dat.sum.symp$intervention_class = droplevels(dat.sum.symp$intervention_class)
p1b2 <- ggplot(dat.sum.symptom) + ylab("cases saved") +  facet_grid(~intervention_class)  + scale_fill_manual(values = fillz) + scale_color_manual(values = fillz) +
  geom_bar(aes(x=iso_lag, y=cases_saved, fill=iso_lag), stat = "identity", position=position_dodge(),  color="black", show.legend = F) + coord_cartesian(ylim=c(-1000,12000)) +
  geom_errorbar(aes(x=iso_lag, ymin=cases_saved_lci, ymax=cases_saved_uci, color=iso_lag), width=.3, position=position_dodge(width = .9), show.legend = F) +
  theme_bw()+ theme(legend.position = c(.8,.8), panel.grid = element_blank(), strip.background = element_rect(fill="white"), # panel.spacing = unit(3, "lines"),
                    axis.text.y = element_blank(), axis.ticks.y = element_blank(), axis.title.y= element_blank(), strip.text = element_blank(),
                    axis.text.x = element_text(size=12), axis.title.x = element_text(size=12), plot.margin = unit(c(.2,.3,.5,.1), "lines")) + xlab("lag to isolation (days)")

print(p1b2)

dat.sum.group$distance_limit = factor(dat.sum.group$distance_limit, levels =c(6,12,16,20,50))

fillz2 <- c("6" = "tomato", "12" = "goldenrod", "16" = "springgreen3", "20"="turquoise3", "50"="royalblue")
p1c2 <- ggplot(dat.sum.group) + xlab("group size limit (# persons)") + ylab("cases saved") +  facet_grid(~intervention_class) + scale_fill_manual(values = fillz2) + scale_color_manual(values = fillz2) +
  geom_bar(aes(x=distance_limit, y=cases_saved, fill=distance_limit), stat = "identity", position=position_dodge(),  color="black", show.legend = F) + coord_cartesian(ylim=c(-1000,12000)) +
  geom_errorbar(aes(x=distance_limit, ymin=cases_saved_lci, ymax=cases_saved_uci, color=distance_limit), width=.3, position=position_dodge(width = .9),  show.legend = F) +
  theme_bw()+ theme(legend.position = c(.9,.85), plot.margin = unit(c(.2,.1,.5,.3), "lines"), panel.grid = element_blank(), strip.background = element_rect(fill="white"), #panel.spacing = unit(3, "lines"), 
                    axis.title = element_text(size=12),  axis.text.y = element_text(size=12), axis.text.x = element_text(size=12), strip.text = element_blank()) 

print(p1c2)



pFig3cases = cowplot::plot_grid(p1c2, p1b2,p1a2,  nrow=1,ncol=3, rel_widths =  c(1.3,1,3))

print(pFig3cases)


#next plot R-effective across each intervention and compare

load("R0.all.int.10.21.Rdata")

dat.test.R0$TAT = factor(dat.test.R0$TAT, levels = c(1,2,3,4,5,10, Inf))
dat.test.R0$LOD = as.character(dat.test.R0$LOD)
dat.test.R0$LOD = factor(dat.test.R0$LOD, levels = c("1e+01", "1e+03", "1e+05", "1e+07"))

dat.test.R0$iso_lag = factor(dat.test.R0$iso_lag, levels=c(1,2,3,4,5, Inf))

dat.test.R0$distance_limit = factor(dat.test.R0$distance_limit, levels=c(6,12,16, 20,50,Inf))

dat.test.R0.testing = subset(dat.test.R0, intervention_class=="testing")
dat.test.R0.symptom = subset(dat.test.R0, intervention_class=="symptom-isolation")
dat.test.R0.group = subset(dat.test.R0, intervention_class=="group-size-limit")

dat.test.R0.testing$test_rotation <- sapply(strsplit(dat.test.R0.testing$test_rotation, "-"), '[', 1)
dat.test.R0.testing$test_rotation[dat.test.R0.testing$test_rotation=="biweekly"] <- "biweekly-testing"
dat.test.R0.testing$test_rotation[dat.test.R0.testing$test_rotation=="weekly"] <- "weekly-testing"
dat.test.R0.testing$test_rotation[dat.test.R0.testing$test_rotation=="two"] <- "every-two-weeks-testing"
dat.test.R0.testing$test_rotation <- factor(dat.test.R0.testing$test_rotation, levels = c("biweekly-testing", "weekly-testing", "every-two-weeks-testing"))


#now plot testing side-by-side with behavior

colorz <- c("1" = "tomato", "2" = "goldenrod", "3" = "springgreen3", "4"="turquoise3", "5"="royalblue", "10" = "magenta", "6" = "tomato", "12" = "goldenrod", "16" = "springgreen3", "20"="turquoise3", "50" ="royalblue")
p.2a2 <- ggplot(data=dat.test.R0.testing) + geom_point(aes(x=TAT, y=mean, color=TAT, alpha=LOD, group=LOD), shape = 16, size=5,show.legend = F,  position=position_dodge(width=.5)) + scale_color_manual(values=colorz) + scale_alpha_manual(values=alphaz) + facet_grid(~test_rotation) + #guides(color="none", shape="none") +
  ylab("reduction in Reff from intervention") + xlab("testing mean turnaround time (days)") +
  geom_errorbar(aes(x=TAT, ymin=lci, max=uci,  color=TAT, alpha=LOD), size=.1, show.legend = F,position=position_dodge(width=.5)) +
  theme_bw() + theme(panel.grid = element_blank(), strip.text = element_text(size=12), axis.text = element_blank(),axis.ticks = element_blank(),axis.title = element_blank(),
                     strip.background = element_rect(fill="white"), plot.margin = unit(c(1.5,.2,0,0), "lines")) #+ geom_label(aes(x=behavioral_int, y=mean, label=behavioral_int))
print(p.2a2)

dat.test.R0.group$split = "group-size-limit"

p.2b2 <- ggplot(data=dat.test.R0.group) + geom_point(aes(x=distance_limit, y=mean, color=distance_limit),  shape = 15, size=5, show.legend = F) + scale_color_manual(values=colorz) + scale_alpha_manual(values=alphaz) + facet_grid(~split) + guides(color="none", shape="none") +
  ylab("reduction in Reff from intervention") +xlab("group size limit (# persons)") +
  geom_errorbar(aes(x=distance_limit, ymin=lci, max=uci,  color=distance_limit, alpha=LOD), size=.1, show.legend = F) +
  theme_bw() +  theme(panel.grid = element_blank(), axis.ticks=element_blank(), axis.text.x = element_blank(), axis.ticks.x=element_blank(),  axis.title.x = element_blank(),  strip.text = element_text(size=12),
                      strip.background = element_rect(fill="white"), plot.margin = unit(c(1.5,0,0,2.1), "lines"),  axis.title.y = element_text(size = 14), axis.text.y = element_text(size=12)) #+ geom_label(aes(x=behavioral_int, y=mean, label=behavioral_int))
print(p.2b2)

dat.test.R0.symptom$split = "symptom-isolation"

p.2c2 <- ggplot(data=dat.test.R0.symptom) + geom_point(aes(x=iso_lag, y=mean, color=iso_lag),  shape = 15, size=5, show.legend = F) + scale_color_manual(values=colorz) + scale_alpha_manual(values=alphaz) + facet_grid(~split) + guides(color="none", shape="none") +
  ylab("reduction in Reff from intervention") + xlab("lag to isolation (days)") +
  geom_errorbar(aes(x=iso_lag, ymin=lci, max=uci,  color=iso_lag, alpha=LOD), size=.1, show.legend = F) +
  theme_bw() +  theme(panel.grid = element_blank(), strip.text = element_text(size=12), axis.text = element_blank(),axis.ticks = element_blank(),axis.title = element_blank(),
                       plot.margin = unit(c(1.5,.25,0,.2), "lines"), strip.background = element_rect(fill="white")) #+ geom_label(aes(x=behavioral_int, y=mean, label=behavioral_int))
print(p.2c2)



pFig3justR0 <- cowplot::plot_grid( p.2b2, p.2c2,p.2a2,ncol=3, nrow=1, rel_widths = c(1.3,1,3))
print(pFig3justR0)


pFig3alt <- cowplot::plot_grid(pFig3justR0, pFig3cases, nrow=2, ncol=1, labels=c("A.", "B."), hjust=-.1)

print(pFig3alt)

ggsave(file =filename,
       units="mm",  
       width=100, 
       height=70, 
       scale=3, 
       dpi=300)
}

make_Fig_3("Fig3.png")


#and save the raw data

make.Table.S1 <- function(filename){
  
load("dat.all.int.10.21.Rdata")
head(dat)
  
  
dat = subset(dat, day <=50)
  
#first, summarize by intervention
dat.sum <- ddply(subset(dat, day==50), .(intervention_class, distance_limit, test_rotation, LOD, TAT, iso_lag), summarize, total_cases=unique(mean_cumulative), lci_cases=unique(lci_cumulative), uci_cases=(uci_cumulative))
  
#and get cases saved
dat.sum$cases_saved = dat.sum$total_cases[dat.sum$intervention_class=="none"] - dat.sum$total_cases 
dat.sum$cases_saved_lci = dat.sum$lci_cases[dat.sum$intervention_class=="none"] - dat.sum$lci_cases 
dat.sum$cases_saved_uci = dat.sum$uci_cases[dat.sum$intervention_class=="none"] - dat.sum$uci_cases 
  
  
dat.sum$TAT = factor(dat.sum$TAT, levels = c(1,2,3,4,5,10, Inf))
dat.sum$LOD = as.character(dat.sum$LOD)
dat.sum$LOD[dat.sum$LOD=="10"] <- "1e+01"
dat.sum$LOD[dat.sum$LOD=="1000"] <- "1e+03"
dat.sum$LOD = factor(dat.sum$LOD, levels = c("1e+01", "1e+03", "1e+05", "1e+07"))
  
dat.sum$iso_lag = factor(dat.sum$iso_lag, levels=c(1,2,3,4,5, Inf))
  
dat.sum$distance_limit = factor(dat.sum$distance_limit, levels=c(6,12,16, 20,50,Inf))


load("R0.all.int.10.21.Rdata")

dat.test.R0$TAT = factor(dat.test.R0$TAT, levels = c(1,2,3,4,5,10, Inf))
dat.test.R0$LOD = as.character(dat.test.R0$LOD)
dat.test.R0$LOD = factor(dat.test.R0$LOD, levels = c("1e+01", "1e+03", "1e+05", "1e+07"))

dat.test.R0$iso_lag = factor(dat.test.R0$iso_lag, levels=c(1,2,3,4,5, Inf))

dat.test.R0$distance_limit = factor(dat.test.R0$distance_limit, levels=c(6,12,16, 20,50,Inf))


#and save the raw data:
head(dat.sum)
head(dat.test.R0)



#and add in Reff
dat.sum$uci_Reff_reduction  <- dat.sum$lci_Reff_reduction <- dat.sum$mean_Reff_reduction <- NA

dat.sum$intervention = dat.sum$test_rotation
dat.sum$intervention[dat.sum$intervention_class=="group-size-limit"] <- paste0("group-lim=",as.character(dat.sum$distance_limit[dat.sum$intervention_class=="group-size-limit"]))
dat.sum$intervention[dat.sum$intervention_class=="symptom-isolation"] <- paste0("symptom-iso=",as.character(dat.sum$iso_lag[dat.sum$intervention_class=="symptom-isolation"]))


dat.test.R0$intervention = dat.test.R0$test_rotation
dat.test.R0$intervention[dat.test.R0$intervention_class=="group-size-limit"] <- paste0("group-lim=",as.character(dat.test.R0$distance_limit[dat.test.R0$intervention_class=="group-size-limit"]))
dat.test.R0$intervention[dat.test.R0$intervention_class=="symptom-isolation"] <- paste0("symptom-iso=",as.character(dat.test.R0$iso_lag[dat.test.R0$intervention_class=="symptom-isolation"]))



for (i in 1:length(dat.test.R0$mean)){
  dat.sum$mean_Reff_reduction[dat.sum$intervention==dat.test.R0$intervention[i]] <- dat.test.R0$mean[i]
  dat.sum$lci_Reff_reduction[dat.sum$intervention==dat.test.R0$intervention[i]] <- dat.test.R0$lci[i]
  dat.sum$uci_Reff_reduction[dat.sum$intervention==dat.test.R0$intervention[i]] <- dat.test.R0$uci[i]
}


dat.sum$test_rotation <- sapply(strsplit(dat.sum$test_rotation, "-"), '[', 1)
dat.sum$test_rotation[dat.sum$test_rotation=="biweekly"] <- "biweekly-testing"
dat.sum$test_rotation[dat.sum$test_rotation=="weekly"] <- "weekly-testing"
dat.sum$test_rotation[dat.sum$test_rotation=="two"] <- "every-two-weeks-testing"
dat.sum$test_rotation <- factor(dat.sum$test_rotation, levels = c("biweekly-testing", "weekly-testing", "every-two-weeks-testing", "none"))

dat.sum$intervention_class = factor(dat.sum$intervention_class, levels=c("testing", "group-size-limit", "symptom-isolation"))
dat.sum <- arrange(dat.sum, intervention_class, test_rotation, LOD, TAT, distance_limit)



write.csv(dat.sum, file = filename)

}
make.Table.S1(filename = "SI-Appendix-Table-S1.csv")
