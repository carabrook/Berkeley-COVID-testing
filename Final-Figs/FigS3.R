rm(list = ls())
library(ggplot2)
library(plyr)
library(dplyr)

#Figure 3 shows gains in testing across a set distance limit - compare frequency/TAT/LOD
#here at a higher proportion asymptomatic

setwd("/Users/caraebrook/Documents/R/R_repositories/Berkeley-Reopening/Dec-2020/all-runs/")
make_FigS3 <- function(filename){
  
load("Final-Figs/dat.all.high.asym.12.29.Rdata")
head(dat)


dat = subset(dat, day <=50)

#first, summarize by intervention
dat.sum <- ddply(subset(dat, day==50), .(intervention_class, distance_limit, test_rotation, LOD, TAT, iso_lag, prop_asym), summarize, total_cases=unique(mean_cumulative), lci_cases=unique(lci_cumulative), uci_cases=(uci_cumulative))

#and get cases saved
dat.sum$cases_saved = dat.sum$total_cases[dat.sum$intervention_class=="none" & dat.sum$prop_asym==.51] - dat.sum$total_cases 
dat.sum$cases_saved_lci = dat.sum$lci_cases[dat.sum$intervention_class=="none" & dat.sum$prop_asym==.51] - dat.sum$lci_cases 
dat.sum$cases_saved_uci = dat.sum$uci_cases[dat.sum$intervention_class=="none" & dat.sum$prop_asym==.51] - dat.sum$uci_cases 


dat.sum$TAT = factor(dat.sum$TAT, levels = c(1,2,3,4,5,10, Inf))
dat.sum$LOD = as.character(dat.sum$LOD)
dat.sum$LOD[dat.sum$LOD=="10"] <- "1e+01"

dat.sum$iso_lag = factor(dat.sum$iso_lag, levels=c(1,2,3,4,5, Inf))

dat.sum$label = NA
dat.sum$label[dat.sum$prop_asym==.32] <- "32% asymptomatic"
dat.sum$label[dat.sum$prop_asym==.51] <- "51% asymptomatic"

dat.sum$label <- factor(dat.sum$label, levels =c("32% asymptomatic", "51% asymptomatic"))

#dat.sum.testing = subset(dat.sum, intervention_class=="testing")
dat.sum.symptom = subset(dat.sum, intervention_class=="symptom-isolation")


#dat.sum.testing$test_rotation <- sapply(strsplit(dat.sum.testing$test_rotation, "-"), '[', 1)
#dat.sum.testing$test_rotation[dat.sum.testing$test_rotation=="biweekly"] <- "biweekly-testing"
#dat.sum.testing$test_rotation[dat.sum.testing$test_rotation=="weekly"] <- "weekly-testing"
#dat.sum.testing$test_rotation[dat.sum.testing$test_rotation=="two"] <- "every-two-weeks-testing"
#dat.sum.testing$test_rotation <- factor(dat.sum.testing$test_rotation, levels = c("biweekly-testing", "weekly-testing", "every-two-weeks-testing"))


#plus cases saved:

fillz <- c("1" = "tomato", "2" = "goldenrod", "3" = "springgreen3", "4"="turquoise3", "5"="royalblue", "10" = "magenta")

#dat.sum.symp$intervention_class = droplevels(dat.sum.symp$intervention_class)
pB <- ggplot(dat.sum.symptom) + ylab("cases saved") +  facet_grid(label~intervention_class)  + scale_fill_manual(values = fillz) + scale_color_manual(values = fillz) +
  geom_bar(aes(x=iso_lag, y=cases_saved, fill=iso_lag), stat = "identity", position=position_dodge(),  color="black", show.legend = F) + coord_cartesian(ylim=c(-1000,18000)) +
  geom_errorbar(aes(x=iso_lag, ymin=cases_saved_lci, ymax=cases_saved_uci, color=iso_lag), width=.3, position=position_dodge(width = .9), show.legend = F) +
  theme_bw()+ theme(legend.position = c(.9,.85), plot.margin = unit(c(.5,.5,.5,.5), "lines"), panel.grid = element_blank(), strip.background = element_rect(fill="white"), #panel.spacing = unit(3, "lines"), 
                    axis.title = element_text(size=14),   axis.text = element_text(size=12), strip.text = element_text(size=11)) + xlab("lag to isolation (days)")

print(pB)



#next plot R-effective across each intervention and compare

load("Final-Figs/R0.all.high.asym.12.29.Rdata")


dat.test.R0$iso_lag = factor(dat.test.R0$iso_lag, levels=c(1,2,3,4,5, Inf))

dat.test.R0$label = NA
dat.test.R0$label[dat.test.R0$prop_asym==.32] <- "32% asymptomatic"
dat.test.R0$label[dat.test.R0$prop_asym==.51] <- "51% asymptomatic"
dat.test.R0$label <- factor(dat.test.R0$label, levels =c("32% asymptomatic", "51% asymptomatic"))


#dat.test.R0.testing = subset(dat.test.R0, intervention_class=="testing")
dat.test.R0.symptom = subset(dat.test.R0, intervention_class=="symptom-isolation")


#now plot testing side-by-side with behavior

colorz <- c("1" = "tomato", "2" = "goldenrod", "3" = "springgreen3", "4"="turquoise3", "5"="royalblue", "10" = "magenta", "6" = "tomato", "12" = "goldenrod", "16" = "springgreen3", "20"="turquoise3", "50" ="royalblue")


p.A <- ggplot(data=dat.test.R0.symptom) + geom_point(aes(x=iso_lag, y=mean, color=iso_lag),  shape = 15, size=5, show.legend = F) +
  scale_color_manual(values=colorz) + facet_grid(label~intervention_class) + guides(color="none", shape="none") +
  ylab("reduction in Reff from intervention") + xlab("lag to isolation (days)") +
  geom_errorbar(aes(x=iso_lag, ymin=lci, max=uci,  color=iso_lag), size=.1, show.legend = F) +
  theme_bw() +  theme(panel.grid = element_blank(), strip.text = element_text(size=12), axis.title = element_text(size=14),   axis.text = element_text(size=12),
                       plot.margin = unit(c(.5,.5,.5,.5), "lines"), strip.background = element_rect(fill="white")) #+ geom_label(aes(x=behavioral_int, y=mean, label=behavioral_int))
print(p.A)

pFig3S1 <- cowplot::plot_grid(p.A, pB, nrow=1, ncol=2, labels=c("A.", "B."), hjust=-.1)

print(pFig3S1)

ggsave(file =filename,
       units="mm",  
       width=80, 
       height=65, 
       scale=3, 
       dpi=300)
}

make_FigS3("Final-Figs/FigS3.png")
