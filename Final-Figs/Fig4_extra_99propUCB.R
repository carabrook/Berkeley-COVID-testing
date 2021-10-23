rm(list = ls())
library(ggplot2)
library(plyr)
library(dplyr)

homewd="/Users/caraebrook/Documents/R/R_repositories/Berkeley-Reopening/Berkeley-COVID-testing/"
setwd(homewd)

#higher proportion of cases in UCB community (99%) - no difference from baseline 90% model

#Figure 4 shows gains in testing across the number of test days per week -- under binom R0

make_Fig4_extra <- function(filename){
load("Final-Figs/dat.test.all.10.20.higher.prop.UCB.Rdata")
head(dat.test.all)

unique(dat.test.all$TAT)
dat.test.all$LOD <- as.character(dat.test.all$LOD)
dat.test.all$LOD[dat.test.all$LOD==10] <- "1e+01"
unique(dat.test.all$LOD)
unique(dat.test.all$TAT)

dat.test.all$test_rotation = as.character(dat.test.all$test_rotation)
dat.test.all$test_rotation[dat.test.all$test_rotation=="semi-weekly"] <- "semi-weekly testing"
dat.test.all$test_rotation[dat.test.all$test_rotation=="weekly"] <- "weekly testing"
dat.test.all$test_rotation[dat.test.all$test_rotation=="two-week"] <- "every two week testing"
dat.test.all$test_rotation = factor(dat.test.all$test_rotation, levels = c("semi-weekly testing", "weekly testing", "every two week testing"))
unique(dat.test.all$test_rotation) 
unique(dat.test.all$n_test_days)



#then, go ahead and take only the first 80 timesteps
#dat = subset(dat, day <50)
dat.test.all = subset(dat.test.all, day<=50)
dat.test = subset(dat.test.all, n_test_days==2)

# test.dat = subset(dat, distance_limit=="12")
# test.dat$test_rotation = factor(test.dat$test_rotation, levels = c("every-other-day", "biweekly-all", "weekly-all", "two-week-all", "no-test"))
# 
leg.dat =rbind(cbind.data.frame(0,0, "every two week testing", "test + trace +\nsymptom-iso + group limit"),cbind.data.frame(0,0,  "every two week testing", "test + trace +\nsymptom-iso + group limit"),cbind.data.frame(0,0,  "every two week testing", "test + trace +\nsymptom-iso + group limit"), cbind.data.frame(0,0,  "every two week testing", "test + trace +\nsymptom-iso + group limit"))
 names(leg.dat) <- c("tmp1", "tmp2", "test_rotation", "intervention_class")
 leg.dat$type <- c("daily exposures", "daily testing isolations", "daily tracing isolations", "daily symptomatic isolations")
# #leg.dat$distance_limit = factor(leg.dat$distance_limit, levels= c("6", "12", "20", "50", "none" ))
leg.dat$type <- factor(leg.dat$type, levels = c("daily exposures", "daily testing isolations", "daily tracing isolations", "daily symptomatic isolations"))
leg.dat$test_rotation = factor(leg.dat$test_rotation, levels = c("biweekly testing", "weekly testing", "every two week testing"))
leg.dat$intervention_class = factor(leg.dat$intervention_class, levels = c("test + trace +\nsymptom-iso + group limit","test + trace + symptom-iso", "test + trace", "testing only"))
# 
 colz =  c("daily exposures" = "black", "daily testing isolations" ="cornflowerblue", "daily tracing isolations" = "seagreen", "daily symptomatic isolations" = "magenta")


p4c <- ggplot(data=dat.test) + scale_color_manual(values = colz) + #scale_fill_manual(values = colz) +
  geom_ribbon(aes(x=day,ymin=lci_exposures, ymax=uci_exposures),  alpha=.3) + #ylim(0,350) +
  geom_line(data=leg.dat, aes(x=tmp1, y=tmp2, color=type)) + 
  #geom_ribbon(data=leg.dat, aes(x=tmp1, ymin=tmp2, ymax=tmp2, fill=type)) +
  #geom_ribbon(aes(x=day,ymin=lci_isolations, ymax=uci_isolations), fill="red", alpha=.3) +
  geom_ribbon(aes(x=day,ymin=lci_testing_iso, ymax=uci_testing_iso), fill="cornflowerblue", alpha=.3) +
  geom_ribbon(aes(x=day,ymin=lci_tracing_iso, ymax=uci_tracing_iso), fill="seagreen", alpha=.3) +
  geom_ribbon(aes(x=day,ymin=lci_symptomatic_iso, ymax=uci_symptomatic_iso), fill="magenta", alpha=.3) +
  #geom_line(aes(x=day, y= mean_isolations), color="red") + 
  geom_line(aes(x=day, y= mean_testing_iso), color="cornflowerblue") + 
  geom_line(aes(x=day, y= mean_tracing_iso), color="seagreen") + 
  geom_line(aes(x=day, y= mean_symptomatic_iso), color="magenta") + 
  geom_line(aes(x=day, y= mean_exposures)) + theme_bw() + coord_cartesian(ylim=c(0,300)) +
  theme(panel.grid  = element_blank(), legend.title = element_blank(), legend.position = c(.85,.9), strip.background = element_blank(),
        strip.text.x = element_text(size=12), strip.text.y = element_text(size=12)) +# ggtitle("no TAT") +
  ylab("daily cases") + xlab("day of epidemic") + facet_grid(intervention_class~test_rotation)

print(p4c)


#and calculate total cases
dat.sum <- ddply(dat.test.all, .(intervention_class, test_rotation, n_test_days), summarise, total_cases=max(mean_cumulative), lci_cases = max(lci_cumulative), uci_cases=max(uci_cumulative))
dat.sum
dat.sum$sd = (dat.sum$uci_cases - dat.sum$total_cases)/1.96
dat.sum$variance = (dat.sum$sd)^2


#and add in no test at all
load("Final-Figs/dat.all.int.12.28.Rdata")
head(dat)
dat.none = subset(dat, intervention_class=="none" & day<=50)
dat.none.sum <-ddply(dat.none, .(intervention_class), summarise, total_cases=max(mean_cumulative), lci_cases = max(lci_cumulative), uci_cases=max(uci_cumulative))

dat.sum$cases_saved <-   dat.none.sum$total_cases - dat.sum$total_cases
dat.sum$cases_saved_lci <- dat.none.sum$lci_cases - dat.sum$total_cases
dat.sum$cases_saved_uci <- dat.none.sum$uci_cases - dat.sum$total_cases


dat.sum.plot = subset(dat.sum, n_test_days ==2)


newcolz = c("test + trace +\nsymptom-iso + group limit" = "purple","test + trace + symptom-iso" = "magenta", "test + trace" = "seagreen", "testing only" = "cornflowerblue")
p4b  <- ggplot(dat.sum.plot) + xlab("") + ylab("cases saved") + scale_fill_manual(values = newcolz) +  scale_color_manual(values = newcolz) +
  geom_bar(aes(x=test_rotation, y=cases_saved, fill=intervention_class), stat = "identity", position=position_dodge(), show.legend = F) +
  geom_errorbar(aes(x=test_rotation, ymin=cases_saved_lci, ymax=cases_saved_uci,color = intervention_class ), width =.3, size=.5, 
                position=position_dodge(width = .9), show.legend = F) + coord_cartesian(ylim=c(0,18000)) +
  theme_bw()+ theme(panel.grid = element_blank(),
                    axis.title.y = element_text(size=14),  axis.text.y = element_text(size=12), axis.text.x = element_text(size=12),
                    plot.margin = unit(c(1.5,.5,.5,.5), "lines"))
print(p4b)  



# and plot R0

load("Final-Figs/R0_test_all_group_10_20_higher_prop_UCB.Rdata")


dat.test.R0$intervention_class <- NA
dat.test.R0$intervention_class[dat.test.R0$iso_lag==Inf  & dat.test.R0$trace_lag==Inf & dat.test.R0$distance_limit==Inf & dat.test.R0$test_rotation=="none"] <- "none"
dat.test.R0$intervention_class[dat.test.R0$iso_lag==Inf  & dat.test.R0$trace_lag==Inf & dat.test.R0$distance_limit==Inf& dat.test.R0$test_rotation!="none"] <- "testing only"
dat.test.R0$intervention_class[dat.test.R0$iso_lag==Inf  & dat.test.R0$trace_lag==1 & dat.test.R0$distance_limit==Inf] <- "test + trace"
dat.test.R0$intervention_class[dat.test.R0$iso_lag==1  & dat.test.R0$trace_lag==1 & dat.test.R0$distance_limit==Inf] <- "test + trace + symptom-iso"
dat.test.R0$intervention_class[dat.test.R0$iso_lag==1  & dat.test.R0$trace_lag==1 & dat.test.R0$distance_limit==12] <- "test + trace + symptom-iso + group limit"

dat.test.R0$intervention_class = factor(dat.test.R0$intervention_class,levels = c("test + trace + symptom-iso + group limit","test + trace + symptom-iso", "test + trace", "testing only", "none"))


dat.test.R0.plot = subset(dat.test.R0, n_test_days==2)



newcolz2 = c("test + trace + symptom-iso + group limit" = "purple","test + trace + symptom-iso" = "magenta", "test + trace" = "seagreen", "testing only" = "cornflowerblue")

p.4a <- ggplot(data=dat.test.R0.plot) + geom_point(aes(x=test_rotation, y=mean, color=intervention_class), shape = 16, size=5,  position=position_dodge(width=.5)) + 
  ylab("reduction in Reff from intervention") + scale_color_manual(values = newcolz2) + coord_cartesian(ylim=c(0,8)) +
  geom_errorbar(aes(x=test_rotation, ymin=lci, max=uci,  color=intervention_class), size=.1, show.legend = F,position=position_dodge(width=.5)) +
  theme_bw() + theme(panel.grid = element_blank(), axis.text= element_text(size = 12), axis.title = element_text(size=14), strip.text = element_text(size=14),
                     strip.background = element_rect(fill="white"), plot.margin = unit(c(1.5,.5,.5,2), "lines"), axis.title.x = element_blank(), 
                     legend.title = element_blank(), legend.position = c(.72,.85), legend.text = element_text(size=8)) 
print(p.4a)




#and all together
Fig4ab <- cowplot::plot_grid(p.4a, p4b, labels = c("A.", "B."), nrow=2, ncol=1)

print(Fig4ab)



Fig4all <- cowplot::plot_grid(Fig4ab, p4c, ncol=2, nrow=1, labels = c("", "C."), rel_widths = c(1, 1.5))
print(Fig4all)


ggsave(file = filename,
       units="mm",  
       width=120, 
       height=80, 
       scale=3, 
       dpi=300)


}

make_Fig4_extra(filename="Final-Figs/Fig4-extra-99propUCB.png")
