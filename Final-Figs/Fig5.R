rm(list = ls())
library(ggplot2)
library(plyr)
library(dplyr)
library(png)
library(grid)


make_Figure_5<- function(filename){
img <- readPNG(source = "Final-Figs/Figure5A.png")

p5A <- ggplot() +annotation_custom(rasterGrob(img, 
                                       width = unit(1,"npc"),
                                       height = unit(1,"npc")), 
                            -Inf, Inf, -Inf, Inf) + theme(plot.margin = unit(c(.5,.5,.5,.5), "lines"))
  
print(p5A)

#Fig 2 is gains by testing at a set distance limit. Supplemental is distance limit by all types of testing
load("Final-Figs/dat_multi_10_12.Rdata")
head(dat.multi)


#unique(dat.multi$TAT)
dat.multi$LOD <- as.character(dat.multi$LOD)
dat.multi$LOD[dat.multi$LOD==10] <- "1e+01"
unique(dat.multi$LOD)
#unique(dat.multi$TAT)

#dat.multi <- subset(dat.multi, day<=50)
unique(dat.multi$test_rotation) #this is correct

dat.plot = subset(dat.multi, test_rotation=="all-no-test" | test_rotation=="two-week-all" | test_rotation=="biweekly-high-three-week-low" | test_rotation=="biweekly-all" | test_rotation=="weekly-all")

dat.plot$test_rotation = factor(dat.plot$test_rotation, levels = c("biweekly-all", "weekly-all", "biweekly-high-three-week-low",  "two-week-all", "all-no-test"))
colz = c("biweekly-all" = "indianred1", "biweekly-high-three-week-low" = "cyan2", "two-week-all" = "springgreen3", "weekly-all" = "purple",  "all-no-test" = "gray50")

p5b <- ggplot(data=dat.plot) + geom_line(aes(x=day, y=mean_cumulative, color=test_rotation)) + scale_color_manual(values=colz) + 
  scale_fill_manual(values = colz) +
  geom_ribbon(aes(x=day, ymin=lci_cumulative, ymax=uci_cumulative, fill=test_rotation), alpha=.3) +
  coord_cartesian(ylim=c(0,20000)) + 
  ylab("cumulative daily exposures") + theme_bw() +  xlab("day of epidemic") +
  theme(panel.grid  = element_blank(), axis.title = element_text(size=14), axis.text = element_text(size=12),
        legend.title = element_blank(), legend.position = c(.8,.23), strip.background = element_blank(),
        strip.text = element_text(face="bold", size=14), plot.margin = unit(c(.5,.5,1.5,.5), "lines"))
print(p5b)  



#and calculate total cases
dat.sum <- ddply(dat.plot, .(test_rotation), summarise, total_cases=max(mean_cumulative), lci_cases = max(lci_cumulative), uci_cases=max(uci_cumulative))
dat.sum
dat.sum$cases_saved = dat.sum$total_cases[dat.sum$test_rotation=="all-no-test"] - dat.sum$total_cases
dat.sum$cases_saved_uci = dat.sum$lci_cases[dat.sum$test_rotation=="all-no-test"] - dat.sum$lci_cases
dat.sum$cases_saved_lci = dat.sum$uci_cases[dat.sum$test_rotation=="all-no-test"] - dat.sum$uci_cases

dat.sum <- arrange(dat.sum, total_cases)
dat.sum = subset(dat.sum, test_rotation!="all-no-test")


dat.sum <- arrange(dat.sum, desc(cases_saved))
dat.sum$test_rotation = as.character(dat.sum$test_rotation)
dat.sum$test_rotation<- factor(dat.sum$test_rotation, levels=c("biweekly-all", "weekly-all", "biweekly-high-three-week-low", "biweekly-high-two-week-low", "two-week-all"))

dat.sum$label <- as.character(dat.sum$test_rotation)
dat.sum$label[dat.sum$label=="biweekly-high-three-week-low"] <- "biweekly-high\nthree-week-low"
dat.sum$label <- factor(dat.sum$label, levels = c("biweekly-all", "weekly-all", "biweekly-high\nthree-week-low", "biweekly-high-two-week-low", "two-week-all") )

p5c  <- ggplot(dat.sum) + xlab("") + ylab("cases saved (2 years)") +
  geom_bar(aes(x=label, y=cases_saved, fill = test_rotation), stat = "identity", 
           position=position_dodge(), show.legend = F, color="black") + scale_color_manual(values = colz) + scale_fill_manual(values=colz) +
  geom_errorbar(aes(x=label, ymin=cases_saved_lci, ymax=cases_saved_uci, color= test_rotation), width =.3, size=.5, 
                position=position_dodge(width = .9), show.legend = F) + #coord_cartesian(ylim=c(0,3500)) +
  theme_bw()+ theme(legend.position = c(.8,.9), legend.title = element_blank(), panel.grid = element_blank(),
                    axis.title.y = element_text(size=14),  axis.text.y = element_text(size=12), axis.text.x = element_text(size=12),
                    plot.margin = unit(c(.2,.5,.5,.5), "lines"))
print(p5c)  

# and plot R0


load("Final-Figs/R0_multi_group_10_12.Rdata")
head(dat.multi.R0)

dat.multi.R0 <- arrange(dat.multi.R0, desc(mean))
dat.multi.R0.plot = subset(dat.multi.R0, test_rotation=="biweekly-all" | test_rotation== "weekly-all" | test_rotation== "biweekly-high-three-week-low" | test_rotation== "two-week-all" )

dat.multi.R0.plot$test_rotation <- factor(dat.multi.R0.plot$test_rotation, levels=c("biweekly-all",  "weekly-all", "biweekly-high-three-week-low",  "two-week-all"))

p.5d <- ggplot(data=dat.multi.R0.plot) + geom_point(aes(x=test_rotation, y=mean, color=test_rotation), shape = 16, size=5,  position=position_dodge(width=.5), show.legend = F) + 
  ylab("reduction in Reff from intervention") + xlab("test rotation")  + scale_color_manual(values = colz) + 
  geom_errorbar(aes(x=test_rotation, ymin=lci, max=uci, color=test_rotation), size=.1, show.legend = F,position=position_dodge(width=.5)) +
  theme_bw() + theme(panel.grid = element_blank(), axis.text= element_text(size = 12), axis.title = element_text(size=14), strip.text = element_text(size=14),
                     strip.background = element_rect(fill="white"), axis.title.x = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank(),
                     plot.margin = unit(c(.5,.5,0,1.2), "lines")) #+ geom_label(aes(x=behavioral_int, y=mean, label=behavioral_int))
print(p.5d)


Fig5CD <- cowplot::plot_grid(p.5d, p5c, nrow=2, ncol=1, labels = c("C.", "D."), rel_heights = c(1,1.2))
print(Fig5CD)

Fig5AB <- cowplot::plot_grid(p5A, p5b, nrow=2, ncol=1, labels = c("A.", "B."), rel_heights = c(1,1.2))
print(Fig5AB)


Fig5 <- cowplot::plot_grid(Fig5AB, Fig5CD, ncol=2, nrow=1, rel_widths = c(1,1))
print(Fig5)

ggsave(file = filename,
       units="mm",  
       width=120, 
       height=65, 
       scale=3, 
       dpi=200)
}
make_Figure_5("Final-Figs/Fig5.png")
