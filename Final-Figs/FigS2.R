rm(list=ls())

library(ggplot2)
library(plyr)
library(dplyr)

#setwd("/Users/caraebrook/Documents/R/R_repositories/Berkeley-Reopening/Dec-2020/all-runs/")
#figure 2
#and plot some distributions to go with it


make_Figure_2_S1 <- function(filename){
  
dat1 = data.frame(x = seq(0, 50, length=1000),
                  y = dlnorm(x = seq(0, 50, length=1000),
                              meanlog = log(1.05),
                              sdlog =log(1.233)))




#setwd("/Users/caraebrook/Documents/R/R_repositories/Berkeley-Reopening/test-sensitivity/")
load("Final-Figs/dat.group.lim.lognorm.12.28.Rdata")
#Fig 1 is just gains by distance limit

head(dat)
group.dat = subset(dat, day<=50)
head(group.dat)
unique(group.dat$distance_limit)
group.dat$distance_limit = factor(group.dat$distance_limit,  levels= c("6", "12", "16", "20", "50", "Inf" ))# "12",
group.dat$label = paste0("group limit = ", group.dat$distance_limit)
group.dat$label[group.dat$label=="group limit = Inf"] <- "no group limit"
group.dat$label = factor(group.dat$label, levels = c("group limit = 6", "group limit = 12",   "group limit = 16",  "group limit = 20", "group limit = 50", "no group limit")) # "group limit = 12",

#and total cases
plot.sum = subset(group.dat, day==50)
plot.sum = ddply(subset(group.dat, day ==max(day)), .(label), summarize, total_cases=unique(mean_cumulative), lci_cases=unique(lci_cumulative), uci_cases =unique(uci_cumulative))
plot.sum$lci_cases[plot.sum$lci_cases<0] <- 0
plot.sum$mean_cases_saved <- plot.sum$total_cases[plot.sum$label=="no group limit"] -plot.sum$total_cases
plot.sum$lci_cases_saved <- plot.sum$lci_cases[plot.sum$label=="no group limit"] -plot.sum$lci_cases
plot.sum$uci_cases_saved <- plot.sum$uci_cases[plot.sum$label=="no group limit"] -plot.sum$uci_cases
plot.sum$lci_cases_saved[plot.sum$lci_cases_saved<0] <- 0


plot.sum$label2 <- as.character(plot.sum$label)
plot.sum$label2 <- sub(pattern= " =", replacement = "\n=", x=plot.sum$label2)
plot.sum$label2[plot.sum$label2=="no group limit"] <- "no group\n limit"

plot.sum$label2 = factor(plot.sum$label2, levels = c("group limit\n= 6", "group limit\n= 12", "group limit\n= 16", "group limit\n= 20","group limit\n= 50","no group\n limit")) #


#and alt

label.dat <- data.frame(rbind(c("group size limits", 27,.4)))
names(label.dat) <- c("label", "x", "y")
label.dat$x = as.numeric(label.dat$x)
label.dat$y = as.numeric(label.dat$y)


p2A2 <- ggplot(data=dat1) + geom_line(aes(x=x,y=y), color="firebrick1") +
  geom_ribbon(aes(x=x,ymin=0, ymax=y), alpha=.3, fill="firebrick1") + theme_bw() +  
  ylab("probability density") + xlab("R") +
  #geom_vline(aes(xintercept=4), linetype=2, size=2) +#color="firebrick1"
  theme(panel.grid = element_blank(), 
        axis.title = element_text(size=14), axis.text = element_text(size=12),
        plot.margin = unit(c(3,1,.7,1), "lines")) +
  geom_vline(aes(xintercept=6), linetype=2,  size=1.5, color="tomato") +
  geom_vline(aes(xintercept=6), linetype=2) +
  geom_vline(aes(xintercept=12), linetype=2, size=1.5, color="goldenrod") +
  geom_vline(aes(xintercept=12), linetype=2) +
  geom_vline(aes(xintercept=16), linetype=2, size=1.5, color="springgreen3") +
  geom_vline(aes(xintercept=16), linetype=2) +
  geom_vline(aes(xintercept=20), linetype=2, size=1.5, color="turquoise3") +
  geom_vline(aes(xintercept=20), linetype=2) +
  geom_vline(aes(xintercept=50), linetype=2, size=1.5, color="royalblue") +
  geom_vline(aes(xintercept=50), linetype=2) +
  geom_label(data=label.dat[1,], aes(label=label, x=x,y=y),  fill="white") +#color="firebrick1",
  geom_label(data=label.dat[2,], aes(label=label, x=x,y=y),  fill="white") #color="turquoise3",
#annotate("text", label="grouping limits", x=15,y=.4) 

print(p2A2)  

group.dat$label2 <- paste0("group size\nlimit=", group.dat$distance_limit)
group.dat$label2[group.dat$label2=="group size\nlimit=Inf"] <- "no group\nsize limit"
group.dat$label2 <- factor(group.dat$label2, levels = c("group size\nlimit=6", "group size\nlimit=12", "group size\nlimit=16", "group size\nlimit=20", "group size\nlimit=50", "no group\nsize limit"))

colorz <- c("6" = "tomato", "12" = "goldenrod", "16" = "springgreen3", "20"="turquoise3", "50"="royalblue", "Inf" = "black")

p2B2 <- ggplot(data=group.dat) + scale_color_manual(values=colorz) + scale_fill_manual(values=colorz) +
  geom_ribbon(aes(x=day,ymin=lci_exposures, ymax=uci_exposures, fill=distance_limit),  alpha=.3, show.legend = F) + #ylim(0,350) +
  geom_line(aes(x=day, y= mean_exposures), color="black", show.legend = F) + theme_bw() + coord_cartesian(ylim=c(0,600)) +
  theme(panel.grid  = element_blank(), legend.title = element_blank(), legend.text = element_blank(), 
        strip.background = element_blank(),
        axis.title.y = element_text(size=14), axis.text.y = element_text(size=12),
        axis.title.x = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank(),
        strip.text = element_text( size=14), plot.margin = unit(c(.5,.5,.2,1.5), "lines")) +# ggtitle("no TAT") +
  ylab("daily cases") + xlab("day of epidemic") + facet_grid(~label2)

print(p2B2)


p2D2 <- ggplot(data=group.dat) + geom_line(aes(x=day, y=mean_cumulative), color="black", show.legend = F) + scale_color_manual(values=colorz) + scale_fill_manual(values=colorz) +
  geom_ribbon(aes(x=day, ymin=lci_cumulative, ymax=uci_cumulative, fill=distance_limit), alpha=.3, show.legend = F) + coord_cartesian(ylim=c(0,20000)) +
  facet_grid(~label2) + ylab("cumulative cases") + theme_bw() +  xlab("day of epidemic") +
  theme(panel.grid  = element_blank(), legend.title = element_blank(), legend.position = "bottom", strip.background = element_blank(),
        strip.text = element_blank(), axis.title = element_text(size=14), axis.text = element_text(size=12),  
        plot.margin = unit(c(0,.5,.5,.5), "lines")) 

print(p2D2)  

Fig2BC <- cowplot::plot_grid(p2B2,p2D2, nrow=2,ncol=1, labels = c("B.", "C."), label_size = 18, rel_heights = c(1,.95))
print(Fig2BC)


Fig2final <- cowplot::plot_grid(p2A2, Fig2BC, nrow=1,ncol=2, labels = c("A.",""), label_size = 18, rel_widths = c(.6,1))
print(Fig2final)

ggsave(file = filename,
       units="mm",  
       width=150, 
       height=60, 
       scale=3, 
       dpi=200)
}

make_Figure_2_S1("Final-Figs/Fig2-S1.png")


