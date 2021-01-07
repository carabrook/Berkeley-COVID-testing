rm(list=ls())


library(deSolve)
library(ggplot2)

#make ode of simple viral dynamics for SARS-CoV-2
#parameters pulled from the simple target cell version of the model in
#Ke et al. 2020 medRxiv 
#Kinetics of SARS-CoV-2 infection in the human upper and lower respiratory tracts and their relationship with infectiousness

#we use only nasal swab surveillance in our campus testing so we focus on equations for the URT only here
#(additionally, in Ke et al. 2020, the authors ignore transfer of virions from LRT to URT anyway and assume infections
# start in the URT and eventually migrate downward, making URT dynamics independent)
simple.target <- function(t,x,parms){
  
  Tc = x[1]
  E = x[2]
  I = x[3]
  V = x[4]
  
  
  # The with() function gives access to the named values of parms within the
  # local environment created by the function
  with(as.list(parms),{
    
    dTcdt <- -parms$beta*Tc*V
    dEdt <- parms$beta*Tc*V - parms$k*E
    dIdt <- parms$k*E - parms$delta*I
    dVdt <- parms$piT*I - parms$c*V
    
    list(c(dTcdt,dEdt, dIdt,dVdt))
  })
}
normal_fn <- function(meanpar=NULL, sdpar=NULL){
  out <- purrr::partial(rnorm,
                        mean = meanpar,
                        sd = sdpar)
  return(out)
} 
lognormal_fn <- function(meanlogpar=NULL, sdlogpar=NULL){
  out <- purrr::partial(rlnorm,
                        meanlog = meanlogpar,
                        sdlog = sdlogpar)
  return(out)
} 
one.model.run <- function(beta, delta, k, c, piT, xstart, times){
  
  #build parameter matrix
  params = list(beta=beta,
                delta=delta,
                c=c,
                k=k,
                piT=piT)
  
  out <- as.data.frame(lsoda(y = xstart, times = times, func = simple.target, parms = params))
  out.rep <- dplyr::select(out, time, V)
  return(out.rep)
}
add.seq.ID <- function(dat, employID){
  dat$employIDs <- employID
  return(dat)
}
add.noise<- function(par){
  par.new <-  rnorm(n=1, mean=par, sd=.3*par)    
  #not allowed to be negative, so if, so, give a small number:
  if(par.new<0){
    par.new = .0000001
  }
  return(par.new)
}


#here is the deterministic model run using the best fit parameters for URT infections from Ke et al. 2020:
out = one.model.run(beta=(1.9*(10^-6)),# swab/day, URT infectivity parameter (FITTED)
                    k=4,#per day; 1/ duration of time (here 6 hrs) for newly infected cells to begin producing virus (FIXED)
                    piT= 51.4,#per swab per day; composite parameter for virus production and sampling (URT) (FITTED)
                    c=10, #per day; clearance rate of free virus (URT) (FIXED)
                    delta = 1.9, #per day; death rate of infected cells (URT) (FITTED)
                    xstart = c(Tc = (4*(10^6)), #free target cells in URT (FIXED)
                               E=0,
                               I=1, #initial number of infected cells in URT (FIXED)
                               V=0), 
                    times=1:20)

#with(out, plot(time, V, type="l", xlim=c(0,20)))
with(out, plot(time, log10(V), type="l", xlim=c(0,20), ylim=c(0,12)))


#wrap model and run many trajectories
wrap.viral.loads <- function(beta.mean, delta.mean, k.mean,
                             c.mean, piT.mean, piT.sd, n.pop, times, xstart){
  
  #shake the parameters with a little noise every time
  beta.list = as.list(rep(beta.mean, length=n.pop))
  beta.list = lapply(beta.list, add.noise)
  
  delta.list = as.list(rep(delta.mean, length=n.pop))
  delta.list = lapply(delta.list, add.noise)
  
  piT.list = as.list(rep(piT.mean, length=n.pop))
  piT.list = lapply(piT.list, add.noise)
  
  k.list = as.list(rep(k.mean, length=n.pop))
  k.list = lapply(k.list, add.noise)
  
  c.list = as.list(rep(c.mean, length=n.pop))
  c.list = lapply(c.list, add.noise)
  
  
  #and run over all the individuals in the pop
  out.list = mapply(one.model.run, beta=beta.list, k = k.list,
                    delta=delta.list, c=c.list, piT=piT.list, 
                    MoreArgs = list(xstart=xstart, times=times), SIMPLIFY = FALSE) 
  
  #each entry in list givens another individual
  employ.ID.list <- as.list(1:n.pop)
  
  out.list <- mapply(add.seq.ID, dat = out.list, employID = employ.ID.list, SIMPLIFY = FALSE)
  
  #and bind
  out.df <- do.call("rbind", out.list)
  head(out.df)
  
  #and attach par 
  out.df$beta <- rep(c(unlist(beta.list)), each = length(times))
  out.df$delta <- rep(c(unlist(delta.list)), each = length(times))
  out.df$c <- rep(c(unlist(c.list)), each = length(times))
  out.df$k <- rep(c(unlist(k.list)), each = length(times))
  out.df$piT <- rep(c(unlist(piT.list)), each = length(times))
  
  
  #and return
  return(out.df)
}


out <- wrap.viral.loads(beta.mean=(1.9*(10^-6)),
                        k.mean=4,
                        piT.mean= 51.4,
                        c.mean=10,
                        delta.mean = 1.9,
                        xstart = c(Tc = (4*(10^6)), E=0,I=1,V=0),
                        times=0:50,
                        n.pop = 100)

head(out)

p1 <- ggplot(out) + geom_line(aes(x=time, y=log10(V), color=piT, group=employIDs)) + coord_cartesian(xlim=c(0,20), ylim =c(0,12)) +
      theme_bw() + geom_hline(yintercept = 3, color="red", linetype=2) + annotate("text", label="LOD=10^3", y = 3.2, x=25) + 
      ylab("log10(viral load)") + xlab("days since exposure") + scale_color_gradient(low="gray87", high="gray48") +
      geom_hline(yintercept = 5, color="red", linetype=2) + annotate("text", label="LOD=10^5", y = 5.2, x=25) +
      #geom_hline(yintercept = 7, color="red", linetype=2) + annotate("text", label="LOD=10^7", y = 7.2, x=25) +
      geom_hline(yintercept = 6, color="green", linetype=2) + 
      annotate("text", label="threshold viral load for infectiousness", y = 6.2, x=25) 
      
      
print(p1)

ggsave(file = "model-sandbox/within-host-model/Individual_Titers_Dec2020.png",
       units="mm",  
       width=60, 
       height=40, 
       scale=6, 
       dpi=300)

#and load and plot from 20000 runs:
load("model-sandbox/within-host-model/titer.dat.20K.Rdata")
head(titer.dat)


p2 <- ggplot(titer.dat) + geom_line(aes(x=time, y=log10(V), color=piT, group=employ_ids)) + coord_cartesian(xlim=c(0,20), ylim =c(0,8)) +
  theme_bw() + geom_hline(yintercept = 3, color="red", linetype=2) + annotate("text", label="LOD=10^3", y = 3.2, x=25) + 
  ylab("log10(viral load)") + xlab("days since exposure") + scale_color_gradient(low="gray87", high="gray48") +
  geom_hline(yintercept = 5, color="red", linetype=2) + annotate("text", label="LOD=10^5", y = 5.2, x=25) +
  #geom_hline(yintercept = 7, color="red", linetype=2) + annotate("text", label="LOD=10^7", y = 7.2, x=25) +
  geom_hline(yintercept = 6, color="green", linetype=2) + 
  annotate("text", label="threshold viral load for infectiousness", y = 6.2, x=25) 


print(p2)

ggsave(file = "model-sandbox/within-host-model/all_titers_20K_Dec2020.png",
       units="mm",  
       width=60, 
       height=40, 
       scale=6, 
       dpi=300)
