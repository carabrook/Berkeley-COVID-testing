rm(list=ls())


library(deSolve)
library(ggplot2)

#make ode of simple viral dynamics
#parameters pulled from the simple target cell version of the model in
#Baccam et al. 2006 J of Virology
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
    dIdt <- parms$k*E - parms$sigma*I
    dVdt <- parms$p*I - parms$c*V
    
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
one.model.run <- function(beta, sigma, k, c, p, xstart, times){
  
  #build parameter matrix
  params = list(beta=beta,
                sigma=sigma,
                c=c,
                k=k,
                p=p)
  
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


out = one.model.run(beta=1*10^-7,
                    k=1,
                    p= 100,
                    c=3.8,
                    sigma = 3.8,
                    xstart = c(Tc = 4*10^8, E=0,I=0,V=4.9), 
                    times=1:20)

#with(out, plot(time, V, type="l", xlim=c(0,5)))
with(out, plot(time, log10(V), type="l", xlim=c(0,20), ylim=c(0,12)))


#wrap model and run many trajectories
wrap.viral.loads <- function(beta.mean, sigma.mean, k.mean,
                             c.mean, p.mean, p.sd, n.pop, times, xstart){
  
  #shake the parameters with a little noise every time
  beta.list = as.list(rep(beta.mean, length=n.pop))
  beta.list = lapply(beta.list, add.noise)
  
  sigma.list = as.list(rep(sigma.mean, length=n.pop))
  sigma.list = lapply(sigma.list, add.noise)
  
  p.list = as.list(rep(p.mean, length=n.pop))
  p.list = lapply(p.list, add.noise)
  
  k.list = as.list(rep(k.mean, length=n.pop))
  k.list = lapply(k.list, add.noise)
  
  c.list = as.list(rep(c.mean, length=n.pop))
  c.list = lapply(c.list, add.noise)
  
  
  #and run over all the individuals in the pop
  out.list = mapply(one.model.run, beta=beta.list, k = k.list,
                    sigma=sigma.list, c=c.list, p=p.list, 
                    MoreArgs = list(xstart=xstart, times=times), SIMPLIFY = FALSE) 
  
  #each entry in list givens another individual
  employ.ID.list <- as.list(1:n.pop)
  
  out.list <- mapply(add.seq.ID, dat = out.list, employID = employ.ID.list, SIMPLIFY = FALSE)
  
  #and bind
  out.df <- do.call("rbind", out.list)
  head(out.df)
  
  #and attach par 
  out.df$beta <- rep(c(unlist(beta.list)), each = length(times))
  out.df$sigma <- rep(c(unlist(sigma.list)), each = length(times))
  out.df$c <- rep(c(unlist(c.list)), each = length(times))
  out.df$k <- rep(c(unlist(k.list)), each = length(times))
  out.df$p <- rep(c(unlist(p.list)), each = length(times))
  
  
  #and return
  return(out.df)
}


out <- wrap.viral.loads(beta.mean=1*10^-7,
                        k.mean=1,
                        p.mean= 100,
                        c.mean=3.8,
                        sigma.mean = 3.8,
                        xstart = c(Tc = 4*10^8, E=0,I=0,V=4.9),
                        times=1:50,
                        n.pop = 100)

head(out)

p1 <- ggplot(out) + geom_line(aes(x=time, y=log10(V), color=p, group=employIDs)) + coord_cartesian(xlim=c(0,30), ylim =c(0,12)) +
      theme_bw() + geom_hline(yintercept = 3, color="red", linetype=2) + annotate("text", label="LOD=10^3", y = 3.2, x=25) + 
      ylab("log10(viral load)") + xlab("days since exposure") + scale_color_gradient(low="gray87", high="gray48") +
      geom_hline(yintercept = 5, color="red", linetype=2) + annotate("text", label="LOD=10^5", y = 5.2, x=25) +
      geom_hline(yintercept = 7, color="red", linetype=2) + annotate("text", label="LOD=10^7", y = 7.2, x=25) +
      geom_hline(yintercept = 8.5, color="green", linetype=2) + 
      annotate("text", label="symptom onset for self-quarantine", y = 8.7, x=25) 
      
      
print(p1)

ggsave(file = "Individual_titers.png",
       units="mm",  
       width=60, 
       height=40, 
       scale=6, 
       dpi=300)





#and run these trajectories and save them

titer.dat <- wrap.viral.loads(beta.mean=1*10^-7,
                        k.mean=1,
                        p.mean= 100,
                        c.mean=3.8,
                        sigma.mean = 3.8,
                        xstart = c(Tc = 4*10^8, E=0,I=0,V=4.9),
                        times=1:50,
                        n.pop = 20000)



#save(titer.dat, file ="titer.dat.20K.Rdata")
