#Simulation of adaptive ES sampling allocation in Kolkata;
#set direction
setwd("~/stat/SurveillanceSimulation/ES sampling sites update design/Simulation/ES simul Kolkata/")
#version
#This version compare adaptive method with fixed method.
version <- "v6"

#packages
library(SSN)
library(sp)
library(calibrate)
library(rgeos)
library(ggplot2)
library(doParallel)
library(doRNG)

#generate the sewage map and lartine locations
#source(paste0("./",version,"/",version,".simulSewage.R"))
#loading data
#if not generate ssn, then load it;
source(paste0("./",version,"/",version,".load.R"))
#loading parameters and settings
source(paste0("./",version,"/",version,".param.R"))
#loading functions
source(paste0("./",version,"/",version,".func.R"))

vec.lambda <- c(10,20,30,40,50)
res.adaptive <- array(NA,dim=c(length(vec.lambda),n.sim,n.update,3))
for (k.lam in 1:length(vec.lambda)){
  lambda <- vec.lambda[k.lam]
  for (k.iter in 1:n.sim){
    #generate the simulation results
    #source(paste0("./",version,"/",version,".simul.R"))
    #source(paste0("./",version,"/",version,".simul.outbreak.R"))
    source(paste0("./",version,"/",version,".simul.MultiOutbreak.R"))
    #adapt
    source(paste0("./",version,"/",version,".adapt.R"))
    #1 is PS, #2 is PSU only, #3 is PS+PSU
    res.adaptive[k.lam,k.iter,,1] <- rep(mean(uf.pos),n.update)
    res.adaptive[k.lam,k.iter,,2] <- system.res.adapt
    res.adaptive[k.lam,k.iter,,3] <- system.res.all.adapt
  }
}
save(res.adaptive,file=paste0("./",version,"/output/adapt_multiOutbreak_",Sys.Date(),".rda"))

for (k.lam in 1:length(vec.lambda)){
  lambda <- vec.lambda[k.lam]
  pdf(file=paste0("./",version,"/output/adapt_multiOutbreak_",lambda,"_",Sys.Date(),".pdf"),width=8,height=8)
  plot(NA,ylim=c(0,100),xlim=c(0,20),
       main="adaptive sampling results",xlab="update",ylab="sensitivity (%)")
  lines(1:n.update,apply(res.adaptive[k.lam,,,1],2,function(x) quantile(x,probs=0.5))*100,col="black",lty=1)
  lines(1:n.update,apply(res.adaptive[k.lam,,,1],2,function(x) quantile(x,probs=0.05))*100,col="black",lty=2)
  lines(1:n.update,apply(res.adaptive[k.lam,,,1],2,function(x) quantile(x,probs=0.95))*100,col="black",lty=2)
  lines(1:n.update,apply(res.adaptive[k.lam,,,2],2,function(x) quantile(x,probs=0.5))*100,col="brown1",lty=1)
  lines(1:n.update,apply(res.adaptive[k.lam,,,2],2,function(x) quantile(x,probs=0.05))*100,col="brown1",lty=2)
  lines(1:n.update,apply(res.adaptive[k.lam,,,2],2,function(x) quantile(x,probs=0.95))*100,col="brown1",lty=2)
  lines(1:n.update,apply(res.adaptive[k.lam,,,3],2,function(x) quantile(x,probs=0.5))*100,col="skyblue",lty=1)
  lines(1:n.update,apply(res.adaptive[k.lam,,,3],2,function(x) quantile(x,probs=0.05))*100,col="skyblue",lty=2)
  lines(1:n.update,apply(res.adaptive[k.lam,,,3],2,function(x) quantile(x,probs=0.95))*100,col="skyblue",lty=2)
  legend("topleft",c("PS","PSU","PS+PSU"),col=c("black","brown1","skyblue"),lty=1,bty="n",lwd=1)
  dev.off()
}
