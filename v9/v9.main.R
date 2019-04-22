#Simulation of adaptive ES sampling allocation in Kolkata;
#set direction
setwd("~/stat/SurveillanceSimulation/ES sampling sites update design/Simulation/ES simul Kolkata/")
#version
#This version compare adaptive method with fixed method.
version <- "v9"

#packages
library(SSN)
library(sp)
library(calibrate)
library(rgeos)
library(ggplot2)
library(doParallel)
library(doRNG)
library(svMisc)

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
#res.adaptive <- array(NA,dim=c(length(vec.lambda),n.sim,n.update,3,6))

###pick up###
  k.lam = 5
  lambda <- vec.lambda[k.lam]
  for (k.iter in 710:n.sim){
    progress(k.iter, n.sim)
    Sys.sleep(0.01)
    
    source(paste0("./",version,"/",version,".param.R"))
    #generate the simulation results
    source(paste0("./",version,"/",version,".simul.MultiOutbreak.R"))
    #adapt
    source(paste0("./",version,"/",version,".adapt.R"))
    #1 is PS, #2 is PSU only, #3 is PS+PSU
    res.adaptive[k.lam,k.iter,,1,1] <- rep(mean(uf.pos),n.update)
    res.adaptive[k.lam,k.iter,,2,1] <- system.res.adapt
    res.adaptive[k.lam,k.iter,,3,1] <- system.res.all.adapt
    
    #parameter changes: double sample size initialization;
    n.row <- 4
    n.col <- 5
    #adapt
    source(paste0("./",version,"/",version,".adapt.R"))
    #1 is PS, #2 is PSU only, #3 is PS+PSU
    res.adaptive[k.lam,k.iter,,1,2] <- rep(mean(uf.pos),n.update)
    res.adaptive[k.lam,k.iter,,2,2] <- system.res.adapt
    res.adaptive[k.lam,k.iter,,3,2] <- system.res.all.adapt
    
    #parameter changes: 3 samples per day;
    n.days.per.sample <- 3
    #adapt
    source(paste0("./",version,"/",version,".adapt.R"))
    #1 is PS, #2 is PSU only, #3 is PS+PSU
    res.adaptive[k.lam,k.iter,,1,3] <- rep(mean(uf.pos),n.update)
    res.adaptive[k.lam,k.iter,,2,3] <- system.res.adapt
    res.adaptive[k.lam,k.iter,,3,3] <- system.res.all.adapt
    
    #parameter changes: 8 weeks per update;
    n.week.per.update <- 8
    #adapt
    source(paste0("./",version,"/",version,".adapt.R"))
    #1 is PS, #2 is PSU only, #3 is PS+PSU
    res.adaptive[k.lam,k.iter,,1,4] <- rep(mean(uf.pos),n.update)
    res.adaptive[k.lam,k.iter,,2,4] <- system.res.adapt
    res.adaptive[k.lam,k.iter,,3,4] <- system.res.all.adapt
    
    #parameter changes: update 1 site;
    n.drop.add <- 1
    #adapt
    source(paste0("./",version,"/",version,".adapt.R"))
    #1 is PS, #2 is PSU only, #3 is PS+PSU
    res.adaptive[k.lam,k.iter,,1,5] <- rep(mean(uf.pos),n.update)
    res.adaptive[k.lam,k.iter,,2,5] <- system.res.adapt
    res.adaptive[k.lam,k.iter,,3,5] <- system.res.all.adapt
    
    #adapt
    source(paste0("./",version,"/",version,".adapt.srs.R"))
    #1 is PS, #2 is PSU only, #3 is PS+PSU
    res.adaptive[k.lam,k.iter,,1,6] <- rep(mean(uf.pos),n.update)
    res.adaptive[k.lam,k.iter,,2,6] <- system.res.adapt
    res.adaptive[k.lam,k.iter,,3,6] <- system.res.all.adapt
    
    if (k.iter==n.sim) cat("Part Done!\n")
  }

########################################## four scenarios ##########################################
# for (k.lam in 4:length(vec.lambda)){
#   k.lam=5
#   lambda <- vec.lambda[k.lam]
#   for (k.iter in 1:n.sim){
#     progress(k.iter, n.sim)
#     Sys.sleep(0.01)
#     
#     source(paste0("./",version,"/",version,".param.R"))
#     #generate the simulation results
#     source(paste0("./",version,"/",version,".simul.MultiOutbreak.R"))
#     #adapt
#     source(paste0("./",version,"/",version,".adapt.R"))
#     #1 is PS, #2 is PSU only, #3 is PS+PSU
#     res.adaptive[k.lam,k.iter,,1,1] <- rep(mean(uf.pos),n.update)
#     res.adaptive[k.lam,k.iter,,2,1] <- system.res.adapt
#     res.adaptive[k.lam,k.iter,,3,1] <- system.res.all.adapt
#     
#     #parameter changes: double sample size initialization;
#     n.row <- 4
#     n.col <- 5
#     #adapt
#     source(paste0("./",version,"/",version,".adapt.R"))
#     #1 is PS, #2 is PSU only, #3 is PS+PSU
#     res.adaptive[k.lam,k.iter,,1,2] <- rep(mean(uf.pos),n.update)
#     res.adaptive[k.lam,k.iter,,2,2] <- system.res.adapt
#     res.adaptive[k.lam,k.iter,,3,2] <- system.res.all.adapt
#     
#     #parameter changes: 3 samples per day;
#     n.days.per.sample <- 3
#     #adapt
#     source(paste0("./",version,"/",version,".adapt.R"))
#     #1 is PS, #2 is PSU only, #3 is PS+PSU
#     res.adaptive[k.lam,k.iter,,1,3] <- rep(mean(uf.pos),n.update)
#     res.adaptive[k.lam,k.iter,,2,3] <- system.res.adapt
#     res.adaptive[k.lam,k.iter,,3,3] <- system.res.all.adapt
#     
#     #parameter changes: 8 weeks per update;
#     n.week.per.update <- 8
#     #adapt
#     source(paste0("./",version,"/",version,".adapt.R"))
#     #1 is PS, #2 is PSU only, #3 is PS+PSU
#     res.adaptive[k.lam,k.iter,,1,4] <- rep(mean(uf.pos),n.update)
#     res.adaptive[k.lam,k.iter,,2,4] <- system.res.adapt
#     res.adaptive[k.lam,k.iter,,3,4] <- system.res.all.adapt
#     
#     #parameter changes: update 1 site;
#     n.drop.add <- 1
#     #adapt
#     source(paste0("./",version,"/",version,".adapt.R"))
#     #1 is PS, #2 is PSU only, #3 is PS+PSU
#     res.adaptive[k.lam,k.iter,,1,5] <- rep(mean(uf.pos),n.update)
#     res.adaptive[k.lam,k.iter,,2,5] <- system.res.adapt
#     res.adaptive[k.lam,k.iter,,3,5] <- system.res.all.adapt
#     
#     #adapt
#     source(paste0("./",version,"/",version,".adapt.srs.R"))
#     #1 is PS, #2 is PSU only, #3 is PS+PSU
#     res.adaptive[k.lam,k.iter,,1,6] <- rep(mean(uf.pos),n.update)
#     res.adaptive[k.lam,k.iter,,2,6] <- system.res.adapt
#     res.adaptive[k.lam,k.iter,,3,6] <- system.res.all.adapt
#     
#     if (k.iter==n.sim) cat("Part Done!\n")
#   }
# }
save(res.adaptive,file=paste0("./",version,"/output/adapt_settings",Sys.Date(),".rda"))

# for (k.lam in 1:length(vec.lambda)){
#   lambda <- vec.lambda[k.lam]
#   pdf(file=paste0("./",version,"/output/adapt_doubleInitial_",lambda,"_",Sys.Date(),".pdf"),width=8,height=8)
#   plot(NA,ylim=c(0,100),xlim=c(0,20),main="adaptive sampling results (PS+PSU)",xlab="update",ylab="sensitivity (%)")
#   lines(1:n.update,apply(res.adaptive[k.lam,,,3,1],2,function(x) quantile(x,probs=0.5))*100,col="skyblue",lty=1)
#   lines(1:n.update,apply(res.adaptive[k.lam,,,3,1],2,function(x) quantile(x,probs=0.05))*100,col="skyblue",lty=2)
#   lines(1:n.update,apply(res.adaptive[k.lam,,,3,1],2,function(x) quantile(x,probs=0.95))*100,col="skyblue",lty=2)
#   lines(1:n.update,apply(res.adaptive[k.lam,,,3,2],2,function(x) quantile(x,probs=0.5))*100,col="brown1",lty=1)
#   lines(1:n.update,apply(res.adaptive[k.lam,,,3,2],2,function(x) quantile(x,probs=0.05))*100,col="brown1",lty=2)
#   lines(1:n.update,apply(res.adaptive[k.lam,,,3,2],2,function(x) quantile(x,probs=0.95))*100,col="brown1",lty=2)
#   legend("bottomright",c("Normal size initialization","Double size initialization"),col=c("skyblue","brown1"),lty=1,bty="n",lwd=1)
#   dev.off()
# }

# for (k.lam in 1:length(vec.lambda)){
#   lambda <- vec.lambda[k.lam]
#   pdf(file=paste0("./",version,"/output/adapt_3daysPerSample_",lambda,"_",Sys.Date(),".pdf"),width=8,height=8)
#   plot(NA,ylim=c(0,100),xlim=c(0,20),main="adaptive sampling results (PS+PSU)",xlab="update",ylab="sensitivity (%)")
#   lines(1:n.update,apply(res.adaptive[k.lam,,,3,1],2,function(x) quantile(x,probs=0.5))*100,col="skyblue",lty=1)
#   lines(1:n.update,apply(res.adaptive[k.lam,,,3,1],2,function(x) quantile(x,probs=0.05))*100,col="skyblue",lty=2)
#   lines(1:n.update,apply(res.adaptive[k.lam,,,3,1],2,function(x) quantile(x,probs=0.95))*100,col="skyblue",lty=2)
#   lines(1:n.update,apply(res.adaptive[k.lam,,,3,2],2,function(x) quantile(x,probs=0.5))*100,col="brown1",lty=1)
#   lines(1:n.update,apply(res.adaptive[k.lam,,,3,2],2,function(x) quantile(x,probs=0.05))*100,col="brown1",lty=2)
#   lines(1:n.update,apply(res.adaptive[k.lam,,,3,2],2,function(x) quantile(x,probs=0.95))*100,col="brown1",lty=2)
#   legend("bottomright",c("7 days per sample","3 days per sample"),col=c("skyblue","brown1"),lty=1,bty="n",lwd=1)
#   dev.off()
# }

# for (k.lam in 1:length(vec.lambda)){
#   lambda <- vec.lambda[k.lam]
#   pdf(file=paste0("./",version,"/output/adapt_8SamplePerUpdate_",lambda,"_",Sys.Date(),".pdf"),width=8,height=8)
#   plot(NA,ylim=c(0,100),xlim=c(0,20),main="adaptive sampling results (PS+PSU)",xlab="update",ylab="sensitivity (%)")
#   lines(1:n.update,apply(res.adaptive[k.lam,,,3,1],2,function(x) quantile(x,probs=0.5))*100,col="skyblue",lty=1)
#   lines(1:n.update,apply(res.adaptive[k.lam,,,3,1],2,function(x) quantile(x,probs=0.05))*100,col="skyblue",lty=2)
#   lines(1:n.update,apply(res.adaptive[k.lam,,,3,1],2,function(x) quantile(x,probs=0.95))*100,col="skyblue",lty=2)
#   lines(1:n.update,apply(res.adaptive[k.lam,,,3,2],2,function(x) quantile(x,probs=0.5))*100,col="brown1",lty=1)
#   lines(1:n.update,apply(res.adaptive[k.lam,,,3,2],2,function(x) quantile(x,probs=0.05))*100,col="brown1",lty=2)
#   lines(1:n.update,apply(res.adaptive[k.lam,,,3,2],2,function(x) quantile(x,probs=0.95))*100,col="brown1",lty=2)
#   legend("bottomright",c("12 samples per update","8 samples per update"),col=c("skyblue","brown1"),lty=1,bty="n",lwd=1)
#   dev.off()
# }

# for (k.lam in 1:length(vec.lambda)){
#   lambda <- vec.lambda[k.lam]
#   pdf(file=paste0("./",version,"/output/adapt_1SitePerUpdate_",lambda,"_",Sys.Date(),".pdf"),width=8,height=8)
#   plot(NA,ylim=c(0,100),xlim=c(0,20),main="adaptive sampling results (PS+PSU)",xlab="update",ylab="sensitivity (%)")
#   lines(1:n.update,apply(res.adaptive[k.lam,,,3,1],2,function(x) quantile(x,probs=0.5))*100,col="skyblue",lty=1)
#   lines(1:n.update,apply(res.adaptive[k.lam,,,3,1],2,function(x) quantile(x,probs=0.05))*100,col="skyblue",lty=2)
#   lines(1:n.update,apply(res.adaptive[k.lam,,,3,1],2,function(x) quantile(x,probs=0.95))*100,col="skyblue",lty=2)
#   lines(1:n.update,apply(res.adaptive[k.lam,,,3,2],2,function(x) quantile(x,probs=0.5))*100,col="brown1",lty=1)
#   lines(1:n.update,apply(res.adaptive[k.lam,,,3,2],2,function(x) quantile(x,probs=0.05))*100,col="brown1",lty=2)
#   lines(1:n.update,apply(res.adaptive[k.lam,,,3,2],2,function(x) quantile(x,probs=0.95))*100,col="brown1",lty=2)
#   legend("bottomright",c("relocate 2 sites per update","relocate 1 site per update"),col=c("skyblue","brown1"),lty=1,bty="n",lwd=1)
#   dev.off()
# }
##############################################################################################################################

########################################## stratified sampling ##########################################;
# for (k.lam in 1:length(vec.lambda)){
#   source(paste0("./",version,"/",version,".param.R"))
#   lambda <- vec.lambda[k.lam]
#   for (k.iter1 in 1:n.sim){
#     progress(k.iter1, n.sim)
#     Sys.sleep(0.01)
#     #generate the simulation results
#     source(paste0("./",version,"/",version,".simul.MultiOutbreak.R"))
#     #adapt
#     source(paste0("./",version,"/",version,".adapt.R"))
#     #1 is PS, #2 is PSU only, #3 is PS+PSU
#     res.adaptive[k.lam,k.iter1,,1,1] <- rep(mean(uf.pos),n.update)
#     res.adaptive[k.lam,k.iter1,,2,1] <- system.res.adapt
#     res.adaptive[k.lam,k.iter1,,3,1] <- system.res.all.adapt
#     if (k.iter1==n.sim) cat("Part 1 Done!\n")
#   }
#   for (k.iter2 in 1:n.sim){
#     progress(k.iter2, n.sim)
#     Sys.sleep(0.01)
#     #generate the simulation results
#     source(paste0("./",version,"/",version,".simul.MultiOutbreak.R"))
#     #adapt
#     source(paste0("./",version,"/",version,".adapt.srs.R"))
#     #1 is PS, #2 is PSU only, #3 is PS+PSU
#     res.adaptive[k.lam,k.iter2,,1,2] <- rep(mean(uf.pos),n.update)
#     res.adaptive[k.lam,k.iter2,,2,2] <- system.res.adapt
#     res.adaptive[k.lam,k.iter2,,3,2] <- system.res.all.adapt
#     if (k.iter2==n.sim) cat("Part 2 Done!\n")
#   }
# }
# save(res.adaptive,file=paste0("./",version,"/output/adapt_stratifiedSampling_",Sys.Date(),".rda"))
# 
# for (k.lam in 1:length(vec.lambda)){
#   lambda <- vec.lambda[k.lam]
#   pdf(file=paste0("./",version,"/output/adapt_stratifiedSampling_",lambda,"_",Sys.Date(),".pdf"),width=8,height=8)
#   plot(NA,ylim=c(0,100),xlim=c(0,20),main="adaptive sampling results (PS+PSU)",xlab="update",ylab="sensitivity (%)")
#   lines(1:n.update,apply(res.adaptive[k.lam,,,3,1],2,function(x) quantile(x,probs=0.5))*100,col="skyblue",lty=1)
#   lines(1:n.update,apply(res.adaptive[k.lam,,,3,1],2,function(x) quantile(x,probs=0.05))*100,col="skyblue",lty=2)
#   lines(1:n.update,apply(res.adaptive[k.lam,,,3,1],2,function(x) quantile(x,probs=0.95))*100,col="skyblue",lty=2)
#   lines(1:n.update,apply(res.adaptive[k.lam,,,3,2],2,function(x) quantile(x,probs=0.5))*100,col="brown1",lty=1)
#   lines(1:n.update,apply(res.adaptive[k.lam,,,3,2],2,function(x) quantile(x,probs=0.05))*100,col="brown1",lty=2)
#   lines(1:n.update,apply(res.adaptive[k.lam,,,3,2],2,function(x) quantile(x,probs=0.95))*100,col="brown1",lty=2)
#   legend("bottomright",c("stratified sampling","simple random sampling"),col=c("skyblue","brown1"),lty=1,bty="n",lwd=1)
#   dev.off()
# }