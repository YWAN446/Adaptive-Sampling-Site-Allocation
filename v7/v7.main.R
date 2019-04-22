#Simulation of adaptive ES sampling allocation in Kolkata;
#set direction
setwd("~/stat/SurveillanceSimulation/ES sampling sites update design/Simulation/ES simul Kolkata/")
#version
#This version compare adaptive method with fixed method.
version <- "v7"

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

#simulate which latrines the shedders will go;
#unequal length
latrine_point <- sort(c(runif(n.latrine-1,0,n.latrine),0,n.latrine))
#equal length
# latrine_point <- 0:n.latrine
latrine_range <- data.frame(start=latrine_point[-length(latrine_point)], end=latrine_point[-1])
latrine_range$length <- latrine_range$end - latrine_range$start

gps.points.range <- iterative.ssn@obspoints@SSNPoints[[1]]@points.bbox
gps.points.range.x <- gps.points.range[1,2]-gps.points.range[1,1]
gps.points.range.y <- gps.points.range[2,2]-gps.points.range[2,1]
gps.points <- iterative.ssn@obspoints@SSNPoints[[1]]@point.coords

outbreak.x <- runif(n.outbreak.area,gps.points.range[1,1],gps.points.range[1,2])
outbreak.y <- runif(n.outbreak.area,gps.points.range[2,1],gps.points.range[2,2])
outbreak.points <- which(((gps.points[,1]-outbreak.x[1])/gps.points.range.x)^2+((gps.points[,2]-outbreak.y[1])/gps.points.range.y)^2<=outbreak.range)
for (j in 1:n.outbreak.area){
  outbreak.points <- unique(c(outbreak.points,which(((gps.points[,1]-outbreak.x[j])/gps.points.range.x)^2+((gps.points[,2]-outbreak.y[j])/gps.points.range.y)^2<=outbreak.range)))
}
latrine_range$outbreak.length <- latrine_range$length
latrine_range$outbreak.length[outbreak.points] <- latrine_range$length[outbreak.points]*outbreak.factor
latrine_range$outbreak.length <- latrine_range$outbreak.length/sum(latrine_range$outbreak.length)*n.latrine #standardize to n.latrine;
latrine_range$outbreak.end <- cumsum(latrine_range$outbreak.length)
latrine_range$outbreak.start <- c(0,latrine_range$outbreak.end[-n.latrine])

vec.lambda <- c(10,20,30,40,50)
res.adaptive <- array(NA,dim=c(length(vec.lambda),n.sim,n.update,n.sample,2))
res.simul.latrine <- array(NA,dim=c(length(vec.lambda),n.latrine,3))
for (k.lam in 1:length(vec.lambda)){
  lambda <- vec.lambda[k.lam]
  #generate the simulation results
  #source(paste0("./",version,"/",version,".simul.R"))
  #source(paste0("./",version,"/",version,".simul.outbreak.R"))
  source(paste0("./",version,"/",version,".simul.sameriskzone.R"))
  #source(paste0("./",version,"/",version,".simul.MultiOutbreak.R"))
  res.simul.latrine[k.lam,,1] <- iterative.ssn@obspoints@SSNPoints[[1]]@point.coords[,1]
  res.simul.latrine[k.lam,,2] <- iterative.ssn@obspoints@SSNPoints[[1]]@point.coords[,2]
  res.simul.latrine[k.lam,,3] <- latrine_range$outbreak.length
  for (k.iter in 1:n.sim){
    progress(k.iter)
    Sys.sleep(0.01)
    #adapt
    source(paste0("./",version,"/",version,".adapt.R"))
    res.adaptive[k.lam,k.iter,,,1] <- matrix(psu.loc.x,nrow = n.update,byrow = TRUE)
    res.adaptive[k.lam,k.iter,,,2] <- matrix(psu.loc.y,nrow = n.update,byrow = TRUE)
    if (k.iter==n.sim) cat("Done!\n")
  }
}
save(res.adaptive, res.simul.latrine, file=paste0("./",version,"/output/adapt_heatmap_",Sys.Date(),".rda"))
#load(paste0("./v7/output/adapt_heatmap_2018-09-19.rda"))

for (k.lam in 1:length(vec.lambda)){
  for (k.update in 1:n.update){
    #heatmap for sampling location;
    locationPSU <- data.frame(long=as.vector(t(res.adaptive[k.lam,,k.update,,1])), lat=as.vector(t(res.adaptive[k.lam,,k.update,,2])))
    latrine.df <- data.frame(long = res.simul.latrine[k.lam,,1], 
                             lat = res.simul.latrine[k.lam,,2],
                             cex = res.simul.latrine[k.lam,,3])
    lambda <- vec.lambda[k.lam]
    
    library(sp)
    png(file=paste0("./",version,"/output/heatmap_update_",k.update,"_lambda_",lambda,"_",Sys.Date(),".png"),width=12,height=10,units = "in",res=300)
    #pdf(file=paste0("./",version,"/output/heatmap_update_",k.update,"_lambda_",lambda,"_",Sys.Date(),".pdf"),width=14,height=10)
    sewage_fortify <- fortify(sewage.lines.df)
    heatmap <- ggplot(data = sewage_fortify) + geom_path(aes(x=long, y=lat, group=group)) + 
      ylim(rev(sewage.lines.df@bbox[2,])) + theme_classic() + ggtitle(paste0("Update ",k.update)) +
      theme(axis.line=element_blank(),axis.text.x=element_blank(),
            axis.text.y=element_blank(),axis.ticks=element_blank(),
            axis.title.x=element_blank(),
            axis.title.y=element_blank(),legend.position="none",
            plot.title = element_text(size = 40, face = "bold"),
            panel.background=element_blank(),panel.border=element_blank(),panel.grid.major=element_blank(),
            panel.grid.minor=element_blank(),plot.background=element_blank())
    heatmap <- heatmap + geom_point(data = latrine.df, aes(x=long, y=lat, size = cex, colour = "red", alpha = 0.6)) + scale_radius(range = c(0.5, 10))
    heatmap <- heatmap + 
      stat_density2d(aes(x = long, y = lat, alpha = 0.55),
                     bins = 5, geom = "polygon",
                     data = locationPSU) +
      scale_fill_gradient(low = "gray70", high = "black")
    plot(heatmap)
    dev.off()
  }
}