#Simulation of adaptive ES sampling allocation in Kolkata;
#set direction
setwd("~/stat/SurveillanceSimulation/ES sampling sites update design/Simulation/ES simul Kolkata/")
#version
version <- "v3"

#packages
library(SSN)
library(sp)
library(calibrate)
library(rgeos)

#loading data
iterative.ssn <- importSSN(paste0("./simul/",version,".SimIterative2018-07-17.ssn"))
#iterative.ssn <- importSSN(paste0("./simul/v2.SimIterative2018-07-17.ssn"))

#for (lambda in c(100,150,200,250,300,400,500)){
  #loading parameters and settings
  source(paste0("./",version,"/",version,".param.R"))
  #loading functions
  source(paste0("./",version,"/",version,".func.R"))
  #generate the sewage map and lartine locations
  #source(paste0("./",version,"/",version,".simulSewage.R"))
  #generate the simulation results
  source(paste0("./",version,"/",version,".simul.outbreak.R"))
  #plotting
  #source(paste0("./",version,"/",version,".plot.R"))
  #sampling
  #source(paste0("./",version,"/",version,".sample.R"))
  #analysis
  #source(paste0("./",version,"/",version,".analysis.R"))
  #exploration
  #source(paste0("./",version,"/",version,".explor.R"))
  #adapt
  source(paste0("./",version,"/",version,".adapt.R"))
#}

