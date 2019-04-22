#analysis
samples.res <- cbind(samples,t(mat.conc[,samples$ID]))
# psu.res <- aggregate(samples.res[,5:length(samples.res[1,])], by=list(samples.res$PSU),FUN=mean)
# names(psu.res)[1] <- "PSU"
# psu.res.wo.id <- as.matrix(psu.res[,-1])
# psu.res.wo.id[which(psu.res.wo.id<LLOD.test)] <- 0
# psu.res.wo.id[which(psu.res.wo.id>=LLOD.test)] <- 1
# psu.res <- cbind(psu.res$PSU,rowMeans(psu.res.wo.id),psu.res.wo.id)

# system.res <- as.numeric(colSums(psu.res[,-c(1,2)]))
# system.res[which(system.res>=1)]=1
# 
# system.res.all <- as.numeric(colSums(rbind(psu.res[,-c(1,2)],uf.pos)))
# system.res.all[which(system.res.all>=1)]=1

#plotting
# #pdf(file=paste0("./",version,"/output/res_pooling_b4adapt",lambda,"_",Sys.Date(),".pdf"),width=14,height=10)
# #pdf(file=paste0("./",version,"/output/res_pooling_af100adapt",lambda,"_",Sys.Date(),".pdf"),width=14,height=10)
# pdf(file=paste0("./",version,"/output/res_lambda_",lambda,"_week_",i,"_",Sys.Date(),".pdf"),width=14,height=10)
# plot(sewage.lines.df, lwd = sqrt(iterative.ssn@data$addfunccol)*5, col = "black", xlab="x-coordinate", ylab="y-coordinate", 
#      main = paste("Peak number of new shedders",lambda,"per day: week ",i), asp=0.25, ylim=rev(sewage.lines.df@bbox[2,]), cex.main=2)
# points(iterative.ssn@obspoints@SSNPoints[[1]]@point.coords[,1],
#        iterative.ssn@obspoints@SSNPoints[[1]]@point.coords[,2],cex=latrine_range$length/8)
# # points(iterative.ssn@obspoints@SSNPoints[[1]]@point.coords[,1], 
# #        iterative.ssn@obspoints@SSNPoints[[1]]@point.coords[,2],cex=mat.latrine.length[(i*7-6),]/8) #only for outbreak;
# points(iterative.ssn@obspoints@SSNPoints[[1]]@point.coords[samples[,2],1], 
#        iterative.ssn@obspoints@SSNPoints[[1]]@point.coords[samples[,2],2],
#        cex=1,pch=4)
# points(0,0,cex=2,pch=5) #pumping station sample location
# sample.x <- unique(samples$x)
# sample.y <- unique(samples$y)
# rect(sample.x-percent.x*gps.points.range.x,sample.y-percent.y*gps.points.range.y,
#      sample.x+percent.x*gps.points.range.x,sample.y+percent.y*gps.points.range.y,
#      lwd = 1.5)
# textxy(0,0,paste0("Pumping Station: ",round(uf.prob.pos*100,1),"%"),cex=1)
# # textxy(unique(samples.res$x),unique(samples.res$y),paste0("  ",psu.res[,2]*100,"%  "),cex=1)
# # legend( x="bottomleft", bty='n', cex=1,
# #         legend=paste0("The sensitivity of surveillance system (pooling samples) is ", mean(system.res)*100, "%\n",
# #                                                "The sensitivity of surveillance system (pooling samples + pumping station) is ", mean(system.res.all)*100, "%"))
# dev.off()