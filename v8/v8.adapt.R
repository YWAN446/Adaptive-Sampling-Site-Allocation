system.res.adapt <- c()
system.res.all.adapt <- c()

i=1
samples <- siteInit(iterative.ssn, n.row, n.col, n.minimum.latrine)
samples <- samples[which(!is.na(samples$ID)),]
samples$PSU <- rep(1:(length(samples$PSU)/n.sample.pooling),each=n.sample.pooling)
#samples <- siteSample(iterative.ssn, 10, n.minimum.latrine)
results <- siteAnalysis(samples, mat.conc, uf.pos)
system.res.adapt[1] <- mean(results[[2]])
system.res.all.adapt[1] <- mean(results[[3]])
psu.loc.x <- unique(samples$x)
psu.loc.y <- unique(samples$y)
source(paste0("./",version,"/",version,".analysis.R"))
for (i in 2:(n.update*n.week.per.update)){
  if ((i-1) %% n.week.per.update == 0){
    if (length(unique(samples$PSU))>n.sample){
      psu.drop <- siteDrop(results[[1]][,c(1,2,2+seq((i-1)*n.days.per.sample+1,(i+n.week.per.update-1)*n.days.per.sample,by=n.days.per.sample))],
                           uf.pos[seq((i-1)*n.days.per.sample+1,(i+n.week.per.update-1)*n.days.per.sample,by=n.days.per.sample)],n.drop=(n.drop + n.drop.add))
      samples <- samples[-which(samples$PSU %in% psu.drop),]
      samples <- siteAdd(n.add=n.drop.add, psu.drop, iterative.ssn, n.minimum.latrine, samples)
    } else {
      psu.drop <- siteDrop(results[[1]][,c(1,2,2+seq((i-1)*n.days.per.sample+1,(i+n.week.per.update-1)*n.days.per.sample,by=n.days.per.sample))],
                           uf.pos[seq((i-1)*n.days.per.sample+1,(i+n.week.per.update-1)*n.days.per.sample,by=n.days.per.sample)],n.drop=n.drop.add)
      samples <- samples[-which(samples$PSU %in% psu.drop),]
      samples <- siteAdd(n.add=n.drop.add, psu.drop, iterative.ssn, n.minimum.latrine, samples)
    }
    psu.loc.x <- c(psu.loc.x,unique(samples$x))
    psu.loc.y <- c(psu.loc.y,unique(samples$y))
    results <- siteAnalysis(samples, mat.conc, uf.pos)
    system.res.adapt[i] <- mean(results[[2]])
    system.res.all.adapt[i] <- mean(results[[3]])
    source(paste0("./",version,"/",version,".analysis.R"))
  }
}

system.res.adapt <- system.res.adapt[(1:n.update)*n.week.per.update-n.week.per.update+1]
system.res.all.adapt <- system.res.all.adapt[(1:n.update)*n.week.per.update-n.week.per.update+1]

# pdf(file=paste0("./",version,"/output/adapt_",lambda,"_",Sys.Date(),".pdf"),width=8,height=8)
# plot(1:n.update,system.res.adapt*100,col="brown1",type="l",ylim=c(0,100),
#      main="adaptive sampling results",xlab="update",ylab="sensitivity (%)")
# lines(1:n.update,system.res.all.adapt*100,col="skyblue")
# legend("topleft",c("PSU","PS+PSU"),col=c("brown1","skyblue"),lty=1,bty="n",lwd=2)
# dev.off()
# 
# #heatmap for sampling location;
# locationPSU <- data.frame(long=psu.loc.x, lat=psu.loc.y)
# latrine.df <- data.frame(long = iterative.ssn@obspoints@SSNPoints[[1]]@point.coords[,1], 
#                          lat = iterative.ssn@obspoints@SSNPoints[[1]]@point.coords[,2],
#                          #cex = std.vec.prob.pos)
#                          cex = latrine_range$outbreak.length)
# 
# library(sp)
# pdf(file=paste0("./",version,"/output/heatmap_lambda_",lambda,"_",Sys.Date(),".pdf"),width=14,height=10)
# sewage_fortify <- fortify(sewage.lines.df)
# heatmap <- ggplot(data = sewage_fortify) + geom_path(aes(x=long, y=lat, group=group)) + 
#   ylim(rev(sewage.lines.df@bbox[2,])) + theme_classic() +
#   theme(axis.line=element_blank(),axis.text.x=element_blank(),
#         axis.text.y=element_blank(),axis.ticks=element_blank(),
#         axis.title.x=element_blank(),
#         axis.title.y=element_blank(),legend.position="none",
#         panel.background=element_blank(),panel.border=element_blank(),panel.grid.major=element_blank(),
#         panel.grid.minor=element_blank(),plot.background=element_blank())
# heatmap <- heatmap + geom_point(data = latrine.df, aes(x=long, y=lat, size = cex, colour = "red", alpha = 0.6)) + scale_radius(range = c(0.5, 10))
# heatmap <- heatmap + 
#   stat_density2d(aes(x = long, y = lat, alpha = 0.55),
#                  bins = 6, geom = "polygon",
#                  data = locationPSU) +
#   scale_fill_gradient(low = "gray70", high = "black")
# heatmap
# dev.off()
