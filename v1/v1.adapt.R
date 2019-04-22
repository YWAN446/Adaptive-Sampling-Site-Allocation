system.res.adapt <- c()
system.res.all.adapt <- c()

samples <- siteSample(iterative.ssn, 10, n.minimum.latrine)
results <- siteAnalysis(samples, mat.conc, uf.pos)
system.res.adapt[1] <- mean(results[[2]])
system.res.all.adapt[1] <- mean(results[[3]])
source(paste0("./",version,"/",version,".analysis_b4.R"))
for (i in 2:n.update){
  psu.drop <- siteDrop(results[[1]],n.drop=2)
  samples <- samples[-which(samples$PSU %in% psu.drop),]
  samples <- siteAdd(n.add=2, psu.drop, iterative.ssn, n.minimum.latrine, samples)
  results <- siteAnalysis(samples, mat.conc, uf.pos)
  system.res.adapt[i] <- mean(results[[2]])
  system.res.all.adapt[i] <- mean(results[[3]])
}
source(paste0("./",version,"/",version,".analysis_af.R"))
pdf(file=paste0("./",version,"/output/adapt_",lambda,"_",Sys.Date(),".pdf"),width=8,height=8)
plot(1:n.update,system.res.adapt*100,col="red",type="l",ylim=c(0,100),
     main="adaptive sampling results",xlab="update",ylab="sensitivity (%)")
lines(1:n.update,system.res.all.adapt*100,col="blue")
dev.off()