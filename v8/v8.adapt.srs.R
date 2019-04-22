system.res.adapt <- c()
system.res.all.adapt <- c()

i=1
#samples <- siteInit(iterative.ssn, n.row, n.col, n.minimum.latrine)
samples <- siteSample(iterative.ssn, 9, n.minimum.latrine)
samples <- samples[which(!is.na(samples$ID)),]
samples$PSU <- rep(1:(length(samples$PSU)/n.sample.pooling),each=n.sample.pooling)
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