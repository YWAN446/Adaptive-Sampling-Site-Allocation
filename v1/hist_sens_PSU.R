res.test <- c()
for (k.iter in 1:300){
  #sampling
  gps.points.range <- iterative.ssn@obspoints@SSNPoints[[1]]@points.bbox
  gps.points.range.x <- gps.points.range[1,2]-gps.points.range[1,1]
  gps.points.range.y <- gps.points.range[2,2]-gps.points.range[2,1]
  gps.points <- iterative.ssn@obspoints@SSNPoints[[1]]@point.coords
  
  samples <- c()
  i=1
  sample.x <- c()
  sample.y <- c()
  while (i <= n.sample){
    sample.x[i] <- runif(1,gps.points.range[1,1],gps.points.range[1,2])
    sample.y[i] <- runif(1,gps.points.range[2,1],gps.points.range[2,2])
    samples.pot <- which((sample.x[i]-percent.x*gps.points.range.x)<=gps.points[,1] & gps.points[,1]<=(sample.x[i]+percent.x*gps.points.range.x) &
                           (sample.y[i]-percent.y*gps.points.range.y)<=gps.points[,2] & gps.points[,2]<=(sample.y[i]+percent.y*gps.points.range.y))
    if (length(samples.pot)>=n.minimum.latrine){
      i=i+1
      samples.cand <- sample(samples.pot,n.sample.pooling,replace=FALSE)
      while (any(samples.cand %in% samples)){
        samples.cand <- sample(samples.pot,n.sample.pooling,replace=FALSE)
      }
      samples <- c(samples,samples.cand)
    }
  }
  samples <- as.data.frame(cbind(rep(1:n.sample,each=n.sample.pooling),as.numeric(samples),rep(sample.x,each=n.sample.pooling),rep(sample.y,each=n.sample.pooling)),col.names=c("PSU","ID","sample.x","sample.y"))
  names(samples) <- c("PSU","ID","x","y")
  
  #analysis
  samples.res <- cbind(samples,t(mat.conc[,samples$ID]))
  psu.res <- aggregate(samples.res[,5:length(samples.res[1,])], by=list(samples.res$PSU),FUN=mean)
  names(psu.res)[1] <- "PSU"
  psu.res.wo.id <- as.matrix(psu.res[,-1])
  psu.res.wo.id[which(psu.res.wo.id<LLOD.test)] <- 0
  psu.res.wo.id[which(psu.res.wo.id>=LLOD.test)] <- 1
  psu.res <- cbind(psu.res$PSU,rowMeans(psu.res.wo.id),psu.res.wo.id)
  
  res.test <- c(res.test,psu.res[,2]*100)
}

par(mfrow=c(2,1))
hist(vec.prob.pos*100, breaks = 50, probability = T,main="Histogram of percent of detection",
     xlab="Percent of detection (%)",ylab="Density",cex.lab=1.3,cex.axis=1.5,cex.main=2,
     xlim=c(0,15),ylim=c(0,0.5))
hist(res.test, breaks = 50, probability = T,main="Histogram of percent of detection (Pooling)",
     xlab="Percent of detection (%)",ylab="Density",cex.lab=1.3,cex.axis=1.5,cex.main=2,
     xlim=c(0,15),ylim=c(0,0.5))
