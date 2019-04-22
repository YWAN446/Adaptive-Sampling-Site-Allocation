n.infected <- rpois(n.days,getLambda(1:n.days,lambda,season = FALSE))

#simulate which latrines the shedders will go;
#unequal length
latrine_point <- sort(c(runif(n.latrine-1,0,n.latrine),0,n.latrine))
#equal length
# latrine_point <- 0:n.latrine
latrine_range <- data.frame(start=latrine_point[-length(latrine_point)], end=latrine_point[-1])
latrine_range$length <- latrine_range$end - latrine_range$start

mat.latrine.pos <- matrix(0,nrow=n.days,ncol=n.latrine)
for (t in which(n.infected>0)){
  vec.loc.pos <- sample(x=1:n.latrine,size=n.infected[t],prob=latrine_range$length,replace = TRUE)
  mat.latrine.pos[t,unique(sort(vec.loc.pos))] <- table(vec.loc.pos)
}

#shedding
mat.count.gen <- matrix(0,nrow=n.days,ncol=n.latrine)
array.count.gen <- array(0, dim=c(n.days,n.latrine,n.days.shed))

k.shed <- sort(unique(as.vector(mat.latrine.pos)))
k.shed <- k.shed[which(k.shed>0)]

Mu.shed <- log(mu.shed)-0.5*sigma.shed^2

for (k in 1:n.days.shed){
  mat.count.tmp <- matrix(0,nrow=n.days,ncol=n.latrine)
  for (i in k.shed){
    latrine.shed.i <- which(mat.latrine.pos==i)
    if (i==1){
      mat.count.tmp[latrine.shed.i] <- round(rlnorm(length(latrine.shed.i),Mu.shed,sigma.shed)*rbinom(length(latrine.shed.i),1,5*dnbinom(k,3,0.4)),0)
    } else {
      mat.count.tmp[latrine.shed.i] <- round(rowSums(matrix(rlnorm(length(latrine.shed.i)*i,Mu.shed,sigma.shed)*rbinom(length(latrine.shed.i)*i,1,5*dnbinom(k,3,0.4)),ncol=i)),0)
    }
  }
  array.count.gen[k:n.days,,k] <- mat.count.tmp[1:(n.days+1-k),]
}

mat.count.gen <- rowSums(array.count.gen,dims=2)

#assume decay and lost as pgamma(x,2.5,0.5);

mat.upstream.dis.points.prop <- 1-pgamma(mat.upstream.dis.points,shape=gamma.shape,rate=gamma.rate)
mat.upstream.dis.points.prop[which(is.na(mat.upstream.dis.points.prop))] <- 0
mat.count <- matrix(NA,nrow=n.days,ncol=n.latrine)

mat.count <- round(mat.count.gen%*%t(mat.upstream.dis.points.prop),0)

#prob of detecting positive
vol.dilution <- total.vol * iterative.ssn@obspoints@SSNPoints[[1]]@point.data$addfunccol

mat.conc <- t(apply(mat.count,1,FUN = function(x) x/vol.dilution))
mat.pos <- mat.conc
mat.pos[which(mat.conc<LLOD.test)] <- 0
mat.pos[which(mat.conc>=LLOD.test)] <- 1
vec.prob.pos <- colMeans(mat.pos)

std.vec.prob.pos <- vec.prob.pos/max(vec.prob.pos)

mat.dis.pump.stat <- as.numeric(iterative.ssn@obspoints@SSNPoints[[1]]@point.data$upDist)
mat.dis.pump.stat.prop <- 1-pgamma(mat.dis.pump.stat,shape=gamma.shape,rate=gamma.rate)
mat.dis.pump.stat.prop[which(is.na(mat.dis.pump.stat.prop))] <- 0

total.count <- round(rowSums(sweep(mat.count.gen,MARGIN=2,mat.dis.pump.stat.prop,`*`)),0)
uf.conc <- total.count/total.vol
uf.pos <- uf.conc
uf.pos[which(uf.conc<LLOD.uf.test)] <- 0
uf.pos[which(uf.conc>=LLOD.uf.test)] <- 1
uf.prob.pos <- mean(uf.pos,na.rm = TRUE)
std.uf.prob.pos <- uf.prob.pos/max(vec.prob.pos)

sample.latrine <- siteInitLatrine(iterative.ssn, n.row, n.col)
result.latrine <- siteAnalysisLatrine(sample.latrine, mat.conc, uf.pos)
res.latrine <- mean(result.latrine[[2]])
res.latrine.PS <- mean(result.latrine[[3]])

sample.PSU <- siteInit(iterative.ssn, n.row, n.col, n.minimum.latrine)
sample.PSU <- sample.PSU[which(!is.na(sample.PSU$ID)),]
sample.PSU$PSU <- rep(1:(length(sample.PSU$PSU)/n.sample.pooling),each=n.sample.pooling)
result.PSU <- siteAnalysis(sample.PSU, mat.conc, uf.pos)
res.PSU <- mean(result.PSU[[2]])
res.PSU.PS <- mean(result.PSU[[3]])
