n.infected <- rpois(n.days,getLambda(1:n.days,lambda,season = FALSE))

#simulate which latrines the shedders will go;
#unequal length
latrine_point <- sort(c(runif(n.latrine-1,0,n.latrine),0,n.latrine))
#equal length
# latrine_point <- 0:n.latrine
latrine_range <- data.frame(start=latrine_point[-length(latrine_point)], end=latrine_point[-1])
latrine_range$length <- latrine_range$end - latrine_range$start

mat.latrine.pos <- matrix(0,nrow=n.days,ncol=n.latrine)
for (t in 1:n.days){
  if (n.infected[t]>0){
    vec.loc.pos <- runif(n.infected[t],0,n.latrine)
    for (k in sapply(vec.loc.pos,FUN = function(x){which(x>=latrine_range$start & x<latrine_range$end)})){
      mat.latrine.pos[t,k] = mat.latrine.pos[t,k] + 1
    }
  }
}

#shedding
mat.count.gen <- matrix(0,nrow=n.days,ncol=n.latrine)
mat.count.tmp <- matrix(0,nrow=n.days,ncol=n.latrine)
array.count.gen <- array(0, dim=c(n.days,n.latrine,n.days.shed))
for (k in 1:n.days.shed){
  for (i in unique(as.vector(mat.latrine.pos))){
    if (i!=0){
      for (j in which(mat.latrine.pos==i)){
        mat.count.tmp[j] <- sum(round(rlnorm(rbinom(1,i,5*dnbinom(k,3,0.4)),log(mu.shed)-0.5*sigma.shed^2,sigma.shed),0))
      }
    }
  }
  array.count.gen[k:n.days,,k] <- mat.count.tmp[1:(n.days+1-k),]
}

for (i in 1:n.days){
  for (j in 1:n.latrine){
    mat.count.gen[i,j] <- sum(array.count.gen[i,j,])
  }
}

# # This is the old function without autocorrelation in temporal;
# # shedding
# mat.count.gen <- matrix(0,nrow=n.days,ncol=n.points)
# for (i in unique(as.vector(mat.latrine.pos))){
#   if (i!=0){
#     for (j in which(mat.latrine.pos==i)){
#       mat.count.gen[j] <- sum(rlnorm(i,log(mu.shed)-0.5*sigma.shed^2,sigma.shed))
#     }
#   }
# }

# #assume no decay and no lost for S. Typhi;
# mat.count <- matrix(NA,nrow=n.days,ncol=n.points)
# for (i in 1:n.days){
#   for (j in 1:n.latrine){
#     mat.count[i,j] <- sum(mat.upstream.points[j,]*mat.count.gen[i,])
#   }
# }

#assume decay and lost as pgamma(x,2.5,0.5);
mat.upstream.dis.points.prop <- 1-pgamma(mat.upstream.dis.points,shape=gamma.shape,rate=gamma.rate)
mat.upstream.dis.points.prop[which(is.na(mat.upstream.dis.points.prop))] <- 0
mat.count <- matrix(NA,nrow=n.days,ncol=n.latrine)
for (i in 1:n.days){
  for (j in 1:n.latrine){
    mat.count[i,j] <- round(sum(mat.upstream.dis.points.prop[j,]*mat.upstream.points[j,]*mat.count.gen[i,]),0)
  }
}

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
