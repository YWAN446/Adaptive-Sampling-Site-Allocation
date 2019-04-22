#exploration of the outcomes;
#Check the correlation of depth of nodes and percent of detecting positive;
#upstream lines;
pdf(file=paste0("./",version,"/output/points_depth_",lambda,"_",Sys.Date(),".pdf"),width=8,height=8)
par(mfrow=c(2,2))
plot(iterative.ssn@obspoints@SSNPoints[[1]]@point.data$shreve,vec.prob.pos,
     xlab="number of upstream lines",ylab="percent of detecting positive")
plot(iterative.ssn@obspoints@SSNPoints[[1]]@point.data$shreve,vec.prob.pos,
     xlab="number of upstream lines",ylab="percent of detecting positive",
     xlim=c(0,50))
lines(lowess(iterative.ssn@obspoints@SSNPoints[[1]]@point.data$shreve,vec.prob.pos), col="red",lwd=2)

#upstream lines;
plot(rowSums(mat.upstream.points),vec.prob.pos,
     xlab="number of upstream latrines",ylab="percent of detecting positive")
plot(rowSums(mat.upstream.points),vec.prob.pos,
     xlab="number of upstream latrines",ylab="percent of detecting positive",
     xlim=c(0,300))
lines(lowess(rowSums(mat.upstream.points),vec.prob.pos), col="red",lwd=2)
dev.off()