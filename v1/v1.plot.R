#plotting
library(rgeos)
sewage.lines <- SpatialLines(iterative.ssn@lines)

## sample data: line lengths
df <- data.frame(len = sapply(1:length(sewage.lines), function(i) gLength(sewage.lines[i, ])))
rownames(df) <- sapply(1:length(sewage.lines), function(i) sewage.lines@lines[[i]]@ID)


## SpatialLines to SpatialLinesDataFrame
sewage.lines.df <- SpatialLinesDataFrame(sewage.lines, data = df)

pdf(file=paste0("./",version,"/output/res_",lambda,"_",Sys.Date(),".pdf"),width=14,height=10)
plot(sewage.lines.df, lwd = sqrt(iterative.ssn@data$addfunccol)*5, col = "black", xlab="x-coordinate", ylab="y-coordinate", 
     main = paste("Peak number of new shedders",lambda,"per day"), asp=0.25, ylim=rev(sewage.lines.df@bbox[2,]), cex.main=2)
points(iterative.ssn@obspoints@SSNPoints[[1]]@point.coords[,1], iterative.ssn@obspoints@SSNPoints[[1]]@point.coords[,2],cex=sqrt(std.vec.prob.pos)*2)
#textxy(iterative.ssn@obspoints@SSNPoints[[1]]@point.coords[,1], iterative.ssn@obspoints@SSNPoints[[1]]@point.coords[,2], as.character(1:n.latrine))
points(0,0,cex=2,pch=5) #pumping station sample location
textxy(0,0,paste0("Pumping Station: ",uf.prob.pos*100,"%"),cex=1)
legend( x="bottomleft", bty='n', title="Percent of Detecting Positive",
        legend=paste(as.character(round(quantile(vec.prob.pos,probs = c(1,0.95,0.75,0.5,0.25))*100,1)),"% - ",c("max","95th percentile","75th percentile","medium","25th percentile")),
        col="black", cex=1, pch=1, pt.cex=sqrt(quantile(std.vec.prob.pos,probs = c(1,0.95,0.75,0.5,0.25)))*2)
dev.off()
# #plot the probability of detecting positive and the number of upstream points;
# plot((rowSums(mat.upstream.points)-1),vec.prob.pos,
#      xlab="number of upstream points",ylab="probability of detecting positive")
# #plot the probability of detecting positive and the number of upstream lines;
# plot((rowSums(mat.upstream.line)-1),vec.prob.pos,
#      xlab="number of upstream points",ylab="probability of detecting positive")
