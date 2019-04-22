iterative.ssn <- importSSN(paste0("./simul/v2.SimIterative2018-07-17.ssn"))
#generate connection matrix for sewage lines;
n.lines <- length(iterative.ssn@lines)
mat.connect.line <- matrix(0,nrow=n.lines,ncol=n.lines)
for (i in 1:(n.lines-1)){
  for (j in (i+1):n.lines){
    if (all(coordinates(iterative.ssn@lines[[i]])[[1]][1,]==coordinates(iterative.ssn@lines[[j]])[[1]][2,])){
      mat.connect.line[i,j] <- 1
    }
  }
}

#generate upstream matrix for sewage lines;
mat.upstream.line <- mat.connect.line
for (i in n.lines:2){
  up.connect.line <- which(mat.connect.line[,i]==1)[which(mat.connect.line[,i]==1)<i]
  vec.upstream.line <- c(i,which(mat.upstream.line[i,]==1))
  mat.upstream.line[up.connect.line,vec.upstream.line] <- 1
}

mat.upstream.line = mat.upstream.line + diag(n.lines)

#generate connection matrix for points;
n.points <- length(iterative.ssn@obspoints@SSNPoints[[1]]@point.data[,1])
mat.connect.points <- matrix(0,nrow=n.points,ncol=n.points)
for (i in 1:n.points){
  for (j in 1:n.points){
    if (mat.upstream.line[as.numeric(iterative.ssn@obspoints@SSNPoints[[1]]@point.data$rid[i]),as.numeric(iterative.ssn@obspoints@SSNPoints[[1]]@point.data$rid[j])]==1 |
        mat.upstream.line[as.numeric(iterative.ssn@obspoints@SSNPoints[[1]]@point.data$rid[j]),as.numeric(iterative.ssn@obspoints@SSNPoints[[1]]@point.data$rid[i])]==1){
      mat.connect.points[i,j] <- 1
    }
  }
}

#generate upstream matrix for points;
mat.upstream.points <- matrix(0,nrow=n.points,ncol=n.points)
for (i in 1:n.points){
  connect.points <- which(mat.connect.points[i,]==1)
  up.connect.points <- connect.points[which(as.numeric(iterative.ssn@obspoints@SSNPoints[[1]]@point.data$upDist)[connect.points]>as.numeric(iterative.ssn@obspoints@SSNPoints[[1]]@point.data$upDist)[i])]
  mat.upstream.points[i,up.connect.points] <- 1
}
mat.upstream.points = mat.upstream.points + diag(n.points)

#generate upstream distance matrix for points;
mat.upstream.dis.points <- matrix(NA,nrow=n.points,ncol=n.points)
for (i in 1:n.points){
  connect.points <- which(mat.connect.points[i,]==1)
  up.connect.dis.points <- connect.points[which(as.numeric(iterative.ssn@obspoints@SSNPoints[[1]]@point.data$upDist)[connect.points]>=as.numeric(iterative.ssn@obspoints@SSNPoints[[1]]@point.data$upDist)[i])]
  mat.upstream.dis.points[i,up.connect.dis.points] <- as.numeric(iterative.ssn@obspoints@SSNPoints[[1]]@point.data$upDist)[up.connect.dis.points]-as.numeric(iterative.ssn@obspoints@SSNPoints[[1]]@point.data$upDist)[i]
}
mat.upstream.dis.points <- mat.upstream.dis.points/max(mat.upstream.dis.points,na.rm = T)*10 #the longest distance is 10km

#load sewage lines;
sewage.lines <- SpatialLines(iterative.ssn@lines)
## sample data: line lengths
df <- data.frame(len = sapply(1:length(sewage.lines), function(i) gLength(sewage.lines[i, ])))
rownames(df) <- sapply(1:length(sewage.lines), function(i) sewage.lines@lines[[i]]@ID)
## SpatialLines to SpatialLinesDataFrame
sewage.lines.df <- SpatialLinesDataFrame(sewage.lines, data = df)
#sampling
gps.points.range <- iterative.ssn@obspoints@SSNPoints[[1]]@points.bbox
gps.points.range.x <- gps.points.range[1,2]-gps.points.range[1,1]
gps.points.range.y <- gps.points.range[2,2]-gps.points.range[2,1]
gps.points <- iterative.ssn@obspoints@SSNPoints[[1]]@point.coords