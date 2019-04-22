pdf(file=paste0("./",version,"/output/AdapativeSites_",Sys.Date(),".pdf"),width=14,height=10)
plot(sewage.lines.df, lwd = sqrt(iterative.ssn@data$addfunccol)*5, col = "black", xlab="Longitude", ylab="Latitude",
     main = "", asp=0.25, ylim=rev(sewage.lines.df@bbox[2,]), cex.main=2, cex.lab=1.5)
points(iterative.ssn@obspoints@SSNPoints[[1]]@point.coords[,1],
       iterative.ssn@obspoints@SSNPoints[[1]]@point.coords[,2],cex=0.5)
points(0,0,cex=3,pch=5) #pumping station sample location

gps.points.range <- iterative.ssn@obspoints@SSNPoints[[1]]@points.bbox
gps.points.range.x <- gps.points.range[1,2]-gps.points.range[1,1]
gps.points.range.y <- gps.points.range[2,2]-gps.points.range[2,1]
gps.points <- iterative.ssn@obspoints@SSNPoints[[1]]@point.coords

abline(v=gps.points.range[1,1]+c(0,1,2,3)/n.col*gps.points.range.x)
abline(h=gps.points.range[2,1]+c(0,1,2,3)/n.row*gps.points.range.y)

samples <- siteInit(iterative.ssn, n.row, n.col, n.minimum.latrine)

sample.x <- unique(samples$x)
sample.y <- unique(samples$y)

points(iterative.ssn@obspoints@SSNPoints[[1]]@point.coords[samples[,2],1],
       iterative.ssn@obspoints@SSNPoints[[1]]@point.coords[samples[,2],2],
       cex=2,pch=4)
rect(sample.x-percent.x*gps.points.range.x,sample.y-percent.y*gps.points.range.y,
     sample.x+percent.x*gps.points.range.x,sample.y+percent.y*gps.points.range.y,
     lwd = 2)
textxy(0,0,"Pumping Station",cex=1)
# textxy(unique(samples.res$x),unique(samples.res$y),paste0("  ",psu.res[,2]*100,"%  "),cex=1)
legend( x="topright", bty='n', cex=2, c('Toilet','PSU',"Toilet Selected","Pumping Station"),pch=c(1,0,4,5))
dev.off()

#sensitivity plots;
vec.lambda <- c(1,2,5,seq(10,200,by=10))
sens.label <- c("low","high")
sens.index <- c(0,5)
for (k.scenario in 1:4){
  load(paste0("./v5/output/scenario",k.scenario,".rda"))
  for (k.sens in 1:2){
    pdf(file=paste0("./v5/output/latrine_PSU_scenario",k.scenario,"sens",sens.label[k.sens],"_",Sys.Date(),".pdf"))
    par(mar=c(5.1,5.1,4.1,1))
    plot(1, type="n",col="black",xlim=c(0,max(vec.lambda)),ylim=c(0,1),xlab="lambda",ylab="sensitivity of detection",cex.main=2,cex.lab=2,cex.axis=2)
    #polygon(c(vec.lambda,rev(vec.lambda)),c(apply(sens[,,1],2,function(x) quantile(x,probs=0.95)),rev(apply(sens[,,1],2,function(x) quantile(x,probs=0.05)))),col="grey")
    lines(vec.lambda, apply(sens[,,sens.index[k.sens]+1],2,function(x) quantile(x,probs=0.5)),col="black",lty=1,lwd=3)
    lines(vec.lambda, apply(sens[,,sens.index[k.sens]+1],2,function(x) quantile(x,probs=0.95)),col="black",lty=2,lwd=3)
    lines(vec.lambda, apply(sens[,,sens.index[k.sens]+1],2,function(x) quantile(x,probs=0.05)),col="black",lty=2,lwd=3)
    
    #polygon(c(vec.lambda,rev(vec.lambda)),c(apply(sens[,,2],2,function(x) quantile(x,probs=0.95)),rev(apply(sens[,,2],2,function(x) quantile(x,probs=0.05)))),col="skyblue")
    lines(vec.lambda, apply(sens[,,sens.index[k.sens]+2],2,function(x) quantile(x,probs=0.5)),col="skyblue",lty=1,lwd=3)
    lines(vec.lambda, apply(sens[,,sens.index[k.sens]+2],2,function(x) quantile(x,probs=0.95)),col="skyblue",lty=2,lwd=3)
    lines(vec.lambda, apply(sens[,,sens.index[k.sens]+2],2,function(x) quantile(x,probs=0.05)),col="skyblue",lty=2,lwd=3)
    
    #polygon(c(vec.lambda,rev(vec.lambda)),c(apply(sens[,,4],2,function(x) quantile(x,probs=0.95)),rev(apply(sens[,,4],2,function(x) quantile(x,probs=0.05)))),col="brown1")
    lines(vec.lambda, apply(sens[,,sens.index[k.sens]+4],2,function(x) quantile(x,probs=0.5)),col="brown1",lty=1,lwd=3)
    lines(vec.lambda, apply(sens[,,sens.index[k.sens]+4],2,function(x) quantile(x,probs=0.95)),col="brown1",lty=2,lwd=3)
    lines(vec.lambda, apply(sens[,,sens.index[k.sens]+4],2,function(x) quantile(x,probs=0.05)),col="brown1",lty=2,lwd=3)
    
    #legend("bottomright",c("PS","Latrine","PSU"),col=c("black","skyblue","brown1"),lty=1,bty="n",lwd=2)
    
    dev.off()
    
    pdf(file=paste0("./v5/output/latrine_PSU_scenario",k.scenario,"sens",sens.label[k.sens],"_",Sys.Date(),".pdf"))
    par(mar=c(5.1,5.1,4.1,1))
    plot(1, type="n",col="black",xlim=c(0,max(vec.lambda)),ylim=c(0,1),xlab="lambda",ylab="sensitivity of detection",cex.main=2,cex.lab=2,cex.axis=2)
    #polygon(c(vec.lambda,rev(vec.lambda)),c(apply(sens[,,1],2,function(x) quantile(x,probs=0.95)),rev(apply(sens[,,1],2,function(x) quantile(x,probs=0.05)))),col="grey")
    lines(vec.lambda, apply(sens[,,sens.index[k.sens]+1],2,function(x) quantile(x,probs=0.5)),col="black",lty=1,lwd=3)
    lines(vec.lambda, apply(sens[,,sens.index[k.sens]+1],2,function(x) quantile(x,probs=0.95)),col="black",lty=2,lwd=3)
    lines(vec.lambda, apply(sens[,,sens.index[k.sens]+1],2,function(x) quantile(x,probs=0.05)),col="black",lty=2,lwd=3)
    
    #polygon(c(vec.lambda,rev(vec.lambda)),c(apply(sens[,,2],2,function(x) quantile(x,probs=0.95)),rev(apply(sens[,,2],2,function(x) quantile(x,probs=0.05)))),col="skyblue")
    lines(vec.lambda, apply(sens[,,sens.index[k.sens]+4],2,function(x) quantile(x,probs=0.5)),col="brown1",lty=1,lwd=3)
    lines(vec.lambda, apply(sens[,,sens.index[k.sens]+4],2,function(x) quantile(x,probs=0.95)),col="brown1",lty=2,lwd=3)
    lines(vec.lambda, apply(sens[,,sens.index[k.sens]+4],2,function(x) quantile(x,probs=0.05)),col="brown1",lty=2,lwd=3)
    
    #polygon(c(vec.lambda,rev(vec.lambda)),c(apply(sens[,,4],2,function(x) quantile(x,probs=0.95)),rev(apply(sens[,,4],2,function(x) quantile(x,probs=0.05)))),col="brown1")
    lines(vec.lambda, apply(sens[,,sens.index[k.sens]+5],2,function(x) quantile(x,probs=0.5)),col="skyblue",lty=1,lwd=3)
    lines(vec.lambda, apply(sens[,,sens.index[k.sens]+5],2,function(x) quantile(x,probs=0.95)),col="skyblue",lty=2,lwd=3)
    lines(vec.lambda, apply(sens[,,sens.index[k.sens]+5],2,function(x) quantile(x,probs=0.05)),col="skyblue",lty=2,lwd=3)
    
    #legend("bottomright",c("PS","PSU","PS+PSU"),col=c("black","brown1","skyblue"),lty=1,bty="n",lwd=2)
    
    dev.off()
  }
}

#legend.pdf;
plot(1, type="n")
legend("topleft",c("Pumping Station","Latrine","PSU"),col=c("black","skyblue","brown1"),title="Sampling Site",lty=1,bty="n",lwd=3,cex=2.5)
legend("bottomleft",c("Median","5th/95th Percentile"),col="black",title="Curve",lty=c(1,2),bty="n",lwd=3,cex=2.5)

plot(1, type="n")
legend("topleft",c("Pumping Station","PSU","Pumping Station + PSU"),col=c("black","brown1","skyblue"),title="Sampling Site",lty=1,bty="n",lwd=3,cex=2.5)
legend("bottomleft",c("Median","5th/95th Percentile"),col="black",title="Curve",lty=c(1,2),bty="n",lwd=3,cex=2.5)


vec.lambda <- c(10,20,30,40,50)
load("./v6/output/adapt_multiOutbreak_2018-09-12.rda")
for (k.lam in 1:length(vec.lambda)){
  lambda <- vec.lambda[k.lam]
  pdf(file=paste0("./v6/output/adapt_multiOutbreak_",lambda,"_",Sys.Date(),".pdf"),width=8,height=8)
  par(mar=c(5.1,5.1,4.1,1))
  plot(NA,ylim=c(0,100),xlim=c(0,20),
       main=paste0("Lambda = ",lambda),xlab="update",ylab="sensitivity (%)",cex.main=3,cex.lab=3,cex.axis=2)
  lines(1:n.update,apply(res.adaptive[k.lam,,,1],2,function(x) quantile(x,probs=0.5))*100,col="black",lty=1,lwd=3)
  lines(1:n.update,apply(res.adaptive[k.lam,,,1],2,function(x) quantile(x,probs=0.05))*100,col="black",lty=2,lwd=3)
  lines(1:n.update,apply(res.adaptive[k.lam,,,1],2,function(x) quantile(x,probs=0.95))*100,col="black",lty=2,lwd=3)
  lines(1:n.update,apply(res.adaptive[k.lam,,,2],2,function(x) quantile(x,probs=0.5))*100,col="brown1",lty=1,lwd=3)
  lines(1:n.update,apply(res.adaptive[k.lam,,,2],2,function(x) quantile(x,probs=0.05))*100,col="brown1",lty=2,lwd=3)
  lines(1:n.update,apply(res.adaptive[k.lam,,,2],2,function(x) quantile(x,probs=0.95))*100,col="brown1",lty=2,lwd=3)
  lines(1:n.update,apply(res.adaptive[k.lam,,,3],2,function(x) quantile(x,probs=0.5))*100,col="skyblue",lty=1,lwd=3)
  lines(1:n.update,apply(res.adaptive[k.lam,,,3],2,function(x) quantile(x,probs=0.05))*100,col="skyblue",lty=2,lwd=3)
  lines(1:n.update,apply(res.adaptive[k.lam,,,3],2,function(x) quantile(x,probs=0.95))*100,col="skyblue",lty=2,lwd=3)
  #legend("topleft",c("PS","PSU","PS+PSU"),col=c("black","brown1","skyblue"),lty=1,bty="n",lwd=3,cex=2)
  dev.off()
}



