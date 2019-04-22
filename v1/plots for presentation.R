plot(sewage.lines.df, lwd = sqrt(iterative.ssn@data$addfunccol)*5, col = "black", xlab="x-coordinate", ylab="y-coordinate", 
     main = "", asp=0.25, ylim=rev(sewage.lines.df@bbox[2,]))
points(iterative.ssn@obspoints@SSNPoints[[1]]@point.coords[,1], 
       iterative.ssn@obspoints@SSNPoints[[1]]@point.coords[,2],
       cex=1)


par(mfrow=c(2,2))
x1 <- seq(0,50,by=1)
y1 <- dpois(x1,20)
plot(x1,y1,type="l",main="Number of new infections",xlab="New infection",ylab="Prob",cex.lab=1.3,cex.axis=1.3,cex.main=2)

x2 <- seq(0,356,by=0.01)
y2 <- ((-sin(x2*pi*2/365)+1)/2)
plot(x2,y2,type="l",main="Seasonality",xlab="Day",ylab="Multiplying fraction",cex.lab=1.3,cex.axis=1.3,cex.main=2)

x3 <- seq(0,14,by=1)
y3 <- 5*dnbinom(x3,3,0.4)
plot(x3,y3,type="l",main="Intermittent shedding",xlab="Day of infection",ylab="Prob of shedding",cex.lab=1.3,cex.axis=1.3,cex.main=2)

mu.shed <- 10^8
sigma.shed <- 1
x4 <- seq(0,10,by=0.01)
y4 <- dlnorm(10^x4,log(mu.shed)-0.5*sigma.shed^2,sigma.shed)
plot(x4,y4,type="l",main="Number of S.Typhi shed",xlab="Log10 number of S.Typhi",ylab="Prob",cex.lab=1.3,cex.axis=1.3,cex.main=2)

par(mfrow=c(1,1))
gamma.shape <- 1
gamma.rate <- 0.25
x5 <- seq(0,20,by=0.01)
y5 <- 1-pgamma(x5,shape=gamma.shape,rate=gamma.rate)
plot(x5,y5,type="l",main="Decay and loss curve",xlab="Length unit",ylab="Proportion of surviving",cex.lab=1.3,cex.axis=1.5,cex.main=2)

hist(vec.prob.pos*100, breaks = 50, probability = T,main="Histogram of percent of detection",xlab="Percent of detection (%)",ylab="Density",cex.lab=1.3,cex.axis=1.5,cex.main=2)
