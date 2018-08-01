set.seed(3000)
xseq<-seq(-4,4,.01)
densities<-dnorm(xseq, 0.5,4)


par(mfrow=c(1,1), mar=c(3,4,4,2))

plot(xseq, densities, col="darkgreen",xlab="", ylab="Density",
     type="l",lwd=2, cex=2, main="PDF of Standard Normal", cex.axis=.8,
     xlim = c(0,1))

y1 = function(x){
  y = 0.1*x*exp(-0.1*x^2)
  return(y)
}
y2 = function(x){
  y = 0.5*x*exp(-0.5*x^2)
  return(y)
}
curve(y1,0,10)
curve(y2,0,10)
