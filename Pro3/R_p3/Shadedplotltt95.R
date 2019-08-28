## S3 plotting method for objects of class 'ltt95'
## written by Liam J. Revell 2014, 2019
plot.ltt95<-function(x,...){
  if(hasArg(lend)) lend<-list(...)$lend
  else lend<-3
  old.lend<-par()$lend
  par(lend=lend)
  if(hasArg(log)) log<-list(...)$log
  else log<-attr(x,"log")
  if(hasArg(xaxis)) xaxis<-list(...)$xaxis
  else xaxis<-"standard"
  if(hasArg(shaded)) shaded<-list(...)$shaded
  else shaded<-FALSE
  if(shaded) bg<-if(hasArg(bg)) list(...)$bg else rgb(0,0,1,0.25)
  if(attr(x,"method")=="times"){
    n<-max(x[,1])
    if(xaxis=="negative"){ 
      x[,2:4]<-x[,2:4]-max(x[,2:4])
    }
    if(xaxis=="flipped"){
      x[,2:4]<-max(x[,2:4])-x[,2:4]
      x.lim<-c(max(x[,2:4]),min(x[,2:4]))
    } else x.lim<-range(x[,2:4])
    plot(x[,3],x[,1],lwd=2,xlim=x.lim,
         type=if(attr(x,"mode")=="median") "s" else "l",main=NULL,
         xlab=if(xaxis=="standard") "time from the oldest root" else 
           if(xaxis=="negative") "time from the present" else 
             if(xaxis=="flipped") "time before the present day",
         ylab="lineages",log=if(log) "y" else "")
    if(!shaded){
      lines(x[,2],x[,1],lty="dashed",type="s")
      lines(x[,4],x[,1],lty="dashed",type="s")
    } else { 
      xx<-c(x[1,2],rbind(x[2:nrow(x),2],x[2:nrow(x),2]),
            rbind(x[nrow(x):2,4],x[nrow(x):2,4]),x[1,4])
      yy<-c(rbind(x[1:(nrow(x)-1),1],x[1:(nrow(x)-1),1]),x[nrow(x),1],
            x[nrow(x),1],rbind(x[(nrow(x)-1):1,1],x[(nrow(x)-1):1,1]))
      polygon(xx,yy,border=NA,col=bg)
      lines(x[,3],x[,1],lwd=2,
            type=if(attr(x,"mode")=="median") "s" else "l")
    }
  } else if(attr(x,"method")=="lineages"){
    if(xaxis=="negative") x[,1]<-x[,1]-max(x[,1])
    if(xaxis=="flipped"){
      x[,1]<-max(x[,1])-x[,1]
      x.lim<-c(max(x[,1]),min(x[,1]))
    } else x.lim<-range(x[,1])
    plot(x[,1],x[,3],xlim=x.lim,ylim=c(min(x[,2]),max(x[,4])),lwd=2,
         type=if(attr(x,"mode")=="median") "s" else "l",main=NULL,
         xlab=if(xaxis=="standard") "time from the oldest root" else 
           if(xaxis=="negative") "time from the present" else 
             if(xaxis=="flipped") "time before the present day",
         ylab="lineages",log=if(log) "y" else "")
    if(!shaded){
      lines(x[,1],x[,2],lty="dashed",type="s")
      lines(x[,1],x[,4],lty="dashed",type="s")
    } else { 
      xx<-c(x[1,1],rbind(x[2:nrow(x),1],x[2:nrow(x),1]),
            rbind(x[nrow(x):2,1],x[nrow(x):2,1]),x[1,1])
      yy<-c(rbind(x[1:(nrow(x)-1),2],x[1:(nrow(x)-1),2]),x[nrow(x),2],
            x[nrow(x),4],rbind(x[(nrow(x)-1):1,4],x[(nrow(x)-1):1,4]))
      polygon(xx,yy,border=NA,col=bg)
      lines(x[,1],x[,3],lwd=2,
            type=if(attr(x,"mode")=="median") "s" else "l")
    }
  }
  par(lend=old.lend)
}