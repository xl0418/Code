## This script creates a "legoplot" similar to those produced by the Broad Institute
## The plot shows the relative abundance of each of the 6 possible mutations in the 
## 16 sequence contexts

## Load packages
library(rgl)

#### START OF FUNCTIONS

## Functions modified from the "demo(hist3d)" examples in the rgl package:
# library(rgl)
# demo(hist3d)
## Note; would it have killed the original author to comment their code?

## Draws a single "column" or "stack".
## X and Y coordinates determine the area of the column
## The Z coordinate determines the height of the column
## We include "lit=FALSE" arguments to remove the nasty shiny surfaces caused by lighting
stackplot.3d<-function(x,y,z,alpha=1,topcol="#078E53",sidecol="#aaaaaa"){

  z.bot = z[1]
  z.q1 = z[2]
  z.mean=z[3]
  z.q3=z[4]
  z.top = z[5]

  ## These lines allow the active rgl device to be updated with multiple changes
  ## This is necessary to draw the sides and ends of the column separately  
  save <- par3d(skipRedraw=TRUE)
  on.exit(par3d(save))


  ## Determine the coordinates of each surface of the column and its edges
  x1=c(rep(c(x[1],x[2],x[2],x[1]),3),rep(x[1],4),rep(x[2],4))
  z1=c(rep(z.bot,4),rep(c(z.bot,z.bot,z.mean,z.mean),4))
  y1=c(y[1],y[1],y[2],y[2],rep(y[1],4),rep(y[2],4),rep(c(y[1],y[2],y[2],y[1]),2))
  
  ## The edges of the columns
  x2=c(rep(c(x[1],x[1],x[2],x[2]),2),rep(c(x[1],x[2],rep(x[1],3),rep(x[2],3)),2),
       rep((x[1]+x[2])/2,2))
  z2=c(rep(c(z.bot,z.mean),4),rep(z.bot,8),rep(z.mean,8), z.q1,z.q3)
  y2=c(rep(y[1],4),rep(y[2],4),rep(c(rep(y[1],3),rep(y[2],3),y[1],y[2]),2),
       rep((y[1]+y[2])/2,2))
  
  ## The error lines 
  x3=c(rep((x[1]+x[2])/2,2))
  z3=c( z.q1,z.q3)
  y3=c(rep((y[1]+y[2])/2,2))
  
  ## These lines create the sides of the column and its coloured top surface
  ## Side surfaces of the main box
  rgl.quads(x1,z1,y1,col=rep(sidecol,each=4),alpha=alpha,lit=FALSE)
  ## Top and bottom surfaces of the main box
  rgl.quads(c(x[1],x[2],x[2],x[1]),rep(z.mean,4),c(y[1],y[1],y[2],y[2]),
            col=rep(topcol,each=4),alpha=1,lit=FALSE) 

  ## This line adds black edges to the column
  rgl.lines(x2,z2,y2,col="#000000",lit=FALSE)
  rgl.lines(x3,z3,y3,col="#000000",lit=FALSE,lwd=2)
  

  
}
# Example:
# stackplot.3d(c(0,1),c(0,1),3,alpha=0.6)

## Calls stackplot.3d repeatedly to create a barplot
## z.top is the heights of the columns and must be an appropriately named vector
barplot3d<-function(z,group.dim=c(15,18),no.column=3,no.row=3,alpha=1,scalexy=10,scalez=1,
                    gap=0.2,gap.sce.mode=TRUE,y.intercept=c(30,150,250,400),
                    barcolors=c("#8CD790","#EFDC05","#30A9DE")){
  ## These lines allow the active rgl device to be updated with multiple changes
  ## This is necessary to add each column sequentially
  dimensions=group.dim
  if(dim(z)[2]!=5){
    return(print('5 columns are expected!'))
  }
  if(length(barcolors)!=no.column*no.row){
    barcolors = rep(barcolors[1],no.column)
    print("WARNING: the number of barcolors doesn't match the number of groups! The color
          is uniformed by the first color element!")
  }
  z=z[,c(1:5)]
  
  
  save <- par3d(skipRedraw=TRUE)
  on.exit(par3d(save))
  
 
  ## Define dimensions of the plot 

  if(group.dim[1]%%no.column!=0){
    return(print("Number of groups deviding dim 1 is not an integer!"))
  }
  ## Scale column area and the gap between columns 
  y=seq(1,dimensions[1]+no.column-1)*scalexy
  x=seq(1,dimensions[2]+no.row-1)*scalexy
  gap=gap*scalexy
  
  z = z*scalez
  
  ## Set up colour palette
  broadcolors=barcolors
  colors1 = rep(broadcolors[1:3],dimensions[1]/no.column)
  colors2 = rep(broadcolors[4:6],dimensions[1]/no.column)
  colors3 = rep(broadcolors[7:9],dimensions[1]/no.column)
  colors = rep(c(colors1,colors2,colors3),each = dimensions[2]/no.row)
  
  # colors=as.vector(sapply(broadcolors,rep,dimensions[2]/no.row))
  dim.each.block = c(dimensions[1]/no.column,dimensions[2]/no.row)
  size.each.block = dim.each.block[1]*dim.each.block[2]
  ## Scale z.top coordinate
  
  ## Reorder data into 6 regions
  color.order = c(1:nrow(z))
  neworder.first.i=NULL
  for(i in 1:dim.each.block[1]){
    for(j in 1:dimensions[2]){
      neworder.temp=(j-1)%/%(dimensions[2]/no.row)*size.each.block + 
        ifelse(j%%dim.each.block[2]==0,dim.each.block[2],j%%dim.each.block[2])+
        dim.each.block[2]*(i-1) 
      neworder.first.i = c(neworder.first.i,neworder.temp)
    }
  }
  neworder = c(neworder.first.i,neworder.first.i+90, neworder.first.i+180)
  it=1
  ## Plot each of the columns
  for(i in 1:dimensions[1]){
    for(j in 1:dimensions[2]){
         # Variable to work out which column to plot; counts from 1:96
      color.it = (i-1)*dimensions[2] + j
      if(gap.sce.mode==TRUE){
        gap.sce.i = (i-1)%/%(dimensions[1]/no.column)*scalexy
        gap.sce.j = (j-1)%/%(dimensions[2]/no.row)*scalexy
        
      }
      else{
        gap.sce.i=0
        gap.sce.j=0
        
      }
      stackplot.3d(c(gap+x[j]+gap.sce.j,x[j]+scalexy+gap.sce.j),
                   c(-gap-y[i]-gap.sce.i,-y[i]-scalexy-gap.sce.i),
                   z[neworder[it],],
                   alpha=alpha,
                   topcol=colors[color.order[color.it]],
                   sidecol=colors[color.order[color.it]]
                   )
      it = it +1
    }
  }
  
  dot.gap = 2
  # Background grid lines
  bg.x = rep(c(rep(scalexy,(y[length(y)]+scalexy-10)/dot.gap+1),seq(scalexy,x[length(x)]+scalexy,dot.gap)),length(y.intercept))
  bg.y = rep(y.intercept,each=length(bg.x)/length(y.intercept))
  bg.z = rep(c(seq(-10,-y[length(y)]-scalexy,-dot.gap),rep(-y[length(y)]-scalexy,x[length(x)]/2+1)),length(y.intercept))
  rgl.points(bg.x,bg.y,bg.z,col="#000000",lit=FALSE,size = 0.3)
  
  

  
  ## Set the viewpoint and add axes and labels
  ## theta: the horizontal angle    phi: the vertical angle
  rgl.viewpoint(theta=70,phi=35,fov=30)
  axes3d("y-+",labels=TRUE,
         at=y.intercept,
         nticks=4,lwd=2)
  # axis for phi
  zlabels <- c('0','0.25','0.5','0.75','1')
  axes3d("z+-", labels=zlabels,nticks=length(zlabels),at=seq(-15,-55,-10)*scalexy/10,lwd=2.5)
  # axis for sigma_phi
  xlabels <- c('1',expression(10^-2),expression(10^-4),expression(10^-6),expression(10^-8),'0')
  axis3d("x-+",nticks=6,at=seq(155,205,10)*scalexy/10,labels=xlabels,lwd=2.5)
  text3d(matrix(c(0,245,210,100,0,0,0,-45,20)*scalexy/10,ncol=3),
         texts=c('Diversity',expression(psi), expression(phi) ),
         cex = 2)
  # spatial distances indicators
  text3d(matrix(c(0,0,0,40,120,180,150,150,150,150,150,150,
                  -40,-90,-150,-200,-200,-200)*scalexy/10,ncol=3),
         texts=c(expression(sigma[disp] == 10),
                 expression(sigma[disp] == 1),
                 expression(sigma[disp] == 0.1),
                 expression(sigma[JC] == 10),
                 expression(sigma[JC] == 1),
                 expression(sigma[JC] == 0.1)),
         cex = 1.5)
  
}
