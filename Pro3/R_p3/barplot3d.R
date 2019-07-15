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
stackplot.3d<-function(x,y,z,alpha=1,topcol="#078E53",sidecol="#aaaaaa",mode='m5'){
  if(mode=='m2'){
    z.bot = z[1]
    z.top = z[2]    
  }else if(mode=='m5'){
    z.bot = z[1]
    z.q1 = z[2]
    z.mean=z[3]
    z.q3=z[4]
    z.top = z[5]
  }
  ## These lines allow the active rgl device to be updated with multiple changes
  ## This is necessary to draw the sides and ends of the column separately  
  save <- par3d(skipRedraw=TRUE)
  on.exit(par3d(save))

  if(mode=='m2'){
    ## Determine the coordinates of each surface of the column and its edges
    x1=c(rep(c(x[1],x[2],x[2],x[1]),3),rep(x[1],4),rep(x[2],4))
    z1=c(rep(z.bot,4),rep(c(z.bot,z.bot,z.top,z.top),4))
    y1=c(y[1],y[1],y[2],y[2],rep(y[1],4),rep(y[2],4),rep(c(y[1],y[2],y[2],y[1]),2))
    x2=c(rep(c(x[1],x[1],x[2],x[2]),2),rep(c(x[1],x[2],rep(x[1],3),rep(x[2],3)),2))
    z2=c(rep(c(z.bot,z.top),4),rep(z.bot,8),rep(z.top,8) )
    y2=c(rep(y[1],4),rep(y[2],4),rep(c(rep(y[1],3),rep(y[2],3),y[1],y[2]),2) )
    
    ## These lines create the sides of the column and its coloured top surface
    rgl.quads(x1,z1,y1,col=rep(sidecol,each=4),alpha=alpha,lit=FALSE)
    rgl.quads(c(x[1],x[2],x[2],x[1]),rep(z.top,4),c(y[1],y[1],y[2],y[2]),
              col=rep(topcol,each=4),alpha=1,lit=FALSE) 
    rgl.quads(c(x[1],x[2],x[2],x[1]),rep(z.bot,4),c(y[1],y[1],y[2],y[2]),
              col=rep(topcol,each=4),alpha=1,lit=FALSE) 
    ## This line adds black edges to the column
    rgl.lines(x2,z2,y2,col="#000000",lit=FALSE)    
  }else if(mode=='m5'){
    ## Determine the coordinates of each surface of the column and its edges
    x1=c(rep(c(x[1],x[2],x[2],x[1]),3),rep(x[1],4),rep(x[2],4))
    z1=c(rep(z.bot,4),rep(c(z.bot,z.bot,z.mean,z.mean),4))
    y1=c(y[1],y[1],y[2],y[2],rep(y[1],4),rep(y[2],4),rep(c(y[1],y[2],y[2],y[1]),2))
    x2=c(rep(c(x[1],x[1],x[2],x[2]),2),rep(c(x[1],x[2],rep(x[1],3),rep(x[2],3)),2),
         rep((x[1]+x[2])/2,2))
    z2=c(rep(c(z.bot,z.mean),4),rep(z.bot,8),rep(z.mean,8), z.q1,z.q3)
    y2=c(rep(y[1],4),rep(y[2],4),rep(c(rep(y[1],3),rep(y[2],3),y[1],y[2]),2),
         rep((y[1]+y[2])/2,2))
    
    ## These lines create the sides of the column and its coloured top surface
    ## Side surfaces of the main box
    rgl.quads(x1,z1,y1,col=rep(sidecol,each=4),alpha=alpha,lit=FALSE)
    ## Top and bottom surfaces of the main box
    rgl.quads(c(x[1],x[2],x[2],x[1]),rep(z.mean,4),c(y[1],y[1],y[2],y[2]),
              col=rep(topcol,each=4),alpha=1,lit=FALSE) 
    # rgl.quads(c(x[1],x[2],x[2],x[1]),rep(z.q1,4),c(y[1],y[1],y[2],y[2]),
    #           col=rep(topcol,each=4),alpha=1,lit=FALSE) 
    ## Max and min surfaces
    # rgl.quads(c(x[1],x[2],x[2],x[1]),rep(z.top,4),c(y[1],y[1],y[2],y[2]),
    #           col=rep(topcol,each=4),alpha=.2,lit=FALSE) 
    # rgl.quads(c(x[1],x[2],x[2],x[1]),rep(z.bot,4),c(y[1],y[1],y[2],y[2]),
    #           col=rep(topcol,each=4),alpha=.2,lit=FALSE) 
    ## This line adds black edges to the column
    rgl.lines(x2,z2,y2,col="#000000",lit=FALSE)
    
    # bg.x = c(rep(c(10,10),4),rep(c(10,70),4))
    # bg.y = c(rep(rep(c(80,120,160,200),each = 2),2))
    # bg.z = c(rep(c(-10,-100),4),rep(c(-100,-100),4))
    # rgl.lines(bg.x,bg.y,bg.z,col="#000000",lit=FALSE)
    
    bg.x = rep(c(rep(10,46),seq(10,70,2)),4)
    bg.y = rep(c(80,120,160,200),each=77)
    bg.z = rep(c(seq(-10,-100,-2),rep(-100,31)),4)
    rgl.points(bg.x,bg.y,bg.z,col="#000000",lit=FALSE,size = 0.3)
    
  }
 
}
# Example:
# stackplot.3d(c(0,1),c(0,1),3,alpha=0.6)

## Calls stackplot.3d repeatedly to create a barplot
## z.top is the heights of the columns and must be an appropriately named vector
barplot3d<-function(z,alpha=1,scalexy=10,scalez=1,gap=0.2,mode='m5'){
  ## These lines allow the active rgl device to be updated with multiple changes
  ## This is necessary to add each column sequentially
  if(mode=='m2'){
    if(dim(z)[2] != 2){
      return(print('2 columns are expected!'))
    }
    z=z[,c(1,2)]

  }else if(mode=='m5'){
    if(dim(z)[2]!=5){
      return(print('5 columns are expected!'))
    }
    z=z[,c(1:5)]
    
  }else{
    return(print('Pls specify mode!'))
  }
  
  
  save <- par3d(skipRedraw=TRUE)
  on.exit(par3d(save))
  
  ## Recreate Broad order
  types=c("Low",'Intermediate','High')
  contexts=c("jc1p1","jc1p2","jc1p3","jc1p4","jc1p5","jc1p6",
             "jc2p1","jc2p2","jc2p3","jc2p4","jc2p5","jc2p6",
             "jc3p1","jc3p2","jc3p3","jc3p4","jc3p5","jc3p6")
  typeorder=c()
  for(type in types){
    typeorder=c(typeorder,paste(type,contexts,sep="_"))
  }
  names(z.top)=typeorder
  names(z.bot)=typeorder
  
  ## Reorder data into 6 regions
  neworder=c(1:length(z.top))
  
  ## Define dimensions of the plot 
  dimensions=c(9,6)
  
  ## Scale column area and the gap between columns 
  y=seq(1,dimensions[1])*scalexy
  x=seq(1,dimensions[2])*scalexy
  gap=gap*scalexy
  
  z = z*scalez
  
  ## Set up colour palette
  broadcolors=c("#74B655","#F2EC3C","#3F4F9D")
  colors=as.vector(sapply(broadcolors,rep,18))
  
  ## Scale z.top coordinate
  if(mode=='m2'){


    ## Plot each of the columns
    for(i in 1:dimensions[1]){
      for(j in 1:dimensions[2]){
        # Variable to work out which column to plot
        it=(i-1)*dimensions[2]+j 
        stackplot.3d(c(gap+x[j],x[j]+scalexy),
                     c(-gap-y[i],-y[i]-scalexy),
                     z[neworder[it],],
                     alpha=alpha,
                     topcol=colors[neworder[it]],
                     sidecol=colors[neworder[it]],
                     mode=mode)
      }
    }
  }else if(mode=='m5'){

    ## Plot each of the columns
    for(i in 1:dimensions[1]){
      for(j in 1:dimensions[2]){
        it=(i-1)*dimensions[2]+j # Variable to work out which column to plot; counts from 1:96
        stackplot.3d(c(gap+x[j],x[j]+scalexy),
                     c(-gap-y[i],-y[i]-scalexy),
                     z[neworder[it],],
                     alpha=alpha,
                     topcol=colors[neworder[it]],
                     sidecol=colors[neworder[it]],
                     mode=mode)
      }
    }
  }
  
  ## Set the viewpoint and add axes and labels
  rgl.viewpoint(theta=50,phi=40,fov=0)
  axes3d("y-+",labels=TRUE,at=seq(80,200,40),nticks=4)
  # axis for phi
  zlabels <- c('0','0.5','1')
  axes3d("z+-", labels=zlabels,nticks=3,at=seq(-15,-35,-10))
  # axis for sigma_phi
  xlabels <- c('0','1e2','1e4','1e6','1e8','Inf')
  axis3d("x-+",nticks=6,at=seq(15,65,10),labels=xlabels)
  text3d(matrix(c(-10,95,40,140,80,80,10,-25,10),ncol=3),texts=c('Abundance', 'Phy', 'JC'))
  
}
