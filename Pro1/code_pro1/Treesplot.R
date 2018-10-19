library(ggplot2)
library(gridExtra)
library(nLTT)
library(DDD)
library(gridGraphics)

library(grid)
source(paste0(getwd(), '/multiplot.R' ))
source(paste0(getwd(), '/event_matrix.R'))
source(paste0(getwd(), 'Plottree_single_Pro1.R'))

# par(mfrow = c(3, 5))
# par(cex = 0.6)
# par(mar = c(0, 0, 0, 0), oma = c(4, 4, 0.5, 0.5))
# par(tcl = -0.25)
# par(mgp = c(2, 0.6, 0))
Treesplot = function(sce){
j = 2
age = 15
i = 1
plist = list()

plotA = matrix(c(1:20), 5, 4)
plotB = rep(21,4)
layout_matrix = rbind(plotA,plotB)
 # layout(matrix(c(1:20), 5, 4))
 layout(layout_matrix, heights=c(2,2,2,2,2,1))
 
 par(mar = rep(0, 4), oma=c(1, 7, 4, 1), las=1)
  # num_group = c(c(101,104,107:109),c(111,114,117:119),c(121,124,127:129),c(131,134,137:139))
 num_group = c(c(101,104,107:109),c(111,114,117:119),c(121,124,127:129),c(131,134,137:139))+40*(sce-1)
 
 for(num in num_group){
    filename = paste(getwd(),j,"Modeltestgroup",age,"age",num,"/out",i,"sim.Rdata",sep = "")
     plottree(file = filename, ls = 0,dropextinct = T,cycle = F)
     i = i+1
  print(paste("Ploting group",num))
  }
  plot_colors <- c("blue","red", "green")

  
  par(mai=c(0,0,0,0))
  plot.new()
  
  legend(x = "center",legend = c("Endemic species", "Endemic species", "Widespread species"), 
         col=plot_colors, lwd=5, cex=1.3, pt.cex = 1, horiz = TRUE, bty="n")
  
  mtext(substitute(paste(mu," = 0")),at=.1,side=3,outer=T,cex=1.2,line = 0)  
  mtext(substitute(paste(mu," = 0.1")),at=.35,side=3,outer=T,cex=1.2,line = 0)  
  mtext(substitute(paste(mu," = 0.2")),at=.6,side=3,outer=T,cex=1.2,line = 0)  
  mtext(substitute(paste(mu," = 0.4")),at=.85,side=3,outer=T,cex=1.2,line = 0)
  
  # mtext(substitute(M[0]),at=1.01,side=2,outer=T,cex=1.2,line = 0.5)  
  
  mtext(substitute(paste(M[0]," = 1000")),at=.18,side=2,outer=T,cex=1.2,line = 7,adj =0)  
  mtext(substitute(paste(M[0]," = 5")),at=.36,side=2,outer=T,cex=1.2,line = 7,adj =0)  
  mtext(substitute(paste(M[0]," = 1")),at=.55,side=2,outer=T,cex=1.2,line = 7,adj =0)  
  mtext(substitute(paste(M[0]," = 0.15")),at=.74,side=2,outer=T,cex=1.2,line = 7,adj =0)
  mtext(substitute(paste(M[0]," = 0")),at=.93,side=2,outer=T,cex=1.2,line = 7,adj =0)
  
}