library(DDD)
library(ggplot2)
library(ggridges)
library(ggthemes)
library(viridis)
library("RColorBrewer")
library(grid)
library(gridExtra)
source(paste0(getwd(),'/g_legend.R'))
source('C:/Liang/Code/Pro3/R_p3/event2L.R', echo=TRUE)
# source('C:/Liang/Code/Pro3/R_p3/barplot3d.R', echo=TRUE)
source('C:/Liang/Code/Pro3/R_p3/multi3dbar.R', echo=TRUE)
moviedir = 'C:/Liang/Googlebox/Research/Project3/replicate_sim_final1/'
dir = 'C:/Liang/Googlebox/Research/Project3/replicate_sim_final1/1e+07/'
scenario = c('Hevent','Mevent','Levent')
sce.short = c('H','M','L')
jclabel = c(0,0.25,0.5,0.75,1)
plabel = c(0,1e2,1e4,1e6,1e8,-1)
diversity.upperlimit = 250
diversity.lowerlimit = 95
lowcol = "#F2F0F7"
highcol = '#54278F'
tas1 = list()
r_df = NULL

for(i_n in c(1:3)){
  sce = scenario[i_n]
  f.name = sce.short[i_n]
  for(i in c(1,4,2,5,3)){
    for(j in c(1:6)){
      ns.vec=NULL
      if(i==1){
        comb=paste0(i,i)
      }else{
        comb=paste0(i,j)
      }
      
      for(rep in c(1:100)){
        
        rname = paste0(dir,'spatialpara1e+07',f.name,comb,'/',sce,comb,'rep',rep,'.csv')
        events = read.csv(rname,header = FALSE)
        ns = tail(events[,2],n=1)
        ns.vec = c(ns.vec, ns)
      }
      quantile1 = quantile(ns.vec)
      r_df = rbind(r_df,c(quantile1,i,j,i_n))
    }
  }
}
z = r_df[,c(1:5)]
z[,1] = 95 # min(z)

# Figure size
r3dDefaults$windowRect <- c(0, 100, 800, 800) 

# Animate the 3d plot and save as a gif
open3d()
barplot3d(z,alpha=0.6)
legend3d("topright", legend = c('Low spatiality','Intermediate spatiality','High spatiality'), pch = 12,
         col = c("#8CD790","#EFDC05","#30A9DE"), cex=1.5, inset=c(0.02))

if (!rgl.useNULL())
  # play3d(spin3d(axis = c(0, 1, 0), rpm = 3), duration = 20)
  movie3d(spin3d(axis = c(0, 1, 0), rpm = 3), duration = 20,dir=moviedir,fps=10)


## Save your images to files if you wish
barplot3d(z,alpha=0.6)
legend3d("topright", legend = c('Low spatiality','Intermediate spatiality','High spatiality'), pch = 12,
         col = c("#8CD790","#EFDC05","#30A9DE"), cex=1.5, inset=c(0.02))

rgl.snapshot(filename=paste0(moviedir,"Species_richness.png"))


testz = z
testz[,1] = 96
# Test 15*18 9 groups
color.plate = c('#DB4D6D','#F596AA','#FEDFE1',
                '#CB4042','#EB7A77','#F19483',
                '#D0104C','#B5495B','#E87A90'
                )
# Figure size
r3dDefaults$windowRect <- c(0, 100, 800, 800) 

# Animate the 3d plot and save as a gif
open3d()
barplot3d(testz,group.dim=c(15,18),alpha=0.7,barcolors = color.plate,
          y.intercept=c(95,120,160,200))
legend3d("topright", legend = c('Low spatiality','Intermediate spatiality','High spatiality'), pch = 12,
         col = c("#8CD790","#EFDC05","#30A9DE"), cex=1.5, inset=c(0.02))

